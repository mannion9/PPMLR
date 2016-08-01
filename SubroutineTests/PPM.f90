module CommonData
implicit none
real,parameter ::  x_min=0.,x_max=5.,t_min=0.,t_max=1.,c=1.,COURANT=.1
integer,parameter :: Nx=64, Nt=256,choice=1,detect=0
!integer,parameter :: Nx=8,Nt=2,choice=0,detect=0
integer,parameter :: i_min=-2,i_max=Nx+3,iMIN=0,iMAX=Nx-1
end module

program main
use CommonData
implicit none
real :: dt,cdt,dx_uni
real,dimension(iMIN:iMAX)   :: a_true
real,dimension(i_min:i_max)   :: a,x,dx
real,dimension(i_min:i_max+1) :: x_12
integer :: i

dx_uni = (x_max-x_min)/real(Nx)

! -------------------------------
! Create cell centers and write to file
! Create cell edges   and write to file
! CALL linspace(x_min-1.*dx_uni,x_max+2.*dx_uni,x,SIZE(x))
! CALL linspace(x_min-1.*(dx_uni/2.),x_max+2.*(dx_uni/2.),x_12,SIZE(x_12)) 
CALL linspace(x_min+.5*dx_uni,x_max+.5*dx_uni,x(iMIN:iMAX),Nx)
CALL linspace(x_min,x_max,x_12(iMIN:iMAX+1),Nx+1)
open(unit=1,file='Output/domain.txt')
open(unit=2,file='Output/cells.txt')
write(1,*) x(iMIN:iMAX)      !(0:Nx-1)
write(2,*) x_12(iMIN:iMAX+1) !(0:Nx)
close(1)
close(1)
!-------------------------------

dx(iMIN:iMAX)  = x_12(iMIN+1:iMAX+1)-x_12(iMIN:iMAX)
dt  = COURANT*MAXVAL(dx(iMIN:iMAX))/c
cdt = c*dt

!------------------------------
! Initialize solutions (This gives all values including ghost cells values.)
CALL initial_data(x(iMIN:iMAX),a(iMIN:iMAX))
CALL initial_data(x(iMIN:iMAX),a_true)
!-----------------------------

!------------------------------
! Time Marching
open(unit=1,file='Output/a.txt')
open(unit=2,file='Output/true.txt')
open(unit=3,file='Output/uL.txt')
open(unit=4,file='Output/uR.txt')
open(unit=7,file='Output/u12.txt')
do i=0,Nt
	write(1,*) a(iMIN:iMAX)
	write(2,*) a_true
	! Apply boundary conditions
	CALL boundaries(dx,a,1)
	CALL initial_data(x(iMIN:iMAX)-i*(c*dt),a_true) ! Create analyitic solution
	CALL forward(cdt,dx,a)              ! Create numerical solution
end do
close(1)
close(2)
close(3)
close(4)
close(7)
end program main 

subroutine forward(cdt,dx,a)
!----------------------------------
! Time march one step forward
! Input:
!   N    - Number of elements in solution
!   cdt  - Speed time dt
!   dx   - Cell widths
! Input/Output: 
!   a    - Current solution 
!   a_first- First order solution
use CommonData
implicit none
real,intent(in)    :: cdt
real,intent(inout),dimension(i_min:i_max) :: a,dx
real,dimension(i_min:i_max) :: aL,aR,dela,a6
real,dimension(i_min:i_max) :: x,fp  
CALL interpol(a,dx,aR,aL,dela,a6)
x = cdt/dx
fp = aR - .5*x*( dela -a6*(1-(2./3.)*x))
a(iMIN:iMAX) = a(iMIN:iMAX) + (cdt/dx(iMIN:iMAX))*(fp(iMIN-1:iMAX-1)-fp(iMIN:iMAX))

write(3,*) aL(iMIN-1:iMAX)
write(4,*) aR(iMIN-1:iMAX)
end subroutine

subroutine interpol(a,dx,aR,aL,dela,a6)
!-------------------------------
! Calculate interpolation value of a_(j+1/2)
! Input:
!   N  : Number of elements in dx
!   af : Solution on extended domain
!   dx : Cell widths on extended domain
! Input/Output:
!   a12    : Interpolated value
! GIVES YOU aR,aL... that are true on iMIN:iMAX+1
!------------------------------
use CommonData
implicit none
real,intent(in),dimension(i_min:i_max)   :: a,dx
real,intent(inout),dimension(i_min:i_max) :: aR,aL,dela,a6
real,dimension(iMIN-1:iMAX) :: one,two,three,four,five
real,dimension(i_min:i_max) :: a12,del,aRd,aLd
integer :: js,je
js = iMIN-1 ! First cell I interpolate is the i=-1  cell
je = iMAX+1   ! Last cell  I interpolate is the i=Nm-1 cell
one   = dx(js:je)/(dx(js:je)+dx(js+1:je+1))
two   = 1./(dx(js-1:je-1)+dx(js:je)+dx(js+1:je+1)+dx(js+2:je+2))
three = (2.*dx(js+1:je+1)*dx(js:je))/(dx(js:je)+dx(js+1:je+1))
four  = (dx(js-1:je-1)+dx(js:je))/(2.*dx(js:je)+dx(js+1:je+1))
five  = (dx(js+2:je+2)+dx(js+1:je+1))/(2.*dx(js+1:je+1)+dx(js:je))

CALL avgSlope(a,dx,del)
a12(js:je) = a(js:je) + one*(a(js+1:je+1)-a(js:je)) &
			+ two*(three*(four-five)*(a(js+1:je+1)-a(js:je)) &
            - dx(js:je)*four*del(js+1:je+1)+dx(js+1:je+1)*five*del(js:je))	
write(7,*) a12		
aR(iMIN:iMAX+1) = a12(iMIN:iMAX+1)
aL(iMIN:iMAX+1) = a12(iMIN-1:iMAX)
if (detect==1) then
	call discont_detect(a,dx,aR,aRd,aL,aLd,del)
end if 
CALL monotone(a,aL,aR)
dela = aR-aL
a6   = 6.*(a-.5*(aL+aR))			
end subroutine

subroutine discont_detect(a,dx,aR,aRd,aL,aLd,del)
use CommonData
implicit none
real,intent(in),dimension(i_min:i_max) :: a,dx,del
real,intent(inout),dimension(i_min:i_max) :: aL,aR
real,dimension(i_min:i_max) :: aLd,aRd,delsquare,nabla
real :: eps = 0.01,nabla1=20.,nabla2=0.05, nabla_temp 
integer :: i,js,je
js = iMIN-1
je = iMAX+2
aRd(iMIN:iMAX+1) = a(iMIN+1:iMAX+2) - .5*del(iMIN+1:iMAX+2)
aLd(iMIN:iMAX+1) = a(iMIN-1:iMAX)   + .5*del(iMIN-1:iMAX)
delsquare(js:je) = 1./(dx(js-1:je-1)+dx(js:je)+dx(js+1:je+1)) &
				 * ((a(js+1:je+1)-a(js:je))/(dx(js+1:je+1)+dx(js:je)) &
				 - (a(js:je)-a(js-1:je-1))/(dx(js:je)+dx(js-1:je-1)))
do i=iMIN,iMAX+1
	if (-1.*delsquare(i+1)*delsquare(i-1).GT.0. .AND. ABS(a(i+1)-a(i-1))-eps*MIN(ABS(a(i+1)),ABS(a(i-1))).GT.0.) then
		nabla_temp = -1.*((delsquare(i+1)-delsquare(i-1))/(dx(i)+dx(i-1)))*((dx(i-1)**(3)+dx(i)**(3))/(a(i+1))-a(i-1))
	else 
		nabla_temp = 0.0
	end if 
	nabla(i) = MAX(0.0,MIN(nabla1*(nabla_temp-nabla2),1.))
end do 
aL = aL*(1.-nabla)+aLd*nabla
aR = aR*(1.-nabla)+aRd*nabla
end subroutine

subroutine avgSlope(af,dxf,delmj)
!------------------------------
! Calculate the average slope in the jth cell of
! the parabola with zone average a j,j+1,j-1
! Input:
!   j  - Index ranging from 0 to 1
!   N  - Number of elements in dx
!   af - Solution on exteneded domain
!   dxf- Cell widths on extended domain
! Input/Output:
!   delj - Average slove 
!-----------------------------
use CommonData
implicit none
real,intent(in),dimension(i_min:i_max) :: af,dxf
real,intent(inout),dimension(i_min:i_max) :: delmj
real,dimension(iMIN-1:iMAX+1) :: one,two,three,diff1,diff2,delj
integer :: st0,stm1,stp1,en0,enm1,enp1,i
real :: sgn
st0 = iMIN-1 ! goes from -1 ghost cell
en0 = iMAX+2 ! up untile the +1 ghost cell 
stm1= st0 -1
enm1= en0 -1
stp1= st0 +1
enp1= en0 + 1
one   = dxf(st0:en0)/(dxf(stm1:enm1)+dxf(st0:en0)+dxf(stp1:enp1))
two   = (2.*dxf(stm1:enm1)+dxf(st0:en0))/(dxf(stp1:enp1)+dxf(st0:en0))
three = (dxf(st0:en0)+2.*dxf(stp1:enp1))/(dxf(stm1:enm1)+dxf(st0:en0)) 
delj  = one*(two*(af(stp1:enp1)-af(st0:en0))+three*(af(st0:en0)-af(stm1:enm1)))
diff1 = af(stp1:enp1)-af(st0:en0)
diff2 = af(st0:en0)-af(stm1:enm1)
do i=st0,en0
  if (diff1(i)*diff2(i) > 0.) then
    if (delj(i)<0.) then
      sgn = -1.
    else 
      sgn = 1.
    end if 
    delmj(i) = SIGN( MIN(ABS(delj(i)),2.*ABS(diff1(i)),2.*ABS(diff2(i))),sgn)
  else 
    delmj(i) = 0. 
  end if 
end do
end subroutine

subroutine monotone(a,aR,aL)
!-------------------------------
! Ensures interpolation is monotone
! Input:
!   N  - Size of a
!   a  - Current solution
! Input/Output:
!   aR - Approx. of a_(j+1/2) from cell to right
!   aL - Approx. of a_(j-1/2) from cell to left
!-----------------------------
use CommonData
implicit none
real,intent(in),dimension(i_min:i_max) :: a
real,intent(inout),dimension(i_min:i_max) :: aR,aL
real,dimension(i_min:i_max) :: aR_temp,aL_temp, cond1,cond2,cond3
integer :: i
aR_temp = aR
aL_temp = aL
cond1 = (aR_temp-a)*(a-aL_temp)
cond2 = (aR_temp-aL_temp)*(a-.5*(aL_temp+aR_temp))
cond3 = ((aR_temp-aL_temp)**(2))/6.
do i=iMIN-1,iMAX
  if (cond1(i).LE.0.) then
    aL(i) = a(i)
    aR(i) = a(i)
  end if
  if (cond2(i)>cond3(i)) then
    aL(i) = 3.*a(i)-2.*aR_temp(i)
  end if 
  if (-1.*cond3(i)>cond2(i)) then
    aR(i) = 3.*a(i)-2.*aL_temp(i)
  end if 
end do 
end subroutine

subroutine boundaries(dx,a,BC)
use CommonData
implicit none
integer,intent(in) :: BC
real,intent(inout),dimension(i_min:i_max) :: dx,a
integer i
! YOU CAN MAKE THIS EASIER TO READ
if (BC==0) then ! Reflective
	 do i = 1,ABS(i_min)
		 dx(iMIN-i) = dx(iMIN+i-1)
		 a(iMIN-i)  = a(iMIN+i-1)
	 end do 
	 do i = iMAX,i_max
		 dx(i+1) = dx(iMAX-(i-iMAX))
		 a(i+1)  = a(iMAX-(i-iMAX))
	end do 
else if (BC==1) then ! Periodic
	 do i = 1,ABS(i_min)
		 dx(iMIN-i) = dx(iMAX-i+1)
		 a(iMIN-i)  = a(iMAX-i+1)
	 end do
	 do i = iMAX+1,i_max
		 dx(i) = dx(iMIN+i-iMAX-1)
		 a(i)  = a(iMIN+i-iMAX-1)
	end do 
end if 
end subroutine

subroutine initial_data(x,y)
!--------------------------------
! Assigns array y values f(x) where
! f is the intial condition function.
! Input:
!   choice - Choice of intial condition
!            Set to 1 -- Gaussian
!            Set to 0 -- Square Wave
!   x      - Domain of cell ceneters
! Input/Output:
!   y      - Array to fill with intial condtion 
!--------------------------------
use CommonData
implicit none
real,intent(in),dimension(iMIN:iMAX) :: x
real,intent(inout),dimension(iMIN:iMAX) :: y
real :: mean=1.5,sigma=.1,x_L=1.,x_R=2.
integer :: i

do i=iMIN,iMAX
  if (choice==0) then
    if (x(i).LE.x_R .AND. x(i).GE.x_L) then
      y(i) = 1.
    else 
      y(i) = 0.
    end if 
  else if (choice==1) then
  	y(i) = EXP(-(x(i)-mean)**(2)/(2.*sigma))
  else if (choice==2) then
	  y(i)= 1.
  end if 
end do
end subroutine

subroutine linspace(xmin,xmax,x,N)
!---------------------------------
! Assigns array x values, that are 
! linearly spaced. For an array with indexing
! [0,N-1] the input N should be N.
! Input:
!   xmin - Minimum value
!   xmax - Maximum value
!   N    - Number of elements in x
! Input/Output:
!   x    - Array with N elements 
!-------------------------------
real,intent(in)    :: xmin,xmax
integer,intent(in) :: N
real,intent(inout),dimension(0:N-1) :: x
real    :: dx
integer :: i
dx = (xmax-xmin)/real(N-1)
do i=0,N-1
  x(i) = xmin + i*dx
end do
end subroutine