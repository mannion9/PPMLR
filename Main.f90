! Ideas to spped up code
! 1) When calling interpol, just put uR as the last input instead of a_12,
! unless you want to write out a_12

module CommonData
! Nm    - Number of cells in true domain
! Nt    - Bound on number of time steps (Needed for plotting against excact)
! i_min - Index of left  most cell of full domain (including ghost cells)
! i_max - Index of right most cell of full domain (including ghost cells)
! iMIN  - Index of left  most cell of true domain
! iMAX  - Index of right most cell of true domain 
! rMIN  - Cell center of left  most real cell
! rMAX  - Cell center of right most real cell
! t_max - Bound of time in domain
! gm    - Ratio of specific heats (ideal gas -> gm=1.4)
! alpha - Dimension parameter
! ap1   - alpha+1
! delta - Machine precision on zero test
! COURANT- Courant number
! dr_uni- Unifrom Eulerian cell widths
! r_min - Left most ghost cell center
! r_max - Right most ghost cell center
implicit none
!integer,parameter :: Nm = 6,Nt = 6
integer,parameter :: Nm=100 ,Nt = 500
integer,parameter :: i_min=-3,i_max=Nm+3 , iMIN=0,iMAX=Nm-1 !Nm is number of interior real
real,parameter    :: rMIN=0.,rMAX=1.,t_max=.2,gm = 1.4,alpha=0.,ap1 = alpha+1.,  &
					 delta=1.0d-30 , COURANT = 0.1 , dr_uni = (rMAX-rMIN)/REAL(Nm-1), &
					 r_min = rMIN -dr_uni*(iMIN-i_min) , &
					 r_max = rMAX +dr_uni*(i_max-iMAX)
end module

program main
use CommonData
implicit none
real,dimension(i_min:i_max)   :: dm,dr,r
real,dimension(i_min:i_max+1) :: r_12_i,r_12_ip
real,dimension(3,i_min:i_max) :: U,V
real,dimension(3,i_min:i_max) :: aL,dela,a6
real,dimension(i_min:i_max) :: itters  ! DO NOT NEED JUST FOR DEBUGGING
real :: dt,total_mass,time=0.0
integer :: i

! Create cell center,edges,and initial conditions


CALL linspace(r_min,r_max,r,SIZE(r)) 						   ! All cell centers (including ghost cells)
CALL linspace(r_min-.5*dr_uni,r_max+.5*dr_uni,r_12_i,SIZE(r_12_i))  ! All cell edges   (including ghost cells)
open(unit=1,file='Output/CellCenter.txt')
write(1,*) r(iMIN:iMAX)
close(1)
! Initialize cell width and initial conditions and total mass inside real cells

dr(iMIN:iMAX) = r_12_i(iMIN+1:iMAX+1)-r_12_i(iMIN:iMAX)
CALL intialCondition(r,U,V,1)
dm = dr*V(1,:)
total_mass = SUM(dm(iMIN:iMAX))  ! rho*dV
print*,'Total Mass:',total_mass

open(unit=1,file='Output/Density.txt')
open(unit=2,file='Output/Velocity.txt')
open(unit=3,file='Output/InternalEnergy.txt')
open(unit=4,file='Output/Pressure.txt')
open(unit=5,file='Output/Energy.txt')
open(unit=7,file='Output/dt.txt')

do i = i_min,i_max
	itters(i) = real(i)
end do 
print*,'itt :',itters(iMIN:iMAX)
i=0

do while (time.LE.t_max .AND. i.LE.NT)
	i=i+1
	CALL boundaries(dr,U,V,0)
	r_12_ip = r_12_i  
	dm = dr*V(1,:)
	CALL TimeStep(U,V,dr,dt)
	IF (ABS(total_mass-SUM(dm(iMIN:iMAX)))/total_mass.GE. .1) THEN
		PRINT*,'Conservation of mass broken at itteration:',i,'at time:',time
		exit
	END IF 
	time = time+dt
	CALL Output(U,V,dm,dt,i,1,0) ! Write out
	CALL LagrangeStep(U,V,dr,dm,dt,r_12_i,r_12_ip) ! updates  tau,u,E
	CALL Output(U,V,dm,dt,i,0,0)
	CALL remap1(dt,dr,r_12_i,r_12_ip,U,V,dm)
	CALL Output(U,V,dm,dt,i,0,0)
end do 
print*,'Finished at time:',time,',after:',i,'itterations:'
close(1)
close(2)
close(3)
end program

subroutine remap1(dt,dr,r_12_i,r_12_ip,U,V,dm)
use CommonData
implicit none
real,intent(in) :: dt
real,intent(in),dimension(i_min:i_max) :: dr
real,intent(in),dimension(i_min:i_max+1) :: r_12_i,r_12_ip
real,intent(inout),dimension(3,i_min:i_max) :: U,V
real,intent(inout),dimension(i_min:i_max) :: dm
real,dimension(i_min:i_max) :: rhoR,rhoL,delrho,rho6,rho12, &
							   uR  ,uL  ,delu  ,u6  ,u12, &
							   pR  ,pL  ,delp  ,p6  ,p12, &
							   ER  ,EL  ,delE  ,E6  ,E12, & 
							   dm_ip , dx_ip , mf=0.d0, &
							   uf=0.d0,Ef=0.d0
real,dimension(i_min:i_max+1):: deltax

real :: y , twth = 2./3.
integer :: i,j
! Calculate cell widths and overlap
dx_ip   = r_12_ip(i_min+1:i_max+1)-r_12_ip(i_min:i_max)  ! New lagrange zone widths
deltax  = r_12_ip-r_12_i				 ! Overlap between new and old zone edges
! Interpolate rho,u,P in lagrange cell to attain an interpolation of E

call interpol(V(1,:),dx_ip,rho12,rhoR,rhoL,delrho,rho6,0)  !!!!!!! SET THE LAST TERM TO 1 
call interpol(U(2,:),dx_ip,u12  ,uR  ,uL  ,delu  ,u6  ,0)
call interpol(V(3,:),dx_ip,p12  ,pR  ,pL  ,delp  ,p6  ,0)
! Interpolate energy, given the P,rho,u interpolation vlaues
E12=P12/((gm-1.)*rho12) + .5*u12**(2)
ER(iMIN-1:iMAX+1) = E12(iMIN-1:iMAX+1)
EL(iMIN-1:iMAX+1) = E12(iMIN-2:iMAX)
delE = ER-EL
E6   = 6.*(U(3,:)-.5*(ER+EL))

do i = iMIN, iMAX+1
	if(deltax(i) >= 0.0) then
		j = i - 1
		y = deltax(i)/dx_ip(j)
		mf(i)   = (rhoL(j) - 0.5*y*(delrho(j) - (1.-twth*y)*rho6(j)) )*deltax(i)
		uf(i)   = (uL(j)   - 0.5*y*(delu(j)   - (1.-twth*y)*u6(j))   )*mf(i)
		Ef(i)   = (EL(j)   - 0.5*y*(delE(j)   - (1.-twth*y)*E6(j))   )*mf(i)
	else
		y = deltax(i)/dx_ip(i)
		mf(i) = (rhoL(i) + 0.5*y*(delrho(i) + (1.-twth*y)*rho6(i))    )*deltax(i)
		uf(i) = (uL(i)   + 0.5*y*(delu(i) 	+ (1.-twth*y)*u6(i))      )*mf(i)
		Ef(i) = (EL(i)   + 0.5*y*(delE(i) 	+ (1.-twth*y)*E6(i))      )*mf(i)
	endif
enddo

do i = iMIN, iMAX
  dm_ip(i) = (dm(i) + mf(i) - mf(i+1))
  V(1,i)   = dm_ip(i)/dr(i)
  V(1,i)   = max(0.00001,V(1,i))
  dm_ip(i) = (V(1,i)*dr(i))
  U(2,i)   = (U(2,i)*dm(i) + uf(i)-uf(i+1))/dm_ip(i)
  U(3,i)   = (U(3,i)*dm(i) + Ef(i)-Ef(i+1))/dm_ip(i)
enddo


U(1,:) = 1./V(1,:)
do i=0,Nm-1
	V(2,i) = MAX(-0.000001,U(3,i) - .5*U(2,i)**(2)) 
end do 
V(3,:) = (gm-1.)*V(1,:)*V(2,:)

dm = dm_ip


end subroutine

subroutine LagrangeStep(U,V,dr,dm,dt,r_12_i,r_12_ip)
!-------------------------------------------------
! Time marches forward one lagrange step and will
! updates r_12_ip,U(1,:),U(2,:),U(3,:)
! Input:
!   dm     - Cell widths
!   dr     - Length of cells
! Input/Output:
!   U  - Array of conservec variables (tau,u,E)
!   V      - Array of primative variables (rho,e,P)
!   r_12_i - Space coordinate before step
!   r_12_ip- Space coordinate after  step 
!--------------------------------------------------
use CommonData
implicit none
real,intent(in) :: dt
real,intent(in),dimension(i_min:i_max)       :: dm,dr
real,intent(inout),dimension(3,i_min:i_max)  :: U,V
real,intent(inout),dimension(i_min:i_max+1)  :: r_12_i,r_12_ip
real,dimension(i_min:i_max)  :: u_star,p_star,Abar  ! Setting to 4. for debugging purposes
real,dimension(i_min:i_max)    :: pL=4.,pR,uL=4.,uR,c_p,c_m
integer :: i
real :: test  ! needed as a real to test if u_Star is zero (fortran formality)

call RiemannEffectiveStates(U,V,r_12_i,pL,pR,uL,uR,c_p,c_m,dt,0)
call RiemannSolver(pL,pR,uL,uR,c_p,c_m,p_star,u_star)

print*,' '
! print*,'pL:',pL(iMIN:iMAX+1)
! print*,'pR:',pR(iMIN:iMAX+1)
! print*,'uL:',uL(iMIN:iMAX+1)
! print*,'uR:',uR(iMIN:iMAX+1)
! print*,'p*:',p_star(iMIN:iMAX+1)
! print*,'u*:',u_star(iMIN:iMAX+1)
!-------------------------------------
! print*,'Rho :',V(1,:)
! print*,'RhoR:',rhoR
! print*,'RhoL:',rhoL
! print*,'u   :',U(2,:)
! print*,'uR  :',uR
! print*,'uL  :',uL
! print*,'P   :',V(3,:)
! print*,'PR  :',pR 
! print*,'PL  :',pL 

r_12_ip(iMIN:iMAX+1)  = r_12_i(iMIN:iMAX+1) + dt*u_star(iMIN:iMAX+1)
do i=i_min,i_max    ! This code is needed b/c if u_star==0 then Abar in floating point goes to zero, but in exact it would go to zero (plug in r_j+12^(n+1) in the def. of Abar you will see goes to 1)
	if (u_star(i).LE.TINY(test)) then
		Abar(i) = 1.
	else 
		Abar(i)  = (r_12_ip(i)**(ap1)-r_12_i(i)**(ap1))/(ap1*u_star(i)*dt)
	end if 
end do
U(1,iMIN:iMAX) = ((r_12_ip(iMIN+1:iMAX+1)**(ap1)-r_12_ip(iMIN:iMAX)**(ap1)))/(ap1*dm(iMIN:iMAX))      ! Tau
U(2,iMIN:iMAX) = U(2,iMIN:iMAX) + &
				.5*(Abar(iMIN+1:iMAX+1)+Abar(iMIN:iMAX))*(dt/dm(iMIN:iMAX))*(p_star(iMIN:iMAX)-p_star(iMIN+1:iMAX+1)) ! U
U(3,iMIN:iMAX) = U(3,iMIN:iMAX) + (dt/dm(iMIN:iMAX))*(Abar(iMIN:iMAX)*u_star(iMIN:iMAX)*p_star(iMIN:iMAX) &
								- Abar(iMIN+1:iMAX+1)*u_star(iMIN+1:iMAX+1)*p_star(iMIN+1:iMAX+1))

! Update primatives
V(1,:) = 1./U(1,:)		     ! Rho
V(2,:) = U(3,:) - .5*U(2,:)**(2) ! Internal energy
V(3,:) = (gm-1.)*V(1,:)*V(2,:)   ! Pressure

end subroutine

subroutine RiemannSolver(pL,pR,uL,uR,c_p,c_m,p_star,u_star)
!---------------------------------------
! Solves RP for two shock waves for p_star 
! and then solves for u_star.
! Input:
!   pL/pR - Pressure approx. array
!   uL/uR - Velocity approx. array
!   c_pm  - Effective speed of sound arrays
! Input/Output:
!   p_star- Pressure between cells
!   u_star- Velocity between cells
! GIVES YOU p_star,u_star on iMIN:iMAX+1
!---------------------------------------
use CommonData
implicit none
real,intent(in),dimension(i_min:i_max)   :: pL,pR,uL,uR,c_m,c_p
real,intent(inout),dimension(i_min:i_max):: p_star,u_star  
real :: TOL=10.**(-6),diff,Temp,Wl,del_u,f,f_var,fprime
real :: gm1=gm-1,gmp1=gm+1
integer :: j, i, itterMAX = 20
p_star = .5*(pL+pR)! initial guess value, the addition of the 10^(-5) helps if your guess is exact 
do j=iMIN,iMAX+1  ! For each cell interface
  diff = 1. ! just intialize doesnt matter yet
  do i=1,itterMAX
    if (diff > TOL) then
      Temp  = p_star(j)
      del_u = uR(j)-uL(j)
      call Func_RP(c_p(j),c_m(j),p_star(j),pL(j),pR(j),del_u,f)
      call Func_RP(c_p(j),c_m(j),p_star(j)+10.**(-6),pL(j),pR(j),del_u,f_var)
      fprime = (f_var-f)/10.**(-6)
      if (fprime.LE.delta .OR. fprime.NE.fprime .OR. p_star(j).LT.delta) then ! ensures no divide by zero or NAN or negative pressure
		p_staR(j) = Temp
		exit
	  end if 
	  p_star(j) = p_star(j) - f/fprime 
      diff = 2.*ABS(p_star(j)+Temp)/(p_star(j)-Temp)
    end if
  end do
Wl = SQRT(ABS(c_p(j)*(1.+((gm+1.)/(2.*gm))*((p_star(j)/pL(j))-1.))))
u_star(j)  = -(p_star(j)-pL(j))/Wl + uL(j)  
end do 

!print*,'Solution to RP'
!print*,'cp:',c_pm(1,:)
! print*,'cm:',c_pm(2,:)
! print*,'pL :',pL
! print*,'pR:',pR
! print*,'GUESS:',.5*(pL+pR)
! print*,'p_star:',p_star
! print*,'u_star:',u_star

end subroutine

subroutine RiemannEffectiveStates(U,V,r_12_i,pp,pm,up,um,c_p,c_m,dt,choice)
use CommonData
implicit none
real,intent(in),dimension(3,i_min:i_max) :: U,V
real,intent(in) :: dt
integer,intent(in) :: choice
real,intent(inout),dimension(i_min:i_max) :: pp,pm,up,um,c_p,c_m
real,dimension(i_min:i_max) :: 	 rhoR,rhoL,delrho,rho6,rho12, &
								 uR  ,uL  ,delu  ,u6  ,u12,&
								 pR  ,pL  ,delp  ,p6  ,p12, &
								 cs  , A , dr , rhop , rhom ,z
real,dimension(i_min:i_max+1) :: r_12_i
integer :: i,j
dr = r_12_i(i_min+1:i_max+1)-r_12_i(i_min:i_max)
A  = (r_12_i(i_min+1:i_max+1)**(ap1)-r_12_i(i_min:i_max)**(ap1))/(ap1*dr)
cs = SQRT(gm*V(3,:)*V(1,:))
z  = dt*cs*A/dr

call interpol(V(1,:),dr,rho12,rhoR,rhoL,delrho,rho6,0)
call interpol(U(2,:),dr,u12  ,uR  ,uL  ,delu  ,u6  ,0)
call interpol(V(3,:),dr,p12  ,pR  ,pL  ,delp  ,p6  ,0)
! prinT*,'U :',U(2,iMIN-1:iMAX+2)
! print*,'u12:',u12(iMIN:iMAX+2)
! print*,'uR:',uR(iMIN:iMAX+1)
! print*,'uL:',uL(iMIN:iMAX+1)
if (choice==1) then
	do i=iMIN,iMAX+1  ! For each interface we need a left and right state 
		j = i-1
		pp(i)  = pL(j) - 0.5*z(j)*(delp(j)-(1.-(2./3.)*z(j))*p6(j))
		up(i)  = uL(j) - 0.5*z(j)*(delu(j)-(1.-(2./3.)*z(j))*u6(j))
		rhop(i)= rhoL(j)-0.5*z(j)*(delrho(j)-1.-(2./3.)*z(j))*rho6(j)
		c_p(i) = gm*pp(i)*rhop(i) ! The square of speed of sound

		pm(i)  = pL(i) + 0.5*z(i)*(delp(i)+(1.-(2./3.)*z(i))*p6(i))
		um(i)  = uL(i) + 0.5*z(i)*(delu(i)+(1.-(2./3.)*z(i))*u6(i))
		rhom(i)= rhoL(i)+0.5*z(i)*(delrho(i)+(1.-(2./3.)*z(i))*rho6(i))
		c_m(i) = gm*pm(i)*rhom(i) ! The square of speed of sound
	end do 
else 
	rhop(iMIN:iMAX+1) =  V(1,iMIN-1:iMAX)!1./rhoL(iMIN-1:iMAX)
	rhom(iMIN:iMAX+1) =  V(1,iMIN:iMAX+1)!1./rhoR(iMIN:iMAX+1)
	pp(iMIN:iMAX+1)   =  V(3,iMIN-1:iMAX)!pL(iMIN:iMAX)
	pm(iMIN:iMAX+1)   =  V(3,iMIN:iMAX+1)!pR(iMIN:iMAX)
	up(iMIN:iMAX+1)   =  U(2,iMIN-1:iMAX)!uL(iMIN:iMAX)
	um(iMIN:iMAX+1)   =  U(2,iMIN:iMAX+1)!uR(iMIN:iMAX)
	c_p(iMIN:iMAX+1) = (gm*pp(iMIN:iMAX+1)*rhop(iMIN:iMAX+1)) ! The square of speed of sound 
	c_m(iMIN:iMAX+1) = (gm*pm(iMIN:iMAX+1)*rhom(iMIN:iMAX+1))
end if 
end subroutine 

subroutine func_RP(cp,cm,p,pl,pr,del_u,eval)
!----------------------------------
! Calculates RP funciton to find root of 
!----------------------------------
use CommonData
implicit none
real,intent(in) :: cp,cm,p,pl,pr,del_u
real,intent(inout) :: eval
real :: Wl,Wr
Wl   = SQRT(ABS(cp*(1.+((gm+1.)/(2.*gm))*((p/pl)-1.))))
Wr   = SQRT(ABS(cm*(1.+((gm+1.)/(2.*gm))*((p/pr)-1.))))
eval = (p-pL)/Wl + (p-pR)/Wr - del_u
end subroutine

subroutine interpol(a,dx,a12,aR,aL,dela,a6,detect)
!-------------------------------
! Calculate interpolation value of a_(j+1/2)
! Input:
!   N  : Number of elements in dx
!   af : Solution on extended domain
!   dx : Cell widths on extended domain
! Input/Output:
!   a12    : Interpolated value
!------------------------------
use CommonData
implicit none
integer,intent(in) :: detect
real,intent(in),dimension(i_min:i_max)   :: a,dx
real,intent(inout),dimension(i_min:i_max) :: a12,aR,aL,dela,a6
real,dimension(iMIN-2:iMAX+1) :: one,two,three,four,five
real,dimension(i_min:i_max) :: del,aRd,aLd
real :: test_1,test_2,K_0=0.1
integer :: js,je
js = iMIN-2   ! First cell I interpolate is the i=-2  cell
je = iMAX+1   ! Last cell  I interpolate is the i=Nm  cell
one   = dx(js:je)/(dx(js:je)+dx(js+1:je+1))
two   = 1./(dx(js-1:je-1)+dx(js:je)+dx(js+1:je+1)+dx(js+2:je+2))
three = (2.*dx(js+1:je+1)*dx(js:je))/(dx(js:je)+dx(js+1:je+1))
four  = (dx(js-1:je-1)+dx(js:je))/(2.*dx(js:je)+dx(js+1:je+1))
five  = (dx(js+2:je+2)+dx(js+1:je+1))/(2.*dx(js+1:je+1)+dx(js:je))

CALL avgSlope(a,dx,del)
a12(js:je) = a(js:je) + one*(a(js+1:je+1)-a(js:je)) &
			+ two*(three*(four-five)*(a(js+1:je+1)-a(js:je)) &
            - dx(js:je)*four*del(js+1:je+1)+dx(js+1:je+1)*five*del(js:je))			
aR(iMIN-1:iMAX+1) = a12(iMIN-1:iMAX+1)  ! This is a statment of treating these arrays on level footing ( the whole poinf of the update 7-11)
aL(iMIN-1:iMAX+1) = a12(iMIN-2:iMAX)	! these are the left and right of the interface iMIN
if (detect==1) then
	call discont_detect(a,dx,aR,aRd,aL,aLd,del)
end if 
CALL monotone(a,aL,aR)
dela = aR-aL
a6   = 6.*(a-.5*(aL+aR))			

end subroutine

subroutine interpolOrig(a,dx,a12,aR,aL,dela,a6,detect)
!-------------------------------
! Calculate interpolation value of a_(j+1/2)
! Input:
!   N  : Number of elements in dx
!   af : Solution on extended domain
!   dx : Cell widths on extended domain
! Input/Output:
!   a12    : Interpolated value
!------------------------------
use CommonData
implicit none
integer,intent(in) :: detect
real,intent(in),dimension(i_min:i_max)   :: a,dx
real,intent(inout),dimension(i_min:i_max) :: a12,aR,aL,dela,a6
real,dimension(iMIN-2:iMAX+1) :: one,two,three,four,five
real,dimension(i_min:i_max) :: del,aRd,aLd
real :: test_1,test_2,K_0=0.1
integer :: js,je
js = iMIN-2   ! First cell I interpolate is the i=-2  cell
je = iMAX+1   ! Last cell  I interpolate is the i=Nm  cell
one   = dx(js:je)/(dx(js:je)+dx(js+1:je+1))
two   = 1./(dx(js-1:je-1)+dx(js:je)+dx(js+1:je+1)+dx(js+2:je+2))
three = (2.*dx(js+1:je+1)*dx(js:je))/(dx(js:je)+dx(js+1:je+1))
four  = (dx(js-1:je-1)+dx(js:je))/(2.*dx(js:je)+dx(js+1:je+1))
five  = (dx(js+2:je+2)+dx(js+1:je+1))/(2.*dx(js+1:je+1)+dx(js:je))

CALL avgSlope(a,dx,del)
a12(js:je) = a(js:je) + one*(a(js+1:je+1)-a(js:je)) &
			+ two*(three*(four-five)*(a(js+1:je+1)-a(js:je)) &
            - dx(js:je)*four*del(js+1:je+1)+dx(js+1:je+1)*five*del(js:je))			
aR(iMIN-1:iMAX+1) = a12(iMIN-1:iMAX+1)  ! This is a statment of treating these arrays on level footing ( the whole poinf of the update 7-11)
aL(iMIN-1:iMAX+1) = a12(iMIN-2:iMAX)	! these are the left and right of the interface iMIN
if (detect==1) then
	call discont_detect(a,dx,aR,aRd,aL,aLd,del)
end if 
CALL monotone(a,aL,aR)
dela = aR-aL
a6   = 6.*(a-.5*(aL+aR))			

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

subroutine TimeStep(U,V,dr,dt)
use CommonData 
implicit none
real,intent(in),dimension(i_min:i_max) :: dr
real,intent(in),dimension(3,i_min:i_max) :: U,V
real,intent(inout) :: dt 
real,dimension(i_min:i_max) :: csound,dx
csound = SQRT(gm*V(3,:)/V(1,:))
dt = COURANT*MINVAL(dr)*MAX(MAXVAL(ABS(U(2,:))),MAXVAL(csound))
end subroutine 

subroutine boundaries(dr,U,V,BC)
use CommonData
implicit none
integer,intent(in) :: BC
real,intent(inout),dimension(i_min:i_max) :: dr
real,intent(inout),dimension(3,i_min:i_max) :: U,V
integer i,j,k
!!
!! Working but make it clearer and easier to read, must be better way to do loop 
!!
IF (BC==0) then ! Reflective
	do i=1,ABS(i_min)
		dr(iMIN-i) = dr(iMIN+i-1)
		U(1,iMIN-i)= U(1,iMIN+i-1)
		U(2,iMIN-i)= U(2,iMIN+i-1)
		U(3,iMIN-i)= U(3,iMIN+i-1)
		V(1,iMIN-i)= V(1,iMIN+i-1)
		V(2,iMIN-i)= V(2,iMIN+i-1)
		V(3,iMIN-i)= V(3,iMIN+i-1)
	end do
	do i=iMAX,i_max-1
		dr(i+1) = dr(iMAX-(i-iMAX))
		U(1,i+1)= U(1,iMAX-(i-iMAX))
		U(2,i+1)= U(2,iMAX-(i-iMAX))
		U(3,i+1)= U(3,iMAX-(i-iMAX))
		V(1,i+1)= V(1,iMAX-(i-iMAX))
		V(2,i+1)= V(2,iMAX-(i-iMAX))
		V(3,i+1)= V(3,iMAX-(i-iMAX))
	end do
ELSE IF (BC==1) then ! Periodic
	do i=1,ABS(i_min)
		dr(iMIN-i) = dr(iMAX-i+1)
		U(1,iMIN-i)= U(1,iMAX-i+1)
		U(2,iMIN-i)= U(2,iMAX-i+1)
		U(3,iMIN-i)= U(3,iMAX-i+1)
		V(1,iMIN-i)= V(1,iMAX-i+1)
		V(2,iMIN-i)= V(2,iMAX-i+1)
		V(3,iMIN-i)= V(3,iMAX-i+1)
	end do
	do i=iMAX+1,i_max
		dr(i) = dr(iMIN+i-iMAX-1)
		U(1,i)= U(1,iMIN+i-iMAX-1)
		U(2,i)= U(2,iMIN+i-iMAX-1)
		U(3,i)= U(3,iMIN+i-iMAX-1)
		
		V(1,i)= V(1,iMIN+i-iMAX-1)
		V(2,i)= V(2,iMIN+i-iMAX-1)
		V(3,i)= V(3,iMIN+i-iMAX-1)
	end do
END IF 
end subroutine

subroutine intialCondition(r,U,V,choice)
!--------------------------------
! Creates intial data
! Input:
!   N : Number of elementsi n U
!   r : Space coordinate
! Input/Output:
!   U : Vector of variables (tau,u,E)
!   V : Vector of variables (rho,e,P)
!-------------------------------
use CommonData
implicit none
integer,intent(in) :: choice
real,intent(in),dimension(i_min:i_max) :: r
real,intent(inout),dimension(3,i_min:i_max) :: U,V
integer :: i
if (choice==1) then
  print*,'Sod Test'
  do i=iMIN,iMAX
    if (r(i) < .5) then
      V(1,i) = 1.
      U(2,i) = 0.0
      V(3,i) = 1.
    else
      V(1,i) = 0.125
      U(2,i) = 0.0
      V(3,i) = 0.1
    end if 
  end do
else if (choice == 2) then
	print*,'Constant Move Right'
	V(1,iMIN:iMAX) = 1.
	U(2,iMIN:iMAX) = 1.
	V(3,iMIN:iMAX) = 1.
else if (choice==3) then 
	print*,' SOFT Sod Test'
	do i=iMIN,iMAX
		if (r(i) < .5) then
		V(1,i) = 1.
		U(2,i) = 0.0
		V(3,i) = 1.
    else
		V(1,i) = .9
		U(2,i) = 0.0
		V(3,i) = .9
    end if 
end do
end if 
U(1,iMIN:iMAX) = 1./V(1,iMIN:iMAX)						 ! tau = 1/rho
V(2,iMIN:iMAX) = V(3,iMIN:iMAX)/((gm-1.)*V(1,iMIN:iMAX)) ! e = P/(gm-1)/rho
U(3,iMIN:iMAX) = V(2,iMIN:iMAX) + .5*U(2,iMIN:iMAX)**(2) ! E = e + .5 u**(2)
end subroutine

subroutine Output(U,V,dm,dt,itter,writeOut,printOut)
use CommonData
implicit none
real,intent(in),dimension(3,i_min:i_max) :: U,V
real,intent(in),dimension(i_min:i_max)   :: dm 
real,intent(in) :: dt
integer,intent(in) :: itter,writeOut,printOut
integer :: st,en,choice = 0
if (choice==0) then ! print out real domain
	st = iMIN
	en = iMAX
else if (choice==1) then ! print out full domain 
	st = i_min
	en = i_max
end if 
!----------------------------------------
if (printOut==1) then
	print*,'------------------------------- ITTERATION:',itter,'-------------------------------'
	print*,'Tau:',U(1,st:en)
	print*,'U  :',U(2,st:en)
	print*,'E  :',U(3,st:en)
	print*,'Rho:',V(1,st:en)
	print*,'e  :',V(2,st:en)
	print*,'P  :',V(3,st:en)
	print*,'dt :',dt
	print*,'M  :',SUM(dm(iMIN:iMAX))
end if 
!----------------------------------------
if (writeOut==1) then
	write(1,*) V(1,iMIN:iMAX) ! Rho 
	write(2,*) U(2,iMIN:iMAX) ! Velocity
	write(3,*) V(2,iMIN:iMAX) ! Internal Energy
	write(4,*) V(3,iMIN:iMAX) ! Pressure
	write(5,*) U(3,iMIN:iMAX) ! Total Energy
	write(7,*) dt  	  ! dt
end if 
end subroutine

subroutine linspace(xmin,xmax,x,N)
!---------------------------------
! Assigns array x values, that are 
! linearly spaced.
! Input:
!   xmin - Minimum value
!   xmax - Maximum value
!   N    - Number of elements in x
! Input/Output:
!   x    - Array with N elements 
!-------------------------------
integer,intent(in) :: N
real,intent(in)    :: xmin,xmax
real,intent(inout),dimension(0:N-1) :: x
real    :: dx
integer :: i
dx = (xmax-xmin)/real(N-1)
do i=0,N-1
    x(i) = xmin + i*dx
end do
end subroutine
