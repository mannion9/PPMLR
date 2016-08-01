program RiemannSolver
use CommonData
implicit none
! Primative W = (rho,u,p)^(T)
real,dimension(0:Nm-1) :: x
real,dimension(0:Nt)   :: dt = 0.d0
real :: g,gm1,gp1,p_star,u_star,f,fp,fl,fr,guess,S,e
real :: newton , time , reader
real,dimension(0:2) :: w_l,w_r,w
real,dimension(0:Nm-1) :: rho,u,P,eng
integer :: i,j,test=1,check ! Test 2 and 4 have p_star issue
common/gammas/ g,gm1,gp1 

if (test==1) then
  open(unit=1,file='Exact-Inputs/Input_1.txt')
  print*,'Test 1'
else if (test==2) then
  open(unit=1,file='Exact-Inputs/Input_2.txt')
  print*,'Test 2'
end if 
open(unit=2,file='Output/dt.txt')


do i=0,2
  read(1,*) w_l(i)
end do
do i=0,2
  read(1,*) w_r(i)
end do
do i=0,Nt
	read(2,*,IOSTAT=check) reader
	if (check>0 .OR. check<0) then ! This statment will stop the loop if we have reached the end of the file
		exit
	else 
		dt(i) = reader
	end if 
end do 
close(1)
close(2)

g = 1.4        ! C_v/C_P = 1.4 (air)
gm1 = g-1.
gp1 = g+1.
! This is a guess for p*
guess = .5*(w_l(2)+w_r(2))

p_star = newton(guess,w_l,w_r)
CALL f_func(p_star,w_l,w_r,f,fp,fl,fr)
u_star = .5*(w_l(1)+w_r(1))+.5*(fr-fl)
! print*,'-----------------------------------------------------------'
! print*,'   P_star           U_star     '
! print*,p_star,u_star
CALL linspace(0.,1.,Nm,x)
open(unit=1,file='Exact-Output/domain.txt')
write(1,*) x
close(1)
open(unit=1,file='Exact-Output/rho.txt')
open(unit=2,file='Exact-Output/velocity.txt')
open(unit=3,file='Exact-Output/pressure.txt')
open(unit=4,file='Exact-Output/energy.txt')
! Time equal zero
do i=0,Nm-1
	if (x(i).LE.0.5) then
		rho(i) = w_l(0)
		u(i)   = w_l(1)
		P(i)   = w_l(2)
		eng(i) = w_l(2)/(gm1*w_l(0))
	else 
		rho(i) = w_r(0)
		u(i)   = w_r(1)
		P(i)   = w_r(2)
		eng(i) = w_r(2)/(gm1*w_r(0))
	end if 
end do 
write(1,*) rho
write(2,*) u
write(3,*) P
write(4,*) eng
! Future time
time = 0.
j = 0
do while(time.LE.t_max)
	do i=0,Nm-1
		!S = (x(i)-.5)/(j*dt(j))
		S = (x(i)-0.5)/time
		CALL sample(S,u_star,p_star,w_l,w_r,w)
		CALL internal_energy(w,e)
		rho(i) = w(0)
		u(i)   = w(1)
		P(i)   = w(2)
		eng(i) = e
	end do 
	time = time + dt(j)
	j = j + 1
	write(1,*) rho
	write(2,*) u 
	write(3,*) P
	write(4,*) eng
end do
close(1)
close(2)
close(3)
close(4)

end program

subroutine internal_energy(w,e)
!--------------------------------------------
! Determies the interal energy of an ideal gas
! INPUTS:
!   w : primative variables
! OUTPUTS:
!   e : internal energy
!-------------------------------------------
implicit none
real,intent(in),dimension(0:2) :: w
real,intent(inout) :: e
real :: g,gm1,gp1
common/gammas/g,gm1,gp1
e = w(2)/(gm1*w(0))  ! Toro Eq. 4.3
end subroutine

subroutine sample(S,U_star,P_star,w_l,w_r,w)
!-------------------------------------------------
! Sample determines for the value of the primates
! in each at a specifc location
! INPUTS: 
!   S     : x/t is the location in x-t space
!   U_star: U in starred region
!   P_star: P in starred region
!   w_L   : primative variables left
!   w_R   : primative variables right
! OUTPUTS:
!   w     : primative variables at x/t
!-------------------------------------------------

implicit none
real,intent(in) :: S,U_star,P_star
real,intent(in),dimension(0:2) :: w_l,w_r
real,intent(inout),dimension(0:2) :: w
real :: a_l,a_r, a_l_star,a_r_star, &
        S_l,S_l_tl,S_l_hd, &
        S_r,S_r_tl,S_r_hd, &
        rho_l_star_shock,rho_r_star_shock, &
        rho_l_star_rare,rho_r_star_rare, &
        rho_l_fan,rho_r_fan, &
        u_l_fan,u_r_fan,p_l_fan,p_r_fan
real,dimension(0:2) :: w_l_star,w_r_star,w_l_fan_star,w_r_fan_star, w_r_fan,w_l_fan
real :: g,gm1,gp1
common/gammas/g,gm1,gp1

a_l = SQRT(g*w_l(2)/w_l(0))   ! Toro 3.6 
a_r = SQRT(g*w_r(2)/w_r(0))  
a_l_star = a_l*(p_star/w_l(2))**(gm1/(2.*g))  ! Toro 4.54/4.61 
a_r_star = a_r*(p_star/w_r(2))**(gm1/(2.*g))
rho_l_star_shock = w_l(0)*(((P_star/w_l(2))+gm1/gp1)/(gm1*P_star/(gp1*w_l(2))+1)) ! Toro 4.50/4.57 
rho_r_star_shock = w_r(0)*(((P_star/w_r(2))+gm1/gp1)/(gm1*P_star/(gp1*w_r(2))+1)) 
rho_l_star_rare  = w_l(0)*(P_star/w_l(2))**(1./g) ! Toro 4.53/4.60
rho_r_star_rare  = w_r(0)*(P_star/w_r(2))**(1./g)
rho_l_fan = w_l(0)*(2./gp1+gm1/(a_l*gp1)*(w_l(1)-S))**(2./gm1) ! Toro 4.56/4.63
rho_r_fan = w_r(0)*(2./gp1-gm1/(a_r*gp1)*(w_r(1)-S))**(2./gm1)
u_l_fan = (2./gp1)*(a_l+(gm1/2.)*w_l(1)+S)   ! Toro 4.56/4.63
u_r_fan = (2./gp1)*(-a_r+(gm1/2.)*w_r(1)+S)
p_l_fan = w_l(2)*(2./gp1+gm1/(gp1*a_l)*(w_l(1)-S) )**(2.*g/gm1)
p_r_fan = w_r(2)*(2./gp1-gm1/(gp1*a_r)*(w_r(1)-S) )**(2.*g/gm1)
S_l = w_l(1) - a_l*SQRT(gp1*P_star/(2.*g*w_l(2))+(gm1/(2.*g))) ! Toro Eq. 4.52
S_r = w_r(1) + a_r*SQRT(gp1*P_star/(2.*g*w_r(2))+(gm1/(2.*g))) ! Toro Eq. 4.52
S_l_hd = w_l(1) - a_l       ! Toro 4.55
S_l_tl = U_star - a_l_star
S_r_hd = w_r(1) + a_r       ! Toro 4.62
S_r_tl = U_star + a_r_star
! Shock
w_l_star = (/rho_l_star_shock,U_star,P_star/)
! Rarefaction (in starred region)
w_l_fan_star = (/rho_l_star_rare,U_star,P_star/)
! Rarefaction (within fan)
w_l_fan = (/rho_l_fan,u_l_fan,p_l_fan/)
! Shock
w_r_star = (/rho_r_star_shock,U_star,P_star/)
! Rarefaction (in starred region)
w_r_fan_star = (/rho_r_star_rare,U_star,P_star/)
! Rarefaction (in fan)
w_r_fan = (/rho_r_fan,u_r_fan,p_r_fan/)

if (S.LE.U_star) then
  ! Left of discontinuity
  if (P_star>w_l(2)) then
    ! Shock
    if (S.LE.S_l) then
      ! To the left of the shock
      w = w_l  ! Will have condition from left
    else
      ! To the right of the shock
      w= w_l_star  ! Will have conditions in star region
    end if 
  else 
    ! Rarefaction
    if (S.LE.S_l_hd) then
      ! To the right of head
      w = w_l
    elseif (S.GE.S_l_hd .AND. S.LE.S_l_tl) then
     ! Inside the fan
      w = w_l_fan
    else 
     ! In the star region
      w = w_l_fan_star
    end if 
  end if 
else if (S.GE.U_star) then
  ! Right of disconinuity
  if (P_star>w_r(2)) then
    ! Shock 
    if (S.GE.S_r) then
      ! To the right of the shock
      w = w_r
    else 
      ! To the left of the shock
      w = w_r_star
    end if 
  else
    ! Rarefaction
    if (S.GE.S_r_tl .AND. S.LE.S_r_hd) then
      ! Inside the fan
      w = w_r_fan
    else if (S.GE.S_r_hd) then
      ! To the right of head
      w = w_r
    else
      ! In the star region
      w = w_r_fan_star
    end if
  end if 
end if 
end subroutine

function newton(guess,w_l,w_r) result(p_star)
!---------------------------------------------
! Newton determines the root of an equation
! using newtons method. 
!
! INPUTS:
!   guess : intial guess of value of root
!   w_l   : primative variables left     (0)
!   w_l   : primative variables right    (0)
!   f_func: the function we seek root of
! OUTPUT
!   p_star: value of p in stared region
!--------------------------------------------
implicit none
real :: guess,p_star
real,dimension(0:2) :: w_l,w_r
real :: f,fp,fl,fr,it, itold, TOL=10.**(-6) , CHA
integer :: i, itter_max=20,break=0
it = guess
CHA =  guess 
do i=1,itter_max 
  itold = it
  if (CHA > TOL) then !.AND. i<itter_max) then
    CALL f_func(it,w_l,w_r,f,fp,fl,fr)
    it = itold - f/fp
    if (it < 0) then  ! Requires that the solution be positive for physical
      it = itold
    end if 
    CHA = 2.*abs(it-itold)/(it+itold)
  else 
    exit 
  end if 
end do
p_star = it
end  function newton 

subroutine f_func(p,w_l,w_r,f,fp,f_l,f_r) 
!---------------------------------------
! The function, whose root returns p* 
! in the star region.
! 
! INPUTS: 
!   p  : pressure in sarted region (0)
!   w_l: primative variables left  (0)
!   w_r: primative variables right (0:2)
!   recall - W = (rho,u,p)^(T)
! OUTPUT:
!   f  : scalar (0)
!   fp : derivative of f
!   fl : function fl
!   fr : function fr
!---------------------------------------
implicit none
real,intent(in),dimension(0:2) :: w_l,w_r
real,intent(in) :: p
real,intent(inout) :: f,fp,f_l,f_r
real :: f_lp,f_rp,del_u,A_l,A_r,B_l,B_r,al,ar
real :: A,B,ak,rho,pk,fk,fkp
integer:: i
real :: g,gm1,gp1
common/gammas/ g,gm1,gp1
del_u  = w_r(1)-w_l(1)
A_l    = (2./gp1)/w_l(0)   ! Toro 4.8     
A_r    = (2./gp1)/w_r(0)                  
B_l    = (gm1/gp1)*w_l(2)                 
B_r    = (gm1/gp1)*w_r(2)                 
al     = SQRT(g*w_l(2)/w_l(0))   ! Toro 3.6 
ar     = SQRT(g*w_r(2)/w_r(0))            
do i=1,2
  if (i==1) then
    A    = A_l
    B    = B_l
    ak    = al
    rho  = w_l(0)
    pk   = w_l(2)
  else
    A    = A_r
    B    = B_r
    ak    = a_r
    rho = w_r(0)
    pk   = w_r(2)
  end if 
  if (p>pk) then
    ! shock
      fk  = (p-pk)*SQRT(A/(p+B))
      fkp = SQRT(A/(B+p))*(1-(p-pk)/(2.*(B+p)))
  else
    ! rarefaction
      fk  = (2.*ak)/gm1*((p/pk)**(gm1/(2.*g))-1.)
      fkp = (1./(rho*ak))*(p/pk)**(-1.*gm1/(2.*g))
  end if
  if (i==1) then
    f_l  = fk
    f_lp = fkp
  else
    f_r  = fk
    f_rp = fkp
  end if
end do
f   = f_l + f_r + del_u
fp  = f_lp + f_rp
end subroutine

subroutine linspace(x_min,x_max,N,x)
!-------------------------------------------
! Creates a linear spaced vector of lenght N
! with indicies 0 to  N-1, with values 
! ranging from x_min to x_max.
!-------------------------------------------
implicit none
real,intent(in) :: x_min , x_max
integer,intent(in) :: N
real,intent(inout),dimension(0:N-1) :: x
integer :: i
real :: dx
dx = (x_max-x_min)/real(N-1)
do i=0,N-1
    x(i) = x_min + i*dx
  end do
end subroutine



