module CommonData
implicit none 
!real,parameter :: rhoL=1.,uL=0.,pL=1.,rhoR=0.125,uR=0.,pR=.1,gm=1.4
real,parameter :: rhoL=1.,uL=0.,pL=1.,rhoR=0.125,uR=0.,pR=1.,gm=1.4, &
					delta=1.-30
real :: cm=SQRT(gm*pL*rhoL),cp=SQRT(gm*pR*rhoR)
end module 

program main
use CommonData
implicit none
real :: newton
real :: guess
real :: p_star,u_star,rho_starL,rho_starR

print*,'cp:',cp
print*,'cm:',cm

guess  = .5*(pL+pR)
print*,'pL:',pL
print*,'pR:',pR
p_star = newton(guess)
call starred(p_star,u_star) 
print*,'p_star:',p_star
print*,'u_star:',u_star

end program main

function newton(guess) result(sol)
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
use CommonData 
implicit none
real :: guess,sol
real :: f,fvar,fp,sol_m1,sol_var,TOL=10.**(-6) , diff 
integer :: i, itter_max=20
sol = guess
diff =  guess 
print*,'GUESS:',guess
do i=1,itter_max 
  sol_m1 = sol
  if (diff > TOL) then !.AND. i<itter_max) then
    CALL f_func(sol,f)
    sol_var = sol + TOL
    CALL f_func(sol_var,fvar)
	fp = (f-fvar)/TOL
    sol = sol_m1 - f/fp
	if (fp.LE.delta .OR. fp.NE.fp .OR. sol.LT.delta) then ! ensures no divide by zero or NAN or negative for physicsal meaning
		sol = sol_m1
		exit
	end if 
  end if
  end do
end  function newton 

subroutine f_func(p_star,eval)
use CommonData
implicit none
real,intent(inout) :: p_star,eval
real :: Wl,Wr,del_u
del_u = uR-uL

Wl   = SQRT(ABS(cp*(1.+((gm+1.)/(2.*gm))*((p_star/Pl)-1.))))
Wr   = SQRT(ABS(cm*(1.+((gm+1.)/(2.*gm))*((p_star/Pr)-1.))))
eval = (p_star-pL)/Wl + (p_star-pR)/Wr - del_u
end subroutine

subroutine starred(p_star,u_star)
use CommonData
implicit none
real,intent(in) :: p_star
real,intent(inout) :: u_star
real :: Wl
Wl   = SQRT(ABS(cp*(1.+((gm+1.)/(2.*gm))*((p_star/Pl)-1.))))
u_star =  uL-(p_star-pL)/Wl 
end subroutine
