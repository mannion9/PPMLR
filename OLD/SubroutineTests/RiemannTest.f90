module CommonData
implicit none 
real,parameter :: gm = 1.4,delta=1.-30 , TOL=1.0d-5
end module 

program main
use CommonData
implicit none
real :: c_p,c_m,p_star,u_star
real :: rhoL=1.,uL=0.,pL=1.,rhoR=0.125,uR=0.,pR=.1  ! Sod Test
!real :: rhoL=1.,uL=1.,pL=1.,rhoR=1.   ,uR=1.,pR=1.


c_p=(gm*pL*rhoL) ! Square of speed of sound
c_m=(gm*pR*rhoR) ! Square of speed of sound

print*,'rhoL                     uL              pL              rhoR              uR         pR'
print*,rhoL,uL,pL,rhoR,uR,pR
CALL RiemannSolver(pL,pR,uL,uR,c_p,c_m,p_star,u_star)
print*,'SoundL:',c_p
print*,'SoundR:',c_m
print*,'p_star:',p_star
print*,'u_star:',u_star

end program main

subroutine RiemannSolver(pL,pR,uL,uR,c_p,c_m,p_star,u_star)
use CommonData
implicit none
real,intent(in)    :: pL,pR,uL,uR,c_m,c_p
real,intent(inout) :: p_star,u_star
real:: diff,Temp,Wl,del_u,f,f_var,fprime,gmp1=gm+1
integer :: j,i,itterMax=20
p_star  = .5*(pL+pR)
! Loop in Main.f90 starts here
diff = 1.
del_u = uR - uL
do i=1,itterMax
	if (diff.GT.TOL) then
		Temp = p_star
		call Func_RP(c_p,c_m,p_star,pL,pR,del_u,f)	! f(x)
		call Func_RP(c_p,c_m,(p_star+TOL),pL,pR,del_u,f_var)	! f(x+dx)
		fprime = (f_var-f)/TOL									! f'(x) ! (f(x)-f(x+dx))/dx
		if (fprime.LE.delta .OR. fprime.NE.fprime .OR. p_star.LT.delta) then ! ensures no divide by zero or NAN or negative pressure
			p_star = Temp
			exit
		end if 
		p_star = p_star - f/fprime  ! Newtons method itteration
		diff   = 2.*ABS((p_star-Temp)/(p_star+Temp))
	else 
		exit
	end if 
end do
print*,'Itterations:',i
Wl = SQRT(ABS(c_p*(1+(gmp1/(2.*gm))*((p_star/pL)-1.))))
u_star = -(p_star-pL)/Wl + uL
end subroutine

subroutine Func_RP(cp,cm,p,pl,pr,del_u,eval)
use CommonData
implicit none
real,intent(in)    :: cp,cm,p,pl,pr,del_u
real,intent(inout) :: eval
real :: Wl,Wr
Wl   = SQRT(ABS(cp*(1.+((gm+1)/(2.*gm))*((p/Pl)-1.))))
Wr   = SQRT(ABS(cm*(1.+((gm+1)/(2.*gm))*((p/Pr)-1.))))
eval = (p-pL)/Wl + (p-pR)/Wr - del_u
end subroutine


