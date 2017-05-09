! The length of the time array and number of spatial cells is "unkown" for this program since this depends on the MAIN
! programs result and t_max. We can bound the time array by itterMAX from the MAIN program.
! To work around this we just make the time array itterMAX long. We make the primative variable arrays have itterMAX length
! aswell. If the MAIN program has more spatial cells than itterMAX and error is thrown.
! We will count from jmin=0,jmax=N-1 where N is the number of grid cells.

program RiemannSolver
use constants
implicit none
real :: ps,us,f,fp,fl,fr,guess,S,center,aR,aL
real,dimension(0:itterMax) :: x,d,u,p,e
real,dimension(0:itterMax) :: time = 0.d0
real :: dl,pl,ul,dr,pr,ur,xmin,xmax
integer :: i,j,Ntmax,N,jmin,jmax

call OpenFiles()
call ReadInputs(dl,dr,pl,pr,ul,ur,center,time,Ntmax,N,xmin,xmax)
call Init(x,jmin,jmax,N,xmin,xmax,d,u,p,e,center,dl,dr,ul,ur,pl,pr)

! Find ps and us
guess = .5*(pl+pr)
aL = SQRT(gamma*pl/dl)
aR = SQRT(gamma*pr/dR)
!guess = .5*(pl+pr)-.125*(ur-ul)*(dl+dr)*(aL+aR)

call Pstar(guess,dl,pl,ul,dr,pr,ur,ps)
CALL f_func(ps,dl,pl,ul,dr,pr,ur,f,fp,fl,fr)
us = .5*(ul+ur)+.5*(fr-fl)

call Output(x,d,u,p,e,jmin,jmax)

do j=1,Ntmax-1
	do i=jmin,jmax
		S = (x(i)-center)/time(j)
		CALL sample(S,us,ps,dl,pl,ul,dr,pr,ur,d(i),u(i),p(i))
		e(i) = p(i)/(gm1*d(i))
	end do
	call Output(x,d,u,p,e,jmin,jmax)
end do

call CloseFiles()
end program

subroutine sample(S,us,ps,dl,pl,ul,dr,pr,ur,d,u,p)
use constants
implicit none
real,intent(in)    :: S,us,ps,dl,pl,ul,dr,pr,ur
real,intent(inout) :: d,u,p
real :: aL,aR, aLs,aRs,SL,SLtl,SLhd,SR,SRtl,S_r_hd, &
		dLss,dRss,dLsr,dRsr,dLfan,dRfan,uLf,uRf,pLf,pRf

aL	  = SQRT(gamma*pl/dl)   								! Toro Eq. 3.6
aR    = SQRT(gamma*pr/dr)
aLs   = aL*(ps/pl)**(gm1/(2.*gamma))  						! Toro Eq. 4.54/4.61
aRs   = aR*(ps/pr)**(gm1/(2.*gamma))
dLss  = dl*(((ps/pl)+gm1/gp1)/(gm1*ps/(gp1*pl)+1)) 			! Toro Eq. 4.50/4.57
dRss  = dr*(((ps/pr)+gm1/gp1)/(gm1*ps/(gp1*pr)+1))
dLsr  = dl*(ps/pl)**(1./gamma)					      		! Toro 4.53/4.60
dRsr  = dr*(ps/pr)**(1./gamma)
dLfan = dl*(2./gp1+gm1/(aL*gp1)*(ul-S))**(2./gm1) 		    ! Toro 4.56/4.63
dRfan = dr*(2./gp1-gm1/(aR*gp1)*(ur-S))**(2./gm1)
uLf   = (2./gp1)*(aL+(gm1/2.)*ul+S)   						! Toro 4.56/4.63
uRf   = (2./gp1)*(-aR+(gm1/2.)*ur+S)
pLf   = pl*(2./gp1+gm1/(gp1*aL)*(ul-S) )**(2.*gamma/gm1)
pRf   = pr*(2./gp1-gm1/(gp1*aR)*(ur-S) )**(2.*gamma/gm1)
SL    = ul - aL*SQRT(gp1*ps/(2.*gamma*pl)+(gm1/(2.*gamma))) ! Toro Eq. 4.52
SR    = ur + aR*SQRT(gp1*ps/(2.*gamma*pr)+(gm1/(2.*gamma))) ! Toro Eq. 4.52
SLhd  = ul - aL      										! Toro 4.55
SLtl  = us - aLs
S_r_hd= ur + aR       										! Toro 4.62
SRtl  = us + aRs

if (S.LE.us) then 	   			! Left of discontinuity
	if (ps.GT.pl) then 			! Shock
		if (S.LE.SL) then 		! To the left of the shock
			d = dl ; u = ul ; p = pl
		else			   		! To the right of the shock
			d = dLss ; u = us ; p = ps
		end if
	else				     	 ! Rarefaction
		if (S.LE.SLhd) then 	 ! To the right of head
			d = dl ; u = ul ; p = pl
		else if (S.GE.SLhd .AND. S.LE.SLtl) then	! Inside the fan
			d = dLfan ; u = uLf ; p = pLf
		else										! In the star region
			d = dLsr ; u = us ;p = ps
		end if
	end if
else if (S.GE.us) then		! Right of disconinuity
	if (ps.GT.pr) then		! Shock
		if (S.GE.SR) then	! To the right of the shock
			d = dr ; u = ur ; p = pr
		else				! To the left of the shock
			d = dRss ; u = us ; p = ps
		end if
	else										! Rarefaction
		if (S.GE.SRtl .AND. S.LE.S_r_hd) then	! Inside the fan
			d = dRfan ; u = uRf ; p = pRf
		else if (S.GE.S_r_hd) then				! To the right of head
			d = dr ; u = ur ; p = pr
		else									! In the star region
			d = dRsr ; u = us ; p = ps
		end if
	end if
end if
end subroutine

subroutine Pstar(guess,dl,pl,ul,dr,pr,ur,ps)
implicit none
real :: guess,ps
real :: dl,pl,ul,dr,pr,ur
real :: f,fp,fl,fr,it, itold, TOL=10.**(-6) , CHA
integer :: i, itter_max=20
it = guess
CHA =  guess
do i=1,itter_max
	itold = it
	if (CHA > TOL) then !.AND. i<itter_max) then
		CALL f_func(it,dl,pl,ul,dr,pr,ur,f,fp,fl,fr)
		it = itold - f/fp
		if (it < 0) then  ! Requires that the solution be positive for physical
			it = itold
		end if
		CHA = 2.*abs(it-itold)/(it+itold)
	else
		exit
	end if
end do
ps = it
end

subroutine f_func(p,dl,pl,ul,dr,pr,ur,f,fp,f_l,f_r)
use constants
implicit none
real,intent(in):: dl,pl,ul,dr,pr,ur
real,intent(in) :: p
real,intent(inout) :: f,fp,f_l,f_r
real :: f_lp,f_rp,del_u,A_l,A_r,B_l,B_r,al,ar
real :: A,B,ak,rho,pk,fk,fkp
integer:: i
del_u  = ur-ul
A_l    = (2./gp1)/dl   	    ! Toro 4.8
A_r    = (2./gp1)/dr
B_l    = (gm1/gp1)*pl
B_r    = (gm1/gp1)*pr
al     = SQRT(gamma*pl/dl)   ! Toro 3.6
ar     = SQRT(gamma*pr/dr)
do i=1,2
  if (i==1) then
	A    = A_l
	B    = B_l
	ak    = al
	rho  = dl
	pk   = pl
  else
	A   = A_r
	B   = B_r
	ak  = a_r
	rho = dr
	pk  = pr
  end if
  if (p>pk) then
	! shock
	  fk  = (p-pk)*SQRT(A/(p+B))
	  fkp = SQRT(A/(B+p))*(1-(p-pk)/(2.*(B+p)))
  else
	! rarefaction
	  fk  = (2.*ak)/gm1*((p/pk)**(gm1/(2.*gamma))-1.)
	  fkp = (1./(rho*ak))*(p/pk)**(-1.*gm1/(2.*gamma))
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
