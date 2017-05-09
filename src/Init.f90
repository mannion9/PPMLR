subroutine InitDomain(r,dv,rc,r0,dv0,rc0)
use Constants
implicit none
real,intent(inout),dimension(imin:imax+1) :: r,r0
real,intent(inout),dimension(imin:imax)  :: dv,dv0,rc,rc0
integer :: i
dv = (rmax-rmin)/REAL(N-1)	! Cell widths
r(imin:imax+1) = (/(rmin+dv(jmin)*(real(i-jmin)),i=imin,imax+1,1)/)
rc  = .5*(r(imin:imax)+r(imin+1:imax+1))				! The i-jmin allows us to start couting at any number
dv(imin:imax) = abs((/(r(i+1)**ap1-r(i)**ap1,i=imin,imax,1)/))/REAL(ap1)	! Volume of cell
r0 = r ; dv0= dv ; rc0 = rc	 									! Stationary cell edges wdith and center
end subroutine

subroutine InitCond(d,e,h,u,m,dv,rc)
use Constants
implicit none
real,intent(inout),dimension(imin:imax)   :: d,e,u,h,m,dv,rc
real,dimension(imin:imax) :: p
real :: dL,dR,uL,uR,pL,pR,x0=rmin+(rmax-rmin)/2.
integer :: i

if (IC.EQ.0) then		! Piston problem
	dL = 1.    ; dR = dL
	pL = 1.    ; pR = pL
	uL = 0.0   ; uR = uL
else if (IC.EQ.1) then  ! Modified sod test
	dL = 1.    ; dR = dL/8.
	pL = 1.    ; pR = pL/10.
	pL = 383703079320.175   ; pR = 38366184029.3157
	uL = 0.0 ; uR = 0.0
else if (IC.EQ.2) then  ! Strong rarefaction waves
	dL = 1.    ; dR = dL
	pL = 0.4   ; pR = pL
	uL = -2.   ; uR = -uL
else if (IC.EQ.3) then
	dL = 1.    ; dR = dL
	pL = 1000. ; pR = 0.01
	uL = 0.0   ; uR = 0.0
else if (IC.EQ.4) then
	dL = 5.99924
	dR = 5.99242
	pL = 460.894
	pR = 46.0950
	uL = 19.5975
	uR = -6.19633
else if (IC.EQ.5) then
	dL = 1.  ; dR = dL
	pL = 1000.
	pR = 0.01
	uL = -19.59745
	uR = uL
else if (IC.EQ.6) then  ! Strong compression inward
	dL = 1.4    ; dR = dL
	pL = .7143  ; pR = pL
	uL = 5.0    ; uR = -uL
end if

do i = imin,imax
	if (rc(i).LT.x0) then
		d(i) = dL ; p(i) = pL ; u(i) = uL
	else
		d(i) = dR ; p(i) = pR ; u(i) = uR
	end if
end do

! Sedov Blast
if (IC.EQ.7) then
	d = 1.
	u   = delta
	p(jmin) = 3.*(5./3.-1.)*d(jmin)
	p(jmin+1:imax) = delta
	! p(imin:INT(N/2)-1) = delta
	! p(INT(N/2))        = 3*(5./3.-1.)*d(3)  ! Give internal energy of .3 to orgin cell (this just tells us the pressure, given .3 internal energy)
	! p(INT(N/2)+1:imax) = delta
end if

! Collela-Woodward Blast
if (IC.EQ.8) then
	d = 1.0
	u = 0.0
	do i=imin,imax
		if (rc(i).LT.0.1*(rmax-rmin)) then
			p(i)   = 1000.
		else if (rc(i).GE.0.1*(rmax-rmin).AND.rc(i).LT.0.9*(rmax-rmin)) then   ! Left state
			p(i) = 0.01
		else  ! Right state
			p(i) = 100.0
		end if
	end do
end if

if (IC.EQ.9) then
	d = 1. ; p = 1. ; u = 1.
end if

! Assumes ideal gas at initial state (makes intial conditions easier to match with others)
e  = p/((5./3-1.)*d)
h  = e + .5*(u)**2
m  = d*dv

! Write this file for Riemann Solver
write(13,*) dL   ; write(13,*) dR
write(13,*) pL   ; write(13,*) pR
write(13,*) uL   ; write(13,*) uR
write(13,*) x0   ; write(15,*) N
write(15,*) rmax ; write(15,*) rmin
close(13)        ; close(15)
end subroutine
