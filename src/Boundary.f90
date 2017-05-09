subroutine Boundary(d,u,p)
use Constants
implicit none
real,intent(inout),dimension(imin:imax) :: d,p,u


if (BC.EQ.'Free') then
	u(imin:jmin-1) = u(jmin+2:jmin:-1)
	u(jmax+1:imax) = u(jmax:jmax-3:-1)
else if (BC.EQ.'Reflective') then
	u(imin:jmin-1) = -u(jmin+2:jmin:-1)
	u(jmax+1:imax) = -u(jmax:jmax-3:-1)
end if

if (IC.EQ.0) u(jmax+1:imax) = -100. !!u(imin:jmin-1)=100.

p(imin:jmin-1) = p(jmin+2:jmin:-1)
p(jmax+1:imax) = p(jmax:jmax-3:-1)
d(imin:jmin-1) = d(jmin+2:jmin:-1)
d(jmax+1:imax) = d(jmax:jmax-3:-1)
end subroutine
