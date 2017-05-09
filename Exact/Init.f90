subroutine Init(x,jmin,jmax,N,xmin,xmax,d,u,p,e,center,dl,dr,ul,ur,pl,pr)
use constants
implicit none
real,intent(inout),dimension(0:itterMax) :: x,d,u,p,e
real,intent(inout) 	  :: xmin,xmax,center,dl,dr,ul,ur,pl,pr
integer,intent(inout) :: N,jmin,jmax
integer :: i
real 	:: dx
! Construct cell centers
jmin = 0 ; jmax = N-1
dx   = (xmax-xmin)/REAL(N-1)
x(jmin:jmax+1) = (/(xmin+real(i-jmin)*dx,i=jmin,jmax+1,1)/)
x(jmin:jmax) = .5*(x(jmin:jmax)+x(jmin+1:jmax+1))
do i=jmin,jmax
	if (x(i).LE.center) then
		d(i) = dl ; u(i) = ul ; p(i) = pl
		e(i) = P(i)/(gm1*d(i))
	else
		d(i) = dr ; u(i) = ur ; p(i) = pr
		e(i) = P(i)/(gm1*d(i))
	end if
end do
end subroutine
