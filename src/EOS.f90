subroutine EOS_ideal_gas(d,e,p,a,gm,igm)
! Will give you pressure and lagrange speed of sound
use constants
real,intent(in)   ,dimension(imin:imax) :: d,e
real,intent(inout),dimension(imin:imax) :: p,a,gm,igm
real :: gamma
gamma = 5./3.
p = (gamma-1.)*d*e
gm   = p/(d*e)+1.
a = SQRT(gm*p/d)     ! speed of sound (Eulerian)
igm  = a**(2)*d/p
end subroutine
