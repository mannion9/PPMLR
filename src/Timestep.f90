subroutine Timestep(u,r,dt,a)
use Constants
implicit none
real,intent(inout) :: dt
real,intent(in),dimension(imin:imax)  :: u,a
real,intent(in),dimension(imin:imax+1):: r
real,dimension(imin:imax) :: dx
real :: dt0,dtLag,dtRemap
dt0 = dt                             ! Previous time step
dx  = r(imin+1:imax+1)-r(imin:imax)  ! Width of cells
if (Eulerian.EQ.'Off') then          ! Pure Lagrange Calculation
  dt = COURANT*MINVAL(dx/(a+abs(u)))
else
  dtRemap = MIN(MINVAL(dx/MAXVAL(ABS(u))),MINVAL(dx/(a+ABS(u))))
  dtLag   = MINVAL(dx/(a+abs(u)))
  dt       = COURANT*MIN(dtLag,dtRemap)
end if
dt = MIN(dt,1.1*dt0)
end subroutine
