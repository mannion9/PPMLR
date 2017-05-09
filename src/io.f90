subroutine Output(d,e,h,p,u,m,time,r,rc,itter,gm)
use Constants
implicit none
real,intent(in),dimension(imin:imax)   :: d,e,u,h,p,m,rc,gm
real,intent(in),dimension(imin:imax+1) :: r
integer,intent(in) :: itter
real,intent(in) :: time
write(1,*)  d(jmin:jmax)  ! Density
write(2,*)  u(jmin:jmax)  ! Velocity
write(3,*)  e(jmin:jmax)  ! Internal Energy
write(4,*)  p(jmin:jmax)  ! Pressure
write(7,*)  r(jmin:jmax+1)! Lagrange/Remapped Cell Edges
write(8,*)  time          ! Current Time
write(9,*)  rc(jmin:jmax) ! Lagrange/Remapped Cell Ceters
write(10,*) SUM(m)        ! Total mass
write(11,*) itter         ! Itteration number
write(12,*) h(jmin:jmax)  ! Total energy
write(14,*) gm(jmin:jmax) ! Gamma
end subroutine

subroutine OpenFiles()
implicit none
! Open Files to write out to
open(unit=1, file='Output/Density.txt'       ,status="replace",action="write")
open(unit=2, file='Output/Velocity.txt'      ,status="replace",action="write")
open(unit=3, file='Output/InternalEnergy.txt',status="replace",action="write")
open(unit=4, file='Output/Pressure.txt'      ,status="replace",action="write")
open(unit=7, file='Output/LagnCellEdge.txt'  ,status="replace",action="write")
open(unit=8, file='Output/CurrentTime.txt'   ,status="replace",action="write")
open(unit=9, file='Output/LagnCellCenter.txt',status="replace",action="write")
open(unit=10,file='Output/TotalMass.txt'     ,status="replace",action="write")
open(unit=11,file='Output/Step.txt'          ,status="replace",action="write")
open(unit=12,file='Output/TotalEnergy.txt'   ,status="replace",action="write")
open(unit=13,file='Output/InitialState.txt'  ,status="replace",action="write")
open(unit=14,file='Output/Gamma.txt'         ,status="replace",action="write")
open(unit=15,file='Output/Domain.txt'        ,status="replace",action="write")
end subroutine

subroutine CloseFiles()
implicit none
! Open Files to write out to
close(1) ; close(2) ; close(3) ; close(4) ; close(7) ; close(8)
close(9) ; close(10); close(11); close(12); close(13); close(14)
close(15)
end subroutine
