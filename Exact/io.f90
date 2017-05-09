subroutine ReadInputs(dl,dr,pl,pr,ul,ur,center,time,Ntmax,N,xmin,xmax)
use constants
implicit none
real,intent(inout),dimension(0:itterMax) :: time
integer,intent(inout) :: Ntmax,N
real,intent(inout) 	  :: center,xmin,xmax
real 	:: dl,pl,ul,dr,pr,ur
integer :: i,check

! Read in left and right state and discontinity location
read(1,*) dl ; read(1,*) dr
read(1,*) pl ; read(1,*) pr
read(1,*) ul ; read(1,*) ur
read(1,*) center

! Read in time steps
do i=0,itterMAX
	read(2,*,IOSTAT=check) time(i)
	if (check>0 .OR. check<0) exit
end do
Ntmax = MAXLOC(time,1)

! Read in cell number and limits
read(3,*) N ; read(3,*) xmax ; read(3,*) xmin

if (itterMax.LT.N) print*,'WARNING!:', 'Exact riemann solver will not have large enough arrays.', &
						  'Increase itterMax in Exact/Constants.f90', &
						  ', shut off exact solver, or decrease number of spacial zones.'
end subroutine

subroutine Output(x,d,u,p,e,jmin,jmax)
use constants
implicit none
real,intent(in),dimension(0:itterMax) :: x,d,u,p,e
integer,intent(in) :: jmin,jmax
write(4,*)  d(jmin:jmax)
write(7,*)  u(jmin:jmax)
write(8,*)  P(jmin:jmax)
write(9,*)  e(jmin:jmax)
write(10,*) x(jmin:jmax)
end subroutine

subroutine OpenFiles()
implicit none
! Open input files
open(unit=1,file='Inputs/InitialState.txt',status='old',action="read")
open(unit=2,file='Inputs/CurrentTime.txt' ,status='old',action="read")
open(unit=3,file='Inputs/Domain.txt'      ,status='old',action="read")
open(unit=4,file='Output/rho.txt'		  ,status='replace',action="write")
open(unit=7,file='Output/velocity.txt'	  ,status='replace',action="write")
open(unit=8,file='Output/pressure.txt'	  ,status='replace',action="write")
open(unit=9,file='Output/energy.txt'	  ,status='replace',action="write")
open(unit=10,file='Output/Position.txt'	  ,status='replace',action="write")
end subroutine

subroutine CloseFiles()
implicit none
close(1) ; close(2) ; close(3)
close(4) ; close(7) ; close(8)
close(9) ; close(10); close(11)
end subroutine
