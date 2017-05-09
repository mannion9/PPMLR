module CommonData
! Nm    - Number of cells in true domain
! Nt    - Bound on number of time steps (Needed for plotting against excact)
! i_min - Index of left  most cell of full domain (including ghost cells)
! i_max - Index of right most cell of full domain (including ghost cells)
! iMIN  - Index of left  most cell of true domain
! iMAX  - Index of right most    cell of true domain 
! rMIN  - Cell center of left  most real cell
! rMAX  - Cell center of right most real cell
! t_max - Bound of time in domain
! gm    - Ratio of specific heats (ideal gas -> gm=1.4)
! alpha - Dimension parameter
! ap1   - alpha+1
! delta - Machine precision on zero test
! COURANT- Courant number
! dr_uni- Unifrom Eulerian cell widths
! r_min - Left most ghost cell center
! r_max - Right most ghost cell center
implicit none
integer,parameter :: BoundaryCon = 0 , & ! (0-Transmisive boundaries 1-Periodic boundaries)
                     InitialCon  = 1 , & ! (1-Sod 2-Flat to right 3-Gaussian  , 4-Traingle)
                     ReZone      = 0 , & ! (1-Remapping on 0- off)
                     RP_Effec    = 0 , & ! (0-Cell average, 1- Domain of dependence )
					 UniSpacing  = 1    ! (0-Uniform spacing, 1- Shock Spacing)
integer,parameter :: Nm = 20,Nt = 10
!integer,parameter :: Nm=100 ,Nt = 10
integer,parameter :: i_min=-3,i_max=Nm+3 , iMIN=0,iMAX=Nm-1,alpha=0,ap1 = alpha+1!Nm is number of interior real
real,parameter    :: rMIN=0.,rMAX=1.,t_max=.25,gm = 1.4,  &
                     delta=1.0d-30 , TOL=1.0d-5 , COURANT = 0.5                  
end module

program main
use CommonData
implicit none
real,dimension(i_min:i_max) :: r,dr
real,dimension(i_min:i_max+1) :: r_12_i=2.d0

call InitDomain(r,r_12_i,dr)

open(unit=1,file='Output_NonUniSpace/Edges.txt')
open(unit=2,file='Output_NonUniSpace/width.txt')
write(1,*) r_12_i(iMIN:iMAX+1)
write(2,*) dr(iMIN:iMAX)
close(1)
close(2)
close(3)
end program

subroutine InitDomain(r,r_12_i,dr)
!-----------------------------------------------------
! Assigns value to cell centers and cell edges and cell widths
! and writes out cell center locations.
! Input/Output :
!    r         - Cell centers
!   r_12_i    - Cell edges
!   dr      - Cell widths
!------------------------------------------------------
use CommonData
implicit none
real,intent(inout),dimension(i_min:i_max)   :: r , dr
real,intent(inout),dimension(i_min:i_max+1) :: r_12_i
real :: dr_uni,r_min,r_max,f,xL,xR,dxS,x_max_L,x_min_L,x_min_R,x_max_R,dx_uni_L,dx_uni_R
integer :: i,NS,Nremain,NL,NR

if (UniSpacing.EQ.0) then 
	! Uniform Spacing
	dr_uni = (rMAX-rMIN)/REAL(Nm-1)  ! uniform spacing
	r_min = rMIN -dr_uni*(iMIN-i_min) ! edge of left  ghost cell
	r_max = rMAX +dr_uni*(i_max-iMAX) ! edge of right ghost cell
	r_12_i = (/(r_min+(i-.5)*dr_uni,i=0,SIZE(r_12_i)-1,1)/)
else 
	f = .50  			! Fraction of all cells, that wil be in refined region
	xL = .4			! Left  boundary of refined region
	xR = .6				! Right boundary of refined region
	NS = INT(f*Nm)	    ! Number of cell in refined region
	Nremain = Nm - NS   ! Number of cells in left and right region
	NL = int(.5*Nremain)! Number of cells in left region
	NR = int(.5*Nremain)! Number of cells in rigt region
	x_max_L = xL - dxS  ! Begining of refined region cell edges
	x_max_R = rMAX		! End 	   of refined region cell edges
	x_min_L = rMIN      ! Begining of left    region cell edges
	x_min_R = xR + 3.*dxS ! Beginign of right region cell edges  ! The 2*dxS ensure the right sidde starts in correct space coord.
	dxS = (xR-xL)/(NS-1)				! Spacing in refined region
	dx_uni_L = (x_max_L-x_min_L)/(NL-1) ! Spacing in left  region
	dx_uni_R = (x_max_R-x_min_R)/(NR-1) ! Spacing in right region
	
	if (dx_uni_L.LT.dxS.OR.dx_uni_R.LT.dxS) then
		print*,'Your spacing in the shock will be worse than in other regions.'
		print*,'Must either increase the fraction of the toal cells to put in the shock region'
		print*,'or increase shock region limits.'
	else if (NL.LE.1) then
		print*,'Numer of cells in refined region less than or equal to one.'
		print*,'Choose a larger fraction or a if you have a small region or'
		prinT*,'choose a smaller fraction if you have a large region.'
	end if 
	! Interior region
	r_12_i(iMIN:NL-1) = (/(x_min_L+i*dx_uni_L ,i=-1,NL-2)/)
	r_12_i(NL:NL+NS-1)  = (/(r_12_i(NL-1)+i*dxS ,i=1,NS)/)
	r_12_i(NL+NS:iMAX+1) = (/(r_12_i(NL+NS-1)+i*dx_uni_R,i=1,NR+1)/)
	! Ghost cells
	r_12_i(i_min:iMIN-1) =(/(r_12_i(iMIN)+i*dx_uni_L,i=-(iMIN-i_min),-1)/)  ! Thiis weird indexing just reverse the ordering
	r_12_i(iMAX+2:i_max+1) = (/(r_12_i(iMAX+1)+i*dx_uni_R,i=1,i_max-iMAX)/)
end if 

! Calculate cell centers and widths
r(i_min:i_max) = .5*(r_12_i(i_min:i_max)+r_12_i(i_min+1:i_max+1))
dr(iMIN:iMAX)= r_12_i(iMIN+1:iMAX+1)-r_12_i(iMIN:iMAX)  ! Cell widths


 print*,r_12_i(iMIN:iMAX+1)
 print*,'  '
 print*,'  '
 print*,dr(iMIN:iMAX)
 print*,'   '
! print*,x_max_L
! print*,x_min_R


! Write out cell center positions
open(unit=1,file='Output_NonUniSpace/CellCenter.txt')
write(1,*) r(iMIN:iMAX)
close(1)

end subroutine