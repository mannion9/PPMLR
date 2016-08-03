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
!integer,parameter :: Nm = 10,Nt = 10
integer,parameter :: Nm=100 ,Nt = 10
integer,parameter :: i_min=-3,i_max=Nm+3 , iMIN=0,iMAX=Nm-1,alpha=0,ap1 = alpha+1!Nm is number of interior real
real,parameter    :: rMIN=0.,rMAX=1.,t_max=.25,gm = 1.4,  &
                     delta=1.0d-30 , TOL=1.0d-5 , COURANT = 0.5                  
end module

program main
use CommonData
implicit none
real,dimension(i_min:i_max) :: r,dr
real,dimension(i_min:i_max+1) :: r_12_i

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
real :: dr_uni,r_min,r_max,refine,StartR,EndR,outter_dr,inner_dr
integer :: nc,n,nl,nr
integer :: i 

dr_uni = (rMAX-rMIN)/REAL(Nm-1)  ! uniform spacing
r_min = rMIN -dr_uni*(iMIN-i_min) ! edge of left  ghost cell
r_max = rMAX +dr_uni*(i_max-iMAX) ! edge of right ghost cell
if (UniSpacing.EQ.0) then 
	! Uniform Spacing
	r      = (/(r_min+i*dr_uni,i=0,SIZE(r)-1,1)/)
	r_12_i = (/(r_min+(i-.5)*dr_uni,i=0,SIZE(r_12_i)-1,1)/)
else 
	refine = 1
	StartR = .4
	EndR   = .6
	outter_dr = dr_uni*refine
	inner_dr  = dr_uni/refine
	r_12_i(iMIN) = r_min-outter_dr
	do i=iMIN+1,iMAX+1
		if (r_12_i(i-1).GT.StartR.AND.r_12_i(i-1).LE.EndR) then 
			! refined
			r_12_i(i) = r_12_i(i-1) + inner_dr
		else
			r_12_i(i) = r_12_i(i-1) + outter_dr
		end if 
	end do 
end if 
r(iMIN:iMAX) = .5*(r_12_i(iMIN:iMAX)+r_12_i(iMIN+1:iMAX+1))
dr(iMIN:iMAX)= r_12_i(iMIN+1:iMAX+1)-r_12_i(iMIN:iMAX)  ! Cell widths
open(unit=3,file='Output_NonUniSpace/CellCenter.txt')
write(3,*) r(iMIN:iMAX)
close(3)

end subroutine