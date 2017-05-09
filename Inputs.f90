module constants
implicit none
!----------------------- USER INPUTS ----------------------------------------- !
integer,parameter :: N 	      = 100   ,& ! Number of cells + 1
					 itterMax = 300  ,& ! Maximum itterations
					 writeStep= 10    ,& ! Step interval between writing out data
					 IC       = 1     ,&   ! Initial condiiton (0-pisition,1-Sod, 6-compression , 7-Sedove, 8- Blast Wave, 9-Constant state moving to righ
					 alpha    = 0         ! Geometric Factor (0-planar,1-cylinderical,2-spherical)
real,parameter    :: Courant  = .4    ,&
					 rmin     = 0.0   ,&
					 rmax     = 1.0   ,&
					 tmax     = 0.4   ,&
					 dtmin    = 1.e-10,&   ! Minimum time step
					 delta    = 1.e-30,&   ! Small number
					 infty    = 1.e30      ! Huge number
character (LEN=*),parameter :: Solver  = 'TwoShock'  ,& ! ('HLL','HLLC','TwoShock')
							   Method  = 'PPM'   ,& ! ('Godunov','MUSCL','PPM')
							   Eulerian= 'On'   ,& ! ('On','Off')
							   Remaper = Method  ,&
							   BC      = 'Free'      !('Free','Reflective')
!----------------------- Program Constants ----------------------------------------- !
integer,parameter :: jmin     = 1          ,& ! Index of left  most LHS real  cell center
					 jmax     = (N-1)+jmin ,& ! Index of right most RHS real  cell center
					 imin     = jmin-3     ,& ! Index of left  most LHS ghost cell center
					 imax     = jmax+4     ,& ! Index of right most RHS ghost cell center
					 ap1	  = alpha+1

end module
