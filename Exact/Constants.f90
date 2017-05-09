module constants
implicit none
real,parameter :: gamma = 5./3.     ,&
				  gm1   = gamma-1.  ,&
				  gp1   = gamma+1.
integer,parameter :: itterMax = 2000
end module
