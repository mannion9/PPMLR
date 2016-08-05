! To do
! - Domain of dependence
! - Discontinuity detection
! - Remapping check
! Ideas to spped up code
! - When calling interpol, just put uR as the last input instead of a_12,
!     unless you want to write out a_12

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
                     RP_Effec    = 4 , & ! (0-Cell average, 1- Domain of dependence , 2-Manually Transl. , 3- Triangle Profile  , 4- Use cell averages)
					 UniSpacing  = 0    ! (0-Uniform spacing, 1- Shock Spacing)
!integer,parameter :: Nm = 10,Nt = 10
integer,parameter :: Nm=100 ,Nt = 100
integer,parameter :: i_min=-3,i_max=Nm+3 , iMIN=0,iMAX=Nm-1,alpha=0,ap1 = alpha+1!Nm is number of interior real
real,parameter    :: rMIN=0.,rMAX=1.,t_max=.25,gm = 1.4,  &
                     delta=1.0d-30 , TOL=1.0d-5 , COURANT = 0.8                  
end module

program main
use CommonData
implicit none
real :: total_mass,current_mass,time=0.0,dt=1.0d+30  ! dt=Large number for first ittertaiton, not important
real,dimension(i_min:i_max+1) :: r_12_i,r_12_ip
real,dimension(i_min:i_max)   :: dm,dr,r
real,dimension(i_min:i_max)   :: itters  ! DO NOT NEED JUST FOR DEBUGGING
real,dimension(3,i_min:i_max) :: U,V
integer :: i 

!!Debugging to show cell integer values
do i = i_min,i_max
    itters(i) = real(i)
 end do 
!print*,'itt :',itters(iMIN:iMAX)

! Tell plotting software if this is pure Lagrange or L+R
open(unit=1,file='Output/RemapChoice.txt')
write(1,*) ReZone
close(1)

! Create cells and initial conditions
CALL InitDomain(r,r_12_i,dr)
CALL intialCondition(r,U,V,InitialCon)

dm = dr*V(1,:)
total_mass = SUM(dm(iMIN:iMAX))  ! Total Mass initially in system
current_mass = total_mass        ! Total Mass in system at each step

CALL OpenFiles()                   ! Open all files to write out to

i=0
do while (time.LE.t_max .AND. i.LT.NT)
    CALL boundaries(dr,dm,U,V,r_12_i,r,BoundaryCon)         ! Impose boundary conditions for ghost cells (0-Transmissive ; 1-Periodic)
	CALL TimeStep(U,V,dr,dt)                       ! Deterimne time step
	if (MOD(i,10).EQ.0) then
		CALL Output(U,V,current_mass,dt,r_12_i,i,1,0)  ! Write out data
	end if 
    CALL LagrangeStep(U,V,dm,dt,r_12_i,r_12_ip)    ! Preform Lagrange step
    if (ReZone.EQ.1) then
        CALL Remap(dr,r_12_i,r_12_ip,U,V,dm)
    end if 
    CALL Update(i,time,dt,r_12_i,r_12_ip,dm,dr,current_mass,total_mass,V(1,:)) ! Update 
end do 

if (i.NE.NT+1) then ! Completed without issues
	print*,' ' 
	if (time.GE.t_max) then 
		print*,'Reached the correct final time:',time,'after',i,'itterations.'
	else
		print*,'Did not reach the correct final time, but did not break conservation laws'
		print*,'Reached a final time of :',time,',after:',i,'itterations:'
	end if 
end if 
CALL CloseFiles()
end program

subroutine LagrangeStep(U,V,dm,dt,r_12_i,r_12_ip)
!-------------------------------------------------
! Time marches forward one Lagrange using CW Eq. 2.10
! Input:
!   dm     - Mass in Lagrange cell
!   dt     - Current time step
! Input/Output:
!   U - Vector of variables (tau,u,E)
!   V - Vector of variables (rho,e,P)
!   r_12_i - Lagrange cell position edge before step
!   r_12_ip- Lagrange cell position edge after step
!--------------------------------------------------
use CommonData
implicit none
real,intent(inout),dimension(i_min:i_max+1)  :: r_12_i,r_12_ip
real,intent(inout),dimension(3,i_min:i_max)  :: U,V
real,intent(in),dimension(i_min:i_max)       :: dm
real,intent(in) :: dt
real,dimension(i_min:i_max) :: u_star,p_star,Abar  
real,dimension(i_min:i_max) :: pL=1.,pR=1.,uL=1.,uR=1.,c_p,c_m
integer :: i

call RiemannEffectiveStates(U,V,r_12_i,pL,pR,uL,uR,c_p,c_m,dt,dm,RP_Effec)
call RiemannSolver(pL,pR,uL,uR,c_p,c_m,p_star,u_star)

if (RP_Effec.EQ.2) then
	! Impose solution to RP
	u_star(iMIN:iMAX+1)  = 1.
	p_star(iMIN:iMAX+1) = 1.
else if (RP_Effec.EQ.3) then 
	! Impose triangular 
	u_star(iMIN:iMAX/2)     =  (/(real(i),i=iMIN,iMAX/2,1)/)
	u_star(iMAX/2+1:iMAX+1) =  (/(-real(i),i=0,iMAX/2+1,1)/) + u_star(iMAX/2) !  -real((iMAX/2)+1:iMAX+1) + r(iMAX/2)
else if (RP_Effec.EQ.4) then
	! Impose cell averages for interface values
	u_star(iMIN:iMAX+1) = .5*(U(2,iMIN-1:iMAX)+U(2,iMIN:iMAX+1))
	p_star(iMIN:iMAX+1) = .5*(V(3,iMIN-1:iMAX)+V(3,iMIN:iMAX+1))
end if 

r_12_ip(iMIN:iMAX+1)  = r_12_i(iMIN:iMAX+1) + dt*u_star(iMIN:iMAX+1)

do i=iMIN,iMAX+1    ! This code is needed b/c if u_star==0 then Abar = 1.
    if (u_star(i).LE.TINY(dt).OR.alpha.EQ.0) then  
		! CHECK THE LIMT OF U_STAR-->0 THAT THIS ACTUALLY IS CORRECT
        ! TINY's argument just needs to be a real variable
        Abar(i) = 1.
    else 
        Abar(i)  = (r_12_ip(i)**(ap1)-r_12_i(i)**(ap1))/(ap1*u_star(i)*dt)
    end if 
end do

! Finite difference
U(1,iMIN:iMAX) = ((r_12_ip(iMIN+1:iMAX+1)**(ap1)-r_12_ip(iMIN:iMAX)**(ap1)))/(ap1*dm(iMIN:iMAX))      ! Tau
U(2,iMIN:iMAX) = U(2,iMIN:iMAX) - (dt/(dm(iMIN:iMAX)))*(V(3,iMIN+1:iMAX+1)-V(3,iMIN:iMAX))
U(3,iMIN:iMAX) = U(3,iMIN:iMAX) - (dt/(dm(iMIN:iMAX)))*(U(2,iMIN+1:iMAX+1)*V(3,iMIN:iMAX) & 
					+U(2,iMIN:IMAX)*(V(3,iMIN+1:iMAX+1)-2.*V(3,iMIN:iMAX)))

! Finite volume with Abar = 1 already
! U(1,iMIN:iMAX) = ((r_12_ip(iMIN+1:iMAX+1)**(ap1)-r_12_ip(iMIN:iMAX)**(ap1)))/(ap1*dm(iMIN:iMAX))      ! Tau
! U(2,iMIN:iMAX) = U(2,iMIN:iMAX) + (dt/dm(iMIN:iMAX))*(p_star(iMIN:iMAX)-p_star(iMIN+1:iMAX+1)) ! U
! U(3,iMIN:iMAX) = U(3,iMIN:iMAX) + (dt/dm(iMIN:iMAX))*(u_star(iMIN:iMAX)*p_star(iMIN:iMAX) & 
									 ! -u_star(iMIN+1:iMAX+1)*p_star(iMIN+1:iMAX+1))
									
									
! Finite volume
! U(1,iMIN:iMAX) = ((r_12_ip(iMIN+1:iMAX+1)**(ap1)-r_12_ip(iMIN:iMAX)**(ap1)))/(ap1*dm(iMIN:iMAX))      ! Tau
! U(2,iMIN:iMAX) = U(2,iMIN:iMAX) + &
!                 .5*(Abar(iMIN+1:iMAX+1)+Abar(iMIN:iMAX))*(dt/dm(iMIN:iMAX))*(p_star(iMIN:iMAX)-p_star(iMIN+1:iMAX+1)) ! U
!U(3,iMIN:iMAX) = U(3,iMIN:iMAX) + (dt/dm(iMIN:iMAX))*(Abar(iMIN:iMAX)*u_star(iMIN:iMAX)*p_star(iMIN:iMAX) &
!                                 - Abar(iMIN+1:iMAX+1)*u_star(iMIN+1:iMAX+1)*p_star(iMIN+1:iMAX+1))


! Update primatives
V(1,:) = 1./U(1,:)               ! Rho
V(2,:) = U(3,:) - .5*U(2,:)**(2) ! Internal energy
V(3,:) = (gm-1.)*V(1,:)*V(2,:)   ! Pressure
end subroutine

! Riemann Solver Subroutines

subroutine RiemannSolver(pL,pR,uL,uR,c_p,c_m,p_star,u_star)
!---------------------------------------
! Solves RP for two shock waves for p_star 
! and then solves for u_star. This is an HLL solver.
! Our initial guess p_star is the average of the
! left and right states values.
! Input:
!   pL/pR - Pressure left and right states
!   uL/uR - Velocity left and right states
!   c_p   - Speed of sound from left  state
!   c_m   - Speed of sound from right state
! Input/Output:
!   p_star- Pressure between cells
!   u_star- Velocity between cells
!---------------------------------------
use CommonData
implicit none
real,intent(in),dimension(i_min:i_max)   :: pL,pR,uL,uR,c_m,c_p
real,intent(inout),dimension(i_min:i_max):: p_star,u_star  
real :: diff,Temp,Wl,del_u,f,f_var,fprime,gmp1=gm+1
integer :: j, i, itterMAX = 20
!p_star = .5*(pL+pR)! initial guess value
do j=iMIN,iMAX+1  ! For each cell interface
	diff = 1.       ! Intialize doesnt matter for first itteration
	del_u = uR(j)-uL(j)
	p_star  = (c_m*pL+c_p*pR-c_p*c_m*del_u)/(c_p+c_m) ! Van Leer Eq. 60 
	do i=1,itterMAX
		if (diff.GT.TOL) then
			Temp  = p_star(j)
			call Func_RP(c_p(j),c_m(j),p_star(j),pL(j),pR(j),del_u,f)        ! f(x)
			call Func_RP(c_p(j),c_m(j),p_star(j)+TOL,pL(j),pR(j),del_u,f_var)! f(x+dx)
			fprime = (f_var-f)/TOL                                           ! f'(x) ~ (f(x)-f(x+dx))/dx    
			if (fprime.LE.delta .OR. fprime.NE.fprime .OR. p_star(j).LT.delta) then ! ensures no divide by zero or NAN or negative pressure
				p_star(j) = Temp
				exit
			end if 
		p_star(j) = p_star(j) - f/fprime                    ! Newtons method
		diff      = 2.*ABS((p_star(j)-Temp)/(p_star(j)+Temp)) ! Difference from last value
		else 
			exit
			end if
	end do
	if (i.EQ.itterMAX) then
		print*,'Warning Riemann solver has reached its maximum # of itteration, may not have converged.'
	end if 
Wl = c_p(j)*SQRT(ABS((1.+(gmp1/(2.*gm))*((p_star(j)/pL(j))-1.))))
u_star(j)  = -(p_star(j)-pL(j))/Wl + uL(j)  
end do 
! print*,'sound:',c_p(iMIN)
! print*,'Diff:',p_star(iMIN)-pL(iMIN)
! print*,'Lagrange Shock speed:', SQRT(ABS(c_p(iMIN)*(1.+(gmp1/(2.*gm))*((p_star(iMIN)/pL(iMIN))-1.))))
! Print RP solution
! print*,'Solution to RP'
! print*,'cp:',(c_p(iMIN:iMAX+1))
! print*,'cm:',(c_m(iMIN:iMAX+1))
! print*,'pL:',pL(iMIN:iMAX+1)
! print*,'pR:',pR(iMIN:iMAX+1)
! print*,'uL:',uL(iMIN:iMAX)
! print*,'uR:',uR(iMIN:iMAX)
! print*,'p_star:',p_star(iMIN:iMAX+1)
! print*,'u_star:',u_star(iMIN:iMAX+1)

end subroutine

subroutine RiemannEffectiveStates(U,V,r_12_i,pp,pm,up,um,c_p,c_m,dt,dm,choice)
!------------------------------------------------------------------
! Calculate the effective right and left states
! by taking into account the domain of dependence
! using CW Eq. 2.6 , which is using the f from Eq. 1.12
! Input:
!   U - Vector of variables (tau,u,E)
!   V - Vector of variables (rho,e,P) 
!   dt - Current time step
!   dm - Mass in Lagrange cell
!   r_12_i - Lagrange cell edges
!   choice - Set to 1 to use correct domain of dependence, set to 0 to use full cell average
! Input/Output:
!   pp/pm - Pressure from left and right state (These values are assigned to pL/pR when outside of the subroutine because we are passing in pL/pR, dont confuse this)
!   up/um - Velocity from left and right state
!   cp/cm - Lagrange speed of sound from left and right states
!------------------------------------------------------------------
use CommonData
implicit none
real,intent(inout),dimension(i_min:i_max) :: pp,pm,up,um,c_p,c_m
real,intent(in),dimension(3,i_min:i_max)  :: U,V
real,intent(in),dimension(i_min:i_max)       :: r_12_i,dm
integer,intent(in) :: choice
real,intent(in) :: dt
real,dimension(i_min:i_max) ::   rhoR,rhoL,delrho,rho6,rho12, &
                                 uR  ,uL  ,delu  ,u6  ,u12,&
                                 pR  ,pL  ,delp  ,p6  ,p12, &
                                 cs  , A , dr , rhop , rhom ,z
real    :: twth = 2./3.
integer :: i,j

! This goes from iMIN-1:iMAX+1 because we need to calculate effective states 
! from the zero interface to the iMAX+1 interface, so we need the first ghost
! cell on left and right.
dr(iMIN-1:iMAX+1) = r_12_i(iMIN:iMAX+2)-r_12_i(iMIN-1:iMAX+1)
cs(iMIN-1:iMAX+1) = SQRT(gm*V(3,iMIN-1:iMAX+1)*V(1,iMIN-1:iMAX+1))

 if (alpha.EQ.0) then    ! If alpha = 1 then A reduces to 1 analyitically, this guarentees that
        A(iMIN-1:iMAX+1) = 1.
    else 
        A(iMIN-1:iMAX+1) = (r_12_i(iMIN:iMAX+2)**(ap1)-r_12_i(iMIN-1:iMAX+1)**(ap1))/(ap1*dr(iMIN-1:iMAX+1))
end if 


call interpol(V(1,:),dm,rho12,rhoR,rhoL,delrho,rho6,0)
call interpol(U(2,:),dm,u12  ,uR  ,uL  ,delu  ,u6  ,0)
call interpol(V(3,:),dm,p12  ,pR  ,pL  ,delp  ,p6  ,0)

z  = dt*cs*A/dm

if (choice==1) then ! Actually do the correct domain of dependence using flux integrals
    do i=iMIN,iMAX+1  ! For each interface we need a left and right state 
        ! data from the state on the left of the interfaces
        j = i -1 
        pp(i)  = pR(j)  - 0.5*z(j)*(delp(j) - (1.-twth*z(j))*p6(j)  )
        up(i)  = uR(j)  - 0.5*z(j)*(delu(j) - (1.-twth*z(j))*u6(j)  )
        rhop(i)= rhoR(j)- 0.5*z(j)*(delrho(j)-(1.-twth*z(j))*rho6(j))
        c_p(i) = SQRT(gm*pp(i)*rhop(i)) ! The lagrange speed of sound
        ! Data from the state on the right of the interfaces
        pm(i)  = pL(i)  + 0.5*z(i)*(delp(i)  +(1.-twth*z(i))*p6(i)  )
        um(i)  = uL(i)  + 0.5*z(i)*(delu(i)  +(1.-twth*z(i))*u6(i)  )
        rhom(i)= rhoL(i)+ 0.5*z(i)*(delrho(i)+(1.-twth*z(i))*rho6(i))
        c_m(i) = SQRT(gm*pm(i)*rhom(i)) ! The lagrange speed of sound
    end do 
else  ! Just use full cell average as left and right states
    rhop(iMIN:iMAX+1) =  V(1,iMIN-1:iMAX)   ! Left state for each interface
    rhom(iMIN:iMAX+1) =  V(1,iMIN:iMAX+1)    ! Right state for each interface 
    pp(iMIN:iMAX+1)   =  V(3,iMIN-1:iMAX)    
    pm(iMIN:iMAX+1)   =  V(3,iMIN:iMAX+1)
    up(iMIN:iMAX+1)   =  U(2,iMIN-1:iMAX)
    um(iMIN:iMAX+1)   =  U(2,iMIN:iMAX+1)
    c_p(iMIN:iMAX+1) = SQRT((gm*pp(iMIN:iMAX+1)*rhop(iMIN:iMAX+1))) ! The Lagrange speed of sound from left cell
    c_m(iMIN:iMAX+1) = SQRT((gm*pm(iMIN:iMAX+1)*rhom(iMIN:iMAX+1))) ! The Lagrange speed of sound from right cell
end if 
! Print effective states
!print*,'diff r-',MAXVAL(rhom(iMIN:iMAX+1)- V(1,iMIN:iMAX+1))
!print*,'diff u-',MAXVAL(um(iMIN:iMAX+1)- U(2,iMIN:iMAX+1))
!print*,'diff P-',MAXVAL(pm(iMIN:iMAX+1)- V(3,iMIN:iMAX+1))
! print*,A(iMIN-1),dm(iMIN-1)
! print*,'coef:',pR(iMIN-1),z(iMIN-1),delp(iMIN-1),p6(iMIN-1)
! print*,'c :',cs(iMIN-1)
! print*,'p-:',pm(iMIN)
! print*,'u-:',um(iMIN)
! print*,'p+:',pp(iMIN)
! print*,'u+:',um(iMIN)
end subroutine 

subroutine func_RP(cp,cm,p,pl,pr,del_u,eval)
!----------------------------------
! Calculates RP funciton to find root of 
! from CW Eq. 2.8 . The acutal equation is found 
! by adding the 2.8 a) and 2.8 c) to one equation
! and inserting 2.8 b) and 2.8 d).
! Input:
!    cp    - Square of Lagrange sound speed from left
!    cm    - Square of Lagrange sound speed from right 
!    p     - Current pressure estimate from Riemann solver
!    pl    - Pressure in left state
!    pr       - Pressure in right state
!    del_u - Difference in right and left state velocity
! Input/Output:
!    eval - Result of function evaulation
!----------------------------------
use CommonData
implicit none
real,intent(in) :: cp,cm,p,pl,pr,del_u
real,intent(inout) :: eval
real :: Wl,Wr
Wl   = SQRT(ABS(cp*(1.+((gm+1.)/(2.*gm))*((p/pl)-1.))))
Wr   = SQRT(ABS(cm*(1.+((gm+1.)/(2.*gm))*((p/pr)-1.))))
eval = (p-pL)/Wl + (p-pR)/Wr - del_u
end subroutine

! Interpolation Subroutines

subroutine interpol(a,dx,a12,aR,aL,dela,a6,detect)
!-------------------------------
! Calculate interpolation value of a_(j+1/2)
! along with the coefficants of the interpolation
! parabola using Eq. 1.5 and 1.6
! Input:
!   a        - Current solution
!   dx       - Cell widths 
!    detect- Set to one for discontinuity detection
! Input/Output:
!   a12   - Interpolated value at cell interfaces
!   aR/aL - Cell interface from right and left cell values
!    dela  - Difference in aR and aL 
!    a6      - Interpolation coefficant
!------------------------------
use CommonData
implicit none
real,intent(inout),dimension(i_min:i_max) :: a12,aR,aL,dela,a6
real,intent(in),dimension(i_min:i_max)    :: a,dx
integer,intent(in) :: detect
real,dimension(iMIN-2:iMAX+1) :: one,two,three,four,five
real,dimension(i_min:i_max)   :: del,aRd,aLd
integer :: js,je
js = iMIN-2   ! First cell we interpolate is the i=-2  cell
je = iMAX+1   ! Last  cell we interpolate is the i=Nm  cell
one   = dx(js:je)/(dx(js:je)+dx(js+1:je+1))
two   = 1./(dx(js-1:je-1)+dx(js:je)+dx(js+1:je+1)+dx(js+2:je+2))
three = (2.*dx(js+1:je+1)*dx(js:je))/(dx(js:je)+dx(js+1:je+1))
four  = (dx(js-1:je-1)+dx(js:je))/(2.*dx(js:je)+dx(js+1:je+1))
five  = (dx(js+2:je+2)+dx(js+1:je+1))/(2.*dx(js+1:je+1)+dx(js:je))

CALL avgSlope(a,dx,del)
a12(js:je) = a(js:je) + one*(a(js+1:je+1)-a(js:je)) &
            + two*(three*(four-five)*(a(js+1:je+1)-a(js:je)) &
            - dx(js:je)*four*del(js+1:je+1)+dx(js+1:je+1)*five*del(js:je))            
aR(iMIN-1:iMAX+1) = a12(iMIN-1:iMAX+1)  
aL(iMIN-1:iMAX+1) = a12(iMIN-2:iMAX)    
if (detect==1) then
    call discont_detect(a,dx,aR,aRd,aL,aLd,del)
end if 
CALL monotone(a,aL,aR)
dela = aR-aL
a6   = 6.*(a-.5*(aL+aR))            

end subroutine

subroutine avgSlope(a,dx,delmj)
!------------------------------
! Calculate the average slope in the jth cell of
! the parabola with zone average a j,j+1,j-1 from
! CW Eq. 1.7 and Eq. 1.8
! Input:
!   a - Current Solution 
!   dx- Cell widths 
! Input/Output:
!   delj - Average slove 
!-----------------------------
use CommonData
implicit none
real,intent(in),dimension(i_min:i_max)    :: a,dx
real,intent(inout),dimension(i_min:i_max) :: delmj
real,dimension(iMIN-1:iMAX+2) :: one,two,three,diff1,diff2,delj
integer :: st0,stm1,stp1,en0,enm1,enp1,i
real :: sgn
! These value are the start (st#) and the end (en#) 
! of the parts of a and dx that we evaluate in Eq. 1.8
! We will assign value to the iMIN-1-iMAX+2 elements of 
! the delmj vector.
st0 = iMIN-1 ! goes from -1 ghost cell
en0 = iMAX+2 ! up until the +1 ghost cell 
stm1= st0 -1
enm1= en0 -1
stp1= st0 +1
enp1= en0 + 1
one   = dx(st0:en0)/(dx(stm1:enm1)+dx(st0:en0)+dx(stp1:enp1))
two   = (2.*dx(stm1:enm1)+dx(st0:en0))/(dx(stp1:enp1)+dx(st0:en0))
three = (dx(st0:en0)+2.*dx(stp1:enp1))/(dx(stm1:enm1)+dx(st0:en0)) 
delj  = one*(two*(a(stp1:enp1)-a(st0:en0))+three*(a(st0:en0)-a(stm1:enm1)))
diff1 = a(stp1:enp1)-a(st0:en0)
diff2 = a(st0:en0)  -a(stm1:enm1)
do i=st0,en0
  if (diff1(i)*diff2(i).GT. 0.) then
    if (delj(i).LT.0.) then
      sgn = -1.
    else 
      sgn = 1.
    end if 
    delmj(i) = SIGN( MIN(ABS(delj(i)),2.*ABS(diff1(i)),2.*ABS(diff2(i))),sgn)
  else 
    delmj(i) = 0. 
  end if 
end do
end subroutine

subroutine discont_detect(a,dx,aR,aRd,aL,aLd,del)
use CommonData
implicit none
real,intent(in),dimension(i_min:i_max) :: a,dx,del
real,intent(inout),dimension(i_min:i_max) :: aL,aR
real,dimension(i_min:i_max) :: aLd,aRd,delsquare,nabla
real :: eps = 0.01,nabla1=20.,nabla2=0.05, nabla_temp 
integer :: i,js,je
js = iMIN-1
je = iMAX+2
aRd(iMIN:iMAX+1) = a(iMIN+1:iMAX+2) - .5*del(iMIN+1:iMAX+2)
aLd(iMIN:iMAX+1) = a(iMIN-1:iMAX)   + .5*del(iMIN-1:iMAX)
delsquare(js:je) = 1./(dx(js-1:je-1)+dx(js:je)+dx(js+1:je+1)) &
                 * ((a(js+1:je+1)-a(js:je))/(dx(js+1:je+1)+dx(js:je)) &
                 - (a(js:je)-a(js-1:je-1))/(dx(js:je)+dx(js-1:je-1)))
do i=iMIN,iMAX+1
    if (-1.*delsquare(i+1)*delsquare(i-1).GT.0. .AND. ABS(a(i+1)-a(i-1))-eps*MIN(ABS(a(i+1)),ABS(a(i-1))).GT.0.) then
        nabla_temp = -1.*((delsquare(i+1)-delsquare(i-1))/(dx(i)+dx(i-1)))*((dx(i-1)**(3)+dx(i)**(3))/(a(i+1))-a(i-1))
    else 
        nabla_temp = 0.0
    end if 
    nabla(i) = MAX(0.0,MIN(nabla1*(nabla_temp-nabla2),1.))
end do 
aL = aL*(1.-nabla)+aLd*nabla
aR = aR*(1.-nabla)+aRd*nabla
end subroutine

subroutine monotone(a,aR,aL)
!-------------------------------
! Ensures interpolation is monotone using CW Eq. 1.10
! Input:
!   a  - Current solution
! Input/Output:
!   aR - Approx. of a_(j+1/2) from cell to right
!   aL - Approx. of a_(j-1/2) from cell to left
!-----------------------------
use CommonData
implicit none
real,intent(inout),dimension(i_min:i_max) :: aR,aL
real,intent(in),dimension(i_min:i_max)       :: a
real,dimension(i_min:i_max) :: aR_temp,aL_temp, cond1,cond2,cond3
integer :: i
aR_temp = aR
aL_temp = aL
cond1 = (aR-a)*(a-aL)
cond2 = (aR-aL)*(a-.5*(aL+aR))
cond3 = ((aR-aL)**(2))/6.
do i=iMIN-1,iMAX+1
  if (cond1(i).LE.0.) then
    aL(i) = a(i)
    aR(i) = a(i)
  end if
  if (cond2(i)>cond3(i)) then
    aL(i) = 3.*a(i)-2.*aR_temp(i)
  end if 
  if (-1.*cond3(i)>cond2(i)) then
    aR(i) = 3.*a(i)-2.*aL_temp(i)
  end if 
end do 
end subroutine

! Misc. Subroutines

subroutine TimeStep(U,V,dr,dt)
!-----------------------------------------
! Calculate the required time step with the given Courant number
! and takes into acound the maximum speeds and distances in the problem
! Recall the CFL condition -> MAX(|u|)*dt/MIN(dx) .LE. COURANT
! So we get the MAX for a given COURANT number is given by
!    dt = COURANT*MIN(dx)/MAX(|u|) 
! Another feauture is that this ensure that the time step
! never grows too quickly, if the proposed time step is twice of the old 
! time step, we choose the time step to be slightly larger than the 
! previous time step. 
! Inputs:
!   U - Vector of variables (tau,u,E)
!   V - Vector of variables (rho,e,P)
!   dr- Mass in Lagrange Cell
! Output:
!    dt - Time step
!-----------------------------------------
use CommonData 
implicit none
real,intent(in),dimension(3,i_min:i_max) :: U,V
real,intent(in),dimension(i_min:i_max) :: dr
real,intent(inout) :: dt 
real :: dt_temp , mc ,mu , mspeed
dt_temp = dt                                              ! Previous time step
mc      = MAXVAL(SQRT(gm*V(3,iMIN:iMAX)*V(1,iMIN:iMAX)))  ! Maximum speed of sound 
mu      = MAXVAL(ABS(U(2,iMIN:iMAX)))                     ! Maximum fluid speed
mspeed  = MAX(mc,mu,delta)                                ! Maximum speed in problem (delta ensure non-zero max speed)
dt = COURANT*MINVAL(dr(iMIN:iMAX))/mspeed                 ! Maximum time step 
dt = MIN(dt,2.*dt_temp)

end subroutine 

subroutine boundaries(dr,dm,U,V,r_12_i,r,BC)
!--------------------------------
! Assigns ghost cells values determined 
! by the chosen boundary conditions.
! Transmisive BC - 
!	Sets the ghost cells equal to real cell values but in reverse order
! 	so that the closest ghost cell is set equal to the closest real cell value, the 
! 	farthest ghost cell is set equal to the farthest real cell value, its a reflection.
! Periodic BC-
!	Sets the ghost cellls equal to the real cells at other end of grid
! Inputs:
!    BC - Chosen boundary condition
! Input/Output:
!   dr - Euler cell widths
!   dm - Mass in lagrange cells
!   U - Vector of variables (tau,u,E)
!   V - Vector of variables (rho,e,P)
!	r_12_i - Edges of cells
!	r  - Cell centers
!----------------------------------
use CommonData
implicit none
integer,parameter :: nl=abs(i_min),nr=i_max-iMAX
real,intent(inout),dimension(3,i_min:i_max) :: U,V
real,intent(inout),dimension(i_min:i_max)   :: dr,dm,r
real,intent(inout),dimension(i_min:i_max+1)  :: r_12_i
integer,intent(in) :: BC
integer :: i
IF (BC==0) then ! Transmissive
    ! Update left ghost cells
    dr(i_min:iMIN-1)  = dr(nl-1:iMIN:-1)
    dm(i_min:iMIN-1)  = dm(nl-1:iMIN:-1)
    U(1,i_min:iMIN-1) = U(1,nl-1:iMIN:-1)
    U(2,i_min:iMIN-1) = U(2,nl-1:iMIN:-1)
    U(3,i_min:iMIN-1) = U(3,nl-1:iMIN:-1)
    V(1,i_min:iMIN-1) = V(1,nl-1:iMIN:-1)
    V(2,i_min:iMIN-1) = V(2,nl-1:iMIN:-1)
    V(3,i_min:iMIN-1) = V(3,nl-1:iMIN:-1)
    ! Update right ghot cells
    dr(Nm:i_max)  = dr(iMAX:Nm-nr:-1)
    dm(Nm:i_max)  = dm(iMAX:Nm-nr:-1)
    U(1,Nm:i_max) = U(1,iMAX:Nm-nr:-1)
    U(2,Nm:i_max) = U(2,iMAX:Nm-nr:-1)
    U(3,Nm:i_max) = U(3,iMAX:Nm-nr:-1)
    V(1,Nm:i_max) = V(1,iMAX:Nm-nr:-1)
    V(2,Nm:i_max) = V(2,iMAX:Nm-nr:-1)
    V(3,Nm:i_max) = V(3,iMAX:Nm-nr:-1)
ELSE IF (BC==1) then ! Periodic
    ! Update left ghost cells
    dr(i_min:iMIN-1) = dr(Nm-nl:iMAX)
    dm(i_min:iMIN-1) = dm(Nm-nl:iMAX)
    U(1,i_min:iMIN-1)= U(1,Nm-nl:iMAX)
    U(2,i_min:iMIN-1)= U(2,Nm-nl:iMAX)
    U(3,i_min:iMIN-1)= U(3,Nm-nl:iMAX)
    V(1,i_min:iMIN-1)= V(1,Nm-nl:iMAX)
    V(2,i_min:iMIN-1)= V(2,Nm-nl:iMAX)
    V(3,i_min:iMIN-1)= V(3,Nm-nl:iMAX)
    ! Update right ghost cells
    dr(Nm:i_max) = dr(iMIN:nr-1)
    dm(Nm:i_max) = dm(iMIN:nr-1)
    U(1,Nm:i_max)= U(1,iMIN:nr-1)
    U(2,Nm:i_max)= U(2,iMIN:nr-1)
    U(3,Nm:i_max)= U(3,iMIN:nr-1)
    V(1,Nm:i_max)= V(1,iMIN:nr-1)
    V(2,Nm:i_max)= V(2,iMIN:nr-1)
    V(3,Nm:i_max)= V(3,iMIN:nr-1)
END IF 
! Update ghost cell ceters and edges
do i = 1,3
	r_12_i(iMIN-i)   = r_12_i(iMIN-i+1) - dr(iMIN-i)
	r_12_i(iMAX+1+i) = r_12_i(iMAX+i)   + dr(iMAX+i)
end do 
r(i_min:i_max) = .5*(r_12_i(i_min:i_max)+r_12_i(i_min+1:i_max+1))
end subroutine

subroutine intialCondition(r,U,V,choice)
!--------------------------------
! Creates intial data
! Input:
!     choice - Chooses initial state
!    r      - Euler cell center locations
! Input/Output:
!   U - Vector of variables (tau,u,E)
!   V - Vector of variables (rho,e,P)
!-------------------------------
use CommonData
implicit none
real,intent(inout),dimension(3,i_min:i_max) :: U,V
real,intent(in),dimension(i_min:i_max) :: r
integer,intent(in) :: choice
integer :: i
if (choice==1) then
  print*,'Sod Test'
  do i=iMIN,iMAX
    if (r(i) < .5) then
      V(1,i) = 1.
      U(2,i) = delta
      V(3,i) = 1.
    else
      V(1,i) = 0.125
      U(2,i) = delta
      V(3,i) = 0.1
    end if 
  end do
else if (choice == 2) then
    print*,'Constant'
    V(1,iMIN:iMAX) = 1.
    U(2,iMIN:iMAX) = 1.
    V(3,iMIN:iMAX) = 1.
else if (choice==3) then 
    print*,' Gaussian Pressure Test'
    V(1,iMIN:iMAX) = 1.
    U(2,iMIN:iMAX) = 0
    V(3,iMIN:iMAX) = EXP(-(r(iMIN:iMAX)-.5)**(2)/.1)
else if (choice==4) then
	print*,"Triangle"
	V(1,iMIN:iMAX)     = 1.
	U(2,iMIN:iMAX/2)   = (/(real(r(i)),i=iMIN,iMAX/2,1)/)
	U(2,iMAX/2+1:iMAX) = U(2,iMAX/2:iMIN:-1)
	V(3,iMIN:iMAX) 	   = 1.
end if 
! Assign other variables
U(1,iMIN:iMAX) = 1./V(1,iMIN:iMAX)                       ! tau = 1/rho
V(2,iMIN:iMAX) = V(3,iMIN:iMAX)/((gm-1.)*V(1,iMIN:iMAX)) ! e = P/(gm-1)/rho
U(3,iMIN:iMAX) = V(2,iMIN:iMAX) + .5*U(2,iMIN:iMAX)**(2) ! E = e + .5 u**(2)
end subroutine

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

! Write out cell center positions
open(unit=1,file='Output/CellCenter.txt')
write(1,*) r(iMIN:iMAX)
close(1)

end subroutine

! I/O Subroutines

subroutine Update(i,time,dt,r_12_i,r_12_ip,dm,dr,current_mass,total_mass,rho)
!--------------------------------------------------------------------
! Updates values after full itteration and also checks that we have 
! conserved mass, and that all cells are disentangled. If the mass
! is not conserved or if the cells get tangled, then this subroutine
! triggers the program to terminate.
! Input/Output:
!    r_12_i  - New Eulerian cell edges after one itteration
!   r_12_ip - Cell Eulerian cell edges after Lagrange step and a remapping 
!    dm        - Mass in eulerian cell after remapping, or Mass in Lagrange cell if remap is off
!   dr         - Eulerian cell widths after one itteration
!   rho     - Density
!    time    - Time in simulation
!   dt         - Time step taken in this itteration
!    current_mass - Mass of system at current step
!    total_mass      - Mass of system at initial condition
!    i        - Itteration counter.
!------------------------------------------------------------------
use CommonData
real,intent(inout),dimension(i_min:i_max) :: r_12_i,r_12_ip,dm,dr,rho
real,intent(inout) :: time,dt,current_mass,total_mass
integer,intent(inout) :: i
i = i + 1                                                ! Update stepper
time = time + dt                                         ! Update time 
r_12_i = r_12_ip                                         ! Update cell edges  (Will only change if in full Lagrange mode) ** Remap sets r_12_ip=startionary cells
dr(iMIN:iMAX) = r_12_i(iMIN+1:iMAX+1)-r_12_i(iMIN:iMAX) ! Update cell widths (Will only change if in full Lagrange mode)
dm = dr*rho                                              ! Update the mass in each cell
current_mass = SUM(dm(iMIN:iMAX))                          ! Update total mass
IF (current_mass/total_mass.LT. .9 .OR. current_mass.NE.current_mass) THEN ! If the mass is less then 90% of the original mass break
    PRINT*,'Conservation of mass broken at itteration:',i,'at time:',time
    print*,'Program Terminated Early'
    i = NT+1 ! This will stop Main program
else if (MINVAL(dr(iMIN:iMAX)).LE.0.0) then ! If a cell has been tangled, and now has negative volume break
    print*,'Cell with negative volume has been created at itteration',i,'at time:',time
    print*,'Program Terminated Early'
    i = NT+1 ! This will stop Main program
end if 
end subroutine

subroutine Output(U,V,M,dt,r_12_i,itter,writeOut,printOut)
!---------------------------------
! Will either write out current state of system
! to the terminal or to a series of .txt files
! that have already been opened before the call.
! Input:
!   U  - Vector of variables (tau,u,E)
!   V  - Vector of variables (rho,e,P)
!   M  - Total mass in system at this time step
!    dt - Current time step
!   r_12_i   - Lagrange cell edges
!    itter       - Current time itteration 
!    writeOut - Set to 1 to write out data to .txt files
!   printOut - Set to 1 to write out data to terminal 
!-------------------------------
use CommonData
implicit none
real,intent(in),dimension(3,i_min:i_max) :: U,V
real,intent(in),dimension(i_min:i_max)   :: r_12_i
integer,intent(in) :: itter,writeOut,printOut
real,intent(in) :: dt, M
integer :: st,en,choice 
choice = 0 ! Set choice to 0 to prin out real domain, 1 to print out entire domain (including ghost cells)
if (choice==0) then ! print out real domain
    st = iMIN
    en = iMAX
else if (choice==1) then ! print out full domain 
    st = i_min
    en = i_max
end if 

if (printOut==1) then
    print*,'------------------------------- ITTERATION:',itter,'-------------------------------'
    print*,'Tau:',U(1,st:en)
    print*,'U  :',U(2,st:en)
    print*,'E  :',U(3,st:en)
    print*,'Rho:',V(1,st:en)
    print*,'e  :',V(2,st:en)
    print*,'P  :',V(3,st:en)
    print*,'dt :',dt
    print*,'L_r:',.5*(r_12_i(iMIN:iMAX)+r_12_i(iMIN+1:iMAX+1)) ! Lagrange Cell Ceters
    print*,'M  :',M
end if 

if (writeOut==1) then
    write(1,*) V(1,iMIN:iMAX) ! Rho 
    write(2,*) U(2,iMIN:iMAX) ! Velocity
    write(3,*) V(2,iMIN:iMAX) ! Internal Energy
    write(4,*) V(3,iMIN:iMAX) ! Pressure
    write(7,*) U(3,iMIN:iMAX) ! Total Energy
    write(8,*) dt             ! dt
    write(9,*) .5*(r_12_i(iMIN:iMAX)+r_12_i(iMIN+1:iMAX+1)) ! Lagrange Cell Ceters
    write(10,*) M               ! Total mass
    write(11,*) dt               ! dt for Exact RP solver
end if 
end subroutine

subroutine OpenFiles()
implicit none
! Open Files to write out to 
open(unit=1,file='Output/Density.txt')
open(unit=2,file='Output/Velocity.txt')
open(unit=3,file='Output/InternalEnergy.txt')
open(unit=4,file='Output/Pressure.txt')
open(unit=7,file='Output/Energy.txt')
open(unit=8,file='Output/dt.txt')
open(unit=9,file='Output/LagnCellCenter.txt')
open(unit=10,file='Output/TotalMass.txt')
open(unit=11,file='RP-Exact/Inputs/dt.txt') ! Needed for Exact riemann solver
end subroutine

subroutine CloseFiles()
implicit none
! Open Files to write out to 
close(1) ! Close density
close(2) ! Close velocity
close(3) ! Close Internal e
close(4) ! Close pressure
close(7) ! Close total energy
close(8) ! Close dt for my solution
close(9) ! Close Lagrage cell ceter
close(10)! Close total mass
close(11)! Close dt for RP-Exact
end subroutine

subroutine Remap(dr,r_12_i,r_12_ip,U,V,dm)
use CommonData
implicit none
real,intent(in),dimension(i_min:i_max) :: dr
real,intent(inout),dimension(i_min:i_max+1) :: r_12_i,r_12_ip
real,intent(inout),dimension(3,i_min:i_max) :: U,V
real,intent(inout),dimension(i_min:i_max) :: dm
real,dimension(i_min:i_max) :: rhoR,rhoL,delrho,rho6,rho12, &
                               uR  ,uL  ,delu  ,u6  ,u12, &
                               pR  ,pL  ,delp  ,p6  ,p12, &
                               ER  ,EL  ,delE  ,E6  ,E12, & 
                               dm_ip , dx_ip , mf , uf , Ef
real,dimension(i_min:i_max+1):: deltax

real :: y , twth = 2./3.
integer :: i,j
! Calculate cell widths and overlap
dx_ip   = r_12_ip(i_min+1:i_max+1)-r_12_ip(i_min:i_max)  ! New lagrange zone widths
deltax  = r_12_ip-r_12_i                                 ! Overlap between new and old zone edges

! Interpolate rho,u,P in space domain of updated lagrange cell
call interpol(V(1,:),dx_ip,rho12,rhoR,rhoL,delrho,rho6,0)  !!!!!!! SET THE LAST TERM TO 1 FOR STEEPEING
call interpol(U(2,:),dx_ip,u12  ,uR  ,uL  ,delu  ,u6  ,0)
call interpol(V(3,:),dx_ip,p12  ,pR  ,pL  ,delp  ,p6  ,0)

! Interpolate energy, given the P,rho,u interpolation vlaues
E12=P12/((gm-1.)*rho12) + .5*u12**(2)
ER(iMIN-1:iMAX+1) = E12(iMIN-1:iMAX+1)
EL(iMIN-1:iMAX+1) = E12(iMIN-2:iMAX)
CALL monotone(U(2,:),EL,ER)
delE = ER-EL
E6   = 6.*(U(3,:)-.5*(ER+EL))

! Calculate mass flux
do i = iMIN, iMAX+1
    if(deltax(i).GE.0.0) then
        j = i - 1
        y = deltax(i)/dx_ip(j)
        mf(i) = (rhoR(j) - 0.5*abs(y)*(delrho(j) - (1.-twth*abs(y))*rho6(j)))*deltax(i)
    else
        y = deltax(i)/dx_ip(i)
        mf(i) = (rhoL(i) + 0.5*abs(y)*(delrho(i) + (1.-twth*abs(y))*rho6(i)))*deltax(i)   
    endif
enddo

! Calculate energy/velocity flux
do i = iMIN, iMAX+1
    if(mf(i).GE.0.0) then
        j = i - 1
        y = mf(i)/(dx_ip(j)*V(1,j)) 
        uf(i) = (uR(j)   - 0.5*abs(y)*(delu(j)   - (1.-twth*abs(y))*u6(j))   )*mf(i)
        Ef(i) = (ER(j)   - 0.5*abs(y)*(delE(j)   - (1.-twth*abs(y))*E6(j))   )*mf(i)
    else
        y = mf(i)/(dx_ip(i)*V(1,i))
        uf(i) = (uL(i)   + 0.5*y*(delu(i)     + (1.-twth*y)*u6(i))      )*mf(i)
        Ef(i) = (EL(i)   + 0.5*y*(delE(i)     + (1.-twth*y)*E6(i))      )*mf(i)
    endif
enddo
! Perform Remapping
do i = iMIN, iMAX
    dm(i)    = V(1,i)*dx_ip(i)                    ! Mass in the lagrange cell after lagrange step
    dm_ip(i) = dm(i) + mf(i) - mf(i+1)            ! Mass in remaped cell
    V(1,i)   = dm_ip(i)/dr(i)                     ! Density in remapped cell (rho = mass/volume)
    V(1,i)   = max(-0.0000001,V(1,i))             ! Lower limit on density
    dm_ip(i) = (V(1,i)*dr(i))                     ! Recalculate mass after density check
    U(2,i)   = (U(2,i)*dm(i) + uf(i)-uf(i+1))/dm_ip(i) ! Velocity in remapped cell
    U(3,i)   = (U(3,i)*dm(i) + Ef(i)-Ef(i+1))/dm_ip(i) ! Energy   in remapped cell
enddo

U(1,:) = 1./V(1,:)               ! Tau in remapped cell
V(2,:) = U(3,:) - .5*U(2,:)**(2) ! Internal energy in remapped cell
do i=iMIN,iMAX
    V(2,i) = MAX(0.000001,V(2,i)) ! Require positive internal energy
end do 
V(3,:) = (gm-1.)*V(1,:)*V(2,:)   ! Pressure in remappend cell

r_12_ip = r_12_i ! Update cell edges (that is they are the same as the stationary eulerian cells
dm = dm_ip          ! Update mass

end subroutine
