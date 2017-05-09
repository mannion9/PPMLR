program main
use constants
implicit none
real :: dt=500. , time=.0 , usdt
real,dimension(imin:imax)   :: d,p,e,u,h,dv,dv0,m,rc,rc0,a,gm,igm ! Located at cell centers
real,dimension(imin:imax+1) :: r,r0,ri,geom  			          ! Located at all cell edges
real,dimension(jmin:jmax+1) :: us,ps  			                   ! Located at physical cell interfaces
integer :: i,j

call OpenFiles()
call InitDomain(r,dv,rc,r0,dv0,rc0)
call InitCond(d,e,h,u,m,dv,rc)
call EOS_ideal_gas(d,e,p,a,gm,igm)
call Output(d,e,h,p,u,m,time,r,rc,0,gm)

do i=0,itterMax-1
	call Boundary(d,u,p)
	call EOS_ideal_gas(d,e,p,a,gm,igm)
	call Timestep(u,r,dt,a)
	call RiemannSolver(u,p,d,us,ps,a,dv,dt,gm,igm,r)
	ri = r  ! Store position
	! Move ghost cells to prevent negative volumes in ghost cells
	r(imin  :jmin  ) = r(imin  :jmin  ) + dt*us(jmin  )
	r(jmax+2:imax+1) = r(jmax+2:imax+1) + dt*us(jmax+1)
	if(alpha.EQ.0) geom(jmin) = 1.										! CW Eq. 2.10 Abar
	if(alpha.EQ.1) geom(jmin) = r(jmin)+.5* dt*us(jmin)
	if(alpha.EQ.2) geom(jmin) = r(jmin)**2+ dt*us(jmin)*(r(jmin)+ dt*us(jmin)/3.)
	do j=jmin,jmax
		usdt  = us(j+1)*dt
		r(j+1)= r(j+1) + usdt
		if(alpha.EQ.0) geom(j+1) = 1.										! CW Eq. 2.10 Abar
		if(alpha.EQ.1) geom(j+1) = r(j+1)+.5*usdt
		if(alpha.EQ.2) geom(j+1) = r(j+1)**2+usdt*(r(j+1)+usdt/3.)
		dv(j)    = (r(j+1)**ap1-r(j)**ap1)/REAL(ap1)
		d(j)     = m(j)/dv(j)
		u(j)     = u(j) + (dt/m(j))*(ps(j)-ps(j+1))*.5*(geom(j+1)+geom(j))
		h(j)     = h(j) + (dt/m(j))*(us(j)*ps(j)*geom(j)-us(j+1)*ps(j+1)*geom(j))
		e(j)     = MAX(h(j) - .5*u(j)**2,delta)
	end do
	call EOS_ideal_gas(d,e,p,a,gm,igm)       ! Calculate new pressure
	rc   = .5*(r(imin:imax)+r(imin+1:imax+1))! Find new cell centers
	dv(imin:jmin-1) = abs(r(imin+1:jmin)**ap1-r(imin:jmin-1)**ap1)/ap1    ! Volume of ghost cells
	dv(jmax+1:imax) = abs(r(jmax+2:imax+1)**ap1-r(jmax+1:imax)**ap1)/ap1
	! dv(imin:jmin-1) = (r(imin+1:jmin)**ap1-r(imin:jmin-1)**ap1)/ap1    ! Volume of ghost cells
	! dv(jmax+1:imax) = (r(jmax+2:imax+1)**ap1-r(jmax+1:imax)**ap1)/ap1
	time = time+ dt
	if (Eulerian.EQ.'On')      call Remap(u,p,d,h,e,us,dv,m,r,r0,dv0,gm,igm,a)
	if (MOD(i,writeStep).EQ.0) call Output(d,e,h,p,u,m,time,r,rc,i,gm)
	if (MINVAL(dv(jmin:jmax)).LT..0) print*,'Negative volume in cell number',MINLOC(dv)
	if (dt.LT.dtmin)                 print*,'Time step has become smaller than dtmin'
	if (dt.LT.dtmin)                 exit
	if (time.GT.tmax)                exit
	if (MINVAL(dv(jmin:jmax)).LT..0) exit
end do

print*,'Terminated at step:',i,'at time:',time
call CloseFiles()
end program
