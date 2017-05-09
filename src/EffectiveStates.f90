subroutine RiemannEffectiveStates(u,p,d,gm,igm,up,um,pp,pm,dp,dm,gmp,gmm,igmp,igmm,dx,dt,a,r)
! Following CW convention, our interpolation quantities are aL,aR,dela,a6, which are all properties of each cell.
! These should not be confused with ap,am, which are the effective right and left states for each interface.
use constants
implicit none
real,intent(in) :: dt
real,intent(in),dimension(imin:imax+1)   :: r
real,intent(in),dimension(imin:imax)     :: p,u,d,dx,a,gm,igm
real,intent(inout),dimension(jmin:jmax+1):: up,um,pp,pm,dp,dm,gmp,gmm,igmp,igmm    ! Left and right states of interfaces
real,dimension(jmin-1:jmax+1) :: pL,pR,p6,delp,& ! p for pressure  ! Interpolated in cells [jmin-1,jmax+1] because we need ghost cell average
								 uL,uR,u6,delu,& ! u for velocity
								 dL,dR,d6,deld   ! d for density
real,dimension(jmin-1:jmax+2) :: p_12,u_12,d_12	 ! Interface interpolated values
real,dimension(imin:imax)     :: z,Ageom,dmc
integer :: i


do i=imin,imax
	if(alpha.EQ.0) Ageom(i) = 1.									! CW Eq 2.7
	if(alpha.EQ.1) Ageom(i) = .5*(r(i+1)+r(i))
	if(alpha.EQ.2) Ageom(i) = (r(i+1)**2+r(i+1)*r(i)+r(i)**2)/3.
end do
dmc = d*dx
if (Method=='Godunov') then
	up(jmin:jmax+1) = u(jmin-1:jmax)
	pp(jmin:jmax+1) = p(jmin-1:jmax)
	dp(jmin:jmax+1) = d(jmin-1:jmax)
	um(jmin:jmax+1) = u(jmin:jmax+1)
	pm(jmin:jmax+1) = p(jmin:jmax+1)
	dm(jmin:jmax+1) = d(jmin:jmax+1)

else if (Method=='MUSCL') then
	if (Eulerian.EQ.'On') then
		call PPM_interpol(p,dx,pL,pR,delp,p6,p_12) ! Note that in the pure lagrangian+remp scheme we interpolate in volume coordinate
		call PPM_interpol(u,dx,uL,uR,delu,u6,u_12)
		call PPM_interpol(d,dx,dL,dR,deld,d6,d_12)
	else
		call PPM_interpol(p,dmc,pL,pR,delp,p6,p_12) 	! Note that in the pure lagrangian scheme we interpolate in mass coordinate
		call PPM_interpol(u,dmc,uL,uR,delu,u6,u_12)
		call PPM_interpol(d,dmc,dL,dR,deld,d6,d_12)
	end if
	z = dt*a*Ageom/dx  ! Note this is equivilant to dt*C/dm where C is the lagrange speed of sound
	pp(jmin:jmax+1) = pR(jmin-1:jmax) - .5*z(jmin-1:jmax)*(delp(jmin-1:jmax)-p6(jmin-1:jmax))
	up(jmin:jmax+1) = uR(jmin-1:jmax) - .5*z(jmin-1:jmax)*(delu(jmin-1:jmax)-u6(jmin-1:jmax))
	dp(jmin:jmax+1) = dR(jmin-1:jmax) - .5*z(jmin-1:jmax)*(deld(jmin-1:jmax)-d6(jmin-1:jmax))
	pm(jmin:jmax+1) = pL(jmin:jmax+1) + .5*z(jmin:jmax+1)*(delp(jmin:jmax+1)+p6(jmin:jmax+1))
	um(jmin:jmax+1) = uL(jmin:jmax+1) + .5*z(jmin:jmax+1)*(delu(jmin:jmax+1)+u6(jmin:jmax+1))
	dm(jmin:jmax+1) = dL(jmin:jmax+1) + .5*z(jmin:jmax+1)*(deld(jmin:jmax+1)+d6(jmin:jmax+1))

else if (Method=='PPM') then
	if (Eulerian.EQ.'On') then
		call PPM_interpol(p,dx,pL,pR,delp,p6,p_12) ! Note that in the pure lagrangian+remp scheme we interpolate in volume coordinate
		call PPM_interpol(u,dx,uL,uR,delu,u6,u_12)
		call PPM_interpol(d,dx,dL,dR,deld,d6,d_12)
	else
		call PPM_interpol(p,dmc,pL,pR,delp,p6,p_12) ! Note that in the pure lagrangian scheme we interpolate in mass coordinate
		call PPM_interpol(u,dmc,uL,uR,delu,u6,u_12)
		call PPM_interpol(d,dmc,dL,dR,deld,d6,d_12)
	end if
	z = dt*a*Ageom/dx  ! Note this is equivilant to dt*C/dm where C is the lagrange speed of sound
	pp(jmin:jmax+1) = pR(jmin-1:jmax) - .5*z(jmin-1:jmax)*(delp(jmin-1:jmax)-(1.-2.*z(jmin-1:jmax)/3.)*p6(jmin-1:jmax))
	up(jmin:jmax+1) = uR(jmin-1:jmax) - .5*z(jmin-1:jmax)*(delu(jmin-1:jmax)-(1.-2.*z(jmin-1:jmax)/3.)*u6(jmin-1:jmax))
	dp(jmin:jmax+1) = dR(jmin-1:jmax) - .5*z(jmin-1:jmax)*(deld(jmin-1:jmax)-(1.-2.*z(jmin-1:jmax)/3.)*d6(jmin-1:jmax))
	pm(jmin:jmax+1) = pL(jmin:jmax+1) + .5*z(jmin:jmax+1)*(delp(jmin:jmax+1)+(1.-2.*z(jmin:jmax+1)/3.)*p6(jmin:jmax+1))
	um(jmin:jmax+1) = uL(jmin:jmax+1) + .5*z(jmin:jmax+1)*(delu(jmin:jmax+1)+(1.-2.*z(jmin:jmax+1)/3.)*u6(jmin:jmax+1))
	dm(jmin:jmax+1) = dL(jmin:jmax+1) + .5*z(jmin:jmax+1)*(deld(jmin:jmax+1)+(1.-2.*z(jmin:jmax+1)/3.)*d6(jmin:jmax+1))
end if
gmp = gm(jmin-1:jmax)
gmm = gm(jmin:jmax+1)
igmp= igm(jmin-1:jmax)
igmm= igm(jmin:jmax+1)

end subroutine
