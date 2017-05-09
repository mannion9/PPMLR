subroutine Remap(u,p,d,h,e,us,dv,m,r,r0,dv0,gm,igm,a)
use constants
implicit none
real,intent(inout),dimension(imin:imax)   :: u,p,d,h,e,a,dv,dv0,m,gm,igm ! Updated values after lagrange step
real,intent(inout),dimension(imin:imax+1) :: r,r0													! r0 is the stationary grid
real,intent(in),dimension(jmin:jmax+1)    :: us
real,dimension(jmin-2:jmax+1) :: p_12,u_12,d_12,h_12,e_12  ! Cells in which we need an interoplation value, for its right interface
real,dimension(imin:imax+1)   :: dfv,dfu,dfh,dfe,dm	 ! Interface fluxes
real,dimension(jmin-1:jmax+1) :: pL,pR,p6,delp  ,& ! p for pressure  ! Interpolated in cells [jmin-1,jmax+1] because we need ghost cell average
								uL,uR,u6,delu  ,& ! u for velocity
								dL,dR,d6,deld  ,& ! d for density
								hL,hR,h6,delh  ,& ! h for total energy (hamiltonian)
								eL,eR,e6,dele
real,dimension(imin:imax)     :: m_temp
integer :: i,j
real :: z

dfv = (r**ap1-r0**ap1)/REAL(ap1)  ! Difference after lagrange step

if (Remaper.EQ.'Godunov') then
	! Calculate fluxes
	do i=jmin,jmax+1
		if (us(i).GE.0.0) then ! Interface has moved to right
			j=i-1
			dm(i) = dfv(i)*d(j)
			dfu(i)= dm(i)*u(j)
			dfh(i)= dm(i)*h(j)
		else 					 ! Interface has moved the left
			j=i
			dm(i) = dfv(i)*d(j)
			dfu(i)= dm(i)*u(j)
			dfh(i)= dm(i)*h(j)
		end if
	end do
else if (Remaper.EQ.'MUSCL') then
	!Calcuate interpolation of d,u,E which require interpolation of pressure
	call PPM_interpol(p,dv,pL,pR,delp,p6,p_12)
	call PPM_interpol(u,dv,uL,uR,delu,u6,u_12)
	call PPM_interpol(d,dv,dL,dR,deld,d6,d_12)
	call PPM_interpol(e,dv,eL,eR,dele,e6,e_12)
	h_12 = p_12/((gm(jmin-2:jmax+1)-1.)*d_12)+.5*u_12**2
	hL   = h_12(jmin-2:jmax)
	hR   = h_12(jmin-1:jmax+1)
	call PPM_monotone(h,hR,hL)
	delh = hR-hL
	h6   = 6.*(h(jmin-1:jmax+1)-.5*(hL+hR))
	! Calculate fluxes
	do i=jmin,jmax+1
		if (dfv(i).GE.0.0) then ! Interface has moved to right
			j = i-1
			z      = dfv(i)/dv(j)
			dm(i)  = dfv(i)*(dR(j)-.5*z*(deld(j)-d6(j)))
			z 	   = dm(i)/(dv(j)*d(j))
			dfu(i) = dm(i)*(uR(j)-.5*z*(delu(j)-u6(j)))
			dfh(i) = dm(i)*(hR(j)-.5*z*(delh(j)-h6(j)))
			dfe(i) = dm(i)*(eR(j)-.5*z*(dele(j)-e6(j)))
		else  					! Interface has moved to left
			j = i
			z      = abs(dfv(i))/dv(j) ! Note CW equations assume y>0
			dm(i)  = dfv(i)*(dL(j)+.5*z*(deld(j)+d6(j)))
			z      = abs(dm(i))/(dv(j)*d(j))
			dfu(i) = dm(i)*(uL(j)+.5*z*(delu(j)+u6(j)))
			dfh(i) = dm(i)*(hL(j)+.5*z*(delh(j)+h6(j)))
			dfe(i) = dm(i)*(eL(j)+.5*z*(dele(j)+e6(j)))
		end if
	end do

else if (Remaper.EQ.'PPM') then
	! Calcuate interpolation of d,u,E which require interpolation of pressure
	call PPM_interpol(p,dv,pL,pR,delp,p6,p_12)
	call PPM_interpol(u,dv,uL,uR,delu,u6,u_12)
	call PPM_interpol(d,dv,dL,dR,deld,d6,d_12)
	call PPM_interpol(e,dv,eL,eR,dele,e6,e_12)
	h_12 = p_12/((gm(jmin-2:jmax+1)-1.)*d_12)+.5*u_12**2
	hL   = h_12(jmin-2:jmax)
	hR   = h_12(jmin-1:jmax+1)
	call PPM_monotone(h,hR,hL)
	delh = hR-hL
	h6   = 6.*(h(jmin-1:jmax+1)-.5*(hL+hR))
	! Calculate fluxes
	do i=jmin,jmax+1
		if (dfv(i).GE.0.0) then ! Interface has moved to right
			j = i-1
			z 	   = dfv(i)/dv(j)
			dm(i)  = dfv(i)*(dR(j)-.5*z*(deld(j)-(1.-2.*z/3.)*d6(j)))
			z      = dm(i)/(dv(j)*d(j))
			dfu(i) = dm(i)*(uR(j)-.5*z*(delu(j)-(1.-2.*z/3.)*u6(j)))
			dfh(i) = dm(i)*(hR(j)-.5*z*(delh(j)-(1.-2.*z/3.)*h6(j)))
			dfe(i) = dm(i)*(eR(j)-.5*z*(dele(j)-(1.-2.*z/3.)*e6(j)))
		else 					! Interface has moved to left
			j = i
			z      = abs(dfv(i))/dv(j) ! Note CW equations assume y>0
			dm(i)  = dfv(i)*(dL(j)+.5*z*(deld(j)+(1.-2.*z/3.)*d6(j)))
			z 	   = abs(dm(i))/(dv(j)*d(j))
			dfu(i) = dm(i)*(uL(j)+.5*z*(delu(j)+(1.-2.*z/3.)*u6(j)))
			dfh(i) = dm(i)*(hL(j)+.5*z*(delh(j)+(1.-2.*z/3.)*h6(j)))
			dfe(i) = dm(i)*(eL(j)+.5*z*(dele(j)+(1.-2.*z/3.)*e6(j)))
		end if
	end do
end if

do j=jmin,jmax
	m_temp(j) = m(j) + dm(j) - dm(j+1)
	d(j) = m_temp(j)/dv0(j)
	u(j) = (u(j)*m(j)+dfu(j)-dfu(j+1))/m_temp(j)
	h(j) = (h(j)*m(j)+dfh(j)-dfh(j+1))/m_temp(j)
	e(j) = MAX(delta,h(j)-.5*u(j)**2)
	! e(j) = (e(j)*m(j)+dfe(j)-dfe(j+1))/m_temm(j)
	! h(j) = e(j)+.5*u(j)**2
end do
call EOS_ideal_gas(d,e,p,a,gm,igm)
m(jmin:jmax) = m_temp(jmin:jmax)
r = r0
dv= dv0
end subroutine
