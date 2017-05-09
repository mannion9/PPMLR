subroutine RiemannSolver(u,p,d,us,ps,a,dx,dt,gm,igm,r)
use constants
implicit none
real,intent(in) :: dt
real,intent(in),dimension(imin:imax)      :: u,p,d,a,dx,gm,igm
real,intent(in),dimension(imin:imax+1)    :: r
real,intent(inout),dimension(jmin:jmax+1) :: us,ps
real,dimension(jmin:jmax+1) :: pL,pR,uL,uR,dL,dR,gmL,gmR,igmL,igmR  ! Left and right state of interfaces
! Determine effective left and right states for Riemann Solver
call RiemannEffectiveStates(u,p,d,gm,igm,uL,uR,pL,pR,dL,dR,gmL,gmR,igmL,igmR,dx,dt,a,r)
if (Solver.EQ.'HLL' )     call RP_HLL_roe(pL,pR,uL,uR,dL,dR,us,ps)
if (Solver.EQ.'HLLC')     call RP_HLLC_roe(pL,pR,uL,uR,dL,dR,us,ps)
if (Solver.EQ.'TwoShock') call TwoShock(ps,us,pL,pR,uL,uR,dL,dR,gmL,gmR,igmL,igmR,gm)
if (Solver.NE.'HLL'.AND.Solver.NE.'HLLC'.AND.Solver.NE.'TwoShock') print*,'Choose a correct Riemann Solver'
end subroutine

subroutine RP_HLL_roe(pL,pR,uL,uR,dL,dR,us,ps)
!!! ASSUMES IDEAL GAS FOR MOMENT
use Constants
implicit none
real,intent(in),dimension(jmin:jmax+1)    :: pL,pR,uL,uR,dL,dR  ! Left and right states of interfaces
real,intent(inout),dimension(jmin:jmax+1) :: ps,us
real :: ar,al,croe,gamma =5./3.
integer :: i
do i=jmin,jmax+1
	croe = SQRT(gamma*(pL(i)+pR(i))/(1./dR(i)+1./dL(i))) ! Eq. 3.26 Munz
	al   = MAX(delta,croe)  !  The only difference between left and right is a sign, but I took that into account in my equation for ps
	ar   = MAX(delta,croe)
	ps(i) = (pL(i)*ar+pR(i)*al-ar*al*(uR(i)-uL(i)))/(ar+al)
	ps(i) = MAX(delta,ps(i)) ! Requires positive pressure
	us(i) =  uL(i) - (ps(i)-pL(i))/al
end do
end subroutine

subroutine RP_HLLC_roe(pL,pR,uL,uR,dL,dR,us,ps)
! ASSUMES IDEAL GAS FOR MOMENT
use Constants
implicit none
real,intent(in),dimension(jmin:jmax+1)    :: pL,pR,uL,uR,dL,dR  ! Left and right states of interfaces
real,intent(inout),dimension(jmin:jmax+1) :: ps,us
real:: ar,al,croe,gamma =5./3.
integer :: i
do i=jmin,jmax+1
	croe=SQRT(gamma*(pL(i)+pR(i))/(1./dR(i)+1./dL(i))) ! Eq. 3.26 Munz
	al = MAX(delta,croe)  ! The only difference between left and right is a sign, but I took that into accound in my equation for ps
	ar = MAX(delta,croe)
	us(i) = ((pL(i)-pR(i))+uL(i)*al+uR(i)*ar)/(al+ar)
	ps(i) =  pR(i) + ar*(us(i)-uR(i))
	ps(i) = MAX(delta,ps(i))  ! Requires positive pressure
end do
end subroutine

subroutine TwoShock(ps,us,pL,pR,uL,uR,dL,dR,gmL,gmR,igmL,igmR,gm)
use Constants
implicit none
real,intent(in),dimension(imin:imax)      :: gm
real,intent(inout),dimension(jmin:jmax+1) ::  us,ps,pL,pR,uL,uR,dL,dR,gmL,gmR,igmL,igmR
real,dimension(jmin:jmax+1) :: cL,cR,wL,wR
real :: pmin=1.0E-4,tol=1.0E-6,gmavg,igmavg,pdf,pdfL,pdfR,udf,gmsL,gmsR,wdenL,wdenR,gmin,gmax,psi,psii,usLi,usRi,usLii,usRii
integer :: i,j,itmax=5
do j=jmin,jmax+1
	gmavg   = .5*(gmL(j)+gmR(j))												! CG Eq. 33
	igmavg  = .5*(igmL(j)+igmR(j))
	cL(j)   = SQRT(igmL(j)*pL(j)*dL(j))   										! CG Eq. 7 Lagrange "speed of sound" (dm/dt)
	cR(j)   = SQRT(igmR(j)*pR(j)*dR(j))
	!  Initial guess for ustar left and right. Note that the secant method requires to intial guesses for pstar
	usLi  = uL(j) ; usRi  = uR(j)
	usLii = cL(j) ; usRii = cR(j)
	psii  = .5*(pL(j)+pR(j))
	psi   = (cR(j)*pL(j)+cL(j)*pR(j)-cL(j)*cR(j)*(uR(j)-uL(j)))/(cL(j)+cR(j)) ! Initial guess assuming weak wave, Van Leer Eq (59)
	if (psi.LT.pmin) psi = pmin
	do i=1,itmax
		! Calculate effective gamma in star region
		gmsL   = gmL(j)+2.*(1.-gmavg/igmavg)*(gmavg-1.)*(psi-pL(j))/(psi+pL(j))	! CG Eq. 31
		gmsR   = gmR(j)+2.*(1.-gmavg/igmavg)*(gmavg-1.)*(psi-pR(j))/(psi+pR(j))
		gmin   = MIN(gm(j-1),gm(j),gm(j+1)) ; gmax = MAX(gm(j-1),gm(j),gm(j+1))
		gmsL   = MAX(gmin,MIN(gmsL,gmax))	; gmsR = MAX(gmin,MIN(gmsR,gmax))   ! CG Eq. 32
		! Calculate wave speeds
		pdfL   = psi-pL(j)					; pdfR = psi-pR(j)
		wdenL  = (psi-pL(j)*(gmsL-1.)/(gmL(j)-1.))					            ! CG Eq. 34 denominator
		wdenR  = (psi-pR(j)*(gmsR-1.)/(gmR(j)-1.))
		if (ABS(wdenL).LT.delta) wdenL = delta									! pathelogical case if wdenL=0.0
		if (ABS(wdenR).LT.delta) wdenR = delta
		wL(j)  = SQRT(ABS(dL(j)*pdfL*(psi+.5*(gmsL-1.)*(psi+pL(j)))/wdenL)) 	! CG Eq. 34
		wR(j)  = SQRT(ABS(dR(j)*pdfR*(psi+.5*(gmsR-1.)*(psi+pR(j)))/wdenR))
		if (ABS(pdfL).LT.delta) wL(j)=cL(j)										! CG Eq. 34 has a pathelogical case if ps-pL is small (weak wave) and should reduce to sound speed
		if (ABS(pdfR).LT.delta) Wr(j)=cR(j)
		! Secant Method
		usLi = uL(j) - pdfL/wL(j)		  							            ! CG Eq. 19
		usRi = uR(j) + pdfR/wR(j)
		pdf  = ABS(psi-psii)
		udf  = ABS(usLi-usLii)+ABS(usRi-usRii)
		if (pdf.LT.delta.OR.udf.LT.delta) then 									! CG Eq. 18 has two possible pathalogical cases
			pdf=0. ; udf=1.
		end if
		ps(j) = psi-(usRi-usLi)*pdf/udf											! CG Eq. 18
		if (ps(j).LT.pmin) ps(j)=pmin
		if (2.*(ps(j)-psi)/(ps(j)+psi).LT.tol) exit								! Convergence has been reached
		usLii = usLi   ; usRii = usRi
		psii  = psi
	end do
	wdenL  = (ps(j)-pL(j)*(gmsL-1.)/(gmL(j)-1.))					            ! CG Eq. 34 denominator
	wdenR  = (ps(j)-pR(j)*(gmsR-1.)/(gmR(j)-1.))
	if (ABS(wdenL).LT.delta) wdenL = delta										! pathelogical case if wdenL=0.0
	if (ABS(wdenR).LT.delta) wdenR = delta
	wL(j)  = SQRT(ABS(dL(j)*(ps(j)-pL(j))*(ps(j)+.5*(gmsL-1.)*(ps(j)+pL(j)))/wdenL))  ! CG Eq. 34
	wR(j)  = SQRT(ABS(dR(j)*(ps(j)-pR(j))*(ps(j)+.5*(gmsR-1.)*(ps(j)+pR(j)))/wdenR))
	if (ABS(ps(j)-pL(j)).LT.delta.OR.wL(j).GT.infty) wL(j)=cL(j)				! Limit of weak wave is a sound wave also there can be floating point error if ps=pL and very large
	if (ABS(ps(j)-pR(j)).LT.delta.OR.wR(j).GT.infty) Wr(j)=cR(j)
	us(j)  = (Wl(j)*uL(j)+Wr(j)*uR(j)-(pR(j)-pL(j)))/(Wl(j)+Wr(j)) 				! Van Leer Eq. 62
end do
end subroutine

subroutine TwoShock_Newton(ps,us,pL,pR,uL,uR,dL,dR)
! Old routine used when assuming and ideal EOS, solves two shock riemann problem by newton method
use Constants
implicit none
real,intent(inout),dimension(jmin:jmax+1) ::  us,ps,pL,pR,uL,uR,dL,dR
real,dimension(jmin:jmax+1) :: cL,cR,wL,wR,zL,zR,usL,usR
real :: gm = 5./3.,pmin=1.0E-4,gmf
integer :: i,j,itmax=2
gmf = (gm+1.)/(2.*gm)
do j=jmin,jmax+1
	cL(j) = SQRT(gm*pL(j)*dL(j))    ! Lagrange "speed of sound" (dm/dt)
	cR(j) = SQRT(gm*pR(j)*dR(j))
	ps(j) = (cR(j)*pL(j)+cL(j)*pR(j)-cL(j)*cR(j)*(uR(j)-uL(j)))/(cL(j)+cR(j)) ! Initial guess assuming weak wave, Van Leer Eq (59)
	if (ps(j).LT.pmin) ps(j) = pmin
	do i=1,itmax
		wL(j)  = cL(j)*SQRT(1.+gmf*(ps(j)/pL(j)-1.))  ! CW Eq. 2.8
		wR(j)  = cR(j)*SQRT(1.+gmf*(ps(j)/pR(j)-1.))
		usL(j) = uL(j) - (ps(j)-pL(j))/wL(j)		  ! Collela & Glaz Eq. 19
		usR(j) = uR(j) + (ps(j)-pR(j))/wR(j)
		zL(j)  = -4.*wL(j)**3/dL(j)/(4.*wL(j)**2/dL(j)-(gm+1.)*(ps(j)-pL(j)))	! PPM report Eq. 64
		zR(j)  =  4.*wR(j)**3/dR(j)/(4.*wR(j)**2/dR(j)-(gm+1.)*(ps(j)-pR(j)))
		ps(j)  = ps(j) + zL(j)*zR(j)*(usR(j)-usL(j))/(zR(j)-zL(j))				! PPM report Eq. 65
		if (ps(j).LT.pmin) ps(j)=pmin
	end do
	us(j) = usL(j) + (usR(j)-usL(j))*zR(j)/(zR(j)-zL(j)) ! PPM report Eq. 66
end do
end subroutine

subroutine TwoShock_Secant(ps,us,pL,pR,uL,uR,dL,dR)
! Old routine used when assuming and ideal EOS, solves two shock riemann problem by secant method
use Constants
implicit none
real,intent(inout),dimension(jmin:jmax+1) ::  us,ps,pL,pR,uL,uR,dL,dR
real,dimension(jmin:jmax+1) :: cL,cR,wL,wR,usL,usR,usLi,usRi
real :: gm = 5./3.,pmin=1.0E-4,gmf,pdf,udf,psi,psii
integer :: i,j,itmax=3
gmf = (gm+1.)/(2.*gm)
do j=jmin,jmax+1
	usLi(j) = uL(j)				! Initial guess for ustar left and right. Note that the secant method requires to intial guesses for pstar
	usRi(j) = uR(j)
	psii    = .5*(pL(j)+pR(j))
	cL(j) = SQRT(gm*pL(j)*dL(j))    ! Lagrange "speed of sound" (dm/dt)
	cR(j) = SQRT(gm*pR(j)*dR(j))
	psi = (cR(j)*pL(j)+cL(j)*pR(j)-cL(j)*cR(j)*(uR(j)-uL(j)))/(cL(j)+cR(j)) ! Initial guess assuming weak wave, Van Leer Eq (59)
	if (psi.LT.pmin) psi = pmin
	do i=1,itmax
		wL(j)  = cL(j)*SQRT(1.+gmf*(psi/pL(j)-1.))  ! CW Eq. 2.8
		wR(j)  = cR(j)*SQRT(1.+gmf*(psi/pR(j)-1.))
		usL(j) = uL(j) - (psi-pL(j))/wL(j)		  ! Collela & Glaz Eq. 19
		usR(j) = uR(j) + (psi-pR(j))/wR(j)
		pdf    = ABS(psi-psii)
		udf    = ABS(usL(j)-usLi(j))+ABS(usR(j)-usRi(j))
		if (pdf.LT.delta.OR.udf.LT.delta) then ! Two pathalogical cases
			pdf=0. ; udf=1.
		end if
		ps(j) = psi-(usR(j)-usL(j))*pdf/udf	! Collela & Glaz Eq. 18			! PPM report Eq. 65
		if (ps(j).LT.pmin) ps(j)=pmin
		usLi(j) = usL(j)
		usRi(j) = usR(j)
		psii    = psi
	end do
	wL(j)  = cL(j)*SQRT(1.+gmf*(ps(j)/pL(j)-1.))  ! CW Eq. 2.8
	wR(j)  = cR(j)*SQRT(1.+gmf*(ps(j)/pR(j)-1.))
	us(j)  = (Wl(j)*uL(j)+Wr(j)*uR(j)-(pR(j)-pL(j)))/(Wl(j)+Wr(j))
	print*,j,ps(j),us(j)
end do
end subroutine
