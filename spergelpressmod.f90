!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Spergel-Press WIMP heat transport module !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Contains the functions used in the Spergel Press section of transgen.f90. These are:
!	-nx_isothermal: Calculates the WIMP density in the Spergel-Press scheme
! 	-Etrans_sp: calculates the WIMP transported energy (eps_x) given the WIMP temperature (Tx)
!	-Tx_integral: to be used in newtons_meth
!	-newtons_meth: solves Tx_integral=0 which defines Tx

! All units are cgs except tab_r and tab_dr
! I apologize for the long function calls.

module spergelpressmod
use capmod
implicit none

double precision, parameter :: kB=1.38064852d-16, mnucg=1.6726219e-24

contains


function nx_isothermal(T_x, Nwimps)
implicit none
double precision, intent(in) :: T_x, Nwimps
double precision :: nx_isothermal(nlines)
double precision :: n_0, mxg
double precision :: R(nlines), phi(nlines)
! Calculates the isothermal wimp number density using eq. (2.25) in https://arxiv.org/pdf/0809.1871.pdf

r = tab_r*Rsun ! cm
phi = -tab_vesc**2/2.d0 ! erg/g
mxg = mdm*1.782662d-24
! open(95,"randphi.dat")
! write(95,*) r , phi
! close(95)
! WIMP number density in isothermal approximation
nx_isothermal = exp(-mxg*(phi-phi(1))/kB/T_x)
n_0 = Nwimps/trapz(r, 4.d0*pi*r**2.d0*nx_isothermal, nlines) ! Normalize so that integral(nx) = Nwimps
nx_isothermal = n_0*nx_isothermal

if (any(isnan(nx_isothermal))) print *, "NAN encountered in nx_isothermal"
return
end function


function Etrans_sp(T_x, sigma_N, Nwimps, niso)
implicit none
! Calculates WIMP transported energy (erg/g/s) using eq. (2.40) in https://arxiv.org/pdf/0809.1871.pdf

integer, intent(in) :: niso
double precision, intent(in) :: T_x, Nwimps
double precision, intent(in) :: sigma_N(niso)
double precision :: n_0, mxg
double precision :: R(nlines), phi(nlines), n_nuc(niso,nlines)
double precision :: n_x(nlines), species_indep(nlines), species_dep(nlines), sigma_nuc(niso)
double precision :: Etrans_sp(nlines)
integer :: i, j
! T_x in K, sigma_N in cm**2,

R = tab_r*Rsun ! R in cm
phi = -tab_vesc**2/2.d0 ! phi in erg/g
mxg = mdm*1.782662d-24 ! WIMP mass in g
! n_nuc in cm**-3
do i=1,niso
	n_nuc(i,:) = tab_mfr(:,i)*tab_starrho/AtomicNumber(i)/mnucg ! tab_starrho in gcm**-3
enddo
sigma_nuc = 2.d0*sigma_N ! Total WIMP-nucleus cross section in cm**2v. Only works for q/v independent cross-sections

! isothermal WIMP number density in cm**-3.
n_x = nx_isothermal(T_x, Nwimps)

! Separate calc into species dependent and independent factors
species_indep = 8.0d0*sqrt(2.d0/pi)*kB**(3.d0/2.d0)*n_x*(T_x-tab_T)/tab_starrho ! The species independent part

! Now sum over species to get the species dependent factor
species_dep=0.d0
do i=1,niso
	species_dep = species_dep + sigma_nuc(i)*n_nuc(i,:)*mxg*mnucg*AtomicNumber(i)/((mxg+mnucg*AtomicNumber(i))**2)* &
		sqrt(tab_T/(mnucg*AtomicNumber(i)) + T_x/mxg)
enddo

Etrans_sp = species_indep*species_dep ! erg/g/s


!! Useful when troubleshooting
!open(55, file="/home/luke/summer_2021/mesa/test_files/Etrans_sp_params.txt")
!write(55,*) "scalar params: T_x=", T_x, "m_x=", mxg, "m_nuc=", mnucg, "sigma_nuc=", sigma_nuc(1), &
!	"nlines=", nlines, "niso=", niso
!do i=1,nlines
!	write(55,*) R(i), tab_T(i), n_x(i), Etrans_sp(i) !n_x(i), tab_starrho(i), n_nuc(1,i), species_indep(i), phi(i)
!enddo
!close(55)

! print*, "************** Inside Etrans_sp ************** "
! print*, "tab_starrho: ", tab_starrho(253)
! print*, "mdm_g: ", mxg
! print*, "targetMass: ", mnucg
! print*, "nx: ", n_x(100)
! print*, "ndensity_target: ", n_nuc(1,253)
! print*, "Tx: ", T_x
! print*, "tab_T: ", tab_T(253)
! print*, "kb: ", kB

return
end function

!SB: Function to calculate the normalization constant of the SP number density function for a given Tx (eq 2.4)
function normalizationConst_mine(Tx, Nwimps)
	double precision, intent(in):: Tx, Nwimps
	double precision :: normalizationConst_mine
	double precision:: mxg
	integer:: i
	double precision:: phi(nlines),norm_integrand(nlines), R(nlines), phi2(nlines)
	phi = -tab_vesc**2/2.d0
	!phi = phi_mine()
	mxg = mdm*1.78d-24 !DM mass in grams
	R = tab_r*Rsun !Radius from core in [cm]
	norm_integrand = 4.d0*pi*R**2*exp(-mxg*(phi-phi(1))/kb/Tx)
	normalizationConst_mine = Nwimps/trapz(R, norm_integrand, nlines)
return
end function

!SB: Function to calculate the SP number density for a given Tx (eq 2.4)
function nxIso_mine(Tx, Nwimps)
	double precision::nxIso_mine(nlines)
	double precision, intent(in):: Tx, Nwimps
	double precision :: normConst
	double precision:: mxg, x
	integer:: i
	double precision:: phi(nlines)
	phi = -tab_vesc**2/2.d0
	!phi = phi_mine()
	mxg = mdm*1.78d-24

	normConst = normalizationConst_mine(Tx, Nwimps)
	nxIso_mine = normConst*exp(-mxg*(phi-phi(1))/kb/Tx)

return
end function


!SB: This is used to define Etrans as done in 2.10 (2111.0695)
!SB: and will be editted to include for general interactions for NREO
function Etrans_sp_mine(nq, nv, sigma_0, targetMass, electron_v_nucleons ,Tx, Nwimps,ndensity_target)
	implicit none
	double precision:: sigma_0, Nwimps, targetMass, Tx ! g
	double precision:: Etrans_sp_mine(nlines), nx(nlines), ndensity_target(nlines), Tcutoff(nlines)
	integer:: nq, nv, electron_v_nucleons, n
	double precision:: Afactor, Bfactor, Qfactor, mdm_g, q0_c0

	mdm_g = mdm*1.782662d-24 ![g]
	q0_c0 = 0.04*c0 ![GeV*cm/s]
	! q0_c0 =  mdm*v0 !0.04*c0 ![GeV*cm/s]
	n = nq+nv
		! targetMass = 1.6726219e-24
	!Determining the B_{2n_q} term

	if (nq.eq.1) then
		Bfactor = 8d0/3d0
	else if (nq.eq.2) then
		Bfactor = 4d0
	else if (nq.eq.-1) then
		Bfactor = 2d0
	else if (nq.eq.3) then
		Bfactor = 32d0/5d0
	end if
	!Determining the A_{2(nq+nv)} term
	if (n.eq.-1) then
		Afactor = 2d0
	else if (n.eq.0) then
		Afactor = 8d0
	else if (n.eq.1) then
		Afactor = 48d0
	else if (n.eq.2) then
		Afactor = 384d0
	end if

	if (nq.eq.0) then
		Qfactor = 2d0*sigma_0/v0**(2.*nv)
	else if (nq.ne.0) then
		Qfactor = Bfactor*2.**nq*(mdm/q0_c0)**(nq*2.)*sigma_0/(1.+mdm_g/targetMass)**(nq*2.)/v0**(2.*nv)
	end if

	nx = nx_isothermal(Tx, Nwimps)
	ETrans_sp_mine = Afactor/tab_starrho*sqrt(2./pi)*mdm_g*targetMass/(mdm_g+targetMass)**2&
										*nx*ndensity_target*Qfactor*kb*(Tx-tab_T)&
										*(kb*tab_T/targetMass+kb*Tx/mdm_g)**(0.5d0+nq+nv)

return
end function

!SB: This is used to define Etrans as done in 2.10 (2111.0695)
!SB: and will be editted to include for general interactions for NREO
function Etrans_sp_nreo(nq, nv, sigma_0, targetMass, electron_v_nucleons ,Tx, Nwimps, ndensity_target)
	implicit none
	double precision:: sigma_0, Nwimps, targetMass, Tx ! g
	double precision:: Etrans_sp_nreo(nlines), nx(nlines), ndensity_target(nlines), Tcutoff(nlines)
	integer:: nq, nv, electron_v_nucleons, n, i
	double precision:: Afactor, Bfactor, Qfactor, mdm_g

	mdm_g = mdm*1.782662d-24 ![g]
	n = nq+nv

	!Determining the B_{2n_q} term
	if (nq.eq.1) then
		Bfactor = 8./3.
	else if (nq.eq.2) then
		Bfactor = 4.
	else if (nq.eq.-1) then
		Bfactor = 2.
	else if (nq.eq.3) then
		Bfactor = 32./5.
	else if (nq.eq.4) then
		Bfactor = 2./3.
	end if


	if (nq.eq.0) then
		Qfactor = 2.*sigma_0/c0**(2.*nv)
	else if (nq.ne.0) then
		Qfactor = Bfactor*2.**nq*(mdm/c0)**(nq*2.)*sigma_0/(1.+mdm_g/targetMass)**(nq*2.)/c0**(2.*nv)
	end if

	nx = nxIso_mine(Tx, Nwimps)

	  if ((nq.eq.0).and.(nv.eq.0)) then
			ETrans_sp_nreo = (8*kb**(5./2.)*nx*Sqrt(2./pi)*Qfactor*(ndensity_target)* &
											  (Tx - (tab_T))*Sqrt((targetMass*(mdm_g)*(targetMass*Tx + &
											      (mdm_g)*(tab_T)))/kb**2))/((targetMass + (mdm_g))**2* &
											  (tab_starrho))
		! print*, "****************************************"
		! print*, "ETrans_sp_nreo(200): ", ETrans_sp_nreo(200)
		! print*, "nx(200): ", nx(200)
		! print*, "kb: ", kb
		! print*, "Qfactor: ", Qfactor
		! print*, "ndensity_target(200): ", ndensity_target(200)
		! print*, "Tx: ", Tx
		! print*, "tab_T(200): ", tab_T(200)
		! print*, "targetMass: ", targetMass
		! print*, "mdm_g: ", mdm_g
		! print*, "tab_starrho(200): ", tab_starrho(200)
		! print*, "q0: ", q0
		! print*, "sigma: ", sigma_0

		! open(unit=10, file='etrans_data_radial_debug.txt', status='replace', action='write')
		! 	 write(10, *) '# ETrans_sp | nx | kb | Qfactor | ndensity_target | Tx | tab_T | targetMass | mdm_g | rho'
		! do i=1, nlines
		! 		write(10, *) ETrans_sp_nreo(i), nx(i), kb, Qfactor, ndensity_target(i), Tx, tab_T(i), targetMass, mdm_g, tab_starrho(i)
		! end do
		!  close(10)

		else if ((nq.eq.0).and.(nv.eq.1)) then !v2
			ETrans_sp_nreo = (-48*kb**(5./2.)*nx*Sqrt(2/Pi)*Qfactor*(ndensity_target)* &
										  (-Tx + (tab_T))*(targetMass*Tx + (mdm_g)*(tab_T))**(3./2.))/ &
										 (Sqrt(targetMass*(mdm_g))*(targetMass + (mdm_g))**2* &
										  (tab_starrho))
		else if ((nq.eq.1).and.(nv.eq.0)) then !q2
			ETrans_sp_nreo = (-48.*kb**(5./2.)*nx*Sqrt(2/Pi)*Qfactor*(ndensity_target)* &
											  (-Tx + (tab_T))*(targetMass*Tx + (mdm_g)*(tab_T))**(3./2.))/ &
											 (Sqrt(targetMass*(mdm_g))*(targetMass + (mdm_g))**2.* &
											  (tab_starrho))

		else if ((nq.eq.1).and.(nv.eq.1)) then !v2 q2
			ETrans_sp_nreo = (-384.*kb**(7./2.)*nx*Sqrt(2/Pi)*Qfactor*(ndensity_target)* &
											  (-Tx + (tab_T))*(targetMass*Tx + (mdm_g)*(tab_T))**(5./2.))/ &
											 ((targetMass*(mdm_g))**(3./2.)*(targetMass + (mdm_g))**2* &
											  (tab_starrho))
		else if ((nq.eq.2).and.(nv.eq.0)) then !q4
			ETrans_sp_nreo = (-384.*kb**(7./2.)*nx*Sqrt(2./Pi)*Qfactor*(ndensity_target)* &
										  (-Tx + (tab_T))*(targetMass*Tx + (mdm_g)*(tab_T))**(5./2.))/ &
										 ((targetMass*(mdm_g))**(3./2.)*(targetMass + (mdm_g))**2* &
										  (tab_starrho))
		else if ((nq.eq.2).and.(nv.eq.1)) then !v2 q4
			 ETrans_sp_nreo = (-15.*nx*Sqrt(2/Pi)*Qfactor*(ndensity_target)*(kb/(targetMass*Tx + (mdm_g)*(tab_T)))**(9./2.)* &
										   (256.*targetMass**8.*Tx**8.*(-Tx + (tab_T)) +  &
											 2048.*targetMass**7.*Tx**7.*(mdm_g)*(tab_T)*(-Tx + (tab_T)) + &
										    7168.*targetMass**6.*Tx**6.*(mdm_g)**2.*(tab_T)**2.*(-Tx + (tab_T)) + &
												 14336.*targetMass**5.*Tx**5.*(mdm_g)**3.*(tab_T)**3.*(-Tx + (tab_T)) + &
										    17920.*targetMass**4.*Tx**4.*(mdm_g)**4.*(tab_T)**4.*(-Tx + (tab_T)) +  &
												7168.*targetMass**2.*Tx**2.*(mdm_g)**6.*(tab_T)**6.*(-Tx + (tab_T)) + &
										    2048.*targetMass*Tx*(mdm_g)**7.*(tab_T)**7.*(-Tx + (tab_T)) +  &
												256.*(mdm_g)**8.*(tab_T)**8.*(-Tx + (tab_T)) + &
										    7.*targetMass**3.*Tx**3.*(mdm_g)**5.*(tab_T)**5.*(-2039.*Tx +  &
												2048.*(tab_T))))/(tab_starrho*(targetMass*(mdm_g))**(5./2.)*(targetMass + (mdm_g))**2.)

		else if ((nq.eq.3).and.(nv.eq.0)) then !q6
				ETrans_sp_nreo = (-15.*nx*Sqrt(2/Pi)*Qfactor*(ndensity_target)*(kb/(targetMass*Tx + (mdm_g)*(tab_T)))**(9./2.)* &
										   (256.*targetMass**8.*Tx**8.*(-Tx + (tab_T)) + 2048.*targetMass**7.* &
											 Tx**7.*(mdm_g)*(tab_T)*(-Tx + (tab_T)) + &
										    7168.*targetMass**6.*Tx**6.*(mdm_g)**2.*(tab_T)**2.*(-Tx + (tab_T))  &
												+ 14336.*targetMass**5.*Tx**5.*(mdm_g)**3.*(tab_T)**3.*(-Tx + (tab_T)) + &
										    17920.*targetMass**4.*Tx**4.*(mdm_g)**4.*(tab_T)**4.*(-Tx + (tab_T)) +  &
												7168.*targetMass**2.*Tx**2.*(mdm_g)**6.*(tab_T)**6.*(-Tx + (tab_T)) + &
										    2048.*targetMass*Tx*(mdm_g)**7.*(tab_T)**7.*(-Tx + (tab_T)) +  &
												256.*(mdm_g)**8.*(tab_T)**8.*(-Tx + (tab_T)) + &
										    7.*targetMass**3.*Tx**3.*(mdm_g)**5.*(tab_T)**5.*(-2039.*Tx +  &
												2048.*(tab_T))))/((targetMass*(mdm_g))**(5./2.)*(targetMass + (mdm_g))**2.*(tab_starrho))

		else if ((nq.eq.3).and.(nv.eq.1)) then !v2 q6
	 			ETrans_sp_nreo = (-45*kb**(11./2.)*nx*Sqrt(2./Pi)*Qfactor*(ndensity_target)* &
												(1024.*targetMass**9.*Tx**9.*(-Tx + (tab_T)) + &
										    9216.*targetMass**8.*Tx**8.*(mdm_g)*(tab_T)*(-Tx + (tab_T)) +  &
												36864.*targetMass**7.*Tx**7.*(mdm_g)**2.*(tab_T)**2.*(-Tx + (tab_T)) + &
										    86016.*targetMass**6.*Tx**6.*(mdm_g)**3.*(tab_T)**3.*(-Tx + (tab_T)) +  &
												129024.*targetMass**5.*Tx**5.*(mdm_g)**4.*(tab_T)**4.*(-Tx + (tab_T)) + &
										    128772.*targetMass**4.*Tx**4.*(mdm_g)**5.*(tab_T)**5.*(-Tx + (tab_T)) +  &
												36864.*targetMass**2.*Tx**2.*(mdm_g)**7.*(tab_T)**7.*(-Tx + (tab_T)) + &
										    9216.*targetMass*Tx*(mdm_g)**8.*(tab_T)**8.*(-Tx + (tab_T)) +  &
												1024.*(mdm_g)**9.*(tab_T)**9.*(-Tx + (tab_T)) + &
										    21.*targetMass**3.*Tx**3.*(mdm_g)**6.*(tab_T)**6.*(-4043.*Tx +  &
												4096.*(tab_T))))/((targetMass*(mdm_g))**(7./2.)*(targetMass + (mdm_g))**2.*(tab_starrho)*  &
										   (targetMass*Tx + (mdm_g)*(tab_T))**(9./2.))

		else if ((nq.eq.4).and.(nv.eq.0)) then !q8
				ETrans_sp_nreo = (-45*kb**(11/2)*nx*Sqrt(2/Pi)*Qfactor*(ndensity_target)*(1024*targetMass**9* &
												Tx**9*(-Tx + (tab_T)) + &
										    9216*targetMass**8*Tx**8*(mdm_g)*(tab_T)*(-Tx + (tab_T)) +  &
												36864*targetMass**7*Tx**7*(mdm_g)**2*(tab_T)**2*(-Tx + (tab_T)) + &
										    86016*targetMass**6*Tx**6*(mdm_g)**3*(tab_T)**3*(-Tx + (tab_T)) +  &
												129024*targetMass**5*Tx**5*(mdm_g)**4*(tab_T)**4*(-Tx + (tab_T)) + &
										    128772*targetMass**4*Tx**4*(mdm_g)**5*(tab_T)**5*(-Tx + (tab_T)) +  &
												36864*targetMass**2*Tx**2*(mdm_g)**7*(tab_T)**7*(-Tx + (tab_T)) + &
										    9216*targetMass*Tx*(mdm_g)**8*(tab_T)**8*(-Tx + (tab_T)) +  &
												1024*(mdm_g)**9*(tab_T)**9*(-Tx + (tab_T)) + &
										    21*targetMass**3*Tx**3*(mdm_g)**6*(tab_T)**6*(-4043*Tx +  &
												4096*(tab_T))))/((targetMass*(mdm_g))**(7/2)*(targetMass +  &
												(mdm_g))**2*(tab_starrho)* &
										   (targetMass*Tx + (mdm_g)*(tab_T))**(9/2))

		 else if ((nq.eq.4).and.(nv.eq.1)) then !v2 q8
 				 ETrans_sp_nreo = (-315.*kb**(13./2.)*nx*Sqrt(2./Pi)*Qfactor*(ndensity_target)* &
				 									(2048.*targetMass**10.*Tx**10.*(-Tx + (tab_T)) + &
											    20480.*targetMass**9.*Tx**9.*(mdm_g)*(tab_T)*(-Tx  &
													+ (tab_T)) + 92160.*targetMass**8.*Tx**8.*(mdm_g)**2.*(tab_T)**2.*(-Tx + (tab_T)) + &
											    245760.*targetMass**7.*Tx**7.*(mdm_g)**3.*(tab_T)**3.*(-Tx + (tab_T)) +  &
													430080.*targetMass**6.*Tx**6.*(mdm_g)**4.*(tab_T)**4.*(-Tx + (tab_T)) + &
											    515592.*targetMass**5.*Tx**5.*(mdm_g)**5.*(tab_T)**5.*(-Tx + (tab_T)) +  &
													427350.*targetMass**4.*Tx**4.*(mdm_g)**6.*(tab_T)**6.*(-Tx + (tab_T)) + &
											    92160.*targetMass**2.*Tx**2.*(mdm_g)**8.*(tab_T)**8.*(-Tx + (tab_T)) +  &
													20480.*targetMass*Tx*(mdm_g)**9.*(tab_T)**9.*(-Tx + (tab_T)) + &
											    2048.*(mdm_g)**10.*(tab_T)**10.*(-Tx + (tab_T)) +  &
													15.*targetMass**3.*Tx**3.*(mdm_g)**7.*(tab_T)**7.*(-15983.*Tx + 16384.*(tab_T))))/ &
											  ((targetMass + (mdm_g))**2.*(tab_starrho)*(targetMass*(mdm_g)*(targetMass*Tx + (mdm_g)*(tab_T)))**(9./2.))

		 else if ((nq.eq.5).and.(nv.eq.0)) then !q10
					ETrans_sp_nreo = (-315.*kb**(13./2.)*nx*Sqrt(2./Pi)*Qfactor*(ndensity_target)*&
														(2048.*targetMass**10.*Tx**10.*(-Tx + (tab_T)) +&
												    20480.*targetMass**9.*Tx**9.*(mdm_g)*(tab_T)*(-Tx + (tab_T)) +&
														 92160.*targetMass**8.*Tx**8.*(mdm_g)**2.*(tab_T)**2.*(-Tx + (tab_T)) + &
												    245760.*targetMass**7.*Tx**7.*(mdm_g)**3.*(tab_T)**3.*(-Tx + (tab_T)) + &
														430080.*targetMass**6.*Tx**6.*(mdm_g)**4.*(tab_T)**4.*(-Tx + (tab_T)) + &
												    515592.*targetMass**5.*Tx**5.*(mdm_g)**5.*(tab_T)**5.*(-Tx + (tab_T)) +  &
														427350.*targetMass**4.*Tx**4.*(mdm_g)**6.*(tab_T)**6.*(-Tx + (tab_T)) + &
												    92160.*targetMass**2.*Tx**2.*(mdm_g)**8.*(tab_T)**8.*(-Tx + (tab_T)) +  &
														20480.*targetMass*Tx*(mdm_g)**9.*(tab_T)**9.*(-Tx + (tab_T)) + &
												    2048.*(mdm_g)**10.*(tab_T)**10.*(-Tx + (tab_T)) +  &
														15.*targetMass**3.*Tx**3.*(mdm_g)**7.*(tab_T)**7.*(-15983.*Tx + 16384.*(tab_T))))/ &
												  ((targetMass + (mdm_g))**2.*(tab_starrho)*(targetMass*(mdm_g)*(targetMass*Tx + (mdm_g)*(tab_T)))**(9./2.))
			else if ((nq.eq.5).and.(nv.eq.1)) then !v2 q10
				ETrans_sp_nreo = (-315.*kb**(15./2.)*nx*Sqrt(2/Pi)	*Qfactor*(ndensity_target)* &
													(32768.*targetMass**11.*Tx**11.*(-Tx + (tab_T)) + &
											    360448.*targetMass**10.*Tx**10.*(mdm_g)*(tab_T)*(-Tx +  &
													(tab_T)) + 1802240.*targetMass**9.*Tx**9.*(mdm_g)**2.*(tab_T)**2.*(-Tx + (tab_T)) + &
											    5406720.*targetMass**8.*Tx**8.*(mdm_g)**3.*(tab_T)**3.*(-Tx + (tab_T)) +  &
													10813440.*targetMass**7.*Tx**7.*(mdm_g)**4.*(tab_T)**4.*(-Tx + (tab_T)) + &
											    15130752.*targetMass**6.*Tx**6.*(mdm_g)**5.*(tab_T)**5.*(-Tx + (tab_T)) +  &
													15087072.*targetMass**5.*Tx**5.*(mdm_g)**6.*(tab_T)**6.*(-Tx + (tab_T)) + &
											    10673520.*targetMass**4.*Tx**4.*(mdm_g)**7.*(tab_T)**7.*(-Tx + (tab_T)) +  &
													1802240.*targetMass**2.*Tx**2.*(mdm_g)**9.*(tab_T)**9.*(-Tx + (tab_T)) + &
											    360448.*targetMass*Tx*(mdm_g)**10.*(tab_T)**10.*(-Tx + (tab_T)) +  &
													32768.*(mdm_g)**11.*(tab_T)**11.*(-Tx + (tab_T)) + &
											    165.*targetMass**3.*Tx**3.*(mdm_g)**8.*(tab_T)**8.*(-31525.*Tx + &
													32768.*(tab_T))))/((targetMass*(mdm_g))**(11./2.)*(targetMass + (mdm_g))**2.*(tab_starrho)* &
											   (targetMass*Tx + (mdm_g)*(tab_T))**(9./2.))
			else if ((nq.eq.6).and.(nv.eq.0)) then !q12
				ETrans_sp_nreo = (-315.*kb**(15./2.)*nx*Sqrt(2./Pi)*Qfactor*(ndensity_target)*&
													(32768.*targetMass**11.*Tx**11.*(-Tx + (tab_T)) +&
											    360448.*targetMass**10.*Tx**10.*(mdm_g)*(tab_T)*(-Tx + (tab_T))&
													 + 1802240.*targetMass**9.*Tx**9.*(mdm_g)**2.*(tab_T)**2.*(-Tx + (tab_T)) + &
											    5406720.*targetMass**8.*Tx**8.*(mdm_g)**3.*(tab_T)**3.*(-Tx + (tab_T)) +  &
													10813440.*targetMass**7.*Tx**7.*(mdm_g)**4.*(tab_T)**4.*(-Tx + (tab_T)) + &
											    15130752.*targetMass**6.*Tx**6.*(mdm_g)**5.*(tab_T)**5.*(-Tx + (tab_T)) +  &
													15087072.*targetMass**5.*Tx**5.*(mdm_g)**6.*(tab_T)**6.*(-Tx + (tab_T)) + &
											    10673520.*targetMass**4.*Tx**4.*(mdm_g)**7.*(tab_T)**7.*(-Tx + (tab_T)) +  &
													1802240.*targetMass**2.*Tx**2.*(mdm_g)**9.*(tab_T)**9.*(-Tx + (tab_T)) + &
											    360448.*targetMass*Tx*(mdm_g)**10.*(tab_T)**10.*(-Tx + (tab_T)) +  &
													32768.*(mdm_g)**11.*(tab_T)**11.*(-Tx + (tab_T)) + &
											    165.*targetMass**3.*Tx**3.*(mdm_g)**8.*(tab_T)**8.*(-31525.*Tx +  &
													32768.*(tab_T))))/((targetMass*(mdm_g))**(11./2.)*(targetMass + (mdm_g))**2.*(tab_starrho)* &
											   (targetMass*Tx + (mdm_g)*(tab_T))**(9./2.))

			 else if ((nq.eq.6).and.(nv.eq.1)) then !v2 q12
	 			ETrans_sp_nreo = 	1d0 !DID NOT INTEGRATE ANALYTICALLY

			else if ((nq.eq.7).and.(nv.eq.0)) then !q14
				ETrans_sp_nreo = 	1d0 !DID NOT INTEGRATE ANALYTICALLY
			end if

return
end function

!SB: calculate the Tx_integral using the energy transfer from "Etrans_sp_mine" function
function Tx_integral_mine(T_x, sigma_0, targetMass, electron_v_nucleons, Nwimps, nabund)
	implicit none

	double precision, intent(in) :: T_x, Nwimps
	double precision, intent(in) :: sigma_0, targetMass
	integer, INTENT(IN):: electron_v_nucleons
	double precision :: R(nlines), integrand(nlines)
	double precision :: Tx_integral_mine, nabund(nlines)
	! integrand units: erg/cm/s
	R = tab_r*Rsun
	integrand = 4.d0*pi*R**2*tab_starrho*Etrans_sp_mine(nq, nv, sigma_0, targetMass, electron_v_nucleons ,T_x, Nwimps, nabund)

	! integral is Etrans_tot (erg/s)
	Tx_integral_mine = trapz(R, integrand, nlines)
	return
end function

!SB: Finds the Tx which satisfies the luminosity condition with the integral "Tx_integral_mine"
function binary_search_mine(f, sigma_0, targetMass, electron_v_nucleons, nabund, Nwimps, guess_1, guess_2, reltolerance)
	implicit none
	integer :: i
	double precision :: f ! Tx_integral_mine
	double precision, intent(in) :: Nwimps, reltolerance, guess_1, guess_2
	double precision, intent(in) :: sigma_0, targetMass,  nabund(nlines)
	integer, INTENT(IN):: electron_v_nucleons
	double precision :: x_1, x_2, x_3, f1, f2, f3, error
	double precision :: binary_search_mine

	! x_1 and x_2 are temperatures (K)
	x_1 = guess_1
	x_2 = guess_2
	error = reltolerance + 1.d0	! So that the first iteration is executed
	! Binary search loop
	i = 0
	do while (error > reltolerance)
		x_3 = (x_1 + x_2)/2.d0
		f1 = f(x_1, sigma_0, targetMass, electron_v_nucleons, Nwimps, nabund)
		f2 = f(x_2, sigma_0, targetMass, electron_v_nucleons, Nwimps, nabund)
		f3 = f(x_3, sigma_0, targetMass, electron_v_nucleons, Nwimps, nabund)
		if (f3 == 0.d0) then
			exit
		else if (f1*f3 .gt. 0) then ! if f1 and f3 have the same sign
			x_1 = x_3
		else if (f2*f3 .gt. 0) then
			x_2 = x_3
		endif
		error = abs(x_2-x_1)/x_2
		i = i + 1
	enddo
	binary_search_mine = x_3

	return
end function

!SB: Finds the Tx which satisfies the luminosity condition with the integral "Tx_integral_nreo"
!SB: This accounts for a linear combination for v-dep and q-dep terms (which occurs for NREO)
function binary_search_nreo(f, sigma_0, targetMass, mreduced, electron_v_nucleons, Nwimps, guess_1, &
				guess_2, reltolerance, vdep)
	implicit none
	integer :: i
	double precision :: f !Tx_integral_nreo
	double precision, intent(in) :: Nwimps, reltolerance, guess_1, guess_2
	double precision, intent(in) :: sigma_0, targetMass, mreduced
	integer, INTENT(IN):: electron_v_nucleons, vdep
	double precision :: x_1, x_2, x_3, f1, f2, f3, error
	double precision :: binary_search_nreo
	integer :: nv_here, nq_here
	! x_1 and x_2 are temperatures (K)
	x_1 = guess_1
	x_2 = guess_2
	error = reltolerance + 1.d0	! So that the first iteration is executed
	! Binary search loop
	i = 0
	do while (error > reltolerance)
		x_3 = (x_1 + x_2)/2.d0
		if (vdep.eq.0) then
			nv_here = 0
			nq_here = nq
			f1 = f(x_1, sigma_0, targetMass, electron_v_nucleons, Nwimps, nq_here, nv_here)
			f2 = f(x_2, sigma_0, targetMass, electron_v_nucleons, Nwimps, nq_here, nv_here)
			f3 = f(x_3, sigma_0, targetMass, electron_v_nucleons, Nwimps, nq_here, nv_here)
		else if (vdep.eq.1) then
			nv_here = 1
			nq_here = nq
			f1 = f(x_1, sigma_0, targetMass, electron_v_nucleons, Nwimps, nq_here, nv_here)/c0**2 &
						- 1/(2*mreduced)**2*f(x_1, sigma_0, targetMass, electron_v_nucleons, Nwimps, nq_here+1, 0)
			f2 = f(x_2, sigma_0, targetMass, electron_v_nucleons, Nwimps, nq_here, nv_here)/c0**2 &
						- 1/(2*mreduced)**2*f(x_2, sigma_0, targetMass, electron_v_nucleons, Nwimps, nq_here+1, 0)
			f3 = f(x_3, sigma_0, targetMass, electron_v_nucleons, Nwimps, nq_here, nv_here)/c0**2 &
						- 1/(2*mreduced)**2*f(x_3, sigma_0, targetMass, electron_v_nucleons, Nwimps, nq_here+1, 0)
		end if
		if (f3 == 0.d0) then
			exit
		else if (f1*f3 .gt. 0) then ! if f1 and f3 have the same sign
			x_1 = x_3
		else if (f2*f3 .gt. 0) then
			x_2 = x_3
		endif
		error = abs(x_2-x_1)/x_2
		i = i + 1
	enddo
	binary_search_nreo = x_3
	return
end function

!SB: calculate the Tx_integral using the energy transfer from "Etrans_sp_mine" function
function Tx_integral_nreo(T_x, sigma_0, targetMass, electron_v_nucleons, Nwimps, nq_in, nv_in, nabund)
	implicit none

	double precision, intent(in) :: T_x, Nwimps
	double precision, intent(in) :: sigma_0, targetMass
	integer, INTENT(IN):: electron_v_nucleons
	double precision :: R(nlines), integrand(nlines), nabund(nlines)
	double precision :: Tx_integral_nreo
	integer :: nq_in, nv_in
	! integrand units: erg/cm/s
	R = tab_r*Rsun
	integrand = 4.d0*pi*R**2*tab_starrho*Etrans_sp_nreo(nq_in, nv_in, sigma_0, targetMass, electron_v_nucleons ,T_x, Nwimps, nabund)

	! integral is Etrans_tot (erg/s)
	Tx_integral_nreo = trapz(R, integrand, nlines)
	return
end function

function Tx_integral(T_x, sigma_N, Nwimps, niso)
	implicit none
	! Calculates the Tx defining integral

	integer, intent(in) :: niso
	double precision, intent(in) :: T_x, Nwimps
	double precision, intent(in) :: sigma_N(niso)
	double precision :: R(nlines), integrand(nlines)
	double precision :: Tx_integral
	! integrand units: erg/cm/s
	R = tab_r*Rsun

	integrand = 4*pi*R**2*tab_starrho*Etrans_sp(T_x, sigma_N, Nwimps, niso)

	! integral is Etrans_tot (erg/s)
	Tx_integral = trapz(R, integrand, nlines)

	return
end function


function newtons_meth(f, sigma_N, Nwimps, niso, guess_1, guess_2, reltolerance)
	! Performs Newton's method to solve Tx_integral(T_x)=0 for T_x (the function returns T_x)
	! The parameter f is the Tx_integral function
	implicit none

	integer, intent(in) :: niso
	double precision :: f ! Tx_integral
	double precision, intent(in) :: Nwimps, reltolerance, guess_1, guess_2
	double precision, intent(in) :: sigma_N(niso)
	double precision :: x_1, x_2, x_3, f1, f2, error
	double precision :: newtons_meth
	! m_x and m_p in grams, T_x, T_star in Kelvin, sigma in cm**2, n_nuc, n_x in cm**-3

	! x_1 and x_2 are temperatures (K)
	x_1 = guess_1
	x_2 = guess_2
	error = reltolerance + 1	! So that the first iteration is executed

	! Newton's method loop
	do while (error > reltolerance)
		! Update x_3 using Newton's method formula
		f1 = f(x_1, sigma_N, Nwimps, niso)
		f2 = f(x_2, sigma_N, Nwimps, niso)
		x_3 = x_2 - f2*(x_2-x_1)/(f2 - f1)
		error = abs(x_3-x_2)/x_2
		x_1 = x_2
		x_2 = x_3
	enddo

	newtons_meth = x_3 ! The solution to the nonlinear equation

	return
end function

function binary_search(f, sigma_N, Nwimps, niso, guess_1, guess_2, reltolerance)
	integer, intent(in) :: niso
	integer :: i
	double precision :: f ! Tx_integral
	double precision, intent(in) :: Nwimps, reltolerance, guess_1, guess_2
	double precision, intent(in) :: sigma_N(niso)
	double precision :: x_1, x_2, x_3, f1, f2, f3, error
	double precision :: binary_search

	! x_1 and x_2 are temperatures (K)
	x_1 = guess_1
	x_2 = guess_2
	error = reltolerance + 1.d0	! So that the first iteration is executed

	! Binary search loop
	i = 0
	do while (error > reltolerance)
		x_3 = (x_1 + x_2)/2.d0
		f1 = f(x_1, sigma_N, Nwimps, niso)
		f2 = f(x_2, sigma_N, Nwimps, niso)
		f3 = f(x_3, sigma_N, Nwimps, niso)
		if (f3 == 0.d0) then
			exit
		else if (f1*f3 .gt. 0) then ! if f1 and f3 have the same sign
			x_1 = x_3
		else if (f2*f3 .gt. 0) then
			x_2 = x_3
		endif
		error = abs(x_2-x_1)/x_2
		i = i + 1
	enddo

	binary_search = x_3

	return
end function

subroutine fourier_smooth(x, y, x_even, y_even, cutoff, noise_indicator, nlines, lensav, ierr)
	! Cuts out the high frequency components of y. E.g. if cutoff=0.05, the top 95% of frequency components are cut
	! Also returns a "noise indicator" - The sum of the frequency components above the cutoff
	integer, intent(in) :: nlines, lensav
	integer :: ierr, i
	double precision, intent(in) :: x(nlines), x_even(nlines), cutoff
	double precision, intent(inout) :: y(nlines)
	double precision, intent(out) :: noise_indicator
	double precision :: y_even(nlines), work(nlines), wsave(lensav), bcoeff(nlines), ccoeff(nlines), dcoeff(nlines)
	double precision :: ispline, denominator

	! Make evenly spaced y array
	call spline(x, y, bcoeff, ccoeff, dcoeff, nlines)
	do i=1,nlines
		y_even(i) = ispline(x_even(i), x, y, bcoeff, ccoeff, dcoeff, nlines)
	enddo

	! Compute FFT of y
	call dfft1i (nlines, wsave, lensav, ierr)  !Initialize (required by fftpack)
	if (ierr /= 0) print *, "FFT initializer 'dfft1i' failed with error ", ierr

	call dfft1f(nlines, 1, y_even, nlines, wsave, lensav, work, nlines, ierr) ! Take FFT
	if (ierr /= 0) print *, "Forward FFT calculator 'dfft1f' failed with error ", ierr
	! dTdr_even is now the array Fourier components of dTdr_even (the way fftpack works)

	noise_indicator = 0.d0
	! Take the ratio of high frequency components to low frequency components as a measure of how noisy the data is
	do i=int(cutoff*nlines),nlines
		noise_indicator = noise_indicator + abs(y_even(i))
	enddo
	do i=1,int(cutoff*nlines)
		denominator = denominator + abs(y_even(i))
	enddo
	noise_indicator = noise_indicator/denominator

	! Cut out top 100*(1-cutoff)% of Fourier components
	do i=1,nlines
		if (i > int(cutoff*nlines)) then
			y_even(i) = 0.d0
		endif
	enddo

	! Rebuild y with high frequency components cut out
	call dfft1b(nlines, 1, y_even, nlines, wsave, lensav, work, nlines, ierr)
	if (ierr /= 0) print *, "Backward FFT calculator 'dfft1b' failed with error ", ierr

	! Evaluate y on original grid (ie go convert y_even --> y)
	call spline(x_even, y_even, bcoeff, ccoeff, dcoeff, nlines)
	do i=1,nlines
		y(i) = ispline(x(i), x_even, y_even, bcoeff, ccoeff, dcoeff, nlines)
	enddo
end subroutine

function rolling_avg(y, nlines)
	! Takes a 1D array f of length N, returns an array of length N whose ith entry
	! is the average of f(i) and its 4 nearest neighbours
	integer, intent(in) ::  nlines
	double precision, intent(in) :: y(nlines)
	integer :: i, j
	double precision :: rolling_avg(nlines)

	do i=5,nlines-4
	    rolling_avg(i) = 0.d0
	    do j=-4,4
	        rolling_avg(i) = rolling_avg(i) + y(i+j)
	    enddo
	    rolling_avg(i) = rolling_avg(i)/9.d0
	enddo
	! do boundary values manually
	rolling_avg(1) = (y(1)+y(2)+y(3))/3.d0
	rolling_avg(2) = (y(1)+y(2)+y(3)+y(4))/4.d0
	rolling_avg(3) = (y(1)+y(2)+y(3)+y(4)+y(5)+y(6))/6.d0
	rolling_avg(4) = (y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)+y(8))/8.d0
	rolling_avg(nlines-1) = (y(nlines-7)+y(nlines-6)+y(nlines-5)+y(nlines-4) &
								+y(nlines-3)+y(nlines-2)+y(nlines-1)+y(nlines))/8.d0
	rolling_avg(nlines-1) = (y(nlines-5)+y(nlines-4)+y(nlines-3)+y(nlines-2)+y(nlines-1)+y(nlines))/6.d0
	rolling_avg(nlines-1) = (y(nlines-3)+y(nlines-2)+y(nlines-1)+y(nlines))/4.d0
	rolling_avg(nlines) = (y(nlines-2)+y(nlines-1)+y(nlines))/3.d0

	return
end function

end module spergelpressmod
