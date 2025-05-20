!!!!!! TRANSPORTER GENERAL !!!!!
!!! Asymmetric dark matter transport routine, check out https://arxiv.org/pdf/1311.2074.pdf
!!! The zetas in Eq. 31 should not be there
!!! for constant, q- and v- dependent cross sections
!!! Uses capture_mod from capgen.f90

!Input:
! sigma_0: DM scattering cross-section
! num_wimps: Total number of DM particles in the star. I know ADM is not WIMPs, stop complaining
! num_isotopes: number of isotopes: 1 = spin-dependent
! nq, nv: v^n, q^n numberwang
! is_spin_dep: spin dependence: 1 = spin-dependent scattering, 0 = spin-independent scattering
! transport_formalism: 1=Gould & Raffelt, 2=Spergel & Press, 3=rescaled Spergel & Press

!dm properties are set when you call init_sun.


!Output
!transported erg/g/s

subroutine transport_energy(sigma_0, num_wimps, num_isotopes, q_pow, v_pow, is_spin_dep, transport_formalism, temperature_dm, &
	noise_indicator, transported, total_luminosity)

! m_dm is stored in capture_mod
! temperature_dm is the output one-zone WIMP temp
use capture_mod
use alpha_kappa_mod
use spergel_press_mod
implicit none
!nlines might be redundant
integer, intent(in) :: transport_formalism
logical splinelog, DCT !for PCHIP
integer, intent(in):: num_isotopes, v_pow, q_pow, is_spin_dep
double precision, intent(in) :: sigma_0, num_wimps
double precision, intent(out) :: noise_indicator
integer, parameter :: decsize = 75 !this should be done a bit more carefully
integer i, j, ri, ierr
integer (kind=4) :: lensav
double precision :: epso,total_luminosity
double precision, parameter :: GN = 6.674d-8, kBeV=8.617e-5 ! kB and mnucg defined in spergelpressmod
double precision :: mxg, q0_cgs, rchi, Tc, rhoc, K, L, integrand
double precision :: capped, capture_maximum !this is the output
double precision :: sigma_SI, sigma_SD, a
double precision :: phi(nlines), Ltrans(nlines),transported(nlines),mfp(nlines),nabund(num_isotopes,nlines),sigma_N(num_isotopes)
double precision :: nxLTE(nlines), thermavg_sigma(nlines), zeta_v(nlines), zeta_q(nlines)
double precision :: nx(nlines),alphaofR(nlines),kappaofR(nlines),cumint(nlines),cumNx,nxIso(nlines),nxIso_func(nlines),cumNxIso
double precision :: r_even(nlines), T_even(nlines), dTdr_even(nlines), work(2*nlines) ! Evenly spaced arrays for Fourier smoothing
! More evenly spaced arrays for Fourier
double precision :: L_even(nlines), dLdr_even(nlines), Etrans_even(nlines), Etrans_test(nlines), Ltrans_cond(nlines)
double precision :: r_double(2*nlines), dTdr_mirror(2*nlines-2), r_mesa(1999)
double precision :: muarray(num_isotopes),alpha(num_isotopes),kappa(num_isotopes),dphidr(nlines),dTdr(nlines)
double precision :: fgoth, hgoth(nlines), ggoth_mesa(1999), ggoth(nlines), dLdR(nlines),isplined1,dLdRscratch(nlines)
double precision :: dggothdr_mesa(1999), dggothdr(nlines), test_array(2000)
double precision :: biggrid(nlines), bcoeff(nlines), ccoeff(nlines), dcoeff(nlines) ! for spline
integer lwk !for pchip
double precision :: pchipScratch(3*nlines) !for pchip
double precision, allocatable :: wsave(:)
double precision :: brcoeff(nlines), crcoeff(nlines), drcoeff(nlines) ! for spline
double precision :: bdcoeff(decsize), cdcoeff(decsize), ddcoeff(decsize) ! for spline
double precision :: smallgrid(decsize), smallR(decsize), smallT(decsize), smallL(decsize),smalldL(decsize),smalldT(decsize),ispline
double precision :: temperature_dm, guess_1, guess_2, reltolerance ! For the Spergel & Press scheme
double precision :: nK_0!(7) ! For the recalibrated Spergel & Press scheme
double precision :: T_eq_Tx_index, r_T, a1, b1, c1, a2, b2, c2, A_MC, x0_MC, sigma_MC, b_MC, chi_MC(nlines), g_MC(nlines)
double precision :: A_LTE, x0_LTE, sigma_LTE, b_LTE, Ltrans_LTE(nlines), chi_LTE(nlines), g_LTE(nlines), T_index_array(1)

lwk = 3*nlines !This is the length of pchipScratch. don't redefine this without also changing pchipScratch
epso = star_r(2)/10.d0 ! small number to prevent division by zero
! smallr = (/((i*1./dble(decsize-1)),i=1,decsize)/) - 1./dble(decsize-1)
smallgrid =  (/((i*1./dble(decsize-1)),i=1,decsize)/) - 1./dble(decsize-1) !(/i, i=1,decsize /)
biggrid =  (/((i*1./dble(nlines-1)),i=1,nlines)/) - 1./dble(nlines-1) !(/i, i=1,nlines/)



mxg = m_dm*1.78d-24
q0_cgs = q0*5.344d-14
Tc = star_temp(1)
rhoc = star_rho(1)
nq = q_pow
nv = v_pow

if (is_spin_dep == 1) then
  sigma_SD = sigma_0
  sigma_SI = 0.d0
  if (num_isotopes .ne. 1) then
  	print *, "Warning: transport_energy does not properly handle spin-dependent scattering on elements that aren't hydrogen."
  	print *, "For heat transport with spin-dependent cross sections, set num_isotopes=1."
  endif
else if (is_spin_dep == 0) then
  sigma_SD = 0.d0
  sigma_SI = sigma_0
end if


if (decsize .ge. nlines) stop "Major problem in transport_energy: your low-res size is larger than the original"
!Check if the stellar parameters have been allocated
if (.not. allocated(star_r)) stop "Error: stellar parameters not allocated in transport_energy"


!set up extra stellar arrays that we need
phi = - star_escape**2/2.d0
dphidr = -star_grav

! smooth T derivative
! some gymnastics are necessary, because the temperature is not smooth at all
! a simple spline -> derivative doesn't help. Fourier method (implemented below with fourier_smooth) works better

! DPCHEZ just outputs dTdr using a cubic spline (it doesn't do any smoothing)
splinelog = .false.
call DPCHEZ( nlines, star_r, star_temp, dTdr, SPLINElog, pchipScratch, LWK, IERR )
if (ierr .lt. 0) then
	print*, 'DPCHEZ interpolant failed with error ', IERR
	return
ENDIF

! Smooth dTdr with FFT
! First build evenly spaced r and dTdr arrays
do i=1,nlines
	r_even(i) =  i*1./dble(nlines)
enddo
lensav = nlines + int(log(real(nlines))) + 4 ! Minimum length required by fftpack
! Cut out high frequency components of dTdr. The subroutine fourier_smooth is located in spergelpressmod.f90
! Keep lowest 5% of components, delete top 95% of frequency components
call fourier_smooth(star_r, dTdr, r_even, dTdr_even, 0.05d0, noise_indicator, nlines, lensav, ierr)
dTdr = dTdr/radius_star

if (any(isnan(dTdr))) print *, "NAN encountered in dT/dr"

! calculate sigma_i, interpolate alpha_i, and kappa_i from tables
! nq and nv are set in darkInputs.txt for DarkMESA or in main.f90 for captngen_test.x
call read_alpha_kappa(nq,nv)
alphaofR(:) = 0.d0
kappaofR(:) = 0.d0
do i = 1,num_isotopes
  a = atomic_nums(i)
  !this is fine for SD as long as it's just hydrogen. Otherwise, spins must be added (use effective operator method)
  muarray(i) = m_dm/a/m_proton
  sigma_N(i) = a**2 * (sigma_SI*a**2 + sigma_SD) * (m_dm+m_proton)**2 / (m_dm+a*m_proton)**2
  nabund(i,:) = star_fractions(:,i)*star_rho(:)/a/mnucg
  !these shouldn't really be done every iteration, can fix later
  call interpolate_1d(dm_target_ratio, alpha_loaded, alpha_kappa_length,muarray(i),alpha(i))
  call interpolate_1d(dm_target_ratio, kappa_loaded, alpha_kappa_length,muarray(i),kappa(i))
end do

!need separate zeta factors for q- and v- dependent interactions
do i = 1,nlines
  zeta_q(i) = q0_cgs/(mxg*sqrt(2.d0*kB*star_temp(i)/mxg))
  zeta_v(i) = v0/(sqrt(2.d0*kB*star_temp(i)/mxg))
end do

! mean free path calcs for each nq,nv case here
! equations from 1311.2074 (eqns 69 to 74 on arxiv copy) with corrections from Hannah Banks
if (nq*nv .ne. 0) then
  stop "Oh no! nq and nv can't both be nonzero."
else if ((nq .eq. 0) .and. (nv .eq. 0)) then
  do i = 1,nlines
    mfp(i) = 1./sum(2.d0*sigma_N*nabund(:,i)) ! Since sigma_tot = 2*sigma_0 for v/q independent scattering
  end do
else if ((nq .eq. 1)) then
  do i = 1,nlines
    mfp(i) = 1./sum(6.*nabund(:,i)*sigma_N/(1.+muarray)/(zeta_q(i)**2))
  end do
else if ((nq .eq. 2)) then
  do i = 1,nlines
    mfp(i) = 1./sum(40.*nabund(:,i)*sigma_N/((1.+muarray)**2)/(zeta_q(i)**4))
  end do
else if ((nq .eq. -1)) then
  do i = 1,nlines
    mfp(i) = 1./sum(nabund(:,i)*sigma_N*(1.+muarray)*zeta_q(i)**2)
  end do
else if ((nv .eq. 1)) then
  do i = 1,nlines
    mfp(i) = 1./sum(nabund(:,i)*sigma_N*(1.+muarray)*3./(zeta_v(i)**2))
  end do
else if ((nv .eq. 2)) then
  do i = 1,nlines
    mfp(i) = 1./sum(2*nabund(:,i)*sigma_N*((1.+muarray)**2)*15./2./(zeta_v(i)**4))
  end do
else if ((nv .eq. -1)) then
  do i = 1,nlines
    mfp(i) = 1./sum(4*nabund(:,i)*sigma_N*zeta_v(i)**2/(1.+muarray))
  end do
end if

rchi = (3.*(kB*Tc)/(2.*pi*GN*rhoc*mxg))**.5;
K = mfp(1)/rchi;
print *, K, "= K"
!K = 1;

! this loop does a number of things: gets alpha and kappa averages (average over isotopes) to get alpha(r), kappa(r),
! and calculates nxLTE

cumint(1) = 0.d0
cumNx = 0.d0
do i = 1,nlines
  !get alpha & kappa averages
  alphaofR(i) = sum(alpha*sigma_N*nabund(:,i))/sum(sigma_N*nabund(:,i))
  if ((nq .eq. 0) .and. (nv .eq. 0)) then
    kappaofR(i) = mfp(i)*sum(2.d0*sigma_N*nabund(:,i)/kappa) ! Since sigma_tot = 2*sigma_0
  else if ((nq .eq. 1)) then
    kappaofR(i) = mfp(i)*sum(6.*nabund(:,i)*sigma_N/(1.+muarray)/(zeta_q(i)**2)/kappa)
  else if ((nq .eq. 2)) then
    kappaofR(i) = mfp(i)*sum(40.*nabund(:,i)*sigma_N/((1.+muarray)**2)/(zeta_q(i)**4)/kappa)
  else if ((nq .eq. -1)) then
    kappaofR(i) = mfp(i)*sum(nabund(:,i)*sigma_N*(1.+muarray)*zeta_q(i)**2/kappa)
  else if ((nv .eq. 1)) then
    kappaofR(i) = mfp(i)*sum(nabund(:,i)*sigma_N*(1.+muarray)*3./2./(zeta_v(i)**2)/kappa)
  else if ((nv .eq. 2)) then
    kappaofR(i) = mfp(i)*sum(nabund(:,i)*sigma_N*((1.+muarray)**2)*15./4./(zeta_v(i)**4)/kappa)
  else if ((nv .eq. -1)) then
    kappaofR(i) = mfp(i)*sum(nabund(:,i)*sigma_N*2*zeta_v(i)**2/(1.+muarray)/kappa)
  end if
  kappaofR(i) = 1./kappaofR(i)

  !perform the integral inside the exponent in nx
  integrand = (kB*alphaofR(i)*dTdr(i) + mxg*dphidr(i))/(kB*star_temp(i))

  if (i > 1) then
  	cumint(i) = cumint(i-1) + integrand*star_dr(i)*radius_star
  end if

  nxLTE(i) = (star_temp(i)/Tc)**(3./2.)*exp(-cumint(i))
  nxIso(i) = num_wimps*exp(-radius_star**2*star_r(i)**2/rchi**2)/(pi**(3./2.)*rchi**3) !normalized correctly

  cumNx = cumNx + 4.*pi*star_dr(i)*star_r(i)**2*nxLTE(i)*radius_star**3.

end do




nxLTE = nxLTE/cumNx*num_wimps !normalize density

! temperature_dm is the Spergel & Press one-zone WIMP temperature in K - calculate it here to use in nxIso
guess_1 = maxval(star_temp)*1.1d0 ! One-zone WIMP temp guesses in K.
guess_2 = maxval(star_temp)/10.d0
reltolerance = 1.0d-6


! newtons_method finds the one-zone wimp temp that gives 0 total transported energy in Spergel-Press scheme
temperature_dm = binary_search(luminosity_dm, sigma_N, num_wimps, num_isotopes, guess_1, guess_2, reltolerance) ! defined in spergelpressmod.f90
! Using Spergel-Press nxIso in Gould-Raffelt scheme gives numerical problems, but ideally we would use it.
!nxIso = isothermal_dm_num_density(temperature_dm, num_wimps) ! Defined in spergelpressmod.f90



!These are the interpolating functions used by G&R for transition to LTE regime
fgoth = 1./(1.+(K/.4)**2)
hgoth = ((star_r*radius_star - rchi)/rchi)**3 +1.
hgoth(1) = 0.d0 !some floating point shenanigans.

nx = fgoth*nxLTE + (1.-fgoth)*nxIso

! 3 options to calculate transported: 1: G&R, 2: S&P, 3: S&P rescaled
select case (transport_formalism)

	case (1) ! transport_formalism=1 -> use Gould & Raffelt
		print*, "GR"
		open(4, file = 'LtransGR.dat')

		Ltrans_LTE = 4.*pi*(star_r+epso)**2.*radius_star**2.*kappaofR*nx*mfp*sqrt(kB*star_temp/mxg)*kB*dTdr;
		Ltrans = fgoth*hgoth*Ltrans_LTE

		if (any(isnan(Ltrans))) print *, "NAN encountered in Ltrans"

		!get derivative of luminosity - also noisy. Fourier method doesn't work as well here
		! There is no smoothing currently implemented here
		splinelog = .false.
		call DPCHEZ( nlines, star_r, Ltrans, dLdR, SPLINElog, pchipScratch, LWK, IERR )
		if (ierr .lt. 0) then
			print*, 'DPCHEZ interpolant failed with error ', IERR
			return
		ENDIF
		dLdr = dLdr/radius_star
		write(4,*) Ltrans
		close(4)

		transported = 1./(4.*pi*(star_r+epso)**2*star_rho)*dLdR/radius_star**2

!		! Useful when troubleshooting
!		! Check Ltrans
!		open(55,file = "scalar_params_gr.dat")
!		write(55,*) fgoth, rchi, radius_star
!		close(55)
!		open(55,file = "etrans_gr.dat")
!		do i=1,nlines
!			write(55,*) star_r(i), transported(i), kappaofR(i), alphaofR(i), mfp(i), star_temp(i), dTdR(i), star_rho(i), nx(i), &
!			dphidr(i), Ltrans(i), dLdr(i), star_fractions(i,1), cumint(i), hgoth(i), phi(i), hgoth(i)
!		end do
!		close(55)


	! case (2) ! transport_formalism=2 -> use Gould & Raffelt rescaled to agree with MC simulations of a realistic star

	! 	skew gaussian rescaling is explained in MCrescaling.pdf

	! 	T_index_array = (minloc(star_temp-temperature_dm)) ! Just a stupid rank mismatch thing
	! 	T_eq_Tx_index = T_index_array(1)
	! 	r_T = star_r(T_eq_Tx_index)! dimensionless

	! 	a1 = -41.64d0
	! 	b1 = -3.26d0
	! 	c1 = 1.1481d0
	! 	a2 = 8.888d15
	! 	b2 = -70.3d0
	! 	c2 = 11.79d0

	! 	A_MC = a1*exp(-((log(K)-b1)/c1)**2.d0) + a2*exp(-((log(K)-b2)/c2)**2.d0)
	! 	x0_MC = 0.18d0
	! 	sigma_MC = 0.34d0
	! 	b_MC = -0.03658d0*K**(-1.818d0) - 3.227d0
	! 	chi_MC = (log10(star_r+epso/r_T) - x0_MC)/sigma_MC
	! 	g_MC = A_MC*exp(-chi_MC**2.d0/2.d0)*(1+erf(b_MC*chi_MC/sqrt(2.d0)))

	! 	A_LTE = 27.17d0*K
	! 	x0_LTE = 0.17d0
	! 	sigma_LTE = 0.35d0
	! 	b_MC = -4.35d0
	! 	chi_LTE = (log10(star_r+epso/r_T) - x0_LTE)/sigma_LTE
	! 	g_LTE = A_LTE*exp(-chi_LTE**2.d0/2.d0)*(1+erf(b_LTE*chi_LTE/sqrt(2.d0)))

	! 	Ltrans_LTE = 4.*pi*(star_r+epso)**2.*radius_star**2.*kappaofR*nxLTE*mfp*sqrt(kB*star_temp/mxg)*kB*dTdr
	! 	Ltrans = (g_MC/g_LTE)*Ltrans_LTE ! g_MC/g_LTE replaces fgoth*hgoth

	! 	if (any(isnan(Ltrans))) print *, "NAN encountered in Ltrans"

	! 	!get derivative of luminosity - also noisy. Fourier method doesn't work as well here
	! 	! There is no smoothing currently implemented here
	! 	splinelog = .false.
	! 	call DPCHEZ( nlines, star_r, Ltrans, dLdR, SPLINElog, pchipScratch, LWK, IERR )
	! 	if (ierr .lt. 0) then
	! 		print*, 'DPCHEZ interpolant failed with error ', IERR
	! 		return
	! 	ENDIF
	! 	dLdr = dLdr/radius_star

	! 	transported = 1./(4.*pi*(star_r+epso)**2*star_rho)*dLdR/radius_star**2

!		! Useful when troubleshooting
!		open(55,file = "scalar_params_gr_skew.dat")
!		write(55,*) fgoth, rchi, radius_star
!		close(55)
!		open(55,file = "etrans_gr_skew.dat")
!		do i=1,nlines
!			write(55,*) star_r(i), transported(i), kappaofR(i), alphaofR(i), mfp(i), star_temp(i), dTdR(i), star_rho(i), nx(i), &
!			dphidr(i), Ltrans(i), dLdr(i), star_fractions(i,1), cumint(i), hgoth(i), phi(i), g_MC(i), g_LTE(i), chi_MC(i), chi_LTE(i)
!		end do
!		close(55)


	case (2) ! transport_formalism=2 -> use Spergel & Press

		print*, "SP"
		! The Spergel-Press heat transport scheme: articles.adsabs.harvard.edu/pdf/1985ApJ...294..663S
		! The functions of interest are in spergelpressmod.f90. These also use https://arxiv.org/pdf/0809.1871.pdf

		! transported in erg/g/s (according to Spergel Press)
		transported = transport_energy_sp(temperature_dm, sigma_N, num_wimps, num_isotopes) ! erg/g/s

		!open a file to write Ltrans data to
		! open(5, file = 'LtransSP.dat')
		! Calculate Ltrans
		do i=1,nlines
			Ltrans(i) = trapezoid(star_r*radius_star, 4.d0*pi*(star_r*radius_star)**2.d0*transported*star_rho, i)
			! write(5,*) star_r(i), Ltrans(i)
		enddo

		!close the Ltrans data file
		close(5)

!		! useful when troubleshooting
!		open(55,file = "etrans_sp.dat")
!		do i=1,nlines
!			write(55,*) star_r(i), Ltrans(i), transported(i), nx(i), star_temp(i), star_grav(i), dTdr(i), nabund(1,i)
!		end do
!		close(55)
	case(3) ! transport_formalism=3 -> use rescaled Spergel & Press

		! The rescaled Spergel & Press transport scheme from Banks et. al https://arxiv.org/abs/2111.06895

		print*, "SP Recalculated"

    !there's a clever lookup way of doing this...
		! nK_0 = [0.40,0.21,1.05,1.72,0.11,0.73,1.20]
    if (nv .eq. -1) then
      nK_0 = 0.11
    else if (nv .eq. 1) then
      nK_0 = 0.73
    else if (nv .eq. 2) then
        nK_0 = 1.20
    else if (nq .eq. -1) then
        nK_0 = 0.21
    else if (nq .eq. 1) then
        nK_0 = 1.05
    else if (nq .eq. -1) then
        nK_0 = 1.72
    else
        nK_0 = 0.4
    end if



		! transported in erg/g/s (according to Spergel Press)
		transported = transport_energy_sp(temperature_dm, sigma_N, num_wimps, num_isotopes) ! erg/g/s
		! open(7, file = 'LtransNewSP.dat')

		do i=2,nlines
			Ltrans(i) = trapezoid(star_r*radius_star, 4.d0*pi*(star_r*radius_star)**2.d0*transported*star_rho, i)
      Ltrans(i) =  0.5*(1/(1+(nK_0/K)**2.))*Ltrans(i)
			! L = 0.5*(1/(1+(nK_0(j)/K)**2.))*Ltrans(i)
			! write(7,*) star_r(i), L
		enddo
		! close(7)

	case default

	stop "Invalid transport formalism. Should be an integer: 1, 2, or 3."

end select

! The total WIMP transported energy (erg/s). In the S&P scheme, this should be 0 by definition of temperature_dm.
total_luminosity = trapezoid(star_r*radius_star, 4.d0*pi*(star_r*radius_star)**2*transported*star_rho, nlines)
! total_luminosity = 1

! This is just to determine how noisy transported is. noise_indicator is the sum of frequency components above the cutoff
Etrans_test = transported
call fourier_smooth(star_r, Etrans_test, r_even, dTdr_even, 0.05d0, noise_indicator, nlines, lensav, ierr)

return

end subroutine transport_energy
