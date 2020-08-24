!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Nonlocal WIMP heat transport module !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Contains the functions used in the Spergel Press section of transgen.f90. These are:
! 	-Etrans_nl: calculates the WIMP transported energy (eps_x) given the WIMP temperature (Tx)
!	-Tx_integral: to be used in newtons_meth
!	-newtons_meth: solves Tx_integral=0 which defines Tx 

! I apologize for the horribly long function calls.

module nonlocalmod
use capmod
implicit none

contains


function nx_func(T_x, Tc, r, T, phi, rho_star, dr, dTdr, dphidr, mfp, n_nuc, sigma_N, sigma_0, alpha, kappa, &
	Nwimps, m_x, K, niso, nlines)
implicit none
integer, intent(in) :: nlines, niso
double precision, intent(in) :: T_x, Tc, Nwimps, m_x, sigma_0, K
double precision, intent(in) :: r(nlines), T(nlines), phi(nlines), rho_star(nlines), dr(nlines), dTdr(nlines)
double precision, intent(in) :: mfp(nlines), dphidr(nlines), n_nuc(niso, nlines)
double precision, intent(in) :: alpha(niso), kappa(niso), sigma_N(niso)
double precision :: n_0, fgoth, rchi, integrand, alphaofR(nlines), kappaofR(nlines), cumint(nlines), cumNx(nlines)
double precision, parameter :: pi=3.1415927d0, kB=1.38064852d-16, GN = 6.674d-8
double precision :: nx(nlines), nx_iso(nlines), nx_func(nlines)
integer :: i
! Calculates the isothermal wimp number density
! I know the function call is disgustingly long. I had to make nx a function in order to be able to use it
! in newton's method

!nx_iso = exp(-m_x*phi/kB/T_x) 
!n_0 = Nwimps/trapz(r, 4.d0*pi*r**2.d0*nx_iso, nlines) ! Normalize so that integral(nx) = Nwimps
!nx_iso = n_0*nx_iso

rchi = (3.*(kB*Tc)/(2.*pi*GN*rho_star(1)*m_x))**.5

cumint(1) = 0.d0
do i = 1,nlines

! 1) get alpha & kappa averages
  alphaofR(i) = sum(alpha*sigma_N*n_nuc(:,i))/sum(sigma_N*n_nuc(:,i))
  kappaofR(i) = mfp(i)*sum(sigma_0*sigma_N*n_nuc(:,i)/kappa)
  kappaofR(i) = 1./kappaofR(i)
  !perform the integral inside the nx integral

  integrand = (kB*alphaofR(i)*dTdr(i) + m_x*dphidr(i))/(kB*T(i))

  if (i > 1) then
  cumint(i) = cumint(i-1) + integrand*dr(i)
  end if

  nx(i) = (T(i)/Tc)**(3./2.)*exp(-cumint(i))

  cumNx = cumNx + 4.*pi*r(i)**2*dr(i)*nx(i)
  
  nx_iso(i) = Nwimps*exp(-r(i)**2/rchi**2)/(pi**(3./2.)*rchi**3) !normalized correctly
end do

nx = nx/cumNx*nwimps !normalize density
fgoth = 1./(1.+(K/.4)**2)

nx_func = fgoth*nx + (1.-fgoth)*nx_iso

if (any(isnan(alphaofR))) print *, "NAN encountered in alpha"
if (any(isnan(kappaofR))) print *, "NAN encountered in kappa"
if (any(isnan(cumint))) print *, "NAN encountered in cumint"
if (any(isnan(cumNx))) print *, "NAN encountered in cumNx"
if (any(isnan(nx_iso))) print *, "NAN encountered in nx_iso"
if (any(isnan(nx))) print *, "NAN encountered in nx"

return
end function


function Etrans_nl(T_x, Tc, r, T_star, phi, rho_star, dr, dTdr, dphidr, mfp, m_x, n_nuc, m_nuc, sigma_N, sigma_0, &
	alpha, kappa, Nwimps, K, nlines, niso)
implicit none
! Calculates WIMP Etrans using eq. (2.40) in https://arxiv.org/pdf/0809.1871.pdf

integer, intent(in) :: nlines, niso
double precision, intent(in) :: T_x, Tc, Nwimps, m_x, sigma_0, K
double precision, intent(in) :: r(nlines), T_star(nlines), rho_star(nlines), phi(nlines), dr(nlines), dTdr(nlines)
double precision, intent(in) :: dphidr(nlines), mfp(nlines), n_nuc(niso, nlines)
double precision, intent(in) :: alpha(niso), kappa(niso), sigma_N(niso), m_nuc(niso)
double precision :: n_0, fgoth, integrand, alphaofR(nlines), kappaofR(nlines), cumint(nlines), cumNx(nlines)
double precision, parameter :: pi=3.1415927d0, kB=1.38064852d-16, Rsun=6.957d10
double precision :: n_x(nlines), sigma_nuc(niso), species_indep(nlines), species_dep(nlines)
double precision :: Etrans_nl(nlines)
integer :: i
! m_x and m_p in grams, T_x, T_star in Kelvin, sigma in cm^2, n_nuc, n_x in cm^-3

sigma_nuc = 2.d0*sigma_N*sigma_0 ! Total WIMP-nucleus cross section

!! A better nxIso estimate than in Transgen
!n_x = exp(-m_x*phi/kB/T_x) 
!n_0 = Nwimps/trapz(r, 4.d0*pi*r**2.d0*n_x, nlines) ! Normalize so that integral(nx) = Nwimps
!n_x = n_0*n_x

n_x = nx_func(T_x, Tc, r, T_star, phi, rho_star, dr, dTdr, dphidr, mfp, n_nuc, sigma_N, sigma_0, alpha, kappa, &
	Nwimps, m_x, K, niso, nlines)

! Separate calc into species dependent and independent factors for convenience
species_indep = 8.0d0*sqrt(2.d0/pi)*kB**(3.d0/2.d0)*n_x*(T_x-T_star)/rho_star ! The species independent part

! Now sum over species to get the species dependent factor
species_dep=0.d0
species_dep = sigma_nuc(1)*n_nuc(1,:)*m_x*m_nuc(1)/((m_x+m_nuc(1))**2)*sqrt(T_star/m_nuc(1) + T_x/m_x)
!print *, "In Etrans_nl:", "m_nuc=", m_nuc(1), "mx=", m_x, "Tx=", T_x, "sigma_nuc", sigma_nuc(1), &
!	"rho_star=", maxval(rho_star), "n_nuc=", n_nuc(1,1)
!do i=1,niso
!	species_dep = species_dep + sigma_nuc(i)*n_nuc(i,:)*m_x*m_nuc(i)/((m_x+m_nuc(i))**2) * & 
!		sqrt(T_star/m_nuc(i) + T_x/m_x)
!enddo

Etrans_nl = species_indep*species_dep ! erg/g/s

open(55, file="/home/luke/summer_2020/mesa/test_files/Etrans_nl_params.txt")
write(55,*) "scalar params: T_x=", T_x, "m_x=", m_x, "m_nuc=", m_nuc(1), "sigma_nuc=", sigma_nuc(1), &
	"nlines=", nlines, "niso=", niso
do i=1,nlines
	write(55,*) r(i), T_star(i), n_x(i), rho_star(i), n_nuc(1,i), species_indep(i), species_dep(i), Etrans_nl(i)
enddo
close(55)

return
end function


function Tx_integral(T_x, Tc, r, T_star, phi, rho_star, dr, dTdr, dphidr, mfp, m_x, n_nuc, m_nuc, sigma_N, sigma_0, &
	alpha, kappa, Nwimps, K, nlines, niso)
implicit none
! Calculates the Tx defining integral 
! T_x, T_star in K, r in cm, phi in erg/g, rho_star in g/cm^3, m_x and m_nuc in g, n_nuc in cm^-3, sigma in cm^2

double precision, intent(in) :: T_x	! The independent variable. All others are params (for Newton's method)
integer, intent(in) :: nlines, niso
double precision, intent(in) :: Tc, Nwimps, m_x, sigma_0, K
double precision, intent(in) :: r(nlines), T_star(nlines), rho_star(nlines), phi(nlines), dr(nlines), dTdr(nlines)
double precision, intent(in) :: dphidr(nlines), mfp(nlines), n_nuc(niso, nlines)
double precision, intent(in) :: alpha(niso), kappa(niso), sigma_N(niso), m_nuc(niso)
double precision :: n_0, fgoth, alphaofR(nlines), kappaofR(nlines), cumint(nlines), cumNx(nlines)
double precision, parameter :: pi=3.1415927d0, kB=1.38064852d-16, Rsun=6.957d10
double precision :: n_x(nlines), integrand(nlines)
double precision :: Tx_integral

! integrand units: erg/cm/s
integrand = 4*pi*r**2*rho_star*Etrans_nl(T_x, Tc, r, T_star, phi, rho_star, dr, dTdr, dphidr, mfp, & 
	m_x, n_nuc, m_nuc, sigma_N, sigma_0, alpha, kappa, Nwimps, K, nlines, niso)

! integral is Etrans_tot (erg/s)
Tx_integral = trapz(r, integrand, nlines)

return
end function


function newtons_meth(f, Tc, r, T_star, phi, rho_star, dr, dTdr, dphidr, mfp, m_x, n_nuc, m_nuc, sigma_N, sigma_0, &
	alpha, kappa, Nwimps, K, nlines, niso, guess_1, guess_2, tolerance)
! Performs Newton's method to solve Tx_integral(T_x)=0 for T_x (the function returns T_x)
! The parameter f is the Tx_integral function
implicit none

double precision :: f ! Tx_integral
double precision, intent(in) :: tolerance, guess_1, guess_2
integer, intent(in) :: nlines, niso
integer :: i, half
double precision, intent(in) :: Tc, Nwimps, m_x, sigma_0, K
double precision, intent(in) :: r(nlines), T_star(nlines), rho_star(nlines), phi(nlines), dr(nlines), dTdr(nlines)
double precision, intent(in) :: dphidr(nlines), mfp(nlines), n_nuc(niso, nlines)
double precision, intent(in) :: alpha(niso), kappa(niso), sigma_N(niso), m_nuc(niso)
double precision :: n_0, fgoth, integrand, alphaofR(nlines), kappaofR(nlines), cumint(nlines), cumNx(nlines)
double precision, parameter :: pi=3.1415927d0, kB=1.38064852d-16, Rsun=6.957d10
double precision :: n_x(nlines), sigma_nuc(niso), species_indep(nlines), species_dep(nlines)
double precision :: x_1, x_2, x_3, error
double precision :: newtons_meth
! m_x and m_p in grams, T_x, T_star in Kelvin, sigma in cm^2, n_nuc, n_x in cm^-3

! x_1 and x_2 are temperatures (K)
x_1 = guess_1
x_2 = guess_2
error = tolerance + 1	! So that the first iteration is executed

! Newton's method loop
do while (error > tolerance)
	! Update x_3 using Newton's method formula
	x_3 = x_2 - f(x_2, Tc, r, T_star, phi, rho_star, dr, dTdr, dphidr, mfp, & 
	m_x, n_nuc, m_nuc, sigma_N, sigma_0, alpha, kappa, Nwimps, K, nlines, niso)*(x_2-x_1) &
		/(f(x_2, Tc, r, T_star, phi, rho_star, dr, dTdr, dphidr, mfp, & 
	m_x, n_nuc, m_nuc, sigma_N, sigma_0, alpha, kappa, Nwimps, K, nlines, niso) - &
		f(x_1, Tc, r, T_star, phi, rho_star, dr, dTdr, dphidr, mfp, & 
	m_x, n_nuc, m_nuc, sigma_N, sigma_0, alpha, kappa, Nwimps, K, nlines, niso))
	error = abs(x_3-x_2)
	x_1 = x_2
	x_2 = x_3
!	print*, "T_x=", x_3, "Etranstot=", trapz(r, 4.d0*pi*r**2*f(x_3, Tc, r, T_star, phi, rho_star, &
!	 dr, dTdr, dphidr, mfp, m_x, n_nuc, m_nuc, sigma_N, sigma_0, alpha, kappa, Nwimps, K, nlines, niso)*rho_star, &
!	 nlines), "error=", error, "tolerance=", tolerance
	open(1, file="/home/luke/summer_2020/mesa/test_files/Etranstot_newtons.dat", access="APPEND")
	write(1,*) x_3, trapz(r, 4.d0*pi*r**2*f(x_3, Tc, r, T_star, phi, rho_star, dr, dTdr, dphidr, &
	 mfp, m_x, n_nuc, m_nuc, sigma_N, sigma_0, alpha, kappa, Nwimps, K, nlines, niso)*rho_star, nlines)
enddo

newtons_meth = x_3 ! The solution to the nonlinear equation

return
end function

subroutine fourier_smooth(x, y, x_even, y_even, cutoff, noise_indicator, nlines, lensav, ierr)
! Cuts out the high frequency components of y
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

end module nonlocalmod
