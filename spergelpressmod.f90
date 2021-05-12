!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Spergel-Press WIMP heat transport module !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Contains the functions used in the Spergel Press section of transgen.f90. These are:
!	-nx_func: Calculates the WIMP density profile using the fgoth interpolation
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


function nx_func(T_x, dTdr, mfp, sigma_N, alpha, Nwimps, niso)
implicit none
integer, intent(in) :: niso
double precision, intent(in) :: T_x, Nwimps
double precision, intent(in) :: dTdr(nlines), mfp(nlines)
double precision, intent(in) :: alpha(niso), sigma_N(niso)
double precision :: n_0, fgoth, Tc, K, rchi, integrand, mxg, alphaofR(nlines)
double precision :: R(nlines), phi(nlines), dR(nlines), dphidr(nlines), n_nuc(niso,nlines)
double precision :: nx_LTE(nlines), nx_iso(nlines), nx_func(nlines), cumint(nlines), cumNx(nlines)
integer :: i, j
! Calculates the isothermal wimp number density
! I haven't succesfully integrated this function into the code
! For some reason, it gave problems with Newton's method when I tried to do that (probably a bug)

Tc = tab_T(1) ! K
R = tab_r*Rsun ! cm
dR = tab_dr*Rsun ! cm
phi = -tab_vesc**2/2.d0 ! erg/g
dphidr = -tab_g ! erg/g/cm
mxg = mdm*1.782662d-24

! WIMP number density in isothermal approximation
nx_iso = exp(-mxg*phi/kB/T_x) 
n_0 = Nwimps/trapz(r, 4.d0*pi*r**2.d0*nx_iso, nlines) ! Normalize so that integral(nx) = Nwimps
nx_iso = n_0*nx_iso

do i=1,niso
	n_nuc(i,:) = tab_mfr(:,i)*tab_starrho/AtomicNumber(i)/mnucg
enddo

rchi = (3.*(kB*Tc)/(2.*pi*GNewt*tab_starrho(1)*mxg))**.5 ! cm
K = mfp(1)/rchi

cumint(1) = 0.d0
do i = 1,nlines

! 1) get alpha average
  alphaofR(i) = sum(alpha*2.d0*sigma_N*n_nuc(:,i))/sum(2.d0*sigma_N*n_nuc(:,i))
  
  !perform the integral inside the nx integral
  integrand = (kB*alphaofR(i)*dTdr(i) + mxg*dphidr(i))/(kB*tab_T(i))

  if (i > 1) then
  cumint(i) = cumint(i-1) + integrand*dR(i)
  end if

  nx_LTE(i) = (tab_T(i)/Tc)**(3./2.)*exp(-cumint(i))

  cumNx = cumNx + 4.*pi*R(i)**2*dR(i)*nx_LTE(i)

end do

nx_LTE = nx_LTE/cumNx*nwimps !normalize density
fgoth = 1./(1.+(K/.4)**2)

nx_func = fgoth*nx_LTE + (1.-fgoth)*nx_iso

if (any(isnan(dphidr))) print *, "NAN encountered in dphidr"
if (any(isnan(alphaofR))) print *, "NAN encountered in alphaofR"
if (any(isnan(dTdr))) print *, "NAN encountered in dTdr"

! If the units are wrong, a NaN will likely show up in one of these
if (any(isnan(cumint))) print *, "NAN encountered in cumint"
if (any(isnan(nx_iso))) print *, "NAN encountered in nx_iso"
if (any(isnan(nx_LTE))) print *, "NAN encountered in nx_LTE"

return
end function


function Etrans_sp(T_x, dTdr, mfp, sigma_N, alpha, Nwimps, niso)
implicit none
! Calculates WIMP transported energy (erg/g/s) using eq. (2.40) in https://arxiv.org/pdf/0809.1871.pdf
! The Spergel Press formalism doesn't actually us dT/dr. 
! I have it here in case we want to use the fgoth nx interpolation with nx_func (which requires dTdr)

integer, intent(in) :: niso
double precision, intent(in) :: T_x, Nwimps
double precision, intent(in) :: dTdr(nlines), mfp(nlines)
double precision, intent(in) :: alpha(niso), sigma_N(niso)
double precision :: n_0, mxg, alphaofR(nlines)
double precision :: R(nlines), phi(nlines), n_nuc(niso,nlines)
double precision :: n_x(nlines), species_indep(nlines), species_dep(nlines), sigma_nuc(niso)
double precision :: Etrans_sp(nlines)
integer :: i, j
! T_x in K, dTdr in K/r, mfp in cm, sigma_N in cm^2, 

R = tab_r*Rsun ! R in cm
phi = -tab_vesc**2/2.d0 ! phi in erg/g
mxg = mdm*1.782662d-24 ! WIMP mass in g
! n_nuc in cm^-3
do i=1,niso
!	do j=1,nlines
		n_nuc(i,:) = tab_mfr(:,i)*tab_starrho/AtomicNumber(i)/mnucg ! tab_starrho in gcm^-3
!	enddo
enddo

sigma_nuc = 2.d0*sigma_N ! Total WIMP-nucleus cross section in cm^2v. Only works for q/v independent cross-sections

! n_x in cm^-3. nx_func interpolates between nx_LTE and nx_iso using fgoth
n_x = nx_func(T_x, dTdr, mfp, sigma_N, alpha, Nwimps, niso)

! Separate calc into species dependent and independent factors
species_indep = 8.0d0*sqrt(2.d0/pi)*kB**(3.d0/2.d0)*n_x*(T_x-tab_T)/tab_starrho ! The species independent part

! Now sum over species to get the species dependent factor
species_dep=0.d0
!species_dep = sigma_nuc(1)*n_nuc(1,:)*mxg*mnucg/((mxg+mnucg)**2)*sqrt(tab_T/mnucg + T_x/mxg)
! Uncomment the next four lines for spin-independent scattering 
do i=1,niso
	species_dep = species_dep + sigma_nuc(i)*n_nuc(i,:)*mxg*mnucg*AtomicNumber(i)/((mxg+mnucg*AtomicNumber(i))**2)* &
		sqrt(tab_T/(mnucg*AtomicNumber(i)) + T_x/mxg)
enddo

Etrans_sp = species_indep*species_dep ! erg/g/s

! Useful when troubleshooting
open(55, file="/home/luke/summer_2021/mesa/test_files/Etrans_sp_params.txt")
write(55,*) "scalar params: T_x=", T_x, "m_x=", mxg, "m_nuc=", mnucg, "sigma_nuc=", sigma_nuc(1), &
	"nlines=", nlines, "niso=", niso
do i=1,nlines
	write(55,*) R(i), tab_T(i), n_x(i) !n_x(i), tab_starrho(i), n_nuc(1,i), species_indep(i), phi(i), Etrans_sp(i)
enddo
close(55)

return
end function


function Tx_integral(T_x, dTdr, mfp, sigma_N, alpha, Nwimps, niso)
implicit none
! Calculates the Tx defining integral 

integer, intent(in) :: niso
double precision, intent(in) :: T_x, Nwimps
double precision, intent(in) :: dTdr(nlines), mfp(nlines)
double precision, intent(in) :: alpha(niso), sigma_N(niso)
double precision :: alphaofR(nlines), R(nlines), integrand(nlines)
double precision :: Tx_integral
! integrand units: erg/cm/s
R = tab_r*Rsun

integrand = 4*pi*R**2*tab_starrho*Etrans_sp(T_x, dTdr, mfp, sigma_N, alpha, Nwimps, niso)

! integral is Etrans_tot (erg/s)
Tx_integral = trapz(R, integrand, nlines)

return
end function


function newtons_meth(f, dTdr, mfp, sigma_N, alpha, Nwimps, niso, guess_1, guess_2, reltolerance)
! Performs Newton's method to solve Tx_integral(T_x)=0 for T_x (the function returns T_x)
! The parameter f is the Tx_integral function
implicit none

integer, intent(in) :: niso
double precision :: f ! Tx_integral
double precision, intent(in) :: Nwimps, reltolerance, guess_1, guess_2
double precision, intent(in) :: dTdr(nlines), mfp(nlines)
double precision, intent(in) :: alpha(niso), sigma_N(niso)
double precision :: x_1, x_2, x_3, f1, f2, error
double precision :: newtons_meth
! m_x and m_p in grams, T_x, T_star in Kelvin, sigma in cm^2, n_nuc, n_x in cm^-3

! x_1 and x_2 are temperatures (K)
x_1 = guess_1
x_2 = guess_2
error = reltolerance + 1	! So that the first iteration is executed

! Newton's method loop
do while (error > reltolerance)
	! Update x_3 using Newton's method formula
	f1 = f(x_1, dTdr, mfp, sigma_N, alpha, Nwimps, niso)
	f2 = f(x_2, dTdr, mfp, sigma_N, alpha, Nwimps, niso)
	x_3 = x_2 - f2*(x_2-x_1)/(f2 - f1)
	error = abs(x_3-x_2)/x_2
	x_1 = x_2
	x_2 = x_3
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

end module spergelpressmod
