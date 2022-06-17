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
! T_x in K, sigma_N in cm^2,

R = tab_r*Rsun ! R in cm
phi = -tab_vesc**2/2.d0 ! phi in erg/g
mxg = mdm*1.782662d-24 ! WIMP mass in g
! n_nuc in cm^-3
do i=1,niso
	n_nuc(i,:) = tab_mfr(:,i)*tab_starrho/AtomicNumber(i)/mnucg ! tab_starrho in gcm^-3
enddo

sigma_nuc = 2.d0*sigma_N ! Total WIMP-nucleus cross section in cm^2v. Only works for q/v independent cross-sections

! isothermal WIMP number density in cm^-3.
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
! m_x and m_p in grams, T_x, T_star in Kelvin, sigma in cm^2, n_nuc, n_x in cm^-3

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
