!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Spergel-Press WIMP heat transport module !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Contains the functions used in the Spergel Press section of transgen.f90. These are:
!	-isothermal_dm_num_density: Calculates the WIMP density in the Spergel-Press scheme
! 	-transport_energy_sp: calculates the WIMP transported energy (eps_x) given the WIMP temperature (Tx)
!	-luminosity_dm: to be used in newtons_method
!	-newtons_method: solves luminosity_dm=0 which defines the isothermal dark matter temperature

! All units are cgs except star_r and star_dr
! I apologize for the long function calls.

module spergel_press_mod
use capture_mod
implicit none

double precision, parameter :: kB=1.38064852d-16, mnucg=1.6726219e-24

contains


function isothermal_dm_num_density(iso_temperature_dm, num_wimps)
implicit none
double precision, intent(in) :: iso_temperature_dm, num_wimps
double precision :: isothermal_dm_num_density(nlines)
double precision :: n_0, mxg
double precision :: R(nlines), phi(nlines)
integer :: i
! Calculates the isothermal wimp number density using eq. (2.25) in https://arxiv.org/pdf/0809.1871.pdf


r = star_r*radius_star ! cm
phi = -star_escape**2/2.d0 ! erg/g
mxg = m_dm*1.782662d-24  ! g


!print*, 'nx_iso here'
! WIMP number density in isothermal approximation

!isothermal_dm_num_density = exp(-mxg*phi/kB/iso_temperature_dm)          !previous calculation that doesn't work above 8GeV
isothermal_dm_num_density = exp(-mxg*(phi-phi(1))/kB/iso_temperature_dm)  !the minus phi(1) lets the code run with a mass above 8 GeV

n_0 = num_wimps/trapezoid(r, 4.d0*pi*r**2.d0*isothermal_dm_num_density, nlines) ! Normalize so that integral(nx) = num_wimps
isothermal_dm_num_density = n_0*isothermal_dm_num_density

if (any(isnan(isothermal_dm_num_density))) print *, "NAN encountered in isothermal_dm_num_density"

return
end function


function transport_energy_sp(iso_temperature_dm, diff_sigma, num_wimps, num_isotopes)
implicit none
! Calculates WIMP transported energy (erg/g/s) using eq. (2.40) in https://arxiv.org/pdf/0809.1871.pdf

integer, intent(in) :: num_isotopes
double precision, intent(in) :: iso_temperature_dm, num_wimps
double precision, intent(in) :: diff_sigma(num_isotopes)
double precision :: n_0, mxg, B, A, initial_q
double precision :: R(nlines), phi(nlines), n_nuc(num_isotopes,nlines)
double precision :: n_x(nlines), species_indep(nlines), species_dep(nlines), sigma_nuc(num_isotopes)
double precision :: transport_energy_sp(nlines)
integer :: i, j, p
! iso_temperature_dm in K, diff_sigma in cm^2,



R = star_r*radius_star ! R in cm
phi = -star_escape**2/2.d0 ! phi in erg/g
mxg = m_dm*1.782662d-24 ! WIMP mass in g
initial_q = q0*5.344d-14 !cgs conversion for q0

! n_nuc in cm^-3
do i=1,num_isotopes
	n_nuc(i,:) = star_fractions(:,i)*star_rho/atomic_nums(i)/mnucg ! star_rho in gcm^-3
end do

sigma_nuc = 2.d0*diff_sigma ! Total WIMP-nucleus cross section in cm^2v. Only works for q/v independent cross-sections

!print*, diff_sigma, sigma_nuc

!print*,'Etrans here'
! isothermal WIMP number density in cm^-3.
n_x = isothermal_dm_num_density(iso_temperature_dm, num_wimps)

p = (nv + nq)
if ((p .eq. 0)) then
	A = 8.d0
	B = 0.d0
else if ((p .eq. 1)) then
	A = 48.d0
	B = 8.d0/3.d0
else if ((p .eq. 2)) then
	A = 384.d0
	B = 4.d0
else if ((p .eq. -1)) then
	A = 2.d0
	B = 2.d0
end if

species_dep=0.d0

if ( (nq .eq. 0) .and. (nv .eq. 0) ) then
	! Separate calc into species dependent and independent factors
	species_indep = A*sqrt(2.d0/pi)*kB**(3.d0/2.d0)*n_x*(iso_temperature_dm-star_temp)/star_rho ! The species independent part
	do i=1,num_isotopes
	species_dep = species_dep + sigma_nuc(i)*n_nuc(i,:)*mxg*mnucg*atomic_nums(i)/((mxg+mnucg*atomic_nums(i))**2)* &
		(star_temp/(mnucg*atomic_nums(i)) + iso_temperature_dm/mxg)**(1.d0/2.d0)
	end do
	transport_energy_sp = species_indep*species_dep ! erg/g/s
else if (nv .ne. 0) then
	! Separate calc into species dependent and independent factors
	species_indep = A*sqrt(2.d0/pi)*kB**(3.d0/2.d0+nv)*n_x*(iso_temperature_dm-star_temp)/star_rho/v0**(2.d0*nv) ! The species independent part
	do i=1,num_isotopes
	species_dep = species_dep + sigma_nuc(i)*n_nuc(i,:)*mxg*mnucg*atomic_nums(i)/((mxg+mnucg*atomic_nums(i))**2)* &
		(star_temp/(mnucg*atomic_nums(i)) + iso_temperature_dm/mxg)**(1.d0/2.d0+nv)
	end do
	transport_energy_sp = species_indep*species_dep
else if (nq .ne. 0) then
	! Separate calc into species dependent and independent factors
	species_indep = A*sqrt(2.d0/pi)*kB**(3.d0/2.d0+nq)*n_x*(iso_temperature_dm-star_temp)/star_rho*B/(initial_q)**(2.d0*nq)* &
		(2.**nq)*mxg**(2.d0*nq) ! The species independent part
	do i=1,num_isotopes
	species_dep = species_dep + diff_sigma(i)*n_nuc(i,:)*mxg*mnucg*atomic_nums(i)/((mxg+mnucg*atomic_nums(i))**2)* &
		(star_temp/(mnucg*atomic_nums(i)) + iso_temperature_dm/mxg)**(1.d0/2.d0+nq)/(1.+mxg/(mnucg*atomic_nums(i)))**(2.d0*nq)
	end do
	transport_energy_sp = species_indep*species_dep
end if


!! Useful when troubleshooting
!open(55, file="/home/luke/summer_2021/mesa/test_files/Etrans_sp_params.txt")
!write(55,*) "scalar params: iso_temperature_dm=", iso_temperature_dm, "m_x=", mxg, "m_nuc=", mnucg, "sigma_nuc=", sigma_nuc(1), &
!	"nlines=", nlines, "num_isotopes=", num_isotopes
!do i=1,nlines
!	write(55,*) R(i), star_temp(i), n_x(i), transport_energy_sp(i) !n_x(i), star_rho(i), n_nuc(1,i), species_indep(i), phi(i)
!end do
!close(55)

return
end function


function luminosity_dm(iso_temperature_dm, diff_sigma, num_wimps, num_isotopes)
implicit none
! Calculates the dark matter temperature defining integral

integer, intent(in) :: num_isotopes
double precision, intent(in) :: iso_temperature_dm, num_wimps
double precision, intent(in) :: diff_sigma(num_isotopes)
double precision :: R(nlines), integrand(nlines)
double precision :: luminosity_dm

! integrand units: erg/cm/s
R = star_r*radius_star

!print*, 'TX here'
integrand = 4*pi*R**2*star_rho*transport_energy_sp(iso_temperature_dm, diff_sigma, num_wimps, num_isotopes)

! integral is Etrans_tot (erg/s)
luminosity_dm = trapezoid(R, integrand, nlines)

return
end function


function newtons_method(func, diff_sigma, num_wimps, num_isotopes, guess_1, guess_2, relative_tolerance)
! Performs Newton's method to solve luminosity_dm(T_x)=0 for T_x (the function returns T_x)
! The parameter func is the luminosity_dm function
implicit none

integer, intent(in) :: num_isotopes
double precision :: func ! luminosity_dm
double precision, intent(in) :: num_wimps, relative_tolerance, guess_1, guess_2
double precision, intent(in) :: diff_sigma(num_isotopes)
double precision :: x_1, x_2, x_3, f1, f2, error
double precision :: newtons_method
! m_x and m_p in grams, T_x, T_star in Kelvin, sigma in cm^2, n_nuc, n_x in cm^-3

! x_1 and x_2 are temperatures (K)
x_1 = guess_1
x_2 = guess_2
error = relative_tolerance + 1	! So that the first iteration is executed

! Newton's method loop
do while (error > relative_tolerance)
	! Update x_3 using Newton's method formula
	f1 = func(x_1, diff_sigma, num_wimps, num_isotopes)
	f2 = func(x_2, diff_sigma, num_wimps, num_isotopes)
	x_3 = x_2 - f2*(x_2-x_1)/(f2 - f1)
	error = abs(x_3-x_2)/x_2
	x_1 = x_2
	x_2 = x_3
end do

newtons_method = x_3 ! The solution to the nonlinear equation

return
end function

function binary_search(func, diff_sigma, num_wimps, num_isotopes, guess_1, guess_2, relative_tolerance)
integer, intent(in) :: num_isotopes
integer :: i
double precision :: func ! luminosity_dm
double precision, intent(in) :: num_wimps, relative_tolerance, guess_1, guess_2
double precision, intent(in) :: diff_sigma(num_isotopes)
double precision :: x_1, x_2, x_3, f1, f2, f3, error
double precision :: binary_search


! x_1 and x_2 are temperatures (K)
x_1 = guess_1
x_2 = guess_2
error = relative_tolerance + 1.d0	! So that the first iteration is executed

! Binary search loop
i = 0

do while (error > relative_tolerance)
	x_3 = (x_1 + x_2)/2.d0
	f1 = func(x_1, diff_sigma, num_wimps, num_isotopes)
	f2 = func(x_2, diff_sigma, num_wimps, num_isotopes)
	f3 = func(x_3, diff_sigma, num_wimps, num_isotopes)
	if (f3 == 0.d0) then
		exit
	else if (f1*f3 .gt. 0) then ! if f1 and f3 have the same sign
		x_1 = x_3
	else if (f2*f3 .gt. 0) then
		x_2 = x_3
	endif
	error = abs(x_2-x_1)/x_2
	i = i + 1
end do

binary_search = x_3

return
end function

subroutine fourier_smooth(x, y, x_even, y_even, cutoff, noise_indicator, num_lines, prime_length, error)
! Cuts out the high frequency components of y. E.g. if cutoff=0.05, the top 95% of frequency components are cut
! Also returns a "noise indicator" - The sum of the frequency components above the cutoff
integer, intent(in) :: num_lines, prime_length
integer :: error, i
double precision, intent(in) :: x(num_lines), x_even(num_lines), cutoff
double precision, intent(inout) :: y(num_lines)
double precision, intent(out) :: noise_indicator
double precision :: y_even(num_lines), work(num_lines), wsave(prime_length), bcoeff(num_lines), ccoeff(num_lines), dcoeff(num_lines)
double precision :: ispline, denominator

! Make evenly spaced y array
call spline(x, y, bcoeff, ccoeff, dcoeff, num_lines)
do i=1,num_lines
	y_even(i) = ispline(x_even(i), x, y, bcoeff, ccoeff, dcoeff, num_lines)
end do

! Compute FFT of y
call dfft1i (num_lines, wsave, prime_length, error)  !Initialize (required by fftpack)
if (error /= 0) print *, "FFT initializer 'dfft1i' failed with error ", error

call dfft1f(num_lines, 1, y_even, num_lines, wsave, prime_length, work, num_lines, error) ! Take FFT
if (error /= 0) print *, "Forward FFT calculator 'dfft1f' failed with error ", error
! dTdr_even is now the array Fourier components of dTdr_even (the way fftpack works)

noise_indicator = 0.d0
! Take the ratio of high frequency components to low frequency components as a measure of how noisy the data is
do i=int(cutoff*num_lines),num_lines
	noise_indicator = noise_indicator + abs(y_even(i))
end do
do i=1,int(cutoff*num_lines)
	denominator = denominator + abs(y_even(i))
end do
noise_indicator = noise_indicator/denominator

! Cut out top 100*(1-cutoff)% of Fourier components
do i=1,num_lines
	if (i > int(cutoff*num_lines)) then
		y_even(i) = 0.d0
	endif
end do

! Rebuild y with high frequency components cut out
call dfft1b(num_lines, 1, y_even, num_lines, wsave, prime_length, work, num_lines, error)
if (error /= 0) print *, "Backward FFT calculator 'dfft1b' failed with error ", error

! Evaluate y on original grid (ie go convert y_even --> y)
call spline(x_even, y_even, bcoeff, ccoeff, dcoeff, num_lines)
do i=1,num_lines
	y(i) = ispline(x(i), x_even, y_even, bcoeff, ccoeff, dcoeff, num_lines)
end do

end subroutine

function rolling_average(y, num_lines)
! Takes a 1D array f of length N, returns an array of length N whose ith entry
! is the average of f(i) and its 4 nearest neighbours
integer, intent(in) ::  num_lines
double precision, intent(in) :: y(num_lines)
integer :: i, j
double precision :: rolling_average(num_lines)

do i=5,num_lines-4
    rolling_average(i) = 0.d0
    do j=-4,4
        rolling_average(i) = rolling_average(i) + y(i+j)
    end do
    rolling_average(i) = rolling_average(i)/9.d0
end do
! do boundary values manually
rolling_average(1) = (y(1)+y(2)+y(3))/3.d0
rolling_average(2) = (y(1)+y(2)+y(3)+y(4))/4.d0
rolling_average(3) = (y(1)+y(2)+y(3)+y(4)+y(5)+y(6))/6.d0
rolling_average(4) = (y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)+y(8))/8.d0
rolling_average(num_lines-1) = (y(num_lines-7)+y(num_lines-6)+y(num_lines-5)+y(num_lines-4) &
							+y(num_lines-3)+y(num_lines-2)+y(num_lines-1)+y(num_lines))/8.d0
rolling_average(num_lines-1) = (y(num_lines-5)+y(num_lines-4)+y(num_lines-3)+y(num_lines-2)+y(num_lines-1)+y(num_lines))/6.d0
rolling_average(num_lines-1) = (y(num_lines-3)+y(num_lines-2)+y(num_lines-1)+y(num_lines))/4.d0
rolling_average(num_lines) = (y(num_lines-2)+y(num_lines-1)+y(num_lines))/3.d0

return
end function

end module spergel_press_mod
