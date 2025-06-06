!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Spergel-Press WIMP heat transport module !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Contains the functions used in the Spergel Press section of transgen.f90. These are:
!	-nx_isothermal: Calculates the WIMP density in the Spergel-Press scheme
! 	-Etrans_sp: calculates the WIMP transported energy (eps_x) given the WIMP temperature (Tx)
!	-Tx_integral: to be used in newtons_meth
!	-newtons_meth: solves Tx_integral=0 which defines Tx

! All units are cgs except tab_r and tab_dr
! I apologize for the long function calls.

module spergelpressmod
use sharedmod
use capmod
implicit none

contains


function nx_isothermal(T_x, Nwimps)
use phys, only : kB, pi
implicit none
double precision, intent(in) :: T_x, Nwimps
double precision :: nx_isothermal(nlines)
double precision :: n_0, mxg
double precision :: R(nlines), phi(nlines)
integer :: i
! Calculates the isothermal wimp number density using eq. (2.25) in https://arxiv.org/pdf/0809.1871.pdf


r = tab_r*Rsun ! cm
phi = -tab_vesc**2/2.d0 ! erg/g
mxg = mdm*1.782662d-24  ! g


!print*, 'nx_iso here'
! WIMP number density in isothermal approximation

!nx_isothermal = exp(-mxg*phi/kB/T_x)          !previous calulation that doesn't work above 8GeV
nx_isothermal = exp(-mxg*(phi-phi(1))/kB/T_x)  !the minus phi(1) lets the code run with a mass above 8 GeV

n_0 = Nwimps/trapz(r, 4.d0*pi*r**2.d0*nx_isothermal, nlines) ! Normalize so that integral(nx) = Nwimps
nx_isothermal = n_0*nx_isothermal

if (any(isnan(nx_isothermal))) print *, "NAN encountered in nx_isothermal"

return
end function


function Etrans_sp(T_x, sigma_N, Nwimps, niso)
use phys, only : mnuc, gev_erg, c0, pi, kB
implicit none
! Calculates WIMP transported energy (erg/g/s) using eq. (2.40) in https://arxiv.org/pdf/0809.1871.pdf

integer, intent(in) :: niso
double precision, intent(in) :: T_x, Nwimps
double precision, intent(in) :: sigma_N(niso)
double precision :: n_0, mxg, B, A, initial_q
double precision :: R(nlines), phi(nlines), n_nuc(niso,nlines)
double precision :: n_x(nlines), species_indep(nlines), species_dep(nlines), sigma_nuc(niso)
double precision :: Etrans_sp(nlines)
double precision, parameter :: mnucg = mnuc / gev_erg / c0**2
integer :: i, j, p
! T_x in K, sigma_N in cm^2,



R = tab_r*Rsun ! R in cm
phi = -tab_vesc**2/2.d0 ! phi in erg/g
mxg = mdm*1.782662d-24 ! WIMP mass in g
initial_q = q0*5.344d-14 !cgs conversion for q0

! n_nuc in cm^-3
do i=1,niso
	n_nuc(i,:) = tab_mfr(:,i)*tab_starrho/AtomicNumber(i)/mnucg ! tab_starrho in gcm^-3
enddo

sigma_nuc = 2.d0*sigma_N ! Total WIMP-nucleus cross section in cm^2v. Only works for q/v independent cross-sections

!print*, sigma_N, sigma_nuc

!print*,'Etrans here'
! isothermal WIMP number density in cm^-3.
n_x = nx_isothermal(T_x, Nwimps)

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
	species_indep = A*sqrt(2.d0/pi)*kB**(3.d0/2.d0)*n_x*(T_x-tab_T)/tab_starrho ! The species independent part
	do i=1,niso
	species_dep = species_dep + sigma_nuc(i)*n_nuc(i,:)*mxg*mnucg*AtomicNumber(i)/((mxg+mnucg*AtomicNumber(i))**2)* &
		(tab_T/(mnucg*AtomicNumber(i)) + T_x/mxg)**(1.d0/2.d0)
	enddo
	Etrans_sp = species_indep*species_dep ! erg/g/s
else if (nv .ne. 0) then
	! Separate calc into species dependent and independent factors
	species_indep = A*sqrt(2.d0/pi)*kB**(3.d0/2.d0+nv)*n_x*(T_x-tab_T)/tab_starrho/v0**(2.d0*nv) ! The species independent part
	do i=1,niso
	species_dep = species_dep + sigma_nuc(i)*n_nuc(i,:)*mxg*mnucg*AtomicNumber(i)/((mxg+mnucg*AtomicNumber(i))**2)* &
		(tab_T/(mnucg*AtomicNumber(i)) + T_x/mxg)**(1.d0/2.d0+nv)
	enddo
	Etrans_sp = species_indep*species_dep
else if (nq .ne. 0) then
	! Separate calc into species dependent and independent factors
	species_indep = A*sqrt(2.d0/pi)*kB**(3.d0/2.d0+nq)*n_x*(T_x-tab_T)/tab_starrho*B/(initial_q)**(2.d0*nq)* &
		(2.**nq)*mxg**(2.d0*nq) ! The species independent part
	do i=1,niso
	species_dep = species_dep + sigma_N(i)*n_nuc(i,:)*mxg*mnucg*AtomicNumber(i)/((mxg+mnucg*AtomicNumber(i))**2)* &
		(tab_T/(mnucg*AtomicNumber(i)) + T_x/mxg)**(1.d0/2.d0+nq)/(1.+mxg/(mnucg*AtomicNumber(i)))**(2.d0*nq)
	enddo
	Etrans_sp = species_indep*species_dep
end if


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

subroutine transport_sp_generic(n, temp_dm, num_dm, m_target, ndensity_target, integral_result)
	!! This gives the result of the Spergel & Press energy transfer with a given target isotope as defined in Eq. 2.10 of
	!! [[arXiv:2111.06895](https://arxiv.org/pdf/2111.06895#equation.2.10)], except for the interaction-dependent terms
	!! \( (1-Q) \sigma_\text{tot} \).
	use phys, only : pi, kB, gev_erg, c0
	implicit none
	integer, intent(in) :: n !! Total number of relative velocity \( z^{2n} \) terms in the integrand [\( 1 \)]
	double precision, intent(in) :: temp_dm !! Isothermal temperature of the dark matter [\( \text{K} \)]
	double precision, intent(in) :: num_dm !! Total number of dark matter particles in the star [\( 1 \)]
	double precision, intent(in) :: m_target !! Mass of the target isotope [\( \text{GeV} \)]
	double precision, intent(in) :: ndensity_target(:) !! Radial profile of the number density of the target isotope [\( \text{cm}^{-3} \)]
	double precision, intent(out) :: integral_result(:) !! [\( (\text{erg} \cdot \text{g}^{-1} \text{s}^{-1}) (\text{cm}^{-2}) {(\text{cm} \cdot \text{s}^{-1})}^{2n} \)]
	double precision :: a_factor

	a_factor = 2.d0**(2+n) * gamma(real(n)+3)
	!* @note
	! These `a_factor`s have been calculated via Sympy and Mathematica analytic solutions to the Eq. 2.10 linked above. They are
	! *very* slow to calculate, and I was unsucessful in convincing either CAS to produce a general expression, so I'm stuck with
	! simply (and slowly) looping over successive powers of \( z^{2n} \). In doing this I discovered that sucessive numerical
	! factors followed the recursive relation \( A_n=(2n+4)A_{n-1} \). The expression for the nth term \( A_n=2^{2+n}\Gamma(n+3) \)
	! matches the values I found using Sympy and Mathematica up to and including \( A_{12} \). @endnote
	!!

	integral_result = a_factor/tab_starrho * sqrt(2/pi) * mdm*m_target/(mdm+m_target)**2 * nx_isothermal(temp_dm, num_dm) &
		* ndensity_target * (temp_dm - tab_t) * kB * sqrt(((tab_t/m_target + temp_dm/mdm) * kB*gev_erg*c0**2)**(1+2*n))

end subroutine transport_sp_generic

subroutine transport_sp_qv(q_pow, v_pow, sigma_0, temp_dm, num_dm, m_target, ndensity_target, epsilon_sp)
	!! \( \epsilon_\text{SP} \) of a given target isotope as defined in Eq. 2.10 of
	!! [[arXiv:2111.06895](https://arxiv.org/pdf/2111.06895#equation.2.10)]. Using a momentum-velocity scaled differential cross
	!! section defined as
	!! \( \frac{\mathrm{d} \sigma}{\mathrm{d} \cos\theta} = \sigma_0 {\frac{q}{q_0}}^{2n_q} {\frac{v}{v_0}}^{2n_v} \).
	use phys, only : c0
	implicit none
	integer, intent(in) :: q_pow !! The number of powers of transfer momentum \( q^{2 q_\text{pow}} \) [\( 1 \)]
	integer, intent(in) :: v_pow !! The number of powers of velocity \( v^{2 v_\text{pow}} \) [\( 1 \)]
	double precision, intent(in) :: sigma_0 !! Reference cross section [\( \text{cm}^2 \)]
	double precision, intent(in) :: temp_dm !! Isothermal temperature of the dark matter [\( \text{K} \)]
	double precision, intent(in) :: num_dm !! Total number of dark matter particles in the star [\( 1 \)]
	double precision, intent(in) :: m_target !! Mass of the target isotope [\( \text{GeV} \)]
	double precision, intent(in) :: ndensity_target(:) !! Radial profile of the number density of the target isotope [\( \text{cm}^{-3} \)]
	double precision, intent(out) :: epsilon_sp(:) !! [\( \text{erg} \cdot \text{g}^{-1} \text{s}^{-1} \)]
	double precision :: sigma_tot
	double precision, allocatable :: integral_result(:)

	if (.not. allocated(integral_result)) then
		allocate(integral_result(size(epsilon_sp)))
	end if

	sigma_tot = sigma_0 * 2/(q_pow+1) * (2*mdm/(c0*(1+mu)*q0))**(2*q_pow) * v0**(-2*v_pow)
	call transport_sp_generic(q_pow+v_pow, temp_dm, num_dm, m_target, ndensity_target, integral_result)

	epsilon_sp = ( 1 - (-q_pow/(q_pow+2))) * sigma_tot * integral_result

end subroutine transport_sp_qv

subroutine transport_sp_nreo(q_pow, w_pow, prefactor, temp_dm, num_dm, m_target, ndensity_target, epsilon_sp)
	!! \( \epsilon_\text{SP} \) of a given target isotope as defined in Eq. 2.10 of
	!! [[arXiv:2111.06895](https://arxiv.org/pdf/2111.06895#equation.2.10)]. Using an NREO differential cross section defined as
	!! \begin{align}
	!! \frac{\mathrm{d} \sigma_T}{\mathrm{d} \cos\theta} &= \frac{\mathrm{d} E_R}{\mathrm{d} \cos\theta} \frac{\mathrm{d} \sigma_T}{\mathrm{d} E_R} \, , \\
	!! \frac{\mathrm{d} \sigma_T}{\mathrm{d} E_R} &= \frac{- P_{T,n_q,n_w} \hbar^2 c^{2(1-n_q)}}{(2J + 1)(1 + n_q)} {\left( \frac{2m_\chi}{1 + \mu} \right)}^{2(1+n_q)} w^{2(n_q+n_w)} \, .
	!! \end{align}
	!! The units of the prefactor are: \([P_{T,n_q,n_w}] = \text{GeV}^{-4} \cdot \text{GeV}^{-2n_q} \cdot (\text{cm} \cdot \text{s}^{-1})^{-2n_w} \).
	use phys, only : hbar, c0
	implicit none
	integer, intent(in) :: q_pow !! The number of powers of transfer momentum \( q^{2 q_\text{pow}} \) [\( 1 \)]
	integer, intent(in) :: w_pow !! The number of powers of velocity \( w^{2 w_\text{pow}} \) [\( 1 \)]
	double precision, intent(in) :: prefactor !! Numerical RW prefactor for the given `q_pow` and `w_pow`, divided by \( (2J+1) \) [\( \text{GeV}^{-4-2n_q} \cdot (\text{cm} \cdot \text{s}^{-1})^{-2n_w} \)]
	double precision, intent(in) :: temp_dm !! Isothermal temperature of the dark matter [\( \text{K} \)]
	double precision, intent(in) :: num_dm !! Total number of dark matter particles in the star [\( 1 \)]
	double precision, intent(in) :: m_target !! Mass of the target isotope [\( \text{GeV} \)]
	double precision, intent(in) :: ndensity_target(:) !! Radial profile of the number density of the target isotope [\( \text{cm}^{-3} \)]
	double precision, intent(out) :: epsilon_sp(:) !! [\( \text{erg} \cdot \text{g}^{-1} \text{s}^{-1} \)]
	double precision :: sigma_tot
	double precision, allocatable :: integral_result(:)

	if (.not. allocated(integral_result)) then
		allocate(integral_result(size(epsilon_sp)))
	end if

	sigma_tot = abs( -prefactor * (hbar * c0**(1-q_pow) * (2*mdm/(1+mu))**(1+q_pow))**2 / (1+q_pow) )
	!* @warning
	! This is *not* complete, I'm still concerned about how we calculate \( \sigma_\text{tot} \) in the NREO formalism. @endwarning
	!!
	call transport_sp_generic(q_pow+w_pow, temp_dm, num_dm, m_target, ndensity_target, integral_result)

	epsilon_sp = ( 1 - (-q_pow/(q_pow+2))) * sigma_tot * integral_result

end subroutine transport_sp_nreo

	! This is *not* complete, I'm still concerned about how we calculate \( \sigma_\text{tot} \) in the NREO formalism. For now we
	! enforce that the total cross section is strictly positive, but is there a more convincing argument beyond that? @endwarning
function Tx_integral(T_x, sigma_N, Nwimps, niso)
use phys, only : pi
implicit none
! Calculates the Tx defining integral

integer, intent(in) :: niso
double precision, intent(in) :: T_x, Nwimps
double precision, intent(in) :: sigma_N(niso)
double precision :: R(nlines), integrand(nlines)
double precision :: Tx_integral

! integrand units: erg/cm/s
R = tab_r*Rsun

!print*, 'TX here'
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
