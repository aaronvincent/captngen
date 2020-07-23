!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Nonlocal WIMP heat transport module !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Contains the functions used in the Spergel Press section of transgen.f90. These are:
! 	-Etrans_nl: calculates the WIMP transported energy (eps_x) given the WIMP temperature (Tx)
!	-Tx_integral: to be used in newtons_meth
!	-newtons_meth: solves Tx_integral=0 which defines Tx 

module nonlocalmod
use capmod
implicit none

contains


function Etrans_nl(T_x, r, T_star, phi, rho_star, m_x, n_nuc, m_nuc, sigma_nuc, Nwimps, nlines, niso)
implicit none
! Calculates WIMP Etrans using eq. (2.40) in https://arxiv.org/pdf/0809.1871.pdf
integer, intent(in) :: nlines, niso
double precision, intent(in) :: T_x, m_x, Nwimps
double precision, intent(in) :: T_star(nlines), phi(nlines), rho_star(nlines), r(nlines)
double precision :: n_nuc(niso, nlines)
double precision, intent(in) :: m_nuc(niso), sigma_nuc(niso)
double precision, parameter :: k=1.38064852d-16, pi=3.1415962d0 ! kB in cgs
double precision :: species_indep(nlines), species_dep(nlines), n_x(nlines), n_0 ! Internal variables
double precision :: Etrans_nl(nlines)
integer :: i
! m_x and m_p in grams, T_x, T_star in Kelvin, sigma in cm^2, n_nuc, n_x in cm^-3

! A better nxIso estimate than in Transgen
n_x = exp(-m_x*phi/k/T_x) 
n_0 = Nwimps/trapz(r, 4.d0*pi*r**2.d0*n_x, nlines) ! Normalize so that integral(nx) = Nwimps
n_x = n_0*n_x

!print *, "In Etrans_nl: n_x = ", n_x(10)

! Separate calc into species dependent and independent factors for convenience
species_indep = 8.0d0*sqrt(2.d0/pi)*k**(3.d0/2.d0)*n_x*(T_x-T_star)/rho_star ! The species independent part

! Now sum over species to get the species dependent factor
species_dep=0.d0
species_dep = sigma_nuc(1)*n_nuc(1,:)*m_x*m_nuc(1)/((m_x+m_nuc(1))**2)*sqrt(T_star/m_nuc(1) + T_x/m_x)
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


function Tx_integral(T_x, r, T_star, phi, rho_star, m_x, n_nuc, m_nuc, sigma, Nwimps, nlines, niso)
implicit none
! Calculates the Tx defining integral 
! T_x, T_star in K, r in cm, phi in erg/g, rho_star in g/cm^3, m_x and m_nuc in g, n_nuc in cm^-3, sigma in cm^2

double precision, intent(in) :: T_x	! The independent variable. All others are params (for Newton's method)
integer, intent(in) :: nlines, niso
double precision, intent(in) :: m_x, Nwimps
double precision, intent(in) :: r(nlines), T_star(nlines), phi(nlines), rho_star(nlines)
double precision, intent(in) :: n_nuc(niso, nlines), m_nuc(niso), sigma(niso)
double precision :: integrand(nlines), n_x(nlines)
double precision :: Tx_integral
double precision :: k=1.38064852d-16, pi=3.14159265 ! Boltzmann constant in cgs

! integrand units: erg/cm/s
integrand = 4*pi*r**2*rho_star*Etrans_nl(T_x, r, T_star, phi, rho_star, m_x, n_nuc, m_nuc, sigma, & 
	Nwimps, nlines, niso)

! integral is Etrans_tot (erg/s)
Tx_integral = trapz(r, integrand, nlines)

return
end function


function newtons_meth(f, r, T_star, phi, rho_star, m_x, n_nuc, m_nuc, sigma_nuc, Nwimps, nlines, niso, &
	guess_1, guess_2, tolerance)
! Performs Newton's method to solve Tx_integral(T_x)=0 for T_x (the function returns T_x)
! The parameter f is the Tx_integral function
implicit none

double precision :: f ! Tx_integral
double precision, intent(in) :: tolerance, guess_1, guess_2
integer, intent(in) :: nlines, niso
integer :: i, half
double precision, intent(in) :: m_x, Nwimps
double precision, intent(in) :: r(nlines), T_star(nlines), phi(nlines), rho_star(nlines)
double precision, intent(in) :: n_nuc(niso, nlines)
double precision, intent(in) :: m_nuc(niso), sigma_nuc(niso)
double precision :: x_1, x_2, x_3, error, n_x(nlines)
double precision :: newtons_meth

! x_1 and x_2 are temperatures (K)
x_1 = guess_1
x_2 = guess_2
error = tolerance + 1	! So that the first iteration is executed

! Newton's method loop
do while (error > tolerance)
	! Update x_3 using Newton's method formula
	x_3 = x_2 - f(x_2,r,T_star,phi,rho_star,m_x,n_nuc,m_nuc,sigma_nuc,Nwimps,nlines,niso)*(x_2-x_1) &
		/(f(x_2,r,T_star,phi,rho_star,m_x,n_nuc,m_nuc,sigma_nuc,Nwimps,nlines,niso) - &
		f(x_1,r,T_star,phi,rho_star,m_x,n_nuc,m_nuc,sigma_nuc,Nwimps,nlines,niso))
	error = abs(x_3-x_2)
	x_1 = x_2
	x_2 = x_3
enddo

newtons_meth = x_3 ! The solution to the nonlinear equation

return
end function

subroutine fourier_smooth(x, y, x_even, y_even, cutoff, nlines, lensav, ierr)
integer, intent(in) :: nlines, lensav
integer :: ierr, i
double precision, intent(in) :: x(nlines), x_even(nlines), cutoff
double precision, intent(inout) :: y(nlines)
double precision :: y_even(nlines), work(nlines), wsave(lensav), bcoeff(nlines), ccoeff(nlines), dcoeff(nlines)
double precision :: ispline

! Make evenly spaced y array
call spline(x, y, bcoeff, ccoeff, dcoeff, nlines)
do i=1,nlines
	y_even(i) = ispline(x_even(i), x, y, bcoeff, ccoeff, dcoeff, nlines)
enddo

! Compute FFT of y
call dfft1i (nlines, wsave, lensav, ierr)  !Initialize (required by fftpack)
if (ierr /= 0) print *, "FFT initializer 'dfft1i' failed with error ", ierr
print *, "FFT initialized"

call dfft1f(nlines, 1, y_even, nlines, wsave, lensav, work, nlines, ierr) ! Take FFT
if (ierr /= 0) print *, "Forward FFT calculator 'dfft1f' failed with error ", ierr
! dTdr_even is now the array Fourier components of dTdr_even (the way fftpack works)
print *, "FFT Taken"

! Cut out top 95% of Fourier components
do i=1,nlines
	if (i > int(cutoff*nlines)) then
		y_even(i) = 0.d0
	endif
enddo

! Rebuild y with high frequency components cut out
call dfft1b(nlines, 1, y_even, nlines, wsave, lensav, work, nlines, ierr)
if (ierr /= 0) print *, "Backward FFT calculator 'dfft1b' failed with error ", ierr
print *, "Inverse FFT Taken"

! Evaluate y on original grid (ie go convert y_even --> y)
call spline(x_even, y_even, bcoeff, ccoeff, dcoeff, nlines)
do i=1,nlines
	y(i) = ispline(x(i), x_even, y_even, bcoeff, ccoeff, dcoeff, nlines)
enddo

end subroutine

end module nonlocalmod
