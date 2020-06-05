! Test program to calculate the WIMP temperature according to Spergel and Press

module nonlocal
implicit none

contains


function Etrans(T_x, T_star, phi, rho_star, m_x, n_nuc, m_nuc, sigma_nuc, nlines, niso)
implicit none
! Calculates Etrans using eq. (2.40) in https://arxiv.org/pdf/0809.1871.pdf
integer, intent(in) :: nlines, niso
double precision, intent(in) :: T_x, m_x
double precision, intent(in) :: T_star(nlines), phi(nlines), rho_star(nlines)
double precision, intent(in) :: n_nuc(niso, nlines)
double precision, intent(in) :: m_nuc(niso), sigma_nuc(niso)
double precision, parameter :: k=1.38064852d-23, pi=3.1415962
double precision :: species_indep(nlines), species_dep(nlines), n_x(nlines), n_0 ! Internal variables
double precision :: Etrans(nlines)
integer :: i

n_0 = 0.4*1.78266192d-27*1.0d6/m_x ! DM number density in m^-3
n_x = n_0*exp(m_x*phi(nlines)/k/T_x)*exp(-m_x*phi/k/T_x) ! At the edge of the star, set n_x=local galactic DM density

species_indep = 8*sqrt(2/pi)*k**(3/2)/rho_star*n_x*(T_star-T_x) ! The species independent part
! Now sum over species to get the species dependent factor
species_dep=0
do i=1,niso
species_dep = species_dep + sigma_nuc(i)*n_nuc(i,:)*m_x*m_nuc(i)/(m_x+m_nuc(i))**2*sqrt(T_star/m_nuc(i)+T_x/m_x)
enddo
!print *, "species_dep = ", species_dep(1)
Etrans = species_indep*species_dep
!print *, "species_indep = ", species_indep(1)
return
end function


function Tx_integral(T_x, r, T_star, phi, rho_star, m_x, n_nuc, m_nuc, sigma_nuc, nlines, niso)
implicit none
! Calculates the Tx defining integral, equation (4.10) in Spergel&Press "Effect of hypothetical, WIMPs on E transport..."

double precision, intent(in) :: T_x	! The independent variable. All others are params (for Newton's method)
integer, intent(in) :: nlines, niso
double precision, intent(in) :: m_x
double precision, intent(in) :: r(nlines), T_star(nlines), phi(nlines), rho_star(nlines)
double precision, intent(in) :: n_nuc(niso, nlines), m_nuc(niso), sigma_nuc(niso)
double precision :: integrand(nlines)
double precision :: Tx_integral
double precision :: k=1.38064852d-23 ! Boltzmann constant

!print *, "integral params:", T_x, r(nlines), T_star(nlines), phi(nlines), rho_star(nlines), m_x, &
!	n_nuc(1,nlines), m_nuc(1), sigma_nuc(2), nlines, niso
!print *, "phi = ", phi(int(nlines)/2)
integrand = r**2*Etrans(T_x, T_star, phi, rho_star, m_x, n_nuc, m_nuc, sigma_nuc, nlines, niso)
!print *, "integrand = ", integrand
Tx_integral = trapz(r, integrand, nlines)

return
end function


function newtons_meth(f, r, T_star, phi, rho_star, m_x, n_nuc, m_nuc, sigma_nuc, nlines, niso, &
	guess_1, guess_2, tolerance)
!Performs gradient descent to solve sp_int=0.
implicit none

double precision :: f ! sp_integral
double precision, intent(in) :: tolerance, guess_1, guess_2
integer, intent(in) :: nlines, niso
double precision, intent(in) :: m_x
double precision, intent(in) :: r(nlines), T_star(nlines), phi(nlines), rho_star(nlines)
double precision, intent(in) :: n_nuc(niso, nlines)
double precision, intent(in) :: m_nuc(niso), sigma_nuc(niso)
double precision :: x_1, x_2, x_3, error
double precision :: newtons_meth

x_1 = guess_1
x_2 = guess_2
error = tolerance + 1

do while (error > tolerance)
	! Update x_3 using Newton's method formula
	x_3 = x_2 - f(x_2,r,T_star,phi,rho_star,m_x,n_nuc,m_nuc,sigma_nuc,nlines,niso)*(x_2-x_1) &
	/(f(x_2,r,T_star,phi,rho_star,m_x,n_nuc,m_nuc,sigma_nuc,nlines,niso) - &
	f(x_1,r,T_star,phi,rho_star,m_x,n_nuc,m_nuc,sigma_nuc,nlines,niso))
	error = abs(x_3-x_2)
	x_1 = x_2
	x_2 = x_3
enddo

newtons_meth = x_3	! The solution to the nonlinear equation
return
end function

end module nonlocal



