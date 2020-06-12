! Test program to calculate the WIMP temperature according to Spergel and Press

module nonlocalmod
use capmod, only:trapz
implicit none

contains


function Etrans_nl(T_x, r, T_star, phi, rho_star, m_x, n_nuc, m_nuc, sigma_nuc, Nwimps, nlines, niso)
implicit none
! Calculates Etrans using eq. (2.40) in https://arxiv.org/pdf/0809.1871.pdf
integer, intent(in) :: nlines, niso
double precision, intent(in) :: T_x, m_x, Nwimps
double precision, intent(in) :: T_star(nlines), phi(nlines), rho_star(nlines), r(nlines)
double precision, intent(in) :: n_nuc(niso, nlines)
double precision, intent(in) :: m_nuc(niso), sigma_nuc(niso)
double precision, parameter :: k=1.38064852d-16, pi=3.1415962 ! kB in cgs
double precision :: species_indep(nlines), species_dep(nlines), n_x(nlines), n_0 ! Internal variables
double precision :: Etrans_nl(nlines)
integer :: i, half
!print *, "In Etrans"
half = int(nlines/2)

n_x = exp(-m_x*phi/k/T_x) 
n_0 = Nwimps/trapz(r, 4*pi*r**2*n_x, nlines) ! Normalize so that integral(nx) = Nwimps
n_x = n_0*n_x

species_indep = 8*sqrt(2/pi)*k**(3/2)/rho_star*n_x*(T_x-T_star) ! The species independent part
! Now sum over species to get the species dependent factor
species_dep=0
!print *, "Entering species loop"
do i=1,niso
	species_dep = species_dep + sigma_nuc(i)*n_nuc(i,:)*m_x*m_nuc(i)/(m_x+m_nuc(i))**2 * & 
		sqrt(T_star/m_nuc(i) + T_x/m_x)
enddo
!print *, "Species loop finished"
Etrans_nl = species_indep*species_dep ! erg/g/s
if (isnan(T_x) .eqv. .false.) then
	open(1, file="/home/luke/summer_2020/mesa/captngen/Etrans.txt")
	write(1,*) "params: T_x=", T_x, "phi=", phi(2), "m_x=", m_x, "kB=", k, &
		"m_nuc=", m_nuc(1), "sigma_nuc=", sigma_nuc(1)
		
	write(1,*) "--------Etrans---------"
	do i=1, nlines
		write(1,*) i, Etrans_nl(i)
	enddo
	
	write(1,*) "--------species_indep---------"
	do i=1, nlines
		write(1,*) i, species_indep(i)
	enddo
	
	write(1,*) "--------species_dep---------"
	do i=1, nlines
		write(1,*) i, species_dep(i)
	enddo
	
	write(1,*) "--------m_x*(phi(nlines)-phi)/(kB*T_x)---------"
	do i=1, nlines
		write(1,*) i, m_x*(phi(nlines)-phi(i))/(k*T_x)
	enddo
	
	write(1,*) "--------n_nuc(1,:)---------"
	do i=1, nlines
		write(1,*) i, n_nuc(1,i)
	enddo
	
	write(1,*) "--------rho_star---------"
	do i=1, nlines
		write(1,*) i, rho_star(i)
	enddo
	close(1)
endif
return
end function


function Tx_integral(T_x, r, T_star, phi, rho_star, m_x, n_nuc, m_nuc, sigma_nuc, Nwimps, nlines, niso)
implicit none
! Calculates the Tx defining integral, equation (4.10) in Spergel&Press "Effect of hypothetical, WIMPs on E transport..."

double precision, intent(in) :: T_x	! The independent variable. All others are params (for Newton's method)
integer, intent(in) :: nlines, niso
double precision, intent(in) :: m_x, Nwimps
double precision, intent(in) :: r(nlines), T_star(nlines), phi(nlines), rho_star(nlines)
double precision, intent(in) :: n_nuc(niso, nlines), m_nuc(niso), sigma_nuc(niso)
double precision :: integrand(nlines)
double precision :: Tx_integral
double precision :: k=1.38064852d-16, pi=3.14159265 ! Boltzmann constant in cgs

!print *, "In Tx_integral"

!print *, "integral params:", T_x, r(nlines), T_star(nlines), phi(nlines), rho_star(nlines), m_x, &
!	n_nuc(1,nlines), m_nuc(1), sigma_nuc(2), nlines, niso
!print *, "phi = ", phi(int(nlines)/2)
! integrand units: erg/cm/s
integrand = 4*pi*r**2*rho_star*Etrans_nl(T_x, r, T_star, phi, rho_star, m_x, n_nuc, m_nuc, sigma_nuc, & 
	Nwimps, nlines, niso)
!print *, "integrand = ", integrand
Tx_integral = trapz(r, integrand, nlines)

return
end function


function newtons_meth(f, r, T_star, phi, rho_star, m_x, n_nuc, m_nuc, sigma_nuc, Nwimps, nlines, niso, &
	guess_1, guess_2, tolerance)
!Performs gradient descent to solve sp_int=0.
implicit none

double precision :: f ! sp_integral
double precision, intent(in) :: tolerance, guess_1, guess_2
integer, intent(in) :: nlines, niso
integer :: i, half
double precision, intent(in) :: m_x, Nwimps
double precision, intent(in) :: r(nlines), T_star(nlines), phi(nlines), rho_star(nlines)
double precision, intent(in) :: n_nuc(niso, nlines)
double precision, intent(in) :: m_nuc(niso), sigma_nuc(niso)
double precision :: x_1, x_2, x_3, error
double precision :: newtons_meth

x_1 = guess_1
x_2 = guess_2
error = tolerance + 1
!print *, "Entering Newton's loop"
i=0
half = int(nlines/2)
do while (error > tolerance)
	! Update x_3 using Newton's method formula
!	print *, "Calculating x_3 ..."
	x_3 = x_2 - f(x_2,r,T_star,phi,rho_star,m_x,n_nuc,m_nuc,sigma_nuc,Nwimps,nlines,niso)*(x_2-x_1) &
	/(f(x_2,r,T_star,phi,rho_star,m_x,n_nuc,m_nuc,sigma_nuc,Nwimps,nlines,niso) - &
	f(x_1,r,T_star,phi,rho_star,m_x,n_nuc,m_nuc,sigma_nuc,Nwimps,nlines,niso))
	error = abs(x_3-x_2)
!	print *, "error = ", error
	x_1 = x_2
!	print *, "x_1 = ", x_1
	x_2 = x_3
!	print *, "x_2 = ", x_2
	i = i + 1
!	print *, "i = ", i
	newtons_meth = f(x_2,r,T_star,phi,rho_star,m_x,n_nuc,m_nuc,sigma_nuc,Nwimps,nlines,niso)
!	print *, "f(x_2) finished"
	newtons_meth = f(x_1,r,T_star,phi,rho_star,m_x,n_nuc,m_nuc,sigma_nuc,Nwimps,nlines,niso)
!	print *, "f(x_1) finished"
!	print *, "x_3=", x_3
!	print *, "i=", i
enddo
!print *, "Newton's loop finished"
newtons_meth = x_3	! The solution to the nonlinear equation
return
end function

end module nonlocalmod



