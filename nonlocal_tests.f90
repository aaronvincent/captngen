! Test program to calculate the WIMP temperature according to Spergel and Press

module nonlocal
implicit none

contains


function Etrans_nl(T_x, T_star, phi, rho_star, m_x, n_nuc, m_nuc, sigma_nuc, nlines, niso)
implicit none
! Calculates Etrans using eq. (2.40) in https://arxiv.org/pdf/0809.1871.pdf
integer, intent(in) :: nlines, niso
double precision, intent(in) :: T_x, m_x
double precision, intent(in) :: T_star(nlines), phi(nlines), rho_star(nlines)
double precision, intent(in) :: n_nuc(niso, nlines)
double precision, intent(in) :: m_nuc(niso), sigma_nuc(niso)
double precision, parameter :: k=1.38064852d-23, pi=3.1415962
double precision :: species_indep(nlines), species_dep(nlines), n_x(nlines), n_0 ! Internal variables
double precision :: Etrans_nl(nlines)
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
Etrans_nl = species_indep*species_dep
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
integrand = r**2*Etrans_nl(T_x, T_star, phi, rho_star, m_x, n_nuc, m_nuc, sigma_nuc, nlines, niso)
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


!Fast trapezoidal integral (copied from gencap.f90)
      function trapz(x,y,flen)
      implicit none
      integer, intent(in) :: flen
      double precision, intent (in) :: x(flen), y(flen)
      double precision trapz
      integer i

      trapz = y(1)*(x(2)-x(1))/2. + y(flen)*(x(flen)-x(flen-1))/2.
      do i = 2,flen-1;
        trapz = trapz + y(i)*(x(i)-x(i-1))
      end do

      return
      end function


!function binary_search(f, nlines, mx, mp, np, Tp, phi, r, guess_1, guess_2, tolerance)
!double precision :: f ! sp_integral
!double precision, intent(in) :: tolerance, guess_1, guess_2
!integer, intent(in) :: nlines
!double precision, intent(in) :: mx, mp ! 2 scalar parameters
!! 4 vector parameters
!double precision, intent(in) :: np(nlines), Tp(nlines), phi(nlines), r(nlines)
!double precision :: x_1, x_2, x_mid, error, f_1, f_2, f_mid 
!double precision :: binary_search

!x_1 = guess_1
!x_2 = guess_2
!error = tolerance + 1

!do while (error > tolerance)
!x_mid = (x_1+x_2)/2
!f_1 = f(x_1,nlines,mx,mp,np,Tp,phi,r)
!f_2 = f(x_2,nlines,mx,mp,np,Tp,phi,r)
!f_mid = f(x_mid,nlines,mx,mp,np,Tp,phi,r)

!!print *, x_1, x_2, f_1, f_2
!if (f_mid*f_1 > 0) then ! If f1, fmid are the same sign, set x1 to xmid
!	x_1 = x_mid
!endif
!if (f_mid*f_2 > 0) then ! If f2, fmid are the same sign, set x2 to xmid
!	x_2 = x_mid
!endif
!error = abs(x_2-x_1) ! error=interval
!enddo

!binary_search = (x_1+x_2)/2

!end function


! The following reads in solar parameters. Copied from gencap, easier than using module
! capmod while I'm testing
! read in solar parameters from Aldo Serenelli-style files, with header removed
    
    subroutine get_nlines(filename, nlines)
    character*300 :: filename
    integer :: nlines, iostatus=0
    
    ! count number of lines
    open(99,file=filename)
    nlines=0
    do
    read(99,*, iostat=iostatus)
    if(iostatus/=0) then ! to avoid end of file error.
    exit
    else
    nlines=nlines+1
    end if
    end do
    close(99)
    nlines = nlines -1
    
    end subroutine
    
    
    subroutine get_solar_params(filename, nlines, tab_r, tab_starrho, tab_T, tab_g, &
	tab_mencl, phi, tab_mfr)
    character*300 :: filename
    double precision :: Temp, Pres, Lumi !these aren't used, but dummies are required
    integer :: i,j, nlines,iostatus
    double precision, parameter :: GMoverR=1.908e11, Rsun = 69.57d9 
    double precision :: tab_r(nlines), tab_starrho(nlines), tab_T(nlines), tab_g(nlines), &
	tab_mencl(nlines), phi(nlines), tab_mfr(nlines, 29)

    open(99,file=filename)
    do i=1,nlines
    read(99,*) tab_mencl(i),tab_r(i), tab_T(i), tab_starrho(i), Pres, Lumi, tab_mfr(i,:)
    end do
    close(99)
	
    phi(nlines) = -GMoverR
    do i = 1,nlines-1
    j = nlines-i !trapezoid integral
    phi(j) = phi(j+1) + GMoverR*(tab_r(j)-tab_r(j+1))/2.*(tab_mencl(j)/tab_r(j)**2+tab_mencl(j+1)/tab_r(j+1)**2)
    end do
    tab_g(nlines) = tab_g(nlines-1)
    return
    end


end module nonlocal



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Testing program
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program tests
use nonlocal
implicit none

double precision :: Tx, mp, mx, integral, xy_int, guess_1, guess_2, tolerance, sigma_x, L_xtot, n_0, sigma_0
double precision, parameter :: kB = 1.38064852d-23
integer :: nlines, niso, i, h
double precision, allocatable :: tab_r(:), tab_starrho(:), tab_T(:), tab_g(:),  tab_mencl(:), phi(:), &
	tab_mfr(:,:), eps(:), tab_np(:), n_nuc(:,:), AtomicNumber(:), m_nuc(:), sigma_N(:), n_x(:)

call get_nlines("/home/luke/summer_2020/mesa/captngen/solarmodels/model_agss09ph_nohead.dat", nlines)

niso = 29

allocate(tab_r(nlines))
allocate(tab_starrho(nlines))
allocate(tab_np(nlines))
allocate(tab_T(nlines))
allocate(tab_g(nlines))
allocate(tab_mencl(nlines))
allocate(phi(nlines))
allocate(tab_mfr(nlines, niso))
allocate(eps(nlines))
allocate(n_nuc(niso, nlines))
allocate(AtomicNumber(niso))
allocate(m_nuc(niso))
allocate(sigma_N(niso))
allocate(n_x(nlines))

call get_solar_params("/home/luke/summer_2020/mesa/captngen/solarmodels/model_agss09ph_nohead.dat", &
nlines, tab_r, tab_starrho, tab_T, tab_g, tab_mencl, phi, tab_mfr)

!print *, "tab_starrho = ", tab_starrho

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Source of error?
! Convert to SI
tab_mencl = tab_mencl*1.98847d30	! Convert to kg
tab_r = tab_r*6.96340d8	! Convert to metres
tab_np = tab_starrho/(1.6726219d-24)*(1.d6)	! Convert to number density in m^-3

! Convert mass fraction to ni
AtomicNumber  = (/ 1., 4., 3., 12., 13., 14., 15., 16., 17., &
                      18., 20.2, 22.99, 24.3, 26.97, 28.1, 30.97,32.06, 35.45, &
                      39.948, 39.098, 40.08, 44.95, 47.86, 50.94, 51.99, &
                      54.93, 55.845, 58.933, 58.693/)
mp = 1.6726219d-27
mx = mp*(5./4.)
m_nuc = AtomicNumber*mp
sigma_0 = 1.d-32
do i=1,29
sigma_N(i) = AtomicNumber(i)**4*(mx+mp)**2/(mx+AtomicNumber(i)*mp)**2 !not yet multiplied by sigma_0
n_nuc(i,:) = tab_mfr(:,i)*tab_starrho(:)/AtomicNumber(i)/mp
enddo
sigma_N = sigma_N*sigma_0

Tx = 8000000
h = int(nlines/2) ! half
!print*, Tx, tab_r(h), tab_T(h), phi(h), tab_starrho(h), mx, n_nuc(1,h), m_nuc(1), sigma_N(1), nlines, niso
!print*, "Etrans = ", Etrans_nl(Tx, tab_T, phi, tab_starrho, mx, n_nuc, m_nuc, sigma_N, nlines, niso)

! Newton's method
guess_1 = 12000000
guess_2 = 12100000
tolerance = 1.0d-8
Tx = newtons_meth(Tx_integral, tab_r, tab_T, phi, tab_starrho, mx, n_nuc, m_nuc, sigma_N, nlines, niso, &
	guess_1, guess_2, tolerance)
print *, "The wimp temperature according to Newton: Tx = ", Tx
print *, "Newton's integral: ", Tx_integral(Tx, tab_r, tab_T, phi, tab_starrho, mx, n_nuc, m_nuc, sigma_N, &
	nlines, niso)

! Calculate luminosity at r (this is what is fed to MESA)
eps = Etrans_nl(Tx, tab_T, phi, tab_starrho, mx, n_nuc, m_nuc, sigma_N, nlines, niso)
L_xtot = trapz(tab_r, eps*tab_r**2, nlines)
print *, "The total wimp luminosity according to Newton: L_xtot = ", L_xtot

n_0 = 0.4*1.78266192d-27*1.0d6/mx ! DM number density in m^-3
n_x = n_0*exp(mx*phi(nlines)/kB/Tx)*exp(-mx*phi/kB/Tx) ! At the edge of the star, set n_x=local galactic

! Plot the function
open(1, file="sp_int.txt")
do i= 1.125e7,1.175e7,1000
	Tx = dble(i)
	eps = Etrans_nl(Tx, tab_T, phi, tab_starrho, mx, n_nuc, m_nuc, sigma_N, nlines, niso)
	L_xtot = trapz(tab_r, eps*tab_r**2, nlines)
	write(1,*) Tx, L_xtot
enddo
close(1)

! Plot the potential, heat tranfer, and DM density
open(2, file="grav_potential.txt")
eps = Etrans_nl(Tx, tab_T, phi, tab_starrho, mx, n_nuc, m_nuc, sigma_N, nlines, niso)
do i=1,nlines
	write(2,*) tab_r(i), phi(i), eps(i), n_x(i)
enddo
close(2)

end program tests
