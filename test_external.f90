program test_external
use capmod, only: trapz, getnlines, get_solar_params
nonlocalmod

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

end program test_external
