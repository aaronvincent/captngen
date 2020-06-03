! Test program to calculate the WIMP temperature according to Spergel and Press

module TWIMP
implicit none

contains


function sp_integral(Tx, nlines, mx, mp, np, Tp, phi, r)
implicit none
! Calculates the Tx defining integral, equation (4.10) in Spergel&Press "Effect of hypothetical, WIMPs on E transport..."

double precision, intent(in) :: Tx	! The independent variable. All others are params
integer, intent(in) :: nlines
double precision, intent(in) :: mx, mp, np(nlines), Tp(nlines), phi(nlines), r(nlines)
double precision :: integrand(nlines)
double precision :: sp_integral
double precision :: k=1.38064852d-23 ! Boltzmann constant

integrand = np*sqrt((mp*Tx+mx*Tp)/(mx*mp))*(Tp-Tx)*exp(-mx*phi/(k*Tx))*r**2
sp_integral = trapz(r, integrand, nlines)

end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Source of error?
function Etrans(T_x, T_p, n_x, n_p, m_x, m_p, sigma_x, nlines)
implicit none
integer, intent(in) :: nlines
double precision, intent(in) :: T_x, m_x, m_p, sigma_x
double precision, intent(in) :: n_x(nlines), n_p(nlines), T_p(nlines)
double precision, parameter :: k=1.38064852d-23, pi=3.1415962
double precision :: Etrans(nlines)

Etrans = 8*sqrt(2/pi)*n_x*n_p*sigma_x*(m_x*m_p)/((m_x+m_p)**2)*sqrt((m_p*T_x+m_x*T_p)/(m_x*m_p))*k**(3/2)*(T_p-T_x)

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

function newtons_meth(f, nlines, mx, mp, np, Tp, phi, r, guess_1, guess_2, tolerance)
!Performs gradient descent to solve sp_int=0.
implicit none

double precision :: f ! sp_integral
double precision, intent(in) :: tolerance, guess_1, guess_2
integer, intent(in) :: nlines
double precision, intent(in) :: mx, mp ! 2 scalar parameters
! 4 vector parameters
double precision, intent(in) :: np(nlines), Tp(nlines), phi(nlines), r(nlines) 
double precision :: x_1, x_2, x_3, error
double precision :: newtons_meth

x_1 = guess_1
x_2 = guess_2
error = tolerance + 1

do while (error > tolerance)
	! Update x_3 using Newton's method formula
	x_3 = x_2 - f(x_2,nlines,mx,mp,np,Tp,phi,r)*(x_2-x_1)/(f(x_2,nlines,mx,mp,np,Tp,phi,r) & 
	-f(x_1,nlines,mx,mp,np,Tp,phi,r))
	error = abs(x_3-x_2)
	x_1 = x_2
	x_2 = x_3
enddo

newtons_meth = x_3	! The solution to the nonlinear equation

end function


function binary_search(f, nlines, mx, mp, np, Tp, phi, r, guess_1, guess_2, tolerance)
double precision :: f ! sp_integral
double precision, intent(in) :: tolerance, guess_1, guess_2
integer, intent(in) :: nlines
double precision, intent(in) :: mx, mp ! 2 scalar parameters
! 4 vector parameters
double precision, intent(in) :: np(nlines), Tp(nlines), phi(nlines), r(nlines) 
double precision :: x_1, x_2, x_mid, error, f_1, f_2, f_mid 
double precision :: binary_search

x_1 = guess_1
x_2 = guess_2
error = tolerance + 1

do while (error > tolerance)
x_mid = (x_1+x_2)/2
f_1 = f(x_1,nlines,mx,mp,np,Tp,phi,r)
f_2 = f(x_2,nlines,mx,mp,np,Tp,phi,r)
f_mid = f(x_mid,nlines,mx,mp,np,Tp,phi,r)

!print *, x_1, x_2, f_1, f_2
if (f_mid*f_1 > 0) then ! If f1, fmid are the same sign, set x1 to xmid
	x_1 = x_mid
endif
if (f_mid*f_2 > 0) then ! If f2, fmid are the same sign, set x2 to xmid
	x_2 = x_mid
endif
error = abs(x_2-x_1) ! error=interval
enddo

binary_search = (x_1+x_2)/2

end function


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


end module TWIMP



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Testing program
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program tests
use TWIMP
implicit none

double precision :: Tx, Tx_b, mp, mx, integral, xy_int, guess_1, guess_2, tolerance, sigma_x, L_xtot, n_0
double precision, parameter :: kB = 1.38064852d-23
integer :: nlines, i
double precision, allocatable :: tab_r(:), tab_starrho(:), nx(:), tab_T(:), tab_g(:), tab_mencl(:), phi(:), &
	tab_mfr(:,:), eps(:), tab_np(:)

call get_nlines("/home/luke/summer_2020/mesa/captngen/solarmodels/model_agss09ph_nohead.dat", nlines)

allocate(tab_r(nlines))
allocate(tab_starrho(nlines))
allocate(tab_np(nlines))
allocate(tab_T(nlines))
allocate(tab_g(nlines))
allocate(tab_mencl(nlines))
allocate(phi(nlines))
allocate(tab_mfr(nlines, 29))
allocate(eps(nlines))
allocate(nx(nlines))

call get_solar_params("/home/luke/summer_2020/mesa/captngen/solarmodels/model_agss09ph_nohead.dat", &
nlines, tab_r, tab_starrho, tab_T, tab_g, tab_mencl, phi, tab_mfr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Source of error?
! Convert to SI
tab_mencl = tab_mencl*1.98847d30	! Convert to kg
tab_r = tab_r*6.96340d8	! Convert to metres
tab_np = tab_starrho/(1.6726219d-24)*(1.d6)	! Convert to number density in m^-3

mp = 1.6726219d-27
mx = mp*(5./4.)
sigma_x = 1.0d-38

! Newton's method
guess_1 = 12000000
guess_2 = 12100000
tolerance = 1.0d-8
Tx = newtons_meth(sp_integral, nlines, mx, mp, tab_np, tab_T, phi, tab_r, guess_1, guess_2, tolerance)
print *, "The wimp temperature according to Newton: Tx = ", Tx
print *, "Newton's integral: ", sp_integral(Tx, nlines, mx, mp, tab_np, tab_T, phi, tab_r)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Source of error?
! Calculate wimp number density in sun
n_0 = 0.4*1.78266192d-27/mx*1.0d6 	! Local wimp density in m^-3
nx = n_0*exp(-mx*phi/kB/Tx)*exp(mx*phi(nlines)/(kB*Tx))

! Calculate luminosity at r (this is what is fed to MESA)
eps = Etrans(Tx, tab_T, nx, tab_np, mx, mp, sigma_x, nlines)
L_xtot = trapz(tab_r, eps*tab_r**2, nlines)
print *, "The total wimp luminosity according to Newton: L_xtot = ", L_xtot

! Plot the function
open(1, file="sp_int.txt")
do i= 1.225e7,1.3e7,1000
	Tx = dble(i)
	eps = Etrans(Tx, tab_T, nx, tab_np, mx, mp, sigma_x, nlines)
	L_xtot = trapz(tab_r, eps*tab_r**2, nlines)
	write(1,*) Tx, L_xtot
enddo
close(1)

! Plot the potential
open(2, file="grav_potential.txt")
do i=1,nlines
	write(2,*) tab_r(i), phi(i)
enddo
close(2)

open(3, file="wimp_density.txt")
do i=1,nlines
	write(3,*) tab_r(i), nx(i)
enddo
close(3)

end program tests
