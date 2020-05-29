! Test program to calculate the WIMP temperature according to Spergel and Press

module TWIMP
!use minpack ! I've included all of the minpack functions in minpack.f90
implicit none

contains

function sp_integral(Tx, nlines, mx, mp, np, Tp, phi, r)
implicit none
! Calculates the Tx defining integral, equation (4.10) in Spergel&Press "Effect of hypothetical, WIMPs on E transport..."
! The parameter arrays are: args_scal=(/mp,mx/), args_vec=(/np,Tp,phi,r/)

double precision, intent(in) :: Tx	! The independent variable. All others are params
integer, intent(in) :: nlines
double precision, intent(in) :: mx, mp, np(nlines), Tp(nlines), phi(nlines), r(nlines)
double precision :: integrand(nlines)
double precision :: sp_integral
double precision :: k=1.0 ! Boltzmann constant
integer :: i

integrand = np*sqrt((mp*Tx+mx*Tp)/(mx*mp))*(Tp-Tx)*exp(-mx*phi/(k*Tx))*r**2
sp_integral = trapz(r, integrand, nlines)

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

!        if (trapz .lt. 0.d0) then
!          print *, "negative encountered in trapz: i = ", i
!        end if
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
	! Update x_3 using gradient descent formula
	x_3 = x_2 - f(x_2,nlines,mp,mx,np,Tp,phi,r)*(x_2-x_1)/(f(x_2,nlines,mp,mx,np,Tp,phi,r) & 
	-f(x_1,nlines,mp,mx,np,Tp,phi,r))
	error = abs(x_3-x_2)
	! Update x_1 and x_2
	x_1 = x_2
	x_2 = x_3
enddo
newtons_meth = x_3	! The solution to the nonlinear equation

end function

function Etrans(T_x, T_p, n_x, n_p, m_x, m_p, sigma_x)
implicit none
double precision, intent(in) :: T_x, T_p, n_x, n_p, m_x, m_p, sigma_x
double precision, parameter :: k=1.0, pi=3.1415962
double precision :: Etrans

Etrans = 8*sqrt(2/pi)*n_x*n_p*(m_x*m_p)/((m_x+m_p)**2)*sqrt((m_p*k*T_x+m_x*k*T_p)/(m_x*m_p))*k*(T_p-T_x)

return
end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    double precision, parameter :: GMoverR=1.908e15, Rsun = 69.57d9 
    double precision :: tab_r(nlines), tab_starrho(nlines), tab_T(nlines), tab_g(nlines), &
	tab_mencl(nlines), phi(nlines), tab_mfr(nlines, 29)
    
    print *, "nlines = ", nlines
    print *, "tab_r(1967) = ", tab_r(1967)
    print *, "tab_r(1968) = ", tab_r(1968)

    !now actually read in the file
    open(99,file=filename)
    do i=1,nlines
    if (i==1) print *, "Loop started"
    if (i>=nlines) print *, "i>=nlines = ", i
    read(99,*) tab_mencl(i),tab_r(i), tab_T(i), tab_starrho(i), Pres, Lumi, tab_mfr(i,:)
    end do
    print *, "Loop ended"
    close(99)
    
	phi(nlines) = -GMoverR
    do i = 1,nlines-1
    j = nlines-i !trapezoid integral
    phi(j) = phi(j+1) + GMoverR*(tab_r(j)-tab_r(j+1))/2.*(tab_mencl(j)/tab_r(j)**2+tab_mencl(j+1)/tab_r(j+1)**2)
    end do
    tab_g(nlines) = tab_g(nlines-1)
    return
    end

!end get_solar_params
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module TWIMP



program tests
use TWIMP
implicit none

double precision :: Tx, mp, mx, integral, xy_int, guess_1, guess_2, tolerance
integer :: nlines, i
double precision, allocatable :: tab_r(:), tab_starrho(:), tab_T(:), tab_g(:), tab_mencl(:), phi(:), tab_mfr(:,:)
!double precision :: np(1000), Tp(1000), phi(1000), r(1000), x(1000), y(1000)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Make sure sp_integral works 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!call get_nlines("/home/luke/summer_2020/mesa/captngen/solarmodels/model_agss09ph_nohead.dat", nlines)

!mp = 
!mx = 0.2
!Tx = 1.0

!guess_1 = 10 
!guess_2 = 20
!tolerance = 1.0e-6


!! Set everything to 1 to make sure function works
!do i=1,nlines
!	np(i) = 1.0/i
!	Tp(i) = 1.0e4
!	phi(i) = 1.0*i
!	r(i) = dble(i)/nlines
!!	print *, i, r(i), Tp(i)
!	cycle
!enddo

!print *, "Values assigned."

!integral = sp_integral(Tx, nlines, mx, mp, np, Tp, phi, r)
!print *, "First integral evaluated. Value: ", integral, ".", "Analytical result:", -1/exp(1.0)/3, "."


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Test Newton's Method 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Plot the function
!open(1, file="sp_int.txt")

!do i=0,1000,1
!	Tx = dble(i)
!	write(1,*) Tx, sp_integral(Tx, nlines, mx, mp, np, Tp, phi, r)
!enddo
!close(1)

!Tx = newtons_meth(sp_integral, nlines, mx, mp, np, Tp, phi, r, guess_1, guess_2, tolerance)
!print *,  "The result for Tx: ", Tx

call get_nlines("/home/luke/summer_2020/mesa/captngen/solarmodels/model_agss09ph_nohead.dat", nlines)
print *, "The number of lines in the solar model file:", nlines

allocate(tab_r(nlines))
allocate(tab_starrho(nlines))
allocate(tab_T(nlines))
allocate(tab_g(nlines))
allocate(tab_mencl(nlines))
allocate(phi(nlines))
allocate(tab_mfr(nlines, 29))

call get_solar_params("/home/luke/summer_2020/mesa/captngen/solarmodels/model_agss09ph_nohead.dat", &
nlines, tab_r, tab_starrho, tab_T, tab_g, tab_mencl, phi, tab_mfr)

tab_r = tab_r*6.96340d8	! Convert to metres
tab_starrho = tab_starrho/(1.6726219d-24)	! Convert to number density
phi = phi/1.d4	! convert to SI
mp = 1.2726219d-27
mx = mp*(5./4.)

guess_1 = 1.d8
guess_2 = 2.d8
Tx = newtons_meth(sp_integral, nlines, mx, mp, tab_starrho, tab_T, phi, tab_r, guess_1, guess_2, tolerance)

print *, "The wimp temperature: T_x = ", Tx

end program tests
