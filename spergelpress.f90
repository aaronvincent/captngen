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

!subroutine s_p(n, Tx, intg, iflag, mp, mx, np, Tp, phi, r, nlines)
! Wierd formatting to agree with minpack
!integer(kind = 4) :: n=1
!real(kind = 8) :: intg(n)
!integer(kind = 4) :: iflag
!real(kind = 8) :: Tx(n)

!fvec = sp_integral(x, 1, 1)

!end subroutine



!function wimp_temp(r, integrand, nlines)

! integrand must be an array dependent only on Tx
! r, integrand are 1d arrays of length nlines
!implicit none

!integer :: nlines 
!double precision :: r(nlines), integrand(nlines)
!double precision :: lhs, Tx=1.d4, tol=1.d-6, info=0, fvec(1)	! Initial guess for Tx is 10^4
!integer :: lwa=8	! Length of wa - required by hybrd1
!double precision :: wa(lwa)	! Also required by hybrd1
!double precision :: wimp_temp

!lhs = trapz(r, integrand, flen)	! We want Tx such that lhs=0

!call hybrd1(lhs, 1, Tx, fvec, tol, info, wa, lwa)	! MINPACK function to solve integral(integrand(Tx))=0
!wimp_temp = Tx

!end function wimp_temp

end module TWIMP

program tests
use TWIMP
implicit none

double precision :: Tx, mp, mx, integral, xy_int, guess_1, guess_2, tolerance
integer :: nlines, i
double precision :: np(1000), Tp(1000), phi(1000), r(1000), x(1000), y(1000)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Make sure sp_integral works
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

nlines = 1000
mp = 0.1
mx = 0.2
!Tx = 1.0

guess_1 = 10 
guess_2 = 20
tolerance = 1.0e-6

print *, "Variables declared."

! Set everything to 1 to make sure function works
do i=1,nlines
	np(i) = 1.0/i
	Tp(i) = 1.0e4
	phi(i) = 1.0*i
	r(i) = dble(i)/nlines
!	print *, i, r(i), Tp(i)
	cycle
enddo

print *, "Values assigned."

!integral = sp_integral(Tx, nlines, mx, mp, np, Tp, phi, r)
!print *, "First integral evaluated. Value: ", integral, ".", "Analytical result:", -1/exp(1.0)/3, "."


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Test Gauss-Newton solver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Plot the function to get an idea
open(1, file="sp_int.txt")

do i=0,1000,1
	Tx = dble(i)
	write(1,*) Tx, sp_integral(Tx, nlines, mx, mp, np, Tp, phi, r)
enddo
close(1)

Tx = newtons_meth(sp_integral, nlines, mx, mp, np, Tp, phi, r, guess_1, guess_2, tolerance)

print *,  "The result for Tx: ", Tx

end program tests
