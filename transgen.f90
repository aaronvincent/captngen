!!!!!! TRANSPORTER GENERAL !!!!!
!!! Asymmetric dark matter transport routine, check out https://arxiv.org/pdf/1311.2074.pdf
!!! The zetas in Eq. 31 should not be there
!!! for constant, q- and v- dependent cross sections
!!! Uses capmod from capgen.f90

!Input:
!nwimps: Total number of DM particles in the star. I know ADM is not WIMPs, stop complaining
!niso: number of isotopes: 1 = spin-dependent
!nq, nv: v^n, q^n numberwang

!dm properties are set when you call capgen.


!Output
!Etrans erg/g/s (I think)

subroutine transgen(sigma_0,Nwimps,niso,nonlocal,Tx,noise_indicator,etrans,EtransTot)
!mdm is stored in capmod
! Tx is unchanged in the LTE scheme, and is the output one-zone WIMP temp in the nonlocal scheme
use capmod
use akmod
use nonlocalmod
implicit none
!nlines might be redundant
logical, intent(in) :: nonlocal
logical splinelog, DCT !for PCHIP
integer, intent(in):: niso
double precision, intent(in) :: sigma_0, Nwimps
double precision, intent(out) :: noise_indicator
integer, parameter :: decsize = 75 !this should be done a bit more carefully
integer i, ri, ierr
integer (kind=4) :: lensav 
double precision :: epso,EtransTot
double precision, parameter :: GN = 6.674d-8, kBeV=8.617e-5 ! kB and mnucg defined in nonlocalmod
double precision :: mxg, rchi, Tc,rhoc,K, integrand
double precision :: capped, maxcap !this is the output
double precision :: phi(nlines), Ltrans(nlines),Etrans(nlines),mfp(nlines),nabund(niso,nlines),sigma_N(niso)
double precision :: nx(nlines),alphaofR(nlines), kappaofR(nlines),cumint(nlines),cumNx,nxIso(nlines),cumNxIso, n_0
double precision :: r_even(nlines), T_even(nlines), dTdr_even(nlines), work(2*nlines) ! Evenly spaced arrays for FFT
! More evenly spaced arrays for FFT
double precision :: L_even(nlines), dLdr_even(nlines), Etrans_even(nlines), Etrans_test(nlines), Ltrans_cond(nlines)
double precision :: r_double(2*nlines), dTdr_mirror(2*nlines-2), r_mesa(1999)
double precision :: muarray(niso),alpha(niso),kappa(niso),dphidr(nlines),dTdr(nlines)
double precision :: fgoth, hgoth(nlines), ggoth_mesa(1999), ggoth(nlines), dLdR(nlines),isplined1,dLdRscratch(nlines)
double precision :: dggothdr_mesa(1999), dggothdr(nlines), test_array(2000)
double precision :: biggrid(nlines), bcoeff(nlines), ccoeff(nlines), dcoeff(nlines) ! for spline
integer lwk !for pchip
double precision :: pchipScratch(3*nlines) !for pchip
double precision, allocatable :: wsave(:)

double precision :: brcoeff(nlines), crcoeff(nlines), drcoeff(nlines) ! for spline
double precision :: bdcoeff(decsize), cdcoeff(decsize), ddcoeff(decsize) ! for spline
double precision :: smallgrid(decsize), smallR(decsize), smallT(decsize), smallL(decsize),smalldL(decsize),smalldT(decsize),ispline
double precision :: Tx, guess_1, guess_2, tolerance ! For the Spergel & Press nonlocal scheme

lwk = 3*nlines !This is the length of pchipScratch. don't redefine this without also changing pchipScratch
epso = tab_r(2)/10.d0 ! small number to prevent division by zero
! smallr = (/((i*1./dble(decsize-1)),i=1,decsize)/) - 1./dble(decsize-1)
smallgrid =  (/((i*1./dble(decsize-1)),i=1,decsize)/) - 1./dble(decsize-1) !(/i, i=1,decsize /)
biggrid =  (/((i*1./dble(nlines-1)),i=1,nlines)/) - 1./dble(nlines-1) !(/i, i=1,nlines/)

mxg = mdm*1.78d-24
Tc = tab_T(1)
rhoc = tab_starrho(1)

if (decsize .ge. nlines) stop "Major problem in transgen: your low-res size is larger than the original"
!Check if the stellar parameters have been allocated
if (.not. allocated(tab_r)) stop "Error: stellar parameters not allocated in transgen"


!set up extra stellar arrays that we need
phi = - tab_vesc**2/2.d0
dphidr = -tab_g

alphaofR(:) = 0.d0
kappaofR(:) = 0.d0

do i = 1,niso
    !this is fine for SD as long as it's just hydrogen. Otherwise, spins must be added
    muarray(i) = mdm/AtomicNumber(i)/mnuc
    sigma_N(i) = AtomicNumber(i)**4*(mdm+mnuc)**2/(mdm+AtomicNumber(i)*mnuc)**2 !not yet multiplied by sigma_0
    nabund(i,:) = tab_mfr(:,i)*tab_starrho(:)/AtomicNumber(i)/mnucg
    !these shouldn't really be done every iteration, can fix later
    call interp1(muVect,alphaVect,nlinesinaktable,muarray(i),alpha(i))
    call interp1(muVect,kappaVect,nlinesinaktable,muarray(i),kappa(i))
end do

! PS: I've commented out the following line -- the array bounds don't match, so it
! creates memory corruption(!) It also looks like it is only here by accident...
!    alphaofR = alphaofR/(sigma_N*sum(nabund,1))

!compute mean free path
if ((nq .eq. 0) .and. (nv .eq. 0)) then
  do i = 1,nlines
    mfp(i) = 1/sum(sigma_N*nabund(:,i))/sigma_0/2. !factor of 2 b/c  sigma_tot = 2 sigma_0
  end do
! else if ((nq .eq. )) !q, v dependence goes here
end if

rchi = (3.*(kB*Tc)/(2.*pi*GN*rhoc*mxg))**.5;

K = mfp(1)/rchi;

!T derivatives
! smooth T
! some gymnastics are necessary, because the temperature is not smooth at all
! a simple spline -> derivative doesn't help. Fourier method works better
! Take dT/dr with pchip
splinelog = .false.
call DPCHEZ( nlines, tab_r, tab_T, dTdr, SPLINElog, pchipScratch, LWK, IERR )
if (ierr .lt. 0) then
	print*, 'DPCHEZ interpolant failed with error ', IERR
	return
ENDIF

if (any(isnan(dTdr))) print *, "NAN encountered in dT/dr"

! Smooth dTdr with FFT
! First build evenly spaced r and dTdr arrays
do i=1,nlines
	r_even(i) =  i*1./dble(nlines)
enddo

lensav = nlines + int(log(real(nlines))) + 4 ! Minimum length required by fftpack

! Cut out high frequency components. The subroutine fourier_smooth is located in nonlocalmod.f90
! Keep lowest 5% of components
call fourier_smooth(tab_r, dTdr, r_even, dTdr_even, 0.05d0, noise_indicator, nlines, lensav, ierr)
dTdr = dTdr/Rsun

!this loop does a number of things
cumint(1) = 0.d0
cumNx = 0.d0

do i = 1,nlines

! 1) get alpha & kappa averages
  alphaofR(i) = sum(alpha*sigma_N*nabund(:,i))/sum(sigma_N*nabund(:,i))
  kappaofR(i) = mfp(i)*sum(sigma_0*sigma_N*nabund(:,i)/kappa)
  kappaofR(i) = 1./kappaofR(i)
  !perform the integral inside the nx integral

  integrand = (kB*alphaofR(i)*dTdr(i) + mxg*dphidr(i))/(kB*tab_T(i))

  ! print*, alphaofR(i),mxg
  if (i > 1) then
  cumint(i) = cumint(i-1) + integrand*tab_dr(i)*Rsun
  end if

  nx(i) = (tab_T(i)/Tc)**(3./2.)*exp(-cumint(i))

  cumNx = cumNx + 4.*pi*tab_r(i)**2*tab_dr(i)*nx(i)*Rsun**3.

  nxIso(i) = Nwimps*exp(-Rsun**2*tab_r(i)**2/rchi**2)/(pi**(3./2.)*rchi**3) !normalized correctly
end do

nx = nx/cumNx*nwimps !normalize density
fgoth = 1./(1.+(K/.4)**2)
hgoth = ((tab_r*Rsun - rchi)/rchi)**3 +1.
hgoth(1) = 0.d0 !some floating point shenanigans.

nx = fgoth*nx + (1.-fgoth)*nxIso


if (nonlocal .eqv. .false.) then ! if nonlocal=false, use Gould & Raffelt regime to calculate transport

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Gould Raffelt section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Ltrans = 4.*pi*(tab_r+epso)**2.*Rsun**2.*kappaofR*fgoth*hgoth*nx*mfp*sqrt(kB*tab_T/mxg)*kB*dTdr;

if (any(isnan(tab_r))) print *, "NAN encountered in tab_r"
if (any(isnan(kappaofR))) print *, "NAN encountered in kappa"
if (any(isnan(hgoth))) print *, "NAN encountered in hgoth"
if (any(isnan(nx))) print *, "NAN encountered in nx"
if (any(isnan(tab_T))) print *, "NAN encountered in tab_T"
if (any(isnan(Ltrans))) print *, "NAN encountered in Ltrans"

!get derivative of luminosity - also needs smoothing. Fourier method doesn't work as well here
!Hermite polynomial interpolation seems to work ok.

!Cubic hermite polynomial to smooth dLdr
splinelog = .false.
call DPCHEZ( nlines, tab_r, Ltrans, dLdR, SPLINElog, pchipScratch, LWK, IERR )
if (ierr .lt. 0) then
  print*, 'DPCHEZ interpolant failed with error ', IERR
  return
ENDIF
dLdr = dLdr/Rsun

Etrans = 1./(4.*pi*(tab_r+epso)**2*tab_starrho)*dLdR/Rsun**2
! Etrans_test is to check how the noise in Etrans. noise_indicator is the sum of the components above the cutoff
Etrans_test = Etrans
call fourier_smooth(tab_r, Etrans_test, r_even, dTdr_even, 0.05d0, noise_indicator, nlines, lensav, ierr)
EtransTot = trapz(tab_r,abs(dLdR)*Rsun,nlines)
print *, "Transgen: total G&R transported energy = ", EtransTot
print *, "fogth=", fgoth

! Useful info to have when troubleshooting, but make sure to change the file path!
! Check Ltrans
open(55,file = "/home/luke/summer_2020/mesa/test_files/Ltrans_gr.dat")
do i=1,nlines
	write(55,*) tab_r(i), Ltrans(i), Etrans(i), kappaofR(i), mfp(i), tab_T(i), dTdR(i), hgoth(i), dLdR(i), nx(i), &
		tab_starrho(i)
end do
close(55)
open(55, file="/home/luke/summer_2020/mesa/test_files/Lmax_gr.dat", access="APPEND")
write(55,*) mfp(1), maxval(-Ltrans)
!close(55)

return


else if (nonlocal) then ! if nonlocal=true, use Spergel & Press regime to calculate heat transport

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Spergel Press section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The nonlocal transport scheme: articles.adsabs.harvard.edu/pdf/1985ApJ...294..663S
! The functions of interest are in nonlocalmod.f90. These also use https://arxiv.org/pdf/0809.1871.pdf

! One-zone WIMP temp guesses in K. They both have to be either greater than or less than the actual
! Tx, so I just hard set them here.
guess_1 = 1.0d7 ! ***** A possible source of error for exotic situations (eg large mass stars) *****
guess_2 = 1.01d7
tolerance = 1.0d-1

! Tx is the Spergel & Press one-zone WIMP temperature in Kelvin

Tx = newtons_meth(Tx_integral, dTdr, mfp, sigma_N, sigma_0, alpha, kappa, Nwimps, niso, guess_1, guess_2, tolerance)

nx = nx_func(Tx, dTdr, mfp, sigma_N, sigma_0, alpha, kappa, Nwimps, niso)

! Etrans in erg/g/s (according to Spergel Press)
Etrans = Etrans_nl(Tx, dTdr, mfp, sigma_N, sigma_0, alpha, kappa, Nwimps, niso) ! erg/g/s

! Calculate Ltrans
do i=1,nlines
	Ltrans(i) = trapz(tab_r*Rsun, 4.d0*pi*(tab_r*Rsun)**2.d0*Etrans*tab_starrho, i)
enddo

!! Gilliland interpolation articles.adsabs.harvard.edu/pdf/1986ApJ...306..703G 
!Ltrans_cond = 4.d0*pi*(tab_r*Rsun)**2.d0*nxIso*sqrt(kB*tab_T/mxg)*kB*dTdr/&
!	(nabund(1,:)*sigma_0*sqrt((3.d0*mxg+mnucg)/(mxg+mnucg)))
!Ltrans = (Ltrans_cond*Ltrans)/(Ltrans_cond+Ltrans+epso)

! If you develop a conversion factor to go from SP --> GR, the following section will implement it.
!------------------------------------------------------------------------------
!! Now convert Etrans_sp to Etrans_gr (by multiplying Ltrans by ggoth, a conversion factor). 
!! Easier to work with Ltrans than Etrans because Ltrans is always positive in both schemes.
!! Load ggoth
!open(55, file="/home/luke/summer_2020/mesa/test_files/ggoth.dat")
!do i = 1, 1999
!	read(55, *) r_mesa(i), ggoth_mesa(i)
!enddo
!! Interpolate ggoth to MESA grid (cubic spline)
!splinelog = .false.
!call DPCHEZ(1999, r_mesa, ggoth_mesa, dggothdr_mesa, SPLINElog, pchipScratch, LWK, IERR )
!if (ierr .lt. 0) then
!	print*, 'DPCHEZ interpolant failed with error ', IERR
!	return
!ENDIF
!call DPCHEV (1999, r_mesa, ggoth_mesa, dggothdr_mesa, nlines, tab_r, ggoth, dggothdr, IERR )
!! Convert Ltrans_sp --> Ltrans_gr with ggoth
!!Ltrans = Ltrans*ggoth
!do i=1,nlines-1
!	dLdr(i) = (Ltrans(i+1)-Ltrans(i))/(tab_r(i+1)-tab_r(i))/Rsun
!enddo
!Etrans = 1./(4.*pi*(tab_r+epso)**2*tab_starrho)*dLdR/Rsun**2
!------------------------------------------------------------------------------

! The total WIMP transported energy (erg/s). In the S&P scheme, this should be 0 by definition of Tx.
EtransTot = trapz(tab_r*Rsun, 4.d0*pi*(tab_r*Rsun)**2*Etrans*tab_starrho, nlines)
print *, "Transgen: total S&P transported energy = ", EtransTot

! Check the noise in Etrans (same as GR scheme) - note that dTdr_even is just a work array, I'm not actually using 
! anything about dTdr in the noise calculation
Etrans_test = Etrans
call fourier_smooth(tab_r, Etrans_test, r_even, dTdr_even, 0.05d0, noise_indicator, nlines, lensav, ierr)

! Write things to file
open(55,file = "/home/luke/summer_2020/mesa/test_files/Ltrans_sp.dat")
do i=1,nlines
	write(55,*) tab_r(i), Ltrans(i), Etrans(i), nx(i) , tab_T(i), tab_g(i), dTdr(i), nxIso(i), nabund(1,i)
end do
close(55)
!open(55, file="/home/luke/summer_2020/mesa/test_files/Lmax_sp.dat", access="APPEND")
!write(55,*) mfp(1), maxval(abs(Ltrans)), sigma_0
!close(55)

return


endif

end subroutine transgen
