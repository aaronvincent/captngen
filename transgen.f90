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

subroutine transgen(sigma_0,Nwimps,niso,nonlocal,etrans,EtransTot)
!mdm is stored in capmod
use capmod
use akmod
use nonlocalmod
implicit none
!nlines might be redundant
logical, intent(in) :: nonlocal
integer, intent(in):: niso
double precision, intent(in) :: sigma_0, Nwimps
integer, parameter :: decsize = 180 !this should be done a bit more carefully
integer i, ri
double precision :: epso,EtransTot
double precision, parameter :: GN = 6.674d-8, kB = 1.3806d-16,kBeV=8.617e-5,mnucg=1.67e-24
double precision :: mxg, rchi, Tc,rhoc,K, integrand
double precision :: capped, maxcap !this is the output
double precision :: phi(nlines), Ltrans(nlines),Etrans(nlines),mfp(nlines),nabund(niso,nlines),sigma_N(niso)
double precision :: nx(nlines),alphaofR(nlines), kappaofR(nlines),cumint(nlines),cumNx,nxIso(nlines),cumNxIso
double precision :: muarray(niso),alpha(niso),kappa(niso),dphidr(nlines),dTdr(nlines)
double precision :: fgoth, hgoth(nlines), dLdR(nlines),isplined1
double precision :: biggrid(nlines), bcoeff(nlines), ccoeff(nlines), dcoeff(nlines) ! for spline
double precision :: brcoeff(nlines), crcoeff(nlines), drcoeff(nlines) ! for spline
double precision :: bdcoeff(decsize), cdcoeff(decsize), ddcoeff(decsize) ! for spline
double precision :: smallgrid(decsize), smallR(decsize), smallT(decsize), smallL(decsize),smalldL(decsize),smalldT(decsize),ispline
double precision :: Tx, guess_1, guess_2, tolerance ! For the Spergel & Press nonlocal schemeF

epso = tab_r(2)/10.d0 ! small number to prevent division by zero
! smallr = (/((i*1./dble(decsize-1)),i=1,decsize)/) - 1./dble(decsize-1)
smallgrid =  (/((i*1./dble(decsize-1)),i=1,decsize)/) - 1./dble(decsize-1) !(/i, i=1,decsize /)
biggrid =  (/((i*1./dble(nlines-1)),i=1,nlines)/) - 1./dble(nlines-1) !(/i, i=1,nlines/)

mxg = mdm*1.78d-24
Tc = tab_T(1)
rhoc = tab_starrho(1)
! print*, "Nwimps in ", Nwimps

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


if (nonlocal .eqv. .false.) then ! if nonlocal=false, use Gould & Raffelt regime to calculate transport

!compute mean free path
if ((nq .eq. 0) .and. (nv .eq. 0)) then
  do i = 1,nlines
    mfp(i) = 1/sum(sigma_N*nabund(:,i))/sigma_0/2. !factor of 2 b/c  sigma_tot = 2 sigma_0
  end do
! else if ((nq .eq. )) !q, v dependence goes here
end if

rchi = (3.*(kB*Tc)/(2.*pi*GN*rhoc*mxg))**.5;

K = mfp(1)/rchi;


!smooth T
!some gymnastics are necessary, because the temperature is not smooth at all
!first build a cubic spline fit
call spline(biggrid, tab_R, brcoeff, crcoeff, drcoeff, nlines)
call spline(tab_R, tab_T, bcoeff, ccoeff, dcoeff, nlines)

!now build a lower resolution array: this effectively smooths to relevant scales
!the smallR is to ensure the adaptive grid is preserved
do i= 1,decsize
smallR(i) = ispline(smallgrid(i),biggrid,tab_R,brcoeff, crcoeff, drcoeff, nlines)
smallT(i) = ispline(smallr(i),tab_R,tab_T,bcoeff,ccoeff,dcoeff,nlines)
end do
call sgolay(smallT,decsize,4,1,smalldT) !differentiate
smalldT(decsize) = 0.d0
smalldT(1) = 0.d0
call spline(smallR, smalldT, bdcoeff, cdcoeff, ddcoeff, decsize) !spline for derivative
!Re-expand to the full array size
do i= 1,nlines
dTdR(i) = ispline(tab_R(i),smallR,smalldT,bdcoeff,cdcoeff,ddcoeff,decsize)
end do
dTdR(1) = 0.d0
dTdR(nlines) = 0.d0
dTdR = dTdR/Rsun*dble(decsize-1)


! call sgolay(tab_T,nlines,3,0,tab_T)
! call spline(tab_r, tab_T, bcoeff, ccoeff, dcoeff, nlines)
! dTdR = bcoeff/Rsun
! call sgolay(dTdR,nlines,3,0,dTdR) !don't ask
! take derivative (for more fun)
! Get derivative of T
! call sgolay(tab_T,nlines,3,1,dTdr)
! dTdr = dTdr/Rsun/tab_dr


! do i = 2,nlines
!   dTdr(i) = (tab_T(i)-tab_T(i-1))/tab_dr(i) !does this kind of indexing work?
! end do
! dTdr(nlines) = 0.d0


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

  ! print*,nx(i)
  cumNx = cumNx + 4.*pi*tab_r(i)**2*tab_dr(i)*nx(i)*Rsun**3.

  nxIso(i) = Nwimps*exp(-Rsun**2*tab_r(i)**2/rchi**2)/(pi**(3./2.)*rchi**3) !normalized correctly
  ! print*,tab_r(i), nxIso(i)
  ! print*,exp(-Rsun**2*tab_r(i)**2/rchi**2)
end do

nx = nx/cumNx*nwimps !normalize density
fgoth = 1./(1.+(K/.4)**2)
hgoth = ((tab_r*Rsun - rchi)/rchi)**3 +1.
hgoth(1) = 0.d0 !some floating point shenanigans.

! Check nx_LTE
!open(55,file = "/home/luke/summer_2020/mesa/captngen/nx_LTE_change.dat")
! do i=1,nlines
! write(55,*) tab_r(i), nx(i), cumint(i), tab_T(i), dTdR(i), tab_g(i), dphidr(i)
! end do
! close(55)

! nx = nxIso
nx = fgoth*nx + (1.-fgoth)*nxIso

Ltrans = 4.*pi*(tab_r+epso)**2.*Rsun**2*kappaofR*fgoth*hgoth*nx*mfp*sqrt(kB*tab_T/mxg)*kB*dTdr;

!get derivative of luminosity - same nonsense as with the temperature
!I'm going to reuse the temperature array, don't get confused :-)
call spline(tab_R, Ltrans, bcoeff, ccoeff, dcoeff, nlines)
do i= 1,decsize
    smallL(i) = ispline(smallr(i),tab_R,Ltrans,bcoeff,ccoeff,dcoeff,nlines)
end do
call sgolay(smallL,decsize,4,1,smalldL) !Take the derivative
! smalldL(1) = 0.d0
! smalldL(1) = smalldL(2)
smalldL(decsize) = 0.d0
call spline(smallR, smalldL, bdcoeff, cdcoeff, ddcoeff, decsize) !spline for derivative
do i= 1,nlines
  dLdR(i) = ispline(tab_R(i),smallR,smalldL,bdcoeff,cdcoeff,ddcoeff,decsize)
end do

dLdR = dLdR/Rsun*dble(decsize-1)

if (any(abs(dLdR) .gt. 1.d100)) then
  open(55,file = "crashsmallarrays.dat")
  do i=1,decsize
    write(55,*) smallR(i), smallT(i), smalldT(i), smallL(i), smalldL(i)
  write(55,*)
  end do
  close(55)
  stop "Infinite luminosity derivative encountered"

end if



! call sgolay(Ltrans,nlines,4,1,Ltrans)
! call sgolay(Ltrans,nlines,3,1,dLdr)
! ! call spline(tab_r, Ltrans, bcoeff, ccoeff, dcoeff, nlines)
! ! dLdr = bcoeff/Rsun
! ! dLdr = dLdr/Rsun/tab_dr
! ! dLdr(1)= 0.d0
! call sgolay(dLdr,nlines,4,0,dLdr)

! do i = 2,nlines
!   dLdr(i) = (Ltrans(i)-Ltrans(i-1))/tab_dr(i) !does this kind of indexing work?
! end do
!
! dLdr(nlines) = 0.d0
! dLdr = dLdr/Rsun


Etrans = 1./(4.*pi*(tab_r+epso)**2*tab_starrho)*dLdR/Rsun**2;

! print*, Etrans

EtransTot = trapz(tab_r,abs(dLdR),nlines)

<<<<<<< HEAD
! Write Etrans to file
! open(55,file = "/home/luke/summer_2020/mesa/captngen/captranstest.dat")
! do i=1,nlines
! write(55,*) tab_r(i), nx(i), tab_T(i), Ltrans(i), Etrans(i), dTdR(i), dLdR(i), tab_starrho(i), tab_g(i), dphidr(i)
! end do
! close(55)

return

!
! open(55,file = "smallarrays.dat")
! do i=1,decsize
!   write(55,*) smallR(i), smallT(i), smalldT(i), smallL(i), smalldL(i)
! write(55,*)
! end do
! close(55)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Spergel Press section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
else if (nonlocal .eqv. .true.) then ! if nonlocal=true, use Spergel & Press regime to calculate heat transport
! The nonlocal transport scheme: http://articles.adsabs.harvard.edu/pdf/1985ApJ...294..663S
! the functions of interest are in nonlocalmod. These also use https://arxiv.org/pdf/0809.1871.pdf
!print *, "calculating spergel press"
guess_1 = 1.0d7 ! Change these to better guesses later?
guess_2 = 1.01d7
tolerance = 1.0d-4

! Tx is the Spergel & Press one-zone WIMP temperature
Tx = newtons_meth(Tx_integral, tab_r*Rsun, tab_T, phi, tab_starrho, mxg, nabund, AtomicNumber*mnucg, & 
	sigma_0*sigma_N, Nwimps, nlines, niso, guess_1, guess_2, tolerance)

! Etrans in erg/g/s
Etrans = Etrans_nl(Tx, tab_r*Rsun, tab_T, phi, tab_starrho, mxg, nabund, AtomicNumber*mnucg, &
	 sigma_0*sigma_N, Nwimps, nlines, niso) ! erg/g/s
print *, "Transgen: Tx = ", Tx

EtransTot = trapz(tab_r*Rsun, 4*pi*(tab_r*Rsun)**2*Etrans*tab_starrho, nlines)
print *, "Transgen: total transported energy = ", EtransTot

! Write Etrans to file
!open(10, file="/home/luke/summer_2020/mesa/captngen/Etrans_sp.txt")
!do i=1,nlines
!	write(10, *) tab_r(i), Etrans(i), tab_starrho(i), tab_T(i), sum(nabund(:,1))*mnucg
!enddo
!close(10)

return

endif

end subroutine transgen
