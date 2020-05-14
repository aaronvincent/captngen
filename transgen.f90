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


!   Updated 2020 to include q- and v- dependent transport
!   Added mean free paths for each nq,nv=-1,0,1,2 case
!   nq and nv can be read in from darkInputs.txt if using DarkMESA
!   if not, just set them in main

subroutine transgen(Nwimps,niso_in,nq_in,nv_in,etrans,EtransTot)
!mdm is stored in capmod
use capmod
use akmod
implicit none
!nlines might be redundant
integer, intent(in):: nq_in, niso_in, nv_in
double precision, intent(in) :: Nwimps
integer, parameter :: decsize = 180 !this should be done a bit more carefully
integer i, ri
double precision :: epso,EtransTot
double precision, parameter :: GN = 6.674d-8, kB = 1.3806d-16,kBeV=8.617e-5,mnucg=1.67e-24
double precision :: mxg, rchi, Tc,rhoc,K, integrand
double precision :: capped, maxcap !this is the output
double precision :: phi(nlines), Ltrans(nlines),Etrans(nlines),mfp(nlines),nabund(niso_in,nlines),sigma_N(niso_in)
double precision :: thermavg_sigma(nlines), zeta_v(nlines), zeta_q(nlines)
double precision :: nx(nlines),alphaofR(nlines), kappaofR(nlines),cumint(nlines),cumNx,nxIso(nlines),cumNxIso
double precision :: muarray(niso_in),alpha(niso_in),kappa(niso_in),dphidr(nlines),dTdr(nlines)
double precision :: fgoth, hgoth(nlines), dLdR(nlines),isplined1
double precision :: biggrid(nlines), bcoeff(nlines), ccoeff(nlines), dcoeff(nlines) ! for spline
double precision :: brcoeff(nlines), crcoeff(nlines), drcoeff(nlines) ! for spline
double precision :: bdcoeff(decsize), cdcoeff(decsize), ddcoeff(decsize) ! for spline
double precision :: smallgrid(decsize), smallR(decsize), smallT(decsize), smallL(decsize),smalldL(decsize),smalldT(decsize),ispline

epso = tab_r(2)/10.d0 ! small number to prevent division by zero
! smallr = (/((i*1./dble(decsize-1)),i=1,decsize)/) - 1./dble(decsize-1)
smallgrid =  (/((i*1./dble(decsize-1)),i=1,decsize)/) - 1./dble(decsize-1) !(/i, i=1,decsize /)
biggrid =  (/((i*1./dble(nlines-1)),i=1,nlines)/) - 1./dble(nlines-1) !(/i, i=1,nlines/)

! niso = niso
mxg = mdm*1.78d-24
Tc = tab_T(1)
rhoc = tab_starrho(1)
niso = niso_in
nq = nq_in
nv = nv_in

if (decsize .ge. nlines) stop "Major problem in transgen: your low-res size is larger than the original"
!Check if the stellar parameters have been allocated
if (.not. allocated(tab_r)) stop "Error: stellar parameters not allocated in transgen"


!set up extra stellar arrays that we need
phi = - tab_vesc**2/2.d0
dphidr = -tab_g

!remember, nq and nv are set in darkInputs.txt for DarkMESA or if not, just in main
call get_alpha_kappa(nq,nv)
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

  !get weighted alpha and kappa vs r
end do
alphaofR = alphaofR/(sigma_N*sum(nabund,1))

!need separate zeta factors for q- and v- dependent interactions
do i = 1,nlines
  zeta_v(i) = v0/(sqrt(2*kB*tab_T(i)/mxg))
  zeta_q(i) = q0/c0/(mxg*sqrt(2*kB*tab_T(i)/mxg))
end do

! mean free path calcs for each nq,nv case here
! equations from 1311.2074 (eqns 69 to 74 on arxiv copy)
if ((nq .eq. 0) .and. (nv .eq. 0)) then
  do i = 1,nlines
    mfp(i) = 1/sum(sigma_N*nabund(:,i))/sigma_0/2. !factor of 2 b/c  sigma_tot = 2 sigma_0
  end do
else if ((nq .eq. 1)) then
  do i = 1,nlines
    mfp(i) = 1/sum(6.*nabund(:,i)*sigma_N/(1.+muarray)/(zeta_q**2))/sigma_0/2.
  end do
else if ((nq .eq. 2)) then
  do i = 1,nlines
    mfp(i) = 1/sum(40.*nabund(:,i)*sigma_N/((1.+muarray)**2)/(zeta_q**4))/sigma_0/2.
  end do
else if ((nq .eq. -1)) then
  do i = 1,nlines
    mfp(i) = 1/sum(nabund(:,i)*sigma_N*(1.+muarray)*zeta_q**2)/sigma_0/2.
  end do
else if ((nv .eq. 1)) then
  do i = 1,nlines
    mfp(i) = 1/sum(nabund(:,i)*sigma_N*(1.+muarray)*3./2./(zeta_v**2))/sigma_0/2.
  end do
else if ((nv .eq. 2)) then
  do i = 1,nlines
    mfp(i) = 1/sum(nabund(:,i)*sigma_N*((1.+muarray)**2)*15./4./(zeta_v**4))/sigma_0/2.
  end do
else if ((nv .eq. -1)) then
  do i = 1,nlines
    mfp(i) = 1/sum(nabund(:,i)*sigma_N*2*zeta_v**2/(1.+muarray))/sigma_0/2.
  end do
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



!this loop does a number of things
cumint(1) = 0.d0
cumNx = 0.d0
do i = 1,nlines

  !get alpha & kappa averages
  alphaofR(i) = sum(alpha*sigma_N*nabund(:,i))/sum(sigma_N*nabund(:,i))
  kappaofR(i) = mfp(i)*sum(sigma_0*sigma_N*nabund(:,i)/kappa)
  kappaofR(i) = 1./kappaofR(i)

  !perform the integral inside the nx integral
  integrand = (kB*alphaofR(i)*dTdr(i) + mxg*dphidr(i))/(kB*tab_T(i))
  if (i > 1) then
    cumint(i) = cumint(i-1) + integrand*tab_dr(i)*Rsun
  end if

  nx(i) = (tab_T(i)/Tc)**(3./2.)*exp(cumint(i))
  cumNx = cumNx + 4.*pi*tab_r(i)**2*tab_dr(i)*nx(i)*Rsun**3.
  nxIso(i) = Nwimps*exp(-Rsun**2*tab_r(i)**2/rchi**2)/(pi**(3./2.)*rchi**3) !normalized correctly
end do
nx = nx/cumNx*nwimps !normalize density

!These are the interpolating functions used by G&R for transition to LTE regime
fgoth = 1./(1.+(K/.4)**2)
hgoth = ((tab_r*Rsun - rchi)/rchi)**3 +1.
hgoth(1) = 0.d0 !some floating point shenanigans.

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

Etrans = 1./(4.*pi*(tab_r+epso)**2*tab_starrho)*dLdR/Rsun**2;

EtransTot = trapz(tab_r,abs(dLdR),nlines)

return

end subroutine transgen