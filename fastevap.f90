!!!!! FAST EVAPORATION ROUTINE !!!
!WORKS WELL FOR SIGMA <~ 10 SIGMA_K
!Don't trust in the LTE limit you will underestimate the evaporation rate
!Also needs to be properly generalized to SD scattering (heavier than H)

subroutine fastevap(Nwimps,niso_in,EvapRate)

  use capmod
  ! use akmod
  implicit none
  integer, intent(in):: niso_in
  double precision, intent(in) :: Nwimps
  double precision, intent(out) :: EvapRate
  double precision :: Earg(nlines), muarray(niso_in),nabund(niso_in,nlines),sigma_N(niso_in)
  double precision :: suparg,Knud,rchi
  double precision :: mdmg,mnucg, Tc, rhoc,mn, Tw,nin,escFrac(nlines),vescc(nlines),nxIso(nlines),mfp(nlines),scatrate(nlines)
  double precision, parameter :: kBeV=8.617d-5, GN =6.674d-8, kB=1.3806d-16

  integer :: i, j, k
  mdmg = mdm*1.78d-24
  mnucg = mnuc*1.78d-24

  Tc = tab_T(1)
  rhoc = tab_starrho(1)
  niso = niso_in
  vescc = tab_vesc/c0

  ! print*,"attempting to evaporate, Nwimps = ", Nwimps


  do i = 1,niso_in
  muarray(i) = mdm/AtomicNumber(i)/mnuc
  sigma_N(i) = AtomicNumber(i)**4*(mdm+mnuc)**2/(mdm+AtomicNumber(i)*mnuc)**2 !not yet multiplied by sigma_0
  nabund(i,:) = tab_mfr(:,i)*tab_starrho(:)/AtomicNumber(i)/mnucg
  end do

  call Twimp(nabund,niso_in,Tw)
  Tw = Tw*Tc*kBeV*1.d-9
  ! print*,"Tw = ", Tw
  !this vastly overestimates the evap rate
  ! nxIso(i) = Nwimps*exp(-Rsun**2*tab_r(i)**2/rchi**2)/(pi**(3./2.)*rchi**3)
  nxIso = exp(mdm*vescc**2/2./Tw);
  nin = 4.d0*pi*trapz(tab_r,tab_r**2.*nxIso,nlines) !%niso norm
  nxIso = nxIso/nin

  ! print*,"norm guy ", nin ! "one: ", 4.d0*pi*trapz(tab_r,tab_r**2.*nxIso,nlines)
  !Fraction of the kinetic distribution above the local escape velocity
  escFrac = sqrt(2.d0/pi)*vescc*sqrt(mdm/Tw)*exp(-mdm*vescc**2/Tw/2.d0) - derf(sqrt(mdm/Tw/2.d0)*vescc) + 1.d0;
  ! print*,"EscFrac = ", escFrac,




  !get scattering rate, mean free path
  !this is copied from transgen, maybe unify?
  if ((nq .eq. 0) .and. (nv .eq. 0)) then
    do i = 1,nlines
      mfp(i) = 1./sum(sigma_N*nabund(:,i))/sigma_0/2. !factor of 2 b/c  sigma_tot = 2 sigma_0
    end do
    scatrate = 1./mfp*sqrt(3.*Tw/mdm)
  ! else if ((nq .eq. )) !q, v dependence goes here
  end if

  rchi = (3.*(kB*Tc)/(2.*pi*GN*rhoc*mdmg))**.5;
  Knud = mfp(1)/rchi;
  ! print*,"Knud", Knud
  if (Knud .lt. 0.1) then
    print*, "WARNING, K = ", knud, " is < 0.1. Approximate evaporation scheme is likely very wrong."
  end if

  suparg = trapz(tab_r,rsun/mfp,nlines)
  Earg = escFrac*scatrate*exp(-suparg)*c0

  open(55,file = "sv.dat")
  do j=1,nlines
  write(55,*) tab_r(j), Earg(j),suparg,nxIso(j),escFrac(j),scatrate(j),exp(-suparg)
  end do
  close(55)


  EvapRate = Nwimps*4.*pi*trapz(tab_r,tab_r**2*nxIso*Earg,nlines)

  if (isnan(EvapRate)) then
    stop "NaN evap rate, check it"
  end if

end subroutine fastevap

!Returns WIMP temperature over central temperature
!knows about nv, nq, niso (which controls SI/SD here) via capmod module
subroutine Twimp(nabund,niso_in,Tw)
  use capmod
  implicit none
  integer, intent(in) :: niso_in
  double precision :: sv(nlines),TGeV(nlines), TcGeV,Tw, tol, Tw_out,dT,TwK,mdmg,mN,beta,sigmaN
  double precision :: Tw_out_num(niso), Tw_out_denom(niso),nxIso(nlines)
  double precision mu,nabund(niso_in,nlines)
  double precision, parameter :: GN = 6.674d-8, kB = 1.3806d-16,kBeV=8.617e-5,mnucg=1.67e-24
  integer i,j
  tol = 1.d-8! tolerance: good enough for evap, not for luminosity calc
  TGeV = tab_T*kBeV*1.d-9
  TcGeV = TGeV(1)
  mdmg = mdm*1.78266e-24
  ! print*,"TcGeV ", TcGeV
  Tw = TcGeV
  Tw_out = 0.d0;
  dT = abs(Tw - Tw_out)/Tw;



  do while (dT .gt. Tol)
    TwK = Tw/(kBeV*1.d-9)

    do i=1,Niso
      mN = AtomicNumber(i)*mnuc
      mu = mdm/mN
      beta = AtomicNumber(i)*mnuc*(mdm + mnuc)/mnuc/(mdm + AtomicNumber(i)*mnuc);

      sigmaN = beta**2.*Atomicnumber(i)**2

      call sigmav(2*nv,2*nq,Tw/mdm,TGeV/mn,nlines,sv)
      ! print*,sv
      if (nq .ne. 0) then
          sv = sv*(2.*mdm**2)**(nq)/(1.+mu)**(2.*nq)/q0**(2*nq)
      elseif (nv .ne. 0) then
          sv = sv/v0**(2*nv)
      end if
      nxIso = exp(mdmg*tab_vesc**2/2./TwK/kB)
      Tw_out_num(i) = trapz(tab_r,tab_r**2*TGeV*sv*nxIso*nabund(i,:),nlines);
      Tw_out_denom(i) = trapz(tab_r,tab_r**2*sv*nxIso*nabund(i,:),nlines);



    end do
            Tw_out = sum(Tw_out_num)/sum(Tw_out_denom);
            dT = abs(Tw - Tw_out)/Tw;
            Tw = Tw_out;
  end do
  Tw = Tw/TcGeV;

  ! open(55,file = "svTw.dat")
  ! do j=1,nlines
  ! write(55,*) tab_r(j),sv(j), nxIso(j), nabund(1,j),TGeV(j),tab_vesc(j),TwK
  ! end do
  ! close(55)




  ! call sigmav(2*nv,2*nq,tab_T/mdm,tab_T/mn,nlines,sv)

end subroutine Twimp


SUBROUTINE sigmav(vpow,qpow,xx,xn,nlines,sv) !dimensionless <sigma v>
  implicit none
  double precision, parameter :: pi=3.141592653
  integer, intent(in) :: vpow,qpow,nlines
  double precision :: n, fofn
  double precision, intent(in) :: xx, xn(nlines)
  double precision, intent(out) :: sv(nlines)
! x = T/m;
!remember to multiply by sigma_i/ q^2n
if (vpow .ne. 0) then
    n = dble(vpow/2)
    sv = 2.d0*2.**(n+3./2.)*dgamma(n+2.)*(xx + xn)**(n+1./2.)/sqrt(pi)
elseif (qpow .ne. 0) then
    n = qpow/2
    sv = 2.d0**(n+3./2.)*dgamma(n+2.)*(xx + xn)**(n+1./2.)/sqrt(pi)*fofn(n) !needs to be multiplied by that other factor
else
    sv = 2.*2.d0**(3./2.)*sqrt(xx+xn)/sqrt(pi)
end if
end subroutine sigmav

double precision function fofn(n)
  double precision f,n
if (n .le. 0.d0) then
    f = 2.d0
else
    f = 2.d0**((n)+1.)/((n)+1.)
end if
fofn = f
return
end function

!now in Capmod
! !Fast trapezoidal integral
!   function trapz(x,y,flen)
!   implicit none
!   integer, intent(in) :: flen
!   double precision, intent (in) :: x(flen), y(flen)
!   double precision trapz
!
!   integer i
!
!
!   trapz = y(1)*(x(2)-x(1))/2. + y(flen)*(x(flen)-x(flen-1))/2.
!   do i = 2,flen-1
!     trapz = trapz + y(i)*(x(i)-x(i-1))
!
!     if (trapz .lt. 0.d0) then
!       print*, "negative encountered in trapz: i = ", i
!     end if
!   end do
!
!
!   return
!   end function
