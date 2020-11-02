!   Capt'n General
!   For a continuum of Q-dependent capture
!   Simplified, general solar DM capture routine
!   Standalone code for q^2n, v^2n
!   Useful stuff is run at the end; beginning is the module that does the heavy lifting
!   Future plans: add form factor handling (a la Catena & Schwabe)
!   Made for GAMBIT, with marginal competence
!   Aaron Vincent 2017
!   all units of distance: cm
!   all units of mass/energy : GeV (or GeV/c^2, don't forget)
!   all units of time: seconds
!   Sticking with notation of 1504.04378. Cite that paper. Or 1605.06502 it's even better.
!   Reference q0 is 40 MeV, and v0 is 220 km/s.

!   Updated 2020 just to handle all nq,nv=-1,0,1,2 cases (integration limits were causing issues).
!   Working for spin-dependent interactions with atomic hydrogen (niso=1)
!   NOTE: removed evaporation calcs - for some reason fastevap was still being used even when
!     option was turned off from DarkMESA side.

    module capmod
      implicit none
      double precision, parameter :: pi=3.141592653, NAvo=6.0221409d23,GMoverR = 1.908e15,GNewt = 6.672d-8
      double precision, parameter :: c0 = 2.99792458d10, mnuc = 0.938, q0 = 0.04,v0 = 220.d5
      !these are now set in captn_init
      double precision :: usun , u0 ,rho0, vesc_halo, Rsun
      !this goes with the Serenelli table format
      double precision :: AtomicNumber(29) !29 is is the number from the Serenelli files; if you have fewer it shouldn't matter
      double precision, allocatable :: tab_mencl(:),tab_starrho(:),tab_mfr(:,:), tab_atomic(:)
      double precision, allocatable :: tab_r(:), tab_vesc(:), tab_dr(:),tab_T(:),tab_g(:)

      ! nq and nv can be -1, 0, 1, 2; this is set in the main program
      integer :: nq, nv, nlines
      double precision :: mdm, vesc_shared, a_shared, mu, muplus

        contains

      !velocity distribution,
      function vdist_over_u(u)
      double precision :: u, vdist_over_u, normfact
      vdist_over_u = (3./2.)**(3./2.)*4.*rho0*u/sqrt(pi)/mdm/u0**3 &
      *exp(-3.*(usun**2+u**2)/(2.*u0**2))*sinh(3.*u*usun/u0**2)/(3.*u*usun/u0**2)
      !normfact = .5*erf(sqrt(3./2.)*(vesc_halo-usun)/u0) + &
      !.5*erf(sqrt(3./2.)*(vesc_halo+usun)/u0)+ u0/(sqrt(6.*pi)*usun) &
      !*(exp(-3.*(usun+vesc_halo)/2./u0**2)-exp(-3.*(usun-vesc_halo)/2./u0**2))
      normfact = 1.
      !print*,normfact
      vdist_over_u = vdist_over_u/normfact
      end function vdist_over_u

      !generalized form factor: hydrogen
      function GFFI_H(w,vesc)
      double precision :: p, w,vesc,u,GFFI_H,G
      p = mdm*w
      u = sqrt(w**2-vesc**2)
      if (nq .ne. -1) then
        G = (p/q0/c0)**(2.d0*dble(nq))*mdm*w**2/(2.d0*mu**dble(nq))*1./(1.+dble(nq)) &
        *((mu/muplus**2)**(dble(nq)+1.)-(u**2/w**2)**(dble(nq)+1.))
      else
        !eps added to make the log finite: the value of eps does not affect the answer
        G = ((p)/q0/c0)**(2.d0*dble(nq))*mdm*w**2/(2.d0*mu**dble(nq))*log(mu/muplus**2*w**2/(u+eps)**2)
      endif
      GFFI_H = G
      end function GFFI_H

      !generalized form factor: other elements
      function GFFI_A(w,vesc,A)
        double precision :: p, w,vesc,u,mN,A,Ei,B
        double precision :: dgamic,GFFI_A
        p = mdm*w
        u = sqrt(w**2-vesc**2)
        mN = A*mnuc
        Ei  = 5.8407d-2/(mN*(0.91*mN**(1./3.)+0.3)**2)
        B = .5*mdm*w**2/Ei/c0**2
        if (nq .eq. 0) then
          GFFI_A = Ei*c0**2*(exp(-mdm*u**2/2/Ei/c0**2)-exp(-B*mu/muplus**2))
        else
          GFFI_A = ((p+eps)/q0/c0)**(2*dble(nq))*Ei*c0**2/(B*mu)**dble(nq)*(dgamic(1.+dble(nq),B*u**2/w**2+eps) &
                  - dgamic(1.+dble(nq),B*mu/muplus**2+eps))
        end if
      end function GFFI_A

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !read in solar parameters from Aldo Serenelli-style files, with header removed
      subroutine get_solar_params(filename,nlines)
        character*300 :: filename
        double precision :: Pres, Lumi !these aren't used, but dummies are required
        double precision, allocatable :: phi(:) !this is used briefly
        integer :: i,j, nlines,iostatus

        Rsun = 69.57d9 !this is set here, for other stars, this sub is not called

        !Get number of lines in the file
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

        !allocate the arrays
        allocate(tab_mencl(nlines))
        allocate(tab_r(nlines))
        allocate(tab_starrho(nlines))
        allocate(tab_mfr(nlines,29)) !we could just allocate niso, but this leads to problems
        allocate(tab_vesc(nlines))
        allocate(phi(nlines))
        allocate(tab_dr(nlines))
        allocate(tab_T(nlines)) !not used in capgen; used for transgen (and anngen? )
        allocate(tab_g(nlines))


        !now actually read in the file
        open(99,file=filename)
        do i=1,nlines
          read(99,*) tab_mencl(i),tab_r(i), tab_T(i), tab_starrho(i), Pres, Lumi, tab_mfr(i,:)
        end do
        close(99)

        !we calculate the escape velocity here since all the ingredients are ready
        phi(nlines) = -GMoverR
        tab_vesc(nlines) = sqrt(-2.d0*phi(nlines))
        tab_dr(nlines) = tab_r(nlines)-tab_r(nlines-1)
        do i = 1,nlines-1
          j = nlines-i !trapezoid integral
          phi(j) = phi(j+1) + GMoverR*(tab_r(j)-tab_r(j+1))/2.*(tab_mencl(j)/tab_r(j)**2+tab_mencl(j+1)/tab_r(j+1)**2)
          tab_vesc(j) = sqrt(-2.d0*phi(j)) !escape velocity in cm/s
          tab_dr(j) = -tab_r(j)+tab_r(j+1) !while we're here, populate dr
          ! tab_g(j) = -(-phi(j)+phi(j+1))/tab_dr(j)
          tab_g(i) = -GMoverR*tab_mencl(i)/tab_r(i)**2/Rsun
        end do
        ! tab_g(nlines) = tab_g(nlines-1)
        tab_g(nlines) = -GMoverR*tab_mencl(nlines)/tab_r(nlines)**2/Rsun

        open(55,file = "tab_serenelli.dat")
        do i=1,nlines
          write(55,*) tab_r(i), tab_starrho(i), tab_vesc(i), tab_mfr(i,:)
          end do
          close(55)

          ! Populate the atomic number tables here (because it relies on a specific format)
        AtomicNumber  = (/ 1., 4., 3., 12., 13., 14., 15., 16., 17., &
                          18., 20.2, 22.99, 24.3, 26.97, 28.1, 30.97,32.06, 35.45, &
                          39.948, 39.098, 40.08, 44.95, 47.86, 50.94, 51.99, &
                          54.93, 55.845, 58.933, 58.693/)


        return
      end subroutine get_solar_params

      !this is to make sure the integrator does what it's supposed to
      function gaussinmod(x)
        double precision :: x,gaussinmod
        gaussinmod = nq*exp(-x**2/2.d0)
      end function gaussinmod


      !Fast trapezoidal integral
      function trapz(x,y,flen)
        implicit none
        integer, intent(in) :: flen
        double precision, intent (in) :: x(flen), y(flen)
        double precision trapz
        integer i

        trapz = y(1)*(x(2)-x(1))/2. + y(flen)*(x(flen)-x(flen-1))/2.
        do i = 2,flen-1
          trapz = trapz + y(i)*(x(i)-x(i-1))

          if (trapz .lt. 0.d0) then
            print*, "negative encountered in trapz: i = ", i
          end if
        end do

        return
      end function

    end module capmod


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    ! Some functions that have to be external, because of the integrator.

    !Just a test for the integrator. Nothing to see here
    function gausstest(x)
      use capmod
      double precision :: x,gausstest
      gausstest = gaussinmod(x)
    end function gausstest


    !The integrand for the integral over u
    function integrand(u,foveru)
      use capmod
      double precision :: u, w, integrand, foveru
      external foveru

      w = sqrt(u**2+vesc_shared**2)

      !Switch depending on whether we are capturing on Hydrogen or not
      if (a_shared .gt. 2.d0) then
        integrand = foveru(u)*GFFI_A(w,vesc_shared,a_shared)
      else
        integrand = foveru(u)*GFFI_H(w,vesc_shared)
      end if

      !Rescale for velocity-dependent cross-sections
      if (nv .ne. 0) then
        integrand = integrand*(w/v0)**(2*nv)
      end if
    end function integrand


    ! This is for use with MESA arrays specifically (though could likely be used for GARSTEC)
    ! star initialization done elsewhere
    ! subroutine captn_mesa(mx_in,sigma_0_in,niso_in,nq_in,nv_in,capped)
    !   subroutine captn_mesa()
    !   use capmod
    !   implicit none
    !   ! integer, intent(in):: nq_in, niso_in, nv_in
    !   integer i, ri
    !   ! double precision, intent(in) :: mx_in, sigma_0_in
    !   double precision :: capped, maxcap !this is the output
    !   double precision :: epsabs, epsrel,limit,result,abserr,neval !for integrator
    !   double precision :: ier,alist,blist,rlist,elist,iord,last!for integrator
    !   double precision, allocatable :: u_int_res(:)
    !
    !
    ! end subroutine captn_mesa


    subroutine captn_general(mx_in,sigma_0,niso,nq_in,nv_in,spin_in,capped)
      use capmod
      implicit none
      integer, intent(in):: nq_in, nv_in, niso, spin_in
      ! integer, intent(in):: spin_in
      integer eli, ri, limit
      double precision, intent(in) :: mx_in, sigma_0
      double precision :: capped !this is the output
      double precision :: sigma_SD, sigma_SI
      double precision :: maxcap, maxcapped, a, muminus, sigma_N, umax, umin, vesc
      double precision :: epsabs, epsrel, abserr, neval  !for integrator
      double precision :: ier,alist,blist,rlist,elist,iord,last!for integrator
      double precision, allocatable :: u_int_res(:)

      dimension alist(1000),blist(1000),elist(1000),iord(1000),   rlist(1000)!for integrator
      external integrand
      external gausstest !this is just for testing

      epsabs=1.d-8
      epsrel=1.d-8
      limit=1000

      mdm = mx_in
      nq = nq_in
      nv = nv_in

      if (spin_in == 1) then
        sigma_SD = sigma_0
        sigma_SI = 0.d0
      else if (spin_in == 0) then
        sigma_SD = 0.d0
        sigma_SI = sigma_0
      end if

      if (nq*nv .ne. 0) then
      print*, "Oh no! nq and nv can't both be nonzero."
      return
      end if

      if (.not. allocated(tab_r)) then
        stop "You haven't yet called captn_init to load the solar model!"
      end if
      allocate(u_int_res(nlines))

      capped = 0.d0

      !Loop over the shells of constant radius in the star
      do ri = 1, nlines

        vesc = tab_vesc(ri)
        vesc_shared = vesc !make accessible via the module

        !Loop over the different elements
        do eli = 1, niso

          a = AtomicNumber(eli)
          a_shared = a !make accessible via the module

          !This is fine for SD as long as it's just hydrogen. Otherwise, spins must be added.
          sigma_N = a**2 * (sigma_SI*a**2 + sigma_SD) * (mx_in+mnuc)**2/(mx_in+a*mnuc)**2

          mu = mx_in/(mnuc*a)
          muplus = (1.+mu)/2.
          muminus = (mu-1.d0)/2.

          ! Bottom part of the integral is always zero -- happy little slow DM particles can always be captured.
          umin = 0.d0
          ! Chop the top of the integral off at the smaller of the halo escape velocity or the minimum velocity required for capture.
          umax = min(vesc * sqrt(mu)/abs(muminus), vesc_halo)

          !Call integrator
          call dsntdqagse(integrand,vdist_over_u,umin,umax, &
          epsabs,epsrel,limit,u_int_res(ri),abserr,neval,ier,alist,blist,rlist,elist,iord,last)
          u_int_res(ri) = u_int_res(ri) * 2.d0 * sigma_N * NAvo * tab_starrho(ri)*tab_mfr(ri,eli) * (muplus/mx_in)**2
          capped = capped + tab_r(ri)**2*u_int_res(ri)*tab_dr(ri)

          if (isnan(capped)) then
            capped = 0.d0
            stop 'NaN encountered whilst trying compute capture rate.'
          end if

        end do

      end do

      capped = 4.d0*pi*Rsun**3*capped

      if (capped .gt. 1.d100) then
        print*,"Capt'n General says: Oh my, it looks like you are capturing an"
        print*,"infinite amount of dark matter in the Sun. Best to look into that."
      end if

      maxcapped = maxcap(mx_in)
      if (capped .gt. maxcapped) then
        capped = maxcapped
      end if
    end subroutine captn_general


    !This is fine as long as the escape velocity is large enough
    function maxcap(mx)
      use capmod
      implicit none
      double precision maxcap
      double precision, intent(in) :: mx

      maxcap = pi/3.d0*rho0/mx*Rsun**2 &
      *(exp(-3./2.*usun**2/u0**2)*sqrt(6.d0/pi)*u0 &
      + (6.d0*GMoverR/usun + (u0**2 + 3.d0*usun**2)/usun)*erf(sqrt(3./2.)*usun/u0))

    end function maxcap


    !! captn_specific calculates the capture rate for constant cross section.
    ! subroutine captn_specific(mx_in,sigma_0,capped_SD,capped_SI)
    !   implicit none
    !   double precision, intent(in) :: mx_in, sigma_0
    !   double precision :: capped_SD,capped_SI

    !   call captn_general(mx_in,sigma_0,1,0,0,1,capped_SD)
    !   call captn_general(mx_in,sigma_0,29,0,0,0,capped_SI)
    ! end subroutine captn_specific

    subroutine captn_specific(mx_in,sigma_0_SD_in,sigma_0_SI_in,capped_SD,capped_SI)
      implicit none
      double precision, intent(in) :: mx_in, sigma_0_SD_in,sigma_0_SI_in
      double precision :: capped_SD,capped_SI

      call captn_general(mx_in,sigma_0_SD_in,1,0,0,1,capped_SD)
      call captn_general(mx_in,sigma_0_SI_in,29,0,0,0,capped_SI)
    end subroutine captn_specific


!------!------!------!------!------INITIALIZATION FCT

    subroutine captn_init(solarmodel,rho0_in,usun_in,u0_in,vesc_in)
      !input velocities in km/s, not cm/s!!!
      use capmod
      use iso_c_binding, only: c_ptr
      implicit none
      character (len=300) solarmodel
      double precision,intent(in) :: rho0_in,usun_in,u0_in,vesc_in
      !common solarmodel
      !external solarmodel

      if  (.not. allocated(tab_r)) then !
          print*,"Capgen initializing from model: ",solarmodel
          call get_solar_params(solarmodel,nlines)
      end if

      usun = usun_in*1.d5
      u0 =  u0_in*1.d5
      rho0 =rho0_in
      vesc_halo = vesc_in*1.d5

    end subroutine captn_init



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! For mesa interface only: allocate arrays.
  subroutine allocate_stellar_arrays(nlines_mesa)
    use capmod
    integer, intent(in) :: nlines_mesa
    nlines = nlines_mesa
    allocate(tab_mencl(nlines))       !M(<r)
    allocate(tab_r(nlines))           !r
    allocate(tab_starrho(nlines))     !rho
    allocate(tab_mfr(nlines,8))       !mass fraction per isotope
    allocate(tab_atomic(8))
    allocate(tab_vesc(nlines))        !local escape velocity
    allocate(tab_T(nlines))           !temperature
    ! allocate(phi(nlines)) !! <--- not needed; computed in wimp_support.f
    allocate(tab_dr(nlines))          !dr (nice)
    allocate(tab_g(nlines))           !local gravitational acceleration, needed for transport

    RETURN
  end subroutine allocate_stellar_arrays

  subroutine deallocate_stellar_arrays()
    use capmod
    deallocate(tab_mencl)
    deallocate(tab_r)
    deallocate(tab_starrho)
    deallocate(tab_mfr) !we could just allocate niso, but this leads to problems
    deallocate(tab_atomic)
    deallocate(tab_vesc)
    deallocate(tab_T)
    deallocate(tab_dr)
    deallocate(tab_g)
    RETURN
  end subroutine deallocate_stellar_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This is called INSTEAD of get_solar_params, for use with MESA interface.
  subroutine get_stellar_params(rmesa,rhomesa,mfrmesa,atomicmesa,mesavesc,Tmesa, &
                                mesag,mesamass,mesaradius,rho0_in,usun_in,u0_in,vesc_in)
    use capmod
    !mesamass & mesaradius unused here but subroutine used in a few other places so I left them
    !in just in case
    double precision :: mesamass, mesaradius
    double precision :: rhomesa(nlines), rmesa(nlines), mfrmesa(8,nlines)
    double precision :: mesavesc(nlines),mesag(nlines),Tmesa(nlines)
    double precision :: atomicmesa(8)
    integer i
    double precision,intent(in) :: rho0_in,usun_in,u0_in,vesc_in

    usun = usun_in*1.d5
    u0 =  u0_in*1.d5
    rho0 =rho0_in
    vesc_halo = vesc_in*1.d5

    Rsun = rmesa(nlines)
    tab_r = rmesa/Rsun
    tab_starrho = rhomesa
    tab_vesc = mesavesc
    tab_T = tmesa
    tab_g = -mesag
    do i= 1,8
      tab_mfr(:,i) = mfrmesa(i,:)
    end do
    tab_atomic = atomicmesa

    do i = 1, nlines-1
      tab_dr(i) = -tab_r(i)+tab_r(i+1) !while we're here, populate dr
    end do
    tab_dr(nlines) = tab_r(nlines)-tab_r(nlines-1)

    AtomicNumber(1:8) = tab_atomic

    RETURN
  end subroutine get_stellar_params


  subroutine getnlines(nlines_out) !a little auxiliary trick
    use capmod
    integer, intent(out) :: nlines_out
    nlines_out = nlines
    return
  end subroutine getnlines
