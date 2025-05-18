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

    module capture_mod

      use shared_mod
      implicit none
      double precision, parameter :: GNewt = 6.672d-8
      double precision, parameter :: q0 = 0.04,v0 = 220.d5

      ! nq and nv can be -1, 0, 1, 2; this is set in the main program
      integer :: nq, nv

        contains

      !generalized form factor: hydrogen
      function gffi_h_qv(vel_dm, vel_esc)
      double precision :: p, vel_dm,vel_esc,u,gffi_h_qv,G
      p = m_dm*vel_dm
      u = sqrt(vel_dm**2-vel_esc**2)
      if (nq .ne. -1) then
        G = (p/q0/c0)**(2.d0*dble(nq))*m_dm*vel_dm**2/(2.d0*mu**dble(nq))*1./(1.+dble(nq)) &
        *((mu/mu_plus**2)**(dble(nq)+1.)-(u**2/vel_dm**2)**(dble(nq)+1.))
      else
        G = ((p)/q0/c0)**(2.d0*dble(nq))*m_dm*vel_dm**2/(2.d0*mu**dble(nq))*log(mu/mu_plus**2*vel_dm**2/(u)**2)
      endif
      gffi_h_qv = G
      end function gffi_h_qv

      !generalized form factor: other elements
      function gffi_a_qv(vel_dm, vel_esc, atomic_number)
        double precision :: p, vel_dm,vel_esc,u,mN,atomic_number,Ei,B
        double precision :: dgamic,gffi_a_qv
        p = m_dm*vel_dm
        u = sqrt(vel_dm**2-vel_esc**2)
        mN = atomic_number*m_proton
        Ei  = 5.8407d-2/(mN*(0.91*mN**(1./3.)+0.3)**2)
        B = .5*m_dm*vel_dm**2/Ei/c0**2
        if (nq .eq. 0) then
          gffi_a_qv = Ei*c0**2*(exp(-m_dm*u**2/2/Ei/c0**2)-exp(-B*mu/mu_plus**2))
        else
          gffi_a_qv = ((p)/q0/c0)**(2*dble(nq))*Ei*c0**2/(B*mu)**dble(nq)*(dgamic(1.+dble(nq),B*u**2/vel_dm**2) &
                  - dgamic(1.+dble(nq),B*mu/mu_plus**2))
        end if
      end function gffi_a_qv

      !Fast trapezoidal integral
      function trapezoid(x, y, length)
        implicit none
        integer, intent(in) :: length
        double precision, intent (in) :: x(length), y(length)
        double precision trapezoid
        integer i

        trapezoid = y(1)*(x(2)-x(1))/2. 

        do i = 2,length-1
          trapezoid = trapezoid + y(i)*(x(i)-x(i-1))
        end do

        trapezoid = trapezoid + y(length)*(x(length)-x(length-1))/2.

        return
      end function

    end module capture_mod


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    ! Some functions that have to be external, because of the integrator.


    !The integrand for the integral over u
    function velocity_integrand_qv(init_velocity, dist_over_vel)
      use capture_mod
      double precision :: init_velocity, w, velocity_integrand_qv, dist_over_vel
      external dist_over_vel

      w = sqrt(init_velocity**2+vesc_shared**2)

      !Switch depending on whether we are capturing on Hydrogen or not
      if (atomic_shared .gt. 2.d0) then
        velocity_integrand_qv = dist_over_vel(init_velocity)*gffi_a_qv(w,vesc_shared,atomic_shared)
      else
        velocity_integrand_qv = dist_over_vel(init_velocity)*gffi_h_qv(w,vesc_shared)
      end if

      !Rescale for velocity-dependent cross-sections
      if (nv .ne. 0) then
        velocity_integrand_qv = velocity_integrand_qv*(w/v0)**(2*nv)
      end if
    end function velocity_integrand_qv


    subroutine capture_rate_qv(mass_dm, sigma_0, num_isotopes, q_pow, v_pow, is_spin_dep, capture_rate)
      use capture_mod
      implicit none
      integer, intent(in):: q_pow, v_pow, num_isotopes, is_spin_dep
      ! integer, intent(in):: is_spin_dep
      integer eli, ri, limit
      double precision, intent(in) :: mass_dm, sigma_0
      double precision :: capture_rate !this is the output
      double precision :: sigma_SD, sigma_SI
      double precision :: capture_maximum, maxcapped, a, muminus, sigma_N, umax, umin, vesc
      double precision :: epsabs, epsrel, abserr, neval  !for integrator
      double precision :: ier,alist,blist,rlist,elist,iord,last!for integrator
      double precision :: int_result

      dimension alist(1000),blist(1000),elist(1000),iord(1000),   rlist(1000)!for integrator
      external velocity_integrand_qv
      external gaussian_test !this is just for testing

      epsabs=1.d-8
      epsrel=1.d-8
      limit=1000

      m_dm = mass_dm
      nq = q_pow
      nv = v_pow

      if (is_spin_dep == 1) then
        sigma_SD = sigma_0
        sigma_SI = 0.d0
      else if (is_spin_dep == 0) then
        sigma_SD = 0.d0
        sigma_SI = sigma_0
      end if

      if (nq*nv .ne. 0) then
        stop "Oh no! nq and nv can't both be nonzero."
      end if

      if (.not. allocated(star_r)) then
        stop "You haven't yet called init_sun to load the solar model!"
      end if

      capture_rate = 0.d0

      !Loop over the shells of constant radius in the star
      do ri = 1, nlines

        vesc = star_escape(ri)
        vesc_shared = vesc !make accessible via the module

        !Loop over the different elements
        do eli = 1, num_isotopes

          a = atomic_nums(eli)
          atomic_shared = a !make accessible via the module

          !This is fine for SD as long as it's just hydrogen. Otherwise, spins must be added.
          sigma_N = a**2 * (sigma_SI*a**2 + sigma_SD) * (mass_dm+m_proton)**2/(mass_dm+a*m_proton)**2

          mu = mass_dm/(m_proton*a)
          mu_plus = (1.+mu)/2.
          muminus = (mu-1.d0)/2.

          ! Bottom part of the integral is always zero -- happy little slow DM particles can always be captured.
          umin = 0.d0
          ! Chop the top of the integral off at the smaller of the halo escape velocity or the minimum velocity required for capture.
          umax = min(vesc * sqrt(mu)/abs(muminus), escape_halo)

          !Call integrator
          call dsntdqagse(velocity_integrand_qv,distribution_over_vel,umin,umax, &
          epsabs,epsrel,limit,int_result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
          int_result = int_result * 2.d0 * sigma_N * avogadro * star_rho(ri)*star_fractions(ri,eli) * (mu_plus/mass_dm)**2
          capture_rate = capture_rate + star_r(ri)**2*int_result*star_dr(ri)

          if (isnan(capture_rate)) then
            capture_rate = 0.d0
            stop 'NaN encountered whilst trying compute capture rate.'
          end if

        end do

      end do

      capture_rate = 4.d0*pi*radius_star**3*capture_rate

      if (capture_rate .gt. 1.d100) then
        print*,"Capt'n General says: Oh my, it looks like you are capturing an"
        print*,"infinite amount of dark matter in the Sun. Best to look into that."
      end if

      maxcapped = capture_maximum(mass_dm)
      if (capture_rate .gt. maxcapped) then
        capture_rate = maxcapped
      end if
    end subroutine capture_rate_qv


    !! capture_rate_constant calculates the capture rate for constant cross section.
    ! subroutine capture_rate_constant(mx_in,sigma_0,capped_SD,capped_SI)
    !   implicit none
    !   double precision, intent(in) :: mx_in, sigma_0
    !   double precision :: capped_SD,capped_SI

    !   call capture_rate_qv(mx_in,sigma_0,1,0,0,1,capped_SD)
    !   call capture_rate_qv(mx_in,sigma_0,29,0,0,0,capped_SI)
    ! end subroutine capture_rate_constant

    subroutine capture_rate_constant(mx_in,sigma_0_SD_in,sigma_0_SI_in,capped_SD,capped_SI)
      implicit none
      double precision, intent(in) :: mx_in, sigma_0_SD_in,sigma_0_SI_in
      double precision :: capped_SD,capped_SI

      call capture_rate_qv(mx_in,sigma_0_SD_in,1,0,0,1,capped_SD)
      call capture_rate_qv(mx_in,sigma_0_SI_in,29,0,0,0,capped_SI)
    end subroutine capture_rate_constant


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! For mesa interface only: allocate arrays.
  subroutine allocate_stellar_arrays(nlines_mesa)
    use capture_mod
    integer, intent(in) :: nlines_mesa
    nlines = nlines_mesa
    allocate(star_enclosed(nlines))       !M(<r)
    allocate(star_r(nlines))           !r
    allocate(star_rho(nlines))     !rho
    allocate(star_fractions(nlines,8))       !mass fraction per isotope
    allocate(tab_atomic(8))
    allocate(star_escape(nlines))        !local escape velocity
    allocate(star_temp(nlines))           !temperature
    ! allocate(phi(nlines)) !! <--- not needed; computed in wimp_support.f
    allocate(star_dr(nlines))          !dr (nice)
    allocate(star_grav(nlines))           !local gravitational acceleration, needed for transport

    RETURN
  end subroutine allocate_stellar_arrays

  subroutine deallocate_stellar_arrays()
    use capture_mod
    deallocate(star_enclosed)
    deallocate(star_r)
    deallocate(star_rho)
    deallocate(star_fractions) !we could just allocate niso, but this leads to problems
    deallocate(tab_atomic)
    deallocate(star_escape)
    deallocate(star_temp)
    deallocate(star_dr)
    deallocate(star_grav)
    RETURN
  end subroutine deallocate_stellar_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This is called INSTEAD of read_solar_params, for use with MESA interface.
  subroutine get_stellar_params(rmesa,rhomesa,mfrmesa,atomicmesa,mesavesc,Tmesa, &
                                mesag,mesamass,mesaradius,rho0_in,usun_in,u0_in,vesc_in)
    use capture_mod
    !mesamass & mesaradius unused here but subroutine used in a few other places so I left them
    !in just in case
    double precision :: mesamass, mesaradius
    double precision :: rhomesa(nlines), rmesa(nlines), mfrmesa(8,nlines)
    double precision :: mesavesc(nlines),mesag(nlines),Tmesa(nlines)
    double precision :: atomicmesa(8)
    integer i
    double precision,intent(in) :: rho0_in,usun_in,u0_in,vesc_in

    vel_sun = usun_in*1.d5
    dispersion_dm =  u0_in*1.d5
    rho_dm =rho0_in
    escape_halo = vesc_in*1.d5

    radius_star = rmesa(nlines)
    star_r = rmesa/radius_star
    star_rho = rhomesa
    star_escape = mesavesc
    star_temp = tmesa
    star_grav = -mesag
    do i= 1,8
      star_fractions(:,i) = mfrmesa(i,:)
    end do
    tab_atomic = atomicmesa
    atomic_nums(1:8) = tab_atomic

    do i = 1, nlines-1
      star_dr(i) = -star_r(i)+star_r(i+1) !while we're here, populate dr
    end do
    star_dr(nlines) = star_r(nlines)-star_r(nlines-1)

    RETURN
  end subroutine get_stellar_params


  subroutine getnlines(nlines_out) !a little auxiliary trick
    use capture_mod
    integer, intent(out) :: nlines_out
    nlines_out = nlines
    return
  end subroutine getnlines
