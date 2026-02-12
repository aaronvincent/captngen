!   Capt'n General subroutines
!   Split from gencap.f90 so that gfortran properly compiles OMP directives
!   (gfortran silently drops OMP in subroutines defined after 'end module' in the same file)

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
      double precision :: int_result

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
        stop "Oh no! nq and nv can't both be nonzero."
      end if

      if (.not. allocated(tab_r)) then
        stop "You haven't yet called captn_init to load the solar model!"
      end if

      capped = 0.d0

      !Loop over the shells of constant radius in the star
      !$OMP parallel &
      !$OMP private(ri, vesc, eli, a, sigma_N, muminus, umin, umax, int_result, &
      !$OMP   abserr, neval, ier, alist, blist, rlist, elist, iord, last) &
      !$OMP shared(nlines, niso, mx_in, sigma_SI, sigma_SD, tab_vesc, tab_starrho, &
      !$OMP   tab_mfr, tab_r, tab_dr, AtomicNumber, mdm, vesc_halo, epsabs, epsrel, limit, nq, nv) &
      !$OMP reduction(+:capped)
      !$OMP do schedule(dynamic)
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
          epsabs,epsrel,limit,int_result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
          int_result = int_result * 2.d0 * sigma_N * NAvo * tab_starrho(ri)*tab_mfr(ri,eli) * (muplus/mx_in)**2
          capped = capped + tab_r(ri)**2*int_result*tab_dr(ri)

        end do

      end do
      !$OMP end do
      !$OMP end parallel

      if (isnan(capped)) then
        capped = 0.d0
        stop 'NaN encountered whilst trying compute capture rate.'
      end if

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


    subroutine captn_specific(mx_in,sigma_0_SD_in,sigma_0_SI_in,capped_SD,capped_SI)
      implicit none
      double precision, intent(in) :: mx_in, sigma_0_SD_in,sigma_0_SI_in
      double precision :: capped_SD,capped_SI

      call captn_general(mx_in,sigma_0_SD_in,1,0,0,1,capped_SD)
      call captn_general(mx_in,sigma_0_SI_in,29,0,0,0,capped_SI)
    end subroutine captn_specific


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
    AtomicNumber(1:8) = tab_atomic

    do i = 1, nlines-1
      tab_dr(i) = -tab_r(i)+tab_r(i+1) !while we're here, populate dr
    end do
    tab_dr(nlines) = tab_r(nlines)-tab_r(nlines-1)

    RETURN
  end subroutine get_stellar_params


  subroutine getnlines(nlines_out) !a little auxiliary trick
    use capmod
    integer, intent(out) :: nlines_out
    nlines_out = nlines
    return
  end subroutine getnlines
