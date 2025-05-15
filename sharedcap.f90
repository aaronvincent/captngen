!   Capt'n Shared
!   Designed as a module to house shared vaiables and functions
!   between both the General and Operator varients
!   Most of this was already written by Aaron Vincent in the older gencap.f90 file
!   Neal Avis Kozar 2020
!   all units of distance: cm
!   all units of mass/energy : GeV (or GeV/c^2, don't forget)
!   all units of time: seconds
!   Sticking with notation of 1504.04378. Cite that paper. Or 1605.06502 it's even better.


module shared_mod
    use omp_lib
    implicit none
    double precision, parameter :: pi=3.141592653, avogadro=6.0221409d23, gm_over_r_sun=1.908e15
    double precision, parameter :: c0=2.99792458d10, m_proton=0.938
    !these are now set in init_sun
    double precision :: vel_sun, dispersion_dm, rho_dm, escape_halo, radius_star
    !tab: means tabulated from file; so as not to be confused with other variables
    double precision, allocatable :: tab_mencl(:), tab_starrho(:), tab_mfr(:,:), tab_r(:), tab_vesc(:), tab_dr(:)
    double precision, allocatable :: tab_mfr_oper(:,:), tab_T(:), tab_g(:), tab_atomic(:), vesc_shared_arr(:)
    !this goes with the Serenelli table format
    double precision :: atomic_nums(29) !29 is is the number from the Serenelli files; if you have fewer it shouldn't matter

    integer :: nlines, shell_index_shared!, ri_for_omega
    double precision :: mdm, vesc_shared, atomic_shared, mu, mu_plus
    !$OMP threadprivate(shell_index_shared, atomic_shared)
    
    contains

    !   this is the function f_sun(u) in 1504.04378 eqn 2.2 divided by u
    !velocity distribution,
    function vdist_over_u(u)
        double precision :: u, vdist_over_u, normfact
        vdist_over_u = (3./2.)**(3./2.)*4.*rho_dm*u/sqrt(pi)/mdm/dispersion_dm**3 &
        *exp(-3.*(vel_sun**2+u**2)/(2.*dispersion_dm**2))*sinh(3.*u*vel_sun/dispersion_dm**2)/(3.*u*vel_sun/dispersion_dm**2)
        !normfact = .5*erf(sqrt(3./2.)*(escape_halo-vel_sun)/dispersion_dm) + &
        !.5*erf(sqrt(3./2.)*(escape_halo+vel_sun)/dispersion_dm)+ dispersion_dm/(sqrt(6.*pi)*vel_sun) &
        !*(exp(-3.*(vel_sun+escape_halo)/2./dispersion_dm**2)-exp(-3.*(vel_sun-escape_halo)/2./dispersion_dm**2))
        normfact = 1.
        !print*,normfact
        vdist_over_u = vdist_over_u/normfact
    end function vdist_over_u

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !read in solar parameters from Aldo Serenelli-style files, with header removed
    subroutine get_solar_params(filename,nlines)
        character*300 :: filename
        double precision :: Pres, Lumi !these aren't used, but dummies are required
        double precision, allocatable :: phi(:) !this is used briefly
        integer :: i,j, nlines,iostatus

        radius_star = 69.57d9 !this is set here, for other stars, this sub is not called

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
        allocate(tab_mfr_oper(nlines,16)) ! for the operator method
        allocate(vesc_shared_arr(nlines)) ! for OMP stuff


        !now actually read in the file
        open(99,file=filename)
        do i=1,nlines
          read(99,*) tab_mencl(i),tab_r(i), tab_T(i), tab_starrho(i), Pres, Lumi, tab_mfr(i,:)
        end do
        close(99)

        !we calculate the escape velocity here since all the ingredients are ready
        phi(nlines) = -gm_over_r_sun
        tab_vesc(nlines) = sqrt(-2.d0*phi(nlines))
        tab_dr(nlines) = tab_r(nlines)-tab_r(nlines-1)
        do i = 1,nlines-1
          j = nlines-i !trapezoid integral
          phi(j) = phi(j+1) + gm_over_r_sun*(tab_r(j)-tab_r(j+1))/2.*(tab_mencl(j)/tab_r(j)**2+tab_mencl(j+1)/tab_r(j+1)**2)
          tab_vesc(j) = sqrt(-2.d0*phi(j)) !escape velocity in cm/s
          tab_dr(j) = -tab_r(j)+tab_r(j+1) !while we're here, populate dr
          ! tab_g(j) = -(-phi(j)+phi(j+1))/tab_dr(j)
          tab_g(i) = -gm_over_r_sun*tab_mencl(i)/tab_r(i)**2/radius_star
        end do
        ! tab_g(nlines) = tab_g(nlines-1)
        tab_g(nlines) = -gm_over_r_sun*tab_mencl(nlines)/tab_r(nlines)**2/radius_star

          ! Populate the atomic number tables here (because it relies on a specific format)
        atomic_nums  = (/ 1., 4., 3., 12., 13., 14., 15., 16., 17., &
                          18., 20.2, 22.99, 24.3, 26.97, 28.1, 30.97,32.06, 35.45, &
                          39.948, 39.098, 40.08, 44.95, 47.86, 50.94, 51.99, &
                          54.93, 55.845, 58.933, 58.693/)


        return
      end subroutine get_solar_params

    ! !this is to make sure the integrator does what it's supposed to
      function gaussinmod(x)
        double precision :: x,gaussinmod
        gaussinmod = 1*exp(-x**2/2.d0)!*nq
      end function gaussinmod
end module shared_mod

!Some functions that have to be external, because of the integrator.

!Just a test for the integrator. Nothing to see here
function gausstest(x)
    use shared_mod
    double precision :: x,gausstest
    gausstest = gaussinmod(x)
end function gausstest

! function dummyf(x)
!     double precision :: x, dummyf
!     dummyf = 1.d0
! end function dummyf

!   this is eqn 2.15 in 1504.04378
!This is fine as long as the escape velocity is large enough
  function capture_maximum(mx)
    use shared_mod
    implicit none
    double precision capture_maximum
    double precision, intent(in) :: mx

    capture_maximum = pi/3.d0*rho_dm/mx*radius_star**2 &
    *(exp(-3./2.*vel_sun**2/dispersion_dm**2)*sqrt(6.d0/pi)*dispersion_dm &
    + (6.d0*gm_over_r_sun/vel_sun + (dispersion_dm**2 + 3.d0*vel_sun**2)/vel_sun)*erf(sqrt(3./2.)*vel_sun/dispersion_dm))

  end function capture_maximum


!------!------!------!------!------INITIALIZATION FCT

  subroutine init_sun(solarmodel,rho0_in,usun_in,u0_in,vesc_in)
    !input velocities in km/s, not cm/s!!!
    use shared_mod
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

    vel_sun = usun_in*1.d5
    dispersion_dm =  u0_in*1.d5
    rho_dm =rho0_in
    escape_halo = vesc_in*1.d5

  end subroutine init_sun
