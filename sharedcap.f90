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
    double precision, allocatable :: star_enclosed(:), star_rho(:), star_fractions(:,:), star_r(:), star_escape(:), star_dr(:)
    double precision, allocatable :: star_fractions_nreo(:,:), star_temp(:), star_grav(:), tab_atomic(:), escape_shared(:)
    !this goes with the Serenelli table format
    double precision :: atomic_nums(29) !29 is is the number from the Serenelli files; if you have fewer it shouldn't matter

    integer :: nlines, shell_index_shared!, ri_for_omega
    double precision :: m_dm, vesc_shared, atomic_shared, mu, mu_plus
    !$OMP threadprivate(shell_index_shared, atomic_shared)
    
    contains

    !   this is the function f_sun(u) in 1504.04378 eqn 2.2 divided by u
    !velocity distribution,
    function distribution_over_vel(velocity)
        double precision :: velocity, distribution_over_vel, normfact
        distribution_over_vel = (3./2.)**(3./2.)*4.*rho_dm*velocity/sqrt(pi)/m_dm/dispersion_dm**3 &
        *exp(-3.*(vel_sun**2+velocity**2)/(2.*dispersion_dm**2))*sinh(3.*velocity*vel_sun/dispersion_dm**2) &
          / (3.*velocity*vel_sun/dispersion_dm**2)
        !normfact = .5*erf(sqrt(3./2.)*(escape_halo-vel_sun)/dispersion_dm) + &
        !.5*erf(sqrt(3./2.)*(escape_halo+vel_sun)/dispersion_dm)+ dispersion_dm/(sqrt(6.*pi)*vel_sun) &
        !*(exp(-3.*(vel_sun+escape_halo)/2./dispersion_dm**2)-exp(-3.*(vel_sun-escape_halo)/2./dispersion_dm**2))
        normfact = 1.
        !print*,normfact
        distribution_over_vel = distribution_over_vel/normfact
    end function distribution_over_vel

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !read in solar parameters from Aldo Serenelli-style files, with header removed
    subroutine read_solar_params(file_name, file_length)
        character*300 :: file_name
        double precision :: Pres, Lumi !these aren't used, but dummies are required
        double precision, allocatable :: phi(:) !this is used briefly
        integer :: i,j, file_length,iostatus

        radius_star = 69.57d9 !this is set here, for other stars, this sub is not called

        !Get number of lines in the file
        open(99,file=file_name)
        file_length=0
        do
          read(99,*, iostat=iostatus)
          if(iostatus/=0) then ! to avoid end of file error.
            exit
          else
            file_length=file_length+1
          end if
        end do
        close(99)
        file_length = file_length -1

        !allocate the arrays
        allocate(star_enclosed(file_length))
        allocate(star_r(file_length))
        allocate(star_rho(file_length))
        allocate(star_fractions(file_length,29)) !we could just allocate niso, but this leads to problems
        allocate(star_escape(file_length))
        allocate(phi(file_length))
        allocate(star_dr(file_length))
        allocate(star_temp(file_length)) !not used in capture_rate_qv; used for transgen (and anngen? )
        allocate(star_grav(file_length))
        allocate(star_fractions_nreo(file_length,16)) ! for the operator method
        allocate(escape_shared(file_length)) ! for OMP stuff


        !now actually read in the file
        open(99,file=file_name)
        do i=1,file_length
          read(99,*) star_enclosed(i),star_r(i), star_temp(i), star_rho(i), Pres, Lumi, star_fractions(i,:)
        end do
        close(99)

        !we calculate the escape velocity here since all the ingredients are ready
        phi(file_length) = -gm_over_r_sun
        star_escape(file_length) = sqrt(-2.d0*phi(file_length))
        star_dr(file_length) = star_r(file_length)-star_r(file_length-1)
        do i = 1,file_length-1
          j = file_length-i !trapezoid integral
          phi(j) = phi(j+1) + gm_over_r_sun*(star_r(j)-star_r(j+1))/2. &
            * (star_enclosed(j)/star_r(j)**2+star_enclosed(j+1)/star_r(j+1)**2)
          star_escape(j) = sqrt(-2.d0*phi(j)) !escape velocity in cm/s
          star_dr(j) = -star_r(j)+star_r(j+1) !while we're here, populate dr
          ! star_grav(j) = -(-phi(j)+phi(j+1))/star_dr(j)
          star_grav(i) = -gm_over_r_sun*star_enclosed(i)/star_r(i)**2/radius_star
        end do
        ! star_grav(num_lines) = star_grav(num_lines-1)
        star_grav(file_length) = -gm_over_r_sun*star_enclosed(file_length)/star_r(file_length)**2/radius_star

          ! Populate the atomic number tables here (because it relies on a specific format)
        atomic_nums  = (/ 1., 4., 3., 12., 13., 14., 15., 16., 17., &
                          18., 20.2, 22.99, 24.3, 26.97, 28.1, 30.97,32.06, 35.45, &
                          39.948, 39.098, 40.08, 44.95, 47.86, 50.94, 51.99, &
                          54.93, 55.845, 58.933, 58.693/)


        return
      end subroutine read_solar_params

    ! !this is to make sure the integrator does what it's supposed to
      function gaussian(x)
        double precision :: x,gaussian
        gaussian = 1*exp(-x**2/2.d0)!*nq
      end function gaussian
end module shared_mod

!Some functions that have to be external, because of the integrator.

!Just a test for the integrator. Nothing to see here
function gaussian_test(x)
    use shared_mod
    double precision :: x,gaussian_test
    gaussian_test = gaussian(x)
end function gaussian_test

! function dummyf(x)
!     double precision :: x, dummyf
!     dummyf = 1.d0
! end function dummyf

!   this is eqn 2.15 in 1504.04378
!This is fine as long as the escape velocity is large enough
  function capture_maximum(mass_dm)
    use shared_mod
    implicit none
    double precision capture_maximum
    double precision, intent(in) :: mass_dm

    capture_maximum = pi/3.d0*rho_dm/mass_dm*radius_star**2 &
    *(exp(-3./2.*vel_sun**2/dispersion_dm**2)*sqrt(6.d0/pi)*dispersion_dm &
    + (6.d0*gm_over_r_sun/vel_sun + (dispersion_dm**2 + 3.d0*vel_sun**2)/vel_sun)*erf(sqrt(3./2.)*vel_sun/dispersion_dm))

  end function capture_maximum


!------!------!------!------!------INITIALIZATION FCT

  subroutine init_sun(file_name, rho0, u_sun, u0, vel_esc)
    !input velocities in km/s, not cm/s!!!
    use shared_mod
    use iso_c_binding, only: c_ptr
    implicit none
    character (len=300) file_name
    double precision,intent(in) :: rho0,u_sun,u0,vel_esc
    !common file_name
    !external file_name

    if  (.not. allocated(star_r)) then !
        print*,"init_sun initializing from model: ",file_name
        call read_solar_params(file_name,nlines)
    end if

    vel_sun = u_sun*1.d5
    dispersion_dm =  u0*1.d5
    rho_dm =rho0
    escape_halo = vel_esc*1.d5

  end subroutine init_sun
