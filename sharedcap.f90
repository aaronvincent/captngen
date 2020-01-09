!   Capt'n Shared
!   Designed as a module to house shared vaiables and functions
!   between both the General and Operator varients
!   Neal Avis Kozar 2020
!   all units of distance: cm
!   all units of mass/energy : GeV (or GeV/c^2, don't forget)
!   all units of time: seconds
!   Sticking with notation of 1504.04378. Cite that paper. Or 1605.06502 it's even better.
!   Reference q0 is 40 MeV, and v0 is 220 km/s.


module sharedcap
    implicit none
    double precision, parameter :: pi=3.141592653, NAvo=6.0221409d23, GMoverR=1.908e15
    double precision, parameter :: Rsun=69.57d9
    double precision, parameter :: c0=2.99792458d10, mnuc=0.938
    double precision, parameter :: eps=1d-10 !stops divisions by zero
    !these are now set in captn_init
    double precision :: usun , u0 ,rho0, vesc_halo
    !tab: means tabulated from file; so as not to be confused with other variables
    double precision, allocatable :: tab_mencl(:), tab_starrho(:), tab_mfr(:,:), tab_r(:), tab_vesc(:), tab_dr(:), tab_mfr_oper(:,:)

    integer :: niso, ri_for_omega, nlines
    double precision :: mdm
    
    contains

    !   this is the function f_sun(u) in 1504.04378 eqn 2.2
    !velocity distribution,
    function get_vdist(u)
        double precision :: u,get_vdist, f, normfact
        f = (3./2.)**(3./2.)*4.*rho0*u**2/sqrt(pi)/mdm/u0**3 &
            *exp(-3.*(usun**2+u**2)/(2.*u0**2))*sinh(3.*u*usun/u0**2)/(3.*u*usun/u0**2)
        !normfact = .5*erf(sqrt(3./2.)*(vesc_halo-usun)/u0) + &
        !.5*erf(sqrt(3./2.)*(vesc_halo+usun)/u0)+ u0/(sqrt(6.*pi)*usun) &
        !*(exp(-3.*(usun+vesc_halo)/2./u0**2)-exp(-3.*(usun-vesc_halo)/2./u0**2))
        normfact = 1.
        !print*,normfact
        f = f/normfact
        get_vdist=f
    end function get_vdist

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !read in solar parameters from Aldo Serenelli-style files, with header removed
    subroutine get_solar_params(filename,nlines)
        character*300 :: filename
        double precision :: Temp, Pres, Lumi !these aren't used, but dummies are required
        double precision, allocatable :: phi(:) !this is used briefly
        integer :: i,j, nlines,iostatus
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
        
        allocate(tab_mfr_oper(nlines,16)) ! for the operator method

        !now actually read in the file
        open(99,file=filename)
        do i=1,nlines
            read(99,*) tab_mencl(i),tab_r(i), Temp, tab_starrho(i), Pres, Lumi, tab_mfr(i,:)
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
        end do
        return
    end subroutine get_solar_params

    !this is to make sure the integrator does what it's supposed to
    function gaussinmod(x)
        double precision :: x,gaussinmod
        gaussinmod = 1*exp(-x**2/2.d0)
    end function gaussinmod

end module sharedcap

!Some functions that have to be external, because of the integrator.
function gausstest(x) !just a test for the integrator. Nothing to see here
    use sharedcap
    double precision :: x,gausstest
    gausstest = gaussinmod(x)
end function gausstest

function dummyf(x)
    double precision :: x, dummyf
    dummyf = 1.d0
end function dummyf

!   this is eqn 2.15 in 1504.04378
!This is fine as long as the escape velocity is large enough
subroutine captn_maxcap(mwimp_in,maxcap)
    use sharedcap
    implicit none
    double precision maxcap
    double precision, intent(in) :: mwimp_in
    mdm = mwimp_in
    maxcap = pi/3.d0*rho0/mdm*Rsun**2 &
        *(exp(-3./2.*usun**2/u0**2)*sqrt(6.d0/pi)*u0 &
        + (6.d0*GMoverR/usun + (u0**2 + 3.d0*usun**2)/usun)*erf(sqrt(3./2.)*usun/u0))
end subroutine captn_maxcap


!------!------!------!------!------INITIALIZATION FCT

subroutine captn_init(solarmodel,rho0_in,usun_in,u0_in,vesc_in)
    !input velocities in km/s, not cm/s!!!
    use sharedcap
    use iso_c_binding, only: c_ptr
    implicit none
    character (len=300) solarmodel
    double precision,intent(in) :: rho0_in,usun_in,u0_in,vesc_in
    !common solarmodel
    !external solarmodel

    if  (.not. allocated(tab_r)) then
        print*,"Capgen initializing from model: ",solarmodel
        call get_solar_params(solarmodel,nlines)
    end if
    !print*,"Capgen tabulons already allocated, you might be overdoing it by calling the init function more than once."
    usun = usun_in*1.d5
    u0 =  u0_in*1.d5
    rho0 =rho0_in
    vesc_halo = vesc_in*1.d5
end subroutine captn_init
