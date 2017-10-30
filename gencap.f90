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


    module capmod
    implicit none
    double precision, parameter :: pi=3.141592653, NAvo=6.0221409d23,GMoverR = 1.908e15
    double precision, parameter :: Rsun = 69.57d9
    double precision, parameter :: c0 = 2.99792458d10, mnuc = 0.938, q0 = 0.04,v0 = 220.d5
    !these are now set in captn_init
    double precision :: usun , u0 ,rho0,vesc_halo
    !this goes with the Serenelli table format
    double precision, parameter :: AtomicNumber(29) = (/ 1., 4., 3., 12., 13., 14., 15., 16., 17., &
                                                        18., 20.2, 22.99, 24.3, 26.97, 28.1, 30.97,32.06, 35.45, &
                                                        39.948, 39.098, 40.08, 44.95, 47.86, 50.94, 51.99, &
                                                        54.93, 55.845, 58.933, 58.693/) !29 is the max niso, corresponding to Ni
    !tab: means tabulated from file; so as not to be confused with other variables
    double precision, allocatable :: tab_mencl(:),tab_starrho(:),tab_mfr(:,:), tab_r(:), tab_vesc(:), tab_dr(:)

    ! nq and nv can be -1, 0, 1, 2; this is set in the main program
    integer :: nq, nv, niso, ri_for_omega, nlines
    double precision :: mdm, sigma_0


    contains


    !velocity distribution,
    function get_vdist(u)
    double precision :: u,get_vdist, f, normfact
    f = (3./2.)**(3./2.)*4.*rho0*u**2/sqrt(pi)/mdm/u0**3 &
    *exp(-3.*(usun**2+u**2)/(2.*u0**2))*sinh(3.*u*usun/u0**2)/(3.*u*usun/u0**2)
!    normfact = .5*erf(sqrt(3./2.)*(vesc_halo-usun)/u0) + &
!    .5*erf(sqrt(3./2.)*(vesc_halo+usun)/u0)+ u0/(sqrt(6.*pi)*usun) &
!    *(exp(-3.*(usun+vesc_halo)/2./u0**2)-exp(-3.*(usun-vesc_halo)/2./u0**2))
    normfact = 1.
!    print*,normfact
    f = f/normfact
    get_vdist=f
    end function get_vdist


     !generalized form factor: hydrogen
    function GFFI_H(w,vesc)
    double precision :: p, mu,w,vesc,u,muplus,GFFI_H,G
    p = mdm*w
    mu = mdm/mnuc
    muplus = (1.+mu)/2.
    u = sqrt(w**2-vesc**2)
    if (nq .ne. -1) then
    G = (p/q0/c0)**(2.d0*nq)*mdm*w**2/(2.d0*mu**nq)*1./(1.+nq)*((mu/muplus**2)**(nq+1)-(u**2/w**2)**(nq+1))
    else
    G = (p/q0/c0)**(2.d0*nq)*mdm*w**2/(2.d0*mu**nq)*log(mu/muplus**2*w**2/u**2)
    endif
    GFFI_H = G
    end function GFFI_H


    !generalized form factor: other elements
    function GFFI_A(w,vesc,A)
    double precision :: p, mu,w,vesc,u,muplus,mN,A,Ei,B
    double precision :: dgamic,GFFI_A
    p = mdm*w
    mu = mdm/mnuc/A
    muplus = (1.+mu)/2.
    u = sqrt(w**2-vesc**2)
    mN = A*mnuc
    Ei  = 5.8407d-2/(mN*(0.91*mN**(1./3.)+0.3)**2)
    B = .5*mdm*w**2/Ei/c0**2
    if (nq .eq. 0) then
    GFFI_A = Ei*c0**2*(exp(-mdm*u**2/2/Ei/c0**2)-exp(-B*mu/muplus**2))
    else
    GFFI_A = (p/q0/c0)**(2*nq)*Ei*c0**2/(B*mu)**nq*(dgamic(1.+dble(nq),B*u**2/w**2) &
    - dgamic(1.+dble(nq),B*mu/muplus**2))
    end if
    end function GFFI_A


    !this is omega/sigma_0
    function OMEGA(rindex,w)
    double precision :: sigma_N, GF,vesc,Omega,mu,muplus,muminus,u,w
    integer i, rindex
    vesc = tab_vesc(rindex)
    u = sqrt(w**2-vesc**2)
    Omega = 0.d0

    do i = 1,niso
    !this is fine for SD as long as it's just hydrogen. Otherwise, spins must be added
    sigma_N = AtomicNumber(i)**4*(mdm+mnuc)**2/(mdm+AtomicNumber(i)*mnuc)**2
    !hydrogen


    mu = mdm/mnuc/AtomicNumber(i)
    muplus = (1.+mu)/2.
    muminus = (mu-1.d0)/2.
    if (mu*vesc**2/muminus**2 .gt. u**2) then

    if (i .eq. 1) then
    GF = GFFI_H(w,vesc)
    else
    GF = GFFI_A(w,vesc,AtomicNumber(i))
    end if
    Omega = Omega+sigma_N*NAvo*tab_starrho(rindex)/AtomicNumber(i)/mnuc*tab_mfr(rindex,i)*muplus**2/mu*GF

    end if
    end do
    Omega = Omega*2.d0/mdm/w

    end function omega


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
    end

! this is to make sure the integrator does what it's supposed to
    function gaussinmod(x)

    double precision :: x,gaussinmod
    gaussinmod = nq*exp(-x**2/2.d0)
    end function gaussinmod

!end get_solar_params
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end module capmod



!   Some functions that have to be external, because of the integrator.
    function gausstest(x) !just a test for the integrator. Nothing to see here
    use capmod
    double precision :: x,gausstest
    gausstest = gaussinmod(x)
    end function gausstest


    !THIS IS THE IMPORTANT FUNCTION: the integrand for the integral over u
    function integrand(u)
    use capmod
    double precision :: u, w, vesc, integrand, int
    vesc = tab_vesc(ri_for_omega)
    w = sqrt(u**2+vesc**2)
!    print*, NAvo, tab_starrho(ri_for_omega), mnuc, tab_mfr(ri_for_omega,1) !, Omega(ri_for_omega,w)
    int = get_vdist(u)/u*w*Omega(ri_for_omega,w)
    if (nv .ne. 0) then
    int = int*(w/v0)**(2*nv)
    end if
    integrand = int
    end function integrand


    function dummyf(x)
    double precision :: x, dummyf
    dummyf = 1.d0
    end function dummyf

    subroutine captn_general(mx_in,sigma_0_in,niso_in,nq_in,nv_in,capped)
    use capmod
    implicit none
    integer :: nq_in, niso_in, nv_in
    integer i, ri
    double precision :: mx_in, sigma_0_in
    double precision :: capped, maxcap !this is the output
    double precision :: epsabs, epsrel,limit,result,abserr,neval !for integrator
    double precision :: ier,alist,blist,rlist,elist,iord,last!for integrator
    double precision, allocatable :: u_int_res(:)

    dimension alist(1000),blist(1000),elist(1000),iord(1000),   rlist(1000)!for integrator
    external gausstest !this is just for testing
    external integrand
    external dummyf
    epsabs=1.d-17
    epsrel=1.d-17
    limit=1000


    mdm = mx_in
    sigma_0 = sigma_0_in
    niso = niso_in
    nq = nq_in
    nv = nv_in



    if (nq*nv .ne. 0) then
    print*, "Oh no! nq and nv can't both be nonzero. "
    return
    end if

    if (.not. allocated(tab_r)) then !
        print*,"Errorface of errors: you haven't called captn_init to load the solar model!"
        return
    end if
        allocate(u_int_res(nlines))

!    As far as I can tell, the second argument (fofuoveru) does nothing in this integrator. I've sent it to an inert dummy just in case.
    capped = 0.d0

    do ri=1,nlines !loop over the star
    result = 0.d0
    ri_for_omega = ri !accessed via the module
    !call integrator
    call dsntdqagse(integrand,dummyf,1.d0,vesc_halo, &
    epsabs,epsrel,limit,result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
    u_int_res(ri) = result*sigma_0
    capped = capped + tab_r(ri)**2*u_int_res(ri)*tab_dr(ri)
    end do
    capped = 4.d0*pi*Rsun**3*capped

    maxcap = pi/3.d0*rho0/mdm*Rsun**2 &
    *(exp(-3./2.*usun**2/u0**2)*sqrt(6.d0/pi)*u0 &
    + (6.d0*GMoverR/usun + (u0**2 + 3.d0*usun**2)/usun)*erf(sqrt(3./2.)*usun/u0))

!    print*,"sigma_0 =", sigma_0, "; m = ", mdm, "; nq = ", nq, "; Capture rate: ", capped, "max = ", maxcap

    if (capped .gt. maxcap) then
        capped = maxcap
    end if


!    open(55,file = "captest.dat")
!    do i=1,nlines
!    write(55,*) tab_r(i), u_int_res(i)
!    end do
!    close(55)
!    end do

    end subroutine captn_general
!


    !! captn_specific calculates the capture rate for constant cross section.
    subroutine captn_specific(mx_in,sigma_0_SD_in,sigma_0_SI_in,capped_SD,capped_SI)
      implicit none
    double precision, intent(in) :: mx_in, sigma_0_SD_in,sigma_0_SI_in
    double precision :: capped_SD,capped_SI

    call captn_general(mx_in,sigma_0_SD_in,1,0,0,capped_SD)
    call captn_general(mx_in,sigma_0_SI_in,29,0,0,capped_SI)

    end subroutine captn_specific

!------!------!------!------!------INITIALIZATION FCT

    subroutine captn_init(solarmodel,rho0_in,usun_in,u0_in,vesc_in)
    use capmod
    use iso_c_binding, only: c_ptr
    implicit none
    character (len=300) solarmodel
    double precision,intent(in) :: rho0_in,usun_in,u0_in,vesc_in
!    common solarmodel
!    external solarmodel


    if  (.not. allocated(tab_r)) then !
        print*,"Capgen initializing from model: ",solarmodel
        call get_solar_params(solarmodel,nlines)
    end if
    ! print*,"Capgen tabulons already allocated, you might be overdoing it by calling the init function more than once."
    usun = usun_in
    u0 =  u0_in
    rho0 =rho0_in
    vesc_halo = vesc_in
    
    end subroutine captn_init
! moved to main.f90:
!    PROGRAM GENCAP
!    implicit none
!    double precision :: mx, sigma_0, capped,capped_sd(250),capped_si(250)
!    integer :: niso, nq, nv, i
!
!    niso = 29
!    nq = 1
!    nv = 0
!    mx = 1000.d0
!    sigma_0 = 1.0d-40 !cm^2
!    do i = 1,250
!    mx = 10**(.02*i - 0.02)
!!    print*,mx
!    call captn_general(mx,sigma_0,29,nq,nv,capped_si(i))
!    call captn_general(mx,sigma_0,1,nq,nv,capped_sd(i))
!
!    end do
!
!    open(55,file = "captest_agss_q2.dat")
!    do i=1,250
!    write(55,*) 10**(.02*i - 0.02), capped_sd(i),capped_si(i)
!    end do
!    close(55)
!
!    END PROGRAM GENCAP
!
