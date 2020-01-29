!   Capt'n General
!   For a continuum of Q-dependent capture
!   Simplified, general solar DM capture routine
!   Standalone code for q^2n, v^2n
!   Useful stuff is run at the end; beginning is the module that does the heavy lifting
!   Future plans: add form factor handling (a la Catena & Schwabe)
!   Made for GAMBIT, with marginal competence
!   Aaron Vincent 2017
!   Neal Avis Kozar 2018-2020
!   This houses the module and functions specific to the v&q dependent capture
!   all units of distance: cm
!   all units of mass/energy : GeV (or GeV/c^2, don't forget)
!   all units of time: seconds
!   Sticking with notation of 1504.04378. Cite that paper. Or 1605.06502 it's even better.
!   Reference q0 is 40 MeV, and v0 is 220 km/s.


module capmod
    use sharedmod
    implicit none
    double precision, parameter :: q0=0.04, v0=220.d5
    !this goes with the Serenelli table format
    double precision, parameter :: AtomicNumber(29) = (/ 1., 4., 3., 12., 13., 14., 15., 16., 17., &
                                                        18., 20.2, 22.99, 24.3, 26.97, 28.1, 30.97,32.06, 35.45, &
                                                        39.948, 39.098, 40.08, 44.95, 47.86, 50.94, 51.99, &
                                                        54.93, 55.845, 58.933, 58.693/) !29 is the max niso, corresponding to Ni
    
    ! nq and nv can be -1, 0, 1, 2; this is set in the main program
    integer :: nq, nv
    double precision :: sigma_0
    integer :: ri_for_omega
    
    contains

    !   this is eqn 2.9 in 1504.04378
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
            !eps added to make the log finite: the value of eps does not affect the answer
            G = ((p)/q0/c0)**(2.d0*nq)*mdm*w**2/(2.d0*mu**nq)*log(mu/muplus**2*w**2/(u+eps)**2)
        endif
        GFFI_H = G
    end function GFFI_H

    !   this is eqn 2.10 in 1504.04378
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
            GFFI_A = ((p+eps)/q0/c0)**(2*nq)*Ei*c0**2/(B*mu)**nq*(dgamic(1.+dble(nq),B*u**2/w**2+eps) &
                - dgamic(1.+dble(nq),B*mu/muplus**2+eps))
        end if
    end function GFFI_A


    !   this is eqn 2.4 in 1504.04378
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
end module capmod

!   this is the integral over R in eqn 2.7 in 1504.04378
!THIS IS THE IMPORTANT FUNCTION: the integrand for the integral over u
function integrand(u)
    use capmod
    double precision :: u, w, vesc, integrand, int
    vesc = tab_vesc(ri_for_omega)
    w = sqrt(u**2+vesc**2)
    !print*, NAvo, tab_starrho(ri_for_omega), mnuc, tab_mfr(ri_for_omega,1) !, Omega(ri_for_omega,w)
    int = get_vdist(u)/u*w*Omega(ri_for_omega,w)
    !print*, "omega: ", Omega(ri_for_omega,w)
    if (nv .ne. 0) then
        int = int*(w/v0)**(2*nv)
    end if
    integrand = int
end function integrand

subroutine captn_general(mx_in,sigma_0_in,niso_in,nq_in,nv_in,capped)
    use capmod
    implicit none
    integer, intent(in):: nq_in, niso_in, nv_in
    integer i, ri
    double precision, intent(in) :: mx_in, sigma_0_in
    double precision :: capped, maxcap !this is the output
    double precision :: epsabs, epsrel,limit,result,abserr,neval !for integrator
    double precision :: ier,alist,blist,rlist,elist,iord,last!for integrator
    double precision, allocatable :: u_int_res(:)

    dimension alist(1000),blist(1000),elist(1000),iord(1000),rlist(1000)!for integrator
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

    !As far as I can tell, the second argument (fofuoveru) does nothing in this integrator. I've sent it to an inert dummy just in case.
    capped = 0.d0

    do ri=1,nlines !loop over the star
        result = 0.d0
        ri_for_omega = ri !accessed via the module
        !call integrator
        call dsntdqagse(integrand,dummyf,1.d0,vesc_halo, &
            epsabs,epsrel,limit,result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
        !print*, "result: ", result
        u_int_res(ri) = result*sigma_0
        capped = capped + tab_r(ri)**2*u_int_res(ri)*tab_dr(ri)
    end do

    capped = 4.d0*pi*Rsun**3*capped

    if (capped .gt. 1.d100) then
      print*,"Capt'n General says: Oh my, it looks like you are capturing an  &
      infinite amount of dark matter in the Sun. Best to look into that."
    end if

    !this now has its own function:
    ! maxcap = pi/3.d0*rho0/mdm*Rsun**2 &
    ! *(exp(-3./2.*usun**2/u0**2)*sqrt(6.d0/pi)*u0 &
    ! + (6.d0*GMoverR/usun + (u0**2 + 3.d0*usun**2)/usun)*erf(sqrt(3./2.)*usun/u0))

    !  print*,"sigma_0 =", sigma_0, "; m = ", mdm, "; nq = ", nq, "; Capture rate: ", capped, "max = ", maxcap

    ! if (capped .gt. maxcap) then
    !     capped = maxcap
    ! end if
end subroutine captn_general

!! captn_specific calculates the capture rate for constant cross section.
subroutine captn_specific(mx_in,sigma_0_SD_in,sigma_0_SI_in,capped_SD,capped_SI)
    implicit none
    double precision, intent(in) :: mx_in, sigma_0_SD_in,sigma_0_SI_in
    double precision :: capped_SD,capped_SI

    call captn_general(mx_in,sigma_0_SD_in,1,0,0,capped_SD)
    call captn_general(mx_in,sigma_0_SI_in,29,0,0,capped_SI)
end subroutine captn_specific
