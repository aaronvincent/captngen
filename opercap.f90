!   Capt'n Oper
!   Module to house everying specific to captn operator
!   Neal Avis Kozar 2020
!   all units of distance: cm
!   all units of mass/energy : GeV (or GeV/c^2, don't forget)
!   all units of time: seconds
!   Sticking with notation of 1504.04378. Cite that paper. Or 1605.06502 it's even better.


module opermod
    use sharedmod
    implicit none
    double precision, parameter :: hbar=6.582d-25 !GeV*s
    !this goes with the Serenelli table format

    double precision, parameter :: AtomicNumber_oper(16) = (/ 1., 3., 4., 12., 14., 16., 20., 23., 24., 27., &
                                                        28., 32., 40., 40., 56., 58./) !the isotopes the catena paper uses
    character (len=4) :: isotopes(16) = [character(len=4) :: "H","He3","He4","C12","N14","O16","Ne20","Na23","Mg24", &
                                                                "Al27", "Si28","S32","Ar40","Ca40","Fe56","Ni58"] !the isotopes in text form to match against the W functions
    double precision, parameter :: AtomicSpin_oper(16) = (/ 0.5, 0.5, 0., 0., 1., 0., 0., 1.5, 0., 2.5, &
                                                        0., 0., 0., 0., 0., 0./) !spins pulled from https://physics.nist.gov/PhysRefData/Handbook/element_name.htm
    double precision :: coupling_Array(14,2)
    double precision :: W_array(8,16,2,2,7)
    double precision :: yConverse_array(16)
    double precision :: coupling_to_sigma(14) !SB: Used to link the coupling constant to a cross section for energy transport calculations

    integer :: q_shared
    logical :: w_shared, w_flag !SB: I am using this to determine if the coupling has a dependence on v2. This will be used when converting cross sections to couplings
    !$OMP threadprivate(q_shared, w_shared)

    contains

    ! having removed the scaling momentum, are the units off here? I'm looking at the p/c0 in particular
    function GFFI_H_oper(w,vesc,mq)
        double precision :: p, mu,w,vesc,u,muplus,GFFI_H_oper,G
        integer mq
        p = mdm*w
        mu = mdm/mnuc
        muplus = (1.+mu)/2.
        u = sqrt(w**2-vesc**2)
        if (mq .ne. -1) then
            G = (p/c0)**(2.d0*mq)*mdm*w**2/(2.d0*mu**mq)*1./(1.+mq)*((mu/muplus**2)**(mq+1)-(u**2/w**2)**(mq+1))
        else
            G = (p/c0)**(2.d0*mq)*mdm*w**2/(2.d0*mu**mq)*log(mu/muplus**2*w**2/(u)**2)
        endif
        GFFI_H_oper = G
    end function GFFI_H_oper

    function GFFI_H_oper_electron(w,vesc,mq)
        double precision :: p, mu,w,vesc,u,muplus,GFFI_H_oper_electron,G
        integer mq
        p = mdm*w
        mu = mdm/melectron
        muplus = (1.+mu)/2.
        u = sqrt(w**2-vesc**2)
        if (mq .ne. -1) then
            G = (p/c0)**(2.d0*mq)*mdm*w**2/(2.d0*mu**mq)*1./(1.+mq)*((mu/muplus**2)**(mq+1)-(u**2/w**2)**(mq+1))
        else
            G = (p/c0)**(2.d0*mq)*mdm*w**2/(2.d0*mu**mq)*log(mu/muplus**2*w**2/(u)**2)
        endif
        GFFI_H_oper_electron = G
    end function GFFI_H_oper_electron

    function GFFI_A_oper(w,vesc,A,mq)
        double precision :: p, mu,w,vesc,u,muplus,mN,A,Ei,B
        double precision :: dgamic,GFFI_A_oper
        integer :: mq
        p = mdm*w
        mu = mdm/mnuc/A
        muplus = (1.+mu)/2.
        u = sqrt(w**2-vesc**2)
        mN = A*mnuc
        Ei = 1./4.d0/mN/264.114*(45.d0*A**(-1./3.)-25.d0*A**(-2./3.))
        B = .5*mdm*w**2/Ei/c0**2
        if (mq .eq. 0) then
            GFFI_A_oper = Ei*c0**2*(exp(-mdm*u**2/2/Ei/c0**2)-exp(-B*mu/muplus**2))
        else
            GFFI_A_oper = ((p)/c0)**(2*mq)*Ei*c0**2/(B*mu)**mq*(dgamic(1.+dble(mq),B*u**2/w**2) &
                - dgamic(1.+dble(mq),B*mu/muplus**2))
        end if
    end function GFFI_A_oper
end module opermod

subroutine captn_init_oper()
    use opermod
    implicit none
    integer :: i, j, k, l, m
    character (len=2) :: terms(7) = [character(len=2) :: "y0", "y1", "y2", "y3", "y4", "y5", "y6"]
    real :: WM, WS2, WS1, WP2, WMP2, WP1, WD, WS1D

    ! tab_mfr_oper is allocated in the get_solar_params subroutine
    ! take the regular array tab_mfr and extract the isotopes used in the 1501.03729 paper (otherwise indices won't match on arrays)
    do i=1,nlines
        tab_mfr_oper(i,1) = tab_mfr(i,1)
        tab_mfr_oper(i,2) = tab_mfr(i,3)
        tab_mfr_oper(i,3) = tab_mfr(i,2)
        tab_mfr_oper(i,4) = tab_mfr(i,4)
        tab_mfr_oper(i,5) = tab_mfr(i,6)
        tab_mfr_oper(i,6) = tab_mfr(i,8)
        tab_mfr_oper(i,7) = tab_mfr(i,11)
        tab_mfr_oper(i,8) = tab_mfr(i,12)
        tab_mfr_oper(i,9) = tab_mfr(i,13)
        tab_mfr_oper(i,10) = tab_mfr(i,14)
        tab_mfr_oper(i,11) = tab_mfr(i,15)
        tab_mfr_oper(i,12) = tab_mfr(i,17)
        tab_mfr_oper(i,13) = tab_mfr(i,19)
        tab_mfr_oper(i,14) = tab_mfr(i,21)
        tab_mfr_oper(i,15) = tab_mfr(i,27)
        tab_mfr_oper(i,16) = tab_mfr(i,29)
    end do

    ! this array stores each of the constants of the W polynomials from paper 1501.03729's appendix individually
    ! array index m handles the 8 varients of the W functions in order [M, S", S', P", MP", P', Delta, S'Delta]
    ! index i handles the 16 isotopes [H, He3, He4, C12, N14, O16, Ne20, Na 23, Mg24, Al27, Si28, S32, Ar40, Ca40, Fe56, Ni58]
    ! index j & k handle the two superscripts for each W function, each taking values of 0 and 1
    ! index L determines the power of each constant ranging from y^0 to y^6
    do m=1,8
        do i=1,16
            do j=1,2
                do k=1,2
                    do l=1,7
                        if (m.eq.1) then
                            W_array(m,i,j,k,l) = WM(j-1,k-1,isotopes(i),terms(l))
                        else if (m.eq.2) then
                            W_array(m,i,j,k,l) = WS2(j-1,k-1,isotopes(i),terms(l))
                        else if (m.eq.3) then
                            W_array(m,i,j,k,l) = WS1(j-1,k-1,isotopes(i),terms(l))
                        else if (m.eq.4) then
                            W_array(m,i,j,k,l) = WP2(j-1,k-1,isotopes(i),terms(l))
                        else if (m.eq.5) then
                            W_array(m,i,j,k,l) = WMP2(j-1,k-1,isotopes(i),terms(l))
                        else if (m.eq.6) then
                            W_array(m,i,j,k,l) = WP1(j-1,k-1,isotopes(i),terms(l))
                        else if (m.eq.7) then
                            W_array(m,i,j,k,l) = WD(j-1,k-1,isotopes(i),terms(l))
                        else
                            W_array(m,i,j,k,l) = WS1D(j-1,k-1,isotopes(i),terms(l))
                        end if
                    end do
                end do
            end do
        end do
    end do

    ! initiate the coupling_Array (full of the coupling constants) with all zeros
    ! populate_array will place the non-zero value into a chosen slot at runtime
    coupling_Array = reshape((/0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, &
                                0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0/), (/14, 2/))
    ! yconv comes from arxiv:1501.03729 page 10, where yconv = (b/{2 hbar c})^2
    do i = 1, 16
        yConverse_array(i) = 264.114/(45.d0*AtomicNumber_oper(i)**(-1./3.)-25.d0*AtomicNumber_oper(i)**(-2./3.))
    end do
end subroutine captn_init_oper

! this is the integral over R in eqn 2.3 in 1501.03729
! note that Omega there is expanded and broken into terms of the form const. * q^2n * exp{E_R/E_i}
! I've doen this so that I can tap into the GFFI functions in eqn 2.9 of 1504.04378
!THIS IS THE IMPORTANT FUNCTION: the integrand for the integral over u
function integrand_oper(u, foveru)
    use opermod
    implicit none
    interface
        function foveru(arg1)
            double precision :: arg1, foveru
        end function foveru
    end interface
    double precision :: u, integrand_oper
    double precision :: w

    w = sqrt(u**2+vesc_shared_arr(rindex_shared)**2)

    !Switch depending on whether we are capturing on Hydrogen or not
    if (a_shared .gt. 2.d0) then
        integrand_oper = foveru(u)*GFFI_A_oper(w,vesc_shared_arr(rindex_shared),a_shared,q_shared)
    else
        integrand_oper = foveru(u)*GFFI_H_oper(w,vesc_shared_arr(rindex_shared),q_shared)
    end if
    if (w_shared) then
        integrand_oper = integrand_oper * w**2
    end if
end function integrand_oper

function integrand_oper_electron(u, foveru)
    use opermod
    use capmod
    implicit none
    interface
        function foveru(arg1)
            double precision :: arg1, foveru
        end function foveru
    end interface
    double precision :: u, integrand_oper_electron
    double precision :: w, x
    w = sqrt(u**2+vesc_shared**2)
    if((nq.eq.0).and.(nv.eq.0)) then
      integrand_oper_electron = foveru(u)*diff_scattering_rate_cst(w, vesc_shared)
    end if
end function integrand_oper_electron

function integrand_oper_Rminus(u, foveru)
    !use opermod
    use capmod
    implicit none
    ! interface
    !     function foveru(arg1)
    !         double precision :: arg1, foveru
    !     end function foveru
    ! end interface

    double precision :: u, integrand_oper_Rminus, foveru
    double precision :: w, x

    external foveru

    w = sqrt(u**2+vesc_shared**2)
    if((nq.eq.0).and.(nv.eq.0)) then
      integrand_oper_Rminus = foveru(u)*diff_scattering_rate_cst(w, vesc_shared)
    else if (nv.eq.1) then
       integrand_oper_Rminus = foveru(u)*diff_scattering_rate_v1(w, vesc_shared)
    else if (nq.eq.1) then
      integrand_oper_Rminus = foveru(u)*diff_scattering_rate_q1(w, vesc_shared)
    else if (nq.eq.2) then
      integrand_oper_Rminus = foveru(u)*diff_scattering_rate_minus_q2(w, vesc_shared)
    end if
end function integrand_oper_Rminus

! call captn_oper to run capt'n with the effective operator method
subroutine captn_oper(mx_in, jx_in, niso, capped)!, isotopeChosen)
    use opermod
    implicit none
    interface
        function integrand_oper(arg1, func1)
            double precision :: arg1, integrand_oper
            interface
                function func1(arg2)
                    double precision :: arg2, func1
                end function func1
            end interface
        end function integrand_oper
    end interface
    integer, intent(in):: niso!, isotopeChosen
    integer ri, eli, limit!, i
    double precision, intent(in) :: mx_in, jx_in
    double precision :: capped !this is the output
    double precision :: maxcap, maxcapped, a, umax, umin, vesc, partialCapped, elementalResult, integrateResult
    double precision :: epsabs, epsrel, abserr, neval !for integrator
    double precision :: ier,alist,blist,rlist,elist,iord,last !for integrator
    ! double precision, allocatable :: u_int_res(:)

    ! specific to captn_oper
    integer :: funcType, tau, taup, term_R, term_W, q_pow, w_pow ! loop indicies
    integer :: q_functype, q_index
    double precision :: J, j_chi, RFuncConst, WFuncConst, mu_T, prefactor_functype, factor_final, prefactor_current
    double precision :: RD, RM, RMP2, RP1, RP2, RS1, RS1D, RS2 !R functions stored in their own source files
    double precision :: prefactor_array(niso,11,2)

    dimension alist(1000),blist(1000),elist(1000),iord(1000),rlist(1000)!for integrator

    epsabs=1.d-6
    epsrel=1.d-6
    limit=1000

    mdm = mx_in
    j_chi = jx_in

    if (.not. allocated(tab_r)) then
        print*,"Errorface of errors: you haven't called captn_init to load the solar model!"
        return
    end if
    ! allocate(u_int_res(nlines))

    do eli = 1, niso
        do q_pow = 1, 11
            do w_pow = 1, 2
                prefactor_array(eli,q_pow,w_pow) = 0.d0
            end do
        end do
    end do

    ! First I set the entries in prefactor_array(niso,11,2)
    ! These are the constants that mulitply the corresonding integral evaluation
    do eli=1,niso !isotopeChosen, isotopeChosen
        ! I'll need mu_T to include in the prefactor when there is a v^2 term
        a = AtomicNumber_oper(eli)
        mu_T = (mnuc*a*mdm)/(mnuc*a+mdm)

        ! the current response function type in order: M, S2, S1, P2, MP2, P1, D, S1D
        do funcType = 1,8

            ! contribution to q^2 count from sum over function types
            q_functype = 0
            prefactor_functype = 1.
            if ( functype.gt.3 ) then
                q_functype = 1
                prefactor_functype = 1./mnuc**2
            end if

            ! the first index on each response function
            do tau=1,2

                ! the second index on each response function
                do taup=1,2

                    ! the possible y-terms for each W function in order: y^0, y^1, y^2, y^3, y^4, y^5, y^6
                    do term_W = 1,7

                        WFuncConst = W_array(funcType,eli,tau,taup,term_W)

                        ! skip if the result gets multiplied by zero in the WFunction
                        if (WFuncConst.ne.0.) then

                            ! the possible terms for each R function in order: c, v2, q2, v2q2, q4, v2q4
                            do term_R = 1,6

                                ! pick out the appropriate term's constant from a given R function of tau, taup, and term_R
                                ! currently passes mnuc, and c0 - these are constants that could be shared to it through the shared module?
                                select case (funcType)
                                case (1)
                                    RFuncConst = RM(mnuc,c0,tau,taup,term_R-1,j_chi,coupling_Array) !!!!!!!!!!!!!!! in the R functions the R term starts at zero, should change it to start at 1 like other Fortran things do for consistency
                                case (2)
                                    RFuncConst = RS2(mnuc,c0,tau,taup,term_R-1,j_chi,coupling_Array)
                                case (3)
                                    RFuncConst = RS1(mnuc,c0,tau,taup,term_R-1,j_chi,coupling_Array)
                                case (4)
                                    RFuncConst = RP2(mnuc,tau,taup,term_R-1,j_chi,coupling_Array)
                                case (5)
                                    RFuncConst = RMP2(mnuc,tau,taup,term_R-1,j_chi,coupling_Array)
                                case (6)
                                    RFuncConst = RP1(mnuc,tau,taup,term_R-1,j_chi,coupling_Array)
                                case (7)
                                    RFuncConst = RD(mnuc,tau,taup,term_R-1,j_chi,coupling_Array)
                                case (8)
                                    RFuncConst = RS1D(tau,taup,term_R-1,j_chi,coupling_Array)
                                case default
                                    RFuncConst = 0.
                                    print*, "Um, I ran out of R functions to choose from?"
                                end select

                                ! skip if the result gets multiplied by zero in the RFunction
                                if (RFuncConst.ne.0.) then

                                    ! calculates the total number of q^2
                                    q_index = 1 + q_functype + term_W - 1 + floor((term_R-1.)/2.)
                                    prefactor_current = prefactor_functype*RFuncConst*WFuncConst*yConverse_array(eli)**(term_W-1)

                                    ! check if term_R is even (in my index convention this corresponds to it having a v^2 in the term)
                                    ! v^2 = w^2 - q^2/(2mu_T)^2
                                    if ( mod(term_R,2).eq.0 ) then
                                        ! this is the -q^2/(2mu_T)^2 contribution
                                        ! it has one extra q^2 contribution compared to the current W & R function contributions
                                        prefactor_array(eli,q_index+1,1) = prefactor_array(eli,q_index+1,1) - prefactor_current * &
                                            (c0**2/(4.*mu_T**2)) ! The Rfunctions are programmed with the 1/c0^2 in their v_perp^2 term (so I need to un-correct it for the- q^2/(2*mu_T)^2, and leave it be for the w^2/c^2)
                                        ! this is the +w^2 contribution
                                        ! it has the same q^2 contribution, but has a v_perp^2 contribution
                                        prefactor_array(eli,q_index,2) = prefactor_array(eli,q_index,2) + prefactor_current
                                    else
                                        prefactor_array(eli,q_index,1) = prefactor_array(eli,q_index,1) + prefactor_current

                                    end if
                                end if
                            end do !term_R
                        end if
                    end do !term_W
                end do !taup
            end do !tau
        end do !functype
    end do !eli

    ! now with all the prefactors computed, any 0.d0 entries in prefactor_array means that we can skip that integral evaluation!
    umin = 0.d0
    capped = 0.d0
    !$OMP parallel default(none) &
    !$OMP private(vesc, elementalResult, a, mu, muplus, muminus, J, umax, integrateResult, factor_final, partialCapped, &
    !$OMP   abserr,neval,ier,alist,blist,rlist,elist,iord,last) &
    !$OMP shared(nlines,niso,mdm,vesc_halo,prefactor_array,tab_vesc,vesc_shared_arr,tab_starrho,tab_mfr_oper,tab_r,tab_dr, capped, &
    !$OMP   umin,limit,epsabs,epsrel)
    partialCapped = 0.d0
    !$OMP do
    do ri=1,nlines
        vesc = tab_vesc(ri)
        rindex_shared = ri !make accessible via the module
        vesc_shared_arr(ri) = vesc !make accessible via the module

        do eli=1,niso !isotopeChosen, isotopeChosen
            ! u_int_res(ri) = 0.d0
            elementalResult = 0.d0
            a = AtomicNumber_oper(eli)
            a_shared = a !make accessible via the module

            mu = mdm/(mnuc*a)
            muplus = (1.+mu)/2.
            muminus = (mu-1.d0)/2.

            J = AtomicSpin_oper(eli)

            ! Chop the top of the integral off at the smaller of the halo escape velocity or the minimum velocity required for capture.
            umax = min(vesc * sqrt(mu)/abs(muminus), vesc_halo)

            do w_pow=1,2
                ! toggles whether we integrate with the w^2 term on
                w_shared = .false.
                if(w_pow.eq.2) then
                    w_shared = .true.
                end if

                do q_pow=1,11
                    if ( prefactor_array(eli,q_pow,w_pow).ne.0. ) then
                        integrateResult = 0.d0
                        q_shared = q_pow - 1
                        !Call integrator
                        call dsntdqagse(integrand_oper,vdist_over_u,umin,umax, &
                            epsabs,epsrel,limit,integrateResult,abserr,neval,ier,alist,blist,rlist,elist,iord,last)

                        elementalResult = elementalResult + integrateResult * prefactor_array(eli,q_pow,w_pow)
                    end if
                end do !q_pow
            end do !w_pow

            factor_final = (2*mnuc*a)/(2*J+1) * NAvo*tab_starrho(ri)*tab_mfr_oper(ri,eli)/(mnuc*a) * &
                tab_r(ri)**2*tab_dr(ri) * (hbar*c0)**2
            partialCapped = partialCapped + elementalResult * factor_final
        end do !eli
    end do !ri
    !$OMP critical
    capped = capped + partialCapped
    !$OMP end critical
    !$OMP end parallel

    capped = 4.d0*pi*Rsun**3*capped

    maxcapped = maxcap(mx_in)
    if (capped .gt. maxcapped) then
      capped = maxcapped
    end if
end subroutine captn_oper

!SB: Only works for Hydrogen + const for now (21-11-2023)
subroutine trans_oper_new(mx_in, jx_in, niso, nwimpsin, K, Tx, etransCum)!, isotopeChosen)
    use opermod
    use spergelpressmod
    implicit none
    integer, intent(in):: niso!, isotopeChosen
    integer ri, eli, limit!, i
    double precision, intent(in) :: mx_in, jx_in, nwimpsin
    double precision :: a

    integer :: funcType, tau, taup, term_R, term_W, q_pow, w_pow ! loop indicies
    integer :: q_functype, q_index
    double precision :: J, j_chi, RFuncConst, WFuncConst, mu_T, prefactor_functype, factor_final, prefactor_current
    double precision :: RD, RM, RMP2, RP1, RP2, RS1, RS1D, RS2 !R functions stored in their own source files
    double precision :: prefactor_array(niso,11,2)

    double precision :: invMFPElemental(nlines), invMFP(nlines), MFP(nlines), MeanFreePathInverseTerm(nlines) !Elemental inverse of the MFP, inverse MFP, MFP
    double precision :: nabund(nlines), etrans(nlines), etransCum(nlines)
    double precision :: K, rchi, Tx, guess_1, guess_2, reltolerance, x_1, x_2, x_3, error, f1, f2, f3
    double precision :: sigma_0, mdm_g, mtarget_g
    double precision :: GeV_cmMinus1_convert = 2.d-14, GN = 6.674d-8

    mdm = mx_in
    mdm_g =mdm*1.782662d-24 ![g]
    j_chi = jx_in

    if (.not. allocated(tab_r)) then
        print*,"Errorface of errors: you haven't called captn_init to load the solar model!"
        return
    end if

    do eli = 1, niso
        do q_pow = 1, 11
            do w_pow = 1, 2
                prefactor_array(eli,q_pow,w_pow) = 0.d0
            end do
        end do
    end do

  !************ looping over to find the prefactor for all nq and nv powers ************
    do eli=1,niso
        ! I'll need mu_T to include in the prefactor when there is a v^2 term
        a = AtomicNumber_oper(eli)
        mu_T = (mnuc*a*mdm)/(mnuc*a+mdm)

        ! the current response function type in order: M, S2, S1, P2, MP2, P1, D, S1D
        do funcType = 1,8

            ! contribution to q^2 count from sum over function types
            q_functype = 0
            prefactor_functype = 1.
            if ( functype.gt.3 ) then
                q_functype = 1
                prefactor_functype = 1./mnuc**2
            end if

            ! the first index on each response function
            do tau=1,2

                ! the second index on each response function
                do taup=1,2

                    ! the possible y-terms for each W function in order: y^0, y^1, y^2, y^3, y^4, y^5, y^6
                    do term_W = 1,7

                        WFuncConst = W_array(funcType,eli,tau,taup,term_W)

                        ! skip if the result gets multiplied by zero in the WFunction
                        if (WFuncConst.ne.0.) then

                            ! the possible terms for each R function in order: c, v2, q2, v2q2, q4, v2q4
                            do term_R = 1,6

                                ! pick out the appropriate term's constant from a given R function of tau, taup, and term_R
                                ! currently passes mnuc, and c0 - these are constants that could be shared to it through the shared module?
                                select case (funcType)
                                case (1)
                                    RFuncConst = RM(mnuc,c0,tau,taup,term_R-1,j_chi,coupling_Array) !!!!!!!!!!!!!!! in the R functions the R term starts at zero, should change it to start at 1 like other Fortran things do for consistency
                                case (2)
                                    RFuncConst = RS2(mnuc,c0,tau,taup,term_R-1,j_chi,coupling_Array)
                                case (3)
                                    RFuncConst = RS1(mnuc,c0,tau,taup,term_R-1,j_chi,coupling_Array)
                                case (4)
                                    RFuncConst = RP2(mnuc,tau,taup,term_R-1,j_chi,coupling_Array)
                                case (5)
                                    RFuncConst = RMP2(mnuc,tau,taup,term_R-1,j_chi,coupling_Array)
                                case (6)
                                    RFuncConst = RP1(mnuc,tau,taup,term_R-1,j_chi,coupling_Array)
                                case (7)
                                    RFuncConst = RD(mnuc,tau,taup,term_R-1,j_chi,coupling_Array)
                                case (8)
                                    RFuncConst = RS1D(tau,taup,term_R-1,j_chi,coupling_Array)
                                case default
                                    RFuncConst = 0.
                                    print*, "Um, I ran out of R functions to choose from?"
                                end select

                                ! skip if the result gets multiplied by zero in the RFunction
                                if (RFuncConst.ne.0.) then
                                    ! calculates the total number of q^2
                                    q_index = 1 + q_functype + term_W - 1 + floor((term_R-1.)/2.)
                                    prefactor_current = prefactor_functype*RFuncConst*WFuncConst*yConverse_array(eli)**(term_W-1)
                                    ! check if term_R is even (in my index convention this corresponds to it having a v^2 in the term)
                                    ! v^2 = w^2 - q^2/(2mu_T)^2
                                    if ( mod(term_R,2).eq.0 ) then
                                        ! this is the -q^2/(2mu_T)^2 contribution
                                        ! it has one extra q^2 contribution compared to the current W & R function contributions
                                        prefactor_array(eli,q_index+1,1) = prefactor_array(eli,q_index+1,1) - prefactor_current * &
                                            (c0**2/(4.*mu_T**2)) ! The Rfunctions are programmed with the 1/c0^2 in their v_perp^2 term
                                        ! it has the same q^2 contribution, but has a v_perp^2 contribution
                                        prefactor_array(eli,q_index,2) = prefactor_array(eli,q_index,2) + prefactor_current*c0**2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                    else
                                        prefactor_array(eli,q_index,1) = prefactor_array(eli,q_index,1) + prefactor_current

                                    end if
                                end if
                            end do !term_R
                        end if
                    end do !term_W
                end do !taup
            end do !tau
        end do !functype
    end do !eli

  !************ Calculating K ************
    invMFP = 0 !inverse of the meanfree path
    do eli=1,niso !isotopeChosen, isotopeChosen
        a = AtomicNumber_oper(eli)
        a_shared = a
        mu_T = (mnuc*a*mdm)/(mnuc*a+mdm)
        invMFPElemental = 0.d0 !inverse mean free path for a given iso

        mu = mdm/(mnuc*a)
        muplus = (1.+mu)/2.
        muminus = (mu-1.d0)/2.

        J = AtomicSpin_oper(eli)

        do w_pow=1,2
            do q_pow=1,11
                if ( prefactor_array(eli,q_pow,w_pow).ne.0. ) then
                    q_shared = q_pow - 1
                    sigma_0 = prefactor_array(eli,q_pow,w_pow)&
                               *(hbar*c0)**2*2*mu_T**2. !just linking sigma 0 to the coupling formalism. Not all sigmas have the same units!!!!!!!
                    call MeanFreePathInverse_calculate(mdm, w_pow-1, q_pow-1, sigma_0, MeanFreePathInverseTerm) !calculate the inverse mfp for a given isotope + one nq and nv pair
                    invMFPElemental = invMFPElemental + MeanFreePathInverseTerm
                end if
            end do !q_pow
        end do !w_pow
        invMFP = invMFP + invMFPElemental
    end do !eli
    MFP = 1./invMFP
    rchi = sqrt(3*kb*tab_T(1)/(2*pi*GN*tab_starrho(1)*mdm_g))
    K = MFP(1)/rchi

    !************ Calculating Tx ************
    guess_1 = maxval(tab_T)*15d0 ! One-zone WIMP temp guesses in K.
    guess_2 = maxval(tab_T)/150.d0
    reltolerance = 1.0d-8

    !starting binary search method
    x_1 = guess_1
  	x_2 = guess_2
    error = reltolerance + 1.d0

    do while (error > reltolerance)
      f1 = 0d0
      f2 = 0d0
      f3 = 0d0
      do eli=1,niso !isotopeChosen, isotopeChosen
          mtarget_g = a*mnuc*1.782662d-24 ![g]
          a = AtomicNumber_oper(eli)
          a_shared = a !make accessible via the module
          mu_T = (mnuc*a*mdm)/(mnuc*a+mdm)

          mu = mdm/(mnuc*a)
          muplus = (1.+mu)/2.
          muminus = (mu-1.d0)/2.
          nabund = tab_mfr(:,eli)*tab_starrho/mtarget_g
          do w_pow=1,2
              do q_pow=1,11
                  if ( prefactor_array(eli,q_pow,w_pow).ne.0. ) then
                      q_shared = q_pow - 1
                      sigma_0 = prefactor_array(eli,q_pow,w_pow)&
                                /GeV_cmMinus1_convert**2*2*mu_T**2. !just linking sigma 0 to the coupling formalism
                      x_3 = (x_1 + x_2)/2.d0
                      nq = q_shared
                      nv = w_pow-1
                      f1 = f1 + Tx_integral_nreo(x_1, sigma_0, mtarget_g, 1, nwimpsin, nq, nv, nabund)
                  		f2 = f2 + Tx_integral_nreo(x_2, sigma_0, mtarget_g, 1, nwimpsin, nq, nv, nabund)
                  		f3 = f3 + Tx_integral_nreo(x_3, sigma_0, mtarget_g, 1, nwimpsin, nq, nv, nabund)
                  end if
              end do !q_pow
          end do !w_pow
      end do !eli
      if (f3 == 0.d0) then
  			exit
  		else if (f1*f3 .gt. 0) then ! if f1 and f3 have the same sign
  			x_1 = x_3
  		else if (f2*f3 .gt. 0) then
  			x_2 = x_3
  		endif
  		error = abs(x_2-x_1)/x_2
    end do
    Tx = x_3
    print*, "Tx: ", Tx

    etransCum = 0d0
    do eli=1,niso !isotopeChosen, isotopeChosen
        a = AtomicNumber_oper(eli)
        a_shared = a !make accessible via the module
        mu_T = (mnuc*a*mdm)/(mnuc*a+mdm)
        mtarget_g = a*mnuc*1.782662d-24

        mu = mdm/(mnuc*a)
        muplus = (1.+mu)/2.
        muminus = (mu-1.d0)/2.

        do w_pow=1,2
            do q_pow=1,11
                if ( prefactor_array(eli,q_pow,w_pow).ne.0. ) then
                    q_shared = q_pow - 1
                    sigma_0 = prefactor_array(eli,q_pow,w_pow)&
                              /GeV_cmMinus1_convert**2*2*mu_T**2. !just linking sigma 0 to the coupling formalism
                    nq = q_shared
                    nv = w_pow-1
                    etrans = Etrans_sp_nreo(nq, nv, sigma_0, mtarget_g, 1 ,Tx, nwimpsin, nabund)
                    etransCum = etransCum + etrans
                end if
            end do !q_pow
        end do !w_pow
    end do !eli
    etransCum = 0.5/(1.d0+(0.4d0/K)**2)*etransCum
end subroutine trans_oper_new

!SB: this is a modified version of "captn_oper" subroutine to consider electrons or hydrogen ONLY NREO capture
!SB: This uses 1308.6288 equation (23) for the matrix elements
! subroutine captn_oper_Rminus(mx_in, jx_in, electron_v_nucleons, capped)!, isotopeChosen)
!     use opermod
!     use capmod
!     implicit none
!     ! interface
!     !     function integrand_oper_Rminus(arg1, func1)
!     !         double precision :: arg1, integrand_oper_Rminus
!     !         interface
!     !             function func1(arg2)
!     !                 double precision :: arg2, func1
!     !             end function func1
!     !         end interface
!     !     end function integrand_oper_Rminus
!     ! end interface
!     !integer, intent(in):: niso!, isotopeChosen
!     integer, intent(in) :: electron_v_nucleons
!     integer ri, eli, limit, x, y, z !x,y,z just iteration variables for testing
!     double precision, intent(in) :: mx_in, jx_in
!     double precision :: capped !this is the output
!     double precision :: maxcap, maxcapped, a, umax, umin, vesc, partialCapped, elementalResult, integrateResult
!     double precision :: epsabs, epsrel, abserr, neval !for integrator
!     double precision :: ier,alist,blist,rlist,elist,iord,last !for integrator
!     double precision :: mtarget, mfr(nlines), mtargetKg, test
!     ! double precision, allocatable :: u_int_res(:)
!
!     ! specific to captn_oper
!     integer :: funcType, tau, taup, term_R, term_W, q_pow, w_pow! loop indicies
!     integer :: q_functype, q_index ! q_functype denotes if RW product is multiplied by q^2 in Ptot
!     double precision :: J, j_chi, RFuncConst, WFuncConst, mu_T, prefactor_functype, factor_final, prefactor_current
!     double precision :: RD, RM, RMP2, RP1, RP2, RS1, RS1D, RS2 !R functions stored in their own source files
!
!     integer :: counter
!
!     external integrand_oper_Rminus
!     external integrand_Rminus
!
!     double precision :: c1, c1p
!     double precision :: c3, c3p
!     double precision :: c4, c4p
!     double precision :: c5, c5p
!     double precision :: c6, c6p
!     double precision :: c7, c7p
!     double precision :: c8, c8p
!     double precision :: c9, c9p
!     double precision :: c10, c10p
!     double precision :: c11, c11p
!     double precision :: c12, c12p
!     double precision :: c13, c13p
!     double precision :: c14, c14p
!     double precision :: c15, c15p
!
!     double precision :: prefactor_array(11,2)
!     character(len=300) :: NREO_prefactor_array_outname
!     character(len=5) :: cplConsts(14) = [character(len=5) :: "c1-0", "c3-0", "c4-0",  "c5-0",  "c6-0",  "c7-0",  "c8-0", &
!                                                             "c9-0", "c10-0", "c11-0", "c12-0", "c13-0", "c14-0", "c15-0"]
!
!     dimension alist(1000),blist(1000),elist(1000),iord(1000),rlist(1000)!for integrator
!
!     epsabs=1.d-6
!     epsrel=1.d-6
!     limit=1000
!
!     mdm = mx_in
!     j_chi = jx_in
!     !v0 = 22000000 !cm/s
!     q0 = v0*mdm
!
!     if (electron_v_nucleons.eq.0) then
!       mtarget = melectron
!       mfr = tab_electron_mfr
!     else if (electron_v_nucleons.eq.1) then
!       mtarget = mnuc
!       mfr = tab_mfr(:,1)
!     end if
!     mtargetKg = mtarget*1.782662d-27 ![kg]
!     if (.not. allocated(tab_r)) then
!         print*,"Errorface of errors: you haven't called captn_init to load the solar model!"
!         return
!     end if
!     ! allocate(u_int_res(nlines))
!
!     c1 = coupling_Array(1,1)
!     c3 = coupling_Array(2,1)
!     c4 = coupling_Array(3,1)
!     c5 = coupling_Array(4,1)
!     c6 = coupling_Array(5,1)
!     c7 = coupling_Array(6,1)
!     c8 = coupling_Array(7,1)
!     c9 = coupling_Array(8,1)
!     c10 = coupling_Array(9,1)
!     c11 = coupling_Array(10,1)
!     c12 = coupling_Array(11,1)
!     c13 = coupling_Array(12,1)
!     c14 = coupling_Array(13,1)
!     c15 = coupling_Array(14,1)
!     do q_pow = 1, 11
!         do w_pow = 1, 2
!             prefactor_array(q_pow,w_pow) = 0.d0
!         end do
!     end do
!
!     ! First I set the entries in prefactor_array(11,2)
!     ! These are the constants that mulitply the corresonding integral evaluation
!
!     ! I'll need mu_T to include in the prefactor when there is a v^2 term
!     a = 1
!     mu_T = (mtarget*mdm)/(mtarget+mdm)
!     !c1 = sqrt(4.*pi/mu_T**2*1d-40/(2.d-14)**2.)
!     !c11 = sqrt(1/q0**2*12.*pi*1d-42/mu_T**2*mtarget**2/(j_chi*(j_chi+1)))/(2.*1d-14)
!     !c6 = sqrt(48.*pi*1d-42/mu_T**2*mtarget**4/(j_chi*(j_chi+1)))/(2.*1d-14)
!     !c3 = sqrt(16.*pi*1d-40/v0**2/q0**2/mu_T**2*mtarget**2.)/(2.*1d-14)
!
!     ! the current response function type in order: M, S2, S1, P2, MP2, P1, D, S1D
!     ! the possible terms for each R function in order: c, v2, q2, v2q2, q4, v2q4
!     do term_R = 1,6
!         select case (term_R)
!         case (1)
!             RFuncConst = (c1**2./8. &
!                           + 1/8.*j_chi*(j_chi+1)/12.*c4**2 &
!                           + j_chi*(j_chi+1)/12./4.*c4**2.)/(pi)
!         case (2)
!               RFuncConst = (j_chi*(j_chi+1)/3./8.*c8**2 &
!                             + 1/2./4.*j_chi*(j_chi+1)/12.*c12**2. &
!                             + 1/4./8.*c7**2.)/(pi)
!         case (3)
!               RFuncConst =  (j_chi*(j_chi+1)/3./8.*c11**2 &
!                             + 1/4./8.*c10**2 &
!                             +j_chi*(j_chi+1)/12./4.*c9**2.)/(pi*mtarget**2)
!         case (4)
!             RFuncConst = (j_chi*(j_chi+1)/3./8.*c5**2 &
!                           + 1/8.*j_chi*(j_chi+1)/12.*c13**2.&
!                           + 1/4./8.*c3**2.&
!                           +j_chi*(j_chi+1)/12./4./2.*c14**2.)/(pi*mtarget**2*2)
!         case (5)
!               RFuncConst = (1/8.*j_chi*(j_chi+1)/12.*c6**2.)/(pi*mtarget**4*2)
!         case (6)
!             RFuncConst = (j_chi*(j_chi+1)/12./4./2*c15**2.)/(4*pi*mtarget**4)
!         case default
!             RFuncConst = 0.
!             print*, "Um, I ran out of R functions to choose from?"
!         end select
!
!         ! skip if the result gets multiplied by zero in the RFunction
!         if (RFuncConst.ne.0.) then
!             ! RFuncConst = 1d-40
!             ! calculates the total number of q^2
!             q_index = 1 + floor((term_R-1.)/2.)
!             prefactor_current = RFuncConst
!             ! check if term_R is even (in my index convention this corresponds to it having a v^2 in the term)
!             ! v^2 = w^2 - q^2/(2mu_T)^2
!             if (mod(term_R,2).eq.0) then
!                 ! this is the -q^2/(2mu_T)^2 contribution
!                 ! it has one extra q^2 contribution compared to the current W & R function contributions
!                 prefactor_array(q_index+1,1) = prefactor_array(q_index+1,1) - &
!                                   (q0/c0)**(2.*q_index)*prefactor_current/(2.*mu_T)**2./(muplus**(2.*q_index)*(1.+(q_index)))
!                 ! this is the +w^2 contribution
!                 ! it has the same q^2 contribution, but has a v_perp^2 contribution
!                 prefactor_array(q_index,2) = prefactor_array(q_index,2) + &
!                                 (q0/c0)**(2.*(q_index-1.))*v0**2*prefactor_current/c0**2/(muplus**(2.*(q_index &
!                                 -1.))*(1.+(q_index-1)))
!             else
!                 prefactor_array(q_index,1) = prefactor_array(q_index,1) + &
!                                           (q0/c0)**(2.*(q_index-1))*prefactor_current/(muplus**(2*(q_index-1.))*(1.+(q_index-1)))
!             end if
!         end if
!     end do !term_R
!
!     ! print*, "Prefactor array: "
!     ! print*, prefactor_array(:, 1)
!     ! print*, prefactor_array(:, 2)
!
!     ! now with all the prefactors computed, any 0.d0 entries in prefactor_array means that we can skip that integral evaluation!
!     umin = 0.d0
!     capped = 0.d0
!     a_shared = 1
!     partialCapped = 0.d0
!     do ri=1, nlines
!         vesc = tab_vesc(ri)
!         rindex_shared = ri !make accessible via the module
!         vesc_shared_arr(ri) = vesc !make accessible via the module
!         vesc_shared = vesc
!         ! u_int_res(ri) = 0.d0
!         elementalResult = 0.d0
!         ue_at_r = sqrt(2*kb_here*tab_T(ri)/mtargetKg)*100
!         mu = mdm/mtarget
!         muplus = (1.+mu)/2.
!         muminus = (mu-1.d0)/2.
!         J = 0.5
!
!         !Chop the top of the integral off at the smaller of the halo escape velocity or the minimum velocity required for capture.
!         umax = min(vesc * sqrt(mu)/abs(muminus), vesc_halo)
!
!         do w_pow=1,2
!             ! toggles whether we integrate with the w^2 term on
!             w_shared = .false.
!             if(w_pow.eq.2) then
!                 w_shared = .true.
!             end if
!
!             do q_pow=1,11
!                 if ( prefactor_array(q_pow,w_pow).ne.0. ) then
!                     integrateResult = 0.d0
!                     q_shared = q_pow - 1
!                     nq = q_shared
!                     nv = w_pow-1
!                     ! print*, "nq: ", nq
!                     ! print*, "nv: ", nv
!                     !Call integrator
!
!                     !test = reset_Rminus(mdm, mtarget, ue_at_r, vesc)
!
!                     call dsntdqagse(integrand_Rminus,vdist_over_u,umin,umax, &
!                         epsabs,epsrel,limit,integrateResult,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
!
!                     !print*, "nq: ", nq, " nv: ", nv, " integrateResult: ", integrateResult,&
!                     !        " prefactor: ", prefactor_array(q_pow,w_pow), " mu_T: ", mu_T
!                     if (integrateResult.lt.0d0) then
!                       print*, "********"
!                       print*, "nq: ", nq
!                       print*, "nv: ", nv
!                       print*, "r: ", ri
!                       print*, integrateResult
!                       print*, elementalResult
!                     end if
!                     elementalResult = elementalResult + integrateResult * prefactor_array(q_pow,w_pow)
!                     ! print*, "nq: ", nq
!                     ! print*, "nv: ", nv
!                     ! print*, "v0: ", v0
!                     ! print*, "q0: ", q0
!                     ! print*, "intResult: ", integrateResult
!                     ! print*, "prefactor: ", prefactor_array(q_pow,w_pow)
!                     ! print*, "integrateResult*prefactor_array", integrateResult*prefactor_array(q_pow,w_pow)
!                     ! print*, "ratio: ", prefactor_array(2,2)/prefactor_array(3,1)
!                     ! print*, "muplus: ", muplus
!                 end if
!             end do !q_pow
!         end do !w_pow
!         factor_final = 2.*mu_T**2.*NAvo*tab_starrho(ri)*mfr(ri)/mtarget* &
!             tab_r(ri)**2*tab_dr(ri)*(2.d-14)**2 !converting to [cm]^2
!         ! factor_final = 2*NAvo*tab_starrho(ri)*mfr(ri)/mtarget* &
!         !     tab_r(ri)**2*tab_dr(ri)*(2.d-14)**2 !converting to [cm]^2
!         partialCapped = partialCapped + elementalResult * factor_final
!     end do !ri
!     capped = capped + partialCapped
!
!     capped = 4.d0*pi*Rsun**3*capped
!     maxcapped = maxcap(mx_in)
!     ! if (capped .gt. maxcapped) then
!     !   capped = maxcapped
!     ! end if
! end subroutine captn_oper_Rminus

!SB: This calculate the inverse mean free path
subroutine MeanFreePathInverse_calculate(mx, w_pow, q_pow, sigma_0, MeanFreePathInverseTerm)
    use opermod
    ! use akmod
    use spergelpressmod
    implicit none
    integer :: i, w_pow, q_pow
    double precision :: mx,sigma_0, mdm_g
    double precision :: MeanFreePathInverseTerm(nlines), nabund(nlines), vTArray(nlines)
    double precision :: mtarget, targetmass_g, mreduced
    double precision:: GN = 6.674d-8

      mtarget = mnuc

      mreduced = mtarget*mx/(mtarget+mx) ![GeV]
      targetmass_g = mtarget*1.782662d-24  ![g]
      mdm_g = mx*1.782662d-24 ![g]

      nabund = tab_mfr(:,1)*tab_starrho/targetmass_g

      nq = q_pow
      nv = w_pow

      ! !******************Begin Inverse Mean Free Path Calc***************
      vTArray = sqrt(2.d0*kb*tab_T/mdm_g)/c0 ![vTArray] = natural units

      if ((nq.eq.0).and.(nv.eq.0)) then !constttt
        MeanFreePathInverseTerm = 2.d0 ! Since sigma_tot = 2*sigma_0 for v/q independent scattering

      else if ((nq.eq.0).and.(nv.eq.1)) then !v2
        MeanFreePathInverseTerm = 3.d0
      else if ((nq.eq.1).and.(nv.eq.0)) then !q2
        MeanFreePathInverseTerm = 6.d0

      else if ((nq.eq.1).and.(nv.eq.1)) then !v2 q2
        MeanFreePathInverseTerm = 15.d0
      else if ((nq.eq.2).and.(nv.eq.0)) then !q4
        MeanFreePathInverseTerm = 40.d0

      else if ((nq.eq.2).and.(nv.eq.1)) then !v2 q4
        MeanFreePathInverseTerm = 140.d0
      else if ((nq.eq.3).and.(nv.eq.0)) then !q6
        MeanFreePathInverseTerm = 420.d0

      else if ((nq.eq.2).and.(nv.eq.1)) then !v2 q6
        MeanFreePathInverseTerm = 1890.d0
      else if ((nq.eq.3).and.(nv.eq.0)) then !q8
        MeanFreePathInverseTerm = 6048.d0

      else if ((nq.eq.2).and.(nv.eq.1)) then !v2 q8
        MeanFreePathInverseTerm = 33264.d0
      else if ((nq.eq.3).and.(nv.eq.0)) then !q10
        MeanFreePathInverseTerm = 110880.d0

      else if ((nq.eq.2).and.(nv.eq.1)) then !v2 q10
        MeanFreePathInverseTerm = 720720.d0
      else if ((nq.eq.3).and.(nv.eq.0)) then !q12
        MeanFreePathInverseTerm = 2471040.d0

      else if ((nq.eq.2).and.(nv.eq.1)) then !v2 q12
        MeanFreePathInverseTerm = 18532800.d0
      else if ((nq.eq.3).and.(nv.eq.0)) then !q14
        MeanFreePathInverseTerm = 64864800.d0

      else if ((nq.eq.2).and.(nv.eq.1)) then !v2 q14
        MeanFreePathInverseTerm = 551350800.d0
      else if ((nq.eq.3).and.(nv.eq.0)) then !q16
        MeanFreePathInverseTerm = 1960358400.d0

      else if ((nq.eq.2).and.(nv.eq.1)) then !v2 q16
        MeanFreePathInverseTerm = 18623404800.d0
      else if ((nq.eq.3).and.(nv.eq.0)) then !q18
        MeanFreePathInverseTerm = 67044257280.d0

      else if ((nq.eq.2).and.(nv.eq.1)) then !v2 q18
        MeanFreePathInverseTerm = 703964701440.d0
      else if ((nq.eq.3).and.(nv.eq.0)) then !q20
        MeanFreePathInverseTerm = 2559871641600.d0

      else if ((nq.eq.2).and.(nv.eq.1)) then !v2 q20
        MeanFreePathInverseTerm = 29438523878400.d0
        print*, "In here bad!!!!!! "
      else if ((nq.eq.3).and.(nv.eq.0)) then !q22
        MeanFreePathInverseTerm = 107941254220800.d0
        print*, "In here bad!!!!!! "
      end if

      MeanFreePathInverseTerm = MeanFreePathInverseTerm*&
                          mx**(2*nq)*sigma_0*vTArray**(2*(nq+nv))*nabund/(1+mu)**(nq-nv)
end subroutine

subroutine populate_array(val, couple, isospin)
    ! in the 1501.03729 paper, the non-zero values chosen were 1.65*10^-8 (represented as 1.65d-8 in the code)
    ! I was trying to directly edit 'couple' and 'isospin' to use in the array indices, but Fortran was throwing segfaults when doing this
    ! might want a way to quit out of subroutine early if error is reached
    use opermod
    implicit none
    integer :: couple, isospin
    double precision :: val
    integer :: cpl, iso

    ! isospin can be 0 or 1
    if ((-1.lt.isospin).and.(isospin.lt.2)) then
        iso = isospin + 1 !fortran arrays start at 1
    else
        stop "Error: isospin can only be 0 or 1!"
    endif


    ! couple can be integer from 1 to 15, BUT 2 IS NOT ALLOWED!
    if (couple.lt.1) then
        stop "Error: you cannot pick a coupling constant lower than one!"
    else if (couple.eq.1) then
        cpl = couple
    else if (couple.eq.2) then
        stop "Error: you cannot use the second coupling constant!"
    else if (couple.gt.2) then
        cpl = couple - 1 !the coupling array doesn't have a slot for 2, so all constants other than one are shifted in row number
    else if (couple.gt.15) then
        stop "Error: you cannot pick a coupling constant past 15!"
    endif

    ! val is the value you want to populate with
    ! set the value picked in the slot chosen
    coupling_Array(cpl,iso) = val
end subroutine populate_array
