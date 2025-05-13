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
    !this goes with the Serenelli table format
    
    double precision, parameter :: AtomicNumber_oper(16) = (/ 1., 3., 4., 12., 14., 16., 20., 23., 24., 27., &
                                                        28., 32., 40., 40., 56., 58./) !! Atomic masses of the isotopes used in [[arXiv:1501.03729](https://arxiv.org/abs/1501.03729)].
    character (len=4) :: isotopes(16) = [character(len=4) :: "H","He3","He4","C12","N14","O16","Ne20","Na23","Mg24", &
                                                                "Al27", "Si28","S32","Ar40","Ca40","Fe56","Ni58"] !! The [[arXiv:1501.03729](https://arxiv.org/abs/1501.03729)] isotopes in text form to match against the W functions.
    double precision, parameter :: AtomicSpin_oper(16) = (/ 0.5, 0.5, 0., 0., 1., 0., 0., 1.5, 0., 2.5, &
                                                        0., 0., 0., 0., 0., 0./) !! Atomic spins of the [[arXiv:1501.03729](https://arxiv.org/abs/1501.03729)] isotopes pulled from [NIST](https://physics.nist.gov/PhysRefData/Handbook/element_name.htm).
    double precision :: coupling_Array(14,2)
    double precision :: W_array(8,16,2,2,7)
    double precision :: yConverse_array(16)

    integer :: q_shared
    logical :: w_shared
    !$OMP threadprivate(q_shared, w_shared)
    
    contains

    function GFFI_H_oper(w,vesc,mq)
        use phys, only : mnuc, c0
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
    
    function GFFI_A_oper(w,vesc,A,mq)
        use phys, only : mnuc, c0
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

    subroutine RW_prefactors(j_chi, total_prefactors)
        !! Populates the `total_prefactors` array with the numerical prefactors \(P_{i,n_q,n_w}\) for each isotope's differential
        !! cross section term corresponding to the R and W response functions as defined by:
        !! \[ \frac{\mathrm{d} \sigma}{\mathrm{d} E_R} = \frac{m_T}{2\pi w^2} \frac{4\pi}{2J+1} \sum_{\tau,\tau^\prime,k}
        !! R^{\tau\tau^\prime}_k\left({v_T^\perp}^2,\frac{q^2}{m_N^2}\right) W^{\tau\tau^\prime}_k\left(y\right) \\
        !! = \frac{2m_T}{w^2(2J+1)} \sum_{i,n_q,n_w} P_{i,n_q,n_w} q^{2n_q} w^{2n_w} \]
        !! following Eq. ([3.26](https://arxiv.org/pdf/1501.03729#equation.3.26)) and
        !! Eq. ([3.23](https://arxiv.org/pdf/1501.03729#equation.3.23)) from [[arxiv:1501.03729](https://arxiv.org/abs/1501.03729)].
        !! The bounds on the terms' powers are \(n_q = [0,8]\) and \(n_w = [0,1]\), with the largest terms
        !! (\(q^{16}w^0, q^{14}w^2\)) arising from \(\frac{q^2}{m_N^2}{v_T^\perp}^2c^{\tau}_{5}c^{\tau^\prime}_{5}\) in
        !! [\(R^{\tau\tau^\prime}_M\)](https://arxiv.org/pdf/1501.03729#equation.A.1) multiplied with
        !! \(W^{\tau\tau^\prime}_M\propto y^6\) (which can occur in isotopes
        !! [\(^{40}\text{Ar}\)](https://arxiv.org/pdf/1501.03729#equation.C.13),
        !! [\(^{40}\text{Ca}\)](https://arxiv.org/pdf/1501.03729#equation.C.14),
        !! [\(^{56}\text{Fe}\)](https://arxiv.org/pdf/1501.03729#equation.C.15), and
        !! [\(^{58}\text{Ni}\)](https://arxiv.org/pdf/1501.03729#equation.C.16)). Here \(q\) is the momentum transferred in the
        !! interaction, and \(w\) is the relative velocity between the dark matter and target nucleus. A prefactor \(P_{i,n_q,n_w}\)
        !! carries units of \(\text{GeV}^{-4-2n_q} {(\text{cm}\cdot\text{s}^{-1})}^{-2n_w}\).
        use phys, only : mnuc, c0
        double precision, intent(in):: j_chi
            !! The spin of the dark matter.
        double precision, intent(out) :: total_prefactors(:,:,:)
            !! The returned array of prefactors. It should be of size \(N_\text{isotopes}, \max(n_q)+1, \max(n_w)+1\) (Fortran
            !! arrays start with 1). This typically means `16,9,2`.

        integer :: eli, func_type, tau, tau_p, term_w, term_r ! loop indices
        integer :: q_func, q_index ! indices used in tracking the powers of momentum transfer q^{2 (q_index-1)}
        double precision :: prefactor_func, r_const, prefactor ! intermediate variables
        double precision :: rd, rm, rmp2, rp1, rp2, rs1, rs1d, rs2 ! DM response R functions stored in their own source files

        total_prefactors = 0.d0
        do eli = 1, size(total_prefactors,dim=1)
            ! I'll need the reduced mass mu to include in the prefactor when there is a v^2 term
            mu = (mnuc*AtomicNumber_oper(eli) * mdm)/(mnuc*AtomicNumber_oper(eli) + mdm)
    
            ! the current response function type in order: M, S2, S1, P2, MP2, P1, D, S1D
            do func_type = 1, 8
    
                ! contribution to q^2 count from sum over function types
                if ( func_type .lt. 4 ) then
                    q_func = 0
                    prefactor_func = 1.
                else
                    q_func = 1
                    prefactor_func = 1./mnuc**2
                end if
    
                ! the first index on each response function
                do tau = 1, 2
    
                    ! the second index on each response function
                    do tau_p = 1, 2
    
                        ! the possible y-terms for each nuclear response W function fit in order: y^0, y^1, y^2, y^3, y^4, y^5, y^6
                        do term_w = 1, 7
    
                            ! skip if the result gets multiplied by zero in the WFunction
                            if ( W_array(func_type,eli,tau,tau_p,term_w) .ne. 0.d0 ) then
    
                                ! the possible terms for each DM response R function in order: c, v2, q2, v2q2, q4, v2q4
                                do term_r = 1, 6
    
                                    ! pick appropriate constant from a given DM response R function with indices (tau,tau_p,term_r)
                                    ! note for possible future change: currently passes mnuc, and c0 - these are constants that could be shared to it through the shared module?
                                    select case (func_type)
                                    case (1)
                                        r_const =   rm(mnuc,c0,tau,tau_p,term_r-1,j_chi,coupling_Array) !!!!!!!!!!!!!!! in the DM response R functions the R term starts at zero, should change it to start at 1 like other Fortran things do for consistency
                                    case (2)
                                        r_const =  rs2(mnuc,c0,tau,tau_p,term_r-1,j_chi,coupling_Array)
                                    case (3)
                                        r_const =  rs1(mnuc,c0,tau,tau_p,term_r-1,j_chi,coupling_Array)
                                    case (4)
                                        r_const =  rp2(mnuc,tau,tau_p,term_r-1,j_chi,coupling_Array)
                                    case (5)
                                        r_const = rmp2(mnuc,tau,tau_p,term_r-1,j_chi,coupling_Array)
                                    case (6)
                                        r_const =  rp1(mnuc,tau,tau_p,term_r-1,j_chi,coupling_Array)
                                    case (7)
                                        r_const =   rd(mnuc,tau,tau_p,term_r-1,j_chi,coupling_Array)
                                    case (8)
                                        r_const = rs1d(tau,tau_p,term_r-1,j_chi,coupling_Array)
                                    case default
                                        print*, "Um, I ran out of DM response R functions to choose from?"
                                        stop
                                    end select
    
                                    ! skip if the result gets multiplied by zero in the RFunction
                                    if (r_const.ne.0.) then
    
                                        ! calculates the total number of q^2, accounting for Fortran indexing at 1
                                        ! i.e. q^{2*(q_index-1)}
                                        q_index = 1 + q_func + term_w - 1 + floor((term_r-1.)/2.)
                                        prefactor = prefactor_func * r_const &
                                            * W_array(func_type,eli,tau,tau_p,term_w) * yConverse_array(eli)**(term_w-1)
    
                                        ! check if term_r is even (in my index convention this corresponds to it having a v_perp^2
                                        ! in the DM response R function), decomposed into v_perp^2 = w^2 - q^2/(2mu)^2
                                        if ( mod(term_r,2).eq.0 ) then
                                            ! this is the -q^2/(2mu)^2 contribution (one extra q^2 compared to current q_index)
                                            total_prefactors(eli,q_index+1,1) = total_prefactors(eli,q_index+1,1) &
                                                + prefactor * (-c0**2/(4.*mu**2)) ! The DM response R functions are programmed with the 1/c0^2 in their v_perp^2 term (so I need to un-correct it for the - q^2/(2*mu_T)^2, and leave it be for the w^2/c^2)
                                            ! this is the +w^2 contribution (same q^2, but has a w^2 contribution)
                                            total_prefactors(eli,q_index,2) = total_prefactors(eli,q_index,2) + prefactor
                                        else
                                            total_prefactors(eli,q_index,1) = total_prefactors(eli,q_index,1) + prefactor
                                        end if
                                    end if
                                end do !term_r
                            end if
                        end do !term_w
                    end do !tau_p
                end do !tau
            end do !func_type
            ! total_prefactors(eli,:,:) = 2*mnuc*AtomicNumber_oper(eli)/(2*AtomicSpin_oper(eli)+1) * total_prefactors(eli,:,:)
        end do !eli
    end subroutine RW_prefactors
end module opermod

subroutine captn_init_oper()
    use opermod
    implicit none
    integer :: i, j, k, l, m
    character (len=2) :: terms(7) = [character(len=2) :: "y0", "y1", "y2", "y3", "y4", "y5", "y6"]
    real :: WM, WS2, WS1, WP2, WMP2, WP1, WD, WS1D
    
    !* @note
    ! `tab_mfr_oper` is allocated in the `get_solar_params` subroutine. Here we take `tab_mfr` and extract the isotopes used in
    ! [[arxiv:1501.03729](https://arxiv.org/abs/1501.03729)] (otherwise indices won't match on arrays). The elements beyond **Ne**
    ! are reported as a sum of all isotopes in the solar model files, so we can calculate individual isotopic abundances from Tab. 9
    ! in [[arxiv:1912.00844](https://arxiv.org/pdf/1912.00844#page=52)]. @endnote
    !!
    do i=1,nlines
        tab_mfr_oper(i,1) = tab_mfr(i,1)                ! H
        tab_mfr_oper(i,2) = tab_mfr(i,3)                ! He3
        tab_mfr_oper(i,3) = tab_mfr(i,2)                ! He4
        tab_mfr_oper(i,4) = tab_mfr(i,4)                ! C12
        tab_mfr_oper(i,5) = tab_mfr(i,6)                ! N14
        tab_mfr_oper(i,6) = tab_mfr(i,8)                ! O16
        tab_mfr_oper(i,7) = 0.931251 * tab_mfr(i,11)    ! Ne20 = 0.931251 * Ne
        tab_mfr_oper(i,8) = tab_mfr(i,12)               ! Na23 -- Na23 is the only isotope in ratio table
        tab_mfr_oper(i,9) = 0.78992 * tab_mfr(i,13)     ! Mg24 = 0.78992 * Mg
        tab_mfr_oper(i,10) = tab_mfr(i,14)              ! Al27 -- Al27 is the only isotope in ratio table
        tab_mfr_oper(i,11) = 0.9223 * tab_mfr(i,15)     ! Si28 = 0.9223 * Si
        tab_mfr_oper(i,12) = 0.9504074 * tab_mfr(i,17)  !  S32 = 0.9504074 * S
        tab_mfr_oper(i,13) = 4.d-5 * tab_mfr(i,19)      ! Ar40 = 0.00004 * Ar
        tab_mfr_oper(i,14) = 0.96941 * tab_mfr(i,21)    ! Ca40 = 0.96941 * Ca
        tab_mfr_oper(i,15) = 0.91754 * tab_mfr(i,27)    ! Fe56 = 0.91754 * Fe
        tab_mfr_oper(i,16) = 0.680769 * tab_mfr(i,29)   ! Ni58 = 0.680769 * Ni
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


! call captn_oper to run capt'n with the effective operator method
subroutine captn_oper(mx_in, jx_in, capped)!, isotopeChosen)
    use phys, only : mnuc, NAvo, hbar, c0, pi
    use opermod
    implicit none
    interface !Required unless these functions are moved to a different module file that gets compiled first
        function integrand_oper(arg1, func1)
            double precision :: arg1, integrand_oper
            interface
                function func1(arg2)
                    double precision :: arg2, func1
                end function func1
            end interface
        end function integrand_oper
    end interface
    integer ri, eli, limit!, i
    double precision, intent(in) :: mx_in, jx_in
    double precision :: capped !this is the output
    double precision :: maxcap, maxcapped, a, muminus, umax, umin, vesc, partialCapped, elementalResult, integrateResult
    double precision :: epsabs, epsrel, abserr, neval !for integrator
    double precision :: ier,alist,blist,rlist,elist,iord,last !for integrator
    ! double precision, allocatable :: u_int_res(:)
    
    ! specific to captn_oper
    integer :: q_pow, w_pow ! loop indicies
    double precision :: J, j_chi, factor_final
    double precision :: prefactor_array(size(tab_mfr_oper,dim=2),9,2)
    
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

    ! Get the prefactors for the q and v terms
    call RW_prefactors(j_chi, prefactor_array)

    ! now with all the prefactors computed, any 0.d0 entries in prefactor_array means that we can skip that integral evaluation!
    umin = 0.d0
    capped = 0.d0
    !$OMP parallel default(none) &
    !$OMP private(vesc, elementalResult, a, mu, muplus, muminus, J, umax, integrateResult, factor_final, partialCapped, &
    !$OMP   abserr,neval,ier,alist,blist,rlist,elist,iord,last) &
    !$OMP shared(nlines,mdm,vesc_halo,prefactor_array,tab_vesc,vesc_shared_arr,tab_starrho,tab_mfr_oper,tab_r,tab_dr, capped, &
    !$OMP   umin,limit,epsabs,epsrel)
    partialCapped = 0.d0
    !$OMP do
    do ri=1,nlines
        vesc = tab_vesc(ri)
        rindex_shared = ri !make accessible via the module
        vesc_shared_arr(ri) = vesc !make accessible via the module

        do eli=1,size(prefactor_array,dim=1)
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

            do w_pow=1,size(prefactor_array,dim=3)
                ! toggles whether we integrate with the w^2 term on
                w_shared = .false.
                if(w_pow.eq.2) then
                    w_shared = .true.
                end if

                do q_pow=1,size(prefactor_array,dim=2)
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
