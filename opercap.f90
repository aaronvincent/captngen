!   Capt'n Oper
!   Module to house everying specific to captn operator
!   Neal Avis Kozar 2020
!   all units of distance: cm
!   all units of mass/energy : GeV (or GeV/c^2, don't forget)
!   all units of time: seconds
!   Sticking with notation of 1504.04378. Cite that paper. Or 1605.06502 it's even better.


module nreo_mod
    use shared_mod
    implicit none
    double precision, parameter :: hbar=6.582d-25 !GeV*s
    !this goes with the Serenelli table format
    
    double precision, parameter :: atomic_nums_nreo(16) = (/ 1., 3., 4., 12., 14., 16., 20., 23., 24., 27., &
                                                        28., 32., 40., 40., 56., 58./) !the isotopes the catena paper uses
    character (len=4) :: isotope_strings(16) = [character(len=4) :: "H","He3","He4","C12","N14","O16","Ne20","Na23","Mg24", &
                                                                "Al27", "Si28","S32","Ar40","Ca40","Fe56","Ni58"] !the isotopes in text form to match against the W functions
    double precision, parameter :: atomic_spins_nreo(16) = (/ 0.5, 0.5, 0., 0., 1., 0., 0., 1.5, 0., 2.5, &
                                                        0., 0., 0., 0., 0., 0./) !spins pulled from https://physics.nist.gov/PhysRefData/Handbook/element_name.htm
    double precision :: couplings_nreo(14,2)
    double precision :: nuclear_responses(8,16,2,2,7)
    double precision :: y_over_q2s(16)

    integer :: q_shared
    logical :: w_shared
    !$OMP threadprivate(q_shared, w_shared)
    
    contains

    ! having removed the scaling momentum, are the units off here? I'm looking at the p/c0 in particular
    function gffi_h_nreo(vel_dm, vel_esc, q_pow)
        double precision :: p, mu,vel_dm,vel_esc,u,muplus,gffi_h_nreo,G
        integer q_pow
        p = mdm*vel_dm
        mu = mdm/m_proton
        muplus = (1.+mu)/2.
        u = sqrt(vel_dm**2-vel_esc**2)
        if (q_pow .ne. -1) then
            G = (p/c0)**(2.d0*q_pow)*mdm*vel_dm**2/(2.d0*mu**q_pow)*1./(1.+q_pow) &
                * ((mu/muplus**2)**(q_pow+1)-(u**2/vel_dm**2)**(q_pow+1))
        else
            G = (p/c0)**(2.d0*q_pow)*mdm*vel_dm**2/(2.d0*mu**q_pow)*log(mu/muplus**2*vel_dm**2/(u)**2)
        endif
        gffi_h_nreo = G
    end function gffi_h_nreo
    
    function gffi_a_nreo(vel_dm, vel_esc, atomic_num, q_pow)
        double precision :: p, mu,vel_dm,vel_esc,u,muplus,mN,atomic_num,Ei,B
        double precision :: dgamic,gffi_a_nreo
        integer :: q_pow
        p = mdm*vel_dm
        mu = mdm/m_proton/atomic_num
        muplus = (1.+mu)/2.
        u = sqrt(vel_dm**2-vel_esc**2)
        mN = atomic_num*m_proton
        Ei = 1./4.d0/mN/264.114*(45.d0*atomic_num**(-1./3.)-25.d0*atomic_num**(-2./3.))
        B = .5*mdm*vel_dm**2/Ei/c0**2
        if (q_pow .eq. 0) then
            gffi_a_nreo = Ei*c0**2*(exp(-mdm*u**2/2/Ei/c0**2)-exp(-B*mu/muplus**2))
        else
            gffi_a_nreo = ((p)/c0)**(2*q_pow)*Ei*c0**2/(B*mu)**q_pow*(dgamic(1.+dble(q_pow),B*u**2/vel_dm**2) &
                - dgamic(1.+dble(q_pow),B*mu/muplus**2))
        end if
    end function gffi_a_nreo

    subroutine init_rw_prefactors(spin_dm, all_prefactors)
        !! Populates the `all_prefactors` array with the numerical prefactors \(P_{i,n_q,n_w}\) for each isotope's differential
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
        double precision, intent(in):: spin_dm
            !! The spin of the dark matter.
        double precision, intent(out) :: all_prefactors(:,:,:)
            !! The returned array of prefactors. It should be of size \(N_\text{isotopes}, \max(n_q)+1, \max(n_w)+1\) (Fortran
            !! arrays start with 1). This typically means `16,9,2`.

        integer :: eli, func_type, tau, tau_p, term_w, term_r ! loop indices
        integer :: q_func, q_index ! indices used in tracking the powers of momentum transfer q^{2 (q_index-1)}
        double precision :: prefactor_func, r_const, prefactor ! intermediate variables
        double precision :: rd, rm, rmp2, rp1, rp2, rs1, rs1d, rs2 ! DM response R functions stored in their own source files

        all_prefactors = 0.d0
        do eli = 1, size(all_prefactors,dim=1)
            ! I'll need the reduced mass mu to include in the prefactor when there is a v^2 term
            mu = (m_proton*atomic_nums_nreo(eli) * mdm)/(m_proton*atomic_nums_nreo(eli) + mdm)
    
            ! the current response function type in order: M, S2, S1, P2, MP2, P1, D, S1D
            do func_type = 1, 8
    
                ! contribution to q^2 count from sum over function types
                if ( func_type .lt. 4 ) then
                    q_func = 0
                    prefactor_func = 1.
                else
                    q_func = 1
                    prefactor_func = 1./m_proton**2
                end if
    
                ! the first index on each response function
                do tau = 1, 2
    
                    ! the second index on each response function
                    do tau_p = 1, 2
    
                        ! the possible y-terms for each nuclear response W function fit in order: y^0, y^1, y^2, y^3, y^4, y^5, y^6
                        do term_w = 1, 7
    
                            ! skip if the result gets multiplied by zero in the WFunction
                            if ( nuclear_responses(func_type,eli,tau,tau_p,term_w) .ne. 0.d0 ) then
    
                                ! the possible terms for each DM response R function in order: c, v2, q2, v2q2, q4, v2q4
                                do term_r = 1, 6
    
                                    ! pick appropriate constant from a given DM response R function with indices (tau,tau_p,term_r)
                                    ! note for possible future change: currently passes m_nuc, and c0 - these are constants that could be shared to it through the shared module?
                                    select case (func_type)
                                    case (1)
                                        r_const =   rm(m_proton,c0,tau,tau_p,term_r-1,spin_dm,couplings_nreo) !!!!!!!!!!!!!!! in the DM response R functions the R term starts at zero, should change it to start at 1 like other Fortran things do for consistency
                                    case (2)
                                        r_const =  rs2(m_proton,c0,tau,tau_p,term_r-1,spin_dm,couplings_nreo)
                                    case (3)
                                        r_const =  rs1(m_proton,c0,tau,tau_p,term_r-1,spin_dm,couplings_nreo)
                                    case (4)
                                        r_const =  rp2(m_proton,tau,tau_p,term_r-1,spin_dm,couplings_nreo)
                                    case (5)
                                        r_const = rmp2(m_proton,tau,tau_p,term_r-1,spin_dm,couplings_nreo)
                                    case (6)
                                        r_const =  rp1(m_proton,tau,tau_p,term_r-1,spin_dm,couplings_nreo)
                                    case (7)
                                        r_const =   rd(m_proton,tau,tau_p,term_r-1,spin_dm,couplings_nreo)
                                    case (8)
                                        r_const = rs1d(tau,tau_p,term_r-1,spin_dm,couplings_nreo)
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
                                            * nuclear_responses(func_type,eli,tau,tau_p,term_w) * y_over_q2s(eli)**(term_w-1)
    
                                        ! check if term_r is even (in my index convention this corresponds to it having a v_perp^2
                                        ! in the DM response R function), decomposed into v_perp^2 = w^2 - q^2/(2mu)^2
                                        if ( mod(term_r,2).eq.0 ) then
                                            ! this is the -q^2/(2mu)^2 contribution (one extra q^2 compared to current q_index)
                                            all_prefactors(eli,q_index+1,1) = all_prefactors(eli,q_index+1,1) &
                                                + prefactor * (-c0**2/(4.*mu**2)) ! The DM response R functions are programmed with the 1/c0^2 in their v_perp^2 term (so I need to un-correct it for the - q^2/(2*mu_T)^2, and leave it be for the w^2/c^2)
                                            ! this is the +w^2 contribution (same q^2, but has a w^2 contribution)
                                            all_prefactors(eli,q_index,2) = all_prefactors(eli,q_index,2) + prefactor
                                        else
                                            all_prefactors(eli,q_index,1) = all_prefactors(eli,q_index,1) + prefactor
                                        end if
                                    end if
                                end do !term_r
                            end if
                        end do !term_w
                    end do !tau_p
                end do !tau
            end do !func_type
        end do !eli
    end subroutine init_rw_prefactors
end module nreo_mod

subroutine init_nreo()
    use nreo_mod
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
                            nuclear_responses(m,i,j,k,l) = WM(j-1,k-1,isotope_strings(i),terms(l))
                        else if (m.eq.2) then
                            nuclear_responses(m,i,j,k,l) = WS2(j-1,k-1,isotope_strings(i),terms(l))
                        else if (m.eq.3) then
                            nuclear_responses(m,i,j,k,l) = WS1(j-1,k-1,isotope_strings(i),terms(l))
                        else if (m.eq.4) then
                            nuclear_responses(m,i,j,k,l) = WP2(j-1,k-1,isotope_strings(i),terms(l))
                        else if (m.eq.5) then
                            nuclear_responses(m,i,j,k,l) = WMP2(j-1,k-1,isotope_strings(i),terms(l))
                        else if (m.eq.6) then
                            nuclear_responses(m,i,j,k,l) = WP1(j-1,k-1,isotope_strings(i),terms(l))
                        else if (m.eq.7) then
                            nuclear_responses(m,i,j,k,l) = WD(j-1,k-1,isotope_strings(i),terms(l))
                        else
                            nuclear_responses(m,i,j,k,l) = WS1D(j-1,k-1,isotope_strings(i),terms(l))
                        end if
                    end do
                end do
            end do
        end do
    end do

    ! initiate the couplings_nreo (full of the coupling constants) with all zeros
    !! [[init_coupling]] will place the non-zero value into a chosen slot at runtime
    couplings_nreo = reshape((/0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, &
                                0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0/), (/14, 2/))

    ! Comes from arxiv:1501.03729 page 10, where yconv = (b/{2 hbar c})^2
    do i = 1, 16
        y_over_q2s(i) = 264.114/(45.d0*atomic_nums_nreo(i)**(-1./3.)-25.d0*atomic_nums_nreo(i)**(-2./3.))
    end do
end subroutine init_nreo

! this is the integral over R in eqn 2.3 in 1501.03729
! note that Omega there is expanded and broken into terms of the form const. * q^2n * exp{E_R/E_i}
! I've doen this so that I can tap into the GFFI functions in eqn 2.9 of 1504.04378
!THIS IS THE IMPORTANT FUNCTION: the integrand for the integral over u
function velocity_integrand_nreo(init_velocity, dist_over_vel)
    use nreo_mod
    implicit none
    interface
        function dist_over_vel(arg1)
            double precision :: arg1, dist_over_vel
        end function dist_over_vel
    end interface
    double precision :: init_velocity, velocity_integrand_nreo
    double precision :: w

    w = sqrt(init_velocity**2+vesc_shared_arr(rindex_shared)**2)

    !Switch depending on whether we are capturing on Hydrogen or not
    if (a_shared .gt. 2.d0) then
        velocity_integrand_nreo = dist_over_vel(init_velocity)*gffi_a_nreo(w,vesc_shared_arr(rindex_shared),a_shared,q_shared)
    else
        velocity_integrand_nreo = dist_over_vel(init_velocity)*gffi_h_nreo(w,vesc_shared_arr(rindex_shared),q_shared)
    end if
    if (w_shared) then
        velocity_integrand_nreo = velocity_integrand_nreo * w**2
    end if

end function velocity_integrand_nreo


! call capture_rate_nreo to run capt'n with the effective operator method
subroutine capture_rate_nreo(m_dm, spin_dm, capture_rate)!, isotopeChosen)
    use nreo_mod
    implicit none
    interface !Required unless these functions are moved to a different module file that gets compiled first
        function velocity_integrand_nreo(arg1, func1)
            double precision :: arg1, velocity_integrand_nreo
            interface
                function func1(arg2)
                    double precision :: arg2, func1
                end function func1
            end interface
        end function velocity_integrand_nreo
    end interface
    integer ri, eli, limit!, i
    double precision, intent(in) :: m_dm, spin_dm
    double precision, intent(out) :: capture_rate !this is the output
    double precision :: capture_maximum, maxcapped, a, muminus, umax, umin, vesc, partialCapped, elementalResult, integrateResult
    double precision :: epsabs, epsrel, abserr, neval !for integrator
    double precision :: ier,alist,blist,rlist,elist,iord,last !for integrator
    ! double precision, allocatable :: u_int_res(:)
    
    ! specific to capture_rate_nreo
    integer :: q_pow, w_pow ! loop indicies
    double precision :: J, factor_final
    double precision :: prefactor_array(size(tab_mfr_oper,dim=2),9,2)
    
    dimension alist(1000),blist(1000),elist(1000),iord(1000),rlist(1000)!for integrator
    
    epsabs=1.d-6
    epsrel=1.d-6
    limit=1000

    mdm = m_dm
    
    if (.not. allocated(tab_r)) then 
        print*,"Errorface of errors: you haven't called init_sun to load the solar model!"
        return
    end if
    ! allocate(u_int_res(nlines))

    ! Get the prefactors for the q and v terms
    call init_rw_prefactors(spin_dm, prefactor_array)

    ! now with all the prefactors computed, any 0.d0 entries in prefactor_array means that we can skip that integral evaluation!
    umin = 0.d0
    capture_rate = 0.d0
    !$OMP parallel default(none) &
    !$OMP private(vesc, elementalResult, a, mu, muplus, muminus, J, umax, integrateResult, factor_final, partialCapped, &
    !$OMP   abserr,neval,ier,alist,blist,rlist,elist,iord,last) &
    !$OMP shared(nlines,mdm,escape_halo,prefactor_array,tab_vesc,vesc_shared_arr,tab_starrho,tab_mfr_oper,tab_r,tab_dr, &
    !$OMP   capture_rate,umin,limit,epsabs,epsrel)
    partialCapped = 0.d0
    !$OMP do
    do ri=1,nlines
        vesc = tab_vesc(ri)
        rindex_shared = ri !make accessible via the module
        vesc_shared_arr(ri) = vesc !make accessible via the module

        do eli=1,size(prefactor_array,dim=1)
            ! u_int_res(ri) = 0.d0
            elementalResult = 0.d0
            a = atomic_nums_nreo(eli)
            a_shared = a !make accessible via the module

            mu = mdm/(m_proton*a)
            muplus = (1.+mu)/2.
            muminus = (mu-1.d0)/2.

            J = atomic_spins_nreo(eli)

            ! Chop the top of the integral off at the smaller of the halo escape velocity or the minimum velocity required for capture.
            umax = min(vesc * sqrt(mu)/abs(muminus), escape_halo)

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
                        call dsntdqagse(velocity_integrand_nreo,vdist_over_u,umin,umax, &
                            epsabs,epsrel,limit,integrateResult,abserr,neval,ier,alist,blist,rlist,elist,iord,last)

                        elementalResult = elementalResult + integrateResult * prefactor_array(eli,q_pow,w_pow)
                    end if
                end do !q_pow
            end do !w_pow

            factor_final = (2*m_proton*a)/(2*J+1) * avogadro*tab_starrho(ri)*tab_mfr_oper(ri,eli)/(m_proton*a) * &
                tab_r(ri)**2*tab_dr(ri) * (hbar*c0)**2
            partialCapped = partialCapped + elementalResult * factor_final
        end do !eli
    end do !ri
    !$OMP critical
    capture_rate = capture_rate + partialCapped
    !$OMP end critical
    !$OMP end parallel

    capture_rate = 4.d0*pi*radius_star**3*capture_rate

    maxcapped = capture_maximum(m_dm)
    if (capture_rate .gt. maxcapped) then
      capture_rate = maxcapped
    end if
end subroutine capture_rate_nreo

subroutine init_coupling(value, coupling_index, isospin_index)
    ! in the 1501.03729 paper, the non-zero values chosen were 1.65*10^-8 (represented as 1.65d-8 in the code)
    ! I was trying to directly edit 'coupling_index' and 'isospin_index' to use in the array indices, but Fortran was throwing segfaults when doing this
    ! might want a way to quit out of subroutine early if error is reached
    use nreo_mod
    implicit none
    integer, intent(in) :: coupling_index, isospin_index
    double precision, intent(in) :: value
    integer :: cpl, iso

    ! isospin_index can be 0 or 1
    if ((-1.lt.isospin_index).and.(isospin_index.lt.2)) then
        iso = isospin_index + 1 !fortran arrays start at 1
    else
        stop "Error: isospin_index can only be 0 or 1!"
    endif


    ! coupling_index can be integer from 1 to 15, BUT 2 IS NOT ALLOWED!
    if (coupling_index.lt.1) then
        stop "Error: you cannot pick a coupling constant lower than one!"
    else if (coupling_index.eq.1) then
        cpl = coupling_index
    else if (coupling_index.eq.2) then
        stop "Error: you cannot use the second coupling constant!"
    else if (coupling_index.gt.2) then
        cpl = coupling_index - 1 !the coupling array doesn't have a slot for 2, so all constants other than one are shifted in row number
    else if (coupling_index.gt.15) then
        stop "Error: you cannot pick a coupling constant past 15!"
    endif

    ! set the value picked in the slot chosen
    couplings_nreo(cpl,iso) = value
end subroutine init_coupling
