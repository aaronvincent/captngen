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
    
    ! integer :: pickIsotope
    ! double precision :: j_chi

    double precision, parameter :: AtomicNumber_oper(16) = (/ 1., 3., 4., 12., 14., 16., 20., 23., 24., 27., &
                                                        28., 32., 40., 40., 56., 58./) !the isotopes the catena paper uses
    character (len=4) :: isotopes(16) = [character(len=4) :: "H","He3","He4","C12","N14","O16","Ne20","Na23","Mg24", &
                                                                "Al27", "Si28","S32","Ar40","Ca40","Fe56","Ni58"]
    double precision, parameter :: AtomicSpin_oper(16) = (/ 0.5, 0.5, 0., 0., 1., 0., 0., 1.5, 0., 2.5, &
                                                        0., 0., 0., 0., 0., 0./) !spins pulled from https://physics.nist.gov/PhysRefData/Handbook/element_name.htm
    double precision :: coupling_Array(14,2)
    double precision :: W_array(8,16,2,2,7)

    integer :: q_shared
    logical :: w_shared
    
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
            !eps added to make the log finite: the value of eps does not affect the answer
            G = (p/c0)**(2.d0*mq)*mdm*w**2/(2.d0*mu**mq)*log(mu/muplus**2*w**2/(u)**2)
        endif
        GFFI_H_oper = G
    end function GFFI_H_oper
    
    function GFFI_A_oper(w,vesc,A,mq)
        double precision :: p, mu,w,vesc,u,muplus,mN,A,Ei,B
        double precision :: dgamic,GFFI_A_oper
        integer :: mq
        p = mdm*w
        mu = mdm/mnuc/A
        muplus = (1.+mu)/2.
        u = sqrt(w**2-vesc**2)
        mN = A*mnuc
        Ei  = 5.8407d-2/(mN*(0.91*mN**(1./3.)+0.3)**2)
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
end subroutine captn_init_oper

!   this is the integral over R in eqn 2.3 in 1501.03729
!THIS IS THE IMPORTANT FUNCTION: the integrand for the integral over u
function integrand_oper(u, foveru)
    use opermod
    double precision :: u, w, integrand_oper, foveru !int, vesc
    external foveru

    ! vesc = tab_vesc(ri_for_omega)

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

    ! int = vdist_over_u(u)/u*w*Omega_oper(ri_for_omega,w)
    ! integrand_oper = int
end function integrand_oper


! call captn_oper to run capt'n with the effective operator method
subroutine captn_oper(mx_in, jx_in, niso, capped)
    use opermod
    implicit none
    integer, intent(in):: niso!_in, isotopeChosen
    integer ri, eli, limit!, i
    double precision, intent(in) :: mx_in, jx_in
    double precision :: capped !this is the output
    double precision :: a, muminus, umax, umin, vesc, result
    double precision :: epsabs, epsrel, abserr, neval !for integrator
    double precision :: ier,alist,blist,rlist,elist,iord,last!for integrator
    double precision, allocatable :: u_int_res(:)

    ! specific to captn_oper
    integer :: funcType, tau, taup, term_R, term_W, q_pow, w_pow ! loop indicies
    integer :: q_functype, q_index
    double precision :: J, j_chi, RFuncConst, WFuncConst, yConverse, mu_T, prefactor_functype
    double precision :: RD, RM, RMP2, RP1, RP2, RS1, RS1D, RS2 !R functions stored in their own source files
    double precision :: prefactor_array(niso,11,2)

    dimension alist(1000),blist(1000),elist(1000),iord(1000),rlist(1000)!for integrator
    !external gausstest !this is just for testing
    external integrand_oper
    external integrand_oper_extrawterm
    ! external dummyf

    epsabs=1.d-17
    epsrel=1.d-17
    limit=1000

    mdm = mx_in
    j_chi = jx_in
    ! niso = niso_in
    
    ! pickIsotope = isotopeChosen

    ! temporary, the user will want to choose their coupling constants to match a model
    !                           c1,  c3,  c4, c5,   c6,  c7,  c8,  c9, c10, c11, c12, c13, c14, c15   
    !coupling_Array = reshape((/1.65d-8, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, &
                                !0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0/), (/14, 2/))
    
    if (.not. allocated(tab_r)) then 
        print*,"Errorface of errors: you haven't called captn_init to load the solar model!"
        return
    end if
    allocate(u_int_res(nlines))

    do eli = 1, niso
        do q_pow = 1, 11
            do w_pow = 1, 2
                prefactor_array(eli,q_pow,w_pow) = 0.d0
            end do
        end do
    end do

    !As far as I can tell, the second argument (fofuoveru) does nothing in this integrator.
    !I've sent it to an inert dummy just in case.
    
    ! First I set the entries in prefactor_array(niso,11,2)
    ! These are the constants that mulitply the corresonding integral evaluation
    do eli=1,niso
        ! I'll need mu_T to include in the prefactor when there is a v^2 term
        a = AtomicNumber_oper(eli)
        mu_T = (mnuc*a*mdm)/(mnuc*a+mdm)
        yConverse = (1.38*10.**-27) * ( ((mnuc*a)**(1./3) + 0.33) / (hbar*c0) )**2

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

                                    ! check if term_R is even (in my index convention this corresponds to it having a v^2 in the term)
                                    ! v^2 = w^2 - q^2/(2mu_T)^2
                                    if ( mod(term_R,2).eq.0 ) then
                                        ! this is the -q^2/(2mu_T)^2 contribution
                                        prefactor_array(eli,q_index+1,1) = prefactor_array(eli,q_index+1,1) - &
                                            prefactor_functype*RFuncConst*WFuncConst*yConverse**(term_W-1) * (c0**2/(4.*mu_T**2)) ! The Rfunctions are programmed with the 1/c0^2 in their v_perp^2 term (so I need to un-correct it for the- q^2/(2*mu_T)^2, and leave it be for the w^2/c^2)
                                        ! this is the +w^2 contribution
                                        prefactor_array(eli,q_index,2) = prefactor_array(eli,q_index,2) + &
                                            prefactor_functype*RFuncConst*WFuncConst*yConverse**(term_W-1)
                                        
                                    else
                                        prefactor_array(eli,q_index,1) = prefactor_array(eli,q_index,1) + &
                                            prefactor_functype*RFuncConst*WFuncConst*yConverse**(term_W-1)

                                    end if
                                end if
                            end do !term_R
                        end if
                    end do !term_W
                end do !taup
            end do !tau
        end do !functype
    end do !eli

    capped = 0.d0
    umin = 0.d0
    do ri=1,nlines
        vesc = tab_vesc(ri)
        rindex_shared = ri !make accessible via the module
        vesc_shared_arr(ri) = vesc !make accessible via the module

        do eli=1,niso
            u_int_res(ri) = 0.d0
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
                        result = 0.d0
                        q_shared = q_pow - 1
                        !Call integrator
                        call dsntdqagse(integrand_oper,vdist_over_u,umin,umax, &
                            epsabs,epsrel,limit,result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)

                        u_int_res(ri) = u_int_res(ri) + result * prefactor_array(eli,q_pow,w_pow)
                    end if
                end do !q_pow
            end do !w_pow
            u_int_res(ri) = u_int_res(ri) * (2*mnuc*a)/(2*J+1)
            u_int_res(ri) = u_int_res(ri) * NAvo*tab_starrho(ri)*tab_mfr(ri,eli)/(mnuc*a)
            u_int_res(ri) = u_int_res(ri) * (hbar*c0)**2
            capped = capped + tab_r(ri)**2*u_int_res(ri)*tab_dr(ri)

            ! capped = capped + u_int_res(ri) * tab_r(ri)**2*tab_dr(ri) * (2*mnuc*a)/(2*J+1) * NAvo*tab_starrho(ri)*tab_mfr(ri,eli)/(mnuc*a) * (hbar*c0)**2
        end do !eli
    end do !ri

    capped = 4.d0*pi*Rsun**3*capped

    if (capped .gt. 1.d100) then
      print*,"Capt'n General says: Oh my, it looks like you are capturing an  &
      infinite amount of dark matter in the Sun. Best to look into that."
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