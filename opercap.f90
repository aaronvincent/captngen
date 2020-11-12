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
    
    ! ! breaks the function W into each term, and sums them with the corresponding GFFI
    ! function sumW(w,vesc,iso,tau,tauprime,Wtype,qOffset)
    !     double precision :: w,vesc,yConverse,sumW,tally
    !     integer :: k,iso,tau,tauprime,Wtype,qOffset
    !     double precision :: G
    !     ! y = yConverse * q^2
    !     ! yConverse is the conversion factor to go from q^2 to y
    !     yConverse = 2/3.*((0.91*(mnuc*AtomicNumber_oper(iso))**(1./3)+0.3)*10**-13)**2/(2*hbar*c0)**2
    !     tally = 0
    !     !$OMP parallel default(shared) private(k, G)
    !     !$OMP do
    !     do k=1,7
    !         if (iso.eq.1) then
    !             G = GFFI_H_oper(w,vesc,(k-1+qOffset))
    !         else
    !             G = GFFI_A_oper(w,vesc,AtomicNumber_oper(iso),(k-1+qOffset))
    !         end if
    !         tally = tally + W_array(Wtype,iso,tau,tauprime,k) * yConverse**(k-1) * G
    !     end do
    !     !$OMP end do
    !     !$OMP end parallel
    !     sumW = tally
    ! end function sumW
    
    ! ! this is eqn 3.23 in 1501.03729
    ! ! large sum handled through expansion of R functions
    ! ! many if statements used to check for terms excluded by choice of constants (c1,c3..c15) = 0
    ! function p_tot(w,vesc,i)
    !     double precision :: w,vesc, p_tot
    !     double precision :: m_T,mu_T,GF
    !     integer :: i,tau,taup
    !     integer :: c, v2, q2, v2q2, q4, v2q4
    !     double precision :: RD, RM, RMP2, RP1, RP2, RS1, RS1D, RS2
        
    !     c = 0
    !     v2 = 1
    !     q2 = 2
    !     v2q2 = 3
    !     q4 = 4
    !     v2q4 = 5

    !     ! get the target nucleus mass m_T and the nucleus-dm reduced mass mu_T
    !     m_T = mnuc * AtomicNumber_oper(i)
    !     mu_T = (m_T*mdm)/(m_T+mdm)

    !     ! use GFFI_H for hydrogen or y~0
    !     ! do a check to see if y can simplify W to a constant
        
    !     ! the sum in eqn 3.23 of paper 1501.03729
    !     ! expanded the R functions into each term of consatnt, v^2n, and q^2n
    !     ! check to see if choice of couping constants has zeroed out any R term before doinng any more work with it
    !     ! if not zerod out, is added to sum (sumW expands the W functions as polynomials as given in paper's appendix)
    !     ! note! coupling constants are 2d array of length 14 (c1,c3,c4...c14,c15) (note absence of c2)
    !     ! this results in array call of index [1] -> c1, but a call of index [2] -> c3 !
    !     p_tot = 0.0
    !     !$OMP  parallel default(shared) private(tau)
    !     !$OMP  do
    !     do tau=1,2
    !         !$OMP  parallel default(shared) private(taup)
    !         !$OMP  do
    !         do taup=1,2
    !             ! RM (c, v2, q2, v2q2)
    !             ! c1,c1
    !             if ((coupling_Array(1,tau).ne.0).and.(coupling_Array(1,taup).ne.0)) then
    !                 p_tot = p_tot + RM(mnuc,c0,tau,taup,c,j_chi,coupling_Array) * sumW(w,vesc,i,tau,taup,1,0)
    !             end if
    !             ! c8,c8
    !             if ((coupling_Array(7,tau).ne.0).and.(coupling_Array(7,taup).ne.0)) then
    !                 p_tot = p_tot + RM(mnuc,c0,tau,taup,v2,j_chi,coupling_Array) * w**2 * sumW(w,vesc,i,tau,taup,1,0) &
    !                             - RM(mnuc,c0,tau,taup,v2,j_chi,coupling_Array) * c0**2/(4.*mu_T**2) * sumW(w,vesc,i,tau,taup,1,1) 
    !             end if
    !             ! c11,c11
    !             if ((coupling_Array(10,tau).ne.0).and.(coupling_Array(10,taup).ne.0)) then
    !                 p_tot = p_tot + RM(mnuc,c0,tau,taup,q2,j_chi,coupling_Array) * sumW(w,vesc,i,tau,taup,1,1)
    !             end if
    !             ! c5,c5
    !             if ((coupling_Array(4,tau).ne.0).and.(coupling_Array(4,taup).ne.0)) then
    !                 p_tot = p_tot + RM(mnuc,c0,tau,taup,v2q2,j_chi,coupling_Array) * w**2 * sumW(w,vesc,i,tau,taup,1,1) &
    !                             - RM(mnuc,c0,tau,taup,v2q2,j_chi,coupling_Array) * c0**2/(4.*mu_T**2) * sumW(w,vesc,i,tau,taup,1,2) 
    !             end if
                
    !             ! RS2 (c, v2, q2, v2q2, q4)
    !             ! c4,c4
    !             if ((coupling_Array(3,tau).ne.0).and.(coupling_Array(3,taup).ne.0)) then
    !                 p_tot = p_tot + RS2(mnuc,c0,tau,taup,c,j_chi,coupling_Array) * sumW(w,vesc,i,tau,taup,2,0) 
    !             end if
    !             ! c12,c12
    !             if ((coupling_Array(11,tau).ne.0).and.(coupling_Array(11,taup).ne.0)) then
    !                 p_tot = p_tot + RS2(mnuc,c0,tau,taup,v2,j_chi,coupling_Array) * w**2 * sumW(w,vesc,i,tau,taup,2,0) &
    !                             - RS2(mnuc,c0,tau,taup,v2,j_chi,coupling_Array) * c0**2/(4.*mu_T**2) * sumW(w,vesc,i,tau,taup,2,1)
    !             end if
    !             ! c10,c10 or c4,c6 or c6,c4
    !             if (((coupling_Array(9,tau).ne.0).and.(coupling_Array(9,taup).ne.0)).or. &
    !                     ((coupling_Array(3,tau).ne.0).and.(coupling_Array(5,taup).ne.0)).or. &
    !                     ((coupling_Array(5,tau).ne.0).and.(coupling_Array(3,taup).ne.0))) then
    !                 p_tot = p_tot + RS2(mnuc,c0,tau,taup,q2,j_chi,coupling_Array) * sumW(w,vesc,i,tau,taup,2,1)
    !             end if
    !             ! c13,c13
    !             if ((coupling_Array(12,tau).ne.0).and.(coupling_Array(12,taup).ne.0)) then
    !                 p_tot = p_tot + RS2(mnuc,c0,tau,taup,v2q2,j_chi,coupling_Array) * w**2 * sumW(w,vesc,i,tau,taup,2,1) &
    !                             - RS2(mnuc,c0,tau,taup,v2q2,j_chi,coupling_Array) * c0**2/(4.*mu_T**2) * sumW(w,vesc,i,tau,taup,2,2)
    !             end if
    !             ! c6,c6
    !             if ((coupling_Array(5,tau).ne.0).and.(coupling_Array(5,taup).ne.0)) then
    !                 p_tot = p_tot + RS2(mnuc,c0,tau,taup,q4,j_chi,coupling_Array) * sumW(w,vesc,i,tau,taup,2,2)
    !             end if
                
    !             ! RS1 (c, v2, q2, v2q2, v2q4)
    !             ! c4,c4
    !             if ((coupling_Array(3,tau).ne.0).and.(coupling_Array(3,taup).ne.0)) then
    !                 p_tot = p_tot + RS1(mnuc,c0,tau,taup,c,j_chi,coupling_Array) * sumW(w,vesc,i,tau,taup,3,0)
    !             end if
    !             ! c7,c7 or c12,c12
    !             if (((coupling_Array(6,tau).ne.0).and.(coupling_Array(6,taup).ne.0)).or. &
    !                     ((coupling_Array(11,tau).ne.0).and.(coupling_Array(11,taup).ne.0))) then
    !                 p_tot = p_tot + RS1(mnuc,c0,tau,taup,v2,j_chi,coupling_Array) * w**2 * sumW(w,vesc,i,tau,taup,3,0) &
    !                             - RS1(mnuc,c0,tau,taup,v2,j_chi,coupling_Array) * c0**2/(4.*mu_T**2) * sumW(w,vesc,i,tau,taup,3,1)
    !             end if
    !             ! c9,c9
    !             if ((coupling_Array(8,tau).ne.0).and.(coupling_Array(8,taup).ne.0)) then
    !                 p_tot = p_tot + RS1(mnuc,c0,tau,taup,q2,j_chi,coupling_Array) * sumW(w,vesc,i,tau,taup,3,1)
    !             end if
    !             ! c3,c3 or c14,c14 or c15,c12 or c12,c15
    !             if (((coupling_Array(2,tau).ne.0).and.(coupling_Array(2,taup).ne.0)).or. &
    !                     ((coupling_Array(13,tau).ne.0).and.(coupling_Array(13,taup).ne.0)).or. &
    !                     ((coupling_Array(14,tau).ne.0).and.(coupling_Array(11,taup).ne.0)).or. &
    !                     ((coupling_Array(11,tau).ne.0).and.(coupling_Array(14,taup).ne.0))) then
    !                 p_tot = p_tot + RS1(mnuc,c0,tau,taup,v2q2,j_chi,coupling_Array)  * w**2 * sumW(w,vesc,i,tau,taup,3,1)&
    !                             -RS1(mnuc,c0,tau,taup,v2q2,j_chi,coupling_Array)  * c0**2/(4.*mu_T**2) * sumW(w,vesc,i,tau,taup,3,2)
    !             end if
    !             ! c15,c15
    !             if ((coupling_Array(14,tau).ne.0).and.(coupling_Array(14,taup).ne.0)) then
    !                 p_tot = p_tot + RS1(mnuc,c0,tau,taup,v2q4,j_chi,coupling_Array) * w**2 * sumW(w,vesc,i,tau,taup,3,2) &
    !                             - RS1(mnuc,c0,tau,taup,v2q4,j_chi,coupling_Array) * c0**2/(4.*mu_T**2) * sumW(w,vesc,i,tau,taup,3,3) 
    !             end if
                
    !             ! RP2 (c, q2, q4)
    !             ! c12,c12
    !             if (((coupling_Array(11,tau).ne.0).and.(coupling_Array(11,taup).ne.0))) then
    !                 p_tot = p_tot + 1./mnuc**2 * RP2(mnuc,tau,taup,c,j_chi,coupling_Array) * sumW(w,vesc,i,tau,taup,4,1)
    !             end if
    !             ! c3,c3 or c12,c15 or c15,c12
    !             if (((coupling_Array(2,tau).ne.0).and.(coupling_Array(2,taup).ne.0)).or. &
    !                     ((coupling_Array(11,tau).ne.0).and.(coupling_Array(14,taup).ne.0)).or. &
    !                     ((coupling_Array(14,tau).ne.0).and.(coupling_Array(11,taup).ne.0))) then
    !                 p_tot = p_tot + 1./mnuc**2 * RP2(mnuc,tau,taup,q2,j_chi,coupling_Array) * sumW(w,vesc,i,tau,taup,4,2)
    !             end if
    !             ! c15,c15
    !             if (((coupling_Array(14,tau).ne.0).and.(coupling_Array(14,taup).ne.0))) then
    !                 p_tot = p_tot + 1./mnuc**2 * RP2(mnuc,tau,taup,q4,j_chi,coupling_Array) * sumW(w,vesc,i,tau,taup,4,3)
    !             end if
                
    !             ! RMP2 (c, q2)
    !             ! c3,c1 or c12,c11
    !             if (((coupling_Array(2,tau).ne.0).and.(coupling_Array(1,taup).ne.0)).or. &
    !                     ((coupling_Array(11,tau).ne.0).and.(coupling_Array(10,taup).ne.0))) then
    !                 p_tot = p_tot + 1./mnuc**2 * RMP2(mnuc,tau,taup,c,j_chi,coupling_Array) * sumW(w,vesc,i,tau,taup,5,1)
    !             end if
    !             ! c15,c11
    !             if (((coupling_Array(14,tau).ne.0).and.(coupling_Array(10,taup).ne.0))) then
    !                 p_tot = p_tot + 1./mnuc**2 * RMP2(mnuc,tau,taup,q2,j_chi,coupling_Array) * sumW(w,vesc,i,tau,taup,5,2)
    !             end if
                
    !             ! RP1 (c, q2)
    !             ! c12,c12
    !             if (((coupling_Array(11,tau).ne.0).and.(coupling_Array(11,taup).ne.0))) then
    !                 p_tot = p_tot + 1./mnuc**2 * RP1(mnuc,tau,taup,c,j_chi,coupling_Array) * sumW(w,vesc,i,tau,taup,6,1) 
    !             end if
    !             ! c13,c13
    !             if (((coupling_Array(12,tau).ne.0).and.(coupling_Array(12,taup).ne.0))) then
    !                 p_tot = p_tot + 1./mnuc**2 * RP1(mnuc,tau,taup,q2,j_chi,coupling_Array) * sumW(w,vesc,i,tau,taup,6,2)
    !             end if
                
    !             ! RD (c, q2)
    !             ! c8,c8
    !             if ((coupling_Array(7,tau).ne.0).and.(coupling_Array(7,taup).ne.0)) then
    !                 p_tot = p_tot + 1./mnuc**2 * RD(mnuc,tau,taup,c,j_chi,coupling_Array) * sumW(w,vesc,i,tau,taup,7,1)
    !             end if
    !             ! c5,c5
    !             if ((coupling_Array(4,tau).ne.0).and.(coupling_Array(4,taup).ne.0)) then
    !                 p_tot = p_tot + 1./mnuc**2 * RD(mnuc,tau,taup,q2,j_chi,coupling_Array) * sumW(w,vesc,i,tau,taup,7,2)
    !             end if
                
    !             !RS1D (c)
    !             ! c5,c4 or c8,c9
    !             if (((coupling_Array(4,tau).ne.0).and.(coupling_Array(3,taup).ne.0)).or. &
    !                     ((coupling_Array(7,tau).ne.0).and.(coupling_Array(8,taup).ne.0))) then
    !                 p_tot = p_tot + 1./mnuc**2 * RS1D(tau,taup,c,j_chi,coupling_Array) * sumW(w,vesc,i,tau,taup,8,1)
    !             end if
    !         end do
    !         !$OMP  end do
    !         !$OMP  end parallel
    !     end do
    !     !$OMP  end do
    !     !$OMP  end parallel
    !     p_tot = p_tot  * hbar**2 * c0**2
    ! end function p_tot


    ! !   this is eqn 2.1 in 1501.03729
    ! function OMEGA_oper(rindex,w)
    !     double precision :: w, vesc,mu,muplus,u,J, OMEGA_oper
    !     integer rindex, i
        
    !     vesc = tab_vesc(rindex)
    !     u = sqrt(w**2-vesc**2)
    !     Omega_oper = 0.d0
        
    !     ! if only one isotope is being run
    !     if ((0 .lt. pickIsotope).and.(pickIsotope .lt. 17)) then
    !         J = AtomicSpin_oper(pickIsotope)
    !         mu = mdm/mnuc/AtomicNumber_oper(pickIsotope)
    !         muplus = (1.+mu)/2.
    !         if (mu/muplus**2 .gt. u**2/w**2) then
    !             ! note the factor [(2*mnuc*AtomicNumber(pickIsotope))/(w**2*(2*J+1))],is the product of
    !             ! the factors in front of eqns 3.26 and 3.23 in paper 1501.03729
    !             Omega_oper = Omega_oper+(NAvo*tab_starrho(rindex)/AtomicNumber_oper(pickIsotope)/mnuc* &
    !                 tab_mfr_oper(rindex,pickIsotope))*w*((2*mnuc*AtomicNumber_oper(pickIsotope))/(w**2*(2*J+1)))* &
    !                 p_tot(w,vesc,pickIsotope)
    !         end if
    !     else ! if all the isotopes are being run
    !         !$OMP  parallel default(none) shared(niso, mdm, u, w, tab_starrho, rindex, tab_mfr_oper, vesc) &
    !         !$OMP  private(i, J, mu, muplus)
    !         !$OMP  do
    !         do i = 1,16!niso
    !             J = AtomicSpin_oper(i)
    !             mu = mdm/mnuc/AtomicNumber_oper(i)
    !             muplus = (1.+mu)/2.
    !             if (mu/muplus**2 .gt. u**2/w**2) then
    !                 ! note the factor [(2*mnuc*AtomicNumber(i))/(w**2*(2*J+1))],is the product of
    !                 ! the factors in front of eqns 3.26 and 3.23 in paper 1501.03729
    !                 Omega_oper = Omega_oper+(NAvo*tab_starrho(rindex)/AtomicNumber_oper(i)/mnuc*tab_mfr_oper(rindex,i))*w* &
    !                                         ((2*mnuc*AtomicNumber_oper(i))/(w**2*(2*J+1)))*p_tot(w,vesc,i)
    !             end if
    !         end do
    !         !$OMP  end do
    !         !$OMP  end parallel
    !     end if
    ! end function OMEGA_oper
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
function integrand_oper(u)
    use opermod
    implicit none
    double precision :: u, w, integrand_oper, foveru !int, vesc
    external foveru

    ! vesc = tab_vesc(ri_for_omega)

    w = sqrt(u**2+vesc_shared**2)

    !Switch depending on whether we are capturing on Hydrogen or not
    if (a_shared .gt. 2.d0) then
        integrand_oper = foveru(u)*GFFI_A_oper(w,vesc_shared,a_shared,q_shared)
    else
        integrand_oper = foveru(u)*GFFI_H_oper(w,vesc_shared,q_shared)
    end if
    if (w_shared) then
        integrand_oper = integrand_oper * w**2
    end if

    ! int = vdist_over_u(u)/u*w*Omega_oper(ri_for_omega,w)
    ! integrand_oper = int
end function integrand_oper


!   this is the integral over R in eqn 2.3 in 1501.03729
!THIS IS THE IMPORTANT FUNCTION: the integrand for the integral over u
function integrand_oper_extrawterm(u,foveru)
    use opermod
    double precision :: u, w, integrand_oper_extrawterm, foveru
    external foveru

    w = sqrt(u**2+vesc_shared**2)

    !Switch depending on whether we are capturing on Hydrogen or not
    if (a_shared .gt. 2.d0) then
        integrand_oper_extrawterm = foveru(u)*GFFI_A_oper(w,vesc_shared,a_shared,q_shared+1)
    else
        integrand_oper_extrawterm = foveru(u)*GFFI_H_oper(w,vesc_shared,q_shared+1)
    end if
end function integrand_oper_extrawterm


! call captn_oper to run capt'n with the effective operator method
subroutine captn_oper(mx_in, jx_in, niso, isotopeChosen, capped)
    use opermod
    implicit none
    integer, intent(in):: niso!_in, isotopeChosen
    integer ri, eli, limit!, i
    double precision, intent(in) :: mx_in, jx_in
    double precision :: capped !this is the output
    double precision :: maxcap, maxcapped, a, muminus, umax, umin, vesc, result
    double precision :: epsabs, epsrel, abserr, neval !for integrator
    double precision :: ier,alist,blist,rlist,elist,iord,last!for integrator
    double precision, allocatable :: u_int_res(:)

    integer :: isotopeChosen ! remove later, i'm just too lazy right now to change any call to captn_oper and remove that parameter

    ! specific to captn_oper
    integer :: funcType, tau, taup, term_R, term_W ! loop indicies
    integer :: qOffset
    double precision :: J, j_chi, RFuncConst, WFuncConst, yConverse, mu_T
    double precision :: RD, RM, RMP2, RP1, RP2, RS1, RS1D, RS2 !R functions stored in their own source files

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

    !As far as I can tell, the second argument (fofuoveru) does nothing in this integrator.
    !I've sent it to an inert dummy just in case.
    
    
    capped = 0.d0
    ! Bottom part of the integral is always zero -- happy little slow DM particles can always be captured.
    umin = 0.d0
    ! the current response function type in order: M, S2, S1, P2, MP2, P1, D, S1D
    do funcType = 1,8
        ! the first index on each response function
        do tau=1,2
            ! the second index on each response function
            do taup=1,2
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

                    !!!!!!!!!!!!!!!!!!!!!! if I check here to see if RFuncConst is zero, I could skip the integrator call, no?
                    !!!!!!!!!!!!!!!!!!!!!! what ever pops out of the integrator is just multiplied by this term, so if it's zero it'll kill the integrator's result off (and move on to the next iteration of the loop)
                    if (RFuncConst.ne.0.) then
                        qOffset = 0 ! qOffset is used to tell the GFFI how many extra q^2 terms come from the current R function term
                        if (funcType.gt.3) then ! note that the functions following SigmaPrime in the P_tot sum are all multiplied by an extra factor of q^2/M_N 
                            qOffset = qOffset + 1
                        end if
                        if ((term_R.eq.3) .or. (term_R.eq.4)) then ! add on the single q^2 term
                            qOffset = qOffset + 1
                        else if ((term_R.eq.5) .or. (term_R.eq.6)) then ! add on the two q^2 terms
                            qOffset = qOffset + 2
                        end if
                        
                        w_shared = .false.
                        if((term_R.eq.2).or.(term_R.eq.4).or.(term_R.eq.6)) then ! if the current R term is one with a v_T^2 in it, this accounts for the v^2/c^2 upon expansion (we deal with the negative just after this do loop)
                            w_shared = .true.
                        end if

                        !Loop over the different elements
                        do eli = 1, niso
                            a = AtomicNumber_oper(eli)
                            a_shared = a !make accessible via the module

                            mu = mx_in/(mnuc*a)
                            muplus = (1.+mu)/2.
                            muminus = (mu-1.d0)/2.

                            J = AtomicSpin_oper(eli)
                            !!!!!!!!!!!!!!!!!!!!!!!!!simplify yConverse mathematical opperation to help loop speed?
                            yConverse = 2/3.*((0.91*(mnuc*a)**(1./3)+0.3)*10**-13)**2/(2*hbar*c0)**2 !the conversion factor between q and y: y = yConverse * q^2
                            ! get the target nucleus-dm reduced mass mu_T
                            mu_T = (mnuc*a*mdm)/(mnuc*a+mdm)

                            ! for all terms, do this base calculation
                            ! the possible y-terms for each W function in order: y^0, y^1, y^2, y^3, y^4, y^5, y^6
                            do term_W = 1,7
                                WFuncConst = W_array(funcType,eli,tau,taup,term_W)
                                ! skip if the result gets multiplied by zero in the WFunction
                                if (WFuncConst.ne.0.) then
                                    q_shared = term_W-1+qOffset ! the GFFI call in the integrand needs to know this

                                    !Loop over the shells of constant radius in the star
                                    do ri = 1, nlines
                                        vesc = tab_vesc(ri)
                                        vesc_shared = vesc !make accessible via the module
                                        ! Chop the top of the integral off at the smaller of the halo escape velocity or the minimum velocity required for capture.
                                        umax = min(vesc * sqrt(mu)/abs(muminus), vesc_halo)

                                        u_int_res(ri) = 0.
                                        !Call integrator
                                        call dsntdqagse(integrand_oper,vdist_over_u,umin,umax, &
                                        epsabs,epsrel,limit,u_int_res(ri),abserr,neval,ier,alist,blist,rlist,elist,iord,last)

                                        ! now we deal with the extra negative term mentioned just earlier
                                        if(w_shared) then 
                                            result = 0.
                                            !Call integrator
                                            call dsntdqagse(integrand_oper_extrawterm,vdist_over_u,umin,umax, &
                                            epsabs,epsrel,limit,result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
                                            u_int_res(ri) = u_int_res(ri) - c0**2/(4.*mu_T**2) * result
                                        end if
                                        !!!!!!!!!!!!!!!!!!!!!!!!!!! the rest of these are multiplicative factors pulled out of the integral (and ideally canceled out to be minimally impactful on performance)
                                        ! extra terms from sumW
                                        u_int_res(ri) = u_int_res(ri) * WFuncConst * yConverse**(term_W-1)
                                        ! extra terms from ptot
                                        if (funcType.gt.3) then ! note that the functions following Sigma` in the P_tot sum are all multiplied by an extra factor of q^2/mnuc^2 
                                            u_int_res(ri) = u_int_res(ri) * 1/mnuc**2
                                        end if
                                        u_int_res(ri) = u_int_res(ri) * RFuncConst * (hbar*c0)**2
                                        ! extra terms from omega_oper
                                        u_int_res(ri) = u_int_res(ri) * (2*mnuc*a)/(2*J+1) !!!!!!!!!!!!!!!!!! we have a rogue 'w' here too, flag it for the integrand!
                                        !!!!!!!!!!!!!!!!!!!!!extra terms from captngeneral - check these, it seems Pat has different terms than what I get from my math, specifically this 2.d0*(muplus/mu)**2
                                        u_int_res(ri) = u_int_res(ri)*2.d0*NAvo*tab_starrho(ri)*tab_mfr(ri,eli)*(muplus/mx_in)**2
                                        ! u_int_res(ri) = u_int_res(ri)*NAvo*tab_starrho(ri)*tab_mfr(ri,eli)/(mnuc*a)

                                        ! !!!! trying to collect all multiplicative factors here
                                        ! u_int_res(ri) = u_int_res(ri) * WFuncConst*yConverse**(term_W-1) &
                                        !     * RFuncConst*(hbar*c0)**2 * 2/(2*J+1) &
                                        !     * NAvo*tab_starrho(ri)*tab_mfr(ri,eli)
                                        ! if (funcType.gt.3) then ! note that the functions following Sigma` in the P_tot sum are all multiplied by an extra factor of q^2/mnuc^2 
                                        !     u_int_res(ri) = u_int_res(ri) * 1/mnuc**2
                                        ! end if

                                        capped = capped + tab_r(ri)**2*u_int_res(ri)*tab_dr(ri)

                                        if (isnan(capped)) then
                                            capped = 0.d0
                                            stop 'NaN encountered whilst trying compute capture rate.'
                                        end if
                                    end do !ri
                                end if !WFuncConst>0
                            end do !term_W
                        end do !eli
                    end if !RFuncConst>0
                end do !term_R
            end do !taup
        end do !tau
    end do !funcType


    ! ! completes integral (2.3) in paper 1501.03729 (gives dC/dV as fn of radius)
    ! do ri=1,nlines !loop over the star
    !     result = 0.d0
    !     ri_for_omega = ri !accessed via the module
    !     !call integrator
    !     call dsntdqagse(integrand_oper,dummyf,1.d0,vesc_halo, &
    !         epsabs,epsrel,limit,result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
    !     u_int_res(ri) = result
    !     capped = capped + tab_r(ri)**2*u_int_res(ri)*tab_dr(ri)
    ! end do

    !completes integral (2.4) of paper 1501.03729
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
        print*,"Error: isospin can only be 0 or 1!"
    endif


    ! couple can be integer from 1 to 15, BUT 2 IS NOT ALLOWED!
    if (couple.lt.1) then
        print*,"Error: you cannot pick a coupling constant lower than one!"
    else if (couple.eq.1) then
        cpl = couple
    else if (couple.eq.2) then
        print*,"Error: you cannot use the second coupling constant!"
    else if (couple.gt.2) then
        cpl = couple - 1 !the coupling array doesn't have a slot for 2, so all constants other than one are shifted in row number
    else if (couple.gt.15) then
        print*,"Error: you cannot pick a coupling constant past 15!"
    endif

    ! val is the value you want to populate with
    ! set the value picked in the slot chosen
    coupling_Array(cpl,iso) = val
end subroutine populate_array