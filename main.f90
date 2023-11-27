! Capt'n General testing program
!
! Main capture routines can be found in gencap.f90 and opercap.f90

PROGRAM GENCAP
    implicit none
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables used for original qv scaled capt'n
    integer :: i, j, num_isotopes, spin_dependency, electron_v_nucleons, sigma_v_couplingCst
    integer :: testIndex !just for testing loops
    double precision :: mx, sigma_0, capped, maxcapture, rho0, usun, u0, vesc, maxcap, evapRate,  evapIso, evapLTE, cutoffFactor
    ! double precision :: capped_si_spec, capped_sd_spec ! Used in captn_specific()
    character(len=300) :: modfile, filename, nq_string, nv_string
    character(len=2) :: spinString(2) = [character(len=2) :: "SI", "SD"]
    character(len=9) :: outfile(7) = [character(len=9) :: "const.dat","qm1.dat","q1.dat","q2.dat","vm1.dat","v1.dat","v2.dat"]
    ! integer :: nq(7) = [integer :: 0, -1, 1, 2,  0, 0, 0] ! Choose velocity and momentum transfer powers in differential cross-section
    ! integer :: nv(7) = [integer :: 0,  0, 0, 0, -1, 1, 2]
    integer :: nq(6) = [integer :: 0, 1, 2, 0, 0, 0] ! Choose velocity and momentum transfer powers in differential cross-section
    integer :: nv(6) = [integer :: 0, 0, 0, 1, 2, -1]
    ! integer :: nq(8) = [integer :: 0, 1, 0, 2, 0, 1, 2, 3] ! Choose velocity and momentum transfer powers in differential cross-section
    ! integer :: nv(8) = [integer :: 0, 0, 1, 0, 2, 1, 1, 0]
    integer:: nq_index, nv_index !just to manually specify the powers for sigma

    double precision :: massesData(320)
    double precision :: sigmaData(320)

!-----------------------------------------------------------------------------------------------------------------------------------
! Variables used for the energy transport calculation
    integer :: transport_formalism, nlines
    double precision :: Tx, nwimpsin, noise_indicator, EtransTot, K, maxLum, EtransTotOld, sigmaPower
    ! double precision :: evapRate ! Used in fastevap()
    double precision, allocatable :: Etrans(:), mfp(:), Etrans_old(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables used in the NREO formalism calculation
    integer :: cpl
    double precision :: jx, couplingVal
    character(len=5) :: cplConsts(14) = [character(len=5) :: "c1-0", "c3-0", "c4-0",  "c5-0",  "c6-0",  "c7-0",  "c8-0", &
                                                            "c9-0", "c10-0", "c11-0", "c12-0", "c13-0", "c14-0", "c15-0"]

!-----------------------------------------------------------------------------------------------------------------------------------
! Choose solar model file:
    ! modfile = "solarmodels/model_gs98_nohead.dat"
    ! modfile = "solarmodels/struct_b16_agss09_nohead.dat"
    ! modfile = "solarmodels/struct_b16_agss09_reduce10_nohead.dat" !only every 10th entry from star radius is included
    modfile = "solarmodels/struct_b16_agss09_nohead_Tsmoothed.dat" !temperature smoothed to not nonsense

!-----------------------------------------------------------------------------------------------------------------------------------
    num_isotopes = 1 ! Number of isotopes original capt'n will loop over in the calculation: up to 29 isotopes
    spin_dependency = 1 ! 0=Spin Independent, 1=Spin Dependent
    electron_v_nucleons = 1 !0=electron-DM, 1 = nucleons-DM

!-----------------------------------------------------------------------------------------------------------------------------------
! Initialise capture calculations
    rho0 = 0.4d1 ! GeV cm^-3
    usun = 220.d0 ! km s^-1
    u0 = 270.d0 ! km s^-1
    vesc = 553.d0 ! km s^-1
    call captn_init(modfile,rho0,usun,u0,vesc)

!-----------------------------------------------------------------------------------------------------------------------------------
! Choose the energy transport formalism: Gould & Raffelt [arxiv:], Rescaled G&R to Monte Carlo [arxiv:], or Spergel & Press [arxiv:]
!       WHICH PAPERS DID THE TRANSPORT FORMALISMS COME FROM?
    transport_formalism = 3 ! 1=G&R, 2=Rescaled G&R, 3=S&P

!-----------------------------------------------------------------------------------------------------------------------------------
! Initialise transport calculations
    call getnlines(nlines)
    allocate(etrans(nlines))
    allocate(mfp(nlines))
    allocate(Etrans_old(nlines))
    call get_alpha_kappa(nq,nv)
!-----------------------------------------------------------------------------------------------------------------------------------
    ! open(77,file="/home/sberam/Documents/Masters/DM_annihilation/iceCube_tau_tauBar/iceCube_tau_tauBar.txt")
    ! open(77,file="/home/sberam/Documents/Masters/DM_annihilation/CRESST/cresst_SD.txt")
!     open(77,file="/home/sberam/Documents/Masters/updated_constraints_paper_reproduce/Fig_4_HeatMap/HeatMap_fig4.txt")
!     do i=1,320
!       read(77,*) sigmaData(i),  massesData(i)
!     end do
!     close(77)
!
!     do j =1, 1
!         nq_index = nq(j)
!         nv_index = nv(j)
!         write(nq_string , "(I1)") nq_index
!         write(nv_string , "(I1)") nv_index
!         print*, "electron_v_nucleons: ", electron_v_nucleons
!         print*, "Calculations for interaction type: nq = ", nq_index, " nv = ", nv_index
!         if (electron_v_nucleons.eq.0) then
!           filename = "Electron_nq_"//trim(nq_string)//"_nv_"//trim(nv_string)//".dat"
!         else if (electron_v_nucleons.eq.1) then
!           filename = "Hydrogen_nq_"//trim(nq_string)//"_nv_"//trim(nv_string)//".dat"
!         end if
!         open(94,file = filename)
!         write(94,*) "Number of Isotopes: ", num_isotopes
!         write(94,*) "Spin Dependency: ", spinString(spin_dependency+1)
!         write(94,*) "nq Power: ", nq_index
!         write(94,*) "nv Power: ", nv_index
!         write(94,*) "Sigma_0 | ", "DM Mass | ", "Captured Dark Matter | ",   "MaxCaptures | ", "Number of WIMPs in | ", &
!                     "Etranstot| ", "Max Lum [erg/s]| ",  "Evap Rate [s-1] | ", " evapIso | " , " evapLTE | ", " Knudsen Number |"
!         do i = 1, 1
!             print*, "------------------------------------------------------------------"
!             ! mx = 1.d1 ** (dble(i-10)/5.)
!             mx = dble(i)
!             ! mx = massesData(i)
!             ! mx = 0.1d0
!             ! sigma_0 = 10**(-31+dble(i))
!             sigma_0 = 1d-40
!             ! if(j.eq.1) then
!             !   sigma_0 = 1d-40 !10d0**(-45+dble(i)/2.)
!             ! else
!             !   sigma_0 = 1d-42 !10d0**(-45+dble(i)/2.)
!             ! end if
!             print*
!             ! call captn_general(mx, sigma_0, num_isotopes, nq_index, nv_index, spin_dependency, capped)
!             call captn_general_Rminus(mx, sigma_0, nq_index, nv_index, electron_v_nucleons, capped)
!             maxcapture = maxcap(mx)
!             print*, "sigma_0: ", sigma_0, "cm^2 ", &
!                     "mx: ", mx, "GeV ", &
!                     "Capture rate: ", capped, "s^-1 ", &
!                     "Geometric limit: ", maxcapture, "s^-1 "
!
!             nwimpsin = capped*3.d7*4.57d9
!             ! nwimpsin = 1d0
!             ! nwimpsin = 1.19d42 !nx/nb=1d-15
!             ! nwimpsin = 1d0
!             !SB: adding this method to implement a generalized way of getting the SP energy transfer from 2111.0695
!             !SB: only works for hydrogen or electrons (->num_isotopes=1)
!             call transgen_SP(sigma_0, electron_v_nucleons , nwimpsin, nq(j), nv(j), j, Tx, &
!                             noise_indicator, Etrans, Etranstot, K, maxLum)
!
!             !evaporation calculation
!             cutoffFactor = 1.1d0 !the cutoff velocity will be cutoffFactor*vesc
!             ! call evap_RPlus_total(mx, sigma_0, nwimpsin, nq_index, nv_index, electron_v_nucleons, cutoffFactor,  evapIso, &
!                                   ! evapLTE, evapRate)
!             print*, "Number of WIMPs in: ", nwimpsin, &
!                     "Energy Transport Total: ", EtransTot, "ergs/g/s"
!             write(94,*) sigma_0, mx, capped, maxcapture, nwimpsin, Etranstot, maxLum, evapRate, evapIso, evapLTE,  K
!             print*, "------------------------------------------------------------------"
!         end do
!         !close(195)
!         close(94)
!
!     end do
! print*, " "

!-----------------------------------------------------------------------------------------------------------------------------------
! Specific set up for the Non-Relativistic Effective Operator calculation [arxiv:1501.03729]
    num_isotopes = 1
    jx = 0.5
    couplingVal = 1d-3/(246.2**2.) ![GeV]-2
    call captn_init_oper()

!-----------------------------------------------------------------------------------------------------------------------------------
! Use the new NREO formalism calculation
    print*, "***** Calculations for NREO *****"
    print*, " "
    do cpl= 1, 1
        filename = "oper_"//trim(cplConsts(cpl))//".dat"
        if (electron_v_nucleons.eq.0) then
          filename = "Electron_oper_"//trim(cplConsts(cpl))//".dat"
        else if (electron_v_nucleons.eq.1) then
          filename = "Hydrogen_oper_"//trim(cplConsts(cpl))//".dat"
        end if
        open(55,file=filename)
        write(55,*) "sigma| ", "Coupling Val | ", "DM_Mass | ", "  Captures | ", "  MaxCaptures", &
                    "  Tx| ", "K| ", " MaxLum| "
        print*, " "
        print*, "Running coupling constant: ", cplConsts(cpl)
        print*
        do i = 1, 1
          ! mx = 1.d1**(dble(i-10)/5.)
          mx = 1

          if (cpl==1) then
              call populate_array(couplingVal, cpl, 0)
          else if (cpl==2) then
              call populate_array(0.d0, cpl-1, 0)
              call populate_array(couplingVal, cpl+1, 0)
          else
              call populate_array(0.d0, cpl, 0)
              call populate_array(couplingVal, cpl+1, 0)
          endif

          call captn_oper(mx, jx, num_isotopes, capped)

          maxcapture = maxcap(mx)
          nwimpsin = capped*3.d7*4.57d9

          call trans_oper_new(mx, jx, 1, nwimpsin, K, Tx, Etrans)
          print*
          print*, "couplingVal: ", couplingVal, "GeV^-2 ", &
                  "DM mass: ", mx, "GeV ", &
                  "Capture rate: ", capped, "s^-1 ", &
                  "Geometric limit: ", maxcapture, "s^-1 "
          write(55,*) sigma_0, couplingVal, mx, capped, maxcapture, Tx, K, maxLum
        end do
        close(55)
    end do

END PROGRAM GENCAP
