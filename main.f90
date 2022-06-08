! Capt'n General testing program
!
! Main capture routines can be found in gencap.f90 and opercap.f90

PROGRAM GENCAP
    implicit none
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables used for original qv scaled capt'n
    integer :: i, j, num_isotopes, spin_dependency
    double precision :: mx, sigma_0, capped, maxcapture, rho0, usun, u0, vesc, maximum_capture
    ! double precision :: capped_si_spec, capped_sd_spec ! Used in capture_rate_constant()
    character(len=300) :: modfile, filename
    character(len=2) :: spinString(2) = [character(len=2) :: "SI", "SD"]
    character(len=9) :: outfile(7) = [character(len=9) :: "const.dat","qm1.dat","q1.dat","q2.dat","vm1.dat","v1.dat","v2.dat"]
    integer :: nq(7) = [integer :: 0, -1, 1, 2,  0, 0, 0] ! Choose velocity and momentum transfer powers in differential cross-section
    integer :: nv(7) = [integer :: 0,  0, 0, 0, -1, 1, 2]
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables used for the energy transport calculation
    integer :: transport_formalism, nlines
    double precision :: Tx, nwimpsin, noise_indicator, EtransTot
    ! double precision :: evapRate ! Used in fastevap()
    double precision, allocatable :: Etrans(:)
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
    num_isotopes = 29 ! Number of isotopes original capt'n will loop over in the calculation: up to 29 isotopes
    spin_dependency = 0 ! 0=Spin Independent, 1=Spin Dependent

!-----------------------------------------------------------------------------------------------------------------------------------
! Initialise capture calculations
    rho0 = 0.4d0 ! GeV cm^-3
    usun = 235.d0 ! km s^-1
    u0 = 235.d0 ! km s^-1
    vesc = 550.d0 ! km s^-1
    call init_sun(modfile,rho0,usun,u0,vesc)

!-----------------------------------------------------------------------------------------------------------------------------------
! Choose the energy transport formalism: Gould & Raffelt [arxiv:], Rescaled G&R to Monte Carlo [arxiv:], or Spergel & Press [arxiv:]
!       WHICH PAPERS DID THE TRANSPORT FORMALISMS COME FROM?
    transport_formalism = 1 ! 1=G&R, 2=Rescaled G&R, 3=S&P

!-----------------------------------------------------------------------------------------------------------------------------------
! Initialise transport calculations
    call getnlines(nlines)
    allocate(etrans(nlines))
    call get_alpha_kappa(nq,nv)

!-----------------------------------------------------------------------------------------------------------------------------------
! Use the original qv scaling capture calculation
    do j = 1,7
        open(94,file = trim(outfile(j)))
        write(94,*) "Number of Isotopes: ", num_isotopes
        write(94,*) "Spin Dependency: ", spinString(spin_dependency+1)
        write(94,*) "Power: ", outfile(j)
        write(94,*) "Sigma_0 | ", "DM Mass | ", "Captured Dark Matter | ",   "MaxCaptures | ", "Number of WIMPs in | ", "Etranstot"
        do i = 1,10
            mx = 1.d1 ** (dble(i)/5.)
            sigma_0 = 10.d0**(-37.d0)!10d0**(-45+dble(i)/2.)

            print*
            call capture_rate(mx, sigma_0, num_isotopes, nq(j), nv(j), spin_dependency, capped)
            maxcapture = maximum_capture(mx)
            print*, "sigma_0: ", sigma_0, "cm^2 ", &
                    "mx: ", mx, "GeV ", &
                    "Capture rate: ", capped, "s^-1 ", &
                    "Geometric limit: ", maxcapture, "s^-1 "

            ! call capture_rate_constant(mx,sigma_0,sigma_0,capped_sd_spec,capped_si_spec)
            ! print*, "Capture rates (SI, SD): (", capped_si_spec, capped_sd_spec, ") s^-1"

            nwimpsin = 5.d44
            ! nwimpsin = capped*3.d7*4.57d9
            call transgen(sigma_0, nwimpsin, num_isotopes, nq(j), nv(j), spin_dependency, transport_formalism, Tx, &
                            noise_indicator, Etrans, Etranstot)
            print*, "Number of WIMPs in: ", nwimpsin, &
                    "Energy Transport Total: ", EtransTot, "UNITS?" !FIXME units?

            ! call fastevap(sigma_0, 1.d0, 28, EvapRate)
            ! print*, "Evap rate: ", EvapRate, "s^-1"

            write(94,*) sigma_0, mx, capped, maxcapture, nwimpsin, Etranstot
        end do
        close(94)
    end do


!-----------------------------------------------------------------------------------------------------------------------------------
! Specific set up for the Non-Relativistic Effective Operator calculation [arxiv:1501.03729]
    num_isotopes = 16
    jx = 0.5
    couplingVal = 1.65d-8
    call init_nreo()

!-----------------------------------------------------------------------------------------------------------------------------------
! Use the new NREO formalism calculation
    do cpl=1, 14
        filename = "oper_"//trim(cplConsts(cpl))//".dat"
        open(55,file=filename)
        write(55,*) "Coupling Val | ", "DM_Mass | ", "  Captures | ", "  MaxCaptures"

        if (cpl==1) then
            call init_couplings(couplingVal, cpl, 0)
        else if (cpl==2) then
            call init_couplings(0.d0, cpl-1, 0)
            call init_couplings(couplingVal, cpl+1, 0)
        else
            call init_couplings(0.d0, cpl, 0)
            call init_couplings(couplingVal, cpl+1, 0)
        endif

        print*, "Running coupling constant: ", cplConsts(cpl)
        do i = 1,10
            mx = 1.d1 ** (dble(i)/5.)
            call capture_rate_nreo(mx, jx, num_isotopes, capped)
            maxcapture = maximum_capture(mx)
            print*, "Coupling Value: ", couplingVal, "GeV^-4 ", &
                    "DM mass: ", mx, "GeV ", &
                    "Capture rate: ", capped, "s^-1 ", &
                    "Geometric limit: ", maxcapture, "s^-1 "
            write(55,*) couplingVal, mx, capped, maxcapture
        end do
        close(55)
    end do

END PROGRAM GENCAP
