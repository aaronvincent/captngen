! Capt'n General testing program
!
! Main capture routines can be found in gencap.f90
!
!

    PROGRAM GENCAP
    implicit none
	integer :: transport_formalism
    character*300 :: modfile, filename
    character*100 :: outfile(7)
    double precision :: mx, Tx, jx, sigma_0,capped_sd,capped_si, maxcapture
    double precision :: capped_si_spec,capped_sd_spec
    double precision :: maximum_capture, nwimpsin, evapRate, noise_indicator
    double precision, allocatable :: Etrans(:)
    double precision :: EtransTot
    integer :: nq(7), nv(7), i, j, k, nlines, num_isotopes, spin_dependency, cpl
    character (len=5) :: cplConsts(14) = [character(len=5) :: "c1-0", "c3-0", "c4-0", "c5-0", "c6-0", "c7-0", &
                        "c8-0", "c9-0", "c10-0", "c11-0", "c12-0", "c13-0", "c14-0", "c15-0"]

	transport_formalism = 1 ! 1=Gould & Raffelt, 2=G&R rescaled to agree with MC, 3=Spergel & Press
	
    ! Choose velocity and momentum transfer powers in differential cross-section
    nq = [0,-1,1,2,0,0,0]
    nv = [0,0,0,0,-1,1,2]

    outfile = ['const.dat','qm1--.dat','q1---.dat','q2---.dat','vm1--.dat','v1---.dat','v2---.dat']

    ! Choose solar model file
    !modfile = "solarmodels/model_gs98_nohead.dat"
    !modfile = "solarmodels/struct_b16_agss09_nohead.dat"
    ! modfile = "solarmodels/struct_b16_agss09_nohead_Tsmoothed.dat" !temperature smoothed to not nonsense
    modfile = "solarmodels/struct_b16_agss09_reduce10_nohead.dat"

    ! number of isotopes capt'n will loop over in the calculation: up to 29 isotopes
    num_isotopes = 1

    ! zero for spin_independent, one for spin_dependent
    spin_dependency = 1

    ! Initialise capture calculations
    call init_sun(modfile,0.4d0,235.d0,235.d0,550.d0)

    ! ! Initialise transport calculations
    ! call getnlines(nlines)
    ! allocate(etrans(nlines))
    ! call get_alpha_kappa(nq,nv)
    

    ! do j = 1,1
    !   open(94,file = outfile(j))
    !   write(94,*) "Number of Isotopes: ", num_isotopes
    !   write(94,*) "Spin Dependency: ", spin_dependency
    !   write(94,*) "Power: ", outfile(j)
    !   write(94,*) "Sigma_0 | ", "DM Mass | ", "Capptured Dark Matter | ", "Etranstot"
    !   do i = 1,1
    !     mx = 5.d0 !+ dble(i)/5.
    !     sigma_0 = 10.d0**(-37.d0)!10d0**(-45+dble(i)/2.)
    !     print*
    !     print*, "mx: ", mx, "sigma_0:", sigma_0, "cm^2"

    !     ! print*, "Geometrical limit on capture rate: ", maximum_capture(mx), "s^-1"

    !     ! print*,"Calling capture_rate for SI scattering."
    !     ! call capture_rate(mx,sigma_0,29,nq,nv,capped_si)
    !     ! print*, "Capture rate", capped_si, "s^-1"

    !     ! print*,"Calling capture_rate for SD scattering."
    !     call capture_rate(mx,sigma_0,num_isotopes,nq(j),nv(j),spin_dependency,capped_sd)
    !     ! print*, "Capture rate", capped_sd, "s^-1"

    !     ! print*,"Calling captn_specific for SI and SD scattering."
    !     ! call captn_specific(mx,sigma_0,sigma_0,capped_sd_spec,capped_si_spec)
    !     ! print*, "Capture rates (SI, SD): (", capped_si_spec, capped_sd_spec, ") s^-1"

    !     nwimpsin = 5.d44
    !     ! nwimpsin = capped_sd*3.d7*4.57d9
    !     ! print*,"Calling transgen, with nwimpsin = ", nwimpsin
    !     call transgen(sigma_0,nwimpsin,num_isotopes,nq(j),nv(j),spin_dependency,transport_formalism, &
    !     	Tx,noise_indicator,Etrans,Etranstot)
    !     ! print*, "Etranstot: ", Etranstot !FIXME units?
    !     ! print*,"Calling fastevap."
    !     ! call fastevap(sigma_0,1.d0,28,EvapRate)
    !     ! print*,"Evap rate: ", EvapRate, "s^-1"

    !     write(94,*) sigma_0, mx, capped_sd, Etranstot
    !   end do
    !   close(94)
    ! end do
    
   call captn_init_oper()
   num_isotopes = 16
   jx = 0.5
   do cpl=1, 14
     filename = "oper_"//trim(cplConsts(cpl))//"_model-OMPedit-parallel.dat"
     open(55,file=filename)
     write(55,*) "DM_Mass | ", "  Captures | ", "  MaxCaptures"

     if (cpl==1) then
       call populate_array(1.65d-8, cpl, 0)
     else if (cpl==2) then
       call populate_array(0.d0, cpl-1, 0)
       call populate_array(1.65d-8, cpl+1, 0)
     else
       call populate_array(0.d0, cpl, 0)
       call populate_array(1.65d-8, cpl+1, 0)
     endif
     
     print*, "Running coupling constant: ", cplConsts(cpl)
     do i = 1,1
       mx = 10. !5.d0 + dble(i)/5.
       call captn_oper(mx,jx,num_isotopes,capped_sd)
       maxcapture = maximum_capture(mx)
       write(55,*) mx, capped_sd, maxcapture
       print*, "mx: ",mx, "capped: ",capped_sd, "max_capture:",maxcapture
     end do
     close(55)
   end do

END PROGRAM GENCAP
