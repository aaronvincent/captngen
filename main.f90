! Capt'n General testing program
!
! Main capture routines can be found in gencap.f90
!
!

    PROGRAM GENCAP
    implicit none
    character*300 :: modfile
    double precision :: mx, sigma_0,capped_sd(250),capped_si(250)
    double precision :: capped_si_spec(250),capped_sd_spec(250)
    double precision :: maxcap, nwimpsin, evapRate(50)
    double precision, allocatable :: Etrans(:)
    double precision :: EtransTot
    integer :: nq, nv, i, nlines

    ! Choose velocity and momentum transfer powers in differential cross-section
    nq = 0
    nv = 0

    ! Choose solar model file
    !modfile = "solarmodels/model_gs98_nohead.dat"
    !modfile = "solarmodels/struct_b16_agss09_nohead.dat"
    modfile = "solarmodels/struct_b16_agss09_nohead_Tsmoothed.dat" !temperature smoothed to not nonsense

    ! Initialise capture calculations
    call captn_init(modfile,1d3,220.d0,220.d0,600.d0)

    ! Initialise transport calculations
    call getnlines(nlines)
    allocate(etrans(nlines))
    call get_alpha_kappa(nq,nv)

    ! do i = 1,10

      ! mx = 10**(.19*dble(i))
      mx = 10.d0
      sigma_0 = 1d-37 !10d0**(-42+dble(i)/5.)
      i = 1
      print*
      print*, "mx: ", mx, "sigma_0:", sigma_0, "cm^2"

      print*, "Geometrical limit on capture rate: ", maxcap(mx), "s^-1"

      print*,"Calling captn_general for SI scattering."
      call captn_general(mx,sigma_0,29,nq,nv,capped_si(i))
      print*, "Capture rate", capped_si(i), "s^-1"

      print*,"Calling captn_general for SD scattering."
      call captn_general(mx,sigma_0,1,nq,nv,capped_sd(i))
      print*, "Capture rate", capped_sd(i), "s^-1"

      print*,"Calling captn_specific for SI and SD scattering."
      call captn_specific(mx,sigma_0,sigma_0,capped_sd_spec(i),capped_si_spec(i))
      print*, "Capture rates (SI, SD): (", capped_si_spec(i), capped_sd_spec(i), ") s^-1"

      nwimpsin = capped_sd(i)*3.15d7*4.57d9
      print*,"Calling transgen, with nwimpsin = ", nwimpsin
      call transgen(nwimpsin,sigma_0,1,Etrans,Etranstot)
      print*, "Etranstot: ", Etranstot !FIXME units?

      print*,"Calling fastevap."
      call fastevap(sigma_0,1.d0,28,EvapRate(i))
      print*,"Evap rate: ", EvapRate(i), "s^-1"

    ! end do

    !Output results to file
    !open(55,file = "gencap.dat")
    !do i=1,50
    !  write(55,*) 10d0**(-42+dble(i)/5.), EvapRate(i)
    !  write(55,*) 10**(.02*i - 0.02), capped_sd(i),capped_si(i),capped_sd_spec(i),capped_si_spec(i)
    !end do
    !close(55)

    END PROGRAM GENCAP
!
