! Capt'n General testing program
!
! Main capture routines can be found in gencap.f90
!
!

    PROGRAM GENCAP
    implicit none
    character*300 :: modfile
    character*10 :: outfile
    double precision :: mx, sigma_0,capped_sd,capped_si
    double precision :: capped_si_spec,capped_sd_spec
    double precision :: maxcap, nwimpsin, evapRate
    double precision, allocatable :: Etrans(:)
    double precision :: EtransTot
    integer :: nq, nv, i, j, nlines

    ! Choose velocity and momentum transfer powers in differential cross-section
    nq = 0
    nv = 0

    ! Choose solar model file
    !modfile = "solarmodels/model_gs98_nohead.dat"
    !modfile = "solarmodels/struct_b16_agss09_nohead.dat"
    modfile = "solarmodels/struct_b16_agss09_nohead_Tsmoothed.dat" !temperature smoothed to not nonsense

    ! Initialise capture calculations
    call captn_init(modfile,0.4d0,220.d0,220.d0,600.d0)

    ! Initialise transport calculations
    call getnlines(nlines)
    allocate(etrans(nlines))
    call get_alpha_kappa(nq,nv)


    ! do j = 1,10
      ! write(outfile,'(a,i2.2,a)') "mx",j,".dat"
      open(55,file = "mx10_sun_const.dat")
      do i = 1,100
        ! mx = 0d0 + dble(i)/5.
        mx = 10.d0
        sigma_0 = 10d0**(-50+dble(i)/5.)
        print*
        print*, "mx: ", mx, "sigma_0:", sigma_0, "cm^2"

        ! print*, "Geometrical limit on capture rate: ", maxcap(mx), "s^-1"

        ! print*,"Calling captn_general for SI scattering."
        ! call captn_general(mx,sigma_0,29,nq,nv,capped_si)
        ! print*, "Capture rate", capped_si, "s^-1"

        ! print*,"Calling captn_general for SD scattering."
        call captn_general(mx,sigma_0,1,nq,nv,capped_sd)
        ! print*, "Capture rate", capped_sd, "s^-1"

        ! print*,"Calling captn_specific for SI and SD scattering."
        ! call captn_specific(mx,sigma_0,sigma_0,capped_sd_spec,capped_si_spec)
        ! print*, "Capture rates (SI, SD): (", capped_si_spec, capped_sd_spec, ") s^-1"

        nwimpsin = 1.2d42
        ! nwimpsin = capped_sd*3.d7*4.57d9
        ! print*,"Calling transgen, with nwimpsin = ", nwimpsin
        call transgen(sigma_0,nwimpsin,1,nq,nv,Etrans,Etranstot)
        ! print*, "Etranstot: ", Etranstot !FIXME units?

        ! print*,"Calling fastevap."
        ! call fastevap(sigma_0,1.d0,28,EvapRate)
        ! print*,"Evap rate: ", EvapRate, "s^-1"

        write(55,*) sigma_0, mx, capped_sd, Etranstot
      end do
      close(55)
    ! end do

    ! Output results to file
    ! open(55,file = "gencap.dat")
    ! do i=1,50
    !  write(55,*) 10d0**(-42+dble(i)/5.), EvapRate(i)
    !  write(55,*) 10**(.02*i - 0.02), capped_sd(i),capped_si(i),capped_sd_spec(i),capped_si_spec(i)
    ! end do
    ! close(55)

    END PROGRAM GENCAP
!
