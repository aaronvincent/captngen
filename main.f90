! Capt'n General testing program
!
! Main capture routines can be found in gencap.f90
!
!

    PROGRAM GENCAP
    implicit none
    character*300 :: modfile
    character*100 :: outfile(7)
    double precision :: mx, sigma_0,capped_sd,capped_si
    double precision :: capped_si_spec,capped_sd_spec
    double precision :: maxcap, nwimpsin, evapRate
    double precision, allocatable :: Etrans(:)
    double precision :: EtransTot
    integer :: nq(7), nv(7), i, j, nlines

    ! Choose velocity and momentum transfer powers in differential cross-section
    nq = [0,-1,1,2,0,0,0]
    nv = [0,0,0,0,-1,1,2]

    outfile = ['const.dat','qm1--.dat','q1---.dat','q2---.dat','vm1--.dat','v1---.dat','v2---.dat']

    ! Choose solar model file
    !modfile = "solarmodels/model_gs98_nohead.dat"
    !modfile = "solarmodels/struct_b16_agss09_nohead.dat"
    modfile = "solarmodels/struct_b16_agss09_nohead_Tsmoothed.dat" !temperature smoothed to not nonsense

    ! Initialise capture calculations
    call captn_init(modfile,0.4d0,235.d0,235.d0,550.d0)

    ! Initialise transport calculations
    call getnlines(nlines)
    allocate(etrans(nlines))
    call get_alpha_kappa(nq,nv)


    do j = 1,7
      open(94,file = outfile(j))
      do i = 1,100
        mx = 0d0 + dble(i)/5.
        sigma_0 = 10d0**(-50+dble(i)/5.)
        print*
        print*, "mx: ", mx, "sigma_0:", sigma_0, "cm^2"

        ! print*, "Geometrical limit on capture rate: ", maxcap(mx), "s^-1"

        ! print*,"Calling captn_general for SI scattering."
        ! call captn_general(mx,sigma_0,29,nq,nv,capped_si)
        ! print*, "Capture rate", capped_si, "s^-1"

        ! print*,"Calling captn_general for SD scattering."
        call captn_general(mx,sigma_0,1,nq(j),nv(j),1,capped_sd)
        ! print*, "Capture rate", capped_sd, "s^-1"

        ! print*,"Calling captn_specific for SI and SD scattering."
        ! call captn_specific(mx,sigma_0,sigma_0,capped_sd_spec,capped_si_spec)
        ! print*, "Capture rates (SI, SD): (", capped_si_spec, capped_sd_spec, ") s^-1"

        nwimpsin = 1.2d42
        ! nwimpsin = capped_sd*3.d7*4.57d9
        ! print*,"Calling transgen, with nwimpsin = ", nwimpsin
        call transgen(sigma_0,nwimpsin,1,nq(j),nv(j),1,Etrans,Etranstot)
        ! print*, "Etranstot: ", Etranstot !FIXME units?

        ! print*,"Calling fastevap."
        ! call fastevap(sigma_0,1.d0,28,EvapRate)
        ! print*,"Evap rate: ", EvapRate, "s^-1"

        write(94,*) sigma_0, mx, capped_sd, Etranstot
      end do
      close(94)
    end do

    END PROGRAM GENCAP
!
