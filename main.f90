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
    double precision :: maxcap, nwimpsin, evapRate(10), sigma(10)
    double precision, allocatable :: Etrans(:), Etrans_all(:,:), tab_mencl(:), r(:), tab_T(:), tab_starrho(:), &
    	tab_mfr(:,:)
    double precision :: EtransTot, Lumi, Pres
    integer :: nq, nv, i, nlines
    logical :: nonlocal

    ! Choose velocity and momentum transfer powers in differential cross-section
    nq = 0
    nv = 0
    nonlocal = .true.

    ! Choose solar model file
    !modfile = "solarmodels/model_gs98_nohead.dat"
    !modfile = "solarmodels/struct_b16_agss09_nohead.dat"
    modfile = "solarmodels/struct_b16_agss09_nohead_Tsmoothed.dat" !temperature smoothed to not nonsense

    ! Initialise capture calculations
    call captn_init(modfile,4.d-1,220.d0,220.d0,600.d0)

    ! Initialise transport calculations
    call getnlines(nlines)
    allocate(etrans(nlines))
    allocate(Etrans_all(nlines, 10))
    allocate(r(nlines))
    allocate(tab_mencl(nlines))
    allocate(tab_T(nlines))
    allocate(tab_starrho(nlines))
    allocate(tab_mfr(nlines, 1))
    call get_alpha_kappa(nq,nv)
    
    ! Read radius from solar model file. Other variables are dummies; just needed to read the file
    open(99,file=modfile)
    do i=1,nlines
    read(99,*) tab_mencl(i), r(i), tab_T(i), tab_starrho(i), Pres, Lumi, tab_mfr(i,:)
    end do
    close(99)

    do i = 1,10

      ! mx = 10**(.19*dble(i))
      mx = 10.d0
      sigma_0 = 10d0**(-42+dble(i))
      sigma(i) = sigma_0
      print*
      print*, "mx: ", mx, "sigma_0:", sigma_0, "cm^2"

!      print*, "Geometrical limit on capture rate: ", maxcap(mx), "s^-1"

      print*,"Calling captn_general for SI scattering."
      call captn_general(mx,sigma_0,29,nq,nv,capped_si(i))
      print*, "Capture rate", capped_si(i), "s^-1"

      print*,"Calling captn_general for SD scattering."
      call captn_general(mx,sigma_0,1,nq,nv,capped_sd(i))
      print*, "Capture rate", capped_sd(i), "s^-1"

      print*,"Calling captn_specific for SI and SD scattering."
      call captn_specific(mx,sigma_0,sigma_0,capped_sd_spec(i),capped_si_spec(i))
      print*, "Capture rates (SI, SD): (", capped_si_spec(i), capped_sd_spec(i), ") s^-1"

!      nwimpsin = capped_sd(i)*3.15d7*4.57d9
	  nwimpsin = 5.d47
      print*,"Calling transgen, with nwimpsin = ", nwimpsin
      call transgen(nwimpsin,1,nonlocal,Etrans,Etranstot)
      print*, "Etranstot: ", Etranstot !FIXME units?
      Etrans_all(:,i) = Etrans(:)

      print*,"Calling fastevap."
      call fastevap(nwimpsin, 1, EvapRate(i))
      print*,"Evap rate: ", EvapRate(i), "s^-1"

    end do

    open(55, file="gentest_m10.dat")
    do i=1,nlines
    	write(55,*) r(i), Etrans_all(i,:)
    enddo
    close(55)

    END PROGRAM GENCAP
!
