!   Capt'n General
!program for testing
! the good stuff is in gencap.f90
    PROGRAM GENCAP
    implicit none
    character*300 :: modfile
    double precision :: mx, sigma_0, capped,capped_sd(250),capped_si(250)
    double precision :: capped_si_spec(250),capped_sd_spec(250)
    double precision :: maxcapped, nwimpsin, evapRate(50)
    double precision, allocatable :: Etrans(:), Etrans_all(:,:), tab_mencl(:), tab_r(:), tab_starrho(:), &
    	tab_mfr(:,:), tab_vesc(:), phi(:), tab_dr(:), tab_T(:), tab_g(:)
    double precision :: Pres, Lumi
    double precision :: EtransTot
    integer :: niso, nq, nv, i,nlines
    logical :: nonlocal
   ! modfile = "solarmodels/struct_b16_agss09_nohead.dat"
   modfile = "solarmodels/struct_b16_agss09_nohead_Tsmoothed.dat" !temeperature smoothed to not nonsense
    ! modfile = "solarmodels/model_gs98_nohead.dat"
    call captn_init(modfile,1.d3,220.d0,220.d0,600.d0)
    call getnlines(nlines)
    
    allocate(Etrans(nlines))
    allocate(Etrans_all(nlines,50))
	allocate(tab_mencl(nlines))
    allocate(tab_r(nlines))
    allocate(tab_starrho(nlines))
    allocate(tab_mfr(nlines,29)) !we could just allocate niso, but this leads to problems
    allocate(tab_vesc(nlines))
    allocate(phi(nlines))
    allocate(tab_dr(nlines))
    allocate(tab_T(nlines)) !not used in capgen; used for transgen (and anngen? )
    allocate(tab_g(nlines))


    !now actually read in the file
    open(99,file=modfile)
    do i=1,nlines
    read(99,*) tab_mencl(i),tab_r(i), tab_T(i), tab_starrho(i), Pres, Lumi, tab_mfr(i,:)
    end do
    close(99)

    niso = 1
    nq = 0
    nv = 0
    mx = 10.d0


    call get_alpha_kappa(nq,nv)
    call captn_maxcap(mx,maxcapped)
     sigma_0 = 1.0d-37 !cm^2
    do i = 1,50
    sigma_0 = 10.d0**(-42.+dble(i)/5.)
!	sigma_0 = 1.0d-36

    ! mx = 10**(.02*i - 0.02)
!    print*,mx
    call captn_general(mx,sigma_0,29,nq,nv,capped_si(i))
    print*, "mx: ", mx, "sigma_SI: ", sigma_0, "capped",capped_si(i)
    call captn_general(mx,sigma_0,1,nq,nv,capped_sd(i))
    print*, "mx: ", mx, "sigma_SD: ", sigma_0, "capped",capped_sd(i)

    print*, "maximum captured", maxcapped
    call captn_specific(mx,sigma_0,sigma_0,capped_sd_spec(i),capped_si_spec(i))


	nonlocal = .true.

    print*,"calling transgen"
    nwimpsin = capped_sd(1)*3.d7*4.57d9
    print*,"passing ", nwimpsin
    print *, "nw = ", nwimpsin, "etrans = ", etrans(10), "Etranstot = ", Etranstot
    call transgen(nwimpsin,1,nonlocal,etrans,Etranstot)
    print*, "Etrans = ", Etrans(10)
    Etrans_all(:, i) = Etrans

    call fastevap(1.d0,28,EvapRate(i))
    print*,"Evap rate: ", EvapRate(i)

    end do

    open(55,file = "evaptest.dat")
    do i=1,50
      write(55,*) 10d0**(-42+dble(i)/5.), EvapRate(i)
    end do
    close(55)
    ! write(55,*) 10**(.02*i - 0.02), capped_sd(i),capped_si(i),capped_sd_spec(i),capped_si_spec(i)
    ! end do
    ! close(55)
    if (nonlocal) then
    	open(56, file="etrans_test_sp.dat")
		do i=1,nlines
			write(56,*) tab_r(i), Etrans_all(i,:)
		enddo
		close(56)
	else if (nonlocal .eqv. .false.) then
		open(56, file="etrans_test_gr.dat")
		do i=1,nlines
			write(56,*) tab_r(i), Etrans_all(i,:)
		enddo
		close(56)
	endif
    END PROGRAM GENCAP
!
