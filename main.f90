!   Capt'n General
!program for testing
! the good stuff is in gencap.f90
    PROGRAM GENCAP
    implicit none
    character*300 :: modfile
    double precision :: mx, sigma_0, capped,capped_sd(250),capped_si(250)
    double precision :: capped_si_spec(250),capped_sd_spec(250)
    double precision :: maxcapped, nwimpsin, evapRate(50)
    double precision, allocatable :: Etrans(:)
    double precision :: EtransTot
    integer :: niso, nq, nv, i,nlines
   ! modfile = "solarmodels/struct_b16_agss09_nohead.dat"
   modfile = "solarmodels/struct_b16_agss09_nohead_Tsmoothed.dat" !temeperature smoothed to not nonsense
    ! modfile = "solarmodels/model_gs98_nohead.dat"
    call captn_init(modfile,0.4d0,220.d0,220.d0,600.d0)
    call getnlines(nlines)
    allocate(Etrans(nlines))

    niso = 1
    nq = 0
    nv = 0
    mx = 3.d0


    call get_alpha_kappa(nq,nv)
    call captn_maxcap(mx,maxcapped)
    ! sigma_0 = 1.0d-37 !cm^2
    do i = 1,50
      sigma_0 = 10d0**(-42+dble(i)/5.)

    ! mx = 10**(.02*i - 0.02)
!    print*,mx
    call captn_general(mx,sigma_0,29,nq,nv,capped_si(i))
    print*, "mx: ", mx, "sigma_SI: ", sigma_0, "capped",capped_si(i)
    call captn_general(mx,sigma_0,1,nq,nv,capped_sd(i))
    print*, "mx: ", mx, "sigma_SD: ", sigma_0, "capped",capped_sd(i)

    print*, "maximum captured", maxcapped
    call captn_specific(mx,sigma_0,sigma_0,capped_sd_spec(i),capped_si_spec(i))


    print*,"calling transgen"
    nwimpsin = capped_sd(1)*3.d7*4.57d9
    print*,"passing ", nwimpsin
    call transgen(nwimpsin,1,etrans,Etranstot)
    ! print*,Etrans

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

    END PROGRAM GENCAP
!
