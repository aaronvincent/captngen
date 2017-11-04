!   Capt'n General
!program for testing
! the good stuff is in gencap.f90
    PROGRAM GENCAP
    implicit none
    character*300 :: modfile
    double precision :: mx, sigma_0, capped,capped_sd(250),capped_si(250)
    double precision :: capped_si_spec(250),capped_sd_spec(250)
    double precision :: maxcapped
    integer :: niso, nq, nv, i
   modfile = "solarmodels/struct_b16_agss09_nohead.dat"
    ! modfile = "solarmodels/model_gs98_nohead.dat"
    call captn_init(modfile,0.4d0,270.d0,220.d0,540.d0)

    niso = 29
    nq = -1
    nv = 0
    mx = 10.d0
    call captn_maxcap(mx,maxcapped)
    sigma_0 = 1.0d-40 !cm^2
    ! do i = 1,50
    i = 1
    mx = 10**(.02*i - 0.02)
!    print*,mx
    call captn_general(mx,sigma_0,29,nq,nv,capped_si(i))
    print*, "mx: ", mx, "sigma_SI: ", sigma_0, "capped",capped_si(i)
    call captn_general(mx,sigma_0,1,nq,nv,capped_sd(i))
    print*, "mx: ", mx, "sigma_SD: ", sigma_0, "capped",capped_sd(i)

    print*, "maximum captured", maxcapped
    call captn_specific(mx,sigma_0,sigma_0,capped_sd_spec(i),capped_si_spec(i))
    ! end do

    ! open(55,file = "captest_agss_q2.dat")
    ! do i=1,250
    ! write(55,*) 10**(.02*i - 0.02), capped_sd(i),capped_si(i),capped_sd_spec(i),capped_si_spec(i)
    ! end do
    ! close(55)

    END PROGRAM GENCAP
!
