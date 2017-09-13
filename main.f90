!   Capt'n General
!program for testing
! the good stuff is in gencap.f90
    PROGRAM GENCAP
    implicit none
    double precision :: mx, sigma_0, capped,capped_sd(250),capped_si(250)
    integer :: niso, nq, nv, i

    niso = 29
    nq = 1
    nv = 0
    mx = 1000.d0
    sigma_0 = 1.0d-40 !cm^2
    do i = 1,250
    mx = 10**(.02*i - 0.02)
!    print*,mx
    call captn_general(mx,sigma_0,29,nq,nv,capped_si(i))
    call captn_general(mx,sigma_0,1,nq,nv,capped_sd(i))

    end do

    open(55,file = "captest_agss_q2.dat")
    do i=1,250
    write(55,*) 10**(.02*i - 0.02), capped_sd(i),capped_si(i)
    end do
    close(55)

    END PROGRAM GENCAP
!
