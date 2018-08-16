!	Capt'n OPER General
!	program for testing
!	the good stuff is in gencap.f90
!   runs the capture code with all isotopes at once, summing them

PROGRAM GENCAP
    implicit none
    character*300 :: modfile
    double precision :: mx, jx, maxcapped
    double precision :: capped(250)
    integer :: niso, i

	modfile = "solarmodels/struct_b16_agss09_nohead.dat"
    !modfile = "solarmodels/model_gs98_nohead.dat"
    call captn_init(modfile, 0.4d0, 270.d0, 220.d0, 540.d0)
    call captn_init_oper()

    niso = 16
    mx = 10.d0
	jx = 0.5
    call captn_maxcap(mx, maxcapped)
	print*, "maximum captured", maxcapped
	print*
    do i = 1,21
    	mx = 10**(.1*i + 0.9)
    	call captn_oper(mx,jx,niso,0,capped(i))
		print*, "mx: ", mx, "capped:",capped(i)
    end do

    open(55,file = "captest_oper_c1-0_test.dat")
    do i=1,21
    	write(55,*) 10**(.1*i + 0.9), capped(i)
    end do
    close(55)
END PROGRAM GENCAP
