!	Capt'n OPER General
!	program for testing
!	the good stuff is in gencap.f90
!	runs the capture code one isotope at a time, letting the user sum the total later
!	to compile this version of main, switch out 'MAIN = main.o' for 'MAIN = mainOper.o' in the make file (line 8)

PROGRAM GENCAP
	implicit none
	character*300 :: modfile, filename
	double precision :: mx, jx, maxcapped
	double precision :: capped
	integer :: niso, i, iso, cpl
	character (len=5) :: cplConsts(14) = [character(len=5) :: "c1-0", "c3-0", "c4-0", "c5-0", "c6-0", "c7-0", &
											"c8-0", "c9-0", "c10-0", "c11-0", "c12-0", "c13-0", "c14-0", "c15-0"]

	
	! modfile = "Aarons_Model_Cuts/struct_b16_agss09_nohead.dat"
	modfile = "solarmodels/struct_b16_agss09_reduce10_nohead.dat"
	niso = 16
	jx = 0.5
	iso = 0 !force captnoper to run and sum all isotopes

	call captn_init(modfile, 0.4d0, 220.d0, 220.d0, 540.d0)
	call captn_init_oper()
	filename = "testfile_oper"
	open(55,file=filename)
	do cpl=1, 14
		write(55,*) cplConsts(cpl)
		if (cpl==1) then
			call populate_array(1.65d-8, cpl, 0)
		else
			call populate_array(1.65d-8, cpl+1, 0)
		endif
		
		print*, "Running coupling constant: ", cplConsts(cpl)
		print*

		mx = 5
		call captn_maxcap(mx, maxcapped)
		call captn_oper(mx,jx,niso,iso,capped)
		print*, "mx: ",mx, "maximum captured", maxcapped
		print*, "GeV, capped: ",capped
		write(55,*) mx, maxcapped, capped

		mx = 50
		call captn_maxcap(mx, maxcapped)
		call captn_oper(mx,jx,niso,iso,capped)
		print*, "mx: ",mx, "maximum captured", maxcapped
		print*, "GeV, capped: ",capped
		write(55,*) mx, maxcapped, capped

		mx = 500
		call captn_maxcap(mx, maxcapped)
		call captn_oper(mx,jx,niso,iso,capped)
		print*, "mx: ",mx, "maximum captured", maxcapped
		print*, "GeV, capped: ",capped
		write(55,*) mx, maxcapped, capped
		write(55,*)
		
		if (cpl==1) then
			call populate_array(0, cpl, 0)
		else
			call populate_array(0, cpl+1, 0)
		endif
	end do
	close(55)
END PROGRAM GENCAP
