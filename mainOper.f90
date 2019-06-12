!	Capt'n OPER General
!	program for testing
!	the good stuff is in gencap.f90
!	runs the capture code one isotope at a time, letting the user sum the total later
!	to compile this version of main, switch out 'MAIN = main.o' for 'MAIN = mainOper.o' in the make file (line 8)

PROGRAM GENCAP
	implicit none
	character*300 :: modfile
	double precision :: mx, jx, maxcapped
	double precision :: capped(16,50)
	integer :: niso, i, iso
	character (len=4) :: isotopes(16) = [character(len=4) :: "H","He3","He4","C12","N14","O16","Ne20", "Na23", "Mg24", &
											"Al27", "Si28","S32","Ar40","Ca40","Fe56","Ni58"]

	!modfile = "solarmodels/struct_b16_agss09_nohead.dat"
	modfile = "solarmodels/model_gs98_nohead.dat"
	call captn_init(modfile, 0.4d0, 270.d0, 220.d0, 540.d0)
	call captn_init_oper()
	call populate_array(1.65d-8, 1, 0)
	
	niso = 16
	!mx = 10.d0
	jx = 0.5
	call captn_maxcap(mx, maxcapped)
	print*, "maximum captured", maxcapped
	print*
	do iso=1,niso
		print*, "Isotope: ", isotopes(iso)
		do i = 1,21
			mx = 10**(.1*i + 0.9)
			call captn_oper(mx,jx,niso,iso,capped(iso,i))
			print*, "mx: ",mx, "capped: ",capped(iso,i)
		end do
		print*
	end do

	open(55,file = "captest_oper_c1-0_alliso-gs98.dat")
	do i=1,21
		write(55,*) 10**(.1*i + 0.9), capped(1,i), capped(2,i), capped(3,i), capped(4,i), &
									capped(5,i), capped(6,i), capped(7,i), capped(8,i), &
									capped(9,i), capped(10,i), capped(11,i), capped(12,i), &
									capped(13,i), capped(14,i), capped(15,i), capped(16, i)
	end do
	close(55)
END PROGRAM GENCAP
