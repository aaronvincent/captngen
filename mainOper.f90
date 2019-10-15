!	Capt'n OPER General
!	program for testing
!	the good stuff is in gencap.f90
!	runs the capture code one isotope at a time, letting the user sum the total later
!	to compile this version of main, switch out 'MAIN = main.o' for 'MAIN = mainOper.o' in the make file (line 8)

PROGRAM GENCAP
	implicit none
	character*300 :: modfile, filename
	double precision :: mx, jx, maxcapped
	double precision :: capped(16,50)
	integer :: niso, i, iso, cpl
	character (len=4) :: isotopes(16) = [character(len=4) :: "H","He3","He4","C12","N14","O16","Ne20", "Na23", "Mg24", &
											"Al27", "Si28","S32","Ar40","Ca40","Fe56","Ni58"]
	character (len=5) :: cplConsts(14) = [character(len=5) :: "c1-0", "c3-0", "c4-0", "c5-0", "c6-0", "c7-0", &
											"c8-0", "c9-0", "c10-0", "c11-0", "c12-0", "c13-0", "c14-0", "c15-0"]

	! modfile = "solarmodels/struct_b16_agss09_nohead.dat"
	! modfile = "solarmodels/model_gs98_nohead.dat"

	! from DarkSUSY:
	modfile = "solarmodels/Models_From_DarkSUSY/Serenelli-model_ags05_nohead.dat"
	! modfile = "solarmodels/Models_From_DarkSUSY/Serenelli-model_agss09_nohead.dat"
	! modfile = "solarmodels/Models_From_DarkSUSY/Serenelli-model_agss09ph_nohead.dat"
	! modfile = "solarmodels/Models_From_DarkSUSY/Serenelli-model_gs98_nohead.dat"
	
	do cpl=1, 14
		call captn_init(modfile, 0.4d0, 220.d0, 220.d0, 540.d0)
		call captn_init_oper()
		if (cpl==1) then
			call populate_array(1.65d-8, cpl, 0)
		else
			call populate_array(1.65d-8, cpl+1, 0)
		endif
		
		niso = 16
		!mx = 10.d0
		jx = 0.5
		call captn_maxcap(mx, maxcapped)
		print*, "Running coupling constant: ", cplConsts(cpl)
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
		filename = "testingCatena/Oper_Usun_update/captest_oper_"//trim(cplConsts(cpl))//"_alliso-ags05.dat"
		open(55,file=filename)!"testingCatena/Oper_data/captest_oper_c1-0_alliso-gs98.dat")
		do i=1,21
			write(55,*) 10**(.1*i + 0.9), capped(1,i), capped(2,i), capped(3,i), capped(4,i), &
										capped(5,i), capped(6,i), capped(7,i), capped(8,i), &
										capped(9,i), capped(10,i), capped(11,i), capped(12,i), &
										capped(13,i), capped(14,i), capped(15,i), capped(16, i)
		end do
		close(55)
		print*
	end do
END PROGRAM GENCAP
