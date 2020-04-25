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
	character (len=4) :: isotopes(16) = [character(len=4) :: "H","He3","He4","C12","N14","O16","Ne20", "Na23", "Mg24", &
											"Al27", "Si28","S32","Ar40","Ca40","Fe56","Ni58"]
	character (len=5) :: cplConsts(14) = [character(len=5) :: "c1-0", "c3-0", "c4-0", "c5-0", "c6-0", "c7-0", &
											"c8-0", "c9-0", "c10-0", "c11-0", "c12-0", "c13-0", "c14-0", "c15-0"]

	
	modfile = "solarmodels/Models_From_DarkSUSY/Serenelli-model_ags05_nohead.dat"
	! modfile = "Cut100_Serenelli-model_ags05.dat"
	
	filename = "accuracytesting/testacc_cut1-ags05.dat"
	open(55,file=filename)
	do cpl=1, 14
		write(55,*) cplConsts(cpl)
		call captn_init(modfile, 0.4d0, 220.d0, 220.d0, 540.d0)
		call captn_init_oper()
		if (cpl==1) then
			call populate_array(1.65d-8, cpl, 0)
		else
			call populate_array(1.65d-8, cpl+1, 0)
		endif
		
		niso = 17
		jx = 0.5
		call captn_maxcap(mx, maxcapped)
		print*, "Running coupling constant: ", cplConsts(cpl)
		print*, "maximum captured", maxcapped
		print*

		mx = 5
		call captn_oper(mx,jx,niso,iso,capped)
		print*, "mx: ",mx, "GeV, capped: ",capped
		write(55,*) mx, capped

		mx = 50
		call captn_oper(mx,jx,niso,iso,capped)
		print*, "mx: ",mx, "GeV, capped: ",capped
		write(55,*) mx, capped
		write(55,*)
	end do
	close(55)
END PROGRAM GENCAP
