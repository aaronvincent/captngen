! Capt'n General testing program
!
! Main capture routines can be found in gencap.f90
!
!

	module trapzmod
	contains
	
	function trapz(x,y,flen)
      implicit none
      integer, intent(in) :: flen
      double precision, intent (in) :: x(flen), y(flen)
      double precision trapz

      integer i


      trapz = y(1)*(x(2)-x(1))/2. + y(flen)*(x(flen)-x(flen-1))/2.
      do i = 2,flen-1
        trapz = trapz + y(i)*(x(i)-x(i-1))

!        if (trapz .lt. 0.d0) then
!          print*, "negative encountered in trapz: i = ", i
!        end if
      end do


      return
      end function
	
	end module

    PROGRAM GENCAP
    use trapzmod
    implicit none
    character*300 :: modfile
    double precision :: mx, sigma_0,capped_sd(250),capped_si(250)
    double precision :: capped_si_spec(250),capped_sd_spec(250), AtomicNumber(29)
    double precision :: maxcap, nwimpsin, evapRate(50), Tx, Nbar
    double precision, allocatable :: Etrans(:), Etrans_all(:,:), msum(:)
    double precision :: EtransTot
    integer :: nq, nv, i, nlines
    logical :: nonlocal
    double precision, allocatable :: tab_mencl(:), tab_r(:), tab_T(:), tab_starrho(:), tab_mfr(:,:)
    double precision :: Pres, Lumi
    double precision, parameter :: mnucg = 1.6726219d-24, pi=3.1415626d0, Rsun=6.9634d10

    ! Choose velocity and momentum transfer powers in differential cross-section
    nq = 0
    nv = 0
	nonlocal = .false.

    ! Choose solar model file
    !modfile = "solarmodels/model_gs98_nohead.dat"
    !modfile = "solarmodels/struct_b16_agss09_nohead.dat"
    modfile = "solarmodels/struct_b16_agss09_nohead_Tsmoothed.dat" !temperature smoothed to not nonsense

    ! Initialise capture calculations
    call captn_init(modfile,1.d3,220.d0,220.d0,600.d0)

    ! Initialise transport calculations
    call getnlines(nlines)
    allocate(etrans(nlines))
    allocate(Etrans_all(nlines,50))
    call get_alpha_kappa(nq,nv)
    
    allocate(tab_mencl(nlines))
    allocate(tab_r(nlines))
    allocate(tab_T(nlines))
    allocate(tab_starrho(nlines))
    allocate(tab_mfr(nlines,29)) !we could just allocate niso, but this leads to problems
    allocate(msum(nlines))

    !now actually read in the file
    open(99,file=modfile)
    do i=1,nlines
    	read(99,*) tab_mencl(i),tab_r(i), tab_T(i), tab_starrho(i), Pres, Lumi, tab_mfr(i,:)
    end do
    close(99)
    
     AtomicNumber  = (/ 1., 4., 3., 12., 13., 14., 15., 16., 17., &
                      18., 20.2, 22.99, 24.3, 26.97, 28.1, 30.97,32.06, 35.45, &
                      39.948, 39.098, 40.08, 44.95, 47.86, 50.94, 51.99, &
                      54.93, 55.845, 58.933, 58.693/)
    
    ! A useful sum:
    do i=1,29
    	msum = msum + tab_mfr(:,i)/AtomicNumber(i)
    enddo
    print *, "max(msum)=", maxval(abs(msum))
    
    ! Total number of baryons:
    Nbar = trapz(tab_r*Rsun, msum*tab_starrho/mnucg*4*pi*(tab_r*Rsun)**2, nlines)

    do i = 1,50
      
      mx = 10.d0 !2.d0+dble(i)
      sigma_0 = 10d0**(-40+dble(i)/5.d0)
      print*
      print *, "i=", i
      print*, "mx: ", mx, "sigma_0:", sigma_0, "cm^2"

      print*, "Geometrical limit on capture rate: ", maxcap(mx), "s^-1"

      print*,"Calling captn_general for SI scattering."
      call captn_general(mx,sigma_0,29,nq,nv,capped_si(i))
      print*, "Capture rate", capped_si(i), "s^-1"

      print*,"Calling captn_general for SD scattering."
      call captn_general(mx,sigma_0,1,nq,nv,capped_sd(i))
      print*, "Capture rate", capped_sd(i), "s^-1"

      print*,"Calling captn_specific for SI and SD scattering."
      call captn_specific(mx,sigma_0,sigma_0,capped_sd_spec(i),capped_si_spec(i))
      print*, "Capture rates (SI, SD): (", capped_si_spec(i), capped_sd_spec(i), ") s^-1"

!      nwimpsin = capped_sd(i)*3.d7*4.57d9
	  nwimpsin = 1.d-15*Nbar
	  print*, "Total number of baryons: Nbar=", Nbar
      print*,"Calling transgen, with nwimpsin = ", nwimpsin
      call transgen(sigma_0,nwimpsin,1,nonlocal,Tx,Etrans,Etranstot)
      print*, "Etranstot: ", Etranstot !FIXME units?
      Etrans_all(:,i) = Etrans

      print*,"Calling fastevap."
      call fastevap(sigma_0,1.d0,28,EvapRate(i))
      print*,"Evap rate: ", EvapRate(i), "s^-1"
      
      if (nonlocal) then
      	open(55, file="/home/luke/summer_2020/mesa/test_files/Ltot_sp.dat", access="APPEND")
      else if (.not. nonlocal) then
      	open(55, file="/home/luke/summer_2020/mesa/test_files/Ltot_gr.dat", access="APPEND")
      endif
      
      write(55,*) sigma_0, Etranstot
      close(55)
      
    end do

    !Output results to file
    if (nonlocal) then 
    	open(55,file = "/home/luke/summer_2020/mesa/test_files/gentest_sp.dat")
    else if (.not. nonlocal) then 
    	open(55,file = "/home/luke/summer_2020/mesa/test_files/gentest_gr.dat")
    endif
	do i=1,nlines
		write(55,*) tab_r(i), Etrans_all(i,:)
	enddo
!    do i=1,50
!      write(55,*) 10d0**(-42+dble(i)/5.), EvapRate(i)
!      write(55,*) 10**(.02*i - 0.02), capped_sd(i),capped_si(i),capped_sd_spec(i),capped_si_spec(i)
!    end do
    close(55)

    END PROGRAM GENCAP
