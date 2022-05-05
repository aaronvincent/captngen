! re-running on data from GAMBIT to find the capture rates
program postProcess
    implicit none
    integer :: num_isotopes, cpl, io
    double precision :: jx, v0, vesc, rho0, mx, couplingStrength
    double precision :: capped, maxcap, maxcapture
    character*300 :: filename, modfile

    ! load these as cmd line args?
    filename = "NREOdata/extractedNREO.dat"
    cpl = 4


    modfile = "solarmodels/struct_b16_agss09_reduce10_nohead.dat"
    jx = 0.5
    num_isotopes = 16
    v0 = 0.
    vesc = 0.
    rho0 = 0.
    mx = 0.
    couplingStrength = 0.
    open(unit=9, file=filename, status="old")
    do
        ! load v0 vesc rho0 mx couplingStrength
        read(9, *, IOSTAT=io) v0, vesc, rho0, mx, couplingStrength
        if ( io == 0 ) then
            print*, v0, vesc, rho0, mx, couplingStrength
            
            call captn_init(modfile,rho0,v0,v0,vesc)
            call captn_init_oper()
            
            call populate_array(couplingStrength, cpl, 0)
            call captn_oper(mx,jx,num_isotopes,capped)
            maxcapture = maxcap(mx)
            print*, capped, maxcapture
        else if ( io > 0) then
            print*, io
            stop "Error on reading data from file"
        else
            exit
        end if
    end do
    close(9)

end program postProcess
