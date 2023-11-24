!   Capt'n General
!   For a continuum of Q-dependent capture
!   Simplified, general solar DM capture routine
!   Standalone code for q**2n, v**2n
!   Useful stuff is run at the end; beginning is the module that does the heavy lifting
!   Future plans: add form factor handling (a la Catena & Schwabe)
!   Made for GAMBIT, with marginal competence
!   Aaron Vincent 2017
!   all units of distance: cm
!   all units of mass/energy : GeV (or GeV/c**2, don't forget)
!   all units of time: seconds
!   Sticking with notation of 1504.04378. Cite that paper. Or 1605.06502 it's even better.
!   Reference q0 is 40 MeV, and v0 is 220 km/s.

!   Updated 2020 just to handle all nq,nv=-1,0,1,2 cases (integration limits were causing issues).
!   Working for spin-dependent interactions with atomic hydrogen (niso=1)
!   NOTE: removed evaporation calcs - for some reason fastevap was still being used even when
!     option was turned off from DarkMESA side.

    module capmod

      use sharedmod
      implicit none
      double precision, parameter :: GNewt = 6.672d-8, kb_here = 1.380649d-23
      double precision, parameter :: v0 = 220.d5
      double precision :: q0, ue_at_r, ue, uChi
      double precision ::  w_shared !to make accesible and modifiable in 2D integration
      double precision :: cutoffIntegrand
      ! nq and nv can be -1, 0, 1, 2; this is set in the main program
      integer :: nq, nv

        contains

      !generalized form factor: hydrogen
      function GFFI_H(w,vesc)
        double precision :: p, w,vesc,u,GFFI_H,G
        p = mdm*w
        u = sqrt(w**2-vesc**2)
        if (nq .ne. -1) then
          G = (p/q0/c0)**(2.d0*dble(nq))*mdm*w**2/(2.d0*mu**dble(nq))*1./(1.+dble(nq)) &
          *((mu/muplus**2)**(dble(nq)+1.)-(u**2/w**2)**(dble(nq)+1.))
        else
          G = ((p)/q0/c0)**(2.d0*dble(nq))*mdm*w**2/(2.d0*mu**dble(nq))*log(mu/muplus**2*w**2/(u)**2)
        endif
        GFFI_H = G
      end function GFFI_H

      !generalized form factor: other elements
      function GFFI_A(w,vesc,A)
        double precision :: p, w,vesc,u,mN,A,Ei,B
        double precision :: dgamic,GFFI_A
        p = mdm*w
        u = sqrt(w**2-vesc**2)
        mN = A*mnuc
        Ei  = 5.8407d-2/(mN*(0.91*mN**(1./3.)+0.3)**2)
        B = .5*mdm*w**2/Ei/c0**2
        if (nq .eq. 0) then
          GFFI_A = Ei*c0**2*(exp(-mdm*u**2/2/Ei/c0**2)-exp(-B*mu/muplus**2))
        else
          GFFI_A = ((p)/q0/c0)**(2*dble(nq))*Ei*c0**2/(B*mu)**dble(nq)*(dgamic(1.+dble(nq),B*u**2/w**2) &
                  - dgamic(1.+dble(nq),B*mu/muplus**2))
        end if
      end function GFFI_A

      !SB: equation 2.19 0809.1871 as another way to get a const cross section capture rate
      function w_func_target(u, targetMass, vesc)
        double precision:: u, targetMass, vesc
        double precision:: w_func_target
        double precision:: B, Xi
        B = 3./2./u0**2.
        Xi = sqrt(pi)*(erf(sqrt(B)*(u-usun))-erf(sqrt(B)*(u+usun)))
        w_func_target = vesc**2*Xi/2/sqrt(B)&
                        +muminus**2/(4.*mu*B**(3./2.))*( 2*sqrt(B)*((usun-u)*exp(-B*(u+usun)**2.)&
                        +(usun+u)*exp(-B*(u+usun)**2+4.*B*usun*u)) -(1+2*B*usun**2)*Xi )
        ! print*, "From inside w_func_target"
        ! print*, muminus
      end function w_func_target

!----------------------------------------------------------------------------------------------------
!SB: Calculating the differential scattering rate using the formalism in appendix A of 1702.02768
    function diff_scattering_rate_cst(w, vesc) !Debugged for electrons, and hydrogen, and NREO :D
        implicit none
        double precision :: w, diff_scattering_rate_cst, vesc
        double precision :: wNat, vescNat, ueNat, v0Nat, x
        logical :: exist
        character(len=300) :: filename
        character ( len = 20 ) :: mdmString

        ! print*, "INSIDE: diff_scattering_rate_cst"

        ue = ue_at_r
        write(mdmString,  '(f10.5)') mdm

        ! w = 64484674.699495271d0
        ! vesc =  61773781.326245911d0
        ! ue = 674964.32622479054d0

        wNat = w/c0
        ueNat = ue/c0
        vescNat = vesc/c0
        v0Nat = v0/c0

        diff_scattering_rate_cst = 0d0

        diff_scattering_rate_cst = -(ueNat**2-2*muplus**2*vescNat**2 + 2*muminus**2*wNat**2)/muplus**2 &
                                    *(erf((muplus*vescNat-muminus*wNat)/ueNat) + erf((muplus*vescNat+muminus*wNat)/ueNat))

        diff_scattering_rate_cst = diff_scattering_rate_cst + 2*ueNat*vescNat/muplus/sqrt(pi) &
                                  *(exp(-(muplus*vescNat-muminus*wNat)**2/ueNat**2) + exp(-(muplus*vescNat+muminus*&
                                  wNat)**2/ueNat**2))

        diff_scattering_rate_cst = diff_scattering_rate_cst +    2*ueNat*muminus*wNat/muplus**2/sqrt(pi) &
                                    *(exp(-(muplus*vescNat-muminus*wNat)**2/ueNat**2) - exp(-(muplus*vescNat+muminus*&
                                    wNat)**2/ueNat**2))

        x =2*ueNat**2/mu*(erf((muminus*vescNat-muplus*wNat)/ueNat) + &
                                  erf((muminus*vescNat+muplus*wNat)/ueNat))&
                                  *exp(mu*(wNat**2-vescNat**2)/ueNat**2)

        if (isnan(x)) then
          diff_scattering_rate_cst = diff_scattering_rate_cst
        else
          diff_scattering_rate_cst = diff_scattering_rate_cst - x
        end if


        diff_scattering_rate_cst =   diff_scattering_rate_cst + 2*ueNat**2*muminus/mu/sqrt(mu+muminus**2)&
                                    *exp(mu*wNat**2*(mu+muminus**2-muplus**2)/(mu+muminus**2)/ueNat**2) &
                                    *(erf(((mu+muminus**2)*vescNat-muminus*muplus*wNat)/sqrt(mu+muminus**2)/ueNat) + &
                                    erf(((mu+muminus**2)*vescNat+muminus*muplus*wNat)/sqrt(mu+muminus**2)/ueNat) )

        diff_scattering_rate_cst = c0**2*(muplus)**2/mu*1./4.*diff_scattering_rate_cst

        if (diff_scattering_rate_cst.lt.0d0) then
          print*, "muplus: ", muplus
          print*, "muminus: ", muminus
          print*, "mu: ", mu
          print*, "w: ", w
          print*, "vesc: ", vesc
          print*, "ue: ", ue
          print*, "x: ", x
          print*, "Integral of Rminus cst: ", diff_scattering_rate_cst
        end if

      filename = "debugging_cap_rate/diff_scattering_rate_cst_mdm_"//trim(mdmString)//".dat"

        ! inquire(file=filename, exist=exist)
        ! if (exist) then
        !   open(12, file=filename, status="old", position="append", action="write")
        ! else
        !   open(12, file=filename, status="new", action="write")
        !   write(12, *) "DM mass [GeV] | ", "mu | ", "w [cm/s] | ", "vesc [cm/s] | ", "ue [cm/s] | ", "Rminus Int | "
        ! end if
        ! write(12, *) mdm, mu, w, vesc, ue, diff_scattering_rate_cst
        ! close(12)
    end function

    function diff_scattering_rate_v1(w, vesc) !Debugged for electrons, and hydrogen, and NREO :D
      implicit none
      double precision :: w, diff_scattering_rate_v1, vesc
      double precision :: x
      ue = ue_at_r
      !print*, "INSIDE: diff_scattering_rate_v1"
      diff_scattering_rate_v1 = 0d0
      diff_scattering_rate_v1 = 2./mu/muplus**3/sqrt(pi)/ue_at_r&
                            *(exp(-(muplus*vesc+muminus*w)**2/ue_at_r**2)&
                            *(2*muplus*ue_at_r**2*(muplus*vesc - muminus*w) + &
                            mu*(-2*muminus*ue_at_r**2*w + muplus**2*(2*vesc*w**2 + ue_at_r**2*(7*vesc + 4*w)) +&
                            muplus*(-2*muminus*w**3 + ue_at_r**2*(2*vesc + 2*w-7*muminus*w))))&
                            +exp(-(muplus*vesc-muminus*w)**2/ue_at_r**2)&
                            *(2*muplus*ue_at_r**2*(muplus*vesc + muminus*w) + &
                            mu*(2*muminus*ue_at_r**2*w + muplus**2*(2*vesc*w**2 + ue_at_r**2*(7*vesc - 4*w)) +&
                            muplus*(2*muminus*w**3 + ue_at_r**2*(2*vesc - 2*w+7*muminus*w)))))

      diff_scattering_rate_v1 = diff_scattering_rate_v1 + &
                            -1./mu/muplus**3/ue_at_r**2&
                            *(erf((muplus*vesc-muminus*w)/ue_at_r)+erf((muplus*vesc+muminus*w)/ue_at_r))&
                            *(2*muplus*ue_at_r**2*(ue_at_r**2 - 2*muplus**2*vesc**2 + 2*muminus**2*w**2) +&
                            mu*((2 + 7*muplus)*ue_at_r**4 - 4*muplus**3*vesc**2*w**2 +&
                            4*muminus**2*muplus*w**4 + 2*ue_at_r**2*(-3*muplus**3*vesc**2 +&
                            2*muminus**2*w**2 + (1 - 2*muminus + 7*muminus**2)*muplus*w**2 -&
                            4*muminus*muplus**2*w**2)))

      x = -2./mu**2*exp(mu*(w**2-vesc**2)/ue_at_r**2)*((4 + 3*mu)*ue_at_r**2 + 2*mu*vesc**2)&
          *(erf((muminus*vesc-muplus*w)/ue_at_r)+erf((muminus*vesc+muplus*w)/ue_at_r))

      if (isnan(x)) then
        diff_scattering_rate_v1 = diff_scattering_rate_v1
      else
        diff_scattering_rate_v1 = diff_scattering_rate_v1 +&
                                -2./mu**2*exp(mu*(w**2-vesc**2)/ue_at_r**2)*((4 + 3*mu)*ue_at_r**2 + 2*mu*vesc**2)&
                                *(erf((muminus*vesc-muplus*w)/ue_at_r)+erf((muminus*vesc+muplus*w)/ue_at_r))
      end if


      diff_scattering_rate_v1 = diff_scattering_rate_v1+&
                                2./mu**2/sqrt(mu+muminus**2)*muminus*ue_at_r**2*(4+3*mu)&
                                *exp(mu*w**2*(mu+muminus**2-muplus**2)/ue_at_r**2/(mu+muminus**2))&
                                *(erf((mu*vesc+muminus*(muminus*vesc-muplus*w))/ue_at_r/sqrt(mu+muminus**2))&
                                +erf((mu*vesc+muminus*(muminus*vesc+muplus*w))/ue_at_r/sqrt(mu+muminus**2)))

      diff_scattering_rate_v1 = diff_scattering_rate_v1+ &
                                4./sqrt(pi)/mu/(mu+muminus**2)**2*muminus*ue_at_r&
                                *(exp(-(mu*(vesc**2-w**2)+(muminus*vesc+muplus*w)**2)/ue_at_r**2)&
                                *(-mu*vesc - muminus**2*vesc + muminus*muplus*w)&
                                -exp(4*muminus*muplus*vesc*w/ue_at_r**2-(mu*(vesc**2-w**2)+(muminus*vesc+muplus*w)**2)/ue_at_r**2)&
                                *((mu + muminus**2)*vesc + muminus*muplus*w))

      diff_scattering_rate_v1 = diff_scattering_rate_v1+&
                                +2./mu/sqrt(mu+muminus**2)**5*muminus&
                                *exp(mu*w**2*(mu+muminus**2-muplus**2)/ue_at_r**2/(mu+muminus**2))&
                                *(erf((mu*vesc+muminus*(muminus*vesc-muplus*w))/ue_at_r/sqrt(mu+muminus**2))&
                                +erf((mu*vesc+muminus*(muminus*vesc+muplus*w))/ue_at_r/sqrt(mu+muminus**2)))&
                                *(mu*ue_at_r**2 + muminus**2*(ue_at_r**2 + 2*muplus**2*w**2))

      diff_scattering_rate_v1 =  (muplus)**2/mu*1./8.*ue_at_r**2/v0**2*diff_scattering_rate_v1
      ! if (diff_scattering_rate_v1.lt.0d0) then
      if (isnan(diff_scattering_rate_v1).or.(diff_scattering_rate_v1.lt.0d0)) then
        print*, "nq: ", nq
        print*, "nv: ", nv
        print*, "muplus: ", muplus
        print*, "muminus: ", muminus
        print*, "mu: ", mu
        print*, "w: ", w
        print*, "vesc: ", vesc
        print*, "ue: ", ue
        print*, "Integral of Rminus v1: ", diff_scattering_rate_v1
      end if
    end function

    function diff_scattering_rate_v2(w, vesc)
      implicit none
        double precision :: w, diff_scattering_rate_v2, vesc, x
        ue = ue_at_r
        ! print*, "INSIDE: diff_scattering_rate_v2"
        ! print*, "Needs Debugging!!!!!!"
        diff_scattering_rate_v2 = 1./(4.*mu*muplus**5)*exp(-(muplus*vesc-muminus*w)**2/ue**2)&
                              *(16*muplus**3*ue**5*(muplus*vesc + muminus*w) +&
                              8*mu*muplus**2*ue**3*(4*muminus*ue**2*w + 2*muminus*muplus*w**3 +&
                              muplus*ue**2*(4*vesc + (-4 + 5*muminus)*w) + muplus**2*vesc*(5*ue**2 + 2*w**2)) +&
                              4*mu**4*ue**3*(3*muplus*ue**2*(vesc - 2*w) + 2*muplus**3*(vesc - w)**3 +&
                              5*muminus*ue**2*w + 2*muminus**2*muplus*(vesc - 3*w)*w**2 +&
                              2*muminus**3*w**3 +&
                              2*muminus*muplus**2*w*(vesc**2 - 3*vesc*w + 3*w**2)) +&
                              mu**5*ue**3*(3*muplus*ue**2*(vesc - 2*w) + 2*muplus**3*(vesc - w)**3 +&
                              5*muminus*ue**2*w + 2*muminus**2*muplus*(vesc - 3*w)*w**2 +&
                              2*muminus**3*w**3 +&
                              2*muminus*muplus**2*w*(vesc**2 - 3*vesc*w + 3*w**2)) +&
                              mu**3*ue**3*(muplus*ue**2*(21*vesc - 34*w) +&
                              2*muminus**2*muplus*(7*vesc - 17*w)*w**2 +&
                              2*muplus**3*(vesc - w)*(14*ue**2 + 7*vesc**2 - 10*vesc*w + 7*w**2) +&
                              2*muminus*muplus**2*w*(14*ue**2 + 7*vesc**2 - 17*vesc*w + 17*w**2) +&
                              7*muminus*w*(5*ue**2 + 2*muminus**2*w**2)) +&
                              2*mu**2*ue*(20*muminus*ue**4*w + 8*muminus**3*ue**2*w**3 +&
                              8*muminus*muplus**2*ue**2*w*(5*ue**2 + vesc**2 - vesc*w + w**2) +&
                              4*muplus*ue**2*(ue**2*(3*vesc - 2*w) + 2*muminus**2*(vesc - w)*w**2) +&
                              muplus**4*vesc*(15*ue**4 + 20*ue**2*w**2 + 4*w**4) +&
                              muplus**3*(4*muminus*w**5 + 5*ue**4*(8*vesc + (-8 + 3*muminus)*w) +&
                              4*ue**2*(2*vesc**3 - 2*vesc**2*w +&
                              2*vesc*w**2 + (-2 + 5*muminus)*w**3))))
        diff_scattering_rate_v2 = diff_scattering_rate_v2 +&
                                1./(4.*mu*muplus**5)*exp(-(muplus*vesc+muminus*w)**2/ue**2)&
                                *(16*muplus**3*ue**5*(muplus*vesc - muminus*w) +&
                                8*mu*muplus**2*ue**3*(-4*muminus*ue**2*w - 2*muminus*muplus*w**3 +&
                                muplus*ue**2*(4*vesc + (4 - 5*muminus)*w) +&
                                muplus**2*vesc*(5*ue**2 + 2*w**2)) +&
                                4*mu**4*ue**3*(-5*muminus*ue**2*w - 2*muminus**3*w**3 +&
                                2*muplus**3*(vesc + w)**3 + 3*muplus*ue**2*(vesc + 2*w) +&
                                2*muminus**2*muplus*w**2*(vesc + 3*w) -&
                                2*muminus*muplus**2*w*(vesc**2 + 3*vesc*w + 3*w**2)) +&
                                mu**5*ue**3*(-5*muminus*ue**2*w - 2*muminus**3*w**3 +&
                                2*muplus**3*(vesc + w)**3 + 3*muplus*ue**2*(vesc + 2*w) +&
                                2*muminus**2*muplus*w**2*(vesc + 3*w) -&
                                2*muminus*muplus**2*w*(vesc**2 + 3*vesc*w + 3*w**2)) +&
                                mu**3*ue**3*(2*muminus**2*muplus*w**2*(7*vesc + 17*w) +&
                                muplus*ue**2*(21*vesc + 34*w) +&
                                2*muplus**3*(vesc + w)*(14*ue**2 + 7*vesc**2 + 10*vesc*w + 7*w**2) -&
                                2*muminus*muplus**2*w*(14*ue**2 + 7*vesc**2 + 17*vesc*w + 17*w**2) -&
                                7*(5*muminus*ue**2*w + 2*muminus**3*w**3)) +&
                                2*mu**2*ue*(-8*muminus*muplus**2*ue**2*w*(5*ue**2 + vesc**2 + vesc*w +&
                                w**2) - 4*muminus*ue**2*w*(5*ue**2 + 2*muminus**2*w**2) +&
                                muplus**4*vesc*(15*ue**4 + 20*ue**2*w**2 + 4*w**4) +&
                                4*muplus*ue**2*(2*muminus**2*w**2*(vesc + w) + ue**2*(3*vesc + 2*w)) +&
                                muplus**3*(-4*muminus*w**5 + 5*ue**4*(8*vesc + (8 - 3*muminus)*w) +&
                                4*ue**2*(2*vesc**3 + 2*vesc**2*w +&
                                2*vesc*w**2 + (2 - 5*muminus)*w**3))))

        diff_scattering_rate_v2  = diff_scattering_rate_v2+&
                              -1./(8.*mu*muplus**5)*sqrt(pi)*(erf((muplus*vesc-muminus*w)/ue)+erf((muplus*vesc+muminus*w)/ue))&
                              *(16*muplus**3*ue**4*(ue**2 - 2*muplus**2*vesc**2 + 2*muminus**2*w**2) + &
                              4*mu**4*(3*ue**6 + 6*(2*muminus**2 - 3*muminus*muplus + muplus**2)*ue**4*w**2 + &
                              4*muminus*(muminus - muplus)**3*ue**2*w**4) + &
                              mu**5*(3*ue**6 + &
                              6*(2*muminus**2 - 3*muminus*muplus + muplus**2)*ue**4*w**2 + &
                              4*muminus*(muminus - muplus)**3*ue**2*w**4) + &
                              mu**3*(7*(3 + 4*muplus**2)*ue**6 + &
                              2*(17*muplus**2 + 14*muminus**2*(3 + 2*muplus**2) - &
                              muminus*muplus*(51 + 28*muplus**2))*ue**4*w**2 + &
                              4*muminus*(7*muminus**3 - 17*muminus**2*muplus + &
                              17*muminus*muplus**2 - 7*muplus**3)*ue**2*w**4) - &
                              8*mu*muplus**2*ue**2*(-((4 + 5*muplus)*ue**4) + 4*muplus**3*vesc**2*w**2 - &
                              4*muminus**2*muplus*w**4 +&
                              2*ue**2*(5*muplus**3*vesc**2 - &
                              4*muminus**2*w**2 - (1 - 4*muminus + 5*muminus**2)*muplus*w**2)) + &
                              2*mu**2*((12 + 40*muplus**2 + 15*muplus**3)*ue**6 - &
                              8*muplus**5*vesc**2*w**4 + 8*muminus**2*muplus**3*w**6 + &
                              ue**4*(-30*muplus**5*vesc**2 + 48*muminus**2*w**2 - &
                              24*muminus*muplus*w**2 + 8*(1 + 10*muminus**2)*muplus**2*w**2 + &
                              10*(2 - 8*muminus + 3*muminus**2)*muplus**3*w**2) - &
                              4*ue**2*(10*muplus**5*vesc**2*w**2 - 4*muminus**4*w**4 + &
                              4*muminus**3*muplus*w**4 - &
                              4*muminus**2*muplus**2*w**4 - (1 - 4*muminus + &
                              10*muminus**2)*muplus**3*w**4)))
          x =  -1./2./mu**2*sqrt(pi)*ue**2*exp(mu*(w**2-vesc**2)/ue**2)&
                                  *(erf((muminus*vesc-muplus*w)/ue)+erf((muminus*vesc+muplus*w)/ue))&
                                  *((24 + 40*mu + 15*mu**2)*ue**4 + 4*mu*(4 + 5*mu)*ue**2*vesc**2 + 4*mu**2*vesc**4)
        if (isnan(x))  then
          x = 0d0
        end if
        diff_scattering_rate_v2 = diff_scattering_rate_v2 + x

        diff_scattering_rate_v2 = diff_scattering_rate_v2+&
                                exp(mu*(mu+muminus**2-muplus**2)*w**2/ue**2/(mu+muminus**2))&
                                *muminus*sqrt(pi)*ue**6/2./mu**2/sqrt(mu+muminus**2)&
                                *(24 + 40*mu + 15*mu**2)&
                                *(erf((mu*vesc+muminus*(muminus*vesc-muplus*w))/ue/sqrt(mu+muminus**2))&
                                +erf((mu*vesc+muminus*(muminus*vesc+muplus*w))/ue/sqrt(mu+muminus**2)))

        diff_scattering_rate_v2 = diff_scattering_rate_v2+&
                                  1./(2.*mu*(mu+muminus**2)**(9./2.))*muminus*ue**2&
                                  *(-2*(exp(-(mu*(vesc-w)*(vesc+w)+(muminus*vesc+muplus*w)**2)/ue**2)&
                                  +exp(-(mu*(vesc-w)*(vesc+w)+(muminus*vesc+muplus*w)**2)/ue**2+4*muplus*muminus*vesc*w/ue**2))&
                                  *((mu + muminus**2)**(3./2.)*ue*vesc*((mu + &
                                  muminus**2)*((mu*(11 + 10*mu) + 2*(4 + 5*mu)*muminus**2)*ue**2 +&
                                  2*mu*(mu + muminus**2)*vesc**2) + 2*mu*muminus**2*muplus**2*w**2)))

        diff_scattering_rate_v2 = diff_scattering_rate_v2+&
                                  1./(2.*mu*(mu+muminus**2)**(9./2.))*muminus*ue**2&
                                  *(-2*(-1*exp(-(mu*(vesc-w)*(vesc+w)+(muminus*vesc+muplus*w)**2)/ue**2)&
                                  +exp(-(mu*(vesc-w)*(vesc+w)+(muminus*vesc+muplus*w)**2)/ue**2+4*muplus*muminus*vesc*w/ue**2))&
                                  *muminus*sqrt(mu + muminus**2)*muplus*ue*w*((mu + &
                                  muminus**2)*((mu*(13 + 10*mu) + 2*(4 + 5*mu)*muminus**2)*ue**2 +&
                                  2*mu*(mu + muminus**2)*vesc**2) + 2*mu*muminus**2*muplus**2*w**2))

        diff_scattering_rate_v2 = diff_scattering_rate_v2+&
                                  1./(2.*mu*(mu+muminus**2)**(9./2.))*muminus*ue**2&
                                  *sqrt(pi)*exp(mu*w**2*(mu+muminus**2-muplus**2)/(mu+muminus**2)/ue**2)&
                                  *((mu + muminus**2)**2*(mu*(11 + 10*mu) + 2*(4 + 5*mu)*muminus**2)*ue**4 + &
                                  4*muminus**2*(mu + muminus**2)*(mu*(7 + 5*mu) + (4 + &
                                  5*mu)*muminus**2)*muplus**2*ue**2*w**2 + 4*mu*muminus**4*muplus**4*w**4) &
                                  *(erf(((mu+muminus**2)*vesc-muminus*muplus*w)/sqrt(mu+muminus**2)/ue)&
                                  +erf(((mu+muminus**2)*vesc+muminus*muplus*w)/sqrt(mu+muminus**2)/ue))


        diff_scattering_rate_v2 =  1./(4.*mu**2*sqrt(pi)*v0**4)*muplus**2*diff_scattering_rate_v2

        x =  -1./2./mu**2*sqrt(pi)*ue**2*exp(mu*(w**2-vesc**2)/ue**2)&
                                *(erf((muminus*vesc-muplus*w)/ue)+erf((muminus*vesc+muplus*w)/ue))&
                                *((24 + 40*mu + 15*mu**2)*ue**4 + 4*mu*(4 + 5*mu)*ue**2*vesc**2 + 4*mu**2*vesc**4)

        if (diff_scattering_rate_v2.lt.0d0) then
        ! if(isnan(x)) then
          print*, "nq: ", nq
          print*, "nv: ", nv
          print*, "muplus: ", muplus
          print*, "muminus: ", muminus
          print*, "mu: ", mu
          print*, "w: ", w
          print*, "vesc: ", vesc
          print*, "ue: ", ue
          print*, "x: ", x
          print*, "Integral of Rminus v2: ", diff_scattering_rate_v2
        end if
    end function

    function diff_scattering_rate_q1(w, vesc) !Debugged for electrons, and hydrogen, and NREO :D
      implicit none
      double precision :: w, diff_scattering_rate_q1, vesc
      double precision :: expr1, expr2, expr3, expr4, expr5, expr6, expr7, expr8
      logical :: exist
      character(len=300) :: filename
      character ( len = 20 ) :: mdmString
      ue = ue_at_r
      write(mdmString,  '(f10.5)') mdm
      diff_scattering_rate_q1 = 0d0

      expr1 = 2./mu/muplus**4/ue_at_r &
              *exp(-(muplus*vesc+muminus*w)**2/ue_at_r**2)&
              *(8.*muplus**2*ue_at_r**2*(muplus*vesc-muminus*w)&
              +mu*(5.*muminus*ue_at_r**2*w+2.*muminus**3*w**3 &
              -2.*muplus**3*(vesc**3-2.*vesc*w**2)&
              -muplus*(2.*muminus**2*vesc*w**2 &
              +ue_at_r**2*(3.*vesc+16.*muminus*w))&
              +2.*muplus**2*(8.*ue_at_r**2*(vesc+w)&
              +muminus*w*(vesc**2-2.*w**2))))

      expr2 = diff_scattering_rate_q1 &
              +1./mu/muplus**4/ue_at_r&
              *exp(-(muplus*vesc-muminus*w)**2/ue_at_r**2)&
              *(16.*muplus**2*ue_at_r**2*(muplus*vesc+muminus*w)&
              -2.*mu*(5.*muminus*ue_at_r**2*w+2.*muminus**3*w**3&
              +2.*muplus**3*(vesc**3-2.*vesc*w**2)&
              +muplus*(2.*muminus**2*vesc*w**2+ue_at_r**2*(3*vesc-16.*muminus*w))&
              -2.*muplus**2*(8.*ue_at_r**2*(vesc-w)+muminus*w*(-vesc**2+2.*w**2))))

      expr3 = diff_scattering_rate_q1 &
              +(erf((muplus*vesc-muminus*w)/ue_at_r)+erf((muplus*vesc+muminus*w)/ue_at_r))&
                 *sqrt(pi)*(3*ue_at_r**2/muplus**4-16.*ue_at_r**2/muplus**3&
                 -8.*ue_at_r**2/mu/muplus**2+12.*muminus**2*w**2/muplus**4&
                 -32.*muminus**2*w**2/muplus**3-4.*w**2/muplus**2&
                 +32.*muminus*w**2/muplus**2-16.*muminus**2*w**2/mu/muplus**2&
                 +4.*muminus**4*w**4/muplus**4/ue_at_r**2-8.*muminus**2*w**4/muplus**2/ue_at_r**2&
                 -4.*vesc**2/mu/ue_at_r**2*((-4.)*ue_at_r**2+mu*(vesc**2-2.*w**2)))

      expr4 = diff_scattering_rate_q1 + &
              -8.*sqrt(pi)/mu**2*exp(mu*(-vesc**2+w**2)/ue_at_r**2)&
              *(3*ue_at_r**2+mu*(vesc**2-w**2))&
              *(erf((muminus*vesc-muplus*w)/ue_at_r)+erf((muminus*vesc+muplus*w)/ue_at_r))

      if(isnan(expr4)) then
        expr4 = 0d0
      end if

      expr5 = diff_scattering_rate_q1 &
                                +exp(mu*w**2*(mu+muminus**2-muplus**2)/(mu+muminus**2)/ue_at_r**2)&
                                *(erf((mu*vesc+muminus*(muminus*vesc-muplus*w))/ue_at_r/sqrt(mu+muminus**2))&
                                +erf((mu*vesc+muminus*(muminus*vesc+muplus*w))/ue_at_r/sqrt(mu+muminus**2)))&
                                *(24.*muminus*sqrt(pi)*ue_at_r**2/mu**2/sqrt(mu+muminus**2)&
                                -8.*sqrt(pi)*muminus*w**2/mu/sqrt(mu+muminus**2))


      expr6 = diff_scattering_rate_q1 &
                          +4./mu/sqrt(mu+muminus**2)**(5.)*muminus&
                          *(-2.)*sqrt(mu+muminus**2)*ue_at_r&
                          *(mu*vesc+muminus**2*vesc)&
                          *(exp(-((muminus*vesc+muplus*w)**2+mu*(vesc**2-w**2))/ue_at_r**2)&
                  +exp(4.*muplus*muminus*vesc*w/ue_at_r**2-((muminus*vesc+muplus*w)**2+mu*(vesc**2-w**2))/ue_at_r**2))


      ! expr7 = diff_scattering_rate_q1 &
      !           +4./mu/sqrt(mu+muminus**2)**(5.)*muminus**2&
      !           *(-2.)*sqrt(mu+muminus**2)*ue_at_r*muplus*w&
      !           *(-exp(-((muminus*vesc+muplus*w)**2+mu*(vesc**2-w**2))/ue_at_r**2)&
      !   +exp(4*muplus*muminus*vesc*w/ue_at_r**2-((muminus*vesc+muplus*w)**2+mu*(vesc**2-w**2))/ue_at_r**2))

      expr7 = diff_scattering_rate_q1 &
                +4./mu/sqrt(mu+muminus**2)**(5.)*muminus**2&
                *(-2.)*sqrt(mu+muminus**2)*ue_at_r*muplus*w &
                *(exp(4*muplus*muminus*vesc*w/ue_at_r**2-((muminus*vesc+muplus*w)**2+mu*(vesc**2-w**2))/ue_at_r**2))&
                +4./mu/sqrt(mu+muminus**2)**(5.)*muminus**2&
                *(-2.)*sqrt(mu+muminus**2)*ue_at_r*muplus*w &
                *(-exp(-((1.+mu)*vesc+(mu-1.)*w)**2./4./ue_at_r**2))

      expr8 = diff_scattering_rate_q1 &
                                +4./mu/sqrt(mu+muminus**2)**5*muminus*sqrt(pi)&
                                *(mu*ue_at_r**2+muminus**2*(ue_at_r**2+2.*muplus**2*w**2))&
                                *(erf((mu*vesc+muminus*(muminus*vesc-muplus*w))/sqrt(mu+muminus**2)/ue_at_r)&
                                + erf((mu*vesc+muminus*(muminus*vesc+muplus*w))/sqrt(mu+muminus**2)/ue_at_r))&
                                *exp(-((muminus*vesc+muplus*w)**2+mu*(vesc**2-w**2))/ue_at_r**2&
                                +(mu*vesc+muminus*(muminus*vesc+muplus*w))**2/(mu+muminus**2)/ue_at_r**2)

      diff_scattering_rate_q1 = expr1 + expr2 + expr3 + expr4 +expr5 + expr6 + expr7 +expr8

      diff_scattering_rate_q1 =  (muplus)**2/mu*1./8./mu/sqrt(pi)/v0**2*muplus**2*ue_at_r**2*diff_scattering_rate_q1
      ! print*, "muplus: ", muplus
      ! print*, "muminus: ", muminus
      ! print*, "mu: ", mu
      ! print*, "w: ", w
      ! print*, "vesc: ", vesc
      ! print*, "ue: ", ue
      ! print*, "expr1: ", expr1
      ! print*, "expr2: ", expr2
      ! print*, "expr3: ", expr3
      ! print*, "expr4: ", expr4
      ! print*, "expr5: ", expr5
      ! print*, "expr6: ", expr6
      ! print*, "expr7: ", expr7
      ! print*, "expr8: ", expr8
      ! print*, "Integral of Rminus: ", diff_scattering_rate_q1

      ! if (diff_scattering_rate_q1.lt.0d0) then
      if(isnan(diff_scattering_rate_q1)) then
        print*, "nq: ", nq
        print*, "nv: ", nv
        print*, "muplus: ", muplus
        print*, "muminus: ", muminus
        print*, "mu: ", mu
        print*, "w: ", w
        print*, "vesc: ", vesc
        print*, "ue: ", ue
        print*, "expr1: ", expr1
        print*, "expr2: ", expr2
        print*, "expr3: ", expr3
        print*, "expr4: ", expr4
        print*, "expr5: ", expr5
        print*, "expr6: ", expr6
        print*, "expr7: ", expr7
        print*, "expr8: ", expr8
        print*, "Integral of Rminus q1: ", diff_scattering_rate_q1
      end if

      filename = "debugging_cap_rate/diff_scattering_rate_q1_mdm_"//trim(mdmString)//".dat"

        ! inquire(file=filename, exist=exist)
        ! if (exist) then
        !   open(12, file=filename, status="old", position="append", action="write")
        ! else
        !   open(12, file=filename, status="new", action="write")
        !   write(12, *) "DM mass [GeV] | ", "mu | ", "w [cm/s] | ", "vesc [cm/s] | ", "ue [cm/s] | ", "Rminus Int | " &
        !                 , "expr1 | ", "expr2 | ", "expr3 | ", "expr4 | ", "expr5 | "
        ! end if
        ! write(12, *) mdm, mu, w, vesc, ue, diff_scattering_rate_q1, expr1, expr2, expr3, expr4, expr5, &
        !               expr6, expr7, expr8
        ! close(12)
    end function

    function diff_scattering_rate_minus_q2(w, vesc) !Debugged for electrons, and hydrogen, and NREO :D
      implicit none
      double precision :: w, diff_scattering_rate_minus_q2, vesc, x, y, z
      logical :: exist
      character(len=300) :: filename
      double precision :: expr1, expr2, expr3, vdiff
      double precision :: v0Nat, vescNat, wNat, ueNat, vdiffNat
      character ( len = 20 ) :: mdmString
      integer :: pres
      double precision :: arrayMax(3)
      ue = ue_at_r
      write(mdmString,  '(f10.5)') mdm

       ! print*, "INSIDE: diff_scattering_rate_q2"
      ! w = 6.1820213723613910d7
      ! vesc = 6.1820210289703123d7
      ! ue = 5.4360700234251425d7

      pres = 3
      ! w = numSetPresicion(w, pres)
      ! vesc = numSetPresicion(vesc, pres)
      ! ue = numSetPresicion(ue, pres)

      wNat = w/(c0)
      vescNat = vesc/(c0)
      ueNat = ue/(c0)
      v0Nat = v0/(c0)
      ! print*, "wNat: ", wNat
      ! print*, "vescNat: ", vescNat
      ! print*, "ueNat: ", ueNat
      ! print*, "v0Nat: ", v0Nat

      ! print*, "vesc: ", vesc

      diff_scattering_rate_minus_q2 = 0d0

      expr1 = (ueNat*((12*ueNat**4*(5*(1 + mu)*(1 + &
              5*mu*(1 + mu)*(1 + mu + mu**2))*vescNat + &
               (5 + mu*(30 + mu*(75 + mu*(100 + 3*mu*(25 + mu)))))*wNat) + &
               2*mu*ueNat**2*((1 + mu)**3*(-3 + 2*mu*(-2 + 3*mu + 6*mu**2))*vescNat**3 + &
               3*(1 + mu)**2*(-1 + mu + 2*mu**2*(7 + mu*(13 +  &
               4*mu)))*vescNat**2*wNat + &
               3*(1 + mu)*(1 + mu*(10 + mu*(35 + 2*mu*(30 +  &
               mu*(11 + 2*mu)))))*vescNat*wNat**2 + &
               (3 + mu*(23 + mu*(75 + mu*(135 - 2*mu*(19 +  &
               3*mu)))))*wNat**3) + &
             mu**2*(vescNat + mu*vescNat + wNat - mu*wNat)*((1 + mu)**4*vescNat**4 -  &
             2*(1 + mu)**2*(1 + mu*(4 + mu))*vescNat**2* &
                wNat**2 + (1 + mu*(8 + mu*(30 + mu*(8 + mu))))*wNat**4))/ &
            exp(((1 + mu)*vescNat + (-1 + mu)*wNat)**2/(4*ueNat**2)) + &
           (12*ueNat**4*(5*(1 + mu)*(1 + 5*mu*(1 + mu)*(1 + mu +  &
           mu**2))*vescNat - &
               (5 + mu*(30 + mu*(75 + mu*(100 + 3*mu*(25 +  &
               mu)))))*wNat) + &
             2*mu*ueNat**2*((1 + mu)**3*(-3 + 2*mu*(-2 + 3*mu + 6*mu**2))*vescNat**3 - &
               3*(1 + mu)**2*(-1 + mu + 2*mu**2*(7 + mu*(13 +  &
               4*mu)))*vescNat**2*wNat + &
               3*(1 + mu)*(1 + mu*(10 + mu*(35 + 2*mu*(30 +  &
               mu*(11 + 2*mu)))))*vescNat*wNat**2 + &
               (-3 + mu*(-23 + mu*(-75 + mu*(-135 + 38*mu +  &
               6*mu**2))))*wNat**3) + &
             mu**2*(vescNat + mu*vescNat - wNat + mu*wNat)*((1 + mu)**4*vescNat**4 -  &
             2*(1 + mu)**2*(1 + mu*(4 + mu))*vescNat**2* &
                wNat**2 + (1 + mu*(8 + mu*(30 + mu*(8 + mu))))*wNat**4))/ &
            exp((vescNat + mu*vescNat + wNat -  &
            mu*wNat)**2/(4*ueNat**2))))/(64*mu**5*sqrt(pi)*v0Nat**4)

        expr2 = 0d0

        expr2 = mu**6*(60*ueNat**6 + 20*vescNat**6 - 60*vescNat**4*wNat**2 + 60*vescNat**2*wNat**4 + &
                44*wNat**6 + 72*ueNat**4*(3*vescNat**2 + 7*wNat**2) - &
                15*ueNat**2*(9*vescNat**4 - 18*vescNat**2*wNat**2 - 23*wNat**4)) !+sub10
        expr2 = expr2 + mu**9*(vescNat**2 - wNat**2)**3 !+sub8
        expr2 = expr2 + 3*mu**8*(-3*ueNat**2 + 2*vescNat**2 - &
                2*wNat**2)*(vescNat**2 - wNat**2)**2 !+sub7
        !if problems, rearrange these!!!!!!!!!!!!!
        expr2 = expr2 + 3*mu**7*(vescNat**2 - wNat**2)*(12*ueNat**4 - &
                18*ueNat**2*(vescNat**2 - wNat**2) + &
                  5*(vescNat**2 - wNat**2)**2) !+sub5
        expr2 = expr2 -60*ueNat**6 !+sub1
        expr2 = expr2 + mu*(-360*ueNat**6 + 36*ueNat**4*vescNat**2 - 36*ueNat**4*wNat**2) !+sub2
        expr2 = expr2 + mu**2*(-900*ueNat**6 + 216*ueNat**4*vescNat**2 - 9*ueNat**2*vescNat**4 - &
                216*ueNat**4*wNat**2 + 18*ueNat**2*vescNat**2*wNat**2 - 9*ueNat**2*wNat**4) !+sub3
        expr2 = expr2 -15*mu**5*(24*ueNat**6 - 36*ueNat**4*(vescNat**2 - wNat**2) + &
                12*ueNat**2*(vescNat**2 - wNat**2)**2 - (vescNat**2 - wNat**2)**3) !+sub4
        expr2 = expr2 -3*mu**4*(300*ueNat**6 - 240*ueNat**4*(vescNat**2 - wNat**2) + &
                45*ueNat**2*(vescNat**2 - wNat**2)**2 - 2*(vescNat**2 - wNat**2)**3)  !+sub6
        expr2 = expr2 -(mu**3*(1200*ueNat**6 - 540*ueNat**4*(vescNat**2 - wNat**2) + &
                54*ueNat**2*(vescNat**2 - wNat**2)**2 - (vescNat**2 - wNat**2)**3)) !+sub9
        expr2 =  (erf((vescNat + mu*vescNat + wNat - mu*wNat)/(2*ueNat)) + &
                  erf((vescNat + mu*vescNat - wNat + mu*wNat)/(2*ueNat)))/ &
                  (128*mu**6*v0Nat**4)*(expr2)
        ! print*, "********************"
        ! print*, "expr2: ", expr2
        ! print*, "********************"

        expr3 = (3*exp((mu*(-vescNat**2 + wNat**2))/ueNat**2)*(1 +  &
                mu)**6*ueNat**2*(20*ueNat**4 + 8*mu*ueNat**2*(vescNat**2 - wNat**2) + &
                mu**2*(vescNat**2 - wNat**2)**2)*(erf((vescNat -  &
                mu*vescNat + wNat + mu*wNat)/(2*ueNat)) - &
                erf(((-1 + mu)*vescNat + (1 +  &
                mu)*wNat)/(2*ueNat))))/(128*mu**6*v0Nat**4)

      if(isnan(expr3)) then
        expr3 = 0
      end if

      diff_scattering_rate_minus_q2 = c0**2*(expr1 + expr2 + expr3)

      vdiff = w - vesc
      !
      ! print*, "FROM INSIDE q2 RMINUS FUNCTION"
      ! print*, "mdm: ", mdm
      ! print*, "muplus: ", muplus
      ! print*, "muminus: ", muminus
      ! print*, "mu: ", mu
      ! print*, "w: ", w
      ! print*, "ue: ", ue
      ! print*, "vesc: ", vesc
      ! print*, "vdiff: ", vdiff
      ! print*, "expr1: ", c0**2*expr1
      ! print*, "expr2: ", c0**2*expr2
      ! print*, "expr3: ", c0**2*expr3
      ! print*, "Integral of Rminus: ", diff_scattering_rate_minus_q2

      if (diff_scattering_rate_minus_q2.lt.0d0) then
      ! if (isnan(diff_scattering_rate_minus_q2)) then
        print*, "nq: ", nq
        print*, "nv: ", nv
        print*, "mdm: ", mdm
        print*, "v0Nat: ", v0Nat
        print*, "muplus: ", muplus
        print*, "muminus: ", muminus
        print*, "mu: ", mu
        print*, "w: ", w
        print*, "vesc: ", vesc
        print*, "velDiff: ", w-vesc
        print*, "ue: ", ue
        print*, "expr1: ", c0**2*expr1
        print*, "expr2: ", c0**2*expr2
        print*, "expr3: ", c0**2*expr3
        print*, "Integral of Rminus q2: ", diff_scattering_rate_minus_q2
      end if

      filename = "debugging_cap_rate/diff_scattering_rate_minus_q2_mdm_"//trim(mdmString)//".dat"
      ! ! if (diff_scattering_rate_minus_q2.lt.0d0) then
        ! inquire(file=filename, exist=exist)
        ! if (exist) then
        !   open(12, file=filename, status="old", position="append", action="write")
        ! else
        !   open(12, file=filename, status="new", action="write")
        !   write(12, *) "DMmass [GeV] | ", "mu | ", "w [cm/s] | ", "vesc [cm/s] | ", "ue [cm/s] | ", &
        !               "expr1 | ", "expr2 | ", "expr3 | ", "Rminus Int | "
        ! end if
        ! write(12, *) mdm, mu, w, vesc, ue, c0**2*expr1, c0**2*expr2, c0**2*expr3, diff_scattering_rate_minus_q2
        ! close(12)
      ! ! end if
    end function

    function diff_scattering_rate_q3(w, vesc)
      implicit none
        double precision :: w, diff_scattering_rate_q3, vesc, x1, x2
        double precision:: expr1, expr2
        double precision wNat, vescNat, ueNat, v0Nat, x, y, z
        double precision :: sqrtPi
        ue = ue_at_r

        sqrtPi = 1.772453850905516d0

        ! w = 137920111.81111661d0
        ! vesc = 137920075.03898156d0
        ! ue = 50304493.749889053d0

        wNat = w/c0
        ueNat = ue/c0
        vescNat = vesc/c0
        v0Nat = v0/c0

        ! print*, "INSIDE: diff_scattering_rate_q3"

        x = exp((mu*(-vescNat**2 + wNat**2))/ueNat**2)*(erf((vescNat - &
            mu*vescNat + wNat + mu*wNat)/(2*ueNat)) -  &
            erf(((-1 + mu)*vescNat + (1 + mu)*wNat)/(2*ueNat)))
        if (isnan(x)) then
          x = 0d0
        end if

        y = erf((vescNat + mu*vescNat + wNat - mu*wNat)/(2*ueNat)) &
            + erf((vescNat + mu*vescNat - wNat + mu*wNat)/(2*ueNat))

        expr1 = 0d0

        expr1 = (840*ueNat**6*(vescNat - wNat))/exp((vescNat + &
                mu*vescNat + wNat - mu*wNat)**2/(4*ueNat**2)) + &
                (840*ueNat**6*(vescNat + wNat))/exp(((1 + mu)*vescNat + &
                (-1 + mu)*wNat)**2/(4*ueNat**2))
        expr1 = expr1 + mu*((60*ueNat**4*(vescNat - wNat)*(112*ueNat**2 - &
                vescNat**2 + wNat**2))/exp((vescNat + mu*vescNat + wNat - &
                mu*wNat)**2/(4*ueNat**2)) + &
                (60*ueNat**4*(vescNat + wNat)*(112*ueNat**2 - vescNat**2 + &
                wNat**2))/exp(((1 + mu)*vescNat + (-1 + mu)*wNat)**2/(4*ueNat**2)))

        expr1 = expr1 + mu**11*((32*ueNat**2*vescNat*(vescNat - &
                wNat)**4)/exp((vescNat + mu*vescNat + &
                wNat - mu*wNat)**2/(4*ueNat**2)) + &
                (32*ueNat**2*vescNat*(vescNat + wNat)**4)/exp(((1 + &
                mu)*vescNat + (-1 + mu)*wNat)**2/(4*ueNat**2)))

        expr1 = expr1 + mu**2*((20*ueNat**2*(vescNat - wNat)*(1176*ueNat**4 + &
                (vescNat**2 - wNat**2)**2 + ueNat**2*(-17*vescNat**2 - &
                14*vescNat*wNat + 31*wNat**2)))/ &
                exp((vescNat + mu*vescNat + wNat - mu*wNat)**2/(4*ueNat**2)) + &
              (20*ueNat**2*(vescNat + wNat)*(1176*ueNat**4 + (vescNat**2 - &
              wNat**2)**2 + ueNat**2*(-17*vescNat**2 + 14*vescNat*wNat + 31*wNat**2)))/ &
               exp(((1 + mu)*vescNat + (-1 + mu)*wNat)**2/(4*ueNat**2)))
        expr1 = expr1 + mu**6*(-(((vescNat + wNat)*(-23520*ueNat**6 - 112*ueNat**2*(vescNat &
                + wNat)**2*(11*vescNat**2 - 18*vescNat*wNat + &
                21*wNat**2) - 280*ueNat**4*(23*vescNat**2 + &
                70*vescNat*wNat + 47*wNat**2) + &
              (vescNat + wNat)**3*(35*vescNat**3 - 135*vescNat**2*wNat +  &
              185*vescNat*wNat**2 - 93*wNat**3)))/ &
              exp(((1 + mu)*vescNat + (-1 + mu)*wNat)**2/(4*ueNat**2))) - &
              ((vescNat - wNat)*(-23520*ueNat**6 - 112*ueNat**2*(vescNat - &
              wNat)**2*(11*vescNat**2 + 18*vescNat*wNat + 21*wNat**2) - &
             280*ueNat**4*(23*vescNat**2 - 70*vescNat*wNat +  &
             47*wNat**2) + (vescNat - wNat)**3* &
            (35*vescNat**3 + 135*vescNat**2*wNat &
             + 185*vescNat*wNat**2 + 93*wNat**3)))/ &
             exp((vescNat + mu*vescNat + wNat -  &
             mu*wNat)**2/(4*ueNat**2)))
        expr1 = expr1 + mu**5*(-(((vescNat + wNat)*(-47040*ueNat**6 - &
                56*ueNat**2*(vescNat + wNat)**2*(17*vescNat**2 - 36*vescNat*wNat +&
                  27*wNat**2) - 280*ueNat**4*(13*vescNat**2 + &
                  56*vescNat*wNat + 43*wNat**2) + &
                (vescNat + wNat)**3*(21*vescNat**3 - &
                75*vescNat**2*wNat + 91*vescNat*wNat**2 - 37*wNat**3)))/ &
              exp(((1 + mu)*vescNat + (-1 + mu)*wNat)**2/(4*ueNat**2))) - &
            ((vescNat - wNat)*(-47040*ueNat**6 - 56*ueNat**2*(vescNat -  &
            wNat)**2*(17*vescNat**2 + 36*vescNat*wNat + 27*wNat**2) - &
               280*ueNat**4*(13*vescNat**2 - 56*vescNat*wNat +  &
               43*wNat**2) + (vescNat - wNat)**3* &
                (21*vescNat**3 + 75*vescNat**2*wNat  + &
                91*vescNat*wNat**2 + 37*wNat**3)))/ &
             exp((vescNat + mu*vescNat + wNat - mu*wNat)**2/(4*ueNat**2)))

      expr1 =  expr1 +  mu**4*(-(((vescNat + wNat)*(-58800*ueNat**6 + &
                (7*vescNat - 9*wNat)*(vescNat - wNat)**2*(vescNat + wNat)**3 - &
                560*ueNat**4*(vescNat**2 + 14*vescNat*wNat +  &
                13*wNat**2) - 2*ueNat**2*(vescNat + wNat)**2* &
                 (247*vescNat**2 - 546*vescNat*wNat + 327*wNat**2)))/exp(((1 +  &
                 mu)*vescNat + (-1 + mu)*wNat)**2/(4*ueNat**2))) - &
                 ((vescNat - wNat)*(-58800*ueNat**6 + (vescNat -  &
                wNat)**3*(vescNat + wNat)**2*(7*vescNat + 9*wNat) - &
               560*ueNat**4*(vescNat**2 - 14*vescNat*wNat + 13*wNat**2) -  &
               2*ueNat**2*(vescNat - wNat)**2* &
                (247*vescNat**2 + 546*vescNat*wNat +  &
                327*wNat**2)))/exp((vescNat + mu*vescNat + wNat - mu*wNat)**2/(4*ueNat**2)))

        expr1 = expr1 +  mu**3*(-(((vescNat + wNat)*(-47040*ueNat**6 + &
                560*ueNat**4*(vescNat**2 - 4*vescNat*wNat - 5*wNat**2) + &
                  (vescNat**2 - wNat**2)**3 - 10*ueNat**2*(vescNat + wNat)**2*(15*vescNat**2 -  &
                  32*vescNat*wNat + 17*wNat**2)))/ &
                exp(((1 + mu)*vescNat + (-1 + mu)*wNat)**2/(4*ueNat**2))) - &
              ((vescNat - wNat)*(-47040*ueNat**6 + 560*ueNat**4*(vescNat**2 +  &
              4*vescNat*wNat - 5*wNat**2) + (vescNat**2 - wNat**2)**3 - &
                 10*ueNat**2*(vescNat - wNat)**2*(15*vescNat**2 + 32*vescNat*wNat +  &
                 17*wNat**2)))/ &
               exp((vescNat + mu*vescNat + wNat - mu*wNat)**2/(4*ueNat**2)))

        expr1 = expr1 + mu**10*(-(((vescNat - wNat)**3*(-64*ueNat**2*vescNat*(3*vescNat - &
                4*wNat) + (vescNat + wNat)**4))/ &
                exp((vescNat + mu*vescNat + wNat -  &
                mu*wNat)**2/(4*ueNat**2))) - &
              ((vescNat + wNat)**3*((vescNat -  &
              wNat)**4 - 64*ueNat**2*vescNat*(3*vescNat + 4*wNat)))/ &
               exp(((1 + mu)*vescNat + (-1 + mu)*wNat)**2/(4*ueNat**2)))

        expr1 = expr1 + mu**9*(-(((vescNat + wNat)**2*(-64*ueNat**4*(10*vescNat + wNat) + &
                & (vescNat - wNat)**4*(7*vescNat + 9*wNat) - &
                4*ueNat**2*(131*vescNat**3 + 295*vescNat**2*wNat +  &
                249*vescNat*wNat**2 - 3*wNat**3)))/ &
              exp(((1 + mu)*vescNat + (-1 + &
               mu)*wNat)**2/(4*ueNat**2))) - &
            ((vescNat - wNat)**2*(64*ueNat**4*(-10*vescNat +  &
            wNat) + (7*vescNat - 9*wNat)*(vescNat + wNat)**4 - &
               4*ueNat**2*(131*vescNat**3 - 295*vescNat**2*wNat +  &
               249*vescNat*wNat**2 + 3*wNat**3)))/ &
             exp((vescNat + mu*vescNat + wNat -  &
             mu*wNat)**2/(4*ueNat**2)))

      expr1 = expr1 + mu**8*(-(((vescNat + wNat)*((vescNat - &
              wNat)**4*(21*vescNat**2 + 54*vescNat*wNat + 37*wNat**2) - &
              4*ueNat**4*(755*vescNat**2 + 1402*vescNat*wNat + 83*wNat**2) - &
              2*ueNat**2*(455*vescNat**4 + 1156*vescNat**3*wNat +  &
              1642*vescNat**2*wNat**2 + 1284*vescNat*wNat**3 - 57*wNat**4)))/ &
              exp(((1 + mu)*vescNat + (-1 + mu)*wNat)**2/(4*ueNat**2))) - &
              ((vescNat - wNat)*((vescNat + wNat)**4*(21*vescNat**2 -  &
              54*vescNat*wNat + 37*wNat**2) - &
             4*ueNat**4*(755*vescNat**2 - 1402*vescNat*wNat +  &
             83*wNat**2) + &
             ueNat**2*(-910*vescNat**4 + 2312*vescNat**3*wNat -  &
             3284*vescNat**2*wNat**2 + 2568*vescNat*wNat**3 + 114*wNat**4)))/ &
             exp((vescNat + mu*vescNat + wNat - mu*wNat)**2/(4*ueNat**2)))

        expr1 = expr1 + mu**7*((24*ueNat**6*(245*vescNat + 59*wNat) - &
                (vescNat - wNat)**4*(35*vescNat**3 + 135*vescNat**2*wNat +  &
                185*vescNat*wNat**2 + 93*wNat**3) +  &
                4*ueNat**4*(1505*vescNat**3 + 5181*vescNat**2*wNat + &
                 5283*vescNat*wNat**2 +  &
                127*wNat**3) + 2*ueNat**2*(595*vescNat**5 +  &
                1409*vescNat**4*wNat + 1854*vescNat**3*wNat**2 + &
                 2882*vescNat**2*wNat**3 +  &
                2463*vescNat*wNat**4 - 243*wNat**5))/exp(((1 + mu)*vescNat + &
                 (-1 + mu)*wNat)**2/(4*ueNat**2)) +  &
                (24*ueNat**6*(245*vescNat - 59*wNat) + 4*ueNat**4*(1505*vescNat**3 -  &
                5181*vescNat**2*wNat + 5283*vescNat*wNat**2 -  &
                127*wNat**3) - (vescNat + wNat)**4*(35*vescNat**3 -  &
                135*vescNat**2*wNat + 185*vescNat*wNat**2 - 93*wNat**3) +  &
                2*ueNat**2*(595*vescNat**5 - 1409*vescNat**4*wNat +  &
                1854*vescNat**3*wNat**2 - 2882*vescNat**2*wNat**3 +  &
                2463*vescNat*wNat**4 + 243*wNat**5))/exp((vescNat + mu*vescNat + wNat - mu*wNat)**2/(4*ueNat**2)))

        expr1 = (ueNat/(256.*Sqrt(pi)))*expr1

        expr2 = ((-840*ueNat**8 - mu**12*(vescNat**2 - wNat**2)**4 + 8*mu**11*(vescNat**2 -  &
              wNat**2)**3*(2*ueNat**2 - vescNat**2 + wNat**2) - 480*mu*ueNat**6*(14*ueNat**2 - vescNat**2 + wNat**2) +  &
              mu**4*(-58800*ueNat**8 + 26880*ueNat**6*(vescNat - wNat)*(vescNat + wNat) -  &
              3360*ueNat**4*(vescNat**2 - wNat**2)**2 + 128*ueNat**2*(vescNat**2 - wNat**2)**3 - (vescNat**2 - wNat**2)**4) -  &
              8*mu**5*(5880*ueNat**8 - 4200*ueNat**6*(vescNat - wNat)*(vescNat + wNat) +  &
              840*ueNat**4*(vescNat**2 - wNat**2)**2 - 56*ueNat**2*(vescNat**2 - wNat**2)**3 +  &
              (vescNat**2 - wNat**2)**4) -  &
              120*mu**2*ueNat**4*(196*ueNat**4 + (vescNat**2 - wNat**2)**2 + 32*ueNat**2*(-vescNat**2 +  &
              wNat**2)) - 4*mu**10*(vescNat - wNat)**2*(vescNat + wNat)**2*(30*ueNat**4 + 7*(vescNat**2 -  &
              wNat**2)**2 + 32*ueNat**2*(-vescNat**2 + wNat**2)) +  &
              8*mu**9*(vescNat - wNat)*(vescNat + wNat)*(60*ueNat**6 + 56*ueNat**2*(vescNat**2 -  &
              wNat**2)**2 - 7*(vescNat**2 - wNat**2)**3 + 120*ueNat**4*(-vescNat**2 + wNat**2)) -  &
              56*mu**7*(120*ueNat**8 + 120*ueNat**4*(vescNat**2 - wNat**2)**2 - 20*ueNat**2*(vescNat**2 -  &
              wNat**2)**3 + (vescNat**2 - wNat**2)**4 + 240*ueNat**6*(-vescNat**2 + wNat**2)) -  &
              16*mu**3*(2940*ueNat**8 + 60*ueNat**4*(vescNat**2 - wNat**2)**2 - ueNat**2*(vescNat**2 -  &
              wNat**2)**3 + 840*ueNat**6*(-vescNat**2 + wNat**2)) -  &
              28*mu**6*(840*ueNat**8 + 300*ueNat**4*(vescNat**2 - wNat**2)**2 - 32*ueNat**2*(vescNat**2 - &
               wNat**2)**3 + (vescNat**2 - wNat**2)**4 + 960*ueNat**6*(-vescNat**2 + wNat**2)) +  &
              2*mu**8*(420*ueNat**8 - 35*vescNat**8 + 140*vescNat**6*wNat**2 - 210*vescNat**4*wNat**4 +  &
              140*vescNat**2*wNat**6 + 93*wNat**8 - 1680*ueNat**4*(vescNat**2 - 3*wNat**2)*(vescNat**2 +  &
              wNat**2) + 960*ueNat**6*(2*vescNat**2 + 5*wNat**2) +  &
              448*ueNat**2*(vescNat**6 - 3*vescNat**4*wNat**2 + 3*vescNat**2*wNat**4 + &
              3*wNat**6)))*(erf((vescNat + mu*vescNat + wNat - mu*wNat)/(2*ueNat)) +  &
              erf((vescNat + mu*vescNat - wNat + mu*wNat)/(2*ueNat))) +  &
              4*x*(1 + mu)**8*ueNat**2*(210*ueNat**6 +  &
              90*mu*ueNat**4*(vescNat - wNat)*(vescNat + wNat) + 15*mu**2*ueNat**2*(vescNat**2 -  &
              wNat**2)**2 + mu**3*(vescNat**2 - wNat**2)**3))/(512.*mu)

      diff_scattering_rate_q3 = c0**2*(expr1+expr2)/(mu**7*v0Nat**6)

        ! print*, "mu: ", mu
        ! print*, "w: ", w
        ! print*, "vesc: ", vesc
        ! print*, "ue: ", ue
        ! print*, "expr1: ", expr1
        ! print*, "expr2: ", expr2
        ! print*, "Integral of Rminus q3: ", diff_scattering_rate_q3
        if (diff_scattering_rate_q3.lt.0d0) then
        ! if (isnan(diff_scattering_rate_q3)) then
          print*, "muplus: ", muplus
          print*, "muminus: ", muminus
          print*, "nq: ", nq
          print*, "nv: ", nv
          print*, "mdm: ", mdm
          print*, "mu: ", mu
          print*, "w: ", w
          print*, "vesc: ", vesc
          print*, "ue: ", ue
          print*, "expr1: ", expr1
          print*, "expr2: ", expr2
          print*, "Integral of Rminus q3: ", diff_scattering_rate_q3
        end if
    end function

    function diff_scattering_rate_q1v1(w, vesc) !Debugged for electrons, and hydrogen, and NREO :D
      implicit none
      double precision :: w, diff_scattering_rate_q1v1, vesc
      logical :: exist
      character(len=300) :: filename
      double precision :: expr1, expr2, expr3, expr4, expr5, x
      double precision :: wNat, vescNat, ueNat, v0Nat
      character ( len = 20 ) :: mdmString
      ue = ue_at_r
      write(mdmString,  '(f10.5)') mdm
      ! print*, "INSIDE: diff_scattering_rate_q1v1"
      diff_scattering_rate_q1v1 = 0d0

      ! w = 136010275.103128520d0
      ! vesc = 135880863.89456472d0
      ! ue = 2094478787.8112397d0

      wNat = w/c0
      vescNat = vesc/c0
      ueNat = ue/c0
      v0Nat = v0/c0

      ! expr1 = sqrt(pi)*(6*(-8 + mu*(-37 + mu*(-68 + mu*(-62 + &
      !                             mu*(-28 + 5*mu)))))*ueNat**6 + &
      !                             4*mu*ueNat**4*((1 + mu)**4*(6 + 5*mu)*vescNat**2 - (9 + &
      !                             mu*(41 + mu*(74 + (-11 + mu)*mu*(-6 + 5*mu))))*wNat**2) - &
      !                             2*mu**3*wNat**2*((1 + mu)**4*vescNat**4 - &
      !                             2*(1 + mu)**4*vescNat**2*wNat**2 + (-1 + mu)**2*(1 + mu*(6 + mu))*wNat**4) - &
      !                             mu**2*ueNat**2*((1 + mu)**4*(4 + 5*mu)*vescNat**4 - &
      !                             2*(1 + mu)**4*(8 + 5*mu)*vescNat**2*wNat**2 + (12 + &
      !                             mu*(53 + mu*(92 + mu*(-162 + mu*(32 + 5*mu)))))*wNat**4)) &
      !                             *(erf((muplus*vescNat-muminus*wNat)/ueNat) + erf((muplus*vescNat+muminus*wNat)/ueNat))
      !
      ! expr2 = 2*mu*ueNat * &
      !         (exp(-(muplus*vescNat+muminus*wNat)**2/ueNat**2)* &
      !         (6*ueNat**4*((1 + mu)*(8 + mu*(29 + mu*(39 + 23*mu)))*vescNat + (8 + &
      !         mu*(37 + mu*(68 + mu*(62 + mu))))*wNat) - &
      !         2*mu**2*wNat**2*(vescNat + mu*vescNat + wNat - &
      !         mu*wNat)*((1 + mu)**2*vescNat**2 - (1 + mu*(6 + mu))*wNat**2) + &
      !         mu*ueNat**2*(mu*(1 + mu)**3*(3 + 8*mu)*vescNat**3 + &
      !         mu*(1 + mu)**2*(19 + mu*(53 + 16*mu))*vescNat**2*wNat + (1 + mu)*(12 + &
      !         mu*(65 + mu*(138 + mu*(45 + 8*mu))))*vescNat*wNat**2 + (12 + &
      !         mu*(61 + mu*(129 - mu*(37 + 5*mu))))*wNat**3)))
      !
      ! expr2 = expr2 + 2*mu*ueNat *(exp(-(muplus*vescNat-muminus*wNat)**2/ueNat**2)&
      !         *(6*ueNat**4*((1 + mu)*(8 + mu*(29 + mu*(39 + 23*mu)))*vescNat - (8 + &
      !         mu*(37 + mu*(68 + mu*(62 + mu))))*wNat) - &
      !         2*mu**2*wNat**2*(vescNat + mu*vescNat - wNat + &
      !         mu*wNat)*((1 + mu)**2*vescNat**2 - (1 + mu*(6 + mu))*wNat**2) + &
      !         mu*ueNat**2*(mu*(1 + mu)**3*(3 + 8*mu)*vescNat**3 - &
      !         mu*(1 + mu)**2*(19 + mu*(53 + 16*mu))*vescNat**2*wNat + (1 + mu)*(12 + &
      !         mu*(65 + mu*(138 + mu*(45 + 8*mu))))*vescNat*wNat**2 + (-12 + &
      !         mu*(-61 + mu*(-129 + mu*(37 + 5*mu))))*wNat**3)))
      !
      ! expr3 = 2**5*exp(mu*(wNat**2-vescNat**2)/ueNat**2) &
      !       *muplus**4*sqrt(pi)*ueNat**2* (3*(8 + 5*mu)*ueNat**4 + &
      !       2*mu**2*vescNat**2*(vescNat - wNat)*(vescNat + wNat) + &
      !       mu*ueNat**2*((12 + 5*mu)*vescNat**2 - (6 + 5*mu)*wNat**2)) &
      !       *(erf((muplus*wNat-muminus*vescNat)/ueNat) - erf((muplus*wNat+muminus*vescNat)/ueNat))


      expr1 = 2*mu*ueNat*((2*ueNat**4*((1 + mu)*(12 + mu*(82 + &
              mu*(223 + mu*(203 + 74*mu))))*vescNat + (12 + mu*(110 + mu*(379 &
              + mu*(380 + (177 - 2*mu)*mu))))*wNat) -  &
              2*mu**2*(1 + mu)*wNat**2*(vescNat + mu*vescNat + wNat - mu*wNat)*((1 +  &
              mu)**2*vescNat**2 - (1 + mu*(6 + mu))*wNat**2) +  &
              mu*ueNat**2*((1 + mu)**3*(-4 + mu*(7 + mu*(11 + 8*mu)))*vescNat**3 + (1 +  &
              mu)**2*(-4 + mu*(27 + mu*(68 + mu*(69 + 16*mu))))*vescNat**2*wNat + (1 +  &
              mu)*(12 + mu*(93 + mu*(187 + mu*(183 + mu*(53 + 8*mu)))))*vescNat*wNat**2 -  &
              (-3 + mu)*(4 + mu*(31 + mu*(63 + mu*(57 + 5*mu))))*wNat**3))/exp(((1 +  &
              mu)*vescNat + (-1 + mu)*wNat)**2/(4*ueNat**2)) +  &
              (2*ueNat**4*((1 + mu)*(12 + mu*(82 + mu*(223 + mu*(203 + 74*mu))))*vescNat  &
              + (-12 + mu*(-110 + mu*(-379 + mu*(-380 + mu*(-177 + 2*mu)))))*wNat) -  &
              2*mu**2*(1 + mu)*wNat**2*(vescNat + mu*vescNat - wNat + mu*wNat)*((1 +  &
              mu)**2*vescNat**2 - (1 + mu*(6 + mu))*wNat**2) +  &
              mu*ueNat**2*((1 + mu)**3*(-4 + mu*(7 + mu*(11 + 8*mu)))*vescNat**3 - (1 +  &
              mu)**2*(-4 + mu*(27 + mu*(68 + mu*(69 + 16*mu))))*vescNat**2*wNat + (1 + &
               mu)*(12 + mu*(93 + mu*(187 + mu*(183 + mu*(53 + 8*mu)))))*vescNat*wNat**2 +  &
              (-3 + mu)*(4 + mu*(31 + mu*(63 + mu*(57 + 5*mu))))*wNat**3))/exp((vescNat  &
              + mu*vescNat + wNat - mu*wNat)**2/(4.*ueNat**2)))

      expr2 = -(Sqrt(pi)*((48 + 2*mu*(111 + mu*(257 + mu*(428 + mu*(304 &
                + (79 - 15*mu)*mu)))))*ueNat**6 + 2*mu*ueNat**4*(-2*(1 + mu)**5*(6 + &
                5*mu)*vescNat**2 + (6 + mu*(87 + mu*(335 + mu*(190 + mu*(8 + 15*(-7 + &
                mu)*mu)))))*wNat**2) + &
                2*mu**3*wNat**2*((1 + mu)**5*vescNat**4 - 2*(1 + mu)**5*vescNat**2*wNat**2 + (-1 + &
                mu)**2*(1 + mu)*(1 + mu*(6 + mu))*wNat**4) + &
                mu**2*ueNat**2*((1 + mu)**5*(4 + 5*mu)*vescNat**4 - 2*(1 + mu)**5*(8 + 5*mu)* &
                vescNat**2*wNat**2 + (12 + mu*(81 + mu*(97 + mu*(-22 + mu*(-146 +  &
                mu*(37 + 5*mu))))))*wNat**4))* &
                (erf((vescNat + mu*vescNat + wNat - mu*wNat)/(2*ueNat)) + erf((vescNat + &
                mu*vescNat - wNat + mu*wNat)/(2.*ueNat))))

      expr3 = 4*(1 - mu)*mu*ueNat**3*((ueNat**2*((12 + 53*mu + 63*mu**2 + 27*mu**3 + &
              5*mu**4)*vescNat + (-12 - 37*mu + 27*mu**2 + 17*mu**3 + 5*mu**4)*wNat) + &
              2*mu*((1 + mu)**3*vescNat**3 + (-1 + mu)*(1 + mu)**2*vescNat**2*wNat - &
              4*mu*(1 + mu)*vescNat*wNat**2 - 4*(-1 + mu)*mu*wNat**3))/ &
              exp((vescNat + mu*vescNat + wNat - mu*wNat)**2/(4*ueNat**2)) + (ueNat**2*((12 + &
              53*mu + 63*mu**2 + 27*mu**3 + 5*mu**4)*vescNat + (12 + 37*mu - 27*mu**2 - 17*mu**3 - &
              5*mu**4)*wNat) + 2*mu*((1 + mu)**3*vescNat**3 - (-1 + mu)*(1 + mu)**2*vescNat**2*wNat - &
              4*mu*(1 + mu)*vescNat*wNat**2 + 4*(-1 + mu)*mu*wNat**3))/ &
              exp(((1 + mu)*vescNat + (-1 + mu)*wNat)**2/(4*ueNat**2)))

      expr4 = 2*(1 - mu)*mu*Sqrt(pi)*ueNat**2*(2*(12 + 41*mu + 22*mu**2 + 5*mu**3)*ueNat**4 + &
              (12 + 25*mu - 80*mu**2 + 10*mu**3 + 12*mu**4 + 5*mu**5)*ueNat**2*wNat**2 - 8*(-1 + &
              mu)**2*mu**2*wNat**4)*(-erf((vescNat + mu*vescNat + wNat - mu*wNat)/(2*ueNat)) - &
              erf((vescNat + mu*vescNat - wNat + mu*wNat)/(2*ueNat)))

      expr5 = 2*exp((mu*(-vescNat**2 + wNat**2))/ueNat**2)*(1 + mu)**5*Sqrt(pi) &
              *ueNat**2*(3*(8 + 5*mu)*ueNat**4 + 2*mu**2*vescNat**2*(vescNat - &
              wNat)*(vescNat + wNat) + mu*ueNat**2*((12 + 5*mu)*vescNat**2 - &
               (6 + 5*mu)*wNat**2))* &
              (-erf(((-1 + mu)*vescNat + wNat + mu*wNat)/(2*ueNat)) + &
              erf((vescNat - mu*vescNat + wNat + mu*wNat)/(2*ueNat)))

      if (isnan(expr5)) then
        expr5 = 0d0
      end if

      diff_scattering_rate_q1v1 =  c0**2*(1./(64.*mu**5*(1 + mu)*Sqrt(pi)*v0Nat**4))*(expr1+expr2+expr3+expr4+expr5)

      if (diff_scattering_rate_q1v1.lt.0d0) then
      ! if (isnan(diff_scattering_rate_q1v1)) then
        print*, "nq: ", nq
        print*, "nv: ", nv
        print*, "mdm: ", mdm
        print*, "muplus: ", muplus
        print*, "muminus: ", muminus
        print*, "mu: ", mu
        print*, "w: ", w
        print*, "vesc: ", vesc
        print*, "ue: ", ue
        print*, "expr1: ", expr1
        print*, "expr2: ", expr2
        print*, "expr3: ", expr3
        print*, "expr4: ", expr4
        print*, "expr5: ", expr5
        print*, "Integral of Rminus q1v1: ", diff_scattering_rate_q1v1
      end if

      filename = "debugging_cap_rate/diff_scattering_rate_q1v1_mdm_"//trim(mdmString)//".dat"

      ! inquire(file=filename, exist=exist)
      ! if (exist) then
      !   open(12, file=filename, status="old", position="append", action="write")
      ! else
      !   open(12, file=filename, status="new", action="write")
      !   write(12, *) "DMmass [GeV] | ", "mu | ", "w [cm/s] | ", "vesc [cm/s] | ", "ue [cm/s] | ", &
      !               "expr1 | ", "expr2 | ", "expr3 | ", "Rminus Int | "
      ! end if
      ! write(12, *) mdm, mu, w, vesc, ue, expr1, expr2, expr3, diff_scattering_rate_q1v1
      ! close(12)
    end function

    function diff_scattering_rate_q2v1(w, vesc)
      implicit none
      double precision :: w, diff_scattering_rate_q2v1, vesc
      double precision :: expr1, expr2, expr3, expr4, expr5, expr6, expr7, expr8, expr9
      double precision:: wNat, v0Nat, ueNat, vescNat, x
      logical :: exist
      character(len=300) :: filename
      character ( len = 20 ) :: mdmString
      ue = ue_at_r
      write(mdmString,  '(f10.5)') mdm

      ! w = 72700220.422896266d0
      ! vesc = 71930608.669404954d0
      ! ue = 17864534.274027329d0

      wNat = w/c0
      ueNat = ue/c0
      vescNat = vesc/c0
      v0Nat = v0/c0

      !print*, "INSIDE: diff_scattering_rate_q2v1"
      diff_scattering_rate_q2v1 = 0d0


     !  expr1 = (120*ueNat**4*vescNat)/exp(((1 + mu)*vescNat + (-1 + mu)*wNat)**2/(4*ueNat**2)) + &
     !          (120*ueNat**4*vescNat)/exp((vescNat + mu*vescNat + wNat - mu*wNat)**2/(4*ueNat**2)) + &
     !          (240*ueNat**4*wNat)/exp(((1 + mu)*vescNat + (-1 + mu)*wNat)**2/(4*ueNat**2)) - &
     !          (240*ueNat**4*wNat)/exp((vescNat + mu*vescNat + wNat - mu*wNat)**2/(4*ueNat**2))
     !
     !  expr1 = expr1 + 12*mu*ueNat**2*((vescNat**3 + 9*ueNat**2*(7*vescNat - 12*wNat) - &
     !            2*vescNat**2*wNat + 3*vescNat*wNat**2 - 4*wNat**3)/exp((vescNat + &
     !            mu*vescNat + wNat - mu*wNat)**2/(4*ueNat**2)) + &
     !            (vescNat**3 + 2*vescNat**2*wNat + 3*vescNat*wNat**2 + 4*wNat**3 + &
     !            9*ueNat**2*(7*vescNat + 12*wNat))/exp(((1 + mu)*vescNat + (-1 + mu)*wNat)**2/(4*ueNat**2)))
     !
     !  expr1 =expr1+ 2*mu**2*((6*ueNat**4*(173*vescNat - 262*wNat) + vescNat*(vescNat - wNat)**3*(vescNat + wNat) + &
     !        2*ueNat**2*(25*vescNat**3 - 57*vescNat**2*wNat + 75*vescNat*wNat**2 - &
     !        73*wNat**3))/exp((vescNat + mu*vescNat + wNat - mu*wNat)**2/(4*ueNat**2)) + &
     !        (vescNat*(vescNat - wNat)*(vescNat + wNat)**3 + 6*ueNat**4*(173*vescNat + &
     !        262*wNat) + 2*ueNat**2*(25*vescNat**3 + 57*vescNat**2*wNat + 75*vescNat*wNat**2 + &
     !        73*wNat**3))/exp(((1 + mu)*vescNat + (-1 + mu)*wNat)**2/(4*ueNat**2)))
     !
     !  expr1 = expr1 + 2*mu**3*((2*ueNat**4*(735*vescNat - 869*wNat) + &
     !          ueNat**2*(157*vescNat**3 - 364*vescNat**2*wNat + 403*vescNat*wNat**2 - &
     !          256*wNat**3) + 2*(3*vescNat**5 - 6*vescNat**4*wNat + vescNat**3*wNat**2 + &
     !          2*vescNat**2*wNat**3 + 4*vescNat*wNat**4 - 8*wNat**5))/ &
     !       exp((vescNat + mu*vescNat + wNat - mu*wNat)**2/(4*ueNat**2)) + &
     !       (2*ueNat**4*(735*vescNat + 869*wNat) + ueNat**2*(157*vescNat**3 + 364*vescNat**2*wNat + &
     !       403*vescNat*wNat**2 + 256*wNat**3) + &
     !        2*(3*vescNat**5 + 6*vescNat**4*wNat + vescNat**3*wNat**2 - 2*vescNat**2*wNat**3 + &
     !        4*vescNat*wNat**4 + 8*wNat**5))/exp(((1 + mu)*vescNat + (-1 + mu)*wNat)**2/(4*ueNat**2)))
     !
     !  expr1 =  expr1 + 2*mu**4*((16*vescNat**5 + 2*ueNat**4*(591*vescNat - 629*wNat) - &
     !          36*vescNat**4*wNat + 26*vescNat**3*wNat**2 - 32*vescNat**2*wNat**3 + &
     !          66*vescNat*wNat**4 - 64*wNat**5 + ueNat**2*(265*vescNat**3 - &
     !          648*vescNat**2*wNat + 721*vescNat*wNat**2 - 510*wNat**3))/ &
     !           exp((vescNat + mu*vescNat + wNat - mu*wNat)**2/(4*ueNat**2)) + (16*vescNat**5 + &
     !           36*vescNat**4*wNat + 26*vescNat**3*wNat**2 + 32*vescNat**2*wNat**3 + &
     !           66*vescNat*wNat**4 + 64*wNat**5 + 2*ueNat**4*(591*vescNat + 629*wNat) + &
     !            ueNat**2*(265*vescNat**3 + 648*vescNat**2*wNat + 721*vescNat*wNat**2 + &
     !            510*wNat**3))/exp(((1 + mu)*vescNat + (-1 + mu)*wNat)**2/(4*ueNat**2)))
     !
     ! expr1 = expr1 + mu**5*((51*vescNat**5 + 4*ueNat**4*(273*vescNat - 230*wNat) - &
     !         142*vescNat**4*wNat + 188*vescNat**3*wNat**2 - 242*vescNat**2*wNat**3 + &
     !         241*vescNat*wNat**4 - 64*wNat**5 + 2*ueNat**2*(275*vescNat**3 - &
     !         727*vescNat**2*wNat + 783*vescNat*wNat**2 - 295*wNat**3))/ &
     !       exp((vescNat + mu*vescNat + wNat - mu*wNat)**2/(4*ueNat**2)) + &
     !       (51*vescNat**5 + 142*vescNat**4*wNat + 188*vescNat**3*wNat**2 + &
     !       242*vescNat**2*wNat**3 + 241*vescNat*wNat**4 + 64*wNat**5 + 4*ueNat**4*(273*vescNat + 230*wNat) + &
     !        2*ueNat**2*(275*vescNat**3 + 727*vescNat**2*wNat + 783*vescNat*wNat**2 +&
     !         295*wNat**3))/exp(((1 + mu)*vescNat + (-1 + mu)*wNat)**2/(4*ueNat**2)))
     !
     !
     !  expr1 = expr1 + mu**6*((4*ueNat**4*(57*vescNat - 14*wNat) + ueNat**2*(362*vescNat**3 - &
     !          986*vescNat**2*wNat + 870*vescNat*wNat**2 - 78*wNat**3) + &
     !          vescNat*(55*vescNat**4 - 192*vescNat**3*wNat + &
     !          322*vescNat**2*wNat**2 - 352*vescNat*wNat**3 + 199*wNat**4))/ &
     !          exp((vescNat + mu*vescNat + wNat - mu*wNat)**2/(4*ueNat**2)) + &
     !          (4*ueNat**4*(57*vescNat + 14*wNat) + ueNat**2*(362*vescNat**3 + &
     !          986*vescNat**2*wNat + 870*vescNat*wNat**2 + 78*wNat**3) + &
     !          vescNat*(55*vescNat**4 + 192*vescNat**3*wNat + 322*vescNat**2*wNat**2 + &
     !          352*vescNat*wNat**3 + 199*wNat**4))/exp(((1 + mu)*vescNat + (-1 + mu)*wNat)**2/(4*ueNat**2)))
     !
     !  expr1 = expr1 + 2*mu**7*(((vescNat - wNat)*(ueNat**2*(70*vescNat**2 - 107*vescNat*wNat + &
     !          9*wNat**2) + vescNat*(21*vescNat**3 - 65*vescNat**2*wNat + 87*vescNat*wNat**2 - &
     !          51*wNat**3)))/exp((vescNat + mu*vescNat + wNat - mu*wNat)**2/(4*ueNat**2)) + &
     !          ((vescNat + wNat)*(ueNat**2*(70*vescNat**2 + 107*vescNat*wNat + 9*wNat**2) + &
     !          vescNat*(21*vescNat**3 + 65*vescNat**2*wNat + 87*vescNat*wNat**2 + &
     !          51*wNat**3)))/exp(((1 + mu)*vescNat + (-1 + mu)*wNat)**2/(4*ueNat**2)))
     !
     !  expr1 = expr1 + 2*mu**8*(((vescNat - wNat)**2*(ueNat**2*(12*vescNat - wNat) + &
     !          vescNat*(11*vescNat**2 - 26*vescNat*wNat + 19*wNat**2)))/exp((vescNat &
     !          + mu*vescNat + wNat - mu*wNat)**2/(4*ueNat**2)) + &
     !          ((vescNat + wNat)**2*(ueNat**2*(12*vescNat + wNat) + &
     !        vescNat*(11*vescNat**2 + 26*vescNat*wNat + &
     !        19*wNat**2)))/exp(((1 + mu)*vescNat + (-1 + mu)*wNat)**2/(4*ueNat**2)))
     !
     !  expr1 = expr1 + mu**9*vescNat*(((7*vescNat - 9*wNat)*(vescNat - wNat)**3)/exp((vescNat + &
     !          mu*vescNat + wNat - mu*wNat)**2/(4*ueNat**2)) + &
     !          ((vescNat + wNat)**3*(7*vescNat + 9*wNat))/exp(((1 + &
     !          mu)*vescNat + (-1 + mu)*wNat)**2/(4*ueNat**2)))
     !
     !  expr1 = expr1 + mu**10*vescNat*((vescNat - wNat)**4/exp((vescNat + mu*vescNat + &
     !            wNat - mu*wNat)**2/(4*ueNat**2)) + &
     !            (vescNat + wNat)**4/exp(((1 + mu)*vescNat + (-1 + mu)*wNat)**2/(4*ueNat**2)))
     !
     !  expr1 = (3.*ueNat**3)/(32.*mu**5*(1 + mu)*Sqrt(pi)*v0Nat**6)*expr1
     !

      expr1 = (3.*ueNat**3* &
              (4.*ueNat**4* &
              (3*(exp(-1/4.*((1 + mu)*vescNat + (-1 + mu)*wNat)**2/ &
              ueNat**2) + exp(-1/4.*(vescNat + mu*vescNat + wNat - &
              mu*wNat)**2/ueNat**2))*(1 + mu)* &
              (10 + mu*(53 + mu*(120 + mu*(125 + mu*(72 + 19*mu)))))* &
              vescNat - &
              (-exp(-1/4.*((1 + mu)*vescNat + (-1 + mu)*wNat)**2/ &
              ueNat**2) + exp(-1/4.*(vescNat + mu*vescNat + wNat - &
              mu*wNat)**2/ueNat**2))* &
              (60 + mu*(324 + mu*(786 + mu*(869 + mu*(629 + 2*mu* &
              (115 + 7*mu))))))*wNat) + &
              2*mu*ueNat**2* &
              ((exp(-1/4.*((1 + mu)*vescNat + (-1 + mu)*wNat)**2/ &
              ueNat**2) + exp(-1/4.*(vescNat + mu*vescNat + wNat - &
              mu*wNat)**2/ueNat**2))*(1 + mu)**3* &
              (6 + mu*(32 + mu*(43 + 2*mu*(17 + 6*mu))))*vescNat**3 - &
              (-exp(-1/4.*((1 + mu)*vescNat + (-1 + mu)*wNat)**2/ &
              ueNat**2) + exp(-1/4.*(vescNat + mu*vescNat + wNat - &
              mu*wNat)**2/ueNat**2))*(1 + mu)**2* &
              (12 + mu*(90 + mu*(172 + mu*(214 + mu*(127 + 25*mu)))))* &
              vescNat**2*wNat + &
              (exp(-1/4.*((1 + mu)*vescNat + (-1 + mu)*wNat)**2/ &
              ueNat**2) + exp(-1/4.*(vescNat + mu*vescNat + wNat - &
              mu*wNat)**2/ueNat**2))*(1 + mu)* &
              (18 + mu*(132 + mu*(271 + mu*(450 + mu*(333 + 2*mu* &
              (51 + 7*mu))))))*vescNat*wNat**2 - &
              (-exp(-1/4.*((1 + mu)*vescNat + (-1 + mu)*wNat)**2/ &
              ueNat**2) + exp(-1/4.*(vescNat + mu*vescNat + wNat - &
              mu*wNat)**2/ueNat**2))* &
              (24 + mu*(146 + mu*(256 + mu*(510 + &
              mu*(295 + mu*(39 + mu*(9 + mu)))))))*wNat**3) + &
              mu**2*((exp(-1/4.*((1 + mu)*vescNat + (-1 + mu)*wNat)**2/ &
              ueNat**2) + exp(-1/4.*(vescNat + mu*vescNat + wNat - &
              mu*wNat)**2/ueNat**2))*(1 + mu)**5* &
              (2 + mu*(2 + mu*(2 + mu)))*vescNat**5 - &
              2*(-exp(-1/4.*((1 + mu)*vescNat + (-1 + mu)*wNat)**2/ &
              ueNat**2) + exp(-1/4.*(vescNat + mu*vescNat + wNat - &
              mu*wNat)**2/ueNat**2))*(1 + mu)**4* &
              (2 + mu*(2 + mu)*(2 + mu*(3 + 2*mu)))*vescNat**4*wNat + &
              2*(exp(-1/4.*((1 + mu)*vescNat + (-1 + mu)*wNat)**2/ &
              ueNat**2) + exp(-1/4.*(vescNat + mu*vescNat + wNat - &
              mu*wNat)**2/ueNat**2))*mu*(1 + mu)**3* &
              (2 + mu*(2 + mu)*(10 + 3*mu*(3 + mu)))*vescNat**3* &
              wNat**2 - &
              2*(-exp(-1/4.*((1 + mu)*vescNat + (-1 + mu)*wNat)**2/ &
              ueNat**2) + exp(-1/4.*(vescNat + mu*vescNat + wNat - &
              mu*wNat)**2/ueNat**2))*(1 + mu)**2* &
              (-2 + mu**2*(2 + mu)*(17 + mu*(18 + mu*(9 + 2*mu))))* &
              vescNat**2*wNat**3 + &
              (exp(-1/4.*((1 + mu)*vescNat + (-1 + mu)*wNat)**2/ &
              ueNat**2) + exp(-1/4.*(vescNat + mu*vescNat + wNat - &
              mu*wNat)**2/ueNat**2))*(1 + mu)* &
              (-2 + mu*(18 + mu*(114 + mu*(127 + &
              mu*(72 + mu*(30 + mu*(8 + mu)))))))*vescNat* &
              wNat**4 - &
              32*(-exp(-1/4.*((1 + mu)*vescNat + (-1 + mu)*wNat)**2/ &
              ueNat**2) + exp(-1/4.*(vescNat + mu*vescNat + wNat - &
              mu*wNat)**2/ueNat**2))*mu*(1 + 2*mu*(2 + mu))* &
              wNat**5)))/(32.*mu**5*(1 + mu)*Sqrt(pi)*v0Nat**6)

      expr2 = diff_scattering_rate_q2v1 + &
                                  (-3*ueNat**2*(3*(10 + 53*mu + 120*mu**2 + 125*mu**3 + 72*mu**4 + 19*mu**5)*ueNat**6 + &
                                  3*(10 + 47*mu + 99*mu**2 + 59*mu**3 + 35*mu**4 - 11*mu**5 - 19*mu**6)*ueNat**4*wNat**2 + &
                                  2*mu*(3 + 15*mu + 16*mu**2 + 46*mu**3 - 11*mu**4 - 23*mu**5)*ueNat**2*wNat**4 + &
                                  4*mu**3*(1 + 3*mu - 2*mu**2 - 2*mu**3)*wNat**6)* &
                                  (erf((vescNat + mu*vescNat + wNat - mu*wNat)/(2*ueNat)) + &
                                  erf((vescNat + mu*vescNat - wNat + mu*wNat)/(2*ueNat))))/(8*mu**5*(1 + mu)*v0Nat**6)

      expr3 = diff_scattering_rate_q2v1 + &
                                  ((2*(6*ueNat**7*((60 + 270*mu + 591*mu**2 + 745*mu**3 + 529*mu**4 + 207*mu**5 + &
                                  42*mu**6)*vescNat + &
                                  (60 + 102*mu + 153*mu**2 + mu**3 - 193*mu**4 - 81*mu**5 - 42*mu**6)*wNat) + &
                                  mu*ueNat**5*(-((1 + mu)**3*(72 + 147*mu + 128*mu**2 + 63*mu**3)*vescNat**3) + &
                                  9*(1 + mu)**2*(-8 - 3*mu + 3*mu**2 + mu**3 + 7*mu**4)*vescNat**2*wNat + &
                                  3*(48 + 333*mu + 747*mu**2 + 824*mu**3 + 592*mu**4 + 251*mu**5 +  &
                                  21*mu**6)*vescNat*wNat**2 + (144 + 591*mu + 403*mu**2 - 396*mu**3 - 108*mu**4 - &
                                  571*mu**5 - 63*mu**6)*wNat**3) + &
                                  2*mu**3*ueNat*wNat**2*((1 + mu)**5*vescNat**5 - (-1 + mu)*(1 + mu)**4*vescNat**4*wNat -  &
                                  2*(1 + mu)**3*(1 + 4*mu + mu**2)*vescNat**3*wNat**2 + 2*(1 + mu)**2*(-1 - 3*mu + &
                                  3*mu**2 + mu**3)*vescNat**2*wNat**3 + &
                                  (1 + 9*mu + 38*mu**2 + 38*mu**3 + 9*mu**4 + mu**5)*vescNat*wNat**4 -  &
                                  (-1 - 7*mu - 22*mu**2 + 22*mu**3 + 7*mu**4 + mu**5)*wNat**5) + &
                                  mu**2*ueNat**3*((1 + mu)**5*(6 + 7*mu)*vescNat**5 - (1 + mu)**4*(-6 - mu +  &
                                  7*mu**2)*vescNat**4*wNat - 2*(1 + mu)**3*(15 + 39*mu + 43*mu**2 + 7*mu**3)*vescNat**3*wNat**2 + &
                                  2*(1 + mu)**2*(-15 - 16*mu - 12*mu**2 + 36*mu**3 + 7*mu**4)*vescNat**2*wNat**3 +  &
                                  (24 + 235*mu + 531*mu**2 + 734*mu**3 + 494*mu**4 + 87*mu**5 + 7*mu**6)*vescNat*wNat**4 + &
                                  (24 + 171*mu + 61*mu**2 + 142*mu**3 - 318*mu**4 - 73*mu**5 -  &
                                  7*mu**6)*wNat**5)))/exp(((1 + mu)*vescNat + (-1 + mu)*wNat)**2/(4*ueNat**2)) + &
                                  (2*(6*ueNat**7*((60 + 270*mu + 591*mu**2 + 745*mu**3 + 529*mu**4 +  &
                                  207*mu**5 + 42*mu**6)*vescNat + (-60 - 102*mu - 153*mu**2 - mu**3 + 193*mu**4 + &
                                  81*mu**5 + 42*mu**6)*wNat) + &
                                  mu*ueNat**5*(-((1 + mu)**3*(72 + 147*mu + 128*mu**2 + 63*mu**3)*vescNat**3) -  &
                                  9*(1 + mu)**2*(-8 - 3*mu + 3*mu**2 + mu**3 + 7*mu**4)*vescNat**2*wNat + &
                                  3*(48 + 333*mu + 747*mu**2 + 824*mu**3 + 592*mu**4 + 251*mu**5 +  &
                                  21*mu**6)*vescNat*wNat**2 + (-144 - 591*mu - 403*mu**2 + 396*mu**3 +  &
                                  108*mu**4 + 571*mu**5 + 63*mu**6)*wNat**3) + &
                                  2*mu**3*ueNat*wNat**2*((1 + mu)**5*vescNat**5 + (-1 + mu)*(1 + mu)**4*vescNat**4*wNat -  &
                                  2*(1 + mu)**3*(1 + 4*mu + mu**2)*vescNat**3*wNat**2 - 2*(1 + mu)**2*(-1 - 3*mu + &
                                  3*mu**2 + mu**3)*vescNat**2*wNat**3 + &
                                  (1 + 9*mu + 38*mu**2 + 38*mu**3 + 9*mu**4 + mu**5)*vescNat*wNat**4 + (-1 -  &
                                  7*mu - 22*mu**2 + 22*mu**3 + 7*mu**4 + mu**5)*wNat**5) + &
                                  mu**2*ueNat**3*((1 + mu)**5*(6 + 7*mu)*vescNat**5 + (1 + mu)**4*(-6 - mu +  &
                                  7*mu**2)*vescNat**4*wNat - 2*(1 + mu)**3*(15 + 39*mu + 43*mu**2 + 7*mu**3)*vescNat**3*wNat**2 - &
                                  2*(1 + mu)**2*(-15 - 16*mu - 12*mu**2 + 36*mu**3 + 7*mu**4)*vescNat**2*wNat**3 +  &
                                  (24 + 235*mu + 531*mu**2 + 734*mu**3 + 494*mu**4 + 87*mu**5 + 7*mu**6)*vescNat*wNat**4 + &
                                  (-24 - 171*mu - 61*mu**2 - 142*mu**3 + &
                                  318*mu**4 + 73*mu**5 +  &
                                  7*mu**6)*wNat**5)))/exp((vescNat + mu*vescNat + &
                                  wNat - mu*wNat)**2/(4*ueNat**2)))/(256*mu**6*sqrt(pi)*v0Nat**6)

      expr4 = diff_scattering_rate_q2v1 + &
                                  (-1/256.*((12*(60 + 210*mu + 381*mu**2 + 364*mu**3 + &
                                  165*mu**4 + 42*mu**5)*ueNat**8 - &
                                  12*ueNat**6*(3*(1 + mu)**6*(10 + 7*mu)*vescNat**2 - (30 + 45*mu + &
                                  168*mu**2 + 155*mu**3 + 84*mu**4 + 171*mu**5 + 30*mu**6 + 21*mu**7)*wNat**2) + &
                                  3*mu*ueNat**4*(3*(1 + mu)**6*(8 + 7*mu)*vescNat**4 - &
                                  6*(1 + mu)**6*(12 + 7*mu)*vescNat**2*wNat**2 + (48 + 165*mu + 78*mu**2 - &
                                  53*mu**3 + 372*mu**4 - 101*mu**5 + 174*mu**6 + 21*mu**7)*wNat**4) - &
                                  2*mu**3*wNat**2*((1 + mu)**6*vescNat**6 - &
                                  3*(1 + mu)**6*vescNat**4*wNat**2 + &
                                  3*(1 + mu)**6*vescNat**2*wNat**4 - &
                                  (-1 + mu)**2*(1 + 8*mu + 30*mu**2 + 8*mu**3 + mu**4)*wNat**6) - &
                                  mu**2*ueNat**2*((1 + mu)**6*(6 + 7*mu)*vescNat**6 - &
                                  3*(1 + mu)**6*(12 + 7*mu)*vescNat**4*wNat**2 + 3*(1 + mu)**6*(18 + 7*mu)*vescNat**2*wNat**4 - &
                                  (24 + 151*mu - 78*mu**2 + 201*mu**3 - 428*mu**4 + 249*mu**5 + &
                                  66*mu**6 + 7*mu**7)*wNat**6))*(erf((vescNat + mu*vescNat + wNat - mu*wNat)/(2*ueNat)) +  &
                                  erf((vescNat + mu*vescNat - wNat + mu*wNat)/(2*ueNat))))/(mu**6*v0Nat**6))

      expr5 = diff_scattering_rate_q2v1 + &
                                  (3*(-1 + mu)*ueNat**3*((2*ueNat**4*(-((1 + mu)*(60 + mu*(1 + mu)*(340 + &
                                  mu*(357 + mu*(165 + 28*mu))))*vescNat) - (-1 + mu)*(60 + mu*(388 +  &
                                  mu*(879 + mu*(598 + mu*(207 + 28*mu)))))*wNat) + &
                                  mu*ueNat**2*(-((1 + mu)**3*(24 + mu*(75 + mu*(38 + 7*mu)))*vescNat**3) -  &
                                  (-1 + mu)*(1 + mu)**2*(24 + mu*(91 + mu*(38 + 7*mu)))*vescNat**2*wNat - &
                                  (-1 + mu)**2*(1 + mu)*(24 + mu*(103 + mu*(38 + 7*mu)))*vescNat*wNat**2 -  &
                                  (-1 + mu)**3*(24 + mu*(111 + mu*(38 + 7*mu)))*wNat**3) - &
                                  2*mu**2*((1 + mu)*vescNat + (-1 + mu)*wNat)*((1 + mu)**2*vescNat**2 + (-1 +  &
                                  mu**2)*vescNat*wNat + (-1 + mu)**2*wNat**2)*((1 + mu)**2*vescNat**2 + &
                                  (-1 + mu)**2*wNat**2 + vescNat*(wNat - mu**2*wNat)))/ &
                                  exp((vescNat + mu*vescNat + wNat - mu*wNat)**2/(4*ueNat**2)) + &
                                  (2*ueNat**4*(-((1 + mu)*(60 + mu*(1 + mu)*(340 + mu*(357 + mu*(165 +  &
                                  28*mu))))*vescNat) + (-1 + mu)*(60 + mu*(388 + mu*(879 + mu*(598 + mu*(207 + 28*mu)))))*wNat) + &
                                  mu*ueNat**2*(-((1 + mu)**3*(24 + mu*(75 + mu*(38 + 7*mu)))*vescNat**3) +  &
                                  (-1 + mu)*(1 + mu)**2*(24 + mu*(91 + mu*(38 + 7*mu)))*vescNat**2*wNat - &
                                  (-1 + mu)**2*(1 + mu)*(24 + mu*(103 + mu*(38 + 7*mu)))*vescNat*wNat**2 +  &
                                  (-1 + mu)**3*(24 + mu*(111 + mu*(38 + 7*mu)))*wNat**3) - &
                                  2*mu**2*(vescNat + mu*vescNat + wNat - mu*wNat)*((1 + mu)**2*vescNat**2 +  &
                                  (-1 + mu**2)*vescNat*wNat + (-1 + mu)**2*wNat**2)*((1 + mu)**2*vescNat**2 +  &
                                  (-1 + mu)**2*wNat**2 + vescNat*(wNat - mu**2*wNat)))/ &
                                  exp(((1 + mu)*vescNat + (-1 + mu)*wNat)**2/(4*ueNat**2))))/(128*mu**6*(1 + mu)*sqrt(pi)*v0Nat**6)

      expr6 = diff_scattering_rate_q2v1 + &
                                  (3*(-1 + mu)*ueNat**2*(4*(60 + 340*mu + 697*mu**2 + 522*mu**3 + &
                                  193*mu**4 + 28*mu**5)*ueNat**6 + 4*(-1 + mu)**2*(30 + 206*mu + 491*mu**2 + &
                                  318*mu**3 + 107*mu**4 + 14*mu**5)*ueNat**4*wNat**2 + &
                                  (-1 + mu)**4*mu*(24 + 115*mu + 38*mu**2 + 7*mu**3)*ueNat**2*wNat**4 + &
                                  2*(-1 + mu)**6*mu**2*wNat**6)*(erf((vescNat + mu*vescNat + wNat - mu*wNat)/(2*ueNat)) +  &
                                  erf((vescNat + mu*vescNat - wNat + mu*wNat)/(2*ueNat))))/ &
                                  (256*mu**6*(1 + mu)*v0Nat**6)
      x = exp((mu*(-vescNat**2 + wNat**2))/ueNat**2)*&
          (-erf(((-1 + mu)*vescNat + wNat + mu*wNat)/(2*ueNat)) &
          + erf((vescNat - mu*vescNat +  &
          wNat + mu*wNat)/(2*ueNat)))

      if (isnan(x)) then
        x = 0d0
      end if

      expr7 = diff_scattering_rate_q2v1 + &
                                  (3*(1 + mu)**5*ueNat**2*(2*(-1 + mu)*ueNat**4*(10*(12 + 7*mu)*ueNat**2 - &
                                  3*mu*(8 + 7*mu)*wNat**2)*(erf((vescNat + mu*vescNat + wNat - mu*wNat)/(2*ueNat)) +  &
                                  erf((vescNat + mu*vescNat - wNat + mu*wNat)/(2*ueNat))) +  &
                                  x*(1 + mu)*(20*(12 + 7*mu)*ueNat**6 +  &
                                  mu**2*(24 + 7*mu)*ueNat**2*vescNat**4 + 2*mu**3*vescNat**6 &
                                  + 2*mu*ueNat**4*(4*(15 + 7*mu)*vescNat**2  &
                                  - 3*(8 + 7*mu)*wNat**2))))/(256*mu**7*v0Nat**6)
      expr8 = diff_scattering_rate_q2v1 + &
                                  (-3*(1 - mu)*(1 + mu)*ueNat**2*wNat**2*(-(sqrt(pi)*(2*(16 +  &
                                  mu*(51 + mu*(30 + 7*mu)))*ueNat**4 + &
                                  (16 + mu*(29 + mu*(-84 + mu*(8 + mu*(16 + 7*mu)))))*ueNat**2*wNat**2 + (-1 +  &
                                  mu)**2*mu*(1 + (-6 + mu)*mu)*wNat**4)*erf((vescNat + mu*vescNat + wNat -  &
                                  mu*wNat)/(2*ueNat))) - &
                                  sqrt(pi)*(2*(16 + mu*(51 + mu*(30 + 7*mu)))*ueNat**4 + &
                                  (16 + mu*(29 + mu*(-84 + mu*(8 + mu*(16 + 7*mu)))))*ueNat**2*wNat**2 +  &
                                  (-1 + mu)**2*mu*(1 + (-6 + mu)*mu)*wNat**4)*erf((vescNat +&
                                   mu*vescNat - wNat + mu*wNat)/(2*ueNat)) + &
                                  4*ueNat*((1 + mu)*vescNat*((16 +  &
                                  mu*(51 + mu*(30 + 7*mu)))*ueNat**2 + 2*mu*(1 + &
                                  mu)**2*vescNat**2 + mu*(1 + (-6 + mu)*mu)*wNat**2)* &
                                  (+(1/2.)*exp(-(((-1 + mu**2)*vescNat*wNat)/(&
                                  2*ueNat**2)) - ((1 + mu)**2*vescNat**2 + &
                                  (-1 + mu)**2*wNat**2)/(4*ueNat**2)) + &
                                   1/2.*exp(((-1 + mu**2)*vescNat*wNat)/( &
                                  2*ueNat**2) - ((1 + mu)**2*vescNat**2 + &
                                  (-1 + mu)**2*wNat**2)/(4*ueNat**2)))+ (-1 + mu)*wNat*((16 + mu*(59 + mu*(30 +  &
                                  7*mu)))*ueNat**2 + 2*mu*(1 + mu)**2*vescNat**2 + mu*(1 + (-6 + mu)*mu)*wNat**2)* &
                                  (-(1/2.)*exp(-(((-1 + mu**2)*vescNat*wNat)/(&
                                  2*ueNat**2)) - ((1 + mu)**2*vescNat**2 + &
                                  (-1 + mu)**2*wNat**2)/(4*ueNat**2)) + &
                                   1/2.*exp(((-1 + mu**2)*vescNat*wNat)/( &
                                  2*ueNat**2) - ((1 + mu)**2*vescNat**2 + &
                                  (-1 + mu)**2*wNat**2)/(4*ueNat**2))))))/(128*mu**5*sqrt(pi)*v0Nat**6)

      expr9 = diff_scattering_rate_q2v1 + &
                                  (3*(1 + mu)**5*ueNat**2*wNat**2*((-1 +  &
                                  mu)*ueNat**2*(-2*(16 + 7*mu)*ueNat**2 +&
                                   mu*(8 + 7*mu)*wNat**2)*erf((vescNat + mu*vescNat + wNat - &
                                   mu*wNat)/(2*ueNat)) + &
                                  (-1 + mu)*ueNat**2*(-2*(16 + 7*mu)*ueNat**2 +  &
                                  mu*(8 + 7*mu)*wNat**2)*erf((vescNat + mu*vescNat - wNat + mu*wNat)/(2*ueNat)) - &
                                  x*(1 + &
                                  mu)*(2*(16 + 7*mu)*ueNat**4 + 2*mu**2*vescNat**2*(2*vescNat**2 -  &
                                  wNat**2) + mu*ueNat**2*(2*(16 + 7*mu)*vescNat**2 - (8 + 7*mu)*wNat**2))))/(256*mu**6*v0Nat**6)
      diff_scattering_rate_q2v1 = c0**2*(expr1 + expr2 + expr3 + expr4 + &
                                  expr5 + expr6 + expr7 + expr8 + expr9)

        ! print*, "nq: ", nq
        ! print*, "nv: ", nv
        ! print*, "mdm: ", mdm
        ! print*, "muplus: ", muplus
        ! print*, "muminus: ", muminus
        ! print*, "mu: ", mu
        ! print*, "w: ", w
        ! print*, "vesc: ", vesc
        ! print*, "ue: ", ue
        ! print*, "expr1: ", c0**2*expr1
        ! print*, "expr2: ", c0**2*expr2
        ! print*, "expr3: ", c0**2*expr3
        ! print*, "expr4: ", c0**2*expr4
        ! print*, "expr5: ", c0**2*expr5
        ! print*, "expr6: ", c0**2*expr6
        ! print*, "expr7: ", c0**2*expr7
        ! print*, "expr8: ", c0**2*expr8
        ! print*, "expr9: ", c0**2*expr9
        ! print*, "Integral of Rminus q2v1: ", diff_scattering_rate_q2v1

      if (diff_scattering_rate_q2v1.lt.0d0) then
      ! if (isnan(diff_scattering_rate_q2v1)) then
        print*, "nq: ", nq
        print*, "nv: ", nv
        print*, "mdm: ", mdm
        print*, "muplus: ", muplus
        print*, "muminus: ", muminus
        print*, "mu: ", mu
        print*, "w: ", w
        print*, "vesc: ", vesc
        print*, "ue: ", ue
        print*, "expr1: ", c0**2*expr1
        print*, "expr2: ", c0**2*expr2
        print*, "expr3: ", c0**2*expr3
        print*, "expr4: ", c0**2*expr4
        print*, "expr5: ", c0**2*expr5
        print*, "expr6: ", c0**2*expr6
        print*, "expr7: ", c0**2*expr7
        print*, "expr8: ", c0**2*expr8
        print*, "expr9: ", c0**2*expr9
        print*, "Integral of Rminus q2v1: ", diff_scattering_rate_q2v1
      end if

      ! filename = "debugging_cap_rate/diff_scattering_rate_q2v1_mdm_"//trim(mdmString)//".dat"
      ! ! if (diff_scattering_rate_q2v1.lt.0d0) then
      !   inquire(file=filename, exist=exist)
      !   if (exist) then
      !     open(12, file=filename, status="old", position="append", action="write")
      !   else
      !     open(12, file=filename, status="new", action="write")
      !     write(12, *) "DMmass [GeV] | ", "mu | ", "w [cm/s] | ", "vesc [cm/s] | ", "ue [cm/s] | ", &
      !                 "expr1 | ", "expr2 | ", "expr3 | ", "expr4 | ", "expr5 | ", "expr6 | ", &
      !                 "expr7 | ", "expr8 | ", "expr9 | ", &
      !                 "Rminus Int | "
      !   end if
      !   write(12, *) mdm, mu, w, vesc, ue, c0**2*expr1, c0**2*expr2, c0**2*expr3, &
      !               c0**2*expr4, c0**2*expr5, c0**2*expr6, c0**2*expr7, c0**2*expr8, &
                      ! c0**2*expr9, &
      !               diff_scattering_rate_q2v1
      !   close(12)
      ! ! end if
    end function

    !SB: I do not know why this is necessary but it seems to be.
    function reset_Rminus(mdm_in, mtarget, ue_at_r_in, vesc_in)
      double precision :: mdm_in, mtarget, ue_at_r_in, vesc_in
      double precision :: reset_Rminus
      mdm = mdm_in
      mu = mdm/mtarget
      muminus = (mu-1.)/2.
      muplus = (mu+1.)/2.
      ue_at_r = ue_at_r_in
      vesc_shared = vesc_in
      reset_Rminus = 0d0
    end function

!----------------------------------------------------------------------------------------------------
      !Fast trapezoidal integral
      function trapz(x,y,flen)
        implicit none
        integer, intent(in) :: flen
        double precision, intent (in) :: x(flen), y(flen)
        double precision trapz
        integer i
        trapz = y(1)*(x(2)-x(1))/2. + y(flen)*(x(flen)-x(flen-1))/2.
        do i = 2,flen-1
          trapz = trapz + y(i)*(x(i)-x(i-1))
        end do
        return
      end function

    end module capmod



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    ! Some functions that have to be external, because of the integrator.


    !The integrand for the integral over u
    function integrand(u,foveru)
      use capmod
      double precision :: u, w, integrand, foveru
      external foveru
      w = sqrt(u**2+vesc_shared**2)

      !Switch depending on whether we are capturing on Hydrogen or not
      if (a_shared .gt. 2.d0) then
        integrand = foveru(u)*GFFI_A(w,vesc_shared,a_shared)
      else
        integrand = foveru(u)*GFFI_H(w,vesc_shared)
      end if

      !Rescale for velocity-dependent cross-sections
      if (nv .ne. 0) then
        integrand = integrand*(w/v0)**(2*nv)
      end if
    end function integrand

    !SB: integrand for electron capture that uses the Rminus functions
    function integrand_Rminus(u,foveru)
      use capmod
      implicit none
      double precision :: u, w, integrand_Rminus, foveru
      external foveru
      w = sqrt(u**2+vesc_shared**2)
      !w_for_all = w
      if((nq.eq.0).and.(nv.eq.0)) then
        integrand_Rminus = foveru(u)*diff_scattering_rate_cst(w, vesc_shared)
      else if((nq.eq.1).and.(nv.eq.1)) then
        integrand_Rminus = foveru(u)*diff_scattering_rate_q1v1(w, vesc_shared)
      else if((nq.eq.2).and.(nv.eq.1)) then
        integrand_Rminus = foveru(u)*diff_scattering_rate_q2v1(w, vesc_shared)
      else if (nv.eq.1) then
         integrand_Rminus = foveru(u)*diff_scattering_rate_v1(w, vesc_shared)
      else if (nq.eq.1) then
         integrand_Rminus = foveru(u)*diff_scattering_rate_q1(w, vesc_shared)
      else if (nv.eq.2) then
        integrand_Rminus = foveru(u)*diff_scattering_rate_v2(w, vesc_shared)
      else if (nq.eq.2) then
         integrand_Rminus = foveru(u)*diff_scattering_rate_minus_q2(w, vesc_shared)
       else if (nq.eq.3) then
          integrand_Rminus = foveru(u)*diff_scattering_rate_q3(w, vesc_shared)
      else if (nv.eq.-1) then
          integrand_Rminus = foveru(u)*diff_scattering_rate_vminus1(w, vesc_shared)
      else
        integrand_Rminus = 0
      end if
    end function integrand_Rminus

    subroutine captn_general(mx_in,sigma_0,niso,nq_in,nv_in,spin_in,capped)
      use capmod
      implicit none
      integer, intent(in):: nq_in, nv_in, niso, spin_in
      ! integer, intent(in):: spin_in
      integer eli, ri, limit, j
      double precision, intent(in) :: mx_in, sigma_0
      double precision :: capped !this is the output
      double precision :: sigma_SD, sigma_SI
      double precision :: maxcap, maxcapped, a, sigma_N, umax, umin, vesc, Hn
      double precision :: epsabs, epsrel, abserr, neval  !for integrator
      double precision :: ier,alist,blist,rlist,elist,iord,last!for integrator
      double precision :: int_result, ur, w
      !Steph: testing agreement with NREO
      double precision :: coupling, mreduced, cappedTest

      dimension alist(1000),blist(1000),elist(1000),iord(1000),   rlist(1000)!for integrator
      external integrand
      external gausstest !this is just for testing

      epsabs=1.d-8
      epsrel=1.d-8
      limit=1000

      mdm = mx_in
      nq = nq_in
      nv = nv_in
      coupling = 1d-3/246.2**2
      if (spin_in == 1) then
        sigma_SD = sigma_0
        sigma_SI = 0.d0
      else if (spin_in == 0) then
        sigma_SD = 0.d0
        sigma_SI = sigma_0
      end if

      if (nq*nv .ne. 0) then
        stop "Oh no! nq and nv can't both be nonzero."
      end if

      if (.not. allocated(tab_r)) then
        stop "You haven't yet called captn_init to load the solar model!"
      end if

      capped = 0.d0
      cappedTest = 0.d0
      ! open(105,file = "debug_cap_rate.dat")
      ! write(105,*) "DM mass [GeV]: ", mx_in
      ! write(105,*)
      ! write(105,*) "Radius/RSolar | ", "Temperature [K] | ", "Escape Velocity [cm/s] | ", "u(r) [cm/s]"
      !Loop over the shells of constant radius in the star
      do ri = 1, nlines

        vesc = tab_vesc(ri)
        vesc_shared = vesc !make accessible via the module

        !Loop over the different elements
        do eli = 1, niso

          a = AtomicNumber(eli)
          a_shared = a !make accessible via the module
          mreduced = mdm*(a*mnuc)/(mdm+a*mnuc)
          !This is fine for SD as long as it's just hydrogen. Otherwise, spins must be added.
          sigma_N = (a**2 * (sigma_SI*a**2 + sigma_SD) * (mx_in+mnuc)**2/(mx_in+a*mnuc)**2)
          mu = mx_in/(mnuc*a)
          !sigma_N = sigma_SI
          muplus = (1.+mu)/2.
          muminus = (mu-1.d0)/2.
          ! Hn = (1+dble(nq))*muplus**(2*dble(nq))
          Hn = 1d0
          ! q0 = mdm*v0/c0
          q0 = 0.04

          !sigma_N = q0**2*1./4.*coupling**2/(a*mnuc)**2/pi*mreduced**2*(2d-14)**2/4.

          ! Bottom part of the integral is always zero -- happy little slow DM particles can always be captured.
          umin = 0.d0
          ! Chop the top of the integral off at the smaller of the halo escape velocity or the minimum velocity required for capture.
          umax = min(vesc * sqrt(mu)/abs(muminus), vesc_halo)

          ur = sqrt(2*kb_here*tab_T(ri)/(mnuc*1.78d-27))*100
          ! write(105, *) tab_r(ri), tab_T(ri), tab_vesc(ri), ur

          !Call integrator
          call dsntdqagse(integrand,vdist_over_u,umin,umax, &
          epsabs,epsrel,limit,int_result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
          int_result = Hn*int_result * 2.d0 * sigma_N * NAvo * tab_starrho(ri)*tab_mfr(ri,eli) * (muplus/mx_in)**2
          capped = capped + tab_r(ri)**2*int_result*tab_dr(ri)

          cappedTest = cappedTest + tab_r(ri)**2*tab_dr(ri)*(w_func_target(umax, mnuc, vesc)-w_func_target(0.d0, mnuc, vesc))&
                  *sigma_0*tab_mfr(ri,1)*NAvo*tab_starrho(ri)/mnuc/mx_in/u0/usun*sqrt(3./2./pi)

          if (isnan(capped)) then
            capped = 0.d0
            stop 'NaN encountered whilst trying compute capture rate.'
          end if

        end do

      end do
      capped = 4.d0*pi*Rsun**3*capped
      cappedTest = 4.d0*pi*Rsun**3*cappedTest
      ! print*, "*********************Testing*********************"
      ! print*, "GFFI: ", capped
      ! print*, "wFunc: ", cappedTest
      ! print*, "Ratio: ", cappedTest/capped

      ! open(107,file = "GFFI_H.dat")
      ! write(107,*) "DM mass [GeV]: ", mx_in
      ! write(107,*) "Escape Velocity: ", tab_vesc(1)
      ! write(107,*)
      ! write(107,*) "w [cm/s] | ", "GFFI | "
      ! w = tab_vesc(1)
      ! do j=1,1000
      !   w = w+1000000.d0
      !   write(107,*) w , GFFI_H(w, tab_vesc(1))
      ! end do
      ! close(107)

      if (capped .gt. 1.d100) then
        print*,"Capt'n General says: Oh my, it looks like you are capturing an"
        print*,"infinite amount of dark matter in the Sun. Best to look into that."
      end if
      maxcapped = maxcap(mx_in)
      ! if (capped .gt. maxcapped) then
      !   capped = maxcapped
      ! end if
    end subroutine captn_general


    subroutine captn_general_electron_GFFI(mx_in,sigma_0,niso,nq_in,nv_in,spin_in,capped)
      use capmod
      implicit none
      integer, intent(in):: nq_in, nv_in, niso, spin_in
      ! integer, intent(in):: spin_in
      integer eli, ri, limit, j
      double precision, intent(in) :: mx_in, sigma_0
      double precision :: capped !this is the output
      double precision :: sigma_SD, sigma_SI
      double precision :: maxcap, maxcapped, a, sigma_N, umax, umin, vesc, Hn
      double precision :: epsabs, epsrel, abserr, neval  !for integrator
      double precision :: ier,alist,blist,rlist,elist,iord,last!for integrator
      double precision :: int_result, ur, w
      !Steph: testing agreement with NREO
      double precision :: coupling, mreduced, cappedTest

      dimension alist(1000),blist(1000),elist(1000),iord(1000),   rlist(1000)!for integrator
      external integrand
      external gausstest !this is just for testing

      epsabs=1.d-8
      epsrel=1.d-8
      limit=1000

      mdm = mx_in
      nq = nq_in
      nv = nv_in

      if (spin_in == 1) then
        sigma_SD = sigma_0
        sigma_SI = 0.d0
      else if (spin_in == 0) then
        sigma_SD = 0.d0
        sigma_SI = sigma_0
      end if

      if (nq*nv .ne. 0) then
        stop "Oh no! nq and nv can't both be nonzero."
      end if

      if (.not. allocated(tab_r)) then
        stop "You haven't yet called captn_init to load the solar model!"
      end if

      capped = 0.d0

      !Loop over the shells of constant radius in the star
      do ri = 1, nlines

        vesc = tab_vesc(ri)
        vesc_shared = vesc !make accessible via the module

        !Loop over the different elements
        do eli = 1, niso

          a = AtomicNumber(eli)
          a_shared = 1 !make accessible via the module
          mreduced = mdm*(melectron)/(mdm+melectron)
          !This is fine for SD as long as it's just hydrogen. Otherwise, spins must be added.
          sigma_N = sigma_0
          mu = mx_in/(melectron)
          muplus = (1.+mu)/2.
          muminus = (mu-1.d0)/2.
          Hn = (1+dble(nq))*muplus**(2*dble(nq))
          ! Hn = 1
          q0 = mdm*v0/c0
          ! q0 = 0.04

          !sigma_N = q0**2*1./4.*coupling**2/(a*mnuc)**2/pi*mreduced**2*(2d-14)**2/4.

          ! Bottom part of the integral is always zero -- happy little slow DM particles can always be captured.
          umin = 0.d0
          ! Chop the top of the integral off at the smaller of the halo escape velocity or the minimum velocity required for capture.
          umax = min(vesc * sqrt(mu)/abs(muminus), vesc_halo)

          ur = sqrt(2*kb_here*tab_T(ri)/(melectron*1.78d-27))*100
          ! write(105, *) tab_r(ri), tab_T(ri), tab_vesc(ri), ur

          !Call integrator
          call dsntdqagse(integrand,vdist_over_u,umin,umax, &
          epsabs,epsrel,limit,int_result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
          int_result = Hn*int_result * 2.d0 * sigma_N * NAvo * tab_starrho(ri)*tab_electron_mfr(ri) * (muplus/mx_in)**2
          capped = capped + tab_r(ri)**2*int_result*tab_dr(ri)

          if (isnan(capped)) then
            capped = 0.d0
            stop 'NaN encountered whilst trying compute capture rate.'
          end if

        end do

      end do
      capped = 4.d0*pi*Rsun**3*capped
      cappedTest = 4.d0*pi*Rsun**3*cappedTest
      ! print*, "*********************Testing*********************"
      ! print*, "GFFI: ", capped
      ! print*, "wFunc: ", cappedTest
      ! print*, "Ratio: ", cappedTest/capped

      ! open(107,file = "GFFI_H.dat")
      ! write(107,*) "DM mass [GeV]: ", mx_in
      ! write(107,*) "Escape Velocity: ", tab_vesc(1)
      ! write(107,*)
      ! write(107,*) "w [cm/s] | ", "GFFI | "
      ! w = tab_vesc(1)
      ! do j=1,1000
      !   w = w+1000000.d0
      !   write(107,*) w , GFFI_H(w, tab_vesc(1))
      ! end do
      ! close(107)

      if (capped .gt. 1.d100) then
        print*,"Capt'n General says: Oh my, it looks like you are capturing an"
        print*,"infinite amount of dark matter in the Sun. Best to look into that."
      end if

      maxcapped = maxcap(mx_in)
      if (capped .gt. maxcapped) then
        capped = maxcapped
      end if
    end subroutine captn_general_electron_GFFI


    subroutine captn_general_electron(mx_in,sigma_0,nq_in,nv_in,capped)
      use capmod
      implicit none
      integer, intent(in):: nq_in, nv_in
      ! integer, intent(in):: spin_in
      integer eli, ri, limit, i
      double precision, intent(in) :: mx_in, sigma_0
      double precision :: capped !this is the output
      double precision :: Hn, w
      double precision :: maxcap, maxcapped, a, sigma_N, umax, umin, vesc
      double precision :: epsabs, epsrel, abserr, neval  !for integrator
      double precision :: ier,alist,blist,rlist,elist,iord,last!for integrator
      double precision :: int_result, x, cappedTest

      dimension alist(1000),blist(1000),elist(1000),iord(1000),   rlist(1000)!for integrator
      external integrand
      external integrand_Rminus
      external gausstest !this is just for testing

      epsabs=1.d-8
      epsrel=1.d-8
      limit=1000

      mdm = mx_in
      nq = nq_in
      nv = nv_in
      a_shared = 1d0
      if (nq*nv .ne. 0) then
        stop "Oh no! nq and nv can't both be nonzero."
      end if

      if (.not. allocated(tab_r)) then
        stop "You haven't yet called captn_init to load the solar model!"
      end if

      capped = 0.d0
      cappedTest = 0.d0

      !Loop over the shells of constant radius in the star
      open(98,file = "/home/sberam/Documents/Masters/org-captngen-mine/Jupyter_Notebooks/debugging_cap_rate/Paper_v_captnGen/&
          capped_captnGen.dat")
      write(98,*) "Radius/RSolar| ", "Cumulative Capped| ", "Integral Result"
      do ri = 1, nlines
        ue_at_r = sqrt(2*kb_here*tab_T(ri)/melectronKg)*100 !in cm/s like the other velocities
        vesc = tab_vesc(ri)
        vesc_shared = vesc !make accessible via the module

        mu = mx_in/(melectron)
        muplus = (1.+mu)/2.
        muminus = (mu-1.d0)/2.
        !x = t_integral(133757026d0, vesc/2.)
        !x = diff_scattering_rate_minus_vm2(133757026d0, vesc)
        !q0 = mdm*v0/c0*sqrt((1+mu)**2/2)**(-1)
        Hn = (1+dble(nq))*muplus**(2*dble(nq))
        ! Hn = 1d0
        q0 = mdm*v0/c0
        ! q0 = 0.04

        ! Bottom part of the integral is always zero -- happy little slow DM particles can always be captured.
        umin = 0.d0
        ! Chop the top of the integral off at the smaller of the halo escape velocity or the minimum velocity required for capture.
        umax = min(vesc * sqrt(mu)/abs(muminus), vesc_halo)
        !Call integrator
        call dsntdqagse(integrand_Rminus,vdist_over_u,umin,umax, &
        epsabs,epsrel,limit,int_result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
        int_result = Hn*int_result * sigma_0 * NAvo * tab_starrho(ri)*tab_electron_mfr(ri)/melectron
        capped = capped + tab_r(ri)**2*int_result*tab_dr(ri)
        cappedTest = cappedTest + tab_r(ri)**2*tab_dr(ri)*(w_func_target(umax, melectron, vesc)-&
                      w_func_target(0.d0, melectron, vesc))&
                *sigma_0*tab_electron_mfr(ri)*NAvo*tab_starrho(ri)/melectron/mx_in/u0/usun*sqrt(3./2./pi)
        write(98,*) tab_r(ri), capped, int_result

        if (isnan(capped)) then
          capped = 0.d0
          stop 'NaN encountered whilst trying compute capture rate.'
        end if
      end do

      capped = 4.d0*pi*Rsun**3*capped
      cappedTest = 4.d0*pi*Rsun**3*cappedTest
      if (capped .gt. 1.d100) then
        print*,"Capt'n General says: Oh my, it looks like you are capturing an"
        print*,"infinite amount of dark matter in the Sun. Best to look into that."
      end if

      maxcapped = maxcap(mx_in)
      ! if (capped .gt. maxcapped) then
      !   capped = maxcapped
      ! end if
      close(98)
    end subroutine captn_general_electron

    !SB: capture rate using Rminus method for electrons or hydrogen ONLY
    subroutine captn_general_Rminus(mx_in,sigma_0,nq_in,nv_in,electron_v_nucleons, capped)
      use capmod
      implicit none
      integer, intent(in):: nq_in, nv_in, electron_v_nucleons
      ! integer, intent(in):: spin_in
      integer eli, ri, limit, i
      double precision, intent(in) :: mx_in, sigma_0
      double precision :: capped !this is the output
      double precision :: Hn, w
      double precision :: maxcap, maxcapped, a, sigma_N, umax, umin, vesc
      double precision :: epsabs, epsrel, abserr, neval  !for integrator
      double precision :: ier,alist,blist,rlist,elist,iord,last!for integrator
      double precision :: int_result, x, cappedTest, mtarget, mfr(nlines), mtargetKg

      dimension alist(1000),blist(1000),elist(1000),iord(1000),   rlist(1000)!for integrator
      external integrand
      external integrand_Rminus
      external gausstest !this is just for testing

      epsabs=1.d-8
      epsrel=1.d-8
      limit=1000

      if (electron_v_nucleons.eq.0) then
        mtarget = melectron
        mfr = tab_electron_mfr
      else
        mtarget = mnuc
        mfr = tab_mfr(:,1)
      end if
      mtargetKg = mtarget*1.782662d-27 ![kg]
      mdm = mx_in
      nq = nq_in
      nv = nv_in
      ! nq = 0
      ! nv = -1
      a_shared = 1d0
      ! if (nq*nv .ne. 0) then
      !   stop "Oh no! nq and nv can't both be nonzero."
      ! end if

      if (.not. allocated(tab_r)) then
        stop "You haven't yet called captn_init to load the solar model!"
      end if

      capped = 0.d0
      cappedTest = 0.d0

      ! ue_at_r = 50523957.665483847d0
      ! mu = mx_in/(mtarget)
      ! x = diff_scattering_rate_vminus1(138276456.67132345d0, 138152250.50688523d0)
      !Loop over the shells of constant radius in the star
      do ri = 1, nlines
        ue_at_r = sqrt(2*kb_here*tab_T(ri)/mtargetKg)*100 !in cm/s like the other velocities
        vesc = tab_vesc(ri)
        vesc_shared = vesc !make accessible via the module

        mu = mx_in/(mtarget)
        muplus = (1.+mu)/2.
        muminus = (mu-1.d0)/2.
        Hn = (1.+nq)*muplus**(2.*nq)
        !q0 = mdm*v0/c0
        q0 = 0.04

        ! Bottom part of the integral is always zero -- happy little slow DM particles can always be captured.
        umin = 0.d0
        ! Chop the top of the integral off at the smaller of the halo escape velocity or the minimum velocity required for capture.
        umax = min(vesc * sqrt(mu)/abs(muminus), vesc_halo)

        !Call integrator

        call dsntdqagse(integrand_Rminus,vdist_over_u,umin,umax, &
        epsabs,epsrel,limit,int_result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)

        int_result = int_result * sigma_0 * NAvo * tab_starrho(ri)*mfr(ri)/mtarget
        if(isnan(int_result)) then
          int_result = 0d0
        end if
        capped = capped + tab_r(ri)**2*int_result*tab_dr(ri)
        ! print*, "ri: ", ri
      end do

      capped = 1./(Hn)*(mdm*v0/c0/q0)**(2.*nq)*&
              4.d0*pi*Rsun**3*capped !!!!!!!!!!!!CORRETING MAKING q0 = 0.04 and removing Hn factor

      if (capped .gt. 1.d100) then
        print*,"Capt'n General says: Oh my, it looks like you are capturing an"
        print*,"infinite amount of dark matter in the Sun. Best to look into that."
      end if
      maxcapped = maxcap(mx_in)
      if (capped .gt. maxcapped) then
        capped = maxcapped
      end if
  end subroutine captn_general_Rminus


  ! subroutine evap_RPlus(mx_in,sigma_0,Nwimps,nq_in,nv_in,electron_v_nucleons, velFactor, evapRate)
  !   use capmod
  !   use spergelpressmod
  !   implicit none
  !   integer, intent(in):: nq_in, nv_in, electron_v_nucleons
  !   ! integer, intent(in):: spin_in
  !   integer eli, ri, limit, i
  !   double precision, intent(in) :: mx_in, sigma_0, Nwimps, velFactor
  !   double precision :: evapRate !this is the output
  !   double precision :: Tx
  !   double precision :: a, sigma_N, wmax, wmin, vesc
  !   double precision :: epsabs, epsrel, abserr, neval  !for integrator
  !   double precision :: ier,alist,blist,rlist,elist,iord,last!for integrator
  !   double precision :: int_result, x, mtarget, mfr(nlines), mtargetKg, nChi(nlines), mfp(nlines)
  !   double precision :: zeta_q(nlines), zeta_v(nlines)
  !   double precision :: tau_r(nlines), eta_ang(nlines), eta_multi(nlines), phi_r(nlines), supp_r(nlines), nabund(nlines)
  !   double precision :: mdmKg, rchi, mdm_g, q0_cgs, targetmass_g
  !   double precision :: guess_1, guess_2, reltolerance
  !   double precision:: GN = 6.674d-8
  !
  !   dimension alist(1000),blist(1000),elist(1000),iord(1000),   rlist(1000)!for integrator
  !   external integrand_RPlus
  !   external gausstest !this is just for testing
  !
  !   epsabs=1.d-8
  !   epsrel=1.d-8
  !   limit=1000
  !
  !   if (electron_v_nucleons.eq.0) then
  !     mtarget = melectron
  !     mfr = tab_electron_mfr
  !   else
  !     mtarget = mnuc
  !     mfr = tab_mfr(:,1)
  !   end if
  !   mtargetKg = mtarget*1.782662d-27 ![kg]
  !   targetmass_g = mtargetKg*1d3
  !   mdm = mx_in
  !   mdmKg = mx_in*1.782662d-27
  !   mdm_g = mdmKg*1d3
  !   nq = nq_in
  !   nv = nv_in
  !   q0_cgs = q0*5.344d-14
  !   a_shared = 1d0
  !   nabund = mfr*tab_starrho/(mtargetKg*1d3)
  !   ! if (nq*nv .ne. 0) then
  !   !   stop "Oh no! nq and nv can't both be nonzero."
  !   ! end if
  !
  !   if (.not. allocated(tab_r)) then
  !     stop "You haven't yet called captn_init to load the solar model!"
  !   end if
  !
  !   evapRate = 0.d0
  !
  !   guess_1 = maxval(tab_T)*15d0 ! One-zone WIMP temp guesses in K.
  !   guess_2 = maxval(tab_T)/150.d0
  !   reltolerance = 1.0d-8
  !
  !   Tx = binary_search_mine(Tx_integral_mine, sigma_0, mtargetKg*1d3, electron_v_nucleons, nabund, 1d0, guess_1, guess_2,&
  !         reltolerance)
  !   ! Tx = 0.525*tab_T(1)
  !   nChi = nxIso_mine(Tx, Nwimps)
  !   !**********************************************************************************************************
  !   !Mean Free Path
  !   rchi = sqrt(3*kb*tab_T(1)/(2*pi*GN*tab_starrho(1)*(mdmKg*1d3)))
  !
  !   do i = 1,nlines
  !     zeta_q(i) = q0_cgs/(mdm_g*sqrt(2.d0*kb*tab_T(i)/mdm_g))
  !     zeta_v(i) = v0/(sqrt(2.*kb*tab_T(i)/mdm_g))
  !   end do
  !
  !   if ((nq .eq. 0) .and. (nv .eq. 0)) then
  !     do i = 1,nlines
  !       mfp(i) = 1./(2.d0*sigma_0*nabund(i))  ! Since sigma_tot = 2*sigma_0 for v/q independent scattering
  !     end do
  !   else if ((nq .eq. 1)) then
  !     do i = 1,nlines
  !       mfp(i) = 1./(6.*nabund(i)*sigma_0/(1.+mdm_g/targetmass_g)/(zeta_q(i)**2))
  !     end do
  !   else if ((nq .eq. 2)) then
  !     do i = 1,nlines
  !       mfp(i) = 1./(40.*nabund(i)*sigma_0/((1.+mdm_g/targetmass_g)**2)/(zeta_q(i)**4))
  !     end do
  !   else if ((nq .eq. -1)) then
  !     do i = 1,nlines
  !       mfp(i) = 1./(nabund(i)*sigma_0*(1.+mdm_g/targetmass_g)*zeta_q(i)**2)
  !       print*, nq
  !       print*, (nq .eq. -1)
  !     end do
  !   else if ((nv .eq. 1)) then
  !     do i = 1,nlines
  !       mfp(i) = 1./(nabund(i)*2*sigma_0*(1.+mdm_g/targetmass_g)*3./2./(zeta_v(i)**2))
  !     end do
  !   else if ((nv .eq. 2)) then
  !     do i = 1,nlines
  !       mfp(i) = 1./(nabund(i)*2*sigma_0*((1.+mdm_g/targetmass_g)**2)*15./4./(zeta_v(i)**4))
  !     end do
  !   else if ((nv .eq. -1)) then
  !     do i = 1,nlines
  !       mfp(i) = 1./(nabund(i)*2*sigma_0*2*zeta_v(i)**2/(1.+mdm_g/targetmass_g))
  !     end do
  !   end if
  !
  !   !**********************************************************************************************************
  !   !For suppression factor
  !
  !   do i = 1, nlines
  !     tau_r(i) = trapz(Rsun*tab_r(i:nlines), 1/(mfp(i:nlines)), (nlines-i))
  !   end do
  !
  !   eta_ang = 7./10.*(1-exp(-10.*tau_r/7.))/tau_r
  !   phi_r = -tab_vesc**2/2.d0
  !   eta_multi = exp(1.5*tau_r/phi_r) !FIXXX!!!!!!!!!!!!!!!!!!!!!!!!
  !   supp_r = exp(-tau_r)*eta_multi*eta_ang
  !   supp_r = 1d0
  !   !**********************************************************************************************************
  !
  !   !Loop over the shells of constant radius in the star
  !   do ri = 1, nlines
  !     ue_at_r = sqrt(2*kb_here*tab_T(ri)/mtargetKg)*100 !in cm/s like the other velocities
  !     uChi = sqrt(2*kb_here*Tx/mdmKg)*100
  !     vesc = tab_vesc(ri)
  !     vesc_shared = vesc !make accessible via the module
  !
  !     vcritShared = vesc*velFactor
  !     mu = mx_in/(mtarget)
  !     muplus = (1.+mu)/2.
  !     muminus = (mu-1.d0)/2.
  !
  !     q0 = 0.04
  !
  !     ! Bottom part of the integral is always zero -- happy little slow DM particles can always be captured.
  !     wmin = 0.d0
  !     ! Chop the top of the integral off at the smaller of the halo escape velocity or the minimum velocity required for capture.
  !     ! wmax = min(vesc*sqrt(mu)/abs(muminus), vesc_halo)
  !     wmax = 3d0*vesc
  !
  !     !Call integrator
  !     call dsntdqagse(integrand_RPlus,fchi,wmin,wmax, &
  !     epsabs,epsrel,limit,int_result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
  !
  !     ! print*, "*************** Debugging intresult ************"
  !     ! print*, "mu: ", mu
  !     ! print*, "wmax: ", wmax
  !     ! print*, "vesc: ", vesc
  !     ! print*, "ue: ", ue_at_r
  !     ! print*, "uChi: ", uChi
  !     ! print*, "int_result: ", int_result
  !
  !     if (isnan(int_result)) then
  !       print*, "ri", ri
  !       int_result = 0d0
  !     end if
  !     evapRate = evapRate + supp_r(ri)*tab_r(ri)**2*int_result*tab_dr(ri)*nChi(ri)&
  !               *sigma_0*tab_starrho(ri)*mfr(ri)/(mtargetKg*1d3)
  !   end do
  !
  !   EvapRate = (4.d0*pi)**2.*Rsun**3*evapRate
  !
  ! end subroutine evap_RPlus

    !! captn_specific calculates the capture rate for constant cross section.
    ! subroutine captn_specific(mx_in,sigma_0,capped_SD,capped_SI)
    !   implicit none
    !   double precision, intent(in) :: mx_in, sigma_0
    !   double precision :: capped_SD,capped_SI
    !   call captn_general(mx_in,sigma_0,1,0,0,1,capped_SD)
    !   call captn_general(mx_in,sigma_0,29,0,0,0,capped_SI)
    ! end subroutine captn_specific

    subroutine captn_specific(mx_in,sigma_0_SD_in,sigma_0_SI_in,capped_SD,capped_SI)
      implicit none
      double precision, intent(in) :: mx_in, sigma_0_SD_in,sigma_0_SI_in
      double precision :: capped_SD,capped_SI

      call captn_general(mx_in,sigma_0_SD_in,1,0,0,1,capped_SD)
      call captn_general(mx_in,sigma_0_SI_in,29,0,0,0,capped_SI)
    end subroutine captn_specific


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! For mesa interface only: allocate arrays.
  subroutine allocate_stellar_arrays(nlines_mesa)
    use capmod
    integer, intent(in) :: nlines_mesa
    nlines = nlines_mesa
    allocate(tab_mencl(nlines))       !M(<r)
    allocate(tab_r(nlines))           !r
    allocate(tab_starrho(nlines))     !rho
    allocate(tab_mfr(nlines,8))       !mass fraction per isotope
    allocate(tab_atomic(8))
    allocate(tab_vesc(nlines))        !local escape velocity
    allocate(tab_T(nlines))           !temperature
    ! allocate(phi(nlines)) !! <--- not needed; computed in wimp_support.f
    allocate(tab_dr(nlines))          !dr (nice)
    allocate(tab_g(nlines))           !local gravitational acceleration, needed for transport
    allocate(tab_electron_mfr(nlines))
    RETURN
  end subroutine allocate_stellar_arrays

  subroutine deallocate_stellar_arrays()
    use capmod
    deallocate(tab_mencl)
    deallocate(tab_r)
    deallocate(tab_starrho)
    deallocate(tab_mfr) !we could just allocate niso, but this leads to problems
    deallocate(tab_atomic)
    deallocate(tab_vesc)
    deallocate(tab_T)
    deallocate(tab_dr)
    deallocate(tab_g)
    deallocate(tab_electron_mfr)
    RETURN
  end subroutine deallocate_stellar_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This is called INSTEAD of get_solar_params, for use with MESA interface.
  subroutine get_stellar_params(rmesa,rhomesa,mfrmesa,atomicmesa,mesavesc,Tmesa, &
                                mesag,mesamass,mesaradius,rho0_in,usun_in,u0_in,vesc_in)
    use capmod
    !mesamass & mesaradius unused here but subroutine used in a few other places so I left them
    !in just in case
    double precision :: mesamass, mesaradius
    double precision :: rhomesa(nlines), rmesa(nlines), mfrmesa(8,nlines)
    double precision :: mesavesc(nlines),mesag(nlines),Tmesa(nlines)
    double precision :: atomicmesa(8)
    integer i
    double precision,intent(in) :: rho0_in,usun_in,u0_in,vesc_in

    usun = usun_in*1.d5
    u0 =  u0_in*1.d5
    rho0 =rho0_in
    vesc_halo = vesc_in*1.d5

    Rsun = rmesa(nlines)
    tab_r = rmesa/Rsun
    tab_starrho = rhomesa
    tab_vesc = mesavesc
    tab_T = tmesa
    tab_g = -mesag
    do i= 1,8
      tab_mfr(:,i) = mfrmesa(i,:)
    end do
    tab_electron_mfr = tab_mfr(:,1)*melectron/mnuc
    tab_atomic = atomicmesa
    AtomicNumber(1:8) = tab_atomic

    do i = 1, nlines-1
      tab_dr(i) = -tab_r(i)+tab_r(i+1) !while we're here, populate dr
    end do
    tab_dr(nlines) = tab_r(nlines)-tab_r(nlines-1)

    RETURN
  end subroutine get_stellar_params


  subroutine getnlines(nlines_out) !a little auxiliary trick
    use capmod
    integer, intent(out) :: nlines_out
    nlines_out = nlines
    return
  end subroutine getnlines
