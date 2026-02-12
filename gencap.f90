!   Capt'n General
!   For a continuum of Q-dependent capture
!   Simplified, general solar DM capture routine
!   Standalone code for q^2n, v^2n
!   Useful stuff is run at the end; beginning is the module that does the heavy lifting
!   Future plans: add form factor handling (a la Catena & Schwabe)
!   Made for GAMBIT, with marginal competence
!   Aaron Vincent 2017
!   all units of distance: cm
!   all units of mass/energy : GeV (or GeV/c^2, don't forget)
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
      double precision, parameter :: GNewt = 6.672d-8
      double precision, parameter :: q0 = 0.04,v0 = 220.d5

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

      !Fast trapezoidal integral
      function trapz(x,y,flen)
        implicit none
        integer, intent(in) :: flen
        double precision, intent (in) :: x(flen), y(flen)
        double precision trapz
        integer i

        trapz = y(1)*(x(2)-x(1))/2.

        do i = 2,flen-1
          trapz = trapz + y(i)*(x(i)-x(i-1))
        end do

        trapz = trapz + y(flen)*(x(flen)-x(flen-1))/2.

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
