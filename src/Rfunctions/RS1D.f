      function response_dm_ds1(tau1, tau2, term, j_chi, couplings_nreo)
      implicit none
      double precision :: response_dm_ds1, j_chi
      double precision :: couplings_nreo(14,2)
      integer :: tau1,tau2
      integer :: term
      integer :: c
      double precision :: c5, c4p
      double precision :: c8, c9p
c
      c = 0
c
      c5 = couplings_nreo(4,tau1)
      c4p = couplings_nreo(3,tau2)
      c8 = couplings_nreo(7,tau1)
      c9p = couplings_nreo(8,tau2)
c
      if (term.eq.c) then
         response_dm_ds1 = (j_chi*(j_chi+1))/3. * (c5*c4p - c8*c9p)
      else
         response_dm_ds1 = 0.
      end if
c
      end function response_dm_ds1
