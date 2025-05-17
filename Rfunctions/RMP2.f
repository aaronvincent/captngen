      function response_dm_p2m(m_proton, tau1, tau2, term, j_chi,
     &   couplings_nreo)
      implicit none
      double precision :: response_dm_p2m, m_proton, j_chi
      double precision :: couplings_nreo(14,2)
      integer :: tau1,tau2
      integer :: term
      integer :: c, q2
      double precision :: c3, c1p
      double precision :: c12, c15, c11p
c
      c = 0
      q2 = 2
c
      c3 = couplings_nreo(2,tau1)
      c1p = couplings_nreo(1,tau2)
      c12 = couplings_nreo(11,tau1)
      c15 = couplings_nreo(14,tau1)
      c11p = couplings_nreo(10,tau2)
c
      if (term.eq.c) then
         response_dm_p2m = c3*c1p + (j_chi*(j_chi+1))/3. * c12*c11p
      else if (term.eq.q2) then
         response_dm_p2m = 1/m_proton**2 * (j_chi*(j_chi+1))/3. *
     &      (-c15*c11p)
      else
         response_dm_p2m = 0.
      end if
c
      end function response_dm_p2m
