      function response_dm_m(m_proton, c0, tau1, tau2, term, j_chi,
     &   couplings_nreo)
      implicit none
      double precision :: response_dm_m, m_proton, c0, j_chi
      double precision :: couplings_nreo(14,2)
      integer :: tau1,tau2
      integer :: term
      integer :: c, v2, q2, v2q2
      double precision :: c1, c1p
      double precision :: c5, c5p
      double precision :: c8, c8p
      double precision :: c11, c11p
c
      c = 0
      v2 = 1
      q2 = 2
      v2q2 = 3
c
      c1 = couplings_nreo(1,tau1)
      c1p = couplings_nreo(1,tau2)
      c5 = couplings_nreo(4,tau1)
      c5p = couplings_nreo(4,tau2)
      c8 = couplings_nreo(7,tau1)
      c8p = couplings_nreo(7,tau2)
      c11 = couplings_nreo(10,tau1)
      c11p = couplings_nreo(10,tau2)
c
      if (term.eq.c) then
         response_dm_m = c1*c1p
      else if (term.eq.v2) then
         response_dm_m = 1/c0**2 * (j_chi*(j_chi+1))/3. * c8*c8p
      else if (term.eq.q2) then
         response_dm_m = 1/m_proton**2 * (j_chi*(j_chi+1))/3. * c11*c11p
      else if (term.eq.v2q2) then
         response_dm_m = 1/c0**2 * 1/m_proton**2 *
     &      (j_chi*(j_chi+1))/3. * c5*c5p
      else
         response_dm_m = 0.
      end if
c
      end function response_dm_m
