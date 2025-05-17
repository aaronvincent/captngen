      function response_dm_s2(m_proton, c0, tau1, tau2, term, j_chi,
     &   couplings_nreo)
      implicit none
      double precision :: response_dm_s2, m_proton, c0, j_chi
      double precision :: couplings_nreo(14,2)
      integer :: tau1,tau2
      integer :: term
      integer :: c, v2, q2, v2q2, q4
      double precision :: c4, c4p
      double precision :: c6, c6p
      double precision :: c10, c10p
      double precision :: c12, c12p
      double precision :: c13, c13p
c
      c = 0
      v2 = 1
      q2 = 2
      v2q2 = 3
      q4 = 4
c
      c4 = couplings_nreo(3,tau1)
      c4p = couplings_nreo(3,tau2)
      c6 = couplings_nreo(5,tau1)
      c6p = couplings_nreo(5,tau2)
      c10 = couplings_nreo(9,tau1)
      c10p = couplings_nreo(9,tau2)
      c12 = couplings_nreo(11,tau1)
      c12p = couplings_nreo(11,tau2)
      c13 = couplings_nreo(12,tau1)
      c13p = couplings_nreo(12,tau2)
c
      if (term.eq.c) then
         response_dm_s2 = (j_chi*(j_chi+1))/12. * c4*c4p
      else if (term.eq.v2) then
         response_dm_s2 = 1/c0**2 * (j_chi*(j_chi+1))/12. * c12*c12p
      else if (term.eq.q2) then
         response_dm_s2 = 1/m_proton**2 * (1/4. * c10*c10p +
     &      (j_chi*(j_chi+1))/12. * (c4*c6p+c6*c4p))
      else if (term.eq.v2q2) then
         response_dm_s2 = 1/c0**2 * 1/m_proton**2 *
     &      (j_chi*(j_chi+1))/12. * c13*c13p
      else if (term.eq.q4) then
         response_dm_s2 = 1/m_proton**4 * (j_chi*(j_chi+1))/12. * c6*c6p
      else
         response_dm_s2 = 0.
      end if
c
      end function response_dm_s2
