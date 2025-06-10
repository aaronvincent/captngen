      function response_dm_s1(m_proton, c0, tau1, tau2, term, j_chi,
     &   couplings_nreo)
      implicit none
      double precision :: response_dm_s1, m_proton, c0, j_chi
      double precision :: couplings_nreo(14,2)
      integer :: tau1,tau2
      integer :: term
      integer :: c, v2, q2, v2q2, v2q4
      double precision :: c3, c3p
      double precision :: c4, c4p
      double precision :: c7, c7p
      double precision :: c9, c9p
      double precision :: c12, c12p
      double precision :: c14, c14p
      double precision :: c15, c15p
c
      c = 0
      v2 = 1
      q2 = 2
      v2q2 = 3
      v2q4 = 5
c
      c3 = couplings_nreo(2,tau1)
      c3p = couplings_nreo(2,tau2)
      c4 = couplings_nreo(3,tau1)
      c4p = couplings_nreo(3,tau2)
      c7 = couplings_nreo(6,tau1)
      c7p = couplings_nreo(6,tau2)
      c9 = couplings_nreo(8,tau1)
      c9p = couplings_nreo(8,tau2)
      c12 = couplings_nreo(11,tau1)
      c12p = couplings_nreo(11,tau2)
      c14 = couplings_nreo(13,tau1)
      c14p = couplings_nreo(13,tau2)
      c15 = couplings_nreo(14,tau1)
      c15p = couplings_nreo(14,tau2)
c
      if (term.eq.c) then
         response_dm_s1 = (j_chi*(j_chi+1))/12. * c4*c4p
      else if (term.eq.v2) then
         response_dm_s1 = 1/c0**2 * (1/8. * c7*c7p +
     &      (j_chi*(j_chi+1))/24. * c12*c12p)
      else if (term.eq.q2) then
         response_dm_s1 = 1/m_proton**2 * (j_chi*(j_chi+1))/12. * c9*c9p
      else if (term.eq.v2q2) then
         response_dm_s1 = 1/c0**2 * 1/m_proton**2 * (1/8. * c3*c3p +
     &      (j_chi*(j_chi+1))/24. * (c14*c14p - c12*c15p - c15*c12p))
      else if (term.eq.v2q4) then
         response_dm_s1 = 1/c0**2 * 1/m_proton**4 *
     &   (j_chi*(j_chi+1))/24. * c15*c15p
      else
         response_dm_s1 = 0.
      end if
c
      end function response_dm_s1
