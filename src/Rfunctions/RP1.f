      function response_dm_p1(m_proton, tau1, tau2, term, j_chi,
     &   couplings_nreo)
      implicit none
      double precision :: response_dm_p1, m_proton, j_chi
      double precision :: couplings_nreo(14,2)
      integer :: tau1,tau2
      integer :: term
      integer :: c, q2
      double precision :: c12, c12p
      double precision :: c13, c13p
c
      c = 0
      q2 = 2
c
      c12 = couplings_nreo(11,tau1)
      c12p = couplings_nreo(11,tau2)
      c13 = couplings_nreo(12,tau1)
      c13p = couplings_nreo(12,tau2)
c
      if (term.eq.c) then
         response_dm_p1 = (j_chi*(j_chi+1))/12. * c12*c12p
      else if (term.eq.q2) then
         response_dm_p1 = 1/m_proton**2 * (j_chi*(j_chi+1))/12. *
     &      c13*c13p
      else
         response_dm_p1 = 0.
      end if
c
      end function response_dm_p1
