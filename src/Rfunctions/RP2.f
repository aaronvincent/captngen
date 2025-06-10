      function response_dm_p2(m_proton, tau1, tau2, term, j_chi,
     &   couplings_nreo)
      implicit none
      double precision :: response_dm_p2, m_proton, j_chi
      double precision :: couplings_nreo(14,2)
      integer :: tau1,tau2
      integer :: term
      integer :: c, q2, q4
      double precision :: c3, c3p
      double precision :: c12, c12p
      double precision :: c15, c15p
c
      c = 0
      q2 = 2
      q4 = 4
c
      c3 = couplings_nreo(2,tau1)
      c3p = couplings_nreo(2,tau2)
      c12 = couplings_nreo(11,tau1)
      c12p = couplings_nreo(11,tau2)
      c15 = couplings_nreo(14,tau1)
      c15p = couplings_nreo(14,tau2)
c
      if (term.eq.c) then
         response_dm_p2 = (j_chi*(j_chi+1))/12. * c12*c12p
      else if (term.eq.q2) then
         response_dm_p2 = 1/m_proton**2 * (1/4. * c3*c3p +
     &      (j_chi*(j_chi+1))/12. * (-c12*c15p-c15*c12p)) 
      else if (term.eq.q4) then
         response_dm_p2 = 1/m_proton**4 * (j_chi*(j_chi+1))/12. *
     &      c15*c15p
      else
         response_dm_p2 = 0.
      end if
c
      end function response_dm_p2
