      function response_dm_d(m_proton, tau1, tau2, term, j_chi,
     &   couplings_nreo)
      implicit none
      double precision :: response_dm_d, m_proton, j_chi
      double precision :: couplings_nreo(14,2)
      integer :: tau1,tau2
      integer :: term
      integer :: c, q2
      double precision :: c5, c5p
      double precision :: c8, c8p
c
      c = 0
      q2 = 2
c
      c5 = couplings_nreo(4,tau1)
      c5p = couplings_nreo(4,tau2)
      c8 = couplings_nreo(7,tau1)
      c8p = couplings_nreo(7,tau2)
c
      if (term.eq.c) then
         response_dm_d = (j_chi*(j_chi+1))/3. * c8*c8p
      else if (term.eq.q2) then
         response_dm_d = 1/m_proton**2 * (j_chi*(j_chi+1))/3. * c5*c5p
      else
         response_dm_d = 0.
      end if
c
      end function response_dm_d
