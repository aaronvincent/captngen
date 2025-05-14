      function RS1(m_N, c0, tau, taup, term, j_chi, couplings_nreo)
      implicit none
      double precision :: RS1, m_N, c0, j_chi
      double precision :: couplings_nreo(14,2)
      integer :: tau,taup
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
      c3 = couplings_nreo(2,tau)
      c3p = couplings_nreo(2,taup)
      c4 = couplings_nreo(3,tau)
      c4p = couplings_nreo(3,taup)
      c7 = couplings_nreo(6,tau)
      c7p = couplings_nreo(6,taup)
      c9 = couplings_nreo(8,tau)
      c9p = couplings_nreo(8,taup)
      c12 = couplings_nreo(11,tau)
      c12p = couplings_nreo(11,taup)
      c14 = couplings_nreo(13,tau)
      c14p = couplings_nreo(13,taup)
      c15 = couplings_nreo(14,tau)
      c15p = couplings_nreo(14,taup)
c
      if (term.eq.c) then
         RS1 = (j_chi*(j_chi+1))/12. * c4*c4p
      else if (term.eq.v2) then
         RS1 = 1/c0**2 * (1/8. * c7*c7p +
     &      (j_chi*(j_chi+1))/24. * c12*c12p)
      else if (term.eq.q2) then
         RS1 = 1/m_N**2 * (j_chi*(j_chi+1))/12. * c9*c9p
      else if (term.eq.v2q2) then
         RS1 = 1/c0**2 * 1/m_N**2 * (1/8. * c3*c3p +
     &      (j_chi*(j_chi+1))/24. * (c14*c14p - c12*c15p - c15*c12p))
      else if (term.eq.v2q4) then
         RS1 = 1/c0**2 * 1/m_N**4 * (j_chi*(j_chi+1))/24. * c15*c15p
      else
         RS1 = 0.
      end if
c
      end function RS1
