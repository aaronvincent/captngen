      function RP2(m_N, tau, taup, term, j_chi, coupling_Array)
      implicit none
      double precision :: RP2, m_N, j_chi
      double precision :: coupling_Array(14,2)
      integer :: tau,taup
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
      c3 = coupling_Array(2,tau)
      c3p = coupling_Array(2,taup)
      c12 = coupling_Array(11,tau)
      c12p = coupling_Array(11,taup)
      c15 = coupling_Array(14,tau)
      c15p = coupling_Array(14,taup)
c
      if (term.eq.c) then
         RP2 = (j_chi*(j_chi+1))/12. * c12*c12p
      else if (term.eq.q2) then
         RP2 = 1/m_N**2 * (1/4. * c3*c3p +
     &      (j_chi*(j_chi+1))/12. * (-c12*c15p-c15*c12p)) 
      else if (term.eq.q4) then
         RP2 = 1/m_N**4 * (j_chi*(j_chi+1))/12. * c15*c15p
      else
         RP2 = 0.
      end if
c
      end function RP2
