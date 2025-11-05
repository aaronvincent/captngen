      function RP1(m_N, tau, taup, term, j_chi, coupling_Array)
      implicit none
      double precision :: RP1, m_N, j_chi
      double precision :: coupling_Array(14,2)
      integer :: tau,taup
      integer :: term
      integer :: c, q2
      double precision :: c12, c12p
      double precision :: c13, c13p
c
      c = 0
      q2 = 2
c
      c12 = coupling_Array(11,tau)
      c12p = coupling_Array(11,taup)
      c13 = coupling_Array(12,tau)
      c13p = coupling_Array(12,taup)
c
      if (term.eq.c) then
         RP1 = (j_chi*(j_chi+1))/12. * c12*c12p
      else if (term.eq.q2) then
         RP1 = 1/m_N**2 * (j_chi*(j_chi+1))/12. * c13*c13p
      else
         RP1 = 0.
      end if
c
      end function RP1
