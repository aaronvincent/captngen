      function RS2(m_N, c0, tau, taup, term, j_chi, coupling_Array)
      implicit none
      double precision :: RS2, m_N, c0, j_chi
      double precision :: coupling_Array(14,2)
      integer :: tau,taup
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
      c4 = coupling_Array(3,tau)
      c4p = coupling_Array(3,taup)
      c6 = coupling_Array(5,tau)
      c6p = coupling_Array(5,taup)
      c10 = coupling_Array(9,tau)
      c10p = coupling_Array(9,taup)
      c12 = coupling_Array(11,tau)
      c12p = coupling_Array(11,taup)
      c13 = coupling_Array(12,tau)
      c13p = coupling_Array(12,taup)
c
      if (term.eq.c) then
         RS2 = (j_chi*(j_chi+1))/12. * c4*c4p
      else if (term.eq.v2) then
         RS2 = 1/c0**2 * (j_chi*(j_chi+1))/12. * c12*c12p
      else if (term.eq.q2) then
         RS2 = 1/m_N**2 * 1/4. * c10*c10p
      else if (term.eq.v2q2) then
         RS2 = 1/c0**2 * 1/m_N**2 * (j_chi*(j_chi+1))/12. * c13*c13p
      else if (term.eq.q4) then
         RS2 = 1/m_N**4 * (j_chi*(j_chi+1))/12. * c6*c6p
      else
         RS2 = 0.
      end if
c
      end function RS2
