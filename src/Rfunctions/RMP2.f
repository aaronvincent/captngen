      function RMP2(m_N, tau, taup, term, j_chi, coupling_Array)
      implicit none
      double precision :: RMP2, m_N, j_chi
      double precision :: coupling_Array(14,2)
      integer :: tau,taup
      integer :: term
      integer :: c, q2
      double precision :: c3, c1p
      double precision :: c12, c15, c11p
c
      c = 0
      q2 = 2
c
      c3 = coupling_Array(2,tau)
      c1p = coupling_Array(1,taup)
      c12 = coupling_Array(11,tau)
      c15 = coupling_Array(14,tau)
      c11p = coupling_Array(10,taup)
c
      if (term.eq.c) then
         RMP2 = c3*c1p + (j_chi*(j_chi+1))/3. * c12*c11p
      else if (term.eq.q2) then
         RMP2 = 1/m_N**2 * (j_chi*(j_chi+1))/3. * (-c15*c11p)
      else
         RMP2 = 0.
      end if
c
      end function RMP2
