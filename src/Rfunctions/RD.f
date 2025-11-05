      function RD(m_N, tau, taup, term, j_chi, coupling_Array)
      implicit none
      double precision :: RD, m_N, j_chi
      double precision :: coupling_Array(14,2)
      integer :: tau,taup
      integer :: term
      integer :: c, q2
      double precision :: c5, c5p
      double precision :: c8, c8p
c
      c = 0
      q2 = 2
c
      c5 = coupling_Array(4,tau)
      c5p = coupling_Array(4,taup)
      c8 = coupling_Array(7,tau)
      c8p = coupling_Array(7,taup)
c
      if (term.eq.c) then
         RD = (j_chi*(j_chi+1))/3. * c8*c8p
      else if (term.eq.q2) then
         RD = 1/m_N**2 * (j_chi*(j_chi+1))/3. * c5*c5p
      else
         RD = 0.
      end if
c
      end function RD
