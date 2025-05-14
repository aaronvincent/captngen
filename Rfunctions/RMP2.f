      function RMP2(m_N, tau, taup, term, j_chi, couplings_nreo)
      implicit none
      double precision :: RMP2, m_N, j_chi
      double precision :: couplings_nreo(14,2)
      integer :: tau,taup
      integer :: term
      integer :: c, q2
      double precision :: c3, c1p
      double precision :: c12, c15, c11p
c
      c = 0
      q2 = 2
c
      c3 = couplings_nreo(2,tau)
      c1p = couplings_nreo(1,taup)
      c12 = couplings_nreo(11,tau)
      c15 = couplings_nreo(14,tau)
      c11p = couplings_nreo(10,taup)
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
