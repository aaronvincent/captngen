      function RM(m_N, c0, tau, taup, term, j_chi, couplings_nreo)
      implicit none
      double precision :: RM, m_N, c0, j_chi
      double precision :: couplings_nreo(14,2)
      integer :: tau,taup
      integer :: term
      integer :: c, v2, q2, v2q2
      double precision :: c1, c1p
      double precision :: c5, c5p
      double precision :: c8, c8p
      double precision :: c11, c11p
c
      c = 0
      v2 = 1
      q2 = 2
      v2q2 = 3
c
      c1 = couplings_nreo(1,tau)
      c1p = couplings_nreo(1,taup)
      c5 = couplings_nreo(4,tau)
      c5p = couplings_nreo(4,taup)
      c8 = couplings_nreo(7,tau)
      c8p = couplings_nreo(7,taup)
      c11 = couplings_nreo(10,tau)
      c11p = couplings_nreo(10,taup)
c
      if (term.eq.c) then
         RM = c1*c1p
      else if (term.eq.v2) then
         RM = 1/c0**2 * (j_chi*(j_chi+1))/3. * c8*c8p
      else if (term.eq.q2) then
         RM = 1/m_N**2 * (j_chi*(j_chi+1))/3. * c11*c11p
      else if (term.eq.v2q2) then
         RM = 1/c0**2 * 1/m_N**2 * (j_chi*(j_chi+1))/3. * c5*c5p
      else
         RM = 0.
      end if
c
      end function RM
