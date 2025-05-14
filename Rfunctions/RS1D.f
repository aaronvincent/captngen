      function RS1D(tau, taup, term, j_chi, couplings_nreo)
      implicit none
      double precision :: RS1D, j_chi
      double precision :: couplings_nreo(14,2)
      integer :: tau,taup
      integer :: term
      integer :: c
      double precision :: c5, c4p
      double precision :: c8, c9p
c
      c = 0
c
      c5 = couplings_nreo(4,tau)
      c4p = couplings_nreo(3,taup)
      c8 = couplings_nreo(7,tau)
      c9p = couplings_nreo(8,taup)
c
      if (term.eq.c) then
         RS1D = (j_chi*(j_chi+1))/3. * (c5*c4p - c8*c9p)
      else
         RS1D = 0.
      end if
c
      end function RS1D
