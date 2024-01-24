      function RS1D(tau, taup, term, j_chi, coupling_Array)
      implicit none
      double precision :: RS1D, j_chi
      double precision :: coupling_Array(14,2)
      integer :: tau,taup
      integer :: term
      integer :: c
      double precision :: c5, c4p
      double precision :: c8, c9p
c
      c = 0
c
      c5 = coupling_Array(4,tau)
      c4p = coupling_Array(3,taup)
      c8 = coupling_Array(7,tau)
      c9p = coupling_Array(8,taup)
c
      if (term.eq.c) then
         RS1D = (j_chi*(j_chi+1))/3. * (c5*c4p - c8*c9p)
      else
         RS1D = 0.
      end if
c
      end function RS1D
