      function WS1D(i,j,target,term)
      implicit none
      real y,WS1D
      integer i,j
      character (len=4) :: target
      character (len=2) :: term
      !include 'dsddcom.h'
c
      if (target.eq."Ni58") then
         WS1D = 0.
      end if
c     
      if (target.eq."Fe56") then
         WS1D = 0.
      end if
c
      if (target.eq."Ca40") then
         WS1D = 0.
      end if
c
      if (target.eq."Ar40") then
         WS1D = 0.
      end if
c
      if (target.eq."S32") then
         WS1D = 0.
      end if
c
      if (target.eq."Si28") then
         WS1D = 0.
      end if
c
      if (target.eq."Al27") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WS1D = -0.08832432609168217
            else if (term.eq."y1") then
               WS1D = 0.18577474668344507
            else if (term.eq."y2") then
               WS1D = -0.10400051271821338
            else if (term.eq."y3") then
               WS1D = 0.016363536090601292
            else
               WS1D = 0.
            end if
         else if ((i.eq.1).and.(j.eq.1)) then
            if (term.eq."y0") then
               WS1D = -0.05010448456928069
            else if (term.eq."y1") then
               WS1D = 0.11484454255964766
            else if (term.eq."y2") then
               WS1D = -0.07298984507998893
            else if (term.eq."y3") then
               WS1D = 0.013131504730800993
            else
               WS1D = 0.
            end if
         else if ((i.eq.1).and.(j.eq.0)) then
            if (term.eq."y0") then
               WS1D = -0.07427314267833456
            else if (term.eq."y1") then
               WS1D = 0.17024154960752716
            else if (term.eq."y2") then
               WS1D = -0.1057443365406475
            else if (term.eq."y3") then
               WS1D = 0.018819732308840584
            else
               WS1D = 0.
            end if
         else
            if (term.eq."y0") then
               WS1D = -0.05958337932351564
            else if (term.eq."y1") then
               WS1D = 0.12532320018926416
            else if (term.eq."y2") then
               WS1D = -0.07172038098440625
            else if (term.eq."y3") then
               WS1D = 0.011397990801132206
            else
               WS1D = 0.
            end if
         end if
      end if
c
      if (target.eq."Mg24") then
         WS1D = 0.
      end if
c
      if (target.eq."Na23") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WS1D = -0.0291634
            else if (term.eq."y1") then
               WS1D = 0.0548817
            else if (term.eq."y2") then
               WS1D = -0.0305345
            else if (term.eq."y3") then
               WS1D = 0.00476387
            else
               WS1D = 0.
            end if         
         else if ((i.eq.1).and.(j.eq.1)) then
            if (term.eq."y0") then
               WS1D = -0.0119052
            else if (term.eq."y1") then
               WS1D = 0.0231539
            else if (term.eq."y2") then
               WS1D = -0.0164035
            else if (term.eq."y3") then
               WS1D = 0.00310235
            else
               WS1D = 0.
            end if
         else if ((i.eq.1).and.(j.eq.0)) then
            if (term.eq."y0") then
               WS1D = -0.024821
            else if (term.eq."y1") then
               WS1D = 0.0482732
            else if (term.eq."y2") then
               WS1D = -0.02884
            else if (term.eq."y3") then
               WS1D = 0.00481368
            else
               WS1D = 0.
            end if
         else 
            if (term.eq."y0") then
               WS1D = -0.013988
            else if (term.eq."y1") then
               WS1D = 0.0263236
            else if (term.eq."y2") then
               WS1D = -0.0171362
            else if (term.eq."y3") then
               WS1D = 0.00306717
            else
               WS1D = 0.
            end if
         end if
      end if
c
      if (target.eq."Ne20") then
         WS1D = 0.
      end if
c
      if (target.eq."O16") then
         WS1D = 0.
      end if
c
      if (target.eq."N14") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WS1D = 0.053445072530347415
            else if (term.eq."y1") then
               WS1D = -0.07553245404282635
            else
               WS1D = 0.
            end if
         else
            WS1D = 0.
         end if
      end if
c
      if (target.eq."C12") then
         WS1D = 0.
      end if
c
      if (target.eq."He4") then
         WS1D = 0.
      end if
c     
      if (target.eq."He3") then
         WS1D = 0.
      end if
c
      if (target.eq."H") then
         WS1D = 0.
      end if
c
      end function WS1D
