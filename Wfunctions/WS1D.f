      function response_nuc_ds1(i,j,target,term)
      implicit none
      real y,response_nuc_ds1
      integer i,j
      character (len=4) :: target
      character (len=2) :: term
      !include 'dsddcom.h'
c
      if (target.eq."Ni58") then
         response_nuc_ds1 = 0.
      end if
c     
      if (target.eq."Fe56") then
         response_nuc_ds1 = 0.
      end if
c
      if (target.eq."Ca40") then
         response_nuc_ds1 = 0.
      end if
c
      if (target.eq."Ar40") then
         response_nuc_ds1 = 0.
      end if
c
      if (target.eq."S32") then
         response_nuc_ds1 = 0.
      end if
c
      if (target.eq."Si28") then
         response_nuc_ds1 = 0.
      end if
c
      if (target.eq."Al27") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               response_nuc_ds1 = -0.08832432609168217
            else if (term.eq."y1") then
               response_nuc_ds1 = 0.18577474668344507
            else if (term.eq."y2") then
               response_nuc_ds1 = -0.10400051271821338
            else if (term.eq."y3") then
               response_nuc_ds1 = 0.016363536090601292
            else
               response_nuc_ds1 = 0.
            end if
         else if ((i.eq.1).and.(j.eq.1)) then
            if (term.eq."y0") then
               response_nuc_ds1 = -0.05010448456928069
            else if (term.eq."y1") then
               response_nuc_ds1 = 0.11484454255964766
            else if (term.eq."y2") then
               response_nuc_ds1 = -0.07298984507998893
            else if (term.eq."y3") then
               response_nuc_ds1 = 0.013131504730800993
            else
               response_nuc_ds1 = 0.
            end if
         else if ((i.eq.1).and.(j.eq.0)) then
            if (term.eq."y0") then
               response_nuc_ds1 = -0.07427314267833456
            else if (term.eq."y1") then
               response_nuc_ds1 = 0.17024154960752716
            else if (term.eq."y2") then
               response_nuc_ds1 = -0.1057443365406475
            else if (term.eq."y3") then
               response_nuc_ds1 = 0.018819732308840584
            else
               response_nuc_ds1 = 0.
            end if
         else
            if (term.eq."y0") then
               response_nuc_ds1 = -0.05958337932351564
            else if (term.eq."y1") then
               response_nuc_ds1 = 0.12532320018926416
            else if (term.eq."y2") then
               response_nuc_ds1 = -0.07172038098440625
            else if (term.eq."y3") then
               response_nuc_ds1 = 0.011397990801132206
            else
               response_nuc_ds1 = 0.
            end if
         end if
      end if
c
      if (target.eq."Mg24") then
         response_nuc_ds1 = 0.
      end if
c
      if (target.eq."Na23") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               response_nuc_ds1 = -0.0291634
            else if (term.eq."y1") then
               response_nuc_ds1 = 0.0548817
            else if (term.eq."y2") then
               response_nuc_ds1 = -0.0305345
            else if (term.eq."y3") then
               response_nuc_ds1 = 0.00476387
            else
               response_nuc_ds1 = 0.
            end if         
         else if ((i.eq.1).and.(j.eq.1)) then
            if (term.eq."y0") then
               response_nuc_ds1 = -0.0119052
            else if (term.eq."y1") then
               response_nuc_ds1 = 0.0231539
            else if (term.eq."y2") then
               response_nuc_ds1 = -0.0164035
            else if (term.eq."y3") then
               response_nuc_ds1 = 0.00310235
            else
               response_nuc_ds1 = 0.
            end if
         else if ((i.eq.1).and.(j.eq.0)) then
            if (term.eq."y0") then
               response_nuc_ds1 = -0.024821
            else if (term.eq."y1") then
               response_nuc_ds1 = 0.0482732
            else if (term.eq."y2") then
               response_nuc_ds1 = -0.02884
            else if (term.eq."y3") then
               response_nuc_ds1 = 0.00481368
            else
               response_nuc_ds1 = 0.
            end if
         else 
            if (term.eq."y0") then
               response_nuc_ds1 = -0.013988
            else if (term.eq."y1") then
               response_nuc_ds1 = 0.0263236
            else if (term.eq."y2") then
               response_nuc_ds1 = -0.0171362
            else if (term.eq."y3") then
               response_nuc_ds1 = 0.00306717
            else
               response_nuc_ds1 = 0.
            end if
         end if
      end if
c
      if (target.eq."Ne20") then
         response_nuc_ds1 = 0.
      end if
c
      if (target.eq."O16") then
         response_nuc_ds1 = 0.
      end if
c
      if (target.eq."N14") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               response_nuc_ds1 = 0.053445072530347415
            else if (term.eq."y1") then
               response_nuc_ds1 = -0.07553245404282635
            else
               response_nuc_ds1 = 0.
            end if
         else
            response_nuc_ds1 = 0.
         end if
      end if
c
      if (target.eq."C12") then
         response_nuc_ds1 = 0.
      end if
c
      if (target.eq."He4") then
         response_nuc_ds1 = 0.
      end if
c     
      if (target.eq."He3") then
         response_nuc_ds1 = 0.
      end if
c
      if (target.eq."H") then
         response_nuc_ds1 = 0.
      end if
c
      end function response_nuc_ds1
