      function WD(i,j,target,term)
      implicit none
      real y,WD
      integer i,j
      character (len=4) :: target
      character (len=2) :: term
      !include 'dsddcom.h'
c
      if (target.eq."Ni58") then
         WD = 0.
      end if
c
      if (target.eq."Fe56") then
         WD = 0.
      end if
c
      if (target.eq."Ca40") then
         WD = 0.
      end if
c
      if (target.eq."Ar40") then
         WD = 0.
      end if
c
      if (target.eq."S32") then
         WD = 0.
      end if
c
      if (target.eq."Si28") then
         WD = 0.
      end if
c
      if (target.eq."Al27") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WD = 0.12604325181306578
            else if (term.eq."y1") then
               WD = -0.10083460145045263
            else if (term.eq."y2") then
               WD = 0.0237576551261129
            else
               WD = 0.
            end if         
         else if ((i.eq.1).and.(j.eq.1)) then
            if (term.eq."y0") then
               WD = 0.05736000705306767
            else if (term.eq."y1") then
               WD = -0.04588800564245413
            else if (term.eq."y2") then
               WD = 0.012102012478740131
            else
               WD = 0.
            end if
         else
            if (term.eq."y0") then
               WD = 0.08502847648281743
            else if (term.eq."y1") then
               WD = -0.06802278118625395
            else if (term.eq."y2") then
               WD = 0.01684504782190263
            else
               WD = 0.
            end if
         end if
      end if
c
      if (target.eq."Mg24") then
         WD = 0.
      end if
c
      if (target.eq."Na23") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WD = 0.0335711
            else if (term.eq."y1") then
               WD = -0.0268568
            else if (term.eq."y2") then
               WD = 0.00656896
            else
               WD = 0.
            end if         
         else if ((i.eq.1).and.(j.eq.1)) then
            if (term.eq."y0") then
               WD = 0.00772326
            else if (term.eq."y1") then
               WD = -0.00617861
            else if (term.eq."y2") then
               WD = 0.0021619
            else
               WD = 0.
            end if
         else 
            if (term.eq."y0") then
               WD = 0.0161021
            else if (term.eq."y1") then
               WD = -0.0128817
            else if (term.eq."y2") then
               WD = 0.00362952
            else
               WD = 0.
            end if
         end if
      end if
c
      if (target.eq."Ne20") then
         WD = 0.
      end if
c
      if (target.eq."O16") then
         WD = 0.
      end if
c
      if (target.eq."N14") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WD = 0.042407469580591495
            else
               WD = 0.
            end if
         else
            WD = 0.
         end if
      end if
c
      if (target.eq."C12") then
         WD = 0.
      end if
c
      if (target.eq."He4") then
         WD = 0.
      end if
c     
      if (target.eq."He3") then
         WD = 0.
      end if
c
      if (target.eq."H") then
         WD = 0.
      end if
c
      end function WD
