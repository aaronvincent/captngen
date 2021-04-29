      function WS2(i,j,target,term)
      implicit none
      real y,WS2
      integer i,j
      character (len=4) :: target
      character (len=2) :: term
      !include 'dsddcom.h'
c
      if (target.eq."Ni58") then
         WS2 = 0.
      end if
c
      if (target.eq."Fe56") then
         WS2 = 0.
      end if
c
      if (target.eq."Ca40") then
         WS2 = 0.
      end if
c
      if (target.eq."Ar40") then
         WS2 = 0.
      end if
c
      if (target.eq."S32") then
         WS2 = 0.
      end if
c
      if (target.eq."Si28") then
         WS2 = 0.
      end if
c
      if (target.eq."Al27") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WS2 = 0.030946466658602678
            else if (term.eq."y1") then
               WS2 = -0.036724195099056736
            else if (term.eq."y2") then
               WS2 = 0.02653470770381152
            else if (term.eq."y3") then
               WS2 = -0.002416057844169939
            else if (term.eq."y4") then
               WS2 = 0.011001067283047518
            else
               WS2 = 0.
            end if
         else if ((i.eq.1).and.(j.eq.1)) then
            if (term.eq."y0") then
               WS2 = 0.021883360052856037
            else if (term.eq."y1") then
               WS2 = -0.009444764771082234
            else if (term.eq."y2") then
               WS2 = 0.011506022741508951
            else if (term.eq."y3") then
               WS2 = 0.0009535368761278121
            else if (term.eq."y4") then
               WS2 = 0.010481297740411958
            else
               WS2 = 0.
            end if
         else
            if (term.eq."y0") then
               WS2 = 0.026023310170958405
            else if (term.eq."y1") then
               WS2 = -0.021056715592501225
            else if (term.eq."y2") then
               WS2 = 0.01586427118299763
            else if (term.eq."y3") then
               WS2 = 0.0006060767741112289
            else if (term.eq."y4") then
               WS2 = 0.010571282241782384
            else
               WS2 = 0.
            end if
         end if
      end if
c
      if (target.eq."Mg24") then
         WS2 = 0.
      end if
c
      if (target.eq."Na23") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WS2 = 0.0126672
            else if (term.eq."y1") then
               WS2 = -0.0262533
            else if (term.eq."y2") then
               WS2 = 0.0401886
            else if (term.eq."y3") then
               WS2 = -0.010514
            else if (term.eq."y4") then
               WS2 = 0.00078605
            else
               WS2 = 0.
            end if         
         else if ((i.eq.1).and.(j.eq.1)) then
            if (term.eq."y0") then
               WS2 = 0.00917577
            else if (term.eq."y1") then
               WS2 = -0.0167053
            else if (term.eq."y2") then
               WS2 = 0.0332751
            else if (term.eq."y3") then
               WS2 = -0.00765719
            else if (term.eq."y4") then
               WS2 = 0.000597676
            else
               WS2 = 0.
            end if
         else 
            if (term.eq."y0") then
               WS2 = 0.0107811
            else if (term.eq."y1") then
               WS2 = -0.020986
            else if (term.eq."y2") then
               WS2 = 0.0360971
            else if (term.eq."y3") then
               WS2 = -0.00876213
            else if (term.eq."y4") then
               WS2 = 0.000626718
            else
               WS2 = 0.
            end if
         end if
      end if
c
      if (target.eq."Ne20") then
         WS2 = 0.
      end if
c
      if (target.eq."O16") then
         WS2 = 0.
      end if
c
      if (target.eq."N14") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WS2 = 0.03367774365
            else if (term.eq."y1") then
               WS2 = 0.05567234823
            else if (term.eq."y2") then
               WS2 = 0.02300785342
            else
               WS2 = 0.
            end if
         else
            WS2 = 0.
         end if
      end if

c
      if (target.eq."C12") then
         WS2 = 0.
      end if

c
      if (target.eq."He4") then
       WS2 = 0.
      end if
c
      if (target.eq."He3") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WS2 = 0.0397887
            else
               WS2 = 0.
            end if
         else if ((i.eq.1).and.(j.eq.1)) then
            if (term.eq."y0") then
               WS2 = 0.0397887
            else
               WS2 = 0.
            end if
         else
            if (term.eq."y0") then
               WS2 = -0.0397887
            else
               WS2 = 0.
            end if
         end if
      end if
c
ccc...
ccc*H
ccc...
      if (target.eq."H") then
         if (term.eq."y0") then
            WS2 = 0.039788735772973836!/exp(2.*y)
         else
            WS2 = 0.
         end if
      end if
c
      end function WS2
