      function WS1(i,j,target,term)
      implicit none
      real y,WS1
      integer i,j
      character (len=4) :: target
      character (len=2) :: term
      !include 'dsddcom.h'
c
      if (target.eq."Ni58") then
         WS1 = 0.
      end if
c
      if (target.eq."Fe56") then
         WS1 = 0.
      end if
c     
      if (target.eq."Ca40") then
         WS1 = 0.
      end if
c
      if (target.eq."Ar40") then
         WS1 = 0.
      end if
c
      if (target.eq."S32") then
         WS1 = 0.
      end if
c     
      if (target.eq."Si28") then
         WS1 = 0.
      end if
c
      if (target.eq."Al27") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WS1 = 0.0618929333172054
            else if (term.eq."y1") then
               WS1 = -0.2108475381697648
            else if (term.eq."y2") then
               WS1 = 0.24446620317090406
            else if (term.eq."y3") then
               WS1 = -0.09426815715239213
            else if (term.eq."y4") then
               WS1 = 0.02437372743452902
            else
               WS1 = 0.
            end if
         else if ((i.eq.1).and.(j.eq.1)) then
            if (term.eq."y0") then
               WS1 = 0.043766720105712115
            else if (term.eq."y1") then
               WS1 = -0.1656221156517662
            else if (term.eq."y2") then
               WS1 = 0.22119270508008307
            else if (term.eq."y3") then
               WS1 = -0.10199091700972805
            else if (term.eq."y4") then
               WS1 = 0.027747661383265603
            else
               WS1 = 0.
            end if
         else
            if (term.eq."y0") then
               WS1 = 0.05204662034191686
            else if (term.eq."y1") then
               WS1 = -0.18712976577516616
            else if (term.eq."y2") then
               WS1 = 0.23300730436550043
            else if (term.eq."y3") then
               WS1 = -0.0985081586005909
            else if (term.eq."y4") then
               WS1 = 0.025932726182971457
            else
               WS1 = 0.
            end if
         end if
      end if
c     
      if (target.eq."Mg24") then
         WS1 = 0.
      end if
c
      if (target.eq."Na23") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WS1 = 0.0253345
            else if (term.eq."y1") then
               WS1 = -0.0750847
            else if (term.eq."y2") then
               WS1 = 0.100235
            else if (term.eq."y3") then
               WS1 = -0.0384261
            else if (term.eq."y4") then
               WS1 = 0.00466396
            else
               WS1 = 0.
            end if         
         else if ((i.eq.1).and.(j.eq.1)) then
            if (term.eq."y0") then
               WS1 = 0.0183515
            else if (term.eq."y1") then
               WS1 = -0.0567009
            else if (term.eq."y2") then
               WS1 = 0.0887794
            else if (term.eq."y3") then
               WS1 = -0.0374699
            else if (term.eq."y4") then
               WS1 = 0.00477955
            else
               WS1 = 0.
            end if
         else 
            if (term.eq."y0") then
               WS1 = 0.0215622
            else if (term.eq."y1") then
               WS1 = -0.0652627
            else if (term.eq."y2") then
               WS1 = 0.0941439
            else if (term.eq."y3") then
               WS1 = -0.0379511
            else if (term.eq."y4") then
               WS1 = 0.00472138
            else
               WS1 = 0.
            end if
         end if
      end if
c
      if (target.eq."Ne20") then
         WS1 = 0.
      end if
c
      if (target.eq."O16") then
         WS1 = 0.
      end if
c     
      if (target.eq."N14") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WS1 = 0.06735548727
            else if (term.eq."y1") then
               WS1 = -0.1903833227
            else if (term.eq."y2") then
               WS1 = 0.1345317622
            else
               WS1 = 0.
            end if
         else
            WS1 = 0.
         end if
      end if
c
      if (target.eq."C12") then
         WS1 = 0.
      end if
c
      if (target.eq."He4") then
         WS1 = 0.
      end if
c
      if (target.eq."He3") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WS1 = 0.0795775
            else
               WS1 = 0.
            end if
         else if ((i.eq.1).and.(j.eq.1)) then
            if (term.eq."y0") then
               WS1 = 0.0795775
            else
               WS1 = 0.
            end if
         else
            if (term.eq."y0") then
               WS1 = -0.0795775
            else
               WS1 = 0.
            end if
         end if
      end if
c
ccc...
ccc*H
ccc...
      if (target.eq."H") then
         if (term.eq."y0") then
            WS1 = 0.07957747154594766!/exp(2.*y)
         else
            WS1 = 0.
         end if
      end if
c
      end function WS1
