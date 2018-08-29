      function WP2(i,j,target,term)
      implicit none
      real y,WP2
      integer i,j
      character (len=4) :: target
      character (len=2) :: term
      !include 'dsddcom.h'
c
      if (target.eq."Ni58") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WP2 = 5.4697
            else if (term.eq."y1") then
               WP2 = -8.75152
            else if (term.eq."y2") then
               WP2 = 4.88454
            else if (term.eq."y3") then
               WP2 = -1.10715
            else if (term.eq."y4") then
               WP2 = 0.0875404
            else
               WP2 = 0.
            end if
         else if ((i.eq.1).and.(j.eq.1)) then
            if (term.eq."y0") then
               WP2 = 0.00977975
            else if (term.eq."y1") then
               WP2 = -0.0156476
            else if (term.eq."y2") then
               WP2 = 0.0136707
            else if (term.eq."y3") then
               WP2 = -0.00592935
            else if (term.eq."y4") then
               WP2 = 0.00140426
            else
               WP2 = 0.
            end if
         else
            if (term.eq."y0") then
               WP2 = -0.231284
            else if (term.eq."y1") then
               WP2 = 0.370054
            else if (term.eq."y2") then
               WP2 = -0.264922
            else if (term.eq."y3") then
               WP2 = 0.0935201
            else if (term.eq."y4") then
               WP2 = -0.0110873
            else
               WP2 = 0.
            end if
         end if
      end if
c
      if (target.eq."Fe56") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WP2 = 4.228719571554119
            else if (term.eq."y1") then
               WP2 = -6.765952099008808
            else if (term.eq."y2") then
               WP2 = 3.7906717589007934
            else if (term.eq."y3") then
               WP2 = -0.8674325849707768
            else if (term.eq."y4") then
               WP2 = 0.06950603470249758
            else
               WP2 = 0.
            end if
         else if ((i.eq.1).and.(j.eq.1)) then
            if (term.eq."y0") then
               WP2 = 0.1433777108801172
            else if (term.eq."y1") then
               WP2 = -0.22940433740818778
            else if (term.eq."y2") then
               WP2 = 0.14460626294197715
            else if (term.eq."y3") then
               WP2 = -0.0422756223829621
            else if (term.eq."y4") then
               WP2 = 0.0048692089588925555
            else
               WP2 = 0.
            end if
         else
            if (term.eq."y0") then
               WP2 = -0.778655335898611
            else if (term.eq."y1") then
               WP2 = 1.2458486096667816
            else if (term.eq."y2") then
               WP2 = -0.7416613747520316
            else if (term.eq."y3") then
               WP2 = 0.1946575349056471
            else if (term.eq."y4") then
               WP2 = -0.01839672272092228
            else
               WP2 = 0.
            end if
         end if
      end if
c
      if (target.eq."Ca40") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WP2 = 0.006481658707
            else if (term.eq."y1") then
               WP2 = 0.008644293538
            else if (term.eq."y2") then
               WP2 = 0.003870408522
            else if (term.eq."y3") then
               WP2 = 0.0006590141866
            else if (term.eq."y4") then
               WP2 = 0.00003767182629
            else
               WP2 = 0.
            end if
         else
            WP2 = 0.
         end if
      end if

c
      if (target.eq."Ar40") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WP2 = 0.2996292793379571
            else if (term.eq."y1") then
               WP2 = -0.37379828709206
            else if (term.eq."y2") then
               WP2 = 0.15489450744016164
            else if (term.eq."y3") then
               WP2 = -0.023898303505580407
            else if (term.eq."y4") then
               WP2 = 0.001224739657053445
            else
               WP2 = 0.
            end if
         else if ((i.eq.1).and.(j.eq.1)) then
            if (term.eq."y0") then
               WP2 = 0.004149993889349473
            else if (term.eq."y1") then
               WP2 = -0.01814740982677668
            else if (term.eq."y2") then
               WP2 = 0.024075498682688812
            else if (term.eq."y3") then
               WP2 = -0.009262635721693066
            else if (term.eq."y4") then
               WP2 = 0.001081153434277832
            else
               WP2 = 0.
            end if
         else
            if (term.eq."y0") then
               WP2 = -0.03526272363733006
            else if (term.eq."y1") then
               WP2 = 0.09909552091617986
            else if (term.eq."y2") then
               WP2 = -0.06834531659689337
            else if (term.eq."y3") then
               WP2 = 0.016156146438124963
            else if (term.eq."y4") then
               WP2 = -0.0011507091232451346
            else
               WP2 = 0.
            end if
         end if
      end if
c
      if (target.eq."S32") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WP2 = 0.4787138511
            else if (term.eq."y1") then
               WP2 = -0.3829708692
            else if (term.eq."y2") then
               WP2 = 0.07659413150
            else
               WP2 = 0.
            end if
         else
            WP2 = 0.
         end if
      end if
c
      if (target.eq."Si28") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WP2 = 0.4619395950
            else if (term.eq."y1") then
               WP2 = -0.3695514682
            else if (term.eq."y2") then
               WP2 = 0.07391025207
            else
               WP2 = 0.
            end if
         else
            WP2 = 0.
         end if
      end if
c
      if (target.eq."Al27") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WP2 = 2.804982677628415
            else if (term.eq."y1") then
               WP2 = -2.2430558399020724
            else if (term.eq."y2") then
               WP2 = 0.45549083810832247
            else
               WP2 = 0.
            end if
         else if ((i.eq.1).and.(j.eq.1)) then
            if (term.eq."y0") then
               WP2 = 0.021492965144796174
            else if (term.eq."y1") then
               WP2 = -0.015615900257270558
            else if (term.eq."y2") then
               WP2 = 0.005968858233116949
            else
               WP2 = 0.
            end if
         else
            if (term.eq."y0") then
               WP2 = -0.18041651507963524
            else if (term.eq."y1") then
               WP2 = 0.1373889421134373
            else if (term.eq."y2") then
               WP2 = -0.023961463225677174
            else
               WP2 = 0.
            end if
         end if
      end if
c
      if (target.eq."Mg24") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WP2 = 0.1630103514
            else if (term.eq."y1") then
               WP2 = -0.1304081576
            else if (term.eq."y2") then
               WP2 = 0.02608160681
            else
               WP2 = 0.
            end if
         else
            WP2 = 0.
         end if
      end if
c
      if (target.eq."Na23") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WP2 = 0.612149
            else if (term.eq."y1") then
               WP2 = -0.49308
            else if (term.eq."y2") then
               WP2 = 0.107832
            else
               WP2 = 0.
            end if         
         else if ((i.eq.1).and.(j.eq.1)) then
            if (term.eq."y0") then
               WP2 = 0.00940911
            else if (term.eq."y1") then
               WP2 = -0.00747826
            else if (term.eq."y2") then
               WP2 = 0.00163204
            else
               WP2 = 0.
            end if
         else 
            if (term.eq."y0") then
               WP2 = -0.075893
            else if (term.eq."y1") then
               WP2 = 0.060682
            else if (term.eq."y2") then
               WP2 = -0.0110124
            else
               WP2 = 0.
            end if
         end if
      end if
c
      if (target.eq."Ne20") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WP2 = 0.02175494505
            else if (term.eq."y1") then
               WP2 = -0.01740391092
            else if (term.eq."y2") then
               WP2 = 0.003480773161
            else
               WP2 = 0.
            end if
         else
            WP2 = 0.
         end if
      end if
c
      if (target.eq."O16") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WP2 = 0.0004372029516
            else if (term.eq."y1") then
               WP2 = -0.0002388726962
            else if (term.eq."y2") then
               WP2 = 0.00003262796188
            else
               WP2 = 0.
            end if
         else
            WP2 = 0.
         end if
      end if
c
      if (target.eq."N14") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WP2 = 0.0905047989234487
            else
               WP2 = 0.
            end if
         else
            WP2 = 0.
         end if
      end if
c
      if (target.eq."C12") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WP2 = 0.04808048103491103
            else
               WP2 = 0.
            end if
         else
            WP2 = 0.
         end if
      end if
c     
      if (target.eq."He4") then
         WP2 = 0.
      end if
c     
      if (target.eq."He3") then
         WP2 = 0.
      end if
c
      if (target.eq."H") then
          WP2 = 0.
      end if
c
      end function WP2
