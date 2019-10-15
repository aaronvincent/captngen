      function WP1(i,j,target,term)
      implicit none
      real y,WP1
      integer i,j
      character (len=4) :: target
      character (len=2) :: term
      !include 'dsddcom.h'
c
      if (target.eq."Ni58") then
         WP1 = 0.
      end if
c
      if (target.eq."Fe56") then
         WP1 = 0.
      end if
c     
      if (target.eq."Ca40") then
         WP1 = 0.
      end if
c
      if (target.eq."Ar40") then
         WP1 = 0.
      end if
c
      if (target.eq."S32") then
         WP1 = 0.
      end if
c
      if (target.eq."Si28") then
         WP1 = 0.
      end if
c
      if (target.eq."Al27") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WP1 = 0.00006807028943376759
            else if (term.eq."y1") then
               WP1 = -0.0003766815038380074
            else if (term.eq."y2") then
               WP1 = 0.00340250547302469
C             else if (term.eq."y3") then
C                WP1 = -4.747215506787354e-20
C             else if (term.eq."y4") then
C                WP1 = 2.4413135656776747e-36
            else
               WP1 = 0.
            end if
         else if ((i.eq.1).and.(j.eq.1)) then
            if (term.eq."y0") then
               WP1 = 0.01496216106271638
            else if (term.eq."y1") then
               WP1 = -0.005633072474649519
            else if (term.eq."y2") then
               WP1 = 0.0044038530074910435
C             else if (term.eq."y3") then
C                WP1 = -7.126454225762893e-20
C             else if (term.eq."y4") then
C                WP1 = 2.3946974273275598e-36
            else
               WP1 = 0.
            end if
         else
            if (term.eq."y0") then
               WP1 = -0.0010091970244178066
            else if (term.eq."y1") then
               WP1 = 0.00298227900889512
            else if (term.eq."y2") then
               WP1 = 0.002815253422852913
C             else if (term.eq."y3") then
C                WP1 = 8.439579246990046e-20
C             else if (term.eq."y4") then
C                WP1 = -2.3946974273275598e-36
            else
               WP1 = 0.
            end if
         end if
      end if
c
      if (target.eq."Mg24") then
         WP1 = 0.
      end if
c
      if (target.eq."Na23") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WP1 = 0.000495589
            else if (term.eq."y1") then
               WP1 = -0.00010394
            else if (term.eq."y2") then
               WP1 = 0.00000544981
            else
               WP1 = 0.
            end if         
         else if ((i.eq.1).and.(j.eq.1)) then
            if (term.eq."y0") then
               WP1 = 0.00000616583
            else if (term.eq."y1") then
               WP1 = 0.00008381
            else if (term.eq."y2") then
               WP1 = 0.0002848
            else
               WP1 = 0.
            end if
         else 
            if (term.eq."y0") then
               WP1 = -0.0000552785
            else if (term.eq."y1") then
               WP1 = -0.000369894
            else if (term.eq."y2") then
               WP1 = 0.0000393968
            else
               WP1 = 0.
            end if
         end if
      end if
c
      if (target.eq."Ne20") then
         WP1 = 0.
      end if
c
      if (target.eq."O16") then
         WP1 = 0.
      end if
c
      if (target.eq."N14") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WP1 = 0.001264316851770931
            else
               WP1 = 0.
            end if
         else
            WP1 = 0.
         end if
      end if
c
      if (target.eq."C12") then
         WP1 = 0.
      end if
c
      if (target.eq."He4") then
         WP1 = 0.
      end if
c
      if (target.eq."H") then
         WP1 = 0.
      end if
c
      end function WP1
