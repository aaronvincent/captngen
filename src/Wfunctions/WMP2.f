      function WMP2(i,j,target,term)
      implicit none
      real y,WMP2
      integer i,j
      character (len=4) :: target
      character (len=2) :: term
      !include 'dsddcom.h'
c
      if (target.eq."Ni58") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WMP2 = -19.1326
            else if (term.eq."y1") then
               WMP2 = 40.3764
            else if (term.eq."y2") then
               WMP2 = -30.3339
            else if (term.eq."y3") then
               WMP2 = 10.0435
            else if (term.eq."y4") then
               WMP2 = -1.4629
            else if (term.eq."y5") then
               WMP2 = 0.0741493
            else
               WMP2 = 0.
            end if
         else if ((i.eq.1).and.(j.eq.1)) then
            if (term.eq."y0") then
               WMP2 = -0.0278969
            else if (term.eq."y1") then
               WMP2 = 0.0781112
            else if (term.eq."y2") then
               WMP2 = -0.0956406
            else if (term.eq."y3") then
               WMP2 = 0.0607914
            else if (term.eq."y4") then
               WMP2 = -0.0211633
            else if (term.eq."y5") then
               WMP2 = 0.00276687
            else
               WMP2 = 0.
            end if
         else if ((i.eq.1).and.(j.eq.0)) then
            if (term.eq."y0") then
               WMP2 = 0.659741
            else if (term.eq."y1") then
               WMP2 = -1.84727
            else if (term.eq."y2") then
               WMP2 = 2.0953
            else if (term.eq."y3") then
               WMP2 = -1.10461
            else if (term.eq."y4") then
               WMP2 = 0.25912
            else if (term.eq."y5") then
               WMP2 = -0.0218459
            else
               WMP2 = 0.
            end if
         else
            if (term.eq."y0") then
               WMP2 = 0.809015
            else if (term.eq."y1") then
               WMP2 = -1.7073
            else if (term.eq."y2") then
               WMP2 = 1.48687
            else if (term.eq."y3") then
               WMP2 = -0.692274
            else if (term.eq."y4") then
               WMP2 = 0.145722
            else if (term.eq."y5") then
               WMP2 = -0.00939131
            else
               WMP2 = 0.
            end if
         end if
      end if
c
      if (target.eq."Fe56") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WMP2 = -16.24267965610378
            else if (term.eq."y1") then
               WMP2 = 33.877585537899655
            else if (term.eq."y2") then
               WMP2 = -25.234189173759155
            else if (term.eq."y3") then
               WMP2 = 8.304707912042549
            else if (term.eq."y4") then
               WMP2 = -1.2033352625021452
            else if (term.eq."y5") then
               WMP2 = 0.06042426665427828
            else
               WMP2 = 0.
            end if
         else if ((i.eq.1).and.(j.eq.1)) then
            if (term.eq."y0") then
               WMP2 = -0.21363139010148274
            else if (term.eq."y1") then
               WMP2 = 0.5981678922841512
            else if (term.eq."y2") then
               WMP2 = -0.6223379309750409
            else if (term.eq."y3") then
               WMP2 = 0.30801402118087706
            else if (term.eq."y4") then
               WMP2 = -0.07352110962374007
            else if (term.eq."y5") then
               WMP2 = 0.006698579711466187
            else
               WMP2 = 0.
            end if
         else if ((i.eq.1).and.(j.eq.0)) then
            if (term.eq."y0") then
               WMP2 = 1.160188852206212
            else if (term.eq."y1") then
               WMP2 = -3.248528893797909
            else if (term.eq."y2") then
               WMP2 = 3.3147296224308356
            else if (term.eq."y3") then
               WMP2 = -1.5426380083289353
            else if (term.eq."y4") then
               WMP2 = 0.3258327814311257
            else if (term.eq."y5") then
               WMP2 = -0.025308405249436387  
            else
               WMP2 = 0.
            end if
         else
            if (term.eq."y0") then
               WMP2 = 2.9908460396840395
            else if (term.eq."y1") then
               WMP2 = -6.238049396962003
            else if (term.eq."y2") then
               WMP2 = 4.814220389520166
            else if (term.eq."y3") then
               WMP2 = -1.7448311733293451
            else if (term.eq."y4") then
               WMP2 = 0.28812798565332887
            else if (term.eq."y5") then
               WMP2 = -0.015992977933668308
            else
               WMP2 = 0.
            end if
         end if
      end if

c
      if (target.eq."Ca40") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WMP2 = -0.45421391592141563
            else if (term.eq."y1") then
               WMP2 = 0.7599758847704111
            else if (term.eq."y2") then
               WMP2 = -0.432313511223993
            else if (term.eq."y3") then
               WMP2 = 0.0971138384810061
            else if (term.eq."y4") then
               WMP2 = -0.0073007938698465245
            else if (term.eq."y5") then
               WMP2 = 0.000025114550863123492
            else
               WMP2 = 0.
            end if
         else
            WMP2 = 0.
         end if
      end if
c
      if (target.eq."Ar40") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WMP2 = -3.0882076222504518
            else if (term.eq."y1") then
               WMP2 = 5.126252058965625
            else if (term.eq."y2") then
               WMP2 = -2.892480914451324
            else if (term.eq."y3") then
               WMP2 = 0.6533856099199644
            else if (term.eq."y4") then
               WMP2 = -0.05265756679056644
            else if (term.eq."y5") then
               WMP2 = 0.0008164931047022968
            else
               WMP2 = 0.
            end if
         else if ((i.eq.1).and.(j.eq.1)) then
            if (term.eq."y0") then
               WMP2 = -0.03634500716389358
            else if (term.eq."y1") then
               WMP2 = 0.14028229828121988
            else if (term.eq."y2") then
               WMP2 = -0.17191662947208944
            else if (term.eq."y3") then
               WMP2 = 0.07704558439012923
            else if (term.eq."y4") then
               WMP2 = -0.013497328446086023
            else if (term.eq."y5") then
               WMP2 = 0.0007207689561852215
            else
               WMP2 = 0.
            end if 
         else if ((i.eq.1).and.(j.eq.0)) then
            if (term.eq."y0") then
               WMP2 = 0.3088255012871975
            else if (term.eq."y1") then
               WMP2 = -0.7093939111200537
            else if (term.eq."y2") then
               WMP2 = 0.5153783625844404
            else if (term.eq."y3") then
               WMP2 = -0.15313415178942127
            else if (term.eq."y4") then
               WMP2 = 0.01856407628805948
            else if (term.eq."y5") then
               WMP2 = -0.0007671394154967563
            else
               WMP2 = 0.
            end if
         else
            if (term.eq."y0") then
               WMP2 = 0.36344449433890386
            else if (term.eq."y1") then
               WMP2 = -1.1712414221110667
            else if (term.eq."y2") then
               WMP2 = 1.0911694222159676
            else if (term.eq."y3") then
               WMP2 = -0.37359193648228745
            else if (term.eq."y4") then
               WMP2 = 0.045276229749805184
            else if (term.eq."y5") then
               WMP2 = -0.0007671394154967566
            else
               WMP2 = 0.
            end if
         end if
      end if
c     
      if (target.eq."S32") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WMP2 = -3.122841616248989
            else if (term.eq."y1") then
               WMP2 = 4.111730986960552
            else if (term.eq."y2") then
               WMP2 = -1.6721045582756733
            else if (term.eq."y3") then
               WMP2 = 0.21082675499009432
            else
               WMP2 = 0.
            end if
         else
            WMP2 = 0.
         end if
      end if
c
      if (target.eq."Si28") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WMP2 = -2.684152471299983
            else if (term.eq."y1") then
               WMP2 = 3.3743369182092318
            else if (term.eq."y2") then
               WMP2 = -1.2809999804387553
            else if (term.eq."y3") then
               WMP2 = 0.1442918726957604
            else
               WMP2 = 0.
            end if
         else
            WMP2 = 0.
         end if
      end if
c
      if (target.eq."Al27") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WMP2 = -15.622754546229208
            else if (term.eq."y1") then
               WMP2 = 19.358850503552464
            else if (term.eq."y2") then
               WMP2 = -7.232344450625174
            else if (term.eq."y3") then
               WMP2 = 0.7970503749890764
            else
               WMP2 = 0.
            end if
         else if ((i.eq.1).and.(j.eq.1)) then
            if (term.eq."y0") then
               WMP2 = -0.03707943785718349
            else if (term.eq."y1") then
               WMP2 = 0.08525450185717358
            else if (term.eq."y2") then
               WMP2 = -0.04492839338037083
            else if (term.eq."y3") then
               WMP2 = 0.008669920653769503
            else
               WMP2 = 0
            end if
         else if ((i.eq.1).and.(j.eq.0)) then
            if (term.eq."y0") then
               WMP2 = 0.5786321764832287
            else if (term.eq."y1") then
               WMP2 = -1.0043777009381476
            else if (term.eq."y2") then
               WMP2 = 0.4912519971422017
            else if (term.eq."y3") then
               WMP2 = -0.07306931335081507
            else
               WMP2 = 0.
            end if 
         else
            if (term.eq."y0") then
               WMP2 = 1.0011246866285106
            else if (term.eq."y1") then
               WMP2 = -1.1593444477399983
            else if (term.eq."y2") then
               WMP2 = 0.40275041726200134
            else if (term.eq."y3") then
               WMP2 = -0.03649522438128564
            else
               WMP2 = 0.
            end if
         end if
      end if
c
      if (target.eq."Mg24") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WMP2 = -1.3667321088909605
            else if (term.eq."y1") then
               WMP2 = 1.6097049910493944
            else if (term.eq."y2") then
               WMP2 = -0.5670722872166495
            else if (term.eq."y3") then
               WMP2 = 0.05674699578152061
            else
               WMP2 = 0.
            end if
         else
            WMP2 = 0.
         end if
      end if
c
      if (target.eq."Na23") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WMP2 = -5.07498
            else if (term.eq."y1") then
               WMP2 = 5.86765
            else if (term.eq."y2") then
               WMP2 = -2.09908
            else if (term.eq."y3") then
               WMP2 = 0.226345
            else
               WMP2 = 0.
            end if         
         else if ((i.eq.1).and.(j.eq.1)) then
            if (term.eq."y0") then
               WMP2 = -0.0273574
            else if (term.eq."y1") then
               WMP2 = 0.0474719
            else if (term.eq."y2") then
               WMP2 = -0.0213121
            else if (term.eq."y3") then
               WMP2 = 0.00280825
            else
               WMP2 = 0.
            end if
         else if ((i.eq.1).and.(j.eq.0)) then
            if (term.eq."y0") then
               WMP2 = 0.220651
            else if (term.eq."y1") then
               WMP2 = -0.382932
            else if (term.eq."y2") then
               WMP2 = 0.17682
            else if (term.eq."y3") then
               WMP2 = -0.0226015
            else
               WMP2 = 0.
            end if
         else 
            if (term.eq."y0") then
               WMP2 = 0.62922
            else if (term.eq."y1") then
               WMP2 = -0.727336
            else if (term.eq."y2") then
               WMP2 = 0.243236
            else if (term.eq."y3") then
               WMP2 = -0.0210943
            else
               WMP2 = 0.
            end if
         end if
      end if
c
      if (target.eq."Ne20") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WMP2 = -0.4160771913844551
            else if (term.eq."y1") then
               WMP2 = 0.44381502814597296
            else if (term.eq."y2") then
               WMP2 = -0.14160010905766138
            else if (term.eq."y3") then
               WMP2 = 0.012258593610350725
            else
               WMP2 = 0.
            end if
         else
            WMP2 = 0.
         end if
      end if
c
      if (target.eq."O16") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WMP2 = -0.04718741575882355
            else if (term.eq."y1") then
               WMP2 = 0.03678312344365378
            else if (term.eq."y2") then
               WMP2 = -0.0066464113401998055
            else if (term.eq."y3") then
               WMP2 = 0.00003262796187666331
            else
               WMP2 = 0.
            end if
         else
            WMP2 = 0.
         end if
      end if
c
      if (target.eq."N14") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WMP2 = -1.0241365716917015
            else if (term.eq."y1") then
               WMP2 = 0.4832668390680021
            else
               WMP2 = 0.
            end if
         else
            WMP2 = 0.
         end if
      end if
c
      if (target.eq."C12") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WMP2 = -0.3711335839642701
            else if (term.eq."y1") then
               WMP2 = 0.16494817126159667
            else
               WMP2 = 0.
            end if
         else
            WMP2 = 0.
         endif
      end if
c
      if (target.eq."He4") then
         WMP2 = 0.
      end if
c     
      if (target.eq."He3") then
         WMP2 = 0.
      end if
c
      if (target.eq."H") then
          WMP2 = 0.
      end if
c
      end function WMP2
