      function WM(i,j,target,term)
      implicit none
      real y,WM
      integer i,j
      character (len=4) :: target
      character (len=2) :: term
      !include 'dsddcom.h'
c
      if (target.eq."Ni58") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WM = 66.9246
            else if (term.eq."y1") then
               WM = -175.389
            else if (term.eq."y2") then
               WM = 169.877
            else if (term.eq."y3") then
               WM = -76.127
            else if (term.eq."y4") then
               WM = 16.6597
            else if (term.eq."y5") then
               WM = -1.6839
            else if (term.eq."y6") then
               WM = 0.0628067
            else
               WM = 0.
            end if
         else if ((i.eq.1).and.(j.eq.1)) then
            if (term.eq."y0") then
               WM = 0.0795762
            else if (term.eq."y1") then
               WM = -0.318305
            else if (term.eq."y2") then
               WM = 0.548985
            else if (term.eq."y3") then
               WM = -0.503018
            else if (term.eq."y4") then
               WM = 0.250492
            else if (term.eq."y5") then
               WM = -0.0603789
            else if (term.eq."y6") then
               WM = 0.00545169
            else
               WM = 0.
            end if
         else
            if (term.eq."y0") then
               WM = -2.30773
            else if (term.eq."y1") then
               WM = 7.63937
            else if (term.eq."y2") then
               WM = -10.3404
            else if (term.eq."y3") then
               WM = 6.95311
            else if (term.eq."y4") then
               WM = -2.30652
            else if (term.eq."y5") then
               WM = 0.350525
            else if (term.eq."y6") then
               WM = -0.0185041
            else
               WM = 0.
            end if
         end if
      end if
c
      if (target.eq."Fe56") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WM = 62.388777015508765
            else if (term.eq."y1") then
               WM = -160.42824560643984
            else if (term.eq."y2") then
               WM = 152.6436768887057
            else if (term.eq."y3") then
               WM = -67.2779398120723
            else if (term.eq."y4") then
               WM = 14.478025411058926
            else if (term.eq."y5") then
               WM = -1.4366495961973593
            else if (term.eq."y6") then
               WM = 0.052529136733736465
            else
               WM = 0.
            end if
         else if ((i.eq.1).and.(j.eq.1)) then
            if (term.eq."y0") then
               WM = 0.31830868659112194
            else if (term.eq."y1") then
               WM = -1.2732347463644877
            else if (term.eq."y2") then
               WM = 1.9918832378895788
            else if (term.eq."y3") then
               WM = -1.5456167618581302
            else if (term.eq."y4") then
               WM = 0.6222644086859521
            else if (term.eq."y5") then
               WM = -0.12227728761714192
            else if (term.eq."y6") then
               WM = 0.009215248417080017
            else
               WM = 0.
            end if
         else
            if (term.eq."y0") then
               WM = -4.456331413823822
            else if (term.eq."y1") then
               WM = 14.642230425862513
            else if (term.eq."y2") then
               WM = -18.257941019835176
            else if (term.eq."y3") then
               WM = 10.891907564433797
            else if (term.eq."y4") then
               WM = -3.2295985806489713
            else if (term.eq."y5") then
               WM = 0.44683638613875426
            else if (term.eq."y6") then
               WM = -0.022001569128954065
            else
               WM = 0.
            end if
         end if
      end if
c
      if (target.eq."Ca40") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WM = 31.82985878
            else if (term.eq."y1") then
               WM = -64.06332409
            else if (term.eq."y2") then
               WM = 45.25265043
            else if (term.eq."y3") then
               WM = -13.14661048
            else if (term.eq."y4") then
               WM = 1.377492267
            else if (term.eq."y5") then
               WM = -0.009441496634
            else if (term.eq."y6") then
               WM = 0.00001674303391
            else
               WM = 0.
            end if
         else
            WM = 0. 
         end if
      end if
c
      if (target.eq."Ar40") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WM = 31.82942047318682
            else if (term.eq."y1") then
               WM = -65.96177301433545
            else if (term.eq."y2") then
               WM = 48.58344700506161
            else if (term.eq."y3") then
               WM = -15.193970945783072
            else if (term.eq."y4") then
               WM = 1.9035952169249877
            else if (term.eq."y5") then
               WM = 0.05958862082938618
            else if (term.eq."y6") then
               WM = 0.0005443287364681978
            else
               WM = 0.
            end if
         else if ((i.eq.1).and.(j.eq.1)) then
            if (term.eq."y0") then
               WM = 0.31830397368380947
            else if (term.eq."y1") then
               WM = -1.0652376196047997
            else if (term.eq."y4") then
               WM = 0.14161793642799148
            else if (term.eq."y5") then
               WM =  - 0.013879710940695546
            else if (term.eq."y6") then
               WM = 0.0004805126374568142
            else
               WM = 0.
            end if
         else
            if (term.eq."y0") then
               WM = -3.1829908917036316
            else if (term.eq."y1") then
               WM = 8.62424562744593
            else if (term.eq."y2") then
               WM = -8.025394612628931
            else if (term.eq."y3") then
               WM = 3.1931619505667515
            else if (term.eq."y4") then
               WM = -0.5544673259783914
            else if (term.eq."y5") then
               WM = 0.03537969449718757
            else if (term.eq."y6") then
               WM = -0.0005114262769978376
            else
               WM = 0.
            end if
         end if
      end if
c
      if (target.eq."S32") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WM = 20.37154292
            else if (term.eq."y1") then
               WM = -37.34770103
            else if (term.eq."y2") then
               WM = 23.99417716
            else if (term.eq."y3") then
               WM = 6.303472873
            else if (term.eq."y4") then
               WM = 0.5803045188
            else
               WM = 0.
            end if
         else
            WM = 0.
         end if
      end if
c
      if (target.eq."Si28") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WM = 15.59657272
            else if (term.eq."y1") then
               WM = -26.73668448
            else if (term.eq."y2") then
               WM = 15.65057555
            else if (term.eq."y3") then
               WM = 3.593209109
            else if (term.eq."y4") then
               WM = 0.2816949468
            else
               WM = 0.
            end if
         else
            WM = 0.
         end if
      end if
c
      if (target.eq."Al27") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WM = 87.01458668430054
            else if (term.eq."y1") then
               WM = -146.0969173747615
            else if (term.eq."y2") then
               WM = 83.53672713691884
            else if (term.eq."y3") then
               WM = -18.598141659854168
            else if (term.eq."y4") then
               WM = 1.4344582399866814
            else
               WM = 0.
            end if
         else if ((i.eq.1).and.(j.eq.1)) then
            if (term.eq."y0") then
               WM = 0.11936637441693228
            else if (term.eq."y1") then
               WM = -0.31831033177848606
            else if (term.eq."y2") then
               WM = 0.33729069642027937
            else if (term.eq."y3") then
               WM = -0.1325262656678264
            else if (term.eq."y4") then
               WM = 0.01815501437926362
            else
               WM= 0.
            end if
         else
            if (term.eq."y0") then
               WM = -3.222827288871286
            else if (term.eq."y1") then
               WM = 7.002655937137329
            else if (term.eq."y2") then
               WM = -4.927561674517203
            else if (term.eq."y3") then
               WM = 1.3358703438729338
            else if (term.eq."y4") then
               WM = -0.11523996342329276
            else
               WM = 0.
            end if
         end if
      end if
c
      if (target.eq."Mg24") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WM = 11.45912908
            else if (term.eq."y1") then
               WM = -17.82529183
            else if (term.eq."y2") then
               WM = 9.310979256
            else if (term.eq."y3") then
               WM = -1.850276129
            else if (term.eq."y4") then
               WM = 0.1234671450
            else
               WM = 0.
            end if
         else
            WM = 0.
         end if
      end if
c
      if (target.eq."Na23") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WM = 42.0965
            else if (term.eq."y1") then
               WM = -63.4498
            else if (term.eq."y2") then
               WM = 32.5913
            else if (term.eq."y3") then
               WM = -6.57878
            else if (term.eq."y4") then
               WM = 0.483166
            else
               WM = 0.
            end if         
         else if ((i.eq.1).and.(j.eq.1)) then
            if (term.eq."y0") then
               WM = 0.0795776
            else if (term.eq."y1") then
               WM = -0.212207
            else if (term.eq."y2") then
               WM = 0.182941
            else if (term.eq."y3") then
               WM = -0.0543892
            else if (term.eq."y4") then
               WM = 0.00523012
            else
               WM = 0.
            end if
         else 
            if (term.eq."y0") then
               WM = -1.83028
            else if (term.eq."y1") then
               WM = 3.81972
            else if (term.eq."y2") then
               WM = -2.50445
            else if (term.eq."y3") then
               WM = 0.597822
            else if (term.eq."y4") then
               WM = -0.04545
            else
               WM = 0.
            end if
         end if
      end if
c
      if (target.eq."Ne20") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WM = 7.957741500
            else if (term.eq."y1") then
               WM = -10.61031393
            else if (term.eq."y2") then
               WM = 4.709038740
            else if (term.eq."y3") then
               WM = -0.7815128203
            else if (term.eq."y4") then
               WM = 0.04317233854
            else
               WM = 0.
            end if
         else
            WM = 0.
         end if
      end if

c     
      if (target.eq."O16") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WM = 5.092948703
            else if (term.eq."y1") then
               WM = -5.157409574
            else if (term.eq."y2") then
               WM = 1.331453166
            else if (term.eq."y3") then
               WM = -0.01305394998
            else if (term.eq."y4") then
               WM = 0.00003262796188
            else
               WM = 0.
            end if
         else
            WM = 0.
         end if
      end if
c
      if (target.eq."N14") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WM = 11.697894775374504
            else if (term.eq."y1") then
               WM = -11.14085547491048
            else if (term.eq."y2") then
               WM = 2.6757362417126016
            else
               WM = 0.
            end if
         else
            WM = 0.
         end if
      end if   
c
      if (target.eq."C12") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WM = 2.864782843
            else if (term.eq."y1") then
               WM = -2.546472275
            else if (term.eq."y2") then
               WM = 0.5658824250
            else
               WM = 0.
            end if
         else
            WM = 0.
         end if
      end if
c
      if (target.eq."He4") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WM = 0.31830988618379064
            else
               WM = 0.
            end if
         else
            WM = 0.
         end if   
      end if
c
      if (target.eq."He3") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               WM = 0.358099
            else
               WM = 0.
            end if
         else if ((i.eq.1).and.(j.eq.1)) then
            if (term.eq."y0") then
               WM = 0.0397887
            else
               WM = 0.
            end if
         else
            if (term.eq."y0") then
               WM = 0.119366
            else
               WM = 0.
            end if
         end if
      end if
c
ccc...
ccc*H
ccc...
      if (target.eq."H") then
         if (term.eq."y0") then
            WM = 0.039788735772973836!/exp(2.*y)
         else
            WM = 0.
         end if
      end if 
c
      end function WM
