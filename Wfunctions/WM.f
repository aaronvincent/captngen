      function response_nuc_m(i,j,target,term)
      implicit none
      real y,response_nuc_m
      integer i,j
      character (len=4) :: target
      character (len=2) :: term
      !include 'dsddcom.h'
c
      if (target.eq."Ni58") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               response_nuc_m = 66.9246
            else if (term.eq."y1") then
               response_nuc_m = -175.389
            else if (term.eq."y2") then
               response_nuc_m = 169.877
            else if (term.eq."y3") then
               response_nuc_m = -76.127
            else if (term.eq."y4") then
               response_nuc_m = 16.6597
            else if (term.eq."y5") then
               response_nuc_m = -1.6839
            else if (term.eq."y6") then
               response_nuc_m = 0.0628067
            else
               response_nuc_m = 0.
            end if
         else if ((i.eq.1).and.(j.eq.1)) then
            if (term.eq."y0") then
               response_nuc_m = 0.0795762
            else if (term.eq."y1") then
               response_nuc_m = -0.318305
            else if (term.eq."y2") then
               response_nuc_m = 0.548985
            else if (term.eq."y3") then
               response_nuc_m = -0.503018
            else if (term.eq."y4") then
               response_nuc_m = 0.250492
            else if (term.eq."y5") then
               response_nuc_m = -0.0603789
            else if (term.eq."y6") then
               response_nuc_m = 0.00545169
            else
               response_nuc_m = 0.
            end if
         else
            if (term.eq."y0") then
               response_nuc_m = -2.30773
            else if (term.eq."y1") then
               response_nuc_m = 7.63937
            else if (term.eq."y2") then
               response_nuc_m = -10.3404
            else if (term.eq."y3") then
               response_nuc_m = 6.95311
            else if (term.eq."y4") then
               response_nuc_m = -2.30652
            else if (term.eq."y5") then
               response_nuc_m = 0.350525
            else if (term.eq."y6") then
               response_nuc_m = -0.0185041
            else
               response_nuc_m = 0.
            end if
         end if
      end if
c
      if (target.eq."Fe56") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               response_nuc_m = 62.388777015508765
            else if (term.eq."y1") then
               response_nuc_m = -160.42824560643984
            else if (term.eq."y2") then
               response_nuc_m = 152.6436768887057
            else if (term.eq."y3") then
               response_nuc_m = -67.2779398120723
            else if (term.eq."y4") then
               response_nuc_m = 14.478025411058926
            else if (term.eq."y5") then
               response_nuc_m = -1.4366495961973593
            else if (term.eq."y6") then
               response_nuc_m = 0.052529136733736465
            else
               response_nuc_m = 0.
            end if
         else if ((i.eq.1).and.(j.eq.1)) then
            if (term.eq."y0") then
               response_nuc_m = 0.31830868659112194
            else if (term.eq."y1") then
               response_nuc_m = -1.2732347463644877
            else if (term.eq."y2") then
               response_nuc_m = 1.9918832378895788
            else if (term.eq."y3") then
               response_nuc_m = -1.5456167618581302
            else if (term.eq."y4") then
               response_nuc_m = 0.6222644086859521
            else if (term.eq."y5") then
               response_nuc_m = -0.12227728761714192
            else if (term.eq."y6") then
               response_nuc_m = 0.009215248417080017
            else
               response_nuc_m = 0.
            end if
         else
            if (term.eq."y0") then
               response_nuc_m = -4.456331413823822
            else if (term.eq."y1") then
               response_nuc_m = 14.642230425862513
            else if (term.eq."y2") then
               response_nuc_m = -18.257941019835176
            else if (term.eq."y3") then
               response_nuc_m = 10.891907564433797
            else if (term.eq."y4") then
               response_nuc_m = -3.2295985806489713
            else if (term.eq."y5") then
               response_nuc_m = 0.44683638613875426
            else if (term.eq."y6") then
               response_nuc_m = -0.022001569128954065
            else
               response_nuc_m = 0.
            end if
         end if
      end if
c
      if (target.eq."Ca40") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               response_nuc_m = 31.82985878
            else if (term.eq."y1") then
               response_nuc_m = -64.06332409
            else if (term.eq."y2") then
               response_nuc_m = 45.25265043
            else if (term.eq."y3") then
               response_nuc_m = -13.14661048
            else if (term.eq."y4") then
               response_nuc_m = 1.377492267
            else if (term.eq."y5") then
               response_nuc_m = -0.009441496634
            else if (term.eq."y6") then
               response_nuc_m = 0.00001674303391
            else
               response_nuc_m = 0.
            end if
         else
            response_nuc_m = 0. 
         end if
      end if
c
      if (target.eq."Ar40") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               response_nuc_m = 31.82942047318682
            else if (term.eq."y1") then
               response_nuc_m = -65.96177301433545
            else if (term.eq."y2") then
               response_nuc_m = 48.58344700506161
            else if (term.eq."y3") then
               response_nuc_m = -15.193970945783072
            else if (term.eq."y4") then
               response_nuc_m = 1.9035952169249877
            else if (term.eq."y5") then
               response_nuc_m = 0.05958862082938618
            else if (term.eq."y6") then
               response_nuc_m = 0.0005443287364681978
            else
               response_nuc_m = 0.
            end if
         else if ((i.eq.1).and.(j.eq.1)) then
            if (term.eq."y0") then
               response_nuc_m = 0.31830397368380947
            else if (term.eq."y1") then
               response_nuc_m = -1.0652376196047997
            else if (term.eq."y4") then
               response_nuc_m = 0.14161793642799148
            else if (term.eq."y5") then
               response_nuc_m =  - 0.013879710940695546
            else if (term.eq."y6") then
               response_nuc_m = 0.0004805126374568142
            else
               response_nuc_m = 0.
            end if
         else
            if (term.eq."y0") then
               response_nuc_m = -3.1829908917036316
            else if (term.eq."y1") then
               response_nuc_m = 8.62424562744593
            else if (term.eq."y2") then
               response_nuc_m = -8.025394612628931
            else if (term.eq."y3") then
               response_nuc_m = 3.1931619505667515
            else if (term.eq."y4") then
               response_nuc_m = -0.5544673259783914
            else if (term.eq."y5") then
               response_nuc_m = 0.03537969449718757
            else if (term.eq."y6") then
               response_nuc_m = -0.0005114262769978376
            else
               response_nuc_m = 0.
            end if
         end if
      end if
c
      if (target.eq."S32") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               response_nuc_m = 20.37154292
            else if (term.eq."y1") then
               response_nuc_m = -37.34770103
            else if (term.eq."y2") then
               response_nuc_m = 23.99417716
            else if (term.eq."y3") then
               response_nuc_m = -6.303472873
            else if (term.eq."y4") then
               response_nuc_m = 0.5803045188
            else
               response_nuc_m = 0.
            end if
         else
            response_nuc_m = 0.
         end if
      end if
c
      if (target.eq."Si28") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               response_nuc_m = 15.59657272
            else if (term.eq."y1") then
               response_nuc_m = -26.73668448
            else if (term.eq."y2") then
               response_nuc_m = 15.65057555
            else if (term.eq."y3") then
               response_nuc_m = -3.593209109
            else if (term.eq."y4") then
               response_nuc_m = 0.2816949468
            else
               response_nuc_m = 0.
            end if
         else
            response_nuc_m = 0.
         end if
      end if
c
      if (target.eq."Al27") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               response_nuc_m = 87.01458668430054
            else if (term.eq."y1") then
               response_nuc_m = -146.0969173747615
            else if (term.eq."y2") then
               response_nuc_m = 83.53672713691884
            else if (term.eq."y3") then
               response_nuc_m = -18.598141659854168
            else if (term.eq."y4") then
               response_nuc_m = 1.4344582399866814
            else
               response_nuc_m = 0.
            end if
         else if ((i.eq.1).and.(j.eq.1)) then
            if (term.eq."y0") then
               response_nuc_m = 0.11936637441693228
            else if (term.eq."y1") then
               response_nuc_m = -0.31831033177848606
            else if (term.eq."y2") then
               response_nuc_m = 0.33729069642027937
            else if (term.eq."y3") then
               response_nuc_m = -0.1325262656678264
            else if (term.eq."y4") then
               response_nuc_m = 0.01815501437926362
            else
               response_nuc_m= 0.
            end if
         else
            if (term.eq."y0") then
               response_nuc_m = -3.222827288871286
            else if (term.eq."y1") then
               response_nuc_m = 7.002655937137329
            else if (term.eq."y2") then
               response_nuc_m = -4.927561674517203
            else if (term.eq."y3") then
               response_nuc_m = 1.3358703438729338
            else if (term.eq."y4") then
               response_nuc_m = -0.11523996342329276
            else
               response_nuc_m = 0.
            end if
         end if
      end if
c
      if (target.eq."Mg24") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               response_nuc_m = 11.45912908
            else if (term.eq."y1") then
               response_nuc_m = -17.82529183
            else if (term.eq."y2") then
               response_nuc_m = 9.310979256
            else if (term.eq."y3") then
               response_nuc_m = -1.850276129
            else if (term.eq."y4") then
               response_nuc_m = 0.1234671450
            else
               response_nuc_m = 0.
            end if
         else
            response_nuc_m = 0.
         end if
      end if
c
      if (target.eq."Na23") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               response_nuc_m = 42.0965
            else if (term.eq."y1") then
               response_nuc_m = -63.4498
            else if (term.eq."y2") then
               response_nuc_m = 32.5913
            else if (term.eq."y3") then
               response_nuc_m = -6.57878
            else if (term.eq."y4") then
               response_nuc_m = 0.483166
            else
               response_nuc_m = 0.
            end if         
         else if ((i.eq.1).and.(j.eq.1)) then
            if (term.eq."y0") then
               response_nuc_m = 0.0795776
            else if (term.eq."y1") then
               response_nuc_m = -0.212207
            else if (term.eq."y2") then
               response_nuc_m = 0.182941
            else if (term.eq."y3") then
               response_nuc_m = -0.0543892
            else if (term.eq."y4") then
               response_nuc_m = 0.00523012
            else
               response_nuc_m = 0.
            end if
         else 
            if (term.eq."y0") then
               response_nuc_m = -1.83028
            else if (term.eq."y1") then
               response_nuc_m = 3.81972
            else if (term.eq."y2") then
               response_nuc_m = -2.50445
            else if (term.eq."y3") then
               response_nuc_m = 0.597822
            else if (term.eq."y4") then
               response_nuc_m = -0.04545
            else
               response_nuc_m = 0.
            end if
         end if
      end if
c
      if (target.eq."Ne20") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               response_nuc_m = 7.957741500
            else if (term.eq."y1") then
               response_nuc_m = -10.61031393
            else if (term.eq."y2") then
               response_nuc_m = 4.709038740
            else if (term.eq."y3") then
               response_nuc_m = -0.7815128203
            else if (term.eq."y4") then
               response_nuc_m = 0.04317233854
            else
               response_nuc_m = 0.
            end if
         else
            response_nuc_m = 0.
         end if
      end if

c     
      if (target.eq."O16") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               response_nuc_m = 5.092948703
            else if (term.eq."y1") then
               response_nuc_m = -5.157409574
            else if (term.eq."y2") then
               response_nuc_m = 1.331453166
            else if (term.eq."y3") then
               response_nuc_m = -0.01305394998
            else if (term.eq."y4") then
               response_nuc_m = 0.00003262796188
            else
               response_nuc_m = 0.
            end if
         else
            response_nuc_m = 0.
         end if
      end if
c
      if (target.eq."N14") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               response_nuc_m = 11.697894775374504
            else if (term.eq."y1") then
               response_nuc_m = -11.14085547491048
            else if (term.eq."y2") then
               response_nuc_m = 2.6757362417126016
            else
               response_nuc_m = 0.
            end if
         else
            response_nuc_m = 0.
         end if
      end if   
c
      if (target.eq."C12") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               response_nuc_m = 2.864782843
            else if (term.eq."y1") then
               response_nuc_m = -2.546472275
            else if (term.eq."y2") then
               response_nuc_m = 0.5658824250
            else
               response_nuc_m = 0.
            end if
         else
            response_nuc_m = 0.
         end if
      end if
c
      if (target.eq."He4") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               response_nuc_m = 0.31830988618379064
            else
               response_nuc_m = 0.
            end if
         else
            response_nuc_m = 0.
         end if   
      end if
c
      if (target.eq."He3") then
         if ((i.eq.0).and.(j.eq.0)) then
            if (term.eq."y0") then
               response_nuc_m = 0.358099
            else
               response_nuc_m = 0.
            end if
         else if ((i.eq.1).and.(j.eq.1)) then
            if (term.eq."y0") then
               response_nuc_m = 0.0397887
            else
               response_nuc_m = 0.
            end if
         else
            if (term.eq."y0") then
               response_nuc_m = 0.119366
            else
               response_nuc_m = 0.
            end if
         end if
      end if
c
ccc...
ccc*H
ccc...
      if (target.eq."H") then
         if (term.eq."y0") then
            response_nuc_m = 0.039788735772973836!/exp(2.*y)
         else
            response_nuc_m = 0.
         end if
      end if 
c
      end function response_nuc_m
