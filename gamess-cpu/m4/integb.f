**== dawerf.f
      function dawerf(y)
      implicit real*8 (a-h,o-z)
c
      real*8 c, h
      integer ifirst, ilast
      common/derfcm/c(246),ifirst(40),ilast(40),h
c
c
c        -----  routine evalues the dawson-error function.      -----
c
      x = abs(y)
      if(x.lt.10.0d+00) then
       xn = x/h
       nx = int(xn)
       nx = nx+1
       if = ifirst(nx)
       il = ilast(nx)
       t = c(il)
       kl = il-if
       do k=1,kl
        t = c(il-k)+x*t
       enddo
       dawerf = t
      else
       txt = 1.0d+00/(2.0d+00*x*x)
       tx = txt*x
       dawerf = tx*(1.0d+00+txt*(1.0d+00+txt*(3.0d+00+txt*(15.0d+00+
     *  txt*(105.0d+00+945.0d+00*txt)))))
      endif
c
      return
      end
**==dawert.f
      subroutine dawert
      implicit real*8 (a-h,o-z)
c        -----  routine allocates parameters specifying the     -----
c        -----  piecewise chebychev polynomial fit to the       -----
c        -----  dawson-error function.                          -----
c
c
      real*8 c, h
      integer ifirst, ilast
      common/derfcm/c(246),ifirst(40),ilast(40),h
c
      h = 0.25d+00
c... 0.00 @ x @ .25 interval no. 1 abs.error = 3.0831115438446d-
      ifirst( 1) =  1
      ilast ( 1) =  8
      c( 1) = -0.2227673357873d-15
      c( 2) = -0.3603041588494d-08
      c( 3)= 0.5641900599881d+00
      c( 4)=-0.1855716761667d-04
      c( 5)=-0.3758028089842d+00
      c( 6)=-0.2937263239086d-02
      c( 7)= 0.1586735773807d+00
      c( 8)=-0.3712362707925d-01
c....25 @ x @ .50 interval no. 2 abs.error=1.6200374375330d-
      ifirst( 2)= 9
      ilast ( 2)=16
      c( 9)= 0.1654923320373d-05
      c(10)=-0.4283673776908d-04
      c(11)= 0.5646780085194d+00
      c(12)=-0.3216632920741d-02
      c(13)=-0.3626661648377d+00
      c(14)=-0.3697874666458d-01
      c(15)= 0.2103001701512d+00
      c(16)=-0.7233313577516d-01
c... .50 @ x @  .75  interval no. 3  abs.error=2.6527224861184d-
      ifirst( 3)=17
      ilast ( 3)=24
      c(17)=-0.4909038422764d-03
      c(18)= 0.6169678713757d-02
      c(19)= 0.5309521767077d+00
      c(20)= 0.9902446632012d-01
      c(21)=-0.5497719612405d+00
      c(22)= 0.1699289458172d+00
      c(23)= 0.8215358267937d-01
      c(24)=-0.3800879631724d-01
c... .75 @ x @ 1.00  interval no. 4  abs.error=2.1188384380366d-
      ifirst( 4)=25
      ilast ( 4)=32
      c(25)=-0.8860184581188d-02
      c(26)= 0.8182679610811d-01
      c(27)= 0.2361899484096d+00
      c(28)= 0.7408443731368d+00
      c(29)=-0.1393552470603d+01
      c(30)= 0.8398504617092d+00
      c(31)=-0.2153157176716d+00
      c(32)= 0.1898293835776d-01
c...1.00 @ x @ 1.25  interval no. 5  abs.error=1.0054179711005d-
      ifirst( 5)=33
      ilast ( 5)=40
      c(33)=-0.2622864777452d-01
      c(34)= 0.2070481292278d+00
      c(35)=-0.1517729522946d+00
      c(36)= 0.1410375745664d+01
      c(37)=-0.2088589931697d+01
      c(38)= 0.1273821623264d+01
      c(39)=-0.3662054027830d+00
      c(40)= 0.4151758125850d-01
c...1.25 @ x @ 1.50  interval no. 6  abs.error=1.2743583965857d-
      ifirst( 6)=41
      ilast ( 6)=48
      c(41)= 0.7652753217902d-01
      c(42)=-0.3496979906426d+00
      c(43)= 0.1142731303255d+01
      c(44)=-0.2641065409391d+00
      c(45)=-0.7870989692331d+00
      c(46)= 0.6659477061060d+00
      c(47)=-0.2082257760423d+00
      c(48)= 0.2389272621700d-01
c...1.50 @ x @ 1.75  interval no. 7  abs.error=3.2651215065016d-
      ifirst( 7)=49
      ilast ( 7)=55
      c(49)= 0.4634935721533d+00
      c(50)=-0.2162775089324d+01
      c(51)= 0.4787745789535d+01
      c(52)=-0.4340108008393d+01
      c(53)= 0.1951017416635d+01
      c(54)=-0.4390388356211d+00
      c(55)= 0.3981846570969d-01
c...1.75 @ x @ 2.00  interval no. 8  abs.error=8.9706020389713d-
      ifirst( 8)=56
      ilast ( 8)=63
      c(56)= 0.9892135872598d+00
      c(57)=-0.4311631180681d+01
      c(58)= 0.8556855041877d+01
      c(59)=-0.8017503990110d+01
      c(60)= 0.4106390132111d+01
      c(61)=-0.1197913275632d+01
      c(62)= 0.1884269501482d+00
      c(63)=-0.1248598098755d-01
c...2.00 @ x @ 2.25  interval no. 9  abs.error=3.1388225352202d-
      ifirst( 9)=64
      ilast ( 9)=71
      c(64)= 0.4508582229747d+00
      c(65)=-0.2475021032343d+01
      c(66)= 0.5870952943693d+01
      c(67)=-0.5834804573043d+01
      c(68)= 0.3041872702481d+01
      c(69)=-0.8863349469112d+00
      c(70)= 0.1377496059452d+00
      c(71)=-0.8952617645264d-02
c...2.25 @ x @ 2.50  interval no. 10  abs.error=3.5207392556913d-
      ifirst(10)=72
      ilast (10)=78
      c(72)=-0.2479150925527d+01
      c(73)= 0.6465442349154d+01
      c(74)=-0.5832144628726d+01
      c(75)= 0.2684557558935d+01
      c(76)=-0.6830155078836d+00
      c(77)= 0.9186854104822d-01
      c(78)=-0.5120121563474d-02
c...2.50 @ x @ 2.75  interval no. 11  abs.error=1.5842438472191d-
      ifirst(11)=79
      ilast (11)=85
      c(79)=-0.2550246682056d+01
      c(80)= 0.6663085650406d+01
      c(81)=-0.6057008239737d+01
      c(82)= 0.2819123791215d+01
      c(83)=-0.7278168845029d+00
      c(84)= 0.9975271609922d-01
      c(85)=-0.5693960934877d-02
c...2.75 @ x @ 3.00  interval no. 12  abs.error=2.1390000881638d-
      ifirst(12)=86
      ilast (12)=92
      c(86)=-0.1437322524459d+01
      c(87)= 0.4244945032011d+01
      c(88)=-0.3866819653025d+01
      c(89)= 0.1760657886042d+01
      c(90)=-0.4399493394691d+00
      c(91)= 0.5797869518089d-01
      c(92)=-0.3166941925883d-02
c...3.00 @ x @ 3.25  interval no. 13  abs.error=1.2176926134089d-
      ifirst(13)=93
      ilast (13)=99
      c(93)= 0.4326982515539d-01
      c(94)= 0.1278541925538d+01
      c(95)=-0.1389457146276d+01
      c(96)= 0.6567745067205d+00
      c(97)=-0.1631568799154d+00
      c(98)= 0.2094831465123d-01
      c(99)=-0.1101922864715d-02
c...3.25 @ x @ 3.50  interval no. 14  abs.error=4.2543746303636d-
      ifirst(14) =100
      ilast (14) =106
      c(100)= 0.1115286202567d+01
      c(101)=-0.7098356509943d+00
      c(102)= 0.1477452918517d+00
      c(103)= 0.2274875774917d-01
      c(104)=-0.1601072077271d-01
      c(105)= 0.2728999281923d-02
      c(106)=-0.1616695274909d-03
c...3.50 @ x @ 3.75  interval no. 15  abs.error=2.3193003073629d-
      ifirst(15) =107
      ilast (15) =112
      c(107)= 0.1346851272697d+01
      c(108)=-0.1125834241251d+01
      c(109)= 0.4583884751028d+00
      c(110)=-0.1007033665315d+00
      c(111)= 0.1153266357142d-01
      c(112)=-0.5427038297057d-03
c...3.75 @ x @  4.00   interval no. 16   abs.error=2.6108004647085d-
      ifirst(16) =113
      ilast (16) =118
      c(113)= 0.1220883963156d+01
      c(114)=-0.9581085557663d+00
      c(115)= 0.3690371210149d+00
      c(116)=-0.7689811239270d-01
      c(117)= 0.8360797545174d-02
      c(118)=-0.3736140672117d-03
c...4.00 @ x @  4.25   interval no. 17   abs.error=1.7570833676928d-
      ifirst(17) =119
      ilast (17) =124
      c(119)= 0.1079863694049d+01
      c(120)=-0.7816343308043d+00
      c(121)= 0.2806819253411d+00
      c(122)=-0.5477512389825d-01
      c(123)= 0.5590565502644d-02
      c(124)=-0.2348302397877d-03
c...4.25 @ x @  4.50   interval no. 18   abs.error=1.0521805648978d-
      ifirst(18) =125
      ilast (18) =130
      c(125)= 0.9600857053756d+00
      c(126)=-0.6404875010594d+00
      c(127)= 0.2141383549717d+00
      c(128)=-0.3908625281056d-01
      c(129)= 0.3740753346938d-02
      c(130)=-0.1475725788623d-03
c...4.50 @ x @  4.75   interval no. 19   abs.error=6.2718719107124d-
      ifirst(19) =131
      ilast (19) =136
      c(131)= 0.8652746664820d+00
      c(132)=-0.5349723088368d+00
      c(133)= 0.1671594028443d+00
      c(134)=-0.2862620440719d-01
      c(135)= 0.2576075398247d-02
      c(136)=-0.9569143876433d-04
c...4.75 @ x @  5.00   interval no. 20   abs.error=3.9639402871217d-
      ifirst(20) =137
      ilast (20) =142
      c(137)= 0.7900679556534d+00
      c(138)=-0.4556957919175d+00
      c(139)= 0.1337277737706d+00
      c(140)=-0.2157594052042d-01
      c(141)= 0.1832563601783d-02
      c(142)=-0.6432286463678d-04
c...5.00 @ x @  5.25   interval no. 21   abs.error=2.4784618801732d-
      ifirst(21) =143
      ilast (21) =148
      c(143)= 0.7290302365904d+00
      c(144)=-0.3945841937680d+00
      c(145)= 0.1092502418252d+00
      c(146)=-0.1667318620030d-01
      c(147)= 0.1341496314853d-02
      c(148)=-0.4464581143111d-04
c...5.25 @ x @  5.50   interval no. 22   abs.error=1.6431300764452d-
      ifirst(22) =149
      ilast (22) =154
      c(149)= 0.6781821897825d+00
      c(150)=-0.3461064436915d+00
      c(151)= 0.9076079662805d-01
      c(152)=-0.1314681150288d-01
      c(153)= 0.1005173704471d-02
      c(154)=-0.3181374631822d-04
c...5.50 @ x @  5.75   interval no. 23   abs.error=1.1222134332911d-
      ifirst(23) =155
      ilast (23) =160
      c(155)= 0.6349563489690d+00
      c(156)=-0.3067739392218d+00
      c(157)= 0.7644326048783d-01
      c(158)=-0.1054063730153d-01
      c(159)= 0.7679506117711d-03
      c(160)=-0.2317563630640d-04
c...5.75 @ x @  6.00   interval no. 24   abs.error=7.8559381222476d-
      ifirst(24) =161
      ilast (24) =166
      c(161)= 0.5976075476333d+00
      c(162)=-0.2742701592908d+00
      c(163)= 0.6512719379650d-01
      c(164)=-0.8570613370193d-02
      c(165)= 0.5964515294181d-03
      c(166)=-0.1720313448459d-04
c...6.00 @ x @  6.25   interval no. 25   abs.error=4.6252335295094d-
      ifirst(25) =167
      ilast (25) =171
      c(167)= 0.4530821605957d+00
      c(168)=-0.1556937202255d+00
      c(169)= 0.2621052097467d-01
      c(170)=-0.2184050212463d-02
      c(171)= 0.7237572572194d-04
c...6.25 @ x @  6.50   interval no. 26   abs.error=3.5389913222161d-
      ifirst(26) =172
      ilast (26) =176
      c(172)= 0.4314398767018d+00
      c(173)=-0.1418342842259d+00
      c(174)= 0.2288200151943d-01
      c(175)=-0.1828741476402d-02
      c(176)= 0.5815166332468d-04
c...6.50 @ x @  6.75   interval no. 27   abs.error=2.7426949600340d-
      ifirst(27) =177
      ilast (27) =181
      c(177)= 0.4119346702408d+00
      c(178)=-0.1298244601295d+00
      c(179)= 0.2010878728234d-01
      c(180)=-0.1544113849832d-02
      c(181)= 0.4719618664240d-04
c...6.75 @ x @  7.00   interval no. 28   abs.error=2.1502799540940d-
      ifirst(28) =182
      ilast (28) =186
      c(182)= 0.3942461802274d+00
      c(183)=-0.1193370928285d+00
      c(184)= 0.1777693657917d-01
      c(185)=-0.1313662248549d-02
      c(186)= 0.3865502003464d-04
c...7.00 @ x @  7.25   interval no. 29   abs.error=1.7034818000639d-
      ifirst(29) =187
      ilast (29) =191
      c(187)= 0.3781177456714d+00
      c(188)=-0.1101165554215d+00
      c(189)= 0.1580007215659d-01
      c(190)=-0.1125279834923d-02
      c(191)= 0.3192276744812d-04
c...7.25 @ x @  7.50   interval no. 30   abs.error=1.3624212868990d-
      ifirst(30) =192
      ilast (30) =196
      c(192)= 0.3633408065001d+00
      c(193)=-0.1019602554247d+00
      c(194)= 0.1411174419289d-01
      c(195)=-0.9699476170226d-03
      c(196)= 0.2656330616446d-04
c...7.50 @ x @  7.75   interval no. 31   abs.error=1.0990319765369d-
      ifirst(31) =197
      ilast (31) =201
      c(197)= 0.3497439218126d+00
      c(198)=-0.9470569276499d-01
      c(199)= 0.1266017939918d-01
      c(200)=-0.8408550675085d-03
      c(201)= 0.2225784919574d-04
c...7.75 @ x @  8.00   interval no. 32   abs.error=8.9381835266522d-
      ifirst(32) =202
      ilast (32) =206
      c(202)= 0.3371845853069d+00
      c(203)=-0.8822105128463d-01
      c(204)= 0.1140456233287d-01
      c(205)=-0.7327946570967d-03
      c(206)= 0.1877023896668d-04
c...8.00 @ x @  8.25   interval no. 33   abs.error=7.3210326689832d-
      ifirst(33) =207
      ilast (33) =211
      c(207)= 0.3255431584494d+00
      c(208)=-0.8239832407718d-01
      c(209)= 0.1031237441518d-01
      c(210)=-0.6417393667562d-03
      c(211)= 0.1592339503986d-04
c...8.25 @ x @  8.50   interval no. 34   abs.error=6.0371707633067d-
      ifirst(34) =212
      ilast (34) =216
      c(212)= 0.3147180811939d+00
      c(213)=-0.7714810368023d-01
      c(214)= 0.9357439877093d-02
      c(215)=-0.5645414814808d-03
      c(216)= 0.1358301324217d-04
c...8.50 @ x @  8.75   interval no. 35   abs.error=5.0113246885530d-
      ifirst(35) =217
      ilast (35) =221
      c(217)= 0.3046231027165d+00
      c(218)=-0.7239608540250d-01
      c(219)= 0.8518560948886d-02
      c(220)=-0.4987218951555d-03
      c(221)= 0.1164632612927d-04
c...8.75 @ x @  9.00   interval no. 36   abs.error=4.1842085352073d-
      ifirst(36) =222
      ilast (36) =226
      c(222)= 0.2951839772317d+00
      c(223)=-0.6807982796391d-01
      c(224)= 0.7778392737741d-02
      c(225)=-0.4423078514719d-03
      c(226)= 0.1003385659715d-04
c...9.00 @ x @  9.25   interval no. 37   abs.error=3.5125236053091d-
      ifirst(37) =227
      ilast (37) =231
      c(227)= 0.2863364933588d+00
      c(228)=-0.6414655846973d-01
      c(229)= 0.7122648855400d-02
      c(230)=-0.3937177713169d-03
      c(231)= 0.8683627129358d-05
c...9.25 @ x @  9.50   interval no. 38   abs.error=2.9645175203541d-
      ifirst(38) =232
      ilast (38) =236
      c(232)= 0.2780246155388d+00
      c(233)=-0.6055132372497d-01
      c(234)= 0.6539470823348d-02
      c(235)=-0.3516734769846d-03
      c(236)= 0.7546893357357d-05
c...9.50 @ x @  9.75   interval no. 39   abs.error=2.5142110615661d-
      ifirst(39) =237
      ilast (39) =241
      c(237)= 0.2701996594546d+00
      c(238)=-0.5725581832360d-01
      c(239)= 0.6018987049956d-02
      c(240)=-0.3151372568482d-03
      c(241)= 0.6585092705791d-05
c...9.75 @ x @ 10.00   interval no. 40   abs.error=2.1422863483167d-
      ifirst(40) =242
      ilast (40) =246
      c(242)= 0.2628187673463d+00
      c(243)=-0.5422707501601d-01
      c(244)= 0.5552907092188d-02
      c(245)=-0.2832594125266d-03
      c(246)= 0.5767453330918d-05
      return
      end
**==dawf.f
      function dawf(y)
      implicit real*8 (a-h,o-z)
c
      real*8 c, h
      integer ifirst, ilast
      common/dawfcm/c(249),ifirst(40),ilast(40),h
c
c
c        -----  routine evalues the dawson function.            -----
c
      x = y
      sign = +1.0d+00
      if(x.ge.0.0d+00) go to 10
      x = -y
      sign = -1.0d+00
   10 if(x.ge.10.0d+00) go to 30
      xn = x/h
      nx = int(xn)
      nx = nx+1
      if = ifirst(nx)
      il = ilast(nx)
      t = c(il)
      kl = il-if
      do 20 k=1,kl
   20 t = c(il-k)+x*t
      dawf = t*sign
      return
   30 txt = 1.0d+00/(2.0d+00*x*x)
      tx = txt*x*sign
      dawf = tx*(1.0d+00+txt*(1.0d+00+txt*(3.0d+00+txt*(15.0d+00+
     * 105.0d+00*txt))))
      return
      end
**==dawt.f
      subroutine dawt
      implicit real*8 (a-h,o-z)
c        -----  routine allocates parameters specifying         -----
c        -----  the piecewise chebyshev polonomial fit to       -----
c        -----  the dawson function.                            -----
c
c
      real*8 c, h
      integer ifirst, ilast
      common/dawfcm/c(249),ifirst(40),ilast(40),h
c
c
      h=0.25d+00
c... 0.00 @ x @ .25 interval no. 1 abs.error = 1.6571632954765d-
      ifirst( 1) =  1
      ilast( 1) =  8
      c( 1) = -0.1587492209085d-15
      c( 2) =  0.1000000001679d+01
      c( 3) = -0.2200140634227d-06
      c( 4) = -0.6666582802511d+00
      c( 5) = -0.1406923341036d-03
      c( 6) =  0.2678610335986d+00
      c( 7) = -0.5192471402032d-02
      c( 8) = -0.6640836596489d-01
c...  .25 @ x @  .50  interval no. 2  abs.error = 3.1585400961376d-
      ifirst( 2) =  9
      ilast ( 2) = 16
      c( 9) = -0.1012320887742d-04
      c(10) =  0.1000229551438d+01
      c(11) = -0.2249318275452d-02
      c(12) = -0.6542497725108d+00
      c(13) = -0.4207802992979d-01
      c(14) =  0.3555122257011d+00
      c(15) = -0.1112488455006d+00
      c(16) = -0.8490860462189d-02
c... .50 @ x @ .75 interval no. 3 abs.error = 1.5242918038894d-
      ifirst( 3) = 17
      ilast ( 3) = 24
      c(17) = -0.4921356497286d-03
      c(18) =  0.1006613917835d+01
      c(19) = -0.3889997394948d-01
      c(20) = -0.5358795704274d+00
      c(21) = -0.2746372044452d+00
      c(22) =  0.6337005279825d+00
      c(23) = -0.2989429214171d+00
      c(24) =  0.4661336966923d-01
c...  .75 @ x @ 1.00  interval no. 4  abs.error = 9.9014130228170d-
      ifirst( 4) = 25
      ilast ( 4) = 32
      c(25) =  0.2895369459824d-03
      c(26) =  0.1001188430399d+01
      c(27) = -0.2450339520101d-01
      c(28) = -0.5518735978257d+00
      c(29) = -0.2745338921003d+00
      c(30) =  0.6506738753191d+00
      c(31) = -0.3141764870712d+00
      c(32) =  0.5101503644671d-01
c... 1.00 @ x @ 1.25  interval no. 5  abs.error = 1.8118839761883d-
      ifirst( 5) = 33
      ilast ( 5) = 40
      c(33) =  0.4002056337725d-01
      c(34) =  0.7312836053563d+00
      c(35) =  0.7635870948288d+00
      c(36) = -0.1834201228325d+01
      c(37) =  0.9813999234300d+00
      c(38) = -0.8984573983721d-01
      c(39) = -0.7076886296272d-01
      c(40) =  0.1660415104457d-01
c... 1.25 @ x @ 1.50  interval no. 6  abs.error = 8.9031004790741d-
      ifirst( 6) = 41
      ilast ( 6) = 48
      c(41) =  0.1852097215132d+00
      c(42) = -0.7878409599634d-01
      c(43) =  0.2705240106371d+01
      c(44) = -0.4425988758096d+01
      c(45) =  0.3062214704197d+01
      c(46) = -0.1094631867084d+01
      c(47) =  0.1994346954993d+00
      c(48) = -0.1461141450065d-01
c... 1.50 @ x @ 1.75  interval no. 7  abs.error = 2.2932766796657d-
      ifirst( 7) = 49
      ilast ( 7) = 56
      c(49) =  0.2566795947043d+00
      c(50) = -0.4340789653964d+00
      c(51) =  0.3460647185102d+01
      c(52) = -0.5316677310217d+01
      c(53) =  0.3691359889577d+01
      c(54) = -0.1360915547303d+01
      c(55) =  0.2619748477425d+00
      c(56) = -0.2090002809252d-01
c... 1.75 @ x @ 2.00  interval no. 8  abs.error = 5.7358562344234d-
      ifirst( 8) = 57
      ilast ( 8) = 64
      c(57) = -0.2953293234440d+00
      c(58) =  0.1741536362797d+01
      c(59) = -0.2179107446162d+00
      c(60) = -0.1857700863729d+01
      c(61) =  0.1737803495672d+01
      c(62) = -0.6982079451638d+00
      c(63) =  0.1369430933680d+00
      c(64) = -0.1077900614057d-01
c... 2.00 @ x @ 2.25  interval no. 9  abs.error = 7.9634077110313d-
      ifirst( 9) = 65
      ilast ( 9) = 71
      c(65) = -0.1683235354314d+01
      c(66) =  0.6579307910802d+01
      c(67) = -0.7452680221211d+01
      c(68) =  0.4159630869064d+01
      c(69) = -0.1268324481886d+01
      c(70) =  0.2038507675752d+00
      c(71) = -0.1359998434782d-01
c... 2.25 @ x @ 2.50  interval no. 10  abs.error = 4.9522164147220d-
      ifirst(10) = 72
      ilast (10) = 78
      c(72) = -0.1172884307439d+01
      c(73) =  0.5243170888632d+01
      c(74) = -0.5994671204809d+01
      c(75) =  0.3310821796507d+01
      c(76) = -0.9902715162510d+00
      c(77) =  0.1552556976676d+00
      c(78) = -0.1006003345052d-01
c... 2.50 @ x @ 2.75  interval no. 11  abs.error = 4.0683012514364d-
      ifirst(11) = 79
      ilast (11) = 85
      c(79) =  0.2453071419396d+00
      c(80) =  0.1843675794804d+01
      c(81) = -0.2597381166786d+01
      c(82) =  0.1499064046402d+01
      c(83) = -0.4464691670097d+00
      c(84) =  0.6815297792976d-01
      c(85) = -0.4243532816569d-02
c... 2.75 @ x @ 3.00  interval no. 12  abs.error = 1.7797319173951d-
      ifirst(12) = 86
      ilast (12) = 92
      c(86) =  0.1711553394682d+01
      c(87) = -0.1366851519128d+01
      c(88) =  0.3331226440003d+00
      c(89) =  0.7176989090082d-01
      c(90) = -0.5525697253082d-01
      c(91) =  0.1093749119900d-01
      c(92) = -0.7553007453680d-03
c... 3.00 @ x @ 3.25  interval no. 13  abs.error = 3.6681768733615d-
      ifirst(13) = 93
      ilast (13) = 99
      c(93) =  0.2517095374277d+01
      c(94) = -0.2990903746754d+01
      c(95) =  0.1697891140418d+01
      c(96) = -0.5401196813139d+00
      c(97) =  0.9911334144514d-01
      c(98) = -0.9840574037905d-02
      c(99) =  0.4103935013215d-03
c... 3.25 @ x @ 3.50  interval no. 14  abs.error = 1.1128875598843d-
      ifirst(14) =100
      ilast (14) =106
      c(100) =  0.2590851813134d+01
      c(101) = -0.3134054423500d+01
      c(102) =  0.1813437074513d+01
      c(103) = -0.5897752524597d+00
      c(104) =  0.1110979078173d+00
      c(105) = -0.1138103745567d-01
      c(106) =  0.4927876094977d-03
c... 3.50 @ x @ 3.75  interval no. 15  abs.error = 1.4752643551219d-
      ifirst(15) =107
      ilast (15) =113
      c(107) =  0.2263419843285d+01
      c(108) = -0.2574296380066d+01
      c(109) =  0.1414605080120d+01
      c(110) = -0.4381747671623d+00
      c(111) =  0.7867466333846d-01
      c(112) = -0.7681612002974d-02
      c(113) =  0.3168638795614d-03
c... 3.75 @ x @ 4.00  interval no. 16  abs.error = 9.2281737806843d-
      ifirst(16) =114
      ilast (16) =120
      c(114) =  0.1858711024930d+01
      c(115) = -0.1925962656395d+01
      c(116) =  0.9817358192129d+00
      c(117) = -0.2839951967696d+00
      c(118) =  0.4777651946157d-01
      c(119) = -0.4378302333256d-02
      c(120) =  0.1696776598692d-03
c... 4.00 @ x @ 4.25  interval no. 17  abs.error = 2.0758506025231d-
      ifirst(17) =121
      ilast (17) =126
      c(121) =  0.1097210930469d+01
      c(122) = -0.8018697789038d+00
      c(123) =  0.2901340972873d+00
      c(124) = -0.5698500686831d-01
      c(125) =  0.5849148664856d-02
      c(126) = -0.2469443250448d-03
c... 4.25 @ x @ 4.50  interval no. 18  abs.error = 1.1124878795954d-
      ifirst(18) =127
      ilast (18) =132
      c(127) =  0.9640388573940d+00
      c(128) = -0.6448585002328d+00
      c(129) =  0.2160731515099d+00
      c(130) = -0.3951479965599d-01
      c(131) =  0.3788248790079d-02
      c(132) = -0.1496796030551d-03
c... 4.50 @ x @ 4.75  interval no. 19  abs.error = 6.3691274476696d-
      ifirst(19) =133
      ilast (19) =138
      c(133) =  0.8660393462895d+00
      c(134) = -0.5357754986697d+00
      c(135) =  0.1674970705868d+00
      c(136) = -0.2869722599326d-01
      c(137) =  0.2583548621624d-02
      c(138) = -0.9600615594536d-04
c... 4.75 @ x @ 5.00  interval no. 20  abs.error = 3.8777869804107d-
      ifirst(20) =139
      ilast (20) =144
      c(139) =  0.7902176622367d+00
      c(140) = -0.4558460385033d+00
      c(141) =  0.1337881166153d+00
      c(142) = -0.2158806335901d-01
      c(143) =  0.1833781835739d-02
      c(144) = -0.6437185220420d-04
c... 5.00 @ x @ 5.25  interval no. 21  abs.error = 2.4793500585929d-
      ifirst(21) =145
      ilast (21) =150
      c(145) =  0.7290485015677d+00
      c(146) = -0.3946016171968d+00
      c(147) =  0.1092568925649d+00
      c(148) = -0.1667445598850d-01
      c(149) =  0.1341617573053d-02
      c(150) = -0.4465044476092d-04
c... 5.25 @ x @ 5.50  interval no. 22  abs.error = 1.6431300764452d-
      ifirst(22) =151
      ilast (22) =156
      c(151) =  0.6781842728186d+00
      c(152) = -0.3461083410440d+00
      c(153) =  0.9076148811696d-01
      c(154) = -0.1314693754430d-01
      c(155) =  0.1005185194663d-02
      c(156) = -0.3181416541338d-04
c... 5.50 @ x @ 5.75  interval no. 23  abs.error = 1.1226575225010d-
      ifirst(23) =157
      ilast (23) =162
      c(157) =  0.6349565075405d+00
      c(158) = -0.3067740750644d+00
      c(159) =  0.7644330701296d-01
      c(160) = -0.1054064526443d-01
      c(161) =  0.7679512928007d-03
      c(162) = -0.2317565958947d-04
c... 5.75 @ x @ 6.00  interval no. 24  abs.error = 7.8470563380506d-
      ifirst(24) =163
      ilast (24) =168
      c(163) =  0.5976070486768d+00
      c(164) = -0.2742697363056d+00
      c(165) =  0.6512705037617d-01
      c(166) = -0.8570589057672d-02
      c(167) =  0.5964494688669d-03
      c(168) = -0.1720306463540d-04
c... 6.00 @ x @ 6.25  interval no. 25  abs.error = 5.6132876125048d-
      ifirst(25) =169
      ilast (25) =174
      c(169) =  0.5649180454987d+00
      c(170) = -0.2470090261752d+00
      c(171) =  0.5603266263223d-01
      c(172) = -0.7053466938578d-02
      c(173) =  0.4698947566794d-03
      c(174) = -0.1297991257161d-04
c... 6.25 @ x @ 6.50  interval no. 26  abs.error = 3.5389913222161d-
      ifirst(26) =175
      ilast (26) =179
      c(175) =  0.4314398767018d+00
      c(176) = -0.1418342842259d+00
      c(177) =  0.2288200151943d-01
      c(178) = -0.1828741476402d-02
      c(179) =  0.5815166332468d-04
c... 6.50 @ x @ 6.75  interval no. 27  abs.error = 2.7426949600340d-
      ifirst(27) =180
      ilast (27) =184
      c(180) =  0.4119346702408d+00
      c(181) = -0.1298244601295d+00
      c(182) =  0.2010878728234d-01
      c(183) = -0.1544113849832d-02
      c(184) =  0.4719618664240d-04
c... 6.75 @ x @ 7.00  interval no. 28  abs.error = 2.1502799540940d-
      ifirst(28) =185
      ilast (28) =189
      c(185) =  0.3942461802274d+00
      c(186) = -0.1193370928285d+00
      c(187) =  0.1777693657917d-01
      c(188) = -0.1313662248549d-02
      c(189) =  0.3865502003464d-04
c... 7.00 @ x @ 7.25  interval no. 29  abs.error = 1.7034818000639d-
      ifirst(29) =190
      ilast (29) =194
      c(190) =  0.3781177456714d+00
      c(191) = -0.1101165554215d+00
      c(192) =  0.1580007215659d-01
      c(193) = -0.1125279834923d-02
      c(194) =  0.3192276744812d-04
c... 7.25 @ x @ 7.50  interval no. 30  abs.error = 1.3624212868990d-
      ifirst(30) =195
      ilast (30) =199
      c(195) =  0.3633408065001d+00
      c(196) = -0.1019602554247d+00
      c(197) =  0.1411174419289d-01
      c(198) = -0.9699476170226d-03
      c(199) =  0.2656330616446d-04
c... 7.50 @ x @ 7.75  interval no. 31  abs.error = 1.0990319765369d-
      ifirst(31) =200
      ilast (31) =204
      c(200) =  0.3497439218126d+00
      c(201) = -0.9470569276499d-01
      c(202) =  0.1266017939918d-01
      c(203) = -0.8408550675085d-03
      c(204) =  0.2225784919574d-04
c... 7.75 @ x @ 8.00 interval no. 32 abs.error = 8.9381835266522d-
      ifirst(32) =205
      ilast (32) =209
      c(205) =  0.3371845853069d+00
      c(206) = -0.8822105128463d-01
      c(207) =  0.1140456233287d-01
      c(208) = -0.7327946570967d-03
      c(209) =  0.1877023896668d-04
c... 8.00 @ x @ 8.25  interval no. 33  abs.error = 7.3210326689832d-
      ifirst(33) =210
      ilast (33) =214
      c(210) =  0.3255431584494d+00
      c(211) = -0.8239832407718d-01
      c(212) =  0.1031237441518d-01
      c(213) = -0.6417393667562d-03
      c(214) =  0.1592339503986d-04
c... 8.25 @ x @ 8.50 interval no. 34 abs.error = 6.0371707633067d-
      ifirst(34) =215
      ilast (34) =219
      c(215) =  0.3147180811939d+00
      c(216) = -0.7714810368023d-01
      c(217) =  0.9357439877093d-02
      c(218) = -0.5645414814808d-03
      c(219) =  0.1358301324217d-04
c... 8.50 @ x @ 8.75  interval no. 35  abs.error = 5.0113246885530d-
      ifirst(35) =220
      ilast (35) =224
      c(220) =  0.3046231027165d+00
      c(221) = -0.7239608540250d-01
      c(222) =  0.8518560948886d-02
      c(223) = -0.4987218951555d-03
      c(224) =  0.1164632612927d-04
c... 8.75 @ x @ 9.00 interval no. 36 abs.error = 4.1842085352073d-
      ifirst(36) =225
      ilast (36) =229
      c(225) =  0.2951839772317d+00
      c(226) = -0.6807982796391d-01
      c(227) =  0.7778392737741d-02
      c(228) = -0.4423078514719d-03
      c(229) =  0.1003385659715d-04
c... 9.00 @ x @ 9.25 interval no. 37 abs.error = 3.5125236053091d-
      ifirst(37) =230
      ilast (37) =234
      c(230) =  0.2863364933588d+00
      c(231) = -0.6414655846973d-01
      c(232) =  0.7122648855400d-02
      c(233) = -0.3937177713169d-03
      c(234) =  0.8683627129358d-05
c... 9.25 @ x @ 9.50 interval no. 38 abs.error = 2.9645175203541d-
      ifirst(38) =235
      ilast (38) =239
      c(235) =  0.2780246155388d+00
      c(236) = -0.6055132372497d-01
      c(237) =  0.6539470823348d-02
      c(238) = -0.3516734769846d-03
      c(239) =  0.7546893357357d-05
c... 9.50 @ x @ 9.75  interval no. 39  abs.error = 2.5142110615661d-
      ifirst(39) =240
      ilast (39) =244
      c(240) =  0.2701996594546d+00
      c(241) = -0.5725581832360d-01
      c(242) =  0.6018987049956d-02
      c(243) = -0.3151372568482d-03
      c(244) =  0.6585092705791d-05
c... 9.75 @ x @ 10.00 interval no. 40 abs.error = 2.1422863483167d-
      ifirst(40) =245
      ilast (40) =249
      c(245) =  0.2628187673463d+00
      c(246) = -0.5422707501601d-01
      c(247) =  0.5552907092188d-02
      c(248) = -0.2832594125266d-03
      c(249) =  0.5767453330918d-05
      return
      end
**==dco.f
      function dco(l,m,kx,ky,kz,lp,mp,fpqr,zlm,
     *   lmf,lmx,lmy,lmz)
      implicit real*8 (a-h,o-z)
      dimension fpqr(25,25,25),zlm(*),lmf(*),lmx(*),lmy(*),lmz(*)
c
c        -----  routine evaluates angular momentum coupling     -----
c        -----  coefficients.                                   -----
c        for a given set of l and m, this routine calculates the type 1
c        (l,m=0) or the type 2 angular integrals. (ie the second line of
c        eq 28 or 29 in md's paper) you still must do the sums.
c
      id = l*(l+1)-m+1
      imn = lmf(id)
      imx = lmf(id+1)-1
      jd = lp*(lp+1)-mp+1
      jmn = lmf(jd)
      jmx = lmf(jd+1)-1
      sumi = 0.0d+00
      do i=imn,imx
        sumj = 0.0d+00
        do j=jmn,jmx
          jx = lmx(i)+kx+lmx(j)+1
          jy = lmy(i)+ky+lmy(j)+1
          jz = lmz(i)+kz+lmz(j)+1
          sumj = sumj+zlm(j)*fpqr(jx,jy,jz)
        enddo
      sumi = sumi+zlm(i)*sumj
      enddo
      dco = sumi
      return
      end
**==eccod1.f
      subroutine eccod1(dcoef1,jfst1,lbecp1,fpqr,zlm,lmf,lmx,lmy,lmz,
     *                  numder)
      implicit real*8 (a-h,o-z)
      dimension dcoef1(*),lbecp1(9,*),fpqr(25,25,25),lmf(*),jfst1(*)
      dimension lmx(*),lmy(*),lmz(*),zlm(*)
c
c        -----  routine finds formula code for the regular      -----
c        -----  three-center one electron integrals.            -----
c        regular meaning type 1 angular integrals
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer ncoef1,ncoef2,j1len,j2len,llim,nlim,ntlim,j4len
      common /ecpdim/ ncoef1,ncoef2,j1len,j2len,llim,nlim,ntlim,j4len
c
c
      integer nfst,nx,ny,nz
      common /gbase/ nfst(8),nx(84),ny(84),nz(84)
c
      dimension binco(28)
      parameter(sqrfpi=3.5449077018110d+00, tol=1.0d-10)
c binco -- indexed combinations (b!/(a!(b-a)!)) indexed by b and a
c          as:  (b,a): (0,0),(1,0),(1,1),(2,0),(2,1),(2,2)...
c aka pascal's triangle, entries through i functions (g hessians)
      data binco/          1.d+00,
     1                  1.d+00, 1.d+00,
     2               1.d+00, 2.d+00, 1.d+00,
     3            1.d+00, 3.d+00, 3.d+00, 1.d+00,
     4         1.d+00, 4.d+00, 6.d+00, 4.d+00, 1.d+00,
     5      1.d+00, 5.d+00,10.d+00,10.d+00, 5.d+00,1.d+00,
     6   1.d+00,6.d+00,15.d+00,20.d+00,15.d+00,6.d+00,1.d+00/
c
      jndx = 0
      jjlst=0
c  loop over possible i shells
      do 190 nn1=1,nlim
        nf1 = nfst(nn1)
        nl1 = nfst(nn1+1)-1
c
        if (numder.eq.0.or.nn1.lt.nlim) then
          nn2lim=nn1
        else
          nn2lim=nn1-numder
        end if
c loop over the i shell angualr momentum functions
        do 180 n1=nf1,nl1
          n1t = (n1*(n1-1))/2
          mx1 = (nx(n1)*(nx(n1)+1))/2+1
          my1 = (ny(n1)*(ny(n1)+1))/2+1
          mz1 = (nz(n1)*(nz(n1)+1))/2+1
c  loop over the possible j shells (note j shell <= i shell)
          do 170 nn2=1,nn2lim
            nf2 = nfst(nn2)
            nl2 = min(nfst(nn2+1)-1,n1)
c  loop over j shell angular momentum
            do 160 n2=nf2,nl2
              indx = n2+n1t
              mx2 = (nx(n2)*(nx(n2)+1))/2+1
              my2 = (ny(n2)*(ny(n2)+1))/2+1
              mz2 = (nz(n2)*(nz(n2)+1))/2+1
              llmax = nn1+nn2-2
c  loop over lamda, mu, kx...
              do 130 lamda=0,llmax
                do 120 mu=-lamda,lamda
                  do 110 kx=0,nx(n1)
                    do 100 ky=0,ny(n1)
                      do 90 kz=0,nz(n1)
                        do 80 kxp=0,nx(n2)
                          do 70 kyp=0,ny(n2)
                            do 60 kzp=0,nz(n2)
                              kappa = kx+ky+kz+kxp+kyp+kzp
c note no nonzero integrals for lamda > kappa
c if kappa + lamda is even go ahead a calculate the angular integral
                              if ((mod((lamda+kappa),2).ne.1).and.
     *                          (lamda.le.kappa)) then
                                dc=dco(0,0,kx+kxp,ky+kyp,kz+kzp,
     *                             lamda,mu,fpqr,zlm,lmf,lmx,lmy,lmz)
c if the integral is non-zero then store it away for future use
                                if(abs(dc).gt.tol) then
                                  jndx = jndx+1
                                  lbecp1(1,jndx)=kappa
                                  lbecp1(2,jndx)=lamda
                                  lbecp1(3,jndx)=mu
                                  lbecp1(4,jndx)=kx
                                  lbecp1(5,jndx)=ky
                                  lbecp1(6,jndx)=kz
                                  lbecp1(7,jndx)=kxp
                                  lbecp1(8,jndx)=kyp
                                  lbecp1(9,jndx)=kzp
c the integral is multiplied by the combination factors and by sqrt(4pi)
c since there is an extra 1/sqrt(4pi) in the integral
                                  dcoef1(jndx) = dc*binco(mx1+kx)*
     *                              binco(my1+ky)*binco(mz1+kz)*
     *                              binco(mx2+kxp)*binco(my2+kyp)*
     *                              binco(mz2+kzp)*sqrfpi
                                end if
                              end if
   60                       continue
   70                     continue
   80                   continue
   90                 continue
  100               continue
  110             continue
  120           continue
  130         continue
              jfst1(indx) = jjlst+1
              jfst1(indx+1) = jndx+1
              jjlst = jndx
  160       continue
  170     continue
  180   continue
  190 continue
c
c      check to be sure that the array has stayed within bounds
c
      if (jndx.gt.ncoef1) then
         write(iwr,1000) jndx, ncoef1
         call caserr('eccod1 out of bounds')
      end if
      ncoef1 = jndx
      return
 1000 format(/'****eccod1 out of bounds, used ',i7,' allowed ',i7)
      end
**==eccod2.f
      subroutine eccod2(dcoef2,jfst2,lbecp2,fpqr,zlm,lmf,lmx,lmy,lmz)
      implicit real*8 (a-h,o-z)
c
      dimension dcoef2(*),lbecp2(6,*),fpqr(25,25,25),lmf(*)
      dimension lmx(*),lmy(*),lmz(*),zlm(*), jfst2(*)
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer ncoef1,ncoef2,j1len,j2len,llim,nlim,ntlim,j4len
      common /ecpdim/ ncoef1,ncoef2,j1len,j2len,llim,nlim,ntlim,j4len
c
c
      integer nfst,nx,ny,nz
      common /gbase/ nfst(8),nx(84),ny(84),nz(84)
c
      dimension binco(28)
      parameter(tol=1.0d-10)
c binco -- indexed combinations (b!/(a!(b-a)!)) indexed by b and a
c          as:  (b,a): (0,0),(1,0),(1,1),(2,0),(2,1),(2,2)...
c aka pascal's triangle, entries through i functions (g hessians)
      data binco/          1.d+00,
     1                  1.d+00, 1.d+00,
     2              1.d+00, 2.d+00, 1.d+00,
     3           1.d+00, 3.d+00, 3.d+00, 1.d+00,
     4        1.d+00, 4.d+00, 6.d+00, 4.d+00, 1.d+00,
     5     1.d+00, 5.d+00,10.d+00,10.d+00, 5.d+00,1.d+00,
     6  1.d+00,6.d+00,15.d+00,20.d+00,15.d+00,6.d+00,1.d+00/
c
c        -----  routine finds formula code for one-electron     -----
c        -----  three-center integrals involving projection     -----
c        -----  operators.                                      -----
c        the type 2 angular integral table generator
c
      jndx = 0
      jjlst=0
      lllim=llim-1
      do 170 l=0,lllim
        do 160 m=l,-l,-1
          lmindx = (l*(l+1)-m)*ntlim
          do 150 nn=1,nlim
            nf = nfst(nn)
            nl = nfst(nn+1)-1
            l2mx = l+nn - 1
            do 140 n=nf,nl
              indx = lmindx+n
              mx = (nx(n)*(nx(n)+1))/2 + 1
              my = (ny(n)*(ny(n)+1))/2 + 1
              mz = (nz(n)*(nz(n)+1))/2 + 1
              do 130 lamda=0,l2mx
                do 120 mu=-lamda,lamda
                  do 110 kx=0,nx(n)
                    do 100 ky=0,ny(n)
                      do 90 kz=0,nz(n)
                        isigma = kx+ky+kz
                        dc = dco(l,m,kx,ky,kz,lamda,mu,
     *                       fpqr,zlm,lmf,lmx,lmy,lmz)
                        if(abs(dc).ge.tol) then
                          jndx = jndx+1
                          lbecp2(1,jndx) = lamda
                          lbecp2(2,jndx) = isigma
                          lbecp2(3,jndx) = mu
                          lbecp2(4,jndx) = kx
                          lbecp2(5,jndx) = ky
                          lbecp2(6,jndx) = kz
                          dcoef2(jndx) = dc*binco(mx+kx)*
     *                       binco(my+ky)*binco(mz+kz)
                        end if
   90                 continue
  100               continue
  110             continue
  120           continue
  130         continue
              jfst2(indx) = jjlst+1
              jfst2(indx+1) = jndx+1
              jjlst = jndx
  140       continue
  150     continue
  160   continue
  170 continue
c
c      check to be sure that the array has stayed within bounds
c
      if (jndx.gt.ncoef2) then
         write(iwr,1000) jndx, ncoef2
         call caserr('eccod2 out of bounds')
      end if
      ncoef2=jndx
      return
 1000 format(/'****eccod2 out of bounds, used ',i6,' allowed ',i5)
      end
**==eccod3.f
      subroutine eccod3(fpqr,dcoef4,zlm,lmf,lmx,lmy,lmz)
      implicit real*8 (a-h,o-z)
      dimension fpqr(25,25,25),zlm(*),lmf(*),lmx(*),lmy(*),lmz(*)
      dimension dcoef4(*)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer ncoef1,ncoef2,j1len,j2len,llim,nlim,ntlim,j4len
      common /ecpdim/ ncoef1,ncoef2,j1len,j2len,llim,nlim,ntlim,j4len
c
c
      integer nfst,nx,ny,nz
      common /gbase/ nfst(8),nx(84),ny(84),nz(84)
c
      data sqrfpi/3.5449077018110d+00/
c
      iiidx=0
      do 70 l=0,llim-1
        do 60 m=l,-l,-1
          do 50 nn=1,nlim
            nf = nfst(nn)
            nl = nfst(nn+1)-1
            do 40 n = nf,nl
              indx = (l*(l+1)-m)*ntlim+n
              iiidx = max(indx, iiidx)
              dcoef4(indx) = dco(0,0,nx(n),ny(n),nz(n),l,m,
     *             fpqr,zlm,lmf,lmx,lmy,lmz)*sqrfpi
   40       continue
   50     continue
   60   continue
   70 continue
      if (iiidx.gt.j4len) then
        write(iwr,*) 'dcoef4 out of bounds in eccod3 ',iiidx
        call caserr('eccod3 out of bounds')
      end if
      return
      end
**==eccodr.f
      subroutine eccodr(dcoef1,jfst1,lbecp1,dcoef2,jfst2,
     *                  lbecp2,fpqr,zlm,lmf,lmx, lmy, lmz, numder)
c
      implicit real*8 (a-h,o-z)
c
      dimension dcoef1(*),lbecp1(9,*),fpqr(25,25,25),zlm(581)
      dimension dcoef2(*),lbecp2(6,*),lmf(122),lmx(581),lmy(581)
      dimension lmz(581), jfst1(*), jfst2(*)
c
c
      integer ncoef1,ncoef2,j1len,j2len,llim,nlim,ntlim,j4len
      common /ecpdim/ ncoef1,ncoef2,j1len,j2len,llim,nlim,ntlim,j4len
c
c
c        -----  routine pre-calculates integral formulas for    -----
c        -----  for calculation of the matrix elements of the   -----
c        -----  effective core potential.                       -----
c        -----  for a cartesian gaussian basis.                 -----
c
c        -----  specify the maximal type of basis function      -----
c        -----  and maximal type of projection oper.            -----
c
c        -----  set up basis angular integrals.                 -----
c
      call ecpini(lmf,lmx,lmy,lmz)
      call ztab(zlm)
      call ftab(fpqr, nlim-1)
c
c        -----  set up the code for the integrals of the local  -----
c        -----  operator potential.                             -----
c
      call eccod1(dcoef1,jfst1,lbecp1,fpqr,zlm,lmf,lmx,lmy,lmz,numder)
c
c        -----  set up the code for the integrals of angular    -----
c        -----  momentum projection operator.                   -----
c
      call eccod2(dcoef2,jfst2,lbecp2,fpqr,zlm,lmf,lmx,lmy,lmz)
c
      return
      end
**==ecpa11.f
      subroutine ecpa11(g,coefi,coefj,fpqr,npnp)
c
      implicit real*8 (a-h,o-z)
c
c  ecp type 1 integrals for the special case of all
c  three functions on the same center <aaa> (ie. icab = 1)
c  g will contain the gradient vector upon exit (not symmetrized)
c  npnp = n+n' gives the maximum angular momentum (kappa) needed
c  note: l shells can not be done with this routine due to the coeff
c  of the prims multiplied in the radial integral section
c
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
      dimension g(*), coefi(*), coefj(*), fpqr(25,25,25)
c
c
      real*8 zetc,cax,cay,caz,ca,xca,yca,zca
      real*8 zetb,bax,bay,baz,ba,xba,yba,zba
      real*8 phase,dax,day,daz,da,xda,yda,zda,xint
      integer kcntr
      common /ecp1/ zetc,cax,cay,caz,ca,xca,yca,zca,
     +              zetb,bax,bay,baz,ba,xba,yba,zba,
     +              phase,dax,day,daz,da,xda,yda,zda,xint,
     +              kcntr
c
      common /ecp2  / clp(400),zlp(400),nlp(400),kfirst(maxat,6),
     *                klast(maxat,6),lmax(maxat),lpskip(maxat),
     *                izcore(maxat)
c
      integer nfst,nx,ny,nz
      common /gbase/ nfst(8),nx(84),ny(84),nz(84)
c
c   iamin/max are the angular function min and max for i
c   ipmin/max are the primitive function min and max for i
c   kf1/kl1 are the first and last primitives forming the lmax ecp
c      and thus the type 1 integral
c
      logical iandj, norm, normi, normj
      common /ecpidx/ q2,iamin,iamax,jamin,jamax,ipmin,ipmax,jpmin,
     *                jpmax,kf1,kl1,llmx,npc,npb,iandj,norm,normi,normj
c
      parameter(zero=0.0d+00)
      parameter(sqrt3=1.73205080756888d+00,sqrt5=2.23606797749979d+00)
      parameter (sqrt7=2.64575131106459d+00)
c
c   nbc = nxc + nyc + nzc + nbx + nby + nbz
      nbc = npnp - 1
c
c   the total n for the two shells must be even or the integral is
c   identically zero
c
      if (mod(nbc,2).eq.1) return
c
c   loop over the primitives summing the radial integrals in fpsum
c
      fpsum = zero
      igii=1
      do 400 ig = ipmin,ipmax
        dum = coefi(igii)
        igii=igii+1
        zetc = ex(ig)
        jgjj=1
        do 300 jg = jpmin,jpmax
          dum2 = dum*coefj(jgjj)
          jgjj=jgjj+1
          zetb = ex(jg)
          zetcb = zetc + zetb
          fptemp = zero
          do 200 k=kf1,kl1
            zeta = zetcb+zlp(k)
            fptemp = fptemp + fa(nbc+nlp(k), zeta)*clp(k)
 200      continue
          fpsum = fpsum + fptemp*dum2
 300    continue
 400  continue
c
c   now retrieve the angular integral and multiply it by the
c   radial integral sum and store in g
c
      nn = 1
      mmax=jamax
      do 1000 i=iamin,iamax
        dumi = q2*fpsum
        if (normi) then
          if (i.ge.8.and.i.le.10) then
            dumi = q2*fpsum*sqrt3
          else if (i.ge.14.and.i.le.20) then
            dumi = q2*fpsum*sqrt5
            if (i.eq.20) dumi=dumi*sqrt3
          else if (i.ge.24) then
            dumi = q2*fpsum*sqrt7
            if (i.ge.30) then
              dumi=dumi*sqrt5/sqrt3
              if (i.ge.33) dumi=dumi*sqrt3
            end if
          end if
        end if
        if (iandj) mmax=i
        do 900 j=jamin,mmax
          dumj = dumi
          if (normj) then
            if (j.ge.8.and.j.le.10) then
              dumj = dumi*sqrt3
            else if (j.ge.14.and.j.le.20) then
              dumj = dumi*sqrt5
              if (j.eq.20) dumj=dumj*sqrt3
            else if (j.ge.24) then
              dumj = dumi*sqrt7
              if (j.ge.30) then
                dumj=dumj*sqrt5/sqrt3
                if (j.ge.33) dumj=dumj*sqrt3
              end if
            end if
          end if
          nxc = nx(i)
          nyc = ny(i)
          nzc = nz(i)
          nxb = nx(j)
          nyb = ny(j)
          nzb = nz(j)
          dpqrt = fpqr(nxc+nxb+1,nyc+nyb+1,nzc+nzb+1)
          g(nn) = g(nn) + dpqrt*dumj
          nn = nn + 1
  900   continue
 1000 continue
      return
      end
**==ecpa12.f
      subroutine ecpa12(fp,jfst1,lbecp1,dcoef1,g,icab,npnp,
     *                  zlm,lmf,lmx,lmy,lmz)
      implicit real*8 (a-h,o-z)
c
c  ecp type 1 angular integrals
c  note the actual integral has already been done, we just
c  need to compute the appropriate factors and multiply by the
c  stored integral.
c  fp will contain all the required radial integrals upon entry
c  g will contain the gradient vector upon exit (not symmetrized)
c  lbecp1 contains the packed indeces for the angular integrals
c  dcoef1 contains the packed angular integral values
c  npnp = n+n' gives the maximum angular momentum (kappa) needed
c
      dimension dcoef1(*), lbecp1(9,*),jfst1(*)
      dimension fp(*), ckl(11,11), g(*)
      dimension zfnlm(121),zlm(*),lmf(*),lmx(*),lmy(*),lmz(*)
c
c
      integer nfst,nx,ny,nz
      common /gbase/ nfst(8),nx(84),ny(84),nz(84)
c
c
      real*8 acx,acy,acz,abx,aby,abz
      common /ecp3/ acx(7),acy(7),acz(7),abx(7),aby(7),abz(7)
c
c   iamin/max are the angular function min and max for i
c   ipmin/max are the primitive function min and max for i
c   kf1/kl1 are the first and last primitives forming the lmax ecp
c   and thus the type 1 integral
c
      logical iandj, norm, normi, normj
      common /ecpidx/ q2,iamin,iamax,jamin,jamax,ipmin,ipmax,jpmin,
     *                jpmax,kf1,kl1,llmx,npc,npb,iandj,norm,normi,normj
c
      parameter(zero=0.0d+00, one=1.0d+00)
      parameter(sqrt3=1.73205080756888d+00,sqrt5=2.23606797749979d+00)
      parameter (sqrt7=2.64575131106459d+00)
      parameter(tol=1.0d-10, fpi=12.566370614359d+00)
c
      call zfn(zfnlm, npnp-1,zlm,lmf,lmx,lmy,lmz)
c   nn indexes the gradient vector array
      nn=1
      mmax=jamax
      do 1000 i=iamin,iamax
        dumi = one
        if (normi) then
          if (i.ge.8.and.i.le.10) then
            dumi = sqrt3
          else if (i.ge.14.and.i.le.20) then
            dumi = sqrt5
            if (i.eq.20) dumi=dumi*sqrt3
          else if (i.ge.24) then
            dumi = sqrt7
            if (i.ge.30) then
              dumi=dumi*sqrt5/sqrt3
              if (i.ge.33) dumi=dumi*sqrt3
            end if
          end if
        end if
        if (iandj) mmax=i
        do 900 j=jamin,mmax
          dumj = q2
          if (normj) then
            if (j.ge.8.and.j.le.10) then
              dumj = q2*sqrt3
            else if (j.ge.14.and.j.le.20) then
              dumj = q2*sqrt5
              if (j.eq.20) dumj=dumj*sqrt3
            else if (j.ge.24) then
              dumj = q2*sqrt7
              if (j.ge.30) then
                dumj=dumj*sqrt5/sqrt3
                if (j.ge.33) dumj=dumj*sqrt3
              end if
            end if
          end if
          if (i.ge.j) then
            indx = j+(i*(i-1))/2
          else
            indx = i+(j*(j-1))/2
          end if
          jf = jfst1(indx)
          jl = jfst1(indx+1)-1
          nxc = nx(i)
          nyc = ny(i)
          nzc = nz(i)
          nxb = nx(j)
          nyb = ny(j)
          nzb = nz(j)
c  zero out sum which is used to collect each total integral
          sum = zero
c  zero out ckl which stores the angular integrals to be used later
          do 100 k=1,npnp
            do 100 l=1,k
  100     ckl(k,l) = zero
c
c    now calculate the angular integral which does not depend on primitives
c
          if (icab.eq.2) then
            nxyz = nxc+nyc+nzc
            do 200 jjj=jf,jl
              kappa = lbecp1(1,jjj)
              lamda = lbecp1(2,jjj)
              mu = lbecp1(3,jjj)
              if(i.ge.j) then
                kx = lbecp1(4,jjj)
                ky = lbecp1(5,jjj)
                kz = lbecp1(6,jjj)
                kxp = lbecp1(7,jjj)
                kyp = lbecp1(8,jjj)
                kzp = lbecp1(9,jjj)
              else
                kx = lbecp1(7,jjj)
                ky = lbecp1(8,jjj)
                kz = lbecp1(9,jjj)
                kxp = lbecp1(4,jjj)
                kyp = lbecp1(5,jjj)
                kzp = lbecp1(6,jjj)
              end if
              if((kx+ky+kz).eq.nxyz) then
                ckl(kappa+1,lamda+1) =ckl(kappa+1,lamda+1)
     *             + dcoef1(jjj)*zfnlm(lamda*(lamda+1)-mu+1)*
     *             abx(nxb-kxp+1)*aby(nyb-kyp+1)*abz(nzb-kzp+1)
              end if
  200       continue
c
          else
            nxyz = nxb+nyb+nzb
            do 300 jjj=jf,jl
              kappa = lbecp1(1,jjj)
              lamda = lbecp1(2,jjj)
              mu = lbecp1(3,jjj)
              if(i.ge.j) then
                kx = lbecp1(4,jjj)
                ky = lbecp1(5,jjj)
                kz = lbecp1(6,jjj)
                kxp = lbecp1(7,jjj)
                kyp = lbecp1(8,jjj)
                kzp = lbecp1(9,jjj)
              else
                kx = lbecp1(7,jjj)
                ky = lbecp1(8,jjj)
                kz = lbecp1(9,jjj)
                kxp = lbecp1(4,jjj)
                kyp = lbecp1(5,jjj)
                kzp = lbecp1(6,jjj)
              end if
              if((kxp+kyp+kzp).eq.nxyz) then
                ckl(kappa+1,lamda+1) =ckl(kappa+1,lamda+1)
     *                 + dcoef1(jjj)*zfnlm(lamda*(lamda+1)-mu+1)*
     *                 acx(nxc-kx+1)*acy(nyc-ky+1)*acz(nzc-kz+1)
              end if
  300       continue
          end if
c
c        -----  combine the radial integrals with structure     -----
c        -----  dependent coeff.                                -----
c
          s00 = dumi*dumj
          n=0
          do 600 k=1,npnp
            do 600 l=1,k
              n = n + 1
              if (abs(ckl(k,l)).gt.tol) then
                sum = sum + fp(n)*ckl(k,l)*s00
              end if
  600     continue
c   store the computed integral into g
          g(nn) = g(nn) + sum*fpi
          nn = nn + 1
c   end of i & j angular loops
  900   continue
 1000 continue
      return
      end
**==ecpa14.f
      subroutine ecpa14(fp,jfst1,lbecp1,dcoef1,g,npnp)
c
      implicit real*8 (a-h,o-z)
c
c  ecp type 1 angular integrals
c  icab = 4 general case with 3 different centers <cab>
c  note the actual integral has already been done, we just
c  need to compute the appropriate factors and multiply by the
c  stored integral.
c  fp will contain all the required radial integrals upon entry
c  g will contain the gradient vector upon exit (not symmetrized)
c  lbecp1 contains the packed indeces for the angular integrals
c  dcoef1 contains the packed angular integral values
c  npnp = n+n' gives the maximum angular momentum (kappa) needed
c
      dimension dcoef1(*),lbecp1(9,*),jfst1(*),fp(*),g(*)
c
      common/junk2/cklu(23,12,12)
c
c
      real*8 acx,acy,acz,abx,aby,abz
      common /ecp3/ acx(7),acy(7),acz(7),abx(7),aby(7),abz(7)
c
c
c   iamin/max are the angular function min and max for i
c   ipmin/max are the primitive function min and max for i
c   kf1/kl1 are the first and last primitives forming the lmax ecp
c   and thus the type 1 integral
c
      logical iandj, norm, normi, normj
      common /ecpidx/ q2,iamin,iamax,jamin,jamax,ipmin,ipmax,jpmin,
     *                jpmax,kf1,kl1,llmx,npc,npb,iandj,norm,normi,normj
c
      integer nfst,nx,ny,nz
      common /gbase/ nfst(8),nx(84),ny(84),nz(84)
c
c
      parameter (zero=0.0d+00, one=1.0d+00)
      parameter (sqrt3=1.73205080756888d+00,sqrt5=2.23606797749979d+00)
      parameter (tol=1.0d-10, fpi=12.566370614359d+00)
      parameter (sqrt7=2.64575131106459d+00)
c
c   nn indexes the gradient vector array
c
      nn=0
      mmax=jamax
      do 1000 i=iamin,iamax
        dumi = one
        if (normi) then
          if (i.ge.8.and.i.le.10) then
            dumi = sqrt3
          else if (i.ge.14.and.i.le.20) then
            dumi = sqrt5
            if (i.eq.20) dumi=dumi*sqrt3
          else if (i.ge.24) then
            dumi = sqrt7
            if (i.ge.30) then
              dumi=dumi*sqrt5/sqrt3
              if (i.ge.33) dumi=dumi*sqrt3
            end if
          end if
        end if
        if (iandj) mmax=i
        do 900 j=jamin,mmax
          dumj = one
          if (normj) then
            if (j.ge.8.and.j.le.10) then
              dumj = sqrt3
            else if (j.ge.14.and.j.le.20) then
              dumj = sqrt5
              if (j.eq.20) dumj=dumj*sqrt3
            end if
            else if (j.ge.24) then
              dumj = sqrt7
              if (j.ge.30) then
                dumj=dumj*sqrt5/sqrt3
                if (j.ge.33) dumj=dumj*sqrt3
              end if
          end if
          s00=dumi*dumj*q2
          if (i.ge.j) then
            indx = j+(i*(i-1))/2
          else
            indx = i+(j*(j-1))/2
          end if
          jf = jfst1(indx)
          jl = jfst1(indx+1)-1
          nxc = nx(i)+1
          nyc = ny(i)+1
          nzc = nz(i)+1
          nxb = nx(j)+1
          nyb = ny(j)+1
          nzb = nz(j)+1
c
c  zero out ckl which stores the angular integrals to be used later
c  these zeroing loops originally matched the compute loops 400/390/380
c  below.  however, the value of lamda+mu in loop 200 can be as big
c  as npnp+1, so to avoid core dumps, we zero more of the array out.
c  original dimension of cklu was (22,11,11), changed to (23,12,12).
c
          do 110 k=1,npnp+1
            do 110 l=1,k
              do 110 mu=1,2*l-1
  110     cklu(mu,l,k) = zero
c
          do 200 k=jf,jl
            kappa = lbecp1(1,k)+1
            lamda = lbecp1(2,k)+1
            mu = lbecp1(3,k)
            if(i.ge.j) then
              kx = lbecp1(4,k)
              ky = lbecp1(5,k)
              kz = lbecp1(6,k)
              kxp = lbecp1(7,k)
              kyp = lbecp1(8,k)
              kzp = lbecp1(9,k)
            else
              kxp = lbecp1(4,k)
              kyp = lbecp1(5,k)
              kzp = lbecp1(6,k)
              kx = lbecp1(7,k)
              ky = lbecp1(8,k)
              kz = lbecp1(9,k)
            end if
****
*           write(6,*) '*** k,kappa.lamda,mu,cklu = ', 
*    +                  k,kappa,lamda,mu,
*    +                  cklu(lamda+mu,lamda,kappa)
*           write(6,*) 'kx,ky,kz = ', kx, ky, kz
*           write(6,*) 'nxc,nyc,nzc = ', nxc,nyc,nzc
*           write(6,*) 'kxp,kyp,kzp = ', kxp, kyp, kzp
*           write(6,*) 'nxb,nyb,nzb = ', nxb,nyb,nzb
**** 
            cklu(lamda+mu,lamda,kappa) = cklu(lamda+mu,lamda,kappa)+
     *          dcoef1(k)*acx(nxc-kx) *acy(nyc-ky) *acz(nzc-kz)
     *                   *abx(nxb-kxp)*aby(nyb-kyp)*abz(nzb-kzp)
  200     continue
c
c  zero out sum which is used to collect each total integral
c
          sum = zero
          n=0
          do 400 k=1,npnp
            do 390 l=1,k
              do 380 mu=1,2*l-1
                ckltem=cklu(mu,l,k)
                n=n+1
                if (abs(ckltem).gt.tol) then
                  sum = sum + fp(n)*ckltem
                end if
  380         continue
  390       continue
  400     continue
c
c   store the computed integral into g
c
          nn = nn + 1
          g(nn) = g(nn) + sum*s00*fpi
  900   continue
 1000 continue
c
c   end of i & j angular loops
c
      return
      end
**==ecpa21.f
      subroutine ecpa21(g,ecoef,coefi,coefj,npnp)
      implicit real*8 (a-h,o-z)
c
c  ecp type 2 integrals for the special case of all
c  three functions on the same center <aaa> (ie. icab = 1)
c  g will contain the gradient vector upon exit (not symmetrized)
c  npnp = n+n' gives the maximum angular momentum (kappa) needed
c  note: l shells can not be done with this routine due to the coeff
c  of the prims multiplied in the radial integral section
c
      dimension g(*), coefi(*), coefj(*), ecoef(*)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
c
      real*8 zetc,cax,cay,caz,ca,xca,yca,zca
      real*8 zetb,bax,bay,baz,ba,xba,yba,zba
      real*8 phase,dax,day,daz,da,xda,yda,zda,xint
      integer kcntr
      common /ecp1/ zetc,cax,cay,caz,ca,xca,yca,zca,
     +              zetb,bax,bay,baz,ba,xba,yba,zba,
     +              phase,dax,day,daz,da,xda,yda,zda,xint,
     +              kcntr
c
      common /ecp2  / clp(400),zlp(400),nlp(400),kfirst(maxat,6),
     *                klast(maxat,6),lmax(maxat),lpskip(maxat),
     *                izcore(maxat)
c
      integer ncoef1,ncoef2,j1len,j2len,llim,nlim,ntlim,j4len
      common /ecpdim/ ncoef1,ncoef2,j1len,j2len,llim,nlim,ntlim,j4len
c
c   iamin/max are the angular function min and max for i
c   ipmin/max are the primitive function min and max for i
c   kf1/kl1 are the first and last primitives forming the lmax ecp
c      and thus the type 1 integral
c
      logical iandj, norm, normi, normj
      common /ecpidx/ q2,iamin,iamax,jamin,jamax,ipmin,ipmax,jpmin,
     *                jpmax,kf1,kl1,llmx,npc,npb,iandj,norm,normi,normj
c
      parameter(zero=0.0d+00, one=1.0d+00)
      parameter(sqrt3=1.73205080756888d+00,sqrt5=2.23606797749979d+00)
      parameter (sqrt7=2.64575131106459d+00)
c
c   nbc = nxc + nyc + nzc + nbx + nby + nbz
      nbc = npnp - 1
c
c   the total n for the two shells must be even or the integral is
c   identically zero
c
      if (mod(nbc,2).eq.1) return
c
      mmax=jamax
c
c   loop over the ecp potentials on center k
c
      do 1500 ll=2,llmx
        kf = kfirst(kcntr,ll)
        kl = klast(kcntr,ll)
        l=ll-2
        nlm1 = (l*(l+1))
c
c   loop over the primitives summing the radial integrals in fpsum
c
        fpsum = zero
        igii=1
        do 400 ig = ipmin,ipmax
          dum = coefi(igii)
          igii=igii+1
          zetc = ex(ig)
          jgjj=1
          do 300 jg = jpmin,jpmax
            dum2 = dum*coefj(jgjj)
            jgjj=jgjj+1
            zetb = ex(jg)
            zetcb = zetc + zetb
            fptemp = zero
            do 200 k=kf,kl
              zeta = zetcb+zlp(k)
              fptemp = fptemp + fa(nbc+nlp(k), zeta)*clp(k)
 200        continue
            fpsum = fpsum + fptemp*dum2
 300      continue
 400    continue
c
c   now loop over the angular functions, retrieving the angular
c   integral, multiplying by the radial sum and store in g
c
        nn = 1
        do 1000 i=iamin,iamax
          dumi = one
          if (normi) then
            if (i.ge.8.and.i.le.10) then
              dumi = sqrt3
            else if (i.ge.14.and.i.le.20) then
              dumi = sqrt5
              if (i.eq.20) dumi=dumi*sqrt3
            else if (i.ge.24) then
              dumi = sqrt7
              if (i.ge.30) then
                dumi=dumi*sqrt5/sqrt3
                if (i.ge.33) dumi=dumi*sqrt3
              end if
            end if
          end if
          dumi = dumi*q2*fpsum
          if (iandj) mmax=i
          do 900 j=jamin,mmax
            dumj = one
            if (normj) then
              if (j.ge.8.and.j.le.10) then
                dumj = sqrt3
              else if (j.ge.14.and.j.le.20) then
                dumj = sqrt5
                if (j.eq.20) dumj=dumj*sqrt3
              else if (j.ge.24) then
                dumj = sqrt7
                if (j.ge.30) then
                  dumj=dumj*sqrt5/sqrt3
                  if (j.ge.33) dumj=dumj*sqrt3
                end if
              end if
            end if
            eco = zero
            do 500 m=-l,l
              nlm = (nlm1-m)*ntlim
              eco=eco+ecoef(nlm+i)*ecoef(nlm+j)
 500        continue
            g(nn) = g(nn) + eco*dumi*dumj
            nn = nn + 1
  900     continue
 1000   continue
 1500 continue
      return
      end
**==ecpa22.f
      subroutine ecpa22(fp,jfst2,lbecp2,dcoef2,ecoef,g,icab,npnp,l,
     *  zlm,lmf,lmx,lmy,lmz)
      implicit real*8 (a-h,o-z)
c
c  ecp type 2 angular integrals
c  note the actual integral has already been done, we just
c  need to compute the appropriate factors and multiply by the
c  stored integral.
c
c  fp will contain all the required radial integrals upon entry
c  g will contain the gradient vector upon exit (not symmetrized)
c  lbecp2 contains the packed indeces for the angular integrals
c  dcoef2 contains the packed angular integral values
c  npnp = n+n' gives the maximum angular momentum (kappa) needed
c
      dimension dcoef2(*), lbecp2(6,*), jfst2(*)
      dimension fp(*), ckl(11,11), g(*), ecoef(*)
      dimension zfnlm(121),zlm(*),lmf(*),lmx(*),lmy(*),lmz(*)
c
c
      integer ncoef1,ncoef2,j1len,j2len,llim,nlim,ntlim,j4len
      common /ecpdim/ ncoef1,ncoef2,j1len,j2len,llim,nlim,ntlim,j4len
c
c
      real*8 acx,acy,acz,abx,aby,abz
      common /ecp3/ acx(7),acy(7),acz(7),abx(7),aby(7),abz(7)
c
c
      integer nfst,nx,ny,nz
      common /gbase/ nfst(8),nx(84),ny(84),nz(84)
c
c   iamin/max are the angular function min and max for i
c   ipmin/max are the primitive function min and max for i
c   kf1/kl1 are the first and last primitives forming the lmax ecp
c      and thus the type 1 integral
c
      logical iandj, norm, normi, normj
      common /ecpidx/ q2,iamin,iamax,jamin,jamax,ipmin,ipmax,jpmin,
     *                jpmax,kf1,kl1,llmx,npc,npb,iandj,norm,normi,normj
c
      parameter(zero=0.0d+00, one=1.0d+00)
      parameter(sqrt3=1.73205080756888d+00,sqrt5=2.23606797749979d+00)
      parameter (sqrt7=2.64575131106459d+00)
      parameter (tol=1.0d-10, fpi=12.566370614359d+00)
c
      l2pl=l*(l+1)
      call zfn(zfnlm, npnp-1,zlm,lmf,lmx,lmy,lmz)
c   nn indexes the gradient vector array
      nn=0
      mmax=jamax
      do 1000 i=iamin,iamax
        dumi = one
        if (normi) then
          if (i.ge.8.and.i.le.10) then
            dumi = sqrt3
          else if (i.ge.14.and.i.le.20) then
            dumi = sqrt5
            if (i.eq.20) dumi=dumi*sqrt3
          else if (i.ge.24) then
            dumi = sqrt7
            if (i.ge.30) then
              dumi=dumi*sqrt5/sqrt3
              if (i.ge.33) dumi=dumi*sqrt3
            end if
          end if
        end if
        if (iandj) mmax=i
        do 900 j=jamin,mmax
          dumj = one
          if (normj) then
            if (j.ge.8.and.j.le.10) then
              dumj = sqrt3
            else if (j.ge.14.and.j.le.20) then
              dumj = sqrt5
              if (j.eq.20) dumj=dumj*sqrt3
            else if (j.ge.24) then
              dumj = sqrt7
              if (j.ge.30) then
                dumj=dumj*sqrt5/sqrt3
                if (j.ge.33) dumj=dumj*sqrt3
              end if
            end if
          end if
          nxc = nx(i)
          nyc = ny(i)
          nzc = nz(i)
          nxb = nx(j)
          nyb = ny(j)
          nzb = nz(j)
c  zero out ckl which stores the angular integrals to be used later
          do 100 k=1,npnp
            do 100 ll=1,k
  100     ckl(k,ll) = zero
c
c    now calculate the angular integral which does not depend on primitives
c
          if (icab.eq.2) then
            kappat = nxc+nyc+nzc+1
            do 200 m=l,-l,-1
              nlm = (l2pl-m)*ntlim
              ct = ecoef(nlm+i)
              if(abs(ct).ge.tol) then
                jf = jfst2(nlm+j)
                jl = jfst2(nlm+j+1)-1
                do 180 jj=jf,jl
                  lamda = lbecp2(1,jj)
                  kappa = lbecp2(2,jj)
                  mu = lbecp2(3,jj)
                  kx = lbecp2(4,jj)
                  ky = lbecp2(5,jj)
                  kz = lbecp2(6,jj)
                  ckl(kappa+kappat,lamda+1)=ckl(kappa+kappat,lamda+1)
     *               +ct*dcoef2(jj)*zfnlm(lamda*(lamda+1)-mu+1)*
     *                abx(nxb-kx+1)*aby(nyb-ky+1)*abz(nzb-kz+1)
  180           continue
              end if
  200       continue
c
          else
            kappat=nxb+nyb+nzb+1
            do 300 m=l,-l,-1
              nlm = (l2pl-m)*ntlim
              ct = ecoef(nlm+j)
              if(abs(ct).ge.tol) then
                jf = jfst2(nlm+i)
                jl = jfst2(nlm+i+1)-1
                do 280 jj=jf,jl
                  lamda = lbecp2(1,jj)
                  kappa = lbecp2(2,jj)
                  mu = lbecp2(3,jj)
                  kx = lbecp2(4,jj)
                  ky = lbecp2(5,jj)
                  kz = lbecp2(6,jj)
                  ckl(kappa+kappat,lamda+1)=ckl(kappa+kappat,lamda+1)
     *               +ct*dcoef2(jj)*zfnlm(lamda*(lamda+1)-mu+1)*
     *                acx(nxc-kx+1)*acy(nyc-ky+1)*acz(nzc-kz+1)
  280           continue
              end if
  300       continue
          end if
c
          sum = zero
          s00 = dumi*dumj*q2
c
c     combine the angular and the radial integrals
c
          n=0
          do 600 k=1,npnp
            do 600 ll=1,k
              n = n + 1
              if (abs(ckl(k,ll)).gt.tol) then
                sum = sum + fp(n)*ckl(k,ll)*s00
              end if
  600     continue
c   store the computed integral into g
          nn = nn + 1
          g(nn) = g(nn) + sum*fpi
c   end of i & j angular loops
  900   continue
 1000 continue
      return
      end
**==ecpa24.f
      subroutine ecpa24(fp,jfst2, lbecp2,dcoef2,g,npnp,l,zlm,lmf,lmx,
     *    lmy,lmz)
      implicit real*8 (a-h,o-z)
c
c  ecp type 2 angular integrals
c  note the actual integral has already been done, we just
c  need to compute the appropriate factors and multiply by the
c  stored integral.
c
c  fp will contain all the required radial integrals upon entry
c  g will contain the gradient vector upon exit (not symmetrized)
c  lbecp2 contains the packed indeces for the angular integrals
c  dcoef2 contains the packed angular integral values
c  npnp = n+n' gives the maximum angular momentum (kappa) needed
c  l = the angular momentum of the core represented by the ecp
c
c
      dimension dcoef2(*), lbecp2(6,*), cpq(11,11,11)
      dimension fp(*),cklc(11,11),cklb(11,11), g(*), jfst2(*)
      dimension zfnlmc(121),zfnlmb(121),zlm(*),lmf(*),lmx(*),lmy(*)
      dimension lmz(*)
c
c
      integer ncoef1,ncoef2,j1len,j2len,llim,nlim,ntlim,j4len
      common /ecpdim/ ncoef1,ncoef2,j1len,j2len,llim,nlim,ntlim,j4len
c
c
      real*8 zetc,cax,cay,caz,ca,xca,yca,zca
      real*8 zetb,bax,bay,baz,ba,xba,yba,zba
      real*8 phase,dax,day,daz,da,xda,yda,zda,xint
      integer kcntr
      common /ecp1/ zetc,cax,cay,caz,ca,xca,yca,zca,
     +              zetb,bax,bay,baz,ba,xba,yba,zba,
     +              phase,dax,day,daz,da,xda,yda,zda,xint,
     +              kcntr
c
c
      real*8 acx,acy,acz,abx,aby,abz
      common /ecp3/ acx(7),acy(7),acz(7),abx(7),aby(7),abz(7)
c
c
      integer nfst,nx,ny,nz
      common /gbase/ nfst(8),nx(84),ny(84),nz(84)
c
c
c   iamin/max are the angular function min and max for i
c   ipmin/max are the primitive function min and max for i
c   kf1/kl1 are the first and last primitives forming the lmax ecp
c      and thus the type 1 integral
c
      logical iandj, norm, normi, normj
      common /ecpidx/ q2,iamin,iamax,jamin,jamax,ipmin,ipmax,jpmin,
     *                jpmax,kf1,kl1,llmx,npc,npb,iandj,norm,normi,normj
c
      common /zfncm / x,y,z
c
      parameter (zero=0.0d+00, one=1.0d+00)
      parameter (sqrt5=2.23606797749979d+00)
      parameter (sqrt3=1.73205080756888d+00,fpisq=157.91367041743d+00)
      parameter (sqrt7=2.64575131106459d+00)
      parameter (tol=1.0d-10)
c
      save zfnlmb,zfnlmc
c
      npcpl = npc+l
      npbpl = npb+l
      np1 = max(npc,npb)
      np1pl = np1 + l
      l1max=max(1,l+1)
      l2pl=l*(l+1)
      x = xca
      y = yca
      z = zca
      call zfn(zfnlmc, npcpl-1,zlm,lmf,lmx,lmy,lmz)
      x = xba
      y = yba
      z = zba
      call zfn(zfnlmb, npbpl-1,zlm,lmf,lmx,lmy,lmz)
c   nn indexes the integral vector array
      nn=0
      mmax=jamax
      do 1000 i=iamin,iamax
        dumi = one
        if (normi) then
          if (i.ge.8.and.i.le.10) then
            dumi = sqrt3
          else if (i.ge.14.and.i.le.20) then
            dumi = sqrt5
            if (i.eq.20) dumi=dumi*sqrt3
          else if (i.ge.24) then
            dumi = sqrt7
            if (i.ge.30) then
              dumi=dumi*sqrt5/sqrt3
              if (i.ge.33) dumi=dumi*sqrt3
            end if
          end if
        end if
        if (iandj) mmax=i
        do 900 j=jamin,mmax
          dumj = one
          if (normj) then
            if (j.ge.8.and.j.le.10) then
              dumj = sqrt3
            else if (j.ge.14.and.j.le.20) then
              dumj = sqrt5
              if (j.eq.20) dumj=dumj*sqrt3
            else if (j.ge.24) then
              dumj = sqrt7
              if (j.ge.30) then
                dumj=dumj*sqrt5/sqrt3
                if (j.ge.33) dumj=dumj*sqrt3
              end if
            end if
          end if
          s00=dumi*dumj*q2
          nxc = nx(i)+1
          nyc = ny(i)+1
          nzc = nz(i)+1
          nxb = nx(j)+1
          nyb = ny(j)+1
          nzb = nz(j)+1
c
c  zero out ckl which stores the angular integrals to be used later
c
          do 100 kk=1,2*np1-1
            do 100 ll=1,npcpl
              do 100 jj=1,npbpl
  100     cpq(jj,ll,kk) = zero
c
c    now calculate the angular integral which does not depend on primitives
c
          do 300 m=l,-l,-1
            do 110 kk=1,np1
              do 110 ll=1,np1pl
                cklc(ll,kk)=zero
                cklb(ll,kk)=zero
  110       continue
            nlm = (l2pl-m)*ntlim
            jf = jfst2(nlm+i)
            jl = jfst2(nlm+i+1)-1
cbb            if (nlm+i.gt.j2len.or.jl.gt.ncoef2) then
cbb               write(6,*) 'ncoef2 out a24',nlm,i,ncoef2,jf,jl
cbb               stop
cbb            end if
            do 150 jj=jf,jl
              lamda = lbecp2(1,jj)+1
              kappa = lbecp2(2,jj)+1
              mu = lbecp2(3,jj)
              kx = lbecp2(4,jj)
              ky = lbecp2(5,jj)
              kz = lbecp2(6,jj)
              cklc(lamda,kappa)=cklc(lamda,kappa)
     *           +dcoef2(jj)*zfnlmc(lamda*(lamda-1)-mu+1)*
     *            acx(nxc-kx)*acy(nyc-ky)*acz(nzc-kz)
  150       continue
            nlm = (l2pl-m)*ntlim
            jf = jfst2(nlm+j)
            jl = jfst2(nlm+j+1)-1
cbb            if (nlm+j.gt.j2len.or.jl.gt.ncoef2) then
cbb               write(6,*) 'ncoef2 out a24',nlm,j,ncoef2,jf,jl
cbb               stop
cbb            end if
            do 180 jj=jf,jl
              lamda = lbecp2(1,jj)+1
              kappa = lbecp2(2,jj)+1
              mu = lbecp2(3,jj)
              kx = lbecp2(4,jj)
              ky = lbecp2(5,jj)
              kz = lbecp2(6,jj)
              cklb(lamda,kappa)=cklb(lamda,kappa)
     *           +dcoef2(jj)*zfnlmb(lamda*(lamda-1)-mu+1)*
     *            abx(nxb-kx)*aby(nyb-ky)*abz(nzb-kz)
  180       continue
c        -----  multiply coeff of the two projection halves     -----
c        -----  and store.                                      -----
            do 280 kk=0,np1-1
              do 280 ll=1,npcpl
                ct1 = cklc(ll,kk+1)
                do 270 kk2=1,np1
                  do 270 ll2=1,npbpl
                    cpq(ll2,ll,kk+kk2) = cpq(ll2,ll,kk+kk2)+
     *                    ct1*cklb(ll2,kk2)
  270           continue
  280       continue
  300     continue
c  zero out sum which is used to collect each total integral
          sum = zero
c
c        -----  combine the radial integrals with structure     -----
c        -----  dependent coeff.                                -----
c
c   do kappa=0 first since it only has values on the diagonal
            do 550 ll=1,l1max
              sum = sum+fp(ll)*cpq(ll,ll,1)
  550       continue
c   if there are kappa values larger than 0 do those now
            if (npnp.gt.1) then
              jstart = 1
              n=6
c   do higher kappa values. note: every other angular integral will be
c   zero by symmetry and is therefor skipped
              do 601 kk=2,npnp
                jstart = 1-jstart
                lstart = jstart
                do 591 jj=1,npcpl
                  lstart = 1-lstart
                  do 581 ll=lstart+1,npbpl,2
                    n=n+1
                    if (abs(cpq(ll,jj,kk)).gt.tol) then
                      sum = sum + fp(n)*cpq(ll,jj,kk)
                    end if
  581             continue
  591           continue
  601         continue
           end if
c   store the computed integral into g
          nn = nn + 1
          g(nn) = g(nn) + sum*s00*fpisq
c   end of i & j angular loops
  900   continue
 1000 continue
      return
      end
**==ecpd14.f
      subroutine ecpd14(fp,fp2,coefi,coefi2,coefj,npnp,zlm,lmf,lmx,lmy,
     *     lmz)
      implicit real*8 (a-h,o-z)
c
c  type 1 radial integrals for the given set of shells
c  shell set types <b|a|c> (ie icab=4)
c
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
      dimension fp(*), fip(78), zfnlm(121), coefi(*),coefj(*)
      dimension coefi2(*), fp2(*),zlm(*),lmf(*),lmx(*),lmy(*),lmz(*)
c
      common /ecp2  / clp(400),zlp(400),nlp(400),kfirst(maxat,6),
     *                klast(maxat,6),lmax(maxat),lpskip(maxat),
     *                izcore(maxat)
c
      logical normi, normj, norm, iandj
      common /ecpidx/ q2,iamin,iamax,jamin,jamax,ipmin,ipmax,jpmin,
     *                jpmax,kf1,kl1,llmx,npc,npb,iandj,norm,normi,normj
c
c
      real*8 bmcx,bmcy,bmcz,bpcx,bpcy,bpcz,cbsq,ax,ay,az
      logical candb
      common /ecp4/ bmcx,bmcy,bmcz,bpcx,bpcy,bpcz,cbsq,ax,ay,az,
     +              candb
c
      common /ficmn / alf,xi,xp0,xp1
      common /zfncm / x,y,z
c
      parameter(zero=0.0d+00, one=1.0d+00, two=2.0d+00,half=0.5d+00)
      parameter(tol6=1.0d-06,onds4p=0.28209479177388d+00)
c
c    maximum number of fp integrals needed for the given set of shells
c
      npnpmx = (npnp*(npnp+1)*(2*npnp+1))/6
c
c   first zero out the integral array
c
      do 100 np=1,npnpmx
       fp2(np) = zero
 100  fp(np) = zero
c
c   loop over the primitives storing the radial integrals in fp
c
      igii=1
      do 700 ig = ipmin,ipmax
        dum=coefi(igii)
        dumd=coefi2(igii)
        igii=igii+1
        zetc = ex(ig)
        jgjj=1
        do 600 jg = jpmin,jpmax
          dum2=dum*coefj(jgjj)
          dumd2=dumd*coefj(jgjj)
          jgjj=jgjj+1
          zetb = ex(jg)
          zetcb = zetc + zetb
c  transform the gaussians to center d actually c and b may be the
c  same as long as neither is the same as a. if c and b are different
c  then a phase factor is also required
          s = one/zetcb
          rat = (zetb-zetc)*s
          dax = half*(bpcx+rat*bmcx) - ax
          day = half*(bpcy+rat*bmcy) - ay
          daz = half*(bpcz+rat*bmcz) - az
          da = sqrt(dax*dax+day*day+daz*daz)
          if (.not.candb) then
            dum2 = dum2*exp(-(zetc*zetb*s)*cbsq)
            dumd2 = dumd2*exp(-(zetc*zetb*s)*cbsq)
          end if
c  if center d and a both lie on the origin then all integrals for
c  lamda<>0 are 0 so special case that code
          if (da.ge.tol6) then
            alfa1 = zetcb*da
            xalfa1 = alfa1*da
            alfi1=one/(two*alfa1)
            alfi = alfi1
            x = dax/da
            y = day/da
            z = daz/da
            call zfn(zfnlm, npnp-1,zlm,lmf,lmx,lmy,lmz)
            xp0 = exp(-xalfa1)
            do 300 k=kf1,kl1
              xi=one/sqrt(zetcb+zlp(k))
              alf=alfa1*xi
              xp1=exp(-xalfa1+alf*alf)
              nlpk = nlp(k)
              clpk = clp(k)*dum2
              clpk2 = clp(k)*dumd2
              call fiecp(fip, alfi, nlpk, npnp-1)
              n=0
              nn=0
              do 200 kk=1,npnp
                do 200 l=1,kk
                  nn=nn+1
                  fiptem = fip(nn)*clpk
                  fipt2 = fip(nn)*clpk2
                  kkll = l*(l-1)+1+l
                  do 200 mu=1,2*l-1
                    n=n+1
                    fp(n) = fp(n) + fiptem*zfnlm(kkll-mu)
                    fp2(n) = fp2(n) + fipt2*zfnlm(kkll-mu)
 200          continue
 300        continue
          else
c    special case when d=a=0. there are only l=0 integrals, but
c    they must be stored the same as above so index them into fp.
c    they also must be multiplied by the 0 spherical hamonic 1/sqrt(4pi)
            do 500 k=kf1,kl1
              zeta = zetcb+zlp(k)
              nlpk = nlp(k)
              clpk = clp(k)*dum2*onds4p
              clpk2 = clp(k)*dumd2*onds4p
              do 400 n=1,npnp
                nk = (n*(n-1)*(2*n-1))/6 + 1
                fatemp = fa(n+nlpk-1, zeta)
                fp(nk) = fp(nk)+fatemp*clpk
                fp2(nk) = fp2(nk)+fatemp*clpk2
 400          continue
 500        continue
          end if
 600    continue
 700  continue
      end
**==ecpini.f
      subroutine ecpini(lmf,lmx,lmy,lmz)
      implicit real*8 (a-h,o-z)
      dimension lmf(122),lmx(581),lmy(581),lmz(581)
      dimension ilmf(122),ilmx(581),ilmy(581),ilmz(581)
      dimension infst(8),inx(84),iny(84),inz(84)
c
      integer nfst,nx,ny,nz
      common /gbase/ nfst(8),nx(84),ny(84),nz(84)
c
c lmf (first) point to the first parts of the 
c spherical harmonic coef table (zlm) and to the corresponding powers
c of x, y, and z stored in lmx, lmy, lmz (good through l=10)
      data ilmf/1, 2,3,4, 5,7,8,10,11, 12,14,16,18,20,22,23, 
     *25,28,30,34,36,39,41,43,45, 47,50,53,57,61,64,67,70,72,76,78,
     *81,85,88,94,98,104,107,111,114,117,121,125,128,
     *131,135,139,145,151,157,163,167,171,175,178,184,188,194,197,
     *201,206,210,218,224,233,239,247,251,256,260,264,270,276,282,
     *288,292, 296,301,306,314,322,331,340,348,356,361,366,371,375,383,
     *389,398,404,412,416, 421,427,432,442,450,462,471,483,491,501,506,
     *512,517,522,530,538,547,556,564,572,577, 582/
      data ilmx/0, 1,0,0, 2,0,1,0,0,0,1, 3,1,2,0,1,1,0,0,0,0,1,2,0, 
     *4,2,0,3,1,2,0,0,2,1,1,0,0,0,0,0,1,1,2,0,3,1, 5,3,1,0,2,4,3,1,1,3,
     *2,0,2,0,1,1,1,0,0,0,0,0,0,1,1,2,2,0,0,3,1,4,2,0, 6,4,2,0,5,3,1,4,
     *2,0,4,2,0,3,1,1,3,2,0,2,0,2,0,1,1,1,0,0,0,0,0,0,0,1,1,1,2,0,2,0,3,
     *1,3,1,4,2,0,5,3,1, 7,3,5,1,6,4,0,2,5,3,1,5,3,1,4,2,0,4,2,0,3,1,3,
     *1,3,1,2,0,2,0,2,0,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,2,0,2,0,2,0,3,1,3,
     *1,4,2,0,4,2,0,5,3,1,0,4,2,6, 8,6,4,2,0,7,5,3,1,6,4,2,0,6,4,2,0,5,
     *3,1,5,3,1,4,0,2,4,0,2,4,0,2,3,1,3,1,3,1,2,0,2,0,2,0,2,0,1,1,1,1,0,
     *0,0,0,0,0,0,0,0,1,1,1,1,2,0,2,0,2,0,3,1,3,1,3,1,4,2,0,4,2,0,5,3,1,
     *5,3,1,6,4,2,0,7,5,3,1, 9,7,5,3,1,8,6,4,2,0,7,5,3,1,7,5,3,1,6,4,2,
     *0,6,4,2,0,5,3,1,5,3,1,5,3,1,4,0,2,4,0,2,4,0,2,3,1,3,1,3,1,3,1,2,0,
     *2,0,2,0,2,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,2,0,2,0,2,0,2,0,
     *3,1,3,1,3,1,4,2,0,4,2,0,4,2,0,5,3,1,5,3,1,6,4,2,0,6,4,2,0,7,5,3,1,
     *8,6,4,2,0, 10,8,6,4,2,0,9,7,5,3,1,8,6,4,2,0,8,6,4,2,0,7,5,3,1,7,5,
     *3,1,6,4,2,0,6,4,2,0,6,4,2,0,5,3,1,5,3,1,5,3,1,4,0,2,4,0,2,4,0,2,4,
     *0,2,3,1,3,1,3,1,3,1,2,0,2,0,2,0,2,0,2,0,1,1,1,1,1,0,0,0,0,0,0,0,
     *0,0,0,0,1,1,1,1,1,2,0,2,0,2,0,2,0,3,1,3,1,3,1,3,1,0,2,4,0,2,4,0,2,
     *4,5,3,1,5,3,1,5,3,1,6,4,2,0,6,4,2,0,7,5,3,1,7,5,3,1,8,6,4,2,0,9,7,
     *5,3,1/
      data ilmy/0, 0,0,1, 0,2,0,0,0,1,1, 0,2,0,2,0,0,0,0,1,1,1,1,3,
     *0,2,4,0,2,0,2,2,0,0,0,0,0,0,1,1,1,1,1,3,1,3, 0,2,4,4,2,0,0,2,2,0,
     *0,2,0,2,0,0,0,0,0,0,1,1,1,1,1,1,1,3,3,1,3,1,3,5, 0,2,4,6,0,2,4,0,
     *2,4,0,2,4,0,2,2,0,0,2,0,2,0,2,0,0,0,0,0,0,0,1,1,1,1,1,1,1,3,1,3,1,
     *3,1,3,1,3,5,1,3,5, 0,4,2,6,0,2,6,4,0,2,4,0,2,4,0,2,4,0,2,4,0,2,0,
     *2,0,2,0,2,0,2,0,2,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,3,1,3,1,3,1,3,1,
     *3,1,3,5,1,3,5,1,3,5,7,3,5,1, 0,2,4,6,8,0,2,4,6,0,2,4,6,0,2,4,6,0,
     *2,4,0,2,4,0,4,2,0,4,2,0,4,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,0,0,0,0,
     *0,0,0,0,1,1,1,1,1,1,1,1,1,3,1,3,1,3,1,3,1,3,1,3,1,3,5,1,3,5,1,3,5,
     *1,3,5,1,3,5,7,1,3,5,7, 0,2,4,6,8,0,2,4,6,8,0,2,4,6,0,2,4,6,0,2,4,
     *6,0,2,4,6,0,2,4,0,2,4,0,2,4,0,4,2,0,4,2,0,4,2,0,2,0,2,0,2,0,2,0,2,
     *0,2,0,2,0,2,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,3,1,3,1,3,1,3,
     *1,3,1,3,1,3,1,3,5,1,3,5,1,3,5,1,3,5,1,3,5,1,3,5,7,1,3,5,7,1,3,5,7,
     *1,3,5,7,9, 0,2,4,6,8,10,0,2,4,6,8,0,2,4,6,8,0,2,4,6,8,0,2,4,6,0,2,
     *4,6,0,2,4,6,0,2,4,6,0,2,4,6,0,2,4,0,2,4,0,2,4,0,4,2,0,4,2,0,4,2,0,
     *4,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,0,0,0,0,0,0,0,0,0,0,1,
     *1,1,1,1,1,1,1,1,1,1,3,1,3,1,3,1,3,1,3,1,3,1,3,1,3,5,3,1,5,3,1,5,3,
     *1,1,3,5,1,3,5,1,3,5,1,3,5,7,1,3,5,7,1,3,5,7,1,3,5,7,1,3,5,7,9,1,3,
     *5,7,9/
      data ilmz/0, 0,1,0, 0,0,1,2,0,1,0, 0,0,1,1,2,0,3,1,2,0,1,0,0,
     *0,0,0,1,1,2,2,0,0,3,1,4,2,0,3,1,2,0,1,1,0,0, 0,0,0,1,1,1,2,2,0,0,
     *3,3,1,1,4,2,0,5,3,1,4,2,0,3,1,2,0,2,0,1,1,0,0,0, 0,0,0,0,1,1,1,2,
     *2,2,0,0,0,3,3,1,1,4,4,2,2,0,0,5,3,1,6,4,2,0,5,3,1,4,2,0,3,3,1,1,2,
     *2,0,0,1,1,1,0,0,0, 0,0,0,0,1,1,1,1,2,2,2,0,0,0,3,3,3,1,1,1,4,4,2,
     *2,0,0,5,5,3,3,1,1,6,4,2,0,7,5,3,1,6,4,2,0,5,3,1,4,4,2,2,0,0,3,3,1,
     *1,2,2,2,0,0,0,1,1,1,0,0,0,0, 0,0,0,0,0,1,1,1,1,2,2,2,2,0,0,0,0,3,
     *3,3,1,1,1,4,4,4,2,2,2,0,0,0,5,5,3,3,1,1,6,6,4,4,2,2,0,0,7,5,3,1,8,
     *6,4,2,0,7,5,3,1,6,4,2,0,5,5,3,3,1,1,4,4,2,2,0,0,3,3,3,1,1,1,2,2,2,
     *0,0,0,1,1,1,1,0,0,0,0, 0,0,0,0,0,1,1,1,1,1,2,2,2,2,0,0,0,0,3,3,3,
     *3,1,1,1,1,4,4,4,2,2,2,0,0,0,5,5,5,3,3,3,1,1,1,6,6,4,4,2,2,0,0,7,7,
     *5,5,3,3,1,1,8,6,4,2,0,9,7,5,3,1,8,6,4,2,0,7,5,3,1,6,6,4,4,2,2,0,0,
     *5,5,3,3,1,1,4,4,4,2,2,2,0,0,0,3,3,3,1,1,1,2,2,2,2,0,0,0,0,1,1,1,1,
     *0,0,0,0,0, 0,0,0,0,0,0,1,1,1,1,1,2,2,2,2,2,0,0,0,0,0,3,3,3,3,1,1,
     *1,1,4,4,4,4,2,2,2,2,0,0,0,0,5,5,5,3,3,3,1,1,1,6,6,6,4,4,4,2,2,2,0,
     *0,0,7,7,5,5,3,3,1,1,8,8,6,6,4,4,2,2,0,0,9,7,5,3,1,10,8,6,4,2,0,9,
     *7,5,3,1,8,6,4,2,0,7,7,5,5,3,3,1,1,6,6,4,4,2,2,0,0,5,5,5,3,3,3,1,1,
     *1,4,4,4,2,2,2,0,0,0,3,3,3,3,1,1,1,1,2,2,2,2,0,0,0,0,1,1,1,1,1,0,0,
     *0,0,0/
      data infst/1,2,5,11,21,36,57,85/
      data inx/0,1,0,0,2,0,0,1,1,0,3,0,0,2,2,1,0,1,0,1,4,0,0,3,3,1,0,
     *1,0,2,2,0,2,1,1, 5,0,0,4,4,1,0,1,0,3,3,2,0,2,0,3,1,1,2,2,1,
     *6,0,0,5,5,1,0,1,0,4,4,2,0,2,0,4,1,1,3,3,0,3,3,2,1,2,1,2/,
     1     iny/0,0,1,0,0,2,0,1,0,1,0,3,0,1,0,2,2,0,1,1,0,4,0,1,0,3,3,
     *0,1,2,0,2,1,2,1, 0,5,0,1,0,4,4,0,1,2,0,3,3,0,2,1,3,1,2,1,2,
     *0,6,0,1,0,5,5,0,1,2,0,4,4,0,2,1,4,1,3,0,3,2,1,3,3,1,2,2/,
     2     inz/0,0,0,1,0,0,2,0,1,1,0,0,3,0,1,0,1,2,2,1,0,0,4,0,1,0,1,
     *3,3,0,2,2,1,1,2, 0,0,5,0,1,0,1,4,4,0,2,0,2,3,3,1,1,3,1,2,2,
     *0,0,6,0,1,0,1,5,5,0,2,0,2,4,4,1,1,4,0,3,3,1,2,1,2,3,3,2/
      do 100 i = 1,122
         lmf(i) = ilmf(i)
  100 continue
c
      do 200 j = 1,581
         lmx(j) = ilmx(j)
         lmy(j) = ilmy(j)
         lmz(j) = ilmz(j)
  200 continue
c
      do 300 k = 1,8
         nfst(k) = infst(k)
  300 continue
c
      do 400 l = 1,84
         nx(l) = inx(l)
         ny(l) = iny(l)
         nz(l) = inz(l)
  400 continue
      return
      end
**==ecpint.f
      subroutine ecpint(core,hecp,h0,iso,nshels,oatomic,out)
      implicit real*8  (a-h,p-z),integer (i-n), logical(o)
c
c     routine controls the calculation of the local potential.  
c
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
c
      integer ncoef1,ncoef2,j1len,j2len,llim,nlim,ntlim,j4len
      common /ecpdim/ ncoef1,ncoef2,j1len,j2len,llim,nlim,ntlim,j4len
c
c
      dimension core(*),hecp(*),h0(*),iso(nshels,*)
c
      nav = lenwrd()
c
         ldcf1  = 0
         ljln   = ldcf1  +  ncoef1
         llb1   = ljln   + (j1len-1)/nav+1
         ldcf4  = llb1   + (ncoef1*9)/nav+1
         ldcf2  = ldcf4  +  j4len
         lj2n   = ldcf2  +  ncoef2
         llb2   = lj2n   + (j2len-1)/nav+1
         lfpqr  = llb2   + (ncoef2*6)/nav
         lzlm   = lfpqr  +  15625
         llmf   = lzlm   + 581
         llmx   = llmf   + 122/nav
         llmy   = llmx   + 582/nav
         llmz   = llmy   + 582/nav
         last   = llmz   + 582/nav
c
         ldcf1 = igmem_alloc(last)
c
         ljln   = ldcf1  +  ncoef1
         llb1   = ljln   + (j1len-1)/nav+1
         ldcf4  = llb1   + (ncoef1*9)/nav+1
         ldcf2  = ldcf4  +  j4len
         lj2n   = ldcf2  +  ncoef2
         llb2   = lj2n   + (j2len-1)/nav+1
         lfpqr  = llb2   + (ncoef2*6)/nav
         lzlm   = lfpqr  +  15625
         llmf   = lzlm   + 581
         llmx   = llmf   + 122/nav
         llmy   = llmx   + 582/nav
         llmz   = llmy   + 582/nav
      
         call ecpint2(hecp,h0,core(ldcf1),
     +   core(ljln),core(llb1),core(ldcf4),core(ldcf2),core(lj2n),
     +   core(llb2),core(lfpqr),core(lzlm),core(llmf),core(llmx),
     +   core(llmy),core(llmz),nx,iso,nshels,oatomic,out)
c
         call gmem_free(ldcf1)
c
      return
      end
**==ecpint2.f
      subroutine ecpint2(hecp,h0,dcoef1,jfst1,lbecp1,dcoef4,dcoef2,
     *                  jfst2,lbecp2,fpqr,zlm,lmf,lmx,lmy,lmz,l2,
     *                  mapshl,nshels,oatomic,dbug)
      implicit real*8 (a-h,o-z)
      logical dbug,canda,aandb
      logical oatomic ! for the atomic start up: calculate one-centre
                      ! terms only
c
      dimension hecp(l2),h0(l2),jfst1(*),jfst2(*),fpqr(25,25,25)
      dimension dcoef1(*),lbecp1(9,*),dcoef2(*),lbecp2(6,*),dcoef4(*)
      dimension zlm(*),lmf(*),lmx(*),lmy(*),lmz(*)
      dimension mapshl(nshels,*)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c      g-integral storage - must be imax*jmax(15*28 for g-hess)
c      7 - for dimension of nlim (for g-hess)
c      13 - dimensioned up to nmax=2*nlim-1
c      need nx etc dimensioned up to 35 for g energy 74 for g hess
c      fp is used for radial integral storage, must be at least 2*11**3
c
c  local storage
      dimension g(420),mi(48),iang(35)
      dimension fp(2662), coefi(mxprms), coefj(mxprms)
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
c
      character *8 zcom,ztitle,zaname,ztag,zsymm,zgroup
      character *10 zbflab
      character *8 zscftp,zruntp,zguess,zconf,zstate,zorb,zpseud
      common /runlab/ zcom(19),ztitle(10),zaname(maxat),zbflab(maxorb),
     +   ztag(maxat),zsymm,zgroup,zscftp,zruntp,zguess,zconf,zstate,
     +   zorb(maxorb),zpseud(maxat)
c
c....  zaname : names of atoms during calculation (after reorder)
c....  ztag   : names of atoms as read in 
c
c
c...ecp parameters...
c
c
      integer ncoef1,ncoef2,j1len,j2len,llim,nlim,ntlim,j4len
      common /ecpdim/ ncoef1,ncoef2,j1len,j2len,llim,nlim,ntlim,j4len
c
c
      real*8 zetc,cax,cay,caz,ca,xca,yca,zca
      real*8 zetb,bax,bay,baz,ba,xba,yba,zba
      real*8 phase,dax,day,daz,da,xda,yda,zda,xint
      integer kcntr
      common /ecp1/ zetc,cax,cay,caz,ca,xca,yca,zca,
     +              zetb,bax,bay,baz,ba,xba,yba,zba,
     +              phase,dax,day,daz,da,xda,yda,zda,xint,
     +              kcntr
c
      common /ecp2  / clp(400),zlp(400),nlp(400),kfirst(maxat,6),
     *                klast(maxat,6),lmax(maxat),lpskip(maxat),
     *                izcore(maxat)
c
      logical iandj,norm,normi,normj
      common /ecpidx/ q2,iamin,iamax,jamin,jamax,ipmin,ipmax,jpmin,
     *                jpmax,kf1,kl1,llmx,npc,npb,iandj,norm,normi,normj
c
      real*8 bmcx,bmcy,bmcz,bpcx,bpcy,bpcz,cbsq,ax,ay,az
      logical candb
      common /ecp4/ bmcx,bmcy,bmcz,bpcx,bpcy,bpcz,cbsq,ax,ay,az,
     +              candb
c
c   ecp common stuff
      common /zfncm / x,y,z
c
c  the following are to make use of symmetry
c
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
c
      parameter (zero=0.0d+00)
c
c    iang gives the maximum angular momentum for each shell
c
      data iang/1,3*2,6*3,10*4,15*5 /
      data m18/18/
c
c     ----- ecp modifications to h integrals -----
c
c        -----  routine controls the calculation of the local    -----
c        -----  potential.  -ecp- program  by l. kahn.           -----
c        -----  gaussian basis integrals.                        -----
c        -----  written at battelle, oct. 1972.                  -----
c        -----  update in  1978.                                 -----
c        -----  rewritten at ims, okazaki japan  1980.           -----
c        -----  gradscf version august 1981.                     -----
c        -----  hondo version april 1984.                        -----
c        -----  gamess version october 1989.                     -----
c        -----  rewritten for gamess june 1995 bmb               -----
c
c
      norm=normf.ne.1.or.normp.ne.1
c  always want shells normalized with regular integrals
      normi=norm
      normj=norm
c
      nav = lenwrd()
c
c     -----  read in ecp formulas and data  -----
c
      ldaf91 = (j1len-1)/nav+1 + (9*ncoef1-1)/nav+1
      ldaf93 = (j2len-1)/nav+1 + (6*ncoef2)/nav
*****
      call secget(isect(479),m18,ibl479)
c
      call rdedx(dcoef1,ncoef1,ibl479,idaf)
      call reads(jfst1,ldaf91,idaf)
      call reads(dcoef2,ncoef2,idaf)
      call reads(jfst2,ldaf93,idaf)
c
      call dawt()
      call errt()
      call dawert()
      call ecpini(lmf,lmx,lmy,lmz)
      call ztab(zlm)
      call ftab(fpqr, nlim-1)
      call eccod3(fpqr,dcoef4,zlm,lmf,lmx,lmy,lmz)
      istart = 1
      iend   = nshell
      jstart = 1
      locij  = 0
      natst  = 1
      nated  = nat
      isave  = 0
c
c
c    initialize hecp---needed when symmetry is used
c
      ndum = (num*(num+1))/2
      call vclr(hecp,1,ndum)
c
c        -----  ishell  -----
c
      do 9010 ii=istart,iend
c
c     see if symmetry can be used to avoid this shell
c
         do it=1,nt
           id=mapshl(ii,it)
           if(id.gt.ii) go to 9010
           mi(it)=id
         enddo
c
       i1 = kstart(ii)
       i2 = i1+kng(ii)-1
       ipmin = i1
       ipmax = i2
       icntr = katom(ii)
       imin = kmin(ii)
       imax = kmax(ii)
       loci = kloc(ii)-imin-locij
       iimax = 1
       if (imin.eq.1.and.imax.eq.4) iimax = 2
         do 9000 iii=1,iimax
           if (imin.eq.1.and.imax.eq.4) then
             if (iii.eq.1) then
               iamin = 1
               iamax = 1
             else
               iamin = 2
               iamax = 4
             end if
           else
             iamin = imin
             iamax = imax
           end if
           npc=iang(iamax)
           igii=1
           do jj=ipmin,ipmax
             if (iamin.eq.1) then
               coefi(igii)=cs(jj)
             else if (iamin.lt.5) then
               coefi(igii)=cp(jj)
             else if (iamin.lt.11) then
               coefi(igii)=cd(jj)
             else if (iamin.le.20) then
               coefi(igii)=cf(jj)
             else if (iamin.le.35) then
               coefi(igii)=cg(jj)
             end if
             igii = igii + 1
           enddo
c
c        -----  jshell  -----
c
          if (oatomic) then
             if (ii.gt.1) then
                if (katom(ii).ne.katom(ii-1)) then
c
c                  the assumption is that all shells on
c                  a single atom are contiquous.
c
                   jstart = ii
                endif
             endif
          endif
          do 8010 jj=jstart,ii
           q2=1
c          check symmetry again
             n2=0
             do 80 it=1,nt
               jd=mapshl(jj,it)
               if(jd.gt.ii) go to 8010
               id=mi(it)
               if(id.ge.jd) go to 70
               nd=id
               id=jd
               jd=nd
 70            if(id.lt.ii) go to 80
               if(jd.lt.jj) go to 80
               if(jd.gt.jj) go to 8010
               n2=n2+1
 80          continue
             q2 = dble(nt)/dble(n2)
c
c
           iandj = ii.eq.jj
           j1 = kstart(jj)
           j2 = j1+kng(jj)-1
           jpmin = j1
           jpmax = j2
           jcntr = katom(jj)
           jmin = kmin(jj)
           jmax = kmax(jj)
           locj = kloc(jj)-jmin
           jjmax = 1
           if (jmin.eq.1.and.jmax.eq.4) jjmax = 2
           do 8000 jjj=1,jjmax
             if (jmin.eq.1.and.jmax.eq.4) then
               if (jjj.eq.1) then
                 jamin = 1
                 jamax = 1
               else
                 if (iandj.and.iamin.eq.1) go to 8000
                 jamin = 2
                 jamax = 4
               end if
             else
               jamin = jmin
               jamax = jmax
             end if
             jgjj=1
             do 112 icc=jpmin,jpmax
               if (jamin.eq.1) then
                 coefj(jgjj)=cs(icc)
               else if (jamin.lt.5) then
                 coefj(jgjj)=cp(icc)
               else if (jamin.lt.11) then
                 coefj(jgjj)=cd(icc)
               else if (jamin.le.20) then
                 coefj(jgjj)=cf(icc)
               else if (jamin.le.35) then
                 coefj(jgjj)=cg(icc)
               end if
               jgjj = jgjj + 1
 112         continue
             npb=iang(jamax)
c     n+n' is the sum of the angular momentum plus 1 to index arrays
             npnp = npc + npb - 1
             ijmax = max((iamax-iamin+1)*(jamax-jamin+1),60)
             do 100 i=1,ijmax
  100        g(i) = zero
             candb = icntr .eq. jcntr
             cx = c(1,icntr)
             cy = c(2,icntr)
             cz = c(3,icntr)
             bx = c(1,jcntr)
             by = c(2,jcntr)
             bz = c(3,jcntr)
             bmcx = bx-cx
             bmcy = by-cy
             bmcz = bz-cz
             bpcx = cx+bx
             bpcy = cy+by
             bpcz = cz+bz
             cbsq = bmcx*bmcx+bmcy*bmcy+bmcz*bmcz
c
c now loop over each center with an ecp potential
c
             if (oatomic) then
                natst = icntr
                nated = icntr
             endif
             do 7000 ikcntr=natst,nated
               kcntr = ikcntr
               ax = c(1,kcntr)
               ay = c(2,kcntr)
               az = c(3,kcntr)
               if(lpskip(kcntr).eq.1) go to 7000
               llmx = lmax(kcntr)+1
               kf1 = kfirst(kcntr,1)
               kl1 = klast(kcntr,1)
               canda = icntr .eq. kcntr
               aandb = kcntr .eq. jcntr
               if (canda) then
                 if (aandb) then
c   special case <a|a|a> only one center
                   iicab = 1
                 else
c   case <a|a|b>
                   cax = zero
                   cay = zero
                   caz = zero
                   ca = zero
                   bax = bx - ax
                   bay = by - ay
                   baz = bz - az
                   ba = sqrt(bax*bax+bay*bay+baz*baz)
                   x = bax/ba
                   y = bay/ba
                   z = baz/ba
                   iicab = 2
                   iipow = 1
                 end if
               else
                 if (aandb) then
c   case <c|a|a>
                   cax = cx - ax
                   cay = cy - ay
                   caz = cz - az
                   ca = sqrt(cax*cax+cay*cay+caz*caz)
                   x = cax/ca
                   y = cay/ca
                   z = caz/ca
                   bax = zero
                   bay = zero
                   baz = zero
                   ba = zero
                   iicab = 3
                   iipow = -1
                 else
c   general case <c|a|b> three-center integral
c   actually c and b may still be equal
                   cax = cx - ax
                   cay = cy - ay
                   caz = cz - az
                   ca = sqrt(cax*cax+cay*cay+caz*caz)
                   bax = bx - ax
                   bay = by - ay
                   baz = bz - az
                   ba = sqrt(bax*bax+bay*bay+baz*baz)
                   xca = cax/ca
                   yca = cay/ca
                   zca = caz/ca
                   xba = bax/ba
                   yba = bay/ba
                   zba = baz/ba
                   iicab = 4
                   iipow = 0
                 end if
               end if
c  set up tables of the powers of the cartesian distances (cax, cay ...)
c  for later use. pass in the maximum angular momentum for i and j
c  use max to index iang to make sure we get the max (l shells)
               if ((icntr.ne.kcntr).or.(kcntr.ne.jcntr)) then
                 call ecppwr(iipow,iang(iamax),iang(jamax))
               end if
               if (iicab.eq.1) then
                 call ecpa11(g,coefi,coefj,fpqr,npnp)
                 if (llmx.gt.1) then
                   call ecpa21(g,dcoef4,coefi,coefj,npnp)
                 end if
               else if (iicab.eq.2.or.iicab.eq.3) then
                 call ecpr12(fp,coefi,coefj,iicab,npnp)
                 call ecpa12(fp,jfst1,lbecp1,dcoef1,g,iicab,npnp,
     *              zlm,lmf,lmx,lmy,lmz)
                 if (llmx.gt.1) then
                   do 200 ll=2,llmx
                     call ecpr22(fp,coefi,coefj,iicab,npnp,ll)
                     call ecpa22(fp,jfst2,lbecp2,dcoef2,dcoef4,g,iicab,
     *                           npnp,ll-2,zlm,lmf,lmx,lmy,lmz)
 200               continue
                 end if
               else if (iicab.eq.4) then
                 call ecpr14(fp,coefi,coefj,npnp,zlm,lmf,lmx,lmy,lmz)
                 call ecpa14(fp,jfst1,lbecp1,dcoef1,g,npnp)
                 if (llmx.gt.1) then
                 ntempp=max(npc,npb)
                 do 300 ll=2,llmx
                   call ecpr24(fp,coefi,coefj,npnp,ntempp,ll)
                   call ecpa24(fp,jfst2,lbecp2,dcoef2,g,npnp,ll-2,
     *                zlm,lmf,lmx,lmy,lmz)
 300             continue
               end if
             end if
 7000      continue
c
c end of kcntr ecp potential center loop
c
c        -----  store ecp integrals in array.                    -----
c
           mmax = jamax
           nn = 0
           do 7500 i=iamin,iamax
             in = iky(loci+i)
             if(iandj) mmax = i
             do 7500 j=jamin,mmax
               jn = in+locj+j
               nn = nn+1
               hecp(jn) = g(nn) 
 7500        continue
c
c              end of both pairs of shell loops
c
 8000      continue
 8010    continue
 9000  continue
 9010 continue
c
c     -----  all done, modify bare nucleus hamiltonian -----
c
      if(dbug) then
         write(iwr,*) 'ecp integral modifications'
         call prtri(hecp,num)
      end if
c
        call vadd(h0,1,hecp,1,h0,1,l2)
c
        if(dbug) then
           write(iwr,*) '1e integrals, with ecp changes'
           call prtri(h0,num)
        end if
c
      return
      end
**==ecppwr.f
      subroutine ecppwr(iswtch, ni, nj)
      implicit real*8 (a-h,o-z)
c        -----  routine finds the   necessary powers of the     -----
c        -----  cartesian component vectors separating the      -----
c        -----  basis function and operator centers.            -----
c  this routine calculates (cax)**(nx-kx) where kx ranges from 0 to nx
c  thus we really need cax**(n) where n ranges from 0 to nx. these powers
c  are computed and stored in the acx... arrays with acx(n) being
c  acx**(n-1) power. ni and nj are the maximum angular momentum for the
c  i and j shell. here we precompute all required powers once for each
c  combination of shells.
c
c
      real*8 zetc,cax,cay,caz,ca,xca,yca,zca
      real*8 zetb,bax,bay,baz,ba,xba,yba,zba
      real*8 phase,dax,day,daz,da,xda,yda,zda,xint
      integer kcntr
      common /ecp1/ zetc,cax,cay,caz,ca,xca,yca,zca,
     +              zetb,bax,bay,baz,ba,xba,yba,zba,
     +              phase,dax,day,daz,da,xda,yda,zda,xint,
     +              kcntr
c
c       7 for g hessians
c
      real*8 acx,acy,acz,abx,aby,abz
      common /ecp3/ acx(7),acy(7),acz(7),abx(7),aby(7),abz(7)
c
      parameter (one=1.0d+00)
c first set the first element to one (0 power)
      acx(1) = one
      acy(1) = one
      acz(1) = one
      abx(1) = one
      aby(1) = one
      abz(1) = one
      if(iswtch.lt.1) then
        if (ni.gt.1) then
          acx(2) = -cax
          acy(2) = -cay
          acz(2) = -caz
          if (ni.gt.2) then
            do 50 ii=3,ni
              acx(ii) = acx(ii-1)*(-cax)
              acy(ii) = acy(ii-1)*(-cay)
              acz(ii) = acz(ii-1)*(-caz)
 50         continue
          end if
        end if
      end if
      if (iswtch.gt.(-1)) then
        if (nj.gt.1) then
          abx(2) = -bax
          aby(2) = -bay
          abz(2) = -baz
          if (nj.gt.2) then
            do 100 ii=3,nj
              abx(ii) = abx(ii-1)*(-bax)
              aby(ii) = aby(ii-1)*(-bay)
              abz(ii) = abz(ii-1)*(-baz)
 100        continue
          end if
        end if
      end if
      return
      end
**==ecpr12.f
      subroutine ecpr12(fp,coefi,coefj,icab,npnp)
c*   type 1 radial integrals for the given set of shells
c*   shell set types <b|a|a> or <a|a|c> (ie icab=(2 or 3))
c        fp will contain the radial integrals upon exiting this routine
c        npnp = n + n' indicates the maximum value for kappa needed
c
      implicit real*8 (a-h,o-z)
c
      dimension fp(*), fip(78), coefi(*), coefj(*)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
c
      real*8 zetc,cax,cay,caz,ca,xca,yca,zca
      real*8 zetb,bax,bay,baz,ba,xba,yba,zba
      real*8 phase,dax,day,daz,da,xda,yda,zda,xint
      integer kcntr
      common /ecp1/ zetc,cax,cay,caz,ca,xca,yca,zca,
     +              zetb,bax,bay,baz,ba,xba,yba,zba,
     +              phase,dax,day,daz,da,xda,yda,zda,xint,
     +              kcntr
c
      common /ecp2  / clp(400),zlp(400),nlp(400),kfirst(maxat,6),
     *                klast(maxat,6),lmax(maxat),lpskip(maxat),
     *                izcore(maxat)
c
      logical norm, normi, normj, iandj
      common /ecpidx/ q2,iamin,iamax,jamin,jamax,ipmin,ipmax,jpmin,
     *                jpmax,kf1,kl1,llmx,npc,npb,iandj,norm,normi,normj
      common /ficmn / alf,xi,xp0,xp1
c
      parameter(zero=0.0d+00, one=1.0d+00, two=2.0d+00)
c
c    maximum number of fp integrals needed for the given set of shells
c    fip must be dimensioned larger than max npnpmx (66 for g-hess)
c
      npnpmx = (npnp*(npnp+1))/2
c
c   first zero out the integral array
c
      do 100 np=1,npnpmx
 100  fp(np) = zero
c
c   loop over the primitives storing the radial integrals in fp
c
      igii=1
      do 700 ig = ipmin,ipmax
        zetc = ex(ig)
        dum = coefi(igii)
        igii=igii+1
        jgjj=1
        do 600 jg = jpmin,jpmax
          dum2 = dum*coefj(jgjj)
          jgjj=jgjj+1
          zetb = ex(jg)
          zetcb = zetc + zetb
          if (icab.eq.2) then
            alfa1 = zetb*ba
            xalfa1= alfa1*ba
          else if (icab.eq.3) then
            alfa1 = zetc*ca
            xalfa1 = alfa1*ca
          end if
          xp0 = exp(-xalfa1)
          alfi = one/(two*alfa1)
c   loop over the ecp primitives
          do 200 k=kf1,kl1
            xi=one/sqrt(zetcb+zlp(k))
            alf=alfa1*xi
            xp1=exp(-xalfa1+alf*alf)
            nlpk = nlp(k)
            clpk = clp(k)*dum2
c   calculate the actual radial integrals with fiecp
            call fiecp(fip, alfi, nlpk, npnp-1)
c   now store them in fp after multiplying by the coefs of each primitive
            do 190 n=1,npnpmx
              fp(n) = fp(n) + fip(n)*clpk
 190        continue
 200      continue
 600    continue
 700  continue
      return
      end
**==ecpr14.f
      subroutine ecpr14(fp,coefi,coefj,npnp,zlm,lmf,lmx,lmy,lmz)
c*   type 1 radial integrals for the given set of shells
c*   shell set types <b|a|c> (ie icab=4)
c
      implicit real*8 (a-h,o-z)
      dimension fp(*), fip(78), zfnlm(121), coefi(*),coefj(*)
      dimension zlm(*),lmf(*),lmx(*),lmy(*),lmz(*)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
      common /ecp2  / clp(400),zlp(400),nlp(400),kfirst(maxat,6),
     *                klast(maxat,6),lmax(maxat),lpskip(maxat),
     *                izcore(maxat)
      logical normi, normj, norm, iandj
      common /ecpidx/ q2,iamin,iamax,jamin,jamax,ipmin,ipmax,jpmin,
     *                jpmax,kf1,kl1,llmx,npc,npb,iandj,norm,normi,normj
c
      real*8 bmcx,bmcy,bmcz,bpcx,bpcy,bpcz,cbsq,ax,ay,az
      logical candb
      common /ecp4/ bmcx,bmcy,bmcz,bpcx,bpcy,bpcz,cbsq,ax,ay,az,
     +              candb
c
      common /ficmn / alf,xi,xp0,xp1
      common /zfncm / x,y,z
c
      parameter(zero=0.0d+00, one=1.0d+00, two=2.0d+00,half=0.5d+00)
      parameter(tol6=1.0d-06,onds4p=0.28209479177388d+00)
c
c    maximum number of fp integrals needed for the given set of shells
c
      npnpmx = (npnp*(npnp+1)*(2*npnp+1))/6
c
c   first zero out the integral array
c
      do 100 np=1,npnpmx
 100  fp(np) = zero
c
c   loop over the primitives storing the radial integrals in fp
c
      igii=1
      do 700 ig = ipmin,ipmax
        dum=coefi(igii)
        igii=igii+1
        zetc = ex(ig)
        jgjj=1
        do 600 jg = jpmin,jpmax
          dum2=dum*coefj(jgjj)
          jgjj=jgjj+1
          zetb = ex(jg)
          zetcb = zetc + zetb
c  transform the gaussians to center d actually c and b may be the
c  same as long as neither is the same as a. if c and b are different
c  then a phase factor is also required
          s = one/zetcb
          rat = (zetb-zetc)*s
          dax = half*(bpcx+rat*bmcx) - ax
          day = half*(bpcy+rat*bmcy) - ay
          daz = half*(bpcz+rat*bmcz) - az
          da = sqrt(dax*dax+day*day+daz*daz)
          if (.not.candb) dum2 = dum2*exp(-(zetc*zetb*s)*cbsq)
c  if center d and a both lie on the origin then all integrals for
c  lamda<>0 are 0 so special case that code
          if (da.ge.tol6) then
            alfa1 = zetcb*da
            xalfa1 = alfa1*da
            alfi1=one/(two*alfa1)
            alfi = alfi1
            x = dax/da
            y = day/da
            z = daz/da
            call zfn(zfnlm, npnp-1,zlm,lmf,lmx,lmy,lmz)
            xp0 = exp(-xalfa1)
            do 300 k=kf1,kl1
              xi=one/sqrt(zetcb+zlp(k))
              alf=alfa1*xi
              xp1=exp(-xalfa1+alf*alf)
              nlpk = nlp(k)
              clpk = clp(k)*dum2
              call fiecp(fip, alfi, nlpk, npnp-1)
              n=0
              nn=0
              do 200 kk=1,npnp
                do 200 l=1,kk
                  nn=nn+1
                  fiptem = fip(nn)*clpk
                  kkll = l*(l-1)+1+l
                  do 200 mu=1,2*l-1
                    n=n+1
                    fp(n) = fp(n) + fiptem*zfnlm(kkll-mu)
 200          continue
 300        continue
          else
c    special case when d=a=0. there are only l=0 integrals, but
c    they must be stored the same as above so index them into fp.
c    they also must be multiplied by the 0 spherical hamonic 1/sqrt(4pi)
            do 500 k=kf1,kl1
              zeta = zetcb+zlp(k)
              nlpk = nlp(k)
              clpk = clp(k)*dum2*onds4p
              do 400 n=1,npnp
                nk = (n*(n-1)*(2*n-1))/6 + 1
                fp(nk) = fp(nk)+fa(n+nlpk-1, zeta)*clpk
 400          continue
 500        continue
          end if
 600    continue
 700  continue
      end
**==ecpr22.f
      subroutine ecpr22(fp,coefi,coefj,icab,npnp,ll)
c*   type 2 radial integrals for the given set of shells
c*   shell set types <b|a|a> or <a|a|c> (ie icab=(2 or 3))
c        fp will contain the radial integrals upon exiting this routine
c        npnp = n + n' indicates the maximum value for kappa needed
c        ll gives the value of l for the <|ulmax-ul|> integral
c
      implicit real*8 (a-h,o-z)
      dimension fp(*), fip(78), coefi(*),coefj(*)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
c
      real*8 zetc,cax,cay,caz,ca,xca,yca,zca
      real*8 zetb,bax,bay,baz,ba,xba,yba,zba
      real*8 phase,dax,day,daz,da,xda,yda,zda,xint
      integer kcntr
      common /ecp1/ zetc,cax,cay,caz,ca,xca,yca,zca,
     +              zetb,bax,bay,baz,ba,xba,yba,zba,
     +              phase,dax,day,daz,da,xda,yda,zda,xint,
     +              kcntr
c
      common /ecp2  / clp(400),zlp(400),nlp(400),kfirst(maxat,6),
     *                klast(maxat,6),lmax(maxat),lpskip(maxat),
     *                izcore(maxat)
c
      logical norm, normi, normj, iandj
      common /ecpidx/ q2,iamin,iamax,jamin,jamax,ipmin,ipmax,jpmin,
     *                jpmax,kf1,kl1,llmx,npc,npb,iandj,norm,normi,normj
      common /ficmn / alf,xi,xp0,xp1
c
      parameter(zero=0.0d+00, one=1.0d+00, two=2.0d+00)
c
c    maximum number of fp integrals needed for the given set of shells
c
      npnpmx = (npnp*(npnp+1))/2
c set up the indeces based on ll
      kf = kfirst(kcntr,ll)
      kl = klast(kcntr,ll)
c
c   first zero out the integral array
c
      do 100 n=1,npnpmx
 100  fp(n) = zero
c
c   loop over the primitives storing the radial integrals in fp
c
      igii=1
      do 700 ig = ipmin,ipmax
        dum=coefi(igii)
        igii=igii+1
        zetc = ex(ig)
        jgjj=1
        do 600 jg = jpmin,jpmax
          dum2=dum*coefj(jgjj)
          jgjj=jgjj+1
          zetb = ex(jg)
          zetcb = zetc + zetb
          if (icab.eq.2) then
            gamma = zetb*ba
            xgamma = gamma*ba
          else if (icab.eq.3) then
            gamma = zetc*ca
            xgamma = gamma*ca
          end if
          gammi = one/(two*gamma)
          xp0 = exp(-xgamma)
c     loop over the ecp primitives
          do 200 k=kf,kl
            xi=one/sqrt(zetcb+zlp(k))
            alf=gamma*xi
            xp1=exp(-xgamma+alf*alf)
            nlpk = nlp(k)
            clpk = clp(k)*dum2
            call fiecp(fip, gammi, nlpk, npnp-1)
            do 190 n=1,npnpmx
              fp(n) = fp(n) + fip(n)*clpk
 190        continue
 200      continue
 600    continue
 700  continue
      return
      end
**==ecpr24.f
      subroutine ecpr24(fp,coefi,coefj,npnp,np1,ll)
c
c*   type 2 radial integrals for the given set of shells
c*   shell set types <b|a|c> (ie icab=4)
c
      implicit real*8 (a-h,o-z)
c since fp is compacted it needs to be dimensioned 11*11*6
      dimension fp(*), coefi(*),coefj(*)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
c
      real*8 zetc,cax,cay,caz,ca,xca,yca,zca
      real*8 zetb,bax,bay,baz,ba,xba,yba,zba
      real*8 phase,dax,day,daz,da,xda,yda,zda,xint
      integer kcntr
      common /ecp1/ zetc,cax,cay,caz,ca,xca,yca,zca,
     +              zetb,bax,bay,baz,ba,xba,yba,zba,
     +              phase,dax,day,daz,da,xda,yda,zda,xint,
     +              kcntr
c
      common /ecp2  / clp(400),zlp(400),nlp(400),kfirst(maxat,6),
     *                klast(maxat,6),lmax(maxat),lpskip(maxat),
     *                izcore(maxat)
      logical normi, normj, norm, iandj
      common /ecpidx/ q2,iamin,iamax,jamin,jamax,ipmin,ipmax,jpmin,
     *                jpmax,kf1,kl1,llmx,npc,npb,iandj,norm,normi,normj
      common /fjcmn / alef,beit,xxi,xpls,xmns,xp
      common /fjnew / xka, xkb, gamma1,gamma2,a1,a2,c
c
      parameter(zero=0.0d+00, one=1.0d+00, two=2.0d+00)
      parameter(ablim=1.0d-01)
c
c  local storage
      dimension fjpq(11,11,11)
      save fjpq
c
c    setup various parameters for loops
c
      l=ll-2
      np1 = max(npc,npb)
      npcpl=npc+l
      npbpl=npb+l
      nmax = 6+(npnp-1)*npcpl*(npbpl/2 +1)
      np1pl = np1+l
      ltmax = max(npnp,np1pl)
      lemx = max(l,ltmax/2)
c set up the indeces based on ll
      kf = kfirst(kcntr,ll)
      kl = klast(kcntr,ll)
c
c   first zero out the integral array
c
      do 100 n=1,nmax
        fp(n) = zero
 100  continue
c
c   loop over the primitives storing the radial integrals in fp
c
      igii=1
      do 700 ig = ipmin,ipmax
        dum = coefi(igii)
        igii=igii+1
        zetc = ex(ig)
        alfa2 = zetc*ca
        xalfa = alfa2*ca
        alfi = one/(two*alfa2)
        jgjj=1
        do 600 jg = jpmin,jpmax
          dum2 = dum*coefj(jgjj)
          jgjj=jgjj+1
          zetb = ex(jg)
          zetcb = zetc + zetb
          beta = zetb*ba
          xbeta = beta*ba
          beti = one/(two*beta)
          alfbet = xalfa+xbeta
          xp=exp(-alfbet)
c
          do 300 k=kf,kl
            xxi=one/sqrt(zetcb+zlp(k))
            alef=alfa2*xxi
            beit = beta*xxi
            if(alef*beit.gt.ablim) then
              xka = 2*alfa2
              xkb = 2*beta
              a1 = alef+beit
              a2 = alef-beit
              c = zetcb+zlp(k)
              gamma1 = 0.25d+00*exp(-alfbet+(alef+beit)**2)
              gamma2 = 0.25d+00*exp(-alfbet+(alef-beit)**2)
              xpls = exp(-alfbet+(alef+beit)**2)
              xmns = exp(-alfbet+(alef-beit)**2)
            else
              xpls = exp(-alfbet+alef*alef)
              xmns = exp(-alfbet+beit*beit)
            end if
            nlpk = nlp(k)
            clpk = clp(k)*dum2
            call fjecp(fjpq,alfi,beti, nlpk, npnp, ltmax,lemx)
c  the first block is special since there are only values on the
c  diagonal so stick into linear memory (max of 6 elements)
            do 170 in=1,lemx+1
              fp(in) = fp(in) + fjpq(in,in,1)*clpk
 170        continue
            if (npnp.gt.1) then
c  do higher blocks here, note: every other element is zero and is
c  therefore zero (actually the angular term will be zero)
              nstart = 1
              n=6
              do 200 in=2,npnp
                nstart = 1-nstart
                lstart = nstart
                do 190 ip=1,npcpl
                  lstart = 1-lstart
                  do 180 iq=1+lstart,npbpl,2
                    n=n+1
                    fp(n) = fp(n) + fjpq(iq,ip,in)*clpk
 180              continue
 190            continue
 200          continue
            end if
 300      continue
 600    continue
 700  continue
      end
**==errf.f
      function errf(y)
      implicit real*8 (a-h,o-z)
c
c        -----  routine evaluates the error function.           -----
c        -----  the error function is an odd-parity function.   -----
c        -----  evaluate   by a piecewise chebyshev polynomial  -----
c        -----  if needed use the asymptotic value.             -----
c
      common /errfcm/ c(142),ifirst(20),ilast(20),h
      x = y
      sign = +1.0d+00
      if(x.ge.0.0d+00)go to 10
      x = -y
      sign = -1.0d+00
   10 if(x.ge.4.86d+00)go to 30
      xn = x/h
      nx = int(xn)
      nx = nx+1
      if = ifirst(nx)
      il = ilast(nx)
      t = c(il)
      kl = il-if
      do 20 k = 1,kl
   20 t = c(il-k)+x*t
      errf = t*sign
      return
   30 errf = sign
      return
      end
**==errt.f
      subroutine errt
      implicit real*8 (a-h,o-z)
      common /errfcm/ c(142),ifirst(20),ilast(20),h
c
c        -----  routine allocates parameters for the polynomial -----
c        -----  fit to the error function.                      -----
c
      h = 0.25d+00
c... 0.00 @ x @ .25 interval no. 1 abs.error = 5.1496584774213d-
      ifirst( 1)=1
      ilast ( 1)=8
      c( 1) =  0.9976905682968d-15
      c( 2) =  0.1128379167616d+01
      c( 3) = -0.6825080083050d-07
      c( 4) = -0.3761237879247d+00
      c( 5) = -0.4362679036164d-04
      c( 6) =  0.1132081528735d+00
      c( 7) = -0.1608829014003d-02
      c( 8) = -0.2383745887450d-01
c...  .25 @ x @  .50  interval no. 2  abs.error= 1.0329515021112d-
      ifirst( 2)=9
      ilast ( 2)=16
      c( 9)  = -0.3265133599788d-05
      c( 10) =  0.1128453123654d+01
      c( 11) = -0.7236769295451d-03
      c( 12) = -0.3721384745633d+00
      c( 13) = -0.1348362056472d-01
      c(14) =   0.1412226404729d+00
      c(15) = -0.3539596498013d-01
      c(16) = -0.5454029355730d-02
c... .50 @ x @ .75 interval no. 3 abs.error=6.1888272284705d-
      ifirst( 3)=17
      ilast ( 3)=24
      c(17) = -0.1801988094647d-03
      c(18) =  0.1130782166249d+01
      c(19) = -0.1400474656066d-01
      c(20) = -0.3295493300819d+00
      c(21) = -0.9653664285517d-01
      c(22) =  0.2398163462058d+00
      c(23) = -0.1014115256923d+00
      c(24) =  0.1378314835685d-01
c... .75 @ x @ 1.00 interval no. 4 abs.error=1.6626700016786d-
      ifirst( 4)=25
      ilast ( 4)=32
      c(25) = -0.4457991000617d-03
      c(26) =  0.1133670008174d+01
      c(27) = -0.2725712863074d-01
      c(28) = -0.2961333291272d+00
      c(29) = -0.1466850258377d+00
      c(30) =  0.2847021965842d+00
      c(31) = -0.1236317583493d+00
      c(32) =  0.1848162923540d-01
c... 1.00 @ x @ 1.25  interval no. 5  abs.error=5.4747317790316d-
      ifirst( 5)=33
      ilast ( 5)=40
      c(33) =  0.1036667435991d-01
      c(34) =  0.1060681450430d+01
      c(35) =  0.1844323965440d+00
      c(36) = -0.6381386282176d+00
      c(37) =  0.1857801550899d+00
      c(38) =  0.9020967608584d-01
      c(39) = -0.6022465654782d-01
      c(40) =  0.9593725204468d-02
c... 1.25 @ x @ 1.50 interval no. 6 abs.error=2.1689317009077d-
      ifirst( 6)=41
      ilast ( 6)=47
      c(41) =  0.5096577645657d-01
      c(42) =  0.8301050352144d+00
      c(43) =  0.7461454169558d+00
      c(44) = -0.1399088183681d+01
      c(45) =  0.8049290742492d+00
      c(46) = -0.2123844940215d+00
      c(47) =  0.2202851573626d-01
c... 1.50 @ x @ 1.75 interval no. 7 abs.error=5.4356519285648d-
      ifirst( 7)=48
      ilast ( 7)=55
      c(48) =  0.1356983343156d+00
      c(49) =  0.4248183601853d+00
      c(50) =  0.1578617255285d+01
      c(51) = -0.2350907279473d+01
      c(52) =  0.1459137006763d+01
      c(53) = -0.4826721286933d+00
      c(54) =  0.8417719176837d-01
      c(55) = -0.6134748458862d-02
c... 1.75 @ x @ 2.00  interval no. 8  abs.error=1.3677947663382d-
      ifirst( 8)=56
      ilast ( 8)=63
      c(56) =  0.5650617287234d-01
      c(57) =  0.7303727164846d+00
      c(58) =  0.1073427393133d+01
      c(59) = -0.1886956498357d+01
      c(60) =  0.1203543744299d+01
      c(61) = -0.3982083398317d+00
      c(62) =  0.6867500288146d-01
      c(63) = -0.4915782383510d-02
c... 2.00 @ x @ 2.25  interval no. 9  abs.error =2.7604585284280d-
      ifirst( 9)=64
      ilast ( 9)=70
      c(64) = -0.7044612080369d+00
      c(65) =  0.3313423632378d+01
      c(66) = -0.2688303271381d+01
      c(67) =  0.1159938873721d+01
      c(68) = -0.2789363703341d+00
      c(69) =  0.3510668004553d-01
      c(70) = -0.1778023938338d-02
c... 2.25 @ x @ 2.50 interval no. 10 abs.error=1.8687273950491d-
      ifirst(10)=71
      ilast (10)=77
      c(71) = -0.8530571392369d+00
      c(72) =  0.3718753089373d+01
      c(73) = -0.3149110592528d+01
      c(74) =  0.1439416506820d+01
      c(75) = -0.3743063133443d+00
      c(76) =  0.5246810708195d-01
      c(77) = -0.3095234433810d-02
c... 2.50 @ x @ 2.75  interval no. 11  abs.error=9.7095664841618d-
      ifirst(11)=78
      ilast (11)=84
      c(78) = -0.6325006110506d+00
      c(79) =  0.3195352274705d+01
      c(80) = -0.2631359977752d+01
      c(81) =  0.1166146516443d+01
      c(82) = -0.2931400812425d+00
      c(83) =  0.3960487060249d-01
      c(84) = -0.2245448529720d-02
c... 2.75 @ x @ 3.00  interval no. 12  abs.error=7.6525452641363d-
      ifirst(12)=85
      ilast (12)=91
      c(85) = -0.1578663656307d+00
      c(86) =  0.2160403160352d+01
      c(87) = -0.1690605397678d+01
      c(88) =  0.7098560625163d+00
      c(89) = -0.1685917225356d+00
      c(90) =  0.2146466169506d-01
      c(91) = -0.1144051551819d-02
c... 3.00 @ x @ 3.25  interval no. 13  abs.error =3.9079850466806d-
      ifirst(13)=92
      ilast (13)=98
      c(92) =  0.3337223693570d+00
      c(93) =  0.1174814163831d+01
      c(94) = -0.8669311383990d+00
      c(95) =  0.3425820853154d+00
      c(96) = -0.7643589524863d-01
      c(97) =  0.9127119866510d-02
      c(98) = -0.4555632670720d-03
c... 3.25 @ x @ 3.50  interval no. 14  abs.error=1.5241141682054d-
      ifirst(14)=99
      ilast (14) =105
      c( 99) = 0.6858505344288d+00
      c(100) = 0.5221032227938d+00
      c(101) = -0.3626430750192d+00
      c(102) = 0.1347162453894d+00
      c(103) = -0.2822370117065d-01
      c(104) = 0.3161233228942d-02
      c(105) =-0.1478642225266d-03
c... 3.50 @ x @ 3.75  interval no. 15  abs.error=9.5496943686157d-
      ifirst(15) =106
      ilast (15) =111
      c(106) =  0.9680982760858d+00
      c(107) =  0.4187846825891d-01
      c(108) = -0.2202918541743d-01
      c(109) =  0.5803572496370d-02
      c(110) = -0.7656564703211d-03
      c(111) =  0.4046317189932d-04
c... 3.75 @ x @ 4.00  interval no. 16  abs.error=2.2239987629291d-
      ifirst(16) =112
      ilast (16) =117
      c(112) = 0.9910262142894d+00
      c(113) = 0.1110389476297d-01
      c(114) =-0.5503022341054d-02
      c(115) = 0.1365301693295d-02
      c(116) =-0.1695606159046d-03
      c(117) = 0.8432380855083d-05
c... 4.00 @ x @ 4.25  interval no. 17  abs.error=4.4764192352886d-
      ifirst(17) =118
      ilast (17) =123
      c(118) = 0.9978511542662d+00
      c(119) = 0.2512562013763d-02
      c(120) =-0.1176282812867d-02
      c(121) = 0.2755980203801d-03
      c(122) =-0.3231375012547d-04
      c(123) = 0.1516751945019d-05
c... 4.25 @ x @ 4.50  interval no. 18  abs.error=9.9475983006414d-
      ifirst(18) =124
      ilast (18) =129
      c(124) = 0.9995577955752d+00
      c(125) = 0.4898449770963d-03
      c(126) =-0.2172069044718d-03
      c(127) = 0.4819046007469d-04
      c(128) =-0.5349400453269d-05
      c(129) = 0.2376735210419d-06
c... 4.50 @ x @ 4.75  interval no. 19  abs.error=6.0396132539609d-
      ifirst(19) =130
      ilast (19) =135
      c(130) =  0.9999195780401d+00
      c(131) =  0.8461205943888d-04
      c(132) = -0.3562696920198d-04
      c(133) =  0.7504334644182d-05
      c(134) = -0.7907161489129d-06
      c(135) =  0.3334134817123d-07
c... 4.75 @ x @ 5.00  interval no. 20  abs.error=7.4429351570870d-
      ifirst(20) =136
      ilast (20) =142
      c(136) =  0.8965357422439d+00
      c(137) =  0.1272754720567d+00
      c(138) = -0.6523024655166d-01
      c(139) =  0.1782843628219d-01
      c(140) = -0.2740698420287d-02
      c(141) =  0.2246825024486d-03
      c(142) = -0.7674098014832d-05
      return
      end
**==fa.f
      function fa(n, zeta)
      implicit real*8 (a-h,o-z)
c        routine evaluates the type 1 radial integral for the special
c        case of k=-2(zeta*ca+zetab*cb)=0. thus, q(n,lamda,kr) reduces
c        to integral from zero to infinity of [r**n e**(-zeta*r**2)]
c        good through for n=0 - n=17
      dimension gammo(9),gamme(9)
      data gammo/0.5d+00,0.5d+00,1.0d+00,3.0d+00,12.0d+00,60.0d+00,
     *           360.0d+00,2520.0d+00,20160.0d+00/
      data gamme/0.5d+00,0.25d+00,0.375d+00,0.9375d+00,3.28125d+00,
     *           14.765625d+00,81.2109375d+00,527.87109375d+00,
     *           3959.033203125d+00/
      data sqrpi/1.772453850905d+00/
c n = odd terms
      if (mod(n,2).ne.0) then
        k = (n-1)/2
        fa= gammo(k+1)/zeta**(k+1)
      else
c n = even terms
        k = n/2
        fa= gamme(k+1)*sqrpi/((zeta**k)*sqrt(zeta))
      end if
      return
      end
**==fiecp.f
      subroutine fiecp(fip, alfi, nlp, npnp)
      implicit real*8 (a-h,o-z)
      common /ficmn / alf,xi,xp0,xp1
      dimension fit(19,2), fip(*)
      parameter(half=0.5d+00, one=1.0d+00)
c
c        -----  routine prepares table of i- integrals by use   -----
c        -----  of a recursion relationship and then it picks   -----
c        -----  out among these the currently required ones.    -----
c        -----  in this manner one calculates only the unique   -----
c        -----  minmal number of i-integrals.                   -----
c
c        -----  find the largest n value occuring in both the   -----
c        -----  even and odd lattices.                          -----
      neff = npnp + nlp
      if (mod(neff,2).eq.1) then
        nemx = neff + 1
        nomx = neff
      else
        nemx = neff
        nomx = neff+1
      end if
      if (neff.eq.0) nomx=-1
c        -----  prepare variables.                              -----
      yi = half*xi
      a1 = xi*alf
      a2 = xi*yi
c        -----  generate a table of elements of even lattice    -----
      if(nemx.ge.0) then
        fit(1,1) = fsi0(0)*yi
        if(nemx.ge.1) then
          fit(2,2) = fsi1(1)*yi*yi
          if(nemx.ge.2) then
            do 40 n=2,nemx,2
              fnm1 = n-1
              fit(n+1,1) = a1*fit(n,2)+fnm1*a2*fit(n-1,1)
              fit(n+2,2) = a1*fit(n+1,1)+(fnm1-one)*a2*fit(n,2)
   40       continue
          end if
        end if
      end if
c        -----  generate table of elements of odd-lattice       -----
      if(nomx.ge.0) then
        fit(1,2) = fsi1(0)*yi
        if(nomx.ge.1) then
          fit(2,1) = fsi0(1)*yi*yi
          if(nomx.ge.2) then
            do 60 n=2,nomx,2
              fn = n
              fit(n+1,2) = a1*fit(n,1)+(fn-3.0d+00)*a2*fit(n-1,2)
              fit(n+2,1) = a1*fit(n+1,2)+fn*a2*fit(n,1)
   60       continue
          end if
        end if
      end if
c        -----  retrieve required integs from tables.           -----
      nk=0
      do 90 n=1,npnp+1
        nk = (n*(n-1))/2
        if (n.eq.1) then
          lmax = 1
        else
          lmax = 2
        end if
        do 70 l=1,lmax
          fip(nk+l) = fit(n+nlp,l)
   70   continue
c   apply modified spherical bessal function recursion relation to
c   get integrals for all higher terms
        if (n.gt.2) then
          nm1k=((n-1)*(n-2))/2
          do 80 l=3,n
            fip(nk+l) = fip(nk+l-2)-(2*l-3)*alfi*fip(nm1k+l-1)
   80     continue
        end if
   90 continue
      return
      end
**==fjecp,f
      subroutine fjecp(fjpq,alfi,beti,nlp, npnp, lmax,lemx)
      implicit real*8 (a-h,o-z)
c        -----  routine prepares table of j-integs by use of    -----
c        -----  a recursion relationship and then it picks up   -----
c        -----  from among these the ones that are currently    -----
c        -----  needed. thus
c        -----  one only calculates the unique number.          -----
c
      common /fjcmn / alef,bet,xi,xpls,xmns,xp
      dimension fjpq(11,11,11),fjt(19,2,2)
      dimension rho(15), sigma(15), sigmab(15), tau(15)
      common /fjnew / xka, xkb, gamma1,gamma2,a1,a2,c
      parameter (ablim=1.0d-01, zero=0.0d+00)
c
c        -----  find the largest n values occuring in both the  -----
c        -----  even and odd lattices.                          -----
c
      jeff = npnp - 1 + nlp
      if(mod(jeff,2).eq.1) then
        nemx = jeff
        nomx = jeff
      else
        nemx = jeff
        nomx = jeff
      end if
      if (jeff.eq.0) nomx = -1
      lomx = lemx
c        -----  if the j-integrals are going to be calcualted   -----
c        -----  by a power series, then prepare a table of      -----
c        -----  i integs.                                       -----
      if(alef*bet.gt.ablim) then
        xkakb=xka*xkb
        onedab = 1.0d+00 / xkakb
        ondab2 = onedab*onedab
        if (2*lemx+3 .gt.15) then
          write(6,*) 'lemx too big, lemx= ',lemx
          go to 300
        end if
        do 10 ii=1,15
          rho(ii) = zero
          sigma(ii) = zero
          sigmab(ii) = zero
          tau(ii) = zero
   10   continue
        call fjform(rho,sigma,sigmab,tau, 2*lemx+3)
        fjt(1,1,1) = rho(3)*onedab
        fjt(2,1,1) = rho(2)*onedab
        fjt(1,2,1) = sigmab(3)*onedab - rho(4)*onedab/xka
        fjt(2,2,1) = sigmab(2)*onedab - rho(3)*onedab/xka
        fjt(1,1,2) = sigma(3)*onedab - rho(4)*onedab/xkb
        fjt(2,1,2) = sigma(2)*onedab - rho(3)*onedab/xkb
        fjt(1,2,2) = ondab2*(rho(5)-xka*sigmab(4)-xkb*sigma(4)+
     *               xkakb*tau(3))
        fjt(2,2,2) = ondab2*(rho(4)-xka*sigmab(3)-xkb*sigma(3)+
     *               xkakb*tau(2))
      else
        call sitabl(lemx,lomx)
        fjt(1,1,1) = fjps(0,0,0)
        fjt(2,1,1) = fjps(1,0,0)
        fjt(1,2,1) = fjps(0,1,0)
        fjt(2,2,1) = fjps(1,1,0)
        fjt(1,1,2) = fjps(0,0,1)
        fjt(2,1,2) = fjps(1,0,1)
        fjt(1,2,2) = fjps(0,1,1)
        fjt(2,2,2) = fjps(1,1,1)
      end if
      xa = xi*alef
      xb = xi*bet
      xx = xi*xi*0.5d+00
c        -----  generate a table of elements of the even        -----
c        -----  lattice of j-integs.                            -----
      if(nemx.ge.0) then
        if(nemx.ge.1) then
          if(nemx.ge.2) then
            do 40 n=2,nemx,2
              np1 = n+1
              fnm2 = n-2
              fjt(np1,1,1) = xa*fjt(n,2,1)+xb*fjt(n,1,2)+
     *               (fnm2+1.0d+00)*xx*fjt(n-1,1,1)
              fjt(np1,2,2) = xa*fjt(n,1,2)+xb*fjt(n,2,1)+
     *               (fnm2-3.0d+00)*xx*fjt(n-1,2,2)
              fjt(n+2,2,1) = xa*fjt(np1,1,1)+xb*fjt(np1,2,2)+
     *                fnm2*xx*fjt(n,2,1)
              fjt(n+2,1,2) = xb*fjt(np1,1,1)+xa*fjt(np1,2,2)+
     *                fnm2*xx*fjt(n,1,2)
   40       continue
          end if
        end if
      end if
c        -----  generate table of elements of the odd lattice   -----
c        -----  j-integs.                                       -----
      if(nomx.ge.0) then
        if(nomx.ge.1) then
          if(nomx.ge.2) then
            do 60 n=2,nomx,2
              np1 = n+1
              fnm3 = n-3
              fjt(np1,2,1) = xa*fjt(n,1,1)+xb*fjt(n,2,2)+fnm3*xx*
     *                       fjt(n-1,2,1)
              fjt(np1,1,2) = xb*fjt(n,1,1)+xa*fjt(n,2,2)+fnm3*xx*
     *                       fjt(n-1,1,2)
              fjt(n+2,1,1) = xa*fjt(np1,2,1)+xb*fjt(np1,1,2)+
     *                       (fnm3+3.0d+00)*xx*fjt(n,1,1)
              fjt(n+2,2,2) = xa*fjt(np1,1,2)+xb*fjt(np1,2,1)+
     *                       (fnm3-1.0d+00)*xx*fjt(n,2,2)
   60       continue
          end if
        end if
      end if
c        -----  retrieve the required j-integs from tables.     -----
c  first take care of kappa=0 block
      do 150 ip = 1,lemx+1
        if (ip.lt.3) then
          fjpq(ip,ip,1) = fjt(nlp+1,ip,ip)
        else if (ip.eq.3) then
          if (alef*bet.gt.ablim) then
            if (nlp.gt.2) then
              write(6,*) 'nlp too big!, ',nlp
              go to 300
            end if
            xka2 = xka*xka
            xkb2 = xkb*xkb
            xkakb2= xkakb*xkakb
            fjpq(ip,ip,1) = ondab2*onedab*(9.0d+00*rho(7-nlp)+3.0d+00*
     *        rho(5-nlp)*(xka2+xkb2)+xkakb2*rho(3-nlp)-3.0d+00*xka*
     *        xkakb*sigma(4-nlp)-9.0d+00*xkb*sigma(6-nlp)-3.0d+00*
     *        xkakb*xkb*sigmab(4-nlp)-9.0d+00*xka*sigmab(6-nlp)+
     *        9.0d+00*xkakb*tau(5-nlp))
          else
            fjpq(ip,ip,1) = fjps(nlp,2,2)
          end if
        else if (ip.eq.4) then
          if (alef*bet.gt.ablim) then
            if (nlp.gt.2) then
              write(6,*) 'nlp too big!, ',nlp
              go to 300
            end if
            fjpq(ip,ip,1) = ondab2*ondab2*(2.25d+02*rho(9-nlp)+
     *          9.0d+01*(xka2+xkb2)*rho(7-nlp)+3.6d+01*xkakb2*
     *          rho(5-nlp)-(2.25d+02*xkb*sigma(8-nlp)+xkb*(1.5d+01*
     *          xkb2+9.0d+01*xka2)*sigma(6-nlp)+6.0d+00*xkakb2*xkb*
     *          sigma(4-nlp))-(2.25d+02*xka*sigmab(8-nlp)+xka*(1.5d+01*
     *          xka2+9.0d+01*xkb2)*sigmab(6-nlp)+6.0d+00*xka*xkakb2*
     *          sigmab(4-nlp))+2.25d+02*xkakb*tau(7-nlp)+1.5d+01*xkakb*
     *          (xka2+xkb2)*tau(5-nlp)+xkakb2*xkakb*tau(3-nlp))
          else
            fjpq(ip,ip,1) = fjps(nlp,3,3)
          end if
        else if (ip.eq.5) then
          if (alef*bet.gt.ablim) then
            if (nlp.gt.2) then
              write(6,*) 'nlp too big!, ',nlp,' ip=5'
              go to 300
            end if
            xka4 = xka2*xka2
            xkb4 = xkb2*xkb2
            fjpq(ip,ip,1) = (1.1025d+04*rho(11-nlp)+4.725d+03*(xka2+
     *          xkb2)*rho(9-nlp)+(1.05d+02*(xka4+xkb4)+2.025d+03*
     *          xkakb2)*rho(7-nlp)+4.5d+01*xkakb2*(xka2+xkb2)*
     *          rho(5-nlp)+xkakb2*xkakb2*rho(3-nlp)-(xkb*(1.1025d+04*
     *          sigma(10-nlp)+(1.05d+03*xkb2+4.725d+03*xka2)*
     *          sigma(8-nlp)+(4.5d+02*xkakb2+1.05d+02*xka4)*
     *          sigma(6-nlp)+1.0d+01*xka2*xkakb2*sigma(4-nlp))+xka*
     *          (1.1025d+04*sigmab(10-nlp)+(1.05d+03*xka2+4.725d+03*
     *          xkb2)*sigmab(8-nlp)+(4.5d+02*xkakb2+1.05d+02*xkb4)*
     *          sigmab(6-nlp)+1.0d+01*xkakb2*xkb2*sigmab(4-nlp)))+
     *          xkakb*(1.1025d+04*tau(9-nlp)+1.05d+03*(xka2+xkb2)*
     *          tau(7-nlp)+1.0d+02*xkakb2*tau(5-nlp)))*(ondab2*ondab2*
     *          onedab)
          else
            fjpq(ip,ip,1) = fjps(nlp,4,4)
          end if
        else if (ip.eq.6) then
          if (alef*bet.gt.ablim) then
            if (nlp.gt.2) then
              write(6,*) 'nlp too big!, ',nlp,' ip=6'
              go to 300
            end if
            fjpq(ip,ip,1) = (8.93025d+05*rho(13-nlp)+3.969d+05*(xka2+
     *        xkb2)*rho(11-nlp)+(1.4175d+04*(xka4+xkb4)+1.764d+05*
     *        xkakb2)*rho(9-nlp)+6.3d+03*xkakb2*(xka2+xkb2)*rho(7-nlp)+
     *        2.25d+02*xkakb2*xkakb2*rho(5-nlp)-(xkb*(8.93025d+05*
     *        sigma(12-nlp)+(3.969d+05*xka2+9.9225d+04*xkb2)*
     *        sigma(10-nlp)+(1.4175d+04*xka4+4.41d+04*xkakb2+9.45d+02*
     *        xkb4)*sigma(8-nlp)+(1.575d+03*xkakb2*xka2+4.2d+02*xkakb2*
     *        xkb2)*sigma(6-nlp)+1.5d+01*xkakb2*xkakb2*sigma(4-nlp))+
     *        xka*(8.93025d+05*sigmab(12-nlp)+(3.969d+05*xkb2+
     *        9.9225d+04*xka2)*sigmab(10-nlp)+(1.4175d+04*xkb4+4.41d+04
     *        *xkakb2+9.45d+02*xka4)*sigmab(8-nlp)+(1.575d+03*xkakb2*
     *        xkb2+4.2d+02*xka2*xkakb2)*sigmab(6-nlp)+1.5d+01*xkakb2*
     *        xkakb2*sigmab(4-nlp)))+xkakb*(8.93025d+05*tau(11-nlp)+
     *        9.9225d+04*(xka2+xkb2)*tau(9-nlp)+(9.45d+02*(xka4+xkb4)+
     *        1.1025d+04*xkakb2)*tau(7-nlp)+1.05d+02*xkakb2*(xka2+xkb2)
     *        *tau(5-nlp)+xkakb2*xkakb2*tau(3-nlp)))*(ondab2*ondab2*
     *        ondab2)
          else
            fjpq(ip,ip,1) = fjps(nlp,5,5)
          end if
        else if (ip.eq.7) then
          if (alef*bet.gt.ablim) then
            if (nlp.gt.2) then
              write(6,*) 'nlp too big!, ',nlp,' ip=7'
              go to 300
            end if
            dum = 1.08056025d+08*rho(15-nlp)+4.9116375d+07*
     *       (xka2+xkb2)*rho(13-nlp)+(2.18295d+06*(xka4+xkb4)+xkakb2*
     *       2.2325625d+07)*rho(11-nlp)+(1.0395d+04*(xka4*xka2+xkb4*
     *       xkb2)+9.9225d+05*xkakb2*(xka2+xkb2))*rho(9-nlp)+(4.41d+04*
     *       xkakb2*xkakb2+4.725d+03*xkakb2*(xka4+xkb4))*rho(7-nlp)+
     *       2.1d+02*xkakb2*xkakb2*(xka2+xkb2)*rho(5-nlp)+xkakb2*xkakb2*
     *       xkakb2*rho(3-nlp)
            dum = dum - (xkb*(1.08056025d+08*sigma(14-nlp)+(xka2*
     *       4.9116375d+07+xkb2*1.30977d+07)*sigma(12-nlp)+(2.18295d+05*
     *       (xkb4+1.0d+01*xka4)+5.9535d+06*xkakb2)*sigma(10-nlp)+(xkb2
     *       *xkakb2*9.9225d+04+2.646d+05*xka2*xkakb2+1.0395d+04*xka4*
     *       xka2)*sigma(8-nlp)+xkakb2*(1.26d+03*xka4+4.41d+03*xkakb2)*
     *       sigma(6-nlp)+2.1d+01*xka2*xkakb2*xkakb2*sigma(4-nlp))+xka*
     *       (1.08056025d+08*sigmab(14-nlp)+(1.30977d+07*xka2+xkb2*
     *       4.9116375d+07)*sigmab(12-nlp)+(2.18295d+05*(xka4+1.0d+01*
     *       xkb4)+5.9535d+06*xkakb2)*sigmab(10-nlp)+(xkakb2*(xka2*
     *       9.9225d+04+2.646d+05*xkb2)+1.0395d+04*xkb2*xkb4)*
     *       sigmab(8-nlp)+xkakb2*(1.26d+03*xkb4+4.41d+03*xkakb2)*
     *       sigmab(6-nlp)+2.1d+01*xkakb2*xkakb2*xkb2*sigmab(4-nlp)))
            dum = dum +
     *       xkakb*(1.08056025d+08*tau(13-nlp)+1.30977d+07*(xka2+xkb2)*
     *       tau(11-nlp)+(2.18295d+05*(xka4+xkb4)+1.5876d+06*xkakb2)*
     *       tau(9-nlp)+2.646d+04*xkakb2*(xka2+xkb2)*tau(7-nlp)+
     *       4.41d+02*xkakb2*xkakb2*tau(5-nlp))
            fjpq(ip,ip,1) = dum * (ondab2*ondab2*ondab2*onedab)
          else
            fjpq(ip,ip,1) = fjps(nlp,6,6)
          end if
        else if (ip.ge.8) then
          write(6,*) 'error ip = ',ip,' is too large in fjecp'
          go to 300
        end if
  150 continue
      if (npnp.gt.1) then
         do 200 ik=2,npnp
           do 170 ip=1,2
             do 160 iq=1,2
               fjpq(iq,ip,ik) = fjt(nlp+ik,ip,iq)
  160        continue
  170      continue
           if (lmax.gt.2) then
             do 175 ip=3,lmax
               do 174 iq=1,2
                 fjpq(iq,ip,ik) = fjpq(iq,ip-2,ik)-(2*ip-3)*alfi*
     *                            fjpq(iq,ip-1,ik-1)
  174          continue
  175        continue
             do 177 ip=1,2
               do 176 iq=3,lmax
                 fjpq(iq,ip,ik) = fjpq(iq-2,ip,ik)-(2*iq-3)*beti*
     *                            fjpq(iq-1,ip,ik-1)
  176          continue
  177        continue
c                 
             do 190 ip=3,lmax
               do 180 iq=3,lmax
                 if (ip.gt.iq) then
                   fjpq(iq,ip,ik) = fjpq(iq,ip-2,ik)-(2*ip-3)*alfi*
     *                              fjpq(iq,ip-1,ik-1)
                 else
                   fjpq(iq,ip,ik) = fjpq(iq-2,ip,ik)-(2*iq-3)*beti*
     *                              fjpq(iq-1,ip,ik-1)
                 end if
  180          continue
  190        continue
           end if
  200    continue
      end if
      return
  300 call caserr('dimensioning problem in fjecp - contact authors')
      return
      end
**==fjform.f
      subroutine fjform(rho, sigma, sigmab, tau, nmax)
      implicit real*8 (a-h,o-z)
      dimension rho(15), sigma(15), sigmab(15), tau(15)
      dimension b1(15), b2(15), c1(15), c2(15), facti(13)
c  routine forms all of the needed radial integrals and
c  stores them in rho...tau so that the final integral may
c  be formed elsewhere
      common /fjcmn / alef,bet,xi,xpls,xmns,xp
      common /fjnew / xka, xkb, gamma1,gamma2,a1,a2,c
      data sqpi/1.772453850905d+00/
c  factorial(i)*i
      data facti/1.0d+00,4.0d+00,1.8d+01,9.6d+01,6.0d+02,4.32d+03,
     *           3.528d+04,3.2256d+05,3.26592d+06,3.6288d+07,
     *           4.390848d+08,5.7480192d+09,8.09512704d+10/
c
      sqpidc = sqpi*xi
      b1(1) = gamma1*sqpidc
      b2(1) = gamma2*sqpidc
      c1(1) = sqpidc*gamma1*errf(a1)
      c2(1) = sqpidc*gamma2*errf(a2)
      b1(2) = 2*sqpi*gamma1*dawerf(a1)
      b2(2) = 2*sqpi*gamma2*dawerf(a2)
      c1(2) = 2*sqpi*gamma1*dawf(a1)
      c2(2) = 2*sqpi*gamma2*dawf(a2)
c
c   use recurrence formulas to generate all of the other
c   b and c integrals. note: gamma cancels out the exp in the
c   recurrence formula to leave xp. errors may accumulate if
c   gamma and the exp are used instead.
c
      xkapkb = xka+xkb
      xkamkb = xka-xkb
      do 100 i=3,nmax
        b1(i) = (1+((-1)**i))*(xkapkb**(i-2))/facti(i-2)*
     *          0.25d+00*xp - 2*c*b1(i-2)/(i-2) +
     *          xkapkb*c1(i-1)/(i-2)
        b2(i) = (1+((-1)**i))*(xkamkb**(i-2))/facti(i-2)*
     *          0.25d+00*xp - 2*c*b2(i-2)/(i-2) +
     *          xkamkb*c2(i-1)/(i-2)
        c1(i) = (1-((-1)**i))*(xkapkb**(i-2))/facti(i-2)*
     *          0.25d+00*xp - 2*c*c1(i-2)/(i-2) +
     *          xkapkb*b1(i-1)/(i-2)
        c2(i) = (1-((-1)**i))*(xkamkb**(i-2))/facti(i-2)*
     *          0.25d+00*xp - 2*c*c2(i-2)/(i-2) +
     *          xkamkb*b2(i-1)/(i-2)
 100  continue
c
c  store the integrals in rho...
c
      do 200 i=1,nmax
        rho(i) = b1(i) - b2(i)
        sigma(i) = c1(i) + c2(i)
        sigmab(i) = c1(i) - c2(i)
        tau(i) = b1(i) + b2(i)
 200  continue
      return
      end
**==fjps.f
      function fjps(n,lalf,lbet)
      implicit real*8 (a-h,o-z)
      common /fjcmn / alf,bet,xi,xpls,xmns,xp
      common /fsicmn/ si(30,7)
      dimension dfctrl(7)
      parameter(limab=10)
c      double factorial
      data dfctrl/1.0d+00,3.0d+00,15.0d+00,105.0d+00,945.0d+00,
     *            1.0395d+04,1.35135d+05/
c
c        -----  routine evaluates j-integrals by use of a       -----
c        -----  power series.                                   -----
c        -----  set up the parameters in the power series       -----
c        -----  depending on the relative values of the         -----
c        -----  alfa and beta variables.                        -----
c
      if ((alf-bet).gt.0.0d+00) then
        x = bet*bet
        t = bet**lbet/dfctrl(lbet+1)
        l1= lalf+1
        l2= lbet+n+1
        l3= lbet+lbet+1
      else
        x = alf*alf
        t = alf**lalf/dfctrl(lalf+1)
        l1= lbet+1
        l2= lalf+n+1
        l3= lalf+lalf+1
      end if
c        -----  start the power series sum.                     -----
      sum = t*si(l2,l1)
      do 40 k=1,limab
        t = t*x/(2*k*(k+k+l3))
        sum = sum+t*si(k+k+l2,l1)
   40 continue
      fjps = sum*(0.5d+00*xi)**(n+1)
      return
      end
**==fsi.f
      function fsi(n)
      implicit real*8 (a-h,o-z)
c        -----  routine evalues the scaled i(n,l) integs        -----
c        -----  from analytic formulas,   or if necessary       -----
c        -----  transfers to a power series approach.           -----
c
      save alfa2,errfs,dawfs
c
      common /ficmn / alfa,xi,xp0,xp1
      parameter(alim=0.317d+00)
      data sqpi/1.772453850905d+00/
      data errfs/0.0d+00/
c
c     -------------
      entry fsi0(n)
c     -------------
      if(alfa.gt.alim) go to 5
      fsi0 = fsips(n,0)
      fsi = fsi0
      return
c
    5 np1 = n+1
      if(np1.eq.2) go to 20
      dawfs = dawf(alfa)
      fsi0 = sqpi*dawfs*xp1/alfa
      fsi = fsi0
      return
c
   20 fsi0 = errfs/alfa
      fsi = fsi0
      return
c
c     -------------
      entry fsi1(n)
c     -------------
      if(alfa.gt.alim) go to 25
      fsi1 = fsips(n,1)
      fsi = fsi1
      return
c
   25 np1 = n+1
      if(np1-2) 30,40,50
   30 errfs = sqpi*errf(alfa)*xp1
      fsi1 = (0.5d+00*errfs/alfa-xp0)/alfa
      fsi = fsi1
      return
c
   40 fsi1 = sqpi*(alfa-dawfs)*xp1/(alfa*alfa)
      fsi = fsi1
      return
c
   50 alfa2 = alfa*alfa
      fsi1 = (2.0d+00*alfa*xp0+(2.0d+00*alfa2-1.0d+00)*errfs)/alfa2
      fsi = fsi1
      return
c
c     -------------
      entry fsi2(n)
c     -------------
      if(alfa.gt.alim) go to 55
      fsi2 = fsips(n,2)
      fsi = fsi2
      return
c
   55 np1 = n+1
      go to (60,70,80,90),np1
   60 alfa2 = alfa*alfa
      fsi2=0.25d+00*sqpi*(3.0d+00*alfa-(2.0d+00*alfa2+3.0d+00)*dawfs)
     * *xp1/(alfa2*alfa)
      fsi = fsi2
      return
c
   70 alfa2 = alfa*alfa
      fsi2 = (0.5d+00*(2.0d+00*alfa2-3.0d+00)*errfs+3.0d+00*alfa*xp0)/
     *      (alfa2*alfa)
      fsi = fsi2
      return
c
   80 alfa2 = alfa*alfa
      fsi2 = sqpi*(alfa*(2.0d+00*alfa2-3.0d+00)+3.0d+00*dawfs)*xp1/
     *      (alfa2*alfa)
      fsi = fsi2
      return
c
   90 alfa2 = alfa*alfa
      fsi2 =((alfa2*(4.0d+00*alfa2-4.0d+00)+3.0d+00)*errfs+2.0d+00*alfa*
     *       (2.0d+00*alfa2-3.0d+00)*xp0)/(alfa2*alfa)
      fsi = fsi2
      return
c
c     -------------
      entry fsi3(n)
c     -------------
      if(alfa.gt.alim) go to 95
      fsi3 = fsips(n,3)
      fsi = fsi3
      return
c
   95 np1 = n+1
      go to (100,110,120,130,140),np1
  100 alfa2 = alfa*alfa
      fsi3 = (0.25d+00*(2.0d+00*alfa2-5.0d+00)*errfs+alfa*(4.0d+00*alfa2
     *       +15.0d+00)*xp0/6.0d+00)/(alfa2*alfa2)
      fsi = fsi3
      return
c
  110 alfa2 = alfa*alfa
      fsi3 = 0.25d+00*sqpi*(alfa*(4.0d+00*alfa2-15.0d+00)+3.0d+00*
     *      (2.0d+00*alfa2+5.0d+00)*dawfs)*xp1/(alfa2*alfa2)
      fsi = fsi3
      return
c
  120 alfa2 = alfa*alfa
      fsi3 = (0.5d+00*(alfa2*(4.0d+00*alfa2-12.0d+00)+15.0d+00)*errfs+
     *      alfa*(2.0d+00*alfa2-15.0d+00)*xp0)/(alfa2*alfa2)
      fsi = fsi3
      return
c
  130 alfa2 = alfa*alfa
      fsi3 = sqpi*(alfa*(alfa2*(4.0d+00*alfa2-10.0d+00)+15.0d+00)-
     *      15.0d+00*dawfs)*xp1/(alfa2*alfa2)
      fsi = fsi3
      return
c
  140 alfa2 = alfa*alfa
      fsi3 = (2.0d+00*alfa*(alfa2*(4.0d+00*alfa2-8.0d+00)+15.0d+00)*
     *       xp0+(alfa2*(alfa2*(8.0d+00*alfa2-12.0d+00)+18.0d+00)
     *       -15.0d+00)*errfs)/(alfa2*alfa2)
      fsi = fsi3
      return
      end
**==fsips
      function fsips(n,l)
      implicit real*8 (a-h,o-z)
c        -----  routine evaluates the scaled i(n,l) integs      -----
c        -----  by a power series method.                       -----
c
      common /ficmn / alfa,xi,xp0,xp1
      dimension fctrl(10),dfctrl(10)
      parameter(lima=10)
      data fctrl/1.d+00,1.d+00,2.d+00,6.d+00,24.d+00,120.d+00,720.d+00,
     1           5040.d+00,40320.d+00,362880.d+00/
      data dfctrl/1.d+00,1.d+00,3.d+00,15.d+00,105.d+00,945.d+00,
     1            10395.d+00,135135.d+00,2027025.d+00,34459425.d+00/
      data sqpi/1.772453850905d+00/
c
c        -----  test whether (n+l) is even or odd integer       -----
      nl = n+l
      if(mod(nl,2).eq.1) go to 30
c        -----  lattice of even (n+l) scaled i(n,l) integs.     -----
      lambda = nl/2
      x = alfa*alfa
      l2= l+l+1
c        -----  sum the power series.                           -----
      t = dfctrl(lambda+1)/dfctrl(l+2)
      sum = t
      do 20 k=1,lima
      t = (t*x*(k+k+nl-1))/(k*(k+k+l2))
      sum = sum+t
   20 continue
      fsips = sum*sqpi*xp0*(2.0d+00**lambda)*(alfa**l)
      return
c
c        -----  lattice of odd (n+l) scaled i(n,l) integs.      -----
   30 lambda = (nl-1)/2
      x = 2.0d+00*alfa*alfa
      l2= l+l+1
      t = fctrl(lambda+1)/dfctrl(l+2)
      sum = t
      do 40 k=1,lima
      t = (t*x*(k+lambda))/(k*(k+k+l2))
      sum = sum+t
   40 continue
      fsips = sum*xp0*(2.0d+00**nl)*(alfa**l)
      return
      end
**==ftab.f
      subroutine ftab(fpqr, nmax)
      implicit real*8 (a-h,o-z)
      dimension fpqr(25,25,25)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      integer p,q,r
c        -----  routine sets a table of f-function values.      -----
c        where f means the angular integral given by:
c        int[dw x^r y^s z^t = o for r, s, or t odd
c                           = (r-1)!!(s-1)!!(t-1)!!/(r+s+t+1)!!
c                               for r, s, and t even
c        needs dimension 4*(max real angular momentum) + 1
c        where 1 is added since we really need the 0th element
c        thus for f gradients need g(4) so 4*4+1=17
      parameter(pi=3.1415926535898d+00, one=1.0d+00)
c
      n4max = 4*nmax + 1
      if (n4max.gt.25) then
        write(iwr,999) n4max
 999    format(1x,'ftab: nmax to large =',i3)
        call caserr('invalid dimensions in ecp routine - tfab')
      end if
c        -----  zero out the table.                             -----
      do 10 p=1,n4max
        do 10 q=1,n4max
          do 10 r=1,n4max
   10 fpqr(r,q,p) = 0.0d+00
c        -----  recursively generate non-zero entries.          -----
      fpqr(1,1,1) = 4.0d+00*pi
      do 80 p=1,n4max,2
        pp = p-1
        do 70 q=1,n4max,2
          qq = q-1
          do 60 r=1,n4max,2
            rr = r-1
            if (p.gt.1) then
              fpqr(p,q,r) = (pp-one)*fpqr(p-2,q,r)/(pp+qq+rr+one)
            else if (q.gt.1) then
              fpqr(p,q,r) = (qq-one)*fpqr(p,q-2,r)/(pp+qq+rr+one)
            else if (r.gt.1) then
              fpqr(p,q,r) = (rr-one)*fpqr(p,q,r-2)/(pp+qq+rr+one)
            end if
   60     continue
   70   continue
   80 continue
      return
      end
**==recur.f
      subroutine recur(nmin,nmax,lmax,x)
      implicit real*8 (a-h,o-z)
c        -----  routine uses recursion relations to generate    -----
c        -----  a table of scaled i(n,l) integrals.             -----
c
      common /fsicmn/ si(30,7)
      tx = x+x
      if (lmax.lt.2) then
c
        do 20 n=nmin,nmax,2
          np1 = n+1
          tfnm2 = n+n-4.0d+00
          si(np1,1) = (tfnm2+2.0d+00)*si(n-1,1)+tx*si(n,2)
          si(n+2,2) = tfnm2*si(n,2)+tx*si(np1,1)
   20   continue
        return
      else if (lmax.eq.3) then
        do 60 n=nmin,nmax,2
          np1 = n+1
          np2 = np1+1
          tfnm2 = n+n-4.0d+00
          si(np1,1) = (tfnm2+2.0d+00)*si(n-1,1)+tx*si(n,2)
          si(np2,2) = tfnm2*si(n,2)+tx*si(np1,1)
          si(n+3,3) = tfnm2*si(np1,3)+tx*si(np2,2)
   60   continue
        return
      else
        do 110 n=nmin,nmax,2
          np1 = n+1
          np2 = np1+1
          np3 = np2+1
          tfnm2 = n+n-4.0d+00
          si(np1,1) = (tfnm2+2.0d+00)*si(n-1,1)+tx*si(n,2)
          si(np2,2) = tfnm2*si(n,2)+tx*si(np1,1)
          si(np3,3) = tfnm2*si(np1,3)+tx*si(np2,2)
          si(n+4,4) = tfnm2*si(np2,4)+tx*si(np3,3)
c the following is a guess
          si(n+5,5) = tfnm2*si(np3,5)+tx*si(n+4,4)
          si(n+6,6) = tfnm2*si(n+4,6)+tx*si(n+5,5)
          si(n+7,7) = tfnm2*si(n+5,7)+tx*si(n+6,6)
  110   continue
      end if
      return
      end
**==sitabl.f
      subroutine sitabl(lemax,lomax)
      implicit real*8 (a-h,o-z)
c
c        -----  routine prepares a table of scaled i integs     -----
c        -----  for use in the evaluation of the j integs       -----
c        -----  by a power series.                              -----
c
      common /fjcmn / alfj,betj,xj,xpls,xmns,xp
      common /ficmn / alfi,xi,xp0,xp1
      common /fsicmn/ si(30,7)
c
c        -----  prepare parameters of scaled integs depending   -----
c        -----  on the sign of the alfa and beta variables.     -----
      if ((alfj-betj).le.0.0d+00) then
        alfi = betj
        x  = betj
        xi = xj
        xp0= xp
        xp1= xmns
      else
        x = alfj
        alfi = alfj
        xi  = xj
        xp0 = xp
        xp1 = xpls
      end if
c        -----  generate lattice of scaled i(n,l) integs for    -----
c        -----  which (n+l) is even.                            -----
      si(1,1) = fsi0(0)
      si(2,2) = fsi1(1)
      if (lemax.ge.2) then
        si(3,3) = fsi2(2)
        if (lemax.ge.3) then
          si(4,4) = fsi3(3)
          if (lemax.ge.4) then
            si(5,5) = fsips(4,4)
            if (lemax.ge.5) then
              si(6,6) = fsips(5,5)
              if (lemax.ge.6) then
                si(7,7) = fsips(6,6)
              end if
            end if
          end if
        end if
      end if
      call recur(2,22,lemax,x)
c        -----  generate lattice of scaled i(n,l) integs for    -----
c        -----  for which (n+l) is odd.                         -----
      si(1,2) = fsi1(0)
      si(2,1) = fsi0(1)
      si(3,2) = fsi1(2)
      if (lomax.ge.2) then
        si(4,3) = fsi2(3)
        if (lomax.ge.3) then
          si(5,4) = fsi3(4)
          if (lomax.ge.4) then
            si(6,5) = fsips(4,5)
            if (lomax.ge.5) then
              si(7,6) = fsips(5,6)
              if (lomax.ge.6) then
                si(8,7) = fsips(6,7)
              end if
            end if
          end if
        end if
      end if
      call recur(3,21,lomax,x)
      return
      end
**==zfn.f
      subroutine zfn(zfnlm,lmax,zlm,lmf,lmx,lmy,lmz)
      implicit real*8 (a-h,o-z)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      dimension zlm(*),lmf(*),lmx(*),lmy(*),
     +          lmz(*), zfnlm(125)
c        -----  this routine evaluates the real spherical
c        -----  harmonics given the value of the three
c        -----  cartesian components of a unit vector.
c
      common /zfncm / x,y,z
      data zero/0.0d+00/
c
      if (lmax.le.10) then
        do 500 l=0,lmax
          do 400 m=-l,l
            id = l*(l+1)-m+1
            imn = lmf(id)
            imx = lmf(id+1)-1
            sum = zero
            do 300 i=imn,imx
              dummy = zlm(i)
              if(lmx(i).gt.0) dummy = dummy * (x**lmx(i))
              if(lmy(i).gt.0) dummy = dummy * (y**lmy(i))
              if(lmz(i).gt.0) dummy = dummy * (z**lmz(i))
              sum = sum + dummy
  300       continue
            zfnlm(id) = sum
  400     continue
  500   continue
      else
c
c    if execution reachs this point then an error has occured so abort
c    since zfn does not implement functions higher than l=10
c
        write(iwr,1000) lmax
 1000   format(1x,'error: zfn max lamda=10, you requested lamda = ',i3)
        call caserr('invalid lamda in spherical harmonics - zfn')
      end if
      return
      end
**==ztab.f
      subroutine ztab(zlm)
      implicit real*8 (a-h,o-z)
      dimension zlm(*)
      parameter(fpi=12.566370614359d+00, one=1.0d+00, half=0.5d+00)
c
c        -----  routine sets up the real spherical harmonics    -----
c        -----  in the form of linear combinations of           -----
c        -----  cartesian products.                             -----
c     in other words these are the coefficients of the spherical harmonics
c     after converting to products of cartesian functions.  these are 
c     indexed by lmf and lml to be matched with the correct cartesian
c     powers specified in lmx, lmy, lmz. the maximum l value required is
c     the sum of the maximum l values for the two basis function involved in
c     the integral. thus to do f-gradients you need (3+4)=7 as the maximum l.
c     these will match the formula's given in deck zfn
c
      oned4p = one/fpi
c  l=0, ml=0
      zlm(1) = sqrt(oned4p)
c  l=1, ml=1,0,-1 powers (x, z, y)
      zlm(2) = sqrt(3.0d+00*oned4p)
      zlm(3) = zlm(2)
      zlm(4) = zlm(2)
c  l=2, ml=2...-2 powers ([5]x^2+[6]y^2,[7]xz,[8]z^2+[9],[10]yz,[11]xy)
      zlm(5) = sqrt(15.0d+00*oned4p)*half
      zlm(6) = -zlm(5)
      zlm(7) = 2.0d+00*zlm(5)
      zlm(9) = -sqrt(5.0d+00*oned4p)*half
      zlm(8) = -zlm(9)*3.0d+00
      zlm(10) = zlm(7)
      zlm(11) = zlm(7)
c  l=3 ml=3...-3 powers ([12]x^3+[13]xy^2,[14]x^2z+[15]y^2z,[16]xz^2+[17]x,
c     [18]z^3+[19]z,[20]yz^2+[21]y,[22]xyz,[23]x^2y+[24]y^3)
      zlm(12) = sqrt(35.0d+00*oned4p/8.0d+00)
      zlm(13) = -3.0d+00*zlm(12)
      zlm(14) = sqrt(105.0d+00*oned4p/4.0d+00)
      zlm(15) = -zlm(14)
      zlm(17) = -sqrt(21.0d+00*oned4p/8.0d+00)
      zlm(16) = -zlm(17)*5.0d+00
      zlm(18) = 5.0d+00*sqrt(7.0d+00*oned4p)*half
      zlm(19) = -3.0d+00*zlm(18)*0.2d+00
      zlm(20) = zlm(16)
      zlm(21) = zlm(17)
      zlm(22) = 2.0d+00*zlm(14)
      zlm(23) = -zlm(13)
      zlm(24) = -zlm(12)
c  l=4 ml=4...-4
      zlm(25) = sqrt(315.0d+00*oned4p/64.0d+00)
      zlm(26) = -6.0d+00*zlm(25)
      zlm(27) = zlm(25)
      zlm(28) = sqrt(315.0d+00/(8.0d+00*fpi))
      zlm(29) = -3.0d+00*zlm(28)
      temp    = sqrt(45.0d+00*oned4p)*0.25d+00
      zlm(30) = 7.0d+00*temp
      zlm(31) = -zlm(30)
      zlm(32) = temp
      zlm(33) = -temp
      temp    = sqrt(45.0d+00/(8.0d+00*fpi))
      zlm(34) = 7.0d+00*temp
      zlm(35) = -3.0d+00*temp
      temp    = sqrt(9.0d+00*oned4p)/8.0d+00
      zlm(36) = 35.0d+00*temp
      zlm(37) = -30.0d+00*temp
      zlm(38) = 3.0d+00*temp
      zlm(39) = zlm(34)
      zlm(40) = zlm(35)
      temp    = sqrt(45.0d+00*oned4p*0.25d+00)
      zlm(41) = 7.0d+00*temp
      zlm(42) = -temp
      zlm(43) = -zlm(29)
      zlm(44) = -zlm(28)
      zlm(45) = sqrt(315.0d+00*oned4p*0.25d+00)
      zlm(46) = -zlm(45)
c  l=5
      zlm(47) = sqrt(693.0d+00/(128.0d+00*fpi))
      zlm(48) = -10.0d+00*zlm(47)
      zlm(49) = 5.0d+00*zlm(47)
      zlm(50) = sqrt(3465.0d+00/(64.0d+00*fpi))
      zlm(51) = -6.0d+00*zlm(50)
      zlm(52) = zlm(50)
      temp    = sqrt(385.0d+00/(128.0d+00*fpi))
      zlm(53) = 9.0d+00*temp
      zlm(54) = -27.0d+00*temp
      zlm(55) = 3.0d+00*temp
      zlm(56) = -temp
      temp    = sqrt(1155.0d+00*oned4p)*0.25d+00
      zlm(57) = 3.0d+00*temp
      zlm(58) = -zlm(57)
      zlm(59) = -temp
      zlm(60) = +temp
      temp    = sqrt(165.0d+00*oned4p)/8.0d+00
      zlm(61) = 21.0d+00*temp
      zlm(62) = -14.0d+00*temp
      zlm(63) = temp
      temp    = sqrt(11.0d+00*oned4p)/8.0d+00
      zlm(64) = 63.0d+00*temp
      zlm(65) = -70.0d+00*temp
      zlm(66) = 15.0d+00*temp
      zlm(67) = zlm(61)
      zlm(68) = zlm(62)
      zlm(69) = zlm(63)
      temp    = sqrt(1155.0d+00*oned4p)*half
      zlm(70) = 3.0d+00*temp
      zlm(71) = -temp
      zlm(72) = -zlm(54)
      zlm(73) = -zlm(55)
      zlm(74) = -zlm(53)
      zlm(75) = -zlm(56)
      zlm(76) = sqrt(3465.0d+00*oned4p)*half
      zlm(77) = -zlm(76)
      zlm(78) = zlm(49)
      zlm(79) = zlm(48)
      zlm(80) = zlm(47)
c l=6
      zlm(81) = sqrt(3003.0d+00/(512.0d+00*fpi))
      zlm(82) = -15.0d+00*zlm(81)
      zlm(83) = -zlm(82)
      zlm(84) = -zlm(81)
      zlm(85) = sqrt(9009.0d+00/(128.0d+00*fpi))
      zlm(86) = -10.0d+00*zlm(85)
      zlm(87) = 5.0d+00*zlm(85)
      temp    = sqrt(819.0d+00/(256.0d+00*fpi))
      zlm(88) = 11.0d+00*temp
      zlm(89) = -66.0d+00*temp
      zlm(90) = zlm(88)
      zlm(91) = -temp
      zlm(92) = 6.0d+00*temp
      zlm(93) = -temp
      temp    = sqrt(1365.0d+00/(128.0d+00*fpi))
      zlm(94) = 11.0d+00*temp
      zlm(95) = -33.0d+00*temp
      zlm(96) = 9.0d+00*temp
      zlm(97) = -3.0d+00*temp
      temp    = sqrt(1365.0d+00/(512.0d+00*fpi))
      zlm(98) = 33.0d+00*temp
      zlm(99) = -zlm(98)
      zlm(100) = -18.0d+00*temp
      zlm(101) = +18.0d+00*temp
      zlm(102) = temp
      zlm(103) = -temp
      temp     = sqrt(273.0d+00*oned4p)/8.0d+00
      zlm(104) = 33.0d+00*temp
      zlm(105) = -30.0d+00*temp
      zlm(106) = 5.0d+00*temp
      temp     = sqrt(13.0d+00*oned4p)/16.0d+00
      zlm(107) = 231.0d+00*temp
      zlm(108) = -315.0d+00*temp
      zlm(109) = 105.0d+00*temp
      zlm(110) = -5.0d+00*temp
      zlm(111) = zlm(104)
      zlm(112) = zlm(105)
      zlm(113) = zlm(106)
      temp     = sqrt(1365.0d+00/(128.0d+00*fpi))
      zlm(114) = 33.0d+00*temp
      zlm(115) = -18.0d+00*temp
      zlm(116) = temp
      zlm(117) = -zlm(95)
      zlm(118) = -zlm(94)
      zlm(119) = -zlm(96)
      zlm(120) = -zlm(97)
      temp     = sqrt(819.0d+00*oned4p)*0.25d+00
      zlm(121) = 11.0d+00*temp
      zlm(122) = -zlm(121)
      zlm(123) = -temp
      zlm(124) = temp
      zlm(125) = zlm(87)
      zlm(126) = zlm(86)
      zlm(127) = zlm(85)
      temp    = sqrt(3003.0d+00/(512.0d+00*fpi))
      zlm(128) = 6.0d+00*temp
      zlm(129) = -20.0d+00*temp
      zlm(130) = zlm(128)
c....l=7
c....ml=-7
      zlm(131)=sqrt(6.435d+03*oned4p/1.024d+03)
      zlm(132)=zlm(131)*35.0d+00
      zlm(133)=-zlm(131)*21.0d+00
      zlm(134)=-zlm(131)*7.0d+00
c...ml=-6
      zlm(135)=sqrt(4.5045d+04*oned4p/5.12d+02)
      zlm(136)=-zlm(135)*15.0d+00
      zlm(137)=-zlm(135)
      zlm(138)=-zlm(136)
c...ml=-5
      temp=sqrt(3.465d+03*oned4p/1.024d+03)
      zlm(139)=13.0d+00*temp
      zlm(140)=-10.0d+00*zlm(139)
      zlm(141)=5.0d+00*zlm(139)
      zlm(142)=-temp
      zlm(143)=-10.0d+00*zlm(142)
      zlm(144)=5.0d+00*zlm(142)
c...ml=-4
      temp=sqrt(3.465d+03*oned4p/2.56d+02)
      zlm(145)=13.0d+00*temp
      zlm(146)=-78.0d+00*temp
      zlm(147)=zlm(145)
      zlm(148)=-3.0d+00*temp
      zlm(149)=18.0d+00*temp
      zlm(150)=zlm(148)
c...ml=-3
      temp=sqrt(3.15d+02*oned4p/1.024d+03)
      zlm(151)=temp*143.0d+00
      zlm(152)=-zlm(151)*3.0d+00
      zlm(153)=-66.0d+00*temp
      zlm(154)=-3.0d+00*zlm(153)
      zlm(155)=3.0d+00*temp
      zlm(156)=-9.0d+00*temp
c...ml=-2
      temp=sqrt(3.15d+02*oned4p/5.12d+02)
      zlm(157)=temp*1.43d+02
      zlm(158)=-zlm(157)
      zlm(159)=-110.0d+00*temp
      zlm(160)=-zlm(159)
      zlm(161)=1.5d+01*temp
      zlm(162)=-zlm(161)
c....ml=-1
      temp=sqrt(1.05d+02*oned4p/1.024d+03)
      zlm(163)=429.0d+00*temp
      zlm(164)=-495.0d+00*temp
      zlm(165)=135.0d+00*temp
      zlm(166)=-5.0d+00*temp
c....ml=0
      temp=sqrt(1.5d+01*oned4p/2.56d+02)
      zlm(167)=temp*429.0d+00
      zlm(168)=-temp*693.0d+00
      zlm(169)=temp*315.0d+00
      zlm(170)=-temp*35.0d+00
c....ml=1
      zlm(171)=zlm(163)
      zlm(172)=zlm(164)
      zlm(173)=zlm(165)
      zlm(174)=zlm(166)
c....ml=2
      zlm(175)=zlm(157)*2.0d+00
      zlm(176)=zlm(159)*2.0d+00
      zlm(177)=zlm(161)*2.0d+00
c....ml=3
      zlm(178)=-zlm(152)
      zlm(179)=-zlm(151)
      zlm(180)=-zlm(154)
      zlm(181)=-zlm(153)
      zlm(182)=-zlm(156)
      zlm(183)=-zlm(155)
c....ml=4
      zlm(184)=zlm(145)*4.0d+00
      zlm(185)=-zlm(184)
      zlm(186)=zlm(148)*4.0d+00
      zlm(187)=-zlm(186)
c....ml=5
      zlm(188)=zlm(141)
      zlm(189)=zlm(140)
      zlm(190)=zlm(139)
      zlm(191)=zlm(144)
      zlm(192)=zlm(143)
      zlm(193)=zlm(142)
c....ml=6
      temp = sqrt(4.5045d+04*oned4p/1.28d+02)
      zlm(194)=3.0d+00*temp
      zlm(195)=-10.0d+00*temp
      zlm(196)=zlm(194)
c....ml=-7
      zlm(197)=-zlm(131)
      zlm(198)=-zlm(132)
      zlm(199)=-zlm(133)
      zlm(200)=-zlm(134)
c l=8, ml=-8
      zlm(201)=(3.0d+00/1.28d+02)*sqrt(1.2155d+04*oned4p)
      zlm(202)=-2.8d+01*zlm(201)
      zlm(203)= 7.0d+01*zlm(201)
      zlm(204)=zlm(202)
      zlm(205)=zlm(201)
c ml=-7
      zlm(206)=(3.0d+00/3.2d+01)*sqrt(1.2155d+04*oned4p)
      zlm(207)=-2.1d+01*zlm(206)
      zlm(208)= 3.5d+01*zlm(206)
      zlm(209)=-7.0d+00*zlm(206)
c ml=-6
      temp=(sqrt(7.293d+03*oned4p*0.5d+00))/3.2d+01
      zlm(210)=1.5d+01*temp
      zlm(211)=-1.5d+01*zlm(210)
      zlm(212)=-zlm(211)
      zlm(213)=-zlm(210)
      zlm(214)=-temp
      zlm(215)=zlm(210)
      zlm(216)=-zlm(215)
      zlm(217)=temp
c ml=-5
      temp=(3.0d+00/3.2d+01)*sqrt(1.7017d+04*oned4p)
      zlm(218)=5.0d+00*temp
      zlm(219)=-1.0d+01*zlm(218)
      zlm(220)=5.0d+00*zlm(218)
      zlm(221)=-temp
      zlm(222)=1.0d+01*temp
      zlm(223)=5.0d+00*zlm(221)
c ml=-4
      temp=(3.0d+00/6.4d+01)*sqrt(1.309d+03*oned4p)
      zlm(224)=6.5d+01*temp
      zlm(225)=zlm(224)
      zlm(226)=-6.0d+00*zlm(224)
      zlm(227)=-26.0d+00*temp
      zlm(228)=zlm(227)
      zlm(229)=-6.0d+00*zlm(227)
      zlm(230)=temp
      zlm(231)=temp
      zlm(232)=-6.0d+00*temp
c ml=-3
      temp=(sqrt(1.9635d+04*oned4p))/3.2d+01
      zlm(233)=3.9d+01*temp
      zlm(234)=-3.0d+00*zlm(233)
      zlm(235)=-2.6d+01*temp
      zlm(236)=-3.0d+00*zlm(235)
      zlm(237)=3.0d+00*temp
      zlm(238)=-3.0d+00*zlm(237)
c ml=-2
      temp=(3.0d+00/3.2d+01)*sqrt(5.95d+02*oned4p*0.5d+00)
      zlm(239)=1.43d+02*temp
      zlm(240)=-zlm(239)
      zlm(241)=zlm(240)
      zlm(242)=zlm(239)
      zlm(243)=3.3d+01*temp
      zlm(244)=-zlm(243)
      zlm(245)=-temp
      zlm(246)=temp
c ml=-1
      temp=(3.0d+00/3.2d+01)*sqrt(1.7d+01*oned4p)
      zlm(247)=7.15d+02*temp
      zlm(248)=-1.001d+03*temp
      zlm(249)=3.85d+02*temp
      zlm(250)=-3.5d+01*temp
c ml=0
      temp=(sqrt(1.7d+01*oned4p))/1.28d+02
      zlm(251)=6.435d+03*temp
      zlm(252)=-1.2012d+04*temp
      zlm(253)=6.930d+03*temp
      zlm(254)=-1.260d+03*temp
      zlm(255)=3.5d+01*temp
c ml=1
      zlm(256)=zlm(247)
      zlm(257)=zlm(248)
      zlm(258)=zlm(249)
      zlm(259)=zlm(250)
c ml=2
      temp=(3.0d+00/1.6d+01)*sqrt(5.95d+02*oned4p*0.5d+00)
      zlm(260)=1.43d+02*temp
      zlm(261)=-zlm(260)
      zlm(262)=3.3d+01*temp
      zlm(263)=-temp
c ml=3
      zlm(264)=-zlm(234)
      zlm(265)=-zlm(233)
      zlm(266)=-zlm(236)
      zlm(267)=-zlm(235)
      zlm(268)=-zlm(238)
      zlm(269)=-zlm(237)
c ml=4
      temp=(3.0d+00/1.6d+01)*sqrt(1.309d+03*oned4p)
      zlm(270)=6.5d+01*temp
      zlm(271)=-zlm(270)
      zlm(272)=-26.0d+00*temp
      zlm(273)=-zlm(272)
      zlm(274)=temp
      zlm(275)=-temp
c ml=5
      temp=(3.0d+00/3.2d+01)*sqrt(1.7017d+04*oned4p)
      zlm(276)=2.5d+01*temp
      zlm(277)=-5.0d+01*temp
      zlm(278)=5.0d+00*temp
      zlm(279)=-5.0d+00*temp
      zlm(280)=1.0d+01*temp
      zlm(281)=-temp
c ml=6
      temp=(sqrt(7.293d+03*oned4p*0.5d+00))/1.6d+01
      zlm(282)=4.5d+01*temp
      zlm(283)=-1.5d+02*temp
      zlm(284)=zlm(282)
      zlm(285)=-3.0d+00*temp
      zlm(286)=1.0d+01*temp
      zlm(287)=zlm(285)
c ml=7
      zlm(288)=-zlm(209)
      zlm(289)=-zlm(208)
      zlm(290)=-zlm(207)
      zlm(291)=-zlm(206)
c ml=8
      zlm(292)=8.0d+00*zlm(201)
      zlm(293)=-7.0d+00*zlm(292)
      zlm(294)=-zlm(293)
      zlm(295)=-zlm(292)
c**** l = 9  ml=-9
      zlm(296)=(sqrt(2.30945d+05*oned4p*0.5d+00))/1.28d+02
      zlm(297)=-3.6d+01*zlm(296)
      zlm(298)=1.26d+02*zlm(296)
      zlm(299)=-8.4d+01*zlm(296)
      zlm(300)=9.0d+00*zlm(296)
c ml=-8
      zlm(301)=(3.0d+00/1.28d+02)*sqrt(2.30945d+05*oned4p)
      zlm(302)=-2.8d+01*zlm(301)
      zlm(303)=7.0d+01*zlm(301)
      zlm(304)=zlm(302)
      zlm(305)=zlm(301)
c ml=-7
      temp=(3.0d+00/1.28d+02)*sqrt(1.3585d+04*oned4p*0.5d+00)
      zlm(306)=1.7d+01*temp
      zlm(307)=-2.1d+01*zlm(306)
      zlm(308)=3.5d+01*zlm(306)
      zlm(309)=-7.0d+00*zlm(306)
      zlm(310)=-temp
      zlm(311)=2.1d+01*temp
      zlm(312)=-3.5d+01*temp
      zlm(313)=7.0d+00*temp
c ml=-6
      temp=(sqrt(4.0755d+04*oned4p*0.5d+00))/3.2d+01
      zlm(314)=1.7d+01*temp
      zlm(315)=-1.5d+01*zlm(314)
      zlm(316)=-zlm(315)
      zlm(317)=-zlm(314)
      zlm(318)=-3.0d+00*temp
      zlm(319)=-1.5d+01*zlm(318)
      zlm(320)=-zlm(319)
      zlm(321)=-zlm(318)
c ml=-5
      temp=(3.0d+00/6.4d+01)*sqrt(2.717d+03*oned4p*0.5d+00)
      zlm(322)=8.5d+01*temp
      zlm(323)=-1.0d+01*zlm(322)
      zlm(324)=5.0d+00*zlm(322)
      zlm(325)=-3.0d+01*temp
      zlm(326)=-1.0d+01*zlm(325)
      zlm(327)=5.0d+00*zlm(325)
      zlm(328)=temp
      zlm(329)=-1.0d+01*temp
      zlm(330)=5.0d+00*temp
c ml=-4
      temp=(3.0d+00/6.4d+01)*sqrt(9.5095d+04*oned4p)
      zlm(331)=1.7d+01*temp
      zlm(332)=zlm(331)
      zlm(333)=-6.0d+00*zlm(331)
      zlm(334)=-1.0d+01*temp
      zlm(335)=zlm(334)
      zlm(336)=-6.0d+00*zlm(334)
      zlm(337)=temp
      zlm(338)=temp
      zlm(339)=-6.0d+00*temp
c ml=-3
      temp=(sqrt(2.1945d+04*oned4p*0.5d+00))/6.4d+01
      zlm(340)=2.21d+02*temp
      zlm(341)=-3.0d+00*zlm(340)
      zlm(342)=-1.95d+02*temp
      zlm(343)=-3.0d+00*zlm(342)
      zlm(344)=3.9d+01*temp
      zlm(345)=-3.0d+00*zlm(344)
      zlm(346)=-temp
      zlm(347)=3.0d+00*temp
c ml=-2
      temp=(3.0d+00/3.2d+01)*sqrt(1.045d+03*oned4p*0.5d+00)
      zlm(348)=2.21d+02*temp
      zlm(349)=-zlm(348)
      zlm(350)=-2.73d+02*temp
      zlm(351)=-zlm(350)
      zlm(352)=9.1d+01*temp
      zlm(353)=-zlm(352)
      zlm(354)=-7.0d+00*temp
      zlm(355)=-zlm(354)
c ml=-1
      temp=(3.0d+00/1.28d+02)*sqrt(9.5d+01*oned4p)
      zlm(356)=2.431d+03*temp
      zlm(357)=-4.004d+03*temp
      zlm(358)=2.002d+03*temp
      zlm(359)=-3.08d+02*temp
      zlm(360)=7.0d+00*temp
c ml=0
      temp=(sqrt(1.9d+01*oned4p))/1.28d+02
      zlm(361)=1.2155d+04*temp
      zlm(362)=-2.5740d+04*temp
      zlm(363)=1.8018d+04*temp
      zlm(364)=-4.620d+03*temp
      zlm(365)=3.15d+02*temp
c ml=1
      zlm(366)=zlm(356)
      zlm(367)=zlm(357)
      zlm(368)=zlm(358)
      zlm(369)=zlm(359)
      zlm(370)=zlm(360)
c ml=2
      temp=(3.0d+00/1.6d+01)*sqrt(1.045d+03*oned4p*0.5d+00)
      zlm(371)=2.21d+02*temp
      zlm(372)=-2.73d+02*temp
      zlm(373)=9.1d+01*temp
      zlm(374)=-7.0d+00*temp
c ml=3
      zlm(375)=-zlm(341)
      zlm(376)=-zlm(340)
      zlm(377)=-zlm(343)
      zlm(378)=-zlm(342)
      zlm(379)=-zlm(345)
      zlm(380)=-zlm(344)
      zlm(381)=-zlm(347)
      zlm(382)=-zlm(346)
c ml=4
      temp=(3.0d+00/1.6d+01)*sqrt(9.5095d+04*oned4p)
      zlm(383)=1.7d+01*temp
      zlm(384)=-zlm(383)
      zlm(385)=-1.0d+01*temp
      zlm(386)=-zlm(385)
      zlm(387)=temp
      zlm(388)=-zlm(387)
c ml=5
      zlm(389)=zlm(324)
      zlm(390)=zlm(323)
      zlm(391)=zlm(322)
      zlm(392)=zlm(327)
      zlm(393)=zlm(326)
      zlm(394)=zlm(325)
      zlm(395)=zlm(330)
      zlm(396)=zlm(329)
      zlm(397)=zlm(328)
c ml=6
      temp=(sqrt(4.0755d+04*oned4p*0.5d+00))/1.6d+01
      zlm(398)=5.1d+01*temp
      zlm(399)=-1.7d+02*temp
      zlm(400)=zlm(398)
      zlm(401)=-9.0d+00*temp
      zlm(402)=3.0d+01*temp
      zlm(403)=zlm(401)
c ml=7
      zlm(404)=-zlm(309)
      zlm(405)=-zlm(308)
      zlm(406)=-zlm(307)
      zlm(407)=-zlm(306)
      zlm(408)=-zlm(313)
      zlm(409)=-zlm(312)
      zlm(410)=-zlm(311)
      zlm(411)=-zlm(310)
c ml=8
      zlm(412)=8.0d+00*zlm(301)
      zlm(413)=-7.0d+00*zlm(412)
      zlm(414)=-zlm(413)
      zlm(415)=-zlm(412)
c ml=9
      zlm(416)=zlm(300)
      zlm(417)=zlm(299)
      zlm(418)=zlm(298)
      zlm(419)=zlm(297)
      zlm(420)=zlm(296)
c*** l=10, ml=-10
      zlm(421)=(sqrt(9.69969d+05*oned4p*0.5d+00))/2.56d+02
      zlm(422)=-4.5d+01*zlm(421)
      zlm(423)=2.10d+02*zlm(421)
      zlm(424)=-zlm(423)
      zlm(425)=-zlm(422)
      zlm(426)=-zlm(421)
c ml=-9
      zlm(427)=(sqrt(4.849845d+06*oned4p*0.5d+00))/1.28d+02
      zlm(428)=-3.6d+01*zlm(427)
      zlm(429)=1.26d+02*zlm(427)
      zlm(430)=-8.4d+01*zlm(427)
      zlm(431)=9.0d+00*zlm(427)
c ml=-8
      temp=(sqrt(2.55255d+05*oned4p))/2.56d+02
      zlm(432)=1.9d+01*temp
      zlm(433)=-2.8d+01*zlm(432)
      zlm(434)=7.0d+01*zlm(432)
      zlm(435)=zlm(433)
      zlm(436)=zlm(432)
      zlm(437)=-temp
      zlm(438)=-2.8d+01*zlm(437)
      zlm(439)=7.0d+01*zlm(437)
      zlm(440)=zlm(438)
      zlm(441)=zlm(437)
c ml=-7
      temp=(3.0d+00/1.28d+02)*sqrt(8.5085d+04*oned4p*0.5d+00)
      zlm(442)=1.9d+01*temp
      zlm(443)=-2.1d+01*zlm(442)
      zlm(444)=3.5d+01*zlm(442)
      zlm(445)=-7.0d+00*zlm(442)
      zlm(446)=-3.0d+00*temp
      zlm(447)=-2.1d+01*zlm(446)
      zlm(448)=3.5d+01*zlm(446)
      zlm(449)=-7.0d+00*zlm(446)
c ml=-6
      temp=(3.0d+00/2.56d+02)*sqrt(5.005d+03*oned4p*0.5d+00)
      zlm(450)=3.23d+02*temp
      zlm(451)=-1.5d+01*zlm(450)
      zlm(452)=-zlm(451)
      zlm(453)=-zlm(450)
      zlm(454)=-1.02d+02*temp
      zlm(455)=-1.5d+01*zlm(454)
      zlm(456)=-zlm(455)
      zlm(457)=-zlm(454)
      zlm(458)=3.0d+00*temp
      zlm(459)=-1.5d+01*zlm(458)
      zlm(460)=-zlm(459)
      zlm(461)=-zlm(458)
c ml=-5
      temp=(3.0d+00/6.4d+01)*sqrt(1.001d+03*oned4p*0.5d+00)
      zlm(462)=3.23d+02*temp
      zlm(463)=-1.0d+01*zlm(462)
      zlm(464)=5.0d+00*zlm(462)
      zlm(465)=-1.70d+02*temp
      zlm(466)=-1.0d+01*zlm(465)
      zlm(467)=5.0d+00*zlm(465)
      zlm(468)=1.5d+01*temp
      zlm(469)=-1.0d+01*zlm(468)
      zlm(470)=5.0d+00*zlm(468)
c ml=-4
      temp=(3.0d+00/1.28d+02)*sqrt(5.005d+03*oned4p)
      zlm(471)=3.23d+02*temp
      zlm(472)=zlm(471)
      zlm(473)=-6.0d+00*zlm(471)
      zlm(474)=-2.55d+02*temp
      zlm(475)=zlm(474)
      zlm(476)=-6.0d+00*zlm(474)
      zlm(477)=4.5d+01*temp
      zlm(478)=zlm(477)
      zlm(479)=-6.0d+00*zlm(477)
      zlm(480)=-temp
      zlm(481)=zlm(480)
      zlm(482)=-6.0d+00*zlm(480)
c ml=-3
      temp=(3.0d+00/6.4d+01)*sqrt(5.005d+03*oned4p*0.5d+00)
      zlm(483)=3.23d+02*temp
      zlm(484)=-3.0d+00*zlm(483)
      zlm(485)=-3.57d+02*temp
      zlm(486)=-3.0d+00*zlm(485)
      zlm(487)=1.05d+02*temp
      zlm(488)=-3.0d+00*zlm(487)
      zlm(489)=-7.0d+00*temp
      zlm(490)=-3.0d+00*zlm(489)
c ml=-2
      temp=(3.0d+00/2.56d+02)*sqrt(3.85d+02*oned4p)
      zlm(491)=4.199d+03*temp
      zlm(492)=-zlm(491)
      zlm(493)=-6.188d+03*temp
      zlm(494)=-zlm(493)
      zlm(495)=2.730d+03*temp
      zlm(496)=-zlm(495)
      zlm(497)=-3.64d+02*temp
      zlm(498)=-zlm(497)
      zlm(499)=7.0d+00*temp
      zlm(500)=-zlm(499)
c ml=-1
      temp=(sqrt(1.155d+03*oned4p))/1.28d+02
      zlm(501)=4.199d+03*temp
      zlm(502)=-7.956d+03*temp
      zlm(503)=4.914d+03*temp
      zlm(504)=-1.092d+03*temp
      zlm(505)=6.3d+01*temp
c ml=0
      temp=(sqrt(2.1d+01*oned4p))/2.56d+02
      zlm(506)=4.6189d+04*temp
      zlm(507)=-1.09395d+05*temp
      zlm(508)=9.0090d+04*temp
      zlm(509)=-3.0030d+04*temp
      zlm(510)=3.465d+03*temp
      zlm(511)=-6.3d+01*temp
c ml=1
      zlm(512)=zlm(501)
      zlm(513)=zlm(502)
      zlm(514)=zlm(503)
      zlm(515)=zlm(504)
      zlm(516)=zlm(505)
c ml=2
      zlm(517)=2.0d+00*zlm(491)
      zlm(518)=2.0d+00*zlm(493)
      zlm(519)=2.0d+00*zlm(495)
      zlm(520)=2.0d+00*zlm(497)
      zlm(521)=2.0d+00*zlm(499)
c ml=3
      zlm(522)=-zlm(484)
      zlm(523)=-zlm(483)
      zlm(524)=-zlm(486)
      zlm(525)=-zlm(485)
      zlm(526)=-zlm(488)
      zlm(527)=-zlm(487)
      zlm(528)=-zlm(490)
      zlm(529)=-zlm(489)
c ml=4
      zlm(530)=4.0d+00*zlm(471)
      zlm(531)=-zlm(530)
      zlm(532)=4.0d+00*zlm(474)
      zlm(533)=-zlm(532)
      zlm(534)=4.0d+00*zlm(477)
      zlm(535)=-zlm(534)
      zlm(536)=4.0d+00*zlm(480)
      zlm(537)=-zlm(536)
c ml=5
      zlm(538)=zlm(462)
      zlm(539)=zlm(463)
      zlm(540)=zlm(464)
      zlm(541)=zlm(465)
      zlm(542)=zlm(466)
      zlm(543)=zlm(467)
      zlm(544)=zlm(468)
      zlm(545)=zlm(469)
      zlm(546)=zlm(470)
c ml=6
      zlm(547)=6.0d+00*zlm(450)
      zlm(548)=-2.0d+01*zlm(450)
      zlm(549)=zlm(547)
      zlm(550)=6.0d+00*zlm(454)
      zlm(551)=-2.0d+01*zlm(454)
      zlm(552)=zlm(550)
      zlm(553)=6.0d+00*zlm(458)
      zlm(554)=-2.0d+01*zlm(458)
      zlm(555)=zlm(553)
c ml=7
      zlm(556)=-zlm(445)
      zlm(557)=-zlm(444)
      zlm(558)=-zlm(443)
      zlm(559)=-zlm(442)
      zlm(560)=-zlm(449)
      zlm(561)=-zlm(448)
      zlm(562)=-zlm(447)
      zlm(563)=-zlm(446)
c ml=8
      zlm(564)=8.0d+00*zlm(432)
      zlm(565)=-5.6d+01*zlm(432)
      zlm(566)=-zlm(565)
      zlm(567)=-zlm(564)
      zlm(568)=8.0d+00*zlm(437)
      zlm(569)=-5.6d+01*zlm(437)
      zlm(570)=-zlm(569)
      zlm(571)=-zlm(568)
c ml=9
      zlm(572)=zlm(431)
      zlm(573)=zlm(430)
      zlm(574)=zlm(429)
      zlm(575)=zlm(428)
      zlm(576)=zlm(427)
c ml=10
      zlm(577)=1.0d+01*zlm(421)
      zlm(578)=-1.20d+02*zlm(421)
      zlm(579)=2.52d+02*zlm(421)
      zlm(580)=zlm(578)
      zlm(581)=zlm(577)
      return
      end
**==ecp1d.f
      subroutine ecp1d(de,datot,l2,dcoef1,jfst1,lbecp1,dcoef4,
     *          dcoef2,jfst2,lbecp2,deloc,fpqr,fp,fp2,g,g2,xin,yin,
     *          zin,zlm,lmf,lmx,lmy,lmz,mapshl,nshels)
      implicit real*8 (a-h,o-z)
      logical out,outall,canda,aandb
c
      dimension datot(l2),deloc(3,*),jfst1(*),jfst2(*),
     *          dcoef1(*),lbecp1(9,*),dcoef2(*),lbecp2(6,*),dcoef4(*)
      dimension xin(225),yin(225),zin(225),g(315),g2(315),
     *          fpqr(25,25,25), de(3,*)
      dimension zlm(*), lmf(*), lmx(*), lmy(*), lmz(*)
      dimension mapshl(nshels,*)
c  fp is used for radial integral storage and must be at least 2*11**3
      dimension fp(2662), fp2(2662)
c
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
c
      character *8 zcom,ztitle,zaname,ztag,zsymm,zgroup
      character *10 zbflab
      character *8 zscftp,zruntp,zguess,zconf,zstate,zorb,zpseud
      common /runlab/ zcom(19),ztitle(10),zaname(maxat),zbflab(maxorb),
     +   ztag(maxat),zsymm,zgroup,zscftp,zruntp,zguess,zconf,zstate,
     +   zorb(maxorb),zpseud(maxat)
c
c....  zaname : names of atoms during calculation (after reorder)
c....  ztag   : names of atoms as read in 
c
c
c...ecp parameters
c
c
      real*8 zetc,cax,cay,caz,ca,xca,yca,zca
      real*8 zetb,bax,bay,baz,ba,xba,yba,zba
      real*8 phase,dax,day,daz,da,xda,yda,zda,xint
      integer kcntr
      common /ecp1/ zetc,cax,cay,caz,ca,xca,yca,zca,
     +              zetb,bax,bay,baz,ba,xba,yba,zba,
     +              phase,dax,day,daz,da,xda,yda,zda,xint,
     +              kcntr
c
      common /ecp2  / clp(400),zlp(400),nlp(400),kfirst(maxat,6),
     *                klast(maxat,6),lmax(maxat),lpskip(maxat),
     *                izcore(maxat)
      logical iandj,norm,normi,normj
      common /ecpidx/ q2,iamin,iamax,jamin,jamax,ipmin,ipmax,jpmin,
     *                jpmax,kf1,kl1,llmx,npc,npb,iandj,norm,normi,normj
c
      real*8 bmcx,bmcy,bmcz,bpcx,bpcy,bpcz,cbsq,ax,ay,az
      logical candb
      common /ecp4/ bmcx,bmcy,bmcz,bpcx,bpcy,bpcz,cbsq,ax,ay,az,
     +              candb
c
c   ecp common stuff
      common /zfncm / x,y,z
c
c  the following are to make use of symmetry
c
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
c
c     storage for derivative fock matrix output
      common/blkin/gout(5109),nword
c
      logical oianj, okanl, oident, omisc, oham, opdipd, omp2
      integer ipos1, ipos2
      common /misc/ oianj,okanl,oident,omisc,
     +            oham,opdipd,omp2,ipos1,ipos2
c
c
c  local storage
      dimension coefi(mxprms),coefj(mxprms),coefip(mxprms),iang(56),
     *          iamina(7)
c
      parameter (zero=0.0d+00)
c
      character *3 dnam
      dimension dnam(3)
      data dnam /'e*x','e*y','e*z'/
      data m5110/5107/
c
c    iang gives the maximum angular momentum for each shell
c
      data iang/1,3*2,6*3,10*4,15*5, 21*6 /
      data iamina/1,2,5,11,21,36,57/
c
c        -----  routine calculates the gradient of the ecp      -----
c        -----  integrals and produces an output vector.        -----
c        -----  routine based on new version of ecpint, by      -----
c        -----  l.    kahn, of battele.                         -----
c        -----  gradient modifications by ims group in japan.   -----
c        -----  gradscf version  august  -1981-                 -----
c        reworked for gamess by brett bode 1998
c
      out = nprint.eq. - 3 .or. nprint.eq. - 10
      outall = nprint.eq. - 10
c
      norm=normf.ne.1.or.normp.ne.1
      nword = 1
c
c  the shell being differentiated will be normalized when the derivative
c  is formed in formdr, so don't norm in the radial routines
      normi = .false.
      normj = norm
      iandj = .false.
c   zero out the gradient vector
      n3x = nat*3
      call vclr(deloc,1,n3x)
c
c        -----  loop over  ishell.                              -----
c        -----  note loops run over full square array           -----
c        -----  of ecp matrix for derivatives.                  -----
c
      do 9000 ii=1,nshell
c
       i1 = kstart(ii)
       i2 = i1+kng(ii)-1
       ipmin = i1
       ipmax = i2
       icntr = katom(ii)
       imin = kmin(ii)
       imax = kmax(ii)
       loci = kloc(ii)-imin
       iimax = 1
       if (imin.eq.1.and.imax.eq.4) iimax = 2
       do 8900 iii=1,iimax
        if (imin.eq.1.and.imax.eq.4) then
         if (iii.eq.1) then
           iamin = 1
           iamax = 1
         else
           iamin = 2
           iamax = 4
         end if
        else
         iamin = imin
         iamax = imax
        end if
        npc=iang(iamax)
c  store the coefs for later use. since we are taking the derivative
c  of the ishell we can combine the coef with the factor of 2ex
        igii=1
        do 111 jj=ipmin,ipmax
         if (iamin.eq.1) then
           coefi(igii)=cs(jj)
           itemp = 1
         else if (iamin.lt.5) then
           coefi(igii)=cp(jj)
           itemp = 2
         else if (iamin.lt.11) then
           coefi(igii)=cd(jj)
           itemp = 3
         else if (iamin.le.20) then
           coefi(igii)=cf(jj)
           itemp = 4
         else if (iamin.le.35) then
           coefi(igii)=cg(jj)
           itemp = 5
         end if
         coefip(igii) = coefi(igii)*(-2.0d+00*ex(jj))
         igii = igii + 1
 111    continue
c
c        -----  jshell  -----
c
        do 8000 jj=1,nshell
c    check symmetry
         n2=0
         ii0 = max(ii,jj)
         jj0 = min(ii,jj)
         do 80 it=1,nt
          id=mapshl(ii,it)
          jd=mapshl(jj,it)
          idd = max(id,jd)
          jdd = min(id,jd)
          if(idd.gt.ii0) go to 8000
          if(idd.lt.ii0) go to 80
          if(jdd.gt.jj0) go to 8000
          if(jdd.lt.jj0) go to 80
          n2=n2+1
 80      continue
         q2 = dble(nt)
         q2 = q2/dble(n2)
c
c
         j1 = kstart(jj)
         j2 = j1+kng(jj)-1
         jpmin = j1
         jpmax = j2
         jcntr = katom(jj)
         jmin = kmin(jj)
         jmax = kmax(jj)
         locj = kloc(jj)-jmin
         jjmax = 1
         if (jmin.eq.1.and.jmax.eq.4) jjmax = 2
         do 7900 jjj=1,jjmax
          if (jmin.eq.1.and.jmax.eq.4) then
            if (jjj.eq.1) then
              jamin = 1
              jamax = 1
            else
              if (iandj.and.iamin.eq.1) go to 7900
              jamin = 2
              jamax = 4
            end if
          else
            jamin = jmin
            jamax = jmax
          end if
          jgjj=1
          do 112 icc=jpmin,jpmax
            if (jamin.eq.1) then
              coefj(jgjj)=cs(icc)
            else if (jamin.lt.5) then
              coefj(jgjj)=cp(icc)
            else if (jamin.lt.11) then
              coefj(jgjj)=cd(icc)
            else if (jamin.le.20) then
              coefj(jgjj)=cf(icc)
            else if (jamin.le.35) then
              coefj(jgjj)=cg(icc)
            end if
            jgjj = jgjj + 1
 112      continue
          npb=iang(jamax)
c     n+n' is the sum of the angular momentum plus 1 to index arrays
          npnp = npc + npb - 1
          ijmax = (iamina(itemp+2)-iamina(itemp+1))*(jamax-jamin+1)
          candb = icntr .eq. jcntr
          cx = c(1,icntr)
          cy = c(2,icntr)
          cz = c(3,icntr)
          bx = c(1,jcntr)
          by = c(2,jcntr)
          bz = c(3,jcntr)
          bmcx = bx-cx
          bmcy = by-cy
          bmcz = bz-cz
          bpcx = cx+bx
          bpcy = cy+by
          bpcz = cz+bz
          cbsq = bmcx*bmcx+bmcy*bmcy+bmcz*bmcz
c now loop over each center with an ecp potential
          do 7250 ikcntr=1,nat
           if (icntr.eq.ikcntr) go to 7250
c zero out the arrays which will collect individual integrals for
c later gradient formation
           do 100 i=1,ijmax
             g(i) = zero
  100        g2(i) = zero
           kcntr = ikcntr
           ax = c(1,kcntr)
           ay = c(2,kcntr)
           az = c(3,kcntr)
           if(lpskip(kcntr).eq.1) go to 7250
           llmx = lmax(kcntr)+1
           kf1 = kfirst(kcntr,1)
           kl1 = klast(kcntr,1)
           canda = icntr .eq. kcntr
           aandb = kcntr .eq. jcntr
           if (canda) then
             if (aandb) then
c   special case <a|a|a> only one center
               iicab = 1
             else
c   case <a|a|b>
               cax = zero
               cay = zero
               caz = zero
               ca = zero
               bax = bx - ax
               bay = by - ay
               baz = bz - az
               ba = sqrt(bax*bax+bay*bay+baz*baz)
               x = bax/ba
               y = bay/ba
               z = baz/ba
               iicab = 2
               iipow = 1
             end if
           else
             if (aandb) then
c   case <c|a|a>
               cax = cx - ax
               cay = cy - ay
               caz = cz - az
               ca = sqrt(cax*cax+cay*cay+caz*caz)
               x = cax/ca
               y = cay/ca
               z = caz/ca
               bax = zero
               bay = zero
               baz = zero
               ba = zero
               iicab = 3
               iipow = -1
             else
c   general case <c|a|b> three-center integral
c   actually c and b may still be equal
               cax = cx - ax
               cay = cy - ay
               caz = cz - az
               ca = sqrt(cax*cax+cay*cay+caz*caz)
               bax = bx - ax
               bay = by - ay
               baz = bz - az
               ba = sqrt(bax*bax+bay*bay+baz*baz)
               xca = cax/ca
               yca = cay/ca
               zca = caz/ca
               xba = bax/ba
               yba = bay/ba
               zba = baz/ba
               iicab = 4
               iipow = 0
             end if
           end if
c  set up tables of the powers of the cartesian distances (cax, cay ...)
c  for later use. pass in the maximum angular momentum for i and j
c  use max to index iang to make sure we get the max (l shells)
           if ((icntr.ne.kcntr).or.(kcntr.ne.jcntr)) then
             call ecppwr(iipow,iang(iamax)+1,iang(jamax))
           end if
c  calculate the integrals for the portion shifted up
           npnp = npnp+1
           npc = npc+1
           iamin = iamina(itemp+1)
           iamax = iamina(itemp+2)-1
           ijmax = (iamax-iamin+1)*(jamax-jamin+1)
           if (iicab.eq.1) then
             call ecpa11(g2,coefip,coefj,fpqr,npnp)
             if (llmx.gt.1) then
               call ecpa21(g2,dcoef4,coefip,coefj,npnp)
             end if
           else if (iicab.eq.2.or.iicab.eq.3) then
             call ecpr12(fp,coefip,coefj,iicab,npnp)
             call ecpa12(fp,jfst1,lbecp1,dcoef1,g2,
     *                   iicab,npnp,zlm,lmf,lmx,lmy,lmz)
             if (llmx.gt.1) then
               do 400 ll=2,llmx
                 call ecpr22(fp,coefip,coefj,iicab,npnp,ll)
                 call ecpa22(fp,jfst2,lbecp2,dcoef2,dcoef4,g2,
     *                    iicab,npnp,ll-2,zlm,lmf,lmx,lmy,lmz)
 400           continue
             end if
           else if (iicab.eq.4) then
             call ecpd14(fp,fp2,coefip,coefi,coefj,npnp,
     *          zlm,lmf,lmx,lmy,lmz)
             call ecpa14(fp,jfst1,lbecp1,dcoef1,g2,npnp)
             if (llmx.gt.1) then
               ntempp=max(npc+1,npb)
               do 500 ll=2,llmx
                 call ecpr24(fp,coefip,coefj,npnp,ntempp,ll)
                 call ecpa24(fp,jfst2,lbecp2,dcoef2,g2,
     *                       npnp,ll-2,zlm,lmf,lmx,lmy,lmz)
 500           continue
             end if
           end if
           npnp = npnp -1
           npc = npc-1
c  and now the part shifted down
           if (npc .gt. 1) then
             npnp = npnp - 1
             iamin = iamina(itemp-1)
             iamax = iamina(itemp)-1
             ijmax = (iamax-iamin+1)*(jamax-jamin+1)
             if (iicab.eq.1) then
                call ecpa11(g,coefi,coefj,fpqr,npnp)
                if (llmx.gt.1) then
                  call ecpa21(g,dcoef4,coefi,coefj,npnp)
                end if
             else if (iicab.eq.2.or.iicab.eq.3) then
               call ecpr12(fp,coefi,coefj,iicab,npnp)
               call ecpa12(fp,jfst1,lbecp1,dcoef1,g,
     *                     iicab,npnp,zlm,lmf,lmx,lmy,lmz)
               if (llmx.gt.1) then
                 do 200 ll=2,llmx
                   call ecpr22(fp,coefi,coefj,iicab,npnp,ll)
                   call ecpa22(fp,jfst2,lbecp2,dcoef2,dcoef4,g,
     *                    iicab,npnp,ll-2,zlm,lmf,lmx,lmy,lmz)
 200             continue
               end if
             else if (iicab.eq.4) then
               call ecpa14(fp2,jfst1,lbecp1,dcoef1,g,npnp)
               if (llmx.gt.1) then
                 ntempp=max(npc-1,npb)
                 do 300 ll=2,llmx
                   call ecpr24(fp,coefi,coefj,npnp,ntempp,ll)
                   call ecpa24(fp,jfst2,lbecp2,dcoef2,g,
     *                         npnp,ll-2,zlm,lmf,lmx,lmy,lmz)
 300             continue
               end if
             end if
             npnp = npnp + 1
c   reset ijmax for use below
             iamin = iamina(itemp+1)
             iamax = iamina(itemp+2)-1
             ijmax = (iamax-iamin+1)*(jamax-jamin+1)
           end if
           iamin = iamina(itemp)
           iamax = iamina(itemp+1)-1
c   form the derivitive integrals from the computed regular integrals
           call formdr(g,g2,xin,yin,zin,iamin,iamax,jamin,jamax,norm)
c   compute the gradient by multiplying the derivitive terms by the
c   corresponding density matrix elements
           n = 0
           do 7400 j=jamin,jamax
            jn = locj+j
            do 7400 i=iamin,iamax
             n = n+1
             in = loci+i
             nn = iky(in)+jn
             if(jn.gt.in) nn = iky(jn)+in
             dum = -datot(nn)*2.0d+00
c       next three lines produce the derivative on the basis function
c       center
             deloc(1,icntr) = deloc(1,icntr)+dum*xin(n)
             deloc(2,icntr) = deloc(2,icntr)+dum*yin(n)
             deloc(3,icntr) = deloc(3,icntr)+dum*zin(n)
c       next three lines produce the gradient due to the derivative
c       on the ecp center
             deloc(1,ikcntr) = deloc(1,ikcntr)-dum*xin(n)
             deloc(2,ikcntr) = deloc(2,ikcntr)-dum*yin(n)
             deloc(3,ikcntr) = deloc(3,ikcntr)-dum*zin(n)
c       the next block stores the 1st derivative elements of the
c       second derivative. 
             if (oham) then
               dum = 1.0d+00
               if (in.eq.jn) dum = 2.0d+00
               gout(nword  ) = - dum * xin(n)
               gout(nword+1) = - dum * yin(n)
               gout(nword+2) = - dum * zin(n)
               gout(nword+3) =   dum * xin(n)
               gout(nword+4) =   dum * yin(n)
               gout(nword+5) =   dum * zin(n)
               nword = nword + 6
               if (nword.ge.m5110) then
                nword = nword - 1
                call wrt3(gout,m5110,ipos1,num8)
                nword = 1
                ipos1 = ipos1 + 10
               end if
             end if
c
 7400      continue
c end of kcntr ecp potential center loop
 7250     continue
 7900    continue
 8000   continue
 8900  continue
 9000 continue
c
      if (oham) then
       nword = nword - 1
       call wrt3(gout,m5110,ipos1,num8)
       ipos1 = ipos1 + 10
       call clredx
      end if
c
c        now add to de
c
      do 9050 j=1,nat
        do 9050 i=1,3
 9050 de(i,j)=de(i,j)+deloc(i,j)
c
      if(outall) then
        write(iwr,9995)
        mmax = 0
 420    mmin = mmax + 1
        mmax = mmax + 8
        if (mmax .gt. nat) mmax = nat
        write (iwr,6020)
        write (iwr,6030) (j,j=mmin,mmax)
        write (iwr,6020)
        do i = 1 , 3
         write (iwr,6040) dnam(i) , (deloc(i,j),j=mmin,mmax)
        enddo
        if (mmax.lt.nat) go to 420
        write(iwr,9999)
      end if
      return
 6020 format (/)
 6030 format (5x,'atom',8(6x,i3,6x))
 6040 format (7x,a3,8e15.7)

 9995 format(/10x,30('*')/10x,'-ecp- contribution to gradient'/
     1        10x,30('*'))
 9999 format(/' ...... end of -ecp- gradient ..... ')
      end
**==ecp2d.f
      subroutine ecp2d(eg,eh,datot,l2,dcoef1,jfst1,lbecp1,dcoef4,
     *          dcoef2,jfst2,lbecp2,fpqr,gl1l1,gl2l,fp,glm2l,gll,
     *          glp2l,glm1lm1,glm1lp1,glp1lm1,glp1lp1,zlm,lmf,
     *          lmx,lmy,lmz,mapshl,nshels)
      implicit real*8 (a-h,o-z)
      logical dbug,canda,aandb
c
c  dimension passed in arrays
c  225 holds formed integrals gxg
c  441 holds intermediate integrals hxh (bigger than gxi)
c  note: some arrays could be smaller, but the code is simpler this way
c
      dimension datot(l2),jfst1(*),jfst2(*),
     *          dcoef1(*),lbecp1(9,*),dcoef2(*),lbecp2(6,*),dcoef4(*)
      dimension gl1l1(9,225), gl2l(9,225),
     *          fpqr(25,25,25),eh(9,*),eg(3,*)
      dimension glm2l(441), gll(441), glp2l(441), glm1lm1(441),
     *          glm1lp1(441), glp1lm1(441), glp1lp1(441)
      dimension zlm(*),lmf(*),lmx(*),lmy(*),lmz(*)
      dimension mapshl(nshels,*)
c  fp is used for radial integral storage and must be at least 2*11**3
      dimension fp(2662)
c
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
c
      character *8 zcom,ztitle,zaname,ztag,zsymm,zgroup
      character *10 zbflab
      character *8 zscftp,zruntp,zguess,zconf,zstate,zorb,zpseud
      common /runlab/ zcom(19),ztitle(10),zaname(maxat),zbflab(maxorb),
     +   ztag(maxat),zsymm,zgroup,zscftp,zruntp,zguess,zconf,zstate,
     +   zorb(maxorb),zpseud(maxat)
c
c....  zaname : names of atoms during calculation (after reorder)
c....  ztag   : names of atoms as read in 
c
c
c...ecp parameters
c
c
      real*8 zetc,cax,cay,caz,ca,xca,yca,zca
      real*8 zetb,bax,bay,baz,ba,xba,yba,zba
      real*8 phase,dax,day,daz,da,xda,yda,zda,xint
      integer kcntr
      common /ecp1/ zetc,cax,cay,caz,ca,xca,yca,zca,
     +              zetb,bax,bay,baz,ba,xba,yba,zba,
     +              phase,dax,day,daz,da,xda,yda,zda,xint,
     +              kcntr
c
      common /ecp2  / clp(400),zlp(400),nlp(400),kfirst(maxat,6),
     *                klast(maxat,6),lmax(maxat),lpskip(maxat),
     *                izcore(maxat)
c
      logical iandj,norm,normi,normj
      common /ecpidx/ q2,iamin,iamax,jamin,jamax,ipmin,ipmax,jpmin,
     *                jpmax,kf1,kl1,llmx,npc,npb,iandj,norm,normi,normj
c
      real*8 bmcx,bmcy,bmcz,bpcx,bpcy,bpcz,cbsq,ax,ay,az
      logical candb
      common /ecp4/ bmcx,bmcy,bmcz,bpcx,bpcy,bpcz,cbsq,ax,ay,az,
     +              candb
c
c   ecp common stuff
      common /zfncm / x,y,z
c
c  the following are to make use of symmetry
c
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
c
c  local storage
      dimension ci(mxprms), cj(mxprms), cip1(mxprms), cip2(mxprms),
     *          cjp1(mxprms),iang(84),iamina(8)
c
      parameter (zero=0.0d+00)
c
c    iang gives the maximum angular momentum for each shell
c
      data iang/1,3*2,6*3,10*4,15*5,21*6,28*7/
      data iamina/1,2,5,11,21,36,57,85/
c
c        -----  routine calculates the hessian of the ecp       -----
c        -----  integrals and produces an output vector.        -----
c        -----  routine based on new version of ecp1d, by       -----
c        -----  brett bode                                      -----
c
      dbug = nprint.eq. - 3 .or. nprint.eq. - 10
      norm=normf.ne.1.or.normp.ne.1
c  the shell being differentiated will be normalized when the derivative
c  is formed in formxx, so don't norm in the radial routines
      normi = .false.
      normj = .false.
      iandj = .false.
c
c        -----  loop over  ishell.                              -----
c        -----  note loops run over full square array           -----
c        -----  of ecp matrix for derivatives.                  -----
c
      do 9000 ii=1,nshell
c
       i1 = kstart(ii)
       i2 = i1+kng(ii)-1
       ipmin = i1
       ipmax = i2
       icntr = katom(ii)
       imin = kmin(ii)
       imax = kmax(ii)
       loci = kloc(ii)-imin
       iimax = 1
       if (imin.eq.1.and.imax.eq.4) iimax = 2
       do 8900 iii=1,iimax
        if (imin.eq.1.and.imax.eq.4) then
         if (iii.eq.1) then
           iamin = 1
           iamax = 1
         else
           iamin = 2
           iamax = 4
         end if
        else
         iamin = imin
         iamax = imax
        end if
        npc=iang(iamax)
c  store the coefs for later use. since we are taking the derivative
c  of the ishell we can combine the coef with the factor of 2ex
        igii=1
        do 111 jj=ipmin,ipmax
         if (iamin.eq.1) then
           ci(igii)=cs(jj)
           itype = 1
         else if (iamin.lt.5) then
           ci(igii)=cp(jj)
           itype = 2
         else if (iamin.lt.11) then
           ci(igii)=cd(jj)
           itype = 3
         else if (iamin.le.20) then
           ci(igii)=cf(jj)
           itype = 4
         else if (iamin.le.35) then
           ci(igii)=cg(jj)
           itype = 5
         end if
         cip1(igii) = ci(igii)*(-2.0d+00*ex(jj))
         cip2(igii) = cip1(igii)*(-2.0d+00*ex(jj))
         igii = igii + 1
 111    continue
c
c        -----  jshell  -----
c
        do 8000 jj=1,nshell
c    check symmetry
         n2=0
         ii0 = max(ii,jj)
         jj0 = min(ii,jj)
         do 80 it=1,nt
          id=mapshl(ii,it)
          jd=mapshl(jj,it)
          idd = max(id,jd)
          jdd = min(id,jd)
          if(idd.gt.ii0) go to 8000
          if(idd.lt.ii0) go to 80
          if(jdd.gt.jj0) go to 8000
          if(jdd.lt.jj0) go to 80
          n2=n2+1
 80      continue
         q2 = nt
         q2 = q2/n2
c
         j1 = kstart(jj)
         j2 = j1+kng(jj)-1
         jpmin = j1
         jpmax = j2
         jcntr = katom(jj)
         jmin = kmin(jj)
         jmax = kmax(jj)
         locj = kloc(jj)-jmin
         jjmax = 1
         if (jmin.eq.1.and.jmax.eq.4) jjmax = 2
         do 7900 jjj=1,jjmax
          if (jmin.eq.1.and.jmax.eq.4) then
            if (jjj.eq.1) then
              jamin = 1
              jamax = 1
            else
              if (iandj.and.iamin.eq.1) go to 7900
              jamin = 2
              jamax = 4
            end if
          else
            jamin = jmin
            jamax = jmax
          end if
          jgjj=1
          do 112 icc=jpmin,jpmax
            if (jamin.eq.1) then
              cj(jgjj)=cs(icc)
              jtype = 1
            else if (jamin.lt.5) then
              cj(jgjj)=cp(icc)
              jtype = 2
            else if (jamin.lt.11) then
              cj(jgjj)=cd(icc)
              jtype = 3
            else if (jamin.le.20) then
              cj(jgjj)=cf(icc)
              jtype = 4
            else if (jamin.le.35) then
              cj(jgjj)=cg(icc)
              jtype = 5
            end if
            cjp1(jgjj) = cj(jgjj)*(-2.0d+00*ex(icc))
            jgjj = jgjj + 1
 112      continue
          npb=iang(jamax)
c     n+n' is the sum of the angular momentum plus 1 to index arrays
          npnp = npc + npb - 1
          candb = icntr .eq. jcntr
          cx = c(1,icntr)
          cy = c(2,icntr)
          cz = c(3,icntr)
          bx = c(1,jcntr)
          by = c(2,jcntr)
          bz = c(3,jcntr)
          bmcx = bx-cx
          bmcy = by-cy
          bmcz = bz-cz
          bpcx = cx+bx
          bpcy = cy+by
          bpcz = cz+bz
          cbsq = bmcx*bmcx+bmcy*bmcy+bmcz*bmcz
c now loop over each center with an ecp potential
          do 7250 ikcntr=1,nat
           if ((icntr.eq.ikcntr).and.(jcntr.eq.ikcntr)) go to 7250
           ijmax = max((iamina(itype+2)-iamina(itype+1))*
     *            (iamina(jtype+2)-iamina(jtype+1)),
     *            (iamina(itype+3)-iamina(itype+2))*
     *            (iamina(jtype+1)-iamina(jtype)))
c zero out the arrays which will collect individual integrals for
c later hessian formation
           do 100 i=1,ijmax
             glm2l(i) = zero
             gll(i) = zero
             glp2l(i) = zero
             glm1lm1(i) = zero
             glm1lp1(i) = zero
             glp1lm1(i) = zero
  100        glp1lp1(i) = zero
           kcntr = ikcntr
           ax = c(1,kcntr)
           ay = c(2,kcntr)
           az = c(3,kcntr)
           if(lpskip(kcntr).eq.1) go to 7250
           llmx = lmax(kcntr)+1
           kf1 = kfirst(kcntr,1)
           kl1 = klast(kcntr,1)
           canda = icntr .eq. kcntr
           aandb = kcntr .eq. jcntr
           if (canda) then
             if (aandb) then
c   special case <a|a|a> only one center
               iicab = 1
             else
c   case <a|a|b>
               cax = zero
               cay = zero
               caz = zero
               ca = zero
               bax = bx - ax
               bay = by - ay
               baz = bz - az
               ba = sqrt(bax*bax+bay*bay+baz*baz)
               x = bax/ba
               y = bay/ba
               z = baz/ba
               iicab = 2
               iipow = 1
             end if
           else
             if (aandb) then
c   case <c|a|a>
               cax = cx - ax
               cay = cy - ay
               caz = cz - az
               ca = sqrt(cax*cax+cay*cay+caz*caz)
               x = cax/ca
               y = cay/ca
               z = caz/ca
               bax = zero
               bay = zero
               baz = zero
               ba = zero
               iicab = 3
               iipow = -1
             else
c   general case <c|a|b> three-center integral
c   actually c and b may still be equal
               cax = cx - ax
               cay = cy - ay
               caz = cz - az
               ca = sqrt(cax*cax+cay*cay+caz*caz)
               bax = bx - ax
               bay = by - ay
               baz = bz - az
               ba = sqrt(bax*bax+bay*bay+baz*baz)
               xca = cax/ca
               yca = cay/ca
               zca = caz/ca
               xba = bax/ba
               yba = bay/ba
               zba = baz/ba
               iicab = 4
               iipow = 0
             end if
           end if
c  set up tables of the powers of the cartesian distances (cax, cay ...)
c  for later use. pass in the maximum angular momentum for i and j
c  use max to index iang to make sure we get the max (l shells)
           if ((icntr.ne.kcntr).or.(kcntr.ne.jcntr)) then
             call ecppwr(iipow,iang(iamax)+2,iang(jamax)+1)
           end if
           if (icntr.ne.ikcntr) then
            normj = norm
c  calculate the <l|u|l> term
            if (iicab.eq.1) then
              call ecpa11(gll,cip1,cj,fpqr,npnp)
              if (llmx.gt.1) then
                call ecpa21(gll,dcoef4,cip1,cj,npnp)
              end if
            else if (iicab.eq.2.or.iicab.eq.3) then
              call ecpr12(fp,cip1,cj,iicab,npnp)
              call ecpa12(fp,jfst1,lbecp1,dcoef1,gll,
     *                    iicab,npnp,zlm,lmf,lmx,lmy,lmz)
              if (llmx.gt.1) then
                do 200 ll=2,llmx
                  call ecpr22(fp,cip1,cj,iicab,npnp,ll)
                  call ecpa22(fp,jfst2,lbecp2,dcoef2,dcoef4,gll,
     *                        iicab,npnp,ll-2,zlm,lmf,lmx,lmy,lmz)
 200            continue
              end if
            else if (iicab.eq.4) then
              call ecpr14(fp,cip1,cj,npnp,zlm,lmf,lmx,lmy,lmz)
              call ecpa14(fp,jfst1,lbecp1,dcoef1,gll,npnp)
              if (llmx.gt.1) then
                ntempp=max(npc,npb)
                do 210 ll=2,llmx
                  call ecpr24(fp,cip1,cj,npnp,ntempp,ll)
                  call ecpa24(fp,jfst2,lbecp2,dcoef2,gll,npnp,ll-2,
     *                        zlm,lmf,lmx,lmy,lmz)
 210            continue
              end if
            end if
c  calculate the <l+2|u|l> term
            npnp = npnp+2
            npc = npc+2
            iamin = iamina(itype+2)
            iamax = iamina(itype+3)-1
            ijmax = (iamax-iamin+1)*(jamax-jamin+1)
            if (iicab.eq.1) then
              call ecpa11(glp2l,cip2,cj,fpqr,npnp)
              if (llmx.gt.1) then
                call ecpa21(glp2l,dcoef4,cip2,cj,npnp)
              end if
            else if (iicab.eq.2.or.iicab.eq.3) then
              call ecpr12(fp,cip2,cj,iicab,npnp)
              call ecpa12(fp,jfst1,lbecp1,dcoef1,glp2l,
     *                    iicab,npnp,zlm,lmf,lmx,lmy,lmz)
              if (llmx.gt.1) then
                do 220 ll=2,llmx
                  call ecpr22(fp,cip2,cj,iicab,npnp,ll)
                  call ecpa22(fp,jfst2,lbecp2,dcoef2,dcoef4,glp2l,
     *                        iicab,npnp,ll-2,zlm,lmf,lmx,lmy,lmz)
 220            continue
              end if
            else if (iicab.eq.4) then
              call ecpr14(fp,cip2,cj,npnp,zlm,lmf,lmx,lmy,lmz)
              call ecpa14(fp,jfst1,lbecp1,dcoef1,glp2l,npnp)
              if (llmx.gt.1) then
                ntempp=max(npc,npb)
                do 230 ll=2,llmx
                  call ecpr24(fp,cip2,cj,npnp,ntempp,ll)
                  call ecpa24(fp,jfst2,lbecp2,dcoef2,glp2l,npnp,ll-2,
     *                        zlm,lmf,lmx,lmy,lmz)
 230            continue
              end if
            end if
c  calculate the <l-2|u|l> term
            npnp = npnp-4
            npc = npc-4
            if (npc .ge. 1) then
              iamin = iamina(itype-2)
              iamax = iamina(itype-1)-1
              ijmax = (iamax-iamin+1)*(jamax-jamin+1)
              if (iicab.eq.1) then
                call ecpa11(glm2l,ci,cj,fpqr,npnp)
                if (llmx.gt.1) then
                  call ecpa21(glm2l,dcoef4,ci,cj,npnp)
                end if
              else if (iicab.eq.2.or.iicab.eq.3) then
                call ecpr12(fp,ci,cj,iicab,npnp)
                call ecpa12(fp,jfst1,lbecp1,dcoef1,glm2l,
     *                      iicab,npnp,zlm,lmf,lmx,lmy,lmz)
                if (llmx.gt.1) then
                  do 240 ll=2,llmx
                    call ecpr22(fp,ci,cj,iicab,npnp,ll)
                    call ecpa22(fp,jfst2,lbecp2,dcoef2,dcoef4,glm2l,
     *                          iicab,npnp,ll-2,zlm,lmf,lmx,lmy,lmz)
 240              continue
                end if
              else if (iicab.eq.4) then
                call ecpr14(fp,ci,cj,npnp,zlm,lmf,lmx,lmy,lmz)
                call ecpa14(fp,jfst1,lbecp1,dcoef1,glm2l,npnp)
                if (llmx.gt.1) then
                  ntempp=max(npc,npb)
                  do 250 ll=2,llmx
                    call ecpr24(fp,ci,cj,npnp,ntempp,ll)
                    call ecpa24(fp,jfst2,lbecp2,dcoef2,glm2l,npnp,ll-2,
     *                          zlm,lmf,lmx,lmy,lmz)
 250              continue
                end if
              end if
            end if
            npnp = npnp + 2
            npc = npc+2
           end if
c  turn of j normalization since it is now also derivatized
           normj = .false.
c  calculate the <l+1|u|l+1> term
           npnp = npnp+2
           npc = npc+1
           npb = npb+1
           iamin = iamina(itype+1)
           iamax = iamina(itype+2)-1
           jamin = iamina(jtype+1)
           jamax = iamina(jtype+2)-1
           ijmax = (iamax-iamin+1)*(jamax-jamin+1)
           if (iicab.eq.1) then
             call ecpa11(glp1lp1,cip1,cjp1,fpqr,npnp)
             if (llmx.gt.1) then
               call ecpa21(glp1lp1,dcoef4,cip1,cjp1,npnp)
             end if
           else if (iicab.eq.2.or.iicab.eq.3) then
             call ecpr12(fp,cip1,cjp1,iicab,npnp)
             call ecpa12(fp,jfst1,lbecp1,dcoef1,glp1lp1,
     *                   iicab,npnp,zlm,lmf,lmx,lmy,lmz)
             if (llmx.gt.1) then
               do 260 ll=2,llmx
                 call ecpr22(fp,cip1,cjp1,iicab,npnp,ll)
                 call ecpa22(fp,jfst2,lbecp2,dcoef2,dcoef4,glp1lp1,
     *                       iicab,npnp,ll-2,zlm,lmf,lmx,lmy,lmz)
 260           continue
             end if
           else if (iicab.eq.4) then
             call ecpr14(fp,cip1,cjp1,npnp,zlm,lmf,lmx,lmy,lmz)
             call ecpa14(fp,jfst1,lbecp1,dcoef1,glp1lp1,npnp)
             if (llmx.gt.1) then
               ntempp=max(npc,npb)
               do 270 ll=2,llmx
                 call ecpr24(fp,cip1,cjp1,npnp,ntempp,ll)
                 call ecpa24(fp,jfst2,lbecp2,dcoef2,glp1lp1,npnp,ll-2,
     *                       zlm,lmf,lmx,lmy,lmz)
 270           continue
             end if
           end if
           npnp = npnp - 2
           npc = npc-1
           npb = npb-1
c  calculate the <l+1|u|l-1> term
           npc = npc+1
           npb = npb-1
           if (npb .ge. 1) then
             iamin = iamina(itype+1)
             iamax = iamina(itype+2)-1
             jamin = iamina(jtype-1)
             jamax = iamina(jtype)-1
             ijmax = (iamax-iamin+1)*(jamax-jamin+1)
             if (iicab.eq.1) then
               call ecpa11(glp1lm1,cip1,cj,fpqr,npnp)
               if (llmx.gt.1) then
                 call ecpa21(glp1lm1,dcoef4,cip1,cj,npnp)
               end if
             else if (iicab.eq.2.or.iicab.eq.3) then
               call ecpr12(fp,cip1,cj,iicab,npnp)
               call ecpa12(fp,jfst1,lbecp1,dcoef1,glp1lm1,
     *                     iicab,npnp,zlm,lmf,lmx,lmy,lmz)
               if (llmx.gt.1) then
                 do 280 ll=2,llmx
                   call ecpr22(fp,cip1,cj,iicab,npnp,ll)
                   call ecpa22(fp,jfst2,lbecp2,dcoef2,dcoef4,glp1lm1,
     *                         iicab,npnp,ll-2,zlm,lmf,lmx,lmy,lmz)
 280             continue
               end if
             else if (iicab.eq.4) then
               call ecpr14(fp,cip1,cj,npnp,zlm,lmf,lmx,lmy,lmz)
               call ecpa14(fp,jfst1,lbecp1,dcoef1,glp1lm1,npnp)
               if (llmx.gt.1) then
                 ntempp=max(npc,npb)
                 do 290 ll=2,llmx
                  call ecpr24(fp,cip1,cj,npnp,ntempp,ll)
                  call ecpa24(fp,jfst2,lbecp2,dcoef2,glp1lm1,npnp,ll-2,
     *                        zlm,lmf,lmx,lmy,lmz)
 290             continue
               end if
             end if
           end if
           npc = npc-1
           npb = npb+1
c  calculate the <l-1|u|l+1> term
           npc = npc-1
           npb = npb+1
           if (npc .ge. 1) then
             iamin = iamina(itype-1)
             iamax = iamina(itype)-1
             jamin = iamina(jtype+1)
             jamax = iamina(jtype+2)-1
             ijmax = (iamax-iamin+1)*(jamax-jamin+1)
             if (iicab.eq.1) then
               call ecpa11(glm1lp1,ci,cjp1,fpqr,npnp)
               if (llmx.gt.1) then
                 call ecpa21(glm1lp1,dcoef4,ci,cjp1,npnp)
               end if
             else if (iicab.eq.2.or.iicab.eq.3) then
               call ecpr12(fp,ci,cjp1,iicab,npnp)
               call ecpa12(fp,jfst1,lbecp1,dcoef1,glm1lp1,
     *                     iicab,npnp,zlm,lmf,lmx,lmy,lmz)
               if (llmx.gt.1) then
                 do 300 ll=2,llmx
                   call ecpr22(fp,ci,cjp1,iicab,npnp,ll)
                   call ecpa22(fp,jfst2,lbecp2,dcoef2,dcoef4,glm1lp1,
     *                         iicab,npnp,ll-2,zlm,lmf,lmx,lmy,lmz)
 300             continue
               end if
             else if (iicab.eq.4) then
               call ecpr14(fp,ci,cjp1,npnp,zlm,lmf,lmx,lmy,lmz)
               call ecpa14(fp,jfst1,lbecp1,dcoef1,glm1lp1,npnp)
               if (llmx.gt.1) then
                 ntempp=max(npc,npb)
                 do 310 ll=2,llmx
                  call ecpr24(fp,ci,cjp1,npnp,ntempp,ll)
                  call ecpa24(fp,jfst2,lbecp2,dcoef2,glm1lp1,npnp,ll-2,
     *                        zlm,lmf,lmx,lmy,lmz)
 310             continue
               end if
             end if
           end if
           npc = npc+1
           npb = npb-1
c  calculate the <l-1|u|l-1> term
           npc = npc-1
           npb = npb-1
           if ((npc .ge. 1).and.(npb .ge. 1)) then
             npnp = npnp - 2
             iamin = iamina(itype-1)
             iamax = iamina(itype)-1
             jamin = iamina(jtype-1)
             jamax = iamina(jtype)-1
             ijmax = (iamax-iamin+1)*(jamax-jamin+1)
             if (iicab.eq.1) then
               call ecpa11(glm1lm1,ci,cj,fpqr,npnp)
               if (llmx.gt.1) then
                 call ecpa21(glm1lm1,dcoef4,ci,cj,npnp)
               end if
             else if (iicab.eq.2.or.iicab.eq.3) then
               call ecpr12(fp,ci,cj,iicab,npnp)
               call ecpa12(fp,jfst1,lbecp1,dcoef1,glm1lm1,
     *                     iicab,npnp,zlm,lmf,lmx,lmy,lmz)
               if (llmx.gt.1) then
                 do 320 ll=2,llmx
                   call ecpr22(fp,ci,cj,iicab,npnp,ll)
                   call ecpa22(fp,jfst2,lbecp2,dcoef2,dcoef4,glm1lm1,
     *                         iicab,npnp,ll-2,zlm,lmf,lmx,lmy,lmz)
 320             continue
               end if
             else if (iicab.eq.4) then
               call ecpr14(fp,ci,cj,npnp,zlm,lmf,lmx,lmy,lmz)
               call ecpa14(fp,jfst1,lbecp1,dcoef1,glm1lm1,npnp)
               if (llmx.gt.1) then
                 ntempp=max(npc,npb)
                 do 330 ll=2,llmx
                  call ecpr24(fp,ci,cj,npnp,ntempp,ll)
                  call ecpa24(fp,jfst2,lbecp2,dcoef2,glm1lm1,npnp,ll-2,
     *                        zlm,lmf,lmx,lmy,lmz)
 330             continue
               end if
             end if
             npnp = npnp + 2
           end if
           npc = npc+1
           npb = npb+1
           iamin = iamina(itype)
           iamax = iamina(itype+1)-1
           jamin = iamina(jtype)
           jamax = iamina(jtype+1)-1
           ijmax = (iamax-iamin+1)*(jamax-jamin+1)
c   form the derivitive integrals from the computed regular integrals
           call formii(glp1lp1,glp1lm1,glm1lp1,glm1lm1,gl1l1,
     *          iamin,iamax,jamin,jamax,norm)
           call formij(glp2l,gll,glm2l,gl2l,iamin,iamax,jamin,jamax,
     *         norm)

c   compute the gradient by multiplying the derivitive terms by the
c   corresponding density matrix elements
           n = 0
           do 7400 j=jamin,jamax
            jn = locj+j
            do 7400 i=iamin,iamax
             n = n+1
             in = loci+i
             nn = iky(in)+jn
             if(jn.gt.in) nn = iky(jn)+in
             dum = datot(nn)*2.0d+00
             do 500 ll=1,9
               if (icntr.ne.kcntr) then
c   2nd derivative of the basis function
                 ncntr = icntr*(icntr-1)/2 + icntr
                 eh(ll,ncntr) = eh(ll,ncntr)+dum*gl2l(ll,n)
c   2nd derivative of the ecp center
                 ncntr = kcntr*(kcntr-1)/2 + kcntr
                 eh(ll,ncntr) = eh(ll,ncntr)+dum*gl2l(ll,n)
c   one derivative on the ecp center
c   there are two such terms which are equal when i=k,
c   and related by exchange of the order of differentiation
c   otherwise. thus only one off diagonal term is stored
                 if (kcntr.eq.icntr) then
                   ncntr = icntr*(icntr-1)/2 + icntr
                   eh(ll,ncntr) = eh(ll,ncntr)-2*dum*gl2l(ll,n)
                 else if (icntr.lt.kcntr) then
                   ncntr = kcntr*(kcntr-1)/2 + icntr
                   eh(ll,ncntr) = eh(ll,ncntr)-dum*gl2l(ll,n)
                 else
                   ncntr = icntr*(icntr-1)/2 + kcntr
                   eh(ll,ncntr) = eh(ll,ncntr)-dum*gl2l(ll,n)
                 end if
               end if
c                term from the two first derivatives on ecp center
               ncntr = kcntr*(kcntr-1)/2 + kcntr
               eh(ll,ncntr) = eh(ll,ncntr)+dum*gl1l1(ll,n)
c                derivatives on each basis function center
               if (jcntr.le.icntr) then
                 ncntr = icntr*(icntr-1)/2 + jcntr
                 eh(ll,ncntr) = eh(ll,ncntr)+dum*gl1l1(ll,n)
               end if
c                first derivative terms with one derivative on the
c                ecp center
               if (kcntr.le.icntr) then
                 ncntr = icntr*(icntr-1)/2 + kcntr
                 eh(ll,ncntr) = eh(ll,ncntr)-dum*gl1l1(ll,n)
               end if
               if (jcntr.le.kcntr) then
                 ncntr = kcntr*(kcntr-1)/2 + jcntr
                 eh(ll,ncntr) = eh(ll,ncntr)-dum*gl1l1(ll,n)
               end if
 500         continue
 7400      continue
c end of kcntr ecp potential center loop
 7250     continue
 7900    continue
 8000   continue
 8900  continue
 9000 continue
c
c
c   do not symmetrize the ecp contribution to the hessian
c
c     call symeh(eh)
c
      if(dbug) then
         write(iwr,9995)
         call hssprt(nat,eg,eh)
         write(iwr,9999)
      end if
      return
 9995 format(/10x,30('-')/10x,'-ecp- contribution to hessian'/
     1        10x,30('-'))
 9999 format(/,' ...... end of -ecp- hessian ..... ')
      end
**==ecdint.f
      subroutine ecdint(x,datot,iso,nshels,oham,out)
      implicit real*8 (a-h,o-z)
      logical out,oham
      dimension datot(*), x(*), iso(nshels,*)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer ncoef1,ncoef2,j1len,j2len,llim,nlim,ntlim,j4len
      common /ecpdim/ ncoef1,ncoef2,j1len,j2len,llim,nlim,ntlim,j4len
c
      common/junk/de(3,maxat)
      data m18/18/
c
      l2 = (num*num+num)/2
      t0 = cpulft(1)
c
c  simple driver to allocate memory and call the routine to
c  compute the ecp modifications to the gradient
c
      nav = lenwrd()
c
      ldcf1  = 0
      ljln   = ldcf1  +  ncoef1
      llb1   = ljln   + (j1len-1)/nav+1
      ldcf4  = llb1   + (9*ncoef1-1)/nav+1
      ldcf2  = ldcf4  +  j4len
      lj2n   = ldcf2  +  ncoef2
      llb2   = lj2n   + (j2len-1)/nav+1
      ldeecp = llb2   + (6*ncoef2)/nav
      lfpqr  = ldeecp + 3*nat
      lfp    = lfpqr  + 15625
      lfp2   = lfp    + 2662
      lgg    = lfp2   + 2662
      lgg2   = lgg    + 315
      lxin   = lgg2   + 315
      lyin   = lxin   + 225
      lzin   = lyin   + 225
      lzlm   = lzin   + 225
      llmf   = lzlm   + 581
      llmx   = llmf   + 122/nav
      llmy   = llmx   + 582/nav
      llmz   = llmy   + 582/nav
      last   = llmz   + 582/nav
c
      ldcf1 = igmem_alloc(last)
c
      ljln   = ldcf1  +  ncoef1
      llb1   = ljln   + (j1len-1)/nav+1
      ldcf4  = llb1   + (9*ncoef1-1)/nav+1
      ldcf2  = ldcf4  +  j4len
      lj2n   = ldcf2  +  ncoef2
      llb2   = lj2n   + (j2len-1)/nav+1
      ldeecp = llb2   + (6*ncoef2)/nav
      lfpqr  = ldeecp + 3*nat
      lfp    = lfpqr  + 15625
      lfp2   = lfp    + 2662
      lgg    = lfp2   + 2662
      lgg2   = lgg    + 315
      lxin   = lgg2   + 315
      lyin   = lxin   + 225
      lzin   = lyin   + 225
      lzlm   = lzin   + 225
      llmf   = lzlm   + 581
      llmx   = llmf   + 122/nav
      llmy   = llmx   + 582/nav
      llmz   = llmy   + 582/nav
c
c     -----  read in ecp formulas and data  -----
c
        ldaf91 = (j1len-1)/nav+1 + (9*ncoef1-1)/nav+1
        ldaf93 = (j2len-1)/nav+1 + (6*ncoef2)/nav
c
         call secget(isect(479),m18,ibl479)
c
         call rdedx(x(ldcf1),ncoef1,ibl479,idaf)
         call reads(x(ljln),ldaf91,idaf)
         call reads(x(ldcf2),ncoef2,idaf)
         call reads(x(lj2n),ldaf93,idaf)
c
c the next calls just fill out tables which have probably already
c been done in the energy step, but just in case this is a restart run
c
        call dawt()
        call errt()
        call dawert()
        call ecpini(x(llmf),x(llmx),x(llmy),x(llmz))
        call ztab(x(lzlm))
        call ftab(x(lfpqr), nlim-1)
        call eccod3(x(lfpqr),x(ldcf4),x(lzlm),x(llmf),x(llmx),x(llmy),
     *              x(llmz))
c
c now compute the gradient
c
         call ecp1d(de,datot,l2,x(ldcf1),x(ljln),x(llb1),x(ldcf4),
     *     x(ldcf2),x(lj2n),x(llb2),x(ldeecp),x(lfpqr),x(lfp),x(lfp2),
     *     x(lgg),x(lgg2),x(lxin),x(lyin),x(lzin),x(lzlm),x(llmf),
     *     x(llmx),x(llmy),x(llmz),iso,nshels)
c
        if(out) then
           t1=cpulft(1)
           tg = t1-t0
           write(iwr,9030) tg
        end if
c
      call gmem_free(ldcf1)
c
      return
 9030 format(1x,'time to compute ecp gradient integrals =',f10.2)
      end
**==dr2ecp.f
      subroutine dr2ecp(dd,isec46)
c------------------------------------------------------------------
c     second derivatives of ECP
c------------------------------------------------------------------
      implicit real*8  (a-h,o-z)
      dimension dd(*)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
      real*8 rcigrd
      integer isecdd, isecll, ifil2d, iblk2d, iword, mnnr, mnc
      integer mnv, mnx, iscigr, isecmo
      integer isecnd, isecsy, irlagr, iadfrc, nfc, intlgr
      integer ncepa, ispaer
      integer nd2mo, ncore, ncact, nvr, ifilh, iblkh, iblk1
      integer ibl222, ntot, nupact, ijr3
      logical cigr, cicv, mpgr, mcgr, cicx, umpgr
      logical lcisd, lcepa, lacpf, lnewci, lsingl
      common /cigrad/ cigr,isecdd,isecll,ifil2d,iblk2d,iword,cicv,
     +                mnnr,mnc,mnv,mnx,mpgr,mcgr,cicx,iscigr,isecmo,
     +                isecnd,isecsy,irlagr,iadfrc,nfc,intlgr,umpgr,
     +                lcisd,lcepa,lacpf,lnewci,lsingl,ncepa,ispaer(20),
     +                nd2mo,ncore,ncact,nvr,ifilh,iblkh,iblk1,
     +                ibl222,ntot,nupact,ijr3,rcigrd(70)
c
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      real*8 pito52, pidiv4, root3, root5, root53, root7
      common /picon/ pito52,pidiv4,root3,root5,root53,root7
c
c
      logical omem, ofirt
      integer maxq, maxqq, iqqoff, mxblock
c
      common/maxlen/maxq,maxqq,omem(maxlfn),iqqoff(maxlfn),
     +              ofirt,mxblock
c
      logical out
c
c     second derivatives of ECP.
c
      L2 = (NUM*NUM+NUM)/2
      LEG  = 1
      nat3 = 3 * nat
      LEH  = LEG    + nat3
      LDA  = LEH    + 9*(NAT*NAT+NAT)/2
      LDB  = LDA    + L2
      ISO  = LDB    + 3*NAT*L2
      LAST = ISO    + lenrel(nw196(5))
c
      negh=3*nat+9*(nat*(nat+1))/2
      call vclr(dd(leg),1,negh)
c
      call rdedx(dd(iso),nw196(5),ibl196(5),idaf)
      call onepdm(dd(lda),dd(ldb))
      out = odebug(6) .or. odebug(7)
c
      call ecphes(dd(last),dd(leg),dd(leh),dd(lda),
     =            dd(iso), nshell,out)
c
c     add to isec46
      nlen = 9 * nat * nat
      call rdedx(dd(ldb),nlen,isec46,idaf)
      call addecp(dd(ldb),dd(leh),nat,nat3)
c
      call wrt3(dd(ldb),nlen,isec46,idaf)
c
      if (out) then
       write (iwr,6020)
       call prnder(dd(ldb),nat3,iwr)
      endif
      return
 6020 format (/1x,'including contribution from second derivatives of',
     +        ' ecp matrix'//)
      end
      subroutine addecp(orig,ecphes,nat,nat3)
      implicit real*8 (a-h,o-z)
      dimension orig(nat3,nat3), ecphes(9,*)
c
      N = 0
      DO IATOM = 1,NAT
      iat3 = (iatom-1)* 3
        DO JATOM = 1,IATOM
           jat3 = (jatom-1)*3
           N = N + 1
           ij = 0
            do i = 1,3
            do j = 1,3
            ij = ij + 1
            orig(iat3+i,jat3+j) = ecphes(ij,n) + orig(iat3+i,jat3+j)
            enddo
            enddo
            if(iatom.ne.jatom) then
            ij = 0
            do i = 1,3
            do j = 1,3
            ij = ij + 1
            orig(jat3+j,iat3+i) = ecphes(ij,n) + orig(jat3+j,iat3+i)
            enddo
            enddo
            endif
          enddo
         enddo
c
       return
      end
**==ecphes.f
      subroutine ecphes(x,eg,eh,dab,iso,nshels,out)
      implicit real*8 (a-h,o-z)
c
      logical out
      dimension eg(3,*), eh(9,*), dab(*), x(*)
      dimension iso(nshels,*)
c
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
c
      character *8 zcom,ztitle,zaname,ztag,zsymm,zgroup
      character *10 zbflab
      character *8 zscftp,zruntp,zguess,zconf,zstate,zorb,zpseud
      common /runlab/ zcom(19),ztitle(10),zaname(maxat),zbflab(maxorb),
     +   ztag(maxat),zsymm,zgroup,zscftp,zruntp,zguess,zconf,zstate,
     +   zorb(maxorb),zpseud(maxat)
c
c....  zaname : names of atoms during calculation (after reorder)
c....  ztag   : names of atoms as read in 
c
c
      integer ncoef1,ncoef2,j1len,j2len,llim,nlim,ntlim,j4len
      common /ecpdim/ ncoef1,ncoef2,j1len,j2len,llim,nlim,ntlim,j4len
c
c
      character*8 rhf, uhf
      data rhf,uhf/'rhf' ,'uhf'    /
      data m18/18/
c
c     out = nprint.eq.-3. or. nprint.eq.-10
      l2 = (num*num+num)/2
      t0 = cpulft(1)
c
c  simple driver to allocate memory and call the routines to
c  compute the gradient and hessian ecp modifications
c
      nav = lenwrd()
c
      ldcf1 = 1
      ljln   = ldcf1  +  ncoef1
      llb1   = ljln   + (j1len-1)/nav+1
      ldcf4  = llb1   + (9*ncoef1-1)/nav+1
      ldcf2  = ldcf4  +  j4len
      lj2n   = ldcf2  +  ncoef2
      llb2   = lj2n   + (j2len-1)/nav+1
      lfpqr  = llb2   + (6*ncoef2)/nav
      lfp    = lfpqr  + 15625
      lfp2   = lfp    + 2662
      lzlm   = lfp2   + 2662
      llmf   = lzlm   + 581
      llmx   = llmf   + 122/nav
      llmy   = llmx   + 582/nav
      llmz   = llmy   + 582/nav
c
c this array is used only by the gradient portion
c
      lgg    = llmz   + 582/nav
      lgg2   = lgg    + 315
      lxin   = lgg2   + 315
      lyin   = lxin   + 225
      lzin   = lyin   + 225
      ldeecp = lzin   + 225
c
c these are used only by the hessian code and thus overwrite
c the gradient only terms
c
      lgl11 = llmz   + 582/nav
      lgl2  = lgl11  + 2025
      lgm2  = lgl2   + 2025
      lgll  = lgm2   + 441
      lglp2 = lgll   + 441
      lgm1m1= lglp2  + 441
      lgm1p1= lgm1m1 + 441
      lgp1m1= lgm1p1 + 441
      lgp1p1= lgp1m1 + 441
c
c     -----  read in ecp formulas and data  -----
c
        ldaf91 = (j1len-1)/nav+1 + (9*ncoef1-1)/nav+1
        ldaf93 = (j2len-1)/nav+1 + (6*ncoef2)/nav
         call secget(isect(479),m18,ibl479)
c
         call rdedx(x(ldcf1),ncoef1,ibl479,idaf)
         call reads(x(ljln),ldaf91,idaf)
         call reads(x(ldcf2),ncoef2,idaf)
         call reads(x(lj2n),ldaf93,idaf)
c
c the next calls just fill out tables which have probably already
c been done in the energy step, but just in case this is a restart run
c
        call dawt()
        call errt()
        call dawert()
        call ecpini(x(llmf),x(llmx),x(llmy),x(llmz))
        call ztab(x(lzlm))
        call ftab(x(lfpqr), nlim-1)
        call eccod3(x(lfpqr),x(ldcf4),x(lzlm),x(llmf),x(llmx),x(llmy),
     *              x(llmz))
c
c compute the 2nd derivative terms
c
        call ecp2d(eg,eh,dab,l2,x(ldcf1),x(ljln),x(llb1),x(ldcf4),
     *    x(ldcf2),x(lj2n),x(llb2),x(lfpqr),x(lgl11),x(lgl2),
     *    x(lfp),x(lgm2),x(lgll),x(lglp2),x(lgm1m1),x(lgm1p1),
     *    x(lgp1m1),x(lgp1p1),x(lzlm),x(llmf),x(llmx),x(llmy),x(llmz),
     *    iso,nshels)
        if(out) then
           t1 = cpulft(1)
           tg = t1-t0
           write(iwr,9010) tg
           t0 = t1
        end if
c
c        symmetrize the fock derivatives (single determinant jobs
c        will catch this later, after 2 electron contributions)
c
      if(zscftp.eq.rhf  .or.  zscftp.eq.uhf) return
c
      return
c
 9010 format(1x,'time to do ecp 2nd derivative integrals=',f10.2)
      end
**==formdr.f
      subroutine formdr(g,g2,xin,yin,zin,mini,maxi,minj,maxj,norm)
      implicit real*8 (a-h,o-z)
      logical norm
      dimension g(*), g2(*), xin(*), yin(*), zin(*)
      parameter (sqrt3=1.73205080756888d+00,sqrt5=2.23606797749979d+00)
      parameter (sqrt7=2.64575131106459d+00)
c
c   form the first derivative integrals
c   this routine is completly general and good through g'
c
      nn=0
      nj = maxj-minj+1
c   s'
        if (mini.eq.1) then
          do 100 j=1,nj
            nn=nn+1
            xin(nn) =  g2(j)
            yin(nn) =  g2(j+nj)
            zin(nn) =  g2(j+2*nj)
 100      continue
          return
        end if
c   p'
        if (mini.le.2) then
          do 200 j=1,nj
            nn = nn+1
            xin(nn) = (g2(j)+g(j))
            yin(nn) =  g2(j+3*nj)
            zin(nn) =  g2(j+4*nj)
            nn = nn+1
            xin(nn) =  g2(j+3*nj)
            yin(nn) = (g2(j+nj)+g(j))
            zin(nn) =  g2(j+5*nj)
            nn = nn+1
            xin(nn) =  g2(j+4*nj)
            yin(nn) =  g2(j+5*nj)
            zin(nn) =  g2(j+2*nj)+g(j)
 200      continue
          return
        end if
c   d'
        if (mini.le.5) then
          do 300 j=1,nj
            nn = nn+1
            xin(nn) = (g2(j)+g(j)+g(j))
            yin(nn) =  g2(j+3*nj)
            zin(nn) =  g2(j+4*nj)
            nn = nn+1
            xin(nn) =  g2(j+5*nj)
            yin(nn) = (g2(j+nj)+g(j+nj)+g(j+nj))
            zin(nn) =  g2(j+6*nj)
            nn = nn+1
            xin(nn) =  g2(j+7*nj)
            yin(nn) =  g2(j+8*nj)
            zin(nn) = (g2(j+2*nj)+g(j+2*nj)+g(j+2*nj))
            nn = nn+1
            dum=1.0d+00
            if(norm) dum=sqrt3
            xin(nn) = dum*(g2(j+3*nj)+g(j+nj))
            yin(nn) = dum*(g2(j+5*nj)+g(j))
            zin(nn) = dum* g2(j+9*nj)
            nn = nn+1
            xin(nn) = dum*(g2(j+4*nj)+g(j+2*nj))
            yin(nn) = dum* g2(j+9*nj)
            zin(nn) = dum*(g2(j+7*nj)+g(j))
            nn = nn+1
            xin(nn) = dum* g2(j+9*nj)
            yin(nn) = dum*(g2(j+6*nj)+g(j+2*nj))
            zin(nn) = dum*(g2(j+8*nj)+g(j+nj))
 300      continue
          return
        end if
c   f'
        if (mini.le.11) then
          do 400 j=1,nj
            nn=nn+1
            xin(nn)=(g2(j)+3.0d+00*g(j))
            yin(nn)= g2(j+3*nj)
            zin(nn)= g2(j+4*nj)
            nn=nn+1
            xin(nn)= g2(j+5*nj)
            yin(nn)=(g2(j+nj)+3.0d+00*g(j+nj))
            zin(nn)= g(j+6*nj)
            nn=nn+1
            xin(nn)= g2(j+7*nj)
            yin(nn)= g2(j+8*nj)
            zin(nn)=(g2(j+2*nj)+3.0d+00*g(j+2*nj))
            nn=nn+1
            dum = 1.0d+00
            if(norm) dum=sqrt5
            xin(nn)=dum* (g2(j+3*nj)+g(j+3*nj)+g(j+3*nj))
            yin(nn)=dum* (g2(j+9*nj)+g(j))
            zin(nn)=dum*  g2(j+12*nj)
            nn=nn+1
            xin(nn)=dum* (g2(j+4*nj)+g(j+4*nj)+g(j+4*nj))
            yin(nn)=dum*  g2(j+12*nj)
            zin(nn)=dum* (g2(j+10*nj)+g(j))
            nn=nn+1
            xin(nn)=dum* (g2(j+9*nj)+g(j+nj))
            yin(nn)=dum* (g2(j+5*nj)+g(j+3*nj)+g(j+3*nj))
            zin(nn)=dum*  g2(j+13*nj)
            nn=nn+1
            xin(nn)=dum*  g2(j+13*nj)
            yin(nn)=dum* (g2(j+6*nj)+g(j+5*nj)+g(j+5*nj))
            zin(nn)=dum* (g2(j+11*nj)+g(j+nj))
            nn=nn+1
            xin(nn)=dum* (g2(j+10*nj)+g(j+2*nj))
            yin(nn)=dum*  g2(j+14*nj)
            zin(nn)=dum* (g2(j+13*nj)+g(j+4*nj)+g(j+4*nj))
            nn=nn+1
            xin(nn)=dum*  g2(j+14*nj)
            yin(nn)=dum* (g2(j+11*nj)+g(j+2*nj))
            zin(nn)=dum* (g2(j+8*nj)+g(j+5*nj)+g(j+5*nj))
            nn=nn+1
            dum=1.0d+00
            if(norm) dum=sqrt3*sqrt5
            xin(nn)=dum* (g2(j+12*nj)+g(j+5*nj))
            yin(nn)=dum* (g2(j+13*nj)+g(j+4*nj))
            zin(nn)=dum* (g2(j+14*nj)+g(j+3*nj))
 400      continue
          return
        end if
c   g'
        if (mini.le.21) then
          do 500 j=1,nj
            nn = nn+1
            xin(nn) = (g2(j)+4.0d+00*g(j))
            yin(nn) =  g2(j+3*nj)
            zin(nn) =  g2(j+4*nj)
            nn = nn+1
            xin(nn) =  g2(j+5*nj)
            yin(nn) = (g2(j+nj)+4.0d+00*g(j+nj))
            zin(nn) =  g2(j+6*nj)
            nn = nn+1
            xin(nn) =  g2(j+7*nj)
            yin(nn) =  g2(j+8*nj)
            zin(nn) = (g2(j+2*nj)+4.0d+00*g(j+2*nj))
            dum = 1.0d+00
            if (norm) dum = sqrt7
            nn = nn+1
            xin(nn) =dum* (g2(j+3*nj)+3.0d+00*g(j+3*nj))
            yin(nn) =dum*  g2(j+9*nj)+g(j)
            zin(nn) =dum*  g2(j+15*nj)
            nn = nn+1
            xin(nn) =dum* (g2(j+4*nj)+3.0d+00*g(j+4*nj))
            yin(nn) =dum*  g2(j+15*nj)
            zin(nn) =dum*  g2(j+10*nj)+g(j)
            nn = nn+1
            xin(nn) =dum* (g2(j+11*nj)+g(j+nj))
            yin(nn) =dum*  g2(j+5*nj)+3.0d+00*g(j+5*nj)
            zin(nn) =dum*  g2(j+16*nj)
            nn = nn+1
            xin(nn) =dum*  g2(j+16*nj)
            yin(nn) =dum*  g2(j+6*nj)+3.0d+00*g(j+6*nj)
            zin(nn) =dum*  g2(j+12*nj)+g(j+nj)
            nn = nn+1
            xin(nn) =dum*  g2(j+13*nj)+g(j+2*nj)
            yin(nn) =dum*  g2(j+17*nj)
            zin(nn) =dum*  g2(j+7*nj)+3.0d+00*g(j+7*nj)
            nn = nn+1
            xin(nn) =dum*  g2(j+17*nj)
            yin(nn) =dum*  g2(j+14*nj)+g(j+2*nj)
            zin(nn) =dum*  g2(j+8*nj)+3.0d+00*g(j+8*nj)
            if (norm) dum = dum*sqrt5/sqrt3
            nn = nn+1
            xin(nn) =dum*  g2(j+9*nj)+2.0d+00*g(j+5*nj)
            yin(nn) =dum*  g2(j+11*nj)+2.0d+00*g(j+3*nj)
            zin(nn) =dum*  g2(j+18*nj)
            nn = nn+1
            xin(nn) =dum*  g2(j+10*nj)+2.0d+00*g(j+7*nj)
            yin(nn) =dum*  g2(j+19*nj)
            zin(nn) =dum*  g2(j+13*nj)+2.0d+00*g(j+4*nj)
            nn = nn+1
            xin(nn) =dum*  g2(j+20*nj)
            yin(nn) =dum*  g2(j+12*nj)+2.0d+00*g(j+8*nj)
            zin(nn) =dum*  g2(j+14*nj)+2.0d+00*g(j+6*nj)
            if (norm) dum = dum*sqrt3
            nn = nn+1
            xin(nn) =dum*  g2(j+15*nj)+2.0d+00*g(j+9*nj)
            yin(nn) =dum*  g2(j+18*nj)+g(j+4*nj)
            zin(nn) =dum*  g2(j+19*nj)+g(j+3*nj)
            nn = nn+1
            xin(nn) =dum*  g2(j+18*nj)+g(j+6*nj)
            yin(nn) =dum*  g2(j+16*nj)+2.0d+00*g(j+9*nj)
            zin(nn) =dum*  g2(j+20*nj)+g(j+5*nj)
            nn = nn+1
            xin(nn) =dum*  g2(j+19*nj)+g(j+8*nj)
            yin(nn) =dum*  g2(j+20*nj)+g(j+7*nj)
            zin(nn) =dum*  g2(j+17*nj)+2.0d+00*g(j+9*nj)
 500      continue
          return
        end if
      write(6,*)'in formdr h attempt',mini,maxi,minj,maxj
      call caserr('dimensioning problem in formdr - contact authors')
      return
      end
**==formii.f
      subroutine formii(gp1p1,gp1m1,gm1p1,gm1m1,gpp,mini,maxi,
     *                  minj,maxj,norm)
      implicit real*8 (a-h,o-z)
      logical norm
      dimension gp1p1(*), gp1m1(*), gm1p1(*), gm1m1(*), gpp(9,*)
      parameter (sqrt3=1.73205080756888d+00,sqrt5=2.23606797749979d+00)
      parameter (sqrt7=2.64575131106459d+00, zero=1.0d+00)
c
c     local storage
      dimension xm1(225),ym1(225),zm1(225),xp1(225),yp1(225),zp1(225)
c
c   form the second derivative integrals of type <d'| |d'>
c   this routine is completly general and good through g',g'
c
      do 50 i=1,225
        xm1(i)=zero
        ym1(i)=zero
        zm1(i)=zero
        xp1(i)=zero
        yp1(i)=zero
        zp1(i)=zero
  50  continue
c
c   first form the integrals due to the derivative on the left
c   since this is just like a first derivative formdr is used
c
      minjm1 = 1
      maxjm1 = 1
      minjp1 = 2
      maxjp1 = 4
      jtype = 1
      if (minj.eq.2) then
        minjp1 = 5
        maxjp1 = 10
        jtype = 2
      else if (minj.eq.5) then
        minjm1 = 2
        maxjm1 = 4
        minjp1 = 11
        maxjp1 = 20
        jtype = 3
      else if (minj.eq.11) then
        minjm1 = 5
        maxjm1 = 10
        minjp1 = 21
        maxjp1 = 35
        jtype = 4
      else if (minj.eq.21) then
        minjm1 = 11
        maxjm1 = 20
        minjp1 = 36
        maxjp1 = 56
        jtype = 5
      end if
c
      if (jtype.ge.2)
     *  call formdr(gm1m1,gp1m1,xm1,ym1,zm1,mini,maxi,minjm1,maxjm1,
     *              norm)
      call formdr(gm1p1,gp1p1,xp1,yp1,zp1,mini,maxi,minjp1,maxjp1,norm)
c
c  now form the total second derivative by forming the derivative on the
c  right side
c
c  note that this code is similar to the regular formdr code, but we are
c  forming the derivative on the right and forming a second derivative
c
      nn=0
      ni = maxi-mini+1
c   s'
      if (minj.eq.1) then
        do 100 i=1,ni
          gpp(1,i) = xp1(i)
          gpp(2,i) = xp1(i+ni)
          gpp(3,i) = xp1(i+2*ni)
          gpp(4,i) = yp1(i)
          gpp(5,i) = yp1(i+ni)
          gpp(6,i) = yp1(i+2*ni)
          gpp(7,i) = zp1(i)
          gpp(8,i) = zp1(i+ni)
          gpp(9,i) = zp1(i+2*ni)
 100    continue
        return
      end if
c   p'
      if (minj.le.2) then
        do 200 i=1,ni
          nn = i
          gpp(1,nn) = xp1(i)+xm1(i)
          gpp(2,nn) = xp1(i+3*ni)
          gpp(3,nn) = xp1(i+4*ni)
          gpp(4,nn) = yp1(i)+ym1(i)
          gpp(5,nn) = yp1(i+3*ni)
          gpp(6,nn) = yp1(i+4*ni)
          gpp(7,nn) = zp1(i)+zm1(i)
          gpp(8,nn) = zp1(i+3*ni)
          gpp(9,nn) = zp1(i+4*ni)
          nn = i+ni
          gpp(1,nn) = xp1(i+3*ni)
          gpp(2,nn) = xp1(i+ni)+xm1(i)
          gpp(3,nn) = xp1(i+5*ni)
          gpp(4,nn) = yp1(i+3*ni)
          gpp(5,nn) = yp1(i+ni)+ym1(i)
          gpp(6,nn) = yp1(i+5*ni)
          gpp(7,nn) = zp1(i+3*ni)
          gpp(8,nn) = zp1(i+ni)+zm1(i)
          gpp(9,nn) = zp1(i+5*ni)
          nn = i+2*ni
          gpp(1,nn) = xp1(i+4*ni)
          gpp(2,nn) = xp1(i+5*ni)
          gpp(3,nn) = xp1(i+2*ni)+xm1(i)
          gpp(4,nn) = yp1(i+4*ni)
          gpp(5,nn) = yp1(i+5*ni)
          gpp(6,nn) = yp1(i+2*ni)+ym1(i)
          gpp(7,nn) = zp1(i+4*ni)
          gpp(8,nn) = zp1(i+5*ni)
          gpp(9,nn) = zp1(i+2*ni)+zm1(i)
 200    continue
        return
      end if
c   d'
      if (minj.le.5) then
        do 300 i=1,ni
          nn = i
          gpp(1,nn) = xp1(i)+2.0d+00*xm1(i)
          gpp(2,nn) = xp1(i+3*ni)
          gpp(3,nn) = xp1(i+4*ni)
          gpp(4,nn) = yp1(i)+2.0d+00*ym1(i)
          gpp(5,nn) = yp1(i+3*ni)
          gpp(6,nn) = yp1(i+4*ni)
          gpp(7,nn) = zp1(i)+2.0d+00*zm1(i)
          gpp(8,nn) = zp1(i+3*ni)
          gpp(9,nn) = zp1(i+4*ni)
          nn = i + ni
          gpp(1,nn) = xp1(i+5*ni)
          gpp(2,nn) = xp1(i+ni)+2.0d+00*xm1(i+ni)
          gpp(3,nn) = xp1(i+6*ni)
          gpp(4,nn) = yp1(i+5*ni)
          gpp(5,nn) = yp1(i+ni)+2.0d+00*ym1(i+ni)
          gpp(6,nn) = yp1(i+6*ni)
          gpp(7,nn) = zp1(i+5*ni)
          gpp(8,nn) = zp1(i+ni)+2.0d+00*zm1(i+ni)
          gpp(9,nn) = zp1(i+6*ni)
          nn = i + 2*ni
          gpp(1,nn) = xp1(i+7*ni)
          gpp(2,nn) = xp1(i+8*ni)
          gpp(3,nn) = xp1(i+2*ni)+2.0d+00*xm1(i+2*ni)
          gpp(4,nn) = yp1(i+7*ni)
          gpp(5,nn) = yp1(i+8*ni)
          gpp(6,nn) = yp1(i+2*ni)+2.0d+00*ym1(i+2*ni)
          gpp(7,nn) = zp1(i+7*ni)
          gpp(8,nn) = zp1(i+8*ni)
          gpp(9,nn) = zp1(i+2*ni)+2.0d+00*zm1(i+2*ni)
          nn = i + 3*ni
          dum=1.0d+00
          if(norm) dum=sqrt3
          gpp(1,nn) = dum*(xp1(i+3*ni)+xm1(i+ni))
          gpp(2,nn) = dum*(xp1(i+5*ni)+xm1(i))
          gpp(3,nn) = dum* xp1(i+9*ni)
          gpp(4,nn) = dum*(yp1(i+3*ni)+ym1(i+ni))
          gpp(5,nn) = dum*(yp1(i+5*ni)+ym1(i))
          gpp(6,nn) = dum* yp1(i+9*ni)
          gpp(7,nn) = dum*(zp1(i+3*ni)+zm1(i+ni))
          gpp(8,nn) = dum*(zp1(i+5*ni)+zm1(i))
          gpp(9,nn) = dum* zp1(i+9*ni)
          nn = i + 4*ni
          gpp(1,nn) = dum*(xp1(i+4*ni)+xm1(i+2*ni))
          gpp(2,nn) = dum* xp1(i+9*ni)
          gpp(3,nn) = dum*(xp1(i+7*ni)+xm1(i))
          gpp(4,nn) = dum*(yp1(i+4*ni)+ym1(i+2*ni))
          gpp(5,nn) = dum* yp1(i+9*ni)
          gpp(6,nn) = dum*(yp1(i+7*ni)+ym1(i))
          gpp(7,nn) = dum*(zp1(i+4*ni)+zm1(i+2*ni))
          gpp(8,nn) = dum* zp1(i+9*ni)
          gpp(9,nn) = dum*(zp1(i+7*ni)+zm1(i))
          nn = i + 5*ni
          gpp(1,nn) = dum* xp1(i+9*ni)
          gpp(2,nn) = dum*(xp1(i+6*ni)+xm1(i+2*ni))
          gpp(3,nn) = dum*(xp1(i+8*ni)+xm1(i+ni))
          gpp(4,nn) = dum* yp1(i+9*ni)
          gpp(5,nn) = dum*(yp1(i+6*ni)+ym1(i+2*ni))
          gpp(6,nn) = dum*(yp1(i+8*ni)+ym1(i+ni))
          gpp(7,nn) = dum* zp1(i+9*ni)
          gpp(8,nn) = dum*(zp1(i+6*ni)+zm1(i+2*ni))
          gpp(9,nn) = dum*(zp1(i+8*ni)+zm1(i+ni))
 300    continue
        return
      end if
c   f'
      if (mini.le.11) then
        do 400 j=1,ni
          nn=i
          gpp(1,nn) = xp1(i)+3.0d+00*xm1(i)
          gpp(2,nn) = xp1(i+3*ni)
          gpp(3,nn) = xp1(i+4*ni)
          gpp(4,nn) = yp1(i)+3.0d+00*ym1(i)
          gpp(5,nn) = yp1(i+3*ni)
          gpp(6,nn) = yp1(i+4*ni)
          gpp(7,nn) = zp1(i)+3.0d+00*zm1(i)
          gpp(8,nn) = zp1(i+3*ni)
          gpp(9,nn) = zp1(i+4*ni)
          nn=i+ni
          gpp(1,nn) = xp1(i+5*ni)
          gpp(2,nn) = xp1(i+ni)+3.0d+00*xm1(i+ni)
          gpp(3,nn) = xp1(i+6*ni)
          gpp(4,nn) = yp1(i+5*ni)
          gpp(5,nn) = yp1(i+ni)+3.0d+00*ym1(i+ni)
          gpp(6,nn) = yp1(i+6*ni)
          gpp(7,nn) = zp1(i+5*ni)
          gpp(8,nn) = zp1(i+ni)+3.0d+00*zm1(i+ni)
          gpp(9,nn) = zp1(i+6*ni)
          nn=i+2*ni
          gpp(1,nn) = xp1(i+7*ni)
          gpp(2,nn) = xp1(i+8*ni)
          gpp(3,nn) = xp1(i+2*ni)+3.0d+00*xm1(i+2*ni)
          gpp(4,nn) = yp1(i+7*ni)
          gpp(5,nn) = yp1(i+8*ni)
          gpp(6,nn) = yp1(i+2*ni)+3.0d+00*ym1(i+2*ni)
          gpp(7,nn) = zp1(i+7*ni)
          gpp(8,nn) = zp1(i+8*ni)
          gpp(9,nn) = zp1(i+2*ni)+3.0d+00*zm1(i+2*ni)
          nn=i+3*ni
          dum = 1.0d+00
          if(norm) dum=sqrt5
          gpp(1,nn) = dum*(xp1(i+3*ni)+2.0d+00*xm1(i+3*ni))
          gpp(2,nn) = dum*(xp1(i+9*ni)+xm1(i))
          gpp(3,nn) = dum* xp1(i+12*ni)
          gpp(4,nn) = dum*(yp1(i+3*ni)+2.0d+00*ym1(i+3*ni))
          gpp(5,nn) = dum*(yp1(i+9*ni)+ym1(i))
          gpp(6,nn) = dum* yp1(i+12*ni)
          gpp(7,nn) = dum*(zp1(i+3*ni)+2.0d+00*zm1(i+3*ni))
          gpp(8,nn) = dum*(zp1(i+9*ni)+zm1(i))
          gpp(9,nn) = dum* zp1(i+12*ni)
          nn=i+4*ni
          gpp(1,nn) = dum*(xp1(i+4*ni)+2.0d+00*xm1(i+4*ni))
          gpp(2,nn) = dum* xp1(i+12*ni)
          gpp(3,nn) = dum*(xp1(i+10*ni)+xm1(i))
          gpp(4,nn) = dum*(yp1(i+4*ni)+2.0d+00*ym1(i+4*ni))
          gpp(5,nn) = dum* yp1(i+12*ni)
          gpp(6,nn) = dum*(yp1(i+10*ni)+ym1(i))
          gpp(7,nn) = dum*(zp1(i+4*ni)+2.0d+00*zm1(i+4*ni))
          gpp(8,nn) = dum* zp1(i+12*ni)
          gpp(9,nn) = dum*(zp1(i+10*ni)+zm1(i))
          nn=i+5*ni
          gpp(1,nn) = dum*(xp1(i+9*ni)+xm1(i+ni))
          gpp(2,nn) = dum*(xp1(i+5*ni)+2.0d+00*xm1(i+3*ni))
          gpp(3,nn) = dum* xp1(i+13*ni)
          gpp(4,nn) = dum*(yp1(i+9*ni)+ym1(i+ni))
          gpp(5,nn) = dum*(yp1(i+5*ni)+2.0d+00*ym1(i+3*ni))
          gpp(6,nn) = dum* yp1(i+13*ni)
          gpp(7,nn) = dum*(zp1(i+9*ni)+zm1(i+ni))
          gpp(8,nn) = dum*(zp1(i+5*ni)+2.0d+00*zm1(i+3*ni))
          gpp(9,nn) = dum* zp1(i+13*ni)
          nn=i+6*ni
          gpp(1,nn) = dum* xp1(i+13*ni)
          gpp(2,nn) = dum*(xp1(i+6*ni)+2.0d+00*xm1(i+5*ni))
          gpp(3,nn) = dum*(xp1(i+11*ni)+xm1(i+ni))
          gpp(4,nn) = dum* yp1(i+13*ni)
          gpp(5,nn) = dum*(yp1(i+6*ni)+2.0d+00*ym1(i+5*ni))
          gpp(6,nn) = dum*(yp1(i+11*ni)+ym1(i+ni))
          gpp(7,nn) = dum* zp1(i+13*ni)
          gpp(8,nn) = dum*(zp1(i+6*ni)+2.0d+00*zm1(i+5*ni))
          gpp(9,nn) = dum*(zp1(i+11*ni)+zm1(i+ni))
          nn=i+7*ni
          gpp(1,nn) = dum*(xp1(i+10*ni)+xm1(i+2*ni))
          gpp(2,nn) = dum* xp1(i+14*ni)
          gpp(3,nn) = dum*(xp1(i+13*ni)+2.0d+00*xm1(i+4*ni))
          gpp(4,nn) = dum*(yp1(i+10*ni)+ym1(i+2*ni))
          gpp(5,nn) = dum* yp1(i+14*ni)
          gpp(6,nn) = dum*(yp1(i+13*ni)+2.0d+00*ym1(i+4*ni))
          gpp(7,nn) = dum*(zp1(i+10*ni)+zm1(i+2*ni))
          gpp(8,nn) = dum* zp1(i+14*ni)
          gpp(9,nn) = dum*(zp1(i+13*ni)+2.0d+00*zm1(i+4*ni))
          nn=i+8*ni
          gpp(1,nn) = dum* xp1(i+14*ni)
          gpp(2,nn) = dum*(xp1(i+11*ni)+xm1(i+2*ni))
          gpp(3,nn) = dum*(xp1(i+8*ni)+2.0d+00*xm1(i+5*ni))
          gpp(4,nn) = dum* yp1(i+14*ni)
          gpp(5,nn) = dum*(yp1(i+11*ni)+ym1(i+2*ni))
          gpp(6,nn) = dum*(yp1(i+8*ni)+2.0d+00*ym1(i+5*ni))
          gpp(7,nn) = dum* zp1(i+14*ni)
          gpp(8,nn) = dum*(zp1(i+11*ni)+zm1(i+2*ni))
          gpp(9,nn) = dum*(zp1(i+8*ni)+2.0d+00*zm1(i+5*ni))
          nn=i+9*ni
          dum=1.0d+00
          if(norm) dum=sqrt3*sqrt5
          gpp(1,nn) = dum*(xp1(i+12*ni)+xm1(i+5*ni))
          gpp(2,nn) = dum*(xp1(i+13*ni)+xm1(i+4*ni))
          gpp(3,nn) = dum*(xp1(i+14*ni)+xm1(i+3*ni))
          gpp(4,nn) = dum*(yp1(i+12*ni)+ym1(i+5*ni))
          gpp(5,nn) = dum*(yp1(i+13*ni)+ym1(i+4*ni))
          gpp(6,nn) = dum*(yp1(i+14*ni)+ym1(i+3*ni))
          gpp(7,nn) = dum*(zp1(i+12*ni)+zm1(i+5*ni))
          gpp(8,nn) = dum*(zp1(i+13*ni)+zm1(i+4*ni))
          gpp(9,nn) = dum*(zp1(i+14*ni)+zm1(i+3*ni))
 400    continue
        return
      end if
c   g'
      if (mini.le.21) then
        do 500 j=1,ni
          nn = i
          gpp(1,nn) = xp1(i)+4.0d+00*xm1(i)
          gpp(2,nn) = xp1(i+3*ni)
          gpp(3,nn) = xp1(i+4*ni)
          gpp(4,nn) = yp1(i)+4.0d+00*ym1(i)
          gpp(5,nn) = yp1(i+3*ni)
          gpp(6,nn) = yp1(i+4*ni)
          gpp(7,nn) = zp1(i)+4.0d+00*zm1(i)
          gpp(8,nn) = zp1(i+3*ni)
          gpp(9,nn) = zp1(i+4*ni)
          nn = i+ni
          gpp(1,nn) = xp1(i+5*ni)
          gpp(2,nn) = xp1(i+ni)+4.0d+00*xm1(i+ni)
          gpp(3,nn) = xp1(i+6*ni)
          gpp(4,nn) = yp1(i+5*ni)
          gpp(5,nn) = yp1(i+ni)+4.0d+00*ym1(i+ni)
          gpp(6,nn) = yp1(i+6*ni)
          gpp(7,nn) = zp1(i+5*ni)
          gpp(8,nn) = zp1(i+ni)+4.0d+00*zm1(i+ni)
          gpp(9,nn) = zp1(i+6*ni)
          nn = i+2*ni
          gpp(1,nn) = xp1(i+7*ni)
          gpp(2,nn) = xp1(i+8*ni)
          gpp(3,nn) = xp1(i+2*ni)+4.0d+00*xm1(i+2*ni)
          gpp(4,nn) = yp1(i+7*ni)
          gpp(5,nn) = yp1(i+8*ni)
          gpp(6,nn) = yp1(i+2*ni)+4.0d+00*ym1(i+2*ni)
          gpp(7,nn) = zp1(i+7*ni)
          gpp(8,nn) = zp1(i+8*ni)
          gpp(9,nn) = zp1(i+2*ni)+4.0d+00*zm1(i+2*ni)
          dum = 1.0d+00
          if (norm) dum = sqrt7
          nn = i+3*ni
          gpp(1,nn) = dum*(xp1(i+3*ni)+3.0d+00*xm1(i+3*ni))
          gpp(2,nn) = dum*(xp1(i+9*ni)+xm1(i))
          gpp(3,nn) = dum* xp1(i+15*ni)
          gpp(4,nn) = dum*(yp1(i+3*ni)+3.0d+00*ym1(i+3*ni))
          gpp(5,nn) = dum*(yp1(i+9*ni)+ym1(i))
          gpp(6,nn) = dum* yp1(i+15*ni)
          gpp(7,nn) = dum*(zp1(i+3*ni)+3.0d+00*zm1(i+3*ni))
          gpp(8,nn) = dum*(zp1(i+9*ni)+zm1(i))
          gpp(9,nn) = dum* zp1(i+15*ni)
          nn = i+4*ni
          gpp(1,nn) = dum*(xp1(i+4*ni)+3.0d+00*xm1(i+4*ni))
          gpp(2,nn) = dum* xp1(i+15*ni)
          gpp(3,nn) = dum*(xp1(i+10*ni)+xm1(i))
          gpp(4,nn) = dum*(yp1(i+4*ni)+3.0d+00*ym1(i+4*ni))
          gpp(5,nn) = dum* yp1(i+15*ni)
          gpp(6,nn) = dum*(yp1(i+10*ni)+ym1(i))
          gpp(7,nn) = dum*(zp1(i+4*ni)+3.0d+00*zm1(i+4*ni))
          gpp(8,nn) = dum* zp1(i+15*ni)
          gpp(9,nn) = dum*(zp1(i+10*ni)+zm1(i))
          nn = i+5*ni
          gpp(1,nn) = dum*(xp1(i+11*ni)+xm1(i+ni))
          gpp(2,nn) = dum*(xp1(i+5*ni)+3.0d+00*xm1(i+5*ni))
          gpp(3,nn) = dum* xp1(i+16*ni)
          gpp(4,nn) = dum*(yp1(i+11*ni)+ym1(i+ni))
          gpp(5,nn) = dum*(yp1(i+5*ni)+3.0d+00*ym1(i+5*ni))
          gpp(6,nn) = dum* yp1(i+16*ni)
          gpp(7,nn) = dum*(zp1(i+11*ni)+zm1(i+ni))
          gpp(8,nn) = dum*(zp1(i+5*ni)+3.0d+00*zm1(i+5*ni))
          gpp(9,nn) = dum* zp1(i+16*ni)
          nn = i+6*ni
          gpp(1,nn) = dum* xp1(i+16*ni)
          gpp(2,nn) = dum*(xp1(i+6*ni)+3.0d+00*xm1(i+6*ni))
          gpp(3,nn) = dum*(xp1(i+12*ni)+xm1(i+ni))
          gpp(4,nn) = dum* yp1(i+16*ni)
          gpp(5,nn) = dum*(yp1(i+6*ni)+3.0d+00*ym1(i+6*ni))
          gpp(6,nn) = dum*(yp1(i+12*ni)+ym1(i+ni))
          gpp(7,nn) = dum* zp1(i+16*ni)
          gpp(8,nn) = dum*(zp1(i+6*ni)+3.0d+00*zm1(i+6*ni))
          gpp(9,nn) = dum*(zp1(i+12*ni)+zm1(i+ni))
          nn = i+7*ni
          gpp(1,nn) = dum*(xp1(i+13*ni)+xm1(i+2*ni))
          gpp(2,nn) = dum* xp1(i+17*ni)
          gpp(3,nn) = dum*(xp1(i+7*ni)+3.0d+00*xm1(i+7*ni))
          gpp(4,nn) = dum*(yp1(i+13*ni)+ym1(i+2*ni))
          gpp(5,nn) = dum* yp1(i+17*ni)
          gpp(6,nn) = dum*(yp1(i+7*ni)+3.0d+00*ym1(i+7*ni))
          gpp(7,nn) = dum*(zp1(i+13*ni)+zm1(i+2*ni))
          gpp(8,nn) = dum* zp1(i+17*ni)
          gpp(9,nn) = dum*(zp1(i+7*ni)+3.0d+00*zm1(i+7*ni))
          nn = i+8*ni
          gpp(1,nn) = dum* xp1(i+17*ni)
          gpp(2,nn) = dum*(xp1(i+14*ni)+xm1(i+2*ni))
          gpp(3,nn) = dum*(xp1(i+8*ni)+3.0d+00*xm1(i+8*ni))
          gpp(4,nn) = dum* yp1(i+17*ni)
          gpp(5,nn) = dum*(yp1(i+14*ni)+ym1(i+2*ni))
          gpp(6,nn) = dum*(yp1(i+8*ni)+3.0d+00*ym1(i+8*ni))
          gpp(7,nn) = dum* zp1(i+17*ni)
          gpp(8,nn) = dum*(zp1(i+14*ni)+zm1(i+2*ni))
          gpp(9,nn) = dum*(zp1(i+8*ni)+3.0d+00*zm1(i+8*ni))
          if (norm) dum = dum*sqrt5/sqrt3
          nn = i+9*ni
          gpp(1,nn) = dum*(xp1(i+9*ni)+2.0d+00*xm1(i+5*ni))
          gpp(2,nn) = dum*(xp1(i+11*ni)+2.0d+00*xm1(i+3*ni))
          gpp(3,nn) = dum* xp1(i+18*ni)
          gpp(4,nn) = dum*(yp1(i+9*ni)+2.0d+00*ym1(i+5*ni))
          gpp(5,nn) = dum*(yp1(i+11*ni)+2.0d+00*ym1(i+3*ni))
          gpp(6,nn) = dum* yp1(i+18*ni)
          gpp(7,nn) = dum*(zp1(i+9*ni)+2.0d+00*zm1(i+5*ni))
          gpp(8,nn) = dum*(zp1(i+11*ni)+2.0d+00*zm1(i+3*ni))
          gpp(9,nn) = dum* zp1(i+18*ni)
          nn = i+10*ni
          gpp(1,nn) = dum*(xp1(i+10*ni)+2.0d+00*xm1(i+7*ni))
          gpp(2,nn) = dum* xp1(i+19*ni)
          gpp(3,nn) = dum*(xp1(i+13*ni)+2.0d+00*xm1(i+4*ni))
          gpp(4,nn) = dum*(yp1(i+10*ni)+2.0d+00*ym1(i+7*ni))
          gpp(5,nn) = dum* yp1(i+19*ni)
          gpp(6,nn) = dum*(yp1(i+13*ni)+2.0d+00*ym1(i+4*ni))
          gpp(7,nn) = dum*(zp1(i+10*ni)+2.0d+00*zm1(i+7*ni))
          gpp(8,nn) = dum* zp1(i+19*ni)
          gpp(9,nn) = dum*(zp1(i+13*ni)+2.0d+00*zm1(i+4*ni))
          nn = i+11*ni
          gpp(1,nn) = dum* xp1(i+20*ni)
          gpp(2,nn) = dum*(xp1(i+12*ni)+2.0d+00*xm1(i+8*ni))
          gpp(3,nn) = dum*(xp1(i+14*ni)+2.0d+00*xm1(i+6*ni))
          gpp(4,nn) = dum* yp1(i+20*ni)
          gpp(5,nn) = dum*(yp1(i+12*ni)+2.0d+00*ym1(i+8*ni))
          gpp(6,nn) = dum*(yp1(i+14*ni)+2.0d+00*ym1(i+6*ni))
          gpp(7,nn) = dum* zp1(i+20*ni)
          gpp(8,nn) = dum*(zp1(i+12*ni)+2.0d+00*zm1(i+8*ni))
          gpp(9,nn) = dum*(zp1(i+14*ni)+2.0d+00*zm1(i+6*ni))
          if (norm) dum = dum*sqrt3
          nn = i+12*ni
          gpp(1,nn) = dum*(xp1(i+15*ni)+2.0d+00*xm1(i+9*ni))
          gpp(2,nn) = dum*(xp1(i+18*ni)+xm1(i+4*ni))
          gpp(3,nn) = dum*(xp1(i+19*ni)+xm1(i+3*ni))
          gpp(4,nn) = dum*(yp1(i+15*ni)+2.0d+00*ym1(i+9*ni))
          gpp(5,nn) = dum*(yp1(i+18*ni)+ym1(i+4*ni))
          gpp(6,nn) = dum*(yp1(i+19*ni)+ym1(i+3*ni))
          gpp(7,nn) = dum*(zp1(i+15*ni)+2.0d+00*zm1(i+9*ni))
          gpp(8,nn) = dum*(zp1(i+18*ni)+zm1(i+4*ni))
          gpp(9,nn) = dum*(zp1(i+19*ni)+zm1(i+3*ni))
          nn = i+13*ni
          gpp(1,nn) = dum*(xp1(i+18*ni)+xm1(i+6*ni))
          gpp(2,nn) = dum*(xp1(i+16*ni)+2.0d+00*xm1(i+9*ni))
          gpp(3,nn) = dum*(xp1(i+20*ni)+xm1(i+5*ni))
          gpp(4,nn) = dum*(yp1(i+18*ni)+ym1(i+6*ni))
          gpp(5,nn) = dum*(yp1(i+16*ni)+2.0d+00*ym1(i+9*ni))
          gpp(6,nn) = dum*(yp1(i+20*ni)+ym1(i+5*ni))
          gpp(7,nn) = dum*(zp1(i+18*ni)+zm1(i+6*ni))
          gpp(8,nn) = dum*(zp1(i+16*ni)+2.0d+00*zm1(i+9*ni))
          gpp(9,nn) = dum*(zp1(i+20*ni)+zm1(i+5*ni))
          nn = i+14*ni
          gpp(1,nn) = dum*(xp1(i+19*ni)+xm1(i+8*ni))
          gpp(2,nn) = dum*(xp1(i+20*ni)+xm1(i+7*ni))
          gpp(3,nn) = dum*(xp1(i+17*ni)+2.0d+00*xm1(i+9*ni))
          gpp(4,nn) = dum*(yp1(i+19*ni)+ym1(i+8*ni))
          gpp(5,nn) = dum*(yp1(i+20*ni)+ym1(i+7*ni))
          gpp(6,nn) = dum*(yp1(i+17*ni)+2.0d+00*ym1(i+9*ni))
          gpp(7,nn) = dum*(zp1(i+19*ni)+zm1(i+8*ni))
          gpp(8,nn) = dum*(zp1(i+20*ni)+zm1(i+7*ni))
          gpp(9,nn) = dum*(zp1(i+17*ni)+2.0d+00*zm1(i+9*ni))
 500    continue
        return
      end if
      write(6,*)'in formii h attempt',mini,maxi,minj,maxj
      call caserr('dimensioning problem in formii - contact authors')
      return
      end
**==formij.f
      subroutine formij(gp2,g,gm2,gpp,mini,maxi,minj,maxj,norm)
      implicit real*8 (a-h,o-z)
      logical norm
      dimension g(*), gp2(*), gm2(*), gpp(9,*)
      parameter (sqrt3=1.73205080756888d+00,sqrt5=2.23606797749979d+00)
      parameter (sqrt7=2.64575131106459d+00)
c
c   form the second derivative integrals of type <d''| |d>
c   this routine is completly general and good through g''
c   remember that order of derivatives is not important so we compute
c   6 terms and equate to the other 3 (ie xy = yx)
c
      nn=0
      nj = maxj-minj+1
c   s''
        if (mini.eq.1) then
          do 100 j=1,nj
            nn=nn+1
            gpp(1,nn) = gp2(j)  + g(j)
            gpp(2,nn) = gp2(j+3*nj)
            gpp(3,nn) = gp2(j+4*nj)
            gpp(4,nn) = gpp(2,nn)
            gpp(5,nn) = gp2(j+nj) + g(j)
            gpp(6,nn) = gp2(j+5*nj)
            gpp(7,nn) = gpp(3,nn)
            gpp(8,nn) = gpp(6,nn)
            gpp(9,nn) = gp2(j+2*nj) + g(j) 
 100      continue
          return
        end if
c   p''
        if (mini.le.2) then
          do 200 j=1,nj
            nn = nn+1
            gpp(1,nn) = gp2(j) + 3.0d+00*g(j)
            gpp(2,nn) = gp2(j+3*nj) + g(j+nj)
            gpp(3,nn) = gp2(j+4*nj) + g(j+2*nj)
            gpp(4,nn) = gpp(2,nn)
            gpp(5,nn) = gp2(j+5*nj) + g(j)
            gpp(6,nn) = gp2(j+9*nj)
            gpp(7,nn) = gpp(3,nn)
            gpp(8,nn) = gpp(6,nn)
            gpp(9,nn) = gp2(j+7*nj) + g(j) 
            nn = nn+1
            gpp(1,nn) = gp2(j+3*nj) + g(j+nj)
            gpp(2,nn) = gp2(j+5*nj) + g(j)
            gpp(3,nn) = gp2(j+9*nj)
            gpp(4,nn) = gpp(2,nn)
            gpp(5,nn) = gp2(j+nj) + 3.0d+00*g(j+nj)
            gpp(6,nn) = gp2(j+6*nj) + g(j+2*nj)
            gpp(7,nn) = gpp(3,nn)
            gpp(8,nn) = gpp(6,nn)
            gpp(9,nn) = gp2(j+8*nj) + g(j+nj) 
            nn = nn+1
            gpp(1,nn) = gp2(j+4*nj) + g(j+2*nj)
            gpp(2,nn) = gp2(j+9*nj)
            gpp(3,nn) = gp2(j+7*nj) + g(j)
            gpp(4,nn) = gpp(2,nn)
            gpp(5,nn) = gp2(j+6*nj) + g(j+2*nj)
            gpp(6,nn) = gp2(j+8*nj) + g(j+nj)
            gpp(7,nn) = gpp(3,nn)
            gpp(8,nn) = gpp(6,nn)
            gpp(9,nn) = gp2(j+2*nj) + 3.0d+00*g(j+2*nj) 
 200      continue
          return
        end if
c   d''
        if (mini.le.5) then
          do 300 j=1,nj
            nn = nn+1
            gpp(1,nn) = gp2(j)  + 5.0d+00*g(j) + 2.0d+00*gm2(j)
            gpp(2,nn) = gp2(j+3*nj) + 2.0d+00*g(j+3*nj)
            gpp(3,nn) = gp2(j+4*nj) + 2.0d+00*g(j+4*nj)
            gpp(4,nn) = gpp(2,nn)
            gpp(5,nn) = gp2(j+9*nj) + g(j)
            gpp(6,nn) = gp2(j+12*nj)
            gpp(7,nn) = gpp(3,nn)
            gpp(8,nn) = gpp(6,nn)
            gpp(9,nn) = gp2(j+10*nj) + g(j) 
            nn = nn+1
            gpp(1,nn) = gp2(j+9*nj) + g(j+nj)
            gpp(2,nn) = gp2(j+5*nj) + 2.0d+00*g(j+3*nj)
            gpp(3,nn) = gp2(j+13*nj)
            gpp(4,nn) = gpp(2,nn)
            gpp(5,nn) = gp2(j+nj) + 5.0d+00*g(j+nj) + 2.0d+00*gm2(j)
            gpp(6,nn) = gp2(j+6*nj) + 2.0d+00*g(j+5*nj)
            gpp(7,nn) = gpp(3,nn)
            gpp(8,nn) = gpp(6,nn)
            gpp(9,nn) = gp2(j+11*nj) + g(j+nj)
            nn = nn+1
            gpp(1,nn) = gp2(j+10*nj) + g(j+2*nj)
            gpp(2,nn) = gp2(j+14*nj)
            gpp(3,nn) = gp2(j+7*nj) + 2.0d+00*g(j+4*nj)
            gpp(4,nn) = gpp(2,nn)
            gpp(5,nn) = gp2(j+11*nj) + g(j+2*nj)
            gpp(6,nn) = gp2(j+8*nj) + 2.0d+00*g(j+5*nj)
            gpp(7,nn) = gpp(3,nn)
            gpp(8,nn) = gpp(6,nn)
            gpp(9,nn) = gp2(j+2*nj)+5.0d+00*g(j+2*nj)+2.0d+00*gm2(j)
            nn = nn+1
            dum=1.0d+00
            if(norm) dum=sqrt3
            gpp(1,nn) = dum*(gp2(j+3*nj) + 3.0d+00*g(j+3*nj))
            gpp(2,nn) = dum*(gp2(j+9*nj)+g(j)+g(j+nj)+gm2(j))
            gpp(3,nn) = dum*(gp2(j+12*nj) + g(j+5*nj))
            gpp(4,nn) = gpp(2,nn)
            gpp(5,nn) = dum*(gp2(j+5*nj) + 3.0d+00*g(j+3*nj))
            gpp(6,nn) = dum*(gp2(j+13*nj) + g(j+4*nj))
            gpp(7,nn) = gpp(3,nn)
            gpp(8,nn) = gpp(6,nn)
            gpp(9,nn) = dum*(gp2(j+14*nj)+g(j+3*nj))
            nn = nn+1
            gpp(1,nn) = dum*(gp2(j+4*nj) + 3.0d+00*g(j+4*nj))
            gpp(2,nn) = dum*(gp2(j+12*nj)+g(j+5*nj))
            gpp(3,nn) = dum*(gp2(j+10*nj)+g(j)+g(j+2*nj)+gm2(j))
            gpp(4,nn) = gpp(2,nn)
            gpp(5,nn) = dum*(gp2(j+13*nj) + g(j+4*nj))
            gpp(6,nn) = dum*(gp2(j+14*nj) + g(j+3*nj))
            gpp(7,nn) = gpp(3,nn)
            gpp(8,nn) = gpp(6,nn)
            gpp(9,nn) = dum*(gp2(j+7*nj)+3.0d+00*g(j+4*nj))
            nn = nn+1
            gpp(1,nn) = dum*(gp2(j+12*nj) + g(j+5*nj))
            gpp(2,nn) = dum*(gp2(j+13*nj)+g(j+4*nj))
            gpp(3,nn) = dum*(gp2(j+14*nj)+g(j+3*nj))
            gpp(4,nn) = gpp(2,nn)
            gpp(5,nn) = dum*(gp2(j+6*nj) + 3.0d+00*g(j+5*nj))
            gpp(6,nn) = dum*(gp2(j+11*nj)+g(j+nj)+g(j+2*nj)+gm2(j))
            gpp(7,nn) = gpp(3,nn)
            gpp(8,nn) = gpp(6,nn)
            gpp(9,nn) = dum*(gp2(j+8*nj)+3.0d+00*g(j+5*nj))
 300      continue
          return
        end if
c   f''
        if (mini.le.11) then
          do 400 j=1,nj
            nn=nn+1
            gpp(1,nn) = gp2(j) + 7.0d+00*g(j) + 6.0d+00*gm2(j)
            gpp(2,nn) = gp2(j+3*nj)+3.0d+00*g(j+3*nj)
            gpp(3,nn) = gp2(j+4*nj)+3.0d+00*g(j+4*nj)
            gpp(4,nn) = gpp(2,nn)
            gpp(5,nn) = gp2(j+9*nj) + g(j)
            gpp(6,nn) = gp2(j+15*nj)
            gpp(7,nn) = gpp(3,nn)
            gpp(8,nn) = gpp(6,nn)
            gpp(9,nn) = gp2(j+10*nj)+g(j)
            nn=nn+1
            gpp(1,nn) = gp2(j+11*nj) + g(j+nj)
            gpp(2,nn) = gp2(j+5*nj)+3.0d+00*g(j+5*nj)
            gpp(3,nn) = gp2(j+16*nj)
            gpp(4,nn) = gpp(2,nn)
            gpp(5,nn) = gp2(j+nj)+7.0d+00*g(j+nj)+6.0d+00*gm2(j+nj)
            gpp(6,nn) = gp2(j+6*nj)+3.0d+00*g(j+6*nj)
            gpp(7,nn) = gpp(3,nn)
            gpp(8,nn) = gpp(6,nn)
            gpp(9,nn) = gp2(j+12*nj)+g(j+nj)
            nn=nn+1
            gpp(1,nn) = gp2(j+13*nj) + g(j+2*nj)
            gpp(2,nn) = gp2(j+17*nj)
            gpp(3,nn) = gp2(j+7*nj)+3.0d+00*g(j+7*nj)
            gpp(4,nn) = gpp(2,nn)
            gpp(5,nn) = gp2(j+14*nj)+g(j+2*nj)
            gpp(6,nn) = gp2(j+8*nj)+3.0d+00*g(j+8*nj)
            gpp(7,nn) = gpp(3,nn)
            gpp(8,nn) = gpp(6,nn)
            gpp(9,nn)=gp2(j+2*nj)+7.0d+00*g(j+2*nj)+6.0d+00*gm2(j+2*nj)
            nn=nn+1
            dum = 1.0d+00
            if(norm) dum=sqrt5
            gpp(1,nn) = dum*(gp2(j+3*nj)+5.0d+00*g(j+3*nj)+
     *                       2.0d+00*gm2(j+nj))
            gpp(2,nn) = dum*(gp2(j+9*nj)+g(j)+2.0d+00*g(j+5*nj)+
     *                       2.0d+00*gm2(j))
            gpp(3,nn) = dum*(gp2(j+15*nj)+2.0d+00*g(j+9*nj))
            gpp(4,nn) = gpp(2,nn)
            gpp(5,nn) = dum*(gp2(j+11*nj) + 3.0d+00*g(j+3*nj))
            gpp(6,nn) = dum*(gp2(j+18*nj)+g(j+4*nj))
            gpp(7,nn) = gpp(3,nn)
            gpp(8,nn) = gpp(6,nn)
            gpp(9,nn) = dum*(gp2(j+19*nj)+g(j+3*nj))
            nn=nn+1
            gpp(1,nn) = dum*(gp2(j+4*nj)+5.0d+00*g(j+4*nj)+
     *                       2.0d+00*gm2(j+2*nj))
            gpp(2,nn) = dum*(gp2(j+15*nj)+2.0d+00*g(j+9*nj))
            gpp(3,nn) = dum*(gp2(j+10*nj)+2.0d+00*g(j+7*nj)+g(j)+
     *                       2.0d+00*gm2(j))
            gpp(4,nn) = gpp(2,nn)
            gpp(5,nn) = dum*(gp2(j+18*nj) + g(j+7*nj))
            gpp(6,nn) = dum*(gp2(j+19*nj)+g(j+3*nj))
            gpp(7,nn) = gpp(3,nn)
            gpp(8,nn) = gpp(6,nn)
            gpp(9,nn) = dum*(gp2(j+13*nj)+3.0d+00*g(j+4*nj))
            nn=nn+1
            gpp(1,nn) = dum*(gp2(j+9*nj)+3.0d+00*g(j+5*nj))
            gpp(2,nn) = dum*(gp2(j+11*nj)+g(j+nj)+2.0d+00*g(j+3*nj)+
     *                       2.0d+00*gm2(j+nj))
            gpp(3,nn) = dum*(gp2(j+18*nj)+g(j+6*nj))
            gpp(4,nn) = gpp(2,nn)
            gpp(5,nn) = dum*(gp2(j+5*nj)+5.0d+00*g(j+5*nj)+
     *                       2.0d+00*gm2(j))
            gpp(6,nn) = dum*(gp2(j+16*nj)+2.0d+00*g(j+9*nj))
            gpp(7,nn) = gpp(3,nn)
            gpp(8,nn) = gpp(6,nn)
            gpp(9,nn) = dum*(gp2(j+20*nj)+g(j+3*nj))
            nn=nn+1
            gpp(1,nn) = dum*(gp2(j+18*nj) + g(j+6*nj))
            gpp(2,nn) = dum*(gp2(j+16*nj)+2.0d+00*g(j+9*nj))
            gpp(3,nn) = dum*(gp2(j+20*nj)+g(j+5*nj))
            gpp(4,nn) = gpp(2,nn)
            gpp(5,nn) = dum*(gp2(j+6*nj)+5.0d+00*g(j+6*nj)+
     *                       2.0d+00*gm2(j+2*nj))
            gpp(6,nn) = dum*(gp2(j+12*nj)+g(j+nj)+2.0d+00*g(j+8*nj)+
     *                       2.0d+00*gm2(j+nj))
            gpp(7,nn) = gpp(3,nn)
            gpp(8,nn) = gpp(6,nn)
            gpp(9,nn) = dum*(gp2(j+14*nj)+3.0d+00*g(j+6*nj))
            nn=nn+1
            gpp(1,nn) = dum*(gp2(j+10*nj)+3.0d+00*g(j+7*nj))
            gpp(2,nn) = dum*(gp2(j+19*nj)+g(j+8*nj))
            gpp(3,nn) = dum*(gp2(j+13*nj)+g(j+2*nj)+2.0d+00*g(j+4*nj)+
     *                       2.0d+00*gm2(j+2*nj))
            gpp(4,nn) = gpp(2,nn)
            gpp(5,nn) = dum*(gp2(j+20*nj)+g(j+7*nj))
            gpp(6,nn) = dum*(gp2(j+17*nj)+2.0d+00*g(j+9*nj))
            gpp(7,nn) = gpp(3,nn)
            gpp(8,nn) = gpp(6,nn)
            gpp(9,nn) = dum*(gp2(j+7*nj)+5.0d+00*g(j+7*nj)+
     *                       2.0d+00*gm2(j))
            nn=nn+1
            gpp(1,nn) = dum*(gp2(j+19*nj)+g(j+8*nj))
            gpp(2,nn) = dum*(gp2(j+20*nj)+g(j+7*nj))
            gpp(3,nn) = dum*(gp2(j+17*nj)+2.0d+00*g(j+9*nj))
            gpp(4,nn) = gpp(2,nn)
            gpp(5,nn) = dum*(gp2(j+12*nj)+3.0d+00*g(j+8*nj))
            gpp(6,nn) = dum*(gp2(j+14*nj)+g(j+2*nj)+2.0d+00*g(j+6*nj)+
     *                       2.0d+00*gm2(j+2*nj))
            gpp(7,nn) = gpp(3,nn)
            gpp(8,nn) = gpp(6,nn)
            gpp(9,nn) = dum*(gp2(j+8*nj)+5.0d+00*g(j+8*nj)+
     *                       2.0d+00*gm2(j+nj))
            nn=nn+1
            dum=1.0d+00
            if(norm) dum=sqrt3*sqrt5
            gpp(1,nn) = dum*(gp2(j+15*nj)+3.0d+00*g(j+9*nj))
            gpp(2,nn) = dum*(gp2(j+18*nj)+g(j+6*nj)+g(j+4*nj)+
     *                       gm2(j+2*nj))
            gpp(3,nn) = dum*(gp2(j+19*nj)+g(j+8*nj)+g(j+3*nj)+
     *                       gm2(j+nj))
            gpp(4,nn) = gpp(2,nn)
            gpp(5,nn) = dum*(gp2(j+16*nj)+3.0d+00*g(j+9*nj))
            gpp(6,nn) = dum*(gp2(j+20*nj)+g(j+7*nj)+g(j+5*nj)+gm2(j))
            gpp(7,nn) = gpp(3,nn)
            gpp(8,nn) = gpp(6,nn)
            gpp(9,nn) = dum*(gp2(j+17*nj)+3.0d+00*g(j+9*nj))
 400      continue
          return
        end if
c   g''
        if (mini.le.21) then
          do 500 j=1,nj
            nn = nn+1
            gpp(1,nn) = gp2(j)+9.0d+00*g(j)+12.0d+00*gm2(j)
            gpp(2,nn) = gp2(j+3*nj)+4.0d+00*g(j+3*nj)
            gpp(3,nn) = gp2(j+4*nj)+4.0d+00*g(j+4*nj)
            gpp(4,nn) = gpp(2,nn)
            gpp(5,nn) = gp2(j+9*nj)+g(j)
            gpp(6,nn) = gp2(j+15*nj)
            gpp(7,nn) = gpp(3,nn)
            gpp(8,nn) = gpp(6,nn)
            gpp(9,nn) = gp2(j+10*nj)+g(j)
            nn = nn+1
            gpp(1,nn) = gp2(j+11*nj)+g(j+nj)
            gpp(2,nn) = gp2(j+5*nj)+4.0d+00*g(j+5*nj)
            gpp(3,nn) = gp2(j+17*nj)
            gpp(4,nn) = gpp(2,nn)
            gpp(5,nn) = gp2(j+nj)+9.0d+00*g(j+nj)+12.0d+00*gm2(j+nj)
            gpp(6,nn) = gp2(j+6*nj)+4.0d+00*g(j+6*nj)
            gpp(7,nn) = gpp(3,nn)
            gpp(8,nn) = gpp(6,nn)
            gpp(9,nn) = gp2(j+12*nj)+g(j+nj)
            nn = nn+1
            gpp(1,nn) = gp2(j+13*nj)+g(j+2*nj)
            gpp(2,nn) = gp2(j+17*nj)
            gpp(3,nn) = gp2(j+7*nj)+4.0d+00*g(j+7*nj)
            gpp(4,nn) = gpp(2,nn)
            gpp(5,nn) = gp2(j+14*nj)+g(j+2*nj)
            gpp(6,nn) = gp2(j+8*nj)+4.0d+00*g(j+8*nj)
            gpp(7,nn) = gpp(3,nn)
            gpp(8,nn) = gpp(6,nn)
            gpp(9,nn) = gp2(j+2*nj)+9.0d+00*g(j+2*nj)+
     *                  12.0d+00*gm2(j+2*nj)
            dum = 1.0d+00
            if (norm) dum = sqrt7
            nn = nn+1
            gpp(1,nn) = dum*(gp2(j+3*nj)+7.0d+00*g(j+3*nj)+
     *                       6.0d+00*gm2(j+3*nj))
            gpp(2,nn) = dum*(gp2(j+9*nj)+g(j)+3.0d+00*g(j+9*nj)+
     *                       3.0d+00*gm2(j))
            gpp(3,nn) = dum*(gp2(j+15*nj)+3.0d+00*g(j+12*nj))
            gpp(4,nn) = gpp(2,nn)
            gpp(5,nn) = dum*(gp2(j+18*nj)+3.0d+00*g(j+3*nj))
            gpp(6,nn) = dum*(gp2(j+21*nj)+g(j+4*nj))
            gpp(7,nn) = gpp(3,nn)
            gpp(8,nn) = gpp(6,nn)
            gpp(9,nn) = dum*(gp2(j+22*nj)+g(j+3*nj))
            nn = nn+1
            gpp(1,nn) = dum*(gp2(j+4*nj)+7.0d+00*g(j+4*nj)+
     *                       6.0d+00*gm2(j+4*nj))
            gpp(2,nn) = dum*(gp2(j+15*nj)+3.0d+00*g(j+12*nj))
            gpp(3,nn) = dum*(gp2(j+10*nj)+g(j)+3.0d+00*g(j+10*nj)+
     *                       3.0d+00*gm2(j))
            gpp(4,nn) = gpp(2,nn)
            gpp(5,nn) = dum*(gp2(j+21*nj)+g(j+4*nj))
            gpp(6,nn) = dum*(gp2(j+22*nj)+g(j+3*nj))
            gpp(7,nn) = gpp(3,nn)
            gpp(8,nn) = gpp(6,nn)
            gpp(9,nn) = dum*(gp2(j+19*nj)+3.0d+00*g(j+4*nj))
            nn = nn+1
            gpp(1,nn) = dum*(gp2(j+18*nj)+3.0d+00*g(j+5*nj))
            gpp(2,nn) = dum*(gp2(j+11*nj)+g(j+nj)+3.0d+00*g(j+9*nj)+
     *                       3.0d+00*gm2(j+nj))
            gpp(3,nn) = dum*(gp2(j+23*nj)+g(j+6*nj))
            gpp(4,nn) = gpp(2,nn)
            gpp(5,nn) = dum*(gp2(j+5*nj)+7.0d+00*g(j+5*nj)+
     *                       6.0d+00*gm2(j+3*nj))
            gpp(6,nn) = dum*(gp2(j+16*nj)+3.0d+00*g(j+13*nj))
            gpp(7,nn) = gpp(3,nn)
            gpp(8,nn) = gpp(6,nn)
            gpp(9,nn) = dum*(gp2(j+24*nj)+g(j+5*nj))
            nn = nn+1
            gpp(1,nn) = dum*(gp2(j+23*nj)+g(j+6*nj))
            gpp(2,nn) = dum*(gp2(j+16*nj)+3.0d+00*g(j+13*nj))
            gpp(3,nn) = dum*(gp2(j+24*nj)+g(j+5*nj))
            gpp(4,nn) = gpp(2,nn)
            gpp(5,nn) = dum*(gp2(j+8*nj)+7.0d+00*g(j+6*nj)+
     *                       6.0d+00*gm2(j+5*nj))
            gpp(6,nn) = dum*(gp2(j+12*nj)+g(j+nj)+3.0d+00*g(j+11*nj)+
     *                       3.0d+00*gm2(j+nj))
            gpp(7,nn) = gpp(3,nn)
            gpp(8,nn) = gpp(6,nn)
            gpp(9,nn) = dum*(gp2(j+20*nj)+3.0d+00*g(j+6*nj))
            nn = nn+1
            gpp(1,nn) = dum*(gp2(j+19*nj)+3.0d+00*g(j+7*nj))
            gpp(2,nn) = dum*(gp2(j+25*nj)+g(j+8*nj))
            gpp(3,nn) = dum*(gp2(j+13*nj)+g(j+2*nj)+3.0d+00*g(j+10*nj)+
     *                       3.0d+00*gm2(j+2*nj))
            gpp(4,nn) = gpp(2,nn)
            gpp(5,nn) = dum*(gp2(j+26*nj)+g(j+7*nj))
            gpp(6,nn) = dum*(gp2(j+17*nj)+3.0d+00*g(j+14*nj))
            gpp(7,nn) = gpp(3,nn)
            gpp(8,nn) = gpp(6,nn)
            gpp(9,nn) = dum*(gp2(j+7*nj)+7.0d+00*g(j+7*nj)+
     *                       6.0d+00*gm2(j+4*nj))
            nn = nn+1
            gpp(1,nn) = dum*(gp2(j+25*nj)+g(j+8*nj))
            gpp(2,nn) = dum*(gp2(j+26*nj)+g(j+7*nj))
            gpp(3,nn) = dum*(gp2(j+17*nj)+3.0d+00*g(j+14*nj))
            gpp(4,nn) = gpp(2,nn)
            gpp(5,nn) = dum*(gp2(j+20*nj)+3.0d+00*g(j+8*nj))
            gpp(6,nn) = dum*(gp2(j+14*nj)+g(j+2*nj)+3.0d+00*g(j+11*nj)+
     *                       3.0d+00*gm2(j+2*nj))
            gpp(7,nn) = gpp(3,nn)
            gpp(8,nn) = gpp(6,nn)
            gpp(9,nn) = dum*(gp2(j+8*nj)+7.0d+00*g(j+8*nj)+
     *                       6.0d+00*gm2(j+5*nj))
            if (norm) dum = dum*sqrt5/sqrt3
            nn = nn+1
            gpp(1,nn) = dum*(gp2(j+9*nj)+5.0d+00*g(j+9*nj)+
     *                       2.0d+00*gm2(j+nj))
            gpp(2,nn) = dum*(gp2(j+18*nj)+2.0d+00*(g(j+3*nj)+
     *                       g(j+5*nj))+4.0d+00*gm2(j+3*nj))
            gpp(3,nn) = dum*(gp2(j+21*nj)+2.0d+00*g(j+13*nj))
            gpp(4,nn) = gpp(2,nn)
            gpp(5,nn) = dum*(gp2(j+11*nj)+5.0d+00*g(j+9*nj)+
     *                       2.0d+00*gm2(j))
            gpp(6,nn) = dum*(gp2(j+23*nj)+2.0d+00*g(j+12*nj))
            gpp(7,nn) = gpp(3,nn)
            gpp(8,nn) = gpp(6,nn)
            gpp(9,nn) = dum*(gp2(j+27*nj)+g(j+9*nj))
            nn = nn+1
            gpp(1,nn) = dum*(gp2(j+10*nj)+5.0d+00*g(j+10*nj)+
     *                       2.0d+00*gm2(j+2*nj))
            gpp(2,nn) = dum*(gp2(j+22*nj)+2.0d+00*g(j+14*nj))
            gpp(3,nn) = dum*(gp2(j+19*nj)+2.0d+00*(g(j+4*nj)+
     *                       g(j+7*nj))+4.0d+00*gm2(j+4*nj))
            gpp(4,nn) = gpp(2,nn)
            gpp(5,nn) = dum*(gp2(j+27*nj)+g(j+10*nj))
            gpp(6,nn) = dum*(gp2(j+25*nj)+2.0d+00*g(j+12*nj))
            gpp(7,nn) = gpp(3,nn)
            gpp(8,nn) = gpp(6,nn)
            gpp(9,nn) = dum*(gp2(j+13*nj)+5.0d+00*g(j+10*nj)+
     *                       2.0d+00*gm2(j))
            nn = nn+1
            gpp(1,nn) = dum*(gp2(j+27*nj)+g(j+11*nj))
            gpp(2,nn) = dum*(gp2(j+24*nj)+2.0d+00*g(j+14*nj))
            gpp(3,nn) = dum*(gp2(j+26*nj)+2.0d+00*g(j+13*nj))
            gpp(4,nn) = gpp(2,nn)
            gpp(5,nn) = dum*(gp2(j+12*nj)+5.0d+00*g(j+11*nj)+
     *                       2.0d+00*gm2(j+2*nj))
            gpp(6,nn) = dum*(gp2(j+20*nj)+2.0d+00*(g(j+6*nj)+
     *                       g(j+8*nj))+4.0d+00*gm2(j+5*nj))
            gpp(7,nn) = gpp(3,nn)
            gpp(8,nn) = gpp(6,nn)
            gpp(9,nn) = dum*(gp2(j+14*nj)+5.0d+00*g(j+11*nj)+
     *                       2.0d+00*gm2(j+nj))
            if (norm) dum = dum*sqrt3
            nn = nn+1
            gpp(1,nn) = dum*(gp2(j+15*nj)+5.0d+00*g(j+12*nj)+
     *                       2.0d+00*gm2(j+5*nj))
            gpp(2,nn) = dum*(gp2(j+21*nj)+2.0d+00*g(j+13*nj)+
     *                       g(j+4*nj)+2.0d+00*gm2(j+4*nj))
            gpp(3,nn) = dum*(gp2(j+22*nj)+2.0d+00*g(j+14*nj)+
     *                       g(j+3*nj)+2.0d+00*gm2(j+3*nj))
            gpp(4,nn) = gpp(2,nn)
            gpp(5,nn) = dum*(gp2(j+23*nj)+3.0d+00*g(j+12*nj))
            gpp(6,nn) = dum*(gp2(j+27*nj)+g(j+10*nj)+
     *                       g(j+9*nj)+gm2(j))
            gpp(7,nn) = gpp(3,nn)
            gpp(8,nn) = gpp(6,nn)
            gpp(9,nn) = dum*(gp2(j+25*nj)+3.0d+00*g(j+12*nj))
            nn = nn+1
            gpp(1,nn) = dum*(gp2(j+21*nj)+3.0d+00*g(j+13*nj))
            gpp(2,nn) = dum*(gp2(j+23*nj)+2.0d+00*g(j+12*nj)+
     *                       g(j+6*nj)+2.0d+00*gm2(j+5*nj))
            gpp(3,nn) = dum*(gp2(j+27*nj)+g(j+11*nj)+
     *                       g(j+9*nj)+gm2(j+nj))
            gpp(4,nn) = gpp(2,nn)
            gpp(5,nn) = dum*(gp2(j+16*nj)+5.0d+00*g(j+13*nj)+
     *                       2.0d+00*gm2(j+4*nj))
            gpp(6,nn) = dum*(gp2(j+24*nj)+2.0d+00*g(j+14*nj)+
     *                       g(j+5*nj)+2.0d+00*gm2(j+3*nj))
            gpp(7,nn) = gpp(3,nn)
            gpp(8,nn) = gpp(6,nn)
            gpp(9,nn) = dum*(gp2(j+26*nj)+3.0d+00*g(j+13*nj))
            nn = nn+1
            gpp(1,nn) = dum*(gp2(j+22*nj)+3.0d+00*g(j+14*nj))
            gpp(2,nn) = dum*(gp2(j+27*nj)+g(j+11*nj)+
     *                       g(j+10*nj)+gm2(j+2*nj))
            gpp(3,nn) = dum*(gp2(j+25*nj)+2.0d+00*g(j+12*nj)+
     *                       g(j+8*nj)+2.0d+00*gm2(j+5*nj))
            gpp(4,nn) = gpp(2,nn)
            gpp(5,nn) = dum*(gp2(j+24*nj)+3.0d+00*g(j+14*nj))
            gpp(6,nn) = dum*(gp2(j+26*nj)+2.0d+00*g(j+13*nj)+
     *                       g(j+7*nj)+2.0d+00*gm2(j+4*nj))
            gpp(7,nn) = gpp(3,nn)
            gpp(8,nn) = gpp(6,nn)
            gpp(9,nn) = dum*(gp2(j+17*nj)+5.0d+00*g(j+14*nj)+
     *                       2.0d+00*gm2(j+3*nj))
 500      continue
          return
        end if
      write(6,*) 'in formij h attempt',mini,maxi,minj,maxj
      call caserr('dimensioning problem in formij - contact authors')
      return
      end
      subroutine hssprt(nat,eg,eh)
c
      implicit real*8 (a-h,o-z)
      logical opg_root
c
      dimension eg(3,nat),eh(9,*)
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common /hsspar/ first,secnd,cphf,both,mfirst,msecnd,mcphf
c
      if(.not.opg_root()) return
c
      n = 0
      do 50 iatom = 1,nat
         do 40 jatom = 1,iatom
            n = n + 1
            write(iwr,110) iatom,jatom,eh(1,n),eh(2,n),eh(3,n)
            write(iwr,120)             eh(4,n),eh(5,n),eh(6,n)
            write(iwr,130)             eh(7,n),eh(8,n),eh(9,n)
   40    continue
   50 continue
      return
c
  110 format(1h ,5x,'atom(', 2i3,')  eh  :     d2/dx2  =',e18.10,
     +           5x,'d2/dydx =',e18.10,5x,'d2/dzdx =',e18.10)
  120 format(1h ,29x,                         'd2/dxdy =',e18.10,
     +           5x,'d2/dy2  =',e18.10,5x,'d2/dzdy =',e18.10)
  130 format(1h ,29x,                         'd2/dxdz =',e18.10,
     +           5x,'d2/dydz =',e18.10,5x,'d2/dz2  =',e18.10)
      end
      subroutine ver_integb(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/integb.m,v $
     +     "/
      data revision /"$Revision: 5774 $"/
      data date /"$Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
