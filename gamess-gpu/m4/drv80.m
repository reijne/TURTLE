c     deck=drv80
c ******************************************************
c ******************************************************
c             =   drv80  =
c ******************************************************
c ******************************************************
**==efill1.f
      subroutine efill1(ishell,jshell,kshell,lshell,la,lb,lc,ld,ax1,
     $                 dm,e,dmax,q4)
c
c***********************************************************************
c     routine to pluck density matrix contributions according
c     to shell numbers (ishell,...,lshell) and store into local
c     array e.  also, the max density matrix element for this block
c     is determined.
c
c     arguments
c
c     ishell to lshell ... shell numbers of the four shells.
c     dm               ... arrays containing density matrix.
c     e                ... array of dimension 256, filled with
c                          combinations of density matrix elements.
c     dmax             ... max combination formed.  useful later on
c                          when deciding just what to keep.
c
c***********************************************************************
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c
      dimension dm(*),e(*)
c
c     connect to the basis common.
c
INCLUDE(common/nshel)
INCLUDE(common/infoa)
INCLUDE(common/mapper)
_IF(ccpdft)
INCLUDE(common/ccpdft.hf77)
_ENDIF
c
      data dzero/0.0d0/,p25/0.25d0/
c
_IF(ccpdft)
      ohf_coul = (.not.CD_active()).or.CD_HF_coulomb_deriv()
      ohf_exch = (.not.CD_active()).or.CD_HF_exchange()
      hf_wght  = CD_HF_exchange_weight()
_ENDIF
c
c     determine starting addresses of basis functions within shells.
c
      iaos = kloc(ishell) - kmin(ishell)
      jaos = kloc(jshell) - kmin(jshell)
      kaos = kloc(kshell) - kmin(kshell)
      laos = kloc(lshell) - kmin(lshell)
c
c     determine starting addresses of basis functions within shells.
c
      dmax = dzero
      do 120 l = 1 , ld
         las = laos + l
c
         do 110 k = 1 , lc
            kas = kaos + k
            lkl = iky(max(kas,las)) + min(kas,las)
            d34 = dm(lkl)
c
            do 100 j = 1 , lb
               id = 16*j + 4*k + l - 84
               jas = jaos + j
               ljk = iky(max(jas,kas)) + min(jas,kas)
               d23 = dm(ljk)
               ljl = iky(max(jas,las)) + min(jas,las)
               d24 = dm(ljl)
c
               do 90 i = 1 , la
                  id = id + 64
                  ias = iaos + i
                  lij = iky(max(ias,jas)) + min(ias,jas)
                  d12 = dm(lij)
                  lik = iky(max(ias,kas)) + min(ias,kas)
                  d13 = dm(lik)
                  lil = iky(max(ias,las)) + min(ias,las)
                  d14 = dm(lil)
_IF(ccpdft)
                  d1234=0.0d0
                  if(ohf_coul)then
                     d1234 = d12*d34*ax1
                  endif
                  if(ohf_exch)then
                     wght = p25*hf_wght*ax1
                     d1234 = d1234 - wght*(d13*d24+d23*d14)
                  endif
_ELSE
                  d1234 = (d12*d34-p25*(d13*d24+d23*d14))*ax1
_ENDIF
c
c     store computed value into e, and check dmax.
c
                  e(id) = d1234*q4
                  d1234 = dabs(d1234)
                  if (d1234.gt.dmax) dmax = d1234
 90            continue
 100        continue
 110     continue
 120  continue
c
c     all done, return.
      return
      end
**==efill2.f
      subroutine efill2(ishell,jshell,kshell,lshell,la,lb,lc,ld,
     +ax1,dm,db,e,dmax,q4)
c
c***********************************************************************
c     routine to pluck density matrix contributions according
c     to shell numbers (ishell,...,lshell) and store into local
c     array e.  also, the max density matrix element for this block
c     is determined.
c
c     arguments
c
c     ishell to lshell ... shell numbers of the four shells.
c     iscf             ... scf mode flag.  see comments in
c     dm and dn        ... arrays containing uhf density matrices.
c     e                ... array of dimension 256, filled with
c                          combinations of density matrix elements.
c     dmax             ... max combination formed.  useful later on
c                          when deciding just what to keep.
c
c***********************************************************************
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c
      dimension dm(*),db(*),e(*)
c
c     connect to the basis common.
c
INCLUDE(common/nshel)
INCLUDE(common/infoa)
INCLUDE(common/mapper)
_IF(ccpdft)
INCLUDE(common/ccpdft.hf77)
_ENDIF
c
      data dzero/0.0d0/,p25/0.25d0/
c
_IF(ccpdft)
      ohf_coul = (.not.CD_active()).or.CD_HF_coulomb_deriv()
      ohf_exch = (.not.CD_active()).or.CD_HF_exchange()
      hf_wght  = CD_HF_exchange_weight()
_ENDIF
c
c     determine starting addresses of basis functions within shells.
c
      iaos = kloc(ishell) - kmin(ishell)
      jaos = kloc(jshell) - kmin(jshell)
      kaos = kloc(kshell) - kmin(kshell)
      laos = kloc(lshell) - kmin(lshell)
c
c     determine starting addresses of basis functions within shells.
c
      dmax = dzero
      do 120 l = 1 , ld
         las = laos + l
c
         do 110 k = 1 , lc
            kas = kaos + k
            lkl = iky(max(kas,las)) + min(kas,las)
            d34 = dm(lkl)
c
            do 100 j = 1 , lb
               id = 16*j + 4*k + l - 84
               jas = jaos + j
               ljk = iky(max(jas,kas)) + min(jas,kas)
               d23 = dm(ljk)
               ljl = iky(max(jas,las)) + min(jas,las)
               d24 = dm(ljl)
               d23b = db(ljk)
               d24b = db(ljl)
c
               do 90 i = 1 , la
                  id = id + 64
                  ias = iaos + i
                  lij = iky(max(ias,jas)) + min(ias,jas)
                  d12 = dm(lij)
                  lik = iky(max(ias,kas)) + min(ias,kas)
                  d13 = dm(lik)
                  lil = iky(max(ias,las)) + min(ias,las)
                  d14 = dm(lil)
c
                  d13b = db(lik)
                  d14b = db(lil)

_IF(ccpdft)
                  d1234=0.0d0
                  if(ohf_coul)then
                     d1234 = d12*d34*ax1
                  endif
                  if(ohf_exch)then
                     d1234 = d1234 - (d13*d24 + d23*d14
     +                     + d13b*d24b + d23b*d14b)*p25*hf_wght*ax1
                  endif
_ELSE
                  d1234 = ((d12)*(d34)
     +                    -p25*(d13*d24+d23*d14+d13b*d24b+d23b*d14b))
     +                    *ax1
_ENDIF
c
c     store computed value into e, and check dmax.
c
                  e(id) = d1234*q4
                  d1234 = dabs(d1234)
                  if (d1234.gt.dmax) dmax = d1234
 90            continue
 100        continue
 110     continue
 120  continue
c
      return
      end
**==efill3.f
      subroutine efill3(ishell,jshell,kshell,lshell,la,lb,lc,ld,
     $                 ax1,dm,db,e,dmax,q4,nconf)
c
c***********************************************************************
c     routine to pluck density matrix contributions according
c     to shell numbers (ishell,...,lshell) and store into local
c     array e.  also, the max density matrix element for this block
c     is determined.
c
c     arguments
c
c     ishell to lshell ... shell numbers of the four shells.
c     iscf             ... scf mode flag.  see comments in
c     dm and dn        ... arrays containing density matrices.
c                          exactly what they contain depends on
c                          iscf.  
c     e                ... array of dimension 256, filled with
c                          combinations of density matrix elements.
c     dmax             ... max combination formed.  useful later on
c                          when deciding just what to keep.
c
c***********************************************************************
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c
      dimension dm(*),db(*),e(*),nconf(*)
c
c     connect to the basis common.
c
INCLUDE(common/nshel)
INCLUDE(common/dmisc)
INCLUDE(common/scfwfn)
INCLUDE(common/infoa)
INCLUDE(common/mapper)
_IF(ccpdft)
INCLUDE(common/ccpdft.hf77)
_ENDIF
c
      data dzero/0.0d0/,h/0.5d0/,two/2.0d0/
c
c     determine starting addresses of basis functions within shells.
c
      iaos = kloc(ishell) - kmin(ishell)
      jaos = kloc(jshell) - kmin(jshell)
      kaos = kloc(kshell) - kmin(kshell)
      laos = kloc(lshell) - kmin(lshell)
c
c     determine starting addresses of basis functions within shells.
c

_IF(ccpdft)
c
c Safety check
c      
      if(CD_2e())then
         call caserr('DFT not coded for DRV80 + chosen SCFTYPE')
      endif
_ENDIF

      dmax = dzero
      do 120 l = 1 , ld
         las = laos + l
c
         do 110 k = 1 , lc
            kas = kaos + k
            lkl = iky(max(kas,las)) + min(kas,las)
            d34 = dm(lkl)
c
            do 100 j = 1 , lb
               id = 16*j + 4*k + l - 84
               jas = jaos + j
               ljk = iky(max(jas,kas)) + min(jas,kas)
               d23 = dm(ljk)
               ljl = iky(max(jas,las)) + min(jas,las)
               d24 = dm(ljl)
c
               do 90 i = 1 , la
                  id = id + 64
                  ias = iaos + i
                  lij = iky(max(ias,jas)) + min(ias,jas)
                  d12 = dm(lij)
                  lik = iky(max(ias,kas)) + min(ias,kas)
                  d13 = dm(lik)
                  lil = iky(max(ias,las)) + min(ias,las)
                  d14 = dm(lil)
c
c     gvb or grhf... generalised valence bond.
c
                  d1234 = dzero
                  if (.not.(onocor)) then
                     d1234 = d1234 + alpha(1)*d12*d34 + h*beta(1)
     +                       *(d13*d24+d14*d23)
                     if (onopen) then
                        d1234 = d1234*ax1*two
                        go to 80
                     else
                        nco1 = nco + 1
                        do 50 iox = nco1 , norb
                           io = num*(iox-1)
                           iojo = nconf(iox)*(nconf(iox)-1)/2 + 1
                           d1234 = d1234 + alpha(iojo)
     +                             *(d12*db(kas+io)*db(las+io)
     +                             +d34*db(ias+io)*db(jas+io))
     +                             + h*beta(iojo)
     +                             *(db(jas+io)*(d13*db(las+io)
     +                             +d14*db(kas+io))+db(ias+io)
     +                             *(d24*db(kas+io)+d23*db(las+io)))
 50                     continue
                     end if
                  end if
                  if (.not.(onopen)) then
                     nco1 = nco + 1
                     do 70 iox = nco1 , norb
                        io = num*(iox-1)
                        do 60 jox = nco1 , norb
                           jo = num*(jox-1)
                           iof = nconf(iox)
                           jof = nconf(jox)
                           iojo = max(iof,jof)*(max(iof,jof)-1)
     +                            /2 + min(iof,jof)
                           d1234 = d1234 + alpha(iojo)*db(ias+io)
     +                             *db(jas+io)*db(kas+jo)*db(las+jo)
     +                             + h*beta(iojo)*db(ias+io)*db(jas+jo)
     +                             *(db(kas+io)*db(las+jo)+db(las+io)
     +                             *db(kas+jo))
 60                     continue
 70                  continue
                  end if
                  d1234 = d1234*ax1*two
c
c     store computed value into e, and check dmax.
c
 80               e(id) = d1234*q4
                  d1234 = dabs(d1234)
                  if (d1234.gt.dmax) dmax = d1234
 90            continue
 100        continue
 110     continue
 120  continue
c
c     all done, return.
      return
      end
**==fpppp.f
      subroutine fpppp
c
      implicit REAL  (a-h,o-z)
INCLUDE(common/ffq)
INCLUDE(common/dindex)
c
      equivalence (gxyoo,gyxoo)
      equivalence (gxzoo,gzxoo)
      equivalence (gyzoo,gzyoo)
      equivalence (gyxxo,gxyxo,gxxyo)
      equivalence (gzxxo,gxzxo,gxxzo)
      equivalence (gzyyo,gyzyo,gyyzo)
      equivalence (gxyyo,gyxyo,gyyxo)
      equivalence (gxzzo,gzxzo,gzzxo)
      equivalence (gyzzo,gzyzo,gzzyo)
      equivalence (gxyzo,gyzxo,gzxyo,gzyxo,gyxzo,gxzyo)
      equivalence (gyyxx,gyxyx,gyxxy,gxyyx,gxyxy,gxxyy)
      equivalence (gzzxx,gzxzx,gzxxz,gxzzx,gxzxz,gxxzz)
      equivalence (gzzyy,gzyzy,gzyyz,gyzzy,gyzyz,gyyzz)
      equivalence (gyxxx,gxyxx,gxxyx,gxxxy)
      equivalence (gzxxx,gxzxx,gxxzx,gxxxz)
      equivalence (gxyyy,gyxyy,gyyxy,gyyyx)
      equivalence (gzyyy,gyzyy,gyyzy,gyyyz)
      equivalence (gxzzz,gzxzz,gzzxz,gzzzx)
      equivalence (gyzzz,gzyzz,gzzyz,gzzzy)
      equivalence (gxyzz,gxzyz,gxzzy,gyxzz,gyzxz,gyzzx,
     *             gzxyz,gzxzy,gzyxz,gzyzx,gzzxy,gzzyx)
      equivalence (gyzxx,gyxzx,gyxxz,gzyxx,gzxyx,gzxxy,
     *             gxyzx,gxyxz,gxzyx,gxzxy,gxxyz,gxxzy)
      equivalence (gzxyy,gzyxy,gzyyx,gxzyy,gxyzy,gxyyz,
     *             gyzxy,gyzyx,gyxzy,gyxyz,gyyzx,gyyxz)
c
      data three/3.0d0/
      data p25/0.25d0/
      data h/0.5d0/
c
      v0001 = qqox*goooo + gooxo
      v0002 = qqoy*goooo + gooyo
      v0003 = qqoz*goooo + goozo
      ve00 = ve00 + (v0001*e(2)+v0002*e(3)+v0003*e(4))*csssp
      temp = v0000*csssp
      ve14 = temp*e(2)
      ve24 = temp*e(3)
      ve34 = temp*e(4)
      ptqq = -a12*a34i
      ptqqs = ptqq**2
      gooxx = gxxoo*ptqqs
      gooyy = gyyoo*ptqqs
      goozz = gzzoo*ptqqs
      gooxy = gxyoo*ptqqs
      gooxz = gxzoo*ptqqs
      gooyz = gyzoo*ptqqs
      v0011 = qqxx*goooo + qqxox*gooxo + gooxx
      v0022 = qqyy*goooo + qqyoy*gooyo + gooyy
      v0033 = qqzz*goooo + qqzoz*goozo + goozz
      v0012 = qqxo*v0002 + qqoy*gooxo + gooxy
      v0013 = qqxo*v0003 + qqoz*gooxo + gooxz
      v0021 = qqyo*v0001 + qqox*gooyo + gooxy
      v0023 = qqyo*v0003 + qqoz*gooyo + gooyz
      v0031 = qqzo*v0001 + qqox*goozo + gooxz
      v0032 = qqzo*v0002 + qqoy*goozo + gooyz
      ve00 = ve00 + (v0011*e(6)+v0012*e(7)+v0013*e(8)+v0021*e(10)
     +       +v0022*e(11)+v0023*e(12)+v0031*e(14)+v0032*e(15)
     +       +v0033*e(16))*csspp
      ve14 = ve14 + (v0010*e(6)+v0020*e(10)+v0030*e(14))*csspp
      ve24 = ve24 + (v0010*e(7)+v0020*e(11)+v0030*e(15))*csspp
      ve34 = ve34 + (v0010*e(8)+v0020*e(12)+v0030*e(16))*csspp
      ve13 = ve13 + (v0001*e(6)+v0002*e(7)+v0003*e(8))*csspp
      ve23 = ve23 + (v0001*e(10)+v0002*e(11)+v0003*e(12))*csspp
      ve33 = ve33 + (v0001*e(14)+v0002*e(15)+v0003*e(16))*csspp
      v0101 = qqox*v0100 + ppox*gooxo + gxoxo
      v0102 = qqoy*v0100 + ppox*gooyo + gxoyo
      v0103 = qqoz*v0100 + ppox*goozo + gxozo
      v0201 = qqox*v0200 + ppoy*gooxo + gxoyo
      v0202 = qqoy*v0200 + ppoy*gooyo + gyoyo
      v0203 = qqoz*v0200 + ppoy*goozo + gyozo
      v0301 = qqox*v0300 + ppoz*gooxo + gxozo
      v0302 = qqoy*v0300 + ppoz*gooyo + gyozo
      v0303 = qqoz*v0300 + ppoz*goozo + gzozo
      ve00 = ve00 + (v0101*e(18)+v0102*e(19)+v0103*e(20)+v0201*e(34)
     +       +v0202*e(35)+v0203*e(36)+v0301*e(50)+v0302*e(51)
     +       +v0303*e(52))*cspsp
      ve14 = ve14 + (v0100*e(18)+v0200*e(34)+v0300*e(50))*cspsp
      ve24 = ve24 + (v0100*e(19)+v0200*e(35)+v0300*e(51))*cspsp
      ve34 = ve34 + (v0100*e(20)+v0200*e(36)+v0300*e(52))*cspsp
      ve12 = ve12 + (v0001*e(18)+v0002*e(19)+v0003*e(20))*cspsp
      ve22 = ve22 + (v0001*e(34)+v0002*e(35)+v0003*e(36))*cspsp
      ve32 = ve32 + (v0001*e(50)+v0002*e(51)+v0003*e(52))*cspsp
      v1001 = qqox*v1000 + ppxo*gooxo + gxoxo
      v1002 = qqoy*v1000 + ppxo*gooyo + gxoyo
      v1003 = qqoz*v1000 + ppxo*goozo + gxozo
      v2001 = qqox*v2000 + ppyo*gooxo + gxoyo
      v2002 = qqoy*v2000 + ppyo*gooyo + gyoyo
      v2003 = qqoz*v2000 + ppyo*goozo + gyozo
      v3001 = qqox*v3000 + ppzo*gooxo + gxozo
      v3002 = qqoy*v3000 + ppzo*gooyo + gyozo
      v3003 = qqoz*v3000 + ppzo*goozo + gzozo
      ve00 = ve00 + (v1001*e(66)+v1002*e(67)+v1003*e(68)+v2001*e(130)
     +       +v2002*e(131)+v2003*e(132)+v3001*e(194)+v3002*e(195)
     +       +v3003*e(196))*cpssp
      ve14 = ve14 + (v1000*e(66)+v2000*e(130)+v3000*e(194))*cpssp
      ve24 = ve24 + (v1000*e(67)+v2000*e(131)+v3000*e(195))*cpssp
      ve34 = ve34 + (v1000*e(68)+v2000*e(132)+v3000*e(196))*cpssp
      ve11 = ve11 + (v0001*e(66)+v0002*e(67)+v0003*e(68))*cpssp
      ve21 = ve21 + (v0001*e(130)+v0002*e(131)+v0003*e(132))*cpssp
      ve31 = ve31 + (v0001*e(194)+v0002*e(195)+v0003*e(196))*cpssp
      v1101 = qqox*v1100 + c1110
      v1102 = qqoy*v1100 + c1120
      v1103 = qqoz*v1100 + c1130
      v1201 = qqox*v1200 + c1210
      v1202 = qqoy*v1200 + c1220
      v1203 = qqoz*v1200 + c1230
      v1301 = qqox*v1300 + c1310
      v1302 = qqoy*v1300 + c1320
      v1303 = qqoz*v1300 + c1330
      v2101 = qqox*v2100 + c2110
      v2102 = qqoy*v2100 + c2120
      v2103 = qqoz*v2100 + c2130
      v2201 = qqox*v2200 + c2210
      v2202 = qqoy*v2200 + c2220
      v2203 = qqoz*v2200 + c2230
      v2301 = qqox*v2300 + c2310
      v2302 = qqoy*v2300 + c2320
      v2303 = qqoz*v2300 + c2330
      v3101 = qqox*v3100 + c3110
      v3102 = qqoy*v3100 + c3120
      v3103 = qqoz*v3100 + c3130
      v3201 = qqox*v3200 + c3210
      v3202 = qqoy*v3200 + c3220
      v3203 = qqoz*v3200 + c3230
      v3301 = qqox*v3300 + c3310
      v3302 = qqoy*v3300 + c3320
      v3303 = qqoz*v3300 + c3330
      ve00 = ve00 + (v1101*e(82)+v1102*e(83)+v1103*e(84)+v1201*e(98)
     +       +v1202*e(99)+v1203*e(100)+v1301*e(114)+v1302*e(115)
     +       +v1303*e(116)+v2101*e(146)+v2102*e(147)+v2103*e(148)
     +       +v2201*e(162)+v2202*e(163)+v2203*e(164)+v2301*e(178)
     +       +v2302*e(179)+v2303*e(180)+v3101*e(210)+v3102*e(211)
     +       +v3103*e(212)+v3201*e(226)+v3202*e(227)+v3203*e(228)
     +       +v3301*e(242)+v3302*e(243)+v3303*e(244))*cppsp
      ve14 = ve14 + (v1100*e(82)+v1200*e(98)+v1300*e(114)+v2100*e(146)
     +       +v2200*e(162)+v2300*e(178)+v3100*e(210)+v3200*e(226)
     +       +v3300*e(242))*cppsp
      ve24 = ve24 + (v1100*e(83)+v1200*e(99)+v1300*e(115)+v2100*e(147)
     +       +v2200*e(163)+v2300*e(179)+v3100*e(211)+v3200*e(227)
     +       +v3300*e(243))*cppsp
      ve34 = ve34 + (v1100*e(84)+v1200*e(100)+v1300*e(116)+v2100*e(148)
     +       +v2200*e(164)+v2300*e(180)+v3100*e(212)+v3200*e(228)
     +       +v3300*e(244))*cppsp
      ve12 = ve12 + (v1001*e(82)+v1002*e(83)+v1003*e(84)+v2001*e(146)
     +       +v2002*e(147)+v2003*e(148)+v3001*e(210)+v3002*e(211)
     +       +v3003*e(212))*cppsp
      ve22 = ve22 + (v1001*e(98)+v1002*e(99)+v1003*e(100)+v2001*e(162)
     +       +v2002*e(163)+v2003*e(164)+v3001*e(226)+v3002*e(227)
     +       +v3003*e(228))*cppsp
      ve32 = ve32 + (v1001*e(114)+v1002*e(115)+v1003*e(116)+v2001*e(178)
     +       +v2002*e(179)+v2003*e(180)+v3001*e(242)+v3002*e(243)
     +       +v3003*e(244))*cppsp
      ve11 = ve11 + (v0101*e(82)+v0102*e(83)+v0103*e(84)+v0201*e(98)
     +       +v0202*e(99)+v0203*e(100)+v0301*e(114)+v0302*e(115)
     +       +v0303*e(116))*cppsp
      ve21 = ve21 + (v0101*e(146)+v0102*e(147)+v0103*e(148)+v0201*e(162)
     +       +v0202*e(163)+v0203*e(164)+v0301*e(178)+v0302*e(179)
     +       +v0303*e(180))*cppsp
      ve31 = ve31 + (v0101*e(210)+v0102*e(211)+v0103*e(212)+v0201*e(226)
     +       +v0202*e(227)+v0203*e(228)+v0301*e(242)+v0302*e(243)
     +       +v0303*e(244))*cppsp
      gxoyz = gxyzo*ptqq
      gxoxx = gxxxo*ptqq
      gyoyy = gyyyo*ptqq
      gzozz = gzzzo*ptqq
      gyoxx = gxxyo*ptqq
      gzoxx = gxxzo*ptqq
      gxoyy = gyyxo*ptqq
      gzoyy = gyyzo*ptqq
      gxozz = gzzxo*ptqq
      gyozz = gzzyo*ptqq
      c1011 = qqxx*gxooo + qqxox*gxoxo + gxoxx
      c1022 = qqyy*gxooo + qqyoy*gxoyo + gxoyy
      c1033 = qqzz*gxooo + qqzoz*gxozo + gxozz
      c2011 = qqxx*gyooo + qqxox*gxoyo + gyoxx
      c2022 = qqyy*gyooo + qqyoy*gyoyo + gyoyy
      c2033 = qqzz*gyooo + qqzoz*gyozo + gyozz
      c3011 = qqxx*gzooo + qqxox*gxozo + gzoxx
      c3022 = qqyy*gzooo + qqyoy*gyozo + gzoyy
      c3033 = qqzz*gzooo + qqzoz*gzozo + gzozz
      c1012 = qqxy*gxooo + qqxo*gxoyo + qqoy*gxoxo + gyoxx
      c1013 = qqxz*gxooo + qqxo*gxozo + qqoz*gxoxo + gzoxx
      c1021 = qqyx*gxooo + qqyo*gxoxo + qqox*gxoyo + gyoxx
      c1023 = qqyz*gxooo + qqyo*gxozo + qqoz*gxoyo + gxoyz
      c1031 = qqzx*gxooo + qqzo*gxoxo + qqox*gxozo + gzoxx
      c1032 = qqzy*gxooo + qqzo*gxoyo + qqoy*gxozo + gxoyz
      c2012 = qqxy*gyooo + qqxo*gyoyo + qqoy*gxoyo + gxoyy
      c2013 = qqxz*gyooo + qqxo*gyozo + qqoz*gxoyo + gxoyz
      c2021 = qqyx*gyooo + qqyo*gxoyo + qqox*gyoyo + gxoyy
      c2023 = qqyz*gyooo + qqyo*gyozo + qqoz*gyoyo + gzoyy
      c2031 = qqzx*gyooo + qqzo*gxoyo + qqox*gyozo + gxoyz
      c2032 = qqzy*gyooo + qqzo*gyoyo + qqoy*gyozo + gzoyy
      c3012 = qqxy*gzooo + qqxo*gyozo + qqoy*gxozo + gxoyz
      c3013 = qqxz*gzooo + qqxo*gzozo + qqoz*gxozo + gxozz
      c3021 = qqyx*gzooo + qqyo*gxozo + qqox*gyozo + gxoyz
      c3023 = qqyz*gzooo + qqyo*gzozo + qqoz*gyozo + gyozz
      c3031 = qqzx*gzooo + qqzo*gxozo + qqox*gzozo + gxozz
      c3032 = qqzy*gzooo + qqzo*gyozo + qqoy*gzozo + gyozz
      v0111 = ppox*v0011 + c1011
      v0112 = ppox*v0012 + c1012
      v0113 = ppox*v0013 + c1013
      v0121 = ppox*v0021 + c1021
      v0122 = ppox*v0022 + c1022
      v0123 = ppox*v0023 + c1023
      v0131 = ppox*v0031 + c1031
      v0132 = ppox*v0032 + c1032
      v0133 = ppox*v0033 + c1033
      v0211 = ppoy*v0011 + c2011
      v0212 = ppoy*v0012 + c2012
      v0213 = ppoy*v0013 + c2013
      v0221 = ppoy*v0021 + c2021
      v0222 = ppoy*v0022 + c2022
      v0223 = ppoy*v0023 + c2023
      v0231 = ppoy*v0031 + c2031
      v0232 = ppoy*v0032 + c2032
      v0233 = ppoy*v0033 + c2033
      v0311 = ppoz*v0011 + c3011
      v0312 = ppoz*v0012 + c3012
      v0313 = ppoz*v0013 + c3013
      v0321 = ppoz*v0021 + c3021
      v0322 = ppoz*v0022 + c3022
      v0323 = ppoz*v0023 + c3023
      v0331 = ppoz*v0031 + c3031
      v0332 = ppoz*v0032 + c3032
      v0333 = ppoz*v0033 + c3033
      ve00 = ve00 + (v0111*e(22)+v0112*e(23)+v0113*e(24)+v0121*e(26)
     +       +v0122*e(27)+v0123*e(28)+v0131*e(30)+v0132*e(31)
     +       +v0133*e(32)+v0211*e(38)+v0212*e(39)+v0213*e(40)
     +       +v0221*e(42)+v0222*e(43)+v0223*e(44)+v0231*e(46)
     +       +v0232*e(47)+v0233*e(48)+v0311*e(54)+v0312*e(55)
     +       +v0313*e(56)+v0321*e(58)+v0322*e(59)+v0323*e(60)
     +       +v0331*e(62)+v0332*e(63)+v0333*e(64))*csppp
      ve14 = ve14 + (v0110*e(22)+v0120*e(26)+v0130*e(30)+v0210*e(38)
     +       +v0220*e(42)+v0230*e(46)+v0310*e(54)+v0320*e(58)
     +       +v0330*e(62))*csppp
      ve24 = ve24 + (v0110*e(23)+v0120*e(27)+v0130*e(31)+v0210*e(39)
     +       +v0220*e(43)+v0230*e(47)+v0310*e(55)+v0320*e(59)
     +       +v0330*e(63))*csppp
      ve34 = ve34 + (v0110*e(24)+v0120*e(28)+v0130*e(32)+v0210*e(40)
     +       +v0220*e(44)+v0230*e(48)+v0310*e(56)+v0320*e(60)
     +       +v0330*e(64))*csppp
      ve13 = ve13 + (v0101*e(22)+v0102*e(23)+v0103*e(24)+v0201*e(38)
     +       +v0202*e(39)+v0203*e(40)+v0301*e(54)+v0302*e(55)
     +       +v0303*e(56))*csppp
      ve23 = ve23 + (v0101*e(26)+v0102*e(27)+v0103*e(28)+v0201*e(42)
     +       +v0202*e(43)+v0203*e(44)+v0301*e(58)+v0302*e(59)
     +       +v0303*e(60))*csppp
      ve33 = ve33 + (v0101*e(30)+v0102*e(31)+v0103*e(32)+v0201*e(46)
     +       +v0202*e(47)+v0203*e(48)+v0301*e(62)+v0302*e(63)
     +       +v0303*e(64))*csppp
      ve12 = ve12 + (v0011*e(22)+v0012*e(23)+v0013*e(24)+v0021*e(26)
     +       +v0022*e(27)+v0023*e(28)+v0031*e(30)+v0032*e(31)
     +       +v0033*e(32))*csppp
      ve22 = ve22 + (v0011*e(38)+v0012*e(39)+v0013*e(40)+v0021*e(42)
     +       +v0022*e(43)+v0023*e(44)+v0031*e(46)+v0032*e(47)
     +       +v0033*e(48))*csppp
      ve32 = ve32 + (v0011*e(54)+v0012*e(55)+v0013*e(56)+v0021*e(58)
     +       +v0022*e(59)+v0023*e(60)+v0031*e(62)+v0032*e(63)
     +       +v0033*e(64))*csppp
      v1011 = ppxo*v0011 + c1011
      v1012 = ppxo*v0012 + c1012
      v1013 = ppxo*v0013 + c1013
      v1021 = ppxo*v0021 + c1021
      v1022 = ppxo*v0022 + c1022
      v1023 = ppxo*v0023 + c1023
      v1031 = ppxo*v0031 + c1031
      v1032 = ppxo*v0032 + c1032
      v1033 = ppxo*v0033 + c1033
      v2011 = ppyo*v0011 + c2011
      v2012 = ppyo*v0012 + c2012
      v2013 = ppyo*v0013 + c2013
      v2021 = ppyo*v0021 + c2021
      v2022 = ppyo*v0022 + c2022
      v2023 = ppyo*v0023 + c2023
      v2031 = ppyo*v0031 + c2031
      v2032 = ppyo*v0032 + c2032
      v2033 = ppyo*v0033 + c2033
      v3011 = ppzo*v0011 + c3011
      v3012 = ppzo*v0012 + c3012
      v3013 = ppzo*v0013 + c3013
      v3021 = ppzo*v0021 + c3021
      v3022 = ppzo*v0022 + c3022
      v3023 = ppzo*v0023 + c3023
      v3031 = ppzo*v0031 + c3031
      v3032 = ppzo*v0032 + c3032
      v3033 = ppzo*v0033 + c3033
      ve00 = ve00 + (v1011*e(70)+v1012*e(71)+v1013*e(72)+v1021*e(74)
     +       +v1022*e(75)+v1023*e(76)+v1031*e(78)+v1032*e(79)
     +       +v1033*e(80)+v2011*e(134)+v2012*e(135)+v2013*e(136)
     +       +v2021*e(138)+v2022*e(139)+v2023*e(140)+v2031*e(142)
     +       +v2032*e(143)+v2033*e(144)+v3011*e(198)+v3012*e(199)
     +       +v3013*e(200)+v3021*e(202)+v3022*e(203)+v3023*e(204)
     +       +v3031*e(206)+v3032*e(207)+v3033*e(208))*cpspp
      ve14 = ve14 + (v1010*e(70)+v1020*e(74)+v1030*e(78)+v2010*e(134)
     +       +v2020*e(138)+v2030*e(142)+v3010*e(198)+v3020*e(202)
     +       +v3030*e(206))*cpspp
      ve24 = ve24 + (v1010*e(71)+v1020*e(75)+v1030*e(79)+v2010*e(135)
     +       +v2020*e(139)+v2030*e(143)+v3010*e(199)+v3020*e(203)
     +       +v3030*e(207))*cpspp
      ve34 = ve34 + (v1010*e(72)+v1020*e(76)+v1030*e(80)+v2010*e(136)
     +       +v2020*e(140)+v2030*e(144)+v3010*e(200)+v3020*e(204)
     +       +v3030*e(208))*cpspp
      ve13 = ve13 + (v1001*e(70)+v1002*e(71)+v1003*e(72)+v2001*e(134)
     +       +v2002*e(135)+v2003*e(136)+v3001*e(198)+v3002*e(199)
     +       +v3003*e(200))*cpspp
      ve23 = ve23 + (v1001*e(74)+v1002*e(75)+v1003*e(76)+v2001*e(138)
     +       +v2002*e(139)+v2003*e(140)+v3001*e(202)+v3002*e(203)
     +       +v3003*e(204))*cpspp
      ve33 = ve33 + (v1001*e(78)+v1002*e(79)+v1003*e(80)+v2001*e(142)
     +       +v2002*e(143)+v2003*e(144)+v3001*e(206)+v3002*e(207)
     +       +v3003*e(208))*cpspp
      ve11 = ve11 + (v0011*e(70)+v0012*e(71)+v0013*e(72)+v0021*e(74)
     +       +v0022*e(75)+v0023*e(76)+v0031*e(78)+v0032*e(79)
     +       +v0033*e(80))*cpspp
      ve21 = ve21 + (v0011*e(134)+v0012*e(135)+v0013*e(136)+v0021*e(138)
     +       +v0022*e(139)+v0023*e(140)+v0031*e(142)+v0032*e(143)
     +       +v0033*e(144))*cpspp
      ve31 = ve31 + (v0011*e(198)+v0012*e(199)+v0013*e(200)+v0021*e(202)
     +       +v0022*e(203)+v0023*e(204)+v0031*e(206)+v0032*e(207)
     +       +v0033*e(208))*cpspp
      qfq4 = (qa1*qa2)**2*fq4
      temp = a1234i**2
      hfq3 = h*qa*fq3*temp
      tfq3 = three*hfq3
      sfq3 = tfq3 + tfq3
      p25fq2 = p25*fq2*temp
      p75fq2 = three*p25fq2
      temp = pqxx*qfq4
      gxxxx = pqxx*(temp-sfq3) + p75fq2
      gyxxx = pqxy*(temp-tfq3)
      gzxxx = pqxz*(temp-tfq3)
      gzyxx = pqyz*(temp-hfq3)
      gyyxx = pqyy*(temp-hfq3) - pqxx*hfq3 + p25fq2
      temp = pqyy*qfq4
      gyyyy = pqyy*(temp-sfq3) + p75fq2
      gxyyy = pqxy*(temp-tfq3)
      gzxyy = pqxz*(temp-hfq3)
      gzyyy = pqyz*(temp-tfq3)
      gzzyy = pqzz*(temp-hfq3) - pqyy*hfq3 + p25fq2
      temp = pqzz*qfq4
      gzzzz = pqzz*(temp-sfq3) + p75fq2
      gxyzz = pqxy*(temp-hfq3)
      gxzzz = pqxz*(temp-tfq3)
      gyzzz = pqyz*(temp-tfq3)
      gxxzz = pqxx*(temp-hfq3) - pqzz*hfq3 + p25fq2
c
      vp4 = (ppxx*v0011+ppxox*c1011+qqxx*gxxoo+qqxox*gxxxo+gxxxx)*e(86)
      vp4 = (ppxx*v0012+ppxox*c1012+qqxy*gxxoo+qqxo*gxxyo+qqoy*gxxxo+
     +      gxxxy)*e(87) + vp4
      vp4 = (ppxx*v0013+ppxox*c1013+qqxz*gxxoo+qqxo*gxxzo+qqoz*gxxxo+
     +      gxxxz)*e(88) + vp4
      vp4 = (ppxx*v0021+ppxox*c1021+qqyx*gxxoo+qqyo*gxxxo+qqox*gxxyo+
     +      gxxyx)*e(90) + vp4
      vp4 = (ppxx*v0022+ppxox*c1022+qqyy*gxxoo+qqyoy*gxxyo+gxxyy)*e(91)
     +      + vp4
      vp4 = (ppxx*v0023+ppxox*c1023+qqyz*gxxoo+qqyo*gxxzo+qqoz*gxxyo+
     +      gxxyz)*e(92) + vp4
      vp4 = (ppxx*v0031+ppxox*c1031+qqzx*gxxoo+qqzo*gxxxo+qqox*gxxzo+
     +      gxxzx)*e(94) + vp4
      vp4 = (ppxx*v0032+ppxox*c1032+qqzy*gxxoo+qqzo*gxxyo+qqoy*gxxzo+
     +      gxxzy)*e(95) + vp4
      vp4 = (ppxx*v0033+ppxox*c1033+qqzz*gxxoo+qqzoz*gxxzo+gxxzz)*e(96)
     +      + vp4
      vp4 = (ppxy*v0011+ppxo*c2011+ppoy*c1011+qqxx*gxyoo+qqxox*gxyxo+
     +      gxyxx)*e(102) + vp4
      vp4 = (ppxy*v0012+ppxo*c2012+ppoy*c1012+qqxy*gxyoo+qqxo*gxyyo+
     +      qqoy*gxyxo+gxyxy)*e(103) + vp4
      vp4 = (ppxy*v0013+ppxo*c2013+ppoy*c1013+qqxz*gxyoo+qqxo*gxyzo+
     +      qqoz*gxyxo+gxyxz)*e(104) + vp4
      vp4 = (ppxy*v0021+ppxo*c2021+ppoy*c1021+qqyx*gxyoo+qqyo*gxyxo+
     +      qqox*gxyyo+gxyyx)*e(106) + vp4
      vp4 = (ppxy*v0022+ppxo*c2022+ppoy*c1022+qqyy*gxyoo+qqyoy*gxyyo+
     +      gxyyy)*e(107) + vp4
      vp4 = (ppxy*v0023+ppxo*c2023+ppoy*c1023+qqyz*gxyoo+qqyo*gxyzo+
     +      qqoz*gxyyo+gxyyz)*e(108) + vp4
      vp4 = (ppxy*v0031+ppxo*c2031+ppoy*c1031+qqzx*gxyoo+qqzo*gxyxo+
     +      qqox*gxyzo+gxyzx)*e(110) + vp4
      vp4 = (ppxy*v0032+ppxo*c2032+ppoy*c1032+qqzy*gxyoo+qqzo*gxyyo+
     +      qqoy*gxyzo+gxyzy)*e(111) + vp4
      vp4 = (ppxy*v0033+ppxo*c2033+ppoy*c1033+qqzz*gxyoo+qqzoz*gxyzo+
     +      gxyzz)*e(112) + vp4
      vp4 = (ppxz*v0011+ppxo*c3011+ppoz*c1011+qqxx*gxzoo+qqxox*gxzxo+
     +      gxzxx)*e(118) + vp4
      vp4 = (ppxz*v0012+ppxo*c3012+ppoz*c1012+qqxy*gxzoo+qqxo*gxzyo+
     +      qqoy*gxzxo+gxzxy)*e(119) + vp4
      vp4 = (ppxz*v0013+ppxo*c3013+ppoz*c1013+qqxz*gxzoo+qqxo*gxzzo+
     +      qqoz*gxzxo+gxzxz)*e(120) + vp4
      vp4 = (ppxz*v0021+ppxo*c3021+ppoz*c1021+qqyx*gxzoo+qqyo*gxzxo+
     +      qqox*gxzyo+gxzyx)*e(122) + vp4
      vp4 = (ppxz*v0022+ppxo*c3022+ppoz*c1022+qqyy*gxzoo+qqyoy*gxzyo+
     +      gxzyy)*e(123) + vp4
      vp4 = (ppxz*v0023+ppxo*c3023+ppoz*c1023+qqyz*gxzoo+qqyo*gxzzo+
     +      qqoz*gxzyo+gxzyz)*e(124) + vp4
      vp4 = (ppxz*v0031+ppxo*c3031+ppoz*c1031+qqzx*gxzoo+qqzo*gxzxo+
     +      qqox*gxzzo+gxzzx)*e(126) + vp4
      vp4 = (ppxz*v0032+ppxo*c3032+ppoz*c1032+qqzy*gxzoo+qqzo*gxzyo+
     +      qqoy*gxzzo+gxzzy)*e(127) + vp4
      vp4 = (ppxz*v0033+ppxo*c3033+ppoz*c1033+qqzz*gxzoo+qqzoz*gxzzo+
     +      gxzzz)*e(128) + vp4
      vp4 = (ppyx*v0011+ppyo*c1011+ppox*c2011+qqxx*gyxoo+qqxox*gyxxo+
     +      gyxxx)*e(150) + vp4
      vp4 = (ppyx*v0012+ppyo*c1012+ppox*c2012+qqxy*gyxoo+qqxo*gyxyo+
     +      qqoy*gyxxo+gyxxy)*e(151) + vp4
      vp4 = (ppyx*v0013+ppyo*c1013+ppox*c2013+qqxz*gyxoo+qqxo*gyxzo+
     +      qqoz*gyxxo+gyxxz)*e(152) + vp4
      vp4 = (ppyx*v0021+ppyo*c1021+ppox*c2021+qqyx*gyxoo+qqyo*gyxxo+
     +      qqox*gyxyo+gyxyx)*e(154) + vp4
      vp4 = (ppyx*v0022+ppyo*c1022+ppox*c2022+qqyy*gyxoo+qqyoy*gyxyo+
     +      gyxyy)*e(155) + vp4
      vp4 = (ppyx*v0023+ppyo*c1023+ppox*c2023+qqyz*gyxoo+qqyo*gyxzo+
     +      qqoz*gyxyo+gyxyz)*e(156) + vp4
      vp4 = (ppyx*v0031+ppyo*c1031+ppox*c2031+qqzx*gyxoo+qqzo*gyxxo+
     +      qqox*gyxzo+gyxzx)*e(158) + vp4
      vp4 = (ppyx*v0032+ppyo*c1032+ppox*c2032+qqzy*gyxoo+qqzo*gyxyo+
     +      qqoy*gyxzo+gyxzy)*e(159) + vp4
      vp4 = (ppyx*v0033+ppyo*c1033+ppox*c2033+qqzz*gyxoo+qqzoz*gyxzo+
     +      gyxzz)*e(160) + vp4
      vp4 = (ppyy*v0011+ppyoy*c2011+qqxx*gyyoo+qqxox*gyyxo+gyyxx)*e(166)
     +      + vp4
      vp4 = (ppyy*v0012+ppyoy*c2012+qqxy*gyyoo+qqxo*gyyyo+qqoy*gyyxo+
     +      gyyxy)*e(167) + vp4
      vp4 = (ppyy*v0013+ppyoy*c2013+qqxz*gyyoo+qqxo*gyyzo+qqoz*gyyxo+
     +      gyyxz)*e(168) + vp4
      vp4 = (ppyy*v0021+ppyoy*c2021+qqyx*gyyoo+qqyo*gyyxo+qqox*gyyyo+
     +      gyyyx)*e(170) + vp4
      vp4 = (ppyy*v0022+ppyoy*c2022+qqyy*gyyoo+qqyoy*gyyyo+gyyyy)*e(171)
     +      + vp4
      vp4 = (ppyy*v0023+ppyoy*c2023+qqyz*gyyoo+qqyo*gyyzo+qqoz*gyyyo+
     +      gyyyz)*e(172) + vp4
      vp4 = (ppyy*v0031+ppyoy*c2031+qqzx*gyyoo+qqzo*gyyxo+qqox*gyyzo+
     +      gyyzx)*e(174) + vp4
      vp4 = (ppyy*v0032+ppyoy*c2032+qqzy*gyyoo+qqzo*gyyyo+qqoy*gyyzo+
     +      gyyzy)*e(175) + vp4
      vp4 = (ppyy*v0033+ppyoy*c2033+qqzz*gyyoo+qqzoz*gyyzo+gyyzz)*e(176)
     +      + vp4
      vp4 = (ppyz*v0011+ppyo*c3011+ppoz*c2011+qqxx*gyzoo+qqxox*gyzxo+
     +      gyzxx)*e(182) + vp4
      vp4 = (ppyz*v0012+ppyo*c3012+ppoz*c2012+qqxy*gyzoo+qqxo*gyzyo+
     +      qqoy*gyzxo+gyzxy)*e(183) + vp4
      vp4 = (ppyz*v0013+ppyo*c3013+ppoz*c2013+qqxz*gyzoo+qqxo*gyzzo+
     +      qqoz*gyzxo+gyzxz)*e(184) + vp4
      vp4 = (ppyz*v0021+ppyo*c3021+ppoz*c2021+qqyx*gyzoo+qqyo*gyzxo+
     +      qqox*gyzyo+gyzyx)*e(186) + vp4
      vp4 = (ppyz*v0022+ppyo*c3022+ppoz*c2022+qqyy*gyzoo+qqyoy*gyzyo+
     +      gyzyy)*e(187) + vp4
      vp4 = (ppyz*v0023+ppyo*c3023+ppoz*c2023+qqyz*gyzoo+qqyo*gyzzo+
     +      qqoz*gyzyo+gyzyz)*e(188) + vp4
      vp4 = (ppyz*v0031+ppyo*c3031+ppoz*c2031+qqzx*gyzoo+qqzo*gyzxo+
     +      qqox*gyzzo+gyzzx)*e(190) + vp4
      vp4 = (ppyz*v0032+ppyo*c3032+ppoz*c2032+qqzy*gyzoo+qqzo*gyzyo+
     +      qqoy*gyzzo+gyzzy)*e(191) + vp4
      vp4 = (ppyz*v0033+ppyo*c3033+ppoz*c2033+qqzz*gyzoo+qqzoz*gyzzo+
     +      gyzzz)*e(192) + vp4
      vp4 = (ppzx*v0011+ppzo*c1011+ppox*c3011+qqxx*gzxoo+qqxox*gzxxo+
     +      gzxxx)*e(214) + vp4
      vp4 = (ppzx*v0012+ppzo*c1012+ppox*c3012+qqxy*gzxoo+qqxo*gzxyo+
     +      qqoy*gzxxo+gzxxy)*e(215) + vp4
      vp4 = (ppzx*v0013+ppzo*c1013+ppox*c3013+qqxz*gzxoo+qqxo*gzxzo+
     +      qqoz*gzxxo+gzxxz)*e(216) + vp4
      vp4 = (ppzx*v0021+ppzo*c1021+ppox*c3021+qqyx*gzxoo+qqyo*gzxxo+
     +      qqox*gzxyo+gzxyx)*e(218) + vp4
      vp4 = (ppzx*v0022+ppzo*c1022+ppox*c3022+qqyy*gzxoo+qqyoy*gzxyo+
     +      gzxyy)*e(219) + vp4
      vp4 = (ppzx*v0023+ppzo*c1023+ppox*c3023+qqyz*gzxoo+qqyo*gzxzo+
     +      qqoz*gzxyo+gzxyz)*e(220) + vp4
      vp4 = (ppzx*v0031+ppzo*c1031+ppox*c3031+qqzx*gzxoo+qqzo*gzxxo+
     +      qqox*gzxzo+gzxzx)*e(222) + vp4
      vp4 = (ppzx*v0032+ppzo*c1032+ppox*c3032+qqzy*gzxoo+qqzo*gzxyo+
     +      qqoy*gzxzo+gzxzy)*e(223) + vp4
      vp4 = (ppzx*v0033+ppzo*c1033+ppox*c3033+qqzz*gzxoo+qqzoz*gzxzo+
     +      gzxzz)*e(224) + vp4
      vp4 = (ppzy*v0011+ppzo*c2011+ppoy*c3011+qqxx*gzyoo+qqxox*gzyxo+
     +      gzyxx)*e(230) + vp4
      vp4 = (ppzy*v0012+ppzo*c2012+ppoy*c3012+qqxy*gzyoo+qqxo*gzyyo+
     +      qqoy*gzyxo+gzyxy)*e(231) + vp4
      vp4 = (ppzy*v0013+ppzo*c2013+ppoy*c3013+qqxz*gzyoo+qqxo*gzyzo+
     +      qqoz*gzyxo+gzyxz)*e(232) + vp4
      vp4 = (ppzy*v0021+ppzo*c2021+ppoy*c3021+qqyx*gzyoo+qqyo*gzyxo+
     +      qqox*gzyyo+gzyyx)*e(234) + vp4
      vp4 = (ppzy*v0022+ppzo*c2022+ppoy*c3022+qqyy*gzyoo+qqyoy*gzyyo+
     +      gzyyy)*e(235) + vp4
      vp4 = (ppzy*v0023+ppzo*c2023+ppoy*c3023+qqyz*gzyoo+qqyo*gzyzo+
     +      qqoz*gzyyo+gzyyz)*e(236) + vp4
      vp4 = (ppzy*v0031+ppzo*c2031+ppoy*c3031+qqzx*gzyoo+qqzo*gzyxo+
     +      qqox*gzyzo+gzyzx)*e(238) + vp4
      vp4 = (ppzy*v0032+ppzo*c2032+ppoy*c3032+qqzy*gzyoo+qqzo*gzyyo+
     +      qqoy*gzyzo+gzyzy)*e(239) + vp4
      vp4 = (ppzy*v0033+ppzo*c2033+ppoy*c3033+qqzz*gzyoo+qqzoz*gzyzo+
     +      gzyzz)*e(240) + vp4
      vp4 = (ppzz*v0011+ppzoz*c3011+qqxx*gzzoo+qqxox*gzzxo+gzzxx)*e(246)
     +      + vp4
      vp4 = (ppzz*v0012+ppzoz*c3012+qqxy*gzzoo+qqxo*gzzyo+qqoy*gzzxo+
     +      gzzxy)*e(247) + vp4
      vp4 = (ppzz*v0013+ppzoz*c3013+qqxz*gzzoo+qqxo*gzzzo+qqoz*gzzxo+
     +      gzzxz)*e(248) + vp4
      vp4 = (ppzz*v0021+ppzoz*c3021+qqyx*gzzoo+qqyo*gzzxo+qqox*gzzyo+
     +      gzzyx)*e(250) + vp4
      vp4 = (ppzz*v0022+ppzoz*c3022+qqyy*gzzoo+qqyoy*gzzyo+gzzyy)*e(251)
     +      + vp4
      vp4 = (ppzz*v0023+ppzoz*c3023+qqyz*gzzoo+qqyo*gzzzo+qqoz*gzzyo+
     +      gzzyz)*e(252) + vp4
      vp4 = (ppzz*v0031+ppzoz*c3031+qqzx*gzzoo+qqzo*gzzxo+qqox*gzzzo+
     +      gzzzx)*e(254) + vp4
      vp4 = (ppzz*v0032+ppzoz*c3032+qqzy*gzzoo+qqzo*gzzyo+qqoy*gzzzo+
     +      gzzzy)*e(255) + vp4
      vp4 = (ppzz*v0033+ppzoz*c3033+qqzz*gzzoo+qqzoz*gzzzo+gzzzz)*e(256)
     +      + vp4
      vpppp = vp4*cpppp
      ve00 = ve00 + vpppp
      ve14 = ve14 + (v1110*e(86)+v1120*e(90)+v1130*e(94)+v1210*e(102)
     +       +v1220*e(106)+v1230*e(110)+v1310*e(118)+v1320*e(122)
     +       +v1330*e(126)+v2110*e(150)+v2120*e(154)+v2130*e(158)
     +       +v2210*e(166)+v2220*e(170)+v2230*e(174)+v2310*e(182)
     +       +v2320*e(186)+v2330*e(190)+v3110*e(214)+v3120*e(218)
     +       +v3130*e(222)+v3210*e(230)+v3220*e(234)+v3230*e(238)
     +       +v3310*e(246)+v3320*e(250)+v3330*e(254))*cpppp
      ve24 = ve24 + (v1110*e(87)+v1120*e(91)+v1130*e(95)+v1210*e(103)
     +       +v1220*e(107)+v1230*e(111)+v1310*e(119)+v1320*e(123)
     +       +v1330*e(127)+v2110*e(151)+v2120*e(155)+v2130*e(159)
     +       +v2210*e(167)+v2220*e(171)+v2230*e(175)+v2310*e(183)
     +       +v2320*e(187)+v2330*e(191)+v3110*e(215)+v3120*e(219)
     +       +v3130*e(223)+v3210*e(231)+v3220*e(235)+v3230*e(239)
     +       +v3310*e(247)+v3320*e(251)+v3330*e(255))*cpppp
      ve34 = ve34 + (v1110*e(88)+v1120*e(92)+v1130*e(96)+v1210*e(104)
     +       +v1220*e(108)+v1230*e(112)+v1310*e(120)+v1320*e(124)
     +       +v1330*e(128)+v2110*e(152)+v2120*e(156)+v2130*e(160)
     +       +v2210*e(168)+v2220*e(172)+v2230*e(176)+v2310*e(184)
     +       +v2320*e(188)+v2330*e(192)+v3110*e(216)+v3120*e(220)
     +       +v3130*e(224)+v3210*e(232)+v3220*e(236)+v3230*e(240)
     +       +v3310*e(248)+v3320*e(252)+v3330*e(256))*cpppp
      ve13 = ve13 + (v1101*e(86)+v1102*e(87)+v1103*e(88)+v1201*e(102)
     +       +v1202*e(103)+v1203*e(104)+v1301*e(118)+v1302*e(119)
     +       +v1303*e(120)+v2101*e(150)+v2102*e(151)+v2103*e(152)
     +       +v2201*e(166)+v2202*e(167)+v2203*e(168)+v2301*e(182)
     +       +v2302*e(183)+v2303*e(184)+v3101*e(214)+v3102*e(215)
     +       +v3103*e(216)+v3201*e(230)+v3202*e(231)+v3203*e(232)
     +       +v3301*e(246)+v3302*e(247)+v3303*e(248))*cpppp
      ve23 = ve23 + (v1101*e(90)+v1102*e(91)+v1103*e(92)+v1201*e(106)
     +       +v1202*e(107)+v1203*e(108)+v1301*e(122)+v1302*e(123)
     +       +v1303*e(124)+v2101*e(154)+v2102*e(155)+v2103*e(156)
     +       +v2201*e(170)+v2202*e(171)+v2203*e(172)+v2301*e(186)
     +       +v2302*e(187)+v2303*e(188)+v3101*e(218)+v3102*e(219)
     +       +v3103*e(220)+v3201*e(234)+v3202*e(235)+v3203*e(236)
     +       +v3301*e(250)+v3302*e(251)+v3303*e(252))*cpppp
      ve33 = ve33 + (v1101*e(94)+v1102*e(95)+v1103*e(96)+v1201*e(110)
     +       +v1202*e(111)+v1203*e(112)+v1301*e(126)+v1302*e(127)
     +       +v1303*e(128)+v2101*e(158)+v2102*e(159)+v2103*e(160)
     +       +v2201*e(174)+v2202*e(175)+v2203*e(176)+v2301*e(190)
     +       +v2302*e(191)+v2303*e(192)+v3101*e(222)+v3102*e(223)
     +       +v3103*e(224)+v3201*e(238)+v3202*e(239)+v3203*e(240)
     +       +v3301*e(254)+v3302*e(255)+v3303*e(256))*cpppp
      ve12 = ve12 + (v1011*e(86)+v1012*e(87)+v1013*e(88)+v1021*e(90)
     +       +v1022*e(91)+v1023*e(92)+v1031*e(94)+v1032*e(95)
     +       +v1033*e(96)+v2011*e(150)+v2012*e(151)+v2013*e(152)
     +       +v2021*e(154)+v2022*e(155)+v2023*e(156)+v2031*e(158)
     +       +v2032*e(159)+v2033*e(160)+v3011*e(214)+v3012*e(215)
     +       +v3013*e(216)+v3021*e(218)+v3022*e(219)+v3023*e(220)
     +       +v3031*e(222)+v3032*e(223)+v3033*e(224))*cpppp
      ve22 = ve22 + (v1011*e(102)+v1012*e(103)+v1013*e(104)+v1021*e(106)
     +       +v1022*e(107)+v1023*e(108)+v1031*e(110)+v1032*e(111)
     +       +v1033*e(112)+v2011*e(166)+v2012*e(167)+v2013*e(168)
     +       +v2021*e(170)+v2022*e(171)+v2023*e(172)+v2031*e(174)
     +       +v2032*e(175)+v2033*e(176)+v3011*e(230)+v3012*e(231)
     +       +v3013*e(232)+v3021*e(234)+v3022*e(235)+v3023*e(236)
     +       +v3031*e(238)+v3032*e(239)+v3033*e(240))*cpppp
      ve32 = ve32 + (v1011*e(118)+v1012*e(119)+v1013*e(120)+v1021*e(122)
     +       +v1022*e(123)+v1023*e(124)+v1031*e(126)+v1032*e(127)
     +       +v1033*e(128)+v2011*e(182)+v2012*e(183)+v2013*e(184)
     +       +v2021*e(186)+v2022*e(187)+v2023*e(188)+v2031*e(190)
     +       +v2032*e(191)+v2033*e(192)+v3011*e(246)+v3012*e(247)
     +       +v3013*e(248)+v3021*e(250)+v3022*e(251)+v3023*e(252)
     +       +v3031*e(254)+v3032*e(255)+v3033*e(256))*cpppp
      ve11 = ve11 + (v0111*e(86)+v0112*e(87)+v0113*e(88)+v0121*e(90)
     +       +v0122*e(91)+v0123*e(92)+v0131*e(94)+v0132*e(95)
     +       +v0133*e(96)+v0211*e(102)+v0212*e(103)+v0213*e(104)
     +       +v0221*e(106)+v0222*e(107)+v0223*e(108)+v0231*e(110)
     +       +v0232*e(111)+v0233*e(112)+v0311*e(118)+v0312*e(119)
     +       +v0313*e(120)+v0321*e(122)+v0322*e(123)+v0323*e(124)
     +       +v0331*e(126)+v0332*e(127)+v0333*e(128))*cpppp
      ve21 = ve21 + (v0111*e(150)+v0112*e(151)+v0113*e(152)+v0121*e(154)
     +       +v0122*e(155)+v0123*e(156)+v0131*e(158)+v0132*e(159)
     +       +v0133*e(160)+v0211*e(166)+v0212*e(167)+v0213*e(168)
     +       +v0221*e(170)+v0222*e(171)+v0223*e(172)+v0231*e(174)
     +       +v0232*e(175)+v0233*e(176)+v0311*e(182)+v0312*e(183)
     +       +v0313*e(184)+v0321*e(186)+v0322*e(187)+v0323*e(188)
     +       +v0331*e(190)+v0332*e(191)+v0333*e(192))*cpppp
      ve31 = ve31 + (v0111*e(214)+v0112*e(215)+v0113*e(216)+v0121*e(218)
     +       +v0122*e(219)+v0123*e(220)+v0131*e(222)+v0132*e(223)
     +       +v0133*e(224)+v0211*e(230)+v0212*e(231)+v0213*e(232)
     +       +v0221*e(234)+v0222*e(235)+v0223*e(236)+v0231*e(238)
     +       +v0232*e(239)+v0233*e(240)+v0311*e(246)+v0312*e(247)
     +       +v0313*e(248)+v0321*e(250)+v0322*e(251)+v0323*e(252)
     +       +v0331*e(254)+v0332*e(255)+v0333*e(256))*cpppp
      return
      end
**==jkdr80.f
      subroutine jkdr80(zscftp,core,prefac,iso,nshels)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      parameter (mxp2 = mxprms*mxprms)
c
INCLUDE(common/specal)
INCLUDE(common/cslosc)
INCLUDE(common/iofile)
INCLUDE(common/nshel)
INCLUDE(common/infoa)
INCLUDE(common/restar)
INCLUDE(common/timez)
INCLUDE(common/grad2)
INCLUDE(common/funct)
INCLUDE(common/symtry)
INCLUDE(common/dmisc)
INCLUDE(common/ijlab)
c
      character*10 charwall
      dimension zscft(4),f(3,maxat),iso(nshels,*)
      dimension core(*)
c
      equivalence (f(1,1),egrad(1))
      data zscft/'rhf','uhf','grhf','gvb'/
c
      if (odscf) then
         do 20 kqkq = 1 , 10
            intcut(kqkq) = 0
 20      continue
      end if
c     tim0 = cpulft(1)
c
c     ---- determine scf type ----
c
      iscf = locatc(zscft,4,zscftp)
      if (iscf.eq.0) call caserr(
     +    'invalid scftype requested in derivatives')
      ntt = (num*(num+1))/2
      nbsq = num*num
      nat3 = 3*nat
c
      ntpdm = 1
      lenb = lensec(nx)
      ndens = 1
      ofock = .false.
      ompir = .false.
      ofokab = .false.
c
c     ---- set up pointers for core ----
c
      length = 9*mxp2 + 5*mxshel + num + ntt 
      if (iscf.eq.2) length = length + ntt
      if (iscf.ge.3) length = length + nbsq

      i10 = igmem_alloc(length)
      i20 = i10 + num
      i30b = i20 + ntt
      i30c = i30b + mxp2
      i30d = i30c + mxp2
      i30e = i30d + mxp2
      i30f = i30e + mxp2
      i30g = i30f + mxp2
      i30h = i30g + mxp2
      i30i = i30h + mxp2
      i30j = i30i + mxp2
      i30k = i30j + mxp2
      i30l = i30k + mxshel
      i30m = i30l + mxshel
      i30n = i30m + mxshel
      i30o = i30n + mxshel
      i30 = i30o + mxshel
c     i40 = i30
c     if (iscf.eq.2) i40 = i30 + ntt
c     if (iscf.ge.3) i40 = i30 + nbsq

      call vclr(core(i10),1,length)

      do 30 ig = 1 , nshels
         core(i30l-1+ig) = c(1,katom(ig))
         core(i30m-1+ig) = c(2,katom(ig))
         core(i30n-1+ig) = c(3,katom(ig))
 30   continue
c
      call ddebut(zscftp,core(i20),core(i30),core(i10),core(i10),
     +            core(i10),num,ntt,nbsq,ist,jst,kst,lst,core)
c
      if (ist.le.nshels) then
c
c     calculate 2-electron contribution to hartree-fock forces.
c     routine returns with vee, and force contributions in
c     fxyz.
c
         call twldrv(iscf,core(i20),core(i30),vee,de,
     +               nprint,iso,nshels,core(i10),core(i30b),core(i30c)
     +               ,core(i30d),core(i30e),core(i30f),core(i30g),
     +               core(i30h),core(i30i),core(i30j),core(i30k),
     +               core(i30l),core(i30m),core(i30n),core(i30o),
     +               prefac,core(1)
     +               )
         if (odscf .and. nprint.ne.-5) then
            write (iwr,6010) intcut(1) , intcut(2)
         end if
         call texit(0,6)
         if (tim.ge.timlim) then
            write (iwr,6030) tim , charwall(), ist , jst , kst , lst
            write (iwr,6040)
            go to 40
         end if
      end if
      if (onocnt) then
         if (nt.ne.1) then
c
c     ----- allocate core memory for symde
c
          isymd = igmem_alloc(nw196(6))
c
          call symde(core(isymd),nat)
c
c     ----- reset core memory from symde
c
          call gmem_free(isymd)
c
      endif
      end if
      call dfinal(core,1,nshels)
 40   continue

      call gmem_free(i10)

      return
 6010 format (/' integral cuts'/1x,13('-')/' on ij   shells ',
     +        i8/' on ijkl shells ',i8)
 6020 format (/)
c9028 format(///'  end of calculation of the energy',
c    *' gradient at ',f8.2,' seconds'/)
 6030 format (//
     +        ' insufficient time to complete evaluation of 2-electron',
     +        ' contribution to gradient'///' job dumped at ',f8.2,
     +        ' seconds',a10,' wall'//' next batch ',4i5)
 6040 format (///,10x,27('*')//10x,'*** warning ***'/10x,
     +        'this job must be restarted'//10x,27('*')//)
      end
**==twldrv.f
_IF(hpux11)
c HP compiler bug JAGae56094
c$HP$ OPTIMIZE LEVEL2
_ENDIF
      subroutine twldrv(iscf,dm,dn,vee,fxyz,
     $                  lprint,iso,nshels,nconf,
     *   qx,qy,qz, pe34,ts3,ts4,ts34,ta34i,csmcd,iatm,p,q,r,ishelt,
     *   prefac,core)
c
c        evaluation of the two electron integral contribution to the
c        forces
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
INCLUDE(common/sizes)
      parameter (mxp2 = mxprms*mxprms)
c
INCLUDE(common/cslosc)
c
INCLUDE(common/restar)
      common/restri/jjfile(63),lds(508),isect(508),ldsect(508)
INCLUDE(common/timez)
INCLUDE(common/iofile)
INCLUDE(common/ffq)
INCLUDE(common/auxvar)
INCLUDE(common/prints)
INCLUDE(common/nshel)
INCLUDE(common/dindex)
      common /junk/ tbaa(400),tbba(400),tbca(400),
     *               tbab(400),tbbb(400),tbcb(400),
     *               tbac(400),tbbc(400),tbcc(400),
     *               tbad(400),tbbd(400),tbcd(400),
     *               tbae(400),tbbe(400),tbce(400),
     *               tbaf(400),tbbf(400),tbcf(400)
INCLUDE(common/symtry)
INCLUDE(common/parallel)
c
      dimension prefac(*)
      dimension dm(*),nconf(*),dn(*),fxyz(*),iso(nshels,*)
      dimension p(*),q(*),r(*),ishelt(*),core(*)
      dimension iatm(*)
      dimension pe34(*),qx(*),qy(*),qz(*),
     $          ta34i(*),ts3(*),ts4(*),ts34(*),csmcd(*)
      dimension m0(48),m1(48),m2(48),m3(48)
c
c
c
      data three/3.0d0/
      data p25/0.25d0/
      data two,h,dzero/2.0d0,0.5d0,0.0d0/
      data done5/1.5d0/
      data twenty/20.0d0/
      data two5,three5,four5/2.5d0,3.5d0,4.5d0/
      data sxtn/16.0d0/
      data done,four/1.0d0,4.0d0/
      data tenm7,tenm4,tenm2/1.0d-7,1.0d-4,1.0d-2/
      data tenm8,tenm12,tenm20/1.0d-8,1.0d-12,1.0d-20/
      data rln10 /2.30258d0/
      data m22,mword/22,7200/
c
      if (nprint.ne.-5 .and. oprint(57)) then
         write (iwr,6020)
         call writel(prefac,nshels)
      end if
c
c     initialize some variables
      css = dzero
      csc = dzero
      cpc = dzero
      cps = dzero
      csp = dzero
      cpp = dzero
      csd = dzero
      cpd = dzero
      ir = ist
c
c      set up f(m,t) interpolation table in common table
c
      call secget(isect(503),m22,iblk22)
      call rdedx(tbaa,mword,iblk22,idaf)
c
c     compute pi related constants.
c
      pi = four*datan(done)
      pito52 = two*(pi**two5)
      pidiv4 = p25*pi
c     trp = two/dsqrt(pi)
c
c     set up correspondence table between shells and atoms.
c
      do 20 i = 1 , nshels
         ishelt(i) = ktype(i) - 1
         iatm(i) = katom(i)
 20   continue
c
c     clear out accumulator vee.
c
      vee = dzero
c
c     set up critical constants.
c
      tol = rln10 * itol
c
c     set accuracy factors (see auxvar)
c
      if (ifastd.eq.2) then
       vtol  = tenm8
       vtol1 = tenm12
       vtol2 = tenm12
       vtols = tenm20
      else if (ifastd.eq.1) then
       vtol  = tenm8
       vtol1 = vtol*tenm4
       vtol2 = vtol*tenm4
       vtols = vtol**2
      else
       vtol = tenm7
       vtol1 = vtol*tenm4
       vtol2 = vtol*tenm2
       vtols = vtol**2
      endif
c
_IF(parallel)
c***   **MPP**
      next = ipg_dlbtask()
c***   **MPP**
_ENDIF
c        begin loops over the shells
c     ---- ishell ----
c
_IF(parallel)
c      do 310 ishell = nshels , ir , -1
      do 310 ishell = ir , nshels
_ELSE
      do 310 ishell = ir , nshels
_ENDIF
         if (ishelt(ishell).le.1) then
            ijshel = ishell*(ishell-1)/2
            do 40 it = 1 , nt
               ids = iso(ishell,it)
               if (ids.gt.ishell) go to 310
               m0(it) = ids
 40         continue
c
c     ---- jshell ----
c
_IF(parallel)
c            do 300 jshell = ishell , 1 , -1
            do 300 jshell = 1 , ishell
_ELSE
            do 300 jshell = 1 , ishell
_ENDIF
               if (ishelt(jshell).le.1) then
                  tolij = dlntol + prefac(ijshel+jshell)
c *** selection on ij parts must not be too heavy. 3.401=ln(30)
                  if (tolij.gt.-3.401d0) then
                     do 60 it = 1 , nt
                        ids = m0(it)
                        jds = iso(jshell,it)
                        if (jds.gt.ishell) go to 300
                        if (ids.lt.jds) then
                           nds = ids
                           ids = jds
                           jds = nds
                        end if
                        if (ids.eq.ishell .and. jds.gt.jshell) go to 300
                        m1(it) = ids
                        m2(it) = jds
 60                  continue
_IF(parallel)
c***   **MPP**
                icount_dlb = icount_dlb + 1
                if(icount_dlb .eq. next) then
c***   **MPP**
_ENDIF
c
c     ---- kshell ----
c
                     do 290 kshell = 1 , ishell
                        if (ishelt(kshell).le.1) then
                           klshel = kshell*(kshell-1)/2
                           do 80 it = 1 , nt
                              kds = iso(kshell,it)
                              if (kds.gt.ishell) go to 290
                              m3(it) = kds
 80                        continue
c
c     ---- lshell ----
c
                           do 280 lshell = 1 , kshell
                              if (ishelt(lshell).le.1) then
                                 tolijk = tolij + prefac(klshel+lshell)
                                 if (tolijk.gt.0.d0) then
                                    if (kshell.eq.ishell .and.
     +                                  lshell.gt.jshell) go to 3001
                                    n4 = 0
                                    do 100 it = 1 , nt
                                       mds = iso(lshell,it)
                                       if (mds.gt.ishell) go to 280
                                       kds = m3(it)
                                       if (kds.lt.mds) then
                                         nds = kds
                                         kds = mds
                                         mds = nds
                                       end if
                                       ids = m1(it)
                                       jds = m2(it)
                                       if (ids.eq.ishell .or.
     +                                    kds.eq.ishell) then
                                         if (kds.ge.ids) then
                                         if (kds.ne.ids .or. mds.gt.jds)
     +                                      then
                                         nds = ids
                                         ids = kds
                                         kds = nds
                                         nds = jds
                                         jds = mds
                                         mds = nds
                                         end if
                                         end if
                                         if (jds.ge.jshell) then
                                         if (jds.gt.jshell) go to 280
                                         if (kds.ge.kshell) then
                                         if (kds.gt.kshell) go to 280
                                         if (mds.ge.lshell) then
                                         if (mds.gt.lshell) go to 280
                                         n4 = n4 + 1
                                         end if
                                         end if
                                         end if
                                       end if
 100                                continue
c
c     ---- calculate q4 for these shells ----
c
                                    q4 = dfloat(nt)/dfloat(n4)
c
                                    ax1 = h
                                    if (ishell.ne.jshell) ax1 = ax1 +
     +                                  ax1
                                    if (kshell.ne.lshell) ax1 = ax1 +
     +                                  ax1
                                    if (ishell.ne.kshell .or.
     +                                  jshell.ne.lshell) ax1 = ax1 +
     +                                  ax1
c
c        put shells into standard order
c
                                    inew = ishell
                                    jnew = jshell
                                    knew = kshell
                                    lnew = lshell
                                    la = ishelt(inew)
                                    lb = ishelt(jnew)
                                    lc = ishelt(knew)
                                    ld = ishelt(lnew)
                                    if (la.lt.lb) then
                                       inew = jshell
                                       jnew = ishell
                                    end if
                                    if (lc.lt.ld) then
                                       knew = lshell
                                       lnew = kshell
                                    end if
                                    if (la+lb-lc.lt.ld) then
                                       id = inew
                                       inew = knew
                                       knew = id
                                       id = jnew
                                       jnew = lnew
                                       lnew = id
                                    end if
                                    la = 3*ishelt(inew) + 1
                                    lb = 3*ishelt(jnew) + 1
                                    lc = 3*ishelt(knew) + 1
                                    ld = 3*ishelt(lnew) + 1
c
c        obtain information about the shells
c
                                    oijsam = inew.eq.jnew
                                    oklsam = knew.eq.lnew
                                    oijkl = (inew.eq.knew) .and.
     +                                 (jnew.eq.lnew)
                                    jtype = (la+lb+lc+lc+ld-2)/3
c
                                    if (jtype.gt.6) call caserr(
     +                'invalid orbital type in derivative routines')
c
                                    ifqmax = ishelt(inew) + ishelt(jnew)
     +                                 + ishelt(knew) + ishelt(lnew) + 1
                                    iat = iatm(inew)
                                    jat = iatm(jnew)
                                    kat = iatm(knew)
                                    lat = iatm(lnew)
                                    if ((iat.ne.jat) .or. (iat.ne.kat)
     +                                  .or. (iat.ne.lat)) then
                                       iatx = 3*(iat-1) + 1
                                       jatx = 3*(jat-1) + 1
                                       katx = 3*(kat-1) + 1
                                       latx = 3*(lat-1) + 1
                                       iaty = iatx + 1
                                       jaty = jatx + 1
                                       katy = katx + 1
                                       laty = latx + 1
                                       iatz = iaty + 1
                                       jatz = jaty + 1
                                       katz = katy + 1
                                       latz = laty + 1
                                       ax = p(inew)
                                       bx = p(jnew)
                                       cx = p(knew)
                                       dx = p(lnew)
                                       ay = q(inew)
                                       by = q(jnew)
                                       cy = q(knew)
                                       dy = q(lnew)
                                       az = r(inew)
                                       bz = r(jnew)
                                       cz = r(knew)
                                       dz = r(lnew)
                                       abx = ax - bx
                                       aby = ay - by
                                       abz = az - bz
                                       cdx = cx - dx
                                       cdy = cy - dy
                                       cdz = cz - dz
                                       r34 = cdx**2 + cdy**2 + cdz**2
                                       r12 = abx**2 + aby**2 + abz**2
                                       istart = kstart(inew)
                                       jstart = kstart(jnew)
                                       kstars = kstart(knew)
                                       lstart = kstart(lnew)
                                       iend = istart + kng(inew) - 1
                                       jend = jstart + kng(jnew) - 1
                                       kend = kstars + kng(knew) - 1
                                       lend = lstart + kng(lnew) - 1
c
c     fill e with the required density matrix contributions,
c     and at the same time, get dmax.
c
                                      if(iscf.eq.1) then
                                       call efill1(inew,jnew,knew,lnew,
     +                                    la,lb,lc,ld,ax1,dm,e,dmax,q4)
                                      else if (iscf.eq.2) then
                                       call efill2(inew,jnew,knew,lnew,
     +                                    la,lb,lc,ld,ax1,dm,dn,e,
     +                                    dmax,q4)
                                      else
                                       call efill3(inew,jnew,knew,lnew,
     +                                    la,lb,lc,ld,ax1,dm,dn,e,
     +                                    dmax,q4,nconf)
                                      endif
c
c     reject current shell case if dmax is .lt. vtol.
c
                                       if (dmax.ge.vtol1) then
c
c     hunt out largest density*contraction, and reject if possible.
c
                                         e34max = dzero
                                         kdzero = -mxprms
                                         do 120 k = kstars , kend
                                         kdzero = kdzero + mxprms
                                         kl = kdzero
                                         a3 = ex(k)
                                         csumc=dabs(cs(k))+dabs(cp(k))
                                         lnd = lend
                                         if (oklsam) lnd = k
                                         do 110 l = lstart , lnd
                                         a4 = ex(l)
                                         kl = kl + 1
                                         a34 = a3 + a4
                                         a34i = done/a34
                                         s3 = a3*a34i
                                         s4 = a4*a34i
                                         s34 = a3*s4
                                         qx(kl) = s3*cx + s4*dx
                                         qy(kl) = s3*cy + s4*dy
                                         qz(kl) = s3*cz + s4*dz
                                         ta34i(kl) = a34i
                                         ts3(kl) = s3
                                         ts4(kl) = s4
                                         ts34(kl) = s34
                                         dum = r34*s34
                                         if(dum.le.tol) then
                                            e34 = dexp(-dum)*a34i
                                            if (oklsam .and. (k.ne.l))
     +                                         e34 = e34 + e34
                                            pe34(kl) = e34
                                            e34 = e34*csumc*(dabs(cs(l))
     +                                         +dabs(cp(l)))
                                         else
                                            e34 = 0.0d0
                                            pe34(kl) = 0.0d0
                                         endif
                                         if (e34.gt.e34max) e34max = e34
                                         csmcd(kl) = e34**2
 110                                     continue
 120                                     continue
                                         if (dmax*e34max.ge.vtol1) then
                                         fxi = dzero
                                         fxj = dzero
                                         fxk = dzero
                                         fxl = dzero
                                         fyi = dzero
                                         fyj = dzero
                                         fyk = dzero
                                         fyl = dzero
                                         fzi = dzero
                                         fzj = dzero
                                         fzk = dzero
                                         fzl = dzero
                                         ve11 = dzero
                                         ve12 = dzero
                                         ve13 = dzero
                                         ve14 = dzero
                                         ve21 = dzero
                                         ve22 = dzero
                                         ve23 = dzero
                                         ve24 = dzero
                                         ve31 = dzero
                                         ve32 = dzero
                                         ve33 = dzero
                                         ve34 = dzero
c
c        loop over the uncontracted gaussians within the shells
c
                                         idzero = -mxprms
c
c..... gaussians at center a.
c
                                         do 270 i = istart , iend
                                         idzero = idzero + mxprms
                                         ij = idzero
                                         a1 = ex(i)
                                         csa = cs(i)
                                         cpa = cp(i)
                                         csuma = (dabs(csa)+dabs(cpa))
     +                                      *dmax
                                         jnd = jend
                                         if (oijsam) jnd = i
c
c..... gaussians at center b.
c
                                         do 260 j = jstart , jnd
                                         ij = ij + 1
                                         a2 = ex(j)
                                         a12 = a1 + a2
                                         a12i = done/a12
                                         s1 = a1*a12i
                                         s2 = a2*a12i
                                         s12 = a1*s2
                                         dum = r12*s12
                                         if (dum .le. tol) then
                                            e12=dexp(-dum)*pito52*a12i
                                            if (oijsam .and. (i.ne.j))
     +                                      e12 = e12 + e12
                                            csb = cs(j)*e12
                                            cpb = cp(j)*e12
                                            csmab = csuma*(dabs(csb)
     +                                           +dabs(cpb))
                                         else
                                            e12 = 0.0d0
                                            csb = 0.0d0
                                            cpb = 0.0d0
                                            csmab = 0.0d0
                                         endif
                                         if (csmab*e34max.ge.vtol2) then
                                         csmab = csmab**2
                                         css = csa*csb
                                         px = s1*ax + s2*bx
                                         py = s1*ay + s2*by
                                         pz = s1*az + s2*bz
                                         if (la.ne.1) then
                                         cps = cpa*csb
                                         ppxo = -s2*abx
                                         ppyo = -s2*aby
                                         ppzo = -s2*abz
                                         if (lb.ne.1) then
                                         csp = csa*cpb
                                         cpp = cpa*cpb
                                         ppox = s1*abx
                                         ppoy = s1*aby
                                         ppoz = s1*abz
                                         ppxox = ppxo + ppox
                                         ppyoy = ppyo + ppoy
                                         ppzoz = ppzo + ppoz
                                         ha12i = h*a12i
                                         ppxx = ppxo*ppox + ha12i
                                         ppyy = ppyo*ppoy + ha12i
                                         ppzz = ppzo*ppoz + ha12i
                                         ppxy = ppxo*ppoy
                                         ppyx = ppxy
                                         ppxz = ppxo*ppoz
                                         ppzx = ppxz
                                         ppyz = ppyo*ppoz
                                         ppzy = ppyz
                                         end if
                                         end if
                                         ve00s = dzero
                                         ve11s = dzero
                                         ve12s = dzero
                                         ve21s = dzero
                                         ve22s = dzero
                                         ve31s = dzero
                                         ve32s = dzero
                                         dvexs = dzero
                                         dveys = dzero
                                         dvezs = dzero
                                         knd = kend
                                         if (oijkl) knd = i
                                         kdzero = -mxprms
c
c..... gaussians at center c.
c
                                         do 250 k = kstars , knd
                                         kdzero = kdzero + mxprms
                                         kl = kdzero
                                         a3 = ex(k)
                                         csc = cs(k)
                                         cpc = cp(k)
                                         csss = css*csc
                                         cssp = css*cpc
                                         cpss = cps*csc
                                         cpsp = cps*cpc
                                         csps = csp*csc
                                         cspp = csp*cpc
                                         cpps = cpp*csc
                                         cppp = cpp*cpc
                                         lnd = lend
                                         if (oklsam) lnd = k
                                         if (oijkl .and. (i.eq.k))
     +                                      lnd = j
c
c..... gaussians at center d.
c
                                         do 240 l = lstart , lnd
                                         kl = kl + 1
                                         a34 = a3 + ex(l)
                                         a1234 = a12 + a34
                                         a1234i = done/a1234
                                         vtest = csmab*csmcd(kl)*a1234i
                                         if (vtest.lt.vtols) go to 240
                                         qa2 = a34*a1234i
                                         qa = a12*qa2
                                         pqx = px - qx(kl)
                                         pqy = py - qy(kl)
                                         pqz = pz - qz(kl)
                                         pqxx = pqx*pqx
                                         pqyy = pqy*pqy
                                         pqzz = pqz*pqz
                                         t = qa*(pqxx+pqyy+pqzz)
                                         a34i = ta34i(kl)
                                         s3 = ts3(kl)
                                         s4 = ts4(kl)
                                         s34 = ts34(kl)
                                         e34 = pe34(kl)*dsqrt(a1234i)
                                         if (oijkl .and. (ij.ne.kl))
     +                                      e34 = e34 + e34
                                         csd = cs(l)*e34
                                         eoooo = e(1)*csss*csd
                                         if (t.lt.sxtn) then
                                         qq = t*twenty
                                         theta = qq - idint(qq)
                                         n = qq - theta
                                         theta2 = theta*(theta-done)
                                         theta3 = theta2*(theta-two)
                                         theta4 = theta2*(theta+done)
                                         go to (170,160,150,140,130) ,
     +                                      ifqmax
                                         else
                                         tj = done/t
                                         if (vtest*tj.lt.vtols)
     +                                      go to 240
                                         fq0 = dsqrt(pidiv4*tj)
                                         fq1 = h*fq0*tj
                                         if (jtype.eq.1) go to 210
                                         fq2 = done5*fq1*tj
                                         if (jtype.eq.2) go to 220
                                         fq3 = two5*fq2*tj
                                         fq4 = three5*fq3*tj
                                         fq5 = four5*fq4*tj
                                         go to 180
                                         end if
 130                                     fq5 = tbaf(n+1)
     +                                      + theta*tbbf(n+1)
     +                                      - theta3*tbcf(n+1)
     +                                      + theta4*tbcf(n+2)
 140                                     fq4 = tbae(n+1)
     +                                      + theta*tbbe(n+1)
     +                                      - theta3*tbce(n+1)
     +                                      + theta4*tbce(n+2)
 150                                     fq3 = tbad(n+1)
     +                                      + theta*tbbd(n+1)
     +                                      - theta3*tbcd(n+1)
     +                                      + theta4*tbcd(n+2)
 160                                     fq2 = tbac(n+1)
     +                                      + theta*tbbc(n+1)
     +                                      - theta3*tbcc(n+1)
     +                                      + theta4*tbcc(n+2)
 170                                     fq1 = tbab(n+1)
     +                                      + theta*tbbb(n+1)
     +                                      - theta3*tbcb(n+1)
     +                                      + theta4*tbcb(n+2)
                                         fq0 = tbaa(n+1)
     +                                      + theta*tbba(n+1)
     +                                      - theta3*tbca(n+1)
     +                                      + theta4*tbca(n+2)
 180                                     if (jtype.eq.1) go to 210
                                         if (jtype.eq.2) go to 220
                                         pqxy = pqx*pqy
                                         pqyz = pqy*pqz
                                         pqxz = pqx*pqz
                                         qa1 = a12*a1234i
                                         cpd = cp(l)*e34
                                         if (lc.ne.1) then
                                         cssps = cssp*csd
                                         cpsps = cpsp*csd
                                         qqxo = -s4*cdx
                                         qqyo = -s4*cdy
                                         qqzo = -s4*cdz
                                         if (ld.ne.1) then
                                         csssp = csss*cpd
                                         csspp = cssp*cpd
                                         cspsp = csps*cpd
                                         cpssp = cpss*cpd
                                         csppp = cspp*cpd
                                         cpspp = cpsp*cpd
                                         cppsp = cpps*cpd
                                         cpppp = cppp*cpd
                                         qqox = s3*cdx
                                         qqoy = s3*cdy
                                         qqoz = s3*cdz
                                         qqxx = qqxo*qqox + h*a34i
                                         qqyy = qqyo*qqoy + h*a34i
                                         qqzz = qqzo*qqoz + h*a34i
                                         qqxy = qqxo*qqoy
                                         qqxz = qqxo*qqoz
                                         qqyx = qqxy
                                         qqyz = qqyo*qqoz
                                         qqzx = qqxz
                                         qqzy = qqyz
                                         qqxox = qqxo + qqox
                                         qqyoy = qqyo + qqoy
                                         qqzoz = qqzo + qqoz
                                         end if
                                         end if
c
c   ******************************************************************
c
c        begin the first pass through the integral calculations
c
c   ********************************************************************
c
                                         oflag = .false.
c
c     integral evaluation section
c
 190                                     goooo = fq0
                                         v0000 = goooo
                                         ve00 = v0000*eoooo
                                         if (jtype.lt.2) go to 200
c
c        (ps,ss) section
c
                                         qfq1 = -qa2*fq1
                                         gxooo = pqx*qfq1
                                         v1000 = ppxo*goooo + gxooo
                                         gyooo = pqy*qfq1
                                         v2000 = ppyo*goooo + gyooo
                                         gzooo = pqz*qfq1
                                         v3000 = ppzo*goooo + gzooo
                                         cpsss = cpss*csd
                                         ve00 = ve00 +
     +                                      (v1000*e(65)+v2000*e(129)
     +                                      +v3000*e(193))*cpsss
                                         temp = v0000*cpsss
                                         ve11 = temp*e(65)
                                         ve21 = temp*e(129)
                                         ve31 = temp*e(193)
                                         if (jtype.lt.3) go to 200
                                         if (jtype.ne.3) then
c
c        (ps,ps) + (ss,ps) section
c
                                         qfq1 = qa1*fq1
                                         gooxo = pqx*qfq1
                                         v0010 = qqxo*goooo + gooxo
                                         gooyo = pqy*qfq1
                                         v0020 = qqyo*goooo + gooyo
                                         goozo = pqz*qfq1
                                         v0030 = qqzo*goooo + goozo
                                         hfq1 = h*fq1*a1234i
                                         qfq2 = -qa*fq2*a1234i
                                         gxoxo = pqxx*qfq2 + hfq1
                                         v1010 = ppxo*v0010 +
     +                                      qqxo*gxooo + gxoxo
                                         gyoyo = pqyy*qfq2 + hfq1
                                         v2020 = ppyo*v0020 +
     +                                      qqyo*gyooo + gyoyo
                                         gzozo = pqzz*qfq2 + hfq1
                                         v3030 = ppzo*v0030 +
     +                                      qqzo*gzooo + gzozo
                                         gxoyo = pqxy*qfq2
                                         v1020 = ppxo*v0020 +
     +                                      qqyo*gxooo + gxoyo
                                         v2010 = ppyo*v0010 +
     +                                      qqxo*gyooo + gxoyo
                                         gxozo = pqxz*qfq2
                                         v1030 = ppxo*v0030 +
     +                                      qqzo*gxooo + gxozo
                                         v3010 = ppzo*v0010 +
     +                                      qqxo*gzooo + gxozo
                                         gyozo = pqyz*qfq2
                                         v2030 = ppyo*v0030 +
     +                                      qqzo*gyooo + gyozo
                                         v3020 = ppzo*v0020 +
     +                                      qqyo*gzooo + gyozo
                                         ve00 = ve00 +
     +                                      (v0010*e(5)+v0020*e(9)
     +                                      +v0030*e(13))*cssps
                                         temp = v0000*cssps
                                         ve13 = temp*e(5)
                                         ve23 = temp*e(9)
                                         ve33 = temp*e(13)
                                         ve00 = ve00 +
     +                                      (v1010*e(69)+v1020*e(73)
     +                                      +v1030*e(77)+v2010*e(133)
     +                                      +v2020*e(137)+v2030*e(141)
     +                                      +v3010*e(197)+v3020*e(201)
     +                                      +v3030*e(205))*cpsps
                                         ve11 = ve11 +
     +                                      (v0010*e(69)+v0020*e(73)
     +                                      +v0030*e(77))*cpsps
                                         ve13 = ve13 +
     +                                      (v1000*e(69)+v2000*e(133)
     +                                      +v3000*e(197))*cpsps
                                         ve21 = ve21 +
     +                                      (v0010*e(133)+v0020*e(137)
     +                                      +v0030*e(141))*cpsps
                                         ve23 = ve23 +
     +                                      (v1000*e(73)+v2000*e(137)
     +                                      +v3000*e(201))*cpsps
                                         ve31 = ve31 +
     +                                      (v0010*e(197)+v0020*e(201)
     +                                      +v0030*e(205))*cpsps
                                         ve33 = ve33 +
     +                                      (v1000*e(77)+v2000*e(141)
     +                                      +v3000*e(205))*cpsps
                                         if (jtype.eq.4) go to 200
                                         end if
c
c        (pp,ss) + (sp,ss) section
c
                                         v0100 = ppox*goooo + gxooo
                                         v0200 = ppoy*goooo + gyooo
                                         v0300 = ppoz*goooo + gzooo
                                         cspss = csps*csd
                                         ve00 = ve00 +
     +                                      (v0100*e(17)+v0200*e(33)
     +                                      +v0300*e(49))*cspss
                                         temp = v0000*cspss
                                         ve12 = temp*e(17)
                                         ve22 = temp*e(33)
                                         ve32 = temp*e(49)
                                         hfq1 = -h*qa2*fq1*a12i
                                         qfq2 = qa2*qa2*fq2
                                         gxxoo = pqxx*qfq2 + hfq1
                                         v1100 = ppxx*goooo +
     +                                      ppxox*gxooo + gxxoo
                                         gyyoo = pqyy*qfq2 + hfq1
                                         v2200 = ppyy*goooo +
     +                                      ppyoy*gyooo + gyyoo
                                         gzzoo = pqzz*qfq2 + hfq1
                                         v3300 = ppzz*goooo +
     +                                      ppzoz*gzooo + gzzoo
                                         gxyoo = pqxy*qfq2
                                         v1200 = ppxo*v0200 +
     +                                      ppoy*gxooo + gxyoo
                                         v2100 = ppyo*v0100 +
     +                                      ppox*gyooo + gxyoo
                                         gxzoo = pqxz*qfq2
                                         v1300 = ppxo*v0300 +
     +                                      ppoz*gxooo + gxzoo
                                         v3100 = ppzo*v0100 +
     +                                      ppox*gzooo + gxzoo
                                         gyzoo = pqyz*qfq2
                                         v2300 = ppyo*v0300 +
     +                                      ppoz*gyooo + gyzoo
                                         v3200 = ppzo*v0200 +
     +                                      ppoy*gzooo + gyzoo
                                         cppss = cpps*csd
                                         ve00 = ve00 +
     +                                      (v1100*e(81)+v1200*e(97)
     +                                      +v1300*e(113)+v2100*e(145)
     +                                      +v2200*e(161)+v2300*e(177)
     +                                      +v3100*e(209)+v3200*e(225)
     +                                      +v3300*e(241))*cppss
                                         ve11 = ve11 +
     +                                      (v0100*e(81)+v0200*e(97)
     +                                      +v0300*e(113))*cppss
                                         ve21 = ve21 +
     +                                      (v0100*e(145)+v0200*e(161)
     +                                      +v0300*e(177))*cppss
                                         ve31 = ve31 +
     +                                      (v0100*e(209)+v0200*e(225)
     +                                      +v0300*e(241))*cppss
                                         ve12 = ve12 +
     +                                      (v1000*e(81)+v2000*e(145)
     +                                      +v3000*e(209))*cppss
                                         ve22 = ve22 +
     +                                      (v1000*e(97)+v2000*e(161)
     +                                      +v3000*e(225))*cppss
                                         ve32 = ve32 +
     +                                      (v1000*e(113)+v2000*e(177)
     +                                      +v3000*e(241))*cppss
                                         if (jtype.ne.3) then
c
c        (pp,ps) + (sp,ps) section
c
                                         v0110 = qqxo*v0100 +
     +                                      ppox*gooxo + gxoxo
                                         v0120 = qqyo*v0100 +
     +                                      ppox*gooyo + gxoyo
                                         v0130 = qqzo*v0100 +
     +                                      ppox*goozo + gxozo
                                         v0210 = qqxo*v0200 +
     +                                      ppoy*gooxo + gxoyo
                                         v0220 = qqyo*v0200 +
     +                                      ppoy*gooyo + gyoyo
                                         v0230 = qqzo*v0200 +
     +                                      ppoy*goozo + gyozo
                                         v0310 = qqxo*v0300 +
     +                                      ppoz*gooxo + gxozo
                                         v0320 = qqyo*v0300 +
     +                                      ppoz*gooyo + gyozo
                                         v0330 = qqzo*v0300 +
     +                                      ppoz*goozo + gzozo
                                         cspps = cspp*csd
                                         ve00 = ve00 +
     +                                      (v0110*e(21)+v0120*e(25)
     +                                      +v0130*e(29)+v0210*e(37)
     +                                      +v0220*e(41)+v0230*e(45)
     +                                      +v0310*e(53)+v0320*e(57)
     +                                      +v0330*e(61))*cspps
                                         ve12 = ve12 +
     +                                      (v0010*e(21)+v0020*e(25)
     +                                      +v0030*e(29))*cspps
                                         ve22 = ve22 +
     +                                      (v0010*e(37)+v0020*e(41)
     +                                      +v0030*e(45))*cspps
                                         ve32 = ve32 +
     +                                      (v0010*e(53)+v0020*e(57)
     +                                      +v0030*e(61))*cspps
                                         ve13 = ve13 +
     +                                      (v0100*e(21)+v0200*e(37)
     +                                      +v0300*e(53))*cspps
                                         ve23 = ve23 +
     +                                      (v0100*e(25)+v0200*e(41)
     +                                      +v0300*e(57))*cspps
                                         ve33 = ve33 +
     +                                      (v0100*e(29)+v0200*e(45)
     +                                      +v0300*e(61))*cspps
                                         qfq3 = qa1*qa2*qa2*fq3
                                         hfq2 = -h*qa2*fq2*a1234i
                                         tfq2 = three*hfq2
                                         gxyzo = pqxy*pqz*qfq3
                                         c1230 = ppxy*goozo +
     +                                      ppxo*gyozo + ppoy*gxozo +
     +                                      gxyzo
                                         v1230 = qqzo*v1200 + c1230
                                         c1320 = ppxz*gooyo +
     +                                      ppxo*gyozo + ppoz*gxoyo +
     +                                      gxyzo
                                         v1320 = qqyo*v1300 + c1320
                                         c2130 = ppyx*goozo +
     +                                      ppyo*gxozo + ppox*gyozo +
     +                                      gxyzo
                                         v2130 = qqzo*v2100 + c2130
                                         c2310 = ppyz*gooxo +
     +                                      ppyo*gxozo + ppoz*gxoyo +
     +                                      gxyzo
                                         v2310 = qqxo*v2300 + c2310
                                         c3120 = ppzx*gooyo +
     +                                      ppzo*gxoyo + ppox*gyozo +
     +                                      gxyzo
                                         v3120 = qqyo*v3100 + c3120
                                         c3210 = ppzy*gooxo +
     +                                      ppzo*gxoyo + ppoy*gxozo +
     +                                      gxyzo
                                         v3210 = qqxo*v3200 + c3210
                                         temp = pqxx*qfq3
                                         gxxxo = pqx*(temp+tfq2)
                                         c1110 = ppxx*gooxo +
     +                                      ppxox*gxoxo + gxxxo
                                         v1110 = qqxo*v1100 + c1110
                                         gxxyo = pqy*(temp+hfq2)
                                         c1120 = ppxx*gooyo +
     +                                      ppxox*gxoyo + gxxyo
                                         v1120 = qqyo*v1100 + c1120
                                         c1210 = ppxy*gooxo +
     +                                      ppxo*gxoyo + ppoy*gxoxo +
     +                                      gxxyo
                                         v1210 = qqxo*v1200 + c1210
                                         c2110 = ppyx*gooxo +
     +                                      ppyo*gxoxo + ppox*gxoyo +
     +                                      gxxyo
                                         v2110 = qqxo*v2100 + c2110
                                         gxxzo = pqz*(temp+hfq2)
                                         c1130 = ppxx*goozo +
     +                                      ppxox*gxozo + gxxzo
                                         v1130 = qqzo*v1100 + c1130
                                         c1310 = ppxz*gooxo +
     +                                      ppxo*gxozo + ppoz*gxoxo +
     +                                      gxxzo
                                         v1310 = qqxo*v1300 + c1310
                                         c3110 = ppzx*gooxo +
     +                                      ppzo*gxoxo + ppox*gxozo +
     +                                      gxxzo
                                         v3110 = qqxo*v3100 + c3110
                                         temp = pqyy*qfq3
                                         gyyyo = pqy*(temp+tfq2)
                                         c2220 = ppyy*gooyo +
     +                                      ppyoy*gyoyo + gyyyo
                                         v2220 = qqyo*v2200 + c2220
                                         gyyxo = pqx*(temp+hfq2)
                                         c2210 = ppyy*gooxo +
     +                                      ppyoy*gxoyo + gyyxo
                                         v2210 = qqxo*v2200 + c2210
                                         c2120 = ppyx*gooyo +
     +                                      ppyo*gxoyo + ppox*gyoyo +
     +                                      gyyxo
                                         v2120 = qqyo*v2100 + c2120
                                         c1220 = ppxy*gooyo +
     +                                      ppxo*gyoyo + ppoy*gxoyo +
     +                                      gyyxo
                                         v1220 = qqyo*v1200 + c1220
                                         gyyzo = pqz*(temp+hfq2)
                                         c2230 = ppyy*goozo +
     +                                      ppyoy*gyozo + gyyzo
                                         v2230 = qqzo*v2200 + c2230
                                         c2320 = ppyz*gooyo +
     +                                      ppyo*gyozo + ppoz*gyoyo +
     +                                      gyyzo
                                         v2320 = qqyo*v2300 + c2320
                                         c3220 = ppzy*gooyo +
     +                                      ppzo*gyoyo + ppoy*gyozo +
     +                                      gyyzo
                                         v3220 = qqyo*v3200 + c3220
                                         temp = pqzz*qfq3
                                         gzzzo = pqz*(temp+tfq2)
                                         c3330 = ppzz*goozo +
     +                                      ppzoz*gzozo + gzzzo
                                         v3330 = qqzo*v3300 + c3330
                                         gzzxo = pqx*(temp+hfq2)
                                         c3310 = ppzz*gooxo +
     +                                      ppzoz*gxozo + gzzxo
                                         v3310 = qqxo*v3300 + c3310
                                         c3130 = ppzx*goozo +
     +                                      ppzo*gxozo + ppox*gzozo +
     +                                      gzzxo
                                         v3130 = qqzo*v3100 + c3130
                                         c1330 = ppxz*goozo +
     +                                      ppxo*gzozo + ppoz*gxozo +
     +                                      gzzxo
                                         v1330 = qqzo*v1300 + c1330
                                         gzzyo = pqy*(temp+hfq2)
                                         c3320 = ppzz*gooyo +
     +                                      ppzoz*gyozo + gzzyo
                                         v3320 = qqyo*v3300 + c3320
                                         c3230 = ppzy*goozo +
     +                                      ppzo*gyozo + ppoy*gzozo +
     +                                      gzzyo
                                         v3230 = qqzo*v3200 + c3230
                                         c2330 = ppyz*goozo +
     +                                      ppyo*gzozo + ppoz*gyozo +
     +                                      gzzyo
                                         v2330 = qqzo*v2300 + c2330
                                         cppps = cppp*csd
                                         ve00 = ve00 +
     +                                      (v1110*e(85)+v1120*e(89)
     +                                      +v1130*e(93)+v1210*e(101)
     +                                      +v1220*e(105)+v1230*e(109)
     +                                      +v1310*e(117)+v1320*e(121)
     +                                      +v1330*e(125)+v2110*e(149)
     +                                      +v2120*e(153)+v2130*e(157)
     +                                      +v2210*e(165)+v2220*e(169)
     +                                      +v2230*e(173)+v2310*e(181)
     +                                      +v2320*e(185)+v2330*e(189)
     +                                      +v3110*e(213)+v3120*e(217)
     +                                      +v3130*e(221)+v3210*e(229)
     +                                      +v3220*e(233)+v3230*e(237)
     +                                      +v3310*e(245)+v3320*e(249)
     +                                      +v3330*e(253))*cppps
                                         ve11 = ve11 +
     +                                      (v0110*e(85)+v0120*e(89)
     +                                      +v0130*e(93)+v0210*e(101)
     +                                      +v0220*e(105)+v0230*e(109)
     +                                      +v0310*e(117)+v0320*e(121)
     +                                      +v0330*e(125))*cppps
                                         ve12 = ve12 +
     +                                      (v1010*e(85)+v1020*e(89)
     +                                      +v1030*e(93)+v2010*e(149)
     +                                      +v2020*e(153)+v2030*e(157)
     +                                      +v3010*e(213)+v3020*e(217)
     +                                      +v3030*e(221))*cppps
                                         ve13 = ve13 +
     +                                      (v1100*e(85)+v1200*e(101)
     +                                      +v1300*e(117)+v2100*e(149)
     +                                      +v2200*e(165)+v2300*e(181)
     +                                      +v3100*e(213)+v3200*e(229)
     +                                      +v3300*e(245))*cppps
                                         ve21 = ve21 +
     +                                      (v0110*e(149)+v0120*e(153)
     +                                      +v0130*e(157)+v0210*e(165)
     +                                      +v0220*e(169)+v0230*e(173)
     +                                      +v0310*e(181)+v0320*e(185)
     +                                      +v0330*e(189))*cppps
                                         ve22 = ve22 +
     +                                      (v1010*e(101)+v1020*e(105)
     +                                      +v1030*e(109)+v2010*e(165)
     +                                      +v2020*e(169)+v2030*e(173)
     +                                      +v3010*e(229)+v3020*e(233)
     +                                      +v3030*e(237))*cppps
                                         ve23 = ve23 +
     +                                      (v1100*e(89)+v1200*e(105)
     +                                      +v1300*e(121)+v2100*e(153)
     +                                      +v2200*e(169)+v2300*e(185)
     +                                      +v3100*e(217)+v3200*e(233)
     +                                      +v3300*e(249))*cppps
                                         ve31 = ve31 +
     +                                      (v0110*e(213)+v0120*e(217)
     +                                      +v0130*e(221)+v0210*e(229)
     +                                      +v0220*e(233)+v0230*e(237)
     +                                      +v0310*e(245)+v0320*e(249)
     +                                      +v0330*e(253))*cppps
                                         ve32 = ve32 +
     +                                      (v1010*e(117)+v1020*e(121)
     +                                      +v1030*e(125)+v2010*e(181)
     +                                      +v2020*e(185)+v2030*e(189)
     +                                      +v3010*e(245)+v3020*e(249)
     +                                      +v3030*e(253))*cppps
                                         ve33 = ve33 +
     +                                      (v1100*e(93)+v1200*e(109)
     +                                      +v1300*e(125)+v2100*e(157)
     +                                      +v2200*e(173)+v2300*e(189)
     +                                      +v3100*e(221)+v3200*e(237)
     +                                      +v3300*e(253))*cppps
                                         if (jtype.eq.6) then
c
c        (pp,pp), (pp,sp), (ps,pp), (sp,pp), (ss,pp), (sp,sp), (ps,sp),
c        and (ss,sp) section
c
                                         call fpppp
                                         end if
                                         end if
 200                                     if (oflag) then
c
c        end of the second pass through the integral section
c
                                         qve00 = qa*(ve00+ve00)
                                         dvex = -(ve11+ve12)
     +                                      *qa2 + (ve13+ve14)
     +                                      *qa1 - pqx*qve00
                                         dvey = -(ve21+ve22)
     +                                      *qa2 + (ve23+ve24)
     +                                      *qa1 - pqy*qve00
                                         dvez = -(ve31+ve32)
     +                                      *qa2 + (ve33+ve34)
     +                                      *qa1 - pqz*qve00
                                         dvexs = dvexs + dvex
                                         dveys = dveys + dvey
                                         dvezs = dvezs + dvez
                                         go to 230
                                         else
c
c   ******************************************************************
c
c        end of the first pass through the integral section
c
c   ******************************************************************
c
                                         oflag = .true.
                                         fq0 = fq1
                                         fq1 = fq2
                                         fq2 = fq3
                                         fq3 = fq4
                                         fq4 = fq5
                                         ve00s = ve00s + ve00
                                         ve11s = ve11s + ve11
                                         ve21s = ve21s + ve21
                                         ve31s = ve31s + ve31
                                         ve12s = ve12s + ve12
                                         ve22s = ve22s + ve22
                                         ve32s = ve32s + ve32
                                         cdve00 = s34*(ve00+ve00)
                                         dx2x = -ve13*s4 + ve14*s3 -
     +                                      cdx*cdve00
                                         dx2y = -ve23*s4 + ve24*s3 -
     +                                      cdy*cdve00
c
c        branch to second pass
c
                                         dx2z = -ve33*s4 + ve34*s3 -
     +                                      cdz*cdve00
                                         go to 190
                                         end if
c
c        special (ss,ss) section
c
 210                                     ve00 = fq0*eoooo
                                         ve00s = ve00s + ve00
                                         cdve00 = s34*(ve00+ve00)
                                         dx2x = -cdx*cdve00
                                         dx2y = -cdy*cdve00
                                         dx2z = -cdz*cdve00
                                         qve00 = qa*fq1*(eoooo+eoooo)
                                         dvex = -pqx*qve00
                                         dvey = -pqy*qve00
                                         dvez = -pqz*qve00
                                         dvexs = dvexs + dvex
                                         dveys = dveys + dvey
                                         dvezs = dvezs + dvez
                                         go to 230
c
c        special (ps,ss) section
c
 220                                     cpsss = cpss*csd
                                         exooo = e(65)*cpsss
                                         eyooo = e(129)*cpsss
                                         ezooo = e(193)*cpsss
                                         temp1 = ppxo*exooo +
     +                                      ppyo*eyooo + ppzo*ezooo
                                         temp2 = -
     +                                      (pqx*exooo+pqy*eyooo+pqz*
     +                                      ezooo)*qa2
                                         ve00 = fq0*(eoooo+temp1)
     +                                      + fq1*temp2
                                         ve00s = ve00s + ve00
                                         ve11s = ve11s + fq0*exooo
                                         ve21s = ve21s + fq0*eyooo
                                         ve31s = ve31s + fq0*ezooo
                                         cdve00 = s34*(ve00+ve00)
                                         dx2x = -cdx*cdve00
                                         dx2y = -cdy*cdve00
                                         dx2z = -cdz*cdve00
                                         ve00 = fq1*(eoooo+temp1)
     +                                      + fq2*temp2
                                         qve00 = qa*(ve00+ve00)
                                         temp = -qa2*fq1
                                         dvex = -pqx*qve00 + temp*exooo
                                         dvey = -pqy*qve00 + temp*eyooo
                                         dvez = -pqz*qve00 + temp*ezooo
                                         dvexs = dvexs + dvex
                                         dveys = dveys + dvey
                                         dvezs = dvezs + dvez
c
c        summation of contributions from the uncontracted gaussians
c
 230                                     fxk = fxk + dx2x - dvex*s3
                                         fyk = fyk + dx2y - dvey*s3
                                         fzk = fzk + dx2z - dvez*s3
 240                                     continue
 250                                     continue
                                         vee = vee + ve00s
                                         abve00 = s12*(ve00s+ve00s)
                                         dx1x = -ve11s*s2 + ve12s*s1 -
     +                                      abx*abve00
                                         dx1y = -ve21s*s2 + ve22s*s1 -
     +                                      aby*abve00
                                         dx1z = -ve31s*s2 + ve32s*s1 -
     +                                      abz*abve00
                                         fxi = fxi + dx1x + dvexs*s1
                                         fyi = fyi + dx1y + dveys*s1
                                         fzi = fzi + dx1z + dvezs*s1
                                         fxj = fxj - dx1x + dvexs*s2
                                         fyj = fyj - dx1y + dveys*s2
                                         fzj = fzj - dx1z + dvezs*s2
                                         end if
 260                                     continue
 270                                     continue
c
c        summation of the contributions from the shells
c
                                         fxl = -(fxi+fxj+fxk)
                                         fyl = -(fyi+fyj+fyk)
                                         fzl = -(fzi+fzj+fzk)
                                         fxyz(iatx) = fxyz(iatx) + fxi
                                         fxyz(jatx) = fxyz(jatx) + fxj
                                         fxyz(katx) = fxyz(katx) + fxk
                                         fxyz(latx) = fxyz(latx) + fxl
                                         fxyz(iaty) = fxyz(iaty) + fyi
                                         fxyz(jaty) = fxyz(jaty) + fyj
                                         fxyz(katy) = fxyz(katy) + fyk
                                         fxyz(laty) = fxyz(laty) + fyl
                                         fxyz(iatz) = fxyz(iatz) + fzi
                                         fxyz(jatz) = fxyz(jatz) + fzj
                                         fxyz(katz) = fxyz(katz) + fzk
                                         fxyz(latz) = fxyz(latz) + fzl
                                         if (lprint.eq.-4)
     +                                      write (iwr,6010) ishell ,
     +                                      jshell , kshell , lshell ,
     +                                      fxi , fxj , fxk , fxl ,
     +                                      fyi , fyj , fyk , fyl ,
     +                                      fzi , fzj , fzk , fzl
                                         end if
                                       end if
                                    end if
                                 end if
                              end if
 280                       continue
                        end if
 290                 continue
 3001                continue
_IF(parallel)
c***   **MPP**
                   next = ipg_dlbtask()
                endif
c***   **MPP**
_ENDIF
                  end if
               end if
 300        continue
            call dfinal(core,0,ishell)
c
c     write(iwr,9001) ishell,jst,kst,lst,tx,tim
c
            if (tim.gt.timlim) then
               ist = ishell + 1
               go to 320
            end if
         end if
 310  continue
 320  continue
_IF(parallel)
      call pg_dlbpush
_ENDIF
c
      return
c
 6010 format (1x,4i2,12f10.6)
 6020 format (/40x,'prefactor matrix in 2e-derivatives'/)
c9001 format(8x,'  ist,jst,kst,lst = ',4(1x,i3),6x,
c    +'del(time) = ',f10.3,2x,'time = ',f10.3)
      end
_IF(hpux11)
c$HP$ OPTIMIZE LEVEL3
_ENDIF
      subroutine ver_drv80(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/drv80.m,v $
     +     "/
      data revision /"$Revision: 6176 $"/
      data date /"$Date: 2010-08-10 16:49:47 +0200 (Tue, 10 Aug 2010) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
