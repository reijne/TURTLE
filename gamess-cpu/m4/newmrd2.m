c 
c  $Author: jmht $
c  $Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
c  $Locker:  $
c  $Revision: 5774 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/newmrd2.m,v $
c  $State: Exp $
c  
c ******************************************************
c ******************************************************
c             =   table = newmrd2    =
c ******************************************************
c ******************************************************
c
c ####################################################################
c #                                                                  #
c #   generator of formulatape for table-ci                          #
c #                                                                  #
c #   cardinput for large table (up to septets)                      #
c #   7    5    5    4    4    3    3    3                           #
c #                                                                  #
c ####################################################################
      subroutine tablegen(debug)
      implicit REAL  (a-h,o-z)
      logical debug
_IF(notused)
      logical*1 car(4),cart(4000)
_ELSE
      integer ot
      common/craypk/ot(4000),ic,i2,i3
      common/bufd/gout(510),ne,nspace
      common/blksi3/nsz,nsz51(6),nsz510
      dimension lout(512)
      equivalence (lout(1),gout(1))
_ENDIF
c
      dimension isq(10),jerk(10),iloc(5),khog(5),jbne(10),nrd(5)
      dimension nplu(5),imar(5)
      dimension imike(4),ijog(10),jjog(10),itim(10),jhog(67),ipar(45)
      dimension imox(90),ibil(3)
c
      common /junk/ a(15876)
INCLUDE(common/iofile)
INCLUDE(common/prints)
      integer  nston,mtape,mdisk
      integer  ideli,ltype,linf,ntab,kfile,kclab
      integer  ntype, mstvt, nf01, nf62, nhead
      integer  nf99, mtapev, nf11
      common /ftap/ nston, mtape, mdisk,
     .              ideli, ltype, linf,  ntab,  kfile,
     .              kclab, ntype, mstvt, nf01, nf62,
     +              nhead, nf99, mtapev, nf11
c
      REAL af, ae
_IF(notused)
      dimension ae(7829)
      common/scrtch/af(8010),icart(1000)
      equivalence (af(26),ae(1))
      integer nsac(5),ndat(5),iaz(5),iez(5),jan(7),jbn(7),idra(5)
      integer idrc(5),jdrc(5),j9(49)
      equivalence (af(1),j9(1)),(j9(1),ndat),(j9(6),nsac),
     1(j9(11),iaz),(j9(16),iez),(j9(21),idra),    
     2(j9(26),idrc),(j9(31),jdrc),(j9(36),jan),(j9(43),jbn)
_ELSE
      integer nsac,ndat,iaz,iez,jan,jbn,idra
      integer idrc,jdrc,ispace
      common/scrtch/ndat(5),nsac(5),iaz(5),iez(5),idra(5),
_IF(cray,ksr,i8)
     +              idrc(5),jdrc(5),jan(7),jbn(7),
_ELSE
     +              idrc(5),jdrc(5),jan(7),jbn(7),ispace,
_ENDIF
     +              ae(7829)
      dimension af(8010)
      equivalence (af(1),ndat(1))
_ENDIF
c
      common/scra/e(5292),imat(504),jmat(1439),icol(126),irow(126)
      dimension iper(4),ihog(48)
c
_IF(notused)
      equivalence (icar,car), (cart(1),icart(1)),
     *  (cart(1),ipos1),(cart(5),ipos2),(cart(9),ipos3)
_ENDIF
c
      write (iwr,12) cpulft(1)
  12  format(/40x,40('*')/
     *40x,'*** MRD-CI V2.0: Table Generation Module'/
     *40x,40('*')//
     + 1x,'***  start of table generation at ',f10.2,' seconds'/)
c
_IF(notused)
      icar=0                                                   
      mason=48
_ENDIF
      ikar=12
      nmax=10
      ien=ikar+1
 900  format(2x,40i3)
      call vclr(af,1,8010)
      ibb=4000
      iwd=ibb/8
      ibc=1-iwd
      imax=5
      kmjx=48
      isrp=126
      ept=1.0d-15
      eps=0.001d0
      write (iwr,4400)
4400  format(/2x,
     +'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'/
     +2x,
     +'+       This is the generator module for the CI Table.       +'/
     +2x,
     +'+ The author of this program is Professor Robert J. Buenker, +'/
     +2x,
     +'+                   University of Wuppertal                  +'/
     +2x,
     +'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     +  //)
c
_IF(notused)
c     bufferlength in bytes
c
      call dwopen(ntab,ibb)
_ELSE
      call setbfc
_ENDIF
      irec=0
c ####################################################################
c     itape=5
c     read(itape,1)mmax,(isq(i),i=1,mmax)
c #####################################################################
       mmax=7
       isq(1)=5    
       isq(2)=5    
       isq(3)=4    
       isq(4)=4    
       isq(5)=3    
       isq(6)=3    
       isq(7)=3
      write(iwr,900) mmax,(isq(i),i=1,mmax)
c #################
1     format(14i5)
c #################
      if(mmax.gt.nmax) 
     +  call caserr('invalid mmax in table generator')
      do 2 nmul=1,mmax
      jerk(nmul)=irec
      iswh=isq(nmul)
      if(iswh.gt.imax) go to 7781
      is=0
      ispin=(nmul-1)/2
      fms=dfloat(ispin)
      if (nmul.eq.2*(nmul/2)) fms=fms+0.5d0
      ibl=0
      nod=nmul-3
      max=0
      do 100 iw=1,iswh
      nod=nod+2
      nrd(iw)=nod
      write (iwr,241) nmul, nod
241   format(/1x,'case : ',i2,5x,'nod = ',i2)
      nmns=iw-1
      nplu(iw)=nod-nmns
      if (nmns.ne.0) go to 101
      ndat(1)=1
      nsac(1)=1
      khog(1)=0
      ae(1)=1.0d0
      ae(2)=1.0d0
      iaz(1)=-1
      iez(1)=1
      jf=2
      imar(1)=1
      ihg=1
      jhog(1)=1
      iloc(1)=0
      if (nod.eq.0) go to 100
      do 80 j=1,nod
      is=is+1
  80  jmat(is)=1
      go to 100
101   ndet=nod
      if (nmns.gt.1) go to 300
      do 301 j=1,nod
  301 imat(j)=j
      go to 125
  300 nmon=1
      nmot=nod
      do 126 j=2,nmns
      nmot=nmot-1
      ndet=ndet*nmot
126   nmon=nmon*j
      ndet=ndet/nmon
      if (ndet.gt.isrp) go to 954
      do 128 j=1,nmns
      iper(j)=j
128   imat(j)=j
      ib=nmns
134   jb=nmns
      jm=nod
129   km=iper(jb)
      if (km.ne.jm) go to 130
      jb=jb-1
      if (jb.eq.0) go to 131
      jm=jm-1
      go to 129
130   ip=iper(jb)+1
      iper(jb)=ip
      if (jb.eq.nmns) go to 132
      jb=jb+1
      do 133 j=jb,nmns
      ip=ip+1
133   iper(j)=ip
132   do 135 j=1,nmns
      ib=ib+1
135   imat(ib)=iper(j)
      go to 134
131   if (ib.ne.ndet*nmns) go to 955
125   nx=nod*ndet
      nd=is+1
      ne=is+nx
      do 531 j=nd,ne
531   jmat(j)=1
      jh=0
      iloc(iw)=is
      do 532 j=1,ndet
      do 533 k=1,nmns
      jh=jh+1
      ly=imat(jh)+is
533   jmat(ly)=-1
532   is=is+nod
      fnod=dfloat(nod)
      khog(iw)=ihg
      fnod=fnod*0.5d0-fms
      nin=ndet+1
      nnam=-ndet
      do 136 j=1,ndet
      nnam=nnam+nin
136   a(nnam)=fnod
      ib=0
      kb=0
      do 137 j=2,ndet
      ib=ib+nmns
      j1=j-1
      jb=-nmns
      kb=kb+ndet
      lb=kb
      ilk=-ndet+j
      do 137 k=1,j1
      ilk=ilk+ndet
      lb=lb+1
      if (nmns.eq.1) go to 139
      mx=0
      jb=jb+nmns
      jn=jb
      do 138 l=1,nmns
      jn=jn+1
      ip=imat(jn)
      in=ib
      do 140 ll=1,nmns
      in=in+1
      if (imat(in)-ip) 140,138,141
140   continue
141   if (mx.eq.1) go to 142
      mx=1
138   continue
139   a(lb)=1.0d0
      go to 137
142   a(lb)=0.0d0
137   a(ilk)=a(lb)
      call dmfgrt(a,ndet,ndet,eps,irank,irow,icol)
      kmj=ndet-irank
      if (kmj.gt.kmjx) go to 950
      ifac=nod+1-nmns
      ifac=(ifac-nmns)*ndet/ifac
      if (ifac.ne.kmj) go to 951
4320  format(2x,40i3)
      if (debug) write (iwr,4320) ndet,kmj
      if (debug) write (iwr,4320) (icol(j),j=1,ndet)
      js=irank*ndet
      in=irank
      do 143 i=1,kmj
      in=in+1
      sum=1.0d0
      do 144 j=1,irank
      ico=(icol(j)-1)*ndet+in
      js=js+1
      ab=a(js)
      sum=sum+ab*ab
  144 a(ico)=ab
      k=irank
      js=js+kmj
      do 145 j=1,kmj
      k=k+1
      ico=(icol(k)-1)*ndet+in
      if (i.eq.j) go to 1460
      a(ico)=0.0d0
      go to 145
1460  a(ico)=1.0d0
145   continue
      k=in-ndet
      sum=1.0d0/ sqrt(sum)
      do 146 j=1,ndet
      k=k+ndet
146   a(k)=a(k)*sum
      if (i.eq.1) go to 143
      i1=i-1
      kb=1-ndet
      do 147 j=1,i1
      kb=kb+ndet
147   a(kb)=0
      jn=in-ndet
      kn=irank-ndet
      do 148 j=1,ndet
      jn=jn+ndet
      cj=a(jn)
      kb=1-ndet
      kn=ndet+kn
      ln=kn
      do 148 k=1,i1
      ln=ln+1
      kb=kb+ndet
148   a(kb)=a(kb)+cj*a(ln)
      kb=1-ndet
      sum=0.0d0
      do 149 k=1,i1
      kb=kb+ndet
      cj=a(kb)
149   sum=sum+cj*cj
      sum=-1.0d0/ sqrt(1.0d0-sum)
      kb=1-ndet
      do 150 j=1,i1
      kb=kb+ndet
150   a(kb)=a(kb)*sum
      kb=kb+ndet
      a(kb)=-sum
      nb=i-kmj
      lb=-kmj
      do 151 j=1,ndet
      sum=0.0d0
      lb=lb+ndet
      jb=lb
      kb=1-ndet
      do 152 k=1,i
      kb=kb+ndet
      jb=jb+1
152   sum=sum+a(jb)*a(kb)
      nb=nb+ndet
151   a(nb)=sum
143   continue
      ir=irank
      do 250 i=1,kmj
      ir=ir+1
250   ihog(i)=icol(ir)
      ir=0
      nb=-kmj
      do 154 i=1,kmj
      nb=nb+1
      kb=nb
      do 154 j=1,ndet
      ir=ir+1
      kb=kb+ndet
154   e(ir)=a(kb)
      jb=kmj*ndet
      if (debug) write (iwr,4322) (e(j),j=1,jb)
4322  format(2x,15f8.4)
      imar(iw)=jb
      ndat(iw)=ndet
      nsac(iw)=kmj
      iaz(iw)=jf-ndet
      do 534 j=1,jb
      jf=jf+1
  534 ae(jf)=e(j)
      iez(iw)=jf
      do 535 j=1,kmj
      ihg=ihg+1
  535 jhog(ihg)=ihog(j)
      jb=-kmj
      do 157 i=1,kmj
      ix=ihog(i)
      ib=ix*ndet-kmj
      jb=jb+1
      kb=jb
      do 157 j=1,kmj
      ib=ib+1
      kb=kb+kmj
157   e(kb)=a(ib)
      ie=-kmj
      jp=-kmj
      il=-kmj
      do 158 i=1,kmj
      il=il+kmj
      jl=il
      jp=jp+1
      kp=jp
      ie=ie+1
      je=-kmj
      do 158 j=1,i
      jl=jl+1
      kp=kp+kmj
      sum=0.0d0
      je=je+1
      jk=je
      ik=ie
      do 159 k=1,kmj
      ik=ik+kmj
      jk=jk+kmj
159   sum=sum+e(ik)*e(jk)
      a(kp)=sum
158   a(jl)=sum
      call dgelgt(e,a,kmj,kmj,ept,ier)
      if (ier.eq.0) go to 160
      write(iwr,161) ier
161   format(10x,'ier error=',i2)
160   jl=kmj*kmj
      do 536 j=1,jl
      jf=jf+1
 536  ae(jf)=e(j)
      if (debug) write (iwr,4322)  (e(j),j=1,jl)
100   continue
_IF(notused)
      nq=(mason +jf  )/iwd+1
_ELSE
      nq=(jf-1)/nsz510+2
_ENDIF
      jbne(nmul)=nq
      irec=irec+nq
      jrec=irec
      do 537 j=1,iswh
      nod=nrd(j)
      kmj=nsac(j)
      if(j.lt.3) go to 538
      ibl=ibl+1
      jan(ibl)=irec
      nd=nod*(nod-1)*(nod-2)*(nod-3)/24
      nd=nd*(1+kmj+kmj)
      nd=(nd-1)/ibb+1
      jbn(ibl)=nd
      irec=irec+nd
538   if(j.eq.1) go to 539
      ibl=ibl+1
      nd=ikar
      if(nod.lt.3) go to 540
      ne=nod*(nod-1)*(nod-2)*(nod-2)/6
      nd=ne*(1+kmj+kmj)+nd
540   ne=nod*(nod-1)/2
      nd=ne*kmj+nd
      nf=(nod+nmul)*kmj+1
      nd=nf*ne+nd
      nd=(nd-1)/ibb+1
      jan(ibl)=irec
      jbn(ibl)=nd
      irec=irec+nd
539   idra(j)=irec
      nps=nplu(j)
      nms=j-1
      if(nps.lt.2) nps=0
      if(nms.lt.2) nms=0
      nps=kmj*(nms+nps)+3*nms*nps*kmj+kmj+ikar
      nd=(nps-1)/ibb+1
      idrc(j)=nd
      irec=irec+nd
      nd=kmj+ikar
      if(nod.lt.2) go to 541
      nps=nod*(nod-1)/2
      nps=nps*(nps+1)/2
      nps=nps*(3*kmj+1)
      nd=nd+nps
541   if(nod.eq.0) go to 542
      nps=nod*(nod+1)/2
      nd=nd+nps*(1+kmj)
      nd=nd+nps*(2*nod*kmj+1)
542   nd=(nd-1)/ibb+1
      jdrc(j)=nd
537   irec=irec+nd
c
_IF(notused)
     ne=ibc
     do 543 j=1,nq
     ne=ne+iwd
     call dwritet(af(ne))
543  continue
_ELSE
      call stopbk3
      jjj=jerk(nmul)
_IF(cray,ksr)
      call fmove(ndat,lout,49)
_ELSE
      call icopy(49,ndat,1,lout,1)
_ENDIF
      call sttout3(jjj)
      call stopbk3
      jjj=jjj+nsz
      ne=0
      do 543 j=1,jf
      ne=ne+1
      gout(ne)=ae(j)
      if(ne.ne.nsz510)go to 543
      call sttout3(jjj)
      call stopbk3
      ne=0
      jjj =jjj +nsz
 543  continue
      if(ne.eq.0)go to 5555
      call sttout3(jjj)
      call stopbk3
_ENDIF
 5555 if (debug) write (iwr,75) ndat,nsac,iaz,iez,idra,idrc,
     +            jdrc,jan,jbn
  75  format(2x,20i6)
      jbl=0
      nmns=-1
      do 551 j=1,iswh
_IF(notused)
_ELSE
      ic=0
      i2=0
      i3=0
      call pack12
_ENDIF
      ihg=khog(j)
      nqns=nmns
      nmns=nmns+1
      nod=nrd(j)
      ndet=ndat(j)
      kmj=nsac(j)
      ii=kmj+kmj+1
      kr=iloc(j)-nod
      jps=nplu(j)
      kmt=kmj*jps
      if(j.lt.3) go to 552
      jbl=jbl+1
      mm=0
      nodj=nod-4
      j8=j-2
      ndj=ndat(j8)
      jr=iloc(j8)-nodj
      do 553 k=1,4
553   imike(k)=k
      nix=0
      if(nodj.eq.0) go to 558
561   ipr=1
      ips=imike(1)
      do 554 k=1,nod
554   itim(k)=k
      do 555 k=1,nod
      k1=k+1
590   ico=itim(k)
      if(ico.ne.ips) go to 555
      nix=nix+nod-k
      do 559 l=k1,nod
559   itim(l-1)=itim(l)
      itim(nod)=ips
      if(ipr.eq.4) go to 558
      ipr=ipr+1
      ips=imike(ipr)
      go to 590
555   continue
558   mm=mm+1
      ico=1
_IF(notused)
      if(nix-2*(nix/2).gt.0) ico=255
      icar=ico
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 560
      call dwritet(cart)
      jrec=jrec+1
_ENDIF
      if(nix-2*(nix/2).gt.0) ico=-1
      ot(mm)=ico
      if(mm.lt.ibb) go to 560
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      jrec=jrec+nsz
      mm=0
560   jhg=ihg
      do 562 k=1,kmj
      mm=mm+1
      jhg=jhg+1
      khg=jhog(jhg)
      nd=khg*nod+kr
      nps=0
      ipr=1
      ico=imike(1)
      do 563 l=1,nod
      nd=nd+1
      if(l.eq.ico) go to 564
      nps=nps+1
      ijog(nps)=jmat(nd)
      go to 563
564   itim(ipr)=jmat(nd)
      if(ipr.eq.4) go to 563
      ipr=ipr+1
      ico=imike(ipr)
563   continue
      if(itim(1).gt.0) go to 565
      if(itim(2).gt.0) go to 566
      if(itim(3).gt.0) go to 567
_IF(notused)
569   icar=0
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 568
      mm=0
      jrec=jrec+1
      call dwritet(cart)
_ELSE
569   ot(mm)=0
      if(mm.lt.ibb) go to 568
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
_ENDIF
568   mm=mm+1
      if(mm.lt.ibb) go to 562
_IF(notused)
      mm=0
      jrec=jrec+1
      call dwritet(cart)
_ELSE
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
_ENDIF
      go to 562
567   if(itim(4).lt.0) go to 569
      ico=1
      go to 570
566   if(itim(3).gt.0) go to 571
      if(itim(4).lt.0) go to 569
      ico=2
      go to 570
571   if(itim(4).gt.0) go to 569
      ico=3
      go to 570
565   if(itim(2).lt.0) go to 572
      if(itim(3).gt.0) go to 569
      if(itim(4).gt.0) go to 569
      ico=1
      go to 570
572   if(itim(3).gt.0) go to 573
      if(itim(4).lt.0) go to 569
      ico=3
      go to 570
573   if(itim(4).gt.0) go to 569
      ico=2
570   if(nodj.gt.0) go to 574
_IF(notused)
      icar=1
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 575
      mm=0
      call dwritet(cart)
      jrec=jrec+1
575   mm=mm+1
      icar=ico
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 562
      mm=0
      jrec=jrec+1
      call dwritet(cart)
_ELSE
      ot(mm)=1
      if(mm.lt.ibb) go to 575
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
575   mm=mm+1
      ot(mm)=ico
      if(mm.lt.ibb) go to 562
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
_ENDIF  
      go to 562
574   nps=jr
      do 576 l=1,ndj
      nps=nps+nodj
      ipr=nps
      do 577 m=1,nodj
      ipr=ipr+1
      if (ijog(m).ne.jmat(ipr)) go to 576
577   continue
      go to 579
576   continue
      go to 7781
_IF(notused)
579   icar=l
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 580
      mm=0
      jrec=jrec+1
      call dwritet(cart)
580   mm=mm+1
      icar=ico
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 562
      mm=0
      jrec=jrec+1
      call dwritet(cart)
_ELSE
579   ot(mm)=l
      if(mm.lt.ibb) go to 580
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
580   mm=mm+1
      ot(mm)=ico
      if(mm.lt.ibb) go to 562
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
_ENDIF
562   continue
      if(nodj.eq.0) go to 616
      nix=0
      ipr=2
      ico=imike(1)+1
170   if(ico.eq.imike(ipr)) go to 581
      imike(ipr-1)=ico
      if (ipr.eq.2) go to 561
      ibd=ipr-2
      do 866 k=1,ibd
 866  imike(k)=k
      go to 561
581   if(ipr.eq.4) go to 582
      ico=ico+1
      ipr=ipr+1
      go to 170
 582  if (ico.eq.nod) go to 616
      imike(4)=ico+1
      do 867 k=1,3
 867  imike(k)=k
      go to 561
616   if(mm.eq.0) go to 584
_IF(notused)
      call dwritet(cart)
      mm=0
      jrec=jrec+1
_ELSE
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      jrec=jrec+nsz
_ENDIF
      go to 584
552   if(j.eq.1) go to 583
584   nodj=nod-2
      j8=j-1
      jr=iloc(j8)-nodj
      ndj=ndat(j8)
      jbl=jbl+1
_IF(notused)
      ipos1=ien-ii
_ELSE
      ic=ien-ii
_ENDIF
      mm=ikar
      nd=ikar
      if(nod.lt.3) go to 585
      ne=nod*(nod-1)*(nod-2)*(nod-2)/6
      nd=ne*ii+nd
_IF(notused)
585   ipos2=nd-kmj
_ELSE
585   i2=nd-kmj
_ENDIF
      jj=nod+nmul
      ne=nod*(nod-1)/2
_IF(notused)
      ipos3=nd+(ne-jj)*kmj
_ELSE
      i3=nd+(ne-jj)*kmj
      if(debug) write(6,99988) ic, i2, i3
99988 format(1x,'**** ic, i2, i3 = ', 3i10)
      call pack12
_ENDIF
      if(nod.lt.3) go to 586
      do 587 k=1,3
587   imike(k)=k
      nix=0
      if(nod.eq.3) go to 588
589   ipr=1
      ips=imike(1)
      do 162 k=1,nod
162   itim(k)=k
      do 592 k=1,nod
      k1=k+1
593   ico=itim(k)
      if(ico.ne.ips) go to 592
      nix=nix+nod-k
      do 163 l=k1,nod
163   itim(l-1)=itim(l)
      itim(nod)=ips
      if(ipr.eq.3) go to 588
      ipr=ipr+1
      ips=imike(ipr)
      go to 593
592   continue
588   ica=1
      if(nix-2*(nix/2).gt.0) ica=-1
      ipr=nodj
      do 594 k=1,nodj
      ipr=ipr-1
      jco=ica
      if(ipr-2*(ipr/2).gt.0) jco=-jco
      mm=mm+1
_IF(notused)
      if(jco.lt.0) jco=jco+256
      icar=jco
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 595
      mm=0
      jrec=jrec+1
      call dwritet(cart)
_ELSE
      ot(mm)=jco
      if(mm.lt.ibb) go to 595
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
_ENDIF
595   jhg=ihg
      do 596 l=1,kmj
      mm=mm+1
      jhg=jhg+1
      khg=jhog(jhg)
      nd=khg*nod+kr
      nps=0
      ipl=1
      ico=imike(1)
      do 597 m=1,nod
      nd=nd+1
      if(m.eq.ico) go to 598
      nps=nps+1
      if(nps.ne.k) go to 599
      nps=nps+1
599   ijog(nps)=jmat(nd)
      go to 597
598   itim(ipl)=jmat(nd)
      if(ipl.eq.3) go to 597
      ipl=ipl+1
      ico=imike(ipl)
597   continue
      if(itim(1).gt.0) go to 600
      if(itim(2).gt.0) go to 601
      if(itim(3).gt.0) go to 602
_IF(notused)
603   icar=0
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if (mm.lt.ibb) go to 604
      mm=0
      jrec=jrec+1
      call dwritet(cart)
604   mm=mm+1
      if(mm.lt.ibb) go to 596
      mm=0
      jrec=jrec+1
      call dwritet(cart)
_ELSE
603   ot(mm)=0
      if (mm.lt.ibb) go to 604
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
604   mm=mm+1
      if(mm.lt.ibb) go to 596
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
_ENDIF
      go to 596
602   ico=-1
      ijog(k)=-1
      go to 605
601   if(itim(3).gt.0) go to 606
      ico=-2
      ijog(k)=-1
      go to 605
606   ico=3
      ijog(k)=1
      go to 605
600   if(itim(2).gt.0) go to 607
      if(itim(3).gt.0) goto 608
      ico=-3
      ijog(k)=-1
      go to 605
608   ico=2
      ijog(k)=1
      go to 605
607   if(itim(3).gt.0) go to 603
      ico=1
      ijog(k)=1
605   nps=jr
      do 609 m=1,ndj
      nps=nps+nodj
      ipl=nps
      do 610 n=1,nodj
      ipl=ipl+1
      if(ijog(n).ne.jmat(ipl)) go to 609
610   continue
      go to 611
609   continue
      go to 7781
_IF(notused)
611   icar=m
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 612
      mm=0
      jrec=jrec+1
      call dwritet(cart)
612   mm=mm+1
      if(ico.lt.0) ico=ico+256
      icar=ico
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 596
      mm=0
      jrec=jrec+1
      call dwritet(cart)
_ELSE
611   ot(mm)=m
      if(mm.lt.ibb) go to 612
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
612   mm=mm+1
      ot(mm)=ico
      if(mm.lt.ibb) go to 596
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
_ENDIF
596   continue
594   continue
      if(nod.eq.3) go to 586
      nix=0
      ipr=2
      ico=imike(1)+1
613   if(ico.eq.imike(ipr)) go to 614
      imike(ipr-1)=ico
      if (ipr.eq.3) imike(1)=1
      go to 589
614   if(ipr.eq.3) go to 615
      ico=ico+1
      ipr=ipr+1
      go to 613
615   if(ico.eq.nod) go to 586
      imike(3)=ico+1
      imike(1)=1
      imike(2)=2
      go to 589
586   imike(1)=1
      imike(2)=2
      max=0
      mix=0
      nix=0
      if(nod.eq.2) go to 627
621   ipr=1
      ips=imike(1)
      do 622 k=1,nod
622   itim(k)=k
      do 624 k=1,nod
      k1=k+1
625   ico=itim(k)
      if(ico.ne.ips) go to 624
      nix=nix+nod-k
      do 626 l=k1,nod
626   itim(l-1)=itim(l)
      itim(nod)=ips
      if(ipr.eq.2) go to 627
      ipr=2
      ips=imike(2)
      go to 625
624   continue
627   ica=1
      if(nix-2*(nix/2).gt.0) ica=-1
      max=max+1
      ipar(max)=ica
      mix=mix+2
      imox(mix)=imike(2)
      imox(mix-1)=imike(1)
      jhg=ihg
      do 628 k=1,kmj
      mm=mm+1
      jhg=jhg+1
      khg=jhog(jhg)
      nd=khg*nod+kr
      nps=0
      ipr=1
      ico=imike(1)
      do 629 l=1,nod
      nd=nd+1
      if(l.eq.ico) go to 166
      nps=nps+1
      ijog(nps)=jmat(nd)
      go to 629
166   itim(ipr)=jmat(nd)
      if(ipr.eq.2) go to 629
      ipr=2
      ico=imike(2)
629   continue
      if(itim(1).gt.0) go to 631
      if(itim(2).gt.0) go to 632
_IF(notused)
634      icar=0
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 628
      mm=0
      jrec=jrec+1
      call dwritet(cart)
_ELSE
634   ot(mm)=0
      if(mm.lt.ibb) go to 628
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
_ENDIF
      go to 628
632   jco=-1
      go to 633
631   if(itim(2).gt.0) go to 634
      jco=1
633   if(ica.eq.-1) jco=-jco
      if(ndj.ne.1) go to 646
      l=1
      go to 637
646   nps=jr
      do 635 l=1,ndj
      nps=nps+nodj
      ipr=nps
      do 636 m=1,nodj
      ipr=ipr+1
      if(ijog(m).ne.jmat(ipr)) go to 635
636   continue
      go to 637
635   continue
      go to 7781
637   if(jco.eq.-1) l=-l
_IF(notused)
      if(l.lt.0) l=l+256
      icar=l
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 628
      mm=0
      jrec=jrec+1
      call dwritet(cart)
_ELSE
      ot(mm)=l
      if(mm.lt.ibb) go to 628
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
_ENDIF
628   continue
      if(nodj.eq.0) go to 630
      nix=0
      ico=imike(1)+1
      if(ico.eq.imike(2)) go to 165
      imike(1)=ico
      go to 621
165   if(ico.eq.nod) go to 630
      imike(2)=ico+1
      imike(1)=1
      go to 621
 630  mix=0
      do 638 k=1,max
      ica=ipar(k)
      mm=mm+1
_IF(notused)
      if(ica.lt.0) ica=ica+256
      icar=ica
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 640
      mm=0
      jrec=jrec+1
      call dwritet(cart)
_ELSE
      ot(mm)=ica
      if(mm.lt.ibb) go to 640
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
_ENDIF
640   nn=mm
      mix=mix+2
      imike(1)=imox(mix-1)
      imike(2)=imox(mix)
      jhg=ihg
      do 680 l=1,kmj
      mm=nn+1
      iba=jrec
      jhg=jhg+1
      khg=jhog(jhg)
      nd=khg*nod+kr
      nps=0
      ipr=1
      ico=imike(1)
      do 641 m=1,nod
      nd=nd+1
      if(m.eq.ico) go to 642
      nps=nps+1
      ijog(nps) =jmat(nd)
      go to 641
642   itim(ipr)=jmat(nd)
      if(ipr.eq.2) go to 641
      ipr=2
      ico=imike(2)
641   continue
      if(itim(1).lt.0) go to 643
      if(itim(2).gt.0) go to 644
_IF(notused)
      icar=1
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
_ELSE
      ot(mm)=1
_ENDIF
      if(ndj.ne.1)    go to 647
      m=1
      go to 648
647   nps=jr
      do 645 m=1,ndj
      nps=nps+nodj
      ipr=nps
      do 649 n=1,nodj
      ipr=ipr+1
      if(ijog(n).ne.jmat(ipr))go to 645
649   continue
      go to 648
645   continue
      go to 7781
648   if(mm.lt.ibb) go to 650
_IF(notused)
      mm=0
      jrec=jrec+1
      call dwritet(cart)
_ELSE
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
_ENDIF
650   mm=mm+1
_IF(notused)
      icar=m
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 651
      mm=0
      jrec=jrec+1
      call dwritet(cart)
_ELSE
      ot(mm)=m
      if(mm.lt.ibb) go to 651
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
_ENDIF
651   if(nodj.eq.0) go to 639
      do 653 m=1,nodj
      if (ijog(m).lt.0) go to 653
      mm=mm+1
_IF(notused)
      icar=m
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 653
      mm=0
      jrec=jrec+1
      call dwritet(cart)
_ELSE
      ot(mm)=m
      if(mm.lt.ibb) go to 653
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
_ENDIF
653   continue
      if(j8.eq.1) go to 639
      do 654 m=1,nodj
      if(ijog(m).gt.0) go to 654
      mm=mm+1
_IF(notused)
      icar=m
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 654
      mm=0
      jrec=jrec+1
      call dwritet(cart)
654   continue
      go to 639
643   if(itim(2).lt.0) go to 655
      icar=2
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 656
      mm=0
      jrec=jrec+1
      call dwritet (cart)
_ELSE
      ot(mm)=m
      if(mm.lt.ibb) go to 654
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
654   continue
      go to 639
643   if(itim(2).lt.0) go to 655
      ot(mm)=2
      if(mm.lt.ibb) go to 656
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
_ENDIF
656   if(ndj.ne.1) go to 657
      m=1
      go to 658
657   nps=jr
      do 659 m=1,ndj
      nps=nps+nodj
      ipr=nps
      do 660 n=1,nodj
      ipr=ipr+1
      if(ijog(n).ne.jmat(ipr)) go to 659
660   continue
      go to 658
659   continue
      go to 7781
658   mm=mm+1
_IF(notused)
      icar=m
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 661
      mm=0
      jrec=jrec+1
      call dwritet(cart)
_ELSE
      ot(mm)=m
      if(mm.lt.ibb) go to 661
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
_ENDIF
661   if(nodj.eq.0) go to 639
      if(j8.eq.1) go to 663
      do 662 m=1,nodj
      if(ijog(m).gt.0) go to 662
      mm=mm+1
_IF(notused)
      icar=m
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 662
      mm=0
      jrec=jrec+1
      call dwritet(cart)
_ELSE
      ot(mm)=m
      if(mm.lt.ibb) go to 662
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
_ENDIF
662   continue
663   do 664 m=1,nodj
      if(ijog(m).lt.0) go to 664
      mm=mm+1
_IF(notused)
      icar=m
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 664
      mm=0
      jrec=jrec+1
      call dwritet(cart)
_ELSE
      ot(mm)=m
      if(mm.lt.ibb) go to 664
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
_ENDIF
664   continue
      go to 639
_IF(notused)
644   icar=3
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 665
      mm=0
      jrec=jrec+1
      call dwritet(cart)
_ELSE
644   ot(mm)=3
      if(mm.lt.ibb) go to 665
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
_ENDIF
665   do 666 m=1,nodj
      if(ijog(m).gt.0) go to 666
      do 667 n=1,nodj
667   jjog(n)=ijog(n)
      jjog(m)=1
      nps=jr
      do 668 n=1,ndj
      nps=nps+nodj
      ipr=nps
      do 669 if=1,nodj
      ipr=ipr+1
      if(jjog(if).ne.jmat(ipr)) go to 668
669   continue
      go to 672
668   continue
      go to 7781
672   mm=mm+1
_IF(notused)
      icar=n
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 671
      mm=0
      jrec=jrec+1
      call dwritet(cart)
671   mm=mm+1
      icar=m
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 666
      mm=0
      jrec=jrec+1
      call dwritet(cart)
666   continue
      go to 639
655      icar=4
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 673
      mm=0
      jrec=jrec+1
      call dwritet(cart)
_ELSE
      ot(mm)=n
      if(mm.lt.ibb) go to 671
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
671   mm=mm+1
      ot(mm)=m
      if(mm.lt.ibb) go to 666
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
666   continue
      go to 639
655   ot(mm)=4
      if(mm.lt.ibb) go to 673
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
_ENDIF
673   do 674 m=1,nodj
      if(ijog(m).lt.0) go to 674
      do 675 n=1,nodj
675   jjog(n)=ijog(n)
      jjog(m)=-1
      nps=jr
      do 676 n=1,ndj
      nps=nps+nodj
      ipr=nps
      do 677 if=1,nodj
      ipr=ipr+1
      if(jjog(if).ne.jmat(ipr)) go to 676
677   continue
      go to 678
676   continue
      go to 7781
678   mm=mm+1
_IF(notused)
      icar=n
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 679
      mm=0
      jrec=jrec+1
      call dwritet(cart)
679   mm=mm+1
      icar=m
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 674
      mm=0
      jrec=jrec+1
      call dwritet(cart)
674   continue
639   nn=nn+jj
      if(nn.lt.ibb) go to 680
      nn=nn-ibb
      if (jrec.gt.iba) go to 680
      jrec=jrec+1
      call dwritet(cart)
680   continue
638   mm=nn
      if(mm.eq.0) go to 681
      jrec=jrec+1
      call dwritet(cart)
      go to 681
583      icar=1
_IF(littleendian)
      cart(ien)=car(1)
_ELSE
      cart(ien)=car(4)
_ENDIF
      if(jps.lt.2) go to 682
      ipos2=1
      mm=ien
      do 683 k=1,jps
      mm=mm+1
      icar=k
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
683   continue
_ELSE
      ot(mm)=n
      if(mm.lt.ibb) go to 679
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
679   mm=mm+1
      ot(mm)=m
      if(mm.lt.ibb) go to 674
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
674   continue
639   nn=nn+jj
      if(nn.lt.ibb) go to 680
      nn=nn-ibb
      if (jrec.gt.iba) go to 680
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      jrec=jrec+nsz
680   continue
638   mm=nn
      if(mm.eq.0) go to 681
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      jrec=jrec+nsz
      go to 681
 583  ot(13)=1
      if(jps.lt.2) go to 682
      i2=1
      mm=13
      if(debug) write(6,99988) ic, i2, i3
      call pack12
      do 683 k=1,jps
      mm=mm+1
683   ot(mm)=k
_ENDIF
      mix=0
      ipr=1
      jpr=2
 798  mix=mix+2
      imox(mix-1)=ipr
      imox(mix)=jpr
      nix=ipr+jpr
      ica=1
      if (nix-2*(nix/2).eq.0) ica=-1
      max=max+1
      ipar(max)=ica
      ipr=ipr+1
      if (ipr.lt.jpr) go to 798
      if (jpr.eq.nod) go to 682
      jpr=jpr+1
      ipr=1
      go to 798
_IF(notused)
681   ipos1=kmj
      i2=kmj+jps*nmns*3*kmj
      ipos2=i2
      ipos3=kmt+i2
      mm=ikar
_ELSE
 681  ic=kmj
      i2=kmj+jps*nmns*3*kmj
      i3=i2+kmt
      if(debug) write(6,99988) ic, i2, i3
      mm=ikar
      call pack12
_ENDIF
      ipr=kr
      icy=0
      jhg=ihg
      do 302 k=1,kmj
      mm=mm+1
      jhg=jhg+1
_IF(notused)
      icar=jhog(jhg)
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
302   continue
_ELSE
 302  ot(mm)=jhog(jhg)
_ENDIF
      do 684 k=1,ndet
      ipr=ipr+nod
      ico=ipr
      do 685 l=1,nod
      ico=ico+1
      if(jmat(ico).gt.0) go to 685
      icy=icy+1
      imat(icy)=l
685   continue
684   continue
      icy=ihg
      do 857 if=1,kmj
      icy=icy+1
      ih=jhog(icy)
      il=(ih-1)*nmns
      jl=-nmns
      do 857 k=1,ndet
      jl=jl+nmns
      if (ih.eq.k) go to 857
      kl=il
      nix=0
      mix=0
      do 859 l=1,nmns
      kl=kl+1
      kzl=imat(kl)
      ll=jl
      do 860 kk=1,nmns
      ll=ll+1
      if(kzl-imat(ll)) 861,862,860
862   mix=mix+1
      ibil(mix)=kk
      go to 859
  860 continue
861   if (nix.eq.1) go to 857
      nix=1
      mal=kzl
859   continue
      if (j.gt.2) go to 791
      nal=imat(jl+1)
      go to 865
 791  do 864 l=1,nqns
      if(ibil(l).eq.l) go to 864
      nal=imat(jl+l)
      go to 865
864   continue
      nal=imat(jl+nmns)
865   mm=mm+1
_IF(notused)
      icar=k
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if (mm.lt.ibb) go to 792
      mm=0
      jrec=jrec+1
      call dwritet(cart)
 792  mm=mm+1
      if (mal.gt.nal) go to 793
      icar=nal
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if (mm.lt.ibb) go to 794
      mm=0
      jrec=jrec+1
      call dwritet(cart)
 794  mm=mm+1
      icar=mal
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if (mm.lt.ibb) go to 857
      mm=0
      jrec=jrec+1
      call dwritet(cart)
      go to 857
793      icar=mal
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if (mm.lt.ibb) go to 795
      mm=0
      jrec=jrec+1
      call dwritet(cart)
 795  mm=mm+1
      icar=nal
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if (mm.lt.ibb) go to 857
      mm=0
      jrec=jrec+1
      call dwritet(cart)
_ELSE
      ot(mm)=k
      if (mm.lt.ibb) go to 792
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
 792  mm=mm+1
      if (mal.gt.nal) go to 793
      ot(mm)=nal
      if (mm.lt.ibb) go to 794
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
 794  mm=mm+1
      ot(mm)=mal
      if (mm.lt.ibb) go to 857
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
      go to 857
 793  ot(mm)=mal
      if (mm.lt.ibb) go to 795
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
 795  mm=mm+1
      ot(mm)=nal
      if (mm.lt.ibb) go to 857
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
_ENDIF
857   continue
      if(jps.lt.2) go to 688
      jhg=ihg
      do 689 k=1,kmj
      jhg=jhg+1
      nd=jhog(jhg)*nod+kr
      do 690 l=1,nod
      nd=nd+1
      if(jmat(nd).lt.0) go to 690
      mm=mm+1
_IF(notused)
      icar=l
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 690
      mm=0
      jrec=jrec+1
      call dwritet(cart)
_ELSE
      ot(mm)=l
      if(mm.lt.ibb) go to 690
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
_ENDIF
690   continue
689   continue
      if(j.eq.2) go to 688
      jhg=ihg
      do 691 k=1,kmj
      jhg=jhg+1
      nd=jhog(jhg)*nod+kr
      do 692 l=1,nod
      nd=nd+1
      if(jmat(nd).gt.0) go to 692
      mm=mm+1
_IF(notused)
      icar=l
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 692
      mm=0
      jrec=jrec+1
      call dwritet(cart)
692   continue
691   continue
688   if(mm.eq.0) go to 693
682   call dwritet(cart)
      jrec=jrec+1
_ELSE
      ot(mm)=l
      if(mm.lt.ibb) go to 692
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
692   continue
691   continue
688   if(mm.eq.0) go to 693
682   call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      jrec=jrec+nsz
_ENDIF
693   mm=ikar
      jhg=ihg
      do 694 k=1,kmj
      jhg=jhg+1
      mm=mm+1
_IF(notused)
      icar=jhog(jhg)
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
694   continue
      if(nod.gt.1) go to 695
      if(nod.eq.0) go to 797
      ipos2=mm-kmj
      jj=nod+nod
      ipos3=mm+kmj+1-jj*kmj
      go to 735
797   call dwritet(cart)
      jrec=jrec+1
      go to 551
695   ii=ii+kmj
      ipos1=mm-ii+1
      nps=max*(max+1)/2
      i2=mm+nps*ii-kmj
      ipos2=i2
      jj=nod+nod
      ipos3=i2+kmj+(nod+max)*(kmj+1)-jj*kmj
_ELSE
694   ot(mm)=jhog(jhg)
      if(nod.gt.1) go to 695
      if(nod.eq.0) go to 797
      i2=mm-kmj
      jj=nod+nod
      i3=mm+kmj+1-jj*kmj
      if(debug) write(6,99988) ic, i2, i3
      call pack12
      go to 735
797   call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      jrec=jrec+nsz
      go to 551
695   ii=ii+kmj
      ic=mm-ii+1
      nps=max*(max+1)/2
      i2=mm+nps*ii-kmj
      jj=nod+nod
      i3=i2+kmj+(nod+max)*(kmj+1)-jj*kmj
      if(debug) write(6,99988) ic, i2, i3
      call pack12
_ENDIF
      mix=0
      do 698 k=1,max
      mix=mix+2
      imike(1)=imox(mix-1)
      imike(2)=imox(mix)
      ica=ipar(k)
      nix=0
      do 698 l=1,k
      nix=nix+2
      jjog(1)=imox(nix-1)
      jjog(2)=imox(nix)
      jca=ipar(l)
      if(ica.gt.0) go to 699
      jca=-jca
699   mm=mm+1
_IF(notused)
      if(jca.lt.0) jca=jca+256
      icar=jca
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 700
      mm=0
      jrec=jrec+1
      call dwritet(cart)
_ELSE
      ot(mm)=jca
      if(mm.lt.ibb) go to 700
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
_ENDIF
700   jhg=ihg
      do 701 m=1,kmj
      mm=mm+1
      jhg=jhg+1
      khg=jhog(jhg)
      nd=khg*nod+kr
      nps=0
      ipr=1
      ico=imike(1)
      jpr=jjog(1)
      do 702 n=1,nod
      nd=nd+1
      if(n.eq.ico) go to 703
      nps=nps+1
705   if(nps.ne.jpr) go to 702
      nps=nps+1
      if(ipr.eq.2) go to 702
      ipr=2
      jpr=jjog(2)
      go to 705
702   ijog(nps)=jmat(nd)
703   ico=imike(2)
      itim(1)=jmat(nd)
      jca=n+1
      do 706 n=jca,nod
      nd=nd+1
      if(n.eq.ico) go to 707
      nps=nps+1
708   if(nps.ne.jpr) go to 706
      nps=nps+1
      if(ipr.eq.2) go to 706
      ipr=2
      jpr=jjog(2)
      go to 708
706   ijog(nps)=jmat(nd)
 707  itim(2)=jmat(nd)
      if (n.eq.nod) go to 710
      jca=n+1
      do 711 n=jca,nod
      nd=nd+1
      nps=nps+1
712   if(nps.ne.jpr) go to 711
      nps=nps+1
      if(ipr.eq.2) go to 711
      ipr=2
      jpr=jjog(2)
      go to 712
711   ijog(nps)=jmat(nd)
710   if(itim(1).lt.0) go to 713
      if(itim(2).lt.0) go to 714
_IF(notused)
      icar=1
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 715
      mm=0
      jrec=jrec+1
      call dwritet(cart)
_ELSE
      ot(mm)=1
      if(mm.lt.ibb) go to 715
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
_ENDIF
715   ipr=jjog(1)
      ijog(ipr)=1
      ipr=jjog(2)
      ijog(ipr)=1
722   nps=kr
      do 716 n=1,ndet
      nps=nps+nod
      ipr=nps
      do 717 if=1,nod
      ipr=ipr+1
      if(ijog(if).ne.jmat(ipr)) go to 716
717   continue
      go to 718
716   continue
      go to 7781
718   mm=mm+1
_IF(notused)
      icar=n
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 719
      mm=0
      jrec=jrec+1
      call dwritet(cart)
719   mm=mm+1
      if(mm.lt.ibb) go to 701
      mm=0
      jrec=jrec+1
      call dwritet(cart)
      go to 701
713   if(itim(2).gt.0) go to 720
      icar=1
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 721
      mm=0
      jrec=jrec+1
      call dwritet(cart)
721   ipr=jjog(1)
      ijog(ipr)=-1
      ipr=jjog(2)
      ijog(ipr)=-1
      go to 722
714      icar=2
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 723
      mm=0
      jrec=jrec+1
      call dwritet(cart)
_ELSE
      ot(mm)=n
      if(mm.lt.ibb) go to 719
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
719   mm=mm+1
      if(mm.lt.ibb) go to 701
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
      go to 701
713   if(itim(2).gt.0) go to 720
      ot(mm)=1
      if(mm.lt.ibb) go to 721
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
721   ipr=jjog(1)
      ijog(ipr)=-1
      ipr=jjog(2)
      ijog(ipr)=-1
      go to 722
714   ot(mm)=2
      if(mm.lt.ibb) go to 723
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
_ENDIF
723   ipr=jjog(1)
      ijog(ipr)=-1
      jco=1
      jpr=jjog(2)
      ijog(jpr)=1
728   nps=kr
      do 724 n=1,ndet
      nps=nps+nod
      kpr=nps
      do 725 if=1,nod
      kpr=kpr+1
      if(ijog(if).ne.jmat(kpr)) go to 724
725   continue
      go to 726
724   continue
      go to 7781
726   mm=mm+1
_IF(notused)
      icar=n
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 727
      mm=0
      jrec=jrec+1
      call dwritet(cart)
_ELSE
      ot(mm)=n
      if(mm.lt.ibb) go to 727
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
_ENDIF
727   if(jco.lt.0) go to 701
      ijog(ipr) =1
      ijog(jpr)=-1
      jco=-1
      go to 728
_IF(notused)
720      icar=2
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 729
      mm=0
      jrec=jrec+1
      call dwritet(cart)
_ELSE
720   ot(mm)=2
      if(mm.lt.ibb) go to 729
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
_ENDIF
729   ipr=jjog(1)
      ijog(ipr)=1
      jpr=jjog(2)
      ijog(jpr)=-1
      jco=1
730   nps=kr
      do 731 n=1,ndet
      nps=nps+nod
      kpr=nps
      do 732 if=1,nod
      kpr=kpr+1
      if(ijog(if).ne.jmat(kpr)) go to 731
732   continue
      go to 733
731   continue
      go to 7781
733   mm=mm+1
_IF(notused)
      icar=n
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 734
      mm=0
      jrec=jrec+1
      call dwritet(cart)
_ELSE
      ot(mm)=n
      if(mm.lt.ibb) go to 734
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
_ENDIF
734   if(jco.lt.0) go to 701
      ijog(ipr)=-1
      ijog(jpr)=1
      jco=-1
      go to 730
701   continue
698   continue
735   do 736 k=1,nod
      ica=k
      do 736 l=1,k
      mm=mm+1
      ica=ica+1
      jca=-1
      if(ica-2*(ica/2).gt.0) jca=1
_IF(notused)
      if(jca.lt.0) jca=jca+256
      icar=jca
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 737
      mm=0
      jrec=jrec+1
      call dwritet(cart)
_ELSE
      ot(mm)=jca
      if(mm.lt.ibb) go to 737
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
_ENDIF
737   if(j.eq.1) go to 746
      jhg=ihg
      do 738 m=1,kmj
      mm=mm+1
      jhg=jhg+1
      khg=jhog(jhg)
      nd=khg*nod+kr
      itim(l)=jmat(nd+k)
      nps=0
      do 740 n=1,nod
      nd=nd+1
      if(n.eq.k) go to 740
      nps=nps+1
      if(nps.ne.l) go to 741
      nps=nps+1
741   itim(nps)=jmat(nd)
740   continue
      nps=kr
      do 744 n=1,ndet
      nps=nps+nod
      ipr=nps
      do 745 if=1,nod
      ipr=ipr+1
      if(itim(if).ne.jmat(ipr)) go to 744
745   continue
      go to 747
744   continue
      go to 7781
_IF(notused)
747   icar=n
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 738
      mm=0
      jrec=jrec+1
      call dwritet(cart)
  738 continue
      go to 736
  746 do 748 m=1,kmj
      mm=mm+1
      icar=1
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if (mm.lt.ibb) go to 748
      mm=0
      jrec=jrec+1
      call dwritet (cart)
_ELSE
747   ot(mm)=n
      if(mm.lt.ibb) go to 738
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
  738 continue
      go to 736
  746 do 748 m=1,kmj
      mm=mm+1
      ot(mm)=1
      if (mm.lt.ibb) go to 748
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
_ENDIF
  748 continue
  736 continue
      do 749 k=1,nod
      ica=k
      do 749 l=1,k
      mm=mm+1
      ica=ica+1
      jca=1
      if (ica-2*(ica/2).gt.0) jca=-1
_IF(notused)
      if(jca.lt.0) jca=jca+256
      icar=jca
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if (mm.lt.ibb) go to 750
      mm=0
      jrec=jrec+1
      call dwritet(cart)
_ELSE
      ot(mm)=jca
      if (mm.lt.ibb) go to 750
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
_ENDIF
  750 if (j.eq.1) go to 752
      jhg=ihg
      do 751 m=1,kmj
      mm=mm+1
      jhg=jhg+1
      khg=jhog(jhg)
      nd=khg*nod+kr
      ipr=jmat(nd+k)
      itim(l)=ipr
      nps=0
      do 753 n=1,nod
      nd=nd+1
      if (n.eq.k) go to 753
      nps=nps+1
      if (nps.ne.l) go to 754
      nps=nps+1
  754 itim(nps)=jmat(nd)
  753 continue
      nps=kr
      do 755 n=1,ndet
      nps=nps+nod
      jpr=nps
      do 756 if=1,nod
      jpr=jpr+1
      if (itim(if).ne.jmat(jpr)) go to 755
  756 continue
      go to 757
  755 continue
      go to 7781
  757 jpr=1
      if (ipr.lt.0) jpr=2
_IF(notused)
      if(jpr.lt.0) jpr=jpr+256
      icar=jpr
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if (mm.lt.ibb) go to 758
      mm=0
      jrec=jrec+1
      call dwritet(cart)
 758  mm=mm+1
      icar=n
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if (mm.lt.ibb) go to 759
      mm=0
      jrec=jrec+1
      call dwritet(cart)
_ELSE
      ot(mm)=jpr
      if (mm.lt.ibb) go to 758
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
 758  mm=mm+1
      ot(mm)=n
      if (mm.lt.ibb) go to 759
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
_ENDIF
  759 if (jpr.eq.2) go to 760
      nps=0
      do 761 n=1,nod
      if (n.eq.l) go to 761
      nps=nps+1
      if (itim(n).lt.0) go to 761
      mm=mm+1
      if (mm.lt.ibb) go to 762
_IF(notused)
      mm=0
      jrec=jrec+1
      call dwritet(cart)
  762 mm=mm+1
      icar=nps
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if (mm.lt.ibb) go to 761
      mm=0
      jrec=jrec+1
      call dwritet(cart)
_ELSE
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
  762 mm=mm+1
      ot(mm)=nps
      if (mm.lt.ibb) go to 761
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
_ENDIF
  761 continue
      nps=0
      do 764 n=1,nod
      if (n.eq.l) go to 764
      nps=nps+1
      if (itim(n).gt.0) go to 764
      do 765 if=1,nod
  765 ijog(if)=itim(if)
      ijog(l)=-1
      ijog(n)=1
      jca=kr
      do 766 if=1,ndet
      jca=jca+nod
      mal=jca
      do 767 nal=1,nod
      mal=mal+1
      if (ijog(nal).ne.jmat(mal)) go to 766
  767 continue
      go to 768
  766 continue
      go to 7781
  768 mm=mm+1
_IF(notused)
      icar=if
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if (mm.lt.ibb) go to 769
      mm=0
      jrec=jrec+1
      call dwritet(cart)
  769 mm=mm+1
      icar=nps
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if (mm.lt.ibb) go to 764
      mm=0
      jrec=jrec+1
      call dwritet(cart)
_ELSE
      ot(mm)=if
      if (mm.lt.ibb) go to 769
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
  769 mm=mm+1
      ot(mm)=nps
      if (mm.lt.ibb) go to 764
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
_ENDIF
  764 continue
      go to 751
  760 if (j.eq.2) go to 784
      nps=0
      do 775 n=1,nod
      if (n.eq.l) go to 775
      nps=nps+1
      if (itim(n).gt.0) go to 775
      mm=mm+1
      if (mm.lt.ibb) go to 776
_IF(notused)
      mm=0
      jrec=jrec+1
      call dwritet(cart)
776   mm=mm+1
      icar=nps
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 775
      mm=0
      jrec=jrec+1
      call dwritet(cart)
_ELSE
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
776   mm=mm+1
      ot(mm)=nps
      if(mm.lt.ibb) go to 775
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
_ENDIF
775   continue
784   nps=0
      do 777 n=1,nod
      if(n.eq.l) go to 777
      nps=nps+1
      if(itim(n).lt.0) go to 777
      do 778 if=1,nod
778   ijog(if)=itim(if)
      ijog(l)=1
      ijog(n)=-1
      jca=kr
      do 779 if=1,ndet
      jca=jca+nod
      mal=jca
      do 780 nal=1,nod
      mal=mal+1
      if(ijog(nal).ne.jmat(mal)) go to 779
780   continue
      go to 781
779   continue
      go to 7781
781   mm=mm+1
_IF(notused)
      icar=if
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 782
      mm=0
      jrec=jrec+1
      call dwritet(cart)
782   mm=mm+1
      icar=nps
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if(mm.lt.ibb) go to 777
      mm=0
      jrec=jrec+1
      call dwritet(cart)
_ELSE
      ot(mm)=if
      if(mm.lt.ibb) go to 782
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
782   mm=mm+1
      ot(mm)=nps
      if(mm.lt.ibb) go to 777
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
_ENDIF
777   continue
751   continue
      go to 749
  752 do 790 m=1,kmj
      mm=mm+1
_IF(notused)
      icar=1
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if (mm.lt.ibb) go to 771
      mm=0
      jrec=jrec+1
      call dwritet(cart)
  771 mm=mm+1
      icar=1
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if (mm.lt.ibb) go to 772
      mm=0
      jrec=jrec+1
      call dwritet(cart)
_ELSE
      ot(mm)=1
      if (mm.lt.ibb) go to 771
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
  771 mm=mm+1
      ot(mm)=1
      if (mm.lt.ibb) go to 772
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
_ENDIF
  772 if (nod.eq.1) go to 790
      nps=0
      do 770 n=1,nod
      if (n.eq.k) go to 770
      nps=nps+1
      mm=mm+1
      if (mm.lt.ibb) go to 774
_IF(notused)
      mm=0
      jrec=jrec+1
      call dwritet(cart)
  774 mm=mm+1
      icar=nps
_IF(littleendian)
      cart(mm)=car(1)
_ELSE
      cart(mm)=car(4)
_ENDIF
      if (mm.lt.ibb) go to 770
      mm=0
      jrec=jrec+1
      call dwritet(cart)
_ELSE
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
  774 mm=mm+1
      ot(mm)=nps
      if (mm.lt.ibb) go to 770
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      mm=0
      jrec=jrec+nsz
_ENDIF
  770 continue
  790 continue
749   continue
      if(mm.eq.0) go to 551
_IF(notused)
      call dwritet(cart)
      jrec=jrec+1
_ELSE
      call stopbk3
      call pack(lout,8,ot,ibb)
      call sttout3(jrec)
      jrec=jrec+nsz
_ENDIF
551   continue
2     write(iwr,783)nmul,irec,jrec
820   format(/10x,10i10/)
_IF(notused)
      call dwclos
_ELSE
      call stopbk3
c
      call clredx
_ENDIF
783   format(/5x,'end case ',3i10/)
      write(iwr,22)cpulft(1)
   22 format(/1x,'***  end of table generation at ',f10.2,' seconds'/)
      return
950   write(iwr,952) kmj,nod
952   format(5x,'kmj too big',2i6)
      go to 7781
951   write(iwr,953) kmj,ifac,nod
953   format(5x,'kmj came out wrong',3i6)
      go to 7781
954   write(iwr,957) ndet,nod
957   format(5x,'ndet too big',2i6)
      go to 7781
955   write(iwr,958) ib,ndet,nmns
958   format(5x,'error in imat',3i6)
7781  call caserr('invalid parameters detected in table generator')
      return
      end
      subroutine dgelgt(r,a,m,n,eps,ier)
      implicit REAL (a-h,o-z)
c
      dimension a(*),r(*)
      if(m)23,23,1
c
c     search for greatest element in matrix a
    1 ier=0
      piv=0.d0
      mm=m*m
      nm=n*m
      do 3 l=1,mm
      tb= abs(a(l))
      if(tb-piv)3,3,2
    2 piv=tb
      i=l
    3 continue
      tol=eps*piv
c     a(i) is pivot element. piv contains the absolute value of a(i).
c
c
c     start elimination loop
      lst=1
      do 17 k=1,m
c
c     test on singularity
      if(piv)23,23,4
    4 if(ier)7,5,7
    5 if(piv-tol)6,6,7
    6 ier=k-1
    7 pivi=1.d0/a(i)
      j=(i-1)/m
      i=i-j*m-k
      j=j+1-k
c     i+k is row-index, j+k column-index of pivot element
c
c     pivot row reduction and row interchange in right hand side r
      do 8 l=k,nm,m
      ll=l+i
      tb=pivi*r(ll)
      r(ll)=r(l)
    8 r(l)=tb
c
c     is elimination terminated
      if(k-m)9,18,18
c
c     column interchange in matrix a
    9 lend=lst+m-k
      if(j)12,12,10
   10 ii=j*m
      do 11 l=lst,lend
      tb=a(l)
      ll=l+ii
      a(l)=a(ll)
   11 a(ll)=tb
c
c     row interchange and pivot row reduction in matrix a
   12 do 13 l=lst,mm,m
      ll=l+i
      tb=pivi*a(ll)
      a(ll)=a(l)
   13 a(l)=tb
c
c     save column interchange information
      a(lst)=j
c
c     element reduction and next pivot search
      piv=0.d0
      lst=lst+1
      j=0
      do 16 ii=lst,lend
      pivi=-a(ii)
      ist=ii+m
      j=j+1
      do 15 l=ist,mm,m
      ll=l-j
      a(l)=a(l)+pivi*a(ll)
      tb= abs(a(l))
      if(tb-piv)15,15,14
   14 piv=tb
      i=l
   15 continue
      do 16 l=k,nm,m
      ll=l+j
   16 r(ll)=r(ll)+pivi*r(l)
   17 lst=lst+m
c     end of elimination loop
c
c
c     back substitution and back interchange
   18 if(m-1)23,22,19
   19 ist=mm+m
      lst=m+1
      do 21 i=2,m
      ii=lst-i
      ist=ist-lst
      l=ist-m
      l=a(l)+.5d0
      do 21 j=ii,nm,m
      tb=r(j)
      ll=j
      do 20 k=ist,mm,m
      ll=ll+1
   20 tb=tb-a(k)*r(ll)
      k=j+l
      r(j)=r(k)
   21 r(k)=tb
   22 return
c
c
c     error return
   23 ier=-1
      return
      end
      subroutine dmfgrt(a,m,n,eps,irank,irow,icol)
      implicit REAL (a-h,o-z)
c
      dimension a(*),irow(*),icol(*)
c
c        test of specified dimensions
      if(m)2,2,1
    1 if(n)2,2,4
    2 irank=-1
    3 return
c        return in case of formal errors
c
c
c        initialize column index vector
c        search first pivot element
    4 irank=0
      piv=0.d0
      jj=0
      do 6 j=1,n
      icol(j)=j
      do 6 i=1,m
      jj=jj+1
      hold=a(jj)
      if( abs(piv)- abs(hold))5,6,6
    5 piv=hold
      ir=i
      ic=j
    6 continue
c
c        initialize row index vector
      do 7 i=1,m
    7 irow(i)=i
c
c        set up internal tolerance
      tol=abs(eps*     piv )
c
c        initialize elimination loop
      nm=n*m
      do 190 ncol=m,nm,m
c
c        test for feasibility of pivot element
      if(abs(piv)-tol)20,20,9
c
c        update rank
    9 irank=irank+1
c
c        interchange rows if necessary
      jj=ir-irank
      if(jj)12,12,10
   10 do 11 j=irank,nm,m
      i=j+jj
      save=a(j)
      a(j)=a(i)
   11 a(i)=save
c
c        update row index vector
      jj=irow(ir)
      irow(ir)=irow(irank)
      irow(irank)=jj
c
c        interchange columns if necessary
   12 jj=(ic-irank)*m
      if(jj)15,15,13
   13 kk=ncol
      do 14 j=1,m
      i=kk+jj
      save=a(kk)
      a(kk)=a(i)
      kk=kk-1
   14 a(i)=save
c
c        update column index vector
      jj=icol(ic)
      icol(ic)=icol(irank)
      icol(irank)=jj
   15 kk=irank+1
      mm=irank-m
      ll=ncol+mm
c
c        test for last row
      if(mm)16,25,25
c
c        transform current submatrix and search next pivot
   16 jj=ll
      save=piv
      piv=0.d0
      do 191 j=kk,m
      jj=jj+1
      hold=a(jj)/save
      a(jj)=hold
      l=j-irank
c
c        test for last column
      if(irank-n) 17,191,191
   17 ii=jj
      do 19 i=kk,n
      ii=ii+m
      mm=ii-l
      a(ii)=a(ii)-hold*a(mm)
      if( abs(a(ii))- abs(piv))19,19,18
   18 piv=a(ii)
      ir=j
      ic=i
   19 continue
  191 continue
  190 continue
c
c        set up matrix expressing row dependencies
   20 if(irank-1)3,25,21
   21 ir=ll
      do 24 j=2,irank
      ii=j-1
      ir=ir-m
      jj=ll
      do 23 i=kk,m
      hold=0.d0
      jj=jj+1
      mm=jj
      ic=ir
      do 22 l=1,ii
      hold=hold+a(mm)*a(ic)
      ic=ic-1
   22 mm=mm-m
   23 a(mm)=a(mm)-hold
   24 continue
c
c        test for column regularity
   25 if(n-irank)3,3,26
c
c        set up matrix expressing basic variables in terms of free
c        parameters (homogeneous solution).
   26 ir=ll
      kk=ll+m
      do 30 j=1,irank
      do 29 i=kk,nm,m
      jj=ir
      ll=i
      hold=0.d0
      ii=j
   27 ii=ii-1
      if(ii)29,29,28
   28 hold=hold-a(jj)*a(ll)
      jj=jj-m
      ll=ll-1
      goto 27
   29 a(ll)=(hold-a(ll))/a(jj)
   30 ir=ir-1
      return
      end
_IF(notused)
**************************************************************************      
      subroutine dwritet(a)                                                      
      integer a
      dimension a(nwords)                                                         
      common /directio/ lu,nbytes,nwords,nrec,irc,iopt                            
INCLUDE(common/iofile)
      nrec=nrec+1                                                               
      if (iopt.eq.1) then                                                       
      write (iwr,*) 'dwrite nrec=',nrec,' nbytes=',nbytes,
     +              ' nwords=', nwords                                                                    
      endif                                                                     
      write (lu,iostat=irc,rec=nrec) a                                          
c@    write (iwr,101) nrec,lu,nbytes                                              
c@101 format(' record no',i10,' written to lu',i2,' - ',i10,' bytes')           
c     write (lu,iostat=irc) a                                                   
      if (irc.ne.0) then                                                        
       write (iwr,100) lu,nrec,nbytes,irc                                         
100    format(' ###### dwrite error on logical unit ',i10/,                     
     *        '        record number is             ',i10/,                     
     *        '        record length (bytes)        ',i10/,                     
     *        '        returncode is                ',z8/)                      
       call caserr('dwrite error')
      endif                                                                     
      return                                                                    
      end                                                                       
c---------------------------------------------------------------------
c     fortran77 dreadt/dwrite routine - bernd hess - 22.5.87
c      2.6.89 - ibm version: status='delete'
c     25.5.88 - drdel
c---------------------------------------------------------------------
c  delcx        convex version
c  delc@        debug
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     drclos only (for hand-linked new mrd-ci package)  (m.c. and g.j.)
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine dwclos
      common /directio/ lu,nbytes,nwords,nrec,irc,iopt
INCLUDE(common/iofile)
      if (iopt.eq.1) then
      write (iwr,*) 'drclos'
      endif
      close (lu)
      return
      end
c---------------------------------------------------------------------
c     fortran77 dreadt/dwrite routine - bernd hess - 22.5.87
c      2.6.89 - ibm version: status='delete'
c     25.5.88 - drdel
c---------------------------------------------------------------------
c  delcx        convex version
c  delc@        debug
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c   dreadt only  (for hand-linked new mrd-ci package)  (g.j. and m.c.)
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine dreadt(a,nr)
      integer a(nwords)
      common /directio/ lu,nbytes,nwords,nrec,irc,iopt
INCLUDE(common/iofile)
      if (iopt.eq.1) then
      write (iwr,*) 'dreadt nrec=',nr,' nbytes=',nbytes,
     *' nwords=',nwords
      endif
      read (lu,rec=nr,iostat=irc) a
      if (irc.ne.0) then
       write (iwr,100) lu,nr,nrec,nbytes,irc
100    format(' ######  dreadt error on logical unit ',i10/,
     *        '        record number is             ',i10/,
     *        '        written by previous dwrite   ',i10/,
     *        '        record length (bytes)        ',i10/,
     *        '        returncode is                ',z8/)
        call caserr('dreadt error in table-ci')
      endif
      return
      end
c---------------------------------------------------------------------
c     fortran77 dreadt/dwrite routine - bernd hess - 22.5.87
c      2.6.89 - ibm version: status='delete'
c     25.5.88 - drdel
c---------------------------------------------------------------------
c  delcx        convex version
c  delc@        debug
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c  dwopen only (for hand-linked new mrd-ci package)  (g.j. and m.c.)
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine dwopen(lun,nbyt)
c
      common /directio/ lu,nbytes,nwords,nrec,irc,iopt
INCLUDE(common/iofile)
c
c     machine dependent parameter: word length
c
      logical open
      data nwlen /4/
c     iopt=1
      if (iopt.eq.1) then
      write (iwr,*) 'dwopen lu=',lun,' nbyt=',nbyt
      endif
      nbytes=(nbyt/nwlen)*nwlen
      if (nbytes.ne.nbyt) then
       write (iwr,101) lun,nbyt,nwlen
101    format(' ###### dwopen error on logical unit ',i10/,
     *        ' record length (in bytes)          : ',i10/,
     *        ' must be divisible by word length  : ',i10)
       call caserr('dwopen error number  1')
      endif
      lu=lun
      nwords=nbyt/nwlen
      inquire (lun,opened=open)
cbe   if (.not.open) then
      open (lu,iostat=irc,access='direct',form='unformatted',
     +      status='unknown',file='table-ci',recl=nbyt)
      if (irc.ne.0) then
       write (iwr,100) lu,nbytes,irc
100    format(' ###### dwopen error on logical unit ',i10/,
     *        '        record length (bytes)        ',i10/,
     *        '        returncode is                ',z8/)
       call caserr('dwopen error number 2')
      endif
c     endif
c@    write (iwr,*) ' dwopen lu',lu,' lrecl=',nbytes
      return
      end
_ENDIF
      subroutine ver_newmrd2(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/newmrd2.m,v $
     +     "/
      data revision /"$Revision: 5774 $"/
      data date /"$Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
