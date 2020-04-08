c ******************************************************
c ******************************************************
c             =   Table-Ci (select module) =
c ******************************************************
c ******************************************************
_EXTRACT(smrd0,mips4,sun)
      subroutine smrd0
      implicit REAL (a-h,o-z), integer(i-n)
      character *8 zcomm,textm,single,double
      character *1 dash
INCLUDE(common/sizes)
      parameter (idim=10970)
c     dimension of array lkon was increased from 2970 to this value
c     this array is in common scrtch and equivalence mkon
c     also b in common/scrtch/ should be changed from 15876 to lower
      common/junkc/zcomm(29),textm
      common/craypk/icall(7),jkonp(mxcsf*mxnshl),npkon(mxcsf)
c  1920 = 24*80 = nshl*80
INCLUDE(common/iofile)
       common/ftape/
     *mtape,kclab,ntype,nf4,mdisk,ltype,kfile,nf22,
     *ntape,nf32,linf,ideli,nf35(8)
      common /tap/ical,m,nko,mxex,nmul,ispace,nprin
      common /scrtch/
     + f(2304),nconf(5),jab(36),isym(8),mj(8),lj(8),iab(maxorb),
     + kj(8),lkon(idim),jcon(maxorb),icuf(12),iduf(15),ideks(9),
     + kcon(mxcrec),nytl(5),ntil(8),nbal(9),isw(10),nj(8),
     + ijoe(mxcsf),ndub(mxcsf),nop(mxcsf),nyzl(mxcsf),
     + ilee(9),inver(8),jkon(mxcsf*mxnshl)
      common/b/ig,ibl,nbxa,nbxb,nbox,jpaw,klx,nshl,kly,ispin,
     + idum1(86+3*mxcsf)
INCLUDE(common/prints)
      dimension mkon(mxcrec),itag(mxcsf)
      equivalence (lkon(1),mkon(1))
      data single,double,dash/
     *' singly',' doubly','-'/
c     nshl=24
      nshl=mxnshl
c     constant maximal number of shells = space for 1 main in field jkon
      nops=5
      ibl=mxcrec
      iswh=5
c     lko=70
      lko=mxcsf
      mmax=40
      call rewftn(ntape)
      call rewftn(mdisk)
      read (ntape) n,nod,jdum,kdum,nbox,mj,kj,lj,nj,nsym,ntil,nbal,
     +             isym,jab,idum
      nbxa=nbox-1
      nbxb=nbox-2
      nnx=(m+nmul-3)/2
      do 87 i=1,iswh
      nnx=nnx+1
 87   nytl(i)=nnx
c     nytl=number of numbers defining configuration in array jkon
      if(.not.oprint(29)) write(iwr,8498)nbox,m,mxex,textm,
     +                    ispace,nko
 8498 format(/
     *' no. of active orbitals      ',i4/
     *' no. of active electrons     ',i4/
     *' excitation level requested  ',i4/
     *' state spin multiplicity  ',a7/
     *' state spatial symmetry      ',i4/
     *' no. of root configurations  ',i4/)
      if(nprin.ne.0)write(iwr,8600)
 8600 format(/
     *' ** print configuration data and symmetry coefficients'
     *//)
      if(nbox.le.0.or.nbox.ge.maxorb) go to 787
      if(mxex.gt.4.or.mxex.lt.0) go to 767
      if(nmul.lt.1.or.nmul.gt.7) go to 783
      do 88 i=1,10
88     isw(i)=0
      if(nmul.gt.1) isw(nmul-1)=1
       isw(nmul+1)=2
       isw(nmul+3)=3
      if(m.lt.2.or.m.gt.mmax) go to 779
c     bad num (too many) of electrons
      if(ispace.lt.1.or.ispace.gt.8) go to 785
c     bad number of irrep
      if(nko.gt.lko) go to 766
c     too many mains
      do 95 i=1,nsym
 95   inver(i)=0
      do 96 i=1,n
      ix=isym(i)
      if (lj(i).eq.0) go to 96
      inver(ix)=i
 96   continue
      ix=0
      do 89 i=1,n
      ilee(i)=ix
      l=lj(i)
      if (l.eq.0) go to 89
      lx=isym(i)
      do 90 j=1,l
      ix=ix+1
   90 iab(ix)=lx
   89 continue
      ilee(n+1)=nbox
      im=0
      do 2 i=1,nko
      nt=im+1
      im=im+nshl
c     where is this main in array with all mains
      do 2000 j=nt,im
2000  jkon(j)=jkonp(j)
      np=npkon(i)
      if (np.gt.nops) go to 775
c     too many open shells
      nyt=1
      if(np.eq.0) go to  94
      nyt=isw(np)
 94   nyt=nytl(nyt)
      nyzl(i)=nyt
      if(i.eq.1) go to 503
      i1=i-1
      mx=-nshl
      do 504 j=1,i1
      mx=mx+nshl
      if(nop(j).ne.np) go to 504
      nv=mx
      la=nt-1
      do 505 k=1,nyt
      la=la+1
      nv=nv+1
      if(jkon(la).ne.jkon(nv)) go to 504
 505  continue
      go to 777
c     at least two ident. mains
 504  continue
 503  na=nmul+np
      if(na-2*(na/2).eq.0) go to 771
      if(np.gt.m) go to 771
      if(np.gt.nmul+3) go to 771
      nop(i)=np
      ndub(i)=(m-np)/2
      if(np.lt.2) go to 500
      nz=jkon(nt)
      if(nz.lt.1.or.nz.gt.nbox) go to 773
c     orbital numbering weired in mains
      do 501 j=2,np
      nt=nt+1
      nv=jkon(nt)
      if(nv.le.nz) go to 769
      if(nv.gt.nbox) go to 773
501   nz=nv
      if(np.eq.m) go to 2
518   mx=np+im-nshl+1
      nv=jkon(mx)
      if(nv.lt.1.or.nv.gt.nbox) go to 773
      lb=im-nshl
      do 506 j=1,np
      lb=lb+1
      if(jkon(lb)-nv) 506,769,507
506   continue
507   jm=np+2
      if(jm.gt.nyt) go to 2
      kp=mx
      do 508 j=jm,nyt
      kg=mx
      mx=mx+1
      nv=jkon(mx)
      if(nv.lt.1.or.nv.gt.nbox) go to 773
      do 509 k=kp,kg
      nz=jkon(k)
      if(nv.le.nz) go to 769
 509   continue
      kg=im-nshl
      do 510 k=1,np
      kg=kg+1
      if(jkon(kg)-nv) 510,769,508
 510   continue
 508   continue
       go to 2
 500   if(np.eq.0) go to 511
       nv=jkon(nt)
       if(nv.lt.1.or.nv.gt.nbox) go to 773
       go to 518
 511   nz=jkon(nt)
       if(nz.lt.1.or.nz.gt.nbox) go to 773
       if(nyt.eq.1) go to 2
       kp=nt
       do 512 j=2,nyt
       kg=nt
       nt=nt+1
       nv=jkon(nt)
       if(nv.lt.1.or.nv.gt.nbox) go to 773
       do 512 k=kp,kg
       nz=jkon(k)
       if(nv.le.nz) go to 769
 512   continue
   2   continue
       kg=nbox*nko
       if (kg .gt. idim)  write(iwr,9996) kg,idim
 9996  format(//1x,'WARNING!'/
     * 1x,'index of array LKON will be out of range!'/
     * 1x,'max index=num.of.orbit*num.of.mains=',i8,' dimension=',i8)
c
       do 513 i=1,kg
 513   lkon(i)=0
c      lkon must be greater or equal then lko*nbox
       kg=-nshl
       kb=-nbox
       do 514 i=1,nko
       kg=kg+nshl
       kb=kb+nbox
       np=nop(i)
       nyt=nyzl(i)
       if(np.eq.0) go to 516
       do 515 j=1,np
       nv=kb+jkon(j+kg)
 515   lkon(nv)=1
       if(np.eq.m) go to 514
 516   np=np+1
       do 517  j=np,nyt
       nv=kb+jkon(j+kg)
 517   lkon(nv)=2
514    continue
      if(oprint(33)) write(iwr,9000)
 9000 format(40x,19('=')/
     *40x,'root configurations'/40x,19('=')/)
      kg=-nshl
      do 521 i=1,nko
      kg=kg+nshl
      np=nop(i)
      nv=kg+nyzl(i)
      if(oprint(33)) write(iwr,9001)i
 9001 format(/
     *' reference configuration no.',i3)
      if(oprint(33)) write(iwr,9002)(dash,j=1,30)
 9002 format(1x,129a1)
      if(np.eq.0)go to 9003
      if(oprint(33)) write(iwr,9004)single
 9004 format(/a7,' occupied orbitals')
      do 9005 k=1,nsym
      ik=0
      do 9006 j=1,np
      l=jkon(j+kg)
      if(iab(l).ne.k)go to 9006
      ik=ik+1
      itag(ik)=l
 9006 continue
      if(ik.eq.0)go to 9005
      if(oprint(33)) 
     + write(iwr,9007)k,(itag(loop),loop=1,ik)
 9007 format(/1x,' irrep. no.',i1,2x,
     *' orbital sequence nos.',2x,30i4)
 9005 continue
 9003 if(nyzl(i).eq.np)go to 521
      l=nyzl(i)-np
      if(oprint(33)) write(iwr,9004)double
      do 9008 k=1,nsym
      ik=0
      do 9009 j=1,l
      ll=jkon(kg+np+j)
      if(iab(ll).ne.k)go to 9009
      ik=ik+1
      itag(ik)=ll
 9009 continue
      if(ik.eq.0)go to 9008
      if(oprint(33)) 
     + write(iwr,9007)k,(itag(loop),loop=1,ik)
 9008 continue 
 521  continue
c
      if(oprint(33)) write(iwr,9002)(dash,i=1,30)
      if(.not.oprint(29)) write(iwr,8499)
 8499 format(  //' orbital symmetry designations'/1x,29('-')//
     *' irrep.no.',5x,'orbital sequence nos.'/)
      do 9010 i=1,nsym
      ik=0
      do 9011 j=1,nbox
      if(iab(j).ne.i)go to 9011
      ik=ik+1
 9011 continue
 9010 itag(i)=ik
      il=0
      do 9012 i=1,nsym
      j=itag(i)
      if(j.eq.0)go to 9012
      ik=il+j
      il=il+1
      if(.not.oprint(29)) write(iwr,9013)i,il,ik
 9013 format(4x,i2,13x,i4,'  to ',i3)
      il=ik
 9012 continue
      if(.not.oprint(29)) write(iwr,9002)(dash,i=1,129)
c
      cpu=cpulft(1)
      if(.not.oprint(29)) write(iwr,9015)cpu
 9015 format(/
     *' commence configuration generation at ',f8.2,' secs.'/)
      ispin=(nmul-1)/2
_IF1()     fms=ispin
_IF1()     if(m.ne.2*(m/2)) fms=fms+0.5d0
      ideks(1)=0
      do 7 i=1,nsym
7     ideks(i+1)=ideks(i)+i
      kg=-nshl
      do 524 i=1,nko
      np=nop(i)
      kg=kg+nshl
      if(np.gt.0) go to 525
      jg=1
      go to 530
 525  jg=jkon(kg+1)
      jg=iab(jg)
      if(np.eq.1) go to 530
      do 527 j=2,np
      nv=jkon(kg+j)
      nv=iab(nv)
      if(nv.lt.jg) go to 528
      jg=ideks(nv)+jg
      go to 527
 528   jg=ideks(jg)+nv
 527  jg=jab(jg)
 530  continue
      nmns=np/2+1-ispin
      if(nmns.le.0) go to 771
 524  ijoe(i)=jg
      do 8 i=1,iswh
8     nconf(i)=0
      ig=0
      klx=0
      llx=0
22    kly=klx
      klx=klx+1
      if (klx.gt.nko) go to 23
      do 24 i=1,nbox
      llx=llx+1
24    jcon(i)=lkon(llx)
      kar=0
      imel=ijoe(klx)
      nbd=1
      nod=nop(klx)
      ncl=ndub(klx)
      nmns=nod/2+1-ispin
      if (klx.eq.1) go to 26
      kg=-nshl
      do 162 i=1,kly
      mx=0
      kg=kg+nshl
      np=nop(i)
      ja=ndub(i)
      kz=kg
      if (np.eq.0) go to 163
      do 164 j=1,np
      kz=kz+1
      jd=jkon(kz)
      if (jcon(jd).gt.0) go to 164
      mx=mx+1
      if (mx.gt.mxex) go to 162
  164 continue
      if (ja.eq.0) go to 28
  163 do 165 j=1,ja
      kz=kz+1
      jd=jkon(kz)
      jd=2-jcon(jd)
      if (jd.eq.0) go to 165
      mx=mx+jd
      if (mx.gt.mxex) go to 162
  165 continue
      go to  28
 162  continue
  26  if (imel.ne.ispace) go to 28
      if (nmns.le.0) go to 28
      ig=ig+1
      kz=kly*nshl
      kcon(ig)=nmns
      nconf(nmns)=nconf(nmns)+1
      if (ig.lt.ibl) go to 166
      ig=0
      write (mdisk) kcon
 166  if (nod.eq.0) go to 167
      do 168 j=1,nod
      kz=kz+1
      ig=ig+1
      kcon(ig)=jkon(kz)
      if (ig.lt.ibl) go to 168
      ig=0
      write (mdisk) kcon
 168  continue
      if (ncl.eq.0) go to 28
 167  do 169 j=1,ncl
      kz=kz+1
      ig=ig+1
      kcon(ig)=jkon(kz)
      if (ig.lt.ibl) go to 169
      ig=0
      write (mdisk) kcon
 169  continue
  28  if (mxex.eq.0) go to 22
      jpaw = min(imel,ispace) + ideks( max(imel,ispace))
      jpaw=jab(jpaw)
38    kar=kar+1
      if (kar.gt.nbox) go to 39
      if (jcon(kar).eq.0) go to 38
      iq=iab(kar)
      kpaw = min(iq,jpaw) + ideks( max(iq,jpaw))
      kpaw=jab(kpaw)
      kpaw=inver(kpaw)
      if (kpaw.eq.0) go to 38
      kas=ilee(kpaw)
      npox=ilee(kpaw+1)
      jcon(kar)=jcon(kar)-1
30    kas=kas+1
      if (kas.eq.kar) go to 30
      if (kas.gt.npox) go to 31
      if (jcon(kas).eq.2) go to 30
      jcon(kas)=jcon(kas)+1
33    if (klx.eq.1) go to 25
      lly=-nshl
      do 34 i=1,kly
      mx=0
      lly=lly+nshl
      np=nop(i)
      ja=ndub(i)
      llz=lly
      if (np.eq.0) go to 173
      do 172 j=1,np
      llz=llz+1
      jd=jkon(llz)
      if (jcon(jd).gt.0) go to 172
      mx=mx+1
      if (mx.gt.mxex) go to 34
 172  continue
      if (ja.eq.0) go to 36
 173  do 174 j=1,ja
      llz=llz+1
      jd=jkon(llz)
      jd=2-jcon(jd)
      if (jd.eq.0) go to 174
      mx=mx+jd
      if (mx.gt.mxex) go to 34
 174  continue
      go to 36
34    continue
      go to 25
31    jcon(kar)=jcon(kar)+1
      go to 38
 3    jcon(kas)=jcon(kas)-1
      go to 30
39    if(mxex.eq.1) go to 22
      kar=0
      nbd=2
      if (jpaw.ne.1) go to 47
46    kar=kar+1
      if (kar.gt.nbox) go to 47
      if (jcon(kar).lt.2) go to 46
      jcon(kar)=0
      kas=0
41    kas=kas+1
      if(kas.eq.kar) go to 41
      if (kas.gt.nbox) go to 44
      if (jcon(kas).gt. 0) go to 41
      jcon(kas) =2
      go to 33
44    jcon(kar)=2
      go to 46
 42   jcon(kas)=0
      go to 41
 47   nbd=3
      kar=0
53    kar=kar+1
      if (kar.gt.nbox) go to 54
      if (jcon(kar).lt.2) go to 53
      jcon(kar)=0
      las=0
55    las=las+1
      if (las.eq.nbox) go to 56
      if (las.eq.kar) go to 55
      if (jcon(las).eq.2)  go to 55
      iq=iab(las)
      kpaw = min(iq,jpaw) + ideks( max(iq,jpaw))
      kpaw=jab(kpaw)
      kpaw=inver(kpaw)
      if (kpaw.eq.0) go to 55
      kas=ilee(kpaw)
      npox=ilee(kpaw+1)
      if (kas.ge.las) go to 184
      if (las.ge.npox) go to 55
      kas=las
 184  jcon(las)=jcon(las)+1
50    kas=kas+1
      if (kas.eq.kar) go to 50
      if (kas.gt.npox) go to 51
      if (jcon(kas).eq.2) go to 50
      jcon(kas)=jcon(kas)+1
      go to 33
 48   jcon(kas)=jcon(kas)-1
      go to 50
51    jcon(las)=jcon(las)-1
      go to 55
56    jcon(kar)=2
      go to 53
54    nbd=4
      kar=0
62    kar=kar+1
      if (kar.eq.nbox) go to 63
      if (jcon(kar).eq.0) go to 62
      iq=iab(kar)
      kpaw = min(iq,jpaw) + ideks( max(iq,jpaw))
      kpaw=jab(kpaw)
      kpaw=inver(kpaw)
      if (kpaw.eq.0) go to 62
      lar=ilee(kpaw)
      npox=ilee(kpaw+1)
      if (lar.ge.kar) go to 185
      if (kar.ge.npox) go to 62
      lar=kar
 185  jcon(kar)=jcon(kar)-1
64    lar=lar+1
      if (lar.gt.npox) go to 65
      if (jcon(lar).eq.0) go to 64
      jcon(lar)=jcon(lar)-1
      kas=0
59    kas=kas+1
      if (kas.eq.lar) go to 59
      if (kas.eq.kar) go to 59
      if (kas.gt.nbox) go to 60
      if (jcon(kas).gt.0)  go to 59
      jcon(kas)=2
      go to 33
60    jcon(lar)=jcon(lar)+1
      go to 64
65    jcon(kar)=jcon(kar)+1
      go to 62
 57   jcon(kas)=0
      go to 59
63    nbd=5
      kar=0
71    kar=kar+1
      if (kar.eq.nbox) go to 581
      if (jcon(kar).eq.0) go to 71
      jcon(kar)=jcon(kar)-1
      iq=iab(kar)
      kpaw = min(iq,jpaw) + ideks( max(iq,jpaw))
      kpaw=jab(kpaw)
      lar=kar
72    lar=lar+1
      if (lar.gt.nbox) go to 73
      if (jcon(lar).eq.0) go to 72
      jcon(lar)=jcon(lar)-1
      iq=iab(lar)
      lpaw = min(iq,kpaw) + ideks( max(iq,kpaw))
      lpaw=jab(lpaw)
      las=0
74    las=las+1
      if (las.eq.nbox) go to 75
      if (las.eq.kar) go to 74
      if (las.eq.lar) go to 74
      if (jcon(las).eq.2) go to 74
      iq=iab(las)
      mpaw = min(iq,lpaw) + ideks( max(iq,lpaw))
      mpaw=jab(mpaw)
      mpaw=inver(mpaw)
      if (mpaw.eq.0) go to 74
      kas=ilee(mpaw)
      npox=ilee(mpaw+1)
      if (kas.ge.las) go to 192
      if (las.ge.npox) go to 74
      kas=las
 192  jcon(las)=jcon(las)+1
68    kas=kas+1
      if (kas.eq.kar) go to 68
      if (kas.eq.lar) go to 68
      if (kas.gt.npox) go to 69
      if (jcon(kas).eq.2) go to 68
      jcon(kas)=jcon(kas)+1
      go to 33
66    jcon(kas)=jcon(kas)-1
      go to 68
69    jcon(las)=jcon(las)-1
      go to 74
75    jcon(lar)=jcon(lar)+1
      go to 72
73    jcon(kar)=jcon(kar)+1
      go to 71
581   if(mxex.eq.2) go to 22
      call smrd01
      go to 22
36    go to (3,42,48,57,66), nbd
25    nod=0
      ncl=0
      do 76 i=1,nbox
      jc=jcon(i)
      if (jc.eq.0) go to 76
      if (jc.eq.1) go to 77
      ncl=ncl+1
      iduf(ncl)=i
      go to 76
77    nod=nod+1
      icuf(nod)=i
76    continue
      nmns=nod/2+1-ispin
      if(nmns.le.0) goto 36
      ig=ig+1
      kcon(ig)=nmns
      nconf(nmns)=nconf(nmns)+1
      if (ig.lt.ibl) go to 85
      ig=0
      write(mdisk)kcon
 85   if (nod.eq.0) go to 80
      do 86 i=1,nod
      ig=ig+1
      kcon(ig)=icuf(i)
      if (ig.lt.ibl) go to 86
      ig=0
      write(mdisk) kcon
86    continue
      if (ncl.eq.0) go to 36
80    do 81 i=1,ncl
      ig=ig+1
      kcon(ig)=iduf(i)
      if (ig.lt.ibl) go to 81
      ig=0
      write(mdisk) kcon
81    continue
      go to 36
   23 j=nmul-3
      do 9020 i=1,5
      j=j+2
      if(.not.oprint(29)) write(iwr,8601)j,nconf(i)
      if(nconf(i).eq.0)go to 9020
      if(j.gt.m) then
         write(iwr,8603) i
 8603    format(//1x,'DANGER!'//
     *1x,'configurations with more open shells than'/
     *1x,'active electrons occured!'/
     *1x,'nconf(',i2,')=0 done, but scratch file contents wrong')
         nconf(i)=0
      endif
 8601 format(' no. of configurations with ',i2,
     *' open shell(s) = ',i7)
 9020 continue
      write(mdisk) kcon
      call rewftn(mtape)
      write (mtape) nbox,m,nmul,ispace,nconf,nytl,mxex,jkon,nko,iswh,
     anyzl
_IF1()      nod=nmul-3
      nerrors=0
      do 100 iw=1,iswh
      nc=nconf(iw)
_IF1()      nod=nod+2
      if (nc.eq.0) go to 100
      jg=0
      nl=nytl(iw)
_IF1()c      write (iwr,241) nc,nod
      call rewftn(mdisk)
      read(mdisk)kcon
      ig=0
      do 201 i=1,nc
202   ig=ig+1
      nad=kcon(ig)
      if ((nmul-3+2*nad) .gt. m) then
        write(iwr,9999) iw,i,nad
9999    format(/1x,'ERROR IN REORDERING !!!'/
     *  1x,'Supercat= ',i2,' conf= ',i7/
     *  'from mdisk was read nad=',i3)
        nerrors=nerrors+1
        if (nerrors.gt.500) then
          write(iwr,9998) nerrors
9998      format('TOO MANY ERRORS:',i4,' - exiting the program')
          call caserr('selection errors determinded')
        endif
      endif
c read number of supercat
      if (ig.lt.ibl) go to 203
      ig=0
      read(mdisk) kcon
203   if(nad.eq.iw) go to 220
c skip conf from different supercat
      ip=nytl(nad)
      ig=ig+ip
      if(ig.lt.ibl) go to 202
      ig=ig-ibl
      read(mdisk)kcon
      go to 202
c when is from treated supercat
220   lx=ig+1
      ig=ig+nl
      if(ig.gt.ibl) go to 221
      do 222 j=lx,ig
      jg=jg+1
c read this conf in field mkon and write to mtape in blocks
      mkon(jg)=kcon(j)
      if (jg.lt.ibl) go to 222
      write(mtape)mkon
      jg=0
222   continue
      if (ig.lt.ibl) go to 201
      ig=0
      read(mdisk)kcon
      go to 201
221   do 224 j=lx,ibl
      jg=jg+1
      mkon(jg)=kcon(j)
      if (jg.lt.ibl) go to 224
      write(mtape) mkon
      jg=0
224   continue
      ig=ig-ibl
      read(mdisk) kcon
      do 225 j=1,ig
      jg=jg+1
      mkon(jg)=kcon(j)
      if (jg.lt.ibl) go to 225
      write(mtape)mkon
      jg=0
225   continue
201   continue
      write(mtape)mkon
100   continue
c
c
      if(nerrors.gt.0) write(iwr,9995) nerrors
9995  format(/1x,'NUMBER OF ERRORS during reordering: ',i4/)
      cpu=cpulft(1)
      if(.not.oprint(29)) write(iwr,240)cpu
240   format(/' end of configuration generation at',
     * f8.2,' secs.')
      return
766   write(iwr,760) nko,lko
760   format(5x,'too many mains',2i6)
      call caserr('invalid number of main configurations')
  767 write (iwr,768)
  768 format(5x, 'excitation class not allowed')
      call caserr('excitation class not allowed')
  779 write(iwr,780)
 780  format(5x,'requested number of electrons troublesome')
      call caserr('invalid number of electrons')
 777  write(iwr,778)
 778  format(5x,'at least two of the main configurations are identical')
      call caserr('identical main configurations detected')
 775  write(iwr,776)
 776  format(5x,'too many open shells in mains for dimensions')
      call caserr('invalid number of open shells')
 773  write(iwr,774)
 774  format(5x,'orbital numbering is weird in mains')
      call caserr('invalid orbital numbering detected')
 771  write(iwr,772)
 772  format(5x,'open shell structure in mains inconsistent with multipl
     1icity')
      call caserr('inconsistent multiplicity and open shell structure')
 769  write(iwr,770) i,nv,nz,jkon(kg),kg,jkon(lb),lb
 770  format(5x,'pauli was right or maybe permutation error',7i5)
      call caserr('possible permutation error detected')
 783  write(iwr,784)
 784   format(5x,'invalid spin multiplicity')
       call caserr('invalid spin multiplicity')
 785   write(iwr,786)
 786   format(5x,'invalid spatial symmetry specified')
       call caserr('invalid spatial symmetry specified')
 787   write(iwr,788)
 788   format(5x,'invalid number of active orbitals')
       call caserr('invalid number of active orbitals')
       return
       end
_ENDEXTRACT
      subroutine smrd01
      implicit REAL (a-h,o-z), integer(i-n)
INCLUDE(common/sizes)
      parameter (idim=10970)
      common/ftape/
     *mtape,kclab,ntype,nf4,mdisk,ltype,kfile,nf22,
     *ntape,nf32,linf,ideli,
     * nf35(8)
      common /tap/ical,m,nko,mxex,nmul,ispace,nprin
      common /scrtch/
     + f(2304),nconf(5),jab(36),isym(8),mj(8),lj(8),iab(maxorb),
     + kj(8),lkon(idim),jcon(maxorb),icuf(12),iduf(15),ideks(9),
     + kcon(mxcrec),nytl(5),ntil(8),nbal(9),isw(10),nj(8),
     + ijoe(mxcsf),ndub(mxcsf),nop(mxcsf),nyzl(mxcsf),
     + ilee(9),inver(8),jkon(mxcsf*mxnshl)
      common/b/ ig,ibl,nbxa,nbxb,nbox,jpaw,klx,nshl,kly,ispin,
     + idum1(86+3*mxcsf)
      dimension mkon(mxcrec)
      equivalence (lkon(1),mkon(1))
      nbd=1
      kar=0
414   kar=kar + 1
      if(kar.gt.nbox) go to 415
      if(jcon(kar).eq.2)  go to 414
      iq=iab(kar)
      kpaw = min(iq,jpaw) + ideks( max(iq,jpaw))
      kpaw= jab(kpaw)
      kpaw=inver(kpaw)
      if(kpaw.eq.0) go to 414
      kas=ilee(kpaw)
      npox= ilee(kpaw+ 1)
      jcon(kar) =jcon(kar)  + 1
418   kas= kas + 1
      if(kas.gt.npox) go to 419
      if(kas.eq.kar) go to 418
      if(jcon(kas).eq.0) go to 418
      jcon(kas) = jcon(kas) - 1
      lar=0
420   lar=lar + 1
      if(lar.eq.kar.or.lar.eq.kas) go to 420
      if(lar.gt.nbox) go to 421
      if(jcon(lar).lt.2) go to 420
      jcon(lar)=0
      las=0
422   las=las + 1
      if(las.eq.kar.or.las.eq.kas) go to 422
      if(las.eq.lar) go to 422
      if(las.gt.nbox) go to 423
      if(jcon(las).gt.0) go to 422
      jcon(las) = 2
      go to 33
423   jcon(lar) =2
      go to 420
421   jcon(kas) =jcon(kas)  + 1
      go to 418
419   jcon(kar)= jcon(kar) -1
      go to 414
400   jcon(las) = 0
      go to 422
415   nbd=2
      kar=0
424   kar = kar + 1
      if(kar.eq.nbxa) go to 425
      if(jcon(kar).eq.0 ) go to 424
      jcon(kar) = jcon(kar) - 1
      iq= iab(kar)
      kpaw = min(iq,jpaw) + ideks( max(iq,jpaw))
      kpaw= jab(kpaw)
      lar  = kar
428   lar = lar + 1
      if(lar.eq.nbox) go to 429
      if(jcon(lar).eq.0) go to 428
      jcon(lar) = jcon(lar) - 1
      iq=iab(lar)
      lpaw = min(iq,kpaw) + ideks( max(iq,kpaw))
      lpaw= jab(lpaw)
      mar= lar
432   mar= mar + 1
      if(mar.gt.nbox) go to 433
      if(jcon(mar).eq.0) go to 432
      iq=iab(mar)
      mpaw = min(iq,lpaw) + ideks( max(iq,lpaw))
      mpaw=jab(mpaw)
      mpaw=inver(mpaw)
      if(mpaw.eq.0) go to 432
      kas=ilee(mpaw)
      npox=ilee(mpaw + 1)
      jcon(mar) =jcon(mar) - 1
436   kas=kas + 1
      if(kas.gt.npox) go to 437
      if(kas.eq.kar.or.kas.eq.lar) go to 436
      if(kas.eq.mar.or.jcon(kas).eq.2) go to 436
      jcon(kas) = jcon(kas) + 1
      las=0
438   las = las + 1
      if(las.gt.nbox) go to 439
      if(las.eq.kar.or.las.eq.lar) go to 438
      if(las.eq.mar.or.las.eq.kas) go to 438
      if(jcon(las).ne.0) go to 438
      jcon(las) =2
      go to 33
439   jcon(kas) =jcon(kas) - 1
      go to 436
401   jcon(las) =0
      go to 438
437   jcon(mar) = jcon(mar) + 1
      go to 432
433   jcon(lar) = jcon(lar) + 1
      go to 428
429   jcon(kar) =jcon(kar) + 1
      go to 424
425   nbd=3
      kar=0
440   kar=kar + 1
      if(kar.eq.nbxa) go to 441
      if(jcon(kar).eq.2) go to 440
      jcon(kar) = jcon(kar) + 1
      iq=iab(kar)
      kpaw = min(iq,jpaw) + ideks( max(iq,jpaw))
      kpaw=jab(kpaw)
      lar=kar
444   lar= lar + 1
      if(lar.eq.nbox) go to 445
      if(jcon(lar).eq.2) go to 444
      jcon(lar) = jcon(lar) + 1
      iq=iab(lar)
      lpaw = min(iq,kpaw) + ideks( max(iq,kpaw))
      lpaw=jab(lpaw)
      mar=lar
448   mar= mar + 1
      if(mar.gt.nbox) go to 449
      if(jcon(mar).eq.2) go to 448
      iq=iab(mar)
      mpaw = min(iq,lpaw) + ideks( max(iq,lpaw))
      mpaw=jab(mpaw)
      mpaw=inver(mpaw)
      if(mpaw.eq.0) go to 448
      kas = ilee(mpaw)
      npox=ilee(mpaw + 1)
      jcon(mar)=jcon(mar)  + 1
452   kas = kas + 1
      if(kas.gt.npox) go to 453
      if(kas.eq.kar.or.kas.eq.lar) go to 452
      if(kas.eq.mar.or.jcon(kas).eq.0) go to 452
      jcon(kas) = jcon(kas) - 1
      las =0
454   las = las + 1
      if(las.gt.nbox) go to 455
      if(las.eq.kar.or.las.eq.lar ) go to 454
      if(las.eq.mar.or.las.eq.kas) go to 454
      if(jcon(las).ne.2) go to 454
      jcon(las) = 0
      go to 33
455   jcon(kas) =jcon(kas) + 1
      go to 452
402   jcon(las) = 2
      go to 454
453   jcon(mar) = jcon(mar) - 1
      go to 448
449   jcon(lar) = jcon(lar) - 1
      go to 444
445   jcon(kar) = jcon(kar) - 1
      go to 440
441   nbd = 4
      kar=0
456   kar = kar + 1
      if(kar.eq.nbxa) go to 582
      if(jcon(kar).eq.2) go to 456
      jcon(kar) = jcon(kar)  + 1
      iq=iab(kar)
      kpaw = min(iq,jpaw) + ideks( max(iq,jpaw))
      kpaw=jab(kpaw)
      lar= kar
460   lar=lar + 1
      if(lar.eq.nbox) go to 461
      if(jcon(lar).eq.2) go to 460
      jcon(lar) = jcon(lar ) + 1
      iq=iab(lar)
      lpaw = min(iq,kpaw) + ideks( max(iq,kpaw))
      lpaw=jab(lpaw)
      mar = lar
464   mar = mar + 1
      if(mar.gt.nbox) go to 465
      if(jcon(mar).eq.2) go to 464
      iq=iab(mar)
      mpaw = min(iq,lpaw) + ideks( max(iq,lpaw))
      mpaw= jab(mpaw)
      jcon(mar) = jcon(mar)  + 1
      kas =0
468   kas = kas + 1
      if(kas.eq.nbxa) go to 469
      if(kas.eq.kar.or.kas.eq.lar) go to 468
      if(kas.eq.mar.or.jcon(kas).eq.0) go to 468
      jcon(kas) = jcon(kas) -1
      iq=iab(kas)
      npaw = min(iq,mpaw) + ideks( max(iq,mpaw))
      npaw=jab(npaw)
      las = kas
472   las = las + 1
      if(las.eq.nbox) go to 473
      if(las.eq.kar.or.las.eq.lar) go to 472
      if(las.eq.mar.or.jcon(las).eq.0) go to 472
      iq=iab(las)
      ipaw = min(iq,npaw) + ideks( max(iq,npaw))
      ipaw=jab(ipaw)
      ipaw=inver(ipaw)
      if(ipaw.eq.0) go to 472
      mom= ilee(ipaw)
      npox=ilee(ipaw + 1)
      if(mom.ge.las) go to 478
      if(las.ge.npox) go to 472
      mom = las
478   jcon(las) = jcon(las)  - 1
476   mom = mom + 1
      if(mom.gt.npox) go to 477
      if(mom.eq.kar.or.mom.eq.lar) go to 476
      if(mom.eq.mar.or.jcon(mom).eq.0) go to 476
      jcon(mom) = jcon(mom) - 1
      go to 33
477   jcon(las) = jcon(las) + 1
      go to 472
403   jcon(mom) = jcon(mom) + 1
      go to 476
473   jcon(kas) = jcon(kas) + 1
      go to 468
469   jcon(mar) = jcon(mar) - 1
      go to 464
465   jcon(lar) = jcon(lar) - 1
      go to 460
461   jcon(kar) = jcon(kar) - 1
      go to 456
582   if(mxex.eq.3) return
      nbd= 5
      if(jpaw.ne.1) go to 479
      kar=0
480   kar=kar + 1
      if(kar.eq.nbox) go to 479
      if(jcon(kar).ne.2) go to 480
      jcon(kar) = 0
      lar =kar
481   lar=lar + 1
      if(lar.gt.nbox) go to 482
      if(jcon(lar).ne.2) go to 481
      jcon(lar) = 0
      kas =0
483   kas = kas + 1
      if(kas.eq.nbox) go to 484
      if(kas.eq.kar.or.kas.eq.lar) go to 483
      if(jcon(kas).ne.0) go to 483
      jcon(kas) = 2
      las = kas
485   las = las + 1
      if(las.gt.nbox) go to 486
      if(las.eq.kar.or.las.eq.lar) go to 485
      if(jcon(las).ne.0)  go to 485
      jcon(las) = 2
      go to 33
486   jcon(kas) = 0
      go to 483
404   jcon(las) = 0
      go to 485
484   jcon(lar) = 2
      go to 481
482   jcon(kar) = 2
      go to 480
479   nbd= 6
      kar=0
487   kar=kar + 1
      if(kar.eq.nbox) go to 488
      if(jcon(kar).eq.0) go to 487
      iq=iab(kar)
      kpaw = min(iq,jpaw) + ideks( max(iq,jpaw))
      kpaw=jab(kpaw)
      kpaw= inver(kpaw)
      if(kpaw.eq.0) go to 487
      lar=ilee(kpaw)
      npox= ilee(kpaw + 1)
      if(lar.ge.kar) go to 491
      if(kar.ge.npox)  go to 487
      lar =kar
491   jcon(kar) = jcon(kar) - 1
492   lar =lar + 1
      if(lar.gt.npox) go to 493
      if(jcon(lar).eq.0) go to 492
      jcon(lar) = jcon(lar) - 1
      mar =0
494   mar = mar + 1
      if(mar.gt.nbox) go to 495
      if(mar.eq.kar.or.mar.eq.lar) go to 494
      if(jcon(mar).ne.2) go to 494
      jcon(mar) = 0
      kas=0
496   kas = kas + 1
      if(kas.eq.nbox) go to 497
      if(kas.eq.kar.or.kas.eq.lar) go to 496
      if(kas.eq.mar.or.jcon(kas).ne.0) go to 496
      jcon(kas) =2
      las = kas
498   las=las + 1
      if(las.gt.nbox) go to 499
      if(las.eq.kar.or.las.eq.lar) go to 498
      if(las.eq.mar.or.jcon(las).ne.0) go to 498
      jcon(las) =2
      go to 33
499   jcon(kas) =0
      go to 496
405   jcon(las) =0
      go to 498
497   jcon(mar) = 2
      go to 494
495   jcon(lar) =jcon(lar)  + 1
      go to 492
493   jcon(kar) =jcon(kar) + 1
      go to 487
488   nbd=7
      kar=0
532   kar = kar + 1
      if(kar.eq.nbox) go to 533
      if(jcon(kar).eq.2) go to 532
      iq=iab(kar)
      kpaw = min(iq,jpaw) + ideks( max(iq,jpaw))
      kpaw=jab(kpaw)
      kpaw= inver(kpaw)
      if(kpaw.eq.0) go to 532
      lar=ilee(kpaw)
      npox= ilee(kpaw + 1)
      if (lar.ge.kar) go to 536
      if (kar.ge.npox) go to 532
      lar = kar
536   jcon(kar) = jcon(kar) + 1
537   lar = lar + 1
      if(lar.gt.npox) go to 538
      if(jcon(lar).eq.2) go to 537
      jcon(lar) =jcon(lar) + 1
      mar=0
539   mar=mar + 1
      if(mar.gt.nbox) go to 540
      if(mar.eq.kar.or.mar.eq.lar) go to 539
      if(jcon(mar).ne.0) go to 539
      jcon(mar) =2
      kas =0
541   kas = kas + 1
      if(kas.eq.nbox) go to 542
      if(kas.eq.kar.or.kas.eq.lar) go to 541
      if(kas.eq.mar.or.jcon(kas).ne.2) go to 541
      jcon(kas) =0
      las = kas
543   las = las + 1
      if(las.gt.nbox) go to 544
      if(las.eq.kar.or.las.eq.lar) go to 543
      if (las.eq.mar.or.jcon(las).ne.2) go to 543
      jcon(las) =0
      go to 33
544   jcon(kas) =2
      go to 541
406   jcon(las) =2
      go to 543
542   jcon(mar) =0
      go to 539
540   jcon(lar) =jcon(lar) - 1
      go to 537
538   jcon(kar) =jcon(kar) -1
      go to 532
533   nbd= 8
      kar =0
545   kar = kar +1
      if (kar.eq.nbox) go to 546
      if (jcon(kar).eq.0) go to 545
      jcon(kar) =jcon (kar) - 1
      iq= iab(kar)
      kpaw = min(iq,jpaw) + ideks( max(iq,jpaw))
      kpaw=jab(kpaw)
      lar = kar
548   lar = lar + 1
      if(lar.gt.nbox) go to 549
      if(jcon(lar).eq.0 ) go to 548
      jcon(lar) =jcon(lar)  - 1
      iq = iab(lar)
      lpaw = min(iq,kpaw) + ideks( max(iq,kpaw))
      lpaw= jab(lpaw)
      mar =0
552   mar = mar  + 1
      if(mar.eq.nbox)  go to 553
      if (mar.eq.kar.or.mar.eq.lar) go to 552
      if(jcon(mar).eq.2)  go to 552
      iq=iab(mar)
      mpaw = min(iq,lpaw) + ideks( max(iq,lpaw))
      mpaw=jab(mpaw)
      mpaw=inver(mpaw)
      if(mpaw.eq.0)  go to 552
      kas = ilee(mpaw)
      npox= ilee(mpaw + 1)
      if (kas.ge.mar)   go to 556
      if (mar.ge.npox)   go to 552
      kas = mar
 556  jcon(mar)  = jcon(mar)  + 1
557   kas  = kas   + 1
      if(kas.gt.npox)   go to 558
      if(kas.eq.kar.or.kas.eq.lar)   go to 557
      if (jcon(kas).eq.2)   go to 557
      jcon(kas)  =jcon(kas)   + 1
      las =0
559   las  = las  + 1
      if(las.gt.nbox)  go to 560
      if(las.eq.kar.or.las.eq.lar)  go to 559
      if(las.eq.mar.or.las.eq.kas)   go to 559
      if (jcon(las).ne.2)   go to 559
      jcon(las)  =0
      mom = 0
561   mom = mom  + 1
      if(mom.gt.nbox)   go to 562
      if(mom.eq.kar.or.mom.eq.lar)   go to 561
      if(mom.eq.mar.or.mom.eq.kas)   go to 561
      if(mom.eq.las.or.jcon(mom).ne.0)   go to 561
      jcon(mom) = 2
      go to 33
562   jcon(las)   =2
      go to 559
407   jcon(mom)   = 0
      go to 561
560   jcon(kas)  =jcon(kas)  - 1
      go to 557
558   jcon(mar)  =jcon(mar)   - 1
      go to 552
553   jcon(lar)   =jcon(lar)  + 1
      go to 548
549   jcon(kar)   =jcon(kar)  + 1
      go to 545
546   nbd= 9
      kar =0
801   kar = kar  + 1
      if(kar.eq.nbxb)  go to 563
      if(jcon(kar).eq.0 )   go to 801
      jcon(kar) = jcon(kar)    - 1
      iq= iab(kar)
      kpaw = min(iq,jpaw) + ideks( max(iq,jpaw))
      kpaw=jab(kpaw)
      lar = kar
566   lar = lar + 1
      if (lar.eq.nbxa)  go to 567
      if (jcon(lar).eq.0)  go to 566
      jcon(lar)  = jcon(lar)   - 1
      iq=iab(lar)
      lpaw = min(iq,kpaw) + ideks( max(iq,kpaw))
      lpaw=jab(lpaw)
      mar = lar
570   mar = mar + 1
      if(mar.eq.nbox)  go to 571
      if (jcon(mar).eq. 0) go to 570
      iq=iab(mar)
      mpaw = min(iq,lpaw) + ideks( max(iq,lpaw))
      mpaw= jab(mpaw)
      mpaw=inver(mpaw)
      if(mpaw.eq.0)  go to 570
      kas = ilee(mpaw)
      npox=ilee(mpaw + 1)
      if(kas.ge.mar)   go to 574
      if(mar.ge.npox)   go to 570
      kas = mar
574   jcon(mar)  = jcon(mar)  - 1
575   kas = kas + 1
      if(kas.gt.npox) go to 576
      if(jcon(kas).eq.0)   go to 575
      jcon(kas)   = jcon(kas)  - 1
      las =0
577   las = las + 1
      if(las.eq.nbox)  go to 578
      if(las.eq.kar.or.las.eq.lar)   go to 577
      if (las.eq.mar.or.las.eq.kas)   go to 577
      if (jcon(las).ne.0) go to 577
      jcon(las)  = 2
      mom= las
579   mom=mom + 1
      if(mom.gt.nbox)  go to 580
      if(mom.eq.kar.or.mom.eq.lar)   go to 579
      if (mom.eq.mar.or.mom.eq.kas)  go to 579
      if (jcon(mom).ne.0)  go to 579
      jcon(mom)   = 2
      go to 33
580   jcon(las)   =0
      go to 577
408   jcon(mom)   =0
      go to 579
578   jcon(kas) = jcon(kas)  + 1
      go to 575
576   jcon(mar)   =jcon(mar)  + 1
      go to 570
571   jcon(lar)  = jcon(lar)  + 1
      go to 566
567   jcon(kar) = jcon(kar)  + 1
      go to 801
563   nbd=10
      kar =0
583   kar = kar + 1
      if(kar.eq.nbxb)  go to 584
      if(jcon(kar).eq.2)  go to 583
      jcon(kar) = jcon(kar)  + 1
      iq= iab(kar)
      kpaw = min(iq,jpaw) + ideks( max(iq,jpaw))
      kpaw= jab(kpaw)
      lar = kar
587   lar = lar +1
      if(lar.eq.nbxa)  go to 588
      if(jcon(lar).eq.2)  go to 587
      jcon(lar) = jcon(lar)  + 1
      iq= iab(lar)
      lpaw = min(iq,kpaw) + ideks( max(iq,kpaw))
      lpaw=jab(lpaw)
      mar = lar
591   mar = mar + 1
      if(mar.eq.nbox)  go to 592
      if(jcon(mar).eq.2)  go to 591
      iq= iab(mar)
      mpaw = min(iq,lpaw) + ideks( max(iq,lpaw))
      mpaw= jab(mpaw)
      mpaw=inver(mpaw)
      if(mpaw.eq.0)    go to 591
      kas = ilee(mpaw)
      npox=ilee(mpaw + 1)
      if(kas.ge.mar)   go to 595
      if(mar.ge.npox)  go to 591
      kas = mar
595   jcon(mar) =jcon(mar) + 1
596   kas = kas + 1
      if(kas.gt.npox)  go to 597
      if(jcon(kas).eq.2)  go to 596
      jcon(kas) = jcon(kas) + 1
      las =0
598   las =las + 1
      if(las.eq.nbox)  go to 599
      if(las.eq.kar.or.las.eq.lar)  go to 598
      if(las.eq.mar.or.las.eq.kas)  go to 598
      if(jcon(las).ne.2)  go to 598
      jcon(las) = 0
       mom = las
600   mom = mom + 1
      if(mom.gt.nbox)  go to 601
      if(mom.eq.kar.or.mom.eq.lar)  go to 600
      if(mom.eq.mar.or.mom.eq.kas)   go to 600
      if(jcon(mom).ne.2)   go to 600
      jcon(mom) = 0
      go to 33
601   jcon(las) = 2
      go to 598
409   jcon(mom) = 2
      go to 600
599   jcon(kas) = jcon(kas) - 1
      go to 596
597   jcon(mar) = jcon(mar)  - 1
      go to 591
592   jcon(lar) =jcon(lar) - 1
      go to 587
588   jcon(kar) = jcon(kar) - 1
      go to 583
584   nbd=11
      kar=0
602   kar = kar  + 1
      if(kar.eq.nbxb)  go to 603
      if(jcon(kar).eq.0)  go to 602
      jcon(kar) = jcon(kar) - 1
      iq=iab(kar)
      kpaw = min(iq,jpaw) + ideks( max(iq,jpaw))
      kpaw=jab(kpaw)
      lar = kar
606   lar = lar + 1
      if(lar.eq.nbxa)  go to 607
      if(jcon(lar).eq.0)  go to 606
      jcon(lar) = jcon(lar) - 1
      iq=iab(lar)
      lpaw = min(iq,kpaw) + ideks( max(iq,kpaw))
      lpaw=jab(lpaw)
      mar = lar
610   mar = mar + 1
      if(mar.eq.nbox)  go to 611
      if(jcon(mar).eq.0)  go to 610
      jcon(mar) = jcon(mar)  - 1
      iq= iab(mar)
      mpaw = min(iq,lpaw) + ideks( max(iq,lpaw))
      mpaw= jab(mpaw)
      kas = mar
614   kas = kas + 1
      if(kas.gt.nbox)  go to 615
      if(jcon(kas).eq.0)   go to 614
      jcon(kas)= jcon(kas) - 1
      iq= iab(kas)
      npaw = min(iq,mpaw) + ideks( max(iq,mpaw))
      npaw= jab(npaw)
      las =0
618   las = las + 1
      if(las.eq.nbox) go to 619
      if( las.eq.kar.or.las.eq.lar)   go to 618
      if (las.eq.mar.or.las.eq.kas)   go to 618
      if(jcon(las).eq.2)   go to 618
      iq= iab(las)
      ipaw = min(iq,npaw) + ideks( max(iq,npaw))
      ipaw= jab(ipaw)
      ipaw= inver(ipaw)
      if(ipaw.eq.0)   go to 618
      mas= ilee(ipaw)
      npox= ilee(ipaw + 1)
      if(mas.ge.las)   go to 622
      if(las.ge.npox)   go to 618
      mas = las
622   jcon(las) = jcon(las)  + 1
623   mas = mas + 1
      if(mas.gt.npox)   go to 624
      if(jcon(mas).eq.2)  go to 623
      if(mas.eq.kar.or.mas.eq.lar)  go to 623
      if(mas.eq.mar.or.mas.eq.kas)   go to 623
      jcon(mas) = jcon(mas) + 1
      mom =0
625   mom = mom + 1
      if(mom.gt.nbox)   go to 626
      if(jcon(mom).ne.0) go to 625
      if(mom.eq.kar.or.mom.eq.lar) go to 625
      if(mom.eq.mar.or.mom.eq.kas) go to 625
      if(mom.eq.las.or.mom.eq.mas) go to 625
      jcon(mom) =2
      go to 33
626   jcon(mas) = jcon(mas)  - 1
      go to 623
410   jcon(mom) = 0
      go to 625
624   jcon(las) = jcon(las) - 1
      go to 618
619   jcon(kas) = jcon(kas) + 1
      go to 614
615   jcon(mar) = jcon(mar)  + 1
      go to 610
611   jcon(lar) = jcon(lar) + 1
      go to 606
607   jcon(kar) =jcon(kar) + 1
      go to 602
603   nbd=12
      kar=0
627   kar = kar + 1
      if(kar.eq.nbxb) go to 628
      if(jcon(kar).eq.2)  go to 627
      jcon(kar) = jcon(kar) + 1
      iq=iab(kar)
      kpaw = min(iq,jpaw) + ideks( max(iq,jpaw))
      kpaw= jab(kpaw)
      lar = kar
631   lar = lar + 1
      if(lar.eq.nbxa) go to 632
      if(jcon(lar).eq.2) go to 631
      jcon(lar) = jcon(lar) + 1
      iq= iab(lar)
      lpaw = min(iq,kpaw) + ideks( max(iq,kpaw))
      lpaw= jab(lpaw)
      mar = lar
635   mar = mar + 1
      if(mar.eq.nbox)  go to 636
      if(jcon(mar).eq.2)   go to 635
      jcon(mar) = jcon(mar)  + 1
      iq=iab(mar)
      mpaw = min(iq,lpaw) + ideks( max(iq,lpaw))
      mpaw= jab(mpaw)
      kas = mar
639   kas = kas + 1
      if(kas.gt.nbox)  go to 640
      if(jcon(kas).eq.2)   go to 639
      jcon(kas) = jcon(kas) + 1
      iq= iab(kas)
      npaw = min(iq,mpaw) + ideks( max(iq,mpaw))
      npaw= jab(npaw)
      las=0
643   las = las + 1
      if(las.eq.nbox) go to 644
      if(las.eq.kar.or.las.eq.lar)   go to 643
      if(las.eq.mar.or.las.eq.kas)   go to 643
      if(jcon(las).eq.0)   go to 643
      iq= iab(las)
      ipaw = min(iq,npaw) + ideks( max(iq,npaw))
      ipaw=jab(ipaw)
      ipaw= inver(ipaw)
      if(ipaw.eq.0)  go to 643
      mas = ilee(ipaw)
      npox= ilee(ipaw + 1)
      if(mas.ge.las)   go to 647
      if(las.ge.npox)   go to 643
      mas = las
647   jcon(las)=jcon(las)  - 1
648   mas = mas + 1
      if(mas.gt.npox)  go to 649
      if(jcon(mas).eq.0)   go to 648
      if(mas.eq.kar.or.mas.eq.lar)   go to 648
      if(mas.eq.mar.or.mas.eq.kas)   go to 648
      jcon(mas) = jcon(mas)   - 1
      mom = 0
650   mom = mom + 1
      if(mom.gt.nbox)   go to 651
      if(jcon(mom).ne.2)   go to 650
      if(mom.eq.kar.or.mom.eq.lar)   go to 650
      if(mom.eq.mar.or.mom.eq.kas)   go to 650
      if (mom.eq.las.or.mom.eq.mas)   go to 650
      jcon(mom)   = 0
      go to 33
651   jcon(mas)   = jcon(mas)   + 1
      go to 648
411   jcon(mom)   = 2
      go to 650
649   jcon(las) = jcon(las) + 1
      go to 643
644   jcon(kas) = jcon(kas)  - 1
      go to 639
640   jcon(mar) = jcon(mar)  - 1
      go to 635
636   jcon(lar) = jcon(lar)   - 1
      go to 631
632   jcon(kar) = jcon(kar)   -1
      go to 627
628   nbd=13
      kar = 0
652   kar = kar + 1
      if(kar.eq.nbxb)   return
      if(jcon(kar).eq.0)  go to 652
      jcon(kar) = jcon(kar) - 1
      iq= iab(kar)
      kpaw = min(iq,jpaw) + ideks( max(iq,jpaw))
      kpaw= jab(kpaw)
      lar = kar
656   lar = lar + 1
      if(lar.eq.nbxa)  go to 657
      if(jcon(lar).eq.0)   go to 656
      jcon(lar) = jcon(lar) -  1
      iq=iab(lar)
      lpaw = min(iq,kpaw) + ideks( max(iq,kpaw))
      lpaw= jab(lpaw)
      mar = lar
660   mar = mar + 1
      if(mar.eq.nbox)   go to 661
       if(jcon(mar).eq.0)  go to 660
      jcon(mar) = jcon(mar)   -   1
      iq= iab(mar)
      mpaw = min(iq,lpaw) + ideks( max(iq,lpaw))
      mpaw= jab(mpaw)
      kas = mar
664   kas = kas + 1
      if(kas.gt.nbox)   go to 665
      if(jcon(kas).eq.0)   go to 664
      jcon(kas) = jcon(kas)   -1
      iq= iab(kas)
      npaw = min(iq,mpaw) + ideks( max(iq,mpaw))
      npaw= jab(npaw)
      las =0
668   las = las + 1
      if(las.eq.nbxb)  go to 669
      if(las.eq.kar.or.las.eq.lar)   go to 668
      if(las.eq.mar.or.las.eq.kas)   go to 668
      if(jcon(las).eq.2)   go to 668
      jcon(las) = jcon(las)   +1
      iq= iab(las)
      ipaw = min(iq,npaw) + ideks( max(iq,npaw))
      ipaw= jab(ipaw)
      mas = las
672   mas = mas + 1
      if(mas.eq.nbxa)   go to 673
      if(jcon(mas).eq.2)  go to 672
      if(mas.eq.kar.or.mas.eq.lar)   go to 672
      if(mas.eq.mar.or.mas.eq.kas)   go to 672
      jcon(mas) = jcon(mas)  +  1
      iq= iab(mas)
      ipop = min(iq,ipaw) + ideks( max(iq,ipaw))
      ipop= jab(ipop)
      mom = mas
676   mom = mom + 1
      if(mom.eq.nbox)   go to 677
      if(jcon(mom).eq.2)   go to 676
      if(mom.eq.kar.or.mom.eq.lar)   go to 676
      if(mom.eq.mar.or.mom.eq.kas) go to 676
      iq= iab(mom)
      idad = min(iq,ipop) + ideks( max(iq,ipop))
      idad= jab(idad)
      idad= inver(idad)
      if(idad.eq.0) go to 676
      mel= ilee(idad)
      npox= ilee(idad + 1)
      if(mel.ge.mom) go to 680
      if(mom.ge.npox)  go to 676
      mel= mom
680   jcon(mom) = jcon(mom)  +  1
681   mel= mel + 1
      if(mel.gt.npox)  go to 682
      if(jcon(mel).eq.2)   go to 681
      if(mel.eq.kar.or.mel.eq.lar)   go to 681
      if(mel.eq.mar.or.mel.eq.kas)   go to 681
      jcon(mel) = jcon(mel) + 1
      go to 33
682   jcon(mom) = jcon(mom)  -  1
      go to 676
412   jcon(mel) = jcon(mel) -  1
      go to 681
677   jcon(mas) = jcon(mas)     -1
      go to 672
673   jcon(las) = jcon(las)   -   1
      go to 668
669   jcon(kas) = jcon(kas)   +   1
      go to 664
665   jcon(mar) = jcon(mar)   +   1
      go to 660
661   jcon(lar) = jcon(lar)   +   1
      go to 656
657   jcon (kar) = jcon(kar)   +   1
      go to 652
33    if (klx.eq.1) go to 25
      lly=-nshl
      do 34 i=1,kly
      mx=0
      lly=lly+nshl
      np=nop(i)
      ja=ndub(i)
      llz=lly
      if (np.eq.0) go to 173
      do 172 j=1,np
      llz=llz+1
      jd=jkon(llz)
      if (jcon(jd).gt.0) go to 172
      mx=mx+1
      if (mx.gt.mxex) go to 34
 172  continue
      if (ja.eq.0) go to 36
 173  do 174 j=1,ja
      llz=llz+1
      jd=jkon(llz)
      jd=2-jcon(jd)
      if (jd.eq.0) go to 174
      mx=mx+jd
      if (mx.gt.mxex) go to 34
 174  continue
      go to 36
34    continue
25    nod=0
      ncl=0
      do 76 i=1,nbox
      jc=jcon(i)
      if (jc.eq.0) go to 76
      if (jc.eq.1) go to 77
      ncl=ncl+1
      iduf(ncl)=i
      go to 76
77    nod=nod+1
      icuf(nod)=i
76    continue
      nmns=nod/2+1-ispin
      if(nmns.le.0) goto 36
      ig=ig+1
      kcon(ig)=nmns
      nconf(nmns)=nconf(nmns)+1
      if (ig.lt.ibl) go to 85
      ig=0
      write(mdisk)kcon
 85   if (nod.eq.0) go to 80
      do 86 i=1,nod
      ig=ig+1
      kcon(ig)=icuf(i)
      if (ig.lt.ibl) go to 86
      ig=0
      write(mdisk) kcon
86    continue
      if (ncl.eq.0) go to 36
80    do 81 i=1,ncl
      ig=ig+1
      kcon(ig)=iduf(i)
      if (ig.lt.ibl) go to 81
      ig=0
      write(mdisk) kcon
81    continue
36    go to (400,401,402,403,404,405,406,407,408,409,410,411,
     *412), nbd
      return
      end
      subroutine smrd10(iot,ideks)
      implicit REAL (a-h,o-z), integer(i-n)
      integer olab
_IF(cray,ksr,i8)
      integer olab8
_ENDIF
INCLUDE(common/sizes)
      common /b/ iwod,nid,lg,irsf,ntape,ntype,ical,ntab,nc,lt,lz,nx,
     + nd,nston,itape,jtape,mdisk,ideli,mtape,iswh,nr,mm,jblk,jto,
     + igmax,nr1,ltype,linf,kfile,nawk,isec,nshl,
     + nsec,lulu,nrootx,lsng,nprin,ipt0,mconf(5),nconf(5),mh,jswh,
     + id(mxcsf),ie(mxcsf),nsc(mxcsf),
     + imo,mc,mz,ndx,nrx,ntx,ny,nt,iht,ihz,nz,mt,nzx,jy,
     + jx,nd1,nd2,md1,md2,nr3,nr4,mr1,mr2,nt1,nt2,mt1,mt2,ih1,ih2,ih3,
     + ih4,jh2,jh3,jh4,iab,nt3,iac,nr2,m,nad,nad1,ifr1,mx,md,mr,idum1(1)
      common/scrtch/jdeks(9),kdeks(8),
     + lj(8),olab(2000),nit(666),nytl(5),ndub(5),
     + icon(5),nod(5),jcon(maxorb),nir(maxorb),loc(maxorb)
      common/junk/ican(3),icbn,olab8(250),jdum(8200),mkon(mxcrec),
     + kc(mxcrec),kd(mxcrec),lab(3),ij(8),jtest(mxnshl),nyzl(mxcsf),
     + jkan(mxcrec),nplu(5),jkon(mxnshl*mxcsf)
      dimension iot(*),ideks(*)
      mm=0
      jto=0
      jblk=0
      do 14 j=1,imo
  14  jcon(j)=0
      do 11 i=1,iswh
      nc=nconf(i)
      if (nc.eq.0) go to 11
      lt=icon(i)
      nx=nytl(i)
      nd=ndub(i)
      nr=nod(i)
      iab=ideks(nr)
      nr1=nr-1
      nr2=nr-2
      lz=lt+nr
      if (i.lt.3) go to 12
      j8=i-2
      mc=nconf(j8)
      if (mc.eq.0) go to 100
      mx=nx-2
      md=nd+2
      mr=nr-4
      nt=lt
      nz=lz
      jt=icon(j8)
      iz=mr+jt
      do 13 j=1,nc
      it=nt
      mt=jt
      mz=iz
      do 15 k=1,nr
      it=it+1
      ll=iot(it)
   15 jcon(ll)=1
      if (nd.eq.0) go to 16
      do 18 k=1,nd
      it=it+1
      ll=iot(it)
   18 jcon(ll)=2
   16 do 17 k=1,mc
      jz=mz
      do 19 l=1,md
      jz=jz+1
      ll=iot(jz)
      if (jcon(ll)-1) 20,21,19
   21 ia=l+1
      go to 22
   19 continue
   20 mm=mm+1
      olab(mm)=0
      if (mm.lt.nid) go to 61
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      mm=0
      go to 61
   22 do 23 l=ia,md
      jz=jz+1
      nn=iot(jz)
      if (jcon(nn)-1) 20,24,23
   24 ja=l+1
      go to 25
   23 continue
   25 if (ja.gt.md) go to 26
      do 27 l=ja,md
      jz=jz+1
      ii=iot(jz)
      if (jcon(ii).ne.2) go to 20
   27 continue
   26 if (mr.eq.0) go to 28
      jz=mt
      do 29 l=1,mr
      jz=jz+1
      ii=iot(jz)
      if (jcon(ii).ne.1) go to 20
   29 continue
   28 jz=nt
      kz=mt+1
      jq=0
      if (mr.gt.0) iq=iot(kz)
      do 30 l=1,nr
      jz=jz+1
      ii=iot(jz)
      if (ii.eq.ll) go to 31
      if (jq.eq.mr.or.iq.ne.ii) go to 32
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   30 continue
   31 ip=l
      ia=l+1
      do 33 l=ia,nr
      jz=jz+1
      ii=iot(jz)
      if (ii.eq.nn) go to 34
      if (jq.eq.mr.or.iq.ne.ii) go to 35
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   33 continue
   32 ip=l
      ia=l+1
      do 36 l=ia,nr
      jz=jz+1
      jj=iot(jz)
      if (jj.eq.ll) go to 37
      if (jq.eq.mr.or.iq.ne.jj) go to 38
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   36 continue
   34 jp=ideks(l-1)
      ia=l+1
      do 39 l=ia,nr
      jz=jz+1
      ii=iot(jz)
      if (jq.eq.mr.or.iq.ne.ii) go to 40
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   39 continue
   40 kp=jdeks(l-2)
      ia=l+1
      do 41 l=ia,nr
      jz=jz+1
      jj=iot(jz)
      if (jq.eq.mr.or.iq.ne.jj) go to 42
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   41 continue
   42 mm=mm+1
      olab(mm)=1
      lp=kdeks(l-3)
      if (mm.lt.nid) go to 43
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      mm=0
      go to 43
   35 jp=ideks(l-1)
      ia=l+1
      do 44l=ia,nr
      jz=jz+1
      jj=iot(jz)
      if (jj.eq.nn) go to 45
      if (jq.eq.mr.or.jj.ne.iq) go to 46
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   44 continue
   37 jp=ideks(l-1)
      ia=l+1
      do 47 l=ia,nr
      jz=jz+1
      jj=iot(jz)
      if (jj.eq.nn) go to 48
      if (jq.eq.mr.or.jj.ne.iq) go to 49
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   47 continue
   38 jp=ideks(l-1)
      ia=l+1
      do 50 l=ia,nr
      jz=jz+1
      kk=iot(jz)
      if (kk.eq.ll) go to 51
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   50 continue
   51 kp=jdeks(l-2)
      ia=l+1
      do 52 l=ia,nr
      jz=jz+1
      kk=iot(jz)
      if (kk.eq.nn) go to 42
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   52 continue
   48 kp=jdeks(l-2)
      ia=l+1
      do 53 l=ia,nr
      jz=jz+1
      jj=iot(jz)
      if (jq.eq.mr.or.jj.ne.iq) go to 54
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   53 continue
   54 lp=kdeks(l-3)
      mm=mm+1
      olab(mm)=3
      if (mm.lt.nid) go to 43
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      mm=0
      go to 43
   49 kp=jdeks(l-2)
      ia=l+1
      do 55 l=ia,nr
      jz=jz+1
      kk=iot(jz)
      if (kk.eq.nn) go to 56
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   55 continue
   56 lp=kdeks(l-3)
      mm=mm+1
      olab(mm)=2
      if (mm.lt.nid) go to 43
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      mm=0
      go to 43
   45 kp=jdeks(l-2)
      ia=l+1
      do 57 l=ia,nr
      jz=jz+1
      jj=iot(jz)
      if (jq.eq.mr.or.jj.ne.iq) go to 56
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   57 continue
   46 kp=jdeks(l-2)
      ia=l+1
      do 58 l=ia,nr
      jz=jz+1
      kk=iot(jz)
      if (kk.eq.nn) go to 54
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   58 continue
 43   mm=mm+1
      olab(mm)=ip+jp+kp+lp
      if (mm.lt.nid) go to 170
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      mm=0
 170  nt1r = nir(ll)
      nl1r = loc (ll)
      nt2r = nir (nn)
      nl2r = loc (nn)
      kix=3
      nb1r = nir(ii)
      nm1r = loc (ii)
      nb2r = nir (jj)
      nm2r = loc (jj)
  120 if (nt1r-nb1r)123,122,199
  199 iax = ideks(nt1r)+nb1r
      if (nt2r.lt.nb2r) go to 124
      ibx = ideks(nt2r)+nb2r
      if (iax-ibx) 125,126,127
  127 icx = ideks(iax) + ibx
      icq = (nl1r-1)*lj(nb1r) + nm1r
  140 ljn = lj(nb2r)
      idq = (nl2r-1)*ljn + nm2r
      icq = (icq-1)*ljn*lj(nt2r)+idq
      go to 129
  126 icx = ideks (iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn+nm1r
  141 idq = (nl2r-1)*ljn+nm2r
      go to 133
  125 icx = ideks(ibx) + iax
      icq = (nl2r-1)*lj(nb2r)+nm2r
  135 ljn = lj(nb1r)
      idq = (nl1r-1)*ljn+nm1r
      icq =(icq-1)*ljn*lj(nt1r) +idq
      go to 129
  124 ibx = ideks(nb2r) + nt2r
      if (iax-ibx) 130,131,132
  132 icx = ideks(iax)+ibx
      icq = (nl1r-1)*lj(nb1r)+nm1r
      ljn = lj(nt2r)
      idq = (nm2r-1)*ljn+nl2r
      icq = (icq-1)*ljn*lj(nb2r) + idq
      go to 129
  131 icx = ideks(iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn+nm1r
      idq = (nm2r-1)*ljn+nl2r
  133 if (icq.lt.idq) go to 134
      icq = ideks(icq)+idq
      go to 129
  134 icq = ideks(idq)+icq
      go to 129
  130 icx=ideks(ibx)+iax
      icq = (nm2r-1)*lj(nt2r)+nl2r
      go to 135
  123 iax = ideks(nb1r)+nt1r
      if (nt2r.lt.nb2r) go to 136
      ibx = ideks(nt2r)+nb2r
      if (iax-ibx) 137,138,139
  139 icx = ideks(iax)+ibx
      icq = (nm1r-1)*lj(nt1r)+nl1r
      go to 140
  138 icx = ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      go to 141
  137 icx = ideks(ibx) + iax
      icq = (nl2r-1)*lj(nb2r) +nm2r
  145 ljn = lj(nt1r)
      idq = (nm1r-1)*ljn         +nl1r
      icq = (icq-1)*ljn*lj(nb1r) +idq
      go to 129
  136 ibx = ideks(nb2r)+nt2r
      if (iax-ibx) 142,143,144
  144 icx = ideks(iax)+ibx
      icq =(nm1r-1) *lj(nt1r)+nl1r
      ljn = lj(nt2r)
      idq = (nm2r-1)*ljn+nl2r
      icq = (icq-1)*ljn*lj(nb2r)+idq
      go to 129
  143 icx=ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      idq = (nm2r-1)*ljn+nl2r
      go to 133
  142 icx = ideks(ibx)+iax
      icq =(nm2r-1)*lj(nt2r)+nl2r
      go to 145
  122 iax=ideks(nt1r+1)
      if (nl1r.lt.nm1r) go to 146
      iay=ideks(nl1r)+nm1r
      go to 147
  146 iay= ideks(nm1r)+nl1r
  147 if (nl2r.lt.nm2r) go to 148
      iby=ideks(nl2r)+nm2r
      go to 149
  148 iby=ideks(nm2r)+nl2r
  149 if (nt1r.eq.nt2r) go to 150
      ibx=ideks(nt2r+1)
      if (iax.lt.ibx) go to 151
      icx=ideks(iax)+ibx
      ljn =lj(nt2r)
      icq = (iay-1)*ideks(ljn+1)+iby
      go to 129
  151 icx=ideks(ibx)+iax
      ljn=lj(nt1r)
      icq =(iby-1)*ideks(ljn+1)+iay
      go to 129
  150 icx=ideks(iax+1)
      if (iay.lt.iby) go to 153
      icq=ideks(iay)+iby
      go to 129
  153 icq=ideks(iby)+iay
 129  icq=nit(icx)  +  icq
      idud=(icq -1)/igmax
      jto=jto + 1
      kc(jto) = idud  + 1
      icq=icq  -  idud*igmax
      kd(jto)  = icq
      if(jto.lt.iwod)  go to 160
      write(ntype)  kc,kd
      jto=0
      jblk=jblk  + 1
  160 if (kix.lt.0) go to 61
      kix=-3
      itr=nb1r
      nb1r=nb2r
      nb2r=itr
      itr=nm1r
      nm1r=nm2r
      nm2r=itr
      go to 120
   61 mz=mz+mx
   17 mt=mt+mx
      do 62 k=1,nx
      nt=nt+1
      ll=iot(nt)
   62 jcon(ll)=0
   13 nz=nt+nr
      go to 100
   12 if (i.eq.1) go to 64
  100 j8=i-1
      mc=nconf(j8)
      if (mc.eq.0) go to 64
      mx=nx-1
      md=nd+1
      mr=nr-2
      jt=icon(j8)
      iz=mr+jt
      nt=lt
      do 65 j=1,nc
      it=nt
      nt1=nt+1
      nz=nt+nr
      mt=jt
      mz=iz
      do 66 k=1,nr
      it=it+1
      ll=iot(it)
   66 jcon(ll)=1
      if (nd.eq.0) go to 67
      do 68 k=1,nd
      it=it+1
      ll=iot(it)
   68 jcon(ll)=2
   67 do 69 k=1,mc
      jz=mz
      do 70 l=1,md
      jz=jz+1
      jj=iot(jz)
      if (jcon(jj)-1) 71,72,70
   70 continue
   71 if (l.eq.md) go to 73
      ia=l+1
      do 74 l=ia,md
      jz=jz+1
      ll=iot(jz)
      if (jcon(ll).ne.2) go to 75
   74 continue
      go to 73
  75  mm=mm+1
      olab(mm)=1
      if (mm.lt.nid) go to 260
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      mm=0
      go to 260
   72 l1=l
      if (l.eq.md) go to 76
      ia=l+1
      do 78 l=ia,md
      jz=jz+1
      ll=iot(jz)
      if (jcon(ll)-1) 75,77,78
  78  continue
      go to 76
  77  l2=l
      if (l.eq.md) go to 79
      ia=l+1
      do 80 l=ia,md
      jz=jz+1
      kk=iot(jz)
      if (jcon(kk).ne.2) go to 75
   80 continue
      go to 79
   73 if (mr.eq.0) go to 81
      jz=mt
      do 82 l=1,mr
      jz=jz+1
      ll=iot(jz)
      if (jcon(ll).ne.1) go to 75
   82 continue
  81  mm=mm+1
      olab(mm)=3
      if (mm.lt.nid) go to 83
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      mm=0
  83  jz=nt
      kz=mt+1
      jq=0
      if (mr.gt.0) iq=iot(kz)
      do 84 l=1,nr
      jz=jz+1
      ll=iot(jz)
      if (jq.eq.mr.or.iq.ne.ll) go to 85
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   84 continue
   85 ia=l+1
      ip=l
      do 86 l=ia,nr
      jz=jz+1
      nn=iot(jz)
      if (jq.eq.mr.or.iq.ne.nn) go to 87
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   86 continue
  87  mm=mm+1
      olab(mm)=ideks(l-1)+ip
      if (mm.lt.nid) go to 88
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      mm=0
  88  nt1r=nir(jj)
      nl1r=loc(jj)
      nb1r=nir(ll)
      nm1r=loc(ll)
      nm2r=loc(nn)
  220 if (nt1r-nb1r) 223,222,299
  299 iax = ideks(nt1r)+nb1r
      icx = ideks (iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn
      idq = icq+nm2r
      icq=ideks(idq)+icq+nm1r
      go to 229
  223 iax = ideks(nb1r)+nt1r
      icx=ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      idq = (nm2r-1)*ljn+nl1r
      icq=ideks(idq)+icq
      go to 229
  222 iax=ideks(nt1r+1)
      icx=ideks(iax+1)
      if (nl1r.lt.nm1r) go to 246
      iay=ideks(nl1r)+nm1r
      go to 247
  246 iay= ideks(nm1r)+nl1r
      iby=ideks(nm2r)+nl1r
      go to 253
  247 if (nl1r.lt.nm2r) go to 248
      iby=ideks(nl1r)+nm2r
      go to 253
  248 iby=ideks(nm2r)+nl1r
  253 icq=ideks(iby)+iay
  229 icq=nit(icx)  +  icq
      idud=(icq -1)/igmax
      jto=jto + 1
      kc(jto) = idud  + 1
      icq=icq  -  idud*igmax
      kd(jto)  = icq
      if(jto.lt.iwod)  go to 260
      write(ntype) kc,kd
      jto=0
      jblk=jblk  + 1
      go to 260
   79 if (mr.eq.0) go to 89
      jz=mt
      do 90 l=1,mr
      jz=jz+1
      ii=iot(jz)
      if (jcon(ii)-1) 75,90,91
  90  continue
      go to 89
  91  l1=l
      if (l.eq.mr) go to 93
      ia=l+1
      do 92 l=ia,mr
      jz=jz+1
      jc=iot(jz)
      if (jcon(jc).ne.1) go to 75
  92  continue
      go to 93
  89  mm=mm+1
      olab(mm)=3
      if (mm.lt.nid) go to 94
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      mm=0
   94 jz=nt
      do 95 l=1,nr
      jz=jz+1
      if (iot(jz).eq.jj) go to 96
   95 continue
   96 ia=l+1
      ip=l
      do 97 l=ia,nr
      jz=jz+1
      if (iot(jz).eq.ll) go to 98
   97 continue
  98  mm=mm+1
      olab(mm)=-ip-ideks(l-1)
      if (mm.lt.nid) go to 99
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      mm=0
   99 jz=nz
      if (l1.eq.1) go to 180
      kz=mz
      ja=l1-1
      do 181 l=1,ja
      kz=kz+1
      jz=jz+1
      kk=iot(jz)
      if (kk.ne.iot(kz)) go to 182
  181 continue
  180 if (l2.eq.l1+1) go to 183
      kz=mz+l1
      ia=l1+1
      ja=l2-1
      do 184 l=ia,ja
      kz=kz+1
      jz=jz+1
      kk=iot(jz)
      if (kk.ne.iot(kz)) go to 182
  184 continue
  183 if (l2.lt.md) go to 185
      kk=iot(jz+1)
      go to 182
  185 kz=mz+l2
      ia=l2+1
      do 186 l=ia,md
      kz=kz+1
      jz=jz+1
      kk=iot(jz)
      if (kk.ne.iot(kz)) go to 182
  186 continue
      kk=iot(jz+1)
  182 nt1r=nir(kk)
      nl1r=loc(kk)
      nb1r=nir(jj)
      nm1r=loc(jj)
      nm2r=loc(ll)
      go to 220
  93  mm=mm+1
      olab(mm)=2
      if (mm.lt.nid) go to 187
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      mm=0
 187  jz=nt
      kz=mt
      ig=0
      if (l1.eq.1) go to 188
      ja=l1-1
      ka=1
      do 189 l=ka,ja
      kz=kz+1
      jp=iot(kz)
 802  jz=jz+1
      if (jp.eq.iot(jz)) go to 189
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.3) go to 190
      go to 802
  189 continue
  188 if (l1.eq.mr) go to 191
      ia=l1+1
      kz=kz+1
      do 192 l=ia,mr
      kz=kz+1
      jp=iot(kz)
 803  jz=jz+1
      if (jp.eq.iot(jz)) go to 192
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.3) go to 190
      go to 803
  192 continue
  191 lab(3)=nr
      if (ig-1) 284,285,190
  284 lab(2)=nr1
      lab(1)=nr2
      go to 190
  285 lab(2)=nr1
  190 jz=lab(1)+nt
      nn=iot(jz)
      mm=mm+1
      if (nn.eq.jj) go to 194
      olab(mm)=1
  195 if (mm.lt.nid) go to 196
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  196 mm=mm+1
      olab(mm)=l1
      if (mm.lt.nid) go to 197
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  197 mm=mm+1
      jp=lab(2)-1
      kp=lab(3)-2
      olab(mm)=lab(1)+ideks(jp)+jdeks(kp)
      if (mm.lt.nid) go to 198
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      go to 198
  194 jz=lab(2)+nt
      nn=iot(jz)
      if (nn.eq.ll) go to 271
      olab(mm)=2
      go to 195
  271 jz=lab(3)+nt
      nn=iot(jz)
      olab(mm)=3
      go to 195
  198 nt1r=nir(jj)
      nt2r=nir(ll)
      nl1r=loc(jj)
      nl2r=loc(ll)
      kix=3
      nb1r=nir(nn)
      nb2r=nir(ii)
      nm1r=loc(nn)
      nm2r=loc(ii)
  320 if (nt1r-nb1r) 323,322,399
  399 iax = ideks(nt1r)+nb1r
      if (nt2r.lt.nb2r) go to 324
      ibx = ideks(nt2r)+nb2r
      if (iax-ibx) 325,326,327
  327 icx = ideks(iax) + ibx
      icq = (nl1r-1)*lj(nb1r) + nm1r
  340 ljn = lj(nb2r)
      idq = (nl2r-1)*ljn + nm2r
      icq = (icq-1)*ljn*lj(nt2r)+idq
      go to 329
  326 icx = ideks (iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn+nm1r
  341 idq = (nl2r-1)*ljn+nm2r
      go to 333
  325 icx = ideks(ibx) + iax
      icq = (nl2r-1)*lj(nb2r)+nm2r
  335 ljn = lj(nb1r)
      idq = (nl1r-1)*ljn+nm1r
      icq =(icq-1)*ljn*lj(nt1r) +idq
      go to 329
  324 ibx = ideks(nb2r) + nt2r
      if (iax-ibx) 330,331,332
  332 icx = ideks(iax)+ibx
      icq = (nl1r-1)*lj(nb1r)+nm1r
      ljn = lj(nt2r)
      idq = (nm2r-1)*ljn+nl2r
      icq = (icq-1)*ljn*lj(nb2r) + idq
      go to 329
  331 icx = ideks(iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn+nm1r
      idq = (nm2r-1)*ljn+nl2r
  333 if (icq.lt.idq) go to 334
      icq = ideks(icq)+idq
      go to 329
  334 icq = ideks(idq)+icq
      go to 329
  330 icx=ideks(ibx)+iax
      icq = (nm2r-1)*lj(nt2r)+nl2r
      go to 335
  323 iax = ideks(nb1r)+nt1r
      if (nt2r.lt.nb2r) go to 336
      ibx = ideks(nt2r)+nb2r
      if (iax-ibx) 337,338,339
  339 icx = ideks(iax)+ibx
      icq = (nm1r-1)*lj(nt1r)+nl1r
      go to 340
  338 icx = ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      go to 341
  337 icx = ideks(ibx) + iax
      icq = (nl2r-1)*lj(nb2r) +nm2r
  345 ljn = lj(nt1r)
      idq = (nm1r-1)*ljn         +nl1r
      icq = (icq-1)*ljn*lj(nb1r) +idq
      go to 329
  336 ibx = ideks(nb2r)+nt2r
      if (iax-ibx) 342,343,344
  344 icx = ideks(iax)+ibx
      icq =(nm1r-1) *lj(nt1r)+nl1r
      ljn = lj(nt2r)
      idq = (nm2r-1)*ljn+nl2r
      icq = (icq-1)*ljn*lj(nb2r)+idq
      go to 329
  343 icx=ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      idq = (nm2r-1)*ljn+nl2r
      go to 333
  342 icx = ideks(ibx)+iax
      icq =(nm2r-1)*lj(nt2r)+nl2r
      go to 345
  322 iax=ideks(nt1r+1)
      if (nl1r.lt.nm1r) go to 346
      iay=ideks(nl1r)+nm1r
      go to 347
  346 iay= ideks(nm1r)+nl1r
  347 if (nl2r.lt.nm2r) go to 348
      iby=ideks(nl2r)+nm2r
      go to 349
  348 iby=ideks(nm2r)+nl2r
  349 if (nt1r.eq.nt2r) go to 350
      ibx=ideks(nt2r+1)
      if (iax.lt.ibx) go to 351
      icx=ideks(iax)+ibx
      ljn =lj(nt2r)
      icq = (iay-1)*ideks(ljn+1)+iby
      go to 329
  351 icx=ideks(ibx)+iax
      ljn=lj(nt1r)
      icq =(iby-1)*ideks(ljn+1)+iay
      go to 329
  350 icx=ideks(iax+1)
      if (iay.lt.iby) go to 353
      icq=ideks(iay)+iby
      go to 329
  353 icq=ideks(iby)+iay
  329 icq=nit(icx)  +  icq
      idud=(icq -1)/igmax
      jto=jto + 1
      kc(jto) = idud  + 1
      icq=icq  -  idud*igmax
      kd(jto)  = icq
      if(jto.lt.iwod)  go to 360
      write(ntype) kc,kd
      jto=0
      jblk=jblk  + 1
 360  if (kix.lt.0) go to 260
      kix=-3
      itr=nb1r
      nb1r=nb2r
      nb2r=itr
      itr=nm1r
      nm1r=nm2r
      nm2r=itr
      go to 320
  76  if (mr.eq.0) go to 273
      jz=mt
      do 274 l=1,mr
      jz=jz+1
      ll=iot(jz)
      if (jcon(ll)-1) 275,274,75
 274  continue
      go to 273
 275  l1=l
      if (l.eq.mr) go to 277
      ia=l+1
      do 276 l=ia,mr
      jz=jz+1
      nn=iot(jz)
      if (jcon(nn).ne.1) go to 75
  276 continue
 277  mm=mm+1
      olab(mm)=2
      if (mm.lt.nid) go to 278
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      mm=0
 278  jz=nt
      kz=mt
      ig=0
      if (l1.eq.1) go to 279
      ja=l1-1
      do 280 l=1,ja
      kz=kz+1
      jp=iot(kz)
 254  jz=jz+1
      if (jp.eq.iot(jz)) go to 280
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.3) go to 281
      go to 254
 280  continue
 279  if (l1.eq.mr) go to 282
      ia=l1+1
      kz=kz+1
      do 283 l=ia,mr
      kz=kz+1
      jp=iot(kz)
 255  jz=jz+1
      if (jp.eq.iot(jz)) go to 283
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.3) go to 281
      go to 255
  283 continue
  282 lab(3)=nr
      if (ig-1) 286,287,281
  286 lab(2)=nr1
      lab(1)=nr2
      go to 281
  287 lab(2)=nr1
  281 jz=lab(1)+nt
      kk=iot(jz)
      mm=mm+1
      if (kk.ne.jj) go to 288
      olab(mm)=-1
      jz=lab(2)+nt
      kk=iot(jz)
      jz=lab(3)+nt
      nn=iot(jz)
 289  if(mm.lt.nid) go to 290
      mm=0
      call pack(olab8,8,olab,2000)
      write(ntape) olab8
 290  mm=mm+1
      olab(mm)=l1
      if(mm.lt.nid) go to 291
      mm=0
      call pack(olab8,8,olab,2000)
      write(ntape) olab8
 291  mm=mm+1
      jp=lab(2)-1
      kp=lab(3)-2
      olab(mm)=lab(1)+ideks(jp)+jdeks(kp)
      if(mm.lt.nid)go to 292
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      go to 292
 288  jz=lab(2)+nt
      nn=iot(jz)
      if (nn.ne.jj) go to 293
      olab(mm)=-2
      jz=lab(3)+nt
      nn=iot(jz)
      go to 289
 293  olab(mm)=-3
      go to 289
 292  nt1r=nir(jj)
      nt2r=nir(ll)
      nl1r=loc(jj)
      nl2r=loc(ll)
      kix=3
      nb1r=nir(kk)
      nb2r=nir(nn)
      nm1r=loc(kk)
      nm2r=loc(nn)
      go to 320
 273  mm=mm+1
      olab(mm)=4
      if(mm.lt.nid) go  to 294
      mm=0
      call pack(olab8,8,olab,2000)
      write(ntape) olab8
 294  jz=nt
      mm=mm+1
      if(mr.eq.0) go to 297
      kz=mt
      do 295 l=1,mr
      kz=kz+1
      jz=jz+1
      ll=iot(jz)
      if(ll.ne.iot(kz)) go to 296
 295  continue
 297  ll=iot(jz+1)
      if(ll.eq.jj) go to 298
      olab(mm)=iab
      if(mm.lt.nid) go to 371
      mm=0
      call pack(olab8,8,olab,2000)
      write(ntape) olab8
      go to 371
 298  ll=iot(jz+2)
      olab(mm)=-iab
      if(mm.lt.nid) go to 371
      mm=0
      call pack(olab8,8,olab,2000)
      write(ntape) olab8
      go to 371
 296  l1 =jz-nt
      if(ll.eq.jj) go to 373
      do 374 l=1,nr
      jz=jz+1
      if(iot(jz).eq.jj) go to 375
 374  continue
 375  olab(mm)=l1+ideks(jz-nt1)
      if(mm.lt.nid) go to 371
      mm=0
      call pack(olab8,8,olab,2000)
      write(ntape) olab8
      go to 371
 373  kz=kz-1
      do 376 ia=l,mr
      kz=kz+1
      jz=jz+1
      ll=iot(jz)
      if(ll.ne.iot(kz)) go to 377
 376  continue
      jz=jz+1
      ll=iot(jz)
 377  olab(mm)=-l1-ideks(jz-nt1)
      if(mm.lt.nid) go to 371
      mm=0
      call pack(olab8,8,olab,2000)
      write(ntape) olab8
 371  mm=mm+1
      iax=nir(jj)
      olab(mm)=iax
      if(mm.lt.nid) go to 372
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  372 iay=ideks(iax+1)
      iaz=ideks(iay+1)
      mq=lj(iax)
      mv=ideks(mq+1)
      iaq=iay-iax
      iaw=iaz-iay
      iat=nit(iaz)
      ii=loc(jj)
      kk=loc(ll)
      km=ideks(kk)
      nm=ideks(ii)
      nn=nm+ii
      if(ii.gt.kk) go to 380
      mm=mm+1
      jb=ii+km
      olab(mm)=kk
      if(mm.lt.nid) go to 378
      mm=0
      call pack(olab8,8,olab,2000)
      write(ntape) olab8
 378  mm=mm+1
      olab(mm)=ii
      if(mm.lt.nid) go to 379
      mm=0
      call pack(olab8,8,olab,2000)
      write(ntape) olab8
 379  ib=ideks(jb) +nn+iat
      go to 381
 380  lb=kk-ii
      jb=nn+lb
      ib=ideks(nn+1)+lb+iat
      mm=mm+1
      olab(mm)=ii
      if(mm.lt.nid) go to 900
      mm=0
      call pack(olab8,8,olab,2000)
      write(ntape) olab8
 900  mm=mm+1
      olab(mm)=kk
      if(mm.lt.nid) go to 381
      mm=0
      call pack(olab8,8,olab,2000)
      write(ntape) olab8
 381  idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if(jto.lt.iwod) go to 382
      jto=0
      write(ntype) kc,kd
      jblk=jblk+1
 382  if(mr.eq.0) go to 383
      jz=mt
      ir=mr
      kix=1
 472  do 384 l=1,ir
      jz=jz+1
      mi=iot(jz)
      mar=nir(mi)
      mlr=loc(mi)
      mlb=ideks(mlr)
      mla=mlb+mlr
      if (mar-iax) 385,395,396
 395  if(mla.lt.jb) go to 386
      ib=iat+ideks(mla) +jb
 388  idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto = jto+1
      kc(jto)=idud+1
      kd(jto) =ib
      if(jto.lt.iwod) go to 387
      jto =0
      write(ntype) kc,kd
      jblk=jblk+1
      go to 387
 386  ib=iat+ideks(jb) +mla
      go to 388
 387  if(mlr.lt.ii) go to 389
      kb=mlb+ii
      go to 390
 389  kb=nm+mlr
 390  if(mlr.lt.kk) go to 391
      lb=mlb+kk
      go to 392
 391  lb=mlr+km
 392  if(kb.lt.lb) go to 393
      ib=iat+ideks(kb) +lb
      go to 394
 393  ib=iat+ideks(lb) +kb
 394  idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if(jto.lt.iwod) go to 384
      jto=0
      jblk=jblk+1
      write(ntype) kc,kd
      go to 384
 385  iby=ideks(mar+1)
      iby=iaw+iby
      ibx=iaq+mar
      ibx=ideks(ibx+1)
      ibx=nit(ibx)
      iby=nit(iby)
      ml=lj(mar)
      ib=ideks(ml+1)*(jb-1) +mla+iby
      idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if(jto.lt.iwod) go to 397
      jto=0
      jblk=jblk+1
      write(ntype) kc,kd
 397  kb=(ii-1)*ml +mlr
      lb=(kk-1)*ml +mlr
      if(kb.lt.lb) go to 398
      ib=ideks(kb)+lb+ibx
      go to 471
 398  ib=ideks(lb) +kb+ibx
 471  idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto) = ib
      if(jto.lt.iwod) go to 384
      jto=0
      jblk=jblk+1
      write(ntype) kc,kd
      go to 384
396   ibx=ideks(mar+1)
      iby=ideks(ibx) +iay
      iby=nit(iby)
      ibx=ibx-mar+iax
      ibx=ideks(ibx+1)
      ibx=nit(ibx)
      ib=(mla-1)*mv +jb+iby
      idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if(jto.lt.iwod) go to 473
      jto=0
      jblk=jblk+1
      write(ntype) kc,kd
 473  kb=(mlr-1)*mq
      lb=kb+ii
      kb=kb+kk
      if(kb.lt.lb) go to 474
      ib=ideks(kb)+lb+ibx
      go to 475
 474  ib=ideks(lb)+kb+ibx
 475  idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if (jto.lt.iwod) go to 384
      jto=0
      jblk=jblk+1
      write (ntype) kc,kd
 384  continue
      if (kix.lt.0) go to 260
 383  if (nd.eq.0) go to 260
      kix=-1
      jz=nz
      ir=nd
      go to 472
 260  mz=mz+mx
  69  mt=mt+mx
      do 476 k=1,nx
      nt=nt+1
      ll=iot(nt)
  476 jcon(ll)=0
   65 continue
   64 call smr100(iot,ideks)
   11 continue
      return
      end
      subroutine smrd12(g,h,lc,ld,mxdim)
      implicit REAL (a-h,o-z), integer(i-n)
INCLUDE(common/sizes)
INCLUDE(common/prints)
      common /b/ iwod,nid,lg,irsf,ntape,ntype,ical,ntab,nc,lt,lz,nx,
     + nd,nston,itape,jtape,mdisk,ideli,mtape,iswh,nr,mm,jblk,jto,
     + igmax,nr1,mtype,linf,kfile,nawk,isec,nshl,
     + nsec,lulu,nrootx,lsng,nprin,ipt0,mconf(5),nconf(5),mh,jswh,
     + id(mxcsf),ie(mxcsf),nsc(mxcsf),
     + imo,mc,mz,ndx,nrx,ntx,ny,nt,iht,ihz,nz,mt,nzx,jy,
     + jx,nd1,nd2,md1,md2,nr3,nr4,mr1,mr2,nt1,nt2,mt1,mt2,ih1,ih2,ih3,
     + ih4,jh2,jh3,jh4,iab,nt3,iac,nr2,m,nad,nad1,ifr1, mx,md,mr
     +,idum1(1)
      common/scrtch/ kc(mxcrec),kd(mxcrec),y(2000)
      dimension g(mxdim),h(mxdim),lc(mxdim),ld(mxdim)
      if (jto.gt.0) jblk=jblk+1
      call rewftn(mtype)
      lgh=(lg-1)/igmax+1
      ib8=(lg-1)/nid+1
      igma=igmax/nid
      if(.not.oprint(29)) write(jtape,10) jblk,jto
 10   format(1x,5i6)
      igy=igmax/iwod
      nres=ib8-(lgh-1)*igma
      lx=lg-(ib8-1)*nid
      iga=(jblk-1)/igy
      jblk=jblk-iga*igy
      iga=iga+1
      igz=igy
      iw=iwod
      if(jto.eq.0) jto=iwod
      do 1 i=1,iga
      if (i.eq.iga) igz=jblk
      ig=0
      do 2 j=1,igz
      read(ntype) kc,kd
      if (i.eq.iga.and.j.eq.igz) iw=jto
      do 2 k=1,iw
      ig=ig+1
      lc(ig)=kc(k)
2     ld(ig)=kd(k)
      call rewftn(nston)
      read(nston)
      read(nston)
      nx=igma
      lz=nid
      do 3 k=1,lgh
      if(k.eq.lgh) nx=nres
      il=0
      do 4 l=1,nx
      read(nston)y
      if (k.eq.lgh.and.l.eq.nres) lz=lx
      do 4 ll=1,lz
      il=il+1
4     g(il)=y(ll)
      do 5 l=1,ig
      if (lc(l).ne.k) go to 5
      kx=ld(l)
      h(l)=g(kx)
5     continue
3     continue
      if (ig.eq.igmax) go to 6
      igy=(ig-1)/iwod+1
6     la=0
      do 7 k=1,igy
      lb=la+1
      la=la+iwod
_IF1()      le=lb+20
7     write(mtype)(h(l),l=lb,la)
1     continue
      write(mtype) h
      call rewftn(mtype)
      read (mtype) (h(l),l=1,iwod)
      call rewftn(mtype)
      if (lg-ib8*nid.eq.0) read (nston)
      return
      end
      subroutine smr133(sx,trwe,trwf,trwg,trwh,trwi,trwj,imsec)
      implicit REAL (a-h,o-z), integer(i-n)
      integer olab,ot
_IF(cray,ksr,i8)
      integer olab8
_ENDIF
INCLUDE(common/sizes)
      common/linkmr/ ic,i2,i3
      common/bufd/ gout(510),nword
      common/scrtch/ndet(5),nsac(5),iaz(5),iez(5),idra(5),
     + idrc(5),jdrc(5),jan(7),jbn(7),ispa2(7),
     + ae(7829),ot(44000),ideks(maxorb),
     + olab(2000),ij(8),nytl(5),nplu(5),ndub(5),nod(5),jkan(mxcrec),
     + nir(maxorb),loc(maxorb),sac(10),
_IF(cray,ksr,i8)
     + olab8(250),f(2304),h(mxcrec),isc(5),bob(3),itest(mxnshl),
_ELSE
     + olab8(250),f(2304),h(mxcrec),isc(5),isd,bob(3),itest(mxnshl),
_ENDIF
     + core,
     + w(mxcsf*(mxcsf+1)/2),
     + ew(mxroot),vect(mxcsf*mxroot),trsum(10,mxroot),
     + wect(mxcsf*mxroot),senk(500),iqr(mxcrec),istm(10)
     +,ht(mxroot),hs(mxcsf),ea(mxroot),iper(10),
     + coff(mxcsf*mxroot),doff(mxcsf),nrec(5),jtrk(mxcsf),
     + isym(mxroot+1),mkon(mxcrec),ikan(mxcrec),hy(mxcsf*mxcsf)
     +,ktrk(mxcsf)
      common/jany/ jerk(10),jbnk(10)
     *,idum11(10)
      common /b/ iwod,nid,lg,irsf,ntape,ntype,ical,ntab,nc,lt,lz,nx,
     + nd,nston,itape,jtape,mdisk,ideli,mtape,iswh,nr,mm,jblk,jto,
     + igmax,nr1,mtype,linf,kfile,
     + nawk,isec,nshl,
     + nsec,lulu,nrootx,lsng,nprin,ipt0,nconf(5),mconf(5),mh,jswh,
     + id(mxcsf),mn(mxcsf),nsc(mxcsf),
     + imo,mc,mz,ndx,nrx,ntx,ny,nt,iht,ihz,nz,mt,nzx,jy,
     + jx,nd1,nd2,md1,md2,nr3,nr4,mr1,mr2,nt1,nt2,mt1,mt2,ih1,ih2,ih3,
     + ih4,jh2,jh3,jh4,iab,nt3,iac,nr2,m,nad,nad1,ifr1, mx,md,mr
     + ipag
_IF1()c      common /g/ trwe,trwf,trwg,trwh,trwi,trwj,imsec
       common /d/ jbab,kml,ifl,ii,nis,iiz,nzer,kmj,ndj,jsac,ie,iej,
     * ja,kc,kmk,ndk,kr,iiq,ic1,icex,i21,i2ex,i31,i3ex,iix,
     * kzer,ndt,iek,j3,j3a,mps,nps,mms,nms,jc,jca,lcl,ndi,kmi,
     * ia,jai,nzei,nisi,lc,kk,ndl,km2,la,lzer,lr,lps,lms,
     * ic2,i22,i32,jcb,j3b
      dimension a(5292),e(2304),lout(510),
     2xemp(500),ihog(48),
     3imap(504),sx(imsec)
      equivalence  (a(1),hy(1)),(a(1),e(1)),(gout(1),lout(1))
_IF(cray,ksr,i8)
     1,(xemp(1),ndet(1)),(ihog(1),a(2809)),(imap(1),a(2305))
_ELSE
     1,(xemp(1),ndet(1)),(ihog(1),a(2557)),(imap(1),a(2305))
_ENDIF
      call rewftn(kfile)
      read (kfile)
      call rewftn(ideli)
      ltype=mtype
      call rewftn(ltype)
      ltape=linf
      call rewftn(ltape)
      read (ltape)
      xemp(5)=trwe
      xemp(6)=trwf
      xemp(7)=trwg
      xemp(8)=trwh
      xemp(9)=trwi
      xemp(10)=trwj
      i=4
  449 i=i+1
      if (i.gt.10) go to 450
      if (istm(i).gt.imsec) go to 449
      trash=xemp(i)
      isp=i
      write(jtape,451) trash
  451 format(/5x,20hcorrect threshold is,f11.8/)
      read(ideli) trwa,trwb,trwc,dumd,tdel
      read(ideli) xemp
      read (kfile) iqr
      kft=0
      imk=0
      imn=0
      do 452 i=1,iswh
      nc=nconf(i)
      if (nc.eq.0) go to 452
      nx=nytl(i)
      read(ltape)ndt,kml,a
      write(ltype) ndt,kml,a
      read(ltape) ihog,imap,e
      write(ltype) ihog,imap,e
      read(ltape) jkan
      lpix=0
      jmk=0
      ibc=1
      do 453 j=1,nc
      kft=kft+1
      kx=iqr(kft)
      if (kx.gt.4) go to 708
      if (kft.lt.iwod) go to 453
      kft=0
      read (kfile) iqr
      go to 453
 708  imk=imk+kml
      if (imk.le.irsf) go to 3273
      imk=imk-irsf
      read(ideli) xemp
3273  pmt=xemp(imk)
      if (kft.lt.iwod) go to 707
      kft=0
      read (kfile) iqr
 707  if (pmt.lt.trash) go to 454
      do 455 k=1,kml
      imn=imn+1
455   sx(imn)=pmt
      do 456 k=1,nx
      lpix=lpix+1
      jmk=jmk+1
      ikan(jmk)=jkan(lpix)
      if (lpix.lt.iwod) go to 457
      lpix=0
      read(ltape) jkan
457   if (jmk.lt.iwod) go to 456
      jmk=0
      ibc=ibc+1
      write(ltype) ikan
456   continue
      go to 453
454   lpix=lpix+nx
      if (lpix.lt.iwod) go to 453
      lpix=lpix-iwod
      read(ltape) jkan
453   continue
      ikan(jmk+1)=0
      write(ltype) ikan
      nrec(i)=ibc
452   continue
      call rewftn(kfile)
      call rewftn(ltape)
      call rewftn(ltype)
      call rewftn(ideli)
      write(ideli)trwa,trwb,trwc,trash,tdel,isp,vect,ew
      i1=1-irsf
      i2=0
  464 i1=i1+irsf
      i2=i2+irsf
      write(ideli)(sx(i),i=i1,i2)
      if (i2.lt.imn) go to 464
      write (jtape,712) nrootx,imn,istm
 712  format(5x,12i10)
      write(ideli)nrootx,((trsum(j,i),j=1,10),i=1,nrootx),istm
      read (ltape)
      do 460 i=1,iswh
      if (nconf(i).eq.0) go to 460
      read (ltype) ndt,kml,a
      write(ltape) ndt,kml,a
      read(ltype) ihog,imap,e
      write(ltape) ihog,imap,e
      ibc=nrec(i)
      do 461 j=1,ibc
      read(ltype) ikan
461   write(ltape) ikan
460   continue
      call clredx
      call closbf3
      return
 450   write(jtape,446) istm(4),imsec
 446   format(/10x,i6,'dimension exceeds maximum of',i5)
      call errors(531)
      end
      subroutine sm1313(bkay,intmax)
      implicit REAL (a-h,o-z), integer(i-n)
      integer olab,ot
_IF(cray,ksr,i8)
      integer olab8
_ENDIF
INCLUDE(common/sizes)
INCLUDE(common/prints)
      common/linkmr/ icex,i2ex,i3ex
      common/scra /res(1)
      common/bufd/ gout(510),nword
      common/scrtch/ndet(5),nsac(5),iaz(5),iez(5),idra(5),
     +idrc(5),jdrc(5),
     +jan(7),jbn(7),ispa2(7),ae(7829),ot(44000),ideks(maxorb),
     +olab(2000),ij(8),nytl(5),nplu(5),ndub(5),nod(5),jkan(mxcrec),
     +nir(maxorb),loc(maxorb),sac(10),
_IF(cray,ksr,i8)
     + olab8(250),f(2304),h(mxcrec),isc(5),bob(3),itest(mxnshl),
_ELSE
     + olab8(250),f(2304),h(mxcrec),isc(5),isd,bob(3),itest(mxnshl),
_ENDIF
     + core,
     + w(mxcsf*(mxcsf+1)/2),
     + ew(mxroot),vect(mxcsf*mxroot),trsum(10,mxroot),
     + wect(mxcsf*mxroot),senk(500),iqr(mxcrec),istm(10)
     +,ht(mxroot),hs(mxcsf),ea(mxroot),iper(10),
     + coff(mxcsf*mxroot),doff(mxcsf),nrec(5),jtrk(mxcsf),
     + isym(mxroot+1),mkon(mxcrec),ikan(mxcrec),hy(mxcsf*mxcsf)
     +,ktrk(mxcsf)
      common/jany/ jerk(10),jbnk(10),idum11(10)
      dimension bkay(intmax)
      dimension hp(4800)
      common /b/ iwod,nid,lg,irsf,ntape,ntype,ical,ntab,nc,lt,lz,nx,
     + nd,nston,itape,jtape,mdisk,ideli,mtape,iswh,nr,mm,jblk,jto,
     + igmax,nr1,mtype,linf,kfile,nawk,isec,nshl,
     + nsec,lulu,nrootx,lsng,nprin,ipt0,
     + mconf(5),nconf(5),mh,jswh,
     + id(mxcsf),mn(mxcsf),nsc(mxcsf),
     + imo,lc,mz,ndx,nrx,ntx,ny,nt,iht,ihz,nz,mt,nzx,jy,
     + jx,nd1,nd2,md1,md2,nr3,nr4,mr1,mr2,nt1,nt2,mt1,mt2,ih1,ih2,ih3,
     + ih4,jh2,jh3,jh4,iab,nt3,iac,nr2,m,nad,nad1,ifr1, mx,md,lr,ipag
       common /d/ jbab,kml,ifl,iiz,nis,ii,lzer,km2,ndl,jsec,ie,iej,
     * la,kc,kmk,ndk,kr,iiq,ic1,ic2,i21,i22,i31,i32,iix,
     * kzer,ndt,iek,j3b,j3a,lps,nps,lms,nms,jcb,jca,lcl,ndi,kmi,
     * ia,jai,nzei,nisi,mc,kk,ndj,kmj,ja,nzer,mr,mps,mms,
     * ic,i2,i3,jc,j3
      equivalence (hp(1),hy(1))
c
      if (oprint(31)) then
       write(6,*) 'mc,kmj,jbab,mm',mc,kmj,jbab,mm
1860   format(1x,14i5)
      endif
      do 29 l=1,mc
      kbab=jbab
      jbab=jbab+kmj
      if (ifl.eq.1) go to 400
      ifl=1
      go to 30
400   mm=mm+1
      kk=olab(mm)
      if (mm.lt.nid) go to 30
      mm=0
      read(ntape) olab8
      call unpack(olab8,8,olab,2000)
30    go to (405,31,32,33), kk
31    mm=mm+1
      ll=olab(mm)
      if (ll.gt.128) ll=ll-256
      if (mm.lt.nid) go to 34
      mm=0
      read(ntape) olab8
      call unpack(olab8,8,olab,2000)
34    mm=mm+1
      jj=olab(mm)
      if (mm.lt.nid) go to 35
      mm=0
      read(ntape) olab8
      call unpack(olab8,8,olab,2000)
35    mm=mm+1
      kk=olab(mm)
      if (mm.lt.nid) go to 36
      mm=0
      read(ntape) olab8
      call unpack(olab8,8,olab,2000)
36    jj=(kk-1)*mr +jj
      jj=jj*ii+ic
      ip=ot(jj)
      if(ip.eq.255) ll=-ll
      jto=jto+1
      if (ll.lt.0) go to 820
      coul = h(jto)
      if (jto.lt.iwod) go to 37
      jto=0
      read(mtype) h
37    jto=jto+1
      exc=-h(jto)
      if (jto.lt.iwod) go to 38
      jto=0
      read(mtype) h
      go to 38
820   coul=-h(jto)
      if (jto.lt.iwod) go to 821
      jto=0
      read(mtype) h
821   jto=jto+1
      ll=-ll
      exc=h(jto)
      if (jto.lt.iwod) go to 38
      jto=0
      read(mtype) h
38    if (ll-2)810,811,812
810   bob(1)=-exc
      bob(2)=-coul
      bob(3)=coul+exc
      go to 49
811   bob(1)=exc
      bob(2)=-coul-exc
      bob(3)=coul
      go to 49
812   bob(1)=coul+exc
      bob(2)=-exc
      bob(3)=-coul
 49   continue
_IF1()c     do 39 m=1,nzer
_IF1()c39   f(m)=0
      call vclr(f,1,nzer)
      do 40 m=1,kml
      jj=jj+1
      my=ot(jj)
      jj=jj+1
      if(my.eq.0) go to 40
_IF1(iv)      kx=ja+my
_IFN1(iv)      kx=ja+my+ndj
      mike=ot(jj)
      if(mike.gt.128) go to 42
      tim=bob(mike)
      go to 43
42    tim=-bob(256-mike)
 43   continue
_IF(ibm)
      lx=m-kml
      do 44 if=1,kmj
      kx=kx+ndj
      lx=lx+kml
44    f(lx)=f(lx)+tim*ae(kx)
_ELSE
      call daxpy(kmj,tim,ae(kx),ndj,f(m),kml)
_ENDIF
40    continue
47    mx=ie
_IF(cray)
      call mxma(ae(ie+1),kml,1,f,1,kml,hp(kbab+1),jsec,1,kml,kml,kmj)
_ELSE
      call mxmaa(ae(ie+1),kml,1,f,1,kml,hp(kbab+1),jsec,1,kml,kml,kmj)
_ENDIF
      go to 29
405   ig7=kbab-jsec
      do 406 m=1,kml
      ig7=ig7+jsec
      in3=ig7
      do 406 if=1,kmj
      in3=in3+1
 406  hp(in3)=0
      go to 29
  32  mm=mm+1
      ll=olab(mm)
      if (mm.lt.nid) go to 51
      mm=0
      read(ntape)olab8
      call unpack(olab8,8,olab,2000)
51    jto=jto+1
      sm=h(jto)
      if(jto.lt.iwod) go to 52
      jto=0
      read(mtype) h
52    if(ll.lt.128) go to 53
      sm=-sm
      ll=256-ll
53    continue
      if (oprint(31)) then
       write(6,5655)kml,ll,i2
5655   format('sm1313 - 53: kml,ll,i2 =', i3, 2i10)
      endif
      jj=ll*kml+i2
_IF1()c     do 55 m=1,nzer
_IF1()c55   f(m)=0
      call vclr(f,1,nzer)
      if (oprint(31)) then
       write(6,5656) kml,jj
5656   format(1x,'sm1313 - do 56:kml,jj = ',i3,2x,i4)
       write(6,5657) (ot(loop+jj),loop=1,kml)
5657   format(1x,'sm1313 - do 56:ot     = ',10i4)
      endif
      do 56 m=1,kml
      jj=jj+1
      my=ot(jj)
      if (my.eq.0) go to 56
      if(my.lt.128) go to 58
_IF1(iv)      kx=ja-my+256
_IFN1(iv)      kx=ja-my+256+ndj
      tim=-sm
      go to 59
58    tim=sm
_IF1(iv)      kx=my+ja
_IFN1(iv)      kx=my+ja+ndj
 59   continue
_IF1(iv)      lx=m-kml
_IF1(iv)      do 60 if=1,kmj
_IF1(iv)      kx=kx+ndj
_IF1(iv)      lx=lx+kml
_IF1(iv) 60   f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv)      call daxpy(kmj,tim,ae(kx),ndj,f(m),kml)
56    continue
      go to 47
  33  mm=mm+1
      ll=olab(mm)
      if (mm.lt.nid) go to 61
      mm=0
      read(ntape)olab8
      call unpack(olab8,8,olab,2000)
61    mm=mm+1
      ntar=olab(mm)
      if (mm.lt.nid) go to 62
      mm=0
      read (ntape) olab8
      call unpack(olab8,8,olab,2000)
62    mm=mm+1
      nb=olab(mm)
      if(mm.lt.nid) go to 63
      mm=0
      read (ntape) olab8
      call unpack(olab8,8,olab,2000)
63    mm=mm+1
      mb=olab(mm)
      if (mm.lt.nid) go to 41
      mm=0
      read (ntape) olab8
      call unpack(olab8,8,olab,2000)
  41  nb=ideks(nb)+mb+ij(ntar)
      sm=bkay(nb)
      jto=jto+1
      sm=sm+h(jto)
      if(jto.lt.iwod) go to 64
      jto=0
      read(mtype) h
64    if (mr.eq.0 ) go to 65
      do 66 m=1,mr
      jto=jto+1
      sm=sm+h(jto)
      if(jto.lt.iwod) go to 68
      jto=0
      read(mtype) h
68    jto=jto+1
      sac(m)=-h(jto)
      if (jto.lt.iwod) go to 66
      jto=0
      read(mtype) h
66    continue
65    if(nd.eq.0) go to 67
      do 69 m=1,nd
      jto=jto+1
      tim=h(jto)
      sm=sm+tim+tim
      if(jto.lt.iwod) go to 70
      jto=0
      read(mtype) h
70    jto=jto+1
      sm=sm-h(jto)
      if(jto.lt.iwod) go to 69
      jto=0
      read(mtype) h
69    continue
67    if(ll.gt.128) go to 86
      jj=ll*j3+i3
      ip=ot(jj)
      if(ip.eq.1) go to 71
      sm=-sm
      if(mr.eq.0) go to 71
      do 72 m=1,mr
72    sac(m)=-sac(m)
71    jj=jj+1
_IF1()c     do 80 m=1,nzer
_IF1()c80   f(m)=0
      call vclr(f,1,nzer)
_IF1(iv)      jx=-kml
      do 73 m=1,kml
_IF1(iv)      jx=jx+1
      in=jj
      im=ot(in)
      go to (74,75,76,77),im
 74   in=in+1
      kx=ot(in)
_IF1(iv)      kx=kx+ja
_IFN1(iv)      kx=kx+ja+ndj
      tim=sm
      if(mps.eq.0) go to 78
_IF1(iv)      do 79 if=1,mps
_IF1(iv)       in=in+1
_IF1(iv)      im=ot(in)
_IF1(iv) 79   tim=tim+sac(im)
_IF1(c)      call gather(mps,res,sac,ot(in+1))
_IFN1(civ)      call dgthr(mps,sac,res,ot(in+1))
_IFN1(iv)      tim=tim+dsum(mps,res,1)
78    continue
_IF1(iv)       lx=jx
_IF1(iv)       do 81 if=1,kmj
_IF1(iv)       kx=kx+ndj
_IF1(iv)       lx=lx+kml
_IF1(iv)  81   f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv)      call daxpy(kmj,tim,ae(kx),ndj,f(m),kml)
      go to 73
 75   in=in+1
      kx=ot(in)
_IF1(iv)      kx=kx+ja
_IFN1(iv)      kx=kx+ja+ndj
      tim=-sm
      if(mms.eq.0) go to 78
_IF1(iv)      do 83 if=1,mms
_IF1(iv)      in=in+1
_IF1(iv)      im=ot(in)
_IF1(iv) 83   tim=tim-sac(im)
_IF1(c)      call gather(mms,res,sac,ot(in+1))
_IFN1(civ)      call dgthr(mms,sac,res,ot(in+1))
_IFN1(iv)      tim=tim-dsum(mms,res,1)
      go to 78
76    do 84 if=1,nms
      in=in+1
      kx=ot(in)
_IF1(iv)      kx=kx+ja
_IFN1(iv)      kx=kx+ja+ndj
      in=in+1
      im=ot(in)
      tim=-sac(im)
_IF1(iv)      lx=jx
_IF1(iv)      do 84 ig=1,kmj
_IF1(iv)      kx=kx+ndj
_IF1(iv)      lx=lx+kml
_IF1(iv) 84   f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv) 84   call daxpy(kmj,tim,ae(kx),ndj,f(m),kml)
      go to 73
77    do 85 if=1,nps
      in=in+1
      kx=ot(in)
_IF1(iv)      kx=kx+ja
_IFN1(iv)      kx=kx+ja+ndj
      in=in+1
      im=ot(in)
      tim=sac(im)
_IF1(iv)      lx=jx
_IF1(iv)      do 85 ig=1,kmj
_IF1(iv)      kx=kx+ndj
_IF1(iv)      lx=lx+kml
_IF1(iv) 85   f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv) 85   call daxpy(kmj,tim,ae(kx),ndj,f(m),kml)
73    jj=jj+jc
      go to 47
86    jj=(256-ll)*j3+i3
      ip=ot(jj)
      if (ip.eq.1) go to 87
      sm=-sm
      if (mr.eq.0) go to 87
      do 88 m=1,mr
88    sac(m)=-sac(m)
87    jj=jj+1
_IF1()c     do 89 m=1,nzer
_IF1()c89   f(m)=0
      call vclr(f,1,nzer)
_IF1(iv)      jx=-kml
      do 90 m=1,kml
_IF1(iv)      jx=jx+1
      in=jj
      im=ot(in)
      go to (91,95,96,97),im
91    in=in+1
      kx=ot(in)
_IF1(iv)      kx=kx+ja
_IFN1(iv)      kx=kx+ja+ndj
      tim=sm
      if(mms.eq.0) go to 92
_IF1(iv)      in=in+mps
_IF1(iv)      do 93 if=1,mms
_IF1(iv)      in=in+1
_IF1(iv)      im=ot(in)
_IF1(iv) 93   tim=tim+sac(im)
_IF1(c)      call gather(mms,res,sac,ot(in+mps+1))
_IFN1(civ)      call dgthr(mms,sac,res,ot(in+mps+1))
_IFN1(iv)      tim=tim+dsum(mms,res,1)
92    continue
_IF1(iv)      lx=jx
_IF1(iv)      do 94 if=1,kmj
_IF1(iv)      kx=kx+ndj
_IF1(iv)      lx=lx+kml
_IF1(iv) 94   f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv)      call daxpy(kmj,tim,ae(kx),ndj,f(m),kml)
      go to 90
95    in=in+1
      kx=ot(in)
_IF1(iv)      kx=kx+ja
_IFN1(iv)      kx=kx+ja+ndj
      tim=-sm
      if(mps.eq.0) go to 92
_IF1(iv)      in=in+mms
_IF1(iv)      do 99 if=1,mps
_IF1(iv)      in=in+1
_IF1(iv)      im=ot(in)
_IF1(iv) 99   tim=tim-sac(im)
_IF1(c)      call gather(mps,res,sac,ot(in+mms+1))
_IFN1(civ)      call dgthr(mps,sac,res,ot(in+mms+1))
_IFN1(iv)      tim=tim-dsum(mps,res,1)
      go to 92
96    do 100 if=1,nms
      in=in+1
      kx=ot(in)
_IF1(iv)      kx=kx+ja
_IFN1(iv)      kx=kx+ja+ndj
      in=in+1
      im=ot(in)
      tim=sac(im)
_IF1(iv)      lx=jx
_IF1(iv)      do 100 ig=1,kmj
_IF1(iv)      kx=kx+ndj
_IF1(iv)      lx=lx+kml
_IF1(iv) 100  f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv) 100  call daxpy(kmj,tim,ae(kx),ndj,f(m),kml)
      go to 90
97    do 101 if=1,nps
      in=in+1
      kx=ot(in)
_IF1(iv)      kx=kx+ja
_IFN1(iv)      kx=kx+ja+ndj
      in=in+1
      im=ot(in)
      tim=-sac(im)
_IF1(iv)      lx=jx
_IF1(iv)      do 101 ig=1,kmj
_IF1(iv)      kx=kx+ndj
_IF1(iv)      lx=lx+kml
_IF1(iv) 101  f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv) 101  call daxpy(kmj,tim,ae(kx),ndj,f(m),kml)
90    jj=jj+jc
      go to 47
29    continue
      return
      end
      subroutine smr112(iot,ideks)
      implicit REAL (a-h,o-z), integer(i-n)
      integer olab
_IF(cray,ksr,i8)
      integer olab8
_ENDIF
INCLUDE(common/sizes)
      common /b/ iwod,nid,lg,irsf,ntape,ntype,ical,ntab,nc,lt,lz,ntx,
     + ndx,nston,itape,jtape,mdisk,ideli,mtape,iswh,nrx,mm,jblk,jto,
     + igmax,nr1,ltype,linf,kfile,
     + nawk,isec,nshl,
     + nsec,lulu,nrootx,lsng,nprin,ipt0,mconf(5),nconf(5),mh,jswh,
     + id(mxcsf),ie(mxcsf),nsc(mxcsf),
     + imo,mc,jh4,nd,nr,nx,nt,ny,iht,ihz,
     + nzx,ih4,nz,jy,jx,nd1,md,md1,md2,nr3,mr,mr1,mr2,nt1,mx,
     + mt1,mt2,ih1,ih2,ih3,mt,jh2,jh3,mz,iab,nt3,iac,nr2,m,
     + nad,nad1,ifr1,nt2,nd2,nr4,idum1(1)
      common/scrtch/jdeks(9),kdeks(8),
     + lj(8),olab(2000),nit(666),nytl(5),ndub(5),
     + icon(5),nod(5),jdum(maxorb),nir(maxorb),loc(maxorb)
      common/junk/ican(3),icbn,olab8(250),jcon(8200),mkon(mxcrec),
     + kc(mxcrec),kd(mxcrec),lab(3)
     + ,ij(8),jtest(mxnshl),nyzl(mxcsf),jkan(mxcrec),nplu(5),
     + jkon(mxnshl*mxcsf)
      dimension iot(*),ideks(*)
      l3=mz
      l4=mt
      do 17 k=1,mc
      jz=l3
      do 19 l=1,md
      jz=jz+1
      ll=iot(jz)
      if (jcon(ll)-1) 20,21,19
   21 ia=l+1
      go to 22
   19 continue
   20 mm=mm+1
      olab(mm)=0
      if (mm.lt.nid) go to 61
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      mm=0
      go to 61
   22 do 23 l=ia,md
      jz=jz+1
      nn=iot(jz)
      if (jcon(nn)-1) 20,24,23
   24 ja=l+1
      go to 25
   23 continue
   25 if (ja.gt.md) go to 26
      do 27 l=ja,md
      jz=jz+1
      ii=iot(jz)
      if (jcon(ii).ne.2) go to 20
   27 continue
   26 if (mr.eq.0) go to 28
      jz=l4
      do 29 l=1,mr
      jz=jz+1
      ii=iot(jz)
      if (jcon(ii).ne.1) go to 20
   29 continue
   28 jz=nt
      kz=l4+1
      jq=0
      if (mr.gt.0) iq=iot(kz)
      do 30 l=1,nr
      jz=jz+1
      ii=iot(jz)
      if (ii.eq.ll) go to 31
      if (jq.eq.mr.or.iq.ne.ii) go to 32
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   30 continue
   31 ip=l
      ia=l+1
      do 33 l=ia,nr
      jz=jz+1
      ii=iot(jz)
      if (ii.eq.nn) go to 34
      if (jq.eq.mr.or.iq.ne.ii) go to 35
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   33 continue
   32 ip=l
      ia=l+1
      do 36 l=ia,nr
      jz=jz+1
      jj=iot(jz)
      if (jj.eq.ll) go to 37
      if (jq.eq.mr.or.iq.ne.jj) go to 38
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   36 continue
   34 jp=ideks(l-1)
      ia=l+1
      do 39 l=ia,nr
      jz=jz+1
      ii=iot(jz)
      if (jq.eq.mr.or.iq.ne.ii) go to 40
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   39 continue
   40 kp=jdeks(l-2)
      ia=l+1
      do 41 l=ia,nr
      jz=jz+1
      jj=iot(jz)
      if (jq.eq.mr.or.iq.ne.jj) go to 42
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   41 continue
   42 mm=mm+1
      olab(mm)=1
      lp=kdeks(l-3)
      if (mm.lt.nid) go to 43
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      mm=0
      go to 43
   35 jp=ideks(l-1)
      ia=l+1
      do 44l=ia,nr
      jz=jz+1
      jj=iot(jz)
      if (jj.eq.nn) go to 45
      if (jq.eq.mr.or.jj.ne.iq) go to 46
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   44 continue
   37 jp=ideks(l-1)
      ia=l+1
      do 47 l=ia,nr
      jz=jz+1
      jj=iot(jz)
      if (jj.eq.nn) go to 48
      if (jq.eq.mr.or.jj.ne.iq) go to 49
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   47 continue
   38 jp=ideks(l-1)
      ia=l+1
      do 50 l=ia,nr
      jz=jz+1
      kk=iot(jz)
      if (kk.eq.ll) go to 51
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   50 continue
   51 kp=jdeks(l-2)
      ia=l+1
      do 52 l=ia,nr
      jz=jz+1
      kk=iot(jz)
      if (kk.eq.nn) go to 42
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   52 continue
   48 kp=jdeks(l-2)
      ia=l+1
      do 53 l=ia,nr
      jz=jz+1
      jj=iot(jz)
      if (jq.eq.mr.or.jj.ne.iq) go to 54
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   53 continue
   54 lp=kdeks(l-3)
      mm=mm+1
      olab(mm)=3
      if (mm.lt.nid) go to 43
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      mm=0
      go to 43
   49 kp=jdeks(l-2)
      ia=l+1
      do 55 l=ia,nr
      jz=jz+1
      kk=iot(jz)
      if (kk.eq.nn) go to 56
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   55 continue
   56 lp=kdeks(l-3)
      mm=mm+1
      olab(mm)=2
      if (mm.lt.nid) go to 43
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      mm=0
      go to 43
   45 kp=jdeks(l-2)
      ia=l+1
      do 57 l=ia,nr
      jz=jz+1
      jj=iot(jz)
      if (jq.eq.mr.or.jj.ne.iq) go to 56
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   57 continue
   46 kp=jdeks(l-2)
      ia=l+1
      do 58 l=ia,nr
      jz=jz+1
      kk=iot(jz)
      if (kk.eq.nn) go to 54
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   58 continue
 43   mm=mm+1
      olab(mm)=ip+jp+kp+lp
      if (mm.lt.nid) go to 170
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      mm=0
 170  nt1r = nir(ll)
      nl1r = loc (ll)
      nt2r = nir (nn)
      nl2r = loc (nn)
      kix=3
      nb1r = nir(ii)
      nm1r = loc (ii)
      nb2r = nir (jj)
      nm2r = loc (jj)
  120 if (nt1r-nb1r)123,122,199
  199 iax = ideks(nt1r)+nb1r
      if (nt2r.lt.nb2r) go to 124
      ibx = ideks(nt2r)+nb2r
      if (iax-ibx) 125,126,127
  127 icx = ideks(iax) + ibx
      icq = (nl1r-1)*lj(nb1r) + nm1r
  140 ljn = lj(nb2r)
      idq = (nl2r-1)*ljn + nm2r
      icq = (icq-1)*ljn*lj(nt2r)+idq
      go to 129
  126 icx = ideks (iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn+nm1r
  141 idq = (nl2r-1)*ljn+nm2r
      go to 133
  125 icx = ideks(ibx) + iax
      icq = (nl2r-1)*lj(nb2r)+nm2r
  135 ljn = lj(nb1r)
      idq = (nl1r-1)*ljn+nm1r
      icq =(icq-1)*ljn*lj(nt1r) +idq
      go to 129
  124 ibx = ideks(nb2r) + nt2r
      if (iax-ibx) 130,131,132
  132 icx = ideks(iax)+ibx
      icq = (nl1r-1)*lj(nb1r)+nm1r
      ljn = lj(nt2r)
      idq = (nm2r-1)*ljn+nl2r
      icq = (icq-1)*ljn*lj(nb2r) + idq
      go to 129
  131 icx = ideks(iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn+nm1r
      idq = (nm2r-1)*ljn+nl2r
  133 if (icq.lt.idq) go to 134
      icq = ideks(icq)+idq
      go to 129
  134 icq = ideks(idq)+icq
      go to 129
  130 icx=ideks(ibx)+iax
      icq = (nm2r-1)*lj(nt2r)+nl2r
      go to 135
  123 iax = ideks(nb1r)+nt1r
      if (nt2r.lt.nb2r) go to 136
      ibx = ideks(nt2r)+nb2r
      if (iax-ibx) 137,138,139
  139 icx = ideks(iax)+ibx
      icq = (nm1r-1)*lj(nt1r)+nl1r
      go to 140
  138 icx = ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      go to 141
  137 icx = ideks(ibx) + iax
      icq = (nl2r-1)*lj(nb2r) +nm2r
  145 ljn = lj(nt1r)
      idq = (nm1r-1)*ljn         +nl1r
      icq = (icq-1)*ljn*lj(nb1r) +idq
      go to 129
  136 ibx = ideks(nb2r)+nt2r
      if (iax-ibx) 142,143,144
  144 icx = ideks(iax)+ibx
      icq =(nm1r-1) *lj(nt1r)+nl1r
      ljn = lj(nt2r)
      idq = (nm2r-1)*ljn+nl2r
      icq = (icq-1)*ljn*lj(nb2r)+idq
      go to 129
  143 icx=ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      idq = (nm2r-1)*ljn+nl2r
      go to 133
  142 icx = ideks(ibx)+iax
      icq =(nm2r-1)*lj(nt2r)+nl2r
      go to 145
  122 iax=ideks(nt1r+1)
      if (nl1r.lt.nm1r) go to 146
      iay=ideks(nl1r)+nm1r
      go to 147
  146 iay= ideks(nm1r)+nl1r
  147 if (nl2r.lt.nm2r) go to 148
      iby=ideks(nl2r)+nm2r
      go to 149
  148 iby=ideks(nm2r)+nl2r
  149 if (nt1r.eq.nt2r) go to 150
      ibx=ideks(nt2r+1)
      if (iax.lt.ibx) go to 151
      icx=ideks(iax)+ibx
      ljn =lj(nt2r)
      icq = (iay-1)*ideks(ljn+1)+iby
      go to 129
  151 icx=ideks(ibx)+iax
      ljn=lj(nt1r)
      icq =(iby-1)*ideks(ljn+1)+iay
      go to 129
  150 icx=ideks(iax+1)
      if (iay.lt.iby) go to 153
      icq=ideks(iay)+iby
      go to 129
  153 icq=ideks(iby)+iay
 129  icq=nit(icx)  +  icq
      idud=(icq -1)/igmax
      jto=jto + 1
      kc(jto) = idud  + 1
      icq=icq  -  idud*igmax
      kd(jto)  = icq
      if(jto.lt.iwod)  go to 160
      write(ntype)  kc,kd
      jto=0
      jblk=jblk  + 1
  160 if (kix.lt.0) go to 61
      kix=-3
      itr=nb1r
      nb1r=nb2r
      nb2r=itr
      itr=nm1r
      nm1r=nm2r
      nm2r=itr
      go to 120
   61 l3=l3+mx
   17 l4=l4+mx
      return
      end
      subroutine smr114(iot,ideks)
      implicit REAL (a-h,o-z), integer(i-n)
      integer olab
_IF(cray,ksr,i8)
      integer olab8
_ENDIF
INCLUDE(common/sizes)
      common /b/ iwod,nid,lg,irsf,ntape,ntype,ical,ntab,nc,lt,lz,ntx,
     + ndx,nston,itape,jtape,mdisk,ideli,mtape,iswh,nrx,mm,jblk,jto,
     + igmax,nr2,ltype,linf,kfile,
     + nawk,isec,nshl,
     + nsec,lulu,nrootx,lsng,nprin,ipt0,mconf(5),nconf(5),mh,jswh,
     + id(mxcsf),ie(mxcsf),nsc(mxcsf),
     + imo,mc,jh3,nd,nr,nx,nt,ny,iht,ihz,
     + nzx,ih3,nz,jy,jx,md,nd2,md1,md2,mr,nr4,mr1,mr2,nt1,nt2,
     + mt1,mt2,ih1,ih2,mt,ih4,jh2,mz,jh4,iac,mx,iab,nr1,m,
     + nad,nad1,ifr1,nt3,nd1,nr3,if
      common/scrtch/jdeks(9),kdeks(8),
     + lj(8),olab(2000),nit(666),nytl(5),ndub(5),
     + icon(5),nod(5),jdum(maxorb),nir(maxorb),loc(maxorb)
      common/junk/ican(3),icbn,olab8(250),jcon(8200),mkon(mxcrec),
     + kc(mxcrec),kd(mxcrec),lab(3)
     + ,ij(8),jtest(mxnshl),nyzl(mxcsf),jkan(mxcrec),nplu(5),
     + jkon(mxnshl*mxcsf)
      dimension iot(*),ideks(*)
      l3=mz
      l4=mt
      do 69 k=1,mc
      jz=l3
      do 70 l=1,md
      jz=jz+1
      jj=iot(jz)
      if (jcon(jj)-1) 71,72,70
   70 continue
   71 if (l.eq.md) go to 73
      ia=l+1
      do 74 l=ia,md
      jz=jz+1
      ll=iot(jz)
      if (jcon(ll).ne.2) go to 75
   74 continue
      go to 73
  75  mm=mm+1
      olab(mm)=1
      if (mm.lt.nid) go to 260
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      mm=0
      go to 260
   72 l1=l
      if (l.eq.md) go to 76
      ia=l+1
      do 78 l=ia,md
      jz=jz+1
      ll=iot(jz)
      if (jcon(ll)-1) 75,77,78
  78  continue
      go to 76
  77  l2=l
      if (l.eq.md) go to 79
      ia=l+1
      do 80 l=ia,md
      jz=jz+1
      kk=iot(jz)
      if (jcon(kk).ne.2) go to 75
   80 continue
      go to 79
   73 if (mr.eq.0) go to 81
      jz=l4
      do 82 l=1,mr
      jz=jz+1
      ll=iot(jz)
      if (jcon(ll).ne.1) go to 75
   82 continue
  81  mm=mm+1
      olab(mm)=3
      if (mm.lt.nid) go to 83
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      mm=0
  83  jz=nt
      kz=l4+1
      jq=0
      if (mr.gt.0) iq=iot(kz)
      do 84 l=1,nr
      jz=jz+1
      ll=iot(jz)
      if (jq.eq.mr.or.iq.ne.ll) go to 85
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   84 continue
   85 ia=l+1
      ip=l
      do 86 l=ia,nr
      jz=jz+1
      nn=iot(jz)
      if (jq.eq.mr.or.iq.ne.nn) go to 87
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   86 continue
  87  mm=mm+1
      olab(mm)=ideks(l-1)+ip
      if (mm.lt.nid) go to 88
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      mm=0
  88  nt1r=nir(jj)
      nl1r=loc(jj)
      nb1r=nir(ll)
      nm1r=loc(ll)
      nm2r=loc(nn)
  220 if (nt1r-nb1r) 223,222,299
  299 iax = ideks(nt1r)+nb1r
      icx = ideks (iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn
      idq = icq+nm2r
      icq=ideks(idq)+icq+nm1r
      go to 229
  223 iax = ideks(nb1r)+nt1r
      icx=ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      idq = (nm2r-1)*ljn+nl1r
      icq=ideks(idq)+icq
      go to 229
  222 iax=ideks(nt1r+1)
      icx=ideks(iax+1)
      if (nl1r.lt.nm1r) go to 246
      iay=ideks(nl1r)+nm1r
      go to 247
  246 iay= ideks(nm1r)+nl1r
      iby=ideks(nm2r)+nl1r
      go to 253
  247 if (nl1r.lt.nm2r) go to 248
      iby=ideks(nl1r)+nm2r
      go to 253
  248 iby=ideks(nm2r)+nl1r
  253 icq=ideks(iby)+iay
  229 icq=nit(icx)  +  icq
      idud=(icq -1)/igmax
      jto=jto + 1
      kc(jto) = idud  + 1
      icq=icq  -  idud*igmax
      kd(jto)  = icq
      if(jto.lt.iwod)  go to 260
      write(ntype) kc,kd
      jto=0
      jblk=jblk  + 1
      go to 260
   79 if (mr.eq.0) go to 89
      jz=l4
      do 90 l=1,mr
      jz=jz+1
      ii=iot(jz)
      if (jcon(ii)-1) 75,90,91
  90  continue
      go to 89
  91  l1=l
      if (l.eq.mr) go to 93
      ia=l+1
      do 92 l=ia,mr
      jz=jz+1
      jc=iot(jz)
      if (jcon(jc).ne.1) go to 75
  92  continue
      go to 93
  89  mm=mm+1
      olab(mm)=3
      if (mm.lt.nid) go to 94
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      mm=0
   94 jz=nt
      do 95 l=1,nr
      jz=jz+1
      if (iot(jz).eq.jj) go to 96
   95 continue
   96 ia=l+1
      ip=l
      do 97 l=ia,nr
      jz=jz+1
      if (iot(jz).eq.ll) go to 98
   97 continue
  98  mm=mm+1
      olab(mm)=-ip-ideks(l-1)
      if (mm.lt.nid) go to 99
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      mm=0
   99 jz=nz
      if (l1.eq.1) go to 180
      kz=l3
      ja=l1-1
      do 181 l=1,ja
      kz=kz+1
      jz=jz+1
      kk=iot(jz)
      if (kk.ne.iot(kz)) go to 182
  181 continue
  180 if (l2.eq.l1+1) go to 183
      kz=l3+l1
      ia=l1+1
      ja=l2-1
      do 184 l=ia,ja
      kz=kz+1
      jz=jz+1
      kk=iot(jz)
      if (kk.ne.iot(kz)) go to 182
  184 continue
  183 if (l2.lt.md) go to 185
      kk=iot(jz+1)
      go to 182
  185 kz=l3+l2
      ia=l2+1
      do 186 l=ia,md
      kz=kz+1
      jz=jz+1
      kk=iot(jz)
      if (kk.ne.iot(kz)) go to 182
  186 continue
      kk=iot(jz+1)
  182 nt1r=nir(kk)
      nl1r=loc(kk)
      nb1r=nir(jj)
      nm1r=loc(jj)
      nm2r=loc(ll)
      go to 220
  93  mm=mm+1
      olab(mm)=2
      if (mm.lt.nid) go to 187
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      mm=0
 187  jz=nt
      kz=l4
      ig=0
      if (l1.eq.1) go to 188
      ja=l1-1
      ka=1
      do 189 l=ka,ja
      kz=kz+1
      jp=iot(kz)
 802  jz=jz+1
      if (jp.eq.iot(jz)) go to 189
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.3) go to 190
      go to 802
  189 continue
  188 if (l1.eq.mr) go to 191
      ia=l1+1
      kz=kz+1
      do 192 l=ia,mr
      kz=kz+1
      jp=iot(kz)
 803  jz=jz+1
      if (jp.eq.iot(jz)) go to 192
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.3) go to 190
      go to 803
  192 continue
  191 lab(3)=nr
      if (ig-1) 284,285,190
  284 lab(2)=nr1
      lab(1)=mr
      go to 190
  285 lab(2)=nr1
  190 jz=lab(1)+nt
      nn=iot(jz)
      mm=mm+1
      if (nn.eq.jj) go to 194
      olab(mm)=1
  195 if (mm.lt.nid) go to 196
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  196 mm=mm+1
      olab(mm)=l1
      if (mm.lt.nid) go to 197
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  197 mm=mm+1
      jp=lab(2)-1
      kp=lab(3)-2
      olab(mm)=lab(1)+ideks(jp)+jdeks(kp)
      if (mm.lt.nid) go to 198
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      go to 198
  194 jz=lab(2)+nt
      nn=iot(jz)
      if (nn.eq.ll) go to 271
      olab(mm)=2
      go to 195
  271 jz=lab(3)+nt
      nn=iot(jz)
      olab(mm)=3
      go to 195
  198 nt1r=nir(jj)
      nt2r=nir(ll)
      nl1r=loc(jj)
      nl2r=loc(ll)
      kix=3
      nb1r=nir(nn)
      nb2r=nir(ii)
      nm1r=loc(nn)
      nm2r=loc(ii)
  320 if (nt1r-nb1r) 323,322,399
  399 iax = ideks(nt1r)+nb1r
      if (nt2r.lt.nb2r) go to 324
      ibx = ideks(nt2r)+nb2r
      if (iax-ibx) 325,326,327
  327 icx = ideks(iax) + ibx
      icq = (nl1r-1)*lj(nb1r) + nm1r
  340 ljn = lj(nb2r)
      idq = (nl2r-1)*ljn + nm2r
      icq = (icq-1)*ljn*lj(nt2r)+idq
      go to 329
  326 icx = ideks (iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn+nm1r
  341 idq = (nl2r-1)*ljn+nm2r
      go to 333
  325 icx = ideks(ibx) + iax
      icq = (nl2r-1)*lj(nb2r)+nm2r
  335 ljn = lj(nb1r)
      idq = (nl1r-1)*ljn+nm1r
      icq =(icq-1)*ljn*lj(nt1r) +idq
      go to 329
  324 ibx = ideks(nb2r) + nt2r
      if (iax-ibx) 330,331,332
  332 icx = ideks(iax)+ibx
      icq = (nl1r-1)*lj(nb1r)+nm1r
      ljn = lj(nt2r)
      idq = (nm2r-1)*ljn+nl2r
      icq = (icq-1)*ljn*lj(nb2r) + idq
      go to 329
  331 icx = ideks(iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn+nm1r
      idq = (nm2r-1)*ljn+nl2r
  333 if (icq.lt.idq) go to 334
      icq = ideks(icq)+idq
      go to 329
  334 icq = ideks(idq)+icq
      go to 329
  330 icx=ideks(ibx)+iax
      icq = (nm2r-1)*lj(nt2r)+nl2r
      go to 335
  323 iax = ideks(nb1r)+nt1r
      if (nt2r.lt.nb2r) go to 336
      ibx = ideks(nt2r)+nb2r
      if (iax-ibx) 337,338,339
  339 icx = ideks(iax)+ibx
      icq = (nm1r-1)*lj(nt1r)+nl1r
      go to 340
  338 icx = ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      go to 341
  337 icx = ideks(ibx) + iax
      icq = (nl2r-1)*lj(nb2r) +nm2r
  345 ljn = lj(nt1r)
      idq = (nm1r-1)*ljn         +nl1r
      icq = (icq-1)*ljn*lj(nb1r) +idq
      go to 329
  336 ibx = ideks(nb2r)+nt2r
      if (iax-ibx) 342,343,344
  344 icx = ideks(iax)+ibx
      icq =(nm1r-1) *lj(nt1r)+nl1r
      ljn = lj(nt2r)
      idq = (nm2r-1)*ljn+nl2r
      icq = (icq-1)*ljn*lj(nb2r)+idq
      go to 329
  343 icx=ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      idq = (nm2r-1)*ljn+nl2r
      go to 333
  342 icx = ideks(ibx)+iax
      icq =(nm2r-1)*lj(nt2r)+nl2r
      go to 345
  322 iax=ideks(nt1r+1)
      if (nl1r.lt.nm1r) go to 346
      iay=ideks(nl1r)+nm1r
      go to 347
  346 iay= ideks(nm1r)+nl1r
  347 if (nl2r.lt.nm2r) go to 348
      iby=ideks(nl2r)+nm2r
      go to 349
  348 iby=ideks(nm2r)+nl2r
  349 if (nt1r.eq.nt2r) go to 350
      ibx=ideks(nt2r+1)
      if (iax.lt.ibx) go to 351
      icx=ideks(iax)+ibx
      ljn =lj(nt2r)
      icq = (iay-1)*ideks(ljn+1)+iby
      go to 329
  351 icx=ideks(ibx)+iax
      ljn=lj(nt1r)
      icq =(iby-1)*ideks(ljn+1)+iay
      go to 329
  350 icx=ideks(iax+1)
      if (iay.lt.iby) go to 353
      icq=ideks(iay)+iby
      go to 329
  353 icq=ideks(iby)+iay
  329 icq=nit(icx)  +  icq
      idud=(icq -1)/igmax
      jto=jto + 1
      kc(jto) = idud  + 1
      icq=icq  -  idud*igmax
      kd(jto)  = icq
      if(jto.lt.iwod)  go to 360
      write(ntype) kc,kd
      jto=0
      jblk=jblk  + 1
 360  if (kix.lt.0) go to 260
      kix=-3
      itr=nb1r
      nb1r=nb2r
      nb2r=itr
      itr=nm1r
      nm1r=nm2r
      nm2r=itr
      go to 320
  76  if (mr.eq.0) go to 273
      jz=l4
      do 274 l=1,mr
      jz=jz+1
      ll=iot(jz)
      if (jcon(ll)-1) 275,274,75
 274  continue
      go to 273
 275  l1=l
      if (l.eq.mr) go to 277
      ia=l+1
      do 276 l=ia,mr
      jz=jz+1
      nn=iot(jz)
      if (jcon(nn).ne.1) go to 75
  276 continue
 277  mm=mm+1
      olab(mm)=2
      if (mm.lt.nid) go to 278
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      mm=0
 278  jz=nt
      kz=l4
      ig=0
      if (l1.eq.1) go to 279
      ja=l1-1
      do 280 l=1,ja
      kz=kz+1
      jp=iot(kz)
 254  jz=jz+1
      if (jp.eq.iot(jz)) go to 280
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.3) go to 281
      go to 254
 280  continue
 279  if (l1.eq.mr) go to 282
      ia=l1+1
      kz=kz+1
      do 283 l=ia,mr
      kz=kz+1
      jp=iot(kz)
 255  jz=jz+1
      if (jp.eq.iot(jz)) go to 283
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.3) go to 281
      go to 255
  283 continue
  282 lab(3)=nr
      if (ig-1) 286,287,281
  286 lab(2)=nr1
      lab(1)=mr
      go to 281
  287 lab(2)=nr1
  281 jz=lab(1)+nt
      kk=iot(jz)
      mm=mm+1
      if (kk.ne.jj) go to 288
      olab(mm)=-1
      jz=lab(2)+nt
      kk=iot(jz)
      jz=lab(3)+nt
      nn=iot(jz)
 289  if(mm.lt.nid) go to 290
      mm=0
      call pack(olab8,8,olab,2000)
      write(ntape) olab8
 290  mm=mm+1
      olab(mm)=l1
      if(mm.lt.nid) go to 291
      mm=0
      call pack(olab8,8,olab,2000)
      write(ntape) olab8
 291  mm=mm+1
      jp=lab(2)-1
      kp=lab(3)-2
      olab(mm)=lab(1)+ideks(jp)+jdeks(kp)
      if(mm.lt.nid)go to 292
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      go to 292
 288  jz=lab(2)+nt
      nn=iot(jz)
      if (nn.ne.jj) go to 293
      olab(mm)=-2
      jz=lab(3)+nt
      nn=iot(jz)
      go to 289
 293  olab(mm)=-3
      go to 289
 292  nt1r=nir(jj)
      nt2r=nir(ll)
      nl1r=loc(jj)
      nl2r=loc(ll)
      kix=3
      nb1r=nir(kk)
      nb2r=nir(nn)
      nm1r=loc(kk)
      nm2r=loc(nn)
      go to 320
 273  mm=mm+1
      olab(mm)=4
      if(mm.lt.nid) go  to 294
      mm=0
      call pack(olab8,8,olab,2000)
      write(ntape) olab8
 294  jz=nt
      mm=mm+1
      if(mr.eq.0) go to 297
      kz=l4
      do 295 l=1,mr
      kz=kz+1
      jz=jz+1
      ll=iot(jz)
      if(ll.ne.iot(kz)) go to 296
 295  continue
 297  ll=iot(jz+1)
      if(ll.eq.jj) go to 298
      olab(mm)=iab
      if(mm.lt.nid) go to 371
      mm=0
      call pack(olab8,8,olab,2000)
      write(ntape) olab8
      go to 371
 298  ll=iot(jz+2)
      olab(mm)=-iab
      if(mm.lt.nid) go to 371
      mm=0
      call pack(olab8,8,olab,2000)
      write(ntape) olab8
      go to 371
 296  l1 =jz-nt
      if(ll.eq.jj) go to 373
      do 374 l=1,nr
      jz=jz+1
      if(iot(jz).eq.jj) go to 375
 374  continue
 375  olab(mm)=l1+ideks(jz-nt1)
      if(mm.lt.nid) go to 371
      mm=0
      call pack(olab8,8,olab,2000)
      write(ntape) olab8
      go to 371
 373  kz=kz-1
      do 376 ia=l,mr
      kz=kz+1
      jz=jz+1
      ll=iot(jz)
      if(ll.ne.iot(kz)) go to 377
 376  continue
      jz=jz+1
      ll=iot(jz)
 377  olab(mm)=-l1-ideks(jz-nt1)
      if(mm.lt.nid) go to 371
      mm=0
      call pack(olab8,8,olab,2000)
      write(ntape) olab8
 371  mm=mm+1
      iax=nir(jj)
      olab(mm)=iax
      if(mm.lt.nid) go to 372
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  372 iay=ideks(iax+1)
      iaz=ideks(iay+1)
      mq=lj(iax)
      mv=ideks(mq+1)
      iaq=iay-iax
      iaw=iaz-iay
      iat=nit(iaz)
      ii=loc(jj)
      kk=loc(ll)
      km=ideks(kk)
      nm=ideks(ii)
      nn=nm+ii
      if(ii.gt.kk) go to 380
      mm=mm+1
      jb=ii+km
      olab(mm)=kk
      if(mm.lt.nid) go to 378
      mm=0
      call pack(olab8,8,olab,2000)
      write(ntape) olab8
 378  mm=mm+1
      olab(mm)=ii
      if(mm.lt.nid) go to 379
      mm=0
      call pack(olab8,8,olab,2000)
      write(ntape) olab8
 379  ib=ideks(jb) +nn+iat
      go to 381
 380  lb=kk-ii
      jb=nn+lb
      ib=ideks(nn+1)+lb+iat
      mm=mm+1
      olab(mm)=ii
      if(mm.lt.nid) go to 900
      mm=0
      call pack(olab8,8,olab,2000)
      write(ntape) olab8
 900  mm=mm+1
      olab(mm)=kk
      if(mm.lt.nid) go to 381
      mm=0
      call pack(olab8,8,olab,2000)
      write(ntape) olab8
 381  idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if(jto.lt.iwod) go to 382
      jto=0
      write(ntype) kc,kd
      jblk=jblk+1
 382  if(mr.eq.0) go to 383
      jz=l4
      ir=mr
      kix=1
 472  do 384 l=1,ir
      jz=jz+1
      mi=iot(jz)
      mar=nir(mi)
      mlr=loc(mi)
      mlb=ideks(mlr)
      mla=mlb+mlr
      if (mar-iax) 385,395,396
 395  if(mla.lt.jb) go to 386
      ib=iat+ideks(mla) +jb
 388  idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto = jto+1
      kc(jto)=idud+1
      kd(jto) =ib
      if(jto.lt.iwod) go to 387
      jto =0
      write(ntype) kc,kd
      jblk=jblk+1
      go to 387
 386  ib=iat+ideks(jb) +mla
      go to 388
 387  if(mlr.lt.ii) go to 389
      kb=mlb+ii
      go to 390
 389  kb=nm+mlr
 390  if(mlr.lt.kk) go to 391
      lb=mlb+kk
      go to 392
 391  lb=mlr+km
 392  ib=iat+min(kb,lb)+ideks(max(kb,lb))
      idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if(jto.lt.iwod) go to 384
      jto=0
      jblk=jblk+1
      write(ntype) kc,kd
      go to 384
 385  iby=ideks(mar+1)
      iby=iaw+iby
      ibx=iaq+mar
      ibx=ideks(ibx+1)
      ibx=nit(ibx)
      iby=nit(iby)
      ml=lj(mar)
      ib=ideks(ml+1)*(jb-1) +mla+iby
      idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if(jto.lt.iwod) go to 397
      jto=0
      jblk=jblk+1
      write(ntype) kc,kd
 397  kb=(ii-1)*ml +mlr
      lb=(kk-1)*ml +mlr
      ib=ibx+min(kb,lb)+ideks(max(kb,lb))
      idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto) = ib
      if(jto.lt.iwod) go to 384
      jto=0
      jblk=jblk+1
      write(ntype) kc,kd
      go to 384
396   ibx=ideks(mar+1)
      iby=ideks(ibx) +iay
      iby=nit(iby)
      ibx=ibx-mar+iax
      ibx=ideks(ibx+1)
      ibx=nit(ibx)
      ib=(mla-1)*mv +jb+iby
      idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if(jto.lt.iwod) go to 473
      jto=0
      jblk=jblk+1
      write(ntype) kc,kd
 473  kb=(mlr-1)*mq
      lb=kb+ii
      kb=kb+kk
      ib=ibx+min(kb,lb)+ideks(max(kb,lb))
      idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if (jto.lt.iwod) go to 384
      jto=0
      jblk=jblk+1
      write (ntype) kc,kd
 384  continue
      if (kix.lt.0) go to 260
 383  if (nd.eq.0) go to 260
      kix=-1
      jz=nz
      ir=nd
      go to 472
 260  l3=l3+mx
  69  l4=l4+mx
      return
      end
      subroutine finmo(coef,numb)
      implicit REAL (a-h,o-z), integer(i-n)
      character *4 itest,space
      common/work/jrec
      dimension coef(*)
      data space/' '/
      nsd=0
 1    if(nsd.ge.numb)go to 2
      call input
 3    call inpa4(itest)
      if(itest.eq.space)go to 1
      jrec=jrec-1
      call inpf(test)
      nsd=nsd+1
      coef(nsd)=test
      go to 3
 2    if(nsd.ne.numb)call errors(562)
      return
      end
_IF(cray)
      subroutine upackx(ihint)
      implicit REAL (a-h,o-z), integer(i-n)
      common/linkmr/ia,ib,ic
      dimension ihint(*)
      ia=shiftl(ihint(1),24).or.shiftl(ihint(2),16).or.
     *   shiftl(ihint(3),8).or.ihint(4)
      ib=shiftl(ihint(5),24).or.shiftl(ihint(6),16).or.
     *   shiftl(ihint(7),8).or.ihint(8)
      ic=shiftl(ihint(9),24).or.shiftl(ihint(10),16).or.
     *   shiftl(ihint(11),8).or.ihint(12)
      return
      end
_ELSEIF(ksr)
      subroutine upackx(ihint)
      implicit REAL (a-h,o-z), integer(i-n)
      common/linkmr/ia,ib,ic
      dimension ihint(*)
      ia = or (shiftl(ihint(1),24),
     +         or (shiftl(ihint(2),16),
     +               or (shiftl(ihint(3),8),ihint(4))))
      ib = or (shiftl(ihint(5),24),
     +         or (shiftl(ihint(6),16),
     +               or (shiftl(ihint(7),8),ihint(8))))
      ic = or (shiftl(ihint(9),24),
     +         or (shiftl(ihint(10),16),
     +               or (shiftl(ihint(11),8),ihint(12))))
      return
      end
_ELSEIF(i88)
      subroutine upackx(ihint)
      implicit REAL (a-h,o-z), integer(i-n)
      common/linkmr/ia,ib,ic
      dimension ihint(*)
      ia = ior(ishft(ihint(1),24),
     +         ior(ishft(ihint(2),16),
     +              ior(ishft(ihint(3),8),ihint(4))))
      ib = ior(ishft(ihint(5),24),
     +        ior(ishft(ihint(6),16),
     +               ior(ishft(ihint(7),8),ihint(8))))
      ic = ior(ishft(ihint(9),24),
     +         ior(ishft(ihint(10),16),
     +               ior(ishft(ihint(11),8),ihint(12))))
      return
      end
_ELSEIF(i8)
      subroutine upackx(olog)
      implicit REAL (a-h,o-z), integer(i-n)
      logical *1 olog,oia(12)
      integer *4 ic4,i24,i34
      common/linkmr/ic,i2,i3,ic4,i24,i34
      dimension olog(*)
      equivalence (oia(1),ic4)
_IF(littleendian)
      i4=1
_ELSE
      i4=8
_ENDIF
      do 1 loop=1,12
      oia(loop) = olog(i4)
  1   i4=i4+8
      ic = ic4
      i2 = i24
      i3 = i34
      return
      end
_ELSE
      subroutine upackx(olog)
      implicit REAL (a-h,o-z), integer(i-n)
      logical *1 olog,oia
      common/linkmr/oia(12)
      dimension olog(*)
_IF(littleendian)
      i4=1
_ELSE
      i4=4
_ENDIF
      do 1 loop=1,12
      oia(loop) = olog(i4)
  1   i4=i4+4
      return
      end
_ENDIF
      subroutine smrdci(qq,iqq,lword)
      implicit REAL (a-h,o-z), integer(i-n)
      character *8 multxt
      character *8 text,zcomm
INCLUDE(common/sizes)
INCLUDE(common/prints)
INCLUDE(common/iofile)
      common/ftape/
     *mtape,kclab,ntype,nf4,mdisk,ltype,kfile,nf22,
     *ntape,nf32,linf,ideli,
     * nf35(8)
      common /tap/ical,m,nko,mxex,nmul,ispace,nprin
      common/junkc/zcomm(29),text
      common/b/issss(96+3*mxcsf)
      common/craypk/icall,nmull,mm,jspace,nkoo,mxexx,nprinn
      dimension qq(*),iqq(*),multxt(6)
      character*10 charwall
      data multxt/
     *'singlet','doublet','triplet','quartet','quintet','sextet'/
      if(.not.oprint(29)) write(iwr,7777)
 7777 format(1x,104('=')//40x,45('*')/
     1 40x,'Table-ci -- configuration selection module --'/
     2 40x,45('*')/)
_IF1(c)      call szero(issss,96+3*mxcsf)
_IFN1(c)      call izero(96+3*mxcsf,issss,1)
      ical=icall
      nmul=nmull
      if(nmul.le.0.or.nmul.ge.7)call caserr(
     *'invalid spin multiplicity specified')
      text=multxt(nmul)
      m=mm
      ispace=jspace
      nko=nkoo
      mxex=mxexx
       if (mxex.le.0) mxex=2
      nprin=nprinn
c     use conventional 8-byte packing value as maxorb (256)
      maxo = 256
       if (ical.eq.0) call smrd0
       intmax = maxo * (maxo+1)/2
       intrel = lenint(intmax)
       intadd = lenrel(intrel)+1
       lwordd = lword - intrel
       call smrd1(qq(intrel+1),iqq(intadd),iqq(intadd),
     +           iqq(intadd+8),iqq(intadd),iqq(1),intmax,lwordd)
       call rewftn(mtape)
       call rewftn(kclab)
       call rewftn(ntype)
       call rewftn(mdisk)
       call rewftn(ltype)
       call rewftn(kfile)
       call rewftn(ntape)
       call rewftn(linf)
       call rewftn(ideli)
       cpu=cpulft(1)
       if(.not.oprint(29)) write (iwr,3) cpu ,charwall()
    3  format(/2x,' end of label generation at ',
     *f8.2,' seconds',a10,' wall'//)
      return
       end
      subroutine smrd1(qq,iot,kj,mj,mike,ideks,intmax,lword)
      implicit REAL (a-h,o-z), integer(i-n)
      integer olab
_IF(cray,ksr,i8)
       integer olab8
_ENDIF
      character *1 star
INCLUDE(common/sizes)
INCLUDE(common/prints)
      common/craypk/icall(7+mxnshl*mxcsf+mxcsf),
     * mrootx,iipt0,nprtn,msng,mulu,mxset,
     * iipt1,ktest(mxcsf)
INCLUDE(common/iofile)
      common/ftape/
     *mtape,ntape,ntype,nf4,mdisk,ltype,kfile,nf22,
     *nston,nf32,linf,ideli,
     * nf35(8)
      common/tap/ical,ma,nko,mxex,nmul,ispace,nprin
      common /b/ iwod,nid,lg,irsf,i1,i2,i3,i4,nc,lt,lz,nx,nd,i5,i6,i7,
     + i8,i9,j1,iswh,nr,mm,jblk,jto,igmax,nr1,j2,j3,j4,nawk,isec,nshl,
     + nsec,lulu,nrootx,lsng,k7,ipt0,mconf(5),nconf(5),mh,jswh,
     + nn(mxcsf),mn(mxcsf),nsc(mxcsf),
     + imo,mc,mz,ndx,nrx,ntx,ny,nt,iht,ihz,nz,mt,nzx,jy,
     + jx,nd1,nd2,md1,md2,nr3,nr4,mr1,mr2,nt1,nt2,mt1,mt2,ih1,ih2,ih3,
     + ih4,jh2,jh3,jh4,iab,nt3,iac,nr2,m,nad,nad1,ifr1,mx,md,mr
     *,idum1(1)
      common /bx/rmin,smin,crus,crud,crub,sch,jsec,kmax
      common/scrtch/ jdeks(9),kdeks(8),
     2lj(8),olab(2000),nit(666),nytl(5),ndub(5),
     3icon(5),nod(5),jcon(maxorb),nir(maxorb),loc(maxorb)
      common/junk/ican(3),icbn,olab8(250),jdum(8200),mkon(mxcrec),
     *kc(mxcrec),kd(mxcrec),lab(3),ij(8),jtest(mxnshl),nyzl(mxcsf),
     *jkan(mxcrec),nplu(5),jkon(mxnshl*mxcsf)
      dimension qq(*),iot(*),kj(8),mj(8),mike(mxcrec),ltest(mxcsf)
      dimension ideks(intmax)
      equivalence (kc(mxcrec),j9)
c     character *1 dash
c     data dash/'-'/
      data star/'*'/
      call rewftn(linf)
      call rewftn(mtape)
      call rewftn(nston)
      call rewftn(ntype)
      call rewftn(ntape)
      i1=ntape
      i2=ntype
      i3=ical
      i5=nston
      i6=ird
      i7=iwr
      i8=mdisk
      i9=ideli
c     nirm=maxorb
c..   see also mrdci6  
c..   this should fit the value of maxo ... (and have 1 to spare)
c,,      rh/jvl 1999
      nirm=255
      nav=lenwrd()
      igmax=maxig
      if(lword*nav.lt.mxtrm) call caserr(
     +                       'insufficient memory in select')
      im3=nirm*(nirm+1)/2
      nshl=mxnshl
      nsec=mxcsf
      j1=mtape
      j2=ltype
      j3=linf
      j4=kfile
      k7=nprin
      read(mtape) ndum,m,nmul,ispace,mconf,nytl,mxex,jkon,nko,iswh,nyzl
      jto=0
      nopmax=10
      crus=999.0d0
      crub=990.0d0
      crud=940.0d0
      sch=0.0005d0
      irsf=500
_IF1()      lko=80
_IF1()      nb42=99
_IF1()      incr=10
_IF1()      msec=nsec+1
      nlec=40
      nawk=nlec/2+4
_IF1()      nrmax=10
      rmin=0.7072d0
      smin=0.90d0
      ideks(1)=0
      do 1 i=1,im3
   1  ideks(i+1)=ideks(i)+i
      np=nopmax-2
      jdeks(1)=0
      do 3 i=1,np
   3  jdeks(i+1)=jdeks(i)+ideks(i+1)
      np=np-1
      kdeks(1)=0
      do 4 i=1,np
   4  kdeks(i+1)=kdeks(i)+jdeks(i+1)
      read (nston) n,iwod,nid,kdum,imo,kj,mj,lj
      read (nston) nit,lg,ij
      call rewftn(nston)
      ix=0
      do 5 i=1,n
      kap=lj(i)
      if (kap.eq.0) go to 5
      do 6 j=1,kap
      ix=ix+1
      nir(ix)=i
   6  loc(ix)=j
   5  continue
      nx=nmul-2
      ndz=nytl(1)-nx
      do 361 i=1,iswh
      nx=nx+1
      ndz=ndz-1
      nplu(i)=nx
      nod(i)=nplu(i)+i-1
 361  ndub(i)=ndz
      if (ical.gt.1) go to 520
      isec=mxset
      nrootx=mrootx
      ipt0=iipt0
      nprin=nprtn
      lsng=msng
      lulu=mulu
      ipt1=iipt1
      if (mxex.ne.1) then
       lulu=nko
      else
       do 20010 i=1,lulu
       ltest(i)=ktest(i)
       if (ltest(i).le.0.or.ltest(i).gt.nko) call errors(561)
20010  continue
       if (lulu.eq.nko) go to 501
       ix=0
       do 502 i=1,lulu
       jt=ltest(i)
       if (jt.eq.i) go to 503
       jx=(jt-1)*nshl
       do 504 j=1,nshl
       jx=jx+1
       ix=ix+1
 504   jkon(ix)=jkon(jx)
       go to 502
 503   ix=ix+nshl
 502   continue
      endif
 501  ny=0
      if(isec.le.0) isec=mxconf/1000
      i=isec*1000
      if(iabs(ipt0).gt.1)call caserr(
     *'invalid root specification option')
      if(.not.oprint(29)) then
       write(iwr,850)lulu,nrootx,nmul,ispace,i
 850   format(/
     * ' no. of reference configurations to be'/
     * ' used in zero order secular problem   ',i3/
     * ' no. of roots to be treated           ',i3/
     * ' state spin multiplicity              ',i3/
     * ' state spatial symmetry               ',i3/
     * ' maximum no. of configurations to be included',i7/)
       if(lsng.eq.0)write(iwr,9010)
 9010  format(/
     * ' no automatic selection of singly excited configurations')
       if(lsng.ne.0)write(iwr,9011)lsng
 9011  format(/
     * ' include all singly excited configurations with'/
     * ' respect to root function no. ',i3)
       if(nprin.eq.0)write(iwr,9012)
 9012  format(/' default print specified')
       if(nprin.eq.1)write(iwr,9013)
 9013  format(/' print all test species selected')
       if(nprin.eq.2)write(iwr,9014)
 9014  format(/' print all test species')
       write(iwr,9015)
 9015  format(/' *** root certifying input specification : ')
       if(ipt0.eq.0)write(iwr,9016)nrootx
 9016  format(/20x,'lowest ',i2,' roots to be considered')
       if(ipt0.eq.1)write(iwr,9017)
 9017  format(/20x,'roots to be user nominated')
       if(ipt0.eq.-1)write(iwr,9018)
 9018  format(/20x,'trial eigenvectors to be user specified')
       if(ipt1.gt.0)write(iwr,9019)
 9019  format(///
     * ' overlap test between input eigenvectors and'/
     * ' zero order counterparts requested')
       if(mxex.ne.1.or.lulu.ne.nko)go to 9020
       write(iwr,9050)
 9050  format(//
     * ' following reference configurations to be used for'/
     * ' zero order secular problem :')
       write(iwr,9051)(ltest(i),i=1,lulu)
 9051  format(/40x,20i4)
       write(iwr,9021)(star,i=1,129)
 9020  continue
      endif
      ix=0
      do 505 i=1,3
      nc=mconf(i)
      if (nc.eq.0) go to 505
      nt=nytl(i)
      if (ix.lt.lulu) go to 506
      nb=(nc*nt)/iwod+1
      do 507 j=1,nb
 507  read (mtape)
      go to 505
 506  read (mtape) mkon
      ig=0
      do 508 j=1,nc
      do 509 k=1,nt
      ig=ig+1
      jtest(k)=mkon(ig)
      if (ig.lt.iwod) go to 509
      ig=0
      read (mtape) mkon
 509  continue
      kg=-nshl
      do 510 k=1,lulu
      kg=kg+nshl
      if (nyzl(k).ne.nt) go to 510
      l7=kg
      do 511 l=1,nt
      l7=l7+1
      if (jtest(l).ne.jkon(l7)) go to 510
 511  continue
      ix=ix+1
      mn(ix)=k
      nn(ix)=j
      nsc(ix)=i
      if (k.eq.lsng) mh=ix*imo
      if (ix.ne.lulu) go to 508
      go to 512
 510  continue
 508  continue
      go to 505
 512  nb=(nc-j)*nt+ig
      nb=nb/iwod
      if (nb.eq.0) go to 505
      do 513 j=1,nb
 513  read (mtape)
 505  continue
      if (ix.ne.lulu) go to 755
      if(.not.oprint(29)) then
        write(iwr,9021)(star,i=1,129)
        write (iwr,514)
 514    format(/5x,'permuted order of reference configurations')
        write (iwr,515)  (mn(i),i=1,lulu)
 515    format(/20x,20i3)
 9021   format(/1x,129a1)
      endif
      do 516 i=1,3
 516  nconf(i)=0
c *****
      call setsto(iwod,0,iot)
c *****
      do 517 i=1,lulu
      ix=nsc(i)
 517  nconf(ix)=nconf(ix)+1
      jswh=ix
      ix=0
      do 518 i=1,3
      icon(i)=ny
      nx=nytl(i)
      nc=nconf(i)
      if (nc.eq.0) go to 518
      do 519 j=1,nc
      ix=ix+1
      jx=mn(ix)
      jx=(jx-1)*nshl
      do 519 k=1,nx
      jx=jx+1
      ny=ny+1
 519  iot(ny)=jkon(jx)
 518  continue
      if(ny.gt.iwod) go to 762
      write(kfile) ij,nconf,mn,jswh,lulu,nsc,nrootx,isec,nprin,
     x ipt0,nplu,ndub,nod,nir,loc,ideks,mike
      call smrd10(iot,ideks)
      call smrd11(iot,ideks)
      kc(jto+1)=0
      kc(iwod)=0
      write (ntype) kc,kd
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      call rewftn(ntype)
      call rewftn(ntape)
      go to 521
 520  jto=0
      do 884 jblk=1,132000
      read (ntype) kc
      if(j9.eq.0) go to 887
  884 continue
 887  do 885 j=1,iwod
      jto=jto+1
      if(kc(jto).eq.0) go to 886
 885  continue
 886  jto=jto-1
      jblk=jblk-1
      call rewftn(ntype)
c
 521  i10 = 1
      i20 = i10 + igmax
      i30 = i20 + igmax
      i40 = i30 + lenint(igmax)
      last= i40 + lenint(igmax)
      if(last.gt.lword) call caserr(
     + 'insufficient memory for smrd12')
      call smrd12(qq(i10),qq(i20),qq(i30),qq(i40),igmax)
c
      i10 = 1
      i20 = i10 + intmax
      i30 = i20 + intmax
      i40 = i30 + intmax
      last= i40 + mxconf
      if(last.gt.lword) call caserr(
     + 'insufficient memory for smrd13')
      call smrd13(qq(i10),qq(i20),qq(i30),qq(i30),intmax,qq(i40))
      return
755   write(iwr,761) ix,lulu
761   format(/10x,'not all mains on configuration tape '/
     *' no. on tape    ',i3/
     *' no. expected   ',i3//)
      call errors(521)
  762 write(iwr,763) ny,iwod
  763 format(/10x,'storage of mains takes up too much space',2i6)
      call errors(531)
      return
      end
      subroutine smr130(bkay,acoul,aexc,intmax)
      implicit REAL (a-h,o-z), integer(i-n)
      integer olab,ot,oiot,oihog
_IF(cray,ksr,i8)
      integer olab8
_ENDIF
INCLUDE(common/sizes)
INCLUDE(common/prints)
      common/scra /res(1)
      common/scrtch/ndet(5),nsac(5),iaz(5),iez(5),idra(5),
     + idrc(5),jdrc(5),jan(7),jbn(7),ispa2(7),
     + ae(7829),ot(44000),ideks(maxorb),
     + olab(2000),ij(8),nytl(5),nplu(5),ndub(5),nod(5),jkan(mxcrec),
     + nir(maxorb),loc(maxorb),sac(10),
_IF(cray,ksr,i8)
     + olab8(250),f(2304),h(mxcrec),isc(5),bob(3),itest(mxnshl),
_ELSE
     + olab8(250),f(2304),h(mxcrec),isc(5),isd,bob(3),itest(mxnshl),
_ENDIF
     + core,
     + w(mxcsf*(mxcsf+1)/2),
     + ew(mxroot),vect(mxcsf*mxroot),trsum(10,mxroot),
     + wect(mxcsf*mxroot),senk(500),iqr(mxcrec),istm(10)
     +,ht(mxroot),hs(mxcsf),ea(mxroot),iper(10),
     + coff(mxcsf*mxroot),doff(mxcsf),nrec(5),jtrk(mxcsf),
     + ixxx(mxroot+1),mkon(mxcrec),ikan(mxcrec),hy(mxcsf*mxcsf)
     +,ktrk(mxcsf)
      common /b/ iwod,nid,lg,irsf,ntape,ntype,ical,ntab,nc,lt,lz,nx,
     + nd,nston,itape,jtape,mdisk,ideli,mtape,jswh,nr,mm,jblk,jto,
     + igmax,nr1,mtype,linf,kfile,nawk,isec,nshl,
     + nsec,lulu,nrootx,lsng,nprin,ipt0,mconf(5),nconf(5),mh,iswh,
     + id(mxcsf),mn(mxcsf),nsc(mxcsf),
     + imo,mc,mz,ndx,nrx,ntx,ny,nt,iht,ihz,nz,mt,nzx,jy,
     + jx,nd1,nd2,md1,md2,nr3,nr4,mr1,mr2,nt1,nt2,mt1,mt2,ih1,ih2,ih3,
     + ih4,jh2,jh3,jh4,iab,nt3,iac,nr2,m,nad,nad1,ifr1, mx,md,mr
     +,idum1(1)
      common/linkmr/ ic,i2,i3
      common/bufd/ gout(510),nword
      common/blksi3/nsz
      common /bx/ rmin,smin,crus,crud,crub,sch,jsec,kmax
      common/lsort/isppp(16),
     * ii,ia,ie,iz,mzer,n3,kml,ndt,nms,nps,
     * np1,np2,nm1,nm2
      dimension bkay(intmax),acoul(intmax),aexc(intmax)
      character*10 charwall
_IF(cray,ksr,i8)
      dimension af(7885),oihog(48),oiot(10962),j9(56)
      dimension t(100000),iy(8),lout(510),kj(8),lj(8),nj(8),
     *ntil(8),nbal(9),isym(8),jsym(36),lsym(2040),ncomp(256),
     *e(7401),c(400),newya(255)
      equivalence (t(1),af(1)),(af(1),j9(1)),(t(1),ndet(1)),(
     *t(17430),c(1)),(gout(1),lout(1)),
     3(ot(13),oihog(1)),(ot(13),oiot(1))
      equivalence (t(7886),iy(1)),(t(7902),lj(1)),(t(7943),jsym(1)),
     *(t(7894),kj(1)),(t(7910),nj(1)),(t(7918),ntil(1)),
     *(t(7926),nbal(1)),(t(7935),isym(1)),(t(7979),lsym(1)),
     *(t(10019),ncomp(1)),(t(10275),e(1))
_ELSE
      dimension af(7857),oihog(48),oiot(10962),j9(56)
      dimension t(100000),iy(8),lout(510),kj(8),lj(8),nj(8),
     *ntil(8),nbal(9),isym(8),jsym(36),lsym(2040),ncomp(256),
     *e(6633),c(400),newya(255)
      equivalence (t(1),af(1)),(af(1),j9(1)),(t(1),ndet(1)),(
     *t(15558),c(1)),(gout(1),lout(1)),
     3(ot(13),oihog(1)),(ot(13),oiot(1))
      equivalence (t(7858),iy(1)),(t(7866),lj(1)),(t(7887),jsym(1)),
     *(t(7862),kj(1)),(t(7870),nj(1)),(t(7874),ntil(1)),
     *(t(7878),nbal(1)),(t(7883),isym(1)),(t(7905),lsym(1)),
     *(t(8925),ncomp(1)),(t(9053),e(1))
_ENDIF
      mult=nplu(1)+1
      call rewftn(nston)
      read(nston)n,iwod,nid,ksum,imo,kj,iy,lj,nj,nsel,ntil,nbal,isym,jsy
     1m,iorbs,knu,newya,lsym,ncomp,e,c,vnuc,zero
      ibl=0
      mq=ideks(jsec+1)
_IF1()c     do 270 i=1,mq
_IF1()c270  w(i)=0
      call vclr(w,1,mq)
      ih=0
_IFN1(c)      call setsto(44000,0,ot)
      do 6 i=1,iswh
      nc=nconf(i)
      if (nc.gt.0) go to 7
      if (i.lt.3) go to 8
      ibl=ibl+2
      go to 6
8     ibl=ibl+i-1
      go to 6
7     ndt=ndet(i)
      ia=iaz(i)
      kml=nsac(i)
      ie=iez(i)
      iz=isc(i)
      mzer=kml*kml
      nd=ndub(i)
      nps=nplu(i)
      nms=i-1
      niw=nps*nms
      nr=nod(i)
      nx=nytl(i)
      n1=nr+1
      n2=nr+2
      n3=nr-1
      np1=nps-1
      np2=np1+np1
      nd1=nd-1
      nm1=i-2
      nm2=nm1+nm1
      ii=kml+kml+1
      nis=1-ii
      if (i.lt.3) go to 9
      ibl=ibl+1
      j8=i-2
      mc=nconf(j8)
      if(mc.eq.0) go to 10
      ndj=ndet(j8)
      kmj=nsac(j8)
      ja=iaz(j8)
      jz=isc(j8)
      nzer=kml*kmj
      mr=nod(j8)
      js=jan(ibl)
      ks=jbn(ibl)
      ix=nad1
      kss=0
  11  if (kss.eq.ks) go to 5560
      ix=ix+nad
      call rdbak3(js)
      js=js+nsz
      call stopbk3
      call unpack(lout,8,ot(ix),nad)
      kss=kss+1
      go to 11
 5560     kz=iz
      do 12 k=1,nc
      lz=jz
      do 13 l=1,mc
      mm=mm+1
      ll=olab(mm)
      if(mm.lt.nid) go to 14
      mm=0
      read(ntape)olab8
      call unpack(olab8,8,olab,2000)
14    if(ll.eq.0) go to 13
      mm=mm+1
      jj=olab(mm)
      if (mm.lt.nid) go to 15
      mm=0
      read(ntape) olab8
      call unpack(olab8,8,olab,2000)
15    jj=jj*ii+nis
      jto=jto+1
      ip=ot(jj)
      if(ip.eq.255) go to 803
      coul=h(jto)
      if (jto.lt.iwod) go to 16
      jto=0
      read(mtype) h
16    jto=jto+1
      exc=-h(jto)
      if (jto.lt.iwod) go to 805
      jto=0
      read(mtype)h
      go to 805
803   coul=-h(jto)
      if (jto.lt.iwod) go to 904
      jto=0
      read(mtype) h
904   jto=jto+1
      exc=h(jto)
      if(jto.lt.iwod) go to 805
      jto=0
      read(mtype) h
 805  if(ll-2) 800,801,802
 800  bob(1)=-coul-exc
      bob(2)=exc
      bob(3)=coul
      go to 804
 801  bob(1)=-exc
      bob(2)=coul+exc
      bob(3)=-coul
      go to 804
802   bob(1)=exc
      bob(2)=coul
      bob(3)=-coul-exc
804   continue
_IF1()c     do 905 m=1,nzer
_IF1()c905  f(m)=0
      call vclr(f,1,nzer)
      do 18 m=1,kml
      jj=jj+1
      my=ot(jj)
      jj=jj+1
      if (my.eq.0) go to 18
_IF1(iv)      kx=ja+my
_IFN1(iv)      kx=ja+my+ndj
      mike=ot(jj)
      tim=bob(mike)
_IF1(iv)      lx=m-kml
_IF1(iv)      do 22 if=1,kmj
_IF1(iv)      kx=kx+ndj
_IF1(iv)      lx=lx+kml
_IF1(iv) 22   f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv)      call daxpy(kmj,tim,ae(kx),ndj,f(m),kml)
18    continue
      mz=kz
_IF1(c)      call mxma(ae(ie+1),kml,1,f,1,kml,res,1,kml,kml,kml,kmj)
_IFN1(c)      call mxmaa(ae(ie+1),kml,1,f,1,kml,res,1,kml,kml,kml,kmj)
      do 23 m=1,kml
_IF1()c     ly=0
      mz=mz+1
_IF1()c     my=mx
      nz=lz
      ifm=0
      do 231 if=1,kmj
      nz=nz+1
_IF1()c      tim=0
_IF1()c      mx=my
_IF1()c      do 24 ig=1,kml
_IF1()c      ly=ly+1
_IF1()c      mx=mx+1
_IF1()c 24   tim=tim+f(ly)*ae(mx)
      tim=res(ifm+m)
      if (dabs(tim).lt.1.0d-7) go to 231
      w(ideks(mz)+nz)=tim
231   ifm=ifm+kml
23    continue
13    lz=lz+kmj
12    kz=kz+kml
      go to 10
9     if (i.eq.1) go to 25
10    ibl=ibl+1
      j8=i-1
      mc=nconf(j8)
      if (mc.eq.0) go to 25
      js=jan(ibl)
      ks=jbn(ibl)
      ndj=ndet(j8)
      kmj=nsac(j8)
      ja=iaz(j8)
      jz=isc(j8)
      nzer=kml*kmj
      mr=nod(j8)
      mps=nplu(j8)
      mms=j8-1
      ix=nad1
      kss=0
 27   if (kss.eq.ks) go to 5570
      ix=ix+nad
      call rdbak3(js)
      js=js+nsz
      call stopbk3
      call unpack(lout,8,ot(ix),nad)
      kss=kss+1
      go to 27
 5570 call upackx(ot)
      if (oprint(31)) then
       write(6,*) 'smr130 - 5570: ic,i2,i3', ic,i2,i3
      endif
      jc=nr+mult
      j3=jc*kml+1
      kz=iz
      do 28 k=1,nc
      lz=jz
      do 29 l=1,mc
      mm=mm+1
      kk=olab(mm)
      if (mm.lt.nid) go to 30
      mm=0
      read(ntape) olab8
      call unpack(olab8,8,olab,2000)
30    go to (29,31,32,33), kk
 31   mm=mm+1
      ll=olab(mm)
      if (ll.gt.128) ll=ll-256
      if (mm.lt.nid) go to 34
      mm=0
      read(ntape) olab8
      call unpack(olab8,8,olab,2000)
34    mm=mm+1
      jj=olab(mm)
      if (mm.lt.nid) go to 35
      mm=0
      read(ntape) olab8
      call unpack(olab8,8,olab,2000)
35    mm=mm+1
      kk=olab(mm)
      if (mm.lt.nid) go to 36
      mm=0
      read(ntape) olab8
      call unpack(olab8,8,olab,2000)
36    jj=(kk-1)*mr +jj
      jj=jj*ii+ic
      ip=ot(jj)
      if(ip.eq.255) ll=-ll
      jto=jto+1
      if (ll.lt.0) go to 820
      coul = h(jto)
      if (jto.lt.iwod) go to 37
      jto=0
      read(mtype) h
37    jto=jto+1
      exc=-h(jto)
      if (jto.lt.iwod) go to 38
      jto=0
      read(mtype) h
      go to 38
820   coul=-h(jto)
      if (jto.lt.iwod) go to 821
      jto=0
      read(mtype) h
821   jto=jto+1
      ll=-ll
      exc=h(jto)
      if (jto.lt.iwod) go to 38
      jto=0
      read(mtype) h
38    if (ll-2)810,811,812
810   bob(1)=-exc
      bob(2)=-coul
      bob(3)=coul+exc
      go to 49
811   bob(1)=exc
      bob(2)=-coul-exc
      bob(3)=coul
      go to 49
812   bob(1)=coul+exc
      bob(2)=-exc
      bob(3)=-coul
49    continue
_IF1()c     do 39 m=1,nzer
_IF1()c39   f(m)=0
      call vclr(f,1,nzer)
      do 40 m=1,kml
      jj=jj+1
      my=ot(jj)
      jj=jj+1
      if(my.eq.0) go to 40
_IF1(iv)      kx=ja+my
_IFN1(iv)      kx=ja+my+ndj
      mike=ot(jj)
      if(mike.gt.128) go to 42
      tim=bob(mike)
      go to 43
42    tim=-bob(256-mike)
43    continue
_IF1(iv)      lx=m-kml
_IF1(iv)      do 44 if=1,kmj
_IF1(iv)      kx=kx+ndj
_IF1(iv)      lx=lx+kml
_IF1(iv) 44   f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv)      call daxpy(kmj,tim,ae(kx),ndj,f(m),kml)
40    continue
47    mz=kz
_IF1(c)      call mxma(ae(ie+1),kml,1,f,1,kml,res,1,kml,kml,kml,kmj)
_IFN1(c)      call mxmaa(ae(ie+1),kml,1,f,1,kml,res,1,kml,kml,kml,kmj)
      do 45 m=1,kml
_IF1()c     ly=0
      mz=mz+1
      nz=lz
_IF1()c     my=mx
      ifm=0
      do 451 if=1,kmj
      nz=nz+1
_IF1()c      tim=0
_IF1()c      mx=my
_IF1()c      do 46 ig=1,kml
_IF1()c      ly=ly+1
_IF1()c      mx=mx+1
_IF1()c 46   tim=tim+f(ly)*ae(mx)
      tim=res(ifm+m)
      if (dabs(tim).lt.1.0d-7) go to 451
      w(ideks(mz)+nz)=tim
451   ifm=ifm+kml
45    continue
      go to 29
  32  mm=mm+1
      ll=olab(mm)
      if (mm.lt.nid) go to 51
      mm=0
      read(ntape)olab8
      call unpack(olab8,8,olab,2000)
51    jto=jto+1
      sm=h(jto)
      if(jto.lt.iwod) go to 52
      jto=0
      read(mtype) h
52    if(ll.lt.128) go to 53
      sm=-sm
      ll=256-ll
53    jj=ll*kml+i2
_IF1()c     do 55 m=1,nzer
_IF1()c55   f(m)=0
      call vclr(f,1,nzer)
      do 56 m=1,kml
      jj=jj+1
      my=ot(jj)
      if (my.eq.0) go to 56
      if(my.lt.128) go to 58
_IF1(iv)      kx=ja-my+256
_IFN1(iv)      kx=ja-my+256+ndj
      tim=-sm
      go to 59
58    tim=sm
_IF1(iv)      kx=my+ja
_IFN1(iv)      kx=my+ja+ndj
59    continue
_IF1(iv)      lx=m-kml
_IF1(iv)      do 60 if=1,kmj
_IF1(iv)      kx=kx+ndj
_IF1(iv)      lx=lx+kml
_IF1(iv) 60   f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv)      call daxpy(kmj,tim,ae(kx),ndj,f(m),kml)
56    continue
      go to 47
  33  mm=mm+1
      ll=olab(mm)
      if (mm.lt.nid) go to 61
      mm=0
      read(ntape)olab8
      call unpack(olab8,8,olab,2000)
61    mm=mm+1
      ntar=olab(mm)
      if (mm.lt.nid) go to 62
      mm=0
      read (ntape) olab8
      call unpack(olab8,8,olab,2000)
62    mm=mm+1
      nb=olab(mm)
      if(mm.lt.nid) go to 63
      mm=0
      read (ntape) olab8
      call unpack(olab8,8,olab,2000)
63    mm=mm+1
      mb=olab(mm)
      if (mm.lt.nid) go to 41
      mm=0
      read (ntape) olab8
      call unpack(olab8,8,olab,2000)
  41  nb=ideks(nb)+mb+ij(ntar)
      sm=bkay(nb)
      jto=jto+1
      sm=sm+h(jto)
      if(jto.lt.iwod) go to 64
      jto=0
      read(mtype) h
64    if (mr.eq.0 ) go to 65
      do 66 m=1,mr
      jto=jto+1
      sm=sm+h(jto)
      if(jto.lt.iwod) go to 68
      jto=0
      read(mtype) h
68    jto=jto+1
      sac(m)=-h(jto)
      if (jto.lt.iwod) go to 66
      jto=0
      read(mtype) h
66    continue
65    if(nd.eq.0) go to 67
      do 69 m=1,nd
      jto=jto+1
      tim=h(jto)
      sm=sm+tim+tim
      if(jto.lt.iwod) go to 70
      jto=0
      read(mtype) h
70    jto=jto+1
      sm=sm-h(jto)
      if(jto.lt.iwod) go to 69
      jto=0
      read(mtype) h
69    continue
67    if(ll.gt.128) go to 86
      jj=ll*j3+i3
      ip=ot(jj)
      if(ip.eq.1) go to 71
      sm=-sm
      if(mr.eq.0) go to 71
      do 72 m=1,mr
72    sac(m)=-sac(m)
71    jj=jj+1
_IF1()c     do 80 m=1,nzer
_IF1()c80   f(m)=0
      call vclr(f,1,nzer)
_IF1(iv)      jx=-kml
      do 73 m=1,kml
_IF1(iv)      jx=jx+1
      in=jj
      im=ot(in)
      go to (74,75,76,77),im
 74   in=in+1
      kx=ot(in)
_IF1(iv)      kx=kx+ja
_IFN1(iv)      kx=kx+ja+ndj
      tim=sm
      if(mps.eq.0) go to 78
      do 79 if=1,mps
      in=in+1
      im=ot(in)
79    tim=tim+sac(im)
78    continue
_IF1(iv)       lx=jx
_IF1(iv)      do 81 if=1,kmj
_IF1(iv)      kx=kx+ndj
_IF1(iv)      lx=lx+kml
_IF1(iv) 81   f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv)      call daxpy(kmj,tim,ae(kx),ndj,f(m),kml)
      go to 73
 75   in=in+1
      kx=ot(in)
_IF1(iv)      kx=kx+ja
_IFN1(iv)      kx=kx+ja+ndj
      tim=-sm
      if(mms.eq.0) go to 78
      do 83 if=1,mms
      in=in+1
      im=ot(in)
83    tim=tim-sac(im)
      go to 78
76    do 84 if=1,nms
      in=in+1
      kx=ot(in)
_IF1(iv)      kx=kx+ja
_IFN1(iv)      kx=kx+ja+ndj
      in=in+1
      im=ot(in)
      tim=-sac(im)
_IF1(iv)      lx=jx
_IF1(iv)      do 84 ig=1,kmj
_IF1(iv)      kx=kx+ndj
_IF1(iv)      lx=lx+kml
_IF1(iv) 84   f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv) 84   call daxpy(kmj,tim,ae(kx),ndj,f(m),kml)
      go to 73
77    do 85 if=1,nps
      in=in+1
      kx=ot(in)
_IF1(iv)      kx=kx+ja
_IFN1(iv)      kx=kx+ja+ndj
      in=in+1
      im=ot(in)
      tim=sac(im)
_IF1(iv)      lx=jx
_IF1(iv)      do 85 ig=1,kmj
_IF1(iv)      kx=kx+ndj
_IF1(iv)      lx=lx+kml
_IF1(iv) 85   f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv) 85   call daxpy(kmj,tim,ae(kx),ndj,f(m),kml)
73    jj=jj+jc
      go to 47
86    jj=(256-ll)*j3+i3
      ip=ot(jj)
      if (ip.eq.1) go to 87
      sm=-sm
      if (mr.eq.0) go to 87
      do 88 m=1,mr
88    sac(m)=-sac(m)
87    jj=jj+1
_IF1()c     do 89 m=1,nzer
_IF1()c89   f(m)=0
      call vclr(f,1,nzer)
_IF1(iv)      jx=-kml
      do 90 m=1,kml
_IF1(iv)      jx=jx+1
      in=jj
      im=ot(in)
      go to (91,95,96,97),im
91    in=in+1
      kx=ot(in)
_IF1(iv)      kx=kx+ja
_IFN1(iv)      kx=kx+ja+ndj
      tim=sm
      if(mms.eq.0) go to 92
      in=in+mps
      do 93 if=1,mms
      in=in+1
      im=ot(in)
93    tim=tim+sac(im)
92    continue
_IF1(iv)      lx=jx
_IF1(iv)      do 94 if=1,kmj
_IF1(iv)      kx=kx+ndj
_IF1(iv)      lx=lx+kml
_IF1(iv) 94   f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv)      call daxpy(kmj,tim,ae(kx),ndj,f(m),kml)
      go to 90
95    in=in+1
      kx=ot(in)
_IF1(iv)      kx=kx+ja
_IFN1(iv)      kx=kx+ja+ndj
      tim=-sm
      if(mps.eq.0) go to 92
      in=in+mms
      do 99 if=1,mps
      in=in+1
      im=ot(in)
99    tim=tim-sac(im)
      go to 92
96    do 100 if=1,nms
      in=in+1
      kx=ot(in)
_IF1(iv)      kx=kx+ja
_IFN1(iv)      kx=kx+ja+ndj
      in=in+1
      im=ot(in)
      tim=sac(im)
_IF1(iv)      lx=jx
_IF1(iv)      do 100 ig=1,kmj
_IF1(iv)      kx=kx+ndj
_IF1(iv)      lx=lx+kml
_IF1(iv) 100  f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv)100   call daxpy(kmj,tim,ae(kx),ndj,f(m),kml)
      go to 90
97    do 101 if=1,nps
      in=in+1
      kx=ot(in)
_IF1(iv)      kx=kx+ja
_IFN1(iv)      kx=kx+ja+ndj
      in=in+1
      im=ot(in)
      tim=-sac(im)
_IF1(iv)      lx=jx
_IF1(iv)      do 101 ig=1,kmj
_IF1(iv)      kx=kx+ndj
_IF1(iv)      lx=lx+kml
_IF1(iv) 101  f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv) 101  call daxpy(kmj,tim,ae(kx),ndj,f(m),kml)
90    jj=jj+jc
      go to 47
29    lz=lz+kmj
28    kz=kz+kml
25    js=idra(i)
      ks=idrc(i)
      ix=nad1
      kss=0
  102 if (kss.eq.ks) go to 5580
      ix=ix+nad
      call rdbak3(js)
      js=js+nsz
      call stopbk3
      call unpack(lout,8,ot(ix),nad)
      kss=kss+1
      go to 102
 5580 call upackx(ot)
      if (oprint(31)) then
       write(6,*) 'smr130 - 5580: ic,i2,i3', ic,i2,i3
      endif
      kz=iz
      do 103 k=1,nc
      do 104 l=1,nx
      ih=ih+1
      itest(l)=jkan(ih)
104   continue
_IF1()c     do 105 l=1,mzer
_IF1()c105  f(l)=0
      call vclr(f,1,mzer)
      care=core
      if (nr.eq.0) go to 106
      do 107 l=1,nr
      mal=itest(l)
      ntar=nir(mal)
      mal=loc(mal)+1
      mal=ideks(mal)+ij(ntar)
107   care=care+bkay(mal)
      if(nd.eq.0) go to 108
106   do 109 l=n1,nx
      mal=itest(l)
      ntar=nir(mal)
      nal=ideks(mal+1)
      mal=loc(mal)+1
      mal=ideks(mal)+ij(ntar)
      sm=bkay(mal)
109   care=care+sm+sm+acoul(nal)
108   if(nr.lt.2) go to 110
      do 111 l=2,nr
      mal=itest(l)
      mal=ideks(mal)
      it=l-1
      do 111 m=1,it
      nal=mal+itest(m)
111   care=care+acoul(nal)
110   if(nd.lt.2) go to 112
      sm=0
      do 113 l=n2,nx
      mal=itest(l)
      mal=ideks(mal)
      it=l-1
      do 113 m=n1,it
      nal=mal+itest(m)
      tim=acoul(nal)
113   sm=sm+tim+tim-aexc(nal)
      care=care+sm+sm
112   if(nr.eq.0.or.nd.eq.0) go to 114
      do 115 l=1,nr
      mal=itest(l)
      do 115 m=n1,nx
      nal=itest(m)
      nal = min(mal,nal) + ideks( max(mal,nal))
      sm=acoul(nal)
115   care=care+sm+sm-aexc(nal)
114   kk=ic
      kn=i2
      ml=i3
_IF1(iv)      jx=-kml
      do 118 l=1,kml
_IF1(iv)      jx=jx+1
_IF1(iv)      lx=jx
      kx=oihog(l)
_IF1(iv)      kx=kx+ia
_IFN1(iv)      kx=kx+ia+ndt
      tim=care
      if(nps.lt.2) go to 122
      mx=kn+1
      do 123 if=2,nps
      mx=mx+1
      it=if-1
      mal=oiot(mx)
      mal=itest(mal)
      mal=ideks(mal)
      mz=kn
      do 123 in=1,it
      mz=mz+1
      nal=oiot(mz)
      nal=itest(nal)+mal
123   tim=tim-aexc(nal)
      kn=mx
      if(nms.lt.2) go to 122
      mx=ml+1
      do 124 if=2,nms
      mx=mx+1
      mal=oiot(mx)
      mal=itest(mal)
      mal=ideks(mal)
      mz=ml
      it=if-1
      do 124 in=1,it
      mz=mz+1
      nal=oiot(mz)
      nal=itest(nal)+mal
124   tim=tim-aexc(nal)
      ml=mx
122   continue
_IF1(iv)      do 121 if=1,kml
_IF1(iv)      lx=lx+kml
_IF1(iv)      kx=kx+ndt
_IF1(iv) 121  f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv)      call daxpy(kml,tim,ae(kx),ndt,f(l),kml)
      if(niw.eq.0) go to 118
      do 120 m=1,niw
      kk=kk+1
      kx=oiot(kk)
_IF1(iv)      kx=kx+ia
_IFN1(iv)      kx=kx+ia+ndt
      kk=kk+1
      mal=oiot(kk)
      kk=kk+1
      nal=oiot(kk)
      if (k.eq.0) write (6,241) kx,mal,nal,kk
      mal=itest(mal)
      nal=ideks(mal)+itest(nal)
      tim=-aexc(nal)
_IF1(iv)      lx=jx
 241  format(2x,21i6)
_IF1(iv)      do 120 if=1,kml
_IF1(iv)      lx=lx+kml
_IF1(iv)      kx=kx+ndt
_IF1(iv) 120  f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv) 120  call daxpy(kml,tim,ae(kx),ndt,f(l),kml)
_IF1(iv)      if (k.eq.0) write (6,239) lx,kx,mal,nal,tim,f(lx),ae(kx)
_IF1(iv)239   format(2x,4i6,3f20.8)
118   continue
      mz=kz
      mx=ie
      do 125 l=1,kml
      ly=1
      nz=kz
      mz=mz+1
      do 1251 m=1,l
      nz=nz+1
      tim=ddot(kml,f(ly),1,ae(mx+1),1)
      if(l.gt.m) go to 127
      jg=ideks(mz)+nz
      w(jg)=tim
      go to 1251
127   if(dabs(tim).lt.1.0d-7) go to 1251
      jg=ideks(mz)+nz
      w(jg)=tim
1251  ly=ly+kml
125   mx=mx+kml
103   kz=mz
      if(nc.eq.1) go to 6
      ks=jdrc(i)
      ix=nad1
      kss=0
 128  if (kss.eq.ks) go to 5590
      ix=ix+nad
      call rdbak3(js)
      js=js+nsz
      call stopbk3
      call unpack(lout,8,ot(ix),nad)
      kss=kss+1
      go to 128
 5590 call upackx(ot)
      if (oprint(31)) then
       write(6,*) 'smr130 - 5590: ic,i2,i3', ic,i2,i3
      endif
      call sm1300(bkay,aexc,intmax)
6     continue
      call clredx
      call closbf3
      cpu=cpulft(1)
      if(.not.oprint(29)) write (jtape,8000) cpu ,charwall()
 8000 format(//' end of zero order ci at ',f8.2,' seconds',a10,' wall')
      return
      end
      subroutine sm1300(bkay,aexc,intmax)
      implicit REAL (a-h,o-z), integer(i-n)
      integer olab,ot,oiot,oihog
_IF(cray,ksr,i8)
      integer olab8
_ENDIF
INCLUDE(common/sizes)
      common/scra /res(1)
      common/scrtch/ndet(5),nsac(5),iaz(5),iez(5),idra(5),
     + idrc(5),jdrc(5),jan(7),jbn(7),ispa2(7),
     + ae(7829),ot(44000),ideks(maxorb),
     + olab(2000),ij(8),nytl(5),nplu(5),ndub(5),nod(5),jkan(mxcrec),
     + nir(maxorb),loc(maxorb),sac(10),
_IF(cray,ksr,i8)
     + olab8(250),f(2304),h(mxcrec),isc(5),bob(3),itest(mxnshl),
_ELSE
     + olab8(250),f(2304),h(mxcrec),isc(5),isd,bob(3),itest(mxnshl),
_ENDIF
     + core,
     + w(mxcsf*(mxcsf+1)/2),
     + ew(mxroot),vect(mxcsf*mxroot),trsum(10,mxroot),
     + wect(mxcsf*mxroot),senk(500),iqr(mxcrec),istm(10)
     +,ht(mxroot),hs(mxcsf),ea(mxroot),iper(10),
     + coff(mxcsf*mxroot),doff(mxcsf),nrec(5),jtrk(mxcsf),
     + ixxx(mxroot+1),mkon(mxcrec),ikan(mxcrec),hy(mxcsf*mxcsf)
     +,ktrk(mxcsf)
      common /b/ iwod,nid,lg,irsf,ntape,ntype,ical,ntab,nc,lt,lz,nx,
     + nd,nston,itape,jtape,mdisk,ideli,mtape,jswh,nr,mm,jblk,jto,
     + igmax,nr1,mtype,linf,kfile,
     + nawk,isec,nshl,
     + nsec,lulu,nrootx,lsng,nprin,ipt0,mconf(5),nconf(5),mh,iswh,
     + id(mxcsf),mn(mxcsf),nsc(mxcsf),
     + imo,mc,mz,ndx,nrx,ntx,ny,nt,iht,ihz,nz,mt,nzx,jy,
     + jx,nd1,nd2,md1,md2,nr3,nr4,mr1,mr2,nt1,nt2,mt1,mt2,ih1,ih2,ih3,
     + ih4,jh2,jh3,jh4,iab,nt3,iac,nr2,m,nad,nad1,ifr1, mx,md,mr
     +,idum1(1)
      common/linkmr/ ic,i2,i3
      common/bufd/ gout(510),nword
      common/blksi3/nsz
      common /bx/ rmin,smin,crus,crud,crub,sch,jsec,kmax
      common/lsort/isppp(16),
     * ii,ia,ie,iz,mzer,n3,kml,ndt,nms,nps,
     * np1,np2,nm1,nm2
      dimension bkay(intmax),aexc(intmax)
_IF(cray,ksr,i8)
      dimension af(7885),oihog(48),oiot(10962),j9(56)
      dimension t(100000),iy(8),lout(510),kj(8),lj(8),nj(8),
     *ntil(8),nbal(9),isym(8),jsym(36),lsym(2040),ncomp(256),
     *e(7401),c(400)
      equivalence (t(1),af(1)),(af(1),j9(1)),(t(1),ndet(1)),(
     *t(17430),c(1)),(gout(1),lout(1)),
     3(ot(13),oihog(1)),(ot(13),oiot(1))
      equivalence (t(7886),iy(1)),(t(7902),lj(1)),(t(7943),jsym(1)),
     *(t(7894),kj(1)),(t(7910),nj(1)),(t(7918),ntil(1)),
     *(t(7926),nbal(1)),(t(7935),isym(1)),(t(7979),lsym(1)),
     *(t(10019),ncomp(1)),(t(10275),e(1))
_ELSE
      dimension af(7857),oihog(48),oiot(10962),j9(56)
      dimension t(100000),iy(8),lout(510),kj(8),lj(8),nj(8),
     *ntil(8),nbal(9),isym(8),jsym(36),lsym(2040),ncomp(256),
     *e(6633),c(400)
      equivalence (t(1),af(1)),(af(1),j9(1)),(t(1),ndet(1)),(
     *t(15558),c(1)),(gout(1),lout(1)),
     3(ot(13),oihog(1)),(ot(13),oiot(1))
      equivalence (t(7858),iy(1)),(t(7866),lj(1)),(t(7887),jsym(1)),
     *(t(7862),kj(1)),(t(7870),nj(1)),(t(7874),ntil(1)),
     *(t(7878),nbal(1)),(t(7883),isym(1)),(t(7905),lsym(1)),
     *(t(8925),ncomp(1)),(t(9053),e(1))
_ENDIF
      ii=ii+kml
      j2=kml+1
      jc=nr+nr
      j3=jc*kml+1
      kz=iz+kml
      j8=0
      do 129 j=2,nc
      lz=iz
      j8=j8+1
      do 130 k=1,j8
      mm=mm+1
      kk=olab(mm)
      if (mm.lt.nid) go to 131
      mm=0
      read(ntape)olab8
      call unpack(olab8,8,olab,2000)
 131  go to (130,132,133,134,135),kk
 135  mm=mm+1
      kk=olab(mm)
      if(mm.lt.nid) go to 136
      mm=0
      read(ntape)olab8
      call unpack(olab8,8,olab,2000)
136   mm=mm+1
      ll=olab(mm)
      if(mm.lt.nid) go to 137
      mm=0
      read(ntape)olab8
      call unpack(olab8,8,olab,2000)
137   kk=min(kk,ll)+ideks(max(kk,ll))
_IF1()c     do 140 l=1,mzer
_IF1()c140  f(l)=0
      call vclr(f,1,mzer)
      tim=aexc(kk)
_IF1(iv)      jx=-kml
      do 141 l=1,kml
      kx=oihog(l)
_IF1(iv)      kx=kx+ia
_IFN1(iv)      kx=kx+ia+ndt
_IF1(iv)      jx=jx+1
_IF1(iv)       lx=jx
_IF1(iv)      do 141 m=1,kml
_IF1(iv)      lx=lx+kml
_IF1(iv)      kx=kx+ndt
_IF1(iv) 141  f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv) 141  call daxpy(kml,tim,ae(kx),ndt,f(l),kml)
146   mz=kz
_IF1(c)      call mxma(ae(ie+1),kml,1,f,1,kml,res,1,kml,kml,kml,kml)
_IFN1(c)      call mxmaa(ae(ie+1),kml,1,f,1,kml,res,1,kml,kml,kml,kml)
      do 145 l=1,kml
_IF1()c     ly=0
      mz=mz+1
_IF1()c     my=mx
      nz=lz
      ifm=0
      do 1451 m=1,kml
      nz=nz+1
_IF1()c      tim=0
_IF1()c      mx=my
_IF1()c      do 143 if=1,kml
_IF1()c      ly=ly+1
_IF1()c      mx=mx+1
_IF1()c 143  tim=tim+f(ly)*ae(mx)
       tim=res(ifm+l)
      if(dabs(tim).lt.1.0d-7) go to 1451
      jg=ideks(mz)+nz
      w(jg)=tim
1451  ifm=ifm+kml
145   continue
      go to 130
_IF1(c) 144  call mxma(ae(ie+1),kml,1,f,1,kml,res,1,kml,kml,kml,kml)
_IFN1(c) 144  call mxmaa(ae(ie+1),kml,1,f,1,kml,res,1,kml,kml,kml,kml)
      nzz=lz
      do 253 l=1,kml
_IF1()c     ly=0
      nzz=nzz+1
_IF1()c     my=mx
      mzz=kz
      ifm=0
      do 2531 m=1,kml
      mzz=mzz+1
_IF1()c      tim=0
_IF1()c      mx=my
_IF1()c      do 254 if=1,kml
_IF1()c      ly=ly+1
_IF1()c      mx=mx+1
_IF1()c 254  tim=tim+f(ly)*ae(mx)
       tim=res(ifm+l)
      if (dabs(tim).lt.1.0d-7) go to 2531
      jg=ideks(mzz)+nzz
      w(jg)=tim
2531  ifm=ifm+kml
 253  continue
      go to 130
 132  mm=mm+1
      ll=olab(mm)
      if(mm.lt.nid) go to 147
      mm=0
      read(ntape)olab8
      call unpack(olab8,8,olab,2000)
147   mm=mm+1
       jj=olab(mm)
      if(mm.lt.nid) go to 148
      mm=0
      read(ntape)olab8
      call unpack(olab8,8,olab,2000)
148   mm=mm+1
      kk=olab(mm)
      if(mm.lt.nid) go to 149
      mm=0
      read(ntape)olab8
      call unpack(olab8,8,olab,2000)
149   if(jj.lt.kk) go to 150
      jj=ideks(jj)+kk
      ibob=0
      go to 151
150   jj=ideks(kk)+jj
      ibob=1
151   jj=jj*ii+ic
      jto=jto+1
      ip=ot(jj)
      if(ip.eq.255) go to 830
      coul=h(jto)
      if(jto.lt.iwod) go to 152
      jto=0
      read(mtype) h
152   jto=jto+1
      exc=-h(jto)
      if(jto.lt.iwod) go to 153
      jto=0
      read(mtype) h
      go to 153
830   coul=-h(jto)
      if (jto.lt.iwod) go to 831
      jto=0
      read (mtype) h
831   jto=jto+1
      exc=h(jto)
      if(jto.lt.iwod) go to 153
      jto=0
      read (mtype) h
  153 if (ll-2) 832,833,834
  832 bob(1)=coul+exc
      bob(2)=exc
      bob(3)=coul
      go to 835
  833 bob(1)=exc
      bob(2)=-coul
      bob(3)=coul+exc
      go to 835
 834  bob(1)=-exc
      bob(2)=-coul-exc
      bob(3)=coul
835   continue
_IF1()c     do 154 l=1,mzer
_IF1()c154  f(l)=0
      call vclr(f,1,mzer)
_IF1(iv)      mx=-kml
      do 155 l=1,kml
      jj=jj+1
_IF1(iv)      mx=mx+1
      ip=ot(jj)
      if(ip.eq.2)  go to 159
      jj=jj+1
      kx=ot(jj)
_IF1(iv)      kx=kx+ia
_IFN1(iv)      kx=kx+ia+ndt
      jj=jj+1
      tim=bob(1)
_IF1(iv)      lx=mx
_IF1(iv)      do 161 m=1,kml
_IF1(iv)      kx=kx+ndt
_IF1(iv)      lx=lx+kml
_IF1(iv) 161  f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv)      call daxpy(kml,tim,ae(kx),ndt,f(l),kml)
      go to 155
159   jj=jj+1
      kx=ot(jj)
_IF1(iv)      kx=kx+ia
_IFN1(iv)      kx=kx+ia+ndt
      tim=bob(2)
_IF1(iv)      lx=mx
_IF1(iv)      do 160 m=1,kml
_IF1(iv)      kx=kx+ndt
_IF1(iv)      lx=lx+kml
_IF1(iv) 160  f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv)      call daxpy(kml,tim,ae(kx),ndt,f(l),kml)
      jj=jj+1
      kx=ot(jj)
_IF1(iv)      kx=kx+ia
_IFN1(iv)      kx=kx+ia+ndt
      tim=bob(3)
_IF1(iv)      lx=mx
_IF1(iv)      do 167 m=1,kml
_IF1(iv)      lx=lx+kml
_IF1(iv)      kx=kx+ndt
_IF1(iv) 167  f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv)      call daxpy(kml,tim,ae(kx),ndt,f(l),kml)
155   continue
      if (ibob.eq.1) go to 144
      go to 146
 133  mm=mm+1
      ll=olab(mm)
      if (mm.lt.nid) go to 168
      mm=0
      read(ntape)olab8
      call unpack(olab8,8,olab,2000)
168   mm=mm+1
      jj=olab(mm)
      if(mm.lt.nid) go to 169
      mm=0
      read(ntape)olab8
      call unpack(olab8,8,olab,2000)
169   if(jj.gt.ll) go to 170
      jj=ideks(ll)+jj
      ibob=0
      go to 171
170   jj=ideks(jj)+ll
      ibob=1
171   jto=jto+1
      tim=h(jto)
      if(jto.lt.iwod) go to 172
      jto=0
      read(mtype)h
172   jj=jj*j2+i2
_IF1()c     do 177 l=1,mzer
_IF1()c177  f(l)=0
      call vclr(f,1,mzer)
_IF1(iv)      jx=-kml
      ip=ot(jj)
      if(ip.eq.255) tim=-tim
      do 173 l=1,kml
_IF1(iv)      jx=jx+1
      jj=jj+1
      kx=ot(jj)
_IF1(iv)      kx=kx+ia
_IFN1(iv)      kx=kx+ia+ndt
_IF1(iv)      lx=jx
_IF1(iv)      do 173 m=1,kml
_IF1(iv)      kx=kx+ndt
_IF1(iv)      lx=lx+kml
_IF1(iv) 173  f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv) 173  call daxpy(kml,tim,ae(kx),ndt,f(l),kml)
      if (ibob.eq.1) go to 144
      go to 146
 134  mm=mm+1
      ll=olab(mm)
      if(mm.lt.nid) go to 179
      mm=0
      read(ntape)olab8
      call unpack(olab8,8,olab,2000)
179   mm=mm+1
      jj=olab(mm)
      if(mm.lt.nid) go to 180
      mm=0
      read(ntape)olab8
      call unpack(olab8,8,olab,2000)
180   mm=mm+1
      ntar=olab(mm)
      if(mm.lt.nid) go to 183
      mm=0
      read(ntape)olab8
      call unpack(olab8,8,olab,2000)
183   mm=mm+1
      nb=olab(mm)
      if(mm.lt.nid) go to 184
      mm=0
      read(ntape)olab8
      call unpack(olab8,8,olab,2000)
184   mm=mm+1
      mb=olab(mm)
      if (mm.lt.nid) go to 48
      mm=0
      read (ntape) olab8
      call unpack(olab8,8,olab,2000)
  48  nb=ideks(nb)+mb+ij(ntar)
      sm=bkay(nb)
      if(ll.gt.128) go to 181
      if(ll.lt.jj) go to 182
      ibob=0
      jj=ideks(ll)+jj
      go to 205
182   jj=ideks(jj)+ll
      ibob=1
205   jj=jj*j3+i3
      if(nr.eq.1) go to 185
      do 186 l=1,n3
      jto=jto+1
      sm=sm+h(jto)
      if(jto.lt.iwod) go to 187
      jto=0
      read(mtype) h
187   jto=jto+1
      sac(l)=-h(jto)
      if(jto.lt.iwod) go to 186
      jto=0
      read(mtype) h
186   continue
185   if(nd.eq.0) go to 188
      do 189 l=1,nd
      jto=jto+1
      tim=h(jto)
      sm=sm+tim+tim
      if(jto.lt.iwod) go to 190
      jto=0
      read(mtype) h
190   jto=jto+1
      sm=sm-h(jto)
      if(jto.lt.iwod) go to 189
      jto=0
      read(mtype) h
189   continue
188   ip=ot(jj)
      if(ip.eq.1) go to 191
      sm=-sm
      if (nr.eq.1) go to 191
      do 192 l=1,n3
192   sac(l)=-sac(l)
191   jj=jj+1
_IF1()c     do 193 l=1,mzer
_IF1()c193  f(l)=0
      call vclr(f,1,mzer)
_IF1(iv)      jx=-kml
      do 194 l=1,kml
      in=jj
      im=ot(in)
_IF1(iv)      jx=jx+1
      in=in+1
      kx=ot(in)
_IF1(iv)      kx=kx+ia
_IFN1(iv)      kx=kx+ia+ndt
_IF1(iv)      lx=jx
      tim=sm
      if(im.eq.2) go to 195
      if(nps.eq.1) go to 196
      do 197 m=1,np1
      in=in+2
      nn=ot(in)
197   tim=tim+sac(nn)
196   continue
_IF1(iv)      do 198 m=1,kml
_IF1(iv)      lx=lx+kml
_IF1(iv)      kx=kx+ndt
_IF1(iv) 198  f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv)      call daxpy(kml,tim,ae(kx),ndt,f(l),kml)
      if(nms.eq.0) go to 194
      do 199 m=1,nms
      in=in+1
      kx=ot(in)
_IF1(iv)      kx=kx+ia
_IFN1(iv)      kx=kx+ia+ndt
_IF1(iv)      lx=jx
      in=in+1
      mx=ot(in)
      tim=sac(mx)
_IF1(iv)      do 199 if=1,kml
_IF1(iv)      lx=lx+kml
_IF1(iv)      kx=kx+ndt
_IF1(iv) 199  f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv)199   call daxpy(kml,tim,ae(kx),ndt,f(l),kml)
      go to 194
195   if(nms.eq.1) go to 200
      do 201 m=1,nm1
      in=in+2
      nn=ot(in)
201   tim=tim+sac(nn)
200   continue
_IF1(iv)      do 202 m=1,kml
_IF1(iv)      lx=lx+kml
_IF1(iv)      kx=kx+ndt
_IF1(iv) 202  f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv)      call daxpy(kml,tim,ae(kx),ndt,f(l),kml)
      do 203 m=1,nps
      in=in+1
      kx=ot(in)
_IF1(iv)      kx=kx+ia
_IFN1(iv)      kx=kx+ia+ndt
_IF1(iv)      lx=jx
      in=in+1
      mx=ot(in)
      tim=sac(mx)
_IF1(iv)      do 203 if=1,kml
_IF1(iv)      lx=lx+kml
_IF1(iv)      kx=kx+ndt
_IF1(iv) 203  f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv)203   call daxpy(kml,tim,ae(kx),ndt,f(l),kml)
194   jj=jj+jc
      if (ibob.eq.1) go to 144
      go to 146
181   ll=256-ll
      if(ll.lt.jj) go to 204
      ibob=0
      jj=ideks(ll)+jj
      go to 206
204   jj=ideks(jj)+ll
      ibob=1
206   jj=jj*j3+i3
      jto=jto+1
      sm=sm+h(jto)
      if(jto.lt.iwod) go to 207
      jto=0
      read(mtype) h
207   jto=jto+1
      sm=sm+h(jto)
      if(jto.lt.iwod) go to 208
      jto=0
      read(mtype) h
208   if(nr.eq.1) go to 209
      do 211 l=1,n3
      jto=jto+1
      sm=sm+h(jto)
      if(jto.lt.iwod) go to 210
      jto=0
      read(mtype) h
210   jto=jto+1
      sac(l)=-h(jto)
      if(jto.lt.iwod) go to 211
      jto=0
      read(mtype) h
211   continue
209   if(nd.eq.1) go to 212
      do 213 l=1,nd1
      jto=jto+1
      tim=h(jto)
      sm=sm+tim+tim
      if(jto.lt.iwod) go to 214
      jto=0
      read(mtype) h
214   jto=jto+1
      sm=sm-h(jto)
      if (jto.lt.iwod) go to 213
      jto=0
      read(mtype) h
213   continue
212   ip=ot(jj)
      if(ip.eq.1) go to 215
      sm=-sm
      if(nr.eq.1) go to 215
      do 216 l=1,n3
216   sac(l)=-sac(l)
215   jj=jj+1
_IF1()c     do 217 l=1,mzer
_IF1()c217  f(l)=0
      call vclr(f,1,mzer)
_IF1(iv)      jx=-kml
      do 218 l=1,kml
      in=jj
      im=ot(in)
_IF1(iv)      jx=jx+1
      in=in+1
      kx=ot(in)
_IF1(iv)      kx=kx+ia
_IFN1(iv)      kx=kx+ia+ndt
_IF1(iv)      lx=jx
      tim=-sm
      if (im.eq.2) go to 219
      if(nms.eq.0) go to 220
      in=in+np2
      jn=in
      do 221 m=1,nms
      jn=jn+2
      nn=ot(jn)
221   tim=tim-sac(nn)
220   continue
_IF1(iv)      do 222 m=1,kml
_IF1(iv)      lx=lx+kml
_IF1(iv)      kx=kx+ndt
_IF1(iv) 222  f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv)      call daxpy(kml,tim,ae(kx),ndt,f(l),kml)
      if(nms.eq.0) go to 218
      do 223 m=1,nms
      in=in+1
      kx=ot(in)
_IF1(iv)      kx=kx+ia
_IFN1(iv)      kx=kx+ia+ndt
_IF1(iv)      lx=jx
      in=in+1
      mx=ot(in)
      tim=sac(mx)
_IF1(iv)      do 223 if =1,kml
_IF1(iv)      lx=lx+kml
_IF1(iv)      kx=kx+ndt
_IF1(iv) 223  f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv)223   call daxpy(kml,tim,ae(kx),ndt,f(l),kml)
      go to 218
219   in=in+nm2
      jn=in
      do 224 m=1,nps
      jn=jn+2
      nn=ot(jn)
224   tim=tim-sac(nn)
_IF1(iv)      do 225 m=1,kml
_IF1(iv)      lx=lx+kml
_IF1(iv)      kx=kx+ndt
_IF1(iv) 225  f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv)      call daxpy(kml,tim,ae(kx),ndt,f(l),kml)
      do 226 m=1,nps
      in=in+1
      kx=ot(in)
_IF1(iv)      kx=kx+ia
_IFN1(iv)      kx=kx+ia+ndt
_IF1(iv)      lx=jx
      in=in+1
      mx=ot(in)
      tim=sac(mx)
_IF1(iv)      do 226 if=1,kml
_IF1(iv)      lx=lx+kml
_IF1(iv)      kx=kx+ndt
_IF1(iv) 226  f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv)226   call daxpy(kml,tim,ae(kx),ndt,f(l),kml)
218   jj=jj+jc
      if (ibob.eq.1) go to 144
      go to 146
130   lz=lz+kml
129   kz=kz+kml
      return
      end
      subroutine smrd13(bkay,acoul,aexc,jdeks,intmax,sx)
      implicit REAL (a-h,o-z), integer(i-n)
       integer olab,ot,oiot,oihog
_IF(cray,ksr,i8)
       integer olab8
_ENDIF
INCLUDE(common/sizes)
INCLUDE(common/prints)
      common/linkmr/ ic,i2,i3
      common/bufd/ gout(510),nword
      common/blksi3/nsz
      common/scrtch/ndet(5),nsac(5),iaz(5),iez(5),idra(5),
     + idrc(5),jdrc(5),jan(7),jbn(7),ispa2(7),ae(7829),
     + ot(44000),ideks(maxorb),
     + olab(2000),ij(8),nytl(5),nplu(5),ndub(5),nod(5),jkan(mxcrec),
     + nir(maxorb),loc(maxorb),sac(10),
_IF(cray,ksr,i8)
     + olab8(250),f(2304),h(mxcrec),isc(5),bob(3),itest(mxnshl),
_ELSE
     + olab8(250),f(2304),h(mxcrec),isc(5),isd,bob(3),itest(mxnshl),
_ENDIF
     + core,
     + w(mxcsf*(mxcsf+1)/2),
     + ew(mxroot),vect(mxcsf*mxroot),trsum(10,mxroot),
     + wect(mxcsf*mxroot),senk(500),iqr(mxcrec),istm(10)
     + ,ht(mxroot),hs(mxcsf),ea(mxroot),iper(10),
     + coff(mxcsf*mxroot),doff(mxcsf),nrec(5),jtrk(mxcsf),
     + isym(mxroot+1),mkon(mxcrec),ikan(mxcrec),hy(mxcsf*mxcsf)
     + ,ktrk(mxcsf)
      common/scra/vred(mxcsf,mxcsf),temp(mxcsf),ilifq(mxcsf)
      common/craypk/icall(14+mxcsf*mxnshl+mxcsf+mxcsf),isymm(mxroot),
     *cofff(mxcsf*mxroot),t0,t1
      common/jany/ jerk(10),jbnk(10),jbun(10)
      common /b/ iwod,nid,lg,irsf,ntape,ntype,ical,ntab,nc,lt,lz,nx,
     + nd,nston,itape,jtape,mdisk,ideli,mtape,iswh,nr,mm,jblk,jto,
     + igmax,nr1,mtype,linf,kfile,
     + nawk,isec,nshl,
     + nsec,lulu,nrootx,lsng,nprin,ipt0,nconf(5),mconf(5),mh,jswh,
     + id(mxcsf),mn(mxcsf),nsc(mxcsf),
     + imo,mc,mz,ndx,nrx,ntx,ny,nt,iht,ihz,nz,mt,nzx,jy,
     + jx,nd1,nd2,md1,md2,nr3,nr4,mr1,mr2,nt1,nt2,mt1,mt2,ih1,ih2,ih3,
     + ih4,jh2,jh3,jh4,iab,nt3,iac,nr2,m,nad,nad1,ifr1, mx,md,mr
     +,ipag
       common /d/ jbab,kml,ifl,ii,nis,iiz,nzer,kmj,ndj,jsac,ie,iej,
     * ja,kc,kmk,ndk,kr,iiq,ic1,ixex,i21,i2ex,i31,i3ex,iix,
     * kzer,ndt,iek,j3,j3a,mps,nps,mms,nms,jc,jca,lcl,ndi,kmi,
     * ia,jai,nzei,nisi,lc,kk,ndl,km2,la,lzer,lr,lps,lms,
     * ic2,i22,i32,jcb,j3b
      common /bx/rmin,smin,crus,crud,crub,sch,jsec,kmax
      common/lsort/isppp(16),
     * trwa,trwb,trwc,trwd,trwe,trwf,trwg,trwh,trwi,trwj,
     * dmin,dmax,jg,ik,mzer,jpag,ica,
     * n1,n2,n3,nm1,nm2,np1,np2,niw,
     * kft,lpix,ifr,ihp,i3a,i2a,j2,i
      dimension bkay(intmax),acoul(intmax),aexc(intmax),jdeks(intmax)
      dimension sx(*),ifront(10)
      dimension a(5292),e(2304),lout(510),
     1hp(4800),xemp(500),jkon(mxnshl*mxcsf),itym(mxroot),
     2ihog(48),imap(504),j9(56),oihog(48),
     3oiot(10962),t(50000)
_IF(cray,ksr,i8)
      dimension af(7885),ibm(10454)
_ELSE
      dimension af(7857),ribm( 8360)
_ENDIF
      equivalence (t(1),j9(1),ndet(1)),(gout(1),lout(1)),
     6 (t(1),af(1)),(ew(1),itym(1)),
_IF(cray,ksr,i8)
     7 (t(7886),jkon(1)),(ibm(1),t(7886)),
_ELSE
     7 (t(7858),jkon(1)),(ribm(1),t(7858)),
_ENDIF
     3 (ot(13),oihog(1)),(oihog(1),oiot(1)),
     4 (t(1),xemp(1)),(a(1),hy(1)),(a(1),e(1)),(hp(1),a(1)),
_IF(cray,ksr,i8)
     5 (a(2305),imap(1)),(a(2809),ihog(1))
_ELSE
     5 (a(2305),imap(1)),(a(2557),ihog(1))
_ENDIF
      data ifront / 1, 2, 6, 20, 70, 1, 1, 2, 5, 14 /
      call rewftn(ntape)
      call rewftn(ideli)
      nrmx = mxcsf
      jdeli=19
      jfile=18
      call rewftn(linf)
      call rewftn(kfile)
      if(ical.gt.1) call rewftn(jdeli)
      if(ical.lt.3) goto 156
      call rewftn(jfile)
      read(jfile) ij,mconf,mn,jswh,lulu,nsc,nrootx,isec,nprin,
     1ipto,nplu,ndub,nod,nir,loc,jdeks,jkan
      write(kfile) ij,mconf,mn,jswh,lulu,nsc,nrootx,isec,nprin,
     1ipto,nplu,ndub,nod,nir,loc,jdeks,jkan
  156 nrmap=10
         imsec=mxconf
         nad=4000
         nad1=-3999
      call rewftn(kfile)
      call rewftn(mtape)
      read (kfile) ij,mconf,mn,jswh,lulu,nsc,nrootx,isec,nprin,
     * ipt0,nplu,ndub,nod,nir,loc,jdeks,jkan
      if(ical.gt.1) ipt0=0
      do 602 i=1,maxorb
 602  ideks(i)=jdeks(i)
      read (mtape) imo,m,nmul,ispace,nconf,nytl,mxex,jkon,nko,iswh
      mult=nplu(1)+1
_IFN(parallel)
      call setbfc
_ELSE
c **** MPP
      call closbf3
      call setbfc(-1)
c **** MPP
_ENDIF
c
c ... now check the contents of the table data base
c ... must exhibit controlled abort if non-valid first record (usually 
c ... due to failure to assign correctly) rather than floating point error!
      call stopbk3
      is = 0
      call rdbak3(is)
      call stopbk3
      do 6612 i = 1,10
      if (lout(i).ne.ifront(i)) then
       write (jtape,6613) i, lout(i), ifront(i)
6613   format(/1x,'*** table data base in error: word ',i3,' = ',i4,
     + ' and should be ',i4/
     +         1x,'*** check assignment of table data set'//)
       call caserr('table data base has incorrect format')
      endif
6612  continue
      if(.not.oprint(29)) write(jtape,6614)
6614  format(/1x,'*** table data base assigned correctly')
c
_IFN1(c)      call setsto(2000,0,olab)
      is=jerk(mult)
      il=jbnk(mult)
      jz=jbun(mult)
      ifr1=1-irsf
      ix=ifr1
      ill=2
      call rdbak3(is)
      call stopbk3
_IF1(c)      call fmove(lout,ndet,49)
_IFN1(c)      call icopy(49,lout,1,ndet,1)
      is=is+nsz
      ik=1
 5551 if (il.lt.ill) go to 5550
      call rdbak3(is)
      is=is+nsz
      call stopbk3
      call dcopy(nword,gout(1),1,ae(ik),1)
      ik=ik+nword
      ill=ill+1
      go to 5551
 5550     jsec=0
      do 634 i=1,jswh
      isc(i)=jsec
 634  jsec=jsec+mconf(i)*nsac(i)
      jsac=jsec
       if(ipt0) 680, 600,681
  680  do 780 i=1,lulu
      do 781 j=1,lulu
      if(mn(j).eq.i) goto 780
 781   continue
780    jtrk(i)=j
      ix=0
       do 782 i=1,lulu
       ktrk(i)=ix
      jj=jtrk(i)
       ns=nsc(jj)
 782  ix=ix+nsac(ns)
      ix=0
      write(jtape,9020)
 9020 format(//
     *' normalised input eigenvectors'/)
      loop=1
      do 545 i=1,nrootx
      call dcopy(jsec,cofff(loop),1,doff,1)
      loop=loop+jsec
      corea=0.0d0
      do 554 j=1,jsec
      butzi=doff(j)
554   corea=corea+butzi*butzi
      corea=1.0d0/dsqrt(corea)
      do 555 j=1,jsec
 555   doff(j)=corea*doff(j)
      write(jtape,9021)i
9021  format(/'vector',i2,3x,'coefficients:'/)
      write(jtape,286)(doff(j),j=1,jsec)
 286  format(1x,14f8.5)
      do 545 j=1,lulu
       jj=mn(j)
      iz=ktrk(jj)
        ns=nsc(j)
      ns=nsac(ns)
      do 545 k=1,ns
       iz=iz+1
      ix=ix+1
 545  coff(ix)=doff(iz)
      go to 600
681   do 20010 i=1,nrootx
      isym(i)=isymm(i)
      if (isym(i).lt.1.or.isym(i).gt.jsec) call caserr(
     *'invalid root specified for selection')
20010 continue
      if(.not.oprint(29)) write (jtape,9010) (isym(i),i=1,nrootx)
9010  format(/' input specified roots for use in selection :
     1 ',3x,20i3)
 600  jto=0
      mm=0
      kft=0
      call rewftn(mtype)
      read (ntape) olab8
      call unpack(olab8,8,olab,2000)
      read (mtype) h
      if (ical.eq.2) read (kfile) iqr
      read (nston) mfg,(bkay(i),i=1,mfg),
     *mfg1,(acoul(i),i=1,mfg1),
     *mfg1,(aexc(i),i=1,mfg1),core
      call rewftn(nston)
_IF(cray,ksr,i8)
      read (nston) ibm,zero
_ELSE
      read (nston) ribm,zero
_ENDIF
      call smr130(bkay,acoul,aexc,intmax)
      if(oprint(31))write(jtape,861)
 861  format(//20x,'root hamiltonian matrix'/20x,23('-')/)
      if(oprint(31))call writel(w,jsec)
      if(oprint(31))write(jtape,9001)zero
 9001 format(/
     *' **** note : diagonal elements scaled by ',f16.8,
     *' hartree')
      if(jsec.gt.1) go to 66
      hy(1)=1.0d0
       go to 67
  66    call eigmrd(w,hy,jsec,0)
  67     ib=0
       do 68 i=1,jsec
         ib=ib+i
  68    hs(i)=w(ib)
      ib=jsec+1
      km=jsec*jsec
      do 76 i=1,jsec
      ib=ib-1
      temp(i)=hs(ib)+zero
       lm=km-jsec
      lq=lm+1
      ilc=0
      do 7000 j=lq,km
      ilc=ilc+1
 7000 vred(ilc,i)=hy(j)
      if (ipt0.ge.0) go to 76
      ilc=0
      ico=0
      do 21 j=1,nrootx
      corea=0
      do 1554 k=lq,km
      ilc=ilc+1
1554  corea=corea+coff(ilc)*hy(k)
      corea=dabs(corea)
      if (corea.gt.rmin) ico=j
21    ht(j)=corea
      write(jtape,1555) i,(j,ht(j),j=1,nrootx)
 1555 format(/
     *' overlap between zero order eigenvector no.',i3,
     *' and input eigenvectors'/
     */6(/i3,5x,f8.4))
      jtrk(i)=ico
   76 km=lm
      do 99999 loop=1,jsec
99999 ilifq(loop)=(loop-1)*nrmx
      if(.not.oprint(29)) then
       write(jtape,9000)
 9000  format(/40x,'root eigenvectors'/40x,17('-')/)
       call writem(vred,ilifq,jsec,jsec)
       write(jtape,9003)
 9003  format(///40x,'root eigen values'/40x,17('-')/)
       write(jtape,9002)(temp(k),k=1,jsec)
 9002  format(/10x,8f14.7)
      endif
      if(ical.lt.2) goto 2533
      read(jdeli) trwa,trwb,trwc,trash,tdel,irom,wect
      trash=trash+(4-irom)*tdel
      ilc=jsec*jsec
      do 2534 i=1,jsec
      ilc=ilc-jsec
      ix=0
      w(i)=0.0d0
      do 2535 j=1,nrootx
      corea=0.0d0
      kmax=ilc
      do 2536 k=1,jsec
      kmax=kmax+1
      ix=ix+1
 2536 corea=corea+hy(kmax)*wect(ix)
 2535 w(i)=w(i)+corea*corea
 2534 continue
      write(jtape,861) (w(i),i=1,jsec)
      do 2537 i=1,nrootx
      corea=0.0d0
      do 2538 j=1,jsec
      if(w(j).lt.corea) goto 2538
      corea=w(j)
      k=j
 2538 continue
      w(k)=-2.0d0
 2537 isym(i)=k
      goto 684
 2533 if (ipt0) 683,643,684
 683  do 556 i=1,nrootx
      do 557 j=1,jsec
      if (jtrk(j).eq.i) go to 556
557   continue
      go to 759
556   isym(i)=j
  684 do 493 i=1,nrootx
      ilc=10000
      do 494 j=1,nrootx
      if(isym(j).ge.ilc) goto 494
      ilc=isym(j)
      ix=j
  494 continue
      itym(i)=ilc
  493 isym(ix)=10001
      do 495 i=1,nrootx
  495 isym(i)=itym(i)
      go to 645
  643 do 644 i=1,nrootx
  644 isym(i)=i
  645 write(jtape,559)(isym(i),i=1,nrootx)
 559  format(/1x,
     *'*** selection based on following roots of zero order problem:'/
     *10x,20i3)
      ilc=0
      do 561 i=1,nrootx
      ix=isym(i)
      ix=jsec-ix+1
      ew(i)=hs(ix)
        lm=(ix-1)*jsec
      do 561 j=1,jsec
      ilc=ilc+1
      lm=lm+1
561   vect(ilc)=hy(lm)
      dmax=-50000.0d0
      do 641 i=1,nrootx
      if (ew(i).gt.dmax) dmax=ew(i)
  641 continue
      dmin=50000.0d0
      do 642 i=1,nrootx
      if (ew(i).lt.dmin) dmin=ew(i)
  642 continue
      dmin=dmin-0.2d0
      dmax=dmax+0.2d0
      if(ical.eq.2) goto 576
      write(linf) jsec,nrootx,nytl,nplu,ndub,vect,ew,nconf
      if(ical.eq.3) goto 290
      goto 575
  576 read(linf)
      goto 290
 575  trash=t0
      tdel=t1
      write (jtape,9300) trash,tdel
      trash=trash* 1.0d-06
 9300 format(/' *** threshold specified ***'//
     *' minimal selection threshold ',f7.2,' microhartree'/
     *' threshold increment for use in selection ',f7.2,
     *' microhartree'/)
      tdel=tdel* 1.0d-06
      trwa=0.05d0*trash
      trwb=0.25d0*trash
      trwc=0.50d0*trash
  290 trwd=trash
      trwe=trash+tdel
      trwf=trwe+tdel
      trwg=trwf+tdel
      trwh=trwg+tdel
      trwi=trwh+tdel
      trwj=trwi+tdel
      write (jtape,1600) trwa,trwb,trwc,trwd,trwe,trwf,trwg
     1,trwh,trwi,trwj
 1600 format( ' following threshold classes to be investigated :'
     *//2x,10f11.6)
      do 400 i=1,10
      istm(i)=0
      do 400 j=1,nrootx
  400 trsum(i,j)=0.0d0
      do 287 j=1,mxnshl
  287 itest(j)=0
      i=4
      write(ideli)trwa,trwb,trwc,trwd,tdel,i,vect,ew
      jpag=1
      ipag=1
      ifr=0
      write(jtape,9302)
      if(oprint(33))write(jtape,9303)
 9302 format(/1x,104('-')
     */40x,20('*')/40x,'results of selection'/
     *40x,20('*') )
 9303 format(/
     *' note : configurations are tagged as follows : '//
     *5x,'type m : main configuration'/
     *5x,'type r : retained configuration'/
     *5x,'type s : selected configuration  - a function of threshold'/
     *5x,'type d : discarded configuration - a function of threshold'/)
      if(oprint(33))write(jtape,9502)
9502  format(1x,104('-')//1x,
     *'sequence no.   type  configuration',8x,
     *'computed energy lowering with respect'/1x,
     *'config.  saf',30x,
     *'to specified roots (micro-hartree)')
      if(oprint(33))write(jtape,9505)(isym(l),l=1,nrootx)
 9505 format(/35x,'roots : ',10i8/(43x,10i8) )
_IFN(parallel)
      call setbfc
_ELSE
c **** MPP
      call closbf3
      call setbfc(-1)
c **** MPP
_ENDIF
      do 79 i=1,iswh
      nc=nconf(i)
      noddd=i-1+nplu(i)
      kml=nsac(i)
      if(nc.ne.0.and.oprint(33))write(jtape,9301)i,noddd
 9301 format(1x,104('=')/1x,
     *'**** SAFs in supercategory ',i3/1x,
     *'     (no. of open shells = ',i2,')',20x,'Energy Lowerings'
     * /1x,104('='))
      ik=i+3
      go to (603,605,606,607,608),i
 603  jsl=idra(1)
      ix=nad1
      kss=0
      if (jz-2) 481,482,483
 483  js=jsl
      ks=idra(3)+idrc(3)+jdrc(3)-js
 485  if (kss.eq.ks) go to 484
      ix=ix+nad
      call rdbak3(js)
      js=js+nsz
      call stopbk3
      call unpack(lout,8,ot(ix),nad)
      kss=kss+1
      go to 485
 482  ks=idra(3)+idrc(3)-jsl
      js=jsl
      go to 485
 481  ks=idra(2)+idrc(2)-jsl
      js=jsl
 486  if (kss.eq.ks) go to 6620
      ix=ix+nad
      call rdbak3(js)
      js=js+nsz
      call stopbk3
      call unpack(lout,8,ot(ix),nad)
      kss=kss+1
      go to 486
6620  js=jan(2)
      ks=jbn(2)
 487  if (kss.eq.ks) go to 6621
      ix=ix+nad
      call rdbak3(js)
      js=js+nsz
      call stopbk3
      call unpack(lout,8,ot(ix),nad)
      kss=kss+1
      go to 487
 6621 js=idra(3)
      ks=idrc(3)
 488  if (kss.eq.ks) go to 484
      ix=ix+nad
      call rdbak3(js)
      js=js+nsz
      call stopbk3
      call unpack(lout,8,ot(ix),nad)
      kss=kss+1
      go to 488
 484  if (nc.eq.0) go to 79
      mcl=mconf(i)
      ihp=0
      call upackx(ot)
      if (oprint(31)) then
       write(6,*) 'smrd13 - 484: ic,i2,i3', ic,i2,i3
      endif
      ic7=ic
      i27=i2
      i37=i3
      if (mcl.eq.0) go to 647
      jqn=idrc(1)*nad+1
      jrn=idrc(1)*nad
      call upackx(ot(jqn))
      if (oprint(31)) then
       write(6,*) 'smrd13 - 647: ic,i2,i3', ic,i2,i3
      endif
      ica=ic+jrn
      i2a=i2+jrn
      i3a=i3+jrn
647   j8=3
      mc=mconf(j8)
      if (mc.eq.0) go to 402
      kmj=nsac(j8)
      iej=iez(j8)
      nis=-kmj-kmj+(jan(2)-jsl)*nad
      nzer=kmj*kml
      iix=kmj+kmj+1
 402  j7=2
      kc=mconf(j7)
      if (kc.eq.0) go to 403
      kmk=nsac(j7)
      iek=iez(j7)
      kzer=kmk*kml
      kr=nod(j7)
      mps=nplu(j7)
      mms=1
      jca=kr+mult
      j3a=jca*kmk+1
      jab=jan(1)-jsl
      jaq=jab*nad+1
      jab=jab*nad
      call upackx(ot(jaq))
      if (oprint(31)) then
       write(6,*) 'smrd13 - 403: ic,i2,i3', ic,i2,i3
      endif
      ic1=ic+jab
      i21=i2+jab
      i31=i3+jab
      iiq=kmk+kmk+1
      md1=ndub(j7)
      go to 403
 605  if (nc.eq.0) go to 79
      mcl=mconf(i)
      ihp=idra(2)-jsl
      jqn=ihp*nad+1
      ihp=ihp*nad
      call upackx(ot(jqn))
      if (oprint(31)) then
       write(6,*) 'smrd13 - 605: ic,i2,i3', ic,i2,i3
      endif
      ic7=ic+ihp
      i27=i2+ihp
      i37=13+ihp
      if (mcl.eq.0) go to 648
      jqn=idra(2)+idrc(2)-jsl
      jrn=jqn*nad
      jqn=jqn*nad+1
      call upackx(ot(jqn))
      if (oprint(31)) then
       write(6,*) 'smrd13 - 648: ic,i2,i3', ic,i2,i3
      endif
      ica=ic+jrn
      i2a=i2+jrn
      i3a=i3+jrn
  648 j7=3
      kc=mconf(j7)
      if (kc.eq.0) go to 609
      kmk=nsac(j7)
      iek=iez(j7)
      kzer=kmk*kml
      iiq=kmk+kmk+1
      md1=ndub(j7)
      kr=nod(j7)
      jca=kr+mult
      j3a=jca*kmk+1
      mps=nplu(j7)
      mms=2
      jab=jan(3)-jsl
      jaq=jab*nad+1
      jab=jab*nad
      call upackx(ot(jaq))
      if (oprint(31)) then
       write(6,*) 'smrd13 - 609: ic,i2,i3', ic,i2,i3
      endif
      ic1=ic+jab
      i21=i2+jab
      i31=i3+jab
 609  j6=1
      lc=mconf(j6)
      if (lc.eq.0) go to 403
      jab=jan(1)-jsl
 615  ndl=ndet(j6)
      km2=nsac(j6)
      la=iaz(j6)
      lzer=kml*km2
      lr=nod(j6)
      jcb=nod(i)+mult
      j3b=jcb*kml+1
      lps=nplu(j6)
      lms=j6-1
      jaq=jab*nad+1
      jab=jab*nad
      call upackx(ot(jaq))
      if (oprint(31)) then
       write(6,*) 'smrd13 - 615: ic,i2,i3'
_IF()
       write(6,92124) ndl,km2,la,lzer,lr,jcb,j3b
92124  format(1x,'ndl,km2,la,lzer,lr,jcb,j3b = ',7i8)
       write(6,92125) lps,lms,jaq,jab
92125  format(1x,'lps,lms,jaq,jab = ', 4i8)
       loop1 = jaq
       do loop=1,12
       write(6,92127) loop,ot(loop1)
92127  format(1x,i3,2x,z16)
       loop1 = loop1 + 1
       enddo
       write(6,92123) ic, i2, i3
92123  format(1x,'ic,i2,i3 = ', 3i8)
_ENDIF
      endif
      ic2=ic+jab
      i22=i2+jab
      i32=i3+jab
      go to 403
 606  if (nc.eq.0) go to 79
      mcl=mconf(i)
      ihp=idra(3)-jsl
      if (jz.eq.1) ihp=ihp-jdrc(2)-jbn(3)
      jqn=ihp*nad+1
      ihp=ihp*nad
      call upackx(ot(jqn))
      if (oprint(31)) then
       write(6,*) 'smrd13 - 606: ic,i2,i3', ic,i2,i3
      endif
      ic7=ic+ihp
      i27=i2+ihp
      i37=i3+ihp
      if (mcl.eq.0) go to 611
      jqn=idra(3)+idrc(3)-jsl
      jrn=jqn*nad
      jqn=jqn*nad+1
      call upackx(ot(jqn))
      if (oprint(31)) then
       write(6,*) 'smrd13 - 611: ic,i2,i3', ic,i2,i3
      endif
      ica=ic+jrn
      i2a=i2+jrn
      i3a=i3+jrn
 611  j5=1
      lcl=mconf(j5)
      if(lcl.eq.0) go to 612
      ndi=ndet(j5)
      kmi=nsac(j5)
      jai=iaz(j5)
      nzei=kmi*kml
      nisi=-kml-kml+(jan(2)-jsl)*nad
      if (jz.eq.1) nisi=nisi-(jdrc(2)*nad)
 612  j6=2
      lc=mconf(j6)
      if (lc.eq.0) go to 403
      jab=jan(3)-jsl
      go to 615
  607 if (nc.eq.0) go to 79
      j5=2
      lcl=mconf(j5)
      jsl=jan(4)
      if (jz.eq.3) go to 489
      ks=jbn(4)
      jqn=ks
      js=jsl
      ix=nad1
      kss=0
 490  if (kss.eq.ks) go to 6222
      ix=ix+nad
      call rdbak3(js)
      js=js+nsz
      call stopbk3
      call unpack(lout,8,ot(ix),nad)
      kss=kss+1
      go to 490
 6222     ks=idrc(4)
      js=idra(4)
      kss=0
 491  if (kss.eq.ks) go to 492
      ix=ix+nad
      call rdbak3(js)
      js=js+nsz
      call stopbk3
      call unpack(lout,8,ot(ix),nad)
      kss=kss+1
      go to 491
 489  ks=idra(4)+idrc(4)-jsl
      js=jsl
      ix=nad1
      kss=0
 613  if (kss.eq.ks) go to 6223
      ix=ix+nad
      call rdbak3(js)
      js=js+nsz
      call stopbk3
      call unpack(lout,8,ot(ix),nad)
      kss=kss+1
      go to 613
 6223     jqn=jbn(4)+jbn(5)
  492 ihp=jqn*nad
      jqn=jqn*nad+1
      call upackx(ot(jqn))
      if (oprint(31)) then
       write(6,*) 'smrd13 - 492: ic,i2,i3', ic,i2,i3
      endif
      ic7=ic+ihp
      i27=i2+ihp
      i37=i3+ihp
      if (lcl.eq.0) go to 614
      ndi=ndet(j5)
      kmi=nsac(j5)
      jai=iaz(j5)
      nzei=kmi*kml
      nisi=-kml-kml
 614  j6=3
      lc=mconf(j6)
      if (lc.eq.0) go to 403
      jab=jbn(4)
      go to 615
  608 if (nc.eq.0) go to 79
      jsl=jan(6)
      ks=jbn(6)
      js=jsl
      ix=nad1
      kss=0
 616  if (kss.eq.ks) go to 6224
      ix=ix+nad
      call rdbak3(js)
      js=js+nsz
      call stopbk3
      call unpack(lout,8,ot(ix),nad)
      kss=kss+1
      go to 616
 6224     ksl=idra(5)
      ks=idrc(5)
      js=ksl
      kss=0
  617  if (kss.eq.ks) go to 6225
      ix=ix+nad
      call rdbak3(js)
      js=js+nsz
      call stopbk3
      call unpack(lout,8,ot(ix),nad)
      kss=kss+1
      go to 617
 6225     j5=3
      lcl=mconf(j5)
      if (lcl.eq.0) go to 618
      ndi=ndet(j5)
      kmi=nsac(j5)
      jai=iaz(j5)
      nzei=kmi*kml
      nisi=-kml-kml
 618  ihp=jbn(6)
      jqn=ihp*nad+1
      ihp=ihp*nad
      call upackx(ot(jqn))
      if (oprint(31)) then
       write(6,*) 'smrd13 - 618: ic,i2,i3', ic,i2,i3
      endif
      ic7=ic+ihp
      i27=i2+ihp
      i37=i3+ihp
 403  ndt=ndet(i)
      ic=ic7
      i2=i27
      i3=i37
      ia=iaz(i)
      ie=iez(i)
      mzer=kml*kml
      nd=ndub(i)
      nps=nplu(i)
      nms=i-1
      niw=nps*nms
      nr=nod(i)
      nx=nytl(i)
      n1=nr+1
      n2=nr+2
      n3=nr-1
      np1=nps-1
      np2=np1+np1
      nd1=nd-1
      nm1=i-2
      nm2=nm1+nm1
      j2=kml+1
      jc=nr+nr
      j3=jc*kml+1
      iiz=kml+kml+1
      ii=iiz+kml
      lpix=0
      read (mtape) mkon
      jg=0
      if (ical.eq.2) go to 619
      np=kml*ndt
      ix=ia+ndt
      do 620 l=1,np
      ix=ix+1
 620  a(l)=ae(ix)
      write (linf) ndt,kml,a
      ix=ie
      do 622 l=1,mzer
      ix=ix+1
 622  e(l)=ae(ix)
      ix=ihp
      do 623 l=1,kml
      ix=ix+1
  623 ihog(l)=oihog(ix)
      if (i.eq.1) go to 624
      if (i.gt.2) go to 625
      do 426 l=1,ndt
 426  imap(l)=l
      go to 624
 625  np=ndt*nms
      do 427 l=1,nms
      iper(l)=l
 427  imap(l)=l
      ib=nms
 628  jb=nms
      jm=nr
 629  km=iper(jb)
      if (km.ne.jm) go to 630
      jb=jb-1
      if (jb.eq.0) go to 631
      jm=jm-1
      go to 629
 630  ip=iper(jb)+1
      iper(jb)=ip
      if (jb.eq.nms) go to 432
      jb=jb+1
      do 633 l=jb,nms
      ip=ip+1
  633 iper(l)=ip
  432 do 435 l=1,nms
      ib=ib+1
  435 imap(ib)=iper(l)
      go to 628
  631 if (ib.ne.np) go to 755
  624 write (linf) ihog,imap,e
      go to 436
 619  read(linf)
      read (linf)
 436  call smr131(bkay,acoul,aexc,intmax)
      if (ical.ne.2) go to 469
      read (linf)
      go to 79
469   lpix=lpix+1
      jkan(lpix)=0
      write(linf) jkan
79    continue
      iqr(kft+1)=0
      if (ical.ne.2) write (kfile) iqr
      immm=ipag-1
      write(jtape,429)
 429  format(//22x,
     *'summary of energy lowerings in various threshold classes'
     1/22x,56('*')/)
      crus=1000.0d0
      bkay(1)=trwa*crus
      bkay(2)=trwb*crus
      bkay(3)=trwc*crus
      bkay(4)=trwd*crus
      bkay(5)=trwe*crus
      bkay(6)=trwf*crus
      bkay(7)=trwg*crus
      bkay(8)=trwh*crus
      bkay(9)=trwi*crus
      bkay(10)=trwj*crus
      write(jtape,430)  (bkay(i),i=1,10)
  430 format(2x,16hthreshold ranges,3x,10f11.8/)
      do 440 i=1,4
  440 immm=immm+istm(i)
      istm(1)=immm-istm(1)
      do 441 i=2,10
  441 istm(i)=istm(i-1)-istm(i)
      write(jtape,431) immm,(istm(i),i=1,10)
  431 format(2x,' saf totals',1x,i7,i10,9i11/)
  433 format(1x,' energy sums root', i3,10f11.5/)
      do 442 i=1,nrootx
      bkay(1)=trsum(1,i)*crus
      do 443 j=2,10
      j1=j-1
      trsum(j,i)=trsum(j1,i)+trsum(j,i)
  443 bkay(j)=trsum(j,i)*crus
  442 write(jtape,433) i,(bkay(j),j=1,10)
      if(ical.ne.2) goto 291
      irom=istm(irom)
      irom=(irom-1)/irsf+1
      do 292 i=1,irom
      read(jdeli) senk
  292 write(ideli) senk
      goto 428
  291 if (ifr.eq.0) go to 447
      write(ideli) senk
  447 if (istm(4).le.imsec) go to 428
      call smr133(sx,trwe,trwf,trwg,trwh,trwi,trwj,imsec)
      return
428   write(ideli) nrootx,((trsum(j,i),j=1,10),i=1,nrootx),istm
      call clredx
      return
 759  write(jtape,765) i,nrmap
  765 format(/10x, 'desired root not found on basis of overlap criterio
     1n',2i6)
      call errors(530)
755   write (jtape,756) ib,ndt,nms
756   format(10x,'error in imap generation',3i6)
      call errors(504)
      return
      end
      subroutine smr131(bkay,acoul,aexc,intmax)
      implicit REAL (a-h,o-z), integer(i-n)
      integer olab,ot,oiot,oihog
_IF(cray,ksr,i8)
      integer olab8
_ENDIF
      character *1 typet,typem,typer,types,typed
      character * 3 itag,iblnk,char3i
INCLUDE(common/sizes)
INCLUDE(common/prints)
      common/linkmr/ ic,i2,i3
      common/bufd/ gout(510),nword
      common/blksi3/nsz
      common/scrtch/ndet(5),nsac(5),iaz(5),iez(5),idra(5),
     + idrc(5),jdrc(5),
     + jan(7),jbn(7),ispa2(7),ae(7829),ot(44000),ideks(maxorb),
     + olab(2000),ij(8),nytl(5),nplu(5),ndub(5),nod(5),jkan(mxcrec),
     + nir(maxorb),loc(maxorb),sac(10),
_IF(cray,ksr,i8)
     + olab8(250),f(2304),h(mxcrec),isc(5),bob(3),itest(mxnshl),
_ELSE
     + olab8(250),f(2304),h(mxcrec),isc(5),isd,bob(3),itest(mxnshl),
_ENDIF
     + core,
     + w(mxcsf*(mxcsf+1)/2),
     + ew(mxroot),vect(mxcsf*mxroot),trsum(10,mxroot),
     * wect(mxcsf*mxroot),senk(500),iqr(mxcrec),istm(10)
     *,ht(mxroot),hs(mxcsf),ea(mxroot),iper(10),
     * coff(mxcsf*mxroot),doff(mxcsf),nrec(5),jtrk(mxcsf),
     * isym(mxroot+1),mkon(mxcrec),ikan(mxcrec),hy(mxcsf*mxcsf)
     *,ktrk(mxcsf)
      common/scra/vred(mxcsf,mxcsf),temp(mxcsf),ilifq(mxcsf)
      common/craypk/icall(14+mxcsf*mxnshl+mxcsf+mxcsf),isymm(mxroot),
     *cofff(mxcsf*mxroot),t0,t1
      common/jany/ jerk(10),jbnk(10),jbun(10)
      common /b/ iwod,nid,lg,irsf,ntape,ntype,ical,ntab,nc,lt,lz,nx,
     + nd,nston,itape,jtape,mdisk,ideli,mtape,iswh,nr,mm,jblk,jto,
     + igmax,nr1,mtype,linf,kfile,
     + nawk,isec,nshl,
     + nsec,lulu,nrootx,lsng,nprin,ipt0,nconf(5),mconf(5),mh,jswh,
     + id(mxcsf),mn(mxcsf),nsc(mxcsf),
     + imo,mc,mz,ndx,nrx,ntx,ny,nt,iht,ihz,nz,mt,nzx,jy,
     + jx,nd1,nd2,md1,md2,nr3,nr4,mr1,mr2,nt1,nt2,mt1,mt2,ih1,ih2,ih3,
     + ih4,jh2,jh3,jh4,iab,nt3,iac,nr2,m,nad,nad1,ifr1, mx,md,mr
     +,ipag
       common /d/ jbab,kml,ifl,ii,nis,iiz,nzer,kmj,ndj,jsac,ie,iej,
     * ja,kc,kmk,ndk,kr,iiq,ic1,ixex,i21,i2ex,i31,i3ex,iix,
     * kzer,ndt,iek,j3,j3a,mps,nps,mms,nms,jc,jca,lcl,ndi,kmi,
     * ia,jai,nzei,nisi,lc,kk,ndl,km2,la,lzer,lr,lps,lms,
     * ic2,i22,i32,jcb,j3b
      common /bx/rmin,smin,crus,crud,crub,sch,jsec,kmax
      common/lsort/isppp(16),
     * trwa,trwb,trwc,trwd,trwe,trwf,trwg,trwh,trwi,trwj,
     * dmin,dmax,jg,ik,mzer,jpag,ica,
     * n1,n2,n3,nm1,nm2,np1,np2,niw,
     * kft,lpix,ifr,ihp,i3a,i2a,j2,i
      dimension bkay(intmax),acoul(intmax),aexc(intmax)
      dimension itag(mxnshl),iea(mxroot)
      dimension a(5292),e(2304),lout(510),
     1hp(4800),xemp(500),itym(mxroot),
     2ihog(48),imap(504),j9(56),oihog(48),
     3oiot(10962),t(50000)
_IF(cray,ksr,i8)
      dimension af(7885)
_ELSE
      dimension af(7857)
_ENDIF
      equivalence (t(1),j9(1)),(j9(1),ndet(1)),(gout(1),lout(1)),
     6 (t(1),af(1)),(ew(1),itym(1)),
     3 (ot(13),oihog(1)),(oihog(1),oiot(1)),
     4 (t(1),xemp(1)),(a(1),hy(1)),(a(1),e(1)),(hp(1),a(1)),
_IF(cray,ksr,i8)
     5 (a(2305),imap(1)),(a(2809),ihog(1))
_ELSE
     5 (a(2305),imap(1)),(a(2557),ihog(1))
_ENDIF
c     character *1 dash,typez
c     data dash,typez/'-','z'/
      data iblnk,fact6/'   ',10.0d6/
      data typet,typem,typer,types,typed/
     *'t','m','r','s','d'/
      do 573 j=1,nc
      do 9500 k=1,mxnshl
 9500 itag(k)=iblnk
      kk = 0
      do 574 k=1,nx
      jg=jg+1
      itest(k)=mkon(jg)
      if(itest(k).eq.0)go to 9501
      kk=kk+1
      itag(kk)=char3i(itest(k))
9501  if (jg.lt.iwod) go to 574
      jg=0
      read (mtape) mkon
  574 continue
      mm=mm+1
      kk=olab(mm)
      if (mm.lt.nid) go to 437
      mm=0
      read (ntape) olab8
      call unpack(olab8,8,olab,2000)
 437  if (kk-6) 438,439,640
 439  corea=crud
702   format(2i7,4x,a1,3x,24a3/22x,24a3)
      if(oprint(33))write(jtape,702) jpag,ipag,typem,(itag(ll),ll=1,nx)
      go to 1
 640  corea=crus
      if(oprint(31))
     *write(jtape,702)jpag,ipag,typet,(itag(ll),ll=1,nx)
      go to 1
 438  ifl=0
      jbab=0
      do 850 k=1,jswh
      md=mconf(k)
      if (md.eq.0) go to 850
      kfc=ik-k
      go to (851,852,853,854,855), kfc
      mz=nsac(k)
      do 856 ix=1,md
      kbab=jbab
      jbab=jbab+mz
      ig7=kbab-jsec
      do 856 l=1,kml
      ig7=ig7+jsec
      in3=ig7
      do 856 ll=1,mz
      in3=in3+1
 856  hp(in3)=0
      go to 850
 853  do 130 lb=1,md
      kbab=jbab
      jbab=jbab+kml
      if (ifl.eq.1) go to 857
      ifl=1
      go to 131
 857  mm=mm+1
      kk=olab(mm)
_IF1()c     if(j.gt.7100) write(6,*) 'smrd13: ',md,kml,jbab,mm,kk,lb
      if (mm.lt.nid) go to 131
      mm=0
      read(ntape)olab8
      call unpack(olab8,8,olab,2000)
 131  go to (858,132,133,134,135),kk
 135  mm=mm+1
      kk=olab(mm)
      if(mm.lt.nid) go to 136
      mm=0
      read(ntape)olab8
      call unpack(olab8,8,olab,2000)
136   mm=mm+1
      ll=olab(mm)
      if(mm.lt.nid) go to 137
      mm=0
      read(ntape)olab8
      call unpack(olab8,8,olab,2000)
137   kk=min(kk,ll)+ideks(max(kk,ll))
      call vclr(f,1,mzer)
_IF1()c     do 140 l=1,mzer
_IF1()c140  f(l)=0
      tim=aexc(kk)
_IF1(iv)      jx=-kml
      do 141 l=1,kml
_IF1(iv)      kx=oihog(l+ihp)+ia
_IFN1(iv)      kx=oihog(l+ihp)+ia+ndt
_IF1(iv)      jx=jx+1
_IF1(iv)      lx=jx
_IF1(iv)      do 141 m=1,kml
_IF1(iv)      lx=lx+kml
_IF1(iv)      kx=kx+ndt
_IF1(iv)141   f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv)141   call daxpy(kml,tim,ae(kx),ndt,f(l),kml)
146   ig7=kbab-jsec
      mx=ie
_IF1()c      do 145 l=1,kml
_IF1()c      ig7=ig7+jsec
_IF1()c      in3=ig7
_IF1()c      ly=0
_IF1()c      my=mx
_IF1()c      do 145 m=1,kml
_IF1()c      tim=0
_IF1()c      mx=my
_IF1()c      do 143 if=1,kml
_IF1()c      ly=ly+1
_IF1()c      mx=mx+1
_IF1()c 143  tim=tim+f(ly)*ae(mx)
_IF1()c      in3=in3+1
_IF1()c 145  hp(in3)=tim
_IF1(c)      call mxma(ae(mx+1),kml,1,f,1,kml,hp(kbab+1),jsec,1,kml,kml,kml)
_IFN1(c)       call mxmaa(ae(mx+1),kml,1,f,1,kml,hp(kbab+1),jsec,1,kml,kml,kml)
      go to 130
 144  mx=ie
_IF1()c      ig7=kbab-jsec
_IF1()c      do 253 l=1,kml
_IF1()c      ly=0
_IF1()c      ig7=ig7+1
_IF1()c      my=mx
_IF1()c      in3=ig7
_IF1()c      do 253 m=1,kml
_IF1()c      in3=in3+jsec
_IF1()c      tim=0
_IF1()c      mx=my
_IF1()c      do 254 if=1,kml
_IF1()c      ly=ly+1
_IF1()c      mx=mx+1
_IF1()c 254  tim=tim+f(ly)*ae(mx)
_IF1()c 253  hp(in3)=tim
_IF1(c)      call mxma(ae(ie+1),kml,1,f,1,kml,hp(kbab+1),1,jsec,kml,kml,kml)
_IFN1(c)      call mxmaa(ae(ie+1),kml,1,f,1,kml,hp(kbab+1),1,jsec,kml,kml,kml)
      go to 130
 858  ig7=kbab-jsec
      do 459 mz=1,kml
      ig7=ig7+jsec
      in3=ig7
      do 459 nz=1,kml
      in3=in3+1
 459  hp(in3)=0
      go to 130
 132  mm=mm+1
      ll=olab(mm)
      if(mm.lt.nid) go to 147
      mm=0
      read(ntape)olab8
      call unpack(olab8,8,olab,2000)
147   mm=mm+1
      jj=olab(mm)
      if(mm.lt.nid) go to 148
      mm=0
      read(ntape)olab8
      call unpack(olab8,8,olab,2000)
148   mm=mm+1
      kk=olab(mm)
      if(mm.lt.nid) go to 149
      mm=0
      read(ntape)olab8
      call unpack(olab8,8,olab,2000)
149   if(jj.lt.kk) go to 150
      jj=ideks(jj)+kk
      ibob=0
      go to 151
150   jj=ideks(kk)+jj
      ibob=1
151   jj=jj*ii+ica
      jto=jto+1
      ip=ot(jj)
      if(ip.eq.255) go to 830
      coul=h(jto)
      if(jto.lt.iwod) go to 152
      jto=0
      read(mtype) h
152   jto=jto+1
      exc=-h(jto)
      if(jto.lt.iwod) go to 153
      jto=0
      read(mtype) h
      go to 153
830   coul=-h(jto)
      if (jto.lt.iwod) go to 831
      jto=0
      read (mtype) h
831   jto=jto+1
      exc=h(jto)
      if(jto.lt.iwod) go to 153
      jto=0
      read (mtype) h
  153 if (ll-2) 832,833,834
  832 bob(1)=coul+exc
      bob(2)=exc
      bob(3)=coul
      go to 835
  833 bob(1)=exc
      bob(2)=-coul
      bob(3)=coul+exc
      go to 835
 834  bob(1)=-exc
      bob(2)=-coul-exc
      bob(3)=coul
835   continue
_IF1()c     do 154 l=1,mzer
_IF1()c154  f(l)=0
      call vclr(f,1,mzer)
_IF1(iv)      mx=-kml
      do 155 l=1,kml
      jj=jj+1
_IF1(iv)      mx=mx+1
      ip=ot(jj)
      if(ip.eq.2)  go to 159
      jj=jj+1
      kx=ot(jj)
_IF1(iv)      kx=kx+ia
_IFN1(iv)      kx=kx+ia+ndt
      jj=jj+1
      tim=bob(1)
_IF1(iv)      lx=mx
_IF1(iv)      do 161 m=1,kml
_IF1(iv)      kx=kx+ndt
_IF1(iv)      lx=lx+kml
_IF1(iv) 161  f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv)      call daxpy(kml,tim,ae(kx),ndt,f(l),kml)
      go to 155
159   jj=jj+1
      kx=ot(jj)
_IF1(iv)      kx=kx+ia
_IFN1(iv)      kx=kx+ia+ndt
      tim=bob(2)
_IF1(iv)      lx=mx
_IF1(iv)      do 160 m=1,kml
_IF1(iv)      kx=kx+ndt
_IF1(iv)      lx=lx+kml
_IF1(iv) 160  f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv)      call daxpy(kml,tim,ae(kx),ndt,f(l),kml)
      jj=jj+1
      kx=ot(jj)
_IF1(iv)      kx=kx+ia
_IFN1(iv)      kx=kx+ia+ndt
      tim=bob(3)
_IF1(iv)      lx=mx
_IF1(iv)      do 167 m=1,kml
_IF1(iv)      lx=lx+kml
_IF1(iv)      kx=kx+ndt
_IF1(iv) 167  f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv)      call daxpy(kml,tim,ae(kx),ndt,f(l),kml)
155   continue
      if (ibob.eq.1) go to 144
      go to 146
 133  mm=mm+1
      ll=olab(mm)
      if (mm.lt.nid) go to 168
      mm=0
      read(ntape)olab8
      call unpack(olab8,8,olab,2000)
168   mm=mm+1
      jj=olab(mm)
      if(mm.lt.nid) go to 169
      mm=0
      read(ntape)olab8
      call unpack(olab8,8,olab,2000)
169   if(jj.gt.ll) goto 170
      jj=ideks(ll)+jj
      ibob=0
      go to 171
170   jj=ideks(jj)+ll
      ibob=1
171   jto=jto+1
      tim=h(jto)
      if(jto.lt.iwod) go to 172
      jto=0
      read(mtype)h
172   jj=jj*j2+i2a
_IF1()c     do 177 l=1,mzer
_IF1()c177  f(l)=0
      call vclr(f,1,mzer)
_IF1(iv)      jx=-kml
      ip=ot(jj)
      if(ip.eq.255) tim=-tim
      do 173 l=1,kml
_IF1(iv)      jx=jx+1
      jj=jj+1
      kx=ot(jj)
_IF1(iv)      kx=kx+ia
_IFN1(iv)      kx=kx+ia+ndt
_IF1(iv)      lx=jx
_IF1(iv)      do 173 m=1,kml
_IF1(iv)      kx=kx+ndt
_IF1(iv)      lx=lx+kml
_IF1(iv) 173  f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv)173   call daxpy(kml,tim,ae(kx),ndt,f(l),kml)
      if (ibob.eq.1) go to 144
      go to 146
 134  mm=mm+1
      ll=olab(mm)
      if(mm.lt.nid) go to 179
      mm=0
      read(ntape)olab8
      call unpack(olab8,8,olab,2000)
179   mm=mm+1
      jj=olab(mm)
      if(mm.lt.nid) go to 180
      mm=0
      read(ntape)olab8
      call unpack(olab8,8,olab,2000)
180   mm=mm+1
      ntar=olab(mm)
      if(mm.lt.nid) go to 302
      mm=0
      read(ntape)olab8
      call unpack(olab8,8,olab,2000)
302   mm=mm+1
      nb=olab(mm)
      if(mm.lt.nid) go to 303
      mm=0
      read(ntape)olab8
      call unpack(olab8,8,olab,2000)
303   mm=mm+1
      mb=olab(mm)
      if (mm.lt.nid) go to 48
      mm=0
      read (ntape) olab8
      call unpack(olab8,8,olab,2000)
  48  nb=ideks(nb)+mb+ij(ntar)
      sm=bkay(nb)
      if(ll.gt.128) go to 181
      if(ll.lt.jj) go to 182
      ibob=0
      jj=ideks(ll)+jj
      go to 205
182   jj=ideks(jj)+ll
      ibob=1
205   jj=jj*j3+i3a
      if(nr.eq.1) go to 185
      do 186 l=1,n3
      jto=jto+1
      sm=sm+h(jto)
      if(jto.lt.iwod) go to 187
      jto=0
      read(mtype) h
187   jto=jto+1
      sac(l)=-h(jto)
      if(jto.lt.iwod) go to 186
      jto=0
      read(mtype) h
186   continue
185   if(nd.eq.0) go to 188
      do 189 l=1,nd
      jto=jto+1
      tim=h(jto)
      sm=sm+tim+tim
      if(jto.lt.iwod) go to 190
      jto=0
      read(mtype) h
190   jto=jto+1
      sm=sm-h(jto)
      if(jto.lt.iwod) go to 189
      jto=0
      read(mtype) h
189   continue
188     ip=ot(jj)
      if(ip.eq.1) go to 191
      sm=-sm
      if (nr.eq.1) go to 191
      do 192 l=1,n3
192   sac(l)=-sac(l)
191   jj=jj+1
_IF1()c     do 193 l=1,mzer
_IF1()c193  f(l)=0
      call vclr(f,1,mzer)
_IF1(iv)      jx=-kml
      do 194 l=1,kml
      in=jj
      im=ot(in)
_IF1(iv)      jx=jx+1
      in=in+1
      kx=ot(in)
_IF1(iv)      kx=kx+ia
_IFN1(iv)      kx=kx+ia+ndt
_IF1(iv)      lx=jx
      tim=sm
      if(im.eq.2) go to 195
      if(nps.eq.1) go to 196
      do 197 m=1,np1
      in=in+2
      nn=ot(in)
197   tim=tim+sac(nn)
196   continue
_IF1(iv)      do 198 m=1,kml
_IF1(iv)      lx=lx+kml
_IF1(iv)      kx=kx+ndt
_IF1(iv) 198  f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv)      call daxpy(kml,tim,ae(kx),ndt,f(l),kml)
      if(nms.eq.0) go to 194
      do 199 m=1,nms
      in=in+1
      kx=ot(in)
_IF1(iv)      kx=kx+ia
_IFN1(iv)      kx=kx+ia+ndt
_IF1(iv)      lx=jx
      in=in+1
      mx=ot(in)
      tim=sac(mx)
_IF1(iv)      do 199 if=1,kml
_IF1(iv)      lx=lx+kml
_IF1(iv)      kx=kx+ndt
_IF1(iv) 199  f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv)199   call daxpy(kml,tim,ae(kx),ndt,f(l),kml)
      go to 194
195   if(nms.eq.1) go to 200
      do 201 m=1,nm1
      in=in+2
      nn=ot(in)
201   tim=tim+sac(nn)
200   continue
_IF1(iv)      do 202 m=1,kml
_IF1(iv)      lx=lx+kml
_IF1(iv)      kx=kx+ndt
_IF1(iv) 202  f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv)      call daxpy(kml,tim,ae(kx),ndt,f(l),kml)
      do 203 m=1,nps
      in=in+1
      kx=ot(in)
_IF1(iv)      kx=kx+ia
_IFN1(iv)      kx=kx+ia+ndt
_IF1(iv)      lx=jx
      in=in+1
      mx=ot(in)
      tim=sac(mx)
_IF1(iv)      do 203 if=1,kml
_IF1(iv)      lx=lx+kml
_IF1(iv)      kx=kx+ndt
_IF1(iv) 203  f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv)203   call daxpy(kml,tim,ae(kx),ndt,f(l),kml)
194   jj=jj+jc
      if (ibob.eq.1) go to 144
      go to 146
181   ll=256-ll
      if(ll.lt.jj) go to 204
      ibob=0
      jj=ideks(ll)+jj
      go to 206
204   jj=ideks(jj)+ll
      ibob=1
206   jj=jj*j3+i3a
      jto=jto+1
      sm=sm+h(jto)
      if(jto.lt.iwod) go to 207
      jto=0
      read(mtype) h
207   jto=jto+1
      sm=sm+h(jto)
      if(jto.lt.iwod) go to 208
      jto=0
      read(mtype) h
208   if(nr.eq.1) go to 209
      do 211 l=1,n3
      jto=jto+1
      sm=sm+h(jto)
      if(jto.lt.iwod) go to 210
      jto=0
      read(mtype) h
210   jto=jto+1
      sac(l)=-h(jto)
      if(jto.lt.iwod) go to 211
      jto=0
      read(mtype) h
211   continue
209   if(nd.eq.1) go to 212
      do 213 l=1,nd1
      jto=jto+1
      tim=h(jto)
      sm=sm+tim+tim
      if(jto.lt.iwod) go to 214
      jto=0
      read(mtype) h
214   jto=jto+1
      sm=sm-h(jto)
      if (jto.lt.iwod) go to 213
      jto=0
      read(mtype) h
213   continue
 212     ip=ot(jj)
      if(ip.eq.1) go to 215
      sm=-sm
      if(nr.eq.1) go to 215
      do 216 l=1,n3
216   sac(l)=-sac(l)
215   jj=jj+1
_IF1()c     do 217 l=1,mzer
_IF1()c217  f(l)=0
      call vclr(f,1,mzer)
_IF1(iv)      jx=-kml
      do 218 l=1,kml
      in=jj
      im=ot(in)
_IF1(iv)      jx=jx+1
      in=in+1
      kx=ot(in)
_IF1(iv)      kx=kx+ia
_IFN1(iv)      kx=kx+ia+ndt
_IF1(iv)      lx=jx
      tim=-sm
      if (im.eq.2) go to 219
      if(nms.eq.0) go to 220
      in=in+np2
      jn=in
      do 221 m=1,nms
      jn=jn+2
      nn=ot(jn)
221   tim=tim-sac(nn)
 220  continue
_IF1(iv)      do 222 m=1,kml
_IF1(iv)      lx=lx+kml
_IF1(iv)      kx=kx+ndt
_IF1(iv) 222  f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv)      call daxpy(kml,tim,ae(kx),ndt,f(l),kml)
      if(nms.eq.0) go to 218
      do 223 m=1,nms
      in=in+1
      kx=ot(in)
_IF1(iv)      kx=kx+ia
_IFN1(iv)      kx=kx+ia+ndt
_IF1(iv)      lx=jx
      in=in+1
      mx=ot(in)
      tim=sac(mx)
_IF1(iv)      do 223 if =1,kml
_IF1(iv)      lx=lx+kml
_IF1(iv)      kx=kx+ndt
_IF1(iv) 223  f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv)223   call daxpy(kml,tim,ae(kx),ndt,f(l),kml)
      go to 218
219   in=in+nm2
      jn=in
      do 224 m=1,nps
      jn=jn+2
      nn=ot(jn)
224   tim=tim-sac(nn)
_IF1(iv)      do 225 m=1,kml
_IF1(iv)      lx=lx+kml
_IF1(iv)      kx=kx+ndt
_IF1(iv) 225  f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv)      call daxpy(kml,tim,ae(kx),ndt,f(l),kml)
      do 226 m=1,nps
      in=in+1
      kx=ot(in)
_IF1(iv)      kx=kx+ia
_IFN1(iv)      kx=kx+ia+ndt
_IF1(iv)      lx=jx
      in=in+1
      mx=ot(in)
      tim=sac(mx)
_IF1(iv)      do 226 if=1,kml
_IF1(iv)      lx=lx+kml
_IF1(iv)      kx=kx+ndt
_IF1(iv) 226  f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv)226   call daxpy(kml,tim,ae(kx),ndt,f(l),kml)
218   jj=jj+jc
      if (ibob.eq.1) go to 144
      go to 146
130   continue
      go to 850
851   call sm1310
      go to 850
855   call sm1311
      go to 850
852   call sm1312(bkay,intmax)
      go to 850
854   call sm1313(bkay,intmax)
850   continue
_IF1()c     do 105 l=1,mzer
_IF1()c105  f(l)=0
      call vclr(f,1,mzer)
      care=core
      if (nr.eq.0) go to 106
      do 107 l=1,nr
      mal=itest(l)
      ntar=nir(mal)
      mal=loc(mal)+1
      mal=ideks(mal)+ij(ntar)
107   care=care+bkay(mal)
      if(nd.eq.0) go to 108
106   do 109 l=n1,nx
      mal=itest(l)
      ntar=nir(mal)
      nal=ideks(mal+1)
      mal=loc(mal)+1
      mal=ideks(mal)+ij(ntar)
      sm=bkay(mal)
109   care=care+sm+sm+acoul(nal)
108   if(nr.lt.2) go to 110
      do 111 l=2,nr
      mal=itest(l)
      mal=ideks(mal)
      it=l-1
      do 111 m=1,it
      nal=mal+itest(m)
111   care=care+acoul(nal)
110   if(nd.lt.2) go to 112
      sm=0
      do 113 l=n2,nx
      mal=itest(l)
      mal=ideks(mal)
      it=l-1
      do 113 m=n1,it
      nal=mal+itest(m)
      tim=acoul(nal)
113   sm=sm+tim+tim-aexc(nal)
      care=care+sm+sm
112   if(nr.eq.0.or.nd.eq.0) go to 114
      do 115 l=1,nr
      mal=itest(l)
      do 115 m=n1,nx
      nal=itest(m)
      nal=min(mal,nal)+ideks(max(mal,nal))
      sm=acoul(nal)
115   care=care+sm+sm-aexc(nal)
114   kk=ic
      kn=i2
      ml=i3
_IF1(iv)      jx=-kml
      do 118 l=1,kml
_IF1(iv)      jx=jx+1
_IF1(iv)      lx=jx
_IF1(iv)      kx=oihog(l+ihp)+ia
_IFN1(iv)      kx=oihog(l+ihp)+ia+ndt
      tim=care
      if(nps.lt.2) go to 122
      mx=kn+1
      do 301 if=2,nps
      mx=mx+1
      it=if-1
      mal=oiot(mx)
      mal=itest(mal)
      mal=ideks(mal)
      mz=kn
      do 301 in=1,it
      mz=mz+1
      nal=oiot(mz)
      nal=itest(nal)+mal
301   tim=tim-aexc(nal)
      kn=mx
      if(nms.lt.2) go to 122
      mx=ml+1
      do 124 if=2,nms
      mx=mx+1
      mal=oiot(mx)
      mal=itest(mal)
      mal=ideks(mal)
      mz=ml
      it=if-1
      do 124 in=1,it
      mz=mz+1
      nal=oiot(mz)
      nal=itest(nal)+mal
124   tim=tim-aexc(nal)
      ml=mx
 122  continue
_IF1(iv)      do 300 if=1,kml
_IF1(iv)      lx=lx+kml
_IF1(iv)      kx=kx+ndt
_IF1(iv) 300  f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv)      call daxpy(kml,tim,ae(kx),ndt,f(l),kml)
      if(niw.eq.0) go to 118
      do 120 m=1,niw
      kk=kk+1
      kx=oiot(kk)+ia
      kk=kk+1
      mal=oiot(kk)
      kk=kk+1
      nal=oiot(kk)
      if (k.eq.0.and.oprint(33)) write (jtape,241) kx,mal,nal,kk
_IFN1(iv)      kx=kx+ndt
      mal=itest(mal)
      nal=ideks(mal)+itest(nal)
      tim=-aexc(nal)
_IF1(iv)      lx=jx
 241  format(2x,21i6)
_IF1(iv)      do 120 if=1,kml
_IF1(iv)      lx=lx+kml
_IF1(iv)      kx=kx+ndt
_IF1(iv) 120  f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv)120   call daxpy(kml,tim,ae(kx),ndt,f(l),kml)
_IF1(iv)      if (k.eq.0.and.oprint(33)) write (6,239) lx,kx,mal,nal,tim,f(lx)
_IF1(iv)     * ,ae(kx)
_IF1(iv)239   format(2x,4i6,3f20.8)
118   continue
      ly=ie
      mxq=0
      mz=0
      do 1558 k=1,nrootx
1558  ea(k)=0
      do 670 k=1,kml
      eigv=0
      do 902 l=1,kml
      ly=ly+1
      mz=mz+1
 902  eigv=eigv+f(mz)*ae(ly)
      if (k.gt.1) go to 184
      if (ical.ne.2) go to 305
      kft=kft+1
      kx=iqr(kft)
      if (kft.lt.iwod) go to 462
      kft=0
      read (kfile) iqr
 462  if (kx.eq.12) go to 381
 305  if (eigv.lt.dmin) go to 183
      if (eigv.lt.dmax) go to 382
      ilse=0
      go to 184
  183 ilse=1
       goto 184
 382  if (ical.ne.2)  go to 465
      ilse=0
      if(oprint(33))
     *write (jtape,466) ipag,jpag,i,kml,nr,eigv,(itest(ll),ll=1,nx)
 466  format(2x,'z',2i5,i2,2i3,3x,f12.8,2x,24i3)
  184 ilc=0
      do 656 l=1,nrootx
      ll=mxq
      dif=0
      do 399 mw=1,jsec
      ll=ll+1
      ilc=ilc+1
 399  dif=dif+hp(ll)*vect(ilc)
 656  ea(l)=ea(l)+(dif*dif)/(eigv-ew(l))
 670  mxq=ll
      dif=0
         if(ilse.eq.0) go to 880
      do 881 l=1,nrootx
       if(ea(l).lt.dif) dif=ea(l)
 881    continue
        dif=-dif
       go to 882
880     do 121 l=1,nrootx
      if (ea(l).gt.dif) dif=ea(l)
121   continue
882   dif=dif/kml
      do 9504 l=1,nrootx
 9504 iea(l)=idint(ea(l)*fact6)
      if (ical.ne.2) go to 463
      if (kx.eq.11) go to 414
      istm(kx)=istm(kx)+kml
      do 306 l=1,nrootx
 306  trsum(kx,l)=trsum(kx,l)+ea(l)
      if (kx.lt.5) go to 416
      go to 414
  463 kft=kft+1
      if(dif.lt.trwa) go to 404
      if (dif.lt.trwb) go to 405
      if (dif.lt.trwc) go to 406
      if (dif.lt.trwd) go to 407
      if (dif.lt.trwe) go to 408
      if (dif.lt.trwf) go to 409
      if (dif.lt.trwg) go to 410
      if (dif.lt.trwh) go to 411
      if (dif.lt.trwi) go to 412
      if (dif.lt.trwj) go to 413
      go to 5
  404 istm(1)=istm(1)+kml
      do 415 l=1,nrootx
  415 trsum(1,l)=trsum(1,l)+ea(l)
      iqr(kft)=1
      if (kft.lt.iwod) go to 416
      kft=0
      write (kfile) iqr
      go to 416
  405 istm(2)=istm(2)+kml
      do 417 l=1,nrootx
  417 trsum(2,l)=trsum(2,l)+ea(l)
      iqr(kft)=2
      if (kft.lt.iwod) go to 416
      kft=0
      write (kfile) iqr
      go to 416
  406 istm(3)=istm(3)+kml
      do 418 l=1,nrootx
  418 trsum(3,l)=trsum(3,l)+ea(l)
      iqr(kft)=3
      if (kft.lt.iwod) go to 416
      kft=0
      write (kfile) iqr
      go to 416
  407 istm(4)=istm(4)+kml
      do 419 l=1,nrootx
  419 trsum(4,l)=trsum(4,l)+ea(l)
      iqr(kft)=4
      if (kft.lt.iwod) go to 416
      kft=0
      write (kfile) iqr
      go to 416
  408 istm(5)=istm(5)+kml
      do 420 l=1,nrootx
  420 trsum(5,l)=trsum(5,l)+ea(l)
      iqr(kft)=5
      if (kft.lt.iwod) go to 414
      kft=0
      write (kfile) iqr
      go to 414
  409 istm(6)=istm(6)+kml
      do 421 l=1,nrootx
  421 trsum(6,l)=trsum(6,l)+ea(l)
      iqr(kft)=6
      if (kft.lt.iwod) go to 414
      kft=0
      write (kfile) iqr
      go to 414
  410 istm(7)=istm(7)+kml
      do 422 l=1,nrootx
  422 trsum(7,l)=trsum(7,l)+ea(l)
      iqr(kft)=7
      if (kft.lt.iwod) go to 414
      kft=0
      write (kfile) iqr
      go to 414
  411 istm(8)=istm(8)+kml
      do 423 l=1,nrootx
  423 trsum(8,l)=trsum(8,l)+ea(l)
      iqr(kft)=8
      if (kft.lt.iwod) go to 414
      kft=0
      write (kfile) iqr
      go to 414
  412 istm(9)=istm(9)+kml
      do 424 l=1,nrootx
  424 trsum(9,l)=trsum(9,l)+ea(l)
      iqr(kft)=9
      if (kft.lt.iwod) go to 414
      kft=0
      write (kfile) iqr
      go to 414
  413 istm(10)=istm(10)+kml
      do 425 l=1,nrootx
  425 trsum(10,l)=trsum(10,l)+ea(l)
      iqr(kft)=10
      if (kft.lt.iwod) go to 414
      kft=0
      write (kfile) iqr
      go to 414
   5  iqr(kft)=11
      if (kft.lt.iwod) go to 414
      kft=0
      write (kfile) iqr
414   corea=dif
      if (.not.oprint(31).and.dif.lt.sch) go to 685
      if(oprint(33)) then
      write(jtape,401)jpag,ipag,types,(itag(loop),loop=1,nx)
      write(jtape,4401)(iea(mo),mo=1,nrootx)
401   format(2i7,4x,a1,3x,24a3/22x,24a3)
4401  format((43x,10(1x,i7)) )
      endif
  685 if(ical.eq.2) goto 2
      go to 700
416   if(.not.oprint(32)) go to 573
      write(jtape,123)j,typed,(itag(loop),loop=1,nx)
123   format(i14,4x,a1,3x,24a3/22x,24a3)
      write(jtape,4401)(iea(mo),mo=1,nrootx)
      go to 573
 465  kft=kft+1
      iqr(kft)=12
      if (kft.lt.iwod) go to 381
      kft=0
      write (kfile) iqr
381   corea=crub
      if(oprint(33)) then
      write(jtape,701)jpag,ipag,typer,(itag(loop),loop=1,nx)
      write(jtape,7701)eigv
701   format(2i7,4x,a1,3x,24a3/22x,24a3)
7701  format(43x,f12.8)
      endif
      if(ical.eq.2) goto 2
      go to 700
  1   kft=kft+1
      if (ical.eq.2) go to 2
      iqr(kft)=11
      if (kft.lt.iwod) go to 700
      kft=0
      write (kfile) iqr
      go to 700
  2   if (kft.lt.iwod) go to 686
      kft=0
      read (kfile) iqr
  686 lpix=lpix+nx
      jpag=jpag+1
      ipag=ipag+kml
      if(lpix.lt.iwod) goto 573
      lpix=lpix-iwod
      read(linf)
      goto 573
700   do 704 ll=1,nx
      lpix=lpix+1
      jkan(lpix)=itest(ll)
      if (lpix.lt.iwod) go to 704
      lpix=0
      write(linf)jkan
704   continue
      jpag=jpag+1
      ipag=ipag+kml
      do 705 ll=1,kml
      ifr=ifr+1
      senk(ifr)=corea
      if(ifr.lt.irsf) go to 705
      ifr=0
      write(ideli) senk
705   continue
573   continue
      return
      end
      subroutine sm1310
      implicit REAL (a-h,o-z), integer(i-n)
      integer olab,ot
_IF(cray,ksr,i8)
      integer olab8
_ENDIF
INCLUDE(common/sizes)
      common /b/ iwod,nid,lg,irsf,ntape,ntype,ical,ntab,nc,lt,lz,nx,
     + nd,nston,itape,jtape,mdisk,ideli,mtape,iswh,nr,mm,jblk,jto,
     + igmax,nr1,mtype,linf,kfile,nawk,isec,nshl,
     + nsec,lulu,nrootx,lsng,nprin,ipt0,
     + mconf(5),nconf(5),mh,jswh,
     + id(mxcsf),mn(mxcsf),nsc(mxcsf),
     + imo,mc,mz,ndx,nrx,ntx,ny,nt,iht,ihz,nz,mt,nzx,jy,
     + jx,nd1,nd2,md1,md2,nr3,nr4,mr1,mr2,nt1,nt2,mt1,mt2,ih1,ih2,ih3,
     + ih4,jh2,jh3,jh4,iab,nt3,iac,nr2,m,nad,nad1,ifr1, mx,md,mr,ipag
       common /d/ jbab,kmj,ifl,iiz,nis,iix,nzer,kml,ndt,jsec,iej,ie,
     * ia,kc,kmk,ndk,kr,iiq,ic1,ic,i21,i2,i31,i3,ii,
     * kzer,ndj,iek,j3,j3a,mps,nps,mms,nms,jc,jca,lcl,ndi,kmi,
     * ja,jai,nzei,nisi,lc,ll,ndl,km2,la,lzer,lr,lps,lms,
     * ic2,i22,i32
     *,idum10(2)
      common/linkmr/ icex,i2ex,i3ex
      common/bufd/ gout(510),nword
      common/scrtch/ndet(5),nsac(5),iaz(5),iez(5),idra(5),
     + idrc(5),jdrc(5),jan(7),jbn(7),ispa2(7),
     + ae(7829),ot(44000),ideks(maxorb),
     + olab(2000),ij(8),nytl(5),nplu(5),ndub(5),nod(5),jkan(mxcrec),
     + nir(maxorb),loc(maxorb),sac(10),
_IF(cray,ksr,i8)
     + olab8(250),f(2304),h(mxcrec),isc(5),bob(3),itest(mxnshl),
_ELSE
     + olab8(250),f(2304),h(mxcrec),isc(5),isd,bob(3),itest(mxnshl),
_ENDIF
     + core,
     + w(mxcsf*(mxcsf+1)/2),
     + ew(mxroot),vect(mxcsf*mxroot),trsum(10,mxroot),
     * wect(mxcsf*mxroot),senk(500),iqr(mxcrec),istm(10)
     *,ht(mxroot),hs(mxcsf),ea(mxroot),iper(10),
     * coff(mxcsf*mxroot),doff(mxcsf),nrec(5),jtrk(mxcsf),
     * isym(mxroot+1),mkon(mxcrec),ikan(mxcrec),hy(mxcsf*mxcsf)
     *,ktrk(mxcsf)
      common/jany/ jerk(10),jbnk(10),idum11(10)
      dimension hp(4800)
      equivalence (hp(1),hy(1))
      do 13 l=1,mc
      kbab=jbab
      jbab=jbab+kml
      if (ifl.eq.1) go to 400
      ifl=1
      go to 14
 400  mm=mm+1
      ll=olab(mm)
      if(mm.lt.nid) go to 14
      mm=0
      read(ntape)olab8
      call unpack(olab8,8,olab,2000)
14    if(ll.eq.0) go to 405
      mm=mm+1
      jj=olab(mm)
      if (mm.lt.nid) go to 15
      mm=0
      read(ntape) olab8
      call unpack(olab8,8,olab,2000)
15    jj=jj*ii+nis
      jto=jto+1
      ip=ot(jj)
      if(ip.eq.255) go to 803
      coul=h(jto)
      if (jto.lt.iwod) go to 16
      jto=0
      read(mtype) h
16    jto=jto+1
      exc=-h(jto)
      if (jto.lt.iwod) go to 805
      jto=0
      read(mtype)h
      go to 805
803   coul=-h(jto)
      if (jto.lt.iwod) go to 904
      jto=0
      read(mtype) h
904   jto=jto+1
      exc=h(jto)
      if(jto.lt.iwod) go to 805
      jto=0
      read(mtype) h
 805  if(ll-2) 800,801,802
 800  bob(1)=-coul-exc
      bob(2)=exc
      bob(3)=coul
      go to 804
 801  bob(1)=-exc
      bob(2)=coul+exc
      bob(3)=-coul
      go to 804
802   bob(1)=exc
      bob(2)=coul
      bob(3)=-coul-exc
804   continue
_IF1()c     do 905 m=1,nzer
_IF1()c905  f(m)=0
      call vclr(f,1,nzer)
      do 18 m=1,kml
      jj=jj+1
      my=ot(jj)
      jj=jj+1
      if (my.eq.0) go to 18
_IF1(iv)      kx=ja+my
_IFN1(iv)      kx=ja+my+ndj
      mike=ot(jj)
      tim=bob(mike)
_IF1(iv)       lx=m-kml
_IF1(iv)      do 22 if=1,kmj
_IF1(iv)      kx=kx+ndj
_IF1(iv)      lx=lx+kml
_IF1(iv) 22    f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv)      call daxpy(kmj,tim,ae(kx),ndj,f(m),kml)
18    continue
_IF1()c      ig7=kbab-jsec
_IF1()c      mx=ie
_IF1()c      do 23 m=1,kml
_IF1()c      ly=0
_IF1()c      ig7=ig7+1
_IF1()c      my=mx
_IF1()c      in3=ig7
_IF1()c      do 23 if=1,kmj
_IF1()c      in3=in3+jsec
_IF1()c      tim=0
_IF1()c      mx=my
_IF1()c      do 24 ig=1,kml
_IF1()c      ly=ly+1
_IF1()c      mx=mx+1
_IF1()c24    tim=tim+f(ly)*ae(mx)
_IF1()c  23  hp(in3)=tim
_IF1(c)      call mxma(ae(ie+1),kml,1,f,1,kml,hp(kbab+1),1,jsec,kml,kml,kmj)
_IFN1(c)      call mxmaa(ae(ie+1),kml,1,f,1,kml,hp(kbab+1),1,jsec,kml,kml,kmj)
      go to 13
 405  ig7=kbab-jsec
      do 406 m=1,kml
      ig7=ig7+1
      in3=ig7
      do 406 if=1,kmj
      in3=in3+jsec
406   hp(in3)=0
13    continue
      return
      end
      subroutine sm1312(bkay,intmax)
      implicit REAL (a-h,o-z), integer(i-n)
      integer olab,ot
_IF(cray,ksr,i8)
      integer olab8
_ENDIF
INCLUDE(common/sizes)
      common/linkmr/ icxx,i2xx,i3xx
      common/bufd/ gout(510),nword
      common/scrtch/ndet(5),nsac(5),iaz(5),iez(5),idra(5),
     + idrc(5),jdrc(5),jan(7),jbn(7),ispa2(7),
     + ae(7829),ot(44000),ideks(maxorb),
     + olab(2000),ij(8),nytl(5),nplu(5),ndub(5),nod(5),jkan(mxcrec),
     + nir(maxorb),loc(maxorb),sac(10),
_IF(cray,ksr,i8)
     + olab8(250),f(2304),h(mxcrec),isc(5),bob(3),itest(mxnshl),
_ELSE
     + olab8(250),f(2304),h(mxcrec),isc(5),isd,bob(3),itest(mxnshl),
_ENDIF
     + core,
     + w(mxcsf*(mxcsf+1)/2),
     + ew(mxroot),vect(mxcsf*mxroot),trsum(10,mxroot),
     + wect(mxcsf*mxroot),senk(500),iqr(mxcrec),istm(10)
     +,ht(mxroot),hs(mxcsf),ea(mxroot),iper(10),
     + coff(mxcsf*mxroot),doff(mxcsf),nrec(5),jtrk(mxcsf),
     + isym(mxroot+1),mkon(mxcrec),ikan(mxcrec),hy(mxcsf*mxcsf)
     +,ktrk(mxcsf)
      common/jany/ jerk(10),jbnk(10)
     *,idum11(10)
      dimension bkay(intmax)
      dimension hp(4800)
      equivalence (hp(1),hy(1))
      common /b/ iwod,nid,lg,irsf,ntape,ntype,ical,ntab,nc,lt,lz,nx,
     + md1,nston,itape,jtape,mdisk,ideli,mtape,iswh,mr,mm,jblk,jto,
     + igmax,nr1,mtype,linf,kfile,
     + nawk,isec,nshl,
     + nsec,lulu,nrootx,lsng,nprin,ipt0,mconf(5),nconf(5),mh,jswh,
     + id(mxcsf),mn(mxcsf),nsc(mxcsf),
     + imo,kc,mz,ndx,nrx,
     + ntx,ny,nt,iht,ihz,nz,mt,nzx,jy,
     + jx,nd1,nd2,nd,md2,nr3,nr4,mr1,mr2,nt1,nt2,mt1,mt2,ih1,ih2,ih3,
     + ih4,jh2,jh3,jh4,iab,nt3,iac,nr2,m,nad,nad1,ifr1, mx,md,nr,ipag
       common /d/ jbab,kmj,ifl,iiq,nis,iiz,kzer,kmk,ndt,jsec,iek,iej,
     * ia,mc,kml,ndk,kr,ii,ic,ic1,i2,i21,i3,i31,iix,
     * nzer,ndj,ie,j3a,j3,nps,mps,nms,mms,jca,jc,lcl,ndi,kmi,
     * ja,jai,nzei,nisi,lc,kk,ndl,km2,la,lzer,lr,lps,lms,
     * ic2,i22,i32
     *,idum10(2)
      do 29 l=1,mc
      kbab=jbab
      jbab=jbab+kml
      if (ifl.eq.1) go to 400
      ifl=1
      go to 30
 400  mm=mm+1
      kk=olab(mm)
_IF1()c     if(jqx.gt.7100) write(6,*) 'sm1312: ',mc,kml,jbab,mm,kk,l
_IF1()c     if (jqx.eq.33) write (6,900) mm,kk,l,mc,jbab
      if (mm.lt.nid) go to 30
      mm=0
      read(ntape) olab8
      call unpack(olab8,8,olab,2000)
30    go to (405,31,32,33), kk
 31   mm=mm+1
      ll=olab(mm)
      if (ll.gt.128) ll=ll-256
      if (mm.lt.nid) go to 34
      mm=0
      read(ntape) olab8
      call unpack(olab8,8,olab,2000)
34    mm=mm+1
      jj=olab(mm)
      if (mm.lt.nid) go to 35
      mm=0
      read(ntape) olab8
      call unpack(olab8,8,olab,2000)
35    mm=mm+1
      kk=olab(mm)
      if (mm.lt.nid) go to 36
      mm=0
      read(ntape) olab8
      call unpack(olab8,8,olab,2000)
36    jj=(kk-1)*mr +jj
      jj=jj*ii+ic
      ip=ot(jj)
      if(ip.eq.255) ll=-ll
      jto=jto+1
      if (ll.lt.0) go to 820
      coul = h(jto)
      if (jto.lt.iwod) go to 37
      jto=0
      read(mtype) h
37    jto=jto+1
      exc=-h(jto)
      if (jto.lt.iwod) go to 38
      jto=0
      read(mtype) h
      go to 38
820   coul=-h(jto)
      if (jto.lt.iwod) go to 821
      jto=0
      read(mtype) h
821   jto=jto+1
      ll=-ll
      exc=h(jto)
      if (jto.lt.iwod) go to 38
      jto=0
      read(mtype) h
38    if (ll-2)810,811,812
810   bob(1)=-exc
      bob(2)=-coul
      bob(3)=coul+exc
      go to 49
811   bob(1)=exc
      bob(2)=-coul-exc
      bob(3)=coul
      go to 49
812   bob(1)=coul+exc
      bob(2)=-exc
      bob(3)=-coul
49    continue
_IF1()c     do 39 m=1,nzer
_IF1()c39   f(m)=0
      call vclr(f,1,nzer)
      do 40 m=1,kml
      jj=jj+1
      my=ot(jj)
      jj=jj+1
      if(my.eq.0) go to 40
_IF1(iv)      kx=ja+my
_IFN1(iv)      kx=ja+my+ndj
      mike=ot(jj)
      if(mike.gt.128) go to 42
      tim=bob(mike)
      go to 43
42    tim=-bob(256-mike)
43    continue
_IF1(iv)       lx=m-kml
_IF1(iv)      do 44 if=1,kmj
_IF1(iv)      kx=kx+ndj
_IF1(iv)      lx=lx+kml
_IF1(iv) 44   f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv)      call daxpy(kmj,tim,ae(kx),ndj,f(m),kml)
40    continue
47    ig7=kbab-jsec
      mx=ie
_IF1()c      do 45 m=1,kml
_IF1()c      ly=0
_IF1()c      ig7=ig7+1
_IF1()c      in3=ig7
_IF1()c      my=mx
_IF1()c      do 45 if=1,kmj
_IF1()c      in3=in3+jsec
_IF1()c      tim=0
_IF1()c      mx=my
_IF1()c      do 46 ig=1,kml
_IF1()c      ly=ly+1
_IF1()c      mx=mx+1
_IF1()c 46   tim=tim+f(ly)*ae(mx)
_IF1()c  45  hp(in3)=tim
_IF1(c)      call mxma(ae(ie+1),kml,1,f,1,kml,hp(kbab+1),1,jsec,kml,kml,kmj)
_IFN1(c)      call mxmaa(ae(ie+1),kml,1,f,1,kml,hp(kbab+1),1,jsec,kml,kml,kmj)
      go to 29
 405  ig7=kbab-jsec
      do 406 m=1,kml
      ig7=ig7+1
      in3=ig7
      do 406 if=1,kmj
      in3=in3+jsec
 406  hp(in3)=0
      go to 29
  32  mm=mm+1
      ll=olab(mm)
      if (mm.lt.nid) go to 51
      mm=0
      read(ntape)olab8
      call unpack(olab8,8,olab,2000)
51    jto=jto+1
      sm=h(jto)
      if(jto.lt.iwod) go to 52
      jto=0
      read(mtype) h
52    if(ll.lt.128) go to 53
      sm=-sm
      ll=256-ll
53    jj=ll*kml+i2
      do 55 m=1,nzer
55    f(m)=0
      do 56 m=1,kml
      jj=jj+1
      my=ot(jj)
      if (my.eq.0) go to 56
      if(my.lt.128) go to 58
      kx=ja-my+256
      tim=-sm
      go to 59
58    tim=sm
      kx=my+ja
59    lx=m-kml
      do 60 if=1,kmj
      kx=kx+ndj
      lx=lx+kml
60    f(lx)=f(lx)+tim*ae(kx)
56    continue
      go to 47
  33  mm=mm+1
      ll=olab(mm)
      if (mm.lt.nid) go to 61
      mm=0
      read(ntape)olab8
      call unpack(olab8,8,olab,2000)
61    mm=mm+1
      ntar=olab(mm)
      if (mm.lt.nid) go to 62
      mm=0
      read (ntape) olab8
      call unpack(olab8,8,olab,2000)
62    mm=mm+1
      nb=olab(mm)
      if(mm.lt.nid) go to 63
      mm=0
      read (ntape) olab8
      call unpack(olab8,8,olab,2000)
63    mm=mm+1
      mb=olab(mm)
      if (mm.lt.nid) go to 41
      mm=0
      read (ntape) olab8
      call unpack(olab8,8,olab,2000)
  41  nb=ideks(nb)+mb+ij(ntar)
      sm=bkay(nb)
      jto=jto+1
      sm=sm+h(jto)
      if(jto.lt.iwod) go to 64
      jto=0
      read(mtype) h
64    if (mr.eq.0 ) go to 65
      do 66 m=1,mr
      jto=jto+1
      sm=sm+h(jto)
      if(jto.lt.iwod) go to 68
      jto=0
      read(mtype) h
68    jto=jto+1
      sac(m)=-h(jto)
      if (jto.lt.iwod) go to 66
      jto=0
      read(mtype) h
66    continue
65    if(nd.eq.0) go to 67
      do 69 m=1,nd
      jto=jto+1
      tim=h(jto)
      sm=sm+tim+tim
      if(jto.lt.iwod) go to 70
      jto=0
      read(mtype) h
70    jto=jto+1
      sm=sm-h(jto)
      if(jto.lt.iwod) go to 69
      jto=0
      read(mtype) h
69    continue
67    if(ll.gt.128) go to 86
      jj=ll*j3+i3
      ip=ot(jj)
      if(ip.eq.1) go to 71
      sm=-sm
      if(mr.eq.0) go to 71
      do 72 m=1,mr
72    sac(m)=-sac(m)
71    jj=jj+1
      do 80 m=1,nzer
80    f(m)=0
      jx=-kml
      do 73 m=1,kml
      jx=jx+1
      in=jj
      im=ot(in)
      go to (74,75,76,77),im
 74   in=in+1
      kx=ot(in)
      kx=kx+ja
      tim=sm
      if(mps.eq.0) go to 78
      do 79 if=1,mps
      in=in+1
      im=ot(in)
79    tim=tim+sac(im)
78    lx=jx
      do 81 if=1,kmj
      kx=kx+ndj
      lx=lx+kml
81    f(lx)=f(lx)+tim*ae(kx)
      go to 73
 75   in=in+1
      kx=ot(in)
      kx=kx+ja
      tim=-sm
      if(mms.eq.0) go to 78
      do 83 if=1,mms
      in=in+1
      im=ot(in)
83    tim=tim-sac(im)
      go to 78
76    do 84 if=1,nms
      in=in+1
      kx=ot(in)
      kx=kx+ja
      in=in+1
      im=ot(in)
      tim=-sac(im)
      lx=jx
      do 84 ig=1,kmj
      kx=kx+ndj
      lx=lx+kml
84    f(lx)=f(lx)+tim*ae(kx)
      go to 73
77    do 85 if=1,nps
      in=in+1
      kx=ot(in)
      kx=kx+ja
      in=in+1
      im=ot(in)
      tim=sac(im)
      lx=jx
      do 85 ig=1,kmj
      kx=kx+ndj
      lx=lx+kml
85    f(lx)=f(lx)+tim*ae(kx)
73    jj=jj+jc
      go to 47
86    jj=(256-ll)*j3+i3
      ip=ot(jj)
      if (ip.eq.1) go to 87
      sm=-sm
      if (mr.eq.0) go to 87
      do 88 m=1,mr
88    sac(m)=-sac(m)
87    jj=jj+1
      do 89 m=1,nzer
89    f(m)=0
      jx=-kml
      do 90 m=1,kml
      jx=jx+1
      in=jj
      im=ot(in)
      go to (91,95,96,97),im
91    in=in+1
      kx=ot(in)
      kx=kx+ja
      tim=sm
      if(mms.eq.0) go to 92
      in=in+mps
      do 93 if=1,mms
      in=in+1
      im=ot(in)
93    tim=tim+sac(im)
92    lx=jx
      do 94 if=1,kmj
      kx=kx+ndj
      lx=lx+kml
94    f(lx)=f(lx)+tim*ae(kx)
      go to 90
95    in=in+1
      kx=ot(in)
      kx=kx+ja
      tim=-sm
      if(mps.eq.0) go to 92
      in=in+mms
      do 99 if=1,mps
      in=in+1
      im=ot(in)
99    tim=tim-sac(im)
      go to 92
96    do 100 if=1,nms
      in=in+1
      kx=ot(in)
      kx=kx+ja
      in=in+1
      im=ot(in)
      tim=sac(im)
      lx=jx
      do 100 ig=1,kmj
      kx=kx+ndj
      lx=lx+kml
100   f(lx)=f(lx)+tim*ae(kx)
      go to 90
97    do 101 if=1,nps
      in=in+1
      kx=ot(in)
      kx=kx+ja
      in=in+1
      im=ot(in)
      tim=-sac(im)
      lx=jx
      do 101 ig=1,kmj
      kx=kx+ndj
      lx=lx+kml
101   f(lx)=f(lx)+tim*ae(kx)
90    jj=jj+jc
      go to 47
29    continue
      return
      end
      subroutine sm1311
      implicit REAL (a-h,o-z), integer(i-n)
      integer olab,ot
_IF(cray,ksr,i8)
      integer olab8
_ENDIF
INCLUDE(common/sizes)
      common/linkmr/ icex,i2ex,i3ex
      common/bufd/ gout(510),nword
      common/scrtch/ndet(5),nsac(5),iaz(5),iez(5),idra(5),
     + idrc(5),jdrc(5),jan(7),jbn(7),ispa2(7),
     + ae(7829),ot(44000),ideks(maxorb),
     + olab(2000),ij(8),nytl(5),nplu(5),ndub(5),nod(5),jkan(mxcrec),
     + nir(maxorb),loc(maxorb),sac(10),
_IF(cray,ksr,i8)
     + olab8(250),f(2304),h(mxcrec),isc(5),bob(3),itest(mxnshl),
_ELSE
     + olab8(250),f(2304),h(mxcrec),isc(5),isd,bob(3),itest(mxnshl),
_ENDIF
     + core,
     + w(mxcsf*(mxcsf+1)/2),
     + ew(mxroot),vect(mxcsf*mxroot),trsum(10,mxroot),
     + wect(mxcsf*mxroot),senk(500),iqr(mxcrec),istm(10)
     +,ht(mxroot),hs(mxcsf),ea(mxroot),iper(10),
     + coff(mxcsf*mxroot),doff(mxcsf),nrec(5),jtrk(mxcsf),
     + isym(mxroot+1),mkon(mxcrec),ikan(mxcrec),hy(mxcsf*mxcsf)
     +,ktrk(mxcsf)
      common/jany/ jerk(10),jbnk(10),idum11(10)
      dimension hp(4800)
      equivalence (hp(1),hy(1))
      common /b/ iwod,nid,lg,irsf,ntape,ntype,ical,ntab,nc,lt,lz,nx,
     + nd,nston,itape,jtape,mdisk,ideli,mtape,iswh,nr,mm,jblk,jto,
     + igmax,nr1,mtype,linf,kfile,nawk,isec,nshl,
     + nsec,lulu,nrootx,lsng,nprin,ipt0,
     + mconf(5),nconf(5),mh,jswh,
     + id(mxcsf),mn(mxcsf),nsc(mxcsf),
     + imo,lcl,mz,ndx,nrx,ntx,
     + ny,nt,iht,ihz,nz,mt,nzx,jy,
     + jx,nd1,nd2,md1,md2,nr3,nr4,mr1,mr2,nt1,nt2,mt1,mt2,ih1,ih2,ih3,
     + ih4,jh2,jh3,jh4,iab,nt3,iac,nr2,m,nad,nad1,ifr1, mx,md,mr
     +,idum1(1)
       common /d/ jbab,kml,ifl,iiz,nisi,ii,nzei,kmi,ndi,jsec,ie,iej,
     * jai,kc,kmk,ndk,kr,iiq,ic1,ic,i21,i2,i31,i3,iix,
     * kzer,ndt,iek,j3,j3a,mps,nps,mms,nms,jc,jca,mc,ndj,kmj,
     * ia,ja,nzer,nis,lc,ll,ndl,km2,la,lzer,lr,lps,lms,
     * ic2,i22,i32
     *,idum10(2)
      do 13 l=1,mc
      kbab=jbab
      jbab=jbab+kmj
      if (ifl.eq.1) go to 400
      ifl=1
      go to 14
400   mm=mm+1
      ll=olab(mm)
      if(mm.lt.nid) go to 14
      mm=0
      read(ntape)olab8
      call unpack(olab8,8,olab,2000)
14    if(ll.eq.0) go to 405
      mm=mm+1
      jj=olab(mm)
      if (mm.lt.nid) go to 15
      mm=0
      read(ntape) olab8
      call unpack(olab8,8,olab,2000)
15    jj=jj*ii+nis
      jto=jto+1
      ip=ot(jj)
      if(ip.eq.255) go to 803
      coul=h(jto)
      if (jto.lt.iwod) go to 16
      jto=0
      read(mtype) h
16    jto=jto+1
      exc=-h(jto)
      if (jto.lt.iwod) go to 805
      jto=0
      read(mtype)h
      go to 805
803   coul=-h(jto)
      if (jto.lt.iwod) go to 904
      jto=0
      read(mtype) h
904   jto=jto+1
      exc=h(jto)
      if(jto.lt.iwod) go to 805
      jto=0
      read(mtype) h
 805  if(ll-2) 800,801,802
 800  bob(1)=-coul-exc
      bob(2)=exc
      bob(3)=coul
      go to 804
 801  bob(1)=-exc
      bob(2)=coul+exc
      bob(3)=-coul
      go to 804
802   bob(1)=exc
      bob(2)=coul
      bob(3)=-coul-exc
804   continue
_IF1()c     do 905 m=1,nzer
_IF1()c905  f(m)=0
      call vclr(f,1,nzer)
      do 18 m=1,kml
      jj=jj+1
      my=ot(jj)
      jj=jj+1
      if (my.eq.0) go to 18
_IF1(iv)      kx=ja+my
_IFN1(iv)      kx=ja+my+ndj
      mike=ot(jj)
      tim=bob(mike)
_IF1(iv)      lx=m-kml
_IF1(iv)      do 22 if=1,kmj
_IF1(iv)      kx=kx+ndj
_IF1(iv)      lx=lx+kml
_IF1(iv) 22   f(lx)=f(lx)+tim*ae(kx)
_IFN1(iv)      call daxpy(kmj,tim,ae(kx),ndj,f(m),kml)
18    continue
      mx=ie
_IF1(c)      call mxma(ae(mx+1),kml,1,f,1,kml,hp(kbab+1),jsec,1,kml,kml,kmj)
_IFN1(c)      call mxmaa(ae(mx+1),kml,1,f,1,kml,hp(kbab+1),jsec,1,kml,kml,kmj)
      go to 13
 405  ig7=kbab-jsec
      do 406 m=1,kml
      ig7=ig7+jsec
      in3=ig7
      do 406 if=1,kmj
      in3=in3+1
406   hp(in3)=0
13    continue
      return
      end
      subroutine smrd11(iot,ideks)
      implicit REAL (a-h,o-z), integer(i-n)
      integer olab
_IF(cray,ksr,i8)
       integer olab8
_ENDIF
INCLUDE(common/sizes)
INCLUDE(common/prints)
      common /b/ iwod,nid,lg,irsf,ntape,ntype,ical,ntab,nc,lt,lz,nx,
     + nd,nston,itape,jtape,ltape,ideli,mtape,iswh,nr,mm,jblk,jto,
     + igmax,nr1,ltype,linf,kfile,nawk,isec,nshl,
     + nsec,lulu,nrootx,lsng,nprin,ipt0,nconf(5),mconf(5),mh,jswh,
     + nn(mxcsf),mn(mxcsf),nsc(mxcsf),
     + imo,mc,mz,ndx,nrx,ntx,ny,nt,iht,ihz,nz,mt,nzx,jy,
     + jx,nd1,nd2,md1,md2,nr3,nr4,mr1,mr2,nt1,nt2,mt1,mt2,ih1,ih2,ih3,
     + ih4,jh2,jh3,jh4,iab,nt3,iac,nr2,m,nad,nad1,ifr1,mx,md,mr,j
      common/scrtch/jdeks(9),kdeks(8),
     + lj(8),olab(2000),nit(666),nytl(5),ndub(5),
     + icon(5),nod(5),jdum(maxorb),nir(maxorb),loc(maxorb)
      common/junk/ican(3),icbn,olab8(250),jcon(8200),mkon(mxcrec),
     + kc(mxcrec),kd(mxcrec),lab(3)
     +,ij(8),jtest(mxnshl),nyzl(mxcsf),jkan(mxcrec),nplu(5),
     + jkon(mxnshl*mxcsf)
      dimension iot(*),ideks(*)
      call rewftn(mtape)
      read (mtape)
      ih=imo*(lulu+1)
      do 500 i=1,ih
  500 jcon(i)=0
      ny=0
      nx=0
      do 501 i=1,jswh
      nc=mconf(i)
      if (nc.eq.0) go to 501
      nr=nod(i)
      nd=ndub(i)
      do 502 j=1,nc
      nx=nx+imo
      if (nr.eq.0) go to 504
      do 503 k=1,nr
      ny=ny+1
      jd=iot(ny)+nx
      if(jd.gt.8200)write(6,5002)jd
 5002 format(1x,' smrd115002:jd= ',i8)
 503  jcon(jd)=1
      if (nd.eq.0) go to 502
 504  do 505 k=1,nd
      ny=ny+1
      jd=iot(ny)+nx
      if(jd.gt.8200)write(6,5003)jd
 5003 format(1x,' smrd115003:jd= ',i8)
 505  jcon(jd)=2
 502  continue
 501  continue
      ky=imo
      do 515 i=1,jswh
      ican(i)=ky
 515  ky=ky+imo*mconf(i)
      nt1=ny+1
      lulb=1
      ksc=nsc(1)
      jnerd=nn(1)
      do 506 i=1,iswh
      nc=nconf(i)
      if (nc.eq.0) go to 506
      read (mtape) mkon
      ig=0
      ndx=ndub(i)
      nrx=nod(i)
      ntx=nytl(i)
      iht=icon(i)
      ihz=iht+nrx
      nzx=ny+nrx
      jy=ican(i+2)
      jx=ican(i+1)
      nd1=ndx+1
      nd2=ndx+2
      md1=ndx-1
      md2=ndx-2
      nr3=nrx-2
      nr4=nrx-4
      mr1=nrx+2
      mr2=nrx+4
      nt3=ntx-1
      nt2=ntx-2
      mt1=ntx+1
      mt2=ntx+2
      ih1=icon(i+1)
      ih2=icon(i+2)
      ih3=icon(i-1)
      ih4=icon(i-2)
      jh2=ih2+mr2
      jh3=ih3+nr3
      jh4=ih4+nr4
      iab=ideks(nrx+2)
      nr1=nrx+1
      iac=ideks(nrx)
      nr2=nrx-1
      do 507 j=1,nc
      if (j.eq.0) write (6,8009) (mkon(iqq),iqq=1,10)
 8009 format(5x,10i10)
      ky=ny
      do 508 k=1,ntx
      ig=ig+1
      ky=ky+1
      if (j.eq.0) write  (6,8003) ig,ky,iwod
 8003 format(5x,3i7,'i am in smrd11 for ig,ky')
 8010 format(5x,3i8,'iot(ky),mkon(ig),mkon(1)')
      iot(ky)=mkon(ig)
       if (j.eq.0) write (6,8010) iot(ky),mkon(ig),mkon(1),ky,ig,iwod
      if (ig.lt.iwod) go to 508
      ig=0
      read (mtape) mkon
 508  continue
_IF1()c     kkq=ny+1
_IF1()c     kkz=ny+ntx
_IF1()c     if(j.eq.7110) write(6,*) 'iot at ',j
_IF1()c     if(j.eq.7110) write(6,1860) (iot(jjjjj),jjjjj=ny+1,ny+ntx)
_IF1()c1860 format(1x,14i5)
_IF1()c     if(j.eq.7111.or.j.eq.7109) write(6,*) 'mm,olab at ',j,mm
_IF1()c     if(j.eq.7111.or.j.eq.7109) write(6,1860) olab
      if (ksc.ne.i.or.jnerd.ne.j) go to 525
      mm=mm+1
      olab(mm)=6
      if (mm.lt.nid) go to 527
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
 527  if (lulb.eq.lulu) go to 528
      lulb=lulb+1
      ksc=nsc(lulb)
      jnerd=nn(lulb)
      go to 507
 528  ksc=0
      go to 507
 525  if (lsng.eq.0) go to 529
      nix=0
      ky=ny
      if (nrx.eq.0) go to 530
      do 531 k=1,nrx
      ky=ky+1
      jd=iot(ky)+mh
      if(jd.gt.8200)write(6,5010)jd
 5010 format(1x,' smrd115010: jd= ',i8)
      if (jcon(jd).gt.0) go to 531
      if (nix.eq.1) go to 529
      nix=1
 531  continue
      if (ndx.eq.0) go to 532
 530  do 533 k=1,ndx
      ky=ky+1
      jd=iot(ky)+mh
      if(jd.gt.8200)write(6,5007)jd
 5007 format(1x,' smrd115007:jd= ',i8)
      if (jcon(jd)-1) 529,534,533
 534  if (nix.eq.1) go to 529
      nix=1
 533  continue
 532  mm=mm+1
      olab(mm)=7
      if (mm.lt.nid) go to 507
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      go to 507
 529  ky=ny
      if (nrx.eq.0) go to 509
      do 510 k=1,nrx
      ky=ky+1
      jd=iot(ky)
      if(jd.gt.8200)write(6,5008)jd
 5008 format(1x,' smrd115008: jd= ',i8)
 510  jcon(jd)=1
      if (ndx.eq.0) go to 511
 509  do 512 k=1,ndx
      ky=ky+1
      jd=iot(ky)
      if(jd.gt.8200)write(6,5009)jd
 5009 format(1x,' smrd115009: ***error*** jd= ',i8)
 512  jcon(jd)=2
 511  do 513 k=1,jswh
      mc=mconf(k)
      if (mc.eq.0) go to 513
      kfc=i-k+3
      go to (520,521,522,523,524),kfc
      go to 513
 522  call smr110(iot,ideks)
      go to 513
 520  call smr111(iot,ideks)
      go to 513
 524  call smr112(iot,ideks)
      go to 513
 521  call smr113(iot,ideks)
      go to 513
 523  call smr114(iot,ideks)
 513  continue
      ky=ny
      do 535 l=1,ntx
      ky=ky+1
      ll=iot(ky)
 535  jcon(ll)=0
 507  continue
 506  continue
      if(.not.oprint(29))write(jtape,431)
431   format(/2x,' end of label generation'/)
      return
      end
      subroutine smr111(iot,ideks)
      implicit REAL (a-h,o-z), integer(i-n)
      integer olab
_IF(cray,ksr,i8)
      integer olab8
_ENDIF
INCLUDE(common/sizes)
      common /b/ iwod,nid,lg,irsf,ntape,ntype,ical,ntab,nc,lt,lz,mt2,
     + md2,nston,itape,jtape,mdisk,ideli,mtape,iswh,mr2,mm,jblk,jto,
     + igmax,nr1,ltype,linf,kfile,nawk,isec,nshl,
     + nsec,lulu,nrootx,lsng,nprin,ipt0,mconf(5),nconf(5),mh,jswh,
     + id(mxcsf),ie(mxcsf),nsc(mxcsf),
     + imo,mc,nzx,md,mr,mx,mt,ih2,iht,ihz,jh2,ny,mz,jy,
     + jx,nd1,nd2,md1,nd,nr3,nr4,mr1,nr,nt1,nt2,mt1,nx,ih1,nt,ih3,
     + ih4,nz,jh3,jh4,iab,nt3,iac,nr2,m,nad,nad1,ifr1,ntx,ndx,nrx
     +,idum1(1)
      common/scrtch/jdeks(9),kdeks(8),
     + lj(8),olab(2000),nit(666),nytl(5),ndub(5),
     + icon(5),nod(5),jdum(maxorb),nir(maxorb),loc(maxorb)
      common/junk/ican(3),icbn,olab8(250),jcon(8200),mkon(mxcrec),
     + kc(mxcrec),kd(mxcrec),lab(3),ij(8),jtest(mxnshl),
     + nyzl(mxcsf),jkan(mxcrec),nplu(5),jkon(mxnshl*mxcsf)
      dimension iot(*),ideks(*)
      l2=jy
_IF1()      l3=nz
      l4=nt
      do 13 j=1,mc
      jz=mz
      do 19 l=1,md
      jz=jz+1
      ll=iot(jz)
      if (jcon(ll+l2)-1) 20,21,19
   21 ia=l+1
      go to 22
   19 continue
   20 mm=mm+1
      olab(mm)=0
      if (mm.lt.nid) go to 61
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      mm=0
      go to 61
   22 do 23 l=ia,md
      jz=jz+1
      nn=iot(jz)
      if (jcon(nn+l2)-1) 20,24,23
   24 ja=l+1
      go to 25
   23 continue
   25 if (ja.gt.md) go to 26
      do 27 l=ja,md
      jz=jz+1
      ii=iot(jz)
      if (jcon(ii+l2).ne.2) go to 20
   27 continue
   26 if (mr.eq.0) go to 28
      jz=mt
      do 29 l=1,mr
      jz=jz+1
      ii=iot(jz)
      if (jcon(ii+l2).ne.1) go to 20
   29 continue
   28 jz=l4
      kz=mt+1
      jq=0
      if (mr.gt.0) iq=iot(kz)
      do 30 l=1,nr
      jz=jz+1
      ii=iot(jz)
      if (ii.eq.ll) go to 31
      if (jq.eq.mr.or.iq.ne.ii) go to 32
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   30 continue
   31 ip=l
      ia=l+1
      do 33 l=ia,nr
      jz=jz+1
      ii=iot(jz)
      if (ii.eq.nn) go to 34
      if (jq.eq.mr.or.iq.ne.ii) go to 35
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   33 continue
   32 ip=l
      ia=l+1
      do 36 l=ia,nr
      jz=jz+1
      jj=iot(jz)
      if (jj.eq.ll) go to 37
      if (jq.eq.mr.or.iq.ne.jj) go to 38
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   36 continue
   34 jp=ideks(l-1)
      ia=l+1
      do 39 l=ia,nr
      jz=jz+1
      ii=iot(jz)
      if (jq.eq.mr.or.iq.ne.ii) go to 40
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   39 continue
   40 kp=jdeks(l-2)
      ia=l+1
      do 41 l=ia,nr
      jz=jz+1
      jj=iot(jz)
      if (jq.eq.mr.or.iq.ne.jj) go to 42
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   41 continue
   42 mm=mm+1
      olab(mm)=1
      lp=kdeks(l-3)
      if (mm.lt.nid) go to 43
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      mm=0
      go to 43
   35 jp=ideks(l-1)
      ia=l+1
      do 44l=ia,nr
      jz=jz+1
      jj=iot(jz)
      if (jj.eq.nn) go to 45
      if (jq.eq.mr.or.jj.ne.iq) go to 46
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   44 continue
   37 jp=ideks(l-1)
      ia=l+1
      do 47 l=ia,nr
      jz=jz+1
      jj=iot(jz)
      if (jj.eq.nn) go to 48
      if (jq.eq.mr.or.jj.ne.iq) go to 49
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   47 continue
   38 jp=ideks(l-1)
      ia=l+1
      do 50 l=ia,nr
      jz=jz+1
      kk=iot(jz)
      if (kk.eq.ll) go to 51
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   50 continue
   51 kp=jdeks(l-2)
      ia=l+1
      do 52 l=ia,nr
      jz=jz+1
      kk=iot(jz)
      if (kk.eq.nn) go to 42
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   52 continue
   48 kp=jdeks(l-2)
      ia=l+1
      do 53 l=ia,nr
      jz=jz+1
      jj=iot(jz)
      if (jq.eq.mr.or.jj.ne.iq) go to 54
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   53 continue
   54 lp=kdeks(l-3)
      mm=mm+1
      olab(mm)=3
      if (mm.lt.nid) go to 43
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      mm=0
      go to 43
   49 kp=jdeks(l-2)
      ia=l+1
      do 55 l=ia,nr
      jz=jz+1
      kk=iot(jz)
      if (kk.eq.nn) go to 56
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   55 continue
   56 lp=kdeks(l-3)
      mm=mm+1
      olab(mm)=2
      if (mm.lt.nid) go to 43
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      mm=0
      go to 43
   45 kp=jdeks(l-2)
      ia=l+1
      do 57 l=ia,nr
      jz=jz+1
      jj=iot(jz)
      if (jq.eq.mr.or.jj.ne.iq) go to 56
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   57 continue
   46 kp=jdeks(l-2)
      ia=l+1
      do 58 l=ia,nr
      jz=jz+1
      kk=iot(jz)
      if (kk.eq.nn) go to 54
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   58 continue
 43   mm=mm+1
      olab(mm)=ip+jp+kp+lp
      if (mm.lt.nid) go to 170
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      mm=0
 170  nt1r = nir(ll)
      nl1r = loc (ll)
      nt2r = nir (nn)
      nl2r = loc (nn)
      kix=3
      nb1r = nir(ii)
      nm1r = loc (ii)
      nb2r = nir (jj)
      nm2r = loc (jj)
  120 if (nt1r-nb1r)123,122,199
  199 iax = ideks(nt1r)+nb1r
      if (nt2r.lt.nb2r) go to 124
      ibx = ideks(nt2r)+nb2r
      if (iax-ibx) 125,126,127
  127 icx = ideks(iax) + ibx
      icq = (nl1r-1)*lj(nb1r) + nm1r
  140 ljn = lj(nb2r)
      idq = (nl2r-1)*ljn + nm2r
      icq = (icq-1)*ljn*lj(nt2r)+idq
      go to 129
  126 icx = ideks (iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn+nm1r
  141 idq = (nl2r-1)*ljn+nm2r
      go to 133
  125 icx = ideks(ibx) + iax
      icq = (nl2r-1)*lj(nb2r)+nm2r
  135 ljn = lj(nb1r)
      idq = (nl1r-1)*ljn+nm1r
      icq =(icq-1)*ljn*lj(nt1r) +idq
      go to 129
  124 ibx = ideks(nb2r) + nt2r
      if (iax-ibx) 130,131,132
  132 icx = ideks(iax)+ibx
      icq = (nl1r-1)*lj(nb1r)+nm1r
      ljn = lj(nt2r)
      idq = (nm2r-1)*ljn+nl2r
      icq = (icq-1)*ljn*lj(nb2r) + idq
      go to 129
  131 icx = ideks(iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn+nm1r
      idq = (nm2r-1)*ljn+nl2r
  133 icq = min(icq,idq) + ideks( max(icq,idq))
      go to 129
  130 icx=ideks(ibx)+iax
      icq = (nm2r-1)*lj(nt2r)+nl2r
      go to 135
  123 iax = ideks(nb1r)+nt1r
      if (nt2r.lt.nb2r) go to 136
      ibx = ideks(nt2r)+nb2r
      if (iax-ibx) 137,138,139
  139 icx = ideks(iax)+ibx
      icq = (nm1r-1)*lj(nt1r)+nl1r
      go to 140
  138 icx = ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      go to 141
  137 icx = ideks(ibx) + iax
      icq = (nl2r-1)*lj(nb2r) +nm2r
  145 ljn = lj(nt1r)
      idq = (nm1r-1)*ljn         +nl1r
      icq = (icq-1)*ljn*lj(nb1r) +idq
      go to 129
  136 ibx = ideks(nb2r)+nt2r
      if (iax-ibx) 142,143,144
  144 icx = ideks(iax)+ibx
      icq =(nm1r-1) *lj(nt1r)+nl1r
      ljn = lj(nt2r)
      idq = (nm2r-1)*ljn+nl2r
      icq = (icq-1)*ljn*lj(nb2r)+idq
      go to 129
  143 icx=ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      idq = (nm2r-1)*ljn+nl2r
      go to 133
  142 icx = ideks(ibx)+iax
      icq =(nm2r-1)*lj(nt2r)+nl2r
      go to 145
  122 iax=ideks(nt1r+1)
      iay = min(nl1r,nm1r) + ideks( max(nl1r,nm1r))
      iby = min(nl2r,nm2r) + ideks( max(nl2r,nm2r))
      if (nt1r.eq.nt2r) go to 150
      ibx=ideks(nt2r+1)
      if (iax.lt.ibx) go to 151
      icx=ideks(iax)+ibx
      ljn =lj(nt2r)
      icq = (iay-1)*ideks(ljn+1)+iby
      go to 129
  151 icx=ideks(ibx)+iax
      ljn=lj(nt1r)
      icq =(iby-1)*ideks(ljn+1)+iay
      go to 129
  150 icx=ideks(iax+1)
      icq = min(iay,iby) + ideks( max(iay,iby))
 129  icq=nit(icx)  +  icq
      idud=(icq -1)/igmax
      jto=jto + 1
      kc(jto) = idud  + 1
      icq=icq  -  idud*igmax
      kd(jto)  = icq
      if(jto.lt.iwod)  go to 160
      write(ntype)  kc,kd
      jto=0
      jblk=jblk  + 1
  160 if (kix.lt.0) go to 61
      kix=-3
      itr=nb1r
      nb1r=nb2r
      nb2r=itr
      itr=nm1r
      nm1r=nm2r
      nm2r=itr
      go to 120
   61 l4=l4+nx
_IF1()      l3=l3+nx
  13  l2=l2+imo
      return
      end
      subroutine smr113(iot,ideks)
      implicit REAL (a-h,o-z), integer(i-n)
      integer olab
_IF(cray,ksr,i8)
      integer olab8
_ENDIF
INCLUDE(common/sizes)
      common /b/ iwod,nid,lg,irsf,ntape,ntype,ical,ntab,nc,lt,lz,mt1,
     + md1,nston,itape,jtape,mdisk,ideli,mtape,iswh,mr1,mm,jblk,jto,
     + igmax,nr1,ltype,linf,kfile,
     + nawk,isec,nshl,
     + nsec,lulu,nrootx,lsng,nprin,ipt0,mconf(5),nconf(5),mh,jswh,
     + id(mxcsf),ie(mxcsf),nsc(mxcsf),
     + imo,mc,nzx,md,mr,mx,mt,ih1,iht,ihz,nz,ny,mz,jx,
     + jy,nd1,nd2,nd,md2,nr3,nr4,nr,mr2,nt1,nt2,nx,mt2,nt,ih2,ih3,
     + ih4,jh2,jh3,jh4,iab,nt3,iac,nr2,m,nad,nad1,ifr1,ntx,ndx,nrx
     +,idum1(1)
      common/scrtch/jdeks(9),kdeks(8),
     + lj(8),olab(2000),nit(666),nytl(5),ndub(5),
     + icon(5),nod(5),jdum(maxorb),nir(maxorb),loc(maxorb)
      common/junk/ican(3),icbn,olab8(250),jcon(8200),mkon(mxcrec),
     + kc(mxcrec),kd(mxcrec),lab(3)
     + ,ij(8),jtest(mxnshl),nyzl(mxcsf),jkan(mxcrec),nplu(5),
     + jkon(mxnshl*mxcsf)
      dimension iot(*),ideks(*)
      l3=jy
      l4=nt
      do 65 j=1,mc
      nt4=l4+1
      nz=l4+nr
      jz=mz
      do 70 l=1,md
      jz=jz+1
      jj=iot(jz)
      if (jcon(jj+l3)-1) 71,72,70
   70 continue
   71 if (l.eq.md) go to 73
      ia=l+1
      do 74 l=ia,md
      jz=jz+1
      ll=iot(jz)
      if (jcon(ll+l3).ne.2) go to 75
   74 continue
      go to 73
  75  mm=mm+1
      olab(mm)=1
      if (mm.lt.nid) go to 260
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      mm=0
      go to 260
   72 l1=l
      if (l.eq.md) go to 76
      ia=l+1
      do 78 l=ia,md
      jz=jz+1
      ll=iot(jz)
      if (jcon(ll+l3)-1) 75,77,78
  78  continue
      go to 76
  77  l2=l
      if (l.eq.md) go to 79
      ia=l+1
      do 80 l=ia,md
      jz=jz+1
      kk=iot(jz)
      if (jcon(kk+l3).ne.2) go to 75
   80 continue
      go to 79
   73 if (mr.eq.0) go to 81
      jz=mt
      do 82 l=1,mr
      jz=jz+1
      ll=iot(jz)
      if (jcon(ll+l3).ne.1) go to 75
   82 continue
  81  mm=mm+1
      olab(mm)=3
      if (mm.lt.nid) go to 83
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      mm=0
  83  jz=l4
      kz=mt+1
      jq=0
      if (mr.gt.0) iq=iot(kz)
      do 84 l=1,nr
      jz=jz+1
      ll=iot(jz)
      if (jq.eq.mr.or.iq.ne.ll) go to 85
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   84 continue
   85 ia=l+1
      ip=l
      do 86 l=ia,nr
      jz=jz+1
      nn=iot(jz)
      if (jq.eq.mr.or.iq.ne.nn) go to 87
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   86 continue
  87  mm=mm+1
      olab(mm)=ideks(l-1)+ip
      if (mm.lt.nid) go to 88
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      mm=0
  88  nt1r=nir(jj)
      nl1r=loc(jj)
      nb1r=nir(ll)
      nm1r=loc(ll)
      nm2r=loc(nn)
  220 if (nt1r-nb1r) 223,222,299
  299 iax = ideks(nt1r)+nb1r
      icx = ideks (iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn
      idq = icq+nm2r
      icq=ideks(idq)+icq+nm1r
      go to 229
  223 iax = ideks(nb1r)+nt1r
      icx=ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      idq = (nm2r-1)*ljn+nl1r
      icq=ideks(idq)+icq
      go to 229
  222 iax=ideks(nt1r+1)
      icx=ideks(iax+1)
      if (nl1r.lt.nm1r) go to 246
      iay=ideks(nl1r)+nm1r
      go to 247
  246 iay= ideks(nm1r)+nl1r
      iby=ideks(nm2r)+nl1r
      go to 253
  247 iby = min(nl1r,nm2r) + ideks( max(nl1r,nm2r))
  253 icq=ideks(iby)+iay
  229 icq=nit(icx)  +  icq
      idud=(icq -1)/igmax
      jto=jto + 1
      kc(jto) = idud  + 1
      icq=icq  -  idud*igmax
      kd(jto)  = icq
      if(jto.lt.iwod)  go to 260
      write(ntype) kc,kd
      jto=0
      jblk=jblk  + 1
      go to 260
   79 if (mr.eq.0) go to 89
      jz=mt
      do 90 l=1,mr
      jz=jz+1
      ii=iot(jz)
      if (jcon(ii+l3)-1) 75,90,91
  90  continue
      go to 89
  91  l1=l
      if (l.eq.mr) go to 93
      ia=l+1
      do 92 l=ia,mr
      jz=jz+1
      jc=iot(jz)
      if (jcon(jc+l3).ne.1) go to 75
  92  continue
      go to 93
  89  mm=mm+1
      olab(mm)=3
      if (mm.lt.nid) go to 94
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      mm=0
   94 jz=l4
      do 95 l=1,nr
      jz=jz+1
      if (iot(jz).eq.jj) go to 96
   95 continue
   96 ia=l+1
      ip=l
      do 97 l=ia,nr
      jz=jz+1
      if (iot(jz).eq.ll) go to 98
   97 continue
  98  mm=mm+1
      olab(mm)=-ip-ideks(l-1)
      if (mm.lt.nid) go to 99
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      mm=0
   99 jz=nz
      if (l1.eq.1) go to 180
      kz=mz
      ja=l1-1
      do 181 l=1,ja
      kz=kz+1
      jz=jz+1
      kk=iot(jz)
      if (kk.ne.iot(kz)) go to 182
  181 continue
  180 if (l2.eq.l1+1) go to 183
      kz=mz+l1
      ia=l1+1
      ja=l2-1
      do 184 l=ia,ja
      kz=kz+1
      jz=jz+1
      kk=iot(jz)
      if (kk.ne.iot(kz)) go to 182
  184 continue
  183 if (l2.lt.md) go to 185
      kk=iot(jz+1)
      go to 182
  185 kz=mz+l2
      ia=l2+1
      do 186 l=ia,md
      kz=kz+1
      jz=jz+1
      kk=iot(jz)
      if (kk.ne.iot(kz)) go to 182
  186 continue
      kk=iot(jz+1)
  182 nt1r=nir(kk)
      nl1r=loc(kk)
      nb1r=nir(jj)
      nm1r=loc(jj)
      nm2r=loc(ll)
      go to 220
  93  mm=mm+1
      olab(mm)=2
      if (mm.lt.nid) go to 187
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      mm=0
 187  jz=l4
      kz=mt
      ig=0
      if (l1.eq.1) go to 188
      ja=l1-1
      ka=1
      do 189 l=ka,ja
      kz=kz+1
      jp=iot(kz)
 802  jz=jz+1
      if (jp.eq.iot(jz)) go to 189
      ig=ig+1
      lab(ig)=jz-l4
      if (ig.eq.3) go to 190
      go to 802
  189 continue
  188 if (l1.eq.mr) go to 191
      ia=l1+1
      kz=kz+1
      do 192 l=ia,mr
      kz=kz+1
      jp=iot(kz)
 803  jz=jz+1
      if (jp.eq.iot(jz)) go to 192
      ig=ig+1
      lab(ig)=jz-l4
      if (ig.eq.3) go to 190
      go to 803
  192 continue
  191 lab(3)=nr
      if (ig-1) 284,285,190
  284 lab(2)=nr1
      lab(1)=mr
      go to 190
  285 lab(2)=nr1
  190 jz=lab(1)+l4
      nn=iot(jz)
      mm=mm+1
      if (nn.eq.jj) go to 194
      olab(mm)=1
  195 if (mm.lt.nid) go to 196
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  196 mm=mm+1
      olab(mm)=l1
      if (mm.lt.nid) go to 197
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  197 mm=mm+1
      jp=lab(2)-1
      kp=lab(3)-2
      olab(mm)=lab(1)+ideks(jp)+jdeks(kp)
      if (mm.lt.nid) go to 198
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      go to 198
  194 jz=lab(2)+l4
      nn=iot(jz)
      if (nn.eq.ll) go to 271
      olab(mm)=2
      go to 195
  271 jz=lab(3)+l4
      nn=iot(jz)
      olab(mm)=3
      go to 195
  198 nt1r=nir(jj)
      nt2r=nir(ll)
      nl1r=loc(jj)
      nl2r=loc(ll)
      kix=3
      nb1r=nir(nn)
      nb2r=nir(ii)
      nm1r=loc(nn)
      nm2r=loc(ii)
  320 if (nt1r-nb1r) 323,322,399
  399 iax = ideks(nt1r)+nb1r
      if (nt2r.lt.nb2r) go to 324
      ibx = ideks(nt2r)+nb2r
      if (iax-ibx) 325,326,327
  327 icx = ideks(iax) + ibx
      icq = (nl1r-1)*lj(nb1r) + nm1r
  340 ljn = lj(nb2r)
      idq = (nl2r-1)*ljn + nm2r
      icq = (icq-1)*ljn*lj(nt2r)+idq
      go to 329
  326 icx = ideks (iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn+nm1r
  341 idq = (nl2r-1)*ljn+nm2r
      go to 333
  325 icx = ideks(ibx) + iax
      icq = (nl2r-1)*lj(nb2r)+nm2r
  335 ljn = lj(nb1r)
      idq = (nl1r-1)*ljn+nm1r
      icq =(icq-1)*ljn*lj(nt1r) +idq
      go to 329
  324 ibx = ideks(nb2r) + nt2r
      if (iax-ibx) 330,331,332
  332 icx = ideks(iax)+ibx
      icq = (nl1r-1)*lj(nb1r)+nm1r
      ljn = lj(nt2r)
      idq = (nm2r-1)*ljn+nl2r
      icq = (icq-1)*ljn*lj(nb2r) + idq
      go to 329
  331 icx = ideks(iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn+nm1r
      idq = (nm2r-1)*ljn+nl2r
  333 icq = min(icq,idq) + ideks( max(icq,idq))
      go to 329
  330 icx=ideks(ibx)+iax
      icq = (nm2r-1)*lj(nt2r)+nl2r
      go to 335
  323 iax = ideks(nb1r)+nt1r
      if (nt2r.lt.nb2r) go to 336
      ibx = ideks(nt2r)+nb2r
      if (iax-ibx) 337,338,339
  339 icx = ideks(iax)+ibx
      icq = (nm1r-1)*lj(nt1r)+nl1r
      go to 340
  338 icx = ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      go to 341
  337 icx = ideks(ibx) + iax
      icq = (nl2r-1)*lj(nb2r) +nm2r
  345 ljn = lj(nt1r)
      idq = (nm1r-1)*ljn         +nl1r
      icq = (icq-1)*ljn*lj(nb1r) +idq
      go to 329
  336 ibx = ideks(nb2r)+nt2r
      if (iax-ibx) 342,343,344
  344 icx = ideks(iax)+ibx
      icq =(nm1r-1) *lj(nt1r)+nl1r
      ljn = lj(nt2r)
      idq = (nm2r-1)*ljn+nl2r
      icq = (icq-1)*ljn*lj(nb2r)+idq
      go to 329
  343 icx=ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      idq = (nm2r-1)*ljn+nl2r
      go to 333
  342 icx = ideks(ibx)+iax
      icq =(nm2r-1)*lj(nt2r)+nl2r
      go to 345
  322 iax=ideks(nt1r+1)
      iay = min(nl1r,nm1r) + ideks( max(nl1r,nm1r))
      iby = min(nl2r,nm2r) + ideks( max(nl2r,nm2r))
      if (nt1r.eq.nt2r) go to 350
      ibx=ideks(nt2r+1)
      if (iax.lt.ibx) go to 351
      icx=ideks(iax)+ibx
      ljn =lj(nt2r)
      icq = (iay-1)*ideks(ljn+1)+iby
      go to 329
  351 icx=ideks(ibx)+iax
      ljn=lj(nt1r)
      icq =(iby-1)*ideks(ljn+1)+iay
      go to 329
  350 icx=ideks(iax+1)
      icq = min(iay,iby) + ideks( max(iay,iby))
  329 icq=nit(icx)  +  icq
      idud=(icq -1)/igmax
      jto=jto + 1
      kc(jto) = idud  + 1
      icq=icq  -  idud*igmax
      kd(jto)  = icq
      if(jto.lt.iwod)  go to 360
      write(ntype) kc,kd
      jto=0
      jblk=jblk  + 1
 360  if (kix.lt.0) go to 260
      kix=-3
      itr=nb1r
      nb1r=nb2r
      nb2r=itr
      itr=nm1r
      nm1r=nm2r
      nm2r=itr
      go to 320
  76  if (mr.eq.0) go to 273
      jz=mt
      do 274 l=1,mr
      jz=jz+1
      ll=iot(jz)
      if (jcon(ll+l3)-1) 275,274,75
 274  continue
      go to 273
 275  l1=l
      if (l.eq.mr) go to 277
      ia=l+1
      do 276 l=ia,mr
      jz=jz+1
      nn=iot(jz)
      if (jcon(nn+l3).ne.1) go to 75
  276 continue
 277  mm=mm+1
      olab(mm)=2
      if (mm.lt.nid) go to 278
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      mm=0
 278  jz=l4
      kz=mt
      ig=0
      if (l1.eq.1) go to 279
      ja=l1-1
      do 280 l=1,ja
      kz=kz+1
      jp=iot(kz)
 254  jz=jz+1
      if (jp.eq.iot(jz)) go to 280
      ig=ig+1
      lab(ig)=jz-l4
      if (ig.eq.3) go to 281
      go to 254
 280  continue
 279  if (l1.eq.mr) go to 282
      ia=l1+1
      kz=kz+1
      do 283 l=ia,mr
      kz=kz+1
      jp=iot(kz)
 255  jz=jz+1
      if (jp.eq.iot(jz)) go to 283
      ig=ig+1
      lab(ig)=jz-l4
      if (ig.eq.3) go to 281
      go to 255
  283 continue
  282 lab(3)=nr
      if (ig-1) 286,287,281
  286 lab(2)=nr1
      lab(1)=mr
      go to 281
  287 lab(2)=nr1
  281 jz=lab(1)+l4
      kk=iot(jz)
      mm=mm+1
      if (kk.ne.jj) go to 288
      olab(mm)=-1
      jz=lab(2)+l4
      kk=iot(jz)
      jz=lab(3)+l4
      nn=iot(jz)
 289  if(mm.lt.nid) go to 290
      mm=0
      call pack(olab8,8,olab,2000)
      write(ntape) olab8
 290  mm=mm+1
      olab(mm)=l1
      if(mm.lt.nid) go to 291
      mm=0
      call pack(olab8,8,olab,2000)
      write(ntape) olab8
 291  mm=mm+1
      jp=lab(2)-1
      kp=lab(3)-2
      olab(mm)=lab(1)+ideks(jp)+jdeks(kp)
      if(mm.lt.nid)go to 292
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
      go to 292
 288  jz=lab(2)+l4
      nn=iot(jz)
      if (nn.ne.jj) go to 293
      olab(mm)=-2
      jz=lab(3)+l4
      nn=iot(jz)
      go to 289
  293 olab(mm)=-3
      go to 289
 292  nt1r=nir(jj)
      nt2r=nir(ll)
      nl1r=loc(jj)
      nl2r=loc(ll)
      kix=3
      nb1r=nir(kk)
      nb2r=nir(nn)
      nm1r=loc(kk)
      nm2r=loc(nn)
      go to 320
 273  mm=mm+1
      olab(mm)=4
      if(mm.lt.nid) go  to 294
      mm=0
      call pack(olab8,8,olab,2000)
      write(ntape) olab8
 294  jz=l4
      mm=mm+1
      if(mr.eq.0) go to 297
      kz=mt
      do 295 l=1,mr
      kz=kz+1
      jz=jz+1
      ll=iot(jz)
      if(ll.ne.iot(kz)) go to 296
 295  continue
 297  ll=iot(jz+1)
      if(ll.eq.jj) go to 298
      olab(mm)=iab
      if(mm.lt.nid) go to 371
      mm=0
      call pack(olab8,8,olab,2000)
      write(ntape) olab8
      go to 371
 298  ll=iot(jz+2)
      olab(mm)=-iab
      if(mm.lt.nid) go to 371
      mm=0
      call pack(olab8,8,olab,2000)
      write(ntape) olab8
      go to 371
 296  l1 =jz-l4
      if(ll.eq.jj) go to 373
      do 374 l=1,nr
      jz=jz+1
      if(iot(jz).eq.jj) go to 375
 374  continue
  375 olab(mm)=l1+ideks(jz-nt4)
_IF1()      itwt=jz-nt1
      if(mm.lt.nid) go to 371
      mm=0
      call pack(olab8,8,olab,2000)
      write(ntape) olab8
      go to 371
 373  kz=kz-1
      do 376 ia=l,mr
      kz=kz+1
      jz=jz+1
      ll=iot(jz)
      if(ll.ne.iot(kz)) go to 377
 376  continue
      jz=jz+1
      ll=iot(jz)
  377 olab(mm)=-l1-ideks(jz-nt4)
      if(mm.lt.nid) go to 371
      mm=0
      call pack(olab8,8,olab,2000)
      write(ntape) olab8
 371  mm=mm+1
      iax=nir(jj)
      olab(mm)=iax
      if(mm.lt.nid) go to 372
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  372 iay=ideks(iax+1)
      iaz=ideks(iay+1)
      mq=lj(iax)
      mv=ideks(mq+1)
      iaq=iay-iax
      iaw=iaz-iay
      iat=nit(iaz)
      ii=loc(jj)
      kk=loc(ll)
      km=ideks(kk)
      nm=ideks(ii)
      nn=nm+ii
      if(ii.gt.kk) go to 380
      mm=mm+1
      jb=ii+km
      olab(mm)=kk
      if(mm.lt.nid) go to 378
      mm=0
      call pack(olab8,8,olab,2000)
      write(ntape) olab8
 378  mm=mm+1
      olab(mm)=ii
      if(mm.lt.nid) go to 379
      mm=0
      call pack(olab8,8,olab,2000)
      write(ntape) olab8
 379  ib=ideks(jb) +nn+iat
      go to 381
 380  lb=kk-ii
      jb=nn+lb
      ib=ideks(nn+1)+lb+iat
      mm=mm+1
      olab(mm)=ii
      if(mm.lt.nid) go to 900
      mm=0
      call pack(olab8,8,olab,2000)
      write(ntape) olab8
 900  mm=mm+1
      olab(mm)=kk
      if(mm.lt.nid) go to 381
      mm=0
      call pack(olab8,8,olab,2000)
      write(ntape) olab8
 381  idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if(jto.lt.iwod) go to 382
      jto=0
      write(ntype) kc,kd
      jblk=jblk+1
 382  if(mr.eq.0) go to 383
      jz=mt
      ir=mr
      kix=1
 472  do 384 l=1,ir
      jz=jz+1
      mi=iot(jz)
      mar=nir(mi)
      mlr=loc(mi)
      mlb=ideks(mlr)
      mla=mlb+mlr
      if (mar-iax) 385,395,396
 395  if(mla.lt.jb) go to 386
      ib=iat+ideks(mla) +jb
 388  idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto = jto+1
      kc(jto)=idud+1
      kd(jto) =ib
      if(jto.lt.iwod) go to 387
      jto =0
      write(ntype) kc,kd
      jblk=jblk+1
      go to 387
 386  ib=iat+ideks(jb) +mla
      go to 388
 387  if(mlr.lt.ii) go to 389
      kb=mlb+ii
      go to 390
 389  kb=nm+mlr
 390  if(mlr.lt.kk) go to 391
      lb=mlb+kk
      go to 392
 391  lb=mlr+km
 392  ib = iat + min(kb,lb) + ideks( max(kb,lb))
      idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if(jto.lt.iwod) go to 384
      jto=0
      jblk=jblk+1
      write(ntype) kc,kd
      go to 384
 385  iby=ideks(mar+1)
      iby=iaw+iby
      ibx=iaq+mar
      ibx=ideks(ibx+1)
      ibx=nit(ibx)
      iby=nit(iby)
      ml=lj(mar)
      ib=ideks(ml+1)*(jb-1) +mla+iby
      idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if(jto.lt.iwod) go to 397
      jto=0
      jblk=jblk+1
      write(ntype) kc,kd
 397  kb=(ii-1)*ml +mlr
      lb=(kk-1)*ml +mlr
      ib = ibx + min(kb,lb) + ideks( max(kb,lb))
      idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto) = ib
      if(jto.lt.iwod) go to 384
      jto=0
      jblk=jblk+1
      write(ntype) kc,kd
      go to 384
396   ibx=ideks(mar+1)
      iby=ideks(ibx) +iay
      iby=nit(iby)
      ibx=ibx-mar+iax
      ibx=ideks(ibx+1)
      ibx=nit(ibx)
      ib=(mla-1)*mv +jb+iby
      idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if(jto.lt.iwod) go to 473
      jto=0
      jblk=jblk+1
      write(ntype) kc,kd
 473  kb=(mlr-1)*mq
      lb=kb+ii
      kb=kb+kk
      ib = ibx + min(kb,lb) + ideks( max(kb,lb))
      idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if (jto.lt.iwod) go to 384
      jto=0
      jblk=jblk+1
      write (ntype) kc,kd
 384  continue
      if (kix.lt.0) go to 260
 383  if (nd.eq.0) go to 260
      kix=-1
      jz=nz
      ir=nd
      go to 472
 260  l4=l4+nx
  65  l3=l3+imo
      return
      end
      subroutine smr110(iot,ideks)
      implicit REAL (a-h,o-z), integer(i-n)
      integer olab
_IF(cray,ksr,i8)
      integer olab8
_ENDIF
INCLUDE(common/sizes)
      common/scrtch/jdeks(9),kdeks(8),
     + lj(8),olab(2000),nit(666),nytl(5),ndub(5),
     + icon(5),nod(5),jdum(maxorb),nir(maxorb),loc(maxorb)
      common/junk/ican(3),icbn,olab8(250),jcon(8200),mkon(mxcrec),
     + kc(mxcrec),kd(mxcrec),lab(3)
     + ,ij(8),jtest(mxnshl),nyzl(mxcsf),jkan(mxcrec),nplu(5),
     + jkon(mxnshl*mxcsf)
      common /b/ iwod,nid,lg,irsf,ntape,ntype,ical,ntab,nc,lt,lz,ntx,
     + ndx,nston,itape,jtape,mdisk,ideli,mtape,iswh,nrx,mm,jblk,jto,
     + igmax,nr2,ltype,linf,kfile,nawk,isec,nshl,
     + nsec,lulu,nrootx,lsng,nprin,ipt0,mconf(5),nconf(5),mh,jswh,
     + id(mxcsf),ie(mxcsf),nsc(mxcsf),
     + imo,mc,ihz,nd,nr,nx,nt,ny,mt,mz,nzx,
     + iht,nz,jy,jx,nd1,nd2,md1,md2,nr3,nr4,mr1,mr2,nt1,nt2,mt1,mt2,
     + ih1,ih2,ih3,ih4,jh2,jh3,jh4,iab,nt3,iac,nr1,m,nad,nad1,ifr1,
     + mx,md,mr,idum1(1)
      dimension iot(*),ideks(*)
      l3=mz
      l4=mt
      do 105 k=1,mc
      if (nd.eq.0) go to 106
      jz=l3
      do 107 l=1,nd
      jz=jz+1
      ll=iot(jz)
      if (jcon(ll)-1) 108,109,107
 107  continue
      go to 106
 108  ip=l
      if (l.eq.nd) go to 110
      ia=l+1
      do 111 l=ia,nd
      jz=jz+1
      jj=iot(jz)
      if (jcon(jj).ne.2) go to 112
  111 continue
      go to 110
  112 mm=mm+1
      olab(mm)=1
      if (mm.lt.nid) go to 113
      mm=0
      call pack(olab8,8,olab,2000)
      write  (ntape) olab8
      goto 113
 109  ip=l
      if (l.eq.nd) go to 114
      ia=l+1
      do 115 l=ia,nd
      jz=jz+1
      nn=iot(jz)
      if (jcon(nn)-1) 112,116,115
115   continue
      go to 114
116   if (l.eq.nd) go to 117
      ia=l+1
      do 118 l=ia,nd
      jz=jz+1
      kk=iot(jz)
      if (jcon(kk).ne.2)  go to 112
 118  continue
      go to 117
 110  if (nr.eq.0) go to 119
      jz=l4
      do 120 l=1,nr
      jz=jz+1
      jj=iot(jz)
      if (jcon(jj)-1) 112,120,121
 120  continue
      go to 119
 121  ip=l
      if (l.eq.nr) go to 122
      ia=l+1
      do 123 l=ia,nr
      jz=jz+1
      kk=iot(jz)
      if (jcon(kk).ne.1) go to 112
 123  continue
      go to 122
 119  mm=mm+1
      olab(mm)=5
      if (mm.lt.nid) go to 124
      mm=0
      call pack(olab8,8,olab,2000)
      write  (ntape) olab8
 124  jz=nz
      if (ip.eq.1) go to 125
      kz=l3
      ja=ip-1
      do 126 l=1,ja
      kz=kz+1
      jz=jz+1
      kk=iot(jz)
      if (kk.ne.iot(kz)) go to 127
 126  continue
 125  if (ip.eq.nd) go to 128
      ja=ip+1
      kz=l3+ip
      do 129 l=ja,nd
      kz=kz+1
      jz=jz+1
      kk=iot(jz)
      if (kk.ne.iot(kz)) go to 127
 129  continue
 128  kk=iot(jz+1)
  127 mm=mm+1
      olab(mm)=ll
      if (mm.lt.nid) go to 130
      mm=0
      call pack(olab8,8,olab,2000)
      write  (ntape) olab8
  130 mm=mm+1
      olab(mm)=kk
      if (mm.lt.nid) go to 113
      mm=0
      call pack(olab8,8,olab,2000)
      write  (ntape) olab8
      go to 113
  122 mm=mm+1
      olab(mm)=3
      if (mm.lt.nid) go to 131
      mm=0
      call pack(olab8,8,olab,2000)
      write  (ntape) olab8
  131 jz=nt
      if (ip.eq.1) go to 132
      kz=l4
      ja=ip-1
      do 133 l=1,ja
      kz=kz+1
      jz=jz+1
      kk=iot(jz)
      if (kk.ne.iot(kz)) go to 134
 133  continue
 132  if (ip.eq.nr) go to 650
      ja=ip+1
      kz=l4+ip
      do 651 l=ja,nr
      kz=kz+1
      jz=jz+1
      kk=iot(jz)
      if (kk.ne.iot(kz)) go to 134
 651  continue
 650  jz=jz+1
      kk=iot(jz)
  134 mm=mm+1
      olab(mm)=jz-nt
      if (mm.lt.nid) go to 135
      mm=0
      call pack(olab8,8,olab,2000)
      write  (ntape) olab8
  135 mm=mm+1
      olab(mm)=ip
      if (mm.lt.nid) go to 136
      mm=0
      call pack(olab8,8,olab,2000)
      write  (ntape) olab8
  136 nt1r=nir(ll)
      nl1r=loc(ll)
      nb1r=nir(jj)
      nm1r=loc(jj)
      nm2r=loc(kk)
   20 if (nt1r-nb1r) 23,22,99
   99 iax = ideks(nt1r)+nb1r
      icx = ideks (iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn
      idq=icq+nm2r
      icq=icq+nm1r
   33 icq =  min(icq,idq) + ideks( max(icq,idq))
      go to 29
   23 iax=ideks(nb1r)+nt1r
      icx=ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      idq = (nm2r-1)*ljn+nl1r
      go to 33
   22 iax=ideks(nt1r+1)
      iay = min(nl1r,nm1r) + ideks( max(nl1r,nm1r))
      iby = min(nl1r,nm2r) + ideks( max(nl1r,nm2r))
      icx=ideks(iax+1)
      icq = min(iay,iby) + ideks( max(iay,iby))
 29   icq=nit(icx)  +  icq
      idud=(icq -1)/igmax
      jto=jto + 1
      kc(jto) = idud  + 1
      icq=icq  -  idud*igmax
      kd(jto)  = icq
      if(jto.lt.iwod)  go to 113
      write(ntype) kc,kd
      jto=0
      jblk=jblk  + 1
      go to 113
 117  jz=l4
      do 137 l=1,nr
      jz=jz+1
      kk=iot(jz)
      if (jcon(kk)-1) 112,137,138
  137 continue
  138 l1=l
      ia=l+1
      do 139 l=ia,nr
      jz=jz+1
      jj=iot(jz)
      if (jcon(jj)-1) 112,139,140
  139 continue
  140 l2=l
      if (l.eq.nr) go to 141
      ia=l+1
      do 142 l=ia,nr
      jz=jz+1
      ii=iot(jz)
      if (jcon(ii).ne.1) go to 112
 142  continue
 141  mm=mm+1
      olab(mm)=2
      if (mm.lt.nid) go to 143
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  143 mm=mm+1
      olab(mm)=1
      if (mm.lt.nid) go to 144
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
 144  jz=nt
      do 145 l=1,nr
      jz=jz+1
      if (iot(jz).eq.ll) go to 146
 145  continue
 146  ip=l
      ia=l+1
      do 870 l=ia,nr
      jz=jz+1
      if (iot(jz).eq.nn) go to 147
 870  continue
 147  mm=mm+1
      olab(mm)=ip+ideks(l-1)
      if (mm.lt.nid) go to 148
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  148 mm=mm+1
      olab(mm)=ideks(l2-1)+l1
 540  if (mm.lt.nid) go to 149
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
 149  nt1r=nir(ll)
      nt2r=nir(nn)
      nl1r=loc(ll)
      nl2r=loc(nn)
 189  nb1r=nir(kk)
      nb2r=nir(jj)
      nm1r=loc(kk)
      nm2r=loc(jj)
      kix=3
  220 if (nt1r-nb1r) 223,222,299
  299 iax = ideks(nt1r)+nb1r
      if (nt2r.lt.nb2r) go to 224
      ibx = ideks(nt2r)+nb2r
      if (iax-ibx) 225,226,227
  227 icx = ideks(iax) + ibx
      icq = (nl1r-1)*lj(nb1r) + nm1r
  240 ljn = lj(nb2r)
      idq = (nl2r-1)*ljn + nm2r
      icq = (icq-1)*ljn*lj(nt2r)+idq
      go to 229
  226 icx = ideks (iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn+nm1r
  241 idq = (nl2r-1)*ljn+nm2r
      go to 233
  225 icx = ideks(ibx) + iax
      icq = (nl2r-1)*lj(nb2r)+nm2r
  235 ljn = lj(nb1r)
      idq = (nl1r-1)*ljn+nm1r
      icq =(icq-1)*ljn*lj(nt1r) +idq
      go to 229
  224 ibx = ideks(nb2r) + nt2r
      if (iax-ibx) 230,231,232
  232 icx = ideks(iax)+ibx
      icq = (nl1r-1)*lj(nb1r)+nm1r
      ljn = lj(nt2r)
      idq = (nm2r-1)*ljn+nl2r
      icq = (icq-1)*ljn*lj(nb2r) + idq
      go to 229
  231 icx = ideks(iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn+nm1r
      idq = (nm2r-1)*ljn+nl2r
  233 icq =  min(icq,idq) + ideks( max(icq,idq))
      go to 229
  230 icx=ideks(ibx)+iax
      icq = (nm2r-1)*lj(nt2r)+nl2r
      go to 235
  223 iax = ideks(nb1r)+nt1r
      if (nt2r.lt.nb2r) go to 236
      ibx = ideks(nt2r)+nb2r
      if (iax-ibx) 237,238,239
  239 icx = ideks(iax)+ibx
      icq = (nm1r-1)*lj(nt1r)+nl1r
      go to 240
  238 icx = ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      go to 241
  237 icx = ideks(ibx) + iax
      icq = (nl2r-1)*lj(nb2r) +nm2r
  245 ljn = lj(nt1r)
      idq = (nm1r-1)*ljn         +nl1r
      icq = (icq-1)*ljn*lj(nb1r) +idq
      go to 229
  236 ibx = ideks(nb2r)+nt2r
      if (iax-ibx) 242,243,244
  244 icx = ideks(iax)+ibx
      icq =(nm1r-1) *lj(nt1r)+nl1r
      ljn = lj(nt2r)
      idq = (nm2r-1)*ljn+nl2r
      icq = (icq-1)*ljn*lj(nb2r)+idq
      go to 229
  243 icx=ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      idq = (nm2r-1)*ljn+nl2r
      go to 233
  242 icx = ideks(ibx)+iax
      icq =(nm2r-1)*lj(nt2r)+nl2r
      go to 245
  222 iax=ideks(nt1r+1)
      iay =  min(nl1r,nm1r) + ideks( max(nl1r,nm1r))
      iby =  min(nl2r,nm2r) + ideks( max(nl2r,nm2r))
      if (nt1r.eq.nt2r) go to 250
      ibx=ideks(nt2r+1)
      if (iax.lt.ibx) go to 251
      icx=ideks(iax)+ibx
      ljn =lj(nt2r)
      icq = (iay-1)*ideks(ljn+1)+iby
      go to 229
  251 icx=ideks(ibx)+iax
      ljn=lj(nt1r)
      icq =(iby-1)*ideks(ljn+1)+iay
      go to 229
  250 icx=ideks(iax+1)
      icq =  min(iay,iby) + ideks( max(iay,iby))
 229  icq=nit(icx)  +  icq
      idud=(icq -1)/igmax
      jto=jto + 1
      kc(jto) = idud  + 1
      icq=icq  -  idud*igmax
      kd(jto)  = icq
      if(jto.lt.iwod)  go to 260
      write(ntype) kc,kd
      jto=0
      jblk=jblk  + 1
260   if(kix.lt.0)  go to 113
      kix=-3
      itr=nb1r
      nb1r=nb2r
      nb2r=itr
      itr=nm1r
      nm1r=nm2r
      nm2r=itr
      go to 220
 114  jz=l4
      do 150 l=1,nr
      jz=jz+1
      nn=iot(jz)
      if (jcon(nn)-1) 151,150,152
  150 continue
  151 jp=l
      if (l.eq.nr) go to 153
      ia=l+1
      do 154 l=ia,nr
      jz=jz+1
      kk=iot(jz)
      if (jcon(kk)-1) 112,154,155
  154 continue
      go to 153
  155 kp=l
      if (l.eq.nr) go to 156
      ia=l+1
      do 157 l=ia,nr
      jz=jz+1
      ii=iot(jz)
      if (jcon(ii).ne.1) go to 112
  157 continue
      go to 156
  152 jp=l
      if (l.eq.nr) go to 158
      ia=l+1
      do 159 l=ia,nr
      jz=jz+1
      kk=iot(jz)
      if (jcon(kk)-1) 160,159,112
  159 continue
      go to 158
  160 kp=l
      if (l.eq.nr) go to 161
      ia=l+1
      do 162 l=ia,nr
      jz=jz+1
      ii=iot(jz)
      if (jcon(ii).ne.1) go to 112
  162 continue
      go to 161
  153 mm=mm+1
      olab(mm)=3
      if (mm.lt.nid) go to 163
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
 163  jz=nt
      do 164 l=1,nr
      jz=jz+1
      if (iot(jz).eq.ll) go to 165
 164  continue
  165 mm=mm+1
      olab(mm)=l
      if (mm.lt.nid) go to 166
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  166 mm=mm+1
      olab(mm)=jp
      if (mm.lt.nid) go to 167
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  167 jz=nz
      if (ip.eq.1) go to 168
      kz=l3
      ja=ip-1
      do 169 l=1,ja
      jz=jz+1
      kz=kz+1
      kk=iot(jz)
      if (kk.ne.iot(kz)) go to 170
  169 continue
  168 if (ip.eq.nd) go to 171
      ja=ip+1
      kz=l3+ip
      do 172 l=ja,nd
      jz=jz+1
      kz=kz+1
      kk=iot(jz)
      if (kk.ne.iot(kz)) go to 170
  172 continue
  171 kk=iot(jz+1)
  170 nt1r=nir(kk)
      nl1r=loc(kk)
      nb1r=nir(ll)
      nm1r=loc(ll)
      nm2r=loc(nn)
      go to 20
  156 mm=mm+1
      olab(mm)=2
      if (mm.lt.nid) go to 173
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  173 jz=nt
      mm=mm+1
      kz=l4
      ig=0
      if (jp.eq.1) go to 174
      ka=1
      ja=jp-1
  176 do 175 l=ka,ja
      jz=jz+1
      kz=kz+1
      if (iot(jz).eq.iot(kz)) go to 175
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.2) go to 177
      go to 800
 175  continue
      go to 174
 800  ka=l
      kz=kz-1
      go to 176
 174  la=kp-1
      if (jp.eq.la) go to 178
      ja=jp+1
      kz=kz+1
 180  do 179 l=ja,la
      kz=kz+1
      jz=jz+1
      if (iot(kz).eq.iot(jz)) go to 179
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.2) go to 177
      go to 801
 179  continue
      go to 178
 801  ja=l
      kz=kz-1
      go to 180
 178  if (kp.eq.nr) go to 181
      ia=kp+1
      kz=l4+kp
 182  do 183 l=ia,nr
      jz=jz+1
      kz=kz+1
      if (iot(jz).eq.iot(kz)) go to 183
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.2) go to 177
      go to 802
 183  continue
      go to 181
 802  ia=l
      kz=kz-1
      go to 182
 181  lab(2)=nr
      if (ig.eq.1) go to 177
      lab(1)=nr1
 177  ia=lab(1)
      jz=ia+nt
      ja=lab(2)
      jj=iot(jz)
      if (jj.eq.ll) go to 184
      olab(mm)=2
 188  if (mm.lt.nid) go to 185
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  185 mm=mm+1
      olab(mm)=ideks(ja-1)+ia
      if (mm.lt.nid) go to 186
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  186 mm=mm+1
      olab(mm)=ideks(kp-1)+jp
      if (mm.lt.nid) go to 187
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
 187  nt1r=nir(nn)
      nt2r=nir(ll)
      nl1r=loc(nn)
      nl2r=loc(ll)
      go to 189
 184  olab(mm)=3
      jz=nt+ja
      jj=iot(jz)
      go to 188
  161 mm=mm+1
      olab(mm)=2
      if (mm.lt.nid) go to 190
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  190 jz=nt
      mm=mm+1
      kz=l4
      ig=0
      if (jp.eq.1) go to 191
      ka=1
      ja=jp-1
 192  do 193 l=ka,ja
      jz=jz+1
      kz=kz+1
      if (iot(jz).eq.iot(kz)) go to 193
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.2) go to 806
      go to 803
 193  continue
      go to 191
 803  ka=l
      kz=kz-1
      go to 192
 191  la=kp-1
      if (jp.eq.la) go to 194
      ja=jp+1
      kz=kz+1
 195  do 196 l=ja,la
      kz=kz+1
      jz=jz+1
      if (iot(kz).eq.iot(jz)) go to 196
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.2) go to 806
      go to 804
196   continue
      go to 194
 804  ja=l
      kz=kz-1
      go to 195
 194  if (kp.eq.nr) go to 197
      ia=kp+1
      kz=l4+kp
 198  do 199 l=ia,nr
      kz=kz+1
      jz=jz+1
      if (iot(jz).eq.iot(kz)) go to 199
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.2) go to 806
      go to 805
  199 continue
      go to 197
  805 ia=l
      kz=kz-1
      go to 198
  197 lab(2)=nr
      if (ig.eq.1) go to 806
      lab(1)=nr-1
  806 ia=lab(1)
      jz=ia+nt
      ja=lab(2)
      jj=iot(jz)
      if (jj.eq.ll) go to 501
      olab(mm)=3
 502  if (mm.lt.nid) go to 503
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  503 mm=mm+1
      olab(mm)=ideks(ja-1)+ia
      if (mm.lt.nid) go to 504
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  504 mm=mm+1
      olab(mm)=ideks(kp-1)+jp
      if (mm.lt.nid) go to 505
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
 505  nt1r=nir(nn)
      nt2r=nir(jj)
      nl1r=loc(nn)
      nl2r=loc(jj)
      kix=3
      nb1r=nir(kk)
      nb2r=nir(ll)
      nm1r=loc(kk)
      nm2r=loc(ll)
      go to 220
 501  olab(mm)=2
      jz=nt+ja
      jj=iot(jz)
      go to 502
 158  mm=mm+1
      olab(mm)=4
      if (mm.lt.nid) go to 506
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
 506  jz=nt
      do 507 l=1,nr
      jz=jz+1
      if (iot(jz).eq.ll) go to 508
 507  continue
 508  mm=mm+1
      olab(mm)=-l
      if (mm.lt.nid) go to 509
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
 509  mm=mm+1
      olab(mm)=jp
      if (mm.lt.nid) go to 510
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
 510  mm=mm+1
      iax=nir(ll)
      iay=ideks(iax+1)
      iaz=ideks(iay+1)
      olab(mm)=iax
      if (mm.lt.nid) go to 893
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
 893  mm=mm+1
      mq=lj(iax)
      mv=ideks(mq+1)
      iaq=iay-iax
      iaw=iaz-iay
      iat=nit(iaz)
      ii=loc(ll)
      kk=loc(nn)
      km=ideks(kk)
      kl=km+kk
      nm=ideks(ii)
      nn=nm+ii
      if (ii.gt.kk) go to 513
      olab(mm)=kk
      if (mm.lt.nid) go to 511
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
 511  mm=mm+1
      olab(mm)=ii
      if (mm.lt.nid) go to 512
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
 512  jb=km+ii
      ib=ideks(jb)+nn+iat
      kb=ideks(kl)+jb+iat
      go to 514
 513  jb=nm+kk
      olab(mm)=ii
      if (mm.lt.nid) go to 551
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
 551  mm=mm+1
      olab(mm)=kk
      if (mm.lt.nid) go to 552
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
 552  ib=ideks(jb)+kl+iat
      kb=ideks(nn)+jb+iat
 514  idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if (jto.lt.iwod) go to 515
      jto=0
      write (ntype) kc,kd
      jblk=jblk+1
 515  idud=(kb-1)/igmax
      kb=kb-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=kb
      nix=1
      if (jto.lt.iwod) go to 516
      jto=0
      write (ntype) kc,kd
      jblk=jblk+1
 516  if (nr.eq.1) go to 517
      jz=l4
      kix=1
      lp=jp
      ir=nr
 518  do 384 l=1,ir
      jz=jz+1
      if (l.eq.lp) go to 384
      mi=iot(jz)
      mar=nir(mi)
      mlr=loc(mi)
      mlb=ideks(mlr)
      mla=mlb+mlr
      if (mar-iax) 385,395,396
 395  if(mla.lt.jb) go to 386
      ib=iat+ideks(mla) +jb
 388  idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto = jto+1
      kc(jto)=idud+1
      kd(jto) =ib
      if(jto.lt.iwod) go to 387
      jto =0
      write(ntype) kc,kd
      jblk=jblk+1
      go to 387
 386  ib=iat+ideks(jb) +mla
      go to 388
 387  if(mlr.lt.ii) go to 389
      kb=mlb+ii
      go to 390
 389  kb=nm+mlr
 390  if(mlr.lt.kk) go to 391
      lb=mlb+kk
      go to 392
 391  lb=mlr+km
 392  ib = iat + min(kb,lb) + ideks( max(kb,lb))
      idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if(jto.lt.iwod) go to 384
      jto=0
      jblk=jblk+1
      write(ntype) kc,kd
      go to 384
 385  iby=ideks(mar+1)
      iby=iaw+iby
      ibx=iaq+mar
      ibx=ideks(ibx+1)
      ibx=nit(ibx)
      iby=nit(iby)
      ml=lj(mar)
      ib=ideks(ml+1)*(jb-1) +mla+iby
      idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if(jto.lt.iwod) go to 397
      jto=0
      jblk=jblk+1
      write(ntype) kc,kd
 397  kb=(ii-1)*ml +mlr
      lb=(kk-1)*ml +mlr
      ib = ibx + min(kb,lb) + ideks( max(kb,lb))
      idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto) = ib
      if(jto.lt.iwod) go to 384
      jto=0
      jblk=jblk+1
      write(ntype) kc,kd
      go to 384
396   ibx=ideks(mar+1)
      iby=ideks(ibx) +iay
      iby=nit(iby)
      ibx=ibx-mar+iax
      ibx=ideks(ibx+1)
      ibx=nit(ibx)
      ib=(mla-1)*mv +jb+iby
      idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if(jto.lt.iwod) go to 473
      jto=0
      jblk=jblk+1
      write(ntype) kc,kd
 473  kb=(mlr-1)*mq
      lb=kb+ii
      kb=kb+kk
      ib = ibx + min(kb,lb) + ideks( max(kb,lb))
      idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if (jto.lt.iwod) go to 384
      jto=0
      jblk=jblk+1
      write (ntype) kc,kd
 384  continue
      if (kix.lt.0) go to 113
      if (nix.lt.0) go to 519
 517  if (nd.eq.1) go to 113
      kix=-1
      jz=l3
      lp=ip
      ir=nd
      go to 518
  519 if (nd.eq.0) go to 113
      kix=-1
      jz=l3
      lp=0
      ir=nd
      go to 518
  106 jz=l4
      do 520 l=1,nr
      jz=jz+1
      ll=iot(jz)
      if (jcon(ll)-1) 521,520,112
  520 continue
  521 ip=l
      if (l.eq.nr) go to 522
      ia=l+1
      do 523 l=ia,nr
      jz=jz+1
      nn=iot(jz)
      if (jcon(nn)-1) 524,523,112
  523 continue
      go to 522
  524 jp=l
      if (l.eq.nr) go to 525
      ia=l+1
      do 526 l=ia,nr
      jz=jz+1
      ii=iot(jz)
      if (jcon(ii).ne.1) go to 112
  526 continue
  525 mm=mm+1
      olab(mm)=2
      if (mm.lt.nid) go to 527
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  527 mm=mm+1
      olab(mm)=1
      if (mm.lt.nid) go to 528
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  528 mm=mm+1
      jz=nt
      kz=l4
      ig=0
      if (ip.eq.1) go to 529
      ka=1
      ja=ip-1
  531 do 530 l=ka,ja
      jz=jz+1
      kz=kz+1
      if (iot(jz).eq.iot(kz)) go to 530
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.2) go to 532
      go to 807
  530 continue
      go to 529
  807 ka=l
      kz=kz-1
      go to 531
  529 la=jp-1
      if (ip.eq.la) go to 533
      ja=ip+1
      kz=kz+1
 534  do 535 l=ja,la
      jz=jz+1
      kz=kz+1
      if (iot(jz).eq.iot(kz)) go to 535
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.2) go to 532
      go to 808
  535 continue
      go to 533
  808 ja=l
      kz=kz-1
      go to 534
  533 if (jp.eq.nr) go to 536
      ia=jp+1
      kz=l4+jp
  537 do 538 l=ia,nr
      kz=kz+1
      jz=jz+1
      if (iot(jz).eq.iot(kz)) go to 538
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.2) go to 532
      go to 809
  538 continue
      go to 536
  809 ia=l
      kz=kz-1
      go to 537
  536 lab(2)=nr
      if (ig.eq.1) go to 532
      lab(1)=nr1
  532 l1=lab(1)
      jz=nt+l1
      kk=iot(jz)
      l2=lab(2)
      jz=nt+l2
      jj=iot(jz)
      olab(mm)=ideks(l2-1)+l1
      if (mm.lt.nid) go to 539
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  539 mm=mm+1
      olab(mm)=ideks(jp-1)+ip
      go to 540
  522 mm=mm+1
      olab(mm)=4
      if (mm.lt.nid) go to 541
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  541 mm=mm+1
      jz=nt
      if (nr.eq.1) go to 546
      kz=l4
      if (ip.eq.1) go to 542
      ja=ip-1
      do 543 l=1,ja
      jz=jz+1
      kz=kz+1
      nn=iot(jz)
      if (nn.ne.iot(kz)) go to 544
  543 continue
      if (ip.eq.nr) go to 546
  542 ja=ip+1
      kz=kz+1
      do 547 l=ja,nr
      jz=jz+1
      kz=kz+1
      nn=iot(jz)
      if (nn.ne.iot(kz)) go to 544
  547 continue
  546 nn=iot(jz+1)
      l1=nr
      olab(mm)=l1
      go to 548
  544 l1=jz-nt
  548 olab(mm)=l1
      if (mm.lt.nid) go to 549
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  549 mm=mm+1
      olab(mm)=ip
      if (mm.lt.nid) go to 550
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  550 mm=mm+1
      iax=nir(ll)
      olab(mm)=iax
      if (mm.lt.nid) go to 897
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
 897  mm=mm+1
      iay=ideks(iax+1)
      iaz=ideks(iay+1)
      mq=lj(iax)
      mv=ideks(mq+1)
      iaq=iay-iax
      iaw=iaz-iay
      iat=nit(iaz)
      ii=loc(ll)
      kk=loc(nn)
      km=ideks(kk)
      nm=ideks(ii)
      if (ii.gt.kk) go to 560
      olab(mm)=kk
      if (mm.lt.nid) go to 561
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  561 mm=mm+1
      olab(mm)=ii
      if (mm.lt.nid) go to 562
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  562 jb=km+ii
      go to 565
  560 olab(mm)=ii
      if (mm.lt.nid) go to 563
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  563 mm=mm+1
      olab(mm)=kk
      if (mm.lt.nid) go to 564
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  564 jb=nm+kk
  565 nix=-1
      if (nr.eq.1) go to 519
      jz=l4
      kix=1
      lp=ip
      ir=nr
      go to 518
  113 l3=l3+nx
  105 l4=l4+nx
      return
      end
      subroutine smr100(iot,ideks)
      implicit REAL (a-h,o-z), integer(i-n)
      integer olab
_IF(cray,ksr,i8)
      integer olab8
_ENDIF
INCLUDE(common/sizes)
      common/scrtch/jdeks(9),kdeks(8),
     + lj(8),olab(2000),nit(666),nytl(5),ndub(5),
     + icon(5),nod(5),jcon(maxorb),nir(maxorb),loc(maxorb)
      common/junk/ican(3),icbn,olab8(250),jdum(8200),mkon(mxcrec),
     + kc(mxcrec),kd(mxcrec),lab(3)
     +,ij(8),jtest(mxnshl),nyzl(mxcsf),jkan(mxcrec),nplu(5),
     + jkon(mxnshl*mxcsf)
      common /b/ iwod,nid,lg,irsf,ntape,ntype,ical,ntab,nc,lt,lz,nx,
     + nd,nston,itape,jtape,mdisk,ideli,mtape,iswh,nr,mm,jblk,jto,
     + igmax,nr1,ltype,linf,kfile,nawk,isec,nshl,
     + nsec,lulu,nrootx,lsng,nprin,ipt0,mconf(5),nconf(5),mh,jswh,
     + id(mxcsf),ie(mxcsf),nsc(mxcsf),
     + imo,mc,mz,ndx,nrx,ntx,ny,nt,iht,ihz,nz,mt,nzx,jy,
     + jx,nd1,nd2,md1,md2,nr3,nr4,mr1,mr2,nt1,nt2,mt1,mt2,ih1,ih2,ih3,
     + ih4,jh2,jh3,jh4,iab,nt3,iac,nr2,m,nad,nad1,ifr1,mx,md,mr
     +,idum1(1)
      dimension iot(*),ideks(*)
      if (nc.eq.1) return
      nt=lt+nx
      nz=lz+nx
      j1=0
      do 100 j=2,nc
      j1=j1+1
      it=nt
      mt=lt
      mz=lz
      if (nr.eq.0) go to 102
      do 101 k=1,nr
      it=it+1
      ll=iot(it)
 101  jcon(ll)=1
      if (nd.eq.0) go to 103
 102  do 104 k=1,nd
      it=it+1
      ll=iot(it)
 104  jcon(ll)=2
 103  do 105 k=1,j1
      if (nd.eq.0) go to 106
      jz=mz
      do 107 l=1,nd
      jz=jz+1
      ll=iot(jz)
      if (jcon(ll)-1) 108,109,107
 107  continue
      go to 106
 108  ip=l
      if (l.eq.nd) go to 110
      ia=l+1
      do 111 l=ia,nd
      jz=jz+1
      jj=iot(jz)
      if (jcon(jj).ne.2) go to 112
  111 continue
      go to 110
  112 mm=mm+1
      olab(mm)=1
      if (mm.lt.nid) go to 113
      mm=0
      call pack(olab8,8,olab,2000)
      write  (ntape) olab8
      goto 113
 109  ip=l
      if (l.eq.nd) go to 114
      ia=l+1
      do 115 l=ia,nd
      jz=jz+1
      nn=iot(jz)
      if (jcon(nn)-1) 112,116,115
115   continue
      go to 114
116   if (l.eq.nd) go to 117
      ia=l+1
      do 118 l=ia,nd
      jz=jz+1
      kk=iot(jz)
      if (jcon(kk).ne.2)  go to 112
 118  continue
      go to 117
 110  if (nr.eq.0) go to 119
      jz=mt
      do 120 l=1,nr
      jz=jz+1
      jj=iot(jz)
      if (jcon(jj)-1) 112,120,121
 120  continue
      go to 119
 121  ip=l
      if (l.eq.nr) go to 122
      ia=l+1
      do 123 l=ia,nr
      jz=jz+1
      kk=iot(jz)
      if (jcon(kk).ne.1) go to 112
 123  continue
      go to 122
 119  mm=mm+1
      olab(mm)=5
      if (mm.lt.nid) go to 124
      mm=0
      call pack(olab8,8,olab,2000)
      write  (ntape) olab8
 124  jz=nz
      if (ip.eq.1) go to 125
      kz=mz
      ja=ip-1
      do 126 l=1,ja
      kz=kz+1
      jz=jz+1
      kk=iot(jz)
      if (kk.ne.iot(kz)) go to 127
 126  continue
 125  if (ip.eq.nd) go to 128
      ja=ip+1
      kz=mz+ip
      do 129 l=ja,nd
      kz=kz+1
      jz=jz+1
      kk=iot(jz)
      if (kk.ne.iot(kz)) go to 127
 129  continue
 128  kk=iot(jz+1)
  127 mm=mm+1
      olab(mm)=ll
      if (mm.lt.nid) go to 130
      mm=0
      call pack(olab8,8,olab,2000)
      write  (ntape) olab8
  130 mm=mm+1
      olab(mm)=kk
      if (mm.lt.nid) go to 113
      mm=0
      call pack(olab8,8,olab,2000)
      write  (ntape) olab8
      go to 113
  122 mm=mm+1
      olab(mm)=3
      if (mm.lt.nid) go to 131
      mm=0
      call pack(olab8,8,olab,2000)
      write  (ntape) olab8
  131 jz=nt
      if (ip.eq.1) go to 132
      kz=mt
      ja=ip-1
      do 133 l=1,ja
      kz=kz+1
      jz=jz+1
      kk=iot(jz)
      if (kk.ne.iot(kz)) go to 134
 133  continue
 132  if (ip.eq.nr) go to 650
      ja=ip+1
      kz=mt+ip
      do 651 l=ja,nr
      kz=kz+1
      jz=jz+1
      kk=iot(jz)
      if (kk.ne.iot(kz)) go to 134
 651  continue
 650  jz=jz+1
      kk=iot(jz)
  134 mm=mm+1
      olab(mm)=jz-nt
      if (mm.lt.nid) go to 135
      mm=0
      call pack(olab8,8,olab,2000)
      write  (ntape) olab8
  135 mm=mm+1
      olab(mm)=ip
      if (mm.lt.nid) go to 136
      mm=0
      call pack(olab8,8,olab,2000)
      write  (ntape) olab8
  136 nt1r=nir(ll)
      nl1r=loc(ll)
      nb1r=nir(jj)
      nm1r=loc(jj)
      nm2r=loc(kk)
   20 if (nt1r-nb1r) 23,22,99
   99 iax = ideks(nt1r)+nb1r
      icx = ideks (iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn
      idq=icq+nm2r
      icq=icq+nm1r
   33 icq =  min(icq,idq) + ideks( max(icq,idq))
      go to 29
   23 iax=ideks(nb1r)+nt1r
      icx=ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      idq = (nm2r-1)*ljn+nl1r
      go to 33
   22 iax=ideks(nt1r+1)
      iay =  min(nl1r,nm1r) + ideks( max(nl1r,nm1r))
      iby =  min(nl1r,nm2r) + ideks( max(nl1r,nm2r))
      icx=ideks(iax+1)
      icq =  min(iay,iby) + ideks( max(iay,iby))
 29   icq=nit(icx)  +  icq
      idud=(icq -1)/igmax
      jto=jto + 1
      kc(jto) = idud  + 1
      icq=icq  -  idud*igmax
      kd(jto)  = icq
      if(jto.lt.iwod)  go to 113
      write(ntype) kc,kd
      jto=0
      jblk=jblk  + 1
      go to 113
 117  jz=mt
      do 137 l=1,nr
      jz=jz+1
      kk=iot(jz)
      if (jcon(kk)-1) 112,137,138
  137 continue
  138 l1=l
      ia=l+1
      do 139 l=ia,nr
      jz=jz+1
      jj=iot(jz)
      if (jcon(jj)-1) 112,139,140
  139 continue
  140 l2=l
      if (l.eq.nr) go to 141
      ia=l+1
      do 142 l=ia,nr
      jz=jz+1
      ii=iot(jz)
      if (jcon(ii).ne.1) go to 112
 142  continue
 141  mm=mm+1
      olab(mm)=2
      if (mm.lt.nid) go to 143
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  143 mm=mm+1
      olab(mm)=1
      if (mm.lt.nid) go to 144
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
 144  jz=nt
      do 145 l=1,nr
      jz=jz+1
      if (iot(jz).eq.ll) go to 146
 145  continue
 146  ip=l
      ia=l+1
      do 870 l=ia,nr
      jz=jz+1
      if (iot(jz).eq.nn) go to 147
 870  continue
 147  mm=mm+1
      olab(mm)=ip+ideks(l-1)
      if (mm.lt.nid) go to 148
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  148 mm=mm+1
      olab(mm)=ideks(l2-1)+l1
 540  if (mm.lt.nid) go to 149
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
 149  nt1r=nir(ll)
      nt2r=nir(nn)
      nl1r=loc(ll)
      nl2r=loc(nn)
 189  nb1r=nir(kk)
      nb2r=nir(jj)
      nm1r=loc(kk)
      nm2r=loc(jj)
      kix=3
  220 if (nt1r-nb1r) 223,222,299
  299 iax = ideks(nt1r)+nb1r
      if (nt2r.lt.nb2r) go to 224
      ibx = ideks(nt2r)+nb2r
      if (iax-ibx) 225,226,227
  227 icx = ideks(iax) + ibx
      icq = (nl1r-1)*lj(nb1r) + nm1r
  240 ljn = lj(nb2r)
      idq = (nl2r-1)*ljn + nm2r
      icq = (icq-1)*ljn*lj(nt2r)+idq
      go to 229
  226 icx = ideks (iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn+nm1r
  241 idq = (nl2r-1)*ljn+nm2r
      go to 233
  225 icx = ideks(ibx) + iax
      icq = (nl2r-1)*lj(nb2r)+nm2r
  235 ljn = lj(nb1r)
      idq = (nl1r-1)*ljn+nm1r
      icq =(icq-1)*ljn*lj(nt1r) +idq
      go to 229
  224 ibx = ideks(nb2r) + nt2r
      if (iax-ibx) 230,231,232
  232 icx = ideks(iax)+ibx
      icq = (nl1r-1)*lj(nb1r)+nm1r
      ljn = lj(nt2r)
      idq = (nm2r-1)*ljn+nl2r
      icq = (icq-1)*ljn*lj(nb2r) + idq
      go to 229
  231 icx = ideks(iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn+nm1r
      idq = (nm2r-1)*ljn+nl2r
  233 icq =  min(icq,idq) + ideks( max(icq,idq))
      go to 229
  230 icx=ideks(ibx)+iax
      icq = (nm2r-1)*lj(nt2r)+nl2r
      go to 235
  223 iax = ideks(nb1r)+nt1r
      if (nt2r.lt.nb2r) go to 236
      ibx = ideks(nt2r)+nb2r
      if (iax-ibx) 237,238,239
  239 icx = ideks(iax)+ibx
      icq = (nm1r-1)*lj(nt1r)+nl1r
      go to 240
  238 icx = ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      go to 241
  237 icx = ideks(ibx) + iax
      icq = (nl2r-1)*lj(nb2r) +nm2r
  245 ljn = lj(nt1r)
      idq = (nm1r-1)*ljn         +nl1r
      icq = (icq-1)*ljn*lj(nb1r) +idq
      go to 229
  236 ibx = ideks(nb2r)+nt2r
      if (iax-ibx) 242,243,244
  244 icx = ideks(iax)+ibx
      icq =(nm1r-1) *lj(nt1r)+nl1r
      ljn = lj(nt2r)
      idq = (nm2r-1)*ljn+nl2r
      icq = (icq-1)*ljn*lj(nb2r)+idq
      go to 229
  243 icx=ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      idq = (nm2r-1)*ljn+nl2r
      go to 233
  242 icx = ideks(ibx)+iax
      icq =(nm2r-1)*lj(nt2r)+nl2r
      go to 245
  222 iax=ideks(nt1r+1)
      iay =  min(nl1r,nm1r) + ideks( max(nl1r,nm1r))
      iby =  min(nl2r,nm2r) + ideks( max(nl2r,nm2r))
      if (nt1r.eq.nt2r) go to 250
      ibx=ideks(nt2r+1)
      if (iax.lt.ibx) go to 251
      icx=ideks(iax)+ibx
      ljn =lj(nt2r)
      icq = (iay-1)*ideks(ljn+1)+iby
      go to 229
  251 icx=ideks(ibx)+iax
      ljn=lj(nt1r)
      icq =(iby-1)*ideks(ljn+1)+iay
      go to 229
  250 icx=ideks(iax+1)
      icq =  min(iay,iby) + ideks( max(iay,iby))
 229  icq=nit(icx)  +  icq
      idud=(icq -1)/igmax
      jto=jto + 1
      kc(jto) = idud  + 1
      icq=icq  -  idud*igmax
      kd(jto)  = icq
      if(jto.lt.iwod)  go to 260
      write(ntype) kc,kd
      jto=0
      jblk=jblk  + 1
260   if(kix.lt.0)  go to 113
      kix=-3
      itr=nb1r
      nb1r=nb2r
      nb2r=itr
      itr=nm1r
      nm1r=nm2r
      nm2r=itr
      go to 220
 114  jz=mt
      do 150 l=1,nr
      jz=jz+1
      nn=iot(jz)
      if (jcon(nn)-1) 151,150,152
  150 continue
  151 jp=l
      if (l.eq.nr) go to 153
      ia=l+1
      do 154 l=ia,nr
      jz=jz+1
      kk=iot(jz)
      if (jcon(kk)-1) 112,154,155
  154 continue
      go to 153
  155 kp=l
      if (l.eq.nr) go to 156
      ia=l+1
      do 157 l=ia,nr
      jz=jz+1
      ii=iot(jz)
      if (jcon(ii).ne.1) go to 112
  157 continue
      go to 156
  152 jp=l
      if (l.eq.nr) go to 158
      ia=l+1
      do 159 l=ia,nr
      jz=jz+1
      kk=iot(jz)
      if (jcon(kk)-1) 160,159,112
  159 continue
      go to 158
  160 kp=l
      if (l.eq.nr) go to 161
      ia=l+1
      do 162 l=ia,nr
      jz=jz+1
      ii=iot(jz)
      if (jcon(ii).ne.1) go to 112
  162 continue
      go to 161
  153 mm=mm+1
      olab(mm)=3
      if (mm.lt.nid) go to 163
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
 163  jz=nt
      do 164 l=1,nr
      jz=jz+1
      if (iot(jz).eq.ll) go to 165
 164  continue
  165 mm=mm+1
      olab(mm)=l
      if (mm.lt.nid) go to 166
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  166 mm=mm+1
      olab(mm)=jp
      if (mm.lt.nid) go to 167
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  167 jz=nz
      if (ip.eq.1) go to 168
      kz=mz
      ja=ip-1
      do 169 l=1,ja
      jz=jz+1
      kz=kz+1
      kk=iot(jz)
      if (kk.ne.iot(kz)) go to 170
  169 continue
  168 if (ip.eq.nd) go to 171
      ja=ip+1
      kz=mz+ip
      do 172 l=ja,nd
      jz=jz+1
      kz=kz+1
      kk=iot(jz)
      if (kk.ne.iot(kz)) go to 170
  172 continue
  171 kk=iot(jz+1)
  170 nt1r=nir(kk)
      nl1r=loc(kk)
      nb1r=nir(ll)
      nm1r=loc(ll)
      nm2r=loc(nn)
      go to 20
  156 mm=mm+1
      olab(mm)=2
      if (mm.lt.nid) go to 173
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  173 jz=nt
      mm=mm+1
      kz=mt
      ig=0
      if (jp.eq.1) go to 174
      ka=1
      ja=jp-1
  176 do 175 l=ka,ja
      jz=jz+1
      kz=kz+1
      if (iot(jz).eq.iot(kz)) go to 175
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.2) go to 177
      go to 800
 175  continue
      go to 174
 800  ka=l
      kz=kz-1
      go to 176
 174  la=kp-1
      if (jp.eq.la) go to 178
      ja=jp+1
      kz=kz+1
 180  do 179 l=ja,la
      kz=kz+1
      jz=jz+1
      if (iot(kz).eq.iot(jz)) go to 179
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.2) go to 177
      go to 801
 179  continue
      go to 178
 801  ja=l
      kz=kz-1
      go to 180
 178  if (kp.eq.nr) go to 181
      ia=kp+1
      kz=mt+kp
 182  do 183 l=ia,nr
      jz=jz+1
      kz=kz+1
      if (iot(jz).eq.iot(kz)) go to 183
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.2) go to 177
      go to 802
 183  continue
      go to 181
 802  ia=l
      kz=kz-1
      go to 182
 181  lab(2)=nr
      if (ig.eq.1) go to 177
      lab(1)=nr1
 177  ia=lab(1)
      jz=ia+nt
      ja=lab(2)
      jj=iot(jz)
      if (jj.eq.ll) go to 184
      olab(mm)=2
 188  if (mm.lt.nid) go to 185
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  185 mm=mm+1
      olab(mm)=ideks(ja-1)+ia
      if (mm.lt.nid) go to 186
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  186 mm=mm+1
      olab(mm)=ideks(kp-1)+jp
      if (mm.lt.nid) go to 187
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
 187  nt1r=nir(nn)
      nt2r=nir(ll)
      nl1r=loc(nn)
      nl2r=loc(ll)
      go to 189
  184 olab(mm)=3
      jz=nt+ja
      jj=iot(jz)
      go to 188
  161 mm=mm+1
      olab(mm)=2
      if (mm.lt.nid) go to 190
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  190 jz=nt
      mm=mm+1
      kz=mt
      ig=0
      if (jp.eq.1) go to 191
      ka=1
      ja=jp-1
 192  do 193 l=ka,ja
      jz=jz+1
      kz=kz+1
      if (iot(jz).eq.iot(kz)) go to 193
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.2) go to 806
      go to 803
 193  continue
      go to 191
 803  ka=l
      kz=kz-1
      go to 192
 191  la=kp-1
      if (jp.eq.la) go to 194
      ja=jp+1
      kz=kz+1
 195  do 196 l=ja,la
      kz=kz+1
      jz=jz+1
      if (iot(kz).eq.iot(jz)) go to 196
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.2) go to 806
      go to 804
196   continue
      go to 194
 804  ja=l
      kz=kz-1
      go to 195
 194  if (kp.eq.nr) go to 197
      ia=kp+1
      kz=mt+kp
 198  do 199 l=ia,nr
      kz=kz+1
      jz=jz+1
      if (iot(jz).eq.iot(kz)) go to 199
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.2) go to 806
      go to 805
  199 continue
      go to 197
  805 ia=l
      kz=kz-1
      go to 198
  197 lab(2)=nr
      if (ig.eq.1) go to 806
      lab(1)=nr-1
  806 ia=lab(1)
      jz=ia+nt
      ja=lab(2)
      jj=iot(jz)
      if (jj.eq.ll) go to 501
      olab(mm)=3
 502  if (mm.lt.nid) go to 503
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  503 mm=mm+1
      olab(mm)=ideks(ja-1)+ia
      if (mm.lt.nid) go to 504
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  504 mm=mm+1
      olab(mm)=ideks(kp-1)+jp
      if (mm.lt.nid) go to 505
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
 505  nt1r=nir(nn)
      nt2r=nir(jj)
      nl1r=loc(nn)
      nl2r=loc(jj)
      kix=3
      nb1r=nir(kk)
      nb2r=nir(ll)
      nm1r=loc(kk)
      nm2r=loc(ll)
      go to 220
  501 olab(mm)=2
      jz=nt+ja
      jj=iot(jz)
      go to 502
 158  mm=mm+1
      olab(mm)=4
      if (mm.lt.nid) go to 506
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
 506  jz=nt
      do 507 l=1,nr
      jz=jz+1
      if (iot(jz).eq.ll) go to 508
 507  continue
 508  mm=mm+1
      olab(mm)=-l
      if (mm.lt.nid) go to 509
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
 509  mm=mm+1
      olab(mm)=jp
      if (mm.lt.nid) go to 510
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
 510  mm=mm+1
      iax=nir(ll)
      iay=ideks(iax+1)
      iaz=ideks(iay+1)
      olab(mm)=iax
      if (mm.lt.nid) go to 893
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
 893  mm=mm+1
      mq=lj(iax)
      mv=ideks(mq+1)
      iaq=iay-iax
      iaw=iaz-iay
      iat=nit(iaz)
      ii=loc(ll)
      kk=loc(nn)
      km=ideks(kk)
      kl=km+kk
      nm=ideks(ii)
      nn=nm+ii
      if (ii.gt.kk) go to 513
      olab(mm)=kk
      if (mm.lt.nid) go to 511
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
 511  mm=mm+1
      olab(mm)=ii
      if (mm.lt.nid) go to 512
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
 512  jb=km+ii
      ib=ideks(jb)+nn+iat
      kb=ideks(kl)+jb+iat
      go to 514
 513  jb=nm+kk
      olab(mm)=ii
      if (mm.lt.nid) go to 551
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
 551  mm=mm+1
      olab(mm)=kk
      if (mm.lt.nid) go to 552
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
 552  ib=ideks(jb)+kl+iat
      kb=ideks(nn)+jb+iat
 514  idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if (jto.lt.iwod) go to 515
      jto=0
      write (ntype) kc,kd
      jblk=jblk+1
 515  idud=(kb-1)/igmax
      kb=kb-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=kb
      nix=1
      if (jto.lt.iwod) go to 516
      jto=0
      write (ntype) kc,kd
      jblk=jblk+1
 516  if (nr.eq.1) go to 517
      jz=mt
      kix=1
      lp=jp
      ir=nr
 518  do 384 l=1,ir
      jz=jz+1
      if (l.eq.lp) go to 384
      mi=iot(jz)
      mar=nir(mi)
      mlr=loc(mi)
      mlb=ideks(mlr)
      mla=mlb+mlr
      if (mar-iax) 385,395,396
 395  if(mla.lt.jb) go to 386
      ib=iat+ideks(mla) +jb
 388  idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto = jto+1
      kc(jto)=idud+1
      kd(jto) =ib
      if(jto.lt.iwod) go to 387
      jto =0
      write(ntype) kc,kd
      jblk=jblk+1
      go to 387
 386  ib=iat+ideks(jb) +mla
      go to 388
 387  if(mlr.lt.ii) go to 389
      kb=mlb+ii
      go to 390
 389  kb=nm+mlr
 390  if(mlr.lt.kk) go to 391
      lb=mlb+kk
      go to 392
 391  lb=mlr+km
 392  ib =  iat + min(kb,lb) + ideks( max(kb,lb))
      idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if(jto.lt.iwod) go to 384
      jto=0
      jblk=jblk+1
      write(ntype) kc,kd
      go to 384
 385  iby=ideks(mar+1)
      iby=iaw+iby
      ibx=iaq+mar
      ibx=ideks(ibx+1)
      ibx=nit(ibx)
      iby=nit(iby)
      ml=lj(mar)
      ib=ideks(ml+1)*(jb-1) +mla+iby
      idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if(jto.lt.iwod) go to 397
      jto=0
      jblk=jblk+1
      write(ntype) kc,kd
 397  kb=(ii-1)*ml +mlr
      lb=(kk-1)*ml +mlr
      ib =  ibx + min(kb,lb) + ideks( max(kb,lb))
      idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto) = ib
      if(jto.lt.iwod) go to 384
      jto=0
      jblk=jblk+1
      write(ntype) kc,kd
      go to 384
396   ibx=ideks(mar+1)
      iby=ideks(ibx) +iay
      iby=nit(iby)
      ibx=ibx-mar+iax
      ibx=ideks(ibx+1)
      ibx=nit(ibx)
      ib=(mla-1)*mv +jb+iby
      idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if(jto.lt.iwod) go to 473
      jto=0
      jblk=jblk+1
      write(ntype) kc,kd
 473  kb=(mlr-1)*mq
      lb=kb+ii
      kb=kb+kk
      ib =  ibx + min(kb,lb) + ideks( max(kb,lb))
      idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if (jto.lt.iwod) go to 384
      jto=0
      jblk=jblk+1
      write (ntype) kc,kd
 384  continue
      if (kix.lt.0) go to 113
      if (nix.lt.0) go to 519
 517  if (nd.eq.1) go to 113
      kix=-1
      jz=mz
      lp=ip
      ir=nd
      go to 518
  519 if (nd.eq.0) go to 113
      kix=-1
      jz=mz
      lp=0
      ir=nd
      go to 518
  106 jz=mt
      do 520 l=1,nr
      jz=jz+1
      ll=iot(jz)
      if (jcon(ll)-1) 521,520,112
  520 continue
  521 ip=l
      if (l.eq.nr) go to 522
      ia=l+1
      do 523 l=ia,nr
      jz=jz+1
      nn=iot(jz)
      if (jcon(nn)-1) 524,523,112
  523 continue
      go to 522
  524 jp=l
      if (l.eq.nr) go to 525
      ia=l+1
      do 526 l=ia,nr
      jz=jz+1
      ii=iot(jz)
      if (jcon(ii).ne.1) go to 112
  526 continue
  525 mm=mm+1
      olab(mm)=2
      if (mm.lt.nid) go to 527
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  527 mm=mm+1
      olab(mm)=1
      if (mm.lt.nid) go to 528
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  528 mm=mm+1
      jz=nt
      kz=mt
      ig=0
      if (ip.eq.1) go to 529
      ka=1
      ja=ip-1
  531 do 530 l=ka,ja
      jz=jz+1
      kz=kz+1
      if (iot(jz).eq.iot(kz)) go to 530
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.2) go to 532
      go to 807
  530 continue
      go to 529
  807 ka=l
      kz=kz-1
      go to 531
  529 la=jp-1
      if (ip.eq.la) go to 533
      ja=ip+1
      kz=kz+1
 534  do 535 l=ja,la
      jz=jz+1
      kz=kz+1
      if (iot(jz).eq.iot(kz)) go to 535
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.2) go to 532
      go to 808
  535 continue
      go to 533
  808 ja=l
      kz=kz-1
      go to 534
  533 if (jp.eq.nr) go to 536
      ia=jp+1
      kz=mt+jp
  537 do 538 l=ia,nr
      kz=kz+1
      jz=jz+1
      if (iot(jz).eq.iot(kz)) go to 538
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.2) go to 532
      go to 809
  538 continue
      go to 536
  809 ia=l
      kz=kz-1
      go to 537
  536 lab(2)=nr
      if (ig.eq.1) go to 532
      lab(1)=nr1
  532 l1=lab(1)
      jz=nt+l1
      kk=iot(jz)
      l2=lab(2)
      jz=nt+l2
      jj=iot(jz)
      olab(mm)=ideks(l2-1)+l1
      if (mm.lt.nid) go to 539
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  539 mm=mm+1
      olab(mm)=ideks(jp-1)+ip
      go to 540
  522 mm=mm+1
      olab(mm)=4
      if (mm.lt.nid) go to 541
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  541 mm=mm+1
      jz=nt
      if (nr.eq.1) go to 546
      kz=mt
      if (ip.eq.1) go to 542
      ja=ip-1
      do 543 l=1,ja
      jz=jz+1
      kz=kz+1
      nn=iot(jz)
      if (nn.ne.iot(kz)) go to 544
  543 continue
      if (ip.eq.nr) go to 546
  542 ja=ip+1
      kz=kz+1
      do 547 l=ja,nr
      jz=jz+1
      kz=kz+1
      nn=iot(jz)
      if (nn.ne.iot(kz)) go to 544
  547 continue
  546 nn=iot(jz+1)
      l1=nr
      olab(mm)=l1
      go to 548
  544 l1=jz-nt
  548 olab(mm)=l1
      if (mm.lt.nid) go to 549
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  549 mm=mm+1
      olab(mm)=ip
      if (mm.lt.nid) go to 550
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  550 mm=mm+1
      iax=nir(ll)
      olab(mm)=iax
      if (mm.lt.nid) go to 897
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
 897  mm=mm+1
      iay=ideks(iax+1)
      iaz=ideks(iay+1)
      mq=lj(iax)
      mv=ideks(mq+1)
      iaq=iay-iax
      iaw=iaz-iay
      iat=nit(iaz)
      ii=loc(ll)
      kk=loc(nn)
      km=ideks(kk)
      nm=ideks(ii)
      if (ii.gt.kk) go to 560
      olab(mm)=kk
      if (mm.lt.nid) go to 561
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  561 mm=mm+1
      olab(mm)=ii
      if (mm.lt.nid) go to 562
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  562 jb=km+ii
      go to 565
 560  olab(mm)=ii
      if (mm.lt.nid) go to 563
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  563 mm=mm+1
      olab(mm)=kk
      if (mm.lt.nid) go to 564
      mm=0
      call pack(olab8,8,olab,2000)
      write (ntape) olab8
  564 jb=nm+kk
  565 nix=-1
      if (nr.eq.1) go to 519
      jz=mt
      kix=1
      lp=ip
      ir=nr
      go to 518
  113 mz=mz+nx
  105 mt=mt+nx
      do 570 k=1,nx
      nt=nt+1
      ll=iot(nt)
  570 jcon(ll)=0
  100 nz=nt+nr
      return
      end
      subroutine eigmrd(a,r,n,mv)
      implicit REAL (a-h,o-z), integer(i-n)
      dimension a(*),r(*)
      if(mv-1) 10,25,10
  10  iq=-n
      do 20 j=1,n
      iq=iq+n
      do 20 i=1,n
      ij=iq+i
      r(ij)=0.0d0
      if(i-j) 20,15,20
  15  r(ij)=1.0d0
  20  continue
  25  anorm=0.0d0
      do 35 i=1,n
      do 35 j=i,n
      if(i-j) 30,35,30
  30  ia=i+(j*j-j)/2
      anorm=anorm+a(ia)*a(ia)
  35  continue
      if(anorm) 165,165,40
  40  anorm=1.414d0*dsqrt(anorm)
      anrmx=anorm*1.0d-14/dfloat(n)
      ind=0
      thr=anorm
  45  thr=thr/dfloat(n)
  50  l=1
  55  m=l+1
  60  mq=((m*m)-m)/2
      lq=((l*l)-l)/2
      lm=l+mq
      if(dabs(a(lm))-thr) 130,65,65
  65  ind=1
      ll=l+lq
      mm=m+mq
      x=0.5d0*(a(ll)-a(mm))
      y=-a(lm)/dsqrt(a(lm)*a(lm)+x*x)
      if(x) 70,75,75
  70  y=-y
  75  if(y .gt. 1.0d0) y=1.0d0
      if(y .lt. -1.0d0) y=-1.0d0
      sinx=y/dsqrt(2.0d0*(1.0d0+(dsqrt(1.0d0-y*y))))
      sinx2=sinx*sinx
      cosx=dsqrt(1.0d0-sinx2)
      cosx2=cosx*cosx
      sincs=sinx*cosx
      ilq=n*(l-1)
      imq=n*(m-1)
      do 125 i=1,n
      iq=(i*i-i)/2
      if(i-l) 80,115,80
  80  if(i-m) 85,115,90
  85  im=i+mq
      go to 95
  90  im=m+iq
  95  if(i-l) 100,105,105
 100  il=i+lq
      go to 110
 105  il=l+iq
 110  x=a(il)*cosx-a(im)*sinx
      a(im)=(a(il)*sinx)+(a(im)*cosx)
      a(il)=x
 115  if(mv-1) 120,125,120
 120  ilr=ilq+i
      imr=imq+i
      x=r(ilr)*cosx-r(imr)*sinx
      r(imr)=(r(ilr)*sinx)+(r(imr)*cosx)
      r(ilr)=x
 125  continue
      x=2.0d0*a(lm)*sincs
      y= (a(ll)*cosx2)+(a(mm)*sinx2)-x
      x=(a(ll)*sinx2)+(a(mm)*cosx2)+x
      a(lm)=(a(ll)-a(mm))*sincs+a(lm)*(cosx2-sinx2)
      a(ll)=y
      a(mm)=x
 130  if(m-n) 135,140,135
 135  m=m+1
      go to 60
 140  if(l-n+1) 145,150,145
 145  l=l+1
      go to 55
 150  if(ind-1) 160,155,160
 155  ind=0
      go to 50
 160  if(thr-anrmx) 165,165,45
 165  iq=-n
      do 185 i=1,n
      iq=iq+n
      ll=i+(i*i-i)/2
      jq=n*(i-2)
      do 185 j=i,n
      jq=jq+n
      mm=j+(j*j-j)/2
      if(a(ll)-a(mm)) 170,185,185
 170  x=a(ll)
      a(ll)=a(mm)
      a(mm)=x
      if(mv-1) 175,185,175
 175  do 180 k=1,n
      ilr=iq+k
      imr=jq+k
      x=r(ilr)
      r(ilr)=r(imr)
 180  r(imr)=x
 185  continue
      return
      end
      subroutine ver_mrdci3(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/mrdci3.m,v $
     +     "/
      data revision /"$Revision: 6176 $"/
      data date /"$Date: 2010-08-10 16:49:47 +0200 (Tue, 10 Aug 2010) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
