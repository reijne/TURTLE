c 
c  $Author: rwa $
c  $Date: 2010-08-10 16:49:47 +0200 (Tue, 10 Aug 2010) $
c  $Locker:  $
c  $Revision: 6176 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/tran4.m,v $
c  $State: Exp $
c  
c     deck=tran4
c ******************************************************
c ******************************************************
c             =   tran4    =
c ******************************************************
c ******************************************************
_IF(parallel)
c***  attempt to adapt for MPP and tools  ***
c***  lines will be flagged with c MPP      ***
c********************************************
_ENDIF
      subroutine adapti(q,otran4)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
INCLUDE(common/scra7)
INCLUDE(common/statis)
      common/craypk/kbcray(1360)
INCLUDE(common/blksiz)
INCLUDE(common/atmblk)
INCLUDE(common/infoa)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/prnprn)
INCLUDE(common/restar)
      common/restri/nfils(63),lda(508),isect(508),ldx(508)
      common/junke/maxt,ires,ipass,nteff,npass1,npass2,lentri,
     *nbuck,mloww,mhi,ntri,iacc
INCLUDE(common/cndx40)
INCLUDE(common/runlab)
INCLUDE(common/timez)
INCLUDE(common/tran)
INCLUDE(common/blocko)
INCLUDE(common/files)
INCLUDE(common/discc)
INCLUDE(common/trntim)
INCLUDE(common/harmon)
INCLUDE(common/newmrd_sort)
INCLUDE(common/fpinfo)
INCLUDE(common/restrj)
      character*10 charwall
c
      character*7 fnm
      character*6 snm
      data fnm,snm/"tran4.m","adapti"/
      dimension q(*)
      data zcas/'mcscf'/
c
      call cpuwal(begin,ebegin)
      nav = lenwrd()
      nmaxim=maxorb-1
      iscftp=0
      iacc=iacc4
c
      npass1=npas41
      npass2=npas42
      oprin4(10)=otran4
      oindx4 = otran4
c
c ---- print flag settings
c
      do 87654 i=1,7
      oprin4(i)=.true.
      if(nprint.eq.3)oprin4(i)=.false.
87654 continue
      oprin4(8)=.true.
      oprin4(9)=.true.
c
_IF(parallel)
c     sort file settings
c
      call setsrtp(nsz,o255i,
     +  nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     +  nsz510,nsz341,nsz342,nsz680)
      nsstat=0
_ENDIF
      iblkq4=ibl7la
      nbas4=num
      ndump4=idaf
      newb4=newbas0
      if(nbas4.le.1.or.nbas4.gt.nmaxim)call caserr(
     *'error in 4-index preprocessor')
      nbb4=nbas4
      nm=nmaxim+1
      do i=1,nm
      ilifq(i)=(i-1)*nbas4
      if (o255i) then
       i4096(i)=i*256
      else
       i4096(i)=i*65536
      endif
      iky(i)=((i-1)*i)/2
      enddo
      lenb4=iky(nbas4+1)
      m12=12/nav
      if(otran4)go to 7000
c
      if(nprint.ne.-5)write(iwr,6800)
 6800 format(/1x,104('=')//
     *40x,26('*')/40x,'integral symmetry adaption'/
     *40x,26('*')/)
      mrest=2
      itime=8
      idsec=isect(483)
      nsa4=newbas0
      ncol4=newbas0
      nblkq4=nbas4*ncol4
      do 5009 i=1,nm
      mapie(i)=i
 5009 ilifm(i)=ilifq(i)
      lena4=iky(nsa4+1)
c
      call symat(q)
c
      go to 5078
c
 7000 if(nprint.ne.-5)write(iwr,6801)
 6801 format(/1x,104('=')//
     *40x,23('*')/
     *40x,'integral transformation'/
     *40x,23('*')/)
      idsec=isect(486)
      mrest=8
      itime=11
      if(zscftp.eq.zcas) iscftp = 1
c
      if (opunch(12))
_IFN1(iv)     & open(ipu,file='moints.ascii',form='formatted',status='unknown')
_IF1(iv)     & call caserr('punched moints file not available')
c
      call trani(q)
      if (odontgened6) goto 6010
c
c ----- restart adaption/transformation
c
 5078 if(irest.lt.mrest)go to 5077
c
       call secget(idsec,500,isecbl)
       call readi(master,m12*nav,isecbl,ndump4)
       if(nprint.ne.-5)write(iwr,5076)idsec,ibl3d,yed(ndump4)
5076   format(1x,'4-index dump information restored from section',
     *i4,' of dumpfile starting at block',i6,' on ',a4)
       goto 5075
 5077 nfiles=1
      nfilef=1
      indxi=1
      indxj=1
      master=0
c     secondary mainfile
      junits=n4tape(1)
      jblkas=n4blk(1)
      jblkrs=jblkas-n4last(1)
c     final mainfile
c
      if(otran4)go to 100
      junitf=n1tape(1)
      jblkaf=n1blk(1)
      jblkrf=jblkaf-n1last(1)
      go to 200
 100  junitf=n6tape(1)
      jblkaf=n6blk(1)
      jblkrf=jblkaf-n6last(1)
 200  call secput(idsec,500,1,isecbl)
      if(nprint.ne.-5)write(iwr,5074)idsec,ibl3d,yed(ndump4)
5074  format(' 4-index information dumped to section',
     *i4,' of dumpfile starting at block',i6,' on ',a4)
      call revind
5075  if(nprint.eq.-5)go to 500
      write(iwr,32)
32    format(/' 2-electron integral files'/1x,25('*'))
      if(iscftp.eq.0)
     * call filprn(n2file,n2blk,n2last,n2tape)
      if(iscftp.eq.1)
     * call filprn(n1file,n1blk,n1last,n1tape)
      write(iwr,33)
33    format(/' secondary mainfile'/1x,18('*'))
      call filprn(n4file,n4blk,n4last,n4tape)
      if(otran4) go to 600
       write(iwr,30)
30     format(/' symmetry adapted mainfile'/1x,26('*'))
       call filprn(n1file,n1blk,n1last,n1tape)
      go to 500
 600  write(iwr,40)
 40   format(/' transformed integral files'/1x,26('*'))
      call filprn(n6file,n6blk,n6last,n6tape)
 500  continue
      icc = igmem_alloc_inf((nsa4+nbas4)*nbb4,fnm,snm,"icc",
     +                      IGMEM_NORMAL)
      call vecget(q,q(icc),otran4)
      call setsto(1360,0,kbcray)
      if (opass11) then
         write(iwr,*) ' '
         write(iwr,*) ' **** bypass 4-index transformation ****'
      else
         call adapt4(q,q(icc),otran4)
      end if
      if(tim.lt.timlim)go to 6000
      irest=mrest
      go to 6010
6000  call adapt2(q,q(icc),otran4)
      call gmem_free_inf(icc,fnm,snm,'icc')
      irest=0
c
c ----- revise section idsec for pseudo restart
c
      nfiles=1
      nfilef=1
      indxi=1
      indxj=1
      master=0
c     secondary mainfile
      junits=n4tape(1)
      jblkas=n4blk(1)
      jblkrs=jblkas-n4last(1)
c     final mainfile
c
      if(otran4)go to 300
      junitf=n1tape(1)
      jblkaf=n1blk(1)
      jblkrf=jblkaf-n1last(1)
      go to 400
 300  junitf=n6tape(1)
      jblkaf=n6blk(1)
      jblkrf=jblkaf-n6last(1)
 400  call wrt3i(master,m12*nav,isecbl,ndump4)
      top=cpulft(1)
      if(nprint.ne.-5)write(iwr,71845)top,charwall()
71845 format(/
     *' end of integral transformation at ',
     *f8.2,' seconds',a10,' wall',/)
c
      call revind
      call clredx
c
c rwah mrdci sort of integrals
c
6010  if (osortmr) then
c     first delete sortfile to minimise space
c     since sortdrver creates its own sortfile
      if (.not.ocifp) call closbf(0)
      call sortdriver(q)
c
c     delete transformed integral file if not
c     to be saved
      call delfil(m6tape(1))
c
      endif
c
      call timana(itime)
c
      return
      end
      subroutine trani(q)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
c
      dimension q(*)
INCLUDE(common/machin)
INCLUDE(common/prints)
INCLUDE(common/iofile)
INCLUDE(common/cndx40)
INCLUDE(common/tran)
INCLUDE(common/blocko)
INCLUDE(common/restar)
INCLUDE(common/harmon)
      common/craypk/nact,nacta(maxorb),nprina,
     *              mact,macta(maxorb),mprina, ihamp,iseca
      common/lsort/iijjnk(4*maxorb),
     *             pop(maxorb),pot(2),nspa(3),
     * ncore,mapcie(maxorb),icorec(maxorb),nvalsp,ixspa(maxorb,2),
     * iqsec,jdsec
      common/atmol3/mina(2),mouta,moutb
      common/junkc/zjob,zdate,ztime,zprog,ztype,zspace(14),ztext(10)
      common/junk/
     * evalue(maxorb),eocc(maxorb+1),
     * nbas,newb,ncoll,ivalue,iocc
      common/multic/radius(40),irad(57+mcfzc),itype(maxorb),isecn
c
      character*7 fnm
      character*5 snm
      data fnm,snm/"tran4.m","trani"/
c
      data m29,isecb/29,470/
      isecv=mouta
      if(moutb.ne.0)isecv=moutb
      if(iscftp.eq.1.and.isecn.gt.0)isecv=isecn
c
c     restore vectors
c
      call secget(isecv,3,iblk)
      call rdchr(zjob,m29,iblk,ndump4)
      call reads(evalue,mach(8),ndump4)
      len2=lensec(mach(9))
      iblkv=iblk+len2+lensec(mach(8))+1
      if(nbas4.ne.nbas) call caserr(
     *'restored eigenvectors are in incorrect format')
      ncol4=ncoll
      if(nprint.ne.-5)
     +   write(iwr,89004)isecv,ztype,ztime,zdate,zjob,ztext
89004 format(/' vectors restored from section ',i4//
     *' header block information : '/
     * 1x,a7,'vectors created at ',a8,' on ',a8,
     *' in the job ',a8/
     *' with the title: ',10a8)
      nav = lenwrd()
      if(iscftp.eq.1)go to 4000
      call readis(ilifc,mach(9)*nav,ndump4)
      if(.not.oprint(46). or .nprint.eq.-5) go to 4000
      write(iwr,98)
 98   format(/14x,'list of lcbf'/14x,12('*')//
     *' coefficient old orbital nterm new orbital'/1x,41('-')/)
c
      do 97 i=1,newb4
      n=ntran(i)
      j=ilifc(i)
      write(iwr,5038) ctran(j+1),itran(j+1),n,i
 5038 format(f12.7,i12,i6,i12)
      if(n.eq.1)go to 97
      do 96 k=2,n
 96   write(iwr,5038)ctran(j+k),itran(j+k)
 97   continue
      write(iwr,5037)newb4
 5037 format(/' no. of lcbf = ',i4)
c
4000  iqsec=isecv
      call dcopy(maxorb,eocc,1,pop,1)
c
c     now restore vectors and copy to ed7
c
      nblkq4=ncol4*nbas4
c...   for harmonic the rdedx might read too much (i.e. nbas4*nbas4)
c     i10 = igmem_alloc_inf(nblkq4,fnm,snm,'i10',IGMEM_NORMAL)
      i10 = igmem_alloc_inf(nbas4*nbas4,fnm,snm,'i10',IGMEM_NORMAL)
      call rdedx(q(i10),nbas4*nbas4,iblkv,ndump4)
c
      if(iscftp.eq.1) then
         call comharm(q(i10),'vectors',ilifq)
         ncol4 =newbas0
         nbas4 = newbas0
         nblkq4=ncol4*nbas4
      end if
c
      call wrt3(q(i10),nblkq4,iblkq4,num8)
c
      call gmem_free_inf(i10,fnm,snm,'i10')
c
c     now restore =active= and =core= specifications
c
      call secget(isecb,1005,iblka)
      call readi(nact,mach(14)*nav,iblka,ndump4)
c
c     active list
c
      nsa4=nact
      if(nsa4-1)1000,1001,1002
 1000 call caserr(
     *'invalid number of active orbitals')
 1002 if(nsa4.gt.ncol4)go to 1000
      do 1003 i=2,nsa4
      if((locat1(nacta,i-1,nacta(i))).ne.0)
     * call caserr(
     *'invalid orbital specified in active list')
 1003 continue
      do 1004 i=1,nsa4
      j=nacta(i)
      mapie(i)=j
 1004 ilifm(i)=ilifq(j)
 1001 if(nprint.ne.-5) then
       write(iwr,1005)
 1005  format(/
     * ' the following orbitals included in active list'/
     * ' =============================================='/
     * 4(' function e/i label',4x)/1x,88('-'))
       write(iwr,1006)(mapie(k),k,k=1,nsa4)
 1006  format(4(i9,i10,4x))
      endif
      lena4=iky(nsa4+1)
      if(nprina.eq.1)oprin4(8)=.false.
c
c     frozen list
c
      ncore=mact
      if(ncore-1)2000,2004,2002
 2000 if(nprint.ne.-5) write(iwr,2003)
 2003 format(/
     *' no functions specified in frozen core list')
      go to 3000
 2002 if(ncore.gt.ncol4)call caserr(
     *'invalid number of functions in core list')
 2004 do 2005 i=1,ncore
      j=macta(i)
      mapcie(i)=j
 2005 icorec(i)=ilifq(j)
      if(ncore.le.1)go to 2007
      do 2006 i=2,ncore
      if((locat1(mapcie,i-1,mapcie(i))).ne.0) call caserr(
     *'invalid function specified in frozen core list')
 2006 continue
 2007 if(nprint.ne.-5) then 
       write(iwr,2008)
 2008  format(/
     * ' the following mos included in frozen core list'/
     * ' ----------------------------------------------'/
     * 4(' function e/i label',4x))
       write(iwr,1006)(mapcie(k),k,k=1,ncore)
      endif
c
 3000 if(mprina.eq.1)oprin4(9)=.false.
      jdsec=iseca
      if(jdsec.gt.350.and.jdsec.ne.466)
     * call caserr(
     *'invalid section specified on onelec/core directive')
      if(jdsec.lt.1)jdsec=466
      if(nprint.ne.-5) write(iwr,3001)jdsec
3001  format(/
     *' route transformed core hamiltonian to section ',i4)
c
      return
      end
      subroutine adaps1(g,nijkl)
c...   sorts mainfile onto lfn sort --- so that for a
c...   given rs comb. al pq combs. available
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
_IF(parallel)
INCLUDE(common/sizes)
_ENDIF
c
_IF1(iv)      integer * 2 nijkl
      dimension g(*)
      dimension nijkl(*)
c
      common/stak/btri,mlow,nstack,iblock,mstack
INCLUDE(common/blksiz)
INCLUDE(common/atmblk)
      common/junke/maxt,ires,ipass,nteff,npass1,npass2,lentri,
     *nbuck,mloww,mhi,ntri
INCLUDE(common/cndx40)
      common/junk/nwbuck(1500),
     * itx(3400),ktx(3400),gtx(3400)
      common/bufb/nwbnwb,lnklnk,gout(5119)
      common/three/mark(1500),ibase(1500),ibasen(1500)
INCLUDE(common/filel)
      common/blkin/gin(1)
INCLUDE(common/iofile)
_IF(parallel)
INCLUDE(common/prints)
      common/nodinf/mpid,minode,mihost,ldim,nnodes,nodscr
     * ,maxnod(maxlfn)
_ENDIF
_IF(linux)
      external fget
_ENDIF
c
      data maxb/9999999/
_IF(parallel)
      call closbf(0)
      call setbfa(-1)
_ENDIF
c...   determine base and limit triangles for this pass
      mloww=master
      mhi=master+nteff
      mlow=mloww+1
      if(mhi.gt.lenb4)mhi=lenb4
      mtri=mhi-mlow
c...   determine minimum no. of bucks.
      nbuck=ires
400   ntri=mtri/nbuck
      if(ntri.ge.maxt)goto 401
_IF(parallel)
      nbuck=nbuck-nnodes
      if (nbuck.lt.0) nbuck=0
_ELSE
      nbuck=nbuck-1
_ENDIF
      if(nbuck)400,401,400
_IF(parallel)
401   nbuck=nbuck+nnodes
_ELSE
401   nbuck=nbuck+1
_ENDIF
      ntri=mtri/nbuck+1
_IF(parallel)
      if(oprint(58)) then
      write(6,33333) minode,nbuck, ires,ntri, mtri
33333 format(1x,'adaps1:minode,nbuck,ires,ntri,mtri = ',5i10)
      endif
_ENDIF
c...   ntri=max. no. of triangles controlled by 1 bucket
c...   nbuck=number of buckets
c
_IF(parallel)
      if(nbuck.gt.ires)then
      write(6,33343)minode,nbuck,ires
33343 format(' minode,nbuck,ires = ',2i5)
       call caserr(
     *'invalid number of buckets in adaps1')
      endif
_ENDIF
      btri=ntri
      btri=1.000000001d0/btri
      call setsto(nbuck,maxb,mark)
      call setsto(nbuck,0,nwbuck)
      nstack=0
      nwbnwb=nsz340
      iblock=0
c...   start loop over mainfile blocks
      do 20113 ifile=1,lfile
      lbl=llblk(ifile)
      if(lbl)10000,20113,10000
10000 iunit=lotape(ifile)
      call search(liblk(ifile),iunit)
20110 call fget(gin,m,iunit)
      if(m)20121,20113,20121
20121 nnn=isort1(itx(nstack+1),ktx(nstack+1),gtx(nstack+1))
      nstack=nstack+nnn
      if(nstack.ge.2721) then
       call stackr(g,nijkl)
      endif
      lbl=lbl+1
      if(lbl)20110,20113,20110
20113 continue
c
c...   mainfile now swept
c...   clear up output
c
      if(nstack.ne.0) then
        call stackr(g,nijkl)
      endif
c
      do 7 ibuck=1,nbuck
      nwb=nwbuck(ibuck)
      if(nwb)999,7,999
999   ib=ibase(ibuck)
      ibn=ibasen(ibuck)
      call stopbk
      nwbnwb=nwb
      lnklnk=mark(ibuck)
      call dcopy(nwb,g(ib+1),1,gout,1)
_IFN1(iv)      call pack(gout(nsz341),lab1632,nijkl(ibn+1),nsz680)
_IF1(iv)      call fmove(nijkl(ibn+1),gout(nsz341),nsz170)
      call sttout
      mark(ibuck)=iblock
      iblock=iblock+nsz
7     continue
      return
      end
      subroutine adapc1(q,r,v,qq)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
c
      dimension r(*),v(*)
       dimension qq(*),q(*)
c...   processes sorted mainfile on lfn sort
c...   to produce integrals of the form (i j/r s)
c...   on secondary mainfile
      common/junke/maxt,ires,ipass,nteff,npass1,npass2,lentri,
     *nbuck,mloww,mhi,ntri,iacc
INCLUDE(common/restar)
INCLUDE(common/cndx40)
      common/bufb/nkk,mkk,gin(5119)
      common/blkin/gout(510),mword
INCLUDE(common/iofile)
INCLUDE(common/mlngth)
      common/three/mark(1500),ibase(1500),ibasen(1500)
INCLUDE(common/blocko)
_IF1(iv)      common/craypk/kbcray(340),macray(340)
_IFN1(iv)       common/craypk/kbcray(680)
INCLUDE(common/atmblk)
c
_IF(parallel)
      common/nodinf/mpid,minode,mihost,ldim,nnodes,nodscr,
     * maxnod(maxlfn)
INCLUDE(common/prints)
_ENDIF
      character*7 fnm
      character*6 snm
      data fnm,snm/"tran4.m","adapc1"/
c
      data maxb/9999999/
c
_IF(parallel)
c MPP
c MPP  all sorted // exchange  #'s of backchains
c MPP
      nbpn = (nbuck-1)/nnodes + 1
      i00 = igmem_alloc_inf((nbpn*nnodes+1)/lenwrd(),fnm,snm,
     &                      'bucket-addr',IGMEM_DEBUG)
_IF(ipsc)
      call excbuc(mark,nbuck,qq(i00),r,nnodes)
_ELSE
      call excbuc(mark,nbuck,qq(i00),r,r(nbpn+1),nnodes)
_ENDIF
      call gmem_free_inf(i00,fnm,snm,'bucket-addr')
c MPP
_ENDIF
      mword=0
_IFN1(iv)      int4=-1
      small=10.0d0**(-iacc-2)
      call stopbk
_IF(parallel)
c MPP
      call closbf(0)
      mloww = mloww + minode*ntri
      master = master + minode*ntri
      do 1 ibpn=1,nbpn
c MPP
_ELSE
      do 1 ibuck=1,nbuck
_ENDIF
      mhigh=mloww+ntri
      if(mhigh.gt.mhi)mhigh=mhi
_IF(parallel)
c MPP
      if(mloww.ge.mhi) go to 1
c MPP
_ENDIF
      mtri= mhigh-mloww
      mloww=mloww+1
      nsize=mtri*lenb4
      i10 = igmem_alloc_inf(nsize,fnm,snm,'i10',IGMEM_NORMAL)
      call vclr(qq(i10),1,nsize)
c...   read in a core load
_IFN(parallel)
      mkk=mark(ibuck)
      go to 999
 888  iblock=mkk
      call rdbak(iblock)
      call stopbk
c...   move block to processing area
      call transc(qq(i10))
999   if(mkk.ne.maxb)go to 888
_ENDIF
_IF(parallel)
c MPP
      do 2 inode=1,nnodes
         ibuck = (ibpn-1)*nnodes+inode
         mkk=mark(ibuck)
         if(oprint(58)) then
         write(6,33335) minode,inode,ibuck,mkk
33335    format(' adapc1, minode,inode,ibuck,mkk = ',4i6)
         endif
         call setbfa(inode-1)
         go to 999
 888     iblock=mkk
         call rdbak(iblock)
         call stopbk
c...   move block to processing area
         call transc(qq(i10))
999      if(mkk.ne.maxb)go to 888
2     call closbf(0)
c MPP
_ENDIF
      map=1
_IF(parallel)
       if(oprint(58)) then
       write(6,33333) minode,map,master,mtri,nsa4,nbas4,lena4
33333 format(
     * ' adapc1: minode,map,master,mtri,nsa4,nbas4,lena4 ',7i7)
c      nnnn=nbas4*nbas4
c      write(6,33334)(q(loop),loop=1,nnnn)
c33334 format(1x,8f13.7)
c      write(6,33334)(q(loop+nanb),loop=1,nnnn)
       endif
_ENDIF
      do 777 itri=1,mtri
      master=master+1
c...   master=iky(r)+s
      call mult3(q,qq(i10+map-1),r,v,nsa4,nbas4)
c...   compute integrals of the form (i j/r s)
      do 1009 kb=1,lena4
      top=r(kb)
c...    (i j/r s) now in top
      if( dabs(top).lt.small)goto 1009
      mword=mword+1
      gout(mword)=top
_IF1(iv)      kbcray(mword)=kb
_IF1(iv)      macray(mword)=master
_IFN1(iv)      int4=int4+2
_IFN1(iv)      kbcray(int4  )=kb
_IFN1(iv)      kbcray(int4+1)=master
      if(mword.eq.nintmx)then
       call blocka
_IFN1(iv)      int4=-1
      endif
1009  continue
777   map=map+lenb4
      call gmem_free_inf(i10,fnm,snm,'i10')
_IF(parallel)
       mloww = mhigh + (nnodes-1)*ntri
       master = master + (nnodes-1)*ntri
 1     continue
_ELSE
1      mloww=mhigh
_ENDIF
      if(mword.ge.1)call blocka
_IF(parallel)
c MPP
c MPP   set mloww correctly for  multi-pass
c MPP
      mloww = mloww - minode*ntri
      master = master - minode*ntri
c MPP
c MPP  synchronize
c MPP
      call pg_synch(3333)
c MPP
_ENDIF
      return
      end
      subroutine blocka
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/files)
INCLUDE(common/cndx40)
INCLUDE(common/atmblk)
      common/blkin/gout(510),mword
      common/craypk/integ(680)
      common/disc/isel(5),ipos(maxlfn)
INCLUDE(common/restar)
_IFN1(civ)      common/prints/oprint(60),ciprnt
c
      if(ipos(junits).ne.jblkas)callsearch(jblkas,junits)

_IF(ibm,vax)
      call pak4v(integ,gout(num2ep+1))
_ELSE
      call pack(gout(num2ep+1),lab1632,integ,numlabp)
_ENDIF
_IF(ipsc)
      if (oprint(58))call pr2es(junits,jblkas)
_ELSEIFN(cray,ibm,vax,t3d)
      if (oprint(58)) then
         call dchksm(mword,gout,1,'blocka')
         call pchksm(mword,gout(num2ep+1),1,'blocka')
         call pr2es(junits,jblkas)
      endif
_ENDIF
      call put(gout(1),511,junits)
      mword=0
      jblkas=jblkas+1
      jblkrs=jblkrs+1
      if(jblkrs)41,40,41
c...   change channel
 40    m4file=nfiles
       m4tape(m4file)=n4tape(nfiles)
       m4blk(m4file)=n4blk(nfiles)
       m4last(m4file)=jblkas
       nfiles=nfiles+1
      if(nfiles.gt.n4file)call caserr(
     *'insufficient no. of secondary mainfile sections allocated')
      junits=n4tape(nfiles)
      jblkas=n4blk(nfiles)
      jblkrs=jblkas-n4last(nfiles)
  41  return
      end
      subroutine blockb
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/files)
INCLUDE(common/iofile)
INCLUDE(common/prnprn)
INCLUDE(common/cndx40)
INCLUDE(common/atmblk)
      common/blkin/gout(510),mword
INCLUDE(common/restar)
      common/disc/isel(5),ipos(maxlfn)
_IFN1(civ)      common/prints/oprint(60),ciprnt
      common/craypk/integ(680)
c
      if(ipos(junitf).ne.jblkaf)call search(jblkaf,junitf)
c
      if (oindx4.and.opunch(12)) then
      if(o255i) then
        ifact = 256
      else
        ifact = 65536
      endif
        do 1011 iq = 1,mword
           ii = integ(2*iq-1)/ifact
           jj = integ(2*iq-1) - ifact*ii
           kk = integ(2*iq)/ifact
           ll = integ(2*iq) - ifact*kk
           write(ipu,*) gout(iq),ii,jj,kk,ll
1011    continue
      endif
c
_IF(ibm,vax)
      call pak4v(integ,gout(num2ep+1))
_ELSE
      call pack(gout(num2ep+1),lab1632,integ,numlabp)
_ENDIF

_IF(ipsc)
      if (oprint(58))call pr2ef(junitf,jblkaf)
_ELSEIFN(cray,ibm,vax,t3d)
      if (oprint(58)) then
         call dchksm(mword,gout,1,'blockb')
         call pchksm(mword,gout(num2ep+1),1,'blockb')
         call pr2ef(junitf,jblkaf)
      endif
_ENDIF
      call put(gout(1),511,junitf)
      mword=0
      jblkaf=jblkaf+1
      jblkrf=jblkrf+1
      if(jblkrf)41,40,41
 40   if(.not.oindx4)go to 50
c ... change channel
      m6file=nfilef
      m6tape(m6file)=n6tape(nfilef)
      m6blk (m6file)=n6blk  (nfilef)
      m6last(m6file)=jblkaf
      nfilef=nfilef+1
      if(nfilef.le.n6file)go to 60
 70   call caserr(
     *'insufficient no. of final mainfile sections allocated')
 60   junitf=n6tape(nfilef)
      jblkaf=n6blk(nfilef)
      jblkrf=jblkaf-n6last(nfilef)
      go to 41
 50   m1file=nfilef
      m1tape(m1file)=n1tape(nfilef)
      m1blk (m1file)=n1blk(nfilef)
      m1last(m1file)=jblkaf
      nfilef=nfilef+1
      if(nfilef.gt.n1file) go to 70
      junitf=n1tape(nfilef)
      jblkaf=n1blk(nfilef)
      jblkrf=jblkaf-n1last(nfilef)
41    return
      end
      subroutine vecget(q,cc,otran4)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c
INCLUDE(common/tran)
       common/junk/ilifx(maxorb)
INCLUDE(common/cndx40)
INCLUDE(common/prints)
INCLUDE(common/blocko)
INCLUDE(common/mlngth)
INCLUDE(common/iofile)
INCLUDE(common/zorac)
      common/lsort/iijjnk(4*maxorb),
     +             commen(maxorb),pot(2),nspa(3),ncore,
     +             mapcie(maxorb),icorec(maxorb)
c
       dimension q(*),cc(*)
c
c      ------ get core memory
c
      nanb=nsa4*nbas4
      nanb2=nanb+nanb
      nana=nsa4*nsa4
      nananb=nana+nanb
      nbasd2=nbas4*0.25d0
c checkout
c     nbasd2=nbas4*0.01
      nbnb=nbb4*nbas4
c
      i10 = 1
      i20=i10+nanb
c
      call rdedx(cc(i20),nblkq4,iblkq4,num8)
c
      if(otran4) then
       if(iscftp.eq.1)otran=.true.
       if(oprin4(9).or.ncore.eq.0) go to 600
       call tdown(cc(i10),ilifq,cc(i20),icorec,ncore)
       write(iwr,610)
 610   format(//1x,104('-')/
     * 40x,18('*')/40x,'core eigen vectors'/40x,18('*'))
       call prsql(cc(i10),ncore,nbas4,nbas4)
 600   call tdown(cc(i10),ilifq,cc(i20),ilifm,nsa4)
       if(iscftp.eq.1)otran=.false.
       if(.not.oprin4(8))then
        write(iwr,630)
 630    format(/1x,104('-')/40x,20('*')/
     *                      40x,'active eigen vectors'/
     *                      40x,20('*'))
        call prsql(cc(i10),nsa4,nbas4,nbas4)
       endif
      else
       do 777 i=1,nsa4
777    ilifx(i)=(i-1)*nbb4
       otran=.true.
       call tdown(cc(i10),ilifx,cc(i20),ilifm,nsa4)
       otran=.false.
      endif
c
      looper=i10
      do 99 loop=1,nanb
      if( dabs(cc(looper)).lt.1.0d-11)cc(looper)=0.0d0
 99   looper=looper+1
c
      if(.not.otran4) then
       if(oprint(46))then
        write(iwr,620)
 620    format(/1x,104('-')/40x,16('*')/
     *                      40x,'s.a.b.f. vectors'/
     *                      40x,16('*'))
        call prsql(cc(i10),nsa4,nbas4,nbas4)
       endif
       call wrt3(cc(i10),nanb,iblkq4,num8)
      endif
c
_IF1(civ)      call dagger(nbas4,nsa4,cc(i10),nbas4,cc(i20),nsa4)
_IFN1(civ)      call rmtran(cc(i10),nbas4,cc(i20),nsa4,nbas4,nsa4)
c
       return
      end
      subroutine adapt4(q,cc,otran4)
c...
c...   control routine for the integral symmetry adaption
c...   and integral transformation ..
c...   adaps1,adapc1,adapc2,adapc2
c...
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
c
_IF(parallel)
      common/nodinf/mpid,minode,mihost,ldim,nnodes,nodscr
     * ,maxnod(maxlfn)
_ENDIF
      common/craypk/ijkl(1360)
INCLUDE(common/cndx40)
INCLUDE(common/blksiz)
INCLUDE(common/atmblk)
      common/junke/maxt,ires,ipass,nteff,npass1,npass2,lentri,
     *nbuck,mloww,mhi,ntri,iacc
INCLUDE(common/blocko)
INCLUDE(common/iofile)
      common/blkin/gout(511)
INCLUDE(common/restar)
INCLUDE(common/filel)
      common/three/mark(1500),ibase(1500),ibasen(1500)
INCLUDE(common/trntim)
INCLUDE(common/timez)
INCLUDE(common/mlngth)
_IF(parallel)
INCLUDE(common/prints)
_ENDIF
      character*10 charwall
c
      character*7 fnm
      character*6 snm
      data fnm,snm/"tran4.m","adapt4"/
      dimension q(*)
c
c     cc holds the wavefunction coefficients as read in vecget.
c
      dimension cc(*) 
      data m1,m2/1,2/
c
      slder=cpulft(1)
      if(nprint.eq.-5)go to 403
      if(.not.otran4 )write(iwr,401)slder,charwall()
      if(otran4) write(iwr,402)slder,charwall()
 401  format(/
     *' start of symmetry adaption at',f8.2,' seconds',a10,' wall')
 402  format(/
     *' start of integral transformation at',f8.2,' seconds',
     * a10,' wall')
 403  iv=lena4
      niqq=iv+nananb
c
c     look up how much memory is available
c
      lword4 = igmem_max_memory()-4*igmem_overhead()
      left=lword4-niqq
c
      if(left.lt.1)call caserr(
     *'insufficient memory for integral transformation')
      call ibasgn(1500,0,nsz340,ibase)
      call ibasgn(1500,0,nsz680,ibasen)
      maxt=left/lenb4
      nword=lword4-nanb2
_IF(ibm,vax)
      nword=(nword*2)/3
      nword=(nword/4)*4
      nword2=nword/4
_ELSE
      if(o255i) then
       nword=nword/3
       nword=(nword/2)*2
       nword2=nword
      else
       nword=nword/2
       nword=(nword/2)*2
       nword2=nword
      endif
_ENDIF
      ires=nword/nsz340
_IF(parallel)
      ires=(ires/nnodes)*nnodes
_ENDIF
      mword=nword
c
      i10=igmem_null()
      i20=igmem_null()
      i30=igmem_null()
      i40=igmem_null()
      i41=igmem_null()
      i50=igmem_null()
c
      if(ires.lt.1.or.maxt.lt.1) then
       call caserr(
     +        'insufficient memory for integral transformation')
      endif
      if(ires.gt.1500)ires=1500
       if(master.eq.lenb4)goto 7766
c...   determine min. no. of passes for adaps1/adapc1
      i=lenb4-master-1
102   nteff=i/npass1+1
      if(((nteff-1)/ires).lt.maxt)goto 103
      npass1=npass1+1
      goto 102
103   if(npass1.gt.i)npass1=i+1
      if(nprint.ne.-5)write(iwr,204)m1,npass1
      if(otran4.and.iscftp.eq.1) then
       lfile=m1file
       do 2020 i=1,lfile
       lotape(i)=m1tape(i)
       liblk(i) =m1blk(i)
2020   llblk(i) =liblk(i)-m1last(i)
      else
       lfile=m2file
       do 1000 i=1,lfile
       lotape(i)=m2tape(i)
       liblk(i)=m2blk(i)
 1000  llblk(i)=liblk(i)-m2last(i)
      endif
      do 1 ipass=1,npass1
         call setsto(1360,0,ijkl)
         i10=igmem_alloc_inf(mword,fnm,snm,'i10',IGMEM_DEBUG)
         i40=igmem_alloc_inf(nword2,fnm,snm,'i40',IGMEM_NORMAL)
         call adaps1(q(i10),q(i40))
         call gmem_free_inf(i40,fnm,snm,'i40')
         call gmem_free_inf(i10,fnm,snm,'i10')
         top=cpulft(1)
         if(nprint.ne.-5)write(iwr,2222)m1,ipass,top,charwall()
2222     format(/' end of sort',i1,' pass',i6,' at ',f8.2,' seconds',
     1           a10,' wall')
_IF(parallel)
c MPP
         call pg_synch(3333)
c MPP
_ENDIF
         i10=igmem_alloc_inf(iv,fnm,snm,'i10',IGMEM_DEBUG)
         i20=igmem_alloc_inf(nananb,fnm,snm,'i20',IGMEM_NORMAL)
         call adapc1(cc,q(i10),q(i20),q)
         call gmem_free_inf(i20,fnm,snm,'i20')
         call gmem_free_inf(i10,fnm,snm,'i10')
         call dumpra(1,ipass,npass1,jblkas,junits)
         if(tim.gt.timlim) go to 2223
 1    continue
7766  if(indxi.gt.nsa4)goto 7755
      lfile=m4file
       do 9998 i=1,lfile
       lotape(i)=m4tape(i)
       liblk(i)=m4blk(i)
 9998  llblk(i)=liblk(i)-m4last(i)
      if(nprint.ne.-5)write(iwr,4007)
4007  format(/' status of secondary mainfile'/
     *        ' ****************************')
      if(nprint.ne.-5)call filprn(m4file,m4blk,m4last,m4tape)
_IFN(parallel)
c ..
c ..  to minimize sortfile space, delete and re-open
c ..
      call closbf(0)
      call setbfa
c
_ENDIF
      lentri=iky(nsa4+1)
c...   determine min. no. of passes for adaps2/adapc2
      i=lentri-iky(indxi)-indxj
202   nteff=i/npass2+1
      if(((nteff-1)/ires).lt.maxt)goto 203
      npass2=npass2+1
      goto 202
203   if(npass2.gt.i)npass2=i+1
      if(nprint.ne.-5)write(iwr,204)m2,npass2
204   format(/' no. of sort',i1,' passes=',i4)
_IF(parallel)
       master=0
_ENDIF
      do 2 ipass=1,npass2
         i10=igmem_alloc_inf(mword,fnm,snm,'i10',IGMEM_DEBUG)
         i40=igmem_alloc_inf(nword2,fnm,snm,'i40',IGMEM_NORMAL)
         call adaps2(q(i10),q(i40))
         call gmem_free_inf(i40,fnm,snm,'i40')
         call gmem_free_inf(i10,fnm,snm,'i10')
         top=cpulft(1)
         if(nprint.ne.-5)write(iwr,2222)m2,ipass,top
_IF(parallel)
c MPP
         call pg_synch(3333)
c MPP
_ENDIF
         i10=igmem_alloc_inf(iv,fnm,snm,'i10',IGMEM_DEBUG)
         i20=igmem_alloc_inf(nananb,fnm,snm,'i20',IGMEM_NORMAL)
         call adapc2(cc,q(i10),q(i20),q)
         call gmem_free_inf(i20,fnm,snm,'i20')
         call gmem_free_inf(i10,fnm,snm,'i10')
         call dumpra(2,ipass,npass2,jblkaf,junitf)
         if(tim.gt.timlim)go to 2223
 2    continue
7755  top=cpulft(1)
      call clredx
      if(nprint.ne.-5) write(iwr,300)top,charwall()
300   format(/' end of 4-index transformation at',
     *f8.2,' seconds',a10,' wall'/)
      if(otran4) then
       if(nprint.ne.-5) then
        write(iwr,302)
 302    format(' status of transformed integral file'/
     *         ' ***********************************')
        call filprn(m6file,m6blk,m6last,m6tape)
       endif
_IFN(parallel)
       if(lfile.eq.1 .and. (m4tape(1).ne.m6tape(1)) ) then
        if(iscftp.ne.1) then
         if(m4tape(1).ne.m2tape(1)) call delfil(m4tape(1))
        else
         if(m4tape(1).ne.m1tape(1)) call delfil(m4tape(1))
        endif
       endif
_ENDIF
      else
       if(nprint.ne.-5) then
        write(iwr,301)
 301    format(
     +  ' status of symmetry adapted mainfile'/
     +  ' ***********************************')
        call filprn(m1file,m1blk,m1last,m1tape)
       endif
_IFN(parallel)
       if(lfile.eq.1 .and. m4tape(1).ne.m1tape(1) .and.
     +    m4tape(1).ne.m2tape(1) ) call delfil(m4tape(1))
_ENDIF
      endif
c
c     ----- reset core memory
c
 2223 call revise
c
      return
      end
      subroutine dumpra(iphase,ipass,npass,jbl,jun)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
INCLUDE(common/sizes)
INCLUDE(common/iofile)
INCLUDE(common/restar)
INCLUDE(common/files)
INCLUDE(common/trntim)
INCLUDE(common/timez)
INCLUDE(common/cndx40)
      common/blkin/gout(511)
      character*10 charwall
c
      data m0/0/
c
      if(ipass.ne.npass)goto 88
      call search(jbl,jun)
      call put(gout,m0,jun)
      jbl=jbl+1
 88   if(iphase-1)89,89,90
 89   m4file=nfiles
      m4tape(m4file)=n4tape(nfiles)
      m4blk (m4file)=n4blk (nfiles)
      m4last(m4file)=jbl
      go to 66
 90   if(.not.oindx4)go to 110
      m6file=nfilef
      m6tape(m6file)=n6tape(nfilef)
      m6blk (m6file)=n6blk (nfilef)
      m6last(m6file)=jbl
      go to 66
 110  m1file=nfilef
      m1tape(m1file)=n1tape(nfilef)
      m1blk (m1file)=n1blk (nfilef)
      m1last(m1file)=jbl
66    call search(isecbl,ndump4)
      m12=12/lenwrd()
      call put(master,m12,ndump4)
      b=cpulft(1)
      call clredx
      if(nprint.ne.-5)write(iwr,100)iphase,ipass,b,charwall()
100   format(/' job dumped after calc',i1,' pass',i6,' at ',f9.2,
     *' seconds',a10,' wall')
      if(ipass.eq.npass.and.iphase.eq.2)goto 99
      if(((b-slder)*1.8d0+b).lt.alimin)goto 77
      if(.not.oindx4)write(iwr,101)
 101  format( /10x,38('*')/
     * 10x,'symmetry adaption incomplete - restart'/ 10x,38('*')/)
      if(oindx4)write(iwr,102)
 102  format(/ 10x,44('*')/
     *10x,'integral transformation incomplete - restart'/10x,
     *44('*'))
      call secsum
      call whtps
      tim=timlim+0.1d0
 77   slder=b
99    return
      end
      subroutine adaps2(g,nijkl)
c
c...   sorts secondary mainfile onto lfn sort---so that
c...   for a given ij comb. all rs combs. available
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c
_IF1(iv)      integer *2 nijkl
      dimension g(*)
      dimension nijkl(*)
c
INCLUDE(common/blksiz)
INCLUDE(common/atmblk)
      common/stak/btri,mlow,nstack,iblock,mstack
INCLUDE(common/blocko)
      common/junke/maxt,ires,ipass,nteff,npass1,npass2,lentri,
     *nbuck,mloww,mhi,ntri,iacc
INCLUDE(common/cndx40)
      common/junk/nwbuck(1500),itx(3400),ktx(3400),gtx(3400)
      common/bufb/nwbnwb,lnklnk,gout(5119)
       common/blkin/gin(510)
      common/three/mark(1500),ibase(1500),ibasen(1500)
INCLUDE(common/filel)
_IF(parallel)
      common/nodinf/mpid,minode,mihost,ldim,nnodes,nodscr,
     * maxnod(maxlfn)
INCLUDE(common/prints)
_ENDIF
c
_IF(linux)
      external fget
_ENDIF
      data maxb/9999999/
_IF(parallel)
      call closbf(0)
      call setbfa(-1)
_ENDIF
c
c...   determine base and limit triangles for this pass
c
       mlow=iky(indxi)+indxj
      mloww=mlow-1
       mhi=mloww+nteff
      if(mhi.gt.lentri)mhi=lentri
      mtri=mhi-mlow
c...   determine minimum no. of bucks.
      nbuck=ires
_IF(parallel)
      if(oprint(58)) then
      write(6,33334) minode,nbuck, mloww,mhi
33334 format(1x,'minode,nbuck,mloww,mhi = ',4i10)
      endif
_ENDIF
400   ntri=mtri/nbuck
      if(ntri.ge.maxt)goto 401
_IF(parallel)
      nbuck=nbuck-nnodes
      if (nbuck.lt.0) nbuck=0
_ELSE
      nbuck=nbuck-1
_ENDIF
      if(nbuck)400,401,400
_IF(parallel)
 401   nbuck=nbuck+nnodes
      if(nbuck.gt.ires)then
      write(6,33343)minode,nbuck,ires
33343 format(' adaps2:minode,nbuck,ires = ',3i5)
       call caserr(
     *'invalid number of buckets in adaps2')
      endif
_ELSE
401   nbuck=nbuck+1
_ENDIF
      ntri=mtri/nbuck+1
_IF(parallel)
      if(oprint(58)) then
      write(6,33333) minode,nbuck, ires,ntri, mtri
33333 format(1x,'adaps2:minode,nbuck,ires,ntri,mtri = ',5i10)
      endif
_ENDIF
c...   ntri=max. no. of triangles controlled by 1 bucket
c...   nbuck=number of buckets
      btri=ntri
      btri=1.000000001d0/btri
      call setsto(nbuck,maxb,mark)
      call setsto(nbuck,0,nwbuck)
      nstack=0
      nwbnwb=nsz340
      iblock=0
c...   start loop over secondary mainfile blocks
      do 20113 ifile=1,lfile
      iunit=lotape(ifile)
      call search(liblk(ifile),iunit)
      lbl=llblk(ifile)
20110 call fget(gin,m,iunit)
      if(m)20121,20113,20121
20121 nnn=isort2(itx(nstack+1),ktx(nstack+1),gtx(nstack+1))
      nstack=nstack+nnn
      if(nstack.ge.3061) then
       call stackr(g,nijkl)
      endif
      lbl=lbl+1
       if(lbl)20110,20113,20110
20113  continue
c
c...   secondary mainfile now swept
c...   clear up output
c
      if(nstack.ne.0) then
        call stackr(g,nijkl)
      endif
c
      do 7 ibuck=1,nbuck
      nwb=nwbuck(ibuck)
      if(nwb)9999,7,9999
9999  ib=ibase(ibuck)
      ibn=ibasen(ibuck)
      call stopbk
      nwbnwb=nwb
      lnklnk=mark(ibuck)
      call dcopy(nwb,g(ib+1),1,gout,1)
_IF(ibm,vax)
      call fmove(nijkl(ib+1),gout(nsz341),nsz170)
_ELSE
      call pack(gout(nsz341),lab1632,nijkl(ibn+1),nsz680)
_ENDIF
      call sttout
      mark(ibuck)=iblock
      iblock=iblock+nsz
7     continue
      return
      end
      subroutine adapc2(q,r,v,qq)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
c...   processes sorted secondary mainfile on lfn sort
c...   to produce integrals of the form (i j/k l)
c...   on final mainfile
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
INCLUDE(common/restar)
c
      dimension v(*),r(*),q(*)
      dimension qq(*)
      common/junke/maxt,ires,ipass,nteff,npass1,npass2,lentri,
     *nbuck,mloww,mhi,ntri,iacc
INCLUDE(common/cndx40)
      common/bufb/nkk,mkk,gin(5119)
INCLUDE(common/blocko)
      common/blkin/gout(510),mword
      common/three/mark(1500),ibase(1500),ibasen(1500)
_IF1(iv)       common/craypk/mwcray(340),kbcray(340)
_IFN1(iv)       common/craypk/mwcray(680)
INCLUDE(common/iofile)
INCLUDE(common/mlngth)
INCLUDE(common/atmblk)
_IF(parallel)
INCLUDE(common/prints)
      common/nodinf/mpid,minode,mihost,ldim,nnodes,nodscr,
     * maxnod(maxlfn)
_ENDIF
      character*7 fnm
      character*6 snm
      data fnm,snm/"tran4.m","adapc2"/
      data maxb/9999999/
c
      acc1=10.0d0**(-iacc)
_IF(parallel)
c MPP
c MPP  all sorted // exchange  #'s of backchains
c MPP
      nbpn = (nbuck-1)/nnodes + 1
      i00 = igmem_alloc_inf((nbpn*nnodes+1)/lenwrd(),fnm,snm,
     &                      'bucket-addr',IGMEM_DEBUG)
_IF(ipsc)
      call excbuc(mark,nbuck,qq(i00),r,nnodes)
_ELSE
      call excbuc(mark,nbuck,qq(i00),r,r(nbpn+1),nnodes)
_ENDIF
      call gmem_free_inf(i00,fnm,snm,'bucket-addr')
_ENDIF
      mword=0
_IFN1(iv)      int4=-1
       call stopbk
_IF(parallel)
      call closbf(0)
      mloww = mloww + minode*ntri
      master = master + minode*ntri
      do 1 ibpn=1,nbpn
_ELSE
      do 1 ibuck=1,nbuck
_ENDIF
      mhigh=mloww+ntri
      if(mhigh.gt.mhi)mhigh=mhi
_IF(parallel)
      if(mloww.ge.mhi) go to 1
_ENDIF
      mtri= mhigh-mloww
      mloww=mloww+1
      nsize=mtri*lenb4
      i10 = igmem_alloc_inf(nsize,fnm,snm,'i10',IGMEM_NORMAL)
      call vclr(qq(i10),1,nsize)
c...   read in a core load
_IF(parallel)
      do 2 inode=1,nnodes
         ibuck = (ibpn-1)*nnodes+inode
         mkk=mark(ibuck)
         if(oprint(58)) then
         write(6,33335) minode,inode,ibuck,mkk
33335    format(' adapc2, minode,inode,ibuck,mkk = ',4i6)
         endif
         call setbfa(inode-1)
888      if(mkk.eq.maxb)go to 999
         iblock=mkk
         call rdbak(iblock)
         call stopbk
c...   move block to processing area
         call transc(qq(i10))
         goto 888
999      call closbf(0)
2     continue
      map=1
       if(oprint(58)) then
       write(6,33336) minode
33336 format(
     * ' adapc2:backchain complete minode',1x,i3)
       write(6,33333) minode,map,master,mtri,nsa4,nbas4,lena4
33333 format(
     * ' adapc2:minode,map,master,mtri,nsa4,nbas4,lena4 ',7i7)
       endif
_ELSE
      mkk=mark(ibuck)
888   if(mkk.eq.maxb)go to 999
      iblock=mkk
      call rdbak(iblock)
      call stopbk
c...  move block to processing area
      call transc(qq(i10))
      goto 888
999   map=1
_ENDIF
      do 777 itri=1,mtri
c...  compute integrals of the form (i j/k l)
_IF(parallel)
      master = master + 1
      call ijij(master,indxi,indxj)
_ENDIF
      mword1=indxj+i4096(indxi)
      call mult3(q,qq(i10+map-1),r,v,indxi,nbas4)
      m=0
      do 1009 k=1,indxi
      last=k
      if(k.eq.indxi)last=indxj
      do 1009 j=1,last
      m=m+1
      top=r(m)
c...    (i j/k l) now in top
      if( dabs(top).lt.acc1)go to 1009
      mword=mword+1
_IF1(iv)      kbcray(mword)=j+i4096(k)
_IF1(iv)      mwcray(mword)=mword1
_IFN1(iv)      int4=int4+2
_IFN1(iv)      mwcray(int4  )=mword1
_IFN1(iv)      mwcray(int4+1)=j+i4096(k)
      gout(mword)=top
_IF1(iv)      if(mword.eq.nintmx)call blockb
_IFN1(iv)      if(mword.lt.nintmx)go to 1009
_IFN1(iv)      call blockb
_IFN1(iv)      int4=-1
1009  continue
_IFN(parallel)
      indxj=indxj+1
      if(indxj.gt.indxi) then
        indxi=indxi+1
        indxj=1
      end if
777   map=map+lenb4
      call gmem_free_inf(i10,fnm,snm,'i10')
1     mloww=mhigh
_ENDIF
_IF(parallel)
777   map=map+lenb4
c
c ****** must reset indxi,indxj for multi-pass adaps2
       mmmm = master + 1
       call ijij(mmmm,indxi,indxj)
c
       mloww = mhigh + (nnodes-1)*ntri
       master = master + (nnodes-1)*ntri
       call gmem_free_inf(i10,fnm,snm,'i10')
1      continue
c MPP
c MPP   set mloww correctly for  multi-pass
c MPP
      mloww = mloww - minode*ntri
      master = master - minode*ntri
c MPP
_ENDIF
      if(mword.ge.1)call blockb
_IF(parallel)
c MPP
c MPP  synchronize
c MPP
      call pg_synch(3333)
_ENDIF
c
      return
      end
      subroutine symat(q)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
INCLUDE(common/prints)
INCLUDE(common/blocko)
INCLUDE(common/cndx40)
INCLUDE(common/tran)
INCLUDE(common/restar)
INCLUDE(common/iofile)
      character*7 fnm
      character*5 snm
      data fnm,snm/"tran4.m","symat"/
c
      dimension q(*)
c
      if(otran)call caserr(
     *'the basis set is not symmetry adapted')
c
      call comharm(flop,'ctrans',flop)
c
      if(nprint.ne.-5)write(iwr,98)
98    format( /14x,'list of lcbf'/  14x,12('*')/
     *' coefficient old orbital nterm new orbital'/1x,41('-')/)
      do 97 i=1,newb4
      n=ntran(i)
      j=ilifc(i)
      if(nprint.ne.-5)write(iwr,5038)ctran(j+1),itran(j+1),n,i
5038  format(f12.7,i12,i6,i12)
      if(n.eq.1)goto 97
      do 96 k=2,n
  96  if(nprint.ne.-5)write(iwr,5038)ctran(j+k),itran(j+k)
97    continue
      if(nprint.ne.-5)write(iwr,5037)newb4
5037  format(/' no. of lcbf=',i4)
c
c ----- get core memory
c
      nbsq=nbas4*nbas4
      i10 = igmem_alloc_inf(nbsq,fnm,snm,'i10',IGMEM_NORMAL)
c
      call vclr(q(i10),1,nbsq)
      do 100 i=1,ncol4
      iii=ilifq(i)+i10-1
 100  q(iii+i)=1.0d0
c
      call tdown(q(i10),ilifq,q(i10),ilifq,ncol4)
c
      if(oprint(46)) then
      write(iwr,110)
 110  format(/1x,104('-')/40x,16('*')/
     *                    40x,'s.a.b.f. vectors'/
     *                    40x,16('*'))
      call prsql(q(i10),ncol4,nbas4,nbas4)
      endif
c     call wrt3(q(i10),nbsq,iblkq4,num8)
c...   next write should be consistent with the one in vecget
      call wrt3(q(i10),nblkq4,iblkq4,num8)
c
c ----- reset memory
c
      call gmem_free_inf(i10,fnm,snm,'i10')
c
      return
      end
      subroutine adapt2(q,qvec,otran4)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
INCLUDE(common/gmempara)
INCLUDE(common/restar)
INCLUDE(common/cndx40)
INCLUDE(common/mlngth)
INCLUDE(common/filel)
INCLUDE(common/iofile)
c
      character*10 charwall
      character*7 fnm
      character*6 snm
      data fnm,snm/"tran4.m","adapt2"/
c
      dimension q(*)
c
      top=cpulft(1)
      if(nprint.ne.-5)write(iwr,2)top,charwall()
 2    format(/
     *' start of 2-index transformation at ',f8.2,' seconds',
     * a10,' wall')
c
      if(.not.otran4)then
         nreq = lenb4+lenb4+nbnb
         i10 = igmem_alloc_inf(nreq,fnm,snm,'i10s',IGMEM_NORMAL)
         i20=i10+lenb4
         i30=i20+lenb4
         call corbls(qvec,q(i10),q(i30),nprint,q)
         call gmem_free_inf(i10,fnm,snm,'i10s')
      else
         if(iscftp.ne.1) then
            lfile=m2file
            do 300 i=1,lfile
            lotape(i)=m2tape(i)
            liblk(i) =m2blk(i)
 300        llblk(i) =liblk(i)-m2last(i)
         else
            lfile=m1file
            do 600 i=1,lfile
            lotape(i)=m1tape(i)
            liblk(i) =m1blk(i)
 600        llblk(i) =liblk(i)-m1last(i)
         endif
         nreq = lenb4 + lenb4 + max(lenb4+lenb4, nananb)
         i10 = igmem_alloc_inf(nreq,fnm,snm,'i10d',IGMEM_NORMAL)
         i20 = i10+lenb4
         i30 = i20+lenb4
         call corbld(qvec,q(i10),q(i30),q(i20),nprint,q)
         call gmem_free_inf(i10,fnm,snm,'i10d')
      endif
      top=cpulft(1)
      if(nprint.ne.-5)write(iwr,4)top,charwall()
 4    format(/ ' end of 2-index transformation at ',f8.2,' seconds',
     * a10,' wall')
      return
      end
      subroutine corbls(q,q1,q2,nprint,qq)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c
      dimension q(*),q1(*),q2(*),qq(*)
      dimension yprin(6),array(6),o1e(6)
INCLUDE(common/cndx40)
      common/blkin/bufa(512),potnuc,pp(3),qone(508)
INCLUDE(common/blocko)
INCLUDE(common/scra7)
INCLUDE(common/iofile)
INCLUDE(common/infoa)
INCLUDE(common/dump3)
INCLUDE(common/mlngth)
INCLUDE(common/zorac)
cjvl     re zora cf. corbld
c
      data yprin/'s','t','t+v','x','y','z'/
      data m2,m511/2,511/
         call check_feature('corbls')
c
      call rdedx(q,nanb,iblkq4,num8)
c
      lenb4=iky(nsa4+1)
      lenblk=lensec(lenb4)
      do loop =1,6
       o1e(loop) = .false.
      enddo
c
c     restore header block
c
      call getmat(q2,q2,q2,q2,q2,q2,potnuc,num,o1e,ionsec)
c
      nbloc = lenblk*6 + 1
      call secput(ions2,m2,nbloc,ib)
      if(nprint.ne.-5)write(iwr,103)ions2
      kblkkk=ib
      call wrt3(potnuc, m511, kblkkk, ndump4 )
      kblkkk = kblkkk + 1
c     nbl=lensec(lbas)
c
      do imat=1,6
      o1e(imat) = .true.
      call getmat(q2,q2,q2,q2,q2,q2,array,num,o1e,ionsec)
c
      if (imat.eq.3.and.ozora) call zora(qq,q1,q2,'read')
c
      call mult2(q,q1,q2,nsa4,nsa4,nbas4)
      call wrt3(q1,lenb4,kblkkk,ndump4)
      kblkkk=kblkkk+lenblk
      if(.not.oprin4(imat))then
       write(iwr,3337)yprin(imat)
       call writel(q1,nsa4)
      endif
      o1e(imat) = .false.
      enddo
c
      return
 103  format(/
     *' *** output symmetry adapted 1e-integrals to section',
     *i4,' of the dumpfile')
 3337  format(/1x,104('*')//
     * 35x,a3,' matrix over symmetry adapted basis functions'/
     * 35x,48('-')/)
      end
      subroutine reftym(q,nprint)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      common/blkin/pppp(512),
_IF1(ck)     *              potnuc,dx,dy,dz,nato,numo,sm(502),kword,
_IF1(ck)     *              pad(3)
_IFN1(ck)     *              potnuc,dx,dy,dz,nato,numo,sm(503),kword,
_IFN1(ck)     *              ipad(7)
INCLUDE(common/machin)
INCLUDE(common/cndx40)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/scra7)
INCLUDE(common/infoa)
c
      dimension q(*), j1e(6)
      data m2/2/
c
c     indexing triangle storage on num8
c
      j1e(5)=0
      j1e(6)=1
      do loop=1,4
       j1e(loop)=loop+1
      enddo
c
      if(nprint.ne.-5)write(iwr,103)ionsv4,potnuc,dx,dy,dz
c
      lbas=nbas4*(nbas4+1)/2
      lenl = lensec(lbas)
      nbloc = lenl*6 + 1
      call secput(ionsv4,m2,nbloc,ib)
      lbas = (nsa4*(nsa4+1))/2
      nato = nat
      numo = nbas4
c     header block 
      call wrt3(potnuc,511,ib,ndump4 )
c
      nbl=lensec(lbas)
c     input order from num8:    T+V, X, Y, Z, S, T
c     output order to dumpfile: S, T, T+V, X, Y, Z
      ii=ibl7la + lensec(nblkq4)
      do i=1,6
       ioff=j1e(i)*lenl + 1 + ib
       call rdedx(q,lbas,ii,num8)
       call wrt3(q,lbas,ioff,ndump4)
       ii=ii+nbl
      enddo
c
      call revind
c
      return
 103  format(//
     *' *** output transformed 1e-integrals (s,t,t+v,x,y,z) to ',
     *'section ',i4,' of the dumpfile'//
     * 4x,'effective nuclear terms'/
     * 4x,'======================='//
     * 4x,'core',1x,f20.12/
     * 4x,' x  ',1x,f20.12/
     * 4x,' y  ',1x,f20.12/
     * 4x,' z  ',1x,f20.12//)
      end
      subroutine corbld(qvec,q1,q2,q3,nprint,q)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      dimension qvec(*),q1(*),q2(*),q3(*),q(*)
INCLUDE(common/scra7)
INCLUDE(common/prnprn)
INCLUDE(common/cndx40)
      common/blkin/corev(512),potdum(4),p(508)
cjvl   *note* due to call's to a fock-matrix builder in zora
cjvl   *note* blkin may be overwritten (intega)
cjvl   ****** now only atomic zora is supported so problem is gone
INCLUDE(common/blocko)
c    + , ipad(3*maxorb+5)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/machin)
      common/junkc/zcomm(19),ztext(10)
      common/lsort/iijjnk(4*maxorb),
     *pop(maxorb),potn,core,ncolo,nbas,newb,ncore,
     *mapcie(maxorb),ilifc(maxorb),nval,mapaie(maxorb),ilifa(maxorb),
     *iqsec,jdsec
INCLUDE(common/filel)
      common/craypk/ijkl(1360)
INCLUDE(common/field)
INCLUDE(common/harmon)
INCLUDE(common/atmblk)
INCLUDE(common/zorac)
       dimension j1e(6), field(3),o1e(6),potnuc(6)
       equivalence (field(1),fieldx)
_IF(linux)
      external fget
_ENDIF
      call check_feature('corbld')
c
_IF1(c)      call fmove(mapie,mapaie,nsa4)
_IF1(c)      call fmove(ilifm,ilifa,nsa4)
_IFN1(c)      call icopy(nsa4,mapie,1,mapaie,1)
_IFN1(c)      call icopy(nsa4,ilifm,1,ilifa,1)
      ncolo=ncol4
c     newb =newb4 (only used to write out to section jdsec) 
      newb =nbas4
      nbas =nbas4
      lenblk=lensec(lena4) + lensec(mach(15))
      if(iscftp.eq.1) ionset = ions2
      if(iscftp.eq.0) ionset = ionsec
c
      lenv=iky(nsa4+1)
      lenvvv=lensec(lenv)
      kblkkk=ibl7la + lensec(nblkq4)
c
c     transform 1-electron integrals in turn
c     these are processed in the order T+V, X, Y, Z, S and T
c     j1e indexes this ordering
c
      do 5000 loop=1,4
5000  j1e(loop)=loop+2
      j1e(5)=1
      j1e(6)=2
c
      imatr = 0
4000  imatr = imatr + 1
      if(imatr.gt.6) go to 4010
c
_IFN1(c)      call setstl(6,.false.,o1e)
_IF1(c)      call setsto(6,.false.,o1e)
      o1e(j1e(imatr)) = .true.
c
      call getmat (q3,q3,q3,q3,q3,q3,potnuc,nbas,o1e,ionset)
c
      o1e(j1e(imatr)) = .false.
      if(imatr.eq.1 .and. ofield) then
c
c     must restore x, y, and z, adding into x3 (1-electron H)
c
       do i = 1, 3
        o1e(i+3) = .true.
        call getmat (q2,q2,q2,q2,q2,q2,potnuc,nbas,o1e,ionset)
        ffield = field(i)
         do j = 1, lenb4
         q3(j) = q3(j) + ffield * q2(j)
         enddo
        o1e(i+3) = .false.
       enddo
      endif
c
      if(imatr.eq.1.and.ozora.and.iscftp.ne.1) call zora(q,q1,q3,'read')
c
      if(imatr.gt.4) then
         potn = 0.0d0
      else
         potn=potnuc(imatr)
      endif
      core=potn
c
      if (ncore.gt.0) then
c
        call rdedx(q2,nblkq4,iblkq4,num8)
        if(iscftp.eq.0)call tdown(q2,ilifq,q2,ilifc,ncore)
        call vclr(q1,1,lenb4)
        do 2 k=1,ncore
        mm=ilifq(k)
        m=0
        do 2 i=1,nbas4
        top=q2(mm+i)
        p(i)=top
        call daxpy(i,top,p,1,q1(m+1),1)
 2      m=m+i
        do 9 i=1,nbas4
        m=iky(i+1)
    9   q1(m)=q1(m)*0.5d0
        top=ddot(lenb4,q1,1,q3,1)
        if(imatr.eq.1) then
_IF1(civ)          call szero(q2,lenb4)
_IFN1(civ)          call vclr(q2,1,lenb4)
          call setsto(1360,0,ijkl)
c
_IF1(iv)          do 5 i=1,nbas4
_IF1(iv)          m=iky(i+1)
_IF1(iv)    5     q1(m)=q1(m)*2.0d0
          do 20113 i=1,lfile
          lbl=llblk(i)
          if(lbl)10000,20113,10000
10000     iunit=lotape(i)
          call search(liblk(i),iunit)
 68       call fget(corev,k,iunit)
          if(k)1,20113,1
_IF(ibm,vax)
    1     call gmake(q2,q1,iky)
_ELSE
    1     if(o255i)then
           call sgmata(q2,q1)
          else
           call sgmata_255(q2,q1)
          endif
_ENDIF
          lbl=lbl+1
          if(lbl)68,20113,68
20113     continue
_IF1(civ)          call addvec(q3,q3,q2,lenb4)
_IF(parallel)
c MPP
c MPP  gather and spread the results
c MPP
             call pg_dgop(6011,q2,lenb4,'+')
c MPP
_ENDIF
_IFN1(civ)        call vadd(q3,1,q2,1,q3,1,lenb4)
_IF1(iv)        do 6 i=1,nbas4
_IF1(iv)        m=iky(i+1)
_IF1(iv)    6   q1(m)=q1(m)*0.5d0
        endif
        top=ddot(lenb4,q1,1,q3,1)+top
        core=top+top+core
      endif
      corev(imatr)=core
      call mult3(qvec,q3,q1,q2,nsa4,nbas4)
      if(imatr.eq.1) then
      if(ofield)
     * write(iwr,230)fieldx,fieldy,fieldz
 230  format(/5x,38('*')/
     *        5x,'one-electron hamiltonian modified with'/
     *        5x,'x *',f11.7/
     *        5x,'y *',f11.7/
     *        5x,'z *',f11.7/
     *        5x,38('*') )
      if(nprint.ne.-5) write(iwr,601)core
 601  format(/
     *' effective nuclear repulsion term',e21.12)
c
      if (iscftp.eq.1) then
c...     reset ilifq,etc if changed by trani for mcscf=> ci
         do i=1,newbas1
            ilifq(i) = (i-1)*newbas1
         end do
         ncolo = newbas1
         nbas  = newbas1
         newb = newbas1
      end if
c
      call secput(jdsec,1004,lenblk,kblxxx)
      nval=nsa4
      call wrt3(pop,mach(15),kblxxx,ndump4)
      call wrt3s(q1,lena4,ndump4)
      endif
      if(.not.oprin4(7)) then
       write(iwr,5003)
 5003  format(/1x,104('*')//
     * 35x,'core matrix over active mos'/35x,27('*')/)
       call writel(q1,nsa4)
      endif
      if (oprin4(10).and.opunch(12)) then
c
c     junk for writing out the one electron integrals for rjh
c
        do 1011 ii = 1,nsa4
          do 1011 jj = 1,ii
            iijj = ii*(ii-1)/2 + jj
            write(7,*) q1(iijj),ii,jj,0,0
1011    continue
        write(ipu,*) core,0,0,0,0
        write(ipu,*) 0.0d0,-1,-1,-1,-1
        close(ipu,status='keep')
      endif
c
      if(ionsv4.eq.0) return
      call wrt3(q1,lena4,kblkkk,num8)
      kblkkk=kblkkk+lenvvv
      go to 4000
c
4010  call reftym(qvec,nprint)
      return
      end
      subroutine pr2ef(iunit,iblock)
      implicit REAL  (a-h,p-w), integer (i-n), logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/iofile)
INCLUDE(common/discc)
INCLUDE(common/atmblk)
      common/blkin/g(510),mword,mspace,
_IFN1(iv)     * i205(4,340)
_IF1(iv)     * i205(340),j205(340),k205(340),l205(340)
      dimension f(4),ia(4),ja(4),ka(4),la(4)
c
      write(iwr,100)iblock,yed(iunit)
      n = 0
      if(mword)33,44,33
 44   write(iwr,112)
      go to 404
 33   continue
_IFN1(c)      call setsto(numlab,0,i205)
_IF1(iv)      call upak8v(g(num2e+1),i205)
_IFN1(iv)      call unpack(g(num2e+1),lab816,i205,numlab)
      do 3 m=1,mword
      n=n+1
      f(n)=g(m)
_IF(ibm,vax)
      ia(n)=i205(m)
      ja(n)=j205(m)
      ka(n)=k205(m)
      la(n)=l205(m)
_ELSEIF(littleendian)
      ia(n)=i205(2,m)
      ja(n)=i205(1,m)
      ka(n)=i205(4,m)
      la(n)=i205(3,m)
_ELSE
      ia(n)=i205(1,m)
      ja(n)=i205(2,m)
      ka(n)=i205(3,m)
      la(n)=i205(4,m)
_ENDIF
      if(n.lt.4 )goto 3
      write(iwr,101)(ia(n),ja(n),ka(n),la(n),f(n),n=1,4)
      n=0
3     continue
      if(n.ne.0 )write(iwr,101)(ia(m),ja(m),ka(m),la(m),f(m),m=1,n)
      write(iwr,102)mword
404    return
112   format(/3x,'endfile block')
102   format(/' no. of integrals =',i4)
101   format(4(1x,4i4,f15.8))
 100  format(/1x,104('*')//45x,'list of block',i6,' from ',a4/
     *45x,27('-')//
     *4(4x,'i',3x,'j',3x,'k',3x,'l',6x,'g',8x)/
     *1x,115('+')/)
      end
      subroutine pr2es(iunit,iblock)
      implicit REAL  (a-h,p-w), integer (i-n), logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/iofile)
INCLUDE(common/discc)
INCLUDE(common/atmblk)
      common/blkin/g(510),mword,mspace,
_IF1(iv)     * ij205(340),kl205(340)
_IFN1(iv)     * ij205(2,340)
      dimension f(4),ia(4),ja(4)
c
      write(iwr,100)iblock,yed(iunit)
      n = 0
      if(mword)33,44,33
 44   write(iwr,112)
      go to 404
 33   continue
_IF1(iv)      call upak4v(g(num2ep+1),ij205)
_IFN1(iv)      call unpack(g(num2ep+1),lab1632,ij205,numlabp)
      do 3 m=1,mword
      n=n+1
      f(n)=g(m)
_IF(ibm,vax)
      ia(n)=ij205(m)
      ja(n)=kl205(m)
_ELSE
      ia(n)=ij205(1,m)
      ja(n)=ij205(2,m)
_ENDIF
      if(n.lt.4 )goto 3
      write(iwr,101)(ia(n),ja(n),f(n),n=1,4)
      n=0
3     continue
      if(n.ne.0 )write(iwr,101)(ia(m),ja(m),f(m),m=1,n)
      write(iwr,102)mword
404    return
112   format(/3x,'endfile block')
102   format(/' no. of secondary mainfile elements =',i4)
101   format(4(1x,2i6,f16.9))
 100  format(/1x,104('*')//45x,'list of block',i6,' from ',a4/
     *45x,27('-')//
     *4(5x,'ij',4x,'kl',6x,'g',9x)/
     *1x,104('+')/)
      end
_IF(parallel)
_IF(ipsc)
      subroutine excbuc(mark,nbuck,immm,ittt,nnod)
c MPP
c MPP  exhange bucket-addresses between nodes
c MPP  nnod = nnodes
c MPP  in immm (which will be moved back to mark) at the end
c MPP  each node will contain the addresses for bucket
c MPP  minode+1  + nnodes*(i-1) (i=1,2,...)
c MPP
      implicit integer (a-z)
INCLUDE(common/sizes)
c
      common/nodinf/mpid,minode,mihost,ldim,nnodes,nodscr
     * ,maxnod(maxlfn)
c
      dimension mark(nbuck),immm(nnod,*),ittt(*)
      parameter (itype=14328)
c
      nbpn = (nbuck-1)/nnodes + 1
c
c..   loop over all other nodes
c
      do 100 i=1,nnodes
         ito = mod(minode+i,nnodes)
c MPP
c MPP     gather start of buckets for this node
c MPP
         call setsto(nbpn,9999999,ittt)
         kk = 0
         do 10 j=ito+1,nbuck,nnodes
            kk = kk + 1
10       ittt(kk) = mark(j)
c MPP
c MPP     now send this to node ito
c MPP
         if (ito.ne.minode) then
          call csend(itype,ittt,nbpn*4,ito,mpid)
         endif
c MPP
c MPP     receive stuff set up by another node
c MPP
         if (ito.ne.minode) then
            call cprob(itype)
            call crecv(itype,ittt,nbpn*4)
            ifr = infonode()
         else
            ifr = ito
         end if
c
         do 20 j=1,nbpn
20       immm(ifr+1,j) = ittt(j)
c
100   continue
c
      call icopy(nbpn*nnodes,immm,1,mark,1)
c MPP
c MPP  this synchronizes automatically
c MPP  all have to wait to get all the backchains
c MPP
      return
      end
_ELSE
      subroutine excbuc(mark,nbuck,immm,ittt,isss,nnod)
c MPP
c MPP  exhange bucket-addresses between nodes
c MPP  nnod = nnodes
c MPP  in immm (which will be moved back to mark) at the end
c MPP  each node will contain the addresses for bucket
c MPP  minode+1  + nnodes*(i-1) (i=1,2,...)
c MPP
c     If communication is synchronous then we must explictly pair up
c     send/receive requests on this process with the matching
c     receive/send operations on other processes.
c
      implicit integer (a-z)
INCLUDE(common/sizes)
c
      common/nodinf/mpid,minode,mihost,ldim,nnodes,nodscr
     * ,maxnod(maxlfn)
INCLUDE(common/iofile)
c
      dimension mark(nbuck),immm(nnod,*),ittt(*),isss(*)
      parameter (itype=14328)
c
      nbpn = (nbuck-1)/nnodes + 1
      nn = nnodes + mod(nnodes,2)
      nbpn4 = nbpn * 8 / lenwrd()
c
        do 10 iter = 1, nn-1
          call pairrup(nn, minode, iter, ip)
          if (ip.lt.nnodes) then
           call setsto(nbpn,9999999,ittt)
           kk = 0
           do 30 j=ip+1,nbuck,nnodes
               kk = kk+1
               ittt(kk) = mark(j)
 30         continue
            if (minode. lt. ip) then
              call pg_snd(itype, ittt, nbpn4, ip, 1)
              call pg_rcv(itype, isss, nbpn4, lenmes, ip, node, 1)
            else if (minode.gt.ip) then
               call pg_rcv(itype, isss, nbpn4, lenmes, ip, node, 1)
               call pg_snd(itype, ittt, nbpn4, ip, 1)
            endif
            do 20 loop = 1,nbpn
               immm(ip+1,loop) = isss(loop)
 20         continue
          endif
 10     continue
c
c..   finally loop over all nodes
c
      do 2000 i=1,nnodes
         ito = mod(minode+i,nnodes)
         if(ito.eq.minode) then
          call setsto(nbpn,9999999,ittt)
          kk = 0
          do 2010 j=ito+1,nbuck,nnodes
               kk = kk + 1
               ittt(kk) = mark(j)
 2010       continue
c
          do 2020 j=1,nbpn
               immm(ito+1,j) = ittt(j)
 2020       continue
         endif
c
 2000 continue
c
      call icopy(nbpn*nnodes,immm,1,mark,1)
c
      return
      end
      subroutine pairrup(n, me, iter, ipair)
      implicit REAL  (a-h,o-z)
c
c     one of many ways of generating maximally overlapped pairs
c     (not all that good on a hypercube though!)
c
      if (iter.eq.1) then
        ipair = mod(n+1-me,n)
      else if (me.eq.0) then
        ipair = iter
      else if (me.eq.iter) then
        ipair = 0
      else
        if (ipair.eq.0) ipair = me
        ipair = ipair + 2
        if (ipair.ge.n) ipair = ipair + 1 - n
      endif
      return
      end
_ENDIF
      subroutine ijij(master,i,j)
c
c...   get i and j out of triangle index master
c
      implicit REAL  (a-h,o-z)
c
      i = (dsqrt(8.0d0*master+1.0d0) - 1.0d0)/2.0d0
      i = i + 1
      j = master - i*(i-1)/2
      if (j.eq.0) then
         i = i - 1
         j = i
      end if
c
      return
      end
_ENDIF
      subroutine ver_tran4(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/tran4.m,v $
     +     "/
      data revision /"$Revision: 6176 $"/
      data date /"$Date: 2010-08-10 16:49:47 +0200 (Tue, 10 Aug 2010) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
