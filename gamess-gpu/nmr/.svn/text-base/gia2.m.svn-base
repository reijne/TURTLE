c
c     NMR code from HONDO 2005
c     taken with kind permission of M. Dupuis
c
      subroutine gia2ei(inext)
      implicit REAL (a-h,o-z)
      logical some
      common/iofile/ir,iw
c
      some=.false.
      call fiflsh(iw)
c
      call cpuwal(tim0,wtim0)
c
c     ----- set up -giao- integral calculation -----
c
      iaddr=inext
      call giaset
c
c     ----- -giao- integrals -----
c
      call giamem(iaddr)
      call giapos
      call giaint
c
      call cpuwal(tim,wtim)
       dtim= tim- tim0
      dwtim=wtim-wtim0
      if(some) write(iw,9999)  dtim, tim
      if(some) write(iw,9998) dwtim,wtim
      return
 9999 format(10x,'elapsed - cpu- time = ',f10.3,
     1        5x,'total - cpu- time = ',f10.3)
 9998 format(10x,'elapsed -wall- time = ',f10.3,
     1        5x,'total -wall- time = ',f10.3)
      end
      subroutine giaint
      implicit REAL (a-h,o-z)
      parameter (mxatom=MAXAT)
      parameter (mxshel=MXSHEL)
      parameter (mxprim=MXPRIM)
      parameter (mxioda=1000) 
      logical   out 
      logical   spskip
INCLUDE(../m4/common/vcore)
INCLUDE(../m4/common/restri)
      common/iofile/ir,iw,ipu(3),idaf
      common/hdafile/idafh,navh,ioda(2,mxioda)
      common/inthnd/spskip
      common/nshel/
     1ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),cf(mxprim),cg(mxprim),
     1kstart(mxshel),katom(mxshel),ktype(mxshel),kng(mxshel),
     2kloc(mxshel),kmin(mxshel),kmax(mxshel),nshell
      common/intmem/iclb,inijg,ibuf,igint,iijklg,igijkl,ignkl,
     1 ignm,idij,idkl,ib00,ib01,ib10,ic00,id00,if00,
     2 isj,isk,isl,idijsi,idijsj,idklsk,idklsl,iabv,icv,irw,ichrg
      data m25/25/
c
      out =.false.         
c
      call cpuwal(tim0,wtim0)
c
c     ----- read -coulomb- integral bounds on -daf- -----
chvd        Schwarz integrals
c
      nclb=(nshell*(nshell+1))/2
      call secget(isect(421),m25,iblk25)
      call rdedx(Q(iclb),nclb,iblk25,idaf)
      do i=0,nclb-1
         Q(iclb+i) = dexp(2*Q(iclb+i))
      enddo
chvd
      if(out) then
         write(iw,*) '-coulomb- integral bounds'
         call prtr(Q(iclb),nshell)
         call fiflsh(iw)
      endif
c
c     ----- -giao- integrals by rys quadrature technique -----
c
      spskip=.false.       
      call giahnd(Q(iclb),Q(inijg))
c
      call cpuwal(tim,wtim)
      dhnd=tim-tim0
      if(out) then
         write(iw,9999) dhnd
      endif
c
      return
 9999 format(' -giahnd- elapsed - cpu- time = ',f10.3)
      end
      subroutine giaset
      implicit REAL (a-h,o-z)
      parameter (mxatom=MAXAT)
      parameter (mxshel=MXSHEL)
      parameter (mxprim=MXPRIM)
      parameter  (mxsym =48)
      parameter  (mxioda=1000)
      character*8 errmsg
      character*8 dsk1,dsk2
      character*8 dsk3,dsk4
      character*8 dsk5,dsk6
      character*8 xio1,xio2
      character*8 xio3,xio4
      character*8 xio5,xio6
      character*8 pck2,pck4
      character*8 uhf
      character*8 semid
      character*8 direct
      character*8 rfield
      character*8 blocki
      character*8 rname
      logical     pack2e
      logical     pk,pandk
      logical     block
      logical     some
      logical     out
      logical     dbug
      logical     sp
      logical     nspdsk,nodisk
      integer     thim,thime,semisp,dirstp
      common/output/nprint
      common/iodisk/mxdisk
INCLUDE(../m4/common/symtry)
      common/iofile/ir,iw,ip,ipu(2),idaf
      common/sqfile/ijk,ipk
      common/hdafile/idafh,navh,ioda(2,mxioda)
      common/restar1/nrest
      common/machin1/isingl,nbits
      common/dirpar/thiz,thim,semisp,dirstp
      common/pcklab/labsiz
      common/pckopt/nhex,ntupl,pack2e
      character*8 fldtyp
      common/drfopt/fldtyp
      character*8 hndtyp,dsktyp,filtyp,pcktyp
      common/hndopt/hndtyp,dsktyp,filtyp,pcktyp
      common/scfopt1/scftyp
      common/intopt/nopk,nok,nosqur
      character*8 blktyp
      common/intsym/blktyp
      common/inttyp/ipople,ihondo
      common/intprt/qint(3),valint(3),jcint(16),out
      common/intfil/nintmx
      common/intdsk/maxrec,ijklnu,ijklsp,ijklno,nspdsk,nodisk
      common/intcal/jkskip,notei
      common/baspar/normf,normp,itol
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,zan(mxatom),
     1                                         c(3,mxatom)
      common/nshel/
     1ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),cf(mxprim),cg(mxprim),
     1kstart(mxshel),katom(mxshel),ktype(mxshel),kng(mxshel),
     2kloc(mxshel),kmin(mxshel),kmax(mxshel),nshell
      common/shlord/modshl(mxshel),invshl(mxshel)
      common/shltol/rtol,dtol
      common/giao_shlint/cut,itot,icount,
     1              irec,locint,istrt,jstrt,kstrt,lstrt
      common/pkfil/pk,pandk,block
      common/symshl/mapshl(mxshel,mxsym),mapctr(mxatom,mxsym)
      common/isopac/indin(mxsym),indout(mxsym/2)
      dimension mi(mxsym)
      dimension rname(8)
      dimension iso(mxshel,mxsym/2)
      dimension errmsg(3)
      equivalence (iso(1,1),mapshl(1,1))
      data errmsg /'program ','stop in ','-jkset -'/
c
      data dsk1,dsk2 /'dsk-jk-u','dsk-jk-o'/
      data dsk3,dsk4 /'dsk-p -u','dsk-p -o'/
      data dsk5,dsk6 /'dsk-pk-u','dsk-pk-o'/
      data xio1,xio2 /'xio-jk-u','xio-jk-o'/
      data xio3,xio4 /'xio-p -u','xio-p -o'/
      data xio5,xio6 /'xio-pk-u','xio-pk-o'/
      data pck2,pck4 /'-2-label','-4-label'/
c
      data iblk   /'    '/
      data uhf    /'uhf     '/
      data semid  /'semi-dir'/
      data direct /'direct  '/
      data rfield /'rfield  '/
      data blocki /'blocki  '/
      data zero   /0.0d+00/
      data tenm6  /1.1d-06/
      data ten    /1.0d+01/
      data rln10  /2.30258d+00/
      data rname  /' hondo  ','integral',' file.  ','        ',
     1             '        ','        ','        ','        '/
      data ntyp   /6/
      data numjk  /255/
      data numpk  /362/
c
      data nhondo /0/
      data npkfil /1/
      data nkfil  /1/
      data labfil /0/
      data cutoff /1.0d-09/
      data intpck /0/
      data intskp /0/
      data intgrp /0/
      data intbuf /65536/
      data itotal /0/
      data nrec   /1/
      data intloc /1/
      data ist    /1/
      data jst    /1/
      data kst    /1/
      data lst    /1/
      data npflg  /0/
      data thime  /5/
      data thize  /1.0d-04/
      data mbytes /-1/
c
      namelist /intgrl/ intskp,intpck,intgrp,intbuf,nhondo,npkfil,nkfil,
     1                  labfil,npflg,cutoff,thize,thime,
     2                  nrec,intloc,ist,jst,kst,lst,itotal,mbytes
c
      dbug=.false.
      some=.false.      
      some=some.or.dbug
c
      if(some) then
         write(iw,9999)
      endif
c
c     ----- define default file format -----
c
      if((filtyp.eq.dsk1.or.filtyp.eq.dsk2.or.
     1    filtyp.eq.xio1.or.filtyp.eq.xio2    )) then
            npkfil=1
             nkfil=1
      endif
      if((filtyp.eq.dsk3.or.filtyp.eq.dsk4.or.
     1    filtyp.eq.xio3.or.filtyp.eq.xio4    )) then
            npkfil=0
             nkfil=1
      endif
      if((filtyp.eq.dsk5.or.filtyp.eq.dsk6.or.
     1    filtyp.eq.xio5.or.filtyp.eq.xio6    )) then
            npkfil=0
             nkfil=0
      endif
c
c     ----- read namelist $intgrl -----
c
chvd  suppress reading integral input for now
chvd
chvd  call rewfil(ir)
chvd  read(ir,intgrl,end=10,err=10)
   10 continue
c
      out=npflg.ne.0
c
c     ----- reaction field integral option -----
c
      if(fldtyp.eq.rfield) nhondo=1
      if(fldtyp.eq.rfield) npkfil=1
      if(fldtyp.eq.rfield)  nkfil=1
c
c     ----- for compatibility -----
c
      nosqur=1
c
c     ----- file format -----
c
      if(labfil.eq.2) then
         pcktyp=pck2
      elseif(labfil.eq.4) then
         pcktyp=pck4
      endif
      if(npkfil.ne.0) then
         nkfil=1
         filtyp=dsk1
         pcktyp=pck4
      else
         if(nkfil.ne.0) then
            filtyp=dsk3
         else
            filtyp=dsk5
         endif
      endif
c
c     ----- for -giao- -----
c
      npkfil=1
      nkfil =1
c
      nopk =npkfil
      nok  =nkfil
      pk   =nopk.eq.0
      pandk=pk.and.(nok.eq.0)
c
      nintmx=intbuf
c
      if(dbug            ) write(iw,9995) nintmx
      if(dbug            ) write(iw,9986) labsiz
      if(dbug.and.     pk) write(iw,9998) pandk
      if(dbug.and..not.pk) write(iw,9997)
c
c     ----- verify file format -----
c
      if((.not.pk).and.(num   .gt.numjk).and.
     1                 (isingl.eq.    2).and.(labsiz.ne.2)) then
         write(iw,9988) isingl,labsiz,num,numjk
         call hnderr(3,errmsg)
      endif
      if((     pk).and.(num   .gt.numpk).and.
     1                 (isingl.eq.    2).and.(labsiz.ne.2)) then
         write(iw,9987) isingl,labsiz,num,numpk
         call hnderr(3,errmsg)
      endif
      if((filtyp.eq.dsk1.or.filtyp.eq.dsk2.or.
     1    filtyp.eq.xio1.or.filtyp.eq.xio2    ).and.
     2   (pcktyp.eq.pck2                      )     ) then
            write(iw,9989) filtyp,pcktyp
            call hnderr(3,errmsg)
      endif
c
      ipople=1
      ihondo=0
      do 20 i=1,nshell
      if(ktype(i).le.2) go to 20
      ihondo=1
      go to 30
   20 continue
   30 continue
      if(nhondo.gt.0) ipople=0
      if(nhondo.gt.0) ihondo=1
      if(dbug.and.ipople.eq.1) write(iw,9992)
      if(dbug.and.ipople.eq.0) write(iw,9991)
      if(dbug.and.intskp.ne.0) write(iw,9993) intskp
c
      rtol=rln10*itol
      dtol=ten**(-itol)
      cut=cutoff
      if(dbug) write(iw,9994) cutoff
      dum=dtol/cutoff
      if(dum.gt.tenm6) then
         write(iw,9990) dtol,cutoff
         call hnderr(3,errmsg)
      endif
c
c     ----- direct and semi-direct algorithm options -----
c           for semi-direct integral selection is based
c           on 'thime and mbytes' or 'thime and thize'
c           but not on all three parameters active at
c           the same time.
c
      thim  =thime
      thiz  =thize
      thiz  =max(thiz,cutoff)
c
      dirstp=0
c
c     ----- figure how many integral buffers will go to disk -----
c
      if(mbytes.lt.0) then
         maxrec=2**30
      else
         thiz  =zero
         mbytes=abs(mbytes)
         lendsk=mbytes*1024*1024
         if(pandk) then
            lendsk=lendsk/2
         endif
         lenrec=nintmx*(8+labsiz*(8/isingl))+(8/isingl)
         maxrec=lendsk/lenrec
         maxrec=maxrec-(maxrec/100+1)
         if(maxrec.lt.0) then
            maxrec=0
         endif
      endif
      if(dbug.and.hndtyp.eq.semid ) then
         write(iw,9985) thim,thiz
         if(mbytes.ge.0) then
            write(iw,9984) mbytes,maxrec
            write(iw,9983)
         endif
      endif
c
      itot  =itotal
      irec  =nrec
      locint=intloc
      istrt =ist
      jstrt =jst
      kstrt =kst
      lstrt =lst
c
      jkskip=intskp
      if(notei.ne.0) jkskip=notei
c
c     ----- symmetry blocking parameter ( not used ) -----
c
      block=blktyp.eq.blocki
c
c     ----- initialize print out parameter -----
c
      jcint(1)=0
c
c     -----             packing parameters                  -----
c           nhex = desired hexadecimal accuracy.
c           ntupl= ' of bytes per integer word to be packed.
c           the negative value of -ntupl- means that the
c           integer labels are not to be packed. initialize
c           packing utilities for keeping statistics.
c
      pack2e=intpck.ne.0
      nhex=(-10* dlog10(cutoff))/12+1
      ntupl=-4
      if(pack2e) then
         if(dbug) write(iw,9996) nhex,ntupl
         if(jkskip.eq.0) then
            call initpk(rname)
            call initpn(rname)
         endif
      endif
c
c     ----- symmetry mapping of shells -----
c
chvd  niso=mxshel*(mxsym/2)/isingl
chvd  call daread(idafh,ioda,iso,niso,5)
chvd  do ii=1,nshell
chvd     do it=1,ntwd
chvd        indout(it)=iso(ii,it)
chvd     enddo
chvd     call isoout(nt)
chvd     do it=1,nt
chvd        iso(ii,it)=indin(it)
chvd     enddo
chvd  enddo
      nav = lenwrd()
      if (nt.gt.1) then
         call readi(iso,nw196(5)*nav,ibl196(5),idaf)
      endif
c
c     ----- define new order of shells according to type -----
c
      do ii=1,nshell
         modshl(ii)=ii
         invshl(ii)=ii
      enddo
c
      iinew=0
      do 270 iityp=1,ntyp
      do 260 ii=1,nshell
      lit=ktype(ii)
      sp =ktype(ii).eq.2.and.kmin(ii).eq.1
      if(sp) lit=ntyp
      if(lit.ne.iityp) go to 260
      do 210 it=1,nt
      id=mapshl(ii,it)
      if( id.gt.ii   ) go to 260
  210 mi(it)=id
c
      if(nt.eq.1) go to 240
      do 230 it=2,nt
      do 220 jt=1,it-1
      if(mi(jt).ne.mi(it)) go to 220
      mi(it)=0
  220 continue
  230 continue
  240 continue
c
c     ----- loop -250- may be activatived only if one is sure  -----
c           that all integral and derivative codes can handle
c           a modified shell order. beware ......
c     ----- the second derivative code does not work well yet  -----
c     ----- also for ordered supermatrices, the original shell -----
c           order must be used.
c
      if(.not.pk) then
         do 250 it=1,nt
            iiold=mi(nt+1-it)
            if(iiold.eq.0) go to 250
               iinew=iinew+1
c....          modshl(iinew)=iiold
c....          invshl(iiold)=iinew
  250       continue
      endif
c
  260 continue
  270 continue
c
      return
 9999 format(/,10x,27('-'),/,10x,'-giao- 2 electron integrals',
     1       /,10x,27('-'))
 9998 format(/,' the -pk- option is on , the integrals are in a',
     1         ' supermatrix form.',' -p- and -k- = ',l3)
 9997 format(/,' the -pk- option is off, the integrals are not in a',
     1         ' supermatrix form.')
 9996 format(  ' the integrals are packed with a.d.mclean',
     1         ' packing utilities (1977).',' nhex = ',i5,
     2         ' ntupl = ',i5,/)
 9995 format(/,' number of integrals per integral record = ',i8)
 9994 format(/,' threshold for keeping the integrals on file = ',e8.1)
 9993 format(/,' no integrals calculated in this run. -intskp- =',i4)
 9992 format(/,' integral evaluation method used = -pople- + -hondo-')
 9991 format(/,' integral evaluation method used = -hondo-          ')
 9990 format(/,' the pre-exponential threshold is too large for the',
     1         ' selected value of -cutoff- . ',/,
     2         ' their ratio must be less than 1.0e-06 .',/,
     3         ' increase -itol- in - $basis - . stop.',/,
     4         ' -dtol- , -cutoff- = ',2e20.12)
 9989 format(/,' -filtyp- and -pcktyp- = ',a8,2x,a8,' is not allowed.')
 9988 format(  ' incorrect label packing format. stop.',/,
     1         ' isingl,labsiz,num,numjk = ',4i5)
 9987 format(  ' incorrect label packing format. stop.',/,
     1         ' isingl,labsiz,num,numpk = ',4i5)
 9986 format(/,' integral index packing format is -labsiz- = ',i2)
 9985 format(/,' cost factors in semi-direct method -thime- = ',i10,
     1       /,'                                    -thize- = ',e10.4)
 9984 format(/,i8,' mbytes of disk reserved for the integral file(s),',
     1       /,   ' with up to ',i10,' record(s) of integrals',
     2            ' on each file.')
 9983 format(/,' - semi-direct - option is now in force.')
      end
      subroutine giapos
      implicit REAL (a-h,o-z)
      parameter (mxshel=MXSHEL)
      parameter (mxprim=MXPRIM)
      character*8 dsk4,dsk6 
      character*8 direct 
      logical pk,pandk
      logical block
      integer thime,semisp,dirstp
INCLUDE(../m4/common/vcore)
      common/sqfile/ijk,ipk
      common/restar1/nrest
      character*8 hndtyp,dsktyp,filtyp,pcktyp
      common/hndopt/hndtyp,dsktyp,filtyp,pcktyp
      common/intfil/nintmx
      common/intcal/intskp,notei
      common/intmem/iclb,inijg,ibuf,igint,iijklg,igijkl,ignkl,
     1 ignm,idij,idkl,ib00,ib01,ib10,ic00,id00,if00,
     2 isj,isk,isl,idijsi,idijsj,idklsk,idklsl,iabv,icv,irw,ichrg
      common/giao_shlint/cut,itotal,icount,
     1              irec,locint,istrt,jstrt,kstrt,lstrt
      common/pkfil/pk,pandk,block
      common/dirpar/thize,thime,semisp,dirstp
      common/nshel/
     1ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),cf(mxprim),cg(mxprim),
     1kstart(mxshel),katom(mxshel),ktype(mxshel),kng(mxshel),
     2kloc(mxshel),kmin(mxshel),kmax(mxshel),nshell
      data dsk4,dsk6 /'dsk-p -o','dsk-pk-o'/
      data direct /'direct  '/
c
      nrec  =irec
      intloc=locint
      ist   =istrt
      jst   =jstrt
      kst   =kstrt
      lst   =lstrt
      if(hndtyp.eq.direct.or.dirstp.eq.1) goto 10
         call rewfil(ijk)
         call rewfil(ipk)
   10 continue
c
c     ----- position integral file -----
c
      if(nrest.eq.0) go to 150
      if(nrec.le.1.or.intloc.le.1) go to 150
c
c     -----                  nrest = 1                   -----
c     ----- position the integral file for a restart job -----
c
  100 continue
      icount=intloc
      n=nrec-1
      if(n.eq.0) go to 120
         do 110 i=1,n
            if(filtyp.eq.dsk4.or.filtyp.eq.dsk6) then
               read(ijk)
               if(pandk) then
                  read(ipk)
               endif
            endif
  110       call advfil(ijk,ipk,pandk)
  120 continue
      if(filtyp.eq.dsk4.or.filtyp.eq.dsk6) then
         read(ijk) nxx
         if(pandk) then
            read(ipk) nxx
         endif
      endif
      if(.not.pandk)
     + call pread (ijk,    Q(ibuf),
     +             Q(ibuf+2*nintmx),
     +             nxx,nintmx)
      if(     pandk)
     + call pkread(ijk,ipk,Q(ibuf),
     +             Q(ibuf+  nintmx),
     +             Q(ibuf+2*nintmx),
     +             nxx,nintmx)
      call rewfil(ijk)
      if(n.eq.0) go to 140
         do 130 i=1,n
            if(filtyp.eq.dsk4.or.filtyp.eq.dsk6) then
               read(ijk)
               if(pandk) then
                  read(ipk)
               endif
            endif
  130       call advfil(ijk,ipk,pandk)
  140 continue
      go to 160
c
c     -----  nrest = 0             -----
c     ----- normal start           -----
c     ----- beware -intskp- option -----
c
  150 continue
      if(ist.le.1) ist=1
      if(jst.le.1) jst=1
      if(kst.le.1) kst=1
      if(lst.le.1) lst=1
c
      if(intskp.ne.0) ist=nshell+1
      if(ist.gt.nshell) go to 160
c
      if( (ist.ne.1).or.(jst.ne.1).or.
     1    (kst.ne.1).or.(lst.ne.1)    ) go to 100
      nrec=1
      intloc=1
  160 continue
      irec=nrec
      locint=intloc
      istrt=ist
      jstrt=jst
      kstrt=kst
      lstrt=lst
      icount=locint
c
      return
      end
      subroutine giamem(iaddr)
      implicit REAL (a-h,o-z)
      parameter (mxatom=MAXAT)
      parameter (mxshel=MXSHEL)
      parameter (mxprim=MXPRIM)
      parameter (lenvec=255)
      character*8 errmsg
      logical sp
      logical pk,pandk,block
      logical norm
      logical vector
      logical out 
      common/output/nprint
      common/iofile/ir,iw
      common/machnv/vector
      common/maxlen/maxcor
      common/baspar/normf,normp,itol
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,zan(mxatom),
     1                                         c(3,mxatom)
      common/nshel/
     1ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),cf(mxprim),cg(mxprim),
     1kstart(mxshel),katom(mxshel),ktype(mxshel),kng(mxshel),
     2kloc(mxshel),kmin(mxshel),kmax(mxshel),nshell
      common/intmem/iclb,inijg,ibuf,igint,iijklg,igijkl,ignkl,
     1 ignm,idij,idkl,ib00,ib01,ib10,ic00,id00,if00,
     2 isj,isk,isl,idijsi,idijsj,idklsk,idklsl,iabv,icv,irw,ichrg
      common/intopt/nopk,nok,nosqur
      common/intfil/nintmx
      common/inttyp/npople,nhondo
      common/intxyz/maxxyz,numxyz
      common/intvec/maxvec,minvec
      common/intpko/numpk
      common/pcklab/labsiz
      common/shlbas/maxtyp,maxnum
      common/shlnrm/pnrm(35)
      common/pkfil/pk,pandk,block
      common/shlgia/npwr,modtyp,modlit,modlkt
      dimension ijkln(5)
      dimension errmsg(3)
      data errmsg /'program ','stop in ','-giamem-'/
      data ijkln  /   1,  4, 10, 20, 35/
      data zero   /0.0d+00/
      data one    /1.00000000000000d+00/
      data sqrt3  /1.73205080756888d+00/
      data sqrt5  /2.23606797749979d+00/
      data sqrt7  /2.64575131106459d+00/
c
      out=.false.        
c
      maxvec=lenvec/3
      if(     vector) minvec=16
      if(.not.vector) minvec=maxvec*3+1
c
c     ----- check maximum angular momemtum -----
c
      sp=.false.
      maxtyp=0
      do i=1,nshell
         sp=sp.or.(ktype(i).eq.2.and.kmin(i).eq.1)
         if(ktype(i).gt.maxtyp) then
            maxtyp=ktype(i)
         endif
      enddo
c
      idum=maxtyp-1
      maxnum=((idum+1)*(idum+2))/2
      if(maxnum.eq.3.and.sp) maxnum=4
c
      minxyz=(4*maxtyp-2)/2
c
c     ----- -giao- parameters -----
c
      npwr=1
      modtyp=maxtyp+npwr
c
c     ----- check number of primitive charge distributions -----
c
      nijg=0
      do ii=1,nshell
         do jj=1,ii
            if(jj.ne.ii) nijg=nijg+ kng(ii)* kng(jj)
            if(jj.eq.ii) nijg=nijg+(kng(ii)*(kng(ii)+1))/2
         enddo
      enddo
c
c     -----  at this point it is good to remember that    -----
c            -maxtyp- = highest shell angular momentum
c            -maxfun- = number of functions with angular
c                       momentum less or equal to -maxtyp-
c            -maxnum- = number of functions with angular
c                       momentum         equal to -maxtyp-
c            -maxxyz- = maximum number of primitive integrals
c                       that can be handled in one -vector-
c            -numxyz- = actual maximum length of one -vector-
c            -maxxyz- = it is numxyz/3 . since the x, y, and z
c                       components are treated as a single vector,
c                       -maxxyz- is the number of (primitive-roots)
c                       combinations which can be treated in one
c                       vector. for -ssss- integrals which require
c                       one rys root, maxxyz happens to coincide with
c                       the number of primitive integrals treated in
c                       one vector. for -dddd- integrals which
c                       require five rys roots, the number of primitive
c                       integrals treated in one vector is -maxxyz-/5 .
c
c
c     ----- set normalization constants -----
c
      maxfun=ijkln(maxtyp)
      do i=1,maxfun
         pnrm(i)=one
      enddo
      norm=normf.ne.1.or.normp.ne.1
      if(norm) then      
         sqrt53=sqrt5/sqrt3
         do i=1,maxfun
            go to(110,110,160,160,110,160,160,120,160,160,
     1            110,160,160,130,160,160,160,160,160,120,
     2            110,160,160,140,160,160,160,160,160,150,
     3            160,160,120,160,160),i
  110       fi=one
            go to 160
  120       fi=sqrt3*fi
            go to 160
  130       fi=sqrt5*fi
            go to 160
  140       fi=sqrt7*fi
            go to 160
  150       fi=sqrt53*fi
  160       continue
            pnrm(i)=fi
         enddo
      endif
c
c     ----- size of gijkl =
c                         1   if s  shells
c                        81   if p  shells
c                      1296   if d  shells
c                     10000   if f  shells
c                     50625   if g  shells
c     ----- this version handles up to -g- shells   -----
c                as well as -sp- shells
c
      ngijkl=(maxnum**4)
      ngint =ngijkl*7
c
c     ----- calculate vector length and set memory pointers -----
c
      lfix=iaddr-1
      lfix=lfix+  (nshell*(nshell+1))/2
      lfix=lfix+( (nshell*(nshell+1))/2 )*2
      lfix=lfix+  nintmx*(1+6+1+labsiz)
      lfix=lfix+  ngint
      lfix=lfix+  ngijkl*4
      lfix=lfix+  nijg*13
      lvar=     ( modtyp**2       * modtyp**2       )*3*3
      lvar=lvar+( modtyp**2       *(modtyp+modtyp-1))*3
      lvar=lvar+((modtyp+modtyp-1)*(modtyp+modtyp-1))*3
      lvar=lvar+( modtyp**2                         )*3
      lvar=lvar+((modtyp+modtyp-1)                  )*3
      lvar=lvar+((modtyp+modtyp-1)* 3               )*3
      lvar=lvar+(  3                                )*3
      lvar=lvar+(  3 )
      lvar=lvar+(  4 )
      lvar=lvar+(  5 )
      lvar=lvar+( 18 )
      lvar=lvar+(  2 )
      maxxyz=(maxcor-lfix)/lvar
      if(maxxyz.lt.minxyz) then
         write(iw,9999) lvar,lfix,maxxyz,minxyz,maxcor
         call hnderr(3,errmsg)
      endif
      if(maxxyz.gt.maxvec) maxxyz=maxvec
      numxyz=3*maxxyz
c
c     Q(iclb  ) = coulomb integral threshold
c     Q(inijg ) = charge distribution pointers
c     Q(ibuf  ) = integral buffers
c     Q(igint ) = electron repulsion integrals
c     Q(iijklg) = indices
c     Q(ichrg ) = charge distribution parameters
c     Q(igijkl) = ( 2-d , 4 centers ) integrals ( 3 sets for -giao- )
c     Q(ignkl ) = ( 2-d , 3 centers ) integrals
c     Q(ignm  ) = ( 2-d , 2 centers ) integrals
c     Q(idij  ) = contraction density for -ij- charge distribution
c     Q(idkl  ) = contraction density for -kl- charge distribution
c     Q(ib00  ) = -b00-
c     Q(ib01  ) = -b01-
c     Q(ib10  ) = -b10-
c     Q(ic00  ) = -c00-
c     Q(id00  ) = -d00-
c     Q(if00  ) = -f00-
c     Q(isj   ) = temporary array when -sp- shells
c     Q(isk   ) = temporary array when -sp- shells
c     Q(isl   ) = temporary array when -sp- shells
c     Q(idijsi) = scaling factor for -s- function of an -sp- ii shell
c     Q(idijsj) = scaling factor for -s- function of an -sp- jj shell
c     Q(idklsk) = scaling factor for -s- function of an -sp- kk shell
c     Q(idklsl) = scaling factor for -s- function of an -sp- ll shell
c     Q(iabv  ) = -ab- vector for primitive integrals
c     Q(icv   ) = -cv- vector for primitive integrals
c     Q(irw   ) = -rw- vector for rys roots and weights
c
      iclb  =iaddr
      inijg =iclb  + (nshell*(nshell+1))/2
      ibuf  =inijg +((nshell*(nshell+1))/2)*2
      igint =ibuf  +  nintmx*(1+6+1+labsiz)
      iijklg=igint +  ngint     
      ichrg =iijklg+  ngijkl*4
      igijkl=ichrg +  nijg  *13
      ignkl =igijkl+( modtyp**2       * modtyp**2       )*maxxyz*3*3
      ignm  =ignkl +( modtyp**2       *(modtyp+modtyp-1))*maxxyz*3
      idij  =ignm  +((modtyp+modtyp-1)*(modtyp+modtyp-1))*maxxyz*3
      idkl  =idij  +( modtyp**2                         )*maxxyz*3
      ib00  =idkl  +((modtyp+modtyp-1)                  )*maxxyz*3
      ib01  =ib00  +((modtyp+modtyp-1)                  )*maxxyz*3
      ib10  =ib01  +((modtyp+modtyp-1)                  )*maxxyz*3
      ic00  =ib10  +((modtyp+modtyp-1)                  )*maxxyz*3
      id00  =ic00  +(  1                                )*maxxyz*3
      if00  =id00  +(  1                                )*maxxyz*3
      isj   =if00  +(  1                                )*maxxyz*3
      isk   =isj   +(  1                                )*maxxyz
      isl   =isk   +(  1                                )*maxxyz
      idijsi=isl   +(  1                                )*maxxyz
      idijsj=idijsi+(  1                                )*maxxyz
      idklsk=idijsj+(  1                                )*maxxyz
      idklsl=idklsk+(  1                                )*maxxyz
      iabv  =idklsl+(  1                                )*maxxyz
      icv   =iabv  +(  5                                )*maxxyz
      irw   =icv   +( 18                                )*maxxyz
      ilast =irw   +(  2                                )*maxxyz
      ineed =ilast -1
c
      if(out) then
         write(iw,9998) numxyz,ineed
      endif
      call fiflsh(iw)
c
      return
 9999 format(/,' not enough memory for -gia2ei-',/,
     1         ' lvar,lfix,maxxyz,minxyz,maxcor = ',5i10)
 9998 format(/,' - int  - , maxvec = ',i3,' ineed = ',i10,/)
      end
      subroutine giahnd(gclb,nijg)
      implicit REAL (a-h,o-z)
      parameter (mxatom=MAXAT)
      parameter (mxshel=MXSHEL)
      parameter (mxprim=MXPRIM)
      parameter (maxorb=MAXORB)
      parameter (mxsym=48)
      parameter (mxgsh =30,mxgsh2=mxgsh*mxgsh)
      logical sym
      logical semi
      logical disk
      logical some
      logical out
      logical dbug
      logical lcap
      logical sp
      logical hndskp
      logical ieqj0,keql0,same0,first
      logical iieqjj,kkeqll,ijeqkl,ijgtkl,ijltkl
      logical outjk
      logical outgia
INCLUDE(../m4/common/vcore)
      common/lcapid/nap,iap
      common/iofile/ir,iw
      common/symtry/invt(mxsym),nt,ntmax,ntwd,nosym
      common/symshl/mapshl(mxshel,mxsym),mapctr(mxatom,mxsym)
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,zan(mxatom),
     1                                         c(3,mxatom)
      common/nshel/
     1ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),cf(mxprim),cg(mxprim),
     1kstart(mxshel),katom(mxshel),ktype(mxshel),kng(mxshel),
     2kloc(mxshel),kmin(mxshel),kmax(mxshel),nshell
      common/shlord/modshl(mxshel),invshl(mxshel)
      common/shlold/numi0,numj0,numk0,numl0,
     1              ieqj0,keql0,same0,first
      common/shlequ/iieqjj,kkeqll,ijeqkl,ijgtkl,ijltkl
      common/shlfac/qq4
      common/shlbas/maxtyp,maxnum
      common/giao_shlint/cutoff,itotal,icount,irec,intloc,
     1                   ist,jst,kst,lst
      common/ijchrg/acharg(11,mxgsh2),dasi(mxgsh2),dasj(mxgsh2),nij
      common/klchrg/bcharg(11,mxgsh2),dbsk(mxgsh2),dbsl(mxgsh2),nkl
      common/intfil/nintmx
      common/intmem/iclb,inijg,ibuf,igint,iijklg,igijkl,ignkl,
     1 ignm,idij,idkl,ib00,ib01,ib10,ic00,id00,if00,
     2 isj,isk,isl,idijsi,idijsj,idklsk,idklsl,iabv,icv,irw,ichrg
      common/intxyz/maxxyz,numxyz
chvd
c     This common block was missing. (iexch set but not passed on)
      common/intxch/norgin(3),iexch
chvd
      common/intdsk/maxrec,ijklnu,ijklsp,ijklno,nspdsk,nodisk
      common/timez/timlim,ti,tx,tim,wti,wtx,wtim
      common/jkstat/lvec,nvec
INCLUDE(../m4/common/mapper)
      dimension mi(mxsym),mj(mxsym),mk(mxsym),m0(mxsym)
      dimension gclb(*),nijg(2,*)
      data zero   /0.0d+00/
c
      lvec=0
      nvec=0
c
c     ----- set some parameters -----
c
      sym=nt.gt.1
      lcap=nap.gt.1
c
      dbug=.false.
      out =.false.
      out =out.or.dbug
      some=.false.         
      some=some.or.out 
      if(out) then
         write(iw,9996)
      endif
      call fiflsh(iw)
c
      call cpuwal(tim,wtim)
      tim0=tim
      tim1=tim
c
c     ----- get one-electron charge distribution -----
c
      call oehnd(nijg)
c
c     ----- two-electron integrals -----
c
      numi0=0
      numj0=0
      numk0=0
      numl0=0
      ieqj0=.false.
      keql0=.false.
      same0=.false.
      first=.true.
c
      if(out) then
         write(iw,9997)
      endif
c
c     ----- ishell -----
c
      ijklsh=0
      ijklnu=0
      do 9000 ii=1,nshell
      iimod=modshl(ii)
c
      intrec=irec
      intloc=icount
c
c     ----- jshell -----
c
      do 8000 jj=1,ii
      jjmod=modshl(jj)
c
c     ----- get -ij- charge distribution -----
c
      iijj=iky(max(iimod,jjmod))+min(iimod,jjmod)
      nij0=nijg(1,iijj)
      nij =nijg(2,iijj)
      if(nij.eq.0) go to 8000
      call oeread(Q(ichrg),nij0,acharg,dasi,dasj,nij)
c
c     ----- kshell -----
c
      do 7000 kk=1,nshell
      kkmod=modshl(kk)
c
c     ----- lshell ----
c
      do 6000 ll=1,kk       
      llmod=modshl(ll)
c
c     ----- get -kl- charge distribution -----
c
      kkll=iky(max(kkmod,llmod))+min(kkmod,llmod)
      nkl0=nijg(1,kkll)
      nkl =nijg(2,kkll)
      if(nkl.eq.0) go to 6000
      call oeread(Q(ichrg),nkl0,bcharg,dbsk,dbsl,nkl)
c
      if(sym) then
         if(ii.gt.kk.or.(ii.eq.kk.and.jj.ge.ll)) then
            ii1 = ii
            jj1 = jj
            kk1 = kk
            ll1 = ll
         else
            ii1 = kk
            jj1 = ll
            kk1 = ii
            ll1 = jj
         endif
c
         n4=0
         do 110 it=1,nt
         id=mapshl(iimod,it)
         id=invshl(id)
         if(id.gt.ii1) go to 6000
         jd=mapshl(jjmod,it)
         jd=invshl(jd)
         if(jd.gt.ii1) go to 6000
         if(jd.gt.id) then
            nd=id
            id=jd
            jd=nd
         endif
         kd=mapshl(kkmod,it)
         kd=invshl(kd)
         if(kd.gt.ii1) go to 6000
         ld=mapshl(llmod,it)
         ld=invshl(ld)
         if(ld.gt.ii1) go to 6000
         if(ld.gt.kd) then
            nd=kd
            kd=ld
            ld=nd
         endif
         if(id.ne.ii1.and.kd.ne.ii1) go to 110
         if(kd.gt.id.or.(kd.eq.id.and.ld.gt.jd)) then
            nd=kd
            kd=id
            id=nd
            nd=ld
            ld=jd
            jd=nd
         endif
         if(jd.gt.jj1) go to 6000
         if(jd.lt.jj1) go to 110
         if(kd.gt.kk1) go to 6000
         if(kd.lt.kk1) go to 110
         if(ld.gt.ll1) go to 6000
         if(ld.lt.ll1) go to 110
         n4=n4+1
  110    continue
      else
         n4=1
      endif
      q4= dble(nt)/ dble(n4)
      qq4=q4
c
c     ----- select shell block if parallel processing -----
c
      ijklsh=ijklsh+1
      if(mod(ijklsh,nap)+1.ne.iap) go to 6000
c
c     ----- check integral bounds and skip conditions for semi-d -----
c
      ijij  =iky(max(iimod,jjmod))+min(iimod,jjmod)
      klkl  =iky(max(kkmod,llmod))+min(kkmod,llmod)
      shize =sqrt(gclb(ijij)*gclb(klkl))
      hndskp=shize.lt.cutoff
      if(hndskp) go to 6000
c
c     ----- get shell parameters and set addressing arrays -----
c
c
      if(dbug) then
c
c     ----- (ii,jj//kk,ll) -----
c
         if((kk.lt.ii).or.(kk.eq.ii.and.ll.le.jj)) then
            iexch=1
            ish=iimod
            jsh=jjmod
            ksh=kkmod
            lsh=llmod
chvd        call jkshel(ish,jsh,ksh,lsh)
chvd        call jkindx(Q(iijklg))
chvd        call jkspdf(Q(igint),Q(iijklg),
chvd 1                  Q(igijkl),Q(ignkl),Q(ignm),
chvd 2                  Q(ib00),Q(ib01),Q(ib10),
chvd +                  Q(ic00),Q(id00),Q(if00),
chvd 3                  Q(idij),Q(idkl),Q(isj ),
chvd +                  Q(isk ),Q(isl ),
chvd 4                  Q(idijsi),Q(idijsj),
chvd +                  Q(idklsk),Q(idklsl),
chvd 5                  Q(iabv),Q(icv ),Q(irw   ),
chvd 6                  maxxyz,outjk)
            if(outjk) then
chvd           call jkout(Q(ibuf),Q(ibuf+nintmx),
chvd +                    Q(igint),
chvd 1                    Q(iijklg),disk)
            endif
         endif
      else
c
c     ----- -giao- (ii,jj//kk,ll) -----
c
         iexch=1
         ish=iimod
         jsh=jjmod
         ksh=kkmod
         lsh=llmod
         call giashl(ish,jsh,ksh,lsh)
         ijeqkl = .false.
         call giaind(Q(iijklg))
         if(out) then
            write(iw,*) ' ... calling -giaspd- '
            write(iw,*) 'ish,jsh,ksh,lsh = ',ish,jsh,ksh,lsh
         endif
         call giaspd(Q(igint),Q(iijklg),
     1               Q(igijkl),Q(ignkl),Q(ignm),
     2               Q(ib00),Q(ib01),Q(ib10),
     +               Q(ic00),Q(id00),Q(if00),
     3               Q(idij),Q(idkl),Q(isj ),
     +               Q(isk ),Q(isl ),
     4               Q(idijsi),Q(idijsj),
     +               Q(idklsk),Q(idklsl),
     5               Q(iabv),Q(icv ),Q(irw   ),
     6               maxxyz,outgia)
         if(out) then
            write(iw,*) ' ... in  to outgia = ',outgia
            call fiflsh(iw)
         endif
         if(outgia) then
            call giaout(Q(ibuf),Q(ibuf+nintmx),
     +                  Q(igint),
     1                  Q(iijklg),disk)
         endif
         if(out) then
            write(iw,*) ' ... out of outgia = ',outgia
            call fiflsh(iw)
         endif
c
      endif
 6000 continue
 7000 continue
 8000 continue
c
c     ----- check cpu time ; exit if needed -----
c
      ncall=0
      irest=1
      call texit(ncall,irest)
      if(tim.ge.timlim) then       
         if(dbug) then
chvd        call jkend(0,ii,jj,kk,ll,
chvd 1                 Q(ibuf),
chvd +                 Q(ibuf+nintmx),
chvd +                 Q(ibuf+2*nintmx),disk)
         else
            call giaend(0,ii,jj,kk,ll,
     1                  Q(ibuf),
     +                  Q(ibuf+nintmx),
     +                  Q(ibuf+2*nintmx),disk)
         endif
         go to 9100
      endif
c
c     ----- print intermediate (restart) data -----
c
      call cpuwal(tim,wtim)
      dt0=tim-tim0
      dt1=tim-tim1
      tim1=tim
      if(out) write(iw,9998) ii,iimod,iimod,iimod,iimod,
     1                       intrec,intloc,dt1,dt0
 9000 continue
c
      if(dbug) then
chvd     call jkend(1,ii,jj,kk,ll,
chvd 1              Q(ibuf),
chvd +              Q(ibuf+nintmx),
chvd +              Q(ibuf+2*nintmx),disk)
      else
         call giaend(1,ii,jj,kk,ll,
     1               Q(ibuf),
     +               Q(ibuf+nintmx),
     +               Q(ibuf+2*nintmx),disk)
      endif
 9100 continue
c
      avgvec=zero
      if(nvec.ne.0) avgvec=dble(lvec)/dble(nvec)
      if(out) write(iw,9995) nvec,avgvec
c
      call timit(0)
      dtim=tim-tim0
      if(lcap.and.some) write(iw,9999) dtim,tim
      return
 9999 format(' -giahnd- elapsed - cpu- time = ',f10.3,
     1                5x,'total - cpu- time = ',f10.3)
 9998 format(i4,4i4,i10,i6,f15.3,f15.3)
 9997 format('     ish jsh ksh lsh      irec intloc      del(time)  ',
     1       '        time',/,1x,66('-'))
 9996 format(1x,'-giaint(hondo)-',/,1x,15('-'))
 9995 format(' # of vectors = ',i10,' average vector length = ',f10.1)
      end
      subroutine giaspd(gint,ijklg,gijkl,gnkl,gnm,
     1 b00,b01,b10,c00,d00,f00,dij,dkl,
     2 sj,sk,sl,dijsi,dijsj,dklsk,dklsl,abv,cv,rwv,maxxyz,
     3 outgia)
      implicit REAL (a-h,o-z)
      parameter (mxgsh =30,mxgsh2=mxgsh*mxgsh)
      logical nmaxs,nmaxp,mmaxs,mmaxp
      logical expndi,expndk
      logical spi,spj,spk,spl,spij,spkl,spijkl
      logical iieqjj,kkeqll,ijeqkl,ijgtkl,ijltkl
      logical pk,pandk,block
      logical first
      logical last
      logical outgia
      common/jkstat/lvec,nvec
      common/ijchrg/acharg(11,mxgsh2),dasi(mxgsh2),dasj(mxgsh2),nij
      common/klchrg/bcharg(11,mxgsh2),dbsk(mxgsh2),dbsl(mxgsh2),nkl
      common/shltyp/spi,spj,spk,spl,spij,spkl,spijkl
      common/shlfac/qq4
      common/shlgnm/nmaxs,nmaxp,mmaxs,mmaxp
      common/shlxpn/expndi,expndk
      common/shlpar/lit,ljt,lkt,llt,loci,locj,lock,locl,
     1              mini,minj,mink,minl,maxi,maxj,maxk,maxl
      common/shlnum/numi,numj,numk,numl,ij,kl,ijkl
      common/shltol/rtol,dtol
      common/shlequ/iieqjj,kkeqll,ijeqkl,ijgtkl,ijltkl
      common/shlgia/npwr,modtyp,modlit,modlkt
      common/intvec/maxvec,minvec
      common/pkfil/pk,pandk,block
      common/hfk/xx,u(9),w(9),nroots
      common/root/yy,t(12),v(12),mroots
      dimension gint(*),ijklg(*),gijkl(*),gnkl(*),gnm(*)
      dimension b00(*),b01(*),b10(*),c00(*),d00(*),f00(*)
      dimension sj(*),sk(*),sl(*),dijsi(*),dijsj(*),dklsk(*),dklsl(*)
      dimension abv(5,*),cv(18,*)
      dimension rwv(2,*)
      data pi252   /34.986836655250d+00/
      data pt5,one /0.5d+00,1.0d+00/
c
      outgia=.false.
c
      nimax=modlit
      njmax=   ljt
      nkmax=modlkt
      nlmax=   llt
      nmax =modlit+ljt-1
      mmax =modlkt+llt-1
      nmaxs=nmax.eq.1
      nmaxp=nmax.le.2
      mmaxs=mmax.eq.1
      mmaxp=mmax.le.2
c
      maxg=maxxyz/nroots
      dtol2=dtol*dtol
      q4=pi252*qq4
c
c     ----- pair of k,l primitives -----
c
      first=.true.
      ng=0
      klg=0
  100 klg=klg+1
      if(klg.gt.nkl) go to 300
      db=bcharg( 1,klg)
      bb=bcharg( 2,klg)
      xb=bcharg( 3,klg)
      yb=bcharg( 4,klg)
      zb=bcharg( 5,klg)
      xd=bcharg( 6,klg)
      yd=bcharg( 7,klg)
      zd=bcharg( 8,klg)
      dxkl=bcharg( 9,klg)
      dykl=bcharg(10,klg)
      dzkl=bcharg(11,klg)
      q4db=q4*db
c
c     ----- pair of i,j primitives -----
c
      ijg=0
  200 ijg=ijg+1
      if(ijg.gt.nij) go to 100
      da=acharg( 1,ijg)
      aa=acharg( 2,ijg)
      xa=acharg( 3,ijg)
      ya=acharg( 4,ijg)
      za=acharg( 5,ijg)
      aandb1=one/(aa+bb)
      q4dbda=q4db*da
      dum   =q4dbda*q4dbda*aandb1
      if(dum.le.dtol2) go to 200
      q4dbda=q4dbda* sqrt(aandb1)
      rho   =aa*bb*aandb1
      xx=rho*((xa-xb)**2+(ya-yb)**2+(za-zb)**2)
c
      ng=ng+1
      abv(1,ng)=aa
      abv(2,ng)=bb
      abv(3,ng)=rho
      abv(4,ng)=q4dbda
      abv(5,ng)=xx
c
      xc=acharg( 6,ijg)
      yc=acharg( 7,ijg)
      zc=acharg( 8,ijg)
      dxij=acharg( 9,ijg)
      dyij=acharg(10,ijg)
      dzij=acharg(11,ijg)
c
      if(mmaxs) go to 210
      cv( 1,ng)=aa*(xa-xd)
      cv( 2,ng)=bb*(xb-xd)
      cv( 3,ng)=aa*(ya-yd)
      cv( 4,ng)=bb*(yb-yd)
      cv( 5,ng)=aa*(za-zd)
      cv( 6,ng)=bb*(zb-zd)
  210 if(nmaxs) go to 220
      cv( 7,ng)=aa*(xa-xc)
      cv( 8,ng)=bb*(xb-xc)
      cv( 9,ng)=aa*(ya-yc)
      cv(10,ng)=bb*(yb-yc)
      cv(11,ng)=aa*(za-zc)
      cv(12,ng)=bb*(zb-zc)
  220 continue
      cv(13,ng)=dxij
      cv(14,ng)=dyij
      cv(15,ng)=dzij
      cv(16,ng)=dxkl
      cv(17,ng)=dykl
      cv(18,ng)=dzkl
      if(spi) dijsi(ng)=dasi(ijg)
      if(spj) dijsj(ng)=dasj(ijg)
      if(spk) dklsk(ng)=dbsk(klg)
      if(spl) dklsl(ng)=dbsl(klg)
c
      if(ng.lt.maxg) go to 200
      last=.false.
      go to 310
c
  300 continue
      last=.true.
  310 continue
      numg=ng
      if(numg.eq.0) go to 700
c
      if(nroots.eq.1) go to 480
      if(.not.spi) go to 420
      do 410 iroot=1,nroots
      do 410 ig=1,numg
      dijsi(ig+numg*(iroot-1))=dijsi(ig)
  410 continue
  420 if(.not.spj) go to 440
      do 430 iroot=1,nroots
      do 430 ig=1,numg
      dijsj(ig+numg*(iroot-1))=dijsj(ig)
  430 continue
  440 if(.not.spk) go to 460
      do 450 iroot=1,nroots
      do 450 ig=1,numg
      dklsk(ig+numg*(iroot-1))=dklsk(ig)
  450 continue
  460 if(.not.spl) go to 480
      do 470 iroot=1,nroots
      do 470 ig=1,numg
      dklsl(ig+numg*(iroot-1))=dklsl(ig)
  470 continue
  480 continue
c
      lvec=lvec+numg*nroots*3
      nvec=nvec+1
c
c     ----- compute roots and weights for quadrature -----
c
      call jkwrys(rwv,abv,numg)
c
c     ----- compute coefficients for recursion formulae -----
c
      call jkbcdf(b00,b01,b10,c00,d00,f00,dij,dkl,
     1            abv,cv,rwv,numg,nroots)
c
c     ----- compute -x- , -y- , -z- integrals ( 2 centers,2-d ) -----
c
      call giagnm(gnm,numg*nroots*3,nmax,mmax,
     1            b00,b01,b10,c00,d00,f00)
c
c     ----- compute -x- , -y- , -z- integrals ( 4 centers,2-d ) -----
c
      call giaxyz(gijkl,gijkl,gnkl,gnkl,gnkl,gnm,gnm,
     1            numg*nroots*3,nmax,mmax,nimax,njmax,nkmax,nlmax,
     2            dij,dkl,expndi,expndk)
c
c     ----- -giao- specific -----
c
      call xyzgia(gijkl,numg*nroots,nimax,njmax,nkmax,nlmax)
c
c     ----- zero out first time around -----
c
      if(.not.first) go to 600
      call giazer(gint,ijklg)
      first=.false.
  600 continue
c
c     ----- compute  - (i,j/k,l) -  integrals -----
c
      call spdgia(numg,nroots,gint,ijklg,gijkl,
     1            sj,sk,sl,dijsi,dijsj,dklsk,dklsl)
c
      if(last) go to 700
      ng=0
      go to 200
  700 continue
      if(numg.ne.0.or.(.not.first)) go to 2000
      call giazer(gint,ijklg)
      return
c
 2000 continue
c
c     ----- ready for -giao- combinations -----
c
      call giamak(gint,ijklg)
c
c     ----- ready to write out -----
c
      outgia=.true.
      return
      end
      subroutine giaout(ix,xx,gijkl,ijklg,disk)
      implicit REAL (a-h,o-z)
      parameter (mxatom=MAXAT)
      parameter (maxorb=MAXORB)
      character*8 direct
      integer shiftl
      logical keqi
      logical iieqjj,kkeqll,ijeqkl,ijgtkl,ijltkl
      logical out
      logical disk
      logical dbllab
      logical keep 
      common/iofile/ir,iw,ip
      common/sqfile/ijk,ipk
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,zan(mxatom),
     1                                         c(3,mxatom)
      common/giao_shlint/cutoff,itotal,icount,irec,intloc,
     1                   ist,jst,kst,lst
      common/shlfac/qq4
      common/shlpar/lit,ljt,lkt,llt,loci,locj,lock,locl,
     1              mini,minj,mink,minl,maxi,maxj,maxk,maxl
      common/shllmn/igxyz(4,35),jgxyz(4,35),kgxyz(4,35),lgxyz(4,35)
      common/shlnrm/pnrm(35)
      common/shlequ/iieqjj,kkeqll,ijeqkl,ijgtkl,ijltkl
      common/intfil/nintmx
      common/intprt/qint(3),valint(3),jcint(16),out
      common/pcklab/labsiz
      dimension xx(7,*),ix(*),gijkl(7,*),ijklg(4,*)
      data direct /'direct  '/
      data pt5    /0.5d+00/
c
      shiftl(iarg,ibit)= ishft(iarg, ibit)
      lor(iarg1,iarg2) = ior(iarg1,iarg2)
      iword(iarg1,iarg2,iarg3,iarg4)=
     1                      lor(shiftl(iarg1,24),lor(shiftl(iarg2,16),
     2                      lor(shiftl(iarg3, 8),    shiftl(iarg4, 0))))
      jword(iarg1,iarg2)= lor(shiftl(iarg1,16),shiftl(iarg2, 0))
c
      dbllab=labsiz.eq.2
      num2=(num *(num +1))/2
      num4=(num2*(num2+1))/2
c
c     ----- pack the 4 indices of integrals into one word -----
c     ----- write labels + integrals on tape (ijk)        -----
c
c     ----- general case -----
c
      ijkln=0
      do i=mini,maxi
      fi=pnrm(i)
      jmax=maxj
      if(iieqjj) jmax=i
      do j=minj,jmax
      fij=fi*pnrm(j)
      ii=max((loci+i),(locj+j))
      jj=   ((loci+i)+(locj+j))-ii
      do k=mink,maxk
      fijk=fij*pnrm(k)
      lmax=maxl
      if(kkeqll) lmax=k
      do l=minl,lmax
      fijkl=fijk*pnrm(l)
c
      kk=max((lock+k),(locl+l))
      ll=   ((lock+k)+(locl+l))-kk
c
      ijkln=ijkln+1
      g  =gijkl(1,ijklg(1,ijkln))*fijkl
      gx1=gijkl(2,ijklg(1,ijkln))*fijkl
      gy1=gijkl(3,ijklg(1,ijkln))*fijkl
      gz1=gijkl(4,ijklg(1,ijkln))*fijkl
      gx2=gijkl(5,ijklg(1,ijkln))*fijkl
      gy2=gijkl(6,ijklg(1,ijkln))*fijkl
      gz2=gijkl(7,ijklg(1,ijkln))*fijkl
      keep= 
     1   abs(g  ).ge.cutoff.or.
     2   abs(gx1).ge.cutoff.or.
     3   abs(gy1).ge.cutoff.or.
     4   abs(gz1).ge.cutoff.or.
     5   abs(gx2).ge.cutoff.or.
     6   abs(gy2).ge.cutoff.or.
     7   abs(gz2).ge.cutoff   
      if(keep) then      
         if(iieqjj.and.j.eq.i) then
            g  =g  *pt5
            gx1=gx1*pt5
            gy1=gy1*pt5
            gz1=gz1*pt5
            gx2=gx2*pt5
            gy2=gy2*pt5
            gz2=gz2*pt5
         endif
         if(kkeqll.and.l.eq.k) then
            g  =g  *pt5
            gx1=gx1*pt5
            gy1=gy1*pt5
            gz1=gz1*pt5
            gx2=gx2*pt5
            gy2=gy2*pt5
            gz2=gz2*pt5
         endif
c
         if(out) then       
            n0=ijklg(1,ijkln)
            g0=gijkl(1,n0)*pnrm(i)*pnrm(j)*pnrm(k)*pnrm(l)
            g1=gijkl(2,n0)*pnrm(i)*pnrm(j)*pnrm(k)*pnrm(l)
            g2=gijkl(3,n0)*pnrm(i)*pnrm(j)*pnrm(k)*pnrm(l)
            g3=gijkl(4,n0)*pnrm(i)*pnrm(j)*pnrm(k)*pnrm(l)
            g4=gijkl(5,n0)*pnrm(i)*pnrm(j)*pnrm(k)*pnrm(l)
            g5=gijkl(6,n0)*pnrm(i)*pnrm(j)*pnrm(k)*pnrm(l)
            g6=gijkl(7,n0)*pnrm(i)*pnrm(j)*pnrm(k)*pnrm(l)
            call intout(ii,jj,kk,ll,qq4,n0,g0,0)
            call intout(ii,jj,kk,ll,qq4,n0,g1,0)
            call intout(ii,jj,kk,ll,qq4,n0,g2,0)
            call intout(ii,jj,kk,ll,qq4,n0,g3,0)
            call intout(ii,jj,kk,ll,qq4,n0,g4,0)
            call intout(ii,jj,kk,ll,qq4,n0,g5,0)
            call intout(ii,jj,kk,ll,qq4,n0,g6,0)
         endif
c
         if(dbllab) then
            xx(1, icount            )=g
            xx(2, icount            )=gx1
            xx(3, icount            )=gy1
            xx(4, icount            )=gz1
            xx(5, icount            )=gx2
            xx(6, icount            )=gy2
            xx(7, icount            )=gz2
            ix(  (icount-1)*labsiz+1)=jword(ii,jj)
            ix(  (icount-1)*labsiz+2)=jword(kk,ll)
            icount=icount+1
            if(icount.le.nintmx) go to 100
            nxx=nintmx
            call giaf10(    xx,ix,nxx,nintmx)
            itotal=itotal+nintmx
            irec  =irec  +1
            icount=1
            go to 100
         else
            xx(1,icount)=g
            xx(2,icount)=gx1
            xx(3,icount)=gy1
            xx(4,icount)=gz1
            xx(5,icount)=gx2
            xx(6,icount)=gy2
            xx(7,icount)=gz2
            ix(  icount)=iword(ii,jj,kk,ll)
            icount=icount+1
            if(icount.le.nintmx) go to 100
            nxx=nintmx
            call giaf10(    xx,ix,nxx,nintmx)
            itotal=itotal+nintmx
            irec  =irec  +1
            icount=1
            go to 100
         endif
c
      endif
  100 continue
      enddo
      enddo
      enddo
      enddo
c
      if(out) then
         call intout(ii,jj,kk,ll,qq4,n0,g0,1)
      endif
c
      return
      end
      subroutine giaf10(xx,ix,nxx,nintmx)
      implicit REAL (a-h,o-z)
      logical  out
      common/dirgia/ifck0,ifckx,ifcky,ifckz,idens,idim
INCLUDE(../m4/common/vcore)
      dimension xx(*),ix(*)
c
c ... idens - triangular density matrix
c ... ifck0 - square     fock    matrix
c ... ifckx - square     fock    matrix
c ... ifcky - square     fock    matrix
c ... ifckz - square     fock    matrix
c
c     ----- square fock-like matrix for -giao- -----
c
      out=.false.
      if(out) then
         call prtr(Q(idens),idim)
         call prsq(Q(ifck0),idim,idim,idim)
         call prsq(Q(ifckx),idim,idim,idim)
         call prsq(Q(ifcky),idim,idim,idim)
         call prsq(Q(ifckz),idim,idim,idim)
      endif
c
      call giafck(Q(ifck0),Q(ifckx),
     +            Q(ifcky),Q(ifckz),
     1            Q(idens),idim,xx,ix,abs(nxx))
c
      if(out) then
         call prsq(Q(ifck0),idim,idim,idim)
         call prsq(Q(ifckx),idim,idim,idim)
         call prsq(Q(ifcky),idim,idim,idim)
         call prsq(Q(ifckz),idim,idim,idim)
      endif
      return
      end
      subroutine giafck(f,fx,fy,fz,d,ndim,xx,ix,nint)
      implicit REAL (a-h,o-z)
      parameter (maxorb=MAXORB)
      parameter (pt5=0.5d+00)
      parameter (two=2.0d+00)
      logical dbllab
      logical out
      integer shiftr
      character*8 errmsg
      common/iofile/ir,iw
      common/pcklab/labsiz
INCLUDE(../m4/common/mapper)
      dimension f(ndim,*),fx(ndim,*),fy(ndim,*),fz(ndim,*)
      dimension d(*)
      dimension xx(7,*),ix(*)
      dimension errmsg(3)    
      data errmsg /'program ','stop in ','-giafck-'/
c
      shiftr(iarg,ibit)= ishft(iarg,-ibit)
      land(iarg1,iarg2)= iand(iarg1,iarg2)
      mask1(iarg)=2**iarg-1
      mask2(iarg)=2**iarg-1
      iword1(iarg)=land(shiftr(iarg,24),mask1(8))
      iword2(iarg)=land(shiftr(iarg,16),mask1(8))
      iword3(iarg)=land(shiftr(iarg, 8),mask1(8))
      iword4(iarg)=land(shiftr(iarg, 0),mask1(8))
      jword1(iarg)=land(shiftr(iarg,16),mask1(16))
      jword2(iarg)=land(shiftr(iarg, 0),mask1(16))
c
      out=.false.
c
      dbllab=labsiz.eq.2
c
      do m=1,nint
         if(dbllab) then
            label1=ix((m-1)*labsiz+1)
            label2=ix((m-1)*labsiz+2)
            i=jword1(label1)
            j=jword2(label1)
            k=jword1(label2)
            l=jword2(label2)
         else
            label=ix(m)
            i=iword1(label)
            j=iword2(label)
            k=iword3(label)
            l=iword4(label)
         endif
         ij=iky(max(i,j))+min(i,j)
         kl=iky(max(k,l))+min(k,l)
         ik=iky(max(i,k))+min(i,k) 
         il=iky(max(i,l))+min(i,l) 
         jk=iky(max(j,k))+min(j,k)
         jl=iky(max(j,l))+min(j,l)
c
         dij=d(ij)
         dkl=d(kl)
         dik=d(ik)
         dil=d(il)
         djk=d(jk)
         djl=d(jl)
c
         val  =xx(1,m)
         valx1=xx(2,m)
         valy1=xx(3,m)
         valz1=xx(4,m)
         valx2=xx(5,m)
         valy2=xx(6,m)
         valz2=xx(7,m)
         valh  = val  *pt5
         valx1h= valx1*pt5
         valy1h= valy1*pt5
         valz1h= valz1*pt5
         valx2h= valx2*pt5
         valy2h= valy2*pt5
         valz2h= valz2*pt5
         if(k.eq.l) then
            valx1c=valx1*two
            valy1c=valy1*two
            valz1c=valz1*two
         else
            valx1c=valx1*two
            valy1c=valy1*two
            valz1c=valz1*two
         endif
         if(out) then
            write(iw,9999) i,j,k,l,
     1                     dij,dkl,dik,djl,dil,djk,
     1                     val,
     1                     valx1,valy1,valz1,valx2,valy2,valz2
         endif
c
c     ----- f00   -----
c
         f(i,j)  = f(i,j) +dkl*( val   +val   )
         f(j,i)  = f(j,i) +dkl*( val   +val   )
         f(i,k)  = f(i,k) -djl*  valh
         f(i,l)  = f(i,l) -djk*  valh
         f(j,k)  = f(j,k) -dil*  valh
         f(j,l)  = f(j,l) -dik*  valh
c
c     ----- f10-x -----
c
         fx(i,j) = fx(i,j)+dkl*(+valx1c       )
         fx(j,i) = fx(j,i)+dkl*(-valx1c       )
         fx(i,k) = fx(i,k)-djl*(+valx1h-valx2h)
         fx(i,l) = fx(i,l)-djk*(+valx1h+valx2h)
         fx(j,k) = fx(j,k)-dil*(-valx1h-valx2h)
         fx(j,l) = fx(j,l)-dik*(-valx1h+valx2h)
c
c     ----- f10-y -----
c
         fy(i,j) = fy(i,j)+dkl*(+valy1c       )
         fy(j,i) = fy(j,i)+dkl*(-valy1c       )
         fy(i,k) = fy(i,k)-djl*(+valy1h-valy2h)
         fy(i,l) = fy(i,l)-djk*(+valy1h+valy2h)
         fy(j,k) = fy(j,k)-dil*(-valy1h-valy2h)
         fy(j,l) = fy(j,l)-dik*(-valy1h+valy2h)
c
c     ----- f10-z -----
c
         fz(i,j) = fz(i,j)+dkl*(+valz1c       )
         fz(j,i) = fz(j,i)+dkl*(-valz1c       )
         fz(i,k) = fz(i,k)-djl*(+valz1h-valz2h)
         fz(i,l) = fz(i,l)-djk*(+valz1h+valz2h)
         fz(j,k) = fz(j,k)-dil*(-valz1h-valz2h)
         fz(j,l) = fz(j,l)-dik*(-valz1h+valz2h)
      enddo       
c
      return
 9999 format(' -giafck- i,j,k,l = ',4i5,/,
     1       '          dij,dkl           = ',2f15.10,/,
     1       '          dik,djl           = ',2f15.10,/,
     1       '          dil,djk           = ',2f15.10,/,
     1       '          val               = ', f15.10,/,
     2       '          valx1,valy1,valz1 = ',3f15.10,/,
     3       '          valx2,valy2,valz2 = ',3f15.10)
      end
      subroutine giaend(index,ii,jj,kk,ll,ixpk,xp,xk,disk)
      implicit REAL (a-h,o-z)
      parameter (mxatom=MAXAT)
      parameter (mxshel=MXSHEL)
      parameter (mxprim=MXPRIM)
      logical some
      logical out
      logical dbug
      logical pk,pandk,block
      logical pack2e
      logical disk
      logical nspdsk,nodisk
      integer thime,semisp,dirstp
      common/output/nprint
      common/restar1/nrest
      common/iofile/ir,iw,ip
      common/sqfile/ijk,ipk
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,zan(mxatom),
     1                                         c(3,mxatom)
      common/nshel/
     1ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),cf(mxprim),cg(mxprim),
     1kstart(mxshel),katom(mxshel),ktype(mxshel),kng(mxshel),
     2kloc(mxshel),kmin(mxshel),kmax(mxshel),nshell
      common/shlnos/qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     1 mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     2 nij,ij,kl,ijkl
      common/giao_shlint/cutoff,itotal,icount,irec,intloc,
     1                   ist,jst,kst,lst
      common/intprt/q(3),v(3),jc,n1(3),j1(3),j2(3),j3(3),j4(3),out
      common/intfil/nintmx
      common/inttyp/npople,nhondo
      common/intopt/nopk,nok,nosqur
      common/intdsk/maxrec,ijklnu,ijklsp,ijklno,nspdsk,nodisk
      common/pkfil /pk,pandk,block
      common/pckopt/nhex,ntupl,pack2e
      common/hndopt/hndtyp,dsktyp,filtyp,pcktyp
      common/dirpar/thize,thime,semisp,dirstp
      dimension xp(*),xk(*),ixpk(*)
c
      dbug=.false.
      some=.false.
      some=some.or.dbug       
c
      if(out) then
         call intout(j1(jc+1),j2(jc+1),j3(jc+1),j4(jc+1),
     1               q (jc+1),n1(jc+1),v (jc+1),1)
      endif
c
      if(index.eq.1) go to 20
c
c     ----- get ready for restart -----
c
      lstmax=kk
      if(.not.pk.and.kk.eq.ii) lstmax=jj
      kstmax=jj
      if(.not.pk             ) kstmax=ii
      jstmax=ii
      istmax=nshell
c
      nrest=1
      ist=ii
      jst=jj
      kst=kk
      lst=ll+1
      if(lst.le.lstmax) go to 10
      lst=1
      kst=kk+1
      if(kst.le.kstmax) go to 10
      kst=1
      jst=jj+1
      if(jst.le.jstmax) go to 10
      jst=1
      ist=ii+1
      if(ist.gt.istmax) go to 20
   10 continue
c
      ihondo=0
      if(npople.eq.0) ihondo=1
      write(iw,9996)
      write(iw,9999) irec,icount,ist,jst,kst,lst,ihondo
      nxx=icount-1
      call giaf10(    xp,ixpk,nxx,nintmx)
      itotal=itotal+nxx
      ntot  =itotal
      nint  =nintmx*(irec-1)+icount-1
      write(iw,9998) ntot,nint
      write(iw,9997) irec,icount,irec
c
      return
c
c     ----- all done -----
c
   20 continue
      nrest=0
      ist=1
      jst=1
      kst=1
      lst=1
c
c     ----- total completion -----
c
      nxx=icount-1
      nxx=-nxx
      call giaf10(    xp,ixpk,nxx,nintmx)
      itotal=itotal+(-nxx)
      intloc=icount
      ntot  =itotal
      nint  =nintmx*(irec-1)+icount-1
c
      if(dbug) write(iw,9998) ntot,nint
      if(dbug) write(iw,9997) irec,icount,irec
      if(some) write(iw,9995)
      if(dbug) write(iw,*)    'end of -giaend-'
c
      return
 9999 format(' irec, intloc, ist, jst, kst, lst, nhondo = ',/,
     1       i10,i6,5i4)
 9998 format(/,16x,' total number of 2e-integrals = ',i15,'  (',i15,')')
 9997 format(6x,i10,' record(s) of 2e-integrals on -ijk-',
     1       ' (intloc = ',i10,' )',/,
     2       ' (also',i10,' record(s) on -ipk- if -pandk-=.t.)',/)
 9996 format(' warning   .......   warning   .......',
     1       ' this job must be restarted ..... ')
 9995 format(' ...... end of -giao- two-electron integrals ......')
      end
      subroutine giaind(ijklg)
      implicit REAL (a-h,o-z)
      logical iieqjj,kkeqll,ijeqkl,ijgtkl,ijltkl
      logical keqi
      logical ieqj0,keql0,same0,first
      logical samtyp,samsym
      logical pk,pandk,block
      common/hfk/xx,u(9),w(9),nroots
      common/root/yy,t(12),v(12),mroots
      common/pkfil/pk,pandk,block
      common/intxch/norgin(3),iexch
      common/shlxch/igi(35,3),igj(35,3),igk(35,3),igl(35,3)
      common/shlbas/maxtyp,maxnum
      common/shlequ/iieqjj,kkeqll,ijeqkl,ijgtkl,ijltkl
      common/shllmn/igxyz(4,35),jgxyz(4,35),kgxyz(4,35),lgxyz(4,35)
      common/shlpar/lit,ljt,lkt,llt,loci,locj,lock,locl,
     1              mini,minj,mink,minl,maxi,maxj,maxk,maxl
      common/shlnum/numi,numj,numk,numl,ij,kl,ijkl
      common/shlold/numi0,numj0,numk0,numl0,
     1              ieqj0,keql0,same0,first
      common/shlgia/npwr,modtyp,modlit,modlkt
      dimension ijklx(35),ijkly(35),ijklz(35)
      dimension ijkln(5)
      dimension ijklg(4,1)
c
      data ijkln /   1,  4, 10, 20, 35/
      data ijklx /   0,  1,  0,  0,  2,  0,  0,  1,  1,  0,
     1               3,  0,  0,  2,  2,  1,  0,  1,  0,  1,
     2               4,  0,  0,  3,  3,  1,  0,  1,  0,  2,
     3               2,  0,  2,  1,  1/
      data ijkly /   0,  0,  1,  0,  0,  2,  0,  1,  0,  1,
     1               0,  3,  0,  1,  0,  2,  2,  0,  1,  1,
     2               0,  4,  0,  1,  0,  3,  3,  0,  1,  2,
     3               0,  2,  1,  2,  1/
      data ijklz /   0,  0,  0,  1,  0,  0,  2,  0,  1,  1,
     1               0,  0,  3,  0,  1,  0,  1,  2,  2,  1,
     2               0,  0,  4,  0,  1,  0,  1,  3,  3,  0,
     3               2,  2,  1,  1,  2/
c
      samtyp=numi.eq.numi0.and.numj.eq.numj0.and.
     1       numk.eq.numk0.and.numl.eq.numl0          .and.(.not.first)
      samsym=(iieqjj.eqv.ieqj0).and.
     1       (kkeqll.eqv.keql0).and.(ijeqkl.eqv.same0).and.(.not.first)
      numi0=numi
      numj0=numj
      numk0=numk
      numl0=numl
      ieqj0=iieqjj
      keql0=kkeqll
      same0=ijeqkl
      first=.false.
      if(samtyp.and.samsym.and..not.pk) return
c
c     ----- prepare indices for pairs of (i,j) functions -----
c
      ni=numl*numk*numj
      nj=numl*numk
      do 10 i=mini,maxi
      igxyz(1,i)  =ni*(i-mini)+1
   10 igi(i,iexch)=igxyz(1,i)
      llkjt=llt*modlkt*ljt
      do 20 i=1,ijkln(lit)
      igxyz(2,i)=ijklx(i)*llkjt+1
      igxyz(3,i)=ijkly(i)*llkjt+1
   20 igxyz(4,i)=ijklz(i)*llkjt+1
c
      do 30 j=minj,maxj
      jgxyz(1,j)  =nj*(j-minj)
   30 igj(j,iexch)=jgxyz(1,j)
      llkt=llt*modlkt
      do 40 j=1,ijkln(ljt)
      jgxyz(2,j)=ijklx(j)*llkt
      jgxyz(3,j)=ijkly(j)*llkt
   40 jgxyz(4,j)=ijklz(j)*llkt
c
c     ----- prepare indices for pairs of (k,l) functions -----
c
      nk=numl
      nl=1
      do 110 k=mink,maxk
      kgxyz(1,k)  =nk*(k-mink)
  110 igk(k,iexch)=kgxyz(1,k)
      do 120 l=minl,maxl
      lgxyz(1,l)  =nl*(l-minl)
  120 igl(l,iexch)=lgxyz(1,l)
c
      do 130 k=1,ijkln(lkt)
      kgxyz(2,k)=ijklx(k)*llt
      kgxyz(3,k)=ijkly(k)*llt
  130 kgxyz(4,k)=ijklz(k)*llt
      do 140 l=1,ijkln(llt)
      lgxyz(2,l)=ijklx(l)
      lgxyz(3,l)=ijkly(l)
  140 lgxyz(4,l)=ijklz(l)
c
c     ----- prepare indices for (ij/kl) -----
c
  200 continue
c
      ijkl=0
      do 240 i=mini,maxi
      jmax=maxj
      if(iieqjj) jmax=i
      do 230 j=minj,jmax
      kmax=maxk
      if(ijeqkl) kmax=i
      do 220 k=mink,kmax
      keqi=ijeqkl.and.k.eq.i
      lmax=maxl
      if(kkeqll) lmax=k
      if(keqi  ) lmax=j
      do 210 l=minl,lmax
      ijkl=ijkl+1
      ijklg(1,ijkl)=((igxyz(1,i)+jgxyz(1,j))+kgxyz(1,k))+lgxyz(1,l)
      ijklg(2,ijkl)=((igxyz(2,i)+jgxyz(2,j))+kgxyz(2,k))+lgxyz(2,l)
      ijklg(3,ijkl)=((igxyz(3,i)+jgxyz(3,j))+kgxyz(3,k))+lgxyz(3,l)
      ijklg(4,ijkl)=((igxyz(4,i)+jgxyz(4,j))+kgxyz(4,k))+lgxyz(4,l)
  210 continue
  220 continue
  230 continue
  240 continue
c
c     ----- set number of quadrature points -----
c
      nroots=(npwr+lit+ljt+lkt+llt-2)/2
      mroots=nroots
      return
      end
      subroutine giashl(ish,jsh,ksh,lsh)
      implicit REAL (a-h,o-z)
      parameter (mxatom=MAXAT)
      parameter (mxshel=MXSHEL)
      parameter (mxprim=MXPRIM)
      logical iieqjj,kkeqll,ijeqkl,ijgtkl,ijltkl
      logical spi,spj,spk,spl,spij,spkl,spijkl
      logical expndi,expndk
      common/infoa/nat,ich,mul,num,nnp,ne,na,nb,zan(mxatom),
     1                                          c(3,mxatom)
      common/nshel/
     1ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),cf(mxprim),cg(mxprim),
     1kstart(mxshel),katom(mxshel),ktype(mxshel),kng(mxshel),
     2kloc(mxshel),kmin(mxshel),kmax(mxshel),nshell
      common/shltyp/spi,spj,spk,spl,spij,spkl,spijkl
      common/shlxpn/expndi,expndk
      common/shlpar/lit,ljt,lkt,llt,loci,locj,lock,locl,
     1              mini,minj,mink,minl,maxi,maxj,maxk,maxl
      common/shlnum/numi,numj,numk,numl,ij,kl,ijkl
      common/shlequ/iieqjj,kkeqll,ijeqkl,ijgtkl,ijltkl
      common/shlgia/npwr,modtyp,modlit,modlkt
      common/atmgia/xi,yi,zi,xj,yj,zj,xk,yk,zk,xl,yl,zl,
     1              tijx,tijy,tijz,tklx,tkly,tklz,
     2              qijx,qijy,qijz,qklx,qkly,qklz 
c
      iieqjj=ish.eq.jsh
      kkeqll=ksh.eq.lsh
      ijeqkl=ish.eq.ksh.and.jsh.eq.lsh
      ijgtkl=max(ish,jsh).gt.max(ksh,lsh)
      ijltkl=max(ish,jsh).lt.max(ksh,lsh)
c
c     ----- ishell -----
c
      iat=katom(ish)
      xi =c(1,iat)
      yi =c(2,iat)
      zi =c(3,iat)
      lit=ktype(ish)
      mini=kmin(ish)
      maxi=kmax(ish)
      numi=maxi-mini+1
      loci=kloc(ish)-mini
      spi=lit.eq.2.and.mini.eq.1
c
c     ----- jshell -----
c
      jat=katom(jsh)
      xj =c(1,jat)
      yj =c(2,jat)
      zj =c(3,jat)
      ljt=ktype(jsh)
      minj=kmin(jsh)
      maxj=kmax(jsh)
      numj=maxj-minj+1
      locj=kloc(jsh)-minj
      spj=ljt.eq.2.and.minj.eq.1
      spij=spi.or.spj
      expndi=lit.ge.ljt
c
c     ----- -giao- factors -----
c
      tijx=xi-xj
      tijy=yi-yj
      tijz=zi-zj
      qijx=yi*zj-zi*yj
      qijy=zi*xj-xi*zj
      qijz=xi*yj-yi*xj
c
c     ----- kshell -----
c
      kat=katom(ksh)
      xk =c(1,kat)
      yk =c(2,kat)
      zk =c(3,kat)
      lkt=ktype(ksh)
      mink=kmin(ksh)
      maxk=kmax(ksh)
      numk=maxk-mink+1
      lock=kloc(ksh)-mink
      spk=lkt.eq.2.and.mink.eq.1
c
c     ----- lshell -----
c
      lat=katom(lsh)
      xl =c(1,lat)
      yl =c(2,lat)
      zl =c(3,lat)
      llt=ktype(lsh)
      minl=kmin(lsh)
      maxl=kmax(lsh)
      numl=maxl-minl+1
      locl=kloc(lsh)-minl
      spl=llt.eq.2.and.minl.eq.1
      spkl=spk.or.spl
      spijkl=spij.or.spkl
      expndk=lkt.ge.llt
c
c     ----- -giao- factors -----
c
      tklx=xk-xl
      tkly=yk-yl
      tklz=zk-zl
      qklx=yk*zl-zk*yl
      qkly=zk*xl-xk*zl
      qklz=xk*yl-yk*xl
c
c     ----- -giao- power -----
c
      modlit=lit+npwr
      modlkt=lkt+npwr
      return
      end
      subroutine giagnm(gnm,ng,nmax,mmax,
     1                  b00,b01,b10,c00,d00,f00)
      implicit REAL (a-h,o-z)
      logical nmaxs,nmaxp,mmaxs,mmaxp
      common/shlgnm/nmaxs,nmaxp,mmaxs,mmaxp
      dimension gnm(ng,nmax,mmax)
      dimension c00(ng),d00(ng),f00(ng)
      dimension b00(ng,1),b01(ng,1),b10(ng,1)
      data zero,one /0.0d+00,1.0d+00/
c
c     ----- g(0,0) -----
c
      do 10 ig=1,ng
      gnm(ig,1,1)=f00(ig)
   10 continue
      if(nmaxs.and.mmaxs) return
      if(nmaxs) go to 30
c
c     ----- g(1,0) = c00 * g(0,0) -----
c
      do 20 ig=1,ng
      gnm(ig,2,1)=c00(ig)*gnm(ig,1,1)
   20 continue
c
   30 continue
      if(mmaxs) go to 60
c
c     ----- g(0,1) = d00 * g(0,0) -----
c
      do 40 ig=1,ng
      gnm(ig,1,2)=d00(ig)*gnm(ig,1,1)
   40 continue
      if(nmaxs) go to 60
c
c     ----- g(1,1) = b00 * g(0,0) + d00 * g(1,0) -----
c
      do 50 ig=1,ng
      gnm(ig,2,2)=b00(ig,1)*gnm(ig,1,1)+d00(ig)*gnm(ig,2,1)
   50 continue
c
   60 continue
      iimax=max(nmax-1,mmax-1)
      do 70 m=2,iimax
      do 70 ig=1,ng
   70 b00(ig,m)=b00(ig,m-1)+b00(ig,1)
c
      if(nmaxp) go to 120
c
c     ----- g(n+1,0) = n * b10 * g(n-1,0) + c00 * g(n,0) -----
c
      do 80 n=2,nmax-1
      do 80 ig=1,ng
      b10(ig,n)=b10(ig,n-1)+b10(ig,1)
   80 continue
      do 90 n=2,nmax-1
      do 90 ig=1,ng
      gnm(ig,n+1,1)=b10(ig,n-1)*gnm(ig,n-1,1)+c00(ig)*gnm(ig,n,1)
   90 continue
      if(mmaxs) go to 110
c
c     ----- g(n,1) = n * b00 * g(n-1,0) + d00 * g(n,0) -----
c
      do 100 n=2,nmax-1
      do 100 ig=1,ng
      gnm(ig,n+1,2)=b00(ig,n)*gnm(ig,n,1)+d00(ig)*gnm(ig,n+1,1)
  100 continue
c
  110 continue
c
  120 continue
      if(mmaxp) go to 170
c
c     ----- g(0,m+1) = m * b01 * g(0,m-1) + d00 * g(o,m) -----
c
      do 130 m=2,mmax-1
      do 130 ig=1,ng
      b01(ig,m)=b01(ig,m-1)+b01(ig,1)
  130 continue
      do 140 m=2,mmax-1
      do 140 ig=1,ng
      gnm(ig,1,m+1)=b01(ig,m-1)*gnm(ig,1,m-1)+d00(ig)*gnm(ig,1,m)
  140 continue
      if(nmaxs) go to 160
c
c     ----- g(1,m) = m * b00 * g(0,m-1) + c00 * g(0,m) -----
c
      do 150 m=2,mmax-1
      do 150 ig=1,ng
      gnm(ig,2,m+1)=b00(ig,m)*gnm(ig,1,m)+c00(ig)*gnm(ig,1,m+1)
  150 continue
c
  160 continue
  170 if(nmaxp.or.mmaxp) return
c
c     ----- g(n+1,m) = n * b10 * g(n-1,m  ) -----
c                    +     c00 * g(n  ,m  )
c                    + m * b00 * g(n  ,m-1)
c
      do 180 m=2,mmax-1
      do 180 n=2,nmax-1
      do 180 ig=1,ng
      gnm(ig,n+1,m+1)=b10(ig,n-1)*gnm(ig,n-1,m+1)+
     1                c00(ig    )*gnm(ig,n  ,m+1)+
     2                b00(ig,m  )*gnm(ig,n  ,m  )
  180 continue
c
      return
      end
      subroutine giaxyz(gijkl,hijkl,gnkl,hnkl,fnkl,gnm,hnm,
     1                  ng,nmax,mmax,nimax,njmax,nkmax,nlmax,
     2                  dij,dkl,expndi,expndk)
      implicit REAL (a-h,o-z)
      logical expndi,expndk
      dimension gijkl(ng,3,nlmax*nkmax,njmax,nimax)
      dimension hijkl(ng,3,nlmax*nkmax*njmax,nimax)
      dimension  gnkl(ng,nlmax,nkmax,nmax)
      dimension  hnkl(ng,nlmax*nkmax,nmax)
      dimension  fnkl(ng,nlmax*nkmax*nmax)
      dimension   gnm(ng,nmax,mmax)
      dimension   dij(ng)
      dimension   dkl(ng)
c
c     ----- g(n,k,l) -----
c
      if(expndk) go to 40
c
      do 30 nk=1,nkmax
      do 10 nl=1,nlmax
      do 10  n=1,nmax
      do 10 ig=1,ng
   10 gnkl(ig,nl,nk,n)=gnm(ig,n,nl)
      if(nk.eq.nkmax) go to 30
      iimax=mmax-nk
      do 20  m=1,iimax
      do 20  n=1,nmax
      do 20 ig=1,ng
   20 gnm(ig,n,m)=dkl(ig)*gnm(ig,n,m)+gnm(ig,n,m+1)
   30 continue
c
      go to 100
   40 continue
c
      do 70 nl=1,nlmax
      do 50 nk=1,nkmax
      do 50  n=1,nmax
      do 50 ig=1,ng
   50 gnkl(ig,nl,nk,n)=gnm(ig,n,nk)
      if(nl.eq.nlmax) go to 70
      iimax=mmax-nl
      do 60  m=1,iimax
      do 60  n=1,nmax
      do 60 ig=1,ng
   60 gnm(ig,n,m)=dkl(ig)*gnm(ig,n,m)+gnm(ig,n,m+1)
   70 continue
c
  100 continue
c
c     ----- g(i,j,k,l) -----
c
      if(expndi) go to 140
c
      do 130 ni=1,nimax
      do 110 lkj=1,nlmax*nkmax*njmax
      do 110 ig=1,ng
  110 hijkl(ig,1,lkj,ni)=fnkl(ig,lkj)
      if(ni.eq.nimax) go to 130
      iimax=nmax-ni
      do 120  n=1,iimax
      do 120 nk=1,nkmax
      do 120 nl=1,nlmax
      do 120 ig=1,ng
  120 gnkl(ig,nl,nk,n)=dij(ig)*gnkl(ig,nl,nk,n)+gnkl(ig,nl,nk,n+1)
  130 continue
c
      return
  140 continue
c
      do 170 nj=1,njmax
      do 150 ni=1,nimax
      do 150 lk=1,nlmax*nkmax
      do 150 ig=1,ng
  150 gijkl(ig,1,lk,nj,ni)=hnkl(ig,lk,ni)
      if(nj.eq.njmax) go to 170
      iimax=nmax-nj
      do 160  n=1,iimax
      do 160 nk=1,nkmax
      do 160 nl=1,nlmax
      do 160 ig=1,ng
  160 gnkl(ig,nl,nk,n)=dij(ig)*gnkl(ig,nl,nk,n)+gnkl(ig,nl,nk,n+1)
  170 continue
c
      return
      end
      subroutine xyzgia(gijkl,ng,nimax,njmax,nkmax,nlmax)
      implicit REAL (a-h,o-z)
      dimension gijkl(ng,3,3,nlmax,nkmax,njmax,nimax)
c
c     ----- shift down for x1,y1,z1 -giao- integrals -----
c
      do ni=1,nimax-1
         do nj=1,njmax
            do nk=1,nkmax
               do nl=1,nlmax
                  do ixyz=1,3
                     do ig=1,ng
      gijkl(ig,ixyz,2,nl,nk,nj,ni)=gijkl(ig,ixyz,1,nl,nk,nj,ni+1)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
c
c     ----- shift down for x2,y2,z2 -giao- integrals -----
c
      do ni=1,nimax
         do nj=1,njmax
            do nk=1,nkmax-1
               do nl=1,nlmax
                  do ixyz=1,3
                     do ig=1,ng
      gijkl(ig,ixyz,3,nl,nk,nj,ni)=gijkl(ig,ixyz,1,nl,nk+1,nj,ni)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
c
      return
      end
      subroutine spdgia(ng,nr,gijkl,ijklg,xyz,
     1                  sj,sk,sl,dijsi,dijsj,dklsk,dklsl)
      implicit REAL (a-h,o-z)
      logical is,js,ks,ls
      logical ijs,ijks,ijkls
      logical spi,spj,spk,spl,spij,spkl,spijkl
      logical iieqjj,kkeqll,ijeqkl,ijgtkl,ijltkl
      logical keqi
      logical nojs,nols
      common/shltyp/spi,spj,spk,spl,spij,spkl,spijkl
      common/shlequ/iieqjj,kkeqll,ijeqkl,ijgtkl,ijltkl
      common/shllmn/igxyz(4,35),jgxyz(4,35),kgxyz(4,35),lgxyz(4,35)
      common/shlpar/lit,ljt,lkt,llt,loci,locj,lock,locl,
     1              mini,minj,mink,minl,maxi,maxj,maxk,maxl
      common/shlnum/numi,numj,numk,numl,ij,kl,ijkl
      dimension gijkl(7,*),ijklg(4,*)
      dimension  xyz(ng*nr,3,3,*)
      dimension sj(*),sk(*),sl(*)
      dimension dijsi(*),dijsj(*),dklsk(*),dklsl(*)
      data zero /0.0d+00/
c
      if(spijkl) go to 50
c
c     ----- no shared exponents ; sum up ( ix * iy * iz ) -----
c
      do 20 n=1,ijkl
         nn=ijklg(1,n)
         nx=ijklg(2,n)
         ny=ijklg(3,n)
         nz=ijklg(4,n)
         dum1=zero
         dum2=zero
         dum3=zero
         dum4=zero
         dum5=zero
         dum6=zero
         dum7=zero
         do igr=1,ng*nr
      dum1=dum1+xyz(igr,1,1,nx)*xyz(igr,2,1,ny)*xyz(igr,3,1,nz)
      dum2=dum2+xyz(igr,1,2,nx)*xyz(igr,2,1,ny)*xyz(igr,3,1,nz)
      dum3=dum3+xyz(igr,1,1,nx)*xyz(igr,2,2,ny)*xyz(igr,3,1,nz)
      dum4=dum4+xyz(igr,1,1,nx)*xyz(igr,2,1,ny)*xyz(igr,3,2,nz)
      dum5=dum5+xyz(igr,1,3,nx)*xyz(igr,2,1,ny)*xyz(igr,3,1,nz)
      dum6=dum6+xyz(igr,1,1,nx)*xyz(igr,2,3,ny)*xyz(igr,3,1,nz)
      dum7=dum7+xyz(igr,1,1,nx)*xyz(igr,2,1,ny)*xyz(igr,3,3,nz)
         enddo     
         gijkl(1,nn)=gijkl(1,nn)+dum1
         gijkl(2,nn)=gijkl(2,nn)+dum2
         gijkl(3,nn)=gijkl(3,nn)+dum3
         gijkl(4,nn)=gijkl(4,nn)+dum4
         gijkl(5,nn)=gijkl(5,nn)+dum5
         gijkl(6,nn)=gijkl(6,nn)+dum6
         gijkl(7,nn)=gijkl(7,nn)+dum7
   20 continue
      return
c
   50 continue
c
c     ----- shared exponents ; form ( ix * iy * iz ) -----
c
      n=0
      do 800 i=mini,maxi
      is=spi.and.i.eq.1
      nojs=iieqjj.and.i.eq.mini
c
      jmax=maxj
      if(iieqjj) jmax=i
      do 700 j=minj,jmax
      js=spj.and.j.eq.1
c
      if(js.and..not.nojs) then
         if(is) then
            do 110 igr=1,ng*nr
  110       sj(igr)=dijsj(igr)*dijsi(igr)
         else
            do 120 igr=1,ng*nr
  120       sj(igr)=dijsj(igr)
         endif
      else
         if(is) then
            do 130 igr=1,ng*nr
  130       sj(igr)=dijsi(igr)
         endif
      endif
      ijs=is.or.js
c
      kmax=maxk
      if(ijeqkl) kmax=i
      do 600 k=mink,kmax
      keqi=ijeqkl.and.k.eq.i
      ks=spk.and.k.eq.1
      nols=kkeqll.and.k.eq.mink
c
      if(ks) then
         if(ijs) then
            do 210 igr=1,ng*nr
  210       sk(igr)=dklsk(igr)*sj(igr)
         else
            do 220 igr=1,ng*nr
  220       sk(igr)=dklsk(igr)
         endif
      else
         if(ijs) then
            do 230 igr=1,ng*nr
  230       sk(igr)=sj(igr)
         endif
      endif
      ijks=ijs.or.ks
c
      lmax=maxl
      if(kkeqll) lmax=k
      if(keqi  ) lmax=j
      do 500 l=minl,lmax
      ls=spl.and.l.eq.1
c
      if(ls.and..not.nols) then
         if(ijks) then
            do 310 igr=1,ng*nr
  310       sl(igr)=dklsl(igr)*sk(igr)
         else
            do 320 igr=1,ng*nr
  320       sl(igr)=dklsl(igr)
         endif
      else
         if(ijks) then
            do 330 igr=1,ng*nr
  330       sl(igr)=sk(igr)
         endif
      endif
      ijkls=ijks.or.ls
c
      n=n+1
      nn=ijklg(1,n)
      nx=ijklg(2,n)
      ny=ijklg(3,n)
      nz=ijklg(4,n)
c
      if(ijkls) then
         dum1=zero
         dum2=zero
         dum3=zero
         dum4=zero
         dum5=zero
         dum6=zero
         dum7=zero
         do 410 igr=1,ng*nr
            dum1=dum1+xyz(igr,1,1,nx)*xyz(igr,2,1,ny)
     1                               *xyz(igr,3,1,nz)*sl(igr)
            dum2=dum2+xyz(igr,1,2,nx)*xyz(igr,2,1,ny)
     1                               *xyz(igr,3,1,nz)*sl(igr)
            dum3=dum3+xyz(igr,1,1,nx)*xyz(igr,2,2,ny)
     1                               *xyz(igr,3,1,nz)*sl(igr)
            dum4=dum4+xyz(igr,1,1,nx)*xyz(igr,2,1,ny)
     1                               *xyz(igr,3,2,nz)*sl(igr)
            dum5=dum5+xyz(igr,1,3,nx)*xyz(igr,2,1,ny)
     1                               *xyz(igr,3,1,nz)*sl(igr)
            dum6=dum6+xyz(igr,1,1,nx)*xyz(igr,2,3,ny)
     1                               *xyz(igr,3,1,nz)*sl(igr)
            dum7=dum7+xyz(igr,1,1,nx)*xyz(igr,2,1,ny)
     1                               *xyz(igr,3,3,nz)*sl(igr)
  410    continue
      else
         dum1=zero
         dum2=zero
         dum3=zero
         dum4=zero
         dum5=zero
         dum6=zero
         dum7=zero
         do 420 igr=1,ng*nr
            dum1=dum1+xyz(igr,1,1,nx)*xyz(igr,2,1,ny)
     1                               *xyz(igr,3,1,nz)
            dum2=dum2+xyz(igr,1,2,nx)*xyz(igr,2,1,ny)
     1                               *xyz(igr,3,1,nz)
            dum3=dum3+xyz(igr,1,1,nx)*xyz(igr,2,2,ny)
     1                               *xyz(igr,3,1,nz)
            dum4=dum4+xyz(igr,1,1,nx)*xyz(igr,2,1,ny)
     1                               *xyz(igr,3,2,nz)
            dum5=dum5+xyz(igr,1,3,nx)*xyz(igr,2,1,ny)
     1                               *xyz(igr,3,1,nz)
            dum6=dum6+xyz(igr,1,1,nx)*xyz(igr,2,3,ny)
     1                               *xyz(igr,3,1,nz)
            dum7=dum7+xyz(igr,1,1,nx)*xyz(igr,2,1,ny)
     1                               *xyz(igr,3,3,nz)
  420    continue
      endif
      gijkl(1,nn)=gijkl(1,nn)+dum1
      gijkl(2,nn)=gijkl(2,nn)+dum2
      gijkl(3,nn)=gijkl(3,nn)+dum3
      gijkl(4,nn)=gijkl(4,nn)+dum4
      gijkl(5,nn)=gijkl(5,nn)+dum5
      gijkl(6,nn)=gijkl(6,nn)+dum6
      gijkl(7,nn)=gijkl(7,nn)+dum7
c
  500 continue
  600 continue
  700 continue
  800 continue
c
      return
      end
      subroutine giazer(gijkl,ijklg)
      implicit REAL (a-h,o-z)
      common/shlnum/numi,numj,numk,numl,ij,kl,ijkl
      dimension gijkl(7,*),ijklg(4,*)
      data zero /0.0d+00/
c
c     ----- zero out -gijkl- -----
c
      do n=1,ijkl
         do i=1,7
            gijkl(i,ijklg(1,n))=zero
         enddo
      enddo
c
      return
      end
      subroutine giamak(gijkl,ijklg)
      implicit REAL (a-h,o-z)
      common/shlnum/numi,numj,numk,numl,ij,kl,ijkl
      common/atmgia/xi,yi,zi,xj,yj,zj,xk,yk,zk,xl,yl,zl,
     1              tijx,tijy,tijz,tklx,tkly,tklz,
     2              qijx,qijy,qijz,qklx,qkly,qklz 
      dimension gijkl(7,*),ijklg(4,*)
      data zero /0.0d+00/
c
c     ----- build -giao- x,y,z terms -----
c
      do n=1,ijkl
         dum1=gijkl(1,ijklg(1,n))
         dum2=gijkl(2,ijklg(1,n))
         dum3=gijkl(3,ijklg(1,n))
         dum4=gijkl(4,ijklg(1,n))
         dum5=gijkl(5,ijklg(1,n))
         dum6=gijkl(6,ijklg(1,n))
         dum7=gijkl(7,ijklg(1,n))
         tmp1=        dum1
         tmp2= ( qijx*dum1 + tijy*dum4 - tijz*dum3 )
         tmp3= ( qijy*dum1 + tijz*dum2 - tijx*dum4 )
         tmp4= ( qijz*dum1 + tijx*dum3 - tijy*dum2 )
         tmp5= ( qklx*dum1 + tkly*dum7 - tklz*dum6 )
         tmp6= ( qkly*dum1 + tklz*dum5 - tklx*dum7 )
         tmp7= ( qklz*dum1 + tklx*dum6 - tkly*dum5 )
         gijkl(1,ijklg(1,n))=tmp1
         gijkl(2,ijklg(1,n))=tmp2
         gijkl(3,ijklg(1,n))=tmp3
         gijkl(4,ijklg(1,n))=tmp4
         gijkl(5,ijklg(1,n))=tmp5
         gijkl(6,ijklg(1,n))=tmp6
         gijkl(7,ijklg(1,n))=tmp7
      enddo
c
      return
      end
