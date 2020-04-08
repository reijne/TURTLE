c     deck=util1
c ******************************************************
c ******************************************************
c             =   util1  =
c ******************************************************
c ******************************************************
      subroutine input
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c.... this routine reads a data card and scans it for non - space fields
c    the number of fields is stored in jump, the starting point of a
c    field in istrt(i) and the number of characters in that field
c    in inumb(i).
      character * 132 char1,char2,char1c
      common/workc/char1,char2,char1c
c
c     mxtoken - the maximum number of tokens on an input line
c               for space separated tokens: 132/2 plus 1 to allow
c               for look ahead (132 is the line-length defined in
c               common/workc)
c
      integer mxtoken
      parameter(mxtoken=67)
      integer jrec, jump, istrt, inumb, iwidth, nend, nstart
      integer nline, noline, jwidth, nerr
      logical oswit, oterm, oflush
      common/work/jrec,jump,istrt(mxtoken),inumb(mxtoken),iwidth,
     + nend(mxtoken),nstart(mxtoken), nline,noline,jwidth,nerr,
     + oswit,oterm,oflush
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
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
      dimension xchar(4),ncol(mxtoken),xcol(2),xprin(132)
      character*132 char2c
      data xblnk,xstp,xcol,xtab/' ','.',':','\\','\t'/
      data nchar,xchar/4,'#','?','>','<'/
      data xcomma,xequal/',','='/
c
      nline=nline+1
      if(nline.le.noline)go to 150
      if(oswit) call caserr2('unexpected end of data file')
 100  read(ird,'(a132)',end=300)char2
c
c     Replace any tab-characters with spaces to ensure proper detection
c     of white-space
c
      do i=1,len(char2)
         if (char2(i:i).eq.xtab) char2(i:i)=xblnk
      enddo
c
c     If the input is to be echoed to the output then print the input
c     line here.
c
      last = lstchr(char2)
      if (oecho.and.opg_root()) then
       do i = 1, last
        xprin(i) = char2(i:i)
       enddo
       write(iwr,400) (xprin(i),i=1,last)
 400   format(1x,'>>>>> ', 132a1)
      endif
      if (last.eq.0) goto 100
c
c     If this line is a comment then decide whether to print it or not
c     and subsequently skip to the next input line (i.e. go to 100).
c
      do i=1,nchar
         if(char2(1:1).eq.xchar(i))go to 110
      enddo
      go to 80
 110  char1(1:131)=char2(2:132)
      if(opg_root().and..not.oecho)write(iwr,90)char1
 90   format(/
     *' comment :-',1x,a79)
      go to 100
 80   k=jwidth
c
c     char2c preserves the case of input strings, this is required
c     for things such as file names which are case sensitive
c
      char2c = char2
      call lcase(char2(1:132))
c
      mark=0

      do 130 i=1,jwidth
        if(char2(i:i).ne.xcol(1).and.
     +     char2(i:i).ne.xcol(2)) go to 130
        mark=mark+1
        if (mark.gt.mxtoken) then
          call caserr("mxtoken exceeded, too many tokens on input line")
        endif
        ncol(mark)=i
 130  continue
      noline=1
      if(mark.ne.0)go to 140
      if (noline.gt.mxtoken) then
        call caserr("mxtoken exceeded, too many tokens on input line")
      endif
      nstart(noline)=1
      nend(noline)=jwidth
      go to 200
 140  i=ncol(mark)+1
      if(i.le.jwidth) then
       do j=i,jwidth
        if(char2(j:j).ne.xblnk) go to 170
       enddo
      endif
      k=ncol(mark)-1
      mark=mark-1
c
 170  noline=mark+1
      if (mark.ge.mxtoken) then
        call caserr("mxtoken exceeded, too many tokens on input line")
      endif
      nstart(1)=1
      do i=1,mark
        j=ncol(i)
        nend(i)=j-1
        nstart(i+1)=j+1
      enddo
      if (noline.gt.mxtoken) then
        call caserr("mxtoken exceeded, too many tokens on input line")
      endif
      nend(noline)=k
 200  nline=1
c
 150  jump=0
      jrec=0
      isw=0
      if (nline.gt.mxtoken) then
        call caserr("mxtoken exceeded, too many tokens on input line")
      endif
      nbegin=nstart(nline)
      nfini=nend  (nline)
      iwidth=nfini-nbegin+1
      char1(1:iwidth)=char2(nbegin:nfini)
      char1c(1:iwidth)=char2c(nbegin:nfini)
c
c...   allow comments preceeded by ! like fortran (surrounded by blanks)
c
      do i=iwidth-1,2,-1
         if (char1(i:i).eq.'!'.and.char2(i+1:i+1).eq.' '
     1                        .and.char2(i-1:i-1).eq.' ') then
            iwidth = i-1
            go to 101
         end if
      end do
101   continue
c
c     pad the line from the last character onwards with spaces
c
      j=iwidth+1
      if (j.lt.132) then
        do i=j,132
          char1c(i:i)=xblnk
          char1(i:i)=xblnk
        enddo
      endif
c
c     find all the space separated tokens
c
      do 40 i = 1,iwidth
        if(char1(i:i).eq.xblnk)go to 30
        if(char1(i:i).eq.xcomma)go to 30
        if(char1(i:i).eq.xequal)go to 30

        mark=ichar(char1(i:i))
        if(mark.eq.13)goto 30

        if (isw.le.0) then
          jump = jump +1
          if (jump.gt.mxtoken) then
            call caserr(
     +      "mxtoken exceeded, too many tokens on input line")
          endif
          istrt(jump) = i
          inumb(jump) = 0
          isw=1
        endif
        if (jump.gt.mxtoken) then
          call caserr("mxtoken exceeded, too many tokens on input line")
        endif
        inumb(jump) = inumb(jump) + 1
        go to 40
30      isw = 0
40    continue
      return
c
 300  oswit=.true.
      jump=0
      jrec=0
      return
      end
      subroutine errout(n)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character * 132 char1,char2,char1c
      common/workc/char1,char2,char1c
c
c     mxtoken - the maximum number of tokens on an input line
c               for space separated tokens: 132/2 plus 1 to allow
c               for look ahead (132 is the line-length defined in
c               common/workc)
c
      integer mxtoken
      parameter(mxtoken=67)
      integer jrec, jump, istrt, inumb, iwidth, nend, nstart
      integer nline, noline, jwidth, nerr
      logical oswit, oterm, oflush
      common/work/jrec,jump,istrt(mxtoken),inumb(mxtoken),iwidth,
     + nend(mxtoken),nstart(mxtoken), nline,noline,jwidth,nerr,
     + oswit,oterm,oflush
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      data xpt,xstp/    '*','.'/
      jrec=-1
      write(iwr,50)char1
50    format(1x,a132)
      do 60 i=1,iwidth
60    char1(i:i)=xstp
      char1(n:n)=xpt
      write(iwr,50)char1
      return
      end
      subroutine outrec
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *132 char
      character * 132 char1,char2,char1c
      common/workc/char1,char2,char1c
c
c     mxtoken - the maximum number of tokens on an input line
c               for space separated tokens: 132/2 plus 1 to allow
c               for look ahead (132 is the line-length defined in
c               common/workc)
c
      integer mxtoken
      parameter(mxtoken=67)
      integer jrec, jump, istrt, inumb, iwidth, nend, nstart
      integer nline, noline, jwidth, nerr
      logical oswit, oterm, oflush
      common/work/jrec,jump,istrt(mxtoken),inumb(mxtoken),iwidth,
     + nend(mxtoken),nstart(mxtoken), nline,noline,jwidth,nerr,
     + oswit,oterm,oflush
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
      char=char1(1:iwidth)
      write(iwr,50)char
50    format(1x,a132)
      return
      end
      subroutine inpa(zguf)
c.... this routine examines the contents of char1  and extracts a
c    character string of 8 chars. this string is stored in iguf .
c    characters beyond the eighth in any field are ignored
c      dimension ibuf(8)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character * 132 char1,char2,char1c
      common/workc/char1,char2,char1c
c
c     mxtoken - the maximum number of tokens on an input line
c               for space separated tokens: 132/2 plus 1 to allow
c               for look ahead (132 is the line-length defined in
c               common/workc)
c
      integer mxtoken
      parameter(mxtoken=67)
      integer jrec, jump, istrt, inumb, iwidth, nend, nstart
      integer nline, noline, jwidth, nerr
      logical oswit, oterm, oflush
      common/work/jrec,jump,istrt(mxtoken),inumb(mxtoken),iwidth,
     + nend(mxtoken),nstart(mxtoken), nline,noline,jwidth,nerr,
     + oswit,oterm,oflush
c
      data xblnk/' '/
c
      jrec = jrec + 1
      zguf = xblnk
      if(jrec .gt. jump) return
      n = inumb(jrec)
      if(n.gt.8)n=8
      zguf = char1(istrt(jrec):istrt(jrec)+n-1)
c
      return
      end
      subroutine lcase(string)
      character*(*) string
c
c...  convert to lower case (use ascii table)
c
      ll = len(string)
      do 125 i=1,ll
         mark = ichar(string(i:i))
         if (mark.ge.65.and.mark.le.90) string(i:i) = char(mark+32)
125   continue
c
      return
      end
      subroutine inpan(zguf)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character*(*) zguf
      character * 132 char1,char2,char1c
      common/workc/char1,char2,char1c
c
c     mxtoken - the maximum number of tokens on an input line
c               for space separated tokens: 132/2 plus 1 to allow
c               for look ahead (132 is the line-length defined in
c               common/workc)
c
      integer mxtoken
      parameter(mxtoken=67)
      integer jrec, jump, istrt, inumb, iwidth, nend, nstart
      integer nline, noline, jwidth, nerr
      logical oswit, oterm, oflush
      common/work/jrec,jump,istrt(mxtoken),inumb(mxtoken),iwidth,
     + nend(mxtoken),nstart(mxtoken), nline,noline,jwidth,nerr,
     + oswit,oterm,oflush
c
      data xblnk/' '/
      jrec = jrec + 1
      nguf=len(zguf)
      zguf = xblnk
      if(jrec .gt. jump) return
      n = inumb(jrec)
      if(n.gt.nguf)n=nguf
      zguf = char1(istrt(jrec):istrt(jrec)+n-1)
c
      return
      end
c
c inpan variant preserving input case
c
      subroutine inpanpc(zguf)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character*(*) zguf
      character * 132 char1,char2,char1c
      common/workc/char1,char2,char1c
c
c     mxtoken - the maximum number of tokens on an input line
c               for space separated tokens: 132/2 plus 1 to allow
c               for look ahead (132 is the line-length defined in
c               common/workc)
c
      integer mxtoken
      parameter(mxtoken=67)
      integer jrec, jump, istrt, inumb, iwidth, nend, nstart
      integer nline, noline, jwidth, nerr
      logical oswit, oterm, oflush
      common/work/jrec,jump,istrt(mxtoken),inumb(mxtoken),iwidth,
     + nend(mxtoken),nstart(mxtoken), nline,noline,jwidth,nerr,
     + oswit,oterm,oflush
c
      data xblnk/' '/
      jrec = jrec + 1
      nguf=len(zguf)
      zguf = xblnk
      if(jrec .gt. jump) return
      n = inumb(jrec)
      if(n.gt.nguf)n=nguf
      zguf = char1c(istrt(jrec):istrt(jrec)+n-1)
c
      return
      end
      subroutine inpa4(yguf)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character * 132 char1,char2,char1c
      common/workc/char1,char2,char1c
c
c     mxtoken - the maximum number of tokens on an input line
c               for space separated tokens: 132/2 plus 1 to allow
c               for look ahead (132 is the line-length defined in
c               common/workc)
c
      integer mxtoken
      parameter(mxtoken=67)
      integer jrec, jump, istrt, inumb, iwidth, nend, nstart
      integer nline, noline, jwidth, nerr
      logical oswit, oterm, oflush
      common/work/jrec,jump,istrt(mxtoken),inumb(mxtoken),iwidth,
     + nend(mxtoken),nstart(mxtoken), nline,noline,jwidth,nerr,
     + oswit,oterm,oflush
c
      data xblnk/' '/
      jrec = jrec + 1
      yguf = xblnk

      if(jrec .gt. jump) return
      n = inumb(jrec)
      if(n.gt.4)n=4
      yguf = char1(istrt(jrec):istrt(jrec)+n-1)
c
      return
      end
      subroutine inpal(zguf)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character*(*) zguf
      character * 132 char1,char2,char1c
      common/workc/char1,char2,char1c
c
c     mxtoken - the maximum number of tokens on an input line
c               for space separated tokens: 132/2 plus 1 to allow
c               for look ahead (132 is the line-length defined in
c               common/workc)
c
      integer mxtoken
      parameter(mxtoken=67)
      integer jrec, jump, istrt, inumb, iwidth, nend, nstart
      integer nline, noline, jwidth, nerr
      logical oswit, oterm, oflush
      common/work/jrec,jump,istrt(mxtoken),inumb(mxtoken),iwidth,
     + nend(mxtoken),nstart(mxtoken), nline,noline,jwidth,nerr,
     + oswit,oterm,oflush
c
      data xblnk/' '/
      jrec = jrec + 1
      nguf=len(zguf)
      zguf = xblnk
      if(jrec .gt. jump) return
c     n = inumb(jrec)
      n = inumb(jump)
      if(n.gt.nguf)n=nguf
      zguf = char1(istrt(jrec):istrt(jump)+n-1)
      jrec = jump 
c
      return
      end
      subroutine incon
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character*8 idir,ivar
c Parameter settings for the Molecular Surface program mscon
c
c maxatm     maximum number of atoms
c maxtyp     maximum number of atom types
c maxnbr     maximum number of neighbors an atom may have
c maxsph     maximum number of surface points on a sphere
c maxcir     maximum number of surface points on a circle
c maxarc     maximum number of surface points on an arc
c maxppp     maximum number of surface points per probe
c maxyon     maximum number of yon probes
c maxvic     maximum number of victim probes
c maxeat     maximum number of eaters of a probe's surface
c maxcub     maximum number of cubes in one direction
c mxconpts   maximum number of surface points
c
c maxatm must be greater than or equal to maxyon
c because they share the same cubing arrays
c
      parameter (maxatm=2000)
      parameter (maxrfa2=12000)
      parameter (maxtyp=100)
      parameter (maxnbr=200)
      parameter (maxsph=1000)
      parameter (maxcir=1000)
      parameter (maxarc=1000)
      parameter (maxppp=1000)
      parameter (maxyon=1000)
      parameter (maxvic=6000)
      parameter (maxeat=1000)
      parameter (maxcub=40)
      parameter (mxconpts=4000)

cxxx  real*8 co,rtype,molnum,d,rp
      real*8 co,rtype,d,rp
      integer*2 molnum
      integer nmaxcon,nmincon,natom
      logical isurdens,iminpoint
      common/coninp/co(3,maxrfa2+1),rtype(maxrfa2+1),d,rp
      common/coninp2/molnum(maxrfa2+1)
      common/coninp3/natom,nmaxcon,nmincon,isurdens,iminpoint
c
c     mxtoken - the maximum number of tokens on an input line
c               for space separated tokens: 132/2 plus 1 to allow
c               for look ahead (132 is the line-length defined in
c               common/workc)
c
      integer mxtoken
      parameter(mxtoken=67)
      integer jrec, jump, istrt, inumb, iwidth, nend, nstart
      integer nline, noline, jwidth, nerr
      logical oswit, oterm, oflush
      common/work/jrec,jump,istrt(mxtoken),inumb(mxtoken),iwidth,
     + nend(mxtoken),nstart(mxtoken), nline,noline,jwidth,nerr,
     + oswit,oterm,oflush
c
c
      real*8 dstgrp
      integer igroup, iunit
      common /drf_in/ dstgrp,igroup,iunit
c
      integer indstg, indstmn, indstmx, insp, inrp, inpro, inprj
      integer inbndr, inbndl, insw, insd, insphr, insr
      common /drf_un/ indstg,indstmn,indstmx,
     + insp,inrp,inpro,inprj,inbndr,inbndl,insw,insd,
     + insphr,insr
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c     common/drf_un/indstg,indstmn,indstmx,
c    1 inbndr,inbndl,ind,inrp
      dimension idir(4)
c
      data idir/'spdens','mxsurpts','mnsurpts','rprobe'/
c defaults
      d=1.0d0
      rp=1.0d0
cxxx  nmaxcon=0
cxxx  isurdens=.true.
      nmaxcon=200
      isurdens=.false.
      iminpoint=.false.
      isetden = 0
      isetmax = 0
10    if (jrec.ge.jump) goto 100
      call inpa(ivar)
      i=locatc(idir,4,ivar)
      if (i.eq.0) call caserr2(
     +    'unknown parameter in Connolly directive')
      goto(20,30,40,50) , i
20    call inpf(d)
      isurdens=.true.
      isetden = 1
      insp = 1
      go to 10
30    call inpi(nmaxcon)
      isetmax = 1
      go to 10
40    call inpi(nmincon)
      iminpoint=.true.
      go to 10
50    call inpf(rp)
      inrp = 1
      inpro = 1
      go to 10
100   continue
      if ((isetden .eq. 1) .and. (isetmax .eq. 1))
     + call caserr2
     + ('either specify surface point density OR maximum no. of points') 
      if (nmincon .ge. nmaxcon) then
        write(iwr,1001) nmaxcon
 1001   format(/,' Input error: you may not set mnsurpts >=',
     1  ' mxsurpts= ',i4)
        call caserr2('Connolly input error detected')
      endif
      return
      end 
      subroutine indiel(rfsurf)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character*8 idir,ivar,surftyp,outdiel
      character*4 dielopt1, dielopt2
      character*10 rfsurf
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
c Parameter settings for the Molecular Surface program mscon
c
c maxatm     maximum number of atoms
c maxtyp     maximum number of atom types
c maxnbr     maximum number of neighbors an atom may have
c maxsph     maximum number of surface points on a sphere
c maxcir     maximum number of surface points on a circle
c maxarc     maximum number of surface points on an arc
c maxppp     maximum number of surface points per probe
c maxyon     maximum number of yon probes
c maxvic     maximum number of victim probes
c maxeat     maximum number of eaters of a probe's surface
c maxcub     maximum number of cubes in one direction
c mxconpts   maximum number of surface points
c
c maxatm must be greater than or equal to maxyon
c because they share the same cubing arrays
c
      parameter (maxatm=2000)
      parameter (maxrfa2=12000)
      parameter (maxtyp=100)
      parameter (maxnbr=200)
      parameter (maxsph=1000)
      parameter (maxcir=1000)
      parameter (maxarc=1000)
      parameter (maxppp=1000)
      parameter (maxyon=1000)
      parameter (maxvic=6000)
      parameter (maxeat=1000)
      parameter (maxcub=40)
      parameter (mxconpts=4000)

c
      integer ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep
      integer neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga
      integer idipcal,iclinte,iclintd,ieffpol,isodis,iclintr
      integer iextdip,itsolv,itermax,isolsav,imomsav,idisadd
      integer ianal,ineqex,ithole,irevdis,ihbond,igrppol
      integer maxci_drf
      common/drfpar1/ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep,
     +               neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga,
     +               idipcal,iclinte,iclintd,ieffpol,isodis,iclintr,
     +               iextdip,itsolv,itermax,isolsav,imomsav,idisadd,
     +               ianal,ineqex,ithole,irevdis,ihbond,igrppol,
     +               maxci_drf
c
      integer nxtpts, npol, nexp, natint, nshint, namb, nspec, ngran
      common/drfpar2/nxtpts,npol,nexp,natint,nshint,namb,nspec,ngran
c
      integer ndima, ndimb, nbem, ndim, nexp4, npol3, nwtr, nwtc, nzfa
      integer lomga, nomga, nchd, ngrnam, ngrnam2, nzfp, nzfn
      integer nneq, nneqrf, maxneq, neqdim, nqdim, nqcls
      common/drfpar3/ndima,ndimb,nbem,ndim,nexp4,npol3,nwtr,nwtc,nzfa,
     +               lomga,nomga,nchd,ngrnam,ngrnam2,nzfp,nzfn,
     +               nneq,nneqrf,maxneq,neqdim,nqdim,nqcls
c
      integer ifldin, ifldout, idrfout, modxza, irepopt
      integer iexpza, icmexp, iadexp, nodpe, igetden, iarfcal, iqmclr
      common/drfpar4/ifldin,ifldout,idrfout,modxza,irepopt,
     +               iexpza,icmexp,iadexp,nodpe,igetden,iarfcal,
     +               iqmclr
c
      real*8 gamdrf, dstmin, dstmax, hbondl, hbondr, afact, rfact
      real*8 cvgrel, agrpe, agrpm, agrpc, scffact, acur
      real*8 conci_drf
      common/drfpar5/gamdrf,dstmin,dstmax,hbondl,hbondr,afact,rfact,
     +               cvgrel,agrpe,agrpm,agrpc,scffact,acur,conci_drf
c
      integer iind, isur1, isur2, ilwt, ilvr, illur, ilindx
      common/drfpar6/iind,isur1,isur2,ilwt,ilvr,illur,ilindx
c
      integer itwoeps,ioponly,itwosur,ixnsurf,ibemout,iuniout
      common /drfbem1/ itwoeps,ioponly,itwosur,ixnsurf,ibemout,iuniout
c
      integer leveli, levelo, nbem1, nbem2, nraw
      common /drfbem2/ leveli,levelo,nbem1,nbem2,nraw
c
      character*32 solnam
      real*8 radmax, eps1, eps2, spherad, rprobe, rprobej
      real*8 solrad, swidth, sdist
      real*8 kappa
      real*8 kappa1, kappa2, kappas
      common /drfbem3/ radmax,eps1,eps2,kappa1,kappa2,spherad,rprobe,
     +                 rprobej,solnam,solrad,swidth,sdist
c
      real*8 radat, cgrav
      common /drfbem4/ radat(maxat),cgrav(3)
c
cxxx  real*8 co,rtype,molnum,d,rp
      real*8 co,rtype,d,rp
      integer*2 molnum
      integer nmaxcon,nmincon,natom
      logical isurdens,iminpoint
      common/coninp/co(3,maxrfa2+1),rtype(maxrfa2+1),d,rp
      common/coninp2/molnum(maxrfa2+1)
      common/coninp3/natom,nmaxcon,nmincon,isurdens,iminpoint
c
c     mxtoken - the maximum number of tokens on an input line
c               for space separated tokens: 132/2 plus 1 to allow
c               for look ahead (132 is the line-length defined in
c               common/workc)
c
      integer mxtoken
      parameter(mxtoken=67)
      integer jrec, jump, istrt, inumb, iwidth, nend, nstart
      integer nline, noline, jwidth, nerr
      logical oswit, oterm, oflush
      common/work/jrec,jump,istrt(mxtoken),inumb(mxtoken),iwidth,
     + nend(mxtoken),nstart(mxtoken), nline,noline,jwidth,nerr,
     + oswit,oterm,oflush
c
c
      real*8 dstgrp
      integer igroup, iunit
      common /drf_in/ dstgrp,igroup,iunit
c
      integer indstg, indstmn, indstmx, insp, inrp, inpro, inprj
      integer inbndr, inbndl, insw, insd, insphr, insr
      common /drf_un/ indstg,indstmn,indstmx,
     + insp,inrp,inpro,inprj,inbndr,inbndl,insw,insd,
     + insphr,insr
c
      dimension idir(14)
c
      data idir/'surface','solvent','dieltyp',
     1          'epsstat','epsopt','kappas','kappao',
     2          'sradius','bemlev','dielout','connolly',
     3          'juffer','solvrad','end'/
c defaults
      ibem=5
      rfsurf='connolly'
      solnam='x'
      eps1=1.0d0
      eps2=1.0d0
      kappa1=0.0d0
      kappa2=0.0d0
      itwoeps = 0
      ioponly = 0
      spherad = 1000.0d0
      levelo = 0
      ibemout = 0
      d=1.0d0
      rp=1.0d0
      solrad=-1.0d0
cxxx  nmaxcon=0
cxxx  isurdens=.true.
      nmaxcon=200
      isurdens=.false.
      isetden = 0
      isetmax = 0
      swidth = -1.0d0
      sdist = -1.0d0
      rprobej = 1.0d0
10    continue
      call input
      if (idrfout .ge. 1) call outrec
cxxxx if (jrec.ge.jump) goto 200
      call inpa(ivar)
      i=locatc(idir,14,ivar)
      if (i.eq.0) 
     1  call caserr2('unknown parameter in Dielectric directive')
      goto(20,30,40,50,60,70,80,90,100,110,120,130,140,150) , i
20    call inpa(surftyp)
      if (surftyp .eq. 'sphere') then
        ibem = 4
      elseif (surftyp .eq. 'juffer') then
        ibem = 3
        rfsurf = 'juffer'
      elseif (surftyp .eq. 'connolly') then
        ibem = 5
        rfsurf = 'connolly'
      else
        call caserr2('unknown option in Surface directive')
      endif
      go to 10
30    call inpan(solnam)
      go to 10
40    call inpa4(dielopt1)
      call inpa4(dielopt2)
      istat = 1
      ioptic = 0
      if ((dielopt1 .ne. 'stat') .and. (dielopt2 .ne. 'stat')) istat=0
      if ((dielopt1 .eq. 'opt') .or. (dielopt2 .eq. 'opt')) ioptic=1
      if ((istat .eq. 0) .and. (ioptic .eq. 1)) ioponly = 1
      if ((istat .eq. 1) .and. (ioptic .eq. 1)) itwoeps = 1
      go to 10
50    call inpf(eps1)
      if (eps1 .lt. 1.0d0) eps1=1.0d0
      go to 10
60    call inpf(eps2)
      if (eps2 .lt. 1.0d0) eps2=1.0d0
      go to 10
70    call inpf(kappa1)
      if (kappa1 .lt. 0.0d0) kappa1=0.0d0
      go to 10
80    call inpf(kappa2)
      if (kappa2 .lt. 0.0d0) kappa2=0.0d0
      go to 10
90    call inpf(spherad)
      insphr = 1
      go to 10
100    call inpi(levelo)
      go to 10
110    call inpa(outdiel)
       if (outdiel .eq. 'standard') then
         ibemout = 0
       elseif (outdiel .eq. 'some') then
         ibemout = 2
       elseif (outdiel .eq. 'moderate') then
         ibemout = 5
       elseif (outdiel .eq. 'extended') then
         ibemout = 10
       elseif (outdiel .eq. 'all') then
         ibemout = 100
       else
         call caserr2('unknown parameter in Dielout directive')
       endif
      go to 10
120    call incon
      go to 10
130    continue
 135   if (jrec.ge.jump) goto 136
         call inpa(outdiel)
       if (outdiel(1:7) .eq. 'cylwdth') then
         call inpf(swidth)
         insw = 1
       else if (outdiel(1:8) .eq. 'surfdist') then
         call inpf(sdist)
         insd = 1
       else if (outdiel(1:6) .eq. 'rprobe') then
         call inpf(rprobej)
         inprj = 1
       else
         call caserr2('unknown parameter in Juffer directive')
       endif
        go to 135
 136   continue
      go to 10
140   call inpf(solrad)
      insr = 1
      go to 10
150   continue
200   continue
      rprobe = rp
      if ((solnam .eq. 'x') .and. (solrad .lt. 0.0d0)) then
        solrad = 0.0d0
      endif
      if ((ioponly .eq. 1).and.(eps2 .ne. 1.0d0)) eps1 = eps2
      if ((ioponly .eq. 1).and.(kappa2 .ne. 0.0d0)) kappa1 = kappa2
      return
      end

      subroutine inmontec
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character*8 idir,ivar
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
cafc
c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     this file contains the parameter settings for piet van duijnen's
c     reactionfield programme together with the sizes needed for the
c     v. duijnen - gamess interface as developed by 
c     nederkoorn and v. lenthe
c
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c     ==================================================================
c
c     parameter setting for reaction field programme van duijnen:
c
c=======================================================================
c
c     
      parameter (mxat = 10000)
      parameter (mxpol = 10000)
      parameter (mxpts = 10000)
      parameter (mxspe = 50)
      parameter (mxamb = 10)
      parameter (mxgrp = 500)
      parameter (mxnum = 1024)
      parameter (mxsh = 512)
      parameter (mxgs = 2048)
      parameter (mxnpts = 20)
      parameter (mxneq = 255)
      parameter (mcwrkmx = 10)
c
c ==================================================================
c
c     parameter setting new drf routines (drf/dimpar):
c
c===================================================================
c
      parameter (mxst = 5)
      parameter (mxgran = 10)
cmw      parameter (mxnpts = 50)
      parameter (mxgrpar = (mxgran*(mxgran+1))/2)
c
c     ==================================================================
c
c     parameter setting for v. duijnen - gamess interface:
c
c=======================================================================
c
      parameter (ncatom=500)
c
c  ncatom is maximal allowed classical atoms in zmat
c  total allowed atoms in zmat = maxnz + ncatom
c
c  maxrfa is maximal number of classical atoms in zmat + coordinates
c  total maximal atoms allowed is maxrfa + maxat
c
      parameter (maxrfa=mxpts)
c
c  maxamb is maximal number of ambigous atoms. it is set equal
c  to maximal z-cards possible in gamess input , i.e. 50
c
      parameter (maxamb=mxamb)
c
c=============================================
c for the pes-routines
c=============================================
c
      parameter (mxgr1=128, mxpes=100)
c
c=============================================
c for the drfexpc routines
c=============================================
c
      parameter (maxex=1000)
c
c=============================================
c for the drfsurf routine
c=============================================
c
      parameter (maxinv=10,maxnrc=10000)
c      parameter (maxnrc=10000)
c
c=============================================
c for the analrf routine
c=============================================
c
c
      real*8 ureft, uold, eqm, eqmol, eclas, eclasol, eclase, ecleol
      real*8 eclasd, ecldol, eclasr, eclrol, eint, eintol, eelst, eelstol
      real*8 edisp, edisol, erep, erepol, extnol, epol, epolol
      real*8 eneq, eneqol, eneqp, eneqpol, edifol
      common /mcener/ ureft(mxst),uold(mxst),eqm(mxst),eqmol(mxst),
     +  eclas,eclasol,eclase,ecleol,eclasd,ecldol,eclasr,eclrol,
     +  eint(mxst),eintol(mxst),eelst(mxst),eelstol(mxst),edisp(mxst),
     +  edisol(mxst),erep(mxst),erepol(mxst),extnol(mxst),epol(mxst),
     +  epolol(mxst),eneq(mxst),eneqol(mxst),eneqp(mxst),eneqpol(mxst),
     +  edifol(mxst*(mxst+1)/2)
c
      real*8 utav, utav2, utavq, utavq2, utavc
      real*8 utavc2, utavce, utavce2, utavcd, utavcd2, utavcr, utavcr2
      real*8 utavi, utavi2, utave, utave2, utavd, utavd2, utavr, utavr2
      real*8 utavp, utavp2, utavn, utavn2, utavnp, utavnp2
      real*8 utavdf, utavdf2, q1tot
      common /mcener2/ utav(mxst),utav2(mxst),utavq(mxst),utavq2(mxst),
     +  utavc,utavc2,utavce,utavce2,utavcd,utavcd2,utavcr,utavcr2,
     +  utavi(mxst),utavi2(mxst),utave(mxst),utave2(mxst),utavd(mxst),
     +  utavd2(mxst),utavr(mxst),utavr2(mxst),utavp(mxst),utavp2(mxst),
     +  utavn(mxst),utavn2(mxst),utavnp(mxst),utavnp2(mxst),
     +  utavdf(mxst*(mxst+1)/2),utavdf2(mxst*(mxst+1)/2),q1tot(mxst)
c
      integer itotacc
      common /mcres/ itotacc
c
      real*8 eclasg, eclsolg, eclaseg, ecleolg, eclasdg, ecldolg
      real*8 eclasrg, eclrolg, eintg, eintolg, eelstg, elstolg
      real*8 edispg, edisolg, erepg, erepolg, extnolg
      common/mcener3/eclasg(mxgrpar),eclsolg(mxgrpar),eclaseg(mxgrpar),
     +  ecleolg(mxgrpar),eclasdg(mxgrpar),ecldolg(mxgrpar),
     +  eclasrg(mxgrpar),eclrolg(mxgrpar),eintg(mxgran,mxst),
     +  eintolg(mxgran,mxst),eelstg(mxgran,mxst),elstolg(mxgran,mxst),
     +  edispg(mxgran,mxst),edisolg(mxgran,mxst),erepg(mxgran,mxst),
     +  erepolg(mxgran,mxst),extnolg(mxgran,mxst)
c
      real*8 utavcg, utavcg2, utavceg, utvceg2, utavcdg, utvcdg2
      real*8 utavcrg, utvcrg2, utavig, utavig2, utaveg, utaveg2
      real*8 utavdg, utavdg2, utavrg, utavrg2
      common/mcener4/utavcg(mxgrpar),utavcg2(mxgrpar),utavceg(mxgrpar),
     +  utvceg2(mxgrpar),utavcdg(mxgrpar),utvcdg2(mxgrpar),
     +  utavcrg(mxgrpar),utvcrg2(mxgrpar),utavig(mxgran,mxst),
     +  utavig2(mxgran,mxst),utaveg(mxgran,mxst),utaveg2(mxgran,mxst),
     +  utavdg(mxgran,mxst),utavdg2(mxgran,mxst),utavrg(mxgran,mxst),
     +  utavrg2(mxgran,mxst)
c
c
      real*8 gmdrf
      integer iguess
      common /mcinp/ gmdrf(mxst),iguess(mxst)
c
      character * 132 char1,char2,char1c
      common/workc/char1,char2,char1c
c
c     mxtoken - the maximum number of tokens on an input line
c               for space separated tokens: 132/2 plus 1 to allow
c               for look ahead (132 is the line-length defined in
c               common/workc)
c
      integer mxtoken
      parameter(mxtoken=67)
      integer jrec, jump, istrt, inumb, iwidth, nend, nstart
      integer nline, noline, jwidth, nerr
      logical oswit, oterm, oflush
      common/work/jrec,jump,istrt(mxtoken),inumb(mxtoken),iwidth,
     + nend(mxtoken),nstart(mxtoken), nline,noline,jwidth,nerr,
     + oswit,oterm,oflush
c
c
      integer nacal
      common /rota1/ nacal
c
      real*8 temprt, symfac
      common /rota2/ temprt(10), symfac
c
c
      integer nsamp, nblock, nmoves, ncheck
      common /mcrunp1/ nsamp,nblock,nmoves,ncheck
c
      real*8 darot, dtrans, temp, onekt, excld, delmax, enmin
      common /mcrunp2/ darot,dtrans,temp,onekt,excld,delmax,enmin
c
      logical mcupdt,secstat,excess
      integer notrans, norot, imcout, iseed, imcst, iactst
      common /mcrunp3/ notrans,norot,imcout,iseed,imcst,iactst,
     +               mcupdt,secstat,excess
c
      character*16 outfor
      common /mcrunp4/ outfor
c
      real*8 ratmin, ratmax, amxrot, amxtrn, amnrot, amntrn
      common /mcrunp5/ ratmin,ratmax,amxrot,amxtrn,amnrot,amntrn
c
      real*8 gammc
      common /mcrunp6/ gammc(5)
c
c
      integer ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep
      integer neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga
      integer idipcal,iclinte,iclintd,ieffpol,isodis,iclintr
      integer iextdip,itsolv,itermax,isolsav,imomsav,idisadd
      integer ianal,ineqex,ithole,irevdis,ihbond,igrppol
      integer maxci_drf
      common/drfpar1/ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep,
     +               neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga,
     +               idipcal,iclinte,iclintd,ieffpol,isodis,iclintr,
     +               iextdip,itsolv,itermax,isolsav,imomsav,idisadd,
     +               ianal,ineqex,ithole,irevdis,ihbond,igrppol,
     +               maxci_drf
c
      integer nxtpts, npol, nexp, natint, nshint, namb, nspec, ngran
      common/drfpar2/nxtpts,npol,nexp,natint,nshint,namb,nspec,ngran
c
      integer ndima, ndimb, nbem, ndim, nexp4, npol3, nwtr, nwtc, nzfa
      integer lomga, nomga, nchd, ngrnam, ngrnam2, nzfp, nzfn
      integer nneq, nneqrf, maxneq, neqdim, nqdim, nqcls
      common/drfpar3/ndima,ndimb,nbem,ndim,nexp4,npol3,nwtr,nwtc,nzfa,
     +               lomga,nomga,nchd,ngrnam,ngrnam2,nzfp,nzfn,
     +               nneq,nneqrf,maxneq,neqdim,nqdim,nqcls
c
      integer ifldin, ifldout, idrfout, modxza, irepopt
      integer iexpza, icmexp, iadexp, nodpe, igetden, iarfcal, iqmclr
      common/drfpar4/ifldin,ifldout,idrfout,modxza,irepopt,
     +               iexpza,icmexp,iadexp,nodpe,igetden,iarfcal,
     +               iqmclr
c
      real*8 gamdrf, dstmin, dstmax, hbondl, hbondr, afact, rfact
      real*8 cvgrel, agrpe, agrpm, agrpc, scffact, acur
      real*8 conci_drf
      common/drfpar5/gamdrf,dstmin,dstmax,hbondl,hbondr,afact,rfact,
     +               cvgrel,agrpe,agrpm,agrpc,scffact,acur,conci_drf
c
      integer iind, isur1, isur2, ilwt, ilvr, illur, ilindx
      common/drfpar6/iind,isur1,isur2,ilwt,ilvr,illur,ilindx
c
      integer nopes, lpes, lpesc, intpesc, numpes, npesgrp, ipesout
      common /pescl/ nopes,lpes,lpesc,intpesc,numpes,npesgrp,ipesout
c
c
      integer itwoeps,ioponly,itwosur,ixnsurf,ibemout,iuniout
      common /drfbem1/ itwoeps,ioponly,itwosur,ixnsurf,ibemout,iuniout
c
      integer leveli, levelo, nbem1, nbem2, nraw
      common /drfbem2/ leveli,levelo,nbem1,nbem2,nraw
c
      character*32 solnam
      real*8 radmax, eps1, eps2, spherad, rprobe, rprobej
      real*8 solrad, swidth, sdist
      real*8 kappa
      real*8 kappa1, kappa2, kappas
      common /drfbem3/ radmax,eps1,eps2,kappa1,kappa2,spherad,rprobe,
     +                 rprobej,solnam,solrad,swidth,sdist
c
      real*8 radat, cgrav
      common /drfbem4/ radat(maxat),cgrav(3)
c

c
c-----  namelist $montec input parameters                 default
c
c       imcout : printflag for mc+rf+bem method               0
c       iseed  : seed for random generator              1234567
c       imcst  : number of states in monte carlo calculation  1
c       (maximum 5; may be altered by setting mxst in drf/dimpar)
c       iguess(mxst) : guess option for different states      5*0
c          options:
c             0: use standard guess (current density etc..
c                present on da10) for state ist
c             1: read guess from da10, collected in previous
c                calculations according to rank number of
c                the state (this is the user's responsibility)
c       ***note***
c       this means that state number ist in the monte carlo was
c       given nset=ist in $guess in a separate calculation on
c       this state!!!!
c
c       gmdrf(mxst): dispersion scaling parameter gamdrf for
c                    different states                         5*gamdrf
c                    cf. gamdrf in rfin
c
c       nsamp  : number of samples                            1
c       nblock : number of blocks in each sample              1
c       nmoves : number of attempted moves in each block      1
c       ncheck : number of moves to check acceptance after    0
c                too low or too high an acceptance will lead
c                to the adaptation of darot and dtrans
c                ncheck = 0 means do not perform this check
c                ncheck < 100 or > 10000 will lead to a warning
c       ratmax : maximum acceptance ratio if ncheck > 0      .5
c       ratmin : minimum acceptance ratio if ncheck > 0      .5
c       darot  : maximum rotation step at beginning of run
c                (in radians)                                .15 (13.5 d
c       dtrans : maximum translation step at beginning of run
c                (in au)                                     .45
c       amxrot : maximum rotation step during run
c                (in radians) (if ncheck > 0)                 1. (90 deg
c       amxtrn : maximum translation step during run
c                (in au) (if ncheck > 0 and notrans = 0)      2.
c       amnrot : minimum rotation step during run
c                (in radians) (if ncheck > 0)                 .01 (0.9 d
c       amntrn : minimum translation step during run
c                (in au) (if ncheck > 0 and notrans = 0)      .001
c       temp   : temperature (k)                              298.15
c       excld  : exclusion distance (only when ibem .ne. 0)   7.3
c                (default is the value for water)
c       delmax : maximum energy drop in monte carlo step     -0.01 hartr
c       enmin  : minumum energy threshold                    -1.0d+10 h
c       (if the computed energy is lower, the move will be rejected)
c
c       the last two flags are ad-hoc measures to prevent spurious
c       configurations
c
      dimension idir(21)
      character*8 errmsg(3)
c
c     data gmdrf/5*-1.0d0/
      data errmsg /'program','stop in','inmontec'/
      data nsmax,nblmax,nmovmax /10,10,10000/
      data drotmax /1.0/
c     data iguess /5*0/
      data idir/'imcout','iseed','imcst','nsamp',
     1          'nblock','nmoves','drot','dtrans',
     2          'ncheck','ratmin','ratmax','amnrot',
     3          'amxrot','amntrn','amxtrn','excld',
     4          'outfor','gmdrf','delmax','enmin',
     5          'end'/
c defaults
      do i=1,mxst
         gmdrf(i) = -1.0d0
         iguess(i) = 0
      end do
      imcout=0
      iseed=1234567
      imcst=1
      nsamp=1
      nblock=1
      nmoves=1
      darot=0.15d0
      dtrans=0.45d0
      ncheck=0
      ratmin=0.4d0
      ratmax=0.5d0
      amnrot=0.01d0
      amxrot=1.0d0
      amntrn=0.001d0
      amxtrn=2.0d0
      excld=7.3d0
      outfor='shorsom none'
      delmax = -1.0d-02
      enmin = -1.0d10
c
      pi = 4.0d0*atan(1.0d0)
c
caleko
c  temp is set equal to the temperature specified in 
c  the DRF-input
c
      temp=temprt(1)
 10   continue
      call input
      if (idrfout.ge.1) call outrec
      call inpa(ivar)
      i=locatc(idir,21,ivar)
      if (i.eq.0)
     1  call caserr2('unknown parameter in Monte Carlo directive')
      goto(20,30,40,50,60,70,80,90,100,110,120,130,140,
     2     150,160,170,180,190,200,210,220),i
20    call inpi(imcout)
      goto 10
30    call inpi(iseed)
      goto 10
40    call inpi(imcst)
      do 45 j=1,imcst
        iguess(j)=0
45    continue
caleko
c
c Some of the options available in HONDO have been
c hard-set to their default:
c  iguess(i), i=1,mcst is 0, that is monte carlo
c guess is done according to the current density 
c present on record 16 of the da10.
c Not possible to use density of a previous calculation.
c
      goto 10 
50    call inpi(nsamp)
      goto 10
60    call inpi(nblock)
      goto 10
70    call inpi(nmoves)
      goto 10
80    call inpf(darot)
      goto 10
90    call inpf(dtrans)
      goto 10
100   call inpi(ncheck)
      goto 10
110   call inpf(ratmin)
      goto 10
120   call inpf(ratmax)
      goto 10
130   call inpf(amnrot)
      goto 10
140   call inpf(amxrot)
      goto 10
150   call inpf(amntrn)
      goto 10
160   call inpf(amxtrn)
      goto 10
170   call inpf(excld)
      goto 10
180   call inpa4(outfor(1:4))
      call inpa4(outfor(5:8))
      call inpa4(outfor(9:12))
      goto 10
190   do 11 j=1,imcst
        call inpf(gmdrf(j))
11    continue
      goto 10
200   call inpf(delmax)
      goto 10
210   call inpf(enmin)
      goto 10
220   continue
      
      if (nopes .eq. 0) then
c 1-----
c  -----  overrule some input for classical pes scan calculation
c
        nsamp = 1
        nblock = 1
        nmoves = numpes
        ncheck = 0
c 1-----
      endif
c
c-----  check whether input values are allowed
c
      if ((nsamp .le. 0) .or. (nsamp .gt. nsmax)) then
        write (iwr,1001) nsmax, nsamp
 1001   format(/' check input: ', i6,' > nsamp >0 ',
     1          ' nsamp given as: ', i6)
        call hnderr(3,errmsg)
      endif
c
      if ((nblock .le. 0) .or. (nblock .gt. nblmax)) then
        write (iwr,1002) nblmax, nblock
 1002   format(/' check input: ', i6,' > nblock >0 ',
     1          ' nblock given as: ', i6)
        call hnderr(3,errmsg)
      endif
c
      if ((nmoves .le. 0) .or. (nmoves .gt. nmovmax)) then
        write (iwr,1003) nmovmax, nmoves
 1003   format(/' check input: ', i6,' > nmoves >0 ',
     1          ' nmoves given as: ', i6)
        call hnderr(3,errmsg)
      endif
c
      if ((darot .lt. 0.0) .or. (darot .gt. drotmax)) then
        write (iwr,1004) drotmax, darot
 1004   format(/' check input: ', f8.4,' > darot >0.0 ',
     1          ' darot given as: ', f8.4)
        call hnderr(3,errmsg)
      endif
c
      if (temp .lt. 0.0) then
        write (iwr,1005) temp
 1005   format(/' check input:  temp >0.0 ',
     1          ' temp given as: ', f8.4)
        call hnderr(3,errmsg)
      endif
c
      if (imcst .gt. mxst) then
        write(iwr,1006) mxst
 1006   format(/,' check input: imcst exceeds mxst =',i2,/,
     1  ' mxst may be set in drf/dimpar')
        call hnderr(3,errmsg)
      endif
cc
      darot = pi*darot
c
      notrans = 1
      if (dtrans .ne. 0.0d0) notrans = 0
      norot = 1
      if (darot .ne. 0.0d0) norot = 0
c
      if (notrans .eq. 1) then
        if (norot .eq. 1) then
          write (iwr, 1011)
 1011     format (/,' this mc run is of little use, since both ',
     2    'rotation and translation of the groups is disallowed')
          call hnderr(3,errmsg)
        endif
        write (iwr, 1012)
 1012   format (/,' translation of the classical groups is disallowed')
      endif
c
      if (norot .eq. 1) then
        write (iwr, 1013)
 1013   format (/,' rotation of the classical groups is disallowed')
      endif
c
      if (ncheck .eq. 0) then
c 1-----
        write (iwr,1021)
 1021   format (/, ' the parameters darot and dtrans will not be ',
     2  'changed during the run')
c 1-----
      else
c 1-----
        write (iwr,1022) ncheck
 1022   format (/, ' the parameters darot and dtrans will be ',
     2  'reviewed and adapted during the run',/,
     3  ' this will be done every ',i6,' moves')
        if (ncheck .lt. 100) write (iwr,1023)
 1023   format
     2 (/, ' ncheck quite low for meaningful check on parameters')
        if (ncheck .gt. 100000) write (iwr,1024)
 1024   format
     2 (/, ' ncheck quite high for meaningful check on parameters')
c
        if (ratmin .gt. ratmax) then
          tempor = ratmin
          ratmin = ratmax
          ratmax = tempor
          write (iwr,1031)
 1031     format (/ ' parameters ratmin and ratmax interchanged !!')
        endif
c 1-----
      endif
c
      if ((ncheck .gt. 0) .and. (dtrans .gt. amxtrn)) then
        write (iwr,1041) amxtrn, dtrans
 1041   format(/' check input: amxtrn=', f8.4,' > dtrans',
     1          ' dtrans given as: ', f8.4)
        call hnderr(3,errmsg)
      endif
c
      if ((ncheck .gt. 0) .and. (darot .gt. amxrot)) then
        write (iwr,1051) amxrot, darot
 1051   format(/' check input: amxrot=', f8.4,' < darot ',
     1          ' darot given as: ', f8.4)
        call hnderr(3,errmsg)
      endif
c
      if ((ncheck .gt. 0) .and. (dtrans .lt. amntrn)) then
c 1-----
        if (notrans .eq. 0) then
          write (iwr,1061) amntrn, dtrans
 1061     format(/' check input: amntrn=', f8.4,' < dtrans',
     1          ' dtrans given as: ', f8.4)
          call hnderr(3,errmsg)
        else
c   2-----
c           set amntrn zero, to avoid faults in translation step
c           parameter updating in mcstepp
c
          amntrn = 0.0d0
        endif
c 1-----
      endif
c
      if ((ncheck .gt. 0) .and. (darot .lt. amnrot)) then
c 1-----
        if (norot .eq. 0) then
          write (iwr,1071) amnrot, darot
 1071     format(/' check input: amnrot=', f8.4,' < darot',
     1          ' darot given as: ', f8.4)
          call hnderr(3,errmsg)
        else
c   2-----
c         set amnrot zero, to avoid faults in rotation step
c         parameter updating in mcstepp
c
          amnrot = 0.0d0
        endif
c 1-----
      endif
c
c-----  set exclusion distance if boundary is present
c       if < 0.0 in input, look it up from solvent information
c
      if (notrans .eq. 0) then
        if (iclintr .ne. 1) then
          iclintr = 1
          write(iwr,1081)
 1081     format(/,' note: classical energy expression forced',
     1    ' to include repulsion because translation is allowed')
        endif
        if ((ibem .ne. 0 ).and.(excld .lt. 0.0d0)) call exclset(excld)
      endif
      continue
      return
      end
      subroutine addup(a,b,c,n)
c------
c      returns array c = a + b of length n
c------
      implicit real*8  (a-h,o-z),integer  (i-n)
c
      dimension a(n), b(n), c(n)
c
      do 100, i = 1, n
c 1-----
        c(i) = a(i) + b(i)
c 1-----
  100 continue
c
      return
      end

      subroutine transform(a,trvect,trmatr)
caleko  
c      This routine translates the coordinates in a
c      using the vector trvect and then rotates them
c      using matrix trmat. These are usually the 
c      standard symmetry-imposed operations of gamess.
c
      implicit real*8  (a-h,o-z),integer  (i-n)
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
cafc
c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     this file contains the parameter settings for piet van duijnen's
c     reactionfield programme together with the sizes needed for the
c     v. duijnen - gamess interface as developed by 
c     nederkoorn and v. lenthe
c
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c     ==================================================================
c
c     parameter setting for reaction field programme van duijnen:
c
c=======================================================================
c
c     
      parameter (mxat = 10000)
      parameter (mxpol = 10000)
      parameter (mxpts = 10000)
      parameter (mxspe = 50)
      parameter (mxamb = 10)
      parameter (mxgrp = 500)
      parameter (mxnum = 1024)
      parameter (mxsh = 512)
      parameter (mxgs = 2048)
      parameter (mxnpts = 20)
      parameter (mxneq = 255)
      parameter (mcwrkmx = 10)
c
c ==================================================================
c
c     parameter setting new drf routines (drf/dimpar):
c
c===================================================================
c
      parameter (mxst = 5)
      parameter (mxgran = 10)
cmw      parameter (mxnpts = 50)
      parameter (mxgrpar = (mxgran*(mxgran+1))/2)
c
c     ==================================================================
c
c     parameter setting for v. duijnen - gamess interface:
c
c=======================================================================
c
      parameter (ncatom=500)
c
c  ncatom is maximal allowed classical atoms in zmat
c  total allowed atoms in zmat = maxnz + ncatom
c
c  maxrfa is maximal number of classical atoms in zmat + coordinates
c  total maximal atoms allowed is maxrfa + maxat
c
      parameter (maxrfa=mxpts)
c
c  maxamb is maximal number of ambigous atoms. it is set equal
c  to maximal z-cards possible in gamess input , i.e. 50
c
      parameter (maxamb=mxamb)
c
c=============================================
c for the pes-routines
c=============================================
c
      parameter (mxgr1=128, mxpes=100)
c
c=============================================
c for the drfexpc routines
c=============================================
c
      parameter (maxex=1000)
c
c=============================================
c for the drfsurf routine
c=============================================
c
      parameter (maxinv=10,maxnrc=10000)
c      parameter (maxnrc=10000)
c
c=============================================
c for the analrf routine
c=============================================
c
c
      real*8 xpts
      common /extpts/ xpts(3,mxpts)
c
      character*16 nxcent
      common /namext/ nxcent(mxpts)
c
      real*8 chrg
      common /extcha/ chrg(mxpts)
c
      real*8 polar
      common /extpol/ polar(mxpol)
c
      real*8 atpol
      common /extpol2/ atpol(mxpts)
c
      real*8 atpolt
      common /extpolt/ atpolt(6,mxpts)
c
      real*8 alfext
      common /extalf/ alfext(mxpts)
c
      integer mpol
      common /ipol/ mpol(mxpol)
c
      integer ncutpt
      common /drfcut/ ncutpt(mxpts)
c
      integer nspecl
      common /drfspe/ nspecl(mxspe)
c
      integer nvalel
      common /drfval/ nvalel(mxpts)
c
      real*8 vale
      common /clasval/ vale(mxpts)
c
      real*8 valat
      common /atval/ valat(mxpts)
c
      real*8 polars
      common /claspol/ polars(mxpts)
c
      real*8 alfat
      common /atalf/ alfat(mxat)
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      dimension a(3), b(3) 
      dimension trvect(3), trmatr(3,3)
      do 100, i=1,3
        b(i)= a(i) + trvect(i)
 100  continue

      call tform2(trmatr,b,a)
      return
      end  
      subroutine errors(n)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character*23 msg
      character * 132 char1,char2,char1c
      common/workc/char1,char2,char1c
c
c     mxtoken - the maximum number of tokens on an input line
c               for space separated tokens: 132/2 plus 1 to allow
c               for look ahead (132 is the line-length defined in
c               common/workc)
c
      integer mxtoken
      parameter(mxtoken=67)
      integer jrec, jump, istrt, inumb, iwidth, nend, nstart
      integer nline, noline, jwidth, nerr
      logical oswit, oterm, oflush
      common/work/jrec,jump,istrt(mxtoken),inumb(mxtoken),iwidth,
     + nend(mxtoken),nstart(mxtoken), nline,noline,jwidth,nerr,
     + oswit,oterm,oflush
c
      data msg /'atmol error number     '/
      if (n.eq.62) call caserr2('atmol i/o error (code 62)')
      write (msg(20:23),'(i4)') n
      nerr=n
      call caserr2(msg)
      return
      end

      subroutine caserr(log)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character*(*) log
c
      integer ERR_NO_CODE
      parameter(ERR_NO_CODE=0)
c
c
c ==== if the error will occur on all nodes, we can ====
c      handle the output more cleanly
c
      integer ERR_SYNC, ERR_ASYNC
      parameter(ERR_ASYNC=0)
      parameter(ERR_SYNC=301)
c
c ==== legal error classes - these control an explanatory  ====
c      message (one per class, see machscf.m).
c
      integer ERR_NO_CLASS, ERR_INCOMPREHENSIBLE,
     &  ERR_INCOMPATIBLE, ERR_DIMENSION, ERR_FIXED_DIMENSION,
     &  ERR_UNIMPLEMENTED, ERR_INTERNAL,  ERR_USER_DIMENSION,
     &  ERR_SWITCHED_OUT, ERR_UNLUCKY

      parameter(ERR_NO_CLASS=0)
c - fix input errors and try again
      parameter(ERR_INCOMPREHENSIBLE=101)
c - incompatible 
      parameter(ERR_INCOMPATIBLE=102)
c - code dimensions exceeded - redimension code
      parameter(ERR_DIMENSION=103)
c - code dimensions exceeded - redimension or contact support
      parameter(ERR_FIXED_DIMENSION=104)
c - request unimplemented option
      parameter(ERR_UNIMPLEMENTED=105)
c - internal error - please report to support
      parameter(ERR_INTERNAL=106)

c - input dimension exceeded - change allocation
c   in input file and re-run
      parameter(ERR_USER_DIMENSION=107)

c - option disabled at configure stage
      parameter(ERR_SWITCHED_OUT=108)

c - some algorithmic or case-specififc problem
      parameter(ERR_UNLUCKY=109)

c
c  System message flag
c
      integer ERR_NO_SYS, ERR_SYS
      parameter(ERR_NO_SYS=0)
      parameter(ERR_SYS=201)
c
      call gamerr(log,ERR_NO_CODE, ERR_NO_CLASS, ERR_ASYNC, ERR_NO_SYS)

      end

      subroutine caserr2(log)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character*(*) log
c
      integer ERR_NO_CODE
      parameter(ERR_NO_CODE=0)
c
c
c ==== if the error will occur on all nodes, we can ====
c      handle the output more cleanly
c
      integer ERR_SYNC, ERR_ASYNC
      parameter(ERR_ASYNC=0)
      parameter(ERR_SYNC=301)
c
c ==== legal error classes - these control an explanatory  ====
c      message (one per class, see machscf.m).
c
      integer ERR_NO_CLASS, ERR_INCOMPREHENSIBLE,
     &  ERR_INCOMPATIBLE, ERR_DIMENSION, ERR_FIXED_DIMENSION,
     &  ERR_UNIMPLEMENTED, ERR_INTERNAL,  ERR_USER_DIMENSION,
     &  ERR_SWITCHED_OUT, ERR_UNLUCKY

      parameter(ERR_NO_CLASS=0)
c - fix input errors and try again
      parameter(ERR_INCOMPREHENSIBLE=101)
c - incompatible 
      parameter(ERR_INCOMPATIBLE=102)
c - code dimensions exceeded - redimension code
      parameter(ERR_DIMENSION=103)
c - code dimensions exceeded - redimension or contact support
      parameter(ERR_FIXED_DIMENSION=104)
c - request unimplemented option
      parameter(ERR_UNIMPLEMENTED=105)
c - internal error - please report to support
      parameter(ERR_INTERNAL=106)

c - input dimension exceeded - change allocation
c   in input file and re-run
      parameter(ERR_USER_DIMENSION=107)

c - option disabled at configure stage
      parameter(ERR_SWITCHED_OUT=108)

c - some algorithmic or case-specififc problem
      parameter(ERR_UNLUCKY=109)

c
c  System message flag
c
      integer ERR_NO_SYS, ERR_SYS
      parameter(ERR_NO_SYS=0)
      parameter(ERR_SYS=201)
c
c  This version must only be called when it is known that
c  the error will occur on all nodes
c  calling it on a single node will cause the job to hang
c
      call gamerr(log,ERR_NO_CODE, ERR_NO_CLASS, ERR_SYNC, ERR_NO_SYS)

      end


      subroutine inpf (buf)
      implicit none
      character * 132 char1,char2,char1c
      common/workc/char1,char2,char1c
c
c     mxtoken - the maximum number of tokens on an input line
c               for space separated tokens: 132/2 plus 1 to allow
c               for look ahead (132 is the line-length defined in
c               common/workc)
c
      integer mxtoken
      parameter(mxtoken=67)
      integer jrec, jump, istrt, inumb, iwidth, nend, nstart
      integer nline, noline, jwidth, nerr
      logical oswit, oterm, oflush
      common/work/jrec,jump,istrt(mxtoken),inumb(mxtoken),iwidth,
     + nend(mxtoken),nstart(mxtoken), nline,noline,jwidth,nerr,
     + oswit,oterm,oflush
c
      intrinsic dble
      integer i,j,i1,i2,ie,ie1,ie2,iexp,isign,ibuff
      logical orep
      character*1 xchar(15)
      data xchar /'0','1','2','3','4','5','6','7','8','9'
     1,'+','-','.','e','d'/
      real*8 buf,ten
      data ten/10.0d0/
      buf=0.0d0
      jrec=jrec+1
      if (jrec.gt.jump) return
      i1=istrt(jrec)
      i2=i1+inumb(jrec)-1
      ie2=i2
c...  sign
      isign=1
      if (char1(i1:i1).eq.xchar(12))isign=-1
      if (char1(i1:i1).eq.xchar(12).or.
     +    char1(i1:i1).eq.xchar(11)) i1=i1+1
c...  exponent
      do 10 ie=i1,i2
      if (char1(ie:ie).eq.xchar(14) .or. 
     +    char1(ie:ie).eq.xchar(15)) go to 20
10    continue
      iexp=0
      go to 50
20    i2=ie-1
      iexp=1
      ie1=ie+1
      if (char1(ie1:ie1).eq.xchar(12))iexp=-1
      if (char1(ie1:ie1).eq.xchar(12).or.
     +    char1(ie1:ie1).eq.xchar(11))
     * ie1=ie1+1
      ibuff=0
      do 40 i=ie1,ie2
      do 30 j=1,10
      if (char1(i:i).eq.xchar(j)) go to 40
30    continue
      goto 100
40    ibuff=ibuff*10+j-1
      iexp=iexp*ibuff
c.... the number itself
 50   orep=.false.
      do 90 i=i1,i2
      if(char1(i:i).eq.xchar(13)) go to 80
      do 60 j=1,10
      if (char1(i:i).eq.xchar(j)) go to 70
60    continue
      goto 100
70    buf=buf*ten+ dble(j-1)
      go to 90
 80   if(orep)go to 100
      iexp=iexp+i-i2
      orep=.true.
90    continue
      buf = buf* dble(isign) * ten**iexp
      return
100   call errout(i)
      call caserr2
     * ('illegal character when reading floating point number')
      return
      end
      subroutine inpi(junke)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character * 132 char1,char2,char1c
      common/workc/char1,char2,char1c
c
c     mxtoken - the maximum number of tokens on an input line
c               for space separated tokens: 132/2 plus 1 to allow
c               for look ahead (132 is the line-length defined in
c               common/workc)
c
      integer mxtoken
      parameter(mxtoken=67)
      integer jrec, jump, istrt, inumb, iwidth, nend, nstart
      integer nline, noline, jwidth, nerr
      logical oswit, oterm, oflush
      common/work/jrec,jump,istrt(mxtoken),inumb(mxtoken),iwidth,
     + nend(mxtoken),nstart(mxtoken), nline,noline,jwidth,nerr,
     + oswit,oterm,oflush
c
      dimension xchar(12)
      data xchar /'0','1','2','3','4','5','6','7','8','9'
     1,'+','-'/
c.... subroutine for reading integers from the array char1,
c     starting at char1(istrt(jrec)) and going on for inumb(jrec))
c     elements. plus signs are ignored, the answer is accumulated
c     in jbuf and transferred to junke
      jbuf = 0
      jrec = jrec + 1
      if(jrec.gt.jump)go to 160
      n = inumb(jrec)
      ifact = 1
      ist=istrt(jrec)
      nstrt = ist + n - 1
      do 150 i = 1,n
      xtemp = char1(nstrt:nstrt)
      do 110 j=1,12
      if(xchar(j).eq.xtemp)go to 130
110   continue
120   call errout(nstrt)
      call caserr2('illegal character when reading integer')
130   if(j.lt.11)go to 140
      if(nstrt.ne.ist)go to 120
      if(j.ge.12)jbuf=-jbuf
      go to 160
140   jbuf=jbuf+(j-1)*ifact
      ifact = ifact * 10
150   nstrt=nstrt-1
160   junke=jbuf
      return
      end
c      subroutine inpwid(iwid)
c      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
c      implicit character *8 (z),character *1 (x)
c      implicit character *4 (y)
cINCLUDE(common/work)
cINCLUDE(common/iofile)
c      if((iwid.lt.1).or.(iwid.gt.132))call caserr2(
c     1 'illegal line width specified')
c      jwidth=iwid
c      write(iwr,10)jwidth
c10    format(///' input line width set to',i4,' characters')
c      return
c      end
      subroutine scann
      implicit real*8 (a-h,o-z),integer(i-n)
c
c     mxtoken - the maximum number of tokens on an input line
c               for space separated tokens: 132/2 plus 1 to allow
c               for look ahead (132 is the line-length defined in
c               common/workc)
c
      integer mxtoken
      parameter(mxtoken=67)
      integer jrec, jump, istrt, inumb, iwidth, nend, nstart
      integer nline, noline, jwidth, nerr
      logical oswit, oterm, oflush
      common/work/jrec,jump,istrt(mxtoken),inumb(mxtoken),iwidth,
     + nend(mxtoken),nstart(mxtoken), nline,noline,jwidth,nerr,
     + oswit,oterm,oflush
c
1      if(jrec.lt.jump)goto 999
      call input
      goto 1
999   return
      end
      subroutine rinpi(i)
      implicit real*8 (a-h,o-z),integer(i-n)
c     this subroutine reads an integer,
c     no matter how many cards it takes
       call scann
       call inpi(i)
      return
      end
      subroutine rchar(a)
      implicit real*8 (a-h,o-z),integer(i-n)
c     this subroutine reads a word of text
c     no matter how many cards it takes
      character *8 a
       call scann
       call inpa(a)
      return
      end
      subroutine rinpf(f)
      implicit real*8 (a-h,o-z),integer(i-n)
c     this subroutine reads an real,
c     no matter how many cards it takes
       call scann
       call inpf(f)
      return
      end
      function ytrunc(ztext)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character*(*) ztext
      ytrunc(1:4)=ztext(1:4)
      return
      end
      subroutine clenup
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      call clredx
      call whtps
      call shut
      return
      end
      subroutine wrt3(q,nword,iblk,num3)
c
c     this routine writes nword(s) of real data, commencing at
c     block iblk on unit num3, from the array q.
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension q(*)
      call search(iblk,num3)
      j=1
      k=nword
 20   if(k)30,30,10
 10   call put(q(j),min(k,511),num3)
      j=j+511
      k=k-511
      go to 20
30    return
      end
      function iposun(iunit)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      integer isel,iselr ,iselw,irep,ichek,ipos,
     *  nam,keep,istat,
     *  len,iop,iwhere,isync,
     *  iposf,keepf,
     *  istatf,lrecf,nrecf
      logical oform
      logical oedpr, oftpr, ogafsmode
      integer ioupd
      logical oadd
      common/disc/isel,iselr ,iselw,irep,ichek,ipos(maxlfn),
     *  nam(maxlfn),keep(maxlfn),istat(maxlfn),
     *  len(maxlfn),iop(maxlfn),iwhere(maxlfn),isync(maxlfn),
     *  iposf(maxfrt),keepf(maxfrt),
     *  istatf(maxfrt),lrecf(maxfrt),nrecf(maxfrt)
     * ,oform(maxfrt)
     *  ,oedpr(maxlfn),oftpr(maxfrt), ogafsmode,
     *  ioupd,oadd(maxlfn)
      if(iunit.le.0.or.iunit.gt.maxlfn)
     + call caserr2('invalid positioning requested')
      iposun = ipos(iunit)
      return
      end
      subroutine rdedx(q,nword,iblk,num3)
c
c     this routine reads nwords of real data, commencing at
c     block iblk on unit num3 into the array q.
c     note that nword has to be exactly what lies on disk
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension q(*)
      logical opg_root
      if (nword .eq. 0) return
      call search(iblk,num3)
      j=1
 20   if(j.gt.nword)go to 30
      call find(num3)
      call get(q(j),l)
      j=j+l
      go to 20
30    if(j-1.ne.nword) then
       write(6,*)
     + 'invalid no of words: requested,present = ', nword, j-1
       call caserr2('invalid number of words in rdedx')
      endif
      return
      end
      subroutine rdedx_less(q,nword,iblk,num3)
c...   allow nword to be less than actually read (overflow)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension q(*)
      if (nword .eq. 0) return
      call search(iblk,num3)
      j=1
 20   if(j.gt.nword)go to 30
      call find(num3)
      call get(q(j),l)
      j=j+l
      go to 20
30    return
      end
      subroutine rdedx_prec(q,nword,iblk,num3)
c...   read precisely nword (use extra buffer to accomplish this)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension q(*)
      dimension qbuff(511)
      if (nword .eq. 0) return
      call search(iblk,num3)
      j=1
 20   if(j.gt.nword)go to 30
      call find(num3)
      if (j+510.gt.nword) then
         call get(qbuff,l)
         l = min(l,nword-j+1)
         call dcopy(l,qbuff,1,q(j),1)
      else
         call get(q(j),l)
      end if
      j=j+l
      go to 20
30    return
      end
      subroutine wrt3s(q ,nword,num3)
c
c     this routine writes nword(s) of real data, commencing at
c     the current position/block on unit num3, from the array q.
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension q(*)
      j=1
      k=nword
20    if(k)30,30,10
10    call put( q(j),min(k,511),num3)
      j=j+511
      k=k-511
      go to 20
30     return
       end
      subroutine reads(q,nword,num3)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension q(*)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      j=1
10    if(j.gt.nword)go to 30
      call find(num3)
      call get(q(j),l)
      j=j+l
      go to 10
30    if(j-1.ne.nword) then
       write(iwr,*)
     + 'invalid no of words: requested,present = ', nword, j-1
       call caserr2('invalid number of words in reads')
      endif
      return
      end
      subroutine readsx(q,nword,nreply,num3)
c
c...  called from dirctb ; the overflow of g is intentional
c...  and is accounted for in david1 (jvl,2000/vrs many moons ago)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension q(*)
      j=1
10    if(j.gt.nword)go to 30
      call find(num3)
      call get(q(j),l)
      j=j+l
      go to 10
30    nreply=j-1
      return
      end

c-----------------------------------------------------------------------
c     DUMPFILE routines all concentrated HERE
c-----------------------------------------------------------------------
      subroutine qsector(mpos,ipos,iclass,ilen,op)
c
c...   query common sector - replacement of upack3  
c...   op = get or put
c...   it would be better to channel all operations on sector
c...   through an interface like this; but that is  bridge too far now
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer numdu, iblkdu, maxb, iblkla
      logical orevis
      real*8 apos
      common/sector/numdu,iblkdu,orevis(2),
     + apos(508),maxb,iblkla
       character*(*) op
c
      if (op.eq.'get') then
         if(mpos.lt.1.or.mpos.gt.508) then
             write(iwr,100) mpos
 100         format(/1x,'**************************',/
     *               1x,'problem with section ',i5,/
     *               1x,'**************************')
          call caserr2('attempting to use invalid dumpfile section get')
      end if
c
          call upack3(apos(mpos),ipos,iclass,ilen)
c
       else if (op.eq.'put') then
c
         orevis(1)=.false.
         if(mpos.lt.1.or.mpos.gt.508) then
             write(iwr,100) mpos
          call caserr2('attempting to use invalid dumpfile section put')
         end if
         if(iclass.lt.0.or.iclass.gt.2500) then
             write(iwr,101) mpos,iclass
 101         format(/1x,'*****************************************',/
     *               1x,'problem with section ',i5,' iclass ',i5,/
     *               1x,'*****************************************')
          call caserr2('attempting to use invalid dumpfile type (put)')
         end if
c
      apos(mpos)=pack3(ipos,iclass,ilen)
      else
         call caserr('wrong operation in qsector')
      end if
c
      return
      end
      subroutine secini(ibl,num)
c...   init dumpfile ; 0,0 => start afresh;
c      cf /sector/
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer numdu, iblkdu, maxb, iblkla
      logical orevis
      real*8 apos
      common/sector/numdu,iblkdu,orevis(2),
     + apos(508),maxb,iblkla
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      data izero/0/
c
      nav = lenwrd()
      mhead = 508 
      orevis(1)=.true.
c
      if (ibl.eq.0.and.num.eq.0) then
         mhead = mhead * nav
         call setsto(mhead,izero,apos)
      else
         mhead = mhead + 2/nav
         numdu=num
         iblkdu=ibl
         call rdedx(apos,mhead,iblkdu,numdu)
      end if
c
      return
      end
      logical function sector_block(ibl,num)
c...   check if current (collection of) blocks qualify as dumpfile block
      integer numdu, iblkdu, maxb, iblkla
      logical orevis
      real*8 apos
      common/sector/numdu,iblkdu,orevis(2),
     + apos(508),maxb,iblkla
c
      sector_block = .false.
      call search(ibl,num)
      call fgett(apos,nword,iunit)
c
      mhead = 508 + 2/lenwrd()
      if (nword.ne.mhead) return
c
      if(maxc.gt.maxb.or.maxc.lt.1) return
      do   i=1,508 
         call qsector(i,m,itype,n,'get')
         if (m.gt.0) then
            if (n.eq.0.or.(m+n).gt.maxc.or.itype.eq.0) return 
         end if
      end do
c
      sector_block = .true.
c
      return
      end

      subroutine secinf(ibl,iunit,max,maxbb)
c...   return dumpfile info
      implicit none
      integer ibl,iunit,max,maxbb
      integer numdu, iblkdu, maxb, iblkla
      logical orevis
      real*8 apos
      common/sector/numdu,iblkdu,orevis(2),
     + apos(508),maxb,iblkla
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      ibl = iblkdu
      iunit = numdu
      max = iblkla
      maxbb = maxb
CMR   write(iwr,*) 'CMR secinf says: ',iblkdu,numdu,iblkla,maxb
c
      return
      end

      subroutine secsum
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
      integer numdu, iblkdu, maxb, iblkla
      logical orevis
      real*8 apos
      common/sector/numdu,iblkdu,orevis(2),
     + apos(508),maxb,iblkla
       call secinf(iblkdu,numdu,iblkla,maxbb)
       write(iwr,600)yed(numdu),iblkdu,iblkla,maxb
 600   format(//' *summary of dumpfile on ',a4,' at block',i8/2h *,/
     * ' *current length=',i18,' blocks'/2h *,/
     * ' *maximum length=',i18,' blocks'/2h *,/
     * ' *section type   block  length')
      do i=1,508
         call qsector(i,ipos,iclass,ilen,'get')
         if(ipos .gt. 0) then
            m=iblkdu+ipos
            write(iwr,500)i,iclass,m,ilen
         endif
      end do
      return
 500  format(2h *,i7,i5,2i8)
      end
      subroutine secloc(mpos,oexist,iblock)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      call  secinf(iblkdu,idum,idum,idum)
      call qsector(mpos,ipos,iclass,ilen,'get')
       oexist=.false.
       if(ipos.gt.0) then
        oexist=.true.
        iblock=ipos+iblkdu
       endif
      return
      end
      subroutine secget(mpos,mtype,iblock)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      call secinf(iblkdu,numdu,idum,idum)
      call qsector(mpos,ipos,iclass,ilen,'get')
c
       if(ipos.lt.1)then
       write(iwr,100) mpos,mtype,iblock
 100   format(/1x,'*********************',/
     *         1x,'problem with section ',/
     *         1x,'*********************',//
     *         1x,'mpos, mtype, block ',3(1x,i5)/)
       call secsum
       call caserr2('attempting to retrieve undefined dumpfile section')
       endif
       if(mtype.eq.0)mtype=iclass
       if(iclass.ne.mtype) then
       write(iwr,100) mpos,mtype,iblock
        call caserr2(
     * 'retrieved dumpfile section of wrong type')
       endif
       iblock=ipos+iblkdu
       call search(iblock,numdu)
      return
      end

      subroutine revind
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer numdu, iblkdu, maxb, iblkla
      logical orevis
      real*8 apos
      common/sector/numdu,iblkdu,orevis(2),
     + apos(508),maxb,iblkla
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical odumpm, ooprnt, opun7, omax, ologf
      integer maxbrd, logf, iofsp
      common /utilc/ odumpm,ooprnt,opun7,omax,maxbrd,ologf,
     +               logf,iofsp(4)
c

      
      if(orevis(1)) then
       call clredx
      else
       nav = lenwrd()
       call search(iblkdu,numdu)
       mhead = 508 + 2/nav
       call wrt3s( apos, mhead, numdu )
       orevis(1)=.true.
       if (ooprnt) then
        write(iwr,10) 
10      format (/1x,'****** output index block of dumpfile ***')
       endif
      endif
      return
      end
      subroutine secput(mpos,mtype,length,iblock)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer numdu, iblkdu, maxb, iblkla
      logical orevis
      real*8 apos
      common/sector/numdu,iblkdu,orevis(2),
     + apos(508),maxb,iblkla
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      parameter (maxdlen=16*65536)
      call qsector(mpos,ipos,iclass,ilen,'get')
c
      if(ipos.le.0 .or. length.gt.ilen)then
c
c  create/expand
c
         if(ipos.gt.0)write(iwr,10) mpos
 10      format(/1x,'**** Warning ****'/
     *        1x,'**** Expanding section ',i5,' on the dumpfile')
c
         m=length+iblkla
         if(m.gt.maxb)call caserr2(
     *   'attempting to expand dumpfile beyond maximum allowed size')
          iblock=iblkdu+iblkla
         call qsector(mpos,iblkla,mtype,length,'put')
         iblkla=m
      else
         call qsector(mpos,ipos,mtype,ilen,'put')
         iblock=ipos+iblkdu
      endif
c
      if(length.gt.maxdlen.or.ipos.gt.maxdlen*16) then
c
         write(iwr,11) mpos,mtype,length
11       format(' trying dumpfile section mpos,mtype,length',3i10)
         call caserr2('dumpfile section-length overflow in secput')
      end if  

      call search(iblock,numdu)
c
      return
      end
      subroutine secdrp(mpos)
c Disable section mpos by setting iclass to zero (see sectst).
c AdM JvL 1991
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      call qsector(mpos,ipos,iclass,ilen,'get')
c
      if(ipos.gt.0) then
         iclass = 0
         call qsector(mpos,ipos,iclass,ilen,'put')
      endif
c
      return
      end
      subroutine sectst(mpos,iretrn)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      call qsector(mpos,ipos,iclass,ilen,'get')
c
      if(ipos.gt.0. and. iclass.gt.0) then
         iretrn=1
      else
         iretrn=0
      endif
c
      return
      end

      subroutine maxset
      implicit real*8  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      integer numdu, iblkdu, maxb, iblkla
      logical orevis
      real*8 apos
      common/sector/numdu,iblkdu,orevis(2),
     + apos(508),maxb,iblkla
      call inpi(maxb)
      if(maxb.lt.iblkla)call caserr('invalid dumpfile size specified')
      write(iwr,300)maxb
 300  format(/' maximum size of dumpfile now',i8,' blocks')
      orevis(1)=.false.
      return  
      end     

      subroutine sumvecs
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character*80 charv
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical odiis, optester, odynamic, omaxcyc,onoor,otestd
      integer maxcyc,mconv,nconv,npunch,icoupl,ifolow,irotb
      integer iter,kcount,iextin,iterv,idiisf
      real*8 accdi1,accdi2,dmpcut,acurcy,en,etot
      real*8 ehf,ehf0,diff,rshift,exttol,dmptol,vshtol
      real*8 damp,damp0,diffd,diffp,de,deavg,diffsp
      real*8 ek, vir,diffpp
      common/scfopt/maxcyc,mconv,nconv,npunch,accdi1,accdi2,odiis,
     +      icoupl,ifolow,irotb,dmpcut,acurcy,en,etot,ehf,ehf0,diff,
     +      iter,kcount,rshift,exttol,dmptol,vshtol,iextin,
     +      iterv,damp,damp0,diffd,diffp,diffpp,de,deavg,diffsp,
     +      ek,vir,idiisf,optester,odynamic,omaxcyc,onoor,otestd
c
      common/junkcXX/zcom(19),ztit(10)
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
c
      data m3,m29/3,29/
c
c ---- read a set of vectors off the dumpfile, and write
c      a summary to the output file
c
      call secinf(idum,numdu,idum,idum)
      write(iwr,100)
100   format(/1x,'*summary of vector sections'/1x,'*'/
     +        1x,'*section',4x,'type',6x,'created:',11x,'title:')
      do i = 1, 508
c
         call qsector(i,ipos,iclass,ilen,'get')
c
       if (ipos.gt.0) then
        if(iclass.eq.m3) then
c
          m=iblkdu + ipos
          call secget(i,m3,k)
          call rdchr(zcom,m29,k,numdu)
c
          j = 1
          do loop = 1,10
           charv(j:j+7) = ztit (loop)
           j = j + 8
          enddo
          last = lstchr(charv)
c
          write(iwr,200)i,zcom(5),zcom(3),zcom(2),charv(1:last)
200       format(1x,'*',i7,4x,a8,1x,a8,' on ', a6, 2x, a)
c
        endif
       endif
      enddo
c
      return
      end
c
c----------------------------------------------------------------------
c

      subroutine tfsqc(v,q,t,m,n,ndim)
c
c     ----- back transform the square matrix q with v -----
c     ----- t is scratch
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension v(ndim,*),q(ndim,*),t(ndim,*)
c
      nmm=ndim*m
      call vclr(t,1,nmm)
      call mxmb(q,1,ndim,v,1,ndim,t,1,ndim,n,m,m)
      call dcopy(nmm,t,1,v,1)
      return
      end
       subroutine mult1a(a,q,ilifq,ncore,h,iky,nbasis)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      dimension a(*),q(*),h(*),ilifq(*),iky(*)
      dimension p(maxorb),r(maxorb)
c...   a=q(transpose) * h * q
c...   a and h stored in triangle form
      do 1 j=1,ncore
      mm=ilifq(j)
      p(1)=h(1)*q(mm+1)
      do 2 i=2,nbasis
      m=iky(i)
      call daxpy(i-1,q(mm+i),h(m+1),1,p,1)
      p(i)=ddot(i,h(m+1),1,q(mm+1),1)
 2    continue
      m=iky(j)
      do 1 i=1,j
      a(m+i)=ddot(nbasis,p,1,q(ilifq(i)+1),1)
 1    continue
      return
      end
       subroutine writel(p,newbas)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
       dimension p(*)
      ind(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
      m5=12
      if(oprint(20)) m5=8
      m=1
      n=m5
   6  if(newbas.lt.m)return
      if(n.gt.newbas)n=newbas
      if(.not.oprint(20))write(iwr,200)(i,i=m,n)
      if(oprint(20))write(iwr,100)(i,i=m,n)
100   format(//3x,8i14)
200   format(//12i9)
      write(iwr,101)
101   format(/)
      do 1 j=1,newbas
      if(oprint(20))write(iwr,102)j,(p(ind(i,j)),i=m,n)
      if(.not.oprint(20))write(iwr,202)j,(p(ind(i,j)),i=m,n)
 1    continue
102   format(7x,i3,8f14.7)
202   format(1x,i3,12f9.4)
      m=m+m5
      n=n+m5
      goto 6
      end
      function locatc(label,nf,itext)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character * (*) label,itext
      dimension label(*)
      do 1 i=1,nf
      if(label(i).eq.itext)go to 2
 1    continue
      locatc=0
      return
 2    locatc=i
      return
      end
      subroutine filprn(nfile,iblk,lblk,notape)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
      dimension iblk(*),lblk(*),notape(*)
      write(iwr,1)(yed(notape(i)),i=1,nfile)
      write(iwr,2)(iblk(i),i=1,nfile)
      write(iwr,3)(lblk(i),i=1,nfile)
1     format(/9x,'ddnames ',20(2x,a4))
2     format(' starting blocks',20i6)
3     format(' terminal blocks',20i6)
      return
      end
      subroutine prtri(d,n)
c
c     ----- print out a triangular matrix -----
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
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
      dimension d(*),dd(12)
      mmax = 12
      if (oprint(20)) mmax=7
      imax = 0
  100 imin = imax+1
      imax = imax+mmax
      if (imax .gt. n) imax = n
      write (iwr,9008)
      if(oprint(20)) write (iwr,9028) (i,i = imin,imax)
      if(.not.oprint(20)) write (iwr,8028) (i,i = imin,imax)
      do 160 j = 1,n
      k = 0
      do 140 i = imin,imax
      k = k+1
      m = max(i,j)*(max(i,j)-1)/2 + min(i,j)
  140 dd(k) = d(m)
      if(oprint(20)) write (iwr,9048) j,(dd(i),i = 1,k)
      if(.not.oprint(20)) write (iwr,8048) j,(dd(i),i = 1,k)
  160 continue
      if (imax .lt. n) go to 100
      return
 9008 format(/)
 9028 format(6x,7(6x,i3,6x))
 9048 format(i5,1x,7f15.10)
 8028 format(6x,12(3x,i3,3x))
 8048 format(i5,1x,12f9.4)
      end
      subroutine prsq(v,m,n,ndim)
c
c     ----- print out a square matrix -----
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      dimension v(ndim,*)
      max = 12
      if (oprint(20)) max=7
      imax = 0
  100 imin = imax+1
      imax = imax+max
      if (imax .gt. m) imax = m
      write (iwr,9008)
      if(.not. oprint(20)) write (iwr,8028) (i,i = imin,imax)
      if( oprint(20)) write (iwr,9028) (i,i = imin,imax)
      write (iwr,9008)
      do 120 j = 1,n
      if( oprint(20)) write (iwr,9048) j,(v(j,i),i = imin,imax)
      if(.not.oprint(20)) write (iwr,8048) j,(v(j,i),i = imin,imax)
120   continue
      if (imax .lt. m) go to 100
      return
 9008 format(1x)
 9028 format(6x,7(6x,i3,6x))
 9048 format(i5,1x,7f15.10)
 8028 format(6x,12(3x,i3,3x))
 8048 format(i5,1x,12f9.5)
      end
      subroutine trianc(r,a,mrowr,n)
      implicit real*8  (a-h,o-z)
      dimension r(mrowr,*),a(*)
c... convert from square r to lower triangle a
      k = 0
      do 30 i = 1 , n
         do 20 j = 1 , i
            k = k + 1
            a(k) = r(i,j)
 20      continue
 30   continue
      return
      end
      subroutine openda
c
c     ----- open a direct access file -----
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
c --- initialise the zero block
c
      call wrrec1(id,id,id,id,c,c,d,d,d,c,c)
      call revind
      return
      end
      subroutine closda
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      character*10 charwall
      cpu=cpulft(1)
      write(iwr,1)cpu,charwall()
 1    format(///
     *' end of  G A M E S S   program at ',f12.2,
     *' seconds',a10,' wall',/)
c
      oprintv = .true.
c
      call revind
      call secsum
      if (oprintv) call sumvecs
      call clenup
      call timout
      return
      end
      subroutine setscm(k)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer nmaxly, ntotly, icurly, isecly
      logical oprintm
      common /segm/ nmaxly,ntotly,icurly,isecly(100),
     +              oprintm
c
c
c -----  core management routines for gamess
c
c --- this routine returns the address in k of the next scm block
c --- of core as an index in the array x
c
      k=ntotly+1
      return
      entry cmem(load)
c
c --- set load to the present core usage
c
      load=ntotly
      return
      entry setc(need)
c
c --- update entries in the common/ovly/ header block
c --- and assign another block of core if necessary
c --- totla amount of core required is need
c --- if need is less than the current amount allocated
c --- then a section of core is deleted
c
 10   if(need.gt.nmaxly) goto 900
      if(need.gt.ntotly) goto 20
      if(need.eq.ntotly) goto 30
c
c --- current requirement is less than that allocated
c --- so delete the last section of core and try again
c
      if(icurly.eq.1) goto 30
      icurly=icurly-1
      ntotly=ntotly-isecly(icurly+1)
      if(opg_root().and.oprintm)
     +   write(iwr,*)'setc: freed to ',ntotly
      goto 10
c
c --- current requirement is greater than that allocated
c -- so add a block of core
c
 20   if(icurly.ge.100) goto 920
      icurly=icurly+1
      isecly(icurly)=need-ntotly
      ntotly=need
      if(opg_root().and.oprintm)
     +   write(iwr,*)'setc: allocated to',ntotly
      go to 30
 900  write(iwr,910) need,nmaxly
 910  format(//1x,'core management :- core overflow'/,
     1       ' need ',i10,'     got',i10,//)
      go to 930
 920  write(iwr,940)
 940  format(//1x,'core management :- too many sections'//)
 930  call caserr2('insufficient memory allocated')
c
c --- can keep the current allocation
c
 30   return
c
c ---eroor messages
c
      end
      function loccm()
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      integer nmaxly, ntotly, icurly, isecly
      logical oprintm
      common /segm/ nmaxly,ntotly,icurly,isecly(100),
     +              oprintm
c
      loccm=ntotly
      return
      end
      subroutine init_segm
      implicit none
c
      integer nmaxly, ntotly, icurly, isecly
      logical oprintm
      common /segm/ nmaxly,ntotly,icurly,isecly(100),
     +              oprintm
c
      integer loop
      icurly = 1
      ntotly = 0
      do 120 loop = 1 , 100
         isecly(loop) = 0
 120  continue
      end
      subroutine setlab
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *3 char3,char3i
      character *4 bfnam
      character *1 atnam1
      character *2 atnam2
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
c from March 97, amass is used only for user-input
c mass values
c czanr contains a copy of czan made in basis.m after reorganisation and
c before pseudo potential; so real charges in real order 
c beware, infob is explicit in direct.m and  dirrpa.m
c
      real*8 czann, czanr, cat, amass, cnew
      integer nonsym, map80
      logical ozmat
      common/infob/czanr(maxat),czann(maxat),cat(3,maxat),amass(maxat),
     +             cnew(maxat3,3),nonsym,map80(maxat),ozmat

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
      dimension bfnam(35),atnam1(104),atnam2(104)
      data bfnam/
     + 's','x','y','z',
     + 'xx','yy','zz','xy','xz','yz',
     + 'xxx','yyy','zzz','xxy','xxz','xyy','yyz','xzz',
     + 'yzz','xyz',
     + 'xxxx','yyyy','zzzz','xxxy','xxxz','yyyx','yyyz',
     + 'zzzx','zzzy','xxyy','xxzz','yyzz','xxyz','yyxz',
     + 'zzxy' /
      data atnam1 /' ','h',
     +             'l','b',' ',' ',' ',' ',' ','n',
     +             'n','m','a','s',' ',' ','c','a',
     +             ' ','c','s','t',' ','c','m','f','c',
     +        'n','c','z','g','g','a','s','b','k',
     +             'r','s',' ','z','n','m','t','r','r',
     +        'p','a','c','i','s','s','t',' ','x',
     +'c','b','l','c','p','n','p','s','e','g','t',
     +'d','h','e','t','y','l','h','t',' ','r','o',
     +'i','p','a','h','t','p','b','p','a','r',
     +'f','r','a','t','p',' ','n','p','a','c','b',
     +'c','e','f','m','n','l',
     +             ' '/
c
      data atnam2 / 'h', 'e',
     +              'i', 'e', 'b', 'c', 'n', 'o', 'f', 'e',
     +              'a', 'g', 'l', 'i', 'p', 's', 'l', 'r',
     +              'k', 'a', 'c', 'i', 'v', 'r', 'n', 'e', 'o',
     +         'i', 'u', 'n', 'a', 'e', 's', 'e', 'r', 'r',
     +              'b', 'r', 'y', 'r', 'b', 'o', 'c', 'u', 'h',
     +         'd', 'g', 'd', 'n', 'n', 'b', 'e', 'i', 'e',
     + 's', 'a', 'a', 'e', 'r', 'd', 'm', 'm', 'u', 'd', 'b',
     + 'y', 'o', 'r', 'm', 'b', 'u', 'f', 'a', 'w', 'e', 's',
     + 'r', 't', 'u', 'g', 'l', 'b', 'i', 'o', 't', 'n',
     + 'r', 'a', 'c', 'h', 'a', 'u', 'p', 'u', 'm', 'm', 'k',
     + 'f', 's', 'm', 'd', 'o', 'w',
     +              ' '/
c
      n = 0
      do 100 ii = 1,nshell
      iat = katom(ii)
      cznuc = czanr(iat)
      j =  nint(cznuc)
      if (j .gt.103. or.j . eq. 0) j = 104
      mini = kmin(ii)
      maxi = kmax(ii)
      if(nat.ge.100) then
      do 110 i = mini,maxi
      n = n+1
      zbflab(n)(1:3)=char3i(iat)
      zbflab(n)(4:4)=atnam1(j)
      zbflab(n)(5:6)=atnam2(j)
      zbflab(n)(7:10)=bfnam(i)
  110 continue
      else
      do 120 i = mini,maxi
      n = n+1
      char3=char3i(iat)
      zbflab(n)(1:2)=char3(1:2)
      zbflab(n)(3:3)=atnam1(j)
      zbflab(n)(4:5)=atnam2(j)
      zbflab(n)(6:9)=bfnam (i)
      zbflab(n)(10:10)= ' '
  120 continue
      endif
  100 continue
      return
      end
      function char3i(i)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *3 c,char3i
      character *1 digit(0:9)
      data c/' '/
      data digit/
     *'0','1','2','3','4','5','6','7','8','9'/
      num=iabs(i)
      char3i='   '
      l=3
      do 1 j=3,1,-1
      new=num/10
      n=num-10*new
      if(n.gt.0)l=j
      c(j:j)=digit(n)
 1    num=new
      goto (2,3,4),l
2     char3i(1:3) = c(1:3)
      go to 5
3     char3i(1:2) = c(2:3)
      go to 5
4     char3i(1:1) = c(3:3)
5     return
      end
      subroutine intr(core)
c
c     ----- calculate atom-atom distances -----
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
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
      real*8 toang, ams
      integer ifau, nosymm, iseczz
      common/phycon/toang(30),ams(54),ifau,nosymm,iseczz
c
      real*8 scalel, vbond1, vbond2
      logical opoint1,ochkdst,ofatal
      integer mxbnd, nbond1, nbond2, mbond1, mbond2
      parameter (mxbnd = 1000)
c
      common/structchk/ vbond1(mxbnd),mbond1(2,mxbnd),
     +                  vbond2(mxbnd),mbond2(2,mxbnd),
     +                  scalel, nbond1, nbond2, opoint1,
     +                  ochkdst, ofatal
c
      dimension r(maxat)
c
      dimension core(*)
      data done/1.0d0/
c
c     write(6,*)'oprint',oprint(21)

      if (nat.eq.1)return
      if (oprint(21).and..not.ochkdst)return
      if (oprint(56)) then
       ipass = 1
       fac = done
       if(nprint.ne.-5)write (iwr,9028)
  100  continue
       max = 0
  120  min = max+1
       max = max+7
       if (max .gt. nat) max = nat
       if(nprint.ne.-5)write (iwr,9048)
       if(nprint.ne.-5)write (iwr,9068)(zaname(j),j= min,max)
       if(nprint.ne.-5)write (iwr,9048)
       do i = 1,nat
        do  j = min,max
        rr = (c(1,i)-c(1,j))**2+(c(2,i)-c(2,j))**2+(c(3,i)-c(3,j))**2
        r(j) = dsqrt(rr)*fac
        enddo
        if(nprint.ne.-5)write(iwr,9088)i,zaname(i),
     *  (r(j),j=min,max)
       enddo
       if (max .lt. nat) go to 120
       if (ipass .lt. 2) then
        ipass = 2
        fac = toang(1)
        if(nprint.ne.-5)write (iwr,9008)
        go to 100
        endif
      endif
c
      if (ochkdst) then
       if(lpseud.eq.0) then
        call struct_check(nat,c,czan,nuct,core)
       else
        call struct_check(nat,c,symz,nuct,core)
       endif
      endif
c
      if(nprint.ne.-5) then
       if(lpseud.eq.0) then
        call struct(nat,c,czan,nuct,core)
       else
        call struct(nat,c,symz,nuct,core)
       endif
      endif
c
      return
 9008 format(/40x,'internuclear distances (angs.)'/40x,30('-'))
 9028 format(/40x,'internuclear distances ( a.u.)'/40x,30('-'))
 9048 format(/)
 9068 format(17x,7(4x,a8,3x))
 9088 format(i3,2x,a8,2x,7f15.7)
      end
      subroutine timit(index)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
      tim=cpulft(1)
      tx = tim-ti
      ti = tim
      if(index.eq.0.or.index.eq.3.or.index.eq.-5)return
      write(iwr,9008)tx,tim
 9008 format(/1x,33('-')/
     *'     elapsed time = ',f8.2,' secs'/
     *'     total time   = ',f8.2,' secs'/1x,33('-')/)
      return
      end
      subroutine texit(ncall,nrest)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
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
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
      call timit(ncall)
      if (tim .lt. timlim) return
      irest = nrest
      itask(mtask) = irest
      call revise
      return
      end
      subroutine gamgen
c        *****  computes and tabulates f0(x) to f5(x)           *****
c        *****  in range x = -0.24 to x = 26.4                  *****
c        *****  in units of x = 0.08                            *****
c        *****  used by the one electron routine auxg and by    *****
c        *****  the two electron integral routines              *****
c        *****  the table is generated only once for each entry *****
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
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
      common/inttab/c(1000,6)
      dimension ppp(350)
      data pt184,pt5/ 0.184d0,0.50d0/
      data six,tenm7/6.0d0,1.0d-20 /
      data four,two,done/4.0d0,2.0d0,1.0d0/
c     data m22,mword/22,6000/
      data pt886/0.8862269254527d0/
      q=-done
      do 30 mm=1,6
      m=mm-1
      q=q+done
      qqq = -0.24d0
      do 20 i=1,340
      qqq = qqq+0.08d0
      a=q
c        *****  change limit of approximate solution.           *****
      if(qqq-15.0d0) 1,1,10
    1 a=a+pt5
      term=done/a
      ptlsum=term
      do 2 l=2,50
      a=a+done
      term=term*qqq/a
      ptlsum=ptlsum+term
      if( dabs(term/ptlsum)-tenm7)3,2,2
    2 continue
    3 ppp(i)=pt5*ptlsum* dexp(-qqq)
      go to 20
   10 b=a+pt5
      a=a-pt5
      approx=pt886/(dsqrt(qqq)*qqq**m)
      if(m.eq.0) go to 13
      do 12 l=1,m
      b=b-done
   12 approx=approx*b
   13 fimult=pt5* dexp(-qqq)/qqq
      fiprop=fimult/approx
      term=done
      ptlsum=term
      notrms=qqq
      notrms=notrms+m
      do 14 l=2,notrms
      term=term*a/qqq
      ptlsum=ptlsum+term
      if( dabs(term*fiprop/ptlsum)-tenm7)15,15,14
   14 a=a-done
   15 ppp(i)=approx-fimult*ptlsum
   20 continue
      do 30 i=1,333
      j=i+2
      c(i,mm)=ppp(j)
      c(i+333,mm)=ppp(j+1)-ppp(j)
      temp1=-two*ppp(j)+ppp(j+1)+ppp(j-1)
      temp2=six*ppp(j)-four*ppp(j+1)-four*ppp(j-1)+ppp(j-2)+ppp(j+2)
   30 c(i+666,mm) = (temp1-pt184*temp2)/six
c        *****  write out interpolation table                   *****
c      call secput(isect(502),m22,lensec(mword),iblk22)
c      call wrt3(c,mword,iblk22,idaf)
      return
      end
      subroutine wrrec1(n1,n2,n3,n4,b,d,d1,d2,d3,e,f)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      dimension h(4*maxat*3+3+4)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
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
      dimension b(*),d(*),e(*),f(*),ih(2)
      equivalence (ih(1),h(1))
      data m15/15/
      nav = lenwrd()
      mxcen = nat * 3
      maxtot = mxcen*4+3
      call dcopy(mxcen,b,1,h,1)
      call dcopy(mxcen,d,1,h(  mxcen+1),1)
      call dcopy(mxcen,e,1,h(2*mxcen+1),1)
      call dcopy(mxcen,f,1,h(3*mxcen+1),1)
      itemp = 4*mxcen+1
      h(itemp  ) = d1
      h(itemp+1) = d2
      h(itemp+2) = d3
      imax = nav*maxtot
      ih(imax+1) = n1
      ih(imax+2) = n2
      ih(imax+3) = n3
      ih(imax+4) = n4
      mach4 = 12*nat +3 +4/nav
      call secput(isect(493),m15,lensec(mach4),iblk15)
      call wrt3(h,mach4,iblk15,idaf)
      return
      end
      subroutine rdrec1(n1,n2,n3,n4,b,d,d1,d2,d3,e,f)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      dimension h(4*maxat*3+3+4)
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
      dimension b(*),d(*),e(*),f(*),ih(2)
      equivalence (ih(1),h(1))
      data m15/15/
      nav = lenwrd()
      mach4 = 12*nat +3 +4/nav
      mxcen = nat * 3
      maxtot = mxcen*4+3
      call secget(isect(493),m15,iblk15)
      call rdedx(h,mach4,iblk15,idaf)
      call dcopy(mxcen,h,1,b,1)
      call dcopy(mxcen,h(  mxcen+1),1,d,1)
      call dcopy(mxcen,h(2*mxcen+1),1,e,1)
      call dcopy(mxcen,h(3*mxcen+1),1,f,1)
      itemp=4*mxcen+1
      d1 = h(itemp  )
      d2 = h(itemp+1)
      d3 = h(itemp+2)
      imax = nav*maxtot
      n1 = ih(imax+1)
      n2 = ih(imax+2)
      n3 = ih(imax+3)
      n4 = ih(imax+4)
      return
      end

      function lensec(nword)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c
c     lensec returns the number of ATMOL blocks corresponding to nword 
c     words of data for all nword greater than zero. For nword equals 
c     zero it returns one.
c     Valid inputs are all non-negative integers.
c
c     Proof:
c
c     Define m = the block size which is 511 words
c
c     Assume n = i*m, i > 0 then
c     lensec = (i*m-1)/m+1
c            = i-1      +1
c            = i
c     Assume n = i*m+j, i >= 0 and 0 < j < m
c     lensec = (i*m+j-1)/m+1
c            = (i*m)/m+(j-1)/m+1
c            = i+1
c     Assume n = 0
c     lensec = -1/m+1
c            = 1
c
      lensec=(nword-1)/511+1
      return
      end

      function lenwrd()
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      parameter (nav = 1)
      lenwrd = nav
      return
      end
      subroutine struct(n, c, az, nuct, bl)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
      integer p
      dimension c(3,*),az(*),nuct(*)
      dimension bl(*)
      dimension v(8),m(4,8),ilifz(maxat)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
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
c approximate covalent radius table in au
      real*8 cov
      common/coval/cov(100)
      character*7 fnm
      character*6 snm
      data fnm,snm/"util1.m","struct"/
      data pi /3.14159265358979d0/
c  covalent radii (mostly from coulson & mcweeny) for atoms up to argon
c      data cov /0.57d0, 2.50d0,
c     * 2.53d0, 2.00d0, 1.53d0, 1.46d0, 1.42d0, 1.38d0, 1.34d0,3.00d0,
c     * 2.91d0, 2.69d0, 2.46d0, 2.23d0, 2.08d0, 1.93d0, 1.87d0,3.50d0,
c     * 85 * 3.0d0 /
      data done,two/1.0d0,2.0d0/
      data zsadd/'saddle'/
      data dsmall/10000.0d0/
c
c     ----- get core memory
c
      length=n*n
      i10 = igmem_alloc_inf(length,fnm,snm,'i10',IGMEM_DEBUG)
      ismall=0
      jsmall=0
c
      do 200 i=1,n
 200  ilifz(i)=(i-1)*n+i10-1
c
      scale = 1.0d0
      toler = 0.5d0
      if(zruntp.eq.zsadd)scale=two
      if(oprint(56)) scale =3.0d0
      call vclr(bl(i10),1,length)
      nm1=n-1
      do 21 i=1,nm1
      iz= nint(az(i))
      if(nuct(i).gt.3)go to 21
      covi=two
      if (iz .gt. 0 .and. iz .le. 100) covi=cov(iz)
      ip1=i+1
      do 20 j=ip1,n
      if(nuct(j).gt.3)go to 20
      ij=i+ilifz(j)
      bl(ij)=
     +  dsqrt((c(1,j)-c(1,i))**2+(c(2,j)-c(2,i))**2+(c(3,j)-c(3,i))**2)
      jz= nint(az(j))
      if (bl(ij).lt.dsmall) then
          ismall=i
          jsmall=j
          dsmall=bl(ij)
      endif
      covj=two
      if (jz .gt. 0 .and. jz .le. 100) covj=cov(jz)
      if (bl(ij) .gt. (covi+covj)*scale+toler) bl(ij)=-bl(ij)
      bl(j+ilifz(i))=bl(ij)
20    continue
21    continue
      write (iwr,1001)
 1001 format(//40x,31('=')/
     *40x,'bond lengths in bohr (angstrom)'/40x,31('='))
c
      max=7
      p=0
      do 30 i=1,nm1
      ip1=i+1
      do 30 j=ip1,n
      ij=i+ilifz(j)
      if (bl(ij) .le. 0.0d0) go to 30
      p=p+1
      if (p .le. max) go to 25
      call prtbuf(n,max,2, m, v)
      p=1
25    m(1,p)=i
      m(2,p)=j
      v(p)=bl(ij)
30    continue
      if (p .gt. 0) call prtbuf(n,p,2, m, v)
      if (p .eq. 0) write (iwr,1009)
1009  format (/5x,'--- none ---')
      if (dsmall.lt.0.5) then
          write(iwr,1200) ismall,jsmall,dsmall
      endif
 1200 format(//10x,55('+')/10x,
     *'WARNING: shortest bond length is smaller than 0.5 bohr!'/10x,
     *'atoms ',I5,2x,I5,F10.5/,10x,55('+'))
      write (iwr,1002)
1002  format (//40x,11('=')/40x,'bond angles'/40x,11('=')/)
c
      max=6
      p=0
      do 60 j=1,nm1
      do 50 i=1,n
      ij=i+ilifz(j)
      if (bl(ij) .le. 0.0d0.or. j .eq. i) go to 50
      jp1=j+1
      do 40 k=jp1,n
      ik=i+ilifz(k)
      jk=j+ilifz(k)
      if (bl(ik) .le. 0.0d0.or. k .eq. i) go to 40
      p=p+1
      if (p .le. max) go to 35
      call prtbuf(n,max,3, m, v)
      p=1
35    m(1,p)=j
      m(2,p)=i
      m(3,p)=k
      cosine=(bl(ij)**2+bl(ik)**2-bl(jk)**2)/(two*bl(ij)*bl(ik))
      if ( dabs(cosine) .gt. done) cosine= dsign(done,cosine)
      v(p)=(180.0d0/pi)*dacos(cosine)
40    continue
50    continue
60    continue
      if (p .gt. 0) call prtbuf(n,p,3, m, v)
      if (p .eq. 0) write (iwr,1009)
      write (iwr,1003)
1003  format (//40x,15('=')/40x,'dihedral angles'/40x,15('='))
c
      max=5
      p=0
      do 100 i=1,nm1
      ip1=i+1
      do 90 j=ip1,n
      ij=i+ilifz(j)
      if (bl(ij) .le. 0.0d0) go to 90
c  i and j are bonded.  construct unit vector from i to j.
      a=bl(ij)
      pij=(c(1,j)-c(1,i))/a
      qij=(c(2,j)-c(2,i))/a
      rij=(c(3,j)-c(3,i))/a
      do 80 k=1,n
      ik=i+ilifz(k)
      if (bl(ik) .le. 0.0d0.or. k .eq. i .or. k .eq. j) go to 80
c  i and k are bonded.  construct unit vector in jik plane, normal
c  to ij.
      pik=c(1,k)-c(1,i)
      qik=c(2,k)-c(2,i)
      rik=c(3,k)-c(3,i)
      dot=pik*pij + qik*qij + rik*rij
      pik=pik-dot*pij
      qik=qik-dot*qij
      rik=rik-dot*rij
      a=dsqrt(pik**2+qik**2+rik**2)
      if (a .lt. 1.0d-4) go to 80
      pik=pik/a
      qik=qik/a
      rik=rik/a
      do 70 l=1,n
      jl=j+ilifz(l)
      if (bl(jl) .le. 0.0d0.or. l .eq. i .or. l .eq. k) go to 70
c  j and l are bonded.  construct unit vector in ijl plane, normal
c  to ij.
      pjl=c(1,l)-c(1,j)
      qjl=c(2,l)-c(2,j)
      rjl=c(3,l)-c(3,j)
      dot=pjl*pij + qjl*qij + rjl*rij
      pjl=pjl-dot*pij
      qjl=qjl-dot*qij
      rjl=rjl-dot*rij
      a=dsqrt(pjl**2+qjl**2+rjl**2)
      if (a .lt. 1.0d-4) go to 80
      pjl=pjl/a
      qjl=qjl/a
      rjl=rjl/a
      p=p+1
      if (p .le. max) go to 65
      call prtbuf(n,max,4, m, v)
      p=1
65    m(1,p)=k
      m(2,p)=i
      m(3,p)=j
      m(4,p)=l
      cosine=pik*pjl+qik*qjl+rik*rjl
      sine=pij*(qjl*rik-rjl*qik)
     *   + qij*(rjl*pik-pjl*rik) + rij*(pjl*qik-qjl*pik)
      v(p)=(180.0d0/pi)* datan2(sine,cosine)
70    continue
80    continue
90    continue
100   continue
      if (p .gt. 0) call prtbuf(n,p,4, m, v)
      if (p .eq. 0) write (iwr,1009)
c
c     ----- reset core memory
c
      call gmem_free_inf(i10,fnm,snm,'i10')
      return
      end
      subroutine struct_check(n, c, az, nuct, bl)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
      dimension c(3,*),az(*),nuct(*)
      dimension bl(*)
      dimension ilifz(maxat)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
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
c approximate covalent radius table in au
      real*8 cov
      common/coval/cov(100)
c
      real*8 scalel, vbond1, vbond2
      logical opoint1,ochkdst,ofatal
      integer mxbnd, nbond1, nbond2, mbond1, mbond2
      parameter (mxbnd = 1000)
c
      common/structchk/ vbond1(mxbnd),mbond1(2,mxbnd),
     +                  vbond2(mxbnd),mbond2(2,mxbnd),
     +                  scalel, nbond1, nbond2, opoint1,
     +                  ochkdst, ofatal
c
      real*8 toang, ams
      integer ifau, nosymm, iseczz
      common/phycon/toang(30),ams(54),ifau,nosymm,iseczz
      character*7 fnm
      character*12 snm
      data fnm,snm/"util1.m","struct_check"/
c
      data two/2.0d0/
      data zsadd/'saddle'/
c
      if(zruntp.eq.zsadd) return
c
c     ----- get core memory
c
      length=n*n
      i10 = igmem_alloc_inf(length,fnm,snm,'i10',IGMEM_DEBUG)
      scale = 1.15d0
      toler = 0.5d0
      obond = .false.
      odiss = .false.
c
      do i=1,n
       ilifz(i)=(i-1)*n+i10-1
      enddo
c
      call vclr(bl(i10),1,length)
      nm1=n-1
      do i=1,nm1
      iz= nint(az(i))
      if(nuct(i).le.3)then
       covi=two
       if (iz .gt. 0 .and. iz .le. 100) covi=cov(iz)
       ip1=i+1
c
       do j=ip1,n
        if(nuct(j).le.3) then
         ij=i+ilifz(j)
         bl(ij)=
     +   dsqrt((c(1,j)-c(1,i))**2+(c(2,j)-c(2,i))**2+(c(3,j)-c(3,i))**2)
         jz= nint(az(j))
         covj=two
         if (jz .gt. 0 .and. jz .le. 100) covj=cov(jz)
         if (bl(ij) .gt. (covi+covj)*scale+toler) bl(ij)=-bl(ij)
         bl(j+ilifz(i))=bl(ij)
c
        endif
       enddo
c
      endif
      enddo
c
c
      if (opoint1) then
c
      nbond1=0
      do i=1,nm1
      ip1=i+1
       do j=ip1,n
       ij=i+ilifz(j)
       if (bl(ij) .gt. 0.0d0) then
        nbond1=nbond1+1
        if (nbond1.gt.mxbnd) then
         write(iwr,1002) nbond1
         go to 200
        endif
        mbond1(1,nbond1)=i
        mbond1(2,nbond1)=j
        vbond1(nbond1)=bl(ij) *toang(1)
       endif
       enddo
      enddo
c
      if(oprint(56)) then
       write (iwr,1001)
       write (iwr,1005)
       do i=1,nbond1
        write(iwr,1003) mbond1(1,i),mbond1(2,i),vbond1(i)
       enddo
       write(iwr,1004)
      endif

      else
c
c     all subsequent points
c
      if(oprint(56)) write (iwr,1001)
c
      nbond2=0
      do i=1,nm1
      ip1=i+1
       do j=ip1,n
       ij=i+ilifz(j)
       if (bl(ij) .gt. 0.0d0) then
        nbond2=nbond2+1
        if (nbond2 .gt. mxbnd) go to 250
         mbond2(1,nbond2)=i
         mbond2(2,nbond2)=j
         vbond2(nbond2)=bl(ij) *toang(1)
       endif
       enddo
      enddo
c
      if(oprint(56)) then
       do i=1,nbond2
        write(iwr,1003) mbond1(1,i),mbond1(2,i),vbond2(i)
       enddo
       write(iwr,1004)
      endif
c
c     now compare with the initial point
c
 250   continue
c
      mprint = 0
      do i=1,nbond1
       mb1 = mbond1(1,i)
       mb2 = mbond1(2,i)
       bl1 = vbond1(i)
c      note that the following check may not function
c      correctly given that nbond can fall to zero ..
c      use change in no. of bonds as a second test
       do j=1,nbond2
        if(mb1.eq.mbond2(1,j).and.mb2.eq.mbond2(2,j)) then
c       examine change in bond length
         if(vbond2(j).gt.scalel*bl1) then
c       trouble ?
         obond = .true.
         mprint = mprint + 1
         if(mprint.eq.1)write(iwr,1240)
c        limit printing to the 1st 20 occurrences
         if(mprint.le.20) then
          write(iwr,1250) mbond2(1,j), mbond2(2,j), vbond2(j)
         endif
         endif
        endif
       enddo
      enddo
c
      endif
c
 200  continue
c
      if (opoint1) then
       opoint1 = .false.
      else
       odiss = nbond1 .ne. nbond2. or. obond
       if (odiss) then
        write (iwr,1009) nbond1, nbond2
       else
       if (.not.obond) write (iwr,1010)
       endif
      endif
c
c     ----- reset core memory
c
      call gmem_free_inf(i10,fnm,snm,'i10')
c
c     should we terminate the optimization ?
c     two consecutive failures do the trick here
c
      ofatal = odiss.and.ofatal
      if (ofatal) then
       write(iwr,1006)
 1006  format(//
     +5x,'***************************************************'/
     +5x,'*** An analysis of two consective points in the ***'/
     +5x,'*** optimisation suggests the molecule is       ***'/
     +5x,'*** dissociating. If this appears unreasonable  ***'/
     +5x,'*** consider redefining the starting geometry   ***'/
     +5x,'*** and using a better starting Hessian.        ***'/
     +5x,'***************************************************'//)
       call caserr2('Molecule appears to be dissociating')
      else if(odiss) then
       ofatal = odiss
      endif
c
      return
c
 1002 format(/10x,' WARNING - maximum of ',i3,
     + ' bonds exceeded - distance checking impeded'/)
 1001 format(10x,60('=')/
     +  10x,'*** distance checking of bond lengths in operation'/
     +  10x,60('='))
 1005 format(10x,'*** First point'/10x,60('='))
 1003 format(21x,i5,2x,i5,3x,f10.5)
 1004 format(10x,60('='))
 1240 format(/10x,68('+'))
 1250 format(10x,'WARNING: bond may be dissociating! ',
     + ' atoms ',I4,2x,I4,1x,F10.5, 'angs.')
 1009 format(//5x,'*** apparent change in bonding detected '/ 
     +5x,'*** no. of bonds detected in starting geometry = ',i3/
     +5x,'*** no. of bonds detected in current geometry  = ',i3/)
 1010 format(/10x,'*** bond length tests validated ***')
      end
      subroutine prtbuf(nat,n,k, m, v)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension v(*),m(4,*)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      real*8 toang, ams
      integer ifau, nosymm, iseczz
      common/phycon/toang(30),ams(54),ifau,nosymm,iseczz
      data xbra,xket/'(',')'/
c
      if(nat.le.99) then
      do 10 i=1,n
      do 10 j=2,k
10    if (m(j,i) .lt. 10) m(j,i)=-m(j,i)
c
      go to (20,20,30,40), k
c
20    write (iwr,1000) (m(1,i), m(2,i), v(i), i=1,n)
1000  format (//1x, 7(i2, '-', i2, f11.7, 2x))
      do 25 i=1,n
25    v(i)=v(i)*toang(1)
      write (iwr,1005) (xbra, v(i), xket, i=1,n)
1005  format (7(7x, a1, f9.7, a1))
      return
c
30    write (iwr,1010) (m(1,i), m(2,i), m(3,i), v(i), i=1,n)
1010  format(/1x,5(i2,2('-',i2),f12.6,2x),i2,2('-',i2),f12.6)
      return
c
40    write (iwr,1020) (m(1,i), m(2,i), m(3,i), m(4,i), v(i), i=1,n)
1020  format (/1x, 5(i2, 3('-', i2), f12.6, 2x))
c
      return
c
      elseif (nat.le.999) then
c
      do 15 i=1,n
      do 15 j=2,k
15    if (m(j,i) .lt. 100) m(j,i)=-m(j,i)
c
      go to (50,50,60,70), k
c
50    write (iwr,1001) (m(1,i), m(2,i), v(i), i=1,n)
1001    format (//1x, 7(i3, '-', i3, f9.5, 2x))
      do 55 i=1,n
55    v(i)=v(i)*toang(1)
      write (iwr,1006) (xbra, v(i), xket, i=1,n)
1006  format (7(9x, a1, f7.5, a1))
      return
c
60    write (iwr,1011) (m(1,i), m(2,i), m(3,i), v(i), i=1,n)
1011  format(/1x,5(i3,2('-',i3),f9.3,2x),i3,2('-',i3),f9.3)
      return
c
70    write (iwr,1021) (m(1,i), m(2,i), m(3,i), m(4,i), v(i), i=1,n)
1021  format (/1x, 5(i3, 3('-', i3), f8.2, 2x))
      return
      else
*     4 digit atom numbers
      do 115 i=1,n
      do 115 j=2,k
115   if (m(j,i) .lt. 1000) m(j,i)=-m(j,i)
c
      go to (150,150,160,170), k
c
150    write (iwr,1101) (m(1,i), m(2,i), v(i), i=1,n)
1101    format (//1x, 7(i4, '-', i4, f9.5, 2x))
      do 155 i=1,n
155    v(i)=v(i)*toang(1)
      write (iwr,1106) (xbra, v(i), xket, i=1,n)
1106  format (7(9x, a1, f7.5, a1))
      return
c
160    write (iwr,1111) (m(1,i), m(2,i), m(3,i), v(i), i=1,n)
1111  format(/1x,5(i4,2('-',i4),f9.3,2x),i4,2('-',i4),f9.3)
      return
c
170    write (iwr,1121) (m(1,i), m(2,i), m(3,i), m(4,i), v(i), i=1,n)
1121  format (/1x, 5(i4, 3('-', i4), f8.2, 2x))
      return
      endif
      end
      subroutine inlist(in,max,ival,ixc)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c     mxtoken - the maximum number of tokens on an input line
c               for space separated tokens: 132/2 plus 1 to allow
c               for look ahead (132 is the line-length defined in
c               common/workc)
c
      integer mxtoken
      parameter(mxtoken=67)
      integer jrec, jump, istrt, inumb, iwidth, nend, nstart
      integer nline, noline, jwidth, nerr
      logical oswit, oterm, oflush
      common/work/jrec,jump,istrt(mxtoken),inumb(mxtoken),iwidth,
     + nend(mxtoken),nstart(mxtoken), nline,noline,jwidth,nerr,
     + oswit,oterm,oflush
c
      dimension in(*)
      data zspace,zto/' ','to'/
c     i=0
 1    call inpa(ztest)
      if(ztest.eq.zspace)go to 2
      jrec=jrec-1
      call inpi(j)
      call inpa(ztest)
      if(ztest.ne.zto)go to 3
      call inpi(k)
      go to 4
 3    k=j
      jrec=jrec-1
 4    if(j.lt.1.or.k.lt.j.or.k.gt.max)call caserr2(
     *'invalid iteration count specified')
      do 5 l=j,k
      if(in(l).ne.ixc)in(l)=ival
 5    continue
      go to 1
 2    return
      end
      subroutine ngdiag(a,q,e,iky,n,ndim,jtype,threshs)
c      version of ligen using f02abf
c
c     note that this version always produces eigenvalues
c    is ascending order
c
      implicit real*8  (a-h,o-z)
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
c...   common for overruling local criteria (e.g. dependency, diag thresholds)
      integer  i_depen_check, i_global_diag, n_depen_check
      logical  o_depen_print
      common/overrule/ i_depen_check,i_global_diag, o_depen_print, 
     1                 n_depen_check
c
c     i_depen_check : criterion for throwing out vectors
c     n_depen_check : number  of vectors to be thrown out 
c     o_depen_print : if true +> extra print 

      dimension a(*),e(*),q(ndim,n),iky(*)
      dimension y(maxorb)
      logical jaco
      common/jacobc/jaco
      thresh = threshs
      if (i_global_diag.ne.-999) thresh = 10.0d0**(-i_global_diag*1.0d0)
c  
      if (jaco) then
         call jacob2(a,q,e,iky,n,ndim,jtype,thresh)
      else
         do 30 i = 1 , n
            do 20 j = 1 , i
               ij = iky(i) + j
               q(i,j) = a(ij)
               q(j,i) = a(ij)
 20         continue
 30      continue
         ifail = 0
         call f02abf(q,ndim,n,e,q,ndim,y,ifail)
      end if
      return
      end
      subroutine jacob2(a,q,e,iky,n,ndim,jtype,thresh)
c     original version of ligen (essentially jacobi routine)
      implicit real*8  (a-h,o-z)
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
      dimension a(*),e(*),q(ndim,n),iky(*)
      dimension y(maxorb)
      data zero,one,half/0.0d0,1.0d0,0.5d0/
      data tenth/0.1d0/
      data m1/1/
      do 30 i = 1 , n
         do 20 j = 1 , n
            q(i,j) = zero
 20      continue
         q(i,i) = one
 30   continue
      if (n.eq.m1) then
         e(1) = a(1)
         return
      else
 40      te = zero
         do 60 i = 2 , n
            i1 = i - 1
            ikyi = iky(i)
            do 50 j = 1 , i1
               temp = dabs(a(j+ikyi))
               if (te.lt.temp) te = temp
 50         continue
 60      continue
         if (te.lt.thresh) then
            do 70 i = 1 , n
               e(i) = a(iky(i)+i)
 70         continue
            go to (140,170,200) , jtype
         else
            te = te*tenth
            do 130 i = 2 , n
               i1 = i - 1
               ip1 = i + m1
               ikyi = iky(i)
               itest = n - ip1
               ii = i + ikyi
               do 80 ir = 1 , n
                  y(ir) = q(ir,i)
 80            continue
               do 110 j = 1 , i1
                  ij = j + ikyi
                  vij = a(ij)
                  if (dabs(vij).ge.te) then
                     vii = a(ii)*half
                     j1 = j - 1
                     jp1 = j + m1
                     ikyj = iky(j)
                     jj = j + ikyj
                     vjj = a(jj)*half
                     temp = vii - vjj
                     tem = dsqrt(temp*temp+vij*vij)
                     if (temp.lt.0) then
                        tem = -tem
                     end if
                     cost = (temp+tem)/vij
                     sint = dsqrt(one/(one+cost*cost))
                     cost = cost*sint
                     temp = vii + vjj
                     a(ii) = temp + tem
                     a(jj) = temp - tem
                     a(ij) = zero
                     if (j1.gt.0) then
                        call drot(j1,a(ikyi+1),1,a(ikyj+1),1,cost,sint)
                     end if
                     if (i1.ge.jp1) then
                        do 90 k = jp1 , i1
                           jj = iky(k) + j
                           vij = a(k+ikyi)
                           a(k+ikyi) = vij*cost + a(jj)*sint
                           a(jj) = a(jj)*cost - vij*sint
 90                     continue
                     end if
                     if (itest.ge.0) then
                        do 100 k = ip1 , n
                           ij = iky(k) + i
                           jj = j + iky(k)
                           vij = a(ij)
                           a(ij) = vij*cost + a(jj)*sint
                           a(jj) = a(jj)*cost - vij*sint
 100                    continue
                     end if
                     call drot(n,y,1,q(1,j),1,cost,sint)
                  end if
 110           continue
               do 120 ir = 1 , n
                  q(ir,i) = y(ir)
 120           continue
 130        continue
            go to 40
         end if
      end if
c
c     sort eigenvalues into increasing order
 140  do 160 min = 1 , n
         jm = min
         em = e(min)
         do 150 j = min , n
            if (e(j).lt.em) then
               em = e(j)
               jm = j
            end if
 150     continue
         if (jm.ne.min) then
            temp = e(jm)
            e(jm) = e(min)
            e(min) = temp
            call dswap(n,q(1,jm),1,q(1,min),1)
         end if
 160  continue
      go to 200
c
c     sort into decreasing order
 170  do 190 max = 1 , n
         jm = max
         em = e(max)
         do 180 j = max , n
            if (e(j).gt.em) then
               jm = j
               em = e(j)
            end if
 180     continue
         if (jm.ne.max) then
            temp = e(jm)
            e(jm) = e(max)
            e(max) = temp
            call dswap(n,q(1,jm),1,q(1,max),1)
         end if
 190  continue
c
 200  return
      end
      subroutine gldiag(l0,ldum,l1,h,ib,eig,vector,ia,iop2)
c
c     ----- general calling  routine for giveis or ligen
c           ligen version -----
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c...   common for overruling local criteria (e.g. dependency, diag thresholds)
      integer  i_depen_check, i_global_diag, n_depen_check
      logical  o_depen_print
      common/overrule/ i_depen_check,i_global_diag, o_depen_print, 
     1                 n_depen_check
c
c     i_depen_check : criterion for throwing out vectors
c     n_depen_check : number  of vectors to be thrown out 
c     o_depen_print : if true +> extra print 

      dimension vector(*),h(*),ib(*),eig(*),ia(*)
      do 10 i=1,l0
 10   ib(i)=(i-1)*ldum
      th=1.0d-15
      if (i_global_diag.ne.-999) th = 10.0d0**(-i_global_diag*1.0d0)
c  

      iop1=2
      call jacobi(h,ia,l0,vector,ib,l1,eig,iop1,iop2,th)
      return
      end
c
c  jacobi - general entry point (includes distributed data diagonaliser 
c           if available) 

      subroutine jacobi(a,iky,newbas,q,ilifq,nrow,e,iop1,
     *iop2,threshs)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)

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
      integer idpdiag, ipdiagmode, ntchnk, limchnk, iparapr,
     & idpdiis, idpmult2, idporth, ipiomode, iptest, ipinvmode,
     & ipdiagif
c
c     ipdiagmode - which diag (IDIAG_PEIGS/_PDSYEV/_PDSYEVD....)
c
c     ipdiagif   - which interface (IDIAG_NO_GAWRAP/IDIAG_GAWRAP)
c
c     idpdiag = dimension for parallel diag
c
c     ipinvmode : what matrix inversion algorithm to use
c                 either the old diag based approach or the more
c                 appropriate ScaLAPACK Cholesky factorisation can 
c                 be used (INV_CHOLESKY or INV_DIAG )
c
      logical odebugp

      common /parcntl/ idpdiag, ipdiagmode, ntchnk, limchnk, iparapr,
     & idpdiis, idpmult2, idporth, ipiomode, iptest, ipinvmode,
     & ipdiagif, odebugp
c
      logical ga_initted
      common/ gainit/ ga_initted

      integer IO_NZ,IO_NZ_S,IO_A
      parameter(IO_NZ=1,IO_NZ_S=2,IO_A=3)

      integer INV_CHOLESKY, INV_DIAG
      parameter(INV_CHOLESKY = 100)
      parameter(INV_DIAG     = INV_CHOLESKY + 1)
c
      integer IDIAG_PEIGS,   IDIAG_PDSYEV, IDIAG_PDSYEVX
      integer IDIAG_PDSYEVD, IDIAG_PDSYEVR
      parameter(IDIAG_PEIGS   = 10)
      parameter(IDIAG_PDSYEV  = IDIAG_PEIGS   + 1)
      parameter(IDIAG_PDSYEVX = IDIAG_PDSYEV  + 1)
      parameter(IDIAG_PDSYEVD = IDIAG_PDSYEVX + 1)
      parameter(IDIAG_PDSYEVR = IDIAG_PDSYEVD + 1)

      integer IDIAG_NO_GAWRAP, IDIAG_GAWRAP
      parameter(IDIAG_NO_GAWRAP=200)
      parameter(IDIAG_GAWRAP=201)

      integer maxlablen
      parameter (maxlablen=20)
      character*20 lab
c
c  change the label definitions (init_time_periodsm in util3 )
c  also
c
      integer TP_ENTIRE
      parameter(TP_ENTIRE=1)
c
c  scf 
c
      integer TP_SCF
      parameter (TP_SCF=TP_ENTIRE+1)
      integer TP_ORFOG
      parameter (TP_ORFOG=TP_SCF+1)
      integer TP_DHSTAR
      parameter (TP_DHSTAR=TP_ORFOG+1)
      integer TP_RDMAT
      parameter (TP_RDMAT=TP_DHSTAR+1)
      integer TP_DIIS
      parameter (TP_DIIS=TP_RDMAT+1)
      integer TP_DIAG
      parameter (TP_DIAG=TP_DIIS+1)
c
c  mp2 code
c
      integer TP_APRDMP2
      parameter (TP_APRDMP2=TP_DIAG+1)

      integer TP_DHSTAR_GOP
      parameter (TP_DHSTAR_GOP=TP_APRDMP2+1)

c spare
      integer TP_APRM1234
      parameter (TP_APRM1234=TP_DHSTAR_GOP+1)
c spare
      integer TP_APRQ34
      parameter (TP_APRQ34=TP_APRM1234+1)

      integer TP_APRQ1
      parameter (TP_APRQ1=TP_APRQ34+1)
      integer TP_APRQ2
      parameter (TP_APRQ2=TP_APRQ1+1)
      integer TP_APRQ2D
      parameter (TP_APRQ2D=TP_APRQ2+1)
      integer TP_APRQ34D
      parameter (TP_APRQ34D=TP_APRQ2D+1)
      integer TP_APRMP2E
      parameter (TP_APRMP2E=TP_APRQ34D+1)

      integer TP_MP1PDM
      parameter (TP_MP1PDM=TP_APRMP2E+1)
      integer TP_MP1PDM_1
      parameter (TP_MP1PDM_1=TP_MP1PDM+1)
      integer TP_MP1PDM_2
      parameter (TP_MP1PDM_2=TP_MP1PDM_1+1)

      integer TP_APR1PDM
      parameter (TP_APR1PDM=TP_MP1PDM_2+1)

      integer TP_MP2HESS
      parameter (TP_MP2HESS=TP_APR1PDM+1)
      integer TP_MP2CHF
      parameter (TP_MP2CHF=TP_MP2HESS+1)
      integer TP_MP2MAKEW
      parameter (TP_MP2MAKEW=TP_MP2CHF+1)
      integer TP_MP2DS
      parameter (TP_MP2DS=TP_MP2MAKEW+1)
      integer TP_MP2BACK_P
      parameter (TP_MP2BACK_P=TP_MP2DS+1)
      integer TP_MP2BACK_F
      parameter (TP_MP2BACK_F=TP_MP2BACK_P+1)
      integer TP_MP2BACKTRAN_2
      parameter (TP_MP2BACKTRAN_2=TP_MP2BACK_F+1)
      integer TP_MP2MCDAB
      parameter (TP_MP2MCDAB=TP_MP2BACKTRAN_2+1)
c
c - parallel functions
c
      integer TP_DGOP
      parameter(TP_DGOP=TP_MP2MCDAB+1)

      integer TP_BCAST
      parameter(TP_BCAST=TP_DGOP+1)

      integer TP_NXTVAL
      parameter(TP_NXTVAL=TP_BCAST+1)

      integer TP_GENRAL
      parameter (TP_GENRAL=TP_NXTVAL+1)

      integer TP_GENRAL_1PDM
      parameter (TP_GENRAL_1PDM=TP_GENRAL+1)

      integer TP_GA_PUT_Q2D
      parameter (TP_GA_PUT_Q2D=TP_GENRAL_1PDM+1)

      integer TP_GA_ACC_Q2D
      parameter (TP_GA_ACC_Q2D=TP_GA_PUT_Q2D+1)

      integer TP_MKT2AO
      parameter (TP_MKT2AO   =TP_GA_ACC_Q2D+1)

      integer TP_APRL2
      parameter (TP_APRL2    =TP_MKT2AO+1)

      integer TP_GA_GET_L2
      parameter (TP_GA_GET_L2=TP_APRL2+1)

      integer TP_APRL34
      parameter (TP_APRL34   =TP_GA_GET_L2+1)

      integer TP_DGENRL
      parameter (TP_DGENRL   =TP_APRL34+1)

      integer TP_MP2
      parameter(TP_MP2       =TP_DGENRL+1)

      integer TP_JKDER
      parameter(TP_JKDER     =TP_MP2+1)

      integer TP_APRL1234
      parameter (TP_APRL1234 =TP_JKDER+1)

      integer TP_APRL1
      parameter (TP_APRL1    =TP_APRL1234+1)

      integer TP_APRDMP2_I
      parameter(TP_APRDMP2_I =TP_APRL1+1)

      integer TP_JKDER_GET
      parameter(TP_JKDER_GET =TP_APRDMP2_I+1)

      integer TP_MP1PDM_3
      parameter(TP_MP1PDM_3  =TP_JKDER_GET+1)

      integer TP_MP1PDM_4
      parameter(TP_MP1PDM_4  =TP_MP1PDM_3+1)

cgdf:  added TP_MKT2MO 14.3.95
      integer TP_MKT2MO
      parameter (TP_MKT2MO   =TP_MP1PDM_4+1)

csecd - time periods for second derivatives
c
      integer TP_2D_AOINTS
      parameter(TP_2D_AOINTS     =TP_MKT2MO+1)

      integer TP_2D_SCF
      parameter(TP_2D_SCF        =TP_2D_AOINTS+1)

      integer TP_2D_HFGRDN
      parameter(TP_2D_HFGRDN     =TP_2D_SCF+1)

      integer TP_2D_INDX2T
      parameter(TP_2D_INDX2T     =TP_2D_HFGRDN+1)

      integer TP_2D_MOINTS
      parameter(TP_2D_MOINTS     =TP_2D_INDX2T+1)

      integer TP_2D_TRNFKD
      parameter(TP_2D_TRNFKD     =TP_2D_MOINTS+1)

      integer TP_2D_CHFNDR
      parameter(TP_2D_CHFNDR     =TP_2D_TRNFKD+1)

      integer TP_2D_QMDER
      parameter(TP_2D_QMDER      =TP_2D_CHFNDR+1)

      integer TP_2D_DMDER
      parameter(TP_2D_DMDER      =TP_2D_QMDER+1)

      integer TP_2D_2D
      parameter(TP_2D_2D         =TP_2D_DMDER+1)

      integer TP_2D_CHF
      parameter(TP_2D_CHF        =TP_2D_2D+1)

      integer TP_2D_NUC
      parameter(TP_2D_NUC        =TP_2D_CHF+1)

      integer TP_2D_OVL
      parameter(TP_2D_OVL        =TP_2D_NUC+1)

      integer TP_2D_KE
      parameter(TP_2D_KE         =TP_2D_OVL+1)

      integer TP_2D_PE
      parameter(TP_2D_PE         =TP_2D_KE+1)

      integer TP_2D_2E
      parameter(TP_2D_2E         =TP_2D_PE+1)

      integer TP_2D_TOTAL
      parameter (TP_2D_TOTAL     =TP_2D_2E+1)

      integer TP_2D_CHFDRV
      parameter (TP_2D_CHFDRV    =TP_2D_TOTAL+1)

      integer TP_2D_PDENS
      parameter (TP_2D_PDENS     =TP_2D_CHFDRV+1)

      integer TP_2D_PFOCK
      parameter (TP_2D_PFOCK     =TP_2D_PDENS+1)

      integer TP_2D_CHFHESS
      parameter (TP_2D_CHFHESS   =TP_2D_PFOCK+1)

      integer TP_2D_CHFRHS
      parameter (TP_2D_CHFRHS    =TP_2D_CHFHESS+1)

      integer TP_2D_SYMMRHS
      parameter (TP_2D_SYMMRHS   =TP_2D_CHFRHS+1)

      integer TP_2D_SYMMU
      parameter (TP_2D_SYMMU     =TP_2D_SYMMRHS+1)

      integer TP_2D_PFOCK_OOOO
      parameter (TP_2D_PFOCK_OOOO=TP_2D_SYMMU+1)

      integer TP_2D_PFOCK_VOOO
      parameter (TP_2D_PFOCK_VOOO=TP_2D_PFOCK_OOOO+1)

      integer TP_2D_PFOCK_VVOO
      parameter (TP_2D_PFOCK_VVOO=TP_2D_PFOCK_VOOO+1)

      integer TP_2D_PFOCK_VOVO
      parameter (TP_2D_PFOCK_VOVO=TP_2D_PFOCK_VVOO+1)

      integer TP_2D_PFOCK_SUM
      parameter (TP_2D_PFOCK_SUM =TP_2D_PFOCK_VOVO+1)

      integer TP_2D_AOGEN
      parameter(TP_2D_AOGEN=TP_2D_PFOCK_SUM+1)

      integer TP_2D_AOUT
      parameter(TP_2D_AOUT =TP_2D_AOGEN+1)

      integer TP_TEST1
      parameter(TP_TEST1   =TP_2D_AOUT+1)

      integer TP_TEST2
      parameter(TP_TEST2   =TP_TEST1+1)

      integer TP_TEST3
      parameter(TP_TEST3   =TP_TEST2+1)

      integer TP_TEST4
      parameter(TP_TEST4   =TP_TEST3+1)

      integer TP_TEST5
      parameter(TP_TEST5   =TP_TEST4+1)

      integer TP_TEST6
      parameter(TP_TEST6   =TP_TEST5+1)

      integer TP_HFGRAD
      parameter(TP_HFGRAD  =TP_TEST6+1)

      integer TP_GAMULT2
      parameter(TP_GAMULT2 =TP_HFGRAD+1)

      integer TP_GAORTHOG
      parameter(TP_GAORTHOG=TP_GAMULT2+1)

      integer TP_PDIAG
      parameter (TP_PDIAG  =TP_GAORTHOG+1)

      integer TP_MULT2
      parameter (TP_MULT2  =TP_PDIAG+1)

      integer TP_INTEG
      parameter (TP_INTEG  =TP_MULT2+1)
c
c =================  I/O timers ========================
c
c find
c	
      integer TP_IOFM1, TP_IOF0, TP_IOF1, TP_IOF2, TP_IOF3,
     &	TP_IOF4, TP_IOF5, TP_IOF6, TP_IOF7

      parameter(TP_IOFM1=TP_INTEG+1)
      parameter (TP_IOF0=TP_IOFM1+1)
      parameter (TP_IOF1=TP_IOF0+1)
      parameter (TP_IOF2=TP_IOF1+1)
      parameter (TP_IOF3=TP_IOF2+1)
      parameter (TP_IOF4=TP_IOF3+1)
      parameter (TP_IOF5=TP_IOF4+1)
      parameter (TP_IOF6=TP_IOF5+1)
      parameter (TP_IOF7=TP_IOF6+1)
c
c get
c
      integer TP_IOGM1, TP_IOG0, TP_IOG1, TP_IOG2, TP_IOG3, 
     &  TP_IOG4, TP_IOG5, TP_IOG6, TP_IOG7
      parameter (TP_IOGM1=TP_IOF7+1)
      parameter (TP_IOG0=TP_IOGM1+1)
      parameter (TP_IOG1=TP_IOG0+1)
      parameter (TP_IOG2=TP_IOG1+1)
      parameter (TP_IOG3=TP_IOG2+1)
      parameter (TP_IOG4=TP_IOG3+1)
      parameter (TP_IOG5=TP_IOG4+1)
      parameter (TP_IOG6=TP_IOG5+1)
      parameter (TP_IOG7=TP_IOG6+1)
c
c put
c
      integer TP_IOPM1,  TP_IOP0, TP_IOP1, TP_IOP2, TP_IOP3,
     & TP_IOP4, TP_IOP5, TP_IOP6, TP_IOP7
      parameter (TP_IOPM1=TP_IOG7+1)
      parameter (TP_IOP0=TP_IOPM1+1)
      parameter (TP_IOP1=TP_IOP0+1)
      parameter (TP_IOP2=TP_IOP1+1)
      parameter (TP_IOP3=TP_IOP2+1)
      parameter (TP_IOP4=TP_IOP3+1)
      parameter (TP_IOP5=TP_IOP4+1)
      parameter (TP_IOP6=TP_IOP5+1)
      parameter (TP_IOP7=TP_IOP6+1)
c
c open
c
      integer TP_IOOM1,TP_IOO0,TP_IOO1,TP_IOO2,TP_IOO3,
     & TP_IOO4,TP_IOO5, TP_IOO6, TP_IOO7

      parameter (TP_IOOM1=TP_IOP7+1)
      parameter (TP_IOO0=TP_IOOM1+1)
      parameter (TP_IOO1=TP_IOO0+1)
      parameter (TP_IOO2=TP_IOO1+1)
      parameter (TP_IOO3=TP_IOO2+1)
      parameter (TP_IOO4=TP_IOO3+1)
      parameter (TP_IOO5=TP_IOO4+1)
      parameter (TP_IOO6=TP_IOO5+1)
      parameter (TP_IOO7=TP_IOO6+1)
c
c delfil (only significant for GA-files dumped to disc
c
      integer TP_IO_GAFILE_READ, TP_IO_GAFILE_DUMP
      parameter (TP_IO_GAFILE_READ      =TP_IOO7+1)
      parameter (TP_IO_GAFILE_DUMP      =TP_IO_GAFILE_READ+1)
c
c Peigs parallel diag
c
      integer TP_PEIGS
      parameter (TP_PEIGS               =TP_IO_GAFILE_DUMP+1)
c
c Scalapack parallel diag
c
      integer TP_PDSYEV
      parameter (TP_PDSYEV               =TP_PEIGS+1)
      integer TP_PDSYEVX
      parameter (TP_PDSYEVX              =TP_PDSYEV+1)
      integer TP_PDSYEVD
      parameter (TP_PDSYEVD              =TP_PDSYEVX+1)
      integer TP_PDSYEVR
      parameter (TP_PDSYEVR              =TP_PDSYEVD+1)
c
c timers for CCP1 DFT module
c
      integer    TP_DFT_JFIT
      parameter (TP_DFT_JFIT            =TP_PDSYEVR+1)
      integer    TP_DFT_JFIT_VFORM
      parameter (TP_DFT_JFIT_VFORM      =TP_DFT_JFIT+1)
      integer    TP_DFT_JFIT_TR
      parameter (TP_DFT_JFIT_TR         =TP_DFT_JFIT_VFORM+1)
      integer    TP_DFT_JFIT_NR
      parameter (TP_DFT_JFIT_NR         =TP_DFT_JFIT_TR+1)
      integer    TP_DFT_JFIT_COEF
      parameter (TP_DFT_JFIT_COEF       =TP_DFT_JFIT_NR+1)
      integer    TP_DFT_JFIT_KSMAT
      parameter (TP_DFT_JFIT_KSMAT      =TP_DFT_JFIT_COEF+1)
      integer    TP_DFT_JFIT_ENERGY
      parameter (TP_DFT_JFIT_ENERGY     =TP_DFT_JFIT_KSMAT+1)
      integer    TP_DFT_EXQUAD
      parameter (TP_DFT_EXQUAD          =TP_DFT_JFIT_ENERGY+1)
      integer    TP_DFT_EXQUAD_INTRO
      parameter (TP_DFT_EXQUAD_INTRO    =TP_DFT_EXQUAD+1)
      integer    TP_DFT_EXQUAD_INTEG
      parameter (TP_DFT_EXQUAD_INTEG    =TP_DFT_EXQUAD_INTRO+1)
      integer    TP_DFT_EXQUAD_DGOP
      parameter (TP_DFT_EXQUAD_DGOP     =TP_DFT_EXQUAD_INTEG+1)
      integer    TP_DFT_EXQUADF
      parameter (TP_DFT_EXQUADF         =TP_DFT_EXQUAD_DGOP+1)
      integer    TP_DFT_EXQUADF_INTRO
      parameter (TP_DFT_EXQUADF_INTRO   =TP_DFT_EXQUADF+1)
      integer    TP_DFT_EXQUADF_INTEG
      parameter (TP_DFT_EXQUADF_INTEG   =TP_DFT_EXQUADF_INTRO+1)
      integer    TP_DFT_EXQUADF_DGOP
      parameter (TP_DFT_EXQUADF_DGOP    =TP_DFT_EXQUADF_INTEG+1)
      integer    TP_DFT_EXQUADLHS
      parameter (TP_DFT_EXQUADLHS       =TP_DFT_EXQUADF_DGOP+1)
      integer    TP_DFT_EXQUADLHS_INTRO
      parameter (TP_DFT_EXQUADLHS_INTRO =TP_DFT_EXQUADLHS+1)
      integer    TP_DFT_EXQUADLHS_INTEG
      parameter (TP_DFT_EXQUADLHS_INTEG =TP_DFT_EXQUADLHS_INTRO+1)
      integer    TP_DFT_EXQUADLHS_DGOP
      parameter (TP_DFT_EXQUADLHS_DGOP  =TP_DFT_EXQUADLHS_INTEG+1)
      integer    TP_DFT_EXQUADHES
      parameter (TP_DFT_EXQUADHES       =TP_DFT_EXQUADLHS_DGOP+1)
      integer    TP_DFT_EXQUADHES_INTRO
      parameter (TP_DFT_EXQUADHES_INTRO =TP_DFT_EXQUADHES+1)
      integer    TP_DFT_EXQUADHES_INTEG
      parameter (TP_DFT_EXQUADHES_INTEG =TP_DFT_EXQUADHES_INTRO+1)
      integer    TP_DFT_EXQUADHES_DGOP
      parameter (TP_DFT_EXQUADHES_DGOP  =TP_DFT_EXQUADHES_INTEG+1)
      integer    TP_DFT_EXQUADRHS
      parameter (TP_DFT_EXQUADRHS       =TP_DFT_EXQUADHES_DGOP+1)
      integer    TP_DFT_EXQUADRHS_INTRO
      parameter (TP_DFT_EXQUADRHS_INTRO =TP_DFT_EXQUADRHS+1)
      integer    TP_DFT_EXQUADRHS_INTEG
      parameter (TP_DFT_EXQUADRHS_INTEG =TP_DFT_EXQUADRHS_INTRO+1)
      integer    TP_DFT_EXQUADRHS_DGOP
      parameter (TP_DFT_EXQUADRHS_DGOP  =TP_DFT_EXQUADRHS_INTEG+1)
      integer    TP_DFT_EXQUADDKSX
      parameter (TP_DFT_EXQUADDKSX      =TP_DFT_EXQUADRHS_DGOP+1)
      integer    TP_DFT_EXQUADDKSX_INTRO
      parameter (TP_DFT_EXQUADDKSX_INTRO=TP_DFT_EXQUADDKSX+1)
      integer    TP_DFT_EXQUADDKSX_INTEG
      parameter (TP_DFT_EXQUADDKSX_INTEG=TP_DFT_EXQUADDKSX_INTRO+1)
      integer    TP_DFT_EXQUADDKSX_DGOP
      parameter (TP_DFT_EXQUADDKSX_DGOP =TP_DFT_EXQUADDKSX_INTEG+1)
      integer    TP_DFT_EXQUADDKS
      parameter (TP_DFT_EXQUADDKS       =TP_DFT_EXQUADDKSX_DGOP+1)
      integer    TP_DFT_EXQUADDKS_INTRO
      parameter (TP_DFT_EXQUADDKS_INTRO =TP_DFT_EXQUADDKS+1)
      integer    TP_DFT_EXQUADDKS_INTEG
      parameter (TP_DFT_EXQUADDKS_INTEG =TP_DFT_EXQUADDKS_INTRO+1)
      integer    TP_DFT_EXQUADDKS_DGOP
      parameter (TP_DFT_EXQUADDKS_DGOP  =TP_DFT_EXQUADDKS_INTEG+1)



      integer    TP_DFT_JMULT
      parameter (TP_DFT_JMULT           =TP_DFT_EXQUADDKS_DGOP+1)
      integer    TP_DFT_JMULT_INTRO
      parameter (TP_DFT_JMULT_INTRO     =TP_DFT_JMULT+1)
      integer    TP_DFT_JMULT_SB
      parameter (TP_DFT_JMULT_SB        =TP_DFT_JMULT_INTRO+1)
      integer    TP_DFT_JMULT_FOCK
      parameter (TP_DFT_JMULT_FOCK      =TP_DFT_JMULT_SB+1)
      integer    TP_DFT_JMULT_CJAT0
      parameter (TP_DFT_JMULT_CJAT0     =TP_DFT_JMULT_FOCK+1)
c
c coulomb fitted gradients
c
      integer    TP_DFT_JFITG
      parameter (TP_DFT_JFITG           =TP_DFT_JMULT_CJAT0+1)

      integer    TP_DFT_JFITG_VFORM
      parameter (TP_DFT_JFITG_VFORM     =TP_DFT_JFITG+1)
      integer    TP_DFT_JFITG_TR
      parameter (TP_DFT_JFITG_TR        =TP_DFT_JFITG_VFORM+1)
      integer    TP_DFT_JFITG_NR
      parameter (TP_DFT_JFITG_NR        =TP_DFT_JFITG_TR+1)
      integer    TP_DFT_JFITG_COEF
      parameter (TP_DFT_JFITG_COEF      =TP_DFT_JFITG_NR+1)
      integer    TP_DFT_JFITG_2C
      parameter (TP_DFT_JFITG_2C        =TP_DFT_JFITG_COEF+1)
      integer    TP_DFT_JFITG_3C
      parameter (TP_DFT_JFITG_3C        =TP_DFT_JFITG_2C+1)
      integer    TP_DFT_JFIT_TR_INIT
      parameter (TP_DFT_JFIT_TR_INIT    =TP_DFT_JFITG_3C+1)
      integer    TP_DFT_JFIT_INV
      parameter (TP_DFT_JFIT_INV        =TP_DFT_JFIT_TR_INIT+1)
c
c VB
c
      integer TP_VB
      parameter (TP_VB                  =TP_DFT_JFIT_INV+1)
      integer TP_VB_STRUC
      parameter (TP_VB_STRUC            =TP_VB+1)
      integer TP_VB_ME
      parameter (TP_VB_ME               =TP_VB_STRUC+1)
      integer TP_VB_DIAG
      parameter (TP_VB_DIAG             =TP_VB_ME+1)
      integer TP_VB_TRAN
      parameter (TP_VB_TRAN             =TP_VB_DIAG+1)
      integer TP_VB_VIRT
      parameter (TP_VB_VIRT             =TP_VB_TRAN+1)
      integer TP_VB_LADM
      parameter (TP_VB_LADM             =TP_VB_VIRT+1)
      integer TP_VB_DTRAN
      parameter (TP_VB_DTRAN            =TP_VB_LADM+1)
c
c VB parallel extra
c
      integer TP_ISEND
      parameter (TP_ISEND               =TP_VB_DTRAN+1)
      integer TP_IRECV
      parameter (TP_IRECV               =TP_ISEND+1)
      integer TP_WAIT
      parameter (TP_WAIT                =TP_IRECV+1)
c
c One electron derivatives
c
      integer TP_STVECP
      parameter (TP_STVECP              =TP_WAIT+1)

      integer TP_TVDER
      parameter (TP_TVDER               =TP_STVECP+1)

      integer TP_SDER
      parameter (TP_SDER                =TP_TVDER+1)

      integer TP_SGRAD
      parameter (TP_SGRAD               =TP_SDER+1)

      integer TP_HELFEY
      parameter (TP_HELFEY              =TP_SGRAD+1)


      !F90 time periods start here
      integer TP_F90_START
      parameter( TP_F90_START           = TP_HELFEY+1 )
      integer TP_F90_SCF
      parameter( TP_F90_SCF             = TP_F90_START+1 )
      integer TP_F90_BUILD
      parameter( TP_F90_BUILD           = TP_F90_SCF+1 )
      integer TP_F90_DIIS
      parameter( TP_F90_DIIS            = TP_F90_BUILD+1 )
      integer TP_F90_SIMIL
      parameter( TP_F90_SIMIL           = TP_F90_DIIS+1 )
      integer TP_F90_DIAG
      parameter( TP_F90_DIAG            = TP_F90_SIMIL+1 )
      integer TP_F90_BACK
      parameter( TP_F90_BACK            = TP_F90_DIAG+1 )
      integer TP_F90_ASSIGN
      parameter( TP_F90_ASSIGN          = TP_F90_BACK+1 )
      integer TP_F90_ORTHOG
      parameter( TP_F90_ORTHOG          = TP_F90_ASSIGN+1 )
      integer TP_F90_MAKE_DENS
      parameter( TP_F90_MAKE_DENS       = TP_F90_ORTHOG+1 )
      integer TP_F90_END
      parameter( TP_F90_END             = TP_F90_MAKE_DENS+1 )
      integer TP_F90_LEV_SHIFT
      parameter( TP_F90_LEV_SHIFT       = TP_F90_END+1 )
      integer TP_F90_TESTER_EVAL
      parameter( TP_F90_TESTER_EVAL     = TP_F90_LEV_SHIFT+1 )
      integer TP_F90_DELTA_EVAL
      parameter( TP_F90_DELTA_EVAL      = TP_F90_TESTER_EVAL+1 )
      integer TP_F90_TDOWN
      parameter( TP_F90_TDOWN           = TP_F90_DELTA_EVAL+1 )
      integer TP_F90_RDMAT
      parameter( TP_F90_RDMAT           = TP_F90_TDOWN+1 )
      integer TP_F90_INTS
      parameter( TP_F90_INTS            = TP_F90_RDMAT+1 )
      integer TP_NEWSCF
      parameter( TP_NEWSCF              = TP_F90_INTS+1 )

      integer TP_DENSCF
      parameter( TP_DENSCF              = TP_NEWSCF+1 )
      integer TP_DENSCF_BUILD
      parameter( TP_DENSCF_BUILD        = TP_DENSCF+1 )
      integer TP_DENSCF_RDMAT
      parameter( TP_DENSCF_RDMAT        = TP_DENSCF_BUILD+1 )
      integer TP_DENSCF_INTS
      parameter( TP_DENSCF_INTS         = TP_DENSCF_RDMAT+1 )
      integer TP_DENSCF_DIAG_S
      parameter( TP_DENSCF_DIAG_S       = TP_DENSCF_INTS+1 )
      integer TP_DENSCF_SIMIL
      parameter( TP_DENSCF_SIMIL        = TP_DENSCF_DIAG_S+1 )
      integer TP_DENSCF_DIAG
      parameter( TP_DENSCF_DIAG         = TP_DENSCF_SIMIL+1 )
      integer TP_DENSCF_BACK
      parameter( TP_DENSCF_BACK         = TP_DENSCF_DIAG+1 )
      integer TP_DENSCF_MAKE_DENS
      parameter( TP_DENSCF_MAKE_DENS    = TP_DENSCF_BACK+1 )
      integer TP_DENSCF_TDOWN
      parameter( TP_DENSCF_TDOWN        = TP_DENSCF_MAKE_DENS+1 )

      integer TP_DRHFCL_GA
      parameter( TP_DRHFCL_GA           = TP_DENSCF_TDOWN+1 )
c
c     RPA module
c
      integer TP_RESPONSE
      parameter( TP_RESPONSE            = TP_DRHFCL_GA + 1)
      integer TP_RPA
      parameter( TP_RPA                 = TP_RESPONSE + 1)
      integer TP_TDA
      parameter( TP_TDA                 = TP_RPA + 1)
      integer TP_RPANAL
      parameter( TP_RPANAL              = TP_TDA + 1)
      integer TP_RPA_MO2AO
      parameter( TP_RPA_MO2AO           = TP_RPANAL + 1)
      integer TP_RPA_INT
      parameter( TP_RPA_INT             = TP_RPA_MO2AO + 1)
      integer TP_RPA_CNTRCT
      parameter( TP_RPA_CNTRCT          = TP_RPA_INT + 1)
      integer TP_RPA_AO2MO
      parameter( TP_RPA_AO2MO           = TP_RPA_CNTRCT + 1)
      integer TP_TDA_MO2AO
      parameter( TP_TDA_MO2AO           = TP_RPA_AO2MO + 1)
      integer TP_TDA_INT
      parameter( TP_TDA_INT             = TP_TDA_MO2AO + 1)
      integer TP_TDA_CNTRCT
      parameter( TP_TDA_CNTRCT          = TP_TDA_INT + 1)
      integer TP_TDA_AO2MO
      parameter( TP_TDA_AO2MO           = TP_TDA_CNTRCT + 1)
c
c     Define the common blocks
c
      integer maxtp
      parameter (maxtp=TP_TDA_AO2MO)
      integer maxtpdepth
      parameter (maxtpdepth = 10)
      integer itpdepth
      integer itpstack
      integer ntpc
      integer parent
      real*8  ttotw, ttotc, tsw, tsc
      real*8  ttotu, ttots, tsu, tss
      real*8  taggc
      common/timeperiods/ttotw(maxtp),ttotc(maxtp),
     &     ttotu(maxtp),ttots(maxtp),
     &     tsw(maxtp),tsc(maxtp),
     &     tsu(maxtp),tss(maxtp),
     &     taggc(maxtp),
     &     ntpc(maxtp),parent(maxtp),
     &     itpstack(0:maxtpdepth),itpdepth
      common/timeperiodsc/lab(maxtp)

      dimension a(*),q(*),e(*),iky(*),ilifq(*)
      dimension ppp(2*maxorb),omask(maxorb),iipt(maxorb),
     +             ipt(maxorb)
c...   common for overruling local criteria (e.g. dependency, diag thresholds)
      integer  i_depen_check, i_global_diag, n_depen_check
      logical  o_depen_print
      common/overrule/ i_depen_check,i_global_diag, o_depen_print, 
     1                 n_depen_check
c
c     i_depen_check : criterion for throwing out vectors
c     n_depen_check : number  of vectors to be thrown out 
c     o_depen_print : if true +> extra print 


      integer iskip, nn, ij, ji, k, i, j, ifail
      dimension w(10)
c
      ierr = 0
      thresh = threshs
      if (i_global_diag.ne.-999) thresh = 10.0d0**(-i_global_diag*1.0d0)
c  

      if(newbas .eq. 1)then

         e(1)=a(1)
         q(1)=1.0d0

      else
c
c  use serial jacobi code
c
            if(iop1.ne.1) then
               do  i=1,newbas
                  ilifi=ilifq(i)
                  call vclr(q(ilifi+1),1,nrow)
                  q(ilifi+i)=1.0d0
               enddo
            endif

 86         te=0.0d0
            do i=2,newbas
               i1=i-1
               ikyi=iky(i)
               do j=1,i1
                  temp= dabs(a(j+ikyi))
                  if(te.lt.temp)te=temp
               enddo
            enddo

            if(te.lt.thresh)goto 99
            te=te*0.1d0
      do 22 i=2,newbas
               i1=i-1
               ip1=i+1
               ikyi=iky(i)
               itest=newbas-ip1
               ii=i+ikyi
               ilifi=ilifq(i)
              do 22 j=1,i1
                 ij=j+ikyi
                 vij=a(ij)
                 if( dabs(vij) .lt. te) go to 22
                 vii=a(ii)*0.5d0
                 j1=j-1
                 jp1=j+1
                 ikyj=iky(j)
                 jj=j+ikyj
                 vjj=a(jj)*0.5d0
                 temp=vii-vjj
                 tem=dsqrt(temp*temp+vij*vij)

                 if(temp.lt.0.0d0) tem=-tem
                 cost=(temp+tem)/vij
                 sint=dsqrt(1.0d0/(1.0d0+cost*cost))
                 cost=cost*sint
                 temp=vii+vjj
                 a(ii)=temp+tem
                 a(jj)=temp-tem
                 a(ij)=0.0d0

                 if(j1.gt.0)then
            call drot(j1,a(ikyi+1),1,a(ikyj+1),1,cost,sint)
                endif
                if(i1.ge.jp1)then
                   do k=jp1,i1
                      jj=iky(k)+j
                      vij=a(k+ikyi)
                      a(k+ikyi)=vij*cost+a(jj)*sint
                      a(jj)=a(jj)*cost-vij*sint
                   enddo
                endif
 5              if(itest.ge.0)then
                   do k=ip1,newbas
                      ij=iky(k)+i
                      jj=j+iky(k)
                      vij=a(ij)
                      a(ij)=vij*cost+a(jj)*sint
                      a(jj)=a(jj)*cost-vij*sint
                   enddo
                endif

      call drot(nrow,q(ilifi+1),1,q(ilifq(j)+1),1,cost,sint)

 22            continue

              goto 86
 99           continue
         do 11 i=1,newbas
            omask(i)=.false.
            e(i)=a(iky(i)+i)
            iipt(i)=i/2
 11      continue

         goto (67,55,55,43),iop2
c... binary sort of e.values to increasing value sequence

 55      ipt(1)=1
         do 19 j=2,newbas
             ia=1
             ib=j-1
             test=e(j)
 53          irm1=ib-ia
             if(irm1)58,50,51
 51          ibp=ia+iipt(irm1)
             if(test.lt.e(ipt(ibp)))goto 52
c...  insert into high half
             ia=ibp+1
             goto 53
c... insert into low half
 52          jj=ib
             do 54 i=ibp,ib
                ipt(jj+1)=ipt(jj)
 54          jj=jj-1
             ib=ibp-1
             goto 53
c...  end point of search
 50          jj=ipt(ia)
             if(test.ge.e(jj))goto 57
             ipt(ia+1)=jj
 58          ipt(ia)=j
             goto 19
 57          ipt(ia+1)=j
 19       continue

          goto (67,68,69,43),iop2
c...   sort by decreasing e.value(invert order)
 69       itest=newbas+1
          ip1=iipt(newbas)
          do 41    i=1,ip1
             j=itest-i
             k=ipt(i)
             ipt(i)=ipt(j)
             ipt(j)=k
 41       continue
 68       do 20 i=1,newbas
             k=ipt(i)
             iipt(k)=i
             ppp(i)=e(k)
 20       continue
 59       continue
c
          call dcopy(newbas,ppp(1),1,e(1),1)
c...  iipt(i)=k   means column i is to move to posn k
c...   ipt(i)=k   means column k is to move to posn i
          call sortq(q,ilifq,iipt,newbas,nrow)
          go to 67
c
c ... locking requested
c
 43       do 31 j=1,newbas
             m=ilifq(j)
             temp=0.0d0
             do 32 i=1,newbas
                vij= dabs(q(i+m))
                if(vij.lt.temp.or.omask(i))goto 32
                temp=vij
                k=i
 32          continue
             iipt(j)=k
             omask(k)=.true.
             ppp(k)=e(j)
 31       continue
          goto 59

 67       continue
      endif
      return
      end
c
c utility routine for Peigs interface
c
      subroutine reconstitute_evecs( order, evecs)
      implicit none
c
      integer          order
      real*8 evecs( 1:order * order )
c
      integer first_length, second_length
      integer second_start
c
      first_length  = order * ( order + 1 ) / 2
      second_length = order * ( order - 1 ) / 2
      second_start = first_length + 1
c
      call pg_dgop( 16467, evecs, first_length, '+' )
      call pg_dgop( 16467, evecs( second_start ), second_length,
     &  '+')
c
      end
c
c  Entry point for symmetrised jacobi
c
      subroutine jacobi_symm( a, iky, newbas, q, ilifq, nrow, e, iop1,
     *                       iop2, threshs, 
     +                       basis_symmetry, use_symmetry)
*
*  Version calling distributed data diagonliser and using
* symmetry to ensure the evecs do not become symmetry 
* contaminated.
*     
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
*
      Integer basis_symmetry( 1:nrow )
      Logical use_symmetry
      Logical force_serial
*
*  Parameter for how much to shift diagonal elements in the
* Hamiltonian if using symmetry to help resolve symmetry
* contamination difficulties.
*
      real*8       shift
      Parameter( shift = 1.0d3 )
*
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
      integer idpdiag, ipdiagmode, ntchnk, limchnk, iparapr,
     & idpdiis, idpmult2, idporth, ipiomode, iptest, ipinvmode,
     & ipdiagif
c
c     ipdiagmode - which diag (IDIAG_PEIGS/_PDSYEV/_PDSYEVD....)
c
c     ipdiagif   - which interface (IDIAG_NO_GAWRAP/IDIAG_GAWRAP)
c
c     idpdiag = dimension for parallel diag
c
c     ipinvmode : what matrix inversion algorithm to use
c                 either the old diag based approach or the more
c                 appropriate ScaLAPACK Cholesky factorisation can 
c                 be used (INV_CHOLESKY or INV_DIAG )
c
      logical odebugp

      common /parcntl/ idpdiag, ipdiagmode, ntchnk, limchnk, iparapr,
     & idpdiis, idpmult2, idporth, ipiomode, iptest, ipinvmode,
     & ipdiagif, odebugp
c
      logical ga_initted
      common/ gainit/ ga_initted

      integer IO_NZ,IO_NZ_S,IO_A
      parameter(IO_NZ=1,IO_NZ_S=2,IO_A=3)

      integer INV_CHOLESKY, INV_DIAG
      parameter(INV_CHOLESKY = 100)
      parameter(INV_DIAG     = INV_CHOLESKY + 1)
c
      integer IDIAG_PEIGS,   IDIAG_PDSYEV, IDIAG_PDSYEVX
      integer IDIAG_PDSYEVD, IDIAG_PDSYEVR
      parameter(IDIAG_PEIGS   = 10)
      parameter(IDIAG_PDSYEV  = IDIAG_PEIGS   + 1)
      parameter(IDIAG_PDSYEVX = IDIAG_PDSYEV  + 1)
      parameter(IDIAG_PDSYEVD = IDIAG_PDSYEVX + 1)
      parameter(IDIAG_PDSYEVR = IDIAG_PDSYEVD + 1)

      integer IDIAG_NO_GAWRAP, IDIAG_GAWRAP
      parameter(IDIAG_NO_GAWRAP=200)
      parameter(IDIAG_GAWRAP=201)

      integer maxlablen
      parameter (maxlablen=20)
      character*20 lab
c
c  change the label definitions (init_time_periodsm in util3 )
c  also
c
      integer TP_ENTIRE
      parameter(TP_ENTIRE=1)
c
c  scf 
c
      integer TP_SCF
      parameter (TP_SCF=TP_ENTIRE+1)
      integer TP_ORFOG
      parameter (TP_ORFOG=TP_SCF+1)
      integer TP_DHSTAR
      parameter (TP_DHSTAR=TP_ORFOG+1)
      integer TP_RDMAT
      parameter (TP_RDMAT=TP_DHSTAR+1)
      integer TP_DIIS
      parameter (TP_DIIS=TP_RDMAT+1)
      integer TP_DIAG
      parameter (TP_DIAG=TP_DIIS+1)
c
c  mp2 code
c
      integer TP_APRDMP2
      parameter (TP_APRDMP2=TP_DIAG+1)

      integer TP_DHSTAR_GOP
      parameter (TP_DHSTAR_GOP=TP_APRDMP2+1)

c spare
      integer TP_APRM1234
      parameter (TP_APRM1234=TP_DHSTAR_GOP+1)
c spare
      integer TP_APRQ34
      parameter (TP_APRQ34=TP_APRM1234+1)

      integer TP_APRQ1
      parameter (TP_APRQ1=TP_APRQ34+1)
      integer TP_APRQ2
      parameter (TP_APRQ2=TP_APRQ1+1)
      integer TP_APRQ2D
      parameter (TP_APRQ2D=TP_APRQ2+1)
      integer TP_APRQ34D
      parameter (TP_APRQ34D=TP_APRQ2D+1)
      integer TP_APRMP2E
      parameter (TP_APRMP2E=TP_APRQ34D+1)

      integer TP_MP1PDM
      parameter (TP_MP1PDM=TP_APRMP2E+1)
      integer TP_MP1PDM_1
      parameter (TP_MP1PDM_1=TP_MP1PDM+1)
      integer TP_MP1PDM_2
      parameter (TP_MP1PDM_2=TP_MP1PDM_1+1)

      integer TP_APR1PDM
      parameter (TP_APR1PDM=TP_MP1PDM_2+1)

      integer TP_MP2HESS
      parameter (TP_MP2HESS=TP_APR1PDM+1)
      integer TP_MP2CHF
      parameter (TP_MP2CHF=TP_MP2HESS+1)
      integer TP_MP2MAKEW
      parameter (TP_MP2MAKEW=TP_MP2CHF+1)
      integer TP_MP2DS
      parameter (TP_MP2DS=TP_MP2MAKEW+1)
      integer TP_MP2BACK_P
      parameter (TP_MP2BACK_P=TP_MP2DS+1)
      integer TP_MP2BACK_F
      parameter (TP_MP2BACK_F=TP_MP2BACK_P+1)
      integer TP_MP2BACKTRAN_2
      parameter (TP_MP2BACKTRAN_2=TP_MP2BACK_F+1)
      integer TP_MP2MCDAB
      parameter (TP_MP2MCDAB=TP_MP2BACKTRAN_2+1)
c
c - parallel functions
c
      integer TP_DGOP
      parameter(TP_DGOP=TP_MP2MCDAB+1)

      integer TP_BCAST
      parameter(TP_BCAST=TP_DGOP+1)

      integer TP_NXTVAL
      parameter(TP_NXTVAL=TP_BCAST+1)

      integer TP_GENRAL
      parameter (TP_GENRAL=TP_NXTVAL+1)

      integer TP_GENRAL_1PDM
      parameter (TP_GENRAL_1PDM=TP_GENRAL+1)

      integer TP_GA_PUT_Q2D
      parameter (TP_GA_PUT_Q2D=TP_GENRAL_1PDM+1)

      integer TP_GA_ACC_Q2D
      parameter (TP_GA_ACC_Q2D=TP_GA_PUT_Q2D+1)

      integer TP_MKT2AO
      parameter (TP_MKT2AO   =TP_GA_ACC_Q2D+1)

      integer TP_APRL2
      parameter (TP_APRL2    =TP_MKT2AO+1)

      integer TP_GA_GET_L2
      parameter (TP_GA_GET_L2=TP_APRL2+1)

      integer TP_APRL34
      parameter (TP_APRL34   =TP_GA_GET_L2+1)

      integer TP_DGENRL
      parameter (TP_DGENRL   =TP_APRL34+1)

      integer TP_MP2
      parameter(TP_MP2       =TP_DGENRL+1)

      integer TP_JKDER
      parameter(TP_JKDER     =TP_MP2+1)

      integer TP_APRL1234
      parameter (TP_APRL1234 =TP_JKDER+1)

      integer TP_APRL1
      parameter (TP_APRL1    =TP_APRL1234+1)

      integer TP_APRDMP2_I
      parameter(TP_APRDMP2_I =TP_APRL1+1)

      integer TP_JKDER_GET
      parameter(TP_JKDER_GET =TP_APRDMP2_I+1)

      integer TP_MP1PDM_3
      parameter(TP_MP1PDM_3  =TP_JKDER_GET+1)

      integer TP_MP1PDM_4
      parameter(TP_MP1PDM_4  =TP_MP1PDM_3+1)

cgdf:  added TP_MKT2MO 14.3.95
      integer TP_MKT2MO
      parameter (TP_MKT2MO   =TP_MP1PDM_4+1)

csecd - time periods for second derivatives
c
      integer TP_2D_AOINTS
      parameter(TP_2D_AOINTS     =TP_MKT2MO+1)

      integer TP_2D_SCF
      parameter(TP_2D_SCF        =TP_2D_AOINTS+1)

      integer TP_2D_HFGRDN
      parameter(TP_2D_HFGRDN     =TP_2D_SCF+1)

      integer TP_2D_INDX2T
      parameter(TP_2D_INDX2T     =TP_2D_HFGRDN+1)

      integer TP_2D_MOINTS
      parameter(TP_2D_MOINTS     =TP_2D_INDX2T+1)

      integer TP_2D_TRNFKD
      parameter(TP_2D_TRNFKD     =TP_2D_MOINTS+1)

      integer TP_2D_CHFNDR
      parameter(TP_2D_CHFNDR     =TP_2D_TRNFKD+1)

      integer TP_2D_QMDER
      parameter(TP_2D_QMDER      =TP_2D_CHFNDR+1)

      integer TP_2D_DMDER
      parameter(TP_2D_DMDER      =TP_2D_QMDER+1)

      integer TP_2D_2D
      parameter(TP_2D_2D         =TP_2D_DMDER+1)

      integer TP_2D_CHF
      parameter(TP_2D_CHF        =TP_2D_2D+1)

      integer TP_2D_NUC
      parameter(TP_2D_NUC        =TP_2D_CHF+1)

      integer TP_2D_OVL
      parameter(TP_2D_OVL        =TP_2D_NUC+1)

      integer TP_2D_KE
      parameter(TP_2D_KE         =TP_2D_OVL+1)

      integer TP_2D_PE
      parameter(TP_2D_PE         =TP_2D_KE+1)

      integer TP_2D_2E
      parameter(TP_2D_2E         =TP_2D_PE+1)

      integer TP_2D_TOTAL
      parameter (TP_2D_TOTAL     =TP_2D_2E+1)

      integer TP_2D_CHFDRV
      parameter (TP_2D_CHFDRV    =TP_2D_TOTAL+1)

      integer TP_2D_PDENS
      parameter (TP_2D_PDENS     =TP_2D_CHFDRV+1)

      integer TP_2D_PFOCK
      parameter (TP_2D_PFOCK     =TP_2D_PDENS+1)

      integer TP_2D_CHFHESS
      parameter (TP_2D_CHFHESS   =TP_2D_PFOCK+1)

      integer TP_2D_CHFRHS
      parameter (TP_2D_CHFRHS    =TP_2D_CHFHESS+1)

      integer TP_2D_SYMMRHS
      parameter (TP_2D_SYMMRHS   =TP_2D_CHFRHS+1)

      integer TP_2D_SYMMU
      parameter (TP_2D_SYMMU     =TP_2D_SYMMRHS+1)

      integer TP_2D_PFOCK_OOOO
      parameter (TP_2D_PFOCK_OOOO=TP_2D_SYMMU+1)

      integer TP_2D_PFOCK_VOOO
      parameter (TP_2D_PFOCK_VOOO=TP_2D_PFOCK_OOOO+1)

      integer TP_2D_PFOCK_VVOO
      parameter (TP_2D_PFOCK_VVOO=TP_2D_PFOCK_VOOO+1)

      integer TP_2D_PFOCK_VOVO
      parameter (TP_2D_PFOCK_VOVO=TP_2D_PFOCK_VVOO+1)

      integer TP_2D_PFOCK_SUM
      parameter (TP_2D_PFOCK_SUM =TP_2D_PFOCK_VOVO+1)

      integer TP_2D_AOGEN
      parameter(TP_2D_AOGEN=TP_2D_PFOCK_SUM+1)

      integer TP_2D_AOUT
      parameter(TP_2D_AOUT =TP_2D_AOGEN+1)

      integer TP_TEST1
      parameter(TP_TEST1   =TP_2D_AOUT+1)

      integer TP_TEST2
      parameter(TP_TEST2   =TP_TEST1+1)

      integer TP_TEST3
      parameter(TP_TEST3   =TP_TEST2+1)

      integer TP_TEST4
      parameter(TP_TEST4   =TP_TEST3+1)

      integer TP_TEST5
      parameter(TP_TEST5   =TP_TEST4+1)

      integer TP_TEST6
      parameter(TP_TEST6   =TP_TEST5+1)

      integer TP_HFGRAD
      parameter(TP_HFGRAD  =TP_TEST6+1)

      integer TP_GAMULT2
      parameter(TP_GAMULT2 =TP_HFGRAD+1)

      integer TP_GAORTHOG
      parameter(TP_GAORTHOG=TP_GAMULT2+1)

      integer TP_PDIAG
      parameter (TP_PDIAG  =TP_GAORTHOG+1)

      integer TP_MULT2
      parameter (TP_MULT2  =TP_PDIAG+1)

      integer TP_INTEG
      parameter (TP_INTEG  =TP_MULT2+1)
c
c =================  I/O timers ========================
c
c find
c	
      integer TP_IOFM1, TP_IOF0, TP_IOF1, TP_IOF2, TP_IOF3,
     &	TP_IOF4, TP_IOF5, TP_IOF6, TP_IOF7

      parameter(TP_IOFM1=TP_INTEG+1)
      parameter (TP_IOF0=TP_IOFM1+1)
      parameter (TP_IOF1=TP_IOF0+1)
      parameter (TP_IOF2=TP_IOF1+1)
      parameter (TP_IOF3=TP_IOF2+1)
      parameter (TP_IOF4=TP_IOF3+1)
      parameter (TP_IOF5=TP_IOF4+1)
      parameter (TP_IOF6=TP_IOF5+1)
      parameter (TP_IOF7=TP_IOF6+1)
c
c get
c
      integer TP_IOGM1, TP_IOG0, TP_IOG1, TP_IOG2, TP_IOG3, 
     &  TP_IOG4, TP_IOG5, TP_IOG6, TP_IOG7
      parameter (TP_IOGM1=TP_IOF7+1)
      parameter (TP_IOG0=TP_IOGM1+1)
      parameter (TP_IOG1=TP_IOG0+1)
      parameter (TP_IOG2=TP_IOG1+1)
      parameter (TP_IOG3=TP_IOG2+1)
      parameter (TP_IOG4=TP_IOG3+1)
      parameter (TP_IOG5=TP_IOG4+1)
      parameter (TP_IOG6=TP_IOG5+1)
      parameter (TP_IOG7=TP_IOG6+1)
c
c put
c
      integer TP_IOPM1,  TP_IOP0, TP_IOP1, TP_IOP2, TP_IOP3,
     & TP_IOP4, TP_IOP5, TP_IOP6, TP_IOP7
      parameter (TP_IOPM1=TP_IOG7+1)
      parameter (TP_IOP0=TP_IOPM1+1)
      parameter (TP_IOP1=TP_IOP0+1)
      parameter (TP_IOP2=TP_IOP1+1)
      parameter (TP_IOP3=TP_IOP2+1)
      parameter (TP_IOP4=TP_IOP3+1)
      parameter (TP_IOP5=TP_IOP4+1)
      parameter (TP_IOP6=TP_IOP5+1)
      parameter (TP_IOP7=TP_IOP6+1)
c
c open
c
      integer TP_IOOM1,TP_IOO0,TP_IOO1,TP_IOO2,TP_IOO3,
     & TP_IOO4,TP_IOO5, TP_IOO6, TP_IOO7

      parameter (TP_IOOM1=TP_IOP7+1)
      parameter (TP_IOO0=TP_IOOM1+1)
      parameter (TP_IOO1=TP_IOO0+1)
      parameter (TP_IOO2=TP_IOO1+1)
      parameter (TP_IOO3=TP_IOO2+1)
      parameter (TP_IOO4=TP_IOO3+1)
      parameter (TP_IOO5=TP_IOO4+1)
      parameter (TP_IOO6=TP_IOO5+1)
      parameter (TP_IOO7=TP_IOO6+1)
c
c delfil (only significant for GA-files dumped to disc
c
      integer TP_IO_GAFILE_READ, TP_IO_GAFILE_DUMP
      parameter (TP_IO_GAFILE_READ      =TP_IOO7+1)
      parameter (TP_IO_GAFILE_DUMP      =TP_IO_GAFILE_READ+1)
c
c Peigs parallel diag
c
      integer TP_PEIGS
      parameter (TP_PEIGS               =TP_IO_GAFILE_DUMP+1)
c
c Scalapack parallel diag
c
      integer TP_PDSYEV
      parameter (TP_PDSYEV               =TP_PEIGS+1)
      integer TP_PDSYEVX
      parameter (TP_PDSYEVX              =TP_PDSYEV+1)
      integer TP_PDSYEVD
      parameter (TP_PDSYEVD              =TP_PDSYEVX+1)
      integer TP_PDSYEVR
      parameter (TP_PDSYEVR              =TP_PDSYEVD+1)
c
c timers for CCP1 DFT module
c
      integer    TP_DFT_JFIT
      parameter (TP_DFT_JFIT            =TP_PDSYEVR+1)
      integer    TP_DFT_JFIT_VFORM
      parameter (TP_DFT_JFIT_VFORM      =TP_DFT_JFIT+1)
      integer    TP_DFT_JFIT_TR
      parameter (TP_DFT_JFIT_TR         =TP_DFT_JFIT_VFORM+1)
      integer    TP_DFT_JFIT_NR
      parameter (TP_DFT_JFIT_NR         =TP_DFT_JFIT_TR+1)
      integer    TP_DFT_JFIT_COEF
      parameter (TP_DFT_JFIT_COEF       =TP_DFT_JFIT_NR+1)
      integer    TP_DFT_JFIT_KSMAT
      parameter (TP_DFT_JFIT_KSMAT      =TP_DFT_JFIT_COEF+1)
      integer    TP_DFT_JFIT_ENERGY
      parameter (TP_DFT_JFIT_ENERGY     =TP_DFT_JFIT_KSMAT+1)
      integer    TP_DFT_EXQUAD
      parameter (TP_DFT_EXQUAD          =TP_DFT_JFIT_ENERGY+1)
      integer    TP_DFT_EXQUAD_INTRO
      parameter (TP_DFT_EXQUAD_INTRO    =TP_DFT_EXQUAD+1)
      integer    TP_DFT_EXQUAD_INTEG
      parameter (TP_DFT_EXQUAD_INTEG    =TP_DFT_EXQUAD_INTRO+1)
      integer    TP_DFT_EXQUAD_DGOP
      parameter (TP_DFT_EXQUAD_DGOP     =TP_DFT_EXQUAD_INTEG+1)
      integer    TP_DFT_EXQUADF
      parameter (TP_DFT_EXQUADF         =TP_DFT_EXQUAD_DGOP+1)
      integer    TP_DFT_EXQUADF_INTRO
      parameter (TP_DFT_EXQUADF_INTRO   =TP_DFT_EXQUADF+1)
      integer    TP_DFT_EXQUADF_INTEG
      parameter (TP_DFT_EXQUADF_INTEG   =TP_DFT_EXQUADF_INTRO+1)
      integer    TP_DFT_EXQUADF_DGOP
      parameter (TP_DFT_EXQUADF_DGOP    =TP_DFT_EXQUADF_INTEG+1)
      integer    TP_DFT_EXQUADLHS
      parameter (TP_DFT_EXQUADLHS       =TP_DFT_EXQUADF_DGOP+1)
      integer    TP_DFT_EXQUADLHS_INTRO
      parameter (TP_DFT_EXQUADLHS_INTRO =TP_DFT_EXQUADLHS+1)
      integer    TP_DFT_EXQUADLHS_INTEG
      parameter (TP_DFT_EXQUADLHS_INTEG =TP_DFT_EXQUADLHS_INTRO+1)
      integer    TP_DFT_EXQUADLHS_DGOP
      parameter (TP_DFT_EXQUADLHS_DGOP  =TP_DFT_EXQUADLHS_INTEG+1)
      integer    TP_DFT_EXQUADHES
      parameter (TP_DFT_EXQUADHES       =TP_DFT_EXQUADLHS_DGOP+1)
      integer    TP_DFT_EXQUADHES_INTRO
      parameter (TP_DFT_EXQUADHES_INTRO =TP_DFT_EXQUADHES+1)
      integer    TP_DFT_EXQUADHES_INTEG
      parameter (TP_DFT_EXQUADHES_INTEG =TP_DFT_EXQUADHES_INTRO+1)
      integer    TP_DFT_EXQUADHES_DGOP
      parameter (TP_DFT_EXQUADHES_DGOP  =TP_DFT_EXQUADHES_INTEG+1)
      integer    TP_DFT_EXQUADRHS
      parameter (TP_DFT_EXQUADRHS       =TP_DFT_EXQUADHES_DGOP+1)
      integer    TP_DFT_EXQUADRHS_INTRO
      parameter (TP_DFT_EXQUADRHS_INTRO =TP_DFT_EXQUADRHS+1)
      integer    TP_DFT_EXQUADRHS_INTEG
      parameter (TP_DFT_EXQUADRHS_INTEG =TP_DFT_EXQUADRHS_INTRO+1)
      integer    TP_DFT_EXQUADRHS_DGOP
      parameter (TP_DFT_EXQUADRHS_DGOP  =TP_DFT_EXQUADRHS_INTEG+1)
      integer    TP_DFT_EXQUADDKSX
      parameter (TP_DFT_EXQUADDKSX      =TP_DFT_EXQUADRHS_DGOP+1)
      integer    TP_DFT_EXQUADDKSX_INTRO
      parameter (TP_DFT_EXQUADDKSX_INTRO=TP_DFT_EXQUADDKSX+1)
      integer    TP_DFT_EXQUADDKSX_INTEG
      parameter (TP_DFT_EXQUADDKSX_INTEG=TP_DFT_EXQUADDKSX_INTRO+1)
      integer    TP_DFT_EXQUADDKSX_DGOP
      parameter (TP_DFT_EXQUADDKSX_DGOP =TP_DFT_EXQUADDKSX_INTEG+1)
      integer    TP_DFT_EXQUADDKS
      parameter (TP_DFT_EXQUADDKS       =TP_DFT_EXQUADDKSX_DGOP+1)
      integer    TP_DFT_EXQUADDKS_INTRO
      parameter (TP_DFT_EXQUADDKS_INTRO =TP_DFT_EXQUADDKS+1)
      integer    TP_DFT_EXQUADDKS_INTEG
      parameter (TP_DFT_EXQUADDKS_INTEG =TP_DFT_EXQUADDKS_INTRO+1)
      integer    TP_DFT_EXQUADDKS_DGOP
      parameter (TP_DFT_EXQUADDKS_DGOP  =TP_DFT_EXQUADDKS_INTEG+1)



      integer    TP_DFT_JMULT
      parameter (TP_DFT_JMULT           =TP_DFT_EXQUADDKS_DGOP+1)
      integer    TP_DFT_JMULT_INTRO
      parameter (TP_DFT_JMULT_INTRO     =TP_DFT_JMULT+1)
      integer    TP_DFT_JMULT_SB
      parameter (TP_DFT_JMULT_SB        =TP_DFT_JMULT_INTRO+1)
      integer    TP_DFT_JMULT_FOCK
      parameter (TP_DFT_JMULT_FOCK      =TP_DFT_JMULT_SB+1)
      integer    TP_DFT_JMULT_CJAT0
      parameter (TP_DFT_JMULT_CJAT0     =TP_DFT_JMULT_FOCK+1)
c
c coulomb fitted gradients
c
      integer    TP_DFT_JFITG
      parameter (TP_DFT_JFITG           =TP_DFT_JMULT_CJAT0+1)

      integer    TP_DFT_JFITG_VFORM
      parameter (TP_DFT_JFITG_VFORM     =TP_DFT_JFITG+1)
      integer    TP_DFT_JFITG_TR
      parameter (TP_DFT_JFITG_TR        =TP_DFT_JFITG_VFORM+1)
      integer    TP_DFT_JFITG_NR
      parameter (TP_DFT_JFITG_NR        =TP_DFT_JFITG_TR+1)
      integer    TP_DFT_JFITG_COEF
      parameter (TP_DFT_JFITG_COEF      =TP_DFT_JFITG_NR+1)
      integer    TP_DFT_JFITG_2C
      parameter (TP_DFT_JFITG_2C        =TP_DFT_JFITG_COEF+1)
      integer    TP_DFT_JFITG_3C
      parameter (TP_DFT_JFITG_3C        =TP_DFT_JFITG_2C+1)
      integer    TP_DFT_JFIT_TR_INIT
      parameter (TP_DFT_JFIT_TR_INIT    =TP_DFT_JFITG_3C+1)
      integer    TP_DFT_JFIT_INV
      parameter (TP_DFT_JFIT_INV        =TP_DFT_JFIT_TR_INIT+1)
c
c VB
c
      integer TP_VB
      parameter (TP_VB                  =TP_DFT_JFIT_INV+1)
      integer TP_VB_STRUC
      parameter (TP_VB_STRUC            =TP_VB+1)
      integer TP_VB_ME
      parameter (TP_VB_ME               =TP_VB_STRUC+1)
      integer TP_VB_DIAG
      parameter (TP_VB_DIAG             =TP_VB_ME+1)
      integer TP_VB_TRAN
      parameter (TP_VB_TRAN             =TP_VB_DIAG+1)
      integer TP_VB_VIRT
      parameter (TP_VB_VIRT             =TP_VB_TRAN+1)
      integer TP_VB_LADM
      parameter (TP_VB_LADM             =TP_VB_VIRT+1)
      integer TP_VB_DTRAN
      parameter (TP_VB_DTRAN            =TP_VB_LADM+1)
c
c VB parallel extra
c
      integer TP_ISEND
      parameter (TP_ISEND               =TP_VB_DTRAN+1)
      integer TP_IRECV
      parameter (TP_IRECV               =TP_ISEND+1)
      integer TP_WAIT
      parameter (TP_WAIT                =TP_IRECV+1)
c
c One electron derivatives
c
      integer TP_STVECP
      parameter (TP_STVECP              =TP_WAIT+1)

      integer TP_TVDER
      parameter (TP_TVDER               =TP_STVECP+1)

      integer TP_SDER
      parameter (TP_SDER                =TP_TVDER+1)

      integer TP_SGRAD
      parameter (TP_SGRAD               =TP_SDER+1)

      integer TP_HELFEY
      parameter (TP_HELFEY              =TP_SGRAD+1)


      !F90 time periods start here
      integer TP_F90_START
      parameter( TP_F90_START           = TP_HELFEY+1 )
      integer TP_F90_SCF
      parameter( TP_F90_SCF             = TP_F90_START+1 )
      integer TP_F90_BUILD
      parameter( TP_F90_BUILD           = TP_F90_SCF+1 )
      integer TP_F90_DIIS
      parameter( TP_F90_DIIS            = TP_F90_BUILD+1 )
      integer TP_F90_SIMIL
      parameter( TP_F90_SIMIL           = TP_F90_DIIS+1 )
      integer TP_F90_DIAG
      parameter( TP_F90_DIAG            = TP_F90_SIMIL+1 )
      integer TP_F90_BACK
      parameter( TP_F90_BACK            = TP_F90_DIAG+1 )
      integer TP_F90_ASSIGN
      parameter( TP_F90_ASSIGN          = TP_F90_BACK+1 )
      integer TP_F90_ORTHOG
      parameter( TP_F90_ORTHOG          = TP_F90_ASSIGN+1 )
      integer TP_F90_MAKE_DENS
      parameter( TP_F90_MAKE_DENS       = TP_F90_ORTHOG+1 )
      integer TP_F90_END
      parameter( TP_F90_END             = TP_F90_MAKE_DENS+1 )
      integer TP_F90_LEV_SHIFT
      parameter( TP_F90_LEV_SHIFT       = TP_F90_END+1 )
      integer TP_F90_TESTER_EVAL
      parameter( TP_F90_TESTER_EVAL     = TP_F90_LEV_SHIFT+1 )
      integer TP_F90_DELTA_EVAL
      parameter( TP_F90_DELTA_EVAL      = TP_F90_TESTER_EVAL+1 )
      integer TP_F90_TDOWN
      parameter( TP_F90_TDOWN           = TP_F90_DELTA_EVAL+1 )
      integer TP_F90_RDMAT
      parameter( TP_F90_RDMAT           = TP_F90_TDOWN+1 )
      integer TP_F90_INTS
      parameter( TP_F90_INTS            = TP_F90_RDMAT+1 )
      integer TP_NEWSCF
      parameter( TP_NEWSCF              = TP_F90_INTS+1 )

      integer TP_DENSCF
      parameter( TP_DENSCF              = TP_NEWSCF+1 )
      integer TP_DENSCF_BUILD
      parameter( TP_DENSCF_BUILD        = TP_DENSCF+1 )
      integer TP_DENSCF_RDMAT
      parameter( TP_DENSCF_RDMAT        = TP_DENSCF_BUILD+1 )
      integer TP_DENSCF_INTS
      parameter( TP_DENSCF_INTS         = TP_DENSCF_RDMAT+1 )
      integer TP_DENSCF_DIAG_S
      parameter( TP_DENSCF_DIAG_S       = TP_DENSCF_INTS+1 )
      integer TP_DENSCF_SIMIL
      parameter( TP_DENSCF_SIMIL        = TP_DENSCF_DIAG_S+1 )
      integer TP_DENSCF_DIAG
      parameter( TP_DENSCF_DIAG         = TP_DENSCF_SIMIL+1 )
      integer TP_DENSCF_BACK
      parameter( TP_DENSCF_BACK         = TP_DENSCF_DIAG+1 )
      integer TP_DENSCF_MAKE_DENS
      parameter( TP_DENSCF_MAKE_DENS    = TP_DENSCF_BACK+1 )
      integer TP_DENSCF_TDOWN
      parameter( TP_DENSCF_TDOWN        = TP_DENSCF_MAKE_DENS+1 )

      integer TP_DRHFCL_GA
      parameter( TP_DRHFCL_GA           = TP_DENSCF_TDOWN+1 )
c
c     RPA module
c
      integer TP_RESPONSE
      parameter( TP_RESPONSE            = TP_DRHFCL_GA + 1)
      integer TP_RPA
      parameter( TP_RPA                 = TP_RESPONSE + 1)
      integer TP_TDA
      parameter( TP_TDA                 = TP_RPA + 1)
      integer TP_RPANAL
      parameter( TP_RPANAL              = TP_TDA + 1)
      integer TP_RPA_MO2AO
      parameter( TP_RPA_MO2AO           = TP_RPANAL + 1)
      integer TP_RPA_INT
      parameter( TP_RPA_INT             = TP_RPA_MO2AO + 1)
      integer TP_RPA_CNTRCT
      parameter( TP_RPA_CNTRCT          = TP_RPA_INT + 1)
      integer TP_RPA_AO2MO
      parameter( TP_RPA_AO2MO           = TP_RPA_CNTRCT + 1)
      integer TP_TDA_MO2AO
      parameter( TP_TDA_MO2AO           = TP_RPA_AO2MO + 1)
      integer TP_TDA_INT
      parameter( TP_TDA_INT             = TP_TDA_MO2AO + 1)
      integer TP_TDA_CNTRCT
      parameter( TP_TDA_CNTRCT          = TP_TDA_INT + 1)
      integer TP_TDA_AO2MO
      parameter( TP_TDA_AO2MO           = TP_TDA_CNTRCT + 1)
c
c     Define the common blocks
c
      integer maxtp
      parameter (maxtp=TP_TDA_AO2MO)
      integer maxtpdepth
      parameter (maxtpdepth = 10)
      integer itpdepth
      integer itpstack
      integer ntpc
      integer parent
      real*8  ttotw, ttotc, tsw, tsc
      real*8  ttotu, ttots, tsu, tss
      real*8  taggc
      common/timeperiods/ttotw(maxtp),ttotc(maxtp),
     &     ttotu(maxtp),ttots(maxtp),
     &     tsw(maxtp),tsc(maxtp),
     &     tsu(maxtp),tss(maxtp),
     &     taggc(maxtp),
     &     ntpc(maxtp),parent(maxtp),
     &     itpstack(0:maxtpdepth),itpdepth
      common/timeperiodsc/lab(maxtp)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c...   common for overruling local criteria (e.g. dependency, diag thresholds)
      integer  i_depen_check, i_global_diag, n_depen_check
      logical  o_depen_print
      common/overrule/ i_depen_check,i_global_diag, o_depen_print, 
     1                 n_depen_check
c
c     i_depen_check : criterion for throwing out vectors
c     n_depen_check : number  of vectors to be thrown out 
c     o_depen_print : if true +> extra print 


      dimension a(*),q(*),e(*),iky(*),ilifq(*)

      dimension ppp(2*maxorb),omask(maxorb),iipt(maxorb),
     +     ipt(maxorb)

      integer  k, i, j

      Integer n_irreps, a_point( 1:maxorb ), q_point( 1:maxorb )
      Integer diag_order( 1:2, 1:8 )
      Integer size_diag, block_start

      If( newbas .NE. nrow ) Then
         Call caserr2( 'NEWBAS not equal to NROW in jacobi' )
      End If

      ierr = 0
      thresh = threshs
      if (i_global_diag.ne.-999) thresh = 10.0d0**(-i_global_diag*1.0d0)
c  

      dumtim=dclock()

      If( .Not. use_symmetry ) Then

c        If( opg_root() ) Then
c           Write( 6, * ) 'serial jacobi:', nrow
c        End If

         Call jacobi( a, iky, newbas, q, ilifq, nrow, e, iop1,
     *                iop2, thresh )
      Else

         n_irreps = 0
         Do i = 1, nrow
            n_irreps = Max( n_irreps, basis_symmetry( i ) )
         End Do
      
         Call classify_diags( nrow, basis_symmetry, 
     +                        n_irreps, diag_order )

         Call block_a( nrow, n_irreps, basis_symmetry, diag_order, 
     +                 ipt, a, a_point, q_point, q )

         Do i = 1, nrow * nrow
            q( i ) = 0.0d0
         End Do

         block_start = 1
         Do i = 1, n_irreps

            size_diag = diag_order( 2, i )

            if(size_diag.gt.0) then

               Do j = 1, size_diag
                  ipt ( j ) = a_point( block_start + j - 1 ) -
     +                 a_point( block_start ) 
                  iipt( j ) = ( j - 1 ) * size_diag 
               End Do

               Call jacobi( a( a_point( block_start ) ), ipt,
     +                       size_diag, 
     +                       q( q_point( block_start ) ), iipt,
     +                       size_diag,
     +                       e( block_start ),
     +                       iop1, 1, thresh )

               Call shift_evecs( nrow, size_diag, 
     +                        q, q_point( block_start ) )

               block_start = block_start + size_diag
            endif
         End Do

         Call restore_orig_symmetry( nrow, n_irreps, q, 
     +                               basis_symmetry, diag_order,
     +                               a )
*
c     
c     construct diagonal triangle, A
         
         do i = 1 , nrow
            a(iky(i) + i) = e(i)
            do j = 1 , i-1
               a(iky(i) + j) = 0.0d0
            enddo
         enddo

         do 11 i=1,newbas
            omask(i)=.false.
            e(i)=a(iky(i)+i)
            iipt(i)=i/2
 11      continue

         goto (67,55,55,43),iop2
c... binary sort of e.values to increasing value sequence
 55      ipt(1)=1
         do 19 j=2,newbas
            ia=1
            ib=j-1
            test=e(j)
 53         irm1=ib-ia
            if(irm1)58,50,51
 51         ibp=ia+iipt(irm1)
            if(test.lt.e(ipt(ibp)))goto 52
c...  insert into high half
            ia=ibp+1
            goto 53
c... insert into low half
 52         jj=ib
            do i=ibp,ib
               ipt(jj+1)=ipt(jj)
               jj=jj-1
            enddo
            ib=ibp-1
            goto 53
c...  end point of search
 50         jj=ipt(ia)
            if(test.ge.e(jj))goto 57
            ipt(ia+1)=jj
 58         ipt(ia)=j
            goto 19
 57         ipt(ia+1)=j
 19      continue
         goto (67,68,69,43),iop2
c...   sort by decreasing e.value(invert order)
 69      itest=newbas+1
         ip1=iipt(newbas)
         do i=1,ip1
            j=itest-i
            k=ipt(i)
            ipt(i)=ipt(j)
            ipt(j)=k
         enddo
 68      do i=1,newbas
            k=ipt(i)
            iipt(k)=i
            ppp(i)=e(k)
         enddo
 59      continue
c     
         call dcopy(newbas,ppp(1),1,e(1),1)
c...  iipt(i)=k   means column i is to move to posn k
c...   ipt(i)=k   means column k is to move to posn i
         if (ilifq(2).ne.newbas)call stretch(q,newbas,ilifq(2),
     1          ilifq,0)
         call sortq(q,ilifq,iipt,newbas,nrow)
         go to 67
c     
c     ... locking requested (vectors close packed dim newbas)
c     
 43      do j=1,newbas
            m = (j-1)*newbas
            temp=0.0d0
            do 32 i=1,newbas
               vij= dabs(q(i+m))
               if(vij.lt.temp.or.omask(i))goto 32
               temp=vij
               k=i
 32         continue
            iipt(j)=k
            omask(k)=.true.
            ppp(k)=e(j)
         enddo
         goto 59
 67      continue
      endif
      return
      end

      subroutine sortq(q,ilifq,iipt,newbas,nrow)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      dimension q(*),ilifq(*),iipt(*)
      dimension ppp(2*maxorb),omask(maxorb)
      juse=1
      jnext=maxorb+1
      call setstl(newbas,.false.,omask)
       do 21 i=1,newbas
      if(omask(i))goto 21
      j=i
      call dcopy(nrow,q(ilifq(j)+1),1,ppp(juse),1)
c...   start a permutation cycle
   23 m=iipt(j)
      ilifi=ilifq(m)
      call dcopy(nrow,q(ilifi+1),1,ppp(jnext),1)
      call dcopy(nrow,ppp(juse),1,q(ilifi+1),1)
      if(m.eq.i)goto 21
      juse=jnext
      jnext=maxorb+2-jnext
      omask(m)=.true.
      j=m
      goto 23
   21 continue
      return
      end

      logical function odpdiag(idim)
      implicit none
      integer idim
c
      integer idpdiag, ipdiagmode, ntchnk, limchnk, iparapr,
     & idpdiis, idpmult2, idporth, ipiomode, iptest, ipinvmode,
     & ipdiagif
c
c     ipdiagmode - which diag (IDIAG_PEIGS/_PDSYEV/_PDSYEVD....)
c
c     ipdiagif   - which interface (IDIAG_NO_GAWRAP/IDIAG_GAWRAP)
c
c     idpdiag = dimension for parallel diag
c
c     ipinvmode : what matrix inversion algorithm to use
c                 either the old diag based approach or the more
c                 appropriate ScaLAPACK Cholesky factorisation can 
c                 be used (INV_CHOLESKY or INV_DIAG )
c
      logical odebugp

      common /parcntl/ idpdiag, ipdiagmode, ntchnk, limchnk, iparapr,
     & idpdiis, idpmult2, idporth, ipiomode, iptest, ipinvmode,
     & ipdiagif, odebugp
c
      logical ga_initted
      common/ gainit/ ga_initted

      integer IO_NZ,IO_NZ_S,IO_A
      parameter(IO_NZ=1,IO_NZ_S=2,IO_A=3)

      integer INV_CHOLESKY, INV_DIAG
      parameter(INV_CHOLESKY = 100)
      parameter(INV_DIAG     = INV_CHOLESKY + 1)
c
      integer IDIAG_PEIGS,   IDIAG_PDSYEV, IDIAG_PDSYEVX
      integer IDIAG_PDSYEVD, IDIAG_PDSYEVR
      parameter(IDIAG_PEIGS   = 10)
      parameter(IDIAG_PDSYEV  = IDIAG_PEIGS   + 1)
      parameter(IDIAG_PDSYEVX = IDIAG_PDSYEV  + 1)
      parameter(IDIAG_PDSYEVD = IDIAG_PDSYEVX + 1)
      parameter(IDIAG_PDSYEVR = IDIAG_PDSYEVD + 1)

      integer IDIAG_NO_GAWRAP, IDIAG_GAWRAP
      parameter(IDIAG_NO_GAWRAP=200)
      parameter(IDIAG_GAWRAP=201)

      odpdiag = .false.
      end

***************************************************************
*     
      Subroutine characterize_mo( n, n_vectors, vectors, symm_ao, 
     +                            n_irreps, symm_mo, error )
*     
*  Given a set of vectors VECTORS and the irreducible representations
* that the basis functions for the vectors span this routine will
* try to identify which of thos irreducible representations each
* of the vectors span.
*     
*  Input:
*  N         : The order of the vectors
*  N_VECTORS : The number of vectors
*  VECTORS   : The vectors
*  SYMM_AO   : The irreps of the basis functions
*     
*  Output:
*  N_IRREPS  : The number of different irreps spanned
*  SYMM_MO   : The irreps that the vectors span
*  ERROR     : Non-zero on failure.
*     
      Implicit None
*     
      Integer n
      Integer n_vectors
      real*8    vectors( 1:n, 1:n_vectors )
      Integer symm_ao( 1:n )
      Integer n_irreps
      Integer symm_mo( 1:n_vectors )
      Integer error
*     
*     In any vector a contribution from a vector of another
*     symmetry up to this value is tolerated. The value is chosen to
*     be consistant with that in the subroutine ANALMO.
*     
      real*8       symmetry_tolerance
      Parameter( symmetry_tolerance = 1.0d-3 )
*     
      external idamax
      integer  idamax
*     
      real*8    biggest_wrong_symm
      Integer symm_this_mo
      Integer count_irrep_ao( 1:8 )
      Integer count_irrep_mo( 1:8 )
c     Integer where_wrong
      Integer i, j
*     
*  Find out how many irreps are spanned by the vectors
*
      n_irreps = 0
      Do i = 1, n_vectors
         n_irreps = Max( n_irreps, symm_ao( i ) )
      End Do
*     
      error = 0
*     
      Do i = 1, n_vectors
*     
*     Find the biggest element of this vector and so the symmetry
*     of the corresponding basis function. If there is no symmetry
*     contamination this will be the symmetry of the vector.
*     
         symm_this_mo = symm_ao( idamax( n, vectors( 1, i ), 1 ) )
*     
*     Check for symmetry contamination.
*     
         biggest_wrong_symm = -1.0d0
*     
         Do j = 1, n
            If( symm_this_mo .NE. symm_ao( j ) ) Then
*               biggest_wrong_symm = Max( Abs( vectors( j, i ) ),
*     +                                   biggest_wrong_symm )
*               where_wrong = j
               If( Abs( vectors( j, i ) ) .GT. biggest_wrong_symm )Then
                  biggest_wrong_symm = Abs( vectors( j, i ) )
c                 where_wrong = j
               End If
            End If
         End Do
*     
         If( biggest_wrong_symm .GT. symmetry_tolerance ) Then
            error = 1
            Go To 10
         Else
            symm_mo( i ) = symm_this_mo
         End If
*     
      End Do
*     
 10   Continue
*
*  If we have the same number of vectors as basis functions
* check that the number of times each irrep is spanned
* by the vectors is the number of times it is spanned by
* the basis.
*
      If( error .EQ. 0 .And. n .EQ. n_vectors ) Then
*
         Do i = 1, n_irreps
            count_irrep_ao( i ) = 0
            count_irrep_mo( i ) = 0
         End Do
*
         Do i = 1, n
            count_irrep_ao( symm_ao( i ) ) = 
     +           count_irrep_ao( symm_ao( i ) ) + 1
            count_irrep_mo( symm_mo( i ) ) = 
     +           count_irrep_mo( symm_mo( i ) ) + 1
         End Do
*
         Do i = 1, n_irreps
            If( count_irrep_ao( i ) .NE. count_irrep_mo( i ) ) Then
               error = 2
            End If
         End Do
*
      End If
*
      End
*
***************************************************************
*     
*
***************************************************************
*
      Subroutine classify_diags( n, symmetry, n_irreps, diag_order )
*
*  This routine examines the symmetry label of each MO to
* find out the size of the diagonalization required for
* each irrep. It then orders the diag list ( in diag
* order ) so that the largest diag is done first, the
* second second etc.
*     
      Implicit None
*     
      Integer n
      Integer symmetry( 1:n )
      Integer n_irreps
      Integer diag_order( 1:2, 1:8 )
*     
      Integer count_irrep_ao( 1:8 )
      Integer this_count
      Integer i, j
*
      n_irreps = 0
      Do i = 1, n
         n_irreps = Max( n_irreps, symmetry( i ) )
      End Do
*
      Do i = 1, n_irreps
         count_irrep_ao( i ) = 0
      End do
*
*  Count how many AOs span each irrep and store that info
*
      Do i = 1, n
         count_irrep_ao( symmetry( i ) ) = 
     +        count_irrep_ao( symmetry( i ) ) + 1
      End Do
*
      Do i = 1, n_irreps
         diag_order( 1, i ) = i
         diag_order( 2, i ) = count_irrep_ao( i )
      End Do
*
*  Simple insertion sort in decreasing order of the
* size of the diag.
*
      Do j = 2, n_irreps
         this_count = diag_order( 2, j )
         Do i = j - 1, 1, -1
            If( diag_order( 2, i ) .GE. this_count ) Then
               Goto 10
            End If
            diag_order( 1, i + 1 ) = diag_order( 1, i )
            diag_order( 2, i + 1 ) = diag_order( 2, i )
         End Do
         i = 0
 10      diag_order( 1, i + 1 ) = j
         diag_order( 2, i + 1 ) = this_count
      End Do
*
      End
*
************************************************************************
*
      subroutine qmat_symm(q,s,v,e,scr,ia,l0,l1,l3,ndim,out,
     +                     isymmo, use_symmetry )
      implicit real*8  (a-h,p-w),integer    (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)

      real*8 q(*)
      Integer isymmo( 1:l1 )
      Logical use_symmetry

      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
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
      integer ibl3s, ibl3t, ibl3f, ibl3x, ibl3y, ibl3z
      integer ibl3rs, ibl3qa, ibl3pa, ibl3ea, ibl3qb, ibl3pb
      integer ibl3eb, ibl3g, ibl3hs, ionsec, ibl3op, ions2
      integer isecqm
      common /dump3/ ibl3s, ibl3t, ibl3f, ibl3x, ibl3y, ibl3z,
     +               ibl3rs,ibl3qa,ibl3pa,ibl3ea,ibl3qb,ibl3pb,
     +               ibl3eb,ibl3g,ibl3hs,ionsec,ibl3op,ions2,
     +               isecqm
c
c
      integer ibl7s, ibl7t, ibl7f, ibl7st,ibl7tt,ibl7ft,ibl3qs
      integer ibl7qa,ibl7pa,ibl7ea,ibl7qb,ibl7pb,ibl7eb
      integer ibl7x, ibl7y, ibl7z
      integer ibl7so,ibl7o, ibl7la,ibl7ec,ibl7ha
      common /scra7/ ibl7s, ibl7t, ibl7f, ibl7st,ibl7tt,ibl7ft,ibl3qs,
     +               ibl7qa,ibl7pa,ibl7ea,ibl7qb,ibl7pb,ibl7eb,
     +               ibl7x, ibl7y, ibl7z,
     +               ibl7so,ibl7o(8), ibl7la, ibl7ec, ibl7ha
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
c...   common for harmonic option
      logical  oharm,opharm,odepen
      integer newbas0, newbas1, nsym0, ilifq0, ielimh
      integer newbash,nsymh
      common/harmon/ oharm,opharm,newbas0,newbas1,nsym0(8),
     1               ilifq0(maxorb),ielimh(maxorb),
     2               newbash,nsymh(8),odepen
c
c...   common for overruling local criteria (e.g. dependency, diag thresholds)
      integer  i_depen_check, i_global_diag, n_depen_check
      logical  o_depen_print
      common/overrule/ i_depen_check,i_global_diag, o_depen_print, 
     1                 n_depen_check
c
c     i_depen_check : criterion for throwing out vectors
c     n_depen_check : number  of vectors to be thrown out 
c     o_depen_print : if true +> extra print 

c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
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
      common/craypk/mmmm(65),isymaos(maxorb),isymmos(maxorb)
      common/restri/jjfile(63),lds(508),isect(508),ldsect(508),
     + iacsct(508)
      Integer itab
      dimension ifreq(31)
      dimension s(*),v(ndim,*),e(*),scr(*),ia(*)
      character*7 fnm
      character*9 snm
      data fnm,snm/'util1.m','qmat_symm'/
      data m167,done,tol1,tol0 /167,1.0d0,1.0d-05,1.0d-07/
c
c...  allow specfied crit for dependency (10**-i_depen_check)
c...  (0 => no checking) ; i_depen_check in common/overrule/
c
      if (i_depen_check.ge.0) then
         tol0 = 10.0d0**(-1.0d0*i_depen_check)
         if (i_depen_check.eq.0) tol0 = 1.0d-20
      end if
c
c     ----- diagonalize overlap matrix -----
c
      if (out) then
       write (iwr,9008)
       call prtril(s,l1)
      endif
      len1=lensec(l1+1)
      l2=l1*(l1+1)/2
      lenv=len1+lensec(l2)
      call sectst(isecqm,itest)
      if(itest.eq.0) go to 300
      call secget(isecqm,m167,iblkv)
      call rdedx_prec(scr,1,iblkv,idaf)
      if(scr(1).ne.dble(l1)) go to 300
      call rdedx(scr,l1+1,iblkv,idaf)
      call dcopy(l1,scr(2),1,e,1)
      call reads(scr,l2,idaf)
      call vsub(s,1,scr,1,scr,1,l2)
      do i=1,l2
      if( dabs(scr(i)).gt.tol1) go to 300
      enddo
      call reads(v,l3,idaf)
      ibl3qs=iblkv+lenv
      l0=newbas0
      go to 400
300   lenq=lenv+lensec(l3)
      call secput(isecqm,m167,lenq,iblkv)
      call wrt3(s,l2,iblkv+len1,idaf)
c
      if (oharm) call comharm(s,'tri',isymmo)
      call gldiag_symm(newbas0,newbas0,newbas0,s,scr,e,v,ia,2, 
     +                 isymmo, use_symmetry )
      if (oharm) then
         nn = newbas1-newbas0
         call expharm(v,'vectors',scr)
c...     move real vectors to the end
         if (nn.ne.0) then
            do i=l0,1,-1
               e(i+nn) = e(i)
               call dcopy(ndim,v(1,i),1,v(1,i+nn),1)
            end do
            call vclr(v,1,ndim*nn)
            call vclr(e,1,nn)
         endif
      end if
c
      if (nprint .eq. 5) then
        write (iwr,9028)
        call prev(v,e,l1,l1,ndim)
      endif
c
c     ----- eliminate eigenvectors for which eigenvalue is less
c           than tol=tol0
c
c...     first make a frequency table of eigenvalues
      if (o_depen_print) then
         write(iwr,9052) (i,e(newbas1-newbas0+i),i=1,newbas0)
      end if 
      rb = log10(e(newbas1-newbas0+1))
      if (rb.lt.0.0d0) rb = rb - 1.0d0
      ib = rb
      ie = log10(e(newbas1))
      ib = max(ib,-20)
      ie = min(ie,10)
      do i=1,ie-ib+1
         ifreq(i) = 0
      end do
      do i=1,newbas0
         rrii = log10(e(i+newbas1-newbas0))
         if (rrii.lt.0.0d0) rrii = rrii - 1.0d0
         ii = rrii
         ii = max(ii,-20)
         ii = min(ii,10)
         ifreq(ii-ib+1) = ifreq(ii-ib+1) + 1
      end do
c
      dum = e(l1-newbas0+1)
      j = 0
      k = 0
      kk=0
      do 180 i = 1,l1
      if(e(i) .lt. tol1) then
      kk=kk+1
      endif
      if (e(i) .lt. tol0 .or. 
     1    i.le. (newbas1-newbas0+n_depen_check) ) go to 160
      j = j+1
      e(j) = done/  dsqrt(e(i))
      if (i.ne.j) call dcopy(l1,v(1,i),1,v(1,j),1)
      go to 180
  160 k = k+1
  180 continue
      l0 = l1-k
      if(kk.ne.l1-newbas0) write(iwr,9047) dum,kk-(l1-newbas0),tol1
      if (l0.ne.newbas0)   write (iwr,9048) dum,k-(l1-newbas0),tol0,
     1                                      n_depen_check,l0
c     if (zruntp.ne.'force') n_depen_check = 0
      if (kk.ne.l1-newbas0.or.l0.ne.newbas0) then
         ii = 0
         do i=ib,ie
            ii = max(ii,ifreq(i-ib+1))
         end do
c...     normalise on 80
         ii = ii/80 + 1
         write(iwr,9049) ii
         do i=ib,ie
            if (ifreq(i-ib+1).lt.ii) then
               if (ifreq(i-ib+1).gt.0) write(iwr,9053) i,ifreq(i-ib+1)
               if (ifreq(i-ib+1).eq.0) write(iwr,9053) i
            else
               write(iwr,9050) i,('#',j=1,ifreq(i-ib+1)/ii)
            end if
         end do
         write(iwr,9051)
9049     format(/' Frequency table of eigenvalues of S-matrix',
     1           ' (scaled by ',i3,')',/,1x,80('='))            
9050     format(i5,2x,80a1)
9051     format(1x,80('='),/1x)
9053     format(i5,2x,i1)
      end if
      if (l0.ne.newbas0) then
cjvl     if(l0.ne.newbas0)call caserr2('linear dependance detected')
         newbas0 = l0
         odepen = .true.
c...  get symmetry numbers right
         do i=1,8
            nsymh(i) = nsym0(i)
            nsym0(i) = 0
         end do
         nav = lenwrd()
         call secget(isect(490),51,iblk51)
         call readi(mmmm,mach(13)*nav,iblk51,idaf)
         do i=1,l0
            ibig=idamax(l1,v(1,i),1)
            ii=isymaos(ibig)
            nsym0(ii) = nsym0(ii) + 1 
         end do
c
c        print sabf definition again
c
         write (iwr,6100)
         do loop = 1 , 8
            if (nsym0(loop).gt.0) then
               write (iwr,6110) loop , nsym0(loop)
            end if
         end do
         write (iwr,6120)
 6100 format (/1x,30('=')/1x,'irrep  no. of symmetry adapted'/1x,
     +        '       basis functions'/1x,30('='))
 6110 format (1x,i3,i12)
 6120 format (1x,30('=')/)
c
      endif

c
c...  because l0 may differ from l1 generate a second ilif in ilifq0
      do i=1,l0
         ilifq0(i)=(i-1)*l0
      end do
c
c     ----- form canonical orthonormal orbitals -----
c
       do 200 j=1,l0
      call dscal(l1,e(j),v(1,j),1)
 200   continue
      scr(1)=dble(l1)
      call dcopy(l1,e,1,scr(2),1)
c     clear "non-existent" vectors
      if (l0.lt.l1) call vclr(v(1,l0+1),1,(l1-l0)*l1)
c
c     ----- write the canonical orbitals on ed7  -----
c
      call wrt3(scr,l1+1,iblkv,idaf)
      ibl3qs=iblkv+lenv
      call wrt3(v,l3,ibl3qs,idaf)
400   if (nprint .ne. 5) return
      write (iwr,9068)
      call prev(v,e,l0,l1,ndim)
      return
 9008 format(/,5x,14(1h-),/,5x,14hoverlap matrix,/,5x,14(1h-))
 9028 format(/,5x,17(1h-),/,5x,17heigenvectors of s,/,5x,17(1h-))
 9047 format(//
     b,' ############################################################',/
     2,' ###        possible linear dependence diagnosed          ###',/
     3,' ### Smallest eigenvalue of overlap matrix is', 1pe12.4,' ###',/
     4,' ### There are',i6,' eigenvalue(s) less than ', 1pe12.4,' ###',/
     c,' ############################################################')
 9048 format(//
     a,' ############################################################',/
     b,' ############################################################',/
     2,' ###        possible linear dependence diagnosed          ###',/
     3,' ### Smallest eigenvalue of overlap matrix is', 1pe12.4,' ###',/
     4,' ###',i5,' eigenvector(s) are eliminated (',1pe9.2,
     4 '/',i5,') ###',/
     5,' ### the number of canonical orbitals kept is',i6,'       ###',/
     6,' ###  Beware .... ONLY (not extensively) TESTED on SCF    ###',/
     c,' ############################################################',/
     d,' ############################################################'/)
 9052 format(//,' **the eigenvalues of the overlap matrix (numbered)**',
     1       /,(2x,5(i5,1x,1pe10.3)))
 9068 format(/,5x,30(1h-),/,5x,30hcanonical orthonormal orbitals,/, 5x,
     +     30(1h-))
      end
c
c support code for symmetry blocking of diagonalisation
c
      subroutine gldiag_symm(l0,ldum,l1,h,ib,eig,vector,ia,iop2,
     +                       isymmo, use_symmetry )
c
c     ----- general calling  routine for giveis or ligen
c           ligen version -----
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)

      Integer isymmo( 1:l0 )
      Logical use_symmetry, force_serial

      integer maxlablen
      parameter (maxlablen=20)
      character*20 lab
c
c  change the label definitions (init_time_periodsm in util3 )
c  also
c
      integer TP_ENTIRE
      parameter(TP_ENTIRE=1)
c
c  scf 
c
      integer TP_SCF
      parameter (TP_SCF=TP_ENTIRE+1)
      integer TP_ORFOG
      parameter (TP_ORFOG=TP_SCF+1)
      integer TP_DHSTAR
      parameter (TP_DHSTAR=TP_ORFOG+1)
      integer TP_RDMAT
      parameter (TP_RDMAT=TP_DHSTAR+1)
      integer TP_DIIS
      parameter (TP_DIIS=TP_RDMAT+1)
      integer TP_DIAG
      parameter (TP_DIAG=TP_DIIS+1)
c
c  mp2 code
c
      integer TP_APRDMP2
      parameter (TP_APRDMP2=TP_DIAG+1)

      integer TP_DHSTAR_GOP
      parameter (TP_DHSTAR_GOP=TP_APRDMP2+1)

c spare
      integer TP_APRM1234
      parameter (TP_APRM1234=TP_DHSTAR_GOP+1)
c spare
      integer TP_APRQ34
      parameter (TP_APRQ34=TP_APRM1234+1)

      integer TP_APRQ1
      parameter (TP_APRQ1=TP_APRQ34+1)
      integer TP_APRQ2
      parameter (TP_APRQ2=TP_APRQ1+1)
      integer TP_APRQ2D
      parameter (TP_APRQ2D=TP_APRQ2+1)
      integer TP_APRQ34D
      parameter (TP_APRQ34D=TP_APRQ2D+1)
      integer TP_APRMP2E
      parameter (TP_APRMP2E=TP_APRQ34D+1)

      integer TP_MP1PDM
      parameter (TP_MP1PDM=TP_APRMP2E+1)
      integer TP_MP1PDM_1
      parameter (TP_MP1PDM_1=TP_MP1PDM+1)
      integer TP_MP1PDM_2
      parameter (TP_MP1PDM_2=TP_MP1PDM_1+1)

      integer TP_APR1PDM
      parameter (TP_APR1PDM=TP_MP1PDM_2+1)

      integer TP_MP2HESS
      parameter (TP_MP2HESS=TP_APR1PDM+1)
      integer TP_MP2CHF
      parameter (TP_MP2CHF=TP_MP2HESS+1)
      integer TP_MP2MAKEW
      parameter (TP_MP2MAKEW=TP_MP2CHF+1)
      integer TP_MP2DS
      parameter (TP_MP2DS=TP_MP2MAKEW+1)
      integer TP_MP2BACK_P
      parameter (TP_MP2BACK_P=TP_MP2DS+1)
      integer TP_MP2BACK_F
      parameter (TP_MP2BACK_F=TP_MP2BACK_P+1)
      integer TP_MP2BACKTRAN_2
      parameter (TP_MP2BACKTRAN_2=TP_MP2BACK_F+1)
      integer TP_MP2MCDAB
      parameter (TP_MP2MCDAB=TP_MP2BACKTRAN_2+1)
c
c - parallel functions
c
      integer TP_DGOP
      parameter(TP_DGOP=TP_MP2MCDAB+1)

      integer TP_BCAST
      parameter(TP_BCAST=TP_DGOP+1)

      integer TP_NXTVAL
      parameter(TP_NXTVAL=TP_BCAST+1)

      integer TP_GENRAL
      parameter (TP_GENRAL=TP_NXTVAL+1)

      integer TP_GENRAL_1PDM
      parameter (TP_GENRAL_1PDM=TP_GENRAL+1)

      integer TP_GA_PUT_Q2D
      parameter (TP_GA_PUT_Q2D=TP_GENRAL_1PDM+1)

      integer TP_GA_ACC_Q2D
      parameter (TP_GA_ACC_Q2D=TP_GA_PUT_Q2D+1)

      integer TP_MKT2AO
      parameter (TP_MKT2AO   =TP_GA_ACC_Q2D+1)

      integer TP_APRL2
      parameter (TP_APRL2    =TP_MKT2AO+1)

      integer TP_GA_GET_L2
      parameter (TP_GA_GET_L2=TP_APRL2+1)

      integer TP_APRL34
      parameter (TP_APRL34   =TP_GA_GET_L2+1)

      integer TP_DGENRL
      parameter (TP_DGENRL   =TP_APRL34+1)

      integer TP_MP2
      parameter(TP_MP2       =TP_DGENRL+1)

      integer TP_JKDER
      parameter(TP_JKDER     =TP_MP2+1)

      integer TP_APRL1234
      parameter (TP_APRL1234 =TP_JKDER+1)

      integer TP_APRL1
      parameter (TP_APRL1    =TP_APRL1234+1)

      integer TP_APRDMP2_I
      parameter(TP_APRDMP2_I =TP_APRL1+1)

      integer TP_JKDER_GET
      parameter(TP_JKDER_GET =TP_APRDMP2_I+1)

      integer TP_MP1PDM_3
      parameter(TP_MP1PDM_3  =TP_JKDER_GET+1)

      integer TP_MP1PDM_4
      parameter(TP_MP1PDM_4  =TP_MP1PDM_3+1)

cgdf:  added TP_MKT2MO 14.3.95
      integer TP_MKT2MO
      parameter (TP_MKT2MO   =TP_MP1PDM_4+1)

csecd - time periods for second derivatives
c
      integer TP_2D_AOINTS
      parameter(TP_2D_AOINTS     =TP_MKT2MO+1)

      integer TP_2D_SCF
      parameter(TP_2D_SCF        =TP_2D_AOINTS+1)

      integer TP_2D_HFGRDN
      parameter(TP_2D_HFGRDN     =TP_2D_SCF+1)

      integer TP_2D_INDX2T
      parameter(TP_2D_INDX2T     =TP_2D_HFGRDN+1)

      integer TP_2D_MOINTS
      parameter(TP_2D_MOINTS     =TP_2D_INDX2T+1)

      integer TP_2D_TRNFKD
      parameter(TP_2D_TRNFKD     =TP_2D_MOINTS+1)

      integer TP_2D_CHFNDR
      parameter(TP_2D_CHFNDR     =TP_2D_TRNFKD+1)

      integer TP_2D_QMDER
      parameter(TP_2D_QMDER      =TP_2D_CHFNDR+1)

      integer TP_2D_DMDER
      parameter(TP_2D_DMDER      =TP_2D_QMDER+1)

      integer TP_2D_2D
      parameter(TP_2D_2D         =TP_2D_DMDER+1)

      integer TP_2D_CHF
      parameter(TP_2D_CHF        =TP_2D_2D+1)

      integer TP_2D_NUC
      parameter(TP_2D_NUC        =TP_2D_CHF+1)

      integer TP_2D_OVL
      parameter(TP_2D_OVL        =TP_2D_NUC+1)

      integer TP_2D_KE
      parameter(TP_2D_KE         =TP_2D_OVL+1)

      integer TP_2D_PE
      parameter(TP_2D_PE         =TP_2D_KE+1)

      integer TP_2D_2E
      parameter(TP_2D_2E         =TP_2D_PE+1)

      integer TP_2D_TOTAL
      parameter (TP_2D_TOTAL     =TP_2D_2E+1)

      integer TP_2D_CHFDRV
      parameter (TP_2D_CHFDRV    =TP_2D_TOTAL+1)

      integer TP_2D_PDENS
      parameter (TP_2D_PDENS     =TP_2D_CHFDRV+1)

      integer TP_2D_PFOCK
      parameter (TP_2D_PFOCK     =TP_2D_PDENS+1)

      integer TP_2D_CHFHESS
      parameter (TP_2D_CHFHESS   =TP_2D_PFOCK+1)

      integer TP_2D_CHFRHS
      parameter (TP_2D_CHFRHS    =TP_2D_CHFHESS+1)

      integer TP_2D_SYMMRHS
      parameter (TP_2D_SYMMRHS   =TP_2D_CHFRHS+1)

      integer TP_2D_SYMMU
      parameter (TP_2D_SYMMU     =TP_2D_SYMMRHS+1)

      integer TP_2D_PFOCK_OOOO
      parameter (TP_2D_PFOCK_OOOO=TP_2D_SYMMU+1)

      integer TP_2D_PFOCK_VOOO
      parameter (TP_2D_PFOCK_VOOO=TP_2D_PFOCK_OOOO+1)

      integer TP_2D_PFOCK_VVOO
      parameter (TP_2D_PFOCK_VVOO=TP_2D_PFOCK_VOOO+1)

      integer TP_2D_PFOCK_VOVO
      parameter (TP_2D_PFOCK_VOVO=TP_2D_PFOCK_VVOO+1)

      integer TP_2D_PFOCK_SUM
      parameter (TP_2D_PFOCK_SUM =TP_2D_PFOCK_VOVO+1)

      integer TP_2D_AOGEN
      parameter(TP_2D_AOGEN=TP_2D_PFOCK_SUM+1)

      integer TP_2D_AOUT
      parameter(TP_2D_AOUT =TP_2D_AOGEN+1)

      integer TP_TEST1
      parameter(TP_TEST1   =TP_2D_AOUT+1)

      integer TP_TEST2
      parameter(TP_TEST2   =TP_TEST1+1)

      integer TP_TEST3
      parameter(TP_TEST3   =TP_TEST2+1)

      integer TP_TEST4
      parameter(TP_TEST4   =TP_TEST3+1)

      integer TP_TEST5
      parameter(TP_TEST5   =TP_TEST4+1)

      integer TP_TEST6
      parameter(TP_TEST6   =TP_TEST5+1)

      integer TP_HFGRAD
      parameter(TP_HFGRAD  =TP_TEST6+1)

      integer TP_GAMULT2
      parameter(TP_GAMULT2 =TP_HFGRAD+1)

      integer TP_GAORTHOG
      parameter(TP_GAORTHOG=TP_GAMULT2+1)

      integer TP_PDIAG
      parameter (TP_PDIAG  =TP_GAORTHOG+1)

      integer TP_MULT2
      parameter (TP_MULT2  =TP_PDIAG+1)

      integer TP_INTEG
      parameter (TP_INTEG  =TP_MULT2+1)
c
c =================  I/O timers ========================
c
c find
c	
      integer TP_IOFM1, TP_IOF0, TP_IOF1, TP_IOF2, TP_IOF3,
     &	TP_IOF4, TP_IOF5, TP_IOF6, TP_IOF7

      parameter(TP_IOFM1=TP_INTEG+1)
      parameter (TP_IOF0=TP_IOFM1+1)
      parameter (TP_IOF1=TP_IOF0+1)
      parameter (TP_IOF2=TP_IOF1+1)
      parameter (TP_IOF3=TP_IOF2+1)
      parameter (TP_IOF4=TP_IOF3+1)
      parameter (TP_IOF5=TP_IOF4+1)
      parameter (TP_IOF6=TP_IOF5+1)
      parameter (TP_IOF7=TP_IOF6+1)
c
c get
c
      integer TP_IOGM1, TP_IOG0, TP_IOG1, TP_IOG2, TP_IOG3, 
     &  TP_IOG4, TP_IOG5, TP_IOG6, TP_IOG7
      parameter (TP_IOGM1=TP_IOF7+1)
      parameter (TP_IOG0=TP_IOGM1+1)
      parameter (TP_IOG1=TP_IOG0+1)
      parameter (TP_IOG2=TP_IOG1+1)
      parameter (TP_IOG3=TP_IOG2+1)
      parameter (TP_IOG4=TP_IOG3+1)
      parameter (TP_IOG5=TP_IOG4+1)
      parameter (TP_IOG6=TP_IOG5+1)
      parameter (TP_IOG7=TP_IOG6+1)
c
c put
c
      integer TP_IOPM1,  TP_IOP0, TP_IOP1, TP_IOP2, TP_IOP3,
     & TP_IOP4, TP_IOP5, TP_IOP6, TP_IOP7
      parameter (TP_IOPM1=TP_IOG7+1)
      parameter (TP_IOP0=TP_IOPM1+1)
      parameter (TP_IOP1=TP_IOP0+1)
      parameter (TP_IOP2=TP_IOP1+1)
      parameter (TP_IOP3=TP_IOP2+1)
      parameter (TP_IOP4=TP_IOP3+1)
      parameter (TP_IOP5=TP_IOP4+1)
      parameter (TP_IOP6=TP_IOP5+1)
      parameter (TP_IOP7=TP_IOP6+1)
c
c open
c
      integer TP_IOOM1,TP_IOO0,TP_IOO1,TP_IOO2,TP_IOO3,
     & TP_IOO4,TP_IOO5, TP_IOO6, TP_IOO7

      parameter (TP_IOOM1=TP_IOP7+1)
      parameter (TP_IOO0=TP_IOOM1+1)
      parameter (TP_IOO1=TP_IOO0+1)
      parameter (TP_IOO2=TP_IOO1+1)
      parameter (TP_IOO3=TP_IOO2+1)
      parameter (TP_IOO4=TP_IOO3+1)
      parameter (TP_IOO5=TP_IOO4+1)
      parameter (TP_IOO6=TP_IOO5+1)
      parameter (TP_IOO7=TP_IOO6+1)
c
c delfil (only significant for GA-files dumped to disc
c
      integer TP_IO_GAFILE_READ, TP_IO_GAFILE_DUMP
      parameter (TP_IO_GAFILE_READ      =TP_IOO7+1)
      parameter (TP_IO_GAFILE_DUMP      =TP_IO_GAFILE_READ+1)
c
c Peigs parallel diag
c
      integer TP_PEIGS
      parameter (TP_PEIGS               =TP_IO_GAFILE_DUMP+1)
c
c Scalapack parallel diag
c
      integer TP_PDSYEV
      parameter (TP_PDSYEV               =TP_PEIGS+1)
      integer TP_PDSYEVX
      parameter (TP_PDSYEVX              =TP_PDSYEV+1)
      integer TP_PDSYEVD
      parameter (TP_PDSYEVD              =TP_PDSYEVX+1)
      integer TP_PDSYEVR
      parameter (TP_PDSYEVR              =TP_PDSYEVD+1)
c
c timers for CCP1 DFT module
c
      integer    TP_DFT_JFIT
      parameter (TP_DFT_JFIT            =TP_PDSYEVR+1)
      integer    TP_DFT_JFIT_VFORM
      parameter (TP_DFT_JFIT_VFORM      =TP_DFT_JFIT+1)
      integer    TP_DFT_JFIT_TR
      parameter (TP_DFT_JFIT_TR         =TP_DFT_JFIT_VFORM+1)
      integer    TP_DFT_JFIT_NR
      parameter (TP_DFT_JFIT_NR         =TP_DFT_JFIT_TR+1)
      integer    TP_DFT_JFIT_COEF
      parameter (TP_DFT_JFIT_COEF       =TP_DFT_JFIT_NR+1)
      integer    TP_DFT_JFIT_KSMAT
      parameter (TP_DFT_JFIT_KSMAT      =TP_DFT_JFIT_COEF+1)
      integer    TP_DFT_JFIT_ENERGY
      parameter (TP_DFT_JFIT_ENERGY     =TP_DFT_JFIT_KSMAT+1)
      integer    TP_DFT_EXQUAD
      parameter (TP_DFT_EXQUAD          =TP_DFT_JFIT_ENERGY+1)
      integer    TP_DFT_EXQUAD_INTRO
      parameter (TP_DFT_EXQUAD_INTRO    =TP_DFT_EXQUAD+1)
      integer    TP_DFT_EXQUAD_INTEG
      parameter (TP_DFT_EXQUAD_INTEG    =TP_DFT_EXQUAD_INTRO+1)
      integer    TP_DFT_EXQUAD_DGOP
      parameter (TP_DFT_EXQUAD_DGOP     =TP_DFT_EXQUAD_INTEG+1)
      integer    TP_DFT_EXQUADF
      parameter (TP_DFT_EXQUADF         =TP_DFT_EXQUAD_DGOP+1)
      integer    TP_DFT_EXQUADF_INTRO
      parameter (TP_DFT_EXQUADF_INTRO   =TP_DFT_EXQUADF+1)
      integer    TP_DFT_EXQUADF_INTEG
      parameter (TP_DFT_EXQUADF_INTEG   =TP_DFT_EXQUADF_INTRO+1)
      integer    TP_DFT_EXQUADF_DGOP
      parameter (TP_DFT_EXQUADF_DGOP    =TP_DFT_EXQUADF_INTEG+1)
      integer    TP_DFT_EXQUADLHS
      parameter (TP_DFT_EXQUADLHS       =TP_DFT_EXQUADF_DGOP+1)
      integer    TP_DFT_EXQUADLHS_INTRO
      parameter (TP_DFT_EXQUADLHS_INTRO =TP_DFT_EXQUADLHS+1)
      integer    TP_DFT_EXQUADLHS_INTEG
      parameter (TP_DFT_EXQUADLHS_INTEG =TP_DFT_EXQUADLHS_INTRO+1)
      integer    TP_DFT_EXQUADLHS_DGOP
      parameter (TP_DFT_EXQUADLHS_DGOP  =TP_DFT_EXQUADLHS_INTEG+1)
      integer    TP_DFT_EXQUADHES
      parameter (TP_DFT_EXQUADHES       =TP_DFT_EXQUADLHS_DGOP+1)
      integer    TP_DFT_EXQUADHES_INTRO
      parameter (TP_DFT_EXQUADHES_INTRO =TP_DFT_EXQUADHES+1)
      integer    TP_DFT_EXQUADHES_INTEG
      parameter (TP_DFT_EXQUADHES_INTEG =TP_DFT_EXQUADHES_INTRO+1)
      integer    TP_DFT_EXQUADHES_DGOP
      parameter (TP_DFT_EXQUADHES_DGOP  =TP_DFT_EXQUADHES_INTEG+1)
      integer    TP_DFT_EXQUADRHS
      parameter (TP_DFT_EXQUADRHS       =TP_DFT_EXQUADHES_DGOP+1)
      integer    TP_DFT_EXQUADRHS_INTRO
      parameter (TP_DFT_EXQUADRHS_INTRO =TP_DFT_EXQUADRHS+1)
      integer    TP_DFT_EXQUADRHS_INTEG
      parameter (TP_DFT_EXQUADRHS_INTEG =TP_DFT_EXQUADRHS_INTRO+1)
      integer    TP_DFT_EXQUADRHS_DGOP
      parameter (TP_DFT_EXQUADRHS_DGOP  =TP_DFT_EXQUADRHS_INTEG+1)
      integer    TP_DFT_EXQUADDKSX
      parameter (TP_DFT_EXQUADDKSX      =TP_DFT_EXQUADRHS_DGOP+1)
      integer    TP_DFT_EXQUADDKSX_INTRO
      parameter (TP_DFT_EXQUADDKSX_INTRO=TP_DFT_EXQUADDKSX+1)
      integer    TP_DFT_EXQUADDKSX_INTEG
      parameter (TP_DFT_EXQUADDKSX_INTEG=TP_DFT_EXQUADDKSX_INTRO+1)
      integer    TP_DFT_EXQUADDKSX_DGOP
      parameter (TP_DFT_EXQUADDKSX_DGOP =TP_DFT_EXQUADDKSX_INTEG+1)
      integer    TP_DFT_EXQUADDKS
      parameter (TP_DFT_EXQUADDKS       =TP_DFT_EXQUADDKSX_DGOP+1)
      integer    TP_DFT_EXQUADDKS_INTRO
      parameter (TP_DFT_EXQUADDKS_INTRO =TP_DFT_EXQUADDKS+1)
      integer    TP_DFT_EXQUADDKS_INTEG
      parameter (TP_DFT_EXQUADDKS_INTEG =TP_DFT_EXQUADDKS_INTRO+1)
      integer    TP_DFT_EXQUADDKS_DGOP
      parameter (TP_DFT_EXQUADDKS_DGOP  =TP_DFT_EXQUADDKS_INTEG+1)



      integer    TP_DFT_JMULT
      parameter (TP_DFT_JMULT           =TP_DFT_EXQUADDKS_DGOP+1)
      integer    TP_DFT_JMULT_INTRO
      parameter (TP_DFT_JMULT_INTRO     =TP_DFT_JMULT+1)
      integer    TP_DFT_JMULT_SB
      parameter (TP_DFT_JMULT_SB        =TP_DFT_JMULT_INTRO+1)
      integer    TP_DFT_JMULT_FOCK
      parameter (TP_DFT_JMULT_FOCK      =TP_DFT_JMULT_SB+1)
      integer    TP_DFT_JMULT_CJAT0
      parameter (TP_DFT_JMULT_CJAT0     =TP_DFT_JMULT_FOCK+1)
c
c coulomb fitted gradients
c
      integer    TP_DFT_JFITG
      parameter (TP_DFT_JFITG           =TP_DFT_JMULT_CJAT0+1)

      integer    TP_DFT_JFITG_VFORM
      parameter (TP_DFT_JFITG_VFORM     =TP_DFT_JFITG+1)
      integer    TP_DFT_JFITG_TR
      parameter (TP_DFT_JFITG_TR        =TP_DFT_JFITG_VFORM+1)
      integer    TP_DFT_JFITG_NR
      parameter (TP_DFT_JFITG_NR        =TP_DFT_JFITG_TR+1)
      integer    TP_DFT_JFITG_COEF
      parameter (TP_DFT_JFITG_COEF      =TP_DFT_JFITG_NR+1)
      integer    TP_DFT_JFITG_2C
      parameter (TP_DFT_JFITG_2C        =TP_DFT_JFITG_COEF+1)
      integer    TP_DFT_JFITG_3C
      parameter (TP_DFT_JFITG_3C        =TP_DFT_JFITG_2C+1)
      integer    TP_DFT_JFIT_TR_INIT
      parameter (TP_DFT_JFIT_TR_INIT    =TP_DFT_JFITG_3C+1)
      integer    TP_DFT_JFIT_INV
      parameter (TP_DFT_JFIT_INV        =TP_DFT_JFIT_TR_INIT+1)
c
c VB
c
      integer TP_VB
      parameter (TP_VB                  =TP_DFT_JFIT_INV+1)
      integer TP_VB_STRUC
      parameter (TP_VB_STRUC            =TP_VB+1)
      integer TP_VB_ME
      parameter (TP_VB_ME               =TP_VB_STRUC+1)
      integer TP_VB_DIAG
      parameter (TP_VB_DIAG             =TP_VB_ME+1)
      integer TP_VB_TRAN
      parameter (TP_VB_TRAN             =TP_VB_DIAG+1)
      integer TP_VB_VIRT
      parameter (TP_VB_VIRT             =TP_VB_TRAN+1)
      integer TP_VB_LADM
      parameter (TP_VB_LADM             =TP_VB_VIRT+1)
      integer TP_VB_DTRAN
      parameter (TP_VB_DTRAN            =TP_VB_LADM+1)
c
c VB parallel extra
c
      integer TP_ISEND
      parameter (TP_ISEND               =TP_VB_DTRAN+1)
      integer TP_IRECV
      parameter (TP_IRECV               =TP_ISEND+1)
      integer TP_WAIT
      parameter (TP_WAIT                =TP_IRECV+1)
c
c One electron derivatives
c
      integer TP_STVECP
      parameter (TP_STVECP              =TP_WAIT+1)

      integer TP_TVDER
      parameter (TP_TVDER               =TP_STVECP+1)

      integer TP_SDER
      parameter (TP_SDER                =TP_TVDER+1)

      integer TP_SGRAD
      parameter (TP_SGRAD               =TP_SDER+1)

      integer TP_HELFEY
      parameter (TP_HELFEY              =TP_SGRAD+1)


      !F90 time periods start here
      integer TP_F90_START
      parameter( TP_F90_START           = TP_HELFEY+1 )
      integer TP_F90_SCF
      parameter( TP_F90_SCF             = TP_F90_START+1 )
      integer TP_F90_BUILD
      parameter( TP_F90_BUILD           = TP_F90_SCF+1 )
      integer TP_F90_DIIS
      parameter( TP_F90_DIIS            = TP_F90_BUILD+1 )
      integer TP_F90_SIMIL
      parameter( TP_F90_SIMIL           = TP_F90_DIIS+1 )
      integer TP_F90_DIAG
      parameter( TP_F90_DIAG            = TP_F90_SIMIL+1 )
      integer TP_F90_BACK
      parameter( TP_F90_BACK            = TP_F90_DIAG+1 )
      integer TP_F90_ASSIGN
      parameter( TP_F90_ASSIGN          = TP_F90_BACK+1 )
      integer TP_F90_ORTHOG
      parameter( TP_F90_ORTHOG          = TP_F90_ASSIGN+1 )
      integer TP_F90_MAKE_DENS
      parameter( TP_F90_MAKE_DENS       = TP_F90_ORTHOG+1 )
      integer TP_F90_END
      parameter( TP_F90_END             = TP_F90_MAKE_DENS+1 )
      integer TP_F90_LEV_SHIFT
      parameter( TP_F90_LEV_SHIFT       = TP_F90_END+1 )
      integer TP_F90_TESTER_EVAL
      parameter( TP_F90_TESTER_EVAL     = TP_F90_LEV_SHIFT+1 )
      integer TP_F90_DELTA_EVAL
      parameter( TP_F90_DELTA_EVAL      = TP_F90_TESTER_EVAL+1 )
      integer TP_F90_TDOWN
      parameter( TP_F90_TDOWN           = TP_F90_DELTA_EVAL+1 )
      integer TP_F90_RDMAT
      parameter( TP_F90_RDMAT           = TP_F90_TDOWN+1 )
      integer TP_F90_INTS
      parameter( TP_F90_INTS            = TP_F90_RDMAT+1 )
      integer TP_NEWSCF
      parameter( TP_NEWSCF              = TP_F90_INTS+1 )

      integer TP_DENSCF
      parameter( TP_DENSCF              = TP_NEWSCF+1 )
      integer TP_DENSCF_BUILD
      parameter( TP_DENSCF_BUILD        = TP_DENSCF+1 )
      integer TP_DENSCF_RDMAT
      parameter( TP_DENSCF_RDMAT        = TP_DENSCF_BUILD+1 )
      integer TP_DENSCF_INTS
      parameter( TP_DENSCF_INTS         = TP_DENSCF_RDMAT+1 )
      integer TP_DENSCF_DIAG_S
      parameter( TP_DENSCF_DIAG_S       = TP_DENSCF_INTS+1 )
      integer TP_DENSCF_SIMIL
      parameter( TP_DENSCF_SIMIL        = TP_DENSCF_DIAG_S+1 )
      integer TP_DENSCF_DIAG
      parameter( TP_DENSCF_DIAG         = TP_DENSCF_SIMIL+1 )
      integer TP_DENSCF_BACK
      parameter( TP_DENSCF_BACK         = TP_DENSCF_DIAG+1 )
      integer TP_DENSCF_MAKE_DENS
      parameter( TP_DENSCF_MAKE_DENS    = TP_DENSCF_BACK+1 )
      integer TP_DENSCF_TDOWN
      parameter( TP_DENSCF_TDOWN        = TP_DENSCF_MAKE_DENS+1 )

      integer TP_DRHFCL_GA
      parameter( TP_DRHFCL_GA           = TP_DENSCF_TDOWN+1 )
c
c     RPA module
c
      integer TP_RESPONSE
      parameter( TP_RESPONSE            = TP_DRHFCL_GA + 1)
      integer TP_RPA
      parameter( TP_RPA                 = TP_RESPONSE + 1)
      integer TP_TDA
      parameter( TP_TDA                 = TP_RPA + 1)
      integer TP_RPANAL
      parameter( TP_RPANAL              = TP_TDA + 1)
      integer TP_RPA_MO2AO
      parameter( TP_RPA_MO2AO           = TP_RPANAL + 1)
      integer TP_RPA_INT
      parameter( TP_RPA_INT             = TP_RPA_MO2AO + 1)
      integer TP_RPA_CNTRCT
      parameter( TP_RPA_CNTRCT          = TP_RPA_INT + 1)
      integer TP_RPA_AO2MO
      parameter( TP_RPA_AO2MO           = TP_RPA_CNTRCT + 1)
      integer TP_TDA_MO2AO
      parameter( TP_TDA_MO2AO           = TP_RPA_AO2MO + 1)
      integer TP_TDA_INT
      parameter( TP_TDA_INT             = TP_TDA_MO2AO + 1)
      integer TP_TDA_CNTRCT
      parameter( TP_TDA_CNTRCT          = TP_TDA_INT + 1)
      integer TP_TDA_AO2MO
      parameter( TP_TDA_AO2MO           = TP_TDA_CNTRCT + 1)
c
c     Define the common blocks
c
      integer maxtp
      parameter (maxtp=TP_TDA_AO2MO)
      integer maxtpdepth
      parameter (maxtpdepth = 10)
      integer itpdepth
      integer itpstack
      integer ntpc
      integer parent
      real*8  ttotw, ttotc, tsw, tsc
      real*8  ttotu, ttots, tsu, tss
      real*8  taggc
      common/timeperiods/ttotw(maxtp),ttotc(maxtp),
     &     ttotu(maxtp),ttots(maxtp),
     &     tsw(maxtp),tsc(maxtp),
     &     tsu(maxtp),tss(maxtp),
     &     taggc(maxtp),
     &     ntpc(maxtp),parent(maxtp),
     &     itpstack(0:maxtpdepth),itpdepth
      common/timeperiodsc/lab(maxtp)

c...   common for overruling local criteria (e.g. dependency, diag thresholds)
      integer  i_depen_check, i_global_diag, n_depen_check
      logical  o_depen_print
      common/overrule/ i_depen_check,i_global_diag, o_depen_print, 
     1                 n_depen_check
c
c     i_depen_check : criterion for throwing out vectors
c     n_depen_check : number  of vectors to be thrown out 
c     o_depen_print : if true +> extra print 

      dimension vector(*),h(*),ib(*),eig(*),ia(*)
      do 10 i=1,l0
 10   ib(i)=(i-1)*ldum
c
c     WARNING: original setting of 1.0d-8 caused real problems 
c     in some high symmetry cases .. reset diaacc from
c     th=1.0d-08
c     the value below is that currently used in gldiag. This may be
c     rather aggressive for some basis sets ..
c
      th=1.0d-15
      if (i_global_diag.ne.-999) th = 10.0d0**(-i_global_diag*1.0d0)
c  
      iop1=2
      Call start_time_period( TP_DIAG )
      call jacobi_symm(h,ia,l0,vector,ib,l1,eig,iop1,iop2,th,
     +                 isymmo, use_symmetry)
      Call end_time_period( TP_DIAG )
      return
      end
*
***************************************************************
*
      Subroutine block_a( n, n_irreps, symmetry, diag_order, 
     +                    orbs_this_symm, a, a_point, q_point,
     +                    work )
*
*  Given the Hamiltonian in A which has been built from
* orbitals that span N_IRREPS irreps and a given orbital
* spans the irrep given by SYMMETRY this routine permutes
* the rows and columns of A such that it becomes block
* diagonal, and also provides pointers for the columns
* of A in the new ordering and where the columns of the
* evecs in the new ordering should go.
* WARNING: there be nasties 'round 'ere !
*
*  The ordering of the blocks in A is such that the
* block corresponding to irrep DIAG_ORDER( 1, 1 )
* of size DIAG_ORDER( 2, 1 ) comes first, then directly
* after that comes the block for DIAG_ORDER( 1, 2 )
* etc. The pointers to the columns of A point
* to the position down the column where the block
* starts ( i.e. leading zeros are ignored ). Similarly
* for the pointers for the evecs.
*
      Implicit None
*
      Integer n
      Integer n_irreps
      Integer symmetry( 1:n )
      Integer diag_order( 1:2, 1:8 )
      Integer orbs_this_symm( 1:n )
      real*8    a( 1: ( n * ( n + 1 ) ) / 2 )
      Integer a_point( 1:n )
      Integer q_point( 1:n )
      real*8    work( 1: ( n * ( n + 1 ) ) / 2 )
*
      Integer size_irrep, which_irrep
      Integer point
      Integer irrep_offset
      Integer column
      Integer a_orig_start
      Integer a_entries
      Integer i, j, k
*
*  Generate the pointers to the starts of the columns
* in the block decomposed A. 
*
      irrep_offset = 0
      column       = 1
*
      Do i = 1, n_irreps
         point = 1
         size_irrep = diag_order( 2, i )
         Do j = 1, size_irrep
            a_point( column ) = point + irrep_offset
            point  = point  + j
            column = column + 1
         End Do
         irrep_offset = irrep_offset + 
     +        ( size_irrep * ( size_irrep + 1 ) ) / 2
      End Do
*
*  Generate the pointers into the eigenvector
* array
*
      irrep_offset = 0
      column       = 1
      point        = 1
*
      Do i = 1, n_irreps
         size_irrep = diag_order( 2, i )
         Do j = 1, size_irrep
            q_point( column ) = point + irrep_offset
            point  = point  + n
            column = column + 1
         End Do
         irrep_offset = irrep_offset + size_irrep
      End Do
*
*  O.K., for each irrep in turn find the elements
* of each column that span that irrep and store
* in the correct order in the workspace.
*
      column    = 1
      a_entries = 0
*
      Do i = 1, n_irreps
*
         which_irrep = diag_order( 1, i )
         size_irrep  = diag_order( 2, i )
*
         point = 1
         Do j = 1, n
            If( symmetry( j ) .EQ. which_irrep ) Then
               orbs_this_symm( point ) = j
               point = point + 1
            End If
         End Do
*
         Do j = 1, size_irrep
            point        = a_point( column )
            a_orig_start = 
     +   ( orbs_this_symm( j ) * ( orbs_this_symm( j ) - 1 ) ) / 2 
            Do k = 1, j
               work( point ) = a( a_orig_start + orbs_this_symm( k ) )
               point         = point  + 1
               a_entries     = a_entries + 1
            End Do
            column = column + 1
         End Do
*
      End Do
*
*  Copy the resulting block diagonal Hamiltonian back into A
*
*
      Call dcopy(a_entries, work, 1, a, 1)
*
      End
*
***************************************************************
*
      Subroutine shift_evecs( n, n_irrep, q, q_point )
*
*  Both the serial and parallel diags assume that where
* the evecs should be stored is one contiguous block
* in memory. This is not true if the symmetry of the
* molecule is being used to block diagonalize
* the Hamiltonian. This routine rectifies this
* by shifting the evecs to the correct places and zeroing
* the parts which should be strictly zero by symmetry.  
*
      Implicit None
*
      Integer n
      Integer n_irrep
      real*8    q( 1:n * n )
      Integer q_point( 1:n_irrep )
*     
      Integer q_from, q_to
      Integer i, j
*
*  Have to be a little careful. Move in backwards
* order to avoid overwrites.
*
      q_from = q_point( 1 ) + n_irrep * n_irrep - 1
*
      Do i = n_irrep, 2, -1
         q_to = q_point( i ) + n_irrep - 1
         Do j = n_irrep, 1, -1
*
            q( q_to   ) = q( q_from )
            q( q_from ) = 0.0d0
*
            q_to   = q_to   - 1
            q_from = q_from - 1
*
         End Do
      End Do
*
      End
*
***************************************************************
*
      Subroutine restore_orig_symmetry( n, n_irreps, q, 
     +                                  basis_symmetry, diag_order,
     +                                  work )
*
*  This undoes the permutations used to turn the
* Hamiltonian into block diagonal symmetry and applies
* them to the evecs to ensure their symmetry is
* consistent with the symmetry of the orbitals
* used to build the Hamiltonian ( the permuting
* made the labelling sequence for the orbitals
* different ).
*
*  One complication is the lack of a n * n workspace
* array, we've only got n * ( n + 1 ) / 2. Hence we
* do it in two goes, n/2 columns at a time.
*
      Implicit None
*
      Integer n
      Integer n_irreps
      real*8    q( 1:n, 1:n )
      Integer basis_symmetry( 1:n )
      Integer diag_order( 1:2, 1:8 )
      real*8    work( 1:n, 1:( n + 1 ) / 2 )
*
      Integer offsets_symm( 1:8 )
      Integer n_shift_1, n_shift_2
      Integer i, j
*
*  First of all work out where in a given column each
* irrep starts.
*
      Do i = 1, n_irreps
         offsets_symm( i ) = 1
         j = 1
         Do While( diag_order( 1, j ) .NE. i )
            offsets_symm( i ) = offsets_symm( i ) + diag_order( 2, j )
            j = j + 1
         End Do
      End Do
*
*  Now do the two shifts.
*
      n_shift_1 = n / 2
      n_shift_2 = n - n / 2
*
      Call shift_to_orig_symmetry( n, n_shift_1, q, 
     +                             basis_symmetry,
     +                             offsets_symm, work )
      Call shift_to_orig_symmetry( n, n_shift_2, q( 1, n_shift_1 + 1 ), 
     +                             basis_symmetry,
     +                             offsets_symm, work )
*
      End
*
***************************************************************
*
      Subroutine shift_to_orig_symmetry( n, n_shift, q, basis_symmetry,
     +                                   offsets_symm, work )
*
*  Shifts a set of blocked by symmetry vectors back to the
* original ordering as described by BASIS_SYMMETRY.
*
      Integer n
      Integer n_shift
      real*8    q( 1:n, 1:n_shift )
      Integer basis_symmetry( 1:n )
      Integer offsets_symm( 1:8 )
      real*8    work( 1:n, 1:n_shift )
*
      Integer count_symm( 1:8 )
      Integer this_symm
      Integer which_row
      Integer i
*
      Do i = 1, 8
         count_symm( i ) = 0
      End Do
*
      Do i = 1, n
         this_symm = basis_symmetry( i )
         which_row = offsets_symm( this_symm ) + 
     +                 count_symm( this_symm )
         call dcopy(n_shift, q( which_row, 1 ), n, work( i, 1 ), n)
         count_symm( this_symm ) = count_symm( this_symm ) + 1
      End Do
*
      Call dcopy(n * n_shift, work, 1, q, 1)
*
      End

      subroutine symvec(q,isymao,isymmo,nbas,ncol,iblk)
      implicit real*8  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
      dimension q(*),isymmo(*),isymao(*)
c     nbsq=nbas*ncol (last vector is there; but zeroed)
      nbsq=nbas*nbas
      iblkv=iblk+1+lensec(mach(8))+lensec(mach(9))
      call rdedx(q,nbsq,iblkv,idaf)
      do 10 i=1,ncol
      ii=ilifq(i)
      do 40 j=1,nbas
      if( dabs(q(ii+j)).gt.1.0d-3)go to 50
 40   continue
 50   isymmo(i)=isymao(j)
 10   continue
      return
      end

      function jsubst(ztext)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension xchar(8)
      do 1 i=1,8
 1    xchar(i)=ztext(i:i)
      jsubst = isubst(xchar)
      return
      end
      subroutine gamgn2
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c        *****  computes and tabulates f0(x) to f5(x)           *****
c        *****  in range x = -0.15 to x = 19.9                  *****
c        *****  in units of x = 0.05                            *****
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
      common/junk/c(1200,6)
      dimension qy(410)
      data pt184,pt5/ 0.184d0,0.50d0/
      data six,tenm7/6.0d0,1.0d-20 /
      data four,two,done/4.0d0,2.0d0,1.0d0 /
      data pt886/0.8862269254527d0 /
      data pt15,pt05/0.15d0,0.05d0/
      data m22,mword/22,7200/
      q=-done
      do 30 mm=1,6
      m=mm-1
      q=q+done
      qx = -pt15
      do 20 i=1,410
      qx = qx+pt05
      a=q
c        *****  change limit of approximate solution.           *****
      if(qx-15.0d0) 1,1,10
    1 a=a+pt5
      term=done/a
      ptlsum=term
      do 2 l=2,50
      a=a+done
      term=term*qx/a
      ptlsum=ptlsum+term
      if(dabs(term/ptlsum)-tenm7)3,2,2
    2 continue
    3 qy(i)=pt5*ptlsum*dexp(-qx)
      go to 20
   10 b=a+pt5
      a=a-pt5
      approx=pt886/(dsqrt(qx)*qx**m)
      if(m.eq.0) go to 13
      do 12 l=1,m
      b=b-done
   12 approx=approx*b
   13 fimult=pt5*dexp(-qx)/qx
      fiprop=fimult/approx
      term=done
      ptlsum=term
      notrms=qx
      notrms=notrms+m
      do 14 l=2,notrms
      term=term*a/qx
      ptlsum=ptlsum+term
      if(dabs(term*fiprop/ptlsum)-tenm7)15,15,14
   14 a=a-done
   15 qy(i)=approx-fimult*ptlsum
   20 continue
      do 30 i=1,400
      j=i+2
      c(i,mm)=qy(j)
      c(i+400,mm)=qy(j+1)-qy(j)
      temp1=-two*qy(j)+qy(j+1)+qy(j-1)
      temp2=six*qy(j)-four*qy(j+1)-four*qy(j-1)+qy(j-2)+qy(j+2)
   30 c(i+800,mm) = (temp1-pt184*temp2)/six
      call secput(isect(503),m22,lensec(mword),iblk22)
      call wrt3(c,mword,iblk22,idaf)
      return
      end
      subroutine qmat(q,s,v,e,scr,ia,l0,l1,l3,ndim,out)
      implicit real*8  (a-h,p-w),integer    (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
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
      integer ibl3s, ibl3t, ibl3f, ibl3x, ibl3y, ibl3z
      integer ibl3rs, ibl3qa, ibl3pa, ibl3ea, ibl3qb, ibl3pb
      integer ibl3eb, ibl3g, ibl3hs, ionsec, ibl3op, ions2
      integer isecqm
      common /dump3/ ibl3s, ibl3t, ibl3f, ibl3x, ibl3y, ibl3z,
     +               ibl3rs,ibl3qa,ibl3pa,ibl3ea,ibl3qb,ibl3pb,
     +               ibl3eb,ibl3g,ibl3hs,ionsec,ibl3op,ions2,
     +               isecqm
c
c
      integer ibl7s, ibl7t, ibl7f, ibl7st,ibl7tt,ibl7ft,ibl3qs
      integer ibl7qa,ibl7pa,ibl7ea,ibl7qb,ibl7pb,ibl7eb
      integer ibl7x, ibl7y, ibl7z
      integer ibl7so,ibl7o, ibl7la,ibl7ec,ibl7ha
      common /scra7/ ibl7s, ibl7t, ibl7f, ibl7st,ibl7tt,ibl7ft,ibl3qs,
     +               ibl7qa,ibl7pa,ibl7ea,ibl7qb,ibl7pb,ibl7eb,
     +               ibl7x, ibl7y, ibl7z,
     +               ibl7so,ibl7o(8), ibl7la, ibl7ec, ibl7ha
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
c...   common for harmonic option
      logical  oharm,opharm,odepen
      integer newbas0, newbas1, nsym0, ilifq0, ielimh
      integer newbash,nsymh
      common/harmon/ oharm,opharm,newbas0,newbas1,nsym0(8),
     1               ilifq0(maxorb),ielimh(maxorb),
     2               newbash,nsymh(8),odepen
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
c...   common for overruling local criteria (e.g. dependency, diag thresholds)
      integer  i_depen_check, i_global_diag, n_depen_check
      logical  o_depen_print
      common/overrule/ i_depen_check,i_global_diag, o_depen_print, 
     1                 n_depen_check
c
c     i_depen_check : criterion for throwing out vectors
c     n_depen_check : number  of vectors to be thrown out 
c     o_depen_print : if true +> extra print 

c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
      common/craypk/mmmm(65),isymaos(maxorb),isymmos(maxorb)
      common/restri/jjfile(63),lds(508),isect(508),ldsect(508),
     + iacsct(508)
      integer itab
      character*7 fnm
      character*4 snm
      data fnm,snm/'util1.m','qmat'/
      dimension q(*),s(*),v(ndim,*),e(*),scr(*),ia(*)
      data m167,done,tol1,tol0 /167,1.0d0,1.0d-05,1.0d-07/
c
c...  allow specfied crit for dependency (10**-i_depen_check)
c...  (0 => no checking) ; i_depen_check in common/overrule/
c
      if (i_depen_check.ge.0) then
         tol0 = 10.0d0**(-1.0d0*i_depen_check)
         if (i_depen_check.eq.0) tol0 = 1.0d-20
      end if
c
c     ----- diagonalize overlap matrix -----
c
      if (out) then
      write (iwr,9008)
      call prtril(s,l1)
      endif
      len1=lensec(l1+1)
      l2=l1*(l1+1)/2
      lenv=len1+lensec(l2)
      call sectst(isecqm,itest)
      if(itest.eq.0) go to 300
      call secget(isecqm,m167,iblkv)
      call rdedx(scr,l1+1,iblkv,idaf)
      if(scr(1).ne.dble(l1)) go to 300
      call dcopy(l1,scr(2),1,e,1)
      call reads(scr,l2,idaf)
      call vsub(s,1,scr,1,scr,1,l2)

      do 310 i=1,l2
      if( dabs(scr(i)).gt.tol1) go to 300
 310  continue
      call reads(v,l3,idaf)
      ibl3qs=iblkv+lenv
      l0=newbas0
      go to 400

300   lenq=lenv+lensec(l3)
      call secput(isecqm,m167,lenq,iblkv)
      call wrt3(s,l2,iblkv+len1,idaf)
c
      if (oharm) call comharm(s,'tri',scr)
      call gldiag(newbas0,newbas0,newbas0,s,scr,e,v,ia,2)
      if (oharm) then
         nn = newbas1-newbas0
         call expharm(v,'vectors',scr)
c...     move real vectors to the end
         do i=l0,1,-1
            e(i+nn) = e(i)
            call dcopy(ndim,v(1,i),1,v(1,i+nn),1)
         end do
         call vclr(v,1,ndim*nn)
         call vclr(e,1,nn)
      end if
c
      if (nprint .ne. 5) go to 120
      write (iwr,9028)
      call prev(v,e,l1,l1,ndim)
  120 continue
c
c     ----- eliminate eigenvectors for which eigenvalue is less
c           than tol=tol0
c
      if (o_depen_print) then
         write(iwr,9052) (i,e(newbas1-newbas0+i),i=1,newbas0)
      end if 
      dum = e(l1-newbas0+1)
      j = 0
      k = 0
      kk=0
      do 180 i = 1,l1
      if(e(i) .lt. tol1) then
      kk=kk+1
      endif
      if (e(i) .lt. tol0 .or. 
     1    i.le. (newbas1-newbas0+n_depen_check) ) go to 160
      j = j+1
      e(j) = done/  dsqrt(e(i))
      if (i.ne.j) call dcopy(l1,v(1,i),1,v(1,j),1)
      go to 180
  160 k = k+1
  180 continue
      l0 = l1-k
      if(kk.ne.l1-newbas0) write(iwr,9047) dum,kk-(l1-newbas0),tol1
      if (l0.ne.newbas0) then
         write (iwr,9048) dum,k-(l1-newbas0),tol0,n_depen_check,l0
cjvl     if(l0.ne.newbas0)call caserr2('linear dependance detected')
         newbas0 = l0
         odepen = .true.
c...  get symmetry numbers right
         do i=1,8
            nsymh(i) = nsym0(i)
            nsym0(i) = 0
         end do
         nav = lenwrd()
         call secget(isect(490),51,iblk51)
         call readi(mmmm,mach(13)*nav,iblk51,idaf)
         do i=1,l0
            ibig=idamax(l1,v(1,i),1)
            ii=isymaos(ibig)
            nsym0(ii) = nsym0(ii) + 1 
         end do
      end if
c     if (zruntp.ne.'force') n_depen_check = 0
c
c...  because l0 may differ from l1 generate a second ilif in ilifq0
c
      do i=1,l0
         ilifq0(i)=(i-1)*l0
      end do
c
c     ----- form canonical orthonormal orbitals -----
c
       do 200 j=1,l0
      call dscal(l1,e(j),v(1,j),1)
 200   continue
      scr(1)=dble(l1)
      call dcopy(l1,e,1,scr(2),1)
c     clear "non-existent" vectors
      if (l0.lt.l1) call vclr(v(1,l0+1),1,(l1-l0)*l1)
c
c     ----- write the canonical orbitals on ed7  -----
c
      call wrt3(scr,l1+1,iblkv,idaf)
      ibl3qs=iblkv+lenv
      call wrt3(v,l3,ibl3qs,idaf)
400   if (nprint .ne. 5) return
      write (iwr,9068)
      call prev(v,e,l0,l1,ndim)
      return
 9008 format(/,5x,14(1h-),/,5x,14hoverlap matrix,/,5x,14(1h-))
 9028 format(/,5x,17(1h-),/,5x,17heigenvectors of s,/,5x,17(1h-))
 9047 format(//
     b,' ############################################################',/
     2,' ###        possible linear dependence diagnosed          ###',/
     3,' ### Smallest eigenvalue of overlap matrix is', 1pe12.4,' ###',/
     4,' ### There are',i6,' eigenvalue(s) less than ', 1pe12.4,' ###',/
     c,' ############################################################')
 9048 format(//
     a,' ############################################################',/
     b,' ############################################################',/
     2,' ###        possible linear dependence diagnosed          ###',/
     3,' ### Smallest eigenvalue of overlap matrix is', 1pe12.4,' ###',/
     4,' ###',i5,' eigenvector(s) are eliminated (',1pe9.2,
     4 '/',i5,') ###',/
     5,' ### the number of canonical orbitals kept is',i6,'       ###',/
     6,' ###  Beware .... ONLY (not extensively) TESTED on SCF    ###',/
     c,' ############################################################',/
     d,' ############################################################'/)
 9052 format(//,' **the eigenvalues of the overlap matrix (numbered)**',
     1       /(2x,5(i5,1x,1pe10.3)))
 9068 format(/,5x,30(1h-),/,5x,30hcanonical orthonormal orbitals,/, 5x,
     +     30(1h-))
      end
      subroutine ortho1(q,s,v,t,ia,n,l0,l1,l2,l3,ndim)
c
c...  this routine in intimately linked with gamess (see rdedx)
c...  vectors in orthogonal basis return in v (complete orthogonal set)
c
c...  ia : iky (i*i-1)/1)
c...  n :  # vectors
c...  l0 : dimension of vectors in orthogonal basis
c...  l1 : dimension of vectors in non-orthogonal basis
c...  l2,l3 triangle/square dimensions corresponding to l1
c...  ndim first dimension of v and q arrays
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension q(ndim,*),s(*),v(ndim,*),t(*),ia(*)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer ibl7s, ibl7t, ibl7f, ibl7st,ibl7tt,ibl7ft,ibl3qs
      integer ibl7qa,ibl7pa,ibl7ea,ibl7qb,ibl7pb,ibl7eb
      integer ibl7x, ibl7y, ibl7z
      integer ibl7so,ibl7o, ibl7la,ibl7ec,ibl7ha
      common /scra7/ ibl7s, ibl7t, ibl7f, ibl7st,ibl7tt,ibl7ft,ibl3qs,
     +               ibl7qa,ibl7pa,ibl7ea,ibl7qb,ibl7pb,ibl7eb,
     +               ibl7x, ibl7y, ibl7z,
     +               ibl7so,ibl7o(8), ibl7la, ibl7ec, ibl7ha
c
c
      call rdedx(s,l2,ibl7st,num8)
      do 160 j = 1,n
      t(1)=s(1)*v(1,j)
      do 140 i=2,l1
      m=ia(i)
      call daxpy(i-1,v(i,j),s(m+1),1,t,1)
      t(i)=ddot(i,s(m+1),1,v(1,j),1)
 140  continue
      call dcopy(l1,t(1),1,v(1,j),1)
 160  continue
      call rdedx(q,l3,ibl3qs,idaf)
      do 200 j = 1,n
      do 180 i = 1,l0
      t(i) = ddot(l1,q(1,i),1,v(1,j),1)
  180 continue
      call dcopy(l0,t(1),1,v(1,j),1)
 200  continue
c
      call schmd(v,n,l0,ndim)
      return
      end
      subroutine schmd(v,m,n,ndim)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      dimension v(ndim,*)
      data done,tol /1.0d0,1.0d-10/
      if (m .eq. 0) go to 180
c
c     ----- orthonormalize first -m- mo's -----
c
      do 160 i = 1,m
c
      dumi=done/dnrm2(n,v(1,i),1)
      call dscal(n,dumi,v(1,i),1)
c
      if (i .eq. m) go to 160
      i1 = i+1
      do 140 j = i1,m
      dum = -ddot(n,v(1,j),1,v(1,i),1)
      call daxpy(n,dum,v(1,i),1,v(1,j),1)
  140 continue
c
  160 continue
      if (m .eq. n) return
  180 continue
c
c     ----- get orthogonal space -----
c
      i = m
      j = 0
  200 i0 = i
      i = i+1
      if (i .gt. n) return
  220 j = j+1
      if (j .gt. n) go to 320
      call vclr(v(1,i),1,n)
      v(j,i) = done
c
      do 260 ii = 1,i0
      dum = -ddot(n,v(1,ii),1,v(1,i),1)
      call daxpy(n,dum,v(1,ii),1,v(1,i),1)
  260 continue
      dumi = ddot(n,v(1,i),1,v(1,i),1)
      if (  dabs(dumi) .lt. tol) go to 220
      dumi = done/ dsqrt(dumi)
      call dscal(n,dumi,v(1,i),1)
      go to 200
  320 write (iwr,9008) i0
      call caserr2('redundant vectors detected')
 9008 format(33h redundant set of vectors. stop. , i5,
     +     26h independant vectors found)
      return
      end
      subroutine  setz(v,d,e,l1,l2,l3)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension d(*),v(*),e(*)
      call vclr(e,1,l1)
      call vclr(d,1,l2)
      call vclr(v,1,l3)
      return
      end
      function tracep(a,b,n)
c
c     ----- trace of product of 2 symmetric matrices
c           -a- and -b- stored linearly -----
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension a(*),b(*)
      ltri=n*(n+1)/2
      t=ddot(ltri,a,1,b,1)
      t=t+t
      k=1
      do 1 i=1,n
      t=t-a(k)*b(k)
1     k=k+i+1
      tracep=t
      return
      end
      subroutine pusql(v,m,n,ndim)
c
c     ----- punch out a square matrix with ordering labels -----
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      dimension v(ndim,*)
      do 120 j = 1,m
      ic = 0
      max = 0
  100 min = max+1
      max = max+5
      ic = ic+1
      if (max .gt. n) max = n
      write (ipu,9008) j,ic,(v(i,j),i = min,max)
      if (max .lt. n) go to 100
  120 continue
      return
 9008 format(i2,i3,5e15.8)
      end
c
c-----------------------------------------------------------------------
c
      subroutine prgeom(coord,chg,ztag,natoms,iwr)
      implicit none
      integer natoms, iwr
c     map80 maps the position of charge in the symz array (passed in as chg)
c     to those in the other atomic information arrays
      integer map80(natoms)
      real*8 coord(3,natoms)
      real*8 chg(natoms)
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
c
      integer ERR_NO_CODE
      parameter(ERR_NO_CODE=0)
c
c
c ==== if the error will occur on all nodes, we can ====
c      handle the output more cleanly
c
      integer ERR_SYNC, ERR_ASYNC
      parameter(ERR_ASYNC=0)
      parameter(ERR_SYNC=301)
c
c ==== legal error classes - these control an explanatory  ====
c      message (one per class, see machscf.m).
c
      integer ERR_NO_CLASS, ERR_INCOMPREHENSIBLE,
     &  ERR_INCOMPATIBLE, ERR_DIMENSION, ERR_FIXED_DIMENSION,
     &  ERR_UNIMPLEMENTED, ERR_INTERNAL,  ERR_USER_DIMENSION,
     &  ERR_SWITCHED_OUT, ERR_UNLUCKY

      parameter(ERR_NO_CLASS=0)
c - fix input errors and try again
      parameter(ERR_INCOMPREHENSIBLE=101)
c - incompatible 
      parameter(ERR_INCOMPATIBLE=102)
c - code dimensions exceeded - redimension code
      parameter(ERR_DIMENSION=103)
c - code dimensions exceeded - redimension or contact support
      parameter(ERR_FIXED_DIMENSION=104)
c - request unimplemented option
      parameter(ERR_UNIMPLEMENTED=105)
c - internal error - please report to support
      parameter(ERR_INTERNAL=106)

c - input dimension exceeded - change allocation
c   in input file and re-run
      parameter(ERR_USER_DIMENSION=107)

c - option disabled at configure stage
      parameter(ERR_SWITCHED_OUT=108)

c - some algorithmic or case-specififc problem
      parameter(ERR_UNLUCKY=109)

c
c  System message flag
c
      integer ERR_NO_SYS, ERR_SYS
      parameter(ERR_NO_SYS=0)
      parameter(ERR_SYS=201)
      character*8 ztag(natoms)
      integer nbqnop
c
c     ----- print the geometry of the molecule.
c
c     The format is chosen such that the printed result can be used
c     directly in an input file without requiring editing.
c
      integer i, j, l, isymz
c
      write(iwr,6000)
      nbqnop = 0
      do i=1,natoms
         if (.not. oprint(31) .or. (ztag(i)(1:2) .ne. 'bq'))then
            write(iwr,6010)(coord(j,i),j=1,3),chg(i),ztag(i)
         else
            nbqnop = nbqnop + 1
         endif
      enddo
      if (nbqnop .gt. 0)then
         write (iwr,6011) nbqnop
      endif
      write(iwr,6020)
c
 6000 format(9x,'x',14x,'y',14x,'z',12x,'chg',2x,'tag'/2x,60('='))
 6010 format(2x,3f15.7,1x,f7.2,2x,a8)
 6011 format(2x,'Output of ',i5,' BQ centres suppressed')
 6020 format(2x,60('='))
      end
c
c-----------------------------------------------------------------------
c
      subroutine prev(v,e,m,n,ndim)
c
c     ----- print out e and v-matrices
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
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
      dimension v(ndim,*),e(*)
      max = 10
      if (oprint(20)) max=7
      imax = 0
  100 imin = imax+1
      imax = imax+max
      if (imax .gt. m) imax = m
      write (iwr,9008)
      if(oprint(20)) go to 130
      write (iwr,8068) (e(i),i = imin,imax)
      write (iwr,9008)
      write (iwr,8028) (i,i = imin,imax)
      write (iwr,9008)
      do 150 j = 1,n
  150 write (iwr,8048) j,zbflab(j),(v(j,i),i = imin,imax)
      go to 140
130   write (iwr,9068) (e(i),i = imin,imax)
      write (iwr,9008)
      write (iwr,9028) (i,i = imin,imax)
      write (iwr,9008)
      do 120 j = 1,n
  120 write (iwr,9048) j,zbflab(j),(v(j,i),i = imin,imax)
 140  if (imax .lt. m) go to 100
      return
 9008 format(/)
 8028 format(17x,10(3x,i3,3x))
 8048 format(i5,2x,a10,10f9.4)
 8068 format(17x,10f9.4)
 9028 format(17x,12(6x,i3,6x))
 9048 format(i5,2x,a10,7f15.10)
 9068 format(17x,7f15.10)
      end
c
c-----------------------------------------------------------------------
c
      subroutine prevg(v,e,m,n,ndim)
c
c     ----- Print out the eigenvalues and eigenvectors of the g-matrix
c           alongside the variable labels from the z-matrix (if any).
c
c           This routine is used in formbgxz if the determinant of
c           the g-matrix is small to hopefully assist the user in 
c           redefining the z-matrix when problems with the g-matrix
c           inversion are suspected.
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
c
      integer ianz, iz, lbl, lalpha, lbeta, nz, nvar
      real*8 bl, alpha, beta
      common /czmat/ ianz(maxnz),iz(maxnz,4),bl(maxnz),
     +               alpha(maxnz),beta(maxnz),lbl(maxnz),
     +               lalpha(maxnz),lbeta(maxnz),nz,nvar
c
      character *8 zvar
      common/csubch/zvar(maxvar)
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
      dimension v(ndim,*),e(*)
      max = 10
      imax = 0
  100 imin = imax+1
      imax = imax+max
      if (imax .gt. m) imax = m
      write (iwr,9008)
      write (iwr,8068) (e(i),i = imin,imax)
      write (iwr,9008)
      write (iwr,8028) (i,i = imin,imax)
      write (iwr,9008)
      if (nz.ge.2) then
         do 20 i = 2 , nz
            ibl = iabs(lbl(i))
            if (ibl.ne.0) then
              write (iwr,8048) i-1,zvar(ibl),(v(i-1,k),k = imin,imax)
            else
              write (iwr,8049) i-1,(v(i-1,k),k = imin,imax)
            end if
 20      continue
c
         if (nz.ge.3) then
            j = nz - 3
            do 30 i = 3 , nz
               ialpha = iabs(lalpha(i))
               if (ialpha.ne.0) then
                 write (iwr,8048) i+j,zvar(ialpha),(v(i+j,k),
     +                                            k = imin,imax)
               else
                 write (iwr,8049) i+j,(v(i+j,k),k = imin,imax)
               end if
 30         continue
c
            if (nz.ge.4) then
               j = nz + nz - 6
               do 40 i = 4 , nz
                  ibeta = iabs(lbeta(i))
                  if (ibeta.ne.0) then
                    write (iwr,8048) i+j,zvar(ibeta),(v(i+j,k),
     +                                              k = imin,imax)
                  else
                    write (iwr,8049) i+j,(v(i+j,k), k = imin,imax)
                  end if
 40            continue
            end if
         end if
      end if
      if (imax .lt. m) go to 100
      return
 9008 format(/)
 8028 format(17x,10(3x,i3,3x))
 8048 format(i5,2x,a10,10f9.4)
 8049 format(i5,12x,10f9.4)
 8068 format(17x,10f9.4)
      end
c
c-----------------------------------------------------------------------
c
      subroutine prtril(d,n)
c
c     ----- print out a triangular matrix with labels
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
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
      dimension d(*),dd(10)
      mmax = 10
      if (oprint(20)) mmax=7
      imax = 0
  100 imin = imax+1
      imax = imax+mmax
      if (imax .gt. n) imax = n
      write (iwr,9008)
      if(.not.oprint(20)) write (iwr,8028) (zbflab(i),i = imin,imax)
      if( oprint(20)) write (iwr,9028) (zbflab(i),i = imin,imax)
      write (iwr,9008)
      do 160 j = 1,n
      k = 0
      do 140 i = imin,imax
      k = k+1
      m =max(i,j)*(max(i,j)-1)/2+min(i,j)
  140 dd(k) = d(m)
      if( oprint(20)) write (iwr,9048) j,zbflab(j),(dd(i),i = 1,k)
      if(.not.oprint(20)) write (iwr,8048) j,zbflab(j),(dd(i),i = 1,k)
  160 continue
      if (imax .lt. n) go to 100
      return
 9008 format(/)
 9028 format(17x,7(3x,a10,2x))
 9048 format(i5,2x,a10,7f15.10)
 8028 format(17x,10(1x,a8))
 8048 format(i5,2x,a10,10f9.4)
      end
      subroutine prsql(v,m,n,ndim)
c
c     ----- print out a square matrix with labels
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
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
      dimension v(ndim,*)
      max = 10
      if (oprint(20)) max=7
      imax = 0
  100 imin = imax+1
      imax = imax+max
      if (imax .gt. m) imax = m
      write (iwr,9008)
      if(oprint(20)) write (iwr,9028) (i,i = imin,imax)
      if(.not.oprint(20)) write (iwr,8028) (i,i = imin,imax)
      write (iwr,9008)
      do 120 j = 1,n
      if(oprint(20)) write (iwr,9048) j,zbflab(j),(v(j,i),i = imin,imax)
      if(.not.oprint(20)) write (iwr,8048) j,zbflab(j),(v(j,i),i =
     * imin, imax)
 120   continue
      if (imax .lt. m) go to 100
      return
 9008 format(1x)
 9028 format(17x,7(6x,i3,6x))
 9048 format(i5,2x,a10,7f15.10)
 8028 format(17x,10(3x,i3,3x))
 8048 format(i5,2x,a10,10f9.4)
      end
      subroutine getq(q,eig,pop,nbas,newb,mode,ieig,ipop,numd,zname)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      parameter (mxorb1=maxorb+1)
      dimension q(*),eig(*),pop(*)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical odiis, optester, odynamic, omaxcyc,onoor,otestd
      integer maxcyc,mconv,nconv,npunch,icoupl,ifolow,irotb
      integer iter,kcount,iextin,iterv,idiisf
      real*8 accdi1,accdi2,dmpcut,acurcy,en,etot
      real*8 ehf,ehf0,diff,rshift,exttol,dmptol,vshtol
      real*8 damp,damp0,diffd,diffp,de,deavg,diffsp
      real*8 ek, vir,diffpp
      common/scfopt/maxcyc,mconv,nconv,npunch,accdi1,accdi2,odiis,
     +      icoupl,ifolow,irotb,dmpcut,acurcy,en,etot,ehf,ehf0,diff,
     +      iter,kcount,rshift,exttol,dmptol,vshtol,iextin,
     +      iterv,damp,damp0,diffd,diffp,diffpp,de,deavg,diffsp,
     +      ek,vir,idiisf,optester,odynamic,omaxcyc,onoor,otestd
c
      integer numdu, iblkdu, maxb, iblkla
      logical orevis
      real*8 apos
      common/sector/numdu,iblkdu,orevis(2),
     + apos(508),maxb,iblkla
      common/junkc/zcom(19),ztit(10)
c      common/junk/ilifd(maxorb),ntrad(maxorb),itrad(mxorb3),
c     *ctrad(mxorb3)
      common/blkqio/deig(maxorb),dpop(mxorb1),
     *nba,new,ncol,jeig,jpop
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
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
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
      data m3,m29/3,29/
      if(numd.gt.508)call caserr2(
     *'invalid section specified for vector input')

      if(nprint.ne.-5)write(iwr,100)zname,numd,iblkdu,yed(numdu)
100   format(/a7,'vectors restored from section',i4,
     *' of dumpfile starting at block',i6,' of ',a4)
      call secget(numd,m3,k)
      call rdchr(zcom(1),m29,k,numdu)
      call reads(deig,mach(8),numdu)
      k=k+1+lensec(mach(8))
      if(nprint.ne.-5)write(iwr,200)(zcom(7-i),i=1,6),ztit
c
      if(new.ne.ncol) then
c           (**harmonic**)
        ncol = new
      end if
c
      if(new.eq.ncol)go to 300
 310  call caserr2(
     *'vectors restored from dumpfile have incorrect format')
 300  nbas=nba
      newb=new
      ieig=jeig
      ipop=jpop
      call dcopy(ncol,dpop,1,pop,1)
      call dcopy(ncol,deig,1,eig,1)
      etot=dpop(mxorb1)
200   format(/' header block information :'/
     *' vectors created under acct. ',a8/1x,
     *a7,'vectors created by ',a8,'  program at ',
     *a8,' on ',a8,' in the job ',a8/
     *' with the title: ',10a8)
      if(mode.lt.3)go to 101
      if(num.lt.nbas.or.num.lt.newb)go to 310
      go to 500
 101  if(num.ne.nbas.or.num.ne.newb)go to 310
      if(mode.le.2)call ctrchk(k)
500   nw=nbas*newb
      k=k+lensec(mach(9))
      call rdedx(q,nw,k,numdu)
      return
      end
      subroutine putdev(q,iseco,ipos,ioff)
c
c      retrieve and store gamess q, eval and d
c
      implicit real*8  (a-h,o-z)
      character *8 grhf
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
      parameter (mxorb1=maxorb+1)
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
      integer master, indxi, indxj, nfiles, junits, jblkrs, jblkas
      integer nfilef, junitf, jblkrf, jblkaf, isecbl, iscftp
      integer lword4, ilow4, ncol4, nsa4, newb4, nbas4, lenb4
      integer ndump4, iblkq4, nblkq4, lena4, nbb4, ionsv4, isecv4
      logical oprin4, oindx4
      integer npas41, npas42, iacc4
c
      common /cndx40/ master,indxi,indxj,nfiles,junits,jblkrs,jblkas,
     + nfilef,junitf,jblkrf,jblkaf,isecbl,iscftp,
     + lword4,ilow4,ncol4,nsa4,newb4,nbas4,lenb4,ndump4,iblkq4,
     + nblkq4,lena4,nbb4,
     + oprin4(10),oindx4,ionsv4,npas41,npas42,iacc4,isecv4
c
       integer len_cndx40
       parameter (len_cndx40=37)
c      used: restre(util1),revise(util1),utyp21(server)
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
c
      character *8 com,dtitle
      common/junkc/com(19),dtitle(10)
      logical iftran
      common/small/
     + ilifcs(maxorb),ntrans(maxorb),itrans(mxorb3),ctrans(mxorb3),
     + iftran,isp,
     + value(maxorb),occ(mxorb1),
     + nbasis,newbas,ncol,ivalue,iocc,ift
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      integer ltask, lcoord, iplot, iprop, imull, idma, inbo
      logical opass1, opass2, opass3, orege, opass4, offsym
      logical opass5, opass6, opass8, opass9, odisc, oirc
      logical orgall, orestj, opass10, opass11
      logical omrdci, ofulci, omech, omopac
      logical orpa, omclr, odirpa, occsd, occsdt, oqcisd, oqcisdt
      logical omrdcid,opark,oclunp,ofill,odiesel
c
      common /restrj/ ltask,lcoord,opass1,opass2,opass3,orege,iplot,
     + iprop,imull,opass4,offsym,opass5,opass6,opass10,idma,inbo,
     + omrdci,ofulci,omech,omopac,opass8,opass9,odisc,
     + oirc,orpa,omclr,odirpa,occsd,occsdt,oqcisd,oqcisdt,orgall,
     + omrdcid,orestj(7),opark,oclunp,ofill(6),odiesel, opass11
c
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
      logical lfield,fixed,lex,ldam12,ldam13,ldam23,ldiis
      common/scfblk/en,etotal,ehf,sh1(2),sh2(2),gap1(2),gap2(2),
     &              d12,d13,d23,canna,cannb,cannc,fx,fy,fz,
     &              lfield,fixed,lex,ldam12,ldam13,ldam23,ldiis,
     &              ncyc,ischm,lock,maxit,nconv,npunch
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
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
c
      real*8 ctran
      integer ilifc, ntran, itran
      logical otran,otri
      common/tran/ilifc(maxorb),ntran(maxorb),itran(mxorb3),
     +            ctran(mxorb3),otran,otri
c
      dimension q(*),ioffs(2)
c
c allocate storage
c
      data m1,m19/1,19/
      data grhf/'grhf'/
c
      if (iacsct(isect(ipos  )).lt.0. or. 
     +    iacsct(isect(ipos+1)).lt.0. or. 
     +    iacsct(isect(ipos+2)).lt.0 ) return
c
      l3 = num*num

      if (scftyp.ne.grhf) then
         ilen = l3 + 2*nx
      else
         ilen = 2*l3 + nx + num
      end if

      i10 = igmem_alloc(ilen)
      i20 = i10 + l3
      i30 = i20 + nx
      if (scftyp.ne.grhf) then
c        last = i30 + nx
      else
         i40 = i30 + l3
c        last = i40 + num
      end if

      len = lensec(num)
      lenx = lensec(nx)
c
c     retrieve and store density matrix
c
      call secget(isect(497),m19,iblok)
      ioffs(1) = iblok
      ioffs(2) = iblok + len + lenx
      call rdedx(q(i20),nx,ioffs(ioff),ifild)
      if (scftyp.eq.grhf) then
         call rdedx(q(i30),nx,ioffs(ioff+1),ifild)
         call vadd(q(i20),1,q(i30),1,q(i20),1,nx)
      end if
c     store
      call secput(isect(ipos),ipos,lenx,iblkdc)
      lds(isect(ipos)) = nx
      call wrt3(q(i20),nx,iblkdc,ifild)
c
c     lagrangian for grhf cases
c
      if (scftyp.eq.grhf) then
         call putlag(q(i20),q(i10),q(i30),q(i40),num)
      end if
c
c     store vectors ( including a 'ctrans' block')
c
      nav = lenwrd()
      nprin = nprint
      nprint = -5
      call getq(q(i10),value,occ,nbasis,newbas,m1,m1,m1,iseco,scftyp)
      if (.not.otran) call tdown(q(i10),ilifq,q(i10),ilifq,num)
      nprint = nprin
c
      lenw = num*ncoorb
      lds(isect(ipos+1)) = lenw
      len = lensec(lenw) + mvadd
      call secput(isect(ipos+1),ipos+1,len,iblkdc)
      iftran = .true.
      ncol = ncoorb
      iocc = ncoorb
      ivalue = ncoorb
      m29 = 29
      call wrtc(com,m29,iblkdc,ifild)
      call wrt3s(value,mach(8),ifild)
      call wrt3is(ilifcs(1),mach(9)*nav,ifild)
      call wrt3s(q(i10),lenw,ifild)
c
c     store eigenvalues
c
      len = lensec(num)
      lds(isect(ipos+2)) = num
      call secput(isect(ipos+2),ipos+2,len,iblkdc)
      call wrt3(value,num,iblkdc,ifild)
c
      call gmem_free(i10)
      return
      end
      subroutine putgrd(eg)
      implicit real*8  (a-h,o-z)
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
      dimension eg(*)
c
c ----- store gradient on isect(14)
c
      if (iacsct(isect(14)). lt. 0) return
      nblock = lensec(maxat*3)
      lds(isect(14)) = maxat*3
      call secput(isect(14),14,nblock,iblok)
      call wrt3(eg,lds(isect(14)),iblok,ifild)
      return
      end
      subroutine putlag(d,v,e,dd,ndim)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      integer ibl3s, ibl3t, ibl3f, ibl3x, ibl3y, ibl3z
      integer ibl3rs, ibl3qa, ibl3pa, ibl3ea, ibl3qb, ibl3pb
      integer ibl3eb, ibl3g, ibl3hs, ionsec, ibl3op, ions2
      integer isecqm
      common /dump3/ ibl3s, ibl3t, ibl3f, ibl3x, ibl3y, ibl3z,
     +               ibl3rs,ibl3qa,ibl3pa,ibl3ea,ibl3qb,ibl3pb,
     +               ibl3eb,ibl3g,ibl3hs,ionsec,ibl3op,ions2,
     +               isecqm
c
c
      integer igrad,ifout1,ifout2,ifmola,if2,if4,iflagr
      integer isecdm,isecda,iseclm,isecla
      integer idmmc,iblkmm,lblmm,ifiled,iwordd,iseclg
      common /dm/ igrad,ifout1,ifout2,ifmola,if2,if4,iflagr,
     +            isecdm,isecda,iseclm,isecla,
     +            idmmc,iblkmm,lblmm,ifiled,iwordd,iseclg
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
      real*8 cicoef, f, alpha, beta
      integer no, nco, nseto, npair, ncores, ibm, nset, nopen
      integer  nope, noe
      logical old
      common/scfwfn/cicoef(2,12),f(25),alpha(325),beta(325),
     + no(10),nco,nseto,npair,ncores,ibm,old,nset,nopen,nope,noe(10)
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
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
      real*8 ctran
      integer ilifc, ntran, itran
      logical otran,otri
      common/tran/ilifc(maxorb),ntran(maxorb),itran(mxorb3),
     +            ctran(mxorb3),otran,otri
c
      dimension v(ndim,*),d(*),e(*),dd(*)
      data m20/20/
      data pt5/0.5d0/
c
      if(iacsct(isect(42)).lt.0) return
      l3 = num*num
c
      norb = nco + npair + npair
      if (nseto.gt.0) then
         do 20 i = 1 , nseto
            norb = norb + no(i)
 20      continue
      end if
      call rdedx(v,l3,ibl3qa,idaf)
      if (.not.otran) call tdown(v,ilifq,v,ilifq,num)
      l3orb = norb*norb
      call secget(iseclg,m20,iblk20)
      call rdedx(e,l3orb,iblk20,idaf)
c
c     ----- zero out weighted density array -----
c
      call vclr(d,1,nx)
c
c     ----- calculate -tr(ce(ct)sa) -----
c
c     ----- note that e(kl) is used exactly twice. divide by
c           two to get the values appropriate for the generalized
c           lagrangian multipliers -----
c
      do 70 i = 1 , num
c
c     ---calculate the half transform first -----
c
         kl = 0
         do 40 l = 1 , norb
            dd(l) = 0.0d0
            do 30 k = 1 , norb
               kl = kl + 1
               dd(l) = dd(l) - v(i,k)*e(kl)
 30         continue
 40      continue
         call dscal(norb,pt5,dd,1)
         do 60 j = 1 , num
            ij = iky(max(i,j)) + min(i,j)
            do 50 l = 1 , norb
               d(ij) = d(ij) + dd(l)*v(j,l)
 50         continue
 60      continue
 70   continue
c
      len = lensec(nx)
      m = 42
      call secput(isect(42),m,len,iblk42)
      call wrt3(d,nx,iblk42,idaf)
      lds(isect(42)) = nx
c
      return
      end
      subroutine putstv(q)
      implicit real*8  (a-h,o-z)
      logical o1e
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
      integer ltask, lcoord, iplot, iprop, imull, idma, inbo
      logical opass1, opass2, opass3, orege, opass4, offsym
      logical opass5, opass6, opass8, opass9, odisc, oirc
      logical orgall, orestj, opass10, opass11
      logical omrdci, ofulci, omech, omopac
      logical orpa, omclr, odirpa, occsd, occsdt, oqcisd, oqcisdt
      logical omrdcid,opark,oclunp,ofill,odiesel
c
      common /restrj/ ltask,lcoord,opass1,opass2,opass3,orege,iplot,
     + iprop,imull,opass4,offsym,opass5,opass6,opass10,idma,inbo,
     + omrdci,ofulci,omech,omopac,opass8,opass9,odisc,
     + oirc,orpa,omclr,odirpa,occsd,occsdt,oqcisd,oqcisdt,orgall,
     + omrdcid,orestj(7),opark,oclunp,ofill(6),odiesel, opass11
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
      integer ibl3s, ibl3t, ibl3f, ibl3x, ibl3y, ibl3z
      integer ibl3rs, ibl3qa, ibl3pa, ibl3ea, ibl3qb, ibl3pb
      integer ibl3eb, ibl3g, ibl3hs, ionsec, ibl3op, ions2
      integer isecqm
      common /dump3/ ibl3s, ibl3t, ibl3f, ibl3x, ibl3y, ibl3z,
     +               ibl3rs,ibl3qa,ibl3pa,ibl3ea,ibl3qb,ibl3pb,
     +               ibl3eb,ibl3g,ibl3hs,ionsec,ibl3op,ions2,
     +               isecqm
c
      dimension o1e(6)
      dimension corev(512),array(6)
c
      dimension q(*)
c
      if (iacsct(isect(5)).lt.0.or.
     +    iacsct(isect(6)).lt.0 ) return
      l2 = num*(num+1)/2
      i10 = igmem_alloc(l2)
      i20 = igmem_alloc(l2)
      do loop =1,6
       o1e(loop) = .false.
      enddo
      o1e(1) = .true.
      o1e(3) = .true.
c
c ----- restore s,t,f from section 192 and store on dumpfile
c
      call getmat(q(i10),q(i10),q(i20),q(i10),q(i10),q(i10),
     +   array,num,o1e,ionsec)
c
c ----- store on isect(5) and isect(6)
c
      nblock = lensec(l2)
      lds(isect(5)) = l2
      call secput(isect(5),5,nblock,iblok)
      call wrt3(q(i10),l2,iblok,idaf)
c
      lds(isect(6)) = l2
      call secput(isect(6),6,nblock,iblok)
      call wrt3(q(i20),l2,iblok,idaf)
c
c ----- reset core
c
      call gmem_free(i20)
      call gmem_free(i10)
c
      return
      end
      subroutine copq_again(mpos)
c
c...  copy (vectors) section to  another section, to keep copies of orbitals
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      dimension buf(512)
      character*(*) test
      logical off,oexist
      save mpos_again,maxpos_again,int_again,int,off
      data off/.true./
c
      if (off) return
c
      int = int + 1 
c
      entry copq_again2(mpos)
      if (off) return
c
      if ((int/int_again)*int_again.eq.int.or.int.eq.1) then
         if (mpos_again.gt.maxpos_again) return
c
         call secinf(iblkdu,numdu,idum,idum)
         call secloc(mpos,oexist,ibli)
         if (.not.oexist) call caserr('no section in copq_again')
         call qsector(mpos,idum,iclass,ilen,'get')
         call secput(mpos_again,iclass,ilen,iblo)
         do i=1,ilen
            call search(ibli,numdu)
            call get(buf,nw,numdu)
            ibli = ibli + 1
            call search(iblo,numdu)
            call put(buf,nw,numdu)
            iblo = iblo + 1
         end do
         write(iwr,'(a,i4,a,i4)') ' wrote vectors(copq) for point ',int,
     1                            ' to section ',mpos_again
         mpos_again = mpos_again + 1
      end if
c
      return      
c
      entry set_copq_again(test)
c
      if (test.eq.'on') then
         off = .false.
      else if (test.eq.'off') then
         off = .true.
      else if (test.eq.'init') then
         off = .false.
         call inpi(int_again)
         if (int_again.eq.0) int_again = 5
         call inpi(mpos_again)
         if (mpos_again.eq.0) mpos_again = 100
         call inpi(maxpos_again)
         if (maxpos_again.eq.0) maxpos_again = 200
c
         write(iwr,'(a,i4,a,a,i4,a,i4,a)') 
     1      ' *** write vectors(copq) every ',int_again,' point',
     2      ' starting with section ',mpos_again,' till ',maxpos_again,
     3      ' *** '
      else
         call caserr('wrong set_putq_again')
      end if
c
      return
      end
      subroutine putq(zcomm,ztit,eig,def,norb,norbn,ncolu,ieig,
     *ideff,q,mpos,iblkq)
c... standard e.vector outputting routine(+ header blocks)
c... allow VB/servec weird dimensioning
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      dimension q(*),zcomm(*),ztit(*),eig(*),def(*)
      integer numdu, iblkdu, maxb, iblkla
      logical orevis
      real*8 apos
      common/sector/numdu,iblkdu,orevis(2),
     + apos(508),maxb,iblkla
      common/junkc/zcom(19),ztitle(10)
      common/blkqio/value(maxorb),pop(maxorb),escf,
     *nbasis,newbas,ncol,ivalue,ipop
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
c
      logical odiis, optester, odynamic, omaxcyc,onoor,otestd
      integer maxcyc,mconv,nconv,npunch,icoupl,ifolow,irotb
      integer iter,kcount,iextin,iterv,idiisf
      real*8 accdi1,accdi2,dmpcut,acurcy,en,etot
      real*8 ehf,ehf0,diff,rshift,exttol,dmptol,vshtol
      real*8 damp,damp0,diffd,diffp,de,deavg,diffsp
      real*8 ek, vir,diffpp
      common/scfopt/maxcyc,mconv,nconv,npunch,accdi1,accdi2,odiis,
     +      icoupl,ifolow,irotb,dmpcut,acurcy,en,etot,ehf,ehf0,diff,
     +      iter,kcount,rshift,exttol,dmptol,vshtol,iextin,
     +      iterv,damp,damp0,diffd,diffp,diffpp,de,deavg,diffsp,
     +      ek,vir,idiisf,optester,odynamic,omaxcyc,onoor,otestd
c
c
      real*8 ctran
      integer ilifc, ntran, itran
      logical otran,otri
      common/tran/ilifc(maxorb),ntran(maxorb),itran(mxorb3),
     +            ctran(mxorb3),otran,otri
c
      data m29/29/
       if(mpos.gt.508)call caserr2(
     *'invalid section specified for vector output')
      len1=lensec(mach(8))
      nbsq=norb*norb
      m3 = 3
c
      if (zcomm(5).eq.'vb'.or.zcomm(5).eq.'servec') then
         nbsq = norb*ncolu
         m3 = 33
      end if
c
      lenv=lensec(nbsq)
      len2=lensec(mach(9))
      j=len2+lenv+len1+1
      call secput(mpos,m3,j,iblk)
      iblkq=iblk+len2+len1+1
      escf=etot
      do 100 i=1,10
 100  ztitle(i)=ztit(i)
      do 200 i=1,19
 200  zcom(i)=zcomm(i)
      call dcopy(ncolu,eig,1,value,1)
      call dcopy(ncolu,def,1,pop,1)
      nbasis=norb
      newbas=norbn
      ncol=ncolu
      ivalue=ieig
      ipop=ideff
      call wrtc(zcom,m29,iblk,numdu)
c...    clear non-existent vectors abd populations and make e.v. large
      if(ncol.ne.norb.and.zcomm(5).ne.'vb'.and.zcomm(5).ne.'servec')then
         call vclr(q(ncol*norb+1),1,(norb-ncol)*norb)
         do i=ncol+1,norb
            value(i) = 9999900.0d0 + i*0.1d0
            pop(i) = 0.0d0
         end do
      end if
c
      call wrt3s(value(1),mach(8),numdu)
      nav = lenwrd()
      call wrt3is(ilifc(1),mach(9)*nav,numdu)
      call wrt3(q(1),nbsq,iblkq,numdu)
      call clredx
      call revind
      return
      end
      subroutine ctrchk(iblk)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
cfu      common/junk/ilifd(maxorb),ntrand(maxorb),itrand(mxorb3),
cfu     * ctrand(mxorb3),
cfu     *ilif(maxorb),ntr(maxorb),itr(mxorb3),ctr(mxorb3),of2
      common/blkctr/ilif(maxorb),ntr(maxorb),itr(mxorb3),ctr(mxorb3),of2
c
      real*8 ctran
      integer ilifc, ntran, itran
      logical otran,otri
      common/tran/ilifc(maxorb),ntran(maxorb),itran(mxorb3),
     +            ctran(mxorb3),otran,otri
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      integer numdu, iblkdu, maxb, iblkla
      logical orevis
      real*8 apos
      common/sector/numdu,iblkdu,orevis(2),
     + apos(508),maxb,iblkla
      equivalence (if1,otran),(if2,of2)
c
      nav = lenwrd()
      call readi(ilif,mach(9)*nav,iblk,numdu)
      if(if1.eq.if2)go to 2
      write(iwr,30)
 30   format(1x,' ***** inconsistent ctrans detected')
 3    call caserr2('retrieved eigenvectors are invalid')
 2    if(otran)return
      do 1 i=1,num
      ic=ilifc(i)
      n=ntran(i)
      if(ic.ne.ilif(i).or.n.ne.ntr(i)) then
       write(iwr,10) i, ilif(i), ntr(i), ic, n
 10   format(/' basis function ',i4,' ilif,ntr (input) = ',2i5/
     * 20x,' ilif,ntr (setup) = ',2i5)
       go to 3
      endif
      do 1 j=1,n
      if(itr(ic+j).ne.itran(ic+j)) then
       write(iwr,20) i,j, itr(ic+j), itran(ic+j)
  20   format(/' basis function ',i4, ' component ',i2,
     * ' input = ',i4,' setup = ',i4/)
      go to 3
      endif
 1    continue
      return
      end
      subroutine dmtx(d,v,p,ia,m,n,ndim)
      implicit real*8  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension ia(*),d(*),v(ndim,*),p(*)
      data dzero /0.0d0/
      call vclr(d,1,ia(n+1))
      if(m.eq.0) return
      ii=0
      do 120 i = 1,n
      do 130 k = 1,m
      dum = p(k)*v(i,k)
      if(dum.ne.dzero) then
      call daxpy(i,dum,v(1,k),1,d(ii+1),1)
         endif
  130 continue
  120 ii=ii+i
      return
      end
      subroutine dmtxp(d,v,p,ia,m,n,ndim)
      implicit real*8  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension ia(*),d(*),v(ndim,*),p(*)
      data dzero /0.0d0/
      ltri = ia(n+1)
      call vclr(d,1,ltri)
      if(m.ne.0) then
         idum = iipsci()
         ii=0
         do 120 i = 1,n
            if(.not. oipsci())then
               do 130 k = 1,m
                  dum = p(k)*v(i,k)
                  if(dum.ne.dzero) then
                     call daxpy(i,dum,v(1,k),1,d(ii+1),1)
                  endif
 130           continue
            endif
            ii=ii+i
 120     continue
         call pg_dgop(1771,d,ltri,'+')
      endif
      return
      end

CMR   renamed function enuc for debugger namespace...
      function enucf(n,cz,c)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      logical indij,indji
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
c---------
c  omtchg - omit 1e- integrals involving the specified charges
c
c   integrals of the form  < bf(a) | q(b) | bf (a) > are deleted
c
c   where a and b are centres indicated as follows
c
c  ztagcc - list of centre names a
c  ztagcq - for each a in ztagcc, the list of centres b
c
c  omtslf - omit self-energy of point charge array (applied to all
c           atoms with centre label of form bq* )
c      CHECK
c      parameter (mxombq=20,mxomat=50)
c---------

c
c obqf  - special treatment of forces .. skip force contributions
c         on non-bq centres
c         for use in qm/mm optimisations of surrounding
c         assumption is that bq's in this case don't have basis
c         functions or ecps
c
c omitted charge specifications
c

      integer mxombq, mxomat
      parameter (mxombq=50,mxomat=20)
      character*8 ztagcc, ztagcq
      logical omtchg, omtslf, obqf
      common /chgcc/ ztagcc(mxomat),ztagcq(mxomat,mxombq),omtchg,
     &     omtslf,obqf
c
c core-core repulsion parameters (see integb_lib for input)
c
      real*8 eta, d
      integer ichgat, indx, jndx
      logical indi
      integer mxpairs
      parameter(mxpairs=10)
      common/repp/eta(mxpairs),d(mxpairs),ichgat(750),npairs,
     +            indx(mxpairs),jndx(mxpairs),indi
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      character*8 field
      character*24 version
      common /drfopt/ field,version
      logical oreact
      integer intdrf
      common /drfopti/ oreact,intdrf
c
cafc
c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     this file contains the parameter settings for piet van duijnen's
c     reactionfield programme together with the sizes needed for the
c     v. duijnen - gamess interface as developed by 
c     nederkoorn and v. lenthe
c
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c     ==================================================================
c
c     parameter setting for reaction field programme van duijnen:
c
c=======================================================================
c
c     
      parameter (mxat = 10000)
      parameter (mxpol = 10000)
      parameter (mxpts = 10000)
      parameter (mxspe = 50)
      parameter (mxamb = 10)
      parameter (mxgrp = 500)
      parameter (mxnum = 1024)
      parameter (mxsh = 512)
      parameter (mxgs = 2048)
      parameter (mxnpts = 20)
      parameter (mxneq = 255)
      parameter (mcwrkmx = 10)
c
c ==================================================================
c
c     parameter setting new drf routines (drf/dimpar):
c
c===================================================================
c
      parameter (mxst = 5)
      parameter (mxgran = 10)
cmw      parameter (mxnpts = 50)
      parameter (mxgrpar = (mxgran*(mxgran+1))/2)
c
c     ==================================================================
c
c     parameter setting for v. duijnen - gamess interface:
c
c=======================================================================
c
      parameter (ncatom=500)
c
c  ncatom is maximal allowed classical atoms in zmat
c  total allowed atoms in zmat = maxnz + ncatom
c
c  maxrfa is maximal number of classical atoms in zmat + coordinates
c  total maximal atoms allowed is maxrfa + maxat
c
      parameter (maxrfa=mxpts)
c
c  maxamb is maximal number of ambigous atoms. it is set equal
c  to maximal z-cards possible in gamess input , i.e. 50
c
      parameter (maxamb=mxamb)
c
c=============================================
c for the pes-routines
c=============================================
c
      parameter (mxgr1=128, mxpes=100)
c
c=============================================
c for the drfexpc routines
c=============================================
c
      parameter (maxex=1000)
c
c=============================================
c for the drfsurf routine
c=============================================
c
      parameter (maxinv=10,maxnrc=10000)
c      parameter (maxnrc=10000)
c
c=============================================
c for the analrf routine
c=============================================
c
c
      integer ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep
      integer neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga
      integer idipcal,iclinte,iclintd,ieffpol,isodis,iclintr
      integer iextdip,itsolv,itermax,isolsav,imomsav,idisadd
      integer ianal,ineqex,ithole,irevdis,ihbond,igrppol
      integer maxci_drf
      common/drfpar1/ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep,
     +               neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga,
     +               idipcal,iclinte,iclintd,ieffpol,isodis,iclintr,
     +               iextdip,itsolv,itermax,isolsav,imomsav,idisadd,
     +               ianal,ineqex,ithole,irevdis,ihbond,igrppol,
     +               maxci_drf
c
      integer nxtpts, npol, nexp, natint, nshint, namb, nspec, ngran
      common/drfpar2/nxtpts,npol,nexp,natint,nshint,namb,nspec,ngran
c
      integer ndima, ndimb, nbem, ndim, nexp4, npol3, nwtr, nwtc, nzfa
      integer lomga, nomga, nchd, ngrnam, ngrnam2, nzfp, nzfn
      integer nneq, nneqrf, maxneq, neqdim, nqdim, nqcls
      common/drfpar3/ndima,ndimb,nbem,ndim,nexp4,npol3,nwtr,nwtc,nzfa,
     +               lomga,nomga,nchd,ngrnam,ngrnam2,nzfp,nzfn,
     +               nneq,nneqrf,maxneq,neqdim,nqdim,nqcls
c
      integer ifldin, ifldout, idrfout, modxza, irepopt
      integer iexpza, icmexp, iadexp, nodpe, igetden, iarfcal, iqmclr
      common/drfpar4/ifldin,ifldout,idrfout,modxza,irepopt,
     +               iexpza,icmexp,iadexp,nodpe,igetden,iarfcal,
     +               iqmclr
c
      real*8 gamdrf, dstmin, dstmax, hbondl, hbondr, afact, rfact
      real*8 cvgrel, agrpe, agrpm, agrpc, scffact, acur
      real*8 conci_drf
      common/drfpar5/gamdrf,dstmin,dstmax,hbondl,hbondr,afact,rfact,
     +               cvgrel,agrpe,agrpm,agrpc,scffact,acur,conci_drf
c
      integer iind, isur1, isur2, ilwt, ilvr, illur, ilindx
      common/drfpar6/iind,isur1,isur2,ilwt,ilvr,illur,ilindx
c
      real*8 unucrep, eoneel, ekin, enua, etwoel, extel, extnuc
      real*8 snucnuc, selel, snua, snucel, selnuc, stwoel, sextel
      real*8 sextnuc, selext, snucext, upolqm 
      real*8 uneqnuc, uneqel, uclas, suclas, uclase, uclasd, uclasr
      real*8 upolcl, udisp, rdisp, repmod, ustanuc, ustael
      common /rfene1/ unucrep,eoneel(mxst),ekin(mxst),enua(mxst),
     +  etwoel(mxst),extel(mxst),extnuc(mxst),snucnuc(mxst),
     +  selel(mxst),snua(mxst),snucel(mxst),selnuc(mxst),stwoel(mxst),
     +  sextel(mxst),
     +  sextnuc(mxst),selext(mxst),snucext(mxst),upolqm(mxst),
     +  uneqnuc(mxst),uneqel(mxst),uclas,suclas,uclase,uclasd,uclasr,
     +  upolcl,udisp(mxst),rdisp(mxst),repmod(mxst),
     +  ustanuc(mxst),ustael(mxst)
c
      real*8 smolnuc, smolel, smolmol, snucmol, selmol, sextmol
      real*8 smolext
      real*8 stotnuc, stotel, stotext, stotmol, upoleq, ucstst, ucstpl
      common/rfene2/smolnuc(mxst),smolel(mxst),smolmol(mxst),
     +  snucmol(mxst),selmol(mxst),sextmol(mxst),smolext(mxst),
     +  stotnuc(mxst),stotel(mxst),stotext(mxst),stotmol(mxst),
     +  upoleq(mxst),ucstst(mxst),ucstpl(mxst)
c
      real*8 uscf, uqm, uelst, uint, suqm, suint, uneq, uneqqm, uneqcl
      real*8 upolneq, usta, ustaqm, ustacl, ustaneq, uens
      common /rfene3/ uscf(mxst),uqm(mxst),uelst(mxst),uint(mxst),
     +  suqm(mxst),suint(mxst),uneq(mxst),uneqqm(mxst),uneqcl,
     +  upolneq(mxst),usta(mxst),ustaqm(mxst),ustacl,ustaneq(mxst),
     +  uens(mxst)
c
      real*8 snucno, selelo, snuao, snucelo, selnuco, stwoelo, sextelo
      real*8 sextno, selexto, snucexo, upolqmo, suclaso, upolclo
      common/rfene4/snucno(mxst),selelo(mxst),snuao(mxst),snucelo(mxst),
     +  selnuco(mxst),stwoelo(mxst),sextelo(mxst),sextno(mxst),
     +  selexto(mxst),snucexo(mxst),upolqmo(mxst),suclaso,upolclo
c
      real*8 smolno, smolelo, smolmo, snucmo, selmolo, sextmo, smolexo
      real*8 stotno, stotelo, stotexo, stotmo, upoleqo, ucstplo
      common/rfene5/smolno(mxst),smolelo(mxst),smolmo(mxst),
     +  snucmo(mxst),selmolo(mxst),sextmo(mxst),smolexo(mxst),
     +  stotno(mxst),stotelo(mxst),stotexo(mxst),stotmo(mxst),
     +  upoleqo(mxst),ucstplo(mxst)
c
      real*8 suqmo, suinto, stabtot, stabto 
      common/rfene6/suqmo(mxst),suinto(mxst),stabtot(mxst),stabto(mxst)
c
      real*8 extelg, extnucg, sextelg, sxtnucg, selextg, snucxtg, uclasg
      real*8 uclaseg, uclasdg, uclasrg, suclasg, upolclg, repmodg
      common /rfene7/ extelg(mxgran,mxst),extnucg(mxgran,mxst),
     +  sextelg(mxgran,mxst),sxtnucg(mxgran,mxst),selextg(mxgran,mxst),
     +  snucxtg(mxgran,mxst),uclasg(mxgrpar),uclaseg(mxgrpar),
     +  uclasdg(mxgrpar),uclasrg(mxgrpar),suclasg(mxgran,mxgran),
     +  upolclg(mxgran),repmodg(mxgran,mxst)
c
      real*8 sxtmolg, smolxtg, stotxtg
      common /rfene8/ sxtmolg(mxgran,mxst),smolxtg(mxgran,mxst),
     +  stotxtg(mxgran,mxst)
c
      real*8 uelstg, uintg, suintg, uneqclg, ustaclg
      common /rfene9/ uelstg(mxgran,mxst),uintg(mxgran,mxst),
     +  suintg(mxgran,mxst),uneqclg(mxgran),ustaclg(mxgran)
c
      real*8 sxtelog, sextnog, selxtog, sncexog, suclsog, uplclog
      common /rfene10/ sxtelog(mxgran,mxst),sextnog(mxgran,mxst),
     +  selxtog(mxgran,mxst),sncexog(mxgran,mxst),
     +  suclsog(mxgran,mxgran),uplclog(mxgran)
c
      real*8 sextmog, smlexog, sttexog
      common/rfene11/sextmog(mxgran,mxst),smlexog(mxgran,mxst),
     +  sttexog(mxgran,mxst)
c
      real*8 suintog
      common /rfene12/ suintog(mxgran,mxst)
c
c
      real*8 unucrepo, eoneelo, ekino, enuao, etwoelo
      real*8 extelo, extnuco
      real*8 uneqno, uneqelo, uclaso, uclaseo, uclasdo, uclasro
      real*8 rdispo, repmo, ustano, ustaelo
      common /rfene21/ unucrepo,eoneelo(mxst),ekino(mxst),enuao(mxst),
     +  etwoelo(mxst),extelo(mxst),extnuco(mxst),
     +  uneqno(mxst),uneqelo(mxst),uclaso,uclaseo,
     +  uclasdo,uclasro,
     +  rdispo(mxst),repmo(mxst),
     +  ustano(mxst),ustaelo(mxst)
c
      real*8 ucststo
      common/rfene22/
     +  ucststo(mxst)
c
      real*8 uscfo, uelsto, uinto, uneqo, uneqqmo, uneqclo
      real*8 upolno, ustao, ustaqmo, ustaclo
      common /rfene23/ uscfo(mxst),uelsto(mxst),uinto(mxst),
     +  uneqo(mxst),uneqqmo(mxst),uneqclo,
     +  upolno(mxst),ustao(mxst),ustaqmo(mxst),ustaclo
c
      real*8 extelgo, extnucgo, uclasgo
      real*8 uclasego, uclasdgo, uclasrgo, repmodgo
      common /rfene27/ extelgo(mxgran,mxst),extnucgo(mxgran,mxst),
     +  uclasgo(mxgrpar),uclasego(mxgrpar),
     +  uclasdgo(mxgrpar),uclasrgo(mxgrpar),
     +  repmodgo(mxgran,mxst)
c
      real*8 uelstgo, uintgo, uneqclgo, ustaclgo
      common /rfene29/ uelstgo(mxgran,mxst),uintgo(mxgran,mxst),
     +  uneqclgo(mxgran),ustaclgo(mxgran)
c
c
      integer nambpt
      common /drfamb/ nambpt(mxat)
c
      real*8 znamb
      common /drfamb2/ znamb(mxat)
c
c
      integer nsamp, nblock, nmoves, ncheck
      common /mcrunp1/ nsamp,nblock,nmoves,ncheck
c
      real*8 darot, dtrans, temp, onekt, excld, delmax, enmin
      common /mcrunp2/ darot,dtrans,temp,onekt,excld,delmax,enmin
c
      logical mcupdt,secstat,excess
      integer notrans, norot, imcout, iseed, imcst, iactst
      common /mcrunp3/ notrans,norot,imcout,iseed,imcst,iactst,
     +               mcupdt,secstat,excess
c
      character*16 outfor
      common /mcrunp4/ outfor
c
      real*8 ratmin, ratmax, amxrot, amxtrn, amnrot, amntrn
      common /mcrunp5/ ratmin,ratmax,amxrot,amxtrn,amnrot,amntrn
c
      real*8 gammc
      common /mcrunp6/ gammc(5)
c
      dimension cz(*),c(3,*)
      data dzero /0.0d0/
      enuc = dzero
      enc  = dzero
      call vdwaals_energy(enc)
      enuc = enuc + enc
cdrf
cdrf  no amibiguous nuclei regarded here
cdrf  code can be found on the hondo8-file scfnew.f
      if ((n .eq. 1) .and. (field(:4) .ne. 'scf') 
     1 .and. (field(5:) .ne. 'scf')) then
          enucf=enuc
          return
      endif
c 
c if (indi)  core-core repulsion included here
c
      enc = dzero
      owarn = .false.
      do i = 2,n
         ni = i-1
         obqi = (omtslf .and. zaname(i)(1:2).eq.'bq')
         do  j = 1,ni
            obqj = (omtslf .and. zaname(j)(1:2).eq.'bq')
            if(.not. (obqi .and. obqj)) then
               rr = dzero
               do k = 1,3
                  rr = rr+(c(k,i)-c(k,j))**2
               enddo
               dsqrr=dsqrt(rr)
               if (indi) then
                 l  = 0
                 do k=1,npairs
                   indij = ichgat(i).eq.indx(k).and.ichgat(j).eq.jndx(k)
                   indji = ichgat(i).eq.jndx(k).and.ichgat(j).eq.indx(k)
                   if (indij.or.indji) l=k
                 enddo
                 if (l.ne.0)enc = enc + d(l)*dexp(-eta(l)*dsqrr)
               endif
               enuct = cz(i)*cz(j)/dsqrr
               enuc = enuc+enuct
            else
               if(abs(cz(i)) .gt. 0.0001d0 .and. 
     &            abs(cz(j)) .gt. 0.0001d0)owarn = .true.
            endif
         enddo
      enddo
      enuc = enuc+enc
c
c informational messages (behaviour has changed in v 6.2)
c
      if(owarn)then
         write(iwr,100)
      else if (.not. omtslf)then
         write(iwr,101)
      endif
 100  format(//1x,'Nuclear energy has been computed in the presence',
     &     ' of a point charge field',/,
     &     1x,'Interaction between bq centres has been omitted.',
     &     ' To include it, include the class 1 directive BQBQ')
 101  format(1x,//'Nuclear energy includes self-energy of point ',
     &     'charge field')

cdrf-------------- unucrep contains vac-repulsion
      unucrep = enuc
cdrf----------------extension for inclusion of nuc-nuc-rep
      if (field(:4) .eq. 'scf') then
        enuc = enuc + extnuc(iactst)
        if (neqsta .eq. 1) enuc = enuc + ustanuc(iactst)
        if (neqrf .eq. 1) enuc = enuc + uneqnuc(iactst)
      endif
      if (field(5:) .eq. 'scf' .and. (iarfcal .eq. 0)) then
        enuc = enuc + snucext(iactst) +
     1                sextnuc(iactst)
      endif
      if ((field(5:) .eq. 'scf') .and. (iarfcal .eq. 0)) then
        enuc = enuc + scffact*snucnuc(iactst) 
      endif
cxxx  if ((field(:4) .eq. 'scf') .or. (field(5:) .eq. 'scf')) then
      if (field(:4) .eq. 'scf') then
        enuc = enuc + repmod(iactst)
      endif
cdrf------------------------
      enucf=enuc
      return
      end

      function crnuc(n,cz,c,etmp)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      dimension cz(*),c(3,*)
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
      data dzero /0.0d0/
      etmp = dzero
      crnuc = dzero
      if (n .eq. 1) return
      do 130 i = 2,n
         ni = i-1
         obqi = (zaname(i)(1:2).eq.'bq')
         do 120 j = 1,ni
            obqj = (zaname(j)(1:2).eq.'bq')
            if(obqi.and.obqj)goto 120
            rr = dzero
            do 100 k = 1,3
               rr = rr+(c(k,i)-c(k,j))**2
 100        continue
            e = cz(i)*cz(j)/dsqrt(rr)
            if(obqi.or.obqj)then
               e = e * 0.5d0
               etmp = etmp + e
            endif
            crnuc = crnuc + e
 120     continue
 130  continue
      write(iwr,160)crnuc
 160  format(1x,'    nuclear energy = ',f20.10)
      write(iwr,170)etmp
 170  format(1x,'    component from crystal field = ',f20.10)
      return
      end
      subroutine orfog(q,qp,b,c,iky,ilifq,newbas,nrow,iop)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension q(*),qp(*),b(*),c(*),iky(*),ilifq(*)
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
      dimension p(maxorb)
c... to orthogonalize the cols of q - result to qp
c... overlap matrix supplied in b - destroyed on exit
c... scratch triangle in c
c... iop=1 normal mode iop=2 qp set as if q was given as identity
c... qp can overwrite q
      top=dsqrt(1.0d0/b(1))
      b(1)=top
      do 1 i=2,newbas
      m=iky(i)
   1  c(m+1)=b(m+1)*top
      do 21 k=2,newbas
      m=iky(k)
      top=b(m+k)
      ll=k-1
      do 22 i=1,ll
       bot=c(m+i)
      p(i)=-bot
   22 top=top-bot*bot
      SAFMIN = DLAMCH( 'Safe minimum' )
      top = max(top,safmin)    !deal with division-by-zero -XJ
      top=dsqrt(1.0d0/top)
      b(m+k)=top
      n=k+1
       if(n.gt.newbas)goto 23
      do 24 l=n,newbas
      nsp=iky(l)
      c(nsp+k)=(ddot(ll,p,1,c(nsp+1),1)  +b(nsp+k))*top
   24 continue
   23 do 21 l=1,ll
      bot=0.0d0
      do 27   i=l,ll
   27 bot=p(i)*b(l+iky(i))+bot
   21 b(l+m)=bot*top
       goto (30,31),iop
   31    do 40 i=1,newbas
      m=ilifq(i)
      n=iky(i)
      do 40 j=1,nrow
      top=0.0d0
      if(j.le.i)top=b(n+j)
   40 qp(m+j)=top
      go to 50
   30 do 41 i=1,nrow
      call dgthr(newbas,q(i+1),p,ilifq)
      m=1
      do 41 j=1,newbas
       qp(i+ilifq(j))=ddot(j,b(m),1,p,1)
   41  m=m+j
   50 continue
      return
      end
       subroutine tranp(p1,p2)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      dimension p1(*),p2(*)
      dimension ct1(mxorb3),it1(mxorb3)
c
      real*8 ctran
      integer ilifc, ntran, itran
      logical otran,otri
      common/tran/ilifc(maxorb),ntran(maxorb),itran(mxorb3),
     +            ctran(mxorb3),otran,otri
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
      l2=ikyp(num)
c
      if(otran)goto 100
      m=0
      do 1 i=1,num
      nt1=ntran(i)
      do 2 k=1,nt1
      l=k+ilifc(i)
      ct1(k)=ctran(l)
2     it1(k)=itran(l)
      do 1 j=1,i
      m=m+1
      p2(m)=0.0d0
      nt2=ntran(j)
      do 1 k=1,nt2
      kk=k+ilifc(j)
      l=itran(kk)
      do 1 ii=1,nt1
      ll=it1(ii)
      top=p1(iky(max(l,ll))+min(ll,l))
1     p2(m)=ct1(ii)*ctran(kk)*top+p2(m)
      return
 100  continue
      call dcopy(l2,p1(1),1,p2(1),1)
      return
      end
      subroutine ciexpr(kcorb,nconf,nham,ia,l1)
c
c     ----- update f, alpha, and beta matrices based on
c           the new ci coefficients -----
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      real*8 cicoef, f, alpha, beta
      integer no, nco, nseto, npair, ncores, ibm, nset, nopen
      integer  nope, noe
      logical old
      common/scfwfn/cicoef(2,12),f(25),alpha(325),beta(325),
     + no(10),nco,nseto,npair,ncores,ibm,old,nset,nopen,nope,noe(10)
c
      dimension ia(l1),kcorb(2,12),nconf(*)
      data dzero,two,tm6/0.0d0,2.0d0,1.0d-06/
      test = tm6
      nocio = 2
c
c     ----- normalize ci coefficients and update f -----
c
      do 140 kpair = 1,npair
      sum = ddot(nocio,cicoef(1,kpair),1,cicoef(1,kpair),1)
      if (sum .lt. test) sum = test
      sum = 1.0d0/dsqrt(sum)
      do 120 korb = 1,nocio
      cicoef(korb,kpair) = cicoef(korb,kpair)*sum
      kk = kcorb(korb,kpair)
      nx = nconf(kk)
  120 f(nx) = cicoef(korb,kpair)**2
  140 continue
c
c     ----- update alpha and beta -----
c
      do 180 kpair = 1,npair
      kone = kcorb(1,kpair)
      ktwo = kcorb(2,kpair)
      lo = nconf(kone)
      lt = nconf(ktwo)
      lolo = ia(lo)+lo
      lolt = ia(lt)+lo
      ltlt = ia(lt)+lt
      do 160 k = 1,nham
      lok = ia(max(lo,k))+min(lo,k)
      ltk = ia(max(lt,k))+min(lt,k)
      alpha(lok) = two*f(lo)*f(k)
      alpha(ltk) = two*f(lt)*f(k)
      beta(lok) = -f(lo)*f(k)
      beta(ltk) = -f(lt)*f(k)
  160 continue
      ann = cicoef(1,kpair)*cicoef(2,kpair)
      alpha(lolo) = f(lo)
      alpha(ltlt) = f(lt)
      beta(lolo) = dzero
      beta(ltlt) = dzero
      beta(lolt) = ann
      alpha(lolt) = dzero
  180 continue
      return
      end
      subroutine dengvb(trans,d,ifock,ia,l1)
c
c     ----- dengvb calculates the density matrix for all molecular
c            orbitals in the ifock'th fock operator -----
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      common/diisd/sta(270),iposa(24),stb(270),iposb(24),
     + ns(15),melsym(15),nconf(maxorb),msympr(maxorb),norb,nham,
     + nsytyp,nsos,
     + ti(25),e1i(25),asym,
     + cilow(12),jpair,kone,ktwo,kcorb(2,12),
     + iojk(49),ioham(25),iojkao(49)
      dimension ia(*),d(*),trans(l1,*)
      data dzero/0.0d0/
      call vclr(d,1,ia(l1+1))
      do 100 k = 1,norb
      nx = nconf(k)
      if (nx .ne. ifock) go to 100
      jbf = msympr(k)
      ij=1
      do 120 i=1,l1
      if(trans(i,jbf).ne.dzero)  then
      call daxpy(i,trans(i,jbf),trans(1,jbf),1,d(ij),1)
      endif
  120 ij=ij+i
  100 continue
      return
      end
      subroutine sto31g(ztext,ibasis,igauss)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension xng(6)
      data xng/
     *'1','2','3','4','5','6'/
      data yblnk,ydunn/' ','dunn' /
c
      ytext=ytrunc(ztext)
      if(ytext.eq.ydunn.or.ytext.eq.yblnk)go to 5
      xmg(1:1)=ztext(1:1)
      do 1 igauss=1,6
      if(xmg.eq.xng(igauss))go to 2
 1    continue
      call caserr2('invalid n specified in n-m1g keyword')
 2    xmg(1:1)=ztext(3:3)
      do 3 ibasis=2,3
      if(xmg.eq.xng(ibasis))go to 4
 3    continue
      call caserr2('invalid m specified in n-m1g')
 4    ibasis=ibasis-2
      return
 5    ibasis=3
c
      return
      end
      subroutine tdown(qnew,ilifn,q,ilifq,nnn)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      dimension qnew(*),q(*),ilifn(*),ilifq(*)
      dimension dum(maxorb)
c
      real*8 ctran
      integer ilifc, ntran, itran
      logical otran,otri
      common/tran/ilifc(maxorb),ntran(maxorb),itran(mxorb3),
     +            ctran(mxorb3),otran,otri
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
      if(otran)go to 60008
      do 731 i=1,nnn
      m=ilifq(i)
      call vclr(dum,1,num)
      do 733 j=1,num
      n=ntran(j)
      do 733 k=1,n
      l=ilifc(j)    +k
733   dum(itran(l))=ctran(l)*q(m+j)+dum(itran(l))
      call dcopy(num,dum(1),1,qnew(ilifn(i)+1),1)
 731  continue
c 
60010 if (nnn.lt.num) then
c...    clear vectors that are extra
         do i=nnn+1,num
            call vclr(qnew(ilifn(i)+1),1,num)
         end do
      end if
c
      return
60008 do 60009 i=1,nnn
      call dcopy(num,q(ilifq(i)+1),1,qnew(ilifn(i)+1),1)
60009 continue
      return
      end

      subroutine wrt3i(ij,nword,iblk,num3)
c
c     this routine writes nword(s) of integer data, commencing at
c     block iblk on unit num3, from the array ij.
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension ij(*)
      common/junkcc/dchar(511)
      dimension iq(1)
      equivalence (dchar(1),iq(1))
      call search(iblk,num3)
      j=1
      k=nword
 20   if(k)30,30,10
 10   call put(ij(j),min(k,511),num3)
      j=j+511
      k=k-511
      go to 20
30    return
      end
      subroutine readi(ij,nword,iblk,num3)
c
c     this routine reads nwords of integer data, commencing at
c     block iblk on unit num3 into the array ij.
c     note that nword has to be exactly what lies on disk
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension ij(*)
      common/junkcc/dchar(511)
      logical odd
      dimension iq(1)
      equivalence (iq(1),dchar(1))
      call search(iblk,num3)
      j=1
 20   if(j.gt.nword)go to 30
      call find(num3)
      call get(ij(j),l)
      j=j+l
      go to 20
30    if(j-1.ne.nword) then
       write(6,*)
     + 'invalid no of words: requested,present = ', nword, j-1
       call caserr2('invalid number of words in readi')
      endif
      return
      end

      subroutine wrt3is(ij,nword,num3)
c
c     this routine writes nword(s) of integer data, commencing at
c     the current position/block on unit num3, from the array ij.
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension ij(*)
      common/junkcc/dchar(511)
      dimension iq(1)
      equivalence (dchar(1),iq(1))
      j=1
      k=nword
 20   if(k)30,30,10
 10   call put(ij(j),min(k,511),num3)
      j=j+511
      k=k-511
      go to 20
30     return
       end
      subroutine readis(ij,nword,num3)
c
c     this routine reads nword(s) of integer data, commencing at
c     the current position/block on unit num3 into the array ij.
c     note that nword has to be exactly what lies on disk.
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension ij(*)
      logical odd
      common/junkcc/dchar(511)
      dimension iq(1)
      equivalence (iq(1),dchar(1))
      j=1
 20   if(j.gt.nword)go to 30
      call find(num3)
      call get(ij(j),l)
      j=j+l
      go to 20
 30   if(j-1.ne.nword) then
       write(6,*)
     + 'invalid no of words: requested,present = ', nword, j-1
       call caserr2('invalid number of words in readis')
      endif
      return
      end

      subroutine restre
      implicit real*8  (a-h,o-z)
c ...
c ... restore restart section on dumpfile
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
      integer master, indxi, indxj, nfiles, junits, jblkrs, jblkas
      integer nfilef, junitf, jblkrf, jblkaf, isecbl, iscftp
      integer lword4, ilow4, ncol4, nsa4, newb4, nbas4, lenb4
      integer ndump4, iblkq4, nblkq4, lena4, nbb4, ionsv4, isecv4
      logical oprin4, oindx4
      integer npas41, npas42, iacc4
c
      common /cndx40/ master,indxi,indxj,nfiles,junits,jblkrs,jblkas,
     + nfilef,junitf,jblkrf,jblkaf,isecbl,iscftp,
     + lword4,ilow4,ncol4,nsa4,newb4,nbas4,lenb4,ndump4,iblkq4,
     + nblkq4,lena4,nbb4,
     + oprin4(10),oindx4,ionsv4,npas41,npas42,iacc4,isecv4
c
       integer len_cndx40
       parameter (len_cndx40=37)
c      used: restre(util1),revise(util1),utyp21(server)
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer ibl3s, ibl3t, ibl3f, ibl3x, ibl3y, ibl3z
      integer ibl3rs, ibl3qa, ibl3pa, ibl3ea, ibl3qb, ibl3pb
      integer ibl3eb, ibl3g, ibl3hs, ionsec, ibl3op, ions2
      integer isecqm
      common /dump3/ ibl3s, ibl3t, ibl3f, ibl3x, ibl3y, ibl3z,
     +               ibl3rs,ibl3qa,ibl3pa,ibl3ea,ibl3qb,ibl3pb,
     +               ibl3eb,ibl3g,ibl3hs,ionsec,ibl3op,ions2,
     +               isecqm
c
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
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
      logical irspl
      common/junko/cspace(30),irspa(700),irspl(40),irspb(1590),
     *             irspc(40),irspd(296),irspe(8)
      data m21/21/
c
      call secget(isect(501),m21,ibl3rs)
      call rdchr(title,ldsect(isect(501)),ibl3rs,idaf)
      call reads(cspace,lds(isect(501)),idaf)
c
      call dcopy(len_restrr,cspace,1,gx,1)
      call icopy(len_restar,irspa,1,nprint,1)
      call icopy(len_restrl,irspl,1,ciopt,1)
      call icopy(len_restri,irspb,1,jjfile,1)
      call icopy(len_cndx40,irspc,1,master,1)
      call icopy(len_cndx41,irspd,1,ncoorb,1)
      call icopy(len_infoa,irspe,1,nat,1)
c
      return
      end
      subroutine revise
      implicit real*8  (a-h,o-z)
c ...
c ... update restart section on dumpfile
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
      integer master, indxi, indxj, nfiles, junits, jblkrs, jblkas
      integer nfilef, junitf, jblkrf, jblkaf, isecbl, iscftp
      integer lword4, ilow4, ncol4, nsa4, newb4, nbas4, lenb4
      integer ndump4, iblkq4, nblkq4, lena4, nbb4, ionsv4, isecv4
      logical oprin4, oindx4
      integer npas41, npas42, iacc4
c
      common /cndx40/ master,indxi,indxj,nfiles,junits,jblkrs,jblkas,
     + nfilef,junitf,jblkrf,jblkaf,isecbl,iscftp,
     + lword4,ilow4,ncol4,nsa4,newb4,nbas4,lenb4,ndump4,iblkq4,
     + nblkq4,lena4,nbb4,
     + oprin4(10),oindx4,ionsv4,npas41,npas42,iacc4,isecv4
c
       integer len_cndx40
       parameter (len_cndx40=37)
c      used: restre(util1),revise(util1),utyp21(server)
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
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
      integer ibl3s, ibl3t, ibl3f, ibl3x, ibl3y, ibl3z
      integer ibl3rs, ibl3qa, ibl3pa, ibl3ea, ibl3qb, ibl3pb
      integer ibl3eb, ibl3g, ibl3hs, ionsec, ibl3op, ions2
      integer isecqm
      common /dump3/ ibl3s, ibl3t, ibl3f, ibl3x, ibl3y, ibl3z,
     +               ibl3rs,ibl3qa,ibl3pa,ibl3ea,ibl3qb,ibl3pb,
     +               ibl3eb,ibl3g,ibl3hs,ionsec,ibl3op,ions2,
     +               isecqm
c
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
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
      logical irspl
      common/junko/cspace(30),irspa(700),irspl(40),irspb(1590),
     *             irspc(40),irspd(296),irspe(8)
      data m21/21/
c
      len501=lensec(ldsect(isect(501))) + lensec(lds(isect(501)))
      call secput(isect(501),m21,len501,ibl3rs)
c
      call dcopy(len_restrr,gx,1,cspace,1)
      call icopy(len_restar,nprint,1,irspa,1)
      call icopy(len_restrl,ciopt,1,irspl,1)
      call icopy(len_restri,jjfile,1,irspb,1)
      call icopy(len_cndx40,master,1,irspc,1)
      call icopy(len_cndx41,ncoorb,1,irspd,1)
      call icopy(len_infoa,nat,1,irspe,1)
c
      call wrtc(title,ldsect(isect(501)),ibl3rs,idaf)
      call wrt3s(cspace,lds(isect(501)),idaf)
c...  check
      if (len_cndx41.gt.296) call caserr('common/cndx41 too big')
c
      return
      end
      function ispchk()
c
c     ----- check for nature of basis -----
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      common/junk/
     *kstart(mxshel),katom(mxshel),ktype(mxshel),
     *kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),nshell,non3(3)
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
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
      logical ciopt, ciforc, mp2, hfgr, bfgs, ump2, lmeth2
      logical ump3, rmp3, ordmo, mp2w, loptor, ladp, lcpf
      logical lopti, lmcscf, lforce, lci, lcart, lmcdat
      logical lfdtrn, unit7, lcontr, lvcd, lgten
      logical ldenom, ignore, ldens, lset, ladapt, lsym, latmol
      logical berny, llibry, limpt, fpres, oss, ldiag, lskip
      logical opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common /restrl/ ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,
     +rmp3,ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
c
      data mone/1/
c
      call secget(isect(491),mone,iblkv)
      m110=10+maxat
      m1420=mxprim*5+maxat
      len=lensec(mxprim)+lensec(m1420)+lensec(m110)
      ibl491 = iblkv + len
      nav = lenwrd()
      call readi(kstart,mach(2)*nav,ibl491,idaf)
c
      itype = 1
      opbas = .false.
      odbas = .false.
      ofbas = .false.
      ogbas = .false.
      do 100 loop = 1,nshell
      if (ktype(loop) .eq. 2) opbas = .true.
      if (ktype(loop) .eq. 3) odbas = .true.
      if (ktype(loop) .eq. 4) ofbas = .true.
      if (ktype(loop) .eq. 5) ogbas = .true.
  100 continue
      if (opbas) then
      itype = 2
      endif
      if (odbas) then
      itype = 3
      endif
      if (ofbas) then
      itype = 4
      endif
      if (ogbas) then
      itype = 5
      endif
      ispchk = itype
      return
      end
      function locats(ztext,ns,ytext)
      implicit real*8  (a-h,p-w),integer (i-n),logical (o)
      implicit character*8 (z),character*1 (x)
      implicit character*4 (y)
      character*(*) ztext,ytext
      dimension ztext(*)
      do 10 i = 1 , ns
            loop = i
            if( index(ztext(i),ytext) ) 10,10,20
  10  continue
      locats = 0
      return
  20  locats = loop
      return
      end
      subroutine linear(a,b,mdim,n)
      implicit real*8  (a-h,p-w), integer (i-n), logical (o)
      implicit character*8 (z), character*1 (x), character*4 (y)
c
c places symmetric square array in a linear form
c
      dimension a(mdim,*),b(*)
      k = 1
      do 10 j = 1 , n
         do 20 i = 1 , j
            b(k) = a(i,j)
            k = k + 1
20    continue
10    continue
      return
      end
      subroutine squar(a,b,mdim,n,key)
      implicit real*8  (a-h,p-w), integer (i-n), logical (o)
      implicit character*8 (z), character*1 (x), character*4 (y)
c
c     places linear array in square form
c     two parts ... key =
c     1  square array to be formed not symmetric
c     0  square array to be formed symmetric
c
      dimension a(*),b(mdim,*)
c
      if(key)1,3,1
c
c     ---- key=1  array not symmetric ----
c
    1 k=n*n
      do  10 j = 1 , n
          jx = n - j + 1
          do 20 i = 1 , n
             ix = n - i + 1
             b(ix,jx) = a(k)
             k = k - 1
20    continue
10    continue
      return
c
c     ---- key=0  array symmetric ----
c
    3 k = n * ( n + 1 ) / 2
      do  30 j = 1 , n
          jx = n - j + 1
          do 40 i = 1 , jx
             ix = jx - i + 1
             b(ix,jx) = a(k)
             k = k - 1
40    continue
30    continue
      do 50 j = 1 , n
         do 60 i = 1 , j
            b(j,i) = b(i,j)
60    continue
50    continue
      return
      end
      subroutine demoao(dmo, dao, coorb, t,nbasis,ncoorb,ndim)
c
c  transform density matrix from mo basis to ao basis
c
c  dmo density matrix in mo basis (triangular form)
c  dao array for density matrix in ao basis (triangular form)
c  coorb mo coefficients
c  t workspace
c  nbasis number of basis functions
c
      implicit real*8  (a-h,o-z)
      dimension dmo(*), dao(*), coorb(ndim,*), t(*)
c
      ij = 0
      do 60 j = 1 , nbasis
         t(1) = dmo(1)*coorb(j,1)
         kk = 1
         do 30 k = 2 , ncoorb
            x = 0.0d0
            km1 = k - 1
            do 20 l = 1 , km1
               x = x + dmo(kk+l)*coorb(j,l)
               t(l) = t(l) + dmo(kk+l)*coorb(j,k)
 20         continue
            t(k) = x + dmo(kk+k)*coorb(j,k)
            kk = kk + k
 30      continue
         do 50 i = 1 , j
            x = 0.0d0
            do 40 k = 1 , ncoorb
               x = x + coorb(i,k)*t(k)
 40         continue
            ij = ij + 1
            dao(ij) = x
 50      continue
 60   continue
      return
c
      end
      subroutine prtris(f,numscf,iw)
      implicit real*8  (a-h,o-z)
      dimension f(*)
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
      dimension n(8)
      mmax = 8
      imax = 0
 20   imin = imax + 1
      imax = imax + mmax
      if (imax.gt.numscf) imax = numscf
      write (iw,6010)
      write (iw,6020) (i,i=imin,imax)
      write (iw,6010)
      do 40 j = 1 , numscf
         nn = 0
         do 30 i = imin , imax
            nn = nn + 1
            n(nn) = iky(max(i,j)) + min(i,j)
 30      continue
         write (iw,6030) j , (f(n(i)),i=1,nn)
 40   continue
      if (imax.lt.numscf) go to 20
      return
 6010 format (/)
 6020 format (5x,8(6x,i3,6x))
 6030 format (i5,8f15.8)
      end
      function lenint(n)
      implicit real*8  (a-h,o-z)
c
c     lenint returns the number of reals corresponding to n integers
c     for all n greater than zero. For n equals zero it returns one.
c     Valid inputs are all non-negative integers.
c
c     Proof:
c
c     Define m = lenwrd() which means that m can be either 1 or 2.
c
c     Assume n = i*m, i > 0 then
c     lenint = (i*m-1)/m+1
c            = i-1      +1
c            = i
c     Assume n = i*m+j, i >= 0 and 0 < j < m
c     lenint = (i*m+j-1)/m+1
c            = (i*m)/m+(j-1)/m+1
c            = i+1
c     Assume n = 0
c     lenint = 0/m+1
c            = 1
c     
c     The max is required to handle the n=0 and m=1 case correctly.
c
      lenint = max(n-1,0)/lenwrd()+1
      return
      end
      function lenrel(n)
c     lenrel is the number of integers corresponding to n reals
      implicit real*8  (a-h,o-z)
      lenrel = n*lenwrd()
      return
      end
      subroutine prsqm(w,nbasis,ncol,ndim,iw)
      implicit real*8  (a-h,o-z)
      dimension w(*)
      m = 1
      n = 6
 20   if (ncol.lt.m) return
      if (n.gt.ncol) n = ncol
      write (iw,6010) (i,i=m,n)
      write (iw,6020)
      do 30 j = 1 , nbasis
         write (iw,6030) j , (w(j+(i-1)*ndim),i=m,n)
 30   continue
      m = m + 6
      n = n + 6
      go to 20
 6010 format (/5x,6i12)
 6020 format (/)
 6030 format (i5,6f12.6)
      end
      subroutine squr(t,sq,n)
      implicit real*8  (a-h,o-z)
c     triangle to square
      dimension sq(2),t(2)
      ij = 0
      ii = 0
      do 30 i = 1 , n
         jj = 0
         do 20 j = 1 , i
            ij = ij + 1
            sq(ii+j) = t(ij)
            sq(jj+i) = t(ij)
            jj = jj + n
 20      continue
         ii = ii + n
 30   continue
      return
      end
      subroutine archiv(option,iwrite)
      implicit real*8  (a-h,o-z)
************************************************************************
*
*   archiv positions the mopac archive file
*   on unit 12 - for data retrieval by gamess.
*   zmatrix, geometry  and vectors defined by option
*
************************************************************************
      character *(*) option
      character *10 string
      character*(132) filatt,zarch
      character *80 line
c
c     mxtoken - the maximum number of tokens on an input line
c               for space separated tokens: 132/2 plus 1 to allow
c               for look ahead (132 is the line-length defined in
c               common/workc)
c
      integer mxtoken
      parameter(mxtoken=67)
      integer jrec, jump, istrt, inumb, iwidth, nend, nstart
      integer nline, noline, jwidth, nerr
      logical oswit, oterm, oflush
      common/work/jrec,jump,istrt(mxtoken),inumb(mxtoken),iwidth,
     + nend(mxtoken),nstart(mxtoken), nline,noline,jwidth,nerr,
     + oswit,oterm,oflush
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      data zarch/'archive'/
c
      iread = 12
      filatt = ' '
      call gtnv(zarch,filatt)
      if (filatt.eq.' ') filatt = zarch
      do loop = len(filatt),1,-1
        if(filatt(loop:loop).ne.' ') goto 30
      enddo
   30 continue
c
      if (option.eq.'zmat') then
       string = 'zmatrix an'
       jrec = 0
      else if (option.eq.'vect') then
       string = 'vectors mo'
       jrec = 0
      else if (option.eq.'geom') then
       string = 'geometry m'
      else
       call caserr2('invalid request for archive data')
      endif
      irdd = 0
c
      close (iread)
      open(unit=iread,file=filatt,form='formatted',status='unknown')
      rewind iread
   10 read(iread,'(a80)',end=20,err=20)line
      irdd = irdd + 1
      if(index(line,string).eq.0) go to 10
      ird = iread
      return
   20 rewind iread
      write(iwrite,'(a)')' archive file missing or empty'
      call caserr2('requested archive file missing or empty')
      return
      end
      subroutine ver_util1(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/util1.m,v $
     +     "/
      data revision /"$Revision: 6246 $"/
      data date /"$Date: 2011-10-31 23:16:26 +0100 (Mon, 31 Oct 2011) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
      subroutine check_feature(module)
c
c...  called to indicate whether a module can handle a feature, so that if one does
c...  not copy all features to an improved routine, one can make sure the user is warned.
c...  Initially set up for ZORA, which was not implemented in direct mode. Thus no ZORA was done
c...  which was easily missed.  
c
      parameter (nmodules=14,ncheck=2)
      character*(*) module
      character*10 modules(0:ncheck,nmodules),type(0:ncheck)
c     
c
c     ZORA common
c
      logical ozora,oscalz,onlyp,onlyd,onlyf,onlyg,opzora,
     1        op2zora,onlys,ocontr,oint_zora,o1e_zora,osmall_zora,
     2        opre_zora,oso,onoext_z,oioraz,oatint_z,oatscf_z,
     3        oscalatz
      real*8 cspeed,critat_z
      integer nbas_zora,ibls_z,iblv_z,nwv_z,ibiso_z,nwiso_z,ibcor_z,
     1        ibsx_z,icoul_z,ibsy_z,ibsz_z,ibsmin_z,ibsplus_z,nwmat_z,
     2        ibtrin_z,ibtrout_z,ibscalf_z,nwtrat_z,ibscale_z,ibcoul_z,
     3        niter_z,int_zora,num_ext,igauge_z,nwcor_z,isexch_z,
     4        nat_z,ibdens_z,nwdens_z,ibscalao_z,ibshat_z,is_z,
     5        irest_z,numdu_z,ibldu_z,iblock_z,lshift_z,ibl7ew_z
      common/zorac/ cspeed,critat_z,ozora,oscalz,
     1              onlyp,onlyd,onlyf,onlyg,opzora,
     1              op2zora,onlys,ocontr,nbas_zora,
     2              oint_zora,o1e_zora,icoul_z,
     3              ibls_z,iblv_z,nwv_z,
     3              ibiso_z,nwiso_z,ibcor_z,nwcor_z,
     4              ibsx_z,ibsy_z,ibsz_z,ibsmin_z,ibsplus_z,nwmat_z,
     5              ibtrin_z,ibtrout_z,ibscalf_z,nwtrat_z,
     6              ibscale_z,ibcoul_z,ibdens_z,nwdens_z,osmall_zora,
     7              niter_z,opre_zora,int_zora,num_ext,isexch_z,
     8              oso,onoext_z,is_z,igauge_z,oioraz,
     9              nat_z,ibscalao_z,oatint_z,oatscf_z,ibshat_z,
     1              oscalatz,irest_z,numdu_z,ibldu_z,iblock_z,lshift_z,
     2              ibl7ew_z
c....    icoul_z : 1 : Full coulomb
c...               2 : atomic coulomb (default)
c...               3 : no coulomb (just 1-electron ints)
c...               4 : small basis coulomb
c...     isexch_z  0 : no exchange in zora correction (default)
c...               1 : exchange in zora correction (like in dft)
c...     niter_z : # iterations between coulomb matrix update
c...     int_zora :internal basis gen. : 0 none, 1, automatic, 2 by hand
c
c...     oscalz  : scaling on or off (true or false)
c...     oioraz  : change metric to perform inf. order scaling 
c...               (Dyall/EvLenthe (not quite))
c...     oscalatz: atoms only scaled zora (together with get atom)
c...     only'x' : internal basis only op to 'x'-type functions
c...     op(2)zora: print flags (op2zora+opzora: prints *everything*)
c...     ocontr  : adds l+1 block to internal basis in contracted way 
c...     nbas_zora: size of internal basis
c...     oint_zora: in internal or external mode???
c...     o1e_zora: one electron integrals have been calculated in *this*
c...               internal basis
c...     osmall_zora: small swap for projected coulomb
c...     opre_zora: pre_zora has been called
c...     oso     : spin orbit calculation 
c...     onoext_z: no external functions copied to internal basis
c...     is_z : section for orbitals used for  coulomb matrix
c...     igauge_z : add so many s-marices to h-matrix
c...                the energy should change by igauge
c...     nat_z  : indicates (if ne 0) that starting density is used
c...              to calc. coulomb; gives # it. to get atomic density
c...     critat_z : accuracy of atomic energy for atomic zora
c...     oatint_z : indicates that integrals to be calculated are atomic
c...     oatscf_z : indicates we are doing a set of atomic SCFs
c...     ibshat_z : block on num8 to store blocked s/h
c...     for disk-adress and lengths (ib.., nw.. ) see subroutine atoms2
c...     irest_z  : .ne.0 indicates atomic zora corrections from dumpfile; 1 not on restart
c...     numdu_z,ibldu_z : if ne 0 other dumpfile to read zora corrections from
c...     iblock_z is the block number on ed3, where the zora corrections (section 425) reside
c...     ibl7ew is where the energy weighted density matrix resides
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      character*60 error
c
      data type   /'module',' zora ', 'general'/
      data modules/'denscf',   'all', ' '       ,
     2             'gvbitr',   'all', ' '       ,
     3             'rhfcld',   'all', ' '       ,
     4             'rhfclm',   'all', ' '       ,
     5             'drhfcl',   'atom',' '       ,
     6             'drhfcl_ga','atom',' '       ,
     7             'dscrf',    'atom',' '       ,
     8             'duhfop',   'atom',' '       ,
     9             'uhfop',    'all', ' '       ,
     x             'drhfgvb',  'atom',' '       ,
     1             'denscf_ga','all', ' '       ,
     2             'indx2t',   'all',' '       ,
     3             'corbls',   'atom',' '       ,
     4             'corbld',   'atomw',' '       /
c
c...   which module ?
c
      do i=1,nmodules
         if (module.eq.modules(0,i)) go to 10
      end do
      error = module//' is not known in check_feature'
      call caserr(error)
c
10    imodule = i
c
      do icheck=1,ncheck
         if (modules(icheck,imodule).eq.' ') return 
         if (icheck.eq.1.and.ozora) then
c...        ZORA
            if (modules(icheck,imodule).eq.'all') return 
            if (modules(icheck,imodule).eq.'atom') then 
               if (nat_z.eq.0.or.oscalatz) then
                  error = module//
     1                  ' does not allow other than simple atomic ZORA'
                  call caserr(error)
               end if
            else if (modules(icheck,imodule).eq.'atomw') then 
               if (nat_z.eq.0.or.oscalatz) then
                  error = module//
     1                  ' does not allow other than simple atomic ZORA'
                  write(iwr,1)  error
1                 format(1x,/1x,75(1h*),/1x,a80,/1x,75(1H*))
               end if
            else if (modules(icheck,imodule).eq.'none') then 
               if (ozora) then
                  error = module//' does not allow for the use of ZORA'
                  if (ozora) call caserr(error)
               end if
            end if
         end if
      end do
c
      return
      end
