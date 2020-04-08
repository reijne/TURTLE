c ******************************************************
c ****** ROUTINES for input and handling of LOCAL and
c ****** non-LOCAL Pseudo-potentials.
c ******************************************************
**==inpseu.f
      subroutine inpseu(iblibb)
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
      character *1 zatm
      common/ecp2/ clp(400),alp(400),nlp(400),kfirst(maxat,6),
     1             klast(maxat,6),lmax(maxat),lpskip(maxat)
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
      character * 132 char1,char2,char1c
      common/workc/char1,char2,char1c
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
      character *8 zrunt, zdd3, zdd4
      character *4 yrunt, ydd
      common/direc/zrunt(50), yrunt(50), ydd(220), zdd3(70) ,zdd4(50)
c
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
      common/bufb/zlib(204),ztagl(10)
      common/blkin/chglib(204),maxlib(204),iblib(204)
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
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
      dimension iblibb(*)
      dimension inda(2),zatm(2)
      dimension ztemp(100)
      dimension zproj(5), zecp(8)
      data m204/204/
      data zproj/ '|s><s|', '|p><p|', '|d><d|', '|f><f|', '|g><g|'/
      data zend /'end'/
      data yrm1/'nonl'/
      data zform/'nwchem'/
      data zecp / 'sbkjc','cep','lanl','lanl2','crenbl',
     +            'crenbs','strlc','strsc' /
      indi = .false.
      diff = 0.0d0
      jrec = jrec - 1
      call inpa(ztext)
      if (ztext.ne.'ecp') then
         call inpa4(ytext)
         if (ytext.eq.yrm1) then
c
c --- nonlocal pseudopotential ---
c
            lpseud = 2
            call setstc(maxat,' ',zpseud)
c
c === library specification
c
            if (jump.gt.2) then
               call inpa4(ytext)
               call inpi(iblki)
               if (iblki.gt.0) iblkpl = iblki
               call filec(ytext,iblkpl,numlib,irep)
               if (irep.ne.0)
     +             call caserr2('invalid library attributes specified')
            end if
            write (iwr,6010) yed(numlib) , iblkpl
            write (iwr,6020)
            nav = lenwrd()
c           mword = 1760 + 8/nav
            m408 = 204 + 408/nav
            call rdchr(zlib,m204,iblkpl,numlib)
            call reads(chglib,m408,numlib)
 20         call input
            call inpa(ztext)
            ytext = ytrunc(ztext)
            i = locatc(ydd(101),limit2,ytext)
            if (i.ne.0) then
               write (iwr,6030)
               do 30 i = 1 , nat
                  if (zpseud(i).ne.' ') then
                     write (iwr,6040) i , zaname(i) , czan(i) , c(1,i) ,
     +                                c(2,i) , c(3,i) , zpseud(i)
                  end if
 30            continue
               go to 130
            else
               numps = locatc(zlib,m204,ztext)
               if (numps.eq.0) call caserr2(
     +         'attempt to locate unknown pseudopotential in library')
               chg = chglib(numps)
c
c --- ztext is the pseudopotential label
c
               jump1 = jump - 1
               do 40 i = 1 , jump1
                  call inpa(ztemp(i))
 40            continue
c
               ktyp = 0
               do 60 i = 1 , jump1
                  do 50 j = 1 , nat
                     if (ztemp(i).eq.zaname(j)) then
                        ktyp = ktyp + 1
c -- sum difference in charges (nonpseudo minus pseudo)
                        diff = diff + (czan(j)-chg)
                        czan(j) = chg
                        zpseud(j) = ztext
                     end if
 50               continue
 60            continue
               if (ktyp.ne.0) go to 20
               go to 100
            end if
         else if (ytext.eq.'ecp') then
         else if (ytext.ne.' ') then
            if (opg_root()) then
               write(iwr,6600)
            endif
            call caserr2('invalid pseudo potential directive')
         end if
      end if
c
c ---- ecp potential
c
      ind = 0
      iecp = 0
      lpseud = 1
      call setsto(maxat,1,lpskip)
      write (iwr,6050)
 70   call input
      oinpu = .false.
      call inpa(ztext)
      ytext = ytrunc(ztext)
      i = locatc(ydd(101),limit2,ytext)
      if (i.ne.0) go to 130
      i = locatc(ydd(1),limit1,ytext)
      if (i.ne.0) go to 130
c
chvd  check for an optional "end" directive to terminate the ECP input
c
      if (ztext.eq.zend) then
        call input
        call inpa(ztext)
        go to 130
      endif
c
c     REPULSION
c
c     Directive controlling core - core repulsion calculation together
c     with ecp calculation (optimization with this is also possible).
c     The energy expression added for a pair of atoms is:
c
c        E_ij = d_ij * exp(-eta_ij*Rij)
c
      if (ytext.eq.'repu') then
       if (dabs(diff).gt.1.0d-10) then
         write(iwr,*)'*** The core-core repulsion pairs must be defined'
         write(iwr,*)'*** in the ECP input block but before any actual'
         write(iwr,*)'*** ECPs, i.e.:'
         write(iwr,*)'***'
         write(iwr,*)'***    PSEUDO ECP'
         write(iwr,*)'***    REPULSION'
         write(iwr,*)'***    <element A> <element B>'
         write(iwr,*)'***    <prefactor> <exponent>'
         write(iwr,*)'***    ...'
         write(iwr,*)'***    END'
         write(iwr,*)'***    <ECP data>'
         write(iwr,*)'***    ...'
         call caserr('define core-core repulsions before ECPs')
       endif
       indi = .true.
       do i=1,nat
        ichgat(i) = nint(czan(i))
       enddo
       npairs = 0
  72   call input
       call inpa(ztext)
       write (iwr,73) char1
  73   format (a80)
       if (ztext.eq.zend) goto 70
       npairs = npairs + 1
       if (npairs.gt.mxpairs) then
         write(iwr,*)' *** The current configuration allows for at most'
         write(iwr,*)' *** ',mxpairs,' core-core repulsion pairs'
         write(iwr,*)' *** Please adjust mxpairs in common/repp/ and',
     +               ' recompile'
         call caserr("too many core-core repulsion pairs")
       endif
       do i=1,2
        zatm(1) = char1(istrt(i):istrt(i))
        if (inumb(i).eq.1) then
           zatm(2) = ' '
        else
           zatm(2)=char1(istrt(i)+1:istrt(i)+1)
        endif
        inda(i)=isubst(zatm)
        ofound = .false.
        do j = 1, nat
          ofound = ofound.or.ichgat(j).eq.inda(i)
        enddo
        if (.not.ofound) then
          write(iwr,*)'*** WARNING: no atom of element ',
     +                char1(istrt(i):istrt(i)+inumb(i)-1),
     +                ' present in molecule'
        endif
       enddo
       indx(npairs) = inda(1)
       jndx(npairs) = inda(2)
       read (ird,*) da,etaa
       write (iwr,75) da,etaa
  75   format (7x,'d',15x,'eta'/2f15.6)
       d(npairs)=da
       eta(npairs)=etaa
       goto 72
      endif
c
c     end of reading core - core repulsion information
c
c     ztext is pseupotential label
c
      if (ztext.eq.'card'.or. ztext.eq.'cards') then
       oinpu = .true.
       call inpa(ztext1)
        if(ztext1.eq.'nwchem') then
         onwchem = .true.
         jump1 = jump-2
        else
         jump1 = jump - 1
         onwchem = .false.
         jrec = jrec -1
        endif
c
      else
c
      oinpu = .false.
c
c     now must deal with more flexible input to cover multiple ecp
c     options over and above simply cep / hw. Must cover alternative
c     specifications:
c     1. atom_type CEP_type [atom-tags from input geometry]
c     2. atom_type  [atom-tags from input geometry] .. old HW default
c     3. atomtype[cep]  [atom-tags from input geometry] .. old CEP specification
c
c     is the second field a recognised ECP type
       call inpa(ztext1)
       iecpt = locatc(zecp,8,ztext1)
       if(iecpt.ne.0) then
c      recognised ECP, so ztext is atom-type
       iecpx = iecpt
       jump1 = jump - 2
       else
c      not recognised as ECP, so must be atom-tag with
c      ECP embedded in first field
c      to make sure, check that it is indeed an atom tag
       icheck = locatc(zaname,nat,ztext1)
       iecpx = 0
       if(icheck.eq.0) then
        call caserr2('attempt to reference invalid ECP or atom')
       endif
       jrec = jrec - 1
       jump1 = jump -1
       endif
c
      endif
c
c     now process list of atom tags ..
c
      natecp = 0
      do loop = 1 , jump1
         call inpa(zchar)
         do j = 1 , nat
            if (zchar.eq.zaname(j)) then
               natecp = natecp + 1
               iblibb(natecp) = j
            end if
         enddo
      enddo
      if (natecp.gt.0 .and. natecp.le.maxat) go to 110
 100  write (iwr,6060)
      call caserr2('error detected in ecp or pseudo potential data')
 110  iecp = iecp + 1
      if (oinpu) then
         call ecpinp(diff,ind,iblibb,natecp,onwchem)
         iecpt = 0
      else
         call ecplib(ztext,iecpt,diff,ind,iblibb,natecp)
      end if
c     print out ecp  data
      if (oinpu.or.iecpt.eq.0) then
       write (iwr,6070) iecp
      else
       write(iwr,6075) zecp(iecpt), iecp
      endif
      write (iwr,6080) (iblibb(i),zaname(iblibb(i)),i=1,natecp)
      do l = 1, natecp
       if (czan(iblibb(l)).lt.0.0) then
c...           ghost
          diff = diff + czan(iblibb(l))
          czan(iblibb(l)) = 0.0
          lpskip(iblibb(l)) = 1
       else
          ipseud(iblibb(l)) = iecpt
       end if
      enddo
      iat = iblibb(1)
      write (iwr,6090) czan(iat) , lmax(iat)
      lmx = lmax(iat) + 1
      do l = 1 , lmx
         ind1 = kfirst(iat,l)
         ind2 = klast(iat,l)
         if (l.eq.1) then
            write (iwr,6100) (nlp(k),clp(k),alp(k),k=ind1,ind2)
         else
            write (iwr,6110) zproj(l-1) ,
     +                       (nlp(k),clp(k),alp(k),k=ind1,ind2)
         end if
      enddo
      go to 70
c
c  adjust number of electrons to take account of pseudopotential
c
 130  jrec = jrec - 1
      idiff = nint(diff)
c
c     Adjust the effective number of electrons for the electrons represented
c     by the ECP. For ECP's with odd numbers of electrons (for example with
c     actinides) we assume that the ECP corresponds to a doublet and remove
c     electrons accordingly from the alpha and beta electrons.
c
      ne = ne - idiff
      nh = idiff/2
      na = na - (idiff - nh)
      nb = nb - nh
c
      write (iwr,6120) ne , na , nb
      write (iwr,6130)
      if (ne.lt.0 . or. na.lt.0  .or. nb.lt.0) then
       write(iwr,6125)
       call caserr2('negative number of effective electrons')
      endif
c
      return
 6010 format (/1x,104('=')///15x,
     +        'n o n l o c a l    p s e u d o p o t e n t i a l   ',
     +        ' o p e r a t o r'/15x,
     +        '---------------------------------------------------',
     +        '-----------------'//15x,
     +        'pseudopotential library restored from ',a4,' at block ',
     +        i5)
 6020 format (/20x,'potential (s) used:'/)
 6030 format (//18x,'nuclear',27x,'coordinates'/9x,'atom',5x,'charge',
     +        13x,'x',19x,'y',19x,'z',10x,'pseudopotential'/)
 6040 format (1x,i4,5x,a4,5x,f5.2,2x,3f20.10,7x,a4)
 6050 format (/
     +     3x,'****************************************************** '/
     +     3x,'*              EFFECTIVE CORE POTENTIALS             * '/
     +     3x,'****************************************************** ')
 6060 format (//1x,'invalid ecp or pseudo-tag specified'/)
 6070 format (/4x,'*** ECP operator # ',i3,' ***'/
     +         4x,'***  for atoms ...')
 6075 format (/4x,'*** ECP operator (internal library: ',a6,') # ',
     +         i3,' ***'/
     +         4x,'*** for atoms ...')
 6080 format (5x,5(i3,'-',a8))
 6090 format (5x,'with effective core charge ',f5.1,'  and lmax ',i2)
 6100 format (5x,'term l+1',5x,i3,2f15.8,20(/,18x,i3,2f15.8))
 6110 format (5x,'term ',a6,2x,i3,2f15.8,20(/,18x,i3,2f15.8))
 6120 format (/5x,'Effective total no. of electrons =',i5/
     +         5x,'Effective no. of occupied orbitals (alpha) =',i5/
     +         5x,'Effective no. of occupied orbitals (beta)  =',i5)
 6125 format (/5x,
     + '**** invalid no. of electrons: check tags and charges ****')
 6130 format (
     +     3x,'******************************************************'/)
 6600 format(' *** Your input was not recognised.',
     +     /,' *** The only valid pseudo potential cards are:',
     +     /,' ***     ECP',
     +     /,' ***     PSEUDO',
     +     /,' ***     PSEUDO ECP',
     +     /,' ***     PSEUDO NONLOCAL [<library file spec>]',//)
      end
**==inpseu2.f
      subroutine inpseu2(core)
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
      common/ecp2/clp(400),alp(400),nlp(400),kfirst(maxat,6),
     1            klast(maxat,6),lmax(maxat),lpskip(maxat)
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
      integer ncoef1,ncoef2,j1len,j2len,llim,nlim,ntlim,j4len
      common /ecpdim/ ncoef1,ncoef2,j1len,j2len,llim,nlim,ntlim,j4len
c
      dimension core(*)
      data m18/18/
c
c   more rigorous checking of requested ecps with data written
c   to section isect(479) for subsequent retrieval in ecp integrals
c   and gradients.
c
c   check for the maximum basis type
c   1 is added here since baschk returns the real angular mometum
c
      call baschk(maxang)
      maxang=maxang+1
c
      call derchk(numder)
      nlim = maxang + numder
c
c     check for the maximum potential type: s=0, p=1, ..., g=4 is ok
c
      llim = 0
      do iat=1,nat
        llim = max(lmax(iat), llim)
      enddo
      nlim = max(nlim,llim)
c
c     tests to see if we can do this run...
c
      if(maxang.ge.6) then
       write(iwr,9100)
       call caserr2('ecp integrals limited to spdfg basis functions')
      end if
c
      if(llim.ge.5) then
       write(iwr,9110)
       call caserr2('ecp integrals limited to spdfg core potentials')
      end if
      if(numder.gt.0  .and.  maxang.eq.6) then
         do 855 jsh=1,nshell
            if(ktype(jsh).ne.4) go to 855
            do 850 iat=1,nat
               if(katom(jsh).eq.iat) go to 850
               if(lpskip(iat).eq.0) then
                  write(iwr,9120)
                  call caserr2('limitations in ecp gradient integrals')
               end if
  850       continue
  855    continue
      end if
c
c  specify the maximal type of basis function
c
      ntlim=(nlim*(nlim+1)*(nlim+2))/6
c
c      initialize ncoef1 based on nlim
c
      ncoef1 = 8520
      ncoef2 = 3424
      if (nlim.eq.5) then
        ncoef2 = 10555
        ncoef1 = 71660
      else if (nlim.eq.6) then
        ncoef1 = 280000
        ncoef2 = 28940
      else if (nlim.eq.7) then
        if (numder.eq.1) then
          ncoef1 = 284868
          ncoef2 = 28940
        else
          ncoef1 = 892584
          ncoef2 = 60382
        end if
      end if
c
c        look at -eccod1/2- to see that -jfst1/2- really do need 1 extra
c
      j1len = (ntlim*ntlim+ntlim)/2 + 1
      j2len = (llim-1)*(llim+1)*ntlim+ntlim + 1
      j4len = ntlim+ntlim*(llim-1)*(llim+1)
c
      nav = lenwrd()
c
      ldcf1  = 0
      ljin   = ldcf1  +  ncoef1
      llb1   = ljin   + (j1len-1)/nav+1
      ldcf2  = llb1   + (ncoef1*9)/nav
      llj2   = ldcf2  +  ncoef2
      llb2   = llj2   + (j2len-1)/nav+1
      lfpqr  = llb2   + (ncoef2*6)/nav
      lzlm   = lfpqr  + 15625
      llmf   = lzlm   + 584
      llmx   = llmf   + 124/nav
      llmy   = llmx   + 584/nav
      llmz   = llmy   + 584/nav
      last   = llmz   + 584/nav
c
      ldcf1 = igmem_alloc(last)
c
      ljin   = ldcf1  +  ncoef1
      llb1   = ljin   + (j1len-1)/nav+1
      ldcf2  = llb1   + (ncoef1*9)/nav
      llj2   = ldcf2  +  ncoef2
      llb2   = llj2   + (j2len-1)/nav+1
      lfpqr  = llb2   + (ncoef2*6)/nav
      lzlm   = lfpqr  + 15625
      llmf   = lzlm   + 584
      llmx   = llmf   + 124/nav
      llmy   = llmx   + 584/nav
      llmz   = llmy   + 584/nav
      last   = llmz   + 584/nav
c
c     ----- encode the ecp angular integrals -----
c
      call eccodr(core(ldcf1),core(ljin),core(llb1),core(ldcf2),
     +            core(llj2),core(llb2),core(lfpqr),core(lzlm),
     +            core(llmf),core(llmx),core(llmy),core(llmz),numder)
c
      ldaf91 = (j1len-1)/nav+1 + (9*ncoef1-1)/nav+1
      ldaf93 = (j2len-1)/nav+1 + (6*ncoef2)/nav
c
      len479 = lensec(ncoef1) + lensec(ldaf91) +
     +         lensec(ncoef2) + lensec(ldaf93)
      call secput(isect(479),m18,len479,ibl479)
c
      call wrt3(core(ldcf1),ncoef1,ibl479,idaf)
      call wrt3s(core(ljin),ldaf91,idaf)
      call wrt3s(core(ldcf2),ncoef2,idaf)
      call wrt3s(core(llj2),ldaf93,idaf)
c
      call gmem_free(ldcf1)
c
      return
c
 9120 format(/1x,'the ecp gradient integrals are limited to spd basis,'/
     *       1x,'functions, except for the following special case:'/
     *       1x,'if the only atom with f functions is also'/
     *       1x,'   the only atom with a core potential,'/
     *       1x,'then spdf ecp gradients are possible.')
 9100 format(/1x,'the ecp integrals are limited to spdfg basis',
     *          ' functions.')
 9110 format(/1x,'the ecp integrals are limited to spdfg core',
     *          ' potentials.')
      end
**==baschk.f
      subroutine baschk(lmax)
c
      implicit real*8 (a-h,o-z)
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
c     return the highest angular momentum present in the basis.
c     note that ktype=1,2,3,4,5 means s, p(l), d, f, g function.
c
      kang = 0
      do  n=1,nshell
       if(ktype(n).gt.kang) kang = ktype(n)
      enddo
      lmax = kang-1
      return
      end
**==derchk.f
      subroutine derchk(maxder)
c
      implicit real*8 (a-h,o-z)
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
      character *8 zrunt, zdd3, zdd4
      character *4 yrunt, ydd
      common/direc/zrunt(50), yrunt(50), ydd(220), zdd3(70) ,zdd4(50)
c
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
c
c     return maximum derivative to be computed by this run
c
      maxder = 0
c
      irun = locatc(zrunt,mxtask,zruntp)
c
      if (irun.ne.0) then
c
      go to (30 ,30 ,30 ,40 ,40 ,40 ,40 ,40 ,30 ,30,
     +       30 ,30 ,30 ,30 ,30 ,30 ,30 ,50 ,30 ,30,
     +       30 ,50 ,50 ,50 ,50 ,50 ,50 ,50 ,30, 30,
     +       50 ,50 ,50 ,50 ,30 ,30 ,30 ,30 ,30 ,30,
     +       30 ,30 ,30 ,30 ,30 ,30 ,30 ,30 ,30 ,30) , irun
c
c     gradient required
 40   maxder = 1
      go to 30
c     assume fully analytic hessian
 50   maxder = 2
c
      endif
c
 30   return
      end
**==ecpinp.f
      subroutine ecpinp (diff,ind,list,nlib,onwchem)
      implicit real*8  (a-h,p-w), integer (i-n), logical(o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *8 textatm, text, ztag1
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
      common/ecp2 / clp(400),alp(400),nlp(400),kfirst(maxat,6),
     1              klast(maxat,6),lmax(maxat),lpskip(maxat)
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
      dimension list(*)
c
      iat = list(1)
      ind2 = ind
c
      if(onwchem) then
c      process intro/label line
       call input
       call input
       call inpa(textatm)
       call inpa(text)
       if(text.ne.'nelec') call caserr2('invalid format')
       call inpi(nelcor)
c      write(6,*)'nelcor = ', i2
c      modify number of electrons and nuclear charge
       if (nelcor.ge.czan(iat))
     +     call caserr2('invalid number of core electrons')
       fnel = dble(nelcor)
       czan(iat) = czan(iat) - fnel
       lmx0 = 0
       ind1 = ind
1210   call input
       call inpa(ztag1)
       if(ztag1.eq.'end') go to 1200
       if(ztag1.eq.textatm)  then
c       write(6,*) ztag1
        lmx0 = lmx0 + 1
        kfirst(iat,lmx0) = ind1 + 1
        klast (iat,lmx0) = ind1
       else
        jrec = jrec -1
        ind1 = ind1 + 1
        call inpi(nlp(ind1))
        call inpf(alp(ind1))
        if (alp(ind1).lt.0.0d0) then
           call caserr2('negative exponent in ecp')
        endif
        call inpf(clp(ind1))
        klast (iat,lmx0) = klast (iat,lmx0) + 1
c       write(6,*) nlp(ind1), alp(ind1), clp(ind1)
       endif
       go to 1210
 1200  lmx = lmx0 -1
       lmax(iat) = lmx
       if(lmx.eq.0 .or. lmx.ge.6)
     +    call caserr2('invalid data in ecp input')
       ind = ind1
      else
       call input
       call inpi(lmx)
       call inpi(nelcor)
c      modify number of electrons and nuclear charge
       if (nelcor.ge.czan(iat))
     +     call caserr2('invalid number of core electrons')
       fnel = dble(nelcor)
       czan(iat) = czan(iat) - fnel
       if (lmx.eq.0 .or. lmx.ge.6)
     +     call caserr2('invalid data in ecp input')
       lmax(iat) = lmx
       lmx = lmx + 1
       do l = 1 , lmx
          call input
          call inpi(nterm)
          ind1 = ind2 + 1
          ind2 = ind2 + nterm
          kfirst(iat,l) = ind1
          klast(iat,l) = ind2
          do i = ind1 , ind2
             call input
             call inpi(nlp(i))
             call inpf(clp(i))
             call inpf(alp(i))
             if (alp(i).lt.0.0d0) then
                call caserr2('negative exponent in ecp')
             endif
          enddo
       enddo
       ind = ind2
      endif
c
      diff = diff + fnel
      lpskip(iat) = 0
c
c
c === load on equivalent centres
c
      do 50 j = 2 , nlib
         kk = list(j)
         czan(kk) = czan(kk) - fnel
         lpskip(kk) = 0
         lmax(kk) = lmax(iat)
         lmx = lmax(iat) + 1
         do 40 l = 1 , lmx
            kfirst(kk,l) = kfirst(iat,l)
            klast(kk,l) = klast(iat,l)
 40      continue
         diff = diff + fnel
 50   continue
c
      return
      end
**==ecplib.f
      subroutine ecplib (zinp,iecpt,diff,ind,list,nlib)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character * 8 (z)
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
      common/ecp2 / clp(400),alp(400),nlp(400),kfirst(maxat,6),
     1              klast(maxat,6),lmax(maxat),lpskip(maxat)
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
      dimension list(*)
      dimension nt(6), n(10,6), cc(10,6), a(10,6)
c
      dimension ztit(103)
      data maxtyp/103/
      data ztit/
     $         'h ', 'he',
     $         'li', 'be', 'b ', 'c ', 'n ', 'o ', 'f ', 'ne',
     $         'na', 'mg', 'al', 'si', 'p ', 's ', 'cl', 'ar',
     $         'k ', 'ca',
     $                     'sc', 'ti', 'v ', 'cr', 'mn',
     $                     'fe', 'co', 'ni', 'cu', 'zn',
     $                     'ga', 'ge', 'as', 'se', 'br', 'kr',
     $ 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag','cd',
     $ 'in','sn','sb','te','i ','xe','cs','ba','la','ce','pr','nd',
     $ 'pm','sm','eu','gd','tb','dy','ho','er','tm','yb','lu','hf',
     $ 'ta','w ','re','os','ir','pt','au','hg','tl','pb','bi','po',
     $ 'at','rn','fr','ra','ac','th','pa','u ','np','pu','am','cm',
     $ 'bk','cf','es','fm','md','no','lw'   /
c
c   this routine contains a library of ecp data taken from
c   various literature sources.
c   the ecps are taken from ....
c   ******************************************************************
c   for li ......ne
c   CEP / SBKJC basis sets
c   and na ......ar from stevens,basch,krauss, jcp 81(1984) 6026
c   for k to rn,
c   w.stevens, p.g.jasien, m.krauss, h.basch, can.j.chem. 70, 612(1992)
c   for the lanthanides,
c   t.r.cundari, w.j.stevens, j.chem.phys. 98, 5555(1993)
c   ******************************************************************
c   and sc,......hg from hay,wadt, jcp 82 (1985) 270
c   and na,......bi from hay,wadt, jcp 82 (1985) 284
c   ******************************************************************
c   CRENBL ECP   -  Christiansen et al. Large Orbital Basis, Small Core Pot.                  
c  Elements  Primitives                         References                        
c  H:       (4s)            T. H. Dunning, Jr. and P. J. Hay, Methods of       
c                           Electronic Structure Theory, Vol. 3, H. F.         
c                           Schaefer III, Ed. Plenum Press (1977).             
c  Li - Ne: (4s,4p)         L. F. Pacios and P. A. Christiansen, J. Chem. Phys.
c                                                                              
c                           82, (1985) 2664.                                   
c  Na - Mg: (6s,4p)                                                            
c  Al - Ar: (4s,4p)                                                            
c  K - Ca:  (5s,4p)         M. M. Hurley et al. J. Chem. Phys. 84, 6840 (1986) 
c  Sc - Zn: (7s,6p,6d)                                                         
c  Ga - Kr: (3s,3p,4d)                                                         
c  Rb - Sr: (5s,5p)         L. A. LaJohn et al. J. Chem. Phys., 87, 2812 (1987)
c                                                                              
c  Y - Cd:  (5s,5p,4d)                                                         
c  In - I:  (3s,3p,4d)                                                         
c  Xe:      (3s,3p,4d)      M. M. Hurley et al., J. Chem. Phys. 84, (1986) 6840
c  .                                                                           
c  Cs:      (5s,5p,4d)      R. B. Ross, W. C. Ermler, P. A. Christiansen et al.
c                                                                              
c                           J. Chem. Phys. 93, 6654 (1990).                    
c  Ba:      (5s,5p,4d)      R.B. ROSS, W.C. ERMLER, P.A. CHRISTIANSEN,         
c                           ET AL. SUB.TO J. CHEM. PHYS.                       
c  La:      (5s,5p,4d)      R.B. ROSS, W.C. ERMLER, P.A. CHRISTIANSEN          
c                           ET AL. J. CHEM. PHYS. 93,(1990)6654.               
c  Ce - Lu: (6s,6p,6d,6f)   R.B. Ross, W.C. Ermler, S. Das, To be published    
c  Hf - Hg: (5s,5p,4d)      R.B. ROSS, W.C. ERMLER, P.A. CHRISTIANSEN          
c                           ET AL. J. CHEM. PHYS. 93,(1990)6654.               
c  Ti - Rn: (3s,3p,4d)                                                         
c  Fr - Ra: (5s,5p,4d)      W.C. ERMLER, R.B. ROSS, P.A. CHRISTIANSEN,         
c                           INT. J. QUANT. CHEM. 40,(1991)829.                 
c  Ac - Pu: (5s,5p,4d,4f)                                                      
c  Am - Uub:  (0s,2p,6d,5f) C.S. NASH, B.E. BURSTEN, W.C. ERMLER,              
c                                         J. CHEM. PHYS. 1997                  
c  Uut - Uuo: (0s,3p,6d,5f)                                                    
**                                                                             
c  These ECPs are sometimes referred to as shape consistent because they          
c  maintain the shape of the atomic orbitals in the valence region.
c   ******************************************************************
c  CRENBS ECP     -   Ermler and Co-workers Small Basis Sets for Use with 
c                     Averaged Relativistic, Large Core ECPs
c  Elements     Contraction                             References                
c  Sc - Co: (4s)             M.M. Hurley et al., J. Chem. Phys., 84, 6840 (1986)  
c  Ni:      (4s)             L.A. LaJohn et al., J. CHEM. PHYS., 87, (1987) 2812. 
c  Cu - Zn: (4s,5d)    -> [1s,1d]    M.M. Hurley et al., J. Chem. Phys., 84, 6840 
c                                    (1986)                                       
c  Y - Cd:  (3s,3p,4d) -> [1s,1p,1d] L.A. LaJohn et al., J. CHEM. PHYS., 87,      
c                                    2812 (1987).                                 
c  La, Hf-Hg: (3s,3p,4d)     R.B. Ross, W.C. Ermler, P.A. CHRISTIANSEN ET AL. J.  
c                            CHEM. PHYS. 93,(1990)6654.                           
c  Ti - Rn: (3s,3p)                                                               
c  Rf - Uut:(0s,5p,6d)       C.S. NASH, B.E. BURSTEN, W.C. ERMLER, J. CHEM.       
c                            PHYS. 1997                                           
c  Uuq - Uuo: (0s,6p,6d)                                                          
c  ** e.g.  for copper only the 4s and 3d orbitals are included in the        
c  valence space.                                                                 
c  These ECPs are sometimes referred to as shape consistent because they          
c  maintain the shape of the atomic orbitals in the valence region.
c   ******************************************************************
c   RSC ECP   Stuttgart Relativistic, Small Core ECP Basis Set                    
c                         Stuttgart                                              
c              Core         Name     Primitives        Contractions              
c     K     : 10[Ne]       ECP10MWB (7s,6p)          -> [5s,4p]                  
c    Ca     : 10           ECP10MWB (6s,6p,5d)       -> [4s,4p,2d]               
c    Sc - Zn: 10           ECP10MDF (8s,7p,6d,1f)    -> [6s,5p,3d,1f]            
c    Rb     : 28[Ar+3d]    ECP28MWB (7s,6p)          -> [5s,4p]                  
c    Sr     : 28           ECP28MWB (6s,6p,5d)       -> [4s,4p,2d]               
c    Y  - Cd: 28           ECP28MHF (8s,7p,6d)       -> [6s,5p,3d]               
c    Cs     : 46[Kr+4d]    ECP46MWB (7s,6p)          -> [5s,4p]                  
c    Ba     : 46           ECP46MWB (6s,6p,5d,1f)    -> [3s,3p,2d,1f]            
c    Ce - Ho: 28[Ar+3d]    ECP28MWB (12s,11p,9d,8f)  -> [5s,5p,4d,3f]            
c    Er - Yb: 28           ECP28MWB (12s,10p,8d,8f)  -> [5s,5p,4d,3f]            
c    Hf - Hg: 60[Kr+4df]   ECP60MWB (8s,7p,6d)       -> [6s,5p,3d]               
c    Ac - Lr: 60           ECP60MWB (12s,11p,10d,8f) -> [8s,7p,6d,4f]            
c    Db     : 92[Xe+4f+5d] ECP92MWB (8s,7p,6d,2f,1g) -> [6s,5p,4d,2f,1g]         
c
c    K      : A. Bergner, M. Dolg, W. Kuechle, H. Stoll, H. Preuss, 
c             Mol. Phys. 80,  1431 (1993).                                                          
c    Ca     : M. Kaupp, P. v. R. Schleyer, H. Stoll, H. Preuss, 
c             J. Chem. Phys. 94, 1360 (1991).                                                          
c    Rf - Db: M. Dolg, H. Stoll, H. Preuss, R.M. Pitzer, 
c             J. Phys. Chem. 97, 5852   (1993).                                                              
c   ******************************************************************
c   RLC ECP   Stuttgart Relativistic, Large Core ECP Basis Set                    
c                        Stuttgart                                                 
c             Core         Name     Primitives        Contractions                 
c   Li - Be:  2[He]       ECP2SDF  (4s,4p)          -> [2s,2p]                     
c    B - N :  2           ECP2MWB  (4s,4p)          -> [2s,2p]                     
c    O - F :  2           ECP2MWB  (4s,5p)          -> [2s,3p]                     
c   Ne     :  2           ECP2MWB  (7s,7p,3d,1f)    -> [4s,4p,3d,1f]               
c   Na - Mg: 10[Ne]       ECP10SDF (4s,4p)          -> [2s,2p]                     
c   Al - P : 10           ECP10MWB (4s,4p)          -> [2s,2p]                     
c    S - Cl: 10           ECP10MWB (4s,5p)          -> [2s,3p]                     
c   Ar     : 10           ECP10MWB (6s,6p,3d,1f)    -> [4s,4p,3d,1f]               
c    K - Ca: 18[Ar]       ECP18SDF (4s,4p)          -> [2s,2p]                     
c   Zn     : 28[Ar+3d]    ECP28MWB (4s,2p)          -> [3s,2p]                     
c   Ga - As: 28           ECP28MWB (4s,2p)          -> [3s,2p]                     
c   Se - Br: 28           ECP28MWB (4s,5p)          -> [2s,3p]                     
c   Kr     : 28           ECP28MWB (6s,6p,3d,1f)    -> [4s,4p,3d,1f]               
c   Rb - Sr: 36[Kr]       ECP36SDF (4s,4p)          -> [2s,2p]                     
c   In - Sb: 46[Kr+4d]    ECP46MWB (4s,4p)          -> [2s,2p]                     
c   Te - I : 46           ECP46MWB (4s,5p)          -> [2s,3p]                     
c   Xe     : 46           ECP46MWB (6s,6p,3d,1f)    -> [4s,4p,3d,1f]               
c   Cs - Ba: 54[Xe]       ECP54SDF (4s,4p)          -> [2s,2p]                     
c   Hg - Rn: 78[Xe+4f+5d] ECP78MWB (4s,4p,1d)       -> [2s,2p,1d]                  
c   Ac - Lr: 78           ECP78MWB (8s,8p,6d,5f,2g) -> [5s,5p,4d,3f,2g]            
c   
c   Li - Be: P. Fuentealba, H. Preuss, H. Stoll, L. v. Szentpaly, 
c   Na       Chem. Phys. Lett.  89, 418 (1982).
c                                                                                  
c   B  - Ne: A. Bergner, M. Dolg, W. Kuechle, H. Stoll, H. Preuss, 
c            Mol. Phys. 80,  1431 (1993).                                                          
c   Mg     : P. Fuentealba, L. v. Szentpaly, H. Preuss, H. Stoll, 
c            J. Phys. B 18,  1287 (1985).                                                          
c   Al     : G. Igel-Mann, H. Stoll, H. Preuss, 
c            Mol. Phys. 65, 1321 (1988).        
c   Hg - Rn: W. Kuechle, M. Dolg, H. Stoll, H. Preuss, 
c            Mol. Phys. 74, 1245 (1991). 
c   Ac - Lr: W. Kuechle, to be published                                           
c   ******************************************************************
      ocep = .false.
      olanl = .false.
      ocrenbl = .false.
      ocrenbs = .false.
      ostrlc = .false.
      ostrsc = .false.
      ohw = .false.
c
c     process iecpt flag
c
      iecpx = iecpt
      if (iecpx.eq.0) then
       iat = index(zinp,'cep')
       if(iat.gt.0) then
        ocep = .true.
        zinp (iat:iat+2) = '   '
        iecpt = 1
       else
        ohw=.true.
        iecpt = 3
        iat = index(zinp,'lanl2')
        if(iat.gt.0) then
         olanl = .true.
         zinp(iat:iat+4) = '      '
         iecpt = 4
        endif
       endif
c
      else if (iecpx.eq.1.or.iecpx.eq.2) then
        ocep = .true.
        iecpt = 1
      else if (iecpx.eq.3.or.iecpx.eq.4) then
        ohw = .true.
        if(iecpx.eq.4) olanl = .true.
      else if (iecpx.eq.5) then
        ocrenbl = .true.
      else if (iecpx.eq.6) then
        ocrenbs = .true.
      else if (iecpx.eq.7) then
        ostrlc = .true.
      else if (iecpx.eq.8) then
        ostrsc = .true.
      endif
c
      ztag = ztit(jsubst(zinp))
      ktype = locatc(ztit,maxtyp,ztag)
      if (ktype.le.0) then
         write (iwr,6010) zinp
         call caserr2('no library available for requested ecp')
      end if
c
      if (ohw) then
c
c     hay,wadt ...
c
       if (ktype.le.18) then
          call ecplb0(ocep,ztag,lmx,nelcor,nt,n,cc,a,olanl)
c
       else if (ktype.le.36) then
          call ecplb1(ztag,lmx,nelcor,nt,n,cc,a,olanl)
c
       else if (ktype.le.54) then
          call ecplb2(ztag,lmx,nelcor,nt,n,cc,a,olanl)
c
       else if (ktype.le.83) then
          call ecplb3(ztag,lmx,nelcor,nt,n,cc,a,olanl)
c
       else if (ktype.ge.92.or.ktype.le.94) then
          call ecplb4(ztag,lmx,nelcor,nt,n,cc,a,olanl)
       else
          call caserr2('no HW ecp for requested element')
       end if
       go to 100
      endif
c
      if (ocep) then
c
c     sjkb CEP potentials
c
       if (ktype.le.18) then
          call ecplb0(ocep,ztag,lmx,nelcor,nt,n,cc,a,olanl)
c
       else if(ktype.ge.19.and.ktype.le.36) then
          if(ktype.ge.21  .and.  ktype.le.31) then
             call sbkpm1(ztag,lmx,nelcor,nt,n,cc,a)
          else
             call sbkp4(ztag,lmx,nelcor,nt,n,cc,a)
          end if
c
       else if(ktype.ge.37.and.ktype.le.54) then
          if(ktype.ge.39  .and.  ktype.le.49) then
             call sbkpm2(ztag,lmx,nelcor,nt,n,cc,a)
          else
             call sbkp5(ztag,lmx,nelcor,nt,n,cc,a)
          end if
c
       else if(ktype.ge.55.and.ktype.le.86) then
         if(ktype.ge.58.and.ktype.le.71) then
            call sbklan(ztag,lmx,nelcor,nt,n,cc,a)
          else if(ktype.ge.57  .and.  ktype.le.81) then
             call sbkpm3(ztag,lmx,nelcor,nt,n,cc,a)
          else
             call sbkp6(ztag,lmx,nelcor,nt,n,cc,a)
          end if
c
       else
           call caserr2('no ecp for requested element')
       end if
       go to 100
      endif
c
      if (ocrenbl) then
c
c     CRENBL
c
       if (ktype.le.18) then
          call crenbl0(ztag,lmx,nelcor,nt,n,cc,a)
c
       else if(ktype.ge.19.and.ktype.le.36) then
          call crenbl1(ztag,lmx,nelcor,nt,n,cc,a)
c
       else if(ktype.ge.37.and.ktype.le.54) then
          call crenbl2(ztag,lmx,nelcor,nt,n,cc,a)
c
       else if(ktype.ge.55.and.ktype.le.71) then
          call crenbl3(ztag,lmx,nelcor,nt,n,cc,a)
c
       else if(ktype.ge.72.and.ktype.le.86) then
          call crenbl4(ztag,lmx,nelcor,nt,n,cc,a)
c
       else if(ktype.ge.87.and.ktype.le.103) then
          call crenbl5(ztag,lmx,nelcor,nt,n,cc,a)
c
       else
          call caserr2('no CRENBL ecp for requested element')
       end if
       go to 100
      endif
c
      if (ocrenbs) then
c
c     CRENBS
c
       if(ktype.ge.21  .and.  ktype.le.30) then
          call crenbs0(ztag,lmx,nelcor,nt,n,cc,a)
c
       else if (ktype.ge.39  .and.  ktype.le.48) then
          call crenbs1(ztag,lmx,nelcor,nt,n,cc,a)
c
       else if (ktype.eq.57. or .
     +         (ktype.ge.72.and.ktype.le.80) .or.
     +         (ktype.ge.82.and.ktype.le.86)   ) then
          call crenbs2(ztag,lmx,nelcor,nt,n,cc,a)
c
       else
         call caserr2('no CRENBS ecp for requested element')
       endif
       go to 100
      endif
c
      if (ostrlc) then
c
c     STRLC
c
       if (ktype.le.18) then
          call strlc0(ztag,lmx,nelcor,nt,n,cc,a)
c
       else if(ktype.ge.19.and.ktype.le.36) then
          call strlc1(ztag,lmx,nelcor,nt,n,cc,a)
c
       else if(ktype.ge.37.and.ktype.le.54) then
          call strlc2(ztag,lmx,nelcor,nt,n,cc,a)
c
       else if(ktype.ge.55.and.ktype.le.71) then
          call strlc3(ztag,lmx,nelcor,nt,n,cc,a)
c
       else if(ktype.ge.72.and.ktype.le.86) then
          call strlc4(ztag,lmx,nelcor,nt,n,cc,a)
c
       else if(ktype.ge.87.and.ktype.le.103) then
          call strlc5(ztag,lmx,nelcor,nt,n,cc,a)
c
       else
          call caserr2('no Stuttgart RLC ecp for requested element')
       end if
       go to 100
      endif
c
      if (ostrsc) then
c
c     STRSC
c
       if(ktype.ge.19.and.ktype.le.36) then
          call strsc1(ztag,lmx,nelcor,nt,n,cc,a)
c
       else if(ktype.ge.37.and.ktype.le.54) then
          call strsc2(ztag,lmx,nelcor,nt,n,cc,a)
c
       else if(ktype.ge.55.and.ktype.le.71) then
          call strsc3(ztag,lmx,nelcor,nt,n,cc,a)
c
       else if(ktype.ge.72.and.ktype.le.86) then
          call strsc4(ztag,lmx,nelcor,nt,n,cc,a)
c
       else if(ktype.ge.87.and.ktype.le.103) then
          call strsc5(ztag,lmx,nelcor,nt,n,cc,a)
c
       else
          call caserr2('no Stuttgart RSC ecp for requested element')
       end if
       go to 100
      endif
c
c     load ecp data
c
100   iat = list(1)
      lpskip(iat) = 0
      fnel = dble(nelcor)
      czan(iat) = czan(iat) - fnel
      diff = diff + fnel
      lmax(iat) = lmx
      lmx = lmx + 1
      do l = 1 , lmx
         nterm = nt(l)
         ind1 = ind + 1
         ind2 = ind + nterm
         kfirst(iat,l) = ind1
         klast(iat,l) = ind2
         kk = 0
         do k = ind1 , ind2
            kk = kk + 1
            nlp(k) = n(kk,l)
            clp(k) = cc(kk,l)
            alp(k) = a(kk,l)
         enddo
         ind = ind2
      enddo
c
c === load on equivalent centres
c
      do j = 2 , nlib
         kk = list(j)
         czan(kk) = czan(kk) - fnel
         lpskip(kk) = 0
         lmax(kk) = lmax(iat)
         lmx = lmax(iat) + 1
         do l = 1 , lmx
            kfirst(kk,l) = kfirst(iat,l)
            klast(kk,l) = klast(iat,l)
         enddo
         diff = diff + fnel
      enddo
c
      return
 6010 format (1x,'*** data not found for ecp type ',a8)
      end
**==ecplb0.f
      subroutine ecplb0 (ocep,zinp,lmx,nelcor,nt,n,cc,a,olanl)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character * 8 (z)
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
      dimension nt(6), n(10,6), cc(10,6), a(10,6)
c
      dimension ztit(18)
      data maxtyp/18/
      data ztit/'h','he',
     *          'li','be','b ','c ','n ','o ','f ','ne',
     *          'na','mg','al','si','p ','s ','cl','ar'/
c
c     this routine contains a library of ecp data taken from
c     various literature sources.
c     the ecps for...
c     na,......ar taken from hay,wadt, jcp 82 (1985) 270
c     cep potentials
c     li,......ar taken from stevens,basch,krauss jcp 81 (1984) 6026
c
      do 20 ityp = 1 , maxtyp
         if (zinp.eq.ztit(ityp)) go to 40
 20   continue
 30   write (iwr,6010) zinp
      call caserr2('no library available for requested ecp')
 40   if(ityp.le.10.and.olanl) go to 310
      go to (30,30,50,60,70,80,90,100,110,120,
     +       130,150,170,190,210,230,250,270) , ityp
c
c     lithium  basch (Li_SBKJC_VDZ_ECP)
c
 50   cc(1,1) = -0.73063d0
      cc(1,2) = 1.80131d0
      cc(2,2) = 3.54971d0
      a(1,1) = 1.34306d0
      a(1,2) = 0.61284d0
      a(2,2) = 1.64881d0
c
      go to 280
c
c     berillium  basch (Be_SBKJC_VDZ_ECP)
c
 60   cc(1,1) = -0.78970d0
      cc(1,2) = 1.84229d0
      cc(2,2) = 6.81169d0
      a(1,1) = 2.79347d0
      a(1,2) = 1.19874d0
      a(2,2) = 3.20738d0
c
      go to 280
c
c     boron  basch (B_SBKJC_VDZ_ECP)
c
 70   cc(1,1) = -0.88355d0
      cc(1,2) = 1.89572d0
      cc(2,2) = 10.97543d0
      a(1,1) = 5.64622d0
      a(1,2) = 1.92061d0
      a(2,2) = 5.55177d0
c
      go to 280
c
c     carbon basch (C_SBKJC_VDZ_ECP)
c
 80   cc(1,1) = -0.89371d0
      cc(1,2) = 1.92926d0
      cc(2,2) = 14.88199d0
      a(1,1) = 8.56468d0
      a(1,2) = 2.81497d0
      a(2,2) = 8.11296d0
c
      go to 280
c
c     nitrogen basch (N_SBKJC_VDZ_ECP)
c
 90   cc(1,1) = -0.91212d0
      cc(1,2) = 1.93565d0
c *** INCORRECT entry
c     cc(2,2) = 21.73359d0
      cc(2,2) = 21.73355d0
      a(1,1) = 11.99686d0
      a(1,2) = 3.83895d0
      a(2,2) = 11.73247d0
c
      go to 280
c
c     oxygen basch (O_SBKJC_VDZ_ECP)
c
 100  cc(1,1) = -0.92550d0
      cc(1,2) = 1.96069d0
      cc(2,2) = 29.13442d0
      a(1,1) = 16.11718d0
      a(1,2) = 5.05348d0
      a(2,2) = 15.95333d0
c
      go to 280
c
c     fluorine basch (F_SBKJC_VDZ_ECP)
c
 110  cc(1,1) = -0.93258d0
      cc(1,2) = 2.03649d0
      cc(2,2) = 27.86279d0
      a(1,1) = 20.73959d0
      a(1,2) = 6.60488d0
      a(2,2) = 18.24092d0
c
      go to 280
c
c     neon basch (Ne_SBKJC_VDZ_ECP)
c
 120  cc(1,1) = -0.94098d0
      cc(1,2) = 2.09952d0
      cc(2,2) = 31.46718d0
      a(1,1) = 26.01666d0
      a(1,2) = 8.28047d0
      a(2,2) = 22.44102d0
c
      go to 280
c
c     sodium
c
 130  if (ocep) then
c     basch (Na_SBKJC_VDZ_ECP)
      cc(1,1) = -2.38446d0
      cc(1,2) = 6.23415d0
      cc(2,2) = 9.08374d0
      cc(1,3) = 3.23971d0
      cc(2,3) = 2.53514d0
      a(1,1) = 0.90009d0
      a(1,2) = 5.37232d0
      a(2,2) = 1.11959d0
      a(1,3) = 1.29158d0
      a(2,3) = 0.65791d0
c
      else
c     hay  (Na_LANL2DZ_ECP)
      cc(1,1) = -10.0000000d0
      cc(2,1) = -47.4902024d0
      cc(3,1) = -17.2283007d0
      cc(4,1) = -6.0637782d0
      cc(5,1) = -0.7299393d0
      cc(1,2) = 3.0000000d0
      cc(2,2) = 36.2847626d0
      cc(3,2) = 72.9304880d0
      cc(4,2) = 23.8401151d0
      cc(5,2) = 6.0123861d0
      cc(1,3) = 5.0000000d0
      cc(2,3) = 117.4495683d0
      cc(3,3) = 423.3986704d0
      cc(4,3) = 109.3247297d0
      cc(5,3) = 31.3701656d0
      cc(6,3) = 7.1241813d0
      a(1,1) = 175.5502590d0
      a(2,1) = 35.0516791d0
      a(3,1) = 7.9060270d0
      a(4,1) = 2.3365719d0
      a(5,1) = 0.7799867d0
      a(1,2) = 243.3605846d0
      a(2,2) = 41.5764759d0
      a(3,2) = 13.2649167d0
      a(4,2) = 3.6797165d0
      a(5,2) = 0.9764209d0
      a(1,3) = 1257.2650682d0
      a(2,3) = 189.6248810d0
      a(3,3) = 54.5247759d0
      a(4,3) = 13.7449955d0
      a(5,3) = 3.6813579d0
      a(6,3) = 0.9461106d0
      endif
      go to 290
c
c     magnesium
c
 150  if (ocep) then
c     basch (Mg_SBKJC_VDZ_ECP)
      cc(1,1) = -2.77296d0
      cc(1,2) = 6.41029d0
      cc(2,2) = 14.12727d0
      cc(1,3) = 3.24762d0
      cc(2,3) = 4.33012d0
      a(1,1) = 1.39220d0
      a(1,2) = 6.96291d0
      a(2,2) = 1.55458d0
      a(1,3) = 1.89177d0
      a(2,3) = 0.98235d0
c
      else
c
c     hay (Mg_LANL2DZ_ECP)
      cc(1,1) = -10.0000000d0
      cc(2,1) = -55.8993968d0
      cc(3,1) = -20.1391957d0
      cc(4,1) = -7.0679107d0
      cc(5,1) = -0.8133109d0
      cc(1,2) = 3.0000000d0
      cc(2,2) = 44.0075660d0
      cc(3,2) = 107.3861344d0
      cc(4,2) = 35.8289088d0
      cc(5,2) = 10.1143435d0
      cc(1,3) = 5.0000000d0
      cc(2,3) = 117.1053672d0
      cc(3,3) = 420.5972073d0
      cc(4,3) = 107.6122959d0
      cc(5,3) = 29.1002576d0
      cc(6,3) = 7.0875570d0
      a(1,1) = 237.5484804d0
      a(2,1) = 47.7520367d0
      a(3,1) = 10.7837852d0
      a(4,1) = 3.1992124d0
      a(5,1) = 1.0636953d0
      a(1,2) = 348.3008631d0
      a(2,2) = 59.4680133d0
      a(3,2) = 19.0767582d0
      a(4,2) = 5.2965613d0
      a(5,2) = 1.3867373d0
      a(1,3) = 1256.8739085d0
      a(2,3) = 189.8608839d0
      a(3,3) = 54.6949631d0
      a(4,3) = 13.8990137d0
      a(5,3) = 3.9597181d0
      a(6,3) = 1.2552787d0
      endif
      go to 290
c
c     aluminium
c
 170  if (ocep) then
c     basch (Al_SBKJC_VDZ_ECP)
      cc(1,1) = -3.03055d0
      cc(1,2) = 6.04650d0
      cc(2,2) = 18.87509d0
      cc(1,3) = 3.29465d0
      cc(2,3) = 6.87029d0
      a(1,1) = 1.95559d0
      a(1,2) = 7.78858d0
      a(2,2) = 1.99025d0
      a(1,3) = 2.83146d0
      a(2,3) = 1.38479d0
c     basch
      else
c     hay (Al_LANL2DZ_ECP)
      cc(1,1) = -10.0000000d0
      cc(2,1) = -63.8079837d0
      cc(3,1) = -22.8972174d0
      cc(4,1) = -8.0063232d0
      cc(5,1) = -0.8829345d0
      cc(1,2) = 3.0000000d0
      cc(2,2) = 50.9705682d0
      cc(3,2) = 143.8716460d0
      cc(4,2) = 48.0055753d0
      cc(5,2) = 14.0114180d0
      cc(1,3) = 5.0000000d0
      cc(2,3) = 93.0595424d0
      cc(3,3) = 266.7686329d0
      cc(4,3) = 70.0725805d0
      cc(5,3) = 17.2290839d0
      cc(6,3) = 0.7105331d0
      a(1,1) = 304.7291926d0
      a(2,1) = 61.5299768d0
      a(3,1) = 13.9259006d0
      a(4,1) = 4.1463626d0
      a(5,1) = 1.3715443d0
      a(1,2) = 467.8437756d0
      a(2,2) = 79.3992216d0
      a(3,2) = 25.4035967d0
      a(4,2) = 6.9954696d0
      a(5,2) = 1.7860129d0
c *** incorrect ENTRY 
c *** a(1,3) = 776.2727190d0
      a(1,3) = 776.2717190d0
      a(2,3) = 118.4992254d0
      a(3,3) = 34.4107276d0
      a(4,3) = 8.7859563d0
      a(5,3) = 2.3406228d0
      a(6,3) = 0.7398386d0
      endif
      go to 290
c
c     silicon
c
 190  if (ocep) then
c     basch (Si_SBKJC_VDZ_ECP)
      cc(1,1) = -2.98240d0
      cc(1,2) = 6.00558d0
      cc(2,2) = 23.77642d0
      cc(1,3) = 3.20345d0
      cc(2,3) = 9.86438d0
      a(1,1) = 2.462720d0
      a(1,2) = 9.46338d0
      a(2,2) = 2.444080d0
      a(1,3) = 3.82585d0
      a(2,3) = 1.83026d0
      else
c     hay (Si_LANL2DZ_ECP)
      cc(1,1) = -10.0000000d0
      cc(2,1) = -84.9236087d0
      cc(3,1) = -30.3299410d0
      cc(4,1) = -12.1049046d0
      cc(5,1) = -1.8945408d0
      cc(1,2) = 3.0000000d0
      cc(2,2) = 60.5206833d0
      cc(3,2) = 201.3086137d0
      cc(4,2) = 65.9399745d0
      cc(5,2) = 19.0300789d0
      cc(1,3) = 5.0000000d0
      cc(2,3) = 6.6413817d0
      cc(3,3) = 247.5972416d0
      cc(4,3) = 129.3715380d0
      cc(5,3) = 47.4617107d0
      cc(6,3) = 11.7376574d0
      a(1,1) = 505.3137693d0
      a(2,1) = 103.2221026d0
      a(3,1) = 23.4569248d0
      a(4,1) = 6.7505693d0
      a(5,1) = 2.1603140d0
      a(1,2) = 689.4910719d0
      a(2,2) = 114.1728510d0
      a(3,2) = 35.7424336d0
      a(4,2) = 9.4529586d0
      a(5,2) = 2.2543590d0
      a(1,3) = 88.9379355d0
      a(2,3) = 76.7773538d0
      a(3,3) = 56.1480987d0
      a(4,3) = 21.1874014d0
      a(5,3) = 6.8277276d0
      a(6,3) = 2.1001192d0
      endif
      go to 290
c
c     phosphorus
c
 210  if (ocep) then
c     basch  (P_SBKJC_VDZ_ECP)
      cc(1,1) = -3.16361d0
      cc(1,2) = 6.20760d0
      cc(2,2) = 29.40024d0
      cc(1,3) = 3.24084d0
      cc(2,3) = 13.10778d0
      a(1,1) = 3.189450d0
      a(1,2) = 11.79328d0
      a(2,2) = 2.953300d0
      a(1,3) = 4.96368d0
      a(2,3) = 2.31865d0
      else
c
c     hay (P_LANL2DZ_ECP)
c
      cc(1,1) = -10.0000000d0
      cc(2,1) = -79.4864658d0
      cc(3,1) = -28.3668251d0
      cc(4,1) = -9.8577589d0
      cc(5,1) = -1.0163783d0
      cc(1,2) = 3.0000000d0
      cc(2,2) = 12.9104154d0
      cc(3,2) = 150.0250298d0
      cc(4,2) = 71.7083146d0
      cc(5,2) = 23.0397012d0
      cc(1,3) = 5.0000000d0
      cc(2,3) = 6.3446507d0
      cc(3,3) = 198.5585104d0
      cc(4,3) = 111.1470820d0
      cc(5,3) = 40.3944144d0
      cc(6,3) = 6.4483233d0
      a(1,1) = 462.1211423d0
      a(2,1) = 93.6863701d0
      a(3,1) = 21.2349094d0
      a(4,1) = 6.3388415d0
      a(5,1) = 2.0620684d0
      a(1,2) = 78.0831823d0
      a(2,2) = 58.9576810d0
      a(3,2) = 36.0571255d0
      a(4,2) = 11.2464453d0
      a(5,2) = 2.6757561d0
      a(1,3) = 75.1617880d0
      a(2,3) = 57.4544041d0
      a(3,3) = 47.9481748d0
      a(4,3) = 18.4588360d0
      a(5,3) = 5.9414190d0
      a(6,3) = 1.8487507d0
      endif
      go to 290
c
c     sulfur 
c
 230  if (ocep) then
c     basch (S_SBKJC_VDZ_ECP)
      cc(1,1) = -3.29839d0
      cc(1,2) = 6.44496d0
      cc(2,2) = 35.62708d0
      cc(1,3) = 3.35704d0
      cc(2,3) = 17.54237d0
      a(1,1) = 3.98930d0
      a(1,2) = 15.58064d0
      a(2,2) = 3.50176d0
      a(1,3) = 6.88540d0
      a(2,3) = 2.90690d0
      else
c     hay (S_LANL2DZ_ECP)
      cc(1,1) = -10.0000000d0
      cc(2,1) = -85.3593846d0
      cc(3,1) = -30.4513290d0
      cc(4,1) = -10.3745886d0
      cc(5,1) = -0.9899295d0
      cc(1,2) = 3.0000000d0
      cc(2,2) = 10.6284036d0
      cc(3,2) = 223.6360469d0
      cc(4,2) = 93.6460845d0
      cc(5,2) = 28.7609065d0
      cc(1,3) = 5.0000000d0
      cc(2,3) = 6.0969842d0
      cc(3,3) = 285.4425500d0
      cc(4,3) = 147.1448413d0
      cc(5,3) = 53.6569778d0
      cc(6,3) = 8.9249559d0
      a(1,1) = 532.6685222d0
      a(2,1) = 108.1342248d0
      a(3,1) = 24.5697664d0
      a(4,1) = 7.3702438d0
      a(5,1) = 2.3712569d0
      a(1,2) = 106.3176781d0
      a(2,2) = 100.8245833d0
      a(3,2) = 53.5858472d0
      a(4,2) = 15.3706332d0
      a(5,2) = 3.1778402d0
      a(1,3) = 101.9709185d0
      a(2,3) = 93.2808973d0
      a(3,3) = 65.1431772d0
      a(4,3) = 24.6347440d0
      a(5,3) = 7.8120535d0
      a(6,3) = 2.3112730d0
      endif
      go to 290
c
c     chlorine
c
 250  if (ocep) then
c     basch (Cl_SBKJC_VDZ_ECP)
      cc(1,1) = -3.40738d0
      cc(1,2) = 6.50966d0
      cc(2,2) = 42.27785d0
      cc(1,3) = 3.42860d0
      cc(2,3) = 22.15256d0
      a(1,1) = 4.87483d0
      a(1,2) = 17.00367d0
      a(2,2) = 4.10380d0
      a(1,3) = 8.90029d0
      a(2,3) = 3.52648d0
      else
c     hay (Cl_LANL2DZ_ECP)
      cc(1,1) = -10.0000000d0
      cc(2,1) = 66.272917d0
      cc(3,1) = -28.968595d0
      cc(4,1) = -12.866337d0
      cc(5,1) = -1.710217d0
      cc(1,2) = 3.0000000d0
      cc(2,2) = 12.852851d0
      cc(3,2) = 275.672398d0
      cc(4,2) = 115.677712d0
      cc(5,2) = 35.060609d0
      cc(1,3) = 5.d0
      cc(2,3) = 7.479486d0
      cc(3,3) = 613.032000d0
      cc(4,3) = 280.800685d0
      cc(5,3) = 107.878824d0
      cc(6,3) = 15.343956d0
      a(1,1) = 94.8130d0
      a(2,1) = 165.6440d0
      a(3,1) = 30.8317d0
      a(4,1) = 10.5841d0
      a(5,1) = 3.7704d0
      a(1,2) = 128.8391d0
      a(2,2) = 120.3786d0
      a(3,2) = 63.5622d0
      a(4,2) = 18.0695d0
      a(5,2) = 3.8142d0
      a(1,3) = 216.5263d0
      a(2,3) = 46.5723d0
      a(3,3) = 147.4685d0
      a(4,3) = 48.9869d0
      a(5,3) = 13.2096d0
      a(6,3) = 3.1831d0
      endif
      go to 290
c
c     argon
c
 270  if (ocep) then
c     basch (Ar_SBKJC_VDZ_ECP)
      cc(1,1) = -3.48806d0
      cc(1,2) = 7.67855d0
      cc(2,2) = 49.92637d0
      cc(1,3) = 3.52500d0
      cc(2,3) = 27.41686d0
      a(1,1) = 5.82303d0
      a(1,2) = 25.93841d0
      a(2,2) = 4.73780d0
      a(1,3) = 11.41178d0
      a(2,3) = 4.20476d0
      else
c     hay (Ar_LANL2DZ_ECP)
      cc(1,1) = -10.0000000d0
      cc(2,1) = -99.0606669d0
      cc(3,1) = -35.2711767d0
      cc(4,1) = -11.8151947d0
      cc(5,1) = -1.0453382d0
      cc(1,2) = 3.0000000d0
      cc(2,2) = 11.2465621d0
      cc(3,2) = 250.7412812d0
      cc(4,2) = 139.1606543d0
      cc(5,2) = 41.2897981d0
      cc(1,3) = 5.0000000d0
      cc(2,3) = 5.2787358d0
      cc(3,3) = 631.7135166d0
      cc(4,3) = 305.1649809d0
      cc(5,3) = 106.4807615d0
      cc(6,3) = 17.0215765d0
      a(1,1) = 711.5242175d0
      a(2,1) = 144.6708689d0
      a(3,1) = 32.9246992d0
      a(4,1) = 9.9103877d0
      a(5,1) = 3.1328926d0
      a(1,2) = 152.8553033d0
      a(2,2) = 82.1424792d0
      a(3,2) = 83.7154800d0
      a(4,2) = 23.8557161d0
      a(5,2) = 4.3128823d0
      a(1,3) = 212.5573747d0
      a(2,3) = 266.4884712d0
      a(3,3) = 139.9222477d0
      a(4,3) = 50.1097659d0
      a(5,3) = 14.5537276d0
      a(6,3) = 3.4828623d0
      endif
      go to 290
 280  lmx = 1
      nelcor = 2
      nt(1) = 1
      nt(2) = 2
      n(1,1) = 1
      n(1,2) = 0
      n(2,2) = 2
      go to 300
 290  lmx = 2
      nelcor = 10
      if (ocep) then
c     basch
       nt(1) = 1
       nt(2) = 2
       nt(3) = 2
       n(1,1) = 1
       n(1,2) = 0
       n(2,2) = 2
       n(1,3) = 0
       n(2,3) = 2
      else
       nt(1) = 5
       nt(2) = 5
       nt(3) = 6
       n(1,1) = 1
       n(2,1) = 2
       n(3,1) = 2
       n(4,1) = 2
       n(5,1) = 2
       n(1,2) = 0
       n(2,2) = 1
       n(3,2) = 2
       n(4,2) = 2
       n(5,2) = 2
       n(1,3) = 0
       n(2,3) = 1
       n(3,3) = 2
       n(4,3) = 2
       n(5,3) = 2
       n(6,3) = 2
      endif
 300  return
 310  write(iwr,6010) zinp
      call caserr2('no library available for requested LANL2 ecp')
 6010 format (1x,'*** data not found for ecp type ',a7)
      return
      end
**==ecplb1.f
      subroutine ecplb1 (zinp,lmx,nelcor,nt,n,cc,a,olanl)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character * 8 (z)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      dimension nt(6), n(10,6), cc(10,6), a(10,6)
c
      dimension ztit(18)
      data maxtyp/18/
      data ztit/'k','ca',
     *          'sc','ti','v ','cr','mn','fe',
     *          'co','ni','cu','zn','ga','ge','as','se','br','kr'/
c
c     this routine contains a library of ecp data taken from
c     various literature sources.
c     the ecps for
c     k,......kr are taken from hay,wadt, jcp 82 (1985) 270
      do 20 ityp = 1 , maxtyp
         if (zinp.eq.ztit(ityp)) go to 30
 20   continue
      write (iwr,6010) zinp
      call caserr2('no library available for requested ecp')
 30   go to (40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,
     +       200,210) , ityp
c     k  hay (K_LANL2DZ_ECP)
 40   lmx = 2
      nelcor = 18
      nt(1) = 5
      nt(2) = 5
      nt(3) = 5
      cc(1,1) = -18.0000000d0
      cc(2,1) = -63.7625207d0
      cc(3,1) = -15.8426437d0
      cc(4,1) = -5.0502393d0
      cc(5,1) = -0.8055766d0
      cc(1,2) = 3.0000000d0
      cc(2,2) = 22.9220976d0
      cc(3,2) = 32.8399986d0
      cc(4,2) = 15.5672939d0
      cc(5,2) = 5.2099815d0
      cc(1,3) = 5.0000000d0
      cc(2,3) = 14.0261987d0
      cc(3,3) = 8.4093418d0
      cc(4,3) = -4.3049715d0
      cc(5,3) = 5.0151840d0
      a(1,1) = 107.3810865d0
      a(2,1) = 19.9413498d0
      a(3,1) = 4.2697333d0
      a(4,1) = 1.0926291d0
      a(5,1) = 0.3594045d0
      a(1,2) = 27.1076482d0
      a(2,2) = 10.9915385d0
      a(3,2) = 5.1237645d0
      a(4,2) = 1.7372696d0
      a(5,2) = 0.4999159d0
      a(1,3) = 6.6617383d0
      a(2,3) = 2.1231073d0
      a(3,3) = 0.8214406d0
      a(4,3) = 0.4472818d0
      a(5,3) = 0.3568265d0
      n(1,1) = 1
      n(2,1) = 2
      n(3,1) = 2
      n(4,1) = 2
      n(5,1) = 2
      n(1,2) = 0
      n(2,2) = 1
      n(3,2) = 2
      n(4,2) = 2
      n(5,2) = 2
      n(1,3) = 0
      n(2,3) = 1
      n(3,3) = 2
      n(4,3) = 2
      n(5,3) = 2
      go to 220
c     ca hay (Ca_LANL2DZ_ECP)
 50   lmx = 2
      nelcor = 18
      nt(1) = 5
      nt(2) = 5
      nt(3) = 4
      cc(1,1) = -18.0000000d0
      cc(2,1) = -72.3806397d0
      cc(3,1) = -18.5797971d0
      cc(4,1) = -5.8911769d0
      cc(5,1) = -0.8339620d0
      cc(1,2) = 3.0000000d0
      cc(2,2) = 19.7879237d0
      cc(3,2) = 44.4890340d0
      cc(4,2) = 22.6988644d0
      cc(5,2) = 7.5593017d0
      cc(1,3) = 5.0000000d0
      cc(2,3) = 14.5671920d0
      cc(3,3) = 11.8395160d0
      cc(4,3) = 3.8072337d0
      a(1,1) = 135.0709504d0
      a(2,1) = 25.1922367d0
      a(3,1) = 5.4940376d0
      a(4,1) = 1.4151874d0
      a(5,1) = 0.4286428d0
      a(1,2) = 29.3419490d0
      a(2,2) = 14.7433814d0
      a(3,2) = 7.9251170d0
      a(4,2) = 2.5306550d0
      a(5,2) = 0.6334168d0
      a(1,3) = 8.0909308d0
      a(2,3) = 3.1248522d0
      a(3,3) = 1.4462879d0
      a(4,3) = 0.4385463d0
      n(1,1) = 1
      n(2,1) = 2
      n(3,1) = 2
      n(4,1) = 2
      n(5,1) = 2
      n(1,2) = 0
      n(2,2) = 1
      n(3,2) = 2
      n(4,2) = 2
      n(5,2) = 2
      n(1,3) = 0
      n(2,3) = 1
      n(3,3) = 2
      n(4,3) = 2
      go to 220
c
c     scandium  hay
c
 60   if(olanl) then
       cc( 1, 1) =        -10.00000000d0
       cc( 2, 1) =        -48.96955530d0
       cc( 3, 1) =         -8.70660440d0
       a ( 1, 1) =        237.31665280d0
       a ( 2, 1) =         42.64890040d0
       a ( 3, 1) =         10.62112290d0
       n ( 1, 1) =   1
       n ( 2, 1) =   2
       n ( 3, 1) =   2
       nt( 1)    =   3
       cc( 1, 2) =          3.00000000d0
       cc( 2, 2) =         17.05845590d0
       cc( 3, 2) =        167.95954120d0
       cc( 4, 2) =        179.83821130d0
       cc( 5, 2) =       -121.15231440d0
       a ( 1, 2) =         69.25523240d0
       a ( 2, 2) =         58.89534160d0
       a ( 3, 2) =         27.26374590d0
       a ( 4, 2) =          5.47373040d0
       a ( 5, 2) =          5.15236070d0
       n ( 1, 2) =   0
       n ( 2, 2) =   1
       n ( 3, 2) =   2
       n ( 4, 2) =   2
       n ( 5, 2) =   2
       nt( 2)    =   5
       cc( 1, 3) =          5.00000000d0
       cc( 2, 3) =          6.10881950d0
       cc( 3, 3) =        192.93849140d0
       cc( 4, 3) =        127.22747550d0
       cc( 5, 3) =         68.09459740d0
       cc( 6, 3) =        -41.64691080d0
       a ( 1, 3) =         77.25129700d0
       a ( 2, 3) =         44.45271550d0
       a ( 3, 3) =         53.68547520d0
       a ( 4, 3) =         18.94211700d0
       a ( 5, 3) =          4.79705730d0
       a ( 6, 3) =          4.41429250d0
       n ( 1, 3) =   0
       n ( 2, 3) =   1
       n ( 3, 3) =   2
       n ( 4, 3) =   2
       n ( 5, 3) =   2
       n ( 6, 3) =   2
       nt( 3)    =   6
       nelcor    =  10
       lmx       =   2
      else
       lmx = 3
       nelcor = 18
       nt(1) = 5
       nt(2) = 4
       nt(3) = 5
       nt(4) = 3
       cc(1,1) = -18.0000000d0
       cc(2,1) = -88.5127611d0
       cc(3,1) = -22.9161516d0
       cc(4,1) = -8.3516042d0
       cc(5,1) = -1.3015863d0
       cc(1,2) = 3.0000000d0
       cc(2,2) = 12.8461779d0
       cc(3,2) = 17.2046793d0
       cc(4,2) = -10.0019149d0
       cc(1,3) = 5.0000000d0
       cc(2,3) = 13.4339672d0
       cc(3,3) = 73.9358072d0
       cc(4,3) = 29.3227672d0
       cc(5,3) = 3.7948832d0
       cc(1,4) = -1.6459835d0
       cc(2,4) = 0.6046650d0
       cc(3,4) = -0.3017156d0
c
       a(1,1) = 192.5195835d0
       a(2,1) = 36.1035075d0
       a(3,1) = 8.2101758d0
       a(4,1) = 2.0455730d0
       a(5,1) = 0.7105354d0
       a(1,2) = 9.6336398d0
       a(2,2) = 2.5946676d0
       a(3,2) = 0.5448709d0
       a(4,2) = 0.4812194d0
       a(1,3) = 22.4571432d0
       a(2,3) = 18.0561916d0
       a(3,3) = 9.7416660d0
       a(4,3) = 2.5427651d0
       a(5,3) = 0.4450257d0
       a(1,4) = 7.6371109d0
       a(2,4) = 1.6990508d0
       a(3,4) = 0.3542965d0
       n(1,1) = 1
       n(2,1) = 2
       n(3,1) = 2
       n(4,1) = 2
       n(5,1) = 2
       n(1,2) = 0
       n(2,2) = 1
       n(3,2) = 2
       n(4,2) = 2
       n(1,3) = 0
       n(2,3) = 1
       n(3,3) = 2
       n(4,3) = 2
       n(5,3) = 2
       n(1,4) = 2
       n(2,4) = 2
       n(3,4) = 2
      endif
      go to 220
c
c     titanium  hay
c
 70   if(olanl) then
       cc( 1, 1) =        -10.00000000d0
       cc( 2, 1) =        -51.84278160d0
       cc( 3, 1) =         -9.14291450d0
       a ( 1, 1) =        265.32639090d0
       a ( 2, 1) =         47.76878150d0
       a ( 3, 1) =         11.89033340d0
       n ( 1, 1) =   1
       n ( 2, 1) =   2
       n ( 3, 1) =   2
       nt( 1)    =   3
       cc( 1, 2) =          3.00000000d0
       cc( 2, 2) =         19.48255790d0
       cc( 3, 2) =        207.33492790d0
       cc( 4, 2) =        235.67445010d0
       cc( 5, 2) =       -166.87843870d0
       a ( 1, 2) =         81.47306960d0
       a ( 2, 2) =         72.64967240d0
       a ( 3, 2) =         31.81282130d0
       a ( 4, 2) =          6.16644680d0
       a ( 5, 2) =          5.82683470d0
       n ( 1, 2) =   0
       n ( 2, 2) =   1
       n ( 3, 2) =   2
       n ( 4, 2) =   2
       n ( 5, 2) =   2
       nt( 2)    =   5
       cc( 1, 3) =          5.00000000d0
       cc( 2, 3) =          5.53488220d0
       cc( 3, 3) =        177.84193840d0
       cc( 4, 3) =        107.42071530d0
       cc( 5, 3) =        -71.90659020d0
       a ( 1, 3) =         50.29669430d0
       a ( 2, 3) =         63.50897540d0
       a ( 3, 3) =         26.09960840d0
       a ( 4, 3) =          5.60225730d0
       a ( 5, 3) =          5.21710690d0
       n ( 1, 3) =   0
       n ( 2, 3) =   1
       n ( 3, 3) =   2
       n ( 4, 3) =   2
       n ( 5, 3) =   2
       nt( 3)    =   5
       nelcor    =  10
       lmx       =   2
      else
       lmx = 3
       nelcor = 18
       nt(1) = 5
       nt(2) = 4
       nt(3) = 5
       nt(4) = 3
       cc(1,1) = -18.0000000d0
       cc(2,1) = -94.1338170d0
       cc(3,1) = -24.4645079d0
       cc(4,1) = -8.8573627d0
       cc(5,1) = -1.2827508d0
       cc(1,2) = 3.0000000d0
       cc(2,2) = 14.1432597d0
       cc(3,2) = 16.8993983d0
       cc(4,2) = -8.5136295d0
       cc(1,3) = 5.0000000d0
       cc(2,3) = 15.5558118d0
       cc(3,3) = 29.0702900d0
       cc(4,3) = 29.8642273d0
       cc(5,3) = 4.3716011d0
       cc(1,4) = -1.6761296d0
       cc(2,4) = 0.6972213d0
       cc(3,4) = -0.3201236d0
c
       a(1,1) = 217.0090397d0
       a(2,1) = 40.7301667d0
       a(3,1) = 9.2567439d0
       a(4,1) = 2.3249337d0
       a(5,1) = 0.7950530d0
       a(1,2) = 10.0762958d0
       a(2,2) = 2.6789265d0
       a(3,2) = 0.6503944d0
       a(4,2) = 0.5646085d0
       a(1,3) = 18.1336424d0
       a(2,3) = 8.3111333d0
       a(3,3) = 8.7336319d0
       a(4,3) = 2.5966642d0
       a(5,3) = 0.5206680d0
       a(1,4) = 7.9748454d0
       a(2,4) = 2.2713857d0
       a(3,4) = 0.3889416d0
       n(1,1) = 1
       n(2,1) = 2
       n(3,1) = 2
       n(4,1) = 2
       n(5,1) = 2
       n(1,2) = 0
       n(2,2) = 1
       n(3,2) = 2
       n(4,2) = 2
       n(1,3) = 0
       n(2,3) = 1
       n(3,3) = 2
       n(4,3) = 2
       n(5,3) = 2
       n(1,4) = 2
       n(2,4) = 2
       n(3,4) = 2
      endif
      go to 220
c
c     vanadium  hay
c
 80   if(olanl) then
       cc( 1, 1) =        -10.00000000d0
       cc( 2, 1) =        -54.89098080d0
       cc( 3, 1) =         -9.70105860d0
       a ( 1, 1) =        296.44875820d0
       a ( 2, 1) =         53.50150520d0
       a ( 3, 1) =         13.32045100d0
       n ( 1, 1) =   1
       n ( 2, 1) =   2
       n ( 3, 1) =   2
       nt( 1)    =   3
       cc( 1, 2) =          3.00000000d0
       cc( 2, 2) =         22.53456890d0
       cc( 3, 2) =        216.56340090d0
       cc( 4, 2) =        252.68225380d0
       cc( 5, 2) =       -176.26716540d0
       a ( 1, 2) =         93.21936540d0
       a ( 2, 2) =         71.10295670d0
       a ( 3, 2) =         33.87669930d0
       a ( 4, 2) =          6.86912000d0
       a ( 5, 2) =          6.48047220d0
       n ( 1, 2) =   0
       n ( 2, 2) =   1
       n ( 3, 2) =   2
       n ( 4, 2) =   2
       n ( 5, 2) =   2
       nt( 2)    =   5
       cc( 1, 3) =          5.00000000d0
       cc( 2, 3) =          5.70035550d0
       cc( 3, 3) =        182.90542740d0
       cc( 4, 3) =        117.27562260d0
       cc( 5, 3) =        -78.58291340d0
       a ( 1, 3) =         52.09539770d0
       a ( 2, 3) =         65.67784110d0
       a ( 3, 3) =         27.47037640d0
       a ( 4, 3) =          6.23520200d0
       a ( 5, 3) =          5.81399090d0
       n ( 1, 3) =   0
       n ( 2, 3) =   1
       n ( 3, 3) =   2
       n ( 4, 3) =   2
       n ( 5, 3) =   2
       nt( 3)    =   5
       nelcor    =  10
       lmx       =   2
      else
       lmx = 3
       nelcor = 18
       nt(1) = 5
       nt(2) = 4
       nt(3) = 5
       nt(4) = 2
       cc(1,1) = -18.0000000d0
       cc(2,1) = -99.0218605d0
       cc(3,1) = -25.6876997d0
       cc(4,1) = -9.2754278d0
       cc(5,1) = -1.2436490d0
       cc(1,2) = 3.0000000d0
       cc(2,2) = 15.5584631d0
       cc(3,2) = 19.1278354d0
       cc(4,2) = -9.8715206d0
       cc(1,3) = 5.0000000d0
       cc(2,3) = 16.4607976d0
       cc(3,3) = 26.5741922d0
       cc(4,3) = 31.0252324d0
       cc(5,3) = 4.8480079d0
       cc(1,4) = 0.3127658d0
       cc(2,4) = -0.3566007d0
c
       a(1,1) = 239.9870839d0
       a(2,1) = 45.0835633d0
       a(3,1) = 10.2253062d0
       a(4,1) = 2.5936435d0
       a(5,1) = 0.8734825d0
       a(1,2) = 10.6463564d0
       a(2,2) = 2.7790805d0
       a(3,2) = 0.7269532d0
       a(4,2) = 0.6351304d0
       a(1,3) = 18.1251385d0
       a(2,3) = 8.2306080d0
       a(3,3) = 8.6771989d0
       a(4,3) = 2.7167571d0
       a(5,3) = 0.5961611d0
       a(1,4) = 1.7439119d0
       a(2,4) = 0.4363769d0
       n(1,1) = 1
       n(2,1) = 2
       n(3,1) = 2
       n(4,1) = 2
       n(5,1) = 2
       n(1,2) = 0
       n(2,2) = 1
       n(3,2) = 2
       n(4,2) = 2
       n(1,3) = 0
       n(2,3) = 1
       n(3,3) = 2
       n(4,3) = 2
       n(5,3) = 2
       n(1,4) = 2
       n(2,4) = 2
      endif
      go to 220
c
c     chromium  hay
c
 90   if(olanl) then
       cc( 1, 1) =        -10.00000000d0
       cc( 2, 1) =        -57.94093980d0
       cc( 3, 1) =        -10.26232850d0
       a ( 1, 1) =        329.27205390d0
       a ( 2, 1) =         59.54875970d0
       a ( 3, 1) =         14.82659660d0
       n ( 1, 1) =   1
       n ( 2, 1) =   2
       n ( 3, 1) =   2
       nt( 1)    =   3
       cc( 1, 2) =          3.00000000d0
       cc( 2, 2) =         25.25345990d0
       cc( 3, 2) =        210.64723660d0
       cc( 4, 2) =        235.28195870d0
       cc( 5, 2) =       -152.53434610d0
       a ( 1, 2) =        100.02606530d0
       a ( 2, 2) =         67.84681240d0
       a ( 3, 2) =         34.67710060d0
       a ( 4, 2) =          7.67570050d0
       a ( 5, 2) =          7.20774270d0
       n ( 1, 2) =   0
       n ( 2, 2) =   1
       n ( 3, 2) =   2
       n ( 4, 2) =   2
       n ( 5, 2) =   2
       nt( 2)    =   5
       cc( 1, 3) =          5.00000000d0
       cc( 2, 3) =          5.63902660d0
       cc( 3, 3) =        164.83621980d0
       cc( 4, 3) =        106.01680180d0
       cc( 5, 3) =        -69.85854680d0
       a ( 1, 3) =         47.97233210d0
       a ( 2, 3) =         59.62248730d0
       a ( 3, 3) =         26.04326000d0
       a ( 4, 3) =          6.55826920d0
       a ( 5, 3) =          6.12682360d0
       n ( 1, 3) =   0
       n ( 2, 3) =   1
       n ( 3, 3) =   2
       n ( 4, 3) =   2
       n ( 5, 3) =   2
       nt( 3)    =   5
       nelcor    =  10
       lmx       =   2
      else
       lmx = 3
       nelcor = 18
       nt(1) = 5
       nt(2) = 4
       nt(3) = 5
       nt(4) = 2
       cc(1,1) = -18.0000000d0
       cc(2,1) = -103.3668518d0
       cc(3,1) = -26.6610566d0
       cc(4,1) = -9.6086551d0
       cc(5,1) = -1.1868727d0
       cc(1,2) = 3.0000000d0
       cc(2,2) = 16.9973360d0
       cc(3,2) = 24.3894723d0
       cc(4,2) = -14.0418751d0
       cc(1,3) = 5.0000000d0
       cc(2,3) = 17.3283897d0
       cc(3,3) = 25.9223919d0
       cc(4,3) = 32.7385948d0
       cc(5,3) = 5.2698153d0
       cc(1,4) = 0.3010003d0
       cc(2,4) = -0.3719276d0
c
       a(1,1) = 261.7636243d0
       a(2,1) = 49.1963290d0
       a(3,1) = 11.1135566d0
       a(4,1) = 2.8486023d0
       a(5,1) = 0.9433796d0
       a(1,2) = 11.5120130d0
       a(2,2) = 2.9984562d0
       a(3,2) = 0.8135381d0
       a(4,2) = 0.7315434d0
       a(1,3) = 18.6651904d0
       a(2,3) = 8.4039778d0
       a(3,3) = 8.8945784d0
       a(4,3) = 2.8776940d0
       a(5,3) = 0.6691169d0
       a(1,4) = 2.1397301d0
       a(2,4) = 0.4773511d0
       n(1,1) = 1
       n(2,1) = 2
       n(3,1) = 2
       n(4,1) = 2
       n(5,1) = 2
       n(1,2) = 0
       n(2,2) = 1
       n(3,2) = 2
       n(4,2) = 2
       n(1,3) = 0
       n(2,3) = 1
       n(3,3) = 2
       n(4,3) = 2
       n(5,3) = 2
       n(1,4) = 2
       n(2,4) = 2
      endif
      go to 220
c
c     manganese  hay
c
 100  if(olanl) then
       cc( 1, 1) =        -10.00000000d0
       cc( 2, 1) =        -60.29902870d0
       cc( 3, 1) =        -10.42178340d0
       a ( 1, 1) =        357.39144690d0
       a ( 2, 1) =         64.64773890d0
       a ( 3, 1) =         16.09608330d0
       n ( 1, 1) =   1
       n ( 2, 1) =   2
       n ( 3, 1) =   2
       nt( 1)    =   3
       cc( 1, 2) =          3.00000000d0
       cc( 2, 2) =         16.25918190d0
       cc( 3, 2) =        276.93739280d0
       cc( 4, 2) =        241.31743420d0
       cc( 5, 2) =       -146.46353290d0
       a ( 1, 2) =        107.41272150d0
       a ( 2, 2) =        111.49589730d0
       a ( 3, 2) =         46.55683460d0
       a ( 4, 2) =          8.36881350d0
       a ( 5, 2) =          7.72374890d0
       n ( 1, 2) =   0
       n ( 2, 2) =   1
       n ( 3, 2) =   2
       n ( 4, 2) =   2
       n ( 5, 2) =   2
       nt( 2)    =   5
       cc( 1, 3) =          5.00000000d0
       cc( 2, 3) =          5.75897560d0
       cc( 3, 3) =        285.29186540d0
       cc( 4, 3) =        143.42226470d0
       cc( 5, 3) =        -88.70318510d0
       a ( 1, 3) =         80.04151030d0
       a ( 2, 3) =        105.60436460d0
       a ( 3, 3) =         40.83004660d0
       a ( 4, 3) =          8.00984570d0
       a ( 5, 3) =          7.33909280d0
       n ( 1, 3) =   0
       n ( 2, 3) =   1
       n ( 3, 3) =   2
       n ( 4, 3) =   2
       n ( 5, 3) =   2
       nt( 3)    =   5
       nelcor    =  10
       lmx       =   2
      else
       lmx = 3
       nelcor = 18
       nt(1) = 5
       nt(2) = 4
       nt(3) = 5
       nt(4) = 2
       cc(1,1) = -18.0000000d0
       cc(2,1) = -107.7292226d0
       cc(3,1) = -27.6423345d0
       cc(4,1) = -9.9330850d0
       cc(5,1) = -1.1330159d0
       cc(1,2) = 3.0000000d0
       cc(2,2) = 18.6308433d0
       cc(3,2) = 29.9037082d0
       cc(4,2) = -18.5369130d0
       cc(1,3) = 5.0000000d0
       cc(2,3) = 17.3564739d0
       cc(3,3) = 57.1846489d0
       cc(4,3) = 47.4013059d0
       cc(5,3) = 5.8126748d0
       cc(1,4) = 0.2968466d0
       cc(2,4) = -0.3897128d0
c
       a(1,1) = 284.6267218d0
       a(2,1) = 53.5171467d0
       a(3,1) = 12.0450769d0
       a(4,1) = 3.1166539d0
       a(5,1) = 1.0160414d0
       a(1,2) = 12.7311355d0
       a(2,2) = 3.2793085d0
       a(3,2) = 0.9046062d0
       a(4,2) = 0.8275742d0
       a(1,3) = 28.8651113d0
       a(2,3) = 13.5526575d0
       a(3,3) = 14.1586483d0
       a(4,3) = 3.9923448d0
       a(5,3) = 0.7265806d0
       a(1,4) = 2.6595345d0
       a(2,4) = 0.5256942d0
       n(1,1) = 1
       n(2,1) = 2
       n(3,1) = 2
       n(4,1) = 2
       n(5,1) = 2
       n(1,2) = 0
       n(2,2) = 1
       n(3,2) = 2
       n(4,2) = 2
       n(1,3) = 0
       n(2,3) = 1
       n(3,3) = 2
       n(4,3) = 2
       n(5,3) = 2
       n(1,4) = 2
       n(2,4) = 2
      endif
      go to 220
c
c     iron  hay
c
 110  if(olanl) then
       cc( 1, 1) =        -10.00000000d0
       cc( 2, 1) =        -63.26675180d0
       cc( 3, 1) =        -10.96133380d0
       a ( 1, 1) =        392.61497870d0
       a ( 2, 1) =         71.17569790d0
       a ( 3, 1) =         17.73202810d0
       n ( 1, 1) =   1
       n ( 2, 1) =   2
       n ( 3, 1) =   2
       nt( 1)    =   3
       cc( 1, 2) =          3.00000000d0
       cc( 2, 2) =         18.17291370d0
       cc( 3, 2) =        339.12311640d0
       cc( 4, 2) =        317.10680120d0
       cc( 5, 2) =       -207.34216490d0
       a ( 1, 2) =        126.05718950d0
       a ( 2, 2) =        138.12642510d0
       a ( 3, 2) =         54.20988580d0
       a ( 4, 2) =          9.28379660d0
       a ( 5, 2) =          8.62890820d0
       n ( 1, 2) =   0
       n ( 2, 2) =   1
       n ( 3, 2) =   2
       n ( 4, 2) =   2
       n ( 5, 2) =   2
       nt( 2)    =   5
       cc( 1, 3) =          5.00000000d0
       cc( 2, 3) =          5.95359300d0
       cc( 3, 3) =        294.26655270d0
       cc( 4, 3) =        154.42446350d0
       cc( 5, 3) =        -95.31642490d0
       a ( 1, 3) =         83.17594900d0
       a ( 2, 3) =        106.05599380d0
       a ( 3, 3) =         42.82849370d0
       a ( 4, 3) =          8.77018050d0
       a ( 5, 3) =          8.03978180d0
       n ( 1, 3) =   0
       n ( 2, 3) =   1
       n ( 3, 3) =   2
       n ( 4, 3) =   2
       n ( 5, 3) =   2
       nt( 3)    =   5
       nelcor    =  10
       lmx       =   2
      else
       lmx = 3
       nelcor = 18
       nt(1) = 5
       nt(2) = 4
       nt(3) = 5
       nt(4) = 2
       cc(1,1) = -18.0000000d0
       cc(2,1) = -111.5658307d0
       cc(3,1) = -28.4011758d0
       cc(4,1) = -10.1669409d0
       cc(5,1) = -1.0631887d0
       cc(1,2) = 3.0000000d0
       cc(2,2) = 17.4313169d0
       cc(3,2) = 29.4117858d0
       cc(4,2) = -16.8563282d0
       cc(1,3) = 5.0000000d0
       cc(2,3) = 18.1533743d0
       cc(3,3) = 54.0296157d0
       cc(4,3) = 48.7957799d0
       cc(5,3) = 6.2143968d0
       cc(1,4) = 0.2738868d0
       cc(2,4) = -0.3947896d0
c
       a(1,1) = 305.9908267d0
       a(2,1) = 57.5401809d0
       a(3,1) = 12.8787208d0
       a(4,1) = 3.3671745d0
       a(5,1) = 1.0778632d0
       a(1,2) = 16.8558207d0
       a(2,2) = 4.4141716d0
       a(3,2) = 0.9390046d0
       a(4,2) = 0.8294794d0
       a(1,3) = 28.7494405d0
       a(2,3) = 13.5008313d0
       a(3,3) = 14.0936860d0
       a(4,3) = 4.1439477d0
       a(5,3) = 0.7948771d0
       a(1,4) = 3.4156612d0
       a(2,4) = 0.5623151d0
       n(1,1) = 1
       n(2,1) = 2
       n(3,1) = 2
       n(4,1) = 2
       n(5,1) = 2
       n(1,2) = 0
       n(2,2) = 1
       n(3,2) = 2
       n(4,2) = 2
       n(1,3) = 0
       n(2,3) = 1
       n(3,3) = 2
       n(4,3) = 2
       n(5,3) = 2
       n(1,4) = 2
       n(2,4) = 2
      endif
      go to 220
c
c     cobalt  hay
c
 120  if(olanl) then
       cc( 1, 1) =        -10.00000000d0
       cc( 2, 1) =        -66.72402760d0
       cc( 3, 1) =        -11.71922430d0
       a ( 1, 1) =        434.10761040d0
       a ( 2, 1) =         78.84027580d0
       a ( 3, 1) =         19.60377490d0
       n ( 1, 1) =   1
       n ( 2, 1) =   2
       n ( 3, 1) =   2
       nt( 1)    =   3
       cc( 1, 2) =          3.00000000d0
       cc( 2, 2) =         20.15531100d0
       cc( 3, 2) =        413.37622440d0
       cc( 4, 2) =        133.87277050d0
       a ( 1, 2) =        151.23075820d0
       a ( 2, 2) =        169.67900260d0
       a ( 3, 2) =         64.99094230d0
       a ( 4, 2) =         12.50762870d0
       n ( 1, 2) =   0
       n ( 2, 2) =   1
       n ( 3, 2) =   2
       n ( 4, 2) =   2
       nt( 2)    =   4
       cc( 1, 3) =          5.00000000d0
       cc( 2, 3) =          6.18923530d0
       cc( 3, 3) =        326.32159190d0
       cc( 4, 3) =         73.87889700d0
       a ( 1, 3) =         92.73199240d0
       a ( 2, 3) =        121.19101830d0
       a ( 3, 3) =         48.62659090d0
       a ( 4, 3) =         11.99669710d0
       n ( 1, 3) =   0
       n ( 2, 3) =   1
       n ( 3, 3) =   2
       n ( 4, 3) =   2
       nt( 3)    =   4
       nelcor    =  10
       lmx       =   2
      else
       lmx = 3
       nelcor = 18
       nt(1) = 5
       nt(2) = 4
       nt(3) = 5
       nt(4) = 2
       cc(1,1) = -18.0000000d0
       cc(2,1) = -114.9386822d0
       cc(3,1) = -28.9823031d0
       cc(4,1) = -10.3160176d0
       cc(5,1) = -0.9816066d0
       cc(1,2) = 3.0000000d0
       cc(2,2) = 18.3122881d0
       cc(3,2) = 27.0944428d0
       cc(4,2) = -13.5591147d0
       cc(1,3) = 5.0000000d0
       cc(2,3) = 18.9879649d0
       cc(3,3) = 54.5504938d0
       cc(4,3) = 51.6244450d0
       cc(5,3) = 6.8211010d0
       cc(1,4) = 0.2655202d0
       cc(2,4) = -0.4143761d0
c
       a(1,1) = 325.9446323d0
       a(2,1) = 61.2895796d0
       a(3,1) = 13.6219044d0
       a(4,1) = 3.6019089d0
       a(5,1) = 1.1293773d0
       a(1,2) = 17.6813526d0
       a(2,2) = 4.6329374d0
       a(3,2) = 1.0443458d0
       a(4,2) = 0.8975469d0
       a(1,3) = 29.7562169d0
       a(2,3) = 13.9792163d0
       a(3,3) = 14.5834185d0
       a(4,3) = 4.4052039d0
       a(5,3) = 0.8947474d0
       a(1,4) = 4.0233194d0
       a(2,4) = 0.6152352d0
       n(1,1) = 1
       n(2,1) = 2
       n(3,1) = 2
       n(4,1) = 2
       n(5,1) = 2
       n(1,2) = 0
       n(2,2) = 1
       n(3,2) = 2
       n(4,2) = 2
       n(1,3) = 0
       n(2,3) = 1
       n(3,3) = 2
       n(4,3) = 2
       n(5,3) = 2
       n(1,4) = 2
       n(2,4) = 2
      endif
      go to 220
c
c     nickel (3d)9 (4s)1  hay ni01
c
 130  if(olanl) then
       cc( 1, 1) =        -10.00000000d0
       cc( 2, 1) =        -69.40848050d0
       cc( 3, 1) =        -12.09510200d0
       a ( 1, 1) =        469.93243310d0
       a ( 2, 1) =         85.42364110d0
       a ( 3, 1) =         21.26749840d0
       n ( 1, 1) =   1
       n ( 2, 1) =   2
       n ( 3, 1) =   2
       nt( 1)    =   3
       cc( 1, 2) =          3.00000000d0
       cc( 2, 2) =         22.02536180d0
       cc( 3, 2) =        443.01810880d0
       cc( 4, 2) =        145.56964110d0
       a ( 1, 2) =        162.16860970d0
       a ( 2, 2) =        176.53332320d0
       a ( 3, 2) =         68.95620100d0
       a ( 4, 2) =         13.57928380d0
       n ( 1, 2) =   0
       n ( 2, 2) =   1
       n ( 3, 2) =   2
       n ( 4, 2) =   2
       nt( 2)    =   4
       cc( 1, 3) =          5.00000000d0
       cc( 2, 3) =          4.98828240d0
       cc( 3, 3) =        256.69458530d0
       cc( 4, 3) =         78.47544500d0
       a ( 1, 3) =         69.01815060d0
       a ( 2, 3) =        275.59555960d0
       a ( 3, 3) =         47.13154530d0
       a ( 4, 3) =         12.98740750d0
       n ( 1, 3) =   0
       n ( 2, 3) =   1
       n ( 3, 3) =   2
       n ( 4, 3) =   2
       nt( 3)    =   4
       nelcor    =  10
       lmx       =   2
      else
       lmx = 3
       nelcor = 18
       nt(1) = 5
       nt(2) = 4
       nt(3) = 5
       nt(4) = 2
       cc(1,1) = -18.0000000d0
       cc(2,1) = -117.9593651d0
       cc(3,1) = -29.4397021d0
       cc(4,1) = -10.3862631d0
       cc(5,1) = -0.8924874d0
       cc(1,2) = 3.0000000d0
       cc(2,2) = 19.2448974d0
       cc(3,2) = 23.9306008d0
       cc(4,2) = -9.3541445d0
       cc(1,3) = 5.0000000d0
       cc(2,3) = 19.8115548d0
       cc(3,3) = 54.3385586d0
       cc(4,3) = 54.0878232d0
       cc(5,3) = 7.3102740d0
       cc(1,4) = 0.2629228d0
       cc(2,4) = -0.4386154d0
c
       a(1,1) = 344.8409988d0
       a(2,1) = 64.8228067d0
       a(3,1) = 14.2847653d0
       a(4,1) = 3.8210097d0
       a(5,1) = 1.1697587d0
       a(1,2) = 18.6423811d0
       a(2,2) = 4.8916063d0
       a(3,2) = 1.1660615d0
       a(4,2) = 0.9523942d0
       a(1,3) = 30.6007023d0
       a(2,3) = 14.3008055d0
       a(3,3) = 15.0330432d0
       a(4,3) = 4.6460085d0
       a(5,3) = 0.9810619d0
       a(1,4) = 4.5600829d0
       a(2,4) = 0.6764731d0
       n(1,1) = 1
       n(2,1) = 2
       n(3,1) = 2
       n(4,1) = 2
       n(5,1) = 2
       n(1,2) = 0
       n(2,2) = 1
       n(3,2) = 2
       n(4,2) = 2
       n(1,3) = 0
       n(2,3) = 1
       n(3,3) = 2
       n(4,3) = 2
       n(5,3) = 2
       n(1,4) = 2
       n(2,4) = 2
      endif
      go to 220
c
c     copper  hay
c
 140  if(olanl) then
       cc( 1, 1) =        -10.00000000d0
       cc( 2, 1) =        -72.55482820d0
       cc( 3, 1) =        -12.74502310d0
       a ( 1, 1) =        511.99517630d0
       a ( 2, 1) =         93.28010740d0
       a ( 3, 1) =         23.22066690d0
       n ( 1, 1) =   1
       n ( 2, 1) =   2
       n ( 3, 1) =   2
       nt( 1)    =   3
       cc( 1, 2) =          3.00000000d0
       cc( 2, 2) =         23.83518250d0
       cc( 3, 2) =        473.89304880d0
       cc( 4, 2) =        157.63458230d0
       a ( 1, 2) =        173.11808540d0
       a ( 2, 2) =        185.24198860d0
       a ( 3, 2) =         73.15178470d0
       a ( 4, 2) =         14.68841570d0
       n ( 1, 2) =   0
       n ( 2, 2) =   1
       n ( 3, 2) =   2
       n ( 4, 2) =   2
       nt( 2)    =   4
       cc( 1, 3) =          5.00000000d0
       cc( 2, 3) =          6.49909360d0
       cc( 3, 3) =        351.46053950d0
       cc( 4, 3) =         85.50160360d0
       a ( 1, 3) =        100.71913690d0
       a ( 2, 3) =        130.83456650d0
       a ( 3, 3) =         53.86837200d0
       a ( 4, 3) =         14.09894690d0
       n ( 1, 3) =   0
       n ( 2, 3) =   1
       n ( 3, 3) =   2
       n ( 4, 3) =   2
       nt( 3)    =   4
       nelcor    =  10
       lmx       =   2
      else
       lmx = 3
       nelcor = 18
       nt(1) = 5
       nt(2) = 4
       nt(3) = 5
       nt(4) = 2
       cc(1,1) = -18.0000000d0
       cc(2,1) = -119.9259397d0
       cc(3,1) = -29.5532867d0
       cc(4,1) = -10.2892433d0
       cc(5,1) = -0.7836363d0
       cc(1,2) = 3.0000000d0
       cc(2,2) = 20.1579275d0
       cc(3,2) = 34.5001906d0
       cc(4,2) = -18.9812003d0
       cc(1,3) = 5.0000000d0
       cc(2,3) = 20.6085328d0
       cc(3,3) = 56.0016888d0
       cc(4,3) = 57.2170107d0
       cc(5,3) = 7.7177878d0
       cc(1,4) = 0.2598616d0
       cc(2,4) = -0.4621680d0
c
       a(1,1) = 359.2137111d0
       a(2,1) = 67.5347369d0
       a(3,1) = 14.7222923d0
       a(4,1) = 3.9975558d0
       a(5,1) = 1.1889410d0
       a(1,2) = 19.6202650d0
       a(2,2) = 5.1604389d0
       a(3,2) = 1.2306099d0
       a(4,2) = 1.0850105d0
       a(1,3) = 31.9385762d0
       a(2,3) = 14.9202125d0
       a(3,3) = 15.6835232d0
       a(4,3) = 4.9311614d0
       a(5,3) = 1.0622167d0
       a(1,4) = 5.1159991d0
       a(2,4) = 0.7396784d0
       n(1,1) = 1
       n(2,1) = 2
       n(3,1) = 2
       n(4,1) = 2
       n(5,1) = 2
       n(1,2) = 0
       n(2,2) = 1
       n(3,2) = 2
       n(4,2) = 2
       n(1,3) = 0
       n(2,3) = 1
       n(3,3) = 2
       n(4,3) = 2
       n(5,3) = 2
       n(1,4) = 2
       n(2,4) = 2
      endif
      go to 220
c
c     zinc  hay (Zn_LANL2DZ_ECP)
c
 150  continue
       lmx = 3
       nelcor = 18
       nt(1) = 5
       nt(2) = 5
       nt(3) = 5
       nt(4) = 3
       cc(1,1) = -18.0000000d0
       cc(2,1) = -124.3527403d0
       cc(3,1) = -30.6601822d0
       cc(4,1) = -10.6358989d0
       cc(5,1) = -0.7683623d0
       cc(1,2) = 3.0000000d0
       cc(2,2) = 22.5234225d0
       cc(3,2) = 48.4465942d0
       cc(4,2) = -44.5560119d0
       cc(5,2) = 12.9983958d0
       cc(1,3) = 5.0000000d0
       cc(2,3) = 20.7435589d0
       cc(3,3) = 90.3027158d0
       cc(4,3) = 74.6610316d0
       cc(5,3) = 9.8894424d0
       cc(1,4) = -4.8490359d0
       cc(2,4) = 3.6913379d0
       cc(3,4) = -0.5037319d0
c
       a(1,1) = 386.7379660d0
       a(2,1) = 72.8587359d0
       a(3,1) = 15.9066170d0
       a(4,1) = 4.3502340d0
       a(5,1) = 1.2842199d0
       a(1,2) = 19.0867858d0
       a(2,2) = 5.0231080d0
       a(3,2) = 1.2701744d0
       a(4,2) = 1.0671287d0
       a(5,2) = 0.9264190d0
       a(1,3) = 43.4927750d0
       a(2,3) = 20.8692669d0
       a(3,3) = 21.7118378d0
       a(4,3) = 6.3616915d0
       a(5,3) = 1.2291195d0
       a(1,4) = 13.5851800d0
       a(2,4) = 9.8373050d0
       a(3,4) = 0.8373113d0
       n(1,1) = 1
       n(2,1) = 2
       n(3,1) = 2
       n(4,1) = 2
       n(5,1) = 2
       n(1,2) = 0
       n(2,2) = 1
       n(3,2) = 2
       n(4,2) = 2
       n(5,2) = 2
       n(1,3) = 0
       n(2,3) = 1
       n(3,3) = 2
       n(4,3) = 2
       n(5,3) = 2
       n(1,4) = 2
       n(2,4) = 2
       n(3,4) = 2
      go to 220
c     ga  hay  (Ga_LANL2DZ_ECP)
 160  lmx = 3
      nelcor = 28
      nt(1) = 5
      nt(2) = 5
      nt(3) = 5
      nt(4) = 5
      cc(1,1) = -28.0000000d0
      cc(2,1) = -161.6248634d0
      cc(3,1) = -49.2659751d0
      cc(4,1) = -16.9305942d0
      cc(5,1) = -1.5149893d0
      cc(1,2) = 3.0000000d0
      cc(2,2) = 93.8597449d0
      cc(3,2) = 308.6495485d0
      cc(4,2) = 94.7290607d0
      cc(5,2) = 27.0931066d0
      cc(1,3) = 5.0000000d0
      cc(2,3) = 22.8495315d0
      cc(3,3) = 43.3974907d0
      cc(4,3) = 48.6203923d0
      cc(5,3) = 14.7365483d0
      cc(1,4) = 3.0000000d0
      cc(2,4) = 24.2994834d0
      cc(3,4) = 35.4415813d0
      cc(4,4) = 9.1055575d0
      cc(5,4) = 0.9601121d0
c
      a(1,1) = 260.9051184d0
      a(2,1) = 50.4244564d0
      a(3,1) = 10.7195520d0
      a(4,1) = 3.1563995d0
      a(5,1) = 0.9502385d0
      a(1,2) = 339.0575995d0
      a(2,2) = 95.9625227d0
      a(3,2) = 30.0989360d0
      a(4,2) = 7.8907039d0
      a(5,2) = 2.1154234d0
      a(1,3) = 30.1202708d0
      a(2,3) = 14.6143825d0
      a(3,3) = 12.4639433d0
      a(4,3) = 5.1039219d0
      a(5,3) = 1.6566033d0
      a(1,4) = 32.2704306d0
      a(2,4) = 10.9444004d0
      a(3,4) = 4.5111791d0
      a(4,4) = 1.5561710d0
      a(5,4) = 0.5550926d0
      n(1,1) = 1
      n(2,1) = 2
      n(3,1) = 2
      n(4,1) = 2
      n(5,1) = 2
      n(1,2) = 0
      n(2,2) = 1
      n(3,2) = 2
      n(4,2) = 2
      n(5,2) = 2
      n(1,3) = 0
      n(2,3) = 1
      n(3,3) = 2
      n(4,3) = 2
      n(5,3) = 2
      n(1,4) = 0
      n(2,4) = 1
      n(3,4) = 2
      n(4,4) = 2
      n(5,4) = 2
      go to 220
c     ge hay  (Ge_LANL2DZ_ECP)
 170  lmx = 3
      nelcor = 28
      nt(1) = 5
      nt(2) = 5
      nt(3) = 5
      nt(4) = 5
      cc(1,1) = -28.0000000d0
      cc(2,1) = -180.9891676d0
      cc(3,1) = -55.0043909d0
      cc(4,1) = -19.7906526d0
      cc(5,1) = -1.8533572d0
      cc(1,2) = 3.0000000d0
      cc(2,2) = 65.2262558d0
      cc(3,2) = 225.2354522d0
      cc(4,2) = 94.0125472d0
      cc(5,2) = 29.9415005d0
      cc(1,3) = 5.0000000d0
      cc(2,3) = 23.4778157d0
      cc(3,3) = 45.0980414d0
      cc(4,3) = 56.3326957d0
      cc(5,3) = 16.6058640d0
      cc(1,4) = 3.0000000d0
      cc(2,4) = 23.7371518d0
      cc(3,4) = 56.4792249d0
      cc(4,4) = 25.8901835d0
      cc(5,4) = 3.0229836d0
c
      a(1,1) = 318.2167583d0
      a(2,1) = 61.5370967d0
      a(3,1) = 13.2986899d0
      a(4,1) = 3.8985215d0
      a(5,1) = 1.2137666d0
      a(1,2) = 205.1886932d0
      a(2,2) = 68.9790278d0
      a(3,2) = 27.9194879d0
      a(4,2) = 8.5481650d0
      a(5,2) = 2.3173734d0
      a(1,3) = 33.2488002d0
      a(2,3) = 15.7777247d0
      a(3,3) = 14.9260722d0
      a(4,3) = 5.8416394d0
      a(5,3) = 1.8349575d0
      a(1,4) = 42.0206343d0
      a(2,4) = 19.2096363d0
      a(3,4) = 9.4133917d0
      a(4,4) = 3.3282907d0
      a(5,4) = 0.8522331d0
      n(1,1) = 1
      n(2,1) = 2
      n(3,1) = 2
      n(4,1) = 2
      n(5,1) = 2
      n(1,2) = 0
      n(2,2) = 1
      n(3,2) = 2
      n(4,2) = 2
      n(5,2) = 2
      n(1,3) = 0
      n(2,3) = 1
      n(3,3) = 2
      n(4,3) = 2
      n(5,3) = 2
      n(1,4) = 0
      n(2,4) = 1
      n(3,4) = 2
      n(4,4) = 2
      n(5,4) = 2
      go to 220
c     as  hay (As_LANL2DZ_ECP)
 180  lmx = 3
      nelcor = 28
      nt(1) = 5
      nt(2) = 5
      nt(3) = 5
      nt(4) = 5
      cc(1,1) = -28.0000000d0
      cc(2,1) = -198.3687357d0
      cc(3,1) = -60.4563798d0
      cc(4,1) = -22.3188274d0
      cc(5,1) = -2.1638290d0
      cc(1,2) = 3.0000000d0
      cc(2,2) = 31.8486777d0
      cc(3,2) = 200.9958215d0
      cc(4,2) = 102.6387604d0
      cc(5,2) = 34.1845359d0
      cc(1,3) = 5.0000000d0
      cc(2,3) = 24.2505747d0
      cc(3,3) = 48.9174963d0
      cc(4,3) = 66.5332259d0
      cc(5,3) = 19.6577666d0
      cc(1,4) = 3.0000000d0
      cc(2,4) = 23.3460712d0
      cc(3,4) = 80.7854919d0
      cc(4,4) = 40.2681347d0
      cc(5,4) = 4.8332790d0
      a(1,1) = 375.4748803d0
      a(2,1) = 72.6554769d0
      a(3,1) = 15.9162677d0
      a(4,1) = 4.6420251d0
      a(5,1) = 1.4789182d0
      a(1,2) = 105.8853870d0
      a(2,2) = 62.5951846d0
      a(3,2) = 34.5017431d0
      a(4,2) = 10.3758594d0
      a(5,2) = 2.5379347d0
      a(1,3) = 37.6471440d0
      a(2,3) = 17.2034665d0
      a(3,3) = 18.9075709d0
      a(4,3) = 6.8153452d0
      a(5,3) = 2.0611853d0
      a(1,4) = 52.0731335d0
      a(2,4) = 28.0037824d0
      a(3,4) = 14.7051208d0
      a(4,4) = 4.7973194d0
      a(5,4) = 1.1062763d0
      n(1,1) = 1
      n(2,1) = 2
      n(3,1) = 2
      n(4,1) = 2
      n(5,1) = 2
      n(1,2) = 0
      n(2,2) = 1
      n(3,2) = 2
      n(4,2) = 2
      n(5,2) = 2
      n(1,3) = 0
      n(2,3) = 1
      n(3,3) = 2
      n(4,3) = 2
      n(5,3) = 2
      n(1,4) = 0
      n(2,4) = 1
      n(3,4) = 2
      n(4,4) = 2
      n(5,4) = 2
      go to 220
c     se hay  (Se_LANL2DZ_ECP)
 190  lmx = 3
      nelcor = 28
      nt(1) = 5
      nt(2) = 5
      nt(3) = 5
      nt(4) = 5
      cc(1,1) = -28.0000000d0
      cc(2,1) = -214.3841762d0
      cc(3,1) = -65.6918782d0
      cc(4,1) = -24.6153932d0
      cc(5,1) = -2.4481497d0
      cc(1,2) = 3.0000000d0
      cc(2,2) = 62.0295390d0
      cc(3,2) = 258.8555984d0
      cc(4,2) = 118.7800153d0
      cc(5,2) = 38.2355279d0
      cc(1,3) = 5.0000000d0
      cc(2,3) = 24.7973458d0
      cc(3,3) = 63.7575640d0
      cc(4,3) = 79.0512831d0
      cc(5,3) = 22.9520183d0
      cc(1,4) = 3.0000000d0
      cc(2,4) = 22.4705907d0
      cc(3,4) = 140.5492887d0
      cc(4,4) = 63.5781835d0
      cc(5,4) = 7.0753614d0
c
      a(1,1) = 433.1931336d0
      a(2,1) = 83.8952157d0
      a(3,1) = 18.5839139d0
      a(4,1) = 5.3955286d0
      a(5,1) = 1.7474326d0
      a(1,2) = 202.8986193d0
      a(2,2) = 78.3820487d0
      a(3,2) = 35.0753037d0
      a(4,2) = 10.8769543d0
      a(5,2) = 2.8005941d0
      a(1,3) = 44.3011875d0
      a(2,3) = 20.3874206d0
      a(3,3) = 23.1889948d0
      a(4,3) = 7.9777664d0
      a(5,3) = 2.2988146d0
      a(1,4) = 73.3628263d0
      a(2,4) = 48.3835618d0
      a(3,4) = 25.6297211d0
      a(4,4) = 7.1705822d0
      a(5,4) = 1.3639538d0
      n(1,1) = 1
      n(2,1) = 2
      n(3,1) = 2
      n(4,1) = 2
      n(5,1) = 2
      n(1,2) = 0
      n(2,2) = 1
      n(3,2) = 2
      n(4,2) = 2
      n(5,2) = 2
      n(1,3) = 0
      n(2,3) = 1
      n(3,3) = 2
      n(4,3) = 2
      n(5,3) = 2
      n(1,4) = 0
      n(2,4) = 1
      n(3,4) = 2
      n(4,4) = 2
      n(5,4) = 2
      go to 220
c     br  hay (Br_LANL2DZ_ECP)
 200  lmx = 3
      nelcor = 28
      nt(1) = 4
      nt(2) = 4
      nt(3) = 5
      nt(4) = 5
      cc(1,1) = -28.0000000d0
      cc(2,1) = -134.9268852d0
      cc(3,1) = -41.9271913d0
      cc(4,1) = -5.9336420d0
      cc(1,2) = 3.0000000d0
      cc(2,2) = 27.3430642d0
      cc(3,2) = 118.8028847d0
      cc(4,2) = 43.4354876d0
      cc(1,3) = 5.0000000d0
      cc(2,3) = 25.0504252d0
      cc(3,3) = 92.6157463d0
      cc(4,3) = 95.8249016d0
      cc(5,3) = 26.2684983d0
      cc(1,4) = 3.0000000d0
      cc(2,4) = 22.5533557d0
      cc(3,4) = 178.1241988d0
      cc(4,4) = 76.9924162d0
      cc(5,4) = 9.4818270d0
c
      a(1,1) = 213.6143969d0
      a(2,1) = 41.0585380d0
      a(3,1) = 8.7086530d0
      a(4,1) = 2.6074661d0
      a(1,2) = 54.1980682d0
      a(2,2) = 32.9053558d0
      a(3,2) = 13.6744890d0
      a(4,2) = 3.0341152d0
      a(1,3) = 54.2563340d0
      a(2,3) = 26.0095593d0
      a(3,3) = 28.2012995d0
      a(4,3) = 9.4341061d0
      a(5,3) = 2.5321764d0
      a(1,4) = 87.6328721d0
      a(2,4) = 61.7373377d0
      a(3,4) = 32.4385104d0
      a(4,4) = 8.7537199d0
      a(5,4) = 1.6633189d0
      n(1,1) = 1
      n(2,1) = 2
      n(3,1) = 2
      n(4,1) = 2
      n(1,2) = 0
      n(2,2) = 1
      n(3,2) = 2
      n(4,2) = 2
      n(1,3) = 0
      n(2,3) = 1
      n(3,3) = 2
      n(4,3) = 2
      n(5,3) = 2
      n(1,4) = 0
      n(2,4) = 1
      n(3,4) = 2
      n(4,4) = 2
      n(5,4) = 2
      go to 220
c     kr hay  (Kr_LANL2DZ_ECP)
 210  lmx = 3
      nelcor = 28
      nt(1) = 4
      nt(2) = 4
      nt(3) = 5
      nt(4) = 5
      cc(1,1) = -28.0000000d0
      cc(2,1) = -143.4259955d0
      cc(3,1) = -44.3723838d0
      cc(4,1) = -6.2911252d0
      cc(1,2) = 3.0000000d0
      cc(2,2) = 37.1893647d0
      cc(3,2) = 132.2973958d0
      cc(4,2) = 47.9346985d0
      cc(1,3) = 5.0000000d0
      cc(2,3) = 26.1564008d0
      cc(3,3) = 96.3176501d0
      cc(4,3) = 103.6632794d0
      cc(5,3) = 31.0779315d0
      cc(1,4) = 3.0000000d0
      cc(2,4) = 24.4356537d0
      cc(3,4) = 85.7229402d0
      cc(4,4) = 91.6970578d0
      cc(5,4) = 11.8300709d0
c
      a(1,1) = 237.5470853d0
      a(2,1) = 45.7406754d0
      a(3,1) = 9.7073237d0
      a(4,1) = 2.9320738d0
      a(1,2) = 70.7606918d0
      a(2,2) = 33.6809745d0
      a(3,2) = 13.8075281d0
      a(4,2) = 3.3188254d0
      a(1,3) = 57.7676804d0
      a(2,3) = 27.6779453d0
      a(3,3) = 30.4757449d0
      a(4,3) = 10.3642639d0
      a(5,3) = 2.8470827d0
      a(1,4) = 85.3902277d0
      a(2,4) = 36.0781993d0
      a(3,4) = 38.0387668d0
      a(4,4) = 10.6069882d0
      a(5,4) = 1.9390148d0
      n(1,1) = 1
      n(2,1) = 2
      n(3,1) = 2
      n(4,1) = 2
      n(1,2) = 0
      n(2,2) = 1
      n(3,2) = 2
      n(4,2) = 2
      n(1,3) = 0
      n(2,3) = 1
      n(3,3) = 2
      n(4,3) = 2
      n(5,3) = 2
      n(1,4) = 0
      n(2,4) = 1
      n(3,4) = 2
      n(4,4) = 2
      n(5,4) = 2
c
c     load ecp data
c
 220  return
 6010 format (1x,'*** data not found for ecp type ',a4)
      end
**==ecplb2.f
      subroutine ecplb2 (zinp,lmx,nelcor,nt,n,cc,a,olanl)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character * 8 (z)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      dimension nt(6), n(10,6), cc(10,6), a(10,6)
c
      dimension ztit(18)
      data maxtyp/18/
      data ztit/'rb','sr',
     *          'y ','zr','nb','mo','tc','ru',
     *          'rh','pd','ag','cd','in','sn','sb','te','i ','xe'/
c
c     this routine contains a library of ecp data taken from
c     various literature sources.
c     the ecps for
c     rb,......xe are taken from hay,wadt, jcp 82 (1985) 270,284
      do 20 ityp = 1 , maxtyp
         if (zinp.eq.ztit(ityp)) go to 30
 20   continue
      write (iwr,6010) zinp
      call caserr2('no library available for requested ecp')
 30   go to (40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,
     +       200,210) , ityp
c     rubidium  hay
 40   lmx = 3
      if(olanl) then
       cc( 1, 1) =         -0.03643620d0
       cc( 2, 1) =        -20.80330520d0
       cc( 3, 1) =       -115.46619100d0
       cc( 4, 1) =        -42.50612930d0
       cc( 5, 1) =         -5.52478890d0
       a ( 1, 1) =        879.59836640d0
       a ( 2, 1) =        142.01886950d0
       a ( 3, 1) =         40.93468820d0
       a ( 4, 1) =          9.90848020d0
       a ( 5, 1) =          3.15560410d0
       n ( 1, 1) =   0
       n ( 2, 1) =   1
       n ( 3, 1) =   2
       n ( 4, 1) =   2
       n ( 5, 1) =   2
       nt( 1)    =   5
       cc( 1, 2) =          2.97479100d0
       cc( 2, 2) =         36.11704350d0
       cc( 3, 2) =        152.27410680d0
       cc( 4, 2) =         53.71519060d0
       a ( 1, 2) =         70.98072150d0
       a ( 2, 2) =         38.69498070d0
       a ( 3, 2) =         15.90820350d0
       a ( 4, 2) =          3.71526590d0
       n ( 1, 2) =   0
       n ( 2, 2) =   1
       n ( 3, 2) =   2
       n ( 4, 2) =   2
       nt( 2)    =   4
       cc( 1, 3) =          4.95510260d0
       cc( 2, 3) =         24.27262360d0
       cc( 3, 3) =        215.23414110d0
       cc( 4, 3) =        151.66876660d0
       cc( 5, 3) =         33.87502800d0
       a ( 1, 3) =         93.30586840d0
       a ( 2, 3) =         47.12213700d0
       a ( 3, 3) =         48.80506400d0
       a ( 4, 3) =         14.37849490d0
       a ( 5, 3) =          3.05478070d0
       n ( 1, 3) =   0
       n ( 2, 3) =   1
       n ( 3, 3) =   2
       n ( 4, 3) =   2
       n ( 5, 3) =   2
       nt( 3)    =   5
       cc( 1, 4) =          3.04652240d0
       cc( 2, 4) =         29.42436240d0
       cc( 3, 4) =         33.44318080d0
       cc( 4, 4) =         11.10294600d0
       a ( 1, 4) =         33.04774300d0
       a ( 2, 4) =          9.96018800d0
       a ( 3, 4) =          3.74035770d0
       a ( 4, 4) =          0.81242080d0
       n ( 1, 4) =   0
       n ( 2, 4) =   1
       n ( 3, 4) =   2
       n ( 4, 4) =   2
       nt( 4)    =   4
       nelcor    =  28
      else
       nelcor = 36
       nt(1) = 6
       nt(2) = 6
       nt(3) = 6
       nt(4) = 5
       cc(1,1) = -0.0354362d0
       cc(2,1) = -23.3062353d0
       cc(3,1) = -41.9708274d0
       cc(4,1) = -8.8054882d0
       cc(5,1) = -3.2762829d0
       cc(6,1) = -0.1823006d0
c
       cc(1,2) = 3.0516375d0
       cc(2,2) = 34.4023818d0
       cc(3,2) = 99.4966874d0
       cc(4,2) = 52.1437811d0
       cc(5,2) = 17.0760321d0
       cc(6,2) = 5.5402103d0
c
       cc(1,3) = 2.0753970d0
       cc(2,3) = 38.2771753d0
       cc(3,3) = 98.7705797d0
       cc(4,3) = 34.3833390d0
       cc(5,3) = 11.6403739d0
       cc(6,3) = 2.4729309d0
c
       cc(1,4) = 3.0303228d0
       cc(2,4) = 28.8442524d0
       cc(3,4) = 104.8105282d0
       cc(4,4) = 43.0034879d0
       cc(5,4) = 10.8173215d0
c
       a(1,1) = 145.0350895d0
       a(2,1) = 21.2749685d0
       a(3,1) = 6.6094585d0
       a(4,1) = 1.8705308d0
       a(5,1) = 0.6177734d0
       a(6,1) = 0.2487102d0
       a(1,2) = 66.0749102d0
       a(2,2) = 29.8443241d0
       a(3,2) = 14.8474854d0
       a(4,2) = 5.6648238d0
       a(5,2) = 1.6674589d0
       a(6,2) = 0.4485012d0
       a(1,3) = 81.0852315d0
       a(2,3) = 29.1588472d0
       a(3,3) = 11.2440746d0
       a(4,3) = 3.7325713d0
       a(5,3) = 1.0373688d0
       a(6,3) = 0.2728727d0
       a(1,4) = 56.8036783d0
       a(2,4) = 29.9914780d0
       a(3,4) = 14.2870365d0
       a(4,4) = 4.4996671d0
       a(5,4) = 0.9210966d0
       n(1,1) = 0
       n(2,1) = 1
       n(3,1) = 2
       n(4,1) = 2
       n(5,1) = 2
       n(6,1) = 2
       n(1,2) = 0
       n(2,2) = 1
       n(3,2) = 2
       n(4,2) = 2
       n(5,2) = 2
       n(6,2) = 2
       n(1,3) = 0
       n(2,3) = 1
       n(3,3) = 2
       n(4,3) = 2
       n(5,3) = 2
       n(6,3) = 2
       n(1,4) = 0
       n(2,4) = 1
       n(3,4) = 2
       n(4,4) = 2
       n(5,4) = 2
      endif
      go to 220
c
c     strontium  hay
c
 50   lmx = 3
      if(olanl) then
       cc( 1, 1) =         -0.03843230d0
       cc( 2, 1) =        -20.61742710d0
       cc( 3, 1) =       -101.17377440d0
       cc( 4, 1) =        -38.77436030d0
       cc( 5, 1) =         -4.64792430d0
       a ( 1, 1) =        782.38046310d0
       a ( 2, 1) =        124.65423380d0
       a ( 3, 1) =         36.98749660d0
       a ( 4, 1) =          9.88288190d0
       a ( 5, 1) =          3.28295880d0
       n ( 1, 1) =   0
       n ( 2, 1) =   1
       n ( 3, 1) =   2
       n ( 4, 1) =   2
       n ( 5, 1) =   2
       nt( 1)    =   5
       cc( 1, 2) =          2.99890220d0
       cc( 2, 2) =         25.65526690d0
       cc( 3, 2) =        183.18185330d0
       cc( 4, 2) =         58.43847390d0
       a ( 1, 2) =         59.32406310d0
       a ( 2, 2) =         55.20384720d0
       a ( 3, 2) =         20.46920920d0
       a ( 4, 2) =          3.95881410d0
       n ( 1, 2) =   0
       n ( 2, 2) =   1
       n ( 3, 2) =   2
       n ( 4, 2) =   2
       nt( 2)    =   4
       cc( 1, 3) =          4.95511890d0
       cc( 2, 3) =         25.44723670d0
       cc( 3, 3) =        203.80027800d0
       cc( 4, 3) =        155.05187400d0
       cc( 5, 3) =         39.36051920d0
       a ( 1, 3) =         92.12019910d0
       a ( 2, 3) =         46.81325590d0
       a ( 3, 3) =         48.65664320d0
       a ( 4, 3) =         14.95032380d0
       a ( 5, 3) =          3.42687850d0
       n ( 1, 3) =   0
       n ( 2, 3) =   1
       n ( 3, 3) =   2
       n ( 4, 3) =   2
       n ( 5, 3) =   2
       nt( 3)    =   5
       cc( 1, 4) =          3.00564510d0
       cc( 2, 4) =         26.70641190d0
       cc( 3, 4) =         74.57569010d0
       cc( 4, 4) =         63.17421210d0
       cc( 5, 4) =         20.29611620d0
       a ( 1, 4) =         65.82913010d0
       a ( 2, 4) =         32.72826210d0
       a ( 3, 4) =         21.11460300d0
       a ( 4, 4) =          9.10712920d0
       a ( 5, 4) =          2.81107540d0
       n ( 1, 4) =   0
       n ( 2, 4) =   1
       n ( 3, 4) =   2
       n ( 4, 4) =   2
       n ( 5, 4) =   2
       nt( 4)    =   5
       nelcor    =  28
      else
       nelcor = 36
       nt(1) = 6
       nt(2) = 6
       nt(3) = 6
       nt(4) = 6
       cc(1,1) = -0.0384323d0
       cc(2,1) = -22.6194918d0
       cc(3,1) = -36.4319079d0
       cc(4,1) = -9.0718265d0
       cc(5,1) = -2.6098384d0
       cc(6,1) = -0.0559578d0
       cc(1,2) = 3.0517388d0
       cc(2,2) = 40.1725046d0
       cc(3,2) = 131.8949569d0
       cc(4,2) = 72.7610900d0
       cc(5,2) = 24.2643629d0
       cc(6,2) = 7.8096640d0
       cc(1,3) = 5.0103024d0
       cc(2,3) = 29.7225305d0
       cc(3,3) = 67.7842103d0
       cc(4,3) = 37.9494334d0
       cc(5,3) = 14.3649831d0
       cc(6,3) = 3.6436641d0
       cc(1,4) = 2.9049721d0
       cc(2,4) = 33.9464015d0
       cc(3,4) = 83.0430063d0
       cc(4,4) = 46.2557229d0
       cc(5,4) = -27.7408480d0
       cc(6,4) = -0.0521022d0
c
       a(1,1) = 111.4938784d0
       a(2,1) = 19.3872731d0
       a(3,1) = 6.5321474d0
       a(4,1) = 1.8697150d0
       a(5,1) = 0.6304149d0
       a(6,1) = 0.2498766d0
       a(1,2) = 85.0969543d0
       a(2,2) = 38.8628138d0
       a(3,2) = 19.4115538d0
       a(4,2) = 7.2565109d0
       a(5,2) = 2.1660364d0
       a(6,2) = 0.5503453d0
       a(1,3) = 32.1182674d0
       a(2,3) = 15.5675759d0
       a(3,3) = 9.8816331d0
       a(4,3) = 4.0615003d0
       a(5,3) = 1.2652503d0
       a(6,3) = 0.3450867d0
       a(1,4) = 70.6203270d0
       a(2,4) = 24.1498704d0
       a(3,4) = 9.1278312d0
       a(4,4) = 1.5891572d0
       a(5,4) = 1.2856014d0
       a(6,4) = 0.2348379d0
       n(1,1) = 0
       n(2,1) = 1
       n(3,1) = 2
       n(4,1) = 2
       n(5,1) = 2
       n(6,1) = 2
       n(1,2) = 0
       n(2,2) = 1
       n(3,2) = 2
       n(4,2) = 2
       n(5,2) = 2
       n(6,2) = 2
       n(1,3) = 0
       n(2,3) = 1
       n(3,3) = 2
       n(4,3) = 2
       n(5,3) = 2
       n(6,3) = 2
       n(1,4) = 0
       n(2,4) = 1
       n(3,4) = 2
       n(4,4) = 2
       n(5,4) = 2
       n(6,4) = 2
      endif
      go to 220
c     y  hay
 60   lmx = 3
      if(olanl) then
       cc( 1, 1) =         -0.04048170d0
       cc( 2, 1) =        -20.61943440d0
       cc( 3, 1) =       -116.75222790d0
       cc( 4, 1) =        -43.79758060d0
       cc( 5, 1) =         -5.42476090d0
       a ( 1, 1) =        578.43103490d0
       a ( 2, 1) =        152.77920040d0
       a ( 3, 1) =         44.93015240d0
       a ( 4, 1) =         11.45879180d0
       a ( 5, 1) =          3.75232670d0
       n ( 1, 1) =   0
       n ( 2, 1) =   1
       n ( 3, 1) =   2
       n ( 4, 1) =   2
       n ( 5, 1) =   2
       nt( 1)    =   5
       cc( 1, 2) =          2.98013390d0
       cc( 2, 2) =         34.78346760d0
       cc( 3, 2) =         28.84532460d0
       cc( 4, 2) =         64.76420880d0
       a ( 1, 2) =         59.47151140d0
       a ( 2, 2) =         17.21735530d0
       a ( 3, 2) =         18.47970930d0
       a ( 4, 2) =          4.32761920d0
       n ( 1, 2) =   0
       n ( 2, 2) =   1
       n ( 3, 2) =   2
       n ( 4, 2) =   2
       nt( 2)    =   4
       cc( 1, 3) =          4.98857830d0
       cc( 2, 3) =         19.65065640d0
       cc( 3, 3) =        194.09431810d0
       cc( 4, 3) =         43.13497690d0
       a ( 1, 3) =         45.72710050d0
       a ( 2, 3) =         49.45958860d0
       a ( 3, 3) =         18.99523730d0
       a ( 4, 3) =          3.66031930d0
       n ( 1, 3) =   0
       n ( 2, 3) =   1
       n ( 3, 3) =   2
       n ( 4, 3) =   2
       nt( 3)    =   4
       cc( 1, 4) =          3.00666470d0
       cc( 2, 4) =         25.98792500d0
       cc( 3, 4) =         85.71728970d0
       cc( 4, 4) =         48.77925680d0
       cc( 5, 4) =         11.45351040d0
       a ( 1, 4) =         62.82684350d0
       a ( 2, 4) =         31.88979040d0
       a ( 3, 4) =         18.36465720d0
       a ( 4, 4) =          7.30624000d0
       a ( 5, 4) =          2.40516350d0
       n ( 1, 4) =   0
       n ( 2, 4) =   1
       n ( 3, 4) =   2
       n ( 4, 4) =   2
       n ( 5, 4) =   2
       nt( 4)    =   5
       nelcor    =  28
      else
       nelcor = 36
       nt(1) = 7
       nt(2) = 5
       nt(3) = 4
       nt(4) = 6
       cc(1,1) = -0.0404817d0
       cc(2,1) = -24.7037153d0
       cc(3,1) = -87.0963659d0
       cc(4,1) = -34.0507769d0
       cc(5,1) = -26.2826968d0
       cc(6,1) = -7.6489315d0
       cc(7,1) = -1.2960907d0
       cc(1,2) = 2.9379959d0
       cc(2,2) = 30.1960922d0
       cc(3,2) = 27.0994409d0
       cc(4,2) = 19.6677401d0
       cc(5,2) = -12.2337233d0
       cc(1,3) = 1.9727650d0
       cc(2,3) = 26.1215711d0
       cc(3,3) = 17.0466642d0
       cc(4,3) = 4.9537127d0
       cc(1,4) = 2.9837545d0
       cc(2,4) = 22.8182359d0
       cc(3,4) = 82.1370962d0
       cc(4,4) = 41.9336967d0
       cc(5,4) = -31.2567263d0
       cc(6,4) = -0.1926108d0
c
       a(1,1) = 286.2468288d0
       a(2,1) = 79.7488604d0
       a(3,1) = 26.9353388d0
       a(4,1) = 12.8167805d0
       a(5,1) = 5.5488178d0
       a(6,1) = 1.5019950d0
       a(7,1) = 0.5336407d0
       a(1,2) = 28.5474197d0
       a(2,2) = 8.4334903d0
       a(3,2) = 2.7123344d0
       a(4,2) = 0.4715347d0
       a(5,2) = 0.4358282d0
       a(1,3) = 22.9190664d0
       a(2,3) = 6.7622487d0
       a(3,3) = 1.8043968d0
       a(4,3) = 0.3843715d0
       a(1,4) = 42.4281313d0
       a(2,4) = 24.1184707d0
       a(3,4) = 10.3891093d0
       a(4,4) = 1.6323830d0
       a(5,4) = 1.5072156d0
       a(6,4) = 0.2580772d0
       n(1,1) = 0
       n(2,1) = 1
       n(3,1) = 2
       n(4,1) = 2
       n(5,1) = 2
       n(6,1) = 2
       n(7,1) = 2
       n(1,2) = 0
       n(2,2) = 1
       n(3,2) = 2
       n(4,2) = 2
       n(5,2) = 2
       n(1,3) = 0
       n(2,3) = 1
       n(3,3) = 2
       n(4,3) = 2
       n(1,4) = 0
       n(2,4) = 1
       n(3,4) = 2
       n(4,4) = 2
       n(5,4) = 2
       n(6,4) = 2
      endif
      go to 220
c
c     zr hay
c
 70   lmx = 3
      if(olanl) then
       cc( 1, 1) =         -0.04258430d0
       cc( 2, 1) =        -20.22224090d0
       cc( 3, 1) =       -101.86951720d0
       cc( 4, 1) =        -41.61957840d0
       cc( 5, 1) =         -4.69861600d0
       a ( 1, 1) =        645.93218730d0
       a ( 2, 1) =        134.75474010d0
       a ( 3, 1) =         42.30746190d0
       a ( 4, 1) =         12.00032270d0
       a ( 5, 1) =          4.12604540d0
       n ( 1, 1) =   0
       n ( 2, 1) =   1
       n ( 3, 1) =   2
       n ( 4, 1) =   2
       n ( 5, 1) =   2
       nt( 1)    =   5
       cc( 1, 2) =          2.79105590d0
       cc( 2, 2) =         41.94894590d0
       cc( 3, 2) =         67.72718660d0
       a ( 1, 2) =        117.62518620d0
       a ( 2, 2) =         22.96460890d0
       a ( 3, 2) =          4.52252980d0
       n ( 1, 2) =   0
       n ( 2, 2) =   1
       n ( 3, 2) =   2
       nt( 2)    =   3
       cc( 1, 3) =          4.99111440d0
       cc( 2, 3) =         20.71931720d0
       cc( 3, 3) =        195.58677580d0
       cc( 4, 3) =         48.28771760d0
       a ( 1, 3) =         47.19531450d0
       a ( 2, 3) =         48.03560330d0
       a ( 3, 3) =         19.45414560d0
       a ( 4, 3) =          4.05128750d0
       n ( 1, 3) =   0
       n ( 2, 3) =   1
       n ( 3, 3) =   2
       n ( 4, 3) =   2
       nt( 3)    =   4
       cc( 1, 4) =          3.00492260d0
       cc( 2, 4) =         25.93779890d0
       cc( 3, 4) =        125.12449340d0
       cc( 4, 4) =         70.76340220d0
       cc( 5, 4) =         15.04928220d0
       a ( 1, 4) =         79.90739830d0
       a ( 2, 4) =         45.82637980d0
       a ( 3, 4) =         26.99035220d0
       a ( 4, 4) =          9.68357180d0
       a ( 5, 4) =          2.79956660d0
       n ( 1, 4) =   0
       n ( 2, 4) =   1
       n ( 3, 4) =   2
       n ( 4, 4) =   2
       n ( 5, 4) =   2
       nt( 4)    =   5
       nelcor    =  28
      else
       nelcor = 36
       nt(1) = 7
       nt(2) = 5
       nt(3) = 4
       nt(4) = 6
       cc(1,1) = -0.0425843d0
       cc(2,1) = -24.6095772d0
       cc(3,1) = -88.1416208d0
       cc(4,1) = -30.5867042d0
       cc(5,1) = -25.5132044d0
       cc(6,1) = -7.5074683d0
       cc(7,1) = -0.9316733d0
       cc(1,2) = 2.9347023d0
       cc(2,2) = 31.7821845d0
       cc(3,2) = 35.4651631d0
       cc(4,2) = 16.6848086d0
       cc(5,2) = -8.1203706d0
       cc(1,3) = 1.9713366d0
       cc(2,3) = 27.0213041d0
       cc(3,3) = 20.2473609d0
       cc(4,3) = 5.9918936d0
       cc(1,4) = 2.9829079d0
       cc(2,4) = 21.2100150d0
       cc(3,4) = 110.0416961d0
       cc(4,4) = 40.6393744d0
       cc(5,4) = -27.1743424d0
       cc(6,4) = -0.2148958d0
       a(1,1) = 258.6887844d0
       a(2,1) = 79.7997982d0
       a(3,1) = 26.3866691d0
       a(4,1) = 12.8834065d0
       a(5,1) = 5.6371132d0
       a(6,1) = 1.5212836d0
       a(7,1) = 0.5288955d0
       a(1,2) = 37.0137467d0
       a(2,2) = 10.7740348d0
       a(3,2) = 3.5716580d0
       a(4,2) = 0.5152143d0
       a(5,2) = 0.4498238d0
       a(1,3) = 28.2944527d0
       a(2,3) = 8.2761325d0
       a(3,3) = 2.3760900d0
       a(4,3) = 0.4378428d0
       a(1,4) = 44.5032208d0
       a(2,4) = 33.8163010d0
       a(3,4) = 13.4383705d0
       a(4,4) = 1.8706775d0
       a(5,4) = 1.6691343d0
       a(6,4) = 0.2966756d0
       n(1,1) = 0
       n(2,1) = 1
       n(3,1) = 2
       n(4,1) = 2
       n(5,1) = 2
       n(6,1) = 2
       n(7,1) = 2
       n(1,2) = 0
       n(2,2) = 1
       n(3,2) = 2
       n(4,2) = 2
       n(5,2) = 2
       n(1,3) = 0
       n(2,3) = 1
       n(3,3) = 2
       n(4,3) = 2
       n(1,4) = 0
       n(2,4) = 1
       n(3,4) = 2
       n(4,4) = 2
       n(5,4) = 2
       n(6,4) = 2
      endif
      go to 220
c
c     nb hay
c
 80   lmx = 3
      if(olanl) then
       cc( 1, 1) =         -0.04474010d0
       cc( 2, 1) =        -20.05351000d0
       cc( 3, 1) =       -103.72440210d0
       cc( 4, 1) =        -41.04591650d0
       cc( 5, 1) =         -4.20462190d0
       a ( 1, 1) =        342.96384050d0
       a ( 2, 1) =        139.93080140d0
       a ( 3, 1) =         43.35234050d0
       a ( 4, 1) =         12.44485330d0
       a ( 5, 1) =          4.32048370d0
       n ( 1, 1) =   0
       n ( 2, 1) =   1
       n ( 3, 1) =   2
       n ( 4, 1) =   2
       n ( 5, 1) =   2
       nt( 1)    =   5
       cc( 1, 2) =          2.79646690d0
       cc( 2, 2) =         42.88457860d0
       cc( 3, 2) =         75.18769660d0
       a ( 1, 2) =        112.86306170d0
       a ( 2, 2) =         22.96786760d0
       a ( 3, 2) =          4.93406730d0
       n ( 1, 2) =   0
       n ( 2, 2) =   1
       n ( 3, 2) =   2
       nt( 2)    =   3
       cc( 1, 3) =          4.94618460d0
       cc( 2, 3) =         24.83163760d0
       cc( 3, 3) =        137.91832920d0
       cc( 4, 3) =         50.97784710d0
       a ( 1, 3) =         63.68019630d0
       a ( 2, 3) =         23.09124290d0
       a ( 3, 3) =         24.42836470d0
       a ( 4, 3) =          4.23100900d0
       n ( 1, 3) =   0
       n ( 2, 3) =   1
       n ( 3, 3) =   2
       n ( 4, 3) =   2
       nt( 3)    =   4
       cc( 1, 4) =          3.00640590d0
       cc( 2, 4) =         25.53474580d0
       cc( 3, 4) =        178.95974520d0
       cc( 4, 4) =         92.95774900d0
       cc( 5, 4) =         18.47442560d0
       a ( 1, 4) =         99.33596360d0
       a ( 2, 4) =         64.05864720d0
       a ( 3, 4) =         37.08910110d0
       a ( 4, 4) =         12.07085740d0
       a ( 5, 4) =          3.15753230d0
       n ( 1, 4) =   0
       n ( 2, 4) =   1
       n ( 3, 4) =   2
       n ( 4, 4) =   2
       n ( 5, 4) =   2
       nt( 4)    =   5
       nelcor    =  28
       else
       nelcor = 36
       nt(1) = 7
       nt(2) = 5
       nt(3) = 4
       nt(4) = 6
       cc(1,1) = -0.0447401d0
       cc(2,1) = -24.8992092d0
       cc(3,1) = -92.9999902d0
       cc(4,1) = -25.9973954d0
       cc(5,1) = -25.2799818d0
       cc(6,1) = -7.0961470d0
       cc(7,1) = -0.6683218d0
       cc(1,2) = 2.9313174d0
       cc(2,2) = 32.4649833d0
       cc(3,2) = 33.7155579d0
       cc(4,2) = 17.7141979d0
       cc(5,2) = -7.6954542d0
       cc(1,3) = 1.9698707d0
       cc(2,3) = 28.0614767d0
       cc(3,3) = 21.5315226d0
       cc(4,3) = 6.7775638d0
       cc(1,4) = 2.9820394d0
       cc(2,4) = 19.7329038d0
       cc(3,4) = 147.4617485d0
       cc(4,4) = 36.6351631d0
       cc(5,4) = -20.3584460d0
       cc(6,4) = -0.1931499d0
c
       a(1,1) = 486.2081260d0
       a(2,1) = 82.6257910d0
       a(3,1) = 26.0707013d0
       a(4,1) = 13.0525293d0
       a(5,1) = 5.6574876d0
       a(6,1) = 1.5285583d0
       a(7,1) = 0.5170244d0
       a(1,2) = 34.0020679d0
       a(2,2) = 10.0551965d0
       a(3,2) = 3.2855386d0
       a(4,2) = 0.5885171d0
       a(5,2) = 0.5033570d0
       a(1,3) = 26.6160240d0
       a(2,3) = 7.8829439d0
       a(3,3) = 2.1744106d0
       a(4,3) = 0.4978172d0
       a(1,4) = 47.7396406d0
       a(2,4) = 52.1900411d0
       a(3,4) = 17.6960849d0
       a(4,4) = 2.0973242d0
       a(5,4) = 1.7470857d0
       a(6,4) = 0.3078355d0
       n(1,1) = 0
       n(2,1) = 1
       n(3,1) = 2
       n(4,1) = 2
       n(5,1) = 2
       n(6,1) = 2
       n(7,1) = 2
       n(1,2) = 0
       n(2,2) = 1
       n(3,2) = 2
       n(4,2) = 2
       n(5,2) = 2
       n(1,3) = 0
       n(2,3) = 1
       n(3,3) = 2
       n(4,3) = 2
       n(1,4) = 0
       n(2,4) = 1
       n(3,4) = 2
       n(4,4) = 2
       n(5,4) = 2
       n(6,4) = 2
      endif
      go to 220
c
c     mo hay
c
 90   lmx = 3
      if(olanl) then
       cc( 1, 1) =         -0.04694920d0
       cc( 2, 1) =        -20.20800840d0
       cc( 3, 1) =       -106.21163020d0
       cc( 4, 1) =        -41.81073680d0
       cc( 5, 1) =         -4.20541030d0
       a ( 1, 1) =        537.96678070d0
       a ( 2, 1) =        147.89829380d0
       a ( 3, 1) =         45.73588980d0
       a ( 4, 1) =         13.29114670d0
       a ( 5, 1) =          4.70599610d0
       n ( 1, 1) =   0
       n ( 2, 1) =   1
       n ( 3, 1) =   2
       n ( 4, 1) =   2
       n ( 5, 1) =   2
       nt( 1)    =   5
       cc( 1, 2) =          2.80637170d0
       cc( 2, 2) =         44.51620120d0
       cc( 3, 2) =         82.77852270d0
       a ( 1, 2) =        110.29917600d0
       a ( 2, 2) =         23.20146450d0
       a ( 3, 2) =          5.35301310d0
       n ( 1, 2) =   0
       n ( 2, 2) =   1
       n ( 3, 2) =   2
       nt( 2)    =   3
       cc( 1, 3) =          4.94208760d0
       cc( 2, 3) =         25.86049760d0
       cc( 3, 3) =        132.47087420d0
       cc( 4, 3) =         57.31497940d0
       a ( 1, 3) =         63.29013970d0
       a ( 2, 3) =         23.33153020d0
       a ( 3, 3) =         24.67594230d0
       a ( 4, 3) =          4.64930400d0
       n ( 1, 3) =   0
       n ( 2, 3) =   1
       n ( 3, 3) =   2
       n ( 4, 3) =   2
       nt( 3)    =   4
       cc( 1, 4) =          3.00545910d0
       cc( 2, 4) =         26.36378510d0
       cc( 3, 4) =        183.38491990d0
       cc( 4, 4) =         98.44530680d0
       cc( 5, 4) =         22.49013770d0
       a ( 1, 4) =        104.48399770d0
       a ( 2, 4) =         66.23072450d0
       a ( 3, 4) =         39.12831760d0
       a ( 4, 4) =         13.11644370d0
       a ( 5, 4) =          3.62802630d0
       n ( 1, 4) =   0
       n ( 2, 4) =   1
       n ( 3, 4) =   2
       n ( 4, 4) =   2
       n ( 5, 4) =   2
       nt( 4)    =   5
       nelcor    =  28
      else
       nelcor = 36
       nt(1) = 7
       nt(2) = 5
       nt(3) = 4
       nt(4) = 6
       cc(1,1) = -0.0469492d0
       cc(2,1) = -24.9754989d0
       cc(3,1) = -104.0253691d0
       cc(4,1) = -18.6448797d0
       cc(5,1) = -26.2555763d0
       cc(6,1) = -6.5246802d0
       cc(7,1) = -0.4606047d0
       cc(1,2) = 2.9278406d0
       cc(2,2) = 34.3483716d0
       cc(3,2) = 35.5827749d0
       cc(4,2) = 17.6822057d0
       cc(5,2) = -6.2317738d0
       cc(1,3) = 1.9683670d0
       cc(2,3) = 30.3775276d0
       cc(3,3) = 25.5509616d0
       cc(4,3) = 7.7414851d0
       cc(1,4) = 2.9811493d0
       cc(2,4) = 21.3037310d0
       cc(3,4) = 136.4686418d0
       cc(4,4) = 36.6796104d0
       cc(5,4) = -18.2198346d0
       cc(6,4) = -0.2024621d0
c
       a(1,1) = 140.4577691d0
       a(2,1) = 89.4739342d0
       a(3,1) = 26.1462922d0
       a(4,1) = 13.3034210d0
       a(5,1) = 5.7109321d0
       a(6,1) = 1.5161163d0
       a(7,1) = 0.4925250d0
       a(1,2) = 33.7771969d0
       a(2,2) = 10.0120020d0
       a(3,2) = 3.2351392d0
       a(4,2) = 0.6672350d0
       a(5,2) = 0.5518063d0
       a(1,3) = 27.8974995d0
       a(2,3) = 8.2676943d0
       a(3,3) = 2.2776224d0
       a(4,3) = 0.5624478d0
       a(1,4) = 50.9491474d0
       a(2,4) = 43.9395130d0
       a(3,4) = 16.9218563d0
       a(4,4) = 2.4215559d0
       a(5,4) = 1.9637047d0
       a(6,4) = 0.3327274d0
       n(1,1) = 0
       n(2,1) = 1
       n(3,1) = 2
       n(4,1) = 2
       n(5,1) = 2
       n(6,1) = 2
       n(7,1) = 2
       n(1,2) = 0
       n(2,2) = 1
       n(3,2) = 2
       n(4,2) = 2
       n(5,2) = 2
       n(1,3) = 0
       n(2,3) = 1
       n(3,3) = 2
       n(4,3) = 2
       n(1,4) = 0
       n(2,4) = 1
       n(3,4) = 2
       n(4,4) = 2
       n(5,4) = 2
       n(6,4) = 2
      endif
      go to 220
c
c     tc hay
c
 100  lmx = 3
      if(olanl) then
       cc( 1, 1) =         -0.04921150d0
       cc( 2, 1) =        -20.14673080d0
       cc( 3, 1) =       -105.57928800d0
       cc( 4, 1) =        -41.59011970d0
       cc( 5, 1) =         -3.94151110d0
       a ( 1, 1) =        501.11994090d0
       a ( 2, 1) =        150.24688440d0
       a ( 3, 1) =         46.70979140d0
       a ( 4, 1) =         13.93390140d0
       a ( 5, 1) =          4.98540690d0
       n ( 1, 1) =   0
       n ( 2, 1) =   1
       n ( 3, 1) =   2
       n ( 4, 1) =   2
       n ( 5, 1) =   2
       nt( 1)    =   5
       cc( 1, 2) =          2.82223740d0
       cc( 2, 2) =         46.73307080d0
       cc( 3, 2) =         90.53186570d0
       a ( 1, 2) =        108.17759950d0
       a ( 2, 2) =         23.65943510d0
       a ( 3, 2) =          5.78149770d0
       n ( 1, 2) =   0
       n ( 2, 2) =   1
       n ( 3, 2) =   2
       nt( 2)    =   3
       cc( 1, 3) =          4.93907820d0
       cc( 2, 3) =         26.92296310d0
       cc( 3, 3) =        128.50797960d0
       cc( 4, 3) =         63.74902500d0
       a ( 1, 3) =         63.30620140d0
       a ( 2, 3) =         23.69355320d0
       a ( 3, 3) =         25.05258540d0
       a ( 4, 3) =          5.07170840d0
       n ( 1, 3) =   0
       n ( 2, 3) =   1
       n ( 3, 3) =   2
       n ( 4, 3) =   2
       nt( 3)    =   4
       cc( 1, 4) =          2.99219020d0
       cc( 2, 4) =         28.18213790d0
       cc( 3, 4) =         89.05484650d0
       cc( 4, 4) =        121.18943250d0
       cc( 5, 4) =         25.58669730d0
       a ( 1, 4) =        102.77213970d0
       a ( 2, 4) =         45.67137560d0
       a ( 3, 4) =         48.26119970d0
       a ( 4, 4) =         15.70891150d0
       a ( 5, 4) =          3.94845760d0
       n ( 1, 4) =   0
       n ( 2, 4) =   1
       n ( 3, 4) =   2
       n ( 4, 4) =   2
       n ( 5, 4) =   2
       nt( 4)    =   5
       nelcor    =  28
      else
       nelcor = 36
       nt(1) = 7
       nt(2) = 5
       nt(3) = 4
       nt(4) = 6
       cc(1,1) = -0.0492115d0
       cc(2,1) = -24.8746952d0
       cc(3,1) = -94.3425609d0
       cc(4,1) = -19.2541442d0
       cc(5,1) = -23.5804838d0
       cc(6,1) = -6.3060155d0
       cc(7,1) = -0.3784024d0
       cc(1,2) = 2.9242712d0
       cc(2,2) = 33.2684200d0
       cc(3,2) = 38.3154372d0
       cc(4,2) = 22.5461021d0
       cc(5,2) = -9.0468255d0
       cc(1,3) = 1.9668255d0
       cc(2,3) = 29.3737383d0
       cc(3,3) = 25.8622165d0
       cc(4,3) = 8.4254815d0
       cc(1,4) = 2.9802373d0
       cc(2,4) = 20.0497296d0
       cc(3,4) = 184.0269963d0
       cc(4,4) = 52.1621823d0
       cc(5,4) = -30.6237885d0
       cc(6,4) = -0.1767003d0
c
       a(1,1) = 303.7613257d0
       a(2,1) = 82.2486843d0
       a(3,1) = 24.8708254d0
       a(4,1) = 13.2618547d0
       a(5,1) = 5.6217429d0
       a(6,1) = 1.5852955d0
       a(7,1) = 0.4794448d0
       a(1,2) = 38.6719811d0
       a(2,2) = 11.7892330d0
       a(3,2) = 3.9915876d0
       a(4,2) = 0.7222026d0
       a(5,2) = 0.5932238d0
       a(1,3) = 32.8241364d0
       a(2,3) = 9.7352606d0
       a(3,3) = 2.8126853d0
       a(4,3) = 0.6071612d0
       a(1,4) = 60.6161543d0
       a(2,4) = 66.0995962d0
       a(3,4) = 22.6740462d0
       a(4,4) = 2.5392036d0
       a(5,4) = 2.1363354d0
       a(6,4) = 0.3435637d0
       n(1,1) = 0
       n(2,1) = 1
       n(3,1) = 2
       n(4,1) = 2
       n(5,1) = 2
       n(6,1) = 2
       n(7,1) = 2
       n(1,2) = 0
       n(2,2) = 1
       n(3,2) = 2
       n(4,2) = 2
       n(5,2) = 2
       n(1,3) = 0
       n(2,3) = 1
       n(3,3) = 2
       n(4,3) = 2
       n(1,4) = 0
       n(2,4) = 1
       n(3,4) = 2
       n(4,4) = 2
       n(5,4) = 2
       n(6,4) = 2
      endif
      go to 220
c
c     ru hay
c
 110  lmx = 3
      if(olanl) then
       cc( 1, 1) =         -0.05152700d0
       cc( 2, 1) =        -20.18165360d0
       cc( 3, 1) =       -105.99669150d0
       cc( 4, 1) =        -42.21667880d0
       cc( 5, 1) =         -3.76750240d0
       a ( 1, 1) =        554.37963030d0
       a ( 2, 1) =        155.10668710d0
       a ( 3, 1) =         48.49762630d0
       a ( 4, 1) =         14.77015940d0
       a ( 5, 1) =          5.20773630d0
       n ( 1, 1) =   0
       n ( 2, 1) =   1
       n ( 3, 1) =   2
       n ( 4, 1) =   2
       n ( 5, 1) =   2
       nt( 1)    =   5
       cc( 1, 2) =          2.95783440d0
       cc( 2, 2) =         25.37487070d0
       cc( 3, 2) =        536.12623720d0
       cc( 4, 2) =       -651.20572210d0
       cc( 5, 2) =        381.38169430d0
       a ( 1, 2) =         66.71180600d0
       a ( 2, 2) =         77.35036320d0
       a ( 3, 2) =         18.35714450d0
       a ( 4, 2) =         11.84047270d0
       a ( 5, 2) =          8.11794790d0
       n ( 1, 2) =   0
       n ( 2, 2) =   1
       n ( 3, 2) =   2
       n ( 4, 2) =   2
       n ( 5, 2) =   2
       nt( 2)    =   5
       cc( 1, 3) =          4.96515570d0
       cc( 2, 3) =         23.88615010d0
       cc( 3, 3) =        464.46313440d0
       cc( 4, 3) =       -714.44517880d0
       cc( 5, 3) =        377.55035940d0
       a ( 1, 3) =         54.99379150d0
       a ( 2, 3) =         13.93992120d0
       a ( 3, 3) =         15.21182460d0
       a ( 4, 3) =         10.54606910d0
       a ( 5, 3) =          7.55394860d0
       n ( 1, 3) =   0
       n ( 2, 3) =   1
       n ( 3, 3) =   2
       n ( 4, 3) =   2
       n ( 5, 3) =   2
       nt( 3)    =   5
       cc( 1, 4) =          3.03529880d0
       cc( 2, 4) =         23.29017230d0
       cc( 3, 4) =        146.09266200d0
       cc( 4, 4) =         28.91297700d0
       a ( 1, 4) =         60.34445950d0
       a ( 2, 4) =         45.21003050d0
       a ( 3, 4) =         19.11900740d0
       a ( 4, 4) =          4.27120900d0
       n ( 1, 4) =   0
       n ( 2, 4) =   1
       n ( 3, 4) =   2
       n ( 4, 4) =   2
       nt( 4)    =   4
       nelcor    =  28
      else
       nelcor = 36
       nt(1) = 7
       nt(2) = 5
       nt(3) = 4
       nt(4) = 6
       cc(1,1) = -0.0515270d0
       cc(2,1) = -24.2111980d0
       cc(3,1) = -90.1324868d0
       cc(4,1) = -18.7412056d0
       cc(5,1) = -21.3351883d0
       cc(6,1) = -5.9996563d0
       cc(7,1) = -0.2994102d0
c
       cc(1,2) = 2.9206085d0
       cc(2,2) = 34.7473053d0
       cc(3,2) = 41.6989147d0
       cc(4,2) = 29.8034158d0
       cc(5,2) = -16.9340057d0
c
       cc(1,3) = 1.9652461d0
       cc(2,3) = 30.4935522d0
       cc(3,3) = 28.0446947d0
       cc(4,3) = 9.0079429d0
c
       cc(1,4) = 2.9793035d0
       cc(2,4) = 20.6159067d0
       cc(3,4) = 175.8361030d0
       cc(4,4) = 64.1834035d0
       cc(5,4) = -40.1987590d0
       cc(6,4) = -0.1880375d0
c
       a(1,1) = 239.2060503d0
       a(2,1) = 77.7081665d0
       a(3,1) = 24.1792246d0
       a(4,1) = 12.9011817d0
       a(5,1) = 5.5165833d0
       a(6,1) = 1.6355769d0
       a(7,1) = 0.4519993d0
       a(1,2) = 40.1276685d0
       a(2,2) = 12.1750285d0
       a(3,2) = 4.0290193d0
       a(4,2) = 0.7120459d0
       a(5,2) = 0.6306855d0
       a(1,3) = 32.6793541d0
       a(2,3) = 9.7452864d0
       a(3,3) = 2.7434318d0
       a(4,3) = 0.6652098d0
       a(1,4) = 57.6689621d0
       a(2,4) = 62.3738213d0
       a(3,4) = 21.9924916d0
       a(4,4) = 2.8136093d0
       a(5,4) = 2.4287338d0
       a(6,4) = 0.3765219d0
       n(1,1) = 0
       n(2,1) = 1
       n(3,1) = 2
       n(4,1) = 2
       n(5,1) = 2
       n(6,1) = 2
       n(7,1) = 2
       n(1,2) = 0
       n(2,2) = 1
       n(3,2) = 2
       n(4,2) = 2
       n(5,2) = 2
       n(1,3) = 0
       n(2,3) = 1
       n(3,3) = 2
       n(4,3) = 2
       n(1,4) = 0
       n(2,4) = 1
       n(3,4) = 2
       n(4,4) = 2
       n(5,4) = 2
       n(6,4) = 2
      endif
      go to 220
c     rh hay
 120  lmx = 3
      if(olanl) then
       cc( 1, 1) =         -0.05389580d0
       cc( 2, 1) =        -20.13162820d0
       cc( 3, 1) =       -105.36541210d0
       cc( 4, 1) =        -42.32743700d0
       cc( 5, 1) =         -3.66540430d0
       a ( 1, 1) =        600.32430320d0
       a ( 2, 1) =        157.69101760d0
       a ( 3, 1) =         49.88419950d0
       a ( 4, 1) =         15.59668950d0
       a ( 5, 1) =          5.50992960d0
       n ( 1, 1) =   0
       n ( 2, 1) =   1
       n ( 3, 1) =   2
       n ( 4, 1) =   2
       n ( 5, 1) =   2
       nt( 1)    =   5
       cc( 1, 2) =          2.97537280d0
       cc( 2, 2) =         25.12303060d0
       cc( 3, 2) =        626.09261450d0
       cc( 4, 2) =       -812.25493850d0
       cc( 5, 2) =        467.37293400d0
       a ( 1, 2) =         59.34425260d0
       a ( 2, 2) =         83.74260610d0
       a ( 3, 2) =         18.45302480d0
       a ( 4, 2) =         12.41946060d0
       a ( 5, 2) =          8.81729130d0
       n ( 1, 2) =   0
       n ( 2, 2) =   1
       n ( 3, 2) =   2
       n ( 4, 2) =   2
       n ( 5, 2) =   2
       nt( 2)    =   5
       cc( 1, 3) =          4.95372130d0
       cc( 2, 3) =         20.48711160d0
       cc( 3, 3) =        598.01201390d0
       cc( 4, 3) =       -718.40590280d0
       cc( 5, 3) =        382.81731510d0
       a ( 1, 3) =         53.43090680d0
       a ( 2, 3) =         65.66718430d0
       a ( 3, 3) =         16.83698620d0
       a ( 4, 3) =         11.30421360d0
       a ( 5, 3) =          8.03124440d0
       n ( 1, 3) =   0
       n ( 2, 3) =   1
       n ( 3, 3) =   2
       n ( 4, 3) =   2
       n ( 5, 3) =   2
       nt( 3)    =   5
       cc( 1, 4) =          3.02795320d0
       cc( 2, 4) =         24.75265160d0
       cc( 3, 4) =        142.68442890d0
       cc( 4, 4) =         32.14068570d0
       a ( 1, 4) =         64.39936530d0
       a ( 2, 4) =         43.46250530d0
       a ( 3, 4) =         19.40203010d0
       a ( 4, 4) =          4.68793280d0
       n ( 1, 4) =   0
       n ( 2, 4) =   1
       n ( 3, 4) =   2
       n ( 4, 4) =   2
       nt( 4)    =   4
       nelcor    =  28
      else
       nelcor = 36
       nt(1) = 7
       nt(2) = 5
       nt(3) = 4
       nt(4) = 6
       cc(1,1) = -0.0538958d0
       cc(2,1) = -23.9307707d0
       cc(3,1) = -89.8969366d0
       cc(4,1) = -21.1724377d0
       cc(5,1) = -18.3052670d0
       cc(6,1) = -5.4081849d0
       cc(7,1) = -0.1924897d0
c
       cc(1,2) = 2.9168518d0
       cc(2,2) = 35.6832341d0
       cc(3,2) = 43.2107769d0
       cc(4,2) = 42.7779660d0
       cc(5,2) = -28.8814407d0
c
       cc(1,3) = 1.9636285d0
       cc(2,3) = 32.0601598d0
       cc(3,3) = 31.4780808d0
       cc(4,3) = 9.5877511d0
c
       cc(1,4) = 2.9783479d0
       cc(2,4) = 21.4067434d0
       cc(3,4) = 167.9848402d0
       cc(4,4) = 81.6117137d0
       cc(5,4) = -55.3483688d0
       cc(6,4) = -0.1868933d0
c
       a(1,1) = 317.3059724d0
       a(2,1) = 77.6959662d0
       a(3,1) = 24.7489753d0
       a(4,1) = 12.0137393d0
       a(5,1) = 5.2096064d0
       a(6,1) = 1.6235995d0
       a(7,1) = 0.3788652d0
       a(1,2) = 39.8340266d0
       a(2,2) = 12.2387614d0
       a(3,2) = 4.0240217d0
       a(4,2) = 0.7718550d0
       a(5,2) = 0.7129208d0
       a(1,3) = 33.1563494d0
       a(2,3) = 9.9603336d0
       a(3,3) = 2.7873484d0
       a(4,3) = 0.7239474d0
       a(1,4) = 58.3209808d0
       a(2,4) = 56.2607768d0
       a(3,4) = 21.4571280d0
       a(4,4) = 3.0690822d0
       a(5,4) = 2.7220803d0
       a(6,4) = 0.4010151d0
       n(1,1) = 0
       n(2,1) = 1
       n(3,1) = 2
       n(4,1) = 2
       n(5,1) = 2
       n(6,1) = 2
       n(7,1) = 2
       n(1,2) = 0
       n(2,2) = 1
       n(3,2) = 2
       n(4,2) = 2
       n(5,2) = 2
       n(1,3) = 0
       n(2,3) = 1
       n(3,3) = 2
       n(4,3) = 2
       n(1,4) = 0
       n(2,4) = 1
       n(3,4) = 2
       n(4,4) = 2
       n(5,4) = 2
       n(6,4) = 2
      endif
      go to 220
c     palladium (4d09 (5s)1  hay  ecp : pd01
 130  lmx = 3
      if(olanl) then
       cc( 1, 1) =         -0.05631770d0
       cc( 2, 1) =        -20.12880360d0
       cc( 3, 1) =       -105.81979230d0
       cc( 4, 1) =        -42.57333450d0
       cc( 5, 1) =         -3.61650860d0
       a ( 1, 1) =        598.33364440d0
       a ( 2, 1) =        162.42982900d0
       a ( 3, 1) =         51.57147710d0
       a ( 4, 1) =         16.48882600d0
       a ( 5, 1) =          5.82876560d0
       n ( 1, 1) =   0
       n ( 2, 1) =   1
       n ( 3, 1) =   2
       n ( 4, 1) =   2
       n ( 5, 1) =   2
       nt( 1)    =   5
       cc( 1, 2) =          3.00036510d0
       cc( 2, 2) =         32.43500930d0
       cc( 3, 2) =        459.08303830d0
       cc( 4, 2) =       -868.06290290d0
       cc( 5, 2) =        514.47260980d0
       a ( 1, 2) =         73.38063040d0
       a ( 2, 2) =         14.75504380d0
       a ( 3, 2) =         17.83502040d0
       a ( 4, 2) =         12.71114770d0
       a ( 5, 2) =          9.32920630d0
       n ( 1, 2) =   0
       n ( 2, 2) =   1
       n ( 3, 2) =   2
       n ( 4, 2) =   2
       n ( 5, 2) =   2
       nt( 2)    =   5
       cc( 1, 3) =          4.95930990d0
       cc( 2, 3) =         21.17110290d0
       cc( 3, 3) =        605.05600920d0
       cc( 4, 3) =       -726.96418460d0
       cc( 5, 3) =        396.32748830d0
       a ( 1, 3) =         55.66898770d0
       a ( 2, 3) =         64.23377710d0
       a ( 3, 3) =         17.62549520d0
       a ( 4, 3) =         11.90581550d0
       a ( 5, 3) =          8.51008320d0
       n ( 1, 3) =   0
       n ( 2, 3) =   1
       n ( 3, 3) =   2
       n ( 4, 3) =   2
       n ( 5, 3) =   2
       nt( 3)    =   5
       cc( 1, 4) =          3.05087450d0
       cc( 2, 4) =         22.25065800d0
       cc( 3, 4) =        674.83576980d0
       cc( 4, 4) =      -1040.85540480d0
       cc( 5, 4) =        505.93751470d0
       a ( 1, 4) =         49.99947280d0
       a ( 2, 4) =         39.74775470d0
       a ( 3, 4) =         11.43213660d0
       a ( 4, 4) =          9.17900800d0
       a ( 5, 4) =          7.56244290d0
       n ( 1, 4) =   0
       n ( 2, 4) =   1
       n ( 3, 4) =   2
       n ( 4, 4) =   2
       n ( 5, 4) =   2
       nt( 4)    =   5
       nelcor    =  28
      else
       nelcor = 36
       nt(1) = 7
       nt(2) = 5
       nt(3) = 4
       nt(4) = 6
       cc(1,1) = -0.0563177d0
       cc(2,1) = -23.8745415d0
       cc(3,1) = -89.5288200d0
       cc(4,1) = -19.2344644d0
       cc(5,1) = -17.9234918d0
       cc(6,1) = -4.8549744d0
       cc(7,1) = -0.1739145d0
       cc(1,2) = 2.9130004d0
       cc(2,2) = 36.7672813d0
       cc(3,2) = 45.7854317d0
       cc(4,2) = 59.8826421d0
       cc(5,2) = -44.8525516d0
       cc(1,3) = 1.9619728d0
       cc(2,3) = 33.9442017d0
       cc(3,3) = 34.9055940d0
       cc(4,3) = 10.1029853d0
       cc(1,4) = 2.9773705d0
       cc(2,4) = 20.2385942d0
       cc(3,4) = 232.4896834d0
       cc(4,4) = 92.9849161d0
       cc(5,4) = -62.9215422d0
       cc(6,4) = -0.1583761d0
c
       a(1,1) = 303.1163008d0
       a(2,1) = 77.5066717d0
       a(3,1) = 24.5130225d0
       a(4,1) = 12.1944317d0
       a(5,1) = 5.1391571d0
       a(6,1) = 1.6455752d0
       a(7,1) = 0.3689377d0
       a(1,2) = 39.8409904d0
       a(2,2) = 12.4879019d0
       a(3,2) = 4.1135057d0
       a(4,2) = 0.8433201d0
       a(5,2) = 0.7989471d0
       a(1,3) = 34.0332507d0
       a(2,3) = 10.1800074d0
       a(3,3) = 2.8740853d0
       a(4,3) = 0.7811346d0
       a(1,4) = 73.7537338d0
       a(2,4) = 88.4201466d0
       a(3,4) = 29.4935149d0
       a(4,4) = 3.2315488d0
       a(5,4) = 2.8402123d0
       a(6,4) = 0.4241456d0
       n(1,1) = 0
       n(2,1) = 1
       n(3,1) = 2
       n(4,1) = 2
       n(5,1) = 2
       n(6,1) = 2
       n(7,1) = 2
       n(1,2) = 0
       n(2,2) = 1
       n(3,2) = 2
       n(4,2) = 2
       n(5,2) = 2
       n(1,3) = 0
       n(2,3) = 1
       n(3,3) = 2
       n(4,3) = 2
       n(1,4) = 0
       n(2,4) = 1
       n(3,4) = 2
       n(4,4) = 2
       n(5,4) = 2
       n(6,4) = 2
      endif
      go to 220
c     ag   hay
 140  lmx = 3
      if(olanl) then
       cc( 1, 1) =         -0.05879300d0
       cc( 2, 1) =        -20.11451460d0
       cc( 3, 1) =       -104.27331140d0
       cc( 4, 1) =        -40.45397870d0
       cc( 5, 1) =         -3.44200090d0
       a ( 1, 1) =        568.70062370d0
       a ( 2, 1) =        162.35790660d0
       a ( 3, 1) =         51.10257550d0
       a ( 4, 1) =         16.92058220d0
       a ( 5, 1) =          6.16695960d0
       n ( 1, 1) =   0
       n ( 2, 1) =   1
       n ( 3, 1) =   2
       n ( 4, 1) =   2
       n ( 5, 1) =   2
       nt( 1)    =   5
       cc( 1, 2) =          2.98615270d0
       cc( 2, 2) =         35.15764600d0
       cc( 3, 2) =        450.18099060d0
       cc( 4, 2) =       -866.02483080d0
       cc( 5, 2) =        523.11101760d0
       a ( 1, 2) =         76.09746580d0
       a ( 2, 2) =         15.33273590d0
       a ( 3, 2) =         18.77153450d0
       a ( 4, 2) =         13.36632940d0
       a ( 5, 2) =          9.82369480d0
       n ( 1, 2) =   0
       n ( 2, 2) =   1
       n ( 3, 2) =   2
       n ( 4, 2) =   2
       n ( 5, 2) =   2
       nt( 2)    =   5
       cc( 1, 3) =          4.96406710d0
       cc( 2, 3) =         21.50282190d0
       cc( 3, 3) =        546.02754530d0
       cc( 4, 3) =       -600.38225560d0
       cc( 5, 3) =        348.29492890d0
       a ( 1, 3) =         56.33180430d0
       a ( 2, 3) =         69.06090980d0
       a ( 3, 3) =         19.27179980d0
       a ( 4, 3) =         12.57706540d0
       a ( 5, 3) =          8.79566700d0
       n ( 1, 3) =   0
       n ( 2, 3) =   1
       n ( 3, 3) =   2
       n ( 4, 3) =   2
       n ( 5, 3) =   2
       nt( 3)    =   5
       cc( 1, 4) =          3.04674860d0
       cc( 2, 4) =         23.36567050d0
       cc( 3, 4) =        777.25401170d0
       cc( 4, 4) =      -1238.86024230d0
       cc( 5, 4) =        608.06771210d0
       a ( 1, 4) =         53.46410780d0
       a ( 2, 4) =         40.19754570d0
       a ( 3, 4) =         11.90860730d0
       a ( 4, 4) =          9.75281830d0
       a ( 5, 4) =          8.17889970d0
       n ( 1, 4) =   0
       n ( 2, 4) =   1
       n ( 3, 4) =   2
       n ( 4, 4) =   2
       n ( 5, 4) =   2
       nt( 4)    =   5
       nelcor    =  28
c
      else
       nelcor = 36
       nt(1) = 6
       nt(2) = 5
       nt(3) = 5
       nt(4) = 6
       cc(1,1) = -0.0587930d0
       cc(2,1) = -23.2562332d0
       cc(3,1) = -92.3272663d0
       cc(4,1) = -20.4374147d0
       cc(5,1) = -5.5519951d0
       cc(6,1) = -0.6687210d0
       cc(1,2) = 2.9090535d0
       cc(2,2) = 37.4068186d0
       cc(3,2) = 113.5661277d0
       cc(4,2) = -142.2229207d0
       cc(5,2) = 81.9775247d0
       cc(1,3) = 1.9602788d0
       cc(2,3) = 29.8658485d0
       cc(3,3) = 85.7999036d0
       cc(4,3) = 32.8323292d0
       cc(5,3) = 8.5669483d0
       cc(1,4) = 2.9763712d0
       cc(2,4) = 24.1109965d0
       cc(3,4) = 153.6079747d0
       cc(4,4) = 77.1773219d0
       cc(5,4) = -46.1904997d0
       cc(6,4) = -0.1981531d0
c
       a(1,1) = 2.3171322d0
       a(2,1) = 73.3524297d0
       a(3,1) = 21.7187945d0
       a(4,1) = 5.9968258d0
       a(5,1) = 2.1258601d0
       a(6,1) = 0.7019195d0
       a(1,2) = 36.5725995d0
       a(2,2) = 11.0710579d0
       a(3,2) = 2.6697805d0
       a(4,2) = 1.9537679d0
       a(5,2) = 1.4890121d0
       a(1,3) = 5.6629108d0
       a(2,3) = 42.3303327d0
       a(3,3) = 12.7053187d0
       a(4,3) = 2.9383016d0
       a(5,3) = 0.7751946d0
       a(1,4) = 63.7701165d0
       a(2,4) = 47.3607192d0
       a(3,4) = 20.8076266d0
       a(4,4) = 3.7440069d0
       a(5,4) = 3.2144028d0
       a(6,4) = 0.4838792d0
       n(1,1) = 0
       n(2,1) = 1
       n(3,1) = 2
       n(4,1) = 2
       n(5,1) = 2
       n(6,1) = 2
       n(1,2) = 0
       n(2,2) = 1
       n(3,2) = 2
       n(4,2) = 2
       n(5,2) = 2
       n(1,3) = 0
       n(2,3) = 1
       n(3,3) = 2
       n(4,3) = 2
       n(5,3) = 2
       n(1,4) = 0
       n(2,4) = 1
       n(3,4) = 2
       n(4,4) = 2
       n(5,4) = 2
       n(6,4) = 2
      endif
      go to 220
c     cd   hay (Cd_LANL2DZ_ECP)
 150  lmx = 3
      nelcor = 36
      nt(1) = 7
      nt(2) = 5
      nt(3) = 4
      nt(4) = 6
      cc(1,1) = -0.0613214d0
      cc(2,1) = -24.9355593d0
      cc(3,1) = -87.7728765d0
      cc(4,1) = -19.0665030d0
      cc(5,1) = -16.1901254d0
      cc(6,1) = -3.9441254d0
      cc(7,1) = -0.3097476d0
      cc(1,2) = 2.9050102d0
      cc(2,2) = 39.6289198d0
      cc(3,2) = 57.4435790d0
      cc(4,2) = 90.2351166d0
      cc(5,2) = -71.7139236d0
      cc(1,3) = 1.9585463d0
      cc(2,3) = 37.5669844d0
      cc(3,3) = 44.0875104d0
      cc(4,3) = 12.5678599d0
      cc(1,4) = 2.9753499d0
      cc(2,4) = 24.7054373d0
      cc(3,4) = 218.2163143d0
      cc(4,4) = 118.9945370d0
      cc(5,4) = -82.4999671d0
      cc(6,4) = -0.2413758d0
c
      a(1,1) = 350.0762490d0
      a(2,1) = 81.1757778d0
      a(3,1) = 24.7747818d0
      a(4,1) = 12.5760330d0
      a(5,1) = 4.9889712d0
      a(6,1) = 1.8449910d0
      a(7,1) = 0.5309494d0
      a(1,2) = 54.5786193d0
      a(2,2) = 16.0082316d0
      a(3,2) = 5.2767462d0
      a(4,2) = 0.9731070d0
      a(5,2) = 0.9305646d0
      a(1,3) = 41.3581766d0
      a(2,3) = 12.1214982d0
      a(3,3) = 3.4478806d0
      a(4,3) = 0.9362175d0
      a(1,4) = 85.3710499d0
      a(2,4) = 73.0500781d0
      a(3,4) = 28.4262887d0
      a(4,4) = 3.9450158d0
      a(5,4) = 3.5349379d0
      a(6,4) = 0.6482194d0
      n(1,1) = 0
      n(2,1) = 1
      n(3,1) = 2
      n(4,1) = 2
      n(5,1) = 2
      n(6,1) = 2
      n(7,1) = 2
      n(1,2) = 0
      n(2,2) = 1
      n(3,2) = 2
      n(4,2) = 2
      n(5,2) = 2
      n(1,3) = 0
      n(2,3) = 1
      n(3,3) = 2
      n(4,3) = 2
      n(1,4) = 0
      n(2,4) = 1
      n(3,4) = 2
      n(4,4) = 2
      n(5,4) = 2
      n(6,4) = 2
      go to 220
c     in  hay  (In_LANL2DZ_ECP)
 160  lmx = 3
      nelcor = 46
      nt(1) = 5
      nt(2) = 5
      nt(3) = 5
      nt(4) = 5
      cc(1,1) = -0.0639031d0
      cc(2,1) = -29.9082901d0
      cc(3,1) = -63.3446824d0
      cc(4,1) = -19.0058428d0
      cc(5,1) = -2.5161921d0
      cc(1,2) = 3.2301607d0
      cc(2,2) = 76.7809140d0
      cc(3,2) = 245.4770760d0
      cc(4,2) = 79.2892721d0
      cc(5,2) = 26.4089460d0
      cc(1,3) = 4.4690290d0
      cc(2,3) = 43.6473399d0
      cc(3,3) = 50.7565406d0
      cc(4,3) = -48.6730017d0
      cc(5,3) = 58.7692653d0
      cc(1,4) = 3.1020770d0
      cc(2,4) = 36.5445090d0
      cc(3,4) = 48.8041808d0
      cc(4,4) = 19.4362361d0
      cc(5,4) = 3.9855647d0
c
      a(1,1) = 0.1802071d0
      a(2,1) = 33.0176528d0
      a(3,1) = 9.7219749d0
      a(4,1) = 2.3481578d0
      a(5,1) = 0.7226853d0
      a(1,2) = 173.1403228d0
      a(2,2) = 58.2752431d0
      a(3,2) = 21.1294304d0
      a(4,2) = 5.9581061d0
      a(5,2) = 1.4830272d0
      a(1,3) = 35.9438033d0
      a(2,3) = 11.2431083d0
      a(3,3) = 3.5340830d0
      a(4,3) = 1.4120937d0
      a(5,3) = 1.2828947d0
      a(1,4) = 42.6276817d0
      a(2,4) = 14.5649330d0
      a(3,4) = 5.5505378d0
      a(4,4) = 1.8351485d0
      a(5,4) = 0.5288122d0
      n(1,1) = 0
      n(2,1) = 1
      n(3,1) = 2
      n(4,1) = 2
      n(5,1) = 2
      n(1,2) = 0
      n(2,2) = 1
      n(3,2) = 2
      n(4,2) = 2
      n(5,2) = 2
      n(1,3) = 0
      n(2,3) = 1
      n(3,3) = 2
      n(4,3) = 2
      n(5,3) = 2
      n(1,4) = 0
      n(2,4) = 1
      n(3,4) = 2
      n(4,4) = 2
      n(5,4) = 2
      go to 220
c     sn   hay  (Sn_LANL2DZ_ECP)
 170  lmx = 3
      nelcor = 46
      nt(1) = 5
      nt(2) = 5
      nt(3) = 5
      nt(4) = 5
      cc(1,1) = -0.0665380d0
      cc(2,1) = -29.7469416d0
      cc(3,1) = -61.8307278d0
      cc(4,1) = -18.9410634d0
      cc(5,1) = -2.2384982d0
      cc(1,2) = 2.9482777d0
      cc(2,2) = 63.7614718d0
      cc(3,2) = 230.5276788d0
      cc(4,2) = 85.0713314d0
      cc(5,2) = 28.8215644d0
      cc(1,3) = 2.6941946d0
      cc(2,3) = 44.8271079d0
      cc(3,3) = 57.6223370d0
      cc(4,3) = -44.6853312d0
      cc(5,3) = 57.1351815d0
      cc(1,4) = 3.0259297d0
      cc(2,4) = 39.5341071d0
      cc(3,4) = 68.4686504d0
      cc(4,4) = 29.1123764d0
      cc(5,4) = 5.8428577d0
c
      a(1,1) = 0.2169021d0
      a(2,1) = 33.3097182d0
      a(3,1) = 9.8846832d0
      a(4,1) = 2.5019593d0
      a(5,1) = 0.7833323d0
      a(1,2) = 155.6062293d0
      a(2,2) = 56.5145464d0
      a(3,2) = 22.6583648d0
      a(4,2) = 6.5837080d0
      a(5,2) = 1.5947702d0
      a(1,3) = 40.4666705d0
      a(2,3) = 13.6218440d0
      a(3,3) = 4.1092688d0
      a(4,3) = 1.5046137d0
      a(5,3) = 1.3646762d0
      a(1,4) = 60.1915052d0
      a(2,4) = 20.4961142d0
      a(3,4) = 8.4438328d0
      a(4,4) = 2.6103790d0
      a(5,4) = 0.6501317d0
      n(1,1) = 0
      n(2,1) = 1
      n(3,1) = 2
      n(4,1) = 2
      n(5,1) = 2
      n(1,2) = 0
      n(2,2) = 1
      n(3,2) = 2
      n(4,2) = 2
      n(5,2) = 2
      n(1,3) = 0
      n(2,3) = 1
      n(3,3) = 2
      n(4,3) = 2
      n(5,3) = 2
      n(1,4) = 0
      n(2,4) = 1
      n(3,4) = 2
      n(4,4) = 2
      n(5,4) = 2
      go to 220
c     sb   hay (Sb_LANL2DZ_ECP)
 180  lmx = 3
      nelcor = 46
      nt(1) = 6
      nt(2) = 5
      nt(3) = 5
      nt(4) = 5
      cc(1,1) = -0.0692261d0
      cc(2,1) = -31.5434712d0
      cc(3,1) = -115.8407306d0
      cc(4,1) = -42.0224212d0
      cc(5,1) = -20.0463938d0
      cc(6,1) = -2.8796028d0
      cc(1,2) = 2.9331303d0
      cc(2,2) = 44.1104819d0
      cc(3,2) = 228.0198050d0
      cc(4,2) = 92.5610162d0
      cc(5,2) = 31.2537244d0
      cc(1,3) = 4.9947369d0
      cc(2,3) = 38.0223860d0
      cc(3,3) = 62.0520654d0
      cc(4,3) = -49.5489389d0
      cc(5,3) = 51.3666000d0
      cc(1,4) = 2.9738470d0
      cc(2,4) = 41.6495825d0
      cc(3,4) = 104.5919562d0
      cc(4,4) = 42.3187857d0
      cc(5,4) = 7.9413503d0
c
      a(1,1) = 227.2723558d0
      a(2,1) = 77.1429882d0
      a(3,1) = 26.0116501d0
      a(4,1) = 9.1421159d0
      a(5,1) = 2.9281082d0
      a(6,1) = 0.9231136d0
      a(1,2) = 116.0010888d0
      a(2,2) = 59.6755065d0
      a(3,2) = 27.6135849d0
      a(4,2) = 7.8045081d0
      a(5,2) = 1.6934418d0
      a(1,3) = 27.8050796d0
      a(2,3) = 10.3929685d0
      a(3,3) = 3.6381067d0
      a(4,3) = 1.9839681d0
      a(5,3) = 1.5546562d0
      a(1,4) = 83.2183030d0
      a(2,4) = 30.3070501d0
      a(3,4) = 13.2014120d0
      a(4,4) = 3.6961611d0
      a(5,4) = 0.7711144d0
      n(1,1) = 0
      n(2,1) = 1
      n(3,1) = 2
      n(4,1) = 2
      n(5,1) = 2
      n(6,1) = 2
      n(1,2) = 0
      n(2,2) = 1
      n(3,2) = 2
      n(4,2) = 2
      n(5,2) = 2
      n(1,3) = 0
      n(2,3) = 1
      n(3,3) = 2
      n(4,3) = 2
      n(5,3) = 2
      n(1,4) = 0
      n(2,4) = 1
      n(3,4) = 2
      n(4,4) = 2
      n(5,4) = 2
      go to 220
c    te  hay ("Te_LANL2DZ_ECP")
 190  lmx = 3
      nelcor = 46
      nt(1) = 6
      nt(2) = 5
      nt(3) = 5
      nt(4) = 5
      cc(1,1) = -0.0719675d0
      cc(2,1) = -31.3740754d0
      cc(3,1) = -112.9711996d0
      cc(4,1) = -41.9599185d0
      cc(5,1) = -19.1274006d0
      cc(6,1) = -2.2917778d0
      cc(1,2) = 2.9286639d0
      cc(2,2) = 68.9626773d0
      cc(3,2) = 276.9637063d0
      cc(4,2) = 104.4701765d0
      cc(5,2) = 34.8384805d0
      cc(1,3) = 4.9732188d0
      cc(2,3) = 38.8983634d0
      cc(3,3) = 137.3139246d0
      cc(4,3) = 69.7370323d0
      cc(5,3) = 22.7487178d0
      cc(1,4) = 3.1548778d0
      cc(2,4) = 33.1153649d0
      cc(3,4) = 210.3056318d0
      cc(4,4) = -321.1963746d0
      cc(5,4) = 152.4853001d0
c
      a(1,1) = 285.0323292d0
      a(2,1) = 75.9105798d0
      a(3,1) = 26.1827859d0
      a(4,1) = 9.0635141d0
      a(5,1) = 2.9543257d0
      a(6,1) = 0.9169942d0
      a(1,2) = 188.3327504d0
      a(2,2) = 68.6306921d0
      a(3,2) = 27.4402195d0
      a(4,2) = 7.8835282d0
      a(5,2) = 1.8443362d0
      a(1,3) = 54.0822795d0
      a(2,3) = 28.5434714d0
      a(3,3) = 16.8660857d0
      a(4,3) = 5.5244176d0
      a(5,3) = 1.4625225d0
      a(1,4) = 34.5312648d0
      a(2,4) = 11.2699366d0
      a(3,4) = 2.2704361d0
      a(4,4) = 1.7766992d0
      a(5,4) = 1.4251666d0
      n(1,1) = 0
      n(2,1) = 1
      n(3,1) = 2
      n(4,1) = 2
      n(5,1) = 2
      n(6,1) = 2
      n(1,2) = 0
      n(2,2) = 1
      n(3,2) = 2
      n(4,2) = 2
      n(5,2) = 2
      n(1,3) = 0
      n(2,3) = 1
      n(3,3) = 2
      n(4,3) = 2
      n(5,3) = 2
      n(1,4) = 0
      n(2,4) = 1
      n(3,4) = 2
      n(4,4) = 2
      n(5,4) = 2
      go to 220
c     i  hay  (I_LANL2DZ_ECP)
 200  lmx = 3
      nelcor = 46
      nt(1) = 5
      nt(2) = 5
      nt(3) = 5
      nt(4) = 5
      cc(1,1) = -0.0747621d0
      cc(2,1) = -30.0811224d0
      cc(3,1) = -75.3722721d0
      cc(4,1) = -22.0563758d0
      cc(5,1) = -1.6979585d0
      cc(1,2) = 2.9380036d0
      cc(2,2) = 41.2471267d0
      cc(3,2) = 287.8680095d0
      cc(4,2) = 114.3758506d0
      cc(5,2) = 37.6547714d0
      cc(1,3) = 2.2222630d0
      cc(2,3) = 39.4090831d0
      cc(3,3) = 177.4075002d0
      cc(4,3) = 77.9889462d0
      cc(5,3) = 25.7547641d0
      cc(1,4) = 7.0524360d0
      cc(2,4) = 33.3041635d0
      cc(3,4) = 186.9453875d0
      cc(4,4) = 71.9688361d0
      cc(5,4) = 9.3630657d0
c
      a(1,1) = 1.0715702d0
      a(2,1) = 44.1936028d0
      a(3,1) = 12.9367609d0
      a(4,1) = 3.1956412d0
      a(5,1) = 0.8589806d0
      a(1,2) = 127.9202670d0
      a(2,2) = 78.6211465d0
      a(3,2) = 36.5146237d0
      a(4,2) = 9.9065681d0
      a(5,2) = 1.9420086d0
      a(1,3) = 13.0035304d0
      a(2,3) = 76.0331404d0
      a(3,3) = 24.1961684d0
      a(4,3) = 6.4053433d0
      a(5,3) = 1.5851786d0
      a(1,4) = 40.4278108d0
      a(2,4) = 28.9084375d0
      a(3,4) = 15.6268936d0
      a(4,4) = 4.1442856d0
      a(5,4) = 0.9377235d0
      n(1,1) = 0
      n(2,1) = 1
      n(3,1) = 2
      n(4,1) = 2
      n(5,1) = 2
      n(1,2) = 0
      n(2,2) = 1
      n(3,2) = 2
      n(4,2) = 2
      n(5,2) = 2
      n(1,3) = 0
      n(2,3) = 1
      n(3,3) = 2
      n(4,3) = 2
      n(5,3) = 2
      n(1,4) = 0
      n(2,4) = 1
      n(3,4) = 2
      n(4,4) = 2
      n(5,4) = 2
      go to 220
c     xe  hay (Xe_LANL2DZ_ECP)
 210  lmx = 3
      nelcor = 46
      nt(1) = 6
      nt(2) = 5
      nt(3) = 5
      nt(4) = 5
      cc(1,1) = -0.0776099d0
      cc(2,1) = -32.4332953d0
      cc(3,1) = -151.0853351d0
      cc(4,1) = -46.4977626d0
      cc(5,1) = -17.8080990d0
      cc(6,1) = -1.6481163d0
      cc(1,2) = 2.9170308d0
      cc(2,2) = 52.0620644d0
      cc(3,2) = 293.2865672d0
      cc(4,2) = 122.8610373d0
      cc(5,2) = 41.6867277d0
      cc(1,3) = 2.0968515d0
      cc(2,3) = 30.0220424d0
      cc(3,3) = 172.8757636d0
      cc(4,3) = 67.7609841d0
      cc(5,3) = 29.8970913d0
      cc(1,4) = 7.0414935d0
      cc(2,4) = 33.3391744d0
      cc(3,4) = 249.2893077d0
      cc(4,4) = 89.3576166d0
      cc(5,4) = 11.6206324d0
c
      a(1,1) = 442.5138313d0
      a(2,1) = 106.0005534d0
      a(3,1) = 32.6123743d0
      a(4,1) = 9.4739667d0
      a(5,1) = 3.1147195d0
      a(6,1) = 0.9194472d0
      a(1,2) = 148.7400724d0
      a(2,2) = 74.7687436d0
      a(3,2) = 34.1548413d0
      a(4,2) = 9.6899675d0
      a(5,2) = 2.1065947d0
      a(1,3) = 8.1495157d0
      a(2,3) = 113.9235640d0
      a(3,3) = 35.1266431d0
      a(4,3) = 8.1449579d0
      a(5,3) = 1.7187138d0
      a(1,4) = 49.4518384d0
      a(2,4) = 38.9480918d0
      a(3,4) = 20.0161598d0
      a(4,4) = 5.0614904d0
      a(5,4) = 1.0672963d0
      n(1,1) = 0
      n(2,1) = 1
      n(3,1) = 2
      n(4,1) = 2
      n(5,1) = 2
      n(6,1) = 2
      n(1,2) = 0
      n(2,2) = 1
      n(3,2) = 2
      n(4,2) = 2
      n(5,2) = 2
      n(1,3) = 0
      n(2,3) = 1
      n(3,3) = 2
      n(4,3) = 2
      n(5,3) = 2
      n(1,4) = 0
      n(2,4) = 1
      n(3,4) = 2
      n(4,4) = 2
      n(5,4) = 2
c
c     load ecp data
c
 220  return
 6010 format (1x,'*** data not found for ecp type ',a6)
      end
**==ecplb3.f
      subroutine ecplb3 (zinp,lmx,nelcor,nt,n,cc,a,olanl)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character * 8 (z)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      dimension nt(6), n(10,6), cc(10,6), a(10,6)
c
      dimension ztit(29)
      data maxtyp/29/
      data ztit/
     $ 'cs','ba','la','ce','pr','nd',
     $ 'pm','sm','eu','gd','tb','dy','ho','er','tm','yb','lu','hf',
     $ 'ta','w ','re','os','ir','pt','au','hg','tl','pb','bi'/
c
c     this routine contains a library of ecp data taken from
c     various literature sources.
c     the ecps for
c     cs,......bi are taken from hay,wadt, jcp 82 (1985) 270,284
      do 20 ityp = 1 , maxtyp
         if (zinp.eq.ztit(ityp)) go to 40
 20   continue
 30   write (iwr,6010) zinp
      call caserr2('no library available for requested ecp')
 40   go to (50,60,70,30,30,30,30,30,30,30,30,30,30,30,30,30,30,80,90,
     +       100,110,120,130,140,150,160,170,180,190) , ityp
c     cs  hay
 50   lmx = 3
      if(olanl) then
       cc( 1, 1) =         -0.08051090d0
       cc( 2, 1) =        -33.43114480d0
       cc( 3, 1) =       -193.92708870d0
       cc( 4, 1) =        -56.83014330d0
       cc( 5, 1) =        -18.65049750d0
       cc( 6, 1) =         -2.54113990d0
       a ( 1, 1) =        636.78401300d0
       a ( 2, 1) =        147.89124460d0
       a ( 3, 1) =         42.57249110d0
       a ( 4, 1) =         11.24278830d0
       a ( 5, 1) =          3.56897070d0
       a ( 6, 1) =          1.31074620d0
       n ( 1, 1) =   0
       n ( 2, 1) =   1
       n ( 3, 1) =   2
       n ( 4, 1) =   2
       n ( 5, 1) =   2
       n ( 6, 1) =   2
       nt( 1)    =   6
       cc( 1, 2) =          3.04999340d0
       cc( 2, 2) =         37.96775140d0
       cc( 3, 2) =         45.23406730d0
       a ( 1, 2) =         40.86480700d0
       a ( 2, 2) =         10.93293290d0
       a ( 3, 2) =          2.22949030d0
       n ( 1, 2) =   0
       n ( 2, 2) =   1
       n ( 3, 2) =   2
       nt( 2)    =   3
       cc( 1, 3) =          4.64069510d0
       cc( 2, 3) =         62.03135330d0
       cc( 3, 3) =        129.59094810d0
       cc( 4, 3) =         31.51868250d0
       a ( 1, 3) =        118.85981710d0
       a ( 2, 3) =         30.68479580d0
       a ( 3, 3) =         10.17629990d0
       a ( 4, 3) =          1.82728500d0
       n ( 1, 3) =   0
       n ( 2, 3) =   1
       n ( 3, 3) =   2
       n ( 4, 3) =   2
       nt( 3)    =   4
       cc( 1, 4) =          7.30544730d0
       cc( 2, 4) =         33.27857600d0
       cc( 3, 4) =         53.12969480d0
       cc( 4, 4) =         16.21188780d0
       a ( 1, 4) =         16.00817040d0
       a ( 2, 4) =          2.64882910d0
       a ( 3, 4) =         14.55556720d0
       a ( 4, 4) =          0.67036040d0
       n ( 1, 4) =   0
       n ( 2, 4) =   1
       n ( 3, 4) =   2
       n ( 4, 4) =   2
       nt( 4)    =   4
       nelcor    =  46
      else
       nelcor = 54
       nt(1) = 6
       nt(2) = 6
       nt(3) = 6
       nt(4) = 5
       cc(1,1) = -0.0805109d0
       cc(2,1) = -34.1328001d0
       cc(3,1) = -52.8122678d0
       cc(4,1) = -15.8154657d0
       cc(5,1) = -4.1663980d0
       cc(6,1) = -0.4762118d0
       cc(1,2) = 3.1442140d0
       cc(2,2) = 56.5854559d0
       cc(3,2) = 128.7603827d0
       cc(4,2) = 55.3017056d0
       cc(5,2) = 18.0902119d0
       cc(6,2) = 5.3472572d0
       cc(1,3) = 4.9674276d0
       cc(2,3) = 44.0083673d0
       cc(3,3) = 76.0877085d0
       cc(4,3) = 37.7375053d0
       cc(5,3) = 11.4693477d0
       cc(6,3) = 2.6430068d0
       cc(1,4) = 7.0615366d0
       cc(2,4) = 38.9607427d0
       cc(3,4) = 169.8032412d0
       cc(4,4) = 72.4738876d0
       cc(5,4) = 16.5343849d0
c
       a(1,1) = 96.3863793d0
       a(2,1) = 22.7391269d0
       a(3,1) = 6.3299947d0
       a(4,1) = 2.0074183d0
       a(5,1) = 0.5792657d0
       a(6,1) = 0.2281691d0
       a(1,2) = 91.9848369d0
       a(2,2) = 30.9158046d0
       a(3,2) = 12.5834123d0
       a(4,2) = 4.1623195d0
       a(5,2) = 1.3342361d0
       a(6,2) = 0.3493508d0
       a(1,3) = 39.1545608d0
       a(2,3) = 16.4708990d0
       a(3,3) = 7.7034744d0
       a(4,3) = 2.8527363d0
       a(5,3) = 0.8534501d0
       a(6,3) = 0.2253916d0
       a(1,4) = 38.5503078d0
       a(2,4) = 25.3708163d0
       a(3,4) = 13.6474569d0
       a(4,4) = 3.8458124d0
       a(5,4) = 0.7760261d0
       n(1,1) = 0
       n(2,1) = 1
       n(3,1) = 2
       n(4,1) = 2
       n(5,1) = 2
       n(6,1) = 2
       n(1,2) = 0
       n(2,2) = 1
       n(3,2) = 2
       n(4,2) = 2
       n(5,2) = 2
       n(6,2) = 2
       n(1,3) = 0
       n(2,3) = 1
       n(3,3) = 2
       n(4,3) = 2
       n(5,3) = 2
       n(6,3) = 2
       n(1,4) = 0
       n(2,4) = 1
       n(3,4) = 2
       n(4,4) = 2
       n(5,4) = 2
      endif
      go to 200
c
c     ba  hay
c
 60   lmx = 3
      if(olanl) then
       cc( 1, 1) =         -0.08346520d0
       cc( 2, 1) =        -33.32576710d0
       cc( 3, 1) =       -190.86072320d0
       cc( 4, 1) =        -55.19841720d0
       cc( 5, 1) =        -18.02363400d0
       cc( 6, 1) =         -2.19782810d0
       a ( 1, 1) =        620.96904880d0
       a ( 2, 1) =        146.36488260d0
       a ( 3, 1) =         42.32071140d0
       a ( 4, 1) =         11.21351510d0
       a ( 5, 1) =          3.69638910d0
       a ( 6, 1) =          1.31695020d0
       n ( 1, 1) =   0
       n ( 2, 1) =   1
       n ( 3, 1) =   2
       n ( 4, 1) =   2
       n ( 5, 1) =   2
       n ( 6, 1) =   2
       nt( 1)    =   6
       cc( 1, 2) =          2.81311600d0
       cc( 2, 2) =         55.30506260d0
       cc( 3, 2) =        149.55134020d0
       cc( 4, 2) =         50.50785530d0
       a ( 1, 2) =        140.36692000d0
       a ( 2, 2) =         39.14365470d0
       a ( 3, 2) =         12.95534930d0
       a ( 4, 2) =          2.43087510d0
       n ( 1, 2) =   0
       n ( 2, 2) =   1
       n ( 3, 2) =   2
       n ( 4, 2) =   2
       nt( 2)    =   4
       cc( 1, 3) =          4.91918120d0
       cc( 2, 3) =         39.38720750d0
       cc( 3, 3) =        335.07535840d0
       cc( 4, 3) =        131.25351530d0
       cc( 5, 3) =         36.31750250d0
       a ( 1, 3) =         99.49222610d0
       a ( 2, 3) =         68.27113090d0
       a ( 3, 3) =         36.15181810d0
       a ( 4, 3) =         10.10515370d0
       a ( 5, 3) =          2.02326480d0
       n ( 1, 3) =   0
       n ( 2, 3) =   1
       n ( 3, 3) =   2
       n ( 4, 3) =   2
       n ( 5, 3) =   2
       nt( 3)    =   5
       cc( 1, 4) =          2.97641400d0
       cc( 2, 4) =         46.05711410d0
       cc( 3, 4) =        117.16585880d0
       cc( 4, 4) =         54.01308150d0
       cc( 5, 4) =         15.57849060d0
       a ( 1, 4) =         95.27525360d0
       a ( 2, 4) =         34.66046080d0
       a ( 3, 4) =         15.48910400d0
       a ( 4, 4) =          5.00158950d0
       a ( 5, 4) =          1.32362660d0
       n ( 1, 4) =   0
       n ( 2, 4) =   1
       n ( 3, 4) =   2
       n ( 4, 4) =   2
       n ( 5, 4) =   2
       nt( 4)    =   5
       nelcor    =  46
      else
       nelcor = 54
       nt(1) = 6
       nt(2) = 6
       nt(3) = 6
       nt(4) = 6
       cc(1,1) = -0.0834652d0
       cc(2,1) = -34.0589206d0
       cc(3,1) = -54.3999771d0
       cc(4,1) = -15.1098299d0
       cc(5,1) = -3.6824878d0
       cc(6,1) = -0.3119913d0
       cc(1,2) = 2.8936199d0
       cc(2,2) = 48.0985554d0
       cc(3,2) = 147.1816847d0
       cc(4,2) = 68.4387389d0
       cc(5,2) = 24.2861704d0
       cc(6,2) = 7.0303066d0
       cc(1,3) = 4.9639857d0
       cc(2,3) = 44.8891930d0
       cc(3,3) = 79.7516094d0
       cc(4,3) = 40.2675214d0
       cc(5,3) = 13.7807639d0
       cc(6,3) = 3.5768992d0
       cc(1,4) = 3.0048979d0
       cc(2,4) = 44.7964631d0
       cc(3,4) = 97.1752741d0
       cc(4,4) = 46.2481775d0
       cc(5,4) = 15.1841059d0
       cc(6,4) = -0.2101319d0
c
       a(1,1) = 80.9449309d0
       a(2,1) = 23.5462009d0
       a(3,1) = 6.5042910d0
       a(4,1) = 1.9680159d0
       a(5,1) = 0.5989001d0
       a(6,1) = 0.2345986d0
       a(1,2) = 109.3469378d0
       a(2,2) = 42.0898722d0
       a(3,2) = 19.5169970d0
       a(4,2) = 6.4269427d0
       a(5,2) = 1.9232296d0
       a(6,2) = 0.4036119d0
       a(1,3) = 41.0644520d0
       a(2,3) = 17.3996358d0
       a(3,3) = 8.2896813d0
       a(4,3) = 3.1897863d0
       a(5,3) = 1.0167215d0
       a(6,3) = 0.2603415d0
       a(1,4) = 81.9326543d0
       a(2,4) = 29.4240958d0
       a(3,4) = 13.0047308d0
       a(4,4) = 4.5577279d0
       a(5,4) = 1.3412763d0
       a(6,4) = 0.2116641d0
       n(1,1) = 0
       n(2,1) = 1
       n(3,1) = 2
       n(4,1) = 2
       n(5,1) = 2
       n(6,1) = 2
       n(1,2) = 0
       n(2,2) = 1
       n(3,2) = 2
       n(4,2) = 2
       n(5,2) = 2
       n(6,2) = 2
       n(1,3) = 0
       n(2,3) = 1
       n(3,3) = 2
       n(4,3) = 2
       n(5,3) = 2
       n(6,3) = 2
       n(1,4) = 0
       n(2,4) = 1
       n(3,4) = 2
       n(4,4) = 2
       n(5,4) = 2
       n(6,4) = 2
      endif
      go to 200
c
c     la  hay
c
 70   lmx = 3
      if(olanl) then
       cc( 1, 1) =         -0.08647270d0
       cc( 2, 1) =        -33.51688030d0
       cc( 3, 1) =       -189.07098290d0
       cc( 4, 1) =        -54.24643670d0
       cc( 5, 1) =        -17.13513420d0
       cc( 6, 1) =         -1.80921850d0
       a ( 1, 1) =        740.91660340d0
       a ( 2, 1) =        146.65714990d0
       a ( 3, 1) =         42.35854130d0
       a ( 4, 1) =         11.16475000d0
       a ( 5, 1) =          3.75563430d0
       a ( 6, 1) =          1.28718270d0
       n ( 1, 1) =   0
       n ( 2, 1) =   1
       n ( 3, 1) =   2
       n ( 4, 1) =   2
       n ( 5, 1) =   2
       n ( 6, 1) =   2
       nt( 1)    =   6
       cc( 1, 2) =          2.88949430d0
       cc( 2, 2) =         61.95215920d0
       cc( 3, 2) =        164.80744900d0
       cc( 4, 2) =         55.69546090d0
       a ( 1, 2) =        125.68626650d0
       a ( 2, 2) =         37.57689550d0
       a ( 3, 2) =         12.69548690d0
       a ( 4, 2) =          2.62760820d0
       n ( 1, 2) =   0
       n ( 2, 2) =   1
       n ( 3, 2) =   2
       n ( 4, 2) =   2
       nt( 2)    =   4
       cc( 1, 3) =          4.59784120d0
       cc( 2, 3) =         71.67197360d0
       cc( 3, 3) =        175.89658960d0
       cc( 4, 3) =         37.91364040d0
       a ( 1, 3) =        179.17822470d0
       a ( 2, 3) =         42.61956610d0
       a ( 3, 3) =         13.41359570d0
       a ( 4, 3) =          2.09503920d0
       n ( 1, 3) =   0
       n ( 2, 3) =   1
       n ( 3, 3) =   2
       n ( 4, 3) =   2
       nt( 3)    =   4
       cc( 1, 4) =          2.97054960d0
       cc( 2, 4) =         46.52399460d0
       cc( 3, 4) =        137.40781760d0
       cc( 4, 4) =         63.49451830d0
       cc( 5, 4) =         18.09986540d0
       a ( 1, 4) =        104.76820640d0
       a ( 2, 4) =         40.19303400d0
       a ( 3, 4) =         18.50507980d0
       a ( 4, 4) =          5.80511850d0
       a ( 5, 4) =          1.45832200d0
       n ( 1, 4) =   0
       n ( 2, 4) =   1
       n ( 3, 4) =   2
       n ( 4, 4) =   2
       n ( 5, 4) =   2
       nt( 4)    =   5
       nelcor    =  46
      else
       nelcor = 54
       nt(1) = 6
       nt(2) = 6
       nt(3) = 6
       nt(4) = 6
       cc(1,1) = -0.0864727d0
       cc(2,1) = -33.6524567d0
       cc(3,1) = -59.5594475d0
       cc(4,1) = -15.1332618d0
       cc(5,1) = -4.8263727d0
       cc(6,1) = -0.8658915d0
       cc(1,2) = 2.8782136d0
       cc(2,2) = 47.7896578d0
       cc(3,2) = 133.8659916d0
       cc(4,2) = 66.5553143d0
       cc(5,2) = 24.1093686d0
       cc(6,2) = 6.7969476d0
       cc(1,3) = 2.0081894d0
       cc(2,3) = 49.6039386d0
       cc(3,3) = 138.8098218d0
       cc(4,3) = 60.6949688d0
       cc(5,3) = 19.5911128d0
       cc(6,3) = 5.2749298d0
       cc(1,4) = 2.9841839d0
       cc(2,4) = 38.7922490d0
       cc(3,4) = 139.0867532d0
       cc(4,4) = 63.6790783d0
       cc(5,4) = 26.8061443d0
       cc(6,4) = -13.3307262d0
       a(1,1) = 0.6403421d0
       a(2,1) = 26.9811967d0
       a(3,1) = 7.4780622d0
       a(4,1) = 2.3364873d0
       a(5,1) = 0.9083284d0
       a(6,1) = 0.2932768d0
       a(1,2) = 86.7321968d0
       a(2,2) = 36.2795884d0
       a(3,2) = 16.8958863d0
       a(4,2) = 5.7907750d0
       a(5,2) = 1.7686257d0
       a(6,2) = 0.4036790d0
       a(1,3) = 107.6316973d0
       a(2,3) = 39.2492604d0
       a(3,3) = 15.8048478d0
       a(4,3) = 4.9233797d0
       a(5,3) = 1.3795238d0
       a(6,3) = 0.2981099d0
       a(1,4) = 76.5526601d0
       a(2,4) = 39.1710741d0
       a(3,4) = 19.3646917d0
       a(4,4) = 5.7777291d0
       a(5,4) = 1.0259657d0
       a(6,4) = 0.8249295d0
       n(1,1) = 0
       n(2,1) = 1
       n(3,1) = 2
       n(4,1) = 2
       n(5,1) = 2
       n(6,1) = 2
       n(1,2) = 0
       n(2,2) = 1
       n(3,2) = 2
       n(4,2) = 2
       n(5,2) = 2
       n(6,2) = 2
       n(1,3) = 0
       n(2,3) = 1
       n(3,3) = 2
       n(4,3) = 2
       n(5,3) = 2
       n(6,3) = 2
       n(1,4) = 0
       n(2,4) = 1
       n(3,4) = 2
       n(4,4) = 2
       n(5,4) = 2
       n(6,4) = 2
      endif
      go to 200
c
c     hf hay
c
 80   lmx = 4
      if(olanl) then
       cc( 1, 1) =        -60.00000000d0
       cc( 2, 1) =       -588.49806640d0
       cc( 3, 1) =       -199.77488200d0
       cc( 4, 1) =        -77.23651210d0
       cc( 5, 1) =        -14.65325960d0
       cc( 6, 1) =         -0.67622040d0
       a ( 1, 1) =        685.55029310d0
       a ( 2, 1) =        151.80804850d0
       a ( 3, 1) =         38.05220990d0
       a ( 4, 1) =         10.20575440d0
       a ( 5, 1) =          3.43762870d0
       a ( 6, 1) =          0.89975390d0
       n ( 1, 1) =   1
       n ( 2, 1) =   2
       n ( 3, 1) =   2
       n ( 4, 1) =   2
       n ( 5, 1) =   2
       n ( 6, 1) =   2
       nt( 1)    =   6
       cc( 1, 2) =          3.00000000d0
       cc( 2, 2) =         31.33955320d0
       cc( 3, 2) =       1251.68661380d0
       cc( 4, 2) =        700.89001060d0
       cc( 5, 2) =        266.04458450d0
       cc( 6, 2) =        465.01961550d0
       cc( 7, 2) =       -369.54045680d0
       a ( 1, 2) =        380.10967250d0
       a ( 2, 2) =        724.62345050d0
       a ( 3, 2) =        281.47646030d0
       a ( 4, 2) =         89.35909880d0
       a ( 5, 2) =         22.97657730d0
       a ( 6, 2) =          3.67000350d0
       a ( 7, 2) =          3.47030000d0
       n ( 1, 2) =   0
       n ( 2, 2) =   1
       n ( 3, 2) =   2
       n ( 4, 2) =   2
       n ( 5, 2) =   2
       n ( 6, 2) =   2
       n ( 7, 2) =   2
       nt( 2)    =   7
       cc( 1, 3) =          2.00000000d0
       cc( 2, 3) =         60.27769720d0
       cc( 3, 3) =        192.21749830d0
       cc( 4, 3) =        339.64600890d0
       cc( 5, 3) =       -266.10665430d0
       a ( 1, 3) =        228.31357710d0
       a ( 2, 3) =         64.08784470d0
       a ( 3, 3) =         20.51392890d0
       a ( 4, 3) =          3.16437360d0
       a ( 5, 3) =          2.98809240d0
       n ( 1, 3) =   0
       n ( 2, 3) =   1
       n ( 3, 3) =   2
       n ( 4, 3) =   2
       n ( 5, 3) =   2
       nt( 3)    =   5
       cc( 1, 4) =          3.00000000d0
       cc( 2, 4) =         54.65767600d0
       cc( 3, 4) =        217.58114270d0
       cc( 4, 4) =        125.13630840d0
       cc( 5, 4) =         37.62031950d0
       a ( 1, 4) =        143.71023990d0
       a ( 2, 4) =         63.18591700d0
       a ( 3, 4) =         32.02274490d0
       a ( 4, 4) =         10.76749220d0
       a ( 5, 4) =          3.10032600d0
       n ( 1, 4) =   0
       n ( 2, 4) =   1
       n ( 3, 4) =   2
       n ( 4, 4) =   2
       n ( 5, 4) =   2
       nt( 4)    =   5
       cc( 1, 5) =          4.00000000d0
       cc( 2, 5) =         49.63583930d0
       cc( 3, 5) =        178.46572220d0
       cc( 4, 5) =         82.00771730d0
       cc( 5, 5) =          6.30639430d0
       a ( 1, 5) =         85.03958530d0
       a ( 2, 5) =         41.81988200d0
       a ( 3, 5) =         19.65249670d0
       a ( 4, 5) =          6.07771200d0
       a ( 5, 5) =          1.47766930d0
       n ( 1, 5) =   0
       n ( 2, 5) =   1
       n ( 3, 5) =   2
       n ( 4, 5) =   2
       n ( 5, 5) =   2
       nt( 5)    =   5
       nelcor    =  60
      else
       nelcor = 68
       nt(1) = 6
       nt(2) = 6
       nt(3) = 5
       nt(4) = 6
       nt(5) = 6
       cc(1,1) = -0.1379731d0
       cc(2,1) = -44.7429491d0
       cc(3,1) = -228.3189889d0
       cc(4,1) = -99.5814027d0
       cc(5,1) = -24.5647636d0
       cc(6,1) = -4.0694921d0
       cc(1,2) = 2.7763702d0
       cc(2,2) = 71.9510823d0
       cc(3,2) = 170.9957954d0
       cc(4,2) = 47.9762258d0
       cc(5,2) = 17.5330190d0
       cc(6,2) = -5.9424080d0
       cc(1,3) = 1.9050102d0
       cc(2,3) = 51.7968625d0
       cc(3,3) = 115.3729030d0
       cc(4,3) = 26.0735876d0
       cc(5,3) = 6.8203642d0
       cc(1,4) = 2.9441877d0
       cc(2,4) = 51.4358799d0
       cc(3,4) = 120.7577481d0
       cc(4,4) = 36.5823417d0
       cc(5,4) = -5.2595142d0
       cc(6,4) = -0.3554933d0
       cc(1,5) = 3.9603546d0
       cc(2,5) = 39.4036737d0
       cc(3,5) = 48.2203489d0
       cc(4,5) = 14.9231704d0
       cc(5,5) = 8.3700837d0
       cc(6,5) = -6.1437920d0
c
       a(1,1) = 404.3025872d0
       a(2,1) = 109.9066721d0
       a(3,1) = 37.5316905d0
       a(4,1) = 10.0447963d0
       a(5,1) = 2.9543812d0
       a(6,1) = 0.9259902d0
       a(1,2) = 98.1486452d0
       a(2,2) = 32.5467191d0
       a(3,2) = 10.6119363d0
       a(4,2) = 3.0444242d0
       a(5,2) = 0.6938215d0
       a(6,2) = 0.5910513d0
       a(1,3) = 84.4468520d0
       a(2,3) = 27.0917584d0
       a(3,3) = 8.5339666d0
       a(4,3) = 2.2452259d0
       a(5,3) = 0.4825892d0
       a(1,4) = 76.5905042d0
       a(2,4) = 27.3537558d0
       a(3,4) = 9.6261795d0
       a(4,4) = 2.5842485d0
       a(5,4) = 1.7747158d0
       a(6,4) = 0.3420756d0
       a(1,5) = 0.4572609d0
       a(2,5) = 12.5661301d0
       a(3,5) = 3.5954593d0
       a(4,5) = 0.8968121d0
       a(5,5) = 0.1598810d0
       a(6,5) = 0.1447955d0
       n(1,1) = 0
       n(2,1) = 1
       n(3,1) = 2
       n(4,1) = 2
       n(5,1) = 2
       n(6,1) = 2
       n(1,2) = 0
       n(2,2) = 1
       n(3,2) = 2
       n(4,2) = 2
       n(5,2) = 2
       n(6,2) = 2
       n(1,3) = 0
       n(2,3) = 1
       n(3,3) = 2
       n(4,3) = 2
       n(5,3) = 2
       n(1,4) = 0
       n(2,4) = 1
       n(3,4) = 2
       n(4,4) = 2
       n(5,4) = 2
       n(6,4) = 2
       n(1,5) = 0
       n(2,5) = 1
       n(3,5) = 2
       n(4,5) = 2
       n(5,5) = 2
       n(6,5) = 2
      endif
      go to 200
c
c     ta hay
c
 90   lmx = 4
      if(olanl) then
       cc( 1, 1) =        -60.00000000d0
       cc( 2, 1) =       -629.71293130d0
       cc( 3, 1) =       -221.17505030d0
       cc( 4, 1) =        -83.47696830d0
       cc( 5, 1) =        -17.83702750d0
       cc( 6, 1) =         -0.81606090d0
       a ( 1, 1) =        768.22675330d0
       a ( 2, 1) =        173.59402180d0
       a ( 3, 1) =         43.70422680d0
       a ( 4, 1) =         11.61813190d0
       a ( 5, 1) =          4.03146180d0
       a ( 6, 1) =          1.09053780d0
       n ( 1, 1) =   1
       n ( 2, 1) =   2
       n ( 3, 1) =   2
       n ( 4, 1) =   2
       n ( 5, 1) =   2
       n ( 6, 1) =   2
       nt( 1)    =   6
       cc( 1, 2) =          3.00000000d0
       cc( 2, 2) =         33.97989850d0
       cc( 3, 2) =       1242.57762440d0
       cc( 4, 2) =        706.16638840d0
       cc( 5, 2) =        276.24645630d0
       cc( 6, 2) =        513.22831320d0
       cc( 7, 2) =       -412.96019160d0
       a ( 1, 2) =        365.41115590d0
       a ( 2, 2) =        697.13009150d0
       a ( 3, 2) =        271.78103790d0
       a ( 4, 2) =         87.30742630d0
       a ( 5, 2) =         22.74923390d0
       a ( 6, 2) =          3.86800820d0
       a ( 7, 2) =          3.67188910d0
       n ( 1, 2) =   0
       n ( 2, 2) =   1
       n ( 3, 2) =   2
       n ( 4, 2) =   2
       n ( 5, 2) =   2
       n ( 6, 2) =   2
       n ( 7, 2) =   2
       nt( 2)    =   7
       cc( 1, 3) =          2.00000000d0
       cc( 2, 3) =         61.81741630d0
       cc( 3, 3) =        194.49754790d0
       cc( 4, 3) =        307.83466710d0
       cc( 5, 3) =       -230.63463810d0
       a ( 1, 3) =        228.21073670d0
       a ( 2, 3) =         62.61854620d0
       a ( 3, 3) =         19.64656700d0
       a ( 4, 3) =          3.36515590d0
       a ( 5, 3) =          3.14972360d0
       n ( 1, 3) =   0
       n ( 2, 3) =   1
       n ( 3, 3) =   2
       n ( 4, 3) =   2
       n ( 5, 3) =   2
       nt( 3)    =   5
       cc( 1, 4) =          3.00000000d0
       cc( 2, 4) =         55.26709780d0
       cc( 3, 4) =        236.29575450d0
       cc( 4, 4) =        134.47550720d0
       cc( 5, 4) =         41.14011830d0
       a ( 1, 4) =        151.27555540d0
       a ( 2, 4) =         67.97949110d0
       a ( 3, 4) =         34.74747060d0
       a ( 4, 4) =         11.56564250d0
       a ( 5, 4) =          3.29301480d0
       n ( 1, 4) =   0
       n ( 2, 4) =   1
       n ( 3, 4) =   2
       n ( 4, 4) =   2
       n ( 5, 4) =   2
       nt( 4)    =   5
       cc( 1, 5) =          4.00000000d0
       cc( 2, 5) =         50.75337600d0
       cc( 3, 5) =        162.03128770d0
       cc( 4, 5) =         72.06192980d0
       cc( 5, 5) =          6.00098950d0
       a ( 1, 5) =         82.83680290d0
       a ( 2, 5) =         38.07198720d0
       a ( 3, 5) =         17.45470340d0
       a ( 4, 5) =          5.67859310d0
       a ( 5, 5) =          1.58140410d0
       n ( 1, 5) =   0
       n ( 2, 5) =   1
       n ( 3, 5) =   2
       n ( 4, 5) =   2
       n ( 5, 5) =   2
       nt( 5)    =   5
       nelcor    =  60
      else
       nelcor = 68
       nt(1) = 6
       nt(2) = 6
       nt(3) = 5
       nt(4) = 6
       nt(5) = 6
       cc(1,1) = -0.1418323d0
       cc(2,1) = -44.6105104d0
       cc(3,1) = -230.6825998d0
       cc(4,1) = -100.8703657d0
       cc(5,1) = -24.4967372d0
       cc(6,1) = -3.9862338d0
       cc(1,2) = 2.7695489d0
       cc(2,2) = 63.8012063d0
       cc(3,2) = 154.0916202d0
       cc(4,2) = 42.8001245d0
       cc(5,2) = 18.2668021d0
       cc(6,2) = -7.6317151d0
       cc(1,3) = 1.9022608d0
       cc(2,3) = 55.8292082d0
       cc(3,3) = 121.9975284d0
       cc(4,3) = 28.5395818d0
       cc(5,3) = 7.5403200d0
       cc(1,4) = 2.9426083d0
       cc(2,4) = 52.5314003d0
       cc(3,4) = 127.4437340d0
       cc(4,4) = 39.0662845d0
       cc(5,4) = -4.8954047d0
       cc(6,4) = -0.2876030d0
       cc(1,5) = 3.9592391d0
       cc(2,5) = 40.7108535d0
       cc(3,5) = 59.7356596d0
       cc(4,5) = 17.1992304d0
       cc(5,5) = 9.0725891d0
       cc(6,5) = -6.4941346d0
c
       a(1,1) = 438.9654138d0
       a(2,1) = 112.8894389d0
       a(3,1) = 39.0300509d0
       a(4,1) = 10.5663985d0
       a(5,1) = 3.1450008d0
       a(6,1) = 0.9962317d0
       a(1,2) = 88.7976994d0
       a(2,2) = 29.9605388d0
       a(3,2) = 10.1122075d0
       a(4,2) = 2.9034296d0
       a(5,2) = 0.6314816d0
       a(6,2) = 0.5171641d0
       a(1,3) = 82.8488496d0
       a(2,3) = 26.2097320d0
       a(3,3) = 8.2687456d0
       a(4,3) = 2.1411005d0
       a(5,3) = 0.5300236d0
       a(1,4) = 80.2483832d0
       a(2,4) = 28.8958219d0
       a(3,4) = 10.2668935d0
       a(4,4) = 2.7589509d0
       a(5,4) = 1.6942558d0
       a(6,4) = 0.3198755d0
       a(1,5) = 0.4782017d0
       a(2,5) = 16.4689620d0
       a(3,5) = 4.7111972d0
       a(4,5) = 1.0862608d0
       a(5,5) = 0.1969840d0
       a(6,5) = 0.1756469d0
       n(1,1) = 0
       n(2,1) = 1
       n(3,1) = 2
       n(4,1) = 2
       n(5,1) = 2
       n(6,1) = 2
       n(1,2) = 0
       n(2,2) = 1
       n(3,2) = 2
       n(4,2) = 2
       n(5,2) = 2
       n(6,2) = 2
       n(1,3) = 0
       n(2,3) = 1
       n(3,3) = 2
       n(4,3) = 2
       n(5,3) = 2
       n(1,4) = 0
       n(2,4) = 1
       n(3,4) = 2
       n(4,4) = 2
       n(5,4) = 2
       n(6,4) = 2
       n(1,5) = 0
       n(2,5) = 1
       n(3,5) = 2
       n(4,5) = 2
       n(5,5) = 2
       n(6,5) = 2
      endif
      go to 200
c
c     w hay
c
 100  lmx = 4
      if(olanl) then
       cc( 1, 1) =        -60.00000000d0
       cc( 2, 1) =       -664.19879200d0
       cc( 3, 1) =       -238.61436510d0
       cc( 4, 1) =        -88.41924070d0
       cc( 5, 1) =        -20.60623260d0
       cc( 6, 1) =         -0.92837920d0
       a ( 1, 1) =        839.44891200d0
       a ( 2, 1) =        192.85324820d0
       a ( 3, 1) =         48.66519740d0
       a ( 4, 1) =         12.92217270d0
       a ( 5, 1) =          4.57488900d0
       a ( 6, 1) =          1.26817960d0
       n ( 1, 1) =   1
       n ( 2, 1) =   2
       n ( 3, 1) =   2
       n ( 4, 1) =   2
       n ( 5, 1) =   2
       n ( 6, 1) =   2
       nt( 1)    =   6
       cc( 1, 2) =          3.00000000d0
       cc( 2, 2) =         39.11929670d0
       cc( 3, 2) =       1180.96929740d0
       cc( 4, 2) =        728.95642100d0
       cc( 5, 2) =        293.55911400d0
       cc( 6, 2) =        562.67314930d0
       cc( 7, 2) =       -457.38071850d0
       a ( 1, 2) =        313.42675180d0
       a ( 2, 2) =        699.31554620d0
       a ( 3, 2) =        259.89237410d0
       a ( 4, 2) =         85.49999800d0
       a ( 5, 2) =         22.76359250d0
       a ( 6, 2) =          4.07643170d0
       a ( 7, 2) =          3.88271620d0
       n ( 1, 2) =   0
       n ( 2, 2) =   1
       n ( 3, 2) =   2
       n ( 4, 2) =   2
       n ( 5, 2) =   2
       n ( 6, 2) =   2
       n ( 7, 2) =   2
       nt( 2)    =   7
       cc( 1, 3) =          2.00000000d0
       cc( 2, 3) =         63.89483930d0
       cc( 3, 3) =        205.89018370d0
       cc( 4, 3) =        312.14271530d0
       cc( 5, 3) =       -231.39612810d0
       a ( 1, 3) =        224.39264240d0
       a ( 2, 3) =         61.67369310d0
       a ( 3, 3) =         19.14690430d0
       a ( 4, 3) =          3.55657100d0
       a ( 5, 3) =          3.32632100d0
       n ( 1, 3) =   0
       n ( 2, 3) =   1
       n ( 3, 3) =   2
       n ( 4, 3) =   2
       n ( 5, 3) =   2
       nt( 3)    =   5
       cc( 1, 4) =          3.00000000d0
       cc( 2, 4) =         55.33152560d0
       cc( 3, 4) =        267.19766530d0
       cc( 4, 4) =        146.84855780d0
       cc( 5, 4) =         44.10552430d0
       a ( 1, 4) =        161.52789580d0
       a ( 2, 4) =         75.58146070d0
       a ( 3, 4) =         38.91158520d0
       a ( 4, 4) =         12.54262710d0
       a ( 5, 4) =          3.46151870d0
       n ( 1, 4) =   0
       n ( 2, 4) =   1
       n ( 3, 4) =   2
       n ( 4, 4) =   2
       n ( 5, 4) =   2
       nt( 4)    =   5
       cc( 1, 5) =          4.00000000d0
       cc( 2, 5) =         50.30655230d0
       cc( 3, 5) =        190.73630980d0
       cc( 4, 5) =         91.76055520d0
       cc( 5, 5) =          8.42473120d0
       a ( 1, 5) =         91.21027270d0
       a ( 2, 5) =         45.41527560d0
       a ( 3, 5) =         22.04529670d0
       a ( 4, 5) =          6.98104130d0
       a ( 5, 5) =          1.83804800d0
       n ( 1, 5) =   0
       n ( 2, 5) =   1
       n ( 3, 5) =   2
       n ( 4, 5) =   2
       n ( 5, 5) =   2
       nt( 5)    =   5
       nelcor    =  60
      else
       nelcor = 68
       nt(1) = 6
       nt(2) = 6
       nt(3) = 5
       nt(4) = 6
       nt(5) = 6
       cc(1,1) = -0.1457448d0
       cc(2,1) = -43.8187612d0
       cc(3,1) = -229.4536668d0
       cc(4,1) = -101.4253295d0
       cc(5,1) = -24.2972893d0
       cc(6,1) = -3.8771269d0
       cc(1,2) = 2.7625960d0
       cc(2,2) = 90.9432757d0
       cc(3,2) = 212.5686798d0
       cc(4,2) = 54.6162565d0
       cc(5,2) = 20.7664520d0
       cc(6,2) = -7.5448433d0
       cc(1,3) = 1.8994679d0
       cc(2,3) = 53.3500287d0
       cc(3,3) = 131.0179317d0
       cc(4,3) = 30.9930029d0
       cc(5,3) = 8.8195349d0
       cc(1,4) = 2.9410060d0
       cc(2,4) = 52.8348829d0
       cc(3,4) = 132.5222840d0
       cc(4,4) = 50.2360265d0
       cc(5,4) = -14.2411625d0
       cc(6,4) = -0.2994311d0
       cc(1,5) = 3.9581079d0
       cc(2,5) = 41.7245004d0
       cc(3,5) = 68.1862225d0
       cc(4,5) = 19.2847162d0
       cc(5,5) = 8.6028742d0
       cc(6,5) = -5.6939432d0
c
       a(1,1) = 329.3529656d0
       a(2,1) = 113.6102932d0
       a(3,1) = 39.9308141d0
       a(4,1) = 11.0046548d0
       a(5,1) = 3.3099176d0
       a(6,1) = 1.0597874d0
       a(1,2) = 120.7193610d0
       a(2,2) = 40.4863153d0
       a(3,2) = 12.2414413d0
       a(4,2) = 3.3536130d0
       a(5,2) = 0.7254557d0
       a(6,2) = 0.5842076d0
       a(1,3) = 97.7263572d0
       a(2,3) = 32.0030370d0
       a(3,3) = 10.3044815d0
       a(4,3) = 2.8403761d0
       a(5,3) = 0.5753263d0
       a(1,4) = 80.9725549d0
       a(2,4) = 29.8854425d0
       a(3,4) = 10.7716934d0
       a(4,4) = 2.7612542d0
       a(5,4) = 2.0639584d0
       a(6,4) = 0.3408971d0
       a(1,5) = 0.4778739d0
       a(2,5) = 19.1813197d0
       a(3,5) = 5.4823066d0
       a(4,5) = 1.2327431d0
       a(5,5) = 0.2348829d0
       a(6,5) = 0.2044616d0
       n(1,1) = 0
       n(2,1) = 1
       n(3,1) = 2
       n(4,1) = 2
       n(5,1) = 2
       n(6,1) = 2
       n(1,2) = 0
       n(2,2) = 1
       n(3,2) = 2
       n(4,2) = 2
       n(5,2) = 2
       n(6,2) = 2
       n(1,3) = 0
       n(2,3) = 1
       n(3,3) = 2
       n(4,3) = 2
       n(5,3) = 2
       n(1,4) = 0
       n(2,4) = 1
       n(3,4) = 2
       n(4,4) = 2
       n(5,4) = 2
       n(6,4) = 2
       n(1,5) = 0
       n(2,5) = 1
       n(3,5) = 2
       n(4,5) = 2
       n(5,5) = 2
       n(6,5) = 2
      endif
      go to 200
c
c     re hay
c
 110  lmx = 4
      if(olanl) then
       cc( 1, 1) =         -0.14971040d0
       cc( 2, 1) =      -1669.25557810d0
       cc( 3, 1) =       -346.66611290d0
       cc( 4, 1) =        -96.66848920d0
       cc( 5, 1) =        -11.07385670d0
       cc( 6, 1) =         -0.57985520d0
       a ( 1, 1) =       1833.81438360d0
       a ( 2, 1) =        414.05894560d0
       a ( 3, 1) =         59.03495400d0
       a ( 4, 1) =         11.91770940d0
       a ( 5, 1) =          3.65317640d0
       a ( 6, 1) =          1.27641840d0
       n ( 1, 1) =   1
       n ( 2, 1) =   2
       n ( 3, 1) =   2
       n ( 4, 1) =   2
       n ( 5, 1) =   2
       n ( 6, 1) =   2
       nt( 1)    =   6
       cc( 1, 2) =          3.14971040d0
       cc( 2, 2) =         45.29699690d0
       cc( 3, 2) =       1116.42795850d0
       cc( 4, 2) =        793.64190650d0
       cc( 5, 2) =        318.20991930d0
       cc( 6, 2) =        586.31050250d0
       cc( 7, 2) =       -475.69041640d0
       a ( 1, 2) =        247.00230850d0
       a ( 2, 2) =        631.71957670d0
       a ( 3, 2) =        249.28004410d0
       a ( 4, 2) =         84.29473150d0
       a ( 5, 2) =         23.04079120d0
       a ( 6, 2) =          4.29971110d0
       a ( 7, 2) =          4.09743900d0
       n ( 1, 2) =   0
       n ( 2, 2) =   1
       n ( 3, 2) =   2
       n ( 4, 2) =   2
       n ( 5, 2) =   2
       n ( 6, 2) =   2
       n ( 7, 2) =   2
       nt( 2)    =   7
       cc( 1, 3) =          2.14971040d0
       cc( 2, 3) =         64.22957360d0
       cc( 3, 3) =        218.85046070d0
       cc( 4, 3) =        344.41891680d0
       cc( 5, 3) =       -260.47272670d0
       a ( 1, 3) =        193.68977260d0
       a ( 2, 3) =         56.93083940d0
       a ( 3, 3) =         18.72276020d0
       a ( 4, 3) =          3.73400200d0
       a ( 5, 3) =          3.51274830d0
       n ( 1, 3) =   0
       n ( 2, 3) =   1
       n ( 3, 3) =   2
       n ( 4, 3) =   2
       n ( 5, 3) =   2
       nt( 3)    =   5
       cc( 1, 4) =          3.14971040d0
       cc( 2, 4) =         46.57148680d0
       cc( 3, 4) =        304.38345950d0
       cc( 4, 4) =        157.15896520d0
       cc( 5, 4) =         48.57012110d0
       a ( 1, 4) =        126.23408240d0
       a ( 2, 4) =         80.70908260d0
       a ( 3, 4) =         42.55297140d0
       a ( 4, 4) =         13.40744940d0
       a ( 5, 4) =          3.69042280d0
       n ( 1, 4) =   0
       n ( 2, 4) =   1
       n ( 3, 4) =   2
       n ( 4, 4) =   2
       n ( 5, 4) =   2
       nt( 4)    =   5
       cc( 1, 5) =          3.95696100d0
       cc( 2, 5) =         52.35867810d0
       cc( 3, 5) =        236.03408240d0
       cc( 4, 5) =        116.57872660d0
       cc( 5, 5) =         11.10505300d0
       a ( 1, 5) =        112.77864590d0
       a ( 2, 5) =         56.78159680d0
       a ( 3, 5) =         28.58319660d0
       a ( 4, 5) =          8.57341580d0
       a ( 5, 5) =          2.08927210d0
       n ( 1, 5) =   0
       n ( 2, 5) =   1
       n ( 3, 5) =   2
       n ( 4, 5) =   2
       n ( 5, 5) =   2
       nt( 5)    =   5
       nelcor    =  60
      else
       nelcor = 68
       nt(1) = 6
       nt(2) = 6
       nt(3) = 5
       nt(4) = 6
       nt(5) = 6
       cc(1,1) = -0.1497104d0
       cc(2,1) = -43.8574528d0
       cc(3,1) = -229.7901906d0
       cc(4,1) = -101.3640312d0
       cc(5,1) = -23.5051700d0
       cc(6,1) = -3.6178541d0
       cc(1,2) = 2.7555091d0
       cc(2,2) = 70.6591389d0
       cc(3,2) = 178.7034878d0
       cc(4,2) = 50.1066112d0
       cc(5,2) = 22.0145213d0
       cc(6,2) = -8.5755035d0
       cc(1,3) = 1.8966314d0
       cc(2,3) = 54.7084064d0
       cc(3,3) = 131.5647971d0
       cc(4,3) = 30.9572952d0
       cc(5,3) = 9.5032921d0
       cc(1,4) = 2.9393809d0
       cc(2,4) = 53.4995776d0
       cc(3,4) = 138.2787174d0
       cc(4,4) = 61.4189098d0
       cc(5,4) = -23.6282276d0
       cc(6,4) = -0.2989651d0
       cc(1,5) = 3.9569610d0
       cc(2,5) = 42.3484376d0
       cc(3,5) = 69.0825777d0
       cc(4,5) = 20.4890160d0
       cc(5,5) = 6.5951586d0
       cc(6,5) = -3.6998106d0
c
       a(1,1) = 341.4416283d0
       a(2,1) = 115.4634850d0
       a(3,1) = 40.5580378d0
       a(4,1) = 11.2759241d0
       a(5,1) = 3.3967771d0
       a(6,1) = 1.1070830d0
       a(1,2) = 101.0857749d0
       a(2,2) = 34.2983377d0
       a(3,2) = 11.5707921d0
       a(4,2) = 3.3540317d0
       a(5,2) = 0.7432483d0
       a(6,2) = 0.6011833d0
       a(1,3) = 94.6617320d0
       a(2,3) = 30.8125403d0
       a(3,3) = 9.9084101d0
       a(4,3) = 2.6363334d0
       a(5,3) = 0.6257650d0
       a(1,4) = 83.9171495d0
       a(2,4) = 31.2351625d0
       a(3,4) = 11.3300093d0
       a(4,4) = 2.8056816d0
       a(5,4) = 2.2303172d0
       a(6,4) = 0.3547500d0
       a(1,5) = 0.4765942d0
       a(2,5) = 19.4095736d0
       a(3,5) = 5.5717284d0
       a(4,5) = 1.2758944d0
       a(5,5) = 0.2493947d0
       a(6,5) = 0.2054839d0
       n(1,1) = 0
       n(2,1) = 1
       n(3,1) = 2
       n(4,1) = 2
       n(5,1) = 2
       n(6,1) = 2
       n(1,2) = 0
       n(2,2) = 1
       n(3,2) = 2
       n(4,2) = 2
       n(5,2) = 2
       n(6,2) = 2
       n(1,3) = 0
       n(2,3) = 1
       n(3,3) = 2
       n(4,3) = 2
       n(5,3) = 2
       n(1,4) = 0
       n(2,4) = 1
       n(3,4) = 2
       n(4,4) = 2
       n(5,4) = 2
       n(6,4) = 2
       n(1,5) = 0
       n(2,5) = 1
       n(3,5) = 2
       n(4,5) = 2
       n(5,5) = 2
       n(6,5) = 2
      endif
      go to 200
c
c     os hay
c
 120  lmx = 4
      if(olanl) then
       cc( 1, 1) =         -0.15372930d0
       cc( 2, 1) =      -1642.96955420d0
       cc( 3, 1) =       -336.79545670d0
       cc( 4, 1) =        -94.98390270d0
       cc( 5, 1) =        -10.35162180d0
       a ( 1, 1) =        848.87432920d0
       a ( 2, 1) =        400.46400040d0
       a ( 3, 1) =         58.03594660d0
       a ( 4, 1) =         12.04546820d0
       a ( 5, 1) =          3.44425120d0
       n ( 1, 1) =   1
       n ( 2, 1) =   2
       n ( 3, 1) =   2
       n ( 4, 1) =   2
       n ( 5, 1) =   2
       nt( 1)    =   5
       cc( 1, 2) =          3.15372930d0
       cc( 2, 2) =         25.99482330d0
       cc( 3, 2) =        859.48657560d0
       cc( 4, 2) =        380.13765040d0
       cc( 5, 2) =        242.70106660d0
       cc( 6, 2) =       -125.17531070d0
       a ( 1, 2) =        216.26722490d0
       a ( 2, 2) =        343.85820980d0
       a ( 3, 2) =        136.05236060d0
       a ( 4, 2) =         35.03894570d0
       a ( 5, 2) =          4.47244520d0
       a ( 6, 2) =          3.76751190d0
       n ( 1, 2) =   0
       n ( 2, 2) =   1
       n ( 3, 2) =   2
       n ( 4, 2) =   2
       n ( 5, 2) =   2
       n ( 6, 2) =   2
       nt( 2)    =   6
       cc( 1, 3) =          2.15372930d0
       cc( 2, 3) =         62.86842750d0
       cc( 3, 3) =        293.13229870d0
       cc( 4, 3) =        264.46483880d0
       cc( 5, 3) =       -172.64203140d0
       a ( 1, 3) =        326.06001340d0
       a ( 2, 3) =         95.63302060d0
       a ( 3, 3) =         32.52610430d0
       a ( 4, 3) =          3.76997300d0
       a ( 5, 3) =          3.38958590d0
       n ( 1, 3) =   0
       n ( 2, 3) =   1
       n ( 3, 3) =   2
       n ( 4, 3) =   2
       n ( 5, 3) =   2
       nt( 3)    =   5
       cc( 1, 4) =          3.15372930d0
       cc( 2, 4) =         45.33992000d0
       cc( 3, 4) =        335.85750320d0
       cc( 4, 4) =        166.30374570d0
       cc( 5, 4) =         49.02378250d0
       a ( 1, 4) =        131.03199290d0
       a ( 2, 4) =         89.24649030d0
       a ( 3, 4) =         46.39769800d0
       a ( 4, 4) =         14.18660530d0
       a ( 5, 4) =          3.76688430d0
       n ( 1, 4) =   0
       n ( 2, 4) =   1
       n ( 3, 4) =   2
       n ( 4, 4) =   2
       n ( 5, 4) =   2
       nt( 4)    =   5
       cc( 1, 5) =          3.95579820d0
       cc( 2, 5) =         51.95219570d0
       cc( 3, 5) =        314.25558660d0
       cc( 4, 5) =        147.18398940d0
       cc( 5, 5) =         13.96229470d0
       a ( 1, 5) =        134.73681480d0
       a ( 2, 5) =         73.87713180d0
       a ( 3, 5) =         38.10304250d0
       a ( 4, 5) =         10.53376950d0
       a ( 5, 5) =          2.33382460d0
       n ( 1, 5) =   0
       n ( 2, 5) =   1
       n ( 3, 5) =   2
       n ( 4, 5) =   2
       n ( 5, 5) =   2
       nt( 5)    =   5
       nelcor    =  60
      else
       nelcor = 68
       nt(1) = 6
       nt(2) = 6
       nt(3) = 5
       nt(4) = 6
       nt(5) = 6
       cc(1,1) = -0.1537293d0
       cc(2,1) = -43.9385285d0
       cc(3,1) = -231.3805876d0
       cc(4,1) = -101.9942766d0
       cc(5,1) = -23.0798574d0
       cc(6,1) = -3.4299171d0
       cc(1,2) = 2.7482860d0
       cc(2,2) = 66.2898847d0
       cc(3,2) = 173.0951085d0
       cc(4,2) = 49.9172313d0
       cc(5,2) = 21.1923935d0
       cc(6,2) = -7.2489733d0
       cc(1,3) = 1.8937508d0
       cc(2,3) = 58.1881981d0
       cc(3,3) = 139.0439799d0
       cc(4,3) = 33.8246611d0
       cc(5,3) = 10.1980590d0
       cc(1,4) = 2.9377328d0
       cc(2,4) = 53.3061878d0
       cc(3,4) = 146.6642182d0
       cc(4,4) = 73.2594224d0
       cc(5,4) = -34.2036412d0
       cc(6,4) = -0.2733322d0
       cc(1,5) = 3.9557982d0
       cc(2,5) = 43.4550088d0
       cc(3,5) = 78.7929165d0
       cc(4,5) = 16.6983488d0
       cc(5,5) = 6.9521119d0
       cc(6,5) = -3.8580472d0
c
       a(1,1) = 383.4894828d0
       a(2,1) = 118.1733466d0
       a(3,1) = 41.5391382d0
       a(4,1) = 11.6367092d0
       a(5,1) = 3.5088244d0
       a(6,1) = 1.1562003d0
       a(1,2) = 101.5361305d0
       a(2,2) = 34.3114388d0
       a(3,2) = 11.8625113d0
       a(4,2) = 3.4755205d0
       a(5,2) = 0.7821672d0
       a(6,2) = 0.6099043d0
       a(1,3) = 95.6090848d0
       a(2,3) = 30.5900495d0
       a(3,3) = 9.8352272d0
       a(4,3) = 2.5958260d0
       a(5,3) = 0.6766517d0
       a(1,4) = 86.9478876d0
       a(2,4) = 33.4368860d0
       a(3,4) = 12.2083862d0
       a(4,4) = 2.8023849d0
       a(5,4) = 2.2970520d0
       a(6,4) = 0.3556000d0
       a(1,5) = 1.3712305d0
       a(2,5) = 23.9889426d0
       a(3,5) = 6.8090456d0
       a(4,5) = 1.3516761d0
       a(5,5) = 0.2925597d0
       a(6,5) = 0.2375690d0
       n(1,1) = 0
       n(2,1) = 1
       n(3,1) = 2
       n(4,1) = 2
       n(5,1) = 2
       n(6,1) = 2
       n(1,2) = 0
       n(2,2) = 1
       n(3,2) = 2
       n(4,2) = 2
       n(5,2) = 2
       n(6,2) = 2
       n(1,3) = 0
       n(2,3) = 1
       n(3,3) = 2
       n(4,3) = 2
       n(5,3) = 2
       n(1,4) = 0
       n(2,4) = 1
       n(3,4) = 2
       n(4,4) = 2
       n(5,4) = 2
       n(6,4) = 2
       n(1,5) = 0
       n(2,5) = 1
       n(3,5) = 2
       n(4,5) = 2
       n(5,5) = 2
       n(6,5) = 2
      endif
      go to 200
c
c     ir hay
c
 130  lmx = 4
      if(olanl) then
       cc( 1, 1) =         -0.15780140d0
       cc( 2, 1) =      -1517.52704460d0
       cc( 3, 1) =       -316.53065290d0
       cc( 4, 1) =        -91.88809410d0
       cc( 5, 1) =         -9.22417730d0
       a ( 1, 1) =        823.58801470d0
       a ( 2, 1) =        364.66133360d0
       a ( 3, 1) =         55.70828010d0
       a ( 4, 1) =         12.04645440d0
       a ( 5, 1) =          3.51206100d0
       n ( 1, 1) =   1
       n ( 2, 1) =   2
       n ( 3, 1) =   2
       n ( 4, 1) =   2
       n ( 5, 1) =   2
       nt( 1)    =   5
       cc( 1, 2) =          3.15780140d0
       cc( 2, 2) =         26.83225770d0
       cc( 3, 2) =        800.42500070d0
       cc( 4, 2) =        369.40506830d0
       cc( 5, 2) =        242.41718990d0
       cc( 6, 2) =       -118.21732820d0
       a ( 1, 2) =        188.04907700d0
       a ( 2, 2) =        340.41947120d0
       a ( 3, 2) =        128.23736730d0
       a ( 4, 2) =         33.86449610d0
       a ( 5, 2) =          4.75600050d0
       a ( 6, 2) =          3.96499740d0
       n ( 1, 2) =   0
       n ( 2, 2) =   1
       n ( 3, 2) =   2
       n ( 4, 2) =   2
       n ( 5, 2) =   2
       n ( 6, 2) =   2
       nt( 2)    =   6
       cc( 1, 3) =          2.15780140d0
       cc( 2, 3) =         61.96786100d0
       cc( 3, 3) =        269.05819860d0
       cc( 4, 3) =        231.16547930d0
       cc( 5, 3) =       -133.69526670d0
       a ( 1, 3) =        289.72911390d0
       a ( 2, 3) =         87.46337890d0
       a ( 3, 3) =         30.43637660d0
       a ( 4, 3) =          4.05534120d0
       a ( 5, 3) =          3.55253410d0
       n ( 1, 3) =   0
       n ( 2, 3) =   1
       n ( 3, 3) =   2
       n ( 4, 3) =   2
       n ( 5, 3) =   2
       nt( 3)    =   5
       cc( 1, 4) =          3.15780140d0
       cc( 2, 4) =         45.93498030d0
       cc( 3, 4) =        359.03446680d0
       cc( 4, 4) =        176.47401190d0
       cc( 5, 4) =         54.51552860d0
       a ( 1, 4) =        136.40171060d0
       a ( 2, 4) =         95.07769250d0
       a ( 3, 4) =         49.22584100d0
       a ( 4, 4) =         15.08741450d0
       a ( 5, 4) =          4.04057640d0
       n ( 1, 4) =   0
       n ( 2, 4) =   1
       n ( 3, 4) =   2
       n ( 4, 4) =   2
       n ( 5, 4) =   2
       nt( 4)    =   5
       cc( 1, 5) =          3.95461970d0
       cc( 2, 5) =         52.97736550d0
       cc( 3, 5) =        274.86433830d0
       cc( 4, 5) =        137.20473380d0
       cc( 5, 5) =         14.86333050d0
       a ( 1, 5) =        127.35079080d0
       a ( 2, 5) =         66.23643740d0
       a ( 3, 5) =         34.42992290d0
       a ( 4, 5) =         10.19957210d0
       a ( 5, 5) =          2.54097020d0
       n ( 1, 5) =   0
       n ( 2, 5) =   1
       n ( 3, 5) =   2
       n ( 4, 5) =   2
       n ( 5, 5) =   2
       nt( 5)    =   5
       nelcor    =  60
      else
       nelcor = 68
       nt(1) = 6
       nt(2) = 6
       nt(3) = 5
       nt(4) = 6
       nt(5) = 6
       cc(1,1) = -0.1578014d0
       cc(2,1) = -44.3754600d0
       cc(3,1) = -235.2707019d0
       cc(4,1) = -102.4292503d0
       cc(5,1) = -22.4183750d0
       cc(6,1) = -3.1794489d0
       cc(1,2) = 2.7409244d0
       cc(2,2) = 92.2686824d0
       cc(3,2) = 237.8253940d0
       cc(4,2) = 64.6371640d0
       cc(5,2) = 22.0167490d0
       cc(6,2) = -4.3711058d0
       cc(1,3) = 1.8908260d0
       cc(2,3) = 55.9632847d0
       cc(3,3) = 150.5666843d0
       cc(4,3) = 37.6873512d0
       cc(5,3) = 11.5570752d0
       cc(1,4) = 2.9360618d0
       cc(2,4) = 54.5844845d0
       cc(3,4) = 153.9981694d0
       cc(4,4) = 75.9710731d0
       cc(5,4) = -32.5058084d0
       cc(6,4) = -0.3246828d0
       cc(1,5) = 3.9546197d0
       cc(2,5) = 43.4409800d0
       cc(3,5) = 79.6338840d0
       cc(4,5) = 29.6292010d0
       cc(5,5) = 5.7419073d0
       cc(6,5) = -1.9369124d0
c
       a(1,1) = 519.6134088d0
       a(2,1) = 122.3129747d0
       a(3,1) = 42.5567718d0
       a(4,1) = 11.9253010d0
       a(5,1) = 3.5816124d0
       a(6,1) = 1.1941008d0
       a(1,2) = 132.6007649d0
       a(2,2) = 45.7455978d0
       a(3,2) = 14.3850395d0
       a(4,2) = 4.0679871d0
       a(5,2) = 0.9498427d0
       a(6,2) = 0.6934006d0
       a(1,3) = 114.1855351d0
       a(2,3) = 37.6106966d0
       a(3,3) = 12.3532050d0
       a(4,3) = 3.4744038d0
       a(5,3) = 0.7246351d0
       a(1,4) = 90.6547115d0
       a(2,4) = 34.9840739d0
       a(3,4) = 12.9003051d0
       a(4,4) = 3.0669900d0
       a(5,4) = 2.4884834d0
       a(6,4) = 0.4041198d0
       a(1,5) = 5.7096181d0
       a(2,5) = 19.3598869d0
       a(3,5) = 5.4073364d0
       a(4,5) = 1.4199283d0
       a(5,5) = 0.3315322d0
       a(6,5) = 0.2340301d0
       n(1,1) = 0
       n(2,1) = 1
       n(3,1) = 2
       n(4,1) = 2
       n(5,1) = 2
       n(6,1) = 2
       n(1,2) = 0
       n(2,2) = 1
       n(3,2) = 2
       n(4,2) = 2
       n(5,2) = 2
       n(6,2) = 2
       n(1,3) = 0
       n(2,3) = 1
       n(3,3) = 2
       n(4,3) = 2
       n(5,3) = 2
       n(1,4) = 0
       n(2,4) = 1
       n(3,4) = 2
       n(4,4) = 2
       n(5,4) = 2
       n(6,4) = 2
       n(1,5) = 0
       n(2,5) = 1
       n(3,5) = 2
       n(4,5) = 2
       n(5,5) = 2
       n(6,5) = 2
      endif
      go to 200
c
c     pt hay
c
 140  lmx = 4
      if(olanl) then
       cc( 1, 1) =         -0.16192680d0
       cc( 2, 1) =      -1320.28738520d0
       cc( 3, 1) =       -298.31781350d0
       cc( 4, 1) =        -87.58370650d0
       cc( 5, 1) =         -8.14932740d0
       a ( 1, 1) =        728.93940560d0
       a ( 2, 1) =        320.65678000d0
       a ( 3, 1) =         52.86801740d0
       a ( 4, 1) =         12.02801280d0
       a ( 5, 1) =          3.52389130d0
       n ( 1, 1) =   1
       n ( 2, 1) =   2
       n ( 3, 1) =   2
       n ( 4, 1) =   2
       n ( 5, 1) =   2
       nt( 1)    =   5
       cc( 1, 2) =          2.73342180d0
       cc( 2, 2) =         59.70243290d0
       cc( 3, 2) =        891.45895500d0
       cc( 4, 2) =        368.44676560d0
       cc( 5, 2) =        238.02630900d0
       cc( 6, 2) =       -107.05564540d0
       a ( 1, 2) =        409.44373580d0
       a ( 2, 2) =        274.54192310d0
       a ( 3, 2) =        127.56585700d0
       a ( 4, 2) =         32.90366310d0
       a ( 5, 2) =          5.05938800d0
       a ( 6, 2) =          4.15065560d0
       n ( 1, 2) =   0
       n ( 2, 2) =   1
       n ( 3, 2) =   2
       n ( 4, 2) =   2
       n ( 5, 2) =   2
       n ( 6, 2) =   2
       nt( 2)    =   6
       cc( 1, 3) =          1.88785680d0
       cc( 2, 3) =         76.01386290d0
       cc( 3, 3) =        343.55111160d0
       cc( 4, 3) =        119.49117860d0
       a ( 1, 3) =        466.17288920d0
       a ( 2, 3) =        120.78882590d0
       a ( 3, 3) =         36.41187910d0
       a ( 4, 3) =          5.69854080d0
       n ( 1, 3) =   0
       n ( 2, 3) =   1
       n ( 3, 3) =   2
       n ( 4, 3) =   2
       nt( 3)    =   4
       cc( 1, 4) =          2.93436780d0
       cc( 2, 4) =         59.33065710d0
       cc( 3, 4) =        452.44451940d0
       cc( 4, 4) =        210.47694790d0
       cc( 5, 4) =         58.62541120d0
       a ( 1, 4) =        249.56507630d0
       a ( 2, 4) =        126.66785850d0
       a ( 3, 4) =         63.14305860d0
       a ( 4, 4) =         17.90594700d0
       a ( 5, 4) =          4.22393730d0
       n ( 1, 4) =   0
       n ( 2, 4) =   1
       n ( 3, 4) =   2
       n ( 4, 4) =   2
       n ( 5, 4) =   2
       nt( 4)    =   5
       cc( 1, 5) =          3.95342530d0
       cc( 2, 5) =         53.85551820d0
       cc( 3, 5) =        247.43051330d0
       cc( 4, 5) =        127.81879760d0
       cc( 5, 5) =         15.37720460d0
       a ( 1, 5) =        121.81587990d0
       a ( 2, 5) =         60.87570300d0
       a ( 3, 5) =         31.47671470d0
       a ( 4, 5) =          9.88117510d0
       a ( 5, 5) =          2.73198740d0
       n ( 1, 5) =   0
       n ( 2, 5) =   1
       n ( 3, 5) =   2
       n ( 4, 5) =   2
       n ( 5, 5) =   2
       nt( 5)    =   5
       nelcor    =  60
      else
       nelcor = 68
       nt(1) = 6
       nt(2) = 6
       nt(3) = 5
       nt(4) = 6
       nt(5) = 6
       cc(1,1) = -0.1619268d0
       cc(2,1) = -44.1056601d0
       cc(3,1) = -225.8208289d0
       cc(4,1) = -99.3462061d0
       cc(5,1) = -20.9272244d0
       cc(6,1) = -2.7509383d0
       cc(1,2) = 2.7334218d0
       cc(2,2) = 76.8788294d0
       cc(3,2) = 208.3220582d0
       cc(4,2) = 60.0312229d0
       cc(5,2) = 25.9812187d0
       cc(6,2) = -8.6418215d0
       cc(1,3) = 1.8878568d0
       cc(2,3) = 57.1033025d0
       cc(3,3) = 150.6780890d0
       cc(4,3) = 37.2985110d0
       cc(5,3) = 12.3081605d0
       cc(1,4) = 2.9343678d0
       cc(2,4) = 52.8691101d0
       cc(3,4) = 170.6411111d0
       cc(4,4) = 99.8353949d0
       cc(5,4) = -54.3170702d0
       cc(6,4) = -0.2872878d0
       cc(1,5) = 3.9534253d0
       cc(2,5) = 44.9889323d0
       cc(3,5) = 85.3111945d0
       cc(4,5) = 28.7619037d0
       cc(5,5) = 6.8855706d0
       cc(6,5) = -2.6742821d0
c
       a(1,1) = 415.2990463d0
       a(2,1) = 117.9901046d0
       a(3,1) = 41.2716533d0
       a(4,1) = 11.8572771d0
       a(5,1) = 3.5452121d0
       a(6,1) = 1.2031453d0
       a(1,2) = 118.3123634d0
       a(2,2) = 40.2269849d0
       a(3,2) = 13.6420645d0
       a(4,2) = 3.9884442d0
       a(5,2) = 0.9246716d0
       a(6,2) = 0.7260083d0
       a(1,3) = 111.7931938d0
       a(2,3) = 36.4526949d0
       a(3,3) = 11.9271531d0
       a(4,3) = 3.2553879d0
       a(5,3) = 0.7783811d0
       a(1,4) = 92.3112664d0
       a(2,4) = 39.3772891d0
       a(3,4) = 14.6918292d0
       a(4,4) = 3.0234056d0
       a(5,4) = 2.5666981d0
       a(6,4) = 0.4019102d0
       a(1,5) = 1.6524570d0
       a(2,5) = 23.2155161d0
       a(3,5) = 6.7336889d0
       a(4,5) = 1.4985371d0
       a(5,5) = 0.3567339d0
       a(6,5) = 0.2699638d0
       n(1,1) = 0
       n(2,1) = 1
       n(3,1) = 2
       n(4,1) = 2
       n(5,1) = 2
       n(6,1) = 2
       n(1,2) = 0
       n(2,2) = 1
       n(3,2) = 2
       n(4,2) = 2
       n(5,2) = 2
       n(6,2) = 2
       n(1,3) = 0
       n(2,3) = 1
       n(3,3) = 2
       n(4,3) = 2
       n(5,3) = 2
       n(1,4) = 0
       n(2,4) = 1
       n(3,4) = 2
       n(4,4) = 2
       n(5,4) = 2
       n(6,4) = 2
       n(1,5) = 0
       n(2,5) = 1
       n(3,5) = 2
       n(4,5) = 2
       n(5,5) = 2
       n(6,5) = 2
      endif
      go to 200
c
c     au   hay
c
 150  lmx = 4
      if(olanl) then
       cc( 1, 1) =        -60.00000000d0
       cc( 2, 1) =       -555.52923120d0
       cc( 3, 1) =       -168.00197850d0
       cc( 4, 1) =        -63.03998750d0
       cc( 5, 1) =         -4.25166810d0
       a ( 1, 1) =        622.62879560d0
       a ( 2, 1) =        136.28436070d0
       a ( 3, 1) =         33.15497810d0
       a ( 4, 1) =          9.98948950d0
       a ( 5, 1) =          3.04813120d0
       n ( 1, 1) =   1
       n ( 2, 1) =   2
       n ( 3, 1) =   2
       n ( 4, 1) =   2
       n ( 5, 1) =   2
       nt( 1)    =   5
       cc( 1, 2) =          3.00000000d0
       cc( 2, 2) =         38.60208800d0
       cc( 3, 2) =        864.83707270d0
       cc( 4, 2) =        374.99355200d0
       cc( 5, 2) =        289.79101000d0
       cc( 6, 2) =       -152.45327730d0
       a ( 1, 2) =        194.73743040d0
       a ( 2, 2) =        351.53274470d0
       a ( 3, 2) =        122.32704020d0
       a ( 4, 2) =         32.09146170d0
       a ( 5, 2) =          5.24518120d0
       a ( 6, 2) =          4.49162230d0
       n ( 1, 2) =   0
       n ( 2, 2) =   1
       n ( 3, 2) =   2
       n ( 4, 2) =   2
       n ( 5, 2) =   2
       n ( 6, 2) =   2
       nt( 2)    =   6
       cc( 1, 3) =          2.00000000d0
       cc( 2, 3) =         73.88856250d0
       cc( 3, 3) =        326.67298720d0
       cc( 4, 3) =        126.58145910d0
       a ( 1, 3) =        420.61588010d0
       a ( 2, 3) =        109.44178150d0
       a ( 3, 3) =         34.17142800d0
       a ( 4, 3) =          5.98797500d0
       n ( 1, 3) =   0
       n ( 2, 3) =   1
       n ( 3, 3) =   2
       n ( 4, 3) =   2
       nt( 3)    =   4
       cc( 1, 4) =          3.00000000d0
       cc( 2, 4) =         55.67931490d0
       cc( 3, 4) =        449.19873350d0
       cc( 4, 4) =        215.02690910d0
       cc( 5, 4) =         64.08409950d0
       a ( 1, 4) =        219.26661580d0
       a ( 2, 4) =        122.72977860d0
       a ( 3, 4) =         63.10633690d0
       a ( 4, 4) =         18.36845200d0
       a ( 5, 4) =          4.49728440d0
       n ( 1, 4) =   0
       n ( 2, 4) =   1
       n ( 3, 4) =   2
       n ( 4, 4) =   2
       n ( 5, 4) =   2
       nt( 4)    =   5
       cc( 1, 5) =          4.00000000d0
       cc( 2, 5) =         51.80653350d0
       cc( 3, 5) =        231.21831130d0
       cc( 4, 5) =        119.00473860d0
       cc( 5, 5) =         15.34241880d0
       a ( 1, 5) =        108.55060370d0
       a ( 2, 5) =         56.47955270d0
       a ( 3, 5) =         29.20691590d0
       a ( 4, 5) =          9.54405430d0
       a ( 5, 5) =          2.89651180d0
       n ( 1, 5) =   0
       n ( 2, 5) =   1
       n ( 3, 5) =   2
       n ( 4, 5) =   2
       n ( 5, 5) =   2
       nt( 5)    =   5
       nelcor    =  60
      else
       nelcor = 68
       nt(1) = 6
       nt(2) = 6
       nt(3) = 5
       nt(4) = 6
       nt(5) = 6
       cc(1,1) = -0.1661054d0
       cc(2,1) = -43.7131382d0
       cc(3,1) = -212.4673199d0
       cc(4,1) = -94.7269807d0
       cc(5,1) = -19.4707222d0
       cc(6,1) = -2.3733091d0
       cc(1,2) = 2.7257756d0
       cc(2,2) = 71.6437344d0
       cc(3,2) = 201.6713462d0
       cc(4,2) = 60.8104382d0
       cc(5,2) = 34.2897810d0
       cc(6,2) = -15.7285421d0
       cc(1,3) = 1.8848427d0
       cc(2,3) = 58.2493046d0
       cc(3,3) = 153.2713054d0
       cc(4,3) = 38.7047152d0
       cc(5,3) = 12.8839047d0
       cc(1,4) = 2.9326507d0
       cc(2,4) = 53.7261877d0
       cc(3,4) = 173.1778236d0
       cc(4,4) = 136.1235752d0
       cc(5,4) = -86.8401404d0
       cc(6,4) = -0.2963681d0
       cc(1,5) = 3.9522151d0
       cc(2,5) = 44.9423254d0
       cc(3,5) = 75.3617449d0
       cc(4,5) = 27.4263473d0
       cc(5,5) = 7.7511742d0
       cc(6,5) = -4.2609308d0
c
       a(1,1) = 433.5024482d0
       a(2,1) = 111.3294046d0
       a(3,1) = 39.2392980d0
       a(4,1) = 11.7036609d0
       a(5,1) = 3.5081904d0
       a(6,1) = 1.2088436d0
       a(1,2) = 117.1237904d0
       a(2,2) = 39.8518385d0
       a(3,2) = 14.0221281d0
       a(4,2) = 4.1944438d0
       a(5,2) = 0.9742593d0
       a(6,2) = 0.8357820d0
       a(1,3) = 109.8661084d0
       a(2,3) = 35.8645320d0
       a(3,3) = 11.8152657d0
       a(4,3) = 3.1776488d0
       a(5,3) = 0.8280106d0
       a(1,4) = 92.3129639d0
       a(2,4) = 39.4155045d0
       a(3,4) = 14.9349958d0
       a(4,4) = 3.1809761d0
       a(5,4) = 2.8278414d0
       a(6,4) = 0.4201538d0
       a(1,5) = 1.4516634d0
       a(2,5) = 20.7062622d0
       a(3,5) = 6.0590141d0
       a(4,5) = 1.4432741d0
       a(5,5) = 0.3276190d0
       a(6,5) = 0.2715330d0
       n(1,1) = 0
       n(2,1) = 1
       n(3,1) = 2
       n(4,1) = 2
       n(5,1) = 2
       n(6,1) = 2
       n(1,2) = 0
       n(2,2) = 1
       n(3,2) = 2
       n(4,2) = 2
       n(5,2) = 2
       n(6,2) = 2
       n(1,3) = 0
       n(2,3) = 1
       n(3,3) = 2
       n(4,3) = 2
       n(5,3) = 2
       n(1,4) = 0
       n(2,4) = 1
       n(3,4) = 2
       n(4,4) = 2
       n(5,4) = 2
       n(6,4) = 2
       n(1,5) = 0
       n(2,5) = 1
       n(3,5) = 2
       n(4,5) = 2
       n(5,5) = 2
       n(6,5) = 2
      endif
      go to 200
c
c     hg   hay
c
 160  if(olanl) then
       go to 300
      else
      lmx = 4
      nelcor = 68
      nt(1) = 6
      nt(2) = 5
      nt(3) = 5
      nt(4) = 6
      nt(5) = 6
      cc(1,1) = -0.1703372d0
      cc(2,1) = -24.5189676d0
      cc(3,1) = -28.2305328d0
      cc(4,1) = -50.2160217d0
      cc(5,1) = -13.8325534d0
      cc(6,1) = -1.3977261d0
      cc(1,2) = 2.7440421d0
      cc(2,2) = 68.8671343d0
      cc(3,2) = 246.8006073d0
      cc(4,2) = 97.7864206d0
      cc(5,2) = 26.5465894d0
      cc(1,3) = 1.9446525d0
      cc(2,3) = 64.9315822d0
      cc(3,3) = 272.3521790d0
      cc(4,3) = 138.5794426d0
      cc(5,3) = 41.2158199d0
      cc(6,3) = 15.1712756d0
      cc(1,4) = 3.1153493d0
      cc(2,4) = 44.5580217d0
      cc(3,4) = 230.2588307d0
      cc(4,4) = 98.5271156d0
      cc(5,4) = -43.3718937d0
      cc(6,4) = -0.3260880d0
      cc(1,5) = 4.0027249d0
      cc(2,5) = 53.9600220d0
      cc(3,5) = 25.2756026d0
      cc(4,5) = 145.9184228d0
      cc(5,5) = 11.7637636d0
      cc(6,5) = -0.6759736d0
c
      a(1,1) = 318.3460385d0
      a(2,1) = 95.6187191d0
      a(3,1) = 17.5652631d0
      a(4,1) = 8.6336796d0
      a(5,1) = 3.0025515d0
      a(6,1) = 1.1004827d0
      a(1,2) = 180.0082018d0
      a(2,2) = 60.9782063d0
      a(3,2) = 22.9444668d0
      a(4,2) = 6.8681457d0
      a(5,2) = 1.3345231d0
      a(1,3) = 260.5481365d0
      a(2,3) = 89.5728959d0
      a(3,3) = 37.1461441d0
      a(4,3) = 12.4960359d0
      a(5,3) = 3.7833475d0
      a(6,3) = 0.9041993d0
      a(1,4) = 72.0586232d0
      a(2,4) = 53.8768738d0
      a(3,4) = 19.6556989d0
      a(4,4) = 3.4386734d0
      a(5,4) = 2.7552072d0
      a(6,4) = 0.4923262d0
      a(1,5) = 86.9008481d0
      a(2,5) = 31.3844809d0
      a(3,5) = 33.5683319d0
      a(4,5) = 10.6046957d0
      a(5,5) = 2.2121230d0
      a(6,5) = 0.7670500d0
      n(1,1) = 0
      n(2,1) = 1
      n(3,1) = 1
      n(4,1) = 2
      n(5,1) = 2
      n(6,1) = 2
      n(1,2) = 0
      n(2,2) = 1
      n(3,2) = 2
      n(4,2) = 2
      n(5,2) = 2
      n(1,3) = 0
      n(2,3) = 1
      n(3,3) = 2
      n(4,3) = 2
      n(5,3) = 2
      n(6,3) = 2
      n(1,4) = 0
      n(2,4) = 1
      n(3,4) = 2
      n(4,4) = 2
      n(5,4) = 2
      n(6,4) = 2
      n(1,5) = 0
      n(2,5) = 1
      n(3,5) = 2
      n(4,5) = 2
      n(5,5) = 2
      n(6,5) = 2
      endif
      go to 200
c     tl  hay (13-electron values)
 170  if(olanl) then
       go to 300
      else
      lmx = 4
      nelcor = 68
      nt(1) = 6
      nt(2) = 5
      nt(3) = 6
      nt(4) = 6
      nt(5) = 6
      cc(1,1) = -0.1746222d0
      cc(2,1) = -24.5858057d0
      cc(3,1) = -28.5214719d0
      cc(4,1) = -47.9977348d0
      cc(5,1) = -13.1733337d0
      cc(6,1) = -1.1972529d0
      cc(1,2) = 2.7960747d0
      cc(2,2) = 74.9062731d0
      cc(3,2) = 264.3323823d0
      cc(4,2) = 102.3427683d0
      cc(5,2) = 30.3054582d0
      cc(1,3) = 1.9719244d0
      cc(2,3) = 70.0368448d0
      cc(3,3) = 286.7645766d0
      cc(4,3) = 141.4003647d0
      cc(5,3) = 45.0490130d0
      cc(6,3) = 16.4655317d0
      cc(1,4) = 3.0725115d0
      cc(2,4) = 43.1809603d0
      cc(3,4) = 302.7696484d0
      cc(4,4) = 116.3041971d0
      cc(5,4) = -59.3066449d0
      cc(6,4) = -0.2462594d0
      cc(1,5) = 3.9888161d0
      cc(2,5) = 54.4481186d0
      cc(3,5) = 89.8756016d0
      cc(4,5) = 187.3744788d0
      cc(5,5) = 15.3552927d0
      cc(6,5) = -0.6439303d0
c
      a(1,1) = 373.7910783d0
      a(2,1) = 99.0222044d0
      a(3,1) = 17.6844640d0
      a(4,1) = 8.6281221d0
      a(5,1) = 3.0353877d0
      a(6,1) = 1.1191950d0
      a(1,2) = 169.8812365d0
      a(2,2) = 58.2434958d0
      a(3,2) = 21.7281484d0
      a(4,2) = 6.6359251d0
      a(5,2) = 1.4615973d0
      a(1,3) = 248.9438809d0
      a(2,3) = 84.6601414d0
      a(3,3) = 33.9674109d0
      a(4,3) = 11.5239104d0
      a(5,3) = 3.5354358d0
      a(6,3) = 0.9812915d0
      a(1,4) = 86.2036829d0
      a(2,4) = 76.0280558d0
      a(3,4) = 26.2846764d0
      a(4,4) = 3.3191097d0
      a(5,4) = 2.7009368d0
      a(6,4) = 0.5492105d0
      a(1,5) = 115.9461675d0
      a(2,5) = 43.7443988d0
      a(3,5) = 46.4937231d0
      a(4,5) = 13.2848999d0
      a(5,5) = 2.5347097d0
      a(6,5) = 0.7670500d0
      n(1,1) = 0
      n(2,1) = 1
      n(3,1) = 1
      n(4,1) = 2
      n(5,1) = 2
      n(6,1) = 2
      n(1,2) = 0
      n(2,2) = 1
      n(3,2) = 2
      n(4,2) = 2
      n(5,2) = 2
      n(1,3) = 0
      n(2,3) = 1
      n(3,3) = 2
      n(4,3) = 2
      n(5,3) = 2
      n(6,3) = 2
      n(1,4) = 0
      n(2,4) = 1
      n(3,4) = 2
      n(4,4) = 2
      n(5,4) = 2
      n(6,4) = 2
      n(1,5) = 0
      n(2,5) = 1
      n(3,5) = 2
      n(4,5) = 2
      n(5,5) = 2
      n(6,5) = 2
      endif
      go to 200
c     pb   hay
 180  if(olanl) then
       go to 300
      else
      lmx = 4
      nelcor = 78
      nt(1) = 6
      nt(2) = 6
      nt(3) = 6
      nt(4) = 6
      nt(5) = 6
      cc(1,1) = -0.1789605d0
      cc(2,1) = -54.3972337d0
      cc(3,1) = -199.7061759d0
      cc(4,1) = -79.1223941d0
      cc(5,1) = -24.9869020d0
      cc(6,1) = -4.4397939d0
      cc(1,2) = 2.8115386d0
      cc(2,2) = 65.0367205d0
      cc(3,2) = 212.7868545d0
      cc(4,2) = 72.1053175d0
      cc(5,2) = 33.0140940d0
      cc(6,2) = -5.7708461d0
      cc(1,3) = 4.8754911d0
      cc(2,3) = 63.9148102d0
      cc(3,3) = 148.1064358d0
      cc(4,3) = 47.3106301d0
      cc(5,3) = 21.0306702d0
      cc(6,3) = -7.0930772d0
      cc(1,4) = 3.2161388d0
      cc(2,4) = 55.7386086d0
      cc(3,4) = 121.4168351d0
      cc(4,4) = 19.3456064d0
      cc(5,4) = 15.3675168d0
      cc(6,4) = 6.1298724d0
      cc(1,5) = 4.1353682d0
      cc(2,5) = 67.5128446d0
      cc(3,5) = 258.7373107d0
      cc(4,5) = 113.2478264d0
      cc(5,5) = 34.1680201d0
      cc(6,5) = -6.5531956d0
c
      a(1,1) = 376.5803786d0
      a(2,1) = 86.4840014d0
      a(3,1) = 26.6784276d0
      a(4,1) = 9.4261986d0
      a(5,1) = 2.7101719d0
      a(6,1) = 0.8792031d0
      a(1,2) = 132.4248796d0
      a(2,2) = 47.2376044d0
      a(3,2) = 17.6312727d0
      a(4,2) = 5.4744712d0
      a(5,2) = 1.2634856d0
      a(6,2) = 0.7651447d0
      a(1,3) = 67.8966454d0
      a(2,3) = 24.9898225d0
      a(3,3) = 10.7052939d0
      a(4,3) = 3.2792568d0
      a(5,3) = 0.8452522d0
      a(6,3) = 0.6416245d0
      a(1,4) = 68.8336005d0
      a(2,4) = 24.2815874d0
      a(3,4) = 9.4532762d0
      a(4,4) = 2.4788185d0
      a(5,4) = 2.4789161d0
      a(6,4) = 0.5551738d0
      a(1,5) = 128.1021322d0
      a(2,5) = 54.8029154d0
      a(3,5) = 24.5529308d0
      a(4,5) = 8.1144792d0
      a(5,5) = 1.6931290d0
      a(6,5) = 0.7670500d0
      n(1,1) = 0
      n(2,1) = 1
      n(3,1) = 2
      n(4,1) = 2
      n(5,1) = 2
      n(6,1) = 2
      n(1,2) = 0
      n(2,2) = 1
      n(3,2) = 2
      n(4,2) = 2
      n(5,2) = 2
      n(6,2) = 2
      n(1,3) = 0
      n(2,3) = 1
      n(3,3) = 2
      n(4,3) = 2
      n(5,3) = 2
      n(6,3) = 2
      n(1,4) = 0
      n(2,4) = 1
      n(3,4) = 2
      n(4,4) = 2
      n(5,4) = 2
      n(6,4) = 2
      n(1,5) = 0
      n(2,5) = 1
      n(3,5) = 2
      n(4,5) = 2
      n(5,5) = 2
      n(6,5) = 2
      endif
      go to 200
c     bi   hay
 190  if(olanl) then
       go to 300
      else
      lmx = 4
      nelcor = 78
      nt(1) = 6
      nt(2) = 6
      nt(3) = 6
      nt(4) = 6
      nt(5) = 6
      cc(1,1) = -0.1833520d0
      cc(2,1) = -53.3389492d0
      cc(3,1) = -201.9117841d0
      cc(4,1) = -76.3115295d0
      cc(5,1) = -24.2578468d0
      cc(6,1) = -4.0319700d0
      cc(1,2) = 2.7916510d0
      cc(2,2) = 71.8553817d0
      cc(3,2) = 252.4327444d0
      cc(4,2) = 89.6045109d0
      cc(5,2) = 61.1744136d0
      cc(6,2) = -33.0471906d0
      cc(1,3) = 4.8462914d0
      cc(2,3) = 65.6883396d0
      cc(3,3) = 166.2619673d0
      cc(4,3) = 54.1393329d0
      cc(5,3) = 28.6204880d0
      cc(6,3) = -12.3452204d0
      cc(1,4) = 3.2183832d0
      cc(2,4) = 54.8952681d0
      cc(3,4) = 116.5540708d0
      cc(4,4) = 69.0760662d0
      cc(5,4) = -72.6410870d0
      cc(6,4) = 38.9400055d0
      cc(1,5) = 3.9864075d0
      cc(2,5) = 63.0262078d0
      cc(3,5) = 92.5259959d0
      cc(4,5) = 240.9973672d0
      cc(5,5) = 62.1844577d0
      cc(6,5) = -1.7251134d0
c
      a(1,1) = 223.3091971d0
      a(2,1) = 86.4699216d0
      a(3,1) = 26.8492203d0
      a(4,1) = 9.3823622d0
      a(5,1) = 2.7614708d0
      a(6,1) = 0.9423923d0
      a(1,2) = 155.9783914d0
      a(2,2) = 54.4591139d0
      a(3,2) = 20.1794558d0
      a(4,2) = 6.1273104d0
      a(5,2) = 1.1805983d0
      a(6,2) = 0.9932464d0
      a(1,3) = 74.5091386d0
      a(2,3) = 28.0487647d0
      a(3,3) = 12.0168508d0
      a(4,3) = 3.6914327d0
      a(5,3) = 0.9080824d0
      a(6,3) = 0.7394617d0
      a(1,4) = 68.0856679d0
      a(2,4) = 23.7250355d0
      a(3,4) = 9.2437053d0
      a(4,4) = 1.8452292d0
      a(5,4) = 1.2209859d0
      a(6,4) = 0.8740861d0
      a(1,5) = 134.3538983d0
      a(2,5) = 51.5812797d0
      a(3,5) = 54.9307284d0
      a(4,5) = 17.0475346d0
      a(5,5) = 3.7702824d0
      a(6,5) = 0.7670500d0
      n(1,1) = 0
      n(2,1) = 1
      n(3,1) = 2
      n(4,1) = 2
      n(5,1) = 2
      n(6,1) = 2
      n(1,2) = 0
      n(2,2) = 1
      n(3,2) = 2
      n(4,2) = 2
      n(5,2) = 2
      n(6,2) = 2
      n(1,3) = 0
      n(2,3) = 1
      n(3,3) = 2
      n(4,3) = 2
      n(5,3) = 2
      n(6,3) = 2
      n(1,4) = 0
      n(2,4) = 1
      n(3,4) = 2
      n(4,4) = 2
      n(5,4) = 2
      n(6,4) = 2
      n(1,5) = 0
      n(2,5) = 1
      n(3,5) = 2
      n(4,5) = 2
      n(5,5) = 2
      n(6,5) = 2
      endif
c
c     load ecp data
c
 200  return
 300  call caserr2('no library available for requested LANL2 ecp')
 6010 format (1x,'*** data not found for ecp type ',a6)
      end
**==ecplb4.f
      subroutine ecplb4 (zinp,lmx,nelcor,nt,n,cc,a,olanl)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character * 8 (z)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      dimension nt(6), n(10,6), cc(10,6), a(10,6)
c
      dimension ztit(3)
      data maxtyp/3/
      data ztit/'u','np','pu' /
c
c     this routine contains a library of ecp data taken from
c     various literature sources.
c     the ecps for
c     U and Np are taken from hay,wadt, jcp 82 (1985) 270,284
c
      if (.not.olanl) then
       call caserr2('no library available for requested LANL ecp')
      endif
c
      do 20 ityp = 1 , maxtyp
         if (zinp.eq.ztit(ityp).and.olanl) go to 40
 20   continue
 30   write (iwr,6010) zinp
      call caserr2('no library available for requested ecp')
 40   go to (92,93,30),ityp
c
c     u  hay
c
 92   continue
c
      cc( 1, 1) =         -0.22527100d0
      cc( 2, 1) =        -53.78057010d0
      cc( 3, 1) =       -516.30520530d0
      cc( 4, 1) =       -202.67382220d0
      cc( 5, 1) =        -62.56440490d0
      cc( 6, 1) =        -21.55895590d0
      cc( 7, 1) =         -1.53905450d0
      a ( 1, 1) =          4.51627250d0
      a ( 2, 1) =        368.68769520d0
      a ( 3, 1) =        112.75767080d0
      a ( 4, 1) =         30.20329150d0
      a ( 5, 1) =         10.29783550d0
      a ( 6, 1) =          3.49105150d0
      a ( 7, 1) =          1.26683910d0
      cc( 1, 2) =          2.61188160d0
      cc( 2, 2) =         81.23566570d0
      cc( 3, 2) =        353.24653640d0
      cc( 4, 2) =        311.92922340d0
      cc( 5, 2) =       -423.55257950d0
      cc( 6, 2) =        293.21608360d0
      cc( 7, 2) =        -21.34425520d0
      cc( 8, 2) =          1.66925480d0
      a ( 1, 2) =        257.35294420d0
      a ( 2, 2) =         83.16054280d0
      a ( 3, 2) =         29.72005130d0
      a ( 4, 2) =          6.24821320d0
      a ( 5, 2) =          4.16346720d0
      a ( 6, 2) =          3.00015110d0
      a ( 7, 2) =          1.56000000d0
      a ( 8, 2) =          0.87000000d0
      cc( 1, 3) =          1.84143880d0
      cc( 2, 3) =         75.70547550d0
      cc( 3, 3) =        288.46016740d0
      cc( 4, 3) =        250.72894290d0
      cc( 5, 3) =       -318.89680610d0
      cc( 6, 3) =        160.08449240d0
      cc( 7, 3) =         33.78018910d0
      cc( 8, 3) =         -2.06999970d0
      a ( 1, 3) =        288.72982500d0
      a ( 2, 3) =         76.78384200d0
      a ( 3, 3) =         24.82520030d0
      a ( 4, 3) =          5.12487350d0
      a ( 5, 3) =          3.76033000d0
      a ( 6, 3) =          2.97617060d0
      a ( 7, 3) =          1.99190170d0
      a ( 8, 3) =          0.97388320d0
      cc( 1, 4) =          2.90820640d0
      cc( 2, 4) =         74.35776280d0
      cc( 3, 4) =        587.11403290d0
      cc( 4, 4) =        317.25436790d0
      cc( 5, 4) =        102.84221170d0
      cc( 6, 4) =         22.42142780d0
      a ( 1, 4) =        350.46994070d0
      a ( 2, 4) =        166.70913010d0
      a ( 3, 4) =         85.82061190d0
      a ( 4, 4) =         26.98422550d0
      a ( 5, 4) =          7.15655920d0
      a ( 6, 4) =          1.28408550d0
      cc( 1, 5) =          3.93503390d0
      cc( 2, 5) =         64.94341200d0
      cc( 3, 5) =        265.18047760d0
      cc( 4, 5) =        389.42840900d0
      cc( 5, 5) =        153.47672210d0
      cc( 6, 5) =         12.08074440d0
      cc( 7, 5) =         -1.35101380d0
      a ( 1, 5) =        251.26707540d0
      a ( 2, 5) =        108.01092180d0
      a ( 3, 5) =        114.80218950d0
      a ( 4, 5) =         39.91125450d0
      a ( 5, 5) =         11.90397840d0
      a ( 6, 5) =          3.38507720d0
      a ( 7, 5) =          1.26776070d0
c
      go to 300
c
c     np  hay
c
 93   continue
c
      cc( 1, 1) =         -0.23019480d0
      cc( 2, 1) =        -45.50012810d0
      cc( 3, 1) =       -476.56537270d0
      cc( 4, 1) =       -195.84771820d0
      cc( 5, 1) =        -57.80336640d0
      cc( 6, 1) =        -19.74346000d0
      cc( 7, 1) =         -1.06696490d0
      a ( 1, 1) =          1.04312730d0
      a ( 2, 1) =        317.02760680d0
      a ( 3, 1) =        105.92003400d0
      a ( 4, 1) =         29.21133280d0
      a ( 5, 1) =          9.93276880d0
      a ( 6, 1) =          3.50434480d0
      a ( 7, 1) =          1.21340410d0
      cc( 1, 2) =          2.60187280d0
      cc( 2, 2) =         85.30283220d0
      cc( 3, 2) =        372.58843850d0
      cc( 4, 2) =        355.03931070d0
      cc( 5, 2) =       -432.40889660d0
      cc( 6, 2) =        302.44980930d0
      cc( 7, 2) =        -66.33727890d0
      cc( 8, 2) =         12.20638630d0
      a ( 1, 2) =        196.01142360d0
      a ( 2, 2) =         73.69599680d0
      a ( 3, 2) =         27.95190760d0
      a ( 4, 2) =          6.09350220d0
      a ( 5, 2) =          4.25840970d0
      a ( 6, 2) =          2.84483580d0
      a ( 7, 2) =          1.78183170d0
      a ( 8, 2) =          1.31267430d0
      cc( 1, 3) =          1.83776320d0
      cc( 2, 3) =         86.75593440d0
      cc( 3, 3) =        326.30120140d0
      cc( 4, 3) =        387.47452740d0
      cc( 5, 3) =       -616.35883850d0
      cc( 6, 3) =        367.04762110d0
      cc( 7, 3) =         -9.15507810d0
      cc( 8, 3) =          1.62181410d0
      a ( 1, 3) =        217.87574760d0
      a ( 2, 3) =         66.65246910d0
      a ( 3, 3) =         22.54723630d0
      a ( 4, 3) =          4.51443610d0
      a ( 5, 3) =          3.42615820d0
      a ( 6, 3) =          2.75484600d0
      a ( 7, 3) =          1.32590980d0
      a ( 8, 3) =          0.99468930d0
      cc( 1, 4) =          2.90616090d0
      cc( 2, 4) =         64.04786880d0
      cc( 3, 4) =        474.28555990d0
      cc( 4, 4) =        299.01152200d0
      cc( 5, 4) =         97.81467960d0
      cc( 6, 4) =         25.03463510d0
      a ( 1, 4) =        250.49997840d0
      a ( 2, 4) =        130.40328280d0
      a ( 3, 4) =         75.75346150d0
      a ( 4, 4) =         25.45525460d0
      a ( 5, 4) =          6.90935500d0
      a ( 6, 4) =          1.36862660d0
      cc( 1, 5) =          3.94737410d0
      cc( 2, 5) =         52.49626330d0
      cc( 3, 5) =        488.64808040d0
      cc( 4, 5) =        462.06653820d0
      cc( 5, 5) =        194.58025030d0
      cc( 6, 5) =         21.99611430d0
      cc( 7, 5) =         -1.83792830d0
      a ( 1, 5) =        283.15616940d0
      a ( 2, 5) =        135.85859240d0
      a ( 3, 5) =        146.72244560d0
      a ( 4, 5) =         51.06453820d0
      a ( 5, 5) =         15.29032080d0
      a ( 6, 5) =          4.11604960d0
      a ( 7, 5) =          1.51189600d0
c
 300  n ( 1, 1) =   0
      n ( 2, 1) =   1
      n ( 3, 1) =   2
      n ( 4, 1) =   2
      n ( 5, 1) =   2
      n ( 6, 1) =   2
      n ( 7, 1) =   2
      nt( 1)    =   7
      n ( 1, 2) =   0
      n ( 2, 2) =   1
      n ( 3, 2) =   2
      n ( 4, 2) =   2
      n ( 5, 2) =   2
      n ( 6, 2) =   2
      n ( 7, 2) =   2
      n ( 8, 2) =   2
      nt( 2)    =   8
      n ( 1, 3) =   0
      n ( 2, 3) =   1
      n ( 3, 3) =   2
      n ( 4, 3) =   2
      n ( 5, 3) =   2
      n ( 6, 3) =   2
      n ( 7, 3) =   2
      n ( 8, 3) =   2
      nt( 3)    =   8
      n ( 1, 4) =   0
      n ( 2, 4) =   1
      n ( 3, 4) =   2
      n ( 4, 4) =   2
      n ( 5, 4) =   2
      n ( 6, 4) =   2
      nt( 4)    =   6
      n ( 1, 5) =   0
      n ( 2, 5) =   1
      n ( 3, 5) =   2
      n ( 4, 5) =   2
      n ( 5, 5) =   2
      n ( 6, 5) =   2
      n ( 7, 5) =   2
      nt( 5)    =   7
      nelcor    =  78
      lmx       =   4
c
      return
 6010 format (1x,'*** data not found for ecp type ',a6)
      end
**==sbkp4.f
      subroutine sbkp4 (zinp,lmx,nelcor,nt,n,c,a)
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
      implicit character * 8 (z)
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
      dimension nt(6), n(10,6), c(10,6), a(10,6)
c
      dimension ztit(8)
      data maxtyp/8/
      data ztit/'k','ca','ga','ge','as','se','br','kr'/
c
c     this routine contains a library of ecp data taken from
c     various literature sources.
c     4th period main group sbk potentials
c     
      do 10 ityp = 1 , maxtyp
         if (zinp.eq.ztit(ityp)) go to 20
 10   continue
      write (iwr,6010) zinp
      call caserr2('no library available for requested ecp')
c
c
 20   go to (100,200,300,400,500,600,700,800), ityp
c
c     potassium (K_SBKJC_VDZ_ECP)
c
 100  nelcor = 18
      lmx = 2
      nt(1) = 1
      nt(2) = 2
      nt(3) = 2
      n(1,1) = 1
      n(1,2) = 0
      n(2,2) = 2
      n(1,3) = 0
      n(2,3) = 2
c
      c(1,1)=-4.30002d+00
      c(1,2)= 6.71530d+00
      c(2,2)= 8.46683d+00
      c(1,3)= 7.06252d+00
      c(2,3)= 2.69868d+00
      a(1,1)= 0.49652d+00
      a(1,2)= 2.48731d+00
      a(2,2)= 0.60190d+00
      a(1,3)= 0.94170d+00
      a(2,3)= 0.33623d+00
      return
c
c     calcium (Ca_SBKJC_VDZ_ECP)
c
  200 nelcor = 18
      lmx = 2
      nt(1) = 2
      nt(2) = 2
      nt(3) = 2
      n(1,1) = 1
      n(2,1) = 1
      n(1,2) = 0
      n(2,2) = 2
      n(1,3) = 0
      n(2,3) = 2
c
      c(1,1)=-5.81970d+00
      c(2,1)=-2.16548d+00
      c(1,2)= 8.00522d+00
      c(2,2)= 8.94751d+00
      c(1,3)= 8.66222d+00
      c(2,3)= 5.57218d+00
      a(1,1)= 2.38045d+00
      a(2,1)= 0.40639d+00
      a(1,2)= 2.10543d+00
      a(2,2)= 0.66865d+00
      a(1,3)= 1.79413d+00
      a(2,3)= 0.49323d+00
      return
c
c        gallium
c     the semi-core potential is stored with the transition metals.
c
  300 call caserr2('gallium potential with transition metals')
c
c     germanium (Ge_SBKJC_VDZ_ECP)
c
  400 continue
      c(1,1)=  -1.26859d+00
      a(1,1)=   0.82751d+00
      c(1,2)=   6.77027d+00
      a(1,2)=  12.24839d+00
      c(2,2)=  27.21907d+00
      a(2,2)=   2.41291d+00
      c(1,3)=   8.31253d+00
      a(1,3)=  15.55556d+00
      c(2,3)=  15.53292d+00
      a(2,3)=   1.82031d+00
      c(1,4)=   1.62360d+00
      a(1,4)=   3.06970d+00
      c(2,4)=   3.20660d+00
      a(2,4)=   0.80525d+00
      go to 250
c
c     arsenic (As_SBKJC_VDZ_ECP)
c
  500 continue
      c(1,1)=   -4.77028d+00
      a(1,1)=    1.64583d+00
      c(1,2)=   13.03398d+00
      a(1,2)=   20.76160d+00
      c(2,2)=   43.12429d+00
      a(2,2)=    2.80084d+00
      c(1,3)=    6.44063d+00
      a(1,3)=   16.05834d+00
      c(2,3)=   28.12117d+00
      a(2,3)=    2.29432d+00
      c(1,4)=    2.23508d+00
      a(1,4)=    1.12086d+00
      c(2,4)=    4.60721d+00
      a(2,4)=    1.12140d+00
      go to 250
c
c     selenium (Se_SBKJC_VDZ_ECP)
c
  600 continue
      c(1,1)=   -5.58498d+00
      a(1,1)=    2.09356d+00
      c(1,2)=   13.28081d+00
      a(1,2)=   20.55644d+00
      c(2,2)=   51.00011d+00
      a(2,2)=    3.15490d+00
      c(1,3)=    1.79640d+00
      a(1,3)=    2.58495d+00
      c(2,3)=   30.43787d+00
      a(2,3)=    2.58526d+00
      c(1,4)=    2.40001d+00
      a(1,4)=    1.38290d+00
      c(2,4)=    6.59948d+00
      a(2,4)=    1.42671d+00
      go to 250
c
c     bromine( Br_SBKJC_VDZ_ECP)
c
  700 continue
      c(1,1)=   -4.78206d+00
      a(1,1)=    2.12116d+00
      c(1,2)=    2.60190d+00
      a(1,2)=    3.21538d+00
      c(2,2)=   44.15165d+00
      a(2,2)=    3.21639d+00
      c(1,3)=    2.09478d+00
      a(1,3)=    2.72704d+00
      c(2,3)=   30.59546d+00
      a(2,3)=    2.72762d+00
      c(1,4)=    1.79375d+00
      a(1,4)=    1.73168d+00
      c(2,4)=   10.04660d+00
      a(2,4)=    1.73246d+00
      go to 250
c
c     krypton (Kr_SBKJC_VDZ_ECP)
c
  800 continue
      c(1,1)=   -5.61220d+00
      a(1,1)=    2.64295d+00
      c(1,2)=    2.36066d+00
      a(1,2)=    3.59276d+00
      c(2,2)=   51.88629d+00
      a(2,2)=    3.57933d+00
      c(1,3)=    2.34982d+00
      a(1,3)=    3.07132d+00
      c(2,3)=   35.67029d+00
      a(2,3)=    3.07238d+00
      c(1,4)=    2.85459d+00
      a(1,4)=    7.08677d+00
      c(2,4)=   16.29386d+00
      a(2,4)=    2.15058d+00
c
250   nelcor = 28
      lmx = 3
      nt(1) = 1
      nt(2) = 2
      nt(3) = 2
      nt(4) = 2
      n(1,1) = 1
      n(1,2) = 0
      n(2,2) = 2
      n(1,3) = 0
      n(2,3) = 2
      n(1,4) = 0
      n(2,4) = 2
      return
 6010 format (1x,'*** data not found for ecp type ',a4)
      end
**==sbkp5.f
      subroutine sbkp5 (zinp,lmx,nelcor,nt,n,c,a)
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
      implicit character * 8 (z)
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
      dimension nt(6), n(10,6), c(10,6), a(10,6)
c
      dimension ztit(7)
      data maxtyp/7/
      data ztit/'rb','sr','sn','sb','te','i','xe'/
c
c     this routine contains a library of ecp data taken from
c     various literature sources.
c     5th period main group sjkb potentials
c     
      do 10 ityp = 1 , maxtyp
         if (zinp.eq.ztit(ityp)) go to 20
 10   continue
      write (iwr,6010) zinp
      call caserr2('no library available for requested ecp')
c
c
 20   go to (100,200,400,500,600,700,800), ityp
c
c     rubidium (Rb_SBKJC_VDZ_ECP)
c
  100 continue
      c(1,1)=  -3.23942d+00
      a(1,1)=   0.36395d+00
      c(1,2)=  -2.49855d+00
      a(1,2)=   0.21048d+00
      c(2,2)=  26.37428d+00
      a(2,2)=   0.18693d+00
      c(3,2)= -21.20813d+00
      a(3,2)=   0.40031d+00
      c(1,3)=   2.48839d+00
      a(1,3)=   0.26539d+00
      c(2,3)=   5.46002d+00
      a(2,3)=   0.69850d+00
      c(1,4)=  -4.38861d+00
      a(1,4)=   0.13507d+00
      c(2,4)=  28.32654d+00
      a(2,4)=   0.08272d+00
      c(3,4)= -26.59458d+00
      a(3,4)=   0.26365d+00
      go to 150
c
c     strontium (Sr_SBKJC_VDZ_ECP)
c
  200 continue
      c(1,1)=  -4.13596d+00
      a(1,1)=   0.51676d+00
      c(1,2)=  -3.79121d+00
      a(1,2)=   0.23231d+00
      c(2,2)=  24.16768d+00
      a(2,2)=   0.18029d+00
      c(3,2)= -17.38542d+00
      a(3,2)=   0.52170d+00
      c(1,3)=   5.46439d+00
      a(1,3)=   0.40227d+00
      c(2,3)=  10.49841d+00
      a(2,3)=   1.30074d+00
      c(1,4)=  -1.73727d+00
      a(1,4)=   0.26949d+00
      c(2,4)=  37.43956d+00
      a(2,4)=   0.25550d+00
      c(3,4)= -36.00393d+00
      a(3,4)=   0.28868d+00
 150  nelcor = 36
      lmx = 3
      nt(1) = 1
      nt(2) = 3
      nt(3) = 2
      nt(4) = 3
      n(1,1) = 1
      n(1,2) = 2
      n(2,2) = 0
      n(3,2) = 0
      n(1,3) = 2
      n(2,3) = 0
      n(1,4) = 2
      n(2,4) = 0
      n(3,4) = 0
      return
c
c     tin (Sn_SBKJC_VDZ_ECP)
c
  400 continue
      c(1,1)=  -5.60682d+00
      a(1,1)=   0.82050d+00
      c(1,2)=   7.18487d+00
      a(1,2)=   2.60504d+00
      c(2,2)= -73.06581d+00
      a(2,2)=   0.95780d+00
      c(3,2)=  91.84771d+00
      a(3,2)=   1.01356d+00
      c(1,3)=   5.07311d+00
      a(1,3)=   0.67947d+00
      c(2,3)= -87.03674d+00
      a(2,3)=   0.74276d+00
      c(3,3)=  95.37048d+00
      a(3,3)=   0.76292d+00
      c(1,4)=   5.10243d+00
      a(1,4)=   1.65430d+00
      c(2,4)=   5.27466d+00
      a(2,4)=   0.60301d+00
c
      nelcor = 46
      lmx = 3
      nt(1) = 1
      nt(2) = 3
      nt(3) = 3
      nt(4) = 2
      n(1,1) = 1
      n(1,2) = 0
      n(2,2) = 2
      n(3,2) = 2
      n(1,3) = 0
      n(2,3) = 2
      n(3,3) = 2
      n(1,4) = 0
      n(2,4) = 2
c
      return
c
c      antimony (Sb_SBKJC_VDZ_ECP)
c
  500 continue
      c(1,1)=  -2.90951d+00
      a(1,1)=   0.75490d+00
      c(2,1)= -12.44921d+00
      a(2,1)=   2.83072d+00
      c(1,2)=   7.55592d+00
      a(1,2)=   5.86819d+00
      c(2,2)= -53.01395d+00
      a(2,2)=   1.10036d+00
      c(3,2)=  79.99053d+00
      a(3,2)=   1.21790d+00
      c(1,3)=   8.30360d+00
      a(1,3)=   3.51430d+00
      c(2,3)= -29.04708d+00
      a(2,3)=   0.86461d+00
      c(3,3)=  45.56610d+00
      a(3,3)=   0.96029d+00
      c(1,4)=   9.22182d+00
      a(1,4)=   2.45981d+00
      c(2,4)=   7.92129d+00
      a(2,4)=   0.76415d+00
      go to 250
c
c     tellurium (Te_SBKJC_VDZ_ECP)
c
  600 continue
      c(1,1)=  -3.95984d+00
      a(1,1)=   0.94847d+00
      c(2,1)= -13.42871d+00
      a(2,1)=   3.94682d+00
      c(1,2)=  12.95665d+00
      a(1,2)=   3.49400d+00
      c(2,2)= -47.45882d+00
      a(2,2)=   1.21477d+00
      c(3,2)=  72.77514d+00
      a(3,2)=   1.33530d+00
      c(1,3)=   8.47921d+00
      a(1,3)=   1.37394d+00
      c(2,3)= -39.40684d+00
      a(2,3)=   1.27761d+00
      c(3,3)=  51.93673d+00
      a(3,3)=   1.28425d+00
      c(1,4)=  11.54663d+00
      a(1,4)=   3.66618d+00
      c(2,4)=  11.23776d+00
      a(2,4)=   0.93189d+00
      go to 250
c
c     iodine (I_SBKJC_VDZ_ECP)
c
  700 continue
      c(1,1)=  -3.69639d+00
      a(1,1)=   0.97628d+00
      c(2,1)= -14.00305d+00
      a(2,1)=   4.33343d+00
      c(1,2)=  12.11123d+00
      a(1,2)=   4.13071d+00
      c(2,2)= -41.09206d+00
      a(2,2)=   1.33375d+00
      c(3,2)=  70.73761d+00
      a(3,2)=   1.49121d+00
      c(1,3)=  10.59271d+00
      a(1,3)=   3.04692d+00
      c(2,3)= -46.02273d+00
      a(2,3)=   1.06339d+00
      c(3,3)=  65.05047d+00
      a(3,3)=   1.14405d+00
      c(1,4)=   9.73089d+00
      a(1,4)=   3.93063d+00
      c(2,4)=  13.98880d+00
      a(2,4)=   1.06920d+00
      go to 250
c
c     xenon (Xe_SBKJC_VDZ_ECP)
c
  800 continue
      c(1,1)=  -5.10538d+00
      a(1,1)=   1.27068d+00
      c(2,1)= -14.39873d+00
      a(2,1)=   5.90469d+00
      c(1,2)=  13.44278d+00
      a(1,2)=   7.15744d+00
      c(2,2)= -45.48918d+00
      a(2,2)=   1.45217d+00
      c(3,2)=  83.21627d+00
      a(3,2)=   1.65478d+00
      c(1,3)=  12.73758d+00
      a(1,3)=   4.36144d+00
      c(2,3)= -49.60137d+00
      a(2,3)=   1.25455d+00
      c(3,3)=  75.57733d+00
      a(3,3)=   1.36439d+00
      c(1,4)=  17.79041d+00
      a(1,4)=   5.89897d+00
      c(2,4)=  17.25358d+00
      a(2,4)=   1.23485d+00
c
  250 continue
      nelcor = 46
      lmx = 3
      nt(1) = 2
      nt(2) = 3
      nt(3) = 3
      nt(4) = 2
      n(1,1) = 1
      n(2,1) = 1
      n(1,2) = 0
      n(2,2) = 2
      n(3,2) = 2
      n(1,3) = 0
      n(2,3) = 2
      n(3,3) = 2
      n(1,4) = 0
      n(2,4) = 2
      return
 6010 format (1x,'*** data not found for ecp type ',a4)
      end
**==sbkp6.f
      subroutine sbkp6 (zinp,lmx,nelcor,nt,n,c,a)
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
      implicit character * 8 (z)
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
      dimension nt(6), n(10,6), c(10,6), a(10,6)
c
      dimension ztit(7)
      data maxtyp/7/
      data ztit/'cs','ba','pb','bi','po','at','rn'/
c
c     this routine contains a library of ecp data taken from
c     various literature sources.
c     5th period main group sjkb potentials
c     
      do 10 ityp = 1 , maxtyp
         if (zinp.eq.ztit(ityp)) go to 20
 10   continue
      write (iwr,6010) zinp
      call caserr2('no library available for requested ecp')
c
c
 20   go to (55,56,82,83,84,85,86), ityp
c
c     cesium (Cs_SBKJC_VDZ_ECP)
c
  55  continue
      nelcor = 54
      lmx = 3
      nt(1) = 1
      nt(2) = 3
      nt(3) = 2
      nt(4) = 3
      n(1,1) = 1
      n(1,2) = 0
      n(2,2) = 2
      n(3,2) = 2
      n(1,3) = 0
      n(2,3) = 2
      n(1,4) = 0
      n(2,4) = 2
      n(3,4) = 1
      c(1,1)=  -5.01874d+00
      a(1,1)=   0.31356d+00
      c(1,2)=   6.79629d+00
      a(1,2)=   0.40093d+00
      c(2,2)= -12.27571d+00
      a(2,2)=   0.65203d+00
      c(3,2)=  13.72442d+00
      a(3,2)=   0.45384d+00
      c(1,3)=  12.30018d+00
      a(1,3)=   0.67287d+00
      c(2,3)=   2.95222d+00
      a(2,3)=   0.22290d+00
      c(1,4)=   1.96793d+00
      a(1,4)=   0.40835d+00
      c(2,4)=  -0.80500d+00
      a(2,4)=   0.20154d+00
      c(3,4)=   2.47096d+00
      a(3,4)=   0.24413d+00
      return
c
c     barium (Ba_SBKJC_VDZ_ECP)
c
   56 continue
      nelcor = 54
      lmx = 3
      nt(1) = 2
      nt(2) = 3
      nt(3) = 2
      nt(4) = 2
      n(1,1) = 1
      n(2,1) = 1
      n(1,2) = 0
      n(2,2) = 2
      n(3,2) = 2
      n(1,3) = 0
      n(2,3) = 2
      n(1,4) = 1
      n(2,4) = 1
      c(1,1)= -4.89584d+00
      a(1,1)=  0.37866d+00
      c(2,1)=-17.47201d+00
      a(2,1)=  3.17423d+00
      c(1,2)= 13.46017d+00
      a(1,2)=  0.40687d+00
      c(2,2)=-14.33523d+00
      a(2,2)=  0.94266d+00
      c(3,2)=  8.12634d+00
      a(3,2)=  0.46605d+00
      c(1,3)= 17.03588d+00
      a(1,3)=  1.17727d+00
      c(2,3)=  4.74706d+00
      a(2,3)=  0.28981d+00
      c(1,4)= 12.78807d+00
      a(1,4)=  1.00450d+00
      c(2,4)= -0.28677d+00
      a(2,4)=  0.14394d+00
      return
c
c     lead (Pb_SBKJC_VDZ_ECP)
c
   82 continue
      n(1,1) = 1
      n(1,2) = 0
      n(2,2) = 2
      n(3,2) = 2
      n(1,3) = 0
      n(2,3) = 2
      n(3,3) = 2
      n(1,4) = 0
      n(2,4) = 2
      n(1,5) = 0
      n(2,5) = 2
      n(3,5) = 2
      c(1,1)=   -2.76464d+00
      a(1,1)=    0.51756d+00
      c(1,2)=    6.94774d+00
      a(1,2)=    2.44363d+00
      c(2,2)=  -17.72106d+00
      a(2,2)=    0.82022d+00
      c(3,2)=   32.39564d+00
      a(3,2)=    0.97662d+00
      c(1,3)=    4.94619d+00
      a(1,3)=    1.49041d+00
      c(2,3)=  -17.16033d+00
      a(2,3)=    0.59861d+00
      c(3,3)=   26.06040d+00
      a(3,3)=    0.66916d+00
      c(1,4)=    4.99762d+00
      a(1,4)=    2.19203d+00
      c(2,4)=    4.52730d+00
      a(2,4)=    0.46937d+00
      c(1,5)=    4.95827d+00
      a(1,5)=    0.69081d+00
      c(2,5)= -124.15758d+00
      a(2,5)=    0.20680d+00
      c(3,5)=  125.32798d+00
      a(3,5)=    0.20773d+00
      go to 100
c
c     bismuth (Bi_SBKJC_VDZ_ECP)
c
   83 continue
      n(1,1) = 1
      n(1,2) = 0
      n(2,2) = 2
      n(3,2) = 2
      n(1,3) = 0
      n(2,3) = 2
      n(3,3) = 2
      n(1,4) = 0
      n(2,4) = 2
      n(1,5) = 0
      n(2,5) = 0
      n(3,5) = 2
      c(1,1)=   -5.34991d+00
      a(1,1)=    0.77941d+00
      c(1,2)=    5.91224d+00
      a(1,2)=    0.58082d+00
      c(2,2)=  -46.52354d+00
      a(2,2)=    0.86480d+00
      c(3,2)=   58.30673d+00
      a(3,2)=    0.95682d+00
      c(1,3)=    7.58772d+00
      a(1,3)=    2.41860d+00
      c(2,3)=   -1.55133d+00
      a(2,3)=    0.51229d+00
      c(3,3)=   17.65062d+00
      a(3,3)=    0.98429d+00
      c(1,4)=    8.20935d+00
      a(1,4)=    4.10212d+00
      c(2,4)=    8.38970d+00
      a(2,4)=    0.63052d+00
      c(1,5)=  -31.22481d+00
      a(1,5)=    0.57108d+00
      c(2,5)=   31.99918d+00
      a(2,5)=    1.17643d+00
      c(3,5)=   16.82287d+00
      a(3,5)=    0.78832d+00
      go to 100
c
c     polonium (Po_SBKJC_VDZ_ECP)
c
   84 continue
      n(1,1) = 1
      n(1,2) = 0
      n(2,2) = 2
      n(3,2) = 2
      n(1,3) = 0
      n(2,3) = 2
      n(3,3) = 2
      n(1,4) = 0
      n(2,4) = 2
      n(1,5) = 0
      n(2,5) = 0
      n(3,5) = 2
      c(1,1)=   -6.91426d+00
      a(1,1)=    0.98236d+00
      c(1,2)=    5.38765d+00
      a(1,2)=    0.39388d+00
      c(2,2)=  -41.35987d+00
      a(2,2)=    0.84405d+00
      c(3,2)=   55.71735d+00
      a(3,2)=    0.99601d+00
      c(1,3)=   11.84449d+00
      a(1,3)=    4.87877d+00
      c(2,3)=   -2.29848d+00
      a(2,3)=    0.84347d+00
      c(3,3)=   27.85431d+00
      a(3,3)=    1.26544d+00
      c(1,4)=    9.51261d+00
      a(1,4)=    2.52803d+00
      c(2,4)=   11.38029d+00
      a(2,4)=    0.75864d+00
      c(1,5)=  -32.32727d+00
      a(1,5)=    0.44923d+00
      c(2,5)=   33.57540d+00
      a(2,5)=    0.71317d+00
      c(3,5)=    7.89213d+00
      a(3,5)=    0.56249d+00
c
  100 continue
      nelcor = 78
      lmx = 4
      nt(1) = 1
      nt(2) = 3
      nt(3) = 3
      nt(4) = 2
      nt(5) = 3
      return
c
c     astatine (At_SBKJC_VDZ_ECP)
c
   85 continue
      c(1,1)=   -8.22821d+00
      a(1,1)=    1.18142d+00
      c(1,2)=    6.65302d+00
      a(1,2)=    0.49276d+00
      c(2,2)=  -15.98654d+00
      a(2,2)=    0.85421d+00
      c(3,2)=   33.10060d+00
      a(3,2)=    1.27555d+00
      c(1,3)=   10.73147d+00
      a(1,3)=    3.37381d+00
      c(2,3)= -128.82729d+00
      a(2,3)=    0.93728d+00
      c(3,3)=  150.80821d+00
      a(3,3)=    0.97488d+00
      c(1,4)=   15.83166d+00
      a(1,4)=    3.88638d+00
      c(2,4)=   14.55150d+00
      a(2,4)=    0.87679d+00
      c(1,5)=    4.27097d+00
      a(1,5)=    1.12284d+00
      c(2,5)=   -3.39500d+00
      a(2,5)=    1.15530d+00
      go to 200
c
c     radon (Rn_SBKJC_VDZ_ECP)
c
   86 continue
      c(1,1)=   -9.36991d+00
      a(1,1)=    1.37920d+00
      c(1,2)=    5.25644d+00
      a(1,2)=    0.47169d+00
      c(2,2)=  -39.27598d+00
      a(2,2)=    1.01268d+00
      c(3,2)=   62.19287d+00
      a(3,2)=    1.26952d+00
      c(1,3)=    7.38631d+00
      a(1,3)=    3.32534d+00
      c(2,3)=  -45.67830d+00
      a(2,3)=    0.99950d+00
      c(3,3)=   71.81385d+00
      a(3,3)=    1.11620d+00
      c(1,4)=    4.90484d+00
      a(1,4)=    0.92214d+00
      c(2,4)=   12.61598d+00
      a(2,4)=    0.92221d+00
      c(1,5)=    2.96924d+00
      a(1,5)=    1.61735d+00
      c(2,5)=   -0.17844d+00
      a(2,5)=    0.32150d+00
c
  200 continue
      nelcor = 78
      lmx = 4
      nt(1) = 1
      nt(2) = 3
      nt(3) = 3
      nt(4) = 2
      nt(5) = 2
      n(1,1) = 1
      n(1,2) = 0
      n(2,2) = 2
      n(3,2) = 2
      n(1,3) = 0
      n(2,3) = 2
      n(3,3) = 2
      n(1,4) = 0
      n(2,4) = 2
      n(1,5) = 0
      n(2,5) = 2
      return
c
 6010 format (1x,'*** data not found for ecp type ',a4)
      end
**==sbklan.f
      subroutine sbklan (zinp,lmx,nelcor,nt,n,c,a)
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
      implicit character * 8 (z)
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
      dimension nt(6), n(10,6), c(10,6), a(10,6)
c
      dimension ztit(14)
      data maxtyp/14/
      data ztit/'ce','pr','nd','pm','sm','eu','gd','tb',
     +          'dy','ho','er','tm','yb','lu' /
c
c     this routine contains a library of ecp data for the lanthanides
c       from T.R. Cundari and W.J. Stevens, 
c       J. Chem. Phys. 98,  5555-5565.
c     
c       the variables in this routine correspond as follows to the
c       aforementioned paper:
c          zinp: the tag of the element
c          nelcor: the number of core electrons used in the ecp
c          n(i,j): this corresponds to the n(lk) parameter
c          c(i,j): A(lk) in the paper (Table I)
c          a(l,k): B(lk) in the paper (Table I)
c          nt(i): number of elements in the i'th part of the ecp
c
c       first set up nt, n, and icore which are the same for all
c       of the lanthanides.
c
      do 10 ityp = 1 , maxtyp
         if (zinp.eq.ztit(ityp)) go to 20
 10   continue
      write (iwr,6010) zinp
      call caserr2('no library available for requested ecp')
 20   nelcor = 46
      lmx = 3
      nt(1) = 2
      nt(2) = 3
      nt(3) = 3
      nt(4) = 2
      n(1,1) = 1
      n(2,1) = 1
      n(1,2) = 2
      n(2,2) = 2
      n(3,2) = 0
      n(1,3) = 0
      n(2,3) = 2
      n(3,3) = 2
      n(1,4) = 2
      n(2,4) = 0
c
      go to (58,59,60,61,62,63,64,65,66,67,68,69,70,71), ityp
c
c   cerium (Ce_SBKJC_VDZ_ECP)
c
 58   continue
      c(1,1) =  -15.34875610d+00
      a(1,1) =    9.20747930d+00
      c(2,1) =   -5.84323953d+00
      a(2,1) =    1.86730116d+00
      c(1,2) = -255.56238300d+00
      a(1,2) =    1.89370134d+00
      c(2,2) =  307.31392800d+00
      a(2,2) =    1.97914859d+00
      c(3,2) =   10.66990170d+00
      a(3,2) =   10.74296970d+00
      c(1,3) =   12.22921090d+00
      a(1,3) =    7.75592979d+00
      c(2,3) =  124.94246600d+00
      a(2,3) =    1.81564126d+00
      c(3,3) =  -84.59998680d+00
      a(3,3) =    1.67164719d+00
      c(1,4) =   24.94467550d+00
      a(1,4) =    1.70642047d+00
      c(2,4) =   10.28614470d+00
      a(2,4) =    6.48933742d+00
      return
c
c  praseodymium
c
 59   continue
      c(1,1) =  -15.44352190d+00
      a(1,1) =    9.82972856d+00
      c(2,1) =   -5.89611285d+00
      a(2,1) =    1.98070077d+00
      c(1,2) = -223.64398200d+00
      a(1,2) =    2.12955579d+00
      c(2,2) =  278.13451500d+00
      a(2,2) =    2.22746546d+00
      c(3,2) =   12.62107180d+00
      a(3,2) =    7.28371963d+00
      c(1,3) =   12.52563900d+00
      a(1,3) =    7.80928039d+00
      c(2,3) =  121.95278200d+00
      a(2,3) =    1.92977861d+00
      c(3,3) =  -79.40475120d+00
      a(3,3) =    1.76478925d+00
      c(1,4) =   26.19266520d+00
      a(1,4) =    1.79797432d+00
      c(2,4) =    9.87391206d+00
      a(2,4) =    6.54805153d+00
      return
c
c   neodymium
c
 60   continue
      c(1,1) =  -15.48300440d+00
      a(1,1) =   10.42978840d+00
      c(2,1) =   -5.94833367d+00
      a(2,1) =    2.09229128d+00
      c(1,2) = -219.47084000d+00
      a(1,2) =    2.22478561d+00
      c(2,2) =  280.52892800d+00
      a(2,2) =    2.34459586d+00
      c(3,2) =   11.76232500d+00
      a(3,2) =   11.85798300d+00
      c(1,3) =   11.43612170d+00
      a(1,3) =    9.44306451d+00
      c(2,3) =  120.33535600d+00
      a(2,3) =    1.98949076d+00
      c(3,3) =  -75.78196430d+00
      a(3,3) =    1.79805446d+00
      c(1,4) =   27.38307130d+00
      a(1,4) =    1.88855052d+00
      c(2,4) =    9.61741049d+00
      a(2,4) =    6.64366554d+00
      return
c
c   promethium
c
 61   continue
      c(1,1) =  -15.58675570d+00
      a(1,1) =   11.10670840d+00
      c(2,1) =   -6.01203457d+00
      a(2,1) =    2.21044787d+00
      c(1,2) = -215.58178700d+00
      a(1,2) =    2.27365312d+00
      c(2,2) =  277.53404200d+00
      a(2,2) =    2.40502952d+00
      c(3,2) =   11.55980120d+00
      a(3,2) =   10.76473100d+00
      c(1,3) =   10.33666460d+00
      a(1,3) =    9.85397700d+00
      c(2,3) =  121.46769300d+00
      a(2,3) =    2.07626890d+00
      c(3,3) =  -74.87664400d+00
      a(3,3) =    1.86493486d+00
      c(1,4) =   28.48611040d+00
      a(1,4) =    1.97798221d+00
      c(2,4) =    9.28090576d+00
      a(2,4) =    6.67329417d+00
      return
c
c   samarium
c
 62   continue
      c(1,1) =  -15.65862010d+00
      a(1,1) =   11.75812390d+00
      c(2,1) =   -6.05376496d+00
      a(2,1) =    2.32808297d+00
      c(1,2) = -206.06726600d+00
      a(1,2) =    2.37776611d+00
      c(2,2) =  270.34259800d+00
      a(2,2) =    2.52471590d+00
      c(3,2) =   11.44490440d+00
      a(3,2) =   10.26438980d+00
      c(1,3) =   10.69510070d+00
      a(1,3) =   10.00542360d+00
      c(2,3) =  121.87722800d+00
      a(2,3) =    2.20124125d+00
      c(3,3) =  -72.78839030d+00
      a(3,3) =    1.96809283d+00
      c(1,4) =   29.52394000d+00
      a(1,4) =    2.06720439d+00
      c(2,4) =    9.03717701d+00
      a(2,4) =    6.71288904d+00
      return
c
c   europium
c
 63   continue
      c(1,1) =  -15.72447990d+00
      a(1,1) =   12.41835390d+00
      c(2,1) =   -6.09107368d+00
      a(2,1) =    2.44755847d+00
      c(1,2) = -196.63773800d+00
      a(1,2) =    2.50382845d+00
      c(2,2) =  264.10344900d+00
      a(2,2) =    2.66926410d+00
      c(3,2) =   11.80698020d+00
      a(3,2) =   10.36738240d+00
      c(1,3) =   10.93806010d+00
      a(1,3) =   10.13418310d+00
      c(2,3) =  128.16026900d+00
      a(2,3) =    2.32439941d+00
      c(3,3) =  -76.51532760d+00
      a(3,3) =    2.08337487d+00
      c(1,4) =   30.46732640d+00
      a(1,4) =    2.15438942d+00
      c(2,4) =    8.74122693d+00
      a(2,4) =    6.69379732d+00
      return
c
c   gadolinium
c
 64   continue
      c(1,1) =  -15.78672600d+00
      a(1,1) =   13.09509080d+00
      c(2,1) =   -6.12777814d+00
      a(2,1) =    2.56973856d+00
      c(1,2) = -137.90212500d+00
      a(1,2) =    2.70226523d+00
      c(2,2) =  208.96454000d+00
      a(2,2) =    2.93347741d+00
      c(3,2) =   16.85869360d+00
      a(3,2) =    9.85080964d+00
      c(1,3) =   11.13083430d+00
      a(1,3) =   10.25620120d+00
      c(2,3) =  131.65132100d+00
      a(2,3) =    2.45273213d+00
      c(3,3) =  -77.44430080d+00
      a(3,3) =    2.19693427d+00
      c(1,4) =   31.46641180d+00
      a(1,4) =    2.24716923d+00
      c(2,4) =    8.34507286d+00
      a(2,4) =    6.43064369d+00
      return
c
c   terbium
c
 65   continue
      c(1,1) =  -15.84080810d+00
      a(1,1) =   13.78498800d+00
      c(2,1) =   -6.16678475d+00
      a(2,1) =    2.69386693d+00
      c(1,2) = -140.10322400d+00
      a(1,2) =    2.61817346d+00
      c(2,2) =  215.28993500d+00
      a(2,2) =    2.89041797d+00
      c(3,2) =   11.68123370d+00
      a(3,2) =   17.04897130d+00
      c(1,3) =   11.12162130d+00
      a(1,3) =   10.12533140d+00
      c(2,3) =  139.65713100d+00
      a(2,3) =    2.56570377d+00
      c(3,3) =  -83.30982210d+00
      a(3,3) =    2.30834841d+00
      c(1,4) =   32.12176790d+00
      a(1,4) =    2.32603332d+00
      c(2,4) =    8.22336256d+00
      a(2,4) =    6.59311802d+00
      return
c
c   dysoprosium
c
 66   continue
      c(1,1) =  -15.89723700d+00
      a(1,1) =   14.49058170d+00
      c(2,1) =   -6.19939399d+00
      a(2,1) =    2.82105135d+00
      c(1,2) = -130.00434800d+00
      a(1,2) =    2.73908105d+00
      c(2,2) =  208.63602400d+00
      a(2,2) =    3.04934714d+00
      c(3,2) =   11.48546660d+00
      a(3,2) =   16.57473300d+00
      c(1,3) =   11.62238130d+00
      a(1,3) =   10.69436010d+00
      c(2,3) =  125.41578200d+00
      a(2,3) =    2.74055834d+00
      c(3,3) =  -65.82038280d+00
      a(3,3) =    2.41123463d+00
      c(1,4) =   32.79772690d+00
      a(1,4) =    2.40961369d+00
      c(2,4) =    7.98023135d+00
      a(2,4) =    6.49588252d+00
      return
c
c   holmium
c
 67   continue
      c(1,1) =  -15.95165880d+00
      a(1,1) =   15.21414080d+00
      c(2,1) =   -6.23162568d+00
      a(2,1) =    2.95124916d+00
      c(1,2) = -110.58301100d+00
      a(1,2) =    2.84745682d+00
      c(2,2) =  192.70173300d+00
      a(2,2) =    3.22421598d+00
      c(3,2) =   11.35468650d+00
      a(3,2) =   16.15618810d+00
      c(1,3) =   10.84274330d+00
      a(1,3) =    9.92608326d+00
      c(2,3) =  116.52006000d+00
      a(2,3) =    2.86230717d+00
      c(3,3) =  -55.66109800d+00
      a(3,3) =    2.47244917d+00
      c(1,4) =   33.34249060d+00
      a(1,4) =    2.49115511d+00
      c(2,4) =    7.73004915d+00
      a(2,4) =    6.35100911d+00
      return
c
c   erbium
c
 68   continue
      c(1,1) =  -16.00435820d+00
      a(1,1) =   15.95965310d+00
      c(2,1) =   -6.26509474d+00
      a(2,1) =    3.08505824d+00
      c(1,2) = -150.36700200d+00
      a(1,2) =    3.03618526d+00
      c(2,2) =  236.19775600d+00
      a(2,2) =    3.35253476d+00
      c(3,2) =   11.41646110d+00
      a(3,2) =   16.03924000d+00
      c(1,3) =   11.89253580d+00
      a(1,3) =   11.04634630d+00
      c(2,3) =  130.54350100d+00
      a(2,3) =    3.01799114d+00
      c(3,3) =  -65.69229880d+00
      a(3,3) =    2.64727181d+00
      c(1,4) =   33.73180760d+00
      a(1,4) =    2.57008755d+00
      c(2,4) =    7.48868913d+00
      a(2,4) =    6.16564070d+00
      return
c
c   thulium
c
 69   continue
      c(1,1) =  -16.05536170d+00
      a(1,1) =   16.72435540d+00
      c(2,1) =   -6.29832935d+00
      a(2,1) =    3.22216196d+00
      c(1,2) = -143.24595300d+00
      a(1,2) =    3.18338303d+00
      c(2,2) =  233.19553800d+00
      a(2,2) =    3.53531467d+00
      c(3,2) =   11.61094440d+00
      a(3,2) =   16.16013920d+00
      c(1,3) =    7.39619033d+00
      a(1,3) =    2.46378238d+00
      c(2,3) =  119.57242900d+00
      a(2,3) =    2.72194361d+00
      c(3,3) =  -84.01744530d+00
      a(3,3) =    2.49552876d+00
      c(1,4) =   33.84116700d+00
      a(1,4) =    2.64361636d+00
      c(2,4) =    7.23350094d+00
      a(2,4) =    5.89426345d+00
      return
c
c   ytterbium
c
 70   continue
      c(1,1) =  -16.10497510d+00
      a(1,1) =   17.51231900d+00
      c(2,1) =   -6.33277963d+00
      a(2,1) =    3.36313057d+00
      c(1,2) = -106.50639900d+00
      a(1,2) =    3.29988442d+00
      c(2,2) =  201.14955200d+00
      a(2,2) =    3.77177093d+00
      c(3,2) =   12.01001670d+00
      a(3,2) =   16.63449170d+00
      c(1,3) =    7.59659144d+00
      a(1,3) =    2.65544200d+00
      c(2,3) =  119.64934700d+00
      a(2,3) =    2.83854581d+00
      c(3,3) =  -82.70469640d+00
      a(3,3) =    2.59880232d+00
      c(1,4) =   34.17132690d+00
      a(1,4) =    2.71026496d+00
      c(2,4) =    7.27486233d+00
      a(2,4) =    6.26424186d+00
      return
c
c   lutetium
c
 71   continue
      c(1,1) =  -16.15200380d+00
      a(1,1) =   18.32043650d+00
      c(2,1) =   -6.36766816d+00
      a(2,1) =    3.50797854d+00
      c(1,2) =  -71.31717050d+00
      a(1,2) =    3.34224801d+00
      c(2,2) =  169.40139300d+00
      a(2,2) =    4.00509413d+00
      c(3,2) =   11.88611540d+00
      a(3,2) =   16.18095210d+00
      c(1,3) =    7.71934058d+00
      a(1,3) =    2.72110668d+00
      c(2,3) =  116.10736300d+00
      a(2,3) =    2.98200454d+00
      c(3,3) =  -78.17862840d+00
      a(3,3) =    2.71480084d+00
      c(1,4) =   33.91192520d+00
      a(1,4) =    2.77571168d+00
      c(2,4) =    7.01871005d+00
      a(2,4) =    5.89380126d+00
c
      return
 6010 format (1x,'*** data not found for ecp type ',a4)
      end
**==sbkpm1.f
      subroutine sbkpm1 (zinp,lmx,nelcor,nt,n,c,a)
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
      implicit character * 8 (z)
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
      dimension nt(6), n(10,6), c(10,6), a(10,6)
c
      dimension ztit(11)
      data maxtyp/11/
      data ztit/'sc','ti','v ','cr','mn','fe',
     *          'co','ni','cu','zn','ga'/
c
c     this routine contains a library of ecp data taken from
c     various literature sources.
c     1st transition series sjkb potentials
c     
      do 10 ityp = 1 , maxtyp
         if (zinp.eq.ztit(ityp)) go to 20
 10   continue
      write (iwr,6010) zinp
      call caserr2('no library available for requested ecp')
 20   nelcor = 10
      lmx = 2
      nt(1) = 1
      nt(2) = 3
      nt(3) = 2
      n(1,1) = 1
      n(1,2) = 0
      n(2,2) = 2
      n(3,2) = 2
      n(1,3) = 0
      n(2,3) = 2
c
      go to (21,22,23,24,25,26,27,28,29,30,31), ityp
c
c     scandium (Sc_SBKJC_VDZ_ECP)
c
  21  continue
      c(1,1)=   -3.66956d+00
      a(1,1)=    9.95712d+00
      c(1,2)=    3.71885d+00
      a(1,2)=    1.46409d+00
      c(2,2)=  178.26066d+00
      a(2,2)=    3.59217d+00
      c(3,2)= -161.78394d+00
      a(3,2)=    3.32583d+00
      c(1,3)=    3.93048d+00
      a(1,3)=   23.53849d+00
      c(2,3)=   45.69140d+00
      a(2,3)=    6.58415d+00
      return
c
c     titanium (Ti_SBKJC_VDZ_ECP)
c
  22  continue
      c(1,1)=   -3.66581d+00
      a(1,1)=   11.20971d+00
      c(1,2)=    3.88067d+00
      a(1,2)=    1.61172d+00
      c(2,2)=  172.25051d+00
      a(2,2)=    4.00960d+00
      c(3,2)= -155.08356d+00
      a(3,2)=    3.66406d+00
      c(1,3)=    4.044327d+00
      a(1,3)=   31.67782d+00
      c(2,3)=   53.00048d+00
      a(2,3)=    7.45793d+00
      return
c
c     vanadium (V_SBKJC_VDZ_ECP)
c
  23  continue
      c(1,1)=   -3.79146d+00
      a(1,1)=   12.79659d+00
      c(1,2)=    3.43258d+00
      a(1,2)=    1.50439d+00
      c(2,2)=  178.03957d+00
      a(2,2)=    4.27321d+00
      c(3,2)= -156.30605d+00
      a(3,2)=    3.84121d+00
      c(1,3)=    3.764807d+00
      a(1,3)=   30.83720d+00
      c(2,3)=   58.96503d+00
      a(2,3)=    8.34454d+00
      return
c
c     chromium (Cr_SBKJC_VDZ_ECP)
c
  24  continue
      c(1,1)=   -3.64362d+00
      a(1,1)=   13.93307d+00
      c(1,2)=    3.59679d+00
      a(1,2)=    1.72136d+00
      c(2,2)=  184.07663d+00
      a(2,2)=    4.75681d+00
      c(3,2)= -161.41996d+00
      a(3,2)=    4.25348d+00
      c(1,3)=    3.68862d+00
      a(1,3)=   32.28223d+00
      c(2,3)=   65.48046d+00
      a(2,3)=    9.28137d+00
      return
c
c     manganese (Mn_SBKJC_VDZ_ECP)
c
  25  continue
      c(1,1)=   -3.85254d+00
      a(1,1)=   15.88592d+00
      c(1,2)=    3.71264d+00
      a(1,2)=    2.05981d+00
      c(2,2)=  181.67442d+00
      a(2,2)=    5.43248d+00
      c(3,2)= -156.33704d+00
      a(3,2)=    4.79887d+00
      c(1,3)=    4.368037d+00
      a(1,3)=   49.76765d+00
      c(2,3)=   75.53242d+00
      a(2,3)=   10.36623d+00
      return
c
c     iron (Fe_SBKJC_VDZ_ECP)
c
  26  continue
      c(1,1)=   -3.89423d+00
      a(1,1)=   17.61910d+00
      c(1,2)=    3.81706d+00
      a(1,2)=    2.33150d+00
      c(2,2)=  172.05349d+00
      a(2,2)=    6.08365d+00
      c(3,2)= -144.70056d+00
      a(3,2)=    5.26659d+00
      c(1,3)=    4.13737d+00
      a(1,3)=   49.40456d+00
      c(2,3)=   82.73696d+00
      a(2,3)=   11.40183d+00
      return
c
c     cobalt (Co_SBKJC_VDZ_ECP)
c
  27  continue
      c(1,1)=   -3.94746d+00
      a(1,1)=   19.44559d+00
      c(1,2)=    3.07705d+00
      a(1,2)=    2.25935d+00
      c(2,2)=  215.48382d+00
      a(2,2)=    6.66665d+00
      c(3,2)= -174.11816d+00
      a(3,2)=    5.76818d+00
      c(1,3)=    3.95321d+00
      a(1,3)=   50.28146d+00
      c(2,3)=   90.90043d+00
      a(2,3)=   12.52610d+00
      return
c
c     nickel (Ni_SBKJC_VDZ_ECP)
c
  28  continue
      c(1,1)=  -3.99337d+00
      a(1,1)=  21.36373d+00
      c(1,2)=   3.42572d+00
      a(1,2)=   2.21817d+00
      c(2,2)=  218.54146d+00
      a(2,2)=    6.71585d+00
      c(3,2)= -184.10961d+00
      a(3,2)=    5.80856d+00
      c(1,3)=    3.87540d+00
      a(1,3)=   56.45353d+00
      c(2,3)=   99.68289d+00
      a(2,3)=   13.66534d+00
      return
c
c     copper (Cu_SBKJC_VDZ_ECP)
c
  29  continue
      c(1,1)=   -4.00276d+00
      a(1,1)=   23.29060d+00
      c(1,2)=    3.32465d+00
      a(1,2)=    2.62059d+00
      c(2,2)=  224.32330d+00
      a(2,2)=    7.81082d+00
      c(3,2)= -180.30470d+00
      a(3,2)=    6.63324d+00
      c(1,3)=    3.42130d+00
      a(1,3)=   51.15734d+00
      c(2,3)=  105.70126d+00
      a(2,3)=   14.73406d+00
      return
c
c     zinc (Zn_SBKJC_VDZ_ECP)
c
  30  continue
      c(1,1)=   -3.91648d+00
      a(1,1)=   24.88705d+00
      c(1,2)=    3.43082d+00
      a(1,2)=    2.79159d+00
      c(2,2)=  244.31634d+00
      a(2,2)=    8.32676d+00
      c(3,2)= -199.42711d+00
      a(3,2)=    7.09651d+00
      c(1,3)=    4.11473d+00
      a(1,3)=   75.70603d+00
      c(2,3)=  118.29936d+00
      a(2,3)=   16.07627d+00
      return
c
c     gallium (Ga_SBKJC_VDZ_ECP)
c
  31  continue
      c(1,1)=   -3.87363d+00
      a(1,1)=   26.74302d+00
      c(1,2)=    4.12472d+00
      a(1,2)=    3.46530d+00
      c(2,2)=  260.73263d+00
      a(2,2)=    9.11130d+00
      c(3,2)= -223.96003d+00
      a(3,2)=    7.89329d+00
      c(1,3)=    4.20033d+00
      a(1,3)=   79.99353d+00
      c(2,3)=  127.99139d+00
      a(2,3)=   17.39114d+00
      return
 6010 format (1x,'*** data not found for ecp type ',a4)
      end
**==sbkpm2.f
      subroutine sbkpm2 (zinp,lmx,nelcor,nt,n,c,a)
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
      implicit character * 8 (z)
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
      dimension nt(6), n(10,6), c(10,6), a(10,6)
c
      dimension ztit(11)
      data maxtyp/11/
      data ztit/'y ','zr','nb','mo','tc','ru',
     *          'rh','pd','ag','cd','in'/
c
c     this routine contains a library of ecp data taken from
c     various literature sources.
c     2nd transition series sjkb potentials
c     
      do 10 ityp = 1 , maxtyp
         if (zinp.eq.ztit(ityp)) go to 20
 10   continue
      write (iwr,6010) zinp
      call caserr2('no library available for requested ecp')
 20   nelcor = 28
      lmx = 3
      nt(1) = 1
      nt(2) = 3
      nt(3) = 3
      nt(4) = 2
      n(1,1) = 1
      n(1,2) = 0
      n(2,2) = 2
      n(3,2) = 2
      n(1,3) = 0
      n(2,3) = 2
      n(3,3) = 2
      n(1,4) = 0
      n(2,4) = 2
c
      go to (39,40,41,42,43,44,45,46,47,48,49), ityp
c
c      yttrium (Y_SBKJC_VDZ_ECP)
c
  39  continue
      c(1,1)=   -3.70827d+00
      a(1,1)=    3.13323d+00
      c(1,2)=    3.52458d+00
      a(1,2)=    0.74734d+00
      c(2,2)= -115.67795d+00
      a(2,2)=    2.11046d+00
      c(3,2)=  141.66551d+00
      a(3,2)=    2.37431d+00
      c(1,3)=    3.45319d+00
      a(1,3)=    0.82924d+00
      c(2,3)=  -67.46192d+00
      a(2,3)=    1.92973d+00
      c(3,3)=   85.09453d+00
      a(3,3)=    2.20379d+00
      c(1,4)=    2.29590d+00
      a(1,4)=    9.36770d+00
      c(2,4)=   20.62182d+00
      a(2,4)=    2.83038d+00
      return
c
c     zirconium (Zr_SBKJC_VDZ_ECP)
c
  40  continue
      c(1,1)=   -6.11923d+00
      a(1,1)=    4.31441d+00
      c(1,2)=    4.57031d+00
      a(1,2)=    0.93035d+00
      c(2,2)= -112.55854d+00
      a(2,2)=    2.33003d+00
      c(3,2)=  139.02972d+00
      a(3,2)=    2.66632d+00
      c(1,3)=    4.04524d+00
      a(1,3)=    0.97753d+00
      c(2,3)=  -66.64768d+00
      a(2,3)=    2.14998d+00
      c(3,3)=   86.37034d+00
      a(3,3)=    2.50892d+00
      c(1,4)=    5.64085d+00
      a(1,4)=   16.57598d+00
      c(2,4)=   26.66093d+00
      a(2,4)=    3.32079d+00
      return
c
c     niobium (Nb_SBKJC_VDZ_ECP)
c
  41  continue
      c(1,1)=   -5.81695d+00
      a(1,1)=    4.64198d+00
      c(1,2)=    4.99937d+00
      a(1,2)=    1.07339d+00
      c(2,2)= -100.17856d+00
      a(2,2)=    2.51174d+00
      c(3,2)=  127.21570d+00
      a(3,2)=    2.93940d+00
      c(1,3)=    3.64794d+00
      a(1,3)=    0.99626d+00
      c(2,3)=  -62.17437d+00
      a(2,3)=    2.28020d+00
      c(3,3)=   85.80845d+00
      a(3,3)=    2.73360d+00
      c(1,4)=    2.77206d+00
      a(1,4)=   10.48803d+00
      c(2,4)=   29.05246d+00
      a(2,4)=    3.64266d+00
      return
c
c     molybdenum (Mo_SBKJC_VDZ_ECP)
c
  42  continue
      c(1,1)=  -7.23009d+00
      a(1,1)=   5.62507d+00
      c(1,2)=   3.83695d+00
      a(1,2)=   1.01907d+00
      c(2,2)= -90.19151d+00
      a(2,2)=   2.68386d+00
      c(3,2)= 129.35214d+00
      a(3,2)=   3.27945d+00
      c(1,3)=   3.83295d+00
      a(1,3)=   1.10957d+00
      c(2,3)= -62.98231d+00
      a(2,3)=   2.49708d+00
      c(3,3)=  89.77478d+00
      a(3,3)=   3.04740d+00
      c(1,4)=   2.73486d+00
      a(1,4)=   9.70400d+00
      c(2,4)=  32.53231d+00
      a(2,4)=   4.03098d+00
      return
c
c     technetium (Tc_SBKJC_VDZ_ECP)
c
  43  continue
      c(1,1)=  -7.54356d+00
      a(1,1)=   6.26067d+00
      c(1,2)=   4.31377d+00
      a(1,2)=   1.17731d+00
      c(2,2)= -80.61847d+00
      a(2,2)=   2.86418d+00
      c(3,2)= 120.82547d+00
      a(3,2)=   3.61432d+00
      c(1,3)=   3.38952d+00
      a(1,3)=   1.08449d+00
      c(2,3)= -57.06568d+00
      a(2,3)=   2.59907d+00
      c(3,3)=  89.03016d+00
      a(3,3)=   3.30452d+00
      c(1,4)=   2.86365d+00
      a(1,4)=  11.62092d+00
      c(2,4)=  37.18014d+00
      a(2,4)=   4.45900d+00
      return
c
c     ruthenium (Ru_SBKJC_VDZ_ECP)
c
  44  continue
      c(1,1)=  -7.98448d+00
      a(1,1)=   6.98304d+00
      c(1,2)=   4.82384d+00
      a(1,2)=   1.29393d+00
      c(2,2)= -81.31936d+00
      a(2,2)=   3.04871d+00
      c(3,2)= 124.36917d+00
      a(3,2)=   3.87163d+00
      c(1,3)=   3.86590d+00
      a(1,3)=   1.20273d+00
      c(2,3)= -72.15048d+00
      a(2,3)=   2.82303d+00
      c(3,3)= 103.67006d+00
      a(3,3)=   3.48391d+00
      c(1,4)=   3.21083d+00
      a(1,4)=  14.33743d+00
      c(2,4)=  42.07321d+00
      a(2,4)=   4.90003d+00
      return
c
c     rhodium (Rh_SBKJC_VDZ_ECP)
c
  45  continue
      c(1,1)=   -8.13225d+00
      a(1,1)=    7.64968d+00
      c(1,2)=    5.20773d+00
      a(1,2)=    1.33777d+00
      c(2,2)=  -50.06846d+00
      a(2,2)=    2.95308d+00
      c(3,2)=   90.36859d+00
      a(3,2)=    4.39131d+00
      c(1,3)=    3.85874d+00
      a(1,3)=    1.19468d+00
      c(2,3)=  -72.74348d+00
      a(2,3)=    2.92434d+00
      c(3,3)=  106.11615d+00
      a(3,3)=    3.67503d+00
      c(1,4)=    3.34825d+00
      a(1,4)=   17.89248d+00
      c(2,4)=   47.13883d+00
      a(2,4)=    5.34874d+00
      return
c
c     palladium (Pd_SBKJC_VDZ_ECP)
c
  46  continue
      c(1,1)=   -8.33466d+00
      a(1,1)=    8.34327d+00
      c(1,2)=    3.97809d+00
      a(1,2)=    1.21315d+00
      c(2,2)=  218.99250d+00
      a(2,2)=    4.03301d+00
      c(3,2)= -168.44404d+00
      a(3,2)=    3.45772d+00
      c(1,3)=    4.14021d+00
      a(1,3)=    1.45607d+00
      c(2,3)=  108.32601d+00
      a(2,3)=    4.05510d+00
      c(3,3)=  -72.83147d+00
      a(3,3)=   3.21616d+00
      c(1,4)=   3.30611d+00
      a(1,4)=   1.954477d+01
      c(2,4)=  52.78802d+00
      a(2,4)=   5.84664d+00
      return
c
c     silver (Ag_SBKJC_VDZ_ECP)
c
  47  continue
      c(1,1)=   -8.44669d+00
      a(1,1)=    9.01696d+00
      c(1,2)=    6.20634d+00
      a(1,2)=    1.51740d+00
      c(2,2)=  -56.82624d+00
      a(2,2)=    3.23088d+00
      c(3,2)=   92.40173d+00
      a(3,2)=    4.82747d+00
      c(1,3)=    3.75637d+00
      a(1,3)=    1.21699d+00
      c(2,3)=  -80.95335d+00
      a(2,3)=    3.18615d+00
      c(3,3)=  118.86801d+00
      a(3,3)=    4.05349d+00
      c(1,4)=    3.40895d+00
      a(1,4)=   22.94020d+00
      c(2,4)=   57.10513d+00
      a(2,4)=    6.26627d+00
      return
c
c     cadmium (Cd_SBKJC_VDZ_ECP)
c
  48  continue
      c(1,1)=   -8.69864d+00
      a(1,1)=    9.81325d+00
      c(1,2)=    6.56934d+00
      a(1,2)=    1.60476d+00
      c(2,2)=  -59.36395d+00
      a(2,2)=    3.37727d+00
      c(3,2)=   93.59610d+00
      a(3,2)=    5.08876d+00
      c(1,3)=    3.67863d+00
      a(1,3)=    1.25171d+00
      c(2,3)= -104.50694d+00
      a(2,3)=    3.42592d+00
      c(3,3)=  145.67316d+00
      a(3,3)=    4.21415d+00
      c(1,4)=    3.33094d+00
      a(1,4)=   22.35337d+00
      c(2,4)=   62.06616d+00
      a(2,4)=    6.75408d+00
      return
c
c     indium (In_SBKJC_VDZ_ECP)
c
  49  continue
      c(1,1)=   -8.84031d+00
      a(1,1)=   10.58987d+00
      c(1,2)=    6.90124d+00
      a(1,2)=    1.73883d+00
      c(2,2)=  -68.48319d+00
      a(2,2)=    3.62806d+00
      c(3,2)=  101.77272d+00
      a(3,2)=    5.30270d+00
      c(1,3)=    3.98554d+00
      a(1,3)=    1.43045d+00
      c(2,3)= -129.78594d+00
      a(2,3)=    3.74447d+00
      c(3,3)=  171.90487d+00
      a(3,3)=    4.47679d+00
      c(1,4)=    3.49244d+00
      a(1,4)=   29.16072d+00
      c(2,4)=   68.27816d+00
      a(2,4)=    7.26566d+00
      return
 6010 format (1x,'*** data not found for ecp type ',a4)
      end
**==sbkpm3.f
      subroutine sbkpm3 (zinp,lmx,nelcor,nt,n,c,a)
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
      implicit character * 8 (z)
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
      dimension nt(6), n(10,6), c(10,6), a(10,6)
c
      dimension ztit(11)
      data maxtyp/11/
      data ztit/
     $ 'la','hf',
     $ 'ta','w ','re','os','ir','pt','au','hg','tl'/
c
c     this routine contains a library of ecp data taken from
c     various literature sources.
c     3rd transition series sjkb potentials
c     
      do 10 ityp = 1 , maxtyp
         if (zinp.eq.ztit(ityp)) go to 20
 10   continue
      write (iwr,6010) zinp
      call caserr2('no library available for requested ecp')
c
 20   go to(57,72,73,74,75,76,77,78,79,80,81), ityp
c
c     lanthanum (La_SBKJC_VDZ_ECP)
c
  57  continue
      nelcor = 46
      lmx = 4
      nt(1) = 1
      nt(2) = 3
      nt(3) = 3
      nt(4) = 2
      nt(5) = 2
      n(1,1) = 1
      n(1,2) = 0
      n(2,2) = 2
      n(3,2) = 2
      n(1,3) = 0
      n(2,3) = 2
      n(3,3) = 2
      n(1,4) = 0
      n(2,4) = 2
      n(1,5) = 1
      n(2,5) = 1
      c(1,1)=  -2.85665d+00
      a(1,1)=   1.42345d+00
      c(1,2)=   8.85347d+00
      a(1,2)=  22.29739d+00
      c(2,2)= -68.41778d+00
      a(2,2)=   1.64320d+00
      c(3,2)= 110.12463d+00
      a(3,2)=   1.85122d+00
      c(1,3)=   2.98675d+00
      a(1,3)=   1.77635d+00
      c(2,3)= -35.51815d+00
      a(2,3)=   3.98941d+00
      c(3,3)=  53.65423d+00
      a(3,3)=   2.47524d+00
      c(1,4)=   1.67298d+00
      a(1,4)=   1.52169d+00
      c(2,4)=  18.46625d+00
      a(2,4)=   1.52182d+00
      c(1,5)=  -4.64430d+00
      a(1,5)=   2.67762d+00
      c(2,5)= -14.43101d+00
      a(2,5)=  10.06608d+00
      return
c
c     hafnium (Hf_SBKJC_VDZ_ECP)
c
  72  continue
      c(1,1)=   -2.59020d+00
      a(1,1)=    1.14090d+00
      c(1,2)=    5.08602d+00
      a(1,2)=   18.98210d+00
      c(2,2)=  -76.34646d+00
      a(2,2)=    3.92413d+00
      c(3,2)=  154.12616d+00
      a(3,2)=    4.32002d+00
      c(1,3)=    3.80263d+00
      a(1,3)=    1.11008d+00
      c(2,3)= -133.12772d+00
      a(2,3)=    2.38920d+00
      c(3,3)=  159.39488d+00
      a(3,3)=    2.58171d+00
      c(1,4)=    5.23486d+00
      a(1,4)=   21.0000d+00
      c(2,4)=   30.09851d+00
      a(2,4)=    2.80487d+00
      c(1,5)=   11.03700d+00
      a(1,5)=    1.11377d+00
      go to 30
c
c     tantalum (Ta_SBKJC_VDZ_ECP)
c
  73  continue
      c(1,1)=   -4.33440d+00
      a(1,1)=    1.65262d+00
      c(1,2)=    4.24022d+00
      a(1,2)=    3.78450d+00
      c(2,2)=  -30.38231d+00
      a(2,2)=    3.17838d+00
      c(3,2)=   89.32811d+00
      a(3,2)=    3.95058d+00
      c(1,3)=    3.33868d+00
      a(1,3)=    1.00580d+00
      c(2,3)= -138.02173d+00
      a(2,3)=    2.35390d+00
      c(3,3)=  171.53349d+00
      a(3,3)=    2.58191d+00
      c(1,4)=    3.04913d+00
      a(1,4)=   24.29129d+00
      c(2,4)=   41.20515d+00
      a(2,4)=    3.18604d+00
      c(1,5)=    7.84634d+00
      a(1,5)=    1.13011d+00
      go to 30
c
c     tungsten (W_SBKJC_VDZ_ECP)
c
  74  continue
      c(1,1)=   -5.26480d+00
      a(1,1)=    2.03828d+00
      c(1,2)=    2.47299d+00
      a(1,2)=    0.96514d+00
      c(2,2)=  -52.07178d+00
      a(2,2)=    2.65502d+00
      c(3,2)=  113.46745d+00
      a(3,2)=    3.62713d+00
      c(1,3)=    3.91085d+00
      a(1,3)=    1.23902d+00
      c(2,3)= -135.31607d+00
      a(2,3)=    2.54122d+00
      c(3,3)=  172.36976d+00
      a(3,3)=    2.80797d+00
      c(1,4)=    3.87349d+00
      a(1,4)=   27.55935d+00
      c(2,4)=   49.52511d+00
      a(2,4)=    3.51973d+00
      c(1,5)=    6.99471d+00
      a(1,5)=    1.20561d+00
      go to 30
c
c     rhenium (Re_SBKJC_VDZ_ECP)
c
  75  continue
      c(1,1)=   -6.62225d+00
      a(1,1)=    2.51731d+00
      c(1,2)=    3.32067d+00
      a(1,2)=    1.00603d+00
      c(2,2)=  -39.91860d+00
      a(2,2)=    2.52065d+00
      c(3,2)=  101.83236d+00
      a(3,2)=    3.86550d+00
      c(1,3)=    3.83735d+00
      a(1,3)=    1.22884d+00
      c(2,3)= -121.99416d+00
      a(2,3)=    2.65033d+00
      c(3,3)=  166.17144d+00
      a(3,3)=    3.01686d+00
      c(1,4)=    4.57819d+00
      a(1,4)=   29.23499d+00
      c(2,4)=   58.72776d+00
      a(2,4)=    3.86422d+00
      c(1,5)=    6.44315d+00
      a(1,5)=    1.27456d+00
      go to 30
c
c     osmium (Os_SBKJC_VDZ_ECP)
c
  76  continue
      c(1,1)=   -7.96753d+00
      a(1,1)=    3.01084d+00
      c(1,2)=    4.03478d+00
      a(1,2)=    1.17006d+00
      c(2,2)=  -84.75261d+00
      a(2,2)=    2.93435d+00
      c(3,2)=  148.94410d+00
      a(3,2)=    3.82663d+00
      c(1,3)=    3.10454d+00
      a(1,3)=    1.02178d+00
      c(2,3)= -120.70117d+00
      a(2,3)=    2.66013d+00
      c(3,3)=  172.77128d+00
      a(3,3)=    3.10435d+00
      c(1,4)=    6.63507d+00
      a(1,4)=   31.00000d+00
      c(2,4)=   66.85726d+00
      a(2,4)=    4.17953d+00
      c(1,5)=    6.33120d+00
      a(1,5)=    1.38742d+00
      go to 30
c
c     iridium (Ir_SBKJC_VDZ_ECP)
c
  77  continue
      c(1,1)=   -8.98208d+00
      a(1,1)=    3.46859d+00
      c(1,2)=    4.83629d+00
      a(1,2)=    1.27403d+00
      c(2,2)=  -98.47624d+00
      a(2,2)=    3.07970d+00
      c(3,2)=  162.58173d+00
      a(3,2)=    3.97153d+00
      c(1,3)=    3.80021d+00
      a(1,3)=    1.08342d+00
      c(2,3)= -151.65910d+00
      a(2,3)=    2.79432d+00
      c(3,3)=  203.96216d+00
      a(3,3)=    3.21008d+00
      c(1,4)=    7.96874d+00
      a(1,4)=   33.00000d+00
      c(2,4)=   73.97844d+00
      a(2,4)=    4.47188d+00
      c(1,5)=    6.28559d+00
      a(1,5)=    1.48524d+00
      go to 30
c
c     platium (Pt_SBKJC_VDZ_ECP)
c
  78  continue
      c(1,1)=   -9.79854d+00
      a(1,1)=    3.90813d+00
      c(1,2)=    5.79843d+00
      a(1,2)=    1.36279d+00
      c(2,2)= -187.51594d+00
      a(2,2)=    3.34227d+00
      c(3,2)=  247.96906d+00
      a(3,2)=    3.91136d+00
      c(1,3)=    3.67055d+00
      a(1,3)=    1.07586d+00
      c(2,3)=  -93.02055d+00
      a(2,3)=    2.80443d+00
      c(3,3)=  151.55915d+00
      a(3,3)=    3.52491d+00
      c(1,4)=   13.09353d+00
      a(1,4)=   42.63116d+00
      c(2,4)=   80.29385d+00
      a(2,4)=    4.744683d+00
      c(1,5)=    6.13870d+00
      a(1,5)=    1.58988d+00
      go to 30
c
c     gold (Au_SBKJC_VDZ_ECP)
c
  79  continue
      c(1,1)=  -10.72358d+00
      a(1,1)=    4.38763d+00
      c(1,2)=    6.35612d+00
      a(1,2)=    1.55636d+00
      c(2,2)= -364.44039d+00
      a(2,2)=    3.71593d+00
      c(3,2)=  428.19753d+00
      a(3,2)=    4.06792d+00
      c(1,3)=    4.41518d+00
      a(1,3)=    1.18798d+00
      c(2,3)= -136.57550d+00
      a(2,3)=    3.01551d+00
      c(3,3)=  194.20535d+00
      a(3,3)=    3.59588d+00
      c(1,4)=    8.88198d+00
      a(1,4)=   35.25000d+00
      c(2,4)=   86.76612d+00
      a(2,4)=    5.02307d+00
      c(1,5)=    6.21603d+00
      a(1,5)=    1.68881d+00
      go to 30
c
c     mercury (Hg_SBKJC_VDZ_ECP)
c
  80  continue
      c(1,1)=  -11.51618d+00
      a(1,1)=    4.85764d+00
      c(1,2)=    4.09547d+00
      a(1,2)=    1.17357d+00
      c(2,2)= -310.08898d+00
      a(2,2)=    3.70966d+00
      c(3,2)=  394.68735d+00
      a(3,2)=    4.18793d+00
      c(1,3)=    3.61202d+00
      a(1,3)=    0.95913d+00
      c(2,3)= -148.42850d+00
      a(2,3)=    2.99268d+00
      c(3,3)=  211.83161d+00
      a(3,3)=    3.58353d+00
      c(1,4)=    8.27074d+00
      a(1,4)=   35.54665d+00
      c(2,4)=   92.85678d+00
      a(2,4)=    5.29318d+00
      c(1,5)=    6.28302d+00
      a(1,5)=    1.80978d+00
      go to 30
c
c     thallium (Tl_SBKJC_VDZ_ECP)
c
  81  continue
      c(1,1)=  -12.13787d+00
      a(1,1)=    5.31294d+00
      c(1,2)=    4.77833d+00
      a(1,2)=    1.29566d+00
      c(2,2)= -253.19592d+00
      a(2,2)=    3.80665d+00
      c(3,2)=  337.05349d+00
      a(3,2)=    4.43233d+00
      c(1,3)=    3.45912d+00
      a(1,3)=    0.94782d+00
      c(2,3)= -168.94091d+00
      a(2,3)=    3.14378d+00
      c(3,3)=  238.10901d+00
      a(3,3)=    3.73562d+00
      c(1,4)=    7.41137d+00
      a(1,4)=   33.43080d+00
      c(2,4)=   98.49883d+00
      a(2,4)=    5.55174d+00
      c(1,5)=    6.27606d+00
      a(1,5)=    1.91385d+00
c
 30   nelcor = 60
      lmx = 4
      nt(1) = 1
      nt(2) = 3
      nt(3) = 3
      nt(4) = 2
      nt(5) = 1
      n(1,1) = 1
      n(1,2) = 0
      n(2,2) = 2
      n(3,2) = 2
      n(1,3) = 0
      n(2,3) = 2
      n(3,3) = 2
      n(1,4) = 0
      n(2,4) = 2
      n(1,5) = 0
      return
 6010 format (1x,'*** data not found for ecp type ',a4)
      end
**==crenbl0.f
      subroutine crenbl0(zinp,lmx,nelcor,nt,n,cc,a)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character * 8 (z)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      dimension nt(6), n(10,6), cc(10,6), a(10,6)
c
      dimension ztit(18)
      data maxtyp/18/
      data ztit/'h','he',
     *          'li','be','b ','c ','n ','o ','f ','ne',
     *          'na','mg','al','si','p ','s ','cl','ar'/
c
c     this routine contains a library of ecp data taken from
c     various literature sources.
c     the CRENBL ecps are due to Christiansen etc. al
c     Large Orbital Basis, Small Core Potential
c     shape consistent ECPs
c     L.F. Pacios and P.A> Christiansen J.Chem. Phys. 82 (1985) 2664
c
      do ityp = 1 , maxtyp
         if (zinp.eq.ztit(ityp)) go to 40
      enddo
 30   write (iwr,6010) zinp
      call caserr2('no library available for requested CRENBL ecp')
c
 40   go to (30,30,50,60,70,80,90,100,110,120,
     +      130,150,170,190,210,230,250,270) , ityp
      go to 30
c
c     lithium  Christiansen et al. (Li_CRENBL_ECP)
c
 50   continue
      cc( 1, 1) =         -0.11866700d0
      cc( 2, 1) =         -1.21779400d0
      cc( 3, 1) =         -1.37580501d0
      a ( 1, 1) =          0.80780000d0
      a ( 2, 1) =          2.55000001d0
      a ( 3, 1) =          7.25349998d0
      cc( 1, 2) =         24.33369899d0
      cc( 2, 2) =        -20.66639709d0
      cc( 3, 2) =         -1.14289200d0
      cc( 4, 2) =          2.99401900d0
      a ( 1, 2) =          0.79980000d0
      a ( 2, 2) =          0.77410000d0
      a ( 3, 2) =          1.19430000d0
      a ( 4, 2) =          2.19990000d0
c
      go to 65
c
c     beryllium  Christiansen et al. (Be_CRENBL_ECP)
c
 60   continue
      cc( 1, 1) =         -0.25828800d0
      cc( 2, 1) =         -2.09550900d0
      cc( 3, 1) =         -1.39610200d0
      a ( 1, 1) =          1.88230000d0
      a ( 2, 1) =          6.05660000d0
      a ( 3, 1) =         17.07800000d0
      cc( 1, 2) =        -28.57470500d0
      cc( 2, 2) =         35.73522700d0
      cc( 3, 2) =         -1.31942100d0
      cc( 4, 2) =          2.96972500d0
      a ( 1, 2) =          1.51470000d0
      a ( 2, 2) =          1.58470000d0
      a ( 3, 2) =          2.30910000d0
      a ( 4, 2) =          4.20010000d0
65    n ( 1, 1) =   2
      n ( 2, 1) =   2
      n ( 3, 1) =   1
      nt( 1)    =   3
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   1
      n ( 4, 2) =   0
      nt( 2)    =   4
      nelcor    =   2
      lmx       =   1
c
      go to 280
c
c     boron  Christiansen et al. (B_CRENBL_ECP)
c
 70   continue
      cc( 1, 1) =         -0.43494600d0
      cc( 2, 1) =         -3.20745206d0
      cc( 3, 1) =         -1.41778803d0
      a ( 1, 1) =          3.56130004d0
      a ( 2, 1) =         11.78320026d0
      a ( 3, 1) =         33.37820053d0
      cc( 1, 2) =        -37.55111694d0
      cc( 2, 2) =         51.02549744d0
      cc( 3, 2) =         -2.38443398d0
      cc( 4, 2) =          2.97974205d0
      a ( 1, 2) =          2.51550007d0
      a ( 2, 2) =          2.68079996d0
      a ( 3, 2) =          3.95619988d0
      a ( 4, 2) =          7.41090012d0
c
      go to 125
c
c     carbon Christiansen et al. (C_CRENBL_ECP)
c
 80   continue
      cc( 1, 1) =         -0.55931300d0
      cc( 2, 1) =         -4.07454997d0
      cc( 3, 1) =         -1.43484600d0
      a ( 1, 1) =          5.35280001d0
      a ( 2, 1) =         18.06680012d0
      a ( 3, 1) =         51.61590004d0
      cc( 1, 2) =        -47.09821510d0
      cc( 2, 2) =         71.58925819d0
      cc( 3, 2) =         -4.67536402d0
      cc( 4, 2) =          3.03797001d0
      a ( 1, 2) =          3.81909999d0
      a ( 2, 2) =          4.17320001d0
      a ( 3, 2) =          6.27069998d0
      a ( 4, 2) =         12.21120000d0
c
      go to 125
c
c     nitrogen Christiansen et al. (N_CRENBL_ECP)
c
 90   continue
      cc( 1, 1) =         -0.67093700d0
      cc( 2, 1) =         -4.87202299d0
      cc( 3, 1) =         -1.45602100d0
      a ( 1, 1) =          7.47939998d0
      a ( 2, 1) =         25.34080005d0
      a ( 3, 1) =         73.08880043d0
      cc( 1, 2) =        -60.38108397d0
      cc( 2, 2) =         69.86420631d0
      cc( 3, 2) =         -9.60874701d0
      cc( 4, 2) =          3.27333000d0
      a ( 1, 2) =          4.17259997d0
      a ( 2, 2) =          4.74250001d0
      a ( 3, 2) =          8.16100001d0
      a ( 4, 2) =          2.38600001d0
c
      go to 125
c
c     oxygen Christiansen et al. (O_CRENBL_ECP)
c
 100  continue
      cc( 1, 1) =         -0.79842000d0
      cc( 2, 1) =         -5.76684701d0
      cc( 3, 1) =         -1.48645601d0
      a ( 1, 1) =         10.02859998d0
      a ( 2, 1) =         34.19799995d0
      a ( 3, 1) =        100.00389957d0
      cc( 1, 2) =         11.21630394d0
      cc( 2, 2) =        -16.34447694d0
      cc( 3, 2) =          1.04294400d0
      cc( 4, 2) =          2.19389099d0
      a ( 1, 2) =          2.24790001d0
      a ( 2, 2) =          2.40490001d0
      a ( 3, 2) =          4.37400001d0
      a ( 4, 2) =          2.18920001d0
c
      go to 125
c
c     fluorine Christiansen et al. (F_CRENBL_ECP)
c
 110  continue
      cc( 1, 1) =         -6.72302401d0
      cc( 2, 1) =         -0.92964900d0
      cc( 3, 1) =         -1.52673399d0
      a ( 1, 1) =         44.51660013d0
      a ( 2, 1) =         12.94869995d0
      a ( 3, 1) =        132.49670029d0
      cc( 1, 2) =         12.68530595d0
      cc( 2, 2) =        -19.30258894d0
      cc( 3, 2) =          1.00217900d0
      cc( 4, 2) =          2.24534899d0
      a ( 1, 2) =          2.88350001d0
      a ( 2, 2) =          3.10769999d0
      a ( 3, 2) =          5.61220002d0
      a ( 4, 2) =          2.81459999d0
c
      go to 125
c
c     neon Christiansen et al. (Ne_CRENBL_ECP)
c
 120  continue
      cc( 1, 1) =         -1.05784200d0
      cc( 2, 1) =         -7.71812600d0
      cc( 3, 1) =         -1.57767500d0
      a ( 1, 1) =         16.20239997d0
      a ( 2, 1) =         56.17350006d0
      a ( 3, 1) =        170.70949936d0
      cc( 1, 2) =         13.66500998d0
      cc( 2, 2) =        -21.93867993d0
      cc( 3, 2) =          0.92151100d0
      cc( 4, 2) =          2.28758100d0
      a ( 1, 2) =          3.57620001d0
      a ( 2, 2) =          3.88209999d0
      a ( 3, 2) =          6.98530000d0
      a ( 4, 2) =          3.48379999d0
125   continue
      n ( 1, 1) =   2
      n ( 2, 1) =   2
      n ( 3, 1) =   1
      nt( 1)    =   3
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   1
      n ( 4, 2) =   0
      nt( 2)    =   4
      nelcor    =   2
      lmx       =   1
c
      go to 280
c
c     sodium
c
 130  continue
c     Christiansen et al. (Na_CRENBL_ECP)
      cc( 1, 1) =         -0.00010200d0
      cc( 2, 1) =         -1.71544300d0
      cc( 3, 1) =         -1.36191100d0
      a ( 1, 1) =          1.60680000d0
      a ( 2, 1) =         22.52309990d0
      a ( 3, 1) =         76.23649979d0
      cc( 1, 2) =        -14.53866100d0
      cc( 2, 2) =         31.91120791d0
      cc( 3, 2) =        -32.32224607d0
      cc( 4, 2) =          3.14094701d0
      cc( 5, 2) =          1.87765500d0
      a ( 1, 2) =          2.47830001d0
      a ( 2, 2) =          3.09900001d0
      a ( 3, 2) =          3.94710001d0
      a ( 4, 2) =          8.21659994d0
      a ( 5, 2) =          1.50080000d0
c
      go to 155
c
c     magnesium
c
 150  continue
c     Christiansen et al. (Mg_CRENBL_ECP)
      cc( 1, 1) =         -0.00005300d0
      cc( 2, 1) =         -1.93249400d0
      cc( 3, 1) =         -1.39384700d0
      a ( 1, 1) =          1.13200000d0
      a ( 2, 1) =         27.30570006d0
      a ( 3, 1) =         93.30560017d0
      cc( 1, 2) =        -15.62671006d0
      cc( 2, 2) =         33.53119421d0
      cc( 3, 2) =        -36.07307196d0
      cc( 4, 2) =          3.12963399d0
      cc( 5, 2) =          1.93843301d0
      a ( 1, 2) =          2.95159999d0
      a ( 2, 2) =          3.73249999d0
      a ( 3, 2) =          4.81180000d0
      a ( 4, 2) =          9.91949999d0
      a ( 5, 2) =          1.81850000d0
155   n ( 1, 1) =   2
      n ( 2, 1) =   2
      n ( 3, 1) =   1
      nt( 1)    =   3
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      n ( 4, 2) =   1
      n ( 5, 2) =   0
      nt( 2)    =   5
      nelcor    =   2
      lmx       =   1
c
      go to 280
c
c     aluminium
c
 170  continue
c     Christiansen et al. (Al_CRENBL_ECP)
      cc( 1, 1) =         -0.53798100d0
      cc( 2, 1) =         -5.45975602d0
      cc( 3, 1) =        -16.65534306d0
      cc( 4, 1) =         -6.47521502d0
      a ( 1, 1) =          1.22110000d0
      a ( 2, 1) =          3.36809999d0
      a ( 3, 1) =          9.75000000d0
      a ( 4, 1) =         29.26929998d0
      cc( 1, 2) =        -56.20521307d0
      cc( 2, 2) =        149.68995476d0
      cc( 3, 2) =        -91.45439434d0
      cc( 4, 2) =          3.72894901d0
      cc( 5, 2) =          3.03799400d0
      a ( 1, 2) =          1.56310000d0
      a ( 2, 2) =          1.77120000d0
      a ( 3, 2) =          2.06230000d0
      a ( 4, 2) =          3.35830000d0
      a ( 5, 2) =          2.13000000d0
      cc( 1, 3) =         93.67560577d0
      cc( 2, 3) =       -189.88896751d0
      cc( 3, 3) =        110.24810410d0
      cc( 4, 3) =          4.19959599d0
      cc( 5, 3) =          5.00335598d0
      a ( 1, 3) =          1.82310000d0
      a ( 2, 3) =          2.12490001d0
      a ( 3, 3) =          2.57049999d0
      a ( 4, 3) =          1.75749999d0
      a ( 5, 3) =          6.76929998d0
      go to 275
c
c     silicon
c
 190  continue
c     Christiansen et al. (Si_CRENBL_ECP)
      cc( 1, 1) =         -0.67587200d0
      cc( 2, 1) =         -6.72859800d0
      cc( 3, 1) =        -19.62386489d0
      cc( 4, 1) =         -6.58675301d0
      a ( 1, 1) =          1.62109999d0
      a ( 2, 1) =          4.53210002d0
      a ( 3, 1) =         13.20760000d0
      a ( 4, 1) =         39.96990013d0
      cc( 1, 2) =        -54.56442690d0
      cc( 2, 2) =        146.37274933d0
      cc( 3, 2) =        -88.33109474d0
      cc( 4, 2) =          6.76354098d0
      cc( 5, 2) =          3.11288401d0
      a ( 1, 2) =          1.86890000d0
      a ( 2, 2) =          2.11390001d0
      a ( 3, 2) =          2.45269999d0
      a ( 4, 2) =          3.88029999d0
      a ( 5, 2) =          2.56439999d0
      cc( 1, 3) =        116.56166363d0
      cc( 2, 3) =       -224.59841347d0
      cc( 3, 3) =        138.22863007d0
      cc( 4, 3) =          2.85250399d0
      cc( 5, 3) =          5.07075799d0
      a ( 1, 3) =          2.47999999d0
      a ( 2, 3) =          2.97729999d0
      a ( 3, 3) =          3.77399999d0
      a ( 4, 3) =          5.56639999d0
      a ( 5, 3) =         10.31330001d0
c
      go to 275
c
c     phosphorus
c
 210  continue
c     Christiansen et al.  (P_CRENBL_ECP)
      cc( 1, 1) =         -0.86027200d0
      cc( 2, 1) =         -8.20137298d0
      cc( 3, 1) =        -22.87043905d0
      cc( 4, 1) =         -6.70789498d0
      a ( 1, 1) =          2.10920000d0
      a ( 2, 1) =          5.96789998d0
      a ( 3, 1) =         17.52020001d0
      a ( 4, 1) =         53.42819977d0
      cc( 1, 2) =        -58.04999590d0
      cc( 2, 2) =        161.87587929d0
      cc( 3, 2) =        -98.19643497d0
      cc( 4, 2) =          6.34047300d0
      cc( 5, 2) =          3.02931401d0
      a ( 1, 2) =          2.28240001d0
      a ( 2, 2) =          2.62020001d0
      a ( 3, 2) =          3.10800001d0
      a ( 4, 2) =          5.00660002d0
      a ( 5, 2) =          3.27329999d0
      cc( 1, 3) =        127.97804737d0
      cc( 2, 3) =       -243.39212990d0
      cc( 3, 3) =        150.05020714d0
      cc( 4, 3) =          3.93873200d0
      cc( 5, 3) =          5.00513297d0
      a ( 1, 3) =          3.10389999d0
      a ( 2, 3) =          3.78320000d0
      a ( 3, 3) =          4.91720003d0
      a ( 4, 3) =          4.28490001d0
      a ( 5, 3) =         13.30770004d0
c
      go to 275
c
c     sulfur 
c
 230  continue
c     Christiansen et al. (S_CRENBL_ECP)
      cc( 1, 1) =         -1.09096500d0
      cc( 2, 1) =         -9.88182199d0
      cc( 3, 1) =        -26.41339898d0
      cc( 4, 1) =         -6.83946300d0
      a ( 1, 1) =          2.69159999d0
      a ( 2, 1) =          7.70099998d0
      a ( 3, 1) =         22.82830000d0
      a ( 4, 1) =         70.11800003d0
      cc( 1, 2) =        -59.91221809d0
      cc( 2, 2) =        182.12236595d0
      cc( 3, 2) =       -115.52312756d0
      cc( 4, 2) =          5.37397599d0
      cc( 5, 2) =          2.95619199d0
      a ( 1, 2) =          2.79040000d0
      a ( 2, 2) =          3.26820001d0
      a ( 3, 2) =          3.99610001d0
      a ( 4, 2) =          6.41450000d0
      a ( 5, 2) =          4.27050000d0
      cc( 1, 3) =        -48.71918821d0
      cc( 2, 3) =        140.84932518d0
      cc( 3, 3) =       -100.90417099d0
      cc( 4, 3) =          4.19983298d0
      cc( 5, 3) =          4.77143699d0
      a ( 1, 3) =          2.23629999d0
      a ( 2, 3) =          2.60190001d0
      a ( 3, 3) =          3.14840001d0
      a ( 4, 3) =          5.23589998d0
      a ( 5, 3) =          3.34340000d0
c
      go to 275
c
c     chlorine
c
 250  continue
c     Christiansen et al. (Cl_CRENBL_ECP)
      cc( 1, 1) =         -1.33688600d0
      cc( 2, 1) =        -11.74370098d0
      cc( 3, 1) =        -30.29555702d0
      cc( 4, 1) =         -6.98328000d0
      a ( 1, 1) =          3.35940000d0
      a ( 2, 1) =          9.73140001d0
      a ( 3, 1) =         29.25079989d0
      a ( 4, 1) =         90.54559994d0
      cc( 1, 2) =        -63.02188396d0
      cc( 2, 2) =        187.00407219d0
      cc( 3, 2) =       -114.41922283d0
      cc( 4, 2) =          7.64637101d0
      cc( 5, 2) =          2.99997699d0
      a ( 1, 2) =          3.14739999d0
      a ( 2, 2) =          3.68259999d0
      a ( 3, 2) =          4.48879999d0
      a ( 4, 2) =          7.18110001d0
      a ( 5, 2) =          4.77090001d0
      cc( 1, 3) =        -52.91793823d0
      cc( 2, 3) =        159.08480644d0
      cc( 3, 3) =       -118.43760872d0
      cc( 4, 3) =          5.19199800d0
      cc( 5, 3) =          4.67078602d0
      a ( 1, 3) =          2.70510000d0
      a ( 2, 3) =          3.20280001d0
      a ( 3, 3) =          3.98570001d0
      a ( 4, 3) =          6.57059997d0
      a ( 5, 3) =          4.28479999d0
c
      go to 275
c
c     argon
c
 270  continue
c     Christiansen et al. (Ar_CRENBL_ECP)
      cc( 1, 1) =         -1.63101099d0
      cc( 2, 1) =        -13.63825095d0
      cc( 3, 1) =        -34.38315582d0
      cc( 4, 1) =         -7.13807201d0
      a ( 1, 1) =          4.16430002d0
      a ( 2, 1) =         12.05250001d0
      a ( 3, 1) =         36.67759991d0
      a ( 4, 1) =        114.51200008d0
      cc( 1, 2) =         44.27679110d0
      cc( 2, 2) =       -119.46989155d0
      cc( 3, 2) =        107.43421745d0
      cc( 4, 2) =          1.31314400d0
      cc( 5, 2) =          3.35814700d0
      a ( 1, 2) =          2.13659999d0
      a ( 2, 2) =          2.42739999d0
      a ( 3, 2) =          2.96259999d0
      a ( 4, 2) =          2.01010001d0
      a ( 5, 2) =         11.53380001d0
      cc( 1, 3) =        -73.12434769d0
      cc( 2, 3) =        168.99321175d0
      cc( 3, 3) =       -117.48934937d0
      cc( 4, 3) =          3.61798900d0
      cc( 5, 3) =          4.84710997d0
      a ( 1, 3) =          2.68390000d0
      a ( 2, 3) =          3.28560001d0
      a ( 3, 3) =          4.10600001d0
      a ( 4, 3) =          7.93089998d0
      a ( 5, 3) =          1.48770000d0
275   n ( 1, 1) =   2
      n ( 2, 1) =   2
      n ( 3, 1) =   2
      n ( 4, 1) =   1
      nt( 1)    =   4
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      n ( 4, 2) =   1
      n ( 5, 2) =   0
      nt( 2)    =   5
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      n ( 3, 3) =   2
      n ( 4, 3) =   1
      n ( 5, 3) =   0
      nt( 3)    =   5
      nelcor    =  10
      lmx       =   2
c
 280  return
 6010 format (1x,'*** data not found for ecp type ',a7)
      end
**==crenbl1.f
      subroutine crenbl1(zinp,lmx,nelcor,nt,n,cc,a)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character * 8 (z)
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      dimension nt(6), n(10,6), cc(10,6), a(10,6)
c
      dimension ztit(18)
      data maxtyp/18/
      data ztit/'k','ca',
     *          'sc','ti','v ','cr','mn','fe',
     *          'co','ni','cu','zn','ga','ge','as','se','br','kr'/
c
c     this routine contains a library of ecp data taken from
c     various literature sources.
c     the CRENBL ecps are due to Christiansen etc. al
c     Large Orbital Basis, Small Core Potential
c     shape consistent ECPs
c     L.F. Pacios and P.A> Christiansen J.Chem. Phys. 82 (1985) 2664

      do 20 ityp = 1 , maxtyp
         if (zinp.eq.ztit(ityp)) go to 30
 20   continue
      write (iwr,6010) zinp
      call caserr2('no library available for requested ecp')
 30   go to (40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,
     +       200,210) , ityp
c     k  Christiansen et al. (K_CRENBL_ECP)
 40   continue
      cc( 1, 1) =         -1.96356797d0
      cc( 2, 1) =        -15.53031635d0
      cc( 3, 1) =        -38.57669830d0
      cc( 4, 1) =         -7.30101395d0
      a ( 1, 1) =          4.99889994d0
      a ( 2, 1) =         14.58220005d0
      a ( 3, 1) =         44.89559937d0
      a ( 4, 1) =        141.48350525d0
      cc( 1, 2) =         33.92499542d0
      cc( 2, 2) =       -117.52191162d0
      cc( 3, 2) =        108.25654602d0
      cc( 4, 2) =          5.80375385d0
      cc( 5, 2) =          3.21023202d0
      a ( 1, 2) =          2.41030002d0
      a ( 2, 2) =          2.77449989d0
      a ( 3, 2) =          3.33690000d0
      a ( 4, 2) =          2.22169995d0
      a ( 5, 2) =         12.19449997d0
      cc( 1, 3) =        -69.97570801d0
      cc( 2, 3) =        178.15260315d0
      cc( 3, 3) =       -136.11946106d0
      cc( 4, 3) =          4.43228197d0
      cc( 5, 3) =          4.74835920d0
      a ( 1, 3) =          3.15459990d0
      a ( 2, 3) =          4.03189993d0
      a ( 3, 3) =          5.16720009d0
      a ( 4, 3) =         11.33360004d0
      a ( 5, 3) =          1.74329996d0
c
      go to 55
c
c     ca Christiansen et al. (Ca_CRENBL_ECP)
c
 50   continue
      cc( 1, 1) =         -2.26729298d0
      cc( 2, 1) =        -17.17546654d0
      cc( 3, 1) =        -42.60916519d0
      cc( 4, 1) =         -7.46766806d0
      a ( 1, 1) =          5.92140007d0
      a ( 2, 1) =         17.20369911d0
      a ( 3, 1) =         53.37039948d0
      a ( 4, 1) =        169.98489380d0
      cc( 1, 2) =         29.85718918d0
      cc( 2, 2) =       -117.43736267d0
      cc( 3, 2) =        102.65641785d0
      cc( 4, 2) =         10.56573391d0
      cc( 5, 2) =          3.14651394d0
      a ( 1, 2) =          2.68330002d0
      a ( 2, 2) =          3.08299994d0
      a ( 3, 2) =          3.70659995d0
      a ( 4, 2) =          2.47650003d0
      a ( 5, 2) =         13.28890038d0
      cc( 1, 3) =        -77.17020416d0
      cc( 2, 3) =        180.64172363d0
      cc( 3, 3) =       -125.80883789d0
      cc( 4, 3) =          5.88094378d0
      cc( 5, 3) =          4.76948595d0
      a ( 1, 3) =          3.47429991d0
      a ( 2, 3) =          4.38199997d0
      a ( 3, 3) =          5.66209984d0
      a ( 4, 3) =          2.20009995d0
      a ( 5, 3) =          6.12169981d0
55    continue
      n ( 1, 1) =   2
      n ( 2, 1) =   2
      n ( 3, 1) =   2
      n ( 4, 1) =   1
      nt( 1)    =   4
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      n ( 4, 2) =   1
      n ( 5, 2) =   0
      nt( 2)    =   5
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      n ( 3, 3) =   2
      n ( 4, 3) =   1
      n ( 5, 3) =   0
      nt( 3)    =   5
      nelcor    =  10
      lmx       =   2
      go to 220
c
c     scandium  Christiansen et al.
c
 60   continue
      cc( 1, 1) =         -2.56915808d0
      cc( 2, 1) =        -19.51248169d0
      cc( 3, 1) =        -47.74710464d0
      cc( 4, 1) =         -7.68808413d0
      a ( 1, 1) =          6.93090010d0
      a ( 2, 1) =         20.47940063d0
      a ( 3, 1) =         64.68340302d0
      a ( 4, 1) =        208.68989563d0
      cc( 1, 2) =        -26.87937546d0
      cc( 2, 2) =         86.15200806d0
      cc( 3, 2) =       -158.35868835d0
      cc( 4, 2) =        153.65249634d0
      cc( 5, 2) =         -0.10035100d0
      cc( 6, 2) =          3.41862202d0
      a ( 1, 2) =          2.16470003d0
      a ( 2, 2) =          2.51990008d0
      a ( 3, 2) =          3.27999997d0
      a ( 4, 2) =          4.30459976d0
      a ( 5, 2) =          6.17159987d0
      a ( 6, 2) =         15.90830040d0
      cc( 1, 3) =         36.96915436d0
      cc( 2, 3) =       -117.41492462d0
      cc( 3, 3) =        213.54290771d0
      cc( 4, 3) =       -163.69456482d0
      cc( 5, 3) =          5.76216221d0
      cc( 6, 3) =          4.66992807d0
      a ( 1, 3) =          3.22720003d0
      a ( 2, 3) =          3.82780004d0
      a ( 3, 3) =          5.22049999d0
      a ( 4, 3) =          7.13719988d0
      a ( 5, 3) =         12.81589985d0
      a ( 6, 3) =          4.66620016d0
c
      go to 155
c
c     titanium  Christiansen et al.
c
 70   continue
      cc( 1, 1) =         -2.92140889d0
      cc( 2, 1) =        -21.51072121d0
      cc( 3, 1) =        -53.50305176d0
      cc( 4, 1) =         -8.13299370d0
      a ( 1, 1) =          8.02050018d0
      a ( 2, 1) =         23.73159981d0
      a ( 3, 1) =         76.42070007d0
      a ( 4, 1) =        255.15179443d0
      cc( 1, 2) =        -51.37191772d0
      cc( 2, 2) =        113.32749939d0
      cc( 3, 2) =       -173.03866577d0
      cc( 4, 2) =        146.56672668d0
      cc( 5, 2) =          6.99775124d0
      cc( 6, 2) =          3.13227606d0
      a ( 1, 2) =          2.09669995d0
      a ( 2, 2) =          2.54419994d0
      a ( 3, 2) =          3.42569995d0
      a ( 4, 2) =          4.62150002d0
      a ( 5, 2) =          1.50059998d0
      a ( 6, 2) =         18.84980011d0
      cc( 1, 3) =         32.56519318d0
      cc( 2, 3) =       -122.77304840d0
      cc( 3, 3) =        217.98690796d0
      cc( 4, 3) =       -168.72882080d0
      cc( 5, 3) =          5.54500103d0
      cc( 6, 3) =          4.74136782d0
      a ( 1, 3) =          3.53789997d0
      a ( 2, 3) =          4.22800016d0
      a ( 3, 3) =          5.74259996d0
      a ( 4, 3) =          7.77629995d0
      a ( 5, 3) =         11.05430031d0
      a ( 6, 3) =          3.26570010d0
c
      go to 155
c
c     vanadium  Christiansen et al.
c
 80   continue
      cc( 1, 1) =         -3.13644791d0
      cc( 2, 1) =        -23.11002541d0
      cc( 3, 1) =        -57.05287170d0
      cc( 4, 1) =         -8.13861370d0
      a ( 1, 1) =          8.96249962d0
      a ( 2, 1) =         26.80929947d0
      a ( 3, 1) =         86.54989624d0
      a ( 4, 1) =        287.72830200d0
      cc( 1, 2) =        -46.46123505d0
      cc( 2, 2) =         98.90882874d0
      cc( 3, 2) =       -175.34150696d0
      cc( 4, 2) =        152.27613831d0
      cc( 5, 2) =          9.72578144d0
      cc( 6, 2) =          3.10772490d0
      a ( 1, 2) =          2.22530007d0
      a ( 2, 2) =          2.78329992d0
      a ( 3, 2) =          3.80530000d0
      a ( 4, 2) =          5.07000017d0
      a ( 5, 2) =          1.67320001d0
      a ( 6, 2) =         18.82760048d0
      cc( 1, 3) =        -17.85890961d0
      cc( 2, 3) =         58.30934525d0
      cc( 3, 3) =       -111.78504944d0
      cc( 4, 3) =        115.28333282d0
      cc( 5, 3) =         -5.63639116d0
      cc( 6, 3) =          5.50943995d0
      a ( 1, 3) =          2.21219993d0
      a ( 2, 3) =          2.62529993d0
      a ( 3, 3) =          3.52760005d0
      a ( 4, 3) =          4.77939987d0
      a ( 5, 3) =          8.18389988d0
      a ( 6, 3) =         17.38699913d0
c
      go to 155
c
c     chromium  Christiansen et al.
c
 90   continue
      cc( 1, 1) =         -3.35951304d0
      cc( 2, 1) =        -24.58956909d0
      cc( 3, 1) =        -61.44187164d0
      cc( 4, 1) =         -8.37935352d0
      a ( 1, 1) =         10.01910019d0
      a ( 2, 1) =         29.92650032d0
      a ( 3, 1) =         97.39040375d0
      a ( 4, 1) =        329.05120850d0
      cc( 1, 2) =        -28.40607834d0
      cc( 2, 2) =         91.17616272d0
      cc( 3, 2) =       -171.88383484d0
      cc( 4, 2) =        181.42510986d0
      cc( 5, 2) =         10.19894505d0
      cc( 6, 2) =          3.05980206d0
      a ( 1, 2) =          3.08360004d0
      a ( 2, 2) =          3.57769990d0
      a ( 3, 2) =          4.65269995d0
      a ( 4, 2) =          6.15100002d0
      a ( 5, 2) =         23.50239944d0
      a ( 6, 2) =         20.16909981d0
      cc( 1, 3) =        -18.25805855d0
      cc( 2, 3) =         59.41814041d0
      cc( 3, 3) =       -113.84180450d0
      cc( 4, 3) =        117.49102783d0
      cc( 5, 3) =         -4.18536377d0
      cc( 6, 3) =          5.40804720d0
      a ( 1, 3) =          2.44000006d0
      a ( 2, 3) =          2.89219999d0
      a ( 3, 3) =          3.88299990d0
      a ( 4, 3) =          5.25540018d0
      a ( 5, 3) =          8.93529987d0
      a ( 6, 3) =         18.83979988d0
c
      go to 155
c
c     manganese  Christiansen et al.
c
 100  continue
      cc( 1, 1) =         -3.68953991d0
      cc( 2, 1) =        -26.47902679d0
      cc( 3, 1) =        -66.81269073d0
      cc( 4, 1) =         -8.68252087d0
      a ( 1, 1) =         11.27060032d0
      a ( 2, 1) =         33.69400024d0
      a ( 3, 1) =        110.90019989d0
      a ( 4, 1) =        381.80200195d0
      cc( 1, 2) =        -28.53082275d0
      cc( 2, 2) =         91.91989136d0
      cc( 3, 2) =       -174.68742371d0
      cc( 4, 2) =        189.56167603d0
      cc( 5, 2) =         12.27451992d0
      cc( 6, 2) =          3.10475302d0
      a ( 1, 2) =          3.39709997d0
      a ( 2, 2) =          3.94370008d0
      a ( 3, 2) =          5.13509989d0
      a ( 4, 2) =          6.81869984d0
      a ( 5, 2) =         24.90229988d0
      a ( 6, 2) =         21.42359924d0
      cc( 1, 3) =        -18.95299530d0
      cc( 2, 3) =         61.78389740d0
      cc( 3, 3) =       -117.40400696d0
      cc( 4, 3) =        120.17333984d0
      cc( 5, 3) =         -2.76534510d0
      cc( 6, 3) =          5.31945705d0
      a ( 1, 3) =          2.69440007d0
      a ( 2, 3) =          3.18889999d0
      a ( 3, 3) =          4.26399994d0
      a ( 4, 3) =          5.76149988d0
      a ( 5, 3) =          9.85439968d0
      a ( 6, 3) =         20.32609940d0
c
      go to 155
c
c     iron  Christiansen et al.
c
 110  continue
      cc( 1, 1) =         -3.92273593d0
      cc( 2, 1) =        -28.14956284d0
      cc( 3, 1) =        -71.93667603d0
      cc( 4, 1) =         -8.99319363d0
      a ( 1, 1) =         12.42240047d0
      a ( 2, 1) =         37.31629944d0
      a ( 3, 1) =        124.06909943d0
      a ( 4, 1) =        435.19369507d0
      cc( 1, 2) =        -29.52211952d0
      cc( 2, 2) =         94.51435852d0
      cc( 3, 2) =       -185.51301575d0
      cc( 4, 2) =        213.25717163d0
      cc( 5, 2) =          1.00094295d0
      cc( 6, 2) =          3.35134506d0
      a ( 1, 2) =          3.44029999d0
      a ( 2, 2) =          4.06059980d0
      a ( 3, 2) =          5.46460009d0
      a ( 4, 2) =          7.53329992d0
      a ( 5, 2) =          4.66279984d0
      a ( 6, 2) =         28.74029922d0
      cc( 1, 3) =        -19.13289261d0
      cc( 2, 3) =         62.18554688d0
      cc( 3, 3) =       -118.63009644d0
      cc( 4, 3) =        122.39747620d0
      cc( 5, 3) =         -1.59960997d0
      cc( 6, 3) =          5.25198078d0
      a ( 1, 3) =          2.94460011d0
      a ( 2, 3) =          3.48340011d0
      a ( 3, 3) =          4.66179991d0
      a ( 4, 3) =          6.30060005d0
      a ( 5, 3) =         10.56540012d0
      a ( 6, 3) =         21.84390068d0
c
      go to 155
c
c     cobalt  Christiansen et al.
c
 120  continue
      cc( 1, 1) =         -4.21805000d0
      cc( 2, 1) =        -30.00989914d0
      cc( 3, 1) =        -77.76090240d0
      cc( 4, 1) =         -9.37224770d0
      a ( 1, 1) =         13.73709965d0
      a ( 2, 1) =         41.34360123d0
      a ( 3, 1) =        139.10049438d0
      a ( 4, 1) =        499.13769531d0
      cc( 1, 2) =        -30.76871681d0
      cc( 2, 2) =         96.04718781d0
      cc( 3, 2) =       -186.83770752d0
      cc( 4, 2) =        214.13160706d0
      cc( 5, 2) =          3.17480111d0
      cc( 6, 2) =          3.30929589d0
      a ( 1, 2) =          3.74440002d0
      a ( 2, 2) =          4.40280008d0
      a ( 3, 2) =          5.94509983d0
      a ( 4, 2) =          8.18770027d0
      a ( 5, 2) =          7.63889980d0
      a ( 6, 2) =         29.83790016d0
      cc( 1, 3) =        -28.57496452d0
      cc( 2, 3) =         90.32809448d0
      cc( 3, 3) =       -157.74783325d0
      cc( 4, 3) =        137.59988403d0
      cc( 5, 3) =          1.21656799d0
      cc( 6, 3) =          5.09922504d0
      a ( 1, 3) =          2.91549993d0
      a ( 2, 3) =          3.64409995d0
      a ( 3, 3) =          4.85970020d0
      a ( 4, 3) =          6.52129984d0
      a ( 5, 3) =          1.70739996d0
      a ( 6, 3) =         24.11779976d0
c
      go to 155 
c
c     nickel (3d)9 (4s)1  Christiansen et al. ni01
c
 130  continue
      cc( 1, 1) =         -4.43806982d0
      cc( 2, 1) =        -32.09807205d0
      cc( 3, 1) =        -84.09606934d0
      cc( 4, 1) =         -9.79445744d0
      a ( 1, 1) =         14.88280010d0
      a ( 2, 1) =         45.52199936d0
      a ( 3, 1) =        155.59590149d0
      a ( 4, 1) =        571.42022705d0
      cc( 1, 2) =        -48.54310989d0
      cc( 2, 2) =        106.50110626d0
      cc( 3, 2) =       -202.67271423d0
      cc( 4, 2) =        218.64447021d0
      cc( 5, 2) =          7.88036299d0
      cc( 6, 2) =          3.17160606d0
      a ( 1, 2) =          3.39499998d0
      a ( 2, 2) =          4.34530020d0
      a ( 3, 2) =          6.15019989d0
      a ( 4, 2) =          8.67730045d0
      a ( 5, 2) =          2.48749995d0
      a ( 6, 2) =         31.87439919d0
      cc( 1, 3) =         54.26760101d0
      cc( 2, 3) =          1.34231699d0
      cc( 3, 3) =       -121.26470184d0
      cc( 4, 3) =        154.69761658d0
      cc( 5, 3) =         -6.76423597d0
      cc( 6, 3) =          5.44540691d0
      a ( 1, 3) =          4.04260015d0
      a ( 2, 3) =          1.78050005d0
      a ( 3, 3) =          5.53210020d0
      a ( 4, 3) =          7.91690016d0
      a ( 5, 3) =          2.22440004d0
      a ( 6, 3) =         30.32270050d0
      go to 155
c
c     copper  Christiansen et al.
c
 140  continue
      cc( 1, 1) =         -4.73227501d0
      cc( 2, 1) =        -34.06355667d0
      cc( 3, 1) =        -90.69224548d0
      cc( 4, 1) =        -10.26460838d0
      a ( 1, 1) =         16.30159950d0
      a ( 2, 1) =         49.98759842d0
      a ( 3, 1) =        173.02969360d0
      a ( 4, 1) =        651.10559082d0
      cc( 1, 2) =        -87.13957977d0
      cc( 2, 2) =        209.05120850d0
      cc( 3, 2) =       -202.30523682d0
      cc( 4, 2) =        154.84190369d0
      cc( 5, 2) =          9.21743488d0
      cc( 6, 2) =          3.18838096d0
      a ( 1, 2) =          3.70869994d0
      a ( 2, 2) =          4.51280022d0
      a ( 3, 2) =          5.53380013d0
      a ( 4, 2) =         10.20059967d0
      a ( 5, 2) =          2.66059995d0
      a ( 6, 2) =         32.17929840d0
      cc( 1, 3) =        -19.07518959d0
      cc( 2, 3) =         63.05695343d0
      cc( 3, 3) =       -127.18070221d0
      cc( 4, 3) =        158.41213989d0
      cc( 5, 3) =         -5.66128206d0
      cc( 6, 3) =          5.39882612d0
      a ( 1, 3) =          3.69499993d0
      a ( 2, 3) =          4.45380020d0
      a ( 3, 3) =          6.17630005d0
      a ( 4, 3) =          8.83930016d0
      a ( 5, 3) =         14.67029953d0
      a ( 6, 3) =         30.43350029d0
      go to 155
c
c     zinc  Christiansen et al. (Zn_CRENBL_ECP)
c
 150  continue
      cc( 1, 1) =         -4.95285416d0
      cc( 2, 1) =        -35.81423950d0
      cc( 3, 1) =        -97.18619537d0
      cc( 4, 1) =        -10.75999069d0
      a ( 1, 1) =         17.70450020d0
      a ( 2, 1) =         54.37659836d0
      a ( 3, 1) =        190.26429749d0
      a ( 4, 1) =        733.55322266d0
      cc( 1, 2) =        -25.24020195d0
      cc( 2, 2) =         76.24803925d0
      cc( 3, 2) =       -184.65945435d0
      cc( 4, 2) =        185.51107788d0
      cc( 5, 2) =         14.91602993d0
      cc( 6, 2) =          3.06773496d0
      a ( 1, 2) =          3.00970006d0
      a ( 2, 2) =          5.60309982d0
      a ( 3, 2) =          6.86539984d0
      a ( 4, 2) =         10.42599964d0
      a ( 5, 2) =          2.56110001d0
      a ( 6, 2) =         36.70050049d0
      cc( 1, 3) =        -28.08492470d0
      cc( 2, 3) =         99.61754608d0
      cc( 3, 3) =       -140.06907654d0
      cc( 4, 3) =        151.88270569d0
      cc( 5, 3) =         -3.59432101d0
      cc( 6, 3) =          5.29285192d0
      a ( 1, 3) =          4.11560011d0
      a ( 2, 3) =          4.90450001d0
      a ( 3, 3) =          6.25229979d0
      a ( 4, 3) =          9.47249985d0
      a ( 5, 3) =          5.86420012d0
      a ( 6, 3) =         34.27420044d0
 155  continue
      n ( 1, 1) =   2
      n ( 2, 1) =   2
      n ( 3, 1) =   2
      n ( 4, 1) =   1
      nt( 1)    =   4
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      n ( 4, 2) =   2
      n ( 5, 2) =   1
      n ( 6, 2) =   0
      nt( 2)    =   6
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      n ( 3, 3) =   2
      n ( 4, 3) =   2
      n ( 5, 3) =   1
      n ( 6, 3) =   0
      nt( 3)    =   6
      nelcor    =  10
      lmx       =   2
c
      go to 220
c
c     ga  Christiansen et al.  (Ga_CRENBL_ECP)
c
 160  continue
      cc( 1, 1) =         -0.77471000d0
      cc( 2, 1) =         -5.18385077d0
      cc( 3, 1) =        -16.36283493d0
      cc( 4, 1) =        -46.29306030d0
      cc( 5, 1) =        -12.36761570d0
      a ( 1, 1) =          0.97359997d0
      a ( 2, 1) =          3.25300002d0
      a ( 3, 1) =          9.12030029d0
      a ( 4, 1) =         30.20359993d0
      a ( 5, 1) =         87.78540039d0
      cc( 1, 2) =         47.45992279d0
      cc( 2, 2) =       -133.31851196d0
      cc( 3, 2) =        198.61817932d0
      cc( 4, 2) =       -111.65273285d0
      cc( 5, 2) =         26.52990341d0
      cc( 6, 2) =          3.20841599d0
      a ( 1, 2) =          1.29789996d0
      a ( 2, 2) =          1.46560001d0
      a ( 3, 2) =          1.80809999d0
      a ( 4, 2) =          2.20980000d0
      a ( 5, 2) =          3.41070008d0
      a ( 6, 2) =         13.75800037d0
      cc( 1, 3) =        -22.17953682d0
      cc( 2, 3) =        116.83462524d0
      cc( 3, 3) =       -173.27267456d0
      cc( 4, 3) =        116.79779053d0
      cc( 5, 3) =         13.53756237d0
      cc( 6, 3) =          5.37063789d0
      a ( 1, 3) =          1.60810006d0
      a ( 2, 3) =          1.83080006d0
      a ( 3, 3) =          2.31450009d0
      a ( 4, 3) =          2.96409988d0
      a ( 5, 3) =          6.94740009d0
      a ( 6, 3) =          9.34539986d0
      go to 215
c
c     ge Christiansen et al.  (Ge_CRENBL_ECP)
c
 170  continue
      cc( 1, 1) =         -1.40196705d0
      cc( 2, 1) =         -9.45403576d0
      cc( 3, 1) =        -23.62409592d0
      cc( 4, 1) =        -62.44671631d0
      cc( 5, 1) =        -12.96366024d0
      a ( 1, 1) =          1.35459995d0
      a ( 2, 1) =          5.02229977d0
      a ( 3, 1) =         15.32830048d0
      a ( 4, 1) =         47.47309875d0
      a ( 5, 1) =        136.24209595d0
      cc( 1, 2) =        -27.65140533d0
      cc( 2, 2) =        157.55535889d0
      cc( 3, 2) =       -208.90888977d0
      cc( 4, 2) =        146.75300598d0
      cc( 5, 2) =         26.79065895d0
      cc( 6, 2) =          2.94686103d0
      a ( 1, 2) =          2.30509996d0
      a ( 2, 2) =          2.67100000d0
      a ( 3, 2) =          3.49440002d0
      a ( 4, 2) =          4.66540003d0
      a ( 5, 2) =         14.92539978d0
      a ( 6, 2) =         44.04529953d0
      cc( 1, 3) =        -55.52495575d0
      cc( 2, 3) =        194.35308838d0
      cc( 3, 3) =       -249.20588684d0
      cc( 4, 3) =        154.47665405d0
      cc( 5, 3) =         13.02904606d0
      cc( 6, 3) =          5.37435579d0
      a ( 1, 3) =          1.72140002d0
      a ( 2, 3) =          1.98930001d0
      a ( 3, 3) =          2.55629992d0
      a ( 4, 3) =          3.29259992d0
      a ( 5, 3) =          8.26930046d0
      a ( 6, 3) =         10.14990044d0
      go to 215
c
c     as  Christiansen et al. (As_CRENBL_ECP)
c
 180  continue
      cc( 1, 1) =         -2.03465605d0
      cc( 2, 1) =        -13.08909798d0
      cc( 3, 1) =        -30.30925941d0
      cc( 4, 1) =        -73.83747101d0
      cc( 5, 1) =        -13.68201923d0
      a ( 1, 1) =          1.74109995d0
      a ( 2, 1) =          6.73740005d0
      a ( 3, 1) =         21.89489937d0
      a ( 4, 1) =         63.28400040d0
      a ( 5, 1) =        186.78129578d0
      cc( 1, 2) =         34.11351776d0
      cc( 2, 2) =       -109.88148499d0
      cc( 3, 2) =        198.61093140d0
      cc( 4, 2) =       -128.75582886d0
      cc( 5, 2) =         32.51686478d0
      cc( 6, 2) =          3.33781409d0
      a ( 1, 2) =          1.60979998d0
      a ( 2, 2) =          1.81799996d0
      a ( 3, 2) =          2.25640011d0
      a ( 4, 2) =          2.81340003d0
      a ( 5, 2) =          3.85339999d0
      a ( 6, 2) =         16.49069977d0
      cc( 1, 3) =        -68.50337982d0
      cc( 2, 3) =        231.13151550d0
      cc( 3, 3) =       -286.99252319d0
      cc( 4, 3) =        179.05894470d0
      cc( 5, 3) =         12.12996483d0
      cc( 6, 3) =          5.37297487d0
      a ( 1, 3) =          1.92700005d0
      a ( 2, 3) =          2.25600004d0
      a ( 3, 3) =          2.96729994d0
      a ( 4, 3) =          3.93980002d0
      a ( 5, 3) =         11.13199997d0
      a ( 6, 3) =         11.22399998d0
      go to 215
c
c     se Christiansen et al.  (Se_CRENBL_ECP)
c
 190  continue
      cc( 1, 1) =         -2.58638692d0
      cc( 2, 1) =        -15.59965801d0
      cc( 3, 1) =        -41.25695801d0
      cc( 4, 1) =        -98.63784790d0
      cc( 5, 1) =        -12.61432457d0
      a ( 1, 1) =          2.11820006d0
      a ( 2, 1) =          8.24320030d0
      a ( 3, 1) =         28.67779922d0
      a ( 4, 1) =         93.85810089d0
      a ( 5, 1) =        239.63980103d0
      cc( 1, 2) =        -76.46389771d0
      cc( 2, 2) =        269.16049194d0
      cc( 3, 2) =       -307.33093262d0
      cc( 4, 2) =        196.38562012d0
      cc( 5, 2) =         29.18940735d0
      cc( 6, 2) =          2.91098595d0
      a ( 1, 2) =          2.58559990d0
      a ( 2, 2) =          3.04340005d0
      a ( 3, 2) =          4.03289986d0
      a ( 4, 2) =          5.39890003d0
      a ( 5, 2) =         18.18700027d0
      a ( 6, 2) =         55.40280151d0
      cc( 1, 3) =        -72.78518677d0
      cc( 2, 3) =        252.06176758d0
      cc( 3, 3) =       -313.29913330d0
      cc( 4, 3) =        200.91072083d0
      cc( 5, 3) =         12.05568886d0
      cc( 6, 3) =          5.31686020d0
      a ( 1, 3) =          2.16829991d0
      a ( 2, 3) =          2.57850003d0
      a ( 3, 3) =          3.48589993d0
      a ( 4, 3) =          4.79790020d0
      a ( 5, 3) =         14.74899960d0
      a ( 6, 3) =         13.46450043d0
      go to 215
c
c     br  Christiansen et al. (Br_CRENBL_ECP)
c
 200  continue
      cc( 1, 1) =         -3.20746708d0
      cc( 2, 1) =        -17.81955910d0
      cc( 3, 1) =        -41.36774826d0
      cc( 4, 1) =       -102.00952911d0
      cc( 5, 1) =        -13.96452522d0
      a ( 1, 1) =          2.54780006d0
      a ( 2, 1) =          9.90429974d0
      a ( 3, 1) =         32.94390106d0
      a ( 4, 1) =         98.63279724d0
      a ( 5, 1) =        287.78079224d0
      cc( 1, 2) =         41.54266739d0
      cc( 2, 2) =       -130.86531067d0
      cc( 3, 2) =        229.11135864d0
      cc( 4, 2) =       -139.15881348d0
      cc( 5, 2) =         34.69166565d0
      cc( 6, 2) =          3.13891602d0
      a ( 1, 2) =          1.81379998d0
      a ( 2, 2) =          2.07699990d0
      a ( 3, 2) =          2.65140009d0
      a ( 4, 2) =          3.39140010d0
      a ( 5, 2) =          5.07940006d0
      a ( 6, 2) =         22.00860023d0
      cc( 1, 3) =        -68.31208801d0
      cc( 2, 3) =        236.41137695d0
      cc( 3, 3) =       -285.20968628d0
      cc( 4, 3) =        192.27485657d0
      cc( 5, 3) =         13.39648247d0
      cc( 6, 3) =          5.33779001d0
      a ( 1, 3) =          2.34380007d0
      a ( 2, 3) =          2.78999996d0
      a ( 3, 3) =          3.77500010d0
      a ( 4, 3) =          5.22399998d0
      a ( 5, 3) =         15.22420025d0
      a ( 6, 3) =         14.20930004d0
c
      go to 215
c
c     kr Christiansen et al.  (Kr_CRENBL_ECP)
c
 210  continue
      cc( 1, 1) =         -3.98946905d0
      cc( 2, 1) =        -21.25818825d0
      cc( 3, 1) =        -50.23062134d0
      cc( 4, 1) =       -113.97364807d0
      cc( 5, 1) =        -14.67047119d0
      a ( 1, 1) =          3.05509996d0
      a ( 2, 1) =         12.09599972d0
      a ( 3, 1) =         41.91159821d0
      a ( 4, 1) =        122.37729645d0
      a ( 5, 1) =        363.47369385d0
      cc( 1, 2) =        -78.36124420d0
      cc( 2, 2) =        278.31964111d0
      cc( 3, 2) =       -296.60931396d0
      cc( 4, 2) =        197.15113831d0
      cc( 5, 2) =         13.34094238d0
      cc( 6, 2) =          3.35239601d0
      a ( 1, 2) =          2.98720002d0
      a ( 2, 2) =          3.60640001d0
      a ( 3, 2) =          4.98120022d0
      a ( 4, 2) =          7.04050016d0
      a ( 5, 2) =         22.72100067d0
      a ( 6, 2) =         21.71269989d0
      cc( 1, 3) =        -73.70785522d0
      cc( 2, 3) =        262.28320312d0
      cc( 3, 3) =       -317.21502686d0
      cc( 4, 3) =        218.15017700d0
      cc( 5, 3) =         13.32087612d0
      cc( 6, 3) =          5.26835585d0
      a ( 1, 3) =          2.61249995d0
      a ( 2, 3) =          3.15750003d0
      a ( 3, 3) =          4.38999987d0
      a ( 4, 3) =          6.30530024d0
      a ( 5, 3) =         19.63229942d0
      a ( 6, 3) =         17.43339920d0
215   continue
      n ( 1, 1) =   2
      n ( 2, 1) =   2
      n ( 3, 1) =   2
      n ( 4, 1) =   2
      n ( 5, 1) =   1
      nt( 1)    =   5
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      n ( 4, 2) =   2
      n ( 5, 2) =   1
      n ( 6, 2) =   0
      nt( 2)    =   6
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      n ( 3, 3) =   2
      n ( 4, 3) =   2
      n ( 5, 3) =   1
      n ( 6, 3) =   0
      nt( 3)    =   6
      nelcor    =  18
      lmx       =   2
c
 220  return
 6010 format (1x,'*** data not found for ecp type ',a4)
      end
**==crenbl2.f
      subroutine crenbl2 (zinp,lmx,nelcor,nt,n,cc,a)
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
      implicit character * 8 (z)
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
      dimension nt(6), n(10,6), cc(10,6), a(10,6)
c
      dimension ztit(18)
      data maxtyp/18/
      data ztit/
     $ 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag','cd',
     $ 'in','sn','sb','te','i ','xe' /
c
c     this routine contains a library of ecp data taken from
c     various literature sources.
c     the CRENBL ecps are due to Christiansen etc. al
c     Large Orbital Basis, Small Core Potential
c     shape consistent ECPs
c     L.F. Pacios and P.A> Christiansen J.Chem. Phys. 82 (1985) 2664
c     
      do 10 ityp = 1 , maxtyp
         if (zinp.eq.ztit(ityp)) go to 20
 10   continue
      write (iwr,6010) zinp
      call caserr2('no library available for requested ecp')
 20   continue
c
      go to (37,38,39,40,41,42,43,44,45,46,47,48,49,
     +       50,51,52,53,54), ityp
c
c      rb (RB_CRENBL_ECP)
c
  37  continue
c
      cc( 1, 1) =         -1.04400003d0
      cc( 2, 1) =        -12.26854706d0
      cc( 3, 1) =        -40.49360657d0
      cc( 4, 1) =        -92.10794830d0
      cc( 5, 1) =        -20.25083160d0
      a ( 1, 1) =          1.96459997d0
      a ( 2, 1) =          5.02349997d0
      a ( 3, 1) =         12.31190014d0
      a ( 4, 1) =         39.43920136d0
      a ( 5, 1) =        116.43070221d0
      cc( 1, 2) =         50.81394196d0
      cc( 2, 2) =       -162.04731750d0
      cc( 3, 2) =        313.81082153d0
      cc( 4, 2) =       -309.75451660d0
      cc( 5, 2) =        216.07606506d0
      cc( 6, 2) =         20.86063194d0
      cc( 7, 2) =          3.36120105d0
      a ( 1, 2) =          2.29809999d0
      a ( 2, 2) =          2.66269994d0
      a ( 3, 2) =          3.50929999d0
      a ( 4, 2) =          4.96980000d0
      a ( 5, 2) =          6.94840002d0
      a ( 6, 2) =         17.70389938d0
      a ( 7, 2) =         25.66029930d0
      cc( 1, 3) =         45.41232300d0
      cc( 2, 3) =       -145.47238159d0
      cc( 3, 3) =        283.18420410d0
      cc( 4, 3) =       -305.10214233d0
      cc( 5, 3) =        207.65396118d0
      cc( 6, 3) =         12.15985012d0
      cc( 7, 3) =          5.39989424d0
      a ( 1, 3) =          2.02160001d0
      a ( 2, 3) =          2.33979988d0
      a ( 3, 3) =          3.07839990d0
      a ( 4, 3) =          4.37570000d0
      a ( 5, 3) =          6.15859985d0
      a ( 6, 3) =         16.77890015d0
      a ( 7, 3) =         16.61680031d0
      cc( 1, 4) =         31.68275070d0
      cc( 2, 4) =       -100.62529755d0
      cc( 3, 4) =        186.52160645d0
      cc( 4, 4) =       -239.76072693d0
      cc( 5, 4) =        170.19052124d0
      cc( 6, 4) =          9.91743755d0
      cc( 7, 4) =          7.41062880d0
      a ( 1, 4) =          1.23380005d0
      a ( 2, 4) =          1.41939998d0
      a ( 3, 4) =          1.83389997d0
      a ( 4, 4) =          2.54550004d0
      a ( 5, 4) =          3.47009993d0
      a ( 6, 4) =         10.62069988d0
      a ( 7, 4) =          9.28610039d0
c
      go to 100
c
c      sr (SR_CRENBL_ECP)
c
  38  continue
c
      cc( 1, 1) =         -1.16187501d0
      cc( 2, 1) =        -13.37399960d0
      cc( 3, 1) =        -43.23659134d0
      cc( 4, 1) =       -100.09903717d0
      cc( 5, 1) =        -20.51813126d0
      a ( 1, 1) =          2.23399997d0
      a ( 2, 1) =          5.67100000d0
      a ( 3, 1) =         13.88179970d0
      a ( 4, 1) =         44.81060028d0
      a ( 5, 1) =        132.90260315d0
      cc( 1, 2) =         53.58986664d0
      cc( 2, 2) =       -172.08218384d0
      cc( 3, 2) =        345.58593750d0
      cc( 4, 2) =       -351.22171021d0
      cc( 5, 2) =        257.34286499d0
      cc( 6, 2) =         12.91709232d0
      cc( 7, 2) =          6.33449411d0
      a ( 1, 2) =          2.44670010d0
      a ( 2, 2) =          2.86780000d0
      a ( 3, 2) =          3.86610007d0
      a ( 4, 2) =          5.66069984d0
      a ( 5, 2) =          8.30790043d0
      a ( 6, 2) =         23.49519920d0
      a ( 7, 2) =         21.03380013d0
      cc( 1, 3) =         49.15122604d0
      cc( 2, 3) =       -158.82582092d0
      cc( 3, 3) =        320.08462524d0
      cc( 4, 3) =       -349.10769653d0
      cc( 5, 3) =        239.32991028d0
      cc( 6, 3) =         11.86419868d0
      cc( 7, 3) =          5.33859777d0
      a ( 1, 3) =          2.20950007d0
      a ( 2, 3) =          2.58439994d0
      a ( 3, 3) =          3.46840000d0
      a ( 4, 3) =          5.06430006d0
      a ( 5, 3) =          7.37960005d0
      a ( 6, 3) =         22.22710037d0
      a ( 7, 3) =         19.91810036d0
      cc( 1, 4) =         32.67572403d0
      cc( 2, 4) =       -103.22133636d0
      cc( 3, 4) =        193.27650452d0
      cc( 4, 4) =       -248.63302612d0
      cc( 5, 4) =        183.39566040d0
      cc( 6, 4) =         10.14631939d0
      cc( 7, 4) =          7.38135815d0
      a ( 1, 4) =          1.40730000d0
      a ( 2, 4) =          1.62419999d0
      a ( 3, 4) =          2.11750007d0
      a ( 4, 4) =          2.97239995d0
      a ( 5, 4) =          4.11800003d0
      a ( 6, 4) =         12.52390003d0
      a ( 7, 4) =         10.88269997d0
c
      go to 100
c
c      yttrium (Y_CRENBL_ECP)
c
  39  continue
c
      cc( 1, 1) =         -2.22248793d0
      cc( 2, 1) =        -23.57674599d0
      cc( 3, 1) =        -59.28044128d0
      cc( 4, 1) =       -164.74046326d0
      cc( 5, 1) =        -22.87406540d0
      a ( 1, 1) =          2.90050006d0
      a ( 2, 1) =          7.94630003d0
      a ( 3, 1) =         21.59399986d0
      a ( 4, 1) =         78.34619904d0
      a ( 5, 1) =        265.06298828d0
      cc( 1, 2) =         55.64851379d0
      cc( 2, 2) =       -177.82318115d0
      cc( 3, 2) =        354.56311035d0
      cc( 4, 2) =       -340.00125122d0
      cc( 5, 2) =        239.03950500d0
      cc( 6, 2) =         14.65989113d0
      cc( 7, 2) =          3.39070010d0
      a ( 1, 2) =          2.68510008d0
      a ( 2, 2) =          3.13739991d0
      a ( 3, 2) =          4.20650005d0
      a ( 4, 2) =          6.10620022d0
      a ( 5, 2) =          8.79300022d0
      a ( 6, 2) =         26.35610008d0
      a ( 7, 2) =         26.61560059d0
      cc( 1, 3) =         50.28510666d0
      cc( 2, 3) =       -160.27282715d0
      cc( 3, 3) =        317.67678833d0
      cc( 4, 3) =       -336.20748901d0
      cc( 5, 3) =        237.31672668d0
      cc( 6, 3) =         12.75121307d0
      cc( 7, 3) =          5.37008190d0
      a ( 1, 3) =          2.37949991d0
      a ( 2, 3) =          2.77360010d0
      a ( 3, 3) =          3.70539999d0
      a ( 4, 3) =          5.37099981d0
      a ( 5, 3) =          7.78319979d0
      a ( 6, 3) =         22.16099930d0
      a ( 7, 3) =         20.47360039d0
      cc( 1, 4) =         33.29671097d0
      cc( 2, 4) =       -109.96274567d0
      cc( 3, 4) =        213.90904236d0
      cc( 4, 4) =       -273.61801147d0
      cc( 5, 4) =        207.30749512d0
      cc( 6, 4) =         10.78762054d0
      cc( 7, 4) =          7.32198191d0
      a ( 1, 4) =          1.58969998d0
      a ( 2, 4) =          1.85490000d0
      a ( 3, 4) =          2.44250011d0
      a ( 4, 4) =          3.47909999d0
      a ( 5, 4) =          4.95599985d0
      a ( 6, 4) =         15.23060036d0
      a ( 7, 4) =         13.37300015d0
c
      go to 100
c
c     zirconium (Zr_CRENBL_ECP)
c
  40  continue
c
      cc( 1, 1) =         -0.27778301d0
      cc( 2, 1) =         -6.66180706d0
      cc( 3, 1) =        -38.94683075d0
      cc( 4, 1) =        -86.41357422d0
      cc( 5, 1) =        -20.21328354d0
      a ( 1, 1) =          1.96560001d0
      a ( 2, 1) =          4.82639980d0
      a ( 3, 1) =         12.26459980d0
      a ( 4, 1) =         38.00030136d0
      a ( 5, 1) =        114.34359741d0
      cc( 1, 2) =         51.15351486d0
      cc( 2, 2) =       -166.17219543d0
      cc( 3, 2) =        354.62121582d0
      cc( 4, 2) =       -347.42883301d0
      cc( 5, 2) =        262.91113281d0
      cc( 6, 2) =         44.85594177d0
      cc( 7, 2) =          2.81595111d0
      a ( 1, 2) =          2.70860004d0
      a ( 2, 2) =          3.23639989d0
      a ( 3, 2) =          4.53420019d0
      a ( 4, 2) =          6.97310019d0
      a ( 5, 2) =         10.91759968d0
      a ( 6, 2) =         42.70869827d0
      a ( 7, 2) =        158.05920410d0
      cc( 1, 3) =         53.14321518d0
      cc( 2, 3) =       -172.25927734d0
      cc( 3, 3) =        357.74844360d0
      cc( 4, 3) =       -387.61001587d0
      cc( 5, 3) =        275.60794067d0
      cc( 6, 3) =         12.65426922d0
      cc( 7, 3) =          5.28926516d0
      a ( 1, 3) =          2.57399988d0
      a ( 2, 3) =          3.03749990d0
      a ( 3, 3) =          4.15030003d0
      a ( 4, 3) =          6.21470022d0
      a ( 5, 3) =          9.40369987d0
      a ( 6, 3) =         29.28240013d0
      a ( 7, 3) =         25.48209953d0
      cc( 1, 4) =         35.59209442d0
      cc( 2, 4) =       -115.19468689d0
      cc( 3, 4) =        225.52505493d0
      cc( 4, 4) =       -291.76715088d0
      cc( 5, 4) =        232.36845398d0
      cc( 6, 4) =         11.23965740d0
      cc( 7, 4) =          7.27583790d0
      a ( 1, 4) =          1.78410006d0
      a ( 2, 4) =          2.08969998d0
      a ( 3, 4) =          2.79329991d0
      a ( 4, 4) =          4.06930017d0
      a ( 5, 4) =          5.94059992d0
      a ( 6, 4) =         18.74440002d0
      a ( 7, 4) =         16.13680077d0
c
      go to 100
c
c     niobium (Nb_CRENBL_ECP)
c
  41  continue
c
      cc( 1, 1) =         -2.12992907d0
      cc( 2, 1) =        -22.97262573d0
      cc( 3, 1) =        -56.91117859d0
      cc( 4, 1) =       -145.99894714d0
      cc( 5, 1) =        -21.96123505d0
      a ( 1, 1) =          3.44460011d0
      a ( 2, 1) =          9.11340046d0
      a ( 3, 1) =         22.98679924d0
      a ( 4, 1) =         76.26629639d0
      a ( 5, 1) =        234.28889465d0
      cc( 1, 2) =         51.21849060d0
      cc( 2, 2) =       -165.87748718d0
      cc( 3, 2) =        353.43283081d0
      cc( 4, 2) =       -339.01223755d0
      cc( 5, 2) =        273.15057373d0
      cc( 6, 2) =         19.10134315d0
      cc( 7, 2) =          3.36632800d0
      a ( 1, 2) =          2.88310003d0
      a ( 2, 2) =          3.43950009d0
      a ( 3, 2) =          4.80849981d0
      a ( 4, 2) =          7.36780024d0
      a ( 5, 2) =         11.54769993d0
      a ( 6, 2) =         37.81430054d0
      a ( 7, 2) =         33.57780075d0
      cc( 1, 3) =         52.35134125d0
      cc( 2, 3) =       -168.72013855d0
      cc( 3, 3) =        349.30529785d0
      cc( 4, 3) =       -368.03564453d0
      cc( 5, 3) =        271.09564209d0
      cc( 6, 3) =         13.24002266d0
      cc( 7, 3) =          5.32180595d0
      a ( 1, 3) =          2.75320005d0
      a ( 2, 3) =          3.24040008d0
      a ( 3, 3) =          4.41370010d0
      a ( 4, 3) =          6.57270002d0
      a ( 5, 3) =          9.91119957d0
      a ( 6, 3) =         29.60670090d0
      a ( 7, 3) =         25.94729996d0
      cc( 1, 4) =         34.02252197d0
      cc( 2, 4) =       -109.34761810d0
      cc( 3, 4) =        211.43179321d0
      cc( 4, 4) =       -270.02197266d0
      cc( 5, 4) =        220.20515442d0
      cc( 6, 4) =         11.01616001d0
      cc( 7, 4) =          7.31999588d0
      a ( 1, 4) =          1.90470004d0
      a ( 2, 4) =          2.22099996d0
      a ( 3, 4) =          2.94840002d0
      a ( 4, 4) =          4.25290012d0
      a ( 5, 4) =          6.14580011d0
      a ( 6, 4) =         18.43989944d0
      a ( 7, 4) =         15.86740017d0
c
      go to 100
c
c     molybdenum (Mo_CRENBL_ECP)
c
  42  continue
c
      cc( 1, 1) =         -2.39912391d0
      cc( 2, 1) =        -26.34414482d0
      cc( 3, 1) =        -62.09140015d0
      cc( 4, 1) =       -163.40112305d0
      cc( 5, 1) =        -22.59148788d0
      a ( 1, 1) =          3.86059999d0
      a ( 2, 1) =         10.36299992d0
      a ( 3, 1) =         26.72909927d0
      a ( 4, 1) =         89.23100281d0
      a ( 5, 1) =        279.54220581d0
      cc( 1, 2) =         48.17292404d0
      cc( 2, 2) =       -156.71144104d0
      cc( 3, 2) =        332.59915161d0
      cc( 4, 2) =       -312.30667114d0
      cc( 5, 2) =        157.88722229d0
      cc( 6, 2) =         29.96883774d0
      cc( 7, 2) =          3.41751909d0
      a ( 1, 2) =          2.99819994d0
      a ( 2, 2) =          3.58340001d0
      a ( 3, 2) =          5.01289988d0
      a ( 4, 2) =          7.67030001d0
      a ( 5, 2) =         11.45119953d0
      a ( 6, 2) =          9.75979996d0
      a ( 7, 2) =         37.12519836d0
      cc( 1, 3) =         50.96643829d0
      cc( 2, 3) =       -164.71308899d0
      cc( 3, 3) =        342.27224731d0
      cc( 4, 3) =       -351.59188843d0
      cc( 5, 3) =        269.40219116d0
      cc( 6, 3) =         14.04963779d0
      cc( 7, 3) =          5.35105991d0
      a ( 1, 3) =          2.94350004d0
      a ( 2, 3) =          3.45880008d0
      a ( 3, 3) =          4.69640017d0
      a ( 4, 3) =          6.96299982d0
      a ( 5, 3) =         10.48050022d0
      a ( 6, 3) =         29.96489906d0
      a ( 7, 3) =         26.68639946d0
      cc( 1, 4) =         36.05135727d0
      cc( 2, 4) =       -117.10117340d0
      cc( 3, 4) =        232.92193604d0
      cc( 4, 4) =       -297.79098511d0
      cc( 5, 4) =        252.18852234d0
      cc( 6, 4) =         11.67608738d0
      cc( 7, 4) =          7.26251984d0
      a ( 1, 4) =          2.12299991d0
      a ( 2, 4) =          2.49699998d0
      a ( 3, 4) =          3.36770010d0
      a ( 4, 4) =          4.96120024d0
      a ( 5, 4) =          7.39949989d0
      a ( 6, 4) =         22.86809921d0
      a ( 7, 4) =         19.53870010d0
c
      go to 100
c
c     technetium (Tc_CRENBL_ECP)
c
  43  continue
c
      cc( 1, 1) =         -2.68787909d0
      cc( 2, 1) =        -28.78056145d0
      cc( 3, 1) =        -64.90345001d0
      cc( 4, 1) =       -172.06602478d0
      cc( 5, 1) =        -22.93936157d0
      a ( 1, 1) =          4.28879976d0
      a ( 2, 1) =         11.57050037d0
      a ( 3, 1) =         29.77589989d0
      a ( 4, 1) =         98.04190063d0
      a ( 5, 1) =        307.48928833d0
      cc( 1, 2) =        -27.43015671d0
      cc( 2, 2) =         89.00428009d0
      cc( 3, 2) =       -170.35900879d0
      cc( 4, 2) =        284.85977173d0
      cc( 5, 2) =       -137.86445618d0
      cc( 6, 2) =         20.88568306d0
      cc( 7, 2) =          3.27177596d0
      a ( 1, 2) =          2.15630007d0
      a ( 2, 2) =          2.51180005d0
      a ( 3, 2) =          3.30170012d0
      a ( 4, 2) =          4.74079990d0
      a ( 5, 2) =          6.52559996d0
      a ( 6, 2) =         10.52480030d0
      a ( 7, 2) =          7.10739994d0
      cc( 1, 3) =         52.82186508d0
      cc( 2, 3) =       -173.63539124d0
      cc( 3, 3) =        377.18582153d0
      cc( 4, 3) =       -400.31878662d0
      cc( 5, 3) =        238.46206665d0
      cc( 6, 3) =         17.92267990d0
      cc( 7, 3) =          5.26544809d0
      a ( 1, 3) =          3.08780003d0
      a ( 2, 3) =          3.68630004d0
      a ( 3, 3) =          5.14750004d0
      a ( 4, 3) =          7.94019985d0
      a ( 5, 3) =         12.23849964d0
      a ( 6, 3) =         10.65450001d0
      a ( 7, 3) =         35.27830124d0
      cc( 1, 4) =         38.29109192d0
      cc( 2, 4) =       -125.35187531d0
      cc( 3, 4) =        257.23062134d0
      cc( 4, 4) =       -330.65811157d0
      cc( 5, 4) =        291.91574097d0
      cc( 6, 4) =         12.69731045d0
      cc( 7, 4) =          7.19261312d0
      a ( 1, 4) =          2.34229994d0
      a ( 2, 4) =          2.78130007d0
      a ( 3, 4) =          3.82170010d0
      a ( 4, 4) =          5.76849985d0
      a ( 5, 4) =          8.93659973d0
      a ( 6, 4) =         28.75009918d0
      a ( 7, 4) =         24.47260094d0
c
      go to 100
c
c     ruthenium (Ru_CRENBL_ECP)
c
  44  continue
c
      cc( 1, 1) =         -3.01680589d0
      cc( 2, 1) =        -32.88051987d0
      cc( 3, 1) =        -73.65652466d0
      cc( 4, 1) =       -202.19105530d0
      cc( 5, 1) =        -24.11733246d0
      a ( 1, 1) =          4.78420019d0
      a ( 2, 1) =         13.04769993d0
      a ( 3, 1) =         35.18500137d0
      a ( 4, 1) =        119.48639679d0
      a ( 5, 1) =        392.59069824d0
      cc( 1, 2) =        -29.59748840d0
      cc( 2, 2) =         95.84513855d0
      cc( 3, 2) =       -184.26242065d0
      cc( 4, 2) =        326.30468750d0
      cc( 5, 2) =       -191.65345764d0
      cc( 6, 2) =         49.91513824d0
      cc( 7, 2) =          2.90960789d0
      a ( 1, 2) =          2.42680001d0
      a ( 2, 2) =          2.81320000d0
      a ( 3, 2) =          3.67639995d0
      a ( 4, 2) =          5.28630018d0
      a ( 5, 2) =          7.51270008d0
      a ( 6, 2) =         11.46560001d0
      a ( 7, 2) =         56.84280014d0
      cc( 1, 3) =        -24.95108986d0
      cc( 2, 3) =         86.32252502d0
      cc( 3, 3) =       -172.20999146d0
      cc( 4, 3) =        306.91882324d0
      cc( 5, 3) =       -195.84736633d0
      cc( 6, 3) =         16.67239571d0
      cc( 7, 3) =          5.11160898d0
      a ( 1, 3) =          2.04450011d0
      a ( 2, 3) =          2.42230010d0
      a ( 3, 3) =          3.22399998d0
      a ( 4, 3) =          4.73470020d0
      a ( 5, 3) =          6.77710009d0
      a ( 6, 3) =         11.11660004d0
      a ( 7, 3) =          7.63590002d0
      cc( 1, 4) =         37.18661118d0
      cc( 2, 4) =       -121.85943604d0
      cc( 3, 4) =        248.04385376d0
      cc( 4, 4) =       -313.83843994d0
      cc( 5, 4) =        282.90338135d0
      cc( 6, 4) =         12.31614971d0
      cc( 7, 4) =          7.23575020d0
      a ( 1, 4) =          2.49099994d0
      a ( 2, 4) =          2.94939995d0
      a ( 3, 4) =          4.02759981d0
      a ( 4, 4) =          6.02890015d0
      a ( 5, 4) =          9.24909973d0
      a ( 6, 4) =         28.46680069d0
      a ( 7, 4) =         24.18630028d0
c
      go to 100
c
c     rhodium (Rh_CRENBL_ECP)
c
  45  continue
c
      cc( 1, 1) =         -3.34650707d0
      cc( 2, 1) =        -36.27990723d0
      cc( 3, 1) =        -79.27825165d0
      cc( 4, 1) =       -219.38746643d0
      cc( 5, 1) =        -24.89464378d0
      a ( 1, 1) =          5.26849985d0
      a ( 2, 1) =         14.52709961d0
      a ( 3, 1) =         39.89099884d0
      a ( 4, 1) =        135.15080261d0
      a ( 5, 1) =        451.73019409d0
      cc( 1, 2) =        -26.99048614d0
      cc( 2, 2) =         88.46970367d0
      cc( 3, 2) =       -175.11592102d0
      cc( 4, 2) =        340.60363770d0
      cc( 5, 2) =       -223.21163940d0
      cc( 6, 2) =         47.58080292d0
      cc( 7, 2) =          2.79638195d0
      a ( 1, 2) =          2.33960009d0
      a ( 2, 2) =          2.77189994d0
      a ( 3, 2) =          3.75909996d0
      a ( 4, 2) =          5.70620012d0
      a ( 5, 2) =          8.72089958d0
      a ( 6, 2) =         13.09829998d0
      a ( 7, 2) =         91.69550323d0
      cc( 1, 3) =        -22.92908096d0
      cc( 2, 3) =         84.58609009d0
      cc( 3, 3) =       -171.19743347d0
      cc( 4, 3) =        316.75189209d0
      cc( 5, 3) =       -225.04205322d0
      cc( 6, 3) =         49.74481964d0
      cc( 7, 3) =          4.39694500d0
      a ( 1, 3) =          2.27749991d0
      a ( 2, 3) =          2.70799994d0
      a ( 3, 3) =          3.56419992d0
      a ( 4, 3) =          5.25099993d0
      a ( 5, 3) =          7.78830004d0
      a ( 6, 3) =         11.98419952d0
      a ( 7, 3) =         70.47000122d0
      cc( 1, 4) =         40.01486969d0
      cc( 2, 4) =       -131.79638672d0
      cc( 3, 4) =        277.17028809d0
      cc( 4, 4) =       -353.21633911d0
      cc( 5, 4) =        332.40838623d0
      cc( 6, 4) =         13.58090496d0
      cc( 7, 4) =          7.15987301d0
      a ( 1, 4) =          2.72009993d0
      a ( 2, 4) =          3.25489998d0
      a ( 3, 4) =          4.54110003d0
      a ( 4, 4) =          6.98999977d0
      a ( 5, 4) =         11.21280003d0
      a ( 6, 4) =         36.28559875d0
      a ( 7, 4) =         30.69669914d0
c
      go to 100
c
c     palladium (Pd_CRENBL_ECP)
c
  46  continue
c
      cc( 1, 1) =         -3.67986703d0
      cc( 2, 1) =        -39.08115387d0
      cc( 3, 1) =        -84.86928558d0
      cc( 4, 1) =       -237.43994141d0
      cc( 5, 1) =        -25.76420593d0
      a ( 1, 1) =          5.79349995d0
      a ( 2, 1) =         15.97259998d0
      a ( 3, 1) =         44.46939850d0
      a ( 4, 1) =        151.29159546d0
      a ( 5, 1) =        516.64819336d0
      cc( 1, 2) =        -27.05408478d0
      cc( 2, 2) =         88.65593719d0
      cc( 3, 2) =       -175.50939941d0
      cc( 4, 2) =        342.64071655d0
      cc( 5, 2) =       -222.09843445d0
      cc( 6, 2) =         51.81674576d0
      cc( 7, 2) =          2.82314610d0
      a ( 1, 2) =          2.48130012d0
      a ( 2, 2) =          2.93619990d0
      a ( 3, 2) =          3.97420001d0
      a ( 4, 2) =          6.02339983d0
      a ( 5, 2) =          9.19719982d0
      a ( 6, 2) =         13.76439953d0
      a ( 7, 2) =         89.98899841d0
      cc( 1, 3) =        -28.74325180d0
      cc( 2, 3) =         94.62859344d0
      cc( 3, 3) =       -189.17919922d0
      cc( 4, 3) =        380.36779785d0
      cc( 5, 3) =       -312.48147583d0
      cc( 6, 3) =         56.35126877d0
      cc( 7, 3) =          4.31076288d0
      a ( 1, 3) =          2.39150000d0
      a ( 2, 3) =          2.83590007d0
      a ( 3, 3) =          3.84990001d0
      a ( 4, 3) =          5.86049986d0
      a ( 5, 3) =          9.15100002d0
      a ( 6, 3) =         13.94559956d0
      a ( 7, 3) =        101.21410370d0
      cc( 1, 4) =         39.28127670d0
      cc( 2, 4) =       -129.23152161d0
      cc( 3, 4) =        269.75354004d0
      cc( 4, 4) =       -337.91909790d0
      cc( 5, 4) =        324.82962036d0
      cc( 6, 4) =         13.18474102d0
      cc( 7, 4) =          7.19891405d0
      a ( 1, 4) =          2.88260007d0
      a ( 2, 4) =          3.44020009d0
      a ( 3, 4) =          4.77540016d0
      a ( 4, 4) =          7.29589987d0
      a ( 5, 4) =         11.60599995d0
      a ( 6, 4) =         36.10689926d0
      a ( 7, 4) =         30.47970009d0
c
      go to 100
c
c     silver (Ag_CRENBL_ECP)
c
  47  continue
c
      cc( 1, 1) =         -3.58370399d0
      cc( 2, 1) =        -41.52972031d0
      cc( 3, 1) =        -90.43763733d0
      cc( 4, 1) =       -254.95243835d0
      cc( 5, 1) =        -26.65418816d0
      a ( 1, 1) =          6.26410007d0
      a ( 2, 1) =         17.34429932d0
      a ( 3, 1) =         49.05220032d0
      a ( 4, 1) =        167.69929504d0
      a ( 5, 1) =        583.98468018d0
      cc( 1, 2) =        -27.86531448d0
      cc( 2, 2) =         89.86640167d0
      cc( 3, 2) =       -176.23782349d0
      cc( 4, 2) =        323.62954712d0
      cc( 5, 2) =       -149.80087280d0
      cc( 6, 2) =         18.86888504d0
      cc( 7, 2) =          3.26756096d0
      a ( 1, 2) =          2.66310000d0
      a ( 2, 2) =          3.13039994d0
      a ( 3, 2) =          4.20359993d0
      a ( 4, 2) =          6.23330021d0
      a ( 5, 2) =          8.93859959d0
      a ( 6, 2) =         14.17240047d0
      a ( 7, 2) =          9.84739971d0
      cc( 1, 3) =        -27.75139236d0
      cc( 2, 3) =         94.34455872d0
      cc( 3, 3) =       -187.94503784d0
      cc( 4, 3) =        377.22006226d0
      cc( 5, 3) =       -298.96502686d0
      cc( 6, 3) =         56.06618118d0
      cc( 7, 3) =          4.33762121d0
      a ( 1, 3) =          2.52550006d0
      a ( 2, 3) =          3.00219989d0
      a ( 3, 3) =          4.04619980d0
      a ( 4, 3) =          6.16489983d0
      a ( 5, 3) =          9.54640007d0
      a ( 6, 3) =         14.71339989d0
      a ( 7, 3) =         99.15979767d0
      cc( 1, 4) =         38.58859253d0
      cc( 2, 4) =       -126.88704681d0
      cc( 3, 4) =        263.18853760d0
      cc( 4, 4) =       -325.30596924d0
      cc( 5, 4) =        319.68847656d0
      cc( 6, 4) =         12.95428467d0
      cc( 7, 4) =          7.22897387d0
      a ( 1, 4) =          3.03929996d0
      a ( 2, 4) =          3.62039995d0
      a ( 3, 4) =          5.00589991d0
      a ( 4, 4) =          7.60449982d0
      a ( 5, 4) =         12.02890015d0
      a ( 6, 4) =         36.31949997d0
      a ( 7, 4) =         30.56290054d0
c
      go to 100
c
c     cadmium (Cd_CRENBL_ECP)
c
  48  continue
c
      cc( 1, 1) =         -4.22736502d0
      cc( 2, 1) =        -45.16851807d0
      cc( 3, 1) =        -97.85390472d0
      cc( 4, 1) =       -277.36898804d0
      cc( 5, 1) =        -27.82247353d0
      a ( 1, 1) =          6.88749981d0
      a ( 2, 1) =         19.10289955d0
      a ( 3, 1) =         55.09080124d0
      a ( 4, 1) =        189.13699341d0
      a ( 5, 1) =        674.74859619d0
      cc( 1, 2) =        -29.60593224d0
      cc( 2, 2) =         96.45623016d0
      cc( 3, 2) =       -189.33224487d0
      cc( 4, 2) =        365.01492310d0
      cc( 5, 2) =       -208.19114685d0
      cc( 6, 2) =         50.42354965d0
      cc( 7, 2) =          2.75293803d0
      a ( 1, 2) =          2.95350003d0
      a ( 2, 2) =          3.45910001d0
      a ( 3, 2) =          4.60519981d0
      a ( 4, 2) =          6.83580017d0
      a ( 5, 2) =         10.11009979d0
      a ( 6, 2) =         15.29839993d0
      a ( 7, 2) =         97.10109711d0
      cc( 1, 3) =        -27.20086288d0
      cc( 2, 3) =         90.01723480d0
      cc( 3, 3) =       -179.77046204d0
      cc( 4, 3) =        347.07955933d0
      cc( 5, 3) =       -229.02343750d0
      cc( 6, 3) =         19.48693848d0
      cc( 7, 3) =          5.02399683d0
      a ( 1, 3) =          2.53099990d0
      a ( 2, 3) =          3.00780010d0
      a ( 3, 3) =          4.08690023d0
      a ( 4, 3) =          6.19070005d0
      a ( 5, 3) =          9.26280022d0
      a ( 6, 3) =         14.80700016d0
      a ( 7, 3) =         10.56299973d0
      cc( 1, 4) =         40.55107880d0
      cc( 2, 4) =       -133.15621948d0
      cc( 3, 4) =        284.73294067d0
      cc( 4, 4) =       -357.67636108d0
      cc( 5, 4) =        299.06478882d0
      cc( 6, 4) =         17.66106606d0
      cc( 7, 4) =          7.21160984d0
      a ( 1, 4) =          3.22659993d0
      a ( 2, 4) =          3.89599991d0
      a ( 3, 4) =          5.54120016d0
      a ( 4, 4) =          8.73200035d0
      a ( 5, 4) =         14.25979996d0
      a ( 6, 4) =         12.15629959d0
      a ( 7, 4) =         40.44020081d0
c
 100  continue
      n ( 1, 1) =   2
      n ( 2, 1) =   2
      n ( 3, 1) =   2
      n ( 4, 1) =   2
      n ( 5, 1) =   1
      nt( 1)    =   5
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      n ( 4, 2) =   2
      n ( 5, 2) =   2
      n ( 6, 2) =   1
      n ( 7, 2) =   0
      nt( 2)    =   7
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      n ( 3, 3) =   2
      n ( 4, 3) =   2
      n ( 5, 3) =   2
      n ( 6, 3) =   1
      n ( 7, 3) =   0
      nt( 3)    =   7
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      n ( 3, 4) =   2
      n ( 4, 4) =   2
      n ( 5, 4) =   2
      n ( 6, 4) =   1
      n ( 7, 4) =   0
      nt( 4)    =   7
      nelcor    =  28
      lmx       =   3
c
      return
c
c     indium (In_CRENBL_ECP)
c
  49  continue
c
      cc( 1, 1) =         -1.94056201d0
      cc( 2, 1) =        -11.47955036d0
      cc( 3, 1) =        -53.80504227d0
      cc( 4, 1) =        -22.34759712d0
      a ( 1, 1) =          1.34599996d0
      a ( 2, 1) =          3.56830001d0
      a ( 3, 1) =         14.11800003d0
      a ( 4, 1) =         41.56969833d0
      cc( 1, 2) =          0.07365500d0
      cc( 2, 2) =        -59.16038895d0
      cc( 3, 2) =        180.88560486d0
      cc( 4, 2) =       -114.79077148d0
      cc( 5, 2) =         15.04667282d0
      cc( 6, 2) =          6.86890602d0
      a ( 1, 2) =          0.19520000d0
      a ( 2, 2) =          1.17530000d0
      a ( 3, 2) =          1.40520000d0
      a ( 4, 2) =          1.70469999d0
      a ( 5, 2) =          2.88599992d0
      a ( 6, 2) =          4.41620016d0
      cc( 1, 3) =        -75.74547577d0
      cc( 2, 3) =        171.16397095d0
      cc( 3, 3) =       -196.06549072d0
      cc( 4, 3) =         87.76480865d0
      cc( 5, 3) =         24.37946129d0
      cc( 6, 3) =          5.66978121d0
      a ( 1, 3) =          0.82660002d0
      a ( 2, 3) =          1.04190004d0
      a ( 3, 3) =          1.35699999d0
      a ( 4, 3) =          1.70620000d0
      a ( 5, 3) =          0.60479999d0
      a ( 6, 3) =          7.74840021d0
      cc( 1, 4) =         -1.61824703d0
      cc( 2, 4) =          2.52116203d0
      cc( 3, 4) =        -42.73102570d0
      cc( 4, 4) =         41.48542786d0
      cc( 5, 4) =          8.59886932d0
      cc( 6, 4) =          7.80331993d0
      a ( 1, 4) =          0.82840002d0
      a ( 2, 4) =          1.47130001d0
      a ( 3, 4) =          2.41720009d0
      a ( 4, 4) =          3.16120005d0
      a ( 5, 4) =          1.45360005d0
      a ( 6, 4) =         13.23980045d0
c
      go to 200
c
c     tin (Sn_CRENBL_ECP)
c
  50  continue
c
      cc( 1, 1) =         -2.08150506d0
      cc( 2, 1) =        -12.32594395d0
      cc( 3, 1) =        -58.33305359d0
      cc( 4, 1) =        -22.59982109d0
      a ( 1, 1) =          1.49080002d0
      a ( 2, 1) =          3.96700001d0
      a ( 3, 1) =         15.61760044d0
      a ( 4, 1) =         46.37760162d0
      cc( 1, 2) =          4.67122316d0
      cc( 2, 2) =       -101.40917206d0
      cc( 3, 2) =        153.43458557d0
      cc( 4, 2) =        -99.63557434d0
      cc( 5, 2) =         42.69211578d0
      cc( 6, 2) =          3.70672011d0
      a ( 1, 2) =          0.87500000d0
      a ( 2, 2) =          1.07930005d0
      a ( 3, 2) =          1.31959999d0
      a ( 4, 2) =          1.61670005d0
      a ( 5, 2) =          0.89330000d0
      a ( 6, 2) =         13.34959984d0
      cc( 1, 3) =        -89.44316101d0
      cc( 2, 3) =        188.95716858d0
      cc( 3, 3) =       -210.16172791d0
      cc( 4, 3) =         98.57840729d0
      cc( 5, 3) =         24.45154953d0
      cc( 6, 3) =          5.68232679d0
      a ( 1, 3) =          0.93419999d0
      a ( 2, 3) =          1.16709995d0
      a ( 3, 3) =          1.57159996d0
      a ( 4, 3) =          1.99930000d0
      a ( 5, 3) =          0.67259997d0
      a ( 6, 3) =          8.51509953d0
      cc( 1, 4) =        -20.61536407d0
      cc( 2, 4) =         39.18358994d0
      cc( 3, 4) =        -76.06874084d0
      cc( 4, 4) =         59.03364944d0
      cc( 5, 4) =          8.33738613d0
      cc( 6, 4) =          7.81764507d0
      a ( 1, 4) =          1.23920000d0
      a ( 2, 4) =          1.67750001d0
      a ( 3, 4) =          2.41880012d0
      a ( 4, 4) =          3.30940008d0
      a ( 5, 4) =          1.00779998d0
      a ( 6, 4) =         14.44050026d0
c
      go to 200
c
c     antimony (Sb_CRENBL_ECP)
c
  51  continue
c
      cc( 1, 1) =         -2.23113608d0
      cc( 2, 1) =        -13.31894398d0
      cc( 3, 1) =        -63.36921310d0
      cc( 4, 1) =        -22.87770844d0
      a ( 1, 1) =          1.64499998d0
      a ( 2, 1) =          4.41020012d0
      a ( 3, 1) =         17.29479980d0
      a ( 4, 1) =         51.86349869d0
      cc( 1, 2) =        -91.23473358d0
      cc( 2, 2) =        234.76853943d0
      cc( 3, 2) =       -302.15496826d0
      cc( 4, 2) =        153.18350220d0
      cc( 5, 2) =         25.74815941d0
      cc( 6, 2) =          6.54234791d0
      a ( 1, 2) =          1.32270002d0
      a ( 2, 2) =          1.76600003d0
      a ( 3, 2) =          2.47189999d0
      a ( 4, 2) =          3.24390006d0
      a ( 5, 2) =          0.97100002d0
      a ( 6, 2) =         11.29850006d0
      cc( 1, 3) =         29.66554260d0
      cc( 2, 3) =        -90.75620270d0
      cc( 3, 3) =        147.41183472d0
      cc( 4, 3) =        -87.66103363d0
      cc( 5, 3) =         10.96317959d0
      cc( 6, 3) =          6.60084820d0
      a ( 1, 3) =          0.75700003d0
      a ( 2, 3) =          0.87760001d0
      a ( 3, 3) =          1.12709999d0
      a ( 4, 3) =          1.42659998d0
      a ( 5, 3) =          1.98370004d0
      a ( 6, 3) =          1.52310002d0
      cc( 1, 4) =        -27.32286644d0
      cc( 2, 4) =         70.00717926d0
      cc( 3, 4) =       -102.09190369d0
      cc( 4, 4) =         90.42972565d0
      cc( 5, 4) =          0.34869000d0
      cc( 6, 4) =          8.32599354d0
      a ( 1, 4) =          1.68490005d0
      a ( 2, 4) =          2.05110002d0
      a ( 3, 4) =          2.84559989d0
      a ( 4, 4) =          3.94829988d0
      a ( 5, 4) =          6.62400007d0
      a ( 6, 4) =         14.89630032d0
c
      go to 200
c
c     tellurium (Te_CRENBL_ECP)
c
  52  continue
c
      cc( 1, 1) =         -2.38358688d0
      cc( 2, 1) =        -14.21988106d0
      cc( 3, 1) =        -67.92436981d0
      cc( 4, 1) =        -23.13225365d0
      a ( 1, 1) =          1.81439996d0
      a ( 2, 1) =          4.86980009d0
      a ( 3, 1) =         18.94720078d0
      a ( 4, 1) =         57.30979919d0
      cc( 1, 2) =       -117.50439453d0
      cc( 2, 2) =        161.63246155d0
      cc( 3, 2) =       -120.28221893d0
      cc( 4, 2) =         25.30402184d0
      cc( 5, 2) =         49.01689529d0
      cc( 6, 2) =          5.21098709d0
      a ( 1, 2) =          1.33440006d0
      a ( 2, 2) =          1.62170005d0
      a ( 3, 2) =          2.21479988d0
      a ( 4, 2) =          3.17429996d0
      a ( 5, 2) =          1.07500005d0
      a ( 6, 2) =         16.05739975d0
      cc( 1, 3) =         33.98709869d0
      cc( 2, 3) =       -106.35176849d0
      cc( 3, 3) =        178.77839661d0
      cc( 4, 3) =       -110.42464447d0
      cc( 5, 3) =         29.88292694d0
      cc( 6, 3) =          5.44751978d0
      a ( 1, 3) =          0.86379999d0
      a ( 2, 3) =          1.00119996d0
      a ( 3, 3) =          1.28160000d0
      a ( 4, 3) =          1.62300003d0
      a ( 5, 3) =          2.40730000d0
      a ( 6, 3) =         12.56849957d0
      cc( 1, 4) =        -25.48574638d0
      cc( 2, 4) =         65.04901886d0
      cc( 3, 4) =        -95.47418976d0
      cc( 4, 4) =         84.64656830d0
      cc( 5, 4) =          2.66751909d0
      cc( 6, 4) =          8.15157413d0
      a ( 1, 4) =          1.86860001d0
      a ( 2, 4) =          2.25830007d0
      a ( 3, 4) =          3.07999992d0
      a ( 4, 4) =          4.19910002d0
      a ( 5, 4) =          7.17479992d0
      a ( 6, 4) =         15.52460003d0
c
      go to 200
c
c     iodine (I_CRENBL_ECP)
c
  53  continue
c
      cc( 1, 1) =         -2.60678697d0
      cc( 2, 1) =        -15.30154037d0
      cc( 3, 1) =        -73.32185364d0
      cc( 4, 1) =        -23.43061638d0
      a ( 1, 1) =          2.00659990d0
      a ( 2, 1) =          5.41200018d0
      a ( 3, 1) =         20.88789940d0
      a ( 4, 1) =         63.82109833d0
      cc( 1, 2) =       -144.29801941d0
      cc( 2, 2) =        252.49822998d0
      cc( 3, 2) =       -200.15640259d0
      cc( 4, 2) =         85.98131561d0
      cc( 5, 2) =         32.24239349d0
      cc( 6, 2) =          6.45168209d0
      a ( 1, 2) =          1.53009999d0
      a ( 2, 2) =          1.85660005d0
      a ( 3, 2) =          2.58430004d0
      a ( 4, 2) =          3.75419998d0
      a ( 5, 2) =          1.13870001d0
      a ( 6, 2) =         12.47010040d0
      cc( 1, 3) =         32.43734360d0
      cc( 2, 3) =       -101.57602692d0
      cc( 3, 3) =        169.25859070d0
      cc( 4, 3) =       -105.27653503d0
      cc( 5, 3) =         12.12060356d0
      cc( 6, 3) =          6.45733500d0
      a ( 1, 3) =          0.91109997d0
      a ( 2, 3) =          1.06579995d0
      a ( 3, 3) =          1.38989997d0
      a ( 4, 3) =          1.79579997d0
      a ( 5, 3) =          2.54789996d0
      a ( 6, 3) =          1.17250001d0
      cc( 1, 4) =          3.41730189d0
      cc( 2, 4) =         -3.97894001d0
      cc( 3, 4) =        -37.31042862d0
      cc( 4, 4) =         77.44996643d0
      cc( 5, 4) =         18.01509857d0
      cc( 6, 4) =          6.86909294d0
      a ( 1, 4) =          0.86530000d0
      a ( 2, 4) =          0.89960003d0
      a ( 3, 4) =          4.58720016d0
      a ( 4, 4) =          5.52820015d0
      a ( 5, 4) =         22.20800018d0
      a ( 6, 4) =         19.14010048d0
c
      go to 200
c
c     xeon (Xe_CRENBL_ECP)
c
  54  continue
c
      cc( 1, 1) =         -2.98907089d0
      cc( 2, 1) =        -16.86190414d0
      cc( 3, 1) =        -80.85955811d0
      cc( 4, 1) =        -23.83045006d0
      a ( 1, 1) =          2.24699998d0
      a ( 2, 1) =          6.12760019d0
      a ( 3, 1) =         23.48080063d0
      a ( 4, 1) =         72.73519897d0
      cc( 1, 2) =         36.93946457d0
      cc( 2, 2) =       -117.99051666d0
      cc( 3, 2) =        217.58146667d0
      cc( 4, 2) =       -132.93127441d0
      cc( 5, 2) =         13.86918831d0
      cc( 6, 2) =          7.20968485d0
      a ( 1, 2) =          1.24269998d0
      a ( 2, 2) =          1.45309997d0
      a ( 3, 2) =          1.92760003d0
      a ( 4, 2) =          2.56839991d0
      a ( 5, 2) =          3.60890007d0
      a ( 6, 2) =          2.81739998d0
      cc( 1, 3) =         38.26582718d0
      cc( 2, 3) =       -120.96846008d0
      cc( 3, 3) =        217.59863281d0
      cc( 4, 3) =       -144.18939209d0
      cc( 5, 3) =         33.80591965d0
      cc( 6, 3) =          5.23532295d0
      a ( 1, 3) =          1.02950001d0
      a ( 2, 3) =          1.20910001d0
      a ( 3, 3) =          1.59679997d0
      a ( 4, 3) =          2.10870004d0
      a ( 5, 3) =          3.15689993d0
      a ( 6, 3) =         17.37170029d0
      cc( 1, 4) =        -17.56745148d0
      cc( 2, 4) =         46.21330261d0
      cc( 3, 4) =        -80.85101318d0
      cc( 4, 4) =         97.02409363d0
      cc( 5, 4) =         -0.88147599d0
      cc( 6, 4) =          8.33844185d0
      a ( 1, 4) =          2.25279999d0
      a ( 2, 4) =          2.73140001d0
      a ( 3, 4) =          3.76239991d0
      a ( 4, 4) =          5.32089996d0
      a ( 5, 4) =          8.73659992d0
      a ( 6, 4) =         18.91500092d0
c
 200  continue
      n ( 1, 1) =   2
      n ( 2, 1) =   2
      n ( 3, 1) =   2
      n ( 4, 1) =   1
      nt( 1)    =   4
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      n ( 4, 2) =   2
      n ( 5, 2) =   1
      n ( 6, 2) =   0
      nt( 2)    =   6
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      n ( 3, 3) =   2
      n ( 4, 3) =   2
      n ( 5, 3) =   1
      n ( 6, 3) =   0
      nt( 3)    =   6
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      n ( 3, 4) =   2
      n ( 4, 4) =   2
      n ( 5, 4) =   1
      n ( 6, 4) =   0
      nt( 4)    =   6
      nelcor    =  36
      lmx       =   3
c
      return
 6010 format (1x,'*** data not found for CRENBL ecp type ',a4)
      end
**==crenbl3.f
      subroutine crenbl3 (zinp,lmx,nelcor,nt,n,cc,a)
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
      implicit character * 8 (z)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      dimension nt(6), n(10,6), cc(10,6), a(10,6)
c
c     this routine contains a library of ecp data taken from
c     various literature sources.
c     the CRENBL ecps are due to Christiansen etc. al
c     Large Orbital Basis, Small Core Potential
c     shape consistent ECPs
c     L.F. Pacios and P.A> Christiansen J.Chem. Phys. 82 (1985) 2664
c
      dimension ztit(17)
      data maxtyp/17/
      data ztit/'cs','ba','la',
     +          'ce','pr','nd','pm','sm','eu','gd','tb',
     +          'dy','ho','er','tm','yb','lu' /
c
c
      do 10 ityp = 1 , maxtyp
         if (zinp.eq.ztit(ityp)) go to 20
 10   continue
      write (iwr,6010) zinp
      call caserr2('no library available for requested ecp')
c
 20   continue
c
      go to (55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71), ityp
c
c   cesium (Cs_CRENBL_ECP)
c
 55   continue
      cc( 1, 1) =         -0.78916699d0
      cc( 2, 1) =         -8.42115784d0
      cc( 3, 1) =        -30.98544312d0
      cc( 4, 1) =        -95.03477478d0
      cc( 5, 1) =        -30.07960320d0
      a ( 1, 1) =          0.93849999d0
      a ( 2, 1) =          2.31629992d0
      a ( 3, 1) =          6.00729990d0
      a ( 4, 1) =         20.37969971d0
      a ( 5, 1) =         59.32889938d0
      cc( 1, 2) =         42.85466766d0
      cc( 2, 2) =       -138.00901794d0
      cc( 3, 2) =        275.99960327d0
      cc( 4, 2) =       -280.45663452d0
      cc( 5, 2) =        199.82038879d0
      cc( 6, 2) =         27.73096657d0
      cc( 7, 2) =          3.76870608d0
      a ( 1, 2) =          1.38530004d0
      a ( 2, 2) =          1.63240004d0
      a ( 3, 2) =          2.20580006d0
      a ( 4, 2) =          3.22149992d0
      a ( 5, 2) =          4.64960003d0
      a ( 6, 2) =         15.15250015d0
      a ( 7, 2) =         19.00049973d0
      cc( 1, 3) =         48.66250992d0
      cc( 2, 3) =       -145.70526123d0
      cc( 3, 3) =        264.46368408d0
      cc( 4, 3) =       -279.85159302d0
      cc( 5, 3) =        184.35585022d0
      cc( 6, 3) =         23.30001831d0
      cc( 7, 3) =          5.76792908d0
      a ( 1, 3) =          1.25950003d0
      a ( 2, 3) =          1.44169998d0
      a ( 3, 3) =          1.87639999d0
      a ( 4, 3) =          2.65750003d0
      a ( 5, 3) =          3.63870001d0
      a ( 6, 3) =         10.65320015d0
      a ( 7, 3) =         14.68060017d0
      cc( 1, 4) =         34.86072540d0
      cc( 2, 4) =       -106.79302979d0
      cc( 3, 4) =        188.23532104d0
      cc( 4, 4) =       -217.63992310d0
      cc( 5, 4) =        137.74559021d0
      cc( 6, 4) =         34.42418671d0
      cc( 7, 4) =          7.19875193d0
      a ( 1, 4) =          0.76810002d0
      a ( 2, 4) =          0.87459999d0
      a ( 3, 4) =          1.10769999d0
      a ( 4, 4) =          1.48230004d0
      a ( 5, 4) =          1.92019999d0
      a ( 6, 4) =          6.21829987d0
      a ( 7, 4) =         17.38809967d0
c
      go to 100
c
c   barium (Ba_CRENBL_ECP)
c
 56   continue
      cc( 1, 1) =         -0.88013703d0
      cc( 2, 1) =        -10.01861763d0
      cc( 3, 1) =        -35.70346451d0
      cc( 4, 1) =       -114.57715607d0
      cc( 5, 1) =        -30.99500656d0
      a ( 1, 1) =          0.97610003d0
      a ( 2, 1) =          2.66910005d0
      a ( 3, 1) =          7.10550022d0
      a ( 4, 1) =         24.84989929d0
      a ( 5, 1) =         75.09850311d0
      cc( 1, 2) =         52.18550110d0
      cc( 2, 2) =       -166.64633179d0
      cc( 3, 2) =        336.98910522d0
      cc( 4, 2) =       -346.60510254d0
      cc( 5, 2) =        229.66429138d0
      cc( 6, 2) =         20.49417496d0
      cc( 7, 2) =          6.64949989d0
      a ( 1, 2) =          1.51549995d0
      a ( 2, 2) =          1.80079997d0
      a ( 3, 2) =          2.46210003d0
      a ( 4, 2) =          3.64910007d0
      a ( 5, 2) =          5.34859991d0
      a ( 6, 2) =         17.15789986d0
      a ( 7, 2) =         16.02389908d0
      cc( 1, 3) =         61.51347351d0
      cc( 2, 3) =       -171.17402649d0
      cc( 3, 3) =        303.61636353d0
      cc( 4, 3) =       -324.12673950d0
      cc( 5, 3) =        210.71342468d0
      cc( 6, 3) =         19.11876488d0
      cc( 7, 3) =          5.84502220d0
      a ( 1, 3) =          1.35119998d0
      a ( 2, 3) =          1.55149996d0
      a ( 3, 3) =          2.06870008d0
      a ( 4, 3) =          3.00740004d0
      a ( 5, 3) =          4.23050022d0
      a ( 6, 3) =         14.21850014d0
      a ( 7, 3) =         13.10439968d0
      cc( 1, 4) =         34.92659378d0
      cc( 2, 4) =       -109.23178864d0
      cc( 3, 4) =        196.23254395d0
      cc( 4, 4) =       -224.63766479d0
      cc( 5, 4) =        146.97143555d0
      cc( 6, 4) =         37.01747894d0
      cc( 7, 4) =          7.09744883d0
      a ( 1, 4) =          0.89740002d0
      a ( 2, 4) =          1.02090001d0
      a ( 3, 4) =          1.28859997d0
      a ( 4, 4) =          1.73850000d0
      a ( 5, 4) =          2.28509998d0
      a ( 6, 4) =          7.33769989d0
      a ( 7, 4) =         20.39410019d0
c
      go to 100
c
c   lanthanum (La_CRENBL_ECP)
c
 57   continue
c
      cc( 1, 1) =         -1.78271198d0
      cc( 2, 1) =        -15.01564884d0
      cc( 3, 1) =        -42.60580063d0
      cc( 4, 1) =       -141.34611511d0
      cc( 5, 1) =        -32.09997559d0
      a ( 1, 1) =          1.30850005d0
      a ( 2, 1) =          3.59450006d0
      a ( 3, 1) =          9.54290009d0
      a ( 4, 1) =         32.37649918d0
      a ( 5, 1) =        100.62180328d0
      cc( 1, 2) =         57.45276260d0
      cc( 2, 2) =       -178.51110840d0
      cc( 3, 2) =        348.32040405d0
      cc( 4, 2) =       -348.10836792d0
      cc( 5, 2) =        238.82337952d0
      cc( 6, 2) =         21.93246651d0
      cc( 7, 2) =          6.68807220d0
      a ( 1, 2) =          1.63999999d0
      a ( 2, 2) =          1.92739999d0
      a ( 3, 2) =          2.60120010d0
      a ( 4, 2) =          3.80599999d0
      a ( 5, 2) =          5.52479982d0
      a ( 6, 2) =         17.20549965d0
      a ( 7, 2) =         16.66570091d0
      cc( 1, 3) =         48.54993820d0
      cc( 2, 3) =       -154.24778748d0
      cc( 3, 3) =        313.21530151d0
      cc( 4, 3) =       -352.77603149d0
      cc( 5, 3) =        255.42198181d0
      cc( 6, 3) =         20.53687859d0
      cc( 7, 3) =          8.26657486d0
      a ( 1, 3) =          1.37989998d0
      a ( 2, 3) =          1.62790000d0
      a ( 3, 3) =          2.21709990d0
      a ( 4, 3) =          3.29920006d0
      a ( 5, 3) =          4.93240023d0
      a ( 6, 3) =         16.82349968d0
      a ( 7, 3) =         14.97049999d0
      cc( 1, 4) =         37.32484817d0
      cc( 2, 4) =       -117.13645172d0
      cc( 3, 4) =        216.44523621d0
      cc( 4, 4) =       -252.23069763d0
      cc( 5, 4) =        168.61718750d0
      cc( 6, 4) =         41.21611023d0
      cc( 7, 4) =          6.93647718d0
      a ( 1, 4) =          0.94440001d0
      a ( 2, 4) =          1.08879995d0
      a ( 3, 4) =          1.40939999d0
      a ( 4, 4) =          1.94959998d0
      a ( 5, 4) =          2.63579988d0
      a ( 6, 4) =          9.00669956d0
      a ( 7, 4) =         25.81089973d0
C
 100  continue
      n ( 1, 1) =   2
      n ( 2, 1) =   2
      n ( 3, 1) =   2
      n ( 4, 1) =   2
      n ( 5, 1) =   1
      nt( 1)    =   5
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      n ( 4, 2) =   2
      n ( 5, 2) =   2
      n ( 6, 2) =   1
      n ( 7, 2) =   0
      nt( 2)    =   7
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      n ( 3, 3) =   2
      n ( 4, 3) =   2
      n ( 5, 3) =   2
      n ( 6, 3) =   1
      n ( 7, 3) =   0
      nt( 3)    =   7
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      n ( 3, 4) =   2
      n ( 4, 4) =   2
      n ( 5, 4) =   2
      n ( 6, 4) =   1
      n ( 7, 4) =   0
      nt( 4)    =   7
      nelcor    =  46
      lmx       =   3
c
      return
c
c   cerium (Ce_CRENBL_ECP)
c
 58   continue
      cc( 1, 1) =         -0.31545401d0
      cc( 2, 1) =         -4.03798485d0
      cc( 3, 1) =        -20.49917410d0
      cc( 4, 1) =        -81.90377040d0
      cc( 5, 1) =       -294.86666900d0
      cc( 6, 1) =        -30.30850030d0
      a ( 1, 1) =          0.18850000d0
      a ( 2, 1) =          0.69250000d0
      a ( 3, 1) =          2.45420003d0
      a ( 4, 1) =          9.75559997d0
      a ( 5, 1) =         50.36529920d0
      a ( 6, 1) =        141.28579700d0
      cc( 1, 2) =         87.62731930d0
      cc( 2, 2) =       -175.92492700d0
      cc( 3, 2) =        173.75613400d0
      cc( 4, 2) =       -118.99214900d0
      cc( 5, 2) =         80.98940270d0
      cc( 6, 2) =         -2.07819796d0
      cc( 7, 2) =         36.81498720d0
      cc( 8, 2) =          6.47752524d0
      a ( 1, 2) =          0.62040001d0
      a ( 2, 2) =          0.75950003d0
      a ( 3, 2) =          1.04569995d0
      a ( 4, 2) =          1.52219999d0
      a ( 5, 2) =          2.23819995d0
      a ( 6, 2) =          8.08870029d0
      a ( 7, 2) =          6.39750004d0
      a ( 8, 2) =         17.84670070d0
      cc( 1, 3) =         73.16195680d0
      cc( 2, 3) =        -98.24731440d0
      cc( 3, 3) =         47.34008030d0
      cc( 4, 3) =        -20.87126540d0
      cc( 5, 3) =         44.86721040d0
      cc( 6, 3) =         18.21543310d0
      cc( 7, 3) =         37.83178710d0
      cc( 8, 3) =          5.46868611d0
      a ( 1, 3) =          0.48100001d0
      a ( 2, 3) =          0.54939997d0
      a ( 3, 3) =          0.83929998d0
      a ( 4, 3) =          1.63699997d0
      a ( 5, 3) =          2.63610005d0
      a ( 6, 3) =          5.01919985d0
      a ( 7, 3) =          8.89319992d0
      a ( 8, 3) =         24.34210010d0
      cc( 1, 4) =         -0.01027300d0
      cc( 2, 4) =         -2.24907804d0
      cc( 3, 4) =         -0.02766900d0
      cc( 4, 4) =         98.89129640d0
      cc( 5, 4) =       -233.11683700d0
      cc( 6, 4) =        202.36198400d0
      cc( 7, 4) =         45.45883560d0
      cc( 8, 4) =          6.81029177d0
      a ( 1, 4) =          0.08600000d0
      a ( 2, 4) =          0.84719998d0
      a ( 3, 4) =          0.17870000d0
      a ( 4, 4) =          1.88349998d0
      a ( 5, 4) =          2.53220010d0
      a ( 6, 4) =          3.19779992d0
      a ( 7, 4) =         10.85770040d0
      a ( 8, 4) =         31.21439930d0
c
      go to 200
c
c  praseodymium
c
 59   continue
      cc( 1, 1) =         -0.28507000d0
      cc( 2, 1) =         -3.78517199d0
      cc( 3, 1) =        -19.60263440d0
      cc( 4, 1) =        -80.29591370d0
      cc( 5, 1) =       -304.50256400d0
      cc( 6, 1) =        -35.10055160d0
      a ( 1, 1) =          0.18340001d0
      a ( 2, 1) =          0.68919998d0
      a ( 3, 1) =          2.42919993d0
      a ( 4, 1) =          9.65989971d0
      a ( 5, 1) =         50.57759860d0
      a ( 6, 1) =        170.35189800d0
      cc( 1, 2) =         14.95254990d0
      cc( 2, 2) =        -43.42506030d0
      cc( 3, 2) =         73.12426760d0
      cc( 4, 2) =        -57.31110000d0
      cc( 5, 2) =         53.03722380d0
      cc( 6, 2) =         50.96369170d0
      cc( 7, 2) =         32.92785650d0
      cc( 8, 2) =          4.06864691d0
      a ( 1, 2) =          0.50809997d0
      a ( 2, 2) =          0.82110000d0
      a ( 3, 2) =          1.07540000d0
      a ( 4, 2) =          1.57959998d0
      a ( 5, 2) =          2.12700009d0
      a ( 6, 2) =          5.46960020d0
      a ( 7, 2) =         10.89529990d0
      a ( 8, 2) =         21.70459940d0
      cc( 1, 3) =         53.97642520d0
      cc( 2, 3) =       -116.57189900d0
      cc( 3, 3) =        126.61888900d0
      cc( 4, 3) =        -89.15218350d0
      cc( 5, 3) =         62.60961530d0
      cc( 6, 3) =        -36.40468220d0
      cc( 7, 3) =         46.55743030d0
      cc( 8, 3) =          5.11101008d0
      a ( 1, 3) =          0.50430000d0
      a ( 2, 3) =          0.64709997d0
      a ( 3, 3) =          0.90820003d0
      a ( 4, 3) =          1.32930005d0
      a ( 5, 3) =          1.97520006d0
      a ( 6, 3) =          6.69859982d0
      a ( 7, 3) =          5.61110020d0
      a ( 8, 3) =         24.73170090d0
      cc( 1, 4) =         -0.03509900d0
      cc( 2, 4) =          4.50882912d0
      cc( 3, 4) =        -52.60649110d0
      cc( 4, 4) =        153.83320600d0
      cc( 5, 4) =       -194.82063300d0
      cc( 6, 4) =        152.34256000d0
      cc( 7, 4) =         22.06917950d0
      cc( 8, 4) =          7.82782602d0
      a ( 1, 4) =          0.15080000d0
      a ( 2, 4) =          0.98780000d0
      a ( 3, 4) =          1.27530003d0
      a ( 4, 4) =          1.64900005d0
      a ( 5, 4) =          2.35910010d0
      a ( 6, 4) =          3.32590008d0
      a ( 7, 4) =          8.52200031d0
      a ( 8, 4) =         12.37780000d0
c
      go to 200
c
c   neodymium
c
 60   continue
      cc( 1, 1) =         -0.27911800d0
      cc( 2, 1) =         -3.64531303d0
      cc( 3, 1) =        -18.82904430d0
      cc( 4, 1) =        -78.47742460d0
      cc( 5, 1) =       -294.55139200d0
      cc( 6, 1) =        -36.00532530d0
      a ( 1, 1) =          0.18470000d0
      a ( 2, 1) =          0.69959998d0
      a ( 3, 1) =          2.42899990d0
      a ( 4, 1) =          9.59739971d0
      a ( 5, 1) =         49.50849920d0
      a ( 6, 1) =        169.03700300d0
      cc( 1, 2) =         10.78615860d0
      cc( 2, 2) =        -32.80133440d0
      cc( 3, 2) =         64.46156310d0
      cc( 4, 2) =        -37.78273010d0
      cc( 5, 2) =         34.32389070d0
      cc( 6, 2) =         40.26707840d0
      cc( 7, 2) =         36.91546250d0
      cc( 8, 2) =          4.29324913d0
      a ( 1, 2) =          0.49700001d0
      a ( 2, 2) =          0.91939998d0
      a ( 3, 2) =          1.18159997d0
      a ( 4, 2) =          1.66310000d0
      a ( 5, 2) =          2.40739989d0
      a ( 6, 2) =          4.56169987d0
      a ( 7, 2) =          9.25269985d0
      a ( 8, 2) =         21.05220030d0
      cc( 1, 3) =         42.02197270d0
      cc( 2, 3) =        -91.85214230d0
      cc( 3, 3) =         96.88398740d0
      cc( 4, 3) =        -43.81736380d0
      cc( 5, 3) =         43.62885670d0
      cc( 6, 3) =          6.82531881d0
      cc( 7, 3) =         37.60051350d0
      cc( 8, 3) =          5.58739185d0
      a ( 1, 3) =          0.50430000d0
      a ( 2, 3) =          0.66560000d0
      a ( 3, 3) =          0.92790002d0
      a ( 4, 3) =          1.31359994d0
      a ( 5, 3) =          2.66960001d0
      a ( 6, 3) =          4.85360003d0
      a ( 7, 3) =          7.76230002d0
      a ( 8, 3) =         23.10460090d0
      cc( 1, 4) =         -0.00342300d0
      cc( 2, 4) =         -0.05969700d0
      cc( 3, 4) =        -29.56184010d0
      cc( 4, 4) =        153.39161700d0
      cc( 5, 4) =       -240.44110100d0
      cc( 6, 4) =        188.30503800d0
      cc( 7, 4) =         47.79472730d0
      cc( 8, 4) =          6.79056502d0
      a ( 1, 4) =          0.07540000d0
      a ( 2, 4) =          0.27680001d0
      a ( 3, 4) =          1.41349995d0
      a ( 4, 4) =          1.88569999d0
      a ( 5, 4) =          2.57579994d0
      a ( 6, 4) =          3.49640012d0
      a ( 7, 4) =         11.78209970d0
      a ( 8, 4) =         34.12229920d0
c
      go to 200
c
c   promethium
c
 61   continue
      cc( 1, 1) =         -0.27581099d0
      cc( 2, 1) =         -3.67671704d0
      cc( 3, 1) =        -19.20204740d0
      cc( 4, 1) =        -80.51749420d0
      cc( 5, 1) =       -295.95886200d0
      cc( 6, 1) =        -32.22702030d0
      a ( 1, 1) =          0.18629999d0
      a ( 2, 1) =          0.72170001d0
      a ( 3, 1) =          2.52950001d0
      a ( 4, 1) =         10.04810050d0
      a ( 5, 1) =         51.84069830d0
      a ( 6, 1) =        154.48420700d0
      cc( 1, 2) =          0.02899800d0
      cc( 2, 2) =         59.00218200d0
      cc( 3, 2) =        -79.18837740d0
      cc( 4, 2) =         43.98767090d0
      cc( 5, 2) =         36.28269960d0
      cc( 6, 2) =         30.94334980d0
      cc( 7, 2) =         29.88560490d0
      cc( 8, 2) =          5.22673082d0
      a ( 1, 2) =          0.12690000d0
      a ( 2, 2) =          0.63349998d0
      a ( 3, 2) =          0.73229998d0
      a ( 4, 2) =          1.05410004d0
      a ( 5, 2) =          3.92600012d0
      a ( 6, 2) =          4.81040001d0
      a ( 7, 2) =         13.42710020d0
      a ( 8, 2) =          6.61800003d0
      cc( 1, 3) =         21.08691790d0
      cc( 2, 3) =        -63.05217740d0
      cc( 3, 3) =         87.14070130d0
      cc( 4, 3) =        -51.46301650d0
      cc( 5, 3) =         54.61180120d0
      cc( 6, 3) =          2.29255104d0
      cc( 7, 3) =         37.44199370d0
      cc( 8, 3) =          5.61524105d0
      a ( 1, 3) =          0.47749999d0
      a ( 2, 3) =          0.74510002d0
      a ( 3, 3) =          1.02049994d0
      a ( 4, 3) =          1.54830003d0
      a ( 5, 3) =          2.59549999d0
      a ( 6, 3) =          5.04430008d0
      a ( 7, 3) =          7.40789986d0
      a ( 8, 3) =         21.90730100d0
      cc( 1, 4) =         -0.40402901d0
      cc( 2, 4) =         30.99323080d0
      cc( 3, 4) =        -96.03415680d0
      cc( 4, 4) =        168.40460200d0
      cc( 5, 4) =       -190.01771500d0
      cc( 6, 4) =        150.01397700d0
      cc( 7, 4) =         36.84980780d0
      cc( 8, 4) =          7.18157816d0
      a ( 1, 4) =          0.43770000d0
      a ( 2, 4) =          1.01680005d0
      a ( 3, 4) =          1.22340000d0
      a ( 4, 4) =          1.65579999d0
      a ( 5, 4) =          2.43820000d0
      a ( 6, 4) =          3.39700008d0
      a ( 7, 4) =          9.34249973d0
      a ( 8, 4) =         21.38209920d0
c
      go to 200
c
c   samarium
c
 62   continue
      cc( 1, 1) =         -0.26969901d0
      cc( 2, 1) =         -3.51133609d0
      cc( 3, 1) =        -18.27518270d0
      cc( 4, 1) =        -77.45976250d0
      cc( 5, 1) =       -285.60421800d0
      cc( 6, 1) =        -36.69469070d0
      a ( 1, 1) =          0.18820000d0
      a ( 2, 1) =          0.72750002d0
      a ( 3, 1) =          2.50729990d0
      a ( 4, 1) =          9.85939979d0
      a ( 5, 1) =         49.69639970d0
      a ( 6, 1) =        169.41960100d0
      cc( 1, 2) =         25.61621090d0
      cc( 2, 2) =        -70.77632140d0
      cc( 3, 2) =         38.24119190d0
      cc( 4, 2) =         52.45853040d0
      cc( 5, 2) =        -49.38473510d0
      cc( 6, 2) =         63.82088850d0
      cc( 7, 2) =         38.89189910d0
      cc( 8, 2) =          6.54664612d0
      a ( 1, 2) =          0.61189997d0
      a ( 2, 2) =          0.91939998d0
      a ( 3, 2) =          1.21570003d0
      a ( 4, 2) =          1.22350001d0
      a ( 5, 2) =          1.94959998d0
      a ( 6, 2) =          3.04640007d0
      a ( 7, 2) =          7.50869990d0
      a ( 8, 2) =         19.97680090d0
      cc( 1, 3) =         29.42872240d0
      cc( 2, 3) =        -76.94017030d0
      cc( 3, 3) =        101.34484900d0
      cc( 4, 3) =        -58.27959440d0
      cc( 5, 3) =         50.80226900d0
      cc( 6, 3) =        -14.50483890d0
      cc( 7, 3) =         41.13980490d0
      cc( 8, 3) =          5.43980885d0
      a ( 1, 3) =          0.51340002d0
      a ( 2, 3) =          0.73739999d0
      a ( 3, 3) =          1.02600002d0
      a ( 4, 3) =          1.46430004d0
      a ( 5, 3) =          2.55809999d0
      a ( 6, 3) =          7.50089979d0
      a ( 7, 3) =          6.67859984d0
      a ( 8, 3) =         21.94910050d0
      cc( 1, 4) =        -25.71501730d0
      cc( 2, 4) =         71.07156370d0
      cc( 3, 4) =       -118.19715900d0
      cc( 4, 4) =        185.11428800d0
      cc( 5, 4) =       -194.87985200d0
      cc( 6, 4) =        145.04307600d0
      cc( 7, 4) =         28.46022030d0
      cc( 8, 4) =          7.62672281d0
      a ( 1, 4) =          0.79869998d0
      a ( 2, 4) =          0.92199999d0
      a ( 3, 4) =          1.20060003d0
      a ( 4, 4) =          1.68820000d0
      a ( 5, 4) =          2.45239997d0
      a ( 6, 4) =          3.46329999d0
      a ( 7, 4) =          8.39239979d0
      a ( 8, 4) =         15.12339970d0
c
      go to 200
c
c   europium
c
 63   continue
      cc( 1, 1) =         -0.26617599d0
      cc( 2, 1) =         -3.49868703d0
      cc( 3, 1) =        -18.43297770d0
      cc( 4, 1) =        -79.26859280d0
      cc( 5, 1) =       -299.62655600d0
      cc( 6, 1) =        -36.39295960d0
      a ( 1, 1) =          0.19030000d0
      a ( 2, 1) =          0.74500000d0
      a ( 3, 1) =          2.58220005d0
      a ( 4, 1) =         10.23359970d0
      a ( 5, 1) =         52.61579900d0
      a ( 6, 1) =        179.52040100d0
      cc( 1, 2) =         21.00108340d0
      cc( 2, 2) =        -57.80246350d0
      cc( 3, 2) =         43.95732500d0
      cc( 4, 2) =         46.48395540d0
      cc( 5, 2) =        -53.55662920d0
      cc( 6, 2) =         61.55118560d0
      cc( 7, 2) =         38.45267870d0
      cc( 8, 2) =          6.59646511d0
      a ( 1, 2) =          0.61290002d0
      a ( 2, 2) =          0.97149998d0
      a ( 3, 2) =          1.30719996d0
      a ( 4, 2) =          1.34809995d0
      a ( 5, 2) =          1.94980001d0
      a ( 6, 2) =          3.19490004d0
      a ( 7, 2) =          7.54370022d0
      a ( 8, 2) =         19.38649940d0
      cc( 1, 3) =         24.56742290d0
      cc( 2, 3) =        -74.12057490d0
      cc( 3, 3) =         97.57618710d0
      cc( 4, 3) =        -63.95512010d0
      cc( 5, 3) =         63.80558400d0
      cc( 6, 3) =          1.47463000d0
      cc( 7, 3) =         38.17859650d0
      cc( 8, 3) =          5.62949991d0
      a ( 1, 3) =          0.51599997d0
      a ( 2, 3) =          0.78479999d0
      a ( 3, 3) =          1.06200004d0
      a ( 4, 3) =          1.66129994d0
      a ( 5, 3) =          2.54099989d0
      a ( 6, 3) =          3.76889992d0
      a ( 7, 3) =          7.31090021d0
      a ( 8, 3) =         21.53230100d0
      cc( 1, 4) =        -33.31368260d0
      cc( 2, 4) =         90.53746030d0
      cc( 3, 4) =       -147.84429900d0
      cc( 4, 4) =        225.43699600d0
      cc( 5, 4) =       -234.19693000d0
      cc( 6, 4) =        171.55630500d0
      cc( 7, 4) =         26.45060160d0
      cc( 8, 4) =          7.73455191d0
      a ( 1, 4) =          0.89880002d0
      a ( 2, 4) =          1.02810001d0
      a ( 3, 4) =          1.32400000d0
      a ( 4, 4) =          1.83800006d0
      a ( 5, 4) =          2.67449999d0
      a ( 6, 4) =          3.81629992d0
      a ( 7, 4) =          9.48239994d0
      a ( 8, 4) =         14.85340020d0
c
      go to 200
c
c   gadolinium
c
 64   continue
      cc( 1, 1) =         -0.26610100d0
      cc( 2, 1) =         -3.53136396d0
      cc( 3, 1) =        -18.72265050d0
      cc( 4, 1) =        -81.53036500d0
      cc( 5, 1) =       -300.74014300d0
      cc( 6, 1) =        -33.50960920d0
      a ( 1, 1) =          0.19360000d0
      a ( 2, 1) =          0.76740003d0
      a ( 3, 1) =          2.67459989d0
      a ( 4, 1) =         10.69419960d0
      a ( 5, 1) =         54.95180130d0
      a ( 6, 1) =        168.31979400d0
      cc( 1, 2) =         19.09932330d0
      cc( 2, 2) =        -26.61323930d0
      cc( 3, 2) =        -26.27995870d0
      cc( 4, 2) =         82.21079250d0
      cc( 5, 2) =        -45.14988710d0
      cc( 6, 2) =         65.37348940d0
      cc( 7, 2) =         41.11628720d0
      cc( 8, 2) =          6.52553797d0
      a ( 1, 2) =          0.62239999d0
      a ( 2, 2) =          1.01269996d0
      a ( 3, 2) =          1.02670002d0
      a ( 4, 2) =          1.38409996d0
      a ( 5, 2) =          2.11680007d0
      a ( 6, 2) =          3.51329994d0
      a ( 7, 2) =          8.39589977d0
      a ( 8, 2) =         22.29739950d0
      cc( 1, 3) =         26.62279320d0
      cc( 2, 3) =        -72.51585390d0
      cc( 3, 3) =         99.13574980d0
      cc( 4, 3) =        -65.45091250d0
      cc( 5, 3) =         63.48794560d0
      cc( 6, 3) =        -16.52910610d0
      cc( 7, 3) =         37.74343490d0
      cc( 8, 3) =          5.63271189d0
      a ( 1, 3) =          0.53479999d0
      a ( 2, 3) =          0.79269999d0
      a ( 3, 3) =          1.10570002d0
      a ( 4, 3) =          1.64199996d0
      a ( 5, 3) =          2.63709998d0
      a ( 6, 3) =          4.77899981d0
      a ( 7, 3) =          6.12379980d0
      a ( 8, 3) =         18.10639950d0
      cc( 1, 4) =        -28.08826640d0
      cc( 2, 4) =         81.27649690d0
      cc( 3, 4) =       -135.07252500d0
      cc( 4, 4) =        202.14030500d0
      cc( 5, 4) =       -205.07235700d0
      cc( 6, 4) =        160.03721600d0
      cc( 7, 4) =         27.62781720d0
      cc( 8, 4) =          7.72167111d0
      a ( 1, 4) =          0.94540000d0
      a ( 2, 4) =          1.08700001d0
      a ( 3, 4) =          1.38240004d0
      a ( 4, 4) =          1.91030002d0
      a ( 5, 4) =          2.79839993d0
      a ( 6, 4) =          4.04769993d0
      a ( 7, 4) =          9.75979995d0
      a ( 8, 4) =         15.58259960d0
c
      go to 200
c
c   terbium
c
 65   continue
      cc( 1, 1) =         -0.26984900d0
      cc( 2, 1) =         -3.62558699d0
      cc( 3, 1) =        -19.32910350d0
      cc( 4, 1) =        -84.83631890d0
      cc( 5, 1) =       -308.04827900d0
      cc( 6, 1) =        -29.42491720d0
      a ( 1, 1) =          0.19800000d0
      a ( 2, 1) =          0.79670000d0
      a ( 3, 1) =          2.80559993d0
      a ( 4, 1) =         11.31449990d0
      a ( 5, 1) =         58.80989840d0
      a ( 6, 1) =        154.28930700d0
      cc( 1, 2) =          0.03446100d0
      cc( 2, 2) =         95.90180200d0
      cc( 3, 2) =       -172.62870800d0
      cc( 4, 2) =        111.80809800d0
      cc( 5, 2) =        -17.09731290d0
      cc( 6, 2) =         72.16942590d0
      cc( 7, 2) =         46.53331760d0
      cc( 8, 2) =          6.30115700d0
      a ( 1, 2) =          0.13300000d0
      a ( 2, 2) =          0.78710002d0
      a ( 3, 2) =          0.93940002d0
      a ( 4, 2) =          1.21809995d0
      a ( 5, 2) =          2.20289993d0
      a ( 6, 2) =          4.50530005d0
      a ( 7, 2) =         11.02239990d0
      a ( 8, 2) =         29.26510050d0
      cc( 1, 3) =         25.06032560d0
      cc( 2, 3) =        -70.25512690d0
      cc( 3, 3) =         96.37406920d0
      cc( 4, 3) =        -61.05334470d0
      cc( 5, 3) =         63.21340180d0
      cc( 6, 3) =         -3.11191893d0
      cc( 7, 3) =         38.19372940d0
      cc( 8, 3) =          5.65808010d0
      a ( 1, 3) =          0.54470003d0
      a ( 2, 3) =          0.82410002d0
      a ( 3, 3) =          1.14540005d0
      a ( 4, 3) =          1.72200000d0
      a ( 5, 3) =          2.80520010d0
      a ( 6, 3) =          4.27790022d0
      a ( 7, 3) =          7.36959982d0
      a ( 8, 3) =         20.52960010d0
      cc( 1, 4) =        -18.67657280d0
      cc( 2, 4) =         59.35923770d0
      cc( 3, 4) =       -102.18101500d0
      cc( 4, 4) =        161.22937000d0
      cc( 5, 4) =       -161.95156900d0
      cc( 6, 4) =        141.52011100d0
      cc( 7, 4) =         30.36456680d0
      cc( 8, 4) =          7.63114309d0
      a ( 1, 4) =          0.98070002d0
      a ( 2, 4) =          1.14110005d0
      a ( 3, 4) =          1.43460000d0
      a ( 4, 4) =          1.99930000d0
      a ( 5, 4) =          2.94009996d0
      a ( 6, 4) =          4.38520002d0
      a ( 7, 4) =         10.22750000d0
      a ( 8, 4) =         17.35890010d0
c
      go to 200
c
c   dysoprosium
c
 66   continue
      cc( 1, 1) =         -0.24610400d0
      cc( 2, 1) =         -3.12704992d0
      cc( 3, 1) =        -16.79536250d0
      cc( 4, 1) =        -73.63502500d0
      cc( 5, 1) =       -264.27121000d0
      cc( 6, 1) =        -37.48455430d0
      a ( 1, 1) =          0.19589999d0
      a ( 2, 1) =          0.76220000d0
      a ( 3, 1) =          2.58400011d0
      a ( 4, 1) =         10.09679990d0
      a ( 5, 1) =         48.73460010d0
      a ( 6, 1) =        164.56170700d0
      cc( 1, 2) =          0.03927000d0
      cc( 2, 2) =        119.06783300d0
      cc( 3, 2) =       -201.71113600d0
      cc( 4, 2) =        230.92334000d0
      cc( 5, 2) =       -175.07632400d0
      cc( 6, 2) =         95.23589320d0
      cc( 7, 2) =         45.16174700d0
      cc( 8, 2) =          6.41383982d0
      a ( 1, 2) =          0.13990000d0
      a ( 2, 2) =          0.83789998d0
      a ( 3, 2) =          0.99059999d0
      a ( 4, 2) =          1.47319996d0
      a ( 5, 2) =          1.91900003d0
      a ( 6, 2) =          3.12389994d0
      a ( 7, 2) =          8.81739998d0
      a ( 8, 2) =         25.33569910d0
      cc( 1, 3) =         27.72136120d0
      cc( 2, 3) =        -73.24303430d0
      cc( 3, 3) =         97.90841670d0
      cc( 4, 3) =        -65.15378570d0
      cc( 5, 3) =         62.79773330d0
      cc( 6, 3) =          0.42207801d0
      cc( 7, 3) =         36.59001160d0
      cc( 8, 3) =          5.73405981d0
      a ( 1, 3) =          0.56580001d0
      a ( 2, 3) =          0.83509999d0
      a ( 3, 3) =          1.17279995d0
      a ( 4, 3) =          1.76390004d0
      a ( 5, 3) =          2.75780010d0
      a ( 6, 3) =          3.79119992d0
      a ( 7, 3) =          7.19640017d0
      a ( 8, 3) =         18.36179920d0
      cc( 1, 4) =        -21.80728910d0
      cc( 2, 4) =         57.27083970d0
      cc( 3, 4) =       -100.43373100d0
      cc( 4, 4) =        175.60334800d0
      cc( 5, 4) =       -185.86096200d0
      cc( 6, 4) =        168.78459200d0
      cc( 7, 4) =         25.64606860d0
      cc( 8, 4) =          7.84217119d0
      a ( 1, 4) =          1.00370002d0
      a ( 2, 4) =          1.14460003d0
      a ( 3, 4) =          1.49460006d0
      a ( 4, 4) =          2.10730004d0
      a ( 5, 4) =          3.17759991d0
      a ( 6, 4) =          4.87099981d0
      a ( 7, 4) =         12.15419960d0
      a ( 8, 4) =         15.88690000d0
c
      go to 200
c
c   holmium
c
 67   continue
      cc( 1, 1) =         -0.24419101d0
      cc( 2, 1) =         -3.11587691d0
      cc( 3, 1) =        -16.65380860d0
      cc( 4, 1) =        -73.37957000d0
      cc( 5, 1) =       -258.84875500d0
      cc( 6, 1) =        -37.32218550d0
      a ( 1, 1) =          0.19880000d0
      a ( 2, 1) =          0.77969998d0
      a ( 3, 1) =          2.63420010d0
      a ( 4, 1) =         10.26840020d0
      a ( 5, 1) =         48.84769820d0
      a ( 6, 1) =        162.01890600d0
      cc( 1, 2) =          0.04005300d0
      cc( 2, 2) =        123.58869200d0
      cc( 3, 2) =       -236.74517800d0
      cc( 4, 2) =        171.58247400d0
      cc( 5, 2) =        -54.66671370d0
      cc( 6, 2) =         80.29673000d0
      cc( 7, 2) =         47.89864350d0
      cc( 8, 2) =          6.27561522d0
      a ( 1, 2) =          0.14200001d0
      a ( 2, 2) =          0.86420000d0
      a ( 3, 2) =          1.02839994d0
      a ( 4, 2) =          1.33249998d0
      a ( 5, 2) =          2.17129994d0
      a ( 6, 2) =          4.00290012d0
      a ( 7, 2) =         10.56019970d0
      a ( 8, 2) =         28.76339910d0
      cc( 1, 3) =         10.00089930d0
      cc( 2, 3) =        -18.15449910d0
      cc( 3, 3) =          4.20426417d0
      cc( 4, 3) =         42.85946280d0
      cc( 5, 3) =        -37.44625860d0
      cc( 6, 3) =         62.10699080d0
      cc( 7, 3) =         42.81314470d0
      cc( 8, 3) =          5.61891890d0
      a ( 1, 3) =          0.48590002d0
      a ( 2, 3) =          0.89020002d0
      a ( 3, 3) =          1.39730001d0
      a ( 4, 3) =          1.40050006d0
      a ( 5, 3) =          2.07459998d0
      a ( 6, 3) =          3.35479999d0
      a ( 7, 3) =          8.65760040d0
      a ( 8, 3) =         25.96680070d0
      cc( 1, 4) =        -17.04135320d0
      cc( 2, 4) =         47.99711230d0
      cc( 3, 4) =        -84.43540190d0
      cc( 4, 4) =        149.89463800d0
      cc( 5, 4) =       -158.44630400d0
      cc( 6, 4) =        160.95439200d0
      cc( 7, 4) =         26.41877370d0
      cc( 8, 4) =          7.84279013d0
      a ( 1, 4) =          1.04910004d0
      a ( 2, 4) =          1.20220006d0
      a ( 3, 4) =          1.54100001d0
      a ( 4, 4) =          2.18689990d0
      a ( 5, 4) =          3.33369994d0
      a ( 6, 4) =          5.20219994d0
      a ( 7, 4) =         12.70359990d0
      a ( 8, 4) =         16.47439960d0
c
      go to 200
c
c   erbium
c
 68   continue
c
      cc( 1, 1) =         -0.24749100d0
      cc( 2, 1) =         -3.12518692d0
      cc( 3, 1) =        -16.71219250d0
      cc( 4, 1) =        -73.93254850d0
      cc( 5, 1) =       -266.99044800d0
      cc( 6, 1) =        -39.54242710d0
      a ( 1, 1) =          0.20420000d0
      a ( 2, 1) =          0.80059999d0
      a ( 3, 1) =          2.70350003d0
      a ( 4, 1) =         10.54500010d0
      a ( 5, 1) =         50.36610030d0
      a ( 6, 1) =        179.59350600d0
      cc( 1, 2) =          0.04168000d0
      cc( 2, 2) =        118.45232400d0
      cc( 3, 2) =       -215.00274700d0
      cc( 4, 2) =        177.66186500d0
      cc( 5, 2) =       -124.16603900d0
      cc( 6, 2) =        118.82530200d0
      cc( 7, 2) =         50.10784910d0
      cc( 8, 2) =          6.37769508d0
      a ( 1, 2) =          0.14760000d0
      a ( 2, 2) =          0.88929999d0
      a ( 3, 2) =          1.06110001d0
      a ( 4, 2) =          1.45819998d0
      a ( 5, 2) =          2.33829999d0
      a ( 6, 2) =          3.34649992d0
      a ( 7, 2) =          9.93130016d0
      a ( 8, 2) =         31.89080050d0
      cc( 1, 3) =         10.00590610d0
      cc( 2, 3) =        -32.43202210d0
      cc( 3, 3) =        282.26181000d0
      cc( 4, 3) =       -224.79307600d0
      cc( 5, 3) =        -39.97660450d0
      cc( 6, 3) =         66.70326230d0
      cc( 7, 3) =         41.70607380d0
      cc( 8, 3) =          5.66391897d0
      a ( 1, 3) =          0.49800000d0
      a ( 2, 3) =          0.99419999d0
      a ( 3, 3) =          1.31599999d0
      a ( 4, 3) =          1.31579995d0
      a ( 5, 3) =          2.20199990d0
      a ( 6, 3) =          3.26659989d0
      a ( 7, 3) =          8.36079979d0
      a ( 8, 3) =         23.72389980d0
      cc( 1, 4) =        -10.52753640d0
      cc( 2, 4) =        -12.90994930d0
      cc( 3, 4) =        102.43170900d0
      cc( 4, 4) =       -124.21231100d0
      cc( 5, 4) =        165.57734700d0
      cc( 6, 4) =        -24.26627920d0
      cc( 7, 4) =         60.33625410d0
      cc( 8, 4) =          6.67504692d0
      a ( 1, 4) =          1.72560000d0
      a ( 2, 4) =          1.74810004d0
      a ( 3, 4) =          2.39429998d0
      a ( 4, 4) =          3.59179998d0
      a ( 5, 4) =          5.86060000d0
      a ( 6, 4) =         10.48060040d0
      a ( 7, 4) =         16.07609940d0
      a ( 8, 4) =         49.35649870d0
c
      go to 200
c
c   thulium
c
 69   continue
c
      cc( 1, 1) =         -0.24459700d0
      cc( 2, 1) =         -3.04991007d0
      cc( 3, 1) =        -16.31619640d0
      cc( 4, 1) =        -71.80751800d0
      cc( 5, 1) =       -219.44210800d0
      cc( 6, 1) =        -32.14739230d0
      a ( 1, 1) =          0.20750000d0
      a ( 2, 1) =          0.81089997d0
      a ( 3, 1) =          2.71910000d0
      a ( 4, 1) =         10.51509950d0
      a ( 5, 1) =         46.68840030d0
      a ( 6, 1) =        119.95919800d0
      cc( 1, 2) =          0.04336000d0
      cc( 2, 2) =        116.32991800d0
      cc( 3, 2) =       -230.19224600d0
      cc( 4, 2) =        228.38797000d0
      cc( 5, 2) =       -148.97734100d0
      cc( 6, 2) =        106.57288400d0
      cc( 7, 2) =         46.43137740d0
      cc( 8, 2) =          6.47452116d0
      a ( 1, 2) =          0.15210000d0
      a ( 2, 2) =          0.91829997d0
      a ( 3, 2) =          1.11269999d0
      a ( 4, 2) =          1.52339995d0
      a ( 5, 2) =          2.18759990d0
      a ( 6, 2) =          3.32069993d0
      a ( 7, 2) =          9.13409996d0
      a ( 8, 2) =         25.98240090d0
      cc( 1, 3) =          0.03231000d0
      cc( 2, 3) =        112.35777300d0
      cc( 3, 3) =       -112.68846100d0
      cc( 4, 3) =        190.72164900d0
      cc( 5, 3) =       -206.54730200d0
      cc( 6, 3) =         99.58874510d0
      cc( 7, 3) =         40.48825460d0
      cc( 8, 3) =          5.48724222d0
      a ( 1, 3) =          0.12780000d0
      a ( 2, 3) =          0.63319999d0
      a ( 3, 3) =          0.66270000d0
      a ( 4, 3) =          2.03090000d0
      a ( 5, 3) =          2.27940011d0
      a ( 6, 3) =          3.92199993d0
      a ( 7, 3) =         10.70800020d0
      a ( 8, 3) =         21.23900030d0
      cc( 1, 4) =        -10.80157280d0
      cc( 2, 4) =         29.62875560d0
      cc( 3, 4) =        -56.37639620d0
      cc( 4, 4) =        117.18654600d0
      cc( 5, 4) =       -119.45823700d0
      cc( 6, 4) =        146.09451300d0
      cc( 7, 4) =         29.72262000d0
      cc( 8, 4) =          7.77486515d0
      a ( 1, 4) =          1.12740004d0
      a ( 2, 4) =          1.28230000d0
      a ( 3, 4) =          1.63870001d0
      a ( 4, 4) =          2.34240007d0
      a ( 5, 4) =          3.57110000d0
      a ( 6, 4) =          5.80569983d0
      a ( 7, 4) =         13.21090030d0
      a ( 8, 4) =         18.70999910d0
c
      go to 200
c
c   ytterbium
c
 70   continue
c
      cc( 1, 1) =         -0.23695500d0
      cc( 2, 1) =         -2.94540310d0
      cc( 3, 1) =        -15.92528720d0
      cc( 4, 1) =        -70.09550470d0
      cc( 5, 1) =       -228.63296500d0
      cc( 6, 1) =        -36.79430390d0
      a ( 1, 1) =          0.20919999d0
      a ( 2, 1) =          0.81500000d0
      a ( 3, 1) =          2.72350001d0
      a ( 4, 1) =         10.49730020d0
      a ( 5, 1) =         46.54999920d0
      a ( 6, 1) =        145.14129600d0
      cc( 1, 2) =          0.03595000d0
      cc( 2, 2) =        114.44560200d0
      cc( 3, 2) =       -183.63880900d0
      cc( 4, 2) =         99.74138640d0
      cc( 5, 2) =         -6.33236122d0
      cc( 6, 2) =         78.99283600d0
      cc( 7, 2) =         42.97451780d0
      cc( 8, 2) =          6.56887198d0
      a ( 1, 2) =          0.14550000d0
      a ( 2, 2) =          0.91399997d0
      a ( 3, 2) =          1.06060004d0
      a ( 4, 2) =          1.37839997d0
      a ( 5, 2) =          3.61319995d0
      a ( 6, 2) =          5.37309981d0
      a ( 7, 2) =         11.38829990d0
      a ( 8, 2) =         24.75140000d0
      cc( 1, 3) =          0.04724300d0
      cc( 2, 3) =        115.76213100d0
      cc( 3, 3) =       -227.08677700d0
      cc( 4, 3) =        213.08325200d0
      cc( 5, 3) =       -122.05501600d0
      cc( 6, 3) =         74.42581940d0
      cc( 7, 3) =         44.72285080d0
      cc( 8, 3) =          5.66238213d0
      a ( 1, 3) =          0.13959999d0
      a ( 2, 3) =          0.73640001d0
      a ( 3, 3) =          0.87620002d0
      a ( 4, 3) =          1.18420005d0
      a ( 5, 3) =          1.63479996d0
      a ( 6, 3) =          2.70499992d0
      a ( 7, 3) =          8.01939964d0
      a ( 8, 3) =         26.23830030d0
      cc( 1, 4) =          0.00776900d0
      cc( 2, 4) =        -25.85315900d0
      cc( 3, 4) =         87.56867220d0
      cc( 4, 4) =        -93.73007200d0
      cc( 5, 4) =        151.94841000d0
      cc( 6, 4) =        -28.00010870d0
      cc( 7, 4) =         60.07494360d0
      cc( 8, 4) =          6.73838902d0
      a ( 1, 4) =          0.11500000d0
      a ( 2, 4) =          1.80429995d0
      a ( 3, 4) =          2.39470005d0
      a ( 4, 4) =          3.80189991d0
      a ( 5, 4) =          6.20260000d0
      a ( 6, 4) =         11.57400040d0
      a ( 7, 4) =         15.53549960d0
      a ( 8, 4) =         46.52880100d0
c
      go to 200
c
c   lutetium
c
 71   continue
c
      cc( 1, 1) =         -0.31990600d0
      cc( 2, 1) =         -3.15664697d0
      cc( 3, 1) =        -16.28186800d0
      cc( 4, 1) =        -69.33083340d0
      cc( 5, 1) =       -226.54129000d0
      cc( 6, 1) =        -38.74679190d0
      a ( 1, 1) =          0.25189999d0
      a ( 2, 1) =          0.93210000d0
      a ( 3, 1) =          2.93959999d0
      a ( 4, 1) =         10.88259980d0
      a ( 5, 1) =         46.68960190d0
      a ( 6, 1) =        154.12980700d0
      cc( 1, 2) =          0.05339400d0
      cc( 2, 2) =        117.31488000d0
      cc( 3, 2) =       -227.90637200d0
      cc( 4, 2) =        207.11132800d0
      cc( 5, 2) =        -99.07708740d0
      cc( 6, 2) =        -36.54268270d0
      cc( 7, 2) =         54.19141770d0
      cc( 8, 2) =          3.64155197d0
      a ( 1, 2) =          0.17810001d0
      a ( 2, 2) =          0.99000001d0
      a ( 3, 2) =          1.18410003d0
      a ( 4, 2) =          1.55540001d0
      a ( 5, 2) =          2.13179994d0
      a ( 6, 2) =          4.83319998d0
      a ( 7, 2) =          3.05800009d0
      a ( 8, 2) =         17.84269910d0
      cc( 1, 3) =         32.67243960d0
      cc( 2, 3) =        -78.01814270d0
      cc( 3, 3) =        100.01547200d0
      cc( 4, 3) =        -69.97026820d0
      cc( 5, 3) =         12.09111500d0
      cc( 6, 3) =         54.91968160d0
      cc( 7, 3) =         41.97269820d0
      cc( 8, 3) =          5.71504688d0
      a ( 1, 3) =          0.60970002d0
      a ( 2, 3) =          0.85030001d0
      a ( 3, 3) =          1.20529997d0
      a ( 4, 3) =          1.77429998d0
      a ( 5, 3) =          2.33669996d0
      a ( 6, 3) =          2.90009999d0
      a ( 7, 3) =          7.69759989d0
      a ( 8, 3) =         21.86129950d0
      cc( 1, 4) =          0.08798900d0
      cc( 2, 4) =        -31.79989050d0
      cc( 3, 4) =         46.60775380d0
      cc( 4, 4) =        -41.07834250d0
      cc( 5, 4) =        150.61451700d0
      cc( 6, 4) =       -159.04402200d0
      cc( 7, 4) =         54.65802380d0
      cc( 8, 4) =          6.87314081d0
      a ( 1, 4) =          0.22400001d0
      a ( 2, 4) =          0.92799997d0
      a ( 3, 4) =          0.99860001d0
      a ( 4, 4) =          1.46200001d0
      a ( 5, 4) =          2.50060010d0
      a ( 6, 4) =          3.25530005d0
      a ( 7, 4) =          4.62179995d0
      a ( 8, 4) =         24.98469930d0
c
 200  n ( 1, 1) =   2
      n ( 2, 1) =   2
      n ( 3, 1) =   2
      n ( 4, 1) =   2
      n ( 5, 1) =   2
      n ( 6, 1) =   1
      nt( 1)    =   6
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      n ( 4, 2) =   2
      n ( 5, 2) =   2
      n ( 6, 2) =   2
      n ( 7, 2) =   1
      n ( 8, 2) =   0
      nt( 2)    =   8
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      n ( 3, 3) =   2
      n ( 4, 3) =   2
      n ( 5, 3) =   2
      n ( 6, 3) =   2
      n ( 7, 3) =   1
      n ( 8, 3) =   0
      nt( 3)    =   8
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      n ( 3, 4) =   2
      n ( 4, 4) =   2
      n ( 5, 4) =   2
      n ( 6, 4) =   2
      n ( 7, 4) =   1
      n ( 8, 4) =   0
      nt( 4)    =   8
      nelcor    =  54
      lmx       =   3
c
      return
 6010 format (1x,'*** data not found for CRENBL ecp type ',a4)
      end
**==crenbl4.f
      subroutine crenbl4 (zinp,lmx,nelcor,nt,n,cc,a)
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
      implicit character * 8 (z)
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
      dimension nt(6), n(10,6), cc(10,6), a(10,6)
c
c     this routine contains a library of ecp data taken from
c     various literature sources.
c     the CRENBL ecps are due to Christiansen etc. al
c     Large Orbital Basis, Small Core Potential
c     shape consistent ECPs
c     L.F. Pacios and P.A> Christiansen J.Chem. Phys. 82 (1985) 2664
c
      dimension ztit(15)
      data maxtyp/15/
      data ztit/
     $ 'hf','ta','w ','re','os','ir','pt','au','hg','tl','pb',
     + 'bi','po','at','rn'/
c
c     this routine contains a library of ecp data taken from
c     various literature sources.
c     3rd transition series sjkb potentials
c     
      do 10 ityp = 1 , maxtyp
         if (zinp.eq.ztit(ityp)) go to 20
 10   continue
 81   write (iwr,6010) zinp
      call caserr2('no library available for requested ecp')
c
 20   go to(72,73,74,75,76,77,78,79,80,81,
     +      82,83,84,85,86), ityp
c
c     hafnium (Hf_CRENBL_ECP)
c
  72  continue
c
      cc( 1, 1) =         -0.33688599d0
      cc( 2, 1) =         -5.75432682d0
      cc( 3, 1) =        -34.38175201d0
      cc( 4, 1) =        -89.85289001d0
      cc( 5, 1) =       -248.74275208d0
      cc( 6, 1) =        -44.91285324d0
      a ( 1, 1) =          0.77789998d0
      a ( 2, 1) =          2.20149994d0
      a ( 3, 1) =          5.79570007d0
      a ( 4, 1) =         14.84070015d0
      a ( 5, 1) =         52.10160065d0
      a ( 6, 1) =        162.34359741d0
      cc( 1, 2) =         -0.64158797d0
      cc( 2, 2) =         55.41609955d0
      cc( 3, 2) =       -176.05253601d0
      cc( 4, 2) =        396.81842041d0
      cc( 5, 2) =       -353.61618042d0
      cc( 6, 2) =        329.38104248d0
      cc( 7, 2) =         27.30104828d0
      cc( 8, 2) =          6.73441696d0
      a ( 1, 2) =          1.08210003d0
      a ( 2, 2) =          2.54130006d0
      a ( 3, 2) =          3.23320007d0
      a ( 4, 2) =          4.75519991d0
      a ( 5, 2) =          7.48320007d0
      a ( 6, 2) =         12.02369976d0
      a ( 7, 2) =         35.47389984d0
      a ( 8, 2) =         30.75819969d0
      cc( 1, 3) =         -0.03051100d0
      cc( 2, 3) =          5.96079922d0
      cc( 3, 3) =       -116.54636383d0
      cc( 4, 3) =        354.66134644d0
      cc( 5, 3) =       -346.44451904d0
      cc( 6, 3) =        286.48968506d0
      cc( 7, 3) =         42.98619461d0
      cc( 8, 3) =          5.40717411d0
      a ( 1, 3) =          0.25839999d0
      a ( 2, 3) =          1.89649999d0
      a ( 3, 3) =          3.19810009d0
      a ( 4, 3) =          4.15119982d0
      a ( 5, 3) =          6.22959995d0
      a ( 6, 3) =          9.52000046d0
      a ( 7, 3) =         24.88139915d0
      a ( 8, 3) =         47.12340164d0
      cc( 1, 4) =         -0.01930300d0
      cc( 2, 4) =         41.46035004d0
      cc( 3, 4) =       -123.71736145d0
      cc( 4, 4) =        228.98324585d0
      cc( 5, 4) =       -253.25607300d0
      cc( 6, 4) =        217.21978760d0
      cc( 7, 4) =         27.25153923d0
      cc( 8, 4) =          7.88987684d0
      a ( 1, 4) =          0.25630000d0
      a ( 2, 4) =          1.74779999d0
      a ( 3, 4) =          2.02530003d0
      a ( 4, 4) =          2.68670011d0
      a ( 5, 4) =          3.86759996d0
      a ( 6, 4) =          5.61810017d0
      a ( 7, 4) =         13.09689999d0
      a ( 8, 4) =         18.43250084d0
      cc( 1, 5) =        -29.21422958d0
      cc( 2, 5) =         84.13317108d0
      cc( 3, 5) =       -109.43501282d0
      cc( 4, 5) =         63.54305649d0
      cc( 5, 5) =         34.68389130d0
      cc( 6, 5) =         94.54371643d0
      cc( 7, 5) =         33.17300034d0
      cc( 8, 5) =          0.93162602d0
      a ( 1, 5) =          0.75269997d0
      a ( 2, 5) =          0.83020002d0
      a ( 3, 5) =          0.98720002d0
      a ( 4, 5) =          1.17200005d0
      a ( 5, 5) =          4.78700018d0
      a ( 6, 5) =         12.48890018d0
      a ( 7, 5) =         33.07509995d0
      a ( 8, 5) =         73.49960327d0
c
      go to 100
c
c     tantalum (Ta_CRENBL_ECP)
c
  73  continue
c
      cc( 1, 1) =         -0.45041099d0
      cc( 2, 1) =         -7.63836813d0
      cc( 3, 1) =        -43.87855530d0
      cc( 4, 1) =        -97.77507019d0
      cc( 5, 1) =       -287.56149292d0
      cc( 6, 1) =        -46.50049973d0
      a ( 1, 1) =          0.97380000d0
      a ( 2, 1) =          2.72210002d0
      a ( 3, 1) =          7.14249992d0
      a ( 4, 1) =         18.33259964d0
      a ( 5, 1) =         63.04660034d0
      a ( 6, 1) =        201.99899292d0
      cc( 1, 2) =        -40.60095978d0
      cc( 2, 2) =        126.83788300d0
      cc( 3, 2) =       -235.24430847d0
      cc( 4, 2) =        445.97167969d0
      cc( 5, 2) =       -389.18530273d0
      cc( 6, 2) =        355.67025757d0
      cc( 7, 2) =         28.51909828d0
      cc( 8, 2) =          6.77081203d0
      a ( 1, 2) =          2.10870004d0
      a ( 2, 2) =          2.46819997d0
      a ( 3, 2) =          3.29310012d0
      a ( 4, 2) =          4.87939978d0
      a ( 5, 2) =          7.62050009d0
      a ( 6, 2) =         12.02050018d0
      a ( 7, 2) =         35.52090073d0
      a ( 8, 2) =         30.82360077d0
      cc( 1, 3) =         -0.09473200d0
      cc( 2, 3) =         59.86339188d0
      cc( 3, 3) =       -192.32298279d0
      cc( 4, 3) =        388.45581055d0
      cc( 5, 3) =       -359.23678589d0
      cc( 6, 3) =        300.74194336d0
      cc( 7, 3) =         34.81594086d0
      cc( 8, 3) =          5.73217821d0
      a ( 1, 3) =          0.57340002d0
      a ( 2, 3) =          2.51550007d0
      a ( 3, 3) =          3.03670001d0
      a ( 4, 3) =          4.18940020d0
      a ( 5, 3) =          6.30660009d0
      a ( 6, 3) =          9.71700001d0
      a ( 7, 3) =         25.63400078d0
      a ( 8, 3) =         34.21419907d0
      cc( 1, 4) =         -0.03177600d0
      cc( 2, 4) =         37.44372940d0
      cc( 3, 4) =       -120.92505646d0
      cc( 4, 4) =        226.60716248d0
      cc( 5, 4) =       -243.61962891d0
      cc( 6, 4) =        210.83656311d0
      cc( 7, 4) =         53.84714890d0
      cc( 8, 4) =          6.98262882d0
      a ( 1, 4) =          0.37210000d0
      a ( 2, 4) =          1.78380001d0
      a ( 3, 4) =          2.09549999d0
      a ( 4, 4) =          2.76589990d0
      a ( 5, 4) =          3.95749998d0
      a ( 6, 4) =          5.68820000d0
      a ( 7, 4) =         15.65960026d0
      a ( 8, 4) =         43.08679962d0
      cc( 1, 5) =         26.30927849d0
      cc( 2, 5) =        -63.06319427d0
      cc( 3, 5) =         55.48675537d0
      cc( 4, 5) =        -12.49588585d0
      cc( 5, 5) =         24.29885674d0
      cc( 6, 5) =         89.15798187d0
      cc( 7, 5) =         32.28345489d0
      cc( 8, 5) =          1.07911205d0
      a ( 1, 5) =          0.76929998d0
      a ( 2, 5) =          0.85350001d0
      a ( 3, 5) =          1.00039995d0
      a ( 4, 5) =          1.13390005d0
      a ( 5, 5) =          3.51349998d0
      a ( 6, 5) =         10.11849976d0
      a ( 7, 5) =         26.51479912d0
      a ( 8, 5) =         58.97980118d0
c
      go to 100
c
c     tungsten (W_CRENBL_ECP)
c
  74  continue
c
      cc( 1, 1) =         -0.54709101d0
      cc( 2, 1) =         -9.24740410d0
      cc( 3, 1) =        -52.48021317d0
      cc( 4, 1) =       -106.44053650d0
      cc( 5, 1) =       -320.76925659d0
      cc( 6, 1) =        -47.99345398d0
      a ( 1, 1) =          1.16069996d0
      a ( 2, 1) =          3.19729996d0
      a ( 3, 1) =          8.38640022d0
      a ( 4, 1) =         22.00390053d0
      a ( 5, 1) =         73.78410339d0
      a ( 6, 1) =        241.43850708d0
      cc( 1, 2) =        -39.29640198d0
      cc( 2, 2) =        124.15521240d0
      cc( 3, 2) =       -228.80580139d0
      cc( 4, 2) =        432.31655884d0
      cc( 5, 2) =       -361.88928223d0
      cc( 6, 2) =        357.38705444d0
      cc( 7, 2) =         30.99152946d0
      cc( 8, 2) =          6.80353212d0
      a ( 1, 2) =          2.27020001d0
      a ( 2, 2) =          2.63599992d0
      a ( 3, 2) =          3.46090007d0
      a ( 4, 2) =          5.08890009d0
      a ( 5, 2) =          7.92360020d0
      a ( 6, 2) =         12.65060043d0
      a ( 7, 2) =         35.13869858d0
      a ( 8, 2) =         32.79159927d0
      cc( 1, 3) =        -31.35513878d0
      cc( 2, 3) =        103.97112274d0
      cc( 3, 3) =       -200.15911865d0
      cc( 4, 3) =        372.58468628d0
      cc( 5, 3) =       -344.03176880d0
      cc( 6, 3) =        300.36007690d0
      cc( 7, 3) =         36.69195175d0
      cc( 8, 3) =          5.78834915d0
      a ( 1, 3) =          1.82889998d0
      a ( 2, 3) =          2.16409993d0
      a ( 3, 3) =          2.87249994d0
      a ( 4, 3) =          4.19990015d0
      a ( 5, 3) =          6.40469980d0
      a ( 6, 3) =          9.77509975d0
      a ( 7, 3) =         25.06100082d0
      a ( 8, 3) =         34.81399918d0
      cc( 1, 4) =        -24.22998619d0
      cc( 2, 4) =         75.06442261d0
      cc( 3, 4) =       -149.58029175d0
      cc( 4, 4) =        267.75427246d0
      cc( 5, 4) =       -291.28921509d0
      cc( 6, 4) =        253.97956848d0
      cc( 7, 4) =         26.81559563d0
      cc( 8, 4) =          7.90078878d0
      a ( 1, 4) =          1.33169997d0
      a ( 2, 4) =          1.55110002d0
      a ( 3, 4) =          2.06049991d0
      a ( 4, 4) =          2.91930008d0
      a ( 5, 4) =          4.37589979d0
      a ( 6, 4) =          6.54969978d0
      a ( 7, 4) =         16.10079956d0
      a ( 8, 4) =         20.45079994d0
      cc( 1, 5) =         37.99953842d0
      cc( 2, 5) =        -73.81892395d0
      cc( 3, 5) =        102.53340149d0
      cc( 4, 5) =        -60.67398453d0
      cc( 5, 5) =         25.93541336d0
      cc( 6, 5) =         93.24456787d0
      cc( 7, 5) =         33.22473145d0
      cc( 8, 5) =          1.05253196d0
      a ( 1, 5) =          0.87849998d0
      a ( 2, 5) =          0.95969999d0
      a ( 3, 5) =          1.22749996d0
      a ( 4, 5) =          1.33510005d0
      a ( 5, 5) =          3.73200011d0
      a ( 6, 5) =         10.76949978d0
      a ( 7, 5) =         28.34440041d0
      a ( 8, 5) =         66.98750305d0
c
      go to 100
c
c     rhenium (Re_CRENBL_ECP)
c
  75  continue
      cc( 1, 1) =         -0.63006300d0
      cc( 2, 1) =        -10.51489067d0
      cc( 3, 1) =        -58.95209122d0
      cc( 4, 1) =       -115.01941681d0
      cc( 5, 1) =       -347.15853882d0
      cc( 6, 1) =        -49.26778412d0
      a ( 1, 1) =          1.34420002d0
      a ( 2, 1) =          3.63630009d0
      a ( 3, 1) =          9.45650005d0
      a ( 4, 1) =         25.35820007d0
      a ( 5, 1) =         83.40769959d0
      a ( 6, 1) =        276.78070068d0
      cc( 1, 2) =        -44.26672363d0
      cc( 2, 2) =        129.19606018d0
      cc( 3, 2) =       -242.92758179d0
      cc( 4, 2) =        501.74584961d0
      cc( 5, 2) =       -445.97048950d0
      cc( 6, 2) =        394.00692749d0
      cc( 7, 2) =         26.49606895d0
      cc( 8, 2) =          6.66976500d0
      a ( 1, 2) =          2.36209989d0
      a ( 2, 2) =          2.72939992d0
      a ( 3, 2) =          3.67540002d0
      a ( 4, 2) =          5.53210020d0
      a ( 5, 2) =          8.98340034d0
      a ( 6, 2) =         15.25590038d0
      a ( 7, 2) =         50.81330109d0
      a ( 8, 2) =         35.46549988d0
      cc( 1, 3) =        -40.03034592d0
      cc( 2, 3) =        126.64355469d0
      cc( 3, 3) =       -228.62715149d0
      cc( 4, 3) =        402.92218018d0
      cc( 5, 3) =       -364.56909180d0
      cc( 6, 3) =        314.11911011d0
      cc( 7, 3) =         41.13771057d0
      cc( 8, 3) =          5.77764893d0
      a ( 1, 3) =          2.03469992d0
      a ( 2, 3) =          2.35809994d0
      a ( 3, 3) =          3.05710006d0
      a ( 4, 3) =          4.39909983d0
      a ( 5, 3) =          6.61850023d0
      a ( 6, 3) =         10.01210022d0
      a ( 7, 3) =         25.28219986d0
      a ( 8, 3) =         38.69940186d0
      cc( 1, 4) =        -19.84150887d0
      cc( 2, 4) =         41.82987595d0
      cc( 3, 4) =       -115.67475128d0
      cc( 4, 4) =        256.13897705d0
      cc( 5, 4) =       -264.05648804d0
      cc( 6, 4) =        231.04666138d0
      cc( 7, 4) =         28.69319344d0
      cc( 8, 4) =          7.87635279d0
      a ( 1, 4) =          1.33980000d0
      a ( 2, 4) =          1.50170004d0
      a ( 3, 4) =          2.26449990d0
      a ( 4, 4) =          3.06690001d0
      a ( 5, 4) =          4.50220013d0
      a ( 6, 4) =          6.76459980d0
      a ( 7, 4) =         15.22159958d0
      a ( 8, 4) =         21.53820038d0
      cc( 1, 5) =        -69.01650238d0
      cc( 2, 5) =        172.23371887d0
      cc( 3, 5) =       -173.36340332d0
      cc( 4, 5) =         97.75544739d0
      cc( 5, 5) =         78.55423737d0
      cc( 6, 5) =        126.52191925d0
      cc( 7, 5) =         37.84556961d0
      cc( 8, 5) =          0.79649299d0
      a ( 1, 5) =          1.48140001d0
      a ( 2, 5) =          1.68589997d0
      a ( 3, 5) =          2.11229992d0
      a ( 4, 5) =          2.64039993d0
      a ( 5, 5) =          9.64130020d0
      a ( 6, 5) =         24.78100014d0
      a ( 7, 5) =         57.28340149d0
      a ( 8, 5) =        142.22329712d0
c
      go to 100
c
c     osmium (Os_CRENBL_ECP)
c
  76  continue
      cc( 1, 1) =         -0.71941298d0
      cc( 2, 1) =        -11.78117943d0
      cc( 3, 1) =        -65.20819855d0
      cc( 4, 1) =       -124.69836426d0
      cc( 5, 1) =       -372.65100098d0
      cc( 6, 1) =        -50.55568314d0
      a ( 1, 1) =          1.53260005d0
      a ( 2, 1) =          4.08620024d0
      a ( 3, 1) =         10.54440022d0
      a ( 4, 1) =         28.97459984d0
      a ( 5, 1) =         93.63200378d0
      a ( 6, 1) =        314.22119141d0
      cc( 1, 2) =        -32.76930618d0
      cc( 2, 2) =        110.55011749d0
      cc( 3, 2) =       -222.00558472d0
      cc( 4, 2) =        473.63671875d0
      cc( 5, 2) =       -411.24047852d0
      cc( 6, 2) =        389.15637207d0
      cc( 7, 2) =         27.66662025d0
      cc( 8, 2) =          6.70742512d0
      a ( 1, 2) =          2.42689991d0
      a ( 2, 2) =          2.84459996d0
      a ( 3, 2) =          3.78419995d0
      a ( 4, 2) =          5.72410011d0
      a ( 5, 2) =          9.26539993d0
      a ( 6, 2) =         15.67440033d0
      a ( 7, 2) =         46.98600006d0
      a ( 8, 2) =         38.68320084d0
      cc( 1, 3) =        -33.74450302d0
      cc( 2, 3) =        107.94266510d0
      cc( 3, 3) =       -207.27871704d0
      cc( 4, 3) =        396.81192017d0
      cc( 5, 3) =       -384.67047119d0
      cc( 6, 3) =        398.63189697d0
      cc( 7, 3) =         26.76161575d0
      cc( 8, 3) =          9.74371719d0
      a ( 1, 3) =          1.98339999d0
      a ( 2, 3) =          2.31750011d0
      a ( 3, 3) =          3.07459998d0
      a ( 4, 3) =          4.54629993d0
      a ( 5, 3) =          7.13530016d0
      a ( 6, 3) =         11.63640022d0
      a ( 7, 3) =         33.21730041d0
      a ( 8, 3) =         28.84840012d0
      cc( 1, 4) =        -24.80944824d0
      cc( 2, 4) =         79.30261230d0
      cc( 3, 4) =       -148.69277954d0
      cc( 4, 4) =        251.07769775d0
      cc( 5, 4) =       -255.17108154d0
      cc( 6, 4) =        224.03025818d0
      cc( 7, 4) =         30.83857727d0
      cc( 8, 4) =          7.83259392d0
      a ( 1, 4) =          1.52030003d0
      a ( 2, 4) =          1.75339997d0
      a ( 3, 4) =          2.25880003d0
      a ( 4, 4) =          3.16960001d0
      a ( 5, 4) =          4.66480017d0
      a ( 6, 4) =          6.84849977d0
      a ( 7, 4) =         14.54829979d0
      a ( 8, 4) =         22.59539986d0
      cc( 1, 5) =         90.39373016d0
      cc( 2, 5) =       -196.83497620d0
      cc( 3, 5) =        215.51077271d0
      cc( 4, 5) =       -153.25781250d0
      cc( 5, 5) =        101.81298828d0
      cc( 6, 5) =        103.32333374d0
      cc( 7, 5) =         29.86117363d0
      cc( 8, 5) =          0.96181101d0
      a ( 1, 5) =          1.94879997d0
      a ( 2, 5) =          2.25889993d0
      a ( 3, 5) =          2.92580009d0
      a ( 4, 5) =          4.02810001d0
      a ( 5, 5) =          5.84280014d0
      a ( 6, 5) =         18.53650093d0
      a ( 7, 5) =         50.93960190d0
      a ( 8, 5) =          4.61829996d0
c
      go to 100
c
c     iridium (Ir_CRENBL_ECP)
c
  77  continue
c
      cc( 1, 1) =         -0.77600998d0
      cc( 2, 1) =        -12.55831528d0
      cc( 3, 1) =        -70.03989410d0
      cc( 4, 1) =       -133.59544373d0
      cc( 5, 1) =       -394.30059814d0
      cc( 6, 1) =        -51.70143127d0
      a ( 1, 1) =          1.70930004d0
      a ( 2, 1) =          4.47060013d0
      a ( 3, 1) =         11.48029995d0
      a ( 4, 1) =         32.14350128d0
      a ( 5, 1) =        102.84519959d0
      a ( 6, 1) =        348.19989014d0
      cc( 1, 2) =        -28.48190117d0
      cc( 2, 2) =         98.51782227d0
      cc( 3, 2) =       -206.84092712d0
      cc( 4, 2) =        457.16946411d0
      cc( 5, 2) =       -384.66342163d0
      cc( 6, 2) =        389.39016724d0
      cc( 7, 2) =         28.70732689d0
      cc( 8, 2) =          6.73599482d0
      a ( 1, 2) =          2.52270007d0
      a ( 2, 2) =          2.95079994d0
      a ( 3, 2) =          3.92389989d0
      a ( 4, 2) =          5.92140007d0
      a ( 5, 2) =          9.53390026d0
      a ( 6, 2) =         16.08679962d0
      a ( 7, 2) =         46.95880127d0
      a ( 8, 2) =         39.62530136d0
      cc( 1, 3) =        -32.93295670d0
      cc( 2, 3) =        108.04138947d0
      cc( 3, 3) =       -213.19577026d0
      cc( 4, 3) =        419.07299805d0
      cc( 5, 3) =       -382.39443970d0
      cc( 6, 3) =        340.84274292d0
      cc( 7, 3) =         51.53657913d0
      cc( 8, 3) =          5.29221916d0
      a ( 1, 3) =          2.18729997d0
      a ( 2, 3) =          2.54640007d0
      a ( 3, 3) =          3.35279989d0
      a ( 4, 3) =          4.92460012d0
      a ( 5, 3) =          7.71029997d0
      a ( 6, 3) =         12.32499981d0
      a ( 7, 3) =         33.08100128d0
      a ( 8, 3) =         63.94940186d0
      cc( 1, 4) =        -32.97980499d0
      cc( 2, 4) =         92.38909149d0
      cc( 3, 4) =       -168.03038025d0
      cc( 4, 4) =        302.14080811d0
      cc( 5, 4) =       -324.87911987d0
      cc( 6, 4) =        295.09466553d0
      cc( 7, 4) =         27.30807304d0
      cc( 8, 4) =          7.89976501d0
      a ( 1, 4) =          1.62119997d0
      a ( 2, 4) =          1.84819996d0
      a ( 3, 4) =          2.43589997d0
      a ( 4, 4) =          3.48830009d0
      a ( 5, 4) =          5.28469992d0
      a ( 6, 4) =          8.07320023d0
      a ( 7, 4) =         20.40609932d0
      a ( 8, 4) =         23.80649948d0
      cc( 1, 5) =        -62.39103699d0
      cc( 2, 5) =        210.25509644d0
      cc( 3, 5) =       -335.25369263d0
      cc( 4, 5) =        374.69042969d0
      cc( 5, 5) =       -257.36740112d0
      cc( 6, 5) =        215.63255310d0
      cc( 7, 5) =         34.77194977d0
      cc( 8, 5) =          4.84190178d0
      a ( 1, 5) =          1.69110000d0
      a ( 2, 5) =          1.98189998d0
      a ( 3, 5) =          2.63700008d0
      a ( 4, 5) =          3.81979990d0
      a ( 5, 5) =          5.81199980d0
      a ( 6, 5) =          9.10079956d0
      a ( 7, 5) =         21.42609978d0
      a ( 8, 5) =         35.98310089d0
c
      go to 100
c
c     platium (Pt_CRENBL_ECP)
c
  78  continue
c
      cc( 1, 1) =         -0.81254202d0
      cc( 2, 1) =        -13.05292892d0
      cc( 3, 1) =        -73.72098541d0
      cc( 4, 1) =       -141.17144775d0
      cc( 5, 1) =       -412.57824707d0
      cc( 6, 1) =        -52.72452545d0
      a ( 1, 1) =          1.86479998d0
      a ( 2, 1) =          4.83080006d0
      a ( 3, 1) =         12.30449963d0
      a ( 4, 1) =         34.86819839d0
      a ( 5, 1) =        110.94640350d0
      a ( 6, 1) =        378.44430542d0
      cc( 1, 2) =        -26.62452316d0
      cc( 2, 2) =         90.01412201d0
      cc( 3, 2) =       -193.86340332d0
      cc( 4, 2) =        432.62850952d0
      cc( 5, 2) =       -343.03793335d0
      cc( 6, 2) =        387.98568726d0
      cc( 7, 2) =         29.76094818d0
      cc( 8, 2) =          6.77873421d0
      a ( 1, 2) =          2.60910010d0
      a ( 2, 2) =          3.02920008d0
      a ( 3, 2) =          4.04290009d0
      a ( 4, 2) =          6.07630014d0
      a ( 5, 2) =          9.83839989d0
      a ( 6, 2) =         16.70409966d0
      a ( 7, 2) =         47.62310028d0
      a ( 8, 2) =         40.06890106d0
      cc( 1, 3) =        -30.37420082d0
      cc( 2, 3) =        102.49530029d0
      cc( 3, 3) =       -205.28015137d0
      cc( 4, 3) =        412.87677002d0
      cc( 5, 3) =       -373.50769043d0
      cc( 6, 3) =        355.23352051d0
      cc( 7, 3) =         34.39381027d0
      cc( 8, 3) =          5.90135813d0
      a ( 1, 3) =          2.31279993d0
      a ( 2, 3) =          2.68169999d0
      a ( 3, 3) =          3.50119996d0
      a ( 4, 3) =          5.13490009d0
      a ( 5, 3) =          8.00559998d0
      a ( 6, 3) =         12.87129974d0
      a ( 7, 3) =         35.53670120d0
      a ( 8, 3) =         36.64830017d0
      cc( 1, 4) =         -3.02787995d0
      cc( 2, 4) =         64.21149445d0
      cc( 3, 4) =       -146.10299683d0
      cc( 4, 4) =        266.61413574d0
      cc( 5, 4) =       -301.47027588d0
      cc( 6, 4) =        278.51043701d0
      cc( 7, 4) =         29.71852303d0
      cc( 8, 4) =          7.86891603d0
      a ( 1, 4) =          1.38300002d0
      a ( 2, 4) =          2.07049990d0
      a ( 3, 4) =          2.49049997d0
      a ( 4, 4) =          3.67370009d0
      a ( 5, 4) =          5.48519993d0
      a ( 6, 4) =          8.20989990d0
      a ( 7, 4) =         18.57180023d0
      a ( 8, 4) =         25.42099953d0
      cc( 1, 5) =        -63.38943481d0
      cc( 2, 5) =        215.72419739d0
      cc( 3, 5) =       -333.95880127d0
      cc( 4, 5) =        368.06796265d0
      cc( 5, 5) =       -242.79386902d0
      cc( 6, 5) =        233.02484131d0
      cc( 7, 5) =         43.59355164d0
      cc( 8, 5) =          4.53264999d0
      a ( 1, 5) =          1.86710000d0
      a ( 2, 5) =          2.20180011d0
      a ( 3, 5) =          2.95210004d0
      a ( 4, 5) =          4.37239981d0
      a ( 5, 5) =          6.78480005d0
      a ( 6, 5) =         11.02280045d0
      a ( 7, 5) =         27.73590088d0
      a ( 8, 5) =         54.84930038d0
c
      go to 100
c
c     gold (Au_CRENBL_ECP)
c
  79  continue
c
      cc( 1, 1) =         -0.90820301d0
      cc( 2, 1) =        -13.97286034d0
      cc( 3, 1) =        -78.04593658d0
      cc( 4, 1) =       -150.67013550d0
      cc( 5, 1) =       -433.79324341d0
      cc( 6, 1) =        -53.92461777d0
      a ( 1, 1) =          2.07069993d0
      a ( 2, 1) =          5.28779984d0
      a ( 3, 1) =         13.28610039d0
      a ( 4, 1) =         38.22499847d0
      a ( 5, 1) =        120.82949829d0
      a ( 6, 1) =        415.27328491d0
      cc( 1, 2) =        -45.11843872d0
      cc( 2, 2) =        109.24782562d0
      cc( 3, 2) =       -184.75865173d0
      cc( 4, 2) =        436.13348389d0
      cc( 5, 2) =       -356.46279907d0
      cc( 6, 2) =        411.84240723d0
      cc( 7, 2) =         31.95258903d0
      cc( 8, 2) =          6.79167795d0
      a ( 1, 2) =          2.81620002d0
      a ( 2, 2) =          3.12030005d0
      a ( 3, 2) =          4.16239977d0
      a ( 4, 2) =          6.38819981d0
      a ( 5, 2) =         10.11999989d0
      a ( 6, 2) =         17.04649925d0
      a ( 7, 2) =         50.50329971d0
      a ( 8, 2) =         39.28239822d0
      cc( 1, 3) =        -30.62415695d0
      cc( 2, 3) =        101.11966705d0
      cc( 3, 3) =       -210.51446533d0
      cc( 4, 3) =        435.37313843d0
      cc( 5, 3) =       -396.14450073d0
      cc( 6, 3) =        368.87313843d0
      cc( 7, 3) =         31.59367371d0
      cc( 8, 3) =          5.83880520d0
      a ( 1, 3) =          2.41569996d0
      a ( 2, 3) =          2.79209995d0
      a ( 3, 3) =          3.67700005d0
      a ( 4, 3) =          5.37360001d0
      a ( 5, 3) =          8.36999989d0
      a ( 6, 3) =         13.32330036d0
      a ( 7, 3) =         37.15980148d0
      a ( 8, 3) =         37.42219925d0
      cc( 1, 4) =        -34.10911942d0
      cc( 2, 4) =        100.58750153d0
      cc( 3, 4) =       -172.30566406d0
      cc( 4, 4) =        290.82238770d0
      cc( 5, 4) =       -292.96429443d0
      cc( 6, 4) =        275.67791748d0
      cc( 7, 4) =         30.37005806d0
      cc( 8, 4) =          7.86203289d0
      a ( 1, 4) =          1.87390006d0
      a ( 2, 4) =          2.11870003d0
      a ( 3, 4) =          2.69359994d0
      a ( 4, 4) =          3.81640005d0
      a ( 5, 4) =          5.69740009d0
      a ( 6, 4) =          8.62959957d0
      a ( 7, 4) =         19.68020058d0
      a ( 8, 4) =         25.98660088d0
      cc( 1, 5) =         40.06137466d0
      cc( 2, 5) =       -128.38630676d0
      cc( 3, 5) =        247.77139282d0
      cc( 4, 5) =       -342.54879761d0
      cc( 5, 5) =        359.27548218d0
      cc( 6, 5) =       -182.27494812d0
      cc( 7, 5) =         53.49130630d0
      cc( 8, 5) =          4.28036213d0
      a ( 1, 5) =          1.43149996d0
      a ( 2, 5) =          1.65450001d0
      a ( 3, 5) =          2.16300011d0
      a ( 4, 5) =          3.08719993d0
      a ( 5, 5) =          4.65740013d0
      a ( 6, 5) =          6.90420008d0
      a ( 7, 5) =         10.14780045d0
      a ( 8, 5) =         47.02579880d0
c
      go to 100
c
c     mercury (Hg_CRENBL_ECP)
c
  80  continue
c
      cc( 1, 1) =         -0.94575602d0
      cc( 2, 1) =        -14.01616859d0
      cc( 3, 1) =        -81.47489166d0
      cc( 4, 1) =       -158.27243042d0
      cc( 5, 1) =       -451.69934082d0
      cc( 6, 1) =        -54.98525620d0
      a ( 1, 1) =          2.26189995d0
      a ( 2, 1) =          5.61199999d0
      a ( 3, 1) =         14.06490040d0
      a ( 4, 1) =         40.92689896d0
      a ( 5, 1) =        129.16789246d0
      a ( 6, 1) =        447.18978882d0
      cc( 1, 2) =        -25.16012001d0
      cc( 2, 2) =         91.19132233d0
      cc( 3, 2) =       -201.98194885d0
      cc( 4, 2) =        454.72329712d0
      cc( 5, 2) =       -358.64822388d0
      cc( 6, 2) =        429.54168701d0
      cc( 7, 2) =         35.77702332d0
      cc( 8, 2) =          6.85280085d0
      a ( 1, 2) =          2.83570004d0
      a ( 2, 2) =          3.29949999d0
      a ( 3, 2) =          4.36730003d0
      a ( 4, 2) =          6.54479980d0
      a ( 5, 2) =         10.36839962d0
      a ( 6, 2) =         17.26539993d0
      a ( 7, 2) =         46.42380142d0
      a ( 8, 2) =         41.91410065d0
      cc( 1, 3) =        -29.68764877d0
      cc( 2, 3) =        101.14541626d0
      cc( 3, 3) =       -208.31994629d0
      cc( 4, 3) =        427.21362305d0
      cc( 5, 3) =       -379.29605103d0
      cc( 6, 3) =        371.70361328d0
      cc( 7, 3) =         33.64100647d0
      cc( 8, 3) =          5.87053490d0
      a ( 1, 3) =          2.49449992d0
      a ( 2, 3) =          2.88730001d0
      a ( 3, 3) =          3.77839994d0
      a ( 4, 3) =          5.53910017d0
      a ( 5, 3) =          8.60260010d0
      a ( 6, 3) =         13.72500038d0
      a ( 7, 3) =         37.17200089d0
      a ( 8, 3) =         38.38690186d0
      cc( 1, 4) =        -29.97778702d0
      cc( 2, 4) =         93.63288879d0
      cc( 3, 4) =       -181.20259094d0
      cc( 4, 4) =        337.58990479d0
      cc( 5, 4) =       -357.07931519d0
      cc( 6, 4) =        342.55819702d0
      cc( 7, 4) =         27.91793060d0
      cc( 8, 4) =          7.87457609d0
      a ( 1, 4) =          1.86319995d0
      a ( 2, 4) =          2.15249991d0
      a ( 3, 4) =          2.83509994d0
      a ( 4, 4) =          4.11339998d0
      a ( 5, 4) =          6.33160019d0
      a ( 6, 4) =         10.02639961d0
      a ( 7, 4) =         26.72430038d0
      a ( 8, 4) =         28.09569931d0
      cc( 1, 5) =         41.19937897d0
      cc( 2, 5) =       -132.28570557d0
      cc( 3, 5) =        259.85302734d0
      cc( 4, 5) =       -360.87664795d0
      cc( 5, 5) =        376.67376709d0
      cc( 6, 5) =       -186.24418640d0
      cc( 7, 5) =         56.03879929d0
      cc( 8, 5) =          4.22023201d0
      a ( 1, 5) =          1.55799997d0
      a ( 2, 5) =          1.81250000d0
      a ( 3, 5) =          2.40109992d0
      a ( 4, 5) =          3.48460007d0
      a ( 5, 5) =          5.36180019d0
      a ( 6, 5) =          8.20409966d0
      a ( 7, 5) =         12.02680016d0
      a ( 8, 5) =         55.55799866d0
c
 100  continue
      nelcor    =  60
c
      go to 30
c
c     lead  (Pb_CRENBL_ECP)
c
  82  continue
c
      cc( 1, 1) =         -1.41434395d0
      cc( 2, 1) =        -13.28021622d0
      cc( 3, 1) =        -54.13752747d0
      cc( 4, 1) =       -138.12661743d0
      cc( 5, 1) =       -376.41168213d0
      cc( 6, 1) =        -54.85811996d0
      a ( 1, 1) =          1.30009997d0
      a ( 2, 1) =          3.18720007d0
      a ( 3, 1) =          9.54880047d0
      a ( 4, 1) =         23.68409920d0
      a ( 5, 1) =         81.04219818d0
      a ( 6, 1) =        260.18218994d0
      cc( 1, 2) =         45.00151825d0
      cc( 2, 2) =       -142.81721497d0
      cc( 3, 2) =        280.33993530d0
      cc( 4, 2) =       -305.93109131d0
      cc( 5, 2) =        275.88833618d0
      cc( 6, 2) =       -126.19586945d0
      cc( 7, 2) =         58.64039612d0
      cc( 8, 2) =          6.27654600d0
      a ( 1, 2) =          1.05130005d0
      a ( 2, 2) =          1.21360004d0
      a ( 3, 2) =          1.59290004d0
      a ( 4, 2) =          2.27550006d0
      a ( 5, 2) =          3.41589999d0
      a ( 6, 2) =          5.04549980d0
      a ( 7, 2) =          8.10179996d0
      a ( 8, 2) =         29.92939949d0
      cc( 1, 3) =          0.19562900d0
      cc( 2, 3) =        -60.52689743d0
      cc( 3, 3) =        201.29090881d0
      cc( 4, 3) =       -254.73730469d0
      cc( 5, 3) =        231.82109070d0
      cc( 6, 3) =       -120.33583069d0
      cc( 7, 3) =         43.35644531d0
      cc( 8, 3) =          5.93558693d0
      a ( 1, 3) =          0.34459999d0
      a ( 2, 3) =          1.01220000d0
      a ( 3, 3) =          1.21469998d0
      a ( 4, 3) =          1.65750003d0
      a ( 5, 3) =          2.40109992d0
      a ( 6, 3) =          3.46830010d0
      a ( 7, 3) =          4.71899986d0
      a ( 8, 3) =         14.24479961d0
      cc( 1, 4) =         -1.35663700d0
      cc( 2, 4) =         49.74739456d0
      cc( 3, 4) =       -160.28633118d0
      cc( 4, 4) =        325.69772339d0
      cc( 5, 4) =       -339.07305908d0
      cc( 6, 4) =        199.77455139d0
      cc( 7, 4) =         40.33429718d0
      cc( 8, 4) =          7.63998222d0
      a ( 1, 4) =          0.88599998d0
      a ( 2, 4) =          2.42969990d0
      a ( 3, 4) =          3.09030008d0
      a ( 4, 4) =          4.45060015d0
      a ( 5, 4) =          6.85960007d0
      a ( 6, 4) =         10.51179981d0
      a ( 7, 4) =          8.92780018d0
      a ( 8, 4) =         30.41259956d0
      cc( 1, 5) =         37.60541153d0
      cc( 2, 5) =        -24.43016624d0
      cc( 3, 5) =        -38.51657867d0
      cc( 4, 5) =         60.97928619d0
      cc( 5, 5) =        -10.03380680d0
      cc( 6, 5) =         52.05616760d0
      cc( 7, 5) =         39.91782379d0
      cc( 8, 5) =          5.47486877d0
      a ( 1, 5) =          0.74299997d0
      a ( 2, 5) =          0.93489999d0
      a ( 3, 5) =          0.93709999d0
      a ( 4, 5) =          1.35049999d0
      a ( 5, 5) =          2.45320010d0
      a ( 6, 5) =          3.99049997d0
      a ( 7, 5) =          8.44530010d0
      a ( 8, 5) =         19.77610016d0
c
      go to 40
c
c     bismuth (Bi_CRENBL_ECP)
c
  83  continue
c
      cc( 1, 1) =         -1.29633200d0
      cc( 2, 1) =        -13.97006321d0
      cc( 3, 1) =        -57.20818710d0
      cc( 4, 1) =       -141.80029297d0
      cc( 5, 1) =       -390.70785522d0
      cc( 6, 1) =        -55.54542160d0
      a ( 1, 1) =          1.38310003d0
      a ( 2, 1) =          3.38969994d0
      a ( 3, 1) =         10.19880009d0
      a ( 4, 1) =         25.27059937d0
      a ( 5, 1) =         86.08190155d0
      a ( 6, 1) =        277.96438599d0
      cc( 1, 2) =         55.90422821d0
      cc( 2, 2) =       -171.26004028d0
      cc( 3, 2) =        297.29568481d0
      cc( 4, 2) =       -320.79031372d0
      cc( 5, 2) =        287.28314209d0
      cc( 6, 2) =       -104.22640228d0
      cc( 7, 2) =         56.34839249d0
      cc( 8, 2) =          6.48098993d0
      a ( 1, 2) =          1.12539995d0
      a ( 2, 2) =          1.28199995d0
      a ( 3, 2) =          1.64119995d0
      a ( 4, 2) =          2.34890008d0
      a ( 5, 2) =          3.33470011d0
      a ( 6, 2) =          4.86899996d0
      a ( 7, 2) =          8.19460011d0
      a ( 8, 2) =         25.78339958d0
      cc( 1, 3) =          0.03434800d0
      cc( 2, 3) =        -56.08747101d0
      cc( 3, 3) =        164.21267700d0
      cc( 4, 3) =       -197.05456543d0
      cc( 5, 3) =        155.65385437d0
      cc( 6, 3) =        -55.68730927d0
      cc( 7, 3) =         50.05516815d0
      cc( 8, 3) =          5.74898911d0
      a ( 1, 3) =          0.20800000d0
      a ( 2, 3) =          1.11570001d0
      a ( 3, 3) =          1.28740001d0
      a ( 4, 3) =          1.77660000d0
      a ( 5, 3) =          2.40019989d0
      a ( 6, 3) =          4.85920000d0
      a ( 7, 3) =          6.06729984d0
      a ( 8, 3) =         20.30240059d0
      cc( 1, 4) =        -49.14233017d0
      cc( 2, 4) =        117.17475128d0
      cc( 3, 4) =       -186.97654724d0
      cc( 4, 4) =        353.28527832d0
      cc( 5, 4) =       -373.14822388d0
      cc( 6, 4) =        379.30261230d0
      cc( 7, 4) =         29.70112038d0
      cc( 8, 4) =          7.84041214d0
      a ( 1, 4) =          1.70360005d0
      a ( 2, 4) =          2.04159999d0
      a ( 3, 4) =          2.95609999d0
      a ( 4, 4) =          4.60300016d0
      a ( 5, 4) =          7.35540009d0
      a ( 6, 4) =         12.00380039d0
      a ( 7, 4) =         30.31959915d0
      a ( 8, 4) =         34.39179993d0
      cc( 1, 5) =        104.19187927d0
      cc( 2, 5) =       -218.21029663d0
      cc( 3, 5) =        220.57954407d0
      cc( 4, 5) =       -170.08029175d0
      cc( 5, 5) =        108.73613739d0
      cc( 6, 5) =         46.86751938d0
      cc( 7, 5) =         42.52275085d0
      cc( 8, 5) =          5.20974779d0
      a ( 1, 5) =          0.59480000d0
      a ( 2, 5) =          0.69040000d0
      a ( 3, 5) =          0.89520001d0
      a ( 4, 5) =          1.24119997d0
      a ( 5, 5) =          1.67270005d0
      a ( 6, 5) =          4.99800014d0
      a ( 7, 5) =          9.57050037d0
      a ( 8, 5) =         20.95170021d0
c
      go to 40
c
c     Po (Po_CRENBL_ECP)
c
  84  continue
c
      cc( 1, 1) =         -1.45787299d0
      cc( 2, 1) =        -14.71360588d0
      cc( 3, 1) =        -62.38053894d0
      cc( 4, 1) =       -147.98703003d0
      cc( 5, 1) =       -416.04962158d0
      cc( 6, 1) =        -56.73604202d0
      a ( 1, 1) =          1.51339996d0
      a ( 2, 1) =          3.67219996d0
      a ( 3, 1) =         11.04430008d0
      a ( 4, 1) =         27.63310051d0
      a ( 5, 1) =         94.03929901d0
      a ( 6, 1) =        308.29211426d0
      cc( 1, 2) =        -30.03038788d0
      cc( 2, 2) =         96.16854858d0
      cc( 3, 2) =       -181.03051758d0
      cc( 4, 2) =        305.07928467d0
      cc( 5, 2) =       -303.26153564d0
      cc( 6, 2) =        193.48907471d0
      cc( 7, 2) =         59.07072830d0
      cc( 8, 2) =          6.37750721d0
      a ( 1, 2) =          0.83230001d0
      a ( 2, 2) =          0.94760001d0
      a ( 3, 2) =          1.20280004d0
      a ( 4, 2) =          1.66719997d0
      a ( 5, 2) =          2.41849995d0
      a ( 6, 2) =          3.38520002d0
      a ( 7, 2) =         11.98690033d0
      a ( 8, 2) =         38.50559998d0
      cc( 1, 3) =        -22.48581886d0
      cc( 2, 3) =         71.35298157d0
      cc( 3, 3) =       -130.56465149d0
      cc( 4, 3) =        214.88987732d0
      cc( 5, 3) =       -221.69078064d0
      cc( 6, 3) =        143.38671875d0
      cc( 7, 3) =         45.49112320d0
      cc( 8, 3) =          6.19466209d0
      a ( 1, 3) =          0.65499997d0
      a ( 2, 3) =          0.73460001d0
      a ( 3, 3) =          0.90689999d0
      a ( 4, 3) =          1.21190000d0
      a ( 5, 3) =          1.68499994d0
      a ( 6, 3) =          2.27480006d0
      a ( 7, 3) =          8.02330017d0
      a ( 8, 3) =         19.25239944d0
      cc( 1, 4) =         -7.33783579d0
      cc( 2, 4) =         25.26849174d0
      cc( 3, 4) =        -52.30345917d0
      cc( 4, 4) =         87.85886383d0
      cc( 5, 4) =       -140.93832397d0
      cc( 6, 4) =        145.23880005d0
      cc( 7, 4) =         28.62387466d0
      cc( 8, 4) =          7.73479319d0
      a ( 1, 4) =          0.76300001d0
      a ( 2, 4) =          0.88559997d0
      a ( 3, 4) =          1.14789999d0
      a ( 4, 4) =          1.61360002d0
      a ( 5, 4) =          2.38429999d0
      a ( 6, 4) =          3.38560009d0
      a ( 7, 4) =         15.32299995d0
      a ( 8, 4) =         13.38679981d0
      cc( 1, 5) =          9.20586109d0
      cc( 2, 5) =        -28.85299873d0
      cc( 3, 5) =         52.88306808d0
      cc( 4, 5) =        -85.14630890d0
      cc( 5, 5) =        124.28939819d0
      cc( 6, 5) =       -116.73908997d0
      cc( 7, 5) =         38.18548203d0
      cc( 8, 5) =          5.00899792d0
      a ( 1, 5) =          0.69440001d0
      a ( 2, 5) =          0.78090000d0
      a ( 3, 5) =          0.96869999d0
      a ( 4, 5) =          1.28680003d0
      a ( 5, 5) =          1.77709997d0
      a ( 6, 5) =          2.38439989d0
      a ( 7, 5) =          2.82629991d0
      a ( 8, 5) =         14.08730030d0
c
      go to 40
c
c     at (At_CRENBL_ECP)
c
  85  continue
c
      cc( 1, 1) =         -0.31176600d0
      cc( 2, 1) =         -0.73071098d0
      cc( 3, 1) =         -2.97282290d0
      cc( 4, 1) =        -32.68309402d0
      cc( 5, 1) =       -164.79605103d0
      cc( 6, 1) =        -51.18941498d0
      a ( 1, 1) =          0.03180000d0
      a ( 2, 1) =          0.18099999d0
      a ( 3, 1) =          1.08060002d0
      a ( 4, 1) =          4.99060011d0
      a ( 5, 1) =         20.39839935d0
      a ( 6, 1) =         81.91079712d0
      cc( 1, 2) =          0.00341100d0
      cc( 2, 2) =          0.30304900d0
      cc( 3, 2) =          0.69335699d0
      cc( 4, 2) =         41.97618484d0
      cc( 5, 2) =         58.09056473d0
      cc( 6, 2) =          3.23042393d0
      cc( 7, 2) =         38.14316940d0
      cc( 8, 2) =          7.18064022d0
      a ( 1, 2) =          0.00200000d0
      a ( 2, 2) =          0.03100000d0
      a ( 3, 2) =          0.17010000d0
      a ( 4, 2) =          1.76409996d0
      a ( 5, 2) =          5.54829979d0
      a ( 6, 2) =          9.53069973d0
      a ( 7, 2) =         10.47630024d0
      a ( 8, 2) =         15.64519978d0
      cc( 1, 3) =          0.31012201d0
      cc( 2, 3) =          0.71845901d0
      cc( 3, 3) =         55.89923096d0
      cc( 4, 3) =        -49.86803818d0
      cc( 5, 3) =         75.92673492d0
      cc( 6, 3) =         -0.51495600d0
      cc( 7, 3) =         51.50535202d0
      cc( 8, 3) =          5.78414392d0
      a ( 1, 3) =          0.03170000d0
      a ( 2, 3) =          0.17670000d0
      a ( 3, 3) =          1.51250005d0
      a ( 4, 3) =          2.18499994d0
      a ( 5, 3) =          3.99130011d0
      a ( 6, 3) =          4.55789995d0
      a ( 7, 3) =         11.17259979d0
      a ( 8, 3) =         31.24390030d0
      cc( 1, 4) =          0.31432700d0
      cc( 2, 4) =          0.79817301d0
      cc( 3, 4) =         12.21004009d0
      cc( 4, 4) =        -38.19544220d0
      cc( 5, 4) =        113.45125580d0
      cc( 6, 4) =        -32.16876984d0
      cc( 7, 4) =         69.49798584d0
      cc( 8, 4) =          7.07968616d0
      a ( 1, 4) =          0.03190000d0
      a ( 2, 4) =          0.19059999d0
      a ( 3, 4) =          1.82609999d0
      a ( 4, 4) =          2.64429998d0
      a ( 5, 4) =          4.57380009d0
      a ( 6, 4) =          8.54640007d0
      a ( 7, 4) =         16.26619911d0
      a ( 8, 4) =         67.92829895d0
      cc( 1, 5) =          0.36805499d0
      cc( 2, 5) =          1.64947903d0
      cc( 3, 5) =         14.43013954d0
      cc( 4, 5) =         17.18401146d0
      cc( 5, 5) =         94.00454712d0
      cc( 6, 5) =       -110.99547577d0
      cc( 7, 5) =         81.44644165d0
      cc( 8, 5) =          9.42089081d0
      a ( 1, 5) =          0.03780000d0
      a ( 2, 5) =          0.18050000d0
      a ( 3, 5) =          0.92019999d0
      a ( 4, 5) =          2.18230009d0
      a ( 5, 5) =          3.96090007d0
      a ( 6, 5) =          5.98220015d0
      a ( 7, 5) =          5.49179983d0
      a ( 8, 5) =          0.08980000d0
c
      go to 40
c
c     radon (Rn_CRENBL_ECP)
c
  86  continue
c
      cc( 1, 1) =         -1.11443496d0
      cc( 2, 1) =        -17.17100143d0
      cc( 3, 1) =        -74.57586670d0
      cc( 4, 1) =       -161.65777588d0
      cc( 5, 1) =       -466.10955811d0
      cc( 6, 1) =        -59.22294235d0
      a ( 1, 1) =          1.65950000d0
      a ( 2, 1) =          4.21479988d0
      a ( 3, 1) =         13.02740002d0
      a ( 4, 1) =         33.21269989d0
      a ( 5, 1) =        111.38510132d0
      a ( 6, 1) =        373.93679810d0
      cc( 1, 2) =          0.00817300d0
      cc( 2, 2) =         23.55170059d0
      cc( 3, 2) =        146.60008240d0
      cc( 4, 2) =       -114.77574158d0
      cc( 5, 2) =        -23.55697250d0
      cc( 6, 2) =         74.31896973d0
      cc( 7, 2) =         31.09328651d0
      cc( 8, 2) =          7.87248802d0
      a ( 1, 2) =          0.02280000d0
      a ( 2, 2) =          1.46410000d0
      a ( 3, 2) =          1.70050001d0
      a ( 4, 2) =          1.52970004d0
      a ( 5, 2) =          3.63899994d0
      a ( 6, 2) =          8.39480019d0
      a ( 7, 2) =          8.56000042d0
      a ( 8, 2) =         20.74449921d0
      cc( 1, 3) =          0.05306900d0
      cc( 2, 3) =        -54.25381851d0
      cc( 3, 3) =        240.12178040d0
      cc( 4, 3) =       -300.78848267d0
      cc( 5, 3) =        261.49554443d0
      cc( 6, 3) =       -114.78035736d0
      cc( 7, 3) =         52.48627472d0
      cc( 8, 3) =          5.68981314d0
      a ( 1, 3) =          0.26650000d0
      a ( 2, 3) =          1.38469994d0
      a ( 3, 3) =          1.69739997d0
      a ( 4, 3) =          2.25480008d0
      a ( 5, 3) =          3.34669995d0
      a ( 6, 3) =          4.87270021d0
      a ( 7, 3) =          8.12670040d0
      a ( 8, 3) =         26.43040085d0
      cc( 1, 4) =        -56.26292419d0
      cc( 2, 4) =        148.12635803d0
      cc( 3, 4) =       -236.51480103d0
      cc( 4, 4) =        425.46905518d0
      cc( 5, 4) =       -427.86010742d0
      cc( 6, 4) =        468.46737671d0
      cc( 7, 4) =         29.30323219d0
      cc( 8, 4) =          7.79323196d0
      a ( 1, 4) =          2.29049993d0
      a ( 2, 4) =          2.72210002d0
      a ( 3, 4) =          3.67820001d0
      a ( 4, 4) =          5.54939985d0
      a ( 5, 4) =          8.97550011d0
      a ( 6, 4) =         15.54349995d0
      a ( 7, 4) =         47.95439911d0
      a ( 8, 4) =         40.16630173d0
      cc( 1, 5) =          2.29281402d0
      cc( 2, 5) =         13.14621830d0
      cc( 3, 5) =          0.83706897d0
      cc( 4, 5) =         49.71719360d0
      cc( 5, 5) =          4.64372921d0
      cc( 6, 5) =        147.32864380d0
      cc( 7, 5) =         99.77916718d0
      cc( 8, 5) =          7.89272690d0
      a ( 1, 5) =          0.15279999d0
      a ( 2, 5) =          0.58080000d0
      a ( 3, 5) =          1.35189998d0
      a ( 4, 5) =          1.98010004d0
      a ( 5, 5) =          2.03579998d0
      a ( 6, 5) =          4.89839983d0
      a ( 7, 5) =         10.49590015d0
      a ( 8, 5) =         22.45910072d0
c
 40   nelcor    =  68
 30   continue
      n ( 1, 1) =   2
      n ( 2, 1) =   2
      n ( 3, 1) =   2
      n ( 4, 1) =   2
      n ( 5, 1) =   2
      n ( 6, 1) =   1
      nt( 1)    =   6
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      n ( 4, 2) =   2
      n ( 5, 2) =   2
      n ( 6, 2) =   2
      n ( 7, 2) =   1
      n ( 8, 2) =   0
      nt( 2)    =   8
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      n ( 3, 3) =   2
      n ( 4, 3) =   2
      n ( 5, 3) =   2
      n ( 6, 3) =   2
      n ( 7, 3) =   1
      n ( 8, 3) =   0
      nt( 3)    =   8
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      n ( 3, 4) =   2
      n ( 4, 4) =   2
      n ( 5, 4) =   2
      n ( 6, 4) =   2
      n ( 7, 4) =   1
      n ( 8, 4) =   0
      nt( 4)    =   8
      n ( 1, 5) =   2
      n ( 2, 5) =   2
      n ( 3, 5) =   2
      n ( 4, 5) =   2
      n ( 5, 5) =   2
      n ( 6, 5) =   2
      n ( 7, 5) =   1
      n ( 8, 5) =   0
      nt( 5)    =   8
      lmx       =   4
c
      return
c
 6010 format (1x,'*** data not found for CRENBL ecp type ',a4)
      end
**==crenbl5.f
      subroutine crenbl5 (zinp,lmx,nelcor,nt,n,cc,a)
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
      implicit character * 8 (z)
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
      dimension nt(6), n(10,6), cc(10,6), a(10,6)
c
c     this routine contains a library of ecp data taken from
c     various literature sources.
c     the CRENBL ecps are due to Christiansen etc. al
c     Large Orbital Basis, Small Core Potential
c     shape consistent ECPs
c     L.F. Pacios and P.A> Christiansen J.Chem. Phys. 82 (1985) 2664
c
      dimension ztit(17)
      data maxtyp/17/
      data ztit/
     +  'fr','ra','ac','th','pa','u ','np','pu','am','cm',
     $  'bk','cf','es','fm','md','no','lw'   /
c
c     this routine contains a library of ecp data taken from
c     various literature sources.
c     3rd transition series sjkb potentials
c     
      do 10 ityp = 1 , maxtyp
         if (zinp.eq.ztit(ityp)) go to 20
 10   continue
      write (iwr,6010) zinp
      call caserr2('no library available for requested ecp')
c
 20   go to(87,88,89,90,91,92,93,94,95,96,97,98,99,100,
     +     101,102,103), ityp
c
c     fr (Fr_CRENBL_ECP)
c
  87  continue
c
      cc( 1, 1) =         -0.89879203d0
      cc( 2, 1) =        -10.09708691d0
      cc( 3, 1) =        -31.60297394d0
      cc( 4, 1) =       -120.95409393d0
      cc( 5, 1) =       -265.95950317d0
      cc( 6, 1) =        -57.15474319d0
      a ( 1, 1) =          0.88870001d0
      a ( 2, 1) =          2.02189994d0
      a ( 3, 1) =          4.90730000d0
      a ( 4, 1) =         15.12580013d0
      a ( 5, 1) =         45.64889908d0
      a ( 6, 1) =        138.01640320d0
      cc( 1, 2) =        -33.04865646d0
      cc( 2, 2) =         99.33948517d0
      cc( 3, 2) =       -183.16125488d0
      cc( 4, 2) =        331.93963623d0
      cc( 5, 2) =       -353.79605103d0
      cc( 6, 2) =        270.63644409d0
      cc( 7, 2) =         30.28211212d0
      cc( 8, 2) =         11.24650669d0
      a ( 1, 2) =          0.91130000d0
      a ( 2, 2) =          1.05540001d0
      a ( 3, 2) =          1.38150001d0
      a ( 4, 2) =          1.99930000d0
      a ( 5, 2) =          3.06850004d0
      a ( 6, 2) =          4.69530010d0
      a ( 7, 2) =         14.47799969d0
      a ( 8, 2) =         12.93210030d0
      cc( 1, 3) =        -78.13166809d0
      cc( 2, 3) =        278.22503662d0
      cc( 3, 3) =       -349.11547852d0
      cc( 4, 3) =        318.36791992d0
      cc( 5, 3) =       -181.44857788d0
      cc( 6, 3) =        183.82826233d0
      cc( 7, 3) =         51.52063751d0
      cc( 8, 3) =          5.67256212d0
      a ( 1, 3) =          1.47309995d0
      a ( 2, 3) =          1.79939997d0
      a ( 3, 3) =          2.51900005d0
      a ( 4, 3) =          3.78139997d0
      a ( 5, 3) =          5.86429977d0
      a ( 6, 3) =          9.54720020d0
      a ( 7, 3) =         20.55660057d0
      a ( 8, 3) =         45.63219833d0
      cc( 1, 4) =        -53.33285141d0
      cc( 2, 4) =        183.93737793d0
      cc( 3, 4) =       -259.36651611d0
      cc( 4, 4) =        260.49545288d0
      cc( 5, 4) =       -166.65827942d0
      cc( 6, 4) =        140.39550781d0
      cc( 7, 4) =         66.17495728d0
      cc( 8, 4) =          7.23944521d0
      a ( 1, 4) =          0.91380000d0
      a ( 2, 4) =          1.06050003d0
      a ( 3, 4) =          1.39680004d0
      a ( 4, 4) =          2.00209999d0
      a ( 5, 4) =          3.00699997d0
      a ( 6, 4) =          4.65070009d0
      a ( 7, 4) =         13.49240017d0
      a ( 8, 4) =         41.66809845d0
      cc( 1, 5) =        -11.07681656d0
      cc( 2, 5) =         46.05147552d0
      cc( 3, 5) =        -97.47564697d0
      cc( 4, 5) =        180.38792419d0
      cc( 5, 5) =       -262.32742310d0
      cc( 6, 5) =        297.45623779d0
      cc( 7, 5) =         27.28302383d0
      cc( 8, 5) =          9.93384647d0
      a ( 1, 5) =          1.07920003d0
      a ( 2, 5) =          1.40100002d0
      a ( 3, 5) =          1.89979994d0
      a ( 4, 5) =          2.86960006d0
      a ( 5, 5) =          4.46210003d0
      a ( 6, 5) =          7.10710001d0
      a ( 7, 5) =         21.98730087d0
      a ( 8, 5) =         19.46669960d0
c
      go to 30
c
c     ra (Ra_CRENBL_ECP)
c
  88  continue
c
      cc( 1, 1) =         -1.14673400d0
      cc( 2, 1) =        -12.07713222d0
      cc( 3, 1) =        -35.34965134d0
      cc( 4, 1) =       -135.56118774d0
      cc( 5, 1) =       -295.28918457d0
      cc( 6, 1) =        -58.24514008d0
      a ( 1, 1) =          1.02110004d0
      a ( 2, 1) =          2.33030009d0
      a ( 3, 1) =          5.71700001d0
      a ( 4, 1) =         17.22739983d0
      a ( 5, 1) =         53.49580002d0
      a ( 6, 1) =        161.11909485d0
      cc( 1, 2) =        -25.06594467d0
      cc( 2, 2) =         81.25407410d0
      cc( 3, 2) =       -164.93606567d0
      cc( 4, 2) =        336.00671387d0
      cc( 5, 2) =       -412.40322876d0
      cc( 6, 2) =        435.14797974d0
      cc( 7, 2) =         39.21112061d0
      cc( 8, 2) =         15.50083256d0
      a ( 1, 2) =          0.78680003d0
      a ( 2, 2) =          0.94749999d0
      a ( 3, 2) =          1.31900001d0
      a ( 4, 2) =          2.06419992d0
      a ( 5, 2) =          3.47370005d0
      a ( 6, 2) =          6.25019979d0
      a ( 7, 2) =         23.50049973d0
      a ( 8, 2) =         19.61120033d0
      cc( 1, 3) =        -73.69863892d0
      cc( 2, 3) =        280.92184448d0
      cc( 3, 3) =       -357.36761475d0
      cc( 4, 3) =        335.35665894d0
      cc( 5, 3) =       -183.73077393d0
      cc( 6, 3) =        195.09538269d0
      cc( 7, 3) =         45.13964462d0
      cc( 8, 3) =          7.23854303d0
      a ( 1, 3) =          1.58039999d0
      a ( 2, 3) =          1.94050002d0
      a ( 3, 3) =          2.70449996d0
      a ( 4, 3) =          4.08209991d0
      a ( 5, 3) =          6.39729977d0
      a ( 6, 3) =         10.23750019d0
      a ( 7, 3) =         18.74060059d0
      a ( 8, 3) =         33.47719955d0
      cc( 1, 4) =        -53.52696228d0
      cc( 2, 4) =        190.84431458d0
      cc( 3, 4) =       -261.43304443d0
      cc( 4, 4) =        258.68441772d0
      cc( 5, 4) =       -159.79736328d0
      cc( 6, 4) =        136.84912109d0
      cc( 7, 4) =         56.20154953d0
      cc( 8, 4) =          7.48619318d0
      a ( 1, 4) =          1.00759995d0
      a ( 2, 4) =          1.17990005d0
      a ( 3, 4) =          1.56289995d0
      a ( 4, 4) =          2.27600002d0
      a ( 5, 4) =          3.43030000d0
      a ( 6, 4) =          5.27850008d0
      a ( 7, 4) =         12.86940002d0
      a ( 8, 4) =         30.29739952d0
      cc( 1, 5) =         -1.44440401d0
      cc( 2, 5) =         34.48142242d0
      cc( 3, 5) =        -96.37539673d0
      cc( 4, 5) =        186.22872925d0
      cc( 5, 5) =       -263.07485962d0
      cc( 6, 5) =        323.78308105d0
      cc( 7, 5) =         27.88471031d0
      cc( 8, 5) =          9.85942745d0
      a ( 1, 5) =          0.85030001d0
      a ( 2, 5) =          1.72609997d0
      a ( 3, 5) =          2.19359994d0
      a ( 4, 5) =          3.27309990d0
      a ( 5, 5) =          5.12659979d0
      a ( 6, 5) =          8.45429993d0
      a ( 7, 5) =         26.79529953d0
      a ( 8, 5) =         22.44179916d0
c
      go to 30
c
c     ac (Ac_CRENBL_ECP)
c
  89  continue
c
      cc( 1, 1) =         -1.20669997d0
      cc( 2, 1) =        -12.14964867d0
      cc( 3, 1) =        -36.28235626d0
      cc( 4, 1) =       -139.70959473d0
      cc( 5, 1) =       -308.30279541d0
      cc( 6, 1) =        -58.80256653d0
      a ( 1, 1) =          1.11780000d0
      a ( 2, 1) =          2.49180007d0
      a ( 3, 1) =          5.99900007d0
      a ( 4, 1) =         18.06710052d0
      a ( 5, 1) =         56.65689850d0
      a ( 6, 1) =        171.80650330d0
      cc( 1, 2) =        -22.83414078d0
      cc( 2, 2) =         85.36679077d0
      cc( 3, 2) =       -184.32901001d0
      cc( 4, 2) =        361.85287476d0
      cc( 5, 2) =       -391.65957642d0
      cc( 6, 2) =        311.78582764d0
      cc( 7, 2) =         31.91579056d0
      cc( 8, 2) =         11.15847778d0
      a ( 1, 2) =          0.97649997d0
      a ( 2, 2) =          1.16690004d0
      a ( 3, 2) =          1.54069996d0
      a ( 4, 2) =          2.26390004d0
      a ( 5, 2) =          3.54769993d0
      a ( 6, 2) =          5.63149977d0
      a ( 7, 2) =         16.93210030d0
      a ( 8, 2) =         15.99460030d0
      cc( 1, 3) =         41.84874725d0
      cc( 2, 3) =       -150.01193237d0
      cc( 3, 3) =        335.78051758d0
      cc( 4, 3) =       -384.57653809d0
      cc( 5, 3) =        335.67944336d0
      cc( 6, 3) =       -153.28929138d0
      cc( 7, 3) =         59.26359940d0
      cc( 8, 3) =          5.42515993d0
      a ( 1, 3) =          1.28699994d0
      a ( 2, 3) =          1.48230004d0
      a ( 3, 3) =          1.94239998d0
      a ( 4, 3) =          2.79509997d0
      a ( 5, 3) =          4.23729992d0
      a ( 6, 3) =          6.32700014d0
      a ( 7, 3) =          9.92879963d0
      a ( 8, 3) =         36.35929871d0
      cc( 1, 4) =         43.26000214d0
      cc( 2, 4) =       -138.38842773d0
      cc( 3, 4) =        263.83312988d0
      cc( 4, 4) =       -327.39425659d0
      cc( 5, 4) =        305.68560791d0
      cc( 6, 4) =       -145.50318909d0
      cc( 7, 4) =         48.99058914d0
      cc( 8, 4) =          7.83579922d0
      a ( 1, 4) =          0.75120002d0
      a ( 2, 4) =          0.86790001d0
      a ( 3, 4) =          1.12769997d0
      a ( 4, 4) =          1.58729994d0
      a ( 5, 4) =          2.32139993d0
      a ( 6, 4) =          3.33249998d0
      a ( 7, 4) =          4.73500013d0
      a ( 8, 4) =         15.27309990d0
      cc( 1, 5) =         -2.29345989d0
      cc( 2, 5) =         34.95658493d0
      cc( 3, 5) =       -103.22103882d0
      cc( 4, 5) =        202.08523560d0
      cc( 5, 5) =       -281.09320068d0
      cc( 6, 5) =        342.86886597d0
      cc( 7, 5) =         27.73361969d0
      cc( 8, 5) =          9.87891960d0
      a ( 1, 5) =          0.99159998d0
      a ( 2, 5) =          1.79879999d0
      a ( 3, 5) =          2.35899997d0
      a ( 4, 5) =          3.49049997d0
      a ( 5, 5) =          5.48509979d0
      a ( 6, 5) =          8.97290039d0
      a ( 7, 5) =         27.43889999d0
      a ( 8, 5) =         23.89990044d0
c
      go to 30
c
c     th (Th_CRENBL_ECP)
c
  90  continue
c
      cc( 1, 1) =         -1.18738997d0
      cc( 2, 1) =        -12.30656910d0
      cc( 3, 1) =        -36.65875244d0
      cc( 4, 1) =       -139.81327820d0
      cc( 5, 1) =       -305.43923950d0
      cc( 6, 1) =        -58.70189667d0
      a ( 1, 1) =          1.20299995d0
      a ( 2, 1) =          2.65759993d0
      a ( 3, 1) =          6.24800014d0
      a ( 4, 1) =         18.62809944d0
      a ( 5, 1) =         57.40290070d0
      a ( 6, 1) =        172.54620361d0
      cc( 1, 2) =        -24.04182243d0
      cc( 2, 2) =         78.38962555d0
      cc( 3, 2) =       -161.65811157d0
      cc( 4, 2) =        328.75106812d0
      cc( 5, 2) =       -374.78945923d0
      cc( 6, 2) =        391.01895142d0
      cc( 7, 2) =         35.13041687d0
      cc( 8, 2) =         15.77670383d0
      a ( 1, 2) =          0.90810001d0
      a ( 2, 2) =          1.07920003d0
      a ( 3, 2) =          1.47860003d0
      a ( 4, 2) =          2.26749992d0
      a ( 5, 2) =          3.71530008d0
      a ( 6, 2) =          6.41309977d0
      a ( 7, 2) =         21.16460037d0
      a ( 8, 2) =         17.85589981d0
      cc( 1, 3) =         45.82766724d0
      cc( 2, 3) =       -167.06277466d0
      cc( 3, 3) =        368.67584229d0
      cc( 4, 3) =       -433.55191040d0
      cc( 5, 3) =        400.28353882d0
      cc( 6, 3) =       -169.12930298d0
      cc( 7, 3) =         85.87816620d0
      cc( 8, 3) =          8.52397919d0
      a ( 1, 3) =          1.26510000d0
      a ( 2, 3) =          1.49430001d0
      a ( 3, 3) =          2.00500011d0
      a ( 4, 3) =          2.97460008d0
      a ( 5, 3) =          4.72440004d0
      a ( 6, 3) =          7.69339991d0
      a ( 7, 3) =         12.33460045d0
      a ( 8, 3) =         51.79029846d0
      cc( 1, 4) =         29.37446976d0
      cc( 2, 4) =       -120.37606049d0
      cc( 3, 4) =        259.13232422d0
      cc( 4, 4) =       -321.52166748d0
      cc( 5, 4) =        304.66244507d0
      cc( 6, 4) =       -143.89175415d0
      cc( 7, 4) =         51.10453033d0
      cc( 8, 4) =          7.74334002d0
      a ( 1, 4) =          0.78079998d0
      a ( 2, 4) =          0.93570000d0
      a ( 3, 4) =          1.21609998d0
      a ( 4, 4) =          1.72160006d0
      a ( 5, 4) =          2.54830003d0
      a ( 6, 4) =          3.69059992d0
      a ( 7, 4) =          5.41249990d0
      a ( 8, 4) =         17.09110069d0
      cc( 1, 5) =        -21.65902138d0
      cc( 2, 5) =         65.11306000d0
      cc( 3, 5) =       -122.52014923d0
      cc( 4, 5) =        218.33807373d0
      cc( 5, 5) =       -291.18447876d0
      cc( 6, 5) =        367.08728027d0
      cc( 7, 5) =         28.42512512d0
      cc( 8, 5) =          9.84213829d0
      a ( 1, 5) =          1.41310000d0
      a ( 2, 5) =          1.74030006d0
      a ( 3, 5) =          2.43449998d0
      a ( 4, 5) =          3.69120002d0
      a ( 5, 5) =          5.90170002d0
      a ( 6, 5) =          9.87790012d0
      a ( 7, 5) =         30.56329918d0
      a ( 8, 5) =         26.41230011d0
c
      go to 30
c
c     pa (Pa_CRENBL_ECP)
c
  91  continue
c
      cc( 1, 1) =         -1.13092303d0
      cc( 2, 1) =        -12.63854885d0
      cc( 3, 1) =        -37.00111771d0
      cc( 4, 1) =       -139.85595703d0
      cc( 5, 1) =       -306.96585083d0
      cc( 6, 1) =        -59.04650879d0
      a ( 1, 1) =          1.27950001d0
      a ( 2, 1) =          2.83470011d0
      a ( 3, 1) =          6.52720022d0
      a ( 4, 1) =         19.20420074d0
      a ( 5, 1) =         58.46409988d0
      a ( 6, 1) =        176.57690430d0
      cc( 1, 2) =        -24.29222107d0
      cc( 2, 2) =         80.66983032d0
      cc( 3, 2) =       -165.70108032d0
      cc( 4, 2) =        329.73165894d0
      cc( 5, 2) =       -360.00137329d0
      cc( 6, 2) =        372.31121826d0
      cc( 7, 2) =         33.68949509d0
      cc( 8, 2) =         15.89577580d0
      a ( 1, 2) =          0.96829998d0
      a ( 2, 2) =          1.14800000d0
      a ( 3, 2) =          1.55799997d0
      a ( 4, 2) =          2.36789989d0
      a ( 5, 2) =          3.83089995d0
      a ( 6, 2) =          6.47319984d0
      a ( 7, 2) =         20.25769997d0
      a ( 8, 2) =         17.12579918d0
      cc( 1, 3) =         57.45444107d0
      cc( 2, 3) =       -192.87774658d0
      cc( 3, 3) =        420.92343140d0
      cc( 4, 3) =       -519.91271973d0
      cc( 5, 3) =        509.81158447d0
      cc( 6, 3) =       -216.13882446d0
      cc( 7, 3) =         96.05680847d0
      cc( 8, 3) =          9.86869431d0
      a ( 1, 3) =          1.38759995d0
      a ( 2, 3) =          1.61860001d0
      a ( 3, 3) =          2.17269993d0
      a ( 4, 3) =          3.20580006d0
      a ( 5, 3) =          5.03020000d0
      a ( 6, 3) =          7.85659981d0
      a ( 7, 3) =         13.11429977d0
      a ( 8, 3) =         54.46849823d0
      cc( 1, 4) =         29.77106667d0
      cc( 2, 4) =       -151.98701477d0
      cc( 3, 4) =        298.67770386d0
      cc( 4, 4) =       -310.10202026d0
      cc( 5, 4) =        295.44259644d0
      cc( 6, 4) =       -152.87329102d0
      cc( 7, 4) =         52.50298691d0
      cc( 8, 4) =          7.72490215d0
      a ( 1, 4) =          0.85420001d0
      a ( 2, 4) =          1.03989995d0
      a ( 3, 4) =          1.29579997d0
      a ( 4, 4) =          1.82229996d0
      a ( 5, 4) =          2.77440000d0
      a ( 6, 4) =          3.91989994d0
      a ( 7, 4) =          5.65369987d0
      a ( 8, 4) =         17.96859932d0
      cc( 1, 5) =        -32.19529724d0
      cc( 2, 5) =         88.89807892d0
      cc( 3, 5) =       -149.84609985d0
      cc( 4, 5) =        249.73576355d0
      cc( 5, 5) =       -320.78536987d0
      cc( 6, 5) =        394.20562744d0
      cc( 7, 5) =         28.40286064d0
      cc( 8, 5) =          9.84354782d0
      a ( 1, 5) =          1.62249994d0
      a ( 2, 5) =          1.94340003d0
      a ( 3, 5) =          2.66190004d0
      a ( 4, 5) =          3.99289989d0
      a ( 5, 5) =          6.34350014d0
      a ( 6, 5) =         10.56330013d0
      a ( 7, 5) =         32.51839828d0
      a ( 8, 5) =         28.10650063d0
c
c     6/8/8/8/8
c
  30  continue
c
      n ( 1, 1) =   2
      n ( 2, 1) =   2
      n ( 3, 1) =   2
      n ( 4, 1) =   2
      n ( 5, 1) =   2
      n ( 6, 1) =   1
      nt( 1)    =   6
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      n ( 4, 2) =   2
      n ( 5, 2) =   2
      n ( 6, 2) =   2
      n ( 7, 2) =   1
      n ( 8, 2) =   0
      nt( 2)    =   8
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      n ( 3, 3) =   2
      n ( 4, 3) =   2
      n ( 5, 3) =   2
      n ( 6, 3) =   2
      n ( 7, 3) =   1
      n ( 8, 3) =   0
      nt( 3)    =   8
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      n ( 3, 4) =   2
      n ( 4, 4) =   2
      n ( 5, 4) =   2
      n ( 6, 4) =   2
      n ( 7, 4) =   1
      n ( 8, 4) =   0
      nt( 4)    =   8
      n ( 1, 5) =   2
      n ( 2, 5) =   2
      n ( 3, 5) =   2
      n ( 4, 5) =   2
      n ( 5, 5) =   2
      n ( 6, 5) =   2
      n ( 7, 5) =   1
      n ( 8, 5) =   0
      nt( 5)    =   8
      go to 200
c
c     u (U_CRENBL_ECP)
c
  92  continue
c
      cc( 1, 1) =         -0.95164698d0
      cc( 2, 1) =        -10.77463818d0
      cc( 3, 1) =        -33.54887009d0
      cc( 4, 1) =       -122.39160919d0
      cc( 5, 1) =       -256.04879761d0
      cc( 6, 1) =       -721.33471680d0
      cc( 7, 1) =        -75.18030548d0
      a ( 1, 1) =          1.22290003d0
      a ( 2, 1) =          2.67100000d0
      a ( 3, 1) =          6.10900021d0
      a ( 4, 1) =         17.91930008d0
      a ( 5, 1) =         49.88119888d0
      a ( 6, 1) =        169.55189514d0
      a ( 7, 1) =        605.90167236d0
      cc( 1, 2) =         86.94699860d0
      cc( 2, 2) =       -324.48245239d0
      cc( 3, 2) =        754.80963135d0
      cc( 4, 2) =       -931.61145020d0
      cc( 5, 2) =        867.38934326d0
      cc( 6, 2) =       -567.78674316d0
      cc( 7, 2) =        467.56951904d0
      cc( 8, 2) =         87.23509216d0
      cc( 9, 2) =          6.00922394d0
      a ( 1, 2) =          2.08200002d0
      a ( 2, 2) =          2.36159992d0
      a ( 3, 2) =          3.04959989d0
      a ( 4, 2) =          4.28889990d0
      a ( 5, 2) =          6.36810017d0
      a ( 6, 2) =          9.73639965d0
      a ( 7, 2) =         15.35929966d0
      a ( 8, 2) =         43.85210037d0
      a ( 9, 2) =        131.51210022d0
      cc( 1, 3) =        109.52915955d0
      cc( 2, 3) =       -372.34240723d0
      cc( 3, 3) =        760.67779541d0
      cc( 4, 3) =       -945.12261963d0
      cc( 5, 3) =        906.44946289d0
      cc( 6, 3) =       -618.62194824d0
      cc( 7, 3) =        434.29833984d0
      cc( 8, 3) =         96.91146851d0
      cc( 9, 3) =          8.63370800d0
      a ( 1, 3) =          1.55610001d0
      a ( 2, 3) =          1.77209997d0
      a ( 3, 3) =          2.25839996d0
      a ( 4, 3) =          3.09450006d0
      a ( 5, 3) =          4.51030016d0
      a ( 6, 3) =          6.71129990d0
      a ( 7, 3) =         10.01770020d0
      a ( 8, 3) =         28.29829979d0
      a ( 9, 3) =         89.10600281d0
      cc( 1, 4) =         98.69554901d0
      cc( 2, 4) =       -348.77438354d0
      cc( 3, 4) =        647.68463135d0
      cc( 4, 4) =       -796.75195312d0
      cc( 5, 4) =        775.70019531d0
      cc( 6, 4) =       -576.68280029d0
      cc( 7, 4) =        328.07855225d0
      cc( 8, 4) =         70.24051666d0
      cc( 9, 4) =          7.22914982d0
      a ( 1, 4) =          1.03709996d0
      a ( 2, 4) =          1.17840004d0
      a ( 3, 4) =          1.44579995d0
      a ( 4, 4) =          1.91680002d0
      a ( 5, 4) =          2.65330005d0
      a ( 6, 4) =          3.74600005d0
      a ( 7, 4) =          5.03719997d0
      a ( 8, 4) =         15.92879963d0
      a ( 9, 4) =         48.05899811d0
      cc( 1, 5) =         -0.77639502d0
      cc( 2, 5) =         64.14081573d0
      cc( 3, 5) =       -301.19494629d0
      cc( 4, 5) =        617.67492676d0
      cc( 5, 5) =       -840.64794922d0
      cc( 6, 5) =       1010.84838867d0
      cc( 7, 5) =       -503.98593140d0
      cc( 8, 5) =        116.52574158d0
      cc( 9, 5) =          8.05050755d0
      a ( 1, 5) =          0.90359998d0
      a ( 2, 5) =          2.84689999d0
      a ( 3, 5) =          3.67079997d0
      a ( 4, 5) =          4.89400005d0
      a ( 5, 5) =          7.39130020d0
      a ( 6, 5) =         11.65960026d0
      a ( 7, 5) =         18.25749969d0
      a ( 8, 5) =         28.92289925d0
      a ( 9, 5) =        137.56460571d0
c
      go to 40
c
c     Np (Np_CRENBL_ECP)
c
  93  continue
c
      cc( 1, 1) =         -0.59611499d0
      cc( 2, 1) =         -8.56580639d0
      cc( 3, 1) =        -30.77449989d0
      cc( 4, 1) =       -101.92440033d0
      cc( 5, 1) =       -227.16531372d0
      cc( 6, 1) =       -631.36932373d0
      cc( 7, 1) =        -71.25617981d0
      a ( 1, 1) =          1.16480005d0
      a ( 2, 1) =          2.50749993d0
      a ( 3, 1) =          5.64090013d0
      a ( 4, 1) =         16.30699921d0
      a ( 5, 1) =         42.56740189d0
      a ( 6, 1) =        143.51150513d0
      a ( 7, 1) =        490.29238892d0
      cc( 1, 2) =        163.53736877d0
      cc( 2, 2) =       -407.81634521d0
      cc( 3, 2) =        733.90496826d0
      cc( 4, 2) =       -892.30194092d0
      cc( 5, 2) =        787.23022461d0
      cc( 6, 2) =       -517.18475342d0
      cc( 7, 2) =        524.23852539d0
      cc( 8, 2) =         98.42464447d0
      cc( 9, 2) =          5.98039103d0
      a ( 1, 2) =          2.21970010d0
      a ( 2, 2) =          2.41639996d0
      a ( 3, 2) =          3.12549996d0
      a ( 4, 2) =          4.42549992d0
      a ( 5, 2) =          6.40640020d0
      a ( 6, 2) =         10.66429996d0
      a ( 7, 2) =         16.00900078d0
      a ( 8, 2) =         49.05659866d0
      a ( 9, 2) =        168.13330078d0
      cc( 1, 3) =        115.25938416d0
      cc( 2, 3) =       -379.42489624d0
      cc( 3, 3) =        740.25213623d0
      cc( 4, 3) =       -907.79357910d0
      cc( 5, 3) =        882.03515625d0
      cc( 6, 3) =       -598.22576904d0
      cc( 7, 3) =        426.77703857d0
      cc( 8, 3) =         97.61193085d0
      cc( 9, 3) =          8.64615154d0
      a ( 1, 3) =          1.61839998d0
      a ( 2, 3) =          1.83500004d0
      a ( 3, 3) =          2.32730007d0
      a ( 4, 3) =          3.19989991d0
      a ( 5, 3) =          4.62799978d0
      a ( 6, 3) =          6.85830021d0
      a ( 7, 3) =         10.20699978d0
      a ( 8, 3) =         28.45260048d0
      a ( 9, 3) =         89.37139893d0
      cc( 1, 4) =        126.47127533d0
      cc( 2, 4) =       -397.60015869d0
      cc( 3, 4) =        716.23443604d0
      cc( 4, 4) =       -891.56164551d0
      cc( 5, 4) =        872.76300049d0
      cc( 6, 4) =       -641.78558350d0
      cc( 7, 4) =        368.73754883d0
      cc( 8, 4) =         76.54608917d0
      cc( 9, 4) =          7.11474991d0
      a ( 1, 4) =          1.11280000d0
      a ( 2, 4) =          1.25590003d0
      a ( 3, 4) =          1.56729996d0
      a ( 4, 4) =          2.11560011d0
      a ( 5, 4) =          3.01139998d0
      a ( 6, 4) =          4.34480000d0
      a ( 7, 4) =          6.03599977d0
      a ( 8, 4) =         19.41950035d0
      a ( 9, 4) =         62.30830002d0
      cc( 1, 5) =         -0.98151898d0
      cc( 2, 5) =         41.70447922d0
      cc( 3, 5) =       -261.12640381d0
      cc( 4, 5) =        560.24981689d0
      cc( 5, 5) =       -803.19482422d0
      cc( 6, 5) =       1050.93518066d0
      cc( 7, 5) =       -540.55218506d0
      cc( 8, 5) =        114.44293976d0
      cc( 9, 5) =          8.04148579d0
      a ( 1, 5) =          0.98909998d0
      a ( 2, 5) =          2.71499991d0
      a ( 3, 5) =          3.74250007d0
      a ( 4, 5) =          4.93669987d0
      a ( 5, 5) =          7.59490013d0
      a ( 6, 5) =         11.83139992d0
      a ( 7, 5) =         17.86820030d0
      a ( 8, 5) =         28.38409996d0
      a ( 9, 5) =        130.11839294d0
c
      go to 40
c
c     Pu (Pu_CRENBL_ECP)
c
  94  continue
c
      cc( 1, 1) =         -0.88483602d0
      cc( 2, 1) =        -10.61913300d0
      cc( 3, 1) =        -33.61550522d0
      cc( 4, 1) =       -120.75714111d0
      cc( 5, 1) =       -256.49505615d0
      cc( 6, 1) =       -725.57281494d0
      cc( 7, 1) =        -75.82495117d0
      a ( 1, 1) =          1.31490004d0
      a ( 2, 1) =          2.85489988d0
      a ( 3, 1) =          6.45160007d0
      a ( 4, 1) =         18.74729919d0
      a ( 5, 1) =         51.15800095d0
      a ( 6, 1) =        173.26420593d0
      a ( 7, 1) =        617.64038086d0
      cc( 1, 2) =         88.92472839d0
      cc( 2, 2) =       -330.37750244d0
      cc( 3, 2) =        750.95227051d0
      cc( 4, 2) =       -908.51391602d0
      cc( 5, 2) =        873.33062744d0
      cc( 6, 2) =       -567.66186523d0
      cc( 7, 2) =        466.55908203d0
      cc( 8, 2) =         83.55030060d0
      cc( 9, 2) =          6.17316198d0
      a ( 1, 2) =          2.21939993d0
      a ( 2, 2) =          2.51889992d0
      a ( 3, 2) =          3.23670006d0
      a ( 4, 2) =          4.51700020d0
      a ( 5, 2) =          6.62480021d0
      a ( 6, 2) =         10.00730038d0
      a ( 7, 2) =         15.32999992d0
      a ( 8, 2) =         39.74290085d0
      a ( 9, 2) =        107.35160065d0
      cc( 1, 3) =        151.63360596d0
      cc( 2, 3) =       -307.76840210d0
      cc( 3, 3) =        556.88562012d0
      cc( 4, 3) =       -860.99200439d0
      cc( 5, 3) =        879.28405762d0
      cc( 6, 3) =       -487.93322754d0
      cc( 7, 3) =        359.29550171d0
      cc( 8, 3) =         99.86724854d0
      cc( 9, 3) =          8.64875412d0
      a ( 1, 3) =          1.64100003d0
      a ( 2, 3) =          1.78100002d0
      a ( 3, 3) =          2.42759991d0
      a ( 4, 3) =          3.38750005d0
      a ( 5, 3) =          4.67859983d0
      a ( 6, 3) =          6.86380005d0
      a ( 7, 3) =         10.95189953d0
      a ( 8, 3) =         29.43799973d0
      a ( 9, 3) =         93.42259979d0
      cc( 1, 4) =        118.47415161d0
      cc( 2, 4) =       -396.36349487d0
      cc( 3, 4) =        706.21350098d0
      cc( 4, 4) =       -875.85577393d0
      cc( 5, 4) =        863.94415283d0
      cc( 6, 4) =       -574.99090576d0
      cc( 7, 4) =        318.59259033d0
      cc( 8, 4) =         79.30906677d0
      cc( 9, 4) =          7.11830807d0
      a ( 1, 4) =          1.15670002d0
      a ( 2, 4) =          1.31159997d0
      a ( 3, 4) =          1.61979997d0
      a ( 4, 4) =          2.19190001d0
      a ( 5, 4) =          3.07089996d0
      a ( 6, 4) =          4.36579990d0
      a ( 7, 4) =          6.33850002d0
      a ( 8, 4) =         20.19569969d0
      a ( 9, 4) =         67.26270294d0
      cc( 1, 5) =         -1.48127198d0
      cc( 2, 5) =        162.76635742d0
      cc( 3, 5) =       -264.43045044d0
      cc( 4, 5) =        369.88125610d0
      cc( 5, 5) =       -704.01580811d0
      cc( 6, 5) =        971.67071533d0
      cc( 7, 5) =       -443.89147949d0
      cc( 8, 5) =        108.53334045d0
      cc( 9, 5) =          8.09724331d0
      a ( 1, 5) =          1.11049998d0
      a ( 2, 5) =          2.80390000d0
      a ( 3, 5) =          3.11529994d0
      a ( 4, 5) =          5.01070023d0
      a ( 5, 5) =          7.77769995d0
      a ( 6, 5) =         11.78960037d0
      a ( 7, 5) =         17.65460014d0
      a ( 8, 5) =         29.48049927d0
      a ( 9, 5) =        120.29889679d0
c
  40  continue
      n ( 1, 1) =   2
      n ( 2, 1) =   2
      n ( 3, 1) =   2
      n ( 4, 1) =   2
      n ( 5, 1) =   2
      n ( 6, 1) =   2
      n ( 7, 1) =   1
      nt( 1)    =   7
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      n ( 4, 2) =   2
      n ( 5, 2) =   2
      n ( 6, 2) =   2
      n ( 7, 2) =   2
      n ( 8, 2) =   1
      n ( 9, 2) =   0
      nt( 2)    =   9
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      n ( 3, 3) =   2
      n ( 4, 3) =   2
      n ( 5, 3) =   2
      n ( 6, 3) =   2
      n ( 7, 3) =   2
      n ( 8, 3) =   1
      n ( 9, 3) =   0
      nt( 3)    =   9
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      n ( 3, 4) =   2
      n ( 4, 4) =   2
      n ( 5, 4) =   2
      n ( 6, 4) =   2
      n ( 7, 4) =   2
      n ( 8, 4) =   1
      n ( 9, 4) =   0
      nt( 4)    =   9
      n ( 1, 5) =   2
      n ( 2, 5) =   2
      n ( 3, 5) =   2
      n ( 4, 5) =   2
      n ( 5, 5) =   2
      n ( 6, 5) =   2
      n ( 7, 5) =   2
      n ( 8, 5) =   1
      n ( 9, 5) =   0
      nt( 5)    =   9
      go to 200
c
c     Am (AM_CRENBL_ECP)
c
  95  continue
c
      cc( 1, 1) =         -1.17866504d0
      cc( 2, 1) =        -13.27001476d0
      cc( 3, 1) =        -39.10317993d0
      cc( 4, 1) =       -155.02056885d0
      cc( 5, 1) =       -354.33300781d0
      cc( 6, 1) =        -61.31810379d0
      a ( 1, 1) =          1.44700003d0
      a ( 2, 1) =          3.21869993d0
      a ( 3, 1) =          7.59910011d0
      a ( 4, 1) =         22.65139961d0
      a ( 5, 1) =         71.29299927d0
      a ( 6, 1) =        220.40229797d0
      cc( 1, 2) =        -28.39175415d0
      cc( 2, 2) =         84.97696686d0
      cc( 3, 2) =       -171.02401733d0
      cc( 4, 2) =        357.06427002d0
      cc( 5, 2) =       -387.68786621d0
      cc( 6, 2) =        441.84158325d0
      cc( 7, 2) =         35.83626938d0
      cc( 8, 2) =         15.82378292d0
      a ( 1, 2) =          1.09730005d0
      a ( 2, 2) =          1.29470003d0
      a ( 3, 2) =          1.80130005d0
      a ( 4, 2) =          2.78760004d0
      a ( 5, 2) =          4.60900021d0
      a ( 6, 2) =          8.16489983d0
      a ( 7, 2) =         26.83819962d0
      a ( 8, 2) =         21.21400070d0
      cc( 1, 3) =         64.26541138d0
      cc( 2, 3) =       -210.60380554d0
      cc( 3, 3) =        442.71279907d0
      cc( 4, 3) =       -492.18869019d0
      cc( 5, 3) =        437.84417725d0
      cc( 6, 3) =       -168.01190186d0
      cc( 7, 3) =         70.56053925d0
      cc( 8, 3) =          6.54860497d0
      a ( 1, 3) =          1.70260000d0
      a ( 2, 3) =          1.94879997d0
      a ( 3, 3) =          2.54539990d0
      a ( 4, 3) =          3.65510011d0
      a ( 5, 3) =          5.55019999d0
      a ( 6, 3) =          8.03429985d0
      a ( 7, 3) =         14.14560032d0
      a ( 8, 3) =         48.91350174d0
      cc( 1, 4) =         28.79652405d0
      cc( 2, 4) =       -105.51664734d0
      cc( 3, 4) =        221.98576355d0
      cc( 4, 4) =       -264.26501465d0
      cc( 5, 4) =        261.34280396d0
      cc( 6, 4) =       -110.51098633d0
      cc( 7, 4) =         57.14907837d0
      cc( 8, 4) =          7.61405420d0
      a ( 1, 4) =          0.92580003d0
      a ( 2, 4) =          1.11629999d0
      a ( 3, 4) =          1.49710000d0
      a ( 4, 4) =          2.19849992d0
      a ( 5, 4) =          3.37220001d0
      a ( 6, 4) =          5.19850016d0
      a ( 7, 4) =          7.96210003d0
      a ( 8, 4) =         24.13820076d0
      cc( 1, 5) =         -2.47670388d0
      cc( 2, 5) =         37.72669983d0
      cc( 3, 5) =       -114.31232452d0
      cc( 4, 5) =        236.72572327d0
      cc( 5, 5) =       -310.35437012d0
      cc( 6, 5) =        466.90872192d0
      cc( 7, 5) =         30.38970375d0
      cc( 8, 5) =          9.76973534d0
      a ( 1, 5) =          1.24930000d0
      a ( 2, 5) =          2.41829991d0
      a ( 3, 5) =          3.24990010d0
      a ( 4, 5) =          4.97840023d0
      a ( 5, 5) =          8.13679981d0
      a ( 6, 5) =         14.51039982d0
      a ( 7, 5) =         44.86399841d0
      a ( 8, 5) =         38.13169861d0
c
      go to 30
c
c     Cm (Cm_CRENBL_ECP)
c
  96  continue
c
      cc( 1, 1) =         -1.27089500d0
      cc( 2, 1) =        -13.83742809d0
      cc( 3, 1) =        -40.88765335d0
      cc( 4, 1) =       -160.03480530d0
      cc( 5, 1) =       -364.42514038d0
      cc( 6, 1) =        -61.78550720d0
      a ( 1, 1) =          1.54879999d0
      a ( 2, 1) =          3.41910005d0
      a ( 3, 1) =          8.07050037d0
      a ( 4, 1) =         23.85969925d0
      a ( 5, 1) =         75.01390076d0
      a ( 6, 1) =        231.42030334d0
      cc( 1, 2) =        -26.33760071d0
      cc( 2, 2) =         85.68217468d0
      cc( 3, 2) =       -174.58818054d0
      cc( 4, 2) =        355.67098999d0
      cc( 5, 2) =       -370.72149658d0
      cc( 6, 2) =        423.92636108d0
      cc( 7, 2) =         35.04350281d0
      cc( 8, 2) =         15.92432594d0
      a ( 1, 2) =          1.15559995d0
      a ( 2, 2) =          1.37430000d0
      a ( 3, 2) =          1.88199997d0
      a ( 4, 2) =          2.89520001d0
      a ( 5, 2) =          4.74020004d0
      a ( 6, 2) =          8.27089977d0
      a ( 7, 2) =         25.20420074d0
      a ( 8, 2) =         21.01460075d0
      cc( 1, 3) =        -35.03964233d0
      cc( 2, 3) =        113.44233704d0
      cc( 3, 3) =       -223.71525574d0
      cc( 4, 3) =        424.26696777d0
      cc( 5, 3) =       -484.90432739d0
      cc( 6, 3) =        375.64773560d0
      cc( 7, 3) =         33.73357773d0
      cc( 8, 3) =         11.89339828d0
      a ( 1, 3) =          1.12440002d0
      a ( 2, 3) =          1.29729998d0
      a ( 3, 3) =          1.68669999d0
      a ( 4, 3) =          2.41960001d0
      a ( 5, 3) =          3.68249989d0
      a ( 6, 3) =          5.54580021d0
      a ( 7, 3) =         17.95700073d0
      a ( 8, 3) =         16.02359962d0
      cc( 1, 4) =         39.77001572d0
      cc( 2, 4) =       -126.65122223d0
      cc( 3, 4) =        255.05882263d0
      cc( 4, 4) =       -305.24252319d0
      cc( 5, 4) =        300.05279541d0
      cc( 6, 4) =       -129.08760071d0
      cc( 7, 4) =         60.35084152d0
      cc( 8, 4) =          7.49130201d0
      a ( 1, 4) =          1.00460005d0
      a ( 2, 4) =          1.18710005d0
      a ( 3, 4) =          1.60969996d0
      a ( 4, 4) =          2.36179996d0
      a ( 5, 4) =          3.64210010d0
      a ( 6, 4) =          5.57789993d0
      a ( 7, 4) =          8.67679977d0
      a ( 8, 4) =         27.02840042d0
      cc( 1, 5) =         -2.70642400d0
      cc( 2, 5) =         37.51979065d0
      cc( 3, 5) =       -111.37210846d0
      cc( 4, 5) =        228.42926025d0
      cc( 5, 5) =       -297.12841797d0
      cc( 6, 5) =        458.37960815d0
      cc( 7, 5) =         30.03723907d0
      cc( 8, 5) =          9.80161095d0
      a ( 1, 5) =          1.30649996d0
      a ( 2, 5) =          2.45950007d0
      a ( 3, 5) =          3.30599999d0
      a ( 4, 5) =          5.06689978d0
      a ( 5, 5) =          8.23499966d0
      a ( 6, 5) =         14.63000011d0
      a ( 7, 5) =         44.62329864d0
      a ( 8, 5) =         37.61240005d0
c
      go to 30
c
c     Bk (Bk_CRENBL_ECP)
c
  97  continue
c
      cc( 1, 1) =         -1.91033494d0
      cc( 2, 1) =        -18.14081764d0
      cc( 3, 1) =        -49.25227356d0
      cc( 4, 1) =       -186.21107483d0
      cc( 5, 1) =       -442.21273804d0
      cc( 6, 1) =        -64.78657532d0
      a ( 1, 1) =          1.77400005d0
      a ( 2, 1) =          4.02939987d0
      a ( 3, 1) =          9.97329998d0
      a ( 4, 1) =         28.40850067d0
      a ( 5, 1) =         94.63580322d0
      a ( 6, 1) =        300.05819702d0
      cc( 1, 2) =         52.19206619d0
      cc( 2, 2) =       -186.36035156d0
      cc( 3, 2) =         92.78298950d0
      cc( 4, 2) =        334.95520020d0
      cc( 5, 2) =       -465.58764648d0
      cc( 6, 2) =        553.61657715d0
      cc( 7, 2) =         29.15875244d0
      cc( 8, 2) =         16.74713707d0
      a ( 1, 2) =          1.79050004d0
      a ( 2, 2) =          2.20539999d0
      a ( 3, 2) =          3.21810007d0
      a ( 4, 2) =          3.21900010d0
      a ( 5, 2) =          5.33370018d0
      a ( 6, 2) =         10.14509964d0
      a ( 7, 2) =         23.30200005d0
      a ( 8, 2) =         31.06459999d0
      cc( 1, 3) =         51.22914124d0
      cc( 2, 3) =       -181.65086365d0
      cc( 3, 3) =        443.10906982d0
      cc( 4, 3) =       -565.48089600d0
      cc( 5, 3) =        596.35418701d0
      cc( 6, 3) =       -236.82652283d0
      cc( 7, 3) =        111.00143433d0
      cc( 8, 3) =          9.82366085d0
      a ( 1, 3) =          1.84089994d0
      a ( 2, 3) =          2.15310001d0
      a ( 3, 3) =          2.89759994d0
      a ( 4, 3) =          4.39190006d0
      a ( 5, 3) =          7.07049990d0
      a ( 6, 3) =         11.38609982d0
      a ( 7, 3) =         19.07290077d0
      a ( 8, 3) =         78.96360016d0
      cc( 1, 4) =         38.53779602d0
      cc( 2, 4) =       -128.37193298d0
      cc( 3, 4) =        265.80664062d0
      cc( 4, 4) =       -312.57681274d0
      cc( 5, 4) =        306.93249512d0
      cc( 6, 4) =       -127.93990326d0
      cc( 7, 4) =         61.20502853d0
      cc( 8, 4) =          7.49189091d0
      a ( 1, 4) =          1.07249999d0
      a ( 2, 4) =          1.27520001d0
      a ( 3, 4) =          1.72510004d0
      a ( 4, 4) =          2.52390003d0
      a ( 5, 4) =          3.89089990d0
      a ( 6, 4) =          5.89279985d0
      a ( 7, 4) =          9.33549976d0
      a ( 8, 4) =         28.41810036d0
      cc( 1, 5) =         -2.38182592d0
      cc( 2, 5) =         37.63102341d0
      cc( 3, 5) =       -112.93108368d0
      cc( 4, 5) =        235.85006714d0
      cc( 5, 5) =       -315.54455566d0
      cc( 6, 5) =        493.09857178d0
      cc( 7, 5) =         30.15993309d0
      cc( 8, 5) =          9.81385136d0
      a ( 1, 5) =          1.33210003d0
      a ( 2, 5) =          2.63820004d0
      a ( 3, 5) =          3.53020000d0
      a ( 4, 5) =          5.43610001d0
      a ( 5, 5) =          8.95960045d0
      a ( 6, 5) =         15.72469997d0
      a ( 7, 5) =         45.85969925d0
      a ( 8, 5) =         41.15259933d0
c
      go to 30
c
c     Cf (Cf_CRENBL_ECP)
c
  98  continue
c
      cc( 1, 1) =          0.06396900d0
      cc( 2, 1) =          0.72026199d0
      cc( 3, 1) =          3.60362005d0
      cc( 4, 1) =       -331.19174194d0
      cc( 5, 1) =        -49.50713348d0
      cc( 6, 1) =        -64.19592285d0
      a ( 1, 1) =          0.11240000d0
      a ( 2, 1) =          0.33410001d0
      a ( 3, 1) =          0.86870003d0
      a ( 4, 1) =         35.57229996d0
      a ( 5, 1) =          5.74980021d0
      a ( 6, 1) =        173.79170227d0
      cc( 1, 2) =         -0.02315000d0
      cc( 2, 2) =         -0.24423000d0
      cc( 3, 2) =         -2.08339190d0
      cc( 4, 2) =       -112.74409485d0
      cc( 5, 2) =        361.11566162d0
      cc( 6, 2) =       -304.97174072d0
      cc( 7, 2) =         91.31772614d0
      cc( 8, 2) =          9.39082432d0
      a ( 1, 2) =          0.08660000d0
      a ( 2, 2) =          0.21210000d0
      a ( 3, 2) =          0.56459999d0
      a ( 4, 2) =          2.45079994d0
      a ( 5, 2) =          3.32960010d0
      a ( 6, 2) =          5.10309982d0
      a ( 7, 2) =          7.63899994d0
      a ( 8, 2) =         45.90340042d0
      cc( 1, 3) =         -0.24964599d0
      cc( 2, 3) =         -4.33974504d0
      cc( 3, 3) =          0.84404802d0
      cc( 4, 3) =         -4.78782701d0
      cc( 5, 3) =         77.99574280d0
      cc( 6, 3) =        -16.03907776d0
      cc( 7, 3) =         46.87809753d0
      cc( 8, 3) =          5.98568678d0
      a ( 1, 3) =          0.17330000d0
      a ( 2, 3) =          0.72009999d0
      a ( 3, 3) =          0.94270003d0
      a ( 4, 3) =          1.73440003d0
      a ( 5, 3) =          2.79830003d0
      a ( 6, 3) =          4.64459991d0
      a ( 7, 3) =          9.47770023d0
      a ( 8, 3) =         22.03520012d0
      cc( 1, 4) =          0.00683800d0
      cc( 2, 4) =         -1.26423001d0
      cc( 3, 4) =         23.25278091d0
      cc( 4, 4) =        -61.92188644d0
      cc( 5, 4) =        115.66687775d0
      cc( 6, 4) =        -23.18271065d0
      cc( 7, 4) =         48.15451050d0
      cc( 8, 4) =          8.24353027d0
      a ( 1, 4) =          0.04830000d0
      a ( 2, 4) =          0.57410002d0
      a ( 3, 4) =          1.29489994d0
      a ( 4, 4) =          2.13150001d0
      a ( 5, 4) =          3.01090002d0
      a ( 6, 4) =          5.09649992d0
      a ( 7, 4) =          7.41559982d0
      a ( 8, 4) =         19.88430023d0
      cc( 1, 5) =         -0.02929600d0
      cc( 2, 5) =         -3.00971389d0
      cc( 3, 5) =         -0.38204300d0
      cc( 4, 5) =         -8.69149971d0
      cc( 5, 5) =         86.56561279d0
      cc( 6, 5) =       -115.70391846d0
      cc( 7, 5) =         34.37308121d0
      cc( 8, 5) =          9.23519802d0
      a ( 1, 5) =          0.09060000d0
      a ( 2, 5) =          0.68760002d0
      a ( 3, 5) =          0.24670000d0
      a ( 4, 5) =          1.96379995d0
      a ( 5, 5) =          5.77269983d0
      a ( 6, 5) =          8.67479992d0
      a ( 7, 5) =         16.13820076d0
      a ( 8, 5) =          9.27400017d0
c
      go to 30
c
c     Es (Es_CRENBL_ECP)
c
  99  continue
c
      cc( 1, 1) =          0.22005400d0
      cc( 2, 1) =          2.60089493d0
      cc( 3, 1) =        -31.24035454d0
      cc( 4, 1) =        -69.07208252d0
      cc( 5, 1) =       -153.30227661d0
      cc( 6, 1) =        -49.51485825d0
      a ( 1, 1) =          0.28020000d0
      a ( 2, 1) =          0.84520000d0
      a ( 3, 1) =          4.92689991d0
      a ( 4, 1) =         15.65719986d0
      a ( 5, 1) =         33.65110016d0
      a ( 6, 1) =         84.97429657d0
      cc( 1, 2) =        -14.49978828d0
      cc( 2, 2) =         -0.08987100d0
      cc( 3, 2) =         46.88783264d0
      cc( 4, 2) =       -134.73852539d0
      cc( 5, 2) =        332.50048828d0
      cc( 6, 2) =       -312.70547485d0
      cc( 7, 2) =        120.46249390d0
      cc( 8, 2) =         13.29056931d0
      a ( 1, 2) =          0.90420002d0
      a ( 2, 2) =          0.21770000d0
      a ( 3, 2) =          1.25109994d0
      a ( 4, 2) =          1.93330002d0
      a ( 5, 2) =          3.12360001d0
      a ( 6, 2) =          5.08769989d0
      a ( 7, 2) =          8.04069996d0
      a ( 8, 2) =         51.81460190d0
      cc( 1, 3) =         -0.08161800d0
      cc( 2, 3) =        -40.43290329d0
      cc( 3, 3) =         42.01081467d0
      cc( 4, 3) =        -96.41140747d0
      cc( 5, 3) =        169.88287354d0
      cc( 6, 3) =        -79.35518646d0
      cc( 7, 3) =         57.57525253d0
      cc( 8, 3) =          7.18364286d0
      a ( 1, 3) =          0.21280000d0
      a ( 2, 3) =          0.76670003d0
      a ( 3, 3) =          0.78860003d0
      a ( 4, 3) =          1.87339997d0
      a ( 5, 3) =          2.26539993d0
      a ( 6, 3) =          4.71110010d0
      a ( 7, 3) =          5.98540020d0
      a ( 8, 3) =         24.30929947d0
      cc( 1, 4) =          0.01627700d0
      cc( 2, 4) =         -4.06413889d0
      cc( 3, 4) =        122.57620239d0
      cc( 4, 4) =       -208.42051697d0
      cc( 5, 4) =        160.10260010d0
      cc( 6, 4) =        -78.40477753d0
      cc( 7, 4) =         59.69111252d0
      cc( 8, 4) =          7.67357492d0
      a ( 1, 4) =          0.09400000d0
      a ( 2, 4) =          0.91619998d0
      a ( 3, 4) =          1.60539997d0
      a ( 4, 4) =          1.99860001d0
      a ( 5, 4) =          2.80159998d0
      a ( 6, 4) =          5.72270012d0
      a ( 7, 4) =          5.63560009d0
      a ( 8, 4) =         21.45039940d0
      cc( 1, 5) =         -0.03291600d0
      cc( 2, 5) =         -0.60518497d0
      cc( 3, 5) =        -17.50750542d0
      cc( 4, 5) =         37.31097031d0
      cc( 5, 5) =        -64.48777771d0
      cc( 6, 5) =         64.95304871d0
      cc( 7, 5) =         34.13174820d0
      cc( 8, 5) =          9.33467960d0
      a ( 1, 5) =          0.17290001d0
      a ( 2, 5) =          0.43399999d0
      a ( 3, 5) =          1.32980001d0
      a ( 4, 5) =          1.73720002d0
      a ( 5, 5) =          2.51580000d0
      a ( 6, 5) =          3.58400011d0
      a ( 7, 5) =         17.79479980d0
      a ( 8, 5) =         15.36159992d0
c
      go to 30
c
c     Fm (Fm_CRENBL_ECP)
c
  100 continue
c
      cc( 1, 1) =          0.20972200d0
      cc( 2, 1) =          2.38716412d0
      cc( 3, 1) =        -23.88536072d0
      cc( 4, 1) =        -41.27284241d0
      cc( 5, 1) =       -195.20947266d0
      cc( 6, 1) =        -52.42905426d0
      a ( 1, 1) =          0.37490001d0
      a ( 2, 1) =          1.00329995d0
      a ( 3, 1) =          4.66419983d0
      a ( 4, 1) =         11.07190037d0
      a ( 5, 1) =         30.67930031d0
      a ( 6, 1) =         94.67870331d0
      cc( 1, 2) =        -34.44987488d0
      cc( 2, 2) =         -0.05397000d0
      cc( 3, 2) =         93.45252228d0
      cc( 4, 2) =       -166.93705750d0
      cc( 5, 2) =        315.45736694d0
      cc( 6, 2) =       -287.31268311d0
      cc( 7, 2) =        121.74865723d0
      cc( 8, 2) =         13.41475391d0
      a ( 1, 2) =          1.12090003d0
      a ( 2, 2) =          0.27300000d0
      a ( 3, 2) =          1.38479996d0
      a ( 4, 2) =          1.96840000d0
      a ( 5, 2) =          3.18589997d0
      a ( 6, 2) =          5.44890022d0
      a ( 7, 2) =          8.25049973d0
      a ( 8, 2) =         53.43569946d0
      cc( 1, 3) =         -0.07136800d0
      cc( 2, 3) =         -4.89569998d0
      cc( 3, 3) =        109.29483032d0
      cc( 4, 3) =       -264.04949951d0
      cc( 5, 3) =        416.31094360d0
      cc( 6, 3) =       -281.26287842d0
      cc( 7, 3) =         56.75691986d0
      cc( 8, 3) =          7.06949806d0
      a ( 1, 3) =          0.29049999d0
      a ( 2, 3) =          0.89389998d0
      a ( 3, 3) =          1.55190003d0
      a ( 4, 3) =          1.84819996d0
      a ( 5, 3) =          2.59710002d0
      a ( 6, 3) =          3.55520010d0
      a ( 7, 3) =          5.33809996d0
      a ( 8, 3) =         23.67040062d0
      cc( 1, 4) =          0.02905100d0
      cc( 2, 4) =        -66.32601166d0
      cc( 3, 4) =        222.19554138d0
      cc( 4, 4) =       -283.73937988d0
      cc( 5, 4) =        282.01510620d0
      cc( 6, 4) =       -125.96237183d0
      cc( 7, 4) =         57.05357742d0
      cc( 8, 4) =          7.82386208d0
      a ( 1, 4) =          0.15090001d0
      a ( 2, 4) =          1.38989997d0
      a ( 3, 4) =          1.65369999d0
      a ( 4, 4) =          2.24569988d0
      a ( 5, 4) =          3.42560005d0
      a ( 6, 4) =          5.16309977d0
      a ( 7, 4) =          7.54209995d0
      a ( 8, 4) =         24.02599907d0
      cc( 1, 5) =         -2.39030695d0
      cc( 2, 5) =         -0.09625900d0
      cc( 3, 5) =        -38.28890610d0
      cc( 4, 5) =        111.31268311d0
      cc( 5, 5) =        -75.48458099d0
      cc( 6, 5) =         88.39505768d0
      cc( 7, 5) =          8.39762402d0
      a ( 1, 5) =          0.85829997d0
      a ( 2, 5) =          0.29679999d0
      a ( 3, 5) =          3.45919991d0
      a ( 4, 5) =          5.04069996d0
      a ( 5, 5) =          7.66330004d0
      a ( 6, 5) =         17.76759911d0
      a ( 7, 5) =         81.46340179d0
c
      go to 50
c
c     Md (Md_CRENBL_ECP)
c
  101 continue
c
      cc( 1, 1) =          0.15976401d0
      cc( 2, 1) =          2.23393393d0
      cc( 3, 1) =        -20.63683128d0
      cc( 4, 1) =        -32.38039780d0
      cc( 5, 1) =       -187.27101135d0
      cc( 6, 1) =        -53.30691528d0
      a ( 1, 1) =          0.43070000d0
      a ( 2, 1) =          1.11600006d0
      a ( 3, 1) =          4.73820019d0
      a ( 4, 1) =          9.46599960d0
      a ( 5, 1) =         28.08659935d0
      a ( 6, 1) =         91.21839905d0
      cc( 1, 2) =         10.10878086d0
      cc( 2, 2) =        -36.09821701d0
      cc( 3, 2) =         74.37078094d0
      cc( 4, 2) =       -152.26487732d0
      cc( 5, 2) =        334.36880493d0
      cc( 6, 2) =       -349.59790039d0
      cc( 7, 2) =        141.49493408d0
      cc( 8, 2) =         13.19567204d0
      a ( 1, 2) =          0.75559998d0
      a ( 2, 2) =          0.90410000d0
      a ( 3, 2) =          1.26919997d0
      a ( 4, 2) =          1.96039999d0
      a ( 5, 2) =          3.35109997d0
      a ( 6, 2) =          6.06659985d0
      a ( 7, 2) =          9.84689999d0
      a ( 8, 2) =         76.16600037d0
      cc( 1, 3) =         -0.65193802d0
      cc( 2, 3) =         27.55850983d0
      cc( 3, 3) =       -140.45182800d0
      cc( 4, 3) =        354.18386841d0
      cc( 5, 3) =       -382.05096436d0
      cc( 6, 3) =        271.13302612d0
      cc( 7, 3) =         65.55524445d0
      cc( 8, 3) =          5.57984114d0
      a ( 1, 3) =          0.61199999d0
      a ( 2, 3) =          1.76080000d0
      a ( 3, 3) =          2.16370010d0
      a ( 4, 3) =          2.94819999d0
      a ( 5, 3) =          4.36439991d0
      a ( 6, 3) =          6.09040022d0
      a ( 7, 3) =         20.19359970d0
      a ( 8, 3) =         59.76610184d0
      cc( 1, 4) =          0.04919700d0
      cc( 2, 4) =        -64.92922211d0
      cc( 3, 4) =        223.80036926d0
      cc( 4, 4) =       -295.63165283d0
      cc( 5, 4) =        295.39566040d0
      cc( 6, 4) =       -127.92171478d0
      cc( 7, 4) =         57.74354553d0
      cc( 8, 4) =          7.82849216d0
      a ( 1, 4) =          0.22090000d0
      a ( 2, 4) =          1.46790004d0
      a ( 3, 4) =          1.75820005d0
      a ( 4, 4) =          2.40960002d0
      a ( 5, 4) =          3.59540010d0
      a ( 6, 4) =          5.42309999d0
      a ( 7, 4) =          7.74739981d0
      a ( 8, 4) =         24.39170074d0
      cc( 1, 5) =         13.13120461d0
      cc( 2, 5) =         -0.21304201d0
      cc( 3, 5) =         -5.23304892d0
      cc( 4, 5) =        -81.41338348d0
      cc( 5, 5) =        203.75080872d0
      cc( 6, 5) =       -217.60070801d0
      cc( 7, 5) =         98.32410431d0
      cc( 8, 5) =          8.20197201d0
      a ( 1, 5) =          2.43199992d0
      a ( 2, 5) =          0.44159999d0
      a ( 3, 5) =          1.28919995d0
      a ( 4, 5) =          4.00540018d0
      a ( 5, 5) =          5.93550014d0
      a ( 6, 5) =          9.76329994d0
      a ( 7, 5) =         15.58800030d0
      a ( 8, 5) =         89.99559784d0
c
      go to 30
c
c     No (No_CRENBL_ECP)
c
  102 continue
c
      cc( 1, 1) =         -1.75959897d0
      cc( 2, 1) =        -19.19492912d0
      cc( 3, 1) =        -50.64709091d0
      cc( 4, 1) =       -200.09509277d0
      cc( 5, 1) =       -486.48617554d0
      cc( 6, 1) =        -67.62830353d0
      a ( 1, 1) =          2.02789998d0
      a ( 2, 1) =          4.57170010d0
      a ( 3, 1) =         11.58650017d0
      a ( 4, 1) =         33.05160141d0
      a ( 5, 1) =        110.28289795d0
      a ( 6, 1) =        356.88180542d0
      cc( 1, 2) =         -5.24792719d0
      cc( 2, 2) =         53.06208420d0
      cc( 3, 2) =       -117.28907776d0
      cc( 4, 2) =        306.65054321d0
      cc( 5, 2) =       -331.41055298d0
      cc( 6, 2) =        128.60786438d0
      cc( 7, 2) =         12.21696186d0
      a ( 1, 2) =          0.94800001d0
      a ( 2, 2) =          1.45229995d0
      a ( 3, 2) =          1.93630004d0
      a ( 4, 2) =          3.56380010d0
      a ( 5, 2) =          6.03779984d0
      a ( 6, 2) =          9.50010014d0
      a ( 7, 2) =         49.94449997d0
      cc( 1, 3) =        -38.58179855d0
      cc( 2, 3) =        127.71244049d0
      cc( 3, 3) =       -256.02554321d0
      cc( 4, 3) =        525.35491943d0
      cc( 5, 3) =       -633.99108887d0
      cc( 6, 3) =        508.19781494d0
      cc( 7, 3) =         36.28185272d0
      cc( 8, 3) =         11.82044506d0
      a ( 1, 3) =          1.48730004d0
      a ( 2, 3) =          1.71399999d0
      a ( 3, 3) =          2.21630001d0
      a ( 4, 3) =          3.18919992d0
      a ( 5, 3) =          4.90679979d0
      a ( 6, 3) =          7.49809980d0
      a ( 7, 3) =         24.33909988d0
      a ( 8, 3) =         21.45919991d0
      cc( 1, 4) =         42.61248398d0
      cc( 2, 4) =       -141.19110107d0
      cc( 3, 4) =        286.28744507d0
      cc( 4, 4) =       -338.34857178d0
      cc( 5, 4) =        343.32589722d0
      cc( 6, 4) =       -140.28236389d0
      cc( 7, 4) =         66.82376099d0
      cc( 8, 4) =          7.40976000d0
      a ( 1, 4) =          1.21459997d0
      a ( 2, 4) =          1.45009995d0
      a ( 3, 4) =          1.97140002d0
      a ( 4, 4) =          2.94039989d0
      a ( 5, 4) =          4.58269978d0
      a ( 6, 4) =          7.18160009d0
      a ( 7, 4) =         11.16399956d0
      a ( 8, 4) =         34.78409958d0
      cc( 1, 5) =          0.17846601d0
      cc( 2, 5) =         -0.33072501d0
      cc( 3, 5) =        -67.09947968d0
      cc( 4, 5) =        216.27134705d0
      cc( 5, 5) =       -237.39591980d0
      cc( 6, 5) =        104.37260437d0
      cc( 7, 5) =          8.08529472d0
      a ( 1, 5) =          1.18110001d0
      a ( 2, 5) =          0.78350002d0
      a ( 3, 5) =          5.19530010d0
      a ( 4, 5) =          6.89720011d0
      a ( 5, 5) =         11.00150013d0
      a ( 6, 5) =         17.62080002d0
      a ( 7, 5) =        102.36190033d0
c
      go to 70

c     Lw (Lw_CRENBL_ECP)
c
  103 continue
c
      cc( 1, 1) =         -1.44170296d0
      cc( 2, 1) =        -15.03919792d0
      cc( 3, 1) =        -39.65955353d0
      cc( 4, 1) =       -161.56266785d0
      cc( 5, 1) =       -360.18804932d0
      cc( 6, 1) =        -62.81684494d0
      a ( 1, 1) =          2.02850008d0
      a ( 2, 1) =          4.39300013d0
      a ( 3, 1) =          9.83810043d0
      a ( 4, 1) =         28.35020065d0
      a ( 5, 1) =         82.61139679d0
      a ( 6, 1) =        250.89920044d0
      cc( 1, 2) =         49.13726425d0
      cc( 2, 2) =       -191.04298401d0
      cc( 3, 2) =        209.44651794d0
      cc( 4, 2) =        260.47091675d0
      cc( 5, 2) =       -488.95681763d0
      cc( 6, 2) =        658.06842041d0
      cc( 7, 2) =         43.20983887d0
      cc( 8, 2) =         15.95556164d0
      a ( 1, 2) =          2.16720009d0
      a ( 2, 2) =          2.72460008d0
      a ( 3, 2) =          4.01079988d0
      a ( 4, 2) =          4.01510000d0
      a ( 5, 2) =          6.69960022d0
      a ( 6, 2) =         13.21790028d0
      a ( 7, 2) =         36.08919907d0
      a ( 8, 2) =         40.45830154d0
      cc( 1, 3) =         -7.35526609d0
      cc( 2, 3) =         28.43832397d0
      cc( 3, 3) =        -66.70023346d0
      cc( 4, 3) =        186.59317017d0
      cc( 5, 3) =       -142.15013123d0
      cc( 6, 3) =        190.02313232d0
      cc( 7, 3) =         43.40534973d0
      cc( 8, 3) =          7.59309292d0
      a ( 1, 3) =          1.34739995d0
      a ( 2, 3) =          1.55369997d0
      a ( 3, 3) =          2.02690005d0
      a ( 4, 3) =          3.07389998d0
      a ( 5, 3) =          5.02370024d0
      a ( 6, 3) =          8.87870026d0
      a ( 7, 3) =         20.10149956d0
      a ( 8, 3) =         27.58690071d0
      cc( 1, 4) =         35.04140091d0
      cc( 2, 4) =       -113.56825256d0
      cc( 3, 4) =        235.13388062d0
      cc( 4, 4) =       -271.48348999d0
      cc( 5, 4) =        282.18148804d0
      cc( 6, 4) =       -105.98326874d0
      cc( 7, 4) =         59.50698853d0
      cc( 8, 4) =          7.74769497d0
      a ( 1, 4) =          1.19560003d0
      a ( 2, 4) =          1.42390001d0
      a ( 3, 4) =          1.93990004d0
      a ( 4, 4) =          2.83879995d0
      a ( 5, 4) =          4.37580013d0
      a ( 6, 4) =          6.63329983d0
      a ( 7, 4) =         10.27820015d0
      a ( 8, 4) =         28.45750046d0
      cc( 1, 5) =         -9.85134220d0
      cc( 2, 5) =         60.02737427d0
      cc( 3, 5) =       -125.03877258d0
      cc( 4, 5) =        224.62619019d0
      cc( 5, 5) =       -251.34989929d0
      cc( 6, 5) =        106.51850891d0
      cc( 7, 5) =          8.20704937d0
      a ( 1, 5) =          2.03099990d0
      a ( 2, 5) =          2.98690009d0
      a ( 3, 5) =          4.02019978d0
      a ( 4, 5) =          6.56209993d0
      a ( 5, 5) =         11.02750015d0
      a ( 6, 5) =         17.35169983d0
      a ( 7, 5) =        109.06250000d0
c
c     6/8/8/8/7
 50   continue
      n ( 1, 1) =   2
      n ( 2, 1) =   2
      n ( 3, 1) =   2
      n ( 4, 1) =   2
      n ( 5, 1) =   2
      n ( 6, 1) =   1
      nt( 1)    =   6
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      n ( 4, 2) =   2
      n ( 5, 2) =   2
      n ( 6, 2) =   2
      n ( 7, 2) =   1
      n ( 8, 2) =   0
      nt( 2)    =   8
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      n ( 3, 3) =   2
      n ( 4, 3) =   2
      n ( 5, 3) =   2
      n ( 6, 3) =   2
      n ( 7, 3) =   1
      n ( 8, 3) =   0
      nt( 3)    =   8
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      n ( 3, 4) =   2
      n ( 4, 4) =   2
      n ( 5, 4) =   2
      n ( 6, 4) =   2
      n ( 7, 4) =   1
      n ( 8, 4) =   0
      nt( 4)    =   8
      n ( 1, 5) =   2
      n ( 2, 5) =   2
      n ( 3, 5) =   2
      n ( 4, 5) =   2
      n ( 5, 5) =   2
      n ( 6, 5) =   1
      n ( 7, 5) =   0
      nt( 5)    =   7
c
      go to 200
c
c     6/7/8/8/7
c
 70   continue
      n ( 1, 1) =   2
      n ( 2, 1) =   2
      n ( 3, 1) =   2
      n ( 4, 1) =   2
      n ( 5, 1) =   2
      n ( 6, 1) =   1
      nt( 1)    =   6
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      n ( 4, 2) =   2
      n ( 5, 2) =   2
      n ( 6, 2) =   1
      n ( 7, 2) =   0
      nt( 2)    =   7
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      n ( 3, 3) =   2
      n ( 4, 3) =   2
      n ( 5, 3) =   2
      n ( 6, 3) =   2
      n ( 7, 3) =   1
      n ( 8, 3) =   0
      nt( 3)    =   8
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      n ( 3, 4) =   2
      n ( 4, 4) =   2
      n ( 5, 4) =   2
      n ( 6, 4) =   2
      n ( 7, 4) =   1
      n ( 8, 4) =   0
      nt( 4)    =   8
      n ( 1, 5) =   2
      n ( 2, 5) =   2
      n ( 3, 5) =   2
      n ( 4, 5) =   2
      n ( 5, 5) =   2
      n ( 6, 5) =   1
      n ( 7, 5) =   0
      nt( 5)    =   7
c
 200  nelcor    =  78
      lmx       =   4
c
      return
c
 6010 format (1x,'*** data not found for CRENBL ecp type ',a4)
      end
**==crenbs0.f
      subroutine crenbs0 (zinp,lmx,nelcor,nt,n,cc,a)
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
      implicit character * 8 (z)
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
      dimension nt(6), n(10,6), cc(10,6), a(10,6)
c
      dimension ztit(10)
      data maxtyp/10/
      data ztit/'sc','ti','v ','cr','mn','fe',
     *          'co','ni','cu','zn'/
c
c     this routine contains a library of ecp data taken from
c     various literature sources.
c     1st transition series CRENBS potentials
c     
      do 10 ityp = 1 , maxtyp
         if (zinp.eq.ztit(ityp)) go to 20
 10   continue
      write (iwr,6010) zinp
      call caserr2('no library available for requested ecp')
 20   n ( 1, 1) =   2
      n ( 2, 1) =   2
      n ( 3, 1) =   2
      n ( 4, 1) =   2
      n ( 5, 1) =   1
      nt( 1)    =   5
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      n ( 4, 2) =   2
      n ( 5, 2) =   1
      n ( 6, 2) =   0
      nt( 2)    =   6
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      n ( 3, 3) =   2
      n ( 4, 3) =   2
      n ( 5, 3) =   1
      n ( 6, 3) =   0
      nt( 3)    =   6
      nelcor    =  18
      lmx       =   2
c
      go to (21,22,23,24,25,26,27,28,29,30), ityp
c
c     scandium (Sc_CRENBS_ECP)
c
  21  continue
      cc( 1, 1) =         -0.27631900d0
      cc( 2, 1) =         -2.71528506d0
      cc( 3, 1) =        -10.29641724d0
      cc( 4, 1) =        -36.81774139d0
      cc( 5, 1) =        -13.80766487d0
      a ( 1, 1) =          0.29840001d0
      a ( 2, 1) =          0.97359997d0
      a ( 3, 1) =          3.00230002d0
      a ( 4, 1) =         12.39669991d0
      a ( 5, 1) =         47.46879959d0
      cc( 1, 2) =          0.03074700d0
      cc( 2, 2) =         98.84182739d0
      cc( 3, 2) =       -176.84428406d0
      cc( 4, 2) =         99.04251862d0
      cc( 5, 2) =          8.55434704d0
      cc( 6, 2) =          3.60357404d0
      a ( 1, 2) =          0.18490000d0
      a ( 2, 2) =          1.01020002d0
      a ( 3, 2) =          1.17219996d0
      a ( 4, 2) =          1.40620005d0
      a ( 5, 2) =          4.10559988d0
      a ( 6, 2) =          5.19129992d0
      cc( 1, 3) =         52.92183685d0
      cc( 2, 3) =       -114.06956482d0
      cc( 3, 3) =        111.66765594d0
      cc( 4, 3) =        -52.40676498d0
      cc( 5, 3) =          5.54833794d0
      cc( 6, 3) =          5.54230881d0
      a ( 1, 3) =          0.74400002d0
      a ( 2, 3) =          0.93519998d0
      a ( 3, 3) =          1.25870001d0
      a ( 4, 3) =          1.63880003d0
      a ( 5, 3) =          2.00399995d0
      a ( 6, 3) =          1.07050002d0
c
      return
c
c     titanium (Ti_CRENBS_ECP)
c
  22  continue
      cc( 1, 1) =         -0.30006000d0
      cc( 2, 1) =         -2.86175990d0
      cc( 3, 1) =        -10.93051720d0
      cc( 4, 1) =        -36.18664169d0
      cc( 5, 1) =        -12.74508858d0
      a ( 1, 1) =          0.33790001d0
      a ( 2, 1) =          1.11440003d0
      a ( 3, 1) =          3.42820001d0
      a ( 4, 1) =         13.60970020d0
      a ( 5, 1) =         45.69940186d0
      cc( 1, 2) =          0.03480600d0
      cc( 2, 2) =         21.71248817d0
      cc( 3, 2) =        -39.98028564d0
      cc( 4, 2) =          1.93583095d0
      cc( 5, 2) =         24.44411659d0
      cc( 6, 2) =          3.18623209d0
      a ( 1, 2) =          0.21480000d0
      a ( 2, 2) =          0.99339998d0
      a ( 3, 2) =          1.17780006d0
      a ( 4, 2) =          1.81480002d0
      a ( 5, 2) =          0.93210000d0
      a ( 6, 2) =         10.55080032d0
      cc( 1, 3) =         21.59829712d0
      cc( 2, 3) =        -72.59468079d0
      cc( 3, 3) =         79.00105286d0
      cc( 4, 3) =        -38.15525436d0
      cc( 5, 3) =          7.75176477d0
      cc( 6, 3) =          5.45599794d0
      a ( 1, 3) =          0.76319998d0
      a ( 2, 3) =          1.06750000d0
      a ( 3, 3) =          1.40939999d0
      a ( 4, 3) =          1.81729996d0
      a ( 5, 3) =          0.79560000d0
      a ( 6, 3) =          0.58350003d0
c
      return
c
c     vanadium (V_CRENBS_ECP)
c
  23  continue
      cc( 1, 1) =         -0.31635800d0
      cc( 2, 1) =         -2.86627412d0
      cc( 3, 1) =        -10.88373375d0
      cc( 4, 1) =        -30.68353462d0
      cc( 5, 1) =        -11.54613400d0
      a ( 1, 1) =          0.37729999d0
      a ( 2, 1) =          1.23520005d0
      a ( 3, 1) =          3.70639992d0
      a ( 4, 1) =         13.44499969d0
      a ( 5, 1) =         37.79880142d0
      cc( 1, 2) =          0.04627600d0
      cc( 2, 2) =        116.58562469d0
      cc( 3, 2) =       -206.26571655d0
      cc( 4, 2) =        116.92503357d0
      cc( 5, 2) =         17.72877121d0
      cc( 6, 2) =          2.97051811d0
      a ( 1, 2) =          0.24120000d0
      a ( 2, 2) =          1.31369996d0
      a ( 3, 2) =          1.55100000d0
      a ( 4, 2) =          1.89559996d0
      a ( 5, 2) =          7.14799976d0
      a ( 6, 2) =         20.98839951d0
      cc( 1, 3) =          0.02383400d0
      cc( 2, 3) =         85.66498566d0
      cc( 3, 3) =       -168.82945251d0
      cc( 4, 3) =        101.25710297d0
      cc( 5, 3) =          6.58652687d0
      cc( 6, 3) =          5.59004593d0
      a ( 1, 3) =          0.20580000d0
      a ( 2, 3) =          0.94440001d0
      a ( 3, 3) =          1.10290003d0
      a ( 4, 3) =          1.31379998d0
      a ( 5, 3) =          3.75839996d0
      a ( 6, 3) =          3.91380000d0
c
      return
c
c     chromium (Cr_CRENBS_ECP)
c
  24  continue
      cc( 1, 1) =         -0.18984801d0
      cc( 2, 1) =         -2.21707392d0
      cc( 3, 1) =        -10.17161274d0
      cc( 4, 1) =        -31.33251953d0
      cc( 5, 1) =        -12.17855644d0
      a ( 1, 1) =          0.33739999d0
      a ( 2, 1) =          1.10440004d0
      a ( 3, 1) =          3.44440007d0
      a ( 4, 1) =         13.09500027d0
      a ( 5, 1) =         41.70100021d0
      cc( 1, 2) =         20.50030136d0
      cc( 2, 2) =        -65.95387268d0
      cc( 3, 2) =        115.60340881d0
      cc( 4, 2) =        -68.62667847d0
      cc( 5, 2) =         13.32377911d0
      cc( 6, 2) =          3.15280890d0
      a ( 1, 2) =          0.60399997d0
      a ( 2, 2) =          0.70929998d0
      a ( 3, 2) =          0.90420002d0
      a ( 4, 2) =          1.09179997d0
      a ( 5, 2) =          1.73950005d0
      a ( 6, 2) =          8.37220001d0
      cc( 1, 3) =          0.01963400d0
      cc( 2, 3) =         96.54696655d0
      cc( 3, 3) =       -186.33720398d0
      cc( 4, 3) =        109.43738556d0
      cc( 5, 3) =          7.51034021d0
      cc( 6, 3) =          5.51086378d0
      a ( 1, 3) =          0.17560001d0
      a ( 2, 3) =          1.02059996d0
      a ( 3, 3) =          1.19340003d0
      a ( 4, 3) =          1.44360006d0
      a ( 5, 3) =          4.16370010d0
      a ( 6, 3) =          4.56899977d0
c
      return
c
c     manganese (Mn_CRENBS_ECP)
c
  25  continue
      cc( 1, 1) =         -0.32237199d0
      cc( 2, 1) =         -2.42131090d0
      cc( 3, 1) =         -9.07335377d0
      cc( 4, 1) =        -23.51095963d0
      cc( 5, 1) =        -11.55800056d0
      a ( 1, 1) =          0.45190001d0
      a ( 2, 1) =          1.38629997d0
      a ( 3, 1) =          3.74930000d0
      a ( 4, 1) =         11.57509995d0
      a ( 5, 1) =         33.64369965d0
      cc( 1, 2) =         50.56281281d0
      cc( 2, 2) =       -139.92572021d0
      cc( 3, 2) =        187.17520142d0
      cc( 4, 2) =        -96.70987701d0
      cc( 5, 2) =         18.29615593d0
      cc( 6, 2) =          3.13673401d0
      a ( 1, 2) =          0.74610001d0
      a ( 2, 2) =          0.83649999d0
      a ( 3, 2) =          1.00870001d0
      a ( 4, 2) =          1.19529998d0
      a ( 5, 2) =          1.95889997d0
      a ( 6, 2) =          8.92039967d0
      cc( 1, 3) =         35.46660995d0
      cc( 2, 3) =       -103.00804138d0
      cc( 3, 3) =        106.40811157d0
      cc( 4, 3) =        -43.28920364d0
      cc( 5, 3) =         10.65728283d0
      cc( 6, 3) =          5.33598995d0
      a ( 1, 3) =          1.09370005d0
      a ( 2, 3) =          1.48570001d0
      a ( 3, 3) =          1.97080004d0
      a ( 4, 3) =          2.59910011d0
      a ( 5, 3) =          1.04449999d0
      a ( 6, 3) =          3.00909996d0
c
      return
c
c     iron (Fe_CRENBS_ECP)
c
  26  continue
      cc( 1, 1) =         -0.33538201d0
      cc( 2, 1) =         -2.39878988d0
      cc( 3, 1) =         -8.88721943d0
      cc( 4, 1) =        -23.24546051d0
      cc( 5, 1) =        -11.52505875d0
      a ( 1, 1) =          0.49520001d0
      a ( 2, 1) =          1.50390005d0
      a ( 3, 1) =          3.96690011d0
      a ( 4, 1) =         11.90670013d0
      a ( 5, 1) =         34.95520020d0
      cc( 1, 2) =         53.32651138d0
      cc( 2, 2) =       -144.99348450d0
      cc( 3, 2) =        193.69319153d0
      cc( 4, 2) =       -101.95478058d0
      cc( 5, 2) =         20.21462250d0
      cc( 6, 2) =          3.22412705d0
      a ( 1, 2) =          0.81730002d0
      a ( 2, 2) =          0.91530001d0
      a ( 3, 2) =          1.10710001d0
      a ( 4, 2) =          1.31669998d0
      a ( 5, 2) =          2.08579993d0
      a ( 6, 2) =          8.91240025d0
      cc( 1, 3) =          0.02352700d0
      cc( 2, 3) =         92.19576263d0
      cc( 3, 3) =       -187.78385925d0
      cc( 4, 3) =        121.00855255d0
      cc( 5, 3) =         12.36842060d0
      cc( 6, 3) =          5.20341396d0
      a ( 1, 3) =          0.22290000d0
      a ( 2, 3) =          1.29540002d0
      a ( 3, 3) =          1.55760002d0
      a ( 4, 3) =          1.89040005d0
      a ( 5, 3) =          4.91750002d0
      a ( 6, 3) =          7.78889990d0
c
      return
c
c     cobalt (Co_CRENBS_ECP)
c
  27  continue
      cc( 1, 1) =         -0.32827699d0
      cc( 2, 1) =         -2.07979298d0
      cc( 3, 1) =         -7.91124582d0
      cc( 4, 1) =        -20.62563133d0
      cc( 5, 1) =        -11.21210480d0
      a ( 1, 1) =          0.53230000d0
      a ( 2, 1) =          1.53649998d0
      a ( 3, 1) =          3.87420011d0
      a ( 4, 1) =         10.96930027d0
      a ( 5, 1) =         32.24440002d0
      cc( 1, 2) =         41.98735046d0
      cc( 2, 2) =       -118.17536163d0
      cc( 3, 2) =        165.50724792d0
      cc( 4, 2) =        -86.30142975d0
      cc( 5, 2) =          6.80684519d0
      cc( 6, 2) =          3.70646691d0
      a ( 1, 2) =          0.86699998d0
      a ( 2, 2) =          0.98699999d0
      a ( 3, 2) =          1.22259998d0
      a ( 4, 2) =          1.48570001d0
      a ( 5, 2) =          2.06960011d0
      a ( 6, 2) =          1.47860003d0
      cc( 1, 3) =          0.02688800d0
      cc( 2, 3) =         97.55522156d0
      cc( 3, 3) =       -178.49310303d0
      cc( 4, 3) =        108.75386047d0
      cc( 5, 3) =         12.92386913d0
      cc( 6, 3) =          5.22732306d0
      a ( 1, 3) =          0.24690001d0
      a ( 2, 3) =          1.40540004d0
      a ( 3, 3) =          1.66750002d0
      a ( 4, 3) =          2.07929993d0
      a ( 5, 3) =          5.26149988d0
      a ( 6, 3) =          8.09350014d0
c
      return
c
c     nickel (Ni_CRENBS_ECP)
c
  28  continue
      cc( 1, 1) =         -0.39820799d0
      cc( 2, 1) =         -3.11549497d0
      cc( 3, 1) =        -11.88599491d0
      cc( 4, 1) =        -32.50999451d0
      cc( 5, 1) =        -11.93502808d0
      a ( 1, 1) =          0.60470003d0
      a ( 2, 1) =          1.93700004d0
      a ( 3, 1) =          5.49459982d0
      a ( 4, 1) =         18.08950043d0
      a ( 5, 1) =         53.18939972d0
      cc( 1, 2) =         43.79937744d0
      cc( 2, 2) =       -122.48828125d0
      cc( 3, 2) =        170.38290405d0
      cc( 4, 2) =        -88.12361908d0
      cc( 5, 2) =          8.00157547d0
      cc( 6, 2) =          3.67468309d0
      a ( 1, 2) =          0.93790001d0
      a ( 2, 2) =          1.06649995d0
      a ( 3, 2) =          1.32009995d0
      a ( 4, 2) =          1.60450006d0
      a ( 5, 2) =          2.26979995d0
      a ( 6, 2) =          1.70940006d0
      cc( 1, 3) =          0.02744700d0
      cc( 2, 3) =         94.02155304d0
      cc( 3, 3) =       -187.11549377d0
      cc( 4, 3) =        122.75196838d0
      cc( 5, 3) =         12.05337429d0
      cc( 6, 3) =          5.36371279d0
      a ( 1, 3) =          0.26530001d0
      a ( 2, 3) =          1.50919998d0
      a ( 3, 3) =          1.81280005d0
      a ( 4, 3) =          2.21199989d0
      a ( 5, 3) =          5.48559999d0
      a ( 6, 3) =          7.51889992d0
c
      return
c
c     copper (Cu_CRENBS_ECP)
c
  29  continue
      cc( 1, 1) =         -0.17287099d0
      cc( 2, 1) =         -1.67033803d0
      cc( 3, 1) =         -8.22487736d0
      cc( 4, 1) =        -23.39514351d0
      cc( 5, 1) =        -11.42614460d0
      a ( 1, 1) =          0.47510001d0
      a ( 2, 1) =          1.42369997d0
      a ( 3, 1) =          4.03959990d0
      a ( 4, 1) =         12.44390011d0
      a ( 5, 1) =         38.74039841d0
      cc( 1, 2) =          6.22903204d0
      cc( 2, 2) =        -48.31729889d0
      cc( 3, 2) =        113.06871033d0
      cc( 4, 2) =        -65.06185150d0
      cc( 5, 2) =         19.98191833d0
      cc( 6, 2) =          3.21701002d0
      a ( 1, 2) =          0.81150001d0
      a ( 2, 2) =          1.10739994d0
      a ( 3, 2) =          1.38030005d0
      a ( 4, 2) =          1.65380001d0
      a ( 5, 2) =          2.90939999d0
      a ( 6, 2) =         10.47480011d0
      cc( 1, 3) =         35.69552994d0
      cc( 2, 3) =       -103.46006775d0
      cc( 3, 3) =        151.97230530d0
      cc( 4, 3) =        -87.70287323d0
      cc( 5, 3) =         20.80924225d0
      cc( 6, 3) =          4.79295111d0
      a ( 1, 3) =          0.73470002d0
      a ( 2, 3) =          0.83590001d0
      a ( 3, 3) =          1.03230000d0
      a ( 4, 3) =          1.25300002d0
      a ( 5, 3) =          2.32179999d0
      a ( 6, 3) =          9.20419979d0
c
      return
c
c     zinc (Zn_CRENBS_ECP)
c
  30  continue
      cc( 1, 1) =         -0.42563301d0
      cc( 2, 1) =         -3.32303691d0
      cc( 3, 1) =        -12.90146923d0
      cc( 4, 1) =        -35.57994080d0
      cc( 5, 1) =        -12.05187607d0
      a ( 1, 1) =          0.70620000d0
      a ( 2, 1) =          2.26999998d0
      a ( 3, 1) =          6.47790003d0
      a ( 4, 1) =         21.37000084d0
      a ( 5, 1) =         63.25559998d0
      cc( 1, 2) =         54.52636337d0
      cc( 2, 2) =       -151.99613953d0
      cc( 3, 2) =        214.51309204d0
      cc( 4, 2) =       -113.67844391d0
      cc( 5, 2) =         20.79129219d0
      cc( 6, 2) =          3.02661109d0
      a ( 1, 2) =          1.13189995d0
      a ( 2, 2) =          1.28470004d0
      a ( 3, 2) =          1.58990002d0
      a ( 4, 2) =          1.93540001d0
      a ( 5, 2) =          3.16650009d0
      a ( 6, 2) =         15.25699997d0
      cc( 1, 3) =         37.36040115d0
      cc( 2, 3) =       -108.31131744d0
      cc( 3, 3) =        166.68620300d0
      cc( 4, 3) =       -102.95578766d0
      cc( 5, 3) =          9.70077801d0
      cc( 6, 3) =          5.59525824d0
      a ( 1, 3) =          0.84560001d0
      a ( 2, 3) =          0.96240002d0
      a ( 3, 3) =          1.19850004d0
      a ( 4, 3) =          1.46300006d0
      a ( 5, 3) =          2.19120002d0
      a ( 6, 3) =          1.61720002d0
c
      return
 6010 format (1x,'*** data not found for CRENBS ecp type ',a4)
      end
**==crenbs1.f
      subroutine crenbs1 (zinp,lmx,nelcor,nt,n,cc,a)
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
      implicit character * 8 (z)
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
      dimension nt(6), n(10,6), cc(10,6), a(10,6)
c
      dimension ztit(10)
      data maxtyp/10/
      data ztit/'y ','zr','nb','mo','tc','ru',
     *          'rh','pd','ag','cd' /
c
c     this routine contains a library of ecp data taken from
c     various literature sources.
c     2nd transition series CRENBS potentials
c     
      do 10 ityp = 1 , maxtyp
         if (zinp.eq.ztit(ityp)) go to 20
 10   continue
      write (iwr,6010) zinp
      call caserr2('no library available for requested ecp')
c
 20   nelcor = 36
      lmx = 3
c
      go to (39,40,41,42,43,44,45,46,47,48), ityp
c
c      yttrium (Y_CRENBS_ECP)
c
  39  continue
      cc( 1, 1) =         -0.38814399d0
      cc( 2, 1) =         -3.85477591d0
      cc( 3, 1) =         -9.36550999d0
      cc( 4, 1) =        -39.12420654d0
      cc( 5, 1) =        -93.54985046d0
      cc( 6, 1) =        -25.21870422d0
      a ( 1, 1) =          0.40310001d0
      a ( 2, 1) =          0.95459998d0
      a ( 3, 1) =          2.47849989d0
      a ( 4, 1) =          7.70930004d0
      a ( 5, 1) =         23.46579933d0
      a ( 6, 1) =         73.73259735d0
      cc( 1, 2) =        -27.56872368d0
      cc( 2, 2) =         80.27277374d0
      cc( 3, 2) =       -101.27806091d0
      cc( 4, 2) =         67.88915253d0
      cc( 5, 2) =         37.90885544d0
      cc( 6, 2) =         34.85948563d0
      cc( 7, 2) =          3.96560407d0
      a ( 1, 2) =          0.52429998d0
      a ( 2, 2) =          0.61140001d0
      a ( 3, 2) =          0.80330002d0
      a ( 4, 2) =          0.99680001d0
      a ( 5, 2) =          3.88000011d0
      a ( 6, 2) =          8.83150005d0
      a ( 7, 2) =         22.29000092d0
      cc( 1, 3) =        127.96752930d0
      cc( 2, 3) =       -265.12304688d0
      cc( 3, 3) =        258.52642822d0
      cc( 4, 3) =       -181.08963013d0
      cc( 5, 3) =         94.94129181d0
      cc( 6, 3) =         16.67901993d0
      cc( 7, 3) =          5.70845318d0
      a ( 1, 3) =          0.68919998d0
      a ( 2, 3) =          0.81430000d0
      a ( 3, 3) =          1.08150005d0
      a ( 4, 3) =          1.51940000d0
      a ( 5, 3) =          2.09270000d0
      a ( 6, 3) =          6.04059982d0
      a ( 7, 3) =          7.35090017d0
      cc( 1, 4) =         -0.18540600d0
      cc( 2, 4) =         22.32200432d0
      cc( 3, 4) =        -70.94950867d0
      cc( 4, 4) =        131.66903687d0
      cc( 5, 4) =       -108.04507446d0
      cc( 6, 4) =         11.55035400d0
      cc( 7, 4) =          7.16109085d0
      a ( 1, 4) =          0.23090000d0
      a ( 2, 4) =          1.29820001d0
      a ( 3, 4) =          1.51409996d0
      a ( 4, 4) =          1.99890006d0
      a ( 5, 4) =          2.56730008d0
      a ( 6, 4) =          4.38630009d0
      a ( 7, 4) =          1.97099996d0
      n ( 1, 1) =   2
      n ( 2, 1) =   2
      n ( 3, 1) =   2
      n ( 4, 1) =   2
      n ( 5, 1) =   2
      n ( 6, 1) =   1
      nt( 1)    =   6
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      n ( 4, 2) =   2
      n ( 5, 2) =   2
      n ( 6, 2) =   1
      n ( 7, 2) =   0
      nt( 2)    =   7
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      n ( 3, 3) =   2
      n ( 4, 3) =   2
      n ( 5, 3) =   2
      n ( 6, 3) =   1
      n ( 7, 3) =   0
      nt( 3)    =   7
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      n ( 3, 4) =   2
      n ( 4, 4) =   2
      n ( 5, 4) =   2
      n ( 6, 4) =   1
      n ( 7, 4) =   0
      nt( 4)    =   7
c
      return
c
c     zirconium (Zr_CRENBS_ECP)
c
  40  continue
      cc( 1, 1) =         -0.65333903d0
      cc( 2, 1) =         -5.42695379d0
      cc( 3, 1) =        -15.36555386d0
      cc( 4, 1) =        -64.25210571d0
      cc( 5, 1) =        -23.12731552d0
      a ( 1, 1) =          0.50849998d0
      a ( 2, 1) =          1.25989997d0
      a ( 3, 1) =          3.85929990d0
      a ( 4, 1) =         12.27550030d0
      a ( 5, 1) =         38.51760101d0
      cc( 1, 2) =        -60.99807739d0
      cc( 2, 2) =        189.74023438d0
      cc( 3, 2) =       -245.02844238d0
      cc( 4, 2) =        219.89311218d0
      cc( 5, 2) =       -110.71615601d0
      cc( 6, 2) =         24.82860374d0
      cc( 7, 2) =          3.39029098d0
      a ( 1, 2) =          0.67790002d0
      a ( 2, 2) =          0.79229999d0
      a ( 3, 2) =          1.04059994d0
      a ( 4, 2) =          1.44630003d0
      a ( 5, 2) =          1.96200001d0
      a ( 6, 2) =          2.62660003d0
      a ( 7, 2) =         10.83559990d0
      cc( 1, 3) =        139.57833862d0
      cc( 2, 3) =       -287.81640625d0
      cc( 3, 3) =        279.62557983d0
      cc( 4, 3) =       -194.22854614d0
      cc( 5, 3) =        104.41959381d0
      cc( 6, 3) =         16.83250046d0
      cc( 7, 3) =          5.69084311d0
      a ( 1, 3) =          0.78560001d0
      a ( 2, 3) =          0.93400002d0
      a ( 3, 3) =          1.24940002d0
      a ( 4, 3) =          1.77209997d0
      a ( 5, 3) =          2.47129989d0
      a ( 6, 3) =          7.06279993d0
      a ( 7, 3) =          8.15489960d0
      cc( 1, 4) =          0.00048200d0
      cc( 2, 4) =         -0.17614000d0
      cc( 3, 4) =        -43.04225159d0
      cc( 4, 4) =        130.94032288d0
      cc( 5, 4) =       -109.38676453d0
      cc( 6, 4) =         13.35779476d0
      cc( 7, 4) =          6.97829580d0
      a ( 1, 4) =          0.09670000d0
      a ( 2, 4) =          0.24150001d0
      a ( 3, 4) =          2.06349993d0
      a ( 4, 4) =          2.50259995d0
      a ( 5, 4) =          3.20359993d0
      a ( 6, 4) =          5.27120018d0
      a ( 7, 4) =          3.53810000d0
c
      go to 120
c
c     niobium (Nb_CRENBS_ECP)
c
  41  continue
      cc( 1, 1) =         -1.18631494d0
      cc( 2, 1) =         -7.80844688d0
      cc( 3, 1) =        -31.78058052d0
      cc( 4, 1) =       -104.86203003d0
      cc( 5, 1) =        -27.51132584d0
      a ( 1, 1) =          0.65390003d0
      a ( 2, 1) =          1.74010003d0
      a ( 3, 1) =          6.60360003d0
      a ( 4, 1) =         22.66510010d0
      a ( 5, 1) =         86.96939850d0
      cc( 1, 2) =        -62.14817429d0
      cc( 2, 2) =        190.85089111d0
      cc( 3, 2) =       -243.21609497d0
      cc( 4, 2) =        221.00381470d0
      cc( 5, 2) =       -114.93121338d0
      cc( 6, 2) =         28.06257629d0
      cc( 7, 2) =          3.42669392d0
      a ( 1, 2) =          0.73820001d0
      a ( 2, 2) =          0.85650003d0
      a ( 3, 2) =          1.11399996d0
      a ( 4, 2) =          1.53129995d0
      a ( 5, 2) =          2.06570005d0
      a ( 6, 2) =          2.50189996d0
      a ( 7, 2) =         10.75329971d0
      cc( 1, 3) =        151.04595947d0
      cc( 2, 3) =       -307.14862061d0
      cc( 3, 3) =        293.11035156d0
      cc( 4, 3) =       -199.33251953d0
      cc( 5, 3) =        112.10955811d0
      cc( 6, 3) =         16.48614311d0
      cc( 7, 3) =          5.68367910d0
      a ( 1, 3) =          0.88760000d0
      a ( 2, 3) =          1.06459999d0
      a ( 3, 3) =          1.44579995d0
      a ( 4, 3) =          2.08999991d0
      a ( 5, 3) =          2.99850011d0
      a ( 6, 3) =          8.67920017d0
      a ( 7, 3) =          8.90999985d0
      cc( 1, 4) =          0.23683000d0
      cc( 2, 4) =         -0.23673899d0
      cc( 3, 4) =         21.64477730d0
      cc( 4, 4) =        -54.17031479d0
      cc( 5, 4) =         41.12133026d0
      cc( 6, 4) =         11.94430923d0
      cc( 7, 4) =          7.14471006d0
      a ( 1, 4) =          0.64829999d0
      a ( 2, 4) =          0.28490001d0
      a ( 3, 4) =          1.36430001d0
      a ( 4, 4) =          1.46120000d0
      a ( 5, 4) =          1.66880000d0
      a ( 6, 4) =          7.64370012d0
      a ( 7, 4) =          6.68480015d0
c
      go to 120
c
c     molybdenum (Mo_CRENBS_ECP)
c
  42  continue
      cc( 1, 1) =         -0.97851598d0
      cc( 2, 1) =         -7.12729788d0
      cc( 3, 1) =        -24.65216255d0
      cc( 4, 1) =        -89.41860199d0
      cc( 5, 1) =        -25.30689621d0
      a ( 1, 1) =          0.67839998d0
      a ( 2, 1) =          1.73590004d0
      a ( 3, 1) =          5.99490023d0
      a ( 4, 1) =         19.34300041d0
      a ( 5, 1) =         66.90019989d0
      cc( 1, 2) =        -67.00737000d0
      cc( 2, 2) =        212.17037964d0
      cc( 3, 2) =       -271.79238892d0
      cc( 4, 2) =        243.67256165d0
      cc( 5, 2) =       -122.49226379d0
      cc( 6, 2) =         27.53581619d0
      cc( 7, 2) =          3.30936599d0
      a ( 1, 2) =          0.83120000d0
      a ( 2, 2) =          0.97170001d0
      a ( 3, 2) =          1.28129995d0
      a ( 4, 2) =          1.79719996d0
      a ( 5, 2) =          2.47049999d0
      a ( 6, 2) =          3.40219998d0
      a ( 7, 2) =         14.03419971d0
      cc( 1, 3) =        -47.97874451d0
      cc( 2, 3) =        172.93215942d0
      cc( 3, 3) =       -230.60501099d0
      cc( 4, 3) =        133.04591370d0
      cc( 5, 3) =        -23.99546814d0
      cc( 6, 3) =         25.51008034d0
      cc( 7, 3) =          5.34964418d0
      a ( 1, 3) =          0.66720003d0
      a ( 2, 3) =          0.76480001d0
      a ( 3, 3) =          0.94209999d0
      a ( 4, 3) =          1.16559994d0
      a ( 5, 3) =          2.48569989d0
      a ( 6, 3) =          2.93589997d0
      a ( 7, 3) =          9.79469967d0
      cc( 1, 4) =         -0.00025400d0
      cc( 2, 4) =         -0.25356400d0
      cc( 3, 4) =        -49.08042908d0
      cc( 4, 4) =        143.46177673d0
      cc( 5, 4) =       -139.96762085d0
      cc( 6, 4) =         13.82697392d0
      cc( 7, 4) =          6.99066401d0
      a ( 1, 4) =          0.08900000d0
      a ( 2, 4) =          0.31630000d0
      a ( 3, 4) =          2.15129995d0
      a ( 4, 4) =          2.93219996d0
      a ( 5, 4) =          3.74819994d0
      a ( 6, 4) =          5.98320007d0
      a ( 7, 4) =          1.19649994d0
c
      go to 120
c
c     technetium (Tc_CRENBS_ECP)
c
  43  continue
      cc( 1, 1) =         -0.79472297d0
      cc( 2, 1) =         -6.37796688d0
      cc( 3, 1) =        -18.39558029d0
      cc( 4, 1) =        -75.56327820d0
      cc( 5, 1) =        -23.71250725d0
      a ( 1, 1) =          0.70050001d0
      a ( 2, 1) =          1.72379994d0
      a ( 3, 1) =          5.29659987d0
      a ( 4, 1) =         16.52239990d0
      a ( 5, 1) =         52.40219879d0
      cc( 1, 2) =        -65.76921082d0
      cc( 2, 2) =        210.01269531d0
      cc( 3, 2) =       -267.63656616d0
      cc( 4, 2) =        243.09695435d0
      cc( 5, 2) =       -123.48133850d0
      cc( 6, 2) =         29.38408089d0
      cc( 7, 2) =          3.36798000d0
      a ( 1, 2) =          0.89370000d0
      a ( 2, 2) =          1.04149997d0
      a ( 3, 2) =          1.36290002d0
      a ( 4, 2) =          1.89919996d0
      a ( 5, 2) =          2.60130000d0
      a ( 6, 2) =          3.41840005d0
      a ( 7, 2) =         13.77740002d0
      cc( 1, 3) =        -62.08471680d0
      cc( 2, 3) =        202.83374023d0
      cc( 3, 3) =       -269.54275513d0
      cc( 4, 3) =        247.41317749d0
      cc( 5, 3) =       -117.69210052d0
      cc( 6, 3) =         25.42036819d0
      cc( 7, 3) =          5.31698513d0
      a ( 1, 3) =          0.76109999d0
      a ( 2, 3) =          0.87000000d0
      a ( 3, 3) =          1.11259997d0
      a ( 4, 3) =          1.52359998d0
      a ( 5, 3) =          1.97570002d0
      a ( 6, 3) =          2.91269994d0
      a ( 7, 3) =         10.04749966d0
      cc( 1, 4) =         -0.14150600d0
      cc( 2, 4) =         -7.07095194d0
      cc( 3, 4) =         38.74472809d0
      cc( 4, 4) =        -74.01949310d0
      cc( 5, 4) =         56.81457138d0
      cc( 6, 4) =         16.01659393d0
      cc( 7, 4) =          6.78246212d0
      a ( 1, 4) =          0.28979999d0
      a ( 2, 4) =          0.90539998d0
      a ( 3, 4) =          1.20210004d0
      a ( 4, 4) =          1.52540004d0
      a ( 5, 4) =          1.98080003d0
      a ( 6, 4) =         11.47700024d0
      a ( 7, 4) =         10.15670013d0
c
      go to 120
c
c     ruthenium (Ru_CRENBS_ECP)
c
  44  continue
      cc( 1, 1) =         -0.84380198d0
      cc( 2, 1) =         -6.69124413d0
      cc( 3, 1) =        -19.53392792d0
      cc( 4, 1) =        -79.45478821d0
      cc( 5, 1) =        -23.90557671d0
      a ( 1, 1) =          0.76889998d0
      a ( 2, 1) =          1.88810003d0
      a ( 3, 1) =          5.83479977d0
      a ( 4, 1) =         18.12579918d0
      a ( 5, 1) =         57.63710022d0
      cc( 1, 2) =        -71.56522369d0
      cc( 2, 2) =        228.00354004d0
      cc( 3, 2) =       -284.92181396d0
      cc( 4, 2) =        254.73445129d0
      cc( 5, 2) =       -129.97424316d0
      cc( 6, 2) =         27.73659134d0
      cc( 7, 2) =          3.39043093d0
      a ( 1, 2) =          0.99110001d0
      a ( 2, 2) =          1.14160001d0
      a ( 3, 2) =          1.46480000d0
      a ( 4, 2) =          1.99300003d0
      a ( 5, 2) =          2.64910007d0
      a ( 6, 2) =          3.19720006d0
      a ( 7, 2) =         12.94449997d0
      cc( 1, 3) =        -29.94395828d0
      cc( 2, 3) =        113.42784882d0
      cc( 3, 3) =       -152.27215576d0
      cc( 4, 3) =        142.09107971d0
      cc( 5, 3) =        -67.73392487d0
      cc( 6, 3) =         27.24661064d0
      cc( 7, 3) =          5.29189920d0
      a ( 1, 3) =          0.79670000d0
      a ( 2, 3) =          0.92159998d0
      a ( 3, 3) =          1.20869994d0
      a ( 4, 3) =          1.71340001d0
      a ( 5, 3) =          2.42820001d0
      a ( 6, 3) =          3.44860005d0
      a ( 7, 3) =         11.52260017d0
      cc( 1, 4) =         -0.01948400d0
      cc( 2, 4) =         -0.35189101d0
      cc( 3, 4) =         20.84566689d0
      cc( 4, 4) =        -59.39392853d0
      cc( 5, 4) =         53.66986465d0
      cc( 6, 4) =         14.55651188d0
      cc( 7, 4) =          6.93632078d0
      a ( 1, 4) =          0.23970000d0
      a ( 2, 4) =          0.42269999d0
      a ( 3, 4) =          1.51569998d0
      a ( 4, 4) =          1.81130004d0
      a ( 5, 4) =          2.25309992d0
      a ( 6, 4) =         11.23980045d0
      a ( 7, 4) =          9.81809998d0
c
      go to 120
c
c     rhodium (Rh_CRENBS_ECP)
c
  45  continue
      cc( 1, 1) =         -0.92449600d0
      cc( 2, 1) =         -7.14187098d0
      cc( 3, 1) =        -21.23577690d0
      cc( 4, 1) =        -85.08799744d0
      cc( 5, 1) =        -24.24301338d0
      a ( 1, 1) =          0.84579998d0
      a ( 2, 1) =          2.08389997d0
      a ( 3, 1) =          6.52729988d0
      a ( 4, 1) =         20.22669983d0
      a ( 5, 1) =         64.99040222d0
      cc( 1, 2) =        -77.58691406d0
      cc( 2, 2) =        240.80317688d0
      cc( 3, 2) =       -278.74105835d0
      cc( 4, 2) =        215.33572388d0
      cc( 5, 2) =       -114.76507568d0
      cc( 6, 2) =         36.79577255d0
      cc( 7, 2) =          3.33553004d0
      a ( 1, 2) =          1.06579995d0
      a ( 2, 2) =          1.21669996d0
      a ( 3, 2) =          1.53530002d0
      a ( 4, 2) =          2.08360004d0
      a ( 5, 2) =          2.96230006d0
      a ( 6, 2) =          3.09660006d0
      a ( 7, 2) =         15.23849964d0
      cc( 1, 3) =        -59.61508942d0
      cc( 2, 3) =        208.17607117d0
      cc( 3, 3) =       -285.84307861d0
      cc( 4, 3) =        259.17666626d0
      cc( 5, 3) =       -118.67604828d0
      cc( 6, 3) =         28.71784782d0
      cc( 7, 3) =          5.17818499d0
      a ( 1, 3) =          0.89630002d0
      a ( 2, 3) =          1.03649998d0
      a ( 3, 3) =          1.34590006d0
      a ( 4, 3) =          1.86530006d0
      a ( 5, 3) =          2.50020003d0
      a ( 6, 3) =          3.80550003d0
      a ( 7, 3) =         13.41059971d0
      cc( 1, 4) =         -0.35453799d0
      cc( 2, 4) =         -0.03887800d0
      cc( 3, 4) =        -83.94521332d0
      cc( 4, 4) =        205.58543396d0
      cc( 5, 4) =       -192.86888123d0
      cc( 6, 4) =         18.63974380d0
      cc( 7, 4) =          6.65249395d0
      a ( 1, 4) =          0.48629999d0
      a ( 2, 4) =          0.26460001d0
      a ( 3, 4) =          2.98650002d0
      a ( 4, 4) =          3.98250008d0
      a ( 5, 4) =          5.47550011d0
      a ( 6, 4) =          8.45230007d0
      a ( 7, 4) =          1.47650003d0
 120  continue
      n ( 1, 1) =   2
      n ( 2, 1) =   2
      n ( 3, 1) =   2
      n ( 4, 1) =   2
      n ( 5, 1) =   1
      nt( 1)    =   5
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      n ( 4, 2) =   2
      n ( 5, 2) =   2
      n ( 6, 2) =   1
      n ( 7, 2) =   0
      nt( 2)    =   7
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      n ( 3, 3) =   2
      n ( 4, 3) =   2
      n ( 5, 3) =   2
      n ( 6, 3) =   1
      n ( 7, 3) =   0
      nt( 3)    =   7
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      n ( 3, 4) =   2
      n ( 4, 4) =   2
      n ( 5, 4) =   2
      n ( 6, 4) =   1
      n ( 7, 4) =   0
      nt( 4)    =   7
c
      return
c
c     palladium (Pd_CRENBS_ECP)
c
  46  continue
      cc( 1, 1) =         -1.46114802d0
      cc( 2, 1) =         -9.50723076d0
      cc( 3, 1) =        -41.77036667d0
      cc( 4, 1) =        -21.68132210d0
      a ( 1, 1) =          0.98809999d0
      a ( 2, 1) =          2.62940001d0
      a ( 3, 1) =         10.25819969d0
      a ( 4, 1) =         29.55920029d0
      cc( 1, 2) =        170.32637024d0
      cc( 2, 2) =       -308.21520996d0
      cc( 3, 2) =        278.71969604d0
      cc( 4, 2) =       -136.46862793d0
      cc( 5, 2) =         32.85005951d0
      cc( 6, 2) =          3.32001805d0
      a ( 1, 2) =          1.56780005d0
      a ( 2, 2) =          1.89289999d0
      a ( 3, 2) =          2.56369996d0
      a ( 4, 2) =          3.51640010d0
      a ( 5, 2) =          4.76160002d0
      a ( 6, 2) =         18.32069969d0
      cc( 1, 3) =        104.99304199d0
      cc( 2, 3) =       -201.67398071d0
      cc( 3, 3) =        191.11369324d0
      cc( 4, 3) =        -89.62584686d0
      cc( 5, 3) =         28.29756546d0
      cc( 6, 3) =          5.31397009d0
      a ( 1, 3) =          1.10739994d0
      a ( 2, 3) =          1.34440005d0
      a ( 3, 3) =          1.82319999d0
      a ( 4, 3) =          2.47819996d0
      a ( 5, 3) =          3.49810004d0
      a ( 6, 3) =         12.45730019d0
      cc( 1, 4) =         -0.06195100d0
      cc( 2, 4) =         -1.10503602d0
      cc( 3, 4) =        -16.16929626d0
      cc( 4, 4) =         17.75333786d0
      cc( 5, 4) =          6.82440615d0
      cc( 6, 4) =          7.95571184d0
      a ( 1, 4) =          0.25200000d0
      a ( 2, 4) =          0.69040000d0
      a ( 3, 4) =          1.75610006d0
      a ( 4, 4) =          2.83060002d0
      a ( 5, 4) =          1.04600000d0
      a ( 6, 4) =         10.49960041d0
c
      go to 130
c
c     silver (Ag_CRENBS_ECP)
c
  47  continue
      cc( 1, 1) =         -1.54570901d0
      cc( 2, 1) =         -9.91327572d0
      cc( 3, 1) =        -44.19367218d0
      cc( 4, 1) =        -21.82627869d0
      a ( 1, 1) =          1.07079995d0
      a ( 2, 1) =          2.85200000d0
      a ( 3, 1) =         11.17080021d0
      a ( 4, 1) =         32.34180069d0
      cc( 1, 2) =        171.23280334d0
      cc( 2, 2) =       -307.32479858d0
      cc( 3, 2) =        280.42736816d0
      cc( 4, 2) =       -153.77388000d0
      cc( 5, 2) =         40.92211151d0
      cc( 6, 2) =          3.48360991d0
      a ( 1, 2) =          1.65400004d0
      a ( 2, 2) =          1.99160004d0
      a ( 3, 2) =          2.69039989d0
      a ( 4, 2) =          3.73289990d0
      a ( 5, 2) =          4.36100006d0
      a ( 6, 2) =         18.48649979d0
      cc( 1, 3) =        123.84740448d0
      cc( 2, 3) =       -238.36637878d0
      cc( 3, 3) =        222.44552612d0
      cc( 4, 3) =       -103.42056274d0
      cc( 5, 3) =         30.84042740d0
      cc( 6, 3) =          5.18727016d0
      a ( 1, 3) =          1.21969998d0
      a ( 2, 3) =          1.48839998d0
      a ( 3, 3) =          2.04329991d0
      a ( 4, 3) =          2.82279992d0
      a ( 5, 3) =          4.09789991d0
      a ( 6, 3) =         15.14739990d0
      cc( 1, 4) =         -0.09446900d0
      cc( 2, 4) =         -3.48698401d0
      cc( 3, 4) =        -38.78794098d0
      cc( 4, 4) =         34.73471451d0
      cc( 5, 4) =          9.89155388d0
      cc( 6, 4) =          7.72656393d0
      a ( 1, 4) =          0.29390001d0
      a ( 2, 4) =          0.92500001d0
      a ( 3, 4) =          2.09789991d0
      a ( 4, 4) =          2.52940011d0
      a ( 5, 4) =          1.06369996d0
      a ( 6, 4) =         10.93939972d0
c
      go to 130
c
c     cadmium (Cd_CRENBS_ECP)
c
  48  continue
      cc( 1, 1) =         -1.77006102d0
      cc( 2, 1) =        -10.56859684d0
      cc( 3, 1) =        -49.03567505d0
      cc( 4, 1) =        -22.09135056d0
      a ( 1, 1) =          1.20790005d0
      a ( 2, 1) =          3.17670012d0
      a ( 3, 1) =         12.61629963d0
      a ( 4, 1) =         36.90489960d0
      cc( 1, 2) =         31.95877838d0
      cc( 2, 2) =       -104.65022278d0
      cc( 3, 2) =        173.25119019d0
      cc( 4, 2) =       -104.91418457d0
      cc( 5, 2) =         23.64106560d0
      cc( 6, 2) =          3.90806794d0
      a ( 1, 2) =          0.80620003d0
      a ( 2, 2) =          0.93199998d0
      a ( 3, 2) =          1.18079996d0
      a ( 4, 2) =          1.52180004d0
      a ( 5, 2) =          2.08349991d0
      a ( 6, 2) =         11.30650043d0
      cc( 1, 3) =        -41.68565369d0
      cc( 2, 3) =        145.67808533d0
      cc( 3, 3) =       -189.84249878d0
      cc( 4, 3) =        125.71853638d0
      cc( 5, 3) =         31.13378716d0
      cc( 6, 3) =          5.27771378d0
      a ( 1, 3) =          0.94749999d0
      a ( 2, 3) =          1.10239995d0
      a ( 3, 3) =          1.44000006d0
      a ( 4, 3) =          1.90450001d0
      a ( 5, 3) =          7.30249977d0
      a ( 6, 3) =         19.80590057d0
      cc( 1, 4) =         -0.52976000d0
      cc( 2, 4) =         -5.00963783d0
      cc( 3, 4) =        -15.28431988d0
      cc( 4, 4) =          6.59910011d0
      cc( 5, 4) =         11.88312626d0
      cc( 6, 4) =          7.84789896d0
      a ( 1, 4) =          0.54339999d0
      a ( 2, 4) =          1.57570004d0
      a ( 3, 4) =          1.82249999d0
      a ( 4, 4) =          2.54139996d0
      a ( 5, 4) =          1.30009997d0
      a ( 6, 4) =         10.56420040d0
c
 130  n ( 1, 1) =   2
      n ( 2, 1) =   2
      n ( 3, 1) =   2
      n ( 4, 1) =   1
      nt( 1)    =   4
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      n ( 4, 2) =   2
      n ( 5, 2) =   1
      n ( 6, 2) =   0
      nt( 2)    =   6
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      n ( 3, 3) =   2
      n ( 4, 3) =   2
      n ( 5, 3) =   1
      n ( 6, 3) =   0
      nt( 3)    =   6
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      n ( 3, 4) =   2
      n ( 4, 4) =   2
      n ( 5, 4) =   1
      n ( 6, 4) =   0
      nt( 4)    =   6
c
      return
c
 6010 format (1x,'*** data not found for CRENBS ecp type ',a4)
      end
**==crenbs2.f
      subroutine crenbs2 (zinp,lmx,nelcor,nt,n,cc,a)
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
      implicit character * 8 (z)
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
      dimension nt(6), n(10,6), cc(10,6), a(10,6)
c
      dimension ztit(15)
      data maxtyp/15/
      data ztit/
     $ 'la','hf',
     $ 'ta','w ','re','os','ir','pt','au','hg','pb',
     + 'bi','po','at','rn'/
c
c     this routine contains a library of ecp data taken from
c     various literature sources.
c     3rd transition series sjkb potentials
c     
      do 10 ityp = 1 , maxtyp
         if (zinp.eq.ztit(ityp)) go to 20
 10   continue
      write (iwr,6010) zinp
      call caserr2('no library available for requested ecp')
c
 20   go to(57,72,73,74,75,76,77,78,79,80,
     +      82,83,84,85,86), ityp
c
c     lanthanum (La_CRENBS_ECP)
c
  57  continue
      cc( 1, 1) =         -0.22267400d0
      cc( 2, 1) =         -2.14405394d0
      cc( 3, 1) =         -8.95630169d0
      cc( 4, 1) =        -31.27042389d0
      cc( 5, 1) =        -80.72589874d0
      cc( 6, 1) =        -34.54819107d0
      a ( 1, 1) =          0.20420000d0
      a ( 2, 1) =          0.53509998d0
      a ( 3, 1) =          1.43079996d0
      a ( 4, 1) =          4.21190023d0
      a ( 5, 1) =         13.20160007d0
      a ( 6, 1) =         39.32419968d0
      cc( 1, 2) =        -57.20286942d0
      cc( 2, 2) =        190.33511353d0
      cc( 3, 2) =       -265.08895874d0
      cc( 4, 2) =        250.97280884d0
      cc( 5, 2) =       -178.21974182d0
      cc( 6, 2) =        103.27313995d0
      cc( 7, 2) =         34.95799255d0
      cc( 8, 2) =          6.50971413d0
      a ( 1, 2) =          0.49239999d0
      a ( 2, 2) =          0.56680000d0
      a ( 3, 2) =          0.73250002d0
      a ( 4, 2) =          1.01779997d0
      a ( 5, 2) =          1.46050000d0
      a ( 6, 2) =          2.04999995d0
      a ( 7, 2) =          6.22079992d0
      a ( 8, 2) =         17.26569939d0
      cc( 1, 3) =        -28.08255768d0
      cc( 2, 3) =        111.27140045d0
      cc( 3, 3) =       -172.92143250d0
      cc( 4, 3) =        176.78236389d0
      cc( 5, 3) =       -131.83255005d0
      cc( 6, 3) =         73.78304291d0
      cc( 7, 3) =         29.54417992d0
      cc( 8, 3) =          6.00691700d0
      a ( 1, 3) =          0.37990001d0
      a ( 2, 3) =          0.42660001d0
      a ( 3, 3) =          0.53280002d0
      a ( 4, 3) =          0.71590000d0
      a ( 5, 3) =          0.99599999d0
      a ( 6, 3) =          1.34940004d0
      a ( 7, 3) =          4.17630005d0
      a ( 8, 3) =         11.68099976d0
      cc( 1, 4) =         -0.03236600d0
      cc( 2, 4) =        -24.12875366d0
      cc( 3, 4) =         72.46633911d0
      cc( 4, 4) =       -125.80812836d0
      cc( 5, 4) =        177.71560669d0
      cc( 6, 4) =       -112.06113434d0
      cc( 7, 4) =         34.28693390d0
      cc( 8, 4) =          7.13268280d0
      a ( 1, 4) =          0.11310000d0
      a ( 2, 4) =          0.59270000d0
      a ( 3, 4) =          0.68229997d0
      a ( 4, 4) =          0.86989999d0
      a ( 5, 4) =          1.17309999d0
      a ( 6, 4) =          1.51779997d0
      a ( 7, 4) =          2.46650004d0
      a ( 8, 4) =         12.28919983d0
      n ( 1, 1) =   2
      n ( 2, 1) =   2
      n ( 3, 1) =   2
      n ( 4, 1) =   2
      n ( 5, 1) =   2
      n ( 6, 1) =   1
      nt( 1)    =   6
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      n ( 4, 2) =   2
      n ( 5, 2) =   2
      n ( 6, 2) =   2
      n ( 7, 2) =   1
      n ( 8, 2) =   0
      nt( 2)    =   8
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      n ( 3, 3) =   2
      n ( 4, 3) =   2
      n ( 5, 3) =   2
      n ( 6, 3) =   2
      n ( 7, 3) =   1
      n ( 8, 3) =   0
      nt( 3)    =   8
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      n ( 3, 4) =   2
      n ( 4, 4) =   2
      n ( 5, 4) =   2
      n ( 6, 4) =   2
      n ( 7, 4) =   1
      n ( 8, 4) =   0
      nt( 4)    =   8
      nelcor    =  54
      lmx       =   3
c
      return
c
c     hafnium (Hf_CRENBS_ECP)
c
  72  continue
      cc( 1, 1) =         -0.63135999d0
      cc( 2, 1) =         -7.10547495d0
      cc( 3, 1) =        -24.14714623d0
      cc( 4, 1) =        -85.58006287d0
      cc( 5, 1) =       -185.18806458d0
      cc( 6, 1) =        -47.48435974d0
      a ( 1, 1) =          0.56110001d0
      a ( 2, 1) =          1.35570002d0
      a ( 3, 1) =          3.52559996d0
      a ( 4, 1) =          9.66549969d0
      a ( 5, 1) =         31.15810013d0
      a ( 6, 1) =         91.57109833d0
      cc( 1, 2) =         12.23064327d0
      cc( 2, 2) =        -59.60324478d0
      cc( 3, 2) =        156.65519714d0
      cc( 4, 2) =       -204.99865723d0
      cc( 5, 2) =        187.90031433d0
      cc( 6, 2) =        -82.80174255d0
      cc( 7, 2) =         40.81014252d0
      cc( 8, 2) =          6.65634918d0
      a ( 1, 2) =          0.62709999d0
      a ( 2, 2) =          0.70029998d0
      a ( 3, 2) =          0.87849998d0
      a ( 4, 2) =          1.20389998d0
      a ( 5, 2) =          1.72860003d0
      a ( 6, 2) =          2.40269995d0
      a ( 7, 2) =          4.02040005d0
      a ( 8, 2) =         12.75230026d0
      cc( 1, 3) =         37.69438171d0
      cc( 2, 3) =       -116.20706177d0
      cc( 3, 3) =        196.35423279d0
      cc( 4, 3) =       -211.90403748d0
      cc( 5, 3) =        113.47141266d0
      cc( 6, 3) =         47.26153183d0
      cc( 7, 3) =         45.12688446d0
      cc( 8, 3) =          5.77881813d0
      a ( 1, 3) =          0.42649999d0
      a ( 2, 3) =          0.47530001d0
      a ( 3, 3) =          0.57969999d0
      a ( 4, 3) =          0.74269998d0
      a ( 5, 3) =          0.91619998d0
      a ( 6, 3) =          3.80509996d0
      a ( 7, 3) =          8.30080032d0
      a ( 8, 3) =         23.66410065d0
      cc( 1, 4) =         -0.24424000d0
      cc( 2, 4) =         11.36239910d0
      cc( 3, 4) =        -36.05837250d0
      cc( 4, 4) =         65.08800507d0
      cc( 5, 4) =        -97.13800049d0
      cc( 6, 4) =         75.79869080d0
      cc( 7, 4) =         19.08870888d0
      cc( 8, 4) =          8.35774994d0
      a ( 1, 4) =          0.25749999d0
      a ( 2, 4) =          0.56739998d0
      a ( 3, 4) =          0.67530000d0
      a ( 4, 4) =          0.88800001d0
      a ( 5, 4) =          1.21080005d0
      a ( 6, 4) =          1.56429994d0
      a ( 7, 4) =          6.02829981d0
      a ( 8, 4) =          5.42290020d0
      cc( 1, 5) =         -6.32602882d0
      cc( 2, 5) =         12.38346767d0
      cc( 3, 5) =        -67.87256622d0
      cc( 4, 5) =         70.92347717d0
      cc( 5, 5) =         36.17158890d0
      cc( 6, 5) =         96.39685059d0
      cc( 7, 5) =         28.85315132d0
      cc( 8, 5) =          1.03963399d0
      a ( 1, 5) =          0.60949999d0
      a ( 2, 5) =          0.68379998d0
      a ( 3, 5) =          1.14269996d0
      a ( 4, 5) =          1.23740005d0
      a ( 5, 5) =          4.84910011d0
      a ( 6, 5) =         12.74800015d0
      a ( 7, 5) =         33.74710083d0
      a ( 8, 5) =         33.60300064d0
      nelcor    =  68
c
      go to 30
c
c     tantalum (Ta_CRENBS_ECP)
c
  73  continue
      cc( 1, 1) =         -0.71255600d0
      cc( 2, 1) =         -7.58646679d0
      cc( 3, 1) =        -25.98602867d0
      cc( 4, 1) =        -91.30485535d0
      cc( 5, 1) =       -202.25379944d0
      cc( 6, 1) =        -48.06631088d0
      a ( 1, 1) =          0.62690002d0
      a ( 2, 1) =          1.50860000d0
      a ( 3, 1) =          3.94729996d0
      a ( 4, 1) =         10.69579983d0
      a ( 5, 1) =         35.08679962d0
      a ( 6, 1) =        103.43869781d0
      cc( 1, 2) =         12.13485432d0
      cc( 2, 2) =        -59.38967896d0
      cc( 3, 2) =        160.96141052d0
      cc( 4, 2) =       -207.03143311d0
      cc( 5, 2) =        187.82365417d0
      cc( 6, 2) =        -85.49249268d0
      cc( 7, 2) =         45.54972458d0
      cc( 8, 2) =          6.45374203d0
      a ( 1, 2) =          0.68239999d0
      a ( 2, 2) =          0.76580000d0
      a ( 3, 2) =          0.97399998d0
      a ( 4, 2) =          1.35800004d0
      a ( 5, 2) =          2.00209999d0
      a ( 6, 2) =          2.90079999d0
      a ( 7, 2) =          4.68419981d0
      a ( 8, 2) =         16.76090050d0
      cc( 1, 3) =         38.24298859d0
      cc( 2, 3) =       -117.84931946d0
      cc( 3, 3) =        203.91978455d0
      cc( 4, 3) =       -223.77775574d0
      cc( 5, 3) =        121.21456146d0
      cc( 6, 3) =         57.04997253d0
      cc( 7, 3) =         46.13226318d0
      cc( 8, 3) =          5.67695284d0
      a ( 1, 3) =          0.46160001d0
      a ( 2, 3) =          0.51789999d0
      a ( 3, 3) =          0.64080000d0
      a ( 4, 3) =          0.83700001d0
      a ( 5, 3) =          1.05239999d0
      a ( 6, 3) =          4.61880016d0
      a ( 7, 3) =          9.97859955d0
      a ( 8, 3) =         28.19230080d0
      cc( 1, 4) =         -0.11714300d0
      cc( 2, 4) =        -28.52329445d0
      cc( 3, 4) =         74.39411926d0
      cc( 4, 4) =       -114.17246246d0
      cc( 5, 4) =        155.09744263d0
      cc( 6, 4) =        -83.56156921d0
      cc( 7, 4) =         47.67310333d0
      cc( 8, 4) =          7.06498909d0
      a ( 1, 4) =          0.21950001d0
      a ( 2, 4) =          1.03610003d0
      a ( 3, 4) =          1.21870005d0
      a ( 4, 4) =          1.62769997d0
      a ( 5, 4) =          2.29290009d0
      a ( 6, 4) =          3.11770010d0
      a ( 7, 4) =          5.81510019d0
      a ( 8, 4) =         23.04089928d0
      cc( 1, 5) =         22.67155838d0
      cc( 2, 5) =        -63.55398560d0
      cc( 3, 5) =         73.78347778d0
      cc( 4, 5) =        -36.86455154d0
      cc( 5, 5) =         16.23764420d0
      cc( 6, 5) =         38.82068634d0
      cc( 7, 5) =         28.87883759d0
      cc( 8, 5) =          1.53800905d0
      a ( 1, 5) =          0.47920001d0
      a ( 2, 5) =          0.53090000d0
      a ( 3, 5) =          0.63300002d0
      a ( 4, 5) =          0.78230000d0
      a ( 5, 5) =          1.46529996d0
      a ( 6, 5) =          4.93750000d0
      a ( 7, 5) =         10.75310040d0
      a ( 8, 5) =         28.31669998d0
      nelcor    =  68
c
      go to 30
c
c     tungsten (W_CRENBS_ECP)
c
  74  continue
      cc( 1, 1) =         -0.75037998d0
      cc( 2, 1) =         -8.05651283d0
      cc( 3, 1) =        -27.77778625d0
      cc( 4, 1) =        -96.63814545d0
      cc( 5, 1) =       -219.25021362d0
      cc( 6, 1) =        -48.63983536d0
      a ( 1, 1) =          0.68370003d0
      a ( 2, 1) =          1.65279996d0
      a ( 3, 1) =          4.38019991d0
      a ( 4, 1) =         11.73820019d0
      a ( 5, 1) =         39.06900024d0
      a ( 6, 1) =        115.77369690d0
      cc( 1, 2) =         -0.00559700d0
      cc( 2, 2) =        -57.14800644d0
      cc( 3, 2) =        168.47732544d0
      cc( 4, 2) =       -202.82611084d0
      cc( 5, 2) =        127.62249756d0
      cc( 6, 2) =       -116.13665771d0
      cc( 7, 2) =         68.52722931d0
      cc( 8, 2) =          4.92766380d0
      a ( 1, 2) =          0.10920000d0
      a ( 2, 2) =          0.87210000d0
      a ( 3, 2) =          1.01670003d0
      a ( 4, 2) =          1.35259998d0
      a ( 5, 2) =          1.81850004d0
      a ( 6, 2) =          2.72230005d0
      a ( 7, 2) =          1.86889994d0
      a ( 8, 2) =         19.64010048d0
      cc( 1, 3) =         38.40269089d0
      cc( 2, 3) =       -116.08489227d0
      cc( 3, 3) =        199.11714172d0
      cc( 4, 3) =       -217.43078613d0
      cc( 5, 3) =        119.61345673d0
      cc( 6, 3) =         54.77985001d0
      cc( 7, 3) =         47.16429138d0
      cc( 8, 3) =          5.72915888d0
      a ( 1, 3) =          0.49710000d0
      a ( 2, 3) =          0.55409998d0
      a ( 3, 3) =          0.68030000d0
      a ( 4, 3) =          0.87900001d0
      a ( 5, 3) =          1.09570003d0
      a ( 6, 3) =          4.52099991d0
      a ( 7, 3) =          9.70989990d0
      a ( 8, 3) =         27.31539917d0
      cc( 1, 4) =        -10.73067856d0
      cc( 2, 4) =         -0.09588000d0
      cc( 3, 4) =         13.17100811d0
      cc( 4, 4) =        -74.29113007d0
      cc( 5, 4) =        194.93284607d0
      cc( 6, 4) =       -128.38191223d0
      cc( 7, 4) =         22.44917870d0
      cc( 8, 4) =          8.12122154d0
      a ( 1, 4) =          0.81410003d0
      a ( 2, 4) =          0.22220001d0
      a ( 3, 4) =          0.88480002d0
      a ( 4, 4) =          2.14289999d0
      a ( 5, 4) =          2.70149994d0
      a ( 6, 4) =          3.51130009d0
      a ( 7, 4) =          5.69469976d0
      a ( 8, 4) =          4.35500002d0
      cc( 1, 5) =         42.52615738d0
      cc( 2, 5) =       -130.23117065d0
      cc( 3, 5) =        176.63513184d0
      cc( 4, 5) =       -173.75531006d0
      cc( 5, 5) =        103.75617981d0
      cc( 6, 5) =         59.37061310d0
      cc( 7, 5) =         47.21155930d0
      cc( 8, 5) =          4.74469805d0
      a ( 1, 5) =          0.44949999d0
      a ( 2, 5) =          0.50650001d0
      a ( 3, 5) =          0.62580001d0
      a ( 4, 5) =          0.85600001d0
      a ( 5, 5) =          1.06809998d0
      a ( 6, 5) =          4.84870005d0
      a ( 7, 5) =         11.44680023d0
      a ( 8, 5) =         37.23880005d0
      nelcor    =  68
c
      go to 30
c
c     rhenium (Re_CRENBS_ECP)
c
  75  continue
      cc( 1, 1) =         -0.80902398d0
      cc( 2, 1) =         -8.50144196d0
      cc( 3, 1) =        -29.58082199d0
      cc( 4, 1) =       -101.49237061d0
      cc( 5, 1) =       -235.57843018d0
      cc( 6, 1) =        -49.22165298d0
      a ( 1, 1) =          0.74460000d0
      a ( 2, 1) =          1.80200005d0
      a ( 3, 1) =          4.82390022d0
      a ( 4, 1) =         12.78439999d0
      a ( 5, 1) =         42.99539948d0
      a ( 6, 1) =        128.28930664d0
      cc( 1, 2) =         12.20130444d0
      cc( 2, 2) =        -58.66179276d0
      cc( 3, 2) =        154.63580322d0
      cc( 4, 2) =       -191.60073853d0
      cc( 5, 2) =        181.25675964d0
      cc( 6, 2) =        -75.58824921d0
      cc( 7, 2) =         44.00109863d0
      cc( 8, 2) =          6.71948910d0
      a ( 1, 2) =          0.76179999d0
      a ( 2, 2) =          0.85130000d0
      a ( 3, 2) =          1.06920004d0
      a ( 4, 2) =          1.46480000d0
      a ( 5, 2) =          2.10380006d0
      a ( 6, 2) =          2.93120003d0
      a ( 7, 2) =          4.98059988d0
      a ( 8, 2) =         14.67980003d0
      cc( 1, 3) =         38.36573791d0
      cc( 2, 3) =       -120.06738281d0
      cc( 3, 3) =        213.54537964d0
      cc( 4, 3) =       -236.23445129d0
      cc( 5, 3) =        130.94828796d0
      cc( 6, 3) =         66.50952911d0
      cc( 7, 3) =         48.34811783d0
      cc( 8, 3) =          5.64982986d0
      a ( 1, 3) =          0.53570002d0
      a ( 2, 3) =          0.60270000d0
      a ( 3, 3) =          0.74870002d0
      a ( 4, 3) =          0.98460001d0
      a ( 5, 3) =          1.24810004d0
      a ( 6, 3) =          5.42570019d0
      a ( 7, 3) =         11.56700039d0
      a ( 8, 3) =         31.93650055d0
      cc( 1, 4) =         -0.13630500d0
      cc( 2, 4) =         76.10538483d0
      cc( 3, 4) =        -25.86570549d0
      cc( 4, 4) =       -131.80529785d0
      cc( 5, 4) =        203.93190002d0
      cc( 6, 4) =       -130.43125916d0
      cc( 7, 4) =         53.73962021d0
      cc( 8, 4) =          6.90621901d0
      a ( 1, 4) =          0.26010001d0
      a ( 2, 4) =          1.42289996d0
      a ( 3, 4) =          1.18320000d0
      a ( 4, 4) =          1.91129994d0
      a ( 5, 4) =          2.73810005d0
      a ( 6, 4) =          3.77220011d0
      a ( 7, 4) =          6.40880013d0
      a ( 8, 4) =         28.94709969d0
      cc( 1, 5) =         42.60373306d0
      cc( 2, 5) =       -132.34577942d0
      cc( 3, 5) =        197.70394897d0
      cc( 4, 5) =       -187.96063232d0
      cc( 5, 5) =         99.53369141d0
      cc( 6, 5) =         50.98129654d0
      cc( 7, 5) =         43.48990631d0
      cc( 8, 5) =          4.84967518d0
      a ( 1, 5) =          0.45539999d0
      a ( 2, 5) =          0.51609999d0
      a ( 3, 5) =          0.64760000d0
      a ( 4, 5) =          0.85759997d0
      a ( 5, 5) =          1.09650004d0
      a ( 6, 5) =          4.66149998d0
      a ( 7, 5) =         10.15660000d0
      a ( 8, 5) =         28.48620033d0
      nelcor    =  68
      go to 30
c
c     osmium (Os_CRENBS_ECP)
c
  76  continue
      cc( 1, 1) =         -0.87657899d0
      cc( 2, 1) =         -8.94104481d0
      cc( 3, 1) =        -31.52030945d0
      cc( 4, 1) =       -106.30191803d0
      cc( 5, 1) =       -252.03785706d0
      cc( 6, 1) =        -49.82308578d0
      a ( 1, 1) =          0.80779999d0
      a ( 2, 1) =          1.95379996d0
      a ( 3, 1) =          5.28529978d0
      a ( 4, 1) =         13.87800026d0
      a ( 5, 1) =         47.06290054d0
      a ( 6, 1) =        141.49729919d0
      cc( 1, 2) =         26.27425575d0
      cc( 2, 2) =        -96.17611694d0
      cc( 3, 2) =        211.95298767d0
      cc( 4, 2) =       -253.41165161d0
      cc( 5, 2) =        228.25082397d0
      cc( 6, 2) =       -102.19163513d0
      cc( 7, 2) =         48.77917480d0
      cc( 8, 2) =          6.43857288d0
      a ( 1, 2) =          0.80199999d0
      a ( 2, 2) =          0.90719998d0
      a ( 3, 2) =          1.15779996d0
      a ( 4, 2) =          1.61629999d0
      a ( 5, 2) =          2.36520004d0
      a ( 6, 2) =          3.38779998d0
      a ( 7, 2) =          5.56780005d0
      a ( 8, 2) =         19.54560089d0
      cc( 1, 3) =         37.43206406d0
      cc( 2, 3) =       -116.31768799d0
      cc( 3, 3) =        207.20248413d0
      cc( 4, 3) =       -229.60458374d0
      cc( 5, 3) =        129.95959473d0
      cc( 6, 3) =         63.95362473d0
      cc( 7, 3) =         49.04730225d0
      cc( 8, 3) =          5.68902302d0
      a ( 1, 3) =          0.56970000d0
      a ( 2, 3) =          0.63859999d0
      a ( 3, 3) =          0.78920001d0
      a ( 4, 3) =          1.02970004d0
      a ( 5, 3) =          1.29589999d0
      a ( 6, 3) =          5.37249994d0
      a ( 7, 3) =         11.33450031d0
      a ( 8, 3) =         31.14130020d0
      cc( 1, 4) =         -0.25271201d0
      cc( 2, 4) =         13.53297424d0
      cc( 3, 4) =        -43.62031174d0
      cc( 4, 4) =         79.10832977d0
      cc( 5, 4) =       -120.06851196d0
      cc( 6, 4) =        100.24253845d0
      cc( 7, 4) =         51.20841217d0
      cc( 8, 4) =          7.07804298d0
      a ( 1, 4) =          0.31650001d0
      a ( 2, 4) =          0.78649998d0
      a ( 3, 4) =          0.92989999d0
      a ( 4, 4) =          1.21500003d0
      a ( 5, 4) =          1.65980005d0
      a ( 6, 4) =          2.16799998d0
      a ( 7, 4) =          8.58160019d0
      a ( 8, 4) =         30.71570015d0
      cc( 1, 5) =         42.01272964d0
      cc( 2, 5) =       -130.60057068d0
      cc( 3, 5) =        195.07528687d0
      cc( 4, 5) =       -184.75950623d0
      cc( 5, 5) =         97.17532349d0
      cc( 6, 5) =         47.79431534d0
      cc( 7, 5) =         44.16305923d0
      cc( 8, 5) =          4.86101103d0
      a ( 1, 5) =          0.42480001d0
      a ( 2, 5) =          0.48190001d0
      a ( 3, 5) =          0.60549998d0
      a ( 4, 5) =          0.80210000d0
      a ( 5, 5) =          1.02779996d0
      a ( 6, 5) =          4.36679983d0
      a ( 7, 5) =          9.82260036d0
      a ( 8, 5) =         27.93250084d0
      nelcor    = 68
      go to 30
c
c     iridium (Ir_CRENBS_ECP)
c
  77  continue
      cc( 1, 1) =         -0.91868699d0
      cc( 2, 1) =         -9.34183121d0
      cc( 3, 1) =        -33.08732224d0
      cc( 4, 1) =       -110.59931183d0
      cc( 5, 1) =       -266.87414551d0
      cc( 6, 1) =        -50.38343811d0
      a ( 1, 1) =          0.86720002d0
      a ( 2, 1) =          2.10159993d0
      a ( 3, 1) =          5.73199987d0
      a ( 4, 1) =         14.90610027d0
      a ( 5, 1) =         50.87039948d0
      a ( 6, 1) =        154.07629395d0
      cc( 1, 2) =         27.13777351d0
      cc( 2, 2) =        -97.93394470d0
      cc( 3, 2) =        212.03878784d0
      cc( 4, 2) =       -248.60784912d0
      cc( 5, 2) =        227.19543457d0
      cc( 6, 2) =        -97.36859131d0
      cc( 7, 2) =         47.97808456d0
      cc( 8, 2) =          6.55392122d0
      a ( 1, 2) =          0.84130001d0
      a ( 2, 2) =          0.95020002d0
      a ( 3, 2) =          1.20710003d0
      a ( 4, 2) =          1.67019999d0
      a ( 5, 2) =          2.41790009d0
      a ( 6, 2) =          3.39159989d0
      a ( 7, 2) =          5.76249981d0
      a ( 8, 2) =         18.47249985d0
      cc( 1, 3) =         57.66228867d0
      cc( 2, 3) =       -154.39692688d0
      cc( 3, 3) =        247.69586182d0
      cc( 4, 3) =       -276.29647827d0
      cc( 5, 3) =        270.38247681d0
      cc( 6, 3) =       -159.55268860d0
      cc( 7, 3) =         41.99956512d0
      cc( 8, 3) =          5.93545818d0
      a ( 1, 3) =          0.65630001d0
      a ( 2, 3) =          0.72520000d0
      a ( 3, 3) =          0.90340000d0
      a ( 4, 3) =          1.23109996d0
      a ( 5, 3) =          1.78830004d0
      a ( 6, 3) =          2.38120008d0
      a ( 7, 3) =          3.08200002d0
      a ( 8, 3) =         11.04860020d0
      cc( 1, 4) =         -0.12088400d0
      cc( 2, 4) =         -1.39111197d0
      cc( 3, 4) =         43.52083206d0
      cc( 4, 4) =       -142.59112549d0
      cc( 5, 4) =        247.56771851d0
      cc( 6, 4) =       -151.51245117d0
      cc( 7, 4) =         24.36254120d0
      cc( 8, 4) =          8.04463005d0
      a ( 1, 4) =          0.28540000d0
      a ( 2, 4) =          0.83289999d0
      a ( 3, 4) =          1.77460003d0
      a ( 4, 4) =          2.23830009d0
      a ( 5, 4) =          3.04940009d0
      a ( 6, 4) =          4.09219980d0
      a ( 7, 4) =          6.55800009d0
      a ( 8, 4) =          5.53870010d0
      cc( 1, 5) =        -19.30400848d0
      cc( 2, 5) =         60.73864746d0
      cc( 3, 5) =       -107.93820953d0
      cc( 4, 5) =        136.65959167d0
      cc( 5, 5) =       -107.83163452d0
      cc( 6, 5) =         58.89172745d0
      cc( 7, 5) =         29.42035294d0
      cc( 8, 5) =          5.93143177d0
      a ( 1, 5) =          0.30829999d0
      a ( 2, 5) =          0.34630001d0
      a ( 3, 5) =          0.42930001d0
      a ( 4, 5) =          0.57349998d0
      a ( 5, 5) =          0.79500002d0
      a ( 6, 5) =          1.15090001d0
      a ( 7, 5) =          4.22580004d0
      a ( 8, 5) =          9.36260033d0
      nelcor    = 68
c
      go to 30
c
c     platium (Pt_CRENBS_ECP)
c
  78  continue
      cc( 1, 1) =         -0.97492802d0
      cc( 2, 1) =         -9.49699306d0
      cc( 3, 1) =        -31.97232437d0
      cc( 4, 1) =       -115.02309418d0
      cc( 5, 1) =       -275.53890991d0
      cc( 6, 1) =        -50.80278015d0
      a ( 1, 1) =          0.91210002d0
      a ( 2, 1) =          2.16540003d0
      a ( 3, 1) =          6.06580019d0
      a ( 4, 1) =         15.52079964d0
      a ( 5, 1) =         52.94639969d0
      a ( 6, 1) =        160.88090515d0
      cc( 1, 2) =        -80.75402069d0
      cc( 2, 2) =        255.40243530d0
      cc( 3, 2) =       -314.58084106d0
      cc( 4, 2) =        297.73495483d0
      cc( 5, 2) =       -197.63934326d0
      cc( 6, 2) =        146.50839233d0
      cc( 7, 2) =         48.66068268d0
      cc( 8, 2) =          4.06491184d0
      a ( 1, 2) =          1.14330006d0
      a ( 2, 2) =          1.32480001d0
      a ( 3, 2) =          1.72399998d0
      a ( 4, 2) =          2.41059995d0
      a ( 5, 2) =          3.49970007d0
      a ( 6, 2) =          5.13000011d0
      a ( 7, 2) =         12.90919971d0
      a ( 8, 2) =         34.22660065d0
      cc( 1, 3) =        -17.78671646d0
      cc( 2, 3) =         96.91400146d0
      cc( 3, 3) =       -146.01452637d0
      cc( 4, 3) =        147.23194885d0
      cc( 5, 3) =        -89.58618164d0
      cc( 6, 3) =         85.01801300d0
      cc( 7, 3) =         50.10876846d0
      cc( 8, 3) =          5.64266300d0
      a ( 1, 3) =          0.85159999d0
      a ( 2, 3) =          0.97509998d0
      a ( 3, 3) =          1.26619995d0
      a ( 4, 3) =          1.80509996d0
      a ( 5, 3) =          2.70499992d0
      a ( 6, 3) =          4.24480009d0
      a ( 7, 3) =         10.74390030d0
      a ( 8, 3) =         31.21459961d0
      cc( 1, 4) =         -0.11851300d0
      cc( 2, 4) =        -17.32895088d0
      cc( 3, 4) =         49.30339050d0
      cc( 4, 4) =       -115.68026733d0
      cc( 5, 4) =        226.66456604d0
      cc( 6, 4) =       -154.83076477d0
      cc( 7, 4) =         59.63909912d0
      cc( 8, 4) =          6.81120205d0
      a ( 1, 4) =          0.28200001d0
      a ( 2, 4) =          1.26719999d0
      a ( 3, 4) =          1.56920004d0
      a ( 4, 4) =          2.30189991d0
      a ( 5, 4) =          3.29360008d0
      a ( 6, 4) =          4.60909987d0
      a ( 7, 4) =          7.73649979d0
      a ( 8, 4) =         36.09999847d0
      cc( 1, 5) =        -28.15009499d0
      cc( 2, 5) =         31.23972702d0
      cc( 3, 5) =         47.82197189d0
      cc( 4, 5) =        -29.59718513d0
      cc( 5, 5) =         68.58837128d0
      cc( 6, 5) =        -27.41664124d0
      cc( 7, 5) =         44.57673645d0
      cc( 8, 5) =          4.99491882d0
      a ( 1, 5) =          0.32130000d0
      a ( 2, 5) =          0.33379999d0
      a ( 3, 5) =          1.39119995d0
      a ( 4, 5) =          1.49779999d0
      a ( 5, 5) =          3.37039995d0
      a ( 6, 5) =          3.54679990d0
      a ( 7, 5) =          8.28349972d0
      a ( 8, 5) =         23.96439934d0
      nelcor    = 68
c
      go to 30
c
c     gold (Au_CRENBS_ECP)
c
  79  continue
      cc( 1, 1) =         -0.98690200d0
      cc( 2, 1) =         -9.68630600d0
      cc( 3, 1) =        -32.57246399d0
      cc( 4, 1) =       -118.53087616d0
      cc( 5, 1) =       -286.25463867d0
      cc( 6, 1) =        -51.27423096d0
      a ( 1, 1) =          0.96569997d0
      a ( 2, 1) =          2.29040003d0
      a ( 3, 1) =          6.43839979d0
      a ( 4, 1) =         16.36459923d0
      a ( 5, 1) =         55.98929977d0
      a ( 6, 1) =        171.10139465d0
      cc( 1, 2) =        -84.96929932d0
      cc( 2, 2) =        280.17764282d0
      cc( 3, 2) =       -354.64355469d0
      cc( 4, 2) =        334.33190918d0
      cc( 5, 2) =       -219.39260864d0
      cc( 6, 2) =        170.85874939d0
      cc( 7, 2) =         56.33872986d0
      cc( 8, 2) =          6.40929985d0
      a ( 1, 2) =          1.22479999d0
      a ( 2, 2) =          1.43359995d0
      a ( 3, 2) =          1.89289999d0
      a ( 4, 2) =          2.69630003d0
      a ( 5, 2) =          3.98990011d0
      a ( 6, 2) =          6.00549984d0
      a ( 7, 2) =         14.58259964d0
      a ( 8, 2) =         39.31069946d0
      cc( 1, 3) =         -6.91190100d0
      cc( 2, 3) =         72.12546539d0
      cc( 3, 3) =       -114.97533417d0
      cc( 4, 3) =        123.61668396d0
      cc( 5, 3) =        -71.13571167d0
      cc( 6, 3) =         79.01660156d0
      cc( 7, 3) =         51.28700638d0
      cc( 8, 3) =          5.65198278d0
      a ( 1, 3) =          0.90640002d0
      a ( 2, 3) =          1.03110003d0
      a ( 3, 3) =          1.32869995d0
      a ( 4, 3) =          1.90369999d0
      a ( 5, 3) =          2.83640003d0
      a ( 6, 3) =          4.58839989d0
      a ( 7, 3) =         11.19810009d0
      a ( 8, 3) =         32.29570007d0
      cc( 1, 4) =         -0.11365700d0
      cc( 2, 4) =         -2.40402508d0
      cc( 3, 4) =         62.58111572d0
      cc( 4, 4) =        255.44717407d0
      cc( 5, 4) =       -161.80647278d0
      cc( 6, 4) =       -175.77076721d0
      cc( 7, 4) =         63.62079239d0
      cc( 8, 4) =          6.75031996d0
      a ( 1, 4) =          0.29310000d0
      a ( 2, 4) =          1.03009999d0
      a ( 3, 4) =          2.01859999d0
      a ( 4, 4) =          3.43919992d0
      a ( 5, 4) =          5.15789986d0
      a ( 6, 4) =          2.51489997d0
      a ( 7, 4) =          8.06239986d0
      a ( 8, 4) =         41.02849960d0
      cc( 1, 5) =        -47.60779953d0
      cc( 2, 5) =        176.31802368d0
      cc( 3, 5) =       -168.67689514d0
      cc( 4, 5) =         54.49244690d0
      cc( 5, 5) =         37.71539307d0
      cc( 6, 5) =        -24.13967514d0
      cc( 7, 5) =         49.16802216d0
      cc( 8, 5) =          4.75958586d0
      a ( 1, 5) =          0.43349999d0
      a ( 2, 5) =          0.50970000d0
      a ( 3, 5) =          0.57959998d0
      a ( 4, 5) =          0.77109998d0
      a ( 5, 5) =          2.59349990d0
      a ( 6, 5) =          8.39070034d0
      a ( 7, 5) =          6.92939997d0
      a ( 8, 5) =         22.25889969d0
      nelcor    = 68
      go to 30
c
c     mercury (Hg_CRENBS_ECP)
c
  80  continue
      cc( 1, 1) =         -1.12559295d0
      cc( 2, 1) =        -10.71220970d0
      cc( 3, 1) =        -37.99799728d0
      cc( 4, 1) =       -124.30549622d0
      cc( 5, 1) =       -312.36462402d0
      cc( 6, 1) =        -52.26305771d0
      a ( 1, 1) =          1.06610000d0
      a ( 2, 1) =          2.55520010d0
      a ( 3, 1) =          7.32770014d0
      a ( 4, 1) =         18.27160072d0
      a ( 5, 1) =         62.89289856d0
      a ( 6, 1) =        194.98759460d0
      cc( 1, 2) =        -84.91705322d0
      cc( 2, 2) =        276.98867798d0
      cc( 3, 2) =       -338.39907837d0
      cc( 4, 2) =        317.14492798d0
      cc( 5, 2) =       -203.94287109d0
      cc( 6, 2) =        165.25131226d0
      cc( 7, 2) =         50.15279770d0
      cc( 8, 2) =          3.97959089d0
      a ( 1, 2) =          1.31060004d0
      a ( 2, 2) =          1.52559996d0
      a ( 3, 2) =          1.99909997d0
      a ( 4, 2) =          2.82430005d0
      a ( 5, 2) =          4.17210007d0
      a ( 6, 2) =          6.28849983d0
      a ( 7, 2) =         15.61279964d0
      a ( 8, 2) =         40.71989822d0
      cc( 1, 3) =        -14.86320400d0
      cc( 2, 3) =         94.56565857d0
      cc( 3, 3) =       -141.04397583d0
      cc( 4, 3) =        143.83683777d0
      cc( 5, 3) =        -80.44802094d0
      cc( 6, 3) =         94.06087494d0
      cc( 7, 3) =         52.95009613d0
      cc( 8, 3) =          5.60394382d0
      a ( 1, 3) =          0.97350001d0
      a ( 2, 3) =          1.12080002d0
      a ( 3, 3) =          1.47000003d0
      a ( 4, 3) =          2.11899996d0
      a ( 5, 3) =          3.24379992d0
      a ( 6, 3) =          5.32800007d0
      a ( 7, 3) =         12.84560013d0
      a ( 8, 3) =         36.97430038d0
      cc( 1, 4) =         64.18244171d0
      cc( 2, 4) =         -0.14345001d0
      cc( 3, 4) =        -47.45405960d0
      cc( 4, 4) =       -103.45832825d0
      cc( 5, 4) =        235.98704529d0
      cc( 6, 4) =       -176.66653442d0
      cc( 7, 4) =         67.14499664d0
      cc( 8, 4) =          6.69803810d0
      a ( 1, 4) =          1.48909998d0
      a ( 2, 4) =          0.35940000d0
      a ( 3, 4) =          1.38450003d0
      a ( 4, 4) =          2.66899991d0
      a ( 5, 4) =          3.74539995d0
      a ( 6, 4) =          5.63350010d0
      a ( 7, 4) =          8.81350040d0
      a ( 8, 4) =         47.66930008d0
      cc( 1, 5) =        -41.44611359d0
      cc( 2, 5) =        136.94160461d0
      cc( 3, 5) =       -188.02015686d0
      cc( 4, 5) =        196.11465454d0
      cc( 5, 5) =       -135.42480469d0
      cc( 6, 5) =         97.13301849d0
      cc( 7, 5) =         44.36238861d0
      cc( 8, 5) =          5.04563379d0
      a ( 1, 5) =          0.72380000d0
      a ( 2, 5) =          0.82499999d0
      a ( 3, 5) =          1.04799998d0
      a ( 4, 5) =          1.42739999d0
      a ( 5, 5) =          2.01430011d0
      a ( 6, 5) =          2.81360006d0
      a ( 7, 5) =          8.31449985d0
      a ( 8, 5) =         22.51950073d0
      nelcor    = 68
c
      go to 30
c
c     lead  (Pb_CRENBS_ECP)
c
  82  continue
      cc( 1, 1) =         -0.74859297d0
      cc( 2, 1) =         -8.60412693d0
      cc( 3, 1) =        -26.17820168d0
      cc( 4, 1) =       -107.49287415d0
      cc( 5, 1) =       -246.23016357d0
      cc( 6, 1) =        -56.31761551d0
      a ( 1, 1) =          0.50800002d0
      a ( 2, 1) =          1.29289997d0
      a ( 3, 1) =          3.43869996d0
      a ( 4, 1) =         11.49810028d0
      a ( 5, 1) =         36.33620071d0
      a ( 6, 1) =        113.20490265d0
      cc( 1, 2) =         -0.14575300d0
      cc( 2, 2) =         38.23120117d0
      cc( 3, 2) =        -96.50691986d0
      cc( 4, 2) =        241.82226562d0
      cc( 5, 2) =       -282.88000488d0
      cc( 6, 2) =        164.92074585d0
      cc( 7, 2) =         52.07824707d0
      cc( 8, 2) =          6.50291586d0
      a ( 1, 2) =          0.22970000d0
      a ( 2, 2) =          0.93140000d0
      a ( 3, 2) =          1.07169998d0
      a ( 4, 2) =          1.55499995d0
      a ( 5, 2) =          2.09069991d0
      a ( 6, 2) =          2.90470004d0
      a ( 7, 2) =          9.83230019d0
      a ( 8, 2) =         27.40800095d0
      cc( 1, 3) =         -0.02626400d0
      cc( 2, 3) =        248.63769531d0
      cc( 3, 3) =        -63.24785995d0
      cc( 4, 3) =       -321.25946045d0
      cc( 5, 3) =        179.96823120d0
      cc( 6, 3) =          1.05947304d0
      cc( 7, 3) =         47.10823059d0
      cc( 8, 3) =          5.81061316d0
      a ( 1, 3) =          0.14990000d0
      a ( 2, 3) =          1.24570000d0
      a ( 3, 3) =          1.03120005d0
      a ( 4, 3) =          1.59749997d0
      a ( 5, 3) =          2.07649994d0
      a ( 6, 3) =          5.17929983d0
      a ( 7, 3) =          8.42840004d0
      a ( 8, 3) =         24.58130074d0
      cc( 1, 4) =         93.85963440d0
      cc( 2, 4) =       -202.12463379d0
      cc( 3, 4) =        226.66917419d0
      cc( 4, 4) =       -178.51530457d0
      cc( 5, 4) =         97.73090363d0
      cc( 6, 4) =         47.09598160d0
      cc( 7, 4) =         44.10550690d0
      cc( 8, 4) =          7.95272112d0
      a ( 1, 4) =          0.61720002d0
      a ( 2, 4) =          0.70740002d0
      a ( 3, 4) =          0.91750002d0
      a ( 4, 4) =          1.22179997d0
      a ( 5, 4) =          1.64769995d0
      a ( 6, 4) =          5.08220005d0
      a ( 7, 4) =          8.89140034d0
      a ( 8, 4) =         17.54299927d0
      cc( 1, 5) =        124.17366028d0
      cc( 2, 5) =       -255.84779358d0
      cc( 3, 5) =        257.96179199d0
      cc( 4, 5) =       -191.85612488d0
      cc( 5, 5) =        113.37288666d0
      cc( 6, 5) =         90.45600128d0
      cc( 7, 5) =         51.46012878d0
      cc( 8, 5) =          4.82522202d0
      a ( 1, 5) =          0.81099999d0
      a ( 2, 5) =          0.94230002d0
      a ( 3, 5) =          1.22370005d0
      a ( 4, 5) =          1.67830002d0
      a ( 5, 5) =          2.25040007d0
      a ( 6, 5) =          8.25269985d0
      a ( 7, 5) =         16.56539917d0
      a ( 8, 5) =         42.04230118d0
      nelcor    = 78
c
      go to 30
c
c     bismuth (Bi_CRENBS_ECP)
c
   83 continue
      cc( 1, 1) =         -0.70747900d0
      cc( 2, 1) =         -8.35064983d0
      cc( 3, 1) =        -25.68899536d0
      cc( 4, 1) =        -99.04219818d0
      cc( 5, 1) =       -227.29902649d0
      cc( 6, 1) =        -55.69332886d0
      a ( 1, 1) =          0.56889999d0
      a ( 2, 1) =          1.39390004d0
      a ( 3, 1) =          3.53730011d0
      a ( 4, 1) =         11.37679958d0
      a ( 5, 1) =         33.82360077d0
      a ( 6, 1) =        104.52100372d0
      cc( 1, 2) =        -34.60086823d0
      cc( 2, 2) =         38.51860809d0
      cc( 3, 2) =        -48.42468262d0
      cc( 4, 2) =        150.18840027d0
      cc( 5, 2) =       -223.60801697d0
      cc( 6, 2) =        192.42153931d0
      cc( 7, 2) =         54.63947296d0
      cc( 8, 2) =          6.56644678d0
      a ( 1, 2) =          0.52440000d0
      a ( 2, 2) =          0.54009998d0
      a ( 3, 2) =          1.12189996d0
      a ( 4, 2) =          1.53910005d0
      a ( 5, 2) =          2.36249995d0
      a ( 6, 2) =          3.00480008d0
      a ( 7, 2) =          9.95259953d0
      a ( 8, 2) =         27.90049934d0
      cc( 1, 3) =         -0.02339300d0
      cc( 2, 3) =        -74.60871887d0
      cc( 3, 3) =        231.86610413d0
      cc( 4, 3) =       -272.94262695d0
      cc( 5, 3) =        163.37461853d0
      cc( 6, 3) =         -7.98877001d0
      cc( 7, 3) =         45.16691589d0
      cc( 8, 3) =          5.93953514d0
      a ( 1, 3) =          0.16890000d0
      a ( 2, 3) =          1.08780003d0
      a ( 3, 3) =          1.27289999d0
      a ( 4, 3) =          1.64979994d0
      a ( 5, 3) =          2.12190008d0
      a ( 6, 3) =          9.05420017d0
      a ( 7, 3) =          7.60179996d0
      a ( 8, 3) =         18.94759941d0
      cc( 1, 4) =        106.07496643d0
      cc( 2, 4) =       -232.36282349d0
      cc( 3, 4) =        250.49478149d0
      cc( 4, 4) =       -186.17285156d0
      cc( 5, 4) =        113.54383087d0
      cc( 6, 4) =         38.24295044d0
      cc( 7, 4) =         43.72408295d0
      cc( 8, 4) =          8.00164509d0
      a ( 1, 4) =          0.74769998d0
      a ( 2, 4) =          0.87339997d0
      a ( 3, 4) =          1.14359999d0
      a ( 4, 4) =          1.58449996d0
      a ( 5, 4) =          2.19569993d0
      a ( 6, 4) =          6.44960022d0
      a ( 7, 4) =          9.12709999d0
      a ( 8, 4) =         18.45389938d0
      cc( 1, 5) =        110.70706177d0
      cc( 2, 5) =       -148.41661072d0
      cc( 3, 5) =        123.52549744d0
      cc( 4, 5) =       -189.74470520d0
      cc( 5, 5) =        140.47151184d0
      cc( 6, 5) =         52.21812058d0
      cc( 7, 5) =         37.94526672d0
      cc( 8, 5) =          5.43266582d0
      a ( 1, 5) =          0.57910001d0
      a ( 2, 5) =          0.63279998d0
      a ( 3, 5) =          0.98130000d0
      a ( 4, 5) =          1.35769999d0
      a ( 5, 5) =          1.64100003d0
      a ( 6, 5) =          5.68540001d0
      a ( 7, 5) =         10.56250000d0
      a ( 8, 5) =         17.24570084d0
      nelcor    = 78
c
      go to 30
c
c     Po (Po_CRENBS_ECP)
c
   84 continue
      cc( 1, 1) =         -0.68128300d0
      cc( 2, 1) =         -8.54319859d0
      cc( 3, 1) =        -26.42345810d0
      cc( 4, 1) =        -99.83397675d0
      cc( 5, 1) =       -228.19705200d0
      cc( 6, 1) =        -55.57450104d0
      a ( 1, 1) =          0.62919998d0
      a ( 2, 1) =          1.52219999d0
      a ( 3, 1) =          3.77830005d0
      a ( 4, 1) =         11.88829994d0
      a ( 5, 1) =         34.90750122d0
      a ( 6, 1) =        106.91159821d0
      cc( 1, 2) =         -0.08444800d0
      cc( 2, 2) =        -80.85897827d0
      cc( 3, 2) =        254.40237427d0
      cc( 4, 2) =       -269.28668213d0
      cc( 5, 2) =        235.26022339d0
      cc( 6, 2) =       -115.93589020d0
      cc( 7, 2) =         65.20384979d0
      cc( 8, 2) =          6.07269716d0
      a ( 1, 2) =          0.25819999d0
      a ( 2, 2) =          1.56389999d0
      a ( 3, 2) =          1.88989997d0
      a ( 4, 2) =          2.68810010d0
      a ( 5, 2) =          4.13819981d0
      a ( 6, 2) =          6.86339998d0
      a ( 7, 2) =          9.40100002d0
      a ( 8, 2) =         36.99940109d0
      cc( 1, 3) =         -0.02554400d0
      cc( 2, 3) =        -71.26758575d0
      cc( 3, 3) =        220.39631653d0
      cc( 4, 3) =       -248.64100647d0
      cc( 5, 3) =        156.54104614d0
      cc( 6, 3) =        140.72520447d0
      cc( 7, 3) =         66.54429626d0
      cc( 8, 3) =          5.46325779d0
      a ( 1, 3) =          0.20640001d0
      a ( 2, 3) =          1.15900004d0
      a ( 3, 3) =          1.36510003d0
      a ( 4, 3) =          1.80060005d0
      a ( 5, 3) =          2.40829992d0
      a ( 6, 3) =         10.21560001d0
      a ( 7, 3) =         23.22529984d0
      a ( 8, 3) =         66.82530212d0
      cc( 1, 4) =        -52.11019135d0
      cc( 2, 4) =        168.50531006d0
      cc( 3, 4) =       -243.10089111d0
      cc( 4, 4) =        250.63352966d0
      cc( 5, 4) =       -185.83685303d0
      cc( 6, 4) =        120.02333069d0
      cc( 7, 4) =         45.58230591d0
      cc( 8, 4) =          7.91772079d0
      a ( 1, 4) =          0.65869999d0
      a ( 2, 4) =          0.74949998d0
      a ( 3, 4) =          0.94880003d0
      a ( 4, 4) =          1.28760004d0
      a ( 5, 4) =          1.81920004d0
      a ( 6, 4) =          2.51469994d0
      a ( 7, 4) =          7.50089979d0
      a ( 8, 4) =         16.72979927d0
      cc( 1, 5) =         -0.15798301d0
      cc( 2, 5) =        -49.45363617d0
      cc( 3, 5) =        159.05288696d0
      cc( 4, 5) =       -265.50646973d0
      cc( 5, 5) =        231.14622498d0
      cc( 6, 5) =        -54.27425766d0
      cc( 7, 5) =         53.80932999d0
      cc( 8, 5) =          4.43020916d0
      a ( 1, 5) =          0.42980000d0
      a ( 2, 5) =          2.06590009d0
      a ( 3, 5) =          2.49589992d0
      a ( 4, 5) =          3.37210011d0
      a ( 5, 5) =          4.42969990d0
      a ( 6, 5) =          9.18620014d0
      a ( 7, 5) =         12.39190006d0
      a ( 8, 5) =         49.19279861d0
      nelcor    = 78
c
      go to 30
c
c     at (At_CRENBS_ECP)
c
   85 continue
      cc( 1, 1) =         -0.77448499d0
      cc( 2, 1) =         -9.04102802d0
      cc( 3, 1) =        -28.50047112d0
      cc( 4, 1) =       -109.59169006d0
      cc( 5, 1) =       -245.29591370d0
      cc( 6, 1) =        -56.34064484d0
      a ( 1, 1) =          0.71289998d0
      a ( 2, 1) =          1.67809999d0
      a ( 3, 1) =          4.16319990d0
      a ( 4, 1) =         13.12460041d0
      a ( 5, 1) =         39.27880096d0
      a ( 6, 1) =        119.85289764d0
      cc( 1, 2) =        -45.67423248d0
      cc( 2, 2) =        132.16868591d0
      cc( 3, 2) =       -217.87092590d0
      cc( 4, 2) =        325.10427856d0
      cc( 5, 2) =       -306.25015259d0
      cc( 6, 2) =        195.50265503d0
      cc( 7, 2) =         40.29927444d0
      cc( 8, 2) =          7.16939592d0
      a ( 1, 2) =          0.83800000d0
      a ( 2, 2) =          0.96259999d0
      a ( 3, 2) =          1.22790003d0
      a ( 4, 2) =          1.69780004d0
      a ( 5, 2) =          2.43869996d0
      a ( 6, 2) =          3.35560012d0
      a ( 7, 2) =          9.58790016d0
      a ( 8, 2) =         15.94320011d0
      cc( 1, 3) =         -0.02650500d0
      cc( 2, 3) =        -90.51212311d0
      cc( 3, 3) =        286.96719360d0
      cc( 4, 3) =       -331.29653931d0
      cc( 5, 3) =        193.59178162d0
      cc( 6, 3) =         -9.22439098d0
      cc( 7, 3) =         47.91704941d0
      cc( 8, 3) =          5.86008310d0
      a ( 1, 3) =          0.23510000d0
      a ( 2, 3) =          1.32900000d0
      a ( 3, 3) =          1.55100000d0
      a ( 4, 3) =          2.01480007d0
      a ( 5, 3) =          2.61159992d0
      a ( 6, 3) =          7.40170002d0
      a ( 7, 3) =          8.69470024d0
      a ( 8, 3) =         22.87330055d0
      cc( 1, 4) =        -48.58103561d0
      cc( 2, 4) =        140.93241882d0
      cc( 3, 4) =       -197.97145081d0
      cc( 4, 4) =        188.84280396d0
      cc( 5, 4) =       -125.92828369d0
      cc( 6, 4) =        112.62237549d0
      cc( 7, 4) =         56.06732941d0
      cc( 8, 4) =          7.51697207d0
      a ( 1, 4) =          0.71300000d0
      a ( 2, 4) =          0.79850000d0
      a ( 3, 4) =          1.02859998d0
      a ( 4, 4) =          1.37629998d0
      a ( 5, 4) =          2.12829995d0
      a ( 6, 4) =          3.01410007d0
      a ( 7, 4) =          9.30930042d0
      a ( 8, 4) =         25.96680069d0
      cc( 1, 5) =         31.65410233d0
      cc( 2, 5) =        -52.62594604d0
      cc( 3, 5) =         48.30789185d0
      cc( 4, 5) =        -20.87113571d0
      cc( 5, 5) =         42.20579529d0
      cc( 6, 5) =         17.82291985d0
      cc( 7, 5) =         47.56859589d0
      cc( 8, 5) =          5.06970882d0
      a ( 1, 5) =          0.51190001d0
      a ( 2, 5) =          0.63190001d0
      a ( 3, 5) =          0.91280001d0
      a ( 4, 5) =          1.45959997d0
      a ( 5, 5) =          2.70810008d0
      a ( 6, 5) =          4.43590021d0
      a ( 7, 5) =          9.46450043d0
      a ( 8, 5) =         23.41760063d0
      nelcor    = 78
      go to 30
c
c     radon (Rn_CRENBS_ECP)
c
   86 continue
c
      cc( 1, 1) =         -0.72139901d0
      cc( 2, 1) =         -8.75459194d0
      cc( 3, 1) =        -28.94036484d0
      cc( 4, 1) =       -107.48715973d0
      cc( 5, 1) =       -240.13769531d0
      cc( 6, 1) =        -56.16814423d0
      a ( 1, 1) =          0.76840001d0
      a ( 2, 1) =          1.77160001d0
      a ( 3, 1) =          4.31479979d0
      a ( 4, 1) =         13.43430042d0
      a ( 5, 1) =         39.14410019d0
      a ( 6, 1) =        118.78569794d0
      cc( 1, 2) =        -41.80408478d0
      cc( 2, 2) =        121.57847595d0
      cc( 3, 2) =       -207.84539795d0
      cc( 4, 2) =        339.97991943d0
      cc( 5, 2) =       -331.86575317d0
      cc( 6, 2) =        214.76388550d0
      cc( 7, 2) =         63.20500565d0
      cc( 8, 2) =          6.29603815d0
      a ( 1, 2) =          0.92170000d0
      a ( 2, 2) =          1.05760002d0
      a ( 3, 2) =          1.35520005d0
      a ( 4, 2) =          1.90330005d0
      a ( 5, 2) =          2.81410003d0
      a ( 6, 2) =          4.04860020d0
      a ( 7, 2) =         14.25360012d0
      a ( 8, 2) =         47.45729828d0
      cc( 1, 3) =        -43.26960754d0
      cc( 2, 3) =        179.17414856d0
      cc( 3, 3) =       -231.13307190d0
      cc( 4, 3) =        229.06947327d0
      cc( 5, 3) =       -126.94128418d0
      cc( 6, 3) =        152.43060303d0
      cc( 7, 3) =         58.18738174d0
      cc( 8, 3) =          5.53178597d0
      a ( 1, 3) =          1.26999998d0
      a ( 2, 3) =          1.60210001d0
      a ( 3, 3) =          2.24740005d0
      a ( 4, 3) =          3.36470008d0
      a ( 5, 3) =          5.19239998d0
      a ( 6, 3) =          8.58419991d0
      a ( 7, 3) =         19.64229965d0
      a ( 8, 3) =         52.35680008d0
      cc( 1, 4) =        -41.88452148d0
      cc( 2, 4) =        142.92504883d0
      cc( 3, 4) =       -203.37625122d0
      cc( 4, 4) =        205.23439026d0
      cc( 5, 4) =       -142.26692200d0
      cc( 6, 4) =        118.22145844d0
      cc( 7, 4) =         58.02330017d0
      cc( 8, 4) =          7.45458698d0
      a ( 1, 4) =          0.79329997d0
      a ( 2, 4) =          0.90460002d0
      a ( 3, 4) =          1.15869999d0
      a ( 4, 4) =          1.60080004d0
      a ( 5, 4) =          2.36409998d0
      a ( 6, 4) =          3.43790007d0
      a ( 7, 4) =         10.24720001d0
      a ( 8, 4) =         28.55220032d0
      cc( 1, 5) =          2.04597211d0
      cc( 2, 5) =         10.23069572d0
      cc( 3, 5) =         26.78542519d0
      cc( 4, 5) =         41.19087219d0
      cc( 5, 5) =        -31.82834816d0
      cc( 6, 5) =         78.47965240d0
      cc( 7, 5) =         61.61168289d0
      cc( 8, 5) =          9.20025444d0
      a ( 1, 5) =          0.14659999d0
      a ( 2, 5) =          0.51359999d0
      a ( 3, 5) =          1.44169998d0
      a ( 4, 5) =          2.17149997d0
      a ( 5, 5) =          1.90530002d0
      a ( 6, 5) =          4.78660011d0
      a ( 7, 5) =         10.76189995d0
      a ( 8, 5) =         22.37039948d0
      nelcor    = 78
c
 30   continue
      n ( 1, 1) =   2
      n ( 2, 1) =   2
      n ( 3, 1) =   2
      n ( 4, 1) =   2
      n ( 5, 1) =   2
      n ( 6, 1) =   1
      nt( 1)    =   6
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      n ( 4, 2) =   2
      n ( 5, 2) =   2
      n ( 6, 2) =   2
      n ( 7, 2) =   1
      n ( 8, 2) =   0
      nt( 2)    =   8
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      n ( 3, 3) =   2
      n ( 4, 3) =   2
      n ( 5, 3) =   2
      n ( 6, 3) =   2
      n ( 7, 3) =   1
      n ( 8, 3) =   0
      nt( 3)    =   8
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      n ( 3, 4) =   2
      n ( 4, 4) =   2
      n ( 5, 4) =   2
      n ( 6, 4) =   2
      n ( 7, 4) =   1
      n ( 8, 4) =   0
      nt( 4)    =   8
      n ( 1, 5) =   2
      n ( 2, 5) =   2
      n ( 3, 5) =   2
      n ( 4, 5) =   2
      n ( 5, 5) =   2
      n ( 6, 5) =   2
      n ( 7, 5) =   1
      n ( 8, 5) =   0
      nt( 5)    =   8
      lmx       =   4
      return
 6010 format (1x,'*** data not found for CRENBS ecp type ',a4)
      end
**==strlc0.f
      subroutine strlc0(zinp,lmx,nelcor,nt,n,cc,a)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character * 8 (z)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      dimension nt(6), n(10,6), cc(10,6), a(10,6)
c
      dimension ztit(18)
      data maxtyp/18/
      data ztit/'h','he',
     *          'li','be','b ','c ','n ','o ','f ','ne',
     *          'na','mg','al','si','p ','s ','cl','ar'/
c
c     this routine contains a library of ecp data taken from
c     various literature sources.
c     Stuttgart Relativistic, Large Core ECP
c
      do ityp = 1 , maxtyp
         if (zinp.eq.ztit(ityp)) go to 40
      enddo
 30   write (iwr,6010) zinp
      call caserr2('no library available for requested STRLC ecp')
c
 40   go to (30,30,50,60,70,80,90,100,110,120,
     +      130,150,170,190,210,230,250,270) , ityp
      go to 30
c
c     lithium   (Li_Stuttgart_RLC_ECP)
c
 50   continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =          5.78600000d0
      a ( 1, 2) =          1.27600000d0
      cc( 1, 3) =         -1.06500000d0
      a ( 1, 3) =          1.60700000d0
c
      go to 65
c
c     beryllium   (Be_Stuttgart_RLC_ECP)
c
 60   continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =         13.32500000d0
      a ( 1, 2) =          2.65300000d0
      cc( 1, 3) =         -1.57400000d0
      a ( 1, 3) =          3.12000000d0
c
 65   continue
      n ( 1, 1) =   2
      nt( 1)    =   1
      n ( 1, 2) =   2
      nt( 2)    =   1
      n ( 1, 3) =   2
      nt( 3)    =   1
      nelcor    =   2
      lmx       =   2
c
      return
c
c     boron   (B_Stuttgart_RLC_ECP)
c
 70   continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =         23.99296000d0
      a ( 1, 2) =          4.50610000d0
      cc( 1, 3) =         -1.76672400d0
      a ( 1, 3) =          5.22950000d0
      cc( 1, 4) =         -0.00175400d0
      a ( 1, 4) =          0.01080000d0
c
      go to 105
c
c     carbon  (C_Stuttgart_RLC_ECP)
c
 80   continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =         33.12163800d0
      a ( 1, 2) =          6.40105200d0
      cc( 1, 3) =         -1.98625700d0
      a ( 1, 3) =          7.30774700d0
      cc( 1, 4) =         -9.45431800d0
      a ( 1, 4) =          5.96179600d0
c
      go to 105
c
c     nitrogen  (N_Stuttgart_RLC_ECP)
c
 90   continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =         38.53383100d0
      a ( 1, 2) =          7.97723200d0
      cc( 1, 3) =         -2.55081000d0
      a ( 1, 3) =         10.18385400d0
      cc( 1, 4) =         -2.99554500d0
      a ( 1, 4) =         11.55994700d0
c
      go to 105
c
c     oxygen  (O_Stuttgart_RLC_ECP)
c
 100  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =         50.77106900d0
      a ( 1, 2) =         10.44567000d0
      cc( 1, 3) =         -4.90355100d0
      a ( 1, 3) =         18.04517400d0
      cc( 1, 4) =         -3.31212400d0
      a ( 1, 4) =          8.16479800d0
c
 105  continue
      n ( 1, 1) =   2
      nt( 1)    =   1
      n ( 1, 2) =   2
      nt( 2)    =   1
      n ( 1, 3) =   2
      nt( 3)    =   1
      n ( 1, 4) =   2
      nt( 4)    =   1
      nelcor    =   2
      lmx       =   3
c
      return
c
c     fluorine  (F_Stuttgart_RLC_ECP)
c
 110  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        102.59795200d0
      cc( 2, 2) =         19.04966300d0
      a ( 1, 2) =         22.35040000d0
      a ( 2, 2) =         11.17520000d0
      cc( 1, 3) =        -15.14396000d0
      cc( 2, 3) =          2.80292100d0
      a ( 1, 3) =         26.47680000d0
      a ( 2, 3) =         13.23840000d0
      cc( 1, 4) =         -0.00186300d0
      a ( 1, 4) =          0.03160000d0
c
      go to 125
c
c     neon  (Ne_Stuttgart_RLC_ECP)
c
 120  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        112.52543566d0
      cc( 2, 2) =         28.30083454d0
      a ( 1, 2) =         31.86016200d0
      a ( 2, 2) =         12.36221900d0
      cc( 1, 3) =        -11.12658543d0
      cc( 2, 3) =          3.38754919d0
      a ( 1, 3) =         21.50803400d0
      a ( 2, 3) =         12.91044700d0
      cc( 1, 4) =         -0.18408921d0
      a ( 1, 4) =          0.85038500d0
c
 125  continue
      n ( 1, 1) =   2
      nt( 1)    =   1
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      nt( 2)    =   2
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      nt( 3)    =   2
      n ( 1, 4) =   2
      nt( 4)    =   1
      nelcor    =   2
      lmx       =   3
c
      return
c
c      (Na_Stuttgart_RLC_ECP)
c
 130  continue
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =         10.83900000d0
      a ( 1, 2) =          1.37800000d0
      cc( 1, 3) =          2.30300000d0
      a ( 1, 3) =          0.66390000d0
      cc( 1, 4) =         -1.77700000d0
      a ( 1, 4) =          0.92490000d0
c
      go to 235
c
c      (Mg_Stuttgart_RLC_ECP)
c
 150  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =         14.67600000d0
      a ( 1, 2) =          1.73200000d0
      cc( 1, 3) =          5.17570000d0
      a ( 1, 3) =          1.11500000d0
      cc( 1, 4) =         -1.81600000d0
      a ( 1, 4) =          1.20300000d0
c
      go to 235
c
c      (Al_Stuttgart_RLC_ECP)
c
 170  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =         20.40981300d0
      a ( 1, 2) =          2.19822500d0
      cc( 1, 3) =          8.98049500d0
      a ( 1, 3) =          1.60139500d0
      cc( 1, 4) =         -1.97041100d0
      a ( 1, 4) =          1.49902600d0
c
      go to 235
c
c      (Si_Stuttgart_RLC_ECP)
c
 190  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =         26.62331865d0
      a ( 1, 2) =          2.71362197d0
      cc( 1, 3) =         10.92995391d0
      a ( 1, 3) =          1.96687987d0
      cc( 1, 4) =         -4.66941200d0
      a ( 1, 4) =          2.71001600d0
c
      go to 235
c
c       (P_Stuttgart_RLC_ECP)
c
 210  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =         26.53226131d0
      a ( 1, 2) =          2.94055986d0
      cc( 1, 3) =         11.49721021d0
      a ( 1, 3) =          2.22771255d0
      cc( 1, 4) =        -16.77278000d0
      a ( 1, 4) =          5.66170600d0
c
      go to 235
c
c      (S_Stuttgart_RLC_ECP)
c
 230  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =         37.97481900d0
      a ( 1, 2) =          3.74389164d0
      cc( 1, 3) =         18.79052931d0
      a ( 1, 3) =          3.08608744d0
      cc( 1, 4) =         -7.83796400d0
      a ( 1, 4) =          4.86241400d0
c
235   continue
      n ( 1, 1) =   2
      nt( 1)    =   1
      n ( 1, 2) =   2
      nt( 2)    =   1
      n ( 1, 3) =   2
      nt( 3)    =   1
      n ( 1, 4) =   2
      nt( 4)    =   1
      nelcor    =  10
      lmx       =   3
c
      return
c
c      (Cl_Stuttgart_RLC_ECP)
c
 250  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =         33.13663200d0
      cc( 2, 2) =         16.27072800d0
      a ( 1, 2) =          6.39430000d0
      a ( 2, 2) =          3.19710000d0
      cc( 1, 3) =         24.41699300d0
      cc( 2, 3) =          7.68305000d0
      a ( 1, 3) =          5.62070000d0
      a ( 2, 3) =          2.81030000d0
      cc( 1, 4) =         -8.58764900d0
      a ( 1, 4) =          5.33810000d0
c
      n ( 1, 1) =   2
      nt( 1)    =   1
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      nt( 2)    =   2
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      nt( 3)    =   2
      n ( 1, 4) =   2
      nt( 4)    =   1
      nelcor    =  10
      lmx       =   3
c
      return
c
c      (Ar_Stuttgart_RLC_ECP)
c
 270  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =         68.66778801d0
      cc( 2, 2) =         24.04276629d0
      a ( 1, 2) =         10.26172100d0
      a ( 2, 2) =          3.95272500d0
      cc( 1, 3) =         27.73076331d0
      cc( 2, 3) =          4.04545904d0
      a ( 1, 3) =          5.39271400d0
      a ( 2, 3) =          2.69996700d0
      cc( 1, 4) =         -8.13747696d0
      cc( 2, 4) =         -1.66452808d0
      a ( 1, 4) =          8.08623500d0
      a ( 2, 4) =          4.01663200d0
      cc( 1, 5) =         -3.40009845d0
      a ( 1, 5) =          5.20845900d0
c
      n ( 1, 1) =   2
      nt( 1)    =   1
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      nt( 2)    =   2
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      nt( 3)    =   2
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      nt( 4)    =   2
      n ( 1, 5) =   2
      nt( 5)    =   1
      nelcor    =  10
      lmx       =   4
c
      return
 6010 format (1x,'*** data not found for STRLC ecp type ',a7)
      end
**==strlc1.f
      subroutine strlc1(zinp,lmx,nelcor,nt,n,cc,a)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character * 8 (z)
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      dimension nt(6), n(10,6), cc(10,6), a(10,6)
c
      dimension ztit(18)
      data maxtyp/18/
      data ztit/'k','ca',
     *          'sc','ti','v ','cr','mn','fe',
     *          'co','ni','cu','zn','ga','ge','as','se','br','kr'/
c
c     this routine contains a library of ecp data taken from
c     various literature sources.
c     Stuttgart Relativistic, Large Core ECP
c
      do 20 ityp = 1 , maxtyp
         if (zinp.eq.ztit(ityp)) go to 30
 20   continue
 60   write (iwr,6010) zinp
      call caserr2('no library available for requested STRLC ecp')
 30   go to (40,50,60,60,60,60,60,60,60,60,60,150,160,170,180,190,
     +       200,210) , ityp
c
c     k   (K_Stuttgart_RLC_ECP)
c
 40   continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =         13.56400000d0
      a ( 1, 2) =          0.85300000d0
      cc( 1, 3) =          2.64800000d0
      a ( 1, 3) =          0.36960000d0
      cc( 1, 4) =         -4.51700000d0
      a ( 1, 4) =          0.66390000d0
c
      go to 55
c
c     ca  (Ca_Stuttgart_RLC_ECP)
c
 50   continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =         12.46600000d0
      a ( 1, 2) =          0.89800000d0
      cc( 1, 3) =          5.14600000d0
      a ( 1, 3) =          0.54800000d0
      cc( 1, 4) =         -7.70900000d0
      a ( 1, 4) =          1.11900000d0

 55    continue
      n ( 1, 1) =   2
      nt( 1)    =   1
      n ( 1, 2) =   2
      nt( 2)    =   1
      n ( 1, 3) =   2
      nt( 3)    =   1
      n ( 1, 4) =   2
      nt( 4)    =   1
      nelcor    =  18
      lmx       =   3
c
      return
c
c     zinc   (Zn_Stuttgart_RLC_ECP)
c
 150  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =         14.93680000d0
      cc( 2, 2) =         -2.62660000d0
      a ( 1, 2) =          1.50000000d0
      a ( 2, 2) =          0.75000000d0
      cc( 1, 3) =         11.43290000d0
      cc( 2, 3) =         -1.49390000d0
      a ( 1, 3) =          1.50000000d0
      a ( 2, 3) =          0.75000000d0
      cc( 1, 4) =          1.16900000d0
      cc( 2, 4) =          0.37070000d0
      a ( 1, 4) =          0.75000000d0
      a ( 2, 4) =          0.37500000d0
c 
      n ( 1, 1) =   2
      nt( 1)    =   1
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      nt( 2)    =   2
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      nt( 3)    =   2
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      nt( 4)    =   2
      nelcor    =  28
      lmx       =   3
c
c
      return
c
c     ga    (Ga_Stuttgart_RLC_ECP)
c
 160  continue
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        203.85397200d0
      a ( 1, 2) =          5.21596000d0
      cc( 1, 3) =        156.10339000d0
      a ( 1, 3) =          4.30890400d0
      cc( 1, 4) =          1.03164700d0
      a ( 1, 4) =          0.49635700d0
      cc( 1, 5) =        -10.67373500d0
      a ( 1, 5) =          1.71517000d0
c
      go to 195
c
c     ge   (Ge_Stuttgart_RLC_ECP)
c
 170  continue
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        149.24657900d0
      a ( 1, 2) =          4.81540900d0
      cc( 1, 3) =        132.84433500d0
      a ( 1, 3) =          4.16951500d0
      cc( 1, 4) =          1.34615400d0
      a ( 1, 4) =          0.59195800d0
      cc( 1, 5) =         -7.04422300d0
      a ( 1, 5) =          1.79177000d0
c
      go to 195
c
c     as   (As_Stuttgart_RLC_ECP)
c
 180  continue
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =         53.96562000d0
      a ( 1, 2) =          3.61262500d0
      cc( 1, 3) =         88.94908800d0
      a ( 1, 3) =          3.90792600d0
      cc( 1, 4) =         22.42028800d0
      a ( 1, 4) =          1.92646700d0
      cc( 1, 5) =         -4.70481500d0
      a ( 1, 5) =          1.77343400d0
c
      go to 195
c
c     se   (Se_Stuttgart_RLC_ECP)
c
 190  continue
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =         79.66334500d0
      a ( 1, 2) =          4.23705700d0
      cc( 1, 3) =         31.56099300d0
      a ( 1, 3) =          2.91033400d0
      cc( 1, 4) =         30.80461000d0
      a ( 1, 4) =          2.33570100d0
      cc( 1, 5) =         -6.54687500d0
      a ( 1, 5) =          2.25463900d0
c
  195 continue
      n ( 1, 1) =   2
      nt( 1)    =   1
      n ( 1, 2) =   2
      nt( 2)    =   1
      n ( 1, 3) =   2
      nt( 3)    =   1
      n ( 1, 4) =   2
      nt( 4)    =   1
      n ( 1, 5) =   2
      nt( 5)    =   1
      nelcor    =  28
      lmx       =   4
c
      return
c
c     br   (Br_Stuttgart_RLC_ECP)
c
 200  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =         61.51372100d0
      cc( 2, 2) =          9.02149300d0
      a ( 1, 2) =          5.02180000d0
      a ( 2, 2) =          2.51090000d0
      cc( 1, 3) =         53.87586400d0
      cc( 2, 3) =          4.62940200d0
      a ( 1, 3) =          4.28140000d0
      a ( 2, 3) =          2.14070000d0
      cc( 1, 4) =         20.84967700d0
      cc( 2, 4) =          2.96544400d0
      a ( 1, 4) =          2.88000000d0
      a ( 2, 4) =          1.44000000d0
      cc( 1, 5) =         -8.16149300d0
      a ( 1, 5) =          2.72070000d0
c
      n ( 1, 1) =   2
      nt( 1)    =   1
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      nt( 2)    =   2
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      nt( 3)    =   2
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      nt( 4)    =   2
      n ( 1, 5) =   2
      nt( 5)    =   1
      nelcor    =  28
      lmx       =   4
c
      return
c
c     kr   (Kr_Stuttgart_RLC_ECP)
c
 210  continue
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =         73.91569390d0
      cc( 2, 2) =         16.16825080d0
      a ( 1, 2) =          5.87771800d0
      a ( 2, 2) =          3.08462200d0
      cc( 1, 3) =         58.51769101d0
      cc( 2, 3) =          8.25910073d0
      a ( 1, 3) =          5.16411000d0
      a ( 2, 3) =          2.35830200d0
      cc( 1, 4) =         33.45822776d0
      cc( 2, 4) =          0.67725331d0
      a ( 1, 4) =          3.21536200d0
      a ( 2, 4) =          1.28500800d0
      cc( 1, 5) =        -15.15869859d0
      cc( 2, 5) =         -0.17408825d0
      a ( 1, 5) =          4.08286900d0
      a ( 2, 5) =          1.19396000d0
c     remove g potentials from expansion
c     cc( 1, 6) =         -6.83315877d0
c     a ( 1, 6) =          3.18077500d0
      n ( 1, 1) =   2
      nt( 1)    =   1
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      nt( 2)    =   2
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      nt( 3)    =   2
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      nt( 4)    =   2
      n ( 1, 5) =   2
      n ( 2, 5) =   2
      nt( 5)    =   2
c     n ( 1, 6) =   2
c     nt( 6)    =   1
      nelcor    =  28
c     lmx       =   5
      lmx       =   4
c
      return
 6010 format (1x,'*** data not found for ecp type ',a4)
      end
**==strlc2.f
      subroutine strlc2 (zinp,lmx,nelcor,nt,n,cc,a)
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
      implicit character * 8 (z)
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
      dimension nt(6), n(10,6), cc(10,6), a(10,6)
c
      dimension ztit(18)
      data maxtyp/18/
      data ztit/
     $ 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag','cd',
     $ 'in','sn','sb','te','i ','xe' /
c
c     this routine contains a library of ecp data taken from
c     various literature sources.
c     Stuttgart Relativistic, Large Core ECP
c
      do 10 ityp = 1 , maxtyp
         if (zinp.eq.ztit(ityp)) go to 20
 10   continue
 39   write (iwr,6010) zinp
      call caserr2('no library available for requested ecp')
 20   continue
c
      go to (37,38,39,39,39,39,39,39,39,39,39,39,49,
     +       50,51,52,53,54), ityp
c
c      rb (RB_Stuttgart_RLC_ECP)
c
  37  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      n ( 1, 1) =   2
      nt( 1)    =   1
      cc( 1, 2) =         45.27200000d0
      a ( 1, 2) =          1.01200000d0
      n ( 1, 2) =   2
      nt( 2)    =   1
      cc( 1, 3) =          2.83000000d0
      a ( 1, 3) =          0.30360000d0
      n ( 1, 3) =   2
      nt( 3)    =   1
      cc( 1, 4) =         -1.97500000d0
      a ( 1, 4) =          0.38380000d0
      n ( 1, 4) =   2
      nt( 4)    =   1
      nelcor    =  36
      lmx       =   3
c
      return
c
c      sr (SR_Stuttgart_RLC_ECP)
c
  38  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      n ( 1, 1) =   2
      nt( 1)    =   1
      cc( 1, 2) =         15.38700000d0
      a ( 1, 2) =          0.79620000d0
      n ( 1, 2) =   2
      nt( 2)    =   1
      cc( 1, 3) =          5.07700000d0
      a ( 1, 3) =          0.42120000d0
      n ( 1, 3) =   2
      nt( 3)    =   1
      cc( 1, 4) =         -2.24800000d0
      a ( 1, 4) =          0.49270000d0
      n ( 1, 4) =   2
      nt( 4)    =   1
      nelcor    =  36
      lmx       =   3
c
      return
c
c     indium (In_Stuttgart_RLC_ECP)
c
  49  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      n ( 1, 1) =   2
      nt( 1)    =   1
      cc( 1, 2) =         29.16521900d0
      cc( 2, 2) =         -4.19080600d0
      a ( 1, 2) =          1.43509100d0
      a ( 2, 2) =          0.69580500d0
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      nt( 2)    =   2
      cc( 1, 3) =         36.99054200d0
      cc( 2, 3) =         -3.36582000d0
      a ( 1, 3) =          1.44083200d0
      cc( 1, 4) =         23.65207854d0
      cc( 2, 4) =          3.25844113d0
      a ( 1, 4) =          2.12260500d0
      a ( 2, 4) =          0.79866900d0
      cc( 1, 5) =        -47.70319876d0
      cc( 2, 5) =         -6.54113991d0
      a ( 1, 5) =          6.16436000d0
      a ( 2, 5) =          1.54237400d0
      cc( 1, 6) =         -7.10585060d0
      a ( 1, 6) =          1.84789200d0
c
      a ( 2, 3) =          0.70139200d0
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      nt( 3)    =   2
      cc( 1, 4) =         20.00053100d0
      a ( 1, 4) =          0.96123600d0
      n ( 1, 4) =   2
      nt( 4)    =   1
      cc( 1, 5) =         -6.01909200d0
      a ( 1, 5) =          0.88436900d0
      n ( 1, 5) =   2
      nt( 5)    =   1
      nelcor    =  46
      lmx       =   4
c
      return
c
c     tin (Sn_Stuttgart_RLC_ECP)
c
  50  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      n ( 1, 1) =   2
      nt( 1)    =   1
      cc( 1, 2) =         67.92534700d0
      cc( 2, 2) =         -7.47854600d0
      a ( 1, 2) =          1.96972500d0
      a ( 2, 2) =          0.97237500d0
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      nt( 2)    =   2
      cc( 1, 3) =         56.60288000d0
      cc( 2, 3) =         -2.16177600d0
      a ( 1, 3) =          1.99921000d0
      a ( 2, 3) =          0.99904200d0
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      nt( 3)    =   2
      cc( 1, 4) =          2.57633600d0
      a ( 1, 4) =          0.50036100d0
      n ( 1, 4) =   2
      nt( 4)    =   1
      cc( 1, 5) =        -10.10925300d0
      a ( 1, 5) =          1.23088000d0
      n ( 1, 5) =   2
      nt( 5)    =   1
      nelcor    =  46
      lmx       =   4
c
      return
c
c     antimony (Sb_Stuttgart_RLC_ECP)
c
  51  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      n ( 1, 1) =   2
      nt( 1)    =   1
      cc( 1, 2) =         68.42793800d0
      cc( 2, 2) =         -4.39863100d0
      a ( 1, 2) =          2.49109100d0
      a ( 2, 2) =          1.34157500d0
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      nt( 2)    =   2
      cc( 1, 3) =         63.96546900d0
      cc( 2, 3) =         -0.57872600d0
      a ( 1, 3) =          2.14386400d0
      a ( 2, 3) =          0.58550300d0
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      nt( 3)    =   2
      cc( 1, 4) =          7.80366100d0
      a ( 1, 4) =          0.79540100d0
      n ( 1, 4) =   2
      nt( 4)    =   1
      cc( 1, 5) =        -14.51768700d0
      a ( 1, 5) =          1.60925100d0
      n ( 1, 5) =   2
      nt( 5)    =   1
      nelcor    =  46
      lmx       =   4
c
      return
c
c     tellurium (Te_Stuttgart_RLC_ECP)
c
  52  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      n ( 1, 1) =   2
      nt( 1)    =   1
      cc( 1, 2) =         50.08380500d0
      cc( 2, 2) =          1.96814000d0
      a ( 1, 2) =          2.92379400d0
      a ( 2, 2) =          1.15275400d0
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      nt( 2)    =   2
      cc( 1, 3) =        119.82070200d0
      cc( 2, 3) =         -2.03904800d0
      a ( 1, 3) =          2.60308600d0
      a ( 2, 3) =          0.98544800d0
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      nt( 3)    =   2
      cc( 1, 4) =         37.75721400d0
      a ( 1, 4) =          1.43501900d0
      n ( 1, 4) =   2
      nt( 4)    =   1
      cc( 1, 5) =        -17.86464100d0
      a ( 1, 5) =          1.93927000d0
      n ( 1, 5) =   2
      nt( 5)    =   1
      nelcor    =  46
      lmx       =   4
c
      return
c
c     iodine (I_Stuttgart_RLC_ECP)
c
  53  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      n ( 1, 1) =   2
      nt( 1)    =   1
      cc( 1, 2) =         83.09814545d0
      cc( 2, 2) =          5.06370919d0
      a ( 1, 2) =          3.50642001d0
      a ( 2, 2) =          1.74736492d0
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      nt( 2)    =   2
      cc( 1, 3) =         27.29481509d0
      cc( 2, 3) =         55.60853601d0
      cc( 3, 3) =          0.77464159d0
      cc( 4, 3) =          1.81386562d0
      a ( 1, 3) =          2.99860773d0
      a ( 2, 3) =          3.01690894d0
      a ( 3, 3) =          1.59415934d0
      a ( 4, 3) =          1.19802939d0
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      n ( 3, 3) =   2
      n ( 4, 3) =   2
      nt( 3)    =   4
      cc( 1, 4) =          2.56052702d0
      cc( 2, 4) =          3.72797296d0
      cc( 3, 4) =          7.64641669d0
      cc( 4, 4) =         11.45079545d0
      a ( 1, 4) =          1.03813792d0
      a ( 2, 4) =          1.01158599d0
      a ( 3, 4) =          2.04193864d0
      a ( 4, 4) =          1.99631017d0
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      n ( 3, 4) =   2
      n ( 4, 4) =   2
      nt( 4)    =   4
      cc( 1, 5) =        -10.62474188d0
      cc( 2, 5) =        -14.27512750d0
      cc( 3, 5) =         -0.11972820d0
      cc( 4, 5) =         -0.40105292d0
      a ( 1, 5) =          2.64971585d0
      a ( 2, 5) =          2.75335574d0
      a ( 3, 5) =          0.49970082d0
      a ( 4, 5) =          0.79638982d0
      n ( 1, 5) =   2
      n ( 2, 5) =   2
      n ( 3, 5) =   2
      n ( 4, 5) =   2
      nt( 5)    =   4
      nelcor    =  46
      lmx       =   4
c
      return
c
c     xeon (Xe_Stuttgart_RLC_ECP)
c
  54  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        122.76382934d0
      cc( 2, 2) =          8.30885115d0
      a ( 1, 2) =          3.94026300d0
      a ( 2, 2) =          2.27726400d0
      cc( 1, 3) =         68.82300437d0
      cc( 2, 3) =          3.64674223d0
      a ( 1, 3) =          3.02837300d0
      a ( 2, 3) =          1.39431900d0
      cc( 1, 4) =         23.65207854d0
      cc( 2, 4) =          3.25844113d0
      a ( 1, 4) =          2.12260500d0
      a ( 2, 4) =          0.79866900d0
      cc( 1, 5) =        -47.70319876d0
      cc( 2, 5) =         -6.54113991d0
      a ( 1, 5) =          6.16436000d0
      a ( 2, 5) =          1.54237400d0
c     remove g potentials from expansion
c     cc( 1, 6) =         -7.10585060d0
c     a ( 1, 6) =          1.84789200d0
c
      n ( 1, 1) =   2
      nt( 1)    =   1
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      nt( 2)    =   2
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      nt( 3)    =   2
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      nt( 4)    =   2
      n ( 1, 5) =   2
      n ( 2, 5) =   2
      nt( 5)    =   2
c     n ( 1, 6) =   2
c     nt( 6)    =   1
      nelcor    =  46
c     lmx       =   5
      lmx       =   4
c
      return
 6010 format (1x,'*** data not found for Stuttgart_RLC ecp type ',a4)
      end
**==strlc3.f
      subroutine strlc3 (zinp,lmx,nelcor,nt,n,cc,a)
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
      implicit character * 8 (z)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      dimension nt(6), n(10,6), cc(10,6), a(10,6)
c
c     this routine contains a library of ecp data taken from
c     various literature sources.
c     Stuttgart Relativistic, Large Core ECP
c
      dimension ztit(17)
      data maxtyp/17/
      data ztit/'cs','ba','la',
     +          'ce','pr','nd','pm','sm','eu','gd','tb',
     +          'dy','ho','er','tm','yb','lu' /
c
c
      do 10 ityp = 1 , maxtyp
         if (zinp.eq.ztit(ityp)) go to 20
 10   continue
 57   write (iwr,6010) zinp
      call caserr2('no library available for requested ecp')
c
 20   continue
c
      go to (55,56,57,57,57,57,57,57,57,57,57,57,57,57,57,57,57), ityp
c
c   cesium (Cs_Stuttgart_RLC_ECP)
c
 55   continue
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        303.50000000d0
      a ( 1, 2) =          1.16000000d0
      cc( 1, 3) =          2.93000000d0
      a ( 1, 3) =          0.24470000d0
      cc( 1, 4) =         -1.13800000d0
      a ( 1, 4) =          0.24270000d0
c
      go to 60
c
c   barium (Ba_Stuttgart_RLC_ECP)
c
 56   continue
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =         16.71000000d0
      a ( 1, 2) =          0.63560000d0
      cc( 1, 3) =          4.87400000d0
      a ( 1, 3) =          0.32140000d0
      cc( 1, 4) =         -0.84270000d0
      a ( 1, 4) =          0.24030000d0
  60  continue
      n ( 1, 1) =   2
      nt( 1)    =   1
      n ( 1, 2) =   2
      nt( 2)    =   1
      n ( 1, 3) =   2
      nt( 3)    =   1
      n ( 1, 4) =   2
      nt( 4)    =   1
      nelcor    =  54
      lmx       =   3
c
      return
 6010 format (1x,'*** data not found for Stuttgart_RLC ecp type ',a4)
      end
**==strlc4.f
      subroutine strlc4 (zinp,lmx,nelcor,nt,n,cc,a)
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
      implicit character * 8 (z)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      dimension nt(6), n(10,6), cc(10,6), a(10,6)
c
c     this routine contains a library of ecp data taken from
c     various literature sources.
c     Stuttgart Relativistic, Large Core ECP
c
      dimension ztit(15)
      data maxtyp/15/
      data ztit/
     $ 'hf','ta','w ','re','os','ir','pt','au','hg','tl','pb',
     + 'bi','po','at','rn'/
c
      do 10 ityp = 1 , maxtyp
         if (zinp.eq.ztit(ityp)) go to 20
 10   continue
 72   write (iwr,6010) zinp
      call caserr2('no library available for requested ecp')
c
 20   go to(72,72,72,72,72,72,72,72,80,81,
     +      82,83,84,85,86), ityp
c
c     mercury (Hg_Stuttgart_RLC_ECP)
c
  80  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      n ( 1, 1) =   2
      nt( 1)    =   1
      cc( 1, 2) =         -0.69617800d0
      cc( 2, 2) =         27.75810500d0
      cc( 3, 2) =         48.78047500d0
      a ( 1, 2) =          0.22721000d0
      a ( 2, 2) =          1.65753000d0
      a ( 3, 2) =         10.00024800d0
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      nt( 2)    =   3
      cc( 1, 3) =         -2.73581100d0
      cc( 2, 3) =          8.57563700d0
      a ( 1, 3) =          0.39837700d0
      a ( 2, 3) =          0.64730700d0
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      nt( 3)    =   2
      cc( 1, 4) =         -0.01311800d0
      cc( 2, 4) =          2.79286200d0
      a ( 1, 4) =          0.21799900d0
      a ( 2, 4) =          0.38605800d0
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      nt( 4)    =   2
      cc( 1, 5) =         -2.63516400d0
      a ( 1, 5) =          0.50000000d0
      n ( 1, 5) =   2
      nt( 5)    =   1
c     remove g potentials from expansion
c     cc( 1, 6) =        -13.39371600d0
c     a ( 1, 6) =          0.80075600d0
c     n ( 1, 6) =   2
c     nt( 6)    =   1
      nelcor    =  78
c     lmx       =   5
      lmx       =   4
c
      go to 30
c
c     Tl (Tl_CRENBL_ECP)
c
  81  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      n ( 1, 1) =   2
      nt( 1)    =   1
      cc( 1, 2) =         -1.01649800d0
      cc( 2, 2) =         51.70795900d0
      cc( 3, 2) =         73.18668300d0
      a ( 1, 2) =          0.32623800d0
      a ( 2, 2) =          1.97754100d0
      a ( 3, 2) =         10.00000000d0
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      nt( 2)    =   3
      cc( 1, 3) =         -2.96267300d0
      cc( 2, 3) =         19.73043100d0
      a ( 1, 3) =          0.54306300d0
      a ( 2, 3) =          1.03214000d0
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      nt( 3)    =   2
      cc( 1, 4) =          2.77269000d0
      cc( 2, 4) =         -3.97943900d0
      a ( 1, 4) =          0.35481700d0
      a ( 2, 4) =          0.70963300d0
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      nt( 4)    =   2
      cc( 1, 5) =         -4.42678600d0
      a ( 1, 5) =          0.68915600d0
      n ( 1, 5) =   2
      nt( 5)    =   1
c     remove g potentials from expansion
c     cc( 1, 6) =        -12.27054000d0
c     a ( 1, 6) =          0.82061700d0
c     n ( 1, 6) =   2
c     nt( 6)    =   1
      nelcor    =  78
c     lmx       =   5
      lmx       =   4
c
      go to 30
c
c     lead  (Pb_Stuttgart_RLC_ECP)
c
  82  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      n ( 1, 1) =   2
      nt( 1)    =   1
      cc( 1, 2) =         -1.87334200d0
      cc( 2, 2) =         20.86079700d0
      cc( 3, 2) =         97.58795500d0
      a ( 1, 2) =          0.52916100d0
      a ( 2, 2) =          1.45672700d0
      a ( 3, 2) =          9.99991100d0
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      nt( 2)    =   3
      cc( 1, 3) =         -7.76820900d0
      cc( 2, 3) =         51.71925400d0
      a ( 1, 3) =          0.67811900d0
      a ( 2, 3) =          1.24901300d0
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      nt( 3)    =   2
      cc( 1, 4) =          1.30076000d0
      cc( 2, 4) =          2.64082200d0
      a ( 1, 4) =          0.30744600d0
      a ( 2, 4) =          0.74493000d0
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      nt( 4)    =   2
      cc( 1, 5) =         -5.70605600d0
      a ( 1, 5) =          0.84869900d0
      n ( 1, 5) =   2
      nt( 5)    =   1
c     remove g potentials from expansion
c     cc( 1, 6) =         -7.48418400d0
c     a ( 1, 6) =          0.99994100d0
c     n ( 1, 6) =   2
c     nt( 6)    =   1
      nelcor    =  78
c     lmx       =   5
      lmx       =   4
c
      go to 30
c
c     bismuth (Bi_Stuttgart_RLC_ECP)
c
   83 continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      n ( 1, 1) =   2
      nt( 1)    =   1
      cc( 1, 2) =         -0.16198800d0
      cc( 2, 2) =         14.03169000d0
      cc( 3, 2) =        122.04740100d0
      a ( 1, 2) =          0.16115200d0
      a ( 2, 2) =          1.50983500d0
      a ( 3, 2) =         10.00000000d0
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      nt( 2)    =   3
      cc( 1, 3) =         -6.18852600d0
      cc( 2, 3) =         51.04586800d0
      a ( 1, 3) =          0.76049000d0
      a ( 2, 3) =          1.42641500d0
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      nt( 3)    =   2
      cc( 1, 4) =         20.53580400d0
      cc( 2, 4) =         -0.13619600d0
      a ( 1, 4) =          0.78022600d0
      a ( 2, 4) =          0.26007500d0
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      nt( 4)    =   2
      cc( 1, 5) =         -6.41422600d0
      a ( 1, 5) =          0.97360800d0
      n ( 1, 5) =   2
      nt( 5)    =   1
c     remove g potentials from expansion
c     cc( 1, 6) =         -6.65606400d0
c     a ( 1, 6) =          1.08819500d0
c     n ( 1, 6) =   2
c     nt( 6)    =   1
      nelcor    =  78
c     lmx       =   5
      lmx       =   4
c
      go to 30
c
c     Po (Po_Stuttgart_RLC_ECP)
c
   84 continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      n ( 1, 1) =   2
      nt( 1)    =   1
      cc( 1, 2) =         -4.15930400d0
      cc( 2, 2) =         33.83035400d0
      cc( 3, 2) =        146.33910100d0
      a ( 1, 2) =          0.92238600d0
      a ( 2, 2) =          1.78191500d0
      a ( 3, 2) =         10.00000000d0
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      nt( 2)    =   3
      cc( 1, 3) =         -4.12531100d0
      cc( 2, 3) =         35.00707800d0
      a ( 1, 3) =          0.72429100d0
      a ( 2, 3) =          1.36386000d0
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      nt( 3)    =   2
      cc( 1, 4) =          1.20651800d0
      cc( 2, 4) =         13.35612500d0
      a ( 1, 4) =          0.47697900d0
      a ( 2, 4) =          0.95395700d0
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      nt( 4)    =   2
      cc( 1, 5) =         -6.77517400d0
      a ( 1, 5) =          1.07545400d0
      n ( 1, 5) =   2
      nt( 5)    =   1
c     remove g potentials from expansion
c     cc( 1, 6) =         -5.51544100d0
c     a ( 1, 6) =          1.12209600d0
c     n ( 1, 6) =   2
c     nt( 6)    =   1
      nelcor    =  78
c     lmx       =   5
      lmx       =   4
c
      go to 30
c
c     at (At_Stuttgart_RLC_ECP)
c
   85 continue
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      n ( 1, 1) =   2
      nt( 1)    =   1
      cc( 1, 2) =         -5.52846100d0
      cc( 2, 2) =         39.56886900d0
      cc( 3, 2) =        170.71138600d0
      a ( 1, 2) =          0.92238600d0
      a ( 2, 2) =          1.78191500d0
      a ( 3, 2) =         10.00000000d0
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      nt( 2)    =   3
      cc( 1, 3) =         -2.29538700d0
      cc( 2, 3) =         25.49292000d0
      a ( 1, 3) =          0.72429100d0
      a ( 2, 3) =          1.36386000d0
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      nt( 3)    =   2
      cc( 1, 4) =          4.86510700d0
      cc( 2, 4) =         14.57941300d0
      a ( 1, 4) =          0.63597200d0
      a ( 2, 4) =          1.27194300d0
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      nt( 4)    =   2
      cc( 1, 5) =         -6.85786700d0
      a ( 1, 5) =          1.15410800d0
      n ( 1, 5) =   2
      nt( 5)    =   1
c     remove g potentials from expansion
c     cc( 1, 6) =         -7.61303900d0
c     a ( 1, 6) =          1.34511500d0
c     n ( 1, 6) =   2
c     nt( 6)    =   1
      nelcor    =  78
c     lmx       =   5
      lmx       =   4
c
      go to 30
c
c     radon (Rn_Stuttgart_RLC_ECP)
c
   86 continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      n ( 1, 1) =   2
      nt( 1)    =   1
      cc( 1, 2) =         -5.01900500d0
      cc( 2, 2) =         37.03679000d0
      cc( 3, 2) =        195.10330800d0
      a ( 1, 2) =          0.92238600d0
      a ( 2, 2) =          1.78191500d0
      a ( 3, 2) =         10.80460100d0
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      nt( 2)    =   3
      cc( 1, 3) =         -1.96648100d0
      cc( 2, 3) =         23.46405900d0
      a ( 1, 3) =          0.72429100d0
      a ( 2, 3) =          1.36386000d0
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      nt( 3)    =   2
      cc( 1, 4) =          7.48345700d0
      cc( 2, 4) =          9.36190000d0
      a ( 1, 4) =          0.76940000d0
      a ( 2, 4) =          1.53880000d0
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      nt( 4)    =   2
      cc( 1, 5) =         -6.76315000d0
      a ( 1, 5) =          1.21389700d0
      n ( 1, 5) =   2
      nt( 5)    =   1
c     remove g potentials from expansion
c     cc( 1, 6) =         -9.91566200d0
c     a ( 1, 6) =          1.57646900d0
c     n ( 1, 6) =   2
c     nt( 6)    =   1
      nelcor    =  78
c     lmx       =   5
      lmx       =   4
c
 30   continue
      return
c
 6010 format (1x,'*** data not found for Stuttgart_RLC ecp type ',a4)
      end
**==strlc5.f
      subroutine strlc5 (zinp,lmx,nelcor,nt,n,cc,a)
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
      implicit character * 8 (z)
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
      dimension nt(6), n(10,6), cc(10,6), a(10,6)
c
c     this routine contains a library of ecp data taken from
c     various literature sources.
c
      dimension ztit(17)
      data maxtyp/17/
      data ztit/
     +  'fr','ra','ac','th','pa','u ','np','pu','am','cm',
     $  'bk','cf','es','fm','md','no','lw'   /
c
c     this routine contains a library of ecp data taken from
c     various literature sources.
c     Actinides Stuttgart Relativistic, Large Core ECP
c     
      do 10 ityp = 1 , maxtyp
         if (zinp.eq.ztit(ityp)) go to 20
 10   continue
 87   write (iwr,6010) zinp
      call caserr2('no library available for requested ecp')
c
 20   go to(87,87,89,90,91,92,93,94,95,96,97,98,99,100,
     +     101,102,103), ityp
c
c     ac (Ac_Stuttgart_RLC_ECP)
c
  89  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      n ( 1, 1) =   2
      nt( 1)    =   1
      cc( 1, 2) =        113.05439700d0
      cc( 2, 2) =         17.08323300d0
      cc( 3, 2) =         -2.66661500d0
      a ( 1, 2) =          4.06365300d0
      a ( 2, 2) =          1.88399500d0
      a ( 3, 2) =          0.88656700d0
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      nt( 2)    =   3
      cc( 1, 3) =        115.93765500d0
      cc( 2, 3) =         18.32997800d0
      cc( 3, 3) =          1.30743400d0
      a ( 1, 3) =          3.98618100d0
      a ( 2, 3) =          2.00016000d0
      a ( 3, 3) =          0.96084100d0
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      n ( 3, 3) =   2
      nt( 3)    =   3
      cc( 1, 4) =         60.86374000d0
      cc( 2, 4) =         29.38228900d0
      cc( 3, 4) =          6.96985000d0
      a ( 1, 4) =          4.14797200d0
      a ( 2, 4) =          2.23456300d0
      a ( 3, 4) =          0.91369500d0
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      n ( 3, 4) =   2
      nt( 4)    =   3
      cc( 1, 5) =         48.64607100d0
      cc( 2, 5) =        -29.41016000d0
      cc( 3, 5) =          2.55310600d0
      a ( 1, 5) =          3.99893800d0
      a ( 2, 5) =          1.99884000d0
      a ( 3, 5) =          0.99564100d0
      n ( 1, 5) =   2
      n ( 2, 5) =   2
      n ( 3, 5) =   2
      nt( 5)    =   3
      cc( 1, 6) =        -25.14020900d0
      cc( 2, 6) =         -0.34197000d0
      a ( 1, 6) =          2.85859900d0
      a ( 2, 6) =          0.88750200d0
      n ( 1, 6) =   2
      n ( 2, 6) =   2
      nt( 6)    =   2
c     lmx       =   5
c
      go to 30
c
c     th (Th_Stuttgart_RLC_ECP)
c
  90  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      n ( 1, 1) =   2
      nt( 1)    =   1
      cc( 1, 2) =        113.32646600d0
      cc( 2, 2) =         15.66375500d0
      cc( 3, 2) =         -2.76590200d0
      a ( 1, 2) =          4.06365300d0
      a ( 2, 2) =          1.88399500d0
      a ( 3, 2) =          0.88656700d0
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      nt( 2)    =   3
      cc( 1, 3) =        115.95348300d0
      cc( 2, 3) =         15.76219000d0
      cc( 3, 3) =          1.37885000d0
      a ( 1, 3) =          3.98618100d0
      a ( 2, 3) =          2.00016000d0
      a ( 3, 3) =          0.96084100d0
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      n ( 3, 3) =   2
      nt( 3)    =   3
      cc( 1, 4) =         59.74781100d0
      cc( 2, 4) =         17.82068300d0
      cc( 3, 4) =          6.91366100d0
      a ( 1, 4) =          4.14797200d0
      a ( 2, 4) =          2.23456300d0
      a ( 3, 4) =          0.91369500d0
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      n ( 3, 4) =   2
      nt( 4)    =   3
      cc( 1, 5) =         49.62437800d0
      cc( 2, 5) =        -26.64186100d0
      cc( 3, 5) =          1.76274800d0
      a ( 1, 5) =          3.99893800d0
      a ( 2, 5) =          1.99884000d0
      a ( 3, 5) =          0.99564100d0
      n ( 1, 5) =   2
      n ( 2, 5) =   2
      n ( 3, 5) =   2
      nt( 5)    =   3
      cc( 1, 6) =        -28.59594100d0
      cc( 2, 6) =         -0.30965300d0
      a ( 1, 6) =          3.16637900d0
      a ( 2, 6) =          0.86010500d0
      n ( 1, 6) =   2
      n ( 2, 6) =   2
      nt( 6)    =   2
c     lmx       =   5
c
      go to 30
c
c     pa (Pa_Stuttgart_RLC_ECP)
c
  91  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      n ( 1, 1) =   2
      nt( 1)    =   1
      cc( 1, 2) =        113.29484400d0
      cc( 2, 2) =         15.51386300d0
      cc( 3, 2) =         -3.14394300d0
      a ( 1, 2) =          4.06365300d0
      a ( 2, 2) =          1.88399500d0
      a ( 3, 2) =          0.88656700d0
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      nt( 2)    =   3
      cc( 1, 3) =        115.89528800d0
      cc( 2, 3) =         15.48629100d0
      cc( 3, 3) =          0.97628300d0
      a ( 1, 3) =          3.98618100d0
      a ( 2, 3) =          2.00016000d0
      a ( 3, 3) =          0.96084100d0
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      n ( 3, 3) =   2
      nt( 3)    =   3
      cc( 1, 4) =         59.66921300d0
      cc( 2, 4) =         17.34372600d0
      cc( 3, 4) =          6.44469800d0
      a ( 1, 4) =          4.14797200d0
      a ( 2, 4) =          2.23456300d0
      a ( 3, 4) =          0.91369500d0
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      n ( 3, 4) =   2
      nt( 4)    =   3
      cc( 1, 5) =         51.67708300d0
      cc( 2, 5) =        -22.94129000d0
      cc( 3, 5) =          0.45025800d0
      a ( 1, 5) =          3.99893800d0
      a ( 2, 5) =          1.99884000d0
      a ( 3, 5) =          0.99564100d0
      n ( 1, 5) =   2
      n ( 2, 5) =   2
      n ( 3, 5) =   2
      nt( 5)    =   3
      cc( 1, 6) =        -31.87331800d0
      cc( 2, 6) =         -0.31139900d0
      a ( 1, 6) =          3.48065000d0
      a ( 2, 6) =          0.86010500d0
      n ( 1, 6) =   2
      n ( 2, 6) =   2
      nt( 6)    =   2
c     lmx       =   5
c
      go to 30
c
c     u (U_Stuttgart_RLC_ECP)
c
  92  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      n ( 1, 1) =   2
      nt( 1)    =   1
      cc( 1, 2) =        112.92010300d0
      cc( 2, 2) =         15.64750000d0
      cc( 3, 2) =         -3.68997100d0
      a ( 1, 2) =          4.06365300d0
      a ( 2, 2) =          1.88399500d0
      a ( 3, 2) =          0.88656700d0
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      nt( 2)    =   3
      cc( 1, 3) =        118.75801600d0
      cc( 2, 3) =         15.07722800d0
      cc( 3, 3) =          0.55672000d0
      a ( 1, 3) =          3.98618100d0
      a ( 2, 3) =          2.00016000d0
      a ( 3, 3) =          0.96084100d0
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      n ( 3, 3) =   2
      nt( 3)    =   3
      cc( 1, 4) =         60.85589200d0
      cc( 2, 4) =         29.28004700d0
      cc( 3, 4) =          4.99802900d0
      a ( 1, 4) =          4.14797200d0
      a ( 2, 4) =          2.23456300d0
      a ( 3, 4) =          0.91369500d0
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      n ( 3, 4) =   2
      nt( 4)    =   3
      cc( 1, 5) =         49.92403500d0
      cc( 2, 5) =        -24.67404200d0
      cc( 3, 5) =          1.38948000d0
      a ( 1, 5) =          3.99893800d0
      a ( 2, 5) =          1.99884000d0
      a ( 3, 5) =          0.99564100d0
      n ( 1, 5) =   2
      n ( 2, 5) =   2
      n ( 3, 5) =   2
      nt( 5)    =   3
      cc( 1, 6) =        -36.04097700d0
      cc( 2, 6) =         -0.09051100d0
      a ( 1, 6) =          3.81742200d0
      a ( 2, 6) =          0.26250100d0
      n ( 1, 6) =   2
      n ( 2, 6) =   2
      nt( 6)    =   2
c     lmx       =   5
c
      go to 30
c
c     Np (Np_Stuttgart_RLC_ECP)
c
  93  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      n ( 1, 1) =   2
      nt( 1)    =   1
      cc( 1, 2) =        112.93486000d0
      cc( 2, 2) =         15.54943200d0
      cc( 3, 2) =         -4.13585700d0
      a ( 1, 2) =          4.06365300d0
      a ( 2, 2) =          1.88399500d0
      a ( 3, 2) =          0.88656700d0
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      nt( 2)    =   3
      cc( 1, 3) =        118.82359400d0
      cc( 2, 3) =         15.28048700d0
      cc( 3, 3) =          0.02369600d0
      a ( 1, 3) =          3.98618100d0
      a ( 2, 3) =          2.00016000d0
      a ( 3, 3) =          0.96084100d0
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      n ( 3, 3) =   2
      nt( 3)    =   3
      cc( 1, 4) =         60.81340100d0
      cc( 2, 4) =         28.99230800d0
      cc( 3, 4) =          4.44499000d0
      a ( 1, 4) =          4.14797200d0
      a ( 2, 4) =          2.23456300d0
      a ( 3, 4) =          0.91369500d0
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      n ( 3, 4) =   2
      nt( 4)    =   3
      cc( 1, 5) =         49.54043100d0
      cc( 2, 5) =        -25.56675700d0
      cc( 3, 5) =          1.88888600d0
      a ( 1, 5) =          3.99893800d0
      a ( 2, 5) =          1.99884000d0
      a ( 3, 5) =          0.99564100d0
      n ( 1, 5) =   2
      n ( 2, 5) =   2
      n ( 3, 5) =   2
      nt( 5)    =   3
      cc( 1, 6) =        -39.83411400d0
      cc( 2, 6) =         -0.14165500d0
      a ( 1, 6) =          4.20153100d0
      a ( 2, 6) =          0.34507800d0
      n ( 1, 6) =   2
      n ( 2, 6) =   2
      nt( 6)    =   2
c     lmx       =   5
c
      go to 30
c
c     Pu (Pu_Stuttgart_RLC_ECP)
c
  94  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      n ( 1, 1) =   2
      nt( 1)    =   1
      cc( 1, 2) =        113.80667800d0
      cc( 2, 2) =         16.00076100d0
      cc( 3, 2) =         -4.83994900d0
      a ( 1, 2) =          4.06365300d0
      a ( 2, 2) =          1.88399500d0
      a ( 3, 2) =          0.88656700d0
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      nt( 2)    =   3
      cc( 1, 3) =        118.76901400d0
      cc( 2, 3) =         15.00546500d0
      cc( 3, 3) =         -0.44582000d0
      a ( 1, 3) =          3.98618100d0
      a ( 2, 3) =          2.00016000d0
      a ( 3, 3) =          0.96084100d0
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      n ( 3, 3) =   2
      nt( 3)    =   3
      cc( 1, 4) =         60.83801500d0
      cc( 2, 4) =         29.14099100d0
      cc( 3, 4) =          3.92557300d0
      a ( 1, 4) =          4.14797200d0
      a ( 2, 4) =          2.23456300d0
      a ( 3, 4) =          0.91369500d0
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      n ( 3, 4) =   2
      nt( 4)    =   3
      cc( 1, 5) =         49.57833400d0
      cc( 2, 5) =        -25.37099200d0
      cc( 3, 5) =          1.92996500d0
      a ( 1, 5) =          3.99893800d0
      a ( 2, 5) =          1.99884000d0
      a ( 3, 5) =          0.99564100d0
      n ( 1, 5) =   2
      n ( 2, 5) =   2
      n ( 3, 5) =   2
      nt( 5)    =   3
      cc( 1, 6) =        -32.53511100d0
      cc( 2, 6) =          0.23591700d0
      a ( 1, 6) =          3.47990400d0
      a ( 2, 6) =          0.16371800d0
      n ( 1, 6) =   2
      n ( 2, 6) =   2
      nt( 6)    =   2
c     lmx       =   5
c
      go to 30
c
c     Am (AM_Stuttgart_RLC_ECP)
c
  95  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      n ( 1, 1) =   2
      nt( 1)    =   1
      cc( 1, 2) =        153.79716900d0
      cc( 2, 2) =        -10.63985500d0
      cc( 3, 2) =          1.27070900d0
      a ( 1, 2) =          4.06365300d0
      a ( 2, 2) =          1.88399500d0
      a ( 3, 2) =          0.88656700d0
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      nt( 2)    =   3
      cc( 1, 3) =        118.77906400d0
      cc( 2, 3) =         14.95569600d0
      cc( 3, 3) =         -1.04628000d0
      a ( 1, 3) =          3.98618100d0
      a ( 2, 3) =          2.00016000d0
      a ( 3, 3) =          0.96084100d0
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      n ( 3, 3) =   2
      nt( 3)    =   3
      cc( 1, 4) =         60.83377800d0
      cc( 2, 4) =         29.09948900d0
      cc( 3, 4) =          3.33361000d0
      a ( 1, 4) =          4.14797200d0
      a ( 2, 4) =          2.23456300d0
      a ( 3, 4) =          0.91369500d0
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      n ( 3, 4) =   2
      nt( 4)    =   3
      cc( 1, 5) =         48.97113700d0
      cc( 2, 5) =        -26.59585400d0
      cc( 3, 5) =          2.47738400d0
      a ( 1, 5) =          3.99893800d0
      a ( 2, 5) =          1.99884000d0
      a ( 3, 5) =          0.99564100d0
      n ( 1, 5) =   2
      n ( 2, 5) =   2
      n ( 3, 5) =   2
      nt( 5)    =   3
      cc( 1, 6) =        -37.95557200d0
      cc( 2, 6) =          0.24965900d0
      a ( 1, 6) =          3.94020500d0
      a ( 2, 6) =          0.17460900d0
      n ( 1, 6) =   2
      n ( 2, 6) =   2
      nt( 6)    =   2
c     lmx       =   5
c
      go to 30
c
c     Cm (Cm_Stuttgart_RLC_ECP)
c
  96  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      n ( 1, 1) =   2
      nt( 1)    =   1
      cc( 1, 2) =        253.67174600d0
      cc( 2, 2) =        -11.65215900d0
      cc( 3, 2) =         -3.39677000d0
      a ( 1, 2) =          4.06365300d0
      a ( 2, 2) =          1.88399500d0
      a ( 3, 2) =          0.88656700d0
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      nt( 2)    =   3
      cc( 1, 3) =        118.74048600d0
      cc( 2, 3) =         14.87157000d0
      cc( 3, 3) =         -1.17650900d0
      a ( 1, 3) =          3.98618100d0
      a ( 2, 3) =          2.00016000d0
      a ( 3, 3) =          0.96084100d0
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      n ( 3, 3) =   2
      nt( 3)    =   3
      cc( 1, 4) =         60.90102300d0
      cc( 2, 4) =         29.50038500d0
      cc( 3, 4) =          3.20672900d0
      a ( 1, 4) =          4.14797200d0
      a ( 2, 4) =          2.23456300d0
      a ( 3, 4) =          0.91369500d0
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      n ( 3, 4) =   2
      nt( 4)    =   3
      cc( 1, 5) =         47.23461400d0
      cc( 2, 5) =        -28.71545300d0
      cc( 3, 5) =          3.78352100d0
      a ( 1, 5) =          3.99893800d0
      a ( 2, 5) =          1.99884000d0
      a ( 3, 5) =          0.99564100d0
      n ( 1, 5) =   2
      n ( 2, 5) =   2
      n ( 3, 5) =   2
      nt( 5)    =   3
      cc( 1, 6) =        -39.54389400d0
      cc( 2, 6) =          0.43072600d0
      a ( 1, 6) =          4.08215300d0
      a ( 2, 6) =          0.27418200d0
      n ( 1, 6) =   2
      n ( 2, 6) =   2
      nt( 6)    =   2
c     lmx       =   5
c
      go to 30
c
c     Bk (Bk_Stuttgart_RLC_ECP)
c
  97  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      n ( 1, 1) =   2
      nt( 1)    =   1
      cc( 1, 2) =        203.63782900d0
      cc( 2, 2) =        -11.66776100d0
      cc( 3, 2) =         -2.22777600d0
      a ( 1, 2) =          4.06365300d0
      a ( 2, 2) =          1.88399500d0
      a ( 3, 2) =          0.88656700d0
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      nt( 2)    =   3
      cc( 1, 3) =        118.71541000d0
      cc( 2, 3) =         14.61461600d0
      cc( 3, 3) =         -2.06267900d0
      a ( 1, 3) =          3.98618100d0
      a ( 2, 3) =          2.00016000d0
      a ( 3, 3) =          0.96084100d0
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      n ( 3, 3) =   2
      nt( 3)    =   3
      cc( 1, 4) =         60.82971100d0
      cc( 2, 4) =         29.05802700d0
      cc( 3, 4) =          2.63108000d0
      a ( 1, 4) =          4.14797200d0
      a ( 2, 4) =          2.23456300d0
      a ( 3, 4) =          0.91369500d0
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      n ( 3, 4) =   2
      nt( 4)    =   3
      cc( 1, 5) =         49.01888600d0
      cc( 2, 5) =        -26.01017100d0
      cc( 3, 5) =          2.50941800d0
      a ( 1, 5) =          3.99893800d0
      a ( 2, 5) =          1.99884000d0
      a ( 3, 5) =          0.99564100d0
      n ( 1, 5) =   2
      n ( 2, 5) =   2
      n ( 3, 5) =   2
      nt( 5)    =   3
      cc( 1, 6) =        -44.64169500d0
      cc( 2, 6) =          0.44015300d0
      a ( 1, 6) =          4.54679200d0
      a ( 2, 6) =          0.29525300d0
      n ( 1, 6) =   2
      n ( 2, 6) =   2
      nt( 6)    =   2
c     lmx       =   5
c
      go to 30
c
c     Cf (Cf_Stuttgart_RLC_ECP)
c
  98  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      n ( 1, 1) =   2
      nt( 1)    =   1
      cc( 1, 2) =        263.52760900d0
      cc( 2, 2) =        -12.59516700d0
      cc( 3, 2) =         -5.52088700d0
      a ( 1, 2) =          4.06365300d0
      a ( 2, 2) =          1.88399500d0
      a ( 3, 2) =          0.88656700d0
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      nt( 2)    =   3
      cc( 1, 3) =        118.66891000d0
      cc( 2, 3) =         14.42269600d0
      cc( 3, 3) =         -2.42664400d0
      a ( 1, 3) =          3.98618100d0
      a ( 2, 3) =          2.00016000d0
      a ( 3, 3) =          0.96084100d0
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      n ( 3, 3) =   2
      nt( 3)    =   3
      cc( 1, 4) =         60.81865300d0
      cc( 2, 4) =         28.98946900d0
      cc( 3, 4) =          2.26659100d0
      a ( 1, 4) =          4.14797200d0
      a ( 2, 4) =          2.23456300d0
      a ( 3, 4) =          0.91369500d0
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      n ( 3, 4) =   2
      nt( 4)    =   3
      cc( 1, 5) =         48.41537600d0
      cc( 2, 5) =        -26.42191600d0
      cc( 3, 5) =          2.89169200d0
      a ( 1, 5) =          3.99893800d0
      a ( 2, 5) =          1.99884000d0
      a ( 3, 5) =          0.99564100d0
      n ( 1, 5) =   2
      n ( 2, 5) =   2
      n ( 3, 5) =   2
      nt( 5)    =   3
      cc( 1, 6) =        -52.58700800d0
      cc( 2, 6) =          0.35204400d0
      a ( 1, 6) =          5.23249100d0
      a ( 2, 6) =          0.24375100d0
      n ( 1, 6) =   2
      n ( 2, 6) =   2
      nt( 6)    =   2
c     lmx       =   5
c
      go to 30
c
c     Es (Es_Stuttgart_RLC_ECP)
c
  99  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      n ( 1, 1) =   2
      nt( 1)    =   1
      cc( 1, 2) =        253.52183000d0
      cc( 2, 2) =        -12.75850100d0
      cc( 3, 2) =         -5.80075300d0
      a ( 1, 2) =          4.06365300d0
      a ( 2, 2) =          1.88399500d0
      a ( 3, 2) =          0.88656700d0
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      nt( 2)    =   3
      cc( 1, 3) =        118.67415200d0
      cc( 2, 3) =         14.34295900d0
      cc( 3, 3) =         -2.99442700d0
      a ( 1, 3) =          3.98618100d0
      a ( 2, 3) =          2.00016000d0
      a ( 3, 3) =          0.96084100d0
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      n ( 3, 3) =   2
      nt( 3)    =   3
      cc( 1, 4) =         60.84486900d0
      cc( 2, 4) =         29.12341200d0
      cc( 3, 4) =          1.81106500d0
      a ( 1, 4) =          4.14797200d0
      a ( 2, 4) =          2.23456300d0
      a ( 3, 4) =          0.91369500d0
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      n ( 3, 4) =   2
      nt( 4)    =   3
      cc( 1, 5) =         48.23279100d0
      cc( 2, 5) =        -27.28767500d0
      cc( 3, 5) =          3.32039600d0
      a ( 1, 5) =          3.99893800d0
      a ( 2, 5) =          1.99884000d0
      a ( 3, 5) =          0.99564100d0
      n ( 1, 5) =   2
      n ( 2, 5) =   2
      n ( 3, 5) =   2
      nt( 5)    =   3
      cc( 1, 6) =        -57.40782000d0
      cc( 2, 6) =          0.39380000d0
      a ( 1, 6) =          5.67453900d0
      a ( 2, 6) =          0.27173000d0
      n ( 1, 6) =   2
      n ( 2, 6) =   2
      nt( 6)    =   2
c     lmx       =   5
c
      go to 30
c
c     Fm (Fm_Stuttgart_RLC_ECP)
c
  100 continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      n ( 1, 1) =   2
      nt( 1)    =   1
      cc( 1, 2) =        253.47684100d0
      cc( 2, 2) =        -12.92052300d0
      cc( 3, 2) =         -6.59105600d0
      a ( 1, 2) =          4.06365300d0
      a ( 2, 2) =          1.88399500d0
      a ( 3, 2) =          0.88656700d0
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      nt( 2)    =   3
      cc( 1, 3) =        118.70685700d0
      cc( 2, 3) =         14.36063100d0
      cc( 3, 3) =         -3.66308700d0
      a ( 1, 3) =          3.98618100d0
      a ( 2, 3) =          2.00016000d0
      a ( 3, 3) =          0.96084100d0
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      n ( 3, 3) =   2
      nt( 3)    =   3
      cc( 1, 4) =         60.81616300d0
      cc( 2, 4) =         28.94595500d0
      cc( 3, 4) =          1.40080400d0
      a ( 1, 4) =          4.14797200d0
      a ( 2, 4) =          2.23456300d0
      a ( 3, 4) =          0.91369500d0
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      n ( 3, 4) =   2
      nt( 4)    =   3
      cc( 1, 5) =         47.93239700d0
      cc( 2, 5) =        -27.41319600d0
      cc( 3, 5) =          3.46457700d0
      a ( 1, 5) =          3.99893800d0
      a ( 2, 5) =          1.99884000d0
      a ( 3, 5) =          0.99564100d0
      n ( 1, 5) =   2
      n ( 2, 5) =   2
      n ( 3, 5) =   2
      nt( 5)    =   3
      cc( 1, 6) =        -62.12683600d0
      cc( 2, 6) =          0.43957200d0
      a ( 1, 6) =          6.11905300d0
      a ( 2, 6) =          0.30355100d0
      n ( 1, 6) =   2
      n ( 2, 6) =   2
      nt( 6)    =   2
c     lmx       =   5
c
      go to 30
c
c     Md (Md_Stuttgart_RLC_ECP)
c
  101 continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      n ( 1, 1) =   2
      nt( 1)    =   1
      cc( 1, 2) =        303.34312900d0
      cc( 2, 2) =        -13.76866300d0
      cc( 3, 2) =         -9.54512900d0
      a ( 1, 2) =          4.06365300d0
      a ( 2, 2) =          1.88399500d0
      a ( 3, 2) =          0.88656700d0
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      nt( 2)    =   3
      cc( 1, 3) =        118.59677000d0
      cc( 2, 3) =         13.97306800d0
      cc( 3, 3) =         -4.30323200d0
      a ( 1, 3) =          3.98618100d0
      a ( 2, 3) =          2.00016000d0
      a ( 3, 3) =          0.96084100d0
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      n ( 3, 3) =   2
      nt( 3)    =   3
      cc( 1, 4) =         60.83176900d0
      cc( 2, 4) =         29.00987900d0
      cc( 3, 4) =          0.98395900d0
      a ( 1, 4) =          4.14797200d0
      a ( 2, 4) =          2.23456300d0
      a ( 3, 4) =          0.91369500d0
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      n ( 3, 4) =   2
      nt( 4)    =   3
      cc( 1, 5) =         47.92292900d0
      cc( 2, 5) =        -27.02507800d0
      cc( 3, 5) =          3.38675800d0
      a ( 1, 5) =          3.99893800d0
      a ( 2, 5) =          1.99884000d0
      a ( 3, 5) =          0.99564100d0
      n ( 1, 5) =   2
      n ( 2, 5) =   2
      n ( 3, 5) =   2
      nt( 5)    =   3
      cc( 1, 6) =        -67.00676500d0
      cc( 2, 6) =          0.48622000d0
      a ( 1, 6) =          6.58211900d0
      a ( 2, 6) =          0.33454200d0
      n ( 1, 6) =   2
      n ( 2, 6) =   2
      nt( 6)    =   2
c     lmx       =   5
c
      go to 30
c
c     No (No_Stuttgart_RLC_ECP)
c
  102 continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      n ( 1, 1) =   2
      nt( 1)    =   1
      cc( 1, 2) =        303.28133500d0
      cc( 2, 2) =        -14.01662400d0
      cc( 3, 2) =         -9.93506500d0
      a ( 1, 2) =          4.06365300d0
      a ( 2, 2) =          1.88399500d0
      a ( 3, 2) =          0.88656700d0
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      nt( 2)    =   3
      cc( 1, 3) =        118.53292600d0
      cc( 2, 3) =         13.75407000d0
      cc( 3, 3) =         -4.28920100d0
      a ( 1, 3) =          3.98618100d0
      a ( 2, 3) =          2.00016000d0
      a ( 3, 3) =          0.96084100d0
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      n ( 3, 3) =   2
      nt( 3)    =   3
      cc( 1, 4) =         60.91964600d0
      cc( 2, 4) =         29.43015900d0
      cc( 3, 4) =          0.73652500d0
      a ( 1, 4) =          4.14797200d0
      a ( 2, 4) =          2.23456300d0
      a ( 3, 4) =          0.91369500d0
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      n ( 3, 4) =   2
      nt( 4)    =   3
      cc( 1, 5) =         45.87535300d0
      cc( 2, 5) =        -29.56441300d0
      cc( 3, 5) =          4.86016100d0
      a ( 1, 5) =          3.99893800d0
      a ( 2, 5) =          1.99884000d0
      a ( 3, 5) =          0.99564100d0
      n ( 1, 5) =   2
      n ( 2, 5) =   2
      n ( 3, 5) =   2
      nt( 5)    =   3
      cc( 1, 6) =        -72.34498500d0
      cc( 2, 6) =          0.53059700d0
      a ( 1, 6) =          7.08177000d0
      a ( 2, 6) =          0.35922300d0
      n ( 1, 6) =   2
      n ( 2, 6) =   2
      nt( 6)    =   2
c     lmx       =   5
c
      go to 30

c     Lw (Lw_Stuttgart_RLC_ECP)
c
  103 continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      n ( 1, 1) =   2
      nt( 1)    =   1
      cc( 1, 2) =        353.22492000d0
      cc( 2, 2) =        -14.71008000d0
      cc( 3, 2) =        -12.35925600d0
      a ( 1, 2) =          4.06365300d0
      a ( 2, 2) =          1.88399500d0
      a ( 3, 2) =          0.88656700d0
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      n ( 3, 2) =   2
      nt( 2)    =   3
      cc( 1, 3) =        118.52595600d0
      cc( 2, 3) =         13.72360000d0
      cc( 3, 3) =         -4.44357700d0
      a ( 1, 3) =          3.98618100d0
      a ( 2, 3) =          2.00016000d0
      a ( 3, 3) =          0.96084100d0
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      n ( 3, 3) =   2
      nt( 3)    =   3
      cc( 1, 4) =         61.13011400d0
      cc( 2, 4) =         30.32674300d0
      cc( 3, 4) =          0.46421900d0
      a ( 1, 4) =          4.14797200d0
      a ( 2, 4) =          2.23456300d0
      a ( 3, 4) =          0.91369500d0
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      n ( 3, 4) =   2
      nt( 4)    =   3
      cc( 1, 5) =         44.05127800d0
      cc( 2, 5) =        -31.82235300d0
      cc( 3, 5) =          6.37356600d0
      a ( 1, 5) =          3.99893800d0
      a ( 2, 5) =          1.99884000d0
      a ( 3, 5) =          0.99564100d0
      n ( 1, 5) =   2
      n ( 2, 5) =   2
      n ( 3, 5) =   2
      nt( 5)    =   3
      cc( 1, 6) =        -76.80026400d0
      cc( 2, 6) =          0.58357600d0
      a ( 1, 6) =          7.53958500d0
      a ( 2, 6) =          0.40042600d0
      n ( 1, 6) =   2
      n ( 2, 6) =   2
      nt( 6)    =   2
c
 30   continue
      nelcor    =  78
c     remove g potentials from expansion
c     lmx       =   5
      lmx       =   4
c
      return
c
 6010 format (1x,'*** data not found for Stuttgart_RLC ecp type ',a4)
      end
**==strsc1.f
      subroutine strsc1(zinp,lmx,nelcor,nt,n,cc,a)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character * 8 (z)
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      dimension nt(6), n(10,6), cc(10,6), a(10,6)
c
      dimension ztit(18)
      data maxtyp/18/
      data ztit/'k','ca',
     *          'sc','ti','v ','cr','mn','fe',
     *          'co','ni','cu','zn','ga','ge','as','se','br','kr'/
c
c     this routine contains a library of ecp data taken from
c     various literature sources.
c     Stuttgart Relativistic Small Core ECP

      do 20 ityp = 1 , maxtyp
         if (zinp.eq.ztit(ityp)) go to 30
 20   continue
 160  write (iwr,6010) zinp
      call caserr2('no library available for requested STRSC ecp')
 30   go to (40,50,60,70,80,90,100,110,120,130,140,150,160,160,160,160,
     +       160,160) , ityp
c
c     k   (K_Stuttgart_RSC_ECP)
c
 40   continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =         92.25263900d0
      a ( 1, 2) =          6.84380700d0
      cc( 1, 3) =         23.49803500d0
      a ( 1, 3) =          4.21785700d0
      cc( 1, 4) =         -7.31542000d0
      a ( 1, 4) =          6.49753700d0
      cc( 1, 5) =         -3.85832000d0
      a ( 1, 5) =          6.12361300d0
c
      n ( 1, 1) =   2
      nt( 1)    =   1
      n ( 1, 2) =   2
      nt( 2)    =   1
      n ( 1, 3) =   2
      nt( 3)    =   1
      n ( 1, 4) =   2
      nt( 4)    =   1
      n ( 1, 5) =   2
      nt( 5)    =   1
      nelcor    =  10
      lmx       =   4
c
      return
c
c     ca  (Ca_Stuttgart_RSC_ECP)
c
 50   continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        138.78517400d0
      cc( 2, 2) =         16.50424400d0
      a ( 1, 2) =         11.23167200d0
      a ( 2, 2) =          4.67196000d0
      cc( 1, 3) =         83.12366400d0
      cc( 2, 3) =         13.50227200d0
      a ( 1, 3) =         11.15690700d0
      a ( 2, 3) =          4.81014100d0
      cc( 1, 4) =        -16.20196500d0
      cc( 2, 4) =         -1.13239000d0
      a ( 1, 4) =         13.75472800d0
      a ( 2, 4) =          4.76247000d0
      cc( 1, 5) =        -26.72817800d0
      a ( 1, 5) =         12.76584600d0
c
      n ( 1, 1) =   2
      nt( 1)    =   1
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      nt( 2)    =   2
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      nt( 3)    =   2
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      nt( 4)    =   2
      n ( 1, 5) =   2
      nt( 5)    =   1
      nelcor    =  10
      lmx       =   4
c
      return
c
c     scandium  
c
 60   continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        138.53815200d0
      cc( 2, 2) =         14.83404200d0
      a ( 1, 2) =         11.50000000d0
      a ( 2, 2) =          5.18400000d0
      cc( 1, 3) =         82.45861400d0
      cc( 2, 3) =          8.56520600d0
      a ( 1, 3) =         10.93000000d0
      a ( 2, 3) =          4.58100000d0
      cc( 1, 4) =        -16.12986200d0
      cc( 2, 4) =         -0.53469000d0
      a ( 1, 4) =         13.47000000d0
      a ( 2, 4) =          4.37500000d0
c
      go to 155
c
c     titanium  
c
 70   continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        158.24159300d0
      cc( 2, 2) =         17.51182400d0
      a ( 1, 2) =         13.01000000d0
      a ( 2, 2) =          5.86200000d0
      cc( 1, 3) =         95.23512700d0
      cc( 2, 3) =         10.04785600d0
      a ( 1, 3) =         12.46000000d0
      a ( 2, 3) =          5.21700000d0
      cc( 1, 4) =        -17.56886100d0
      cc( 2, 4) =         -0.58725600d0
      a ( 1, 4) =         15.35000000d0
      a ( 2, 4) =          4.98000000d0
c
      go to 155
c
c     vanadium  
c
 80   continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        178.44797100d0
      cc( 2, 2) =         19.83137500d0
      a ( 1, 2) =         14.49000000d0
      a ( 2, 2) =          6.52400000d0
      cc( 1, 3) =        109.52976300d0
      cc( 2, 3) =         12.57031000d0
      a ( 1, 3) =         14.30000000d0
      a ( 2, 3) =          6.02100000d0
      cc( 1, 4) =        -19.21965700d0
      cc( 2, 4) =         -0.64277500d0
      a ( 1, 4) =         17.48000000d0
      a ( 2, 4) =          5.70900000d0
c
      go to 155
c
c     chromium  
c
 90   continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        201.57888700d0
      cc( 2, 2) =         24.20574100d0
      a ( 1, 2) =         16.39000000d0
      a ( 2, 2) =          7.40200000d0
      cc( 1, 3) =        125.02277400d0
      cc( 2, 3) =         16.47906600d0
      a ( 1, 3) =         16.45000000d0
      a ( 2, 3) =          6.96200000d0
      cc( 1, 4) =        -20.82742100d0
      cc( 2, 4) =         -0.83436800d0
      a ( 1, 4) =         19.93000000d0
      a ( 2, 4) =          6.59800000d0
c
      go to 155
c
c     manganese  
c
 100  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        226.43090200d0
      cc( 2, 2) =         30.35907200d0
      a ( 1, 2) =         18.52000000d0
      a ( 2, 2) =          8.37300000d0
      cc( 1, 3) =        142.15470500d0
      cc( 2, 3) =         21.53650900d0
      a ( 1, 3) =         18.92000000d0
      a ( 2, 3) =          8.01700000d0
      cc( 1, 4) =        -22.56811900d0
      cc( 2, 4) =         -1.20581000d0
      a ( 1, 4) =         22.72000000d0
      a ( 2, 4) =          7.64000000d0
c
      go to 155
c
c     iron  
c
 110  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        253.74958800d0
      cc( 2, 2) =         37.92284500d0
      a ( 1, 2) =         20.93000000d0
      a ( 2, 2) =          9.44500000d0
      cc( 1, 3) =        161.03681200d0
      cc( 2, 3) =         27.65129800d0
      a ( 1, 3) =         21.76000000d0
      a ( 2, 3) =          9.17800000d0
      cc( 1, 4) =        -24.43127600d0
      cc( 2, 4) =         -1.43425100d0
      a ( 1, 4) =         25.90000000d0
      a ( 2, 4) =          8.83500000d0
c
      go to 155
c
c     cobalt  
c
 120  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        283.96056600d0
      cc( 2, 2) =         47.15684600d0
      a ( 1, 2) =         23.66000000d0
      a ( 2, 2) =         10.61000000d0
      cc( 1, 3) =        182.21223600d0
      cc( 2, 3) =         35.23335200d0
      a ( 1, 3) =         25.04000000d0
      a ( 2, 3) =         10.44000000d0
      cc( 1, 4) =        -26.47533300d0
      cc( 2, 4) =         -1.82578700d0
      a ( 1, 4) =         29.54000000d0
      a ( 2, 4) =         10.18000000d0
c
      go to 155 
c
c     nickel (3d)9 (4s)1   ni01
c
 130  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        317.68227200d0
      cc( 2, 2) =         58.25539100d0
      a ( 1, 2) =         26.74000000d0
      a ( 2, 2) =         11.86000000d0
      cc( 1, 3) =        252.47436600d0
      cc( 2, 3) =         36.08150300d0
      a ( 1, 3) =         28.80000000d0
      a ( 2, 3) =         11.79000000d0
      cc( 1, 4) =        -18.52295500d0
      cc( 2, 4) =         -4.55766800d0
      a ( 1, 4) =         33.70000000d0
      a ( 2, 4) =         11.66000000d0
c
      go to 155
c
c     copper  
c
 140  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        355.77015800d0
      cc( 2, 2) =         70.86535700d0
      a ( 1, 2) =         30.22000000d0
      a ( 2, 2) =         13.19000000d0
      cc( 1, 3) =        233.89197600d0
      cc( 2, 3) =         53.94729900d0
      a ( 1, 3) =         33.13000000d0
      a ( 2, 3) =         13.22000000d0
      cc( 1, 4) =        -31.27216500d0
      cc( 2, 4) =         -2.74110400d0
      a ( 1, 4) =         38.42000000d0
      a ( 2, 4) =         13.26000000d0
c
      go to 155
c
c     zinc   (Zn_Stuttgart_RSC_ECP)
c
 150  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        399.98728200d0
      cc( 2, 2) =         85.48565500d0
      a ( 1, 2) =         34.15000000d0
      a ( 2, 2) =         14.59000000d0
      cc( 1, 3) =        277.14896000d0
      cc( 2, 3) =         69.05220500d0
      a ( 1, 3) =         39.78000000d0
      a ( 2, 3) =         14.95000000d0
      cc( 1, 4) =        -34.14934900d0
      cc( 2, 4) =         -3.29183100d0
      a ( 1, 4) =         43.80000000d0
      a ( 2, 4) =         14.98000000d0
c
 155  continue
c
      n ( 1, 1) =   2
      nt( 1)    =   1
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      nt( 2)    =   2
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      nt( 3)    =   2
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      nt( 4)    =   2
      nelcor    =  10
      lmx       =   3
c
      return
 6010 format (1x,'*** data not found for ecp type ',a4)
      end
**==strsc2.f
      subroutine strsc2 (zinp,lmx,nelcor,nt,n,cc,a)
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
      implicit character * 8 (z)
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
      dimension nt(6), n(10,6), cc(10,6), a(10,6)
c
      dimension ztit(18)
      data maxtyp/18/
      data ztit/
     $ 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag','cd',
     $ 'in','sn','sb','te','i ','xe' /
c
c     this routine contains a library of ecp data taken from
c     various literature sources.
c     Stuttgart Relativistic Small Core ECP
c     
      do 10 ityp = 1 , maxtyp
         if (zinp.eq.ztit(ityp)) go to 20
 10   continue
 49   write (iwr,6010) zinp
      call caserr2('no library available for requested ecp')
 20   continue
c
      go to (37,38,39,40,41,42,43,44,45,46,47,48,49,
     +       49,49,49,49,49), ityp
c
c      rb (RB_Stuttgart_RSC_ECP)
c
  37  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =         89.50019822d0
      cc( 2, 2) =          0.49376198d0
      a ( 1, 2) =          5.03655181d0
      a ( 2, 2) =          1.97084963d0
      cc( 1, 3) =         58.56897492d0
      cc( 2, 3) =          0.43179108d0
      a ( 1, 3) =          4.25834119d0
      a ( 2, 3) =          1.47070961d0
      cc( 1, 4) =         26.22489809d0
      cc( 2, 4) =          0.96283950d0
      a ( 1, 4) =          3.02312774d0
      a ( 2, 4) =          0.65038346d0
      cc( 1, 5) =        -12.31690065d0
      a ( 1, 5) =          3.84311471d0
c
      go to 50
c
c      sr (SR_Stuttgart_RSC_ECP)
c
  38  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        135.47943000d0
      cc( 2, 2) =         17.53446300d0
      a ( 1, 2) =          7.40007400d0
      a ( 2, 2) =          3.60637900d0
      cc( 1, 3) =         88.35970900d0
      cc( 2, 3) =         15.39437200d0
      a ( 1, 3) =          6.48486800d0
      a ( 2, 3) =          3.28805300d0
      cc( 1, 4) =         29.88898700d0
      cc( 2, 4) =          6.65941400d0
      a ( 1, 4) =          4.62284100d0
      a ( 2, 4) =          2.24690400d0
      cc( 1, 5) =        -15.80599200d0
      a ( 1, 5) =          4.63397500d0
c
  50  continue
      n ( 1, 1) =   2
      nt( 1)    =   1
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      nt( 2)    =   2
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      nt( 3)    =   2
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      nt( 4)    =   2
      n ( 1, 5) =   2
      nt( 5)    =   1
      nelcor    =  28
      lmx       =   4
c
      return
c
c      yttrium (Y_Stuttgart_RSC_ECP)
c
  39  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        135.15384400d0
      cc( 2, 2) =         15.55244100d0
      a ( 1, 2) =          7.48804900d0
      a ( 2, 2) =          3.74402500d0
      cc( 1, 3) =         87.78499200d0
      cc( 2, 3) =         11.56406600d0
      a ( 1, 3) =          6.44537700d0
      a ( 2, 3) =          3.22268900d0
      cc( 1, 4) =         29.70100100d0
      cc( 2, 4) =          5.53996800d0
      a ( 1, 4) =          4.65844700d0
      a ( 2, 4) =          2.32922400d0
      cc( 1, 5) =        -19.12219800d0
      cc( 2, 5) =         -2.43637500d0
      a ( 1, 5) =          6.58421200d0
      a ( 2, 5) =          3.29210600d0
c
      go to 60
c
c     zirconium (Zr_Stuttgart_RSC_ECP)
c
  40  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        150.26759100d0
      cc( 2, 2) =         18.97621600d0
      a ( 1, 2) =          8.20000000d0
      a ( 2, 2) =          4.08972800d0
      cc( 1, 3) =         99.62212400d0
      cc( 2, 3) =         14.16873300d0
      a ( 1, 3) =          7.11000000d0
      a ( 2, 3) =          3.59679800d0
      cc( 1, 4) =         35.04512400d0
      cc( 2, 4) =          6.11125900d0
      a ( 1, 4) =          5.35000000d0
      a ( 2, 4) =          2.49182100d0
      cc( 1, 5) =        -21.09377600d0
      cc( 2, 5) =         -3.08069400d0
      a ( 1, 5) =          7.54000000d0
      a ( 2, 5) =          3.77000000d0
c
      go to 60
c
c     niobium (Nb_Stuttgart_RSC_ECP)
c
  41  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        165.17914300d0
      cc( 2, 2) =         21.99297400d0
      a ( 1, 2) =          8.90000000d0
      a ( 2, 2) =          4.43000000d0
      cc( 1, 3) =        111.79441400d0
      cc( 2, 3) =         16.63348300d0
      a ( 1, 3) =          7.77000000d0
      a ( 2, 3) =          3.96000000d0
      cc( 1, 4) =         38.11224900d0
      cc( 2, 4) =          8.03916700d0
      a ( 1, 4) =          6.05000000d0
      a ( 2, 4) =          2.84000000d0
      cc( 1, 5) =        -22.92955000d0
      cc( 2, 5) =         -3.66631000d0
      a ( 1, 5) =          8.49000000d0
      a ( 2, 5) =          4.25000000d0
c
      go to 60
c
c     molybdenum (Mo_Stuttgart_RSC_ECP)
c
  42  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        180.10310800d0
      cc( 2, 2) =         24.99722800d0
      a ( 1, 2) =          9.71459400d0
      a ( 2, 2) =          4.68050000d0
      cc( 1, 3) =        123.77275200d0
      cc( 2, 3) =         19.53022800d0
      a ( 1, 3) =          8.14213700d0
      a ( 2, 3) =          4.62598600d0
      cc( 1, 4) =         48.37502200d0
      cc( 2, 4) =          8.89205300d0
      a ( 1, 4) =          6.61841500d0
      a ( 2, 4) =          3.24875200d0
      cc( 1, 5) =        -24.80517700d0
      cc( 2, 5) =         -4.15378200d0
      a ( 1, 5) =          9.45000000d0
      a ( 2, 5) =          4.72000000d0
c
      go to 60
c
c     technetium (Tc_Stuttgart_RSC_ECP)
c
  43  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        195.15916600d0
      cc( 2, 2) =         28.09260300d0
      a ( 1, 2) =         10.42234600d0
      a ( 2, 2) =          5.03651600d0
      cc( 1, 3) =        135.28456600d0
      cc( 2, 3) =         21.80650400d0
      a ( 1, 3) =          8.95044900d0
      a ( 2, 3) =          4.85443900d0
      cc( 1, 4) =         54.32972900d0
      cc( 2, 4) =         11.15506800d0
      a ( 1, 4) =          6.94569700d0
      a ( 2, 4) =          3.97058500d0
      cc( 1, 5) =        -26.56244700d0
      cc( 2, 5) =         -4.58568100d0
      a ( 1, 5) =         10.40000000d0
      a ( 2, 5) =          5.20000000d0
c
      go to 60
c
c     ruthenium (Ru_Stuttgart_RSC_ECP)
c
  44  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        209.82297100d0
      cc( 2, 2) =         30.65472600d0
      a ( 1, 2) =         11.10526900d0
      a ( 2, 2) =          5.41474500d0
      cc( 1, 3) =        146.33618200d0
      cc( 2, 3) =         24.12787700d0
      a ( 1, 3) =          9.77127100d0
      a ( 2, 3) =          5.07399100d0
      cc( 1, 4) =         67.51589700d0
      cc( 2, 4) =          9.87010400d0
      a ( 1, 4) =          7.67142300d0
      a ( 2, 4) =          4.13656500d0
      cc( 1, 5) =        -28.34061600d0
      cc( 2, 5) =         -4.94462900d0
      a ( 1, 5) =         11.36000000d0
      a ( 2, 5) =          5.68000000d0
c
      go to 60
c
c     rhodium (Rh_Stuttgart_RSC_ECP)
c
  45  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        225.34775400d0
      cc( 2, 2) =         32.82318900d0
      a ( 1, 2) =         11.72000000d0
      a ( 2, 2) =          5.82000000d0
      cc( 1, 3) =        158.70941200d0
      cc( 2, 3) =         26.44410000d0
      a ( 1, 3) =         10.42000000d0
      a ( 2, 3) =          5.45000000d0
      cc( 1, 4) =         62.75862600d0
      cc( 2, 4) =         10.97871900d0
      a ( 1, 4) =          8.82000000d0
      a ( 2, 4) =          3.87000000d0
      cc( 1, 5) =        -30.09345600d0
      cc( 2, 5) =         -5.21848200d0
      a ( 1, 5) =         12.31000000d0
      a ( 2, 5) =          6.16000000d0
c
      go to 60
c
c     palladium (Pd_Stuttgart_RSC_ECP)
c
  46  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        240.22904000d0
      cc( 2, 2) =         35.17194300d0
      a ( 1, 2) =         12.43000000d0
      a ( 2, 2) =          6.17075900d0
      cc( 1, 3) =        170.41727600d0
      cc( 2, 3) =         28.47213300d0
      a ( 1, 3) =         11.08000000d0
      a ( 2, 3) =          5.82955400d0
      cc( 1, 4) =         69.01384500d0
      cc( 2, 4) =         11.75086200d0
      a ( 1, 4) =          9.51000000d0
      a ( 2, 4) =          4.13978100d0
      cc( 1, 5) =        -31.92955400d0
      cc( 2, 5) =         -5.39821700d0
      a ( 1, 5) =         13.27000000d0
      a ( 2, 5) =          6.63000000d0
c
      go to 60
c
c     silver (Ag_Stuttgart_RSC_ECP)
c
  47  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        255.13936500d0
      cc( 2, 2) =         36.86612200d0
      a ( 1, 2) =         13.13000000d0
      a ( 2, 2) =          6.51000000d0
      cc( 1, 3) =        182.18186900d0
      cc( 2, 3) =         30.35775100d0
      a ( 1, 3) =         11.74000000d0
      a ( 2, 3) =          6.20000000d0
      cc( 1, 4) =         73.71926100d0
      cc( 2, 4) =         12.50211700d0
      a ( 1, 4) =         10.21000000d0
      a ( 2, 4) =          4.38000000d0
      cc( 1, 5) =        -33.68992000d0
      cc( 2, 5) =         -5.53112000d0
      a ( 1, 5) =         14.22000000d0
      a ( 2, 5) =          7.11000000d0
c
      go to 60
c
c     cadmium (Cd_Stuttgart_RSC_ECP)
c
  48  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        270.00948300d0
      cc( 2, 2) =         38.76730800d0
      a ( 1, 2) =         13.83586900d0
      a ( 2, 2) =          6.85727000d0
      cc( 1, 3) =        193.82962900d0
      cc( 2, 3) =         31.89652500d0
      a ( 1, 3) =         12.40497100d0
      a ( 2, 3) =          6.56779900d0
      cc( 1, 4) =         79.19364700d0
      cc( 2, 4) =         13.23082700d0
      a ( 1, 4) =         10.89692500d0
      a ( 2, 4) =          4.64116500d0
      cc( 1, 5) =        -35.47662600d0
      cc( 2, 5) =         -5.61767700d0
      a ( 1, 5) =         15.18479600d0
      a ( 2, 5) =          7.59239800d0
c
 60   continue
      n ( 1, 1) =   2
      nt( 1)    =   1
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      nt( 2)    =   2
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      nt( 3)    =   2
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      nt( 4)    =   2
      n ( 1, 5) =   2
      n ( 2, 5) =   2
      nt( 5)    =   2
      nelcor    =  28
      lmx       =   4
c
      return
 6010 format (1x,'*** data not found for Stuttgart_RSC ecp type ',a4)
      end
**==strsc3.f
      subroutine strsc3 (zinp,lmx,nelcor,nt,n,cc,a)
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
      implicit character * 8 (z)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      dimension nt(6), n(10,6), cc(10,6), a(10,6)
c
c     this routine contains a library of ecp data taken from
c     various literature sources.
c     Stuttgart Relativistic Small Core ECP
c
      dimension ztit(17)
      data maxtyp/17/
      data ztit/'cs','ba','la',
     +          'ce','pr','nd','pm','sm','eu','gd','tb',
     +          'dy','ho','er','tm','yb','lu' /
c
c
      do 10 ityp = 1 , maxtyp
         if (zinp.eq.ztit(ityp)) go to 20
 10   continue
 57   write (iwr,6010) zinp
      call caserr2('no library available for requested ecp')
c
 20   continue
c
      go to (55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,57), ityp
c
c   cesium (Cs_Stuttgart_RSC_ECP)
c
 55   continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =         84.54773062d0
      cc( 2, 2) =         16.65417395d0
      a ( 1, 2) =          4.07975066d0
      a ( 2, 2) =          2.41740635d0
      cc( 1, 3) =        157.04905939d0
      cc( 2, 3) =         26.42330704d0
      a ( 1, 3) =          5.51408016d0
      a ( 2, 3) =          2.16031682d0
      cc( 1, 4) =         13.17275381d0
      cc( 2, 4) =          3.34283398d0
      a ( 1, 4) =          1.80741000d0
      a ( 2, 4) =          0.85818223d0
      cc( 1, 5) =        -28.88430917d0
      a ( 1, 5) =          3.12326906d0
c
      go to 30
c
c   barium (Ba_Stuttgart_RSC_ECP)
c
 56   continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        427.84581600d0
      cc( 2, 2) =        204.41753000d0
      a ( 1, 2) =          9.52698600d0
      a ( 2, 2) =          4.48751000d0
      cc( 1, 3) =        293.60586400d0
      cc( 2, 3) =        294.19331600d0
      a ( 1, 3) =          8.31593000d0
      a ( 2, 3) =          4.29221700d0
      cc( 1, 4) =        112.55040200d0
      cc( 2, 4) =        181.78262100d0
      a ( 1, 4) =          5.91610800d0
      a ( 2, 4) =          2.87484200d0
      cc( 1, 5) =        -33.47317400d0
      a ( 1, 5) =          3.58946500d0
c
30    continue
      n ( 1, 1) =   2
      nt( 1)    =   1
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      nt( 2)    =   2
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      nt( 3)    =   2
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      nt( 4)    =   2
      n ( 1, 5) =   2
      nt( 5)    =   1
      nelcor    =  46
      lmx       =   4
c
      return
c
c   cerium 
c
 58   continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        580.08345700d0
      a ( 1, 2) =         20.13782900d0
      cc( 1, 3) =        310.30283300d0
      a ( 1, 3) =         15.99848200d0
      cc( 1, 4) =        167.81394400d0
      a ( 1, 4) =         14.97418700d0
      cc( 1, 5) =        -49.39022900d0
      a ( 1, 5) =         23.40245500d0
      cc( 1, 6) =        -21.33187900d0
      a ( 1, 6) =         16.57055300d0
c
      go to 40
c
c  praseodymium
c
 59   continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        577.57312200d0
      a ( 1, 2) =         20.76627800d0
      cc( 1, 3) =        295.78584600d0
      a ( 1, 3) =         16.07844800d0
      cc( 1, 4) =        150.86705500d0
      a ( 1, 4) =         14.70508900d0
      cc( 1, 5) =        -48.73676600d0
      a ( 1, 5) =         23.37896900d0
      cc( 1, 6) =        -22.32948800d0
      a ( 1, 6) =         17.44713800d0
c
      go to 40
c
c   neodymium
c
 60   continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        574.37098000d0
      a ( 1, 2) =         21.35226700d0
      cc( 1, 3) =        280.94644000d0
      a ( 1, 3) =         16.11926500d0
      cc( 1, 4) =        138.67062700d0
      a ( 1, 4) =         14.49410300d0
      cc( 1, 5) =        -47.52266800d0
      a ( 1, 5) =         23.18386000d0
      cc( 1, 6) =        -23.34458700d0
      a ( 1, 6) =         18.34417400d0
c
      go to 40
c
c   promethium
c
 61   continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        575.39574900d0
      a ( 1, 2) =         21.94286500d0
      cc( 1, 3) =        281.70451400d0
      a ( 1, 3) =         16.55516100d0
      cc( 1, 4) =        123.52473700d0
      a ( 1, 4) =         13.96030800d0
      cc( 1, 5) =        -50.74151100d0
      a ( 1, 5) =         24.03354600d0
      cc( 1, 6) =        -24.37251000d0
      a ( 1, 6) =         19.26024500d0
c
      go to 40
c
c   samarium
c
 62   continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        572.98533200d0
      a ( 1, 2) =         22.34447100d0
      cc( 1, 3) =        272.35914500d0
      a ( 1, 3) =         16.69459000d0
      cc( 1, 4) =        115.29390000d0
      a ( 1, 4) =         13.72770500d0
      cc( 1, 5) =        -51.10839200d0
      a ( 1, 5) =         24.05909200d0
      cc( 1, 6) =        -25.42188500d0
      a ( 1, 6) =         20.19724900d0
c
      go to 40
c
c   europium
c
 63   continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        607.65933100d0
      a ( 1, 2) =         23.47138400d0
      cc( 1, 3) =        264.38547600d0
      a ( 1, 3) =         16.77247900d0
      cc( 1, 4) =        115.38137500d0
      a ( 1, 4) =         13.98134300d0
      cc( 1, 5) =        -49.40079400d0
      a ( 1, 5) =         23.96288800d0
      cc( 1, 6) =        -26.74827300d0
      a ( 1, 6) =         21.23245800d0
c
      go to 40
c
c   gadolinium
c
 64   continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        637.20086900d0
      a ( 1, 2) =         24.60215100d0
      cc( 1, 3) =        261.68960100d0
      a ( 1, 3) =         16.88925000d0
      cc( 1, 4) =        106.85653300d0
      a ( 1, 4) =         13.64335800d0
      cc( 1, 5) =        -50.68359000d0
      a ( 1, 5) =         24.12691700d0
      cc( 1, 6) =        -27.57963000d0
      a ( 1, 6) =         22.13188700d0
c
      go to 40
c
c   terbium
c
 65   continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        668.59715500d0
      a ( 1, 2) =         24.95295600d0
      cc( 1, 3) =        266.98047500d0
      a ( 1, 3) =         17.61089900d0
      cc( 1, 4) =         97.50659600d0
      a ( 1, 4) =         12.97600900d0
      cc( 1, 5) =        -52.17575700d0
      a ( 1, 5) =         24.24886900d0
      cc( 1, 6) =        -28.69426800d0
      a ( 1, 6) =         23.13067200d0
c
      go to 40
c
c   dysoprosium
c
 66   continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        705.67122100d0
      a ( 1, 2) =         26.42958600d0
      cc( 1, 3) =        254.86698900d0
      a ( 1, 3) =         17.31703400d0
      cc( 1, 4) =         95.04518700d0
      a ( 1, 4) =         12.91359900d0
      cc( 1, 5) =        -54.57409300d0
      a ( 1, 5) =         24.90787800d0
      cc( 1, 6) =        -29.82827700d0
      a ( 1, 6) =         24.14875300d0
c
      go to 40
c
c   holmium
c
 67   continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        755.70313600d0
      a ( 1, 2) =         28.39725700d0
      cc( 1, 3) =        253.55199800d0
      a ( 1, 3) =         17.43863300d0
      cc( 1, 4) =         89.63567700d0
      a ( 1, 4) =         12.43421200d0
      cc( 1, 5) =        -55.48203600d0
      a ( 1, 5) =         25.38701000d0
      cc( 1, 6) =        -30.99112500d0
      a ( 1, 6) =         25.18850100d0
c
      go to 40
c
c   erbium
c
 68   continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        800.95287600d0
      a ( 1, 2) =         29.79859200d0
      cc( 1, 3) =        262.01986900d0
      a ( 1, 3) =         18.11423700d0
      cc( 1, 4) =         80.17055200d0
      a ( 1, 4) =         11.36958700d0
      cc( 1, 5) =        -42.33628500d0
      a ( 1, 5) =         21.82123300d0
      cc( 1, 6) =        -32.18527800d0
      a ( 1, 6) =         26.25073500d0
c
      go to 40
c
c   thulium
c
 69   continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        845.51074300d0
      a ( 1, 2) =         31.14412200d0
      cc( 1, 3) =        258.58523900d0
      a ( 1, 3) =         18.09235300d0
      cc( 1, 4) =         80.72905900d0
      a ( 1, 4) =         11.46915900d0
      cc( 1, 5) =        -48.70126600d0
      a ( 1, 5) =         23.60554400d0
      cc( 1, 6) =        -33.39549600d0
      a ( 1, 6) =         27.32978100d0
c
      go to 40
c
c   ytterbium
c
 70   continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        891.01377700d0
      a ( 1, 2) =         32.42448400d0
      cc( 1, 3) =        264.03695300d0
      a ( 1, 3) =         18.65623200d0
      cc( 1, 4) =         73.92391900d0
      a ( 1, 4) =         10.49022200d0
      cc( 1, 5) =        -39.59217300d0
      a ( 1, 5) =         20.77418300d0
      cc( 1, 6) =        -34.63863800d0
      a ( 1, 6) =         28.43102800d0
c
  40  continue
      n ( 1, 1) =   2
      nt( 1)    =   1
      n ( 1, 2) =   2
      nt( 2)    =   1
      n ( 1, 3) =   2
      nt( 3)    =   1
      n ( 1, 4) =   2
      nt( 4)    =   1
      n ( 1, 5) =   2
      nt( 5)    =   1
      n ( 1, 6) =   2
      nt( 6)    =   1
c     remove g potentials from expansion
c     lmx       =   5
      lmx       =   4
      nelcor    =  28
c
      return
c
 6010 format (1x,'*** data not found for Stuttgart_RSC ecp type ',a4)
      end
**==strsc4.f
      subroutine strsc4 (zinp,lmx,nelcor,nt,n,cc,a)
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
      implicit character * 8 (z)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      dimension nt(6), n(10,6), cc(10,6), a(10,6)
c
c     this routine contains a library of ecp data taken from
c     various literature sources.
c     Stuttgart Relativistic, Small Core ECP
c
      dimension ztit(15)
      data maxtyp/15/
      data ztit/
     $ 'hf','ta','w ','re','os','ir','pt','au','hg','tl','pb',
     + 'bi','po','at','rn'/
c
c     this routine contains a library of ecp data taken from
c     various literature sources.
c     3rd transition series sjkb potentials
c     
      do 10 ityp = 1 , maxtyp
         if (zinp.eq.ztit(ityp)) go to 20
 10   continue
 81   write (iwr,6010) zinp
      call caserr2('no library available for requested ecp')
c
 20   go to(72,73,74,75,76,77,78,79,80,81,
     +      81,81,81,81,81), ityp
c
c     hafnium (Hf_Stuttgart_RSC_ECP)
c
  72  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =       1499.28471100d0
      cc( 2, 2) =         40.28210100d0
      a ( 1, 2) =         14.76995900d0
      a ( 2, 2) =          7.38497900d0
      cc( 1, 3) =        397.73300500d0
      cc( 2, 3) =         19.31640600d0
      a ( 1, 3) =          9.84949000d0
      a ( 2, 3) =          4.92474500d0
      cc( 1, 4) =        101.32980500d0
      cc( 2, 4) =          5.87343800d0
      a ( 1, 4) =          6.09675600d0
      a ( 2, 4) =          3.04837800d0
      cc( 1, 5) =         10.04672300d0
      a ( 1, 5) =          1.78577000d0
      cc( 1, 6) =         -9.55824400d0
      a ( 1, 6) =          2.63240000d0
c
      go to 30
c
c     tantalum (Ta_Stuttgart_RSC_ECP)
c
  73  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =       1345.88064700d0
      cc( 2, 2) =         36.76680600d0
      a ( 1, 2) =         14.54640800d0
      a ( 2, 2) =          7.27320400d0
      cc( 1, 3) =        378.42530100d0
      cc( 2, 3) =         22.29309100d0
      a ( 1, 3) =          9.93556500d0
      a ( 2, 3) =          4.96778200d0
      cc( 1, 4) =        104.88395600d0
      cc( 2, 4) =          8.75584800d0
      a ( 1, 4) =          6.34737700d0
      a ( 2, 4) =          3.17368800d0
      cc( 1, 5) =         12.01796100d0
      a ( 1, 5) =          2.01788100d0
      cc( 1, 6) =        -11.72893300d0
      a ( 1, 6) =          3.04033000d0
c
      go to 30
c
c     tungsten (W_Stuttgart_RSC_ECP)
c
  74  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =       1192.39588200d0
      cc( 2, 2) =         32.52293300d0
      a ( 1, 2) =         14.32285600d0
      a ( 2, 2) =          7.16142800d0
      cc( 1, 3) =        359.03196700d0
      cc( 2, 3) =         24.03038000d0
      a ( 1, 3) =         10.02164100d0
      a ( 2, 3) =          5.01082000d0
      cc( 1, 4) =        108.30134900d0
      cc( 2, 4) =         10.98252800d0
      a ( 1, 4) =          6.59799700d0
      a ( 2, 4) =          3.29899900d0
      cc( 1, 5) =         14.15257900d0
      a ( 1, 5) =          2.25888800d0
      cc( 1, 6) =        -14.05643500d0
      a ( 1, 6) =          3.46411000d0
c
      go to 30
c
c     rhenium (Re_Stuttgart_RSC_ECP)
c
  75  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =       1038.95157200d0
      cc( 2, 2) =         29.56173800d0
      a ( 1, 2) =         14.09930500d0
      a ( 2, 2) =          7.04965300d0
      cc( 1, 3) =        339.54351000d0
      cc( 2, 3) =         24.91369600d0
      a ( 1, 3) =         10.10771700d0
      a ( 2, 3) =          5.05385800d0
      cc( 1, 4) =        111.69965300d0
      cc( 2, 4) =         12.62432900d0
      a ( 1, 4) =          6.84861800d0
      a ( 2, 4) =          3.42430900d0
      cc( 1, 5) =         16.44985200d0
      a ( 1, 5) =          2.50865100d0
      cc( 1, 6) =        -16.50112000d0
      a ( 1, 6) =          3.90124500d0
c
      go to 30
c
c     osmium (Os_Stuttgart_RSC_ECP)
c
  76  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        885.40571900d0
      cc( 2, 2) =         25.96704000d0
      a ( 1, 2) =         13.87575400d0
      a ( 2, 2) =          6.93787700d0
      cc( 1, 3) =        320.08390200d0
      cc( 2, 3) =         26.14876500d0
      a ( 1, 3) =         10.19379300d0
      a ( 2, 3) =          5.09689600d0
      cc( 1, 4) =        115.04484300d0
      cc( 2, 4) =         13.62257500d0
      a ( 1, 4) =          7.09923800d0
      a ( 2, 4) =          3.54961900d0
      cc( 1, 5) =         18.90945700d0
      a ( 1, 5) =          2.76707500d0
      cc( 1, 6) =        -19.02759500d0
      a ( 1, 6) =          4.34990500d0
c 
      go to 30
c
c     iridium (Ir_Stuttgart_RSC_ECP)
c
  77  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        732.26920000d0
      cc( 2, 2) =         26.48472100d0
      a ( 1, 2) =         13.65220300d0
      a ( 2, 2) =          6.82610100d0
      cc( 1, 3) =        299.48947400d0
      cc( 2, 3) =         26.46623400d0
      a ( 1, 3) =         10.27986800d0
      a ( 2, 3) =          5.13993400d0
      cc( 1, 4) =        124.45759500d0
      cc( 2, 4) =         14.03599500d0
      a ( 1, 4) =          7.34985900d0
      a ( 2, 4) =          3.67492900d0
      cc( 1, 5) =         21.53103100d0
      a ( 1, 5) =          3.03407200d0
      cc( 1, 6) =        -21.60759700d0
      a ( 1, 6) =          4.80885700d0
c
      go to 30
c
c     platium (Pt_Stuttgart_RSC_ECP)
c
  78  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        579.22386100d0
      cc( 2, 2) =         29.66949100d0
      a ( 1, 2) =         13.42865100d0
      a ( 2, 2) =          6.71432600d0
      cc( 1, 3) =        280.86077400d0
      cc( 2, 3) =         26.74538200d0
      a ( 1, 3) =         10.36594400d0
      a ( 2, 3) =          5.18297200d0
      cc( 1, 4) =        120.39644400d0
      cc( 2, 4) =         15.81092100d0
      a ( 1, 4) =          7.60047900d0
      a ( 2, 4) =          3.80024000d0
      cc( 1, 5) =         24.31437600d0
      a ( 1, 5) =          3.30956900d0
      cc( 1, 6) =        -24.21867500d0
      a ( 1, 6) =          5.27728900d0
c
      go to 30
c
c     gold (Au_Stuttgart_RSC_ECP)
c
  79  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        426.70984000d0
      cc( 2, 2) =         35.93882400d0
      a ( 1, 2) =         13.20510000d0
      a ( 2, 2) =          6.60255000d0
      cc( 1, 3) =        261.16102300d0
      cc( 2, 3) =         26.62628400d0
      a ( 1, 3) =         10.45202000d0
      a ( 2, 3) =          5.22601000d0
      cc( 1, 4) =        124.75683100d0
      cc( 2, 4) =         15.77226000d0
      a ( 1, 4) =          7.85110000d0
      a ( 2, 4) =          3.92555000d0
      cc( 1, 5) =         30.56847500d0
      cc( 2, 5) =          5.18377400d0
      a ( 1, 5) =          4.78980000d0
      a ( 2, 5) =          2.39491000d0
c
      n ( 1, 1) =   2
      nt( 1)    =   1
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      nt( 2)    =   2
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      nt( 3)    =   2
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      nt( 4)    =   2
      n ( 1, 5) =   2
      n ( 2, 5) =   2
      nt( 5)    =   2
      nelcor    =  60
      lmx       =   4
c
      return
c
c     mercury (Hg_Stuttgart_RSC_ECP)
c
  80  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        274.53216900d0
      cc( 2, 2) =         49.18219200d0
      a ( 1, 2) =         12.98154900d0
      a ( 2, 2) =          6.49077400d0
      cc( 1, 3) =        237.39577000d0
      cc( 2, 3) =         28.12158400d0
      a ( 1, 3) =         10.53809600d0
      a ( 2, 3) =          5.26904800d0
      cc( 1, 4) =        114.25203400d0
      cc( 2, 4) =         18.49563800d0
      a ( 1, 4) =          8.10172100d0
      a ( 2, 4) =          4.05086000d0
      cc( 1, 5) =         30.36499600d0
      a ( 1, 5) =          3.88579100d0
      cc( 1, 6) =        -29.47311800d0
      a ( 1, 6) =          6.24100000d0
c
  30  n ( 1, 1) =   2
      nt( 1)    =   1
      n ( 1, 2) =   2
      n ( 2, 2) =   2
      nt( 2)    =   2
      n ( 1, 3) =   2
      n ( 2, 3) =   2
      nt( 3)    =   2
      n ( 1, 4) =   2
      n ( 2, 4) =   2
      nt( 4)    =   2
      n ( 1, 5) =   2
      nt( 5)    =   1
      n ( 1, 6) =   2
      nt( 6)    =   1
      nelcor    =  60
c     remove g potentials from expansion
c     lmx       =   5
      lmx       =   4
c
      return
c
 6010 format (1x,'*** data not found for Stuttgart_RSC ecp type ',a4)
      end
**==strsc5.f
      subroutine strsc5 (zinp,lmx,nelcor,nt,n,cc,a)
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
      implicit character * 8 (z)
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
      dimension nt(6), n(10,6), cc(10,6), a(10,6)
c
c     this routine contains a library of ecp data taken from
c     various literature sources.
c     Stuttgart Relativistic, Small Core ECP
c
      dimension ztit(17)
      data maxtyp/17/
      data ztit/
     +  'fr','ra','ac','th','pa','u ','np','pu','am','cm',
     $  'bk','cf','es','fm','md','no','lw'   /
c
c     this routine contains a library of ecp data taken from
c     various literature sources.
c     3rd transition series sjkb potentials
c     
      do 10 ityp = 1 , maxtyp
         if (zinp.eq.ztit(ityp)) go to 20
 10   continue
 87   write (iwr,6010) zinp
      call caserr2('no library available for requested ecp')
c
 20   go to(87,87,89,90,91,92,93,94,95,96,97,98,99,100,
     +     101,102,103), ityp
c
c     ac (Ac_Stuttgart_RSC_ECP)
c
  89  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        476.14692192d0
      a ( 1, 2) =         14.65187951d0
      cc( 1, 3) =        185.89552068d0
      a ( 1, 3) =          9.18034600d0
      cc( 1, 4) =        204.17603473d0
      a ( 1, 4) =          9.52316367d0
      cc( 1, 5) =        100.81825474d0
      a ( 1, 5) =          7.59420962d0
      cc( 1, 6) =        -51.61524181d0
      a ( 1, 6) =         11.02773233d0
c
      go to 30
c
c     th (Th_Stuttgart_RSC_ECP)
c
  90  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        455.05596763d0
      a ( 1, 2) =         14.77642413d0
      cc( 1, 3) =        177.69179070d0
      a ( 1, 3) =          9.08318434d0
      cc( 1, 4) =        197.23270511d0
      a ( 1, 4) =          9.65385771d0
      cc( 1, 5) =         82.90872723d0
      a ( 1, 5) =          7.38076589d0
      cc( 1, 6) =        -55.12641906d0
      a ( 1, 6) =         11.60902462d0
c
      go to 30
c
c     pa (Pa_Stuttgart_RSC_ECP)
c
  91  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        475.47573283d0
      a ( 1, 2) =         15.26917367d0
      cc( 1, 3) =        207.66485962d0
      a ( 1, 3) =          9.95773699d0
      cc( 1, 4) =        156.67411988d0
      a ( 1, 4) =          8.94314126d0
      cc( 1, 5) =         67.98223685d0
      a ( 1, 5) =          7.01046662d0
      cc( 1, 6) =        -57.63059989d0
      a ( 1, 6) =         12.20105601d0
c
      go to 30
c
c     u (U_Stuttgart_RSC_ECP)
c
  92  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        536.51662778d0
      a ( 1, 2) =         16.41403869d0
      cc( 1, 3) =        169.54492465d0
      a ( 1, 3) =          9.06055606d0
      cc( 1, 4) =        142.61559837d0
      a ( 1, 4) =          8.83183198d0
      cc( 1, 5) =         60.39307602d0
      a ( 1, 5) =          7.01851629d0
      cc( 1, 6) =        -60.12998958d0
      a ( 1, 6) =         12.80408844d0
c
      go to 30
c
c     Np (Np_Stuttgart_RSC_ECP)
c
  93  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        541.53853851d0
      a ( 1, 2) =         16.77194770d0
      cc( 1, 3) =        191.08457137d0
      a ( 1, 3) =          9.78355277d0
      cc( 1, 4) =        146.40258923d0
      a ( 1, 4) =          9.03422660d0
      cc( 1, 5) =         60.39433314d0
      a ( 1, 5) =          7.12184827d0
      cc( 1, 6) =        -62.62618508d0
      a ( 1, 6) =         13.41832253d0
c
      go to 30
c
c     Pu (Pu_Stuttgart_RSC_ECP)
c
  94  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        662.75624198d0
      a ( 1, 2) =         18.55968916d0
      cc( 1, 3) =        179.41528979d0
      a ( 1, 3) =          9.47351697d0
      cc( 1, 4) =        139.70591217d0
      a ( 1, 4) =          9.03069065d0
      cc( 1, 5) =         56.04590914d0
      a ( 1, 5) =          7.13879779d0
      cc( 1, 6) =        -65.12121611d0
      a ( 1, 6) =         14.04399655d0
c
      go to 30
c
c     Am (AM_Stuttgart_RSC_ECP)
c
  95  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        672.42690848d0
      a ( 1, 2) =         19.13009648d0
      cc( 1, 3) =        180.98363004d0
      a ( 1, 3) =          9.55914207d0
      cc( 1, 4) =        139.03756389d0
      a ( 1, 4) =          9.16126617d0
      cc( 1, 5) =         51.27355314d0
      a ( 1, 5) =          7.03575463d0
      cc( 1, 6) =        -67.61694208d0
      a ( 1, 6) =         14.68133289d0
c
      go to 30
c
c     Cm (Cm_Stuttgart_RSC_ECP)
c
  96  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        648.91058192d0
      a ( 1, 2) =         19.25034634d0
      cc( 1, 3) =        185.64663794d0
      a ( 1, 3) =          9.74566041d0
      cc( 1, 4) =        137.84026176d0
      a ( 1, 4) =          9.26702314d0
      cc( 1, 5) =         48.18099636d0
      a ( 1, 5) =          7.06371412d0
      cc( 1, 6) =        -70.11497550d0
      a ( 1, 6) =         15.33052565d0
c
      go to 30
c
c     Bk (Bk_Stuttgart_RSC_ECP)
c
  97  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        696.47512220d0
      a ( 1, 2) =         20.28537610d0
      cc( 1, 3) =        188.76026572d0
      a ( 1, 3) =          9.86839647d0
      cc( 1, 4) =        124.93583585d0
      a ( 1, 4) =          8.91655740d0
      cc( 1, 5) =         45.21650778d0
      a ( 1, 5) =          7.02500314d0
      cc( 1, 6) =        -72.61718880d0
      a ( 1, 6) =         15.99180528d0
c
      go to 30
c
c     Cf (Cf_Stuttgart_RSC_ECP)
c
  98  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =        823.74657643d0
      a ( 1, 2) =         22.27163201d0
      cc( 1, 3) =        173.07043952d0
      a ( 1, 3) =          8.80079915d0
      cc( 1, 4) =        125.56875116d0
      a ( 1, 4) =          9.31115746d0
      cc( 1, 5) =         34.95259623d0
      a ( 1, 5) =          6.62957784d0
      cc( 1, 6) =        -75.12525493d0
      a ( 1, 6) =         16.66536973d0
c
      go to 30
c
c     Es (Es_Stuttgart_RSC_ECP)
c
  99  continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =       1109.05253371d0
      a ( 1, 2) =         25.23943442d0
      cc( 1, 3) =        175.31877724d0
      a ( 1, 3) =          8.57450440d0
      cc( 1, 4) =        124.58681414d0
      a ( 1, 4) =          9.45460730d0
      cc( 1, 5) =         29.55192708d0
      a ( 1, 5) =          6.27339051d0
      cc( 1, 6) =        -77.64107478d0
      a ( 1, 6) =         17.35145618d0
c
      go to 30
c
c     Fm (Fm_Stuttgart_RSC_ECP)
c
  100 continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =       1201.22279983d0
      a ( 1, 2) =         26.49971356d0
      cc( 1, 3) =        179.58021904d0
      a ( 1, 3) =          8.45509135d0
      cc( 1, 4) =        122.51530386d0
      a ( 1, 4) =          9.52611052d0
      cc( 1, 5) =         29.49963960d0
      a ( 1, 5) =          6.48364925d0
      cc( 1, 6) =        -80.16520155d0
      a ( 1, 6) =         18.05016051d0
c
      go to 30
c
c     Md (Md_Stuttgart_RSC_ECP)
c
  101 continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =       1099.41031834d0
      a ( 1, 2) =         26.25702941d0
      cc( 1, 3) =        192.25142777d0
      a ( 1, 3) =          9.14467744d0
      cc( 1, 4) =        129.19158391d0
      a ( 1, 4) =          9.76984978d0
      cc( 1, 5) =         34.32932133d0
      a ( 1, 5) =          6.96753140d0
      cc( 1, 6) =        -82.70015285d0
      a ( 1, 6) =         18.76179172d0
c
      go to 30
c
c     No (No_Stuttgart_RSC_ECP)
c
  102 continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =       1099.41255114d0
      a ( 1, 2) =         26.78877434d0
      cc( 1, 3) =        203.31904971d0
      a ( 1, 3) =          9.67194189d0
      cc( 1, 4) =        127.62744557d0
      a ( 1, 4) =          9.65954537d0
      cc( 1, 5) =         31.66186279d0
      a ( 1, 5) =          6.60473647d0
      cc( 1, 6) =        -85.24736741d0
      a ( 1, 6) =         19.48652544d0
c
      go to 30

c     Lw (Lw_Stuttgart_RSC_ECP)
c
  103 continue
c
      cc( 1, 1) =          0.00000000d0
      a ( 1, 1) =          1.00000000d0
      cc( 1, 2) =       4933.05020435d0
      a ( 1, 2) =         41.62309750d0
      cc( 1, 3) =        208.55935232d0
      a ( 1, 3) =          9.66701176d0
      cc( 1, 4) =        105.29895096d0
      a ( 1, 4) =          8.45559601d0
      cc( 1, 5) =         36.22570362d0
      a ( 1, 5) =          7.39101231d0
      cc( 1, 6) =        -87.80798252d0
      a ( 1, 6) =         20.22453342d0
c
 30   continue
      n ( 1, 1) =   2
      nt( 1)    =   1
      n ( 1, 2) =   2
      nt( 2)    =   1
      n ( 1, 3) =   2
      nt( 3)    =   1
      n ( 1, 4) =   2
      nt( 4)    =   1
      n ( 1, 5) =   2
      nt( 5)    =   1
      n ( 1, 6) =   2
      nt( 6)    =   1
      nelcor    =  60
c     remove g potentials from expansion
c     lmx       =   5
      lmx       =   4
      return
c
 6010 format (1x,'*** data not found for Stuttgart_RSC ecp type ',a4)
      end
      subroutine ver_integb_lib(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/integb_lib.m,v $
     +     "/
      data revision /"$Revision: 5814 $"/
      data date /"$Date: 2009-01-16 10:46:55 +0100 (Fri, 16 Jan 2009) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
