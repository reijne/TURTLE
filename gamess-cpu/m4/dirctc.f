      BLOCK data cci
      common/mr_comcci/icci,iclass(3,3,3),istcc,iselcc
c...  Contracted CI & CEPA:
c...    icci   -- flag that switches contracted CI on or off.
c...    iclass -- contains the offsets for the various excitation
c...              classes in the contracted formalism. The offsets are
c...              stored as:
c...
c...                 iclass(nhole+1,npart+1,nexlvl+1)
c...
c...              where iclass is the offset for the excitation class
c...              with nhole holes, npart particles and an excitation
c...              level of nexlvl. The excitation levels may be
c...              0, 1 or 2 corresponding to reference, single or
c...              double respectively.
c...
      data icci/0/
      data iclass/ 0,-1,-1,-1,-1,-1,-1,-1,-1,
     *             1, 2,-1, 3, 4,-1,-1,-1,-1,
     *             5, 6, 7, 8, 9,10,11,12,13/
c...  -1 means that this entry is not used.
c...  (the corresponding class is not defined).
      end
      subroutine mpin(iwr)
c
      implicit real*8 (a-h,o-z), integer (i-n)
c
c...  input routine for (MR) Moller-Plesset
c...     J.H. van Lenthe, P.Pulay 1990
c...  - K.Wolinski, H.L.Sellers, P.Pulay C.P.L. 140(1987) 225
c...
c...  called from inpro
c...
c...   = directive format = (all on 1 line)
c...
c...   mp  (2/3) (model -13/-12/-11/-3/-2/-1/1/2/3/11/12/13) (print 1) 
c...       (crit (absolute/relative) 1.0 6) (level 0.0 (5 0.0)) 
c...       (maxcyc 50) (fignore 1. 8) (fock ed3 1 50)  (nodoc) (nodiag)
c...       (schmidt/schmidtp/schmidt4/lowdin/householder 1. 8)
c...
c...   2/3 : chooses between mp2 and mp3 ; mp3 is more expensive
c...         since all matrix-elements need be evaluated and
c...         all integrals are sorted, and is not implemented
c...   print : more output for a higher number (part. >2000)
c...   crit : stop-criterion for mp2-iteration  (max. residue)
c...     absolute : use stop-criterion as is.
c...     relative : use stop-criterion relative to the norm of the
c...                righthand side of linear system. (only in MPn)
c...   level : levelshift1, # iterations, levelshift2
c...           after the specified # iterations shift2 is used
c...   maxcyc : maximum # cycles in mp2-iteration-process
c...   fignore : absolute smaller FOCK-matrix elements are ignored
c...   fock : specifies dumpfile,block and section for FOCK-integrals
c...          (by default these are calculated by the program)
c...          alternatively the orbital energies corresponding to
c...          the vectors are used (if 'vectors' is specified)
c...   nodoc : tells program to treat doc as voc (for debugging)
c...   nodiag : tells the program to assume the external FOCK-
c...            integral-matrix to be non-diagonal (for debugging)
c...   schmidt : sets the convergence criterion for schmidt
c...             orthogonalisation.
c...   qr : sets the convergence criterion for the QR
c...        orthogonalisation.
c...   = only the first 4 chars of a keyword are read =
c...   = usually the following directive suffices :  mp 2
c...
c...   (nodoc is transmitted to mpini etc by iwnmp(1,1) = 1000)
c...
c...   == ROUTINES involved in MP2 : ==
c...   MPIN  (called from INPRO)
c...   MPSRCC  -  sort configuration by symmetry/occupation
c...     MPSRCS,MPSLDC,NTCSF
c...   MPGIJK,MPHBIL  - get only ints and matrix-el needed for mp2
c...     MPIJKL
c...   MPSTAR (called from program DIRECT) - symbolic/sorting
c...     MPINI
c...     MPGIJ  - fock-sort
c...       MPFOCK,SGMAT - fock-calc
c...     MPPRCC - generate contration-matrices to excited states
c...       MPPRO
c...         MPSING,GENSIN,MPEXPR  - generate/print single excitation
c...           IANNMP,ICREMP,ICHMP
c...         MPFPR  -  expand,orthogonalise and select excitations
c...           MPEXEX,MPSCHM,MPQR,MPPACT
c...     CMPFIJ,CMPFIA,CMPFAB,CMPAAB - fock-symbolic
c...       MPFIJ,MPFIA,MPFAB.MPFAAB
c...     WHERMP
c...   MPMR - control iterative MP2-solver
c...     BASMP - core-partitioning / various constants
c...     CCNTMP - contract HC0 to excited-state basis
c...       CONTMP,SIN1MP,AAGP
c...     CCXPMP - expand excited-state-function vector to
c...              configuration-state-function vector.
c...       CEXPMP,SIN1MP,AAGP
c...     MPDIAG,AADIAG - generate fock-matrix diagonal
c...     MPGEN  - (f-e0)c
c...       MPIJ,DVMPIA,SDMPIA,DDMPAB,SSMPAB
c...         SYMMPA,SYMMP,SQUAMP,MULSSP,multt
c...       MPGEND
c...         MPGIJ,AAMPAA,AAMPAB,AAMPIA
c...     DIISMP  - convergence acceleration
c...       OSINV
c...     MPSUM  - calculate various contributions to MP2-energy
c...
c...  various routines from DIRECT and the atmol-library are also used
c...  routines from DIRECT really changed :
c...     block-data,direct,inpro,pass
c...
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
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
c
      logical malt
      real*8 potnuc,totnuc,e11,eigen,h11,z11,h22,s22,h12,s12,
     *     tester,eshift,thresh,aaa,bbb,scream,gprin,z12,z22,
     *     vdstpr,qester
      integer maxcyc,moddia,mprin,icyc
      common/erg/potnuc,totnuc,e11,eigen(nd200),h11
     *,z11,h22,s22,h12,s12,tester,eshift,thresh,aaa,bbb,scream
     *,gprin,z12,z22,vdstpr(3),qester,maxcyc,moddia,mprin,malt,icyc
c
      common/diisdc/ diiscr,ntdiis
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
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
c
      logical latorb, mspin1, ispind
      integer ispacg, isping, iprwh, iprcon
      common /natorb/ ispacg,isping,latorb,mspin1,ispind,
     +                iprwh,iprcon
c
      logical oisint
      character*8 test,ztext
      character*1 cnum(9)
      parameter (nmodel=14)
      character*43 descr(nmodel)
      data descr /"       ** no projection applied **       ",
     *            "  |ref>F<ref|+|sin>F<sin|+|doub>F<doub|  ",
     *            "      |ref>F<ref| + |doub>F<doub|        ",
     *            "   |ref>F<ref| + |sin+doub>F<sin+doub|   ",
     *            "  |ref>F<ref|+|sin>F<sin|+|rest>F<rest|  ",
     *            "      |ref>F<ref| + |rest>F<rest|        ",
     *            "  |ref>F<ref| + |doub+more>F<doub+more|  ",
     *            "|ref>F<ref|+|sin^>F<sin^|+|doub^>F<doub^|",
     *            "      |ref>F<ref| + |doub^>F<doub^|      ",
     *            " |ref>F<ref| + |sin^+doub^>F<sin^+doub^| ",
     *            " ",
     *            " ",
     *            " ",
     *            " "/
c
      data cnum/'1','2','3','4','5','6','7','8','9'/
c
      data ncall/0/
      save ncall
c
      ncall = ncall + 1
c
c...  if this subroutine has been called once before then 
c...  skip the initialisations and start reading subdirectives
c
      if (ncall.ne.1) goto 50
c
      mp = 2
      isec_mp = 0
      nummp = 0
      maxcyc = 0
      critmp = 1.0d-6
      maxcmp = 50
      iprimp = 1
      mpmod = -1
      iwnmp(1,1) = 0
      indmp(6) = 0
      plevel(1) = .0d0
      plevel(2) = plevel(1)
      fignor = 1.0d-8
      cortho = 1.0d-8
c...  default orthogonalisation method is Householder QR
      iortho = 5
      ifdiaf = 1
c...  disable configuration printing (use prconf after to switch on)
      iprcon = -1
c
c...  set default diis (start straight away)
c
      if (diiscr.eq.0.0d0) then
         diiscr = 100.0d0
         ntdiis = -1
      end if
c
50    continue
c
      call inpa(test)
      if (test.eq.'2') then
         mp = 2
      else if (test.eq.'3') then
         mp = 3
      else if (test.eq.'n') then
         call inpi(mp)
      else
         jrec = jrec -1
      end if
c
10    call inpa(test)
c
      if (test(1:4).eq.'mode') then
         call inpa(ztext)
         call intval(ztext,mpmod,oisint)
         if (.not.oisint) then
            if (ztext.eq.'ruttink') then
               mpmod = -1
            else if (ztext.eq.'pulay') then
               mpmod = 1
            else if (ztext.eq.'andersso') then
               mpmod = 11
            else
               write(*,*)'No such model ',ztext,' !!!'
               write(*,*)'Know only: Ruttink, Pulay and Andersson'
               call caserr('No such model')
            endif
         endif
c...     save sign of mpmod before the checks.
         i = isign(1,mpmod)
         mpmod = abs(mpmod)
         if (mpmod.lt.0.or.mpmod.gt.nmodel-1) call caserr('wrong model')
         if ((mpmod.ne. 1).and.(mpmod.ne. 2).and.(mpmod.ne. 3).and.
     *       (mpmod.ne.11).and.(mpmod.ne.12).and.(mpmod.ne.13)) 
     *       call caserr('no such model')
c...     restore the sign.
         mpmod = isign(mpmod,i)
      else if (test(1:4).eq.'prin') then
         call inpi(iprimp)
      else if (test(1:4).eq.'crit') then
         call inpa(test)
         if (test(1:4).eq.'abso') then
           call inpf(bilbo)
           call inpi(j)
           critmp = bilbo*(0.1d0**j)
         else if (test(1:4).eq.'rela') then
           call inpf(bilbo)
           call inpi(j)
           critmp = - bilbo*(0.1d0**j)
         else
           jrec = jrec - 1
           call inpf(bilbo)
           call inpi(j)
           critmp = bilbo*(0.1d0**j)
         endif
         if(mp.le.3) then
           critmp = dabs(critmp)
         endif
      else if (test(1:4).eq.'fign') then
         call inpf(bilbo)
         call inpi(j)
         fignor = bilbo*(0.1d0**j)
      else if (test(1:8).eq.'schmidtp') then
         iortho = 2
         call inpf(bilbo)
         call inpi(j)
         cortho = bilbo*(0.1d0**j)
      else if (test(1:8).eq.'schmidt4') then
         iortho = 3
         call inpf(bilbo)
         call inpi(j)
         cortho = bilbo*(0.1d0**j)
      else if (test(1:7).eq.'schmidt') then
         iortho = 1
         call inpf(bilbo)
         call inpi(j)
         cortho = bilbo*(0.1d0**j)
      else if (test(1:6).eq.'lowdin') then
         iortho = 4
         call inpf(bilbo)
         call inpi(j)
         cortho = bilbo*(0.1d0**j)
      else if (test(1:5).eq.'house') then
         iortho = 5
         call inpf(bilbo)
         call inpi(j)
         cortho = bilbo*(0.1d0**j)
      else if (test(1:4).eq.'maxc') then
         call inpi(maxcmp)
      else if (test(1:4).eq.'leve') then
         call inpf(plevel(1))
         call inpa(test)
         jrec = jrec - 1
         if (locatc(cnum,9,test(1:1)).ne.0) then
            call inpi(levelp)
            call inpf(plevel(2))
         else
            levelp = 0
            plevel(2) = plevel(1)
         end if
      else if (test(1:4).eq.'nodo') then
         iwnmp(1,1) = 1000
         print *,' ***debug*** no doublys are recognised'
      else if (test(1:4).eq.'nodi') then
         ifdiaf = 0
         print *,' FOCK matrix assumed **non**-diagonal'
      else if (test(1:4).eq.'noor') then
         print *,' no noorth option !!'
      else if (test(1:4).eq.'fock') then
         call inpa(test)
         if (test(1:4).eq.'vect') then
            nummp = -999
         else
           nummp = locatc(yed,maxlfn,test(1:4))
           if (nummp.eq.0) call caserr('fock-subdirective unrecognised')
           call inpi(iblk_mp)
           call inpi(isec_mp)
         end if
      else if (test(1:1).eq.' ') then
         go to 20
      else
         call caserr(' subdirective not recognised ')
      end if
c
      go to 10
c
c...  print input
c
20    write(iwr,1) mp
      if (abs(mpmod).ge.1.and.abs(mpmod).le.3) then
         write(iwr,2) mpmod,descr(abs(mpmod)+1)
      else if (abs(mpmod).ge.11.and.abs(mpmod).le.13) then
         write(iwr,2) mpmod,descr(abs(mpmod)-3)
      endif
      if (mpmod.lt.0) then
         write(iwr,69)
      endif
      write(iwr,3) iprimp
      write(iwr,4) maxcmp,critmp
      if (plevel(1).ne.0.0d0.or.plevel(2).ne.0.0d0)
     *   write(iwr,6) plevel(1),levelp,plevel(2)
      write(iwr,8) fignor
      if (iortho.eq.1) then
        write(iwr,9) 'Schmidt ...................',cortho
      else if (iortho.eq.2) then
        write(iwr,9) 'Pivoted Repeated Schmidt ..',cortho
      else if (iortho.eq.3) then
        write(iwr,9) 'Quadruple Precision Schmidt',cortho
      else if (iortho.eq.4) then
        write(iwr,9) 'Lowdin ....................',cortho
      else if (iortho.eq.5) then
        write(iwr,9) 'Householder QR ............',cortho
      endif
      if (nummp.gt.0) then
         write(iwr,5) yed(nummp),iblk_mp,isec_mp
      else if (nummp.eq.-999) then
         write(iwr,601)
      else
         write(iwr,7)
      end if
c....
1     format(//,
     1 '  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$',/,
     2 '  $ Multi Reference MP ',i1,' option requested CONTRACT      $')
2     format(
     1 '  $ model',i3,1x,a43,                                       '$')
3     format(
     1 '  $    print requested  ',i6,     '                          $')
4     format(
     1 '  $    # cycles ',i6,'  crit ', e9.3,     '                  $')
6     format(
     1 '  $    level-shift ',e8.3,' ; ',i6,'  ,  ',e8.3,     '       $')
8     format(
     1 '  $    ignore FOCK-integrals < ',e8.3,     '                 $')
9     format(
     1 '  $    ',a27, ' upto ',e8.3,                        '        $')
5     format(
     1 '  $    FOCK-matrix from ',a3,i5,i4,'                        $',
     1/,'  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
7     format(
     1 '  $    FOCK-matrix will be calculated                   $',/,
     1 '  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
601   format(
     1 '  $    FOCK-matrix from vectors section                 $',/,
     1 '  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
69    format(
     1 '  $ H0 block diagonal in excitation classes !!!         $')
c
      return
      end
      subroutine mpmask(conf,nref,nw,nword)
      implicit real*8 (a-h,o-z)
c
c...  create the mask that only contains the doubly occupied orbitals
c...  (the doc orbitals). This mask is needed in mpsrcc and should be
c...  generated before the configuration selection. The latter is
c...  because reference configurations for 'legal' excited
c...  configurations may have the wrong (spin)symmetry and may be
c...  discarded in the selection step.
c
      integer conf,dmsk
      dimension conf(nword,*),dmsk(3)
c
      integer mms,ms,nd32,n32m,mmms
      common /cdcryz/ mms,ms,nd32,n32m,mmms
c
c
       call setsto(3,0,dmsk)
c     call vclr(dmsk,1,3)
      do 20 i=1,nw
         dmsk(i) = n32m
         do 10 j=1,nref
10       dmsk(i) = iand(dmsk(i),conf(i,j))
20    dmsk(i) = ior(dmsk(i),ishft(dmsk(i),-1))
      call smpmsk(dmsk)
      end
      subroutine mpsrcc(conf,resul)
c
c...  sort configurations according to the occupation of
c...  the (reference) doubly occupied orbitals and
c...  sort out addressing in conf array's
c...
c...  the  ncmph  array's contain info for the sorted subdivisions :
c...  (also results from the excitation generator)
c...   ncmphn : vacuum ; ncmph1 : n-1 ; ncmph2 : n-2 ; ncmphd : n-2-aa
c...   1st index : # conf's (1,...) / #  csf's (2,...) / # sets (3....)
c...               # singlets (4,..) / # triplets (5,..)
c...               # refs (6,..) / # singles (7,..) / # doubles (8,..)
c...               # excited states (9,..) /total # exc. confs (10,..)
c...   2nd index : # holes in doubly occ. :
c...                0 holes (..,1,..)
c...                1 hole (..,2 - nirr+1,..)
c...                2 holes in 1 orbital (..,2+nirr,..)
c...                2 holes in 2 orbitals (..,2+nirr+1 - 2+nirr+nirr,..)
c...   3d  index : symmetry of conf (only for n-1 and n-2)
c...   clearly not all info is relevant to each type
c...   ncmphd is almost identical to ncmph2 (only symmetry 1 though)
c..
c...
c...  sorting will be done within symmetry within vacuum/n-1/n-2 :
c...   order : no holes
c...           single holes / doc symmetry (= in doc order)
c...           double holes (all  1 symmetry)
c...           two single holes / doc symmetry
c...
c
      implicit real*8 (a-h,o-z)
      integer conf,resul,dmsk,rmsk
      dimension conf(5,*),resul(5,*),dmsk(3)
c
c
      integer nint, next, ninnex, nwm1, nw, nwp1, nintwo
      integer nint21, nintp1
      common /orb/ nint,next,ninnex,nwm1,nw,nwp1,nintwo,
     +             nint21,nintp1
c
c
      integer nref, nnst, nmin1, nmin2, mspin, nelec
      integer nstrt1, nlast1, nstrt2, nlast2, nmin12, nminv2, nminv1
      integer nlst21, nlastx
      common /ccntl/ nref,nnst,nmin1,nmin2,mspin,nelec,
     +               nstrt1,nlast1,nstrt2,nlast2,nmin12,nminv2,nminv1,
     +               nlst21,nlastx
c
c
      integer nirr,mult,ifrep,nirr8,jfrep
      common /table/ nirr,mult(8,8),ifrep,nirr8,jfrep
c
c
      integer nconst, ncond, maxpa, maxpat, maxsin, maxepa
      integer mxdvpa, mxddpa, nspips, nspipt, mxsspa, mxttpa
      integer mxstpa, nsitot, ntrtot
      integer mx12pa, mxvvpa, mxv2pa, maxpsq, max12
      common /auxh/ nconst(8),ncond(8),maxpa(4),maxpat,maxsin(3),
     +              maxepa(4),mxdvpa,mxddpa,nspips(8),nspipt(8),
     +              mxsspa,mxttpa,mxstpa,nsitot(8),ntrtot(8),
     +              mx12pa,mxvvpa,mxv2pa,maxpsq,max12
c
c
      integer n24, magic, n32, n64, n8, n16, n48
      common /cdcryi/ n24,magic,n32,n64,n8,n16,n48
c
c
      integer mms,ms,nd32,n32m,mmms
      common /cdcryz/ mms,ms,nd32,n32m,mmms
c
c
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
      common/mr_mp2l/ ivoc(62),nvoc,idoc(62),ndoc,kvcc(4),nvcc(4),
     *                ncmphn(10,18),ncmph1(10,18,8),ncmph2(10,18,8),
     *                ncmphd(10,18)
c
      dimension iocc(62)
c
c       if nodo-directive => doc treated as voc
c
      call gmpmsk(dmsk)
      if (iwnmp(1,1).eq.1000) then
       call vclr(dmsk,1,3)
      endif

c
c...  make list of double /  variable occupied orbitals
c
      call upack(dmsk,nint,iocc)
      ndoc = 0
      nvoc = 0
      do 30 i=1,nint
         if (iocc(i).eq.2) then
            ndoc = ndoc + 1
            idoc(ndoc) = i
         else
            nvoc = nvoc + 1
            ivoc(nvoc) = i
         end if
30     continue
c
      nresul = 0
c
c...  ==  vacuum ==
c
      call mpsrcs(conf(1,1),nnst,dmsk,resul,nresul,ncmphn)
c
c...  === doublet === //per symmetry\\
c
      kc = nnst + 1
      do 40 ir=1,nirr
         call mpsrcs(conf(1,kc),ncond(ir),dmsk,resul,nresul,
     *               ncmph1(1,1,ir))
         kc = kc + ncond(ir)
40    continue
c
c...  === N-2 === //per symmetry\\
c
      do 50 ir=1,nirr
         call mpsrcs(conf(1,kc),nconst(ir),dmsk,resul,nresul,
     *               ncmph2(1,1,ir))
         kc = kc + nconst(ir)
50    continue
c
      call icopy(10*18,ncmph2,1,ncmphd,1)
c
c...  get n-2 (aa) right (no triplets!!)
c
      do 60 j=1,18
            ncmphd(2,j) = ncmphd(4,j)
            ncmphd(5,j) = 0
60    continue
c
      if (kc.ne.nlast2+1) call caserr('count1 error in mpsorc')
      if (nresul.ne.nlast2) call caserr('count2 error in mpsorc')
c
c...  move everything back
c
      call dcopy(nresul*5,resul,1,conf,1)
c
      return
      end
      subroutine mpsrcs(conf,nconf,dmsk,resul,nresul,nconh)
c
c...  sort configurations according to the occupation of
c...  the (reference) doubly occupied orbitals and
c...  sort out addressing in conf array's
c...  all for the config's 1 .. nconf
c
c...  **c alled by mpsrcc **
c...
c...   order : no holes
c...           single holes / doc symmetry
c...           double holes (all  1 symmetry)
c...           two single holes / doc symmetry
c...   dimensions of config's ,csf's and sets  are returned
c...   (see mpsrcc)
c
      implicit real*8 (a-h,o-z)
      integer conf,dmsk,resul
      dimension conf(5,nconf),resul(5,*),dmsk(3)
      dimension nconh(10,18)
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
      integer nirr,mult,ifrep,nirr8,jfrep
      common /table/ nirr,mult(8,8),ifrep,nirr8,jfrep
c
c
      integer iro, iky, niky, m1e, m2e, m1ia, m2coul
      integer m2exch, m2fock, m2psup, m2qsup
      integer mbase, nbase, npairi, norbi, mbasi, mnbasi
      integer iap, mbfock, mbexch, mbcoul, mbpsup, mbqsup
      integer norbas, m2exab, m2em1e
      common /helpr/ iro(nd200),iky(nd200),niky(8),m1e,m2e,m1ia,m2coul,
     +   m2exch,m2fock,m2psup,m2qsup,
     +   mbase(288),nbase(8),npairi(36),norbi(8),mbasi(8),mnbasi(8),
     +   iap(62),mbfock(62),mbexch(36),mbcoul(36),mbpsup(36),mbqsup(36),
     +   norbas(nd200),m2exab,m2em1e
c
c
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
      common/mr_mp2l/ ivoc(62),nvoc,idoc(62),ndoc,kvcc(4),nvcc(4),
     *                ncmphn(10,18),ncmph1(10,18,8),ncmph2(10,18,8),
     *                ncmphd(10,18)
c
      call izero(9*18,nconh,1)
c..    set # excited states to -1 (cf. mppro)
      do 1 i=1,18
       do 1 j=6,9
1     nconh(j,i) = -1
c
c..   no holes
c
      nn = nresul
      call mpsldc(conf,nconf,0,0,dmsk,resul,nresul)
      nset = nresul - nn
      nconh(1,1) = nset
      nconh(2,1) = ntcsf(resul(1,nn+1),nset,nconh(4,1),nconh(5,1))
      if (nset.ne.0) nconh(3,1) = 1

c
c...  single holes
c
      do 10 i=1,ndoc
         nn = nresul
         call mpsldc(conf,nconf,0,idoc(i),dmsk,resul,nresul)
         ir = iro(idoc(i))
         nset = nresul - nn
         if (nconh(1,ir+1).eq.0) then
             nconh(1,ir+1) = nset
             nconh(2,ir+1) = ntcsf(resul(1,nn+1),nset,nconh(4,ir+1),
     1                                                nconh(5,ir+1))
         end if
         if (nset.ne.0) nconh(3,ir+1) = nconh(3,ir+1) + 1
         if (nconh(1,ir+1).ne.nset) call caserr('1 hole mpsrcs')
10    continue
c
c...  double holes in one doc
c
      do 20 i=1,ndoc
         nn = nresul
         call mpsldc(conf,nconf,idoc(i),idoc(i),dmsk,resul,nresul)
         nset = nresul - nn
         if (i.eq.1) then
            nconh(1,2+nirr) = nset
            nconh(2,2+nirr) = ntcsf(resul(1,nn+1),nset,nconh(4,2+nirr),
     1                                                 nconh(5,2+nirr))
         end if
         if (nset.ne.0) nconh(3,2+nirr) = nconh(3,2+nirr) + 1
         if (nconh(1,2+nirr).ne.nset) call caserr('d2 hole mpsrcs')
20    continue
c
c...  double holes in different docs
c...
c
      do 50 ir=1,nirr
       do 40 i=2,ndoc
        do 30 j=1,i-1
         if (mult(iro(idoc(i)),iro(idoc(j))).eq.ir) then
           nn = nresul
           call mpsldc(conf,nconf,idoc(i),idoc(j),dmsk,resul,nresul)
           nset = nresul - nn
           if (nconh(1,2+nirr+ir).eq.0) then
              nconh(1,2+nirr+ir) = nset
              nconh(2,2+nirr+ir) = ntcsf(resul(1,nn+1),nset,
     1                        nconh(4,2+nirr+ir),nconh(5,2+nirr+ir))
           end if
           if (nset.ne.0) nconh(3,2+nirr+ir) = nconh(3,2+nirr+ir) + 1
           if (nconh(1,2+nirr+ir).ne.nset) call caserr('11 hole mpsrcs')
         end if
30      continue
40     continue
50    continue
c
      return
      end
      subroutine mpsldc(conf,nconf,iorb,jorb,dmsk,resul,nresul)
c
      implicit real*8 (a-h,o-z)
c
cap*  beware of silly action on compares
c
      integer conf,resul,dmsk,holem,temp
c
      dimension conf(5,nconf),resul(5,*),dmsk(3)
c
c...  select from the set in conf(1..nconf) the ones having
c...  holes in iorb and jorb (with respect to dmask = doubly-mask)
c...  write them to resul and return # written in nresul (accumulative)
c...  iorb/jorb may be 0 indicating no hole in doubly space
c
c
      integer nint, next, ninnex, nwm1, nw, nwp1, nintwo
      integer nint21, nintp1
      common /orb/ nint,next,ninnex,nwm1,nw,nwp1,nintwo,
     +             nint21,nintp1
c
c
      dimension holem(3),temp(3)
c
c..   set up detection-mask for this set in holem
c
      if (iorb.ne.0) then
         if (iannmp(dmsk,iorb,temp).eq.0) call caserr(' mpsorc 1 ')
      else
         call dcopy(nwm1,dmsk,1,temp,1)
      end if
      if (jorb.ne.0) then
         if (iannmp(temp,jorb,holem).eq.0) call caserr(' mpsorc 2 ')
      else
         call dcopy(nwm1,temp,1,holem,1)
      end if
c
      do 10 i=1,nconf
         do 5 j=1,nwm1
           temp(j) = iand(dmsk(j),conf(j,i))
           if (temp(j).ne.holem(j)) go to 10
cjvl          temp(j) = dand(dmsk(j),conf(j,i))
cjvl          if (locate(temp(j),1,1,holem(j)).eq.0) go to 10
5        continue
c
c..      found one
c
         nresul = nresul + 1
         call dcopy(5,conf(1,i),1,resul(1,nresul),1)
c
10    continue
c
      return
      end
      integer function ntcsf(conf,nconf,nsing,ntrip)
c
      implicit real*8 (a-h,o-z)
c
      integer conf
      dimension conf(5,nconf)
c
c...  count te total number of (internal csf's associated with the
c...  configurations conf(1..nconf)
c
      nsing = 0
      ntrip = 0
      do 10 i=1,nconf
         call upack2(conf(3,i),ncolt,ncols)
         nsing = nsing + ncols
         ntrip = ntrip + ncolt
10    continue
c
      ntcsf = nsing + ntrip
c
      return
      end
      subroutine mpgijk(cr,icr)
c
c...   modified version of gijkl for mp2
c...   we do'nt need all integrals
c...   ** we tried to adapt to GAMESS **
c
      implicit real*8  (a-h,o-z),integer  (i-n)
      external fget
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
      dimension cr(*),icr(*)
      common/craypk/i205(340),j205(340),k205(340),l205(340)
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
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
c
      integer nref, nnst, nmin1, nmin2, mspin, nelec
      integer nstrt1, nlast1, nstrt2, nlast2, nmin12, nminv2, nminv1
      integer nlst21, nlastx
      common /ccntl/ nref,nnst,nmin1,nmin2,mspin,nelec,
     +               nstrt1,nlast1,nstrt2,nlast2,nmin12,nminv2,nminv1,
     +               nlst21,nlastx
c
c
       integer mapei, mapie
       logical modei
       common/mapp/mapei(nd200),modei,mapie(nd200)
c
      common/stoctl/khbill,khbm1e,khbm2e,ijklno
c
      integer iro, iky, niky, m1e, m2e, m1ia, m2coul
      integer m2exch, m2fock, m2psup, m2qsup
      integer mbase, nbase, npairi, norbi, mbasi, mnbasi
      integer iap, mbfock, mbexch, mbcoul, mbpsup, mbqsup
      integer norbas, m2exab, m2em1e
      common /helpr/ iro(nd200),iky(nd200),niky(8),m1e,m2e,m1ia,m2coul,
     +   m2exch,m2fock,m2psup,m2qsup,
     +   mbase(288),nbase(8),npairi(36),norbi(8),mbasi(8),mnbasi(8),
     +   iap(62),mbfock(62),mbexch(36),mbcoul(36),mbpsup(36),mbqsup(36),
     +   norbas(nd200),m2exab,m2em1e
c
c
      integer nirr,mult,ifrep,nirr8,jfrep
      common /table/ nirr,mult(8,8),ifrep,nirr8,jfrep
c
c
      integer nint, next, ninnex, nwm1, nw, nwp1, nintwo
      integer nint21, nintp1
      common /orb/ nint,next,ninnex,nwm1,nw,nwp1,nintwo,
     +             nint21,nintp1
c
      common/blkin/value(510),mmword,mspace,mult8(8,8)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common/bufb /nwbnwb,lnklnk,buf2(1)
      common/stak/btri,mlowww,nscrp,iblock
c
      real*8 valp
      integer link, nwbuf, lnumb, iblkp, indxi, indxj, indxk
      integer indxl, istakp
      common /presrt/ link(78),nwbuf(78),lnumb(78),ibasen(78),
     +                iblkp,indxi,indxj,indxk,indxl,istakp,valp
c
c
      integer nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl
      integer nsz510,nsz341,nsz342,nsz680
      integer nsstat,nslen,nsmax,nsort, nstid
      logical oasync
      common/blksiz/nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     * nsz510,nsz341,nsz342,nsz680,nsstat,nslen,nsmax,nsort,
     * nstid(2),oasync
c
      common/scra/vscr(3680),istscr(3680)
     *,ijscr(3680),klscr(3680)
      common/mapper/iky2(5,maxorb),i4096(1)
      common/junk/icyb(3412),jcyb(3412),kcyb(3412),lcyb(3412)
c
      integer isec3, nsect2, notap2, iblk2, lblk2
      common /file12/ isec3,nsect2,notap2(20),iblk2(20),lblk2(20)
c
      common/lsort/occ(maxorb),potn,core,
     +    ilsort(4*maxorb+7),itype(nd200)
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
      common/erg/potnuc,totnuc,energy
c
      integer nuse, ngot
      common /corctl/ nuse, ngot
c
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
c
      ihint(i,j) = iky(norbas(max(i,j))) + min(i,j)
     +             + nbase(iro(max(i,j)))
c
      joop = 0
      nin = 0
      nscrp = 0
      ivmax = 0
c...
c... 1-e ij in core
c... 1-e ia bucket 1
c... 1-e ab bucket 2  (not needed)
c... 2-e ijkl bucket 8
c... 2-e ijka bucket 3
c... 2-e ijab bucket 4   not needed
c... 2-e iiab bucket 5   not needed
c... 2-e iajb bucket 6
c... 2-e iaib bucket 7
c... 2-e iabc bucket 8+nirr+i  not needed
c... 2-e abcd bucket 8+iksym   not needed
c...
      iblock = iblkp  
      call cpuwal(top,bot)
      khbm1e=khbill+m1e
      khbm2e=khbm1e+m2e
      nirr8=nirr+8
      nint9=nint+nirr8
      nuse = 0
      nusei = 0
      do 20 loop = 1 , nint9
        lnumb(loop) = nuse
        nuse = nuse + nsz340
        ibasen(loop) = nusei
        nusei = nusei + nsz680
 20   continue
      do 30 loop = 1 , ijklno
        mbase(loop) = mbase(loop) + khbm1e
 30   continue
      do 40 loop = 1 , nirr
        joop = iky(norbi(loop)+norbe(loop)+1) + joop
        nbase(loop) = nbase(loop) + khbill
 40   continue
      nuse = max(khbm2e,nint9*nsz340*2+khbm1e)
      call corfai('mpsort')
      nav = lenwrd()
      i10 = khbm1e + 1
      i20 = i10 + nint9*nsz340
      i20 = (i20-1)*nav + 1
      call setsto(nint9*nsz340*2,0,cr(i10))
      call setsto(nd200,-1000,itype)
      call setsto(ninnex,1,itype)
      call setsto(nint,0,itype)
      call setsto(1360,0,i205)   
      call setsto(3680,0,klscr)
c...
c... sort 1-elec integrals
c...
      call secget(isec3,1004,jblkk)
      call rdedx(occ,mach(15),jblkk,idaf)    
      if(.not.oprint(28)) write (iwr,6010) top , isec3 , ibl3d ,
     +     yed(idaf)
 6010 format (/1x,80('-')
     +        //' **** transformed integral pre-sort called at',f9.2,
     +        ' secs'//
     +        ' transformed 1-electron integrals restored from section',
     +        i4,' of dumpfile starting at block',i8,' of ',a4/,
     *' *** restore only MP2 integrals *** ')
      potnuc=core
      totnuc=potn
      jblkk=jblkk+ lensec(mach(15))
      call search(jblkk,idaf)
      valmax = 0.0d0
      call fget(value,kword,idaf)
      nwbnwb=nsz340
      do 101 iii=1,nd200
      iiii=mapei(iii)
      do 101 jjj=1,iii
       j=mapei(jjj)
       i=iiii
       if(i.ge.j)goto 8332
       i=j
       j=iiii
8332   nin=nin+1
      if(nin.le.kword)goto 102
      call fget(value,kword,idaf)
            if(nscrp.ge.3170)
     *       call stack(cr(i10),icr(i20))
      nin=1
102   val=value(nin)
      if(iro(j).ne.iro(i))goto 511
      istack=itype(i)+itype(j)
      if(istack)101,404,402
c... (i/i)
404   cr(ihint(i,j))=val
      goto 10101
c...   symmetry forbidden integral
511   if (dabs(val).le.dabs(valmax)) go to 101
      valmax = val
      ivmax = i
      jvmax = j
      go to 101
c... (e/i) or (e/e)
402   nscrp=nscrp+1
      vscr(nscrp)=val
      istscr(nscrp)=istack
      ijscr(nscrp) = j + i4096(i)
10101  joop=joop-1
       if(joop.eq.0) go to 10102
101   continue
c...
c... sort 2-elec integrals
c...
10102   if (ivmax.ne.0) write(iwr,512) ivmax,jvmax,valmax
512     format(/' *** symmetry check / largest forbidden h-integral :',
     1         2i4,e12.5)
      do loop=1,nirr
        do ioopn =1,nirr
          mult8(ioopn,loop) = mult(ioopn,loop) + 8
        end do
      end do
      do 20113 ioopn=1,nsect2
      ibl=iblk2(ioopn)
      lbl=ibl-lblk2(ioopn)
      if(lbl)66889,20113,66889
66889 iunit=notap2(ioopn)
      call search(ibl,iunit)
66    call fget(value,joop,iunit)
        if(nscrp.ge.3001)
     *      call stack(cr(i10),icr(i20))
      if (joop.eq.0) go to 20113
c... process block
      if(modei) then
         call unpack(value(num2e+1),lab816,i205,numlab)
      else
         call upak8w(value(num2e+1),i205,mapei)
      end if
c
            int4 = 1  
      do 33 iwword=1,mmword
              j = i205(int4  )
              i = i205(int4+1)
              l = i205(int4+2)
              k = i205(int4+3)
      val=value(iwword)
      iroi=iro(i)
      iroki=mult8(iro(k),iroi)
      irol=iro(l)
      if(mult8(iro(j),irol).ne.iroki)goto 33
      iti=itype(i)+itype(j)+itype(k)+itype(l)
      if(iti)33,34,35
c... (ii/ii)
34     istack=8
       goto 333
35    goto (501,502,503,504),iti
c... (ai/ii)
501   if(nmin1)5011,33,5011
5011  istack=3
      goto 333
c... (aa/ii) or (ai/ai)
502   if(j.le.nint)goto 5502
c... (aa/ij)
      goto 33
c... (aa/ii)
c... (aj/ai)
5502  if(j.eq.l)goto 6502
      istack=6
      goto 333
c... (ai/ai)
6502  istack=7
      goto 333
c... (aa/ai)-standard order   or   (ai/aa)
503   go to 33
c... (aa/aa)
504   go to 33
c.......................
333   nscrp=nscrp+1
      vscr(nscrp)=val
      istscr(nscrp)=istack
      ijscr(nscrp) = j + i4096(i)
      klscr(nscrp) = l + i4096(k)
 33          int4 = int4 + 4
      lbl=lbl+1
      if(lbl)66,20113,66
20113 continue
      if (nscrp.ne.0)
     * call stack(cr(i10),icr(i20))
c...
c... clear up bucketing
c...
      lbl=nsz340+1
      do 7777 joop=1,nint9
      nwb=nwbuf(joop)
      if(nwb)7776,7777,7776
7776  call stopbk
      nwbnwb=nwb
      iiii=lnumb(joop)
      lnklnk=link(joop)
      jjjj = ibasen(joop)
      call dcopy(nwb,cr(i10+iiii),1,buf2,1)
      call pack(buf2(nsz341),lab1632,icr(i20+jjjj),nsz680)
      call sttout
      link(joop)=iblock
      iblock = iblock + nsz
7777  continue
      if(.not.oprint(28)) write (iwr,4100) iblock
      call stopbk
      iblkp = iblock                             
4100  format(//' sort file',i8,' blocks long')
c...
c... get (ij/kl) integrals into core
c...
      nuse = khbm2e
      call corfai('mpgijk')
      call vclr(cr(khbm1e+1),1,m2e)
      call setsto(nsz340,0,icyb)
      call setsto(nsz340,0,jcyb)
      call setsto(nsz340,0,kcyb)
      call setsto(nsz340,0,lcyb)
      iblok=link(8)
 200  if (iblok.ge.0) then
        call rdbak(iblok)
        call stopbk
        call upak8s(buf2,icyb)
        if (odebug(40)) call wrtsor('ijkl  ',icyb,jcyb,kcyb,lcyb,buf2,
     +                              nwbnwb)
        do 210 loop = 1 , nwbnwb
          cr(irint(icyb(loop),jcyb(loop),kcyb(loop),lcyb(loop)))
     +      = buf2(loop)
 210    continue
        iblok = lnklnk
        go to 200
      else
        return
      end if
c
      end
      subroutine mphbil(cr,icr)
cmp
cmp..  special version for mp2  (2-only !!)
cmp..  only calculate interactions with reference configs
cmp
c... model=1 end of coupling coefficients
c... model=2 ijkl vacuum/vacuum
c... model=5 ijka+fock(ia) doublet/vacuum
c... model=11 iajb p+q supermatrix  n-2/vacuum
c... model=12 iaib n-2/vacuum
c
      implicit real*8 (a-h,o-z)
      integer  icr
      dimension cr(*),icr(*)
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
      integer nuse, ngot
      common /corctl/ nuse, ngot
c
c
      integer nref, nnst, nmin1, nmin2, mspin, nelec
      integer nstrt1, nlast1, nstrt2, nlast2, nmin12, nminv2, nminv1
      integer nlst21, nlastx
      common /ccntl/ nref,nnst,nmin1,nmin2,mspin,nelec,
     +               nstrt1,nlast1,nstrt2,nlast2,nmin12,nminv2,nminv1,
     +               nlst21,nlastx
c
      common/symcon/norbbb(7),ndiff,iiii(18),ndimax
c
      integer nconst, ncond, maxpa, maxpat, maxsin, maxepa
      integer mxdvpa, mxddpa, nspips, nspipt, mxsspa, mxttpa
      integer mxstpa, nsitot, ntrtot
      integer mx12pa, mxvvpa, mxv2pa, maxpsq, max12
      common /auxh/ nconst(8),ncond(8),maxpa(4),maxpat,maxsin(3),
     +              maxepa(4),mxdvpa,mxddpa,nspips(8),nspipt(8),
     +              mxsspa,mxttpa,mxstpa,nsitot(8),ntrtot(8),
     +              mx12pa,mxvvpa,mxv2pa,maxpsq,max12
c
c
      integer nirr,mult,ifrep,nirr8,jfrep
      common /table/ nirr,mult(8,8),ifrep,nirr8,jfrep
c
c
      common/stoctl/khbill,khbm1e,khbm2e
c     common/spew/model,ic,ir,inti,intj,intk
      common/spew/rmodel(5)
      equivalence (rmodel(1),model),(rmodel(2),ic),(rmodel(3),ir)
c
      integer nxblk, index, numci, mblkci, ispbig, ispsma
      common/disktl/nxblk,index(199),numci,mblkci,ispbig,ispsma
      common /bufb/ ifi,ifo,iiiblk,jjjblk
c
      data m1/1/
c... spin coupling coefficient evaluation control routine
      call cpuwal(top,bot)
      write(6,444)top,bot
444   format(//1x,80('-')//
     *' ****hamiltonian builder vmp2 called at',f12.1,' wall',
     *f12.1,' secs')
      call brtsb(nxblk)
      index(1)=nxblk
c...
c... (ij/kl) and (ij)
c...
c.... vacuum-vacuum
      kcj = khbm2e+1
      nbuf=maxpat*maxpat
      kck=kcj+nbuf
      kbuf=kck+nbuf
      lbuf=kbuf+nbuf
      nuse=lbuf+nbuf-1
      call corfai('mphbill')
      model=2
      call mpijkl(cr(kcj),cr(kck),cr(kbuf),cr(lbuf),1,nnst,nref,icr,cr)
      call whersb(1)
c...   to get genz to clear doublet/n-2 vectors in out store mode
      model = 3
      ic = 0
      ir = 0
      call wrtsb(rmodel(1),3)
      call whersb(2)
      model = 4
      call wrtsb(rmodel(1),3)
      call whersb(3)
c...
c... (ij/ka) and (ia)
c...
      call basmak
      kcj=khbill+1
      kck=kcj+nbuf
      ifo=kck+nbuf-1
      ifi=ifo
      iiiblk=1
c... doublet - vacuum
      if (nmin1.eq.0) go to 7777
      call ijka(icr,cr(kcj),nstrt1,nlast1,1,nref)
c... singlet/triplet - doublet
      call whersb(5)
c...
c... (ij/ab) (ia/jb) and (ab)
c...
c... doublet - doublet
      call whersb(6)
c.. n-2  -  n-2
7777  if(nmin2.eq.0) go to 6666
      call whersb(7)
6666  call whersb(9)
       if (nmin2.eq.0) go to 2
c... singlet/triplet - vacuum
1     ifi=ifo
      iiiblk=1
      iflop = nnst
      nnst = nref
      call iajb(cr,icr,cr(kcj))
      nnst = iflop
      call ertsc(cr)
c...
2     model = 1
      call wrtsb(rmodel(1),1)
      call ertsb
c...
c... coupling coefficients for natural orbital generator
c...
c... model=1 end of coefficients
c...
c... spinfree natural orbitals
c... model=2 vacuum/vacuum  (ij)
c... model=3 vacuum/vacuum diagonal cases  (ii)
c... model=4 doublet/doublet (ij)
c... model=5 doublet/doublet diagonal cases (ii)
c... model=6 n-2/n-2 (ij)
c... model=7 n-2/n-2 diagonal cases (ii)
c... model=8 doublet/vacuum (ia)
c... model=9 n-2/doublet (ia)
c...
c... spin natural orbitals
c... model=10 vacuum/vacuum (ii) + (ij)
c... model=12 doublet/doublet (ii) + (ij)
c... model=14 n-2/n-2 (ii) + (ij)
c... model=15 doublet/vacuum (ia)
c... model=16 n-2/doublet (ia)
c... doublet/doublet and n-2/n-2 (ab) not labelled by model
      ndimax=2
      index(197)=iposun(numci)
      call brtsb(index(197))
c... vacuum/vacuum (ij)
      call natij(icr,cr(kcj),1,nnst,2)
c... doublet/doublet (ij)
      nn1=nnst
      do 44 isym=1,nirr
      ns=nn1+1
      nn1=nn1+ncond(isym)
      call natij(icr,cr(kcj),ns,nn1,4)
44    continue
c... n-2/n-2 (ij)
      do 55 isym=1,nirr
      ns=nn1+1
      nn1=nn1+nconst(isym)
      call natij(icr,cr(kcj),ns,nn1,6)
55    continue
c... doublet/vacuum (ia)
      call natia(cr,cr(kcj),nstrt1,nlast1,1,nnst,8)
c... n-2/doublet (ia)
      call natia(cr,cr(kcj),nstrt2,nlast2,nstrt1,nlast1,9)
      model = 1
      call wrtsb(rmodel(1),4)
c... doublet/doublet and n-2/n-2 (ab) spin natorb
      call natdub(icr,cr(kcj))
      call ertsb
      return
      end
      subroutine mpijkl(cj,ck,cx,cy,nstart,nlast,nref, conf,cr)
cmp
cmp...  only calculate matrix-elements involving reference confs.
cmp
c*ap*
      implicit real*8 (a-h,o-z)
      integer conf
      dimension cr(*),cj(*),ck(*),cx(*),cy(*), conf(5,*)
      dimension xnumb(2)
c     compute the h-matrix between vacuum states
c     **     this subroutine computes the internal parts for     **
c     ** the (n-1/h/n-1) and (n-2/h/n-2) matrix elements as well **
      common/erg/potnuc,totnuc,energy
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
      integer iro, iky, niky, m1e, m2e, m1ia, m2coul
      integer m2exch, m2fock, m2psup, m2qsup
      integer mbase, nbase, npairi, norbi, mbasi, mnbasi
      integer iap, mbfock, mbexch, mbcoul, mbpsup, mbqsup
      integer norbas, m2exab, m2em1e
      common /helpr/ iro(nd200),iky(nd200),niky(8),m1e,m2e,m1ia,m2coul,
     +   m2exch,m2fock,m2psup,m2qsup,
     +   mbase(288),nbase(8),npairi(36),norbi(8),mbasi(8),mnbasi(8),
     +   iap(62),mbfock(62),mbexch(36),mbcoul(36),mbpsup(36),mbqsup(36),
     +   norbas(nd200),m2exab,m2em1e
c
c     common/spew/model,ic,ir
      common/spew/rmodel(5)
      equivalence (rmodel(1),model),(rmodel(2),ic),(rmodel(3),ir)
      common/symcon/norbbb(7),nd,ii,jj,kk,ll,konte(14)
      common/three/n1l(1024),n2l(1024),minlhs(64),maxlhs(64),iocc(64),
     *n1r(1024),n2r(1024),minrhs(64),maxrhs(64),jocc(64),modelr
c
      integer nint, next, ninnex, nwm1, nw, nwp1, nintwo
      integer nint21, nintp1
      common /orb/ nint,next,ninnex,nwm1,nw,nwp1,nintwo,
     +             nint21,nintp1
c
c     equivalence (conf(1,1),cr(1))
      data xnumb/1.0d0,2.0d0/
      ihint(i,j) = iky(norbas(max(i,j))) + min(i,j)
     +             + nbase(iro(max(i,j)))
c
      if(nlast.lt.nstart)goto 99999
      do 100 ic=nstart,nlast
       call expan1(conf(1,ic))
       call upack2(conf(3,ic),ncolt,ncols)
cmp
         irend = min(ic,nref)
         do 100 ir=nstart,irend
cmp
      call expan2(conf(1,ir))
c.... establish number of differences
c.... if more then 4 => no future in it
c...    (nd comess from expan2 in /symcon/)
            if (nd.gt.4) go to 100
      call upack2(conf(3,ir),nrowt,nrows)
      ndims=nrows*ncols
      ndimt=nrowt*ncolt
      call vclr(cx,1,ndims)
      call vclr(cy,1,ndimt)
            if (nd-2) 20,60,90
c.... no differences => diagonal element
c....  ( h is the common  coulomb-exchange/1-el contribution )
20          h=potnuc
            do 50 i=1,nint
      iocci=iocc(i)
      if(iocci)51,50,51
51     iocci2=iocci-2
      pocci=xnumb(iocci)
               je = i - 1
      if(je)41,40,41
41            do 30 j=1,je
      ioccj=iocc(j)
      if(ioccj)31,30,31
31    coul=cr(irint(i,i,j,j))
      exch=cr(irint(i,j,i,j))
      ioccij=iocci2+ioccj
      if(ioccij)24,24,25
c.... coefficient  (both i and j are singly occupied)
24              call symb2(j,i)
      h=h+coul
                  call symbtr(cj,ck,ifl,0,0)
      if(ifl)7001,7000,7001
7001  call daxpy(ndims,exch,ck,1,cx,1)
7000  call symbtr(cj,ck,ifl,1,1)
      if(ifl)7006,30,7006
7006  call daxpy(ndimt,exch,ck,1,cy,1)
                  go to 30
c.... either one or both doubly occupied => (2j-k)*occ contribution
25    h=(coul+coul-exch)*xnumb(ioccij)+h
30             continue
c.... (i/i) and (ii/ii) contributions
40             h = h + cr(ihint(i,i))*pocci
      if(iocci2)50,7002,50
7002          h = h + cr(irint(i,i,i,i))
50          continue
c.... add coulomb/1-el contribution to diagonals
      if(ncols)7004,7003,7004
7004  call ihk(cx,h,ncols)
7003  if(ncolt)7005,400,7005
7005  call ihk(cy,h,ncolt)
            go to 400
c---- end of diagonal contributions --------------------------------
c.... single difference
60         h =cr(ihint(ii,jj))
c.... first do the orbitals with differences
      if(jocc(ii).eq.2)h=h+cr(irint(ii,ii,ii,jj))
      if(iocc(jj).eq.2)h=h+cr(irint(ii,jj,jj,jj))
c.... fock-type contributions
7016      do 80 i=1,nint
      iocci=iocc(i)
         if(iocci.eq.0.or.i.eq.ii.or.i.eq.jj) go to 80
      coul=cr(irint(ii,jj,i,i))
      exch=cr(irint(ii,i,jj,i))
      if(iocci-1)70,70,75
c.... i singly occupied => interesting coupling
70             call symb2(i,0)
               call symbtr(cj,ck,ifl,0,0)
       h=h+coul
      if(ifl)9001,9000,9001
9001  call daxpy(ndims,exch,ck,1,cx,1)
9000  call symbtr(cj,ck,ifl,1,1)
      if(ifl)9002,80,9002
9002  call daxpy(ndimt,exch,ck,1,cy,1)
      goto 80
c.... i doubly occupied => add contribution to h
75    h=h+coul+coul-exch
80          continue
c.... now add accumulated h + 2j - k contributions
      call symb1
            call sym1tr(ck,ifl,0,0)
      if(ifl)9010,9011,9010
9010  call daxpy(ndims,h,ck,1,cx,1)
9011  call sym1tr(ck,ifl,1,1)
      if(ifl)9012,400,9012
9012  call daxpy(ndimt,h,ck,1,cy,1)
      goto 400
c---- end of delta = 1 -------------------------------------------------
c.... 2 orthogonalities
90     call symb2(0,0)
            coul =cr(irint(ii,jj,kk,ll))
            exch =cr(irint(ii,kk,jj,ll))
      call symbtr(cj,ck,ifl,0,0)
      if(ifl)7041,7040,7041
7041  call ihj(cx,cj,ck,coul,exch,ndims)
7040  call symbtr(cj,ck,ifl,1,1)
      if(ifl)7042,400,7042
7042  call ihj(cy,cj,ck,coul,exch,ndimt)
400   call wrtsb(rmodel(1),3)
      call wrtsb(cx,ndims)
      call wrtsb(cy,ndimt)
c---- end of delta = 2 -------------------------------------------------
100   continue
99999 return
      end

      subroutine mpstar(conf,cr,icr)
c
      implicit real*8 (a-h,o-z)
c
      integer conf
      dimension conf(5,*),cr(*),icr(*)
c
c...  control routine for MR-MP2/MP3 symbolic and integral sorting
c
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
      common/mr_mp2l/ ivoc(62),nvoc,idoc(62),ndoc,kvcc(4),nvcc(4),
     *                ncmphn(10,18),ncmph1(10,18,8),ncmph2(10,18,8),
     *                ncmphd(10,18)
c
      real*8 critcep, crycep
      integer ifcepa, icepa , lcepaul, lsinsh, itcepa, lprcep, indcep
      integer nrfcsf, mcycep, ibrst,icorr_mr,mrdoc
      common/ccepa/critcep,crycep,ifcepa,icepa,lcepaul,lsinsh,itcepa,
     1                lprcep,indcep,nrfcsf,mcycep,ibrst,icorr_mr,mrdoc
c
c
      logical latorb, mspin1, ispind
      integer ispacg, isping, iprwh, iprcon
      common /natorb/ ispacg,isping,latorb,mspin1,ispind,
     +                iprwh,iprcon
c
      common/symcon/nornor(26),ndimax
c
c..   MR-MP model's used
c... model=1 end of coupling coefficients
c... model=2 ijkl vacuum/vacuum
c... model=3 ijkl doublet/doublet
c... model=4 ijkl n-2/n-2
c... model=35 fock(ia) doublet/vacuum
c... model=36 fock(ia) n-2/doublet
c... model=37 fock(ab) doublet/doublet  diagonal/off-diagonal (as ijab)
c... model=38 fock(ab) n-2/n-2 diagonal (i.e. as iaib)
c
c..  model 41 end of aa-couplings
c..  model 42 ijkl - aa
c..  model 43 fock(ia) aa/doublet
c..  model 44 fock(ab) aa-ab
c..  model 45 fock(ab) aa-aa
c
c     indmp : 1 : coupling coefficients / 2(ia) / 3 (ab)-triangle /
c             4 (ab)-square / 5 (bb)-diagonal / 6 aa coupling coeff.
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
      integer nirr,mult,ifrep,nirr8,jfrep
      common /table/ nirr,mult(8,8),ifrep,nirr8,jfrep
c
c
      integer nconst, ncond, maxpa, maxpat, maxsin, maxepa
      integer mxdvpa, mxddpa, nspips, nspipt, mxsspa, mxttpa
      integer mxstpa, nsitot, ntrtot
      integer mx12pa, mxvvpa, mxv2pa, maxpsq, max12
      common /auxh/ nconst(8),ncond(8),maxpa(4),maxpat,maxsin(3),
     +              maxepa(4),mxdvpa,mxddpa,nspips(8),nspipt(8),
     +              mxsspa,mxttpa,mxstpa,nsitot(8),ntrtot(8),
     +              mx12pa,mxvvpa,mxv2pa,maxpsq,max12
c
c
      integer nref, nnst, nmin1, nmin2, mspin, nelec
      integer nstrt1, nlast1, nstrt2, nlast2, nmin12, nminv2, nminv1
      integer nlst21, nlastx
      common /ccntl/ nref,nnst,nmin1,nmin2,mspin,nelec,
     +               nstrt1,nlast1,nstrt2,nlast2,nmin12,nminv2,nminv1,
     +               nlst21,nlastx
c
c
      integer nval, ndoub, nsingl, ntripl, ntotal, ibig, jbig, kbig
      integer nvd, nvds, nsitri, lbig, mbig, nbig, nref0
      common /cntrl/ nval,ndoub,nsingl,ntripl,ntotal,ibig,jbig,kbig,
     +               nvd,nvds,nsitri,lbig,mbig,nbig,nref0
c
c
      integer nint, next, ninnex, nwm1, nw, nwp1, nintwo
      integer nint21, nintp1
      common /orb/ nint,next,ninnex,nwm1,nw,nwp1,nintwo,
     +             nint21,nintp1
c
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
c
      integer iro, iky, niky, m1e, m2e, m1ia, m2coul
      integer m2exch, m2fock, m2psup, m2qsup
      integer mbase, nbase, npairi, norbi, mbasi, mnbasi
      integer iap, mbfock, mbexch, mbcoul, mbpsup, mbqsup
      integer norbas, m2exab, m2em1e
      common /helpr/ iro(nd200),iky(nd200),niky(8),m1e,m2e,m1ia,m2coul,
     +   m2exch,m2fock,m2psup,m2qsup,
     +   mbase(288),nbase(8),npairi(36),norbi(8),mbasi(8),mnbasi(8),
     +   iap(62),mbfock(62),mbexch(36),mbcoul(36),mbpsup(36),mbqsup(36),
     +   norbas(nd200),m2exab,m2em1e
c
c     common/spew/model,ic,ir
      common/spew/rmodel(5)
      equivalence (rmodel(1),model),(rmodel(2),ic),(rmodel(3),ir)

c
c
      real*8 tvec
      integer ivec, nvec, mxvec, iselvc, nselvc, istvc
      integer iprvc
      common /trial/ tvec(20),ivec(20),nvec,mxvec,
     +               iselvc,nselvc,istvc,iprvc
c
      common/stoctl/khbill,khbz,khbc,khbb
c
      integer nxblk, index, numci, mblkci, ispbig, ispsma
      common/disktl/nxblk,index(199),numci,mblkci,ispbig,ispsma
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer nuse, ngot
      common /corctl/ nuse, ngot
c
c
      dimension ncm(9)
c
      call cpuwal(top,bot)
      write(iwr,444)top,bot
444   format(//1x,80('-')//
     *' ****MR-MP prep. v1.2 called at',f12.1,' wall',f12.1,' secs')
c
      if (iprimp.ge.5) call corpr
c...  set iprcon just in case a config needs printing
      iprcon = 1
c
c...    make sure that nwm1 is ok for present purpose
c
      nnnnn = (nint+2-1)/32 + 1
      if (nnnnn.ne.nwm1) then
         write(iwr,1) nint,nnnnn,nwm1
1        format(/' *************************************************',/
     *          ,' **         little tested problem               **',/
     *          ,' **   nint ',i6,'  nw either ',i6,' or ',i6,'   **',/
     *          ,' *************************************************')
         nwm1 = nnnnn
      end if
c
c...   reserve space for fock-matrix and conf auxilaries
c...   put reference index array in front
c...   indexing for coefficients/conf (iadv/iadvt) before that
c...   **note** iadvt is only used for n-2
c
c...  Note on the core partitioning:
c...  Both array iadv and iadtt are (nlast2-1)/lenwrd()+1 integers
c...  long. However, iadtt contains information about triplets only.
c...  This means that the first nlast1/lenwrd() elements that are about
c...  the vacuum and the doublets are unused. To save memory iadv and
c...  iadtt are chosen to overlap such that iadtt(nstrt2) comes right
c...  after iadv(nlast2). 
c...  The array indr starts right after the last element of iadtt.
c
      kiadv  = khbill
      kiadvt = kiadv  + (nlast2-1)/lenwrd()+1
      kiadtt = kiadvt - nlast1/lenwrd()
      kindr  = kiadvt + (nlast2-nstrt2)/lenwrd()+2
c
      kfock = kindr + (nref+1)*3/lenwrd()+1
      kcon2 = kfock + ninnex*(ninnex+1)/2
      kcon3 = kcon2 + nnst*5
      kcr   = kcon3 + nnst*5
c
c...   set up various indexing / and set conf(4,..)
c
      call mpini(conf,cr(kcon2+1),cr(kcon3+1),
     *           cr(kiadv+1),cr(kiadtt+1),cr(kindr+1))
c
c...   get reference ci-vector  (setc1d)
c...   but give it all the core
c...   reference vector starts at kcr+1 (i.e. end of conf+ fock(ij)+..)
c...   use conf with adressing in place of occupation-number conf
c
      iselvc = 3
      nselvc = nrfcsf
      call dcopy(nnst*5,cr(kcon2+1),1,conf,1)
      call setc1d(cr(kcr+1),conf,e0ref,ngot-kcr,iwr)
c
c..   save energy and vector (krefmp is reserved in direct)
c
      call dcopy(nrfcsf,cr(kcr+1),1,cr(krefmp+1),1)
c
      write(iwr,445) e0ref,nrfcsf,(cr(krefmp+i),i=1,nrfcsf)
445   format(/' The reference Energy is ',e16.10,/,
     1        ' The reference Vector is : (',i6,' states)',/(t5,6f12.6))
c
c..  integral sort // no core problems expected
c
c..    read fock-matrix / and write ia/ab matrices
c
      kfia   = kcr
      kfabtr = kfia   + m1ia
      kfabsq = kfabtr + norbsi(1)
      nuse   = kfabsq + norbsq(1)
      call corfai('mp')
c
      call mpgij(cr(kfock+1),cr(kfia+1),cr(kfabtr+1),cr(kfabsq+1),
     *           conf,cr(krefmp+1),cr(nuse+1),ngot-nuse)
c
c...   the adressing (instead of configs) is kept because mpfock
c...   needs the addressing // now put conf back for mpsrcc and foll.
c
      call dcopy(nnst*5,cr(kcon3+1),1,conf,1)
c
c...    do'nt need extra conf arrays anymore // reset nuse
c...    only internal fock-matrix needed
c
      nuse = kfock + nint*(nint+1)/2
c
      call cpuwal(top,bot)
      write(iwr,447)top,bot
447   format(//1x,80('-')//
     *' ****MR-MP symbolic v1.3 called at',f12.1,' wall',f12.1,' secs')
c
c...  build contractors for various spaces
c...  depending on the model one or more may be actually absent
c...  in principle :  contractors onto singles,doubles for
c...  vacuum/n-1 and n-2
c...  amount of core taken by contractors in nuse  (start in kvcc())
c
      call mpprcc(conf,cr,cr(kindr+1),cr(kiadv+1),cr(kiadtt+1))
      if (iprimp.ge.5) call corpr
c
c...  --- symbolic for MR-MP ---
c
      kcx    = nuse
      kcy    = kcx  + maxpat**2
      kcxx   = kcy  + maxpat**2
      kcyy   = kcxx + maxpat**2
      kresul = kcyy + maxpat**2
c...   determine max # csfs / excited states and max. # conf's   per set
      maxcsf = 0
      maxexs = 0
      maxcon = 0
      do 10 i=1,2+nirr*2
         maxcsf = max(maxcsf,ncmphn(2,i))
         maxexs = max(maxexs,ncmphn(6,i),ncmphn(7,i),ncmphn(8,i))
         maxcon = max(maxcon,ncmphn(1,i))
         do 10 j=1,nirr
            maxcsf = max(maxcsf,ncmph1(2,i,j),ncmph2(2,i,j))
            maxexs = max(maxexs,ncmph1(6,i,j),ncmph1(7,i,j),
     *                    ncmph1(8,i,j),ncmph2(6,i,j),
     *                    ncmph2(7,i,j),ncmph2(8,i,j))
            maxcon = max(maxcon,ncmph1(1,i,j),ncmph2(1,i,j))
10    continue
c
      ktemp  = kresul + maxcsf**2
      kfinal = ktemp  + maxcsf*maxexs
      nuse   = kfinal + 2*maxexs**2
      knpp = nuse
      kipair = knpp + (nint+1)/lenwrd()
      nuse = kipair + (2*nint*maxcon+1)/lenwrd()
      call corfai('mp')
c
      nxblk = iposun(numci)
      ndimax = 2
c
c...  calculate FOCK - reference energy    <ref|F|ref>
c...  fill in dummy ncm
c
      ncm(1) = nref
      ncm(2) = nrfcsf
      ncm(3) = 1
      ncm(4) = nrfcsf
      ncm(5) = 0
      ncm(6) = 1
      ncm(7) = 0
      ncm(8) = 0
      ncm(9) = 1
c
      e0mp  = 0.0d0
      model = 2
      call brtsb(nxblk)
      call  mpfij(conf,1,nref,1,nref,1,1,
     *            cr(kiadv+1),cr(kiadtt+1),
     *            cr(krefmp+1),ncm,cr(krefmp+1),ncm,cr(kfock+1),
     *            cr(kresul+1),nrfcsf,nrfcsf,cr(kcx+1),cr(kcy+1),
     *            cr(ktemp+1),cr(kfinal+1))
      call ertsb
      call brtsb(nxblk)
      call getsb(rmodel(1),3)
      if (model.ne.2.or.ic.ne.1.or.ir.ne.1) call caserr('screw e0mp')
      call getsb(e0mp,1)
      write(iwr,446) e0mp
446   format(' the reference FOCK energy is ',e16.10)
      if (iprimp.ge.5) call corpr
c
      indmp(1) = nxblk
      call brtsb(indmp(1))
      icc = 1
c
      call cpuwal(top,bot)
      write(iwr,4471)top,bot
4471  format(' **** Fock IJ called at',f12.1,' wall',f12.1,' secs')
c
c...    vacuum   fock(ij)
c
      model = 2
      call  cmpfij(conf,1,icc,cr(kiadv+1),cr(kiadtt+1),
     *             cr(kvcc(1)),ncmphn,nirr,1,cr(kfock+1),cr(kresul+1),
     *             cr(kcx+1),cr(kcy+1),cr(ktemp+1),cr(kfinal+1))
c
      call whermp(iwnmp(1,1),iwnmp(2,1))
c
c...    doublet  fock(ij)
c
      model = 3
      icc1  = icc
      call cmpfij(conf,nstrt1,icc,cr(kiadv+1),cr(kiadtt+1),
     *            cr(kvcc(2)),ncmph1,nirr,nirr,cr(kfock+1),cr(kresul+1),
     *            cr(kcx+1),cr(kcy+1),cr(ktemp+1),cr(kfinal+1))
c
      call whermp(iwnmp(1,2),iwnmp(2,2))
c
c...    n-2     fock(ij)
c
      model = 4
      icc2  = icc
      call cmpfij(conf,nstrt2,icc,cr(kiadv+1),cr(kiadtt+1),
     *            cr(kvcc(3)),ncmph2,nirr,nirr,cr(kfock+1),cr(kresul+1),
     *            cr(kcx+1),cr(kcy+1),cr(ktemp+1),cr(kfinal+1))
c
      call whermp(iwnmp(1,3),iwnmp(2,3))
c
      iccd = icc
c
      call cpuwal(top,bot)
      write(iwr,4472)top,bot
4472  format(' **** Fock IA called at',f12.1,' wall',f12.1,' secs')
c
c...   doublet-vacuum fock(ia)
c
      model = 35
      call cmpfia(conf,nstrt1,1,icc1,1,cr(kiadv+1),cr(kiadtt+1),
     *            cr(kvcc(2)),ncmph1,cr(kvcc(1)),ncmphn,nirr,nirr,1,
     *            cr(kresul+1),cr(kcx+1),cr(ktemp+1),cr(kfinal+1),
     *            cr(knpp+1),cr(kipair+1),maxcon)
c
      call whermp(iwnmp(1,4),iwnmp(2,4))
c
c...   n-2 - doublet fock(ia)
c
      model = 36
      call cmpfia(conf,nstrt2,nstrt1,icc2,icc1,
     *            cr(kiadv+1),cr(kiadtt+1),
     *            cr(kvcc(3)),ncmph2,cr(kvcc(2)),ncmph1,nirr,nirr,nirr,
     *            cr(kresul+1),cr(kcx+1),cr(ktemp+1),cr(kfinal+1),
     *            cr(knpp+1),cr(kipair+1),maxcon)
c
      call whermp(iwnmp(1,5),iwnmp(2,5))
c
      call cpuwal(top,bot)
      write(iwr,4473)top,bot
4473  format(' **** Fock AB called at',f12.1,' wall',f12.1,' secs')
c
c...   doublet - doublet fock(ab)
c
      model = 37
      call cmpfab(icc1,cr(kvcc(2)),ncmph1,nirr,nirr,cr(kfinal+1))
c
      call whermp(iwnmp(1,6),iwnmp(2,6))
c
c...   n-2 - n-2  fock(ab)  (off diagonal)
c
      model = 38
      call cmpfab(icc2,cr(kvcc(3)),ncmph2,nirr,nirr,cr(kfinal+1))
c
      model = 1
      call wrtsb(rmodel(1),1)
      call ertsb
c
      call cpuwal(top,bot)
      write(iwr,4475)top,bot
4475  format(' **** Fock AA-DIAG called at',f12.1,' wall',f12.1,' secs')
c
c...   do aa diagonals //
c...   redo fij / ia and ab with aa diagonals involved
c...   could be made more efficient  // nlastx = last of symm1 n-2
c...   indmp(6) is switchable / could also depend on model
c...   note only n-2 symmetry 1 is involved
c
      if (nsingl.eq.0) then
         indmp(6) = 0
         go to 99999
      end if
c
      nxblk    = iposun(numci)
      indmp(6) = nxblk
      call brtsb(indmp(6))
c
c...    n-2     fock(ij)  (aa)
c
      model = 42
      icc   = iccd
      call cmpfij(conf,nstrt2,icc,cr(kiadv+1),cr(kiadtt+1),
     *            cr(kvcc(4)),ncmphd,nirr,1,cr(kfock+1),cr(kresul+1),
     *            cr(kcx+1),cr(kcy+1),cr(ktemp+1),cr(kfinal+1))
c
c     ncc = icc - 1
c
c...   n-2 - n-2  fock(ab)  aa/aa no code needed
c
      model = 45
      call wrtsb(rmodel(1),1)
c
c...   (aa/ab) mixed diagonal / off-diagonal
c
      if (ifdiaf.ne.1) then
        model = 44
        call cmpaab(iccd,cr(kvcc(4)),ncmphd,icc2,cr(kvcc(3)),ncmph2,
     *              nirr,cr(kfinal+1))
      end if
c
c...   n-2 - doublet fock(ia)
c
      model = 43
      call cmpfia(conf,nstrt2,nstrt1,iccd,icc1,
     *            cr(kiadv+1),cr(kiadtt+1),
     *            cr(kvcc(4)),ncmphd,cr(kvcc(2)),ncmph1,nirr,1,nirr,
     *            cr(kresul+1),cr(kcx+1),cr(ktemp+1),cr(kfinal+1),
     *            cr(knpp+1),cr(kipair+1),maxcon)
c
      model = 41
      call wrtsb(rmodel(1),1)
      call ertsb
c
99999 continue
c
c...  dump the contraction matrices to ed7
c
      do 99989 i=1,4
         nxblk = iposun(numci)
         call wrt3(cr(kvcc(i)),nvcc(i),nxblk,numci)
         kvcc(i) = nxblk
99989 continue
c
      call whtps
      nxblk = iposun(numci)
c
      if (iprimp.ge.5) call corpr
      call cpuwal(top,bot)
      write(iwr,448)top,bot
448   format(/,
     *' ****MR-MP symbolic v1.3 ended at',f12.1,' wall',f12.1,' secs')
      ndimax = 4
c
      return
      end
      subroutine mpini(conf,con2,con3,iadv,iadvt,indr)
c
      implicit real*8 (a-h,o-z)
c
      integer conf, con2, con3
c
      dimension conf(5,*)
      dimension con2(5,*),con3(5,*),iadv(*),iadvt(*),indr(3,*)
c
c...  initialise various arrays for mp2
c...  idoc : indices of orbitals doubly occupied in ref (ndoc)
c...  ivoc : indices of other internal orbitals  (nvoc)
c...  con2 : copy of vacuum conf with diabas adressing
c..        : for  vvijkl/setc1d
c...  con3 : save-keeping copy of vacuum conf
c...  indr : indexing of reference vector
c...  addressing vector in iadv :
c...  -- for n-2 singlet in iadv / triplet in iadvt
c
c...  also calculate few constants
c...    max. length of contractors (replaced by real length later)
c...    # reference/vacuum/n-1 and n-2 csf
c
c...   store fixed length data in /mp2l/
c
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
      common/mr_mp2l/ ivoc(62),nvoc,idoc(62),ndoc,kvcc(4),nvcc(4),
     *                ncmphn(10,18),ncmph1(10,18,8),ncmph2(10,18,8),
     *                ncmphd(10,18)
c
c
      integer nirr,mult,ifrep,nirr8,jfrep
      common /table/ nirr,mult(8,8),ifrep,nirr8,jfrep
c
c
      integer nref, nnst, nmin1, nmin2, mspin, nelec
      integer nstrt1, nlast1, nstrt2, nlast2, nmin12, nminv2, nminv1
      integer nlst21, nlastx
      common /ccntl/ nref,nnst,nmin1,nmin2,mspin,nelec,
     +               nstrt1,nlast1,nstrt2,nlast2,nmin12,nminv2,nminv1,
     +               nlst21,nlastx
c
c
      real*8 critcep, crycep
      integer ifcepa, icepa , lcepaul, lsinsh, itcepa, lprcep, indcep
      integer nrfcsf, mcycep, ibrst,icorr_mr,mrdoc
      common/ccepa/critcep,crycep,ifcepa,icepa,lcepaul,lsinsh,itcepa,
     1                lprcep,indcep,nrfcsf,mcycep,ibrst,icorr_mr,mrdoc
c
c
      integer nint, next, ninnex, nwm1, nw, nwp1, nintwo
      integer nint21, nintp1
      common /orb/ nint,next,ninnex,nwm1,nw,nwp1,nintwo,
     +             nint21,nintp1
c
c
c
c...  list of double and variable occupied orbitals made in mpsrcc
c
      if (iprimp.gt.10) then
        write(6,'(/,1x,/,2x,a)') '== orbital-classification =='
        if (ndoc.gt.0) write(6,601) 'doubly occupied',(idoc(i),i=1,ndoc)
        if (nvoc.gt.0) write(6,601) ' variably occ. ',(ivoc(i),i=1,nvoc)
601      format(2x,a15,(t20,10i5),1x)
      end if
c
c...  print configuration division
c
      if (iprimp.ge.30) then
         write(6,610)
610      format(/,5x,' **MP division of sorted configurations MP**',/)
         write(6,611) 'voc',('doc',i,i=1,nirr),'ddd',0,
     *               ('d-d',i,i=1,nirr)
         write(6,612) ' # conf     ',(ncmphn(1,i),i=1,2*nirr+2)
         write(6,613) ' # csf      ',(ncmphn(2,i),i=1,2*nirr+2)
         write(6,612) ' # sets     ',(ncmphn(3,i),i=1,2*nirr+2)
         write(6,612) ' # singlets ',(ncmphn(4,i),i=1,2*nirr+2)
         write(6,612) ' # triplets ',(ncmphn(5,i),i=1,2*nirr+2)
         write(6,612) ' # refs     ',(ncmphn(6,i),i=1,2*nirr+2)
         write(6,612) ' # singles  ',(ncmphn(7,i),i=1,2*nirr+2)
         write(6,612) ' # doubles  ',(ncmphn(8,i),i=1,2*nirr+2)
         write(6,612) ' # exc.state',(ncmphn(9,i),i=1,2*nirr+2)
         write(6,612) ' # exc.conf ',(ncmphn(10,i),i=1,2*nirr+2)
         write(6,611) '   '
         write(6,611) '   '
         do 80 ir=1,nirr
            write(6,612) ' # conf     ',(ncmph1(1,i,ir),i=1,2*nirr+2)
            write(6,614) ' n-1 ',ir,' # csf      ',
     *                   (ncmph1(2,i,ir),i=1,2*nirr+2)
            write(6,612) ' # sets     ',(ncmph1(3,i,ir),i=1,2*nirr+2)
            write(6,612) ' # singlets ',(ncmph1(4,i,ir),i=1,2*nirr+2)
            write(6,612) ' # triplets ',(ncmph1(5,i,ir),i=1,2*nirr+2)
            write(6,612) ' # refs     ',(ncmph1(6,i,ir),i=1,2*nirr+2)
            write(6,612) ' # singles  ',(ncmph1(7,i,ir),i=1,2*nirr+2)
            write(6,612) ' # doubles  ',(ncmph1(8,i,ir),i=1,2*nirr+2)
            write(6,612) ' # exc.state',(ncmph1(9,i,ir),i=1,2*nirr+2)
            write(6,612) ' # exc.conf ',(ncmph1(10,i,ir),i=1,2*nirr+2)
            write(6,611) '   '
80       continue
         write(6,611) '   '
         do 90 ir=1,nirr
            write(6,612) ' # conf     ',(ncmph2(1,i,ir),i=1,2*nirr+2)
            write(6,614) ' n-2 ',ir,' # csf      ',
     *                   (ncmph2(2,i,ir),i=1,2*nirr+2)
            write(6,612) ' # sets     ',(ncmph2(3,i,ir),i=1,2*nirr+2)
            write(6,612) ' # singlets ',(ncmph2(4,i,ir),i=1,2*nirr+2)
            write(6,612) ' # triplets ',(ncmph2(5,i,ir),i=1,2*nirr+2)
            write(6,612) ' # refs     ',(ncmph2(6,i,ir),i=1,2*nirr+2)
            write(6,612) ' # singles  ',(ncmph2(7,i,ir),i=1,2*nirr+2)
            write(6,612) ' # doubles  ',(ncmph2(8,i,ir),i=1,2*nirr+2)
            write(6,612) ' # exc.state',(ncmph2(9,i,ir),i=1,2*nirr+2)
            write(6,612) ' # exc.conf ',(ncmph2(10,i,ir),i=1,2*nirr+2)
            write(6,611) '   '
90       continue
         write(6,611) '   '
         write(6,612) ' # conf     ',(ncmphd(1,i),i=1,2*nirr+2)
         write(6,615) ' # csf      ',(ncmphd(2,i),i=1,2*nirr+2)
         write(6,612) ' # sets     ',(ncmphd(3,i),i=1,2*nirr+2)
         write(6,612) ' # singlets ',(ncmphd(4,i),i=1,2*nirr+2)
         write(6,612) ' # triplets ',(ncmphd(5,i),i=1,2*nirr+2)
         write(6,612) ' # refs     ',(ncmphd(6,i),i=1,2*nirr+2)
         write(6,612) ' # singles  ',(ncmphd(7,i),i=1,2*nirr+2)
         write(6,612) ' # doubles  ',(ncmphd(8,i),i=1,2*nirr+2)
         write(6,612) ' # exc.state',(ncmphd(9,i),i=1,2*nirr+2)
         write(6,612) ' # exc.conf ',(ncmphd(10,i),i=1,2*nirr+2)
         write(6,611) '   '
      end if
c
611   format(23x,a3,18(2x,a3,i1))
612   format(8x,a12,18i6)
613   format(' VACUUM ',a12,18i6)
614   format(1x,a5,i1,1x,a12,18i6)
615   format(' n-2-aa ',a12,18i6)
c
c...  run through conf to set up iadv adressing
c...  part of the work already done for sorting in mpsrcc
c...  so follow the indexing of the ncmph-arrays
c...  lproj : length of projection-matrix arrays
c...  also set up vacuum adressing and reference index-vector
c...
c...  since mpfia needs symmetry sorted out correctly check here
c...  one again (may be commented after a few month (iisym checks)
c
      nnrcsf = 0
      nncsf = 0
      nnst5 = nnst * 5
      call dcopy(nnst5,conf,1,con2,1)
      call dcopy(nnst5,conf,1,con3,1)
c...    ---vacuum---
      nvcc(1) = 0
      ic = 0
      do 10 in=1,2+nirr*2
         do 9 iset=1,ncmphn(3,in)
            nn = 0
            iisym = conf(4,ic+1)
            do 8 icon=1,ncmphn(1,in)
               ic = ic + 1
               if (iisym.ne.conf(4,ic)) call caserr('mpini-s1')
               call upack2(conf(3,ic),ncolt,ncols)
               if (ic.le.nref) then
                  indr(1,ic) = ic
                  indr(2,ic) = nnrcsf + 1
                  indr(3,ic) = ncols
                  nnrcsf = nnrcsf + ncols
               end if
               con2(1,ic) = nncsf
               iadv(ic) = nn
               nncsf = nncsf + ncols
               nn = nn + ncols
8           continue
c     debug check
         if (nn.ne.ncmphn(2,in)) call caserr('mpini1')
9        continue
         nvcc(1) = nvcc(1) + ncmphn(3,in)*(ncmphn(2,in)**2)
10    continue
c
      indr(1,nref+1) = 0
      indr(2,nref+1) = indr(2,nref) + indr(3,nref)
      indr(3,nref+1) = 0
c
c... doublets
c
      nvcc(2) = 0
      do 21 irconf=1,nirr
         do 20 in=1,2+nirr*2
            do 19 iset=1,ncmph1(3,in,irconf)
               nn = 0
               iisym = conf(4,ic+1)
               do 18 icon=1,ncmph1(1,in,irconf)
                  ic = ic + 1
                  if (iisym.ne.conf(4,ic)) call caserr('mpinis2')
                  call upack2(conf(3,ic),ncolt,ncols)
                  iadv(ic) = nn
                  nn = nn + ncols
18             continue
c     debug check
           if (nn.ne.ncmph1(2,in,irconf)) call caserr('mpini2')
19          continue
            nvcc(2) = nvcc(2) + ncmph1(3,in,irconf)
     *                        *(ncmph1(2,in,irconf)**2)
20       continue
21    continue
c
c...  n-2
c
      nvcc(3) = 0
      nvcc(4) = 0
      do 31 irconf=1,nirr
         do 30 in=1,2+nirr*2
            do 29 iset=1,ncmph2(3,in,irconf)
               nn = 0
               mm = ncmph2(4,in,irconf)
c...  singlet/triplet
               iisym = conf(4,ic+1)
               do 27 icon=1,ncmph2(1,in,irconf)
                  ic = ic + 1
                  if (iisym.ne.conf(4,ic)) call caserr('mpinis3')
                  call upack2(conf(3,ic),ncolt,ncols)
                  iadv(ic) = nn
                  nn = nn + ncols
                  iadvt(ic) = mm
                  mm = mm + ncolt
27             continue
c     debug check
            if (nn.ne.ncmph2(4,in,irconf)) call caserr('mpini3')
            if (mm.ne.ncmph2(2,in,irconf)) call caserr('mpini4')
29          continue
            nvcc(3) = nvcc(3) + ncmph2(3,in,irconf)
     *                        *(ncmph2(2,in,irconf)**2)
            nvcc(4) = nvcc(4) + ncmph2(3,in,irconf)
     *                        *(ncmph2(4,in,irconf)**2)
30       continue
31    continue
c
c...   debug
c
      if (ic.ne.nlast2) call caserr('mpini5')
c
      if (nnrcsf.ne.nrfcsf) call caserr('nrcsf ne nrfcsf')
c
      return
      end
      subroutine prmpcon
c
      implicit real*8 (a-h,o-z)
c
c...  print various arrays for mp2
c...  idoc : indices of orbitals doubly occupied in ref (ndoc)
c...  ivoc : indices of other internal orbitals  (nvoc)
c...  con2 : copy of vacuum conf with diabas adressing
c..        : for  vvijkl/setc1d
c...  con3 : save-keeping copy of vacuum conf
c...  indr : indexing of reference vector
c...  addressing vector in iadv :
c...  -- for n-2 singlet in iadv / triplet in iadvt
c
c...  the  ncmph  array's contain info for the sorted subdivisions :
c...  (also results from the excitation generator)
c...   ncmphn : vacuum ; ncmph1 : n-1 ; ncmph2 : n-2 ; ncmphd : n-2-aa
c...   1st index : # conf's (1,...) / #  csf's (2,...) / # sets (3....)
c...               # singlets (4,..) / # triplets (5,..)
c...               # refs (6,..) / # singles (7,..) / # doubles (8,..)
c...               # excited states (9,..) /total # exc. confs (10,..)
c...   2nd index : # holes in doubly occ. :
c...                0 holes (..,1,..)
c...                1 hole (..,2 - nirr+1,..)
c...                2 holes in 1 orbital (..,2+nirr,..)
c...                2 holes in 2 orbitals (..,2+nirr+1 - 2+nirr+nirr,..)
c...   3d  index : symmetry of conf (only for n-1 and n-2)
c...   clearly not all info is relevant to each type
c...   ncmphd is almost identical to ncmph2 (only symmetry 1 though)
c
c...   store fixed length data in /mp2l/
c
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
      common/mr_mp2l/ ivoc(62),nvoc,idoc(62),ndoc,kvcc(4),nvcc(4),
     *                ncmphn(10,18),ncmph1(10,18,8),ncmph2(10,18,8),
     *                ncmphd(10,18)
c
      integer nirr,mult,ifrep,nirr8,jfrep
      common /table/ nirr,mult(8,8),ifrep,nirr8,jfrep
c
c
      integer nref, nnst, nmin1, nmin2, mspin, nelec
      integer nstrt1, nlast1, nstrt2, nlast2, nmin12, nminv2, nminv1
      integer nlst21, nlastx
      common /ccntl/ nref,nnst,nmin1,nmin2,mspin,nelec,
     +               nstrt1,nlast1,nstrt2,nlast2,nmin12,nminv2,nminv1,
     +               nlst21,nlastx
c
c
      integer nint, next, ninnex, nwm1, nw, nwp1, nintwo
      integer nint21, nintp1
      common /orb/ nint,next,ninnex,nwm1,nw,nwp1,nintwo,
     +             nint21,nintp1
c
c
c...  list of double and variable occupied orbitals made in mpsrcc
c
c     if (iprimp.gt.10) then
        write(6,'(/,1x,/,2x,a)') '== orbital-classification =='
        if (ndoc.gt.0) write(6,601) 'doubly occupied',(idoc(i),i=1,ndoc)
        if (nvoc.gt.0) write(6,601) ' variably occ. ',(ivoc(i),i=1,nvoc)
601      format(2x,a15,(t20,10i5),1x)
c     end if
c
c...  print configuration division
c
c     if (iprimp.ge.30) then
         write(6,610)
610      format(/,5x,' **MP division of sorted configurations MP**',/)
         write(6,611) 'voc',('doc',i,i=1,nirr),'ddd',0,
     *               ('d-d',i,i=1,nirr)
         write(6,612) ' # conf     ',(ncmphn(1,i),i=1,2*nirr+2)
         write(6,613) ' # csf      ',(ncmphn(2,i),i=1,2*nirr+2)
         write(6,612) ' # sets     ',(ncmphn(3,i),i=1,2*nirr+2)
         write(6,612) ' # singlets ',(ncmphn(4,i),i=1,2*nirr+2)
         write(6,612) ' # triplets ',(ncmphn(5,i),i=1,2*nirr+2)
         write(6,612) ' # refs     ',(ncmphn(6,i),i=1,2*nirr+2)
         write(6,612) ' # singles  ',(ncmphn(7,i),i=1,2*nirr+2)
         write(6,612) ' # doubles  ',(ncmphn(8,i),i=1,2*nirr+2)
         write(6,612) ' # exc.state',(ncmphn(9,i),i=1,2*nirr+2)
         write(6,612) ' # exc.conf ',(ncmphn(10,i),i=1,2*nirr+2)
         write(6,611) '   '
         write(6,611) '   '
         do 80 ir=1,nirr
            write(6,612) ' # conf     ',(ncmph1(1,i,ir),i=1,2*nirr+2)
            write(6,614) ' n-1 ',ir,' # csf      ',
     *                   (ncmph1(2,i,ir),i=1,2*nirr+2)
            write(6,612) ' # sets     ',(ncmph1(3,i,ir),i=1,2*nirr+2)
            write(6,612) ' # singlets ',(ncmph1(4,i,ir),i=1,2*nirr+2)
            write(6,612) ' # triplets ',(ncmph1(5,i,ir),i=1,2*nirr+2)
            write(6,612) ' # refs     ',(ncmph1(6,i,ir),i=1,2*nirr+2)
            write(6,612) ' # singles  ',(ncmph1(7,i,ir),i=1,2*nirr+2)
            write(6,612) ' # doubles  ',(ncmph1(8,i,ir),i=1,2*nirr+2)
            write(6,612) ' # exc.state',(ncmph1(9,i,ir),i=1,2*nirr+2)
            write(6,612) ' # exc.conf ',(ncmph1(10,i,ir),i=1,2*nirr+2)
            write(6,611) '   '
80       continue
         write(6,611) '   '
         do 90 ir=1,nirr
            write(6,612) ' # conf     ',(ncmph2(1,i,ir),i=1,2*nirr+2)
            write(6,614) ' n-2 ',ir,' # csf      ',
     *                   (ncmph2(2,i,ir),i=1,2*nirr+2)
            write(6,612) ' # sets     ',(ncmph2(3,i,ir),i=1,2*nirr+2)
            write(6,612) ' # singlets ',(ncmph2(4,i,ir),i=1,2*nirr+2)
            write(6,612) ' # triplets ',(ncmph2(5,i,ir),i=1,2*nirr+2)
            write(6,612) ' # refs     ',(ncmph2(6,i,ir),i=1,2*nirr+2)
            write(6,612) ' # singles  ',(ncmph2(7,i,ir),i=1,2*nirr+2)
            write(6,612) ' # doubles  ',(ncmph2(8,i,ir),i=1,2*nirr+2)
            write(6,612) ' # exc.state',(ncmph2(9,i,ir),i=1,2*nirr+2)
            write(6,612) ' # exc.conf ',(ncmph2(10,i,ir),i=1,2*nirr+2)
            write(6,611) '   '
90       continue
         write(6,611) '   '
         write(6,612) ' # conf     ',(ncmphd(1,i),i=1,2*nirr+2)
         write(6,615) ' # csf      ',(ncmphd(2,i),i=1,2*nirr+2)
         write(6,612) ' # sets     ',(ncmphd(3,i),i=1,2*nirr+2)
         write(6,612) ' # singlets ',(ncmphd(4,i),i=1,2*nirr+2)
         write(6,612) ' # triplets ',(ncmphd(5,i),i=1,2*nirr+2)
         write(6,612) ' # refs     ',(ncmphd(6,i),i=1,2*nirr+2)
         write(6,612) ' # singles  ',(ncmphd(7,i),i=1,2*nirr+2)
         write(6,612) ' # doubles  ',(ncmphd(8,i),i=1,2*nirr+2)
         write(6,612) ' # exc.state',(ncmphd(9,i),i=1,2*nirr+2)
         write(6,612) ' # exc.conf ',(ncmphd(10,i),i=1,2*nirr+2)
         write(6,611) '   '
c     end if
c
611   format(23x,a3,18(2x,a3,i1))
612   format(8x,a12,18i6)
613   format(' VACUUM ',a12,18i6)
614   format(1x,a5,i1,1x,a12,18i6)
615   format(' n-2-aa ',a12,18i6)
      return
      end
c*******************MR-MP symbolic *************************************
      subroutine mpprcc(conf,cr,indr,iadv,iadvt)
c
      implicit real*8 (a-h,o-z)
      integer conf
c
      dimension conf(5,*),cr(*),indr(*),iadv(*),iadvt(*)
c
c...   control routine for building contracters onto
c...   reference,singles,and doubles spaces
c...   kcr : start (-1)  of free area in cr
c...   kvv/nvv return  adresses/length of contraction arrays
c...   nuse passes end of space for contractors
c
c       projection models   (cf. mpini)  for the FOCK operator
c       0   ** no projection applied **
c       1   |ref>F<ref|+|sin>F<sin|+|doub>F<doub|
c       2   |ref>F<ref| + |doub>F<doub|
c       3   |ref>F<ref| + |sin+doub>F<sin+doub|
c       4   |ref>F<ref|+|sin>F<sin|+|rest>F<rest|
c       5   |ref>F<ref| + |rest>F<rest|
c       6   |ref>F<ref| + |doub+more>F<doub+more|
c       11  |ref>F<ref|+|sin'>F<sin'|+|doub'>F<doub'|
c       12  |ref>F<ref| + |doub'>F<doub'|
c       13  |ref>F<ref| + |sin'+doub'>F<sin'+doub'|
c
c       only models 1-3 and 11-13 implemented
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
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
      common/mr_mp2l/ ivoc(62),nvoc,idoc(62),ndoc,kvcc(4),nvcc(4),
     *                ncmphn(10,18),ncmph1(10,18,8),ncmph2(10,18,8),
     *                ncmphd(10,18)
c
      real*8 critcep, crycep
      integer ifcepa, icepa , lcepaul, lsinsh, itcepa, lprcep, indcep
      integer nrfcsf, mcycep, ibrst,icorr_mr,mrdoc
      common/ccepa/critcep,crycep,ifcepa,icepa,lcepaul,lsinsh,itcepa,
     1                lprcep,indcep,nrfcsf,mcycep,ibrst,icorr_mr,mrdoc
c
c
      integer nirr,mult,ifrep,nirr8,jfrep
      common /table/ nirr,mult(8,8),ifrep,nirr8,jfrep
c
c
      integer iro, iky, niky, m1e, m2e, m1ia, m2coul
      integer m2exch, m2fock, m2psup, m2qsup
      integer mbase, nbase, npairi, norbi, mbasi, mnbasi
      integer iap, mbfock, mbexch, mbcoul, mbpsup, mbqsup
      integer norbas, m2exab, m2em1e
      common /helpr/ iro(nd200),iky(nd200),niky(8),m1e,m2e,m1ia,m2coul,
     +   m2exch,m2fock,m2psup,m2qsup,
     +   mbase(288),nbase(8),npairi(36),norbi(8),mbasi(8),mnbasi(8),
     +   iap(62),mbfock(62),mbexch(36),mbcoul(36),mbpsup(36),mbqsup(36),
     +   norbas(nd200),m2exab,m2em1e
c
c
      integer nconst, ncond, maxpa, maxpat, maxsin, maxepa
      integer mxdvpa, mxddpa, nspips, nspipt, mxsspa, mxttpa
      integer mxstpa, nsitot, ntrtot
      integer mx12pa, mxvvpa, mxv2pa, maxpsq, max12
      common /auxh/ nconst(8),ncond(8),maxpa(4),maxpat,maxsin(3),
     +              maxepa(4),mxdvpa,mxddpa,nspips(8),nspipt(8),
     +              mxsspa,mxttpa,mxstpa,nsitot(8),ntrtot(8),
     +              mx12pa,mxvvpa,mxv2pa,maxpsq,max12
c
c
      integer nref, nnst, nmin1, nmin2, mspin, nelec
      integer nstrt1, nlast1, nstrt2, nlast2, nmin12, nminv2, nminv1
      integer nlst21, nlastx
      common /ccntl/ nref,nnst,nmin1,nmin2,mspin,nelec,
     +               nstrt1,nlast1,nstrt2,nlast2,nmin12,nminv2,nminv1,
     +               nlst21,nlastx
c
c
      integer nval, ndoub, nsingl, ntripl, ntotal, ibig, jbig, kbig
      integer nvd, nvds, nsitri, lbig, mbig, nbig, nref0
      common /cntrl/ nval,ndoub,nsingl,ntripl,ntotal,ibig,jbig,kbig,
     +               nvd,nvds,nsitri,lbig,mbig,nbig,nref0
c
c
      integer nint, next, ninnex, nwm1, nw, nwp1, nintwo
      integer nint21, nintp1
      common /orb/ nint,next,ninnex,nwm1,nw,nwp1,nintwo,
     +             nint21,nintp1
c
c
c
      integer nuse, ngot
      common /corctl/ nuse, ngot
c
c
c
c...  get an acceptable minimum size for maxpat (gensin will warn)
c
      maxpat = max(maxpat,mspin+1)
c
c..   treat all models
c
      if (abs(mpmod).lt.1.or.(abs(mpmod).gt.3.and.abs(mpmod).lt.11).or.
     *    abs(mpmod).gt.13) 
     *    call caserr('wrong model in mpprcc')
c
      if (iprimp.ge.500) call mpexpr(1,cr(krefmp+1),nrfcsf,indr,nref+1,
     *                               1,1,1,1,1,1,'r')
c
c...    1 |ref>F<ref|+|sin>F<sin|+|doub>F<doub|
c...    2 |ref>F<ref| + |doub>F<doub|
c...    3 |ref>F<ref| + |sin+doub>F<sin+doub|
c
c...    space for ioksi (max # singles // nvoc**2 or 1 assigned /type
c
c...   ------------------- vacuuum----------------------------
c
c...    singles/doubles in nvoc  (max # singles = nvoc**2)
c...    leave space for reference at beginning of ind/v
c...    excited states array starts at kvcn
c
      kvcn = nuse+1
      nuse = kvcn + nvcc(1)
      kioksi = nuse
      nuse = kioksi + max(nvoc*nvoc,1)/lenwrd() + 1
c
      if (iprimp.ge.500) write(6,601) 'vacuum','voc-voc',1
601   format(/' ++++++++ excite ',a7,3x,a7,' symmetry ',i1,' ++++++++')
c
      icst = 1
      kcp = 0
      call mppro(conf,icst,iadv,iadvt,ncmphn(1,1),dum,
     *           cr,cr,cr(krefmp+1),indr,cr(kioksi+1),
     *           ivoc,nvoc,ivoc,nvoc,ivoc,nvoc,ivoc,nvoc,
     *           cr(kvcn+kcp),dum,1,iabs(mpmod),.true.)
      icst = icst + ncmphn(1,1)
      kcp = kcp + ncmphn(2,1)*ncmphn(9,1)
c
c...    vacuum single holes in doc  (max # doubles nvoc**3*nref)
c
      if (iprimp.ge.500) write(6,601) 'vacuum','doc/voc',1
      do 10 i=1,ndoc
         irr = iro(idoc(i))+1
         if (ncmphn(1,irr).gt.0)
     *   call mppro(conf,icst,iadv,iadvt,ncmphn(1,irr),dum,
     *              cr,cr,cr(krefmp+1),indr,cr(kioksi+1),
     *              idoc(i),1,ivoc,nvoc,ivoc,nvoc,ivoc,nvoc,
     *              cr(kvcn+kcp),dum,0,iabs(mpmod),.true.)
         icst = icst + ncmphn(1,irr)
         kcp = kcp + ncmphn(2,irr)*ncmphn(9,irr)
10    continue
c
c...    vacuum double holes in  1 doc  // no singles needed
c
      if (iprimp.ge.500) write(6,601) 'vacuum','ddddoc ',1
      do 20 i=1,ndoc
         irr = nirr+2
         if (ncmphn(1,irr).gt.0)
     *   call mppro(conf,icst,iadv,iadvt,ncmphn(1,irr),dum,
     *              cr,cr,cr(krefmp+1),indr,cr(kioksi+1),
     *              idoc(i),1,ivoc,nvoc,idoc(i),1,ivoc,nvoc,
     *              cr(kvcn+kcp),dum,0,iabs(mpmod),.false.)
         icst = icst + ncmphn(1,irr)
         kcp = kcp + ncmphn(2,irr)*ncmphn(9,irr)
20    continue
c
c...    vacuum double holes in different  doc  // no singles needed
c...    because of sorting // handle per symmetry of doc-doc
c
      do 30 ir=1,nirr
        if (iprimp.ge.500) write(6,601) 'vacuum','doc-doc',ir
         irr = ir+nirr+2
         if (ncmphn(1,irr).eq.0)  go to 30
          do 25 i=2,ndoc
           do 25 j=1,i-1
             if (mult(iro(idoc(i)),iro(idoc(j))).ne.ir) go to 25
             call mppro(conf,icst,iadv,iadvt,ncmphn(1,irr),dum,
     *                  cr,cr,cr(krefmp+1),indr,cr(kioksi+1),
     *                  idoc(i),1,ivoc,nvoc,idoc(j),1,ivoc,nvoc,
     *                  cr(kvcn+kcp),dum,0,iabs(mpmod),.false.)
             icst = icst + ncmphn(1,irr)
             kcp = kcp + ncmphn(2,irr)*ncmphn(9,irr)
25        continue
30    continue
c
c...    debug check
c
      if (icst.ne.nstrt1) call caserr('nstrt1 problem in mpprcc')
c
c...    -------------------------- (n-1) ---------------------------
c
c...     as vacuum except that we do it per symmetry of the config
c...     start of excited vectors : kvc1
c
      kvc1 = kvcn + kcp
      nuse = kvc1 + nvcc(2)
      kioksi = nuse
      nuse = kioksi + max(nvoc*nvoc,1)/lenwrd() + 1
      kcp = 0
c
      do 100 irconf=1,nirr
c
         if (iprimp.ge.500) write(6,602) irconf
602      format(/' ############ n-1 symmetry ',i1,' ##############')
c
c...     voc => a // voc => voc
c
         if (iprimp.ge.500) write(6,601) '  n-1 ','voc-voc',1
         if (ncmph1(1,1,irconf).gt.0)
     *   call mppro(conf,icst,iadv,iadvt,ncmph1(1,1,irconf),dum,
     *              cr,cr,cr(krefmp+1),indr,cr(kioksi+1),
     *              ivoc,nvoc,nint+2,1,ivoc,nvoc,ivoc,nvoc,
     *              cr(kvc1+kcp),dum,0,iabs(mpmod),.true.)
         icst = icst + ncmph1(1,1,irconf)
         kcp = kcp + ncmph1(2,1,irconf)*ncmph1(9,1,irconf)
c
c...     doc => a // voc  => voc
c
         if (iprimp.ge.500) write(6,601) '  n-1 ','doc/voc',1
         do 40 i=1,ndoc
          irr = iro(idoc(i))+1
          if (ncmph1(1,irr,irconf).gt.0)
     *    call mppro(conf,icst,iadv,iadvt,ncmph1(1,irr,irconf),dum,
     *               cr,cr,cr(krefmp+1),indr,cr(kioksi+1),
     *               idoc(i),1,nint+2,1,ivoc,nvoc,ivoc,nvoc,
     *               cr(kvc1+kcp),dum,0,iabs(mpmod),.true.)
          icst = icst + ncmph1(1,irr,irconf)
          kcp = kcp + ncmph1(2,irr,irconf)*ncmph1(9,irr,irconf)
40       continue
c
c...     doc => a // same doc => voc
c
         if (iprimp.ge.500) write(6,601) '  n-1 ','dddddoc',1
         do 50 i=1,ndoc
          irr = nirr+2
          if (ncmph1(1,irr,irconf).gt.0)
     *    call mppro(conf,icst,iadv,iadvt,ncmph1(1,irr,irconf),dum,
     *               cr,cr,cr(krefmp+1),indr,cr(kioksi+1),
     *               idoc(i),1,nint+2,1,idoc(i),1,ivoc,nvoc,
     *               cr(kvc1+kcp),dum,0,iabs(mpmod),.false.)
          icst = icst + ncmph1(1,irr,irconf)
          kcp = kcp + ncmph1(2,irr,irconf)*ncmph1(9,irr,irconf)
50       continue
c
c...     doc => a // other doc => voc  (per doc-symmetry)
c
         do 60 ir=1,nirr
          irr = ir+nirr+2
          if (ncmph1(1,irr,irconf).eq.0)  go to 60
          if (iprimp.ge.500) write(6,601) '  n-1 ','doc/doc',ir
            do 55 i=2,ndoc
             do 55 j=1,i-1
                if (mult(iro(idoc(i)),iro(idoc(j))).ne.ir) go to 55
                call mppro(conf,icst,iadv,iadvt,
     *                     ncmph1(1,irr,irconf),dum,
     *                     cr,cr,cr(krefmp+1),indr,cr(kioksi+1),
     *                     idoc(i),1,nint+2,1,idoc(j),1,ivoc,nvoc,
     *                     cr(kvc1+kcp),dum,0,iabs(mpmod),.false.)
                icst = icst + ncmph1(1,irr,irconf)
                kcp = kcp + ncmph1(2,irr,irconf)*ncmph1(9,irr,irconf)
55          continue
60       continue
100   continue
c
c...    debug check
c
      if (icst.ne.nstrt2) call caserr('nstrt2 problem in mpprcc')
c
c...    -------------------------- (n-2) ---------------------------
c...      start of n-2 : kvc2 / of n-2-aa kvc2d
c
      kvc2 = kvc1 + kcp
      kvc2d = kvc2 + nvcc(3)
      nuse = kvc2d + nvcc(4)
      kioksi = nuse
      nuse = kioksi + max(nvoc*nvoc,1)/lenwrd() + 1
c
c...      start of n-2 : kvc2 / of n-2-aa kvc2d
c...      for easier debug modified to nint+1/nint+2
c
      kcp = 0
      kcpd = 0
c
      do 200 irconf=1,nirr
c
         if (iprimp.ge.500) write(6,603) irconf
603      format(/' ############ n-2 symmetry ',i1,' ##############')
c
c...     voc => a // voc => b
c
         if (iprimp.ge.500) write(6,601) '  n-2 ','voc-voc',1
         if (ncmph2(1,1,irconf).gt.0)
     *   call mppro(conf,icst,iadv,iadvt,
     *              ncmph2(1,1,irconf),ncmphd(1,1),
     *              cr,cr,cr(krefmp+1),indr,cr(kioksi+1),
     *              ivoc,nvoc,nint+1,1,ivoc,nvoc,nint+2,1,
     *              cr(kvc2+kcp),cr(kvc2d+kcpd),0,2,.false.)
         icst = icst + ncmph2(1,1,irconf)
         kcp = kcp + ncmph2(2,1,irconf)*ncmph2(9,1,irconf)
         if (irconf.eq.1) kcpd = kcpd + ncmphd(2,1)*ncmphd(9,1)
c
c...     doc => a // voc => b
c
         if (iprimp.ge.500) write(6,601) '  n-2 ','voc/doc',1
         do 110 i=1,ndoc
            irr = iro(idoc(i))+1
            if (ncmph2(1,irr,irconf).gt.0)
     *      call mppro(conf,icst,iadv,iadvt,
     *                 ncmph2(1,irr,irconf),ncmphd(1,irr),
     *                 cr,cr,cr(krefmp+1),indr,cr(kioksi+1),
     *                 idoc(i),1,nint+1,1,ivoc,nvoc,nint+2,1,
     *                 cr(kvc2+kcp),cr(kvc2d+kcpd),0,2,.false.)
            icst = icst + ncmph2(1,irr,irconf)
            kcp = kcp + ncmph2(2,irr,irconf)*ncmph2(9,irr,irconf)
            if (irconf.eq.1) kcpd = kcpd + ncmphd(2,irr)*ncmphd(9,irr)
110      continue
c
c...     doc => a // same doc => b
c
         if (iprimp.ge.500) write(6,601) '  n-2 ','dddddoc',1
         do 120 i=1,ndoc
            irr = nirr+2
            if (ncmph2(1,irr,irconf).gt.0)
     *      call mppro(conf,icst,iadv,iadvt,
     *                 ncmph2(1,irr,irconf),ncmphd(1,irr),
     *                 cr,cr,cr(krefmp+1),indr,cr(kioksi+1),
     *                 idoc(i),1,nint+1,1,idoc(i),1,nint+2,1,
     *                 cr(kvc2+kcp),cr(kvc2d+kcpd),0,2,.false.)
            icst = icst + ncmph2(1,irr,irconf)
            kcp = kcp + ncmph2(2,irr,irconf)*ncmph2(9,irr,irconf)
            if (irconf.eq.1) kcpd = kcpd + ncmphd(2,irr)*ncmphd(9,irr)
120      continue
c
c...     doc => a // other doc => b   (per symmetry of the docs)
c
         do 130 ir=1,nirr
          irr = ir+nirr+2
          if (ncmph2(1,irr,irconf).eq.0)  go to 130
          if (iprimp.ge.500) write(6,601) '  n-2 ','doc-doc',ir
          do 125 i=2,ndoc
             do 125 j=1,i-1
                if (mult(iro(idoc(i)),iro(idoc(j))).ne.ir) go to 125
                call mppro(conf,icst,iadv,iadvt,
     *                     ncmph2(1,irr,irconf),ncmphd(1,irr),
     *                     cr,cr,cr(krefmp+1),indr,cr(kioksi+1),
     *                     idoc(i),1,nint+1,1,idoc(j),1,nint+2,1,
     *                     cr(kvc2+kcp),cr(kvc2d+kcpd),0,2,.false.)
                icst = icst + ncmph2(1,irr,irconf)
                kcp = kcp + ncmph2(2,irr,irconf)*ncmph2(9,irr,irconf)
                if(irconf.eq.1)kcpd = kcpd+ncmphd(2,irr)*ncmphd(9,irr)
125        continue
130      continue
c
200   continue
c
c...    debug check
c
         if (icst.ne.nlast2+1) call caserr('nslast2 problem in mpprcc')
c
c...    --------------------------  end  ---------------------------
c
c
c...  compact the vectors
c
      kk = kvc2 + kcp
      do 310 i=1,kcpd
310   cr(kk+i-1) = cr(kvc2d+i-1)
      kvc2d = kk
c
c..   store present addresses and lengths of the contraction matrices
c
      kvcc(1) = kvcn
      kvcc(2) = kvc1
      kvcc(3) = kvc2
      kvcc(4) = kvc2d
      nvcc(1) = kvc1 - kvcn
      nvcc(2) = kvc2 - kvc1
      nvcc(3) = kcp
      nvcc(4) = kcpd
c
      nuse =  kvcc(4) + nvcc(4) - 1
c
c...  all the contraction matrices are produced
c
      return
      end
      subroutine mppro(conf,icst,iadv,iadvt,ncmph,ncmphd,
     *                 cr,icr,ref,indr,ioksi,
     *                 ioc1,noc1,ioc2,noc2,ioc3,noc3,ioc4,noc4,
     *                 vcont,vcontd,iref,mode,osin)
c
      implicit real*8 (a-h,o-z)
c
      integer  conf
c
      dimension conf(5,*),cr(*),icr(*),ref(*),indr(*),ioksi(noc1*noc2)
      dimension ioc1(noc1),ioc2(noc2),ioc3(noc3),ioc4(noc4)
      dimension vcont(*),vcontd(*)
      dimension iadv(*),iadvt(*),ncmph(10),ncmphd(10)
      logical osin
c
c...  generate excited states corresponding to excitations
c...  from the ioc1 set to the ioc2 set and then ioc3 to ioc4
c...  icst : 1st conf of the ncst conf the excited state may contain
c...  iref (1,0) indicates if reference-vector should be in
c...             set of vectors for schmidt
c...  osin .true. generate singles / .false. all singles are forbidden
c
c...  results are return as orthogonalised excited states in vcont
c...  and the number of refererence,singles and doubles in ncmph
c...  vcontd,ncmphd is for the n-2-diagonals
c
c...   cr and icr are the same array (cr(1)/icr(1),icr(2) etc.
c...   depending on lenwrd() )
c
c
      integer nref, nnst, nmin1, nmin2, mspin, nelec
      integer nstrt1, nlast1, nstrt2, nlast2, nmin12, nminv2, nminv1
      integer nlst21, nlastx
      common /ccntl/ nref,nnst,nmin1,nmin2,mspin,nelec,
     +               nstrt1,nlast1,nstrt2,nlast2,nmin12,nminv2,nminv1,
     +               nlst21,nlastx
c
c
      integer nint, next, ninnex, nwm1, nw, nwp1, nintwo
      integer nint21, nintp1
      common /orb/ nint,next,ninnex,nwm1,nw,nwp1,nintwo,
     +             nint21,nintp1
c
c
      integer nconst, ncond, maxpa, maxpat, maxsin, maxepa
      integer mxdvpa, mxddpa, nspips, nspipt, mxsspa, mxttpa
      integer mxstpa, nsitot, ntrtot
      integer mx12pa, mxvvpa, mxv2pa, maxpsq, max12
      common /auxh/ nconst(8),ncond(8),maxpa(4),maxpat,maxsin(3),
     +              maxepa(4),mxdvpa,mxddpa,nspips(8),nspipt(8),
     +              mxsspa,mxttpa,mxstpa,nsitot(8),ntrtot(8),
     +              mx12pa,mxvvpa,mxv2pa,maxpsq,max12
c
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
      common/mr_mp2l/ ivoc(62),nvoc,idoc(62),ndoc,kvcc(4),nvcc(4)
c
      real*8 critcep, crycep
      integer ifcepa, icepa , lcepaul, lsinsh, itcepa, lprcep, indcep
      integer nrfcsf, mcycep, ibrst,icorr_mr,mrdoc
      common/ccepa/critcep,crycep,ifcepa,icepa,lcepaul,lsinsh,itcepa,
     1                lprcep,indcep,nrfcsf,mcycep,ibrst,icorr_mr,mrdoc
c
c...    only partially ....
c
c
      integer nuse, ngot
      common /corctl/ nuse, ngot
c
c
      logical twice
c
c...  ic does numbering of excited states for print
c...  icd does the same for aa-diagonals
c
      save ic,icd
      data ic,icd/1,1/
      ncst = ncmph(1)
      nstate = ncmph(2)
c
      if (noc4.gt.0.and.noc2.gt.0.and.
     *    ioc4(noc4).gt.ioc2(noc2).and.ioc4(noc4).ne.nint+2)
     *    call caserr('unexpected mppro-call')
c
c...  check for ioc1 ne ioc3 => then do doubles twice
c
      twice = .false.
      if (noc1.ne.noc3) twice = .true.
      if (.not.twice) then
         do 5 i=1,noc1
5        if (ioc1(i).ne.ioc3(i)) twice = .true.
      end if
c
      if (iprimp.ge.500) then
         write(6,601)
601      format(/' -- generate excitations  --')
         if (noc1.gt.0) write(6,602) ioc1
602      format(3x,' annihilate ',(t18,15i4))
         if (noc2.gt.0) write(6,603) ioc2
603      format(3x,' **create** ',(t18,15i4))
         if (noc3.gt.0) write(6,604) ioc3
604      format(3x,' then annihilate ',(t23,15i4))
         if (noc4.gt.0) write(6,605) ioc4
605      format(3x,'  and **create** ',(t23,15i4))
         if (.not.osin) write(6,607)
607      format(3x,' -- ignore singles --')
         if (twice) then
            if (noc3.gt.0) write(6,602) ioc3
            if (noc2.gt.0) write(6,603) ioc2
            if (noc1.gt.0) write(6,604) ioc1
            if (noc4.gt.0) write(6,605) ioc4
         end if
      end if
c
c...   maxexc = max # excitations
c
      maxexc = noc1*noc2 * (1 + noc3*noc4) + iref
      if (twice) maxexc = maxexc + noc1*noc2*noc3*noc4
c
      kch = nuse
      nref1 = nref+1
      nvac = max(maxpa(1)+maxpa(2),maxpat)*nref
      krconf = kch + maxpat**2
      krind = krconf + nref*5
      krvec = krind + nref1*3/lenwrd() + 1
      kind = krvec + nvac
      kv = kind + 3*(nref1*maxexc)/lenwrd()+1
      kvs = kv + iref*nvac
      nuse = kvs + nvac*noc1*noc2
c..     we used quite save estimates  (i hope)
      call corfai('mp')
c...  switch kind to integer addressing
      kind = kind*lenwrd()
      kinds = kind + iref*nref1*3
c
c...  clear krconf with integer 0 
c
      call setsto(nref*5*2,0,cr(krconf+1))
c
c...   if needed get reference vector as first state
      if (iref.ne.0) then
         if (iref.ne.1) call caserr(' wrong iref in mppro')
         call dcopy(nrfcsf,ref,1,cr(kv+1),1)
         call icopy(nref1*3,indr,1,icr(kind+1),1)
      end if
c...       loop over  singles
c...     store in ioksi if we produced an allowed state
      iiex = 1
      nsvec = 0
      kvv = kvs
      kin = kinds
c
      if (osin) then
        do 10 i=1,noc1
          do 10 j=1,noc2
             call mpsing(conf,conf,icst,ncst,ref,indr,nref,cr(kch+1),
     *                   cr(kvv+1),nvac,icr(kin+1),nref1,ioksi(iiex),
     *                   ioc1(i),1,ioc2(j),1)
             nsvec = nsvec + ioksi(iiex)
             kin = kin + ioksi(iiex)*nref1*3
             kvv = kvv + ioksi(iiex)*nvac
             iiex = iiex + 1
10      continue
c
      else
c...   all singles are forbidden
        call setsto(noc1*noc2,0,ioksi)
      end if
c
c...   print
c
      if (nsvec.gt.0.and.iprimp.ge.500)
     *  call mpexpr(iref+1,cr(kvs+1),nvac,icr(kinds+1),nref1,
     *              ioc1,noc1,ioc2,noc2,ioksi,nsvec,'s')
c
c...       doubles // keep on using kin/kvv as addresses
c
      ndvec = 0
c
c...    loop over all singles (allowed or not)
c
      kks = 0
      ifk = 0
      do 20 i=1,noc1
         do 20 j=1,noc2
            ifk = ifk + 1
            nuse = kvv + nvac*noc3*noc4
            call corfai('mp')
            if (ioksi(ifk).eq.1) then
c
c...    symmetry allowed single // is there
c
               call mpsing(conf,conf,icst,ncst,cr(kvs+kks*nvac+1),
     *                     icr(kinds+kks*nref1*3+1),
     *                     nref,cr(kch+1),cr(kvv+1),nvac,
     *                     icr(kin+1),nref1,nn,ioc3,noc3,ioc4,noc4)
               kks = kks + 1
c
            else
c
c...     symmetry forbidden single // make it
c
               call gensin(conf,ref,indr,nref,cr(kch+1),
     *                     cr(krconf+1),nn,cr(krind+1),nref1,
     *                     cr(krvec+1),nvac,ioc1(i),ioc2(j))
c
               if (nn.eq.0) go to 20
               if (iprimp.ge.900) then
c...    print   single just made
                  write(6,610) ioc1(i),ioc2(j)
610               format(/' ..orphaned-single (',i2,' => ',i2,')..')
                  do 15 l=1,nn
15                call prcon(cr(krconf+(l-1)*5+1),l,0)
                  call mpexpr(1,cr(krvec+1),nvac,
     *                        cr(krind+1),nref1,
     *                        ioc1(i),1,ioc2(j),1,1,1,'s')
               end if
c
c        now do doubles from this one
c
               call mpsing(cr(krconf+1),conf,icst,ncst,cr(krvec+1),
     *                     cr(krind+1),
     *                     nref,cr(kch+1),cr(kvv+1),nvac,
     *                     icr(kin+1),nref1,nn,ioc3,noc3,ioc4,noc4)
c
            endif
c
            if (iprimp.ge.500)
     *      call mpexpr(nsvec+iref+ndvec+1,cr(kvv+1),nvac,
     *                  icr(kin+1),nref1,
     *                  ioc1(i),1,ioc2(j),1,nn,nn,'d')
c
            ndvec = ndvec + nn
            kvv = kvv + nn*nvac
            kin = kin + nn*nref1*3
20    continue
c
c...   if twice do doubles for ioc1 and ioc3 interchanged
c
      if (.not.twice) go to 100
c
c...    loop over all singles from ioc3 (none are stored)
c
      do 40 i=1,noc3
         do 40 j=1,noc2
            nuse = kvv + nvac*noc1*noc4
            call corfai('mp')
c
            call gensin(conf,ref,indr,nref,cr(kch+1),
     *                  cr(krconf+1),nn,cr(krind+1),nref1,
     *                  cr(krvec+1),nvac,ioc3(i),ioc2(j))
c
            if (nn.eq.0) go to 40
            if (iprimp.ge.900) then
c...    print   single just made
               write(6,611) ioc3(i),ioc2(j)
611            format(/' ..ignored-single (',i2,' => ',i2,')..')
               do 35 l=1,nn
35             call prcon(cr(krconf+(l-1)*5+1),l,0)
               call mpexpr(1,cr(krvec+1),nvac,
     *                     cr(krind+1),nref1,
     *                     ioc3(i),1,ioc2(j),1,1,1,'s')
            end if
c
c        now do doubles from this one
c
            call mpsing(cr(krconf+1),conf,icst,ncst,cr(krvec+1),
     *                  cr(krind+1),
     *                  nref,cr(kch+1),cr(kvv+1),nvac,
     *                  icr(kin+1),nref1,nn,ioc1,noc1,ioc4,noc4)
c
            if (iprimp.ge.500)
     *      call mpexpr(nsvec+iref+ndvec+1,cr(kvv+1),nvac,
     *                  icr(kin+1),nref1,
     *                  ioc3(i),1,ioc2(j),1,nn,nn,'d')
c
            ndvec = ndvec + nn
            kvv = kvv + nn*nvac
            kin = kin + nn*nref1*3
40    continue
c
100   nuse = kvv
c
c...    we now have nsvec singles and ndvec doubles (and possibly 1 ref)
c...    (nstate = maximum # of configs in space
c
c...      schmidt
      nn = nsvec+ndvec+iref
      nv = nn
      if (nn.gt.maxexc) call caserr('silly in mppro-exc')
      kkept = nuse
      kvr = kkept + nn/lenwrd()+1
      nuse = kvr + nn*nstate
      call corfai('mp')
c
c...  check if excitations are actually symmetry allowed
c
      isym = conf(4,icst)
      if (icst.ge.nstrt1) then
         if (icst.ge.nstrt2) then
            nn = norbtr(isym)
         else
            nn = norbe(isym)
         endif
         if (nn.le.0) then
            if (iprimp.ge.1000) write(6,613) isym
613         format(' == eliminated states for symmetry',i2,' ==')
            do 110 i=6,10
110         ncmph(i) = 0
            go to 120
         end if
      end if
c
      call mpfpr(conf,cr(kv+1),nvac,icr(kind+1),nref1,
     *           cr(kvr+1),nstate,iadv,iadvt,cr(kkept+1),
     *           iref,nsvec,ndvec,mode,vcont,ncmph,ic,
     *           cr)
c
      if (iprimp.ge.2000) then
         write(6,606) ic,ic+ncmph(10)-1,
     *                ncmph(6),ncmph(7),ncmph(8),nstate
606      format(/' --',i5,' ..',i5,' -- # ref',i2,
     *           ' # sing',i3,' # doub',i3,' dim',i3)
         if (ncmph(9).gt.0) call prvc4(vcont,nstate,ncmph(9),1,1,0,10)
         ic = ic + ncmph(10)
      end if
c
c...   for nm2 states do aa-states seperately ( no external triplets)
c...   only symmetry 1 has aa states !!
c

 120   if (icst.lt.nstrt2.or.conf(4,icst).ne.1) return
c
c...   check symmetry, just in case
c
      if (norbsi(1).le.0) then
         write(6,613) 11
         do 130 i=6,10
130      ncmphd(i) = 0
         call caserr(' nonexisting aa in mppro')
      end if
c
      if (iref.gt.0.or.nsvec.gt.0.or.mode.eq.1)
     *     call caserr('strange nm2')
c
      call mpfpr(conf,cr(kv+1),nvac,icr(kind+1),nref1,
     *           cr(kvr+1),ncmphd(2),iadv,iadvt,cr(kkept+1),
     *           iref,nsvec,ndvec,mode,vcontd,ncmphd,ic,
     *           cr)
c
c...  that's it for this set
c
      if (iprimp.ge.2000) then
         write(6,612)
612      format(' // n-2 diagonals \\ ')
         write(6,606) icd,icd+ncmphd(9)-1,
     *                ncmphd(6),ncmphd(7),ncmphd(8),ncmphd(2)
         call prvc4(vcontd,ncmphd(2),ncmphd(9),1,1,0,10)
         icd = icd + ncmphd(10)
      end if
c
      return
      end
      subroutine mpexpr(kstart,v,nvac,ind,nref1,
     *                  ioc1,noc1,ioc2,noc2,ioksi,nsvec,type)
c
      implicit real*8 (a-h,o-z)
c
c...   print excitations obtained
c...  type is either 'r' (reference), 's' (singles), 'd' (doubles)
c
      dimension v(nvac,nsvec),ind(3,nref1,nsvec)
      dimension ioc1(noc1),ioc2(noc2),ioksi(*)
      character*1 type
c
      write(6,'(1x)')
c
      if (type.eq.'r') then
         write(6,101)
101      format(' ====  REFERENCE  ====== ')
         call mp1pr(ind(1,1,1),nref1,v(1,1),0,0,0)
      else if (type.eq.'s') then
         write(6,201) nsvec
201      format(' ==== ',i5,' SINGLES ======')
         kks = 0
         kioksi = 0
         do 250 i=1,noc1
           do 250 j=1,noc2
             kioksi = kioksi + 1
             if (ioksi(kioksi).eq.0) go to 250
             kks = kks + 1
             call mp1pr(ind(1,1,kks),nref1,v(1,kks),kstart+kks-1,
     *                  ioc1(i),ioc2(j))
250      continue
c
      else if (type.eq.'d') then
         write(6,301) nsvec,ioc1(1),ioc2(1)
301      format(' ==== ',i5,' DOUBLES ====== (',i3,' =>',i3,
     *          '  // *=>* )')
         do 350 i=1,nsvec
            call mp1pr(ind(1,1,i),nref1,v(1,i),kstart+i-1,0,0)
350      continue
c
      else
         call caserr(' wrong call to mpexpr')
      end if
c
      return
      end
      subroutine mp1pr(ind,nind,v,num,ioc1,ioc2)
c
c...  print a single excited state
c...  called by mpexpr
c
      implicit real*8 (a-h,o-z), integer(i-n)
c
      dimension ind(3,nind),v(*)
      logical first
      dimension ip(5),rp(5)
c
      first = .true.
      n = 0
      do 240 k=1,nind
         if (ind(1,k).eq.0) go to 240
         do 220 kk=1,ind(3,k)
            n = n + 1
            ip(n) = ind(1,k)
            rp(n) = v(ind(2,k)-1+kk)
            if (n.eq.5) then
               if (first) then
                  if (ioc1.ne.0) then
                     write(6,601) num,ioc1,ioc2,
     *                            (ip(l),rp(l),l=1,n)
                  else
                     write(6,603) num,(ip(l),rp(l),l=1,n)
                  end if
                  first = .false.
               else
                  write(6,602) (ip(l),rp(l),l=1,n)
               end if
               n = 0
            end if
220      continue
240   continue
c
      if (n.ne.0) then
         if (first) then
            if (ioc1.ne.0) then
               write(6,601) num,ioc1,ioc2,
     *                      (ip(l),rp(l),l=1,n)
            else
               write(6,603) num,(ip(l),rp(l),l=1,n)
            end if
         else
            write(6,602) (ip(l),rp(l),l=1,n)
         end if
      end if
c
601     format(1x,'#',i4,2x,i3,' =>',i3,3x,5(i4,f7.4))
602     format(17x,3x,5(i4,f7.4))
603     format(1x,'#',i4,12x,3x,5(i4,f7.4))
c
      return
      end
      subroutine mpsing(rconf,conf,icst,ncst,ref,indref,nref,ch,svec,
     *                  nvac,insvec,nref1,nsvec,ifr,nfr,ito,nto)
c
c     Arkansas MR-MP  (1989)
c     sort out  states for mp2/3 calculation
c     produce single excitations // we only need 1-electron interactions
c     nsvec vectors for single excitations returned in svec (packed)
c     **note** for n-2 singlet/triplet are consecutive/configuration
c
c      rconf : occupation-patterns for reference configs
c      conf  : occupation-patterns for excited configs (5,*)
c              icst : start of conf range / ncst # conf
c      ref   : reference CI-vector (packed) indref : index / nref #confs
c      ch    : scratch space to hold 1-electron couplings
c              length :maxpat**2
c      svec  : single-excitation vectors (result)  (packed)
c      nvac  : max. length needed/vector (say nref*maxpat)
c      insvec: index vectors to svec (#conf / position /length in svec)
c            ; last insvec(1,..) = 0
c      nsvec : # vectors produced
c      ifr(nfr) : orbitals to excite from
c      ito(nto) : orbitals to excite to (may be identical to ifr)
c
c..    we need seperate icre/iann routines since no restrictions and
c..    conf(3-5,..) is already in use (# spin-paths/ classification)
c
c...  **note** since i=>i excitation is allowed nd=0 can occur
c
      implicit real*8 (a-h,o-z)
      integer conf,rconf,cannil
c
      dimension ch(*),rconf(5,*),conf(5,*),ref(*),indref(3,nref)
      dimension ifr(nfr),ito(nto)
      dimension svec(nvac,*),insvec(3,nref1,*)
c
c==== vacuum - internal excitation interactions =======only==========
c
      common/symcon/norbbb(7),nd,ii,jj,kk,ll,konte(14)
c
      dimension cannil(5),ccreat(5)
c     dimension ior(4)
chvd
c...  skipped : stores is a configuration has been skipped before
c...            reason is we don't want to print thousands of
c...            warnings. 
c...            this hack is required for using restrict in 
c...            contracted CI calculations.
      logical skipped
      save skipped
      data skipped/.false./
chvd
c
      nsvec = 0
c
c...   loop over "from" excitations
c
      do 10000 ifrom = 1,nfr
c
c...   loop over "to" excitations
c
         do 1000 itorb=1,nto
c
            nsvec = nsvec + 1
            call vfill(0.0d0,svec(1,nsvec),1,nvac)
            call setsto(nref1*3,0,insvec(1,1,nsvec))
            kks = 1
            kki = 1
c
c..     loop over reference configurations
c..     reference configs are not consecutive, neither are coefficients
c
            do 500 iref=1,nref
               ir = indref(1,iref)
               if (ir.eq.0) go to 500
               if (iannmp(rconf(1,ir),ifr(ifrom),cannil).eq.0) go to 500
               if (icremp(cannil,ito(itorb),ccreat).eq.0) go to 500
               kref = indref(2,iref)
c
c...   look for the configuration
c...   if it's not there the whole thing is symmetry forbidden
c
               ic = ichmp(conf(1,icst),ncst,ccreat,5) + icst-1
               if (ic.eq.icst-1) then
c
c...   if not found check if it's because of impossible spin
c...    (if so do'nt quit)
c
                  nn = nspip(ccreat,1)
                  if (nn.le.0) go to 500
c
c...   now it must be symmetry forbidden so it must be first one
c
                  nsvec = nsvec - 1
                  if (kki.gt.1.and..not.skipped) then
                   write(6,601) ifr(ifrom),ito(itorb)
601    format(/' a configuration was not found for excitation ',
     *         i3,' =>',i3,/,' symmetry is probably not the reason',
     *         '  (if you used restrict this is okay)',
     *         ', so check configuration generation ',
     *         /' the missing configuration is .... ')
                   call prcon(ccreat,kki,0)
c                  call caserr('symmetry elimination not at first conf')
c ??? conflicts with restrict option !!!
c     so we skip the configuration and the crash.
                   skipped = .true.
                   goto 500
                  else if (kki.gt.1.and.skipped) then
                   goto 500
                  end if
                  go to 1000
               end if
c
c...   we found both configurations E(ifr=>ito) conf(ir) = conf(ic)
c
               call expan1(rconf(1,ir))
               call upack2(rconf(3,ir),ncolt,ncols)
c
               call expan2(conf(1,ic))
               call upack2(conf(3,ic),nrowt,nrows)
c
c..   store the single excitation
c
               insvec(1,kki,nsvec) = ic
               insvec(2,kki,nsvec) = kks
               insvec(3,kki,nsvec) = nrows+nrowt
               kki = kki + 1
c
c.... check  number of differences
c
               if (nd.eq.0) then
c
c...   i => i excitation (only vacuum-vacuum // nvoc-nvoc)
c...   set up unit matrix and jump over symb1/sym1tr
c
                  if (ifr(ifrom).ne.ito(itorb))
     *            call caserr(' no difference for no i=> i')
                  if (nrowt.ne.0)
     *            call caserr(' triplet couplings for i=>i  ')
                  call vclr(ch,1,ncols*ncols)
                  do 4 i=1,ncols
4                 ch((i-1)*ncols+i) = 1.0d0
                  go to 5
               end if
c
               if (nd.ne.2) call caserr(' very funny in mpsing ')
c
c...     check orbitals to be absolutely sure (**debug**)
c
               if ((ii.ne.min(ifr(ifrom),ito(itorb))).or.
     *             (jj.ne.max(ifr(ifrom),ito(itorb))))
     *         call caserr('silly excitation failure in mpsing')
c
c.... single difference
c...  coefficients of single : <br/d>  (ch's) (over sigma)
c...  normalization of br : sum <br/d><d/br> = sum  ch.ch
c...  single excitation in ii,jj
c
c     *** note the second conf(ic)  defines most rapid varying index ***
c
               call symb1
c
               if (nrows.eq.0) go to 20
c
               call sym1tr(ch,ifl,0,0)
               if (ifl.eq.0) call caserr('ifl 0 0 = 0 (mpsing)')
c
c..    loop over all reference csf's and store resulting linear comb
c
5              do 10 i=1,ncols
               call daxpy(nrows
     +         ,ref(kref+i-1),ch((i-1)*nrows+1),1,svec(kks,nsvec),1)
10             continue
c
c...  normalise
c
               kks = kks + nrows
c
c...           for (n-2) do triplets as well
c
20             if (nrowt.eq.0) go to 500
c
               call sym1tr(ch,ifl,0,1)
c
c..    loop over all reference csf's and store resulting linear comb
c
               if (ifl.eq.0) call caserr('ifl 0 1 = 0 (mpsing)')
               do 30 i=1,ncols
               call daxpy(nrowt
     +         ,ref(kref+i-1),ch((i-1)*nrowt+1),1,svec(kks,nsvec),1)
30             continue
c
c...  normalise
c
               kks = kks + nrowt
c
500         continue
c
c...  last element of insvec : 0,kks(=length+1),0
c
            insvec(2,nref1,nsvec) = kks
            if (kki.gt.nref1) call caserr('insvec overflow mpsing')
c
c...  if nothing was produced // reset nsvec
c
         if (kki.eq.1) nsvec = nsvec - 1
c
1000     continue
c
10000 continue
c
      return
      end
      subroutine gensin(rconf,ref,indref,nref,ch,
     *                  sconf,ic,inds,nref1,
     *                  svec,nvac,iiii,jjjj)
c
c     Arkansas MR-MP  (1989)
c
c     produce a straight single excitation from iiii to jjjj
c     conf,ref,indref : starting vector
c     sconf,svec,inds : result
c     only used to produce n and n-1 states a / uses full occupations
c
c..    we need seperate icre/iann routines since no restrictions and
c..    conf(3-5,..) is already in use (# spin-paths/ classification)
c
      implicit real*8 (a-h,o-z)
      integer  sconf,rconf,cannil
c
      dimension ch(*),rconf(5,*),sconf(5,*),ref(*),indref(3,nref)
      dimension svec(nvac),inds(3,nref1)
c
c==== vacuum - internal excitation interactions =======only==========
c
c
      integer nconst, ncond, maxpa, maxpat, maxsin, maxepa
      integer mxdvpa, mxddpa, nspips, nspipt, mxsspa, mxttpa
      integer mxstpa, nsitot, ntrtot
      integer mx12pa, mxvvpa, mxv2pa, maxpsq, max12
      common /auxh/ nconst(8),ncond(8),maxpa(4),maxpat,maxsin(3),
     +              maxepa(4),mxdvpa,mxddpa,nspips(8),nspipt(8),
     +              mxsspa,mxttpa,mxstpa,nsitot(8),ntrtot(8),
     +              mx12pa,mxvvpa,mxv2pa,maxpsq,max12
c
      common/symcon/norbbb(7),nd,ii,jj,kk,ll,konte(14)
      common/three/n1l(1024),n2l(1024),minlhs(64),maxlhs(64),iocc(64),
     *n1r(1024),n2r(1024),minrhs(64),maxrhs(64),jocc(64),modelr
c
      dimension cannil(5)
c     dimension ior(4)
c
      call vfill(0.0d0,svec,1,nvac)
      call setsto(nref1*3,0,inds)
      kks = 1
      ic = 1
c
c..  loop over reference configurations
c..  reference configs are not consecutive, neither are coefficients
c
      do 500 iref=1,nref
         ir = indref(1,iref)
         if (ir.eq.0) go to 500
         if (iannmp(rconf(1,ir),iiii,cannil).eq.0) go to 500
         if (icremp(cannil,jjjj,sconf(1,ic)).eq.0) go to 500
         kref = indref(2,iref)
c
c...     establish # spin-possibilities / and check if not over maxpat
c...     if spin forbidden skip
c
         nn = nspip(sconf(1,ic),1)
         if (nn.le.0) go to 500
c
         if (nn.gt.maxpat) call caserr('spin > maxpat in gensin')
         sconf(3,ic) = nn
c
         call expan1(rconf(1,ir))
         call upack2(rconf(3,ir),ncolt,ncols)
c
         call expan2(sconf(1,ic))
         call upack2(sconf(3,ic),nrowt,nrows)
c
c..   store the single excitation
c
         inds(1,ic) = ic
         inds(2,ic) = kks
         inds(3,ic) = nrows+nrowt
         ic = ic + 1
c
c.... check  number of differences
c
         if (nd.eq.0) then
c
c...   i => i excitation (only vacuum-vacuum // nvoc-nvoc)
c...   set up unit matrix and jump over symb1/sym1tr
c
            if (iiii.ne.jjjj)
     *      call caserr(' no difference for no i=> i')
            if (ncols.ne.nrows.or.ncolt.ne.0)
     *      call caserr('spin foul-up in gensin')
            call vclr(ch,1,ncols*ncols)
            do 4 i=1,ncols
4           ch((i-1)*ncols+i) = 1.0d0
            go to 5
         end if
c
         if (nd.ne.2) call caserr('gensin produced double diff')
c
c...     check orbitals to be absolutely sure (**debug**)
         if ((ii.ne.min(iiii,jjjj)).or.
     *       (jj.ne.max(iiii,jjjj)))
     *      call caserr(' silly excitation failure in gensin')
c
c.... single difference
c...  coefficients of single : <br/d>  (ch's) (over sigma)
c...  normalization of br : sum <br/d><d/br> = sum  ch.ch
c...  single excitation in ii,jj
c
c     *** note the second conf(ic)  defines most rapid varying index ***
c
         call symb1
c
         call sym1tr(ch,ifl,0,0)
         if (ifl.eq.0) call caserr('ifl 0 in gensin')
c
c..    loop over all reference csf's and store resulting linear comb
c
5       do 10 i=1,ncols
        call daxpy(nrows
     +  ,ref(kref+i-1),ch((i-1)*nrows+1),1,svec(kks),1)
10      continue
c
c...  normalise
c
         kks = kks + nrows
c
500   continue
c
c...  last element of insvec : 0,kks(=length+1),0
c
      inds(2,nref1) = kks
      if (ic.gt.nref1) call caserr('inds overflow gensin')
c...    ic return # confs actually produced
      ic = ic - 1
c
      return
      end
      function iannmp(ref,i,conf)
c*ap*
c***     MR-MP special version that just does the job and no tango
c***     only  occupation pattern is set (cf icremp,ichmp)
c
c     annihilate an electron in orbital i in configuration ref
c     to produce configuration conf, both using nw words
c     conf contains the occupation-scheme (2 bits per orbital)
c     left-justified , to be read from left to right ,
c     while the number of 1-bits gives the occupation (11b = 2)
c     so each word has room for  32 orbitals
c     upon return : iannmp = 1 => success  / iannmp = 0 => failure
      implicit real*8 (a-h,o-z)
      integer pad1,pad2,pad3,conf,ref,pad
      dimension ref(*),conf(*)
      common/dopey/nshif
c
      integer nint, next, ninnex, nwm1, nw, nwp1, nintwo
      integer nint21, nintp1
      common /orb/ nint,next,ninnex,nwm1,nw,nwp1,nintwo,
     +             nint21,nintp1
c
c
      integer n24, magic, n32, n64, n8, n16, n48
      common /cdcryi/ n24,magic,n32,n64,n8,n16,n48
c
c
      common/apollo/pad1,pad2,pad3
c
c.... see in which word and which position  i is
      kc = (i-1)/n32
      ks = (i-kc*n32)*2
      kc = kc + 1
c.... check occupation of i
      iocc = iand(ishftc(ref(kc),ks,64),pad3)
      if (iocc-1) 10,30,20
c.... no electrons in i => failure
10    iannmp = 0
      return
c.... 2 electrons in i
20    pad = pad2
      go to 40
c.... 1 electron in i
30    pad = pad1
40    iannmp = 1
      call dcopy(nwm1,ref,1,conf,1)
      conf(kc) = ieor(conf(kc),ishftc(pad,n64-ks,64))
c
      return
      end
      integer function icremp(ref,i,conf)
c*ap*
c     create an electron in orbital i in configuration conf
c***     MR-MP special version that just does the job and no tango
c***     only  occupation pattern is set (cf iannmp,ichmp)
c
      implicit real*8 (a-h,o-z)
      integer pad1,pad2,pad3,conf,ref,pad
      dimension ref(*),conf(*)
c
      integer nint, next, ninnex, nwm1, nw, nwp1, nintwo
      integer nint21, nintp1
      common /orb/ nint,next,ninnex,nwm1,nw,nwp1,nintwo,
     +             nint21,nintp1
c
c
      integer n24, magic, n32, n64, n8, n16, n48
      common /cdcryi/ n24,magic,n32,n64,n8,n16,n48
c
      common/apollo/pad1,pad2,pad3
c
c.... determine which word and which position
      kc = (i-1)/32
      ks = (i-kc*32)*2
      kc = kc + 1
c.... check occupation of i
      iocc = iand(ishftc(ref(kc),ks,64),pad3)
      if (iocc-1) 20,30,10
c.... 2 electrons in i => failure
10    icremp = 0
      return
c.... no electrons in orbital i
20    pad = pad1
      go to 40
c.... 1 electron in orbital i
30    pad = pad2
40    icremp = 1
      call dcopy(nwm1,ref,1,conf,1)
      conf(kc) = ieor(conf(kc),ishftc(pad,n64-ks,64))
c
      return
      end
      function ichmp(ref,nref,conf,nword)
c     compare the  occ.-scheme words of conf to those
c     of the reference set to see which one it is
c     if so : ** return the number in ichmp **  // else return 0
c***     MR-MP special version that just does the job and no tango
c***      (cf icremp,iannmp)
      implicit real*8 (a-h,o-z)
      integer ref,conf
      dimension ref(nword,*),conf(*)
c
      integer nint, next, ninnex, nwm1, nw, nwp1, nintwo
      integer nint21, nintp1
      common /orb/ nint,next,ninnex,nwm1,nw,nwp1,nintwo,
     +             nint21,nintp1
c
c
      if (nwm1.eq.1) then
         ichmp = locate(ref(1,1),nref,nword,conf(1))
      else
         k=0
         jr=0
4        ir=locate(ref(1,k+1),nref-k,nword,conf(1))
         if (ir.eq.0) go to 25
3        if(ir.eq.jr)goto 31
         k=k+ir-1
         ir=locate(ref(2,k+1),nref-k,nword,conf(2))
         if (ir.eq.0) go to 25
         if (ir.eq.1) goto 31
         k=k+ir-1
         jr=1
         goto 4
c.... found identical configuration
31       ichmp = ir+k
      end if
c
      return
c         .. no match
25    ichmp = 0
c
      return
      end
      subroutine mpfpr(conf,v,lv,ind,nind,vr,lvr,iadv,iadvt,kept,
     *                 nref,nsing,ndoub,mode,vcont,ncmph,ic,cr)
c
cmp   changed to contraction
c     produce real orthogonalised excited vectors
c     n-2 : singlets first / triplets next
c           for aa- lvr is # singlets => triplets are ignored
c
      implicit real*8 (a-h,o-z)
c
      dimension v(lv,*),ind(3,nind,*),vr(lvr,*),iadv(*),iadvt(*)
      dimension kept(*),vcont(lvr,*),ncmph(10)
      dimension conf(5,*),cr(*)
c
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
c
c
      integer nuse, ngot
      common /corctl/ nuse, ngot
c
c
      dimension nexmp(3)
c
      nv = nref+nsing+ndoub
c
c...   finish of building of contraction vectors
c...   we have a set of vectors representing singly and
c...   doubly excited states
c...   expand them to full vectors
c
      call  mpexex(conf,v,lv,nv,ind,nind,vr,lvr,iadv,iadvt)
c
c...    print
c
      if (iprimp.ge.2000) then
20       write(6,600)
600      format(/' == (too many)  excited states  ==')
         call prvc4 (vr,lvr,nv,dum,dum,0,10)
      end if
c
c
c...  Orthogonalise states according to mode:
c
c
c...  The projection operators for the Pulay models 1, 2, and 3
c...  are defined such that  I = P0 + P1 + P2  .
c...  Where  P0  projects onto the reference function, P1 projects onto
c...  the space of singly substituted states, and P2 projects onto the
c...  the space of doubly excited states.
c
c...  The projection operators for the Andersson models 11, 12, and 13
c...  are defined such that  I = P0 + P'k + P'1 + P'2  .
c...  Where the projection operators basically mean the same as in the
c...  Pulay models, except that the singles and doubles in the 
c...  orthogonal complement of the reference function in the CAS 
c...  reference space are contained in  P'k  . So  P'1  and  P'2  lack
c...  the internal excited stated.
c
c...  The states in  P'k  are assumed not to contribute to the MRMP2/3
c...  energies and are therefore removed in models 11, 12, and 13.
c
      if (mode.eq.1) then
c
c...     H0 = |P0>F<P0| + |P1>F<P1| + |P2>F<P2|
c
         call ibasgn(nv,1,0,kept)
      else if (mode.eq.2) then
c
c...     H0 = |P0>F<P0| + |P2>F<P2| (i.e. singles eliminated)
c
         call ibasgn(nv,1,0,kept)
      else if (mode.eq.3) then
c
c...     H0 = |P0>F<P0| + |P1+2>F<P1+2|
c
         call ibasgn(nv,1,0,kept)
         ndoub = ndoub + nsing
         nsing = 0
      else if (mode.eq.11) then
c
c...     H0 = |P0>F<P0| + |P'1>F<P'1| + |P'2>F<P'2|
c
         call ibasgn(nv,1,0,kept)
         if (nref.ne.0) then
c...        These are the internal excited states
            call ibasgn(nsing+ndoub,0,0,kept(nref+1))
         endif
      else if (mode.eq.12) then
c
c...     H0 = |P0>F<P0| + |P'2>F<P'2| (i.e. singles eliminated)
c
         call ibasgn(nv,1,0,kept)
         if (nref.ne.0) then
c...        These are the internal excited states
            call ibasgn(ndoub+nsing,0,0,kept(nref+1))
         endif
      else if (mode.eq.13) then
c
c...     H0 = |P0>F<P0| + |P'1+2>F<P'1+2|
c
         call ibasgn(nv,1,0,kept)
         ndoub = ndoub + nsing
         nsing = 0
         if (nref.ne.0) then
c...        These are the internal excited states
            call ibasgn(ndoub,0,0,kept(nref+1))
         endif
      else
         call caserr(' wrong mode in mpfpr ')
      endif
c
c
c...  And now for the actual orthogonalisation
c
c
      nused = nuse
      if (iortho.eq.1) then
c
c...    schmidt option (Modified Gramm-Schmidt)
c
        call mpschm(1,nv,1,nv,vr,lvr,nv,kept,cortho,iprimp)
      else if (iortho.eq.2) then
c
c...    schmidtp option (Pivoted Repeated Modified Gramm-Schmidt)
c
        i10  = nuse + 1
        i20  = i10  + nv
        nuse = i20  + nv
        call corfai('mp')
        call prmpschm(nref,nsing,ndoub,vr,lvr,nv,kept,cortho,iprimp,
     *                cr(i10),cr(i20))
      else if (iortho.eq.3) then
c
c...    schmidt4 option (Quadruple Precision Modified Gramm-Schmidt)
c
        i10  = nuse + 1
        nuse = i10  + lvr * nv * 2
        call corfai('mp')
        call caserr('no quadruple precision schmidt available')
      else if (iortho.eq.4) then
c
c...    lowdin option (Lowdin)
c
        i00  = max(nref,nsing,ndoub)
        i10  = nuse + 1
        i20  = i10  + i00 * (i00 + 1) / 2
        i30  = i20  + i00 * (i00 + 1)
        i40  = i30  + lvr * i00 * 2
        j10  = nuse + 1
        nuse = max(i40  + i00 * 2 / lenwrd(), j10 + lvr * nv * 2)
        call corfai('mp')
        call mpsymm(nref,nsing,ndoub,vr,lvr,nv,kept,cortho,iprimp,
     *              cr(i10),cr(i20),cr(i30),cr(i40),cr(j10))
      else if (iortho.eq.5) then
c
c...    house option (Householder QR Factorisation)
c
        i10  = nuse + 1
        i20  = i10  + lvr * lvr
        i30  = i20  + lvr
        i40  = i30  + max(lvr,nv)
        i50  = i40  + nv
        nuse = i50  + (nv + 1) / lenwrd()
        call corfai('mp')
        call mphouse(nref,nsing,ndoub,vr,lvr,nv,kept,cortho,iprimp,
     *               cr(i10),cr(i20),cr(i30),cr(i40),cr(i50))
      else
        call caserr(' strange orthogonalisation option')
      endif
      nuse = nused
c
c...  Fix model 2 after orthogonalisation. This can't be done
c...  before orthogonalisation because the doubles are not yet
c...  uniquely defined.
c
      if (mode.eq.1) then
c
c...     H0 = |P0>F<P0| + |P1>F<P1| + |P2>F<P2|
c
      else if (mode.eq.2) then
c
c...     H0 = |P0>F<P0| + |P2>F<P2| (i.e. singles eliminated)
c
         call ibasgn(nsing,0,0,kept(nref+1))
      else if (mode.eq.3) then
c
c...     H0 = |P0>F<P0| + |P1+2>F<P1+2|
c
      else if (mode.eq.11) then
c
c...     H0 = |P0>F<P0| + |P'1>F<P'1| + |P'2>F<P'2|
c
      else if (mode.eq.12) then
c
c...     H0 = |P0>F<P0| + |P'2>F<P'2| (i.e. singles eliminated)
c
         call ibasgn(nsing,0,0,kept(nref+1))
      else if (mode.eq.13) then
c
c...     H0 = |P0>F<P0| + |P'1+2>F<P'1+2|
c
      else
         call caserr(' wrong mode in mpfpr ')
      endif
c
c
c...  get final vectors neatly alligned
c...  fill info-array nexmp(3) with : # ref, # sing, # doub
c...  a complete-set => unit matrix
c
      call mppact(vr,lvr,kept,nref,nexmp(1),vcont)
      call mppact(vr(1,nref+1),lvr,kept(nref+1),nsing,nexmp(2),
     *            vcont(1,nexmp(1)+1))
      call mppact(vr(1,nref+nsing+1),lvr,kept(nref+nsing+1),ndoub,
     *            nexmp(3),vcont(1,nexmp(1)+nexmp(2)+1))
c
c...  check
c
      if (nexmp(1)+nexmp(2)+nexmp(3).gt.lvr)
     *   call caserr(' dimensionality fault in mpfpro ')
c
c...  if all states are in one set => use csf's
c
      if (nexmp(2).eq.lvr.or.nexmp(3).eq.lvr) then
         call vclr(vcont,1,lvr*lvr)
         do 40 i=1,lvr
40       vcont(i,i) = 1.0d0
      end if
c
c...  either move nexmp info to ncmph or check
c
c     write(6,*)'nexmp = ', nexmp
c     write(6,*)'ncmph = ', ncmph
      do 50 i=1,3
         if (ncmph(5+i).lt.0) then
            ncmph(5+i) = nexmp(i)
         else
            if (ncmph(5+i).ne.nexmp(i)) call caserr('nexmp in mpfpr')
         end if
50    continue
c
c...  ncmph(9) counts total # excited states (for dimensioning)
c...  ncmph(10) counts total # excited configs
c
      ncmph(9) = ncmph(6)+ncmph(7)+ncmph(8)
      ncmph(10) = min(ncmph(6),1)+min(ncmph(7),1)+min(ncmph(8),1)
c
      return
      end
      subroutine mpexex(conf,v,lv,nv,ind,nind,vr,lvr,iadv,iadvt)
c
      implicit real*8 (a-h,o-z), integer (i-n)
      integer conf
      dimension conf(5,*)
c
      dimension v(lv,nv),ind(3,nind,nv),vr(lvr,nv),iadv(*),iadvt(*)
c
c...  expand a packed vector to full form using addressing
c...  in iadv,iadvt
c...  first all singlets then all triplets
c...  if a triplet address is beyond lvr => triplets are ignored
c...  this is used to remove triplets in processing the diagonals.
c...
c...  The reasoning is: 
c...  Singlets and triplets are processed by the same code, however
c...  triplets don't have diagonals in contrast to the singlets. 
c...  Therefore when the diagonals are being processed the length
c...  of the vectors <lvr> is set to the number of singlets which 
c...  are in front of the triplets. 
c
c     write(6,*)'mpexex: lvr, nv = ', lvr,nv
      call vclr(vr,1,lvr*nv)
c
      do 100 iv=1,nv
         do 50 j=1,nind
            iconf = ind(1,j,iv)
            if (iconf.eq.0) go to 100
            call upack2(conf(3,iconf),np3,np1)
            if (np1.gt.0) then
               ifr = ind(2,j,iv)
               ito = iadv(iconf)+1
               call dcopy(np1,v(ifr,iv),1,vr(ito,iv),1)
            end if
            if (np3.gt.0.and.iadvt(iconf).lt.lvr) then
               ifr = ind(2,j,iv) + np1
               ito = iadvt(iconf)+1
               call dcopy(np3,v(ifr,iv),1,vr(ito,iv),1)
            end if
50      continue
100   continue
c
      return
      end
      subroutine mppact(vr,lvr,kept,nold,nnew,v)
c
c..   produce a compact set of vectors in v from the vectors in vr
c..   according to kept (1=keep/0=loose)
c..   return # kept in nnew
c
      implicit real*8 (a-h,o-z)
c
      dimension vr(lvr,*),kept(*),v(lvr,*)
c
      nnew = 0
      do 10 i=1,nold
         if (kept(i).eq.1) then
            nnew = nnew + 1
            call dcopy(lvr,vr(1,i),1,v(1,nnew),1)
         end if
10    continue
c
      return
      end
      subroutine cmpfij(conf,istart,ic,iadv,iadvt,vcc,ncmph,nirr,nir1,
     *                  fock,result,x,y,temp,final)
c
      implicit real*8 (a-h,o-z)
c
      integer conf
      dimension conf(5,*)
      dimension iadv(*),iadvt(*),vcc(*),ncmph(10,18,nir1)
      dimension fock(*),result(*),x(*),y(*),final(*)
c
c...  routine to handle symmetry/doc blocked calling of mpfij
c...  (cf mpsrcc for blocking conventions)
c...  istart is start of conf's to be treated
c...  icc is start of excited states
c
      logical inclass
c
      mlast = istart - 1
      kvcc = 1
      do 100 ijrr=1,nir1
         nnlast = mlast
         irr = ic
         kvcrr = kvcc
         do 90 is=1,2+nirr*2
            nconi = ncmph(1,is,ijrr)
            ncsfi = ncmph(2,is,ijrr)
            nseti = ncmph(3,is,ijrr)
            do 80 iset=1,nseti
               mstart = mlast + 1
               mlast = mlast + nconi
               nlast = nnlast
               ir = irr
               kvcr = kvcrr
               do 70 js=1,is
                  nconj = ncmph(1,js,ijrr)
                  ncsfj = ncmph(2,js,ijrr)
                  nsetj = ncmph(3,js,ijrr)
                  if (js.eq.is) nsetj = iset
                  do 60 jset=1,nsetj
                    nstart = nlast + 1
                    nlast = nlast + nconj
                    if (inclass(is,js,mstart,mlast,nstart,nlast))
     *              call mpfij(conf,mstart,mlast,nstart,nlast,
     *                         ic,ir,iadv,iadvt,
     *                         vcc(kvcc),ncmph(1,is,ijrr),
     *                         vcc(kvcr),ncmph(1,js,ijrr),   fock,
     *                         result,ncsfi,ncsfj,x,y,temp,final)
                    kvcr = kvcr + ncmph(9,js,ijrr)*ncsfj
                    ir = ir + ncmph(10,js,ijrr)
60                continue
70             continue
               kvcc = kvcc + ncmph(9,is,ijrr)*ncsfi
               ic = ic + ncmph(10,is,ijrr)
80          continue
90       continue
100   continue
c
      return
      end
      subroutine cmpfia(conf,istarc,istarr,icc,irr,iadv,iadvt,
     *                  vcc,ncmphc,vcr,ncmphr,nirr,nirrc,nirrr,
     *                  result,cx,temp,final,
     *                  npp,ipair,mxpair)
c
      implicit real*8 (a-h,o-z)
c
      integer conf
      dimension conf(5,*),iadv(*),iadvt(*)
      dimension vcc(*),vcr(*)
      dimension ncmphc(10,18,nirrc),ncmphr(10,18,nirrr)
      dimension result(*),cx(*),final(*)
      dimension npp(*),ipair(2,mxpair,*)
c
c...  routine to handle symmetry/doc blocked calling of mpfia
c...  (cf mpsrcc for blocking conventions)
c...  istart/jstart are start of conf's to be treated
c
      logical inclass
c
      mlast = istarc - 1
      kvcc = 1
      ic = icc
c
c...  loop over ic (doublets/n-2)
c
      do 100 iirr=1,nirrc
         do 90 is=1,2+nirr*2
            nconi = ncmphc(1,is,iirr)
            ncsfi = ncmphc(2,is,iirr)
            nseti = ncmphc(3,is,iirr)
            do 80 iset=1,nseti
               mstart = mlast + 1
               mlast = mlast + nconi
c
               nlast = istarr - 1
               kvcr = 1
               ir = irr
c
c...          loop over ir (vacuum/doublets)
c
               do 75 jjrr=1,nirrr
                  do 70 js=1,2+nirr*2
                     nconj = ncmphr(1,js,jjrr)
                     ncsfj = ncmphr(2,js,jjrr)
                     nsetj = ncmphr(3,js,jjrr)
                     do 60 jset=1,nsetj
                        nstart = nlast + 1
                        nlast = nlast + nconj
                        if (inclass(is,js,mstart,mlast,nstart,nlast))
     *                  call mpfia(conf,mstart,mlast,nstart,nlast,
     *                             ic,ir,iadv,iadvt,
     *                             vcc(kvcc),ncmphc(1,is,iirr),
     *                             vcr(kvcr),ncmphr(1,js,jjrr),
     *                             result,ncsfi,ncsfj,cx,temp,final,
     *                             npp,ipair,mxpair)
                        kvcr = kvcr + ncmphr(9,js,jjrr)*ncsfj
                        ir = ir + ncmphr(10,js,jjrr)
60                   continue
70                continue
75             continue
               kvcc = kvcc + ncmphc(9,is,iirr)*ncsfi
               ic = ic + ncmphc(10,is,iirr)
80          continue
90       continue
100   continue
c
      return
      end
      subroutine cmpfab(icc,vcc,ncmph,nirr,nirr1,final)
c
      implicit real*8 (a-h,o-z)
c
      dimension vcc(*),ncmph(10,18,nirr1),final(*)
c
c...  routine to handle symmetry/doc blocked calling of mpfab
c...  only diagonal blocks need be considered
c...  (cf mpsrcc for blocking conventions)
c
      kvcc = 1
      ic = icc
c
      do 100 ijrr=1,nirr1
        do 90 is=1,2+nirr*2
          ncsfi = ncmph(2,is,ijrr)
          nseti = ncmph(3,is,ijrr)
          do 80 iset=1,nseti
            call mpfab(ic,vcc(kvcc),ncsfi,ncmph(1,is,ijrr),final)
            kvcc = kvcc + ncsfi*ncmph(9,is,ijrr)
            ic = ic + ncmph(10,is,ijrr)
80        continue
90      continue
100   continue
c
      return
      end
      subroutine cmpaab(iccd,vcd,ncmphd,icc2,vc2,ncmph2,nirr,final)
c
      implicit real*8 (a-h,o-z)
c
      dimension vcd(*),ncmphd(10,18),vc2(*),ncmph2(10,18),final(*)
c
c...  routine to handle symmetry/doc blocked calling of mpfaab
c...  only diagonal blocks need be considered (and symmetry 1)
c...  (cf mpsrcc for blocking conventions)
c
      kvc2 = 1
      kvcd = 1
      ic2 = icc2
      icd = iccd
c
      do 90 is=1,2+nirr*2
         ncsfi = ncmphd(2,is)
         ncsfj = ncmph2(2,is)
         nseti = ncmph2(3,is)
         do 80 iset=1,nseti
            call mpfaab(icd,ic2,vcd(kvcd),ncsfi,ncmphd(1,is),
     *                 vc2(kvc2),ncsfj,ncmph2(1,is),final)
            kvcd = kvcd + ncsfi*ncmphd(9,is)
            kvc2 = kvc2 + ncsfj*ncmph2(9,is)
            icd = icd + ncmphd(10,is)
            ic2 = ic2 + ncmph2(10,is)
80       continue
90    continue
c
      return
      end
      subroutine mpfij(conf,mstart,mlast,nstart,nlast,icc,irr,
     *                 iadv,iadvt,vcc,ncmphc,vcr,ncmphr, fock,
     *                 result,nstatc,nstatr,cx,cy,temp1,final)
c
      implicit real*8 (a-h,o-z)
      integer conf
      dimension conf(5,*),iadv(*),iadvt(*)
      dimension vcc(nstatc,*),ncmphc(10),vcr(nstatr,*),ncmphr(10)
      dimension fock(*),result(nstatr,nstatc)
      dimension cx(*),cy(*),temp1(*),final(*)
      dimension xnumb(2)
      logical diag,contr,notri
c
c     compute the fock-matrix between vacuum states (in a sub-block)
c     for MR-MP theory
c     **     this subroutine computes the internal parts for     **
c     ** the (n-1/h/n-1) and (n-2/h/n-2) matrix elements as well **
c
c     icc,icr : numbering of excited states (starting at 1)
c     iadv/iadvt : addressing of internal space
c     vcc,vcr : contraction-matrices with excited states
c     ncmphc,ncmphr : excited state info
c     fock : internal-fock-matrix (full-triangle)
c     result : temp array for coupling coefficients
c     temp : scratch-array's (max(nstatc**2,nstatr**2)
c     final : final array for matrix-elements (max # excited states**2)
c     cx,cy : matrix-elements storage (maxpa**2)
c     nstate : # internal csf's
c
c...  //note// first all singlets then triplets
c...           second conf defines most rapidly varying index
c...  ** no need for seperate singlet-triplet in final matrix-elements
c
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
c
c
      integer nref, nnst, nmin1, nmin2, mspin, nelec
      integer nstrt1, nlast1, nstrt2, nlast2, nmin12, nminv2, nminv1
      integer nlst21, nlastx
      common /ccntl/ nref,nnst,nmin1,nmin2,mspin,nelec,
     +               nstrt1,nlast1,nstrt2,nlast2,nmin12,nminv2,nminv1,
     +               nlst21,nlastx
c
c     common/spew/model,ic,ir
      common/spew/rmodel(5)
      equivalence (rmodel(1),model),(rmodel(2),ic),(rmodel(3),ir)
      common/symcon/norbbb(7),nd,ii,jj,kk,ll,konte(14)
      common/three/n1l(1024),n2l(1024),minlhs(64),maxlhs(64),iocc(64),
     *n1r(1024),n2r(1024),minrhs(64),maxrhs(64),jocc(64),modelr
c
      integer nint, next, ninnex, nwm1, nw, nwp1, nintwo
      integer nint21, nintp1
      common /orb/ nint,next,ninnex,nwm1,nw,nwp1,nintwo,
     +             nint21,nintp1
c
c
      data xnumb/1.0d0,2.0d0/
c
      ind(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
c
c...   determine if it is a diagonal block
c...   and if any triplets are involved
c
      diag = mstart.eq.nstart
      notri = ncmphc(5).eq.0.or.ncmphr(5).eq.0
c
c...  contr indicates if there is any contribution to the fock-matrix
c...  this could be sorted out beforehand actually
c
      contr = .false.
c
      if (nlast.lt.nstart.or.mlast.lt.mstart) return
c
      call vclr(result,1,nstatr*nstatc)
c
c..    loop over states considered
c
      do 200 ic=mstart,mlast
         call expan1(conf(1,ic))
         call upack2(conf(3,ic),ncolt,ncols)
         if (notri) ncolt = 0
c
         do 100 ir=nstart,nlast
            call upack2(conf(3,ir),nrowt,nrows)
            call expan2(conf(1,ir))
c
c...  establish number of differences
c...  if more then 2 => no future in it
c...    (nd comess from expan2 in /symcon/)
c
            if (nd.gt.2) go to 100
            if (nd.eq.2) go to 60
c
c.... no differences => diagonal element
c
            contr = .true.
            if (ic.ne.ir) call caserr('impossible in mpfij')
            h = 0.0d0
            do 50 i=1,nint
               if (iocc(i).ne.0) h = h + fock(ind(i,i))*xnumb(iocc(i))
50          continue
c....  add 1-el contribution to diagonals
            do 51 i=1,ncols
51          result(iadv(ic)+i,iadv(ic)+i) = h
            do 52 i=1,ncolt
52          result(iadvt(ic)+i,iadvt(ic)+i) = h
            go to 100
c---- end of diagonal contributions --------------------------------
c
c.... single difference
c
60          h = fock(ind(ii,jj))
            if (h.eq.0.0d0) go to 100
            ndims=nrows*ncols
            ndimt=nrowt*ncolt
            contr = .true.
c.... now make h contributions  (into full square)
            call symb1
c...   singlets
            call sym1tr(cx,ifl,0,0)
            if (ifl.ne.0) then
               call dscal(ndims,h,cx,1)
               kk = 0
               do 61 kc=1,ncols
                do 61 kr=1,nrows
                  kk = kk + 1
                  result(iadv(ir)+kr,iadv(ic)+kc) = cx(kk)
61             continue
            end if
c...   triplets
            call sym1tr(cy,ifl,1,1)
            if (ifl.ne.0.and..not.notri) then
               call dscal(ndimt,h,cy,1)
               kk = 0
               do 62 kc=1,ncolt
                do 62 kr=1,nrowt
                  kk = kk + 1
                  result(iadvt(ir)+kr,iadvt(ic)+kc) = cy(kk)
62             continue
            end if
c
c---- end of delta = 1 -------------------------------------------------
c
100      continue
c
200   continue
c
c...   if fock-matrix empty return
c
       if (.not.contr) return
c
c...    we now have a full internal fock-matrix in result
c
       if (iprimp.ge.1000) then
          write(6,600) model,mstart,mlast,nstart,nlast
600    format(/' internal FOCK matrix (model',i2,') before contraction'
     *        /'  ic from',i6,' to',i6,' / ir from',i6,' to',i6,'    ')
          call  prvc4 (result,nstatr,nstatc,idum,dum,0,10)
       end if
c
c...  apply contractions // no interactions between sets
c...  each set will count as one (1) seperate config in mpgen
c...  if this gives core-problems we can do otherwise
c...  3 sets are always assumed dimensions in ncmph(6..8)
c
      ic = icc
      ir = irr
c
      ncols = ncmphc(4) + ncmphc(5)
      nrows = ncmphr(4) + ncmphr(5)
      ivcc = 1
      ivcr = 1
c
      do 500 iex=6,8
c
         nexc = ncmphc(iex)
         nexr = ncmphr(iex)
c
         if (nexc.le.0.or.nexr.le.0) go to 490
c
         call wrtsb(rmodel(1),3)
c
c...     multiply  vcr+ result vcc
c
         call mxmd(vcr(1,ivcr),nstatr,1,result,1,nstatr,
     *             temp1,1,nexr,nexr,nrows,nstatc)
c
         call mxmd(temp1,1,nexr,vcc(1,ivcc),1,nstatc,
     *             final,1,nexr,nexr,ncols,nexc)
c
c..       substract e0mp for diagonal
c
         if (diag) then
            ii = 1
            do 410 i=1,nexr
               final(ii) = final(ii) - e0mp
               ii = ii + nexr + 1
410         continue
         end if
c
c...    print coefficients (**debug**)
c
         if (iprimp.ge.1000) then
            write(6,603) ic,ir
603         format(/' internal FOCK matrix aft contraction - E0mp'
     *             ,' ic,ir .. ',2i5)
            call  prvc4 (final,nexr,nexc,idum,dum,0,10)
         end if
c
c...   write  matrix-elements to symbolic file
c
         call wrtsb(final(1),nexr*nexc)
c
490      continue
c
c...   adjust ic,ir // ivcc,ivcr
c
         if (nexc.gt.0) then
            ic = ic + 1
            ivcc = ivcc + nexc
         end if
         if (nexr.gt.0) then
            ir = ir + 1
            ivcr = ivcr + nexr
         end if
c
500   continue
c
      return
      end
      subroutine mpfia(conf,mstart,mlast,nstart,nlast,icc,irr,
     *                 iadv,iadvt,vcc,ncmphc,vcr,ncmphr,
     *                 result,nstatc,nstatr,cj,temp,final,
     *                 npp,ipair,mxpair)
c
      implicit real*8 (a-h,o-z)
c
      integer  conf,conner
c
      dimension conf(5,*),iadv(*),iadvt(*)
      dimension vcc(nstatc,*),ncmphc(10),vcr(nstatr,*),ncmphr(10)
      dimension result(nstatr,nstatc),cj(*),temp(*),final(*)
      dimension ipair(2,mxpair,*),npp(*)
c
c...  to evaluate spin coupling coeffs. for ia(fock) integrals
c...  for MR-MP theory.  / analagous to mpfij
c...  (temp) = max(nstatr**2,nstatc**2)
c...  (final) = 2*(maximum # excited states)**2
c...
c...  ipair(2,mxpair,nint) is used to gather the config's involved
c...  this may be made still more efficient (cf also mpfini)
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
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
c
      logical kltj
c
c
      integer nref, nnst, nmin1, nmin2, mspin, nelec
      integer nstrt1, nlast1, nstrt2, nlast2, nmin12, nminv2, nminv1
      integer nlst21, nlastx
      common /ccntl/ nref,nnst,nmin1,nmin2,mspin,nelec,
     +               nstrt1,nlast1,nstrt2,nlast2,nmin12,nminv2,nminv1,
     +               nlst21,nlastx
c
      common/dopey/ndimsd,ncols,ncold,isingl(62)
c     common/spew/model,ic,ir,ninti
      common/spew/rmodel(5)
      equivalence (rmodel(1),model),(rmodel(2),ic),(rmodel(3),ir)
      equivalence (rmodel(4),ninti)
      common/symcon/norbbb(7),ndiff,ii,jj,kk,ll,kount,kounte
c
      integer nint, next, ninnex, nwm1, nw, nwp1, nintwo
      integer nint21, nintp1
      common /orb/ nint,next,ninnex,nwm1,nw,nwp1,nintwo,
     +             nint21,nintp1
c
      common/expanc/conner(10),iocc(64)
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
c
      integer nirr,mult,ifrep,nirr8,jfrep
      common /table/ nirr,mult(8,8),ifrep,nirr8,jfrep
c
c
      integer iro, iky, niky, m1e, m2e, m1ia, m2coul
      integer m2exch, m2fock, m2psup, m2qsup
      integer mbase, nbase, npairi, norbi, mbasi, mnbasi
      integer iap, mbfock, mbexch, mbcoul, mbpsup, mbqsup
      integer norbas, m2exab, m2em1e
      common /helpr/ iro(nd200),iky(nd200),niky(8),m1e,m2e,m1ia,m2coul,
     +   m2exch,m2fock,m2psup,m2qsup,
     +   mbase(288),nbase(8),npairi(36),norbi(8),mbasi(8),mnbasi(8),
     +   iap(62),mbfock(62),mbexch(36),mbcoul(36),mbpsup(36),mbqsup(36),
     +   norbas(nd200),m2exab,m2em1e
c
c
      common/stoctl/khbill,khbz,khbc,khbb
      dimension jap(62)
c
      if(mlast.lt.mstart.or.nlast.lt.nstart) return
c
c....  sort out needed pairs
c
      mx = 0
      call izero(nint,npp,1)
      do 20 ic = mstart,mlast
         do 10 ir=nstart,nlast
            ndiff =  ndif(conf(1,ic),conf(1,ir),ii)
            if (ndiff.eq.2) then
               npp(ii) = npp(ii) + 1
               mx = max(npp(ii),mx)
               ipair(1,npp(ii),ii) = ic
               ipair(2,npp(ii),ii) = ir
            end if
10       continue
20    continue
c
      if (mx.eq.0) return
      if (mx.gt.mxpair) call caserr('pair overflow in MPFIA')
c
c...   get all the common symmetry info
c
         isym=conf(4,mstart)
         jsym= conf(4,nstart)
         ksym=mult(jsym,isym)
         nek=norbe(ksym)
         if (nek.eq.0) return
c
c...     for normal (non-diagonal n-2 skip if only 1 orbital available
c
         if (model.eq.36.and.isym.eq.1.and.nek.le.1) then
            return
         end if
c
c...  set up internal orbital addressing
c
      do i=1,nint
       jap(i) = iap(i) - iap(1)
      enddo
c
c...    we need to treat each (ia) fock operator seperately
c...    so loop over internal orbitals
c
      do 99999 iiiint=1,nint
c
      if (npp(iiiint).eq.0) go to 99999
c
c... loop over internally least occupied states   (doublets/n-2)
c
      call vclr(result,1,nstatr*nstatc)
c
      do 200 ip = 1,npp(iiiint)
         ic = ipair(1,ip,iiiint)
c
         call expan1(conf(1,ic))
         call upack2(conf(3,ic),ncolt,ncolss)
         ncols=ncolss+ncolt
c
c... loop over internally more occupied state   (vacuum/doublet))
c
            ir = ipair(2,ip,iiiint)
c
            ncold= conf(3,ir)
            call expan2(conf(1,ir))
c... if no. of differences ne 2 not interested (and now impossible)
            if (ndiff.ne.2) then
               print *,' error ',iiiint,ic,ir,ip,nint,ndiff
               print *,' npp ',(npp(i),i=1,nint)
               call caserr('ndiff ne 2 ')
            end if
c
c...    is it the right fock-operator ???
c
            if (ii.ne.iiiint) call caserr(' ii ne iiiint ')
c
            ndimsd=ncold*ncols
            mubas=ncolss*ncold
            kltj=ksym.lt.jsym.or.ncolt.eq.0
            nubas=ndimsd-mubas
c
c... single orthogonality-----------------------------------------------
c... evaluate fock contributions
c
            call symb1
            call symbt(cj,0,kount)
            if (.not.kltj) call dscal(nubas,-1.0d0,cj(mubas+1),1)
c
c..         coefficients in cj (first singlet then triplet)
c
c...        write them to result
c
            kk = 0
            do kc=1,ncolss
             do kr=1,ncold
                kk = kk + 1
                result(iadv(ir)+kr,iadv(ic)+kc) = cj(kk)
             enddo
            enddo
            do kc=1,ncolt
             do kr=1,ncold
                kk = kk + 1
                result(iadv(ir)+kr,iadvt(ic)+kc) = cj(kk)
             enddo
            enddo
c
100      continue
c
200   continue
c
c...    we now have all fock-couplings in result
c
       if (iprimp.ge.1000) then
          write(6,600) model,iiiint,mstart,mlast,nstart,nlast
600       format(/' Fia matrix (model',i2,' orbital ',i2,
     *           ') before contraction',
     *           /'  ic from',i6,' to',i6,' / ir from',i6,' to',i6)
          call  prvc4 (result,nstatr,nstatc,idum,dum,0,10)
       end if
c
c...  apply contractions // no interactions between sets
c...  each set will count as one (1) seperate config in mpgen
c...  if this gives core-problems we can do otherwise
c...  3 sets are always assumed dimensions in ncmph(6..8)
c...  internally least occupied (n-2/doublet) is ic
c...  since ic is the second index : produce first singlets and write
c...                                 triplets next
c
      ic = icc
      ir = irr
c
      ncols = ncmphc(4)
      nrows = ncmphr(4)
      ncolt = ncmphc(5)
      ivcc = 1
      ivcr = 1
c
      do 500 iex=6,8
c
         nexc = ncmphc(iex)
         nexr = ncmphr(iex)
c
         if (nexc.le.0.or.nexr.le.0) go to 490
c
         kfin = 1
c
c...     multiply  vcr+ result vcc
c
         call mxmd(vcr(1,ivcr),nstatr,1,result,1,nstatr,
     *             temp,1,nexr,nexr,nrows,nstatc)
c
c...   singlet first (if available)
c
         if (ncols.gt.0) then
c
            call mxmd(temp,1,nexr,vcc(1,ivcc),1,nstatc,
     *                final(kfin),1,nexr,nexr,ncols,nexc)
c
            kfin = kfin + nexr*nexc
c
         end if
c
c...     now same for triplets
c
         if (ncolt.gt.0) then
c
c...     multiply   result vcc
c
            call mxmd(temp(ncols*nexr+1),1,nexr,vcc(ncols+1,ivcc),
     *                1,nstatc,final(kfin),1,nexr,nexr,ncolt,nexc)
            kfin = kfin + nexr*nexc
         end if
c
c...     print coefficients (**debug**)
c
         if (iprimp.ge.1000) then
            nn = (kfin-1)/nexr
            write(6,603) model,iiiint,ic,ir
603         format(/' Fia matrix (model',i2,' orbital ',i2,
     *       ') ic/ir ',2i5)
            call  prvc4 (final,nexr,nn,idum,dum,0,10)
         end if
c
c...   write  matrix-elements to symbolic file
c
c
         ninti = jap(iiiint)
         call wrtsb(rmodel(1),4)
         call wrtsb(final(1),kfin-1)
c
490      continue
c...   adjust ic/ir
         if (nexc.gt.0) then
            ic = ic + 1
            ivcc = ivcc + nexc
         end if
         if (nexr.gt.0) then
            ir = ir + 1
            ivcr = ivcr + nexr
         end if
c
c...   end of loop over contractors
c
500   continue
c
c...   end of loop over fock-operators
c
99999 continue
c
      return
      end
      subroutine mpfab(icc,vcc,nstatc,ncmphc,final)
c
      implicit real*8 (a-h,o-z)
c
      dimension vcc(nstatc,*),ncmphc(10),final(*)
c
c...  to evaluate symbolic for ab(fock) integrals
c...  for MR-MP theory.  cf. mpfij  (only diagonal interactions here)
c...  the coupling coefficients for the individual csf's form a unit
c...  matrix : get couplings for excited states by v+v
c...  no singlet-triplet interactions can arise !!
cjvl  *** in definite program doublets should not be done here '
c...
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
c
c
      integer nref, nnst, nmin1, nmin2, mspin, nelec
      integer nstrt1, nlast1, nstrt2, nlast2, nmin12, nminv2, nminv1
      integer nlst21, nlastx
      common /ccntl/ nref,nnst,nmin1,nmin2,mspin,nelec,
     +               nstrt1,nlast1,nstrt2,nlast2,nmin12,nminv2,nminv1,
     +               nlst21,nlastx
c
c     common/spew/model,ic,ir
      common/spew/rmodel(5)
      equivalence (rmodel(1),model),(rmodel(2),ic),(rmodel(3),ir)
      common/stoctl/khbill,khbz,khbc,khbb
c
c...  apply contractions // no interactions between sets
c...  each set will count as one (1) seperate config in mpgen
c...  if this gives core-problems we can do otherwise
c...  3 sets are always assumed dimensions in ncmph(6..8)
c
      ic = icc
c
      ncols = ncmphc(4)
      ncolt = ncmphc(5)
      nstate = ncols+ncolt
      ivcc = 1
c
      do 500 iex=6,8
c
         nexc = ncmphc(iex)
c
         if (nexc.le.0) go to 500
c
         call wrtsb(rmodel(1),2)
c
         if (model.eq.37.and.iprimp.lt.1000) go to 499
         if (ifdiaf.eq.1) go to 499
c
c...   singlet first (if available)
c
         if (ncols.gt.0) then
c
c...     multiply  vcc+ vcc
c
            call mxmd(vcc(1,ivcc),nstate,1,vcc(1,ivcc),1,nstate,
     *                final,1,nexc,nexc,ncols,nexc)
c
c...    print coefficients (**debug**)
c
            if (iprimp.ge.1000) then
               write(6,600) model,ic
600            format(/' Fab matrix (model',i2,' ic= ',i2,')')
               call  prvc4 (final,nexc,nexc,idum,dum,0,10)
            end if
c
c...   write  matrix-elements to symbolic file
c
            if (model.ne.37) then
               call wrtsb(final(1),nexc*nexc)
            end if
         end if
c
c...     now same for triplets
c
         if (ncolt.gt.0) then
c
c...     multiply  vcr+ result(triplet-part) vcc
c
            call mxmd(vcc(ncols+1,ivcc),nstate,1,
     *                vcc(ncols+1,ivcc),1,nstate,
     *                final,1,nexc,nexc,ncolt,nexc)
c
c
c...    print coefficients (**debug**)
c
            if (iprimp.ge.1000) then
               write(6,601) model,iex
601            format(/' triplet Fab matrix (model',i2,' # ',i2,')')
               call  prvc4 (final,nexc,nexc,idum,dum,0,10)
            end if
c
c...   write  matrix-elements to symbolic file
c
            call wrtsb(final(1),nexc*nexc)
         end if
c...   adjust ic,ivcc
499      ic = ic + 1
         ivcc = ivcc + nexc
c
500   continue
c
      return
      end
      subroutine mpfaab(icc,irr,vcc,nstatc,ncmphd,vcr,nstatr,ncmph2,
     *                  final)
c
      implicit real*8 (a-h,o-z)
c
      dimension vcc(nstatc,*),ncmphd(10),vcr(nstatr,*),ncmph2(10)
      dimension final(*)
c
c...  to evaluate symbolic for ab(fock) integrals
c...  for MR-MP theory. cf. mpfij for (aa/ab) interactions
c...  on the aa side we only have singlets
c...  so only singlet interactions !!
c...  (note that (aa/aa)-type interactions are done by mpfab)
c
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
c
c     common/spew/model,ic,ir
      common/spew/rmodel(5)
      equivalence (rmodel(1),model),(rmodel(2),ic),(rmodel(3),ir)
c
c...     ab (ir/left)  <=>  aa (ic/right)
c
      ic = icc
      ir = irr
      ivcc = 1
      ivcr = 1
c
      if (nstatc.ne.ncmphd(4).or.ncmphd(4).ne.ncmph2(4)) then
         call prmpcon
         print*,'** nstatc = ',nstatc
         print*,'** ncmphd = ',ncmphd
         print*,'** nstatr = ',nstatr
         print*,'** ncmph2 = ',ncmph2
         call caserr('program fault in mpfaab')
      endif
      nsingl = ncmphd(4)
c
      do 500 iex=6,8
c
         nexc = ncmphd(iex)
         nexr = ncmph2(iex)
c
         if (nexc.le.0.or.nexr.le.0) go to 490
c
c...     multiply
c
         call mxmd(vcr(1,ivcr),nstatr,1,vcc(1,ivcc),1,nstatc,
     *             final,1,nexr,nexr,nsingl,nexc)
c
c...    print coefficients (**debug**)
c
         if (iprimp.ge.1000) then
            write(6,601) model
601   format(' faab fock-matrix model ',i2)
            call  prvc4 (final,nexr,nexc,idum,dum,0,10)
         end if
c
c...   write  matrix-elements to symbolic file
c...   ** matrix is not symmetrical / nor is it used that way **
c
         call wrtsb(rmodel(1),3)
         call wrtsb(final(1),nexr*nexc)
c
490      continue
c
         if (nexr.gt.0) then
            ir = ir +1
            ivcr = ivcr + nexr
         end if
         if (nexc.gt.0) then
            ic = ic +1
            ivcc = ivcc + nexc
         end if
c
500   continue
c
      return
      end
      subroutine prvc4 (c,n,m,ioc,eps,ich,nc)
      implicit real*8 (a-h,o-z)
c
c
c   routine for printing of column-vectors and eigenvalues.
c
c            c  : matrix of coefficients
c            n  : dimension of c(n,m) and eps(n)
c            m  : x columns to be printed
c            eps: vector of eigenvalues
c            ioc: vector of (integer) (occupation)-numbers
c            n  : dimension of c(n,n) and eps(n)/ioc(n)
c            ich: 1 => ioc / 0 => column-numbers / -1 => eig
c            nc : maximum number of columns on one row
c
c

      dimension c(n,m), eps(m) , ioc(m)
      character*5 sub(2)
      data sub /'-----','-----'/
      ncc = nc
      if ( ncc.gt.10 ) ncc=10
      nbl = (m-1)/ncc
      nlast = m-ncc*nbl
      nbl = nbl+1
      nc1 = 0
      do 60 ib=1,nbl
      if ( ib .eq. nbl ) ncc=nlast
      nc0 = nc1+1
      nc1 = nc1+ncc
      if (ich.lt.0) write(6,10) ( eps(ic), ic=nc0,nc1 )
      if (ich.gt.0) write(6,11) ( ioc(ic), ic=nc0,nc1 )
      if (ich.eq.0) write(6,11) ( ic, ic=nc0,nc1 )
   10 format(/8x,10f12.6)
   11 format(/8x,10(4x,i4,4x))
      write(6,20) ( sub, i=1,ncc )
   20 format(8x,10(2x,2a5))
c      write(6,30)
c   30 format(1h )
      do 50 ia=1,n
      write(6,40) ia, ( c(ia,ic), ic=nc0,nc1 )
   40 format(3x,i3,2x,10f12.6)
   50 continue
   60 continue
c      write(6,30)
c
      return
      end
      subroutine whermp(iw,nw)
c...   copied from whersb but more general / used in mp2
      implicit real*8 (a-h,o-z)
      common/disksa/bof(511),no,iblo
      nw=no
      if (no)888,999,888
999   iw = iblo
      goto 777
888   iw = iblo+1
777   return
      end
c
c===== ORTHOGONALISATION SCHEMES FOR PERTURBATION THEORY ===============
c
c     It was found that multi reference perturbation theory is 
c     extremely sensitive to the accuracy with which the basis of 
c     excited states is orthogonalised.
c
c     The researches have led us past the following orthonalisation
c     methods:
c
c        modified Gramm-Schmidt
c        modified Gramm-Schmidt in quadruple precision (real*8*16)
c        pivoted repeated modified Gramm-Schmidt
c        Lowdin
c        Householder
c
c     The only method that solved our problems was the Householder
c     QR factorisation method. All other methods were to inaccurate. 
c     This was shown using CASMP calculations of the size consistency
c     error in the Berylium dimer using the 6-311G* basisset. The
c     CAS was constructed from the 2S and 2P orbitals and 2 electrons
c     per Berylium atom. In the dimer the interatomic distance was
c     set to 1000 bohr.
c
c     Next, the implementation of the various orthogonalisation schemes
c     is given.
c
c                                                           Huub van Dam
c                                                         March 28, 1997
c
c-----------------------------------------------------------------------
c
      subroutine mpsymm(nref,nsing,ndoub,v,lv,nv,kept,cqr,ipr,
     *                  a,q,t,iky,qv)
c                                   
      implicit real*8 (a-h,o-z)
c
c     Symmetric-orthonormalisation for MR-MP
c     ======================================
c
c     Purpose:
c     --------
c
c     For efficient implementation of multi reference Moller-Plesset
c     perturbation theory (MRMP) orthonormalisation of the excited 
c     state space is required. 
c     If the space is strongly linear dependent then modified Gramm-
c     Schmidt is not accurate enough and Lowdin-orthonormalisation
c     may perform better.
c
c     In the MRMP approach formulated by Wolinski et al. the 
c     reference function, the singly substituted, and the doubly
c     substituted spaces should be mutually orthogonal. This means
c     that we can not apply Lowdin-orthogonalisation to the total
c     space. Lowdin has to be applied to the subspaces which are
c     subsequently orthogonalised to each other through projections.
c
c     nref  :
c        The number of reference functions in <v>.
c
c     nsing :
c        The number of singly substituted functions in <v>.
c
c     ndoub :
c        The number of double substituted functions in <v>.
c
c     v     :
c        Space spanned by <nv> column-vectors of length <lv>.
c        The basis of this space is returned in this variable.
c
c     kept  :
c        Vector of length <nv> indicating which vectors are
c        linearly independent and should be kept.
c        If kept(i) = 1 the i-th vector will be kept.
c        If kept(i) = 0 the i-th vector will be dropped.
c
c     cqr   :
c        Criterion for orthogonality, i.e. two vectors are orthogonal
c        if their inproduct is less than <cqr>.
c
c     ipr   :
c        Print flag.
c
c
c     Workspace
c     ---------
c
c     define  ns = max(nref, nsing, ndoub)  is the maximum dimension
c     of a subspace, then
c
c     a     :
c        Triangular matrix of dimension ns*(ns+1)/2 to store the
c        overlap matrix.
c
c     q     :
c        Rectangular matrix of dimension ns*(ns+1) to store the
c        coefficient matrix.
c
c     t     :
c        Rectangular matrix of dimension lv*ns*2 to store the
c        vectors under consideration and the transformed vectors
c        before the back substitution.
c
c     iky   :
c        Vector of dimension ns*2 to store offsets for the diag
c        subroutine.
c
c     qv   :
c        Matrix in QUADRUPLE PRECISION of dimension lv*nv to store
c        the column vectors. This matrix is not used at the same
c        time as the workspace variables a, q, t, and iky. Therefore
c        the memory areas of qv may overlap with that of the other
c        variables.
c
c
c     Theory:
c     -------
c
c     Define the overlap matrix  S = <V|V>  then there is a unitary
c     transformation  U  such that  L = <U|S|U>  where  L  is a 
c     diagonal matrix. For the non-zero eigenvalues of  L  we have
c     I = <U L^(-1/2)|S|U L^(-1/2)> = <V U L^(-1/2)|V U L^(-1/2)>
c     i.e. V U L^(-1/2) is an orthonormal space.
c
c     Algorithm:
c     ----------
c
c     Lowdin all reference functions
c
c     Lowdin all singly substituted functions.
c     Gramm-Schmidt singles on reference, i.e. project reference out.
c     Lowdin all singly substituted functions.
c
c     Lowdin all doubly substituted functions.
c     Gramm-Schmidt doubles on reference and singles, i.e. project
c     reference and singles out.
c     Lowdin all doubly substituted functions.
c
      call mplowd(1,nref,v,lv,nv,kept,cqr,ipr,a,q,t,iky)
c
c     call mplowd(nref+1,nref+nsing,v,lv,nv,kept,cqr,ipr,a,q,t,iky)
c     call qmpschm(1,nref,nref+1,nref+nsing,v,lv,nv,kept,cqr,ipr,qv)
      call mpschm(1,nref,nref+1,nref+nsing,v,lv,nv,kept,cqr,ipr)
      call mplowd(nref+1,nref+nsing,v,lv,nv,kept,cqr,ipr,a,q,t,iky)
c
c     call mplowd(nref+nsing+1,nv,v,lv,nv,kept,cqr,ipr,a,q,t,iky)
c     call qmpschm(1,nref+nsing,nref+nsing+1,nv,v,lv,nv,kept,cqr,ipr,qv)
      call mpschm(1,nref+nsing,nref+nsing+1,nv,v,lv,nv,kept,cqr,ipr)
      call mplowd(nref+nsing+1,nv,v,lv,nv,kept,cqr,ipr,a,q,t,iky)
c
      end 
c
c-----------------------------------------------------------------------
c
      subroutine mplowd(ia,ib,v,lv,nv,kept,clowd,ipr,a,q,t,iky)
c                                   
      implicit real*8 (a-h,o-z)
c
c     Lowdin-orthonormalisation for MR-MP
c     ===================================
c
c     Purpose
c     -------
c
c     Orthonormalises the column vectors in the subspace 
c     v[ia..ib] using Lowdin orthogonalisation.
c
c     ia    , ib    :
c        Lower and Upper limits of the current subspace.
c
c     v     : [1..lv,1..nv]
c        Space spanned by <nv> column-vectors of length <lv>.
c        The basis of this space is returned in this variable.
c
c     kept  : [1..nv]
c        Vector of length <nv> indicating which vectors are
c        linearly independent and should be kept.
c        If kept(i) = 1 the i-th vector will be kept.
c        If kept(i) = 0 the i-th vector will be dropped.
c
c     clowd :
c        Orthogonalisation criterion. In practice it is used as
c        criterion for the diagonalisation routine and as lower
c        limit for non-zero eigenvalues of the overlap matrix.
c
c     ipr   :
c        Print flag.
c        If ipr > 1 the accuracy of the orthonormalisation will be
c        checked and reported.
c        If ipr > 2 the eigenvalue of every vector will be reported
c        as well as the predicate kept or dropped.
c
c  
c     Workspace
c     ---------
c
c     define ns = ib - ia + 1 is the dimension of the current 
c     subspace, then
c
c     a     : 
c        Triangular matrix of dimension ns*(ns+1)/2 to store the
c        overlap matrix.
c
c     q     : 
c        Rectangular matrix of dimension ns*(ns+1) to store the
c        coefficient matrix.
c
c     t     : 
c        Rectangular matrix of dimension lv*ns*2 to store the
c        vectors under consideration and the transformed vectors
c        before the back substitution.
c
c     iky   : 
c        Vector of dimension ns*2 to store offsets for the diag
c        subroutine.
c
c
c     Theory
c     ------
c
c     Define the overlap matrix  S = <V|V>  then there is a unitary
c     transformation  U  such that  L = <U|S|U>  where  L  is a 
c     diagonal matrix. For the non-zero eigenvalues of  L  we have
c     I = <U L**(-1/2)|S|U L**(-1/2)> = <V U L**(-1/2)|V U L**(-1/2)>
c     i.e. V U L**(-1/2) is an orthonormal space.
c
c
      real*8 dabs
      real*8 ddot
      external ddot
      parameter (dzero =1.0d-15)
      dimension v(lv,nv),kept(nv)
c
c...  itri: indexing for lower triangular matrix including diagonal
c...  irec: indexing for rectangular matrix.
c...        columns length  n
c...        column number   j
c...        row number      i
c
      dimension a(*),q(*),t(*)
      dimension iky(*)
c
      itri(i,j)   = max(i,j)*(max(i,j)-1)/2+min(i,j)
      irec(n,i,j) = (j-1)*n+i
c
      depend = clowd
      xnorm  = 0.0d0
c
c...  ipe: pointer to the array of eigenvalues on q.
c...  ipt2: pointer to temp array for transformed vectors on t.
c...  ipil: pointer to index array for q on iky.
      ns   = ib - ia + 1
      ipe  = ns ** 2
      ipt2 = lv * ns
      ipil = ns
c
c...  Copy input vectors to temporary variables
c
      j = 0
      do 17 i = ia, ib
        if (kept(i).eq.1) then
          j = j + 1
          call dcopy(lv,v(1,i),1,t(irec(lv,1,j)),1)
        endif
17    continue
c
c...  Count all vectors still in the race.
c
      nkept = j
c
c...  Calculate overlap matrix
c
      do 300 i = 1, nkept
        it = irec(lv,1,i)
        do 310 j = 1, i
          jt = irec(lv,1,j)
          a(itri(i,j)) = ddot(lv,t(it),1,t(jt),1)
310     continue
        iky(i)      = itri(i,0)
        iky(ipil+i) = irec(nkept,0,i)
300   continue
c
c...  Diagonalise overlap matrix
c
      crit = 1.0d0 / dsqrt(1.0d0*nkept) * clowd
      call jacobi(a,iky,nkept,q,iky(ipil+1),nkept,q(ipe+1),2,1,crit)
c
c...  jacobi call parameter one before last determines eigenvalue order
c...  0 nop / 1 unsorted / 2 increasing / 3 decreasing / 4 lock
c
c...  diag call parameter two before last determines initialisation
c...    1.   use vectors in q as such
c...    2.   initialise q to the unit matrix
c...  diag call parameter one before last determines eigenvalue order
c...  0 nop / 1 unsorted / 2 increasing / 3 decreasing / 4 lock
c...    1.   unsorted
c...    2.   increasing
c...    3.   decreasing
c...    4.   lock
c
c...  Select proper vectors. I.e.
c...  If norm        L[i] < depend  then throw away. 
c...  If eigenvalue  L[i] < zero    then throw away.
c
      do 400 i = 1,nkept
        da = q(ipe+i)
        dd = dsqrt(dabs(da))
        if (da.lt.dzero.or.da.lt.depend) then
          q(ipe+i) = 0.0d0
          call vclr(q(irec(nkept,1,i)),1,nkept)
          if (ipr.gt.2) write(*,900)i,da
          if (da.gt.xnorm) xnorm = da
        else
          if (ipr.gt.2) write(*,910)i,da
          dd = 1.0d0/dd
          call dscal(nkept,dd,q(irec(nkept,1,i)),1)
        endif
400   continue
900   format('** vector ',i6,' with eigenvalue ',e17.5,' is dropped')
910   format('** vector ',i6,' with eigenvalue ',e17.5,' is kept   ')
c
c...  Compute new orthogonal vectors
c
      call mxmd(t,1,lv,q,1,nkept,t(ipt2+1),1,lv,lv,nkept,nkept)
c
c...  Write back the selected vectors
c
      do 500 itmp = 1, nkept
        i = nkept-itmp+1
        ikept = i
        do k = ia, ib
          if (kept(k).eq.1) ikept = ikept - 1
          if (ikept.eq.0) goto 530
        enddo
530     if (q(ipe+i).eq.0.0d0) then
          kept(k) = 0
          call dcopy(lv,t(ipt2+irec(lv,1,i)),1,v(1,k),1)
        else 
          kept(k) = 1
          call dcopy(lv,t(ipt2+irec(lv,1,i)),1,v(1,k),1)
        endif
500   continue
c
c...  Just to be save: normalise again
c
      do 600 i = ia, ib
        if (kept(i).eq.1) then
          dd = dsqrt(ddot(lv,v(1,i),1,v(1,i),1))
          dd = 1.0d0/dd
          call dscal(lv,dd,v(1,i),1)
        endif
600   continue
c 
c...  if printing on then find and print maximum overlap
c
      imax = 0
      jmax = 0
      smax = 0.0d0
      do 170 i = ia, ib
        if (kept(i).eq.1) then
          do 180 j = ia, i-1
            if(kept(j).eq.1) then
              dd = dabs(ddot(lv,v(1,i),1,v(1,j),1))
              if(dd.gt.smax) then
                smax = dd
                imax = i
                jmax = j
              endif
            endif
180       continue
        endif
170   continue
      if (ipr.gt.1.or.xnorm.lt.1.0d-10) then
         write(*,*)'*** LOWDIN orthonormalisation'
         write(*,940)xnorm
         write(*,920)smax
940      format('*** Maximum eigenvalue dropped = ',e17.5)
920      format('*** Maximum overlap            = ',e17.5)
      endif
      if (smax.gt.clowd) then
        write(*,930)smax,clowd
930     format(' ** Maximum overlap ',e12.5,' exceeds threshold ',e12.5,
     *          ' in Lowdin !!!')
        call caserr("Maximum overlap exceeds threshold !!!")
      endif
c
      end 
c
c-----------------------------------------------------------------------
c
      subroutine mpschm(ia,ib,ja,jb,v,lv,nv,kept,cschm,ipr)
c
      implicit real*8 (a-h,o-z)
c
c     Repeated Modified Gramm-Schmidt-orthonormalisation for MR-MP
c     ============================================================
c
c     Purpose
c     -------
c
c     Orthonormalises the column vectors  v[ja..jb]  to the
c     columns vectors  v[ia..ib]  that are markted  kept(i) = 1
c     using the modified Gramm-Schmidt method.
c     It is essential that  ia =< ja  and  ib =< jb and that the
c     vectors  v[ia..ja]  are normalized!
c
c     ia    , ib     :
c        The lower and upper index of the vector set  i  to which the
c        other vectors are to be orthogonalised.
c
c     ja    , jb     :
c        The lower and upper index of the vector set  j  which is to be
c        orthogonalised to set  i.
c
c     v     :
c        Space spanned by <nv> column-vectors of length <lv>.
c        The basis of this space is returned in this variable.
c
c     kept  :
c        Vector of length <nv> indicating which vectors are
c        linearly independent and should be kept.
c        If kept(i) = 1 the i-th vector will be kept.
c        If kept(i) = 0 the i-th vector will be dropped.
c
c     cschm :
c        Criterion for orthogonality, i.e. two vectors are orthogonal
c        if their inproduct is less than <cschm>.
c
c     ipr   :
c        Print flag.
c
c     Theory
c     ------
c
c     Assume a vector  A  and  B  where the inproduct  <A|A> = 1  then
c     then the vector  B' = B - <A|B> A  is orthogonal to  A.
c     However, beware of numerical instability due to these projections
c     and subsequent normalisation.
c
c     It seems reasonable to assume that a vector is linear dependent
c     on a set of vectors if the selfoverlap upon orthogonalisation 
c     is reduced by a factor in the order of the machine precision.
c     So this may be used as a criterion to throw a vector away.
c
c
      real*8 dabs
      real*8 ddot
      external ddot
c
      dimension v(lv,nv),kept(nv)
c
      parameter(zero = 1.0d-15)
c
      depend = cschm
      ortho  = cschm
      smax   = 0.0d0
c
      do 100 j = ja, jb
c
         if (kept(j).eq.0) goto 100
c
c...     calculate and store current self overlap
c
         rho1 = ddot(lv,v(1,j),1,v(1,j),1)
c
         ncyc = 0
c
c...     Orthogonalise vector  j  on the subspace  v[ia..min(ib,j-1)]
c
110      continue
         do 50 i = ia, min(ib,j-1)
            if (kept(i).eq.1) then
               dd = -ddot(lv,v(1,j),1,v(1,i),1)
               call daxpy(lv,dd,v(1,i),1,v(1,j),1)
            endif
50       continue
c
c...     If the vector is linear dependent then throw away
c...     else normalise and update  rho1  and check orthogonality.
c...     If needed recycle on the current vector.
c
         ncyc = ncyc + 1
         rho2 = ddot(lv,v(1,j),1,v(1,j),1)
         if (ncyc.gt.7) then
            write(*,*)'*** problem with vector ',j,rho2/rho1
            call caserr('*** orthogonalisation problem in MPSCHM')
         else if (ncyc.gt.3) then
            write(*,*)'*** problem with vector ',j,rho2/rho1,ncyc
         else if (rho2.le.rho1*depend) then
            kept(j) = 0
         else
            dd   = 1.0d0 / dsqrt(rho2)
            rho1 = rho1  / rho2
            call dscal(lv,dd,v(1,j),1)
c
c...        Test orthogonality
c
            do 120 i = ia, min(ib,j-1)
               if(kept(i).eq.0) go to 120
               dd = dabs(ddot(lv,v(1,i),1,v(1,j),1))
               if(dd.gt.ortho) then
                  if(ipr.gt.2) write(6,620)j,i,dd,ncyc,ortho,j
620               format('Overlap of vector ',i5,' with vector ',i5,
     *                   ' is ',e12.5,' in cycle ',i2,
     *                   ' this exceeds ',e12.5,
     *                   ' so recycle on vector ',i5)
                  goto 110
               endif
120         continue
c
            kept(j) = 1
         end if
c
100   continue
c
c...  test orthogonality
c
      do 200 j = ja, jb
         if (kept(j).eq.1) then
            do 210 i = ia, min(ib,j-1)
               if (kept(i).eq.1) then
                  dd = dabs(ddot(lv,v(1,i),1,v(1,j),1))
                  if (dd.gt.smax) smax = dd
               endif
210         continue
         endif
200   continue
c
      if (ipr.gt.1) then
         write(*,*)'*** Schmidt orthonormalisation'
         write(*,920)smax
920      format('*** Maximum overlap = ',e17.5)
      endif
      if (smax.gt.ortho) then
        write(*,900)smax,ortho
900     format(' ** Maximum overlap ',e12.5,' exceeds threshold ',e12.5,
     *          ' in Schmidt !!!')
        call caserr("Maximum overlap exceeds threshold !!!")
      endif
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine prmpschm(nref,nsing,ndoub,v,lv,nv,kept,cschm,ipr,
     1                    piv,rho)
c
      implicit real*8 (a-h,o-z)
c
c     Repeated Modified Gramm-Schmidt-orthonormalisation for MR-MP
c     ============================================================
c
c     Reference
c     ---------
c
c     Gerard L.G. Sleijpen
c
c
c     Purpose
c     -------
c
c     Orthonormalises the column vectors  v[1..nv] that are markted  
c     kept(i) = 1  using the modified Gramm-Schmidt method with 
c     pivotting.
c
c     nref  :
c        The number of reference functions
c
c     nsing :
c        The number of singly substituted functions
c
c     ndoub :
c        The number of doubly substituted functions
c
c     v     : 
c        Space spanned by <nv> column-vectors of length <lv>.
c        The basis of this space is returned in this variable.
c
c     kept  :
c        Vector of length <nv> indicating which vectors are
c        linearly independent and should be kept. 
c        If kept(i) = 1 the i-th vector will be kept.
c        If kept(i) = 0 the i-th vector will be dropped.
c
c     cschm : 
c        Criterion for orthogonality, i.e. two vectors are orthogonal 
c        if their inproduct is less than <cschm>.
c
c     ipr   : 
c        Print flag.
c        If ipr > 1 the accuracy of the orthonormalisation will be
c        checked and reported.
c
c     piv   :
c        Workspace to store the pivots.
c
c     rho   :
c        Workspace to store the original norms.
c
c
c     Theory
c     ------
c
c     Assume a vector  A  and  B  where the inproduct  <A|A> = 1  then
c     then the vector  B' = B - <A|B> A  is orthogonal to  A.
c     However, beware of numerical instability due to these projections
c     and subsequent normalisation.
c
c     When using pivots some efficiency can be gained through updating
c     these things. If  piv(B) = <B|B>  and  <A|A> = 1  then the 
c     piv(B') = < B-<A|B>A | B-<A|B>A >
c             = piv(B) - <A|B>**2
c
c
c     Algorithm by Gerard Sleijpen
c     ----------------------------
c
c     w te orthogonalizeren op V  (V is n bij k, kolom vektoren norm 1)
c     met resultaat v: w=V*beta+v, V'*v=0  (V' is de transpose van V)
c
c     kies kappa (bv kappa=0.5)
c
c
c     rho=||w||
c 
c     beta=V'*w
c     tw=w-V*beta   (dit is klassiek GS, modGS is ook prima)
c 
c     trho=||tw||
c 
c     if trho < 0.5 rho
c        tbeta=V'*tw
c        ttw=tw-V*tbeta
c        beta=beta+tbeta
c 
c        ttrho = ||ttw||
c 
c        if ttrho < 0.5 trho
c           v=0
c        else
c           v=ttw
c        end
c     else
c        v=tw
c     end
c 
c     Ervaringsfeit: hiermee worden fouten gemaakt van orde
c                    machine precisie/kappa
c 
c
c     Algorithm adapted to pivotting
c     ------------------------------
c
c     for i=ia, ib 
c        calculate and store the pivot of v[i]
c
c     for i=ia, ib
c        find vector v[i] in v[i..ib] with maximum pivot and 
c          move it to position i.
c        else if ( norm < 0.5 * stored norm ) then
c          repeat orthogonalisation to vectors v[ia..i-1]
c          calculate the norm again in norm1
c          if ( norm1 < 0.5 * norm ) then
c            drop vector v[i]
c          else
c            normalise vector[i]
c          endif
c        else
c          normalise vector[i]
c        endif
c
c        for j=i+1, ib
c           orthogonalise vector v[j] to vector v[i]
c
c     Notes
c     -----
c
c        - dmax1 calls required to avoid negative pivots due to limited
c          precision.
c
      real*8 dabs
      real*8 ddot
      external ddot
      parameter (zero   = 1.0d-15)
      parameter (rkappa = 0.5d0  )
c
      dimension v(lv,nv),kept(nv),piv(nv),rho(nv)
c
      xnorm  = 0.0d0
      ifase  = 0
c
c...  piv(i) = <v(i)|v(i)>
c
      do 100 i = 1, nv
         if (kept(i).eq.1) then
            piv(i) = dmax1(0.0d0, ddot(lv,v(1,i),1,v(1,i),1))
            rho(i) = dsqrt(piv(i))
         else 
            piv(i) = 0.0d0
            rho(i) = 0.0d0
         endif
100   continue
c
200   continue
      ifase = ifase + 1
      if (ifase .eq. 1) then
         ia = 1
         ib = nref
      else if (ifase .eq. 2) then
         ia = nref+1
         ib = nref+nsing
      else if (ifase .eq. 3) then
         ia = nref+nsing+1
         ib = nv
      else
         call caserr("subroutine prmpschm: This is insane!!!")
      endif
c
c...  v(j) = v(j) - <v(i)|v(j)> v(i) 
c...  using pivotting
c
      do 300 i = ia, ib
c
c...     find vector with max. pivot in current subspace
c
c...     pmax set to -1 for safety. Insures that if there is a kept
c...     vector left, it will be found.
c
         pmax = -1.0d0
         imax =  0
         do 310 j = i, ib
            if (kept(j).eq.1) then
               if (piv(i).gt.pmax) then
                 pmax = piv(i)
                 imax = j
               endif
            endif
310      continue
         if (imax.eq.0) then
c...        no more kept vectors in this subspace so leave loop
            goto 350
         endif
c
c...     move max. norm vector to current position
c
         kpttmp     = kept(i)
         kept(i)    = kept(imax)
         kept(imax) = kpttmp
c
         tmps       = piv(i)
         piv(i)     = piv(imax)
         piv(imax)  = tmps
c
         tmps       = rho(i)
         rho(i)     = rho(imax)
         rho(imax)  = tmps
c
         do 320 k = 1, lv
            tmps      = v(k,i)
            v(k,i)    = v(k,imax)
            v(k,imax) = tmps
320      continue
c
c...     if needed repeat orthogonalisation on current vector
c
         rho1 = dsqrt(piv(i))
         if (piv(i) .le. zero * rho(i) * rho(i)) then
c...        if pivot reduced by more than the machine precision then
c...        consider it zero.
            kept(i) = 0
         else if (rho1 .le. rkappa * rho(i)) then
            do 330 k = 1, i-1
               if (kept(k).eq.1) then
                  s = -ddot(lv,v(1,k),1,v(1,i),1)
                  call daxpy(lv,s,v(1,k),1,v(1,i),1)
                  piv(i) = dmax1(0.0d0, piv(i) - s * s)
               endif
330         continue
            rho2 = dsqrt(piv(i))
            if (piv(i) .le. zero * rho(i) * rho(i)) then
               kept(i) = 0
            else if (rho2 .le. rkappa * rho1) then
               kept(i)  = 0
               if (rho2.gt.xnorm) xnorm = rho2
            else
               kept(i)  = 1
            endif
         else
            kept(i)  = 1
         endif
c
c...     normalise current vector
c
         if (kept(i).eq.1) then
            dd = dsqrt(ddot(lv,v(1,i),1,v(1,i),1))
            if (dd.gt.zero) then
               dd = 1.0d0 / dd
               call dscal(lv,dd,v(1,i),1)
            else
              kept(i) = 0
            endif
         endif
c
c...     orthogonalise all remaining vectors on current vector
c
         if (kept(i).eq.1) then
            do 340 j = i+1, nv
               if (kept(j).eq.1) then
                  s = -ddot(lv,v(1,i),1,v(1,j),1)
                  call daxpy(lv,s,v(1,i),1,v(1,j),1)
                  piv(j) = dmax1(0.0d0, piv(j) - s * s)
               endif
340         continue
         endif
c
300   continue
c
350   continue
      if (ifase .lt. 3) goto 200
c
c...  find maximum overlap
c
      smax  = 0.0d0
      do 500 i = 1, nv
         if (kept(i).eq.1) then
            do 510 j = 1, i-1
               if (kept(j).eq.1) then
                  s = dabs(ddot(lv,v(1,j),1,v(1,i),1))
                  if (s.gt.smax) then
                     smax = s
                  endif
               endif
510         continue
         endif
500   continue
      if (ipr.gt.1.or.xnorm.lt.1.0d-10) then
         write(*,*)'*** PIVOTED REPEATED MODIFIED GRAMM-SCHMIDT'
         write(*,*)'Maximum norm dropped = ',xnorm
         write(*,*)'Maximum overlap      = ',smax
      endif
      if (smax.gt.cschm) then
        write(*,900)smax,cschm
900     format(' ** Maximum overlap ',e12.5,' exceeds threshold ',e12.5,
     *         ' in pivoted repeated modified Gramm-Schmidt !!!')
        call caserr("Maximum overlap exceeds threshold !!!")
      endif
c
      return
      end
c
c-----------------------------------------------------------------------
c
c
c-----------------------------------------------------------------------
c
      subroutine mphouse(nref,nsing,ndoub,v,lv,nv,kept,crit,ipr,q,h,w,
     1                   cpiv,ipiv)
c
      implicit real*8 (a-h,o-z)
c
c     Pivoted Householder QR orthonormalisation
c     =========================================
c
c     Reference
c     ---------
c
c     "Matrix Computations"
c      Gene H. Golub and Charles F. van Loan
c      Second Edition
c      The Johns Hopkins University Press, 1989
c      ISBN 0-8018-3772-3
c      ISBN 0-8018-3792-1 (pbk.)
c      Algorithm 5.4.1 (Householder QR with column pivoting)
c      Section   5.1.6 (Backward accumulation of qq(ivoff+)) 
c
c     Purpose
c     -------
c
c     Orthonormalises the column vectors  v[1..nv] that are markted  
c     kept(i) = 1  using Householder transformations with pivoting.
c
c     The algorithm given by Golub et al. has been adapted to 
c     guarantee that the space spanned by the singles will be 
c     orthogonalised to the reference functions, and the space spanned
c     by the doubles will be orthogonalised to the space of the 
c     references and singles. This is important because of the 
c     projection operators in the zero order Hamiltonian.
c
c     nref  :
c        The number of reference functions
c
c     nsing :
c        The number of singly substituted functions
c
c     ndoub :
c        The number of doubly substituted functions
c
c     v     : [1..lv,1..nv]
c        Space spanned by <nv> column-vectors of length <lv>.
c        The basis of this space is returned in this variable.
c
c     kept  : [1..nv]
c        Vector of length <nv> indicating which vectors are
c        linearly independent and should be kept. 
c        If kept(i) = 1 the i-th vector will be kept.
c        If kept(i) = 0 the i-th vector will be dropped.
c
c     ipr   : 
c        Print flag. 
c        If ipr > 1 the accuracy of the orthonormalisation will be
c        checked and reported.
c
c     q     : [1..lv,1..lv]
c        Workspace to accumulate the orthogonal space.
c
c     h     : [1..lv]
c        Workspace to store Householder vectors.
c
c     w     : [1..max(lv,nv)]
c        Workspace to store temporary vectors.
c
c     cpiv  : [1..nv]
c        Workspace to store the pivots.
c
c     ipiv  : [1..nv]
c        Workspace to store the vector numbers to administer the
c        permutations due to the pivoting.
c
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      real*8 dabs
      real*8 ddot
      external ddot
      real*8  dd, s, smax
      real*8  cpiv(nv)
      real*8  h(lv), w(nv),q(lv,lv)
      integer ipiv(nv)
c
      parameter (zero   = 1.0d-15)
c
      dimension v(lv,nv),kept(nv)
c
c...  Needed to figure out sign of reference vector.
c...  Householder may change the sign, because it only does an 
c...  orthonormalisation. However, a change of sign leads to 
c...  problems in the MP3.
c
      crtmp  = 1.0d0
      irtmp  = 0
      smin   = 1.0d15
c
c...  Smin will become the minimal selected pivot. In case of
c...  trouble this can point out if the criterion was set to
c...  thight.
c
      if (nref.ne.0) then
         do i = 1, lv
            if (dabs(v(i,1)).gt.zero) then
               irtmp = i
               crtmp = v(i,1)
               goto 20
            endif
         enddo
20       continue
      endif
c
c...  inproduct(i) = <v(i)|v(i)> setting up the pivots
c
      do i = 1, nv
         if (kept(i).eq.1) then
            cpiv(i) = ddot(lv,v(1,i),1,v(1,i),1)
         else
            cpiv(i) = 0.0d0
         endif
      enddo
      call ibasgn(nv,1,1,ipiv)
c
c...  Starting the Householder QR factorisation.
c...  No linear independent vectors found yet (nr.eq.0)
c
      ifase  = 0
      nr     = 0
c
200   continue
      ifase = ifase + 1
      if (ifase .eq. 1) then
         ia = 1
         ib = nref
      else if (ifase .eq. 2) then
         ia = nref+1
         ib = nref+nsing
      else if (ifase .eq. 3) then
         ia = nref+nsing+1
         ib = nv
      else
         call caserr("subroutine mphouse: This is insane!!!")
      endif
c
c...  find vector with max. pivot in current subspace
c
c...  smax set to -1 for safety. Insures that if there is a kept
c...  vector left, it will be found.
c
      smax = -1.0d0
      imax =  0
      do i = ia, ib
         if (cpiv(i).gt.smax) then
            smax = cpiv(i)
            imax = i
         endif
      enddo
c
c...  begin while loop of the orthogonalisation
c
400   if (smax.gt.max(zero,crit)) then
         nr       = nr + 1
         if (smax.lt.smin) smin = smax
c
c...     swap pivots and vectors
c
         itmp       = ipiv(nr)
         ipiv(nr)   = ipiv(imax)
         ipiv(imax) = itmp
c
         itmp       = kept(nr)
         kept(nr)   = kept(imax)
         kept(imax) = itmp
c
         tmp        = cpiv(nr)
         cpiv(nr)   = cpiv(imax)
         cpiv(imax) = tmp
c
         do it = 1, lv
            tmp        = v(it,nr)
            v(it,nr)   = v(it,imax)
            v(it,imax) = tmp
         enddo
c
c...     compute householder vector in h(i)
c
         call house(nr,lv,h(1),v(1,nr))
c
c...     perform householder transformation
c
         call rowhouse(lv,nv,nr,v,lv,h,w)
c
c...     copy householder vector into matrix
c
         do it = nr+1, lv
            v(it,nr) = h(it)
         enddo
c
c...     update the pivots
c
         do it = nr+1,  nv
            if (cpiv(it).gt.0.0d0) then
               cpiv(it) = cpiv(it) - v(nr,it)*v(nr,it)
            endif
         enddo
c
c...     find vector with max. pivot in current subspace
c
         if (max(nr,ia) .lt. ib) then
            smax = -1.0d0
            imax =  0
            do i = max(nr+1,ia), ib
               if (cpiv(i).gt.smax) then
                  smax = cpiv(i)
                  imax = i
               endif
            enddo
            goto 400
         endif
      else
         do i = max(nr+1,ia), ib
            cpiv(i) = 0.0d0
         enddo
      endif
c
c...  if not all subspaces have been treated,
c...  go for the next subspace.
c
      if (ifase .lt. 3) goto 200
c
c...  Othogonalisation is done now accumulate qq(ivoff+) from Householder vectors
c
      do 500 it = 1, lv
         do 510 jt = 1, lv
            q(jt,it) = 0.0d0
510      continue
         q(it,it) = 1.0d0
500   continue
      do 520 it = 1, nr
         i = nr+1-it
         do 530 k = i+1, lv
            h(k) = v(k,i)
530      continue
         h(i) = 1.0d0
         call rowhouse(lv,lv,i,q,lv,h,w)
520   continue
c
c...  copy back qq(ivoff+) into V
c
      do 600 i = 1, nv
         kept(i) = 0
600   continue
      do 610 i = 1, nr
         ic = ipiv(i)
         do 620 it = 1, lv
            v(it,ic) = q(it,i)
620      continue
         kept(ic) = 1
610   continue
c
c...  normalize vectors
c...  if the are reference vectors then restore the proper sign
c
      if (nref.ne.0) then
        crtmp = dsign(1.0d0,crtmp*v(irtmp,1))
      endif
      do 700 i = 1, nv
        if (kept(i).eq.1) then
           dd = dsqrt(ddot(lv,v(1,i),1,v(1,i),1))
           dd = crtmp / dd
           call dscal(lv,dd,v(1,i),1)
        endif
700   continue
c
c...  if printing on then find and print maximum overlap
c
      smax  = 0.0d0
      do 800 i = 1, nv
         if (kept(i).eq.1) then
            do 810 j = 1, i-1
               if (kept(j).eq.1) then
                  s = dabs(ddot(lv,v(1,j),1,v(1,i),1))
                  if (s.gt.smax) then
                     smax = s
                  endif
               endif
810         continue
         endif
800   continue
      if (ipr.gt.100.or.smin.lt.1.0d-10) then
         write(iwr,*) ' PIVOTED HOUSEHOLDER orthonormalisation'
         write(iwr,*) ' Minimal pivot selected = ',smin
         write(iwr,*) ' Maximum overlap        = ',smax
      endif
      if (smax.gt.crit) then
         write(iwr,*) ' Minimal pivot selected = ',smin
         write(iwr,900)smax,crit
900      format(1x,'*** Maximum overlap ',e12.5,' exceeds threshold ',
     *           e12.5,
     *          ' in Householder QR orthonormalisation !!!')
         call caserr("Maximum overlap exceeds threshold !!!")
      endif
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine house(ia,n,v,x)
c
      implicit real*8 (a-h,o-z), integer(i-n)
      dimension v(n), x(n)
c
c     Computes a Householder vector
c     =============================
c
c     Given an vector x(ia:n), this function computes an vector 
c     with v(ia) = 1 such that 
c     (I(ia:n,ia:n) - 2 v(ia:n)v(ia:n)T/(v(ia:n)Tv(ia:n)))x(ia:n) 
c     is zero in all but the first component.
c
c     Reference
c     ---------
c
c     "Matrix Computations"
c      Gene H. Golub, Charles F. van Loan
c      Second Edition
c      The Johns Hopkins University Press, 1989
c      ISBN 0-8018-3772-3
c      ISBN 0-8018-3792-1 (pbk.)
c      Algorithm 5.1.1
c
      u = 0.0d0
      do 90  i = 1,  ia-1
        v(i) = 0.0d0
90    continue
      do 100 i = ia, n
        u    = u + x(i) * x(i)
        v(i) = x(i)
100   continue
      u = dsqrt(u)
c
      if (u.ne.0.0d0) then
         b = x(ia) + dsign(u,x(ia))
         b = 1.0d0 / b
         call dscal(n-ia,b,v(ia+1),1)
      endif
      v(ia) = 1
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine rowhouse(m,n,ir,a,ncdim,v,w)
c
      implicit real*8 (a-h,o-z), integer(i-n)
      dimension a(ncdim,n), v(m), w(n)
c
c     Performs a Householder transformation
c     =====================================
c
c     "Matrix Computations"
c      Gene H. Golub, Charles F. van Loan
c      Second Edition
c      The Johns Hopkins University Press, 1989
c      ISBN 0-8018-3772-3
c      ISBN 0-8018-3792-1 (pbk.)
c      Algorithm 5.1.2
c
      b = -2.0d0 / ddot(m-ir+1,v(ir),1,v(ir),1)
      do 10 i = ir, n
         w(i) = 0.0d0
10    continue
c
c...  w = b * transpose(A) * v
c
      do 110 j = ir, n
         do 100 i = ir, m
            w(j) = w(j) + A(i,j) * v(i)
100      continue
         w(j) = w(j) * b
110   continue
c
c...  A = A + v * transpose(w)
c
      do 120 i = ir, n
         do 130 j = ir, m
            A(j,i) = A(j,i) + v(j) * w(i)
130      continue
120   continue
c
      end
c
c=======================================================================
c
      subroutine prfmat(a,ndim,n,m,kept)
      implicit none
      integer n,m,ndim
      integer kept(m)
      real*8  a(ndim,m)
c
c...  prints the columns of a rectangular matrix a that are marked 
c...  with kept = 1. It is assumed that all other columns are marked
c...  with kept = 0.
c
      integer mxcol
      parameter (mxcol = 10)
c
      integer nkept, ncblk, icblk, irow, icol, i, itmp
c
c     print*
c     print *,'** ndim,n,m = ',ndim,n,m
c     print*,kept
c
      nkept = 0
      do 10 i = 1, m
        if (kept(i).eq.1) then
          nkept = nkept + 1
          kept(nkept) = i
        endif
10    continue
c
      ncblk = (nkept + mxcol - 1) / mxcol
c
      write(*,*)
      do 20 icblk = 1, ncblk-1
        write(*,900)(kept(i),i=(icblk-1)*mxcol+1,icblk*mxcol)
        write(*,*)
        do 30 irow = 1, n
          write(*,910)irow,(a(irow,kept(i)),
     *                      i=(icblk-1)*mxcol+1,icblk*mxcol)
30      continue
        write(*,*)
20    continue
c
      write(*,900)(kept(i),i=max(1,(ncblk-1)*mxcol+1),nkept)
      write(*,*)
      do 40 irow = 1, n
        write(*,910)irow,(a(irow,kept(i)),
     *                    i=max(1,(ncblk-1)*mxcol+1),nkept)
40    continue
      write(*,*)
c
      do 50 itmp = 1, nkept
        i = nkept + 1 - itmp
        icol       = kept(i)
        kept(i)    = 0
        kept(icol) = 1
50    continue
c     print*,kept
c
900   format(6x,10i14)
910   format(i6,10f14.6)
c
      end
c
c=======================================================================
c
      logical function inclass(ms,ns,mstart,mlast,nstart,nlast)
c
c...  sorts out if two configurations belong to the same excitation 
c...  class.
c
      implicit real*8 (a-h,o-z)
c
c
      integer nref, nnst, nmin1, nmin2, mspin, nelec
      integer nstrt1, nlast1, nstrt2, nlast2, nmin12, nminv2, nminv1
      integer nlst21, nlastx
      common /ccntl/ nref,nnst,nmin1,nmin2,mspin,nelec,
     +               nstrt1,nlast1,nstrt2,nlast2,nmin12,nminv2,nminv1,
     +               nlst21,nlastx
c
c
      integer nval, ndoub, nsingl, ntripl, ntotal, ibig, jbig, kbig
      integer nvd, nvds, nsitri, lbig, mbig, nbig, nref0
      common /cntrl/ nval,ndoub,nsingl,ntripl,ntotal,ibig,jbig,kbig,
     +               nvd,nvds,nsitri,lbig,mbig,nbig,nref0
c
chvd
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
      common/mr_mp2l/ ivoc(62),nvoc,idoc(62),ndoc,kvcc(4),nvcc(4),
     *                ncmphn(10,18),ncmph1(10,18,8),ncmph2(10,18,8),
     *                ncmphd(10,18)
      common /mr_mp2con/ nsmpv,nsmp1,nsmp2,nsmp2d,ntotmp,
     *                   nmpv,nmp1,nmp2,nmp2d
c
      integer nirr,mult,ifrep,nirr8,jfrep
      common /table/ nirr,mult(8,8),ifrep,nirr8,jfrep
c
c
      integer nint, next, ninnex, nwm1, nw, nwp1, nintwo
      integer nint21, nintp1
      common /orb/ nint,next,ninnex,nwm1,nw,nwp1,nintwo,
     +             nint21,nintp1
c
chvd
c...  find out number of holes.
c
      if (ms.eq.1) then
         mhole = 0
      else if (max(2,ms).eq.min(ms,nirr+1)) then
         mhole = 1
      else if (max(nirr+2,ms).eq.min(ms,2*nirr+2)) then
         mhole = 2
      else
         call caserr('illegal ms in inclass')
      endif
c
      if (ns.eq.1) then
         nhole = 0
      else if (max(2,ns).eq.min(ns,nirr+1)) then
         nhole = 1
      else if (max(nirr+2,ns).eq.min(ns,2*nirr+2)) then
         nhole = 2
      else
         call caserr('illegal ns in inclass')
      endif
c
c...  find out number of particles.
c
      if (max(0,mstart).eq.min(mstart,nstrt1-1).and.
     *    max(0,mlast) .eq.min(mlast,nstrt1-1)) then
         mpart = 0
      else if (max(nstrt1,mstart).eq.min(mstart,nlast1).and.
     *         max(nstrt1,mlast) .eq.min(mlast,nlast1)) then
         mpart = 1
      else if (max(nstrt2,mstart).eq.min(mstart,nlast2).and.
     *         max(nstrt2,mlast) .eq.min(mlast,nlast2)) then
         mpart = 2
      else
         call caserr('illegal mstart and/or mlast in inclass')
      endif
c
      if (max(0,nstart).eq.min(nstart,nstrt1-1).and.
     *    max(0,nlast) .eq.min(nlast,nstrt1-1)) then
         npart = 0
      else if (max(nstrt1,nstart).eq.min(nstart,nlast1).and.
     *         max(nstrt1,nlast) .eq.min(nlast,nlast1)) then
         npart = 1
      else if (max(nstrt2,nstart).eq.min(nstart,nlast2).and.
     *         max(nstrt2,nlast) .eq.min(nlast,nlast2)) then
         npart = 2
      else
         call caserr('illegal nstart and/or nlast in inclass')
      endif
c
c...  find out if the configurations are in the same class
c
      inclass = (mpmod.ge.0).or.(mhole.eq.nhole.and.mpart.eq.npart)
c
      return
      end

c*****************MR-MP FOCK-matrix SORT*****************************
      subroutine mpgij(fock,fia,fabtr,fabsq,
     *                 conf,ref,cr,lcr)
c
      implicit real*8 (a-h,o-z)
      dimension fock(*),fia(*),fabtr(*),fabsq(*)
      integer conf
      dimension conf(5,*),ref(*),cr(lcr)
      external fget
c
c...   read fock-matrix from section isec_mp on dumpfile (nummp,iblk_mp)
c...   if nummp = -999 use vectors section orbital energies as fock
c...   derived from gijkl
c...   write ia and ab parts to ed7
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
      integer iro, iky, niky, m1e, m2e, m1ia, m2coul
      integer m2exch, m2fock, m2psup, m2qsup
      integer mbase, nbase, npairi, norbi, mbasi, mnbasi
      integer iap, mbfock, mbexch, mbcoul, mbpsup, mbqsup
      integer norbas, m2exab, m2em1e
      common /helpr/ iro(nd200),iky(nd200),niky(8),m1e,m2e,m1ia,m2coul,
     +   m2exch,m2fock,m2psup,m2qsup,
     +   mbase(288),nbase(8),npairi(36),norbi(8),mbasi(8),mnbasi(8),
     +   iap(62),mbfock(62),mbexch(36),mbcoul(36),mbpsup(36),mbqsup(36),
     +   norbas(nd200),m2exab,m2em1e
c
c
       integer mapei, mapie
       logical modei
       common/mapp/mapei(nd200),modei,mapie(nd200)
c
c
      integer nint, next, ninnex, nwm1, nw, nwp1, nintwo
      integer nint21, nintp1
      common /orb/ nint,next,ninnex,nwm1,nw,nwp1,nintwo,
     +             nint21,nintp1
c
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
c
      integer nxblk, index, numci, mblkci, ispbig, ispsma
      common/disktl/nxblk,index(199),numci,mblkci,ispbig,ispsma
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
c
      integer nirr,mult,ifrep,nirr8,jfrep
      common /table/ nirr,mult(8,8),ifrep,nirr8,jfrep
c
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer isec3, nsect2, notap2, iblk2, lblk2
      common /file12/ isec3,nsect2,notap2(20),iblk2(20),lblk2(20)
c
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
      common/junkc/zjob,zdate,ztime,zprog,ztype,zspace(14),ztext(10)
      common/junk/
     * evalue(maxorb),eocc(maxorb+1),
     * nbas_1,newb_1,ncoll,ivalue,iocc
      common/lsort/ilifx(maxorb),ilifs(maxorb),ilifq(maxorb),
     + lelim(maxorb),comman(maxorb),pottt(2),ncolo,nbas,newb,ncore,
     + mapcie(maxorb),ilifc(maxorb),nsact,mapaie(maxorb),
     + ilifa(maxorb),iqsec 
      common/blkin/value(510),mmword,mdum,mult8(8,8)
      dimension commen(5),occ(255),tit(10)
      dimension rjunk(512)
      equivalence(commen(1),rjunk(1))
      equivalence(occ(1),rjunk(6))
      equivalence(tit(1),rjunk(261))
      equivalence(potn,rjunk(271))
      equivalence(core,rjunk(272))
      itr(i) = i*(i+1)/2
c
c... 1-e ij in core
c... 1-e ia in core     straight triangle
c... 1-e ab in core
c
c     indmp : 1 coupling coeff. / 2(ia) /  / 4 (ab) square / 5 (aa) diag
c
      if (nummp.eq.0) then
c
c...   no fock-directive presented so calculate fock-matrix
c...   from the density matrix computed from the reference function ref
c
         call mpfock(conf,ref,fock,cr,lcr,valmax,ivmax,jvmax)
         go to 500
      else if (nummp.eq.-999) then
c...     pick orbital info up from isec3
c...     cf. symsap
         call secget(isec3,1004,jblkk)
         call rdedx(comman,mach(15),jblkk,idaf)
c....    now go to vectors section
         call secget(iqsec,3,jblkk)
         call rdchr(zjob,29,jblkk,idaf)
         call reads(evalue,mach(8),idaf)
         write(iwr,1101) iqsec,ztype,zdate,ztime,zjob
1101  format(
     *   //' Diagonal FOCK-matrix  restored from vectors section',i4,
     *   /,1x,a8,' section created on ',a8,' at ',a8,' in the job ',a8)
         if (ivalue.ne.1) call caserr('no eigenvals on vectors section')
c
c...     get active eigenvalues
c
         do i=1,nsact
           eocc(i) = evalue(mapaie(i))
         end do
c
c...     get in DIRECT order
c
         call vclr(fock,1,ninnex*(ninnex+1)/2)
         do i=1,ninnex
           fock(itr(mapei(i))) = eocc(i)
         end do
c
         go to 500
      end if
c
      call cpuwal(top,bot)
c
      valmax = 0.0d0
      ivmax = 0
      nwb = 0
      do 2 loop=1,nirr
2     nwb=iky(norbi(loop)+norbe(loop)+1)+nwb
c
      call vclr(fock,1,ninnex*(ninnex+1)/2)
c...
c... get 1-elec FOCK integrals
c...
      call secini(iblk_mp,nummp)
      call secget(isec_mp,1004,jblkk)
      call rdedx(commen,272,jblkk,nummp)
      write(6,1100)top,bot,isec_mp,iblk_mp,yed(nummp),commen(2),
     *commen(3),commen(4),commen(1),tit
1100  format(//1x,80('-')//
     *' **** FOCK matrix sorter v1.1 called at',f12.1,' wall'
     *,f12.1,' secs'
     *////' FOCK-matrix integrals restored from section',i4,
     *' of dumpfile starting at block',i8,' of ',a4/
     *' section created on ',a8,' at ',a8,' in the job ',a8,
     *' under acct ',a8/' with the title ',10a8)
      jblkk=jblkk+3
c
      call search(jblkk,nummp)
      call fget(value,kword,nummp)
      nin = 0
c
      do 101 iii=1,200
         iiii=mapei(iii)
         do 101 jjj=1,iii
            j=mapei(jjj)
            i=iiii
            if(i.lt.j) then
               i=j
               j=iiii
            end if
            nin=nin+1
            if (nin.gt.kword) then
               call fget(value,kword,nummp)
               nin=1
            end if
c
            if (i.gt.ninnex.or.j.gt.ninnex) go to 101
c
            if(iro(j).ne.iro(i)) then
c
c...   symmetry forbidden integral
c
               if (dabs(value(nin)).gt.dabs(valmax)) then
                  valmax = value(nin)
                  ivmax = i
                  jvmax = j
               end if
            else
c
c... (i/i) // (i/a) // (a/a) in fock  no symmetry used (cf. gijkl)
c
               fock(i*(i-1)/2+j)=value(nin)
               nwb=nwb-1
               if (nwb.eq.0) go to 10102
            end if
c
101   continue
c
c...   reset dumpfile
c
      call secini(ibl3d,idaf)
c
500   continue
c...
c...   symmetry-check print
c...
10102   if (ivmax.ne.0) write(6,512) ivmax,jvmax,valmax
512     format(/' *** symmetry check / largest forbidden F-integral :',
     1         2i4,e12.5)
      if (iprimp.ge.50) then
         write(6,601)
601      format(/' === FOCK-matrix as read ===')
         call prtri(fock,ninnex)
      end if
c
c...  put too small fock-integrals to 0.0
c
      kk = 0
      valmax = 0.0d0
      ifdiaf = 1
      do 600 i=1,ninnex
        do 600 j=1,i
          kk = kk + 1
          if (fock(kk).ne.0.0d0.and.dabs(fock(kk)).lt.fignor) then
             if (dabs(fock(kk)).gt.dabs(valmax)) then
                valmax = fock(kk)
                ivmax = i
                jvmax = j
             end if
             fock(kk) = 0.0d0
          end if
          if (j.gt.nint.and.i.ne.j.and.fock(kk).ne.0.0d0) ifdiaf = 0
600   continue
      if (valmax.ne.0.0d0) write(6,513) ivmax,jvmax,valmax
513     format(/' *** largest ignored  F-integral :',
     1         2i4,e12.5,' ***')
      if (ifdiaf.eq.0) then
         write(6,514)
514      format(' *** external fock-matrix is **not** diagonal ***')
      else
         write(6,515)
515      format(' *** external fock-matrix is **diagonal** ***')
      end if
c
c...  ****** process FOCK(ia) *******  => indmp(2)
c
      call search(nxblk,numci)
c
c...  gather (ia) integral
cs
      kki = 0
      kke = nint
      kkf = 0
      do 1000 irr=1,nirr
         ni = norbi(irr)
         ne = norbe(irr)
         do 900 j=kki+1,kki+ni
            do 800 i=kke+1,kke+ne
               kkf = kkf + 1
               fia(kkf) = fock(i*(i-1)/2+j)
800         continue
900      continue
         kki = kki + ni
         kke = kke + ne
1000  continue
c
c...   write
c
      indmp(2) = iposun(numci)
      call wrt3(fia,m1ia,indmp(2),numci)
c
      if (iprimp.ge.1500) write(6,602) (fia(i),i=1,m1ia)
602      format(/' Fock(ia) integrals ',/,(2x,10f11.5))
c
c...  ****** process FOCK(ab) *******  => indmp(3) / indmp(4)
c
c
c...  to form fock(ab) operators for doublets +  n-2 states
c...   adressing : norbsi(1) : length of triangular vector
c...               ibassi(irep,1) : pos. in triangular (start at 0)
c...               norbsq(1) : length of squared out vector
c...               ibassq(irep,1) : pos. in square (start at 0)
c...               norbe : dimension/symmetry (of course)
c
c...  gather triangle ab integrals
c
      kkf = 0
      kke = nint
      do 2000 irr=1,nirr
         ne = norbe(irr)
         do 1900 i=kke+1,kke+ne
          do 1900 j=kke+1,i
            kkf = kkf + 1
            fabtr(kkf) = fock(i*(i-1)/2+j)
1900     continue
         kke = kke + ne
2000  continue
c...
c... fock(ab)  doublets/diagonal n-2 (square for ddmpab/tri for ssiaib)
c...
c      if (iprimp.ge.1000) write(6,603) (fabtr(i),i=1,norbsi(1))
c603      format(/' Fock(ab) triangles ',/,(2x,10f11.5))
c
        indmp(3) = -1
c       indmp(3) = iposun(numci)
c       call wrt3(fabtr,norbsi(1),indmp(3),numci)
c...
c... fock(ab)  for off-diagobnal n-2  / squared most convenient
c...
c...  square the operator
c
      do 3000 irr=1,nirr
         call square(fabsq(ibassq(irr,1)+1),fabtr(ibassi(irr,1)+1),
     *               norbe(irr),norbe(irr))
3000  continue
c
c...   write
c
       if (iprimp.ge.1000) write(6,604) (fabsq(i),i=1,norbsq(1))
604    format(/' Fock(ab) squares ',/,(2x,10f11.5))
c
      indmp(4) = iposun(numci)
      call wrt3(fabsq,norbsq(1),indmp(4),numci)
c
c...  gather fock(aa) diagonal
c
      do 4000 i=nint+1,ninnex
4000  fabtr(i-nint) = fock(i*(i+1)/2)
c
      indmp(5) = iposun(numci)
      call wrt3(fabtr,next,indmp(5),numci)
c
      nxblk = iposun(numci)
c

      return
      end
      subroutine mpfock(conf,ref,fock,cr,lcr,valmax,ivmax,jvmax)
c...
c... build fock-operator for MRMP reference function
c... first part derived form nosas (build density matrix)
c... second part from tran (build fock-operator)
c...
      implicit real*8 (a-h,o-z)
c
      integer  conf
c
      dimension conf(5,*),ref(*),fock(*),cr(lcr)
      external fget
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
      integer isec3, nsect2, notap2, iblk2, lblk2
      common /file12/ isec3,nsect2,notap2(20),iblk2(20),lblk2(20)
c
      common/craypk/i205(340),j205(340),k205(340),l205(340) 
      common/blkin/value(510),mmword,mdum,mult8(8,8)
c
       integer mapei, mapie
       logical modei
       common/mapp/mapei(nd200),modei,mapie(nd200)
c
c
      integer iro, iky, niky, m1e, m2e, m1ia, m2coul
      integer m2exch, m2fock, m2psup, m2qsup
      integer mbase, nbase, npairi, norbi, mbasi, mnbasi
      integer iap, mbfock, mbexch, mbcoul, mbpsup, mbqsup
      integer norbas, m2exab, m2em1e
      common /helpr/ iro(nd200),iky(nd200),niky(8),m1e,m2e,m1ia,m2coul,
     +   m2exch,m2fock,m2psup,m2qsup,
     +   mbase(288),nbase(8),npairi(36),norbi(8),mbasi(8),mnbasi(8),
     +   iap(62),mbfock(62),mbexch(36),mbcoul(36),mbpsup(36),mbqsup(36),
     +   norbas(nd200),m2exab,m2em1e
c
c
      integer nconst, ncond, maxpa, maxpat, maxsin, maxepa
      integer mxdvpa, mxddpa, nspips, nspipt, mxsspa, mxttpa
      integer mxstpa, nsitot, ntrtot
      integer mx12pa, mxvvpa, mxv2pa, maxpsq, max12
      common /auxh/ nconst(8),ncond(8),maxpa(4),maxpat,maxsin(3),
     +              maxepa(4),mxdvpa,mxddpa,nspips(8),nspipt(8),
     +              mxsspa,mxttpa,mxstpa,nsitot(8),ntrtot(8),
     +              mx12pa,mxvvpa,mxv2pa,maxpsq,max12
c
c
      integer nint, next, ninnex, nwm1, nw, nwp1, nintwo
      integer nint21, nintp1
      common /orb/ nint,next,ninnex,nwm1,nw,nwp1,nintwo,
     +             nint21,nintp1
c
c
      integer nuse, ngot
      common /corctl/ nuse, ngot
c
c     common/spew/model,ic,ir,ninti
      common/spew/rmodel(5)
      equivalence (rmodel(1),model),(rmodel(2),ic),(rmodel(3),ir),
     1            (rmodel(4),ninti)
c
      integer nxblk, index, numci, mblkci, ispbig, ispsma
      common/disktl/nxblk,index(199),numci,mblkci,ispbig,ispsma
c
      integer nref, nnst, nmin1, nmin2, mspin, nelec
      integer nstrt1, nlast1, nstrt2, nlast2, nmin12, nminv2, nminv1
      integer nlst21, nlastx
      common /ccntl/ nref,nnst,nmin1,nmin2,mspin,nelec,
     +               nstrt1,nlast1,nstrt2,nlast2,nmin12,nminv2,nminv1,
     +               nlst21,nlastx
c
      common/erg/potnuc,totnuc,energy
      equivalence (bilbo,ir)
c
c
c...   triangle-adressing   iijj
c
      iijj(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)
c
      call cpuwal(top,bot)
      write(iwr,1100)top,bot
1100  format(//1x,80('-')//
     *' **** FOCK matrix calc. v1.1 called at',f12.1,' wall'
     *,f12.1,' secs')
c
c...  densitry matrix is at start of cr
c
      icl = -1
      nd = nint*(nint+1)/2
      call vclr(cr,1,nd)
      khe = nd
      khf = nd + maxpa(4)**2
      if (khf+maxpa(4).gt.lcr) call caserr('core-overflow1 in mpfock')
c
      call brtsb(index(197))
1     call getsb(rmodel(1),4)
      goto (99,2,3,99,99,99,99,99,99,10,99,99,99,99,99,99),model
c
c...  we only need vacuum-vacuum  (within reference set)
c
c...
c...  vacuum/vacuum   (ij)
c...
2     continue
      lamc = conf(3,ic)
      lamr = conf(3,ir)
c
      if (ic.gt.nref.or.ir.gt.nref) then
         call skipsb(lamc*lamr)
         go to 1
      end if
c
      call getsb(cr(khe+1),lamc*lamr)
c
      icp = conf(1,ic)
      irp = conf(1,ir)
c
      call mxmd(cr(khe+1),1,lamr,ref(icp+1),1,1,
     *          cr(khf+1),1,1,lamr,lamc,1)
      cr(ninti) = ddot(lamr,cr(khf+1),1,ref(irp+1),1) + cr(ninti)
c
      goto 1
c...
c... vacuum/vacuum   (ii)
c...
3     if (ic.gt.nref) go to 1
c
      if  (ic.ne.icl) then
         jc = conf(1,ic)
         nc = conf(3,ic)
         frodo = ddot(nc,ref(jc+1),1,ref(jc+1),1)
         icl=ic
      end if
c
      cr(ninti) = frodo*bilbo + cr(ninti)
c
      goto 1
c
c...  spin vacuum ij/ii ... skip
c
10    continue
      lamc = conf(3,ic)
      lamr = conf(3,ir)
      call skipsb(lamc*lamr)
      go to 1
c
c...  symbolic gave corrections with respect to closed-shell
c
99    kk = 0
      do 100 i=1,nint
         kk = kk + i
         cr(kk) = cr(kk) + 2.0d0
100   continue
c
c     translate (permute) matrix back to original mo-order
c     first find the maximum index of a non-zero d-matrix element
c     and the maximum size of the fock-matrix in original order
c
      maxdi = 0
      maxde = 0
      do 110 i=1,ninnex
         if (i.le.nint) maxdi = max(maxdi,mapie(i))
         maxde = max(maxde,mapie(i))
110   continue
c
      lenbas = maxde*(maxde+1)/2
      kh = lenbas
      kf = kh + lenbas
      if (kf+lenbas.gt.lcr) call caserr('core overflow2 in mpfock')
c
      call dcopy(nd,cr,1,cr(lenbas+1),1)
      call vclr(cr,1,lenbas)
      do 120 i=1,nint
        do 120 j=1,nint
           ifr = iijj(i,j)
           ito = iijj(mapie(i),mapie(j))
           cr(ito) = cr(lenbas+ifr)
120   continue
c
c...  read the h (1-electron) matrix  (straight)
c
      call vclr(cr(kh+1),1,lenbas)
c...
      call secget(isec3,1004,jblkk)
      jblkk=jblkk+lensec(mach(15))
c
      call search(jblkk,idaf)
c
c...  integrals are on disk as straight triangle so read straight
c
      kk = 0
130      call fget(value,kword,idaf)
         nn = min(kword,lenbas-kk)
         call dcopy(nn,value,1,cr(kh+kk+1),1)
         kk = kk + nn
      if (kk.lt.lenbas) go to 130
c
c...  now build fock-matrix (using ed6)
c
c...  scale density matrix *.5 and diagonal again *.5
c
      call dscal(lenbas,0.5d0,cr,1)
      kk = 0
      do 140 i=1,maxdi
         kk = kk + i
         cr(kk) = cr(kk)*0.5d0
140   continue
c
c...   top will contain d*(h+f) (commented out)
c
c      top=ddot(lenbas,cr,1,cr(kh+1),1)
      call setsto(1360,0,i205) 
      call vclr(cr(kf+1),1,lenbas)
      do 200 i=1,nsect2
         lbl = iblk2(i) - lblk2(i)
         if (lbl.eq.0) go to 200
           iunit = notap2(i)
           call search(iblk2(i),iunit)
150         call fget(value,k,iunit)
            if (k.eq.0) go to 200
            call sgmat_mp(cr(kf+1),cr,maxdi,maxde)
            lbl = lbl+1
            if (lbl.ne.0) go to 150
200    continue
c
      call vadd(cr(kf+1),1,cr(kh+1),1,cr(kh+1),1,lenbas)
c      top = ddot(lenbas,cr,1,cr(kh+1),1) + top
c      top = top + top + potnuc
c
c     we now have the fock-matrix in cr(kh)
c     symmetry sort it back
c
      call vclr(fock,1,ninnex*(ninnex+1)/2)
      valmax = 0.0d0
      ivmax = 0
c
      do 220 i=1,ninnex
        do 220 j=1,ninnex
           ifr = iijj(mapie(i),mapie(j))
           if (iro(i).ne.iro(j)) then
              if (dabs(cr(kh+ifr)).gt.dabs(valmax)) then
                 valmax = cr(kh+ifr)
                 ivmax = i
                 jvmax = j
              end if
           else
              ito = iijj(i,j)
              fock(ito) = cr(kh+ifr)
           end if
220   continue
c
      return
      end
      subroutine sgmat_mp(fock,p,maxi,maxe)
c*******************************************************************
c**             sgmat      taken from scf                         **
c**     requires a few changes (back) in calling routines         **
c**    originally crayscs (old atmol cray scf - 1982)             **
c**    modernised and adapted to  i1,j1,k1,l1,i2,j2,k2,l2 pack    **
c**    for nos/ve  / joop van lenthe   january 1987               **
c**    adapted for speciasl use in DIRECTMP jvl Jun89 (maxi/maxe) **
c*******************************************************************
      implicit real*8 (a-h,o-z)
      dimension p(*),fock(*)
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      common/blkin/gg(510),mword,ndum
c*********
      common /craypk/ i205(340),j205(340),k205(340),l205(340)
c
      iky(i) = i*(i-1)/2
c
      call unpack(gg(num2e+1),lab816,i205,numlab)
c*********
      do 6000 iw=1,mword
              int4 = (iw-1)*4+1
              j = i205(int4  )
              i = i205(int4+1)
              l = i205(int4+2)
              k = i205(int4+3)
cmp
cmp.     direct // we know there are'nt any d-matrix elemnts
cmp.     for orbitals outside of internal set
cmp
         if ((l.gt.maxi.and.j.gt.maxi).or.i.gt.maxe) go to 6000
cmp
c*********
         gik=gg(iw)
         g2=gik+gik
         g4=g2+g2
         ikyi=iky(i)
         ikyj=iky(j)
         ikyk=iky(k)
         ik=ikyi+k
         il=ikyi+l
         ij=ikyi+j
         jk=ikyj+k
         jl=ikyj+l
         kl=ikyk+l
         aij=g4*p(kl)+fock(ij)
         fock(kl)=g4*p(ij)+fock(kl)
         fock(ij)=aij
c... exchange
         gil=gik
         if(i.eq.k.or.j.eq.l)gik=g2
         if(j.eq.k)gil=g2
         if(j.ge.k)goto 1
         jk=ikyk+j
         if(j.ge.l)goto 1
         jl=iky(l)+j
1        ajk=fock(jk)-gil*p(il)
         ail=fock(il)-gil*p(jk)
         aik=fock(ik)-gik*p(jl)
         fock(jl)=fock(jl)-gik*p(ik)
         fock(jk)=ajk
         fock(il)=ail
         fock(ik)=aik
c
6000  continue
c
      return
      end
      subroutine diismp(c,z,indix,ntotal,indc,iprn)
c
      implicit real*8 (a-h,o-z), integer(i-n)
c
      dimension c(ntotal),z(ntotal),indix(3,2)
c
c...  perform 3*3 diis convergence acceleration for ci and mp2/3
c...   c,z   :  work-space
c...   indix :  diskaddresses on ed7 for vector-storage
c...            (1/2/3,1) : c-vectors of iteration i/i-1/i-2
c...            (1/2/3,2) : residue-vectors of iteration i/i-1/i-2
c...   in DIRECT indix is probably index(31-36)
c...   indc    : address for c vector
c...   iprn    : print flag.
c...   on entrance and exit : z = current residuum
c...                          c = scratch-space
c...   // peter pulay, Fayetteville Arkansas 1989 \\
c
c
      integer nxblk, index, numci, mblkci, ispbig, ispsma
      common/disktl/nxblk,index(199),numci,mblkci,ispbig,ispsma
      dimension diism(3,3),lll(3),mmm(3)
      save ncyc
      data ncyc/0/
c
      ncyc = ncyc + 1
c
c...     swap residu and c addresse
c
         do 10 i=1,2
            ind = indix(3,i)
            indix(3,i) = indix(2,i)
            indix(2,i) = indix(1,i)
            indix(1,i) = ind
10       continue
c
c...     dump current c/z vectors (z-vector is in core)
c
            call wrt3(z,ntotal,indix(1,2),numci)
            call rdedx(c,ntotal,indc,numci)
            call wrt3(c,ntotal,indix(1,1),numci)
c
      if (ncyc.gt.2) then
c
c...    dot product of current residuum with  c's
c
c        call rdedx(c,ntotal,indix(1,1),8)  is already in core
c
         diism(1,1) = ddot(ntotal,c,1,z,1)
         call rdedx(c,ntotal,indix(2,1),numci)
         diism(2,1) = ddot(ntotal,c,1,z,1)
         call rdedx(c,ntotal,indix(3,1),numci)
         diism(3,1) = ddot(ntotal,c,1,z,1)
c
c...   dot product of previous residuum with c's
c
         call rdedx(z,ntotaL,indix(2,2),numci)
         diism(3,2) = ddot(ntotal,c,1,z,1)
         call rdedx(c,ntotal,indix(2,1),numci)
         diism(2,2) = ddot(ntotal,c,1,z,1)
         call rdedx(c,ntotaL,indix(1,1),numci)
         diism(1,2) = ddot(ntotal,c,1,z,1)
c
c...   dot product of oldest residuum with c's
c
         call rdedx(z,ntotal,indix(3,2),numci)
         diism(1,3) = ddot(ntotal,c,1,z,1)
         call rdedx(c,ntotal,indix(2,1),numci)
         diism(2,3) = ddot(ntotal,c,1,z,1)
         call rdedx(c,ntotal,indix(3,1),numci)
         diism(3,3) = ddot(ntotal,c,1,z,1)
c
c...  Now the matrix diism(3,3) contains <c(i)|z(j)> where i,j=1
c...  corresponds to the current value, 2 to the previous one,
c...  and 3 to the value two before. Then row 2 is subtracted from
c...  row 1, and row 3 from row 2. This is equivalent of orthogonalizing
c...  the interpolated residue to (C1-C2) and (C2-C3). The third
c...  row of the diism matrix is all 1's: the sum of the three
c...  coefficients is 1. This arrangement minimizes numerical errors
c
c...     according to Hestenes' theory, the interpolated Z must
c...     be orthogonal to delta(present-C) and delta(previous-C)
c...     a good summary of this is given by Wormer,Paldus IJQC ???
c..      see also : ????
c
        do 635 iii=1,3
          diism(1,iii) = diism(1,iii)-diism(2,iii)
          diism(2,iii) = diism(2,iii)-diism(3,iii)
          diism(3,iii) = 1.0d0
 635    continue
c          write(*,987) diism(1,1),diism(1,2),diism(1,3),diism(2,1),
c     1      diism(2,2),diism(2,3),diism(3,1),diism(3,2),diism(3,3)
c
c...    normalize the equations
c
        s1 = diism(1,1)
        if (abs(s1).lt.1.0d-15) s1=1.0d-15
        s2 = diism(2,2)
        if (abs(s2).lt.1.0d-15) s2=1.0d-15
        do 695 iii=1,3
           diism(1,iii)=diism(1,iii)/s1
           diism(2,iii)=diism(2,iii)/s2
 695    continue
c          write(*,987) diism(1,1),diism(1,2),diism(1,3),diism(2,1),
c     1      diism(2,2),diism(2,3),diism(3,1),diism(3,2),diism(3,3)
 987      format(' diis',5x, 3e13.6,/,10x,3e13.6,/,10x,3e13.6)
c
c...    invert the diism matrix. Its last column will give coefficients
c
        call osinv1(diism,3,determ,1.0d-14 ,lll,mmm)
        if(iprn.ge.2) then
           write(6,988) determ,(diism(kk,3),kk=1,3)
        endif
 988    format(11x,'DIIS det ',e11.4,' coeff',3F13.6)
c
c...    calculate extrapolated c as linear combination of the c's
c...    cex = c(1)*diism(1,3) + c(2)*diism(2,3) + c(3)*diism(3,3)
c
        call rdedx(c,ntotal,indix(3,1),numci)
        call dscal(ntotal,diism(3,3),c,1)
        call rdedx(z,ntotal,indix(2,1),numci)
        call daxpy(ntotal,diism(2,3),z,1,c,1)
        call rdedx(z,ntotal,indix(1,1),numci)
        call daxpy(ntotal,diism(1,3),z,1,c,1)
c
c...    replace c by it's extrapolated counterpart
c...    also dump it (for whoever)
c
        call wrt3(c,ntotal,indix(1,1),numci)
        call wrt3(c,ntotal,indc,numci)
c
c...    calculate extrapolated residue as linear combination of the z's
c...    zex = z(1)*diism(1,3) + z(2)*diism(2,3) + z(3)*diism(3,3)
c...    exactly like c's / but leave result in z (to return)
c
        call rdedx(z,ntotal,indix(3,2),numci)
        call dscal(ntotal,diism(3,3),z,1)
        call rdedx(c,ntotal,indix(2,2),numci)
        call daxpy(ntotal,diism(2,3),c,1,z,1)
        call rdedx(c,ntotal,indix(1,2),numci)
        call daxpy(ntotal,diism(1,3),c,1,z,1)
c
c...    replace z by it's extrapolated friend
c...     dumping not needed (return param.)
c
        call wrt3(z,ntotal,indix(1,2),numci)
c
c...   check - calculate the scalar products with the c's
c        call getcz(cr,cr(kmpz+1),indix(30))
c        s1=ddot(ntotal,cr(kmpc+1),1,cr(kmpz+1),1)
c        call getcz(cr,cr(kmpz+1),indix(31))
c        s2=ddot(ntotal,cr(kmpc+1),1,cr(kmpz+1),1)
c        write(*,*) 's1,s2', s1,s2
c
      endif
c
      return
c
      entry dinimp
c
      ncyc = 0
c
      return
      end
c************ MR-MP gen ***********************************************
      subroutine mpmr(cr,conf,energy)
c
      implicit real*8 (a-h,o-z), integer(i-n)
c
      external fget
      integer  conf
      dimension conf(5,*)
      dimension cr(*)
c
c...  control-routine for mr-mp2/3 calculation
c
c...  core usage has been reduced by using better  contraction/expansion
c...  routines are ccntmp2 and ccxpmp2; These require now one space/spin
c...  symmetry part of c/z in core; This may be reduced further,if need be.
c...  Various routines are prepared for this
c...  **note** mrmpn and contracted ci have not been changed
c...           they may be adapted similarly if the need arises
c...  JvL 1998
c
c     disk -layout
c     R, HR (CSF) : index(26),index(27)
c     contracted (excited state) basis :
c     (H-F)R      : indexp(1)
c     Fock-diag   : indexp(2)
c     x and Fx    : indexp(3),indexp(4)
c     diis storage : indexp(5-10)
c
c     iteration process :
c     C0 = -[HR]/(D-E0)
c     Cn = Cn-1 - [HR + (F-E0)Cn-1]/(D-E0)
c     this is equivalent to
c     Cn = -[HR + (F-D)Cn-1]/(D-E0)
c
c     the MP2 energy is calculated in Hylleraas form (more stable)
c        emp2 = ehyll = <C1|F-E0|C1> + 2<C1|H|C0>
c        the 'first order' mp2-energy emp21 = <C1|H|C0> printed evry it
c        as a consequence the last residue is not used in updating !
c
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
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
      common/mr_mp2l/ ivoc(62),nvoc,idoc(62),ndoc,kvcc(4),nvcc(4),
     *                ncmphn(10,18),ncmph1(10,18,8),ncmph2(10,18,8),
     *                ncmphd(10,18)
      common /mr_mp2con/ nsmpv,nsmp1,nsmp2,nsmp2d,ntotmp,
     *                   nmpv,nmp1,nmp2,nmp2d
      common/diisdc/ diiscr,ntdiis
c
c
      integer nref, nnst, nmin1, nmin2, mspin, nelec
      integer nstrt1, nlast1, nstrt2, nlast2, nmin12, nminv2, nminv1
      integer nlst21, nlastx
      common /ccntl/ nref,nnst,nmin1,nmin2,mspin,nelec,
     +               nstrt1,nlast1,nstrt2,nlast2,nmin12,nminv2,nminv1,
     +               nlst21,nlastx
c
c
      integer nval, ndoub, nsingl, ntripl, ntotal, ibig, jbig, kbig
      integer nvd, nvds, nsitri, lbig, mbig, nbig, nref0
      common /cntrl/ nval,ndoub,nsingl,ntripl,ntotal,ibig,jbig,kbig,
     +               nvd,nvds,nsitri,lbig,mbig,nbig,nref0
c
c
      integer nint, next, ninnex, nwm1, nw, nwp1, nintwo
      integer nint21, nintp1
      common /orb/ nint,next,ninnex,nwm1,nw,nwp1,nintwo,
     +             nint21,nintp1
c
c
      integer nxblk, index, numci, mblkci, ispbig, ispsma
      common/disktl/nxblk,index(199),numci,mblkci,ispbig,ispsma
      common/stoctl/khbill,khbz,khbc,khbb,khbzd,khbzs,khbzt,khbcd,
     *khbcs,khbct,kh,khb1,khb2,khb3,khb4,khb5,khb6,khb7,khb8
     *,khb44,khb55,khb66,khb88
c
      integer lcall,  mrstor
      logical modev, moded, modest, instor, loijka, msort, loijab
      logical loiajb, lvar
      common/diactl/modev,moded,modest,instor,
     +              loijka,msort,loijab,loiajb,lcall,mrstor,lvar
c
      logical malt
      real*8 potnuc,totnuc,e11,eigen,h11,z11,h22,s22,h12,s12,
     *     tester,eshift,thresh,aaa,bbb,scream,gprin,z12,z22,
     *     vdstpr,qester
      integer maxcyc,moddia,mprin,icyc
      common/erg/potnuc,totnuc,e11,eigen(nd200),h11
     *,z11,h22,s22,h12,s12,tester,eshift,thresh,aaa,bbb,scream
     *,gprin,z12,z22,vdstpr(3),qester,maxcyc,moddia,mprin,malt,icyc
c
      integer nirr,mult,ifrep,nirr8,jfrep
      common /table/ nirr,mult(8,8),ifrep,nirr8,jfrep
c
c
      integer nconst, ncond, maxpa, maxpat, maxsin, maxepa
      integer mxdvpa, mxddpa, nspips, nspipt, mxsspa, mxttpa
      integer mxstpa, nsitot, ntrtot
      integer mx12pa, mxvvpa, mxv2pa, maxpsq, max12
      common /auxh/ nconst(8),ncond(8),maxpa(4),maxpat,maxsin(3),
     +              maxepa(4),mxdvpa,mxddpa,nspips(8),nspipt(8),
     +              mxsspa,mxttpa,mxstpa,nsitot(8),ntrtot(8),
     +              mx12pa,mxvvpa,mxv2pa,maxpsq,max12
c
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
c
      integer nuse, ngot
      common /corctl/ nuse, ngot
c
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
c
      dimension epair_mp(2,9),indexp(10)
      logical diisc,trip
c
      parameter (nmodel=7)
      character*39 descr(nmodel)
      data descr /'      ** no projection applied **     ',
     *            ' |ref>F<ref|+|sin>F<sin|+|doub>F<doub|',
     *            '     |ref>F<ref| + |doub>F<doub|      ',
     *            '  |ref>F<ref| + |sin+doub>F<sin+doub| ',
     *            ' |ref>F<ref|+|sin>F<sin|+|rest>F<rest|',
     *            '     |ref>F<ref| + |rest>F<rest|      ',
     *            ' |ref>F<ref| + |doub+more>F<doub+more|'/
c
      diisc = .false.
      pplev = plevel(1)
c
      call cpuwal(top,bot)
      write(6,444)top,bot
444   format(//1x,80('-')//
     *' ****MR-MP 1st order v1.3 called at',f12.1,' wall',f12.1,' secs')
       maxtot=max(nval,ndoub)
       do i=1,nirr
          maxtot = max(maxtot,nspips(i)*norbsi(i),nspipt(i)*norbtr(I))
       end do
chvd
c...  save config and common blocks for later use in MP3
c
      call cistor
      call gciblk(mpcblk)
      call wrt3(cr,5*(nnst+nmin1+nmin2),mpcblk,numci)
c
      if (iprimp.ge.2000) then
         if (khbill+ntotal.gt.ngot) then
            write(6,608)
618         format(' - HC0 (CSF) - not printed due to lack of core')
         else
            call getcz(conf,cr(khbill+1),index(27))
            write(6,603) (cr(khbill+i),i=1,ntotal)
603         format(' - HC0 (CSF) -',/,(2x,10e12.5))
         end if
      end if
chvd
c
c...   switch to new contracted addressing (khbill changes)
c
      call basmp(conf)
c
      nbl = (ntotmp-1)/511+1
      indexp(1) = index(28)
      do 5 i=2,10
5     indexp(i) = indexp(i-1) + nbl
c
c...   contract
c
      kmpz = khbill
      kmpc = kmpz + ntotmp
      kcc  = kmpz
      kdc  = kcc  + nsmpv+nsmp1+nsmp2
      kconf = kmpc
      kc   = kconf + nlast2*5
      kvc  = kc   + maxtot
      kscr = kvc  + max(nvcc(1),nvcc(2),nvcc(3),nvcc(4))
      nuse = kscr + max(nspips(1)*next,maxpat**2)
      call corfai('contract')
c
      if (iprimp.ge.2100) then
c       print transformation matrix.
        call vvpri(cr,kc,kcc,kdc,kvc,kscr,ntotal,ntotmp,1.0d-8)
      endif
c
      call ccntmp2(cr(kc+1),cr(kcc+1),cr(kdc+1),cr(kvc+1),cr(kscr+1),
     1            cr(kconf+1),index(27))
c
      if (iprimp.ge.2000) write(6,602) (cr(kcc+i),i=1,ntotmp)
602   format(' - HC0 (excited states) -',/,(2x,10e12.5))
c602   format(' - HC0 (excited states) -',/,(2x,5e22.15))
c
c...  check if first element is indeed e0ref  (**debug**)
c..
      tttt = cr(kcc+1) - e0ref
      if (tttt.ne.0.0d0) write(6,606) tttt
606   format(' #@#@#@ devation from e0ref :',e14.5,' #@#@#')
      if (dabs(tttt).gt.1.0d-9) call caserr('mp-contraction') 
      cr(kcc+1) = 0.0d0
      call wrt3(cr(kcc+1),ntotmp,indexp(1),numci)
c
c...  calculate diagonal of fock-matrix
c
      call vclr(cr(kmpz+1),1,ntotmp)
c
      if(odebug(40)) then
c     debug print
       write(6,*) 'before mpdiag: conf'
       call confprt(conf)
      endif
c
      call mpdiag(cr,cr(kmpz+1),cr(kscr+1))
c
      if(odebug(40)) then
c      debug print
       write(6,*) 'after mpdiag: conf'
       call confprt(conf)
      endif
c
      if (iprimp.ge.2000) write(6,605) (cr(kmpz+i),i=1,ntotmp)
c605   format(' - diagonal (F-E0) -',/,(2x,10e12.5))
605   format(' - diagonal (F-E0) -',/,(2x,5e22.15))
c
c...   invert sign of diagonal (F-E0) => (E0-F)  / dump
c...   put diagonal for reference arbitrary to 1.0d0
c
      call dscal(ntotmp,-1.0d0,cr(kmpz+1),1)
      cr(kmpz+1) = 1.0d0
      call wrt3(cr(kmpz+1),ntotmp,indexp(2),numci)
c
c...   calculate first psi1(1) =  HR/D
c
      call rdedx(cr(kmpc+1),ntotmp,indexp(1),numci)
c
c...    level-shift
c
      if (pplev.ne.0.0d0)
     *   call vsadd(cr(kmpz+1),1,-pplev,cr(kmpz+1),1,ntotmp)
c
      call vvdv(ntotmp,cr(kmpc+1),cr(kmpc+1),cr(kmpz+1))
c
      call wrt3(cr(kmpc+1),ntotmp,indexp(3),numci)
c
c...  calculate start-mp2 energy  e1h0
c
      call rdedx(cr(kmpz+1),ntotmp,indexp(1),numci)
      e1h0 = ddot(ntotmp,cr(kmpc+1),1,cr(kmpz+1),1)
c
      if (iprimp.ge.100) then
         if (iprimp.ge.2000) write(6,604) (cr(kmpc+i),i=1,ntotmp)
c604      format(' - current psi(1) -',/,(2x,10e12.5))
604      format(' - current psi(1) -',/,(2x,5e22.15))
      end if
c
c...  calculate (F-E0)*C
c
      ncyc = 0
c
c...   start of fock-matrix iteration
c
10    continue
c
      call mpgen(conf,cr)
c
c..   calculate <C1|(F-E0)|C1> for hylleraas mp2
c
      call wrt3(cr(kmpz+1),ntotmp,indexp(4),numci)
c
      call rdedx(cr(kmpc+1),ntotmp,indexp(3),numci)
      e1f1 = ddot(ntotmp,cr(kmpc+1),1,cr(kmpz+1),1)
c
c...   ehyll = <C1|F-E0|C1> + 2<C1|H|C0>
c
      ehyll = e1f1 + 2.0d0*e1h0
      if(iprimp.ge.1000) then
         write(6,668) ehyll,e1h0
      endif
668   format(10x,'  E2-Hylleraas ',f16.13,'   E2(1) ',f16.13)
c
c..   build r = HR+(F-E0)C   in z
c
      call rdedx(cr(kmpc+1),ntotmp,indexp(1),numci)
c
      if (iprimp.ge.2000) write(6,607) (cr(kmpz+i),i=1,ntotmp)
c607   format(' - (F-E0) psi(1) -',/,(2x,10e12.5))
607   format(' - (F-E0) psi(1) -',/,(2x,5e22.15))
c
      call vadd(cr(kmpz+1),1,cr(kmpc+1),1,cr(kmpz+1),1,ntotmp)
c
c...  residue calculated (in 'z')
c...  check convergence (on residue)
c
      if (iprimp.ge.2000) write(6,608) (cr(kmpz+i),i=1,ntotmp)
608   format(' - residual-vector -',/,(2x,10e12.5))
      residu = 0.0d0
      do i=1,ntotmp
         residu = max(residu,dabs(cr(kmpz+i)))
      end do
c
c...    if converged skip last updating step
c
      if (residu.lt.critmp) go to 15
c
c...  Peter Pulay's famous DIIS
c
      if (residu.lt.diiscr.and.ncyc.ge.ntdiis) diisc = .true.
      if (diisc)
     * call diismp(cr(kmpc+1),cr(kmpz+1),indexp(5),ntotmp,indexp(3),
     *             iprimp)
c
      if (iprimp.ge.2000) write(6,609) (cr(kmpz+i),i=1,ntotmp)
609   format(' - residual-vector (after DIIS) -',/,(2x,10e12.5))
c
c...  check extrapolated residue
c
      residd = 0.0d0
      do i=1,ntotmp
         residd = max(residd,dabs(cr(kmpz+i)))
      end do
c
c...    r/D
c
      call rdedx(cr(kmpc+1),ntotmp,indexp(2),numci)
c
c...    level-shift
c
      if (pplev.ne.0.0d0)
     *   call vsadd(cr(kmpc+1),1,-pplev,cr(kmpc+1),1,ntotmp)
c
      call vvdv(ntotmp,cr(kmpz+1),cr(kmpz+1),cr(kmpc+1))
c
c...  check largest change in first-order wavefunction (not used)
c
      change = 0.0d0
      do i=1,ntotmp
         change = max(change,dabs(cr(kmpz+i)))
      end do
c
      call cpuwal(top,bot)
      ncyc = ncyc + 1
      if(iprimp.ge.1000) then
         write(6,666) ncyc,top,bot,change,residu,residd
      else
         write(6,670) ncyc,top,bot,ehyll,residu
      endif
666   format(' MP-cycle',i4,' at',f11.1,' wall',f11.1,' secs',
     *       ' change',e11.4,' residue',e11.4,' diis residue',e11.4)
670   format(' MP-cycle',i4,' at',f11.1,' wall',f11.1,' secs',
     *       '  E2-Hylleraas ',f17.13,' residue',e10.3)
      if (ncyc.ge.levelp) pplev = plevel(2)
c
c...   update c :  x = x + r/d
c
      call rdedx(cr(kmpc+1),ntotmp,indexp(3),numci)
      call vadd(cr(kmpz+1),1,cr(kmpc+1),1,cr(kmpc+1),1,ntotmp)
c
c     write(6,603) cr(kmpc+1)
      if (iprimp.ge.2000) write(6,604) (cr(kmpc+i),i=1,ntotmp)
c
      call wrt3(cr(kmpc+1),ntotmp,indexp(3),numci)
c
c...  calculate mp2 energy  e1h0
c
      call rdedx(cr(kmpz+1),ntotmp,indexp(1),numci)
      e1h0 = ddot(ntotmp,cr(kmpc+1),1,cr(kmpz+1),1)
c
      if (ncyc.lt.maxcmp) go to 10
c
c...   no convergence // print out best results // skip fancy print
c
c     write(6,603) cr(kmpc+1)
      cc00 = ddot(ntotmp,cr(kmpc+1),1,cr(kmpc+1),1)
      energy = 0.0d0
      write(6,667) e1h0,ehyll,cc00
667   format(' *&*&*&*&*&*&*&* no convergence *&*&*&*&*&*&*&*',/,
     *       ' &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&',/,
     *       ' &&          Psi(1) is updated               &&',/,
     *       ' &&     new E2(1) ',f16.13,       '          &&',/,
     *       ' &&     old E2(h) ',f16.13,       '          &&',/,
     *       ' &&      <C1|C1>  ',f8.5,'                   &&',/,
     *       ' &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'/)
      go to 100
c
c...  we converged the first-order wavefunction
c
15    call cpuwal(top,bot)
      ncyc = ncyc + 1
      write(6,610) ncyc,top,bot
610   format(' MP-convergence at cycle',i4,' at',f11.1,' wall',f11.1,
     *       ' secs ')
c
c...   get norm of  c-vector
c
      call rdedx(cr(kmpc+1),ntotmp,indexp(3),numci)
      cc00 = ddot(ntotmp,cr(kmpc+1),1,cr(kmpc+1),1)
c
c...  get overlap with reference function
c
      overla = cr(kmpc+1)
c
c...  calculate pair-energies // in mp-mode
c...  first order part :
c...   <ref|H|psi1> =   <c | z28)
c
      call rdedx(cr(kmpz+1),ntotmp,indexp(1),numci)
c
      call vclr(epair_mp,1,2*9)
c
      call mpsum(conf,cr(kmpc+1),cr(kmpz+1),epair_mp)
c
c...  now multiply by 2.0 and get <C1|F-E0|C1> part of pair energies
c...   Emp2 (hylleraas) = <C1|F-E0|C1> + 2<C1|H|C0>
c
      call dscal(2*9,2.0d0,epair_mp,1)
c
      call rdedx(cr(kmpz+1),ntotmp,indexp(4),numci)
c
      call mpsum(conf,cr(kmpc+1),cr(kmpz+1),epair_mp)
c
c...  now epair(i,j) contains following pair-energies
c...  i = 1 singles / i = 2 doubles
c...  j= 1 ..3 : # holes + 1 in vacuum
c...  j= 4 ..6 : # holes + 4 in n-1
c...  j= 7 ..9 : # holes + 7 in n-2
c
      e1vac = epair_mp(1,1)+epair_mp(1,2)+epair_mp(1,3)
      e2vac = epair_mp(2,1)+epair_mp(2,2)+epair_mp(2,3)
      e1n1  = epair_mp(1,4)+epair_mp(1,5)+epair_mp(1,6)
      e2n1  = epair_mp(2,4)+epair_mp(2,5)+epair_mp(2,6)
      e1n2  = epair_mp(1,7)+epair_mp(1,8)+epair_mp(1,9)
      e2n2  = epair_mp(2,7)+epair_mp(2,8)+epair_mp(2,9)
      ecorr = e1vac+e2vac+e1n1+e2n1+e1n2+e2n2
c
      energy = e0ref+ecorr
      write(6,601) mpmod,descr(abs(mpmod)+1),nref,mspin,
     *             e0ref,ecorr,e0ref+ecorr,cc00,residu,
     *             e1vac,e2vac,e1n1,e2n1,e1n2,e2n2
      if (mpmod.lt.0) write(6,*)
     *         ' ** H0 block diagonal in excitation classes         **'
601   format(//' $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$',
     *        /' $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$',
     *        /' $$             Multi-reference MP2                 $$',
     *        /' $$                 model',i6,'                     $$',
     *        /' $$   ',              a39                  ,'       $$',
     *        /' $$          Nref ',i6,'   Spin ',i6,'              $$',
     *        /' $$                                                 $$',
     *        /' $$      E-reference   ',f18.10,        '           $$',
     *        /' $$      E-correlation ',f18.10,        '           $$',
     *        /' $$      --- total --- ',f18.10,        '           $$',
     *        /' $$                                                 $$',
     *        /' $$      <psi1|psi1>   ',e18.5,         '           $$',
     *        /' $$        residue     ',e18.5,         '           $$',
     *        /' $$                                                 $$',
     *        /' $$             singles        doubles              $$',
     *        /' $$  vacuum ',e13.6,3x,e13.6,10x,                 ' $$',
     *        /' $$  n-1    ',e13.6,3x,e13.6,10x,                 ' $$',
     *        /' $$  n-2    ',e13.6,3x,e13.6,10x,                 ' $$',
     *        /' $$                                                 $$',
     *        /' $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$',
     *        /' $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$',
     *        /)
c
c...   if (overlap with psi(r) really big) print warning // error??
c
      if (dabs(overla).gt.critmp) write(6,611) overla
611   format(' #@#@#@#@#@# <psi(1)|psi(r)> =',e11.3,' #@#@#@#@#@#')
100   continue
c
c...  expand the excitated-state vector to configuration-state vector
c
      call rdedx(cr(kcc+1),ntotmp,indexp(3),numci)
      cr(kcc+1) = 1.0d0
      if (iprimp.gt.1) then
        write(6,621)
621     format('MP2 vector')
        call ccvpri(cr,cr(kcc+1),1.0d-4)
      endif
c
      call ccxpmp2(cr(kc+1),cr(kcc+1),cr(kdc+1),cr(kvc+1),cr(kscr+1),
     1             cr(kconf+1),index(26))
c
c...  restore original information in the common blocks to enable
c...  the use of genz.
c
      call ciretr
c
c...  retrieve config information to do the mp3
c
      call gciblk(mpcblk)
      call rdedx(conf,5*(nnst+nmin1+nmin2),mpcblk,numci)
c
c...  get  current vector for instore mode
c
      if (instor) call getcz(conf,cr(khbc+1),index(26))
cmp3
c...  calculate  Z = H*C
c
      if (ncyc.lt.maxcmp.and.mp.ge.3) then
         call corfai('mp3')
         call genz(cr,conf,index(26),h11)
         if(instor) call dumpcz(conf,cr(khbz+1),index(27))
c
c...     E(3) = <C0+C1|H-E1|C0+C1> + E1 = <C0+C1|H|C0+C1> - E1*<C1|C1>
c...                                           - 2*E1*<C0|C1>
c
         emp3 = h11 - e0ref*(cc00+2.0d0*overla)
         h11 = h11/(1.0d0+cc00+2.0d0*overla)
c
c...     output MP3 results
c
         energy = emp3
         write(6,612) emp3-e0ref,emp3
      endif
612   format( ' $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$',
     *       /' $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$',
     *       /' $$             Multi-reference MP3                 $$',
     *       /' $$              Ecorr',f18.10,        '            $$',
     *       /' $$              E    ',f18.10,        '            $$',
     *       /' $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$',
     *       /' $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$',
     *       //)
cmp3
c
c..   if (maxcyc.ne.0) there is more to come
c..   so then normalise c and z-vectors
c..   normalise them anyway print, a warning though what 
c..   the natural orbitals mean
c
      cc00 = dsqrt(1.0d0+cc00+2.0d0*overla)
      cc00 = 1.0d0/cc00
      write (6,613) cc00
613   format(/' $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$',
     *       /' $$      Normalising c/z with ',e12.5,   '          $$',
     *       /' $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'/)
c
      if (maxcyc.le.0) write(6,614)
614   format(/' $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$',
     *       /' $$    All the following (natorb/vector-printer)    $$',
     *       /' $$     is based on this renormalised CI-vector     $$',
     *       /' $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'/)

      ibl = index(26)
      trip = .false.
c
1000  call search(ibl,numci)
      n = 0
1001  call fget(cr(khbc+1),m,numci)
      call dscal(m,cc00,cr(khbc+1),1)
      call wrt3(cr(khbc+1),m,ibl,numci)
      ibl = ibl + 1
      n = n + m
      if (n.ne.ntotal) go to 1001
c
      if (.not.trip) then
         ibl = index(27)
         trip = .true.
         go to 1000
      end if
c
      return
      end
      subroutine mpmr2(cr,conf)
c
      implicit real*8 (a-h,o-z), integer(i-n)
c
      integer  conf
      dimension cr(*),conf(5,1)
c
c...  control-routine for mr-mpn calculation
c
c     disk -layout
c     R, HR (CSF) : index(26),index(27)
c     contracted (excited state) basis :
c     (H-F)R      : indexp(1)
c     Fock-diag   : indexp(2)
c     x and Fx    : indexp(3),indexp(4)
c     diis storage : indexp(5-10)
c     storage for contracted c,z vectors : indexp(11-110)
c     indexp(ic(i)) : disk-address for (i-1)-th order correction to the
c                     wavefunction c(i-1).
c     indexp(iz(i)) : disk-address for z(i-1) = H*c(i-1)
c
c     iteration process :
c     C0 = -[HR]/(D-E0)
c     Cn = Cn-1 - [HR + (F-E0)Cn-1]/(D-E0)
c     this is equivalent to
c     Cn = -[HR + (F-D)Cn-1]/(D-E0)
c
c     the MP2 energy is calculated in Hylleraas form (more stable)
c        emp2 = ehyll = <C1|F-E0|C1> + 2<C1|H|C0>
c        the 'first order' mp2-energy emp21 = <C1|H|C0> printed evry it
c        as a consequence the last residue is not used in updating !
c
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
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
      common/mr_mp2l/ ivoc(62),nvoc,idoc(62),ndoc,kvcc(4),nvcc(4),
     *                ncmphn(10,18),ncmph1(10,18,8),ncmph2(10,18,8),
     *                ncmphd(10,18)
      common /mr_mp2con/ nsmpv,nsmp1,nsmp2,nsmp2d,ntotmp,
     *                   nmpv,nmp1,nmp2,nmp2d
      common/diisdc/ diiscr,ntdiis
c
c
      integer nref, nnst, nmin1, nmin2, mspin, nelec
      integer nstrt1, nlast1, nstrt2, nlast2, nmin12, nminv2, nminv1
      integer nlst21, nlastx
      common /ccntl/ nref,nnst,nmin1,nmin2,mspin,nelec,
     +               nstrt1,nlast1,nstrt2,nlast2,nmin12,nminv2,nminv1,
     +               nlst21,nlastx
c
c
      integer nval, ndoub, nsingl, ntripl, ntotal, ibig, jbig, kbig
      integer nvd, nvds, nsitri, lbig, mbig, nbig, nref0
      common /cntrl/ nval,ndoub,nsingl,ntripl,ntotal,ibig,jbig,kbig,
     +               nvd,nvds,nsitri,lbig,mbig,nbig,nref0
c
c
      integer nint, next, ninnex, nwm1, nw, nwp1, nintwo
      integer nint21, nintp1
      common /orb/ nint,next,ninnex,nwm1,nw,nwp1,nintwo,
     +             nint21,nintp1
c
c
c
      integer nxblk, index, numci, mblkci, ispbig, ispsma
      common/disktl/nxblk,index(199),numci,mblkci,ispbig,ispsma
      common/stoctl/khbill,khbz,khbc,khbb,khbzd,khbzs,khbzt,khbcd,
     *khbcs,khbct,kh,khb1,khb2,khb3,khb4,khb5,khb6,khb7,khb8
     *,khb44,khb55,khb66,khb88
c
      integer lcall,  mrstor
      logical modev, moded, modest, instor, loijka, msort, loijab
      logical loiajb, lvar
      common/diactl/modev,moded,modest,instor,
     +              loijka,msort,loijab,loiajb,lcall,mrstor,lvar
c
      logical malt
      real*8 potnuc,totnuc,e11,eigen,h11,z11,h22,s22,h12,s12,
     *     tester,eshift,thresh,aaa,bbb,scream,gprin,z12,z22,
     *     vdstpr,qester
      integer maxcyc,moddia,mprin,icyc
      common/erg/potnuc,totnuc,e11,eigen(nd200),h11
     *,z11,h22,s22,h12,s12,tester,eshift,thresh,aaa,bbb,scream
     *,gprin,z12,z22,vdstpr(3),qester,maxcyc,moddia,mprin,malt,icyc
c
      integer nirr,mult,ifrep,nirr8,jfrep
      common /table/ nirr,mult(8,8),ifrep,nirr8,jfrep
c
c
      integer nconst, ncond, maxpa, maxpat, maxsin, maxepa
      integer mxdvpa, mxddpa, nspips, nspipt, mxsspa, mxttpa
      integer mxstpa, nsitot, ntrtot
      integer mx12pa, mxvvpa, mxv2pa, maxpsq, max12
      common /auxh/ nconst(8),ncond(8),maxpa(4),maxpat,maxsin(3),
     +              maxepa(4),mxdvpa,mxddpa,nspips(8),nspipt(8),
     +              mxsspa,mxttpa,mxstpa,nsitot(8),ntrtot(8),
     +              mx12pa,mxvvpa,mxv2pa,maxpsq,max12
c
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
c
      integer nuse, ngot
      common /corctl/ nuse, ngot
c
c
      integer maxord
      parameter(maxord = 50)
c...  maxord = maximum order in perturbation theory
      dimension indexp(10+2*maxord), ic(maxord), iz(maxord)
      dimension epert(maxord)
c
      logical diisc
c
      parameter (nmodel=7)
      character*39 descr(nmodel)
      data descr /'      ** no projection applied **     ',
     *            ' |ref>F<ref|+|sin>F<sin|+|doub>F<doub|',
     *            '     |ref>F<ref| + |doub>F<doub|      ',
     *            '  |ref>F<ref| + |sin+doub>F<sin+doub| ',
     *            ' |ref>F<ref|+|sin>F<sin|+|rest>F<rest|',
     *            '     |ref>F<ref| + |rest>F<rest|      ',
     *            ' |ref>F<ref| + |doub+more>F<doub+more|'/
c
      character*8 zblank
      data zblank / ' '/
c
      diisc = .false.
      pplev = plevel(1)
c
      call cpuwal(top,bot)
      write(6,444)top,bot
444   format(//1x,80('-')//
     *' ****MR-MP 1st order v1.3 called at',f12.1,' wall',f12.1,' secs')
chvd
c...  save config and common blocks for later use in MP3
c
      call cistor
      call gciblk(mpcblk)
      call wrt3(cr,5*(nnst+nmin1+nmin2),mpcblk,numci)
chvd
c...  contract the Z-vector
c
      nuse = khbill + ntotal
      call corfai('contract')
c...  get z-vector
      kmpz = khbill
      if (instor) then
         if (khbz.ne.khbill) call caserr('khbz misplaced in mpgen')
      else
         call getcz(cr,cr(kmpz+1),index(27))
      end if
c...  dump back to disk  at index(27) straight
      call wrt3(cr(kmpz+1),ntotal,index(27),numci)
c
c...  switch to new contracted addressing (khbill changes)
c
      call basmp(conf)
c
c...  build disc addresses
c
      nbl = (ntotmp-1)/511+1
      indexp(1) = index(28)
      do 5 i=2,10+2*maxord
5     indexp(i) = indexp(i-1) + nbl
      ic(1) = 10 + 1
      iz(1) = 10 + 2
      do 6 i=2,maxord
        ic(i) = ic(i-1) + 2
        iz(i) = iz(i-1) + 2
6     continue
c
c...  contract
c
      kmpz = khbill
      kmpc = kmpz + ntotmp
      kmpb = kmpc + ntotmp
      kcc  = kmpz
      kdc  = kcc  + nsmpv+nsmp1+nsmp2
      kc   = kmpz + ntotmp
      kvc  = kc   + ntotal
      kscr = kvc  + max(nvcc(1),nvcc(2),nvcc(3),nvcc(4))
      nuse = kscr + max(nspips(1)*next,norbtr(1)*nspips(1))
      call corfai('contract')
chvd
      write(6,446)
446   format(//3x,89('-')/'   order  time(wall)  time(secs)',
     *10x,'energy(au) perturbation(au)'
     */3x,89('-'))
chvd
c...  save config and common blocks for later use in contracted
c...  stuff
c
      call mpstor
      call gmpblk(mpcblk)
      call wrt3(cr,5*(nmpv+nmp1+nmp2+nmp2d),mpcblk,numci)
chvd
      if (iprimp.ge.2100) then
c       print transformation matrix.
        call vvpri(cr,kc,kcc,kdc,kvc,kscr,ntotal,ntotmp,1.0d-8)
      endif
chvd
c...  retrieve, contract, and store zero-th order c-vector
c
      call rdedx(cr(kc+1),ntotal,index(26),numci)
      if (iprimp.ge.2000) write(6,616) (cr(kc+i),i=1,ntotal)
616   format(' - C0 (CSF) -',/,(2x,10e12.5))
      call ccntmp(cr(kc+1),cr(kcc+1),cr(kdc+1),cr(kvc+1),cr(kscr+1))
      if (iprimp.ge.2000) write(6,614) (cr(kcc+i),i=1,ntotmp)
614   format(' - C0 (excited states) -',/,(2x,5e22.15))
      call wrt3(cr(kcc+1),ntotmp,indexp(ic(1)),numci)
c
c...  retrieve, contract, and store zero-th order z-vector
c
      call rdedx(cr(kc+1),ntotal,index(27),numci)
      if (iprimp.ge.2000) write(6,603) (cr(kc+i),i=1,ntotal)
603   format(' - HC0 (CSF) -',/,(2x,10e12.5))
      call ccntmp(cr(kc+1),cr(kcc+1),cr(kdc+1),cr(kvc+1),cr(kscr+1))
      if (iprimp.ge.2000) write(6,602) (cr(kcc+i),i=1,ntotmp)
c602   format(' - HC0 (excited states) -',/,(2x,10e12.5))
602   format(' - HC0 (excited states) -',/,(2x,5e22.15))
      call wrt3(cr(kcc+1),ntotmp,indexp(iz(1)),numci)
c
chvd
c...  check if first element is indeed e0ref  (**debug**)
c..
      tttt = cr(kcc+1) - e0ref
      if (tttt.ne.0.0d0) write(6,606) tttt
606   format(' #@#@#@ devation from e0ref :',e14.5,' #@#@#')
c
c...  calculate diagonal of fock-matrix
c
      call vclr(cr(kmpz+1),1,ntotmp)
      call mpdiag(cr,cr(kmpz+1),cr(kscr+1))
c
      if (iprimp.ge.2000) write(6,605) (cr(kmpz+i),i=1,ntotmp)
c605   format(' - diagonal (F-E0) -',/,(2x,10e12.5))
605   format(' - diagonal (F-E0) -',/,(2x,5e22.15))
c
c...   invert sign of diagonal (F-E0) => (E0-F)  / dump
c...   put diagonal for reference arbitrary to 1.0d0
c
      call dscal(ntotmp,-1.0d0,cr(kmpz+1),1)
      cr(kmpz+1) = 1.0d0
      call wrt3(cr(kmpz+1),ntotmp,indexp(2),numci)
c
chvd  SETUP START
c
c...  calculate E(0) = c(0)H0c(0)
c...  calculate E(1) = c(0)Hc(0)-E(0)
c
      call rdedx(cr(kmpc+1),ntotmp,indexp(ic(1)),numci)
      epert(1) = e0mp
      call rdedx(cr(kmpz+1),ntotmp,indexp(iz(1)),numci)
      epert(2) = ddot(ntotmp,cr(kmpc+1),1,cr(kmpz+1),1) - epert(1)
c
      call cpuwal(top,bot)
      write(6,601)0,top,bot,epert(1),0.0d0
c
c     set starting order
c
      iord = 2
c
chvd  BEGIN MPN
c
30    continue
c
c...  calculate E(i) = c(0)Hc(i-1)
c
      if(iord.gt.2) then
        call rdedx(cr(kmpc+1),ntotmp,indexp(ic(1)),numci)
        call rdedx(cr(kmpz+1),ntotmp,indexp(iz(iord-1)),numci)
        epert(iord) = ddot(ntotmp,cr(kmpc+1),1,cr(kmpz+1),1)
      endif
c
      eref  = epert(1) + epert(2)
      ecorr = 0.0d0
      do 250 i = 3, iord
        ecorr = ecorr + epert(i)
250   continue
c     write(6,601) mp,mpmod,descr(mpmod+1),nref,mspin,
c    *             e0mp,ecorr,e0mp+ecorr
      write(6,601) iord-1,top,bot,eref+ecorr,epert(iord)
c
c...  calculate z = (F - H)c(i-1)
c
      call rdedx(cr(kmpc+1),ntotmp,indexp(ic(iord-1)),numci)
      call mpgen(cr,cr)
      call daxpy(ntotmp,e0mp,cr(kmpc+1),1,cr(kmpz+1),1)
      call rdedx(cr(kmpc+1),ntotmp,indexp(iz(iord-1)),numci)
      call daxpy(ntotmp,-1.0d0,cr(kmpc+1),1,cr(kmpz+1),1)
c
c...  calculate z = z + (E(j)*c(i-j),j=1,i)
c
      do 20 i = 2, iord
        call rdedx(cr(kmpc+1),ntotmp,indexp(ic(iord-i+1)),numci)
        call daxpy(ntotmp,epert(i),cr(kmpc+1),1,cr(kmpz+1),1)
20    continue
      call dscal(ntotmp,-1.0d0,cr(kmpz+1),1)
c
c...  Now only to solve  (F-E0)c(i) = z
c
      call wrt3(cr(kmpz+1),ntotmp,indexp(1),numci)
      rnorm = dsqrt(ddot(ntotmp,cr(kmpz+1),1,cr(kmpz+1),1))
c
chvd  BEGIN LINEAR SYSTEM SOLVER
c
c...   calculate first psi1(1) =  HR/D
c
      call rdedx(cr(kmpc+1),ntotmp,indexp(1),numci)
      call rdedx(cr(kmpz+1),ntotmp,indexp(2),numci)
c
c...    level-shift
c
      if (pplev.ne.0.0d0)
     *   call vsadd(cr(kmpz+1),1,-pplev,cr(kmpz+1),1,ntotmp)
c
      call vvdv(ntotmp,cr(kmpc+1),cr(kmpc+1),cr(kmpz+1))
c
      call wrt3(cr(kmpc+1),ntotmp,indexp(3),numci)
c
c...  calculate start-mp2 energy  e1h0
c
      call rdedx(cr(kmpz+1),ntotmp,indexp(1),numci)
      e1h0 = ddot(ntotmp,cr(kmpc+1),1,cr(kmpz+1),1)
c
      if (iprimp.ge.100) then
         if (iprimp.ge.2000) write(6,604) (cr(kmpc+i),i=1,ntotmp)
c604      format(' - current psi(1) -',/,(2x,10e12.5))
604      format(' - current psi(1) -',/,(2x,5e22.15))
      end if
c
c...  calculate (F-E0)*C
c
      ncyc = 0
      call dinimp
c
c...  start of fock-matrix iteration
c
10    continue
c
      call mpgen(conf,cr)
c
c..   calculate <C1|(F-E0)|C1> for hylleraas mp2
c
      call wrt3(cr(kmpz+1),ntotmp,indexp(4),numci)
c
      call rdedx(cr(kmpc+1),ntotmp,indexp(3),numci)
      e1f1 = ddot(ntotmp,cr(kmpc+1),1,cr(kmpz+1),1)
c
c...  ehyll = <C1|F-E0|C1> + 2<C1|H|C0>
c
      ehyll = e1f1 + 2.0d0*e1h0
      write(6,668) ehyll,e1h0
668   format(10x,'  E2-Hylleraas ',f16.13,'   E2(1) ',f16.13)
c
c..   build r = HR+(F-E0)C   in z
c
      call rdedx(cr(kmpc+1),ntotmp,indexp(1),numci)
c
      if (iprimp.ge.2000) write(6,607) (cr(kmpz+i),i=1,ntotmp)
c607   format(' - (F-E0) psi(1) -',/,(2x,10e12.5))
607   format(' - (F-E0) psi(1) -',/,(2x,5e22.15))
c
      call vadd(cr(kmpz+1),1,cr(kmpc+1),1,cr(kmpz+1),1,ntotmp)
c
c...  residue calculated (in 'z')
c...  check convergence (on residue)
c
      if (iprimp.ge.2000) write(6,608) (cr(kmpz+i),i=1,ntotmp)
608   format(' - residual-vector -',/,(2x,10e12.5))
      residu = 0.0d0
      do i=1,ntotmp
         residu = max(residu,dabs(cr(kmpz+i)))
      end do
c
c...  if converged skip last updating step
c
      if(critmp.le.0.0d0) then
        if (residu.lt.min(dabs(critmp),dabs(critmp)*rnorm)) go to 15
      else
        if (residu.lt.critmp) go to 15
      endif
c
c...  Peter Pulay's famous DIIS
c
      if (residu.lt.diiscr.and.ncyc.ge.ntdiis) diisc = .true.
      if (diisc)
     * call diismp(cr(kmpc+1),cr(kmpz+1),indexp(5),ntotmp,indexp(3),
     *             iprimp)
c
      if (iprimp.ge.2000) write(6,609) (cr(kmpz+i),i=1,ntotmp)
609   format(' - residual-vector (after DIIS) -',/,(2x,10e12.5))
c
c...  check extrapolated residue
c
      residd = 0.0d0
      do i=1,ntotmp
         residd = max(residd,dabs(cr(kmpz+i)))
      end do
c
c...  r/D
c
      call rdedx(cr(kmpc+1),ntotmp,indexp(2),numci)
c
c...  level-shift
c
      if (pplev.ne.0.0d0)
     *   call vsadd(cr(kmpc+1),1,-pplev,cr(kmpc+1),1,ntotmp)
c
      call vvdv(ntotmp,cr(kmpz+1),cr(kmpz+1),cr(kmpc+1))
c
c...  check largest change in first-order wavefunction (not used)
c
      change = 0.0d0
      do i=1,ntotmp
         change = max(change,dabs(cr(kmpz+i)))
      end do
c
      call cpuwal(top,bot)
      ncyc = ncyc + 1
      write(6,666) ncyc,top,bot,change,residu,residd
666   format(' MP-cycle',i4,' at',f8.0,' wall',f8.3,' secs',
     *       ' change',e11.4,' residue',e11.4,' diis residue',e11.4)
      if (ncyc.ge.levelp) pplev = plevel(2)
c
c...   update c :  x = x + r/d
c
      call rdedx(cr(kmpc+1),ntotmp,indexp(3),numci)
      call vadd(cr(kmpz+1),1,cr(kmpc+1),1,cr(kmpc+1),1,ntotmp)
c
c     write(6,603) cr(kmpc+1)
      if (iprimp.ge.2000) write(6,604) (cr(kmpc+i),i=1,ntotmp)
c
      call wrt3(cr(kmpc+1),ntotmp,indexp(3),numci)
c
c...  calculate mp2 energy  e1h0
c
      call rdedx(cr(kmpz+1),ntotmp,indexp(1),numci)
      e1h0 = ddot(ntotmp,cr(kmpc+1),1,cr(kmpz+1),1)
c
      if (ncyc.lt.maxcmp) go to 10
c
c...   no convergence // print out best results // skip fancy print
c
c     write(6,603) cr(kmpc+1)
      cc00 = ddot(ntotmp,cr(kmpc+1),1,cr(kmpc+1),1)
c     write(6,667) e1h0,ehyll,cc00
c667   format(' *&*&*&*&*&*&*&* no convergence *&*&*&*&*&*&*&*',/,
c     *       ' &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&',/,
c     *       ' &&          Psi(1) is updated               &&',/,
c     *       ' &&     new E2(1) ',f16.13,       '          &&',/,
c     *       ' &&     old E2(h) ',f16.13,       '          &&',/,
c     *       ' &&      <C1|C1>  ',f8.5,'                   &&',/,
c     *       ' &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'/)
      write(6,667)
667   format(' *&*&*&*&*&*&*&* no convergence *&*&*&*&*&*&*&*')
      go to 100
c
c...  we converged the first-order wavefunction
c
15    call cpuwal(top,bot)
      ncyc = ncyc + 1
      write(6,610) ncyc,top,bot
610   format(' MP-convergence at cycle',i4,' at',f8.3,' wall',f8.3,
     *       ' secs ')
c
c...   get norm of  c-vector
c
      call rdedx(cr(kmpc+1),ntotmp,indexp(3),numci)
      cc00 = ddot(ntotmp,cr(kmpc+1),1,cr(kmpc+1),1)
c
c...  get overlap with reference function
c
      overla = cr(kmpc+1)
c
c...   if (overlap with psi(r) really big) print warning // error??
c
      write(6,611)overla
      write(*,*)'*** cc00 ',cc00
      if (dabs(overla).gt.dabs(critmp)) write(6,611) overla
611   format(' #@#@#@#@#@# <psi(1)|psi(r)> =',e11.3,' #@#@#@#@#@#')
100   continue
c
c...  expand the excitated-state vector to configuration-state vector
c
      call rdedx(cr(kcc+1),ntotmp,indexp(3),numci)
c...  the vector should have been orthogonal to C0.
      cr(kcc+1) = 0.0d0
      call wrt3(cr(kcc+1),ntotmp,indexp(ic(iord)),numci)
      call ccxpmp(cr(kc+1),cr(kcc+1),cr(kdc+1),cr(kvc+1),cr(kscr+1))
c
c...  restore original information in the common blocks to enable
c...  the use of genz.
c
      call ciretr
      nuse = max(khbz,khbc,khbb)+ntotal
      call corfai(zblank)
c
c...  move the csf-vector to proper position for genz.
c
      if(kc.gt.khbc) then
        do 80 i = 1, ntotal
          cr(khbc+i) = cr(kc+i)
80      continue
      else if(kc.lt.khbc) then
        do 90 i = 0, ntotal-1
          cr(khbc+ntotal-i) = cr(kc+ntotal-i)
90      continue
      endif
c
c...  retrieve config information to do genz
c
      call gciblk(mpcblk)
      call rdedx(conf,5*(nnst+nmin1+nmin2),mpcblk,numci)
c
c...  dump current vector to synchronize dumpfile contents with
c...  core contents
c
      call dumpcz(cr,cr(khbc+1),index(26))
c
c...  calculate  Z = H*C
c
      call genz(cr,conf,index(26),h11)
      if(instor) call dumpcz(cr,cr(khbz+1),index(27))
c
c...  get c,z-vector in proper format
c
      call getcz(cr,cr(khbc+1),index(26))
      call wrt3(cr(khbc+1),ntotal,index(26),numci)
      call getcz(cr,cr(khbz+1),index(27))
      call wrt3(cr(khbz+1),ntotal,index(27),numci)
c
c...  return to the contracted MP universe
c
      call mpretr
      call gmpblk(mpcblk)
      call rdedx(cr,5*(nmpv+nmp1+nmp2+nmp2d),mpcblk,numci)
c
c...  retrieve, contract, and store zero-th order z-vector
c
      call rdedx(cr(kc+1),ntotal,index(27),numci)
      if (iprimp.ge.2000) write(6,603) (cr(kc+i),i=1,ntotal)
c603   format(' - HC0 (CSF) -',/,(2x,10e12.5))
      call ccntmp(cr(kc+1),cr(kcc+1),cr(kdc+1),cr(kvc+1),cr(kscr+1))
      if (iprimp.ge.2000) write(6,602) (cr(kcc+i),i=1,ntotmp)
c602   format(' - HC0 (excited states) -',/,(2x,10e12.5))
      call wrt3(cr(kcc+1),ntotmp,indexp(iz(iord)),numci)
c
c...  increase order
c
      iord = iord + 1
c
      if(iord.le.min(mp+1,maxord)) go to 30
c
c...  END
c
c...  calculate the last order energy-term
c
      call rdedx(cr(kmpc+1),ntotmp,indexp(ic(1)),numci)
      call rdedx(cr(kmpz+1),ntotmp,indexp(iz(iord-1)),numci)
      epert(iord) = ddot(ntotmp,cr(kmpc+1),1,cr(kmpz+1),1)
c
c...  calculate final c,z-vectors
c
      call rdedx(cr(kmpc+1),ntotmp,indexp(ic(1)),numci)
      call rdedx(cr(kmpz+1),ntotmp,indexp(iz(1)),numci)
      cnorm = dsqrt(ddot(ntotmp,cr(kmpc+1),1,cr(kmpc+1),1))
      write(6,618)0,epert(1),epert(1),cnorm
      cnorm = dsqrt(ddot(ntotmp,cr(kmpc+1),1,cr(kmpc+1),1))
      write(6,620)1,cr(kmpz+1),epert(2),cnorm
618   format(//,3x,59('-'),/,
     *'   order       perturbation energy    correction',
     *'     corr.norm',
     */,3x,59('-'),//,3x,i4,3x,e24.15,3x,e13.7,3x,e9.3)
      do 40 i = 2,iord-1
        call rdedx(cr(kmpb+1),ntotmp,indexp(ic(i)),numci)
        cnorm = dsqrt(ddot(ntotmp,cr(kmpb+1),1,cr(kmpb+1),1))
        call daxpy(ntotmp,1.0d0,cr(kmpb+1),1,cr(kmpc+1),1)
        call rdedx(cr(kmpb+1),ntotmp,indexp(iz(i)),numci)
        call daxpy(ntotmp,1.0d0,cr(kmpb+1),1,cr(kmpz+1),1)
        write(6,620)i,cr(kmpz+1),epert(i+1),cnorm
40    continue
      call wrt3(cr(kmpc+1),ntotmp,indexp(ic(iord)),numci)
      call wrt3(cr(kmpz+1),ntotmp,indexp(iz(iord)),numci)
      if (iprimp.gt.1) then
        write(6,621)
621     format('MPn vector')
        call ccvpri(cr,cr(kmpc+1),1.0d-4)
      endif
      if (iprimp.gt.4) then
        dd = ddot(ntotmp,cr(kmpc+1),1,cr(kmpc+1),1)
        dd = 1.0d0/dsqrt(dd)
        call dscal(ntotmp,dd,cr(kmpc+1),1)
        write(6,622)
622     format('MPn vector RENORMALISED')
        call ccvpri(cr,cr(kmpc+1),1.0d-4)
      endif
620   format(3x,i4,3x,e24.15,3x,e13.7,3x,e9.3)
c
c...  calculate final energy
c
      call rdedx(cr(kmpc+1),ntotmp,indexp(ic(1)),numci)
      call rdedx(cr(kmpz+1),ntotmp,indexp(iz(iord)),numci)
      h11 = ddot(ntotmp,cr(kmpc+1),1,cr(kmpz+1),1)
c
c...  Expand c-vector and store it in right format
c
      call rdedx(cr(kcc+1),ntotmp,indexp(ic(iord)),numci)
      call ccxpmp(cr(kc+1),cr(kcc+1),cr(kdc+1),cr(kvc+1),cr(kscr+1))
      call ciretr
      if(kc.gt.khbc) then
        do 200 i = 1, ntotal
          cr(khbc+i) = cr(kc+i)
200     continue
      else if(kc.lt.khbc) then
        do 210 i = 0, ntotal-1
          cr(khbc+ntotal-i) = cr(kc+ntotal-i)
210      continue
      endif
      call gciblk(mpcblk)
      call rdedx(cr,5*(nnst+nmin1+nmin2),mpcblk,numci)
      call dumpcz(cr,cr(khbc+1),index(26))
c
c...  Return to contracted MP universe
c
      call mpretr
      call gmpblk(mpcblk)
      call rdedx(cr,5*(nmpv+nmp1+nmp2+nmp2d),mpcblk,numci)
c
c...  Expand z-vector and store it in right format
c
      call rdedx(cr(kcc+1),ntotmp,indexp(iz(iord)),numci)
      call ccxpmp(cr(kc+1),cr(kcc+1),cr(kdc+1),cr(kvc+1),cr(kscr+1))
      call ciretr
      if(kc.gt.khbc) then
        do 220 i = 1, ntotal
          cr(khbc+i) = cr(kc+i)
220     continue
      else if(kc.lt.khbc) then
        do 230 i = 0, ntotal-1
          cr(khbc+ntotal-i) = cr(kc+ntotal-i)
230      continue
      endif
      call gciblk(mpcblk)
      call rdedx(cr,5*(nnst+nmin1+nmin2),mpcblk,numci)
      call dumpcz(cr,cr(khbc+1),index(27))
c
c...  Calculate total energies
c
      eref  = epert(1) + epert(2)
      ecorr = 0.0d0
      do 240 i = 3, iord
        ecorr = ecorr + epert(i)
240   continue
      write(6,640) mp,abs(mpmod),descr(abs(mpmod)+1),nref,mspin,
     *             eref,ecorr,eref+ecorr
      if (mpmod.lt.0) write(6,*)
     *         ' ** H0 block diagonal in excitation classes         **'
640   format(//' $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$',
     *        /' $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$',
     *        /' $$             Multi-reference MP',i6,'            $$',
     *        /' $$                 model',i6,'                     $$',
     *        /' $$   ',              a39                  ,'       $$',
     *        /' $$          Nref ',i6,'   Spin ',i6,'              $$',
     *        /' $$                                                 $$',
     *        /' $$      E-reference   ',f18.10,        '           $$',
     *        /' $$      E-correlation ',f18.10,        '           $$',
     *        /' $$      --- total --- ',f18.10,        '           $$',
     *        /' $$                                                 $$',
     *        /' $$                                                 $$',
     *        /' $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$',
     *        /' $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$',
     *        /)
c      call cpuwal(top,bot)
c      write(6,601) iord-1,top,bot,e0mp+ecorr,epert(iord)
601   format(i8,f12.0,f12.3,2f18.10)
c     if (maxcyc.eq.0) return
c
c..   if (maxcyc.ne.0) there is more to come
c..   so then normalise c and z-vectors
c
c     if (.not.instor) then
         call getcz(cr,cr(khbc+1),index(26))
         call getcz(cr,cr(khbz+1),index(27))
c     end if
c
      cc00 = dsqrt(ddot(ntotal,cr(khbc+1),1,cr(khbc+1),1))
      cc00 = 1.0d0/cc00
      call dscal(ntotal,cc00,cr(khbc+1),1)
      call dscal(ntotal,cc00,cr(khbz+1),1)
c
      call dumpcz(cr,cr(khbc+1),index(26))
      call dumpcz(cr,cr(khbz+1),index(27))
c
      return
      end
      subroutine basmp(conf)
c
c...  set up addressing for mpgen etc.
c...  results in conf and possibly changed commons
c...  in store only
c...  note that aa-diagonals are seperate states
c...  the singlet/triplet n-2 are combined in c/z vector
c...  in conf(5,*) is config-classification :
c...  i   = 1 for singles and 2 for double excitation
c...  j = 1 + # holes in doubly occ. for vacuum
c...      4 + # holes ""   ""        for n-1
c...      7 + # holes in doubly occ. for n-2
c
      implicit real*8 (a-h,o-z), integer (i-n)
c
      integer conf,pack2
      dimension conf(5,*)
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
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
      common/mr_mp2l/ ivoc(62),nvoc,idoc(62),ndoc,kvcc(4),nvcc(4),
     *                ncmphn(10,18),ncmph1(10,18,8),ncmph2(10,18,8),
     *                ncmphd(10,18)
      common /mr_mp2con/ nsmpv,nsmp1,nsmp2,nsmp2d,ntotmp,
     *                   nmpv,nmp1,nmp2,nmp2d
c
c
      integer iro, iky, niky, m1e, m2e, m1ia, m2coul
      integer m2exch, m2fock, m2psup, m2qsup
      integer mbase, nbase, npairi, norbi, mbasi, mnbasi
      integer iap, mbfock, mbexch, mbcoul, mbpsup, mbqsup
      integer norbas, m2exab, m2em1e
      common /helpr/ iro(nd200),iky(nd200),niky(8),m1e,m2e,m1ia,m2coul,
     +   m2exch,m2fock,m2psup,m2qsup,
     +   mbase(288),nbase(8),npairi(36),norbi(8),mbasi(8),mnbasi(8),
     +   iap(62),mbfock(62),mbexch(36),mbcoul(36),mbpsup(36),mbqsup(36),
     +   norbas(nd200),m2exab,m2em1e
c
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
c
      integer nint, next, ninnex, nwm1, nw, nwp1, nintwo
      integer nint21, nintp1
      common /orb/ nint,next,ninnex,nwm1,nw,nwp1,nintwo,
     +             nint21,nintp1
c
      common/stoctl/khbill,khbz,khbc,khbb,khbzd,khbzs,khbzt,khbcd,
     *khbcs,khbct,kh,khb1,khb2,khb3,khb4,khb5,khb6,khb7,khb8
     *,khb44,khb55,khb66,khb88,khb111,khb222,khb333
c
      integer nconst, ncond, maxpa, maxpat, maxsin, maxepa
      integer mxdvpa, mxddpa, nspips, nspipt, mxsspa, mxttpa
      integer mxstpa, nsitot, ntrtot
      integer mx12pa, mxvvpa, mxv2pa, maxpsq, max12
      common /auxh/ nconst(8),ncond(8),maxpa(4),maxpat,maxsin(3),
     +              maxepa(4),mxdvpa,mxddpa,nspips(8),nspipt(8),
     +              mxsspa,mxttpa,mxstpa,nsitot(8),ntrtot(8),
     +              mx12pa,mxvvpa,mxv2pa,maxpsq,max12
c
c
      integer nirr,mult,ifrep,nirr8,jfrep
      common /table/ nirr,mult(8,8),ifrep,nirr8,jfrep
c
c
      integer nuse, ngot
      common /corctl/ nuse, ngot
c
      common/mr_comcci/icci,iclass(3,3,3),istcc,iselcc
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
c
      character*9 texi(3)
      character*6 texj(4)
      character*2 ts
c
      pack2 (i1,i2)  = ior(ishft(i1,32),i2)
      data texi/'reference',' single  ',' double  '/
      data texj/'vacuum',' n-1  ',' n-2  ',' n-2aa'/
c
c...  singlets are lower triangles as well ...
c...  singlets and triplets are the same vector !!
c
      do 1 i=1,nirr
1     ibassi(i,1) = ibastr(i,1)
c
      maxpat = 0
c
c... vacuum
c
      kk = 0
      loop = 0
      do 30 is=1,2+nirr*2
         if (is.eq.1) then
            nhole = 0
         else if (is.le.nirr+1) then
            nhole = 1
         else
            nhole = 2
         end if
         nset = ncmphn(3,is)
         maxpat = max(maxpat,ncmphn(6,is),ncmphn(7,is),ncmphn(8,is))
         do 20 iset=1,nset
            do 10 iex=6,8
               nex = ncmphn(iex,is)
               if (nex.gt.0) then
                  loop = loop + 1
                  conf(1,loop) = kk
                  conf(3,loop) = nex
                  conf(4,loop) = 1
                  conf(5,loop) = pack2(iex-6,1+nhole)
                  kk = kk + nex
               end if
10          continue
20       continue
30    continue
c
      nsmpv = kk
      nmpv = loop
c
c...  doublets
c
      ll = 0
c
      do 70 ijrr=1,nirr
         norex = norbe(ijrr)
         do 60 is=1,2+nirr*2
            if (is.eq.1) then
               nhole = 0
            else if (is.le.nirr+1) then
               nhole = 1
            else
               nhole = 2
            end if
            nset = ncmph1(3,is,ijrr)
            maxpat = max(maxpat,ncmph1(6,is,ijrr),
     *                    ncmph1(7,is,ijrr),ncmph1(8,is,ijrr))
            do 50 iset=1,nset
               do 40 iex=6,8
                  nex = ncmph1(iex,is,ijrr)
                  if (nex.gt.0) then
                     loop = loop + 1
                     conf(1,loop) = ll
                     conf(3,loop) = nex
                     conf(4,loop) = ijrr
                     conf(5,loop) = pack2(iex-6,4+nhole)
                     ll = ll + nex*norex
                  end if
40             continue
50          continue
60       continue
70    continue
c
      nsmp1 = ll
      nmp1 = loop - nmpv
c
c...  (n-2)  singlet+triplets  (diagonal seperate)
c
      kkk = 0
c
      do 110 ijrr=1,nirr
         norsi = norbtr(ijrr)
         do 100 is=1,2+nirr*2
            if (is.eq.1) then
               nhole = 0
            else if (is.le.nirr+1) then
               nhole = 1
            else
               nhole = 2
            end if
            nset = ncmph2(3,is,ijrr)
            maxpat = max(maxpat,ncmph2(6,is,ijrr),
     *                    ncmph2(7,is,ijrr),ncmph2(8,is,ijrr))
            do 90 iset=1,nset
               do 80 iex=6,8
                  nex = ncmph2(iex,is,ijrr)
                  if (nex.gt.0) then
                     nsex = nex
                     ntex = nex
                     if (ncmph2(4,is,ijrr).eq.0) nsex = 0
                     if (ncmph2(5,is,ijrr).eq.0) ntex = 0
                     loop = loop + 1
                     conf(1,loop) = pack2(kkk,kkk)
                     conf(3,loop) = pack2(ntex,nsex)
                     conf(4,loop) = ijrr
                     conf(5,loop) = pack2(iex-6,7+nhole)
                     kkk = kkk + nex*norsi
                  end if
80             continue
90          continue
100      continue
110   continue
c
      nsmp2 = kkk
      nmp2 = loop - nmpv - nmp1
c
c...  aa-(n-2)-diagonal seperate
c...  addressing/csf in conf(2,loop)
c
      kkk = 0
      lll = 0
      do 140 is=1,2+nirr*2
         if (is.eq.1) then
            nhole = 0
         else if (is.le.nirr+1) then
            nhole = 1
         else
            nhole = 2
         end if
         nset = ncmphd(3,is)
         maxpat = max(maxpat,ncmphd(6,is),ncmphd(7,is),ncmphd(8,is))
         do 130 iset=1,nset
            do 120 iex=6,8
               nex = ncmphd(iex,is)
               if (nex.gt.0) then
                  loop = loop + 1
                  conf(1,loop) = pack2(0,kkk)
                  conf(2,loop) = pack2(0,lll)
                  conf(3,loop) = pack2(0,nex)
                  conf(4,loop) = 1
                  conf(5,loop) = pack2(iex-6,7+nhole)
                  kkk = kkk + nex*next
                  lll = lll + nex
               end if
120         continue
130      continue
140   continue
c
      nmp2d = loop - (nmpv+nmp1+nmp2)
      nsmp2d = kkk
c
      ntotmp = nsmpv + nsmp1 + nsmp2 + nsmp2d
c
c     nblock = (ntotmp-1)/511 + 1
c
c...  we now know the length of conf => adapt khbill
c
      khbill = loop*5
c
      khbz=khbill
      khbzd=khbz+nsmpv
      khbzs=khbzd+nsmp1
      khbzt=khbzs
c      khddz = khbzs + nsmp2
      khbc = khbz + ntotmp
      khbcd=khbc+nsmpv
      khbcs=khbcd+nsmp1
      khbct=khbcs
c      khddc = khbcs + nsmp2
c
      khbb = khbc + ntotmp
c...   vvijkl,ddijkl,ssijkl
      kend = khbb + 2*maxpat**2
c...   dvmpia
      kend = max(kend,khbb+m1ia+maxpat**2)
c...   sdmpia
      khb111 = khbb + m1ia
      khb222 = khb111 + 2*maxpat**2
      khb333 = khb222 + next**2
      kend = max(kend,khb333+next*2*maxpat)
c...   ddmpab
      kend = max(kend,khbb+norbsq(1))
c...   ssmpab
      khb7 = khbb
      khb3 = khb7 + norbsq(1)
      khb4 = khb3 + maxpat**2
      khb8 = khb4 + maxpat**2
      kend = max(kend,khb8+nubsq+mubsq+mubsq)
      kend = max(kend,khb8+next)
c
      nuse = kend
      call corfai('mp')
c
      if(odebug(40)) then
c      debug print
       write(6,*) 'basmp: conf'
       call confprt(conf)
      endif
c
      if (iprimp.lt.-1) return
c
      write(6,601) loop,ntotmp,nmpv,nmp1,nmp2,nmp2d,
     *             nsmpv,nsmp1,nsmp2,nsmp2d
601   format(//,' // excited state classification for',i7,' models '
     *         ,i9,' states \\',
     *       /,3x,'  models  # vac',i5,' # n-1',i6,' # n-2',i8,
     *         ' # n-2d',i6,
     *       /,3x,'  states  # vac',i5,' # n-1',i6,' # n-2',i8,
     *         ' # n-2d',i6,/)
      if (iprimp.lt.100) return
      write(6,602)
602   format('     #        type        # holes symmetry  # model-s ',
     *       ' # states ')
      idumm = 1
      do 600  i=1,loop
         call upack2(conf(5,i),it,jj)
         jt = (jj-1)/3+1
         jj = jj - (jt-1)*3 - 1
         call upack2(conf(3,i),ntex,nsex)
         nex = max(ntex,nsex)
         call upack2(conf(4,i),idum,is)
         if (i.gt.nmpv+nmp1+nmp2) jt = 4
         if (jt.eq.1) nnn = nex
         if (jt.eq.2) nnn = nex*norbe(is)
         if (jt.eq.3) nnn = nex*norbtr(is)
         if (jt.eq.4) nnn = nex*next
         ts = ' '
         if (jt.eq.3.and.nsex.eq.0) ts = '_t'
         if (jt.eq.3.and.ntex.eq.0) ts = '_s'
         write(6,603) i,texi(it+1),texj(jt),ts,jj,is,nex,nnn
     *               ,idumm,idumm+nnn-1
603      format(i6,2x,a9,1x,a6,a2,i6,2x,i6,4x,i6,4x,i6,
     *          2x,'(',i7,' -',i7,')')
         idumm = idumm + nnn
600   continue
c
      return
      end
      subroutine ccntmp(c,cc,dc,vc,scr)
c
c...  transform csf-c/z-vector to excited state c/z-vector
c...  by calling contmp symmetry blocked
c...  for n-2 : add singlet+triplet and get aa-diagonal in dc
c...  results in cc/dc
c
      implicit real*8 (a-h,o-z), integer (i-n)
c
      dimension c(*),cc(*),vc(*),dc(*),scr(*)
c
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
      common/mr_mp2l/ ivoc(62),nvoc,idoc(62),ndoc,kvcc(4),nvcc(4),
     *                ncmphn(10,18),ncmph1(10,18,8),ncmph2(10,18,8),
     *                ncmphd(10,18)
c
      integer nval, ndoub, nsingl, ntripl, ntotal, ibig, jbig, kbig
      integer nvd, nvds, nsitri, lbig, mbig, nbig, nref0
      common /cntrl/ nval,ndoub,nsingl,ntripl,ntotal,ibig,jbig,kbig,
     +               nvd,nvds,nsitri,lbig,mbig,nbig,nref0
c
c
      integer nint, next, ninnex, nwm1, nw, nwp1, nintwo
      integer nint21, nintp1
      common /orb/ nint,next,ninnex,nwm1,nw,nwp1,nintwo,
     +             nint21,nintp1
c
c
      integer nirr,mult,ifrep,nirr8,jfrep
      common /table/ nirr,mult(8,8),ifrep,nirr8,jfrep
c
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
c
      integer nconst, ncond, maxpa, maxpat, maxsin, maxepa
      integer mxdvpa, mxddpa, nspips, nspipt, mxsspa, mxttpa
      integer mxstpa, nsitot, ntrtot
      integer mx12pa, mxvvpa, mxv2pa, maxpsq, max12
      common /auxh/ nconst(8),ncond(8),maxpa(4),maxpat,maxsin(3),
     +              maxepa(4),mxdvpa,mxddpa,nspips(8),nspipt(8),
     +              mxsspa,mxttpa,mxstpa,nsitot(8),ntrtot(8),
     +              mx12pa,mxvvpa,mxv2pa,maxpsq,max12
c
c
      integer nxblk, index, numci, mblkci, ispbig, ispsma
      common/disktl/nxblk,index(199),numci,mblkci,ispbig,ispsma
c
c... vacuum
c
      kc = 1
      kcc = 1
      call rdedx(vc,nvcc(1),kvcc(1),numci)
      kvc = 1
      do 10 is=1,2+nirr*2
         ncsf = ncmphn(2,is)
         nset = ncmphn(3,is)
         nst = ncmphn(9,is)
         call vclr(cc(kcc),1,nst*nset)
         call contmp(c(kc),cc(kcc),1,vc(kvc),ncsf,ncsf,nst,nset)
         kc = kc + ncsf*nset
         kcc = kcc + nst*nset
         kvc = kvc + ncsf*nst*nset
10    continue
c
c...  doublet (n-1)
c
      call rdedx(vc,nvcc(2),kvcc(2),numci)
      kvc = 1
      do 30 ijrr=1,nirr
         nex = norbe(ijrr)
         do 20 is=1,2+nirr*2
            ncsf = ncmph1(2,is,ijrr)
            nset = ncmph1(3,is,ijrr)
            nst = ncmph1(9,is,ijrr)
            call vclr(cc(kcc),1,nst*nset*nex)
            call contmp(c(kc),cc(kcc),nex,vc(kvc),ncsf,ncsf,nst,nset)
            kc = kc + ncsf*nset*nex
            kcc = kcc + nst*nset*nex
            kvc = kvc + ncsf*nst*nset
20       continue
30    continue
c
c...  (n-2)  singlet and triplets   add
c
c...  get symmetry 1 diagonal out  and contract
c
      call rdedx(vc,nvcc(4),kvcc(4),numci)
      kvc = 1
      if (nspips(1).gt.0) call aagp(c(kc),nspips(1),scr,next,'get')
      ksc = 1
      kdc = 1
      do 40 is=1,2+nirr*2
         ncsf = ncmphd(2,is)
         nset = ncmphd(3,is)
         nst = ncmphd(9,is)
         call vclr(dc(kdc),1,nst*nset*next)
         call contmp(scr(ksc),dc(kdc),next,vc(kvc),ncsf,ncsf,nst,nset)
         ksc = ksc + ncsf*nset*next
         kdc = kdc + nst*nset*next
         kvc = kvc + ncsf*nst*nset
40    continue
c
c...  general (n-2)
c
      call rdedx(vc,nvcc(3),kvcc(3),numci)
      kvc = 1
c
c...  symmetry 1 off diagonal   (first compress than contract)
c...  do singlet/triplet and then add
c...  originally triplets start nsingl after singlets
c
      if (nspips(1).ne.0)
     *call sin1mp(c(kc),norbsi(1),scr,norbtr(1),nspips(1),'com')
      kc = kc + nspips(1)*next
      call dcopy(norbtr(1)*nspips(1),scr,1,c(kc),1)
      nstot = nsingl - nspips(1)*next
      kctr = kc + nstot
c
      do 60 ijrr=1,nirr
         nex = norbtr(ijrr)
         do 50 is=1,2+nirr*2
            ncsf = ncmph2(2,is,ijrr)
            nset = ncmph2(3,is,ijrr)
            nsin = ncmph2(4,is,ijrr)
            ntri = ncmph2(5,is,ijrr)
            nst = ncmph2(9,is,ijrr)
c..    contract
            call vclr(cc(kcc),1,nset*nst*nex)
            if (nsin.gt.0)
     *      call contmp(c(kc),cc(kcc),nex,vc(kvc),ncsf,nsin,nst,nset)
            if (ntri.gt.0)
     *      call contmp(c(kctr),cc(kcc),nex,vc(kvc+nsin),
     *                  ncsf,ntri,nst,nset)
            kc = kc + nsin*nset*nex
            kctr = kctr + ntri*nset*nex
            kcc = kcc + nst*nset*nex
            kvc = kvc + ncsf*nst*nset
50       continue
60    continue
c
      return
      end
      subroutine gtczmp(c,ns,conf,indcz,act)
c
c...  a flexible cz (csf basis) vectors getter/putter for mp
c...  disk store relative addresses for c-z in conf(2,i)
c...  c   : the c-vector to be returned
c...  ns    : number of elements to be retrieved (0=all)
c...          for symsin/symtri ns give the symmetry to be retrieved
c...  act   : action to be taken (see if construction below)
c...  indcz : start block on disk of vector
c...  conf  : the conf array with symbolic info
c...          conf has to be in its direct ci state
c
c...  mixing gtczmp and ptczmp is not safe
c
      implicit real*8 (a-h,o-z)
c
      integer nval, ndoub, nsingl, ntripl, ntotal, ibig, jbig, kbig
      integer nvd, nvds, nsitri, lbig, mbig, nbig, nref0
      common /cntrl/ nval,ndoub,nsingl,ntripl,ntotal,ibig,jbig,kbig,
     +               nvd,nvds,nsitri,lbig,mbig,nbig,nref0
c
c
      integer nref, nnst, nmin1, nmin2, mspin, nelec
      integer nstrt1, nlast1, nstrt2, nlast2, nmin12, nminv2, nminv1
      integer nlst21, nlastx
      common /ccntl/ nref,nnst,nmin1,nmin2,mspin,nelec,
     +               nstrt1,nlast1,nstrt2,nlast2,nmin12,nminv2,nminv1,
     +               nlst21,nlastx
c
c
      integer nxblk, index, numci, mblkci, ispbig, ispsma
      common/disktl/nxblk,index(199),numci,mblkci,ispbig,ispsma
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
c
      logical get
      character*(*) act
      integer conf
      dimension c(*),conf(5,*)
      save isconf,itconf
c
      get = .true.
c
      go to 1
c
      entry ptczmp(c,ns,conf,indcz,act)
c...  put corresponding to gtczmp get
      get = .false.
c
1     if (act.eq.'vacuum') then
         if (get) then
            call rdedx(c,nval,indcz,numci)
         else
            call wrt3(c,nval,indcz,numci)
         end if
         else if (act.eq.'doublet') then
         ibl = indcz + conf(2,1)
         if (get) then
            call rdedx(c,ndoub,ibl,numci)
         else
            call wrt3(c,ndoub,ibl,numci)
         end if
      else if (act.eq.'singlet'.or.act.eq.'symsin') then
c...     retrieve all singlet states (up to ns or of symmetry ns)
         kk = 1
         do i=isconf,nlast2
            isi = conf(4,i)
            if (act.eq.'symsin') then
               if (isi.ne.ns) go to 10
            else 
               if (nn.eq.ns) go to 10
               if (nn.gt.ns) call caserr('trouble in singlet gtczmp')
            end if
            call upack2(conf(3,i),lamt,lams)
            nn=lams*norbsi(isi)
            call upack2(conf(2,i),ibt,ibs)
            if (get) then
               call rdedx(c(kk),nn,ibs+indcz,numci)
            else
               call wrt3(c(kk),nn,ibs+indcz,numci)
            end if
            kk = kk + nn
         end do
10       isconf = i 
      else if (act.eq.'triplet'.or.act.eq.'symtri') then
c...     retrieve all triplet states (up to ns or of symmetry ns)
         kk = 1
         do i=itconf,nlast2
            isi = conf(4,i)
            if (act.eq.'symtri') then
               if (isi.ne.ns) go to 20
            else 
               if (nn.eq.ns) go to 20
               if (nn.gt.ns) call caserr('trouble in triplet gtczmp')
            end if
            call upack2(conf(3,i),lamt,lams)
            nn=lamt*norbtr(isi)
            call upack2(conf(2,i),ibt,ibs)
            if (get) then
               call rdedx(c(kk),nn,ibt+indcz,numci)
            else
               call wrt3(c(kk),nn,ibt+indcz,numci)
            end if
            kk = kk + nn
         end do
20       itconf = i
      else if (act.eq.'all') then
         call getcz(conf,c,indcz)
      else if (act.eq.'reset') then
         isconf = nstrt2
         itconf = nstrt2
c...     read conf
         call gciblk(iblk)
         call rdedx(conf,5*nlast2,iblk,numci)
      else
         call caserr(' unrecognised action in gtczmp/ptczmp ')
      end if
c
      end
c
      subroutine ccntmp2(c,cc,dc,vc,scr,conf,indcz)
c
c...  transform csf-c/z-vector to excited state c/z-vector
c...  by calling contmp symmetry blocked
c...  for n-2 : add singlet+triplet and get aa-diagonal in dc
c...  results in cc/dc
c...  reduced core usage; requires only 1 space/spin symmetry in core
c...  uses gtczmp,aagp,sin1mp2,contmp
c
      implicit real*8 (a-h,o-z), integer (i-n)
c
      dimension c(*),cc(*),vc(*),dc(*),scr(*)
      integer conf
      dimension conf(5,*)
c
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
      common/mr_mp2l/ ivoc(62),nvoc,idoc(62),ndoc,kvcc(4),nvcc(4),
     *                ncmphn(10,18),ncmph1(10,18,8),ncmph2(10,18,8),
     *                ncmphd(10,18)
c
c
      integer nval, ndoub, nsingl, ntripl, ntotal, ibig, jbig, kbig
      integer nvd, nvds, nsitri, lbig, mbig, nbig, nref0
      common /cntrl/ nval,ndoub,nsingl,ntripl,ntotal,ibig,jbig,kbig,
     +               nvd,nvds,nsitri,lbig,mbig,nbig,nref0
c
c
      integer nxblk, index, numci, mblkci, ispbig, ispsma
      common/disktl/nxblk,index(199),numci,mblkci,ispbig,ispsma
c
      integer nint, next, ninnex, nwm1, nw, nwp1, nintwo
      integer nint21, nintp1
      common /orb/ nint,next,ninnex,nwm1,nw,nwp1,nintwo,
     +             nint21,nintp1
c
c
      integer nirr,mult,ifrep,nirr8,jfrep
      common /table/ nirr,mult(8,8),ifrep,nirr8,jfrep
c
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
c
      integer nconst, ncond, maxpa, maxpat, maxsin, maxepa
      integer mxdvpa, mxddpa, nspips, nspipt, mxsspa, mxttpa
      integer mxstpa, nsitot, ntrtot
      integer mx12pa, mxvvpa, mxv2pa, maxpsq, max12
      common /auxh/ nconst(8),ncond(8),maxpa(4),maxpat,maxsin(3),
     +              maxepa(4),mxdvpa,mxddpa,nspips(8),nspipt(8),
     +              mxsspa,mxttpa,mxstpa,nsitot(8),ntrtot(8),
     +              mx12pa,mxvvpa,mxv2pa,maxpsq,max12
c
c
      logical trip
c
c... vacuum
c
      call gtczmp(c,0,conf,indcz,'reset')
c
      kc = 1
      kcc = 1
      call rdedx(vc,nvcc(1),kvcc(1),numci)
      call gtczmp(c,0,conf,indcz,'vacuum')
      kvc = 1
      do 10 is=1,2+nirr*2
         ncsf = ncmphn(2,is)
         nset = ncmphn(3,is)
         nst = ncmphn(9,is)
         call vclr(cc(kcc),1,nst*nset)
         call contmp(c(kc),cc(kcc),1,vc(kvc),ncsf,ncsf,nst,nset)
         kc = kc + ncsf*nset
         kcc = kcc + nst*nset
         kvc = kvc + ncsf*nst*nset
10    continue
c
c...  doublet (n-1)
c
      call rdedx(vc,nvcc(2),kvcc(2),numci)
      call gtczmp(c,0,conf,indcz,'doublet')
      kc = 1
      kvc = 1
      do 30 ijrr=1,nirr
         nex = norbe(ijrr)
         do 20 is=1,2+nirr*2
            ncsf = ncmph1(2,is,ijrr)
            nset = ncmph1(3,is,ijrr)
            nst = ncmph1(9,is,ijrr)
            call vclr(cc(kcc),1,nst*nset*nex)
            call contmp(c(kc),cc(kcc),nex,vc(kvc),ncsf,ncsf,nst,nset)
            kc = kc + ncsf*nset*nex
            kcc = kcc + nst*nset*nex
            kvc = kvc + ncsf*nst*nset
20       continue
30    continue
c
c...  (n-2)  singlet and triplets   add
c
c...  get symmetry 1 diagonal out  and contract
c
      call rdedx(vc,nvcc(4),kvcc(4),numci)
      call gtczmp(c,1,conf,indcz,'symsin')
      kc = 1
      kvc = 1
      if (nspips(1).gt.0) call aagp(c(kc),nspips(1),scr,next,'get')
      ksc = 1
      kdc = 1
      do 40 is=1,2+nirr*2
         ncsf = ncmphd(2,is)
         nset = ncmphd(3,is)
         nst = ncmphd(9,is)
         call vclr(dc(kdc),1,nst*nset*next)
         call contmp(scr(ksc),dc(kdc),next,vc(kvc),ncsf,ncsf,nst,nset)
         ksc = ksc + ncsf*nset*next
         kdc = kdc + nst*nset*next
         kvc = kvc + ncsf*nst*nset
40    continue
c
c...  general (n-2)
c
      trip = .false.
      kcctr = kcc
      call rdedx(vc,nvcc(3),kvcc(3),numci)
      kvc = 1
c
c...  symmetry 1 off diagonal   (first compress than contract)
c...  do singlet/triplet and then add
c...  originally triplets start nsingl after singlets
c
      if (nspips(1).ne.0) call sin1mp2(c,nspips(1),'com')
c
100   continue
c
      do 60 ijrr=1,nirr
         if (trip) then
            call gtczmp(c,ijrr,conf,indcz,'symtri')
         else
            if (ijrr.ne.1) call gtczmp(c,ijrr,conf,indcz,'symsin')
         end if
         kc = 1
         nex = norbtr(ijrr)
         do 50 is=1,2+nirr*2
            ncsf = ncmph2(2,is,ijrr)
            nset = ncmph2(3,is,ijrr)
            nsin = ncmph2(4,is,ijrr)
            nst = ncmph2(9,is,ijrr)
c..    contract
            if (.not.trip) then
               call vclr(cc(kcc),1,nset*nst*nex)
               if (nsin.gt.0) 
     *         call contmp(c(kc),cc(kcc),nex,vc(kvc),ncsf,nsin,nst,nset)
               kc = kc + nsin*nset*nex
            else
               ntri = ncmph2(5,is,ijrr)
               if (ntri.gt.0) 
     *         call contmp(c(kc),cc(kcc),nex,vc(kvc+nsin),
     *                     ncsf,ntri,nst,nset)
               kc = kc + ntri*nset*nex
            end if
            kcc = kcc + nst*nset*nex
            kvc = kvc + ncsf*nst*nset
50       continue
60    continue
c
      if (.not.trip) then
         trip = .true.
         kvc = 1
         kcc = kcctr
         go to 100
      end if
c
      return
      end
      subroutine contmp(c,cc,nex,vcc,nvcc,ncsf,nst,nset)
c
c...  contract vectors c  using vcc
c...  yielding vectors cc (external dimension nex)
c...  (singlet/triplet together)
c...  loop over sets here
c
      implicit real*8 (a-h,o-z), integer(i-n)
c
      dimension c(nex,ncsf,nset),cc(nex,nst,nset),vcc(nvcc,nst,nset)
c
      do 100 i=1,nset
         call mxmb(c(1,1,i),1,nex,vcc(1,1,i),1,nvcc,
     *             cc(1,1,i),1,nex,nex,ncsf,nst)
100   continue
c
      return
      end
      subroutine ccxpmp(c,cc,dc,vc,scr)
c
c...  transform excited state c/z-vector to csf c/z-vector
c...  by calling cexpmp symmetry blocked
c...  for n-2 : add singlet+triplet and get aa-diagonal in dc
c
c     c   : resulting csf-vector
c     cc  : input excited state vector
c     vc  : storage for contraction matrices
c     scr : scratch (workspace)
c
c...  derived from ccntmp
c
      implicit real*8 (a-h,o-z), integer (i-n)
c
      dimension c(*),cc(*),vc(*),dc(*),scr(*)
c
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
      common/mr_mp2l/ ivoc(62),nvoc,idoc(62),ndoc,kvcc(4),nvcc(4),
     *                ncmphn(10,18),ncmph1(10,18,8),ncmph2(10,18,8),
     *                ncmphd(10,18)
c
c
      integer nval, ndoub, nsingl, ntripl, ntotal, ibig, jbig, kbig
      integer nvd, nvds, nsitri, lbig, mbig, nbig, nref0
      common /cntrl/ nval,ndoub,nsingl,ntripl,ntotal,ibig,jbig,kbig,
     +               nvd,nvds,nsitri,lbig,mbig,nbig,nref0
c
c
      integer nint, next, ninnex, nwm1, nw, nwp1, nintwo
      integer nint21, nintp1
      common /orb/ nint,next,ninnex,nwm1,nw,nwp1,nintwo,
     +             nint21,nintp1
c
c
      integer nirr,mult,ifrep,nirr8,jfrep
      common /table/ nirr,mult(8,8),ifrep,nirr8,jfrep
c
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
c
      integer nconst, ncond, maxpa, maxpat, maxsin, maxepa
      integer mxdvpa, mxddpa, nspips, nspipt, mxsspa, mxttpa
      integer mxstpa, nsitot, ntrtot
      integer mx12pa, mxvvpa, mxv2pa, maxpsq, max12
      common /auxh/ nconst(8),ncond(8),maxpa(4),maxpat,maxsin(3),
     +              maxepa(4),mxdvpa,mxddpa,nspips(8),nspipt(8),
     +              mxsspa,mxttpa,mxstpa,nsitot(8),ntrtot(8),
     +              mx12pa,mxvvpa,mxv2pa,maxpsq,max12
c
c
      integer nxblk, index, numci, mblkci, ispbig, ispsma
      common/disktl/nxblk,index(199),numci,mblkci,ispbig,ispsma
c
c... vacuum
c
      kc = 1
      kcc = 1
      call rdedx(vc,nvcc(1),kvcc(1),numci)
      kvc = 1
      do 10 is=1,2+nirr*2
         ncsf = ncmphn(2,is)
         nset = ncmphn(3,is)
         nst  = ncmphn(9,is)
         call vclr(c(kc),1,ncsf*nset)
         call cexpmp(c(kc),cc(kcc),1,vc(kvc),ncsf,ncsf,nst,nset)
         kc  = kc  + ncsf*nset
         kcc = kcc + nst*nset
         kvc = kvc + ncsf*nst*nset
10    continue
c
c...  doublet (n-1)
c
      call rdedx(vc,nvcc(2),kvcc(2),numci)
      kvc = 1
      do 30 ijrr=1,nirr
         nex = norbe(ijrr)
         do 20 is=1,2+nirr*2
            ncsf = ncmph1(2,is,ijrr)
            nset = ncmph1(3,is,ijrr)
            nst  = ncmph1(9,is,ijrr)
            call vclr(c(kc),1,ncsf*nset*nex)
            call cexpmp(c(kc),cc(kcc),nex,vc(kvc),ncsf,ncsf,nst,nset)
            kc  = kc + ncsf*nset*nex
            kcc = kcc + nst*nset*nex
            kvc = kvc + ncsf*nst*nset
20       continue
30    continue
c
c...  general (n-2)
c
      call rdedx(vc,nvcc(3),kvcc(3),numci)
      kvc = 1
c
c...  symmetry 1 off diagonal   (first compress than contract)
c...  do singlet/triplet and then add
c...  originally triplets start nsingl after singlets
c
      kc    = kc + nspips(1)*next
      nstot = nsingl - nspips(1)*next
      kctr  = kc + nstot
c...  <kcs> temporary variable, <kc> needed later.
      kcs   = kc
c
      do 60 ijrr=1,nirr
         nex = norbtr(ijrr)
         do 50 is=1,2+nirr*2
            ncsf = ncmph2(2,is,ijrr)
            nset = ncmph2(3,is,ijrr)
            nsin = ncmph2(4,is,ijrr)
            ntri = ncmph2(5,is,ijrr)
            nst  = ncmph2(9,is,ijrr)
c..    expand
            call vclr(c(kcs),1,nset*nsin*nex)
            call vclr(c(kctr),1,nset*ntri*nex)
            if (nsin.gt.0)
     *      call cexpmp(c(kcs),cc(kcc),nex,vc(kvc),ncsf,nsin,nst,nset)
            if (ntri.gt.0)
     *      call cexpmp(c(kctr),cc(kcc),nex,vc(kvc+nsin),
     *                  ncsf,ntri,nst,nset)
            kcs  = kcs + nsin*nset*nex
            kctr = kctr + ntri*nset*nex
            kcc  = kcc + nst*nset*nex
            kvc  = kvc + ncsf*nst*nset
50       continue
60    continue
      call dcopy(norbtr(1)*nspips(1),c(kc),1,scr,1)
      kc = kc - nspips(1)*next
      if (nspips(1).ne.0)
     *call sin1mp(c(kc),norbsi(1),scr,norbtr(1),nspips(1),'exp')
c
c...  (n-2)  singlet and triplets   add
c
c...  get symmetry 1 diagonal out  and contract
c
      call rdedx(vc,nvcc(4),kvcc(4),numci)
      kvc = 1
      ksc = 1
      kdc = 1
      do 40 is=1,2+nirr*2
         ncsf = ncmphd(2,is)
         nset = ncmphd(3,is)
         nst  = ncmphd(9,is)
         call vclr(scr(ksc),1,ncsf*nset*next)
         call cexpmp(scr(ksc),dc(kdc),next,vc(kvc),ncsf,ncsf,nst,nset)
         ksc = ksc + ncsf*nset*next
         kdc = kdc + nst*nset*next
         kvc = kvc + ncsf*nst*nset
40    continue
      if (nspips(1).gt.0) call aagp(c(kc),nspips(1),scr,next,'put')
      return
      end
      subroutine ccxpmp2(c,cc,dc,vc,scr,conf,indcz)
c
c...  transform excited state c/z-vector to csf c/z-vector
c...  by calling cexpmp symmetry blocked
c...  requiring only one space/spin type wordking using getcz format
c
c     c   : resulting csf-vector
c     cc  : input excited state vector
c     vc  : storage for contraction matrices
c     scr : scratch (workspace)
c
c...  derived from ccntmp2
c
      implicit real*8 (a-h,o-z), integer (i-n)
c
      dimension c(*),cc(*),vc(*),dc(*),scr(*)
c
      integer conf
      dimension conf(5,*)
c
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
      common/mr_mp2l/ ivoc(62),nvoc,idoc(62),ndoc,kvcc(4),nvcc(4),
     *                ncmphn(10,18),ncmph1(10,18,8),ncmph2(10,18,8),
     *                ncmphd(10,18)
c
c
      integer nval, ndoub, nsingl, ntripl, ntotal, ibig, jbig, kbig
      integer nvd, nvds, nsitri, lbig, mbig, nbig, nref0
      common /cntrl/ nval,ndoub,nsingl,ntripl,ntotal,ibig,jbig,kbig,
     +               nvd,nvds,nsitri,lbig,mbig,nbig,nref0
c
c
      integer nint, next, ninnex, nwm1, nw, nwp1, nintwo
      integer nint21, nintp1
      common /orb/ nint,next,ninnex,nwm1,nw,nwp1,nintwo,
     +             nint21,nintp1
c
c
      integer nirr,mult,ifrep,nirr8,jfrep
      common /table/ nirr,mult(8,8),ifrep,nirr8,jfrep
c
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
c
      integer nconst, ncond, maxpa, maxpat, maxsin, maxepa
      integer mxdvpa, mxddpa, nspips, nspipt, mxsspa, mxttpa
      integer mxstpa, nsitot, ntrtot
      integer mx12pa, mxvvpa, mxv2pa, maxpsq, max12
      common /auxh/ nconst(8),ncond(8),maxpa(4),maxpat,maxsin(3),
     +              maxepa(4),mxdvpa,mxddpa,nspips(8),nspipt(8),
     +              mxsspa,mxttpa,mxstpa,nsitot(8),ntrtot(8),
     +              mx12pa,mxvvpa,mxv2pa,maxpsq,max12
c
c
      integer nxblk, index, numci, mblkci, ispbig, ispsma
      common/disktl/nxblk,index(199),numci,mblkci,ispbig,ispsma
c
      logical trip
c
c... vacuum
c
      call ptczmp(c,0,conf,indcz,'reset')
c
      kc = 1
      kcc = 1
      call rdedx(vc,nvcc(1),kvcc(1),numci)
      kvc = 1
      do 10 is=1,2+nirr*2
         ncsf = ncmphn(2,is)
         nset = ncmphn(3,is)
         nst  = ncmphn(9,is)
         call vclr(c(kc),1,ncsf*nset)
         call cexpmp(c(kc),cc(kcc),1,vc(kvc),ncsf,ncsf,nst,nset)
         kc  = kc  + ncsf*nset
         kcc = kcc + nst*nset
         kvc = kvc + ncsf*nst*nset
10    continue
      call ptczmp(c,0,conf,indcz,'vacuum')
c
c...  doublet (n-1)
c
      kc = 1
      call rdedx(vc,nvcc(2),kvcc(2),numci)
      kvc = 1
      do 30 ijrr=1,nirr
         nex = norbe(ijrr)
         do 20 is=1,2+nirr*2
            ncsf = ncmph1(2,is,ijrr)
            nset = ncmph1(3,is,ijrr)
            nst  = ncmph1(9,is,ijrr)
            call vclr(c(kc),1,ncsf*nset*nex)
            call cexpmp(c(kc),cc(kcc),nex,vc(kvc),ncsf,ncsf,nst,nset)
            kc  = kc + ncsf*nset*nex
            kcc = kcc + nst*nset*nex
            kvc = kvc + ncsf*nst*nset
20       continue
30    continue
      call ptczmp(c,0,conf,indcz,'doublet')
c
c...  symmetry 1 diagonal saved in scr
c
      call rdedx(vc,nvcc(4),kvcc(4),numci)
      kvc = 1
      ksc = 1
      kdc = 1
      do 40 is=1,2+nirr*2
         ncsf = ncmphd(2,is)
         nset = ncmphd(3,is)
         nst  = ncmphd(9,is)
         call vclr(scr(ksc),1,ncsf*nset*next)
         call cexpmp(scr(ksc),dc(kdc),next,vc(kvc),ncsf,ncsf,nst,nset)
         ksc = ksc + ncsf*nset*next
         kdc = kdc + nst*nset*next
         kvc = kvc + ncsf*nst*nset
40    continue
c
c...  general (n-2) (symmetry 1 needs special treatment)
c
      trip = .false.
      kcctr = kcc
      call rdedx(vc,nvcc(3),kvcc(3),numci)
      kvc = 1
c
100   continue
c

      do 60 ijrr=1,nirr
         kc = 1
         nex = norbtr(ijrr)
         do 50 is=1,2+nirr*2
            ncsf = ncmph2(2,is,ijrr)
            nset = ncmph2(3,is,ijrr)
            nst  = ncmph2(9,is,ijrr)
            nsin = ncmph2(4,is,ijrr)
c..    expand
            if (.not.trip) then
               if (nsin.gt.0) then
                  call vclr(c(kc),1,nset*nsin*nex)
                  call cexpmp(c(kc),cc(kcc),nex,vc(kvc),
     1                        ncsf,nsin,nst,nset)
                  kc  = kc + nsin*nset*nex
               end if
            else
               ntri = ncmph2(5,is,ijrr)
               if (ntri.gt.0) then
                  call vclr(c(kc),1,nset*ntri*nex)
                  call cexpmp(c(kc),cc(kcc),nex,vc(kvc+nsin),
     1                        ncsf,ntri,nst,nset)
                  kc = kc + ntri*nset*nex
               end if
            end if
            kcc  = kcc + nst*nset*nex
            kvc  = kvc + ncsf*nst*nset
50       continue
         if (ijrr.eq.1.and..not.trip.and.nspips(1).ne.0) then
            call sin1mp2(c,nspips(1),'exp')
            call aagp(c,nspips(1),scr,next,'put')
         end if
c
         if (.not.trip) then
            call ptczmp(c,ijrr,conf,indcz,'symsin')
         else
            call ptczmp(c,ijrr,conf,indcz,'symtri')
         end if
c
60    continue
c
      if (.not.trip) then
         trip = .true.
         kvc = 1
         kcc = kcctr
         go to 100
      end if
c
      return
      end
      subroutine cexpmp(c,cc,nex,vcc,nvcc,ncsf,nst,nset)
c
c...  expand   vectors cc using vcc
c...  yielding vectors c  (external dimension nex)
c...  (singlet/triplet together)
c...  loop over sets here
c
c...  derived from contmp
c
c     c = c + cc * transpose(vcc)
c
      implicit real*8 (a-h,o-z), integer(i-n)
c
      dimension c(nex,ncsf,nset),cc(nex,nst,nset),vcc(nvcc,nst,nset)
c
      do 100 i=1,nset
chvd     call mxmb(c(1,1,i),1,nex,vcc(1,1,i),1,nvcc,
chvd *             cc(1,1,i),1,nex,nex,ncsf,nst)
         call mxmb(cc(1,1,i),1,nex,vcc(1,1,i),nvcc,1,
     *             c(1,1,i),1,nex,nex,nst,ncsf)
100   continue
c
      return
      end
      subroutine aagp(c,nspips,d,next,gp)
c
c...   either extract or  insert diagonal (aa) elements
c...   from/in a (n-2) symmetry 1 singlet vector
c
      implicit real*8 (a-h,o-z)
      dimension c(*),d(next,nspips)
      character*3 gp
c
c
      integer nirr,mult,ifrep,nirr8,jfrep
      common /table/ nirr,mult(8,8),ifrep,nirr8,jfrep
c
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
c
      if (nspips.eq.0) call caserr('no symmetry 1 n-2 vector')
c
      mu = 0
      md = 0
      do 100 irr=1,nirr
         ne = norbe(irr)
         if (ne.eq.0) go to 100
         do 50 ie=1,ne
            mu = mu + ie
            md = md + 1
            if (gp.eq.'get') then
               do 10 i=1,nspips
10             d(md,i) = c((i-1)*norbsi(1)+mu)
            else if (gp.eq.'put') then
               do 20 i=1,nspips
20             c((i-1)*norbsi(1)+mu) = d(md,i)
            else
               call caserr('illegal mpgpaa call')
            end if
50       continue
100   continue
c
      return
      end
      subroutine sin1mp2(cc,nspips,op)
c
c...   either expand or compress (take out aa) elements
c...   in a (n-2) symmetry 1 singlet vector
c...   ce : expanded ; cc : compressed
c...   NEW version not requiring extra space
c...   nspips is the number of c-vectors to treat
c
      implicit real*8 (a-h,o-z)
      dimension cc(*)
      character*3 op
c
c
      integer nirr,mult,ifrep,nirr8,jfrep
      common /table/ nirr,mult(8,8),ifrep,nirr8,jfrep
c
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
c
      if (nspips.eq.0) call caserr('no symmetry 1 n-2 vector')
c
      if (op.eq.'exp') then
c...     move cc forwards => expand works like compress
         mf = (norbsi(1)-norbtr(1))*nspips
         do i=norbtr(1)*nspips,1,-1
            cc(mf+i) = cc(i)
         end do
      else if (op.eq.'com') then
         mf = 0
      else
         call caserr('incorrect sin1mp call')
      end if
c
      mt = 0
      do 200 k=1,nspips
         do 100 irr=1,nirr
            do 50 i=1,norbe(irr)
               do 40 j=1,i-1
                  mf = mf + 1
                  mt = mt + 1
                  cc(mt) = cc(mf)
40             continue
               if (op.eq.'exp') then
                  mt = mt + 1
                  cc(mt) = 0.0d0
               else
                  mf = mf + 1
               end if
50          continue
100      continue
200   continue
c
      return
      end
      subroutine sin1mp(ce,ne,cc,nc,nspips,op)
c
c...   either expand or compress (take out aa) elements
c...   in a (n-2) symmetry 1 singlet vector
c...   ce : expanded ; cc : compressed
c
      implicit real*8 (a-h,o-z)
      dimension ce(ne,nspips),cc(nc,nspips)
      character*3 op
c
c
      integer nirr,mult,ifrep,nirr8,jfrep
      common /table/ nirr,mult(8,8),ifrep,nirr8,jfrep
c
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
c
      if (nspips.eq.0) call caserr('no symmetry 1 n-2 vector')
c
      if (op.eq.'exp') then
         call vclr(ce,1,ne*nspips)
      else if (op.eq.'com') then
         call vclr(cc,1,nc*nspips)
      else
         call caserr('incorrect sin1mp call')
      end if
c
      mc = 0
      me = 0
      do 100 irr=1,nirr
         nne = norbe(irr)
         do 50 i=1,nne
            do 40 j=1,i-1
               mc = mc + 1
               me = me + 1
               if (op.eq.'exp') then
                  do 10 k=1,nspips
10                ce(me,k) = cc(mc,k)
               else
                  do 20 k=1,nspips
20                cc(mc,k) = ce(me,k)
               end if
40          continue
            me = me + 1
50       continue
100   continue
c
      return
      end
      subroutine mpdiag(conf,diag,scr)
c...
c... to build diagonal of mp-hamiltonian in diag
c... scr is scratch-space
c...
c...  instore mode only
c...
      implicit real*8 (a-h,o-z)
c
      integer conf
c
      dimension conf(5,*),diag(*),scr(*)
c
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
      common /mr_mp2con/ nsmpv,nsmp1,nsmp2,nsmp2d,ntotmp,
     *                   nmpv,nmp1,nmp2,nmp2d
c
c     common/spew/model,ic,ir
      common/spew/rmodel(5)
      equivalence (rmodel(1),model),(rmodel(2),ic),(rmodel(3),ir)
c
c
      integer nirr,mult,ifrep,nirr8,jfrep
      common /table/ nirr,mult(8,8),ifrep,nirr8,jfrep
c
      common/loco/irebas(8),multt(8),multtt(8,8)
c
      integer nint, next, ninnex, nwm1, nw, nwp1, nintwo
      integer nint21, nintp1
      common /orb/ nint,next,ninnex,nwm1,nw,nwp1,nintwo,
     +             nint21,nintp1
c
c
      integer nconst, ncond, maxpa, maxpat, maxsin, maxepa
      integer mxdvpa, mxddpa, nspips, nspipt, mxsspa, mxttpa
      integer mxstpa, nsitot, ntrtot
      integer mx12pa, mxvvpa, mxv2pa, maxpsq, max12
      common /auxh/ nconst(8),ncond(8),maxpa(4),maxpat,maxsin(3),
     +              maxepa(4),mxdvpa,mxddpa,nspips(8),nspipt(8),
     +              mxsspa,mxttpa,mxstpa,nsitot(8),ntrtot(8),
     +              mx12pa,mxvvpa,mxv2pa,maxpsq,max12
c
c
      integer nxblk, index, numci, mblkci, ispbig, ispsma
      common/disktl/nxblk,index(199),numci,mblkci,ispbig,ispsma
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
c
      ksmp1 = nsmpv
      ksmp2 = nsmpv + nsmp1
      irebas(1) = 0
      do 1 i=2,nirr
1     irebas(i) = irebas(i-1) + norbe(i-1)
c
      call brtsb(indmp(1))
111   call getsb(rmodel(1),1)
cmp
888   if (model.gt.4) model = model - 30
cmp
      goto (999,2,3,4,5,5,7,8,999,999),model
      call caserr('unknown model in mpdiag')
c...
c... vacuum / vacuum (ij/kl)
c...
2     call getsb(rmodel(2),2)
      lams = conf(3,ic)
      if(ic.ne.ir) then
        call skipsb(lams*conf(3,ir))
         goto 22
      end if
      call getsb(scr,lams*lams)
      icc = conf(1,ic) + 1
      call dcopy(lams,scr,lams+1,diag(icc),1)
22    call getsb(rmodel(1),1)
      if (model.eq.2) go to 2
      goto 888
c...
c... doublet / doublet   (ij/kl)
c...
3     call getsb(rmodel(2),2)
       lams = conf(3,ic)
       if (ic.ne.ir) then
          call skipsb(lams*conf(3,ir))
          goto 33
       end if
       call getsb(scr,lams*lams)
       icbas = ksmp1 + conf(1,ic)
       ne = norbe(conf(4,ic))
       lams1 = lams+1
       idbas = 0
       do 43 loop=1,lams
          call vfill(scr(idbas+1),diag(icbas+1),1,ne)
          idbas=idbas+lams1
          icbas=icbas+ne
43     continue
33    call getsb(rmodel(1),1)
      if (model.eq.3) go to 3
      goto 888
c...
c...  n-2 / n-2   (ij/kl)  (cf mpij)
c...
4     call getsb(rmodel(2),2)
      call upack2(conf(3,ic),lamt,lams)
      lams = max(lams,lamt)
      if (ic.ne.ir) then
         call upack2(conf(3,ir),joop,loop)
         loop = max(joop,loop)
         call skipsb(lams*loop)
         goto 84
      end if
      isi = conf(4,ic)
      call upack2(conf(1,ic),kc,jc)
      norsi=norbtr(isi)
      lams1=lams+1
      idbas = 0
      icbas = ksmp2 + jc
      call getsb(scr,lams*lams)
      do 64 loop=1,lams
         call vfill(scr(idbas+1),diag(icbas+1),1,norsi)
         idbas=idbas+lams1
         icbas=icbas+norsi
64    continue
c
84    call getsb(rmodel(1),1)
      if(model.eq.4) go to 4
5      continue
       call possb(iwnmp(2,5),iwnmp(1,5))
       goto 111
c...
c...  doublet / doublet fock(aa)
c...  fock  operator contribution
c...
7     continue
c...    read the one needed diagonal fock-operator
      call rdedx(scr,next,indmp(5),numci)
c
27     call getsb(rmodel(2),1)
       lams = conf(3,ic)
       isi = conf(4,ic)
       ne=norbe(isi)
       idbas = ksmp1 + conf(1,ic)
c... (a/a) contribution   idbas : addres of doublet vector
       do 77 loop=1,lams
          call vadd(scr(irebas(isi)+1)
     +   ,1,diag(idbas+1),1,diag(idbas+1),1,ne)
          idbas=idbas+ne
77     continue
97     call getsb(rmodel(1),1)
       if (model.eq.37) go to 27
       go to 888
c...
c...   n-2 / n-2 fock(aa)
c...
8     continue
c...    read the one needed diagonal fock-operator
      kfock = maxpat*2
      call rdedx(scr(kfock+1),next,indmp(5),numci)
c
18    continue
c     kcoup = kfock + next
      do 1008 isi=1,nirr
      do 1008 loop=1,nirr
      isf=mult(loop,isi)
      nf=norbe(isf)
      if (nf.eq.0) isf = 0
1008  multtt(loop,isi)=isf
c
c...  no code necessary for diagonal aa
c
28    call getsb(rmodel(2),1)
      call upack2(conf(3,ic),lamt,lams)
      lamss=lams*lams
      lamtt=lamt*lamt
      if (ifdiaf.ne.1) call skipsb(lamss+lamtt)
      lams = max(lams,lamt)
      lamspt = lams
      call upack2(conf(1,ic),kc,jc)
      isi = conf(4,ic)
      icbas = ksmp2+jc
      if (isi.ne.1) goto 708
c... totally symmetric external pair
      do 748 loop=1,lams
        idbas = kfock
        do 748 moop=1,nirr
         ne = norbe(moop)
         do 448 joop=2,ne
            ioop=joop-1
            call vadd(scr(idbas+1)
     +     ,1,diag(icbas+1),1,diag(icbas+1),1,ioop)
            call vsadd(diag(icbas+1),1,scr(idbas+joop),
     1                 diag(icbas+1),1,ioop)
448      icbas=icbas+ioop
748   idbas=idbas+ne
      goto 738
c... external pair not totally symmetric
708   call icopy(nirr,multtt(1,isi),1,multt,1)
      do 788 loop=1,lamspt
         idbas = kfock
         do 508 moop=1,nirr
            ne=norbe(moop)
            if (ne.eq.0) go to 508
            isf = multt(moop)
            if (isf.lt.moop) goto 508
            nf = norbe(isf)
            iee = irebas(moop) + idbas
            iff = irebas(isf) + idbas
            do 548 joop=1,nf
               call vadd(scr(iee+1)
     +        ,1,diag(icbas+1),1,diag(icbas+1),1,ne)
               call vsadd(diag(icbas+1),1,scr(iff+joop),
     1                    diag(icbas+1),1,ne)
548         icbas=icbas+ne
508      continue
 788   continue
c
738   call getsb(rmodel(1),1)
      if (model.eq.38) go to 28
      go to 888
c
c...    model 40 : has gone
c
999   continue
c
c...  now do aa-diagonals
c
      if (indmp(6).le.0) return
c
      kaa = nsmpv+nsmp1+nsmp2+1
      call aadiag(conf,diag(kaa),next,scr)
c
      return
      end
      subroutine aadiag(conf,daa,ne,cr)
c...
c...   (n-2-aa) diagonals   (for MP)  CONTRACTED
c...
      implicit real*8 (a-h,o-z)
      integer conf
      dimension conf(5,*),daa(ne,*),cr(*)
c
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
      common /mr_mp2con/ nsmpv,nsmp1,nsmp2,nsmp2d,ntotmp,
     *                   nmpv,nmp1,nmp2,nmp2d
c     common/spew/model,ic,ir,idum1,idum2
      common/spew/rmodel(5),dum(2)
      equivalence (rmodel(1),model),(rmodel(2),ic),(rmodel(3),ir)
c
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
c
      integer nref, nnst, nmin1, nmin2, mspin, nelec
      integer nstrt1, nlast1, nstrt2, nlast2, nmin12, nminv2, nminv1
      integer nlst21, nlastx
      common /ccntl/ nref,nnst,nmin1,nmin2,mspin,nelec,
     +               nstrt1,nlast1,nstrt2,nlast2,nmin12,nminv2,nminv1,
     +               nlst21,nlastx
c
c
      integer nxblk, index, numci, mblkci, ispbig, ispsma
      common/disktl/nxblk,index(199),numci,mblkci,ispbig,ispsma
c
      call brtsb(indmp(6))
      call getsb(rmodel(1),1)
      if (model.ne.42) call caserr(' wrong model for aadiag')
c
c...  ijkl part
c
5     call getsb(rmodel(2),2)
c
      if (ic.ne.ir) then
         call upack2(conf(3,ic),ncolt,ncols)
         call upack2(conf(3,ir),nrowt,nrows)
         call skipsb(ncols*nrows)
         go to 20
      end if
c
      call upack2(conf(3,ic),ncolt,ncols)
      call getsb(cr,ncols*ncols)
      icbas = conf(2,ic)
c
      ib = 1
      do 10 i=1,ncols
         call vfill(cr(ib),daa(1,icbas+i),1,ne)
         ib = ib + ncols+1
10    continue
c
20    call getsb(rmodel(1),1)
      if (model.eq.42) go to 5
c
c...   aa/aa part
c
      if (model.ne.45) call caserr('no aa after ijkl in aadiag')
c
c...     read diagonal fock-matrix
c
      call rdedx(cr,ne,indmp(5),numci)
c
c...   scale diagonals *2 (  check)
c
      if (iprimp.ge.500) Print *,' scaling aa elements'
      call dscal(ne,2.0d0,cr,1)
c
c...  add diagonals
c
      nn = nsmp2d/ne
      do 60 i=1,nn
         call daxpy(ne,1.0d0,cr,1,daa(1,i),1)
60    continue
c
100   call getsb(rmodel(1),1)
      if (model.eq.45) call caserr(' 2 model 45 ')
c
      return
      end
      subroutine mpgen(conf,cr)
c
      implicit real*8 (a-h,o-z)
c
      integer conf
c
      dimension conf(5,*),cr(*)
c
c... to compute z=fock*c where c is a trial ci vector
c... derived from genz      (jvl 1989)
c...    **may be stripped further**
c...    uses singlet and triplet-vectors both without diagonal
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
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
      common /mr_mp2con/ nsmpv,nsmp1,nsmp2,nsmp2d,ntotmp,
     *                   nmpv,nmp1,nmp2,nmp2d
c
      common/stoctl/khbill,khbz,khbc,khbb,khbzd,khbzs,khbzt,khbcd,
     *khbcs,khbct,kh,khb1,khb2,khb3,khb4,khb5,khb6,khb7,khb8
     *,khb44,khb55,khb66,khb88,khb111,khb222,khb333
c     common/spew/model,ic,ir
      common/spew/rmodel(5)
      equivalence (rmodel(1),model),(rmodel(2),ic),(rmodel(3),ir)
c
      common/moco/ns,nt,lamsc,lamtc,lamsr,lamtr,lamss,lamtt,
     *khcsc,khzsc,khctc,khztc,isc,khcsr,khzsr,khctr,khztr,isr,
     *norsir,nortrr,lamstc
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
c
      integer iro, iky, niky, m1e, m2e, m1ia, m2coul
      integer m2exch, m2fock, m2psup, m2qsup
      integer mbase, nbase, npairi, norbi, mbasi, mnbasi
      integer iap, mbfock, mbexch, mbcoul, mbpsup, mbqsup
      integer norbas, m2exab, m2em1e
      common /helpr/ iro(nd200),iky(nd200),niky(8),m1e,m2e,m1ia,m2coul,
     +   m2exch,m2fock,m2psup,m2qsup,
     +   mbase(288),nbase(8),npairi(36),norbi(8),mbasi(8),mnbasi(8),
     +   iap(62),mbfock(62),mbexch(36),mbcoul(36),mbpsup(36),mbqsup(36),
     +   norbas(nd200),m2exab,m2em1e
c
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
c
      dimension norbev(8)
      data norbev/8*1/
c
      call brtsb(indmp(1))
      call getsb(rmodel(1),1)
c
c-----------------out of store version----------------------------------
c      call caserr(' out store disabled')
c..
c..   check old version for the general idea
c
c--------------------in store version-----------------------------------
887   kh=khbb
      call vclr(cr(khbz+1),1,ntotmp)
888   if (model.gt.4) model = model - 30
       if (model.gt.10) call caserr(' model out of bounds ')
      nodel=model
      if (odebug(40)) then
c      debug print
       write(6,*) 'mpgen: conf, nodel', nodel
       call confprt(conf)
      endif
c
      goto (9999,2,3,4,5,6,7,8,9,10),model
c...
c...      vacuum / vacuum (ij/kl) interactions
c...      uses khbb : matrix (maxpat**2)
c...           khbc/khbz
c...
2     call mpij(conf,cr(khbc+1),cr(khbz+1),cr(khbb+1),norbev,nodel)
c
c...  reposition
c
4747  if (nodel-1.gt.6) call caserr('illegal possb for instore mp')
      call possb(iwnmp(2,nodel-1),iwnmp(1,nodel-1))
      call getsb(rmodel(1),1)
      goto 888
c...
c...       doublet / doublet (ij/kl) interactions
c...       uses khbb : matrix (maxpat**2)
c...            khbcd,khbzd
c...
3      call mpij(conf,cr(khbcd+1),cr(khbzd+1),cr(khbb+1),norbe,nodel)
      goto 888
c...
c...       n-2 / n-2 (ij/kl) interactions
c...       uses khbb : matrix (maxpat**2)
c...            khbzs,khbcs
c...
4      call mpij(conf,cr(khbcs+1),cr(khbzs+1),cr(khbb+1),norbtr,nodel)
      go to 888
c...
c...      doublet / vacuum (ij/ka)-like  fock(ia) interactions
c...      uses : khbb : integrals + coupl. : m1ia+maxpat**2
c...             khbc,khbz,khbcd,khbzd
c...
5      call dvmpia(conf,cr)
       goto 888
c...
c...     n-2 / doublet (ij/ka)-like fock(ia) interactions
c...         khb111  coupling coefficients (2*maxpat**2)
c...         khb222  squared c or z matrix sub block (next(sym)**2)
c...         khb333  m or l    (2*maxpat*next(sym))
c...         integrals at ninti (symbolic) + kh=khbb
c...
6      if (nmp2.eq.0) goto 4747
c
        call sdmpia(conf,cr(khbcd+1),cr(khbzd+1),
     *              cr(khbcs+1),cr(khbzs+1),cr(khbb+1),
     *              cr(khb111+1),cr(khb222+1),cr(khb333+1))
        go to 888
c...
c...     doublet / doublet fock(ab)
c...     uses khbb : squared fock-matrix (norbsq(1))
c...          khbcd/khbzd
c...
7      if (nmp1.eq.0) goto 4747
        call ddmpab(conf,cr)
        goto 888
c...
c...    n-2 / n-2 (ia/ib)-like fock(ab)   (diagonal)
c...    uses khb7 (= khbb) fock-matrix  (squared : norbsq(1))
c..          khb3 b(singlet-singlet) (maxpat**2)
c..          khb4 b(triplet-triplet) (maxpat**2)
c..          khb8 scratch space (nubsq+mubsq+mubsq)
c...
c
8       call ssmpab(conf,cr)
c
       go to 888
c...
c...     doublet / doublet (ia/jb) interactions // model 39 not needed
c...
9       call caserr(' model 39 void ')
        goto 888
c...
c...     n-2 / n-2 (ia/jb)-type off-diagonal FOCK interactions
c...
10     call caserr(' model not needed in contracted mp2')
c
c...   end of fock*c
c
9999  continue
c
c...  do diagonal aa part
c
      if (indmp(6).le.0) return
c
      kaa = nsmpv+nsmp1+nsmp2+1
      call mpgend(conf,cr(khbc+1),cr(khbz+1),
     *            cr(khbc+kaa),cr(khbz+kaa),cr,khbb)
c
      return
      end
      subroutine mpij(conf,c,z,bb,norbe,nodel)
c...
c...  fock(ij) interactions for doublet, n-2 and vacuum
c...
c...  conf : symbolic information (global array)
c...  c,z : c and z vectors of this type
c...  bb  : space for matrix-elements maxpat**2)
c...  norbe : external dimension / symmetry
c...  nodel : current model
c...
      implicit real*8 (a-h,o-z)
c
      integer conf
c
      dimension conf(5,*),c(*),z(*),bb(*),norbe(8)
c
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
c     common/spew/model,ic,ir
      common/spew/rmodel(5)
      equivalence (rmodel(1),model),(rmodel(2),ic),(rmodel(3),ir)
c
c
      if(odebug(40)) write(6,*) 'mpij: nodel = ', nodel
999   call getsb(rmodel(2),2)
      if(odebug(40)) write(6,*)' mpij: getsb(999), ic, ir = ', ic, ir
c
       call upack2(conf(3,ic),na,nc)
       if(odebug(40)) write(6,*)' mpij: upack2, na, nc', na, nc
       nc = max(na,nc)
       call upack2(conf(3,ir),na,nr)
       if(odebug(40)) write(6,*)' mpij: upack2, na, nr', na, nr
       nr = max(na,nr)
       call getsb(bb(1),nc*nr)
      if(odebug(40)) then
c      debug print
       write(6,*) 'mpij: conf after getsb:ic,ir =', ic, ir
       call confprt(conf)
      endif
       call upack2(conf(1,ic),itmp,icbas)
       call upack2(conf(1,ir),itmp,irbas)
       if(odebug(40)) write(6,*)'mpij: icbas, irbas', icbas, irbas
       ne = norbe(conf(4,ic))
       if(odebug(40)) write(6,*)'mpij: ne = ', ne
c
       call mxmb(c(irbas+1),1,ne,bb,1,nr,z(icbas+1),1,ne,ne,nr,nc)
       if (ic.ne.ir)
     * call mxmb(c(icbas+1),1,ne,bb,nr,1,z(irbas+1),1,ne,ne,nc,nr)
c
       call getsb(rmodel(1),1)
c
      if (model.eq.nodel) goto 999
c
      return
      end
      subroutine dvmpia(conf,cr)
c
      implicit real*8 (a-h,o-z)
c
      integer conf
c
      dimension conf(5,*),cr(*)
c
c...  doublet/vacuum mp2 fock(ia) contributions
c...  derived from dvijka
c...  process ia integrals as ijka ints
c...  ninti as usual gives addressing in ia array
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
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
c
      common/loco/integr,idum,scr(201)
c     common/spew/model,ic,ir,ninti,icfock
      common/spew/rmodel(5)
      equivalence (rmodel(1),model),(rmodel(2),ic),(rmodel(3),ir)
      equivalence (rmodel(4),ninti)
c
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
c
      integer iro, iky, niky, m1e, m2e, m1ia, m2coul
      integer m2exch, m2fock, m2psup, m2qsup
      integer mbase, nbase, npairi, norbi, mbasi, mnbasi
      integer iap, mbfock, mbexch, mbcoul, mbpsup, mbqsup
      integer norbas, m2exab, m2em1e
      common /helpr/ iro(nd200),iky(nd200),niky(8),m1e,m2e,m1ia,m2coul,
     +   m2exch,m2fock,m2psup,m2qsup,
     +   mbase(288),nbase(8),npairi(36),norbi(8),mbasi(8),mnbasi(8),
     +   iap(62),mbfock(62),mbexch(36),mbcoul(36),mbpsup(36),mbqsup(36),
     +   norbas(nd200),m2exab,m2em1e
c
      common/stoctl/khbill,khbz,khbc,khbb,khbzd,khbzs,khbzt,khbcd,
     *khbcs,khbct,kh
c
      integer nxblk, index, numci, mblkci, ispbig, ispsma
      common/disktl/nxblk,index(199),numci,mblkci,ispbig,ispsma
       icl=0
       khbx = khbb + m1ia
c
c... khbb    (fock-i/a) integrals
c... khbx     coupling coefficients
c... read (fock-i/a) integrals
      call rdedx(cr(khbb+1),m1ia,indmp(2),numci)
55    call getsb(rmodel(2),3)
      if(ic.eq.icl)goto 500
      icl=ic
      itop=conf(1,ic)
      khzdub=khbzd+itop
      khcdub=khbcd+itop
      ne=norbe(conf(4,ic))
      lamd=conf(3,ic)
c
c...   the fock elements are treated as ijka and not fock-matrices
c
500   continue
      lamv=conf(3,ir)
      itop=conf(1,ir)
      khi=ninti+kh
      call getsb(cr(khbx+1),lamv*lamd)
c... construct l
      call mxmd(cr(khbx+1),lamv,1,cr(khbc+itop+1),1,1,scr(1),1,1,
     *lamd,lamv,1)
c... update z doublet
      call mxmb(cr(khi+1),1,1,scr(1),1,1,cr(khzdub+1),1,ne,
     *ne,1,lamd)
c... construct m
      call mxmd(cr(khcdub+1),ne,1,cr(khi+1),1,1,scr(1),1,1,lamd,ne,1)
c... update z vacuum
      call mxmb(cr(khbx+1),1,lamv,scr(1),1,1,cr(khbz+itop+1),1,1,
     *lamv,lamd,1)
      call getsb(rmodel(1),1)
c
      if(model.eq.35) goto 55
c
      return
      end
      subroutine sdmpia(conf,cd,zd,cs,zs,bb,coup,cz2,cml)
c
c...  n-2 / doublet fock(ia) contributions
c...  derived from sdijka
c...  rewritten out of frustration march 1990 (jvl)
c...
c...  lamr=1 special case removed
c...  standard f77 do-loops assumed
c
c... khbcd+1 doublet (n-1) c-vectors             (cd)
c... khbzd+1 doublet (n-1) z-vectors             (cd)
c... khbcs+1 n-2 c-vectors                       (cs)
c... khbzs+1 n-2 z-vectors                       (zs)
c... khbb+1  fock (ia) matrix                    (bb)
c... khb111+1  coupling coefficients            (coup)
c... khb222+1  squared c or z matrix sub block  (cz2)
c... khb333+1  m or l                           (cml)
c
      implicit real*8 (a-h,o-z)
c
      integer conf
c
      dimension conf(5,*)
      dimension cd(*),zd(*),cs(*),zs(*)
      dimension bb(*),coup(*),cz2(*),cml(*)
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
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
c
c
      integer nirr,mult,ifrep,nirr8,jfrep
      common /table/ nirr,mult(8,8),ifrep,nirr8,jfrep
c
c
      integer iro, iky, niky, m1e, m2e, m1ia, m2coul
      integer m2exch, m2fock, m2psup, m2qsup
      integer mbase, nbase, npairi, norbi, mbasi, mnbasi
      integer iap, mbfock, mbexch, mbcoul, mbpsup, mbqsup
      integer norbas, m2exab, m2em1e
      common /helpr/ iro(nd200),iky(nd200),niky(8),m1e,m2e,m1ia,m2coul,
     +   m2exch,m2fock,m2psup,m2qsup,
     +   mbase(288),nbase(8),npairi(36),norbi(8),mbasi(8),mnbasi(8),
     +   iap(62),mbfock(62),mbexch(36),mbcoul(36),mbpsup(36),mbqsup(36),
     +   norbas(nd200),m2exab,m2em1e
c
c
      integer nxblk, index, numci, mblkci, ispbig, ispsma
      common/disktl/nxblk,index(199),numci,mblkci,ispbig,ispsma
      common/loco/integr,nwww,scr(nd200),neii(8),
     *mcolcz(8),mrowcz(8),ibassc(8),ibassz(8)
c     common/spew/model,ic,ir,ninti,icfock
      common/spew/rmodel(5)
      equivalence (rmodel(1),model),(rmodel(2),ic),(rmodel(3),ir)
      equivalence (rmodel(4),ninti),(rmodel(5),icfock)
c
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
c
      icl = 0
c
c...     always read fock-matrices
c
        call rdedx(bb,m1ia,indmp(2),numci)
c
10    call getsb(rmodel(2),3)
      if (ic.eq.icl) go to 30
      icl = ic
c
      call upack2(conf(1,ic),kc,jc)
      isc=conf(4,ic)
      ns=norbtr(isc)
      nt=norbtr(isc)
      call upack2(conf(3,ic),lamtc,lamsc)
c
c...  set up symmetry info
c
      do 20 isr=1,nirr
         isi=mult(isr,isc)
         nei=norbe(isi)
         neii(isr)=nei
         if (isi.lt.isr) then
            iss=isi
            mcolcz(isr)=nei
            mrowcz(isr)=1
         else
            iss=isr
            mcolcz(isr)=1
            mrowcz(isr)=norbe(isr)
         end if
         ibassc(isr)=ibastr(iss,isc)+jc
         ibassz(isr)=ibastr(iss,isc)+jc
20    continue
c
c...   sort out doublet
c
30    continue
      lamr=conf(3,ir)
      isr=conf(4,ir)
      itop=conf(1,ir)
      ner=norbe(isr)
      if (isc.eq.1.and.ner.eq.1) call caserr('symm 1 + ne=1 sdmpia')
      lamlam=lamsc+lamtc
c
      call getsb(coup,lamlam*lamr)
c
      nei=neii(isr)
      mcol=mcolcz(isr)
      mrow=mrowcz(isr)
c
      call mxmd(cd(itop+1),1,ner,coup,1,lamr,
     *          cml,1,ner,ner,lamr,lamlam)
      icbas=ibassc(isr)
      icbat=ibassc(isr)
      imbas= 1
c
c...  ==================== update z(n-2) =======================
c
      if (isc.eq.1) then
c... z(n-2) totally symmetric
         do 9003 lam=1,lamsc
            call mxmd(cml(imbas),1,1,bb(ninti+1),1,1,
     *                cz2,1,ner,ner,1,nei)
            imbas=imbas+ner
            call symmpa(zs(icbas+1),cz2,ner)
9003     icbas=icbas+ns
         do 9004 lam=1,lamtc
            call mxmd(cml(imbas),1,1,bb(ninti+1),1,1,
     *                cz2,1,ner,ner,1,nei)
            imbas=imbas+ner
            call anti1a(zs(icbat+1),cz2,ner)
9004     icbat=icbat+nt
      else
c
c... z(n-2) not totally symmetric
c
         do 9005 lam=1,lamsc
            if(ner.lt.nei) then
               call mxmb(bb(ninti+1),1,1,cml(imbas),1,1,
     *                   zs(icbas+1),mrow,mcol,nei,1,ner)
            else
               call mxmb(cml(imbas),1,1,bb(ninti+1),1,1,
     *                   zs(icbas+1),mcol,mrow,ner,1,nei)
            end if
            imbas=imbas+ner
            icbas=icbas+ns
9005     continue
         do 9006 lam=1,lamtc
            if (ner.lt.nei) then
               call mxmb(bb(ninti+1),1,1,cml(imbas),1,1,
     *                   zs(icbat+1),mrow,mcol,nei,1,ner)
            else
               call mxmb(cml(imbas),1,1,bb(ninti+1),1,1,
     *                   zs(icbat+1),mcol,mrow,ner,1,nei)
            end if
            imbas=imbas+ner
            icbat=icbat+nt
9006     continue
      end if
c
c... ================== update z(doublet) ======================
c
      icbas=ibassc(isr)
      icbat=ibassc(isr)
      imbas= 1
      if (isc.eq.1) then
c
c... c(n-2) totally symmetric
c
         do 7003 lam=1,lamsc
            call squamp(cz2,cs(icbas+1),ner,ner)
            icbas=icbas+ns
            call mxmd(cz2,1,ner,bb(ninti+1),1,1,
     *                cml(imbas),1,1,ner,nei,1)
            imbas=imbas+ner
7003     continue
         do 7004 lam=1,lamtc
            call sqtrip(cz2,cs(icbat+1),ner)
            icbat=icbat+nt
            call mxmd(cz2,1,ner,bb(ninti+1),1,1,
     *                cml(imbas),1,1,ner,nei,1)
            imbas=imbas+ner
7004     continue
      else
c
c...  c(n-2) not totally symmetric
c
         do 7005 lam=1,lamsc
            call mxmd(cs(icbas+1),mcol,mrow,bb(ninti+1),1,1,
     *                cml(imbas),1,1,ner,nei,1)
            icbas=icbas+ns
            imbas=imbas+ner
7005     continue
         do 7006 lam=1,lamtc
            call mxmd(cs(icbat+1),mcol,mrow,bb(ninti+1),1,1,
     *                cml(imbas),1,1,ner,nei,1)
            icbat=icbat+nt
            imbas=imbas+ner
7006     continue
      end if
c
      call mxmb(cml,1,ner,coup,lamr,1,
     *          zd(itop+1),1,ner,ner,lamlam,lamr)
c
       call getsb(rmodel(1),1)
       if (model.eq.36) go to 10
c
      return
      end
      subroutine ddmpab(conf,cr)
c
c... doublet / doublet fock(ab) interactions
c...   simple diagonal fock-operator
c...  derived from ddijab
c
      implicit real*8 (a-h,o-z)
c
      integer conf
c
      dimension conf(5,*),cr(*)
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
      integer nirr,mult,ifrep,nirr8,jfrep
      common /table/ nirr,mult(8,8),ifrep,nirr8,jfrep
c
c
      integer nint, next, ninnex, nwm1, nw, nwp1, nintwo
      integer nint21, nintp1
      common /orb/ nint,next,ninnex,nwm1,nw,nwp1,nintwo,
     +             nint21,nintp1
c
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
c
      integer nconst, ncond, maxpa, maxpat, maxsin, maxepa
      integer mxdvpa, mxddpa, nspips, nspipt, mxsspa, mxttpa
      integer mxstpa, nsitot, ntrtot
      integer mx12pa, mxvvpa, mxv2pa, maxpsq, max12
      common /auxh/ nconst(8),ncond(8),maxpa(4),maxpat,maxsin(3),
     +              maxepa(4),mxdvpa,mxddpa,nspips(8),nspipt(8),
     +              mxsspa,mxttpa,mxstpa,nsitot(8),ntrtot(8),
     +              mx12pa,mxvvpa,mxv2pa,maxpsq,max12
c
c
      integer nxblk, index, numci, mblkci, ispbig, ispsma
      common/disktl/nxblk,index(199),numci,mblkci,ispbig,ispsma
      common/stoctl/khbill,khbz,khbc,khbb,khbzd,khbzs,khbzt,khbcd,
     *khbcs,khbct,kh
c
      integer iro, iky, niky, m1e, m2e, m1ia, m2coul
      integer m2exch, m2fock, m2psup, m2qsup
      integer mbase, nbase, npairi, norbi, mbasi, mnbasi
      integer iap, mbfock, mbexch, mbcoul, mbpsup, mbqsup
      integer norbas, m2exab, m2em1e
      common /helpr/ iro(nd200),iky(nd200),niky(8),m1e,m2e,m1ia,m2coul,
     +   m2exch,m2fock,m2psup,m2qsup,
     +   mbase(288),nbase(8),npairi(36),norbi(8),mbasi(8),mnbasi(8),
     +   iap(62),mbfock(62),mbexch(36),mbcoul(36),mbpsup(36),mbqsup(36),
     +   norbas(nd200),m2exab,m2em1e
c
c
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
c
c     common/spew/model,ic,ir,kkkbas
      common/spew/rmodel(5)
      equivalence (rmodel(1),model),(rmodel(2),ic),(rmodel(3),ir)
      equivalence (rmodel(4),kkkbas)
c
      dimension ibe(8)
c
c..   since only one symmetric operator
c..   use ibasssq for addressing for squared operator
c..   use symbolic only as driver (coupling is unit-matrix)
c...  if the fock-operator is diagonal in the orbitals (ifdiaf=1)
c...  fock-matrix is compressed and special simpler code is used
c
c... khbb   fock (ab)  integrals (squared)
c
      if (ifdiaf.ne.1) then
         call rdedx(cr(khbb+1),norbsq(1),indmp(4),numci)
      else
c
c...  if external fock-matrix is really diagonal compress and
c...  set up addressing in ibe
c
         call rdedx(cr(khbb+1),next,indmp(5),numci)
         ibe(1) = 0
         do 10 k=1,nirr-1
10       ibe(k+1) = ibe(k) + norbe(k)
      end if
c
55    call getsb(rmodel(2),1)
       itop=conf(1,ic)
       khzduc=khbzd+itop
       khcduc=khbcd+itop
       lamc=conf(3,ic)
       isc=conf(4,ic)
       nec=norbe(isc)
       if (ifdiaf.eq.1) then
c...    diagonal F
          joke =ibe(isc) + khbb
          do 70 j=1,lamc
             do 60 i=1,nec
60           cr(khzduc+i) = cr(khzduc+i) + cr(khcduc+i)*cr(joke+i)
             khcduc = khcduc + nec
70       khzduc = khzduc + nec
       else
c...    off-diagonal F
          joke=  ibassq(isc,1) + khbb
c... operator is totally symmetric
          call mxmb(cr(joke+1),1,nec,cr(khcduc+1),1,nec,
     *              cr(khzduc+1),1,nec,nec,nec,lamc)
       end if
c
      call getsb(rmodel(1),1)
      if (model.eq.37)goto 55
c
      return
      end
      subroutine ssmpab(conf,cr)
c
c... n-2 / n-2 fock(ab) interactions
c... derived from ssiaib/ssfock
c... only singlet/singlet and triplet/triplet interactions
c... this will work in store only !!!
c
      implicit real*8 (a-h,o-z)
c
      integer  conf
c
      dimension conf(5,*),cr(*)
c
c... khb7 (= khbb) fock-matrix  (squared)
c... khb3 b(singlet-singlet)
c... khb4 b(triplet-triplet)
c... khb8 scratch space
c
c... if the external fock-matrix is diagonal (ifdiaf=1) :
c...   khb8 : storage for fock-diagonals (next)
c...   khb7 : triangle fock-matrix // external diagonal
c
c
      integer nint, next, ninnex, nwm1, nw, nwp1, nintwo
      integer nint21, nintp1
      common /orb/ nint,next,ninnex,nwm1,nw,nwp1,nintwo,
     +             nint21,nintp1
c
      common/moco/ns,nt,lamsc,lamtc,lamsr,lamtr,lamss,lamtt,
     *khcsc,khzsc,khctc,khztc,isc,khcsr,khzsr,khctr,khztr,isr,
     *norsir,nortrr
      common/stoctl/khbill,khbz,khbc,khbb,khbzd,khbzs,khbzt,khbcd,
     *khbcs,khbct,kh,khb1,khb2,khb3,khb4,khb5,khb6,khb7,khb8
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
c
      integer nirr,mult,ifrep,nirr8,jfrep
      common /table/ nirr,mult(8,8),ifrep,nirr8,jfrep
c
      common/expanc/ibasso(8),mcolo(8),mrowo(8)
      common/loco/irebas(8),multtt(8)
c
      integer nxblk, index, numci, mblkci, ispbig, ispsma
      common/disktl/nxblk,index(199),numci,mblkci,ispbig,ispsma
c     common/spew/model,ic,ir,kkkbas
      common/spew/rmodel(4)
      equivalence (rmodel(1),model),(rmodel(2),ic),(rmodel(3),ir)
      equivalence (rmodel(4),kkkbas)
c
c
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
c
      if (ifdiaf.eq.1) go to 100
c
c...  always read fock-operators (ab)
c...   non-diagonal FOCK
c
      call rdedx(cr(khb7+1),norbsq(1),indmp(4),numci)
c
c...   setup ibasso,mcolo,mrowo for multt
c
      do 1 loop=1,nirr
         nec=norbe(loop)
         if (nec.eq.0) go to 1
          ibasso(loop)=ibassq(loop,1)
          mcolo(loop)=1
          mrowo(loop)=nec
1     continue
c
10    call getsb(rmodel(2),1)
      call upack2(conf(1,ic),kc,jc)
      khzsc=khbzs+jc
      khcsc=khbcs+jc
      call upack2(conf(3,ic),lamtc,lamsc)
      isc=conf(4,ic)
      ns=norbtr(isc)
c
       call getsb(cr(khb3+1),lamsc*lamsc)
       call getsb(cr(khb4+1),lamtc*lamtc)
c
c... z(c) singlet (special routine simpler call)
      if (lamsc.gt.0)
     *call mulssp(cr(khb7+1),cr(khcsc+1),cr(khzsc+1),lamsc,ns,isc,
     *            cr(khb3+1),cr(khb8+1))
c... z(c) triplet
      if (lamtc.gt.0)
     *call multt(cr(khb7+1),cr(khcsc+1),lamtc,ns,isc,
     *     cr(khzsc+1),lamtc,ns,isc,cr(khb4+1),1,lamtc,cr(khb8+1))
c
      call getsb(rmodel(1),1)
      if (model.eq.38) go to 10
c
      return
c
100   continue
c...  always read fock-operators (ab)
c...   diagonal FOCK
c
      call rdedx(cr(khb8+1),next,indmp(5),numci)
      iscc = 0
      irebas(1) = 0
      do 101 i=2,nirr
101   irebas(i) = irebas(i-1) + norbe(i-1)
c
105   call getsb(rmodel(2),1)
      isc=conf(4,ic)
c
      if (iscc.ne.isc) then
c...        set up diagonal
          idbas = khb8
          icbas = khb7
         if (isc.eq.1) then
c...          totally symmetric
            do 120 moop=1,nirr
              do 110 joop=2,norbe(moop)
                call dcopy(joop-1,cr(idbas+1),1,cr(icbas+1),1)
                call vsadd(cr(icbas+1),1,cr(idbas+joop),
     1                     cr(icbas+1),1,joop-1)
110           icbas = icbas+joop-1
120         idbas = idbas + norbe(moop)
         else
c... external pair not totally symmetric
            do 130 loop=1,nirr
               isf=mult(loop,isc)
               nf=norbe(isf)
               if (nf.eq.0) isf = 0
130         multtt(loop)=isf
            do 150 moop=1,nirr
               ne=norbe(moop)
               if (ne.eq.0) go to 150
               isf = multtt(moop)
               if (isf.lt.moop) goto 150
               nf = norbe(isf)
               iee = irebas(moop) + idbas
               iff = irebas(isf) + idbas
               do 140 joop=1,nf
                  call dcopy(ne,cr(iee+1),1,cr(icbas+1),1)
                  call vsadd(cr(icbas+1),1,cr(iff+joop),
     1                       cr(icbas+1),1,ne)
140            icbas=icbas+ne
150         continue
         end if
      end if
c
      call upack2(conf(1,ic),kc,jc)
      khzsc=khbzs+jc
      khcsc=khbcs+jc
      call upack2(conf(3,ic),lamtc,lamsc)
      ns=norbtr(isc)
c
c...     diagonal fock-matrix / no symbolic needed
c
c...  multiply
c
      lamsc = max(lamsc,lamtc)
      do 170 j=1,lamsc
         do 160 i=1,ns
160      cr(khzsc+i) = cr(khzsc+i) + cr(khcsc+i)*cr(khb7+i)
         khcsc = khcsc + ns
170   khzsc = khzsc + ns
c
      call getsb(rmodel(1),1)
      if (model.eq.38) go to 105
c
      return
      end
      subroutine symmpa(r,a,n)
c...   symmetrize and add, ignoring diagonal
      implicit real*8 (a-h,o-z)
      dimension r(*),a(*)
      common/scra/temp(1)
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
      integer iro, iky, niky, m1e, m2e, m1ia, m2coul
      integer m2exch, m2fock, m2psup, m2qsup
      integer mbase, nbase, npairi, norbi, mbasi, mnbasi
      integer iap, mbfock, mbexch, mbcoul, mbpsup, mbqsup
      integer norbas, m2exab, m2em1e
      common /helpr/ iro(nd200),iky(nd200),niky(8),m1e,m2e,m1ia,m2coul,
     +   m2exch,m2fock,m2psup,m2qsup,
     +   mbase(288),nbase(8),npairi(36),norbi(8),mbasi(8),mnbasi(8),
     +   iap(62),mbfock(62),mbexch(36),mbcoul(36),mbpsup(36),mbqsup(36),
     +   norbas(nd200),m2exab,m2em1e
c
      call symmp(temp,a,n)
      call vadd(temp,1,r,1,r,1,iky(n))
      return
      end
      subroutine symmp(r,a,n)
c...   symmetrize ignoring diagonal
      implicit real*8 (a-h,o-z)
      dimension r(*),a(n,*)
      m=1
      do 1 i=1,n
      do 1 j=1,i-1
      r(m)=a(i,j)+a(j,i)
1     m=m+1
      return
      end
      subroutine squamp(r,a,mrowr,n)
c... convert lower triangle a to symmetric square r
      implicit real*8 (a-h,o-z)
      dimension r(mrowr,n),a(*)
      k = 1
      do 40 i=1,n
cdir$ ivdep
      do 30 j=1,i-1
         r(i,j) = a(k)
         r(j,i) = a(k)
         k=k+1
30    continue
40    r(i,i) = 0.0d0
      return
      end
      subroutine mulssp(oper,c,z,lamc,norsic,isc,b,r)
c
c...   mulss for mp2  singlets without diagonal
c
      implicit real*8 (a-h,o-z)
      dimension b(*),oper(*),c(*),z(*),r(*)
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
      integer nirr,mult,ifrep,nirr8,jfrep
      common /table/ nirr,mult(8,8),ifrep,nirr8,jfrep
c
c
      integer iro, iky, niky, m1e, m2e, m1ia, m2coul
      integer m2exch, m2fock, m2psup, m2qsup
      integer mbase, nbase, npairi, norbi, mbasi, mnbasi
      integer iap, mbfock, mbexch, mbcoul, mbpsup, mbqsup
      integer norbas, m2exab, m2em1e
      common /helpr/ iro(nd200),iky(nd200),niky(8),m1e,m2e,m1ia,m2coul,
     +   m2exch,m2fock,m2psup,m2qsup,
     +   mbase(288),nbase(8),npairi(36),norbi(8),mbasi(8),mnbasi(8),
     +   iap(62),mbfock(62),mbexch(36),mbcoul(36),mbpsup(36),mbqsup(36),
     +   norbas(nd200),m2exab,m2em1e
c
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
c
c... z(singlet)  -  c(singlet)   fock(ab) interactions
      ibbb=1
      izzz=1
      do 100 lam=1,lamc
         call mxmd(c,1,norsic,b(ibbb),1,lamc,r,1,1,norsic,lamc,1)
         do 200 iooz=1,nirr
            niz=norbe(iooz)
            if (niz.eq.0) go to 200
            jooz=mult(iooz,isc)
            njz=norbe(jooz)
            if (njz.eq.0) go to 200
            if (njz.eq.1.and.isc.eq.1) go to 200
            jooo=mult(jooz,isc)
c...   z(iooz,jooz) = oper(iooz,jooo) * c(jooo,jooz)
            njo=norbe(jooo)
            if (njo.eq.0) go to 200
            lbas=norsic
            if (jooz-jooo) 211,210,212
c... c vector totally symmetric
210          call squamp(r(lbas+1),r(ibassi(jooo,1)+1),njo,njo)
             kbas=lbas
             lbas=lbas+njo*njo
             goto 213
c... c transposed
211          kbas=ibassi(jooz,isc)
             mmmc=njz
             nnnc=1
             goto 214
c... c not transposed
212          kbas=ibassi(jooo,isc)
213          nnnc=njo
             mmmc=1
214         if(jooz-iooz)301,300,302
c...  z vector totally symmetric
300          call mxmd(oper(ibassq(iooz,1)+1),1,niz,
     *                 r(kbas+1),mmmc,nnnc,r(lbas+1),1,niz,niz,njo,njz)
             call symmpa(z(izzz+ibassi(iooz,1)),r(lbas+1),niz)
             goto 200
c... z vector transposed
301          nnnz=1
             mmmz=njz
             jump=ibassi(jooz,isc)+izzz
             goto 888
c... z vector not transposed
302          mmmz=1
             nnnz=niz
             jump=ibassi(iooz,isc)+izzz
888          if (niz.lt.njz) then
c... evaluate z (transpose)
               call mxmb(r(kbas+1),nnnc,mmmc,oper(ibassq(iooz,1)+1),
     *                   niz,1,z(jump),nnnz,mmmz,
     *                   njz,njo,niz)
             else
c... evaluate z
               call mxmb(oper(ibassq(iooz,1)+1),1,niz,
     *                   r(kbas+1),mmmc,nnnc,z(jump),mmmz,nnnz,
     *                   niz,njo,njz)
             end if
200      continue
         ibbb=ibbb+lamc
100      izzz=izzz+norsic
c
      return
      end
      subroutine mpgend(conf,c,z,caa,zaa,cr,kcr)
c
      implicit real*8 (a-h,o-z)
c
      integer  conf
c
      dimension conf(5,*),c(*),z(*),caa(*),zaa(*),cr(*)
c
c...   control routine for aa contributions to (f-e0)c
c...   c/z : normal c/z vectors (only till sym1-n-2 used)
c...   caa/zaa : aa parts of n-2 vectors
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
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
      common /mr_mp2con/ nsmpv,nsmp1,nsmp2,nsmp2d,ntotmp,
     *                   nmpv,nmp1,nmp2,nmp2d
c
c
      integer nxblk, index, numci, mblkci, ispbig, ispsma
      common/disktl/nxblk,index(199),numci,mblkci,ispbig,ispsma
      common/symcon/root2,root2i
      common/stoctl/khbill,khbz,khbc,khbb,khbzd,khbzs,khbzt,khbcd,
     *khbcs,khbct,kh,khb1,khb2,khb3,khb4,khb5,khb6,khb7,khb8
     *,khb44,khb55,khb66,khb88
c
      integer lcall,  mrstor
      logical modev, moded, modest, instor, loijka, msort, loijab
      logical loiajb, lvar
      common/diactl/modev,moded,modest,instor,
     +              loijka,msort,loijab,loiajb,lcall,mrstor,lvar
c
c
      integer iro, iky, niky, m1e, m2e, m1ia, m2coul
      integer m2exch, m2fock, m2psup, m2qsup
      integer mbase, nbase, npairi, norbi, mbasi, mnbasi
      integer iap, mbfock, mbexch, mbcoul, mbpsup, mbqsup
      integer norbas, m2exab, m2em1e
      common /helpr/ iro(nd200),iky(nd200),niky(8),m1e,m2e,m1ia,m2coul,
     +   m2exch,m2fock,m2psup,m2qsup,
     +   mbase(288),nbase(8),npairi(36),norbi(8),mbasi(8),mnbasi(8),
     +   iap(62),mbfock(62),mbexch(36),mbcoul(36),mbpsup(36),mbqsup(36),
     +   norbas(nd200),m2exab,m2em1e
c
c
      integer nconst, ncond, maxpa, maxpat, maxsin, maxepa
      integer mxdvpa, mxddpa, nspips, nspipt, mxsspa, mxttpa
      integer mxstpa, nsitot, ntrtot
      integer mx12pa, mxvvpa, mxv2pa, maxpsq, max12
      common /auxh/ nconst(8),ncond(8),maxpa(4),maxpat,maxsin(3),
     +              maxepa(4),mxdvpa,mxddpa,nspips(8),nspipt(8),
     +              mxsspa,mxttpa,mxstpa,nsitot(8),ntrtot(8),
     +              mx12pa,mxvvpa,mxv2pa,maxpsq,max12
c
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
c
      integer nint, next, ninnex, nwm1, nw, nwp1, nintwo
      integer nint21, nintp1
      common /orb/ nint,next,ninnex,nwm1,nw,nwp1,nintwo,
     +             nint21,nintp1
c
c
      integer nuse, ngot
      common /corctl/ nuse, ngot
c
c     common/spew/model,ic,ir,ninti,icfock
      common/spew/rmodel(5)
      equivalence (rmodel(1),model),(rmodel(2),ic),(rmodel(3),ir)
      equivalence (rmodel(4),ninti),(rmodel(5),icfock)
c
c..     for mpij
      dimension nextt(2)
      nextt(1) = next
c
c...  all scratch space is in cr starting at kcr
c
      call vclr(zaa,1,nsmp2d)
c..     renorm
      call dscal(nsmp2d,root2,caa,1)
c
      kbb = kcr
      kf  = kbb + maxpat**2
      kscr = kf + max(next,m1ia,norbsq(1))
      nuse = kscr + max(maxpat*next,norbsi(1))
      call corfai('mp')
c
      call brtsb(indmp(6))
      call getsb(rmodel(1),1)
100   continue
c
      if (model.eq.42) then
c
c..   (ij)
c
         call mpij(conf,caa,zaa,cr(kbb+1),nextt,42)
c
      else if (model.eq.45) then
c
c...   aa - aa  (aa)
c
         call aampaa(caa,zaa,next,cr(kf+1))
c
      else if (model.eq.44) then
c
c...   aa - ab (aa)
c
         kcz = nsmpv+nsmp1+1
         call  aampab(conf,caa,zaa,next,c(kcz),z(kcz),
     *                cr(kbb+1),cr(kf+1),cr(kscr+1))
c
      else if (model.eq.43) then
c
c...   aa - doublet (ia)
c
        kcz = nsmpv + 1
        call  aampia(conf,caa,zaa,next,c(kcz),z(kcz),
     *               cr(kbb+1),cr(kf+1),cr(kscr+1))
c
      else if (model.eq.41) then
c
c...  end of contributions
c
         go to 200
c
      else
         call caserr('unknown model in mpgend')
      end if
c
      go to 100
c
c...  cleanup
c
200   continue
c
c..     renorm    check !!!!!
c
      call dscal(nsmp2d,root2i,zaa,1)
      call dscal(nsmp2d,root2i,caa,1)
c
      return
      end
      subroutine aampia(conf,c,z,ne,cd,zd,bb,f,scr)
c
c...  singlet doublet ijaa - ia contributions
c
      implicit real*8 (a-h,o-z)
c
      integer  conf
      dimension conf(5,*)
      dimension c(ne,*),z(ne,*),cd(*),zd(*),bb(*),f(*),scr(*)
c
c     c,z : n-2 singlet c and z-vectors (diagonal only / renorm I)
c     ne : # external orbitals
c     cd/zd : doublet c/z-vector space
c     bb : space for coupling coefficients(maxpat**2)
c     f : space for f-matrix
c     scr : scratch array for intermediates (maxpat*norbe or so)
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
      integer nirr,mult,ifrep,nirr8,jfrep
      common /table/ nirr,mult(8,8),ifrep,nirr8,jfrep
c
c     common/spew/model,ic,ir,ninti,icfock
      common/spew/rmodel(5)
      equivalence (rmodel(1),model),(rmodel(2),ic),(rmodel(3),ir)
      equivalence (rmodel(4),ninti),(rmodel(5),icfock)
c
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
c
      integer iro, iky, niky, m1e, m2e, m1ia, m2coul
      integer m2exch, m2fock, m2psup, m2qsup
      integer mbase, nbase, npairi, norbi, mbasi, mnbasi
      integer iap, mbfock, mbexch, mbcoul, mbpsup, mbqsup
      integer norbas, m2exab, m2em1e
      common /helpr/ iro(nd200),iky(nd200),niky(8),m1e,m2e,m1ia,m2coul,
     +   m2exch,m2fock,m2psup,m2qsup,
     +   mbase(288),nbase(8),npairi(36),norbi(8),mbasi(8),mnbasi(8),
     +   iap(62),mbfock(62),mbexch(36),mbcoul(36),mbpsup(36),mbqsup(36),
     +   norbas(nd200),m2exab,m2em1e
c
c
      integer nxblk, index, numci, mblkci, ispbig, ispsma
      common/disktl/nxblk,index(199),numci,mblkci,ispbig,ispsma
c
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
c
       dimension ibase(8)
c...    adressing in each (aa) diagonal / symmetry
       ibase(1) = 1
       do 1 i=2,nirr
1      ibase(i) = ibase(i-1) + norbe(i-1)
c
c...     read fock-matrices
c
       call rdedx(f,m1ia,indmp(2),numci)
c
10    call getsb(rmodel(2),3)
c
c...    ir : doublets / ic singlets
c
      lamr = conf(3,ir)
      isr = conf(4,ir)
      itop = conf(1,ir)
      ner=norbe(isr)
c
      lams = conf(3,ic)
      ict = conf(2,ic)
c
      call getsb(bb(1),lams*lamr)
c
c...    cd*f = scr
c
      do i=1,lamr
        do j=1,ner
          scr((i-1)*ner+j) = f(ninti+j)*cd(itop+(i-1)*ner+j)
        end do
      end do
c
c...   scale scr by 2 to get proper sqrt(2) for z-vector
c...   combined with the renormalisation
c...   since c has a sqrt(2) this is not needed in the other direction
c
       call dscal(ner*lamr,2.0d0,scr,1)
c
c..     scr*bb  => z(n-2)
c
      call mxmb(scr,1,ner,bb,1,lamr,z(ibase(isr),ict+1),1,ne,
     *          ner,lamr,lams)
c
c...    c*f = scr
c
      do i=1,lams
         do j=1,ner
            scr((i-1)*ner+j) = f(ninti+j)*c(ibase(isr),ict+j)
         end do
      end do
c
c...   bb*scr => zd
c
      call mxmb(scr,1,ner,bb,lamr,1,zd(itop+1),1,ner,ner,lams,lamr)
c
      call getsb(rmodel(1),1)
      if (model.eq.43) goto 10
c
      return
      end
      subroutine aampaa(c,z,ne,f)
c
c...  singlet singlet aa/aa (model 45) (diagonal only)
c...  no real symbolic needed
c
      implicit real*8 (a-h,o-z)
c
      dimension c(ne,*),z(ne,*),f(*)
c
c     c,z : n-2 singlet c and z-vectors (renorm I)
c     ne : # external orbitals
c     f : space for f-matrix
c
c     common/spew/model,ic,ir,ninti,icfock
      common/spew/rmodel(5)
      equivalence (rmodel(1),model),(rmodel(2),ic),(rmodel(3),ir)
      equivalence (rmodel(4),ninti),(rmodel(5),icfock)
c
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
c
      integer nref, nnst, nmin1, nmin2, mspin, nelec
      integer nstrt1, nlast1, nstrt2, nlast2, nmin12, nminv2, nminv1
      integer nlst21, nlastx
      common /ccntl/ nref,nnst,nmin1,nmin2,mspin,nelec,
     +               nstrt1,nlast1,nstrt2,nlast2,nmin12,nminv2,nminv1,
     +               nlst21,nlastx
c
c
      integer nxblk, index, numci, mblkci, ispbig, ispsma
      common/disktl/nxblk,index(199),numci,mblkci,ispbig,ispsma
c
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
      common /mr_mp2con/ nsmpv,nsmp1,nsmp2,nsmp2d,ntotmp,
     *                   nmpv,nmp1,nmp2,nmp2d
c
c...     read diagonal fock-matrix
c
       call rdedx(f,ne,indmp(5),numci)
cd      print *,' scaling faa by 2.0'
      call dscal(ne,2.0d0,f,1)
c
      nn = nsmp2d/ne
c
      do 20 j=1,nn
         do 20 i=1,ne
20       z(i,j) = z(i,j) + c(i,j)*f(i)
c
      call getsb(rmodel(1),1)
c
      return
      end
      subroutine aampab(conf,caa,zaa,ne,c,z,bb,f,scr)
c
      implicit real*8 (a-h,o-z)
c
      integer  conf
c
      dimension conf(5,*),caa(ne,*),zaa(ne,*),c(*),z(*)
      dimension bb(*),f(*),scr(*)
c...
c...     control n-2 / n-2   aa(ic) - ab(ir) interactions  (model 44)
c...
c     common/spew/model,ic,ir
      common/spew/rmodel(5)
      equivalence (rmodel(1),model),(rmodel(2),ic),(rmodel(3),ir)
c
c
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
c
      integer nxblk, index, numci, mblkci, ispbig, ispsma
      common/disktl/nxblk,index(199),numci,mblkci,ispbig,ispsma
c
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
c
c...  read the squared fock-matrix
c
       call rdedx(f,norbsq(1),indmp(4),numci)
c
c...  loop over interactions
c
10    call getsb(rmodel(2),2)
            lamc = conf(3,ic)
            ict = conf(2,ic)
      call upack2(conf(3,ir),lamtr,lamsr)
      call upack2(conf(1,ir),kc,jc)
      call getsb(bb(1),lamc*lamsr)
c
      call mulaap(f,caa(1,ict+1),zaa(1,ict+1),lamc,ne,
     *            c(jc+1),z(jc+1),lamsr,bb,scr)
c
      call getsb(rmodel(1),1)
      if (model.eq.44) go to 10
c
      return
      end
      subroutine mulaap(oper,cd,zd,lamd,nord,c,z,lamr,b,r)
c
c...   mulss for symmetry 1 (only)  ab(ir) - aa(ic) interactions
c
      implicit real*8 (a-h,o-z)
      dimension b(lamr,lamd),oper(*),r(*)
      dimension cd(*),zd(*),c(*),z(*)
c
c
      integer nirr,mult,ifrep,nirr8,jfrep
      common /table/ nirr,mult(8,8),ifrep,nirr8,jfrep
c
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
c
c...   cd(ic) ..  z (ir)
c
c..    cd * couplings => r
c
      call mxmd(cd,1,nord,b,lamr,1,r,1,nord,nord,lamd,lamr)
c
      irrr = 1
      lbas = lamr*nord
      jump = 1
      do 30 lam=1,lamr
         ioper = 1
         do 20 iooz=1,nirr
            niz=norbe(iooz)
            if (niz.eq.0) go to 20
c... r vector only a diagonal // oper symmetrical
            call vclr(r(lbas+1),1,niz*niz)
            it = 1
            do 10 i=1,niz
               call daxpy(niz,r(irrr),oper(ioper),1,r(lbas+it),1)
               irrr = irrr + 1
               ioper = ioper + niz
               it = it + niz
10          continue
c..   add result to z-vector
            call symmpa(z(jump),r(lbas+1),niz)
            jump = jump + niz*(niz-1)/2
20       continue
30    continue
c
c...   c(ir) ..zd (ic)
c
      irrr = 1
      lbas = lamr*nord
      jump = 1
      do 60 lam=1,lamr
         ioper = 1
         do 50 iooz=1,nirr
            niz=norbe(iooz)
            if (niz.eq.0) go to 50
c...  from c*integrals vectors in r
            call squamp(r(lbas+1),c(jump),niz,niz)
            it = 1
            do 40 i=1,niz
               r(irrr) = ddot(niz,oper(ioper),1,r(lbas+it),1)*2.0d0
               irrr = irrr + 1
               ioper = ioper + niz
               it = it + niz
40          continue
            jump = jump + niz*(niz-1)/2
50       continue
60    continue
c
c...  we have now the nord*lamr c*integrals products
c...  multiply them by the couplings and add
c
      call mxmb(r,1,nord,b,1,lamr,zd,1,nord,nord,lamr,lamd)
c
      return
      end
      subroutine mpsum(conf,c,z,epair_mp)
c
c...  sum mp2 vectors into pair-energies
c
      implicit real*8 (a-h,o-z)
c
      integer  conf
      dimension conf(5,*),c(*),z(*),epair_mp(2,9)
c
      common /mr_mp2con/ nsmpv,nsmp1,nsmp2,nsmp2d,ntotmp,
     *                   nmpv,nmp1,nmp2,nmp2d
c
      integer nint, next, ninnex, nwm1, nw, nwp1, nintwo
      integer nint21, nintp1
      common /orb/ nint,next,ninnex,nwm1,nw,nwp1,nintwo,
     +             nint21,nintp1
c
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
c
c...    vacuum
c
      icc = 2
      do 20 ic=2,nmpv
         call upack2(conf(3,ic),ncolt,ncols)
         call upack2(conf(5,ic),ii,jj)
         epair_mp(ii,jj) = epair_mp(ii,jj)
     *                + ddot(ncols,c(icc),1,z(icc),1)
20    icc = icc + ncols
c
c...    doublet
c
      do 30 ic=nmpv+1,nmpv+nmp1
         call upack2(conf(3,ic),ncolt,ncols)
         call upack2(conf(5,ic),ii,jj)
         call upack2(conf(4,ic),idum,isi)
         nn = norbe(isi)*ncols
         epair_mp(ii,jj) = epair_mp(ii,jj)
     *                + ddot(nn,c(icc),1,z(icc),1)
30    icc = icc + nn
c
c...    n-2
c
      do 40 ic=nmpv+nmp1+1,nmpv+nmp1+nmp2
         call upack2(conf(5,ic),ii,jj)
         call upack2(conf(4,ic),idum,isi)
         call upack2(conf(3,ic),lamt,lams)
         lams = max(lams,lamt)
         nn = lams*norbtr(isi)
         epair_mp(ii,jj) = epair_mp(ii,jj)
     *                + ddot(nn,c(icc),1,z(icc),1)
40    icc = icc + nn
c
c...    n-2 aa-diagonal
c
      do 50 ic=nmpv+nmp1+nmp2+1,nmpv+nmp1+nmp2+nmp2d
         call upack2(conf(5,ic),ii,jj)
         call upack2(conf(3,ic),lamt,lams)
         nn = lams*next
         epair_mp(ii,jj) = epair_mp(ii,jj)
     *                + ddot(nn,c(icc),1,z(icc),1)
50    icc = icc + nn
c
      return
      end
c=====================additional-routines============================
      subroutine cistor
      implicit real*8 (a-h,o-z), integer(i-n)
c
c...  Stores current contents of common blocks <stoctl>, <auxh>,
c...  and <sym> in a help common block <ci3>. The stored information
c...  is needed for <genz> to do an uncontracted matrix-vector
c...  multiplication.
c...  However, the contents of the common blocks are modified in the
c...  contracted MP2 calculation. Therefore the data has to be
c...  restored upon transforming from the contracted basis to the
c...  configuration state function basis.
c
c...  The information should be stored before <basmp> is called.
c
      common/stoctl/kh(26)
c
      integer nconst, ncond, maxpa, maxpat, maxsin, maxepa
      integer mxdvpa, mxddpa, nspips, nspipt, mxsspa, mxttpa
      integer mxstpa, nsitot, ntrtot
      integer mx12pa, mxvvpa, mxv2pa, maxpsq, max12
      common /auxh/ nconst(8),ncond(8),maxpa(4),maxpat,maxsin(3),
     +              maxepa(4),mxdvpa,mxddpa,nspips(8),nspipt(8),
     +              mxsspa,mxttpa,mxstpa,nsitot(8),ntrtot(8),
     +              mx12pa,mxvvpa,mxv2pa,maxpsq,max12
c
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
      integer cicblk,ci,cixpat,cibssi
      common/ci3/cicblk,ci(26),cixpat,cibssi(8,8)
c
      do 10 i = 1, 26
        ci(i) = kh(i)
10    continue
      cixpat = maxpat
      do 20 j = 1, 8
        do 30 i = 1, 8
          cibssi(i,j) = ibassi(i,j)
30      continue
20    continue
c
      end
      subroutine ciretr
      implicit real*8 (a-h,o-z), integer(i-n)
c
c...  Retrieves contents of common blocks <stoctl>, <auxh> and <sym>
c...  from a help common block <ci3>. The retrieved information
c...  is needed for <genz> to do an uncontracted matrix-vector
c...  multiplication.
c
c...  The information should be retrieved before <genz> is called.
c
      common/stoctl/kh(26)
c
      integer nconst, ncond, maxpa, maxpat, maxsin, maxepa
      integer mxdvpa, mxddpa, nspips, nspipt, mxsspa, mxttpa
      integer mxstpa, nsitot, ntrtot
      integer mx12pa, mxvvpa, mxv2pa, maxpsq, max12
      common /auxh/ nconst(8),ncond(8),maxpa(4),maxpat,maxsin(3),
     +              maxepa(4),mxdvpa,mxddpa,nspips(8),nspipt(8),
     +              mxsspa,mxttpa,mxstpa,nsitot(8),ntrtot(8),
     +              mx12pa,mxvvpa,mxv2pa,maxpsq,max12
c
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
      integer cicblk,ci,cixpat,cibssi
      common/ci3/cicblk,ci(26),cixpat,cibssi(8,8)
c
      do 10 i = 1, 26
        kh(i) = ci(i)
10    continue
      maxpat = cixpat
      do 20 j = 1, 8
        do 30 i = 1, 8
          ibassi(i,j) = cibssi(i,j)
30      continue
20    continue
c
      end
      subroutine sciblk(iblk)
      implicit real*8 (a-h,o-z), integer(i-n)
c
c...  Sets the block number at which <iconf(2,5,nconf)> may be
c...  stored. The configuration table should be stored before <basmp>
c...  is called. This table is needed for <genz> to do an
c...  uncontracted matrix-vector multiplication.
c
      integer cicblk,ci,cixpat,cibssi
      common/ci3/cicblk,ci(26),cixpat,cibssi(8,8)
c
      cicblk = iblk
c
      end
      subroutine gciblk(iblk)
      implicit real*8 (a-h,o-z), integer(i-n)
c
c...  Gets the block number at which <iconf(2,5,nconf)> may be
c...  stored. The configuration table should be stored before <basmp>
c...  is called. This table is needed for <genz> to do an
c...  uncontracted matrix-vector multiplication.
c
      integer cicblk,ci,cixpat,cibssi
      common/ci3/cicblk,ci(26),cixpat,cibssi(8,8)
c
      iblk = cicblk
c
      end
      subroutine mpstor
      implicit real*8 (a-h,o-z), integer(i-n)
c
c...  Stores current contents of common blocks <stoctl>, <auxh>,
c...  and <sym> in a help common block <mp3>. The stored information
c...  is needed for <genz> to do a MP3 calculation. However, the
c...  contents of the common blocks are modified in the contracted
c...  MP2 calculation. Therefore the data has to be restored upon
c...  transforming from the contracted basis to the configuration
c...  state function basis.
c
c...  The information should be stored before <basmp> is called.
c
      common/stoctl/kh(26)
c
      integer nconst, ncond, maxpa, maxpat, maxsin, maxepa
      integer mxdvpa, mxddpa, nspips, nspipt, mxsspa, mxttpa
      integer mxstpa, nsitot, ntrtot
      integer mx12pa, mxvvpa, mxv2pa, maxpsq, max12
      common /auxh/ nconst(8),ncond(8),maxpa(4),maxpat,maxsin(3),
     +              maxepa(4),mxdvpa,mxddpa,nspips(8),nspipt(8),
     +              mxsspa,mxttpa,mxstpa,nsitot(8),ntrtot(8),
     +              mx12pa,mxvvpa,mxv2pa,maxpsq,max12
c
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
      integer dmsk
      common/mp3/dmsk(3),mpcblk,mp(26),mpxpat,mpbssi(8,8)
c
      do 10 i = 1, 26
        mp(i) = kh(i)
10    continue
      mpxpat = maxpat
      do 20 j = 1, 8
        do 30 i = 1, 8
          mpbssi(i,j) = ibassi(i,j)
30      continue
20    continue
c
      end
      subroutine mpretr
      implicit real*8 (a-h,o-z), integer(i-n)
c
c...  Retrieves contents of common blocks <stoctl>, <auxh> and <sym>
c...  from a help common block <mp3>. The retrieved information
c...  is needed for <genz> to do a MP3 calculation.
c
c...  The information should be retrieved before <genz> is called.
c
      common/stoctl/kh(26)
c
      integer nconst, ncond, maxpa, maxpat, maxsin, maxepa
      integer mxdvpa, mxddpa, nspips, nspipt, mxsspa, mxttpa
      integer mxstpa, nsitot, ntrtot
      integer mx12pa, mxvvpa, mxv2pa, maxpsq, max12
      common /auxh/ nconst(8),ncond(8),maxpa(4),maxpat,maxsin(3),
     +              maxepa(4),mxdvpa,mxddpa,nspips(8),nspipt(8),
     +              mxsspa,mxttpa,mxstpa,nsitot(8),ntrtot(8),
     +              mx12pa,mxvvpa,mxv2pa,maxpsq,max12
c
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
      integer dmsk
      common/mp3/dmsk(3),mpcblk,mp(26),mpxpat,mpbssi(8,8)
c
      do 10 i = 1, 26
        kh(i) = mp(i)
10    continue
      maxpat = mpxpat
      do 20 j = 1, 8
        do 30 i = 1, 8
          ibassi(i,j) = mpbssi(i,j)
30      continue
20    continue
c
      end
      subroutine smpblk(iblk)
      implicit real*8 (a-h,o-z), integer(i-n)
c
c...  Sets the block number at which <iconf(2,5,nconf)> may be
c...  stored. The configuration table should be stored before <basmp>
c...  is called. This table is needed for <genz> to do a MP3
c...  calculation.
c
      integer dmsk
      common/mp3/dmsk(3),mpcblk,mp(26),mpxpat,mpbssi(8,8)
c
      mpcblk = iblk
c
      end
      subroutine gmpblk(iblk)
      implicit real*8 (a-h,o-z), integer(i-n)
c
c...  Gets the block number at which <iconf(2,5,nconf)> may be
c...  stored. The configuration table should be stored before <basmp>
c...  is called. This table is needed for <genz> to do a MP3
c...  calculation.
c
      integer dmsk
      common/mp3/dmsk(3),mpcblk,mp(26),mpxpat,mpbssi(8,8)
c
      iblk = mpcblk
c
      end
      subroutine smpmsk(rmsk)
      implicit real*8 (a-h,o-z), integer(i-n)
      integer rmsk
      dimension rmsk(3)
      integer dmsk
      common/mp3/dmsk(3),mpcblk,mp(26),mpxpat,mpbssi(8,8)
c
c...  Sets the mask for double occupied orbitals
c
      do i=1,3
       dmsk(i) = rmsk(i)
      enddo
c
      end
      subroutine gmpmsk(rmsk)
      implicit real*8 (a-h,o-z), integer(i-n)
      integer rmsk
      dimension rmsk(3)
c
c...  Sets the mask for double occupied orbitals
c
      integer dmsk
      common/mp3/dmsk(3),mpcblk,mp(26),mpxpat,mpbssi(8,8)
c
      do i=1,3
       rmsk(i) = dmsk(i)
      enddo
c
      end
c...  CONTRACTED CI AND CEPA
c...
c...  This file contains additional subroutines to run contracted
c...  multi-reference ci and cepa calculations.
c...
c...  The file contains the following subroutines:
c...
c...  subroutine cciin
c...  subroutine ccepin
c...  subroutine davidc(cr,cc)
c...  subroutine genzcc(cr,h11o,indc,indz,kcc,kc,kdc,kvc,kscr)
c...  subroutine ccvpri(conf,c,crit)
c...  subroutine cntrcc(conf,cr1,cr2,cr,ic1,ic12,sum12)
c...  subroutine hc2cep(cc,nc,h,hs,h2,h2s,ndim)
c...  subroutine ccpvds(conf,z,c)
c...  subroutine vvpri(cr,kc,kcc,kdc,kvc,kscr,ntotal,ntotmp,crit)
c...
c...  The subroutines are straight forward generalisations of
c...  their uncontracted counter parts. Except for "vvpri",
c...  more information on this subroutine is in the subroutine
c...  itself.
c...
c...  hvd, 1995
c
      subroutine cciin
c
      implicit real*8 (a-h,o-z), integer (i-n)
c
c...  input routine for Contracted CI
c...
c...  called from inpro
c...
c...   = directive format = (all on 1 line)
c...
c...   conci (print 1) (crit (absolute/relative) 1.0 6)
c...       (level 0.0 (5 0.0)) (maxcyc 50) (fignore 1. 8)
c...       (nodoc) (nodiag) (schmidt 1. 8) (trial ref/vac 1)
c...
c...   print : more output for a higher number (part. >2000)
c...   crit : stop-criterion for preconditioner-iteration
c...          (max. residue)
c...     absolute : use stop-criterion as is.
c...     relative : use stop-criterion relative to the norm of the
c...                righthand side of linear system. (only in MPn)
c...   level : levelshift1, # iterations, levelshift2
c...           after the specified # iterations shift2 is used
c...   maxcyc : maximum # cycles in precondition-iteration-process
c...   fignore : absolute smaller FOCK-matrix elements are ignored
c...   cschm   : criterion for orthogonality and linear dependency in
c...             modified Gramm-Schmidt procedure.
c...   nodoc : tells program to treat doc as voc (for debugging)
c...   nodiag : tells the program to assume the external FOCK-
c...            integral-matrix to be non-diagonal (for debugging)
c...   trial : to generate a trial vector for the contracted ci.
c...           this is NOT the vector used to construct the contracted
c...           ci space.
c...   = only the first 4 chars of a keyword are read =
c...   = usually the following directive suffices :  mp 2
c...
c...   (nodoc is transmitted to mpini etc by iwnmp(1,1) = 1000)
c...
c...
c...  various routines from DIRECT and the atmol-library are also used
c...  routines from DIRECT really changed :
c...     block-data,direct,inpro,pass
c...
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
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
      common/mr_comcci/icci,iclass(3,3,3),istcc,iselcc
c
      real*8 critcep, crycep
      integer ifcepa, icepa , lcepaul, lsinsh, itcepa, lprcep, indcep
      integer nrfcsf, mcycep, ibrst,icorr_mr,mrdoc
      common/ccepa/critcep,crycep,ifcepa,icepa,lcepaul,lsinsh,itcepa,
     1                lprcep,indcep,nrfcsf,mcycep,ibrst,icorr_mr,mrdoc
c
      logical malt
      real*8 potnuc,totnuc,e11,eigen,h11,z11,h22,s22,h12,s12,
     *     tester,eshift,thresh,aaa,bbb,scream,gprin,z12,z22,
     *     vdstpr,qester
      integer maxcyc,moddia,mprin,icyc
      common/erg/potnuc,totnuc,e11,eigen(nd200),h11
     *,z11,h22,s22,h12,s12,tester,eshift,thresh,aaa,bbb,scream
     *,gprin,z12,z22,vdstpr(3),qester,maxcyc,moddia,mprin,malt,icyc
c
      common/diisdc/ diiscr,ntdiis
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
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
      character*8 test
      character*1 cnum(9)
      parameter (nmodel=14)
      character*43 descr(nmodel)
      data descr /"       ** no projection applied **       ",
     *            "  |ref>F<ref|+|sin>F<sin|+|doub>F<doub|  ",
     *            "      |ref>F<ref| + |doub>F<doub|        ",
     *            "   |ref>F<ref| + |sin+doub>F<sin+doub|   ",
     *            "  |ref>F<ref|+|sin>F<sin|+|rest>F<rest|  ",
     *            "      |ref>F<ref| + |rest>F<rest|        ",
     *            "  |ref>F<ref| + |doub+more>F<doub+more|  ",
     *            "|ref>F<ref|+|sin^>F<sin^|+|doub^>F<doub^|",
     *            "      |ref>F<ref| + |doub^>F<doub^|      ",
     *            " |ref>F<ref| + |sin^+doub^>F<sin^+doub^| ",
     *            " ",
     *            " ",
     *            " ",
     *            " "/
      data cnum/'1','2','3','4','5','6','7','8','9'/
c
      mp = 0
      icci = 1
      isec_mp = 0
      nummp = 0
c     maxcyc = 0
      critmp = 1.0d-6
      maxcmp = 50
      iprimp = 1
      mpmod = 1
      iwnmp(1,1) = 0
      indmp(6) = 0
      plevel(1) = .0d0
      plevel(2) = plevel(1)
      fignor = 1.0d-8
      cortho = 1.0d-8
      iortho = 1
      ifdiaf = 1
      ifcepa = 0
      icepa  = -1
c...  disable configuration printing (use prconf after to switch on)
      iprcon = -1
c...  select the reference function as trial function
      istcc = 1
      iselcc = 3
c
c...  set default diis (start straight away)
c
      if (diiscr.eq.0.0d0) then
         diiscr = 100.0d0
         ntdiis = -1
      end if
c
10    call inpa(test)
c
      if (test(1:4).eq.'mode') then
         call inpi(mpmod)
         mpmod = abs(mpmod)
         if (mpmod.lt.0.or.mpmod.gt.nmodel-1) call caserr('wrong model')
         if ((mpmod.ne. 1).and.(mpmod.ne. 2).and.(mpmod.ne. 3).and.
     *       (mpmod.ne.11).and.(mpmod.ne.12).and.(mpmod.ne.13))
     *       call caserr('no such model')
      else if (test(1:4).eq.'prin') then
         call inpi(iprimp)
      else if (test(1:4).eq.'crit') then
         call inpa(test)
         if (test(1:4).eq.'abso') then
           call inpf(bilbo)
           call inpi(j)
           critmp = bilbo*(0.1d0**j)
         else if (test(1:4).eq.'rela') then
           call inpf(bilbo)
           call inpi(j)
           critmp = - bilbo*(0.1d0**j)
         else
           jrec = jrec - 1
           call inpf(bilbo)
           call inpi(j)
           critmp = bilbo*(0.1d0**j)
         endif
         if(mp.le.3) then
           critmp = dabs(critmp)
         endif
      else if (test(1:4).eq.'fign') then
         call inpf(bilbo)
         call inpi(j)
         fignor = bilbo*(0.1d0**j)
      else if (test(1:8).eq.'schmidtp') then
         iortho = 2
         call inpf(bilbo)
         call inpi(j)
         cortho = bilbo*(0.1d0**j)
      else if (test(1:8).eq.'schmidt4') then
         iortho = 3
         call inpf(bilbo)
         call inpi(j)
         cortho = bilbo*(0.1d0**j)
      else if (test(1:7).eq.'schmidt') then
         iortho = 1
         call inpf(bilbo)
         call inpi(j)
         cortho = bilbo*(0.1d0**j)
      else if (test(1:6).eq.'lowdin') then
         iortho = 4
         call inpf(bilbo)
         call inpi(j)
         cortho = bilbo*(0.1d0**j)
      else if (test(1:5).eq.'house') then
         iortho = 5
         call inpf(bilbo)
         call inpi(j)
         cortho = bilbo*(0.1d0**j)
      else if (test(1:4).eq.'maxc') then
         call inpi(maxcmp)
      else if (test(1:4).eq.'leve') then
         call inpf(plevel(1))
         call inpa(test)
         jrec = jrec - 1
         if (locatc(cnum,9,test(1:1)).ne.0) then
            call inpi(levelp)
            call inpf(plevel(2))
         else
            levelp = 0
            plevel(2) = plevel(1)
         end if
      else if (test(1:4).eq.'nodo') then
         iwnmp(1,1) = 1000
         print *,' ***debug*** no doublys are recognised'
      else if (test(1:4).eq.'nodi') then
         ifdiaf = 0
         print *,' FOCK matrix assumed **non**-diagonal'
      else if (test(1:4).eq.'noor') then
         print *,' no noorth option !!'
      else if (test(1:4).eq.'fock') then
         call inpa(test)
         nummp = locatc(yed,maxlfn,test(1:4))
         if (nummp.eq.0) call caserr('fock-subdirective unknown')
         call inpi(iblk_mp)
         call inpi(isec_mp)
      else if (test(1:4).eq.'tria') then
         call inpa4(test)
         if (test(1:4).eq.'ref ') then
           iselcc = 3
         else if (test(1:4).eq.'vac ') then
           iselcc = 4
         endif
         call inpi(istcc)
         istcc = max(1,istcc)
      else if (test(1:1).eq.' ') then
         go to 20
      else
         call caserr(' subdirective not recognised ')
      end if
c
      go to 10
c
c...  print input
c
20    write(6,1)
      write(6,2) mpmod,descr(mpmod+1)
      write(6,3) iprimp
      write(6,4) maxcmp,critmp
      if (plevel(1).ne.0.0d0.or.plevel(2).ne.0.0d0)
     *   write(6,6) plevel(1),levelp,plevel(2)
      write(6,8) fignor
      if (iortho.eq.1) then
        write(6,9) 'Schmidt ...................',cortho
      else if (iortho.eq.2) then
        write(6,9) 'Pivoted Repeated Schmidt ..',cortho
      else if (iortho.eq.3) then
        write(6,9) 'Quadruple Precision Schmidt',cortho
      else if (iortho.eq.4) then
        write(6,9) 'Lowdin ....................',cortho
      else if (iortho.eq.5) then
        write(6,9) 'Householder QR ............',cortho
      endif
      if (nummp.gt.0) then
         write(6,5) yed(nummp),iblk_mp,isec_mp
      else
         write(6,7)
      end if
c....
1     format(//,
     1 '  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$',/,
     2 '  $ Multi Reference CI ',i1,' option requested CONTRACT      $')
2     format(
     1 '  $ model',i3,1x,a43,                                       '$')
3     format(
     1 '  $    print requested  ',i6,     '                          $')
4     format(
     1 '  $    # cycles ',i6,'  crit ', e9.3,     '                  $')
6     format(
     1 '  $    level-shift ',e8.3,' ; ',i6,'  ,  ',e8.3,     '       $')
8     format(
     1 '  $    ignore FOCK-integrals < ',e8.3,     '                 $')
9     format(
     1 '  $    ',a27, ' upto ',e8.3,                        '        $')
5     format(
     1 '  $    FOCK-matrix from ',a3,i5,i4,'                        $',
     1/,'  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
7     format(
     1 '  $    FOCK-matrix will be calculated                   $',/,
     1 '  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
69    format(
     1 '  $ H0 block diagonal in excitation classes !!!         $')
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine ccepin
c
      implicit real*8 (a-h,o-z), integer (i-n)
c
c...  input routine for Contracted MRCEPA
c...
c...  called from inpro
c...
c...   = directive format = (all on 1 line)
c...
c...   ccepa (model) (print 1) (precrit (absolute/relative) 1.0 6)
c...       (level 0.0 (5 0.0)) (maxcyc 50) (fignore 1. 8)
c...       (nodoc) (nodiag) (sin/nosin) (it 3) (crit 1.0d-2)
c...       (micro 3 0.01d0) (brst) (schmidt 1. 8)
c...
c...   model : one of : mrd, mr1,mr0,aqcc,acpf 
c...   print : more output for a higher number (part. >2000)
c...   crit : stop-criterion for preconditioner-iteration
c...          (max. residue)
c...     absolute : use stop-criterion as is.
c...     relative : use stop-criterion relative to the norm of the
c...                righthand side of linear system. (only in MPn)
c...   level : levelshift1, # iterations, levelshift2
c...           after the specified # iterations shift2 is used
c...   maxcyc : maximum # cycles in precondition-iteration-process
c...   fignore : absolute smaller FOCK-matrix elements are ignored
c...   cortho  : criterion for orthogonality and linear dependency in
c...             modified Gramm-Schmidt procedure and QR procedure.
c...   iortho  : Orthogonalisation method 1=Schmidt 2=QR.
c...   nodoc : tells program to treat doc as voc (for debugging)
c...   nodiag : tells the program to assume the external FOCK-
c...            integral-matrix to be non-diagonal (for debugging)
c...   = only the first 4 chars of a keyword are read =
c...
c...   (nodoc is transmitted to mpini etc by iwnmp(1,1) = 1000)
c...
c...
c...  various routines from DIRECT and the atmol-library are also used
c...  routines from DIRECT really changed :
c...     block-data,direct,inpro,pass
c...
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
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
      common/mr_comcci/icci,iclass(3,3,3),istcc,iselcc
c
      real*8 critcep, crycep
      integer ifcepa, icepa , lcepaul, lsinsh, itcepa, lprcep, indcep
      integer nrfcsf, mcycep, ibrst,icorr_mr,mrdoc
      common/ccepa/critcep,crycep,ifcepa,icepa,lcepaul,lsinsh,itcepa,
     1                lprcep,indcep,nrfcsf,mcycep,ibrst,icorr_mr,mrdoc
c
      common/mr_mp2l/ ivoc(62),nvoc,idoc(62),ndoc,kvcc(4),nvcc(4),
     *                ncmphn(10,18),ncmph1(10,18,8),ncmph2(10,18,8),
     *                ncmphd(10,18)
      logical malt
      real*8 potnuc,totnuc,e11,eigen,h11,z11,h22,s22,h12,s12,
     *     tester,eshift,thresh,aaa,bbb,scream,gprin,z12,z22,
     *     vdstpr,qester
      integer maxcyc,moddia,mprin,icyc
      common/erg/potnuc,totnuc,e11,eigen(nd200),h11
     *,z11,h22,s22,h12,s12,tester,eshift,thresh,aaa,bbb,scream
     *,gprin,z12,z22,vdstpr(3),qester,maxcyc,moddia,mprin,malt,icyc
c
      common/diisdc/ diiscr,ntdiis
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
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
      character*8 test, test1
      character*25 text
      character*1 cnum(9)
      parameter (nmodel=14)
      character*43 descr(nmodel)
      data descr /"       ** no projection applied **       ",
     *            "  |ref>F<ref|+|sin>F<sin|+|doub>F<doub|  ",
     *            "      |ref>F<ref| + |doub>F<doub|        ",
     *            "   |ref>F<ref| + |sin+doub>F<sin+doub|   ",
     *            "  |ref>F<ref|+|sin>F<sin|+|rest>F<rest|  ",
     *            "      |ref>F<ref| + |rest>F<rest|        ",
     *            "  |ref>F<ref| + |doub+more>F<doub+more|  ",
     *            "|ref>F<ref|+|sin^>F<sin^|+|doub^>F<doub^|",
     *            "      |ref>F<ref| + |doub^>F<doub^|      ",
     *            " |ref>F<ref| + |sin^+doub^>F<sin^+doub^| ",
     *            " ",
     *            " ",
     *            " ",
     *            " "/
      data cnum/'1','2','3','4','5','6','7','8','9'/
cci
      icci = 1
cci
cmp
      mp = 0
      isec_mp = 0
      nummp = 0
c     maxcyc = 0
      critmp = 1.0d-6
      maxcmp = 50
      iprimp = 1
      mpmod = 1
      iwnmp(1,1) = 0
      indmp(6) = 0
      plevel(1) = .0d0
      plevel(2) = plevel(1)
      fignor = 1.0d-8
      cortho = 1.0d-8
      iortho = 1
      ifdiaf = 1
cmp
ccepa
      icepa = 10
      text = ' mrdcepa variant '
      lcepaul = 0
      lsinsh = -1
      itcepa = 3
      critcep = 1.0d-2
      lprcep = 0
      mcycep = 3
      crycep = 0.01d0
      ibrst = 1
ccepa
c...  disable configuration printing (use prconf after to switch on)
      iprcon = -1
c
c...  set default diis (start straight away)
c
      if (diiscr.eq.0.0d0) then
         diiscr = 100.0d0
         ntdiis = -1
      end if
c
100   call inpa(test)
c
      if (test(1:3).eq.'mrd') then
         icepa = 10
         text = ' mrdcepa variant '
      else if (test(1:3).eq.'mr1') then
         icepa = 12
         text = ' mrcepa1-all variant '
      else if (test(1:3).eq.'mr0') then
         text = ' mrcepa0 variant '
         icepa = 100
      else if (test(1:4).eq.'aqcc') then
         text = ' aqcc variant '
         icepa = 101
      else if (test(1:4).eq.'acpf') then
         text = ' acpf variant '
         icepa = 102
      else if (test(1:3).eq.'doc') then
         text = ' mrcepa1-doc variant '
         if (icepa.ne.11.and.icepa.ne.12)
     1   call caserr('doc subdirective only for mrcepa1')
         icepa = 11
      else if (test(1:3).eq.'all') then
         text = ' mrcepa1-all variant '
         if (icepa.ne.11.and.icepa.ne.12)
     1   call caserr('all subdirective only for mrcepa1')
         icepa = 12
      else if (test(1:4).eq.'prin') then
         call inpi(iprimp)
         lprcep = 1
      else if (test(1:4).eq.'prec') then
         call inpa(test)
         if (test(1:4).eq.'abso') then
           call inpf(bilbo)
           call inpi(j)
           critmp = bilbo*(0.1d0**j)
         else if (test(1:4).eq.'rela') then
           call inpf(bilbo)
           call inpi(j)
           critmp = - bilbo*(0.1d0**j)
         else
           jrec = jrec - 1
           call inpf(bilbo)
           call inpi(j)
           critmp = bilbo*(0.1d0**j)
         endif
         if(mp.le.3) then
           critmp = dabs(critmp)
         endif
      else if (test(1:4).eq.'fign') then
         call inpf(bilbo)
         call inpi(j)
         fignor = bilbo*(0.1d0**j)
      else if (test(1:4).eq.'schm') then
         call inpf(bilbo)
         call inpi(j)
         iortho = 1
         cortho = bilbo*(0.1d0**j)
      else if (test(1:2).eq.'house') then
         call inpf(bilbo)
         call inpi(j)
         iortho = 5
         cortho = bilbo*(0.1d0**j)
      else if (test(1:4).eq.'maxc') then
         call inpi(maxcmp)
      else if (test(1:4).eq.'leve') then
         call inpf(plevel(1))
         call inpa(test)
         jrec = jrec - 1
         if (locatc(cnum,9,test(1:1)).ne.0) then
            call inpi(levelp)
            call inpf(plevel(2))
         else
            levelp = 0
            plevel(2) = plevel(1)
         end if
      else if (test(1:4).eq.'nodo') then
         iwnmp(1,1) = 1000
         print *,' ***debug*** no doublys are recognised'
      else if (test(1:4).eq.'nodi') then
         ifdiaf = 0
         print *,' FOCK matrix assumed **non**-diagonal'
      else if (test(1:4).eq.'noor') then
         print *,' no noorth option !!'
      else if (test(1:4).eq.'fock') then
         call inpa(test)
         nummp = locatc(yed,maxlfn,test(1:4))
         if (nummp.eq.0) call caserr('fock-subdirective a mystery')
         call inpi(iblk_mp)
         call inpi(isec_mp)
      else if (test(1:3).eq.'sin') then
         lsinsh = 1
      else if (test(1:4).eq.'nosi') then
         lsinsh = 0
      else if (test(1:2).eq.'it') then
         call inpi(itcepa)
      else if (test(1:4).eq.'crit') then
         call inpf(critcep)
      else if (test(1:4).eq.'micr') then
         call inpi(mcycep)
         call inpf(crycep)
         mcycep = max(mcycep,1)
         if (crycep.eq.0.0d0) crycep = 1.0d-15
      else if (test(1:4).eq.'brst') then
         ibrst = 1
         print *,' brst is the only reasonable option'
      else if (test(1:1).eq.' ') then
         go to 200
      else
         call caserr('unrecognised option for contracted (mr)cepa')
      end if
c
      go to 100
c
c...  print input
c
200   continue
      if (lsinsh.lt.0) lsinsh = 1
c...  always single shift  / no consistent way to turn it off yet
      write(6,1) text
      if (lsinsh.gt.0) write(6,10)
      if (lsinsh.le.0) write(6,11)
      if (ibrst.gt.0)  write(6,14)
      write(6,12) mcycep,crycep
      write(6,13) itcepa,critcep
      write(6,2) mpmod,descr(mpmod+1)
      write(6,3) iprimp
      write(6,4) maxcmp,critmp
      if (plevel(1).ne.0.0d0.or.plevel(2).ne.0.0d0)
     *   write(6,6) plevel(1),levelp,plevel(2)
      write(6,8) fignor
      if (nummp.gt.0) then
         write(6,5) yed(nummp),iblk_mp,isec_mp
      else
         write(6,7)
      end if
c....
1     format(//,
     1       '  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$',/,
     2       '  $ Multi Reference Contracted CEPA requested       $',/,
     3       '  $    ',a25,                  '                    $')    
10    format('  $ include shifting of singles                     $')
11    format('  $ exclude shifting of singles                     $')
12    format('  $ number of micro iterations     ',i6,'           $',/,
     1       '  $    until relative difference < ',f8.5,'         $')
13    format('  $ start cepa after it ',i6,'                      $',/,
     1       '  $    and tester      <',f8.5,'                    $')
14    format('  $ include brillioun states                        $')
2     format('  $ model',i2,1x,a39,                             ' $')
3     format('  $    print requested  ',i6,'                      $')
4     format('  $    # cycles ',i6,'  crit ', e9.3,'              $')
6     format('  $    level-shift ',e8.3,' ; ',i6,'  ,  ',e8.3,'   $')
8     format('  $    ignore FOCK-integrals < ',e8.3,'             $')
5     format('  $    FOCK-matrix from ',a3,i5,i4,'                $',/,
     1       '  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
7     format('  $    FOCK-matrix will be calculated               $',/,
     1       '  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine davidc(cr,icr,cc)
c
      implicit real*8 (a-h,o-z), integer(i-n)
c
      integer icr
      dimension cr(*), icr(*), cc(*)
c
c...  contracted Davidson diagonalisation.
c
c     disk -layout
c     R, HR (CSF) : index(26),index(27)
c     contracted (excited state) basis :
c     (H-F)R      : indexp(1)
c     Fock-diag   : indexp(2)
c     x and Fx    : indexp(3),indexp(4)
c     diis storage : indexp(5-10) (not used currently)
c     storage for contracted c,z vectors : indexp(11-110)
c     indexp(ic(i)) : disk-address for (i-1)-th order correction to the
c                     wavefunction c(i-1).
c     indexp(iz(i)) : disk-address for z(i-1) = H*c(i-1)
c
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
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
      common/mr_mp2l/ ivoc(62),nvoc,idoc(62),ndoc,kvcc(4),nvcc(4),
     *                ncmphn(10,18),ncmph1(10,18,8),ncmph2(10,18,8),
     *                ncmphd(10,18)
      common /mr_mp2con/ nsmpv,nsmp1,nsmp2,nsmp2d,ntotmp,
     *                   nmpv,nmp1,nmp2,nmp2d
      common/diisdc/ diiscr,ntdiis
c
c
      integer nref, nnst, nmin1, nmin2, mspin, nelec
      integer nstrt1, nlast1, nstrt2, nlast2, nmin12, nminv2, nminv1
      integer nlst21, nlastx
      common /ccntl/ nref,nnst,nmin1,nmin2,mspin,nelec,
     +               nstrt1,nlast1,nstrt2,nlast2,nmin12,nminv2,nminv1,
     +               nlst21,nlastx
c
c
      integer nval, ndoub, nsingl, ntripl, ntotal, ibig, jbig, kbig
      integer nvd, nvds, nsitri, lbig, mbig, nbig, nref0
      common /cntrl/ nval,ndoub,nsingl,ntripl,ntotal,ibig,jbig,kbig,
     +               nvd,nvds,nsitri,lbig,mbig,nbig,nref0
c
c
      integer nint, next, ninnex, nwm1, nw, nwp1, nintwo
      integer nint21, nintp1
      common /orb/ nint,next,ninnex,nwm1,nw,nwp1,nintwo,
     +             nint21,nintp1
c
c
      integer nxblk, index, numci, mblkci, ispbig, ispsma
      common/disktl/nxblk,index(199),numci,mblkci,ispbig,ispsma
      common/stoctl/khbill,khbz,khbc,khbb,khbzd,khbzs,khbzt,khbcd,
     *khbcs,khbct,kh,khb1,khb2,khb3,khb4,khb5,khb6,khb7,khb8
     *,khb44,khb55,khb66,khb88
c
      integer lcall,  mrstor
      logical modev, moded, modest, instor, loijka, msort, loijab
      logical loiajb, lvar
      common/diactl/modev,moded,modest,instor,
     +              loijka,msort,loijab,loiajb,lcall,mrstor,lvar
c
      logical malt
      real*8 potnuc,totnuc,e11,eigen,h11,z11,h22,s22,h12,s12,
     *     tester,eshift,thresh,aaa,bbb,scream,gprin,z12,z22,
     *     vdstpr,qester
      integer maxcyc,moddia,mprin,icyc
      common/erg/potnuc,totnuc,e11,eigen(nd200),h11
     *,z11,h22,s22,h12,s12,tester,eshift,thresh,aaa,bbb,scream
     *,gprin,z12,z22,vdstpr(3),qester,maxcyc,moddia,mprin,malt,icyc
c
      integer nirr,mult,ifrep,nirr8,jfrep
      common /table/ nirr,mult(8,8),ifrep,nirr8,jfrep
c
c
      integer nconst, ncond, maxpa, maxpat, maxsin, maxepa
      integer mxdvpa, mxddpa, nspips, nspipt, mxsspa, mxttpa
      integer mxstpa, nsitot, ntrtot
      integer mx12pa, mxvvpa, mxv2pa, maxpsq, max12
      common /auxh/ nconst(8),ncond(8),maxpa(4),maxpat,maxsin(3),
     +              maxepa(4),mxdvpa,mxddpa,nspips(8),nspipt(8),
     +              mxsspa,mxttpa,mxstpa,nsitot(8),ntrtot(8),
     +              mx12pa,mxvvpa,mxv2pa,maxpsq,max12
c
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
c
      integer nuse, ngot
      common /corctl/ nuse, ngot
c
c
      common/junk/vdum(510),ivdum(66),q(60,60),f(1830),eig(60),
     +            fsloc(820),f2sloc(820)
      common/dialoc/floc(820),f2loc(820),xn(40),xns(40),ilifq(40)
c
      integer iro, iky, niky, m1e, m2e, m1ia, m2coul
      integer m2exch, m2fock, m2psup, m2qsup
      integer mbase, nbase, npairi, norbi, mbasi, mnbasi
      integer iap, mbfock, mbexch, mbcoul, mbpsup, mbqsup
      integer norbas, m2exab, m2em1e
      common /helpr/ iro(nd200),iky(nd200),niky(8),m1e,m2e,m1ia,m2coul,
     +   m2exch,m2fock,m2psup,m2qsup,
     +   mbase(288),nbase(8),npairi(36),norbi(8),mbasi(8),mnbasi(8),
     +   iap(62),mbfock(62),mbexch(36),mbcoul(36),mbpsup(36),mbqsup(36),
     +   norbas(nd200),m2exab,m2em1e
c
ccepa
c
      real*8 critcep, crycep
      integer ifcepa, icepa , lcepaul, lsinsh, itcepa, lprcep, indcep
      integer nrfcsf, mcycep, ibrst,icorr_mr,mrdoc
      common/ccepa/critcep,crycep,ifcepa,icepa,lcepaul,lsinsh,itcepa,
     1                lprcep,indcep,nrfcsf,mcycep,ibrst,icorr_mr,mrdoc
c
      common/mr_cepctl/ kcpill,kcpv,kcpz,kcpx,kcpcd,kcpd,kcp2,kcpc2,
     1                  kcepc,kcepz,kcpcc,ncpcc,kcpbl,ncpbl,icpsto
     1                  ,kdavcz
      logical icpsto
      integer maxpair
      parameter (maxpair=40)
      real*8 c0,h11ci,ecorr_mr epair
      common /mr_pairs/ c0,h11ci,ecorr_mr,epair(maxpair,maxpair)
c
      real*8 eshifc eshifd
      common /shifts/ eshifc(maxpair,maxpair),eshifd(maxpair)
ccepa
      character*8 text
      integer maxord
      parameter(maxord = 40)
c...  maxord = maximum order in perturbation theory
      dimension indexp(10+2*maxord+2), ic(maxord+1), iz(maxord+1)
c
      data text/'      '/
c
      pplev = plevel(1)
c
      call cpuwal(top,bot)
      write(6,444)top,bot
444   format(//1x,80('-')//
     *' ****MR-CCI/CCEPA v1.3 called at',f12.1,' wall',f12.1,' secs')
chvd
c...  save config and common blocks for later use in matrix multiplication
c
      call cistor
      call gciblk(mpcblk)
      call wrt3(cr,5*(nnst+nmin1+nmin2),mpcblk,numci)
chvd
c...  contract the Z-vector
c
      nuse = khbill + ntotal
      call corfai('contract')
c...  get z-vector
      kmpz = khbill
      if (instor) then
         if (khbz.ne.khbill) call caserr('khbz misplaced in mpgen')
      else
         call getcz(cr,cr(kmpz+1),index(27))
      end if
c...  dump back to disk  at index(27) straight
      call wrt3(cr(kmpz+1),ntotal,index(27),numci)
c
c...  switch to new contracted addressing (khbill changes)
c
      call basmp(icr)
c
c...  build disc addresses
c
      nbl = (ntotmp-1)/511+1
      if (icepa.ge.0) then
c...    cepa uses disc space through dumpcc and getcc
        indexp(1) = kcpbl + (maxord*(maxord+1)/2)*ncpbl
      else
        indexp(1) = index(28)
      endif
      do 5 i=2,10+2*maxord+2
5     indexp(i) = indexp(i-1) + nbl
      ic(1) = 10 + 1
      iz(1) = 10 + 2
      do 6 i=2,maxord+1
        ic(i) = ic(i-1) + 2
        iz(i) = iz(i-1) + 2
6     continue
      ibestc = indexp(ic(maxord+1))
      ibestz = indexp(iz(maxord+1))
c
c...  contract
c
c     kmpz = khbill
      if (icepa.ge.0) then
c...    cepa uses core to calculate <cc>
        kmpz = kcpcc+5*ncpcc+1
      else
        kmpz = khbill
      endif
      kmpc = kmpz + ntotmp
      kmpb = kmpc + ntotmp
      kmpy = kmpb + ntotmp
      kcc  = kmpz
      kdc  = kcc  + nsmpv+nsmp1+nsmp2
      kc   = kmpz + ntotmp
      kvc  = kc   + ntotal
      kscr = kvc  + max(nvcc(1),nvcc(2),nvcc(3),nvcc(4))
      nuse = kscr + max(nspips(1)*next,norbtr(1)*nspips(1))
      call corfai('contract')
chvd
      write(6,446)
446   format(//3x,89('-')/'   order  time(wall)  time(secs)',
     *10x,'energy(au) perturbation(au)'
     */3x,89('-'))
chvd
c...  save config and common blocks for later use in contracted
c...  stuff
c
      call mpstor
      call gmpblk(mpcblk)
      call wrt3(cr,5*(nmpv+nmp1+nmp2+nmp2d),mpcblk,numci)
chvd
chvd
c...  retrieve, contract, and store zero-th order c-vector
c
      call rdedx(cr(kc+1),ntotal,index(26),numci)
      if (iprimp.ge.2000) write(6,616) (cr(kc+i),i=1,ntotal)
616   format(' - C0 (CSF) -',/,(2x,10e12.5))
      call ccntmp(cr(kc+1),cr(kcc+1),cr(kdc+1),cr(kvc+1),cr(kscr+1))
      if (iprimp.ge.2000) write(6,614) (cr(kcc+i),i=1,ntotmp)
614   format(' - C0 (excited states) -',/,(2x,5e22.15))
      call wrt3(cr(kcc+1),ntotmp,indexp(ic(1)),numci)
      call wrt3(cr(kcc+1),ntotmp,ibestc,numci)
c
c...  retrieve, contract, and store zero-th order z-vector
c
      call rdedx(cr(kc+1),ntotal,index(27),numci)
      if (iprimp.ge.2000) write(6,603) (cr(kc+i),i=1,ntotal)
603   format(' - HC0 (CSF) -',/,(2x,10e12.5))
      call ccntmp(cr(kc+1),cr(kcc+1),cr(kdc+1),cr(kvc+1),cr(kscr+1))
      if (iprimp.ge.2000) write(6,602) (cr(kcc+i),i=1,ntotmp)
602   format(' - HC0 (excited states) -',/,(2x,5e22.15))
      call wrt3(cr(kcc+1),ntotmp,indexp(iz(1)),numci)
      call wrt3(cr(kcc+1),ntotmp,ibestz,numci)
c
c...  calculate diagonal of fock-matrix
c
      call vclr(cr(kmpz+1),1,ntotmp)
      call mpdiag(icr,cr(kmpz+1),cr(kscr+1))
c
      if (iprimp.ge.2000) write(6,605) (cr(kmpz+i),i=1,ntotmp)
605   format(' - diagonal (F-E0) -',/,(2x,5e22.15))
c
c...  dump diagonal (F-E0)
c
      call wrt3(cr(kmpz+1),ntotmp,indexp(2),numci)
c
chvd  SETUP START
c
      ispace = maxord
      call ibasgn(ispace,0,60,ilifq)
9999  call rdedx(cr(kmpc+1),ntotmp,indexp(ic(1)),numci)
      call rdedx(cr(kmpz+1),ntotmp,indexp(iz(1)),numci)
      floc(1)   = ddot(ntotmp,cr(kmpc+1),1,cr(kmpz+1),1)
      var       = ddot(ntotmp,cr(kmpz+1),1,cr(kmpz+1),1)
      f2loc(1)  = var
      fsloc(1)  = floc(1)
      h11       = floc(1)
      f2sloc(1) = f2loc(1)
ccepa
      if (icepa.ge.0) then
        if (moddia.gt.0) then
           call cntrcc(cr,cr(kmpc+1),cr(kmpz+1),cc,1,2,temp)
        else
           call vclr(cc,1,ncpcc*2)
           call cntrcc(cr,cr(kmpc+1),cr(kmpc+1),cc,0,1,temp)
        endif
        call dumpcc(cc,1,1)
        if (ifcepa.eq.1) then
          call hc2cep(cc,ncpcc,h11,h11,var,var,1)
        endif
      endif
c
      ncycep = 0
      h11o   = h11
ccepa
c
c     BEGIN DAVIDSON ITERATONS
c
      do 6666 joop=2,ispace
c
      var = var - h11*h11
      if (ifcepa.gt.0) then
        call ccpvds(cr,cr(kmpz+1),cr(kmpc+1))
        var = ddot(ntotmp,cr(kmpz+1),1,cr(kmpz+1),1) - h11*h11
      endif
      call gtrian(ntotmp,h11,cr(kmpc+1),cr(kmpz+1),cr(kmpc+1))
      tester = 0.0d0
      do i=1,ntotmp
         tester = max(tester,dabs(cr(kmpc+i)))
      end do
      call rdedx(cr(kmpz+1),ntotmp,indexp(2),numci)
c     call vsav(ntotmp,eshift-h11,cr(kmpz+1),cr(kmpz+1))
c...  mpdiag contructs diag(F-E0) already.
      call mrgle(ntotmp,0.05d0,cr(kmpz+1))
      call vvdv(ntotmp,cr(kmpc+1),cr(kmpc+1),cr(kmpz+1))
      call mrgge(ntotmp,0.2d0,cr(kmpc+1))
      call cpuwal(top,bot)
      write(6,111)icyc,top,bot,h11,tester,var,text
111   format(i8,f12.0,f12.3,e23.14,e18.9,e19.10,a6)
      if(tester.lt.thresh) goto 777
      if(icyc.ge.maxcyc)   goto 777
c
c...  orthogonalize pert vector
c
      x = ddot(ntotmp,cr(kmpc+1),1,cr(kmpc+1),1)
8000  x = 1.0d0/dsqrt(x)
      call dscal(ntotmp,x,cr(kmpc+1),1)
      do 8010 i = 1, joop-1
        call rdedx(cr(kmpb+1),ntotmp,indexp(ic(i)),numci)
        x = ddot(ntotmp,cr(kmpc+1),1,cr(kmpb+1),1)
        call daxpy(ntotmp,-x,cr(kmpb+1),1,cr(kmpc+1),1)
8010  continue
      x = ddot(ntotmp,cr(kmpc+1),1,cr(kmpc+1),1)
      if(x.lt.0.05d0) go to 8000
      x = 1.0d0/dsqrt(x)
      call dscal(ntotmp,x,cr(kmpc+1),1)
c
      call wrt3(cr(kmpc+1),ntotmp,indexp(ic(joop)),numci)
      call genzcc(cr,h11t,indexp(ic(joop)),indexp(iz(joop)),
     *            kcc,kc,kdc,kvc,kscr)
      call rdedx(cr(kmpc+1),ntotmp,indexp(ic(joop)),numci)
      call rdedx(cr(kmpz+1),ntotmp,indexp(iz(joop)),numci)
c
      m  = iky(joop)
      mm = m + joop
      floc(mm)  = h11t
      f2loc(mm) = ddot(ntotmp,cr(kmpz+1),1,cr(kmpz+1),1)
ccepa
      if (icepa.ge.0) then
c...     c.c and c.z for current c/z vectors
         if (moddia.gt.0) then
           call cntrcc(cr,cr(kmpc+1),cr(kmpz+1),cc,1,2,temp)
         else
           call vclr(cc,1,ncpcc*2)
           call cntrcc(cr,cr(kmpc+1),cr(kmpc+1),cc,0,1,temp)
         endif
         call dumpcc(cc,joop,joop)
      endif
ccepa
c
      do 8020 i = 1, joop-1
        call rdedx(cr(kmpb+1),ntotmp,indexp(ic(i)),numci)
        call rdedx(cr(kmpy+1),ntotmp,indexp(iz(i)),numci)
ccepa
        if (icepa.lt.0) then
          floc(m+i)  = ddot(ntotmp,cr(kmpc+1),1,cr(kmpy+1),1)
        else
          call cntrcc(cr,cr(kmpc+1),cr(kmpb+1),cc,0,1,temp)
          if (moddia.gt.0) then
            call cntrcc(cr,cr(kmpc+1),cr(kmpy+1),cc,0,2,temp)
            call cntrcc(cr,cr(kmpz+1),cr(kmpb+1),cc,0,3,floc(m+i))
          else
            call vclr(cc(ncpcc+1),1,ncpcc*2)
            floc(m+i) = ddot(ntotmp,cr(kmpc+1),1,cr(kmpy+1),1)
          endif
          call dumpcc(cc,joop,i)
        endif
        f2loc(m+i) = ddot(ntotmp,cr(kmpz+1),1,cr(kmpy+1),1)
ccepa
8020  continue
c
      call dcopy(joop,floc(m+1),1,fsloc(m+1),1)
      call dcopy(joop,f2loc(m+1),1,f2sloc(m+1),1)
c
8025  if (moddia.eq.0) then
c...     lock mode
         call dcopy(mm,floc,1,f,1)
         call jacobi(f,iky,joop,q,ilifq,joop,eig,2,0,1.0d-9)
         do 8030 loop=1,joop
            eig(loop)=-dabs(q(1,loop))
8030     continue
      else if (moddia.eq.-1) then
c...     emin mode
         call dcopy(mm,floc,1,f,1)
         call jacobi(f,iky,joop,q,ilifq,joop,eig,2,1,1.0d-9)
      else if (moddia.eq.1) then
c...     vmin mode
         do i=1,mm
            f(i) = f2loc(i) - (h11+h11)*floc(i)
         end do
         call jacobi(f,iky,joop,q,ilifq,joop,eig,2,1,1.0d-9)
         call vsadd(eig,1,h11*h11,eig,1,joop)
      endif
      ivrs = 0
      do i=2,joop
        if (eig(i).lt.eig(ivrs+1)) ivrs = i-1
      end do
      var  = expect(f2loc,q(1,ivrs+1),joop)
      h11  = expect(floc,q(1,ivrs+1),joop)
c
c...  calculate new c,z-vectors
c
      call dcopy(joop,q(1,ivrs+1),1,eig,1)
      call rdedx(cr(kmpc+1),ntotmp,indexp(ic(1)),numci)
      call dscal(ntotmp,eig(1),cr(kmpc+1),1)
      call rdedx(cr(kmpz+1),ntotmp,indexp(iz(1)),numci)
      call dscal(ntotmp,eig(1),cr(kmpz+1),1)
      do 40 i = 2, joop
        call rdedx(cr(kmpb+1),ntotmp,indexp(ic(i)),numci)
        call daxpy(ntotmp,eig(i),cr(kmpb+1),1,cr(kmpc+1),1)
        call rdedx(cr(kmpb+1),ntotmp,indexp(iz(i)),numci)
        call daxpy(ntotmp,eig(i),cr(kmpb+1),1,cr(kmpz+1),1)
40    continue
      call wrt3(cr(kmpc+1),ntotmp,ibestc,numci)
      call wrt3(cr(kmpz+1),ntotmp,ibestz,numci)
ccepa
      if (icepa.ge.0.and.tester.le.critcep.and.joop.ge.itcepa) then
        ifcepa = 1
        text   = ' cepa '
      endif
      if (ifcepa.eq.1.and.ncycep.le.mcycep) then
c
        h11dif = h11 - h11o
        call cpuwal(top,bot)
        if (lprcep.gt.0.and.ncycep.gt.0) then
          write(6,114)ncycep,top,bot,h11o,temp
114       format(' cepa micro it.',i3,' at ',f12.3,' / ',f12.3,
     *       ' last e ',e20.13,' delta ',e10.3,'  cepa ',/)
        endif
        h11dif = dabs(h11dif)
        h11o   = h11
c
c...    expand the excitated-state vector to configuration-state vector
c
        call rdedx(cr(kcc+1),ntotmp,ibestc,numci)
        call ccxpmp(cr(kc+1),cr(kcc+1),cr(kdc+1),cr(kvc+1),cr(kscr+1))
c
c...    restore original information in the common blocks to enable
c...    the use of genz.
c
        call ciretr
        nuse = max(khbz,khbc,khbb)+ntotal
        call corfai('contract')
c
c...    move the csf-vector to proper position for genz.
c
        if(kc.gt.khbc) then
          do 80 i = 1, ntotal
            cr(khbc+i) = cr(kc+i)
80        continue
        else if(kc.lt.khbc) then
          do 90 i = 0, ntotal-1
            cr(khbc+ntotal-i) = cr(kc+ntotal-i)
90        continue
        endif
c
c...    retrieve config information to do genz
c
        call gciblk(mpcblk)
        call rdedx(cr,5*(nnst+nmin1+nmin2),mpcblk,numci)
c
c...    dump current vector to synchronize dumpfile contents with
c...    core contents
c
        call dumpcz(cr,cr(khbc+1),index(26))
c
c...    compute pair energies
c
        call cepair(index(26),cr,cr)
c
c...    what was the ci-energy again ?
c
        h11ci = expect(fsloc,q(1,ivrs+1),joop)
c
c...    compute cepa shifts
c
        call cepac(cr(khbc+1),cr)
c
c...    modify h (floc) and h**2 (f2loc) matrices
c
        call hc2cep(cc,ncpcc,floc,fsloc,f2loc,f2sloc,joop)
c
c...    return to the contracted MP universe
c
        call mpretr
        call gmpblk(mpcblk)
        call rdedx(cr,5*(nmpv+nmp1+nmp2+nmp2d),mpcblk,numci)
c
c...    get current c,z-vectors again
c
        call rdedx(cr(kmpc+1),ntotmp,ibestc,numci)
        call rdedx(cr(kmpz+1),ntotmp,ibestz,numci)
c
        ncycep = ncycep + 1
        go to 8025
      else
        ncycep = 0
      endif
ccepa
c
6666  continue
c
c     restart Davidson
c
      call cpuwal(top,bot)
      write(6,1001)icyc,top
1001  format(i8,f12.3,' space contracted')

      call wrt3(cr(kmpc+1),ntotmp,indexp(ic(1)),numci)
      call wrt3(cr(kmpz+1),ntotmp,indexp(iz(1)),numci)
      go to 9999
c
c     END OF DAVIDSON
c     get every thing right before exit.
c
777   continue
      if (tester.lt.thresh) then
c...     Hurrah, we're converged.
         write(6,101)icyc
101      format(//20x,26('*')//20x,'convergence at cycle',
     *          i6//20x,26('*'))
      else
c...     Boehoe, no convergence.
         write(6,102)
102      format(///20x,14('*')//20x,'no convergence'//20x,14('*'))
      endif
c
c...  Expand c-vector and store it in right format
c
      call rdedx(cr(kcc+1),ntotmp,ibestc,numci)
      if(iprimp.ge.1) call ccvpri(cr,cr(kcc+1),1.0d-3)
      call ccxpmp(cr(kc+1),cr(kcc+1),cr(kdc+1),cr(kvc+1),cr(kscr+1))
      call ciretr
      if(kc.gt.khbc) then
        do 200 i = 1, ntotal
          cr(khbc+i) = cr(kc+i)
200     continue
      else if(kc.lt.khbc) then
        do 210 i = 0, ntotal-1
          cr(khbc+ntotal-i) = cr(kc+ntotal-i)
210      continue
      endif
      call gciblk(mpcblk)
      call rdedx(cr,5*(nnst+nmin1+nmin2),mpcblk,numci)
      call dumpcz(cr,cr(khbc+1),index(26))
c
c...  Return to contracted MP universe
c
      call mpretr
      call gmpblk(mpcblk)
      call rdedx(cr,5*(nmpv+nmp1+nmp2+nmp2d),mpcblk,numci)
c
c...  Expand z-vector and store it in right format
c
      call rdedx(cr(kcc+1),ntotmp,ibestz,numci)
      call ccxpmp(cr(kc+1),cr(kcc+1),cr(kdc+1),cr(kvc+1),cr(kscr+1))
      call ciretr
      if(kc.gt.khbc) then
        do 220 i = 1, ntotal
          cr(khbc+i) = cr(kc+i)
220     continue
      else if(kc.lt.khbc) then
        do 230 i = 0, ntotal-1
          cr(khbc+ntotal-i) = cr(kc+ntotal-i)
230      continue
      endif
      call gciblk(mpcblk)
      call rdedx(cr,5*(nnst+nmin1+nmin2),mpcblk,numci)
      call dumpcz(cr,cr(khbc+1),index(27))
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine genzcc(cr,h11o,indc,indz,kcc,kc,kdc,kvc,kscr)
c
      implicit real*8 (a-h,o-z), integer(i-n)
c
      dimension cr(*)
c
c...  contracted matrix-vector multiplication
c...  (contracted genz)
c
c     disk -layout
c     R, HR (CSF)            : index(26),index(27)
c     R, HR (excited states) : indc, indz
c
c     kcc : memory address excited state vector
c     kc  : memory address CSF-vector
c
c     procedure:
c
c       uncontract excited state c-vector.
c       load CSF-configuration table.
c       CSF matrix-vector multiplication.
c       load excitated state configuration table.
c       contract CSF z-vector.
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
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
      common/mr_mp2l/ ivoc(62),nvoc,idoc(62),ndoc,kvcc(4),nvcc(4),
     *                ncmphn(10,18),ncmph1(10,18,8),ncmph2(10,18,8),
     *                ncmphd(10,18)
      common /mr_mp2con/ nsmpv,nsmp1,nsmp2,nsmp2d,ntotmp,
     *                   nmpv,nmp1,nmp2,nmp2d
      common/diisdc/ diiscr,ntdiis
c
c
      integer nref, nnst, nmin1, nmin2, mspin, nelec
      integer nstrt1, nlast1, nstrt2, nlast2, nmin12, nminv2, nminv1
      integer nlst21, nlastx
      common /ccntl/ nref,nnst,nmin1,nmin2,mspin,nelec,
     +               nstrt1,nlast1,nstrt2,nlast2,nmin12,nminv2,nminv1,
     +               nlst21,nlastx
c
c
      integer nval, ndoub, nsingl, ntripl, ntotal, ibig, jbig, kbig
      integer nvd, nvds, nsitri, lbig, mbig, nbig, nref0
      common /cntrl/ nval,ndoub,nsingl,ntripl,ntotal,ibig,jbig,kbig,
     +               nvd,nvds,nsitri,lbig,mbig,nbig,nref0
c
c
      integer nint, next, ninnex, nwm1, nw, nwp1, nintwo
      integer nint21, nintp1
      common /orb/ nint,next,ninnex,nwm1,nw,nwp1,nintwo,
     +             nint21,nintp1
c
c
      integer nxblk, index, numci, mblkci, ispbig, ispsma
      common/disktl/nxblk,index(199),numci,mblkci,ispbig,ispsma
      common/stoctl/khbill,khbz,khbc,khbb,khbzd,khbzs,khbzt,khbcd,
     *khbcs,khbct,kh,khb1,khb2,khb3,khb4,khb5,khb6,khb7,khb8
     *,khb44,khb55,khb66,khb88
c
      integer lcall,  mrstor
      logical modev, moded, modest, instor, loijka, msort, loijab
      logical loiajb, lvar
      common/diactl/modev,moded,modest,instor,
     +              loijka,msort,loijab,loiajb,lcall,mrstor,lvar
c
      logical malt
      real*8 potnuc,totnuc,e11,eigen,h11,z11,h22,s22,h12,s12,
     *     tester,eshift,thresh,aaa,bbb,scream,gprin,z12,z22,
     *     vdstpr,qester
      integer maxcyc,moddia,mprin,icyc
      common/erg/potnuc,totnuc,e11,eigen(nd200),h11
     *,z11,h22,s22,h12,s12,tester,eshift,thresh,aaa,bbb,scream
     *,gprin,z12,z22,vdstpr(3),qester,maxcyc,moddia,mprin,malt,icyc
c
c
      integer nirr,mult,ifrep,nirr8,jfrep
      common /table/ nirr,mult(8,8),ifrep,nirr8,jfrep
c
c
      integer nconst, ncond, maxpa, maxpat, maxsin, maxepa
      integer mxdvpa, mxddpa, nspips, nspipt, mxsspa, mxttpa
      integer mxstpa, nsitot, ntrtot
      integer mx12pa, mxvvpa, mxv2pa, maxpsq, max12
      common /auxh/ nconst(8),ncond(8),maxpa(4),maxpat,maxsin(3),
     +              maxepa(4),mxdvpa,mxddpa,nspips(8),nspipt(8),
     +              mxsspa,mxttpa,mxstpa,nsitot(8),ntrtot(8),
     +              mx12pa,mxvvpa,mxv2pa,maxpsq,max12
c
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
c
      integer nuse, ngot
      common /corctl/ nuse, ngot
c
c
c     integer maxord
c     parameter(maxord = 50)
c...  maxord = maximum order in perturbation theory
c     dimension indexp(10+2*maxord), ic(maxord), iz(maxord)
c     dimension epert(maxord)
c
c     dimension epair_mp(2,9)
c     logical diisc, dlev
c
c     parameter (nmodel=7)
c     character*39 descr(nmodel)
c     data descr /'      ** no projection applied **     ',
c    *            ' |ref>F<ref|+|sin>F<sin|+|doub>F<doub|',
c    *            '     |ref>F<ref| + |doub>F<doub|      ',
c    *            '  |ref>F<ref| + |sin+doub>F<sin+doub| ',
c    *            ' |ref>F<ref|+|sin>F<sin|+|rest>F<rest|',
c    *            '     |ref>F<ref| + |rest>F<rest|      ',
c    *            ' |ref>F<ref| + |doub+more>F<doub+more|'/
c
c     diisc = .false.
c
c...  expand the excitated-state vector to configuration-state vector
c
      call rdedx(cr(kcc+1),ntotmp,indc,numci)
      call ccxpmp(cr(kc+1),cr(kcc+1),cr(kdc+1),cr(kvc+1),cr(kscr+1))
c
c...  restore original information in the common blocks to enable
c...  the use of genz.
c
      call ciretr
      nuse = max(khbz,khbc,khbb)+ntotal
      call corfai('contract')
c
c...  move the csf-vector to proper position for genz.
c
      if(kc.gt.khbc) then
        do 80 i = 1, ntotal
          cr(khbc+i) = cr(kc+i)
80      continue
      else if(kc.lt.khbc) then
        do 90 i = 0, ntotal-1
          cr(khbc+ntotal-i) = cr(kc+ntotal-i)
90      continue
      endif
c
c...  retrieve config information to do genz
c
      call gciblk(mpcblk)
      call rdedx(cr,5*(nnst+nmin1+nmin2),mpcblk,numci)
c
c...  dump current vector to synchronize dumpfile contents with
c...  core contents
c
      call dumpcz(cr,cr(khbc+1),index(26))
c
c...  calculate  Z = H*C
c
      call genz(cr,cr,index(26),h11o)
      if(instor) call dumpcz(cr,cr(khbz+1),index(27))
c
c...  get c,z-vector in proper format
c
      call getcz(cr,cr(khbc+1),index(26))
      call wrt3(cr(khbc+1),ntotal,index(26),numci)
      call getcz(cr,cr(khbz+1),index(27))
      call wrt3(cr(khbz+1),ntotal,index(27),numci)
c
c...  return to the contracted MP universe
c
      call mpretr
      call gmpblk(mpcblk)
      call rdedx(cr,5*(nmpv+nmp1+nmp2+nmp2d),mpcblk,numci)
c
c...  retrieve, contract, and store zero-th order z-vector
c
      call rdedx(cr(kc+1),ntotal,index(27),numci)
      if (iprimp.ge.2000) write(6,603) (cr(kc+i),i=1,ntotal)
603   format(' - H*C (CSF) -',/,(2x,10e12.5))
      call ccntmp(cr(kc+1),cr(kcc+1),cr(kdc+1),cr(kvc+1),cr(kscr+1))
      if (iprimp.ge.2000) write(6,602) (cr(kcc+i),i=1,ntotmp)
602   format(' - H*C (excited states) -',/,(2x,10e12.5))
      call wrt3(cr(kcc+1),ntotmp,indz,numci)
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine ccvpri(conf,c,crit)
c
c...  print contracted ci-vector
c
      implicit real*8 (a-h,o-z)
      integer  conf
      dimension conf(5,*)
      dimension c(*)
c
chvd
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
      common/mr_mp2l/ ivoc(62),nvoc,idoc(62),ndoc,kvcc(4),nvcc(4),
     *                ncmphn(10,18),ncmph1(10,18,8),ncmph2(10,18,8),
     *                ncmphd(10,18)
      common /mr_mp2con/ nsmpv,nsmp1,nsmp2,nsmp2d,ntotmp,
     *                   nmpv,nmp1,nmp2,nmp2d
c
      integer nirr,mult,ifrep,nirr8,jfrep
      common /table/ nirr,mult(8,8),ifrep,nirr8,jfrep
c
c
      integer nint, next, ninnex, nwm1, nw, nwp1, nintwo
      integer nint21, nintp1
      common /orb/ nint,next,ninnex,nwm1,nw,nwp1,nintwo,
     +             nint21,nintp1
c
chvd
ccepa
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
      common/mr_cepctl/ kcpill,kcpv,kcpz,kcpx,kcpcd,kcpd,kcp2,kcpc2,
     1                  kcepc,kcepz,kcpcc,ncpcc,kcpbl,ncpbl,icpsto
     1                  ,kdavcz
      logical icpsto
ccepa
      common/mr_comcci/icci,iclass(3,3,3),istcc,iselcc
      character*12 texi(3)
      character*6 texj(4)
      character*2 ts
      data texi/'   reference','     single ','     double '/
      data texj/'vacuum',' n-1  ',' n-2  ',' n-2aa'/
c
      write(6,630)crit
630   format(////,' list of excited state coefficients greater than',
     1   e14.4,///,60('-'),/,
     1   ' state   coefficient   excitation class   range of states',
     1   /,60('-'))
      nconf = nmpv+nmp1+nmp2+nmp2d
      idumm = 1
      do 600 i = 1, nconf
         call upack2(conf(5,i),it,jj)
         jt = (jj-1)/3+1
         jj = jj - (jt-1)*3 - 1
         call upack2(conf(3,i),ntex,nsex)
         nex = max(ntex,nsex)
         call upack2(conf(4,i),idum,is)
         if (i.gt.nmpv+nmp1+nmp2) jt = 4
         if (jt.eq.1) nnn = nex
         if (jt.eq.2) nnn = nex*norbe(is)
         if (jt.eq.3) nnn = nex*norbtr(is)
         if (jt.eq.4) nnn = nex*next
         ts = ' '
         if (jt.eq.3.and.nsex.eq.0) ts = '_t'
         if (jt.eq.3.and.ntex.eq.0) ts = '_s'
         if (jt.eq.4) jt = 3
c        write(6,603) i,texi(it+1),texj(jt),ts,jj,is,nex,nnn
c    *               ,idumm,idumm+nnn-1
603      format(i6,2x,a12,1x,a6,a2,i6,2x,i6,4x,i6,4x,i6,
     *          2x,'(',i7,' -',i7,')')
         j=idumm
         write(6,610)j,c(j),texi(it+1),jj,jt-1,idumm,idumm+nnn-1
610      format(i6,2x,f12.8,a12,1x,'(',i1,',',i1,')   (',i7,' -',i7,')')
         do 700 j = idumm+1, idumm+nnn-1
           if(dabs(c(j)).lt.crit) go to 700
           write(6,620) j,c(j),texi(it+1),jj,jt-1
620        format(i6,2x,f12.8,a12,1x,'(',i1,',',i1,')')
700      continue
         idumm = idumm + nnn
600   continue
690   format(i5,f22.12,a)
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine cntrcc(conf,cr1,cr2,cr,ic1,ic12,sum12)
c
c...  contract a combination of c/z vectors to contributions
c...  per model configuration  for mrdcepa / i.e. all refs together
c...  note for n-2 contraction per spin-type
c     cr1  : vector 1   (iblk1 block-number if out-store)
c     cr2  : vector 2   (iblk2 block-number if out-store)
c     cr   : start of result-array for contractesd products
c     ic1  : index of array for the contracted <c/c>-vector (if>0)
c     ic12 : index of array to receive the contracted <c/z>-vector
c            length nnst-nref+1 + nmin1 + 2*nmin2
c     icpsto if .true. no reading needs to be done
c
      implicit real*8 (a-h,o-z)
      dimension cr1(*),cr2(*),cr(*)
      integer  conf
      dimension conf(5,*)
c
chvd
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
      common/mr_mp2l/ ivoc(62),nvoc,idoc(62),ndoc,kvcc(4),nvcc(4),
     *                ncmphn(10,18),ncmph1(10,18,8),ncmph2(10,18,8),
     *                ncmphd(10,18)
      common /mr_mp2con/ nsmpv,nsmp1,nsmp2,nsmp2d,ntotmp,
     *                   nmpv,nmp1,nmp2,nmp2d
c
      integer nirr,mult,ifrep,nirr8,jfrep
      common /table/ nirr,mult(8,8),ifrep,nirr8,jfrep
c
c
      integer nint, next, ninnex, nwm1, nw, nwp1, nintwo
      integer nint21, nintp1
      common /orb/ nint,next,ninnex,nwm1,nw,nwp1,nintwo,
     +             nint21,nintp1
c
chvd
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
      common/mr_cepctl/ kcpill,kcpv,kcpz,kcpx,kcpcd,kcpd,kcp2,kcpc2,
     1                  kcepc,kcepz,kcpcc,ncpcc,kcpbl,ncpbl,icpsto
     1                  ,kdavcz
      logical icpsto
c
      common/mr_comcci/icci,iclass(3,3,3),istcc,iselcc
chvd
      kc   = 0
      kc1  = (ic1 -1)*ncpcc + 1
      kc12 = (ic12-1)*ncpcc + 1
c
      call vclr(cr(kc12),1,ncpcc)
      if (ic1.gt.0) call vclr(cr(kc1),1,ncpcc)
c
      nconf = nmpv+nmp1+nmp2+nmp2d
      idumm = 1
      do 100 i = 1, nconf
         call upack2(conf(5,i),iexlvl,jj)
         jt = (jj-1)/3+1
         nhole = jj - (jt-1)*3 - 1
         call upack2(conf(3,i),ntex,nsex)
         nex = max(ntex,nsex)
         call upack2(conf(4,i),idum,is)
         if (i.gt.nmpv+nmp1+nmp2) jt = 4
         if (jt.eq.1) nnn = nex
         if (jt.eq.2) nnn = nex*norbe(is)
         if (jt.eq.3) nnn = nex*norbtr(is)
         if (jt.eq.4) nnn = nex*next
         if (jt.ne.4) then
           npart  = jt-1
         else
           npart  = 2
           iexlvl = 2
         endif
         ioffst = iclass(nhole+1,npart+1,iexlvl+1)
         cr(kc12+ioffst) = cr(kc12+ioffst) +
     *      ddot(nnn,cr1(idumm),1,cr2(idumm),1)
         if (ic1.gt.0) cr(kc1+ioffst) = cr(kc1+ioffst) +
     *      ddot(nnn,cr1(idumm),1,cr1(idumm),1)
         idumm = idumm + nnn
100   continue
c
c...  calculate the sum
c
      sum12 = 0.0d0
      do 120 i = 1, ncpcc
        sum12 = sum12 + cr(kc12-1+i)
120   continue
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine hc2cep(cc,nc,h,hs,h2,h2s,ndim)
c
c...  produce the modified <h> and <h**2> matrices by applying
c...  the (mrd)cepa diagonal shifts. use the saved original <hs>
c...  and <h2s> matrices and the contracted <c.c> and <c.z> vectors
c...  (from disc) / dimension of h-matrices is ndim
c...   normalisation-factors for the c-vectors are stored in xn
c
      implicit real*8 (a-h,o-z)
      dimension cc(nc,5)
      dimension h(*),hs(*),h2(*),h2s(*)
c
      common/mr_cepctl/ kcpill,kcpv,kcpz,kcpx,kcpcd,kcpd,kcp2,kcpc2,
     1                  kcepc,kcepz,kcpcc,ncpcc,kcpbl,ncpbl,icpsto
     1                  ,kdavcz
      logical icpsto
      integer maxpair
      parameter (maxpair=40)
      real*8 c0,h11ci,ecorr_mr epair
      common /mr_pairs/ c0,h11ci,ecorr_mr,epair(maxpair,maxpair)
c
      real*8 eshifc eshifd
      common /shifts/ eshifc(maxpair,maxpair),eshifd(maxpair)
      common/mr_comcci/icci,iclass(3,3,3),istcc,iselcc
c
      integer nxblk, index, numci, mblkci, ispbig, ispsma
      common/disktl/nxblk,index(199),numci,mblkci,ispbig,ispsma
c
c...  set up shifts and shift**2 arrays
c
c     call expdsh(iconf,cc(1,1))
c     call vvtv(ncpcc,cc(1,2),cc(1,1),cc(1,1))
      do 10 iexlvl = 0, 2
        do 20 nhole = 0, iexlvl
          do 30 npart = 0, iexlvl
            ioffst = iclass(nhole+1,npart+1,iexlvl+1)
c           cc(1+ioffst,1) = eshifc(npart+1,nhole+1)
            cc(1+ioffst,1) = eshifc(1+ioffst,1)
30        continue
20      continue
10    continue
c
      kh = 1
      kbl = kcpbl
c
c...  go over the h-matrix-triangle
c
      do 50 i=1,ndim
         do 40 j=1,i
c...  read contracted vectors  (ci.cj / ci.zj / zi.cj)
            call rdedx(cc(1,3),ncpcc*3,kbl,numci)
c...  h = hs + <ci/diag/cj>
            h(kh) = hs(kh) + ddot(ncpcc,cc(1,1),1,cc(1,3),1)
c...  hs = h2s + <ci/diag/zj> + <zi/diag/cj> + <ci/diag**2/cj>
            h2(kh) = h2s(kh)  + ( ddot(ncpcc,cc(1,1),1,cc(1,4),1)
     *                        +   ddot(ncpcc,cc(1,1),1,cc(1,5),1)
     *                        +   ddot(ncpcc,cc(1,2),1,cc(1,3),1) )
c
            kh = kh + 1
            kbl = kbl + ncpbl
c
40       continue
50    continue
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine ccpvds(conf,z,c)
c
c...  adjust z-vector for Davidson.
c
      implicit real*8 (a-h,o-z)
      integer  conf
      dimension conf(5,*)
      dimension c(*),z(*)
c
chvd
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
      common/mr_mp2l/ ivoc(62),nvoc,idoc(62),ndoc,kvcc(4),nvcc(4),
     *                ncmphn(10,18),ncmph1(10,18,8),ncmph2(10,18,8),
     *                ncmphd(10,18)
      common /mr_mp2con/ nsmpv,nsmp1,nsmp2,nsmp2d,ntotmp,
     *                   nmpv,nmp1,nmp2,nmp2d
c
      integer nirr,mult,ifrep,nirr8,jfrep
      common /table/ nirr,mult(8,8),ifrep,nirr8,jfrep
c
c
      integer nint, next, ninnex, nwm1, nw, nwp1, nintwo
      integer nint21, nintp1
      common /orb/ nint,next,ninnex,nwm1,nw,nwp1,nintwo,
     +             nint21,nintp1
c
chvd
ccepa
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
      common/mr_cepctl/ kcpill,kcpv,kcpz,kcpx,kcpcd,kcpd,kcp2,kcpc2,
     1                  kcepc,kcepz,kcpcc,ncpcc,kcpbl,ncpbl,icpsto
     1                  ,kdavcz
      logical icpsto
      integer maxpair
      parameter (maxpair=40)
      real*8 c0,h11ci,ecorr_mr epair
      common /mr_pairs/ c0,h11ci,ecorr_mr,epair(maxpair,maxpair)
c
      real*8 eshifc eshifd
      common /shifts/ eshifc(maxpair,maxpair),eshifd(maxpair)
ccepa
      common/mr_comcci/icci,iclass(3,3,3),istcc,iselcc
      nconf = nmpv+nmp1+nmp2+nmp2d
      idumm = 1
      do 100 i = 1, nconf
         call upack2(conf(5,i),iexlvl,jj)
         jt = (jj-1)/3+1
         nhole = jj - (jt-1)*3 - 1
         call upack2(conf(3,i),ntex,nsex)
         nex = max(ntex,nsex)
         call upack2(conf(4,i),idum,is)
         if (i.gt.nmpv+nmp1+nmp2) jt = 4
         if (jt.eq.1) nnn = nex
         if (jt.eq.2) nnn = nex*norbe(is)
         if (jt.eq.3) nnn = nex*norbtr(is)
         if (jt.eq.4) nnn = nex*next
         if (jt.ne.4) then
           npart = jt-1
         else
           npart = 2
         endif
         icl = iclass(nhole+1,npart+1,iexlvl+1)+1
         call daxpy(nnn,eshifc(icl,1),c(idumm),1,z(idumm),1)
         idumm = idumm + nnn
100   continue
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine vvpri(cr,kc,kcc,kdc,kvc,kscr,ntotal,ntotmp,crit)
      implicit real*8 (a-h,o-z), integer (i-n)
c
c     This subroutine prints the transformations matrices from
c     csf-space to excited state-space and vice versa. The printed
c     transformation is the transformation that is actually used.
c
c     The printed coefficients are obtained from transforming a
c     unit-vector in one space to a vector in the other space.
c
c     Note that the subroutine is rather inefficient and is meant for
c     debugging purposes only.
c
      dimension cr(*)
c
      crit = dabs(crit)
      write(*,*)'*** PRINTING TRANSFORMATION MATRIX'
      write(*,*)'*** from CSF to EXCITED STATES    '
      write(*,*)
      write(*,*)'     csf  excited  coefficient'
      write(*,*)'           state              '
      write(*,*)
      call vclr(cr(kc+1),1,ntotal)
      do i = 1, ntotal
        if (i.gt.1) cr(kc+i-1) = 0.0d0
        cr(kc+i) = 1.0d0
        call ccntmp(cr(kc+1),cr(kcc+1),cr(kdc+1),cr(kvc+1),cr(kscr+1))
        sum = 0.0d0
        do j = 1, ntotmp
          if(dabs(cr(kcc+j)).gt.crit) then
            write(*,10)i,j,cr(kcc+j)
10          format(1x,i7,2x,i7,2x,e12.5)
            sum = sum + dabs(cr(kcc+j))
          endif
        enddo
        if (sum.gt.0.0d0) then
          write(*,*)
        endif
      enddo
      write(*,*)
      write(*,*)'*** from EXCITED STATES to CSF   '
      write(*,*)
      write(*,*)' excited      csf  coefficient'
      write(*,*)'  state                       '
      write(*,*)
      call vclr(cr(kcc+1),1,ntotmp)
      do i = 1, ntotmp
        if (i.gt.1) cr(kcc+i-1) = 0.0d0
        cr(kcc+i) = 1.0d0
        call ccxpmp(cr(kc+1),cr(kcc+1),cr(kdc+1),cr(kvc+1),cr(kscr+1))
        sum = 0.0d0
        do j = 1, ntotal
          if(dabs(cr(kc+j)).gt.crit) then
            write(*,10)i,j,cr(kc+j)
            sum = sum + dabs(cr(kc+j))
          endif
        enddo
        if (sum.gt.0.0d0) then
          write(*,*)
        endif
      enddo
      write(*,*)
      write(*,*)'*** END OF TRANSFORMATION MATRIX'
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine sttria
      implicit real*8 (a-h,o-z)
c
c     Sets the parameters for the trial vector as requested for the
c     contracted CI davidson. The trial vector for the contracted CI
c     is not necessarily the same as the reference vector for the
c     construction of the CI space.
c
c     Note that it is not clear what the relation is between a state
c     calculated in the contracted space of some other state and
c     and physical reality. But one may want to do these kind of
c     calculations for debugging purposes.
c
c
      real*8 tvec
      integer ivec, nvec, mxvec, iselvc, nselvc, istvc
      integer iprvc
      common /trial/ tvec(20),ivec(20),nvec,mxvec,
     +               iselvc,nselvc,istvc,iprvc
c
      common/mr_comcci/icci,iclass(3,3,3),istcc,iselcc
c
      istvc  = istcc
      iselvc = iselcc
c
      return
      end
      subroutine corfai(text)
c.... core-failure imitation (extended)
      implicit real*8  (a-h,o-z),integer (i-n)
      character*(*) text
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer nuse, ngot
      common /corctl/ nuse, ngot
c
      parameter (mseg=30)
      character*10 frag(mseg)
      integer iuse(mseg)
c
      save maxuse,nseg,frag,iuse
      data maxuse/0/,nseg/0/,iuse/mseg*0/
c
      ll = 1
      if (text.ne.' ') then
         ll = min(len(text),10)
         do 10 i=1,nseg
            if (frag(i)(1:ll).eq.text(1:ll)) then
               iuse(i) = max(iuse(i),nuse)
               go to 20
            end if
10       continue
         nseg = nseg + 1
         frag(nseg) = text(1:ll)
         iuse(nseg) = nuse
      end if
c
20    maxuse=max(maxuse,nuse)
c
      if (nuse.le.ngot) return
      write(iwr,1) text(1:ll),nuse,ngot
1     format(
     + '*****************************************************'/
     + '************* ',a,'   *******************************'/
     + '*** core requested ',i10,' (words) *******************'/
     + '*** core available ',i10,' (words) *******************'/
     + '*****************************************************'/)
      intry = 0
      go to 100
c
      entry corpr
c
      intry = 1
c
c...  print maximum core used sofar
c
100   write(iwr,601) maxuse
      do 50 i=1,nseg
50    write(iwr,602) frag(i),iuse(i)
      write(iwr,603)
601   format(/' ****************************************************',
     1       /' ***** max core used'       ,i20,      ' words  *****')
602   format( '   *****  ',a10,i20,                 '      *****')
603   format( ' ****************************************************')
c
      if (intry.eq.0) call caserr('Insufficient memory in Direct CI')
c
      return
      end
      subroutine dumpcc(cc,i,j)
c
c...  dump the 3 contracted vectors  for h(i,j)
c...  if (i.eq.j) the last two are identical / see to that
c...  ***   ccepa  ***
c
      implicit real*8 ( a-h,o-z)
      dimension cc(*)
      common/mr_cepctl/ kcpill,kcpv,kcpz,kcpx,kcpcd,kcpd,kcp2,kcpc2,
     1                  kcepc,kcepz,kcpcc,ncpcc,kcpbl,ncpbl,icpsto
     1                  ,kdavcz
      logical icpsto
c
      integer nxblk, index, numci, mblkci, ispbig, ispsma
      common/disktl/nxblk,index(199),numci,mblkci,ispbig,ispsma
c
      kbl = kcpbl + ( i*(i-1)/2 + j - 1)*ncpbl
c       copy ci.zj to zi.cj   for i=j
      if (i.eq.j) call dcopy(ncpcc,cc(ncpcc+1),1,cc(ncpcc*2+1),1)
c
      call wrt3(cc,ncpcc*3,kbl,numci)
c
      return
      end
      subroutine inicepa(cr,iconf)
c
c...  routine to start off cepa.
c...  for now really only psi(0) is calculated for e0
c
      implicit real*8 (a-h,o-z), integer (i-n)
      integer iconf
c
c
      real*8 critcep, crycep
      integer ifcepa, icepa , lcepaul, lsinsh, itcepa, lprcep, indcep
      integer nrfcsf, mcycep, ibrst,icorr_mr,mrdoc
      common/ccepa/critcep,crycep,ifcepa,icepa,lcepaul,lsinsh,itcepa,
     1                lprcep,indcep,nrfcsf,mcycep,ibrst,icorr_mr,mrdoc
c
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer nuse, ngot
      common /corctl/ nuse, ngot
c
c
      real*8 tvec
      integer ivec, nvec, mxvec, iselvc, nselvc, istvc
      integer iprvc
      common /trial/ tvec(20),ivec(20),nvec,mxvec,
     +               iselvc,nselvc,istvc,iprvc
c
      common/stoctl/khbill,khbz,khbc
      dimension cr(*),iconf(*)
c
      iselvc_save = iselvc
      nselvc_save = nselvc
      iselvc = 1
      nselvc = nrfcsf
      call setc1d(cr(khbill+1),iconf,e0ref,ngot-khbill,iwr)
      iselvc = iselvc_save
      nselvc =  nselvc_save 
c
      return
      end
      subroutine cepac(cr,conf)
c
c...  compute cepa shifts
c...
c...   closed shell  cepa 0,1,2 :
c...   formulas from c.zirz and r.ahlrichs
c...                 recent developments in coupled pair theories
c...                 proc. of daresbury weekend on electron correlation
c...                 report dl/sci/r14 (theory) november 1979
c...                 ed. by m.f.guest, s.wilson
c...
c...   mrd-cepa (0) version for closed/open shell reference functions
c...   mrcepa(1) and mrcepa-1 versions 
c...
c...   please note
c...   no changes are made to c or z vectors on disk
c...   so they are changed whenever they are used
c...   cepa- changes are surrounded by ccepa  comment cards
c...   iconf(5,...) is used to store excitation information
c...   ** if this info changes the following routines are affected
c...   **  cepinf,expdsh,cepadvv,cepadst,cepaddd,cepvds  **
ccepa  affected subroutines are :
ccepa  a. specially written for cepa :
ccepa     cepa,cepin,cepair,cciaib,cciajb,ccijkl,ccijka,cepinf,
ccepa     cepadvv,cepaddd,cepadst,cepvds,prcepa,cepijkl,cntrci,expdsh,
ccepa     dumpcc,hh2cep
ccepa  b. original ci routines
ccepa     block data,inpro,pass,direct,vvc2,ddc2,ssc2,hbill,
ccepa      diabas,davids
cci
cci    Contracted cepa:
cci    a. new subroutine involved are:
cci       ccepin, davidc, genzcc, cntrcc, hc2cep, ccpvds.
cci    b. modified subroutines are:
cci       cepac.
c
      implicit real*8 (a-h,o-z)
      integer conf
cci
      common/mr_comcci/icci,iclass(3,3,3),istcc,iselcc
cci
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
      real*8 critcep, crycep
      integer ifcepa, icepa , lcepaul, lsinsh, itcepa, lprcep, indcep
      integer nrfcsf, mcycep, ibrst,icorr_mr,mrdoc
      common/ccepa/critcep,crycep,ifcepa,icepa,lcepaul,lsinsh,itcepa,
     1                lprcep,indcep,nrfcsf,mcycep,ibrst,icorr_mr,mrdoc
c
c
      integer nref, nnst, nmin1, nmin2, mspin, nelec
      integer nstrt1, nlast1, nstrt2, nlast2, nmin12, nminv2, nminv1
      integer nlst21, nlastx
      common /ccntl/ nref,nnst,nmin1,nmin2,mspin,nelec,
     +               nstrt1,nlast1,nstrt2,nlast2,nmin12,nminv2,nminv1,
     +               nlst21,nlastx
c
      integer maxpair
      parameter (maxpair=40)
      real*8 c0,h11ci,ecorr_mr epair
      common /mr_pairs/ c0,h11ci,ecorr_mr,epair(maxpair,maxpair)
c
      real*8 eshifc eshifd
      common /shifts/ eshifc(maxpair,maxpair),eshifd(maxpair)
c...
c...  epair contains pair-energies :
c...      epair(i,j)  i >= j  : singlet  / i  < j  triplet for cepa(0-3)
c...      for MRCEPA's no singlet triplet division is made
c...      adddressing completely determined by conf(5,i)
c...  eshifc like epair
c..   eshifd not used in MRCEPA variants
c...
c
      integer nint, next, ninnex, nwm1, nw, nwp1, nintwo
      integer nint21, nintp1
      common /orb/ nint,next,ninnex,nwm1,nw,nwp1,nintwo,
     +             nint21,nintp1
c
      common/mr_mp2l/ ivoc(62),nvoc,idoc(62),ndoc,kvcc(4),nvcc(4),
     *                ncmphn(10,18),ncmph1(10,18,8),ncmph2(10,18,8),
     *                ncmphd(10,18)
      common/mr_mp2/mp,mpmod,indmp(6),isec_mp,nummp,iblk_mp,iprimp,
     *           krefmp,maxcmp,critmp,e0mp,e0ref,iwnmp(2,6),
     *           plevel(2),levelp,ifdiaf,fignor,cortho,iortho
c
c
      integer nval, ndoub, nsingl, ntripl, ntotal, ibig, jbig, kbig
      integer nvd, nvds, nsitri, lbig, mbig, nbig, nref0
      common /cntrl/ nval,ndoub,nsingl,ntripl,ntotal,ibig,jbig,kbig,
     +               nvd,nvds,nsitri,lbig,mbig,nbig,nref0
c
c
      integer nxblk, index, numci, mblkci, ispbig, ispsma
      common/disktl/nxblk,index(199),numci,mblkci,ispbig,ispsma
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      logical malt
      real*8 potnuc,totnuc,e11,eigen,h11,z11,h22,s22,h12,s12,
     *     tester,eshift,thresh,aaa,bbb,scream,gprin,z12,z22,
     *     vdstpr,qester
      integer maxcyc,moddia,mprin,icyc
      common/erg/potnuc,totnuc,e11,eigen(nd200),h11
     *,z11,h22,s22,h12,s12,tester,eshift,thresh,aaa,bbb,scream
     *,gprin,z12,z22,vdstpr(3),qester,maxcyc,moddia,mprin,malt,icyc
      dimension conf(5,*),cr(*)
      real*8 Ai,Bi,Ci,Aa,Ba,Ca,ph_pair(3,3)
      integer Ni,Na
c
      character*6 txt(3)
      data txt/'   ref','single','double'/
c...  correlation energy
      call vclr(eshifc,1,maxpair*maxpair)
      call vclr(eshifd,1,maxpair)
c
c... select cepa variant
c
      if (icepa.gt.9) go to 1001
c
      ecorr = 0.0d0
      do 10 i=1,nint
         do 10 j=1,nint
10    ecorr = ecorr + epair(i,j)
c
      go to (100,200,300),icepa+1
c...  cepa0
100   do 110 i=nstrt2,nlast2
         call upack2(conf(5,i),ia,ib)
         eshifc(ia,ib) = ecorr
         eshifc(ib,ia) = ecorr
110   continue
      if (lsinsh.eq.0) go to 999
c...  doublets (according to ahlrichs not shifted in cepa0)
      do 120 i=nstrt1,nlast1
         call upack2(conf(5,i),ia,ib)
120   eshifd(ib) = ecorr
      go to 999
c
c...  cepa1
c
200   do 260 i=nstrt2,nlast2
         call upack2(conf(5,i),ia,ib)
         ee = 0.0d0
         do 210 k=1,nint
210      ee = ee + epair(ia,k) + epair(ib,k) + epair(k,ia) + epair(k,ib)
         eshifc(ia,ib) = ecorr - ee/4.0d0
         if (ia.ne.ib) eshifc(ib,ia) = eshifc(ia,ib)
c
260   continue
      if (lsinsh.eq.0) go to 999
c...  doublets
      do 290 i=nstrt1,nlast1
         call upack2(conf(5,i),ia,ib)
         ee = 0.0d0
         do 270 k=1,nint
270      ee = ee + epair(ib,k) + epair(k,ib)
290   eshifd(ib) = ecorr - ee/2.0d0
c
      go to 999
c
c...  cepa2
c
300   do 350 i=nstrt2,nlast2
         call upack2(conf(5,i),ia,ib)
         if (lcepaul.ne.0) go to 330
c...  cepa2  a variant  singlet epv = singlet pair etc.
            eshifc(ia,ib) = ecorr - epair(ia,ib)
            if (ia.ne.ib) eshifc(ib,ia) = ecorr - epair(ib,ia)
         go to 350
c...  cepa2  b variant  singlet epv = triplet pair
330         eshifc(ia,ib) = ecorr - epair(ib,ia)
            if (ia.ne.ib) eshifc(ib,ia) = ecorr - epair(ia,ib)
350   continue
c
      if (lsinsh.eq.0) go to 999
c...  doublets   (according to ahlrichs not shifted in cepa2)
      do 390 i=nstrt1,nlast1
         call upack2(conf(5,i),ia,ib)
         ee = 0.0d0
         do 370 k=1,nint
370      ee = dmin1(ee,(epair(ib,k)+epair(k,ib))/2.0d0)
390   eshifd(ib) = ecorr - ee
c
999   if (lprcep.le.0) return
c
      write(iwr,1) ecorr,c0,h11ci,e0ref
      if (lsinsh.gt.0) write(iwr,2) (eshifd(i),i=1,nint)
      write(iwr,3)
      do 998 i=1,nint
998   write(iwr,6) (eshifc(j,i),j=1,nint)
      write(iwr,5)
      do 1000 i=1,nint
1000  write(iwr,6) (epair(j,i),j=1,nint)
      write(iwr,7)
c
      return
c
c********************************************************************
c
c        ***********    mrd  cepa    ***********
c
c     calculate the shifts per excitation class
c     as characterized by the hole and particle numbers
c
c     -> mrdcepa(0)
c     P.J.A. Ruttink, J.H. van Lenthe, R.Zwaans, G.C.Groenenboom
c     J.Chem.Phys. 94, 7212 (1991)
c
c     ecorr = the total correlation energy - hf-like terms
c     shift(j,i) = ecorr - vi-terms
c     i = # of inactive orbital holes  + 1
c     j = # of external orbital particles   + 1 
c     **note** The compatibility with closed shell cepa has been abandoned
c              epair(j,i) contains 'singlet+triplet' pair-energy 
c              shifts are treated similarly  (singlet = triplet)
c              **note** doublet shifts not any more in seperate array
c
c     the vi-terms are given by :
c
c     nh / np    0             1              2
c
c     0        ecort     10/20/01/11/21     10/20
c
c     1   10/01/11/02/12    10/01/11         10
c
c     2        01/02           01         no vi-terms
c
c     the (00) contribution to ecort should vanish
c     since the ref conf set is taken to be complete
c     within the active orbital subset
c     if it does not,it should be neglected in the calculation of epair
c     the vi-terms are the pair energies to be subtracted from the shifts
c
c     add singlet triplet contributions
c
c=======================================================================
c
c        *********** MR-CEPA variants   ***********
c        Now the shifts themselves are given in the table
c
c     The vi-epv table is given by : (these are contributions to shifts)
c     shift(k1,l1) = table(l1,k1;l2,k2) * epair(k2,l2)
c     the epairs and shifts (as in program) are indicated as ( , )
c     shifts : (1,..) = vacuum, (2,..) = eshifd() , (3,..) = doublet
c
c     (k2,l2)    00    10    20    01    11    21    02    12    22
c     (k1,l1)   (1,1) (1,2) (1,3) (2,1) (2,2) (2,3)  (3,1) (3,2) (3,3)
c     00  (1,1)   0     0     0     0     0     0     0     0     0
c     10  (1,2)   0     0     Bi    0     0     Bi    0     0     Bi
c     20  (1,3)   0     Bi    Ai    0     Bi    Ai    0     Bi    Ai
c     01  (2,1)   0     0     0     0     0     0     Aa    Ba    1
c     11  (2,2)   0     0     Bi    0     0     Bi    Ba    Ci*Ca Bi
c     21  (2,3)   0     Bi    Ai    0     Bi    Ai    1     Bi    Ai
c     02  (3,1)   0     0     0     Aa    Ba    1     Aa    Ba    1
c     12  (3,2)   0     0     Bi    Ba    Ci*Ca Bi    Ba    Ci*Ca Bi
c     22  (3,3)   0     Bi    Ai    1     Bi    Ai    1     Bi    Ai
c
c A = (n-2)*(n-3)(n*(n-1)) = 1-2(2n-3)/(n(n-1)) 
c B = (n-2)/n = 1-2/n 
c C = (n-1)/n = 1-1/n
c Ai,Bi,Ci refer to inactive (doubly occupied orbitals)
c ni = # electrons in the doubly occupied orbitals (= 2*ndoc)
c Aa,Ba,Ca are for the active (variably occupied orbitals)
c na = # electrons in active space
c na+ni = n (total # electrons)
c no EPV terms : Ai = Bi = Ci = Aa = Ba = Ca = 1 (=> MR CEPA)
c
c Contribution to shift for class (k1,l1) = (table entry)*epair(k2,l2)
c k = # of holes, l = # of particles
c
c Versions :
c MRCEPA0(icepa=100): (linear) shift = ecorr
c    CEPA MRD0
c AQCC   (icepa=101): (averaged EPV acc. to Szalay) : shift = A*ecorr
c    CEPA AQCC      (P.G.Szalay,R.J.Bartlett, CPL 214,481(1993))
c ACPF   (icepa=102): (averaged EPV according to Ahlrichs) : shift = B*ecorr
c    CEPA ACPF      (R.Gdanitz,R.Ahlrichs, CPL 143,413(1988))
c MRCEPA (icepa= 10): (only VI) : Ai=Bi=Ci=Aa=Ba=Ca=1
c    CEPA MRD       (P.J.A. Ruttink,J.H.van Lenthe,R.Zwaans,G.C.Groenenboom,
c                    J.Chem.Phys.94,7212(1991))
c MRCEPA(1) (icepa= 11): (VI/EPV for inactive(doc)) : Aa=Ba=Ca=1 
c    CEPA MR(1 DOC Cf. L.Fusti-Molnar,P.G.Szalay,J.Phys.Chem.100,6288(1996)
c MRCEPA(1)(icepa= 12): (VI/EPV for all) : as in Table 
c    CEPA MR(1 ALL P.J.A.Ruttink et al, under investigation (2002) 
c MRCEPA-1 (icepa= 20): uses doc holes tp make mrceopa(1) better  
c    CEPA MR-1 ALL P.J.A.Ruttink et al, under investigation (2002) 
c
c In addition the way the correlation energy is calculated may be specified
c for the first 3 variants. Choices are PROJECT (icorr_mr=1) 
c (the expectation value of c projected 'CI' function on the reference space, 
c or VARIA (icorr_mr=2), the energy
c of the variational optimised reference space. 
c Default is VA IA, (Ahlrichs choice)
c
c
c      ADDRESSING ......
c      for  all MRCEPA's the distinction between singlet/triplet is gone (also in cepair)
c
c      **mrd**  cepa(0)
c     i,j are :   1 / # holes + 1 in doc for vacuum
c                 2 / # holes + 1 in doc for nmin1
c                 3 / # holes + 1 in doc for nmin2  (singlet/triplet)
c      **mrd** cepa-1 
c     np = # partcles ; nh = # holes ; ib first doc hole (or 0) ; ia second doc hole (or 0)
c     mrdoc  = # doubly occupied orbitals (in reference function)
c     i,j are :   ip / jp 
c                 all                 jp = 1 + ib
c                 np=0,1,2,nh=0,1     ip = np+1
c                 np=0,1,2,nh=2       ip = 3+mrdoc*np+ia 
c     ndoc numbering is relative (done by routine docfind)
c
1001  if(icci.ne.0) goto 1100
c
      ecorr = 0.0d0
      if (icepa.ge.20.and.icepa.lt.30) then
c...     translate to mrdcepa(0) pair energies
         do i=1,3
            do j=1,3
               ph_pair(j,i) = 0.0d0
            end do
         end do
         do i=1,3
           ph_pair(i,1) = epair(i,1)
           do j=1,mrdoc
             ph_pair(i,2) = ph_pair(i,2) + epair(i,j+1)
             do k=1,mrdoc
               ph_pair(i,3) = ph_pair(i,3) + epair(3+mrdoc*(i-1)+k,j+1)
             end do
           end do
           ecorr = ecorr + ph_pair(i,1) + ph_pair(i,2) + ph_pair(i,3)
         end do
      else
c
        do i=1,3
           ecorr = ecorr + epair(3,i) + epair(2,i) + epair(1,i)
        end do
      end if
c
      e0_used = 0.0d0
c
      if (icorr_mr.eq.2) then
         ecorr = h11 - e0ref
         e0_used = e0ref
      end if
c
      if (icorr_mr.eq.1) then
c
c...   calculated ref-part expectation
c
         kc = 1
         kz = kc + nrfcsf
         kw = kz + nrfcsf
         call rdedx(cr(kc),nval,index(26),numci)
         call vclr(cr(kz),1,nrfcsf)
         call vrijkl(index,nrfcsf,cr(kz),cr(kc),cr(kw),conf)
         e0_used = ddot(nrfcsf,cr(kc),1,cr(kz),1) /
     1             ddot(nrfcsf,cr(kc),1,cr(kc),1)
         ecorr = h11 - e0_used 
      end if
      ecorr_mr = ecorr
c
c...   first the simple models
c
      if(icepa.ge.100) then
        if (icepa.eq.100) then
c...       Plain MRDCEPA 0 
           cepa_shift =ecorr
        else if (icepa.eq.101) then
c...          AQCC
           if (nelec.eq.1) call caserr('aqcc with 1 electron')
           A = (nelec-2.0d0)*(nelec-3.0d0)/(nelec*(nelec-1.0d0))
           cepa_shift = ecorr*A
        else if (icepa.eq.102) then
c...          ACPF
           B = (nelec-2.0d0)/nelec
           cepa_shift = ecorr*B
        end if
c...    Shift all non-reference states with (scaled) correlation energy
        do i = 1, 3
          do j = 1, 3
            eshifc(i,j) = cepa_shift
          enddo
        enddo
        eshifc(1,1) = 0.0d0
        if (lsinsh.le.0) then
           do j=1,3
              eshifc(2,j) = 0.0d0
           end do
        end if
c
      else if (icepa.eq.10.or.icepa.eq.11.or.icepa.eq.12) then
c
         Ni = 2 * ndoc
         Na = nelec - Ni
c
            Ai = 1.0d0
            Bi = 1.0d0
            Ci = 1.0d0
            Aa = 1.0d0
            Ba = 1.0d0
            Ca = 1.0d0
c        
         if (icepa.eq.10) then
c...     MRCEPA - no EPV
         end if
         if ((icepa.eq.11.or.icepa.eq.12).and.Ni.gt.0) then
c...     MR-CEPA(1) (inactive only)
            Ai =  1.0d0-2.0d0*(2.0d0*Ni-3.0d0)/(Ni*(Ni-1))
            Bi =  1.0d0-2.0d0/Ni
            Ci =  1.0d0-1.0d0/Ni
         end if
         if (icepa.eq.12) then
c...     MR-CEPA(1) (all)
            if (Na.le.1) then
               Aa = 0.0d0
               Ba = 0.0d0
               Ca = 0.0d0
            else
               Aa =  1.0d0-2.0d0*(2.0d0*Na-3.0d0)/(Na*(Na-1))
               Ba =  1.0d0-2.0d0/Na
               Ca =  1.0d0-1.0d0/Na
            end if
         end if
c
c...   0 particles  0,1,2 holes       vacuum
cc      eshifc(1,1) = 0.0d0
cc      eshifc(1,2) = ecorr - epair(2,1) - epair(3,1)
cc     1                    - epair(1,2) - epair(2,2)
cc     2                    - epair(3,2)
cc      eshifc(1,3) = ecorr - epair(2,1) - epair(3,1)
c...   1 particle  0,1,2 holes       doublet
cc      eshifd(1) = ecorr - epair(2,1) - epair(1,2)
cc     1                  - epair(2,2) - epair(1,3)
cc     2                  - epair(2,3)
cc      eshifd(2) = ecorr - epair(2,1) - epair(1,2) - epair(2,2)
cc      eshifd(3) = ecorr - epair(2,1)
cc      if (lsinsh.le.0) call vclr(eshifd,1,3)
ccc...   2 particle  0,1,2 holes       singlet/triplet
cc      eshifc(3,1) = ecorr - epair(1,2) - epair(1,3)
cc      eshifc(3,2) = ecorr - epair(1,2)
cc      eshifc(3,3) = ecorr
cc      do 1010 i=1,3
cc         eshifc(4,i) = eshifc(3,i)
cc         eshifc(i,4) = eshifc(3,i)
cc         eshifc(2,i) = eshifd(i)
cc1010  continue
c
c....   combined mrcepa mrcepa1
c....   **note** previously corrections, now contributions
c     
c         0 particles 0,1,2,holes  vacuum
        eshifc(1,1) = 0.0d0
        eshifc(1,2) = Bi * (epair(1,3) + epair(2,3) + epair(3,3))
        eshifc(1,3) = Bi * (epair(1,2) + epair(2,2) + epair(3,2))
     +              + Ai * (epair(1,3) + epair(2,3) + epair(3,3))
c...   1 particle  0,1,2 holes     doublet (was eshifd)
        eshifc(2,1) = Aa *  epair(3,1) + Ba * epair(3,2) + epair(3,3)
        eshifc(2,2) = Bi * epair(1,3)  + Bi * epair(2,3) 
     +              + Ba * epair(3,1)  + Ci*Ca * epair(3,2)
     +              + Bi * epair(3,3)
        eshifc(2,3) = Bi * epair(1,2) + Ai * epair(1,3)
     +              + Bi * epair(2,2) + Ai * epair(2,3)
     +              + epair(3,1)     + Bi * epair(3,2) + Ai * epair(3,3)
c
        if (lsinsh.le.0) then
           do j=1,3
              eshifc(2,j) = 0.0d0
           end do
        end if
c
c...   2 particle  0,1,2 holes       singlet/triplet
        eshifc(3,1) =  Aa * epair(2,1) + Ba * epair(2,2)
     +              +       epair(2,3) + Aa * epair(3,1)
     +              +  Ba * epair(3,2) +      epair(3,3)
        eshifc(3,2) =  Bi * epair(1,3) + Ba * epair(2,1)
     +              +  Ci*Ca * epair(2,2) + Bi * epair(2,3)
     +              +  Ba * epair(3,1) + Ci*Ca * epair(3,2)
     +              +  Bi * epair(3,3)
        eshifc(3,3) =  Bi * epair(1,2) + Ai * epair(1,3)
     +              +       epair(2,1) + Bi * epair(2,2)
     +              +  Ai * epair(2,3) +      epair(3,1)
     +              +  Bi * epair(3,2) + Ai * epair(3,3)
c
      else if (icepa.eq.20) then
c...   MRCEPA-1  --- try to get EPV as correctly as possible
         Ni = 2 * ndoc
         Na = nelec - Ni
c
c...     if Ni or Na are 0, make factors 1.0 
c
         if (Ni.eq.0) then
            Ai = 1.0d0
            Bi = 1.0d0
            Ci = 1.0d0
         else
c           Ai =  1.0d0-2.0d0*(2.0d0*Ni-3.0d0)/(Ni*(Ni-1))
c           Bi =  1.0d0-2.0d0/Ni
            Ci =  1.0d0-1.0d0/Ni
         end if
         if (Na.eq.0) then
            Aa = 1.0d0
            Ba = 1.0d0
            Ca = 1.0d0
         else
            Aa =  1.0d0-2.0d0*(2.0d0*Na-3.0d0)/(Na*(Na-1))
            Ba =  1.0d0-2.0d0/Na
            Ca =  1.0d0-1.0d0/Na
         end if
c
c...     2-hole pair-energy => 2-hole shift
c
         do ia=1,mrdoc
           do ib=1,ia
             eshifcc = 0.0d0
             do np=0,2
               corr = 0.0d0
               do k=1,mrdoc
                  corr = corr + epair(3+np*mrdoc+k,1+ia)
     1                        + epair(3+np*mrdoc+k,1+ib)
               end do
                  corr = corr + epair(3+np*mrdoc+ia,1+ia)
     1                        + epair(3+np*mrdoc+ib,1+ib)
               eshifcc = eshifcc + ph_pair(np+1,3) - corr/4.0d0 
             end do
             eshifc(3+mrdoc*2+ia,1+ib) = eshifcc
             eshifc(3+mrdoc  +ia,1+ib) = eshifcc
             eshifc(3        +ia,1+ib) = eshifcc
           end do
         end do
c
c...     2-hole pair-energy => 1-hole shift
c...     + 0,1 hole => 1-hole shift
c
         do ib=1,mrdoc
           eshifcc = 0.0
           do np=0,2
             corr = 0.0d0
             do k=1,mrdoc
               corr = corr + epair(3+np*mrdoc+k,1+ib)
             end do
             eshifcc = eshifcc + ph_pair(np+1,3) - corr/2.0d0 
           end do
           eshifc(3,1+ib) = eshifcc
     +                    +  Ba *    (ph_pair(2,1) + ph_pair(3,1))
     +                    +  Ci*Ca * (ph_pair(2,2) + ph_pair(3,2))
           eshifc(2,1+ib) = eshifcc
     +                    + Ba * ph_pair(3,1)  + Ci*Ca * ph_pair(3,2)
           eshifc(1,1+ib) = eshifcc
         end do
c
c...     1-hole pair-energy => 2-hole shift (2-hole contribution done above)
c...     + 0-hole contributions 
c
         do ia=1,mrdoc
           do ib=1,ia
             eshifcc = 0.0d0
             do np=0,2
               corr = epair(np+1,1+ia) + epair(np+1,1+ib)
               eshifcc = eshifcc + ph_pair(np+1,2) - corr/2.0d0 
             end do
             eshifc(3+mrdoc*2+ia,1+ib)=eshifc(3+mrdoc*2+ia,1+ib)+eshifcc
     1                                + ph_pair(3,1) + ph_pair(2,1)
             eshifc(3+mrdoc  +ia,1+ib)=eshifc(3+mrdoc  +ia,1+ib)+eshifcc
     1                                + ph_pair(3,1)
             eshifc(3        +ia,1+ib)=eshifc(3        +ia,1+ib)+eshifcc
           end do
         end do
c
c...    all contributions to 0-hole shifts (like MRCEPA(1))
c
         eshifc(1,1) =  0.0d0
         eshifc(2,1) =  Aa*ph_pair(3,1) + Ba*ph_pair(3,2) + ph_pair(3,3)
         eshifc(3,1) =  Aa * ph_pair(2,1) + Ba * ph_pair(2,2)
     +               +       ph_pair(2,3) + Aa * ph_pair(3,1)
     +               +  Ba * ph_pair(3,2) +      ph_pair(3,3)
c
        if (lsinsh.le.0) 
     1              call caserr('no singles not defined for MRCEPA-1')
c
      else
         call caserr('unrecognised icepa in cepac')
      endif
c
      if (lprcep.le.0) return
c
      write(iwr,1) ecorr,c0,h11ci,e0_used
      if (icepa.ge.20.and.icepa.lt.30) then
         write(iwr,8)
         do i=1,3
            write(iwr,6) (ph_pair(i,j),j=1,3)
         end do
         write(iwr,'(a)') '      --- detailed pair-energies ---'
         write(iwr,'(a)') '        0-hole => particles '
         write(iwr,'(5x,3(i4,f10.5))') (i-1,epair(i,1),i=1,3)
         write(iwr,'(a)') '        1-hole => holes '
         do ip=1,3
            write(iwr,'(7x,i2,9h-particle,/(5x,5(i4,f10.5)))')
     1      ip-1,(i,epair(ip,1+i),i=1,mrdoc)
         end do
         do ip=1,3
            write(iwr,'(7x,i2,17h-particle 2 holes)') ip
            write(iwr,'((5x,5(2i4,f10.5)))') 
     1      ((ih,jh,epair(3+mrdoc*(ip-1)+ih,jh+1),jh=1,ih),ih=1,mrdoc) 
         end do
         write(iwr,'(a)') '      --- detailed SHIFTS ---'
         write(iwr,'(a)') '        0-hole => particles'
         write(iwr,'(5x,3(i4,f10.5))') (i-1,eshifc(i,1),i=1,3)
         write(iwr,'(a)') '        1-hole => holes '
         do ip=1,3
            write(iwr,'(7x,i2,9h-particle,/(5x,5(i4,f10.5)))')
     1      ip-1,(i,eshifc(ip,1+i),i=1,mrdoc)
         end do
         do ip=1,3
            write(iwr,'(7x,i2,17h-particle 2 holes)') ip
            write(iwr,'((5x,5(2i4,f10.5)))') 
     1      ((ih,jh,eshifc(3+mrdoc*(ip-1)+ih,jh+1),jh=1,ih),ih=1,mrdoc) 
         end do
c
      else
         if (lsinsh.gt.0) write(iwr,2) (eshifc(2,j),j=1,3)
         write(iwr,4)
         do i=1,3
            write(iwr,6) (eshifc(i,j),j=1,3)
         end do
         write(iwr,8)
         do i=1,3
            write(iwr,6) (epair(i,j),j=1,3)
         end do
      end if
      write(iwr,7)
c
c***********************************************************************
c
1     format(1x/1x,' correlation energy ',f14.9,'  c0 ',f12.9,/,
     1       '  ci-energy ',e22.14,' E0(used)',e22.14)
2     format(' doublet shifts ',(t17,8f12.7))
3     format(' singlet/triplet shifts (triplet under diagonal)')
4     format(' particle - hole shifts => holes ')
8     format(' particle - hole pair-energies => holes ')
5     format(' singlet/triplet pair-energies ')
6     format(5x,10f12.7)
7     format(1x)
9     format('   shifts     excitation class')
11    format(f12.7,3x,a6,' (',i1,',',i1,')')
c
      return
cci
cci---------------------------------------------------------------------
cci   CONTRACTED CEPA
c...  including new cepa variants
c...  the iclass array is used for adressing of the shifts
c...  perhaps it shoud (at some time) be used for addressing epair
c
1100  continue
c
c...  calculate correlation energy
c
      ecorr = 0.0d0
      do 1102 i=1,3
         epair(3,i) = epair(4,i) + epair(i,4)
         epair(4,i) = 0.0d0
         epair(i,4) = 0.0d0
1102  ecorr = ecorr + epair(3,i) + epair(2,i) + epair(1,i)
c
c...   first the simple models
c
      if(icepa.ge.100) then
        if (icepa.eq.100) then
c...       Plain MRDCEPA 0
           cepa_shift =ecorr
        else if (icepa.eq.101) then
c...          AQCC
           if (nelec.eq.1) call caserr('aqcc with 1 electron')
           A = (nelec-2.0d0)*(nelec-3.0d0)/(nelec*(nelec-1.0d0))
           cepa_shift = ecorr*A
        else if (icepa.eq.102) then
c...          ACPF
           B = (nelec-2.0d0)/nelec
           cepa_shift = ecorr*B
        end if
c
c...    Shift all non-reference states with (scaled) correlation energy
c
        do iexlvl = 0, 2
          do nhole = 0, iexlvl
            do npart = 0, iexlvl
              icl = iclass(nhole+1,npart+1,iexlvl+1)+1
              eshifc(icl,1) = cepa_shift
            end do
          end do
        end do
        icl = iclass(0+1,0+1,0+1)+1
        eshifc(icl,1) = 0.0d0
        if (lsinsh.le.0) then
          iexlvl = 1
          do nhole = 0, iexlvl
            do npart = 0, iexlvl
              icl = iclass(nhole+1,npart+1,iexlvl+1)+1
              eshifc(icl,1) = 0.0d0
            end do
          end do
        endif
c
      else if (icepa.eq.10.or.icepa.eq.11.or.icepa.eq.12) then
c
c...  allowed contributions for MR(D)CEPA(1)
c
        con22i = 1.d0
        con12i = 1.d0
        con11i = 1.d0
        con22a = 1.d0
        con12a = 1.d0
        con11a = 1.d0
        Ni = 2*ndoc
        Na = nelec - Ni
        if ((icepa.eq.11.or.icepa.eq.12).and.Ni.gt.0) then
c...  mrdcepa inactive space EPV terms
c...  L Fusti-Molnar, P Szalay, JPC 100 6288 1996
           con22i = (Ni-2.d0)*(Ni-3.d0)/(Ni*(Ni-1.d0))
           con12i = (Ni-2.d0)/Ni
           con11i = (Ni-1.d0)/Ni
        endif
c active epv
        if (icepa.eq.12) then
           if (Na.le.1) then
             con22a = 0.d0
             con12a = 0.d0
             con11a = 0.d0
           else
             con22a = (Na-2.d0)*(Na-3.d0)/(Na*(Na-1.d0))
             con12a = (Na-2.d0)/Na
             con11a = (Na-1.d0)/Na
           endif
        endif
c
c...    0 particles  0,1,2 holes       vacuum
c
        npart = 0
        do iexlvl = 0, 2
          nhole = 0
          icl   = iclass(nhole+1,npart+1,iexlvl+1)+1
          eshifc(icl,1) = 0.0d0
c
          nhole = 1
          icl   = iclass(nhole+1,npart+1,iexlvl+1)+1
          if (icl.gt.0)
     1    eshifc(icl,1) = con12i*(epair(1,3)+epair(2,3)+epair(3,3))
c
          nhole = 2
          icl   = iclass(nhole+1,npart+1,iexlvl+1)+1
          if (icl.gt.0)
     1    eshifc(icl,1) = con12i*(epair(1,2)+epair(2,2)+epair(3,2))
     1                   +con22i*(epair(1,3)+epair(2,3)+epair(3,3))
        end do
c
c...    1 particle  0,1,2 holes       doublet
c
        npart = 1
        do iexlvl = 1, 2
          nhole = 0
          icl   = iclass(nhole+1,npart+1,iexlvl+1)+1
          eshifc(icl,1) = con22a*epair(3,1)+con12a*epair(3,2)+epair(3,3)
c
          nhole = 1
          icl   = iclass(nhole+1,npart+1,iexlvl+1)+1
          eshifc(icl,1) = con12i*(epair(1,3)+epair(2,3)+epair(3,3))
     1                   +con12a*epair(3,1)+con11i*con11a*epair(3,2)
c
          nhole = 2
          icl   = iclass(nhole+1,npart+1,iexlvl+1)+1
          if (icl.gt.0)
     1    eshifc(icl,1) = con12i*(epair(1,2)+epair(2,2)+epair(3,2))
     1                   +con22i*(epair(1,3)+epair(2,3)+epair(3,3))
     1                   +epair(3,1)
c
        end do
c
c...    2 particle  0,1,2 holes       singlet/triplet
c
        iexlvl = 2
        npart  = 2
        nhole  = 0
        icl    = iclass(nhole+1,npart+1,iexlvl+1)+1
        eshifc(icl,1) = con22a*(epair(2,1)+epair(3,1))
     1                 +con12a*(epair(2,2)+epair(3,2))
     1                 +(epair(2,3)+epair(3,3))
c
        nhole = 1
        icl   = iclass(nhole+1,npart+1,iexlvl+1)+1
        eshifc(icl,1) = con12i*(epair(1,3)+epair(2,3)+epair(3,3))
     1                 +con12a*(epair(2,1)+epair(3,1))
     1                 +con11i*con11a*(epair(2,2)+epair(3,2))
c
        nhole = 2
        icl   = iclass(nhole+1,npart+1,iexlvl+1)+1
        eshifc(icl,1) = con12i*(epair(1,2)+epair(2,2)+epair(3,2))
     1                 +con22i*(epair(1,3)+epair(2,3)+epair(3,3))
     1                 +(epair(2,1)+epair(3,1))
c
        if (lsinsh.le.0) then
          iexlvl = 1
          do nhole = 0, iexlvl
            do npart = 0, iexlvl
              icl = iclass(nhole+1,npart+1,iexlvl+1)+1
              eshifc(icl,1) = 0.0d0
            end do
          end do
        endif
c
      else
c
        call caserr('unimplemented contracted cepa option')
c
      endif
c
      if (lprcep.le.0) return
c
      write(iwr,1) ecorr,c0,h11ci
      write(iwr,8)
      do 1421 i=1,3
1421  write(iwr,6) (epair(i,j),j=1,3)
      write(iwr,9)
      do 1410 iexlvl = 0, 2
        do 1420 nhole = 0, iexlvl
          do 1430 npart = 0, iexlvl
            icl = iclass(nhole+1,npart+1,iexlvl+1)+1
            write(iwr,11)eshifc(icl,1),txt(iexlvl+1),nhole,npart
1430      continue
1420    continue
1410  continue
c
      return
      end
      subroutine confprt(conf)
      implicit integer (a-z)
      integer conf
      common /mr_mp2con/ nsmpv,nsmp1,nsmp2,nsmp2d,ntotmp,
     *                   nmpv,nmp1,nmp2,nmp2d
      dimension conf(5,*)
c
c...  print new configurations
c
      loop = nmpv+nmp1+nmp2+nmp2d
      do i=1,loop
       call upack2(conf(1,i), it,jj)
       write(6,*) 'confprt: i,it,jj = ', i, it, jj
      enddo
c
      return
      end
      subroutine ver_dirctc(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/dirctc.m,v $
     +     "/
      data revision /"$Revision: 6178 $"/
      data date /"$Date: 2010-08-10 16:57:58 +0200 (Tue, 10 Aug 2010) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
