c
c  Auxilliary basis set processing
c
c  This version written in terms of GAMESS-UK input processing,
c  and assumes that the DFT code will pick up options.dft on node 0
c  at the start of the DFT run.
c
c  DGauss4.0 basis files can also be read, but at the moment library
c  and input basis sets cant be mixed.
c
c  for ease of allocation in the calling routine, data is stored
c  linearly in a integer and real*8 arrays (currently a static
c  array is used
c
      subroutine store_jbasis
      implicit none
c
c use in basis set processor
c   
c   maximum number of atom types
c
      integer maxtyp
      parameter(maxtyp=10)
c
c number of specific library requests
c
      integer mxbasreq
      parameter(mxbasreq=100)
c
c symbolic names for shell types
c  (spd converted internally to sp+d)
c
      integer S,P,D,F,G,SP,SPD
      parameter(S=1,P=2,D=3,F=4,G=5,SP=6,SPD=7)
      integer nshl(maxtyp)
      character*8 ctyp(maxtyp)
      integer istart(maxtyp)
      integer iend(maxtyp)
      integer idim, fdim
c
      integer ibuff
      real*8  rbuff
      common/junkx/rbuff(10000),ibuff(10000)
c
      character*4 yform
      character*1 cform
      integer ntyp
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

      integer reqz(mxbasreq)
      character*8  reqt(mxbasreq), reqn(mxbasreq)
      character*8 name, zlib, zlab
      integer ireq, nreq, iz, iat
      logical hit, more, oignore, oprint

      external isubst
      integer isubst

      idim = 10000
      fdim = 10000

      call inpa4(yform)
      if(yform .eq. 'igno')then
         oignore = .true.
      else if(yform .eq. ' ')then
         oignore = .false.
      else
         oignore = .false.
         jrec = jrec - 1
      endif

      oprint =.false.
      call inpa4(yform)
      if(yform .eq. 'nwch')then
c
c...     load contracted basis from input stream
c
         call read_basis(ntyp,ctyp,nshl,istart,
     &        ibuff,rbuff,idim,fdim,'n')
      else if(yform .eq. 'game'.or.yform .eq. '    ')then
c
c...     load contracted basis from input stream
c
         call read_basis(ntyp,ctyp,nshl,istart,
     &        ibuff,rbuff,idim,fdim,'g')
      else 
         jrec = jrec-1
c
c...     load contracted basis from internal tables
c
         call read_basis_intern(ntyp,ctyp,nshl,istart,
     &        ibuff,rbuff,idim,fdim,oprint)
      endif
c
c...  echo fit basis set to standard output
c
      call printbasis(ntyp,ctyp,nshl,istart,ibuff,rbuff,oprint)
c
c...  normalise basis functions
c
      call norm_basis(ntyp,ctyp,nshl,istart,ibuff,rbuff)
c
c...  write out in options.dft format
c
      call write_dftbasis(ntyp,ctyp,nshl,istart,ibuff,rbuff)
c
c...  check that all atoms have a basis set
c
      call chkfitbas(oignore,ctyp,ntyp,zaname,nat)
      end
c
c-----------------------------------------------------------------------
c
      subroutine norm_basis(ntyp,ctyp,nshl,istart,index,value)
      implicit none
c
c...  Takes the (contracted) fitting basis as input and returns with
c...  the contraction coefficients modified in such a way as to make
c...  each contracted function normalised to 1.
c
c...  Input variables
c
c...     ntyp      number of atom types
c...     ctyp      character tag for atom types
c...     nshl(i)   number of shells on atom type i
c...     istart(i) start indexing point for shell of atom type i
c...     index     integer indexing array
c
c...  Input/Output variables
c
c...     values    floating point data
c
c...  Data on index array:
c
c...     The first maxtyp elements are a copy of istart. The numbers
c...     passed that contain at offsets, relative to istart(i)
c
c...     0   number of contracted function
c...     1   shell type
c...     2   address of floating values
c...     3   address of next shell on this atom
c
c...  -------------
c
c
c use in basis set processor
c   
c   maximum number of atom types
c
      integer maxtyp
      parameter(maxtyp=10)
c
c number of specific library requests
c
      integer mxbasreq
      parameter(mxbasreq=100)
c
c symbolic names for shell types
c  (spd converted internally to sp+d)
c
      integer S,P,D,F,G,SP,SPD
      parameter(S=1,P=2,D=3,F=4,G=5,SP=6,SPD=7)
      integer     ntyp
      character*8 ctyp(maxtyp)
      integer     nshl(*)
      integer     istart(maxtyp)
      integer     index(*)
      real*8        value(*)
c
c...  Local variables
c
      integer ityp, ishl, iprim, jprim
      integer ik, ik1, ik2
      integer shell_type
      integer ix 
      real*8    facs, facp, facd, facf, facg, ee, ee1, c11, c12
      real*8    fac1, fac2, dum1, dum2, fac
c
      real*8 pt5, pt75, pt187, pt6562
      data pt5,pt75,pt187,pt6562/0.5d0,0.75d0,1.875d0,6.5625d0/
c
      real*8 pi32
      data pi32/5.56832799683170d0/
c
      ix = maxtyp+1-4
      ik = 0
      do ityp = 1, ntyp
         do ishl = 1, nshl(ityp)
            ix = ix + 4
            shell_type = index(ix+1)
c
c...        Account for the normalisation constant of each primitive
c
            ik1 = 0
            do iprim = 1, index(ix+0)
               ik1  = ik1 + 1
               ee   = 2.0d0*value(ik+ik1)
               facs = pi32/(ee*sqrt(ee))
               facp = pt5*facs/ee
               facd = pt75*facs/(ee*ee)
               facf = pt187*facs/(ee**3)
               facg = pt6562*facs/(ee**4)
               ik1  = ik1 + 1
               if (shell_type.eq.S) then
                  value(ik+ik1) = value(ik+ik1)/sqrt(facs)
               else if (shell_type.eq.P) then
                  value(ik+ik1) = value(ik+ik1)/sqrt(facp)
               else if (shell_type.eq.D) then
                  value(ik+ik1) = value(ik+ik1)/sqrt(facd)
               else if (shell_type.eq.F) then
                  value(ik+ik1) = value(ik+ik1)/sqrt(facf)
               else if (shell_type.eq.G) then
                  value(ik+ik1) = value(ik+ik1)/sqrt(facg)
               else if (shell_type.eq.SP) then
                  value(ik+ik1) = value(ik+ik1)/sqrt(facs)
                  ik1 = ik1 + 1
                  value(ik+ik1) = value(ik+ik1)/sqrt(facp)
               endif
            enddo
c
c...        Compute the normalisation constant for contracted function
c
            ik1  = 0
            fac1 = 0
            fac2 = 0
            do iprim = 1, index(ix+0)
               ik1 = ik1 + 1
               ee1 = value(ik+ik1)
               ik1 = ik1 + 1
               c11 = value(ik+ik1)
               if (shell_type.eq.SP) then
                  ik1 = ik1 + 1
                  c12 = value(ik+ik1)
               endif
               ik2 = 0
               do jprim = 1, iprim
                  ik2  = ik2 + 1
                  ee   = ee1 + value(ik+ik2)
                  fac  = ee*sqrt(ee)
                  ik2  = ik2 + 1
                  dum2 = 0.0d0
                  if (shell_type.eq.S) then
                     dum1 = c11*value(ik+ik2)/fac
                  else if (shell_type.eq.P) then
                     dum1 = pt5*c11*value(ik+ik2)/(ee*fac)
                  else if (shell_type.eq.D) then
                     dum1 = pt75*c11*value(ik+ik2)/(ee*ee*fac)
                  else if (shell_type.eq.F) then
                     dum1 = pt187*c11*value(ik+ik2)/(ee**3*fac)
                  else if (shell_type.eq.G) then
                     dum1 = pt6562*c11*value(ik+ik2)/(ee**4*fac)
                  else if (shell_type.eq.SP) then
                     dum1 = c11*value(ik+ik2)/fac
                     ik2  = ik2 + 1
                     dum2 = pt5*c12*value(ik+ik2)/(ee*fac)
                  endif
                  if (jprim.ne.iprim) then
                     dum1 = dum1 + dum1
                     dum2 = dum2 + dum2
                  endif
                  fac1 = fac1 + dum1
                  fac2 = fac2 + dum2
               enddo
            enddo
c
c...        Now adjust the contraction coefficients 
c
            ik1 = 0
            do iprim = 1, index(ix+0)
               ik1 = ik1 + 2
               value(ik+ik1) = value(ik+ik1)/sqrt(fac1*pi32)
               if (shell_type .eq. SP) then
                  ik1 = ik1 + 1
                  value(ik+ik1) = value(ik+ik1)/sqrt(fac2*pi32)
               endif
            enddo
c
            ik = ik + ik1
         enddo
      enddo
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine chkfitbas(oignore,reqt,nreq,zaname,nat)
      implicit none
c
c     This routine checks whether each element in zaname has a fitting
c     basis set. If oignore is false a missing basis set will result
c     in the calculation being terminated.
c
      logical oignore
      integer nreq, nat
      character*8 reqt(nreq), zaname(nat)
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
c     local variables
c
      integer iat, ireq, iz
      logical ofound, oallfound
c
c     functions
c
      integer isubst
      external isubst
c
      oallfound = .true.
      do iat = 1, nat
         iz = isubst(zaname(iat))
         if (iz.gt.0) then
c
c           zaname(iat) is the name of a true chemical element
c
            ireq = 0
            ofound = .false.
 10         if (ireq.lt.nreq.and..not.ofound) then
               ireq = ireq + 1
               ofound = ofound.or.(zaname(iat).eq.reqt(ireq))
               goto 10
            endif
            oallfound = oallfound.and.ofound
            if (.not.ofound) then
               if (oignore) then
                  write(iwr,60)"WARNING:",iat,zaname(iat)
               else
                  write(iwr,60)"ERROR:  ",iat,zaname(iat)
               endif
            endif
         endif
      enddo
      if (.not.oignore.and..not.oallfound) then
         write(iwr,70)
         call caserr('No fitting basis defined for some atoms')
      endif
 60   format(a8,'no fitting basis defined for atom',i6,' called ',a8)
 70   format('No fitting basis set was defined for some atoms',/,
     +       'Use the "ignore" subdirective if you want to run',/,
     +       'this calculation anyway.')
      end
c
c
c input variables
c
c   idim       available length of integer array
c   fdim       available length of float array
c
c  output variables
c
c   ntyp      number of atom types
c   ctyp      character tag for atom types
c   nshl(i)   number of shells on type i
c   istart(i) start indexing point for shell i
c   index     integer indexing array
c   values    floating point data
c
c  offsets into the index array, relative to istart(i)
c
c      0   number of contracted function
c      1   shell type 
c      2   address of floating values
c      3   address of next shell on this atom
c
      subroutine read_basis(
     &     ntyp,ctyp,nshl,istart,index,value,
     &     idim, fdim, format)

      implicit none

c
c use in basis set processor
c   
c   maximum number of atom types
c
      integer maxtyp
      parameter(maxtyp=10)
c
c number of specific library requests
c
      integer mxbasreq
      parameter(mxbasreq=100)
c
c symbolic names for shell types
c  (spd converted internally to sp+d)
c
      integer S,P,D,F,G,SP,SPD
      parameter(S=1,P=2,D=3,F=4,G=5,SP=6,SPD=7)

      integer ntyp
      character*8 ctyp(maxtyp)
      integer nshl(*)
      integer istart(maxtyp)
      integer index(*)
      real*8  value(*)
      character *4 ytmp, y1, y2
      integer ix, iftop, itop
      real*8 c1, c2, c3, ex
      logical onew
      integer i
      integer atom_type, shell_type
      character *1 format
      integer idim, fdim

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
c pointers
c
      itop = maxtyp + 1
      iftop = 1
c
c initialisation
c
      ntyp=0

 100  continue
c
c load a record from input stream
c
      call input
c
c try and read two strings, to check if it is another component 
c of the current shell, or the end of the basis

      call inpa4(y1)
      call inpa4(y2)

      if(format .eq. 'n')then
c
c nwchem format is c s
c
         ytmp  = y2
      else
c
c gamess is of the form s c
c
         ytmp  = y1
      endif
c
      onew = .true.
      if(y1.eq.'end ')then
c     
c all done, return
c
         return
      else if(ytmp .eq. 's   ')then
         shell_type = S
      else if(ytmp .eq. 'p   ')then
         shell_type = P
      else if(ytmp .eq. 'd   ')then
         shell_type = D
      else if(ytmp .eq. 'f   ')then
         shell_type = F
      else if(ytmp .eq. 'g   ')then
         shell_type = G
      else if(ytmp .eq. 'sp   ' .or. ytmp.eq.'l')then
         shell_type = SP
      else
         onew = .false.
      endif
c
c  new shell
c
      if(onew)then
c
c it is a new shell, locate the correct atom
c
         if(format .eq. 'n')then
            ytmp = y1
         else
            ytmp = y2
         endif

         ix=0
         do i =1,ntyp
            if(ytmp .eq. ctyp(i))then
c
c its an existing type, set address based on
c leading scratch area. This will be the 
c address of the existing shell data segment
c so set chain pointer to the new shell.
c
               atom_type = i
               ix = index(atom_type)
               index(ix + 3)  = itop
               ix = index(ix + 3)
            endif
         enddo
         if(ix .eq. 0)then
c
c its a new atom type
c
            ntyp = ntyp + 1
            atom_type = ntyp
            ctyp(atom_type) = ytmp
            nshl(atom_type) = 0
            istart(atom_type) = itop
            ix = istart(atom_type)
         endif
c
c  set integer information for this type
c
         if(ix+3.gt.idim)call caserr('basis i buffer too small')

         index(ix + 0) = 0
         index(ix + 1) = shell_type
         index(ix + 2) = iftop
         index(ix + 3) = 0
c work pointer for current shell
         index(atom_type) = ix
         nshl(atom_type) = nshl(atom_type) + 1
c
c location for next shell
c
         itop = itop + 4
c
c load the next record containing the first set of exponents and coefficients
c
         call input
      else
c
c it's another component of an existing shell
c backspace read routine so we get to read 
c the data again
c blank or unrecognised lines will end up here also
c
         jrec = jrec - 2

      endif
c
c quietly skip blank lines
c
      if ( jump .eq. 0 ) goto 100
c
c read the contraction info
c
      if(format .eq. 'n')then
         call inpf(ex)
         call inpf(c1)
         if(shell_type .eq. SP)then
            if(jrec .ge. jump)call caserr('missing p coef on sp basis')
            call inpf(c2)
         else if(shell_type .eq. SPD)then
            if(jrec .ge. jump)call caserr('missing p coef on sp basis')
            call inpf(c2)
            if(jrec .ge. jump)call caserr('missing d coef on spd basis')
            call inpf(c3)
         endif
      else
         call inpf(c1)
         call inpf(ex)
         if(shell_type .eq. SP)then
            if(jrec .ge. jump)call caserr('missing p coef on sp basis')
            call inpf(c2)
         else if(shell_type .eq. SPD)then
            if(jrec .ge. jump)call caserr('missing p coef on sp basis')
            call inpf(c2)
            if(jrec .ge. jump)call caserr('missing d coef on spd basis')
            call inpf(c3)
         endif
      endif
c
c update the indexing info
c
      index(ix + 0) = index(ix + 0) + 1

      if(iftop+1.gt.fdim)call caserr('basis f buffer too small')

      value(iftop) = ex
      iftop = iftop + 1

      value(iftop) = c1
      iftop = iftop + 1

      if(shell_type .eq. SP)then
         if(iftop.gt.fdim)call caserr('basis f buffer too small')
         value(iftop) = c2
         iftop = iftop + 1
      else if(shell_type .eq. SPD)then
         if(iftop+1.gt.fdim)call caserr('basis f buffer too small')
         value(iftop) = c2
         iftop = iftop + 1
         value(iftop) = c3
         iftop = iftop + 1
      endif

      goto 100

      end
c
c basis set printing routine
c
      subroutine printbasis(
     &     ntyp,ctyp,nshl,istart,index,value,oprint)
      implicit none
      integer ntyp
      logical oprint
      integer istart(*)
c
c use in basis set processor
c   
c   maximum number of atom types
c
      integer maxtyp
      parameter(maxtyp=10)
c
c number of specific library requests
c
      integer mxbasreq
      parameter(mxbasreq=100)
c
c symbolic names for shell types
c  (spd converted internally to sp+d)
c
      integer S,P,D,F,G,SP,SPD
      parameter(S=1,P=2,D=3,F=4,G=5,SP=6,SPD=7)
      integer nshl(*)
      character*8 ctyp(*)
      integer index(*)
      real*8 value(*)
      integer ix, i,j,inext,ncontr
      integer ic, ip, shell_type

      ic = 0
      do i=1,ntyp
         ix = istart(i)
         if(oprint)  write(6,*)'centre type',ctyp(i),nshl(i),ix
         do while (ix .ne. 0)
            ncontr     = index(ix+0)
            shell_type = index(ix+1)
            ip         = index(ix+2)
            inext      = index(ix+3)
            if(oprint)write(6,*)'shell type',shell_type
            do j=1,ncontr
               if(shell_type .eq. SP)then
                  if(oprint)  write(6,*)value(ip),value(ip+1),
     +                        value(ip+2)
                  ip = ip + 3
               else
                  if(oprint)write(6,*)value(ip),value(ip+1)
                  ip = ip + 2
               endif
            enddo
            ix = inext
         enddo
      enddo
      end

      subroutine write_dftbasis(
     &     ntyp,ctyp,nshl,istart,index,value)
      implicit none

c
c use in basis set processor
c   
c   maximum number of atom types
c
      integer maxtyp
      parameter(maxtyp=10)
c
c number of specific library requests
c
      integer mxbasreq
      parameter(mxbasreq=100)
c
c symbolic names for shell types
c  (spd converted internally to sp+d)
c
      integer S,P,D,F,G,SP,SPD
      parameter(S=1,P=2,D=3,F=4,G=5,SP=6,SPD=7)

      integer ntyp
      integer istart(*)
      integer nshl(*)
      character*8 ctyp(*)
      integer index(*)
      real*8 value(*)

      integer  isubst
      external isubst

      integer kkk, z
      integer ix, i,j,inext,ncontr
      integer ic, ip, shell_type
      logical opg_root
      integer ip0, nsp, nspd
      character * 2 iel
c
      dimension iel(105)
c
      data iel/'x ', 'bq',
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


      if(.not.opg_root())return


      open(50,file='options.dft',form='formatted',
     &     status='unknown')

      write(50,1000)'jbasis'
1000  format(a4)

      write(50,1001)ntyp
1001  format(i3)

      ic = 0
      do i=1,ntyp

         z = isubst(ctyp(i))
c
c at the moment we can't process basis on bqs etc
c
         if(z.le.0)then
            write(6,*)ctyp(i)
            call caserr('unrec atom type')
         endif
c
c Dummy loop to determine number of sp or spd shells to correct
c expanded shell count
c
         nsp = 0
         nspd = 0
         ix = istart(i)
         do while (ix .ne. 0)
            shell_type = index(ix+1)
            inext      = index(ix+3)
            if(shell_type .eq. SP)nsp = nsp+1
            if(shell_type .eq. SPD)nspd = nspd+1
            ix = inext
         enddo

         write(50,*)z
         write(50,*)nshl(i) + nsp + 2*nspd

         ix = istart(i)
         do while (ix .ne. 0)
            ncontr     = index(ix+0)
            shell_type = index(ix+1)
            ip         = index(ix+2)
            inext      = index(ix+3)

            if(shell_type .eq. SP)then
               ip0=ip
               write(50,*)ncontr
               write(50,1001)1
               write(50,1001)1
               do j=1,ncontr
                  write(50,1002)value(ip),value(ip+1)
                  ip = ip + 3
               enddo
               ip = ip0
               write(50,*)ncontr
               write(50,1001)2
               write(50,1001)2
               do j=1,ncontr
                  write(50,1002)value(ip),value(ip+2)
                  ip = ip + 3
               enddo
            elseif(shell_type .eq. SPD)then
               ip0=ip
               write(50,*)ncontr
               write(50,1001)1
               write(50,1001)1
               do j=1,ncontr
                  write(50,1002)value(ip),value(ip+1)
                  ip = ip + 4
               enddo
               ip = ip0
               write(50,*)ncontr
               write(50,1001)2
               write(50,1001)2
               do j=1,ncontr
                  write(50,1002)value(ip),value(ip+2)
                  ip = ip + 4
               enddo
               ip = ip0
               write(50,*)ncontr
               write(50,1001)3
               write(50,1001)3
               do j=1,ncontr
                  write(50,1002)value(ip),value(ip+3)
                  ip = ip + 4
               enddo
            else
               write(50,*)ncontr
               write(50,1001)shell_type
               write(50,1001)shell_type
               do j=1,ncontr
                  write(50,1002)value(ip),value(ip+1)
                  ip = ip + 2
               enddo
            endif

c            write(50,*)ncontr
c            if(shell_type .eq. SP)then
c               write(50,1001)2
c               write(50,1001)1
c            else
c               write(50,1001)shell_type
c               write(50,1001)shell_type
c            endif
c
c            if(shell_type .eq. SP)then
c               do j=1,ncontr
c                  write(50,1002)value(ip),value(ip+1),value(ip+2)
c                  ip = ip + 3
c               enddo
c            else
c               do j=1,ncontr
c                  write(50,1002)value(ip),value(ip+1)
c                  ip = ip + 2
c               enddo
c            endif

c 1002       format(3e20.14)
1002       format(3f25.8)
            ix = inext
         enddo
      enddo

      close(50)

      end

      subroutine read_basis_intern(
     &     ntyp,ctyp,nshl,istart,index,value,
     &     idim, fdim,oprint)

      implicit none
c
c input variables
c
c   idim       available length of integer array
c   fdim       available length of float array
c
c  output variables
c
c   ntyp      number of atom types
c   ctyp      character tag for atom types
c   nshl(i)   number of shells on type i
c   istart(i) start indexing point for shell i
c   index     integer indexing array
c   values    floating point data
c
c  offsets into the index array, relative to istart(i)
c
c      0   number of contracted function
c      1   shell type 
c      2   address of floating values
c      3   address of next shell on this atom
c
c
c
c use in basis set processor
c   
c   maximum number of atom types
c
      integer maxtyp
      parameter(maxtyp=10)
c
c number of specific library requests
c
      integer mxbasreq
      parameter(mxbasreq=100)
c
c symbolic names for shell types
c  (spd converted internally to sp+d)
c
      integer S,P,D,F,G,SP,SPD
      parameter(S=1,P=2,D=3,F=4,G=5,SP=6,SPD=7)
c
      integer ntyp
      logical oprint
      character*8 ctyp(maxtyp)
      integer nshl(*)
      integer istart(maxtyp)
      integer index(*)
      real*8  value(*)
      integer ix, iftop, itop
      real*8 c1, c2, ex
      logical onew
      logical opg_root
      integer i
      integer atom_type, shell_type
      integer idim, fdim
      character*10 ztagbs,libname
      character*8 zatagi
      integer libnames
      parameter (libnames=4)
      dimension ztagbs(libnames)
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
c
      real*8 stos1, stos2, stos3, stos4, stos5, stos6, stos7
      real*8 rleesf
      real*8 c11al, c11be, c31al, c31be, c33al, c33be
      integer nbfs, minf, maxf, nangm
c
      common /datbas/ stos1(54),stos2(54),stos3(54),stos4(54),
     + stos5(54),stos6(54),stos7(54),
     + rleesf(36,2),
     + nbfs(24),minf(24),maxf(24),nangm(24),
     + c11al(12),c11be(12),c31al(12),c31be(12),c33al(12),c33be(12)
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
      integer ngauss, ishell, k, k1, loop, nati, natj
      integer nucz, jsubst, nshell, ityp
      integer jumpj,j,igauss
      integer locatc
c
      integer ibuff, intyp, kng, kstart
      real*8  rbuff, eex, cs, cp, cd, cf, cg
      common/junkx/rbuff(10000),ibuff(10000),
     +  eex(1000), cs(1000), cp(1000), cd(1000), 
     +  cf(1000), cg(1000), 
     +  intyp(1000),kng(1000),kstart(1000)
c
      data ztagbs / 'a1-dgauss', 'a2-dgauss', 'demon', 'ahlrichs' /
c 
c pointers
c
      itop = maxtyp + 1
      iftop = 1
c
c initialisation
c
      ntyp=0
c
c determine internal basis set name .. default to dgauss_a1
c
      call inpan(libname)
      jumpj = locatc(ztagbs,libnames,libname)
      if (jumpj.eq.0) then
       jumpj = 1
       libname = ztagbs(1)
      endif
      if(opg_root()) write(iwr,9569)libname
9569  format(1x,
     *'* fitting basis selected is ',a10,' *'/)
c
c ... drive through the nat centres
c ... eliminating those without basis functions
c
      do 3005 nati=1,nonsym
      zatagi=ztag(nati)
      j=nati-1
      if(j.eq.0)go to 3040
      natj=locatc(ztag,j,zatagi)
      if(natj)3005,3040,3005
3040  nucz = jsubst(zatagi)
      if(nucz.le.0)go to 3005
c
c     pick up appropriate internal library routine
c
      go to (13000,13100,13200,13300), jumpj
c
c     DGauss A1
13000 call dgauss_a1_fitl(eex,cs,cp,cd,
     +  intyp,kng,kstart,
     +  nangm,nbfs,minf,maxf,ngauss,nshell,nucz)
      go to 3010
c     DGauss A2
13100 call dgauss_a2_fitl(eex,cs,cp,cd,
     +  intyp,kng,kstart,
     +  nangm,nbfs,minf,maxf,ngauss,nshell,nucz,iwr)
      go to 3010
c     Demon
13200 call demon_fitl(eex,cs,cp,cd,
     +  intyp,kng,kstart,
     +  nangm,nbfs,minf,maxf,ngauss,nshell,nucz)
      go to 3010
c     Ahlrichs
13300 call ahlrichs_fitl(eex,cs,cp,cd,cf,cg,
     +  intyp,kng,kstart,
     +  nangm,nbfs,minf,maxf,ngauss,nshell,nucz)
c
3010  continue
c
c its a new atom type
c
      ntyp = ntyp + 1
      atom_type = ntyp
      ctyp(atom_type) = zatagi
      nshl(atom_type) = 0
      istart(atom_type) = itop
      ix = istart(atom_type)
c
      if(oprint) write(6,*)'nucz,ntyp,nshell = ', nucz,ntyp,nshell
c
      do ishell = 1 , nshell
c
      ityp = intyp(ishell)
      if (ityp.eq.16) then
       shell_type = S
      else if (ityp.eq.17) then
       shell_type = P
      else if (ityp.eq.18) then
       shell_type = D
      else if (ityp.eq.19) then
       shell_type = F
      else if (ityp.eq.20) then
       shell_type = G
      else if (ityp.eq.22) then
       shell_type = SP
      else
       call caserr('invalid shell type')
      endif
      if(ishell.gt.1) then
       ix = index(atom_type)
       index(ix + 3)  = itop
       ix = index(ix + 3)
      endif
c
c  set integer information for this type
c
         if(ix+3.gt.idim)call caserr('basis i buffer too small')

         index(ix + 0) = 0
         index(ix + 1) = shell_type
         index(ix + 2) = iftop
         index(ix + 3) = 0
c work pointer for current shell
         index(atom_type) = ix
         nshl(atom_type) = nshl(atom_type) + 1
c
c location for next shell
c
         itop = itop + 4
c
         igauss = kng(ishell)
         k1 = kstart(ishell)
         if(oprint) write(6,*)'shell_type,igauss,ityp,k1, = ', 
     +              shell_type,igauss,ityp,k1
         do loop = 1, igauss
         k = k1+loop-1
         ex = eex(k)
         if(shell_type.eq.S) then
          c1 = cs(k)
         else if(shell_type.eq.P) then
          c1 = cp(k)
         else if(shell_type.eq.D) then
          c1 = cd(k)
         else if(shell_type.eq.F) then
          c1 = cf(k)
         else if(shell_type.eq.G) then
          c1 = cg(k)
         else if(shell_type.eq.SP) then
          c1 = cs(k)
          c2 = cp(k)
         else
          call caserr('invalid shell type')
         endif
c
c update the indexing info
c
         index(ix + 0) = index(ix + 0) + 1

         if(iftop+1.gt.fdim)call caserr('basis f buffer too small')

         value(iftop) = ex
         iftop = iftop + 1

         value(iftop) = c1
         iftop = iftop + 1
         if(oprint) write(6,*)' loop, ex,c1 ', loop, ex, c1
         if(shell_type .eq. SP) then
          value(iftop) = c2
          iftop = iftop + 1
          if(oprint) write(6,*)' loop, ex,c2 ', loop, ex, c2
         endif

         enddo

      enddo

 3005  continue

      end

      subroutine dgauss_a1_fitl(ex,cs,cp,cd,
     + intyp,kng,kstart,
     + nangm,nbfs,minf,maxf,
     + ngauss,nshell,nucz)
c
      implicit none
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
      real*8 ex,cs,cp,cd
      integer intyp, nangm, kstart, kng
      integer nbfs, minf, maxf
      integer ngauss,nshell,nucz
c
      dimension ex(*),cs(*),cp(*),cd(*)
      dimension intyp(*),nangm(*)
      dimension kstart(*),kng(*)
      dimension nbfs(*),minf(*),maxf(*)
c
      real*8 eex, ccs, ccp, ccd
      common/blkin/eex(50),ccs(50),ccp(50),ccd(50)
c
      integer ng, igauss, ityp, ipass
      integer i, k, k1
      logical odone
c
c     DGAUSS_A1_DFT Coulomb Fitting basis sets
c
      ng = 0
      igauss = 0
      ityp =  0
      nshell = 0
      ngauss = 0
      odone = .false.
      call vclr(eex,1,200)
c
c     ----- hydrogen to helium -----
c
      if (nucz .le. 2) then
          call dgauss_a1_fit0(eex,ccs,ccp,nucz)
c
c     ----- lithium to neon -----
c
      else if (nucz .le. 10) then
          call dgauss_a1_fit1(eex,ccs,ccp,ccd,nucz)
c
c     ----- sodium to argon -----
c
      else if (nucz .le. 18) then
          call dgauss_a1_fit2(eex,ccs,ccp,ccd,nucz)
c
c     ----- potassium to zinc -----
c
      else if(nucz.le.30) then
          call dgauss_a1_fit3(eex,ccs,ccp,ccd,nucz)
c
c     ----- gallium to krypton -----
c
      else if(nucz.le.36) then
          call dgauss_a1_fit4(eex,ccs,ccp,ccd,nucz)
c
c     ----- rubidium to cadmium
c
      else if (nucz .le. 48) then
          call dgauss_a1_fit5(eex,ccs,ccp,ccd,nucz)
c
c     ----- indium to xenon
c
c      else if (nucz .le. 54) then
          call dgauss_a1_fit6(eex,ccs,ccp,ccd,nucz)
c
c
c     ----- past xenon does not exist
c
      else
        call caserr(
     +   'attempting to site dgauss_a1 function on invalid centre')
      endif
c
c     ----- loop over each shell -----
c
      ipass = 0
  210 ipass = ipass+1
        call dgauss_a1_fitsh(nucz,ipass,ityp,igauss,ng,odone)
      if(odone) go to 220
c
c     ----- define the current shell -----
c
      nshell = nshell+1
      kstart(nshell) = ngauss+1
      intyp(nshell) = ityp
      kng(nshell) = igauss
      ngauss = ngauss+igauss
      k1 = kstart(nshell)
      do i = 1,igauss
         k = k1+i-1
         ex(k) = eex(ng+i)
         if(ityp.eq.16) then
          cs(k) = ccs(ng+i)
         else if(ityp.eq.17) then
          cp(k) = ccp(ng+i)
         else if(ityp.eq.18) then
          cd(k) = ccd(ng+i)
         else if(ityp.eq.22) then
          cs(k) = ccs(ng+i)
          cp(k) = ccp(ng+i)
         else
          call caserr('invalid shell type')
         endif
      enddo
c
      go to 210
c
  220 continue
      return
      end
      subroutine dgauss_a1_fitsh(nucz,ipass,itype,
     +                  igauss,ng,odone)
      implicit none
c
      integer nucz,ipass,itype,igauss,ng
      logical odone
c
      integer igau, itypes, kind
      integer ind, mxpass, ityp, loop
c
      dimension igau(16,8), itypes(16,8)
      dimension kind(8)
      data kind/4, 8,10,12,13,14,15,16/
c     note all shells are single primitives in a1-fitting basis
      data igau /
     + 1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,
     + 1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,
     + 1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,
     + 1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,
     + 1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,
     + 1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
     + 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,
     + 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/
      data itypes /
     + 1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,
     + 1,1,1,1,1,7,7,3,0,0,0,0,0,0,0,0,
     + 1,1,1,1,7,7,7,3,3,3,0,0,0,0,0,0,
     + 1,1,1,1,1,7,7,7,7,3,3,3,0,0,0,0,
     + 1,1,1,1,1,7,7,7,7,3,3,3,3,0,0,0,
     + 1,1,1,1,1,7,7,7,7,7,3,3,3,3,0,0,
     + 1,1,1,1,1,7,7,7,7,7,3,3,3,3,3,0,
     + 1,1,1,1,1,7,7,7,7,7,7,3,3,3,3,3 /
c
c     set values for the current dgauss_a1 shell
c
c     return igauss = number of gaussians in current shell
c            ityp   = 1,2,7,3 for s,p,sp,d shell
c            ng     = offset in e,cs,cp,cd arrays for current shell
c
c
      ind=1
      if(nucz.gt. 2) ind=2
      if(nucz.gt. 4) ind=3
      if(nucz.gt.10) ind=4
      if(nucz.gt.12) ind=5
      if(nucz.gt.18) ind=6
      if(nucz.gt.30) ind=7
      if(nucz.gt.36) ind=8
      if(nucz.gt.38) ind=7
      if(nucz.gt.54) then
         call caserr('dgauss_a1 basis sets only extend to xeon')
      end if
c
      mxpass=kind(ind)
c
      if(ipass.gt.mxpass) odone=.true.
      if(.not.odone) then
c
       igauss = igau(ipass,ind)
       ng =0
       do loop = 1, ipass-1
        ng = ng + igau(loop,ind)
       enddo
       ityp = itypes(ipass,ind)
c
      endif
c
      itype = ityp + 15
      return
      end
      subroutine dgauss_a2_fitl(ex,cs,cp,cd,
     + intyp,kng,kstart,nangm,nbfs,minf,maxf,
     + ngauss,nshell,nucz,iwr)
c
      implicit none
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
      real*8 ex,cs,cp,cd
      integer intyp, nangm, kstart, kng
      integer nbfs, minf, maxf
      integer ngauss,nshell,nucz,iwr
c
      dimension ex(*),cs(*),cp(*),cd(*)
      dimension intyp(*),nangm(*)
      dimension kstart(*),kng(*)
      dimension nbfs(*),minf(*),maxf(*)

      real*8 eex, ccs, ccp, ccd
      common/blkin/eex(50),ccs(50),ccp(50),ccd(50)
c
      integer ng, igauss, ityp, ipass
      integer i, k, k1
      logical odone
c
c     DGAUSS_A2_DFT Coulomb Fitting basis sets
c     N. Godbout, D. R. Salahub, J. Andzelm, and E. Wimmer,   
c     Can. J. Chem. 70, 560 (1992).                            
c     Elements             
c     H  - He: ( 4s,1p,1d)  
c     Li - Be: ( 7s,2p,1d)  
c     B  - Ne: ( 8s,4p,4d)                                                           
c     Na - Mg: ( 9s,4p,3d)                                                           
c     Al - Ar: ( 9s,4p,4d)                                                           
c     Sc - Zn: (10s,5p,5d)                                                           
**                                                                             
c     DGauss basis sets provided courtesy of Cray Research, Inc.                     
c
      ng = 0
      igauss = 0
      ityp = 0
      nshell = 0
      ngauss = 0
      odone = .false.
      call vclr(eex,1,200)
c
c     ----- hydrogen to helium -----
c
      if (nucz .le. 2) then
          call dgauss_a2_fit0(eex,ccs,ccp,ccd,nucz)
c
c     ----- lithium to neon -----
c
      else if (nucz .le. 10) then
          call dgauss_a2_fit1(eex,ccs,ccp,ccd,nucz,iwr)
c
c     ----- sodium to argon -----
c
      else if (nucz .le. 18) then
          call dgauss_a2_fit2(eex,ccs,ccp,ccd,nucz,iwr)
c
c     ----- potassium to zinc -----
c
      else if(nucz.le.30) then
          call dgauss_a2_fit3(eex,ccs,ccp,ccd,nucz,iwr)
c
c     ----- past zinc does not exist
c
      else
        call caserr(
     +   'attempting to site dgauss_a2 function on invalid centre')
      endif
c
c     ----- loop over each shell -----
c
      ipass = 0
  210 ipass = ipass+1
        call dgauss_a2_fitsh(nucz,ipass,ityp,igauss,ng,odone)
      if(odone) go to 220
c
c     ----- define the current shell -----
c
      nshell = nshell+1
      kstart(nshell) = ngauss+1
      intyp(nshell) = ityp
      kng(nshell) = igauss
      ngauss = ngauss+igauss
      k1 = kstart(nshell)
      do i = 1,igauss
         k = k1+i-1
         ex(k) = eex(ng+i)
         if(ityp.eq.16) then
          cs(k) = ccs(ng+i)
         else if(ityp.eq.17) then
          cp(k) = ccp(ng+i)
         else if(ityp.eq.18) then
          cd(k) = ccd(ng+i)
         else if (ityp.eq.22) then
          cs(k) = ccs(ng+i)
          cp(k) = ccp(ng+i)
         else
          call caserr('invalid shell type')
         endif
       enddo
c
       go to 210
c
  220 continue
      return
      end
      subroutine dgauss_a2_fitsh(nucz,ipass,itype,
     +                           igauss,ng,odone)
      implicit none
c
      integer nucz,ipass,itype,igauss,ng
      logical odone
c
      integer igau, itypes, kind
      integer ind, mxpass, ityp, loop
c
      dimension igau(15,5), itypes(15,5)
      dimension kind(5)
      data kind/5, 8,12,13,15/
c     note all shells are single primitives in a1-fitting basis
      data igau /
     + 1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,
     + 1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,
     + 1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,
     + 1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
     + 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/
      data itypes /
     + 1,1,1,7,3,0,0,0,0,0,0,0,0,0,0,
     + 1,1,1,1,1,7,7,3,0,0,0,0,0,0,0,
     + 1,1,1,1,7,7,7,7,3,3,3,3,0,0,0,
     + 1,1,1,1,1,7,7,7,7,3,3,3,3,0,0,
     + 1,1,1,1,1,7,7,7,7,7,3,3,3,3,3 /
c
c    pointers to the a1-fitting basis
c    set values for the current dgauss_a2 shell
c
c     return igauss = number of gaussians in current shell
c            ityp   = 1,2,3,7 for s,p,d,sp shell
c            ng     = offset in e,cs,cp,cd arrays for current shell
c
      ind=1
      if(nucz.gt. 2) ind=2
      if(nucz.gt. 4) ind=3
      if(nucz.gt.10) ind=4
      if(nucz.gt.18) ind=5
      if(nucz.gt.30) then
         call caserr('dgauss_a2 basis sets only extend to zinc')
      end if
c
      mxpass=kind(ind)
c
      if(ipass.gt.mxpass) odone=.true.
      if(.not.odone) then
c
       igauss = igau(ipass,ind)
       ng =0
       do loop = 1, ipass-1
        ng = ng + igau(loop,ind)
       enddo
       ityp = itypes(ipass,ind)
c
      endif
c
      itype = ityp + 15
      return
      end
      subroutine demon_fitl(ex,cs,cp,cd,
     + intyp,kng,kstart,
     + nangm,nbfs,minf,maxf,
     + ngauss,nshell,nucz)
c
      implicit none
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
      real*8 ex,cs,cp,cd
      integer intyp, nangm, kstart, kng
      integer nbfs, minf, maxf
      integer ngauss,nshell,nucz
c
      dimension ex(*),cs(*),cp(*),cd(*)
      dimension intyp(*),nangm(*)
      dimension kstart(*),kng(*)
      dimension nbfs(*),minf(*),maxf(*)
c
      real*8 eex, ccs, ccp, ccd
      common/blkin/eex(50),ccs(50),ccp(50),ccd(50)
c
      integer ng, igauss, ityp, ipass
      integer i, k, k1
      logical odone
c
c     DeMon DFT Coulomb Fitting basis sets
************************************************************
c      N. Godbout, D. R. Salahub, J. Andzelm, and E. Wimmer,    
c      Can. J. Chem. 70, 560 (1992).                            
************************************************************
c      H - He : (4s,1p)     
c      Li     : (7s,3p,3d)
c      Be     : (7s,2p,1d)
c      B  - Ne: (7s,3p,3d)
c      Na     : (9s,3p,3d)
c      Mg     : (9s,4p,3d)
c      Al - Ar: (9s,4p,4d)
c      K- Ca  : (10s,5p,4d)
c      Sc - Kr: (10s,5p,5d)
c      Rb - Sr: (11s,6p,5d)
c      Y  - Xe: (10s,5p,5d)
************************************************************
c
      ng = 0
      igauss = 0
      ityp = 0
      nshell = 0
      ngauss = 0
      odone = .false.
      call vclr(eex,1,200)
c
c     ----- hydrogen to helium -----
c
      if (nucz .le. 2) then
          call demon_fit0(eex,ccs,ccp,nucz)
c
c     ----- lithium to neon -----
c
      else if (nucz .le. 10) then
          call demon_fit1(eex,ccs,ccp,ccd,nucz)
c
c     ----- sodium to argon -----
c
      else if (nucz .le. 18) then
          call demon_fit2(eex,ccs,ccp,ccd,nucz)
c
c     ----- potassium to zinc -----
c
      else if(nucz.le.30) then
c     note that K-Zn are same as DGauss A1
          call dgauss_a1_fit3(eex,ccs,ccp,ccd,nucz)
c
c     ----- gallium to krypton -----
c
      else if(nucz.le.36) then
c     note that Ga-Kr are same as DGauss A1
          call dgauss_a1_fit4(eex,ccs,ccp,ccd,nucz)
c
c     ----- rubidium to cadmium
c
      else if (nucz .le. 48) then
c     note that Rb-Cd are same as DGauss A1
          call dgauss_a1_fit5(eex,ccs,ccp,ccd,nucz)
c
c     ----- indium to xenon
c
c      else if (nucz .le. 54) then
c     note that In-Xe are same as DGauss A1
          call dgauss_a1_fit6(eex,ccs,ccp,ccd,nucz)
c
c     ----- past xenon does not exist
c
      else
        call caserr(
     +   'attempting to site demon function on invalid centre')
      endif
c
c     ----- loop over each shell -----
c
      ipass = 0
  210 ipass = ipass+1
        call demon_fitsh(nucz,ipass,ityp,igauss,ng,odone)
      if(odone) go to 220
c
c     ----- define the current shell -----
c
      nshell = nshell+1
      kstart(nshell) = ngauss+1
      intyp(nshell) = ityp
      kng(nshell) = igauss
      ngauss = ngauss+igauss
      k1 = kstart(nshell)
      do i = 1,igauss
         k = k1+i-1
         ex(k) = eex(ng+i)
         if(ityp.eq.16) then
          cs(k) = ccs(ng+i)
         else if(ityp.eq.17) then
          cp(k) = ccp(ng+i)
         else if(ityp.eq.18) then
          cd(k) = ccd(ng+i)
         else if (ityp.eq.22) then
          cs(k) = ccs(ng+i)
          cp(k) = ccp(ng+i)
         else
          call caserr('invalid shell type')
         endif
      enddo
      go to 210
c
  220 continue
      return
      end
      subroutine demon_fitsh(nucz,ipass,itype,
     +                        igauss,ng,odone)
      implicit none
c
      integer nucz,ipass,itype,igauss,ng
      logical odone
c
      integer igau, itypes, kind
      integer ind, mxpass, ityp, loop
c
      dimension igau(16,9), itypes(16,9)
      dimension kind(9)
      data kind/4, 10, 8, 12, 12, 13, 14, 15, 16/
c     note all shells are single primitives in demon-fitting basis
      data igau /
     + 1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,
     + 1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,
     + 1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,
     + 1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,
     + 1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,
     + 1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,
     + 1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
     + 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,
     + 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/
      data itypes /
     + 1,1,1,7,0,0,0,0,0,0,0,0,0,0,0,0,
     + 1,1,1,1,7,7,7,3,3,3,0,0,0,0,0,0,
     + 1,1,1,1,1,7,7,3,0,0,0,0,0,0,0,0,
     + 1,1,1,1,1,1,7,7,7,3,3,3,0,0,0,0,
     + 1,1,1,1,1,7,7,7,7,3,3,3,0,0,0,0,
     + 1,1,1,1,1,7,7,7,7,3,3,3,3,0,0,0,
     + 1,1,1,1,1,7,7,7,7,7,3,3,3,3,0,0,
     + 1,1,1,1,1,7,7,7,7,7,3,3,3,3,3,0,
     + 1,1,1,1,1,7,7,7,7,7,7,3,3,3,3,3 /
c
c    pointers to the demon-fitting basis
c
c     return igauss = number of gaussians in current shell
c            ityp   = 1,2,3,7 for s,p,d,sp shell
c            ng     = offset in e,cs,cp,cd arrays for current shell
c
      ind=1
      if(nucz.gt. 2) ind=2
      if(nucz.gt. 3) ind=3
      if(nucz.gt. 4) ind=2
      if(nucz.gt.10) ind=4
      if(nucz.gt.11) ind=5
      if(nucz.gt.12) ind=6
      if(nucz.gt.18) ind=7
      if(nucz.gt.20) ind=8
      if(nucz.gt.36) ind=9
      if(nucz.gt.38) ind=8
      if(nucz.gt.54) then
         call caserr('demon basis sets only extend to xeon')
      end if
c
      mxpass=kind(ind)
c
      if(ipass.gt.mxpass) odone=.true.
      if(.not. odone) then
c
       igauss = igau(ipass,ind)
       ng =0
       do loop = 1, ipass-1
        ng = ng + igau(loop,ind)
       enddo
       ityp = itypes(ipass,ind)
c
      endif
c
      itype = ityp + 15
      return
      end
      subroutine ahlrichs_fitl(ex,cs,cp,cd,cf,cg,
     + intyp,kng,kstart,nangm,nbfs,minf,maxf,ngauss,
     + nshell,nucz)
c
      implicit none
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
      real*8 ex,cs,cp,cd, cf, cg
      integer intyp, nangm, kstart, kng
      integer nbfs, minf, maxf
      integer ngauss,nshell,nucz
c
      dimension ex(*),cs(*),cp(*),cd(*),cf(*),cg(*)
      dimension intyp(*),nangm(*)
      dimension kstart(*),kng(*)
      dimension nbfs(*),minf(*),maxf(*)

      real*8 eex, ccs, ccp, ccd, ccf, ccg
      common/blkin/eex(50),ccs(50),ccp(50),ccd(50),ccf(50),ccg(50)
c
      real*8 pi32, pt5, pt75, pt187, pt6562
      real*8 facs, facp, facd, facf, facg, ee, tol, done, fac
      real*8 dums, dump, dumd, dumf, dumg
      integer ng, igauss, ityp, ipass
      integer i, k, k1, k2, ig, jg
      logical odone
c
c     Ahlrichs DFT Coulomb Fitting Basis Sets                               
************************************************************
c     K. Eichkorn, O. Treutler, H. Ohm, M. Haser,R. Ahlrichs, 
c     Chem. Phys. Lett. 240, 283 (1995)    
c     K. Eichkorn, F. Weigend, O. Treutler, R. Ahlrichs, 
c     Theor. Chim. Acc. 97,    119 (1997)               
************************************************************
c    Contractions                                                              
c   H:  (4s2p1d)     -> [3s2p1d]      {211/11/1}
c   **
c   He: (4s2p)       -> [2s2p]        {31/11}
c   **
c   Li: (9s2p2d1f)   -> [7s2p2d1f]    {3111111/11/11/1}
c   Be: (9s2p2d1f)   -> [7s2p2d1f]    {3111111/11/11/1}
c   B:  (9s3p3d1f)   -> [7s3p3d1f]    {3111111/111/111/1}
c   **                                                                               
c   C:  (9s3p3d1f)   -> [7s3p3d1f]    {3111111/111/111/1}
c   N:  (9s3p3d1f)   -> [7s3p3d1f]    {3111111/111/111/1}
c   O:  (9s3p3d1f)   -> [7s3p3d1f]    {3111111/111/111/1}
c   F:  (9s3p4d1f)   -> [7s3p4d1f]    {3111111/111/1111/1}
c   Ne: (8s3p3d1f)   -> [6s3p3d1f]    {311111/111/111/1}
c   **                                                                             
c   Na: (12s4p4d1f)   / [5s2p2d1f]    {81111/31/31/1}
c   Mg: (12s4p4d1f)   / [5s2p2d1f]    {81111/31/31/1}
c   Al: (12s4p5d1f)   / [5s3p2d1f]    {81111/211/41/1}
c   Si: (12s6p5d1f1g) / [5s3p2d1f1g]  {81111/411/41/1/1}
c   P:  (12s6p5d1f1g) / [5s3p2d1f1g]  {81111/411/41/1/1}
c   S:  (12s6p5d1f1g) / [5s3p2d1f1g]  {81111/411/41/1/1}
c   Cl: (12s6p5d1f1g) / [5s3p2d1f1g]  {81111/411/41/1/1} 
c   Ar: (12s6p5d1f)   / [5s3p2d1f]    {81111/411/41/1}
c   **                                                                             
c   K:  (16s4p4d1f)   / [6s2p2d1f]    {10 21111/31/31/1}
c   Ca: (16s4p4d1f)   / [6s2p2d1f]    {10 21111/31/31/1}
c   Sc: (16s4p4d3f4g) / [6s4p3d3f2g]  {10 21111/1111/211/111/31}
c   Ti: (16s4p4d3f4g) / [6s4p2d3f2g]  {10 21111/1111/31/111/31}
c   V:  (16s4p4d3f4g) / [6s4p2d3f2g]  {10 21111/1111/31/111/31}
c   Cr: (16s4p4d3f4g) / [6s4p2d3f2g]  {10 21111/1111/31/111/31}
c   Mn: (16s4p4d3f4g) / [6s4p2d3f2g]  {10 21111/1111/31/111/31}
c   Fe: (16s4p4d3f4g) / [6s4p2d3f2g]  {10 21111/1111/31/111/31}
c   Co: (16s4p4d3f4g) / [6s4p3d3f2g]  {10 21111/1111/211/111/31} 
c   Ni: (17s4p4d3f4g) / [7s4p2d3f2g]  {10 211111/1111/31/111/31} 
c   Cu: (17s4p4d3f4g) / [7s4p2d3f2g]  {10 211111/1111/31/111/31}
c   Zn: (17s4p4d3f4g) / [7s4p2d3f2g]  {10 211111/1111/31/111/31} 
c   Ga: (16s4p4d1f)   / [6s4p2d1f]    {10 21111/1111/31/1}
c   Ge: (16s4p5d1f1g) / [6s4p2d1f1g]  {10 21111/1111/41/1/1}
c   As: (16s4p5d1f1g) / [6s4p2d1f1g]  {10 21111/1111/41/1/1}
c   Se: (16s4p5d1f1g) / [6s4p2d1f1g]  {10 21111/1111/41/1/1}
c   Br: (16s4p5d1f1g) / [6s4p2d1f1g]  {10 21111/1111/41/1/1}
c   Kr: (17s4p5d1f1g) / [7s4p2d1f1g]  {10 211111/1111/41/1/1}
c   **                                                                             
c   Rb: (8s4p4d)      / [5s2p2d]      {41111/31/31}
c   Sr: (9s3p3d2f1g)  / [6s3p3d2f1g]  {411111/111/111/11/1}
c   Y:  (8s4p4d3f4g)  / [6s4p2d3f2g]  {311111/1111/31/111/31}
c   Zr: (10s4p4d3f4g) / [6s4p2d3f2g]  {511111/1111/31/111/31}
c   Nb: (8s4p4d3f4g)  / [6s4p3d3f2g]  {311111/1111/211/111/31}
c   Mo: (8s4p4d3f4g)  / [7s4p3d3f2g]  {311111/1111/211/111/31}
c   Tc: (8s4p4d3f4g)  / [6s4p3d3f2g]  {3 11111/1111/211/111/31}
c   Ru: (9s4p4d3f4g)  / [7s4p3d3f2g]  {3111111/1111/211/111/31}
c   Rh: (9s4p4d3f4g)  / [7s4p4d3f2g]  {3111111/1111/1111/111/31}
c   Pd: (9s4p4d3f4g)  / [7s4p2d3f2g]  {3111111/1111/31/111/31}
c   Ag: (9s4p4d3f4g)  / [7s4p3d3f2g]  {3111111/1111/211/111/31}
c   Cd: (9s4p4d3f4g)  / [7s4p2d3f2g]  {3111111/1111/31/111/31}
c   In: (5s3p3d1f1g)  / [3s3p2d1f1g]  {311/111/21/1/1}
c   Sn: (5s3p3d1f1g)  / [4s3p3d1f1g]  {2111/111/111/1/1}
c   Sb: (5s3p3d1f1g)  / [4s3p3d1f1g]  {2111/111/111/1/1}
c   Te: (5s3p3d1f1g)  / [4s3p3d1f1g]  {2111/111/111/1/1}
c   I:  (5s3p3d1f1g)  / [3s3p2d1f1g]  {311/111/21/1/1}
c   Xe: (8s3p3d1f1g)  / [5s3p2d1f1g]  {41111/111/21/1/1}
c   **                                                                             
c   Cs: (9s4p4d)      / [5s2p2d]      {51111/31/31}
c   Ba: (9s3p3d2f1g)  / [6s3p3d2f1g]  {411111/111/111/11/1}
c   Hf: (9s4p4d3f4g)  / [7s4p4d3f2g]  {3111111/1111/1111/111/31}
c   Ta: (8s4p4d3f4g)  / [6s4p3d3f2g]  {311111/1111/211/111/31}
c   W:  (9s4p4d3f4g)  / [7s4p3d3f2g]  {3111111/1111/211/111/31}
c   Re: (9s4p4d3f4g)  / [7s4p3d3f2g]  {3111111/1111/211/111/31}
c   Os: (9s4p4d3f4g)  / [7s4p3d3f2g]  {3111111/1111/211/111/31}
c   Ir: (9s4p4d3f4g)  / [7s4p3d3f2g]  {3111111/1111/211/111/31}
c   Pt: (10s4p5d3f4g) / [8s4p4d3f2g]  {31111111/1111/2111/111/31}
c   Au: (9s4p4d3f4g)  / [8s4p3d3f2g]  {21111111/1111/211/111/31}
c   Hg: (12s4p3d2f2g) / [7s4p3d2f2g]  {6111111/1111/111/11/11}
c   Tl: (5s3p3d1f1g)  / [3s3p2d1f1g]  {311/111/21/1/1}
c   Pb: (5s3p3d1f1g)  / [4s3p2d1f1g]  {2111/111/21/1/1}
c   Bi: (5s3p3d1f1g)  / [4s3p3d1f1g]  {2111/111/111/1/1} 
c   Po: (5s3p3d1f1g)  / [4s3p2d1f1g]  {2111/111/21/1/1}
c   At: (5s3p3d1f1g)  / [3s3p2d1f1g]  {311/111/21/1/1}
************************************************************
c
      data pt5,pt75/0.5d+00,0.75d+00/
      data done/1.0d+00/
      data tol /1.0d-10/
 
      data pi32/5.56832799683170d+00/
      data pt187,pt6562 /1.875d+00,6.5625d+00/
c
      ng = 0
      igauss = 0
      ityp = 0
      nshell = 0
      ngauss = 0
      odone = .false.
      call vclr(eex,1,300)
c
c     ----- hydrogen to helium -----
c
      if (nucz .le. 2) then
          call ahlrichs_fit0(eex,ccs,ccp,ccd,nucz)
c
c     ----- lithium to neon -----
c
      else if (nucz .le. 10) then
          call ahlrichs_fit1(eex,ccs,ccp,ccd,ccf,nucz)
c
c     ----- sodium to argon -----
c
      else if (nucz .le. 18) then
          call ahlrichs_fit2(eex,ccs,ccp,ccd,ccf,ccg,nucz)
c
c     ----- potassium to zinc -----
c
      else if(nucz.le.30) then
          call ahlrichs_fit3(eex,ccs,ccp,ccd,ccf,ccg,nucz)
c
c     ----- gallium to krypton -----
c
      else if(nucz.le.36) then
          call ahlrichs_fit4(eex,ccs,ccp,ccd,ccf,ccg,nucz)
c
c     ----- rubidium to cadmium
c
      else if (nucz .le. 48) then
          call ahlrichs_fit5(eex,ccs,ccp,ccd,ccf,ccg,nucz)
c
c     ----- indium to xenon
c
c      else if (nucz .le. 54) then
          call ahlrichs_fit6(eex,ccs,ccp,ccd,ccf,ccg,nucz)
c
c     ----- past xenon does not exist
c
      else
        call caserr(
     +   'attempting to site ahlrichs function on invalid centre')
      endif
c
c     ----- loop over each shell -----
c
      ipass = 0
  210 ipass = ipass+1
      call ahlrichs_fitsh(nucz,ipass,ityp,igauss,ng,odone)
      if(odone) go to 220
c
c     ----- define the current shell -----
c
      nshell = nshell+1
      kstart(nshell) = ngauss+1
      intyp(nshell) = ityp
      kng(nshell) = igauss
      ngauss = ngauss+igauss
      k1 = kstart(nshell)
      k2 = k1+kng(nshell)-1
      do i = 1,igauss
        k = k1+i-1
        ex(k) = eex(ng+i)
        cs(k) = 0.0d0
        cp(k) = 0.0d0
        cd(k) = 0.0d0
        cf(k) = 0.0d0
        cg(k) = 0.0d0
        if(ityp.eq.16) cs(k) = ccs(ng+i)
        if(ityp.eq.17) cp(k) = ccp(ng+i)
        if(ityp.eq.18) cd(k) = ccd(ng+i)
        if(ityp.eq.19) cf(k) = ccf(ng+i)
        if(ityp.eq.20) cg(k) = ccg(ng+i)
      enddo
      go to 210
c
  220 continue
      return
      end
      subroutine ahlrichs_fitsh(nucz,ipass,itype,igauss,ng,odone)
      implicit none
c
      integer nucz,ipass,itype,igauss,ng
      logical odone
c
      integer igau, itypes, kind
      integer ind, mxpass, ityp, loop
c
      dimension igau(20,18), itypes(20,18)
      dimension kind(18)
      data kind/6, 4, 12, 14, 15, 13, 10, 11, 12, 11,
     +         11,18, 17, 18, 18, 13, 14, 15 /
c
      data igau /
     + 2,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     + 3,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     + 3,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,
     + 3,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,
     + 3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,
     + 3,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,
     + 8,1,1,1,1,3,1,3,1,1,0,0,0,0,0,0,0,0,0,0,
     + 8,1,1,1,1,2,1,1,4,1,1,0,0,0,0,0,0,0,0,0,
     + 8,1,1,1,1,4,1,1,4,1,1,1,0,0,0,0,0,0,0,0,
     + 8,1,1,1,1,4,1,1,4,1,1,0,0,0,0,0,0,0,0,0,
     +10,2,1,1,1,1,3,1,3,1,1,0,0,0,0,0,0,0,0,0,
     +10,2,1,1,1,1,1,1,1,1,2,1,1,1,1,1,3,1,0,0,
     +10,2,1,1,1,1,1,1,1,1,3,1,1,1,1,3,1,0,0,0,
     +10,2,1,1,1,1,1,1,1,1,2,1,1,1,1,1,3,1,0,0,
     +10,2,1,1,1,1,1,1,1,1,1,3,1,1,1,1,3,1,0,0,
     +10,2,1,1,1,1,1,1,1,1,3,1,1,0,0,0,0,0,0,0,
     +10,2,1,1,1,1,1,1,1,1,4,1,1,1,0,0,0,0,0,0,
     +10,2,1,1,1,1,1,1,1,1,1,4,1,1,1,0,0,0,0,0/
c
      data itypes /
     + 1,1,1,2,2,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     + 1,1,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     + 1,1,1,1,1,1,1,2,2,3,3,4,0,0,0,0,0,0,0,0,
     + 1,1,1,1,1,1,1,2,2,2,3,3,3,4,0,0,0,0,0,0,
     + 1,1,1,1,1,1,1,2,2,2,3,3,3,3,4,0,0,0,0,0,
     + 1,1,1,1,1,1,2,2,2,3,3,3,4,0,0,0,0,0,0,0,
     + 1,1,1,1,1,2,2,3,3,4,0,0,0,0,0,0,0,0,0,0,
     + 1,1,1,1,1,2,2,2,3,3,4,0,0,0,0,0,0,0,0,0,
     + 1,1,1,1,1,2,2,2,3,3,4,5,0,0,0,0,0,0,0,0,
     + 1,1,1,1,1,2,2,2,3,3,4,0,0,0,0,0,0,0,0,0,
     + 1,1,1,1,1,1,2,2,3,3,4,0,0,0,0,0,0,0,0,0,
     + 1,1,1,1,1,1,2,2,2,2,3,3,3,4,4,4,5,5,0,0,
     + 1,1,1,1,1,1,2,2,2,2,3,3,4,4,4,5,5,0,0,0,
     + 1,1,1,1,1,1,2,2,2,2,3,3,3,4,4,4,5,5,0,0,
     + 1,1,1,1,1,1,1,2,2,2,2,3,3,4,4,4,5,5,0,0,
     + 1,1,1,1,1,1,2,2,2,2,3,3,4,0,0,0,0,0,0,0,
     + 1,1,1,1,1,1,2,2,2,2,3,3,4,5,0,0,0,0,0,0,
     + 1,1,1,1,1,1,1,2,2,2,2,3,3,4,5,0,0,0,0,0/
c
c    pointers to the a1-fitting basis
c    set values for the current ahlrichs  shell
c
c  return igauss = number of gaussians in current shell
c         ityp   = 1,2,3,4,5 for s,p,d,f,g shell
c         ng     = offset in e,cs,cp,cd,cf and cg arrays for current shell
c
      ind=1
      if(nucz.gt. 1) ind=2
      if(nucz.gt. 2) ind=3
      if(nucz.gt. 4) ind=4
      if(nucz.gt. 8) ind=5
      if(nucz.gt. 9) ind=6
      if(nucz.gt.10) ind=7
      if(nucz.gt.12) ind=8
      if(nucz.gt.13) ind=9
      if(nucz.gt.17) ind=10
      if(nucz.gt.18) ind=11
      if(nucz.gt.20) ind=12
      if(nucz.gt.21) ind=13
      if(nucz.gt.26) ind=14
      if(nucz.gt.27) ind=15
      if(nucz.gt.30) ind=16
      if(nucz.gt.31) ind=17
      if(nucz.gt.35) ind=18
      if(nucz.gt.36) then
         call caserr('ahlrichs basis sets only extend to krypton')
      end if
c
      mxpass=kind(ind)
c
      if(ipass.gt.mxpass) odone=.true.
      if(odone) go to 100
c
       igauss = igau(ipass,ind)
       ng =0
       do loop = 1, ipass-1
        ng = ng + igau(loop,ind)
       enddo
       ityp = itypes(ipass,ind)
c
100   itype = ityp + 15
      return
      end
      subroutine ver_dft_readbasis(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/dft/readbasis.m,v $
     +     "/
      data revision /
     +     "$Revision: 5872 $"
     +      /
      data date /
     +     "$Date: 2009-03-02 12:31:23 +0100 (Mon, 02 Mar 2009) $"
     +     /
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
