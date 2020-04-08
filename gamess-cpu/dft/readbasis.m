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
INCLUDE(common/readbasis)
      integer nshl(maxtyp)
      character*8 ctyp(maxtyp)
      integer istart(maxtyp)
      integer iend(maxtyp)
      integer idim, fdim
c
      integer ibuff
      REAL  rbuff
      common/junkx/rbuff(10000),ibuff(10000)
c
      character*4 yform
      character*1 cform
      integer ntyp
INCLUDE(../m4/common/work)
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/infoa)
INCLUDE(../m4/common/runlab)

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
_IF(dft_basis_lib)
      else if(yform .eq. 'libr')then
c     
         nreq = 0
         more = .true.
         do while (more)

            call input
            call inpa(name)

            if(name(1:3) .eq. 'end') then
               more = .false.
            else
               call inpa(zlib)
cc               write(6,*)'jrec,jump', jrec,jump

               if(jrec .ne. jump)then
                  call inpa(zlab)
c     
c     use for a specified tag type
c     check this tag matches one of the atoms
c
                  hit = .false.
                  do iat=1,nat
                     if(ztag(iat) .eq. zlab)then
                        hit = .true.
                        iz = isubst(ztag(iat))
                     endif
                  enddo

                  if(.not. hit)
     &                 call caserr('unrecognised atom tag in jbas')

c
c  Now check if we already have a specification, if so 
c  replace it
c                     
                  hit = .false.
                  do ireq=1,nreq
                     if(reqz(ireq).eq.iz .and.
     &                    zlab .eq. reqt(ireq))then
                        hit=.true.
                        reqn(ireq) = name
                     endif
                  enddo
c
c Otherwise add it
c
                  if(.not.hit)then
                     nreq=nreq+1
                     reqz(nreq) = iz
                     reqn(nreq) = name
                     reqt(nreq) = zlab
                  endif

               else
c     
c     apply for all atom types that don't already have an 
c     explicit basis stored
c
                  do iat=1,nat
                     iz = isubst(ztag(iat))
                     hit=.false.
                     do ireq=1,nreq
                        if(reqz(ireq).eq.iz .and.
     &                       ztag(iat) .eq. reqt(ireq))hit=.true.
                     enddo
                     if(.not.hit)then
                        nreq=nreq+1
                        reqz(nreq) = iz
                        reqn(nreq) = name
                        reqt(nreq) = ztag(iat)
                     endif
                  enddo
               endif
            endif
         enddo
         if(nreq .gt. mxbasreq)call caserr('jbas limit')

         write(6,*)'the following have been requested from library :'
     $        ,zlib
         do ireq = 1,nreq
            write(6,*)reqz(ireq),reqn(ireq),reqt(ireq)
         enddo
c
c...     load basis sets from library
c
         call store_dgauss_basis(nreq, reqz, reqn, reqt,
     $        ntyp,ctyp,nshl,istart,
     &        ibuff,rbuff,idim,fdim,oprint,
     &        '/h/psh/gamess_dft/dft/data/dgauss.dat')
_ENDIF
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
INCLUDE(common/readbasis)
      integer     ntyp
      character*8 ctyp(maxtyp)
      integer     nshl(*)
      integer     istart(maxtyp)
      integer     index(*)
      REAL        value(*)
c
c...  Local variables
c
      integer ityp, ishl, iprim, jprim
      integer ik, ik1, ik2
      integer shell_type
      integer ix 
      REAL    facs, facp, facd, facf, facg, ee, ee1, c11, c12
      REAL    fac1, fac2, dum1, dum2, fac
c
      REAL pt5, pt75, pt187, pt6562
      data pt5,pt75,pt187,pt6562/0.5d0,0.75d0,1.875d0,6.5625d0/
c
      REAL pi32
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
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/iofile)
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

INCLUDE(common/readbasis)

      integer ntyp
      character*8 ctyp(maxtyp)
      integer nshl(*)
      integer istart(maxtyp)
      integer index(*)
      REAL  value(*)
      character *4 ytmp, y1, y2
      integer ix, iftop, itop
      REAL c1, c2, c3, ex
      logical onew
      integer i
      integer atom_type, shell_type
      character *1 format
      integer idim, fdim

INCLUDE(../m4/common/work)
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
INCLUDE(common/readbasis)
      integer nshl(*)
      character*8 ctyp(*)
      integer index(*)
      REAL value(*)
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

INCLUDE(common/readbasis)
_IF(taskfarm)
      character*10 ch_gnodeid
      character*132,dftopt_filename
      integer global_nodeid, len
_ENDIF

      integer ntyp
      integer istart(*)
      integer nshl(*)
      character*8 ctyp(*)
      integer index(*)
      REAL value(*)

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


_IF(taskfarm)
      write(ch_gnodeid,'(i10)')global_nodeid()

      do i=1,10
         if(ch_gnodeid(i:i).ne.' ') then
            len=i
            go to 999
         endif
      enddo

 999  dftopt_filename='options.dft.'//ch_gnodeid(len:10)

      open(50,file=dftopt_filename,form='formatted',
     &     status='unknown')
_ELSE
      open(50,file='options.dft',form='formatted',
     &     status='unknown')
_ENDIF

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

_IF(dft_basis_library)
      subroutine store_dgauss_basis(nreq, reqz, reqn, reqt,
     $     ntyp,ctyp,nshl,istart,
     &     index,value,idim,fdim,oprint,libname)
c
c
c input variables
c
c
c
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
c  for the first shell on an atom
c
c      0   number of contracted function
c      1   shell type 
c      2   address of floating values
c      3   address of next shell on this atom
c

      implicit none
INCLUDE(common/readbasis)

      integer ntyp
      logical oprint
      character*8 ctyp(maxtyp)
      integer nshl(*)
      integer istart(maxtyp)
      integer index(*)
      REAL  value(*)
      character *4 ytmp, y1, y2
      integer ix, iftop, itop
      REAL c1, c2, ex
      logical onew
      integer i
      integer atom_type, shell_type
      integer idim, fdim
c
c input args specifying which basis sets we need
c
      integer nreq, reqz(*)
      character*8 reqt(*), reqn(*)
      character*(*) libname
c
c communication with free format routine
c
      character*132 col,str
      integer ierr, nstr, nnum
      REAL anum
      common/fr/anum(50),nstr,nnum,ierr,col,str(50)
c
c  functions
c
      external isubst
      integer isubst
c
c local variables
c
      logical more, more2, hit
      integer nskips, nskipc, ireq, isave
      integer otyp, ifo, io, ishl, iprim
      integer iz
c
c initial pointers into storage arrays
c
      itop = maxtyp + 1
      iftop = 1
c
c initialisation
c
      ntyp=0
c
c Open the basis file
c      
      open(unit=97,form='formatted',file=libname,err=101)
c
c Main loop - load the header for a basis set
c
      more = .true.
      do while (more)

         call freef(97)

         if(ierr .ne. 0)then
            more = .false.
         else if(nstr+nnum .eq. 0)then
         else if(nstr .ge. 1 .and. str(1)(1:1) .eq. '#')then
         else

            hit = .false.

cc            write(6,*)'checking header',col(1:30)
            
            if(nnum .ne. 0 .or. nstr .ne. 3)then
               write(6,*)col
               call caserr('dgaus library error')
            endif

c
c check Z and basis name against list of requested
c library basis sets
c
            iz = isubst(str(1))
            do ireq = 1, nreq
               if(reqz(ireq) .eq. iz)then
                  if(str(2) .eq. reqn(ireq))then
                     hit = .true.
                     isave=ireq
                  endif
               endif
            enddo

            if(hit)write(6,*)'hit',col
c
c now scan this basis entry
c
            more2 = .true.
            nskips = -1
            do while (more2)

               call freef(97)

               if(ierr .ne. 0)then
                  call caserr('lib file error')
               else if(nstr+nnum .eq. 0)then
c ignore blanks
               else if(nstr .ge. 1 .and. str(1)(1:1) .eq. '#')then
c ignore comments
               else

                  if(nskips .eq. -1)then
                     nskips = nint(anum(1))
                     if(hit)then
c
c Create the entry
c
                        ntyp = ntyp + 1
                        atom_type = ntyp
                        ctyp(atom_type) = reqt(isave)
                        nshl(atom_type) = nskips
                        istart(atom_type) = itop
                        ix = istart(atom_type)

                     endif

c                     write(6,*)'shells to skip',nskips
                     nskipc = -1

                  else if(nskipc .eq. -1)then
c                     write(6,*)'coefs to skip',nskipc

                     nskipc = nint(anum(1))

                     if(hit)then

                        if(str(1) .eq. 's')then
                           shell_type = S
                        else if(str(1) .eq. 'P')then
                           shell_type = P
                        else if(str(1) .eq. 'D')then
                           shell_type = D
                        else if(str(1) .eq. 'F')then
                           shell_type = F
                        else if(str(1) .eq. 'G')then
                           shell_type = G
                        else if(str(1) .eq. 'SP')then
                           shell_type = SP
                        else if(str(1) .eq. 'SPD')then
                           shell_type = SPD
                        endif

                        if(ix+3.gt.idim)call caserr
     $                       ('basis i buffer too small')

                        index(ix + 0) = nskipc
                        index(ix + 1) = shell_type
                        index(ix + 2) = iftop
                        itop = itop + 4
                        if(nskips.ne.1)then
                           index(ix + 3) = itop
                           ix = itop
                        else
                           index(ix + 3) = 0
                        endif
                     endif
                  else
                     nskipc  = nskipc - 1

c     write(6,*)'add shell',ex,c1,iftop,itop,index(ix + 0)

                     if(hit)then

                        if(iftop+1.gt.fdim)call caserr
     $                       ('basis f buffer too small')
                        value(iftop) = anum(1)
                        iftop = iftop + 1

                        if(nnum .ge. 2)then
                           value(iftop) = anum(2)
                        else
                           value(iftop) = 1.0d0
                        endif
                        iftop = iftop + 1

                        if(shell_type .eq. SP .or. shell_type .eq. SPD
     $                       )then
                           if(iftop.gt.fdim)call caserr
     $                          ('basis f buffer too small')
                           if(nnum .ge. 3)then
                              value(iftop) = anum(3)
                           else
                              value(iftop) = 1.0d0
                           endif
                           iftop = iftop + 1
                        endif

                        if(shell_type .eq. SPD)then
                           if(iftop.gt.fdim)call caserr
     $                          ('basis f buffer too small')
                           if(nnum .ge. 4)then
                              value(iftop) = anum(4)
                           else
                              value(iftop) = 1.0d0
                           endif
                           iftop = iftop + 1
                        endif

                     endif

                     if(nskipc .eq. 0)then
                        nskips = nskips -1
                        if(nskips .eq. 0)then
                           more2 = .false.
                        else
                           nskipc = -1
                        endif
                     endif
                  endif
               endif
            enddo
c
c  In some cases, we may be storing the same basis for
c  more than one tag. clone the lists here
c            
            if(hit)then
               otyp = ntyp
               do ireq = isave+1,nreq
                  if(reqz(ireq) .eq. reqz(isave) .and.
     &                 reqn(ireq) .eq. reqn(isave)) then

                     write(6,*)'copying type'

                     ntyp = ntyp + 1
                     atom_type = ntyp
                     ctyp(atom_type) = reqt(ireq)
                     nshl(atom_type) = nshl(otyp)
                     istart(atom_type) = itop
                     ix = istart(atom_type)
                     io = istart(otyp)

                     do ishl = 1,nshl(atom_type)
                        itop = ix + 4
                        index(ix + 0) = index(io + 0)
                        index(ix + 1) = index(io + 1)
                        index(ix + 2) = iftop
                        index(ix + 3) = itop
                        if(ishl .eq. nshl(atom_type))index(ix + 3) = 0
                        ix = itop
                        ifo = index(io + 2)
                        io = io + 4
                        do iprim = 1, index(ix + 0)
                           if(iftop.gt.fdim)call caserr
     $                          ('basis f buffer too small')
                           value(iftop) = value(ifo)
                           iftop = iftop + 1
                           if(shell_type .eq. SPD)then
                              if(iftop.gt.fdim)call caserr
     $                             ('basis f buffer too small')
                              value(iftop) = value(ifo)
                              iftop = iftop + 1
                           endif
                           if(shell_type .eq. SP .or. shell_type .eq.
     $                          SPD)then
                              if(iftop.gt.fdim)call caserr
     $                             ('basis f buffer too small')
                              value(iftop) = value(ifo)
                              iftop = iftop + 1
                           endif
                        enddo
                     enddo
                  endif
               enddo
            endif
         endif
      enddo

      call printbasis(ntyp,ctyp,nshl,istart,index,value,oprint)

      return

 101  call caserr('problem opening library file')

      end

      integer function add_atom_type(label,
     &     istart,iend,nshl,ctyp)
c
c add_atom_type: 
c NOT FINISHED
c
      implicit none
INCLUDE(common/readbasis)

      character label*(*)
      character*8 ctyp(maxtyp)
      integer nshl(*),iend(*), istart(*)
      integer atom_type

      atom_type = 9999

      ctyp(atom_type) = label
      nshl(atom_type) = 0
      istart(atom_type) = 0
      iend(atom_type) = 0

      add_atom_type = atom_type
      end

      subroutine add_shell(atom_type, exp, coef, shell_type, ncont,
     &  index,value,idim,fdim,itop,iftop)

ccc      implicit none

c
c add_shell
c
c This routine appends a shell to the specified basis set
c definition
c
c NOT FINISHED
c

INCLUDE(common/readbasis)

      integer atom_type
      REAL exp(30)
      REAL coef(30,3)
      integer shell_type
      integer ncont
      integer index(*)
      integer nshl(maxtyp)
      integer istart(maxtyp)
      integer iend(maxtyp)
      REAL value(*)

      nshl(atom_type) = nshl(atom_type) + 1

      istart(atom_type) = itop
      ix = istart(atom_type)
c
c locate previous shell
c

c
c store shell info
c
      itop = ix + 4
      index(ix + 0) = index(io + 0)
      index(ix + 1) = index(io + 1)
      index(ix + 2) = iftop
      index(ix + 3) = 0
      if(ishl .eq. nshl(atom_type))index(ix + 3) = 0
      ix = itop

      do iprim = 1, index(ix + 0)

         if(iftop.gt.fdim)call caserr
     $        ('basis f buffer too small')
         value(iftop) = value(ifo)
         iftop = iftop + 1
         if(shell_type .eq. SPD)then
            if(iftop.gt.fdim)call caserr
     $           ('basis f buffer too small')
            value(iftop) = value(ifo)
            iftop = iftop + 1
         endif
         if(shell_type .eq. SP .or. shell_type .eq.
     $        SPD)then
            if(iftop.gt.fdim)call caserr
     $           ('basis f buffer too small')
            value(iftop) = value(ifo)
            iftop = iftop + 1
         endif
      enddo
      end

      subroutine freef(iunit)
c
c  * this copy of the free-format routine is used for reading outputs
c
      implicit REAL (a-h,o-z)
      character*132 col,str,s,s1,s2
      common/fr/anum(50),nstr,nnum,ierr,col,str(50)
      common/iounit/ninp,nout,ngms
      logical numeric,point,cont
      dimension itype(0:255)
c comment out all data statements except that relevant to your machine.
c if you machine uses a different character set you may get a lot of
c bad character errors (ierr = 1). set up a new data statement.
c ascii (vax,sun,macintosh)
      data itype/9*8,6,22*8,6,1,9,8*1,3,6,3,4,1,10*2,3*1,1,7*1,5,57*1
     $     ,129*8/
c ascii+128  (prime)
c      data itype/128*8,32*8,6,1,9,8*1,3,6,3,4,1,10*2,3*1,7,7*1,5,
c     & 57*1,1*8/
c ebcdic (ibm) still to add quote
c      data itype/64*8,6,9*8,1,4,1,1,3,1,1,9*8,6*1,3,1,8*8,1,6,4*1,
c     & 9*8,5*1,7,1,8,9*1,7*8,9*1,7*8,9*1,3*8,1,15*8,1,8,8,5*1,5,4*1,
c     & 6*8,10*1,6*8,1,8,8*1,6*8,10*2,6*8/
c *** internal character codes:
c     all print characters    1
c     except          0 to 9  2
c                       +-    3
c                     .       4
c                         e   5
c            <tab> <space> ,  6
c     control chars etc       8
c     double quotes        "  9
c
c *** error return codes
c     no error detected       0
c     unprintable character   1
c     too many (>50) numbers  2
c     too many (>50) literals 3
c     end of input file       4
c     other read error        5
c     bad numeric string      6
      maxchar=132
      nnum=0
      nstr=0
      ierr=0
      cont=.false.
      do 5 i=1,50
        str(i)=' '
5       anum(i)=0.d0
c **** this point to read a line
10    ipt=1
      if(iunit.ne.0)then
         read(iunit,1000,end=998,err=999)col
      endif
      cont=.false.
c
20    if(ipt.gt.maxchar)return
      it1=itype(ichar(col(ipt:ipt)))
      if(it1.eq.9)then
c quote
         l=1
 25      if(ipt+l.gt.maxchar)then
            write(6,*)'warning - missing double quote'
            goto 40
         endif
         it=itype(ichar(col(ipt+l:ipt+l)))
         if(it.ne.9)then
            l=l+1
            goto 25
         endif
         if(l.eq.1)then
            write(6,*)'empty quotes'
            goto 40
         endif
         nstr=nstr+1
         if(nstr.gt.50)then
            nstr=50
            ierr=3
            return
         endif
         str(nstr)=col(ipt+1:ipt+l-1)
         ipt=ipt+l
      else if(it1.le.5)then
         iepos=0
         numeric=.true.
         point=.false.
         if(it1.eq.4)point=.true.
         if(it1.eq.1)numeric=.false.
         if(it1.eq.5)numeric=.false.
         if(it1.eq.9)numeric=.false.
         l=1
30       if(ipt+l.gt.maxchar)goto 40
         it=itype(ichar(col(ipt+l:ipt+l)))
         if(it.eq.2)then
c           no action
         else if(it.eq.1)then
            numeric=.false.
         else if(it.eq.6)then
            goto 40
c         else if(it.eq.7)then
c            cont=.true.
c            goto 40
c         else if(it.eq.9)then
c            goto 40
         else if(it.eq.3.and.l.ne.iepos)then
            numeric=.false.
         else if(it.eq.4)then
            if(point)numeric=.false.
            if(iepos.ne.0)numeric=.false.
            point=.true.
         else if(it.eq.5)then
            if(iepos.ne.0)numeric=.false.
            if(.not.point)numeric=.false.
            if(l.le.1)numeric=.false.
            iepos=l+1
         else if(it.eq.8)then
            ierr=1
            return
         endif
         l=l+1
         goto 30
40       if(l.eq.1)then
            if(it1.eq.3)numeric=.false.
            if(it1.eq.4)numeric=.false.
         endif
         s=col(ipt:ipt+l-1)
         ipt=ipt+l-1
         if(numeric)then
            nnum=nnum+1
            if(nnum.gt.50)then
               nnum=50
               ierr=2
               return
            endif
            if(iepos.ne.0)then
c e format number
              il1=iepos-1
              il2=l-iepos
              s1=s(1:il1)
              s2=s(il1+2:80)
              s2(il2+1:il2+1)='.'
              read(s1,1020,err=997)part1
              read(s2,1020,err=997)part2
              if(part2+log10(abs(part1)+1.).gt.36.0d0)goto 997
              anum(nnum)=part1*10.0d0**part2
            else
c *** a dp is added if not present, and the number decoded
              if(.not.point)s(l+1:l+1)='.'
              read(s,1020,err=997)anum(nnum)
            endif
         else
            nstr=nstr+1
            if(nstr.gt.50)then
              nstr=50
              ierr=3
              return
            endif
c            if(iunit.eq.ninp)then
              do 45 i=1,l
                j=ichar(s(i:i))
c ibm case conversion
c                if(j.ge.193.and.j.le.201)then
c                  str(nstr)(i:i)=char(j-64)
c                else if(j.ge.209.and.j.le.217)then
c                  str(nstr)(i:i)=char(j-64)
c                else if(j.ge.226.and.j.le.233)then
c                  str(nstr)(i:i)=char(j-64)
c ibm ends
c ascii case conversion
                if(j.ge.65.and.j.le.90)then
                  str(nstr)(i:i)=char(j+32)
c ascii ends
                else
                  str(nstr)(i:i)=s(i:i)
                endif
45            continue
c            else
c              str(nstr)=s
c            endif
         endif
      else if(it1.eq.6)then
         goto 50
c      else if(it1.eq.7)then
c         cont=.true.
      else if(it1.eq.8)then
         ierr=1
         return
      endif
      if(cont)goto 10
50    ipt=ipt+1
      goto 20
c error traps
997   ierr=6
      nnum=nnum-1
      return
998   ierr=4
      return
999   ierr=5
      return
1000  format(a132)
1020  format(f80.0)
      end
_ENDIF
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
INCLUDE(common/readbasis)
c
      integer ntyp
      logical oprint
      character*8 ctyp(maxtyp)
      integer nshl(*)
      integer istart(maxtyp)
      integer index(*)
      REAL  value(*)
      integer ix, iftop, itop
      REAL c1, c2, ex
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
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/work)
INCLUDE(../m4/common/infoa) 
INCLUDE(../m4/common/infob)
INCLUDE(../m4/common/runlab) 
INCLUDE(../m4/common/datbas) 
INCLUDE(../m4/common/iofile) 
c
      integer ngauss, ishell, k, k1, loop, nati, natj
      integer nucz, jsubst, nshell, ityp
      integer jumpj,j,igauss
      integer locatc
c
      integer ibuff, intyp, kng, kstart
      REAL  rbuff, eex, cs, cp, cd, cf, cg
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
INCLUDE(../m4/common/sizes)
c
      REAL ex,cs,cp,cd
      integer intyp, nangm, kstart, kng
      integer nbfs, minf, maxf
      integer ngauss,nshell,nucz
c
      dimension ex(*),cs(*),cp(*),cd(*)
      dimension intyp(*),nangm(*)
      dimension kstart(*),kng(*)
      dimension nbfs(*),minf(*),maxf(*)
c
      REAL eex, ccs, ccp, ccd
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
_IFN1(civu)      call vclr(eex,1,200)
_IFN(civu)      call szero(eex,200)
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
INCLUDE(../m4/common/sizes)
c
      REAL ex,cs,cp,cd
      integer intyp, nangm, kstart, kng
      integer nbfs, minf, maxf
      integer ngauss,nshell,nucz,iwr
c
      dimension ex(*),cs(*),cp(*),cd(*)
      dimension intyp(*),nangm(*)
      dimension kstart(*),kng(*)
      dimension nbfs(*),minf(*),maxf(*)

      REAL eex, ccs, ccp, ccd
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
_IFN1(civu)      call vclr(eex,1,200)
_IFN(civu)      call szero(eex,200)
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
INCLUDE(../m4/common/sizes)
c
      REAL ex,cs,cp,cd
      integer intyp, nangm, kstart, kng
      integer nbfs, minf, maxf
      integer ngauss,nshell,nucz
c
      dimension ex(*),cs(*),cp(*),cd(*)
      dimension intyp(*),nangm(*)
      dimension kstart(*),kng(*)
      dimension nbfs(*),minf(*),maxf(*)
c
      REAL eex, ccs, ccp, ccd
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
_IFN1(civu)      call vclr(eex,1,200)
_IFN(civu)      call szero(eex,200)
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
INCLUDE(../m4/common/sizes)
c
      REAL ex,cs,cp,cd, cf, cg
      integer intyp, nangm, kstart, kng
      integer nbfs, minf, maxf
      integer ngauss,nshell,nucz
c
      dimension ex(*),cs(*),cp(*),cd(*),cf(*),cg(*)
      dimension intyp(*),nangm(*)
      dimension kstart(*),kng(*)
      dimension nbfs(*),minf(*),maxf(*)

      REAL eex, ccs, ccp, ccd, ccf, ccg
      common/blkin/eex(50),ccs(50),ccp(50),ccd(50),ccf(50),ccg(50)
c
      REAL pi32, pt5, pt75, pt187, pt6562
      REAL facs, facp, facd, facf, facg, ee, tol, done, fac
      REAL dums, dump, dumd, dumf, dumg
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
_IFN1(civu)      call vclr(eex,1,300)
_IFN(civu)      call szero(eex,300)
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
