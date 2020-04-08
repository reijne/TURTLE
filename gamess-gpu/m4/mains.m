c 
c  $Author: jmht $
c  $Date: 2015-03-13 21:56:17 +0100 (Fri, 13 Mar 2015) $
c  $Locker:  $
c  $Revision: 6317 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/mains.m,v $
c  $State: Exp $
c  
_IF(unix,t3d,cray,mips4)
c
c use unix memory allocation from C
c
_IF(chemshell,charmm,taskfarm)
c
c CHARMM/GAMESS-UK and ChemShell linked versions
c
c Currently we assumes that GAMESS-UK is built either 
c for serial execution or uses MPI
c
c init has a slightly different meaning in each case
c charmm - set up system but no energy calc
c chemsh - this is the first gamess point
c
_IF(chemshell,taskfarm)
      subroutine gamess(init,iret)
_ELSEIF(charmm_if4)
      subroutine gamess(init,iret,startup,iver)
_ELSEIF(charmm)
      subroutine gamess(init,iret,startup,lqmewd,iver)
_ENDIF
_ELSE
      program gamess
_ENDIF
      implicit none
      integer ierr
INCLUDE(common/mpidata)
INCLUDE(common/segm)
INCLUDE(common/vcore)
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
INCLUDE(common/maxlen)
INCLUDE(common/errcodes)
_IF(ga)
INCLUDE(common/disc)
_ENDIF

_IF(ma)
_IF(64bitpointers)
      common /mbc_dbl/          dbl_mb(2)
      double precision          dbl_mb

      integer*8 iq_8, type_8, size_8, itag_8

c iq64 is the offset relative to dbl_mb of the 
c start of gamess core. relative to this, all addresses
c allocated by MA should be representable in 32 bits
      integer*8 iq64
      common/pointer64/iq64
      logical ostat 
      character*16 tag
      integer itest, igmem_alloc
_ENDIF
_ENDIF

_IF(64bitpointers)
      integer*8 ioff
_ELSE
      integer ioff
_ENDIF
      integer idum
INCLUDE(common/restar)
INCLUDE(common/iofile)
      integer iret
_IF(dataout)
      logical opg_root
_ENDIF(dataout)
      logical opg_root
      integer ipg_nodeid
_IF(chemshell,charmm,taskfarm)
      integer init, i, ii, l, ioerr
      character*132 infile2
_IF(charmm)
INCLUDE(common/chmgms)
      logical startup,lqmewd
      integer iver
_IF(charmm_if4)
      logical oappdum
_ENDIF
_ENDIF
INCLUDE(common/work)
_IF(taskfarm)
INCLUDE(common/taskfarm)
      character*132 infile,outfile
INCLUDE(common/iotrp)
      logical oapp
_ELSE
INCLUDE(common/chemshell)
_ENDIF

_IF(chemshell)
      external chemshelldata
      logical oapp
_ENDIF
_IF(parallel)
      external fildat
_ENDIF
      external iodata
      external all_data
      external savscf
_ENDIF

_IF(_AND(ma,64bitpointers))
      logical ma_free_heap
      logical ma_alloc_get
_ENDIF

_IF(win95)
      integer iii,growhandles
c...    set control word for FPU
      call setcw
_IFN(ifc)
c...    allow 70 open files
      iii=growhandles(70)
_ENDIF
_ENDIF

_IF(vampir)
      call vttraceoff(ierr)
_ENDIF(vampir)
_IF(chemshell,charmm,taskfarm)
c
c this is done here, but with the PBEGINF/MPI_INIT missed out
c the controlling program will have done it, via the parallel_init
c routine below
c
      if (init .eq. 1)call pg_begin
_ELSE
      call pg_begin
_ENDIF
      call reset_time_periods
_IF(vampir)
      call vttraceon(ierr)
_ENDIF(vampir)

_IF(chemshell,charmm,taskfarm)

_IF(taskfarm)
      ocaser = .false.
      itaskret = 0
_ENDIF

_IF(charmm)
c
c  check that we are using consistent 
c  versions. This source supports version 4 and
c  version 5 (c34b1 release onwards)
c
_IF(charmm_if4)
      if (iver .ne. 4) then
         call caserr('CHARMM/GAMESS-UK Version Mismatch')
      endif
_ELSE
      if (iver .ne. 5) then
         call caserr('CHARMM/GAMESS-UK Version Mismatch')
      endif
_ENDIF
      iver = -iver
c
c Flag for CHARMM/GAMESS-UK initialisation
c
      qinigm = (init .eq. 1)
c these are controlled by the chm directive
      ochmdbg = .false.
      onoatpr = .false.
      oqmewd = lqmewd
c which QM charges to return      
      ichchty = CHM_MULLIKEN
      nchmlab = 0    
      if(startup)then
         charmm_cycle = 0
      endif

      eoff = 0.0d0
c
c Switch off append flag for first cycle only
c to create a new output file
c
      if(qinigm) oapp = .false.
c
_ELSE
      oapp = .false.
_ENDIF
      oswit = .false.

      if(opg_root())then
_IF(charmm_if4)
         call getfilenames(infile, outfile, 132 , oappdum)
_ELSE
         call getfilenames(infile, outfile, 132 )
_ENDIF
         ird=1 
         if(outfile .eq. "stdout")then
            iwr=6
         else
            iwr=2
            if(oapp)then
               open(unit=iwr,file=outfile,form='formatted',
_IF(fio_access_append)
     &              status='old',access='append',iostat=ierr)
_ELSE
     &              status='old',position='append',iostat=ierr)
_ENDIF
            if(ierr.ne.0 .and. opg_root())then
               write(6,*)
     &     'warning: open of output with append setting failed'
            endif
_IF(fio_use_endfile)
c may still be needed for some f77 compilers, hopefully redundant
               if(ierr.eq.0)then
                  ENDFILE (UNIT=iwr)
               endif
_ENDIF
            else
               ierr = 999
            endif

            if(ierr.ne.0)then
               open(unit=iwr,file=outfile,form='formatted',
     &              status='unknown',iostat=ierr)
            endif

            if(ierr.ne.0)then
               infile = 'couldnt open listing: '//outfile
               call gamerr(infile,
     &              ERR_NO_CODE, ERR_NO_CLASS, ERR_ASYNC, ERR_NO_SYS)
            endif
         endif

         open(unit=ird,file=infile,form='formatted',
     &        status='old',iostat=ierr)

         infile2=infile
         if(ierr .ne. 0) then
_IF(charmm)
c Try and remove any existing replica number
            l=132
            CALL TRIMA(infile2,l)
            infile = infile2
            DO ii = 1,3
               IF(infile2(l-ii:l-ii).eq.'_')THEN
                  infile= infile2(1:l-ii-1)
               ENDIF
            ENDDO
            if(infile2 .ne. infile) then

               CALL TRIMA(infile2,l)
               write(iwr,*)'Cant find input: ',infile2(1:l)
               CALL TRIMA(infile,l)
               write(iwr,*)'Trying: ',infile(1:l)

               open(unit=ird,file=infile,form='formatted',
     &              status='old',iostat=ierr)
            endif
_ENDIF
            if(ierr.ne.0)then
               outfile = 'couldnt open input: '//infile2
               call gamerr(infile,
     &              ERR_NO_CODE, ERR_NO_CLASS, ERR_ASYNC, ERR_NO_SYS)
            endif
         endif
      endif
_ENDIF

_IF(dataout)
c     jmht - hack to give a named output file rather then using stdout
      if(opg_root())then
        iwr=2
        open(unit=iwr,file='gamessuk.out',form='formatted',
     &        status='replace',iostat=ierr)
        if(ierr.ne.0)then
          call gamerr('could not open output gamess.out',
     &                 ERR_NO_CODE, ERR_NO_CLASS, ERR_ASYNC, ERR_NO_SYS)
        endif
      endif
_ENDIF(dataout)

_IF(parallel)
CMR   Make a Glue_uk_us enabled version
      call glue()
_ENDIF

       iret=0
       call initj(nmaxly)
_IF(_AND(ga,ma))

      call pg_ga_begin(nmaxly)
c
c if we are running with GA files, we want to
c open ed3/7 at the start to avoid memory fragmentation
c
      if( iwhere( 4  ) .EQ. 6 .Or. iwhere( 4 ) .EQ. 7 ) Then 
         call search(1,4)
      endif
      if( iwhere( 8  ) .EQ. 6 .Or. iwhere( 8 ) .EQ. 7 ) Then 
         call search(1,8)
      endif
c
      ioff = 0
c place values guard here
      qq(1) = 123.456d0
      qq(2) = 123.456d0

_ELSEIF(ma)
c Using the ma library for allocation but no GAs
c
c initialise the MA library
      call guk_ma_init(nmaxly)

c  ioff is not used - we pass in qq(1) to gmem_set which then
c  calculates the offset of dlb_mb(1) (the array relative to which
c  the ma library returns its indexes) to qq(1) and stores as iqoff
      ioff = 0

c place values guard here
      qq(1) = 123.456d0
      qq(2) = 123.456d0
_ENDIF

c
c print out Node information
c placed here so that GA data server processes don't get listed
c
      call pg_pproc
c
_IFN(ma)
c
c If we are not using MA tools for memory, allocate
c
      call mallocc(qq(1),qq(2),nmaxly,ioff,memhandle)
      if(ioff .eq. 0)call caserr('memory problem')
_ENDIF

      maxq=nmaxly


_IF(_AND(ma,64bitpointers))
c
c  The objective here is to pass an address close to
c  that returned by MA, so that the offsets can be 
c  stored in 32 bit quantities
c
      type_8 = 1013
      size_8 = 1
      tag='gamess init'
      itag_8 = 0
      iq_8 = 0
      
      ostat = ma_alloc_get(
     &     type_8,
     &     size_8,
     &     tag,
     &     itag_8,
     &     iq_8)
      if (.not.ostat) then
        write(iwr,*)"*** ma_alloc_get: size    = ",size_8
        write(iwr,*)"*** ma_alloc_get: address = ",iq_8
        call caserr("ma_alloc_get: could not get memory")
      endif

      iq64 = iq_8
      call gmem_set(dbl_mb(iq_8))
      call init_segm
      call driver(dbl_mb(iq_8))        
      ostat = ma_free_heap(itag_8)
_ELSE
      call gmem_set(qq(1 + ioff))
      call init_segm
      call driver(qq(1 + ioff))

_ENDIF
      call gmem_highwater(IGMEM_QUIET)
      call gmem_summarize('mains.m','gamess','after_driver',IGMEM_DEBUG)
_IF(parallel)
c parallel tests if requested
cc       call pg_test
_ENDIF

_IF(ma)
      if(qq(1) .ne. 123.456d0 .or. qq(2) .ne. 123.456d0)then
         call gamerr('Guard values at start of core corrupted',
     &        ERR_NO_CODE, ERR_NO_CLASS, ERR_SYNC, ERR_NO_SYS)
      endif
_ENDIF

c
c flag unconverged scf
c     write(6,*)'irest', irest
c
_IF(taskfarm)
      iret = itaskret
_ENDIF
      if (irest .ne. 0) iret = 1

_IF(ma)
cccc      call ma_summarize_allocated_blocks
_ENDIF
      call clenup2
_IF(chemshell,charmm,taskfarm)
c pg_end calls this has been moved into the charmm/chemshell parallel
c shutdown (see parallel_end)
_ELSE
      call pg_end(iret)
      write(iwr,*)'pg_end should not return'
_ENDIF
      end
c
      integer function icoresize()
_IF(ga)
c
c  memory will be assigned to MA tools 
c
_IF(ipsc)
c default for ipsc860 at DL
      icoresize = 1100000
_ELSEIF(t3e)
c default for 512 MByte nodes
c Cray T3E/1200 at Minnesota
c     icoresize = 48000000
c default for 256 MByte nodes
c Cray T3E/1200 at Manchester
      icoresize = 26000000
c default for 128 MByte nodes
c Cray T3E at Berlin and Orsay
c     icoresize = 9500000
_ELSEIF(t3d)
c default for Cray T3D at Edinburgh and ZIB
      icoresize = 3500000
_ELSEIF(bluegene)
c default for 512 MByte nodes (BlueGene/L)
      icoresize = 20000000
_ELSEIF(rs6000)
c default for Daresbury SP2/TN2 (64MBytes)
c     icoresize = 4000000
c default for Daresbury SP2/P2SC (256 MBytes)
c     icoresize = 20000000
c default for 512 MByte nodes (Minnesota)
c     icoresize = 48000000
c default for 1GByte nodes (HPCx)
      icoresize = 50000000
_ELSEIF(r10000)
c default for TERAS Origin 3000 (1GB RAM/CPU)
      icoresize = 25000000
_ELSEIF(altix)
c default for SGI Altix at Manchester (2GB RAM/CPU)
      icoresize = 60000000
_ELSEIF(hpux11)
c default for Exemplar (512MB RAM/CPU)
      icoresize = 20000000
_ELSEIF(ev5,ev6)
c default for Sierra systems (Columbus ES40s)
      icoresize = 20000000
_ELSEIF(ultra)
c default for Sun MP systems (EPCC)
      icoresize = 20000000
_ELSEIF(linux)
c default for Linux Clusters
      icoresize = 20000000
_ELSE
c default memory for worstation cluster server model
      nn = nnodes()
      id = nodeid()
      if(id .ge. nn/2)then
c memory server proces
         icoresize = 4000000
      else	
c gamess worker
         icoresize = 4000000
      endif
_ENDIF
_ELSEIF(ipsc)
      icoresize=500000
_ELSEIF(sv1)
c default for Cray SV1 at Groningen
      icoresize = 5000000
_ELSEIF(c90,j90)
c default for C90/J90
      icoresize=1000000
_ENDIF
_IF(newscf)
c consider allocations for MPI/scalapack code
_IF(rs6000,opteron,xeon)
c default for 1GByte nodes (HPCx,scarf)
      icoresize = 50000000
_ELSE
c other machines need specifying here
c at the moment take 250% of GA allocation
      icoresize = 5000000
_ENDIF
_ELSEIFN(ga,newscf)
c
c default for workstation version
      icoresize=4000000
_ENDIF
      end
_IF(charmm)
c
c this is intended to be called at the beginning of the charmm run
c (from machdep/paral2.src)
c
      subroutine parallel_init(iret)
_IF(ga)
c assume tcgmsg-mpi is in use, so this will call mpi_init
      call pbeginf
_ELSEIF(mpi,charmmpar)
_IF(i8)
      integer*4 status
_ELSE
      integer status
_ENDIF
      call mpi_init(status)
_ELSE
      call caserr('parallel_init unimplemented')
_ENDIF
      iret=0
      end
c
c this is intended to be called at the end of the charmm run
c (from machdep/paral1.src)
c
      subroutine parallel_end(iret)
      call pg_end(iret)
      iret=0
      end     
_ENDIF
c
_ELSE
c
c - static memory version
c
      program gamess
      implicit none
INCLUDE(common/segm)
INCLUDE(common/vcore)
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
INCLUDE(common/maxlen)

      integer icoresize
      integer ierr

_IF(vampir)
      call vttraceoff(ierr)
_ENDIF(vampir)
      call pg_begin
      call reset_time_periods
_IF(vampir)
      call vttraceon(ierr)
_ENDIF(vampir)
      call initj(nmaxly)

      if(nmaxly.gt.icoresize())call caserr(
     &     'static memory is smaller than requested')
c
      call gmem_set(qq(1))
      call init_segm
      maxq = nmaxly

c uncomment to check memory allocation code before job
c      call memtest(qq(1),nmaxly)

_IF(ga)
      call pg_ga_begin
_ENDIF
_IF(syslog)
      call logsys('Start Job'//char(0))
_ENDIF
_IF(glibc)
      call installsegfaulthandler()
      call installctrlchandler()
_ENDIF
      call driver(qq)
      call gmem_highwater(IGMEM_QUIET)
      call gmem_summarize('mains.m','gamess','after_driver',IGMEM_DEBUG)
_IF(glibc)
      call uninstallsegfaulthandler()
      call uninstallctrlchandler()
_ENDIF
_IF(syslog)
      call logsys('Finalize Job'//char(0))
_ENDIF
      call pg_end(0)
      end
c
c static memory icoresize
c
      integer function icoresize()
      implicit none
c
c explicit dimension for common/vcore/
c must be loaded first
c
      integer mxcore, ivoff
      parameter(mxcore=10000000)
      REAL qq
      common/vcore/qq(mxcore)
      common/vcoreoff/ivoff
      icoresize=mxcore
      ivoff=1
      end

_ENDIF      
      subroutine clenup2
      implicit none
      integer code
_IF(mpi)
      include 'mpif.h'
INCLUDE(common/mpidata)
      integer ierr
      integer icode, i
_ENDIF
_IF(taskfarm,chemshell,charmm)

_IF(taskfarm)
INCLUDE(common/taskfarm)
      character*132 infile,outfile
INCLUDE(common/iotrp)
      logical oapp
_ELSE
INCLUDE(common/chemshell)
_ENDIF

      logical opg_root
      external opg_root
INCLUDE(common/iofile)
_ENDIF
INCLUDE(common/vcore)
      integer idum
_IF(xml)
c Need to include to see if we need to close the xml file
INCLUDE(common/blocks)
INCLUDE(common/xmlout)
_ENDIF(xml)
_IF(taskfarm,chemshell,charmm)
c
c close input and output streams
c
_IF(charmm_if4)
      logical oapp
      call getfilenames(infile, outfile, 132 , oapp)
_ELSE
      call getfilenames(infile, outfile, 132 )
_ENDIF
      if(opg_root())then
         close(unit=1)
         close(unit=ird)
         if(outfile .ne. "stdout")then
            close(unit=iwr)
         endif
      endif
_ENDIF
_IF(xml)
      call xmlStopWriter()
_ENDIF
c
c flush standard out
      call flushn(6)
_IFN(ma,hpux11)
c release core memory
c note problem on C200
      call freec(qq(1), memhandle, idum)
_ENDIF
      end
c
_IF(chemshell)
      block data chemshelldata
      implicit none
INCLUDE(common/chemshell)
      data infile/"gamess.in"/
      data outfile/"gamess.out"/
      end
_ENDIF




**==io.f
      block data iodata
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
INCLUDE(common/sizes)
INCLUDE(common/work)
INCLUDE(common/workc)
c
      common/disc/isel,iselr ,iselw,irep,ichek,ipos(maxlfn),
     *  nam(maxlfn)
_IF(cray)
     * ,keep(maxlfn),keepf(3,maxfrt)
_ELSEIF(ipsc,t3d)
     * ,keep(maxlfn),istat(maxlfn),
     *  len(maxlfn),iop(maxlfn),iwhere(maxlfn),isync(maxlfn),
     *  iposf(maxfrt),keepf(maxfrt),
     *  istatf(maxfrt),lrec(maxfrt),nrec(maxfrt)
     * ,oform(maxfrt)
_ELSE
     * ,keep(maxlfn),istat(maxlfn),llll(maxlfn),iop(maxlfn),
     *  iposf(maxfrt),keepf(maxfrt),istatf(maxfrt),
     *  lrec(maxfrt),nrec(maxfrt),oform(maxfrt)
_ENDIF
c
c /bufa/
c
_IF(parallel)
_IF(ipsc)
c === tpf for 16 files
      parameter(nbuf=16)
_ELSE
      parameter(nbuf=maxlfn)
_ENDIF
      common/bufa/aaaaa(512*nbuf)
_ELSEIF(cray,ibm,vax)
_IF1(i)      common/bufa/aaaaa(8190)
_IF1(v)      common/bufa/aaaaa(4760)
c
_IF1()  c logic inverted, used to be
_IF1()  c
_IF1()  c_ELSEIF(convex,alliant,titan,sun,apollo,sgi,rs6000,dec,hp700,hpux11,ksr)
_IF1()
_ELSE
c
_IF(fortio)
      common/bufnew/bufa(512,maxblo,maxbuf)
      common/discne/idiscn(452)
_ENDIF
_IF(cio,fcio)
      common/bufa/aaaaa(8190)
_ENDIF
_ENDIF
_IF1(c)      common/bufnew/ibufa(512,maxblo,maxbuf)
_IF1(c)      common/discne/idiscn(maxlfn*44+maxbuf*8+9)
      common/blksiz/nsz,nsz512(8),nsznav,nsstat,
     *  nslen,nsmax
_IF1(iv)     * ,nsort
_IF1(c)     * ,nsort(40)
_IF(convex,alliant,titan,sun,apollo,sgi,ipsc,rs6000,dec,hp700,hpux11,ksr)
_IF(t3d)
     * ,nsort,nstid(2),oasync,lensrt
_ELSEIF(cio,fcio)
     * ,nsort
_ENDIF
_IF(fortio)
     *, nsort, nstid
_ENDIF
_ENDIF
      common/blksiz2/ns2,ns2512(10),npznav,npstat,
     * n2len,n2max
_IF1(iv)     *,nport
_IF1(c)     *,nport(40)
_IF(convex,alliant,titan,sun,apollo,sgi,ipsc,rs6000,dec,hp700,hpux11,ksr)
_IF(t3d)
     * ,nport,nptid(2),opsync,lenprt
_ELSEIF(cio,fcio)
     * ,nport
_ENDIF
_IF(fortio)
     *, nport, nptid
_ENDIF
_ENDIF
      common/blksi3/ns3,ns3512(8),npdnav,npdtat,
     * n3len,n3max
_IF1(iv)     *,ntort
_IF1(c)     *,ntort(40)
_IF(convex,alliant,titan,sun,apollo,sgi,ipsc,rs6000,dec,hp700,hpux11,ksr)
_IF(t3d)
     * ,ntort,nttid(2),otsync,lentrt
_ELSEIF(cio,fcio)
     * ,ntort
_ENDIF
_IF(fortio)
     *, ntort, nttid
_ENDIF
_ENDIF
INCLUDE(common/sector)
_IF1(iv)      common/linkio/mode(2),iselx(2),kp(2),mp(2),kpp(2),
_IF1(i)     * ipbn(2),nwor(10),istart(10),ib,nb,iposn,nseg(3),
_IF1(v)     * ipbn(2),nwor(6),istart(6),ib,nb,iposn,nseg(3),
_IF1(iv)     * ojunkk,ojunkl,ojunkm,ojunkn
c
_IF(GIGA_DUMP)
      data maxb,iblkla/999999999,4/
_ELSE
      data maxb,iblkla/9999999,1/
_ENDIF
      data jrec/-1/
      data jwidth/132/,oswit/.false./,oterm/.false./,oflush/.true./
      data nerr/999/
      data nline,noline/0,0/
      data iwidth/132/
      end
_IF(parallel)
c
c these initialisations have been removed from the parallel version
c of inita. 
c
      block data fildat
      implicit none
c
INCLUDE(common/sizes)
INCLUDE(common/disc)
INCLUDE(common/discc)
c
      data yed/'ed0','ed1' ,'ed2' ,'ed3' ,'ed4' ,'ed5' ,'ed6' ,
     *          'ed7' ,'ed8' ,'ed9' ,'ed10','ed11','ed12','ed13',
     *          'ed14','ed15','ed16','ed17','ed18','ed19',
     *          'mt0' ,'mt1' ,'mt2' ,'mt3' ,'mt4' ,'mt5' ,'mt6' ,
     *          'mt7' ,'mt8' ,'mt9' ,'mt10','mt11','mt12','mt13',
     *          'mt14','mt15','mt16','mt17','mt18','mt19'
     *         /
      data yft/
     *'ft01','ft02','ft03','ft04','ft05','ft06','ft07','ft08',
     *'ft09','ft10','ft11',
     *'ft12','ft13','ft14','ft15','ft16','ft17','ft18','ft19',
     *'ft20','ft21','ft22','ft23','ft24','ft25','ft26','ft27',
     *'ft28','ft29','ft30','ft31',
     *'ft32','ft33','ft34','ft35','ft36','ft37','ft38','ft39',
     *'ft40','ft41','ft42','ft43','ft44','ft45','ft46','ft47',
     *'ft48','ft49','ft50','ft51',
     *'ft52','ft53','ft54','ft55','ft56','ft57','ft58','ft59',
     *'ft60'/
      data zftn/
     *'ftn001','ftn002','ftn003','ftn004','ftn005','ftn006','ftn007',
     *'ftn008','ftn009','ftn010','ftn011',
     *'ftn012','ftn013','ftn014','ftn015','ftn016','ftn017','ftn018',
     *'ftn019','ftn020','ftn021','ftn022','ftn023','ftn024','ftn025',
     *'ftn026','ftn027','ftn028','ftn029','ftn030','ftn031',
     *'ftn032','ftn033','ftn034','ftn035','ftn036','ftn037','ftn038',
     *'ftn039',
     *'ftn040','ftn041','ftn042','ftn043','ftn044','ftn045',
     *'ftn046','ftn047','ftn048','ftn049','ftn050','ftn051',
     *'ftn052','ftn053','ftn054','ftn055','ftn056','ftn057',
     *'ftn058','ftn059','ftn060'/
      data oedpr/maxlfn*.false./
      data oftpr/maxfrt*.false./
      end
_ENDIF
**==all.f
      block data all_data
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
INCLUDE(common/sizes)
INCLUDE(common/drfopt)
_IF(drf)
      common/nottwi/obeen,obeen2,obeen3,obeen4
_ENDIF
      parameter (ncmax=65)
      parameter (mcen33 = 9 * maxat)
      parameter (ntr  = 1900)
      common/junk2/pcap(13,ntr)
INCLUDE(common/modj)
INCLUDE(common/modfov)
c
_IF(j90,c90)
      common/bufb/rrbufb(1500+mcen33),iibufb(1400+mxgrps)
_ELSE
      common/bufb/rrbufb(1500+mcen33),iibufb(200+mxgrps)
_ENDIF
      common/maskc/rmaskc(64)
INCLUDE(common/jinfo)
c     common/cslosc/closc(7),icslo(1084)
      common/craypk/icall(14+mxcsf*mxnshl+mxcsf+mxcsf),isymm(mxroot),
     *cofff(mxcsf*mxroot),t0,t1
c = 256 bfn / 100 centres
      common/blkin/ctxxx(maxat*12+3),itcxxx(4)
_IFN1(c)      common/junk/ccxyz(36776),iicxyz(101)
_IF1(c)      common/junk/ccxyz(22500),iicxyz(22500)
c     common/linkmc/
c = 256  bfn / 1000 centres
c
INCLUDE(common/fpinfo)
INCLUDE(common/rcombd)
c     dimension of ydd2
      common/direc/zrunt(50),yrunt(50),ydd1(100),ydd2(120),
     *             zdd3(70) ,zdd4(50)
INCLUDE(common/atmol3)
INCLUDE(common/scra7)
INCLUDE(common/hermit)
INCLUDE(common/wermit)
INCLUDE(common/blockc)
INCLUDE(common/phycon)
INCLUDE(common/mapper)
INCLUDE(common/gjs)
INCLUDE(common/filel)
INCLUDE(common/trntim)
c (common/timanb)
c
      REAL timbeg,  tmlast,  timsec2
      REAL wtimbeg, wtmlast, wtimsec
      REAL cpwa
      integer nodel
      common/timanb/timbeg,tmlast,timsec2(14),
     *wtimbeg,wtmlast,wtimsec(14),cpwa,nodel
c
c
      common/junkc/ylabel(26),ztype(35),zlist(300),znumb(100),zhead(90),
     * zheadr(20),zunit(33),ztittt(100),ztitl(120),zspin(30),
     * zsta,zhz,zgauss,zcont(10),zconta(maxat),zmul(2,maxorb)
INCLUDE(common/runlab)
INCLUDE(common/infob)
INCLUDE(common/nshel)
      common/symtry/invt(206)
INCLUDE(common/funct)
INCLUDE(common/foropt)
INCLUDE(common/vibrtn)
INCLUDE(common/transf)
INCLUDE(common/frame)
      common/czmat/ianz(maxnz,5),bla(maxnz,3),lbla(maxnz,3),nz,nnvar
INCLUDE(common/csubst)
INCLUDE(common/csubch)
INCLUDE(common/ijlab)
      common/shlnos/qq4,litt(21)
INCLUDE(common/shlt)
INCLUDE(common/misc)
INCLUDE(common/denss)
INCLUDE(common/indez)
_IF(vector)
INCLUDE(common/imisc)
_ELSE
      common/setint/ikn(24),bp01(16),ij12(4)
_ENDIF
INCLUDE(common/shlinf)
      common/pkfil/opkfl(3)
INCLUDE(common/bshell)
      common/root/suw(85),nroots(2)
      common/flips/ib123(13+4*mxprms)
INCLUDE(common/shlg70)
INCLUDE(common/picon)
      common/geom/axyz(34)
      common/pqgeom/abcd(24)
INCLUDE(common/ginf)
      common/type/jtype(2)
      common/pgeom/gep(8*mxprms*mxprms)
      common/qgeom/acx(6)
      common/maxc/cmax(mxprim),cmaxx(4*mxprms+2),ismlp(mxprms*mxprms+2)
INCLUDE(common/auxvar)
INCLUDE(common/astore)
INCLUDE(common/miscg)
      common/lt/labcdt(4)
INCLUDE(common/const)
      common/shllfo/shllf(12*mxprms+8)
INCLUDE(common/cntl2)
      common/limy/mpppp(23+mxroot)
      common/reso/iresco(400)
      common/prpspc/centx(3),mopic(3)
INCLUDE(common/datbas)
INCLUDE(common/datgue)
INCLUDE(common/datplt)
INCLUDE(common/datprp)
      common /miscop/ cmin(16,maxvar),
     2 spppp(31),fspa(3,maxat),fspaa(maxvar),ia(maxat+1)
INCLUDE(common/copt)
      common /sver  / tver  (18)
INCLUDE(common/linkan)
      common /dshlno/ dshlno(22)
INCLUDE(common/grad2)
      common /dshlnf/ dshlnf(24*mxprms + mxprms*mxprms + 15),idshln(4)
INCLUDE(common/dmisc)
      common /setd  / setd  (ncmax,10),setdd(6),isetd(34)
INCLUDE(common/ffq)
c...   for UHF natorbs
INCLUDE(common/unocas)
c
INCLUDE(common/inxblk)
INCLUDE(common/rtdata)
c
INCLUDE(common/restar)
      common/restrz/zgam(50)
      common/restrr/realg(30)
INCLUDE(common/restri)
INCLUDE(common/restrj)
      common/restrl/logg(50)
c
INCLUDE(common/cndx40)
INCLUDE(common/cndx41)
INCLUDE(common/field)
c ...
c     simons & jorgensen optimisation.
c ...
INCLUDE(common/zjorgs)
      common /jorgs/maxj,irecal,opowel,obfgs,obfgx,onrfo,ocut,
     *outjor,outdeb,maxsta
c
INCLUDE(common/cntl1)
INCLUDE(common/tol)
INCLUDE(common/iofile)
INCLUDE(common/utilc)
INCLUDE(common/prints)
INCLUDE(common/blocks)
c
c symmetry assignment
c
      common/atmol3o/osim1,osim2
INCLUDE(common/fsymas)
INCLUDE(common/seerch)
      common/l/lockk
INCLUDE(common/statis)
INCLUDE(common/chgcc)
      common/scfwfn/cicoef(2,12),f(25),alpha(325),beta(325),
     * no(10),nco,nseto,npair,ncores,ibm,old,nset,nopen,nope,noe(10)
INCLUDE(common/files)
INCLUDE(common/segm)
      common/scfopt/maxit,mconv,nconv,npunch,
     * accdi1,accdi2,odiis,icoupl(3),dmpcut,acccc(6),iaccv(2),
     * rshift(4),iextn,junks,dampoo(7)
INCLUDE(common/machin)
INCLUDE(common/timez)
_IF(parallel)
      common/scftim/tdiag(10)
_ENDIF
      common/tran/ilifc(maxorb),ntran(maxorb),itran(mxorb3),
     * ctran(mxorb3),otran,otri
INCLUDE(common/tranao)
INCLUDE(common/infoa)
INCLUDE(common/molsym)
      common/dnfnw/oden(3),iiden(3),vdnfnw(8,200),rdnfnw(1275),
     *  tdnfnw(50,3),ndnfnw(50),eltotxx,nomoxx
c
c     ------ scrf common blocks
c
INCLUDE(common/scrf)
      common/dplrfc/tele(3),tnuc(3),tmol(3),dtot
INCLUDE(common/gvalue)
c
c ----- Pisa solvation common blocks
c
      common/psscrf/ptspac,delect,iconv,opssc,
     +              facnrm,inkblk,itotbq,iscf,ispace,msecp
      common/surfgn/numbq(maxat),irejec(maxat),cavsze(maxat)
     * ,nocent,itmax,damper
c
INCLUDE(common/vdwrad)
INCLUDE(common/coval)
INCLUDE(common/indxsv)
      common/orb/nint(9)
c
c ***** hondo common blocks ****
c
      common/blk1/realx(510),irealx(4)
      common/sortpk/spk(3412)
      common/tapes/iko(750)
cjmht does not appear to be needed. It's now in the common directory
c     as an include file if it should be.
cjmht      common/funcon/bohr(9)
      common/cigrad/eds(130)
      common/specal/ipqrs(14)
INCLUDE(common/atmblk)
INCLUDE(common/jspec)
      common/maxlen/maxq,maxqq,omem(2*maxlfn+2)
INCLUDE(common/tdhfx)
INCLUDE(common/tdhf)
INCLUDE(common/ghfblk)
INCLUDE(common/qercom)
INCLUDE(common/qercmx)
      common/out/aj(512)
      common/scfblk/al(34)
INCLUDE(common/lattic)
INCLUDE(common/indsyx)
_IF1(c)c  === to placate unicos cft77
      common/symmos/imos(43,maxorb),immos(162)
INCLUDE(common/defunk)
INCLUDE(common/worksp)
INCLUDE(common/zlabs)
_IFN(charmm)
      common/small/dipd(3,maxat*3),pold(3,3,maxat*3) ,vcd(maxat*3,6)
_ENDIF
_IF1(c)c  === to placate unicos cft77
c     common/mcgrad/av(8)
      common/incrs/bc(60)
INCLUDE(common/crnams)
INCLUDE(common/crnamx)
_IF(cio,fcio)
c flag to avoid recursive calls to caserr
      common/iotrp/ocaser
_ENDIF

c *************************************
c ***** end of hondo common blocks ***
c *************************************
_IF(drf)
cdrf start
      common/c_of_m/pcm,qcm,rcm
cdrf end
_ENDIF
c
INCLUDE(common/periodic)
INCLUDE(common/harmon)
INCLUDE(common/overrule)
INCLUDE(common/zorac)
INCLUDE(common/canon_nat)
INCLUDE(common/fiddle_hesc)
_IF(xml)
INCLUDE(common/xmlin)
INCLUDE(common/xmlout)
_ENDIF(xml)
c
      data field/' '/
      data version/' '/
_IF(drf)
      data pcm,qcm,rcm/3*0.0d0/
      data obeen,obeen2,obeen3 /3*.false./
_ENDIF
c
      data ac6/'c6','c7','c8','c9','c10','cn'/
      data iangs/1,1,1,  2,2,2,2,2,   3,3,3,3,3,3,3,
     1   4,4,4,4,4,4,4,4,4/
      data iangc/1,1,1,  2,2,2,2,2,2,  3,3,3,3,3,3,3,3,3,3,
     1   4,4,4,4, 4,4,4,4,4, 4,4,4,4,4,4   ,0,0,0,0,0,0/
      data pnams/'q10','q11c','q11s','q20','q21c','q21s','q22c','q22s',
     1    'q30','q31c','q31s','q32c','q32s','q33c','q33s',
     2 'q40','q41c','q41s','q42c','q42s','q43c','q43s','q44c','q44s'/
      data pnamc/'x','y','z','xx','yy','zz','xy','xz','yz','xxx',
     1         'yyy','zzz','xyy','yzz','xxz','xzz','xxy','yyz',
     2         'xyz','xxxx','yyyy','zzzz','xxyy','yyzz','xxzz',
     3         'xxxy','xxxz','yyyz','xyyy','xzzz','yzzz','xyzz',
     4         'xyyz', 'xxyz','lx','ly','lz','px','py','pz'/
c
c     runlab
c
      data zbflab/maxorb*'        '/
c
c
c     set data in /rtdata/ used by integration routines
c
      data madd/17,17,18,18,19,19,20,20,21,21,22,22,
     &   23,23,24,24,25,25,26,26/
c
c     roots and weights generated using nag routine d01bcf
c     rhi contains squares of gauss-hermite roots of order
c     2,4,6.......18
c     whi are corresponding weights
c     rlow contains squares of gauss-legendre roots of order
c     2,4,6,......,18
c     wlow are corresponding weights
c     rfac contains ratios rlow/rhi
c
      data rhi/
     &0.50000000000000d+00,0.27247448713916d+01,0.27525512860841d+00,
     &0.55253437422632d+01,0.17844927485433d+01,0.19016350919349d+00,
     &0.85886356890120d+01,0.39269635013583d+01,0.13390972881264d+01,
     &0.14530352150332d+00,0.11807189489972d+02,0.64147297336620d+01,
     &0.30859374437175d+01,0.10745620124369d+01,0.11758132021178d+00,
     &0.15129959781108d+02,0.91242480375311d+01,0.51961525300544d+01,
     &0.25525898026682d+01,0.89830283456962d+00,0.98747014068481d-01,
     &0.18528277495852d+02,0.11989993039824d+02,0.75540913261018d+01,
     &0.43897928867310d+01,0.21805918884504d+01,0.77213792004278d+00,
     &0.85115442997594d-01,0.21984272840963d+02,0.14972627088426d+02,
     &0.10093323675221d+02,0.64831454286271d+01,0.38094763614849d+01,
     &0.19051136350314d+01,0.67724908764929d+00,0.74791882596818d-01,
     &0.25485979166099d+02,0.18046505467729d+02,0.12771825354869d+02,
     &0.87697567302685d+01,0.56944233429577d+01,0.33691762702432d+01,
     &0.16923950797932d+01,0.60323635708175d+00,0.66702230958194d-01/
      data whi/
     &0.88622692545276d+00,0.81312835447246d-01,0.80491409000551d+00,
     &0.45300099055090d-02,0.15706732032286d+00,0.72462959522439d+00,
     &0.19960407221138d-03,0.17077983007414d-01,0.20780232581489d+00,
     &0.66114701255824d+00,0.76404328552331d-05,0.13436457467813d-02,
     &0.33874394455481d-01,0.24013861108232d+00,0.61086263373533d+00,
     &0.26585516843565d-06,0.85736870435883d-04,0.39053905846291d-02,
     &0.51607985615884d-01,0.26049231026416d+00,0.57013523626248d+00,
     &0.86285911681258d-08,0.47164843550191d-05,0.35509261355194d-03,
     &0.78500547264582d-02,0.68505534223466d-01,0.27310560906425d+00,
     &0.53640590971209d+00,0.26548074740116d-09,0.23209808448653d-06,
     &0.27118600925380d-04,0.93228400862421d-03,0.12880311535510d-01,
     &0.83810041398987d-01,0.28064745852853d+00,0.50792947901661d+00,
     &0.78281997721164d-11,0.10467205795793d-07,0.18106544810936d-05,
     &0.91811268679300d-04,0.18885226302685d-02,0.18640042387545d-01,
     &0.97301747641316d-01,0.28480728566998d+00,0.48349569472546d+00/
      data rlow/
     &0.33333333333333d+00,0.74155574714581d+00,0.11558710999705d+00,
     &0.86949939491826d+00,0.43719785275109d+00,0.56939115967007d-01,
     &0.92215660849206d+00,0.63467747623464d+00,0.27618431387246d+00,
     &0.33648268067507d-01,0.94849392628837d+00,0.74833462838728d+00,
     &0.46159736149627d+00,0.18783156765245d+00,0.22163568807218d-01,
     &0.96346127870282d+00,0.81742801326687d+00,0.59275012773154d+00,
     &0.34494237942742d+00,0.13530001165525d+00,0.15683406607401d-01,
     &0.97275575129749d+00,0.86199133320339d+00,0.68426201565315d+00,
     &0.47237153700448d+00,0.26548115726894d+00,0.10183270400277d+00,
     &0.11675871940146d-01,0.97891421016235d+00,0.89222197421380d+00,
     &0.74931737854740d+00,0.57063582016217d+00,0.38177105339712d+00,
     &0.20977936861551d+00,0.79300559811486d-01,0.90273770256471d-02,
     &0.98320148322563d+00,0.91359942257427d+00,0.79673916319752d+00,
     &0.64594166107702d+00,0.47843096553757d+00,0.31334338332122d+00,
     &0.16953901896600d+00,0.63446670693112d-01,0.71868028362264d-02/
      data wlow/
     &0.10000000000000d+01,0.34785484513745d+00,0.65214515486255d+00,
     &0.17132449237917d+00,0.36076157304814d+00,0.46791393457269d+00,
     &0.10122853629038d+00,0.22238103445337d+00,0.31370664587789d+00,
     &0.36268378337836d+00,0.66671344308689d-01,0.14945134915058d+00,
     &0.21908636251598d+00,0.26926671931000d+00,0.29552422471475d+00,
     &0.47175336386513d-01,0.10693932599532d+00,0.16007832854335d+00,
     &0.20316742672307d+00,0.23349253653835d+00,0.24914704581340d+00,
     &0.35119460331752d-01,0.80158087159761d-01,0.12151857068790d+00,
     &0.15720316715819d+00,0.18553839747794d+00,0.20519846372130d+00,
     &0.21526385346316d+00,0.27152459411755d-01,0.62253523938648d-01,
     &0.95158511682493d-01,0.12462897125553d+00,0.14959598881658d+00,
     &0.16915651939500d+00,0.18260341504492d+00,0.18945061045507d+00,
     &0.21616013526485d-01,0.49714548894970d-01,0.76425730254889d-01,
     &0.10094204410629d+00,0.12255520671148d+00,0.14064291467065d+00,
     &0.15468467512627d+00,0.16427648374583d+00,0.16914238296314d+00/
      data rfac/
     &0.66666666666667d+00,0.27215603006790d+00,0.41992718021791d+00,
     &0.15736566546394d+00,0.24499839134006d+00,0.29942188282334d+00,
     &0.10736939391571d+00,0.16162041638918d+00,0.20624663818033d+00,
     &0.23157228207121d+00,0.80331896688366d-01,0.11665879303696d+00,
     &0.14958091987121d+00,0.17479825778177d+00,0.18849566212812d+00,
     &0.63679037660486d-01,0.89588534847421d-01,0.11407481291265d+00,
     &0.13513427777031d+00,0.15061737139021d+00,0.15882410982598d+00,
     &0.52501143266839d-01,0.71892563268415d-01,0.90581644583619d-01,
     &0.10760679357614d+00,0.12174729195090d+00,0.13188408619684d+00,
     &0.13717689210026d+00,0.44527932183337d-01,0.59590208781963d-01,
     &0.74238913033865d-01,0.88018358749513d-01,0.10021614971993d+00,
     &0.11011383507948d+00,0.11709216189089d+00,0.12069995716395d+00,
     &0.38578132580971d-01,0.50624727552267d-01,0.62382560132156d-01,
     &0.73655596266152d-01,0.84017456504924d-01,0.93002965172432d-01,
     &0.10017697462623d+00,0.10517713322195d+00,0.10774456465678d+00/
c
c
      data amps/24.0d0,30.0d0,35.0d0,39.0d0,43.0d0,47.0d0,51.0d0,55.0d0,
     c60.0d0/
      data ipoint/0,1,3,6,10,15,21,28,36/
      data ix/0,
     +        1,0,0,
     +        2,0,0,1,1,0,
     +        3,0,0,2,2,1,0,1,0,1,
     +        4,0,0,3,3,1,0,1,0,2,2,0,2,1,1/
      data iy/0,
     +        0,1,0,
     +        0,2,0,1,0,1,
     +        0,3,0,1,0,2,2,0,1,1,
     +        0,4,0,1,0,3,3,0,1,2,0,2,1,2,1/
      data iz/0,
     +        0,0,1,
     +        0,0,2,0,1,1,
     +        0,0,3,0,1,0,1,2,2,1,
     +        0,0,4,0,1,0,1,3,3,0,2,2,1,1,2/
c ...
      data  jsym/1,2,1,3,4,1,4,3,2,1,5,6,7,8,1,6,5,8,7,2,1,7,8
     *          ,5,6,3,4,1,8,7,6,5,4,3,2,1/
      data nix/0,
     +         1,0,0,
     +         2,1,1,0,0,0,
     +         3,2,2,1,1,0,0,1,0,0,
     +         4,0,0,3,3,1,0,1,0,2,2,0,2,1,1/
      data niy/0,
     +         0,1,0,
     +         0,1,0,2,1,0,
     +         0,1,0,1,2,3,2,0,1,0,
     +         0,4,0,1,0,3,3,0,1,2,0,2,1,2,1/
      data niz/0,
     +         0,0,1,
     +         0,0,1,0,1,2,
     +         0,0,1,1,0,0,1,2,2,3,
     +         0,0,4,0,1,0,1,3,3,0,2,2,1,1,2/
      data isx /1,0,0,1,1,0,1/
      data isy /0,1,0,1,0,1,1/
      data isz /0,0,1,0,1,1,1/
      data irrep/0,0,0,1,0,0,0,1,0,1,1,0,0,0,1,1,0,1,0,1,1,1,1,1/
      data idh /-1,1,-1,1,-1,1,-1,1,-1,-1,1,1,-1,-1,-1,-1,1,1,-1,-1,1,
     11,1,1,-1,-1,-1,-1,-1,1,-1,-1,1,-1,1,1,-1,-1,-1,-1,1,1,-1,-1,1,
     2-1,1,1,-1/
      data jdh /-1,1,-1,1,-1,-1,-1,-1,1/
c
c      list of all possible runtype's
c
      data zrunt/
     *'scf','trans','ci','saddle','optimize','optxyz',
     *'gradient','force','adapt','analyse','gf','tda',
     *'jorgense','branch','integral','nmr','********',
     *'irc','rpagrad','rpaoptim',
     *'index4', 'polariza', 'fockder','hessian','dipder',
     *'polder','magnetiz', 'intensit','test', 'lagrange',
     *'hyperpol','lazzeret','raman'  ,'infrared',6*'--------',
     *'mopac','response','montec',7*'********'/
      data yrunt/
     *'scf ','tran','ci  ','sadd','opti','optx','grad',
     *'forc','adap','anal','gf'  ,'tda' ,'jorg','bran',
     *'inte','nmr' ,'****','irc ','rpag','rpao',
     *'inde','pola','fock','hess','dipd',
     *'pold','magn','inte','test','lagr',
     *'hype','lazz','rama','infr',6*'----',
     *'mopa','resp','mont',7*'****'/
c
c-weights
c ams data deleted, all atomic masses should now be handled by
c amass_get function
c

_IF1()c
_IF1()c    after some discussion, it has been decided to replace the 
_IF1()c    averaged isotope at. weights with values for the 
_IF1()c    leading isotope .. should have minor impact on calc.
_IF1()c    vib. frequencies (chap2.vtab to be modified ..)
_IF1()      data ams/1.00797d0,4.0026d0,6.939d0,9.0122d0,10.811d0,
_IF1()     *12.01115d0,14.0067d0,15.9994d0,18.9984d0,20.183d0,
_IF1()     *22.9898d0,24.312d0,26.9815d0,28.086d0,30.9738d0,
_IF1()     *32.064d0,35.453d0,39.948d0,39.102d0,40.08d0,44.956d0,
_IF1()     *47.90d0,50.942d0,51.996d0,54.9380d0,55.85d0,58.9332d0,
_IF1()     *58.71d0,63.54d0,65.37d0,69.72d0,72.59d0,74.9216d0,
_IF1()     *78.96d0,79.909d0,83.80d0,85.47d0,87.62d0,88.905d0,
_IF1()     *91.22d0,92.906d0,95.94d0,99.0d0,101.07d0,102.905d0,
_IF1()     *106.4d0,107.87d0,112.40d0,114.82d0,118.69d0,121.75d0,
_IF1()     *127.60d0,126.9044d0,131.30d0/
_IF1()
_IF1() see new mass code 
_IF1()
_IF1()      data ams/
_IF1()     + 1.0078252d0, 4.0026036d0,7.016005d0,9.012186d0,11.0093051d0,
_IF1()     + 12.000000d0,14.0030744d0,15.9949149d0,18.998405d0,19.9924404d0,
_IF1()     + 22.989773d0,23.985045d0,26.981535d0,27.976927d0,30.973763d0,
_IF1()     + 31.972074d0,34.968854d0,39.962384d0,38.963714d0,39.962589d0,
_IF1()     + 44.955919d0,47.947948d0,50.943978d0,51.940514d0,54.938054d0,
_IF1()     + 55.93493d0,58.933189d0, 57.93534d0,62.92959d0,63.929145d0,
_IF1()     + 68.9257d0,73.921d0,74.921d0,79.916d0,78.914d0,83.911d0,
_IF1()     + 84.911d0,87.905d0,88.905d0, 89.904d0, 92.906d0,97.906d0,
_IF1()     + 98.906d0,101.903d0,102.904d0,105.903d0,106.905d0,113.903d0,
_IF1()     + 114.904d0,119.902d0,120.903d0,127.904d0,126.9043d0,131.904d0/
c
cjvl
cjvl  add 7,8,9 points hermite quadrature
cjvl   (hfff is dimensioned 7+8+9=24)
cjvl  allows up to g-functions in dma analysis
cjvl  M.Abramowitz, I.A.Stegun, Handbook of Mathematical Fucntions
cjvl  9th printing, (Dover publications, New York) (1972) page 924
c
      data h /0.0d0,
     +-.707106781186548d0,.707106781186548d0,
     +-1.22474487139159d0,0.0d0,
     + 1.22474487139159d0,-1.65068012388578d0,-.524647623275290d0,
     +.524647623275290d0,1.65068012388578d0,
     +-2.02018287045609d0,-.958572464613819d0,0.0d0,
     + .958572464613819d0,2.02018287045609d0,-2.350604973674d0,
     + -1.335849074014d0,-0.436077411928d0,0.436077411928d0,
     +  1.335849074014d0, 2.350604973674d0,
     + -2.651961356835d0, -1.673551628767d0,
     + -0.81628788285896d0,   0.0d0,  0.81628788285896d0,
     +  1.673551628767d0 , 2.651961356835d0,
     + -2.930637420257d0, -1.981656756696d0,
     + -1.157193712447d0, -0.3811869902073d0,
     +  0.3811869902073d0, 1.157193712447d0,
     +  1.981656756696d0,  2.930637420257d0,
     + -3.190993201781528d0, -2.266580584531843d0,
     + -1.468553289216668d0, -0.723551018752838d0, 0.0d0,
     +  0.723551018752838d0,  1.468553289216668d0,
     +  2.266580584531843d0,  3.190993201781528d0  /
c
      data w /1.77245385090552d0,
     +.8862269254528d0,.8862269254528d0,
     +.2954089751509d0,1.181635900604d0,
     + .2954089751509d0, 8.131283544725d-02,8.049140900055d-01,
     + 8.049140900055d-01,8.131283544725d-02,
     + 1.995324205905d-02,3.936193231522d-01,
     + 9.453087204829d-01,3.936193231522d-01,1.995324205905d-02,
     + 4.530009905509d-03,
     + 1.570673203229d-01, 7.246295952244d-01,7.246295952244d-01,
     + 1.570673203229d-01, 4.530009905509d-03,
     + 0.9717812450995d-03, 0.5451558281913d-01,
     + 0.4256072526101d0,   0.8102646175568d0, 0.4256072526101d0,
     + 0.5451558281913d-01, 0.9717812450995d-03, 0.1996040722114d-03, 
     + 0.1707798300741d-01, 0.2078023258149d0,   0.6611470125582d0,
     + 0.6611470125582d0,   0.2078023258149d0,   0.1707798300741d-01, 
     + 0.1996040722114d-03, 0.3960697726326d-04, 0.4943624275537d-02,
     + 0.8847452739438d-01, 0.4326515590026d0, 0.7202352156061d0,
     + 0.4326515590026d0,   0.8847452739438d-01,
     + 0.4943624275537d-02, 0.3960697726326d-04  /
c
      data ydd1 /
     *'safe','minb','bypa','nosy','zmat','ipri',
     *'nopr','angs','ang ','adap','****','main',
     *'mfil','seco','sfil','ffil','tfil','dmfi',
     *'twop','cifi','afil','maxb','loop','rloo',
     *'dump','inte','titl','prin','accu','supe',
     *'atom','geom','basi','gtos','stos','gaus',
     *'symm','char','norm','mult','rest','pseu',
     *'ecp' ,'form','fiel','comb','mapp','free',
     *'hfil','punc','denf','scrf','pssc','orie',
     *'symt','elec','dege','grou','vers','xfie',
     *'moro','reac','harm','zora','bqbq','blur',
     *'news','chm ','bqfo','cres','vb  ','depe',
     *'serv','ghos','disp','keys','xmli','xmlo',
     *'copq','cpwa','****','****','****','****',
     *'****','****','****','****','****','****',
     *'****','****','****','****','****','****',
     *'****','****','****','****'
     * /
c     do add keywords to ydd2, three other places
c     have to be changed: file m4/common/direc, 
c     this file, common/direc/zrun ...ydd2(120),
c     this file, data ... limit2 ...
      data ydd2 /
     *'maxc','thre','npun','damp','lock',
     *'swap','conv','leve','valu','xtol',
     *'step','lsea','tolm','tols','tans','minm',
     *'diis','orbs','noci','stat','simu','sort',
     *'augm','onep','hess','cive','cano','newt',
     *'line','lagr','jkop','ciac','supa','trac',
     *'casp','supe','pass','cipa','stop','open',
     *'popl','upda','nucl','prop','cent','sele',
     *'maxj','reca','powe','bfgs','bfgx','rfo ',
     *'cuto','optp','optd','mine','maxe','pref',
     *'inne','less',
     *'ente','runt','conf','loca','boys','lmo',
     *'grap','core','irc' ,'acti','dire','i.p.',
     *'ip'  ,'mull','mult','mcsc','dma', 'mrdc',
     *'full','mode','rpa' ,'sacc','mopa','vect',
     *'scft','nato','potf','omit','ccsd','qcis',
     *'dpa', 'save','hmat','nbo', 'cdft','dft' ,
     *'smea','aver','dies','vdwa','dlfi','mass',
     *'****','****','****','****','****','****',
     *'****','****','****','****','****','****',
     *'****','****','****','****','****','****' /
      data zdd3 /
     1'output','cards','vcd','gtens','field','skip','imptprin','open',
     2'ncoorb','keepfock','jacobi','oscf','grhf','contract','weights',
     7'uhf','canonica','open-sin','sections','density',
     8'mos','order','epsnes','ignore','salvage','shells','orbcore',
     b 'gauge','copyfile','reorder','restrict','irest2',
     c'irest3','irest4','irest5','mp2','seczero','closed',
     d'spdoc', 'spuoc','helfeygr','mprest','mp3','pwftol'
     e,'origin','fcm','pump2','ci','icflag'
     f,'frozen','minfc','adaptint','cpf',
     g'optorb','mpflag','nonuclea','frequenc',
     h'onelec','cutoff','chfconv','perturba','mixed',
     i'noprint','irestp','irest1','debug','mpasses','irest6','ordermos',
     j'********'/
cjvl  disabled zdd4 'check', 'intermol', 'ghost1', 'ghost2', 'supermol',
cjvl  disabled zdd4 'dispersi', 'multipol',
      data zdd4 / 'hfscf','gradient','gradone','optimize','force',
     1    '*****', '*****', '*****', '*****', '*****',
     2    '*****', '*****', 'index4', 'polariza',
     3    '       ','      ','magnetiz', 'saddle','integral',
     4    'fockder','hessian','dipder','polder','adapt',
     5    'scf','potentia','cichf','property','intensit',
     6    'optimise','energy','moguess','test','drt',
     7    'polarisa','magnetis','cisort',
     8    'lagrange','lazzeret','hyperpol',
     9     'allinfo','raman','infrared',7*'        '/
c
      data cpwa/0.0d0/
      data nbfs  /1,1,3,4,1,3,4, 6,1,3,4, 6,1,3, 6,1,3, 6,10,15,
     +   1,4, 6,4/
      data minf  /1,1,2,1,1,2,1, 5,1,2,1, 5,1,2, 5,1,2, 5,11,21,
     +   1,1, 5,1/
      data maxf  /1,1,4,4,1,4,4,10,1,4,4,10,1,4,10,1,4,10,20,35,
     +   1,4,10,4/
      data nangm /1,1,2,2,1,2,2, 3,1,2,2, 3,1,2, 3,1,2, 3, 4, 5,
     +   1,2, 3,2/
c     1s
      data stos1 /
     +1.24d0,1.69d0,
     +2.69d0,3.68d0,4.68d0,5.67d0,6.67d0,7.66d0,8.65d0,9.64d0,
     +10.61d0,11.59d0,12.56d0,13.53d0,14.50d0,15.47d0,16.43d0,17.40d0,
     +18.61d0,19.58d0,
     +20.56d0,21.54d0,22.53d0,23.52d0,24.50d0,25.49d0,26.47d0,27.46d0,
     +28.44d0,29.43d0,
     +30.42d0,31.40d0,32.39d0,33.37d0,34.36d0,35.34d0,
     +36.32d0,37.31d0,38.29d0,39.27d0,40.26d0,41.24d0,
     +42.22d0,43.21d0,44.19d0,45.17d0,46.15d0,47.14d0,
     +48.12d0,49.10d0,50.08d0,51.07d0,52.05d0,53.03d0/
c     2sp
      data stos2 /
     +0.00d0,0.00d0,
     +0.80d0,1.15d0,1.50d0,1.72d0,1.95d0,2.25d0,2.55d0,2.88d0,
     +3.48d0,3.90d0,4.36d0,4.83d0,5.31d0,5.79d0,6.26d0,6.74d0,
     +7.26d0,7.74d0,
     +8.22d0,8.70d0,9.18d0,9.66d0,10.13d0,10.61d0,11.09d0,11.56d0,
     +12.04d0,12.52d0,
     +12.99d0,13.47d0,13.94d0,14.40d0,14.87d0,15.34d0,
     +15.81d0,16.28d0,16.72d0,17.19d0,17.66d0,18.12d0,
     +18.59d0,19.05d0,19.51d0,19.97d0,20.43d0,20.88d0,
     +21.33d0,21.79d0,22.25d0,22.71d0,23.17d0,23.63d0/
c     3sp
      data stos3 /
     +10*0.0d0,
     +1.75d0,1.70d0,1.70d0,1.75d0,1.90d0,2.05d0,2.10d0,2.33d0,
     +2.75d0,3.01d0,
     +3.21d0,3.44d0,3.67d0,3.89d0,4.11d0,4.33d0,4.55d0,4.76d0,
     +4.98d0, 5.19d0,
     +5.26d0,5.58d0,5.90d0,6.22d0,6.54d0,6.86d0,
     +7.18d0, 7.49d0, 7.97d0, 8.21d0, 8.51d0, 8.82d0,
     +9.14d0, 9.45d0, 9.77d0,10.09d0,10.41d0,10.74d0,
     +11.08d0,11.39d0,11.71d0,12.03d0,12.35d0,12.66d0/
c     3d
      data stos4 /
     +20*0.0d0,
     +1.10d0,1.90d0,2.55d0,3.05d0,3.45d0,3.75d0,4.10d0,4.35d0,
     +4.60d0, 4.90d0,
     +5.26d0,5.58d0,5.90d0,6.22d0,6.54d0,6.86d0,
     +7.18d0, 7.49d0, 7.97d0, 8.21d0, 8.51d0, 8.82d0,
     +9.14d0, 9.45d0, 9.77d0,10.09d0,10.41d0,10.74d0,
     +11.08d0,11.39d0,11.71d0,12.03d0,12.35d0,12.66d0/
c      4sp
      data stos5 /
     +18*0.0d0,
     +1.43d0,1.36d0,
     +1.60d0,1.70d0,1.70d0,1.75d0,1.65d0,1.55d0,1.55d0,1.60d0,
     +1.60d0,1.90d0,
     +1.80d0,2.00d0,2.12d0,2.22d0,2.38d0,2.54d0,
     +3.02d0, 3.16d0, 3.29d0, 3.48d0, 3.67d0, 3.87d0,
     +4.05d0, 4.24d0, 4.41d0, 4.59d0, 4.76d0, 4.93d0,
     +4.65d0, 4.89d0, 5.12d0, 5.36d0, 5.59d0, 5.82d0/
c     4d
      data stos6 /
     +36*0.0d0,
     +0.00d0, 0.00d0, 1.40d0, 1.95d0, 2.40d0, 2.70d0,
     +3.00d0, 3.20d0, 3.45d0, 3.60d0, 3.75d0, 3.95d0,
     +4.65d0, 4.89d0, 5.12d0, 5.36d0, 5.59d0, 5.82d0/
c     5sp
      data stos7 /
     +36*0.0d0,
     +1.90d0, 1.80d0, 1.80d0, 1.90d0, 1.90d0, 1.95d0,
     +1.85d0, 1.75d0, 1.75d0, 1.80d0, 1.80d0, 2.10d0,
     +2.05d0, 2.15d0, 2.20d0, 2.28d0, 2.42d0, 2.57d0/
c
      data rleesf/1.20d+00,1.00d+00,1.03d+00,1.03d+00,1.03d+00,1.00d+00,
     1            0.99d+00,0.99d+00,1.00d+00,1.00d+00,1.00d+00,1.00d+00,
     2            1.00d+00,1.00d+00,1.00d+00,1.00d+00,1.00d+00,1.00d+00,
     3         12*1.00d+00,6*0.00d+00,
     4            1.15d+00,1.00d+00,1.12d+00,1.12d+00,1.12d+00,1.04d+00,
     5            0.98d+00,0.98d+00,1.00d+00,1.00d+00,1.00d+00,1.00d+00,
     6            1.00d+00,1.00d+00,1.02d+00,1.01d+00,1.01d+00,1.00d+00,
     7         12*1.00d+00,6*0.00d+00/
      data c11al/2*0.375d0,2*-0.125d0,2*0.125d0,2*-0.125d0,2*0.375d0,
     c2*0.125d0/
      data c11be/-0.25d0,2*0.25d0,4*-0.25d0,0.25d0,-0.75d0,0.75d0,0.25d0
     c,0.75d0/
      data c31al/2*0.125d0,2*0.625d0,2*0.375d0,2*0.625d0,2*0.125d0,
     c2*0.375d0/
      data c31be/-0.25d0,0.75d0,5*-0.25d0,-1.25d0,3*-0.25d0,0.75d0/
      data c33al/2*1.375d0,2*0.875d0,2*1.125d0,2*0.875d0,2*1.375d0,
     c2*1.125d0/
      data c33be/-0.75d0,2*-0.25d0,2*-0.75d0,-0.25d0,-0.75d0,-0.25d0,
     *-1.25d0,0.25d0,-0.25d0,-0.75d0/
      data ncom/1,6,3,6,3,6,6,7,10,10,3,10,10,9,1,3,3,1,1,6,3,3/
      data lin/0,1,0,0,2,0,0,1,1,0,3,0,0,2,2,1,0,1,0,1,4,0,0,3,3,1,0,1,
     +         0,2,2,0,2,1,1,5,0,0,4,4,1,0,1,0,3,3,2,0,2,0,3,1,1,2,2,1,
     +         6,0,0,5,5,1,0,1,0,4,4,2,0,2,0,4,1,1,3,3,0,3,3,2,1,2,1,2/
      data kin/0,0,1,0,0,2,0,1,0,1,0,3,0,1,0,2,2,0,1,1,0,4,0,1,0,3,3,0,
     +         1,2,0,2,1,2,1,0,5,0,1,0,4,4,0,1,2,0,3,3,0,2,1,3,1,2,1,2,
     +         0,6,0,1,0,5,5,0,1,2,0,4,4,0,2,1,4,1,3,0,3,2,1,3,3,1,2,2/
      data nin/0,0,0,1,0,0,2,0,1,1,0,0,3,0,1,0,1,2,2,1,0,0,4,0,1,0,1,3,
     +         3,0,2,2,1,1,2,0,0,5,0,1,0,1,4,4,0,2,0,2,3,3,1,1,3,1,2,2,
     +         0,0,6,0,1,0,1,5,5,0,2,0,2,4,4,1,1,4,0,3,3,1,2,1,2,3,3,2/
      data ngmx,nfmx,nomx,ncmx,ntmx,ntcol,icon,npmx/
     *1000,1024,1024,100,20,20,24*0,100/
      data cconv/
     *9.07618d0,17.7497d0,.171524d0,.324123d0,2.54154d0,
     *1.344911d0,.79198d0,.2800116d0,.711688d0,.14818d0,.14818d0/
c
c-weights
c atwt data deleted, all atomic masses should now be handled by
c amass_get function
c
_IF1()      data atwt/
_IF1()     *1.0078252d0,4.0026036d0,7.016005d0,9.012186d0,11.0093051d0,
_IF1()     *12.0d0,14.0030744d0,15.9949149d0,18.998405d0,19.9924404d0,
_IF1()     *22.989773d0,23.985045d0,26.981535d0,27.976927d0,30.973763d0,
_IF1()     *31.972074d0,34.968854d0,39.962384d0,38.963714d0,39.962589d0,
_IF1()     *44.955919d0,47.947948d0,50.943978d0,52.940514d0,54.938054d0,
_IF1()     *55.93493d0,58.933189d0,57.93534d0,62.92959d0,63.929145d0/
      data itps,maxnt/5,22/
      data cspin/
     *2.792780d0,-2.1275d0,1.085467d0,-0.392533d0,0.896167d0,
     *0.7024d0,0.2018d0,-0.37874d0,2.6288d0,-0.2206d0,
     *.739167d0,-.17106d0,.72828d0,-.5553d0,1.1317d0,
     *.2144333d0,.273943d0,-.185714d0,.1304667d0,-.188143d0,
     *.679486d0,-.15766d0,.7355714d0,-.1581333d0,.6888d0,
     *.0902d0,.6641429d0,-.2495667d0,.742d0,.1751d0/
c
      data bias,nvar/
     *5.0d0,.0001d0,.001d0,.1d0,.000001d0,.001d0,.0d0,.0d0,
     *0,10,10,50,0,0,0,0,1,1,1,0,1,1/
      data ncs/
     *0,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,
     *0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,
     *0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1/
      data nos/
     *1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
     *0,0,0,0,1,1,1,1,1,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
     *0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1/
      data npopul/
     *1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
     *0,0,0,0,1,2,3,4,5,0,0,0,1,2,3,4,5,0,0,0,0,0,0,0,0,0,0,0,0,0,
     *0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,2,3,4,5,6,7,8,9,0 /
      data ajm/
     *-1.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,
     * 0.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,
     *-1.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,
     * 0.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,
     *.0d0,.0d0,-1.6666666667d0,.1333333333d0,.0d0,.0d0,.0d0,.0d0,.0d0,
     *.0d0,.0d0,.0d0,.0d0,-.6666666667d0,-0.0666666667d0,.0d0,.0d0,.0d0,
     *.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,-.33333333333d0,-.1333333333d0,.0d0,
     *.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,-.1666666667d0,
     *-0.0166666667d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,
     *-0.0666666667d0,0.0053333333d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,
     *.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0/
      data ajm1/
     *-1.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,
     * 0.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,
     * 0.0d0,0.0d0,-1.6666667d0, 0.1333333d0,0.0d0,0.0d0,0.0d0,0.0d0,
     *0.0d0,0.0d0,0.0d0,
     * 0.0d0,0.0d0,-0.6666667d0,-0.0666667d0,0.0d0,0.0d0,0.0d0,0.0d0,
     *0.0d0,0.0d0,0.0d0,
     * 0.0d0,0.0d0,-0.3333333d0,-0.1333333d0,0.0d0,0.0d0,0.0d0,0.0d0,
     *0.0d0,0.0d0,0.0d0,
     * 0.0d0,0.0d0,-0.1666667d0,-0.0166667d0,0.0d0,0.0d0,0.0d0,0.0d0,
     *0.0d0,0.0d0,0.0d0,
     * 0.0d0,0.0d0,-0.0666667d0, 0.0053333d0,0.0d0,0.0d0,0.0d0,0.0d0,
     *0.0d0,0.0d0,0.0d0,
     * 0.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,
     *-1.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,
     *.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0/
      data ajm2/
     *.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,-1.8d0,.05714286d0,.05714286d0,
     *0.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,-.8d0,-.10612245d0,
     *.0367347d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,
     *-.4666666667d0,-.07891156d0,-.0154195d0,.0d0,.0d0,.0d0,.0d0,.0d0,
     *.0d0,.0d0,.0d0,-.3d0,-0.05d0,-0.05d0,.0d0,.0d0,.0d0,.0d0,.0d0,
     *.0d0,.0d0,0.0d0,-.2d0,-0.05714286d0,-.05714286d0,0.0d0,
     $0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     *-0.13333333d0,-0.022222222d0,-0.022222222d0,0.0d0,
     $0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     *-0.08571429d0,-0.01449396d0,-0.002832153d0,0.0d0,
     $0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     *-0.05d0,-0.00663265d0,0.002295918d0,0.0d0,
     $0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     *-0.022222222d0,0.000705467d0,0.000705467d0,0.0d0,
     * 0.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0/
      data ird,iwr,ipu/ 5,6,7/
      data ntapes/1,2,3,4,8/
      data ooprnt,opun7/.false.,.false./
      data omax,maxbrd/.true.,0/
      data ologf,logf/.false.,0/
      data ti,tx/ 0.0d0,0.0d0/
c     limit2= dimension of ydd2
      data limit1,limit2,limit3,mxtask/
     *   80,120, 120, 50/
      data prmoms/3*-1.0d0/
      data indx/0/
      data cov/
     * 0.47d0,3.80d0,    
     * 2.76d0,1.99d0, 1.62d0,1.33d0,1.23d0,1.14d0,0.95d0,3.80d0,
     * 3.42d0,2.85d0, 2.38d0,2.09d0,1.90d0,1.90d0,1.90d0,3.80d0,
     * 4.18d0,3.42d0,
     * 3.04d0,2.66d0,2.57d0,2.66d0,2.66d0,
     * 2.66d0,2.57d0,2.57d0,2.57d0,2.57d0,
     * 2.47d0,2.38d0,2.18d0,2.18d0,2.18d0,3.80d0,
     * 4.46d0,3.80d0,               
     * 3.42d0,2.94d0,2.76d0,2.76d0,2.57d0,
     * 2.47d0,2.57d0,2.66d0,3.04d0,2.94d0,
     * 2.94d0,2.76d0,2.76d0,2.66d0,2.66d0,3.80d0,
     * 4.94d0,4.09d0,
     * 3.71d0,3.52d0,3.52d0,3.52d0,3.52d0,
     * 3.52d0,3.52d0,3.42d0,3.33d0,3.33d0,
     * 3.33d0,3.33d0,3.33d0,3.33d0,
     * 3.33d0,2.94d0,2.76d0,2.57d0,2.57d0,
     * 2.47d0,2.57d0,2.57d0,2.57d0,2.85d0,
     * 3.61d0,3.42d0,3.04d0,3.61d0,3.61d0,3.80d0,
     * 4.94d0,4.09d0,
     * 3.71d0,3.42d0,3.42d0,3.33d0,3.33d0,
     * 3.33d0,3.33d0,3.23d0,3.13d0,3.13d0,
     * 3.13d0,3.13d0 /
     
c      data cov /0.57d0, 2.50d0,
c     * 2.53d0, 2.00d0, 1.53d0, 1.46d0, 1.42d0, 1.38d0, 1.34d0,3.00d0,
c     * 2.91d0, 2.69d0, 2.46d0, 2.23d0, 2.08d0, 1.93d0, 1.87d0,3.50d0,
c     * 2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0,2.00d0,
c     * 2.00d0, 3.00d0, 2.00d0, 3.00d0,
c     * 2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0,
c     * 64*2.0d0/
c
c ----- van der waals radii data in au, atomic number ordering
c
      data radvdw/
     *     2.26767d0,1.87083d0, 2*0.0d0, 4.10071d0, 3.495995d0, 
     *     2.83459d0, 2.645617d0, 2.551131d0, 3*0.0d0, 4.78101d0, 
     *     4.23299d0, 3.590481d0, 3.495995d0, 3.40151d0, 14*0.0d0,
     *     3.81725d0, 3.77945d0,3.77945d0,3.684967d0, 15*0.0d0,
     *     4.15740d0, 4.15740d0,  4.062913d0, 26*0.0d0,
     *     2.83459d0, 25*0.0d0 /
_IF1()c  ---  angstrom
_IF1()      data radvdw/1.2d0,4*0.0d0,1.85d0,1.5d0,1.4d0,1.35d0,5*0.0d0,
_IF1()     * 1.9d0,1.85d0,1.8d0,17*0.0d0,1.95d0,17*0.0d0,2.15d0,52*0.0d0/
_IF(parallel)
c imode flag - how to divide up grids over processors, set default (serial)
      common/pgraph/imode,idum(2)
      data imode/1/
_ENDIF
_IF(cio,fcio)
      data ocaser/.false./
_ENDIF
_IF(xml)
      data oxml/.false./
      data nixfile, ixfcoord, ixfaobas, ixfcdbas, ixfvect, 
     & ixfpop, ixfdens/7*0/
_ENDIF
      data nseto,ncores,npair,nset,nopen,nope/6*0/
      data nfblck/-1/
      data timsec,walsec,tstart,estart/102*0.0d0/
c
c common/periodic
c
      data zelem/
     $         'bq',
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
c  common/harmon
c
      data oharm,opharm/.false.,.false./,odepen/.false./
c
c  common/overrule
c
      data i_depen_check,i_global_diag/-1,-999/
      data n_depen_check,o_depen_print/0,.false./
c
c  common/zorac
c
      data ozora,onlys,onlyp,onlyd,onlyf,onlyg,oscalz,opzora,
     1     op2zora,ocontr,oioraz
     2    /5*.false.,         .true., 5*.false./
      data oint_zora,o1e_zora,osmall_zora/3*.false./
      data icoul_z,isexch_z/2,0/
      data niter_z/100/,opre_zora/.false./,int_zora/1/
      data onoext_z,oso/2*.false./
      data cspeed/137.0359895d0/
      data is_z/0/,igauge_z/0/,nat_z/2/
      data critat_z/1.0d-6/
      data oatint_z,oatscf_z,oscalatz/2*.false.,.false./
      data irest_z,numdu_z/2*0/
c
c  common canon_nat
c
      data o_canon,opr_canon/2*.false./
c
c  common datgue
c
      data natconf/0/
c
c  common fiddle_hesc
c
      data icall_h,ifreq_h/0,987654/
c
      end
      block data savscf
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c...  block data to initialise variables that could be better
c...  handled by save statements in subroutines
c...  unfortunately some compilers have trouble with those
c...  (variables have been checked not to cause interfeence
c...   in any of the other routines mentioned below)
c
c...  subroutines affected :
c     randd   (anala)   : u
c     adapt   (guess)   : idadap
c     chekpp  (input)   : ifirst
c     sinfo   (integc)  : ijato,klato,ijshlo
c     fj      (integb)  : common/fjsav/dp,dm,tab,hm,ep,em
c     fm      (integb)  : common/fmsav/t,est,ect
c     fsi     (integb)  : common/fsisav/errfs,dawfs
c     extrpd/extrpm(scf): common/extsav/odampr,oextpr

c
c...  machine specific uses are left alone in:
c     mxmb (FPS-MAX) (machscf) : ienter,ndots,lenofp,lenofl,ifun
c     ioerr (Apollo) (machscf) : igain
c     timef (atpghs) (machscf) : ienter,tst
c     get (io for ?) (machscf) : irbuf,iselb
c     oipsci (ipsc)  (machscf) : icount
c     proc21 (FPS)   (util3)   : arg
c     jksup1 (FPS)   (util3)   : arg
c
c..   in model are saved but not pre-initialised :
c     model (model) : nr,nr3,len173
c
c     for fixnag (util3/scf) seperate commons nagsav/nagsat are used
c
      common /saveco/ u,idadap,ifirst,ijato,klato,ijshlo
c
      parameter (maxf=10)
      character*20 fail(maxf)
      integer nfail(maxf),ifails(maxf)
      common/nagsat/ fail
      common/nagsav/ nfail,nf,ifails
c
      common/fjsav/dp,dm,tab,hm,ep,em
      common/fmsav/t,est,ect
      common/fsisav/errfs,dawfs
c
INCLUDE(common/extsav)
c
      data nfail/maxf*0/,nf/0/,ifails/maxf*0/
c
      end
      block data  fock_critical_data
      implicit none
c
INCLUDE(common/fock_critic)
      data ofipc,oshmem_inuse,nwfc/.false.,.false.,0/
      data ofipc_debug,ofipc_shmp,ofipc_nofock,ofipc_root/4*.false./
c
      end

      subroutine ver_mains(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/mains.m,v $
     +     "/
      data revision /"$Revision: 6317 $"/
      data date /"$Date: 2015-03-13 21:56:17 +0100 (Fri, 13 Mar 2015) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
_IF(hacked_mpich)
_IFN(charmm)
C----------------------------------------------------------------------
C     This is some stupid & dirty trick
      SUBROUTINE GETARG(I,S)      
      integer I
      character*(*) s
C      call getarg(i,s)
      call get_command_argument(i,s)
      return
      end
_ENDIF      
_ENDIF      
