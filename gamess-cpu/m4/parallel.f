
c  $Author: jmht $
c  $Date: 2013-05-27 14:27:44 +0200 (Mon, 27 May 2013) $
c  $Locker:  $
c  $Revision: 6287 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/parallel.m,v $
c  $State: Exp $
c  
c ======================================================================
c
c  ** pg_begin : initialise parallel processing
c
      subroutine pg_begin
      implicit none

c - other node parameters

      call pg_init

      return
      end

c
c ======================================================================
c
c subroutine pg_ga_begin : intialise GA storage
c
c  allocate memory for globals arrays. For machines without data servers
c  this is postponed till after initj to allow the memory to be set in
c  the input file
c
      subroutine  pg_ga_begin(size)
      implicit none
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
      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
      integer size
      igamem_totsize = 0
      return
      end
c
c this is to evade name clash between nnodes in common/nodinf
c and the tcgmsg function
c
      subroutine set_nodinf
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
c see also routine pg_init in parallel.m
c
      integer mpid, minode, mihost, ldim, nnodes, nodscr, maxnod
      common/nodinf/mpid,minode,mihost,ldim,nnodes,nodscr,maxnod(maxlfn)
c
      integer ipg_nodeid, ipg_nnodes
      mpid = 0
      minode = ipg_nodeid()
      mihost = 0
      nnodes = ipg_nnodes()
      nodscr = 1
      return
      end
c
c ========================================================================
c
c   pg_init: set parameters needed for  node process
c 

      subroutine pg_init
      implicit none
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
      integer icount_dlb, icount_slb
      common /par_dlb/ icount_dlb, icount_slb
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

      integer ipg_nnodes
c
      call set_nodinf
c
      call debugp('pg_init')
c
c use parallel diag for n > num nodes
      idpdiag = max(ipg_nnodes(),200)
c use PeIGS by default (except when the newscf F90 module
c is available with ScaLAPACK)
      ipdiagmode = IDIAG_PEIGS
c
c
c   Inversion algorithm
chvd
c     choose the old method by default for now until the new method 
c     has been suitably tested.
c     ipinvmode  = INV_CHOLESKY
      ipinvmode  = INV_DIAG
chvd
c
c use parallel linear algebra for n > 200
      idpdiis = 200
c
c    note that the ga orthog routine causes problems in
c    geometry optimisations .. this is now disabled in
c    start as a function of runtype until this is rationalised
       idporth = 200
c
       idpmult2 = 200
c aim for 40 chunks per scf
       ntchnk = 40
c limit on chunk size
       limchnk = -1
c dont do parallel tests
       iptest = 0
c
c I/O mode default is node zero i/o 
c screening out of parts of SCF to reduce comms
c is supressed pending further work
c
      ipiomode = IO_NZ
c
c dynamic lb counter
      icount_dlb = 0
c
c print parameter
c 
c 1 = root times
c 2 = all node times
c 3 = ga sizes
c
      iparapr = 1
c
      return
      end
c
c =================================================================
c
c  pg_pproc : print node information
c
c serial stub
      subroutine pg_pproc
      end
c
c========================================================================
c
c  ** pg_end(code)
c
c  orderly closedown of the GAMESS-UK process, returning the code
c  to the OS
c
c  also handles input and out units if needed
c
c  in standalone mode this is called at the end of the job
c  in the chemshell/charmm cases this is deferred until the 
c  shutdown of the controlling processes
c
      subroutine pg_end(code)
      implicit none
      integer code
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c    nothing specific needed for serial code
c
c Default behaviour is to use the C wrapper for exit
c
c      write(6,*)'Exiting code',code
      close(iwr)
      call exitc(code)
      end
c
c========================================================================
c
c  ** opg_root() : return .true. for root process
c
      logical function opg_root()
      opg_root = ipg_nodeid() .eq. 0
      return
      end

      logical function oroot()
      oroot = ipg_nodeid() .eq. 0
      return
      end
c
c=======================================================================
c
c
c========================================================================
c
c ** ipg_nodeid() : return index (0 - (nnodes-1)) for the 
c                   current process
c
      integer function ipg_nodeid()
      implicit none
      ipg_nodeid = 0
      return
      end

c
c========================================================================
c
c ** ipg_nnodes() : return number of nodes
c

      integer function ipg_nnodes()
      implicit none
      ipg_nnodes = 1
      return
      end
c
c========================================================================
c
c  static load balancing functions
c
      function oipsci()
      implicit real*8  (a-h,o-z)
      logical oipsci
       oipsci = .false.
       end
c
c========================================================================
c
       function  iipsci()
c
c...   initialise
c
       iipsci = 0
       return
       end


c
c=======================================================================
c
      subroutine maxabs(iin,iinout,n_mpi,itype_mpi)
      implicit none
c
c     This routine computes
c
c        dinout(i) = max(abs(dinout(i)),abs(din(i))), i=1,n
c
c     This routine is used in MPI_ALLREDUCE to implement the 'absmax'
c     function in the MPI implementation of pg_dgop.
c
      integer*4 n_mpi, itype_mpi
      integer*4 icomm_world_mpi
      integer iin(n_mpi), iinout(n_mpi)
      call caserr("unexpected use of subroutine maxabs")
      return
      end
c
c=======================================================================
c
      subroutine dmaxabs(din,dinout,n_mpi)
      implicit none
c
c     This routine computes
c
c        dinout(i) = max(abs(dinout(i)),abs(din(i))), i=1,n
c
c     This routine is used in MPI_ALLREDUCE to implement the 'absmax'
c     function in the MPI implementation of pg_dgop.
c
      integer*4 n_mpi
      real*8 din(n_mpi), dinout(n_mpi)
c
      integer i
      do i = 1, n_mpi
        dinout(i) = max(abs(dinout(i)),abs(din(i)))
      enddo
      end
c
c=======================================================================
c
      subroutine imaxabs4(iin,iinout,n_mpi)
      implicit none
c
c     This routine computes
c
c        iinout(i) = max(abs(iinout(i)),abs(iin(i))), i=1,n
c
c     This routine is used in MPI_ALLREDUCE to implement the 'absmax'
c     function in the MPI implementation of pg_igop.
c
      integer*4 n_mpi
      integer*4 iin(n_mpi), iinout(n_mpi)
c
      integer i
      do i = 1, n_mpi
        iinout(i) = max(abs(iinout(i)),abs(iin(i)))
      enddo
      end
c
c=======================================================================
c
      subroutine imaxabs8(iin,iinout,n_mpi)
      implicit none
c
c     This routine computes
c
c        iinout(i) = max(abs(iinout(i)),abs(iin(i))), i=1,n
c
c     This routine is used in MPI_ALLREDUCE to implement the 'absmax'
c     function in the MPI implementation of pg_igop.
c
      integer*4 n_mpi
      integer*8 iin(n_mpi), iinout(n_mpi)
c
      integer i
      do i = 1, n_mpi
        iinout(i) = max(abs(iinout(i)),abs(iin(i)))
      enddo
      end
c
c=======================================================================
c
c  ** subroutine pg_dgop : double precision global sum
c
c     Double Global OPeration.
c     x(1:n) is a vector present on each process. ggop 'sums'
c     x accross all nodes using the commutative operator op.
c     The result is broadcast to all nodes. Supported operations
c
c     include '+', '*', 'max', 'min', 'absmax', 'absmin'.
c
      subroutine pg_dgop(TYPE, X, N, OP)
      implicit none
      integer TYPE, N
      real*8 X(N)
      character*(*) OP     
      character*10 fnm
      character*7  snm
      data fnm/'parallel.m'/
      data snm/'pg_dgop'/
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
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
c Use simple binary tree
      integer ipg_nnodes
      if (ipg_nnodes().le.1) return
      call dbintree(x,n,op)
      return
      end
      integer function memreq_pg_dgop(N,OP)
      implicit none
c
c     work out the memory requirements for pg_dgop
c
      integer  icount, igmem_push_estimate, igmem_pop_estimate
      integer  igmem_incr
      external igmem_incr
      integer N, iwork
      character*(*) OP
      character*10 fnm
      character*14 snm
      data fnm/'parallel.m'/
      data snm/'memreq_pg_dgop'/
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
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
c Use simple binary tree
      integer ipg_nnodes
      memreq_pg_dgop = 0
      if (ipg_nnodes().le.1) return
      icount = igmem_push_estimate()
c     call dbintree(x,n,op)
      iwork = igmem_incr(n)
      call gmem_decr(iwork)
      memreq_pg_dgop = igmem_pop_estimate(icount)
      return
      end
c
c  Binary Tree implementation of global ops
c  Real*8 version
c
      subroutine dbintree(x,n,op)
      implicit none
      character*(*) OP
      integer n
      real*8 x(*)

      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
      real*8 qq
      integer ivoff, memhandle
      common/vcore/qq(2)
      common/vcoreoff/ivoff, memhandle
      integer n1, n2, n3, me, nproc, iwork, type, nret, ndfm
      integer igmem_alloc_inf, ipg_nodeid, ipg_nnodes
      external igmem_alloc_inf, ipg_nodeid, ipg_nnodes

      character*10 fnm
      character*8 snm

      data fnm/'parallel.m'/
      data snm/'dbintree'/

      type = 1021

      me = ipg_nodeid()
      nproc = ipg_nnodes()
      iwork = igmem_alloc_inf(n,fnm,snm,'iwork',IGMEM_DEBUG)
      n1 = 2*me+1
      n2 = 2*me+2
      n3 = (me-1)/2
   
      if (n2.lt.nproc) then
         call pg_rcv(type,qq(ivoff+iwork),n*8,nret,n2,ndfm,1) 
         if(nret.ne.n*8)call pg_errmsg('dbl gop msg err',-1)
         call ddoop(n, op, x, qq(ivoff+iwork))
      endif
      if (n1.lt.nproc) then
         call pg_rcv(type,qq(ivoff+iwork),n*8,nret,n1,ndfm,1)
         if(nret.ne.n*8)call pg_errmsg('dbl gop msg err',-1)
         call ddoop(n, op, x, qq(ivoff+iwork))
      endif
      if (me.ne.0) call pg_snd(type, x, n*8, n3, 1)
      call pg_brdcst(type, x, 8*n, 0)
   
      call gmem_free_inf(iwork,fnm,snm,'iwork')
      end

      subroutine ddoop(n, op, x, work)
      implicit none
      integer n
      real*8 x(n), work(n)
      character *(*) op
      integer i
c
c     real*8  Do Op ... do the operation for pg_dgop
c
      if (op .eq. '+') then
        do 10 i = 1,n
          x(i) = x(i) + work(i)
10      continue
      else if (op .eq. '*') then
	do 20 i = 1,n
	  x(i) = x(i) * work(i)
20      continue
      else if (op .eq. 'max') then
	do 30 i = 1,n
	  x(i) = max(x(i), work(i))
30      continue
      else if (op .eq. 'min') then
	do 40 i = 1,n
	  x(i) = min(x(i), work(i))
40      continue
      else if (op .eq. 'absmax') then
	do 50 i = 1,n
	  x(i) = max(abs(x(i)), abs(work(i)))
50      continue
      else if (op .eq. 'absmin') then
	do 60 i = 1,n
	  x(i) = min(abs(x(i)), abs(work(i)))
60      continue
      else
         call pg_errmsg('bad dgop',-1)
      endif
      end
c========================================================================
c
c  ** subroutine pg_igop : integer global sum
c
c     Integer Global OPeration.
c     x(1:n) is a vector present on each process. ggop 'sums'
c     x accross all nodes using the commutative operator op.
c     The result is broadcast to all nodes. Supported operations
c     include '+', '*', 'max', 'min', 'absmax', 'absmin'.
c
      subroutine pg_igop(TYPE, X, N, OP)
      implicit none
      integer TYPE, N
      integer X(N)
      character*(*) OP     
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
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
      character*10 fnm
      character*7  snm
      data fnm/'parallel.m'/
      data snm/'pg_igop'/
      integer ipg_nnodes
      if (ipg_nnodes().le.1) return
      call ibintree(x,n,op)
      return
      end
      integer function memreq_pg_igop(N, OP)
      implicit none
c
c     Work out the memory requirements for pg_igop
c     integer TYPE
c     integer X(N)
      integer N
      character*(*) OP     
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
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
      character*10 fnm
      character*14 snm
      data fnm/'parallel.m'/
      data snm/'memreq_pg_igop'/
      integer icount, iwork
      integer igmem_pop_estimate, igmem_push_estimate, igmem_incr
      integer lenwrd
      external lenwrd
      integer ipg_nnodes
      memreq_pg_igop = 0
      if (ipg_nnodes().le.1) return
      icount = igmem_push_estimate()
c     call ibintree(x,n,op)
      iwork=igmem_incr((n-1)/lenwrd() + 1)
      call gmem_decr(iwork)
      memreq_pg_igop = igmem_pop_estimate(icount)
      return
      end
c
c
c simple binary tree implementation
c Integer version
c
      subroutine ibintree(x,n,op)
      implicit none
      integer x(*), n
      character*(*) op

      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
      real*8 qq
      integer ivoff, memhandle
      common/vcore/qq(2)
      common/vcoreoff/ivoff, memhandle
      integer n1, n2, n3, me, nproc, liw, iwork, type, nret, ndfm

      integer igmem_alloc_inf, ipg_nodeid, ipg_nnodes
      external igmem_alloc_inf, ipg_nodeid, ipg_nnodes

      integer lenwrd
      external lenwrd

      character*10 fnm
      character*8  snm

      data fnm/'parallel.m'/
      data snm/'ibintree'/
c
      me = ipg_nodeid()
      nproc = ipg_nnodes()
      n1 = 2*me+1
      n2 = 2*me+2
      n3 = (me-1)/2
      type = 1022

      liw = 8/lenwrd()

      iwork = igmem_alloc_inf((n-1)/lenwrd() + 1,fnm,snm,'iwork',
     &                        IGMEM_DEBUG)

      if (n2.lt.nproc) then
         call pg_rcv(type,qq(ivoff+iwork),n*liw,nret,n2,ndfm,1)
         if(nret.ne.n*liw)call pg_errmsg('int gop msg err',-1)
         call ddoop(n, op, x, qq(ivoff+iwork))
      endif
      if (n1.lt.nproc) then
         call pg_rcv(type,qq(ivoff+iwork),n*liw,nret,n1,ndfm,1)
         if(nret.ne.n*liw)call pg_errmsg('int gop msg err',-1)
         call idoop(n, op, x, qq(ivoff+iwork))
      endif
      if (me.ne.0) call pg_snd(type, x, n*liw, n3, 1)
      call pg_brdcst(type, x, n*liw, 0)
      call gmem_free_inf(iwork,fnm,snm,'iwork')
      end

      subroutine idoop(n, op, x, work)
      implicit none
      integer n
      integer x(n), work(n)
      character *(*) op
      integer i
c
c     integer Do Op ... do the operation for pg_igop
c
      if (op .eq. '+') then
        do 10 i = 1,n
          x(i) = x(i) + work(i)
10      continue
      else if (op .eq. '*') then
        do 20 i = 1,n
          x(i) = x(i) * work(i)
20      continue
      else if (op .eq. 'max') then
	do 30 i = 1,n
	  x(i) = max(x(i), work(i))
30      continue
      else if (op .eq. 'min') then
	do 40 i = 1,n
	  x(i) = min(x(i), work(i))
40      continue
      else if (op .eq. 'absmax') then
	do 50 i = 1,n
	  x(i) = max(iabs(x(i)), iabs(work(i)))
50      continue
      else if (op .eq. 'absmin') then
	do 60 i = 1,n
	  x(i) = min(iabs(x(i)), iabs(work(i)))
60      continue
      else
         call pg_errmsg('bad igop',-1)
      endif
      end
c========================================================================
c
c ** pg_brdcst : byte-wise broadcast
c
      subroutine pg_brdcst(TYPE, BUF, LENBUF, IFROM)
      implicit none
      INTEGER TYPE    
      INTEGER LENBUF  
      INTEGER IFROM   
c no op
      integer buf(*)
      return 
      end
c
c========================================================================
c
c     subroutine pg_synch : synchronisation
c
      subroutine pg_synch(code)
      implicit none
      integer code
      return
      end
c
c ======================================================================
c
c  ** taskinfodata : load balancing initialisation
c
      block data taskinfodata
c
      integer ichi, ichl, ioff, nleft, itaskl
      integer g_nxtval ! global array handle for nxtval counter
      common /taskinfo/ g_nxtval,ichi,ichl,ioff,nleft,itaskl
c
      data ichi /1/
      data ichl /1/
      data nleft /0/
      data itaskl /0/
      end
c
c =====================================================================
c
      subroutine pg_eig_solve( nrow, a , q, e )
      call caserr
     $     ("pg_scalapack_eig: no GAs so no GA ScaLAPACK interfaces")
      end
c
c =====================================================================
c
c  ** pg_dlbcreate : create the global array for the nxtval counter and
c                    initialise the counter.
c
      subroutine pg_dlbcreate
      implicit none
      end
c
c =====================================================================
c
c  ** pg_dlbdestroy : destroy the global array for the nxtval counter
c
      subroutine pg_dlbdestroy
      implicit none
      end
c
c =====================================================================
c
c  ** pg_dlbchunk : set chunk size
c
      subroutine pg_dlbchunk(ichunk,prnt)
c
c  set up initial chunck size 
c  must be called with the same value from all nodes.
c  this determines the first index of the loop counter
c
      implicit none
c
      integer ichi, ichl, ioff, nleft, itaskl
      integer g_nxtval ! global array handle for nxtval counter
      common /taskinfo/ g_nxtval,ichi,ichl,ioff,nleft,itaskl
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      integer ichunk
      logical prnt
      ichi = ichunk
      ichl = ichunk
      end
c
c ====================================================================
c
c  ** pg_dlbreset : reset global and local counters
c
      subroutine pg_dlbreset
c
c     reset the task allocation code
c
      implicit none
c
      integer ichi, ichl, ioff, nleft, itaskl
      integer g_nxtval ! global array handle for nxtval counter
      common /taskinfo/ g_nxtval,ichi,ichl,ioff,nleft,itaskl
c
      integer icount_dlb, icount_slb
      common /par_dlb/ icount_dlb, icount_slb
      nleft=0
      itaskl=0
      icount_dlb = 0
      return
      end
c
c ====================================================================
c
c  ** pg_dlbfin : end dlb section
c
      subroutine pg_dlbfin
      implicit none
c
      integer ichi, ichl, ioff, nleft, itaskl
      integer g_nxtval ! global array handle for nxtval counter
      common /taskinfo/ g_nxtval,ichi,ichl,ioff,nleft,itaskl
c
      call pg_synch(999)
      call debugp('pg_dlbfin')
      end
c
c ====================================================================
c
c  ** ipg_dlbtask() :  get a global index
c
      integer function ipg_dlbtask()
      implicit none
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
      integer ichi, ichl, ioff, nleft, itaskl
      integer g_nxtval ! global array handle for nxtval counter
      common /taskinfo/ g_nxtval,ichi,ichl,ioff,nleft,itaskl
c
c
c serial version
c
      if(itaskl.eq.0)then
         itaskl = 1
      else
         itaskl = itaskl + 1
      endif
      ipg_dlbtask = itaskl
      end
c
c ====================================================================
c
c   ** ipg_dblpush :  get a global index
c
c            push the last index back for re-use, for use when a task 
c            is not used (as at the end of a load-balanced loop)
c
      subroutine pg_dlbpush
      implicit none
c
      integer ichi, ichl, ioff, nleft, itaskl
      integer g_nxtval ! global array handle for nxtval counter
      common /taskinfo/ g_nxtval,ichi,ichl,ioff,nleft,itaskl
c
      integer ipg_nnodes
      call debugp('pg_dlbpush *')
c
c tcgmsg and serial 
c
      nleft = nleft + 1
      itaskl = itaskl - 1
      end

c
      subroutine pg_dlbtest
      implicit none
      integer icount_dlb, icount_slb
      common /par_dlb/ icount_dlb, icount_slb
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      real*8 tot, tester
      integer ipg_nodeid
      integer ipg_dlbtask
      integer i, niter, next

      tot = 0.0d0
      niter=10

      call pg_dlbreset
      next = ipg_dlbtask()

      do i = 1,niter
         icount_dlb = icount_dlb+1
         if(icount_dlb .eq. next)then
            write(6,*)'task',next,'on node',ipg_nodeid()
            tot = tot + dble(icount_dlb)
c
c check push of indices
c
            next = ipg_dlbtask() 
            call pg_dlbpush

            next = ipg_dlbtask()
         endif
      enddo

      call pg_dgop(1010,tot,1,'+')
      tester= dble(niter*(niter+1)/2)
      if(dabs(tot - tester).lt.1.0d-10)then
         write(iwr,*)'test passed',tot
      else
         write(iwr,*)'test failed',tot,tester
      endif
      call pg_dlbpush

      write(iwr,*)'all done',ipg_nodeid()

      end
c
c ======================================================================
c
c ** pg_err : parallel error handling 
c   
c   should only get here on a parallel system with 
c   an asynchronous error
c
c   code : numerical error code
c
      subroutine pg_err(code)
      implicit none
      integer code
c
c serial implementations, generally use C exit routine
c
      call exitc(code)
      end
      subroutine pg_errmsg(s,i)
      implicit none
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      integer ipg_nodeid, i, me
      character s*(*)
      character tmp*100
      me=ipg_nodeid()
      write(6,*)'******* fatal error on node',me,' code =',i
      tmp = 'Fatal parallel error: '//s
      call caserr(tmp)
c      write(6,*)'       ',s
c      if (iwr .ne. 6)then
c         write(iwr,*)'******* fatal error on node',me,' code =',i
c         write(iwr,*)'       ',s
c      endif
c      call pg_err(i)
      return
      end
c     **just mpi-sync for now** 
c     **either using ga (tcgmsg-mpi) or straight mpi
c ======================================================================
c
c ** pg_snd_rcv : snd/rcv bytewise message
c
      subroutine pg_sndrcv(TYPES, BUFS, LENBUFS, NODE, 
     1                     TYPE, BUF, LENBUF, LENMES, NODESEL, NODEFROM)
      implicit none
      INTEGER TYPES, TYPE      
      INTEGER LENBUFS, LENBUF    
      INTEGER NODE      
      INTEGER BUFS(*), BUF(*)   
      INTEGER LENMES, NODESEL, NODEFROM, ifrom   

      INTEGER ierr
      integer*4 ierr_4,types_4,lenbufs_4,node_4,type_4,lenbuf_4
      integer*4 ifrom_4
      logical test_verb
c
c     fix up as we don't get the true length yet
c
      lenmes = lenbuf_4
c
      return
      end
c
c ======================================================================
c
c ** pg_snd : snd bytewise message
c
      subroutine pg_snd(TYPE, BUF, LENBUF, NODE, SYNC)
      implicit none
      INTEGER TYPE      
      INTEGER LENBUF    
      INTEGER NODE      
      INTEGER SYNC      
      integer buf(*)
      call caserr('pg_snd')
      return
      end
c
c ======================================================================
c
c ** pg_rcv : receive bytewise message
c
      SUBROUTINE pg_rcv(TYPE, BUF, LENBUF, LENMES, NODESEL, 
     &     NODEFROM, SYNC)
      implicit none
      INTEGER TYPE      
      INTEGER LENBUF     
c      BYTE BUF(LENBUF)   
      integer BUF(*)   
      INTEGER LENMES     
      INTEGER NODESEL    
      INTEGER NODEFROM   
      INTEGER SYNC       
      
      call caserr('pg_rcv')
      return
      end
c========================================================================
c
c ** pg_wait  : checks completion of asynchronous send and receive. See notes!
c       notes : - only asynchronous for tcgmsg if supported
c               - mpi checks communcation by status information,
c                 while tcgmsg checks by node. This may give
c                 problems when implementing.
      subroutine pg_wait(inode,imode)
      implicit none
      integer inode,imode
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
c trap not implemented
      call caserr('pg_wait')
      return
      end



c
c>>>>>>>>>>>>>>>>> old file below here <<<<<<<<<<<<<<<<<<<<<
c
c  should discontinue using these stubs, so as to 
c  help restrict the calls to the intel code, but the
c  parallel i/o "uses" them so heavily they are kept 
c
c - except dclock, as it is so widely used, now in machscf.m
c
c
c  utility routines
c
c test cpu availability - this seems to work quite well for 
c IBM SP2 nodes - poorer for HP700. probably some tuning
c of the sample period is required
c
      subroutine getload(fac)
      implicit none
      real*8 dumdum
      common/foolery/dumdum
      real*8 buf(3),w0,c0,w,c, fac
      integer i
      call gms_cputime(buf)
      c0 = buf(1)
      call walltime(w0)
      dumdum=0.0d0
      do i = 1,10000000
         dumdum=dumdum+sqrt(dble(i))
      enddo
      call gms_cputime(buf)
      c = buf(1)-c0
      call walltime(w)
      w = w -w0
      fac = c / w
      return
      end

      subroutine chtoi(i,s)
      integer i(*)
      character s*(*)
      ilen=len(s)
      do 10 ii = 1,ilen
        i(ii)=ichar(s(ii:ii))
10    continue
      return
      end
      subroutine itoch(i,s)
      integer i(*)
      character s*(*)
      ilen=len(s)
      do 10 ii = 1,ilen
        s(ii:ii)=char(i(ii))
10    continue
      return
      end

c  simple internal write fix (needed for apollo)
      subroutine intwrt(s,ip,iv,ic,otrunc)
      character cc(8)*(1), num(9)*(1), s*(*)
      logical otrunc
      otrunc=.false.
      do 10 iii=1,8
 10         cc(iii)='0'
      do 11 iii=1,9
 11      num(iii)=char(iii+ichar('1')-1)
      itmp=iv
      im=100000000
      itest=itmp/im
      if(itest.ne.0)then
         otrunc=.true.
         itmp=itmp-im*itest
      endif
      im=im/10
      do 20 iii=1,8
         itest=itmp/im
         if(itest.ne.0)then
            if(iii.le.8-ic)then
               otrunc=.true.
            else
               cc(iii)=num(itest)
            endif
         endif
         itmp=itmp-itest*im
 20      im=im/10
      do 30 iii=1,ic
 30      s(ip+iii:ip+iii)=cc(8-ic+iii)
      ip=ip+ic
      return
      end
c
c ======================================================================
c  control of error messages
c

      block data verbodat
      implicit none
      integer iverb, ilevel
      common/verbo/iverb(100),ilevel
      data iverb/100*0/
      data ilevel/0/
      end

      subroutine push_verb(i)
      implicit none
      integer iverb, ilevel
      common/verbo/iverb(100),ilevel
      integer i
      ilevel = ilevel + 1
      if(ilevel.eq.101)then
         call pg_errmsg('recursion gone mad',-1)
      endif
      iverb(ilevel)=i
      return
      end

      subroutine pop_verb
      implicit none
      integer iverb, ilevel
      common/verbo/iverb(100),ilevel
      ilevel = ilevel -1
      if(ilevel.le.0)then
         call pg_errmsg('pop_verb gone mad',-1)
      endif
      return
      end
c
      logical function test_verb(i)
      implicit none
      integer i
      integer iverb, ilevel
      common/verbo/iverb(100),ilevel
      if(ilevel.eq.0)then
         call pg_errmsg('bad initialisation or pop_verb error',-1)
      endif
      test_verb = i.le.iverb(ilevel)
      return
      end
c
      integer function get_verb()
      implicit none
      integer iverb, ilevel
      common/verbo/iverb(100),ilevel
      get_verb = iverb(ilevel)
      return
      end
c
c                        Memory allocation routines
c
c  These routines return offsets to dynamically allocated 
c  memory addresses. 
c
c  the addresses returned by igmem_alloc(size) may be used in one
c  of two ways:
c
c     i) they may be resolved as references into the core array
c        directly (core is as passed from mains to master):
c
c            subroutine junk(core,....)
c            real*8 core(*)
c            ....
c            i10 = igmem_alloc(n*n)
c            do i=1,n*n
c               core(i10 + i - 1) = funct(i)
c            enddo
c
c    ii) they may be considered as references into the common/vcore/qq
c        array when offset by ivoff (stored in common/vcoreoff/)
c
c        This offset is incorporated into the source during m4 processing
c        if the macro qq(ivoff+) is used.
c
c   gmem_set(q) 
c
c   gmem_set is called once at the start of the run to establish the
c   relationship between the core array (generally denoted q) which is
c   passed to the driver routines, and specific arrays in common 
c
c   The exact function on the routine depends on the underlying mechanism
c   used to obtain the memory.
c
c   GA-tools
c
c   When the GA tools are used the offset between the q(1) and the dbl_mb
c   array is stored (in real*8 words) as iqoff.
c
c   UNIX-memory:
c
c   When memory is being allocated dynamically from the OS, the 
c
c   When the gmem_routines are simply allocating memory from the GAMESS-UK 
c   stack, via getscm the offset between the qq array in vcore
c
      block data gmem_null_blockdata
      implicit none
c
c     This common block is to be used by the routines
c        gmem_null_blockdata
c        gmem_null_initialise
c        igmem_null
c        gmem_null_finalise
c     only. It holds the index in the main floating point array
c     of a chunk of memory of length 1 that can be passed into dummy
c     arguments when they are not being used. It needs to be converted
c     for use with integer arrays.
c
      integer igmem_null_pt
      integer*8 igmem_nan
      common/gmemnull/igmem_nan,igmem_null_pt
c
      integer null
      parameter(null=2**27+2**26+2**25)
      data igmem_null_pt/null/
      end

c     ****f* memory/gmem_null_initialise
c
c     NAME
c
c       gmem_null_initialise - initialises the null "pointer"
c
c     SOURCE
c
      subroutine gmem_null_initialise(q)
      implicit none
      real*8 q(*)
c
c     This common block is to be used by the routines
c        gmem_null_blockdata
c        gmem_null_initialise
c        igmem_null
c        gmem_null_finalise
c     only. It holds the index in the main floating point array
c     of a chunk of memory of length 1 that can be passed into dummy
c     arguments when they are not being used. It needs to be converted
c     for use with integer arrays.
c
      integer igmem_null_pt
      integer*8 igmem_nan
      common/gmemnull/igmem_nan,igmem_null_pt
c
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
      integer igmem_alloc_inf
      character *10 fnm
      character *20 snm
      integer*8 inan
      real*8 dnan
      equivalence(inan,dnan)
      data fnm/"parallel.m"/
      data snm/"gmem_null_initialise"/
      igmem_null_pt = igmem_alloc_inf(1,fnm,snm,"null",IGMEM_DEBUG)
      inan = -1
      q(igmem_null_pt) = dnan
      dnan = q(igmem_null_pt)
      igmem_nan = inan
      end
c     ******

      integer function igmem_null()
      implicit none
c
c     This common block is to be used by the routines
c        gmem_null_blockdata
c        gmem_null_initialise
c        igmem_null
c        gmem_null_finalise
c     only. It holds the index in the main floating point array
c     of a chunk of memory of length 1 that can be passed into dummy
c     arguments when they are not being used. It needs to be converted
c     for use with integer arrays.
c
      integer igmem_null_pt
      integer*8 igmem_nan
      common/gmemnull/igmem_nan,igmem_null_pt
c
      igmem_null = igmem_null_pt
      return
      end

      subroutine gmem_null_finalise(q)
      implicit none
      real*8 q(*)
c
c     This common block is to be used by the routines
c        gmem_null_blockdata
c        gmem_null_initialise
c        igmem_null
c        gmem_null_finalise
c     only. It holds the index in the main floating point array
c     of a chunk of memory of length 1 that can be passed into dummy
c     arguments when they are not being used. It needs to be converted
c     for use with integer arrays.
c
      integer igmem_null_pt
      integer*8 igmem_nan
      common/gmemnull/igmem_nan,igmem_null_pt
c
      integer igmem_alloc_inf
      character *10 fnm
      character *18 snm
      integer*8 inan
      real*8 dnan
      equivalence(inan,dnan)
      data fnm/"parallel.m"/
      data snm/"gmem_null_finalise"/
      dnan = q(igmem_null_pt)
      if (inan.ne.igmem_nan) then
         write(*,*)'*** inan      = ',inan
         write(*,*)'*** igmem_nan = ',igmem_nan
         call caserr("gmem_null_finalise: data at null pointer changed")
      endif
      call gmem_free_inf(igmem_null_pt,fnm,snm,"null")
      end

      subroutine gmem_set(q)
      implicit none
      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
      integer numchunk
      integer numestimate
      integer ichunksize
      integer itotsize
      integer imaxsize
      integer igamatsize
      integer igatotsize
      integer igamaxsize
      common/gmemcount/
     &     numchunk, numestimate,
     &     ichunksize(mxheap),  igamatsize(mxheap),
     &     itotsize(0:mxcount), igatotsize,
     &     imaxsize(0:mxcount), igamaxsize
      real*8 qq
      integer ivoff, memhandle
      common/vcore/qq(2)
      common/vcoreoff/ivoff, memhandle
      real*8 q(*)
      integer i
c
c
c
c now compute ivoff
c
      call gmem_c_pointer_diff(q(1),qq(1),ivoff)
      ogmem_alloced_all = .false.
      numheap = 0
      igmem_count = 0
      igmem_totsize(0) = 0
      igmem_maxsize(0) = 0
      igamem_totsize = 0
      igamem_maxsize = 0
      numcount = 0
      itotsize(0) = 0
      imaxsize(0) = 0
      numestimate = 0
      numchunk = 0
      do i=1,mxheap
         igamatsize(i)=0
      enddo
      return
      end


      subroutine gmem_usage(itotsize,num)
      implicit none
      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
      integer itotsize, num
      num = numheap
      itotsize = igmem_totsize(0)
      return
      end

      integer function igmem_overhead()
      implicit none
c
c     The number of extra words of memory used to allocate a chunk
c     of memory.
c
      igmem_overhead = 4
      return
      end
c
c     A bunch of utility routines for the purpose of calculating
c     memory usage of the code.
c
      integer function igmem_count_max()
      implicit none
      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
      integer numchunk
      integer numestimate
      integer ichunksize
      integer itotsize
      integer imaxsize
      integer igamatsize
      integer igatotsize
      integer igamaxsize
      common/gmemcount/
     &     numchunk, numestimate,
     &     ichunksize(mxheap),  igamatsize(mxheap),
     &     itotsize(0:mxcount), igatotsize,
     &     imaxsize(0:mxcount), igamaxsize
      igmem_count_max = imaxsize(0)
      return
      end

      integer function igmem_count_current()
      implicit none
      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
      integer numchunk
      integer numestimate
      integer ichunksize
      integer itotsize
      integer imaxsize
      integer igamatsize
      integer igatotsize
      integer igamaxsize
      common/gmemcount/
     &     numchunk, numestimate,
     &     ichunksize(mxheap),  igamatsize(mxheap),
     &     itotsize(0:mxcount), igatotsize,
     &     imaxsize(0:mxcount), igamaxsize
      igmem_count_current = itotsize(0)
      return
      end

      block data gmem_data_count
      implicit none
      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
      data igmem_count /0/
      data numheap /0/
      data igmem_size /mxheap*0/
      data igmem_totsize /mxcount*0,0/
      data igmem_maxsize /mxcount*0,0/
      data igamem_size /mxheap*0/
      data igamem_totsize /0/
      data igamem_maxsize /0/
      end

      block data gmem_data_estimate
      implicit none
      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
      integer numchunk
      integer numestimate
      integer ichunksize
      integer itotsize
      integer imaxsize
      integer igamatsize
      integer igatotsize
      integer igamaxsize
      common/gmemcount/
     &     numchunk, numestimate,
     &     ichunksize(mxheap),  igamatsize(mxheap),
     &     itotsize(0:mxcount), igatotsize,
     &     imaxsize(0:mxcount), igamaxsize
      data numchunk /0/
      data numestimate /0/
      data ichunksize /mxheap*0/
      data itotsize /mxcount*0,0/
      data imaxsize /mxcount*0,0/
      data igamatsize /mxheap*0/
      data igatotsize /0/
      data igamaxsize /0/
      end

      integer function igmem_incr(nsize)
      implicit none
      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
      integer numchunk
      integer numestimate
      integer ichunksize
      integer itotsize
      integer imaxsize
      integer igamatsize
      integer igatotsize
      integer igamaxsize
      common/gmemcount/
     &     numchunk, numestimate,
     &     ichunksize(mxheap),  igamatsize(mxheap),
     &     itotsize(0:mxcount), igatotsize,
     &     imaxsize(0:mxcount), igamaxsize
      integer igmem_overhead
      integer nsize, i
c
      if (numchunk.lt.0) then
         call caserr('igmem_incr: numchunk corrupted')
      endif
      if (numchunk.ge.mxheap) then
         write(*,*)'Too many memory chunk: you want to use ',
     &             numchunk+1,' but the maximum is ',mxheap
         call caserr('Too many memory chunks')
      endif
      numchunk = numchunk + 1
      ichunksize(numchunk) = nsize + igmem_overhead()
      do i=0,numestimate
         itotsize(i) = itotsize(i) + ichunksize(numchunk)
         imaxsize(i) = max(imaxsize(i),itotsize(i))
      enddo
      igmem_incr = numchunk
c
      return
      end

      subroutine gmem_decr(ichunk)
      implicit none
      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
      integer numchunk
      integer numestimate
      integer ichunksize
      integer itotsize
      integer imaxsize
      integer igamatsize
      integer igatotsize
      integer igamaxsize
      common/gmemcount/
     &     numchunk, numestimate,
     &     ichunksize(mxheap),  igamatsize(mxheap),
     &     itotsize(0:mxcount), igatotsize,
     &     imaxsize(0:mxcount), igamaxsize
      integer igmem_overhead
      integer ichunk, i
c
      if (ichunk.ne.numchunk) then
         write(*,*)'Non stack memory: you want to decrease ',ichunk,
     &             'but the current chunk is ',numchunk
         call caserr('Non stack memory counting')
      endif
      do i=0,numestimate
         itotsize(i) = itotsize(i) - ichunksize(numchunk)
      enddo
      numchunk = numchunk - 1
      if (numchunk.lt.0) then
         call caserr('gmem_decr: numchunk corrupted')
      endif
c
      return
      end
c
c **********************************************************************
c * igmem_push_estimate: create a new memory estimate counter
c *
c *   Return the handle of a new memory estimate counter
c *
c **********************************************************************

      integer function igmem_push_estimate()
      implicit none
      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
      integer numchunk
      integer numestimate
      integer ichunksize
      integer itotsize
      integer imaxsize
      integer igamatsize
      integer igatotsize
      integer igamaxsize
      common/gmemcount/
     &     numchunk, numestimate,
     &     ichunksize(mxheap),  igamatsize(mxheap),
     &     itotsize(0:mxcount), igatotsize,
     &     imaxsize(0:mxcount), igamaxsize
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      if (numestimate.lt.mxcount) then
         numestimate=numestimate+1
         itotsize(numestimate)=0
         imaxsize(numestimate)=0
         igmem_push_estimate=numestimate
      else
         write(iwr,*)'numestimate is ',numestimate,' mxcount is ',
     &               mxcount
         call caserr(
     &'igmem_push_estimate: no more memory estimate counters available')
      endif
      return
      end
c
c **********************************************************************
c * igmem_pop_estimate(ihandle): destroy a memory estimate counter
c *
c *   ihandle: the handle of the toplevel estimate counter
c *
c *   Return the maximum memory estimate counter
c *
c **********************************************************************

      integer function igmem_pop_estimate(ihandle)
      implicit none
      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
      integer numchunk
      integer numestimate
      integer ichunksize
      integer itotsize
      integer imaxsize
      integer igamatsize
      integer igatotsize
      integer igamaxsize
      common/gmemcount/
     &     numchunk, numestimate,
     &     ichunksize(mxheap),  igamatsize(mxheap),
     &     itotsize(0:mxcount), igatotsize,
     &     imaxsize(0:mxcount), igamaxsize
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      integer ipg_nodeid
      integer ihandle
c
      if (numestimate.eq.ihandle) then
         if (itotsize(numestimate).ne.0) then
            write(iwr,*)ipg_nodeid(),' WARNING: current memory ',
     &      ' estimate counter ',numestimate,' is not zero!!!'
            write(iwr,*)ipg_nodeid(),' WARNING:',
     &      ' estimate counter ',numestimate,' is ',
     &      itotsize(numestimate)
         endif
         igmem_pop_estimate=imaxsize(numestimate)
         itotsize(numestimate)=0
         imaxsize(numestimate)=0
         numestimate=numestimate-1
      else
         write(iwr,*)'numestimate is ',numestimate,' ihandle is ',
     &               ihandle
         call caserr(
     &        'igmem_pop_estimate: out of stack order destruction')
      endif
      return
      end
c
c **********************************************************************
c * igmem_get_estimate(ihandle): get the current value of the estimated
c *                              memory usage
c *
c *   ihandle: the handle of the estimate counter
c *
c *   Return the current memory estimate counter
c *
c **********************************************************************

      integer function igmem_get_estimate(ihandle)
      implicit none
      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
      integer numchunk
      integer numestimate
      integer ichunksize
      integer itotsize
      integer imaxsize
      integer igamatsize
      integer igatotsize
      integer igamaxsize
      common/gmemcount/
     &     numchunk, numestimate,
     &     ichunksize(mxheap),  igamatsize(mxheap),
     &     itotsize(0:mxcount), igatotsize,
     &     imaxsize(0:mxcount), igamaxsize
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      integer ipg_nodeid
      integer ihandle
c
      if (numestimate.ge.ihandle.and.ihandle.ge.0) then
         if (itotsize(numestimate).lt.0) then
            write(iwr,*)ipg_nodeid(),' WARNING: current memory ',
     &      ' estimate counter ',ihandle,' is ',itotsize(ihandle)
         endif
         igmem_get_estimate=itotsize(ihandle)
      else
         write(iwr,*)'numestimate is ',numestimate,' ihandle is ',
     &               ihandle
         call caserr(
     &        'igmem_get_estimate: ihandle out of range')
      endif
      return
      end
c
c **********************************************************************
c * igmem_max_estimate(ihandle): get the maximum value of the estimated
c *                              memory usage
c *
c *   ihandle: the handle of the estimate counter
c *
c *   Return the maximum memory estimate counter
c *
c **********************************************************************

      integer function igmem_max_estimate(ihandle)
      implicit none
      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
      integer numchunk
      integer numestimate
      integer ichunksize
      integer itotsize
      integer imaxsize
      integer igamatsize
      integer igatotsize
      integer igamaxsize
      common/gmemcount/
     &     numchunk, numestimate,
     &     ichunksize(mxheap),  igamatsize(mxheap),
     &     itotsize(0:mxcount), igatotsize,
     &     imaxsize(0:mxcount), igamaxsize
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      integer ipg_nodeid
      integer ihandle
c
      if (numestimate.ge.ihandle.and.ihandle.ge.0) then
         if (imaxsize(numestimate).lt.0) then
            write(iwr,*)ipg_nodeid(),' WARNING: current memory ',
     &      ' estimate counter ',ihandle,' is ',imaxsize(ihandle)
         endif
         igmem_max_estimate=imaxsize(ihandle)
      else
         write(iwr,*)'numestimate is ',numestimate,' ihandle is ',
     &               ihandle
         call caserr(
     &        'igmem_max_estimate: ihandle out of range')
      endif
      return
      end

c
c **********************************************************************
c * igmem_push_usage: create a new memory usage counter
c *
c *   Return the handle of a new memory counter
c *
c **********************************************************************

      integer function igmem_push_usage()
      implicit none
      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      if (numcount.lt.mxcount) then
         numcount=numcount+1
         igmem_totsize(numcount)=0
         igmem_maxsize(numcount)=0
         igmem_push_usage=numcount
      else
         write(iwr,*)'numcount is ',numcount,' mxcount is ',mxcount
         call caserr(
     &        'igmem_push_usage: no more memory counters available')
      endif
      return
      end
c
c **********************************************************************
c * igmem_pop_usage(ihandle): destroy a memory usage counter
c *
c *   ihandle: the handle of the toplevel counter
c *
c *   Return the maximum memory usage counter
c *
c **********************************************************************

      integer function igmem_pop_usage(ihandle)
      implicit none
      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      integer ipg_nodeid
      integer ihandle
c
      if (numcount.eq.ihandle) then
         if (igmem_totsize(numcount).ne.0) then
            write(iwr,*)ipg_nodeid(),' WARNING: current memory usage ',
     &      ' counter ',numcount,' is not zero!!!'
         endif
         igmem_pop_usage=igmem_maxsize(numcount)
         igmem_totsize(numcount)=0
         igmem_maxsize(numcount)=0
         numcount=numcount-1
      else
         write(iwr,*)'numcount is ',numcount,' ihandle is ',ihandle
         call caserr(
     &        'igmem_pop_usage: out of stack order destruction')
      endif
      return
      end
c
c **********************************************************************
c * igmem_get_usage(ihandle): get the current value of the memory usage
c *
c *   ihandle: the handle of the usage counter
c *
c *   Return the current memory usage counter
c *
c **********************************************************************

      integer function igmem_get_usage(ihandle)
      implicit none
      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      integer ipg_nodeid
      integer ihandle
c
      if (numcount.ge.ihandle.and.ihandle.ge.0) then
         if (igmem_totsize(numcount).lt.0) then
            write(iwr,*)ipg_nodeid(),' WARNING: current memory ',
     &      ' usage counter ',ihandle,' is ',igmem_totsize(ihandle)
         endif
         igmem_get_usage=igmem_totsize(ihandle)
      else
         write(iwr,*)'numcount is ',numcount,' ihandle is ',
     &               ihandle
         call caserr(
     &        'igmem_get_usage: ihandle out of range')
      endif
      return
      end
c
c **********************************************************************
c * igmem_max_usage(ihandle): get the maximum value of the memory usage
c *
c *   ihandle: the handle of the usage counter
c *
c *   Return the maximum memory usage counter
c *
c **********************************************************************

      integer function igmem_max_usage(ihandle)
      implicit none
      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      integer ipg_nodeid
      integer ihandle
c
      if (numcount.ge.ihandle.and.ihandle.ge.0) then
         if (igmem_maxsize(numcount).lt.0) then
            write(iwr,*)ipg_nodeid(),' WARNING: maximum memory ',
     &      ' usage counter ',ihandle,' is ',igmem_maxsize(ihandle)
         endif
         igmem_max_usage=igmem_maxsize(ihandle)
      else
         write(iwr,*)'numcount is ',numcount,' ihandle is ',
     &               ihandle
         call caserr(
     &        'igmem_max_usage: ihandle out of range')
      endif
      return
      end
c
c **********************************************************************
c * igmem_alloc_inf(size,filename,subname,varid,priority) :  allocate 
c *    a segment of memory
c *
c *    New style interface to the memory subsystem. 
c *    This routine delegates the real work to ogmem_alloc_kernel.
c *
c *        size     - size of memory segment                     (input)
c *        filename - name of the file from which 
c *                   the call was made                          (input)
c *        subname  - name of the subroutine from which 
c *                   the call is made                           (input)
c *        varid    - mnemonic identifier of the data structure
c *                   for which the memory is allocated          (input)
c *        priority - how important is it for a user to know 
c *                   the size of this memory segment            (input)
c *
c **********************************************************************

      integer function igmem_alloc_inf(size,filename,subname,varid,
     &                                 priority)
c
      implicit none
      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      integer size
      logical ogmem_alloc_kernel
      integer igmem_max_memory
      integer igmem_overhead
      integer i
      logical opg_root
      logical osuccess
      integer result
      integer maxmem
      character*(*) filename,subname,varid
      integer priority
c
      igmem_alloc_inf = 0
      osuccess = ogmem_alloc_kernel(size,varid,result)
      igmem_size(numheap) = size + igmem_overhead()
      do i=0,numcount
         igmem_totsize(i) = igmem_totsize(i) + igmem_size(numheap)
         igmem_maxsize(i) = max(igmem_maxsize(i),igmem_totsize(i))
      enddo
c
      if (opg_root().and.priority.le.igmem_print) then
         write(iwr,600)size,filename,subname,varid,result,
     &                 igmem_totsize(0)
 600     format(1x,'allocating ',i16,' words in file ',a16,
     &            ' routine ',a16,' for ',a16,' address ',2i12)
      endif
      if (opg_root().and.priority.ge.0) then
         call caserr('Illegal value for priority: are you using common/g
     &mempara allright?')
      endif
c
      if (.not.osuccess) then
         maxmem = igmem_max_memory()
         if (opg_root()) then
            call ma_summarize_allocated_blocks
            write(iwr,610)size,maxmem
         endif
 610     format(1x,'Memory allocation of ',i16,' words failed.',/,
     &          1x,'Only ',i16,' words available.')
         call caserr('memory ran out')
      endif
c
      igmem_priority(numheap)=priority
      zgmem_varid(numheap)   =varid
c
      igmem_alloc_inf = result
      return
      end


c **********************************************************************
c * igmem_alloc_all_inf(size,filename,subname,varid,priority) :  
c *    allocate all memory available
c *
c *    New style interface to the memory subsystem. 
c *    This routine delegates the real work to ogmem_alloc_all_kernel.
c *
c *        size     - size of memory segment                    (output)
c *        filename - name of the file from which
c *                   the call was made                          (input)
c *        subname  - name of the subroutine from which
c *                   the call is made                           (input)
c *        varid    - mnemonic identifier of the data structure
c *                   for which the memory is allocated          (input)
c *        priority - how important is it for a user to know
c *                   the size of this memory segment            (input)
c *
c **********************************************************************

      integer function igmem_alloc_all_inf(size,filename,subname,varid,
     &                                     priority)
c
      implicit none
      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      integer size,i
      logical ogmem_alloc_all_kernel
      integer igmem_max_memory
      integer igmem_overhead
      logical opg_root
      logical osuccess
      integer result
      integer maxmem
      character*(*) filename,subname,varid
      integer priority
c
      ogmem_alloced_all = .true.
      igmem_alloc_all_inf = 0
c
      if (opg_root().and.priority.le.igmem_print) then
         write(iwr,600)filename,subname,varid
 600     format(1x,'allocating    all remaining words in file ',a16,
     &            ' routine ',a16,' for ',a16)
      endif
c
      osuccess = ogmem_alloc_all_kernel(size,varid,result)
      igmem_size(numheap) = size + igmem_overhead()
      do i=0,numcount
         igmem_totsize(i) = igmem_totsize(i) + igmem_size(numheap)
c        igmem_maxsize(i) = max(igmem_maxsize(i),igmem_totsize(i))
      enddo
      if (opg_root().and.priority.ge.0) then
         call caserr('Illegal value for priority: are you using common/g
     &mempara allright?')
      endif
c
      if (.not.osuccess.or.size.le.0) then
         maxmem = igmem_max_memory()
         if (opg_root()) then
            call ma_summarize_allocated_blocks
            write(iwr,610)maxmem
         endif
 610     format(1x,'Memory allocation of all remaining words failed.',/,
     &          1x,'Only ',i16,' words available.')
         call caserr('memory ran out')
      endif
c
      igmem_priority(numheap)=priority
      zgmem_varid(numheap)   =varid
c
      igmem_alloc_all_inf = result
      return
      end


c **********************************************************************
c *
c *   gmem_free_set(ilow,ihigh) : free all memory chunks between and
c *                               including addresses ilow and ihigh.
c *
c *    This routine is meant to lift the burden to deallocate all
c *    memory chunks explicitly. This burden turned out to be a 
c *    significant obstacle to the take up of the dynamic memory 
c *    management.
c *
c *        ilow    - the lowest address of the memory chunks you
c *                  want to free                                (input)
c *        ihigh   - the highest address of the memory chunks you
c *                  want to free                                (input)
c *
c **********************************************************************

      subroutine gmem_free_set(ilow,ihigh)
      implicit none
      integer ilow, ihigh
      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
      integer ix, ixhigh, ixlow, iq_h
c
      ixhigh = numheap
      ixlow  = numheap
      do ix = numheap,1,-1
         if (iq_heap(ix).ge.ihigh) ixhigh = ix
         if (iq_heap(ix).ge.ilow)  ixlow  = ix
      enddo
      call gmem_free(ihigh)
      do ix = ixhigh-1, ixlow+1, -1
         iq_h =iq_heap(ix)
         call gmem_free(iq_h)
      enddo
      if (ihigh.ne.ilow) call gmem_free(ilow)
      end

c **********************************************************************
c *
c *  gmem_free_inf(iq,filename,subname,varid) : free memory at offset iq
c *
c *    New style interface to the memory subsystem. 
c *    This routine delegates the real work to gmem_free_kernel.
c *
c *        iq       - the offset of the memory chunk             (input)
c *        filename - name of the file from which
c *                   the call was made                          (input)
c *        subname  - name of the subroutine from which
c *                   the call is made                           (input)
c *        varid    - mnemonic identifier of the data structure
c *                   for which the memory is allocated          (input)
c *
c *    Although varid is stored when the allocate is done the user
c *    expected to enter it again. This ensures that varid is a string
c *    the user can search for in the source to find matching allocates
c *    and frees.
c *
c **********************************************************************

      subroutine gmem_free_inf(iq,filename,subname,varid)
      implicit none
c
      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      integer iq
      integer igmem_overhead
      logical opg_root
      integer ix, priority, len, size, i
      character*(*) filename,subname,varid
c
      priority=igmem_priority(numheap)
      size    =igmem_size(numheap)
      call strtrm(varid,len)
      len = min(len,16)
      if (varid(1:len).ne.zgmem_varid(numheap)(1:len)) then
         if (opg_root()) then
            write(iwr,*)'WARNING: freeing ',varid,' but was expecting ',
     &                  zgmem_varid(numheap)
         endif
      endif
c
      if (opg_root().and.priority.le.igmem_print) then
         write(iwr,600)size-igmem_overhead(),filename,subname,varid,
     &                 igmem_totsize(0)
 600     format(1x,'releasing  ',i16,' words in file ',a16,
     &            ' routine ',a16,' for ',a16,i12)
      endif
      do i=0,numcount
         igmem_totsize(i) = igmem_totsize(i) - size
      enddo
c
c     if (opg_root().and.ogmem_debug) then
c        do ix = numheap,1,-1
c           if (iq.eq.iq_heap(ix)) then
c              write(iwr,*)'free memory handle= ',itag_heap(ix)
c              goto 100
c           endif
c        enddo
c     endif
 100  continue
      call gmem_free_kernel(iq)
      end


c
c **********************************************************************
c * igmem_alloc(size) :  allocate a segment of memory
c *
c *    Old style interface to the memory subsystem. All calls to this
c *    routine should eventually be replaced by calls to 
c *    igmem_alloc_inf. This routine delegates the real work to 
c *    ogmem_alloc_kernel.
c *
c *                size - size of memory segment (input)
c *
c **********************************************************************

      integer function igmem_alloc(size)
c
      implicit none
      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      integer size, i
      logical ogmem_alloc_kernel
      integer igmem_max_memory
      integer igmem_overhead
      logical opg_root
      logical osuccess
      integer result
      integer maxmem
      character*30 tag
c
      igmem_alloc = 0
      write(tag,610)numheap+1
 610  format('gamess_',i3.3)
      osuccess = ogmem_alloc_kernel(size,tag,result)
      igmem_size(numheap) = size + igmem_overhead()
      do i=0,numcount
         igmem_totsize(i) = igmem_totsize(i) + igmem_size(numheap)
         igmem_maxsize(i) = max(igmem_maxsize(i),igmem_totsize(i))
      enddo
      zgmem_varid(numheap) = ' '
c
      if (.not.osuccess) then
         maxmem = igmem_max_memory()
         if (opg_root()) then
            call ma_summarize_allocated_blocks
            write(iwr,600)size,maxmem
         endif
 600     format(1x,'Memory allocation of ',i16,' words failed.',/,
     &          1x,'Only ',i16,' words available.')
         call caserr('memory ran out')
      endif
c
      if (opg_root().and.IGMEM_DEBUG.le.igmem_print) then
         write(iwr,620)tag,size,itag_heap(numheap),iq_heap(numheap),
     &                 igmem_totsize(0)
 620     format(1x,'allocate ',a10,' size=',i8,' handle= ',i8,
     &          ' gamess address=',i12,i12)
      endif
c
      igmem_alloc = result
      return
      end


c **********************************************************************
c * igmem_alloc_all(size) :  allocate all memory available
c *
c *    Old style interface to the memory subsystem. All calls to this
c *    routine should eventually be replaced by calls to 
c *    igmem_alloc_all_inf. This routine delegates the real work to 
c *    ogmem_alloc_all_kernel.
c *
c *                size - size of memory segment (output)
c *
c **********************************************************************

      integer function igmem_alloc_all(size)
c
      implicit none
      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      integer size, i
      logical ogmem_alloc_all_kernel
      integer igmem_max_memory
      integer igmem_overhead
      logical opg_root
      logical osuccess
      integer result
      integer maxmem
      integer lentag
      character*30 tag
      tag = 'all_remaining_memory'
c
      ogmem_alloced_all = .true.
      igmem_alloc_all = 0
      osuccess = ogmem_alloc_all_kernel(size,tag,result)
      igmem_size(numheap) = size + igmem_overhead()
      do i=0,numcount
         igmem_totsize(i) = igmem_totsize(i) + igmem_size(numheap)
      enddo
      zgmem_varid(numheap) = ' '
c
      if (.not.osuccess.or.size.le.0) then
         maxmem = igmem_max_memory()
         if (opg_root()) then
            call ma_summarize_allocated_blocks
            write(iwr,600)maxmem
         endif
 600     format(1x,'Memory allocation of all remaining words failed.',/,
     &          1x,'Only ',i16,' words available.')
         call caserr('memory ran out')
      endif
c
      call strtrm(tag,lentag)
      if (opg_root().and.IGMEM_DEBUG.le.igmem_print) then
         write(iwr,610)tag(1:lentag),size,itag_heap(numheap),
     &                 iq_heap(numheap)
 610     format(1x,'allocate ',a20,' size=',i14,' handle=',i8,
     &               ' gamess address=',i8)
      endif
c
      igmem_alloc_all = result
      return
      end


c **********************************************************************
c *
c *   gmem_free(iq) : free memory at offset iq
c *
c *    Old style interface to the memory subsystem. All calls to this
c *    routine should eventually be replaced by calls to
c *    gmem_free_inf. This routine delegates the real work to
c *    gmem_free_kernel.
c *
c **********************************************************************

      subroutine gmem_free(iq)
      implicit none
c
      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      integer iq
      logical opg_root
      integer ix
      integer i
      character *30 tag
c
      do i=0,numcount
         igmem_totsize(i) = igmem_totsize(i) - igmem_size(numheap)
      enddo
      if (opg_root().and.IGMEM_DEBUG.le.igmem_print) then
         do ix = numheap,1,-1
            if (iq.eq.iq_heap(ix)) then
               write(tag,610)ix
 610           format('gamess_',i3.3)
               write(iwr,600)tag,itag_heap(ix),iq
 600           format(1x,'free     ',a10,' handle= ',i8,
     &                ' gamess address=',i12,i12)
               goto 100
            endif
         enddo
      endif
 100  continue
      call gmem_free_kernel(iq)
      end


c **********************************************************************
c * ogmem_alloc_kernel(size,tag,offset) :  allocate a segment of memory 
c *
c *    This routine does the real work and never prints any 
c *    diagnostics. This routine should NEVER be called directly
c *    from an application.
c *    The return value is .true. if the allocation was successful
c *
c *                size - size of memory segment (input)
c *                tag  - name of memory segment (input)
c *                offset - offset to memory allocated (output)
c *
c **********************************************************************

      logical function ogmem_alloc_kernel(size,tag,offset)

      implicit none
      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
      real*8 qq
      integer ivoff, memhandle
      common/vcore/qq(2)
      common/vcoreoff/ivoff, memhandle
c
      integer nmaxly, ntotly, icurly, isecly
      logical oprintm
      common /segm/ nmaxly,ntotly,icurly,isecly(100),
     +              oprintm
c
       integer i10, need, loc10, loccm,  loadcm
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      logical osuccess, opg_root
      integer ipg_nodeid
      external ipg_nodeid
      character tag*(*)
      integer size, offset
      integer igmem_null

c
      ogmem_alloc_kernel = .false.
      osuccess = .false.
      numheap = numheap + 1

      if (numheap.gt.mxheap) then
         write(iwr,600)ipg_nodeid(),numheap,mxheap
 600     format(1x,'PROC ',i6,': current number of memory chunks ',i6,
     &             ' exceeds maximum of ',i6,/,
     &          1x,'Please change mxheap in gmemdata and recompile.')
         call caserr("numheap exceeds mxheap in memory allocation")
      endif

c
c implementation in terms of GAMESS-UK stack
c
c  based on model code...
c      length=n*n           
c      call cmem(loadcm  sets loadcm to top of core  (ntotly)
c      call setscm(i10)  get address of next scm bloc
c            k of core (ntotly + 1)
c      last=i10+length      !  
c      loc10=loccm()        !   ( returns ntotly )
c      need=loc10+length    !   
c      call setc(need)      !
c
c
c get current top of core
      call cmem(loadcm) 
c
c return address of the nect scm block of core
      call setscm(i10)
c
c loccm currently always returns ntotly
      loc10 = loccm()
c
c new top address, allowing 4 words for guards
      need = loc10 + size + 4

      if(oprintm) then
         write(iwr,*) 'loadcm ',loadcm,
     &        ' loc10 ',loc10,' i10 ',i10
      endif

      if(need.gt.nmaxly)then
         numheap = numheap - 1
         return
c        write(iwr,200) size, nmaxly-loc10-4
c200     format(1x,'memory allocation error: need',i10,' avail ',i10)
c        call caserr('memory ran out')
      endif

      call setc(need) 
c
c  determine stored offset (following MA convention)
c
      iq_heap(numheap) = i10 + 2
c
c  allocate a tag for tracing 
c
      igmem_count = igmem_count+1
      itag_heap(numheap) = igmem_count
c
c the gamess allocator as currently implemented uses a single
c physical stack - so the address allocated is 1 greater than 
c the "loadcm" value used to free the data. We use this
c assumption, so check for changes that may invalidate it
c
      if(i10 .ne. loadcm+1)call caserr('addressing error')
c
c define values of guards
c
      iqoff=0

      qq(ivoff+iq_heap(numheap)-iqoff-2) = size
      qq(ivoff+iq_heap(numheap)-iqoff-1) = 111.0d0
      qq(ivoff+iq_heap(numheap)-iqoff+size) = 222.0d0
      qq(ivoff+iq_heap(numheap)-iqoff+size+1) = size
c
      osuccess = .true.
c
c
      if (osuccess) then
         offset = iq_heap(numheap) 
         if (ogmem_nanify) then
            call gmem_nanify(qq(ivoff+offset),size)
         endif
         if (ogmem_debug) then
            call gmem_check_guards('gmem_alloc')
         endif
      else 
         numheap = numheap - 1
         offset = igmem_null()
      endif
      ogmem_alloc_kernel = osuccess
c
      return
      end

      function igmem_max_memory ()
c
      implicit none
      integer igmem_max_memory
      logical ogmem_alloc_all_kernel
      integer nmaxly, ioff
      logical osuccess
      logical osave_nanify
      character*30 tag
      data tag/"igmem_max_memory"/
      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
c
c     determine total memory currently available, and use
c     this in deriving data storage. Note that the allocation
c     and subsequent freeing may be inefficient.
c
      osave_nanify = ogmem_nanify
      ogmem_nanify = .false.
      osuccess = ogmem_alloc_all_kernel(nmaxly,tag,ioff)
      ogmem_nanify = osave_nanify
      if (osuccess) then
         igmem_max_memory = nmaxly
         call gmem_free_kernel(ioff)
      else
         igmem_max_memory = 0
      endif
c
      return
      end


c **********************************************************************
c *
c * ogmem_alloc_all_kernel(size,tag,offset) :  allocate all memory 
c *                                            available
c *
c *    This routine does the real work and never prints any 
c *    diagnostics. This routine should NEVER be called directly
c *    from an application.
c *    The return value is .true. if the allocation was successful
c *
c *                size - size of memory segment (output)
c *                tag  - name of memory segment (input)
c *                offset - offset to memory allocated (output)
c *
c **********************************************************************

      logical function ogmem_alloc_all_kernel(size,tag,offset)

      implicit none
      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
      integer size, offset
      character tag*(*)
      integer lentag
      logical osuccess
      logical opg_root
      integer igmem_null

c
c  GAMESS-UK stack implementation
c
      real*8 qq
      integer ivoff, memhandle
      common/vcore/qq(2)
      common/vcoreoff/ivoff, memhandle
c
      integer nmaxly, ntotly, icurly, isecly
      logical oprintm
      common /segm/ nmaxly,ntotly,icurly,isecly(100),
     +              oprintm
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      integer loadcm, i10, loc10, loccm, need
c
c  based on model code...
c      length=n*n           
c      call cmem(loadcm  sets loadcm to top of core  (ntotly)
c      call setscm(i10)  get address of next scm bloc
c            k of core (ntotly + 1)
c      last=i10+length      !  
c      loc10=loccm()        !   ( returns ntotly )
c      need=loc10+length    !   
c      call setc(need)      !
c
      osuccess = .false.
      ogmem_alloc_all_kernel = osuccess
      size = 0
      offset = igmem_null()

      numheap = numheap + 1

c get current top of core
      call cmem(loadcm)   

c return address of the nect scm block of core
      call setscm(i10)

c loccm currently always returns ntotly
      loc10 = loccm()
c
c work out available size allowing 4 words for guards
      size = nmaxly - i10 - 4 + 1
      need = nmaxly

      if(size.lt.0)then
c        write(iwr,200)
c200     format(1x,'memory allocation error: gmem_alloc_all but no',
c    +        ' core remains')
c        call caserr('memory ran out')
         numheap = numheap - 1
         return
      endif

      call setc(need) 
c
c  determine stored offset (following MA convention
c
c      write(6,*)'gamess addr',i10
      iq_heap(numheap) = i10 + 2
c
c  allocate a tag for tracing 
c
      igmem_count = igmem_count+1
      itag_heap(numheap) = igmem_count
c
c the gamess allocator as currently implemented uses a single
c physical stack - so the address allocated is 1 greater than 
c the "loadcm" value used to free the data. We use this
c assumption, so check for changes that may invalidate it
c
      if(i10 .ne. loadcm+1) call caserr('addressing error')
c
c define values of guards
c
      iqoff=0

      qq(ivoff+iq_heap(numheap)-iqoff-2) = size
      qq(ivoff+iq_heap(numheap)-iqoff-1) = 111.0d0
      qq(ivoff+iq_heap(numheap)-iqoff+size) = 222.0d0
      qq(ivoff+iq_heap(numheap)-iqoff+size+1) = size
c
      osuccess = .true.
c
c
      if (osuccess) then
         offset = iq_heap(numheap)
         if (ogmem_nanify) then
            call gmem_nanify(qq(ivoff+offset),size)
         endif
      else
         numheap = numheap - 1
      endif
      ogmem_alloc_all_kernel = osuccess
c
      return
      end

c **********************************************************************
c *
c *   gmem_free_kernel(iq) : free memory at offset iq
c *
c *    This routine does the real work and does not print any 
c *    information. The routine does trigger an error exit if it
c *    finds the memory corrupted. This routine should NEVER be called 
c *    directly from an application.
c *
c **********************************************************************

      subroutine gmem_free_kernel(iq)
      implicit none
c
      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      integer igmem_null
      real*8 qq
      integer ivoff, memhandle
      common/vcore/qq(2)
      common/vcoreoff/ivoff, memhandle
      integer loadcm, size
c
      logical ostat, opg_root
      integer iq, ix
      ostat = .true.
      if (ogmem_debug) then
         call gmem_check_guards('gmem_free')
      endif
      do  ix = numheap,1,-1
         if(iq.eq. iq_heap(ix))then
            if(ix .ne. numheap) then
               write(iwr,600)ix,numheap
               write(iwr,610)iq_heap(ix),iq_heap(numheap)
               call caserr('non-stack memory')
            endif
600         format('You want to free memory chunk ',i3,
     +             ' but the top of stack is memory chunk ',i3)
610         format('The indeces are ',i8,' and ',i8,' respectively.')
c
c implementation in terms of GAMESS-UK stack
c
c check guards
c
            size = nint(qq(ivoff+iq_heap(ix)-2))
            if(ogmem_debug)then
               call gmem_check_guards(' gmem_free ')

               write(6,*)size,qq(ivoff+iq_heap(ix)-1)
               write(6,*)nint(qq(ivoff+iq_heap(ix)+size+1)),
     &              qq(ivoff+iq_heap(ix)+size)
            endif

            if(size.lt.0.or.
     &           qq(ivoff+iq_heap(ix)-1).ne.111.0d0)then
               write(iwr,100)itag_heap(ix),iq_heap(ix)
 100           format(1x,'problem with lower guard, memory handle',i7,
     &                ' address',i7)
               call caserr('problem with memory guard')
            endif

            if(nint(qq(ivoff+iq_heap(ix)+size+1)) .ne. size .or.
     &           qq(ivoff+iq_heap(ix)+size).ne.222.0d0)then
               write(iwr,101)itag_heap(ix),iq_heap(ix)
 101           format(1x,'problem with upper guard, memory handle',i7,
     &                ' address',i7)
               call caserr('problem with memory guard')
            endif

            loadcm = iq_heap(ix) - 3
            call setc( loadcm )
c           if(opg_root() .and. ogmem_debug)
c    &           write(6,*)'free memory handle= ',itag_heap(ix)
            numheap = numheap - 1
            if (.not.ostat)call caserr('problem freeing memory ')
            iq = igmem_null()
            return
         endif
      enddo

      if(opg_root())write(6,*)'attempt to free address= ',iq

      call caserr('gmem_free - bad address')
      end
c
c  - allocate a reserved region
c
      integer function igmem_reserve_region(size)
      implicit none
      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
      real*8 qq
      integer ivoff, memhandle
      common/vcore/qq(2)
      common/vcoreoff/ivoff, memhandle
      integer size
      logical opg_root
      character tag*30

      nresreg = nresreg + 1
      if(nresreg.gt.mxresreg)call 
     &     caserr('too many reserved memory regions requested')

      numheap = numheap + 1
      write(tag,100)nresreg
 100  format('reserve_',i3.3)

      call caserr('gmem func unavailable')
c
      iadr(0,nresreg) = iq_heap(numheap) - 1
      isize(nresreg) = size
      nres(nresreg) = 0
      igmem_reserve_region = nresreg
      if(opg_root()  .and. ogmem_debug)then
         write(6,*)'reserve region tag=',tag(1:10),
     &        'base address =',iadr(0,nresreg) + 1,' size=',size
      endif
      return
      end

      subroutine gmem_free_region(ires)
      implicit none

      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      real*8 qq
      integer ivoff, memhandle
      common/vcore/qq(2)
      common/vcoreoff/ivoff, memhandle
      integer ires
      logical opg_root
      integer iq, ix
      if(ires.ne.nresreg)call 
     &     caserr('bad gmem_free_region call')

      iq = iadr(0,nresreg) + 1

      do  ix = 1,numheap
         if(iq.eq. iq_heap(ix) )then
            if(ix .ne. numheap)then
               write(iwr,600)ix,numheap
               write(iwr,610)iq_heap(ix),iq_heap(numheap)
               call caserr('non-stack memory')
            endif
600         format('You want to free memory chunk ',i3,
     +             ' but the top of stack is memory chunk ',i3)
610         format('The indeces are ',i8,' and ',i8,' respectively.')
      call caserr('gmem func unavailable')
            if(opg_root()  .and. ogmem_debug)
     &           write(6,*)'free reserver handle= ',itag_heap(ix)
            numheap = numheap - 1
cccc
            nresreg = nresreg - 1
cccc
            return
         endif
      enddo
      call caserr('gmem_free - bad address')
      end
c
      integer function igmem_alloc_reserved(ires,size)
      implicit none
      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
      integer size, nleft, ires
      logical opg_root
c
c check space..
c
      if(ires.lt.0.or.ires.gt.nresreg)call 
     &     caserr('bad reserved memory region index')

      nleft = isize(ires) - ( iadr(nres(ires),ires) - 
     &     iadr(0,ires) )
      if(nleft . lt . size)call 
     &     caserr('reserved memory region too small')
      nres(ires) = nres(ires) + 1
      if(nres(ires) .gt. mxressect)call
     &     caserr('too many sections from reserved region')
      iadr(nres(ires),ires) = iadr(nres(ires) - 1,ires) + size
      igmem_alloc_reserved = iadr(nres(ires) - 1,ires) + 1
      if(opg_root()  .and. ogmem_debug)then
        write(6,*)'mem from reserved region=',ires,'index=',nres(ires),
     &        'address=',igmem_alloc_reserved
      endif
      return
      end

      subroutine gmem_free_reserved(ires,i)
      implicit none
      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
      logical opg_root
      integer ires, i
      if(ires.lt.0.or.ires.gt.nresreg)call 
     &     caserr('bad reserved memory region index')
      if(opg_root()  .and. ogmem_debug)then
         write(6,*)'free mem from reserved region=',ires,
     &        'index=',nres(ires),
     &        'address=',iadr(nres(ires) - 1,ires) + 1
         if(iadr(nres(ires) - 1,ires) .ne. i - 1)call
     &        caserr('non-stack gmem_free_reserved call')
      endif
      nres(ires) = nres(ires) - 1
      return
      end
      subroutine ma_summarize_allocated_blocks
      return
      end

c
c **********************************************************************
c *
c *   gmem_nanify(core,size) : set contents of core to NaN.
c *
c *   This routine sets the contents of the a real memory segment to
c *   NaN. This is to be used only for debugging purposes so as to test
c *   for failures to initialise data properly.
c *
c **********************************************************************
c
      subroutine gmem_nanify(core,size)
      implicit none
c
c     Input/Output
c
      integer size
      real*8 core(size)
c
c     Local
c
      integer iunk(2),i
      real*8 junk
      equivalence(junk,iunk(1))
c
      iunk(1) = -1
      iunk(2) = -1
      do i = 1, size
         core(i) = junk
      enddo
      end

c ********************************************************************
c *
c * gmem_check_guards :  checks guards without freeing memory
c *
c *  
c *
c ********************************************************************

      subroutine gmem_check_guards(text)
      implicit none
      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
      real*8 qq
      integer ivoff, memhandle
      common/vcore/qq(2)
      common/vcoreoff/ivoff, memhandle
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      character*(*) text

      integer ix, size
      write(iwr,*)'checking memory guards numheap=',numheap,text
      do  ix = 1,numheap

         size = nint(qq(ivoff+iq_heap(ix)-2))
         if(ogmem_debug)then
cc            write(6,*)size,qq(ivoff+iq_heap(ix)-1)
cc            write(6,*)nint(qq(ivoff+iq_heap(ix)+size+1)),
cc     &           qq(ivoff+iq_heap(ix)+size)
         endif
 
         if(size.lt.0.or.
     &        qq(ivoff+iq_heap(ix)-1).ne.111.0d0)then
            write(iwr,100)itag_heap(ix),iq_heap(ix)
 100        format(1x,'problem with lower guard, memory handle',i7,
     &             ' address',i7)
         endif
 
         if(nint(qq(ivoff+iq_heap(ix)+size+1)) .ne. size .or.
     &        qq(ivoff+iq_heap(ix)+size).ne.222.0d0)then
            write(iwr,101)itag_heap(ix),iq_heap(ix)
 101        format(1x,'problem with upper guard, memory handle',i7,
     &             ' address',i7)
         endif
         
      enddo
      return
      end

      subroutine gmem_summarize(filename,subrname,tag,ipriority)
      implicit none
      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      integer ipriority
      character *(*) filename,subrname,tag
      integer i
      logical opg_root
      external opg_root
c
      if (opg_root().and.ipriority.le.igmem_print) then
         write(iwr,600)filename,subrname,tag,numheap,igmem_totsize(0)
 600     format(//,
     &          1x,'gmem_summarize',/,
     &          1x,'==============',/,
     &          1x,'in file=',a16,' routine=',a16,' tag=',a16,/,
     &          1x,'on stack are ',i8,' chunks totalling ',i8,' words:')
         do i=numheap,1,-1
            write(iwr,630)i,iq_heap(i),itag_heap(i),igmem_size(i),
     &                    zgmem_varid(i)
 630        format(1x,'chunk=',i3,' address=',i8,' handle=',i8,
     &                ' size=',i8,' name=',a16)
         enddo

         call gmem_highwater(ipriority)
      endif
      end

      subroutine gmem_highwater(ipriority)
      implicit none
      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      integer ipriority
      integer i
      logical opg_root
      external opg_root
c
      if (opg_root().and.ipriority.le.igmem_print) then
         write(iwr,600)
 600     format(/,
     &          1x,'Memory high water mark',/,
     &          1x,'======================',/)
         if (ogmem_alloced_all) then
            write(iwr,605)
         endif
         write(iwr,610)igmem_maxsize(0)
         if (igmem_totsize(0).ne.0) then
            write(iwr,630)igmem_totsize(0)
         endif
 605     format(1x,'*** accurate heap memory high water mark not ',
     &             'available',
     &        /,1x,'*** allocation of ALL memory occurred at ',
     &             'some point',/,
     &          1x,'*** the high water mark below excludes such ',
     &             'allocations',/)
 610     format(1x,'heap memory high water mark = ',i12,' words')
 620     format(1x,'GA   memory high water mark = ',i12,
     &             ' words overall')
 630     format(1x,'*** heap not cleared ',i12,' words remain',
     &             ' allocated')
 640     format(1x,'*** GAs  not cleared ',i12,' words remain',
     &             ' allocated')

      endif
      end
c
c     Set the nanify flag
c
      subroutine gmem_set_nanify(flag)
      implicit none
      logical flag
      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
      ogmem_nanify = flag
      end
c
c controls memory debug print
c
      subroutine gmem_set_debug(flag)
      implicit none
      logical flag
      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
      ogmem_debug = .false.
      if (flag) then
        igmem_print = IGMEM_DEBUG
      else 
        igmem_print = IGMEM_NORMAL
      endif
      end
c
      subroutine gmem_set_quiet(flag)
      implicit none
      logical flag
      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
      ogmem_debug = .false.
      if (flag) then
        igmem_print = IGMEM_QUIET
      else 
        igmem_print = IGMEM_NORMAL
      endif
      end
c
      subroutine memtest(q,nmaxly)
      implicit none
      real*8 q(*)
      real*8 small
      integer i, iq, iq2, nmaxly, itest

      integer igmem_alloc_inf
      external igmem_alloc_inf

      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
      real*8 qq
      integer ivoff, memhandle
      common/vcore/qq(2)
      common/vcoreoff/ivoff, memhandle
      character*10 fnm
      character*7  snm

      data fnm/'parallel.m'/
      data snm/'memtest'/

      small = 0.0001d0
c
c test passed q
c
c      write(6,*)'test q'
c      do i=1, nmaxly
c         q(i) = 99.0d0
c      enddo
c      do i=1, nmaxly
c         if(dabs(q(i) - 99.0d0) .gt. small)write(6,*)q(i)
c      enddo
c

      call gmem_set_debug(.true.)

      itest = (nmaxly - 100) / 2
      write(6,*)'allocate',itest
      iq = igmem_alloc_inf(itest,fnm,snm,'iq',IGMEM_DEBUG)
      write(6,*)'allocate',itest
      iq2 = igmem_alloc_inf(itest,fnm,snm,'iq2',IGMEM_DEBUG)

      write(6,*)'test q(iq)', iq, iq2
      do i=1, itest
         q(iq+i-1) = 99.0d0
      enddo
      do i=1, itest
         if(dabs(q(iq+i-1) - 99.0d0) .gt. small)write(6,*)q(i)
      enddo
          
      call gmem_free_inf(iq2,fnm,snm,'iq2')
      call gmem_free_inf(iq,fnm,snm,'iq')

      itest = nmaxly - 100

      write(6,*)'test qq(ivoff+iq)',itest

      iq = igmem_alloc_inf(itest,fnm,snm,'iq',IGMEM_DEBUG)

      do i=1, itest
         qq(ivoff+iq+i-1) = 99.0d0
      enddo
      do i=1, itest
         if(dabs(qq(ivoff+iq+i-1) - 99.0d0) .gt. small)
     &        write(6,*)qq(ivoff+iq+i-1)
      enddo

      call gmem_free_inf(iq,fnm,snm,'iq')

      end
      subroutine debugp(s)
      implicit none
      character s*(*)
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

      integer ipg_nodeid
      if(odebugp)then
         write(6,*)'PAR:',ipg_nodeid(),s
      endif
      end
      subroutine ver_parallel(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/parallel.m,v $
     +     "/
      data revision /"$Revision: 6287 $"/
      data date /"$Date: 2013-05-27 14:27:44 +0200 (Mon, 27 May 2013) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
