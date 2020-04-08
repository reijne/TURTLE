














































c
c     GDF  June 21 2007  glue_uk_us.m
c     'glue' to connect GAMESS(US) code to GAMESS-UK program
c
c     NOTE: need CPP for 'mafdecls.fh' (see ddi_create, below)
c
      subroutine cp_uk_us
      implicit real*8 (a-h,o-z)
c
c     subroutine below copies data from G-UK common blocks into 
c     G(US) ones, some G(US) blocks must be renamed (add '_f')
c
c     -----------------  G-UK common blocks  ------------------
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
c  global array pointers
      integer g_oooo, g_vooo, g_vvoo, g_vovo, g_vvvo, g_ovaa, g_u
     &,  il_oooo, ih_oooo, jl_oooo, jh_oooo
     &,  il_vooo, ih_vooo, jl_vooo, jh_vooo
     &,  il_vvoo, ih_vvoo, jl_vvoo, jh_vvoo
     &,  il_vovo, ih_vovo, jl_vovo, jh_vovo
     &,  il_vvvo, ih_vvvo, jl_vvvo, jh_vvvo
     &,  il_ovaa, ih_ovaa, jl_ovaa, jh_ovaa
     &,  il_u, ih_u, jl_u, jh_u
      common/gms_global/g_oooo, g_vooo, g_vvoo
     &,  g_vovo, g_vvvo, g_ovaa, g_u
     &,  il_oooo, ih_oooo, jl_oooo, jh_oooo
     &,  il_vooo, ih_vooo, jl_vooo, jh_vooo
     &,  il_vvoo, ih_vvoo, jl_vvoo, jh_vvoo
     &,  il_vovo, ih_vovo, jl_vovo, jh_vovo
     &,  il_vvvo, ih_vvvo, jl_vvvo, jh_vvvo
     &,  il_ovaa, ih_ovaa, jl_ovaa, jh_ovaa
     &,  il_u, ih_u, jl_u, jh_u
      integer g_hess, il_hess,ih_hess,
     &  jl_hess,jh_hess
      common/extraglobal/g_hess,il_hess,ih_hess,
     &  jl_hess,jh_hess
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
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
c
      real*8 tr, trx, try, trz, prmoms, aprev
      integer indmx, jaxis, igroup, igrp80
      logical oprmom
      common /molsym/ tr(3,3),trx,try,trz,indmx,jaxis,igroup,igrp80,
     +                prmoms(3),aprev(maxat3,3),oprmom
c
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
c     -----------------  G(US) common blocks  ------------------
c
cNotes:
c     In this subroutine some common block variables are renamed 
c     to avoid clashes, after this they get their original names.
c     Simplest to give G(US) common blocks their own lengths for now.
c
      parameter (mxatm=2000, mxrt=100, mxao=8192, mxsh=5000)
      parameter (mxgtot=20000)
      common /infoa_f/nat_f,ich_f,mul_f,num_f,nqmt,ne_f,na_f,nb_f,
     *                zan(mxatm),c_f(3,mxatm),ian(mxatm)
      common /nshel_f/ex_f(mxgtot),cs_f(mxgtot),cp_f(mxgtot),
     *cd_f(mxgtot),
     *                cf_f(mxgtot),cg_f(mxgtot),ch(mxgtot),
     *ci(mxgtot),
     *                kstart_f(mxsh),katom_f(mxsh),ktype_f(mxsh),
     *kng_f(mxsh),
     *                kloc_f(mxsh),kmin_f(mxsh),kmax_f(mxsh),
     *nshell_f
c
c     these G(US) common blocks are independent of G-UK data
c     but best to set some parameters at the start
c
      character*8 runtyp,exetyp,scftyp,cityp
      common /runopt_f/ runtyp,exetyp,nevals,nglevl,nhlevl
      common /wfnopt_f/ scftyp,cityp,dftype,cctyp,mplevl,mpctyp
      common /iofile_f/ ir,iw,ip,ijko,ijkt,idaf,nav,ioda(950)
      integer         d_oooo,d_vooo,d_vvoo,d_vovo
      common /trfdms/ d_oooo,d_vooo,d_vvoo,d_vovo
      logical                                     goparr,dskwrk,maswrk
      common /par   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
      common /fmoinf/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody
      common /ijpair/ ia(mxao)
      common /enrgys/ enucr,eelct,etot_f,sz,szz,ecore,escf,eerd,e1,e2,
     *                ven,vee,epot,ekin,estate(mxrt),statn,edft(2)
      logical                                           abel
      common /symmol/ group,complex,igroup_f,naxis,ilabmo,abel
      common /machin_f/ nwdvar,maxfm,maxsm,limfm,limsm
      common /symtry_f/ mapshl(mxsh,48),mapctr(mxatm,48),
     *                t(432),invt_f(48),nt_f
      common /output_f/ nprint_f,itol,icut,normf,normp,nopk
c
c     simulate g(us) exe
c
c     for /runopt/
c
      character*8 energy,run
      data energy,run/'energy  ','run     '/
c
c     for /wfnopt/
c
      character*8 rmc,ormas
      data rmc,ormas /'mcscf   ','ormas   '/
c
c     /runopt/ 
c
      runtyp = energy
      exetyp = run
c
c     /wfnopt/ 
c
      scftyp = rmc
      cityp  = ormas
      mplevl = 0
c
c     /iofile/
c
      ir   = 5
      iw   = 6
c
c     /par   /
c
      call ddi_nproc(nproc,me)
      master = 0
      ibtyp  = 1             ! nxtval
      iptim  = 0             ! =1 for // timing output
      goparr = .true.
      dskwrk = .false.       ! flag for master,slave i/o
      maswrk = me.eq.master
c 
c     /infoa/ stuff
c
      nat_f = nat
      ich_f = ich
      mul_f = mul
      num_f = num
      nqmt  = num    ! = /harmon/...newbas0
      ne_f  = ne
      na_f  = na
      nb_f  = nb
      do i = 1, nat
        c_f(1,i) = c(1,i)
        c_f(2,i) = c(2,i)
        c_f(3,i) = c(3,i)
      end do
c
c     /nshel/ stuff
c
      nshell_f = nshell
      do i = 1, nshell
        kstart_f(i) = kstart(i)
         katom_f(i) =  katom(i)
         ktype_f(i) =  ktype(i)
           kng_f(i) =    kng(i)
          kloc_f(i) =   kloc(i)
          kmin_f(i) =   kmin(i)
          kmax_f(i) =   kmax(i)
      end do
      nprim = kstart(nshell) + kng(nshell) -1
      do i = 1, nprim
        ex_f(i) = ex(i)
        cs_f(i) = cs(i)
        cp_f(i) = cp(i)
        cd_f(i) = cd(i)
        cf_f(i) = cf(i)
        cg_f(i) = cg(i)
      end do
c
c     initialize the g(us)-style memory management
c
      call setfm
c
c     when available, use existing GA handles
c
      d_oooo = g_oooo
      d_vooo = g_vooo
      d_vvoo = g_vvoo
      d_vovo = g_vovo
c
c     finally, ensure the pair-index is initialized 
c
      ii = 0
      do i = 1, mxao
        ia(i) = ii
        ii = ii + i
      end do
c
c     /fmoinf/
c
      nfg = 0
c
c     /enrgys/, get enucr from /scfopt/
c
      enucr =  enucf(nat,czan,c)
c
c     /symmol/, from /molsym/
c
      igroup_f = igroup
      naxis    = jaxis
      abel     = .true.
c
c     /symtry/
c
      nt_f = nt
c
c     /machin/, nwdvar saves space with i4
c
      nwdvar = 1
c
c     /output/
c
      icut = 9
      nprint_f = 7
c
      end
c
cgdf  parallel interface from gamess(us) ddi to gamess-uk (pg_)
c     append to parallel.m 
c   ? parallel.m = only file in gamess-uk to be cpp pre-processed
c     needed for 'mafdecls.fh' stuff (e.g. see 6 lines down)
c
c     assume init,end done elsewhere
c
      subroutine ddi_create(handle,nrows,ncols)
      implicit   none




!
!     $Id: mafdecls.fh,v 1.2 2007-08-19 18:54:00 mrdj Exp $
!

!
!     Public header file for a portable dynamic memory allocator.
!
!     This file may be included by internal and external FORTRAN files.
!
































!
!     The guard ends here instead of at the end of the file because we only
!     need the cpp constants (stuff above) defined once per FORTRAN file,
!     but need the declarations (stuff below) to be defined each time this
!     file is included in a FORTRAN file.
!



!
!     constants
!

!     type declarations for datatype constants
      integer	MT_BYTE		! byte
      integer	MT_INT		! integer
      integer	MT_LOG		! logical
      integer	MT_REAL		! real
      integer	MT_DBL		! double precision
      integer	MT_SCPL		! single precision complex
      integer	MT_DCPL		! double precision complex

      integer	MT_F_FIRST	! first type
      integer	MT_F_LAST	! last type

!     parameter declarations for datatype constants
      parameter	(MT_BYTE = (1000 + 9))
      parameter	(MT_INT = (1000 + 10))
      parameter	(MT_LOG = (1000 + 11))
      parameter	(MT_REAL = (1000 + 12))
      parameter	(MT_DBL = (1000 + 13))
      parameter	(MT_SCPL = (1000 + 14))
      parameter	(MT_DCPL = (1000 + 15))

      parameter	(MT_F_FIRST = MT_BYTE)
      parameter	(MT_F_LAST = MT_DCPL)

!
!     function types
!


      logical MA_alloc_get
      logical MA_allocate_heap
      logical MA_chop_stack
      logical MA_free_heap
      logical MA_free_heap_piece
      logical MA_get_index
      logical MA_get_next_memhandle
      logical MA_get_numalign
      logical MA_init
      logical MA_initialized
      logical MA_init_memhandle_iterator
      integer MA_inquire_avail
      integer MA_inquire_heap
      integer MA_inquire_heap_check_stack
      integer MA_inquire_heap_no_partition
      integer MA_inquire_stack
      integer MA_inquire_stack_check_heap
      integer MA_inquire_stack_no_partition
      logical MA_pop_stack
!     subroutine MA_print_stats
      logical MA_push_get
      logical MA_push_stack
      logical MA_set_auto_verify
      logical MA_set_error_print
      logical MA_set_hard_fail
      logical MA_set_numalign
      integer MA_sizeof
      integer MA_sizeof_overhead
!     subroutine MA_summarize_allocated_blocks
!     subroutine MA_trace
      logical MA_verify_allocator_stuff

      external MA_alloc_get
      external MA_allocate_heap
      external MA_chop_stack
      external MA_free_heap
      external MA_free_heap_piece
      external MA_get_index
      external MA_get_next_memhandle
      external MA_get_numalign
      external MA_init
      external MA_initialized
      external MA_init_memhandle_iterator
      external MA_inquire_avail
      external MA_inquire_heap
      external MA_inquire_heap_check_stack
      external MA_inquire_heap_no_partition
      external MA_inquire_stack
      external MA_inquire_stack_check_heap
      external MA_inquire_stack_no_partition
      external MA_pop_stack
!     external MA_print_stats
      external MA_push_get
      external MA_push_stack
      external MA_set_auto_verify
      external MA_set_error_print
      external MA_set_hard_fail
      external MA_set_numalign
      external MA_sizeof
      external MA_sizeof_overhead
!     external MA_summarize_allocated_blocks
!     external MA_trace
      external MA_verify_allocator_stuff


!
!     variables
!





!     common blocks



      common /mbc_byte/		byte_mb(2)
      character*1		byte_mb



      common /mbc_int/		int_mb(2)
      integer			int_mb



      common /mbc_log/		log_mb(2)
      logical			log_mb



      common /mbc_real/		real_mb(2)
      real			real_mb



      common /mbc_dbl/		dbl_mb(2)
      double precision		dbl_mb



      common /mbc_scpl/		scpl_mb(2)
      complex			scpl_mb



      common /mbc_dcpl/		dcpl_mb(2)
      double complex		dcpl_mb


      integer    handle,nrows,ncols
      integer dditype
      logical    success, pg_create
      character*8  name
c
c     'dditype' refers to integer, double - see mafdecls.fh
c
      success=.FALSE.
      end
c
      subroutine ddi_destroy(handle)
      implicit   none
      integer    handle
      logical    success, pg_destroy
c
      success =.FALSE.
      end
c
      subroutine ddi_distrib(handle,me,ilo,ihi,jlo,jhi)
      implicit   none
      integer    handle,me,ilo,ihi,jlo,jhi
c
      end
c
      subroutine ddi_get(handle,ilo,ihi,jlo,jhi,buff)
      implicit   none
      integer    handle,ilo,ihi,jlo,jhi
      real*8     buff(*)
c
      end
c
      subroutine ddi_put(handle,ilo,ihi,jlo,jhi,buff)
      implicit   none
      integer    handle,ilo,ihi,jlo,jhi
      real*8     buff(*)
c
      end
c
      subroutine ddi_acc(handle,ilo,ihi,jlo,jhi,buff)
      implicit   none
      integer    handle,ilo,ihi,jlo,jhi
      real*8     buff(*)
c
      end
c
      subroutine ddi_gsumf(tag,buff,len)
      implicit   none
      integer    tag,len
      real*8     buff(*)
c
      end
c
      subroutine ddi_gsumi(tag,buff,len)
      implicit   none
      integer    tag,len,buff(*)
c
      end
c
      subroutine ddi_bcast(tag,dditype,buff,len,from)
      implicit   none
      integer    tag,buff(*),len,from
      character*1 dditype
c
cgdf  02/08/06  need support for i8
c
      end
c
      subroutine ddi_sync(tag)
      implicit   none
      integer    tag
c
      end
c
      subroutine ddi_nproc(ddi_np,ddi_me)
      implicit   none
      integer    ddi_np,ddi_me,ipg_nnodes,ipg_nodeid
c
      ddi_np=1
      ddi_me=1
      end 
c
      subroutine ddi_dlbnext(dlb_counter)
      implicit   none
      integer    dlb_counter,ipg_dlbtask
c
c     note that the DDI counter starts from 0 
c     while the pg_ one starts from 1
c     easy to modify the quantum chemistry to fit
c     in so leave this as is
c
      end 
c
      subroutine ddi_dlbreset
      implicit   none
c
      end
c
c     memory management interface from gamess(us) to MA tool
c     MA returns word addresses into the q() array
c
      subroutine valfm(loadfm)
      implicit none
      integer  loadfm, igmem_alloc_inf, addr
c
      integer mxaddr,addrss,addrcnt
      parameter (mxaddr=100)
      common /memusuk/ addrss(mxaddr),addrcnt
c
c     In gamess-uk all routines pass the work array, q
c     igmem_alloc_inf provides addresses into q().
c     Give address of next free word then free it.
c
      addr = igmem_alloc_inf(1,'fullnr','valfm','work',-1)
c
c     let gamess(us) do its address arithmetic based on 'addr'
c     then allocate the full amount again in getfm.
c     (the call to gmem_free_inf below corrupts addr)
c     subtract 1 word because valfm returns the top of the 
c     current stack
c
      loadfm = addr - 1
      call gmem_free_inf(addr,'ormas','valfm','work')
      end
c
      subroutine getfm(need)
      implicit none
      integer  loadfm, igmem_alloc_inf, need
      integer  igmem_alloc
c
      integer mxaddr,addrss,addrcnt
      parameter (mxaddr=100)
      common /memusuk/ addrss(mxaddr),addrcnt
c
c     have to assume that valfm,getfm are always called together
c     but not necessarily followed by retfm
c
      addrcnt = addrcnt + 1
      loadfm = igmem_alloc_inf(need,'ormas','getfm','work',-1)
c
c     store this address in a new location for later use 
c     when freeing memory in gmem_free
c
      addrss(addrcnt) = loadfm
      end
c
      subroutine retfm(need)
      implicit none
      integer  need, loadfm
c
      integer mxaddr,addrss,addrcnt
      parameter (mxaddr=100)
      common /memusuk/ addrss(mxaddr),addrcnt
c
c     'need' is not needed here, 
c     gmem_free frees the -address-, like malloc
c
      loadfm = addrss(addrcnt)
      call gmem_free_inf(loadfm,'ormas','retfm','work')
c
c     done with this address
c
      addrcnt = addrcnt - 1
      end
c
      subroutine setfm
      implicit none
c
      integer mxaddr,addrss,addrcnt
      parameter (mxaddr=100)
      common /memusuk/ addrss(mxaddr),addrcnt
c
c     initialize address counter to zero, this 
c     must be done before any g(us)-style memory 
c     management is done
c
      addrcnt = 0
      end
c
c     miscellaneous stubs
c     feel free to add what you like here
c
      subroutine abrt
      call caserr('abrt called')
      end
c
      subroutine flshbf
      call flush(6)
      end
c
      subroutine tsecnd(cpu)
      implicit none
      real*8 cpu, cpulft
c
c     returns the elapsed cpu time
c
      cpu = cpulft(1)
      end

      subroutine glue()
CMR   intialize us-uk glue
CMR   currently a stub
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
CMR     write(iwr,1)
1       format('Glue present')
        return
      end
