c
c     GDF  June 21 2007  glue_uk_us.m
c     'glue' to connect GAMESS(US) code to GAMESS-UK program
c
c     NOTE: need CPP for 'mafdecls.fh' (see ddi_create, below)
c
      subroutine cp_uk_us
      implicit REAL (a-h,o-z)
c
c     subroutine below copies data from G-UK common blocks into 
c     G(US) ones, some G(US) blocks must be renamed (add '_f')
c
c     -----------------  G-UK common blocks  ------------------
c
INCLUDE(common/sizes)
INCLUDE(common/global)
INCLUDE(common/scfopt)
INCLUDE(common/symtry)
INCLUDE(common/molsym)
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
#include "../g/ma/mafdecls.fh" 

      integer    handle,nrows,ncols
      MA_INTEGER dditype
      logical    success, pg_create
      character*8  name
c
c     'dditype' refers to integer, double - see mafdecls.fh
c
_IF(ga)
      dditype = MT_DBL
      success = pg_create(dditype,nrows,ncols,name,1,1,handle)
_ELSE
      success=.FALSE.
_ENDIF
      end
c
      subroutine ddi_destroy(handle)
      implicit   none
      integer    handle
      logical    success, pg_destroy
c
_IF(ga)
      success = pg_destroy(handle)
_ELSE
      success =.FALSE.
_ENDIF
      end
c
      subroutine ddi_distrib(handle,me,ilo,ihi,jlo,jhi)
      implicit   none
      integer    handle,me,ilo,ihi,jlo,jhi
c
_IF(ga)
      call pg_distribution(handle,me,ilo,ihi,jlo,jhi)
_ENDIF
      end
c
      subroutine ddi_get(handle,ilo,ihi,jlo,jhi,buff)
      implicit   none
      integer    handle,ilo,ihi,jlo,jhi
      real*8     buff(*)
c
_IF(ga)
      call pg_get(handle,ilo,ihi,jlo,jhi,buff,1)
_ENDIF
      end
c
      subroutine ddi_put(handle,ilo,ihi,jlo,jhi,buff)
      implicit   none
      integer    handle,ilo,ihi,jlo,jhi
      real*8     buff(*)
c
_IF(ga)
      call pg_put(handle,ilo,ihi,jlo,jhi,buff,1)
_ENDIF
      end
c
      subroutine ddi_acc(handle,ilo,ihi,jlo,jhi,buff)
      implicit   none
      integer    handle,ilo,ihi,jlo,jhi
      real*8     buff(*)
c
_IF(ga)
      call pg_acc(handle,ilo,ihi,jlo,jhi,buff,1)
_ENDIF
      end
c
      subroutine ddi_gsumf(tag,buff,len)
      implicit   none
      integer    tag,len
      real*8     buff(*)
c
_IF(ga)
      call pg_dgop(tag,buff,len,'+')
_ENDIF
      end
c
      subroutine ddi_gsumi(tag,buff,len)
      implicit   none
      integer    tag,len,buff(*)
c
_IF(ga)
      call pg_dgop(tag,buff,len,'+')
_ENDIF
      end
c
      subroutine ddi_bcast(tag,dditype,buff,len,from)
      implicit   none
      integer    tag,buff(*),len,from
      character*1 dditype
c
cgdf  02/08/06  need support for i8
c
_IF(ga)
      if (dditype.eq.'f'.or.dditype.eq.'F') then
        call pg_brdcst(tag,buff,len*8,from)
      else
        call pg_brdcst(tag,buff,len*4,from)
      end if
_ENDIF
      end
c
      subroutine ddi_sync(tag)
      implicit   none
      integer    tag
c
_IF(ga)
      call pg_synch(tag)
_ENDIF
      end
c
      subroutine ddi_nproc(ddi_np,ddi_me)
      implicit   none
      integer    ddi_np,ddi_me,ipg_nnodes,ipg_nodeid
c
_IF(ga)
      ddi_np = ipg_nnodes()
      ddi_me = ipg_nodeid()
_ELSE
      ddi_np=1
      ddi_me=1
_ENDIF
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
_IF(ga)
      dlb_counter = ipg_dlbtask()
_ENDIF
      end 
c
      subroutine ddi_dlbreset
      implicit   none
c
_IF(ga)
      call pg_dlbreset
_ENDIF
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
      REAL cpu, cpulft
c
c     returns the elapsed cpu time
c
      cpu = cpulft(1)
      end

      subroutine glue()
CMR   intialize us-uk glue
CMR   currently a stub
INCLUDE(common/iofile)
CMR     write(iwr,1)
1       format('Glue present')
        return
      end
