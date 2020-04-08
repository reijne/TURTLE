*
*     second-order scf module for nonorthogonal orbitals
*     author Zahid Rashid and Joop H. van Lenthe
*     Version 1.0, Revision 0 DATE: 2013/02/06
*     Revision 1.1, fixed parallel version, DATE: 2014/03/15 
*
************************************************************************
      subroutine vbqcscf(v,maxvec,ci,q)
************************************************************************
*
      implicit real*8 (a-h,o-z), integer (i-n)
*
*     control routine for orbital optimisation
*     using newton-raphson scheme or super-ci
*     adapted from 'subroutine vbscf'
*     v  : orbitals
*     ci : ci vector 
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
c...   turtle parameters
c
      integer melec,ncorq,maxex,maxact,mbonds,maxato,
     &        maxopa,ncases,maxcon,mmsocc,mxorbvb,maxspvb
c
      parameter (melec=300)
      parameter (ncorq=2000000)
      parameter (maxex=20000,maxact=150,mbonds=150,
     &           maxato=750,maxopa=1500)
      parameter (ncases=35)
c      ... maxcon in vbcrestr ...
c92      parameter (maxcon = 100)
      parameter (maxcon = 3000)
c      ... msocc max. # singly occuppieds ; was always enough (2009) ...
      parameter (mmsocc = 30)
c..   maxspvb  maximum # spinfunctions per configuration 
c..   (42 for 10 singlys and the ultimate answer)
      parameter (maxspvb=42) 
c
c     mxorbvb may be as large as allowed by n_8_16 (32767)
c
      parameter (mxorbvb=1000)
*
*     to signal which scf procedure to be used, e.g., super-ci or Newton-Raphson
*     default is super-ci
*
      real*8 scalhs, sci2nr, delta, crmax
      integer nits, nitsci, nege, nneg
      logical oqcscf, ocong, ospci, oinve, ograd, ocorr, oauto,
     +        oiter, oqcof
      common / quadratic / nitsci,  nits,  nege, nneg, scalhs, sci2nr,
     +                     delta, crmax, oqcscf, ocong, ospci, oinve, 
     +                     ograd, ocorr, oauto, oiter, oqcof

c
c .. common for  diis in vb /vbdiist/
c
      real*8 evbpre
      integer iwdiis,ntdiis,madiis,midiis,kvb7_vnew
      logical ovbdiis
      common /vbdiist/ evbpre,iwdiis,ntdiis,madiis,midiis,
     1                 kvb7_vnew,ovbdiis
c     8note* tractlt noo included in vbdens.m mkd1vb
      integer iprint,nbasis,ncol,nsa,ncurt,lenact,lenbas,nactiv
      integer num3,iblkq,ionsec,nbsort,isecdu,n2int,setintpos
      logical incsort,oprs
      integer nnsort,isort,ksort,maxsort
      real*8 scri
c
      common /tractlt/ scri,iprint,nbasis,ncol,nsa,ncurt,lenact,lenbas,
     &                 num3,iblkq,ionsec,nbsort,isecdu,n2int,nactiv,
     &                 incsort,nnsort,isort,ksort,maxsort,oprs
c
c...  incsort : indicates if incore sorting is on
c...  nnsort : # sortbuffers
c...  isort : sortbuffer done
c...  ksort : adres of sortbuffers in qq (igmem_alloc)
c...  maxsort : maximum sort buffer used
      common /bypass/ index4,ihmat,idavid
      logical ohmat
      integer nstruc,iprinv,ihfile,ihbloc,ipbloc,idumhs
      common/hsinfo/ohmat(3),nstruc,iprinv,ihfile,ihbloc,ipbloc,idumhs
c...   common meant for all kind of (fixed) vb criteria
c...   critor .. criterium to determine orthogonality
c...   crilow .. criterium  for orthogonalization (e.g.jacoby in lowdin)
c...   cridep .. criterium for dependency (on c**2)
c...   symcri .. criterium for symmetry equivalence
c...   remcri .. criterium for including excitation on basis of 1-e ints
c...   crihyb .. max hybrid contam. allowed in matrix and (*100) vector
c...   cribri .. min. size of bcoeff to avoid annihilation
c...   crinor .. min. norm of a vector
c...   cribub .. a meaningful difference (in a bubble sort)
c...   criign .. size of matrix elements that really may be ignored
c...   cripop .. criterium for pipek/mizek localization
c...   normally critor = 1.0d-14, crilow = critor/10, cridep = critor*10
c...            symcri = critor*100 remcri = 1.0d-10 crihyb = 1.0d-14
c...            cribri = 1.0d-8 crinor = critor*100  cribub = 1.0d-10
c...            criign = 1.0d-44 cripop 1.0d-10
      real*8 critor,crilow,cridep,symcri,remcri,crihyb,cribri,
     1     crinor,cribub,criign,cripop
      common/vbcri/ critor,crilow,cridep,symcri,remcri,crihyb,cribri,
     1              crinor,cribub,criign,cripop
*
*     for dynamic memory allocation debugging purpose, contains 
*     IGMEM_ vars.
*
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
*
      common /brill/ nex(maxact) ,iex (maxex,maxact),
     &               nexw(maxact),iexw(maxex,maxact),
     &                            ieqsig(maxex,maxact),
     &               iequi(maxact),nequi,iredund(mxorbvb),nredund,
     &               inteq(maxact+1,maxact),ninteq(maxact),ninter,
     &               nscf,nullity,ortscf,schmi,reduce,autscf,
     &               exact(maxact),nbrill,idequi,
     &               ieqs(maxact),nexeq(maxact),iexeq(maxact*maxact)
      logical ortscf,schmi,reduce,autscf,exact
      integer nequi,nredund,ninter,nscf,nullity,nbrill,nex,nexw,
     &        iex,iexw,ieqsig,iequi,iredund,inteq,ninteq,idequi,
     &        ieqs,nexeq,iexeq
*
      real*8  shiscf,criscf,evb,ebi,swave,percmax,perccur,brmax,
     &      dampvb,fixb0
      integer optcri,iprins,nitscf,maxscf,isecv,nitdiis,maxdiis,ntscf,
     &        ipcmt,maxscf_s
      logical perall,scfconv,p0core,unitary
      common /scftvb/dampvb,fixb0,shiscf,criscf,evb,ebi,swave,percmax,
     &               perccur,brmax,
     &               optcri,iprins,nitscf,maxscf,isecv,nitdiis,
     &               maxdiis,ntscf,ipcmt,maxscf_s,perall,scfconv,
     &               p0core,unitary
c
c     core :  frozen core energy for the whole job
c     core_i :  frozen core energy  in makeps0 
c
      real*8 core,pot_nuc,core_i
      integer nfil,nofil,iblvb,lblvb,iblhs
      common /ffile/ core,core_i,pot_nuc,nfil,
     1               nofil(20),iblvb(20),lblvb(20),iblhs
c
c...  frozen core definition for TURTLE
c
      real*8 dxyz
      integer ncore,mapcie,mapiee,ncore_i,nsplice
      parameter (nsplice=512)
      common/splice/dxyz(3),ncore,ncore_i,
     1              mapcie(nsplice),mapiee(nsplice)
*
      common /davcom/ eshift,thresj,thresd,maxcyc,maxdav,maxsel,iprind,
     +                alter,firsts
      logical alter,firsts
*
C...
c...  atonam for names of atoms (for mulliken population analysis,
c...  and automatic atomic (strictly localised) scf-procedure)
c...  atom is used her to denote fragments
c...  nopa contains the # orbitals per atom, iopa contains them
c...  mullik indicates if the analysis is requested, hybry indicates
c...  wether the scf-procedure just should used atom-centered ao's
c...  per atomic scf orbital. iacat contains the active (scf) orbital
c...  per atom, nacat contains the number of the per atom.
c...  atomao : fragment numbers / ao (may be more)
c...           then stored as frag1, maxato*frag2, maxato**2*frag3 etc.
c     ifrzat : hybrid handling
c              0 : normal, 1 frozen, 2 super
c...
      character*8 atoms
      common /atonam/ atoms(maxato)
      integer atomao,nopa,iopa,natom,iacat,nacat,ifrzat
      common /infato/ nopa(maxato),iopa(maxopa,maxato),natom,
     &                iacat(maxact,maxato),nacat(maxato),
     &                atomao(mxorbvb),ifrzat(maxato),
     &                mullik,hybry,clean,hguess,guesao
      logical mullik,hybry,clean,hguess,guesao
c     common for equivalences in VB
c     eqvmo/eqvao list of equivalent MO's and AO's in that order.
      integer eqvmo,eqvao,neqmo,neqao,meqmo,meqao
      logical equiv,equivir
      parameter (meqmo=10,meqao=100)
      common/vbequiv/eqvmo(meqmo),eqvao(meqao,meqmo),neqmo,neqao,
     1               equiv,equivir
*
      common/vbpert/crivarpert(4),shiftp,fpdoc,fpuoc,
     &              nvarp,ifirsp,idelp,npert,
     &              ntype,itype(mxorbvb),
     &              minfock,maxfock,minpert,maxpert,minvar,maxvarp,
     &              minfockd,maxfockd,
     &              fckp,ovbper,fockdiagonal
      INTEGER nvarp,ifirsp,idelp,npert,ntype,itype
      INTEGER minfock,maxfock,minpert,maxpert,minvar,maxvarp
      INTEGER minfockd,maxfockd
      logical ovbper,fckp,fockdiagonal
      real*8 crivarpert,shiftp,fpdoc,fpuoc
c description:
c   ntype -> 1 - specify optimisation per brillioun structure
c            2 - specify optimisation per orbital
c            3 - specify optimisation per kind (type excitation)
c   itype ->
c
c          ntype = 1
c            itype(1..maxex) - optimisation type for each (brillouin) structure
c
c          ntype = 2
c            itype(1..mxorbvb) - optimisation type for each orbital
c            
c          ntype = 3 
c            itype(1) - optimisation type for doubly to variably occupied
c            itype(2) - optimisation type for doubly to unoccupied
c            itype(3) - optimisation type for variably to variably occupied
c            itype(4) - optimisation type for variably to unoccupied
c            
c          ntype = 4
c            itype(1..mxorbvb) - optimisation type for each orbital (to be
c                               specified by clsbril
c
c   contents of itype array:
c            itype(x)=1 for Super CI optimalisation
c            itype(x)=2 for perturbation theory optimalisation
c            itype(x)=4 perturbation theory  determines choice between Super and Pert
c
c   crivarpert ~: critaria for itype=4 case per class if applicable

*
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
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
*
c
c     twice contains info on pauli-principle related restrictions
c     ("doubly","redorb","iflip"), and super, which indicates wether
c     internal excitations are allowed for (not so if super-mcscf)
c     alos ifreez/nfreez orbitals to be frozen
      integer inforb,ninfor,idoubly,ndoubly,isingly,nsingly,kaart,
     &        ioror,igno_sh,igno_sel
      logical super,super_hybrid,ofrzmo,super_cas,super_act
      common /twice/ inforb(maxact*5,maxact),ninfor(maxact),
     &               idoubly(maxact),ndoubly,isingly(maxact),nsingly,
     &               kaart(mxorbvb+maxact),ioror(maxact,maxact),super,
     &               super_hybrid,igno_sh(maxato),ofrzmo(maxact),
     &               super_cas,super_act,igno_sel
c...
c...  common to signal first point of (e.g.) geometry optimisation
      common/first_vb/ofirst,opr_first
      logical ofirst,opr_first
      character*10 hsmatt
      common /hsmattype/hsmatt
*
c
c...  Atoms in Valence Bond (AiVB) directive, common block, work in progress
c...  Marcin Zielinski 08-09/2007, 04-11/2008, 03-11/2009, 04-06,09-10/2010
c
c...  _multname = allowed multiplicity names, in input
c...  _multrdin = multiplicity names readed from input
c...  _atstname = allowed atomic term symbols, in input
c...  _atstrdin = atomic term symbols readed from input
c...  _atnmrdin = atom's names readed from input
c...  _atindex  = atom indexes
c...            ,1: indicates atom number from geometry input
c...            ,2: x index from _statedef
c...            ,3: y index from _statedef
c...            ,4: index from aivb_multname
c...  _atcntr   = atom counter, gives number of atoms affected by input definition per structure
c...  _allatst  = allowed atomic states (terms) for P,S block atoms - (x,y) indexing
c...            x: 1-6 numerates one 'p' to six 'p' electron atoms, 7 and 8 indexes are one and two 's' electron atoms, 
c...               9 represents 's' with no electrons
c...            y: lists possible atomic states
c...  _atstcom  = boundary limits for _confdef/_statedef tables
c...            (x,1)    - x=1..9 - number of atomic states per given number of 'p'/'s' electrons
c...            (1..9,y) - number of projections of each of the atomic state

c...  _confdef  = all possible configuration determinants for given 'x' 'p' electron atoms
c...            x: 1-5 numerates one 'p' to five 'p' electron atoms, 6 and 7 numerates one and two 's' electron atoms
c...            y: indexing all possible determinants involved in the given 's'/'p' situation
c...            z: stores given determinant

c...  _statedef = definitions of atomic states
c...            x: 1-6 indexes are one 'p' to six 'p' electron atoms, 7 and 8 indexes are one and two 's' electron atoms, 
c...               9 represents 's' with no electrons
c...            y: index of the given atomic state withing given 'x'
c...            z: index of the given projection of the atomic state within give 'x'
c...            u: stores the atomic state definition in terms: (i,j1,k1,j2,k2,j3,k3,.....), where
c...               'i' - how many configurations of x,y,z orb. (from _confdef) are involved for given atomic state
c...               'j1' - constant in front of the given configuration involved, 'k1' - which determinant is that (from _confdef)
c...               if 'i' > 1 then 'j2' and 'k2' is the next pair (next determinant from _confdef)

c...  _statelab = string describing given atomic state
c...            x: 1-5 numerates one 'p' to five 'p' electron atoms, 6 and 7 numerates one and two 's' electron atoms
c...            y: numerates given atomic state withing given 'x'
c...            z: index = 1 - stores general name of the atomic state (e.g. trip P)
c...               index 2....n - store the atomic state's descriptions in terms of L, M and cos/sin (e.g. |L=x M=y; cos>)


c...  _detstruc = number of determinants per structure. for now 1000 structures is top number
c...              keep in mind for future !!
c...  _atord    = stores given 'elconfvaria' configuration in <singly><doubly> manner, for 'spinef'
c...  _elgrps   = 
c...         ,1 - stores the group (number of the atom), to which given singly occ. orbital belongs to
c...         ,2 - stores the given spin pattern for the given singly occ. orbitals
c...  _conford  = keeps the number of possible atomic state variations per atom
c...              say: _conford[1] = 3, _conford[2] = 1, _conford[3] = 1 means: 
c...              3 variations of atomic state for atom 1, 1 varia. of atomic state for atom 2 and 3
c...  _confgrps = gives configuration of atomic state variations in group-like scheme, based od _conford
c...              say: for previous example 3,1,1 _confgrps would look like: 1 1 1 2 3
c...  _coreorb  = core orbitals per atom (closed shells only)
c...  _confvaria= number of variations of given configuration (singly occupied orbs. only) - 1000 max. for now
c...  _confpr   = number of configurations per readnr - helpful for 'aimo diff' directive
c...  _reduce   = remove the number of redundant configurations with structure coef. less then 10E-08
c...  _aimo     = indicate that we're dealing with aimo directive, helpful in various places outside vbcrestr.m
c...  _diff     = indicate that we want to use different atomic orbitals for different atomic states of atoms, in various confs.
c...  _trconf   = truncate configurations, contains numbers of configurations per atom per readnr, which we want to keep
c...            x: readnr
c...            y: 1= how many configurations after truncation
c...               2= number indexes of the configurations to keep
c...            z: (x,2,i...max) explanation above
c...  _debugpr  = debug prints flags:
c...         1  = vbcrest.m, cresti subroutine, aivb part printouts -
c...         2  = vbaivb.m, read_aivb_input subroutine -
c...         3  = vbaivb.m, generate_aivb_drv subroutine, coeff and vectors/matrices -
c...         4  = vbaivb.m, generate_aivb_drv subroutine, generate kconf and kconfdetdef -
c...         5  = vbaivb.m, generate_aivb_drv subroutine, generate kconfdef -
c...         6  = vbaivb.m, generate_aivb_drv subroutine, generate elconfvaria -
c...         7  = vbaivb.m, generate_aivb_drv subroutine, generate determinants - print coeffs, vectors and matrices -
c...         8  = vbaivb.m, gen_aivbconf_drv subroutine, generate couplings between different atomic states -
c...         9  = vbcrestr.m, spinef subroutine, permute pairs and assign pairs printing
c...        10  = vbcrestr.m, spinef subroutine, define structures printing
c...        11  = vbcrestr.m, spinef subroutine, debugging prints + the rest
c...        12  = vbaivb.m, memory allocation/deallocation -
c...        13  = vbaivb.m, gen_aivb_descr subroutine -
c...        14  = vbaivb.m, compperm subroutine -
c...  _detnr    = variable used in 'spinef' subroutine = detnrold
c...  _cpats    = configuration nr per given atomic state definition 
c...              indicates whether we're dealing with the first configuration from the definition or further
c...              helps with checks
c... _lmc       = last main configuration, for generating configuration's determinants
c... _nconfs    = how many configurations make up given at. state (1 if >1, 0 if =1)
c... _confnr    = stores the number of the given configuration in the given readnr (regardless whether it's stored or not)
c
      character(len=1  ) aivb_atstname(4),aivb_atstrdin(maxact,maxact)
      character(len=2  ) aivb_zchr(maxact)
      character(len=20 ) aivb_atnmrdin(maxact,maxact)
      character(len=4  ) aivb_multname(6),aivb_multrdin(maxact,maxact)
      character(len=6  ) aivb_allatst(9,10)
      character(len=100) aivb_confdescr(100,5)
      character(len=20 ) aivb_statelab(9,4,6)
      character(len=70 ) aivb_conflab(9,4,5)
c
      common/aivbchar/aivb_atstrdin,aivb_multrdin,aivb_atnmrdin,
     &      aivb_allatst,aivb_multname,aivb_atstname,
     &      aivb_confdescr,aivb_statelab,aivb_conflab,
     &      aivb_zchr
c
      integer :: aimo_atindex(maxact,4),aimo_atcntr(maxact),
     &        aimo_atstcom(10,5),aimo_nnd(maxact),
     &        aimo_detnr,aimo_cpats,
     &        aimo_atord(mmsocc),aimo_detstruc(maxact,2,maxcon),
     &        aimo_elgrps(mmsocc,2),aimo_coreorb(maxact),
     &        aimo_confgrps(maxact*4),aimo_conford(maxact),
     &        aimo_confvaria(maxact),aimo_confpr(maxact),
     &        aimo_trconf(maxact,2,mmsocc),aimo_debugpr(20),
     &        aimo_confdef(9,13,10),aimo_statedef(9,4,5,11),
     &        aimo_confnr(maxact),aivb_ndocc,aimo_lterm(mmsocc),
     &        aimo_nconfs,aivb_diffadd,aivb_mo2perm(maxact),aivb_nmo2p,
     &        aivb_nlterm,aimo_lmc,aivb_inpspp(100)
c
      common/aivbint/aimo_atindex,aimo_atcntr,aimo_atstcom,
     &      aimo_detstruc,aimo_atord,aimo_elgrps,
     &      aimo_coreorb,aimo_conford,aimo_confgrps,aimo_confvaria,
     &      aimo_confpr,aimo_trconf,aimo_debugpr,aimo_lterm,
     &      aimo_confdef,aimo_statedef,aimo_confnr,aimo_nnd,
     &      aimo_detnr,aimo_cpats,aimo_lmc,aimo_nconfs,aivb_ndocc,
     &      aivb_diffadd,aivb_mo2perm,aivb_nmo2p,aivb_nlterm,aivb_inpspp
c
      logical :: aivb_reduce,aivb_set,aivb_diff,aivb_moperm,aivb_optorb
      common /aivblog/ aivb_reduce,aivb_set,aivb_diff,aivb_moperm,
     &       aivb_optorb

      real*8 :: aivb_threshold, aivb_nms(maxact)
      common /aivbreal/ aivb_threshold, aivb_nms
      common /vbtimess/ tvirtu,t4indx0,t4indx,tmatre0,tmatre,tdavid,
     1                  tlagr,tgtr,tfock, 
     2                  n4indx0,n4indx,nmatre0,nmatre,nfock
      real*8  tvirtu,t4indx0,t4indx,tmatre0,tmatre,tdavid,tlagr,tgtr,
     3      tfock
      integer n4indx0,n4indx,nmatre0,nmatre,nfock

*
      common /txema/ struc,ivec
      logical struc
*
      common/checkblock/ icheckblock
      common/gjs   /nirr,mult(8,8),isymao(maxorb),isymmo(maxorb)
     +              ,irr,iss,nstart(8),nfin(8),nmc
*
      character*10 zstr
      logical oroot
      integer idumpr(5,2)
*
      dimension v(*),q(*),ci(*)
*
*
*     ncol2 is # of vectors before redundant ones are deleted
*
      save ncol2
*
      icheckblock = 0
      call init_vbdiis
      swave = 0.0d0                    
      evbpre = 0.0d0
      crmax = 50.0d0
*
      ntscf = ntscf + 1
*
      tvirtu = 0.0d0
      t4indx0 = 0.0d0
      n4indx0 = 0
      t4indx = 0.0d0
      n4indx = 0
      tmatre0 = 0.0d0
      nmatre0 = 0
      tmatre = 0.0d0
      nmatre = 0
      tdavid = 0.0d0
      tlagr  = 0.0d0
      tgtr   = 0.0d0
      minfock = 999999999
      maxfock = 0
      minfockd = 999999999
      maxfockd = 0
      minpert = 999999999
      maxpert = 0
      minvar = 999999999
      maxvarp = 0
*
      kmaxmem = igmem_max_memory()
      kmemscr = kmaxmem/10
      kmemleft = kmaxmem - kmemscr
      ks   = igmem_alloc_inf(nbasis*(nbasis+1)/2,'vbqcscf.m','vbqcscf', 
     +                       'ks',IGMEM_DEBUG)
      kscr = igmem_alloc_inf(kmemscr,'vbqcscf.m','vbqcscf','kscr',
     +                           IGMEM_DEBUG)

      call get1e(q(ks),dummy,'s',q(kscr))

      if (ofirst.and.iprins.gt.50) then
         if (iprins.gt.10000) then
            write(iwr,10)
10          format(/,' the ao metric :')
            call tripri(q(ks),nbasis)
         end if
         call sigorb(q(ks),nbasis)
      end if
*
*     nitscf  = iteration # in the scf process
*     scfconv = true if the optimisation has converged
*     nactiv  = # active, i.e. occupied, mo's
*     ncol    = # orbitals (core,active,virtual), without redundant
*     ncol2   = # orbitals (core,active,virtual), with redundant
*                 (redorb will reduce ncol)
*     nsa     = # orbitals in 4-index (active+virtual)
*
      call getqvb(v,nbasis,ncol,isecv,'nopr')
      call normt(v,nbasis,ncol,q(ks),q(kscr))
*
*     make sure equivalence restrictions are obeyed 
*
      if (equiv) call primequiv(v,nbasis,1,'occ','at start')
*
      call putqvb(v,nbasis,ncol)
      call vbfirst(v,maxvec)
      call gmem_free_inf(kscr,'vbqcscf.m','vbqcscf','kscr')
      call gmem_free_inf(ks,'vbqcscf.m','vbqcscf','ks')
*
      nitscf = 0
      nits   = 0
      scfconv = .false.
      npert = 0
      evb = 0.0d0
*
      nactiv = nsa
*
*     (this nsa was the # active orbitals in the active statement)
*
      nsa = ncol - ncore
      nvirt = nsa - nactiv
      ncol2 = ncol
*
*     start of iterative loop
*

42    nitscf = nitscf + 1

*
*     perform allowed orthogonalisations (if ortscf)
*

      if (ortscf) 
     +   call vbcanon(v,nbasis,ncore+nactiv,100,'occ',q)
*
*     make sure equivalence restrictions are obeyed 
*
      if (equiv) call primequiv(v,nbasis,1,'occ','in scf iteration')
*
*     calculate psi(0)
*       
*
      evbpre = evb
*
20    call makepsi0(v,ci)
*
      if (ovbdiis) then
         if (evbpre.lt.evb) then

            write(iwr,*) ' '
            write(iwr,30)
30          format(100('-'))
            write(iwr,40) evbpre,evb
40          format(4x,'==> DIIS failed  previous E = ',f22.14,
     +             ' last E = ',f22.15,' <==')
            write(iwr,30)

            kvb7_vnew = kscra7vb('kvb7_vnew',nbasis*nactiv,'r','r')
            call rdedx(v(ncore*nbasis+1),nbasis*nactiv,kvb7_vnew,num8)
            ovbdiis = .false.
            call init_vbdiis
            go to 20
         end if
      end if
*
*     delta can also be used as convergence controler
*
      delta = evb - evbpre
*
*
*     perform orthogonalisations as defined in forbex
*
       print*, 'restror '
         print*, ' oqcscf ',oqcscf
      call restror(v(ncore*nbasis+1),nsa,nbasis,q)
       print*, 'after restror '
      a1 = cpulft(1)
      call start_time_period(TP_VB_VIRT)
*
*
*     construct virtual space (see virtual)
*
*     use special techniques to get a properly shaped virtual
*     space (see "virtual")
*
      kocc = 1
      maxvirt = maxvec - (nscf+ncore)
      call virtual(v(kocc),v(kocc + (nscf+ncore)*nbasis),
     +               nbasis,ncore,ncore+nscf,nvirt,maxvirt)
      if (super_hybrid) nsa = ncol - ncore
*
*     remove redundant virtuals
*
*     the vectors can be changed by redorb ! (redundant
*     virtuals are thrown away, iex is adapted) nredund
*     sits in /brill/. ncol is saved in ncol2
*

      if (reduce) then
        if (equiv.and.nitscf.le.0) 
     +         write(iwr,*) ' Remco : redorb is fatal'
        if (nactiv.ne.nsingly+ndoubly) 
     +     call caserr('nactiv ne nsingly+ndoubly')
*
*
        ksao = igmem_alloc_inf(nbasis*(nbasis+1)/2,'vbqcscf.m',
     +                         'vbqcscf','ksao',IGMEM_DEBUG)
        kvcopy = igmem_alloc_inf(nbasis*(nbasis+nactiv+ncore),
     +                           'vbqcscf.m','vbqcscf','kvcopy',
     +                           IGMEM_DEBUG)
        kiocvi = igmem_alloc_inf((nactiv+ncore-1)/nipw()+1,'vbqcscf.m',
     +                           'vbqcscf','kiocvi',IGMEM_DEBUG)
        kiset = igmem_alloc_inf(nbasis,'vbqcscf.m','vbqcscf','kiset',
     +                          IGMEM_DEBUG)
        kiact = igmem_alloc_inf(nactiv,'vbqcscf.m','vbqcscf','kiact',
     +                          IGMEM_DEBUG)
        kidoc = igmem_alloc_inf(nactiv,'vbqcscf.m','vbqcscf','kidoc',
     +                          IGMEM_DEBUG)
        kisoc = igmem_alloc_inf(nactiv,'vbqcscf.m','vbqcscf','kisoc',
     +                          IGMEM_DEBUG)
        kivir = igmem_alloc_inf(nbasis,'vbqcscf.m','vbqcscf','kivir',
     +                           IGMEM_DEBUG)
        nam = nbasis + nscf
        khmoao = igmem_alloc_inf(nam*(nam+1)/2,'vbqcscf.m','vbqcscf',
     +                           'khmoao',IGMEM_DEBUG)
        kiex2 = igmem_alloc_inf(nactiv*(nbasis+1),'vbqcscf.m','vbqcscf',
     +                          'kiex2',IGMEM_DEBUG)
*
*
        call get1e(q(ksao),head,'s',q(ksao))
*
        call redorb(q(ksao),v(ncore*nbasis+1),q(kvcopy),v,q(kiocvi),
     +               q(kiset),q(kiact),q(kidoc),q(kisoc),q(kivir),
     +               q(khmoao),q(kiex2),nbasis,q)

        ncol = ncol2 - nredund
        nsa = ncol - ncore
        nvirt = nsa - nactiv
        call gmem_free_inf(kiex2,  'vbqcscf.m', 'vbqcscf', 'kiex2' )
        call gmem_free_inf(khmoao, 'vbqcscf.m', 'vbqcscf', 'khmoao')
        call gmem_free_inf(kivir,  'vbqcscf.m', 'vbqcscf', 'kivir' )
        call gmem_free_inf(kisoc,  'vbqcscf.m', 'vbqcscf', 'kisoc' )
        call gmem_free_inf(kidoc,  'vbqcscf.m', 'vbqcscf', 'kidoc' )
        call gmem_free_inf(kiact,  'vbqcscf.m', 'vbqcscf', 'kiact' )
        call gmem_free_inf(kiset,  'vbqcscf.m', 'vbqcscf', 'kiset' )
        call gmem_free_inf(kiocvi, 'vbqcscf.m', 'vbqcscf', 'kiocvi')
        call gmem_free_inf(kvcopy, 'vbqcscf.m', 'vbqcscf', 'kvcopy')
        call gmem_free_inf(ksao,   'vbqcscf.m', 'vbqcscf', 'ksao'  )        
*


      end if

*
*     dump vectors (including virtuals) for orthopt
*
      nn = ncore + nactiv + nvirt
      kvb7_vv = kscra7vb('kvb7_vv',nn*nbasis,'r','w')
      call wrt3(v,nn*nbasis,kvb7_vv,num8)
*
*     print the orbitals / excitation patterns the first iter.
*
      if (nitscf.eq.1) then
         nn = nactiv
         if (iprins.gt.1) nn = ncore + nactiv
         if (iprins.gt.10) nn = ncore + nactiv + nvirt
      else
         nn = 0
         if (iprins.ge.5) nn = nactiv
         if (iprins.ge.50) nn = ncore + nactiv + nvirt
      end if
      if (nn.eq.nactiv) then
         kvv = ncore*nbasis + 1
         zstr = 'excl. core'
         kk = 0
      else
         kvv = 1
         zstr = 'incl. core'
         kk = ncore
      end if
*
      if (nn.ge.nactiv) then
         write(iwr,609) ncore,'frozen','fzc',(i,i=1,ncore)
         write(iwr,609) ndoubly,'doubly','doc',
     +                 (idoubly(i)+ncore,i=1,ndoubly)
         write(iwr,609) nsingly,'singly','voc',
     +                 (isingly(i)+ncore,i=1,nsingly)
         write(iwr,610) nvirt,nbasis
609      format(1x,'+++++',i6,1x,a6,' occupieds (',a3,'):',(t36,15i4))
610      format(1x,'+++++',i6,1x,' virtuals / ',i6,' basis functions')
         if (nsingly+ndoubly.ne.nactiv) 
     +       call caserr('nsingly-ndoubly - nactiv  clash')
*
         write(iwr,601) zstr
601      format(1x,/,'         ============================= ',
     +             /,'          real excitations ',a10,
     +             /,'         ============================= ')
         do i=1,nscf
           if (nex(i).gt.0) write(iwr,602) (iex(j,i)+kk,j=1,nex(i)+1)
602        format('       mo',i3,' ==>',15i4,/,(16x,15i4))
         end do
*
*    check
*
         kk = 0
         do i=1,nscf
            do j=1,nex(i)+1
               do k=j+1,nex(i)+1
                  if (iex(j,i).eq.iex(k,i)) then
                     kk = kk + 1
                  end if
               end do
            end do
         end do
         if (kk.gt.0) then
            write(iwr,'(a,i5)') ' *** found double arbitals # ',kk
            call caserr(' program error - double excitations')
         end if
*
         if (super_hybrid.and.nn.le.ncore+nactiv.and.igno_sel.ne.99)then
            write(iwr,606)
606         format(1x,/,'         ============================= ',
     +                /,'          super_hybrid - virt vs. ao ',
     +                /,'         ============================= ')
            kkkv = (ncore+nactiv)*nbasis + 1
            kkk = 1
            do i=1,nvirt
               call abmax(v(kkkv),nbasis,1,'sq',1,am,im)
               if (am.ne.1.0d0) call caserr('super_hybrid funny')
               idumpr(kkk,1) = i+nactiv+kk
               idumpr(kkk,2) = im
               if (kkk.eq.5) then
                  write(iwr,607) (idumpr(k,1),idumpr(k,2),k=1,kkk)
                  kkk = 0
               end if
               kkk = kkk + 1
               kkkv = kkkv + nbasis
            end do
            kkk = kkk - 1
            if(kkk.ne.0)write(iwr,607) (idumpr(k,1),idumpr(k,2),k=1,kkk)
607         format(7x,10(:' (',i4,' - ',i4,')  '))
         end if
*
         if (nn.eq.ncore+nactiv) write(iwr,603) ncore
         if (nn.lt.ncore+nactiv) write(iwr,604)
         if (nn.gt.ncore+nactiv) write(iwr,605) ncore,nvirt
603      format(1x,/,'             ===================== ',
     +             /,'                  real vectors',
     +             /'                 incl',i3,' core',
     +             /,'             ===================== ')
604      format(1x,/,'             ===================== ',
     +             /,'                  real vectors',
     +             /,'             ===================== ')
605      format(1x,/,'             ===================== ',
     +             /,'                  real vectors',
     +             /'            incl',i3,' core and',i3,' virt',
     +             /,'             ===================== ')
         call prvc(v(kvv),nn,nbasis,v,'o','l')
      end if
*
*     check hybrids
*


      if (hybry) call clvecvb(v(ncore*nbasis+1),nactiv+nvirt,nbasis,
     +                        'check','null','in vbqcscf')
*
*     save orbitals to disk as these are the vbqcscf orbitals
*
      call putqvb(v,nbasis,nactiv+ncore)
*
      call flushbuffer
*
*
*     the stuff before had to be done for consistent orthog etc.

      if (scfconv) goto 100

*
      a1 = cpulft(1) - a1
      tvirtu = tvirtu + a1
      call end_time_period(TP_VB_VIRT)
*
*     transform the integrals (2nd transformation including virtuals)
*
      a1 = cpulft(1)
      call start_time_period(TP_VB_TRAN)
      nmc = nsa + ncore 
      if (scfconv) nsa = nactiv
      lenact = nsa*(nsa+1)/2
*
      ks = igmem_alloc_inf(lenact,'vbqcscf.m','vbqcscf','ks',
     +                     IGMEM_DEBUG)
      kh = igmem_alloc_inf(lenact,'vbqcscf.m','vbqcscf','kh',
     +                     IGMEM_DEBUG)
      ng = lenact*(lenact+1)/2+1
*
*
*     (first integral = (00/00) = 0.0 => + 1)
*     don't use to much core in tran as virtual io may interfere
*     with normal io; If the user wants to commit suicide, let him
*     (jvl,2002)
*
      call transformvb(q(ks),q(kh),v)
*

      if ( nitscf.eq.1 .and. n2int+1 .ne. lenact*(lenact+1)/2+1) then
         write(iwr,'(a,I10,a,I10,a)')
     +    ' Truncating nr of 2-el integrals (Brillouin) to ',
     +    n2int+1,' from ',lenact*(lenact+1)/2 + 1,'.'
      end if
*
      kg = igmem_alloc_inf(n2int+1,'vbqcscf.m','vbqcscf',
     +                       'kg',IGMEM_DEBUG)
      call getin2(q(kg))
      call clredx
*
*     print integrals over orbitals if requested
*   
      if (iprint.ge.50) then
         write(iwr,*) ' 1-electron integrals over orbitals'
         call tripr14(q(kh),nsa)
         write(iwr,*) ' overlap matrix between orbitals'
         call tripr14(q(ks),nsa)
         if (iprint.ge.100000) then
            write(iwr,*) ' 2-electron integrals'
            call tripri(q(kg+1),lenact)
         end if
      end if
*
      a1 = cpulft(1) - a1
      t4indx = t4indx + a1
      n4indx = n4indx + 1
      call end_time_period(TP_VB_TRAN)
      call flushbuffer()
*
*     calculate hessian matrix
*
      a1 = cpulft(1)
      call start_time_period(TP_VB_ME)
*

      call wr15vb(ndum1,ndum2,ndum3,ndum4,nstruc,ndum6,ndum7,
     +            ndum8,ndum9,ndum10,ndum11,ndum12,'read')

      mbril = 0
      do i = 1, nscf
         mbril = mbril + nex(i)
      end do
      nbril = mbril + 1
*
      nnd = mbril + nstruc - 1
*
      khes = igmem_alloc_inf(nnd*nnd,'vbqcscf.m','vbqcscf','khes',
     +                        IGMEM_DEBUG)
      kgrad = igmem_alloc_inf(nnd+1,'vbqcscf.m','vbqcscf','kgrad',
     +                        IGMEM_DEBUG)
*
      call vclr(q(khes),1,nnd*nnd)
      call vclr(q(kgrad),1,nnd+1)
*
      call vbquad(q(ks),q(kh),q(kg),q(khes),q(kgrad),ci,nbasis,
     +            nactiv,nvirt,q)
*
      a1 = cpulft(1) - a1
      tmatre = tmatre + a1
      nmatre = nmatre + 1
      call end_time_period(TP_VB_ME)
*
      kupda = igmem_alloc_inf(nnd+1,'vbqcscf.m','vbqcscf','kupda',
     +                         IGMEM_DEBUG)
      if (ospci.and.oqcof) then
*
*        a super-ci iteration
*
         call fmove(q(kgrad),q(kupda),nbril)
*
      else
*
*        a newton-raphson iteration
*
         q(kupda) = 1.0d0
         if (ocong) then
*
*           conjugate gradient method to solve Ax = b
*           later will be switched on so that user have the choice between
*           conjugate gradient and inverse hessian to get the update
*           coefficient
*
            call vbgetupdate(q(khes),q(kgrad),q(kupda+1),nnd,q)
*
         else
*
*           invert the hessian matrix and get update coefficient by
*           multiplying 1/hessian with gradient
*
*
*           'subroutine osinv_vb' destroys the original matrix 
*           if hessian is needed (for storage or something else) then
*           make a copy (zahid) 
*
            ksv1 = igmem_alloc_inf(nnd,'vbin.m','scfina','ksv1',
     +                             IGMEM_DEBUG)
            ksv2 = igmem_alloc_inf(nnd,'vbin.m','scfina','ksv2',
     +                             IGMEM_DEBUG)
*
            fcrit = 1.0d-20
            call osinv_vb(q(khes),nnd,detval,fcrit,q(ksv1),q(ksv2))
            do i = 1, nnd 
               xx = 0.0d0
               do j = 1, nnd
                  xx = xx + q(khes-1+(i-1)*nnd+j) * q(kgrad-1+j)
               end do
               q(kupda+i) = - xx
            end do
            call gmem_free_inf(ksv2,'vbqcscf.m','vbqcscf','ksv2')
            call gmem_free_inf(ksv1,'vbqcscf.m','vbqcscf','ksv1')
         end if
      end if
*
*
*     update orbitals and check convergence
*
      knewv = igmem_alloc_inf(nbasis*nactiv,'vbqcscf.m','vbqcscf',
     +                        'knewv',IGMEM_DEBUG)
      kmaxmem  = igmem_max_memory()
      kmemscr  = kmaxmem/10
      kmemleft = kmaxmem - kmemscr
      kscr     = igmem_alloc_inf(kmemscr,'vbqcscf.m','vbqcscf',
     +                           'kscr',IGMEM_DEBUG)
*
*     change orbitals
*
      mstruc = nstruc
      nstruc = nbril
      call qcorbopt(v,nbasis,nactiv,q(knewv),q(kupda),q(kscr),kmemscr)
      nstruc = mstruc
*
      call gmem_free_inf(kscr,  'vbqcscf.m', 'vbqcscf', 'kscr'  )
      call gmem_free_inf(knewv, 'vbqcscf.m', 'vbqcscf', 'knewv' )
      call gmem_free_inf(kupda, 'vbqcscf.m', 'vbqcscf', 'kupda' )
      call gmem_free_inf(kgrad, 'vbqcscf.m', 'vbqcscf', 'kgrad' )
      call gmem_free_inf(khes,  'vbqcscf.m', 'vbqcscf', 'khes' )
      call gmem_free_inf(kg,    'vbqcscf.m', 'vbqcscf', 'kg'    )
      call gmem_free_inf(kh,    'vbqcscf.m', 'vbqcscf', 'kh'    )
      call gmem_free_inf(ks,    'vbqcscf.m', 'vbqcscf', 'ks'    )
*  
      if (unitary) then
         ksao = igmem_alloc_inf(nbasis*(nbasis+1)/2,'vbscf.m','vbscf', 
     +                          'ksao',IGMEM_DEBUG)
         kscr = igmem_alloc_inf(nbasis*(nbasis+1)/2,'vbscf.m','vbscf',
     +                        'kscr',IGMEM_DEBUG) 

         call get1e(q(ksao),dummy,'s',q(kscr))

         call gmem_free_inf(kscr,'vbscf.m','vbscf','kscr')

         kscr    = igmem_alloc_inf(nbasis*2,'vbscf.m','vbscf','kscr',
     +                              IGMEM_DEBUG)

         call normvc(v(ncore*nbasis+1),q(ksao),q(kscr),nbasis,
     +               nactiv,cridep)  

         call gmem_free_inf(kscr,'vbscf.m','vbscf','kscr')
         call gmem_free_inf(ksao,'vbscf.m','vbscf','ksao')
      end if
*
*
*     end of iterative loop. start on top 
*
      goto 42
*

100   continue

      icheckblock = 1
*
      ks = igmem_alloc_inf(lenact,'vbqcscf.m','vbqcscf','ks',
     +                     IGMEM_DEBUG)
      kh = igmem_alloc_inf(lenact,'vbqcscf.m','vbqcscf','kh',
     +                     IGMEM_DEBUG)

      nsa = nactiv
*
*     read info from tape
*
      call wr15vb(nelec,nalfa,nconf,ndets,nstruc,maxdet,maxstr,
     +            nwpack,ncoeff,imax,ndum11,ndoub,'read')
      norb = nsa + nvirt
c
      ksao = igmem_alloc_inf(nbasis*(nbasis+1)/2,'vbqcscf.m','vbqcscf',
     +                       'ksao',IGMEM_DEBUG)
      ksao2 = igmem_alloc_inf(nbasis*(nbasis+1)/2,'vbqcscf.m','vbqcscf',
     +                        'ksao2',IGMEM_DEBUG)
      kmaxmem  = igmem_max_memory()
      kmemscr  = kmaxmem/10
      kmemleft = kmaxmem - kmemscr
      kscr = igmem_alloc_inf(kmemscr,'vbqcscf.m','vbqcscf','kscr',
     +                       IGMEM_DEBUG)
c
c...  calculate s-matrix over orbitals (all but virtual orbitals)
c
      call get1e(q(ksao),dummy,'s',q(ksao2))
      call fmos(q(ks),q(ksao),v(ncore*nbasis+1),q(kscr),nsa,
     +          nbasis,crilow)

      call gmem_free_inf(kscr,'vbqcscf.m','vbqcscf','kscr')
      call gmem_free_inf(ksao2,'vbqcscf.m','vbqcscf','ksao2')
      call gmem_free_inf(ksao,'vbqcscf.m','vbqcscf','ksao')
c
c...  keep certain parameters
c
      hsmatt = 'vbfunction'
      oldncol = ncol
      ncol = ncore + nscf
      oldncore = ncore
      ncore = ncore + ndoub
      oldnsa = nsa
      nsa = ncol - ncore
      lenact = nsa*(nsa+1)/2
      oldncurt = ncurt
      ncurt = -1
      struc=.true.
      kmaxmem  = igmem_max_memory()
      kmemscr  = kmaxmem/10
      kmemleft = kmaxmem - kmemscr
      kscr = igmem_alloc_inf(kmemscr,'vbqcscf.m','vbqcscf','kscr',
     +                       IGMEM_DEBUG)
      call scfina(v,ci,q(ks),q(kscr),q)


      call gmem_free_inf(kscr,'vbqcscf.m','vbqcscf','kscr')
c
c...  find the orthogonality classes in the integrals
c
c...  why do we even do this??
c
      kmaxmem  = igmem_max_memory()
      kmemscr  = kmaxmem/10
      kmemleft = kmaxmem - kmemscr
      kiort = igmem_alloc_inf(norb,'vbqcscf.m','vbqcscf','kiort',
     +                       IGMEM_DEBUG)
      call sym1(q(ks),norb,q(kiort),northo,'noprint')


      kigro = igmem_alloc_inf(5*nconf,'vbqcscf.m','vbqcscf','kigro',
     +                        IGMEM_DEBUG)
      kpacd = igmem_alloc_inf(nwpack,'vbqcscf.m','vbqcscf','kpacd',
     +                        IGMEM_DEBUG)
      kndet = igmem_alloc_inf(nstruc,'vbqcscf.m','vbqcscf','kndet',
     +                        IGMEM_DEBUG)
      kidps = igmem_alloc_inf(ndets*nstruc,'vbqcscf.m','vbqcscf',
     +                       'kidps',IGMEM_DEBUG)
      kcoef = igmem_alloc_inf(ncoeff*ndets,'vbqcscf.m','vbqcscf',
     +                       'kcoef',IGMEM_DEBUG)
      kidet = igmem_alloc_inf(nelec*ndets,'vbqcscf.m','vbqcscf',
     +                       'kidet',IGMEM_DEBUG)
      kjdet = igmem_alloc_inf(nelec*ndets,'vbqcscf.m','vbqcscf',
     +                       'kjdet',IGMEM_DEBUG)
      kmaxmem  = igmem_max_memory()
      kmemscr  = kmaxmem/10
      kmemleft = kmaxmem - kmemscr
      kscr  = igmem_alloc_inf(kmemscr,'vbqcscf.m','vbqcscf','kscr',
     +                        IGMEM_DEBUG)
c
c...  read information from datatape
c
      call inform(nelec,nalfa,ndets,nstruc,q(kpacd),q(kndet),
     +            q(kidps),q(kcoef),ncoeff,q(kigro),ngroup,q(kiort),
     +            northo,q(kscr),kmemscr)
c
c...  generate psi0 on determinant basis
c
      call psi0det(ci,nstruc,q(kigro),ngroup,q(kcoef),ncoeff,
     +             ndets,q(kndet),q(kidps),q(kidet),q(kjdet),
     +             q(kpacd),nelec,nalfa,ndettot,q(kscr),kmemscr,
     +             'dontsaveondisk')
c
c...  move original parameters back
c
      struc=.false.
      hsmatt = 'full      '
      ncore = oldncore
      ncol = oldncol
      nsa = oldnsa
      lenact = nsa*(nsa+1)/2
      ncurt = oldncurt
      call gmem_free_inf(kscr,'vbqcscf.m','vbqcscf','kscr')
      call gmem_free_inf(kjdet,'vbqcscf.m','vbqcscf','kjdet')
      call gmem_free_inf(kidet,'vbqcscf.m','vbqcscf','kidet')
      call gmem_free_inf(kcoef,'vbqcscf.m','vbqcscf','kcoef')
      call gmem_free_inf(kidps,'vbqcscf.m','vbqcscf','kidps')
      call gmem_free_inf(kndet,'vbqcscf.m','vbqcscf','kndet')
      call gmem_free_inf(kpacd,'vbqcscf.m','vbqcscf','kpacd')
      call gmem_free_inf(kigro,'vbqcscf.m','vbqcscf','kigro')
      call gmem_free_inf(kiort,'vbqcscf.m','vbqcscf','kiort')
      call gmem_free_inf(kh,'vbqcscf.m','vbqcscf','kh')
      call gmem_free_inf(ks,'vbqcscf.m','vbqcscf','ks')
      return
      end

*
************************************************************************
      subroutine vbquad(supers,superh,superg,hes,grad,ci,nbasis,nsa,
     +                  nvirt,q)
************************************************************************
*
      implicit real*8 (a-h,o-z), integer (i-n)
*
*     dynamical core allocation and brillouin control routine
*     for dynamic memory allocation debugging purpose, contains i
*     IGMEM_ vars.
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
c...   turtle parameters
c
      integer melec,ncorq,maxex,maxact,mbonds,maxato,
     &        maxopa,ncases,maxcon,mmsocc,mxorbvb,maxspvb
c
      parameter (melec=300)
      parameter (ncorq=2000000)
      parameter (maxex=20000,maxact=150,mbonds=150,
     &           maxato=750,maxopa=1500)
      parameter (ncases=35)
c      ... maxcon in vbcrestr ...
c92      parameter (maxcon = 100)
      parameter (maxcon = 3000)
c      ... msocc max. # singly occuppieds ; was always enough (2009) ...
      parameter (mmsocc = 30)
c..   maxspvb  maximum # spinfunctions per configuration 
c..   (42 for 10 singlys and the ultimate answer)
      parameter (maxspvb=42) 
c
c     mxorbvb may be as large as allowed by n_8_16 (32767)
c
      parameter (mxorbvb=1000)
c
c... common to  sitmch from 8 to 16 bits orb ital packing
c...  n8_16 is # bits used for packing; i.e. 8 (255 orbitals) 16 (32000)
c...  n16_32 is the correspnding one for 16 vs 32 bits packin 
c...  GAMESS  parameters cannot be used as VB may have more MO's than AO's
c...  parameters are set in vbcrestr
c
      integer n8_16,n16_32,n340
      common/c8_16vb/n8_16,n16_32,n340
      common /brill/ nex(maxact) ,iex (maxex,maxact),
     &               nexw(maxact),iexw(maxex,maxact),
     &                            ieqsig(maxex,maxact),
     &               iequi(maxact),nequi,iredund(mxorbvb),nredund,
     &               inteq(maxact+1,maxact),ninteq(maxact),ninter,
     &               nscf,nullity,ortscf,schmi,reduce,autscf,
     &               exact(maxact),nbrill,idequi,
     &               ieqs(maxact),nexeq(maxact),iexeq(maxact*maxact)
      logical ortscf,schmi,reduce,autscf,exact
      integer nequi,nredund,ninter,nscf,nullity,nbrill,nex,nexw,
     &        iex,iexw,ieqsig,iequi,iredund,inteq,ninteq,idequi,
     &        ieqs,nexeq,iexeq
      logical ohmat
      integer nstruc,iprinv,ihfile,ihbloc,ipbloc,idumhs
      common/hsinfo/ohmat(3),nstruc,iprinv,ihfile,ihbloc,ipbloc,idumhs
*  
      common /icaselist/ icaselist(35)
      common /array/ narray,iarray(100)
      common /arnam/ anames(100)
*
      character*8 anames
*
      dimension supers(*),superh(*),superg(*),hes(*),grad(*),
     +          q(*),ci(*)
*
*
      call wr15vb(nelec,nalfa,nconf,ndets,nstruc,maxdet,maxstr,
     +            nwpack,ncoeff,imax,ndum11,ndum12,'read')
*
      npword = nwpack
      nextot = 0
      maxequ = 0
      do 10 i = 1, nequi
10    maxequ = max(iequi(i),maxequ)

      maxdet = max(ndets,2*ncoeff*maxequ)
      do 20 i = 1, nsa
         do 30 j = 1, nex(i)
            npword = npword + (4*ndets*nelec-1) / (64/n8_16) + 1
30       continue
         nextot = nextot + nex(i)
20    continue
      narray = 0
*
*     pacdet will contain the packed slater determinants. 
*
      kpacde = igmem_alloc_inf(nwpack,'vbqcscf.m','vbquad','kpacde',
     +                         IGMEM_DEBUG)
*
*     idet and jdet contain the unpacked slater determinants 
*     idet and cdet will contain all the determinants and their
*     coefficients, respectively, that make up the reference function.
*
      kidet  = igmem_alloc_inf((ndets+2)*nelec,'vbqcscf.m','vbquad',
     +                         'kidet',IGMEM_DEBUG)
      kcdet = igmem_alloc_inf(ncoeff+2,'vbqcscf.m','vbquad',
     +                        'kcdet',IGMEM_DEBUG)
*
*     idetall and coeff contains the determinants and their 
*     coefficients, respectively, on structure basis. 
*
      kidetall = igmem_alloc_inf(nelec*ncoeff,'vbqcscf.m','vbquad',
     +                        'kidetall',IGMEM_DEBUG)
      kcoeff = igmem_alloc_inf(ncoeff,'vbqcscf.m','vbquad',
     +                        'kcoeff',IGMEM_DEBUG)
*
*     ndetps contains the number of determinants per structure
*
      kndetp = igmem_alloc_inf(nstruc+max((nextot+1),nstruc),
     +                        'vbqcscf.m','vbquad','kndetp',IGMEM_DEBUG)
*
*     idetps contains the determinant numbers that define the structures
*
      kidetp = igmem_alloc_inf(ncoeff,'vbqcscf.m','vbquad','kidetp',
     +                         IGMEM_DEBUG)
*
*     detcomb is scratch for the calculation of the brill. matrix-
*     elements, dettot is scratch for the corresponding overlaps
*
      kdetco = igmem_alloc_inf(16*maxdet*ncoeff,'vbqcscf.m','vbquad',
     +                         'kdetco',IGMEM_DEBUG)
      kdetto = igmem_alloc_inf(16*maxdet*ncoeff,'vbqcscf.m','vbquad',
     +                         'kdetto',IGMEM_DEBUG)
*
*     the matrix element routine (matre3) is given the determinants via
*     icp and jcp (mind parity changes !).
*
      kicp   = igmem_alloc_inf(nelec,'vbqcscf.m','vbquad','kicp',
     +                         IGMEM_DEBUG)
      kjcp   = igmem_alloc_inf(nelec,'vbqcscf.m','vbquad','kjcp',
     +                         IGMEM_DEBUG)
*
*     ig will contain the blocking information per matrix element
*     (5 numbers per block). but at first (in symblo) it contains
*     5 numbers per orthogonality group => at most norb * 10 numbers !
*     (spin-orthogonality on top of spatial orthogonality => 2 * 5)
*     ig(1 to 4,0) is put to zero, ig(5,0)=1, this must be done to deal
*     with the case no alpha block occur in the overlap matrix (ialfa=0)
*      => at most norb * 10 + 5 numbers

      kig    = igmem_alloc_inf(10*(nsa+nvirt)+5,'vbqcscf.m','vbquad',
     +                        'kig',IGMEM_DEBUG)
*
*     igroup contains the number of determinants per bs, the number
*     of structures per bs (=1) and the starting position of idetps per
*     group.
*
      kigrou = igmem_alloc_inf(max(nconf,(nextot+1))*5,'vbqcscf.m',
     +                         'vbquad','kigrou',IGMEM_DEBUG)
*
*     ipos is the position array that is used to gather the integrals
*     per matrix element. it is also used to reorder the matrix
*     elements before they are written to disc (in writh)
*
*     weight contains the weight-factors of the integrals (cofactors)
*     per matrix element
*
*     first calculate the maximum number of positions in ipos/weight
*
                                nbeta = nelec - nalfa
*
*                               one electron part
*
                                nnnnn = nalfa ** 2 + nbeta ** 2
                                nwwww = nnnnn
*   
*                               second order alfa/beta terms apart
*
                                nxxxx = (nalfa * (nalfa-1) / 2)**2 +
     +                                  (nbeta * (nbeta-1) / 2)**2
                                nnnnn = nnnnn + nxxxx * 2
                                nwwww = nwwww + nxxxx
*
*                               mixed terms
*
                                nxxxx = (nalfa**2)*(nbeta**2)
                                nnnnn = nnnnn + nxxxx * 2
                                nwwww = nwwww + nxxxx
*
      kipos  = igmem_alloc_inf(max(nnnnn,nextot+1),'vbqcscf.m','vbquad',
     +                         'kipos',IGMEM_DEBUG)
      kweigh = igmem_alloc_inf(nwwww+1,'vbqcscf.m','vbquad','kweigh',
     +                         IGMEM_DEBUG)
*
*     g contains the integrals that make up the matrix element
*
      kg     = igmem_alloc_inf(nnnnn,'vbqcscf.m','vbquad','kg',
     +                         IGMEM_DEBUG)
*
*     iortho contains the orthogonality number per orbital
*
      kiorth = igmem_alloc_inf(nsa+nvirt,'vbqcscf.m','vbquad','kiorth',
     +                         IGMEM_DEBUG)
*
*     s contains the (biorthogonalised) overlap-matrix per matrix-
*     element
*
      ks     = igmem_alloc_inf(nalfa**2+nbeta**2+1,'vbqcscf.m','vbquad',
     +                        'ks',IGMEM_DEBUG)
*
*     hamil and overl contain the matrix representation of the
*     hamiltonian on a structure basis
*
      nnd = nstruc*(nstruc+1)/2
      khamil = igmem_alloc_inf(nnd,'vbqcscf.m','vbquad','khamil',
     +                         IGMEM_DEBUG)
      koverl = igmem_alloc_inf(nnd,'vbqcscf.m','vbquad','koverl',
     +                         IGMEM_DEBUG)
*
*     scr1 and scr2 are used in hamilt during the transformation of
*     hamil and overl. scr1/scr2/scr3 are used in matre3 using
*     always less than (nelec**2) words. scr1 is used in hamilt to
*     gather (reorder) matrix elements before they are written to disc.
*
      nnnnn  = max( nelec**2 , maxdet * 1 )
      mmmmm  = max( nnnnn    , nextot + 1 )
      kscr1  = igmem_alloc_inf(mmmmm,'vbqcscf.m','vbquad','kscr1',
     +                         IGMEM_DEBUG)
      kscr2  = igmem_alloc_inf(nnnnn,'vbqcscf.m','vbquad','kscr2',
     +                         IGMEM_DEBUG)
      kscr3  = igmem_alloc_inf(nelec*nelec,'vbqcscf.m','vbquad','kscr3',
     +                         IGMEM_DEBUG)
      kmaxmem = igmem_max_memory()
      kmemscr = kmaxmem/10
      kscr4  = igmem_alloc_inf(kmemscr,'vbqcscf.m','vbquad','kscr4',
     +                         IGMEM_DEBUG)
*
      call vbjanuuni(supers,superh,superg(2),hes,grad,ci,q(kpacde),
     +               q(kidet),q(kcdet),q(kidetall),q(kcoeff),q(kndetp),
     +               q(kidetp),q(kdetco),q(kdetto),q(kicp),q(kjcp),
     +               q(kig),q(kigrou),q(kipos),q(kweigh),q(kg),
     +               q(kiorth),q(ks),q(khamil),q(koverl),q(kscr1),
     +               q(kscr2),q(kscr3),q(kscr4),kmemscr)
*
      call gmem_free_inf(kscr4,    'vbqcscf.m', 'vbquad', 'kscr4'     )
      call gmem_free_inf(kscr3,    'vbqcscf.m', 'vbquad', 'kscr3'     )
      call gmem_free_inf(kscr2,    'vbqcscf.m', 'vbquad', 'kscr2'     )
      call gmem_free_inf(kscr1,    'vbqcscf.m', 'vbquad', 'kscr1'     )
      call gmem_free_inf(koverl,   'vbqcscf.m', 'vbquad', 'koverl'    )
      call gmem_free_inf(khamil,   'vbqcscf.m', 'vbquad', 'khamil'    )
      call gmem_free_inf(ks,       'vbqcscf.m', 'vbquad', 'ks'        )
      call gmem_free_inf(kiorth,   'vbqcscf.m', 'vbquad', 'kiorth'    )
      call gmem_free_inf(kg,       'vbqcscf.m', 'vbquad', 'kg'        )
      call gmem_free_inf(kweigh,   'vbqcscf.m', 'vbquad', 'kweigh'    )
      call gmem_free_inf(kipos,    'vbqcscf.m', 'vbquad', 'kipos'     )
      call gmem_free_inf(kigrou,   'vbqcscf.m', 'vbquad', 'kigrou'    )
      call gmem_free_inf(kig,      'vbqcscf.m', 'vbquad', 'kig'       )
      call gmem_free_inf(kjcp,     'vbqcscf.m', 'vbquad', 'kjcp'      )
      call gmem_free_inf(kicp,     'vbqcscf.m', 'vbquad', 'kicp'      )
      call gmem_free_inf(kdetto,   'vbqcscf.m', 'vbquad', 'kdetto'    )
      call gmem_free_inf(kdetco,   'vbqcscf.m', 'vbquad', 'kdetco'    )
      call gmem_free_inf(kidetp,   'vbqcscf.m', 'vbquad', 'kidetp'    )
      call gmem_free_inf(kndetp,   'vbqcscf.m', 'vbquad', 'kndetp'    )
      call gmem_free_inf(kcoeff,   'vbqcscf.m', 'vbquad', 'kcoeff'    )
      call gmem_free_inf(kidetall, 'vbqcscf.m', 'vbquad', 'kidetall'  )
      call gmem_free_inf(kcdet,    'vbqcscf.m', 'vbquad', 'kcdet'     )
      call gmem_free_inf(kidet,    'vbqcscf.m', 'vbquad', 'kidet'     )
      call gmem_free_inf(kpacde,   'vbqcscf.m', 'vbquad', 'kpacde'    )
*
      return
      end
*
************************************************************************
      subroutine vbjanuuni(supers,superh,superg,hes,grad,ci,pacdet,idet,
     +                     cdet,idetall,coeff,ndetps,idetps,detcomb,
     +                     dettot,icp,jcp,ig,igroup,ipos,weight,g,
     +                     iortho,s,hamil,overl,scr,scr2,scr3,q,lword)
************************************************************************
*
      implicit real*8 (a-h,o-z), integer (i-n)
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
c... common to  sitmch from 8 to 16 bits orb ital packing
c...  n8_16 is # bits used for packing; i.e. 8 (255 orbitals) 16 (32000)
c...  n16_32 is the correspnding one for 16 vs 32 bits packin 
c...  GAMESS  parameters cannot be used as VB may have more MO's than AO's
c...  parameters are set in vbcrestr
c
      integer n8_16,n16_32,n340
      common/c8_16vb/n8_16,n16_32,n340
c     8note* tractlt noo included in vbdens.m mkd1vb
      integer iprint,nbasis,ncol,nsa,ncurt,lenact,lenbas,nactiv
      integer num3,iblkq,ionsec,nbsort,isecdu,n2int,setintpos
      logical incsort,oprs
      integer nnsort,isort,ksort,maxsort
      real*8 scri
c
      common /tractlt/ scri,iprint,nbasis,ncol,nsa,ncurt,lenact,lenbas,
     &                 num3,iblkq,ionsec,nbsort,isecdu,n2int,nactiv,
     &                 incsort,nnsort,isort,ksort,maxsort,oprs
c
c...  incsort : indicates if incore sorting is on
c...  nnsort : # sortbuffers
c...  isort : sortbuffer done
c...  ksort : adres of sortbuffers in qq (igmem_alloc)
c...  maxsort : maximum sort buffer used
c
c     core :  frozen core energy for the whole job
c     core_i :  frozen core energy  in makeps0 
c
      real*8 core,pot_nuc,core_i
      integer nfil,nofil,iblvb,lblvb,iblhs
      common /ffile/ core,core_i,pot_nuc,nfil,
     1               nofil(20),iblvb(20),lblvb(20),iblhs
      logical ohmat
      integer nstruc,iprinv,ihfile,ihbloc,ipbloc,idumhs
      common/hsinfo/ohmat(3),nstruc,iprinv,ihfile,ihbloc,ipbloc,idumhs
      real*8  shiscf,criscf,evb,ebi,swave,percmax,perccur,brmax,
     &      dampvb,fixb0
      integer optcri,iprins,nitscf,maxscf,isecv,nitdiis,maxdiis,ntscf,
     &        ipcmt,maxscf_s
      logical perall,scfconv,p0core,unitary
      common /scftvb/dampvb,fixb0,shiscf,criscf,evb,ebi,swave,percmax,
     &               perccur,brmax,
     &               optcri,iprins,nitscf,maxscf,isecv,nitdiis,
     &               maxdiis,ntscf,ipcmt,maxscf_s,perall,scfconv,
     &               p0core,unitary
c
c...   turtle parameters
c
      integer melec,ncorq,maxex,maxact,mbonds,maxato,
     &        maxopa,ncases,maxcon,mmsocc,mxorbvb,maxspvb
c
      parameter (melec=300)
      parameter (ncorq=2000000)
      parameter (maxex=20000,maxact=150,mbonds=150,
     &           maxato=750,maxopa=1500)
      parameter (ncases=35)
c      ... maxcon in vbcrestr ...
c92      parameter (maxcon = 100)
      parameter (maxcon = 3000)
c      ... msocc max. # singly occuppieds ; was always enough (2009) ...
      parameter (mmsocc = 30)
c..   maxspvb  maximum # spinfunctions per configuration 
c..   (42 for 10 singlys and the ultimate answer)
      parameter (maxspvb=42) 
c
c     mxorbvb may be as large as allowed by n_8_16 (32767)
c
      parameter (mxorbvb=1000)
*
*     to signal which scf procedure to be used, e.g., super-ci or Newton-Raphson
*     default is super-ci
*
      real*8 scalhs, sci2nr, delta, crmax
      integer nits, nitsci, nege, nneg
      logical oqcscf, ocong, ospci, oinve, ograd, ocorr, oauto,
     +        oiter, oqcof
      common / quadratic / nitsci,  nits,  nege, nneg, scalhs, sci2nr,
     +                     delta, crmax, oqcscf, ocong, ospci, oinve, 
     +                     ograd, ocorr, oauto, oiter, oqcof
      logical oprop,orunprop,opol,omag,onosym,odebug,oignore,ovbci
      logical ocomplex,ouncoupled,ovbci2
      integer nex_prop,iex_prop,ndubbel
      common/vbproper/oprop,orunprop,opol,omag,onosym,odebug,ovbci,
     &oignore,ocomplex,ouncoupled,ovbci2,ndubbel,
     &nex_prop(maxact),iex_prop(maxex,maxact)
      common /brill/ nex(maxact) ,iex (maxex,maxact),
     &               nexw(maxact),iexw(maxex,maxact),
     &                            ieqsig(maxex,maxact),
     &               iequi(maxact),nequi,iredund(mxorbvb),nredund,
     &               inteq(maxact+1,maxact),ninteq(maxact),ninter,
     &               nscf,nullity,ortscf,schmi,reduce,autscf,
     &               exact(maxact),nbrill,idequi,
     &               ieqs(maxact),nexeq(maxact),iexeq(maxact*maxact)
      logical ortscf,schmi,reduce,autscf,exact
      integer nequi,nredund,ninter,nscf,nullity,nbrill,nex,nexw,
     &        iex,iexw,ieqsig,iequi,iredund,inteq,ninteq,idequi,
     &        ieqs,nexeq,iexeq
c
c     twice contains info on pauli-principle related restrictions
c     ("doubly","redorb","iflip"), and super, which indicates wether
c     internal excitations are allowed for (not so if super-mcscf)
c     alos ifreez/nfreez orbitals to be frozen
      integer inforb,ninfor,idoubly,ndoubly,isingly,nsingly,kaart,
     &        ioror,igno_sh,igno_sel
      logical super,super_hybrid,ofrzmo,super_cas,super_act
      common /twice/ inforb(maxact*5,maxact),ninfor(maxact),
     &               idoubly(maxact),ndoubly,isingly(maxact),nsingly,
     &               kaart(mxorbvb+maxact),ioror(maxact,maxact),super,
     &               super_hybrid,igno_sh(maxato),ofrzmo(maxact),
     &               super_cas,super_act,igno_sel
c...
C...
c...  atonam for names of atoms (for mulliken population analysis,
c...  and automatic atomic (strictly localised) scf-procedure)
c...  atom is used her to denote fragments
c...  nopa contains the # orbitals per atom, iopa contains them
c...  mullik indicates if the analysis is requested, hybry indicates
c...  wether the scf-procedure just should used atom-centered ao's
c...  per atomic scf orbital. iacat contains the active (scf) orbital
c...  per atom, nacat contains the number of the per atom.
c...  atomao : fragment numbers / ao (may be more)
c...           then stored as frag1, maxato*frag2, maxato**2*frag3 etc.
c     ifrzat : hybrid handling
c              0 : normal, 1 frozen, 2 super
c...
      character*8 atoms
      common /atonam/ atoms(maxato)
      integer atomao,nopa,iopa,natom,iacat,nacat,ifrzat
      common /infato/ nopa(maxato),iopa(maxopa,maxato),natom,
     &                iacat(maxact,maxato),nacat(maxato),
     &                atomao(mxorbvb),ifrzat(maxato),
     &                mullik,hybry,clean,hguess,guesao
      logical mullik,hybry,clean,hguess,guesao
c
c...  frozen core definition for TURTLE
c
      real*8 dxyz
      integer ncore,mapcie,mapiee,ncore_i,nsplice
      parameter (nsplice=512)
      common/splice/dxyz(3),ncore,ncore_i,
     1              mapcie(nsplice),mapiee(nsplice)
c...   common meant for all kind of (fixed) vb criteria
c...   critor .. criterium to determine orthogonality
c...   crilow .. criterium  for orthogonalization (e.g.jacoby in lowdin)
c...   cridep .. criterium for dependency (on c**2)
c...   symcri .. criterium for symmetry equivalence
c...   remcri .. criterium for including excitation on basis of 1-e ints
c...   crihyb .. max hybrid contam. allowed in matrix and (*100) vector
c...   cribri .. min. size of bcoeff to avoid annihilation
c...   crinor .. min. norm of a vector
c...   cribub .. a meaningful difference (in a bubble sort)
c...   criign .. size of matrix elements that really may be ignored
c...   cripop .. criterium for pipek/mizek localization
c...   normally critor = 1.0d-14, crilow = critor/10, cridep = critor*10
c...            symcri = critor*100 remcri = 1.0d-10 crihyb = 1.0d-14
c...            cribri = 1.0d-8 crinor = critor*100  cribub = 1.0d-10
c...            criign = 1.0d-44 cripop 1.0d-10
      real*8 critor,crilow,cridep,symcri,remcri,crihyb,cribri,
     1     crinor,cribub,criign,cripop
      common/vbcri/ critor,crilow,cridep,symcri,remcri,crihyb,cribri,
     1              crinor,cribub,criign,cripop
c
c .. common for  diis in vb /vbdiist/
c
      real*8 evbpre
      integer iwdiis,ntdiis,madiis,midiis,kvb7_vnew
      logical ovbdiis
      common /vbdiist/ evbpre,iwdiis,ntdiis,madiis,midiis,
     1                 kvb7_vnew,ovbdiis
*
*
      dimension supers(*), superh(*), superg(*), hes(*), grad(*),
     +          ci(*), pacdet(*), idet(*), cdet(*), idetall(*),
     +          coeff(*), ndetps(*), idetps(*), detcomb(*),
     +          dettot(*), icp(*), jcp(*), ig(*), igroup(*), ipos(*),
     +          weight(*), g(*), iortho(*), s(*), hamil(*), overl(*),
     +          scr(*), scr2(*), scr3(*), q(*)
*
      norb = nsa
*
*     find the orthogonality classes in the integrals
*
      call sym1(supers,norb,iortho,northo,'noprint')
*
*     read information from datatape
*
      call inform(nelec,nalfa,ndets,nstruc,pacdet,ndetps,idetps,cdet,
     +            ncoeff,igroup,ngroup,iortho,northo,q,lword)
*       
*     get wave function 
*     coeff(*)   = coefficients of determinants on structure basis
*     idetall(*) = all determinants on structure basis
*     cdet(*)    = coefficients of wave function as one structure 
*     idet(*)    = all determinants of wave functions as one structure
*
      kjdet  = igmem_alloc_inf((ndets+2)*nelec,'vbqcscf.m','vbjanuuni',
     +                         'kjdet',IGMEM_DEBUG)

      call getpsi0(ci,nstruc,igroup,ngroup,cdet,coeff,ndetps,idetps,
     +             idet,idetall,q(kjdet),pacdet,nelec,nalfa,ndet,ndtot,
     +             ncoeff,q,lword)
      call gmem_free_inf(kjdet, 'vbqcscf.m', 'vbjanuuni', 'kjdet' )
*
      if (ndtot.ne.ncoeff) call 
     +    vberr('no. of determinants not equal to no. of coefficients')
*
*     check overlap of old wavefunction with new one
*
      call overlon(nelec,nalfa,igroup,ngroup,idet,q,lword)
*
*     allocate space for excited determinants and their coefficients
*
*     for doubly occupied orbitals single excitations (i->j) will result
*     in a maximum of 2*ndet excited determinants while double 
*     excitations (i->j,k->l) will result in a maximum of 4*ndet excited
*     determinants. (in case of unitary max 8*ndet)
*     + scratch space for subroutine match( )
*
      memval = 8*ndet*nelec + (8*ndet+ndet+2)*nelec + 
     +         8*ndet*2 + (8*ndet+2)*nelec * 2 + 2
      if (memval.ge.lword)
     + call corfait(memval,lword,' subroutine vbjanuuni before kiscr')
*
      kiscr = igmem_alloc_inf(8*ndet*nelec,'vbqcscf.m',
     +                        'vbjanuuni','kiscr', IGMEM_DEBUG)
      kiscr2 = igmem_alloc_inf((8*ndet+ndet+2)*nelec,'vbqcscf.m',
     +                         'vbjanuuni','kiscr2', IGMEM_DEBUG)
      kiscr3 = igmem_alloc_inf(2,'vbqcscf.m',
     +                         'vbjanuuni','kiscr3', IGMEM_DEBUG)
      kcdet2 = igmem_alloc_inf(8*ndet,'vbqcscf.m','vbjanuuni',
     +                        'kcdet2', IGMEM_DEBUG)
      kcdet3 = igmem_alloc_inf(8*ndet,'vbqcscf.m','vbjanuuni',
     +                        'kcdet3', IGMEM_DEBUG)
      kidet2 = igmem_alloc_inf((8*ndet+2)*nelec,'vbqcscf.m',
     +                        'vbjanuuni','kidet2', IGMEM_DEBUG)
      kidet3 = igmem_alloc_inf((8*ndet+2)*nelec,'vbqcscf.m',
     +                         'vbjanuuni','kidet3', IGMEM_DEBUG)
*
      mbril = 0
      do i = 1, nscf
        mbril = mbril + nex(i)
      end do
*
      nbril = mbril + 1
      nd = mbril + nstruc - 1
      ndim = mbril * (mbril + 1) / 2
      mci = nstruc
*

      ig(1) = 0
      ig(2) = 0
      ig(3) = 0
      ig(4) = 0
      ig(5) = 1
*
*
      call hathad(dd,ds,icp,jcp,supers,superh,superg,ig(6),ipos,weight,
     +            nalfa,scr,scr2,scr3,s,g,nelec,iortho,northo,nbody,
     +            ndet,idet,cdet,detcomb,dettot)
*
*
      e0 = dd/ds
      dnorm = ds
*
      memval = memval + ndim + mci*(mci+1)/2 + mbril*mci + mbril*2 +
     +         (nbril*(nbril+1)/2)* 2
      if (memval.ge.lword)
     + call corfait(memval,lword,' subroutine vbjanuuni before kcihes')
*
      kcihes = igmem_alloc_inf(mci*(mci+1)/2,'vbqcscf.m',
     +                        'vbjanuuni','kcihes', IGMEM_DEBUG)
      kciorbh = igmem_alloc_inf(mbril*mci,'vbqcscf.m',
     +                          'vbjanuuni','kciorbh', IGMEM_DEBUG)
*
      korbh = igmem_alloc_inf(mbril,'vbqcscf.m','vbjanuuni','korbh',
     +                        IGMEM_DEBUG)
      kbhamil = igmem_alloc_inf(nbril*(nbril+1)/2,'vbqcscf.m',
     +                         'vbjanuuni','kbhamil', IGMEM_DEBUG)
      kboverl = igmem_alloc_inf(nbril*(nbril+1)/2,'vbqcscf.m',
     +                         'vbjanuuni','kboverl', IGMEM_DEBUG)
*
      call vclr(q(kbhamil),1,nbril*(nbril+1)/2)
      call vclr(q(kboverl),1,nbril*(nbril+1)/2)
*
      q(kbhamil) = dd
      q(kboverl) = ds
*
      call vborbder(hes,q(korbh),icp,jcp,supers,
     +              superh,superg,grad,ig(6),ipos,s,g,iortho,weight,
     +              scr,scr2,scr3,detcomb,dettot,q(kiscr),cdet,
     +              q(kcdet2),q(kcdet3),idet,q(kidet2),q(kidet3),
     +              ndet,northo,nelec,nalfa,mbril,q(kbhamil),
     +              q(kboverl),q(kiscr2),q(kiscr3),nd)
*
      if (ospci.and.oqcof) then
*
         if (oiter) nits = nits + 1
*
*        get the orbital updates from super-ci
*        gradient vector is not required any more
*        return the update coefficient in grad(*) 
*
*
         memval = memval + nbril*nbril*2 + nbril
         if (memval.ge.lword) call corfait(memval,lword,
     +                    ' subroutine vbjanuuni before kvec')
*
         kvec = igmem_alloc_inf(nbril*nbril,'vbqcscf.m',
     +                          'vbjanuuni','kvec', IGMEM_DEBUG)
         ke = igmem_alloc_inf(nbril,'vbqcscf.m',
     +                         'vbjanuuni','ke', IGMEM_DEBUG)
         kscr = igmem_alloc_inf(nbril*nbril,'vbqcscf.m',
     +                          'vbjanuuni','kscr', IGMEM_DEBUG)
*
         call jacobs(q(kbhamil),nbril,q(kboverl),q(kvec),q(ke),
     +               2, 1.0d-20, q(kscr) )
*
         do i = 1, nbril
            if (q(kvec).lt.0.0d0 ) then
               grad(i) = - q(kvec-1+i)
            else
               grad(i) = q(kvec-1+i)
            end if
         end do
*
         crmax = 0.0d0
         yy = 0.0d0
         do i = 2, nbril
            yy  = dabs(grad(i))/dabs(grad(1))
            if (yy.gt.crmax) crmax = yy
         end do
*
         call gmem_free_inf(kscr,   'vbqcscf.m', 'vbjanuuni', 'kscr'   )
         call gmem_free_inf(ke,     'vbqcscf.m', 'vbjanuuni', 'ke'     )
         call gmem_free_inf(kvec,   'vbqcscf.m', 'vbjanuuni', 'kvec'   )
         call gmem_free_inf(kboverl,'vbqcscf.m', 'vbjanuuni', 'kboverl')
         call gmem_free_inf(kbhamil,'vbqcscf.m', 'vbjanuuni', 'kbhamil')
         call gmem_free_inf(korbh,  'vbqcscf.m', 'vbjanuuni', 'korbh'  )
         call gmem_free_inf(kciorbh,'vbqcscf.m', 'vbjanuuni', 'kciorbh')
         call gmem_free_inf(kcihes, 'vbqcscf.m', 'vbjanuuni', 'kcihes' )
         call gmem_free_inf(kidet3, 'vbqcscf.m', 'vbjanuuni', 'kidet3' )
         call gmem_free_inf(kidet2, 'vbqcscf.m', 'vbjanuuni', 'kidet2' )
         call gmem_free_inf(kcdet3, 'vbqcscf.m', 'vbjanuuni', 'kcdet3' )
         call gmem_free_inf(kcdet2, 'vbqcscf.m', 'vbjanuuni', 'kcdet2' )
         call gmem_free_inf(kiscr3, 'vbqcscf.m', 'vbjanuuni', 'kiscr3' )
         call gmem_free_inf(kiscr2, 'vbqcscf.m', 'vbjanuuni', 'kiscr2' )
         call gmem_free_inf(kiscr,  'vbqcscf.m', 'vbjanuuni', 'kiscr'  )
*
         return
      end if
*
*     proceed with quadratic scf
*     super-ci h and s are not needed anymore.
*
      call gmem_free_inf(kboverl, 'vbqcscf.m', 'vbjanuuni', 'kboverl' )
      call gmem_free_inf(kbhamil, 'vbqcscf.m', 'vbjanuuni', 'kbhamil' )
*
      if (nstruc.gt.1) then
*        for single structure wave function, we don't need to 
*        calculate the ci and mixed part of the hessian
*     
         memval = memval - (nbril*(nbril+1)/2) * 2 + mci * 2
         if (memval.ge.lword) call corfait(memval,lword,
     +                       ' subroutine vbjanuuni before kcih2')
*     
         kcis2 = igmem_alloc_inf(mci,'vbqcscf.m',
     +                             'vbjanuuni','kcis2', IGMEM_DEBUG)
*     
         call vbcider(q(kcihes),q(kcis2),hamil,overl,ndetps,coeff,
     +                idetall,icp,jcp,supers,superh,superg,grad(nbril),
     +                ig(6),ipos,s,g,iortho,weight,scr,scr2,scr3,
     +                detcomb,dettot,q(kiscr),cdet,q(kcdet2),q(kcdet3),
     +                idet,q(kidet2),q(kidet3),idetps,ndet,northo,nelec,
     +                nalfa,mci,e0,dnorm,q(kiscr2),q(kiscr3))
*     
         call vbciorbder(q(kciorbh),q(korbh),q(kcis2),ndetps,coeff,
     +                   idetall,icp,jcp,supers,superh,superg,ig(6),
     +                   ipos,s,g,iortho,weight,scr,scr2,scr3,detcomb,
     +                   dettot,q(kiscr),cdet,q(kcdet2),q(kcdet3),
     +                   idet,q(kidet2),q(kidet3),idetps,ndet,northo,
     +                   nelec,nalfa,mci,e0,dnorm,mbril,q(kiscr2),
     +                   q(kiscr3))
*     
         call gmem_free_inf(kcis2, 'vbqcscf.m', 'vbjanuuni', 'kcis2' )
*                 
         memval = memval - mbril * 2 - mci + mci * mci * 2 
         if (memval.ge.lword) call corfait(memval,lword,
     +                       ' subroutine vbjanuuni before kvec')
*           
         kvec = igmem_alloc_inf(mci*mci,'vbqcscf.m',
     +                           'vbjanuuni','kvec', IGMEM_DEBUG)
         ke   = igmem_alloc_inf(mci,'vbqcscf.m',
     +                           'vbjanuuni','ke', IGMEM_DEBUG)
         kscr = igmem_alloc_inf(mci*mci,'vbqcscf.m',
     +                          'vbjanuuni','kscr', IGMEM_DEBUG)
*     
         if (iprint.gt.10000) then
            write(iwr,*)' hamiltonian matrix '
            call tripr14(hamil,mci)
            write(iwr,*)' overlap matrix '
            call tripr14(overl,mci)
         end if
*     
         call jacobs(hamil,mci,overl,q(kvec),q(ke),2,1.0d-20,q(kscr))
*     
c         write(iwr,80) 
c80       format(/1x,' ** eigenvalues and eigenvectors of Hessian **',/)
c         call prvc(q(kvec),mci,mci,q(ke),'v','a')
c         evbnew = q(ke) + core
c         print*, ' vb energy ',evbnew
*     
         call gmem_free_inf(kscr, 'vbqcscf.m', 'vbjanuuni', 'kscr' )
         call gmem_free_inf(ke,   'vbqcscf.m', 'vbjanuuni', 'ke'   )
*     
*        transform ci hessian
*     
         memval = memval - mci - mci * mci + mci*(mci-1) + 
     +            mci*(mci-1)/2 + mci*(mci-1)+(mci-1)*(mci-1)
         if (memval.ge.lword) call corfait(memval,lword,
     +                       ' subroutine vbjanuuni before kdgr')
*     
         kdgr = igmem_alloc_inf(mci*(mci-1),'vbqcscf.m',
     +                           'vbjanuuni','kdgr', IGMEM_DEBUG)
         kscr = igmem_alloc_inf(mci*(mci-1)/2,'vbqcscf.m',
     +                           'vbjanuuni','kscr', IGMEM_DEBUG)
         kscr2 = igmem_alloc_inf(mci*(mci-1)+(mci-1)*(mci-1),
     +                           'vbqcscf.m','vbjanuuni','kscr2',
     +                           IGMEM_DEBUG)
*     
*        kdgr => dagger
*     
         call formdagger(q(kvec+mci),q(kdgr),mci,mci-1)
*     
*     
         call mult11(q(kcihes),q(kscr),q(kscr2),mci-1,mci,q(kvec+mci),
     +               q(kdgr))
*     
         do i = 1, mci*(mci-1)/2
            q(kcihes+i-1) = q(kscr+i-1)
         end do
*     
         call gmem_free_inf(kscr2,'vbqcscf.m','vbjanuuni','kscr2')
         call gmem_free_inf(kscr,'vbqcscf.m','vbjanuuni','kscr')
*     
*        transform ci/orbital part 
*     
         memval = memval - mci*(mci-1)/2 - mci*(mci-1)+(mci-1)*(mci-1) + 
     +            (mci-1)*mbril
         if (memval.ge.lword) call corfait(memval,lword,
     +                       ' subroutine vbjanuuni before kscr')
*     
         kscr = igmem_alloc_inf((mci-1)*mbril,'vbqcscf.m',
     +                           'vbjanuuni','kscr', IGMEM_DEBUG)
*     
*     
         do i = 1, mbril
            call transformmix(q(kciorbh+mci*(i-1)),q(kscr+(mci-1)*
     +                        (i-1)),q(kdgr),mci-1,mci)
         end do
*     
         call vclr(q(kciorbh),1,mbril*mci)
*     
         do i = 1, mbril*(mci-1)
            q(kciorbh+i-1) = q(kscr+i-1)
         end do
*
         call gmem_free_inf(kscr,'vbqcscf.m','vbjanuuni','kscr')
         call gmem_free_inf(kdgr,'vbqcscf.m','vbjanuuni','kdgr' )
         call gmem_free_inf(kvec,'vbqcscf.m','vbjanuuni','kvec' )
      end if
*
      memval = memval - (mci-1)*mbril - mci*(mci-1) - mci*mci
*
      mci = mci - 1
*
*
*     form complete hessian-matrix
*
      ndim = mbril + mci

      call fillhes(hes,ndim,mbril,q(kcihes),mci,q(kciorbh))
*
      if (iprint.gt.10000) then
         write(iwr,*) ' complete hessian-matrix '
         call prsq14(hes,ndim,ndim,ndim)
      end if
*
*     diagonalise hessian
*
      memval = memval + ndim*(ndim+1)/2 + ndim*ndim + ndim*4 + 1
      if (memval.ge.lword) call corfait(memval,lword,
     +                    ' subroutine vbjanuuni before ktri')
*
      ktri = igmem_alloc_inf(ndim*(ndim+1)/2,'vbqcscf.m',
     +                        'vbjanuuni','ktri', IGMEM_DEBUG)
*
*     copy hes to q(ktri) triangle matrix
*
      call fillmatrix(hes,q(ktri),ndim)
*
      khvec = igmem_alloc_inf(ndim*ndim,'vbqcscf.m',
     +                        'vbjanuuni','khvec', IGMEM_DEBUG)
      ke = igmem_alloc_inf(ndim,'vbqcscf.m',
     +                        'vbjanuuni','ke', IGMEM_DEBUG)
      kscr = igmem_alloc_inf(ndim+1,'vbqcscf.m',
     +                        'vbjanuuni','kscr', IGMEM_DEBUG)
      kscr2 = igmem_alloc_inf(ndim*2,'vbqcscf.m',
     +                        'vbjanuuni','kscr2', IGMEM_DEBUG)
*
100   continue
      call jacodiag(q(ktri),ndim,q(khvec),q(ke),q(kscr),q(kscr2))
c
c      write(iwr,'(/,1x,a,/)')
c     +          'eigen values and vectors of hessian matrix'
c      call prvc(q(khvec),ndim,ndim,q(ke),'v','a')
c
c      print *,' eigenvalues of hessian',ndim
c      write(iwr,'(1p10e13.5)') (q(ke+i-1),i=1,ndim)
c      write(iwr,'(a,f10.5)') ' eigenvalue of hessian ',q(ke)

      nneg = 0
      do i = 1, ndim
         if (q(ke-1+i).lt.0.00001) nneg = nneg + 1
      end do
*
      if (nneg.gt.0) then
*
         call makehpositive(hes,grad,q(ke),ndim)
*
*        check once again just to make sure hes is postive-definite
*
c         call fillmatrix(hes,q(ktri),ndim)
c         go to 100
      end if
*
      call gmem_free_inf(kscr2,   'vbqcscf.m', 'vbjanuuni', 'kscr2'   )
      call gmem_free_inf(kscr,    'vbqcscf.m', 'vbjanuuni', 'kscr'    )
      call gmem_free_inf(ke,      'vbqcscf.m', 'vbjanuuni', 'ke'      )
      call gmem_free_inf(khvec,   'vbqcscf.m', 'vbjanuuni', 'khvec'   )
      call gmem_free_inf(ktri,    'vbqcscf.m', 'vbjanuuni', 'ktri'    )
      call gmem_free_inf(korbh,   'vbqcscf.m', 'vbjanuuni', 'korbh'   )
      call gmem_free_inf(kciorbh, 'vbqcscf.m', 'vbjanuuni', 'kciorbh' )
      call gmem_free_inf(kcihes,  'vbqcscf.m', 'vbjanuuni', 'kcihes'  )
      call gmem_free_inf(kidet3,  'vbqcscf.m', 'vbjanuuni', 'kidet3'  )
      call gmem_free_inf(kidet2,  'vbqcscf.m', 'vbjanuuni', 'kidet2'  )
      call gmem_free_inf(kcdet3,  'vbqcscf.m', 'vbjanuuni', 'kcdet3'  )
      call gmem_free_inf(kcdet2,  'vbqcscf.m', 'vbjanuuni', 'kcdet2'  )
      call gmem_free_inf(kiscr3,  'vbqcscf.m', 'vbjanuuni', 'kiscr3'  )
      call gmem_free_inf(kiscr2,  'vbqcscf.m', 'vbjanuuni', 'kiscr2'  )
      call gmem_free_inf(kiscr,   'vbqcscf.m', 'vbjanuuni', 'kiscr'   )
*
c      call gmem_check_guards(' subroutine vbjanuuni end')
*
      return
      end
*
************************************************************************
      subroutine vborbder(hes,orbh,icp,jcp,smat,hmat,gmat,
     +                    grad,ig,ipos,s,g,iortho,weight,scr,scr2,scr3,
     +                    scr4,scr5,iscr,cdet,cdet2,cdet3,idet,idet2,
     +                    idet3,ndet,northo,nelec,nalfa,mbril,bhamil,
     +                    boverl,iscr2,iscr3,nd)
************************************************************************
*
      implicit real*8 (a-h,o-z), integer (i-n)
*
*     generate orbital hessian matrix
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
c
c...   turtle parameters
c
      integer melec,ncorq,maxex,maxact,mbonds,maxato,
     &        maxopa,ncases,maxcon,mmsocc,mxorbvb,maxspvb
c
      parameter (melec=300)
      parameter (ncorq=2000000)
      parameter (maxex=20000,maxact=150,mbonds=150,
     &           maxato=750,maxopa=1500)
      parameter (ncases=35)
c      ... maxcon in vbcrestr ...
c92      parameter (maxcon = 100)
      parameter (maxcon = 3000)
c      ... msocc max. # singly occuppieds ; was always enough (2009) ...
      parameter (mmsocc = 30)
c..   maxspvb  maximum # spinfunctions per configuration 
c..   (42 for 10 singlys and the ultimate answer)
      parameter (maxspvb=42) 
c
c     mxorbvb may be as large as allowed by n_8_16 (32767)
c
      parameter (mxorbvb=1000)
      common /brill/ nex(maxact) ,iex (maxex,maxact),
     &               nexw(maxact),iexw(maxex,maxact),
     &                            ieqsig(maxex,maxact),
     &               iequi(maxact),nequi,iredund(mxorbvb),nredund,
     &               inteq(maxact+1,maxact),ninteq(maxact),ninter,
     &               nscf,nullity,ortscf,schmi,reduce,autscf,
     &               exact(maxact),nbrill,idequi,
     &               ieqs(maxact),nexeq(maxact),iexeq(maxact*maxact)
      logical ortscf,schmi,reduce,autscf,exact
      integer nequi,nredund,ninter,nscf,nullity,nbrill,nex,nexw,
     &        iex,iexw,ieqsig,iequi,iredund,inteq,ninteq,idequi,
     &        ieqs,nexeq,iexeq
      real*8  shiscf,criscf,evb,ebi,swave,percmax,perccur,brmax,
     &      dampvb,fixb0
      integer optcri,iprins,nitscf,maxscf,isecv,nitdiis,maxdiis,ntscf,
     &        ipcmt,maxscf_s
      logical perall,scfconv,p0core,unitary
      common /scftvb/dampvb,fixb0,shiscf,criscf,evb,ebi,swave,percmax,
     &               perccur,brmax,
     &               optcri,iprins,nitscf,maxscf,isecv,nitdiis,
     &               maxdiis,ntscf,ipcmt,maxscf_s,perall,scfconv,
     &               p0core,unitary
c
c     twice contains info on pauli-principle related restrictions
c     ("doubly","redorb","iflip"), and super, which indicates wether
c     internal excitations are allowed for (not so if super-mcscf)
c     alos ifreez/nfreez orbitals to be frozen
      integer inforb,ninfor,idoubly,ndoubly,isingly,nsingly,kaart,
     &        ioror,igno_sh,igno_sel
      logical super,super_hybrid,ofrzmo,super_cas,super_act
      common /twice/ inforb(maxact*5,maxact),ninfor(maxact),
     &               idoubly(maxact),ndoubly,isingly(maxact),nsingly,
     &               kaart(mxorbvb+maxact),ioror(maxact,maxact),super,
     &               super_hybrid,igno_sh(maxato),ofrzmo(maxact),
     &               super_cas,super_act,igno_sel
c...
*
*     to signal which scf procedure to be used, e.g., super-ci or Newton-Raphson
*     default is super-ci
*
      real*8 scalhs, sci2nr, delta, crmax
      integer nits, nitsci, nege, nneg
      logical oqcscf, ocong, ospci, oinve, ograd, ocorr, oauto,
     +        oiter, oqcof
      common / quadratic / nitsci,  nits,  nege, nneg, scalhs, sci2nr,
     +                     delta, crmax, oqcscf, ocong, ospci, oinve, 
     +                     ograd, ocorr, oauto, oiter, oqcof
*
*
      dimension cdet(*),cdet2(*),cdet3(*),idet(*),idet2(*),
     +          idet3(*),icp(*),jcp(*),smat(*),hmat(*),gmat(*),
     +          iortho(*),weight(*),scr(*),scr2(*),scr3(*),scr4(*),
     +          scr5(*),ig(*),ipos(*),s(*),g(*),
     +          iscr(*),grad(*), orbh(*),
     +          bhamil(*),boverl(*),iscr2(*),iscr3(2),hes(nd,nd)
*
*
      ind(i,j) = max0(i,j)*(max0(i,j)-1)/2+min0(i,j)
*
      e0 = bhamil(1)/boverl(1)
      dnorm = boverl(1)
*
      if (shiscf.ne.0.0d0) bhamil(1) = bhamil(1) - shiscf * boverl(1)
*
      nbril = mbril + 1
      nbdim = nbril * (nbril + 1) / 2 - 1
*
      call vclr(orbh,1,mbril)
*
c      call printpsi(ndet,cdet,idet,nelec)
*
c     tstart = cpulft(1)
*
*     first calculate the gradient part
*
*     < Psi(0) | H - E0 | Psi(ij) > 
*
*
c      write(iwr,*) 'orbital gradient < Psi(0) | H - E0 | Psi(ij) > '

      ij = 0
      do i = 1, nscf
         ifrm = iex(1,i)
         do j = 2, nex(i) + 1
            itoj = iex(j,i)
            ij = ij + 1
            jj = ij * ( ij + 1 ) / 2 + 1
            call vbexcite(ndet,cdet,idet,cdet2,idet2,nelec,ifrm,itoj,
     +                    ndet2,nalfa)
c            call printpsi(ndet2,cdet2,idet2,nelec)
            call breorb(idet2,ndet2,nelec,nalfa,iortho,cdet2)
*
            call match(idet,idet(nelec*ndet+1),ndet,nelec,nalfa,
     +                 iscr3(1),iscr2(1))

            call match(idet2,idet2(nelec*ndet2+1),ndet2,nelec,nalfa,
     +                 iscr3(2),iscr2(ndet+1))
*
            if (iscr3(1).eq.1.and.iscr3(2).eq.1) then
*
               call hathab(dd,ds,icp,jcp,smat,hmat,gmat,ig,ipos,weight,
     +                     nalfa,scr,scr2,scr3,s,g,nelec,iortho,northo,
     +                     nbody,ndet,ndet2,idet,idet2,iscr2(1),
     +                     iscr2(ndet+1),cdet,cdet2,scr4,scr5)
            else
*
               call hatham(dd,ds,icp,jcp,smat,hmat,gmat,ig,ipos,weight,
     +                     nalfa,scr,scr2,scr3,s,g,nelec,iortho,northo,
     +                     nbody,ndet,ndet2,idet,idet2,cdet,cdet2,scr4,
     +                     scr5)
            end if
*
            bhamil(jj) = dd
            boverl(jj) = ds
            orbh(ij) = dd - e0 * ds
            grad(ij)  = 2.0d0 * ( dd - e0 * ds ) / dnorm
*
            if (shiscf.ne.0.0d0)  then
               bhamil(jj+1) = bhamil(jj+1) - shiscf * boverl(jj+1)
            end if
*
c         write(iwr,11) ifrm,itoj,ij, orbh(ij), dd, ds
c11       format(1x,' gradient h and s ',3i5,5x,3f20.14)
*
         end do
      end do
*
*
c      tgend = cpulft(1)
c      tg = tgend - tstart
c      write(iwr,*) 'time  < Psi(0) | H - E0 | Psi(ij) >  ', tg
*
*     double excitations 2.0d0 * < Psi(ij) | H - E0 * S | Psi(kl) >
*
      ij = 0
      do i = 1, nscf
         ifrm = iex(1,i)
         do j = 2, nex(i) + 1
            itoj = iex(j,i)
            ij = ij + 1
            call vbexcite(ndet,cdet,idet,cdet2,idet2,nelec,ifrm,itoj,
     +                    ndet2,nalfa)
            kl = 0
            do k = 1, nscf
               kfrm = iex(1,k)
               do l = 2, nex(k) + 1
                  ktol = iex(l,k)
                  kl = kl + 1
                  if (kl.le.ij) then
*
                     call breorb(idet2,ndet2,nelec,nalfa,iortho,cdet2)
*
                     if (ij.eq.kl) then
*
                        call hathad(dd,ds,icp,jcp,smat,hmat,gmat,ig,
     +                              ipos,weight,nalfa,scr,scr2,scr3,
     +                              s,g,nelec,iortho,northo,nbody,
     +                              ndet2,idet2,cdet2,scr4,scr5)
                     else
*
                        call vbexcite(ndet,cdet,idet,cdet3,idet3,nelec,
     +                                 kfrm,ktol,ndet3,nalfa)
*
                        call breorb(idet3,ndet3,nelec,nalfa,iortho,
     +                              cdet3)
                        call match(idet2,idet2(nelec*ndet2+1),ndet2,
     +                             nelec,nalfa,iscr3(1),iscr2(1))
                        call match(idet3,idet3(nelec*ndet3+1),ndet3,
     +                             nelec,nalfa,iscr3(2),iscr2(ndet2+1))

                        if (iscr3(1).eq.1.and.iscr3(2).eq.1) then
*
                           call hathab(dd,ds,icp,jcp,smat,hmat,gmat,ig,
     +                                 ipos,weight,nalfa,scr,scr2,scr3,
     +                                 s,g,nelec,iortho,northo,nbody,
     +                                 ndet2,ndet3,idet2,idet3,
     +                                 iscr2(1),iscr2(ndet2+1),cdet2,
     +                                 cdet3,scr4,scr5)
                        else
*
                           call hatham(dd,ds,icp,jcp,smat,hmat,gmat,ig,
     +                                 ipos,weight,nalfa,scr,scr2,scr3,
     +                                 s,g,nelec,iortho,northo,nbody,
     +                                 ndet2,ndet3,idet2,idet3,cdet2,
     +                                 cdet3,scr4,scr5)
                        end if
                     end if
*
c       write(iwr,'(a,i3,a,i3,a,i3,a,i3,3f20.15,a,i3)')
c     + '< Psi(ij) |H-E0| Psi(kl) >',ifrm,' -> ',itoj,',',kfrm,' ->',
c     +    ktol,dd, ds, 2.0d0*(dd-e0*ds)/dnorm,' index ',ndet3
*
                     bhamil(ind(1+ij,1+kl)) = dd
                     boverl(ind(1+ij,1+kl)) = ds
                     hes(ij,kl) = 2.0d0 * (dd-e0*ds)/dnorm
*
                  end if
               end do
            end do
         end do
      end do
*
*
      brmax = 0.0d0
      do i = 1, mbril
         brmax = max(brmax,dabs(grad(i)))
      end do
      brmax = brmax/2.0d0
*
c      tsend = cpulft(1)
c      tsingle = tsend - tgend
c      write(iwr,*) 'time  < Psi(ij) | H - E0 | Psi(kl) >  ', tsingle
*
*     various switches to control the next iteration
*
      if (ospci) then
         if (ograd.and.brmax.le.sci2nr) oqcof = .false.
         if (ocorr.and.crmax.le.sci2nr) oqcof = .false.
         if (oauto) then
            if (brmax.le.sci2nr.or.crmax.le.sci2nr) then
                oqcof = .false.
            else
                oqcof = .true.
            end if
         end if
         if (oiter.and.nits.ge.nitsci) oqcof = .false.
      end if
*
      if (oqcof) return
*
*     double excitations 2.0d0 * < Psi(0) | H - E0 | Psi(ij,kl) >
*
      ij = 0
      do i = 1, nscf
         ifrm = iex(1,i)
         do j = 2, nex(i) + 1
            itoj = iex(j,i)
            ij = ij + 1
            kl = 0
            do k = 1, nscf
               kfrm = iex(1,k)
               do l = 2, nex(k) + 1
                  ktol = iex(l,k)
                  kl = kl + 1
                  if (kl.le.ij.and.(ifrm.ne.ktol.or.itoj.ne.kfrm)) then 
*
*                    ifrm = ktol and itoj = kfrm (i.e., 1 -> 2, 2 -> 1) 
*                    produces no brillouin state
*
                     call vbexcite(ndet,cdet,idet,cdet2,idet2,nelec,
     +                             ifrm,itoj,ndet2,nalfa)
*
                     call vbexcite(ndet2,cdet2,idet2,cdet3,idet3,
     +                             nelec,kfrm,ktol,ndet3,nalfa)
*
c                     call printpsi(ndet3,cdet3,idet3,nelec)

                     call breorb(idet3,ndet3,nelec,nalfa,iortho,cdet3)
*
                     call match(idet,idet(nelec*ndet+1),ndet,
     +                          nelec,nalfa,iscr3(1),iscr2(1))

                     call match(idet3,idet3(nelec*ndet3+1),ndet3,
     +                          nelec,nalfa,iscr3(2),iscr2(ndet+1))
*
                     if (iscr3(1).eq.1.and.iscr3(2).eq.1) then
*
                        call hathab(dd,ds,icp,jcp,smat,hmat,gmat,ig,
     +                              ipos,weight,nalfa,scr,scr2,scr3,
     +                              s,g,nelec,iortho,northo,nbody,
     +                              ndet,ndet3,idet,idet3,iscr2(1),
     +                              iscr2(ndet+1),cdet,cdet3,scr4,scr5)
                     else
*
                        call hatham(dd,ds,icp,jcp,smat,hmat,gmat,ig,
     +                              ipos,weight,nalfa,scr,scr2,scr3,
     +                              s,g,nelec,iortho,northo,nbody,ndet,
     +                              ndet3,idet,idet3,cdet,cdet3,scr4,
     +                              scr5)
                     end if
*
                     hes(ij,kl) = hes(ij,kl) + 2.0d0*(dd-e0*ds)/dnorm
*
c       write(iwr,'(a,i3,a,i3,a,i3,a,i3,3f20.15,a,i3)') 
c     + '< Psi(0) |H-E0| Psi(ij,kl) >',ifrm,' -> ',itoj,',',kfrm,' ->',
c     +    ktol,dd, ds, 2.0d0*(dd-e0*ds)/dnorm,' index ',ndet3
c        else
c          write(iwr,'(a,i3,a,i3,a,i3,a,i3)') 
c     +    ' from ',ifrm,' -> ',itoj,',',kfrm,' ->',ktol
c         end if
                  end if
               end do
            end do
         end do
      end do
*
*
*     generate final 2nd derivative matrix of orbitals
*
*     2.0d0 * < Psi(ij) | H - E0 | Psi(kl) > + 
*     2.0d0 * < Psi(0) | H - E0 | Psi(ij,kl) > -
*     4.0d0 * < Psi(0) | H - E0 | Psi(ij) > < Psi(0) | Psi(kl) > -
*     4.0d0 * < Psi(0) | H - E0 | Psi(kl) > < Psi(0) | Psi(ij) >
*
*
      ij = 0
      do i = 1, nscf
         ifrm = iex(1,i)
         do j = 2, nex(i) + 1
            itoj = iex(j,i)
            ij = ij + 1
            kl = 0
            sij = boverl(ij*(ij+1)/2+1)
            hij = bhamil(ij*(ij+1)/2+1) - e0 * sij
            do k = 1, nscf
               kfrm = iex(1,k)
               do l = 2, nex(k) + 1
                  ktol = iex(l,k)
                  kl   = kl + 1
                  if (kl.le.ij) then
                     skl = boverl(kl*(kl+1)/2+1)
                     hkl = bhamil(kl*(kl+1)/2+1) - e0 * skl
                     hes(ij,kl) = hes(ij,kl) -  4.0d0 * (
     +                                     hij * skl + hkl * sij )
                  end if
               end do
            end do
         end do
      end do
*
c      call prsq14(hes,nd,nd,nd)
c      call gmem_check_guards(' after orbital hessian complete')
*
      return
      end
*
************************************************************************
      subroutine vbcider(cihes,cis2,hamil,overl,indetps,coeff,
     +                   jidets,icp,jcp,smat,hmat,gmat,grad,ig,ipos,
     +                   s,g,iortho,weight,scr,scr2,scr3,scr4,scr5,
     +                   iscr,cdet,cdet2,cdet3,idet,idet2,idet3,
     +                   idetps,ndet,northo,nelec,nalfa,mstruc,e0,
     +                   dnorm,iscr2,iscr3)
************************************************************************
*
      implicit real*8 (a-h,o-z), integer (i-n)
*
*     generate ci hessian matrix
*
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
c
c...   turtle parameters
c
      integer melec,ncorq,maxex,maxact,mbonds,maxato,
     &        maxopa,ncases,maxcon,mmsocc,mxorbvb,maxspvb
c
      parameter (melec=300)
      parameter (ncorq=2000000)
      parameter (maxex=20000,maxact=150,mbonds=150,
     &           maxato=750,maxopa=1500)
      parameter (ncases=35)
c      ... maxcon in vbcrestr ...
c92      parameter (maxcon = 100)
      parameter (maxcon = 3000)
c      ... msocc max. # singly occuppieds ; was always enough (2009) ...
      parameter (mmsocc = 30)
c..   maxspvb  maximum # spinfunctions per configuration 
c..   (42 for 10 singlys and the ultimate answer)
      parameter (maxspvb=42) 
c
c     mxorbvb may be as large as allowed by n_8_16 (32767)
c
      parameter (mxorbvb=1000)
      common /brill/ nex(maxact) ,iex (maxex,maxact),
     &               nexw(maxact),iexw(maxex,maxact),
     &                            ieqsig(maxex,maxact),
     &               iequi(maxact),nequi,iredund(mxorbvb),nredund,
     &               inteq(maxact+1,maxact),ninteq(maxact),ninter,
     &               nscf,nullity,ortscf,schmi,reduce,autscf,
     &               exact(maxact),nbrill,idequi,
     &               ieqs(maxact),nexeq(maxact),iexeq(maxact*maxact)
      logical ortscf,schmi,reduce,autscf,exact
      integer nequi,nredund,ninter,nscf,nullity,nbrill,nex,nexw,
     &        iex,iexw,ieqsig,iequi,iredund,inteq,ninteq,idequi,
     &        ieqs,nexeq,iexeq
      logical ohmat
      integer nstruc,iprinv,ihfile,ihbloc,ipbloc,idumhs
      common/hsinfo/ohmat(3),nstruc,iprinv,ihfile,ihbloc,ipbloc,idumhs
*
*
      dimension cihes(*),cis2(*),
     +          hamil(*),overl(*),indetps(*),coeff(*),jidets(*),
     +          cdet(*),cdet2(*),cdet3(*),idet(*),idet2(*),
     +          idet3(*),icp(*),jcp(*),smat(*),hmat(*),gmat(*),
     +          iortho(*),weight(*),scr(*),scr2(*),scr3(*),scr4(*),
     +          scr5(*),ig(*),ipos(*),s(*),g(*),
     +          idetps(*),iscr(*),grad(*),iscr2(*),iscr3(2)
*
*     < phi_k | H - E0 | phi_l >
*
*     jidets (*) =  all determinants of psi0 
*     coeff(*)   =  coefficients of determinants
*     indetps(*) =  number of determinant per structure
*
      call vclr(cihes,1,nstruc*(nstruc+1)/2)
*
      call geths(hamil,nstruc,ihfile,ihbloc,'h')
      call geths(overl,nstruc,ihfile,ihbloc,'s')
*
      do i = 1, nstruc*(nstruc+1)/2
         cihes(i) = 2.0d0 * ( hamil(i) - e0 * overl(i) )
      end do
*     
      do k = 1, nstruc
*
         call getphik(k,indetps,coeff,jidets,ndet2,idet2,cdet2,nelec)

         call breorb(idet2,ndet2,nelec,nalfa,iortho,cdet2)

         call onlys(ds,icp,jcp,smat,hmat,gmat,ig,ipos,weight,nalfa,
     +              scr,scr2,scr3,s,g,nelec,iortho,northo,nbody,
     +              ndet,ndet2,idet,idet2,cdet,cdet2,scr4,scr5)
*
         cis2(k) = ds
         grad(k) = 0.0d0

c         write(iwr,12) k, dd, ds, ds*e0
c12       format(1x,'gradient h and s ',i4,5x,3f20.12)
*
      end do
*
*
      return
      end

************************************************************************
      subroutine vbciorbder(ciorbh,orbh,cis2,indetps,
     +                      coeff,jidets,icp,jcp,smat,hmat,gmat,ig,
     +                      ipos,s,g,iortho,weight,scr,scr2,scr3,
     +                      scr4,scr5,iscr,cdet,cdet2,cdet3,idet,
     +                      idet2,idet3,idetps,ndet,northo,
     +                      nelec,nalfa,nstruc,e0,dnorm,mbril,
     +                      iscr2,iscr3)

************************************************************************
*     
      implicit real*8 (a-h,o-z), integer (i-n)
*
*     generate orbital-ci mixed hessian 
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
c
c...   turtle parameters
c
      integer melec,ncorq,maxex,maxact,mbonds,maxato,
     &        maxopa,ncases,maxcon,mmsocc,mxorbvb,maxspvb
c
      parameter (melec=300)
      parameter (ncorq=2000000)
      parameter (maxex=20000,maxact=150,mbonds=150,
     &           maxato=750,maxopa=1500)
      parameter (ncases=35)
c      ... maxcon in vbcrestr ...
c92      parameter (maxcon = 100)
      parameter (maxcon = 3000)
c      ... msocc max. # singly occuppieds ; was always enough (2009) ...
      parameter (mmsocc = 30)
c..   maxspvb  maximum # spinfunctions per configuration 
c..   (42 for 10 singlys and the ultimate answer)
      parameter (maxspvb=42) 
c
c     mxorbvb may be as large as allowed by n_8_16 (32767)
c
      parameter (mxorbvb=1000)
      common /brill/ nex(maxact) ,iex (maxex,maxact),
     &               nexw(maxact),iexw(maxex,maxact),
     &                            ieqsig(maxex,maxact),
     &               iequi(maxact),nequi,iredund(mxorbvb),nredund,
     &               inteq(maxact+1,maxact),ninteq(maxact),ninter,
     &               nscf,nullity,ortscf,schmi,reduce,autscf,
     &               exact(maxact),nbrill,idequi,
     &               ieqs(maxact),nexeq(maxact),iexeq(maxact*maxact)
      logical ortscf,schmi,reduce,autscf,exact
      integer nequi,nredund,ninter,nscf,nullity,nbrill,nex,nexw,
     &        iex,iexw,ieqsig,iequi,iredund,inteq,ninteq,idequi,
     &        ieqs,nexeq,iexeq
*
*
      dimension ciorbh(*), orbh(*),cis2(*),
     +          indetps(*),coeff(*),jidets(*),
     +          cdet(*),cdet2(*),cdet3(*),idet(*),idet2(*),
     +          idet3(*),icp(*),jcp(*),smat(*),hmat(*),gmat(*),
     +          iortho(*),weight(*),scr(*),scr2(*),scr3(*),scr4(*),
     +          scr5(*),ig(*),ipos(*),s(*),g(*),
     +          idetps(*),iscr(*),iscr2(*),iscr3(2)
*
*
      call vclr(ciorbh,1,nstruc*mbril)
      ij = 0
      do i = 1, nscf
         ifrm = iex(1,i)
         do j = 2, nex(i) + 1
            itoj = iex(j,i)
            ij = ij + 1
            call vbexcite(ndet,cdet,idet,cdet2,idet2,nelec,ifrm,itoj,
     +                    ndet2,nalfa)
            do k = 1, nstruc
               call getphik(k,indetps,coeff,jidets,ndet3,idet3,
     +                        cdet3,nelec)

               call breorb(idet2,ndet2,nelec,nalfa,iortho,cdet2)
               call breorb(idet3,ndet3,nelec,nalfa,iortho,cdet3)
*
               call match(idet2,idet2(nelec*ndet2+1),ndet2,
     +                    nelec,nalfa,iscr3(1),iscr2(1))

               call match(idet3,idet3(nelec*ndet3+1),ndet3,
     +                    nelec,nalfa,iscr3(2),iscr2(ndet2+1))
*
               if (iscr3(1).eq.1.and.iscr3(2).eq.1) then
*
                  call hathab(dd,ds,icp,jcp,smat,hmat,gmat,ig,
     +                        ipos,weight,nalfa,scr,scr2,scr3,
     +                        s,g,nelec,iortho,northo,nbody,
     +                        ndet2,ndet3,idet2,idet3,
     +                        iscr2(1),iscr2(ndet2+1),cdet2,
     +                        cdet3,scr4,scr5)
               else
                  call hatham(dd,ds,icp,jcp,smat,hmat,gmat,ig,ipos,
     +                        weight,nalfa,scr,scr2,scr3,s,g,nelec,
     +                        iortho,northo,nbody,ndet2,ndet3,idet2,
     +                        idet3,cdet2,cdet3,scr4,scr5)
               end if
*
               ciorbh(nstruc*(ij-1)+k) = (2.0d0*(dd-e0*ds))/dnorm
            end do
         end do
      end do
*
      do k = 1, nstruc
         ij = 0
         do i = 1, nscf
             ifrm = iex(1,i)
             do j = 2, nex(i) + 1
                itoj = iex(j,i)
                ij = ij + 1

                call getphik(k,indetps,coeff,jidets,ndet3,idet3,
     +                         cdet3,nelec)
                call vbexcite(ndet3,cdet3,idet3,cdet2,idet2,nelec,
     +                        ifrm,itoj,ndet2,nalfa)

                call breorb(idet2,ndet2,nelec,nalfa,iortho,cdet2)
*
                call match(idet,idet(nelec*ndet+1),ndet,
     +                     nelec,nalfa,iscr3(1),iscr2(1))

                call match(idet2,idet2(nelec*ndet2+1),ndet2,
     +                     nelec,nalfa,iscr3(2),iscr2(ndet+1))
*
                if (iscr3(1).eq.1.and.iscr3(2).eq.1) then
*
                   call hathab(dd,ds,icp,jcp,smat,hmat,gmat,ig,
     +                         ipos,weight,nalfa,scr,scr2,scr3,
     +                         s,g,nelec,iortho,northo,nbody,
     +                         ndet,ndet2,idet,idet2,
     +                         iscr2(1),iscr2(ndet+1),cdet,
     +                         cdet2,scr4,scr5)
                else
*
                   call hatham(dd,ds,icp,jcp,smat,hmat,gmat,ig,ipos,
     +                         weight,nalfa,scr,scr2,scr3,s,g,nelec,
     +                         iortho,northo,nbody,ndet,ndet2,idet,
     +                         idet2,cdet,cdet2,scr4,scr5)
                end if
*
                ciorbh(nstruc*(ij-1)+k) = ciorbh(nstruc * (ij-1)+k) +
     +                                     (2.0d0*(dd-e0*ds))/dnorm
             end do
         end do
      end do
*
*
      do k = 1, nstruc
         ij = 0
         do i = 1, nscf
            ifrm = iex(1,i)
            do j = 2, nex(i) + 1
               itoj = iex(j,i)
               ij = ij + 1
               ciorbh(nstruc*(ij-1)+k) = ciorbh(nstruc*(ij-1)+k)
     +                                   - 4.0d0 * orbh(ij) * cis2(k) 
            end do
         end do
      end do
*
      return
      end
*
************************************************************************
      subroutine vbgetupdate(hes,grad,x,ndim,q)
************************************************************************
*              
      implicit real*8 (a-h,o-z), integer (i-n)
*              
*     prepares to invoke the conjugate gradient method to solve Ax = b
*    
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
*                
*    
      dimension hes(*), grad(*), x(*), q(*)
*    
*    
*       
      kp = igmem_alloc_inf(ndim, 'vbqcscf.m',
     +                          'vbgetupdate', 'kp', IGMEM_DEBUG )
      kq = igmem_alloc_inf(ndim, 'vbqcscf.m',
     +                          'vbgetupdate', 'kq', IGMEM_DEBUG )
      kr = igmem_alloc_inf(ndim, 'vbqcscf.m',
     +                          'vbgetupdate', 'kr', IGMEM_DEBUG )
*
      call conjgrad(hes,grad,x,q(kp),q(kq),q(kr),ndim)
*
      call gmem_free_inf(kr, 'vbqcscf.m', 'vbgetupdate', 'kr' )
      call gmem_free_inf(kq, 'vbqcscf.m', 'vbgetupdate', 'kq' )
      call gmem_free_inf(kp, 'vbqcscf.m', 'vbgetupdate', 'kp' )
*
      return
      end
*
************************************************************************
      subroutine conjgrad(a,g,x,p,q,r,nd)
************************************************************************
*
*     conjugate gradient method to solve Ax = b    where
*
*     A = a (ndim, ndim)  input
*     b = g (ndim)        input
*     x = x (ndim)        output
*
*
      implicit real*8 (a-h,o-z), integer (i-n)
*
*
      dimension a(nd,nd), g(nd), x(nd), p(nd), q(nd), r(nd)
*
*
      iter = 0
      crit = 1.0d-15
*
*     initial guess for x is 0.0d0
*
      do i = 1, nd
         x(i) = 0.0d0
      end do
*
*     r = g - Ax  
*
      do i = 1, nd 
         xx = 0.0d0
         do j = 1, nd
            xx = xx + a(i,j) * x(j)
         end do
*
         r(i) = - g(i) - xx
         p(i) = - g(i) - xx
      end do
*
      rho  = ddot( nd,r,1,r,1 )
*
      do i = 1, nd
         xx = 0.0d0
         do j = 1, nd
            xx = xx + a(i,j) * p(j)
         end do
         q(i) = xx
      end do
*
      alpha = rho / ddot( nd,p,1,q,1)
*
      do i = 1, nd
         x(i) = x(i) + alpha * p(i)
         r(i) = r(i) - alpha * q(i)
      end do
*
10    iter = iter + 1  
*
      rho0 = rho
      rho  = ddot( nd,r,1,r,1 )
      do j = 1, nd
         p(j) = r(j) + ( rho / rho0 ) * p(j)
      end do
      do j = 1,  nd
         xx = 0.0d0
         do k = 1, nd
            xx = xx + a(j,k) * p(k)
         end do
         q(j) = xx
      end do
      alpha = rho / ddot( nd,p,1,q,1 )
      do j = 1, nd
         x(j) = x(j) + alpha * p(j)
         r(j) = r(j) - alpha * q(j)
      end do
      check = dsqrt( ddot(nd,r,1,r,1) )
*
*     normally nd iterations are enough to converge but
*     we are using too tight criterion.
*
      if (check.le.crit.or.iter.gt.2*nd) go to 20
*
      go to 10
*
20    continue
*
*
      if (iter.gt.2*nd.and.check.gt.crit) then
         write(*,30) iter, check, x
         write(iwr,*) ' '
         call vberr( 'conjugate gradient method failed to converge' )
      end if
*
30    format(1x,'iteration',i4,f20.15,/,(t15,5f20.15))
*
      return
      end
*
*************************************************************************
      subroutine makehpositive(hes,grad,eig,nd)
*************************************************************************
*
      implicit real*8 (a-h,o-z), integer (i-n)
*
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
*
*     to signal which scf procedure to be used, e.g., super-ci or Newton-Raphson
*     default is super-ci
*
      real*8 scalhs, sci2nr, delta, crmax
      integer nits, nitsci, nege, nneg
      logical oqcscf, ocong, ospci, oinve, ograd, ocorr, oauto,
     +        oiter, oqcof
      common / quadratic / nitsci,  nits,  nege, nneg, scalhs, sci2nr,
     +                     delta, crmax, oqcscf, ocong, ospci, oinve, 
     +                     ograd, ocorr, oauto, oiter, oqcof
*
      dimension hes(nd,nd), grad(nd), eig(nd)
*
      if (scalhs.le.0.0d0) scalhs = 1.5d0
      shift = dabs(eig(1)) + scalhs * ddot(nd,grad,1,grad,1)
      if (nneg.gt.0) then
         write(iwr,10) nneg
10       format(/1x,' => Hessian matrix has ',i3,
     +           ' negative or (near) zero eigenvalues <=')
c         write(iwr,'(a,f10.5)') 'scaling factor ', shift
      end if
*
      do i = 1, nd
         hes(i,i) = hes(i,i) + shift
      end do
*
      return
      end
*
************************************************************************
      subroutine matrixmul(a,u,c,n)
************************************************************************
*
      implicit real*8 (a-h,o-z), integer (i-n)
*
      dimension a(n,n), u(n,n), c(n,n)
*
      do i = 1, n
         do j = 1, n
            xx = 0.0d0
            do k = 1, n
               xx = xx + a(i,k) * u(k,j)
            end do
            c(i,j) = xx
         end do
      end do
*
      return
      end
*
************************************************************************
      subroutine qcorbopt(v,nb,nact,vnew,bi,q,lword)
************************************************************************
*
      implicit real*8 (a-h,o-z), integer (i-n)
*
*     use info in bi-vector to update orbitals
*     if converged return scfconv = .true.
*     v are all the orbitals (including core)
*     bi update-vector
*     q scratch space
*     nstruc (in hsinfo) now have the dimension of bi(*)
*
c
c .. common for  diis in vb /vbdiist/
c
      real*8 evbpre
      integer iwdiis,ntdiis,madiis,midiis,kvb7_vnew
      logical ovbdiis
      common /vbdiist/ evbpre,iwdiis,ntdiis,madiis,midiis,
     1                 kvb7_vnew,ovbdiis
c     8note* tractlt noo included in vbdens.m mkd1vb
      integer iprint,nbasis,ncol,nsa,ncurt,lenact,lenbas,nactiv
      integer num3,iblkq,ionsec,nbsort,isecdu,n2int,setintpos
      logical incsort,oprs
      integer nnsort,isort,ksort,maxsort
      real*8 scri
c
      common /tractlt/ scri,iprint,nbasis,ncol,nsa,ncurt,lenact,lenbas,
     &                 num3,iblkq,ionsec,nbsort,isecdu,n2int,nactiv,
     &                 incsort,nnsort,isort,ksort,maxsort,oprs
c
c...  incsort : indicates if incore sorting is on
c...  nnsort : # sortbuffers
c...  isort : sortbuffer done
c...  ksort : adres of sortbuffers in qq (igmem_alloc)
c...  maxsort : maximum sort buffer used
      logical ohmat
      integer nstruc,iprinv,ihfile,ihbloc,ipbloc,idumhs
      common/hsinfo/ohmat(3),nstruc,iprinv,ihfile,ihbloc,ipbloc,idumhs
c
c...   turtle parameters
c
      integer melec,ncorq,maxex,maxact,mbonds,maxato,
     &        maxopa,ncases,maxcon,mmsocc,mxorbvb,maxspvb
c
      parameter (melec=300)
      parameter (ncorq=2000000)
      parameter (maxex=20000,maxact=150,mbonds=150,
     &           maxato=750,maxopa=1500)
      parameter (ncases=35)
c      ... maxcon in vbcrestr ...
c92      parameter (maxcon = 100)
      parameter (maxcon = 3000)
c      ... msocc max. # singly occuppieds ; was always enough (2009) ...
      parameter (mmsocc = 30)
c..   maxspvb  maximum # spinfunctions per configuration 
c..   (42 for 10 singlys and the ultimate answer)
      parameter (maxspvb=42) 
c
c     mxorbvb may be as large as allowed by n_8_16 (32767)
c
      parameter (mxorbvb=1000)
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
      common/vbpert/crivarpert(4),shiftp,fpdoc,fpuoc,
     &              nvarp,ifirsp,idelp,npert,
     &              ntype,itype(mxorbvb),
     &              minfock,maxfock,minpert,maxpert,minvar,maxvarp,
     &              minfockd,maxfockd,
     &              fckp,ovbper,fockdiagonal
      INTEGER nvarp,ifirsp,idelp,npert,ntype,itype
      INTEGER minfock,maxfock,minpert,maxpert,minvar,maxvarp
      INTEGER minfockd,maxfockd
      logical ovbper,fckp,fockdiagonal
      real*8 crivarpert,shiftp,fpdoc,fpuoc
c description:
c   ntype -> 1 - specify optimisation per brillioun structure
c            2 - specify optimisation per orbital
c            3 - specify optimisation per kind (type excitation)
c   itype ->
c
c          ntype = 1
c            itype(1..maxex) - optimisation type for each (brillouin) structure
c
c          ntype = 2
c            itype(1..mxorbvb) - optimisation type for each orbital
c            
c          ntype = 3 
c            itype(1) - optimisation type for doubly to variably occupied
c            itype(2) - optimisation type for doubly to unoccupied
c            itype(3) - optimisation type for variably to variably occupied
c            itype(4) - optimisation type for variably to unoccupied
c            
c          ntype = 4
c            itype(1..mxorbvb) - optimisation type for each orbital (to be
c                               specified by clsbril
c
c   contents of itype array:
c            itype(x)=1 for Super CI optimalisation
c            itype(x)=2 for perturbation theory optimalisation
c            itype(x)=4 perturbation theory  determines choice between Super and Pert
c
c   crivarpert ~: critaria for itype=4 case per class if applicable

c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common /brill/ nex(maxact) ,iex (maxex,maxact),
     &               nexw(maxact),iexw(maxex,maxact),
     &                            ieqsig(maxex,maxact),
     &               iequi(maxact),nequi,iredund(mxorbvb),nredund,
     &               inteq(maxact+1,maxact),ninteq(maxact),ninter,
     &               nscf,nullity,ortscf,schmi,reduce,autscf,
     &               exact(maxact),nbrill,idequi,
     &               ieqs(maxact),nexeq(maxact),iexeq(maxact*maxact)
      logical ortscf,schmi,reduce,autscf,exact
      integer nequi,nredund,ninter,nscf,nullity,nbrill,nex,nexw,
     &        iex,iexw,ieqsig,iequi,iredund,inteq,ninteq,idequi,
     &        ieqs,nexeq,iexeq
      real*8  shiscf,criscf,evb,ebi,swave,percmax,perccur,brmax,
     &      dampvb,fixb0
      integer optcri,iprins,nitscf,maxscf,isecv,nitdiis,maxdiis,ntscf,
     &        ipcmt,maxscf_s
      logical perall,scfconv,p0core,unitary
      common /scftvb/dampvb,fixb0,shiscf,criscf,evb,ebi,swave,percmax,
     &               perccur,brmax,
     &               optcri,iprins,nitscf,maxscf,isecv,nitdiis,
     &               maxdiis,ntscf,ipcmt,maxscf_s,perall,scfconv,
     &               p0core,unitary
c...   common meant for all kind of (fixed) vb criteria
c...   critor .. criterium to determine orthogonality
c...   crilow .. criterium  for orthogonalization (e.g.jacoby in lowdin)
c...   cridep .. criterium for dependency (on c**2)
c...   symcri .. criterium for symmetry equivalence
c...   remcri .. criterium for including excitation on basis of 1-e ints
c...   crihyb .. max hybrid contam. allowed in matrix and (*100) vector
c...   cribri .. min. size of bcoeff to avoid annihilation
c...   crinor .. min. norm of a vector
c...   cribub .. a meaningful difference (in a bubble sort)
c...   criign .. size of matrix elements that really may be ignored
c...   cripop .. criterium for pipek/mizek localization
c...   normally critor = 1.0d-14, crilow = critor/10, cridep = critor*10
c...            symcri = critor*100 remcri = 1.0d-10 crihyb = 1.0d-14
c...            cribri = 1.0d-8 crinor = critor*100  cribub = 1.0d-10
c...            criign = 1.0d-44 cripop 1.0d-10
      real*8 critor,crilow,cridep,symcri,remcri,crihyb,cribri,
     1     crinor,cribub,criign,cripop
      common/vbcri/ critor,crilow,cridep,symcri,remcri,crihyb,cribri,
     1              crinor,cribub,criign,cripop
c
c...  frozen core definition for TURTLE
c
      real*8 dxyz
      integer ncore,mapcie,mapiee,ncore_i,nsplice
      parameter (nsplice=512)
      common/splice/dxyz(3),ncore,ncore_i,
     1              mapcie(nsplice),mapiee(nsplice)
C...
c...  atonam for names of atoms (for mulliken population analysis,
c...  and automatic atomic (strictly localised) scf-procedure)
c...  atom is used her to denote fragments
c...  nopa contains the # orbitals per atom, iopa contains them
c...  mullik indicates if the analysis is requested, hybry indicates
c...  wether the scf-procedure just should used atom-centered ao's
c...  per atomic scf orbital. iacat contains the active (scf) orbital
c...  per atom, nacat contains the number of the per atom.
c...  atomao : fragment numbers / ao (may be more)
c...           then stored as frag1, maxato*frag2, maxato**2*frag3 etc.
c     ifrzat : hybrid handling
c              0 : normal, 1 frozen, 2 super
c...
      character*8 atoms
      common /atonam/ atoms(maxato)
      integer atomao,nopa,iopa,natom,iacat,nacat,ifrzat
      common /infato/ nopa(maxato),iopa(maxopa,maxato),natom,
     &                iacat(maxact,maxato),nacat(maxato),
     &                atomao(mxorbvb),ifrzat(maxato),
     &                mullik,hybry,clean,hguess,guesao
      logical mullik,hybry,clean,hguess,guesao
c     common for equivalences in VB
c     eqvmo/eqvao list of equivalent MO's and AO's in that order.
      integer eqvmo,eqvao,neqmo,neqao,meqmo,meqao
      logical equiv,equivir
      parameter (meqmo=10,meqao=100)
      common/vbequiv/eqvmo(meqmo),eqvao(meqao,meqmo),neqmo,neqao,
     1               equiv,equivir
*
*     to signal which scf procedure to be used, e.g., super-ci or Newton-Raphson
*     default is super-ci
*
      real*8 scalhs, sci2nr, delta, crmax
      integer nits, nitsci, nege, nneg
      logical oqcscf, ocong, ospci, oinve, ograd, ocorr, oauto,
     +        oiter, oqcof
      common / quadratic / nitsci,  nits,  nege, nneg, scalhs, sci2nr,
     +                     delta, crmax, oqcscf, ocong, ospci, oinve, 
     +                     ograd, ocorr, oauto, oiter, oqcof
*
      dimension v(nb,*),bi(*),q(lword),vnew(nb,nact)
*
      character*4 lafz
      character*10 charwall
*
      dlog10(aa) = dlog(aa) / dlog(10.d0)
*
*     check convergence
*
      scfconv = .false.
      lafz = 'on'
      if (oqcof) lafz = 'off'
*
      div = 0.0d0
      do 10 i = 2, nstruc
10    div = dmax1(dabs(bi(i)),div)
*
      if (optcri.eq.2) then
*
*        use overlap of current wavefunction with previous one as crit.
*
         if (swave.ge.criscf) scfconv = .true.
*
      else if (optcri.eq.1) then

*        use largest change in mo coefficient as criterion
*
         if (div.le.criscf) scfconv = .true.
*
      else
*
*        use real brillouin theorem
*
         if (brmax.le.criscf) scfconv = .true.
      end if
*
      if (nitscf.ge.maxscf) then
         scfconv = .true.
         write(iwr,20) nitscf
20       format(/,' iterations ended due to excessive iterations, viz.',
     +           i4)
      end if
*
      del = dabs(1.d0-swave)
*
*     1.0d-14 is machine precision .. we may take critor
*
      if (del.gt.critor) then
         del = dlog10(del) - .5d0
         idel=idint(del)
      else
         idel = -15
      end if
*
c      write(iwr,'(a,f35.28)') 'delta E ',delta
c      if (dabs(delta).lt.1.0d-10) scfconv = .true.

      write(iwr,30) nitscf,cpulft(1),charwall(),evb,delta,div,idel,
     +              brmax,lafz
30    format(/,' it.',i5,' at',f9.2,a10,' evb',f17.10,' del-e',f17.10,
     +            ' div ',1pe7.1,'/del',i3,'/brm ',1pe7.1,' qc ',a3)
*
*     build transformation matrix
*
      call vclr(q,1,nsa*nact+nsa*nact)
*
      iq = 1
      jq = 1
      iibi = 1
      maxoc = 0
      do i = 1, nscf
         ifrom = iex(1,i)
*
*        ifrom = 0 means that that orbital is frozen
*
         if (ifrom.ne.0) then
            maxoc = max(maxoc,ifrom)
            q( (ifrom-1)*nsa + ifrom ) = bi(1)
            do j = 2, nex(i) + 1
               iibi = iibi + 1
               ito = iex(j,i)
               q( (ifrom-1)*nsa + ito ) = bi(iibi) * iflip(ifrom,ito)
     +                                             * dampvb
            end do
            jq = jq + 1
            if (jq.le.iequi(iq)) then
               iibi = iibi - nex(i)
            else
               iq = iq + 1
               jq = 1
            end if
         end if
      end do
*
*     update orbitals
*
      if (equiv) call primequiv(v,nbasis,1,'occ','beg orbopt')
      if (iprins.gt.20) then
         write(iwr,60)
60       format(/,' orbitals before update, mos =>  ')
         call prvc(v(1,ncore+1),nact,nbasis,v,'o','l')
         write(iwr,70)
70       format(/,' transformation matrix for orbital update')
         call prvc14(q,maxoc,nsa,q,'o','a')
      end if
*
*     dump old vectors to tape33 for oldnew (determines overlap between
*     old and new wavefunctions) and natorb (1-electron density matrix)
*
      kvb7_vo = kscra7vb('kvb7_vo',nb*(ncore+nact),'r','w')
      kvb7_vn = kscra7vb('kvb7_vn',nb*(ncore+nact),'r','w')
      kvb7_vnew = kscra7vb('kvb7_vnew',nb*nact,'r','w')
      call wrt3(v,nb*(ncore+nact),kvb7_vo,num8)
*
*     if maxscf=0 we have done the work but do not change the orbitals
*
      if (maxscf.eq.0) then
         write(iwr,80)
80       format(' ****** Orbitals are unchanged ******')
         call fmove(v(1,ncore+1),vnew,nbasis*nact)
      else
       call mxmd(v(1,ncore+1),1,nbasis,q,1,nsa,vnew,1,nbasis,
     +           nbasis,nsa,nact)
      end if
*
*
      call fmove(vnew,v(1,ncore+1),nbasis*nact)
*
c      write(iwr,84)
c84    format(/,' orbitals after update, unnormalised =>  ')
c      call prvc(vnew,nact,nbasis,vnew,'o','l')

*
*     normalise orbitals and vnew as well / read s into q
*
*     diis (normalises q itself)
*
      if (idel.lt.iwdiis.and.nitscf.ge.ntdiis.and..not.scfconv) then
         ovbdiis = .true.
         kscr = 2*nsa*nact+1
         call vbdiis(q,v(1,ncore+1),vnew,nbasis,nact,
     +               q(kscr),lword-kscr,madiis,midiis)
*
*        the diis transformations may ruin hybrids, so re-clear
*
         if (clean.and.hybry) then
            call clvecvb(v(1,ncore+1),nact,nbasis,'active','small',
     +                   'after diis')
         end if
      else
         call fmove(vnew,v(1,ncore+1),nbasis*nact)
      end if

      call get1e(q,dummy,'s',q)
      call normt(v(1,ncore+1),nbasis,nact,q,q(nbasis*(nbasis+1)/2+1))
      call normt(vnew,nbasis,nact,q,q(nbasis*(nbasis+1)/2+1))
*
c      write(iwr,85)
c85    format(/,' orbitals after update, normalised =>  ')
c      call prvc(vnew,nact,nbasis,vnew,'o','l')
*
*     dump new vectors to tape33 for oldnew (determines overlap between
*     old and new wavefunctions) and natorb (1-electron density matrix)
*     info : orbitals updated but singles not orthogonal to doubles yet
*     (caused bug in overlap old/new and natural orbitals) now in vbscf
*     also the previous vnew vectors in case diis fails
*
      kvb7_vnew = kscra7vb('kvb7_vnew',nb*nact,'r','w')
      call wrt3(vnew,nb*nact,kvb7_vnew,num8)
*
*     make sure equivalence restrictions are obeyed
*
      if (equiv) call primequiv(v,nbasis,1,'occ','end orbopt')
*
      kvb7_vn = kscra7vb('kvb7_vn',nb*(ncore+nact),'r','w')
      call wrt3(v,nb*(ncore+nact),kvb7_vn,num8)
*
      if (iprins.ge.10) then
         write(iwr,90) (bi(i),i=1,nstruc)
90       format(/,' brillouin-state vector:',/,(10e12.4))
         write(iwr,100)
100      format(//,' --- new vb orbitals ---')
         call prvc(v(1,ncore+1),nact,nbasis,v,'o','l')
      end if
*
      return
      end

************************************************************************
      subroutine vbexcite(ndet,coeff,idet,bcoeff,ibdet,nelec,ifrom,itoj,
     +                    nb,nalfa)
************************************************************************
*
      implicit real*8 (a-h,o-z), integer (i-n)
*
*     create all excited determinants by replacing orbital 'ifr' with
*     orbital 'ito' in the reference determinants
*     ndet  = number of reference determinants
*     coeff = coefficients of reference determinants (input)
*     idet  = reference determinants                 (input) 
*     nb    = number of excited determinants         (output)
*     bcoeff= coefficients of excited determinants   (output)
*     ibdet = excited determinants                   (output)
*
      real*8  shiscf,criscf,evb,ebi,swave,percmax,perccur,brmax,
     &      dampvb,fixb0
      integer optcri,iprins,nitscf,maxscf,isecv,nitdiis,maxdiis,ntscf,
     &        ipcmt,maxscf_s
      logical perall,scfconv,p0core,unitary
      common /scftvb/dampvb,fixb0,shiscf,criscf,evb,ebi,swave,percmax,
     &               perccur,brmax,
     &               optcri,iprins,nitscf,maxscf,isecv,nitdiis,
     &               maxdiis,ntscf,ipcmt,maxscf_s,perall,scfconv,
     &               p0core,unitary
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
*
      dimension idet(nelec,*), ibdet(nelec,*), bcoeff(*), coeff(*)
*
*     we don't want to change 'ifrom' and 'ito' if 'unitary' is
*     used, so copy

      ifr = ifrom
      ito = itoj
      nb = 0
      it = 0
      ione = 1
      isign = + 1
*
10    continue
*
      do k = 1, ndet
          call bcreat(idet(1,k),nelec,nalfa,ibdet(1,nb+1),nbb,ifr,ito)
*
*         nbb = # new brillouin dets. which results by replacing 
*         orbital 'ifrom' by orbital 'itoj' in the 
*         referece det. idet(1,k)
*
          l = nb
20        if (l.lt.nb+nbb) then
             l = l + 1
*
*            check if new determinant is already present in ibdet(*)
*
             if (unitary) then
                ipos = 1
             else
                ipos = iflip(ifrom,itoj)
             end if
*
             do 30 m = 1, nb
                isignq = isame(ibdet(1,m),ibdet(1,l),
     +                   ibdet(1,nb+nbb+1),nelec,nalfa)
                if (isignq.ne.0) then
*
*                  we already had this determinant
*                  so add this one (l-th det.) to the 
*                  one present before
*                  
                   bcoeff(m) = bcoeff(m) + isignq * coeff(k)  
     +                         * ipos * isign
*
                   l = l - 1
                   nbb = nbb - 1
                   if (l.lt.nb+nbb) then
*                  
*                     copy the next one onto l-th det. and check again
*
                      do 40 n = 1, nelec
40                    ibdet(n,l+1) = ibdet(n,l+2)
                   end if
                   go to 20
               end if
30           continue
             it = it + 1
             bcoeff(it) = coeff(k) * ipos * isign
             go to 20
          end if
          nb = nb + nbb
      end do 
*
      if (unitary.and.ione.eq.1) then
         ione = 2
         ifr = itoj
         ito = ifrom
         isign = - 1
         go to 10
      end if
c      write(iwr,'(a,i3,a,i3)') 'after unitary ',ifr,' -> ',ito
c      call printpsi(nb-ipt,bcoeff(ipt+1),ibdet(1,ipt+1),nelec)
*
      return
      end
*
************************************************************************
      subroutine fillhes(hes,ndim,mbril,cih,nstruc,orbci)
************************************************************************
*
      implicit real*8 (a-h,o-z), integer (i-n)
*
*     fill the square matrix hes(*,*) symmetrically with orb(*) and 
*     cih(*) triangles and orbci(*) rectangular array
*
      dimension hes(ndim,ndim), cih(*), orbci(nstruc,mbril)
*
      ind(i,j) = max0(i,j)*(max0(i,j)-1)/2+min0(i,j)
*
      do i = 1, nstruc
         do j = 1, i
            hes(i+mbril,j+mbril) = cih(ind(i,j))
         end do
      end do
      do i = 1, nstruc
         do j = 1, mbril
            hes(i+mbril,j) = orbci(i,j)
         end do
      end do
      do i = 1, ndim
         do j = i+1, ndim
            hes(i,j) = hes(j,i)
         end do
      end do
*
      return
      end
*
************************************************************************
      subroutine getpsi0(ci,nstruc,igroup,ngroup,cdet,coeff,ndetps,
     +                   idetps,idet,idetall,ibdet,pacdet,nelec,nalfa,
     +                   ldet,ndtot,ncoeff,bcoeff,nword)
************************************************************************
*
      implicit real*8 (a-h,o-z),integer (i-n)
*
*     make psi0 on determinant basis 
*
*     ci(*)   = ci coefficients                                 (input)
*
*     idet    = all unique determinant in psi0                  (output)
*     cdet    = coefficients of all unique determinant          (output)
*     ldet    = number of unique determinants in psi0           (output)
*
*     idetall = complete set of determinants for all structures (output)
*     coeff   = coefficients of determinants for all structures (output)
*     ndtot   = total number of determinants for all structures (output)
*    
*     note: 
*          on entry cdet(*) contains the coefficients of all
*          determinants as input
*     
*
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
c... common to  sitmch from 8 to 16 bits orb ital packing
c...  n8_16 is # bits used for packing; i.e. 8 (255 orbitals) 16 (32000)
c...  n16_32 is the correspnding one for 16 vs 32 bits packin 
c...  GAMESS  parameters cannot be used as VB may have more MO's than AO's
c...  parameters are set in vbcrestr
c
      integer n8_16,n16_32,n340
      common/c8_16vb/n8_16,n16_32,n340
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
*
      dimension ci(*),igroup(5,*),coeff(*),cdet(*),idetps(*),
     +          idet(nelec,*),idetall(nelec,*),
     +          pacdet(*),bcoeff(*),ndetps(*),ibdet(nelec,*)
*
*
      id = 1
      ip = 1
      it = 1
      ik = 0
      ldet = 0
      xx = 0.0d0
*
      call vclr(coeff,1,ncoeff)
      call izero(ncoeff*nelec,idetall,1)
*
      do 30 i = 1, ngroup
*
         call vclr(bcoeff(ldet+1),1,igroup(1,i))
         call izero(igroup(1,i)*nelec,idet(1,ldet+1),1)
         call unpack(pacdet(ip),n8_16,idet(1,ldet+1),igroup(1,i)*nelec)
*
         ip = ip + (igroup(1,i)*nelec-1) / (64/n8_16) + 1
*
*        loop over structures in configuration
*
         do 20 j = 1, igroup(2,i)

            xx = ci(it)
*
*           get combined coefficient of determinants in structure
*           i.e.  spin-coefficient * ci-coefficient
*
            do 10 k = 1, ndetps(it)
               ibb = idetps(id)
               ik = ik + 1
               do l = 1, nelec
                  idetall(l,ik) = idet(l,ibb+ldet)
               end do
               coeff(ik) = cdet(id)
               bcoeff(ibb+ldet) = bcoeff(ibb+ldet) + cdet(id) * xx
10          id = id + 1
20       it = it + 1
         ldet = ldet + igroup(1,i)
30    continue
*
      ndtot = ik
*
      if (nword.lt.ldet+nelec*ldet/2+1) then
         call corfait(2*ldet,nword,'need more in psi0det/atmol')
      end if
*
      call checkpsi0(ldet,bcoeff,idet,bcoeff(ldet+1),nelec,nalfa,
     +               lldet,cdet,ibdet)
      do i = 1, lldet
         bcoeff(i) = cdet(i)
         do j = 1, nelec
            idet(j,i) = ibdet(j,i)
         end do
      end do
*
      ldet = lldet
*
      igroup(1,1) = ldet
      igroup(2,1) = 0
      igroup(3,1) = ldet
*
      call abmin(ci,it-1,rmin,imin)
      small = rmin
      call pack(pacdet,n8_16,idet(1,1),ldet*nelec)
      if (ldet.ne.igroup(1,1)) call caserr(' ldet ne igroup')
      kl = (ldet*nelec-1)/(64/n8_16)+1
      nl = lensec(3)+lensec(kl)+lensec(igroup(1,1))
      call secput(isect(79),79,nl,ib)
      igroup(2,1) = nelec
      call wrt3(igroup(1,1),3,ib,idaf)
      igroup(2,1) = 0
      ib = ib+1
      call wrt3(bcoeff,igroup(1,1),ib,idaf)
      ib = ib + lensec(igroup(1,1))
      call wrt3(pacdet,kl,ib,idaf)
*
      return
      end
*
************************************************************************
      subroutine getphik(kstruc,indetps,coeff,idets,nkdet,ikdet,ckdet,
     +                   nelec)
************************************************************************
*
      implicit real*8 (a-h,o-z) , integer (i-n)
*
*     copy the determinats of kstruc-th structure in ikdet(*) and
*     their coefficients in ckdet(*)
* 
      dimension idets(nelec,*), coeff(*), ikdet(nelec,*), ckdet(*)
      dimension indetps(*)
*     
      nkdet = indetps(kstruc)
      nskip = 0
      do i = 1, kstruc - 1
         nskip = nskip + indetps(i)
      end do
      do i = 1, indetps(kstruc)
         ckdet(i) = coeff(nskip+i)
         do j = 1, nelec
            ikdet(j,i) = idets(j,nskip+i)
         end do 
      end do
*
      return
      end
*
*************************************************************************
      subroutine printpsi(n,coeff,idet,nelec)
************************************************************************
*
      implicit real*8 (a-h,o-z) , integer (i-n)
*
*     print the determinant and their coefficients
*     upto 14th decimal places
*
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
*
      dimension coeff(*),idet(nelec,*)
*
      do i = 1, n
         write(iwr,10) i, coeff(i), (idet(j,i),j = 1, nelec)
10       format(1x,'det. ',i3,' coeff. ',f18.14,' : ', (t40,30i4))
      end do
*
      return
      end 
*
************************************************************************
      subroutine tripr14(a,ndim)
************************************************************************
*
      implicit real*8 (a-h,o-z) , integer (i-n)
*
*     print (a maximum of 5 columns per row) the triangular matrix a
*     upto 14th decimal places
*
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
*     
      dimension a(*)
*     
      write(iwr,'(/)')
      kend = ndim / 5 + 1
      jstart = 1
      it = 0
      is = 0
      nextra = 0
      do 30 k = 1, kend
         n = 0
         do 20 j = jstart, ndim
            n = n + 1
            if (n.gt.5) then
               iend = 5
            else
               iend = n
               is = is + n + nextra
            end if
            write(iwr,10) j,(a(it+i),i=1,iend)
10          format(2x,i3,2x,5f20.14)
            it = it + n + nextra
20       continue
         write(iwr,'(/)')
         nextra = nextra + 5
         it = is + nextra
         jstart = jstart + 5
30    continue
*
      return
      end
*
************************************************************************
      subroutine prsq14(v,m,n,ndim)
************************************************************************
*
      implicit real*8 (a-h, p-w), integer (i-n), logical (o)
*
*     print out a square matrix (maximum 5 columns per row)
*     upto 14th decimal places
*     m = columns, n = rows
*
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
*
      dimension v(ndim,*)
*
      max = 5
      imax = 0
10    imin = imax + 1
      imax = imax + max
      if (imax .gt. m) imax = m
      write (iwr,20)
      write (iwr,30) (i,i = imin,imax)
      write (iwr,20)
      do j = 1, n
         write (iwr,40) j,(v(j,i),i = imin,imax) 
      end do
      if (imax .lt. m) go to 10 
      write(iwr,20)
      return
20    format(1x)
30    format(6x,5(8x,i4,8x))
40    format(i5,1x,5f20.14)
*
      end
*
************************************************************************
      subroutine prvc14(c,m,n,eps,ch,test)
************************************************************************
*     
      implicit real*8 (a-h, o-z), integer (i-n)
*
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
*
*     routine for printing of column-vectors and eigenvalues.
*     now heavily changed with respect to original (esp. parameters)
*     # columns fixed to 5 // 14 decimal places
*   
*            c  : matrix of coefficients
*            m  : x columns to be printed
*            n  : dimension of c(n,m) and eps(n)
*            eps: vector of eigenvalues/occupations
*            ch : 'value' or 'v' => values; 'order' or 'o' => columns
*           test: 'label' ior 'l' => gamess labels if possible
*           
*
      dimension c(n,m), eps(m)
      character*(*) test,ch
      dimension sub(2)
      character*5 sub
      data sub /'-----','-----'/
*     
      if (ch.eq.'value'.or.ch.eq.'v') then
         if (test.eq.'label'.or.test.eq.'l') then
            call prev(c,eps,m,n,n)
         else  
            ich = 1
            go to 10
         end if
      else if (ch.eq.'order'.or.ch.eq.'o') then
         if (test.eq.'label'.or.test.eq.'l') then
            call prsql14(c,m,n,n)
         else
            ich = 0
            go to 10
         end if
      else
         call vberr(' wrong prvc call')
      endif

      return
*
10    ncc = 5
      nbl = (m-1) / ncc
      nlast = m - ncc * nbl
      nbl = nbl + 1
      nc1 = 0
      do 20 ib = 1, nbl
         if ( ib .eq. nbl ) ncc = nlast
         nc0 = nc1 + 1
         nc1 = nc1 + ncc
         if (ich.ne.0) write(iwr,30) ( eps(ic), ic=nc0,nc1 )
         if (ich.eq.0) write(iwr,40) ( ic, ic=nc0,nc1 )
30       format(/,7x,10f12.6)
40       format(/,7x,5(8x,i4,8x))
         write(iwr,50) ( sub, i = 1, ncc )
50       format(7x,5(6x,a5,a5,5x))
         write(iwr,60)
60       format(1h )
         do 70 ia = 1, n
            write(iwr,80) ia, ( c(ia,ic), ic = nc0, nc1 )
80          format(2x,i3,2x,5f20.14)
70       continue
20    continue
      write(iwr,30)
*
      return
      end
*
************************************************************************
      subroutine prsql14(v,m,n,ndim)
************************************************************************
*
      implicit real*8 (a-h, p-w), integer (i-n), logical (o)
      implicit character*8 (z), character*1 (x)
      implicit character*4 (y)
*
*     print out a square matrix with labels (maximum 5 columns per row)
*     upto 14th decimal places
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
*
      dimension v(ndim,*)
*
      max = 5
      imax = 0
10    imin = imax+1  
      imax = imax+max
      if (imax .gt. m) imax = m
      write (iwr,20)
      write (iwr,30) (i,i = imin,imax)
      write (iwr,20)
      do j = 1,n
         write (iwr,40) j,zbflab(j),(v(j,i),i = imin,imax)
      end do
      if (imax .lt. m) go to 10
      write (iwr,20)
      return
20    format(1x)
30    format(17x,5(8x,i4,8x))
40    format(i5,2x,a10,5f20.14)
      end
*
************************************************************************
      subroutine triprsq14(a,n)
************************************************************************
*
      implicit real*8 (a-h, p-w), integer (i-n), logical (o)
      implicit character*8 (z), character*1 (x)
      implicit character*4 (y)
*
*     print out a triangular matrix in square form 
*     a maximum of 5 columns per row
*     upto 14th decimal places
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
*
      dimension a(*), dd(5)
*
      mmax = 5
      imax = 0
10    imin = imax + 1
      imax = imax + mmax
      if (imax .gt. n) imax = n
      write (iwr,20)
      write (iwr,30) (i,i = imin,imax)
      write (iwr,20)
      do 40 j = 1, n
         k = 0
         do 50 i = imin, imax
            k = k + 1
            m = max(i,j) * (max(i,j)-1)/2 + min(i,j)
50       dd(k) = a(m)
         write (iwr,60) j, (dd(i),i = 1,k)
40    continue
      if (imax .lt. n) go to 10
      write (iwr,20)
      return
20    format(/)
30    format(6x,5(8x,i4,8x))
60    format(i5,1x,5f20.14)
*
      end
*
************************************************************************
