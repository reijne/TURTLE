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
      implicit REAL (a-h,o-z), integer (i-n)
*
*     control routine for orbital optimisation
*     using newton-raphson scheme or super-ci
*     adapted from 'subroutine vbscf'
*     v  : orbitals
*     ci : ci vector 
*
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/iofile)
INCLUDE(common/turtleparam)
INCLUDE(common/vbqc)

INCLUDE(common/vbdiist)
INCLUDE(common/tractlt)
      common /bypass/ index4,ihmat,idavid
INCLUDE(common/hsinfo)
INCLUDE(common/vbcri)
*
*     for dynamic memory allocation debugging purpose, contains 
*     IGMEM_ vars.
*
INCLUDE(../m4/common/gmemdata)
INCLUDE(../m4/common/gmempara) 
*
INCLUDE(common/brill)
*
INCLUDE(common/scftvb)
INCLUDE(common/ffile)
INCLUDE(common/splice)
*
      common /davcom/ eshift,thresj,thresd,maxcyc,maxdav,maxsel,iprind,
     +                alter,firsts
      logical alter,firsts
*
INCLUDE(common/infato)
INCLUDE(common/vbequiv)
*
INCLUDE(common/vbpert)
*
INCLUDE(../m4/common/mapper)
INCLUDE(../m4/common/restri)
INCLUDE(../m4/common/timeperiods)
*
INCLUDE(common/twice)
INCLUDE(common/first_vb)
INCLUDE(common/hsmattype)
*
INCLUDE(common/aivb)
INCLUDE(common/vbtimess)
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
_IF(debug)
      write(iwr,*) '(1) vbqcscf'
_ENDIF
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
_IF(debug)
      write(iwr,*) '(2) vbqcscf'
_ENDIF
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
_IF(debug)
      write(iwr,*) '(4) vbqcscf'
_ENDIF
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
_IF(debug)
        write(iwr,*) '(5) vbqcscf - reduce'
_ENDIF
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
_IF(debug)
        write(iwr,*) '(55) vbqcscf - reduce'
_ENDIF
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
_IF(debug)
        write(iwr,*) '(555) vbqcscf - reduce'
_ENDIF


      end if

_IF(debug)
      write(iwr,*) '(6) vbqcscf'
_ENDIF
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
_IF(debug)
        write(iwr,*) '(7) vbqcscf'
_ENDIF
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

_IF(debug)
        write(iwr,*) '(8) vbqcscf'
_ENDIF

      if (hybry) call clvecvb(v(ncore*nbasis+1),nactiv+nvirt,nbasis,
     +                        'check','null','in vbqcscf')
_IF(parallel)
*
*     to be sure broadcast orbitals en info; to get all noses aligned
*
      call pg_brdcst(8,v,8*nbasis*(ncore+nactiv+nvirt),0)
      nn = maxact*(1+maxex+1+maxex+maxex+1)+1+mxorbvb+1+
     +     (maxact+1)*maxact+maxact+7+maxact
      nn = nn*8/nipw()
      call pg_brdcst(9,nex,nn,0)
_ENDIF
*
*     save orbitals to disk as these are the vbqcscf orbitals
*
      call putqvb(v,nbasis,nactiv+ncore)
*
      call flushbuffer
*
_IF(peigss)
      call cleargas()
_ENDIF
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
_IF(debug)
      write(iwr,*) '(9) vbqcscf'
_ENDIF
*
*     (first integral = (00/00) = 0.0 => + 1)
*     don't use to much core in tran as virtual io may interfere
*     with normal io; If the user wants to commit suicide, let him
*     (jvl,2002)
*
      call transformvb(q(ks),q(kh),v)
*
_IF(debug)
      write(iwr,*) '(10) vbqcscf'
_ENDIF

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
_IF(debug)
            write(iwr,*) '(11) vbqcscf'
_ENDIF
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
_IF(debug)
      write(iwr,*) '(14) vbqcscf'
_ENDIF
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
_IF(debug)
      write(iwr,*) '(15) vbqcscf'
_ENDIF
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

_IF(debug)
      write(iwr,*) '(16) vbqcscf'
_ENDIF
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

_IF(debug)
      write(iwr,*) '(17) vbqcscf'
_ENDIF

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

_IF(debug)
      write(iwr,*) '(18) vbqcscf'
_ENDIF

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
_IF(debug)
        write(iwr,*) '(19) vbqcscf'
_ENDIF
c
c...  generate psi0 on determinant basis
c
      call psi0det(ci,nstruc,q(kigro),ngroup,q(kcoef),ncoeff,
     +             ndets,q(kndet),q(kidps),q(kidet),q(kjdet),
     +             q(kpacd),nelec,nalfa,ndettot,q(kscr),kmemscr,
     +             'dontsaveondisk')
_IF(debug)
        write(iwr,*) '(20) vbqcscf'
_ENDIF
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
_IF(debug)
        write(iwr,*) '(21) vbqcscf'
_ENDIF
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
_IF(debug)
        write(iwr,*) '(KHATAM) vbqcscf'
_ENDIF
      return
      end

*
************************************************************************
      subroutine vbquad(supers,superh,superg,hes,grad,ci,nbasis,nsa,
     +                  nvirt,q)
************************************************************************
*
      implicit REAL (a-h,o-z), integer (i-n)
*
*     dynamical core allocation and brillouin control routine
*     for dynamic memory allocation debugging purpose, contains i
*     IGMEM_ vars.
*
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/gmemdata)
INCLUDE(../m4/common/gmempara) 
INCLUDE(common/turtleparam)
INCLUDE(common/c8_16vb)
INCLUDE(common/brill)
INCLUDE(common/hsinfo)
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
_IF(atmol)
      rewind 25
      read(25) nelec,nalfa,nconf,ndets,nstruc,maxdet,maxstr,
     +         nwpack,ncoeff,imax
_ELSE
      call wr15vb(nelec,nalfa,nconf,ndets,nstruc,maxdet,maxstr,
     +            nwpack,ncoeff,imax,ndum11,ndum12,'read')
_ENDIF
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
      implicit REAL (a-h,o-z), integer (i-n)
*
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/gmemdata)
INCLUDE(../m4/common/gmempara)
INCLUDE(../m4/common/restri)
INCLUDE(../m4/common/iofile)
INCLUDE(../m4/common/mapper)
INCLUDE(common/c8_16vb)
INCLUDE(common/tractlt)
INCLUDE(common/ffile)
INCLUDE(common/hsinfo)
INCLUDE(common/scftvb)
INCLUDE(common/turtleparam)
INCLUDE(common/vbqc)
INCLUDE(common/vbproper)
INCLUDE(common/brill)
INCLUDE(common/twice)
INCLUDE(common/infato)
INCLUDE(common/splice)
INCLUDE(common/vbcri)
INCLUDE(common/vbdiist)
*
_IF(parallel)
INCLUDE(common/parinf)
_ENDIF
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
_IF(parallel)
      call pg_dgop(1111,dd,1,'+')
      call pg_dgop(1112,ds,1,'+')
_ENDIF
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
      implicit REAL (a-h,o-z), integer (i-n)
*
*     generate orbital hessian matrix
*
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/gmemdata)
INCLUDE(../m4/common/gmempara)
INCLUDE(../m4/common/iofile)
INCLUDE(common/turtleparam)
INCLUDE(common/brill)
INCLUDE(common/scftvb)
INCLUDE(common/twice)
INCLUDE(common/vbqc)
*
_IF(parallel)
INCLUDE(common/parinf)
_ENDIF
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
_IF(parallel)
      call pg_dgop(1113,grad,mbril,'+')
      call pg_dgop(1114,orbh,mbril,'+')
      call pg_dgop(1115,bhamil(2),nbdim,'+')
      call pg_dgop(1116,boverl(2),nbdim,'+')
_ENDIF
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
_IF(parallel)
      call pg_dgop(1118,hes,nd*nd,'+')
_ENDIF
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
      implicit REAL (a-h,o-z), integer (i-n)
*
*     generate ci hessian matrix
*
*
INCLUDE(../m4/common/sizes) 
INCLUDE(../m4/common/gmemdata)
INCLUDE(../m4/common/gmempara)
INCLUDE(../m4/common/iofile)
INCLUDE(common/turtleparam)
INCLUDE(common/brill)
INCLUDE(common/hsinfo)
*
_IF(parallel)
INCLUDE(common/parinf)
_ENDIF
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
_IF(parallel)
      call pg_dgop(1119,cis2,nstruc,'+')
_ENDIF
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
      implicit REAL (a-h,o-z), integer (i-n)
*
*     generate orbital-ci mixed hessian 
*
INCLUDE(../m4/common/sizes) 
INCLUDE(../m4/common/gmemdata)
INCLUDE(../m4/common/gmempara)
INCLUDE(../m4/common/iofile)
INCLUDE(common/turtleparam)
INCLUDE(common/brill)
*
_IF(parallel)
INCLUDE(common/parinf)
_ENDIF
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
_IF(parallel)
      call pg_dgop(1120,ciorbh,nstruc*mbril,'+')
_ENDIF
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
      implicit REAL (a-h,o-z), integer (i-n)
*              
*     prepares to invoke the conjugate gradient method to solve Ax = b
*    
*        
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/gmemdata)
INCLUDE(../m4/common/gmempara)
INCLUDE(../m4/common/iofile)
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
      implicit REAL (a-h,o-z), integer (i-n)
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
      implicit REAL (a-h,o-z), integer (i-n)
*
INCLUDE(../m4/common/iofile)
INCLUDE(common/vbqc)
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
      implicit REAL (a-h,o-z), integer (i-n)
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
      implicit REAL (a-h,o-z), integer (i-n)
*
*     use info in bi-vector to update orbitals
*     if converged return scfconv = .true.
*     v are all the orbitals (including core)
*     bi update-vector
*     q scratch space
*     nstruc (in hsinfo) now have the dimension of bi(*)
*
INCLUDE(common/vbdiist)
INCLUDE(common/tractlt)
INCLUDE(common/hsinfo)
INCLUDE(common/turtleparam)
INCLUDE(../m4/common/sizes)
INCLUDE(common/vbpert)
INCLUDE(../m4/common/iofile)
INCLUDE(common/brill)
INCLUDE(common/scftvb)
INCLUDE(common/vbcri)
INCLUDE(common/splice)
INCLUDE(common/infato)
INCLUDE(common/vbequiv)
INCLUDE(common/vbqc)
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
_IF(atmol)
      write(33) nb,ncore+nact,((v(ii,jj),ii=1,nb),jj=1,ncore+nact)
_ELSE
      kvb7_vo = kscra7vb('kvb7_vo',nb*(ncore+nact),'r','w')
      kvb7_vn = kscra7vb('kvb7_vn',nb*(ncore+nact),'r','w')
      kvb7_vnew = kscra7vb('kvb7_vnew',nb*nact,'r','w')
      call wrt3(v,nb*(ncore+nact),kvb7_vo,num8)
_ENDIF
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
_IF(atmol)
      call vberr(' vnew not properly  dumped')
      write(33) nb,ncore+nact,((v(ii,jj),ii=1,nb),jj=1,ncore+nact)
_ELSE
      kvb7_vn = kscra7vb('kvb7_vn',nb*(ncore+nact),'r','w')
      call wrt3(v,nb*(ncore+nact),kvb7_vn,num8)
_ENDIF
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
      implicit REAL (a-h,o-z), integer (i-n)
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
INCLUDE(common/scftvb)
INCLUDE(../m4/common/iofile)
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
      implicit REAL (a-h,o-z), integer (i-n)
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
      implicit REAL (a-h,o-z),integer (i-n)
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
INCLUDE(../m4/common/iofile)
INCLUDE(common/c8_16vb)
INCLUDE(../m4/common/restri)
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
      implicit REAL (a-h,o-z) , integer (i-n)
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
      implicit REAL (a-h,o-z) , integer (i-n)
*
*     print the determinant and their coefficients
*     upto 14th decimal places
*
INCLUDE(../m4/common/iofile)
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
      implicit REAL (a-h,o-z) , integer (i-n)
*
*     print (a maximum of 5 columns per row) the triangular matrix a
*     upto 14th decimal places
*
INCLUDE(../m4/common/iofile)
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
      implicit REAL (a-h, p-w), integer (i-n), logical (o)
*
*     print out a square matrix (maximum 5 columns per row)
*     upto 14th decimal places
*     m = columns, n = rows
*
INCLUDE(../m4/common/iofile)
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
      implicit REAL (a-h, o-z), integer (i-n)
*
_IFN(atmol)
INCLUDE(../m4/common/iofile)
_ENDIF
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
_IFN(atmol) 
            call prev(c,eps,m,n,n)
         else  
_ENDIF   
            ich = 1
            go to 10
         end if
      else if (ch.eq.'order'.or.ch.eq.'o') then
         if (test.eq.'label'.or.test.eq.'l') then
_IFN(atmol)
            call prsql14(c,m,n,n)
         else
_ENDIF
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
      implicit REAL (a-h, p-w), integer (i-n), logical (o)
      implicit character*8 (z), character*1 (x)
      implicit character*4 (y)
*
*     print out a square matrix with labels (maximum 5 columns per row)
*     upto 14th decimal places
*     
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/prints)
INCLUDE(../m4/common/iofile)
INCLUDE(../m4/common/runlab)
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
      implicit REAL (a-h, p-w), integer (i-n), logical (o)
      implicit character*8 (z), character*1 (x)
      implicit character*4 (y)
*
*     print out a triangular matrix in square form 
*     a maximum of 5 columns per row
*     upto 14th decimal places
*
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/prints)
INCLUDE(../m4/common/iofile)
INCLUDE(../m4/common/runlab)
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
