c  $Author: jvl $
c  $Date: 2010-08-10 16:49:47 +0200 (Tue, 10 Aug 2010) $
c  $Locker:  $
c  $Revision: 6176 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/secd_parallel.m,v $
c  $State: Exp $
c
      subroutine scfdd(q)
c
c     driving routine for scf second derivative calculations
c
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
      dimension q(*)
c
      common/maxlen/maxq
INCLUDE(common/segm)
INCLUDE(common/vectrn)
INCLUDE(common/common)
INCLUDE(common/iofile)
INCLUDE(common/cndx40)
INCLUDE(common/cndx41)
INCLUDE(common/infoa)
INCLUDE(common/mapper)
INCLUDE(common/statis)
INCLUDE(common/symtry)
      common/scfopt/maxcyc,mconv,nconv
c
      common/junke/maxt,ires,ipass,nteff,
     1     npass1,npass2,lentri,nbuck,mloww,mhi,ntri,iacc,iontrn
c
INCLUDE(common/restrj)
      logical lsave
INCLUDE(common/prnprn)
INCLUDE(common/timeperiods)
INCLUDE(common/mp2grad_pointers)
INCLUDE(common/nshel)
INCLUDE(common/sortp)
c
c control common affects mp2 transformation
INCLUDE(common/mp2grd)
INCLUDE(common/secd_pointers)
c
      logical pg_create_inf
      logical pg_destroy_inf
      character*10 charwall
c
      character *8 fkd,blank,closed,oscf,grhf,dipd
      character *15 fnm
      character *5  snm
      data fkd,blank /'fockder','        '/
      data closed,oscf,grhf/'closed','oscf','grhf'/
      data dipd /'dipder'/
      data fnm/'secd_parallel.m'/
      data snm/'scfdd'/
c
c     evaluate integrals
c
      if (opass2) then
         if(nopk.eq.0.or.nopk.ne.nopkr.or.iofsym.eq.1.or.
     +      iofrst.ne.iofsym) then
           write (iwr,130)
           opass2 = .false.
         endif
      endif

      call start_time_period(TP_2D_TOTAL)

      nopk = 1
      iofsym = 0
      isecvv = isect(8)
      itypvv = 8
      nconv = max(nconv,7)
      mmaxq = maxq
      if (irest5.ge.6) then
         go to 30
      else
         if (irest5.lt.1) then
            call start_time_period(TP_2D_AOINTS)
c           call jobstep('integ')
            call integ(q)
            call end_time_period(TP_2D_AOINTS)
            irest5 = 1
            call revise
         end if
         if (irest5.lt.2) then
            call start_time_period(TP_2D_SCF)
c           call jobstep('scfrun')
            call scfrun(q)
            call end_time_period(TP_2D_SCF)
            irest5 = 2
            irest = 5
            call revise
         end if

         if (irest5.lt.3) then
c
c           test code for allocating all of the stashed data into memory
c           this should be performed in situ, and freed up when complete
c           but for the moment do it all up front
c
            if (omem2nd) then
               nat3 = nat * 3
               len1 = nat3 * nx
               if (omem2nd_ga) then
c dh
                  if (.not.pg_create_inf(0,len1,1,'dham-stash',0,0,
     &               iof2nd(1),fnm,snm,IGMEM_NORMAL)) then
                     call caserr('failed to allocate dham-stash')
                  endif
c ds
                  if (.not.pg_create_inf(0,len1,1,'doverlap-stash',0,0,
     &               iof2nd(2),fnm,snm,IGMEM_NORMAL)) then
                     call caserr('failed to allocate doverlap-stash')
                  endif
               else
c dh
                  iof2nd(1) = igmem_alloc_inf(len1,fnm,snm,
     &                     'dham-stash',IGMEM_NORMAL)
c ds
                  iof2nd(2) = igmem_alloc_inf(len1,fnm,snm,
     &                      'doverlap-stash',IGMEM_NORMAL)
               endif

c pdens (now in chfndr)

               len3 = nat3 * ikyp(ncoorb)

c soln (now in chfndr)

               len2 = nat3*(noccb*nvirta + (noccb-nocca)*nocca)

               itotal = 2 * len1 + len3 + max((len2+len2),len3)
               write(iwr,7040) itotal
            endif

            fkder = fkd
            call start_time_period(TP_2D_HFGRDN)
c           call jobstep('gradient')
            call hfgrdn(q)
            call end_time_period(TP_2D_HFGRDN)
            irest5 = 3
            call revise
         end if
         if (opass6) then
            if ((scftyp.eq.grhf .and. iscftp.lt.5) .or.
     +          (scftyp.ne.grhf .and. iscftp.lt.3)) call caserr(
     +          'less restricted 4-index transformation required ')
            go to 20
         end if
c
c  restrict the transformation
c
         if (scftyp.eq.grhf) then
            iscftp = max(iscftp,5)
         else
            iscftp = max(iscftp,3)
         end if
c
         npass1 = max(nps1,1)
         npass2 = max(nps2,1)
         iconvv = max(iconvv,9)
         lword4 = maxq
c
         call revise
      end if
c
c   do the 4-index transformation and set
c
 20   oprn(4) = .false.
c
c     itmp = nw196(5) + num*num + num + nx + 2 * (num * nt)
c     ibase2 = igmem_alloc_inf(itmp,fnm,snm,'core_indxsy',IGMEM_DEBUG)
      ibase2 = igmem_alloc_all_inf(itmp,fnm,snm,'core_indxsy',
     &                             IGMEM_DEBUG)
      call indxsy(q(ibase2))
      call gmem_free_inf(ibase2,fnm,snm,'core_indxsy')
c
      if (irest5.lt.4) then
         call start_time_period(TP_2D_INDX2T)
c        call jobstep('1el tran')
         write(iwr,7030) cpulft(1) ,charwall()
         itmp =  num*num + nx*4 + num + nw196(5)
         ibase2 = igmem_alloc_inf(itmp,fnm,snm,'core_indx2t',
     &                            IGMEM_DEBUG)
c
c   do the 2-index transformation
c
         call indx2t(q(ibase2))
         call gmem_free_inf(ibase2,fnm,snm,'core_indx2t')
         call end_time_period(TP_2D_INDX2T)
c 
         irest5 = 4
         call revise
         call timana(11)
      end if
      if (iscftp.lt.3 .or. (scftyp.eq.grhf .and. iscftp.lt.5))
     +    call caserr('less restricted 4-index transformation required')

c
c  load global arrays with transformed integrals
c  may need this when transformation is disabled
c
      call set41
      call start_time_period(TP_2D_MOINTS)
c     call jobstep('2el tran')
      opg_grad = .true.
      opg_alloc = .false.
c ...  schwarz inequality test
      mp2grad_schw = igmem_alloc_inf( nshell*(nshell+1)/2,fnm,snm,
     &               'mp2grad_schw',IGMEM_DEBUG )
c     switchable schwarz inequality 
      if (oschw) then
        if (nopk.ne.1) then
          i10 = igmem_alloc_inf(30000,fnm,snm,'gout',IGMEM_DEBUG)
        else
          i10 = igmem_alloc_inf(10000,fnm,snm,'gout',IGMEM_DEBUG)
        end if
        call coulmb(q(mp2grad_schw),q(i10))
        call gmem_free_inf(i10,fnm,snm,'gout')
      else
         call dcopy( nshell*(nshell+1)/2
     &        ,0.0d0, 0, q(mp2grad_schw), 1)
      end if
c ...  GA-based transformation
      call aprdmp2(q, emp2,mp2)
      call end_time_period(TP_2D_MOINTS)

      if (irest5.lt.5) then

c
         call start_time_period(TP_2D_TRNFKD)
         itmp = num*num+nx+nx+1
         ibase2 = igmem_alloc_inf(itmp,fnm,snm,'core_trnfkd',
     &                            IGMEM_DEBUG)
c        call jobstep('trnfkd')
         write(iwr,7000)cpulft(1) ,charwall()
         call trnfkd(q(ibase2))
         call gmem_free_inf(ibase2,fnm,snm,'core_trnfkd')
         call end_time_period(TP_2D_TRNFKD)
c
         call dksm_exp(q,q)
c
         call start_time_period(TP_2D_CHFNDR)
c        call jobstep('chfndr')
         call chfndr(q)
         call end_time_period(TP_2D_CHFNDR)
         fkder = blank
c
c     deletion now conducted in aprdmp2
c     call delete_ga(g_oooo)
c     call delete_ga(g_vvoo)
c     call delete_ga(g_vooo)
c     call delete_ga(g_vovo)
c
         call start_time_period(TP_2D_DMDER)
         itmp = num + 4* nx
         ibase2 = igmem_alloc_inf(itmp,fnm,snm,'core_dmder',IGMEM_DEBUG)
c        call jobstep('dmder')
         write(iwr,7010) cpulft(1) ,charwall()
         call dmder(q(ibase2))
         call gmem_free_inf(ibase2,fnm,snm,'core_dmder')
         call end_time_period(TP_2D_DMDER)
c
         call start_time_period(TP_2D_QMDER)
         itmp = nw196(5) + lenint(nshell*nt) + 1
         itmp = max(4 * nx, itmp)
         ibase2 = igmem_alloc_inf(itmp,fnm,snm,'core_qmderi',
     &                            IGMEM_DEBUG)
c        call jobstep('qmderi')
         call qmderi(q(ibase2))
         call gmem_free_inf(ibase2,fnm,snm,'core_qmderi')
         call end_time_period(TP_2D_QMDER)
c
         irest5 = 5
         irest = 7
         call revise
         if (runtyp.eq.dipd) go to 30
      end if
      if (irest5.lt.6) then
         lsave = oprn(5)
         oprn(5) = .false.

         call start_time_period(TP_2D_2D)
c        call jobstep('2nd D ints')
         write(iwr,7020) cpulft(1) ,charwall()
         call dertwo(q,q)
         call end_time_period(TP_2D_2D)

         oprn(5) = lsave
         irest5 = 6
         irest = 0
         call revise
      end if
c
c     ----- reset core allocation
c
 30   continue
      maxq = mmaxq
      if (omem2nd) then
        if (omem2nd_ga) then
c         free memory used to hold "stashed" data
c         fock 
          if (.not.pg_destroy_inf(iof2nd(4),'fock-stash',fnm,snm))
     &        call caserr('failed to destroy fock-stash')
c         pdens
          if (.not.pg_destroy_inf(iof2nd(5),'pdens-stash',fnm,snm))
     &        call caserr('failed to destroy pdens-stash')
        else
c         free memory used to hold "stashed" data
c         fock 
          call gmem_free_inf(iof2nd(4),fnm,snm,'fock-stash')
c         pdens
          call gmem_free_inf(iof2nd(5),fnm,snm,'pdens-stash')
        endif
      endif
c
      call gmem_free_inf(mp2grad_schw,fnm,snm,'mp2grad_schw')
c
      if (omem2nd) then
        if (omem2nd_ga) then
c         free memory used to hold "stashed" data
c         ds 
          if (.not.pg_destroy_inf(iof2nd(2),'doverlap-stash',fnm,snm))
     &        call caserr('failed to destroy doverlap-stash')
c         dh
          if (.not.pg_destroy_inf(iof2nd(1),'dham-stash',fnm,snm))
     &        call caserr('failed to destroy dham-stash')
        else
c         free memory used to hold "stashed" data
c         rhs
c            call gmem_free_inf(iof2nd(3),fnm,snm,'rhs-stash')
c         soln
c            call gmem_free_inf(iof2nd(6),fnm,snm,'soln-stash')
c         ds
             call gmem_free_inf(iof2nd(2),fnm,snm,'doverlap-stash')
c         dh
             call gmem_free_inf(iof2nd(1),fnm,snm,'dham-stash')
c
        endif
      endif
c
      call end_time_period(TP_2D_TOTAL)

      return
 130  format(/
     + 1x,'= = = = = = = = = = = = = = = = = = = = = = ='/
     + 1x,'= integrals must NOT be in supermatrix form ='/
     + 1x,'=        requested bypass is ignored        ='/
     + 1x,'= = = = = = = = = = = = = = = = = = = = = = ='/)
 7000 format(/1x,
     + 'construct transformed 1-electron matrix derivatives at ',
     +  f9.2,' seconds',a10,' wall')
 7010 format(/1x,
     + 'construct dipole and quadrupole derivatives at ',
     +  f9.2,' seconds',a10,' wall')
 7020 format(/1x,
     + 'commence analytic second derivatives integrals at ',
     +  f9.2,' seconds',a10,' wall')
 7030 format(/1x,
     + 'commence 2-index integral transformation at ', f9.2,
     + ' seconds',a10,' wall')
 7040 format(/1x,
     +'**** Node memory required for storing temporary matrices = ',
     +      i10,' Words'/)
      end

      subroutine trnfkd(q)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      logical out
INCLUDE(common/infoa)
      common/small/y(maxorb)
INCLUDE(common/common)
INCLUDE(common/iofile)
INCLUDE(common/cndx40)
INCLUDE(common/cndx41)
      dimension q(*)
INCLUDE(common/mapper)
INCLUDE(common/atmblk)
INCLUDE(common/ghfblk)
INCLUDE(common/prnprn)
      character *8 grhf,oscf
      data grhf/'grhf'/
      data oscf/'oscf'/
c
c
      if (lfdtrn .or. irest.eq.1) then
        if ( odebug(30)) then
         write (iwr,6050)
         if (lfdtrn) write (iwr,6060)
         if (irest.eq.1) write (iwr,6070)
         write (iwr,6030)iochf(13),iochf(14)
        endif
       return
      end if
      if(odebug(30)) write (iwr,6040)
      i1 = num*num + 1
      i2 = i1 + nx
      lennew = iky(ncoorb+1)
      lenblk = lensec(nx)
      newblk = lensec(lennew)
      out = odebug(3) .or. odebug(4)
      nat3 = nat*3
cc      ib = iochf(11)
cc      ib2 = iochf(11)
cc      iochf(13) = ib2
      nblkq = num*ncoorb
      call secget(isect(8),8,iblkq)
      iblkq = iblkq + mvadd
c
c     derivatives of fock matrices (no wavefunction derivatives)
c     a.o. basis at section 11
c     m.o. basis at section 13 of fockfile
c    ******** section 13 overwrites section 11 *********
c
      if (out) then
         write (iwr,6010)
      end if
      call rdedx(q(1),nblkq,iblkq,ifild)
      maxn3 = nat3
      if (scftyp.eq.oscf) maxn3 = nat3 + nat3
      if (scftyp.eq.grhf) maxn3 = (njk+njk+1)*nat3
      do 20 i = 1 , maxn3

csecd         call rdedx(q(i2),nx,ib,ifockf)
         call fetch ( q(i2), nx, 'dh', i )

         call qhq1(q(i1),q(1),ilifq,ncoorb,q(i2),iky,num)
cc         call wrt3(q(i1),lennew,ib2,ifockf)
csecd
         call stash ( q(i1), nx, 'dh', i )
cc         ib = ib + lenblk
cc         ib2 = ib2 + newblk
         if (out) then
            call prtris(q(i1),ncoorb,iwr)
         end if
 20   continue
      if (out) then
         write (iwr,6020)
      end if
cc      ib = iochf(12)
cc      ib2 = iochf(12)
cc      iochf(14) = ib2
c
c     derivatives of overlap matrix
c     a.o. basis at section 12
c     m.o. basis at section 14
c     ******** section 14 overwrites section 12
c
      do 30 i = 1 , nat3

csecd         call rdedx(q(i2),nx, ib, ifockf)
         call fetch ( q(i2), nx, 'ds', i )

         call qhq1(q(i1),q(1),ilifq,ncoorb,q(i2),iky,num)
cc         call wrt3(q(i1),lennew,ib2,ifockf)
csecd
         call stash ( q(i1), nx, 'ds', i )
cc         ib = ib + lenblk
cc         ib2 = ib2 + newblk
         if (out) then
            call prtris(q(i1),ncoorb,iwr)
         end if
 30   continue
      lfdtrn = .true.
      irest = 1
      call revise
      call clredx
      if(odebug(30)) write (iwr,6030)iochf(13),iochf(14)
      return
 6010 format (///5x,'transformed one-electron matrix derivatives'//)
 6020 format (///5x,'transformed overlap matrix derivatives'//)
 6030 format (/1x,'hamfile summary'/
     +         1x,'section 13 at block ',i5/
     +         1x,'section 14 at block ',i5)
 6040 format(/1x,'calling trnfkd')
 6050 format(/1x,'omitting call of trnfkd')
 6060 format(1x,'because lfdtrn is true')
 6070 format(1x,'because irest = 1')
      end
c
      subroutine dksm_exp(q,iq)
      implicit none
c
c     Construct and add the Kohn-Sham contributions to the explicit
c     derivatives of the Fock matrixes. I.e. derivatives with respect
c     to nuclear coordinates and no wavefunction contributions.
c
c     Parameters
c
INCLUDE(common/sizes)
c
c     Input
c
INCLUDE(common/common)
INCLUDE(common/mapper)
INCLUDE(common/debug)
INCLUDE(common/infoa)
INCLUDE(common/iofile)
INCLUDE(common/cndx41)
INCLUDE(common/drive_dft)
INCLUDE(common/atmblk)
INCLUDE(common/dump3)
c
c     Workspace:
c
      REAL q(*)
      integer iq(*)
c
c     Functions:
c
INCLUDE(common/ccpdft.hf77)
INCLUDE(common/gmempara)
      logical opg_root
      integer igmem_alloc_inf
      integer igmem_null
      integer lensec
c
c     Local:
c
      integer inull
      character *8 grhf,oscf,zfock
      character *15 fnm
      character *8  snm
      data fnm/'secd_parallel.m'/
      data snm/'dksm_exp'/
      data grhf/'grhf'/
      data oscf/'oscf'/
      data zfock/'fockder'/
      logical ofock,out
      integer i0,i1,i3,lennew,lenblk,newblk,nat3,ib2,i,nfok,ifok
      integer nblkq,iblkq,nincr
      integer ierror
_IF(ccpdft)
      ofock = fkder.eq.zfock
      inull = igmem_null()
      if (CD_active().and.ofock) then
         lennew = iky(ncoorb+1)
         lenblk = lensec(nx)
         newblk = lensec(lennew)
         out = odebug(3) .or. odebug(4)
         nat3 = nat*3
c
c        Compute DFT contributions
c
         if (ks_dx_bas.eq.KS_DX_AO) then
c
c           create the fock matrices
c
            nincr = nx
            nfok  = nat3
            ifok  = igmem_alloc_inf(nincr*nfok,fnm,
     &                              snm,'der-fock',IGMEM_NORMAL)
            call vclr(q(ifok),1,nincr*nfok)
c
c           load the density matrix
c
            i1 = igmem_alloc_inf(nx,fnm,snm,
     &                           'alpha-dens-mat',IGMEM_NORMAL)
            call rdedx(q(i1),nx,ibl3pa,idaf)
c
c           calculate DFT contributions
c
            ierror = CD_dksm_exp_ao(iq,q,nfok,q(i1),q(inull),
     &                              q(ifok),q(inull),.false.,iwr)
c
c           load the KS vectors
c
            i0 = igmem_alloc_inf(num*num,fnm,snm,
     &                           'alpha-vectors',IGMEM_NORMAL)
            nblkq = num*ncoorb
            call secget(isect(8),8,iblkq)
            iblkq = iblkq + mvadd
            call rdedx(q(i0),nblkq,iblkq,ifild)
c
c           transform DFT contributions to MO-basis
c
            do i = 1, nfok
               i3 = ifok + (i-1)*nincr
               call qhq1(q(i1),q(i0),ilifq,ncoorb,q(i3),iky,num)
               call dcopy(lennew,q(i1),1,q(i3),1)
            enddo
c
            call gmem_free_inf(i0,fnm,snm,'alpha-vectors')
            call gmem_free_inf(i1,fnm,snm,'alpha-dens-mat')
c
         else if (ks_dx_bas.eq.KS_DX_MO) then
c
c           create the fock matrices
c
            nincr = lennew
            nfok  = nat3
            ifok  = igmem_alloc_inf(nincr*nfok,fnm,
     &                              snm,'der-fock',IGMEM_NORMAL)
            call vclr(q(ifok),1,nincr*nfok)
c
c           load the KS vectors
c
            i0 = igmem_alloc_inf(num*num,fnm,snm,
     &                           'alpha-vectors',IGMEM_NORMAL)
            nblkq = num*ncoorb
            call secget(isect(8),8,iblkq)
            iblkq = iblkq + mvadd
            call rdedx(q(i0),nblkq,iblkq,ifild)
c
c           calculate DFT contributions
c
            ierror = CD_dksm_exp_mo(iq,q,nfok,ncoorb,na,0,
     &               q(i0),q(inull),q(ifok),q(inull),.false.,iwr)
c
            call gmem_free_inf(i0,fnm,snm,'alpha-vectors')
c
         else
            if (opg_root()) then
               write(iwr,*)'dksm_exp: ks_dx_bas = ',ks_dx_bas
            endif
            call caserr('dksm_exp: ks_dx_bas has an illegal value!')
         endif
c
c        Add explicit derivative KS matrices onto stored quantities
c        derivatives of fock matrices (no wavefunction derivatives)
c        m.o. basis at section 13 of fockfile
c
         i1  = igmem_alloc_inf(lennew,fnm,snm,
     &                         'tmp-fock',IGMEM_NORMAL)
         ib2 = iochf(13)
         i3 = ifok
         do i = 1 , nfok
c           call rdedx(q(i1),lennew,ib2,ifockf)
            call fetch(q(i1),lennew,'dh',i)
            call daxpy(lennew,1.0d0,q(i3),1,q(i1),1)
c           call wrt3(q(i1),lennew,ib2,ifockf)
            call stash(q(i1),lennew,'dh',i)
c           ib2 = ib2 + newblk
            i3  = i3  + nincr
         enddo
c
         call revise
         call clredx
c
         call gmem_free_inf(i1,fnm,snm,'tmp-fock')
         call gmem_free_inf(ifok,fnm,snm,'der-fock')
      endif
      return
_ENDIF
_IF(old-junk)
      ofock = fkder.eq.zfock
      if (CD_active().and.ofock) then
         i0 = igmem_alloc(num*num)
         i1 = igmem_alloc(nx)
         lennew = iky(ncoorb+1)
         lenblk = lensec(nx)
         newblk = lensec(lennew)
         out = odebug(3) .or. odebug(4)
         nat3 = nat*3
         ib2 = iochf(13)
c
c        load vectors
c
         nblkq = num*ncoorb
         call secget(isect(8),8,iblkq)
         iblkq = iblkq + mvadd
         call rdedx(q(i0),nblkq,iblkq,ifild)
c
c        create the fock matrices
c
         nfok = nat3
         ifok = igmem_alloc(nx*nfok)
         call vclr(q(ifok),1,nx*nfok)
c
c        calculate fock matrix contributions
c
         if (.true.) then
            call rdedx(q(i1),nx,ibl3pa,idaf)
            ierror = CD_dksm_exp_ao(iq,q,nfok,q(i1),q(inull),
     &                              q(ifok),q(inull),.false.,iwr)
            do i = 1, nfok
               i3 = ifok + (i-1)*nx
               call qhq1(q(i1),q(i0),ilifq,ncoorb,q(i3),iky,num)
               call dcopy(lennew,q(i1),1,q(i3),1)
            enddo
         endif
c
c        derivatives of fock matrices (no wavefunction derivatives)
c        m.o. basis at section 13 of fockfile
c
         i3 = ifok
         do i = 1 , nfok
c           call rdedx(q(i1),lennew,ib2,ifockf)
            call fetch(q(i1),lennew,'dh',i)
            call daxpy(lennew,1.0d0,q(i3),1,q(i1),1)
c           call wrt3(q(i1),lennew,ib2,ifockf)
            call stash(q(i1),lennew,'dh',i)
c           ib2 = ib2 + newblk
            i3  = i3  + nx
         enddo
c
         call revise
         call clredx
c
         call gmem_free(ifok)
         call gmem_free(i1)
         call gmem_free(i0)
      endif
      return
_ENDIF
      end
c
      subroutine jobstep(string)
      integer itop, iout
      real*8 cpusec, cpulft, wall
      character*(*) string
      common/maxlen/maxq
INCLUDE(common/iofile)
INCLUDE(common/common)
      logical opg_root
      iout = 6
      cpusec=cpulft(1)
      call walltime(wall)
      if (nprint.ne.-5.and.opg_root()) write(iout,2) cpusec,
     &     wall, string
      return
 2    format(/,5x,'>>> cpu=',f8.2,' wall=',f8.2,' ',a60,/)
      end

      subroutine chfndr(q)
c
c    driving routine for nuclear displacement chf routines
c    -----------------------------------------------------
c
      implicit REAL  (a-h,o-z)
      logical lstop,skipp
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
      common/lsort/skipp(3*maxat)
      common/maxlen/maxq
INCLUDE(common/infoa)
INCLUDE(common/nshel)
INCLUDE(common/common)
INCLUDE(common/iofile)
INCLUDE(common/cndx40)
INCLUDE(common/cndx41)
INCLUDE(common/ghfblk)
      logical lfield,fixed,lex,ldam12,ldam13,ldam23,ldiis
      common/scfblk/en,etot,ehf,sh1(2),sh2(2),gap1(2),gap2(2),
     1              d12,d13,d23,canna,cannb,cannc,fx,fy,fz,
     2              lfield,fixed,lex,ldam12,ldam13,ldam23,ldiis,
     3              ncyc,ischm,lock,maxit,nconv,npunch,lokcyc
      common/mpshl/ns(maxorb)
INCLUDE(common/symtry)
INCLUDE(common/statis)
INCLUDE(common/global)
INCLUDE(common/vcore)
INCLUDE(common/timeperiods)
INCLUDE(common/mp2grd)
INCLUDE(common/mapper)
INCLUDE(common/secd_pointers)
INCLUDE(common/drive_dft)
_IF(ccpdft)
INCLUDE(common/ccpdft.hf77)
_ENDIF
      logical ogeompert_save
      logical pg_create_inf
      logical pg_destroy_inf
c
      dimension q(*)
      character*10 charwall
      character *8  blank,closed
      character *15 fnm
      character *6  snm
      data closed/'closed'/
      data blank/'        '/
      data fnm/'secd_parallel.m'/
      data snm/'chfndr'/
      call cpuwal(begin,ebegin)
c
cc      lenx = lensec(nx)*nat*15
cc      call wrt3z(iblks,ifils,lenx)
      write (iwr,6010) cpulft(1) ,charwall()
      if (scftyp.ne.closed) then
       call caserr('parallel hessian for closed shells only')
      end if
_IF(ccpdft)
      ierror = CD_set_2e()
_ENDIF
      np = nat*3
         mn = noccb*nvirta + (noccb-nocca)*nocca
         call grhfbl(scftyp)
         call bfnshl(ns,nsa4)
c        iso = lenrel(mn) + 1
c        i0 = iso + lenrel(nw196(5))
c        i01 = mn + lenint(nx) + nw196(5) + 1
c        i2 = i01 + mn
c        i3 = i2 + mn
c        i4 = i3 + mn
c        i5 = i4 + num
c        i6 = i5 + nx
c        i7 = i6 + nx
c
c soln and pdens stash memory
         len3 = np * ikyp(ncoorb)
         len2 = np * mn
         if (omem2nd) then
           if (omem2nd_ga) then
c pdens stash memory
             if (.not.pg_create_inf(0,len3,1,'pdens-stash',0,0,
     &                iof2nd(5),fnm,snm,IGMEM_NORMAL)) then
               call caserr('failed to allocate pdens-stash')
             endif
c soln stash memory
             if (.not.pg_create_inf(0,len2,1,'soln-stash',0,0,
     &                iof2nd(6),fnm,snm,IGMEM_NORMAL)) then
               call caserr('failed to allocate soln-stash')
             endif
           else
c pdens stash memory
              iof2nd(5) = igmem_alloc_inf(len3,fnm,snm,'pdens-stash',
     &                                    IGMEM_NORMAL)
c soln stash memory
              iof2nd(6) = igmem_alloc_inf(len2,fnm,snm,'soln-stash',
     &                                    IGMEM_NORMAL)
           endif
         endif
         iso= igmem_alloc_inf(nw196(5),fnm,snm,'iso',IGMEM_DEBUG)
         i0 = igmem_alloc_inf(mn,fnm,snm,'eps',IGMEM_DEBUG)
         if (omem2nd) then
           if (omem2nd_ga) then
c rhs stash memory
             if (.not.pg_create_inf(0,len2,1,'rhs-stash',0,0,
     &                iof2nd(3),fnm,snm,IGMEM_NORMAL)) then
               call caserr('failed to allocate rhs-stash')
             endif
           else
c rhs stash memory
              iof2nd(3) = igmem_alloc_inf(len2,fnm,snm,'rhs-stash',
     &                                    IGMEM_NORMAL)
           endif
         endif
chvd     i00 = igmem_alloc_inf(nx,fnm,snm,'mapnr',IGMEM_DEBUG)
c
         i1 = igmem_alloc_inf(mn,fnm,snm,'bx',IGMEM_DEBUG)
         i2 = igmem_alloc_inf(mn,fnm,snm,'by',IGMEM_DEBUG)
         i3 = igmem_alloc_inf(mn,fnm,snm,'bz',IGMEM_DEBUG)
         i4 = igmem_alloc_inf(num,fnm,snm,'eval',IGMEM_DEBUG)
         i5 = igmem_alloc_inf(nx,fnm,snm,'sx',IGMEM_DEBUG)
         i6 = igmem_alloc_inf(nx,fnm,snm,'sy',IGMEM_DEBUG)
         i7 = igmem_alloc_inf(nx,fnm,snm,'sz',IGMEM_DEBUG)
c
c      sort out a-matrix
c
c     r.h.s of equations
c
         call start_time_period(TP_2D_CHFRHS)

         ibuff = igmem_alloc_inf(nvirta*nocca,fnm,snm,'buf',IGMEM_DEBUG)

c        call jobstep('rhs')
         write(iwr,7010) cpulft(1) ,charwall()
         call rhscl_ga(q(i0),q(i1),q(i2),q(i3),q(i4),
     +        q(i5),q(i6),q(i7), 
     +        q(ibuff))
c
         call gmem_free_inf(ibuff,fnm,snm,'buf')
c
         call gmem_free_inf(i7,fnm,snm,'sz')
         call gmem_free_inf(i6,fnm,snm,'sy')
         call gmem_free_inf(i5,fnm,snm,'sx')
         call gmem_free_inf(i4,fnm,snm,'eval')
         call gmem_free_inf(i3,fnm,snm,'bz')
         call gmem_free_inf(i2,fnm,snm,'by')
         call gmem_free_inf(i1,fnm,snm,'bx')
c
         call rhscl_dft(q,q)
c
         call end_time_period(TP_2D_CHFRHS)
c
         call start_time_period(TP_2D_SYMMRHS)
c        call jobstep('symrhs')
         write(iwr,7000)cpulft(1) ,charwall()
         i00 = igmem_alloc_inf(nx,fnm,snm,'mapnr',IGMEM_DEBUG)
         itmp = (nat * 3 + 1) * mn
         icore = igmem_alloc_inf(itmp,fnm,snm,'core_symrhs',IGMEM_DEBUG)
cc         call symrhs(q(icore),iblks,skipp,q(i00),q(iso),nshell)
         call symrhs(q(icore),skipp,q(i00),q(iso),nshell)
         call gmem_free_inf(icore,fnm,snm,'core_symrhs')
         call gmem_free_inf(i00,fnm,snm,'mapnr')
         call end_time_period(TP_2D_SYMMRHS)
c
c       lenblk = lensec(mn)
c       iblku = iblks + np*lenblk
      lstop = .false.
c
c   solve linear equations
c
c     opg_pertbns_sep = .true.     
c     now set by "runtype hessian separate"
c
      call start_time_period(TP_2D_CHFDRV)
      write(iwr,7020) cpulft(1)
c     call jobstep('chfdrv')
      ogeompert_save = ogeompert
      ogeompert = .true.
      if (opg_pertbns_sep) then
c
c   old ga solver - perturbations done separately
         call chfdrv(q(i0),lstop,skipp)
      else
c
c   new ga solver - perturbations done together
         call chfdrv2(q(i0),lstop,skipp)
      end if
      ogeompert = ogeompert_save
c   free stashed 'rhs' memory
c
      if (omem2nd) then
        if (omem2nd_ga) then
          if (.not.pg_destroy_inf(iof2nd(3),'rhs-stash',fnm,snm))
     &             call caserr('failed to destroy rhs-stash')
        else
          call gmem_free_inf(iof2nd(3),fnm,snm,'rhs-stash')
        endif
      endif
      call gmem_free_inf(i0,fnm,snm,'eps')
      call end_time_period(TP_2D_CHFDRV)

      if (lstop) then
         call revise
         write (iwr,6040)
         call timana(22)
         call clenms('stop')
      else
c
         call start_time_period(TP_2D_SYMMU)
c        call jobstep('symmu')
         write(iwr,7030)cpulft(1) ,charwall()
         itmp = (nat + 1) * mn * 3
         i01 = igmem_alloc_inf(itmp,fnm,snm,'core_symmu',IGMEM_DEBUG)
cc         call symmu(q(i01),iblku,skipp,q(iso),nshell)
         call symmu(q(i01),skipp,q(iso),nshell)
         call gmem_free_inf(i01,fnm,snm,'core_symmu')
         call end_time_period(TP_2D_SYMMU)
c
         call gmem_free_inf(iso,fnm,snm,'iso')
c
c     perturbed density matrices
c
         call start_time_period(TP_2D_PDENS)
c        call jobstep('pdens')
         write(iwr,7040)cpulft(1) ,charwall()
         itmp = nx + nx + mn
         i01 = igmem_alloc_inf(itmp,fnm,snm,'core_pdens',IGMEM_DEBUG)
         call pdens(q(i01),lstop)
         call gmem_free_inf(i01,fnm,snm,'core_pdens')
c
         call end_time_period(TP_2D_PDENS)
c
c   free stashed 'soln' memory
         if (omem2nd) then
           if (omem2nd_ga) then
             if (.not.pg_destroy_inf(iof2nd(6),'soln-stash',fnm,snm))
     &                call caserr('failed to destroy soln-stash')
           else
             call gmem_free_inf(iof2nd(6),fnm,snm,'soln-stash')
           endif
         endif
         mmaxq = maxq
c
         call start_time_period(TP_2D_PFOCK)
c        call jobstep('pfock')
         write(iwr,7050)cpulft(1) ,charwall()
         if (omem2nd) then
           if (omem2nd_ga) then
c fock stash memory
             if (.not.pg_create_inf(0,len3,1,'fock-stash',0,0,
     &                iof2nd(4),fnm,snm,IGMEM_NORMAL)) then
                call caserr('failed to allocate fock-stash')
             endif
           else
c fock stash memory
             iof2nd(4) = igmem_alloc_inf(len3,fnm,snm,'fock-stash',
     &                                   IGMEM_NORMAL)
           endif
         endif
c     determine core available

         no = nocca
         nv = nvirta
         noo2=no*(no+1)/2
         nvv = nv*nv
         ibuf =  max(noo2, nvv)
         ibase2 = igmem_alloc_inf(ibuf,fnm,snm,'buf_pfockc',IGMEM_DEBUG)
         maxqq = igmem_max_memory() 
     &         - memreq_pg_dgop(3*nat*nx,'+') 
     &         - 2*igmem_overhead()
         ibase1 = igmem_alloc_inf(maxqq,fnm,snm,'core_pfockc',
     &                            IGMEM_DEBUG)

         call pfockc(q(ibase1),q(ibase2),maxqq)

         call gmem_free_inf(ibase1,fnm,snm,'core_pfockc')
         call gmem_free_inf(ibase2,fnm,snm,'buf_pfockc')

         call pdksmc(q,q)

         call end_time_period(TP_2D_PFOCK)
c
         maxq = mmaxq
c
         call revise
         call clredx
         call delfil(nofile(1))
         call timana(22)
      end if
_IF(ccpdft)
      ierror = CD_reset_2e()
_ENDIF
 6010 format (/1x,
     + 'commence solution of chf equations for nuclear motions at ',
     +  f9.2,' seconds',a10,' wall')
 6040 format (//1x,'insufficient time to finish chf equations'//1x,
     +        'restart job'//)
 7000 format(/1x,'commence symmetrisation of r.h.s. at ', 
     +       f9.2,' seconds',a10,' wall'/)
 7010 format(/1x,'commence GA-based r.h.s. construction at ', 
     +       f9.2,' seconds',a10,' wall')
 7020 format(/1x,'commence solution of linear equations at ', 
     +       f9.2,' seconds',a10,' wall'/)
 7030 format(/1x,'commence symmu at ',f9.2,' seconds',a10,' wall')
 7040 format(/1x,
     + 'commence construction of perturbed density matrices at ',
     +  f9.2,' seconds',a10,' wall')
 7050 format(/1x,
     +'commence final construction of derivative fock operator at ',
     + f9.2,' seconds',a10,' wall')
      end
      subroutine chfdrv(eps,lstop,skipp)
      implicit REAL  (a-h,o-z)
c
c  driving routine for cphf - perturbations treated separately
c
      logical lstop,skipp
      dimension skipp(*), eps(*)
INCLUDE(common/iofile)
INCLUDE(common/common)
INCLUDE(common/cndx41)
INCLUDE(common/global)
INCLUDE(common/vcore)
INCLUDE(common/gmempara)
      logical pg_create_inf, opg_root
      logical pg_destroy_inf, succeeded
      dimension iatms(1)
      character *15 fnm
      character *6 snm
      data fnm/"secd_parallel.m"/
      data snm/"chfdrv"/

      npstar = 0
      npfin = np
      nov = nocca*nvirta
      nvv = nvirta*nvirta
      maxc = 50
      mynode = ipg_nodeid()
c
c  allocate distributed and replicated workspaces
c                                         
      if (pg_create_inf(0,nov,maxc,'u',nov,1,g_u,fnm,snm,IGMEM_DEBUG)) 
     &then
         call pg_distribution(g_u, mynode, ilo, ihi, jlo, jhi)
         write(*,918)  g_u, mynode, ilo, ihi, jlo, jhi
        il_u = ilo
        ih_u = ihi
        jl_u = jlo
        jh_u = jhi
        call pg_zero(g_u)
      else
        print *,'**GA error** failed to create GA ',g_u 
      end if

      i00 = igmem_alloc_inf(ncoorb,fnm,snm,"e",IGMEM_DEBUG)
      i10 = igmem_alloc_inf(nov,fnm,snm,"u",IGMEM_DEBUG)
      i20 = igmem_alloc_inf(nov,fnm,snm,"unxt",IGMEM_DEBUG)
      i30 = igmem_alloc_inf(nov,fnm,snm,"prhs",IGMEM_DEBUG)
      i40 = igmem_alloc_inf(maxc,fnm,snm,"b",IGMEM_DEBUG)
      i50 = igmem_alloc_inf(maxc,fnm,snm,"cc",IGMEM_DEBUG)
      i60 = igmem_alloc_inf(maxc,fnm,snm,"uu",IGMEM_DEBUG)
      i70 = igmem_alloc_inf(maxc*maxc,fnm,snm,"uau",IGMEM_DEBUG)
      i80 = igmem_alloc_inf(max( nvv , max( 2*maxc , nov )),fnm,snm,
     &                      "buf",IGMEM_DEBUG)
      i90 = igmem_alloc_inf(nvv,fnm,snm,"buf2",IGMEM_DEBUG)
      iz  = igmem_alloc_inf(nov,fnm,snm,"rhs",IGMEM_DEBUG)

c  loop over perturbations
      ip = 0
      do i = npstar + 1 , npfin
         ip = ip + 1
         if (skipp(i)) then
          if(opg_root())write(iwr,*)
     +    'Skip CPHF for perturbation: ',ip
            call vclr(Q(iz),1,nov)
            call stash(Q(iz),nov,'soln',ip)
         else
           if(opg_root())write(iwr,*)
     +        'CPHF for perturbation: ',ip 
c  using fetch
            call fetch(Q(iz),nov,'rhs',ip)
c  solver
            iatms(1) = (ip-1)/3+1
            call secdchf_ga(Q(i00),Q(iz),
     &           Q(i10),Q(i20),Q(i30),
     &           Q(i40),Q(i50),
     &           Q(i60),Q(i70),Q(i80),
     &           Q(i90),maxc,nov,eps,iatms,Q(1))
            call stash(Q(iz),nov,'soln',ip)
        end if
      end do

c  deallocate distributed and replicated workspaces
      call gmem_free_inf(iz,fnm,snm,"rhs")
      call gmem_free_inf(i90,fnm,snm,"buf2")
      call gmem_free_inf(i80,fnm,snm,"buf")
      call gmem_free_inf(i70,fnm,snm,"uau")
      call gmem_free_inf(i60,fnm,snm,"uu")
      call gmem_free_inf(i50,fnm,snm,"cc")
      call gmem_free_inf(i40,fnm,snm,"b")
      call gmem_free_inf(i30,fnm,snm,"prhs")
      call gmem_free_inf(i20,fnm,snm,"unxt")
      call gmem_free_inf(i10,fnm,snm,"u")
      call gmem_free_inf(i00,fnm,snm,"e")

      succeeded=pg_destroy_inf(g_u,'u',fnm,snm)
      if (.not.succeeded)
     &        call caserr('failed to destroy u')

      return
 6010 format (//1x,'insufficient store for chf equations'//1x,
     +        'store available ',i8/1x,'required - at least ',i8,
     +        ' and preferably ',i8)
 918  format(' pg_distribution of GA [',i6,'] to node ',1i4,':'
     &    ,4x,'ilo =',1i6,3x,'ihi =',1i6,3x,'jlo =',1i6,3x,'jhi =',1i6)
      end

      subroutine secdchf_ga(
     &     e,       ! ncoorb       i00
     &     rhs,     ! nov          iz
     &     u,       ! nov          i10
     &     unxt,    ! nov          i20
     &     prhs,    ! nov          i30
     &     b,       ! maxc         i40
     &     cc,      ! maxc         i50
     &     uu,      ! maxc         i60
     &     uau,     ! maxc,maxc    i70  
     &     buf,     !              i80
     &     buf2,    !              i90
     &     maxc,    !      #iters
     &     nov,
     &     eps,
     &     iatms,
     &     q)            

      implicit REAL  (a-h,o-z)
INCLUDE(common/timez)
INCLUDE(common/common)
INCLUDE(common/cndx41)
INCLUDE(common/prnprn)
INCLUDE(common/iofile)
INCLUDE(common/timeperiods)
      integer a
      common/small/alpha(50,50),aa(50,50),wk1(50),wk2(50)
      dimension e(*), rhs(nov)
      dimension u(nov),unxt(nov),prhs(nov)
      dimension b(maxc),cc(maxc),uu(maxc),uau(maxc,maxc)
      dimension buf(*), buf2(*)
      dimension eps(nov)
      dimension iatms(1)
      dimension q(*)
      character*10 charwall
INCLUDE(common/global)
      data smal/1.0d-13/,tich/1.0d-24/
      data four/4.0d0/,zero/0.0d0/
      no = nocca
      nv = nvirta
      n = ncoorb
      no1 = no + 1
      nvv = nv*nv

      call start_time_period(TP_MP2CHF)
c
c  iterative solution of simultaneous equations to give
c  coupled hartree-fock first order wavefunction
c
      nnodes = ipg_nnodes()
      mynode = ipg_nodeid()
      if (maxc.gt.50) maxc = 50
      uconv = 10.0d0**(-iconvv)
      iaa = 50
      ifail = 0
      if (oprn(12)) write (iwr,6010)
      if (oprn(13)) write (iwr,6020)
      call vclr(uau,1,maxc*maxc)
      call vclr(b,1,maxc)
      call vclr(uu,1,maxc)

      write(iwr,6000) cpulft(1) ,charwall()

c      write(6,*)'zero order.. RHF'
c      write(6,99)(rhs(iii),iii=1,nov)        

c  get zeroth order estimate
      do i = 1 , nov
         v = -rhs(i)*eps(i)
         if (dabs(v).le.tich) v = 0.0d0
         b(1) = b(1) + v*v
         u(i) = v
      end do

c
c  copy U into 1st locn of global array
c
c      write(6,*)'put u:',1,nov,g_u
      if (1.ge.jl_u.and.1.le.jh_u)
     &     call pg_put(g_u, 1, nov, 1, 1, u, nov)

c      call wrtvec_ga('rhs   ',nov,u)

      if (oprn(6)) write(iwr,6100)
c
c  start of iterative solution of chf equations
c  50 iterations are allowed ---  usually less than 10 are necessary
c
      do 1 itr = 1 , maxc
        call vclr(unxt,1,nov)
        
c        write(6,*)'enter MXM : estimate'
c        write(6,99)(u(iii),iii=1,nov)        
 99      format(1x,8f10.5)

        call mxmau_ga_secd(u,unxt,buf,buf2,no,nv,nov,nvv)
        call pg_dgop(101+itr,unxt,nov,'+')
_IF(ccpdft)
c
c       Add the DFT contributions
c
        call au_dft(q,u,unxt,1,iatms)
_ENDIF
c        write(6,*)'after MXM and gop: unxt'
c        write(6,99)(unxt(iii),iii=1,nov)        

c  scale unxt by difference of eigenvalues

        do 70 i = 1 , mn
           unxt(i) = unxt(i)*eps(i)
 70     continue

        uu(itr) = ddot(nov,u,1,u,1)

c        write(6,*)'uu:',uu(itr)

        do i = 2 , itr - 1
           if (i.ge.jl_u.and.i.le.jh_u) then
              call pg_get(g_u, 1, nov, i, i, buf, nov)
              uau(itr,i-1) = ddot(nov,u,1,buf,1)
           end if
        end do

        if (mynode.eq.nnodes-1.and.itr.gt.1)
     &       uau(itr,itr - 1) = uu(itr)
        if (mynode.eq.nnodes-1)
     &       uau(itr,itr) = ddot(nov,u,1,unxt,1)
        call dcopy(itr,uau(itr,1),maxc,buf,1)
        call pg_dgop(202+itr,buf,itr,'+')
        call dcopy(itr,buf,1,uau(itr,1),maxc)

        do i = 1 , itr - 1
           if (i.ge.jl_u.and.i.le.jh_u) then
              call pg_get(g_u, 1, nov, i, i, buf, nov)
              uau(i,itr) = ddot(nov,buf,1,unxt,1)
           end if
        end do
        call dcopy(itr-1,uau(1,itr),1,buf,1)
        call pg_dgop(303+itr,buf,itr-1,'+')
        call dcopy(itr-1,buf,1,uau(1,itr),1)

c        write(6,*)'uau:',((uau(iii,jjj),iii=1,itr),jjj=1,itr)

c  nag routine to solve a small set of simultaneous equations
        do i = 1 , itr
          do j = 1 , itr
            alpha(j,i) = uau(j,i)
          end do
        end do
        do i = 1 , itr
          alpha(i,i) = alpha(i,i) + uu(i)
        end do
        call f04atf(alpha,iaa,b,itr,cc,aa,iaa,wk1,wk2,ifail)
c  form new solution vector
        call dcopy(nov,rhs,1,prhs,1)
        call vclr(rhs,1,nov)

c        write(6,*)'cc:',(cc(iii),iii=1,itr)

        do i = 1 , itr
          ccjn = cc(i)
          if (dabs(ccjn).gt.tich) then
            if (i.ge.jl_u.and.i.le.jh_u) then
              call pg_get(g_u, 1, nov, i, i, buf, nov)
              call daxpy(nov,ccjn,buf,1,rhs,1)
            end if
          end if
        end do
        call pg_dgop(404+itr,rhs,nov,'+')

c  check for convergence of solution vectors
        if (itr.ne.1) then
          gmax = zero
          call vsub(rhs,1,prhs,1,prhs,1,nov)
          gnorm = ddot(nov,prhs,1,prhs,1)/dfloat(nov)
          gnorm = dsqrt(gnorm)
          gmax = dmax1(gmax,gnorm)
          if (oprn(6)) write (iwr,6050) itr , gmax
          if (gmax.le.uconv) then
            write (iwr,6090)
            its = itr
            go to 3
          end if
        end if

c  update next expansion vector
        call vclr(prhs,1,nov)
        do i = 1 , itr
          if (i.ge.jl_u.and.i.le.jh_u) then
            fac = -uau(i,itr)/uu(i)
            call pg_get(g_u, 1, nov, i, i, buf, nov)
            call daxpy(nov,fac,buf,1,prhs,1)
          end if
        end do
        call pg_dgop(505+itr,prhs,nov,'+')
        call vadd(prhs,1,unxt,1,unxt,1,nov)

        gmax = zero
        do i = 1 , nov
          if (dabs(unxt(i)).le.tich) unxt(i) = zero
        end do
c  check for convergence of next expansion vector
        gnorm = ddot(nov,unxt,1,unxt,1)/dfloat(nov)
        gnorm = dsqrt(gnorm)
        gmax = dmax1(gmax,gnorm)
        if (oprn(6)) write (iwr,6060) gmax
        if (gmax.le.smal) then
          write (iwr,*)'converged - new expansion vector negligible '
          its = itr
          go to 3
        end if

        nxtr = itr + 1

        if (nxtr.ge.jl_u.and.nxtr.le.jh_u)then
c           write(6,*)'put u:',1,nov,nxtr,nxtr
           call pg_put(g_u, 1, nov, nxtr, nxtr, unxt, nov)
        endif

        call dcopy(nov,unxt,1,u,1)
1     continue
c
c  end of loop
c
      write (iwr,*) ' no full convergence after ',maxc,' iterations '
      write (iwr,*) ' this will require changes to the program ! '
3     call timit(3)
      write (iwr,6070) itr , cpulft(1) ,charwall()
      if (oprn(13)) write (iwr,6080) 1,(rhs(i),i=1,nov)

      call end_time_period(TP_MP2CHF)

      return
 6000 format(/1x,
     +'commence iterative solution of chf equations at ',f8.2,
     +' seconds',a10,' wall')
 6010 format (//1x,'print right-hand-side to chf equations')
 6020 format (//1x,'print solutions to chf equations')
 6030 format (1x,'perturbation',i5,' omitted')
 6040 format (//1x,'perturbation  ',i4//(5x,5f16.8))
 6050 format (i10,5x,f15.10)
 6060 format (30x,f20.15)
 6070 format (/1x,
     + 'chf converged at iteration',i4/1x,
     + 'chf complete at ',f8.2,' seconds',a10,' wall')
 6080 format (//1x,'solution  ',i4//(5x,5f16.8))
 6090 format(/1x,'chf converged - wavefunctions stationary')
 6100 format(/
     +  6x,'iteration',9x,'tester',2x,'expansion vector norm'/
     +  6x,47('=')/)
      end
      subroutine mxmau_ga_secd(a1,a2,buf1,buf2,no,nv,nov,nvv)
      implicit REAL  (a-h,o-z)

      dimension a1(nov), a2(nov), buf1(*), buf2(nv,*)
INCLUDE(common/global)
      integer i,j,k,l,kl,no,nv,nvv,ij,o1,o2,v1,v2
      integer nvt

      nvt = nv*(nv+1)/2
      nvv = nv*nv

_IF(ccpdft)
      hf_wght = CD_HF_exchange_weight()
_ENDIF

      do i = 1,nov
         a2(i) = 0.0d0
      enddo

      ij = 0
      do o1=1,no
         do o2=1,o1
            ij = ij + 1
            if (ij.ge.jl_vvoo.and.ij.le.jh_vvoo) then

c               call pg_get(g_hess,1,nvv,ij,ij,buf2,nvv)
c  orbital hessian generated on-the-fly
          call pg_get(g_vovo,1,nvv,ij,ij,buf1,nvv)
          do k=1,nv
           do l=1,nv
_IF(ccpdft)
            buf2(l,k)=4.0d0*buf1((l-1)*nv+k)-hf_wght*buf1((k-1)*nv+l)
_ELSE
            buf2(l,k)=4.0d0*buf1((l-1)*nv+k)-buf1((k-1)*nv+l)
_ENDIF
           end do
          end do
          call pg_get(g_vvoo,1,nvt,ij,ij,buf1,nvt)
          kl=0
          do k=1,nv
           do l=1,k-1
            kl=kl+1
_IF(ccpdft)
            buf2(l,k) = buf2(l,k)-hf_wght*buf1(kl)
            buf2(k,l) = buf2(k,l)-hf_wght*buf1(kl)
_ELSE
            buf2(l,k) = buf2(l,k)-buf1(kl)
            buf2(k,l) = buf2(k,l)-buf1(kl)
_ENDIF
           end do
           kl=kl+1
_IF(ccpdft)
           buf2(l,l) = buf2(l,l)-hf_wght*buf1(kl)
_ELSE
           buf2(l,l) = buf2(l,l)-buf1(kl)
_ENDIF
          end do
       
               do v1 = 1,nv
                  do v2 = 1,nv
                     ix1 = (v1-1)*no + o1
                     ix2 = (v2-1)*no + o2
c                     a2(ix2) = a2(ix2) + buf1((v1-1)*nv+v2) * a1(ix1)
                     a2(ix2) = a2(ix2) + buf2(v2,v1) * a1(ix1)
                     if( o1 .ne. o2 )then
c                        a2(ix1) = a2(ix1) + buf1((v1-1)*nv+v2) * a1(ix2)
                        a2(ix1) = a2(ix1) + buf2(v2,v1) * a1(ix2)
                     endif
                  enddo
               enddo

            endif	
         enddo
      enddo
      return
      end

      subroutine pfockc(q,buff,maxq)
c
c     adds derivative wavefunction term to derivative integral
c     term to make complete derivative of fock operator
c
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      common/blkin/g(510),nword
      common/craypk/labs(1)
INCLUDE(common/infoa)
INCLUDE(common/nshel)
INCLUDE(common/symtry)
      dimension q(*),buff(*)
c
INCLUDE(common/common)
INCLUDE(common/iofile)
INCLUDE(common/cndx40)
INCLUDE(common/cndx41)
INCLUDE(common/prnprn)
INCLUDE(common/mapper)
INCLUDE(common/vcore)
      nat3 = nat*3
      ltri = ikyp(ncoorb)
cc      ibm1 = kblk(1)
cc      ibm2 = nblk(1)
cc      idevm = nufile(1)
cc      length = lensec(ltri)
cc      ibt = iochf(1)
c
c     perturbed density matrices at section 15
c     derivatives of integrals at section 13
c     complete fock matrix output at section 16 (m.o. basis)
c
cc      ibh = iochf(13)
cc      ibd = iochf(15)
cc      iochf(16) = ibt
c
c      modifications for multipass version
c
      mfok = maxq/(ltri+ltri)
c     mfok is maximum no. of fock matrices per pass , therefore
      npass = (nat3-1)/mfok + 1
      nadd = min(nat3,mfok)
      mi = 1
      ma = nadd
      mhalf = nadd*ltri
      mtwo = mhalf + mhalf
      mthalf = mtwo + mhalf
c     write(iwr,*)' nadd,mi,ma,mhalf',nadd,mi,ma,mhalf
      call setsto(1360,0,labs)


      ipert  = 1
      ipert2 = 1
      ipert3 = 1


      do 150 ipass = 1 , npass
         i = mhalf + 1
         if(npass.gt.1) write(iwr,*)'pfockc: ipass,npass',
     +      ipass,npass
c    read in batch of density matrices
         do 20 n = mi , ma

cc            call rdedx(q(i),ltri,ibd,ifockf)
            call fetch(q(i),ltri,'pdens',ipert)

            ipert = ipert + 1
            i = i + ltri
cc            ibd = ibd + length
 20      continue
c
c     map to active only
c
         ija = 0
         do 50 i = 1 , nsa4
            do 40 j = 1 , i
               ija = ija + 1
               ij = iky(mapie(i)) + mapie(j)
               k = 0
               l = mhalf
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
               do 30 n = mi , ma
                  q(ija+k) = 0.5d0*q(ij+l)
                  k = k + ltri
                  l = l + ltri
 30            continue
 40         continue
 50      continue
c
c     multiply diagonal of density matrix by 0.5
c
         do 70 i = 1 , nsa4
            ii = iky(i+1)
            k = 0
            do 60 n = mi , ma
               q(ii+k) = 0.5d0*q(ii+k)
               k = k + ltri
 60         continue
 70      continue

c
c     construct fock matrix in q(mhalf+1)
c
c  - integral version
c
csecd
cc         call vclr(q(mhalf+1),1,ltri*nadd)
cc         call search(ibm1,idevm)
cc         call find(idevm)
cc         do 80 ib = ibm1 , ibm2
cc            call get(g,nw)
cc            if (nword.eq.0) go to 90
cc            if (nw.eq.0) go to 90
cc            call find(idevm)
cc            call sgmatm(q(mhalf+1),q(1),nadd,ltri)
cc 80      continue
cc
cc 90      continue
cc
cc         do i = 1,nadd
cc            write(6,*)'fock (ints)..',i
cc            call prtri(q(mhalf +(i-1)*ltri + 1), nsa4)
cc         enddo
c
c  - global array implementation
c
         call vclr(q(mhalf+1),1,ltri*nadd)


         no = nocca
         nv = nvirta

         noo2=no*(no+1)/2
         nvv = nv*nv

c        ibuf=igmem_alloc( max(noo2, max(nvv, nadd*ltri)))

         call pfock_ga(q(mhalf+1),q(1),nadd,ltri, 
     &        buff)

c        call gmem_free(ibuf)

c
c     get term involving integral derivatives
c     and form complete derivative fock matrix
c
         i = 1
         do 100 n = mi , ma
cc
cc            call rdedx(q(i),ltri,ibh,ifockf)
            call fetch(q(i),ltri,'dh',ipert2)

            i = i + ltri
cc            ibh = ibh + length
            ipert2 = ipert2 + 1
 100     continue
         ija = 0
         do 130 i = 1 , nsa4
            do 120 j = 1 , i
               ija = ija + 1
               ij = iky(mapie(i)) + mapie(j)
               k = 0
               l = mhalf
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
               do 110 n = mi , ma
                  q(ij+k) = q(ij+k) + q(ija+l)
                  k = k + ltri
                  l = l + ltri
 110           continue
 120        continue
 130     continue
         i = 1
         do 140 n = mi , ma
cc
cc            call wrt3(q(i),ltri,ibt,ifockf)
            call stash(q(i), ltri, 'fock',ipert3)

            if (odebug(4)) then
               write (iwr,6010) n
               call prtris(q(i),ncoorb,iwr)
            end if
            i = i + ltri
cc            ibt = ibt + length
            ipert3 = ipert3 + 1
 140     continue
         mi = mi + nadd
         ma = ma + nadd
         ma = min(ma,nat3)
 150  continue
cc      iochf(1) = ibt
      call revise
      call clredx
      if(odebug(30)) then
       write (iwr,6020) iochf(13),iochf(14),iochf(15),iochf(16)
      endif
      if (.not.mp2) return
      call symfck(q(nw196(5)+1),q,nshell)
      i1 = 1 + ltri
      i2 = i1 + ltri
      i3 = i2 + ncoorb
      i4 = i3 + ncoorb*ncoorb*nat3
      call umat(q(1),q(i1),q(i2),q(i3),ltri,ncoorb,nat3)
      call revise
      return
 6010 format (//5x,'perturbed fock matrix in m.o. basis -- perturbation'
     +        ,i6/)
 6020  format(/1x,'hamfile summary'/
     +         1x,'section 13 at block ',i5/
     +         1x,'section 14 at block ',i5/
     +         1x,'section 15 at block ',i5/
     +         1x,'section 16 at block ',i5/)
      end
c
c
      subroutine pdksmc(q,iq)
      implicit none
c
c     Add the DFT wavefunction contribution to the perturbed Fock
c     matrix
c
c     Parameters:
c
INCLUDE(common/sizes)
      integer m8
      parameter (m8=8)
c
c     Input:
c
INCLUDE(common/mapper)
INCLUDE(common/infoa)
INCLUDE(common/atmblk)
INCLUDE(common/common)
INCLUDE(common/iofile)
INCLUDE(common/cndx41)
INCLUDE(common/cndx40)
INCLUDE(common/mp2grd)
c
c     Workspace:
c
      integer iq(*)
      REAL q(*)
c
c     Local:
c
      integer ibt, ibd, ltri, iatms, iblok
      integer nat3, i, j, k, l, ierror, n, m, lennew
      integer ida, ifa, icc, ida_t
      integer ij, ija
      integer ntot, nnow
      integer inull
      character *15 fnm
      character *6  snm
c
c     Functions:
c
INCLUDE(common/ccpdft.hf77)
INCLUDE(common/gmempara)
      integer lensec,lenrel,iiatms
      integer igmem_alloc_inf, igmem_null
      data fnm/'secd_parallel.m'/
      data snm/'pdksmc'/
c
c     Code:
c
      if (.not.CD_active()) return
c
      inull = igmem_null()
      if (opg_pertbns_sep) then
c
c        Do everything in batches of 1 coordinate at a time
c        if this code works we should always use it just with 
c        different batch sizes...
c
         ltri = ikyp(ncoorb)
         nat3 = 3*nat
         ntot = 1
         ida = igmem_alloc_inf(ltri*ntot,fnm,snm,
     &                         'pert-density',IGMEM_NORMAL)
         ifa = igmem_alloc_inf(ltri*ntot,fnm,snm,
     &                         'pert-fock',IGMEM_NORMAL)
         icc = igmem_alloc_inf(num*num,fnm,snm,
     &                         'alpha-vectors',IGMEM_NORMAL)
         iatms = igmem_alloc_inf(nat3,fnm,snm,
     &                           'iatms',IGMEM_NORMAL)
         lennew = lensec(ltri)
         iiatms = lenrel(iatms-1)+1
c
         call secget(isect(8),m8,iblok)
         iblok = iblok + mvadd
         call rdedx(q(icc),num*ncoorb,iblok,ifild)
c
         nnow  = 0
         n     = 0
         m     = 0
         ida_t = ida
         do i = 1, nat3
c
c           collect upto ntot perturbed density matrices
c
            n    = n + 1
            nnow = nnow + 1
            iq(iiatms+nnow-1) = (i-1)/3+1
            call fetch(q(ida_t),ltri,'pdens',n)
c
            call dscal(ltri,0.5d0,q(ida_t),1)
            do j = 1, nocca
               q(ida_t-1+j*(j+1)/2) = 0.5d0*q(ida_t-1+j*(j+1)/2)
            enddo
            ida_t = ida_t + ltri
c
            if (nnow.eq.ntot.or.i.eq.nat3) then
c
c              got enough perturbed density matrices so do the work
c
               call vclr(q(ifa),1,ltri*nnow)
               ierror = CD_chf_dksm_mo(iq,q,nnow,iq(iiatms),ncoorb,
     &                  nocca,0,q(icc),q(inull),q(ida),
     &                  q(inull),q(ifa),q(inull),.false.,iwr)
               do j = 1, nnow
                  m = m + 1
                  call fetch(q(ida),ltri,'fock',m)
                  ija = 0
                  do k = 1, nsa4
                     do l = 1, k
                        ij  = iky(mapie(k)) + mapie(l) - 1
                        q(ida+ij) = q(ida+ij) + q(ifa+(j-1)*ltri+ija)
                        ija = ija + 1
                     enddo
                  enddo
                  call stash(q(ida),ltri,'fock',m)
               enddo
c
c              reset counters for the next batch
c
               nnow  = 0
               ida_t = ida
            endif
         enddo
         call gmem_free_inf(iatms,fnm,snm,'iatms')
         call gmem_free_inf(icc,fnm,snm,'alpha-vectors')
         call gmem_free_inf(ifa,fnm,snm,'pert-fock')
         call gmem_free_inf(ida,fnm,snm,'pert-density')
         call revise
         call clredx
      else
         ltri = ikyp(ncoorb)
         nat3 = 3*nat
         ida = igmem_alloc_inf(ltri*nat3,fnm,snm,
     &                         'pert-density',IGMEM_NORMAL)
         ifa = igmem_alloc_inf(ltri*nat3,fnm,snm,
     &                         'pert-fock',IGMEM_NORMAL)
         icc = igmem_alloc_inf(num*num,fnm,snm,
     &                         'alpha-vectors',IGMEM_NORMAL)
         iatms = igmem_alloc_inf(nat3,fnm,snm,
     &                           'iatms',IGMEM_NORMAL)
         iiatms = lenrel(iatms-1)+1
         lennew = lensec(ltri)
c
         call secget(isect(8),m8,iblok)
         iblok = iblok + mvadd
         call rdedx(q(icc),num*ncoorb,iblok,ifild)
c
         call vclr(q(ifa),1,ltri*nat3)
c
c        ibd   = iochf(15)
         ida_t = ida
         n = 0
         do i = 1, nat
            do j = 1, 3
               n = n + 1
               iq(iiatms+n-1) = i
c              call rdedx(q(ida_t),ltri,ibd,ifockf)
               call fetch(q(ida_t),ltri,'pdens',n)
               ida_t = ida_t + ltri
c              ibd   = ibd   + lennew
            enddo
         enddo
         call dscal(ltri*nat3,0.5d0,q(ida),1)
         do i = 1, nat3
            do j = 1, nocca
               q(ida-1+(i-1)*ltri+j*(j+1)/2)
     &         = 0.5d0*q(ida-1+(i-1)*ltri+j*(j+1)/2)
            enddo
         enddo
c
         ierror = CD_chf_dksm_mo(iq,q,nat3,iq(iiatms),ncoorb,nocca,0,
     &                           q(icc),q(inull),q(ida),q(inull),
     &                           q(ifa),q(inull),.false.,iwr)
c
c        ibt = iochf(16)
         do k = 1, nat3
c           call rdedx(q(icc),ltri,ibt,ifockf)
            call fetch(q(icc),ltri,'fock',k)
            ija = 0
            do i = 1, nsa4
               do j = 1, i
                  ij  = iky(mapie(i)) + mapie(j) - 1
                  q(icc+ij) = q(icc+ij) + q(ifa+(k-1)*ltri+ija)
                  ija = ija + 1
               enddo
            enddo
            call stash(q(icc),ltri,'fock',k)
c           call wrt3(q(icc),ltri,ibt,ifockf)
c           ibt = ibt + lennew
         enddo
c
         call gmem_free_inf(iatms,fnm,snm,'iatms')
         call gmem_free_inf(icc,fnm,snm,'alpha-vectors')
         call gmem_free_inf(ifa,fnm,snm,'pert-fock')
         call gmem_free_inf(ida,fnm,snm,'pert-density')
         call revise
         call clredx
      endif
      end
c
c
c
_IF(notused)
      subroutine rhscl_old(eps,bx,by,bz,eval,sx,sy,sz)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
c
c     r h s of chf equations ( nuclear displacements )
c
INCLUDE(common/mapper)
      common/maxlen/maxq
INCLUDE(common/common)
INCLUDE(common/cndx40)
INCLUDE(common/cndx41)
c
INCLUDE(common/infoa)
_IFN1(iv)      common/craypk/labs(1360)
_IF1(iv)      common/craypk/i205(340),j205(340),k205(340),l205(340)
_IFN1(iv)      common/small/labout(1360)
_IF1(iv)      common/small/labij(340),labkl(340)
      common/blkin/gin(510),nword
      common/out/gout(510),nw
INCLUDE(common/atmblk)
      dimension eps(*),bx(*),by(*),bz(*),eval(*),sx(*),sy(*),sz(*)
      data smal/1.0d-20/
      data zero,one,half,four/0.0d0,1.0d0,0.5d0,4.0d0/
c
c     sorts out scf eigenvalues and perturbation matrix elements
c
_IFN1(iv)      do 20 iiii = 1 , 1360
_IFN1(iv)         labout(iiii) = 0
_IFN1(iv) 20   continue
_IF1(iv)      call izero(680,labij,1)
      nplus1 = nocca + 1
      call secget(isect(9),9,isec9)
c
c     read in eigenvalues
c
      call rdedx(eval,lds(isect(9)),isec9,ifild)
      do 40 ii = 1 , nocca
         i1 = mapie(ii)
         do 30 jj = nplus1 , nsa4
            j1 = mapie(jj)
            it = (jj-nocca-1)*nocca + ii
            eps(it) = one/(eval(j1)-eval(i1))
 30      continue
 40   continue
c
c     sort out integrals involving one virtual m.o.
c
      i = 0
      j = 0
      k = 0
      l = 0
      ifili = nufile(1)
      ifilo = ifils

      write(6,*)'RHS ifili=',ifili,' ifilo=',ifilo

      ib1 = kblk(1)
      ib2 = nblk(1)
      nw = 0
      nat3 = nat*3
      iblll = lensec(mn)
      iblkb = iblks
      ib3 = iblkb + iblll*nat3
      call search(ib1,ifili)
      call wrt3z(iblkb,ifilo,ib3)
      call search(ib3,ifilo)
      do 70 ibl = ib1 , ib2
         call find(ifili)
         call get(gin,nn)
c
c     complete list of two-electron integrals coming
c     in from transformed mainfile.
c     those of form <aj/kl> going out
c     on scratchfile (ed7)
c
         if (nn.eq.0) go to 80
_IFN1(iv)         call unpack(gin(num2e+1),lab816,labs,numlab)
_IF1(iv)         call upak8v(gin(num2e+1),i205)
         do 60 n = 1 , nword
_IF(ibm,vax,fps,cyber205)
            i = i205(n)
            j = j205(n)
            k = k205(n)
            l = l205(n)
_ELSEIF(alliant,dec,ipsc,littleendian)
            n4 = (n+n) + (n+n)
            i = labs(n4-2)
            j = labs(n4-3)
            k = labs(n4  )
            l = labs(n4-1)
_ELSE
            n4 = (n+n) + (n+n)
            i = labs(n4-3)
            j = labs(n4-2)
            k = labs(n4-1)
            l = labs(n4)
_ENDIF
            if (i.gt.nocca) then
               if (j.le.nocca .and. k.le.nocca .and. l.le.nocca) then
                  nw = nw + 1
                  gout(nw) = gin(n)
_IFN1(iv)                  n4 = (nw+nw) + (nw+nw)
_IFN1(iv)                  labout(n4-3) = i
_IFN1(iv)                  labout(n4-2) = j
_IFN1(iv)                  labout(n4-1) = k
_IFN1(iv)                  labout(n4) = l
_IF1(iv)                  labij(nw) = j + i4096(i)
_IF1(iv)                  labkl(nw) = l + i4096(k)
                  if (nw.ge.340) then
c
c     writing out a block
c
_IFN1(iv)                     call pack(gout(num2e+1),lab816,labout,numlab)
_IF1(iv)                    call pak4v(labij,gout(num2e+1))
                     call put(gout,m511,ifilo)
_IFN1(iv)                     do 50 iiii = 1 , 1360
_IFN1(iv)                        labout(iiii) = 0
_IFN1(iv) 50                  continue
_IF1(iv)                    call izero(680,labij,1)
                     ib3 = ib3 + 1
                     nw = 0
                  end if
               end if
            end if
 60      continue
 70   continue
 80   if (nw.ne.0) then
_IFN1(iv)         call pack(gout(num2e+1),lab816,labout,numlab)
_IF1(iv)         call pak4v(labij,gout(num2e+1))
         call put(gout,m511,ifilo)
_IFN1(iv)         do 90 iiii = 1 , 1360
_IFN1(iv)            labout(iiii) = 0
_IFN1(iv) 90      continue
_IF1(iv)         call izero(680,labij,1)
         ib3 = ib3 + 1
      end if
      ib4 = ib3 - 1
c
      write(6,*)'rhs sorting complete'
c
c     sorting completed
c
c
      ib3 = iblkb + iblll*nat3
c
c      derivatives of overlap matrix section 14 of fockfile
c
      ibs = iochf(14)
      lennew = iky(ncoorb+1)
      newblk = lensec(lennew)
      np = nat3
      i = 0
      j = 0
      k = 0
      l = 0
      do 130 n = 1 , nat

         write(6,*)'atom',n
         ncount = 0

         do 100 jj = 1 , mn
            bx(jj) = zero
            by(jj) = zero
            bz(jj) = zero
 100     continue

         call fetch(sx, lennew, 'ds', (n-1)*3+1)
         call fetch(sy, lennew, 'ds', (n-1)*3+2)
         call fetch(sz, lennew, 'ds', (n-1)*3+3)

cc         call rdedx(sx,lennew,ibs,ifockf)
cc         call reads(sy,lennew,ifockf)
cc         call reads(sz,lennew,ifockf)

         ibs = ibs + newblk*3
c
c     have just read overlap derivatives sx,sy,sz for
c     one atom.  now scan list of integrals <aj/kl>
c
         call search(ib3,ifilo)
         do 120 ibl = ib3 , ib4
            call find(ifilo)
            call get(gin,nn)
_IFN1(iv)            call unpack(gin(num2e+1),lab816,labs,numlab)
_IF1(iv)            call upak8v(gin(num2e+1),i205)
            do 110 int = 1 , nword
_IFN1(iv)               n4 = (int+int) + (int+int)
_IFN1(iv)               i = labs(n4-3)
_IFN1(iv)               j = labs(n4-2)
_IFN1(iv)               k = labs(n4-1)
_IFN1(iv)               l = labs(n4)
_IF1(iv)               i = i205(int)
_IF1(iv)               j = j205(int)
_IF1(iv)               k = k205(int)
_IF1(iv)               l = l205(int)
                     ncount = ncount + 1
                     write(6,*)i,j,k,l,gin(int)
               gg = -gin(int)
               if (k.eq.l) gg = gg*half
               gg4 = gg*four
               ii = (i-nocca-1)*nocca
               nij = ii + j
               nik = ii + k
               nil = ii + l
               j1 = mapie(j)
               k1 = mapie(k)
               l1 = mapie(l)
               nkl = iky(k1) + l1
               njk = iky(j1) + k1
               njl = iky(j1) + l1
               if (k1.gt.j1) njk = iky(k1) + j1
               if (l1.gt.j1) njl = iky(l1) + j1
               bx(nij) = bx(nij) + gg4*sx(nkl)
               by(nij) = by(nij) + gg4*sy(nkl)
               bz(nij) = bz(nij) + gg4*sz(nkl)
               bx(nik) = bx(nik) - gg*sx(njl)
               by(nik) = by(nik) - gg*sy(njl)
               bz(nik) = bz(nik) - gg*sz(njl)
               bx(nil) = bx(nil) - gg*sx(njk)
               by(nil) = by(nil) - gg*sy(njk)
               bz(nil) = bz(nil) - gg*sz(njk)
 110        continue
 120     continue
c
c     have just formed that part of r.h.s. involving
c'    product of s' with <ij/kl>
c     store this on stratchfile
c
         write(6,*)'partial RHS x',n
         write(6,99)(bx(iii),iii=1,mn)
         write(6,*)'partial RHS y',n
         write(6,99)(by(iii),iii=1,mn)
         write(6,*)'partial RHS z',n
         write(6,99)(bz(iii),iii=1,mn)


         call wrt3(bx,mn,iblkb,ifils)
         call wrt3s(by,mn,ifils)
         call wrt3s(bz,mn,ifils)
         iblkb = iblkb + iblll*3

         write(6,*)'ncount',ncount
 130  continue
c
c
c
      write(6,*)'RHS - read perturbed fock'
      ibh = iochf(13)
      ibs = iochf(14)
      iblkb = iblks
      do 180 l = 1 , nat3
c
c     read perturbed fock matrix elements (integral contribution
c     only) from fockfile at section iochf(13)
c
cc         call rdedx(sx,lennew,ibh,ifockf)
         call fetch(sx, lennew, 'dh', l)
c
c     and the overlap derivatives again
c
cc         call rdedx(sy,lennew,ibs,ifockf)
         call fetch(sy, lennew, 'ds', l)

         call rdedx(bx,mn,iblkb,ifils)
         ibh = ibh + newblk
         ibs = ibs + newblk
         do 150 j = 1 , nocca
            j1 = mapie(j)
            do 140 i = nplus1 , nsa4
               i1 = mapie(i)
               it = iky(i1) + j1
               mt = (i-nocca-1)*nocca + j
               bx(mt) = bx(mt) + sx(it) - eval(j1)*sy(it)
 140        continue
 150     continue
         do 170 j = 1 , nocca
            do 160 i = nplus1 , nsa4
               mt = (i-nocca-1)*nocca + j
               if (dabs(bx(mt)).lt.smal) bx(mt) = zero
 160        continue
 170     continue
c
c
c     write complete r.h.s. to scratchfile
c
         write(6,*)'RHS ',l
         write(6,99)(bx(iii),iii=1,mn)
 99      format(1x,8f10.5)

         call wrt3(bx,mn,iblkb,ifils)
         iblkb = iblkb + iblll
 180  continue
      call clredx

      write(6,*)'RHSides complete'

      return
      end
_ENDIF(notused)
c
c  RHS - version reading integrals from
c        global arrays and 1e matrices using
c        stash/fetch
c
      subroutine rhscl_ga(eps,bx,by,bz,eval,sx,sy,sz,buf)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/mapper)
      common/maxlen/maxq
INCLUDE(common/common)
INCLUDE(common/cndx40)
INCLUDE(common/cndx41)
INCLUDE(common/infoa)
INCLUDE(common/atmblk)
INCLUDE(common/global)
_IF(ccpdft)
INCLUDE(common/ccpdft.hf77)
_ENDIF

      dimension eps(*),bx(*),by(*),bz(*),eval(*),sx(*),sy(*),sz(*)
      data smal/1.0d-20/
      data zero,one,half,four/0.0d0,1.0d0,0.5d0,4.0d0/

      REAL buf(*)
c
_IF(ccpdft)
      hf_wght = CD_HF_exchange_weight()
_ENDIF
c
c     sorts out scf eigenvalues and perturbation matrix elements
c
      nplus1 = nocca + 1
      nv = nsa4 - nocca
      no = nocca
ccc      write(6,*)'rhscl_ga: orbs',no,nv
c
c     read in eigenvalues
c
      call secget(isect(9),9,isec9)
      call rdedx(eval,lds(isect(9)),isec9,ifild)
      do 40 ii = 1 , nocca
         i1 = mapie(ii)
         do 30 jj = nplus1 , nsa4
            j1 = mapie(jj)
            it = (jj-nocca-1)*nocca + ii
            eps(it) = one/(eval(j1)-eval(i1))
 30      continue
 40   continue
c
c
c      derivatives of overlap matrix section 14 of fockfile
c
cc      ibs = iochf(14)
      lennew = iky(ncoorb+1)
cc      newblk = lensec(lennew)
      nat3 = nat*3
      np = nat3
cc      iblll = lensec(mn)
cc      iblkb = iblks

      do 130 n = 1 , nat
         do 100 jj = 1 , mn
            bx(jj) = zero
            by(jj) = zero
            bz(jj) = zero
 100     continue

         call fetch(sx, lennew, 'ds', (n-1)*3+1)
         call fetch(sy, lennew, 'ds', (n-1)*3+2)
         call fetch(sz, lennew, 'ds', (n-1)*3+3)
c
c     have just read overlap derivatives sx,sy,sz for
c     one atom.  now scan list of integrals <aj/kl>
c
c  vooo
         jkl = 0
         do k = 1 , no
            do l = 1 , k
               do j = 1 , no
                  jkl = jkl + 1
                  if (jkl.ge.jl_vooo.and.jkl.le.jh_vooo) then
                     call pg_get(g_vooo,1,nv,jkl,jkl,buf,1)
                     do iv = 1 , nv
                        i = iv + no
c
c nb k,l are triangulated
c
c                        write(6,*)i,j,k,l,buf(iv)
                        gg = -buf(iv)
                        if (k.eq.l) gg = gg*half
                        gg4 = gg*four

                        ii = (i-nocca-1)*nocca
                        nij = ii + j
                        nik = ii + k
                        nil = ii + l

                        j1 = mapie(j)
                        k1 = mapie(k)
                        l1 = mapie(l)

                        nkl = iky(k1) + l1
                        njk = iky(j1) + k1
                        njl = iky(j1) + l1

                        if (k1.gt.j1) njk = iky(k1) + j1
                        if (l1.gt.j1) njl = iky(l1) + j1
                        bx(nij) = bx(nij) + gg4*sx(nkl)
                        by(nij) = by(nij) + gg4*sy(nkl)
                        bz(nij) = bz(nij) + gg4*sz(nkl)
_IF(ccpdft)
                        gg = gg*hf_wght
_ENDIF
                        bx(nik) = bx(nik) - gg*sx(njl)
                        by(nik) = by(nik) - gg*sy(njl)
                        bz(nik) = bz(nik) - gg*sz(njl)
                        bx(nil) = bx(nil) - gg*sx(njk)
                        by(nil) = by(nil) - gg*sy(njk)
                        bz(nil) = bz(nil) - gg*sz(njk)
                     end do
                  end if
               end do
            end do
         end do

         call pg_dgop(1004,bx,mn,'+')
         call pg_dgop(1004,by,mn,'+')
         call pg_dgop(1004,bz,mn,'+')
c
c     have just formed that part of r.h.s. involving
c'    product of s' with <ij/kl>
c     store this on scratchfile
c
         if(.false.)then
         write(6,*)'partial RHS x',n
         write(6,99)(bx(iii),iii=1,mn)
         write(6,*)'partial RHS y',n
         write(6,99)(by(iii),iii=1,mn)
         write(6,*)'partial RHS z',n
         write(6,99)(bz(iii),iii=1,mn)
         endif


         call stash(bx,mn,'rhs',(n-1)*3+1)
         call stash(by,mn,'rhs',(n-1)*3+2)
         call stash(bz,mn,'rhs',(n-1)*3+3)

cc         call wrt3(bx,mn,iblkb,ifils)
cc         call wrt3s(by,mn,ifils)
cc         call wrt3s(bz,mn,ifils)
cc         iblkb = iblkb + iblll*3

 130  continue
c
c      write(6,*)'RHS - read perturbed fock'

cc      iblkb = iblks
      do 180 l = 1 , nat3
c
c     read perturbed fock matrix elements (integral contribution
c     only) from fockfile at section iochf(13)
c
cc         call rdedx(sx,lennew,ibh,ifockf)
         call fetch(sx, lennew, 'dh', l)
c
c     and the overlap derivatives again
c
cc         call rdedx(sy,lennew,ibs,ifockf)
         call fetch(sy, lennew, 'ds', l)
c
cc         call rdedx(bx,mn,iblkb,ifils)
         call fetch(bx,mn,'rhs',l)
c
cc        ibh = ibh + newblk
cc         ibs = ibs + newblk
         do 150 j = 1 , nocca
            j1 = mapie(j)
            do 140 i = nplus1 , nsa4
               i1 = mapie(i)
               it = iky(i1) + j1
               mt = (i-nocca-1)*nocca + j
               bx(mt) = bx(mt) + sx(it) - eval(j1)*sy(it)
 140        continue
 150     continue
         do 170 j = 1 , nocca
            do 160 i = nplus1 , nsa4
               mt = (i-nocca-1)*nocca + j
               if (dabs(bx(mt)).lt.smal) bx(mt) = zero
 160        continue
 170     continue
c
c     write complete r.h.s. to scratchfile
c
cdbg
c         write(6,*)'RHS ',l
c         write(6,99)(bx(iii),iii=1,mn)
 99      format(1x,8f10.5)

         call stash(bx,mn,'rhs',l)
cc         call wrt3(bx,mn,iblkb,ifils)
cc         iblkb = iblkb + iblll
 180  continue
      call clredx

c      write(6,*)'RHSides complete'

      return
      end
c
      subroutine rhscl_dft(q,iq)
      implicit none
c
c     Adds the DFT contributions onto the right-hand-sides.
c
c     Parameters:
c
INCLUDE(common/sizes)
      integer m8
      parameter(m8=8)
c
c     Input:
c
INCLUDE(common/infoa)
INCLUDE(common/cndx41)
INCLUDE(common/atmblk)
INCLUDE(common/common)
INCLUDE(common/mapper)
INCLUDE(common/iofile)
c
c     Workspace:
c
      REAL q(*)
      integer iq(*)
c
c     Functions:
c
INCLUDE(common/ccpdft.hf77)
INCLUDE(common/gmempara)
      integer igmem_alloc_inf, igmem_null
      integer lensec, lenrel
c
c     Local:
c
      integer nat3,lennew,newblk,mnblk,is,ib,ic,iblok,iatms,iiatms
      integer ibs,ics,icb,i,j,ierror,iblkb
      integer inull
      character *15 fnm
      character *9  snm
      data fnm/'secd_parallel.m'/
      data snm/'rhscl_dft'/
_IF(ccpdft)
      if (CD_active()) then
         inull = igmem_null()
         nat3 = nat*3
         lennew = iky(ncoorb)+ncoorb
         newblk = lensec(lennew)
         mnblk  = lensec(mn)
         is = igmem_alloc_inf(lennew*nat3,fnm,snm,
     &                        'pert-overlap',IGMEM_NORMAL)
         ib = igmem_alloc_inf(mn*nat3,fnm,snm,
     &                        'right-hand-sides',IGMEM_NORMAL)
         ic = igmem_alloc_inf(num*num,fnm,snm,
     &                        'alpha-vectors',IGMEM_NORMAL)
c
c        Get the MO-coefficients
c
         call secget(isect(8),m8,iblok)
         iblok = iblok + mvadd
         call rdedx(q(ic),num*ncoorb,iblok,ifild)
c
c        Get the derivative overlap matrices
c
c        ibs = iochf(14)
         ics = is
         do i = 1, nat3
c           call rdedx(q(ics),lennew,ibs,ifockf)
            call fetch(q(ics),lennew,'ds',i)
            ics = ics + lennew
c           ibs = ibs + newblk
         enddo
c
c        Calculate the RHS contributions
c
         iatms = igmem_alloc_inf(nat3,fnm,snm,
     &                           'perturbations',IGMEM_DEBUG)
         iiatms = lenrel(iatms-1)+1
         do i = 1, nat
            iq(iiatms+3*(i-1)+0) = i
            iq(iiatms+3*(i-1)+1) = i
            iq(iiatms+3*(i-1)+2) = i
         enddo
         call vclr(q(ib),1,mn*nat3)
         ierror = CD_chf_rhs_mo(iq,q,nat3,iq(iiatms),ncoorb,nocca,0,
     &                          q(ic),q(inull),q(is),q(inull),
     &                          q(ib),q(inull),.false.,iwr)
         call gmem_free_inf(iatms,fnm,snm,
     &                      'perturbations')
c
c        Add the DFT RHS contributions onto the Hartree-Fock parts
c
c        iblkb = iblks
         icb = ib
         do i = 1, nat3
c           call rdedx(q(ic),mn,iblkb,ifils)
            call fetch(q(ic),mn,'rhs',i)
            do j = 0, mn-1
               q(ic+j)=q(ic+j)+q(icb+j)
            enddo
c           call wrt3(q(ic),mn,iblkb,ifils)
c           iblkb = iblkb + mnblk
            call stash(q(ic),mn,'rhs',i)
            icb   = icb   + mn
         enddo
         call clredx
         call gmem_free_inf(ic,fnm,snm,'alpha-vectors')
         call gmem_free_inf(ib,fnm,snm,'right-hand-sides')
         call gmem_free_inf(is,fnm,snm,'pert-overlap')
      endif
_ENDIF
      end
c
      subroutine ovlcl(qq,iqq)
c
c   assemble overlap contribution to closed shell scf second
c   derivatives
c   term involving derivative of lagrangian
c   closed shell case
c
      implicit REAL  (a-h,o-z)
      logical out
INCLUDE(common/sizes)
INCLUDE(common/infoa)
INCLUDE(common/nshel)
      common/bufb/e(maxorb),e1(maxorb)
      dimension qq(*),iqq(*)
      common/maxlen/maxq
c
INCLUDE(common/common)
INCLUDE(common/iofile)
INCLUDE(common/cndx41)
INCLUDE(common/mapper)
INCLUDE(common/symtry)
c
INCLUDE(common/prnprn)
c
      out = odebug(6)
      call secget(isect(9),9,isec9)
      call rdedx(e,lds(isect(9)),isec9,ifild)
      nat3 = nat*3
      nlen = nat3*nat3
      iof = lenrel(nw196(5))
      iofs = nw196(5) + lenint(nat*nt)
      i10 = iofs + 1
      i11 = i10 + nlen
      ioff = i11 + nx
      noc1 = na + 1
c
c     perturbed fock matrix at section 16 , m.o. basis
c     derivatives of overlap at section 14
c
cc      ibf = iochf(16)
cc      ibs = iochf(14)
      ltri = ikyp(ncoorb)
cc      length = lensec(ltri)
      if (odebug(6)) then
         write (iwr,6010)
         do 30 n = 1 , nat3

cc            call rdedx(qq(i11),ltri,ibf,ifockf)
cc            call rdedx(qq(ioff),ltri,ibs,ifockf)

            call fetch(qq(i11),ltri,'fock',n)
            call fetch(qq(ioff),ltri,'ds',n)

cc            ibf = ibf + length
cc            ibs = ibs + length
            do 20 i = 1 , ncoorb
               ii = ikyp(i) - 1
c
c     perturbed eigenvalues ( used as check only)
c
               e1(i) = qq(i11+ii) - e(i)*qq(ioff+ii)
 20         continue
            write (iwr,6020) (e1(i),i=1,ncoorb)
 30      continue
      end if
cc      ibf = iochf(16)
c
c     perturbed density matrix at section 15
c
cc      ibd = iochf(15)
      do 90 n = 1 , nat3
         n10 = i10 - 1 + n
c
c     complete derivative of fock matrix (only elements
c     with two occupied orbitals are needed)
c
cc         call rdedx(qq(i11),ltri,ibf,ifockf)
         call fetch(qq(i11),ltri,'fock',n)
c
c     perturbed density matrix
c
cc         call rdedx(qq(ioff),ltri,ibd,ifockf)
         call fetch(qq(ioff),ltri,'pdens',n)
c
c     derivative lagrangian elements for two occupied orbitals
c
         do 50 i = 1 , na
            ii = iky(i) - 1
            do 40 j = 1 , i
               ij = ii + j
               qq(ioff+ij) = qq(ioff+ij)*(e(i)+e(j)) + qq(i11+ij)
     +                       + qq(i11+ij)
 40         continue
 50      continue
c
c     derivative lagrangian elements for one occupied and one virtual m.
c
         do 70 i = noc1 , ncoorb
            ii = iky(i) - 1
            do 60 j = 1 , na
               ij = ii + j
               qq(ioff+ij) = qq(ioff+ij)*e(j)
 60         continue
 70      continue
c
c
cc         ibf = ibf + length
cc         ibd = ibd + length
cc         ibs = iochf(14)
         do 80 m = 1 , nat3
c
c     derivative overlap matrix
c
cc            call rdedx(qq(i11),ltri,ibs,ifockf)
cc            ibs = ibs + length
            call fetch(qq(i11),ltri,'ds',m)
c
c     take product with overlap derivatives
c
            qq(n10+(m-1)*nat3) = -tracep(qq(i11),qq(ioff),ncoorb)
 80      continue
 90   continue
      call rdedx(qq(1),nw196(5),ibl196(5),ifild)
      if (out) then
         call dr2sym(qq(i10),qq(i11),iqq(1),iqq(iof+1),nat,nat3,
     +               nshell)
         write (iwr,6030)
         call prnder(qq(i10),nat3,iwr)
      end if
c
c#accum derivative of lagrangian
c
      call secget(isect(60),60,isec46)
      call rdedx(qq(i11),nlen,isec46,ifild)
      call vadd(qq(i11),1,qq(i10),1,qq(i11),1,nlen)
      call wrt3(qq(i11),nlen,isec46,ifild)
      return
 6010 format (//5x,'perturbed eigenvalues'//)
 6020 format (//(5x,6f16.8))
 6030 format (//' contribution from derivative of lagrangian')
      end
c
c  GA-based fock builder - replaces pfock/sgmatm
c
      subroutine pfock_ga(fock,p,nfok,ltri,buff)
      implicit none

INCLUDE(common/global)
INCLUDE(common/timeperiods)
INCLUDE(common/cndx41)

      integer i, j, k, l
      integer ij, jkl, kl, ik, jl, kk, ii
      integer iv, noo2, nvv2
      integer no, noo, nov, nv, nvv

      REAL p(*), fock(*), buff(*)
      REAL thresh
      integer ltri, nfok

      data thresh /1.0d-10/
      no = nocca
      nv = nvirta
      noo = no*no
      nvv = nv*nv
      nov = no*nv
      noo2 = no*(no+1)/2
      nvv2 = nv*(nv+1)/2
c
c oooo
c
c      write(6,*)'pfock_ga nfok=',nfok

      call start_time_period(TP_2D_PFOCK_OOOO)

      ij=1
      do i = 1, no
         do j = 1, i
            if(ij .ge. jl_oooo .and. ij .le. jh_oooo)then
c
               call pg_get(g_oooo,1,noo2,ij,ij,buff,1)
               kl = 1
               do k = 1 , no
                  do l = 1 , k
c
c NB GA contains more integrals of oooo class than ed6,
c discard some
c
                    if (ij .ge. kl.and.
     +              abs(buff(kl)).gt.thresh) then
                        call pfock_bld(i,j,k,l,
     &                       buff(kl),
     &                       p,fock,nfok,ltri)
                        
                    endif
                     kl = kl + 1
                  end do
               end do
            endif
            ij = ij + 1
         end do
      end do

      call end_time_period(TP_2D_PFOCK_OOOO)

      call start_time_period(TP_2D_PFOCK_VOOO)

c
c vooo
c
      jkl = 0
      do k = 1 , no
         do l = 1 , k
            do j = 1 , no
               jkl = jkl + 1
               if (jkl.ge.jl_vooo.and.jkl.le.jh_vooo) then
                  call pg_get(g_vooo,1,nv,jkl,jkl,buff,1)
                  do iv = 1 , nv
                     i = iv + no
                     if(abs(buff(iv)).gt.thresh) then
                     call pfock_bld(i,j,k,l,
     &                    buff(iv),
     &                    p,fock,nfok,ltri)
                     endif
                  enddo
               endif
            enddo
         enddo
      enddo
      call end_time_period(TP_2D_PFOCK_VOOO)

      call start_time_period(TP_2D_PFOCK_VVOO)

c
c vvoo
c
c  oo triangle of vv triangles
c
      kl = 1
      do k = 1, no
         do l = 1, k
            if(kl .ge. jl_vvoo .and. kl .le. jh_vvoo)then
               call pg_get(g_vvoo,1,nvv2,kl,kl,buff,1)
               ij = 1
               do i = 1, nv
                  do j = 1, i

                     if(abs(buff(ij)).gt.thresh) then
                     call pfock_bld(no+i,no+j,k,l,
     &                    buff(ij),
     &                    p,fock,nfok,ltri)
                     endif
                     ij = ij + 1
                  enddo
               enddo
            endif
            kl = kl + 1
         enddo
      enddo

      call end_time_period(TP_2D_PFOCK_VVOO)
      call start_time_period(TP_2D_PFOCK_VOVO)
c
c vovo class
c
      jl = 1
      do j = 1, no
         do l = 1, j
            if(jl .ge. jl_vovo .and. jl .le. jh_vovo)then
               call pg_get(g_vovo,1,nvv,jl,jl,buff,1)
               ik = 1
               do kk = 1, nv
                  do ii = 1, nv

                     k = kk + no
                     i = ii + no

                     ij = j + (i*(i-1))/2
                     kl = l + (k*(k-1))/2

                     if(j.eq.l .and. k .gt. i) then
c
c  <ij|kj>  integrals are duplicated
c                        
                     else if(ij .ge. kl)then
c
c already canonical ordering
c
cc                        write(6,*)'ijkl',i,j, k, l, buff(ik)
                        if(abs(buff(ik)).gt.thresh) then
                        call pfock_bld(i,j, k,l,
     &                       buff(ik),
     &                       p,fock,nfok,ltri)
                        endif
                     else
c
c switch ij / kl
c
cc                        write(6,*)'ijkl',k,l, i, j, buff(ik)
                        if(abs(buff(ik)).gt.thresh) then
                        call pfock_bld(k,l,i,j,
     &                       buff(ik),
     &                       p,fock,nfok,ltri)
                        endif
                     endif
                     ik = ik + 1
                  enddo
               enddo
            endif
            jl = jl + 1
         enddo
      enddo
      call end_time_period(TP_2D_PFOCK_VOVO)

      call start_time_period(TP_2D_PFOCK_SUM)
      call pg_dgop(1001,fock,nfok*ltri,'+')
      call end_time_period(TP_2D_PFOCK_SUM)

cdbg
c      do i = 1,nfok
c         write(6,*)'fock (ga)..',i
c         call prtri(fock( (i-1)*ltri + 1), ncoorb)
c      enddo

      return

      end
      subroutine pfock_bld(i,j,k,l,g,p,fock,nfok,ltri)

      implicit REAL  (a-h,o-z)      

INCLUDE(common/sizes)
INCLUDE(common/mapper)
_IF(ccpdft)
INCLUDE(common/ccpdft.hf77)
      logical ohf_exch
_ENDIF

      REAL p(*), fock(*)

_IF(ccpdft)
      hf_wght  = CD_HF_exchange_weight()
      ohf_exch = (.not.CD_active()).or.CD_HF_exchange()
_ENDIF
      gik = g
      g2 = gik + gik
      g4 = g2 + g2
      ikyi = iky(i)
      ikyj = iky(j)
      ikyk = iky(k)
      ik = ikyi + k
      il = ikyi + l
      ij = ikyi + j
      jk = ikyj + k
      jl = ikyj + l
      kl = ikyk + l
      ioff = 0

_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
      do 20 n = 1 , nfok
         aij = g4*p(kl+ioff) + fock(ij+ioff)
         fock(kl+ioff) = g4*p(ij+ioff) + fock(kl+ioff)
         fock(ij+ioff) = aij

cc         write(6,*)'wrt',kl+ioff, g4*p(ij+ioff) + fock(kl+ioff)
cc         write(6,*)'wrt',ij+ioff,aij

         ioff = ioff + ltri
 20   continue

c... exchange
_IF(ccpdft)
      if (.not.ohf_exch) then
         return
      else
         gik = hf_wght*gik
         g2  = hf_wght*g2
      endif
_ENDIF
      gil = gik
      if (i.eq.k .or. j.eq.l) gik = g2
      if (j.eq.k) gil = g2
      if (j.lt.k) then
         jk = ikyk + j
         if (j.lt.l) then
            jl = iky(l) + j
         end if
      end if
      ioff = 0

_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
      do 30 n = 1 , nfok
         ajk = fock(jk+ioff) - gil*p(il+ioff)
         ail = fock(il+ioff) - gil*p(jk+ioff)
         aik = fock(ik+ioff) - gik*p(jl+ioff)
         fock(jl+ioff) = fock(jl+ioff) - gik*p(ik+ioff)
         fock(jk+ioff) = ajk
         fock(il+ioff) = ail
         fock(ik+ioff) = aik
         ioff = ioff + ltri
 30   continue

      end
c
c  chf solver - from GA hessian  temporary in-core solution
c
      subroutine chfeqs_ga(eps,b,cc,wks1,wks2,alpha,aa,
     +     ndim,skipp)
      implicit REAL  (a-h,o-z)
c
c     simultaneous equations - small case
c
      logical skipp
      dimension skipp(100)
c
      dimension eps(ndim),b(ndim),cc(ndim),wks1(ndim),wks2(ndim),
     &    alpha(ndim,ndim),aa(ndim,ndim)
INCLUDE(common/cigrad)
_IFN1(iv)      common/craypk/labs(1360)
_IF1(iv)      common/craypk/ii(340),jj(340)
INCLUDE(common/prnprn)
      common/blkin/g(510),nword
c
INCLUDE(common/common)
INCLUDE(common/iofile)
INCLUDE(common/cndx41)
INCLUDE(common/global)
c
      data zero,one/0.0d0,1.0d0/

      integer o1, o2, v1, v2
c
cc      lenblk = lensec(mn)
c
c     rhs of equations input from scratchfile
c
cc      iblkb = iblks
cc      iblku = iblkb + lenblk*np

      if (oprn(12)) write (iwr,6010)
      if (oprn(13)) write (iwr,6020)
      ifail = 0
c
      do 30 i = 1 , mn
         do 20 j = 1 , mn
            alpha(j,i) = zero
 20      continue
 30   continue

      no = nocca
      nv = nvirta
      nvo = nv*no
      nvv = nv*nv

      ij = 0

      do o1=1,no
         do o2=1,o1

            ij = ij + 1

            call pg_get(g_hess,1,nvv,ij,ij,aa,nvv)

            do v1 = 1,nv
               do v2 = 1,nv

                  ix1 = (v1-1)*no + o1
                  ix2 = (v2-1)*no + o2

                  ix = (ix1-1) * nvo + ix2
                  alpha(ix1,ix2) = aa((v1-1)*nv+v2,1)

                  ix = (ix2-1) * nvo + ix1
                  alpha(ix2,ix1) = aa((v1-1)*nv+v2,1)

               enddo
            enddo
         enddo
      enddo

      do 40 i = 1 , mn
         alpha(i,i) = alpha(i,i) + one/eps(i)
 40   continue

c
c     loop over the perturbations
c
      if (odebug(2)) write (iwr,6030)
      if (odebug(2)) call prsqm(alpha,mn,mn,mn,iwr)

      do 80 i = 1 , mn
         do 70 j = 1 , mn
            alpha(i,j) = alpha(i,j)*eps(i)
 70      continue
 80   continue
c
      do 110 ko = 1 , np
c
c     get rhs
c
         call fetch(b,mn,'rhs',ko)
cc         call rdedx(b,mn,iblkb,ifils)

c      write(iwr,978)(b(i),i=1,mn)
c978   format(' b in chfeqs = ',5f15.10)

         if (oprn(12)) write (iwr,6040) ko , (b(i),i=1,mn)
cc         iblkb = iblkb + lenblk
         do 90 i = 1 , mn
            cc(i) = 0.0d0
 90      continue
         if (odebug(2).and.skipp(ko)) write (iwr,6050) ko
         if (.not.(skipp(ko))) then
            do 100 i = 1 , mn
               b(i) = -b(i)*eps(i)
 100        continue
c
c    use nag routine to solve
c
            call f04atf(alpha,mn,b,mn,cc,aa,mn,wks1,wks2,ifail)
         end if
c
c    write solution onto scratchfile
c
cc         call wrt3(cc,mn,iblku,ifils)
cc         iblku = iblku + lenblk

         call stash(cc,mn,'soln',ko)

         if (oprn(12) .or. oprn(13)) then
            write (iwr,6060) ko , (cc(i),i=1,mn)
         end if
 110  continue
      return
 6010 format (//1x,'print right-hand-side of chf equations')
 6020 format (//1x,'print solution to chf equations')
 6030 format (//1x,'a-matrix routine chfeqs')
 6040 format (//1x,'perturbation  ',i4//(5x,5f16.8))
 6050 format (/1x,'perturbation',i5,' omitted')
 6060 format (//1x,'solution  ',i4//(5x,5f16.8))
      end
c**** stash.M
c
c  stash and fetch routines - designed as a generic
c  interface for storing/retrieving matrices that
c  will be implemented via disk, local memory and GA
c  depending on architecture
c
c
c  simple file-based version for testing
c
c
      subroutine stash(a,n,tag,ix)
      implicit none
      integer n, ix
      REAL a
      dimension a(n)
      INCLUDE(common/vcore)
INCLUDE(common/secd_pointers)
INCLUDE(common/iofile)
      integer iloc, locatc,ioffs
      integer ipg_nodeid
      integer ilo,ihi,jlo,jhi,indx,iofft
      character *5 stashed(6)
      data stashed /
     +  'dh', 'ds', 'rhs', 'fock', 'pdens', 'soln' /
      character tag*(*), filename*30
c
      if (omem2nd) then
c
c        locate index of "stashed" datan and store into memory
c
         iloc = locatc(stashed,6,tag)
         if (iloc.le.0 .or. iloc.gt.6) then
            call caserr('invalid matrix tag in stash')
         endif
c
c        define memory location
c
         ioffs = (ix-1) * n
         iofft =  ix    * n
         if (omem2nd_ga) then
c
c           invoke pg_put such that only local copying is performed:
c
c           step 1) workout which part of the data is held locally 
c           on this processor.
c
            call pg_distribution(iof2nd(iloc),ipg_nodeid(),
     +      ilo,ihi,jlo,jhi)
            if (jlo.ne.1.or.jhi.ne.1) call caserr('stash distribution?')
            jlo = max(ioffs+1,ilo)
            jhi = min(iofft  ,ihi)
            indx = jlo-ioffs
c
c           step 2) if the size of the local part is greater than zero
c               call pg_put to store it.
c

            if (jhi-jlo.ge.0) then
                call pg_put(iof2nd(iloc),jlo,jhi,1,1,a(indx),1)
            endif
c
c           step 3) synchronise for sanity
c
            call pg_synch(103)
         else
            ioffs = ioffs + iof2nd(iloc)
            call dcopy(n,a,1,Q(ioffs),1)
         endif
c
      else
      
         call tagname(tag,len(tag),ix,filename)

         open(unit=20,file=filename,status='unknown',
     +     form='unformatted')
         write(20)n
         write(20)a
c
      endif
      return
      end

      subroutine fetch(a,n,tag,ix)
      implicit none
      integer n, ix, ntmp
      REAL a
      INCLUDE(common/vcore)
INCLUDE(common/secd_pointers)
INCLUDE(common/iofile)
      integer iloc, locatc, ioffs
      character *5 stashed(6)
      data stashed /
     +  'dh', 'ds', 'rhs', 'fock', 'pdens', 'soln' /
      character tag*(*), filename*30
      integer ipg_nodeid, ipg_nnodes
      integer ilo, ihi, jlo, jhi, n0, ip1, n1, n2, iofft
      dimension a(n)
c
      if (omem2nd) then
c
c     locate index of "stashed" datan and retrieve from memory
c
      iloc = locatc(stashed,6,tag)
      if (iloc.le.0 .or. iloc.gt.6) then
       call caserr('invalid matrix tag in fetch')
      endif
c
c     define memory location
c
      ioffs = (ix-1) * n
      iofft =  ix    * n
      if (omem2nd_ga) then
c
c       invoke pg_get such that as most local copying a possible is
c       done:
c
c       assumptions: the data is stored as a linear array of which each 
c                    processor holds the same sized fraction except
c                    for the last processor.
c
c       step 1) find the processor which holds the first element of
c               the data we are interested in.
c
        call pg_distribution(iof2nd(iloc),0,ilo,ihi,jlo,jhi)
        if (jlo.ne.1.or.jhi.ne.1) call caserr("fetch distribution?")
        n0 = ihi-ilo+1
        ip1 = (ioffs+1)/n0
c
c       step 2) check whether the next processor holds a larger portion
c               of the data we are interested than the processor which
c               holds the first element.
c
        if (ip1.lt.ipg_nnodes()-1) then
          call pg_distribution(iof2nd(iloc),ip1,ilo,ihi,jlo,jhi)
          n1 = ihi-ilo+1
          call pg_distribution(iof2nd(iloc),ip1+1,ilo,ihi,jlo,jhi)
          n2 = ihi-ilo+1
          if (n2.gt.n1) ip1 = ip1+1
        endif
c
c       step 3) ip1 is now the processor which holds the largest 
c               fraction of the array. If I am on that processor call
c               pg_get.
c
        if (ip1.eq.ipg_nodeid()) then
           call pg_get(iof2nd(iloc),ioffs+1,iofft,1,1,a,1)
        endif
c
c       now simply broadcast the data to everyone.
c
        call pg_brdcst(22,a,8*n,ip1)
      else
        ioffs = ioffs + iof2nd(iloc)
        call dcopy(n,Q(ioffs),1,a,1)
      endif
c
      else
      
      call tagname(tag,len(tag),ix,filename)
cdbg
c      write(6,*)'fetch: ',filename

      open(unit=20,file=filename, status='old',
     &     form='unformatted',err=80)

      read(20,err=90)ntmp
      if(ntmp .ne. n)call caserr('bad length in fetch')
      read(20,err=90)a
      close(20)
      return
 80   continue
      write(6,*)'problem with file=',tag
      call caserr('bad open in fetch')
 90   continue
      call caserr('bad read in fetch')
      endif
      end

      subroutine tagname(tag,len,ix,file)
      implicit none
      character tag*(*), file*30
      character znum*4, znod*4
      integer ix, i, j, len
      integer ipg_nodeid	
      write(znum,'(i4)') ix
      write(znod,'(i4)') ipg_nodeid()

      i = 1
      j = 1
      file = ' '
      do while (tag(i:i).ne.' ' .and. i .le. len)
         file(i:i) = tag(i:i)
         i = i + 1
         j = j + 1
      enddo

      do i= 1, 4
         if(znum(i:i) .eq. ' ')then
            file(j:j) = '0'
         else
            file(j:j) = znum(i:i)
         endif
         j = j + 1
      enddo

      do i= 2, 4
         if(znod(i:i) .eq. ' ')then
            file(j:j) = '0'
         else
            file(j:j) = znod(i:i)
         endif
         j = j + 1
      enddo

      return

      end
c
c  =============== slightly modified routines
c
      subroutine pdens(ss,lstop)
c
c    perturbed density matrices (nuclear motions) in m.o. basis
c
      implicit REAL  (a-h,o-z)
      logical lstop
INCLUDE(common/sizes)
INCLUDE(common/common)
INCLUDE(common/iofile)
INCLUDE(common/cndx40)
INCLUDE(common/cndx41)
      common/mpshl/ns(maxorb)
INCLUDE(common/ghfblk)
INCLUDE(common/infoa)
      dimension ss(*)
      character *8 grhf,oscf
INCLUDE(common/mapper)
INCLUDE(common/prnprn)
      logical open
      data grhf/'grhf'/
      data oscf/'oscf'/
      data zero/0.0d0/
      open = scftyp.eq.oscf
c
c
      if (odebug(3)) write (iwr,6010)
cc      iblll = lensec(mn)
cc      iblku = iblks + iblll*np
      if (lstop) then
c
c     dump because running out of time
c
         iochf(15) = iochf(1)
cc         iposd = iochf(1)
cc         do 20 n = 1 , npfin
cc            call rdedx(ss,mn,iblku,ifils)
cc            call wrt3(ss,mn,iposd,ifockf)
cc            iblku = iblku + iblll
cc            iposd = iposd + iblll
cc 20      continue
         npstar = npfin
         irest = 1
         return
      else
c
c
         nsoc = noccb - nocca
         ndpls1 = nocca + 1
         ntpls1 = noccb + 1
         nw = iky(ncoorb+1)
c
c
cc         iposs = iochf(14)
cc         iblls = lensec(nw)
c
c     output to section 15 of fockfile.
c     overlap matrix derivatives (m.o. basis) on section 14
c     perturbed wavefunctions on scratchfile
c
cc         iochf(15) = iochf(1)
cc         iposd = iochf(15)
cc         iposdo = iposd + 3*nat*iblls
cc
cc - write zeros onto section 15??
cc
cc         call wrt3z(iposd,ifockf,3*nat*iblls)
cc         if (open) call wrt3z(iposd,ifockf,3*nat*iblls)

         nu = nx + nx
         nopen = nu + mn
         do 170 i = 1 , nat
            do 160 j = 1 , 3

cc  overlap derivatives
cc               call rdedx(ss,nw,iposs,ifockf)
               call fetch(ss,nw,'ds',(i-1)*3+j)

cc cphf solutions
cc               call rdedx(ss(nu+1),mn,iblku,ifils)
               call fetch(ss(nu+1),mn,'soln',(i-1)*3+j)

cc               iposs = iposs + iblls
cc               iblku = iblku + iblll
               do 30 k = 1 , nw
                  if (open) ss(nopen+k) = zero
                  ss(nx+k) = zero
 30            continue
               if (scftyp.eq.grhf) then
                  nr = 0
                  kl = 0
                  do 50 k = 1 , nsa4
                     do 40 l = 1 , k
                        kl = iky(mapie(k)) + mapie(l)
                        ss(nx+kl) = -0.5d0*ss(kl)*
     +                              (fjk(ns(k))+fjk(ns(l)))
                        if (ns(k).ne.ns(l)) then
                           nr = nr + 1
                           ssss = ss(nx+kl)
                           ss(nx+kl) = ss(nx+kl) + 0.5d0*ss(nu+nr)
     +                                  *(fjk(ns(l))-fjk(ns(k)))
                        end if
 40                  continue
 50               continue
               else
                  if (nocca.ne.0) then
                     do 70 k = 1 , nocca
                        do 60 l = 1 , k
                           kl = iky(mapie(k)) + mapie(l)
                           ss(nx+kl) = -2.0d0*ss(kl)
 60                     continue
 70                  continue
                  end if
                  if (nsoc.ne.0 .and. nocca.ne.0) then
                     do 90 k = ndpls1 , noccb
                        do 80 l = 1 , nocca
                           kl = iky(mapie(k)) + mapie(l)
                           mt = (k-nocca-1)*nocca + l
                           if (open) ss(nopen+kl) = -ss(kl) - ss(nu+mt)
                           ss(nx+kl) = -ss(kl) + ss(nu+mt)
 80                     continue
 90                  continue
                  end if
                  if (nsoc.ne.0) then
                     do 110 k = ndpls1 , noccb
                        do 100 l = ndpls1 , k
                           kl = iky(mapie(k)) + mapie(l)
                           if (open) ss(nopen+kl) = -ss(kl)
                           ss(nx+kl) = -ss(kl)
 100                    continue
 110                 continue
                  end if
                  if (nocca.ne.0) then
                     do 130 k = ntpls1 , num
                        do 120 l = 1 , nocca
                           kl = iky(mapie(k)) + mapie(l)
                           mt = (k-nocca-1)*nocca + l
                           ss(nx+kl) = 2.0d0*ss(nu+mt)
 120                    continue
 130                 continue
                  end if
                  if (nsoc.ne.0) then
                     do 150 k = ntpls1 , num
                        do 140 l = ndpls1 , noccb
                           kl = iky(mapie(k)) + mapie(l)
                           mt = nocca*nvirta + (k-nsoc-1)*nsoc + l -
     +                          nocca
                           if (open) ss(nopen+kl) = ss(nu+mt)
                           ss(nx+kl) = ss(nu+mt)
 140                    continue
 150                 continue
                  end if
               end if
               if (odebug(3)) call prtris(ss(nx+1),ncoorb,iwr)

cc write perturbed density
cc               call wrt3(ss(nx+1),nw,iposd,ifockf)
               call stash(ss(nx+1),nw,'pdens',(i-1)*3+j)

cc               if (open) call wrt3(ss(nopen+1),nw,iposdo,ifockf)
cc              iposd = iposd + iblls
cc              if (open) iposdo = iposdo + iblls
 160        continue
 170     continue
cc         iochf(1) = max(iposd,iposdo)
         call clredx
         irest = 0
c
         return
      end if
 6010 format (//1x,'total perturbed density matrices in mo basis')
      end

      subroutine chfcla(qq,iqq)
c
c    assemble coupled hartree fock contribution to closed shell
c    scf second derivatives
c
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      logical out
c
INCLUDE(common/common)
INCLUDE(common/iofile)
INCLUDE(common/cndx41)
c
      common/maxlen/maxq
INCLUDE(common/infoa)
INCLUDE(common/nshel)
INCLUDE(common/mapper)
      common/bufb/ioffn(maxat*3),icol(maxorb)
INCLUDE(common/symtry)
      dimension qq(*),iqq(*)
INCLUDE(common/prnprn)
c
c     contibution from first derivative of density matrix
c     with first derivative of integrals
c
c     write(6,*)'start chfcla'
      nat3 = nat*3
      nlen = nat3*nat3
      iof = lenrel(nw196(5))
      iofs = nw196(5) + lenint(nat*nt)
      i10 = iofs + 1
      i11 = i10 + nlen
      maxw = maxq - nlen - nx - iofs
      if (maxw.lt.(nx*3)) call caserr(' insufficient core')
      do 20 i = 1 , nat3
         icol(i) = (i-1)*nat3 + iofs
 20   continue
      call vclr(qq(i10),1,nlen)
      ltri = ikyp(ncoorb)
cc      lenblk = lensec(ltri)
      out = odebug(6)
c     density matrix derivatives section 15
cc      iposd = iochf(15)
      npass = 1
      maxnuc = 0
      nadd = nat
      ntot = nx*nat3
 30   if (ntot.le.maxw) then
         do 50 ipass = 1 , npass
            minnuc = maxnuc + 1
            maxnuc = maxnuc + nadd
            if (maxnuc.gt.nat) maxnuc = nat
            nuc = maxnuc - minnuc + 1
            nuc3 = nuc*3
            k = (minnuc-1)*3 + 1
            ioffn(k) = nlen + nx + i10
            do 40 i = 1 , nuc3
               ioffn(k+1) = ioffn(k) + nx
               k = k + 1
 40         continue
 50      continue
         maxnuc = 0

         do 90 ipass = 1 , npass
c
c     derivatives of fock matrix (m.o. basis, no wavefunction term)
c     at section 13
c
cc            iposf = iochf(13)
            minnuc = maxnuc + 1
            maxnuc = maxnuc + nadd
            if (maxnuc.gt.nat) maxnuc = nat
            nuc = maxnuc - minnuc + 1
            nuc3 = nuc*3
            do 60 k = 1 , nuc3
               ioff = ioffn(k)
c
c     read perturbed density matrices
c
c#13
cc               call rdedx(qq(ioff),ltri,iposd,ifockf)
               call fetch(qq(ioff),ltri,'pdens',k)

cc               iposd = iposd + lenblk
 60         continue
            do 80 nn = 1 , nat3
c
c     read derivative fock matrix
c
cc               call rdedx(qq(i11),ltri,iposf,ifockf)
               call fetch(qq(i11),ltri,'dh',nn)

cc               iposf = iposf + lenblk
               k = (minnuc-1)*3
               do 70 kk = 1 , nuc3
                  ioff = ioffn(k+kk)
c
c     contribution to second derivative
c
                  dum = tracep(qq(i11),qq(ioff),ncoorb)
                  qq(k+kk+icol(nn)) = dum
 70            continue
 80         continue
 90      continue
         call rdedx(qq(1),nw196(5),ibl196(5),ifild)
         if (out) then
            call dr2sym(qq(i10),qq(i11),iqq(1),iqq(iof+1),nat,nat3,
     +      nshell)
            write (iwr,6010)
            call prnder(qq(i10),nat3,iwr)
         end if
c
c        write(6,*)'summing into derivative'
c
c#accum CHF term (NB performs dr2sym) - for a given pass
c
         call secget(isect(60),60,isec46)
         call rdedx(qq(i11),nlen,isec46,ifild)
         call vadd(qq(i11),1,qq(i10),1,qq(i11),1,nlen)
         call dr2sym(qq(i11),qq(i10),iqq(1),iqq(iof+1),nat,nat3,
     +               nshell)
         call wrt3(qq(i11),nlen,isec46,ifild)
         if (out) then
            write (iwr,6020)
            call prnder(qq(i11),nat3,iwr)
         end if
         return
      else
         npass = npass + 1
         nadd = nat/npass + 1
         ntot = nadd*3*nx
         go to 30
      end if
 6010 format (//' coupled hartree-fock contribution')
 6020 format (//' total so far')

      end

      subroutine symrhs(qq,skipp,mapnr,iso,nshels)
c     subroutine symrhs(qq,ibstar,skipp,mapnr,iso,nshels)
c
c    symmetrise r.h.s. ( d2h and subgroups only )
c
      implicit REAL  (a-h,o-z)
      logical skipp
INCLUDE(common/sizes)
      dimension skipp(3,nat),mapnr(*)
      dimension qq(*),iso(nshels,*)
c
INCLUDE(common/cigrad)
      common/mpshl/ns(maxorb)
INCLUDE(common/nshel)
INCLUDE(common/common)
INCLUDE(common/cndx40)
INCLUDE(common/cndx41)
INCLUDE(common/infoa)
INCLUDE(common/symtry)
INCLUDE(common/runopt)
c
      common/bufb/ptr(3,144),ict(maxat,8)
      common/symmos/imos(8,maxorb)
      data one/1.0d0/
c
c
      call rdedx(ptr(1,1),nw196(1),ibl196(1),ifild)
      nav = lenwrd()
      call readi(iso,nw196(5)*nav,ibl196(5),ifild)
      do 40 ii = 1 , nshell
         ic = katom(ii)
         do 30 it = 1 , nt
            id = iso(ii,it)
            ict(ic,it) = katom(id)
 30      continue
 40   continue
c
      do 80 n = 1 , nat
         do 50 i = 1 , 3
            skipp(i,n) = .false.
 50      continue
         do 70 nop = 1 , nt
            if (ict(n,nop).gt.n.or.  zopti(n).eq.'no') then
               do 60 i = 1 , 3
                  skipp(i,n) = .true.
 60            continue
               go to 80
            end if
 70      continue
 80   continue
c
      do 110 n = 1 , nat
         do 90 nop = 1 , nt
            if (ict(n,nop).ne.n) go to 110
 90      continue
         do 100 i = 1 , 3
            skipp(i,n) = .true.
 100     continue
         nuniq = n
         go to 120
 110  continue
 120  ntpls1 = noccb + 1
      nplus1 = nocca + 1
      nsoc = noccb - nocca
      nvirta = nsa4 - noccb
cc     ibll = lensec(mn)
cc     iblku = ibstar
      an = one/dfloat(nt)
      ioff = mn + 1
      nat3 = nat*3
c     read in u vectors
      do 130 n = 1 , nat3
         call fetch(qq(ioff),mn,'rhs',n)
cc         call rdedx(qq(ioff),mn,iblku,ifils)
cc        iblku = iblku + ibll
         ioff = ioff + mn
 130  continue
c     loop over vectors
cc     iblku = ibstar
      do 370 n = 1 , nat
         do 360 nc = 1 , 3
            ioff = ((n-1)*3+nc)*mn
c     copy vector for atom n, component nc into work area
            do 140 i = 1 , mn
               qq(i) = qq(ioff+i)
 140        continue
c     work along the elements of this vector
c loop over double-single and double-virtual
c N.B. closed-shell only
               if (nocca.ne.0) then
                  do 260 i = 1 , nocca
                     do 250 ia = nplus1 , nsa4
                        ij = (ia-nocca-1)*nocca + i
c     loop over symmetry operations
c     except identity
                        do 240 iop = 2 , nt
                           isign = imos(iop,i)*imos(iop,ia)
                           sign = dfloat(isign)
                           niop = ict(n,iop)
c     niop is the atom equivalent to n under operation
                           ioff = (niop-1)*3*mn
                           npnc = (iop-1)*3 + nc
                           do 230 k = 1 , 3
                              ioff = ioff + mn
                              qq(ij) = ptr(k,npnc)*sign*qq(ioff+ij)
     +                                 + qq(ij)
 230                       continue
 240                    continue
 250                 continue
 260              continue
               end if
c
            do 350 i = 1 , mn
               qq(i) = an*qq(i)
 350        continue

            call stash(qq(1),mn,'rhs',(n-1)*3+nc)
cc            call wrt3(qq(1),mn,iblku,ifils)
cc           iblku = iblku + ibll
c
 360     continue
 370  continue
      return
      end

      subroutine symmu(qq,skipp,iso,nshels)
c     subroutine symmu(qq,ibstar,skipp,iso,nshels)
c
c    symmetrise u-vectors ( chf solutions )
c    d2h and subgroups only
c
      implicit REAL  (a-h,o-z)
      logical skipp
INCLUDE(common/sizes)
      dimension skipp(3,nat)
      dimension qq(*),iso(nshels,*)
c
INCLUDE(common/nshel)
INCLUDE(common/common)
INCLUDE(common/cndx40)
INCLUDE(common/cndx41)
c
      common/mpshl/ns(maxorb)
INCLUDE(common/infoa)
INCLUDE(common/symtry)
      common/bufb/ptr(3,144),ict(maxat,8)
      common/symmos/imos(8,maxorb)
      character *8 grhf
      data grhf/'grhf'/
      data one/1.0d0/
c
      nav = lenwrd()
      call readi(iso,nw196(5)*nav,ibl196(5),ifild)
      call rdedx(ptr(1,1),nw196(1),ibl196(1),ifild)
c
      do 40 ii = 1 , nshell
         ic = katom(ii)
         do 30 it = 1 , nt
            id = iso(ii,it)
            ict(ic,it) = katom(id)
 30      continue
 40   continue
c
      nuniq = 0
      do 60 n = 1 , nat
         do 50 nop = 1 , nt
            if (ict(n,nop).ne.n) go to 60
 50      continue
         nuniq = n
         go to 70
 60   continue
 70   ntpls1 = noccb + 1
      nplus1 = nocca + 1
      nsoc = noccb - nocca
      nvirta = nsa4 - noccb
      ibll = lensec(mn)
cc     iblku = ibstar
      ioff = mn*3 + 1
      nat3 = nat*3
c     read in u vectors
      do 80 n = 1 , nat3
ccsecd
cc         call rdedx(qq(ioff),mn,iblku,ifils)
         call fetch(qq(ioff),mn,'soln',n)

cc        iblku = iblku + ibll
         ioff = ioff + mn
 80   continue
      mn3 = mn*3
c
c
c     loop over vectors
      do 270 n = 1 , nat
         if (.not.(skipp(1,n))) then
            ioff = n*mn3
c     copy vectors for atom n into work area
c
            do 90 i = 1 , mn3
               qq(i) = qq(ioff+i)
 90         continue
c
c     zero all elements related to components of atom n
c     by symmetry
c
            nsame = 0
            do 110 iop = 1 , nt
               niop = ict(n,iop)
               if (niop.eq.n) nsame = nsame + 1
               ioff = niop*mn3
               do 100 i = 1 , mn3
                  qq(ioff+i) = 0.0d0
 100           continue
 110        continue
            an = one/dfloat(nsame)
c
c     work along the elements of this vector
c loop over double-single and double-virtual
c
            if (scftyp.eq.grhf) then
               ij = 0
               do 160 i = 1 , nsa4
                  do 150 ia = 1 , i
                     if (ns(i).ne.ns(ia)) then
                        ij = ij + 1
c     loop over symmetry operations
                        do 140 iop = 1 , nt
                           isign = imos(iop,i)*imos(iop,ia)
                           sign = dfloat(isign)*an
                           niop = ict(n,iop)
c     niop is the atom equivalent to n under operation
                           do 130 nc = 1 , 3
                              ioff = (nc-1)*mn
                              npnc = (iop-1)*3 + nc
                              do 120 k = 1 , 3
                                 iof2 = (niop*3+k-1)*mn
                                 qq(iof2+ij) = qq(iof2+ij)
     +                              + sign*ptr(k,npnc)*qq(ioff+ij)
 120                          continue
 130                       continue
 140                    continue
                     end if
 150              continue
 160           continue
            else
               if (nocca.ne.0) then
                  do 210 i = 1 , nocca
                     do 200 ia = nplus1 , nsa4
                        ij = (ia-nocca-1)*nocca + i
c     loop over symmetry operations
                        do 190 iop = 1 , nt
                           isign = imos(iop,i)*imos(iop,ia)
                           sign = dfloat(isign)*an
                           niop = ict(n,iop)
c     niop is the atom equivalent to n under operation
                           do 180 nc = 1 , 3
                              ioff = (nc-1)*mn
                              npnc = (iop-1)*3 + nc
                              do 170 k = 1 , 3
                                 iof2 = (niop*3+k-1)*mn
                                 qq(iof2+ij) = qq(iof2+ij)
     +                              + sign*ptr(k,npnc)*qq(ioff+ij)
 170                          continue
 180                       continue
 190                    continue
 200                 continue
 210              continue
               end if
               if (noccb.ne.nocca) then
c open shell only - loop over single-virtual
                  do 260 i = nplus1 , noccb
                     do 250 ia = ntpls1 , nsa4
                        ij = nvirta*nocca + (ia-nsoc-1)*nsoc + i - nocca
c     loop over symmetry operations
                        do 240 iop = 1 , nt
                           isign = imos(iop,i)*imos(iop,ia)
                           sign = dfloat(isign)*an
                           niop = ict(n,iop)
c     niop is the atom equivalent to n under operation
                           do 230 nc = 1 , 3
                              ioff = (nc-1)*mn
                              npnc = (iop-1)*3 + nc
                              do 220 k = 1 , 3
                                 iof2 = (niop*3+k-1)*mn
                                 qq(iof2+ij) = qq(iof2+ij)
     +                              + sign*ptr(k,npnc)*qq(ioff+ij)
 220                          continue
 230                       continue
 240                    continue
 250                 continue
 260              continue
               end if
            end if
         end if
 270  continue
c
c
c     translational invariance
c
      iof2 = nuniq*mn3
      if (nuniq.ne.0) then
         do 280 i = 1 , mn3
            qq(iof2+i) = 0.0d0
 280     continue
         do 300 n = 1 , nat
            if (n.ne.nuniq) then
               ioff = n*mn3
               do 290 i = 1 , mn3
                  qq(iof2+i) = qq(iof2+i) - qq(ioff+i)
 290           continue
            end if
 300     continue
      end if
c
c     write it all out again
c
cc     iblku = ibstar
      ioff = mn3 + 1
      lenuu = ibll*nat3
      call secput(isect(66),66,lenuu,iblko)
      lds(isect(66)) = mn*nat3
      call revind
      do 310 n = 1 , nat3
cc         call wrt3(qq(ioff),mn,iblku,ifils)
         call stash(qq(ioff),mn,'soln',n)
c also copy to dumpfile
         call wrt3(qq(ioff),mn,iblko,ifild)

cc        iblku = iblku + ibll
         iblko = iblko + ibll
         ioff = ioff + mn
 310  continue
      return
      end

      subroutine anairr(q)
c
c  ------------------------------------------------
c    analysis of infra-red , raman and vcd effects
c    using previously calculated dipole derivatives
c    polarizability derivatives, vcd parameters and
c    cartesian force constants
c    ----------------------------------------------
c
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      common/maxlen/maxq
INCLUDE(common/infoa)
INCLUDE(common/iofile)
INCLUDE(common/prnprn)
      dimension q(*)
c
      nat3 = nat*3
      nnc = nat3*(nat3+1)/2
      itmp =  nnc + nnc + nat3 * (3*nat3) + 2 * nat3 + maxat*3
c
      i1 = igmem_alloc(itmp)
      i2 = i1 + nnc
      i3 = i2 + nnc
      i4 = i3 + nat3*nat3
      i5 = i4 + nat3
      i6 = i5 + nat3
      i7 = i6 + nat3*nat3
      i8 = i7 + nat3*nat3
      i9 = i8 + maxat*3
c
c  loop over isotopic substitution patterns
c
      nmv = mass_numvec()

      do imv=1,nmv
         
         if(nmv.ne.1)then
             write(iwr,100)imv
 100         format(/,1x,'Considering set of nuclear masses no.',i4,/)
         endif
         call rotcon(q(i1), q(nat+i1), q(2*nat+i1), q(3*nat+i1), nat,
     +        oprn(20),imv)
         oprn(20) = .false.
         call iranal(q(i1),q(i2),q(i3),q(i4),q(i5),q(i6),q(i7),q(i8),
     +        nat3,imv)
      enddo
c
c     ----- reset core allocation
c
      call gmem_free(i1)
c
      return
      end

      subroutine ver_secd_parallel(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/secd_parallel.m,v $
     +     "/
      data revision /"$Revision: 6176 $"/
      data date /"$Date: 2010-08-10 16:49:47 +0200 (Tue, 10 Aug 2010) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
