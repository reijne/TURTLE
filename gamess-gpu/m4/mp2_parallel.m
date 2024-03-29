c   
c  ***********  PARALLEL MP2 ENERGY AND GRADIENT MODULE  ******
c
c  this module uses square (num*num) orbital sets with the
c  superfluous (newbas0) vectors blanked out with *large* eigenvalues
c  jvl (2003)
c
      subroutine emp23(q,total_energy)
c
c  driving routine for parallel mp2 energy 
c
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
cfix
      logical orest
      dimension q(*),array(10)
c
INCLUDE(common/segm)
      common/maxlen/maxq
INCLUDE(common/timez)
INCLUDE(common/cigrad)
c
      logical lfield,fixed,lex,ldam12,ldam13,ldam23,ldiis
      common/scfblk/engy,etot,ehf,sh1(2),sh2(2),gap1(2),gap2(2),
     1              d12,d13,d23,canna,cannb,cannc,fx,fy,fz,
     2              lfield,fixed,lex,ldam12,ldam13,ldam23,ldiis,
     3              ncyc,ischm,lock,maxit,nconv,npunch,lokcyc
c
INCLUDE(common/vectrn)
INCLUDE(common/common)
INCLUDE(common/iofile)
INCLUDE(common/cndx41)
INCLUDE(common/infoa)
INCLUDE(common/mapper)
c
      common/junke/maxt,ires,ipass,nteff,
     1     npass1,npass2,lentrixx,nbuck,mloww,mhi,ntri,iacc,iontrn
c
INCLUDE(common/prnprn)
INCLUDE(common/restrj)
INCLUDE(common/cslosc)
INCLUDE(common/mp2grad_pointers)
INCLUDE(common/timeperiods)
INCLUDE(common/cntl2)
INCLUDE(common/global)
INCLUDE(common/nshel)
INCLUDE(common/sortp)
c
INCLUDE(common/mp2grd)
INCLUDE(common/cntl1)
INCLUDE(common/disp)
INCLUDE(common/errcodes)
      character*10 charwall

_IF(chemshell, charmm)
      external  mp2grad_pointers_init
      external  global_init
      external  mp2grd_init
_ENDIF

      logical opg_root
      character *8 hfsc
      data m1,m10,m13,m16/1,10,13,16/
      data hfsc/'hfscf'/

cjmht Should have been caught in the input but just to be safe...
      if (scftyp .eq. 'open') call gamerr(
     &     'UHF MP2 calculation not available in parallel',
     &     ERR_NO_CODE, ERR_INCOMPATIBLE, ERR_SYNC, ERR_NO_SYS)

      if ( runtyp .eq. 'gradient' .or.
     +     runtyp .eq. 'saddle'   .or. 
     +     runtyp .eq. 'optxyz'   .or.
     +     runtyp. eq. 'force'    .or.
     +     runtyp .eq. 'optimize' ) then
       if (isadle.ne.2) then
         opg_grad = .true.
       endif
      endif
c
      if ( runtyp .eq. 'force') go to 1000
c
c
c  evaluate integrals
c
      if (opass2) then
       if(nopk.eq.0.or.nopk.ne.nopkr.or.iofsym.eq.1.or.
     +    iofrst.ne.iofsym) then
        write (iwr,6020)
        opass2 = .false.
       endif
      endif
      nopk = 1
      iofsym = 0
      isecvv = isect(8)
      itypvv = 8
      nconv = max(nconv,7)
cfix
      orest = .false.
      if (mp2 .or. mp3) then
         call timit(3)
         t1 = tim
      end if
      if (.not.((mp2 .or. mp3) .and. mprest.ge.1)) then
c
c  calculate hf energy
c
         call integ(q)
*****
c     if(opg_alloc) call chk_ga_tools
*****
         call revise
         if (odisp) then
            call dispvec(q)
         else
            call scfrun(q)
         end if
*****
c     if(opg_alloc) call chk_ga_tools
*****

c  allow for incomplete scf i.e. maxcyc exceeded
         if (irest.ge.2) then
          orest = .true.
         else
          irest = 5
         endif

         call revise
      end if

      call timit(3)
      timscf = tim - t1
      mprest = max(1,mprest)
      call revise
cfix
      if (cpulft(0).lt.timscf.or.orest) then
         write (iwr,6010) mprest
         call parclenms('restart job')
      endif
c
 1000 continue
      emp2 = 0.0d0
      if(opg_root()) then 
       write(iwr,6040) cpulft(1) ,charwall()
      endif
c ...  schwarz inequality test
      if(.not.opg_alloc_sch) then
       opg_alloc_sch = .true.
       mp2grad_schw = igmem_alloc( nshell*(nshell+1)/2 )
      endif
c gdf:   switchable schwarz inequality 
      if (oschw) then
       if (nopk.ne.1) then
         i10 = igmem_alloc(151875)
       else
         i10 = igmem_alloc(50625)
       end if
       call coulmb(q(mp2grad_schw),q(i10))
       call gmem_free(i10)
      else
         call dcopy( nshell*(nshell+1)/2 
     &        ,0.0d0, 0, q(mp2grad_schw), 1)
         dlnmxs = 0.0d0
      end if
c
c ...  schwarz inequality test
c
      if (opg_grad) then
c  for optimisation or single-point mp2 gradient
        if (opg_grad_com ) then
c executed on first inovation of emp23, or if call to emp23 
c was preceeded by a call to grmp23
            call start_time_period(TP_MP2)
            if(.not.opg_alloc_dmat) then
             mp2grad_wmat = igmem_alloc(ncoorb*ncoorb)
             mp2grad_pmat = igmem_alloc(ncoorb*ncoorb)
             mp2grad_vecs = igmem_alloc(ncoorb*ncoorb)
             mp2grad_dens = igmem_alloc(ncoorb*ncoorb)
             mp2grad_vals = igmem_alloc(ncoorb)
             opg_alloc_dmat = .true.
            endif
            opg_alloc = .false.
         else
c executed if the last call was to emp23, and was not
c preceeded by a call to grmp23
            call end_time_period(TP_MP2)
            call pg_zero(g_oooo)
            call pg_zero(g_vooo)
            call pg_zero(g_vvoo)
            call pg_zero(g_vovo)
            opg_alloc = .true.
            call start_time_period(TP_MP2)
         end if
      else 
         if(opg_grad_com) then
c           write(6,*) 'opg_grad_com ',opg_grad_com
            opg_alloc = .false.
         else
c           write(6,*) opg_grad_com, g_vovo
            call pg_zero(g_vovo)
            opg_alloc = .true.
         endif
         call start_time_period(TP_MP2)
      end if
      call set41
      call aprdmp2(q, emp2,mp2)

      if (.not.opg_grad) then
c  delete GA's if single MP2 energy run only
         if(isadle.ne.2) then
          call delete_ga_inf(g_vovo,'vovo',fnm,snm)
         endif
         if(opg_alloc_sch) then
          call gmem_free(mp2grad_schw)
          opg_alloc_sch = .false.
         endif
         call end_time_period(TP_MP2)
      end if
      total_energy = etot + emp2
      array(1) = engy
      array(2) = ehf
      array(3) = total_energy
      array(4) = emp2
      do loop = 5,10
       array(loop) = 0.0d0
      enddo
      if (.not.odisp) write(iwr,6030) emp2, total_energy
      call secput(isect(494),m16,m1,iblk9)
      call wrt3(array,m10,iblk9,ifild)
c         mprest = 3
c
c ====  this restart needs attention ====
c
      irest = 0
      if (runtyp.eq.hfsc) mprest = 0
      if (runtyp .ne. 'force') then
c
c    read energy from dumpfile
c
         call secget(isect(13),m13,isec13)
         call rdedx(engy,lds(isect(13)),isec13,ifild)
         call revise
      endif
      opg_grad_com = .false.
      return
 6010 format (//'insufficient time , restart parameter =',i5)
 6020 format(/
     + 1x,'= = = = = = = = = = = = = = = = = = = = = = ='/
     + 1x,'= integrals must NOT be in supermatrix form ='/
     + 1x,'=        requested bypass is ignored        ='/
     + 1x,'= = = = = = = = = = = = = = = = = = = = = = ='/)
 6030 format(/10x,42('*')/10x,'mp2 calculation'/10x,42('*')
     +       /10x,'mp2 correlation energy     ',f15.8
     +       /10x,'total energy (mp2)         ',f15.8
     +       /10x,42('*'))
 6040 format(1x,'commence direct mp2 energy evaluation at', 
     +       f10.2,' seconds',a10,' wall')
      end
c
c  driving routine for parallel mp2 gradient
c
      subroutine grmp23(q,eg)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      dimension eg(*)
      dimension q(*)
INCLUDE(common/segm)
      common/maxlen/maxq
c
INCLUDE(common/common)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/cndx41)
c
INCLUDE(common/cigrad)
INCLUDE(common/infoa)
INCLUDE(common/mp2grad_pointers)
INCLUDE(common/timeperiods)
INCLUDE(common/nshel)
INCLUDE(common/cslosc)
INCLUDE(common/ijlab)
INCLUDE(common/machin)
INCLUDE(common/cntl2)
INCLUDE(common/mp2grd)
INCLUDE(common/statis)
INCLUDE(common/symtry)
INCLUDE(common/vcore)
c
      data m17 /17/
      mp2gcycl = jpoint
      call vclr(eg,1,nat*3)
      call mpgrad(q)

c     allocate gradient section on dumpfile
      ncoord = 3*nat
      nc2 = ncoord*ncoord
      isize = lensec(nc2) + lensec(mach(7))
      call secput(isect(495),m17,isize,ibl3g)
      ibl3hs = ibl3g + lensec(mach(7))

      iprefa = igmem_alloc( nshell*(nshell+1)/2 )
      iiso   = igmem_alloc( nw196(5) )
      call rdmake(q(iprefa))
c restore iso array
      call rdedx(q(iiso),nw196(5),ibl196(5),idaf)
      call cpuwal(begin,ebegin)
      call stvder(scftyp,q,q(iprefa),q(iiso),nshell)
      call timana(6)
      call gmem_free(iiso)
      call gmem_free(iprefa)
c
      istd = 1
      do i = 1 , nshell
        kad(i) = -1
      end do
      nindmx = 1
      call cpuwal(begin,ebegin)
      onocnt = .true.
      call jkder_ga(q,q(mp2grad_schw))
c
***   i00 = igmem_alloc((nshell*nt-1)/lenwrd()+1)
***   call rdedx(Q(i00),nw196(5),ibl196(5),idaf)
***   call jkder_ga(q,q(mp2grad_schw),Q(i00))
***   call gmem_free(i00)

      call gmem_free(mp2grad_vals)
      call gmem_free(mp2grad_dens)
      call gmem_free(mp2grad_vecs)
      call gmem_free(mp2grad_pmat)
      call gmem_free(mp2grad_wmat)
      call gmem_free(mp2grad_schw)
      opg_alloc_sch = .false.
      opg_alloc_dmat = .false.
      opg_grad_com = .true.
      if (mp2 .or. mp3) mprest = 0
      call end_time_period(TP_MP2)
c
      call timana(7)
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MO INTEGRAL CLASSES SECTION
c
c
c  initialisation routine for parallel transformation
c
      subroutine aprdmp2(q, emp2, omp2)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      REAL q
      dimension q(*)
INCLUDE(common/sizes)
INCLUDE(common/fsymas)
INCLUDE(common/modj)
INCLUDE(common/cslosc)
INCLUDE(common/symtry)
INCLUDE(common/nshel)
INCLUDE(common/atmol3)
INCLUDE(common/prints)
INCLUDE(common/scra7)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/mapper)
INCLUDE(common/timez)
INCLUDE(common/infoa)
INCLUDE(common/scfopt)
INCLUDE(common/restar)
INCLUDE(common/segm)
      common/diisd/st(210),cdiis(20),rdiis(19),derror,sdiis(20),
     + iposit(20),nsti(2),ondiis,junkj
INCLUDE(common/runlab)
      common/restri/jjfile(63),lds(508),isect(508),ldsect(508)
INCLUDE(common/psscrf)
INCLUDE(common/timeperiods)
INCLUDE(common/global)
INCLUDE(common/vcore)
INCLUDE(common/ijlab)
INCLUDE(common/mp2grad_pointers)
INCLUDE(common/tran)
INCLUDE(common/machin)
INCLUDE(common/harmon)
INCLUDE(common/direc)
      common/bufb/mmmm(65),isymao(maxorb),isymmo(maxorb)
      common/junkc/zjob(29)
_IF(parallel)
      common/nodeio/ids48(maxlfn),idr48(maxlfn),oputpp(maxlfn)
     + , maxio(maxlfn),oswed3(maxlfn)
_ENDIF
INCLUDE(common/mp2grd)
INCLUDE(common/disp)
      character*10 charwall
c
      data zrhf,zscf/'rhf','scf'/
      data dzero,two,pt5/0.0d0,2.0d0,0.5d0/
      data yblnk,yav/' ',' av'/
      data done,pt2,twopt2 /1.0d0,0.2d0,2.2d0/
      data dmptlc/1.0d-2/
      data igs/5/
      data m29,m51/29,51/
      out = nprint .eq. 5
      nav = lenwrd()
_IF(parallel)
      oswed3(4) = .true.
      oswed3(8) = .true.
_ENDIF
c gdf:
      call start_time_period(TP_APRDMP2)

      call start_time_period(TP_APRDMP2_I)
      jrun=locatc(yrunt,mxtask,ytrunc(zruntp))
      nocc=na
      nvir=num-nocc
      nvec=num
      if(nprint.ne.-5) then
       if(omp2) then
         write (iwr,9348)
       else
         write(iwr,9349) cpulft(1) ,charwall()
       endif
      endif
c
c allocate if not already
      if (.not.opg_alloc)  call alloc_ga2(na, num-na, num )
c
c     - the memory chunk needs to be num*num for this module
c
      l10 = num*num
      i10 = igmem_alloc(l10)
cjvl      l10 = num*newbas0  GAMESS writes square blocks
cjvl      next bit is the case for all
c      if (jrun.ge.4.and.jrun.le.8.or.jrun.eq.13.or.
c     +    jrun.eq.14. or.jrun.eq.18) then
c         l10=num*num
c         nvir=num-nocc
c         nvec=num
c      endif
      l20 = num
      i20 = igmem_alloc(l20)
      nav = lenwrd()
      l30 = nw196(5) 
      i30 = igmem_alloc(l30)
      mxshl=4
      do i=1,nshell
c        mxshl = max(mxshl,kmax(i)-kmin(i)+1)
         if(ktype(i).eq.3) then
          mxshl = max(mxshl,6)
         else if(ktype(i).eq.4) then
          mxshl = max(mxshl,10)
         else if(ktype(i).eq.5) then
          mxshl = max(mxshl,15)
         else
         endif
      enddo
      l50 = mxshl**4 
      i50 = igmem_alloc(l50)
      l60 = mxshl*mxshl*num*nocc 
      i60 = igmem_alloc(l60)
      l80 = max(num*num,num*mxshl*mxshl)
      i80 = igmem_alloc(l80)
      l90 = max(num*num,num*mxshl*mxshl)
      i90 = igmem_alloc(l90)
c
c generate mo symmetries
c
      call secget(isect(490),m51,iblk51)
      call readi(mmmm,mach(13)*nav,iblk51,idaf)
      call secget(mouta,3,iblk51)
      call rdchr(zjob,m29,iblk51,idaf)
      call symvec(Q(i10),isymao,isymmo,num,num,iblk51)
c
c     restore vectors and eigenvalues
c
      call get_mo_coeffs(Q(i10),l10)
      call rdedx(Q(i20),l20,ibl3ea,idaf)
      if (newbas0.ne.newbas1) then
c
c     clear non-existent vectors and make their eigenvalues biggg
c
         call vclr(q(i10+newbas0*num),1,(newbas1-newbas0)*num)
         do i=newbas0+1,newbas1
            q(i20+i-1) = 1.0d99
         end do
      end if
c
c flag shells to compute - should hardwire independent of input..
      call spchck
      call debut(zrhf)
c restore iso array
      call rdedx(Q(i30),nw196(5),ibl196(5),idaf)
      mxshlt=mxshl*mxshl

      call end_time_period(TP_APRDMP2_I)
c
      call aprm1234( q(mp2grad_schw) ,Q(i30),
     &     Q(i50),nshell,
     &     Q(i10),Q(i60),Q(i80),
     &     Q(i90),
     &     na,nvir,num,nvec,mxshl)
c
c
c  debug check of global arrays
c      call chk_ga2(na,num-na,num)
      if (.not.odisp) then
         call apremp2(na,nvir,Q(i20),Q(i80),
     &        Q(i90),emp2)
      else
         call aprdisp(na,nvir,Q(i20),Q(i80),
     &        Q(i90),emp2)
      end if
c
      call gmem_free(i90)
      call gmem_free(i80)
      call gmem_free(i60)
      call gmem_free(i50)
      call gmem_free(i30)
      call gmem_free(i20)
      call gmem_free(i10)
      call end_time_period(TP_APRDMP2)
****
c     call chk_ga_tools
c     write (6,*) 'end of aprdmp2'
c     call wrtmat_ga('(VO|OO)  ',g_vooo)
****
      return
c
 9348 format(//40x,'**************************'/
     +         40x,'* Direct MP2 Calculation *'/
     +         40x,'**************************'/)
 9349 format(/1x,
     + 'commence GA-based integral transformation at ',
     +  f9.2,' seconds',a10,' wall')
 9148 format(//1x,104('-')//
     *  50x,12('-')/50x,'eigenvectors'/50x,12('-'))
 9308 format(' core assignement '/(10i8))
      end
c
c   code from mptran to set up nocca etc
c
      subroutine set41
      implicit none
INCLUDE(common/sizes)
INCLUDE(common/cndx41)
INCLUDE(common/cndx40)
INCLUDE(common/mapper)
INCLUDE(common/prnprn)
INCLUDE(common/iofile)
INCLUDE(common/infoa)
INCLUDE(common/common)
c
      logical uhf
      character *8 open
      data open/'open'/

      integer i, j
c
      uhf = scftyp.eq.open

      nocc = nb
      nvirt = ncoorb - nb
      nocca = 0
      noccb = 0
      if (.not.uhf) then
         do 20 i = 1 , nsa4
            j = mapie(i)
            if (j.le.na) nocca = nocca + 1
            if (j.gt.na .and. j.le.nb) noccb = noccb + 1
 20      continue
         noccb = noccb + nocca
         nvirta = nsa4 - noccb
         if (oprn(6) .and. nprint.ne.-5) then
            write (iwr,6050) num , ncoorb , nsa4 , nocc , noccb
         end if
      else
         nsb = nsa4
         do 30 i = 1 , nsb
            j = mapie(i)
            if (j.le.na) nocca = nocca + 1
            if (j.le.nb) noccb = noccb + 1
 30      continue
         nvirta = nsa4 - nocca
         nvirtb = nsb - noccb
         if (oprn(6) .and. nprint.ne.-5) then
            write (iwr,6060) num , ncoorb , nsa4 , nsb , nocca , 
     +           noccb
         end if
      end if
c
 6050 format (/1x,'number of basis functions          ',i4/1x,
     +        'number of molecular orbitals       ',i4/1x,
     +        'number of active molecular orbitals',i4/1x,
     +        'number of occupied orbitals        ',i4/1x,
     +        'number of occupied active orbitals ',i4)
 6060 format (/1x,'number of basis functions                 ',i4/1x,
     +        'number of molecular orbitals              ',i4/1x,
     +        'number of active alpha molecular orbitals ',i4/1x,
     +        'number of active beta molecular orbitals  ',i4/1x,
     +        'number of occupied alpha orbitals         ',i4/1x,
     +        'number of occupied beta orbitals          ',i4)
      return
      end
c
c  parallel transformation driver
c
      subroutine aprm1234(schwa,iso,gout,nshels,cmo,trn1,
     &                    tmp1,tmp2,nocc,nvir,nbas,nvec,mxshl)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      logical sym_pair
      external sym_pair
      logical sym_quartet
      external sym_quartet
      logical osym, osym_pair
      dimension iso(nshels,*),gout(*),schwa(*)
      REAL cmo(nbas,nvec),trn1(nbas,mxshl*mxshl,nocc),
     &     tmp1(nbas*nbas),tmp2(nbas*nbas)
INCLUDE(common/sizes)
      parameter (mxp2 = mxprms * mxprms)
INCLUDE(common/cslosc)
INCLUDE(common/mapper)
INCLUDE(common/restar)
INCLUDE(common/iofile)
INCLUDE(common/utilc)
INCLUDE(common/nshel)
c
INCLUDE(common/shlnos)
INCLUDE(common/shlt)
INCLUDE(common/symtry)
INCLUDE(common/ijlab)
INCLUDE(common/shlg70)
INCLUDE(common/picon)
INCLUDE(common/misc)
c
INCLUDE(common/timez)
INCLUDE(common/pkfil)
INCLUDE(common/flips)
INCLUDE(common/sortp)
INCLUDE(common/parallel)
c
INCLUDE(common/global)
INCLUDE(common/timeperiods)
c
c     ----- size of gout -
c                         1   if s or k shells
c                        81   if p      shells
c                       256   if      l shells
c                      1296   if d or m shells
c                     10000   if f shells
c                     50625   if g shells
c
c     ----- this version can handle g shells   -----
c
      dimension ib(4,4)
      data ib/64,16,4,1,216,36,6,1,1000,100,10,1,
     +         3375,225,15,1 /
      data done/1.0d0/,two/2.0d0/,twopt5/2.5d0/,four/4.0d0/
c
      call start_time_period(TP_APRM1234)

c  debug write mo coeffs
c      print *,'mo coeffs'
c      do i99 = 1 , nvec
c        write(*,912)(cmo(j99,i99),j99=1,nbas)
c912     format(10f12.8)
c      end do

c
c     ----- initialize parameters -----
c
      dlncutoff = dlog(cutoff)
      q4 = 1.0d0
      osym = .true.
      osym_pair = .true.
      qq4 = 1.0d0
      ii = 1
      if (opdbas) then
         ii = 2
      end if
      if (opfbas) then
         ii = 3
      end if
      if (opgbas) then
         ii = 4
      end if
      do 30 loop = 1 , 3
        igt(loop) = ib(1,ii)
        jgt(loop) = ib(2,ii)
        kgt(loop) = ib(3,ii)
        lgt(loop) = ib(4,ii)
 30   continue
      pidiv4 = datan(done)
      pi = four*pidiv4
      pito52 = two*pi**twopt5
      if (intg76.ne.0) call sinset
      nschwz = 0
      time = cpulft(1)
      tim0 = time
      tim1 = time
      ist0 = 1
      jst0 = 1
      kst0 = 1
      lst0 = 1
ccccc             are the rotated axis integrals to be used only?
      ogauss = intg76.eq.0. or. opdbas .or. opfbas .or. opgbas
      noffij=0
      nproc= ipg_nnodes()
      me=ipg_nodeid()
cgdf  18.03.08  .true. -> nprint.ne.-5
      call pg_dlbchunk(1, nprint.ne.-5 )
      call pg_dlbreset()
      mytask=ipg_dlbtask()
      mxshl2 = mxshl*mxshl
      ltrn1 = nocc*nbas*mxshl2
      if (intg76.ne.0) call filmax
c
c     ----- ishell -----
c
      do 150 ii = ist0 , nshels
       ishell = ii
       mini = kmin(ishell)
       maxi = kmax(ishell)
       loci = kloc(ishell)-mini
       nishl = maxi-mini+1
       kadi = kad(ii)
       dt0 = time - tim0
       dt1 = time - tim1
       tim1 = time
_IF(newints)
       oisp = nishl.eq.4
_ENDIF
c
c     ----- jshell -----
c
       do 140 jj = jst0 , ii
        jst0 = 1
        jshell = jj
        minj = kmin(jj)
        maxj = kmax(jj)
        locj = kloc(jshell)-minj
        njshl = maxj - minj + 1
        nijshl = nishl*njshl
        oianj = ishell .eq. jshell
        if (oianj) nijshl = nishl*(nishl+1)/2
        kadij = kad(jj) + kadi
        itrij = iky(max(ii,jj)) + min(jj,ii)
c
c  symmetry (IJ| 
c
        if (nt.gt.1) then
         osym_pair = sym_pair(iso,nshels,nt,ii,jj)
        endif
        if (osym_pair) then
        if (schwa(itrij)+dlnmxs.lt.dlncutoff) then
           nschw = nschw + nshels*(nshels+1)/2
           noffij = noffij + nijshl
           go to 140
        endif
c
c  parallel load balancer
c
        icount_dlb = icount_dlb + 1
        if (icount_dlb.eq.mytask)then
          mytask=ipg_dlbtask()
c
c     ----- get information about i-shell and j-shell -----
c     ----- for gauss-rys routine
c
         if (ogauss) then
           call shells(gout,1,ishell,jshell,ishell,jshell,1)
           call ijprim
           if (nij.eq.0) go to 170
         endif
         call dcopy(ltrn1,0.0d00,0,trn1,1)
_IF(newints)
         oijsp = oisp.or.(njshl.eq.4)
_ENDIF
c
c     ----- kshell -----
c
         do 130 kk = kst0 , nshels
          kst0 = 1
          kshell = kk
          kadijk = kad(kk) + kadij
          mink = kmin(kshell)
          maxk = kmax(kshell)
          lock = kloc(kshell)-mink
          nkshl = maxk - mink + 1
_IF(newints)
          oijksp = oijsp.or.(nkshl.eq.4)
_ENDIF
c
c     ----- lshell ----
c
          do 120 ll = lst0 , kk
           lst0 = 1
           lshell = ll
           minl = kmin(lshell)
           maxl = kmax(lshell)
           locl = kloc(lshell)-minl
           nlshl = maxl - minl + 1
           nklshl = nkshl*nlshl
           okanl = kshell .eq. lshell
           if (okanl) nklshl= nkshl*(nkshl+1)/2
           itrkl = iky(max(ll,kk)) + min(ll,kk)
c
c  schwarz inequality test
c
           test = schwa(itrij) + schwa(itrkl)
           if (test.lt.dlncutoff) then
             nschwz = nschwz + 1
           else
c
c  symmetry (IJ|KL)
c
             if (nt.gt.1) then
              osym = sym_quartet(iso,nshels,nt, ii, jj, kk, ll, q4)
             endif
             if (osym) then
c
c             print *,'unique quartet'
c             print *,'   ',ii,jj,kk,ll,'  q4 =',nint(q4)
c             print *,'shell  ',ii,jj,kk,ll
              qq4 = q4
              if (kadijk+kad(lshell).lt.0)then
               call aprshel(gout,2,ishell,jshell,kshell,lshell,1)
c              call start_time_period(TP_GENRAL)
c
c  compute (IJ|KL)
c
               call genral(gout)
c              call end_time_period(TP_GENRAL)
c
c  first transformation
c
               call aprq1(gout,cmo,trn1,nbas,nvec,nocc,mxshl2) 
              else
_IF(newints)
               oijklsp = oijksp.or.(nlshl.eq.4)
               if (oijklsp) then
                  isptype = 0
               else
                  isptype = 6
               endif
_ENDIF
c              call start_time_period(TP_GENRAL)
               call genr70(gout,1,.false.)
c              call end_time_period(TP_GENRAL)
               call aprq1s(gout,cmo,trn1,nbas,nvec,nocc,mxshl2) 
              endif
c             call aprint(gout)
             end if
          end if
 120      continue
 130     continue
_IF()
c  debug print first transformation integrals 
c        print *,'first transformation (',ii,jj,'| lamda I)'
c        do i99 = 1 , nijshl
c          do j99 = 1 , nocc
c            write(*,913)j99,(trn1(k99,i99,j99),k99=1,nbas)
c913          format(1i2,3x,10f12.8)
c          end do
c        end do
_ENDIF
c
c  second (& third) transformations
c
         if (g_vvvo.ge.0)then
          call aprq2(cmo,trn1,tmp1
     &,   nbas,nvec,nocc,nvir,mxshl*mxshl,nijshl,noffij)
         else
          call aprq2d(cmo,trn1,tmp1,tmp2
     &,   nbas,nvec,nocc,nvir,mxshl,nijshl,noffij)
         endif
 170     continue
         end if
         end if
         noffij=noffij+nijshl
 140   continue
       time = cpulft(1)
 150  continue
      call pg_dlbreset()

c  debug print third transformation for vovo integrals
c     print *,'third transformation for vovo integrals'
c     call wrtmat_ga('(VO|VO)  ',g_vovo)

c
c  (third &) fourth transformations
c
      call pg_synch(101)

      if (g_vvvo.ge.0)then
       call aprq34(cmo,tmp1,tmp2,nbas,nvec,nocc,nvir,nshels)
      else
       call aprq34d(cmo,tmp1,tmp2,nbas,nvec,nocc,nvir,nshels)
      endif
c 
c  remove mo integrals zero-by-symmetry
c 
****
*     write (6,*) 'prior to sym_mo_ints'
*     call wrtmat_ga('(VO|OO)  ',g_vooo)
****
      if (nt.gt. 1) then
         call pg_synch(99)
         if (g_vvvo.ge.0)then
          call sym_mo_ints_d(tmp1,nbas,nocc,nvir,nprint,iwr)
         else
          call sym_mo_ints(tmp1,nbas,nocc,nvir,nprint,iwr)
         endif
      endif
      call pg_synch(99)
      if(nprint.ne.-5) write(iwr,6040) nschwz
_IF(newints)
      isptype = -1
_ENDIF
      call end_time_period(TP_APRM1234)

      return
6040  format(/1x,'schwarz inequality test skipped ',i10,
     + ' integral blocks')
      end
c
c  determine highest-index unique pairs
c
      logical function sym_pair(iso,nshels,nops,iat,jat)
      implicit none
      integer iat, jat          ! shell indices [input]
      integer nshels
      integer iso(nshels,*)
      integer iat_new, jat_new
      integer nops, op, ij, ij_new
INCLUDE(common/sizes)
INCLUDE(common/mapper)

      sym_pair = .false.
      if (iat .lt. jat) return  ! Labels must be in canonical order
      ij   = iky(max(iat,jat))+min(iat,jat)
c
c     Loop thru operations in the group and map to new pairs
c
      do op = 2, nops
c
c     Map centers
c
         iat_new = iso( iat, op)
         jat_new = iso( jat, op)
c
c     Compare canonical pair indices
c
         ij_new   = iky(max(iat_new,jat_new))
     +                + min(iat_new,jat_new)
         if (ij .lt. ij_new) return
      end do
      sym_pair = .true.
      end
c
c  determine unique quartets
c
      logical function sym_quartet(iso,nshels,nops
     $,     iat, jat, kat, lat, q4)
      implicit none
      integer iat, jat, kat, lat ! shell indices [input]
      REAL q4       ! Constituency number [output]
      integer nshels
      integer iso(nshels,*)
      integer iat_new, jat_new, kat_new, lat_new
      integer nops, op, ij, kl, ij_new, kl_new
      integer n                 ! Counts no. of equivalent pairs
INCLUDE(common/sizes)
INCLUDE(common/mapper)
c
      q4 = 0.0d0
      sym_quartet = .false.
c
      if (iat .lt. jat) return  ! Labels must be in canonical order
      if (kat .lt. lat) return
c
      ij   = iky(max(iat,jat))+min(iat,jat)
      kl   = iky(max(kat,lat))+min(kat,lat)
c
c     Loop thru operations in the group and map to new pairs
c
      n = 1                     ! Identity always counts
      do op = 2, nops
c
c     Map centers
c
         iat_new = iso( iat, op)
         jat_new = iso( jat, op)
         kat_new = iso( kat, op)
         lat_new = iso( lat, op)
c
c     Compare canonical pair indices
c
         ij_new   = iky(max(iat_new,jat_new))+min(iat_new,jat_new)
         kl_new   = iky(max(kat_new,lat_new))+min(kat_new,lat_new)
c
         if (ij .lt. ij_new) return
         if (ij.eq.ij_new .and. kl.lt.kl_new) return
         if (ij.eq.ij_new .and. kl.eq.kl_new) then
            n = n + 1
         endif
      end do
c
      q4 = dfloat(nops) / dfloat(n)
      if (abs(q4-nint(q4)).gt.1d-12) call caserr
     $     ('sym_quartet: not divisible')
      sym_quartet = .true.
c
      end
c  
c  remove mo integrals zero-by-symmetry
c
      subroutine sym_mo_ints(buf,nbas,nocc,nvir,nprint,iwr)
      implicit REAL  (a-h,o-z)
      integer i,j,k,l,ij,kl,ijk,is,js,ks,ls, nbas,nocc,nvir
      REAL buf(*)
INCLUDE(common/sizes)
INCLUDE(common/global)
INCLUDE(common/mp2grd)
      common/bufb/mmmm(65),isymao(maxorb),isymmo(maxorb)
      data zero/0.0d0/

      if(nprint.ne.-5)write(iwr,70)
70    format(/38x,'mo symmetries'/38x,13('*')/)
      if(nprint.ne.-5)write(iwr,80)(isymmo(i),i=1,nbas)
80    format(10x,14i5)

c  vovo class
      ij = 0
      len = nvir*nvir
      do i = 1 , nocc
        do j = 1 , i
          ij = ij + 1
          if (ij.ge.jl_vovo .and. ij.le.jh_vovo) then
            call pg_get(g_vovo,1,len,ij,ij,buf,1)
            kl = 0
            do k = nocc+1 , nbas
              do l = nocc+1 , nbas
                kl = kl + 1
                is = isymmo(i)-1
                js = isymmo(j)-1
                ks = isymmo(k)-1
                ls = isymmo(l)-1
                if ( ieor (ieor ( ieor(is,js), ks), ls ).eq.0 ) then
c  integral is non-vanishing by symmetry
                else
                  buf(kl) = zero
                end if
              end do
            end do
            call pg_put(g_vovo,1,len,ij,ij,buf,1)
          end if
        end do
      end do
      if (opg_grad) then 
c  vvoo class
        ij = 0
        len = nvir*(nvir+1)/2
        do i = 1 , nocc
          do j = 1 , i
            ij = ij + 1
            if (ij.ge.jl_vvoo .and. ij.le.jh_vvoo) then
              call pg_get(g_vvoo,1,len,ij,ij,buf,1)
              kl = 0
              do k = nocc+1 , nbas
                do l = nocc+1 , k
                  kl = kl + 1
                  is = isymmo(i)-1
                  js = isymmo(j)-1
                  ks = isymmo(k)-1
                  ls = isymmo(l)-1
                  if ( ieor (ieor ( ieor(is,js), ks), ls ).eq.0 ) then
c  integral is non-vanishing by symmetry
                  else
                    buf(kl) = zero
                  end if
                end do
              end do
              call pg_put(g_vvoo,1,len,ij,ij,buf,1)
            end if
          end do
        end do
c  vooo class
        ijk = 0
        len = nvir
        do i = 1 , nocc
          do j = 1 , i
            do k = 1 , nocc
              ijk = ijk + 1
              if (ijk.ge.jl_vooo .and. ijk.le.jh_vooo) then
                call pg_get(g_vooo,1,len,ijk,ijk,buf,1)
                kl = 0
                do l = nocc+1 , nbas
                  kl = kl + 1
                  is = isymmo(i)-1
                  js = isymmo(j)-1
                  ks = isymmo(k)-1
                  ls = isymmo(l)-1
                  if ( ieor (ieor ( ieor(is,js), ks), ls ).eq.0 ) then
c  integral is non-vanishing by symmetry
                  else
                    buf(kl) = zero
                  end if
                end do
                call pg_put(g_vooo,1,len,ijk,ijk,buf,1)
              end if
            end do
          end do
        end do
c  oooo class
        ij = 0
        len = nocc*(nocc+1)/2
        do i = 1 , nocc
          do j = 1 , i
            ij = ij + 1
            if (ij.ge.jl_oooo .and. ij.le.jh_oooo) then
              call pg_get(g_oooo,1,len,ij,ij,buf,1)
              kl = 0
              do k = 1 , nocc
                do l = 1 , k
                  kl = kl + 1
                  is = isymmo(i)-1
                  js = isymmo(j)-1
                  ks = isymmo(k)-1
                  ls = isymmo(l)-1
                  if ( ieor (ieor ( ieor(is,js), ks), ls ).eq.0 ) then
c  integral is non-vanishing by symmetry
                  else
                    buf(kl) = zero
                  end if
                end do
              end do
              call pg_put(g_oooo,1,len,ij,ij,buf,1)
            end if
          end do
        end do
      end if
      return
      end
c  
c  remove mo integrals zero-by-symmetry
c  in-core vvvo version
c
      subroutine sym_mo_ints_d(buf,nbas,nocc,nvir,nprint,iwr)
      implicit REAL  (a-h,o-z)
      integer i,j,k,l,ij,kl,ijk,is,js,ks,ls, nbas,nocc,nvir
      REAL buf(*)
INCLUDE(common/sizes)
INCLUDE(common/global)
INCLUDE(common/mp2grd)
      common/bufb/mmmm(65),isymao(maxorb),isymmo(maxorb)
      data zero/0.0d0/

c  vovo class
      ij = 0
      len = nvir*nvir
      do i = 1 , nocc
        do j = 1 , i
          ij = ij + 1
          if (ij.ge.jl_vovo .and. ij.le.jh_vovo) then
            call pg_get(g_vovo,1,len,ij,ij,buf,1)
            kl = 0
            do k = nocc+1 , nbas
              do l = nocc+1 , nbas
                kl = kl + 1
                is = isymmo(i)-1
                js = isymmo(j)-1
                ks = isymmo(k)-1
                ls = isymmo(l)-1
                  if ( ieor (ieor ( ieor(is,js), ks), ls ).eq.0 ) then
c  integral is non-vanishing by symmetry
                else
                  buf(kl) = zero
                end if
              end do
            end do
            call pg_put(g_vovo,1,len,ij,ij,buf,1)
          end if
        end do
      end do
      if (opg_grad) then 
c  vvvo class
        ij = 0
        len = nvir*(nvir+1)/2
        do i = 1 , nocc
          do j = nocc+1 , nbas
            ij = ij + 1
            if (ij.ge.jl_vvvo .and. ij.le.jh_vvvo) then
              call pg_get(g_vvvo,1,len,ij,ij,buf,1)
              kl = 0
              do k = nocc+1 , nbas
                do l = nocc+1 , k
                  kl = kl + 1
                  is = isymmo(i)-1
                  js = isymmo(j)-1
                  ks = isymmo(k)-1
                  ls = isymmo(l)-1
                  if ( ieor (ieor ( ieor(is,js), ks), ls ).eq.0 ) then
c  integral is non-vanishing by symmetry
                  else
                    buf(kl) = zero
                  end if
                end do
              end do
              call pg_put(g_vvvo,1,len,ij,ij,buf,1)
            end if
          end do
        end do
c  vvoo class
        ij = 0
        len = nvir*(nvir+1)/2
        do i = 1 , nocc
          do j = 1 , i
            ij = ij + 1
            if (ij.ge.jl_vvoo .and. ij.le.jh_vvoo) then
              call pg_get(g_vvoo,1,len,ij,ij,buf,1)
              kl = 0
              do k = nocc+1 , nbas
                do l = nocc+1 , k
                  kl = kl + 1
                  is = isymmo(i)-1
                  js = isymmo(j)-1
                  ks = isymmo(k)-1
                  ls = isymmo(l)-1
                  if ( ieor (ieor ( ieor(is,js), ks), ls ).eq.0 ) then
c  integral is non-vanishing by symmetry
                  else
                    buf(kl) = zero
                  end if
                end do
              end do
              call pg_put(g_vvoo,1,len,ij,ij,buf,1)
            end if
          end do
        end do
c  vooo class
        ijk = 0
        len = nvir
        do i = 1 , nocc
          do j = 1 , i
            do k = 1 , nocc
              ijk = ijk + 1
              if (ijk.ge.jl_vooo .and. ijk.le.jh_vooo) then
                call pg_get(g_vooo,1,len,ijk,ijk,buf,1)
                kl = 0
                do l = nocc+1 , nbas
                  kl = kl + 1
                  is = isymmo(i)-1
                  js = isymmo(j)-1
                  ks = isymmo(k)-1
                  ls = isymmo(l)-1
                  if ( ieor (ieor ( ieor(is,js), ks), ls ).eq.0 ) then
c  integral is non-vanishing by symmetry
                  else
                    buf(kl) = zero
                  end if
                end do
                call pg_put(g_vooo,1,len,ijk,ijk,buf,1)
              end if
            end do
          end do
        end do
c  oooo class
        ij = 0
        len = nocc*(nocc+1)/2
        do i = 1 , nocc
          do j = 1 , i
            ij = ij + 1
            if (ij.ge.jl_oooo .and. ij.le.jh_oooo) then
              call pg_get(g_oooo,1,len,ij,ij,buf,1)
              kl = 0
              do k = 1 , nocc
                do l = 1 , k
                  kl = kl + 1
                  is = isymmo(i)-1
                  js = isymmo(j)-1
                  ks = isymmo(k)-1
                  ls = isymmo(l)-1
                  if ( ieor (ieor ( ieor(is,js), ks), ls ).eq.0 ) then
c  integral is non-vanishing by symmetry
                  else
                    buf(kl) = zero
                  end if
                end do
              end do
              call pg_put(g_oooo,1,len,ij,ij,buf,1)
            end if
          end do
        end do
      end if
      return
      end
c
c  first transformation
c
      subroutine aprq1(g,cmo,trn1,nbas,nvec,nocc,mxshlt)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension g(*),cmo(nbas,nvec),trn1(nbas,mxshlt,nocc)
INCLUDE(common/misc)
INCLUDE(common/shlnos)
INCLUDE(common/indez)
INCLUDE(common/timeperiods)
c
c       call start_time_period(TP_APRQ1)

c  debug print integrals
c       print *,'shell integrals (weighted by q4)'
c       print *,'shell integrals'

       ijn = 0
       jmax = maxj
       do i = mini,maxi
        if (oianj) jmax = i
        do j = minj,jmax
         ijn = ijn+1
         kln = 0
         do k = mink,maxk
          lmax = maxl
          if (okanl) lmax = k
          i3 = lock + k
          do l = minl,lmax
           kln = kln+1
           nn = ijgt(ijn)+klgt(kln)
           i4 = locl + l
           val = g(nn)

c           write(*,912) loci+i,locj+j,lock+k,locl+l,val
c912        format(4x,4i2,1f12.8)

           if (i3.eq.i4)val=val*0.5d00
           do iocc=1,nocc
            trn1(i3,ijn,iocc)=trn1(i3,ijn,iocc)+val*cmo(i4,iocc)
            trn1(i4,ijn,iocc)=trn1(i4,ijn,iocc)+val*cmo(i3,iocc)
           enddo
          enddo
         enddo
        enddo
       enddo
c      call end_time_period(TP_APRQ1)
      return
      end
c
c  first transformation fast integrals version
c
      subroutine aprq1s(g,cmo,trn1,nbas,nvec,nocc,mxshlt)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension g(*),cmo(nbas,nvec),trn1(nbas,mxshlt,nocc)
INCLUDE(common/misc)
INCLUDE(common/shlnos)
INCLUDE(common/indez)
INCLUDE(common/timeperiods)
c
      common/flip70/ib1,jb1,kb1,lb1,ib2,jb2,kb2,lb2,ib3,jb3,kb3,lb3
      dimension ib(4,4)
      data ib/0,0,0,0,64,16,4,1,128,32,8,2,192,48,12,3/
c      call start_time_period(TP_APRQ1)
       ijn = 0
       jmax = maxj
       do i = mini,maxi
        if (oianj) jmax = i
        do j = minj,jmax
         ijn = ijn+1
         n1 = ib(ib1,i)+ib(jb1,j)+1
         do k = mink,maxk
          lmax = maxl
          if (okanl) lmax = k
          i3 = lock + k
          do l = minl,lmax
           nn = n1+ib(kb1,k)+ib(lb1,l)
           i4 = locl + l
           val = g(nn)
           if (i3.eq.i4)val=val*0.5d0
           do iocc=1,nocc
            trn1(i3,ijn,iocc)=trn1(i3,ijn,iocc)+val*cmo(i4,iocc)
            trn1(i4,ijn,iocc)=trn1(i4,ijn,iocc)+val*cmo(i3,iocc)
           enddo
          enddo
         enddo
        enddo
       enddo
c      call end_time_period(TP_APRQ1)
      return
      end
c
c  second transformation
c
      subroutine aprq2d(cmo,trn1,tmp1,tmp2
     &,                 nbas,nvec,nocc,nvir,mxshl,nijshl,noffij)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension cmo(nbas,nbas),trn1(nbas*mxshl*mxshl,nocc),
     &          tmp1(nijshl,nbas),tmp2(nbas,nbas)
INCLUDE(common/global)
INCLUDE(common/timeperiods)
INCLUDE(common/mp2grd)
      character*1 xn, xt
      data xn,xt/'n','t'/
c
      call start_time_period(TP_APRQ2D)
      do iocc=1,nocc
c  second transformation
        call dgemm(xt,xn,nijshl,nvec,nbas
     &,            1.0d00,trn1(1,iocc),nbas
     &,            cmo,nbas,0.0d00,tmp1,nijshl)

c  debug print second transformation integrals
c       print *,'second transformation, mo: ',iocc
c       do i99 = 1 , nijshl
c         write(*,912)(tmp1(i99,j99),j99=1,nbas)
c912       format(10f12.8)
c       end do

        if (opg_grad) then
c         call start_time_period(TP_GA_PUT_Q2D)
          do jocc=1,iocc
            kl=iocc*(iocc-1)/2+jocc
            call pg_put(g_vvoo,noffij+1,noffij+nijshl,kl,kl,
     &                  tmp1(1,jocc),nijshl)
          enddo
c         call end_time_period(TP_GA_PUT_Q2D)
        end if
        call aprq3x(cmo,tmp1,tmp2,trn1(1,iocc),
     &              nbas,nvec,nocc,nvir,mxshl,nijshl,iocc)
c
      end do
      call end_time_period(TP_APRQ2D)
      return
      end
c
c  third transformation for vovo class
c
      subroutine aprq3x(cmo,tmp1,tmp2,tmp3,
     &                  nbas,nvec,nocc,nvir,mxshl,nijshl,iocc)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension cmo(nbas,nvec),tmp1(nijshl,nbas),tmp2(nvir,mxshl),
     &          tmp3(nvir,mxshl)
INCLUDE(common/global)
INCLUDE(common/misc)
INCLUDE(common/shlnos)
INCLUDE(common/indez)
INCLUDE(common/timeperiods)
c
      ltmp2=nvir*(maxi-mini+1)
      ltmp3=nvir*(maxj-minj+1)
      do jocc=1,iocc
        call dcopy(ltmp2,0.0d00,0,tmp2,1)
        call dcopy(ltmp3,0.0d00,0,tmp3,1)
        ijn = 0
        jmax = maxj
        ii=0
        do i = mini,maxi
          ii=ii+1
          if (oianj) jmax = i
          jj=0
          do j = minj,jmax
            jj=jj+1
            ijn = ijn+1
            if (i+loci.ne.j+locj)then
              do ibas=1,nvir
               tmp2(ibas,ii)=tmp2(ibas,ii)+
     &                       tmp1(ijn,nocc+ibas)*cmo(j+locj,jocc)
               tmp3(ibas,jj)=tmp3(ibas,jj)+
     &                       tmp1(ijn,nocc+ibas)*cmo(i+loci,jocc)
              enddo
            else
              do ibas=1,nvir
               tmp2(ibas,ii)=tmp2(ibas,ii)+
     &                       tmp1(ijn,nocc+ibas)*cmo(j+locj,jocc)
              enddo
            endif
          enddo
        enddo
c       call start_time_period(TP_GA_ACC_Q2D)
        io=(loci+mini-1)*nvir+1
        jo=(locj+minj-1)*nvir+1
        kl=iocc*(iocc-1)/2+jocc
        call pg_acc(g_vovo,io,io+ltmp2-1,kl,kl,tmp2,1,1.0d00)
        call pg_acc(g_vovo,jo,jo+ltmp3-1,kl,kl,tmp3,1,1.0d00)
c       call end_time_period(TP_GA_ACC_Q2D)
      enddo
      return
      end
c
c  second transformation in-core vvvo version
c
      subroutine aprq2(cmo,trn1,tmp1
     &,                nbas,nvec,nocc,nvir,mxshlt,nijshl,noffij)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension cmo(nbas,nvec),trn1(nbas,mxshlt,nocc),tmp1(nijshl,nbas)
INCLUDE(common/global)
INCLUDE(common/timeperiods)
INCLUDE(common/mp2grd)
      character*1 xn, xt
      data xn,xt/'n','t'/
      call start_time_period(TP_APRQ2)
      do iocc=1,nocc
       call dgemm(xt,xn,nijshl,nvec,nbas
     &   ,1.0d00,trn1(1,1,iocc),nbas
     &   ,cmo,nbas,0.0d00,tmp1,nijshl)
c
       if (opg_grad) then
         do ibas=1,nocc
          kl=(iocc-1)*nocc+ibas
          call pg_put(g_vvoo,noffij+1,noffij+nijshl,kl,kl,
     &                tmp1(1,ibas),nijshl)
         enddo
         do ibas=1,nvir
          kl=(iocc-1)*nvir+ibas
          call pg_put(g_vvvo,noffij+1,noffij+nijshl,kl,kl,
     &                tmp1(1,nocc+ibas),nijshl)
         end do
       end if
c
      end do
      call end_time_period(TP_APRQ2)
      return
      end
c
c  third transformation for vovo class in-core vvvo version
c
      subroutine oaprq3x(cmo,tmp1,tmp2,tmp3,nbas,nocc,nvir,mxshl,
     &                  nijshl,iocc)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension cmo(nbas,nbas),tmp1(nijshl,nbas),tmp2(nocc,mxshl),
     &          tmp3(nocc,mxshl)
INCLUDE(common/global)
INCLUDE(common/misc)
INCLUDE(common/shlnos)
INCLUDE(common/indez)
INCLUDE(common/timeperiods)
c
      ltmp2=nocc*(maxi-mini+1)
      ltmp3=nocc*(maxj-minj+1)
      do ibas=1,nvir
       call dcopy(ltmp2,0.0d00,0,tmp2,1)
       call dcopy(ltmp3,0.0d00,0,tmp3,1)
       ijn = 0
       jmax = maxj
       ii=0
       do i = mini,maxi
        ii=ii+1
        if (oianj) jmax = i
        jj=0
        do j = minj,jmax
         jj=jj+1
         ijn = ijn+1
         val = tmp1(ijn,nocc+ibas)
         if (i+loci.eq.j+locj)val=val*0.5d00
         do jocc=1,nocc
          tmp2(jocc,ii)=tmp2(jocc,ii)+val*cmo(j+locj,jocc)
          tmp3(jocc,jj)=tmp3(jocc,jj)+val*cmo(i+loci,jocc)
         enddo
        enddo
       enddo
       io=(loci+mini-1)*nocc+1
       jo=(locj+minj-1)*nocc+1
       kl=(iocc-1)*nvir+ibas
c       print *,'tmp2 ',kl,io,ltmp2,tmp2
c       print *,'tmp3 ',kl,jo,ltmp3,tmp3
c       call start_time_period(TP_GA_ACC_Q2D)
       call pg_acc(g_vovo,io,io+ltmp2-1,kl,kl,tmp2,1,1.0d00)
       call pg_acc(g_vovo,jo,jo+ltmp3-1,kl,kl,tmp3,1,1.0d00)
c       call end_time_period(TP_GA_ACC_Q2D)
      enddo
c
      return
      end
_EXTRACT(aprq34d,linux)
c
c  fourth transformation for vovo class + 
c  third and fourth transformations for vvoo, vooo, oooo classes
c
      subroutine aprq34d(cmo,tmp1,tmp2,nbas,nvec,nocc,nvir,nshels)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension cmo(nbas,nvec),tmp1(*),tmp2(nbas*nbas)
INCLUDE(common/sizes)
INCLUDE(common/nshel)
INCLUDE(common/global)
INCLUDE(common/timeperiods)
INCLUDE(common/mp2grd)
      character*1 xn, xt
      data xn,xt/'n','t'/
c
      call start_time_period(TP_APRQ34D)
c
c - [nn|oo]
c      print *,'q34 nnoo'
      nbt=nbas*(nbas+1)/2
      do iocc=1,nocc
       do jocc=1,iocc
        kl=iocc*(iocc-1)/2+jocc
        if(kl.ge.jl_vovo .and. kl.le.jh_vovo)then

        if (opg_grad) then 

        call pg_get(g_vvoo,1,nbt,kl,kl,tmp1,nbt)
        ijn = 0
        do ii = 1 , nshels
         do jj = 1 , ii
          do i = kloc(ii),kloc(ii)-kmin(ii)+kmax(ii)
           jmax = kloc(jj)-kmin(jj)+kmax(jj)
           if (ii.eq.jj) jmax = i
           do j = kloc(jj),jmax
            ijn = ijn+1
            tmp2((i-1)*nbas+j)=tmp1(ijn)
            tmp2((j-1)*nbas+i)=tmp1(ijn)
           enddo
          enddo
         enddo
        enddo
        call dgemm(xn,xn,nbas,nvec,nbas
     &    ,1.0d00,tmp2,nbas
     &    ,cmo,nbas,0.0d00,tmp1,nbas)
        call dgemm(xt,xn,nbas,nbas,nvec
     &    ,1.0d00,cmo,nbas
     &    ,tmp1,nbas,0.0d00,tmp2,nbas)
        ijn=0
        do i=1,nocc
         do j=1,i
          ijn=ijn+1
          tmp1(ijn)=tmp2((i-1)*nbas+j)
         enddo
        enddo
        call pg_put(g_oooo,1,ijn,kl,kl,tmp1,ijn)
c
        ijn=0
        do i=nocc+1,nvec
         do j=nocc+1,i
          ijn=ijn+1
          tmp1(ijn)=tmp2((i-1)*nbas+j)
         enddo
        enddo
        call pg_put(g_vvoo,1,ijn,kl,kl,tmp1,ijn)
c
        ijn=0
        do j=1,nocc
         do i=nocc+1,nvec
          ijn=ijn+1
          tmp1(ijn)=tmp2((i-1)*nbas+j)
         enddo
        enddo
        call pg_put(g_vooo,1,nvir,(kl-1)*nocc+1,kl*nocc,tmp1,nvir)
c
        end if

        call pg_get(g_vovo,1,nbas*nvir,kl,kl,tmp1,nbas*nvir)
        call dgemm(xn,xn,nvir,nvir,nbas
     &    ,1.0d00,tmp1,nvir
     &    ,cmo(1,nocc+1),nbas,0.0d00,tmp2,nvir)
        call pg_put(g_vovo,1,nvir*nvir,kl,kl,tmp2,nvir*nvir)

        endif
       enddo
      enddo
c
      call end_time_period(TP_APRQ34D)

      return
      end
_ENDEXTRACT
c
c  fourth transformation for vovo class +
c  third and fourth transformations for vvoo, vooo, oooo classes
c  in-core vvvo version
c
      subroutine aprq34(cmo,tmp1,tmp2,nbas,nvec,nocc,nvir,nshels)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension cmo(nbas,nvec),tmp1(*),tmp2(nbas,nbas)
INCLUDE(common/sizes)
INCLUDE(common/nshel)
INCLUDE(common/global)
INCLUDE(common/timeperiods)
      character*1 xn, xt
      data xn,xt/'n','t'/
c

      call start_time_period(TP_APRQ34)

c - [nn|oo]
      nbt=nbas*(nbas+1)/2
      do iocc=1,nocc
       do jocc=1,iocc
        kl=iocc*(iocc-1)/2+jocc

        if(kl.ge.jl_vvoo .and. kl.le.jh_vvoo)then
        call pg_get(g_vvoo,1,nbt,kl,kl,tmp1,nbt)
        ijn = 0
        do ii = 1 , nshels
         do jj = 1 , ii
          do i = kloc(ii),kloc(ii)-kmin(ii)+kmax(ii)
           jmax = kloc(jj)-kmin(jj)+kmax(jj)
           if (ii.eq.jj) jmax = i
           do j = kloc(jj),jmax
            ijn = ijn+1
            tmp2(i,j)=tmp1(ijn)
            tmp2(j,i)=tmp1(ijn)
           enddo
          enddo
         enddo
        enddo
        call dgemm(xn,xn,nbas,nvec,nbas
     &  ,1.0d00,tmp2,nbas
     &  ,cmo,nbas,0.0d00,tmp1,nbas)
        call dgemm(xt,xn,nbas,nbas,nvec
     &   ,1.0d00,cmo,nbas
     &   ,tmp1,nbas,0.0d00,tmp2,nbas)
        ijn=0
        do i=1,nocc
         do j=1,i
          ijn=ijn+1
          tmp1(ijn)=tmp2(i,j)
         enddo
        enddo
        call pg_put(g_oooo,1,ijn,kl,kl,tmp1,ijn)
c
        ijn=0
        do i=nocc+1,nvec
         do j=nocc+1,i
          ijn=ijn+1
          tmp1(ijn)=tmp2(i,j)
         enddo
        enddo
        call pg_put(g_vvoo,1,ijn,kl,kl,tmp1,ijn)
c
        endif
       enddo
      enddo
c
      do iocc=1,nocc
       do jocc=1,iocc
        kl=(iocc-1)*nvir+jocc
        lk=(jocc-1)*nocc+iocc

        if(kl.ge.jl_vvvo .and. kl.le.jh_vvvo)then

        call pg_get(g_vvvo,1,nbt,kl,kl,tmp1,nbt)
        ijn = 0
        do ii = 1 , nshels
         do jj = 1 , ii
          do i = kloc(ii),kloc(ii)-kmin(ii)+kmax(ii)
           jmax = kloc(jj)-kmin(jj)+kmax(jj)
           if (ii.eq.jj) jmax = i
           do j = kloc(jj),jmax
            ijn = ijn+1
            tmp2(i,j)=tmp1(ijn)
            tmp2(j,i)=tmp1(ijn)
           enddo
          enddo
         enddo
        enddo
        call dgemm(xn,xn,nbas,nvec,nbas
     &   ,1.0d00,tmp2,nbas
     &   ,cmo,nbas,0.0d00,tmp1,nbas)
        call dgemm(xt,xn,nbas,nbas,nvec
     &   ,1.0d00,cmo,nbas
     &   ,tmp1,nbas,0.0d00,tmp2,nbas)
        ijn=0
        do i=nocc+1,nvec
         do j=nocc+1,nvec
          ijn=ijn+1
          tmp1(ijn)=tmp2(i,j)
         enddo
        enddo
        call pg_put(g_vvvo,1,ijn,kl,kl,tmp1,ijn)
c
        ijn=0
        do i=1,nocc
         do j=nocc+1,nvec
          ijn=ijn+1
          tmp1(ijn)=tmp2(i,j)
         enddo
        enddo
        call pg_put(g_vovo,1,ijn,kl,kl,tmp1,ijn)
c
        ijn=0
        do i=1,nocc
         do j=1,nocc
          ijn=ijn+1
          tmp1(ijn)=tmp2(i,j)
         enddo
        enddo
        call pg_put(g_vooo,1,ijn,lk,lk,tmp1,ijn)
c
        endif
       enddo
      enddo
c
      call end_time_period(TP_APRQ34)
c
      return
      end
c
c  compute mp2 energy 
c
      subroutine apremp2(nocc,nvir,eig,tmp1,tmp2,emp2)
      implicit REAL (a-h,o-z)
      integer nocc,nvir
      REAL eig(nocc+nvir),tmp1(nvir,nvir),tmp2(nvir,nocc)
INCLUDE(common/global)
INCLUDE(common/timeperiods)
      REAL emp2
      integer i,j,a,b
      logical opg_root
      character*10 charwall
c
      call start_time_period(TP_APRMP2E)
c
      if(opg_root())then
       write(*,10) cpulft(1) ,charwall()
      endif
      nvv=nvir*nvir
      emp2=0.0d00
c      emp2x=0.0d00
      ijn=0
      do i=1,nocc
         do j=1,i
            ijn=ijn+1
            if(ijn.ge.jl_vovo .and. ijn.le.jh_vovo)then
               call pg_get(g_vovo, 1, nvv, ijn, ijn, tmp1, nvv)
c      print *,'i,j,ijn ',i,j,ijn
c      print *,'tmp1 ',tmp1
               do a=1,nvir
                  do b=1,nvir
c      xiajb=aprmoint(i,nocc+a,j,nocc+b)
c      xibja=aprmoint(i,nocc+b,j,nocc+a)
c      print *,a,b,xiajb,xibja,tmp(a,b),tm
                     if (i.ne.j)then
c                       emp2x=emp2x+2.0d00*xiajb*(xiajb+xiajb-xibja)/
c     &                 (eig(i)+eig(j)-eig(nocc+a)-eig(nocc+b))
                        emp2=emp2+2.0d00*tmp1(a,b)*(tmp1(a,b)+tmp1(a,b)
     &              -tmp1(b,a))/(eig(i)+eig(j)-eig(nocc+a)-eig(nocc+b))
                     else
c                       emp2x=emp2x+xiajb*(xiajb+xiajb-xibja)/
c     &                 (eig(i)+eig(j)-eig(nocc+a)-eig(nocc+b))
                        emp2=emp2+tmp1(a,b)*(tmp1(a,b)+tmp1(a,b)
     &              -tmp1(b,a))/(eig(i)+eig(j)-eig(nocc+a)-eig(nocc+b))
                     endif
                  enddo
               enddo
            endif
         enddo
      enddo

      call pg_dgop(999,emp2,1,'+')
      if(opg_root())then
       write(*,20) emp2
      endif
      call end_time_period(TP_APRMP2E)
c
      return
 10   format (/1x,'commence evaluation of mp2 energy at ',
     +         f10.2,' seconds',a10,' wall')
 20   format (/1x,'mp2 energy = ',f15.8)
      end
      subroutine aprdisp(nocc,nvir,eig,tmp1,tmp2,emp2)
c
c  compute mp2 coulomb dispersion energy 
c  assumes vector set non-orthogonal and ordered occa,occb,virta,virtb
c
      implicit REAL (a-h,o-z)
      integer nocc,nvir
      REAL eig(nocc+nvir),tmp1(nvir,nvir),tmp2(nvir,nocc)
INCLUDE(common/global)
INCLUDE(common/timeperiods)
INCLUDE(common/disp)
INCLUDE(common/iofile)
      REAL emp2
      integer i,j,a,b
      logical opg_root
      character*10 charwall
c
      call start_time_period(TP_APRMP2E)
c
      nvirt_a = nmo_a - nocc_a 
      if(opg_root())then
       write(iwr,10) cpulft(1) ,charwall()
      endif
      nvv=nvir*nvir
      emp2=0.0d00
c      emp2x=0.0d00
      ijn=0
      do i=1,nocc
      do j=1,i
      ijn=ijn+1
      if(ijn.ge.jl_vovo .and. ijn.le.jh_vovo)then
      call pg_get(g_vovo, 1, nvv, ijn, ijn, tmp1, nvv)
c      print *,'i,j,ijn ',i,j,ijn
c      print *,'tmp1 ',tmp1
      do a=1,nvir
      do b=1,nvir
c      xiajb=aprmoint(i,nocc+a,j,nocc+b)
c      xibja=aprmoint(i,nocc+b,j,nocc+a)
c      print *,a,b,xiajb,xibja,tmp(a,b),tm
      if (a.gt.nvirt_a.and.b.le.nvirt_a.and.
     1    i.gt.nocc_a.and.j.le.nocc_a) then
         fac = 2.0d0
         if (i.eq.j) fac = 1.0d0
c      emp2x=emp2x+2.0d00*xiajb*(xiajb+xiajb)/
c     &          (eig(i)+eig(j)-eig(nocc+a)-eig(nocc+b))
         emp2=emp2+fac*tmp1(a,b)*(tmp1(a,b)+tmp1(a,b))/
     &             (eig(i)+eig(j)-eig(nocc+a)-eig(nocc+b))
      end if
      enddo
      enddo
      endif
      enddo
      enddo

      call pg_dgop(999,emp2,1,'+')
      if(opg_root())then
         write (iwr,6021) nocc_a,nmo_a-nocc_a,
     1                    nocc_b,nmo_b-nocc_b, emp2,
     1                    emp2*2625.562d0,emp2*627.52d0,emp2*315777.d0
      endif
      call end_time_period(TP_APRMP2E)
c
      return
 10   format (/1x,'commence evaluation of mp2 energy at ',
     +         f10.2,' seconds',a10,' wall')
 6021 format (/10x,62('*')/10x,'mp2 Dispersion','  occ_a ',i3,
     +        ' virt_a ',i5,' occ_b ',i3,' virt_b ',i4,/10x,62('*'),
     +   /10x,'Polarisation dispersion  energy ',e15.8,' a.u.',
     +   /10x,'                                ',f15.10,' Kjoule/mole',
     +   /10x,'                                ',f15.10,' Kcal/mole',
     +   /10x,'                                ',f15.10,' K',
     +   /10x,62('*'))
      end
c
      subroutine aprshel(gout,nelec,ish,jsh,ksh,lsh,iexch)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      parameter (mxp2 = mxprms * mxprms)
INCLUDE(common/infoa)
INCLUDE(common/nshel)
INCLUDE(common/root)
      common/junk/cxyz(3,5625),aaa(21*mxp2),ijaaa(225)
INCLUDE(common/shlnos)
INCLUDE(common/indez)
INCLUDE(common/shlinf)
INCLUDE(common/misc)
INCLUDE(common/flips)
      dimension gout(*)
      dimension ix(35),iy(35),iz(35),jx(35),jy(35),jz(35),
     +          kx(35),ky(35),kz(35),lx(35),ly(35),lz(35)
      data lx /   0,  1,  0,  0,  2,  0,  0,  1,  1,  0,
     +            3,  0,  0,  2,  2,  1,  0,  1,  0,  1,
     +            4,  0,  0,  3,  3,  1,  0,  1,  0,  2,
     +            2,  0,  2,  1,  1/
      data kx /   0,  5,  0,  0, 10,  0,  0,  5,  5,  0,
     +           15,  0,  0, 10, 10,  5,  0,  5,  0,  5,
     +           20,  0,  0, 15, 15,  5,  0,  5,  0, 10,
     +           10,  0, 10,  5,  5/
      data jx /   0, 25,  0,  0, 50,  0,  0, 25, 25,  0,
     +           75,  0,  0, 50, 50, 25,  0, 25,  0, 25,
     +          100,  0,  0, 75, 75, 25,  0, 25,  0, 50,
     +           50,  0, 50, 25, 25/
      data ix /   1,126,  1,  1,251,  1,  1,126,126,  1,
     +          376,  1,  1,251,251,126,  1,126,  1,126,
     +          501,  1,  1,376,376,126,  1,126,  1,251,
     +          251,  1,251,126,126/
      data ly /   0,  0,  1,  0,  0,  2,  0,  1,  0,  1,
     +            0,  3,  0,  1,  0,  2,  2,  0,  1,  1,
     +            0,  4,  0,  1,  0,  3,  3,  0,  1,  2,
     +            0,  2,  1,  2,  1/
      data ky /   0,  0,  5,  0,  0, 10,  0,  5,  0,  5,
     +            0, 15,  0,  5,  0, 10, 10,  0,  5,  5,
     +            0, 20,  0,  5,  0, 15, 15,  0,  5, 10,
     +            0, 10,  5, 10,  5/
      data jy /   0,  0, 25,  0,  0, 50,  0, 25,  0, 25,
     +            0, 75,  0, 25,  0, 50, 50,  0, 25, 25,
     +            0,100,  0, 25,  0, 75, 75,  0, 25, 50,
     +            0, 50, 25, 50, 25/
      data iy /   1,  1,126,  1,  1,251,  1,126,  1,126,
     +            1,376,  1,126,  1,251,251,  1,126,126,
     +            1,501,  1,126,  1,376,376,  1,126,251,
     +            1,251,126,251,126/
      data lz /   0,  0,  0,  1,  0,  0,  2,  0,  1,  1,
     +            0,  0,  3,  0,  1,  0,  1,  2,  2,  1,
     +            0,  0,  4,  0,  1,  0,  1,  3,  3,  0,
     +            2,  2,  1,  1,  2/
      data kz /   0,  0,  0,  5,  0,  0, 10,  0,  5,  5,
     +            0,  0, 15,  0,  5,  0,  5, 10, 10,  5,
     +            0,  0, 20,  0,  5,  0,  5, 15, 15,  0,
     +           10, 10,  5,  5, 10/
      data jz /   0,  0,  0, 25,  0,  0, 50,  0, 25, 25,
     +            0,  0, 75,  0, 25,  0, 25, 50, 50, 25,
     +            0,  0,100,  0, 25,  0, 25, 75, 75,  0,
     +           50, 50, 25, 25, 50/
      data iz /   1,  1,  1,126,  1,  1,251,  1,126,126,
     +            1,  1,376,  1,126,  1,126,251,251,126,
     +            1,  1,501,  1,126,  1,126,376,376,  1,
     +          251,251,126,126,251/
c
      if (nelec.eq.2) then
         okanl = ksh.eq.lsh
capr         oident = ish.eq.ksh .and. jsh.eq.lsh
         oident = .false.
         ngtk = kgt(iexch)
         ngtl = lgt(iexch)

c
c     ----- kshell
c
         k = katom(ksh)
         pk = c(1,k)
         qk = c(2,k)
         rk = c(3,k)
         k1 = kstart(ksh)
         k2 = k1 + kng(ksh) - 1
         lkt = ktype(ksh)
         mink = kmin(ksh)
         maxk = kmax(ksh)
         lock = kloc(ksh) - mink
         ngc = 0
         do 20 k = k1 , k2
            ngc = ngc + 1
            cgg(ngc) = ex(k)
            csc(ngc) = cs(k)
            cpc(ngc) = cp(k)
            cdc(ngc) = cd(k)
            cfc(ngc) = cf(k)
            cgc(ngc) = cg(k)
 20      continue
c
c     ----- lshell
c
         l = katom(lsh)
         pl = c(1,l)
         ql = c(2,l)
         rl = c(3,l)
         l1 = kstart(lsh)
         l2 = l1 + kng(lsh) - 1
         llt = ktype(lsh)
         minl = kmin(lsh)
         maxl = kmax(lsh)
         locl = kloc(lsh) - minl
         ngd = 0
         do 30 l = l1 , l2
            ngd = ngd + 1
            dg(ngd) = ex(l)
            csd(ngd) = cs(l)
            cpd(ngd) = cp(l)
            cdd(ngd) = cd(l)
            cfd(ngd) = cf(l)
            cgd(ngd) = cg(l)
 30      continue
         nroots = (lit+ljt+lkt+llt-4)/2 + 1
         rrk = ((pk-pl)**2+(qk-ql)**2+(rk-rl)**2)
c
c     ----- prepare indices for pairs of (k,l) functions
c
         kl = 0
         lmax = maxl
         do 50 k = mink , maxk
            nnx = kx(k)
            nny = ky(k)
            nnz = kz(k)
            if (okanl) lmax = k
            do 40 l = minl , lmax
               kl = kl + 1
               klx(kl) = nnx + lx(l)
               kly(kl) = nny + ly(l)
               klz(kl) = nnz + lz(l)
               klgt(kl) = ngtk*(k-mink) + ngtl*(l-minl)
 40         continue
 50      continue
         maxkl = kl
         do 60 i = 1 , ij
            if (oident) maxkl = i
            ik(i) = maxkl
 60      continue
         ijkl = ij*kl
         if (oident) ijkl = ij*(ij+1)/2
c
c     zero integral storage
c
         do 65 i = 1,ij
         ngij = ijgt(i)
         do 65 k = 1,ik(i)
         gout(ngij + klgt(k)) = 0.0d0
  65     continue

         return
      else
         oianj = ish.eq.jsh
         ngti = igt(iexch)
         ngtj = jgt(iexch)
c
c     ----- ishell
c
         i = katom(ish)
         pi = c(1,i)
         qi = c(2,i)
         ri = c(3,i)
         i1 = kstart(ish)
         i2 = i1 + kng(ish) - 1
         lit = ktype(ish)
         mini = kmin(ish)
         maxi = kmax(ish)
         loci = kloc(ish) - mini
         nga = 0
         do 70 i = i1 , i2
            nga = nga + 1
            ag(nga) = ex(i)
            csa(nga) = cs(i)
            cpa(nga) = cp(i)
            cda(nga) = cd(i)
            cfa(nga) = cf(i)
            cga(nga) = cg(i)
 70      continue
c
c     ----- jshell
c
         j = katom(jsh)
         pj = c(1,j)
         qj = c(2,j)
         rj = c(3,j)
         j1 = kstart(jsh)
         j2 = j1 + kng(jsh) - 1
         ljt = ktype(jsh)
         minj = kmin(jsh)
         maxj = kmax(jsh)
         locj = kloc(jsh) - minj
         ngb = 0
         do 80 j = j1 , j2
            ngb = ngb + 1
            bg(ngb) = ex(j)
            csb(ngb) = cs(j)
            cpb(ngb) = cp(j)
            cdb(ngb) = cd(j)
            cfb(ngb) = cf(j)
            cgb(ngb) = cg(j)
 80      continue

         rri = ((pi-pj)**2+(qi-qj)**2+(ri-rj)**2)
c
c     ----- prepare indices for pairs of (i,j) functions
c
         ij = 0
         jmax = maxj
         do 100 i = mini , maxi
            nnx = ix(i)
            nny = iy(i)
            nnz = iz(i)
            if (oianj) jmax = i
            do 90 j = minj , jmax
               ij = ij + 1
               ijx(ij) = nnx + jx(j)
               ijy(ij) = nny + jy(j)
               ijz(ij) = nnz + jz(j)
               ijgt(ij) = ngti*(i-mini) + ngtj*(j-minj) + 1 
 90         continue
 100     continue

         return
      end if
      end
c
c  create global arrays
c
      subroutine alloc_ga2(no, nv, nb)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
INCLUDE(common/global)
INCLUDE(common/iofile)
INCLUDE(common/parcntl)
INCLUDE(common/mp2grd)
INCLUDE(common/gmempara)
      logical pg_create_inf
      character *14 fnm
      character *9 snm
      data fnm/'mp2_parallel.m'/
      data snm/'alloc_ga2'/
c
      nov = no*nv
      nbo = nb*no
      noo = no*(no+1)/2
      nooo = noo*no
      nvv = nv*nv
      ntri = nb*(nb+1)/2
      nsq = nb*nb
      mynode = ipg_nodeid()
      g_vvvo=-1
      call pg_synch(99)
c
capr - change here
c      if(opg_root())then
c        write(6,*)'need ',nsq*noo2
c        call ma_summarize_allocated_blocks
c      endif
c  [VO|VO]
      if (pg_create_inf(0,nsq,noo,'vovo',nsq,1,g_vovo,fnm,snm,
     +                  IGMEM_NORMAL)) then
         call pg_distribution(g_vovo, mynode, ilo, ihi, jlo, jhi)
         if(iparapr.eq.3)
     +      write(iwr,918) g_vovo, mynode, ilo, ihi, jlo, jhi
        il_vovo = ilo
        ih_vovo = ihi
        jl_vovo = jlo
        jh_vovo = jhi
      else
         call pg_error('**GA error** failed to create GA ',g_vovo)
      end if
c      if(opg_root())then
c        write(6,*)'need ',noo2*nov
c        call ma_summarize_allocated_blocks
c      endif

      if (opg_grad) then
c  [VO|OO]
         if (pg_create_inf(0,nv,nooo,'vooo',nv,1,g_vooo,fnm,snm,
     +                     IGMEM_NORMAL)) then
            call pg_distribution(g_vooo, mynode, ilo, ihi, jlo, jhi)
            if(iparapr.eq.3)
     +      write(iwr,918)  g_vooo, mynode, ilo, ihi, jlo, jhi
            il_vooo = ilo
            ih_vooo = ihi
            jl_vooo = jlo
            jh_vooo = jhi
         else
            call pg_error('**GA error** failed to create GA ',g_vooo)
         end if
c      if(opg_root())then
c        write(6,*)'need ',noo2*ntri
c        call ma_summarize_allocated_blocks
c      endif
c  [VV|OO]
         if (pg_create_inf(0,ntri,noo,'vvoo',ntri,1,g_vvoo,fnm,snm,
     +                     IGMEM_NORMAL)) then
            call pg_distribution(g_vvoo, mynode, ilo, ihi, jlo, jhi)
            if(iparapr.eq.3)
     +      write(iwr,918)  g_vvoo, mynode, ilo, ihi, jlo, jhi
            il_vvoo = ilo
            ih_vvoo = ihi
            jl_vvoo = jlo
            jh_vvoo = jhi
         else
            call pg_error('**GA error** failed to create GA ',g_vvoo)
         end if
c      if(opg_root())then
c        write(6,*)'need',noo2*noo2
c        call ma_summarize_allocated_blocks
c      endif
c  [OO|OO]
         if (pg_create_inf(0,noo,noo,'oooo',noo,1,g_oooo,fnm,snm,
     +                     IGMEM_NORMAL)) then
            call pg_distribution(g_oooo, mynode, ilo, ihi, jlo, jhi)
            if(iparapr.eq.3)
     +      write(iwr,918)  g_oooo, mynode, ilo, ihi, jlo, jhi
            il_oooo = ilo
            ih_oooo = ihi
            jl_oooo = jlo
            jh_oooo = jhi
         else
            call pg_error('**GA error** failed to create GA ',g_oooo)
         end if
c  [VV|VO]
         idvvvo=1
         if (idvvvo.eq.0)then
            write(iwr,*)
     &       ' ****vvvo**** integrals stored in global array '
            if (pg_create_inf(0,nsq,nov,'vvvo',nsq,1,g_vvvo,fnm,snm,
     +                        IGMEM_NORMAL)) then
               call pg_distribution(g_vvvo, mynode, ilo, ihi, jlo, jhi)
               if(iparapr.eq.3)
     +         write(iwr,918)  g_vvvo, mynode, ilo, ihi, jlo, jhi
               il_vvvo = ilo
               ih_vvvo = ihi
               jl_vvvo = jlo
               jh_vvvo = jhi
            else
              call pg_error(
     &         '**GA error** failed to create GA ',g_vvvo)
            end if
         else
           write(iwr,*)' ****vvvo**** term done direct '
         endif
      endif

c  initialize global arrays to zero
      call pg_synch(99)
      call zero_ga(g_vovo,il_vovo,ih_vovo,jl_vovo,jh_vovo)
      if (opg_grad) then
        call zero_ga(g_vooo,il_vooo,ih_vooo,jl_vooo,jh_vooo)
        call zero_ga(g_vvoo,il_vvoo,ih_vvoo,jl_vvoo,jh_vvoo)
        call zero_ga(g_oooo,il_oooo,ih_oooo,jl_oooo,jh_oooo)
        if (idvvvo.eq.0)then
          call zero_ga(g_vvvo,il_vvvo,ih_vvvo,jl_vvvo,jh_vvvo)
        end if
      end if
      call pg_synch(99)
      return
918   format(1x,'GA: distribution of GA [',1i6,'] to node ',1i4/
     +       1x,'GA: ilo =',1i6,3x,'ihi =',1i6,3x,
     +              'jlo =',1i6,3x,'jhi =',1i6)
      end
c  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   LAGRANGIAN THREE-VIRTUAL TERMS SECTION
c
c  initialisation routine for mp2 lagrangian three-virtual terms
c
      subroutine apr1pdm(wmo,pmo,cmo,eorb,schwa)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      REAL wmo(num,num),eorb(num),pmo(num,num),cmo(num,num)
      REAL schwa(*)
INCLUDE(common/sizes)
INCLUDE(common/fsymas)
INCLUDE(common/modj)
INCLUDE(common/cslosc)
INCLUDE(common/symtry)
INCLUDE(common/nshel)
INCLUDE(common/atmol3)
INCLUDE(common/prints)
INCLUDE(common/scra7)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/mapper)
INCLUDE(common/timez)
INCLUDE(common/infoa)
INCLUDE(common/scfopt)
INCLUDE(common/restar)
INCLUDE(common/segm)
      common/diisd/st(210),cdiis(20),rdiis(19),derror,sdiis(20),
     + iposit(20),nsti(2),ondiis,junkj
INCLUDE(common/runlab)
      common/restri/jjfile(63),lds(508),isect(508),ldsect(508)
INCLUDE(common/psscrf)
INCLUDE(common/timeperiods)
_IF(parallel)
      common/nodeio/ids48(maxlfn),idr48(maxlfn),oputpp(maxlfn)
     + , maxio(maxlfn),oswed3(maxlfn)
_ENDIF
INCLUDE(common/global)
INCLUDE(common/vcore)
c      logical opg_root
      data zscf/'scf'/
      data dzero,two,pt5/0.0d0,2.0d0,0.5d0/
      data yblnk,yav/' ',' av'/
      data done,pt2,twopt2 /1.0d0,0.2d0,2.2d0/
      data dmptlc/1.0d-2/
      data igs/5/
      out = nprint .eq. 5
      nav = lenwrd()
_IF(parallel)
      oswed3(4) = .true.
      oswed3(8) = .true.
_ENDIF
c gdf:
      call start_time_period(TP_APR1PDM)
c
      nbas=num
      nocc=na
      nvir=nbas-nocc
      if(nprint.ne.-5)write (iwr,9348)
cps
c      write (iwr,9148)
c      call prev(cmo,eorb,nbas,nbas,nbas)
c
      idebug=0
      if (idebug.ne.0)call dcopy(nbas*nbas,0.0d00,0,wmo,1)
      if (idebug.gt.10)then
       do i=1,nbas
        do j=1,nbas
         cmo(i,j)=0.0d00
        enddo
        cmo(i,i)=1.0d00
       enddo
      endif
c
c ----------------------------------------
c  Back transform t2 amplitudes
c ----------------------------------------
      noo=nocc*nocc
      ntri=nbas*(nbas+1)/2

      i00 = igmem_alloc(nbas*nbas)
      i10 = igmem_alloc(nbas*nbas)
      i20 = igmem_alloc(nbas*nbas)

      call mkt2ao(Q(i00),Q(i10),
     &            cmo,eorb,nocc, nvir, nbas,pmo,idebug)

      call gmem_free(i20)
      call gmem_free(i10)

c ------------------------------------------------
c Contract AO t2 amplitudes with [vvvo] integrals
c ------------------------------------------------
      ostat = .true.
c     mxshl=0
      mxshl=4
      do i=1,nshell
c        mxshl = max(mxshl,kmax(i)-kmin(i)+1)
         if(ktype(i).eq.3) then
          mxshl = max(mxshl,6)
         else if(ktype(i).eq.4) then
          mxshl = max(mxshl,10)
         else if(ktype(i).eq.5) then
          mxshl = max(mxshl,15)
         else
         endif
      enddo
      mxshl=max(mxshl,4)

      i40 = igmem_alloc(nw196(5))
      i50 = igmem_alloc(mxshl**4)
      i60 = igmem_alloc(mxshl*mxshl*nbas*nocc*2)
      i80 = igmem_alloc(nocc*nocc)
      i90 = igmem_alloc(nocc*nbas)
      i95 = igmem_alloc(nocc*(nocc+1)/2)
      call rdedx(Q(i40),nw196(5),ibl196(5),idaf)

      call aprl1234(Q(i40),schwa,Q(i50),
     &              nshell,cmo,Q(i60),
     &              Q(i80),Q(i90),nocc,nvir,nbas,
     &              mxshl,eorb,Q(i00),wmo,idebug,
     &	            Q(i95))
      call gmem_free(i95)
      call gmem_free(i90)
      call gmem_free(i80)
      call gmem_free(i60)
      call gmem_free(i50)
      call gmem_free(i40)
      call gmem_free(i00)
c
c ------------------------------------------------
c  transform amplitudes to MO basis
c ------------------------------------------------

      i00 = igmem_alloc(nbas*nbas)
      i10 = igmem_alloc(nbas*nbas)
      i20 = igmem_alloc(nbas*nbas)
 
      call mkt2mo(Q(i00),Q(i10),Q(i20),
     &          cmo,eorb,nocc,nvir,nbas)

      call gmem_free(i20)
      call gmem_free(i10)
      call gmem_free(i00)

      call end_time_period(TP_APR1PDM)
      return
c
 9348 format(/40x,32('*')/
     +        40x,'direct lagrangian formation'/
     +        40x,32('*'))
 9148 format(//1x,104('-')//
     *  50x,12('-')/50x,'eigenvectors'/50x,12('-'))
      end
c
c  parallel integral driver
c
      subroutine aprl1234(iso,schwa,gout,nshels,cmo,trn1,t2ao,lag,nocc,
     &                    nvir,nbas,mxshl,eorb,pao,wmo,idebug,buf)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension gout(*),schwa(*),iso(nshels,*)
      REAL  cmo(nbas*nbas),trn1(nbas,nocc,mxshl,mxshl,2),
     +       t2ao(nocc,nocc),lag(nbas*nocc),eorb(nbas),pao(nbas,nbas),
     +       wmo(nbas,nbas),buf(*)
      logical sym_pair, sym_quartet_l
      logical osym, osym_pair
c
INCLUDE(common/sizes)
      parameter (mxp2 = mxprms * mxprms)
INCLUDE(common/cslosc)
INCLUDE(common/mapper)
INCLUDE(common/restar)
INCLUDE(common/iofile)
INCLUDE(common/utilc)
INCLUDE(common/nshel)
INCLUDE(common/shlnos)
INCLUDE(common/shlt)
INCLUDE(common/symtry)
INCLUDE(common/shlg70)
INCLUDE(common/picon)
INCLUDE(common/ijlab)
c
c     ----- size of gout -
c                         1   if s or k shells
c                        81   if p      shells
c                       256   if      l shells
c                      1296   if d or m shells
c                     10000   if f shells
c                     50625   if g shells
c
c     ----- this version can handle g shells   -----
c
INCLUDE(common/timez)
INCLUDE(common/pkfil)
INCLUDE(common/flips)
INCLUDE(common/sortp)
INCLUDE(common/parallel)
INCLUDE(common/timeperiods)
      dimension mi(48),mj(48),mk(48),m0(48)
c
      dimension ib(4,4)
      data ib/64,16,4,1,216,36,6,1,1000,100,10,1,
     +         3375,225,15,1 /
      data done/1.0d0/,two/2.0d0/,twopt5/2.5d0/,four/4.0d0/
c
      call start_time_period(TP_APRL1234)
c
      dlncutoff = dlog(cutoff)
      q4 = 1.0d0
      osym = .true.
      osym_pair = .true.
      pidiv4 = datan(done)
      pi = four*pidiv4
      pito52 = two*pi**twopt5
c
c     ----- two-electron integrals -----
c
c     zero lagrangian to start with
      call dcopy(nbas*nocc,0.0d00,0,lag,1)
      ii = 1
      if (opdbas) then
         ii = 2
      end if
      if (opfbas) then
         ii = 3
      end if
      if (opgbas) then
         ii = 4
      end if
      do 30 loop = 1 , 3
        igt(loop) = ib(1,ii)
        jgt(loop) = ib(2,ii)
        kgt(loop) = ib(3,ii)
        lgt(loop) = ib(4,ii)
 30   continue
c
      nschwz = 0
c
      ist0 = 1
      jst0 = 1
      kst0 = 1
      lst0 = 1
c
      if (intg76.ne.0) call sinset
c
      ltrn1=nocc*nbas*mxshl*mxshl*2
      call pg_dlbreset()
      mytask=ipg_dlbtask()
      nproc= ipg_nnodes()
      me=ipg_nodeid()
c
      if (intg76.ne.0) call filmax
c
c     are the rotated axis integrals to be used only?
c
      ogauss = intg76.eq.0. or. opdbas .or. opfbas
c
c ----- ishell -----
c
      do 150 ii = ist0 , nshels
       ishell = ii
       mini = kmin(ishell)
       maxi = kmax(ishell)
       loci = kloc(ishell)-mini
       kadi = kad(ii)
c
c ----- kshell -----
c
       do 130 kk = kst0 , ii
         kst0 = 1
         kshell = kk
         mink = kmin(kshell)
         maxk = kmax(kshell)
         lock = kloc(kshell)-mink
         kadik = kadi+kad(kk)
c
c  symmetry (IK| 
c
        if (nt.gt.1 ) then
         osym_pair = sym_pair(iso,nshels,nt,ii,kk)
        endif
        if (osym_pair) then
c
c  load balancer
c
        icount_dlb = icount_dlb + 1
        if (icount_dlb.eq.mytask)then
         mytask=ipg_dlbtask()
c
         call dcopy(ltrn1,0.0d00,0,trn1,1)
c
c ----- jshell -----
c
         do 140 jj = jst0 , nshels
          jst0 = 1
          jshell = jj
          minj = kmin(jshell)
          maxj = kmax(jshell)
          locj = kloc(jshell)-minj
          kadijk = kadik+kad(jj)
          if (ogauss) then
           call aprshel2(gout,1,ishell,jshell,kshell,lshell,1)
           call ijprim
           if (nij.eq.0) go to 170
          endif
c
c ----- lshell ----
c
           do 120 ll = lst0 , nshels
            lst0 = 1
            lshell = ll
            minl = kmin(lshell)
            maxl = kmax(lshell)
            locl = kloc(lshell)-minl
            itrij = iky(max(ii,jj)) + min(ii,jj)
            itrkl = iky(max(kk,ll)) + min(kk,ll)
c
c  schwarz inequality test
c
            test = schwa(itrij) + schwa(itrkl)
            if (test.lt.dlncutoff) then
              nschwz = nschwz + 1
            else
c
c  symmetry (IJ|KL)
c
             if (nt.gt.1) then
              osym = sym_quartet_l(iso,nshels,nt,ii,jj,kk,ll,q4)
             endif
             if (osym) then
              qq4 = q4
              if (kadijk+kad(ll).lt.0)then
               call aprshel2(gout,2,ishell,jshell,kshell,lshell,1)
c        call start_time_period(TP_GENRAL_1PDM)
c
c  compute (IJ|KL)
c
               call genral(gout)
c        call end_time_period(TP_GENRAL_1PDM)
               call aprl1(gout,cmo,trn1,nbas,nocc,mxshl,idebug)
              else
               call genr70(gout,1,.false.)
               call aprl1s(gout,cmo,trn1,nbas,nocc,mxshl,idebug)
              endif
             endif
            endif
 120       continue
 170      continue
 140     continue
         call aprl2(t2ao,trn1,lag,nbas,nocc,mxshl,pao,idebug,buf)
         endif
        endif
 130   continue
 150  continue
      call pg_dlbreset()
      call aprl34(cmo,lag,wmo,nocc,nvir,nbas)
c  ... schwarz inequality log
      if(nprint.ne.-5) write(iwr,6040) nschwz
      call end_time_period(TP_APRL1234)
      return
6040  format(/1x,'schwarz inequality test skipped ',i10,
     + ' integral blocks')
      end
      logical function sym_quartet_l(iso,nshels,nops
     +,     i, j, k, l, q4)
c  constituency numbers adjusted for permutational symmetries
c  employed in vvvo terms
      implicit REAL  (a-h,o-z)
      integer i, j, k, l        ! shell indices [input]
      REAL q4                   ! Constituency number [output]
      integer iso(nshels,*), nshels
      integer i_new, j_new, k_new, l_new
      integer nops, op , ns2, indx, indx_new
      integer n                 ! Counts no. of equivalent pairs
INCLUDE(common/sizes)
INCLUDE(common/mapper)
      q4 = 0.0d0
      sym_quartet_l = .false.
      ns2 = nshels*nshels
      if (k.gt.i) return                            
      indx = (iky(i) + k -1)*ns2 + (j-1)*nshels + l      
c     Loop thru operations in the group and map to new pairs
c
      n = 1                     ! Identity always counts
      do op = 2, nops
c
c     Map centers
c
       i_new = iso( i, op)
       j_new = iso( j, op)
       k_new = iso( k, op)
       l_new = iso( l, op)
       if (i_new.ge.k_new) then
         indx_new = (iky(i_new) + k_new - 1)*ns2    
     &              + (j_new - 1)*nshels + l_new          
       else 
         indx_new = (iky(k_new) + i_new - 1)*ns2    
     &              + (l_new - 1)*nshels + j_new 
       end if
       if (indx_new.gt.indx) return               
       if (indx_new.eq.indx) n = n + 1           
      end do
c
      q4 = dfloat(nops) / dfloat(n)
      if (dabs(q4-nint(q4)).gt.1d-12) call caserr
     &     ('sym_quartet_l: not divisible')
      sym_quartet_l = .true.
c
        return
        end
c
c  transformation of lagrangian
c
      subroutine aprl34(cmo,lag,wmo,nocc,nvir,nbas)
      implicit REAL (a-h,o-z)
      REAL lag
INCLUDE(common/timeperiods)
      dimension lag(nbas*nocc),wmo(nbas,nbas),cmo(nbas,nbas)
c
      character*1 xn, xt
      data xn,xt/'n','t'/

      call start_time_period(TP_APRL34)
      call dgemm(xt,xn,nvir,nocc,nbas
     &    ,1.0d00,cmo(1,nocc+1),nbas
     &    ,lag,nbas,1.0d00,wmo(nocc+1,1),nbas)
      call end_time_period(TP_APRL34)
c
      return
      end
c
c  first transformation
c
      subroutine aprl1(g,cmo,trn1,nbas,nocc,mxshl,idebug)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension g(*),cmo(nbas,nbas),trn1(nbas,nocc,mxshl,mxshl,2)
INCLUDE(common/misc)
INCLUDE(common/shlnos)
INCLUDE(common/indez)
INCLUDE(common/timeperiods)
c      REAL xtst1
c
cc      call start_time_period(TP_APRL1)
      ijn = 0
      ix=0
      do i = loci+mini, loci+maxi
       ix=ix+1
       do j = locj+minj, locj+maxj
        ijn = ijn+1
        kx=0
        kln=0
        do k = lock+mink,lock+maxk
         kx=kx+1
         do l = locl+minl,locl+maxl
          kln = kln+1
          nn = ijgt(ijn)+klgt(kln)
          val = g(nn)
c          xtst1=apraoint(i,j,k,l)
c          if (dabs(val-xtst1).gt.1.0d-8)then
c           write(*,843)i,j,k,l,nn,val,xtst1,val-xtst1
c  843      format(5i4,3e16.6)
c          endif
          do iocc=1,nocc
           trn1(j,iocc,ix,kx,1)=trn1(j,iocc,ix,kx,1)+val*cmo(l,iocc)
           trn1(l,iocc,kx,ix,2)=trn1(l,iocc,kx,ix,2)+val*cmo(j,iocc)
          enddo
         enddo
        enddo
       enddo
      enddo
c
cc      call end_time_period(TP_APRL1)
      return
      end
c
c  first transformation fast integrals version 
c
      subroutine aprl1s(g,cmo,trn1,nbas,nocc,mxshl,idebug)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension g(*),cmo(nbas,nbas),trn1(nbas,nocc,mxshl,mxshl,2)
INCLUDE(common/sizes)
INCLUDE(common/nshel)
INCLUDE(common/misc)
INCLUDE(common/shlnos)
INCLUDE(common/shlg70)
      common/flip70/ib1,jb1,kb1,lb1,ib2,jb2,kb2,lb2,ib3,jb3,kb3,lb3
INCLUDE(common/timeperiods)
      dimension ib(4,4)
      data ib/0,0,0,0,64,16,4,1,128,32,8,2,192,48,12,3/
c
      ix = 0
      do ii = mini, maxi
       ix = ix + 1
       do jj = minj, maxj
        j = jj + locj
        n1 = ib(ib1,ii)+ib(jb1,jj)+1
        kx = 0
        do kk = mink,maxk
         kx = kx + 1
         do ll = minl,maxl
          l = ll + locl
          nn = n1 + ib(kb1,kk)+ib(lb1,ll)
          val = g(nn)
          do iocc=1,nocc
           trn1(j,iocc,ix,kx,1)=trn1(j,iocc,ix,kx,1)+val*cmo(l,iocc)
           trn1(l,iocc,kx,ix,2)=trn1(l,iocc,kx,ix,2)+val*cmo(j,iocc)
          enddo
         enddo
        enddo
       enddo
      enddo
c
      return
      end
c
c  second transformation
c
      subroutine aprl2(t2ao,trn1,lag,nbas,nocc,mxshl,pao,
     &                 idebug,buf)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      REAL t2ao(nocc,nocc),trn1(nbas,nocc,mxshl,mxshl,2),
     &       lag(nbas,nocc),pao(nbas,nbas),buf(*)
INCLUDE(common/misc)
INCLUDE(common/shlnos)
INCLUDE(common/global)
INCLUDE(common/timeperiods)
c      REAL xtst1,xtst2,xtst3,xtst4
      character*1 xn, xt
      data xn,xt/'n','t'/
c
      if (idebug.ne.0)then
       print *,'cant check q1 transformation '
cccc       ix=0
cccc       do i = loci+mini, loci+maxi
cccc        ix=ix+1
cccc        kx=0
cccc        do k = lock+mink, lock+maxk
cccc         kx=kx+1
cccc         do lo=1,nocc
cccc          do j=1,nbas
cccc           xtst1=trn1(j,lo,ix,kx,1)
cccc           xtst2=aprq1int(i,j,k,lo)
cccc           xtst3=trn1(j,lo,kx,ix,2)
cccc           xtst4=aprq1int(k,j,i,lo)
cccc           if (dabs(xtst1-xtst2).gt.1.0d-8.or.
cccc     &         dabs(xtst3-xtst4).gt.1.0d-8)then
cccc           write(*,43)j,lo,i,k,xtst1,xtst2,xtst1-xtst2,
cccc     &                xtst3,xtst4,xtst3-xtst4
cccc   43      format(4i5,3e16.6/,20x,3e16.6)
cccc           endif
cccc          enddo
cccc         enddo
cccc        enddo
cccc       enddo
cccc       print *,'end q1 check '
      endif

        call start_time_period(TP_APRL2)

      ix=0
      do i = loci+mini, loci+maxi
       ix=ix+1
       kx=0
       do k = lock+mink, lock+maxk
        kx=kx+1
        if (i.ge.k)then
           call start_time_period(TP_GA_GET_L2)
           call get_half_trans(t2ao,buf,nbas,nocc,i,k)
           call end_time_period(TP_GA_GET_L2)
         call dgemm(xn,xt,nbas,nocc,nocc
     & ,1.0d00,trn1(1,1,ix,kx,1)
     & ,nbas,t2ao,nocc,1.0d00,lag,nbas)
         do lo=1,nocc
          do j=1,nbas
           lag(k,lo)=lag(k,lo)+4.0d00*trn1(j,lo,ix,kx,1)*pao(i,j)
           lag(i,lo)=lag(i,lo)-trn1(j,lo,ix,kx,1)*pao(k,j)
           lag(j,lo)=lag(j,lo)-trn1(j,lo,ix,kx,1)*pao(i,k)
          enddo
         enddo
         if (i.ne.k)then
         call dgemm(xn,xn,nbas,nocc,nocc
     &     ,1.0d00,trn1(1,1,kx,ix,2)
     &     ,nbas,t2ao,nocc,1.0d00,lag,nbas)
         do lo=1,nocc
          do j=1,nbas
           lag(i,lo)=lag(i,lo)+4.0d00*trn1(j,lo,kx,ix,2)*pao(k,j)
           lag(k,lo)=lag(k,lo)-trn1(j,lo,kx,ix,2)*pao(i,j)
           lag(j,lo)=lag(j,lo)-trn1(j,lo,kx,ix,2)*pao(k,i)
          enddo
         enddo
         endif
        endif
       enddo
      enddo

        call end_time_period(TP_APRL2)

      return
      end
c
      subroutine aprshel2(gout,nelec,ish,jsh,ksh,lsh,iexch)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      parameter (mxp2 = mxprms * mxprms)
INCLUDE(common/infoa)
INCLUDE(common/nshel)
INCLUDE(common/root)
      common/junk/cxyz(3,5625),aaa(21*mxp2),ijaaa(225)
INCLUDE(common/shlnos)
INCLUDE(common/indez)
INCLUDE(common/shlinf)
INCLUDE(common/misc)
INCLUDE(common/flips)
      dimension gout(*)
      dimension ix(35),iy(35),iz(35),jx(35),jy(35),jz(35),
     +          kx(35),ky(35),kz(35),lx(35),ly(35),lz(35)
      data lx /   0,  1,  0,  0,  2,  0,  0,  1,  1,  0,
     +            3,  0,  0,  2,  2,  1,  0,  1,  0,  1,
     +            4,  0,  0,  3,  3,  1,  0,  1,  0,  2,
     +            2,  0,  2,  1,  1/
      data kx /   0,  5,  0,  0, 10,  0,  0,  5,  5,  0,
     +           15,  0,  0, 10, 10,  5,  0,  5,  0,  5,
     +           20,  0,  0, 15, 15,  5,  0,  5,  0, 10,
     +           10,  0, 10,  5,  5/
      data jx /   0, 25,  0,  0, 50,  0,  0, 25, 25,  0,
     +           75,  0,  0, 50, 50, 25,  0, 25,  0, 25,
     +          100,  0,  0, 75, 75, 25,  0, 25,  0, 50,
     +           50,  0, 50, 25, 25/
      data ix /   1,126,  1,  1,251,  1,  1,126,126,  1,
     +          376,  1,  1,251,251,126,  1,126,  1,126,
     +          501,  1,  1,376,376,126,  1,126,  1,251,
     +          251,  1,251,126,126/
      data ly /   0,  0,  1,  0,  0,  2,  0,  1,  0,  1,
     +            0,  3,  0,  1,  0,  2,  2,  0,  1,  1,
     +            0,  4,  0,  1,  0,  3,  3,  0,  1,  2,
     +            0,  2,  1,  2,  1/
      data ky /   0,  0,  5,  0,  0, 10,  0,  5,  0,  5,
     +            0, 15,  0,  5,  0, 10, 10,  0,  5,  5,
     +            0, 20,  0,  5,  0, 15, 15,  0,  5, 10,
     +            0, 10,  5, 10,  5/
      data jy /   0,  0, 25,  0,  0, 50,  0, 25,  0, 25,
     +            0, 75,  0, 25,  0, 50, 50,  0, 25, 25,
     +            0,100,  0, 25,  0, 75, 75,  0, 25, 50,
     +            0, 50, 25, 50, 25/
      data iy /   1,  1,126,  1,  1,251,  1,126,  1,126,
     +            1,376,  1,126,  1,251,251,  1,126,126,
     +            1,501,  1,126,  1,376,376,  1,126,251,
     +            1,251,126,251,126/
      data lz /   0,  0,  0,  1,  0,  0,  2,  0,  1,  1,
     +            0,  0,  3,  0,  1,  0,  1,  2,  2,  1,
     +            0,  0,  4,  0,  1,  0,  1,  3,  3,  0,
     +            2,  2,  1,  1,  2/
      data kz /   0,  0,  0,  5,  0,  0, 10,  0,  5,  5,
     +            0,  0, 15,  0,  5,  0,  5, 10, 10,  5,
     +            0,  0, 20,  0,  5,  0,  5, 15, 15,  0,
     +           10, 10,  5,  5, 10/
      data jz /   0,  0,  0, 25,  0,  0, 50,  0, 25, 25,
     +            0,  0, 75,  0, 25,  0, 25, 50, 50, 25,
     +            0,  0,100,  0, 25,  0, 25, 75, 75,  0,
     +           50, 50, 25, 25, 50/
      data iz /   1,  1,  1,126,  1,  1,251,  1,126,126,
     +            1,  1,376,  1,126,  1,126,251,251,126,
     +            1,  1,501,  1,126,  1,126,376,376,  1,
     +          251,251,126,126,251/
c
      if (nelec.eq.2) then
capr         okanl = ksh.eq.lsh
capr         oident = ish.eq.ksh .and. jsh.eq.lsh
         okanl = .false.
         oident = .false.
         ngtk = kgt(iexch)
         ngtl = lgt(iexch)
c
c     ----- kshell
c
         k = katom(ksh)
         pk = c(1,k)
         qk = c(2,k)
         rk = c(3,k)
         k1 = kstart(ksh)
         k2 = k1 + kng(ksh) - 1
         lkt = ktype(ksh)
         mink = kmin(ksh)
         maxk = kmax(ksh)
         lock = kloc(ksh) - mink
         ngc = 0
         do 20 k = k1 , k2
            ngc = ngc + 1
            cgg(ngc) = ex(k)
            csc(ngc) = cs(k)
            cpc(ngc) = cp(k)
            cdc(ngc) = cd(k)
            cfc(ngc) = cf(k)
            cgc(ngc) = cg(k)
 20      continue
c
c     ----- lshell
c
         l = katom(lsh)
         pl = c(1,l)
         ql = c(2,l)
         rl = c(3,l)
         l1 = kstart(lsh)
         l2 = l1 + kng(lsh) - 1
         llt = ktype(lsh)
         minl = kmin(lsh)
         maxl = kmax(lsh)
         locl = kloc(lsh) - minl
         ngd = 0
         do 30 l = l1 , l2
            ngd = ngd + 1
            dg(ngd) = ex(l)
            csd(ngd) = cs(l)
            cpd(ngd) = cp(l)
            cdd(ngd) = cd(l)
            cfd(ngd) = cf(l)
            cgd(ngd) = cg(l)
 30      continue
         nroots = (lit+ljt+lkt+llt-4)/2 + 1
         rrk = ((pk-pl)**2+(qk-ql)**2+(rk-rl)**2)
c
c     ----- prepare indices for pairs of (k,l) functions
c
         kl = 0
         lmax = maxl
         do 50 k = mink , maxk
            nnx = kx(k)
            nny = ky(k)
            nnz = kz(k)
            if (okanl) lmax = k
            do 40 l = minl , lmax
               kl = kl + 1
               klx(kl) = nnx + lx(l)
               kly(kl) = nny + ly(l)
               klz(kl) = nnz + lz(l)
               klgt(kl) = ngtk*(k-mink) + ngtl*(l-minl)
 40         continue
 50      continue
         maxkl = kl
         do 60 i = 1 , ij
            if (oident) maxkl = i
            ik(i) = maxkl
 60      continue
         ijkl = ij*kl
         if (oident) ijkl = ij*(ij+1)/2
c
c     zero integral storage
c
         do 65 i = 1,ij
         ngij = ijgt(i)
         do 65 k = 1,ik(i)
         gout(ngij + klgt(k)) = 0.0d0
  65     continue
         return
      else
capr         oianj = ish.eq.jsh
         oianj = .false.
         ngti = igt(iexch)
         ngtj = jgt(iexch)
c
c     ----- ishell
c
         i = katom(ish)
         pi = c(1,i)
         qi = c(2,i)
         ri = c(3,i)
         i1 = kstart(ish)
         i2 = i1 + kng(ish) - 1
         lit = ktype(ish)
         mini = kmin(ish)
         maxi = kmax(ish)
         loci = kloc(ish) - mini
         nga = 0
         do 70 i = i1 , i2
            nga = nga + 1
            ag(nga) = ex(i)
            csa(nga) = cs(i)
            cpa(nga) = cp(i)
            cda(nga) = cd(i)
            cfa(nga) = cf(i)
            cga(nga) = cg(i)
 70      continue
c
c     ----- jshell
c
         j = katom(jsh)
         pj = c(1,j)
         qj = c(2,j)
         rj = c(3,j)
         j1 = kstart(jsh)
         j2 = j1 + kng(jsh) - 1
         ljt = ktype(jsh)
         minj = kmin(jsh)
         maxj = kmax(jsh)
         locj = kloc(jsh) - minj
         ngb = 0
         do 80 j = j1 , j2
            ngb = ngb + 1
            bg(ngb) = ex(j)
            csb(ngb) = cs(j)
            cpb(ngb) = cp(j)
            cdb(ngb) = cd(j)
            cfb(ngb) = cf(j)
            cgb(ngb) = cg(j)
 80      continue
         rri = ((pi-pj)**2+(qi-qj)**2+(ri-rj)**2)
c
c     ----- prepare indices for pairs of (i,j) functions
c
         ij = 0
         jmax = maxj
         do 100 i = mini , maxi
            nnx = ix(i)
            nny = iy(i)
            nnz = iz(i)
            if (oianj) jmax = i
            do 90 j = minj , jmax
               ij = ij + 1
               ijx(ij) = nnx + jx(j)
               ijy(ij) = nny + jy(j)
               ijz(ij) = nnz + jz(j)
               ijgt(ij) = ngti*(i-mini) + ngtj*(j-minj) + 1 
 90         continue
 100     continue
         return
      end if
      end
c
c  construct half-back-transformed amplitudes
c
      subroutine mkt2ao(tmp1,tmp3,cmo,eorb,nocc,nvir,
     &                  nbas,pmo,idebug)
      implicit integer (i-n)
      REAL tmp1(nbas*nbas),tmp3(nbas*nbas),
     &       cmo(nbas,nbas),eorb(nbas),pmo(nbas,nbas)
INCLUDE(common/global)
INCLUDE(common/timeperiods)
      integer av,bv,ab,ba
      character*1 xn, xt
      data xn,xt/'n','t'/

c
c - make t2 with virtual indices transformed to ao basis
      if (idebug.ne.0)print *,'mkt2ao '

      call start_time_period(TP_MKT2AO)

      nsq = nbas*nbas
      nov=nocc*nvir
      nvv=nvir*nvir
      ij = 0
      do i = 1 , nocc
       do j = 1 , i
        ij = ij + 1
       if (ij.ge.jl_vovo.and.ij.le.jh_vovo) then
        call pg_get(g_vovo, 1, nvv, ij, ij, tmp1, nvv)
        iabv=0
        do av = nocc+1 , nocc+nvir
         do bv = nocc+1, nocc+nvir
          ab = (bv-nocc-1)*nvir+av-nocc
          ba = (av-nocc-1)*nvir+bv-nocc
          iabv=iabv+1
          tmp3(iabv)=4.0d00*(tmp1(ab)+tmp1(ab)-tmp1(ba))/
     &               (eorb(av)+eorb(bv)-eorb(i)-eorb(j))
         enddo
        enddo
        call dgemm(xn,xn,nbas,nvir,nvir
     &    ,1.0d00,cmo(1,nocc+1),nbas
     &    ,tmp3,nvir,0.0d00,tmp1,nbas)
        call dgemm(xn,xt,nbas,nbas,nvir
     &  ,1.0d00,tmp1,nbas
     &  ,cmo(1,nocc+1),nbas,0.0d00,tmp3,nbas)
        call pg_put(g_vovo, 1, nsq, ij, ij, tmp3, 1)
       end if
       enddo
      enddo
c
c - transform pmo to ao basis
      call dgemm(xn,xn,nbas,nvir,nvir
     &    ,1.0d00,cmo(1,nocc+1),nbas
     &    ,pmo(nocc+1,nocc+1),nbas,0.0d00,tmp3,nbas)
      call dgemm(xn,xt,nbas,nbas,nvir
     &     ,1.0d00,tmp3,nbas
     &     ,cmo(1,nocc+1),nbas,0.0d00,tmp1,nbas)
c
        call end_time_period(TP_MKT2AO)
c
      return
      end
c
c  reverse mkt2ao - restore vovo integrals
c
      subroutine mkt2mo(tmp1,tmp2,cm1,cmo,eorb,nocc,nvir,nbas)
      implicit integer (i-n)
      REAL tmp1(nbas,nbas),tmp2(nbas*nbas),cm1(nbas,nbas),
     &  cmo(nbas,nbas),eorb(nbas)
INCLUDE(common/global)
INCLUDE(common/timeperiods)
INCLUDE(common/scra7)
INCLUDE(common/iofile)
      integer a,b,ab,ba
      REAL third, twothirds , one, zero, tmp2ab, tmp2ba
      data one,zero,third,twothirds/1.0d0,0.0d0
     &, 0.3333333333333333d0, 0.6666666666666667d0/
      character*1 xn, xt
      data xn,xt/'n','t'/
c
c - restore vovo integrals with ao indices transformed to virtual indices

      call start_time_period(TP_MKT2MO)

      nsq = nbas*nbas
      ntri = nbas*(nbas+1)/2
      nov=nocc*nvir
      nvv=nvir*nvir
      nocc1 = nocc + 1
c  form inverse MO coeff matrix:
c   get overlap matrix 
      call rdedx(tmp1,ntri,ibl7s,num8)
      call tr2sq_ga(tmp1,tmp2,nbas)
c   c(-1) = c(t)*s
      call dgemm(xt,xn,nbas,nbas,nbas
     &   ,one,cmo,nbas
     &   ,tmp2,nbas,zero,cm1,nbas)

      ij = 0 
      do i = 1 , nocc
       do j = 1 , i
        ij = ij + 1
       if (ij.ge.jl_vovo.and.ij.le.jh_vovo) then
        eij = eorb(i) + eorb(j)
        call pg_get(g_vovo, 1, nsq, ij, ij, tmp1, 1)
c  transform amplitudes to mo basis
        call dgemm(xn,xn,nbas,nbas,nbas
     &   ,one,cm1,nbas
     &   ,tmp1,nbas,zero,tmp2,nbas)
        call dgemm(xn,xt,nbas,nbas,nbas
     &   ,one,tmp2,nbas
     &   ,cm1,nbas,zero,tmp1,nbas)
c  convert amplitudes back to integrals
        iabv=0
        do a = nocc1, nbas
         do b = nocc1, nbas
          ab = (b-nocc1)*nvir + a-nocc
          ba = (a-nocc1)*nvir + b-nocc
          iabv=iabv+1
          denom = (eorb(a)+eorb(b)-eij)*0.25d0
          tmp2ab = tmp1(a,b)*denom
          tmp2ba = tmp1(b,a)*denom
          tmp2(ba) = twothirds*tmp2ab + third*tmp2ba
          tmp2(ab) = twothirds*tmp2ba + third*tmp2ab
         enddo
        enddo
        call pg_put(g_vovo, 1, nvv, ij, ij, tmp2, 1)
        endif
       enddo
      enddo
c
        call end_time_period(TP_MKT2MO)
c
      return
      end
c  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   1-PARTICLE DENSITY AND CPHF SECTION
c
c
c  construct terms of W and P
c
      subroutine aprwp(pmat,wmat,eorb,nocc,nvir,nbas,buf1,buf2,buf3)
      implicit none
      integer nocc,nvir,nbas
      REAL pmat(nbas,nbas),wmat(nbas,nbas),eorb(nbas),
     &       buf1(nbas*nbas),buf2(nbas*nbas),buf3(nbas*nbas)
INCLUDE(common/global)
INCLUDE(common/timeperiods)
c
      integer i,j,k,l,a,b,nop1,nvv,ab,ba,ij,ik,jk,ijk,ia,nov,
     &        noo,noo2,nvv2,kl,mytask,ipg_dlbtask
      integer ipg_nnodes,ipg_nodeid,nproc,me
      REAL xiajb,tiajb,diajb
      character*1 xn, xt
      data xn,xt/'n','t'/

c
      call start_time_period(TP_MP1PDM)
c
      nop1=nocc+1
      nvv=nvir*nvir
      nov=nvir*nocc
      noo=nocc*nocc
      noo2=nocc*(nocc+1)/2
      nvv2=nvir*(nvir+1)/2
      call vclr(wmat,1,nbas*nbas)
      call vclr(pmat,1,nbas*nbas)
c --------------------------------------------------
c 1st parallel segment
c --------------------------------------------------
c
c -- pvv1 and wvv1
      call start_time_period(TP_MP1PDM_1)
      ij=0
      do i=1,nocc
       do j=1,i
        ij=ij+1
        if (ij.ge.jl_vovo.and.ij.le.jh_vovo) then
         call pg_get(g_vovo,1,nvv,ij,ij,buf2,nvv)
         ba=0
         do a=1,nvir
          ab=a-nvir
          do b=1,nvir
           ba=ba+1
           ab=ab+nvir
           buf1(ba)=(buf2(ab)+buf2(ab)-buf2(ba))/
     &              (eorb(nocc+a)+eorb(nocc+b)-eorb(i)-eorb(j))
          enddo
         enddo
         call dgemm(xn,xn,nvir,nvir,nvir
     &    ,4.0d00,buf2,nvir
     &    ,buf1,nvir,1.0d00,wmat(nop1,nop1),nbas)
         if (i.ne.j)then
          call dgemm(xt,xt,nvir,nvir,nvir
     &   ,4.0d00,buf2,nvir
     &   ,buf1,nvir,1.0d00,wmat(nop1,nop1),nbas)
         endif
         ba=0
         do a=1,nvir
          ab=a-nvir
          do b=1,nvir
           ba=ba+1
           ab=ab+nvir
           buf2(ba)=buf2(ba)/
     &              (eorb(nocc+a)+eorb(nocc+b)-eorb(i)-eorb(j))
          enddo
         enddo
         call dgemm(xt,xt,nvir,nvir,nvir
     &   ,-2.0d00,buf1,nvir
     &   ,buf2,nvir,1.0d00,pmat(nop1,nop1),nbas)
         if (i.ne.j)then
          call dgemm(xn,xn,nvir,nvir,nvir
     &   ,-2.0d00,buf1,nvir
     &   ,buf2,nvir,1.0d00,pmat(nop1,nop1),nbas)
         endif
        endif
       enddo
      enddo
      call end_time_period(TP_MP1PDM_1)
c
c --------------------------------------------------
c 2nd parallel segment
c --------------------------------------------------
c -- poo1 and woo1
      call start_time_period(TP_MP1PDM_2)
      ij=0
      nproc=ipg_nnodes()
      me=ipg_nodeid()
      call pg_dlbreset()
      mytask=ipg_dlbtask()
      ijk = 0
      do i=1,nocc
       do j=1,i
       ij=ij+1
        do k=1,j
         ijk=ijk+1
         if(ijk.eq.mytask)then
          mytask=ipg_dlbtask()
          ik=i*(i-1)/2+k
          jk=j*(j-1)/2+k
          call pg_get(g_vovo,1,nvv,ij,ij,buf1,nvv)
          call pg_get(g_vovo,1,nvv,ik,ik,buf2,nvv)
          call pg_get(g_vovo,1,nvv,jk,jk,buf3,nvv)
c
c         call wpij(pmat,wmat,eorb,nbas,nocc,nvir,i,j,k)
          ba=0
          do a=1,nvir
           ab=a-nvir
           do b=1,nvir
            ba=ba+1
            ab=ab+nvir
            xiajb=buf1(ab)
            tiajb=(buf2(ab)+buf2(ab)-buf2(ba))/
     &            (eorb(nocc+a)+eorb(nocc+b)-eorb(i)-eorb(k))
            diajb=eorb(nocc+a)+eorb(nocc+b)-eorb(i)-eorb(j)
            pmat(k,j)=pmat(k,j)+2.0d0*xiajb*tiajb/diajb
            wmat(k,j)=wmat(k,j)+4.0d00*xiajb*tiajb
           enddo
          enddo
c 
          if (j.ne.k)then
c          call wpij(pmat,wmat,eorb,nbas,nocc,nvir,i,k,j)
           ba=0
           do a=1,nvir
            ab=a-nvir
            do b=1,nvir
             ba=ba+1
             ab=ab+nvir
             xiajb=buf2(ab)
             tiajb=(buf1(ab)+buf1(ab)-buf1(ba))/
     &             (eorb(nocc+a)+eorb(nocc+b)-eorb(i)-eorb(j))
             diajb=eorb(nocc+a)+eorb(nocc+b)-eorb(i)-eorb(k)
             pmat(j,k)=pmat(j,k)+2.0d0*xiajb*tiajb/diajb
             wmat(j,k)=wmat(j,k)+4.0d00*xiajb*tiajb
            enddo
           enddo
          endif
c 
          if (i.ne.j)then
c           call wpij(pmat,wmat,eorb,nbas,nocc,nvir,j,i,k)
           ba=0
           do a=1,nvir
            ab=a-nvir
            do b=1,nvir
             ba=ba+1
             ab=ab+nvir
             xiajb=buf1(ba)
             tiajb=(buf3(ab)+buf3(ab)-buf3(ba))/
     &             (eorb(nocc+a)+eorb(nocc+b)-eorb(j)-eorb(k))
             diajb=eorb(nocc+a)+eorb(nocc+b)-eorb(i)-eorb(j)
             pmat(k,i)=pmat(k,i)+2.0d0*xiajb*tiajb/diajb
             wmat(k,i)=wmat(k,i)+4.0d00*xiajb*tiajb
            enddo
           enddo
c
            if (i.ne.k)then
c             call wpij(pmat,wmat,eorb,nbas,nocc,nvir,j,k,i)
             ba=0
             do a=1,nvir
              ab=a-nvir
              do b=1,nvir
               ba=ba+1
               ab=ab+nvir
               xiajb=buf3(ab)
               tiajb=(buf1(ba)+buf1(ba)-buf1(ab))/
     &               (eorb(nocc+a)+eorb(nocc+b)-eorb(j)-eorb(i))
               diajb=eorb(nocc+a)+eorb(nocc+b)-eorb(k)-eorb(j)
               pmat(i,k)=pmat(i,k)+2.0d0*xiajb*tiajb/diajb
               wmat(i,k)=wmat(i,k)+4.0d00*xiajb*tiajb
              enddo
             enddo
            endif
          endif
c 
          if (k.ne.i.and.k.ne.j)then
c           call wpij(pmat,wmat,eorb,nbas,nocc,nvir,k,j,i)
           ba=0
           do a=1,nvir
            ab=a-nvir
           do b=1,nvir
             ba=ba+1
             ab=ab+nvir
             xiajb=buf3(ba)
             tiajb=(buf2(ba)+buf2(ba)-buf2(ab))/
     &             (eorb(nocc+a)+eorb(nocc+b)-eorb(k)-eorb(i))
             diajb=eorb(nocc+a)+eorb(nocc+b)-eorb(k)-eorb(j)
             pmat(i,j)=pmat(i,j)+2.0d0*xiajb*tiajb/diajb
             wmat(i,j)=wmat(i,j)+4.0d00*xiajb*tiajb
            enddo
           enddo
           if (i.ne.j)then
c            call wpij(pmat,wmat,eorb,nbas,nocc,nvir,k,i,j)
            ba=0
            do a=1,nvir
             ab=a-nvir
             do b=1,nvir
              ba=ba+1
              ab=ab+nvir
              xiajb=buf2(ba)
              tiajb=(buf3(ba)+buf3(ba)-buf3(ab))/
     &              (eorb(nocc+a)+eorb(nocc+b)-eorb(k)-eorb(j))
              diajb=eorb(nocc+a)+eorb(nocc+b)-eorb(k)-eorb(i)
              pmat(j,i)=pmat(j,i)+2.0d0*xiajb*tiajb/diajb
              wmat(j,i)=wmat(j,i)+4.0d00*xiajb*tiajb
             enddo
            enddo
           endif
          endif
c
         endif
        enddo
       enddo
      enddo
      call pg_dlbreset()
      call pg_dgop(1001,pmat,nbas*nbas,'+')
      call end_time_period(TP_MP1PDM_2)
c
c --------------------------------------------------
c 3rd parallel segment
c --------------------------------------------------
c -- wov1
      call start_time_period(TP_MP1PDM_3)
      ij=0
      do k=1,nocc
       do j=1,k
        ij=ij+1
        if (ij.ge.jl_vovo.and.ij.le.jh_vovo) then
         call pg_get(g_vovo,1,nvv,ij,ij,buf2,nvv)
         ba=0
         do a=1,nvir
          ab=a-nvir
          do b=1,nvir
           ba=ba+1
           ab=ab+nvir
           buf1(ba)=(buf2(ab)+buf2(ab)-buf2(ba))/
     &              (eorb(nocc+a)+eorb(nocc+b)-eorb(k)-eorb(j))
          enddo
         enddo
c
         ia=1
         do i=1,nocc
          if (i.ge.k)then
           jk=(i*(i-1)/2+k-1)*nocc+j
          else
           jk=(k*(k-1)/2+i-1)*nocc+j
          endif
          call pg_get(g_vooo,1,nvir,jk,jk,buf2(ia),nvir)
          ia=ia+nvir
         enddo
c -- wov1
         call dgemm(xt,xn,nocc,nvir,nvir
     &     ,4.0d00,buf2,nvir
     &     ,buf1,nvir,1.0d00,wmat(1,nop1),nbas)
c -- wvo2
         ia=0
         do i=1,nocc
          do a=nop1,nbas
           ia=ia+1
           xiajb=buf2(ia)
           wmat(a,j)=wmat(a,j)+xiajb*pmat(i,k)*4.0d00
           wmat(a,i)=wmat(a,i)-xiajb*pmat(j,k)
           wmat(a,k)=wmat(a,k)-xiajb*pmat(j,i)
          end do
         end do
         if (j.ne.k)then
          ia=1
          do i=1,nocc
           if (i.ge.j)then
            jk=(i*(i-1)/2+j-1)*nocc+k
           else
            jk=(j*(j-1)/2+i-1)*nocc+k
           endif
           call pg_get(g_vooo,1,nvir,jk,jk,buf2(ia),nvir)
           ia=ia+nvir
          enddo
c -- wov1
          call dgemm(xt,xt,nocc,nvir,nvir
     &   ,4.0d00,buf2,nvir
     &   ,buf1,nvir,1.0d00,wmat(1,nop1),nbas)
c -- wvo2
          ia=0
          do i=1,nocc
           do a=nop1,nbas
            ia=ia+1
            xiajb=buf2(ia)
            wmat(a,k)=wmat(a,k)+xiajb*pmat(i,j)*4.0d00
            wmat(a,i)=wmat(a,i)-xiajb*pmat(k,j)
            wmat(a,j)=wmat(a,j)-xiajb*pmat(k,i)
           end do
          end do
         endif
c
        endif
       enddo
      enddo
      call end_time_period(TP_MP1PDM_3)
c
c --------------------------------------------------
c 4th parallel segment
c --------------------------------------------------
      call start_time_period(TP_MP1PDM_4)
      ij=0
      do i=1,nocc
       do j=1,i
        ij=ij+1
        if (ij.ge.jl_oooo.and.ij.le.jh_oooo) then
         call pg_get(g_oooo,1,noo2,ij,ij,buf1,noo2)
         call pg_get(g_vvoo,1,nvv2,ij,ij,buf2,nvv2)
         call pg_get(g_vovo,1,nvv,ij,ij,buf3,nvv)
         do k=1,nocc
          do l=1,nocc
           kl=k*(k-1)/2+l
           if (k.lt.l)kl=l*(l-1)/2+k
           xiajb=buf1(kl)
           wmat(i,j)=wmat(i,j)+xiajb*pmat(k,l)*4.0d00
           wmat(i,k)=wmat(i,k)-xiajb*pmat(j,l)
           wmat(i,l)=wmat(i,l)-xiajb*pmat(j,k)
          enddo
        enddo
        do a=1,nvir
          do b=1,nvir
           ab=a*(a-1)/2+b
           if (a.lt.b)ab=b*(b-1)/2+a
           xiajb=buf2(ab)
           wmat(i,j)=wmat(i,j)+xiajb*pmat(a+nocc,b+nocc)*4.0d00
           ab=(b-1)*nvir+a
           xiajb=buf3(ab)
           wmat(i,j)=wmat(i,j)-xiajb*pmat(a+nocc,b+nocc)
           wmat(j,i)=wmat(j,i)-xiajb*pmat(a+nocc,b+nocc)
          enddo
        enddo
c
        if (i.ne.j)then
         do k=1,nocc
           do l=1,nocc
            kl=k*(k-1)/2+l
            if (k.lt.l)kl=l*(l-1)/2+k
            xiajb=buf1(kl)
            wmat(j,i)=wmat(j,i)+xiajb*pmat(k,l)*4.0d00
            wmat(j,k)=wmat(j,k)-xiajb*pmat(i,l)
            wmat(j,l)=wmat(j,l)-xiajb*pmat(i,k)
           enddo
         enddo
c -- woo3
         do a=1,nvir
          do b=1,nvir
            ab=a*(a-1)/2+b
            if (a.lt.b)ab=b*(b-1)/2+a
            xiajb=buf2(ab)
            wmat(j,i)=wmat(j,i)+xiajb*pmat(a+nocc,b+nocc)*4.0d00
            ab=(a-1)*nvir+b
            xiajb=buf3(ab)
            wmat(j,i)=wmat(j,i)-xiajb*pmat(a+nocc,b+nocc)
            wmat(i,j)=wmat(i,j)-xiajb*pmat(a+nocc,b+nocc)
           enddo
          enddo
         endif

        endif
       enddo
      enddo
      call end_time_period(TP_MP1PDM_4)
      call end_time_period(TP_MP1PDM)
c
      return
      end
c
c  mp2 1-particle density driving routine
c
      subroutine mpgrad(q)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
      dimension q(*)
      integer mpass,ipass,mpst,mpfi,mpadd,iaoxf1,iaoxf2,moderf
      common/mpases/mpass,ipass,mpst,mpfi,mpadd,iaoxf1,iaoxf2,moderf
INCLUDE(common/nshel)
INCLUDE(common/infoa)
INCLUDE(common/cigrad)
INCLUDE(common/symtry)
      common/maxlen/maxq
INCLUDE(common/common)
INCLUDE(common/iofile)
INCLUDE(common/cndx41)
INCLUDE(common/atmblk)
INCLUDE(common/dump3)
INCLUDE(common/global)
INCLUDE(common/mp2grad_pointers)
INCLUDE(common/parcntl)
INCLUDE(common/prnprn)
      integer i10, i20, i30, i40, i50, i60, i70, i80, i90, iz
      integer no, nv, n, nov, ntri, nsq, nvv, noo
      integer i, a, jwr, jw, iw, kz, itop, maxc
      integer isec9, iblok
      integer maxq
      integer igmem_alloc, ilo, ihi, jlo, jhi, ipg_nodeid, mynode
      logical pg_create
      integer lensec, isdd, ityp, len
      character *14 fnm
      character *6  snm
      data fnm/'mp2_parallel.m'/
      data snm/'mpgrad'/

      iconvv = 9
      mynode = ipg_nodeid()
      no = nocca
      nv = nvirta
      n = ncoorb
      nov = no*nv
      noo = no*no
      nvv = nv*nv
      ntri = n*(n+1)/2
      nsq = n*n

c  read in eigenvalues
      call secget(isect(9),9,isec9)
      call rdedx(q(mp2grad_vals),n,isec9,ifild)
c  read in eigenvectors 
      call secget(isect(8),8,iblok)
      iblok = iblok + mvadd
      call rdedx(q(mp2grad_vecs),nsq,iblok,ifild)
      if (oprn(15)) then
c       note that wrtmat_ga is limited to TINY examples
c       allocate scratch array for wrtmat_ga
c
        nspace = 100000
        i10 = igmem_alloc(nspace)
        call wrtmat_ga('(VO|VO)  ',g_vovo,q(i10))
        call wrtmat_ga('(VO|OO)  ',g_vooo,q(i10))
        call wrtmat_ga('(VV|OO)  ',g_vvoo,q(i10))
        call wrtmat_ga('(OO|OO)  ',g_oooo,q(i10))
c       call wrtmat_ga('(VV|VO)  ',g_vvvo,q(i10))
c       free scratch space
        call gmem_free(i10)
      end if
c
c       ----  construct Lagrangian rhs and blocks of P(2) and W(2)  ----
c
******
c     call pg_synch(99)
******
c
      i10 = igmem_alloc(nsq)
      i20 = igmem_alloc(nsq)
      i30 = igmem_alloc(nsq)
      call aprwp(q(mp2grad_pmat),q(mp2grad_wmat),q(mp2grad_vals),
     &             no,nv,n,q(i10),q(i20),q(i30))
*     if (oprn(15)) then
*      write(iwr,*) 'W(2) matrix after aprwp'
*      call wrtsqm_ga('W(2) matrix',n,q(mp2grad_wmat))
*     endif
      call gmem_free(i30)
      call gmem_free(i20)
      call gmem_free(i10)
****
**** debug mp2grad_wmat and mp2grad_pmat
****
c     i10 = igmem_alloc(nsq)
c     call dcopy(nsq,q(mp2grad_wmat),1,q(i10),1)
c     call pg_dgop(1002,q(i10),nsq,'+')
c     sssum = dsum(nsq,q(i10),1)
c     write(6,1289)sssum
c1289  format(' DSUM_wmat after aprwp = ',f20.10)
c     call dcopy(nsq,q(mp2grad_pmat),1,q(i10),1)
c     call pg_dgop(1002,q(i10),nsq,'+')
c     sssum = dsum(nsq,q(i10),1)
c     write(6,1287)sssum
c1287  format(' DSUM_pmat after aprwp = ',f20.10)
c     call gmem_free(i20)
c     call gmem_free(i10)
****
      call delete_ga_inf(g_oooo,'oooo',fnm,snm)
*     if (oprn(15)) call wrtsqm_ga('P(2) matrix',n,q(mp2grad_pmat))
      if (g_vvvo.ge.0)then
       call delete_ga_inf(g_vvvo,'vvvo',fnm,snm)
      else if (g_vvvo.lt.0)then
c  direct 3 virtual terms
       call apr1pdm(    q(mp2grad_wmat) 
     &,                 q(mp2grad_pmat)
     &,                 q(mp2grad_vecs)
     &,                 q(mp2grad_vals)
     &,                 q(mp2grad_schw))
      endif

c  global sum of W(2)
      call pg_dgop(1002,q(mp2grad_wmat),nsq,'+')
****
c     sssum = dsum(nsq,q(mp2grad_wmat),1)
c     write(6,1288)sssum
c1288  format(' DSUM after apr1pdm = ',f20.10)
****
*     if (oprn(15))call wrtsqm_ga('W(2) matrix',n,q(mp2grad_wmat))
c
c                      ----  solve for z-matrix  ----
c
c  iz => lagrangian rhs
c
      iz = igmem_alloc(nov)
      kz = iz
      iw = mp2grad_wmat + no
      jwr = mp2grad_wmat + no*n
      do i = 1 , no
        jw = jwr
        do a = 1 , nv
          q(kz) = q(iw) - q(jw)
          kz = kz + 1
          iw = iw + 1
          jw = jw + n
        end do
        iw = iw + no
        jwr = jwr + 1
      end do
c
c  zero blocks of lagrangian if vanishing by symmetry
c
      call sym_lagrangian(q(iz),no,n)
c
      maxc = 50
      if (pg_create(0,nov,maxc,'u',nov,1,g_u)) then
        call pg_distribution(g_u, mynode, ilo, ihi, jlo, jhi)
        if(iparapr.eq.3)
     +   write(iwr,918)  g_u, mynode, ilo, ihi, jlo, jhi
        il_u = ilo
        ih_u = ihi
        jl_u = jlo
        jh_u = jhi
        call zero_ga(g_u,ilo, ihi, jlo, jhi)
      else
        call pg_error('**GA error** failed to create GA ',g_u)
      end if
      i10 = igmem_alloc(nov)
      i20 = igmem_alloc(nov)
      i30 = igmem_alloc(nov)
      i40 = igmem_alloc(maxc)
      i50 = igmem_alloc(maxc)
      i60 = igmem_alloc(maxc)
      i70 = igmem_alloc(maxc*maxc)
      i80 = igmem_alloc(max( nvv , max( 2*maxc , nov )))
      i90 = igmem_alloc(nvv)

      call mp2chf_ga(   q(mp2grad_vals),  q(iz)
     &,q(i10),q(i20),q(i30),q(i40),q(i50),q(i60),q(i70),q(i80),q(i90)
     &,maxc,nov)

      call gmem_free(i90)
      call gmem_free(i80)
      call gmem_free(i70)
      call gmem_free(i60)
      call gmem_free(i50)
      call gmem_free(i40)
      call gmem_free(i30)
      call gmem_free(i20)
      call gmem_free(i10)
      call delete_ga(g_u)
      call delete_ga_inf(g_vvoo,'vvoo',fnm,snm)
c
c              ----  complete the W(2) and P(2) matrices  ----
c
      i10 = igmem_alloc(noo)
      i20 = igmem_alloc(max(noo,nv))
c
      call makew_ga(    q(mp2grad_pmat)
     &,                 q(iz)
     &,                 q(mp2grad_wmat)
     &,                 q(mp2grad_vals)
     &, q(i10),q(i20),no,nv,n)

*     if (oprn(15))call wrtsqm_ga('W(2) matrix',n,q(mp2grad_wmat))

      call gmem_free(i20)
      call gmem_free(i10)
      call gmem_free(iz)
      call delete_ga_inf(g_vooo,'vooo',fnm,snm)

      i10 = igmem_alloc(nsq)
c  back-transform P(2) to AO-basis
      call vtamv(       q(mp2grad_pmat)
     &,                 q(mp2grad_vecs),  q(i10),n)
c  symmetrise P(2) = q(i10) triangle ...
      call trsqsq(      q(mp2grad_pmat),  q(i10),n)
c  write out triangular p(2) for MP2 analysis 
c  mp2nat in anala
      ityp = 0      
      len = lensec(ntri)
      call secput(isecdd,ityp,len,isdd)
      call wrt3(q(i10),ntri,isdd,ifild)
c 
      call squr(q(i10), q(mp2grad_pmat),  n)
c  back-transform W(2) to AO-basis
      call vtamv(       q(mp2grad_wmat)
     &,                 q(mp2grad_vecs),  q(i10),n)
c  symmetrise W(2) = q(i10) triangle ...
      call trsqsq(      q(mp2grad_wmat),  q(i10),n)


      call squr(q(i10), q(mp2grad_wmat),  n)
c  make P(HF)
      call vclr(        q(mp2grad_dens),  1,nsq)
      call mxmb(        q(mp2grad_vecs)
     &,  1,n,           q(mp2grad_vecs)
     &,  n,1,           q(mp2grad_dens),  1,n,n,no,n)
      call gmem_free(i10)

*     if (oprn(15)) then
*       call wrtsqm_ga('W(2) matrix',n,q(mp2grad_wmat))
*       call wrtsqm_ga('P(2) matrix',n,q(mp2grad_pmat))
*       call wrtsqm_ga('eigenvectors',n,q(mp2grad_vecs))
*       call wrtsqm_ga('density matrix',n,q(mp2grad_dens))
*     end if
      return
918   format(1x,'GA: distribution of GA [',1i2,'] to node ',1i4/
     +       1x,'GA: ilo =',1i6,3x,'ihi =',1i6,3x,
     +              'jlo =',1i6,3x,'jhi =',1i6)
      end
c
c  zero blocks of lagrangian vanishing by symmetry
c   (do this for now prior to blocking cphf 
c    - keep as safeguard at least)
c
      subroutine sym_lagrangian(lag,no,n)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      REAL lag(*)
      integer i,a,ia, no,n, isymi, isyma
      character*1 zclass
*     common/crap_96/zclass(maxorb)
      integer mmmm, isymao, isymmo
      common/bufb/mmmm(65),isymao(maxorb),isymmo(maxorb)
      ia = 0
      do i = 1 , no
      isymi = isymmo(i) - 1
        do a = no + 1 , n
        isyma = isymmo(a) - 1
          ia = ia + 1
*         if (zclass(i).ne.zclass(a)) lag(ia) = 0.0d0
          if (ieor(isymi,isyma).eq.0 ) then
          else
          lag(ia) = 0.0d0
          endif
        end do
      end do
      return
      end
c
c  cphf solver
c  in-core global array version of chfeqv 
c
      subroutine mp2chf_ga(e, rhs
     &,     u,unxt,prhs
     &,     b,cc,uu,uau,buf,buf2,   maxc,nov)
      implicit REAL  (a-h,o-z)
INCLUDE(common/timez)
INCLUDE(common/common)
INCLUDE(common/cndx41)
INCLUDE(common/prnprn)
INCLUDE(common/iofile)
      integer a
      common/small/alpha(50,50),aa(50,50),wk1(50),wk2(50)
      dimension e(*), rhs(nov)
      dimension u(nov),unxt(nov),prhs(nov)
      dimension b(maxc),cc(maxc),uu(maxc),uau(maxc,maxc)
      dimension buf(*), buf2(*)
      data smal/1.0d-13/,tich/1.0d-24/
      data four/4.0d0/,zero/0.0d0/
INCLUDE(common/global)
INCLUDE(common/timeperiods)
      character*10 charwall
      no = nocca
      nv = nvirta
      n = ncoorb
      no1 = no + 1

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
c  get zeroth order estimate
      ia = 0
      do i = 1 , no
        do a = no1 , n
          ia = ia + 1
          v = rhs(ia)/(e(i) - e(a))
          if (dabs(v).le.tich) v = zero 
          b(1) = b(1) + v*v
          u(ia) = v
        end do
      end do

      if (1.ge.jl_u.and.1.le.jh_u)
     &  call pg_put(g_u, 1, nov, 1, 1, u, 1)

c      call wrtvec_ga('rhs   ',nov,u)

      if (oprn(6)) write(iwr,6100)
c
c  start of iterative solution of chf equations
c  50 iterations are allowed ---  usually less than 10 are necessary
c
      do 1 itr = 1 , maxc
        call vclr(unxt,1,nov)
        call mxmau_ga(u,unxt,buf,buf2,no,nv)
        call pg_dgop(101+itr,unxt,nov,'+')
c        call wrtvec_ga('unxt   ',nov,unxt)
c  scale unxt by difference of eigenvalues
        ia = 0
        do i = 1 , no
          do a = no1 , n
            ia = ia + 1
            unxt(ia) = unxt(ia)/(e(a) - e(i))
          end do
        end do
        uu(itr) = ddot(nov,u,1,u,1)

        do i = 2 , itr - 1
          if (i.ge.jl_u.and.i.le.jh_u) then
            call pg_get(g_u, 1, nov, i, i, buf, 1)
            uau(itr,i-1) = ddot(nov,u,1,buf,1)
          end if
        end do
        if (mynode.eq.nnodes-1.and.itr.gt.1)
     &     uau(itr,itr - 1) = uu(itr)
        if (mynode.eq.nnodes-1)
     &     uau(itr,itr) = ddot(nov,u,1,unxt,1)
        call dcopy(itr,uau(itr,1),maxc,buf,1)
        call pg_dgop(202+itr,buf,itr,'+')
        call dcopy(itr,buf,1,uau(itr,1),maxc)

        do i = 1 , itr - 1
          if (i.ge.jl_u.and.i.le.jh_u) then
            call pg_get(g_u, 1, nov, i, i, buf, 1)
            uau(i,itr) = ddot(nov,buf,1,unxt,1)
          end if
        end do
        call dcopy(itr-1,uau(1,itr),1,buf,1)
        call pg_dgop(303+itr,buf,itr-1,'+')
        call dcopy(itr-1,buf,1,uau(1,itr),1)

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

        do i = 1 , itr
          ccjn = cc(i)
          if (dabs(ccjn).gt.tich) then
            if (i.ge.jl_u.and.i.le.jh_u) then
              call pg_get(g_u, 1, nov, i, i, buf, 1)
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
            call pg_get(g_u, 1, nov, i, i, buf, 1)
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
        if (nxtr.ge.jl_u.and.nxtr.le.jh_u)
     &     call pg_put(g_u, 1, nov, nxtr, nxtr, unxt, 1)
        call dcopy(nov,unxt,1,u,1)
1     continue
c
c  end of loop
c
      write (iwr,*) ' no full convergence after 50 iterations '
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
c
c  noddy routine to do partial multiply on each node
c
      subroutine mxmau_ga(v1,v2,buf,buf2,no,nv)
      implicit REAL  (a-h,o-z)
      integer no,nv,nvt,nvv,i,j,ij,ic,jc,a,b
      dimension v1(*), v2(*), buf(nv,*), buf2(nv,*)
INCLUDE(common/global)
      data four/4.0d0/
      nvt = nv*(nv+1)/2
      nvv = nv*nv

      ij = 0
      do i = 1 , no
        ic = (i-1)*nv + 1
        do j = 1 , i
          jc = (j-1)*nv + 1
          ij = ij + 1
          if (ij.ge.jl_vvoo.and.ij.le.jh_vvoo) then
            call pg_get(g_vvoo, 1, nvt, ij, ij, buf2, 1)
c  hessian generated otf
            call tr2sq_ga(buf2,buf,nv)
            call pg_get(g_vovo, 1, nvv, ij, ij, buf2, nv)
            do a = 1 , nv
              do b = 1 , a
                x = buf2(a,b)
                y = buf2(b,a)
                z = buf(a,b)
                buf(a,b) = four*x - y - z
                buf(b,a) = four*y - x - z
              end do
            end do
c        call wrtsqm_ga('hessian block',nv,buf)
            call mxmb(buf,nv,1,v1(ic),1,nv,v2(jc),1,nv, nv,nv,1)
            if (i.ne.j) 
     &      call mxmb(buf,1,nv,v1(jc),1,nv,v2(ic),1,nv, nv,nv,1)
          end if
        end do
      end do
      return
      end
c
c  construct remaining terms of W and P
c
      subroutine makew_ga(p, z, w, e, t, buf, no,nv,n) 
      implicit REAL  (a-h,o-z)
      dimension p(n,n), w(n,n), z(*), e(*)
      dimension t(no,no),buf(*)
      integer a,b,av
INCLUDE(common/global)
INCLUDE(common/timeperiods)
      data four/4.0d0/, half/0.5d0/
      noo = no*(no+1)/2
      nosq = no*no
      no1 = no + 1

      call start_time_period(TP_MP2MAKEW)

c  change sign of P(2)
      do l = 1 , n
        do m = 1 , n
          p(l,m) = -p(l,m)
        end do
      end do
c  fill in VO-block of P(2)
      ia = 1
      do i = 1 , no
        do a = no1 , n
          p(a,i) = -z(ia)
          ia = ia + 1
        end do
      end do
c  2nd term rhs (14)
      do i = 1 , no
        do a = no1 , n
          w(a,i) = -p(a,i)*e(i)
        end do
      end do
c
c  final contribution to 3rd term rhs (12):
      call vclr(t,1,nosq)
c  vooo
      ijk = 0
      do i = 1 , no
        do j = 1 , i
          do k = 1 , no
            ijk = ijk + 1
            if (ijk.ge.jl_vooo.and.ijk.le.jh_vooo) then
              call pg_get(g_vooo,1,nv,ijk,ijk,buf,1)
              do a = 1 , nv
                av = a + no
                t(i,j) = t(i,j) + buf(a)*p(av,k)*four
                t(i,k) = t(i,k) - buf(a)*p(av,j)
                t(j,k) = t(j,k) - buf(a)*p(av,i)
              end do
              if (i.ne.j) then
                do a = 1 , nv
                  av = a + no
                  t(j,i) = t(j,i) + buf(a)*p(av,k)*four
                  t(k,i) = t(k,i) - buf(a)*p(av,j)
                  t(k,j) = t(k,j) - buf(a)*p(av,i)
                end do
              end if
            end if
          end do
        end do
      end do
c  global sum vooo term
      call pg_dgop(1004,t,nosq,'+')
c  combine terms of (12)
      do i = 1 , no
        do j = 1 , no
          w(i,j) = half*( w(i,j) - p(i,j)*( e(i) + e(j) ) - t(i,j) )
        end do
      end do
c  combine terms of (13)
      do a = no1 , n
        do b = no1 , n
          w(a,b) = half*( w(a,b) - p(a,b)*( e(a) + e(b) ) )
        end do
      end do

      call end_time_period(TP_MP2MAKEW)

      return
      end
c  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MP2 2PDM AND DERIVATIVE INTEGRALS SECTION
c
c  ga version of jkder - DLB task distribution from here on
c
      subroutine jkder_ga(q,schwa)
***   subroutine jkder_ga(q,schwa,iso)
c
c       ----- two electron contribution to the gradient -----
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
      parameter (mxp2 = mxprms*mxprms)
INCLUDE(common/prnprn)
      common/restrl/ociopt,ocifor,omp2,ohf(7),omp2w
INCLUDE(common/cndx41)
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
INCLUDE(common/restar)
INCLUDE(common/restri)
c
INCLUDE(common/specal)
      common/craypk/ijkl(1360)
      common/sortpk/ijkl1(1360)
      common/three/ijkl2(1360)
      common/lsort/ijkl3(1360)
      common/dipblk/dipd(3,3,maxat),dipn(3,3,maxat),dipi(3,3,maxat)
INCLUDE(common/atmblk)
      common/dbuf/icnt,mxtr,val(5118)
      common/dlabs/iao(5,3120)
INCLUDE(common/ghfblk)
INCLUDE(common/cigrad)
      parameter (ncmax=65)
INCLUDE(common/iofile)
INCLUDE(common/infoa)
INCLUDE(common/nshel)
INCLUDE(common/dshlno)
INCLUDE(common/symtry)
INCLUDE(common/timez)
      common/tgrad/dgout(9)
      common/junk/ddij(ncmax,225),ddkl(ncmax,225),aei(ncmax),
     + aej(ncmax),aek(ncmax),ael(ncmax),
     + aa(mxp2),r(mxp2),p1(mxp2),q1(mxp2),r1(mxp2),
     + dd(4*mxp2),ijden(225),ik(225),
     + ijx(225),ijy(225),ijz(225),klx(225),kly(225),klz(225),
     + dij(225),dkl(225),ijgt(225),klgt(225)
INCLUDE(common/ijlab)
INCLUDE(common/dmisc)
INCLUDE(common/segm)
INCLUDE(common/cslosc)
INCLUDE(common/prints)
INCLUDE(common/incrd)
INCLUDE(common/parallel)
      dimension q(*), schwa(*)
c
      dimension dbufs(5120)
      equivalence(dbufs(1),icnt)
c
***   dimension iso(nshell,*)
      character *14 fnm
      character *8 snm
      data fnm/'mp2_parallel.m'/
      data snm/'jkder_ga'/
      data zrhf,zuhf,zgvb,zgrhf/
     * 'rhf','uhf','gvb','grhf' /
INCLUDE(common/global)
INCLUDE(common/timeperiods)
INCLUDE(common/mp2grad_pointers)
INCLUDE(common/mapper)
      logical out, outv
***   logical sym_pair, sym_quartet_l
      common/shlt/tol,cutoff,icount, ic999, out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
c
      call start_time_period(TP_JKDER)
c
***   call priso(iso,nshell,nt)
c
      call setsto(1360,0,ijkl)
      if (nprint.ne.-5 .and. oprint(57)) then
         write (iwr,6010)
         call writel(q(iprefa),nshell)
      end if
      dlncutoff = dlog(cutoff)

c      call wrtsqm_ga('jkder P(2)   ',ncoorb,q(mp2grad_pmat))
c      call wrtsqm_ga('jkder P(HF)  ',ncoorb,q(mp2grad_dens))
c      call wrtsqm_ga('jkder evecs  ',ncoorb,q(mp2grad_vecs))

c
c     ----- set some parameters -----
c
      no = nocca
      nv = nvirta
      nsq = ncoorb*ncoorb
      ntri = ncoorb*(ncoorb+1)/2
      nov = nocca*nvirta
      noo = no*(no+1)/2
      nosq = no*no
      nn = nshell*(nshell+1)/2
c     mxshl=0
      mxshl=4
      do i=1,nshell
c        mxshl = max(mxshl,kmax(i)-kmin(i)+1)
         if(ktype(i).eq.3) then
          mxshl = max(mxshl,6)
         else if(ktype(i).eq.4) then
          mxshl = max(mxshl,10)
         else if(ktype(i).eq.5) then
          mxshl = max(mxshl,15)
         else
         endif
      enddo
      nnodes = ipg_nnodes()
      nc = ncoorb
      nschwz = 0
      ndgout = 9
      nav = lenwrd()
      ntpdm = 1
      oeof = .false.
      l1 = num
      l2 = (num*(num+1))/2
      l3 = num*num
      nat3 = nat*3
      nbsq = num*num
      lenb = lensec(nx)
      odpres = .false.
      ofpres = .false.
      ogpres = .false.
      kt_max=0
      do 20 i = 1 , nshell
         if (ktype(i).eq.3) odpres = .true.
         if (ktype(i).eq.4) ofpres = .true.
         if (ktype(i).eq.5) ogpres = .true.
         kt_max=max(kt_max,ktype(i))
 20   continue
c
c ps - limits for dgenrl buffer in static allocation version
c
      inc1_max=(kt_max+1)**4
      ncmmm_max = 32
c
      m = ntpdm + 9

      lnddm = 257
      if (odpres)lnddm = 1297
      if (ofpres)lnddm = 10001
      if (ogpres)lnddm = 50625
c
      ifok=99999999

c
c??? any chance of ndens != 1???
c
c
c remnants of original memory area
c
      i10 = 0
      i20 = i10 + maxorb
      i30 = i20 + l2
      len1 = i30 + l2
      len2 = i20 + ndens*l2
      i10 = igmem_alloc(max(len1,len2))

      i20 = i10 + maxorb
      ida = i20 - 1
      ndens = 1
      i30 = i20 + l2
c
c  allocate separate sections for buffers
c    - nb watch for offsets on ic7, iabd
c
      ic7 = igmem_alloc(lnddm*9 + inc1_max*ncmmm_max*6)
      ic7 = ic7 - 1

      iabd = igmem_alloc(lnddm*ntpdm)
      iabd = iabd - 1
cc
cc  separate allocation of ic1 reverted, because of
cc  problems with mp2 second deriv. Size is included
cc  in ic7
cc
cc      ic1 = igmem_alloc(inc1_max*ncmmm_max*6)

      i00 = igmem_alloc(lnddm*(3/nav+1))

      tim0 = cpulft(1)
c
c  pass 2*P(HF) to ddebut
c
      call sq2tr_ga( q(mp2grad_dens), q(i20), ncoorb)
      call vadd(q(i20),1,q(i20),1,q(i20),1,l2)

      zscftp = 'rhf'
      call ddebut(zscftp,q(i20),q(i30),q(i10),q(i10),q(ifok+1),
     +               l1,l2,l3,ist,jst,kst,lst)
c      call search(iblk2d,ifil2d)
      mword = 0
      iword = 0
      kloc(nshell+1) = num + 1
      icnt = 0
      ib1 = 1
c
c   ----  first half-back-transformation of amplitudes  ----
c
      i2 = igmem_alloc(nsq)
      i3 = igmem_alloc(nsq)
      call backtran_ga(q(mp2grad_vals), ncoorb,nocca,nvirta
     &,                q(mp2grad_vecs),  q(i2),q(i3))
      call gmem_free(i3)
      call gmem_free(i2)
c
c   -----  two-electron derivative integrals ----
c
c  workspace for 2nd back-transformation
c
      igns = igmem_alloc(mxshl*mxshl*nosq)
      nonsep = igmem_alloc(mxshl*mxshl*mxshl*ncoorb)
      idum = igmem_alloc(max(no*mxshl,noo))
c
c  **** MPP
      call pg_dlbreset
      next = ipg_dlbtask()
c  **** MPP
c
c     ----- ishell -----
c
      do 140 ii = 1 , nshell
        nmini = kloc(ii)
        nmaxi = kloc(ii+1) - 1
        ishl = nmaxi - nmini + 1
        iceni = katom(ii)
c
c     ----- kshell -----
c
        do 130 kk = 1 , ii
          nmink = kloc(kk)
          nmaxk = kloc(kk+1) - 1
          kshl = nmaxk - nmink + 1
          olab = katom(kk).eq.iceni
c
c  symmetry (IK|
c
***       if (sym_pair(iso,nshell,nt,ii,kk)) then
c
c***   **MPP** DLB counter
c
          icount_dlb = icount_dlb + 1
          if(icount_dlb . eq. next) then
c
c  get all A(nu,i|lamda,j) for nu,lamda spanning (IK)th shell
c
            call start_time_period(TP_JKDER_GET)

            i9oo = igns
            do nu = nmini , nmaxi
              do lamda = nmink , nmaxk
                call get_half_trans(q(i9oo),q(idum),ncoorb
     &,                             nocca,nu,lamda)
                i9oo = i9oo + nosq
               end do 
           end do

            call end_time_period(TP_JKDER_GET)
c
c     ----- jshell -----
c
            do 120 jj = 1 , nshell
              nminj = kloc(jj)
              nmaxj = kloc(jj+1) - 1
              jshl = nmaxj - nminj + 1
              olabc = olab .and. katom(jj).eq.iceni
c
c  store information about the pair (ij)
c
                call dshell(1,ii,jj,kk,ll)
                call dprim
                if (nij.ne.0) then
c
c  second half back-transformation of amplitudes
c
                  call start_time_period(TP_MP2BACKTRAN_2)

                  i9oo = igns
                  i7ns = nonsep
                  do nu = nmini , nmaxi
                    do lamda = nmink , nmaxk
                      call tran2_ga(ncoorb,nocca,nvirta
     &,                             q(mp2grad_vecs), q(i9oo)
     &,                             q(i7ns), nminj,jshl,q(idum))
                      i9oo = i9oo + nosq
                      i7ns = i7ns + ncoorb*jshl
                    end do
                  end do

                  call end_time_period(TP_MP2BACKTRAN_2)
c
c     ----- lshell -----
c
                  do 110 ll = 1 , nshell
                    nminl = kloc(ll)
                    nmaxl = kloc(ll+1) - 1
                    lshl = nmaxl - nminl + 1
                    lshl = nmaxl - nminl + 1
c  ... schwarz inequality test
                    itrij = iky(max(ii,jj)) + min(ii,jj)
                    itrkl = iky(max(kk,ll)) + min(kk,ll)
                    test = schwa(itrij) + schwa(itrkl)
                    if (test.lt.dlncutoff) then
                       nschwz = nschwz + 1
                    else
                      olabcd = olabc .and. katom(ll).eq.iceni
                      if (.not.(olabcd)) then
c
c  symmetry (IJ|KL)
c
***                   if (sym_quartet_l(iso,nshell,nt,
***  &                    ii, jj, kk, ll, q4)) then
c
c     ----- check for redundant combInations -----
c
                        call redund_ga(ii,jj,kk,ll,iwr)
                        if (npass.ne.0) then
                          call vclr(dgout,1,ndgout)
c
c     ----- store information about the pair (kl)
c
                          call dshell(2,ii,jj,kk,ll)
                          call vclr(q(iabd+1),1,lendd*ntpdm)
c
c     ----- combine terms of the 2-particle density matrix -----
c
c                          call start_time_period(TP_MP2MCDAB)

                          call mcdab_ga(q(iabd+1),ii,jj,kk,ll,
     &                                  q(mp2grad_pmat),
     &                                  q(mp2grad_dens),
     &                    nocca, nvirta, ncoorb, q(nonsep))

c                          call end_time_period(TP_MP2MCDAB)
c
c  ---- mess about with the density matrix by eliminating zero elements
c
                          call delim(q(iabd+1),q(ic7+1),
     +                               ijgt,klgt,q(i00),ijd,kld,
     +                               ijkld,lendd,abmax)
c
c  ---- calc derivative integrals
c
                          if (ijkld.ne.0) then
c                            call start_time_period(TP_DGENRL)

                             call dgenrl(q(1),q(i00),q(i00),abmax)

c                            call end_time_period(TP_DGENRL)
c     -----  combine all 4 partial contributions to the gradient ----
                            call formeg
                          end if
                        end if
                      end if
                    end if
***               end if
110             continue
              end if
120         continue
c***   **MPP**
            next = ipg_dlbtask()
          end if
c***   **MPP**
***      end if
130     continue
c
c     ----- save gradient and restart data -----
c
        call dfinal(q,0,ii)
        if (tim.ge.timlim) go to 155
140   continue
c
155   continue
_IF(parallel)
      call pg_dlbpush
_ENDIF
c
c   ----- release memory ------ 
c
      call gmem_free(idum)
      call gmem_free(nonsep)
      call gmem_free(igns)
      call gmem_free(i00)
cc reverted
cc      call gmem_free(ic1)
      iabd = iabd+1
      call gmem_free(iabd)
      ic7 = ic7+1
      call gmem_free(ic7)
      call gmem_free(i10)
c
      if (tim.ge.timlim) go to 150
c
c   ----- end of *shell* loops -----
c
      if (nt.ne.1) then
c
c   ----- symmetrise gradient vector
c
         isymd = igmem_alloc(nw196(6))
         call symde(q(isymd),nat)
         call gmem_free(isymd)

      endif
c
      call dfinal(q,1,ii)
      nindmx = 0
 150  call timit(0)
c  ... schwarz inequality log
      if(nprint.ne.-5) write(iwr,6040) nschwz
      dtim = tim - tim0
c  delete vovo GA
      call delete_ga_inf(g_vovo,'vovo',fnm,snm)
      call end_time_period(TP_JKDER)
      return
6010  format (/40x,'prefactor matrix in 2e-derivatives'/)
6020  format (/1x,'derivative integrals to be output')
6030  format (/1x,'derivative integrals written')
6040  format(/1x,'schwarz inequality test skipped ',i10,
     + ' integral blocks')
      end
c
c  first half-back-transformation of amplitudes
c  in this version, the vovo GA is overwritten
c
      subroutine backtran_ga(e, n,no,nv, vec, buf, dum)
      implicit REAL  (a-h,o-z)
      dimension e(*), vec(n,*), buf(nv,*), dum(*)
INCLUDE(common/global)
INCLUDE(common/timeperiods)
      data two/2.0d0/

      call start_time_period(TP_MP2BACK_F)

      no1 = no + 1
      nsq = n*n
      nvv = nv*nv
c  loop over local blocks
      ij = 0
      do i = 1 , no
        do j = 1 , i
          ij = ij + 1
          if (ij.ge.jl_vovo.and.ij.le.jh_vovo) then
            call pg_get(g_vovo, 1, nvv, ij, ij, buf, nv)
c  form block of double substitution amplitudes
            do k = 1 , nv
              do l = 1 , k
                x = buf(k,l)
                y = buf(l,k)
                denom = two/(e(i)+e(j)-e(k+no)-e(l+no))
                buf(k,l) = (two*x - y)*denom
                buf(l,k) = (two*y - x)*denom
              end do
            end do
c  back-transformation
            call vclr(dum,1,nsq)
            call mxmb(buf,1,nv,vec(1,no1),n,1,dum,1,nv, nv,nv,n)
            call vclr(buf,1,nsq)
            call mxmb(vec(1,no1),1,n,dum,1,nv,buf,1,n,   n,nv,n)
            call pg_put(g_vovo, 1, nsq, ij, ij, buf, n)
          end if
        end do
      end do

      call end_time_period(TP_MP2BACK_F)
      return
      end
c
c  get a block of half-transformed amplitudes
c
      subroutine get_half_trans(aoao,buf,n,no,nu,lamda) 
      integer n,no,nu,lamda
      REAL aoao,buf
      dimension aoao(no,no),buf(*)
INCLUDE(common/global)
      noo = no*(no+1)/2
      ik = (nu - 1)*n + lamda
      call pg_get(g_vovo, ik, ik, 1, noo, buf, 1)
      call tr2sq_ga(buf,aoao,no)
      if (nu.ne.lamda) then
        ki = (lamda - 1)*n + nu
        call pg_get(g_vovo, ki, ki, 1, noo, buf, 1)
        ij = 0
        do i = 1 , no 
          do j = 1 , i
            ij = ij + 1
            aoao(j,i) = buf(ij)
          end do
        end do
      end if
      return
      end
c
c  second half back-transformation of amplitudes
c
      subroutine tran2_ga(n,no,nv, vec, oaoa, gns
     &,  minj,jshl, dum)
      implicit REAL  (a-h,o-z)
      integer n,no,nv,minj,jshl
      dimension vec(n,*), oaoa(no,no), gns(*), dum(*)
INCLUDE(common/global)
INCLUDE(common/timeperiods)

      no1 = no + 1
c  back-transform occupied indices to shell AO indices
      call vclr(dum,1,no*jshl)
      call mxmb(oaoa,1,no,vec(minj,1),n,1,dum,1,no,  no,no,jshl)
c  back-transform occupied indices to all AO indices
      call vclr(gns,1,n*jshl)
      call mxmb(vec,1,n,dum,1,no,gns,1,n,  n,no,jshl)
      return
      end
_EXTRACT(mcdab_ga,linux)
c
c  macdab using partial back-transformed amplitudes
c
      subroutine mcdab_ga(abdens,ii,jj,kk,ll
     &,  pmp2, hfden, nocca, nvirta, ncoorb, gns)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/nshel)
INCLUDE(common/cigrad)
INCLUDE(common/incrd)
INCLUDE(common/dmisc)
      dimension abdens(*), pmp2(*), hfden(*), gns(*)
      data half,one,two/0.5d0,1.0d0,2.0d0/
      logical iieqjj, iieqkk, iieqll, jjeqkk, jjeqll, kkeqll 
c
      iieqjj = ii.eq.jj
      iieqkk = ii.eq.kk
      iieqll = ii.eq.ll
      jjeqkk = jj.eq.kk
      jjeqll = jj.eq.ll
      kkeqll = kk.eq.ll
      mini = kloc(ii)
      minj = kloc(jj)
      mink = kloc(kk)
      minl = kloc(ll)
      maxi = kloc(ii+1) - 1
      maxj = kloc(jj+1) - 1
      maxk = kloc(kk+1) - 1
      maxl = kloc(ll+1) - 1
      ishl = maxi - mini + 1
      jshl = maxj - minj + 1
      kshl = maxk - mink + 1
      lshl = maxl - minl + 1

      fac = one
      if ( iieqkk .and. jj.ne.ll )  fac = half
c
c  form 2pdm over primitives
c
      do i = mini , maxi
         do j = minj , maxj
            j11 = (j-1)*ncoorb
            imn = j11 + i
            do k = mink , maxk
               k11 = (k-1)*ncoorb
               iml = k11 + i
               inl = k11 + j
               do l = minl , maxl
                  l11 = (l-1)*ncoorb
                  ils = l11 + k
                  ims = l11 + i
                  ins = l11 + j
c  SCF term of 2pdm -           ij         kl
                  dscf = hfden(imn)*hfden(ils)*4.0d0
     +                      - hfden(iml)*hfden(ins)
     +                      - hfden(ims)*hfden(inl)
c  MP2 separable term  - 
                  dsep =  pmp2(imn)
     +                   *hfden(ils) + pmp2(ils)*hfden(imn)
     +                   - (pmp2(iml)*hfden(ins)+pmp2(ims)
     +                   *hfden(inl)+pmp2(inl)*hfden(ims)
     +                   +pmp2(ins)*hfden(iml))*0.25d0
c  MP2 non-separable term -
                  igns = (i-mini)*kshl*jshl*ncoorb + 
     &                   (j-minj)*ncoorb + 
     &                   (k-mink)*jshl*ncoorb + l
                  dnon = gns(igns)
c  total MP2 gradient density -
                  gamma = dscf + two*(dsep + dnon)
c  permutational symmetry weight - 
                  gamma = gamma*fac
c gdf:                  gamma = gamma*q4
                  nn = (i-mini)*inc5 + (j-minj)*inc4 + 
     &                 (k-mink)*inc3 + (l-minl) + inc2
                  abdens(nn) = abdens(nn) + gamma
               end do
            end do
         end do
      end do
c
c  symmetrize density matrix 
c
      if ( iieqjj .or. kkeqll ) then 
      do i = mini , maxi
        jjmx = maxj
        if ( iieqjj ) jjmx = i
        do j = minj , jjmx
          do k = mink , maxk
            llmx = maxl
            if ( kkeqll ) llmx = k
            do l = minl , llmx
              ni = (i-mini)*inc5 + (j-minj)*inc4 + 
     &             (k-mink)*inc3 + (l-minl) + inc2
              if ( iieqjj .and. i.ne.j ) then
                nj = (j-minj)*inc5 + (i-mini)*inc4 + 
     &               (k-mink)*inc3 + (l-minl) + inc2
                di = abdens(ni)
                dj = abdens(nj)
                abdens(ni) = di + dj
                abdens(nj) = abdens(ni) 
                if ( kkeqll .and. k.ne.l ) then
                  nk = (i-mini)*inc5 + (j-minj)*inc4 + 
     &                 (l-minl)*inc3 + (k-mink) + inc2
                  nl = (j-minj)*inc5 + (i-mini)*inc4 + 
     &                 (l-minl)*inc3 + (k-mink) + inc2
                  dk = abdens(nk)
                  dl = abdens(nl)
                  abdens(nk) = dk + dl  
                  abdens(nl) = abdens(nk) 
                  abdens(ni) = half*( abdens(ni) + 
     &               abdens(nj) + abdens(nk) + abdens(nl) )
                  abdens(nj) = abdens(ni)
                  abdens(nk) = abdens(ni)
                  abdens(nl) = abdens(ni)
                end if
              else if ( kkeqll .and. k.ne.l ) then
                nj = (i-mini)*inc5 + (j-minj)*inc4 + 
     &               (l-minl)*inc3 + (k-mink) + inc2
                di = abdens(ni)
                dj = abdens(nj)
                abdens(ni) = di + dj 
                abdens(nj) = abdens(ni) 
              end if
            end do
          end do
        end do
      end do
c
      else if ( ( iieqkk .and. jjeqll ) .or. 
     &          ( iieqll .and. jjeqkk ) ) then 
      do i = mini , maxi
        do j = minj , maxj
          do k = mink , maxk
            do l = minl , maxl
              ni = (i-mini)*inc5 + (j-minj)*inc4 + 
     &             (k-mink)*inc3 + (l-minl) + inc2
              if ( iieqkk .and. jjeqll ) then 
                nj = (k-mink)*inc5 + (l-minl)*inc4 + 
     &               (i-mini)*inc3 + (j-minj) + inc2
              else if ( iieqll .and. jjeqkk ) then 
                nj = (l-minl)*inc5 + (k-mink)*inc4 + 
     &               (j-minj)*inc3 + (i-mini) + inc2
              end if
              di = abdens(ni)
              dj = abdens(nj)
              abdens(ni) =  ( di + dj )*half 
              abdens(nj) = abdens(ni) 
            end do
          end do
        end do
      end do
        if ( ii.ne.kk ) 
     &   call dscal(lendd,two,abdens,1)
      end if
      return
      end
_ENDEXTRACT
c
c gdf:   25.1.95 - modifications to REDUND (indicated below)
c   accounting for different permutational symmetries of the 
c   two-electron integrals in the parallel MP2 gradient 
c   calculation 
c
      subroutine redund_ga(ii,jj,kk,ll,iw)
      implicit REAL  (a-h,o-z)
      logical inej,inek,inel,jnek,jnel,knel
INCLUDE(common/sizes)
INCLUDE(common/dmisc)
INCLUDE(common/nshel)
c
      dimension lll(4)
      equivalence (lll(1),lit)
      oskip(1) = .true.
      oskip(2) = .true.
      oskip(3) = .true.
      oskip(4) = .true.
      npass = 0
      do 20 i = 1 , 4
         natomd(i) = 0
 20   continue
      lit = ktype(ii)
      ljt = ktype(jj)
      lkt = ktype(kk)
      llt = ktype(ll)
      iat = katom(ii)
      jat = katom(jj)
      kat = katom(kk)
      lat = katom(ll)
      inej = iat.ne.jat
      inek = iat.ne.kat
      inel = iat.ne.lat
      jnek = jat.ne.kat
      jnel = jat.ne.lat
      knel = kat.ne.lat
      if (inej) then
         if (.not.(inek)) then
            if (.not.(inel)) go to 40
            if (jnel) go to 50
c      iat=kat    jat=lat
            if (ii.ne.kk .or. jj.ne.ll) then
               n1 = (lit+1)*(lkt+1)*ljt*llt
               n2 = lit*lkt*(ljt+1)*(llt+1)
               if (n1.ge.n2) go to 50
               go to 70
            else
c gdf:   30.1.95 - start
c               if (ljt.le.lit) go to 60
               if (ljt.lt.lit) go to 60
               if (ljt.eq.lit.and.jat.lt.iat) goto 60
c gdf:   30.1.95 - end
               go to 40
            end if
         else if (jnek) then
            if (.not.(jnel)) go to 70
            if (.not.(knel)) go to 80
c gdf:   25.1.95 - start
            if (.not.(inel)) go to 85
c gdf:   25.1.95 - end 
c     iat # jat # kat # lat  -- omit one centre
            mmin = lit
            imin = 1
            do 30 iper = 2 , 4
               if (lll(iper).lt.mmin) then
                  mmin = lll(iper)
                  imin = iper
               end if
 30         continue
            go to (90,100,110,120) , imin
            go to 90
         else
            if (.not.(jnel)) go to 60
c gdf:   25.1.95 - start
c     iat=lat  jat=kat   derivative (i'j/kl')
            if ((.not.jnek).and.(.not.inel)) then 
               if (ii.ne.ll .or. jj.ne.kk) then
                  n1 = (lit+1)*(llt+1)*ljt*lkt
                  n2 = lit*llt*(ljt+1)*(lkt+1)
                  if (n1.ge.n2) go to 85
                  go to 35
               else
               if (ljt.le.lit) go to 60
                  go to 40
               end if
            end if
c gdf:   25.1.95 - end 
c     ----- jat = kat # iat # lat -----
35             oskip(1) = .false.
               oskip(4) = .false.
               natomd(1) = iat
               natomd(2) = lat
               natomd(3) = jat
               npass = 2
               go to 130
         end if
      else if (inek) then
         if (.not.(knel)) then
c     iat=jat  ,  kat=lat   differentiate one pair
            n1 = lit*ljt*(lkt+1)*(llt+1)
            n2 = (lit+1)*(ljt+1)*lkt*llt
            if (n2.lt.n1) go to 80
c gdf:   25.1.95 - start
         else if (knel.and.(.not.inel)) then
c     iat=jat=lat   derivative (ij/kl)
            oskip(3) = .false.
            natomd(1) = kat
            natomd(2) = iat
            npass = 1
            go to 130
c gdf:   25.1.95 - end 
         end if
      else
         if (inel) then
c     iat=jat=kat   derivative (ij/kl)
            oskip(4) = .false.
            natomd(1) = lat
            natomd(2) = iat
            npass = 1
         end if
         go to 130
      end if
c     iat=jat   derivatives (ij/k'l) and (ij/kl')
      oskip(3) = .false.
      oskip(4) = .false.
      natomd(1) = kat
      natomd(2) = lat
      natomd(3) = iat
      npass = 2
      go to 130
c     iat=kat=lat   derivative (ij/kl)
 40   oskip(2) = .false.
      natomd(1) = jat
      natomd(2) = iat
      npass = 1
      go to 130
c     iat=kat   derivatives (ij'/kl) and (ij/kl')
 50   oskip(2) = .false.
      oskip(4) = .false.
      natomd(1) = jat
      natomd(2) = lat
      natomd(3) = iat
      npass = 2
      go to 130
c     jat=kat=lat    (ij/kl)
 60   oskip(1) = .false.
      natomd(1) = iat
      natomd(2) = jat
      npass = 1
      go to 130
c      jat=lat    derivatives (i'j/kl) and (ij/k'l)
 70   oskip(1) = .false.
      oskip(3) = .false.
      natomd(1) = iat
      natomd(2) = kat
      natomd(3) = jat
      npass = 2
      go to 130
c     kat=lat   derivatives (i'j/kl) and (ij'/kl)
 80   oskip(1) = .false.
      oskip(2) = .false.
      natomd(1) = iat
      natomd(2) = jat
      natomd(3) = kat
      npass = 2
      go to 130
c gdf:   25.1.95 - start
c     iat = lat # jat # kat -----
 85   oskip(2) = .false.
      oskip(3) = .false.
      natomd(1) = jat
      natomd(2) = kat
      natomd(3) = iat
      npass = 2
      go to 130
c gdf:   25.1.95 - end 
 90   natomd(1) = jat
      natomd(2) = kat
      natomd(3) = lat
      natomd(4) = iat
      npass = 3
      oskip(2) = .false.
      oskip(3) = .false.
      oskip(4) = .false.
      go to 130
 100  natomd(1) = iat
      natomd(2) = kat
      natomd(3) = lat
      natomd(4) = jat
      npass = 3
      oskip(1) = .false.
      oskip(3) = .false.
      oskip(4) = .false.
      go to 130
 110  natomd(1) = iat
      natomd(2) = jat
      natomd(3) = lat
      natomd(4) = kat
      npass = 3
      oskip(1) = .false.
      oskip(2) = .false.
      oskip(4) = .false.
      go to 130
 120  oskip(1) = .false.
      oskip(2) = .false.
      oskip(3) = .false.
      natomd(1) = iat
      natomd(2) = jat
      natomd(3) = kat
      natomd(4) = lat
      npass = 3
c     -----
 130  if (.not.outd) return
      write (iw,6010) ii , jj , kk , ll , 
     +    oskip(1) , oskip(2) , oskip(3) , oskip(4) , 
     +    npass , (natomd(i),i=1,4)
      return
 6010 format (/,' ********  ii,jj,kk,ll =',4i3,' skip1,2,3,4 =',4l3,
     +        ' npass =',i2,' centers =',4i5,/)
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   UTILITIES SECTION
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  routine to delete a GA
c
      subroutine delete_ga(handle)
      integer handle
      logical pg_destroy, opg_root
INCLUDE(common/parcntl)
      if (pg_destroy(handle)) then
        if (opg_root().and.iparapr.eq.3) write(*,10) handle
      else
        call pg_error('**GA error** failed to destroy GA ',handle)
      end if
      return
 10   format(/1x,
     +   'GA: pg_destroy - double precision handle [',i10,']')
      end 
      subroutine delete_ga_inf(handle,anm,fnm,snm)
      integer handle
      logical pg_destroy_inf, opg_root
      character *(*) anm,fnm,snm
INCLUDE(common/parcntl)
      if (pg_destroy_inf(handle,anm,fnm,snm)) then
        if (opg_root().and.iparapr.eq.3) write(*,10) handle
      else
        call pg_error('**GA error** failed to destroy GA ',handle)
      end if
      return
 10   format(/1x,
     +   'GA: pg_destroy_inf - double precision handle [',i10,']')
      end 
c
c  extra layer for read 
c
      subroutine get_mo_coeffs(coeffs,len)
INCLUDE(common/restar)
INCLUDE(common/restri)
INCLUDE(common/atmblk)
      integer len,iblok
      REAL coeffs(*)
      call secget(isect(8),8,iblok)
      iblok = iblok + mvadd
      call rdedx(coeffs,len,iblok,ifild)
      end
c
c  zero a GA
c
      subroutine zero_ga(handle,ilo,ihi,jlo,jhi)
      integer handle,ilo,ihi,jlo,jhi,i
      do i = ilo , ihi
         call pg_put(handle,i,i,jlo,jhi,0.0d0,0)
      end do
      end
c
c  write a vector, comments and node number for parallel debugging
c
      subroutine wrtvec_ga(string,n,a)
      integer n,j,mynode
      REAL a(*)
      character*(*) string
      mynode = ipg_nodeid()
      write(*,1) string, mynode
      write(*,999)(a(j),j=1,n)
      return
1     format(/,5x,'[ ',a20,' ]  node = ',1i3)
999   format(12f14.8)
      end
c
c  write a square matrix, comments and node number for parallel debugging
c
      subroutine wrtsqm_ga(string,n,a)
      implicit REAL  (a-h,o-z)
INCLUDE(common/iofile)
      dimension a(n,n)
      character*(*) string
      mynode = ipg_nodeid()
      write(iwr,1) string, mynode
      do i = 1 , n 
        write(iwr,999)(a(i,j),j=1,n)
      end do
      return
1     format(/,5x,'[ ',a20,' ]  node = ',1i3)
999   format(6f12.8)
      end
c
c  write a GA, comments and node number for parallel debugging
c
      subroutine wrtmat_ga(string,handle,dum)
      implicit REAL  (a-h,o-z)
INCLUDE(common/iofile)
      integer handle, ilo, ihi, jlo, jhi, n, i,j,mynode
      dimension dum(100000)
      character*(*) string
      mynode = ipg_nodeid()
      write(iwr,1) string, mynode
      call pg_distribution(handle, mynode, ilo, ihi, jlo, jhi)
      n = jhi - jlo + 1
      sum2 = 0.0d0
      do i = ilo , ihi
        call pg_get(handle, i, i, jlo, jhi, dum, 1)
*       write(iwr,999)(dum(j),j=1,n)
      sum2 = sum2 + dsum(n,dum,1)
      end do
      call pg_dgop(999,sum2,1,'+')
      write(iwr,*)string,sum2
      return
1     format(/,5x,'[ ',a20,' ]  node = ',1i3)
999   format(6f12.8)
      end
c
c  copy symmetric square matrix to triangle
c
      subroutine sq2tr_ga(sq, tr, n)
      implicit none
      integer n,i,j,k
      REAL sq, tr
      dimension sq(n,n), tr(*)
      k = 0 
      do i = 1 , n
        do j = 1 , i
          k = k + 1  
          tr(k) = sq(i,j)
        end do
      end do
      return
      end
c
c  copy triangle to square
c
      subroutine tr2sq_ga(tr, sq, n)
      implicit none
      integer n,i,j,k
      REAL sq, tr
      dimension sq(n,n), tr(*)
      k = 0
      do i = 1 , n
        do j = 1 , i
          k = k + 1
          sq(i,j) = tr(k)
          sq(j,i) = tr(k)
        end do
      end do
      return
      end
c
c   debug write iso array routine
c
      subroutine priso(iso,nshels,nt)
      implicit REAL  (a-h,o-z)
      dimension iso(nshels,*)
      write(*,2)
      do i = 1 , nshels
       write(*,1)(iso(i,j),j=1,nt)
      end do
      write(*,2)
1     format(1i4,')  ',20i4)
2     format(//,
     &' ****************** ISO ARRAY ******************',//)
      return
      end
c  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  BLOCK DATA INITIALISATIONS
c
c  initialize mp2 gradient core pointers 
c
      block data mp2grad_pointers_init
      implicit none
INCLUDE(common/mp2grad_pointers)
      data mp2grad_wmat /0/
      data mp2grad_pmat /0/
      data mp2grad_vecs /0/
      data mp2grad_dens /0/
      data mp2grad_vals /0/
      data mp2grad_schw /0/
      end
c
c  initialize global array pointers
c
      block data global_init
      implicit none
INCLUDE(common/global)
      data g_oooo, g_vooo, g_vvoo, g_vovo, g_vvvo, g_ovaa, g_u
     &,  il_oooo, ih_oooo, jl_oooo, jh_oooo
     &,  il_vooo, ih_vooo, jl_vooo, jh_vooo
     &,  il_vvoo, ih_vvoo, jl_vvoo, jh_vvoo
     &,  il_vovo, ih_vovo, jl_vovo, jh_vovo
     &,  il_vvvo, ih_vvvo, jl_vvvo, jh_vvvo
     &,  il_ovaa, ih_ovaa, jl_ovaa, jh_ovaa
     &,  il_u, ih_u, jl_u, jh_u
     &  /0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
     &,  0,0,0,0,0,0,0,0,0,0,0,0,0,0/
      end
c
c  initialise gradient search switches
c
      block data mp2grd_init
INCLUDE(common/mp2grd)
      data mp2gcycl/-99/
      data opg_grad/.false./
      data opg_alloc/.false./
      data opg_alloc_sch / .false./ 
      data opg_alloc_dmat / .false./
      data opg_grad_com / .true./
c
      end 
c  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  TEST SECTION
c
_IF(aprtest)
      subroutine aprint(g)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/ijlab)
INCLUDE(common/cslosc)
INCLUDE(common/mapper)
INCLUDE(common/shlt)
INCLUDE(common/restar)
INCLUDE(common/shlnos)
INCLUDE(common/misc)
INCLUDE(common/indez)
      dimension g(*)
      REAL aprscr,aprmo
      integer aprbas,lapr
      parameter (lapr=1)
      common/apr_wrk/aprscr(lapr),aprmo(100,100),aprbas
c
c      print *,'putting data into aprint '
      ijn = 0
      jmax = maxj
      do 260 i = mini,maxi
      if (oianj) jmax = i
      int1 = loci + i
      do 240 j = minj,jmax
      ijn = ijn+1
      n1 = ijgt(ijn)
      int2 = locj + j
      i1 = max(int1,int2)
      i2 = min(int1,int2)
      lmax = maxl
      kln = 0
      do 220 k = mink,maxk
      if (okanl) lmax = k
      int3 = lock + k
      do 200 l = minl,lmax
      kln = kln+1
      nn = n1+klgt(kln)
      val = g(nn)
        int4 = locl + l
        i3 = max(int3,int4)
        i4 = min(int3,int4)
        itr12=iky(i1)+i2
        itr34=iky(i3)+i4
        i1234 = max(itr12,itr34) * (max(itr12,itr34)-1)/2 +
     &          min(itr12,itr34)
        if (i1234.le.lapr)then
         aprscr(i1234)=val
        else
         print *,'exceed dimension of aprscr ',i1234,lapr
         stop 20
        endif
c        write(*,435)i1,i2,i3,i4,val
c  435 format(' int ',4i5,e15.5)
  200 continue
  220 continue
  240 continue
  260 continue
      return
      end
      subroutine aprpmo(nbas,cmo)
      implicit none
      integer nbas,i,j
      REAL cmo(nbas,nbas)
      REAL aprscr,aprmo
      integer aprbas,lapr
      parameter (lapr=1)
      common/apr_wrk/aprscr(lapr),aprmo(100,100),aprbas

      aprbas=nbas
      do i=1,nbas
      do j=1,nbas
c      cmo(i,j)=0.0d00
c      if (i.eq.j)cmo(i,j)=1.0d00
      aprmo(i,j)=cmo(i,j)
      enddo
      enddo
      return
      end
      double precision function aprmoint(ii,jj,kk,ll)
      implicit none
      integer ii,jj,kk,ll,i,j,k,l
      integer ij,kl,ijkl
      REAL xxxx

      REAL aprscr,aprmo
      integer aprbas,lapr
      parameter (lapr=1)
      common/apr_wrk/aprscr(lapr),aprmo(100,100),aprbas
      integer mxidx
      mxidx(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)

c      print *,'ii,jj,kk,ll ',ii,jj,kk,ll
      xxxx=0.0d00
      do i=1,aprbas
      do j=1,aprbas
      ij=mxidx(i,j)
      do k=1,aprbas
      do l=1,aprbas
      kl=mxidx(k,l)
      ijkl=mxidx(ij,kl)
c      print *,i,j,k,l,aprscr(ijkl),aprmo(i,ii),aprmo(j,jj),
c     &          aprmo(k,kk),aprmo(l,ll)
      xxxx=xxxx+aprscr(ijkl)*aprmo(i,ii)*aprmo(j,jj)*
     &          aprmo(k,kk)*aprmo(l,ll)
      enddo
      enddo
      enddo
      enddo
      aprmoint=xxxx
c
      return
      end
      double precision function aprq1int(ii,jj,kk,llm)
      implicit none
      integer ii,jj,kk,llm,i,j,k,l
      integer ij,kl,ijkl
      REAL xxxx

      REAL aprscr,aprmo
      integer aprbas,lapr
      parameter (lapr=1)
      common/apr_wrk/aprscr(lapr),aprmo(100,100),aprbas
      integer mxidx
      mxidx(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)

c      print *,'ii,jj,kk,ll ',ii,jj,kk,ll
      xxxx=0.0d00
      ij=mxidx(ii,jj)
      do l=1,aprbas
      kl=mxidx(kk,l)
      ijkl=mxidx(ij,kl)
      xxxx=xxxx+aprscr(ijkl)*aprmo(l,llm)
      enddo
      aprq1int=xxxx
c
      return
      end
      double precision function apraoint(ii,jj,kk,ll)
      implicit none
      integer ii,jj,kk,ll,i,j
      integer ij,kl,ijkl

      REAL aprscr,aprmo
      integer aprbas,lapr
      parameter (lapr=1)
      common/apr_wrk/aprscr(lapr),aprmo(100,100),aprbas
      integer mxidx
      mxidx(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
c
      ij=mxidx(ii,jj)
      kl=mxidx(kk,ll)
      ijkl=mxidx(ij,kl)
      apraoint=aprscr(ijkl)
c
      return
      end

      subroutine chk_ga2(no,nv,nb)
      implicit REAL  (a-h,o-z)
INCLUDE(common/global)
c
      write(*,919)
919   format(/,1x,'checking global arrays',/)
c  (VV|VO)
      if (g_vvvo.ge.0)then
       print *,'vvvo'
       do i=1,nv
        iv=i+no
        do j=1,nv
         jv=j+no
         do k=1,nv
          kv=k+no
          do l=1,no
           gg1=aprmoint(iv,jv,kv,l)
           ij = (i-1)*nv + j
           lk = (l-1)*nv + k
           call pg_get(g_vvvo, ij, ij, lk, lk, gg2, 1)
           if (dabs(gg1-gg2).gt.1.0d-10)then
            write(*,453)'vvvo1',i,j,k,l,gg1,gg2,gg1-gg2
           endif
  453     format(a5,4i4,3e15.5)
          enddo
         enddo
        enddo
       enddo
      endif
 
c  (VV|OO)
      print *,'vvoo'
      do i=1,nv
       iv=i+no
       do j=1,i
        jv=j+no
        do k=1,no
         do l=1,k
          gg1=aprmoint(iv,jv,k,l)
          ij = i*(i-1)/2+j
          kl = k*(k-1)/2+l
          call pg_get(g_vvoo,ij,ij,kl,kl,gg2,1)
          if (dabs(gg1-gg2).gt.1.0d-10)then
           write(*,453)'vvoo1',i,j,k,l,gg1,gg2,gg1-gg2
          endif
         enddo
        enddo
       enddo
      enddo
c
c  (VO|VO)
      print *,'vovo'
      do i=1,nv
       iv=i+no
       do j=1,no
        do k=1,nv
         kv=k+no
         do l=1,j
          gg1=aprmoint(iv,j,kv,l)
          kl = j*(j-1)/2+l
          ij = (k-1)*nv + i
          call pg_get(g_vovo,ij,ij,kl,kl,gg2,1)
          if (dabs(gg1-gg2).gt.1.0d-10)then
           write(*,453)'vovo1',i,j,k,l,gg1,gg2,gg1-gg2
          endif
         enddo
        enddo
       enddo
      enddo
c
c  (VO|OO)
      print *,'vooo'
      do i=1,nv
       iv=i+no
       do j=1,no
        do k=1,no
         do l=1,k
          gg1=aprmoint(iv,j,k,l)
          ij = i
          kl =(k*(k-1)/2+l-1)*no+j
          call pg_get(g_vooo,ij,ij,kl,kl,gg2,1)
          if (dabs(gg1-gg2).gt.1.0d-10)then
           write(*,453)'vooo1',i,j,k,l,gg1,gg2,gg1-gg2
          endif
         enddo
        enddo
       enddo
      enddo

c  (OO|OO)
      print *,'oooo'
      do i=1,no
       do j=1,i
        do k=1,no
         do l=1,k
          gg1=aprmoint(i,j,k,l)
          ij = i*(i-1)/2+j
          kl = k*(k-1)/2+l
          call pg_get(g_oooo,ij,ij,kl,kl,gg2,1)
          if (dabs(gg1-gg2).gt.1.0d-10)then
           write(*,453)'oooo1',i,j,k,l,gg1,gg2,gg1-gg2
          endif
         enddo
        enddo
       enddo
      enddo
c
      return
      end
      subroutine ckaoint(no,nb)
      implicit REAL  (a-h,o-z)
      integer g_bbbo,ilbbbo,ihbbbo,jlbbbo,jhbbbo
      common/aprglob/g_bbbo,ilbbbo,ihbbbo,jlbbbo,jhbbbo
c
      print *,'ao ints'
      do i=1,nb
       do j=1,i
        do k=1,no
         do l=1,nb
          xx1=aprmoint(i,j,k,l)
          ij = i*(i-1)/2+j
          kl = (k-1)*nb+l
          call pg_get(g_bbbo,ij,ij,kl,kl,xx2,1)
          if (dabs(xx1-xx2).gt.1.0d-10)then
          write(*,843)i,j,k,l,xx1,xx2,xx1-xx2
  843     format('diff ',4i5,3e15.5)
          endif
         enddo
        enddo
       enddo
      enddo
      return
      end
      subroutine chklag(wmo,tmp1,tmp2,nocc,nvir,nbas,eorb,pmo)
      implicit none
      integer nocc,nvir,nbas
      REAL wmo(nbas,nbas),tmp1(nbas*nocc),tmp2(nbas*nocc),
     &       eorb(nbas),pmo(nbas,nbas)
      REAL xai,xacjb(255),xibjc,xicjb,ticjb,aprmoint
      integer av,i,cv,j,bv,iia,joff
c
      print *,' check wmo '
      iia=0
      do i=1,nocc
      do av=nocc+1,nocc+nvir
      iia=iia+1
      tmp2(iia)=wmo(av,i)
      enddo
      enddo
      call pg_dgop(4321,tmp2,nocc*nbas,'+')
c     
      call dcopy(nvir*nocc,0.0d00,0,tmp1,1)
      do j=1,nocc
      joff=(j-1)*nvir-nocc
      do cv=nocc+1,nocc+nvir
      print *,' j,cv ',j,cv
      do bv=nocc+1,nocc+nvir
      do av=nocc+1,nocc+nvir
       xacjb(av)=aprmoint(av,cv,j,bv)
       tmp1(joff+bv)=tmp1(joff+bv)+4.0d00*xacjb(av)*pmo(av,cv)
       tmp1(joff+av)=tmp1(joff+av)-xacjb(av)*pmo(bv,cv)
       tmp1(joff+cv)=tmp1(joff+cv)-xacjb(av)*pmo(av,bv)
      enddo
      iia=0
      do i=1,nocc
      xicjb=aprmoint(i,cv,j,bv)
      xibjc=aprmoint(i,bv,j,cv)
      ticjb=(8.0d00*xicjb-4.0d00*xibjc)/
     &      (eorb(bv)+eorb(cv)-eorb(i)-eorb(j))
      do av=nocc+1,nocc+nvir
       iia=iia+1
       tmp1(iia)=tmp1(iia)+xacjb(av)*ticjb
      enddo
      enddo
      enddo

      enddo
      enddo
      iia=0
      do i=1,nocc
      do av=nocc+1,nocc+nvir
       iia=iia+1
       write(*,843)av,i,tmp2(iia),tmp1(iia),tmp2(iia)-tmp1(iia)
  843  format(2i5,3e16.6)
      enddo
      enddo
      return
      end
      subroutine ckt2ao(g_t2ao,tmp1,cmo,eorb,nocc,nvir,nbas,buf)
      implicit none
      integer g_t2ao,nocc,nvir,nbas
      REAL tmp1(nocc,nocc),cmo(nbas,nbas),eorb(nbas),buf(*)
      integer iabv,a,b,i,j,av,bv
      REAL t2mo,xtst,aprmoint
c
      print *,'checking mkt2ao '
c
c code to check t2ao
      iabv=0
      do a=1,nbas
       do b=1,a
        iabv=iabv+1
c        call ga_get(g_t2ao,1,nocc*nocc,iabv,iabv,tmp1,nocc*nocc)
        call get_half_trans(tmp1,buf,nbas,nocc,nvir,a,b)
        do i = 1 , nocc
         do j = 1 , nocc
          xtst=0.0d00
          do av=nocc+1,nocc+nvir
           do bv=nocc+1,nocc+nvir
            t2mo=(2.0d00*aprmoint(i,av,j,bv)-aprmoint(i,bv,j,av))/
     &           (eorb(av)+eorb(bv)-eorb(i)-eorb(j))
c            print *,'a,b,i,j,av,bv,t2mo',a,b,i,j,av,bv,t2mo
c            print *,aprmoint(i,av,j,bv),aprmoint(i,bv,j,av),
c     &           eorb(av),eorb(bv),eorb(i),eorb(j)
            xtst=xtst+cmo(a,av)*cmo(b,bv)*4.0d00*t2mo
           enddo
          enddo
          write(*,84)a,b,i,j,tmp1(j,i),xtst,tmp1(j,i)-xtst
  84      format(4i5,3e16.6)
         enddo
        enddo
       enddo
      enddo
      return
      end
      subroutine wpij(pmat,wmat,eorb,nbas,nocc,nvir,i,j,k)
      implicit none
      integer i,j,k,a,b,nbas,nocc,nvir
      REAL pmat(nbas,nbas),wmat(nbas,nbas),eorb(nbas)
      REAL xiajb,tiajb,diajb,gtint,gt2
      do a=nocc+1,nbas
       do b=nocc+1,nbas
        xiajb=gtint(nocc,nvir,i,a,j,b)
        tiajb=gt2(nocc,nvir,i,a,k,b,eorb)
        diajb=eorb(a)+eorb(b)-eorb(i)-eorb(j)
        pmat(k,j)=pmat(k,j)+2.0d00*xiajb*tiajb/diajb
        wmat(k,j)=wmat(k,j)+4.0d00*xiajb*tiajb
       enddo
      enddo
      return
      end
      subroutine tst_aprwp(pmat,wmat,eorb,nocc,nvir,nbas,buf1,buf2,buf3)
      implicit none
      integer nocc,nvir,nbas
      REAL pmat(nbas,nbas),wmat(nbas,nbas),eorb(nbas),
     &       buf1(nbas*nbas),buf2(nbas*nbas),buf3(nbas*nbas)
INCLUDE(common/global)
c
      integer i,j,k,l,a,b,c,nop1,nvv,ab,ba,ij,ik,jk,ijk,ia,nov,
     &        noo,noo2,nvv2,kl
      integer ipg_nnodes,ipg_nodeid,nproc,me
      integer iskip1,iskip2,iskip3,iskip4
      REAL xiajb,tiajb,diajb,gtint,gt2
      character*1 xn, xt
      data xn,xt/'n','t'/
c
      nop1=nocc+1
      nvv=nvir*nvir
      nov=nvir*nocc
      noo=nocc*nocc
      noo2=nocc*(nocc+1)/2
      nvv2=nvir*(nvir+1)/2
      call vclr(wmat,1,nbas*nbas)
      call vclr(pmat,1,nbas*nbas)
      iskip1=0
      iskip2=0
      iskip3=0
      iskip4=0
c
c -- pvv1 and wvv1
      if (iskip1.ne.0)goto 123
      ij=0
      do i=1,nocc
       do j=1,i
        ij=ij+1
        if (ij.ge.jl_vovo.and.ij.le.jh_vovo) then
         call pg_get(g_vovo,1,nvv,ij,ij,buf2,nvv)
         ba=0
         do a=1,nvir
          ab=a-nvir
          do b=1,nvir
           ba=ba+1
           ab=ab+nvir
           buf1(ba)=(buf2(ab)+buf2(ab)-buf2(ba))/
     &              (eorb(nocc+a)+eorb(nocc+b)-eorb(i)-eorb(j))
          enddo
         enddo
         call dgemm(xn,xn,nvir,nvir,nvir
     &   ,4.0d00,buf2,nvir
     &   ,buf1,nvir,1.0d00,wmat(nop1,nop1),nbas)
         if (i.ne.j)then
          call dgemm(xt,xt,nvir,nvir,nvir
     &  ,4.0d00,buf2,nvir
     &  ,buf1,nvir,1.0d00,wmat(nop1,nop1),nbas)
         endif
         ba=0
         do a=1,nvir
          ab=a-nvir
          do b=1,nvir
           ba=ba+1
           ab=ab+nvir
           buf2(ba)=buf2(ba)/
     &              (eorb(nocc+a)+eorb(nocc+b)-eorb(i)-eorb(j))
          enddo
         enddo
         call dgemm(xt,xt,nvir,nvir,nvir
     &    ,-2.0d00,buf1,nvir
     &    ,buf2,nvir,1.0d00,pmat(nop1,nop1),nbas)
         if (i.ne.j)then
          call dgemm(xn,xn,nvir,nvir,nvir
     &    ,-2.0d00,buf1,nvir
     &    ,buf2,nvir,1.0d00,pmat(nop1,nop1),nbas)
         endif
        endif
       enddo
      enddo
  123 continue
c
      if (iskip1.ne.0)then
      do i=1,nocc
       do j=1,nocc
        do a=nop1,nbas
         do b=nop1,nbas
          do c=nop1,nbas
           xiajb=gtint(nocc,nvir,i,a,j,b)
           tiajb=gt2(nocc,nvir,i,a,j,c,eorb)
           diajb=eorb(a)+eorb(b)-eorb(i)-eorb(j)
           pmat(c,b)=pmat(c,b)-2.0d00*xiajb*tiajb/diajb
           wmat(b,c)=wmat(b,c)+4.0d00*xiajb*tiajb
          enddo
         enddo
        enddo
       enddo
      enddo
      endif
c      call wrtsqm_ga('W(2)',nbas, wmat)
c
c -- poo1 and woo1
      if (iskip2.ne.0)goto 124
      ij=0
      ijk=0
      nproc=ipg_nnodes()
      me=ipg_nodeid()
c*****SHOULD BE LOAD BALANCED
      do i=1,nocc
       do j=1,i
       ij=ij+1
        do k=1,j
         ijk=ijk+1
         if(mod(ijk,nproc).eq.me)then
          ik=i*(i-1)/2+k
          jk=j*(j-1)/2+k
          call pg_get(g_vovo,1,nvv,ij,ij,buf1,nvv)
          call pg_get(g_vovo,1,nvv,ik,ik,buf2,nvv)
          call pg_get(g_vovo,1,nvv,jk,jk,buf3,nvv)
c
c         call wpij(pmat,wmat,eorb,nbas,nocc,nvir,i,j,k)
          ba=0
          do a=1,nvir
           ab=a-nvir
           do b=1,nvir
            ba=ba+1
            ab=ab+nvir
            xiajb=buf1(ab)
            tiajb=(buf2(ab)+buf2(ab)-buf2(ba))/
     &            (eorb(nocc+a)+eorb(nocc+b)-eorb(i)-eorb(k))
            diajb=eorb(nocc+a)+eorb(nocc+b)-eorb(i)-eorb(j)
            pmat(k,j)=pmat(k,j)+2.0d0*xiajb*tiajb/diajb
            wmat(k,j)=wmat(k,j)+4.0d00*xiajb*tiajb
           enddo
          enddo
c 
          if (j.ne.k)then
c          call wpij(pmat,wmat,eorb,nbas,nocc,nvir,i,k,j)
           ba=0
           do a=1,nvir
            ab=a-nvir
            do b=1,nvir
             ba=ba+1
             ab=ab+nvir
             xiajb=buf2(ab)
             tiajb=(buf1(ab)+buf1(ab)-buf1(ba))/
     &             (eorb(nocc+a)+eorb(nocc+b)-eorb(i)-eorb(j))
             diajb=eorb(nocc+a)+eorb(nocc+b)-eorb(i)-eorb(k)
             pmat(j,k)=pmat(j,k)+2.0d0*xiajb*tiajb/diajb
             wmat(j,k)=wmat(j,k)+4.0d00*xiajb*tiajb
            enddo
           enddo
          endif
c 
          if (i.ne.j)then
c           call wpij(pmat,wmat,eorb,nbas,nocc,nvir,j,i,k)
           ba=0
           do a=1,nvir
            ab=a-nvir
            do b=1,nvir
             ba=ba+1
             ab=ab+nvir
             xiajb=buf1(ba)
             tiajb=(buf3(ab)+buf3(ab)-buf3(ba))/
     &             (eorb(nocc+a)+eorb(nocc+b)-eorb(j)-eorb(k))
             diajb=eorb(nocc+a)+eorb(nocc+b)-eorb(i)-eorb(j)
             pmat(k,i)=pmat(k,i)+2.0d0*xiajb*tiajb/diajb
             wmat(k,i)=wmat(k,i)+4.0d00*xiajb*tiajb
            enddo
           enddo
c
            if (i.ne.k)then
c             call wpij(pmat,wmat,eorb,nbas,nocc,nvir,j,k,i)
             ba=0
             do a=1,nvir
              ab=a-nvir
              do b=1,nvir
               ba=ba+1
               ab=ab+nvir
               xiajb=buf3(ab)
               tiajb=(buf1(ba)+buf1(ba)-buf1(ab))/
     &               (eorb(nocc+a)+eorb(nocc+b)-eorb(j)-eorb(i))
               diajb=eorb(nocc+a)+eorb(nocc+b)-eorb(k)-eorb(j)
               pmat(i,k)=pmat(i,k)+2.0d0*xiajb*tiajb/diajb
               wmat(i,k)=wmat(i,k)+4.0d00*xiajb*tiajb
              enddo
             enddo
            endif
          endif
c 
          if (k.ne.i.and.k.ne.j)then
c           call wpij(pmat,wmat,eorb,nbas,nocc,nvir,k,j,i)
           ba=0
           do a=1,nvir
            ab=a-nvir
           do b=1,nvir
             ba=ba+1
             ab=ab+nvir
             xiajb=buf3(ba)
             tiajb=(buf2(ba)+buf2(ba)-buf2(ab))/
     &             (eorb(nocc+a)+eorb(nocc+b)-eorb(k)-eorb(i))
             diajb=eorb(nocc+a)+eorb(nocc+b)-eorb(k)-eorb(j)
             pmat(i,j)=pmat(i,j)+2.0d0*xiajb*tiajb/diajb
             wmat(i,j)=wmat(i,j)+4.0d00*xiajb*tiajb
            enddo
           enddo
           if (i.ne.j)then
c            call wpij(pmat,wmat,eorb,nbas,nocc,nvir,k,i,j)
            ba=0
            do a=1,nvir
             ab=a-nvir
             do b=1,nvir
              ba=ba+1
              ab=ab+nvir
              xiajb=buf2(ba)
              tiajb=(buf3(ba)+buf3(ba)-buf3(ab))/
     &              (eorb(nocc+a)+eorb(nocc+b)-eorb(k)-eorb(j))
              diajb=eorb(nocc+a)+eorb(nocc+b)-eorb(k)-eorb(i)
              pmat(j,i)=pmat(j,i)+2.0d0*xiajb*tiajb/diajb
              wmat(j,i)=wmat(j,i)+4.0d00*xiajb*tiajb
             enddo
            enddo
           endif
          endif
c
         endif
        enddo
       enddo
      enddo
c
  124 continue
      if (iskip2.ne.0)then
      do i=1,nocc
       do j=1,nocc
        do a=nop1,nbas
         do b=nop1,nbas
          do k=1,nocc
           xiajb=gtint(nocc,nvir,i,a,j,b)
           tiajb=gt2(nocc,nvir,i,a,k,b,eorb)
           diajb=eorb(a)+eorb(b)-eorb(i)-eorb(j)
           pmat(k,j)=pmat(k,j)+2.0d00*xiajb*tiajb/diajb
           wmat(k,j)=wmat(k,j)+4.0d00*xiajb*tiajb
          enddo
         enddo
        enddo
       enddo
      enddo
      endif
c
c -- wov1
      if (iskip3.ne.0)goto 125
      ij=0
      do k=1,nocc
       do j=1,k
        ij=ij+1
        if (ij.ge.jl_vovo.and.ij.le.jh_vovo) then
         call pg_get(g_vovo,1,nvv,ij,ij,buf2,nvv)
         ba=0
         do a=1,nvir
          ab=a-nvir
          do b=1,nvir
           ba=ba+1
           ab=ab+nvir
           buf1(ba)=(buf2(ab)+buf2(ab)-buf2(ba))/
     &              (eorb(nocc+a)+eorb(nocc+b)-eorb(k)-eorb(j))
          enddo
         enddo
c
         ia=1
         do i=1,nocc
          if (i.ge.k)then
           jk=(i*(i-1)/2+k-1)*nocc+j
          else
           jk=(k*(k-1)/2+i-1)*nocc+j
          endif
          call pg_get(g_vooo,1,nvir,jk,jk,buf2(ia),nvir)
          ia=ia+nvir
         enddo
c -- wov1
         call dgemm(xt,xn,nocc,nvir,nvir
     &    ,4.0d00,buf2,nvir
     &    ,buf1,nvir,1.0d00,wmat(1,nop1),nbas)
c -- wvo2
         ia=0
         do i=1,nocc
          do a=nop1,nbas
           ia=ia+1
           xiajb=buf2(ia)
           wmat(a,j)=wmat(a,j)+xiajb*pmat(i,k)*4.0d00
           wmat(a,i)=wmat(a,i)-xiajb*pmat(j,k)
           wmat(a,k)=wmat(a,k)-xiajb*pmat(j,i)
          end do
         end do
         if (j.ne.k)then
          ia=1
          do i=1,nocc
           if (i.ge.j)then
            jk=(i*(i-1)/2+j-1)*nocc+k
           else
            jk=(j*(j-1)/2+i-1)*nocc+k
           endif
           call pg_get(g_vooo,1,nvir,jk,jk,buf2(ia),nvir)
           ia=ia+nvir
          enddo
c -- wov1
          call dgemm(xt,xt,nocc,nvir,nvir
     &    ,4.0d00,buf2,nvir
     &    ,buf1,nvir,1.0d00,wmat(1,nop1),nbas)
c -- wvo2
          ia=0
          do i=1,nocc
           do a=nop1,nbas
            ia=ia+1
            xiajb=buf2(ia)
            wmat(a,k)=wmat(a,k)+xiajb*pmat(i,j)*4.0d00
            wmat(a,i)=wmat(a,i)-xiajb*pmat(k,j)
            wmat(a,j)=wmat(a,j)-xiajb*pmat(k,i)
           end do
          end do
         endif
c
        endif
       enddo
      enddo
  125 continue
c
c -- wov1
      if (iskip3.ne.0)then
      do i=1,nocc
       do j=1,nocc
        do a=nop1,nbas
         do b=nop1,nbas
          do k=1,nocc
           xiajb=gtint(nocc,nvir,i,k,j,b)
           tiajb=gt2(nocc,nvir,k,a,j,b,eorb)
           wmat(i,a)=wmat(i,a)+4.0d00*xiajb*tiajb
          enddo
         enddo
        enddo
       enddo
      enddo
c
c -- wvo2
      do i=1,nocc
       do j=1,nocc
        do k=1,nocc
         do a=nop1,nbas
          xiajb=gtint(nocc,nvir,a,i,j,k)
          wmat(a,i)=wmat(a,i)+xiajb*pmat(j,k)*4.0d00
          wmat(a,j)=wmat(a,j)-xiajb*pmat(i,k)
          wmat(a,k)=wmat(a,k)-xiajb*pmat(i,j)
         end do
        end do
       end do
      end do
      endif
c
c -- wvo3
c      do i=1,nocc
c       do j=1,nocc
c        do a=nop1,nbas
c         do b=nop1,nbas
c          do c=nop1,nbas
c           xiajb=gtint(nocc,nvir,a,c,b,j)
c           tiajb=gt2(nocc,nvir,i,c,j,b,eorb)
c           wmat(a,i)=wmat(a,i)+4.0d00*xiajb*tiajb
c          enddo
c         enddo
c        enddo
c       enddo
c      enddo
c
c -- wvo4
c      do i=1,nocc

c       do a=nop1,nbas
c        do b=nop1,nbas
c         do c=nop1,nbas
c          xiajb=gtint(nocc,nvir,a,i,b,c)
c          wmat(a,i)=wmat(a,i)+xiajb*pmat(b,c)*4.0d00
c          wmat(b,i)=wmat(b,i)-xiajb*pmat(a,c)
c          wmat(c,i)=wmat(c,i)-xiajb*pmat(a,b)
c         end do
c        end do
c       end do
c      end do
c
      if (iskip4.ne.0)goto 126
      ij=0
      do i=1,nocc
       do j=1,i
        ij=ij+1
        if (ij.ge.jl_oooo.and.ij.le.jh_oooo) then
         call pg_get(g_oooo,1,noo2,ij,ij,buf1,noo2)
         call pg_get(g_vvoo,1,nvv2,ij,ij,buf2,nvv2)
         call pg_get(g_vovo,1,nvv,ij,ij,buf3,nvv)
         do k=1,nocc
          do l=1,nocc
           kl=k*(k-1)/2+l
           if (k.lt.l)kl=l*(l-1)/2+k
           xiajb=buf1(kl)
           wmat(i,j)=wmat(i,j)+xiajb*pmat(k,l)*4.0d00
           wmat(i,k)=wmat(i,k)-xiajb*pmat(j,l)
           wmat(i,l)=wmat(i,l)-xiajb*pmat(j,k)
          enddo
        enddo
        do a=1,nvir
          do b=1,nvir
           ab=a*(a-1)/2+b
           if (a.lt.b)ab=b*(b-1)/2+a
           xiajb=buf2(ab)
           wmat(i,j)=wmat(i,j)+xiajb*pmat(a+nocc,b+nocc)*4.0d00
           ab=(b-1)*nvir+a
           xiajb=buf3(ab)
           wmat(i,j)=wmat(i,j)-xiajb*pmat(a+nocc,b+nocc)
           wmat(j,i)=wmat(j,i)-xiajb*pmat(a+nocc,b+nocc)
          enddo
        enddo
c
        if (i.ne.j)then
         do k=1,nocc
           do l=1,nocc
            kl=k*(k-1)/2+l
            if (k.lt.l)kl=l*(l-1)/2+k
            xiajb=buf1(kl)
            wmat(j,i)=wmat(j,i)+xiajb*pmat(k,l)*4.0d00
            wmat(j,k)=wmat(j,k)-xiajb*pmat(i,l)
            wmat(j,l)=wmat(j,l)-xiajb*pmat(i,k)
           enddo
         enddo
c -- woo3
         do a=1,nvir
          do b=1,nvir
            ab=a*(a-1)/2+b
            if (a.lt.b)ab=b*(b-1)/2+a
            xiajb=buf2(ab)
            wmat(j,i)=wmat(j,i)+xiajb*pmat(a+nocc,b+nocc)*4.0d00
            ab=(a-1)*nvir+b
            xiajb=buf3(ab)
            wmat(j,i)=wmat(j,i)-xiajb*pmat(a+nocc,b+nocc)
            wmat(i,j)=wmat(i,j)-xiajb*pmat(a+nocc,b+nocc)
           enddo
          enddo
         endif

        endif
       enddo
      enddo
  126 continue
c
c
c -- woo2
      if (iskip4.ne.0)then
      do i=1,nocc
       do j=1,nocc
        do k=1,nocc
          do l=1,nocc
           xiajb=gtint(nocc,nvir,i,j,k,l)
           wmat(i,j)=wmat(i,j)+xiajb*pmat(k,l)*4.0d00
           wmat(i,k)=wmat(i,k)-xiajb*pmat(j,l)
           wmat(i,l)=wmat(i,l)-xiajb*pmat(j,k)
          enddo
        enddo
       enddo
      enddo
c
c -- woo3
      do i=1,nocc
       do j=1,nocc
        do a=nop1,nbas
          do b=nop1,nbas
           xiajb=gtint(nocc,nvir,i,j,a,b)
           wmat(i,j)=wmat(i,j)+xiajb*pmat(a,b)*4.0d00
           xiajb=gtint(nocc,nvir,i,a,j,b)
           wmat(i,j)=wmat(i,j)-xiajb*pmat(a,b)
           wmat(j,i)=wmat(j,i)-xiajb*pmat(a,b)
          enddo
        enddo
       enddo
      enddo
      endif
c
c
cwoo4
      do i=1,nocc
       do j=1,nocc
c        wmat(j,i)=wmat(j,i)+pmat(j,i)*(eorb(i)+eorb(j))
       enddo
      enddo
c
cwvv4
      do a=nop1,nbas
       do b=nop1,nbas
c        wmat(b,a)=wmat(b,a)+pmat(b,a)*(eorb(a)+eorb(b))
       enddo
      enddo
c
      return
      end
      double precision function gtint(nocc,nvir,i,j,k,l)
      implicit none
INCLUDE(common/global)
      integer nocc,nvir,i,j,k,l,ii,jj,kk,ll,iijj,kkll,itmp
      integer ij,kl
      REAL gg2
      integer ibt
      ibt(i,j)=i*(i-1)/2+j
c
      ii=i
      jj=j
      if (j.gt.i)then
       ii=j
       jj=i
      endif
      kk=k
      ll=l
      if (l.gt.k)then
       kk=l
       ll=k
      endif
      iijj=ii*(ii-1)/2+jj
      kkll=kk*(kk-1)/2+ll
c      print *,'ii,jj,kk,ll',ii,jj,kk,ll
      if (kkll.gt.iijj)then
       itmp=ii
       ii=kk
       kk=itmp
       itmp=jj
       jj=ll
       ll=itmp
      endif
c      print *,'ijkl ',i,j,k,l,ii,jj,kk,ll
c
      if (ii.le.nocc)then
c oooo
c       gg1=aprmoint(i,j,k,l)
       ij = ibt(ii,jj)
       kl = ibt(kk,ll)
c       print *,'oooo ',ij,kl
       call pg_get(g_oooo, ij, ij, kl, kl, gg2, 1)
      else
       if (jj.le.nocc)then
        if (kk.le.nocc)then
c vooo
c         gg1=aprmoint(iv,j,k,l)
         ij = ii-nocc
         kl = (ibt(kk,ll)-1)*nocc+jj
c         print *,'vooo ',ij,kl
         call pg_get(g_vooo, ij,ij, kl, kl, gg2, 1)
        else
         if (ll.le.nocc)then
c vovo
c          gg1=aprmoint(iv,j,kv,l)
          if (jj.ge.ll)then
           ij = (kk-nocc-1)*nvir+ii-nocc
           kl = ibt(jj,ll)
          else
           ij = (ii-nocc-1)*nvir+kk-nocc
           kl = ibt(ll,jj)
          endif
c          print *,'vovo ',ij,kl
          call pg_get(g_vovo, ij, ij, kl, kl, gg2, 1)
         else
c vovv
c          gg1=aprmoint(iv,jv,kv,l)
          print *,'vvvo ',ij,kl
          stop
         endif
        endif
       else
        if (kk.le.nocc)then
c vvoo
c          gg1=aprmoint(iv,jv,k,l)
          ij = ibt(ii-nocc,jj-nocc)
          kl = ibt(kk,ll)
c          print *,'vvoo ',ij,kl
          call pg_get(g_vvoo, ij, ij, kl, kl, gg2, 1)
        else
         if (ll.le.nocc)then
c vvvo
          print *,'vvvo ',ij,kl
          stop
c          call pg_get(g_vvvo, ij, ij, lk, lk, gg2, 1)
         else
c vvvv
         endif
        endif
       endif
      endif
      gtint=gg2
      return
      end
      double precision function gt2(nocc,nvir,i,a,j,b,denom)
      implicit none
INCLUDE(common/global)
      integer nocc,nvir,i,a,j,b
      REAL denom(nocc+nvir)
      integer ij,ab,ba
      REAL xiajb,xibja
c vovo
      if (i.ge.j)then
       ij = i*(i-1)/2+j
       ab = (b-nocc-1)*nvir+a-nocc
       ba = (a-nocc-1)*nvir+b-nocc
      else
       ij = j*(j-1)/2+i
       ab = (a-nocc-1)*nvir+b-nocc
       ba = (b-nocc-1)*nvir+a-nocc
      endif
c      print *,'t2 ',i,a,j,b,ij,ab,ba
      call pg_get(g_vovo, ab, ab, ij, ij, xiajb, 1)
      call pg_get(g_vovo, ba, ba, ij, ij, xibja, 1)
      gt2=(xiajb+xiajb-xibja)/(denom(a)+denom(b)-denom(i)-denom(j))
      return
      end
_ENDIF

      subroutine parclenms(messge)
c
      implicit REAL  (a-h,o-z)
      character messge*(*)
c
INCLUDE(common/iofile)
c
      l = len(messge)
      write (iwr,'(/1x,a)') messge(1:l)
      call clenup
      call pg_end(0)
      end
c
      subroutine chk_ga_tools
      implicit REAL  (a-h,o-z)
INCLUDE(common/global)
c
      write(*,919)
919   format(/,1x,'checking global arrays',/)
c  (VO|VO)
      i = 1
      j = 1
      k = 1
      l = 1
      call pg_get(g_vovo,i,j,k,l,gg2,1)
      print *,'vovo'
      write(6,*)'vovo1',i,j,k,l,gg2
c
      return
      end

      subroutine ga_mp2_memchk(iwr)
      implicit none
INCLUDE(common/errcodes)
INCLUDE(common/restar)

      integer igmem_max_memory
      external igmem_max_memory

      integer ipg_nnodes
      external ipg_nnodes

      integer iwr, np, memrqd, nmax, itest, i,nprhold

      np = ipg_nnodes()

      memrqd = 0

      call ga_mp2_memory(iwr,np,memrqd)
c
c  check requirement against memory allocated
c  this check is performed on all nodes as 
c  a precaution against problems with future distributed
c  data schemes. 
c
      if (memrqd.eq.0 ) return

      nmax = igmem_max_memory()

      write(iwr,100) memrqd
      write(iwr,200) nmax

      if (memrqd.le.nmax) then
         itest = 1
      else
         itest = 0
      endif

      call pg_igop(1010,itest,1,'+')

      if (itest.eq.np) then
         return
      else
c
c before aborting, give summary for different numbers of nodes
c
        nprint = -5

        write(iwr,*)
        write(iwr,*)'Insufficient core available for the parallel MP2'
        write(iwr,*)'The following may help revise allocation:'
        write(iwr,*)
        write(iwr,*)'nodes    requirement (MW)'
        np = 1
        do i = 1, 11
           call ga_mp2_memory(iwr,np,memrqd)
           write(iwr,300)np,memrqd
           np = np*2
        enddo

        call gamerr('insufficient core for parallel MP2',
     &       ERR_NO_CODE, ERR_USER_DIMENSION, ERR_SYNC, ERR_NO_SYS)

      endif
 100  format(/1x,'**** Memory required: ',i10)
 200  format(/1x,'**** Memory available: ',i10)
 300  format(1x,i5,i9)
      end


      subroutine ga_mp2_memory(iwr,nproc,memrqd)
c
c  check if enough memory has been allocated for global arrays
c  to procede with parallel mp2 energy or gradient, or scf hessian
c
c  9/12/96: now that all GA and gamess core memory are handled together
c  this test is not so reliable as core used by optimiser etc will 
c  reduce available core
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/cndx41)
      logical skipp
      common/lsort/skipp(3*maxat)
INCLUDE(common/parcntl)
INCLUDE(common/nshel)
INCLUDE(common/symtry)
INCLUDE(common/runlab)
INCLUDE(common/restrl)
INCLUDE(common/cntl1)

INCLUDE(common/segm)
c
      if (((zruntp.eq.'gradient'.or.zruntp.eq.'optimize'.or.
     +      zruntp.eq.'force'   .or.zruntp.eq.'optxzy'  .or.
     +      zruntp.eq.'scf').and.mp2) .or. 
     +    (zruntp.eq.'hessian'.and.zscftp.eq.'rhf')) then
        call set41
        no = nocca
        nv = nvirta
        n = ncoorb
        nn = n*(n+1)/2
        noo = no*(no+1)/2
        nov = no*nv
        nsq = n*n
        ns = nshell

        maxl=4
        do i = 1 , ns
           if(ktype(i).eq.3) then
            maxl = max(maxl,6)
           else if(ktype(i).eq.4) then
            maxl = max(maxl,10)
           else if(ktype(i).eq.5) then
            maxl = max(maxl,15)
           endif
        enddo
      else
        return
      end if
c
c  for mp2 gradient case 
c
      if ((zruntp.eq.'gradient'.or.
     +    (zruntp.eq.'optimize'.and.isadle.ne.2). or.
     +     zruntp.eq.'force'   .or.zruntp.eq.'optxyz').
     +    and.mp2) then

c  workspace allocated before entering driver + MA overhead
        memws = ns*(ns+1)/2 
     &+         nsq 
     &+         nsq 
     &+         nsq 
     &+         nsq 
     &+         n 
     &+         8*6
c
c  evaluate requirements at memory-critical steps of calc.
c
c  1. transformation
c
        memga = noo*noo                                !  oooo
     &+         nn*noo                                 !  vvoo
     &+         nov*noo                                !  vooo
     &+         nsq*noo                                !  vovo
        memloc = nsq 
     &+          (ns*nt-1)/lenwrd()+1 
     &+          maxl**4 
     &+          maxl*maxl*n*no 
     &+          n 
     &+          max( nsq , n*maxl*maxl ) 
     &+          max( nsq , n*maxl*maxl ) 
        memtrn = (memga + 1)/nproc 
     &+          memloc
c  8-word overhead required by MA tools for each allocation
     &+          8*11

c
c  2. VVV - integrals
c
        memga = nn*noo                                 !  vvoo
     &+         nov*noo                                !  vooo
     &+         nsq*noo                                !  vovo
        memloc1 = nsq 
     &+           (ns*nt-1)/lenwrd()+1 
     &+           maxl**4 
     &+           2*maxl*maxl*n*no 
     &+           no*no 
     &+           no*n 
     &+           noo
        memloc2 = nsq
     &+           nsq
     &+           nsq
        memloc  = max( memloc1 , memloc2 ) 
        memvvv  = (memga + 1)/nproc 
     &+           memloc
c  8-word MA overhead for each allocation
     &+           8*10
c
c  maximum memory requirement for mp2 gradient/optimisation
c
        memrqd = max( memtrn , memvvv ) + memws
c
c  scf hessian run
c
      else if (zruntp.eq.'hessian'.and.zscftp.eq.'rhf') then

c  workspace allocated before entering driver + MA overhead
         memws = ns*(ns+1)/2 
     &+          nsq 
     &+          nsq 
     &+          nsq 
     &+          nsq 
     &+          n 
     &+          8*6
c  find number of perturbations
         npi = 0
         do i = npstar + 1 , npfin
           if (.not.skipp(i)) npi = npi + 1
         end do

         memga = noo*noo                           ! oooo
     &+          nn*noo                            ! vvoo
     &+          nov*noo                           ! vooo
     &+          nsq*noo                           ! vovo
     &+          nov*50*npi                        ! trial vecs
     &+          2*nov*npi                         ! soln vecs
         memloc1 = nsq 
     &+            (ns*nt-1)/lenwrd()+1 
     &+            maxl**4 
     &+            maxl*maxl*n*no 
     &+            n 
     &+            max( nsq , n*maxl*maxl ) 
     &+            max( nsq , n*maxl*maxl ) 
         memloc2 = n 
     &+            nov*npi 
     &+            nov*npi 
     &+            nov*npi 
     &+            nov*npi 
     &+            nov*npi 
     &+            2650 
     &+            nv*nv
         memloc  = max( memloc1 , memloc2 ) 
c
c  maximum memory required for scf hessian calc.
c
         memrqd = (memga + 1)/nproc 
     &+           memloc
c  plus MA overhead for each array
     &+           8*15
c
c  mp2 energy run 
c
      else if ((zruntp.eq.'scf '.or.
     +         (zruntp.eq.'optimize'.and.isadle.eq.2))
     +          .and.mp2) then

        memws = ns*(ns+1)/2 
     &+         nsq 
     &+         nsq 
     &+         nsq 
     &+         nsq 
     &+         n 
     &+         8*6
        memga = nsq*noo                                !  vovo
        memloc = nsq 
     &+          (ns*nt-1)/lenwrd()+1 
     &+          maxl**4 
     &+          maxl*maxl*n*no 
     &+          n 
     &+          max( nsq , n*maxl*maxl ) 
     &+          max( nsq , n*maxl*maxl ) 
        memtrn = (memga + 1)/nproc 
     &+          memloc
c  8-word overhead required by MA tools for each allocation
     &+          8*8
        memrqd = memtrn 
      end if
      end  

      subroutine ver_mp2_parallel(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/mp2_parallel.m,v $
     +     "/
      data revision /"$Revision: 6176 $"/
      data date /"$Date: 2010-08-10 16:49:47 +0200 (Tue, 10 Aug 2010) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
