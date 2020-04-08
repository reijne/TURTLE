c 
c  $Author: jmht $
c  $Date: 2008-10-06 23:36:09 +0200 (Mon, 06 Oct 2008) $
c  $Locker:  $
c  $Revision: 5733 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/dft/integ_te3c_rep.m,v $
c  $State: Exp $
c ********
c ********
c *** these integral routines have been extracted from GAMESS-UK,
c *** and modified to work as far as possible independently of the 
c *** parent code.
c *** As they now stand, they function (i) only in direct-mode.
c *** (ii) assuming that symmetry has been disabled (i.e must be run 
c *** with nosym), and supermatrices have also been disabled.
c *******
      subroutine te3c_rep_drv(icd_tag,iao_tag,gout,amatrix,
     &                        oVform_sw,
     &                        eri_scr,dvec_scr,rmatrix_out)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
***************************************************************
* EXTERNAL routines invoked
**************************************************************
*
* Parallel load balancing etc
*
*     pg_dlbchunk, pg_dlbreset, pg_dlbpush, ipg_nodes
*
*
* Memory allocation routines
*
*     setscm, setc, cmem, loccm
*
* BLAS - linear algebra
*
*      daxpy, vlcr, vadd, setsto
*
* ERROR monitoring
*
*      caserr
*
* TIMING + statistics + ulilities
*
*      cpuwal, timana, cpulft, writel
***************************************************************
***************************************************************
* common blocks used
**************************************************************
**** from main-line code
*********************************************************************
* cslosc * files * infoa *  parallel * parcntl 
* prints * restar * sizes * statis * timez
***************************************************************
* local common blocks
*
* tabinx, junkx, mapperx, flipsx, flip7x, ijlabx
* pqgeomx, qgeomx, pgeomx, geomx, shlg70x, piconx, miscgx, astorex
* ginfx, bshellx, constx, shllfox, maxcx, savecx, typex
* shlinfx
* shltx, miscx, setintx, rootx, denssx, shlnosx, indezx, rtdatx
* inxblkx
********************************************************************
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/cslosc)
INCLUDE(../m4/common/infoa)
INCLUDE(../m4/common/statis)
INCLUDE(common/dft_iofile)
INCLUDE(../m4/common/restar)
INCLUDE(../m4/common/timez)
c
INCLUDE(../m4/common/parallel)
INCLUDE(../m4/common/parcntl)
c
INCLUDE(common/dft_mbasis)
INCLUDE(common/dft_shlt)
INCLUDE(common/dft_mapper)
INCLUDE(common/dft_ijlab)
INCLUDE(common/dft_picon)
INCLUDE(common/dft_root)
c
      dimension gout(*)
      dimension amatrix(*),dvec_scr(*)
      dimension rmatrix_out(*),eri_scr(*)
c
      data done/1.0d0/,two/2.0d0/,twopt5/2.5d0/,four/4.0d0/
c
      call cpuwal(begin,ebegin)
c ***
c *** establish arrays normally set in main line code
c ***
c
c ... /mapperx/
c
      do 30 i = 1 , maxorb
         k = i*(i-1)/2
         iky(i) = k
 30   continue
c
c ... /piconx/
c
      pidiv4 = datan(done)
      pi = four*pidiv4
      pito52 = two*pi**twopt5
      root3 = dsqrt(3.0d0)
      root5 = dsqrt(5.0d0)
      root53= root5/root3
      root7 = dsqrt(7.0d0)
c
c ... /rootx/
c
      do 20 loop = 1 , 60
         dji(loop) = 1.0d0/(loop+loop-1)
 20   continue
c
c ... /auxvarx/
c
c      ofast = .false.
c
c
c /cslosc/
c
      call setsto(10,0,intcut)
      call setsto(1060,0,intmag)
      nopk = 1
c
_IF(parallel)
c***   ***node-MPP***
c .. estimate load-balancing factor
      l2 = nshell(iao_tag)*(nshell(iao_tag)+1)/2
      if(ntchnk.le.0)then
         write(6,*)'chunk factor not initialised',ipg_nnodes(), 
     &        ntchnk, limchnk
         call caserr('problem with parallel set-up')
      endif
      nodet = (l2-1)/ipg_nnodes() + 1
      if (nodet.ge.ntchnk) then
         ichunk = nodet / ntchnk
      else
         ichunk = 1
      endif
      if(limchnk.gt.0)ichunk = min(ichunk,limchnk)
      call pg_dlbchunk(ichunk,(nprint.ne.-5))
      call pg_dlbreset

c***   ***node-MPP***
_ENDIF
c***
c     check for pure sp basis. if so, do gaussian integrals.
      call spchck_3c
      call debut_3c
      nopkr = nopk
      iofrst = iofsym
c
      nindmx = 1

      call aux_find(iao_tag)

      call jkint_gamess_3c(iao_tag,icd_tag,gout,
     &                     oVform_sw,amatrix,dvec_scr,
     &                     eri_scr,rmatrix_out)
      call final_3c
      cpu = cpulft(1)
c ***
         if (outv) then
            write (iwr,6010) (intcut(i),i=1,3)
            if (oimag) then
               do 50 k = 1 , 11
                  isum = 0
                  do 40 kk = (k-1)*4 + 1 , (k-1)*4 + 4
                     isum = isum + intmag(kk)
 40               continue
                  intmag(k) = isum
 50            continue
               write (iwr,6020) (k-34,k=1,42,4)
               write (iwr,6030) (intmag(k),k=1,11)
            end if
         end if
c     if (nprint.ne.-5) write (iwr,6040) cpu

ccc      call setc(loadcm)

      call timana(4)
      return
 6010 format (/,' integral test counts'/1x,30('=')
     +        /' on ij shell         ',i10/' on ijkl shells      ',
     +        i10/' on ijkl shells & den',i10/1x,30('='))
 6020 format (' magnitudes of computed integrals '//' 2**',4x,11i8)
 6030 format (' ',8x,11i8)
 6040 format (/' end of 2-electron integral evaluation at ',f8.2,
     +        ' seconds')
      end
      subroutine jkint_gamess_3c(iao_tag,icd_tag,
     &                           gout,oVform_sw,amatrix,dvec_scr,
     &                           eri_scr,rmatrix_out)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(../m4/common/sizes)
      parameter (mxp2 = mxprms * mxprms)
INCLUDE(../m4/common/cslosc)
INCLUDE(../m4/common/restar)
INCLUDE(common/dft_iofile)
INCLUDE(common/dft_mbasis)
c
INCLUDE(common/dft_mapper)
INCLUDE(common/dft_shlnos)
INCLUDE(common/dft_shlt)
INCLUDE(common/dft_ijlab)
INCLUDE(common/dft_shlg70)
INCLUDE(common/dft_picon)
INCLUDE(common/dft_misc)
INCLUDE(common/dft_auxg)
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
      common/junkx/cxyz(3,5625),aaa(21*mxp2),ijaaa(225)
c
INCLUDE(../m4/common/timez)
INCLUDE(../m4/common/parallel)
c
INCLUDE(common/dft_flips)
INCLUDE(common/dft_shlinf)
c
      dimension ib(4,4)
      dimension gout(*)
      dimension amatrix(*),dvec_scr(*)
      dimension eri_scr(225,*),rmatrix_out(*)
      REAL ovtest_3c,testlim,IL_shlove_tol
      character*1 xn,xt
      data xn,xt/'n','t'/

      data ib/64,16,4,1,216,36,6,1,1000,100,10,1,
     +         3375,225,15,1 /
c
c     ----- two-electron integrals -----
c
      testlim=IL_shlove_tol()
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
c     ----- set some parameters -----
c
c
c
      l2 = iky(numorb(1)+1)
c
      time = cpulft(1)
      tim0 = time
      tim1 = time
c
      ist0 = ist
      jst0 = jst
      kst0 = kst
      lst0 = lst

      ogauss = .true.
_IF(parallel)
c***   **MPP**
      next = ipg_dlbtask()
c      write(6,*)'next', next
c***   **MPP**
_ENDIF
c
c     ----- ishell -----
c
      locij_dm=1
      locij_km=0
_IF(parallel)
      do 150 ii = nshell(iao_tag), ist0 , -1
_ELSE
      do 150 ii = 1 , nshell(iao_tag)
_ENDIF
c
c     ----- print intermediate restart data -----
c
         ishell = ii
         mini = kmin(iao_tag,ishell)
         maxi = kmax(iao_tag,ishell)
         loci = kloc(iao_tag,ishell)-mini
         ibas_num=(maxi-mini)+1
c        kadi = kad(iao_tag,ii)
         dt0 = time - tim0
         dt1 = time - tim1
         tim1 = time
         ikyii = iky(ii)
c
c     ----- jshell -----
c
         j0 = jst0
_IF(parallel)
         do 140 jj = ii, j0, -1
_ELSE
         do 140 jj = 1 , ii
_ENDIF
c           write(6,*) ' '
c           write(6,*) 'New Shell',ii,jj
            jst0 = 1
            minj = kmin(iao_tag,jj)
            maxj = kmax(iao_tag,jj)
            itrij = ikyii + jj
c           kadij = kad(iao_tag,jj) + kadi
_IF(parallel)
c***   **MPP**
            icount_dlb = icount_dlb + 1
c            write(6,*)'icount_dlb,next',icount_dlb,next
            if(icount_dlb . eq. next) then

c              write(6,*)'task', icount_dlb,' on node',ipg_nodeid(),
c     &              ' ij ',ii,jj

c***   **MPP**
_ENDIF
c
c     ----- get information about i-shell and j-shell -----
c
            jshell = jj
            locj = kloc(iao_tag,jshell)-minj
            jbas_num=(maxj-minj)+1
            call shells_3c(gout,iao_tag,icd_tag,
     &                      1,ishell,jshell,ishell,jshell,1)
c
c switch 170 to 140
c
            if(ovtest_3c(ii,jj,rri).lt.testlim)goto 170
            call ijprim_3c
c           write(6,*)'nij',nij

            if (nij.eq.0) go to 170


            iaobas_num=ibas_num*jbas_num
            imc=0
c
c     ----- kshell -----
c
               k0 = kst0
               ipos=0
               do 130 kk = 1, nshell(icd_tag)
                  kst0 = 1
c                 kadijk = kad(icd_tag,kk) + kadij
                  ikykk = iky(kk)
                  kshell = kk
                  mink = kmin(icd_tag,kshell)
                  maxk = kmax(icd_tag,kshell)
                  icdbas_num=(maxk-mink)+1
                  lock = kloc(icd_tag,kshell)-mink
                  itrjk = iky(max(jj,kk)) + min(jj,kk)
c
                  q4 = 1.0d0
                  lshell = 0
                  qq4 = q4
C *
C *For shell quartet with no sp shell(s)
C *
c
c     ----- get information about ksh and lsh -----
c

                  call shells_3c(gout,iao_tag,icd_tag,
     &                           2,ishell,jshell,kshell,lshell,1)

                  call genral_3c(gout)

                              write(6,*)'mat',ii,jj,kk,mini,maxi,
     &                             minj,maxj,mink,maxk,
     &                         iaobas_num,nbasfn(icd_tag),imc

                  call mat_form_3c(gout,mini,maxi,minj,maxj,mink,maxk,
     &                             iaobas_num,nbasfn(icd_tag),imc,
     &                             eri_scr)
                  imc=imc+icdbas_num
 130           continue

               if(oVform_sw) then
                 call aclear_dp(dvec_scr,iaobas_num,0.0d0)
                 call dgemv(xn,iaobas_num,nbasfn(icd_tag),1.0d0
     &                      ,eri_scr,225,amatrix,1,1.0d0
     &                      ,dvec_scr,1)
                 locij_t=locij_km
                 call km_fill(locij_t,ii,jj,
     &                        rmatrix_out,dvec_scr)
                 locij_km=locij_km+locij_t
               else

                 fac=2.0d0
                 if(ii.eq.jj) fac=1.0d0
                 call dvec_fill(locij_dm,
     &                          amatrix,ii,jj,mini,maxi,minj,maxj,
     &                          dvec_scr)

                 call dgemv(xt,iaobas_num,nbasfn(icd_tag),fac
     &                      ,eri_scr,225
     &                      ,dvec_scr,1,1.0d0,rmatrix_out,1)

               endif
 170           continue

_IF(parallel)
            next = ipg_dlbtask()
            endif
_ENDIF
 140     continue
         time = cpulft(1)
 150     continue
_IF(parallel)

         call pg_dlbpush
c
c  globally sum accumulated quantities
c
         if(oVform_sw) then
            nlen =  ((nbasfn(iao_tag)+1)*nbasfn(iao_tag))/2
         else
            nlen =  nbasfn(icd_tag)
         endif
c         write(6,*)'dgop',nlen
         call pg_dgop(2001,rmatrix_out,nlen,'+')
_ENDIF
c

      return
 6020 format (i4,3i5,1x,i10,9x,f11.2,f9.2)
      end
      subroutine spchck_3c
c
c     ----- check for pure sp basis -----
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/restar)
INCLUDE(common/dft_mbasis)
c
INCLUDE(common/dft_ijlab)
c
      ospbas = .true.
      opdbas = .false.
      opfbas = .false.
      opgbas = .false.
      do lbas=1,2
        do 20 i = 1 , nshell(lbas)
           kad(lbas,i) = 0
           ii = i
           if (ktype(lbas,i).gt.2) then
             ospbas = .false.
             if (ktype(lbas,i).eq.3) opdbas = .true.
             if (ktype(lbas,i).eq.4) opfbas = .true.
             if (ktype(lbas,i).eq.5) opgbas = .true.
             kad(lbas,i) = -1
           end if
 20     continue
        if (intg76.eq.0) then
          do 30 i = 1 , nshell(lbas)
             kad(lbas,i) = -1
 30       continue
        end if
      enddo
      return
      end
      subroutine debut_3c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/cslosc)
INCLUDE(../m4/common/prints)
INCLUDE(../m4/common/infoa)
INCLUDE(../m4/common/restar)
INCLUDE(common/dft_iofile)
INCLUDE(../m4/common/files)
INCLUDE(../m4/common/timez)
c
INCLUDE(common/dft_mbasis)
INCLUDE(common/dft_shlt)
INCLUDE(common/dft_bshell)
INCLUDE(common/dft_intctl)
      data done,ten,e /1.0d0,1.0d+01,2.30258d0/
c
      outv = oprint(59)
      if (nprint.eq.-5) outv = .false.
      out = nprint.eq.4
c
      omaxb = .false.

c      icut0 = iabs(icut)
c      if (icut0.eq.0) icut0 = 9
c      cutoff = done/(ten**icut0)
c      if (itol.eq.0) itol = 20
c      tol = e*itol

      icut0 = iabs(icut_3c)
      cutoff = done/(ten**icut0)
      tol = e*itol_3c
      if (nprint.ne.-5 .and. opg_root()) then
        write(iwr,*)'Using cut,tol',icut_3c,itol_3c
      endif

c
c     -----  irest = 0   -----
c     ----- normal start -----
c
      mfilep = 1
      mainp = n2tape(1)
      iblkmp = n2blk(1)
      mblp = iblkmp - n2last(1)
      m2file = 1
      do 30 i = 2 , 20
         m2blk(i) = -1
 30   continue
      m2tape(1) = mainp
      m2blk(1) = iblkmp
      m2last(1) = -1
c
      do 40 i = 1 , nshell(1)
         icc = katom(1,i)
         p(i) = c(1,icc)
         q(i) = c(2,icc)
         r(i) = c(3,icc)
 40   continue
      ist = 1
      jst = 1
      kst = 1
      lst = 1
      nindmx = 0
      nrec = 1
      icount = 1
      ic4 = 1
      cpu = cpulft(1)
      isti = ist
      jsti = jst
      ksti = kst
      lsti = lst
      lastb = iblkmp
      lastu = mfilep
      if (nprint.ne.-5 .and. opg_root()) then
         write (iwr,6060) cpu
         if (outv) write (iwr,6070)
         else
      end if
      return
 6060 format (/' commence 2-e 3-c integral evaluation at ',f8.2,
     +        ' seconds')
 6070 format (/1x,'ist',2x,'jst',2x,'kst',2x,'lst',7x,'nrec',3x,
     +        'intloc',5x,'del(t)',5x,'time'/1x,58('-'))
      end
      subroutine final_3c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(../m4/common/sizes)
INCLUDE(common/dft_iofile)
INCLUDE(../m4/common/restar)
c
INCLUDE(common/dft_shlt)
INCLUDE(common/dft_shlnos)
INCLUDE(common/dft_misc)
c ***
c
      if (icount.gt.1) then
          nrec = nrec + 1
      end if
c
      m2last(m2file) = m2last(m2file) + 1
      irest = 0
      nrec = 1
      ist = 1
      jst = 1
      kst = 1
      lst = 1
c
      cpu = cpulft(1)
      if (nprint.ne.-5 .and. opg_root()) then
         write(iwr,100)cpu
      endif

 100  format(' 3-c 2-e integral evaluation complete at ',f8.2,
     +        ' seconds')
      return
      end
      subroutine indxa_3c(ijx,ijy,ijz,ij,mini,maxi,
     &minj,maxj,iandj,inc1,inc2,inc3)
      implicit REAL  (a-h,o-z)
      logical iandj
      dimension ijx(225),ijy(225),ijz(225)
INCLUDE(common/dft_inxblk)
      dimension jx(35),jy(35),jz(35)
c
_IF1(ct)cdir$ novector
_IF1(a)cvd$  novector
_IF1(x)c$dir scalar
      do 20 j = minj , maxj
         jx(j) = ix(j)*inc2
         jy(j) = iy(j)*inc2
         jz(j) = iz(j)*inc2
 20   continue
      ij = 0
      jmax = maxj
      do 40 i = mini , maxi
         nx = ix(i)*inc1 + inc3
         ny = iy(i)*inc1 + inc3
         nz = iz(i)*inc1 + inc3
         if (iandj) jmax = i
         do 30 j = minj , jmax
            ij = ij + 1
            ijx(ij) = nx + jx(j)
            ijy(ij) = ny + jy(j)
            ijz(ij) = nz + jz(j)
 30      continue
 40   continue
_IF1(ct)cdir$ vector
      return
      end
      subroutine genral_3c(gout)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
INCLUDE(../m4/common/sizes)
      parameter (mxp2 = mxprms * mxprms)
INCLUDE(../m4/common/restar)
c
INCLUDE(common/dft_shlinf)
      common/junkx/pin(5625), qin(5625), rin(5625),
     + aa(mxp2),r(mxp2),p1(mxp2),q1(mxp2),r1(mxp2),dd(16*mxp2),
     + ijd(225)
INCLUDE(common/dft_shlnos)
INCLUDE(common/dft_root)
INCLUDE(common/dft_shlt)
INCLUDE(common/dft_misc)
INCLUDE(common/dft_setint)
INCLUDE(common/dft_denss)
      dimension gout(*),in1(9)
c
      data ijn1,ijn2,kln1,kln2 /125,25,5,1/
      data sqrt3 /1.73205080756888d0/
      data sqrt5 /2.23606797749979d0/
      data sqrt7 /2.64575131106459d0/
      data pi252 /34.986836655250d0/
      data pt5,done /0.5d0,1.0d0/
c
      if (ijkl.eq.1) then
c
c     ----- (s,s//s,s) -----
c
         call s0000_3c(gout)
      else if (ij.eq.1) then
c
c     ----- (s,s//k,l) -----
c
         call sskl_3c(gout)
      else if (kl.eq.1) then
         call ijss_3c(gout)
      else
         factor = pi252*qq4
         onorm = normf.ne.1 .or. normp.ne.1
         onormk = .false.
         onorml = .false.
c
c     ----- select expansion centre for -xyz- integrals -----
c
         if (lit.lt.ljt) then
            ni = ljt - 1
            nj = lit - 1
            ij1 = ijn2
            ij2 = ijn1
            pc = pj
            qc = qj
            rc = rj
            dxij = pj - pi
            dyij = qj - qi
            dzij = rj - ri
         else
            ni = lit - 1
            nj = ljt - 1
            ij1 = ijn1
            ij2 = ijn2
            pc = pi
            qc = qi
            rc = ri
            dxij = pi - pj
            dyij = qi - qj
            dzij = ri - rj
         end if
         if (lkt.lt.llt) then
            nk = llt - 1
            nl = lkt - 1
            kl1 = kln2
            kl2 = kln1
            pd = pl
            qd = ql
            rd = rl
            dxkl = pl - pk
            dykl = ql - qk
            dzkl = rl - rk
         else
            nk = lkt - 1
            nl = llt - 1
            kl1 = kln1
            kl2 = kln2
            pd = pk
            qd = qk
            rd = rk
            dxkl = pk - pl
            dykl = qk - ql
            dzkl = rk - rl
         end if
         nmax = ni + nj
         mmax = nk + nl
         max = nmax + 1
         do 20 i = 1 , max
            n = i - 1
            if (n.le.ni) in1(i) = ij1*n + 1
            if (n.gt.ni) in1(i) = ij1*ni + ij2*(n-ni) + 1
 20      continue
         max = mmax + 1
         do 30 k = 1 , max
            n = k - 1
            if (n.le.nk) kn(k) = kl1*n
            if (n.gt.nk) kn(k) = kl1*nk + kl2*(n-nk)
 30      continue
c
c     ----- k primitive
c
         lgmax = ngd
         do 270 kg = 1 , ngc
            ak = cgg(kg)
            brrk = ak*rrk
            akxk = ak*pk
            akyk = ak*qk
            akzk = ak*rk
            csk = csc(kg)*factor
            cpk = cpc(kg)*factor
            cdk = cdc(kg)*factor
            cfk = cfc(kg)*factor
            cgk = cgc(kg)*factor
c
c     ----- l primitive
c
            if (okanl) lgmax = kg
            do 260 lg = 1 , 1
               al = dg(lg)
               b = ak + al
               b1 = done/b
               bbrrk = al*brrk*b1
               if (bbrrk.le.tol) then
                  csl = csd(lg)
                  cpl = cpd(lg)
                  cdl = cdd(lg)
                  cfl = cfd(lg)
                  cgl = cgd(lg)
                  pb = (akxk+al*pl)*b1
                  qb = (akyk+al*ql)*b1
                  rb = (akzk+al*rl)*b1
                  bxbd = b*(pb-pd)
                  bybd = b*(qb-qd)
                  bzbd = b*(rb-rd)
                  bxbc = b*(pb-pc)
                  bybc = b*(qb-qc)
                  bzbc = b*(rb-rc)
c
c     ----- density factor
c
                  odoub = okanl .and. kg.gt.lg
                  n = 0
                  max = maxl
                  do 210 k = mink , maxk
                     go to (40,50,110,110,
     +                 60,110,110,70,110,110,
     +                 80,110,110,90,110,110,110,110,110,100,
     +                 92,110,110,94,110,110,110,110,110,96,
     +                110,110,98,110,110), k
 40                  dum1 = csk*b1
                     go to 110
 50                  dum1 = cpk*b1
                     go to 110
 60                  dum1 = cdk*b1
                     go to 110
 70                  if (onormk) dum1 = dum1*sqrt3
                     go to 110
 80                  dum1 = cfk*b1
                     go to 110
 90                  if (onormk) dum1 = dum1*sqrt5
                     go to 110
 100                 if (onormk) dum1 = dum1*sqrt3
                     go to 110
 92                  dum1 = cgk*b1
                     go to 110
 94                  if (onormk) dum1 = dum1*sqrt7
                     go to 110
 96                  if (onormk) dum1 = dum1*sqrt5/sqrt3
                     go to 110
 98                  if (onormk) dum1 = dum1*sqrt3
 110                 if (okanl) max = k
                     do 200 l = 1,1
                        go to (120,130,190,190,
     +                    140,190,190,150,190,190,
     +                    160,190,190,170,190,190,190,190,190,180,
     +                    182,190,190,184,190,190,190,190,190,186,
     +                    190,190,188,190,190), l
c
 120                    dum2 = dum1*csl
                        if (odoub) then
                           if (k.gt.1) then
                              dum2 = dum2 + csk*cpl*b1
                           else
                              dum2 = dum2 + dum2
                           end if
                        end if
                        go to 190
 130                    dum2 = dum1*cpl
                        if (odoub) dum2 = dum2 + dum2
                        go to 190
 140                    dum2 = dum1*cdl
                        if (odoub) dum2 = dum2 + dum2
                        go to 190
 150                    if (onorml) dum2 = dum2 * sqrt3
                        go to 190
 160                    dum2 = dum1*cfl
                        if (odoub) dum2 = dum2 + dum2
                        go to 190
 170                    if (onorml) dum2 = dum2 * sqrt5
                        go to 190
 180                    if (onorml) dum2 = dum2 * sqrt3
                        go to 190
 182                    dum2 = dum1*cgl
                        if (odoub) dum2 = dum2+dum2
                        go to 190
 184                   if (onorml) dum2 = dum2 * sqrt7
                        go to 190
 186                   if (onorml) dum2 = dum2 * sqrt5/sqrt3
                        go to 190
 188                   if (onorml) dum2 = dum2 * sqrt3
 190                    n = n + 1
                        dkl(n) = dum2
 200                 continue
 210              continue
c
c     ----- pair of i,j primitives
c
                  nn = 0
                  do 250 n = 1 , nij
                     dum = bbrrk + r(n)
                     if (dum.le.tol) then
                        do 220 i = 1 , ij
                           dij(i) = dd(ijd(i)+nn)
 220                    continue
                        a = aa(n)
                        ab = a*b
                        aandb = a + b
                        expe = dexp(-dum)/dsqrt(aandb)
                        rho = ab/aandb
                        pa = p1(n)
                        qa = q1(n)
                        ra = r1(n)
                        pp = rho*((pa-pb)**2+(qa-qb)**2+(ra-rb)**2)
                        axad = a*(pa-pd)
                        ayad = a*(qa-qd)
                        azad = a*(ra-rd)
                        axac = a*(pa-pc)
                        ayac = a*(qa-qc)
                        azac = a*(ra-rc)
                        c1x = bxbd + axad
                        c2x = a*bxbd
                        c3x = bxbc + axac
                        c4x = b*axac
                        c1y = bybd + ayad
                        c2y = a*bybd
                        c3y = bybc + ayac
                        c4y = b*ayac
                        c1z = bzbd + azad
                        c2z = a*bzbd
                        c3z = bzbc + azac
                        c4z = b*azac
c
c     ----- roots and weights for quadrature
c
                        if (nroots.le.3) call rt123_dft
                        if (nroots.eq.4) call roots4_dft
                        if (nroots.eq.5) call roots5_dft
                        if (nroots.ge.6) call roots_dft
                        mm = 0
                        max = nmax + 1
                        do 240 m = 1 , nroots
                           u2 = u(m)*rho
                           f00 = expe*w(m)
                           do 230 i = 1 , max
                              in(i) = in1(i) + mm
 230                       continue
                           dum = done/(ab+u2*aandb)
                           dum2 = pt5*dum
                           bp01 = (a+u2)*dum2
                           b00 = u2*dum2
                           b10 = (b+u2)*dum2
                           pcp00 = (u2*c1x+c2x)*dum
                           pc00 = (u2*c3x+c4x)*dum
                           qcp00 = (u2*c1y+c2y)*dum
                           qc00 = (u2*c3y+c4y)*dum
                           rcp00 = (u2*c1z+c2z)*dum
                           rc00 = (u2*c3z+c4z)*dum
                           call xyzint_3c
                           mm = mm + 625
 240                    continue
c
c     ----- form (i,j//k,l) integrals over functions
c
                        call spdint_3c(gout)
                     end if
                     nn = nn + 16
 250              continue
               end if
 260        continue
 270     continue
      end if
      return
      end
      subroutine ijprim_3c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(../m4/common/sizes)
      parameter (mxp2 = mxprms * mxprms)
INCLUDE(../m4/common/restar)
c
INCLUDE(common/dft_shlt)
INCLUDE(common/dft_shlinf)
INCLUDE(common/dft_misc)
INCLUDE(common/dft_shlnos)
      common/junkx/cxyz(3,5625),
     + a(mxp2),r(mxp2),p1(mxp2),q1(mxp2),r1(mxp2),dij(16*mxp2),
     + ijd(225)
c
      data sqrt3 /1.73205080756888d0/
      data sqrt5 /2.23606797749979d0/
      data sqrt7 /2.64575131106459d0/
      data done /1.0d0/
c
      if (ij.eq.1) then
c
c     ----- (s,s// -----
c
         call ssprm_3c
         return
      else
         onorm = normf.ne.1 .or. normp.ne.1
         max = maxj
         n = 0
         nn = 0
         do 70 i = mini , maxi
            go to (20,20,30,30,
     +             20,30,30,20,30,30,
     +             20,30,30,20,30,30,30,30,30,20,
     +             20,30,30,20,30,30,30,30,30,20,
     +             30,30,20,30,30) , i
 20         nm = nn
 30         nn = nm
            if (oianj) max = i
            do 60 j = minj , max
               go to (40,40,50,50,
     +                40,50,50,40,50,50,
     +                40,50,50,40,50,50,50,50,50,40,
     +                40,50,50,40,50,50,50,50,50,40,
     +                50,50,40,50,50) , j
 40            nn = nn + 1
 50            n = n + 1
               ijd(n) = nn
 60         continue
 70      continue
c
c     ----- i primitive
c
         nij = 0

         write(6,*)'ijprim',nga,ngb,oianj

         jbmax = ngb
         do 310 ia = 1 , nga
            ai = ag(ia)
            arri = ai*rri
            axi = ai*pi
            ayi = ai*qi
            azi = ai*ri
            csi = csa(ia)
            cpi = cpa(ia)
            cdi = cda(ia)
            cfi = cfa(ia)
            cgi = cga(ia)
c
c     ----- j primitive
c
            if (oianj) jbmax = ia
            do 300 jb = 1 , jbmax
               aj = bg(jb)
               aa = ai + aj
               aa1 = done/aa
               dum = aj*arri*aa1
               if (dum.gt.tol) go to 300
               csj = csb(jb)
               cpj = cpb(jb)
               cdj = cdb(jb)
               cfj = cfb(jb)
               cgj = cgb(jb)
               nm = nij*16
               nn = nm
               nij = nij + 1
               r(nij) = dum
               a(nij) = aa
               p1(nij) = (axi+aj*pj)*aa1
               q1(nij) = (ayi+aj*qj)*aa1
               r1(nij) = (azi+aj*rj)*aa1
c
c     ----- density factor
c
               do 250 i = mini , maxi
                  go to (80,90,250,250,
     +                   100,250,250,110,250,250,
     +                   120,250,250,130,250,250,250,250,250,140,
     +                   142,250,250,144,250,250,250,250,250,146,
     +                   250,250,148,250,250) , i
c
 80               dum1 = csi*aa1
                  go to 150
 90               dum1 = cpi*aa1
                  go to 150
 100              dum1 = cdi*aa1
                  go to 150
 110              if (onorm) dum1 = dum1*sqrt3
                  go to 150
 120              dum1 = cfi*aa1
                  go to 150
 130              if (onorm) dum1 = dum1*sqrt5
                  go to 150
 140              if (onorm) dum1 = dum1*sqrt3
                  go to 150
 142              dum1 = cgi*aa1
                  go to 150
 144              if (onorm) dum1 = dum1*sqrt7
                  go to 150
 146              if (onorm) dum1 = dum1*sqrt5/sqrt3
                  go to 150
 148              if (onorm) dum1 = dum1*sqrt3
 150              if (oianj) max = i
                  do 240 j = minj , max
                     go to (160,170,240,240,
     +                      180,240,240,190,240,240,
     +                      200,240,240,210,240,240,240,240,240,220,
     +                      222,240,240,224,240,240,240,240,240,226,
     +                      240,240,228,240,240),j
 160                 dum2 = dum1*csj
                     go to 230
 170                 dum2 = dum1*cpj
                     go to 230
 180                 dum2 = dum1*cdj
                     go to 230
 190                 if (onorm) dum2 = dum2*sqrt3
                     go to 230
 200                 dum2 = dum1*cfj
                     go to 230
 210                 if (onorm) dum2 = dum2*sqrt5
                     go to 230
 220                 if (onorm) dum2 = dum2*sqrt3
                     go to 230
 222                 dum2 = dum1*cgj
                     go to 230
 224                 if (onorm) dum2 = dum2*sqrt7
                     go to 230
 226                 if (onorm) dum2 = dum2*sqrt5/sqrt3
                     go to 230
 228                 if (onorm) dum2 = dum2*sqrt3
 230                 nn = nn + 1
                     dij(nn) = dum2
 240              continue
 250           continue
               if (.not.oianj) go to 300
               if (ia.eq.jb) go to 300
               go to (290,260,280,270,275) , lit
 260           if (mini.ne.2) then
                  dij(nm+2) = dij(nm+2) + csi*cpj*aa1
                  dij(nm+3) = dij(nm+3) + dij(nm+3)
               end if
               go to 290
 275           dij(nm+10)= dij(nm+10)+ dij(nm+10)
               dij(nm+9) = dij(nm+9) + dij(nm+9)
               dij(nm+8) = dij(nm+8) + dij(nm+8)
               dij(nm+7) = dij(nm+7) + dij(nm+7)
 270           dij(nm+6) = dij(nm+6) + dij(nm+6)
               dij(nm+5) = dij(nm+5) + dij(nm+5)
               dij(nm+4) = dij(nm+4) + dij(nm+4)
 280           dij(nm+2) = dij(nm+2) + dij(nm+2)
               dij(nm+3) = dij(nm+3) + dij(nm+3)
 290           dij(nm+1) = dij(nm+1) + dij(nm+1)
 300        continue
 310     continue
         return
      end if
      end
      subroutine ijss_3c(gout)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(../m4/common/sizes)
      parameter (mxp2 = mxprms * mxprms)
INCLUDE(../m4/common/restar)
c
INCLUDE(common/dft_shlinf)
      common/junkx/pint(5625),qint(5625),rint(5625),
     + aa(mxp2),r(mxp2),p1(mxp2),q1(mxp2),r1(mxp2),dd(16*mxp2),
     + ijd(225)
INCLUDE(common/dft_shlnos)
INCLUDE(common/dft_root)
INCLUDE(common/dft_shlt)
INCLUDE(common/dft_misc)
INCLUDE(common/dft_indez)
INCLUDE(common/dft_setint)
INCLUDE(common/dft_denss)
      dimension gout(*)
      data ijn1,ijn2 /125,25/
      data pi252 /34.986836655250d0/
      data dzero,pt5,done /0.0d0,0.5d0,1.0d0/
      factor = pi252*qq4
c
c     ----- select expansion centre for -xyz- integrals -----
c
      if (lit.lt.ljt) then
         ni = ljt - 1
         nj = lit - 1
         ij1 = ijn2
         ij2 = ijn1
         pc = pj
         qc = qj
         rc = rj
         dxij = pj - pi
         dyij = qj - qi
         dzij = rj - ri
      else
         ni = lit - 1
         nj = ljt - 1
         ij1 = ijn1
         ij2 = ijn2
         pc = pi
         qc = qi
         rc = ri
         dxij = pi - pj
         dyij = qi - qj
         dzij = ri - rj
      end if
      nmax = ni + nj
      max = nmax + 1
      do 20 i = 1 , max
         n = i - 1
         if (n.le.ni) in(i) = ij1*n + 1
         if (n.gt.ni) in(i) = ij1*ni + ij2*(n-ni) + 1
 20   continue
c
c     ----- k primitive
c
      lgmax = ngd
      do 200 kg = 1 , ngc
         ak = cgg(kg)
         brrk = ak*rrk
         akxk = ak*pk
         akyk = ak*qk
         akzk = ak*rk
         csk = csc(kg)*factor
c
c     ----- l primitive
c
         if (okanl) lgmax = kg
         do 190 lg = 1 , lgmax
            al = dg(lg)
            b = ak + al
            b1 = done/b
            bbrrk = al*brrk*b1
            if (bbrrk.le.tol) then
               csl = csd(lg)
               pb = (akxk+al*pl)*b1
               qb = (akyk+al*ql)*b1
               rb = (akzk+al*rl)*b1
               bxbc = b*(pb-pc)
               bybc = b*(qb-qc)
               bzbc = b*(rb-rc)
c
c     ----- density factor
c
               d2 = csk*csl*b1
               if (okanl .and. (kg.gt.lg)) d2 = d2 + d2
c
c     ----- pair of i,j primitives
c
               nn = 0
               do 180 n = 1 , nij
                  dum = bbrrk + r(n)
                  if (dum.gt.tol) then
                     nn = nn + 16
                     go to 180
                  else
                     do 30 i = 1 , ij
                        dij(i) = dd(ijd(i)+nn)
 30                  continue
                     a = aa(n)
                     ab = a*b
                     aandb = a + b
                     expe = d2*dexp(-dum)/dsqrt(aandb)
                     rho = ab/aandb
                     pa = p1(n)
                     qa = q1(n)
                     ra = r1(n)
                     pp = rho*((pa-pb)**2+(qa-qb)**2+(ra-rb)**2)
                     axac = a*(pa-pc)
                     ayac = a*(qa-qc)
                     azac = a*(ra-rc)
                     c3x = bxbc + axac
                     c4x = b*axac
                     c3y = bybc + ayac
                     c4y = b*ayac
                     c3z = bzbc + azac
                     c4z = b*azac
c
c     ----- roots and weights for quadrature
c
                     if (nroots.le.3) then
                        call rt123_dft
                     else if (nroots.eq.4) then
                        call roots4_dft
                     else
                        call roots5_dft
                     end if
                     mm = 0
c
c     compute two-electron  integrals for each root
c
                     do 90 m = 1 , nroots
                        u2 = u(m)*rho
                        f00 = expe*w(m)
                        dum = done/(ab+u2*aandb)
                        b10 = (b+u2)*pt5*dum
                        pc00 = (u2*c3x+c4x)*dum
                        qc00 = (u2*c3y+c4y)*dum
                        rc00 = (u2*c3z+c4z)*dum
c
c     ----- i(0,0) -----
c
                        i1 = in(1) + mm
                        pint(i1) = done
                        qint(i1) = done
                        rint(i1) = f00
                        if (nmax.eq.0) then
                           mm = mm + 625
                           go to 90
                        else
c
c     ----- i(1,0) -----
c
                           i2 = in(2) + mm
                           pint(i2) = pc00
                           qint(i2) = qc00
                           rint(i2) = rc00*f00
                           if (nmax.le.1) then
                              mm = mm + 625
                              go to 90
                           else
c
c     ----- i(iii,0) -----
c
                              c10 = dzero
                              i3 = i1
                              i4 = i2
                              do 40 iii = 2 , nmax
                                 c10 = c10 + b10
                                 i5 = in(iii+1) + mm
                                 pint(i5) = c10*pint(i3) + pc00*pint(i4)
                                 qint(i5) = c10*qint(i3) + qc00*qint(i4)
                                 rint(i5) = c10*rint(i3) + rc00*rint(i4)
                                 i3 = i4
                                 i4 = i5
 40                           continue
                              if (nj.eq.0) then
                                 mm = mm + 625
                                 go to 90
                              else
c
c     ----- i(iii,jjj,0,0) -----
c
                                 i5 = in(nmax+1) + mm
                                 min = ni
                              end if
                           end if
                        end if
 50                     iii = nmax
                        i3 = i5
 60                     i4 = in(iii) + mm
                        pint(i3) = pint(i3) + dxij*pint(i4)
                        qint(i3) = qint(i3) + dyij*qint(i4)
                        rint(i3) = rint(i3) + dzij*rint(i4)
                        i3 = i4
                        iii = iii - 1
                        if (iii.gt.min) go to 60
                        min = min + 1
                        if (min.lt.nmax) go to 50
                        if (ni.ne.0) then
                           i3 = ij2 + i1
                           do 80 jjj = 1 , nj
                              i4 = i3
                              do 70 iii = 1 , ni
                                 pint(i4) = pint(i4+ij1-ij2)
     +                              + dxij*pint(i4-ij2)
                                 qint(i4) = qint(i4+ij1-ij2)
     +                              + dyij*qint(i4-ij2)
                                 rint(i4) = rint(i4+ij1-ij2)
     +                              + dzij*rint(i4-ij2)
                                 i4 = i4 + ij1
 70                           continue
                              i3 = i3 + ij2
 80                        continue
                        end if
                        mm = mm + 625
 90                  continue
c
c     ----- form (i,j//k,l) integrals over functions
c
                     go to (100,120,140,160,165) , nroots
                  end if
 100              do 110 iii = 1 , ij
                     nx = ijx(iii)
                     ny = ijy(iii)
                     nz = ijz(iii)
                     nnn = ijgt(iii)
                     gout(nnn) = (pint(nx)*qint(ny)*rint(nz))
     +                           *dij(iii) + gout(nnn)
 110              continue
                  nn = nn + 16
                  go to 180
 120              do 130 iii = 1 , ij
                     nx = ijx(iii)
                     ny = ijy(iii)
                     nz = ijz(iii)
                     nnn = ijgt(iii)
                     gout(nnn) = (pint(nx)*qint(ny)*rint(nz)+
     +                       pint(nx+625)*qint(ny+625)*rint(nz+625))
     +                       *dij(iii) + gout(nnn)
 130              continue
                  nn = nn + 16
                  go to 180
 140              do 150 iii = 1 , ij
                     nx = ijx(iii)
                     ny = ijy(iii)
                     nz = ijz(iii)
                     nnn = ijgt(iii)
                     gout(nnn) = (pint(nx)*qint(ny)*rint(nz)+
     +                      pint(nx+625)*qint(ny+625)*rint(nz+625)+
     +                      pint(nx+1250)*qint(ny+1250)*rint(nz+1250))
     +                      *dij(iii) + gout(nnn)
 150              continue
                  nn = nn + 16
                  go to 180
 160              do 170 iii = 1 , ij
                     nx = ijx(iii)
                     ny = ijy(iii)
                     nz = ijz(iii)
                     nnn = ijgt(iii)
                     gout(nnn) = (pint(nx)*qint(ny)*rint(nz)+
     +                       pint(nx+625)*qint(ny+625)*rint(nz+625)+
     +                       pint(nx+1250)*qint(ny+1250)*rint(nz+1250)+
     +                       pint(nx+1875)*qint(ny+1875)*rint(nz+1875))
     +                       *dij(iii) + gout(nnn)
 170              continue
                  nn = nn + 16
                  go to 180
 165              do 175 iii = 1 , ij
                     nx = ijx(iii)
                     ny = ijy(iii)
                     nz = ijz(iii)
                     nnn = ijgt(iii)
                     gout(nnn) = (pint(nx)*qint(ny)*rint(nz)+
     +                       pint(nx+625)*qint(ny+625)*rint(nz+625)+
     +                       pint(nx+1250)*qint(ny+1250)*rint(nz+1250)+
     +                       pint(nx+1875)*qint(ny+1875)*rint(nz+1875)+
     +                       pint(nx+2500)*qint(ny+2500)*rint(nz+2500))
     +                       *dij(iii) + gout(nnn)
 175              continue
                  nn = nn + 16
 180           continue
            end if
 190     continue
 200  continue
      return
      end
      subroutine s0000_3c(g)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(../m4/common/sizes)
      parameter (mxp2 = mxprms * mxprms)
INCLUDE(common/dft_shlinf)
INCLUDE(common/dft_shlt)
      common/junkx/cxyz(3,5625),
     + a(mxp2),r(mxp2),p1(mxp2),q1(mxp2),r1(mxp2),dij(16*mxp2),
     + ijd(225)
INCLUDE(common/dft_shlnos)
INCLUDE(common/dft_misc)
      dimension g(*)
      data pi252 /34.986836655250d0/
      data pie4 /7.85398163397448d-01/
      data dzero,done /0.0d0,1.0d0/
      gout = dzero
      lgmax = ngd
      do 40 kg = 1 , ngc
         bk = cgg(kg)
         brrk = bk*rrk
         bxk = bk*pk
         byk = bk*qk
         bzk = bk*rk
         csk = csc(kg)
         if (okanl) lgmax = kg
         do 30 lg = 1 , lgmax
            bl = dg(lg)
            bb = bk + bl
            bb1 = done/bb
            dum = bl*brrk*bb1
            if (dum.le.tol) then
               bbrrk = dum
               d2 = csd(lg)*csk*bb1
               if (okanl .and. lg.ne.kg) d2 = d2 + d2
               bbx = (bxk+bl*pl)*bb1
               bby = (byk+bl*ql)*bb1
               bbz = (bzk+bl*rl)*bb1
               sum = dzero
               nn = 1
               do 20 n = 1 , nij
                  dum = bbrrk + r(n)
                  if (dum.le.tol) then
                     expe = dexp(-dum)
                     aa = a(n)
                     ab = done/(aa+bb)
                     dum = p1(n) - bbx
                     pp = dum*dum
                     dum = q1(n) - bby
                     pp = dum*dum + pp
                     dum = r1(n) - bbz
                     pp = dum*dum + pp
                     p = pp*aa*bb*ab
                     if (p.gt.5.0d0) then
                        if (p.le.15.0d0) then
                           e = dexp(-p)
                           if (p.gt.10.0d0) then
                              ww1 = (((-1.8784686463512d-01/p+
     +                              2.2991849164985d-01)
     +                              /p-4.9893752514047d-01)
     +                              /p-2.1916512131607d-05)
     +                              *e + dsqrt(pie4/p)
                           else
                              ww1 = ((((((4.6897511375022d-01/p-
     +                              6.9955602298985d-01)
     +                              /p+5.3689283271887d-01)
     +                              /p-3.2883030418398d-01)
     +                              /p+2.4645596956002d-01)
     +                              /p-4.9984072848436d-01)
     +                              /p-3.1501078774085d-06)
     +                              *e + dsqrt(pie4/p)
                           end if
                        else if (p.gt.33.0d0) then
                           ww1 = dsqrt(pie4/p)
                        else
                           e = dexp(-p)
                           ww1 = ((1.9623264149430d-01/p-
     +                           4.9695241464490d-01)
     +                           /p-6.0156581186481d-05)
     +                           *e + dsqrt(pie4/p)
                        end if
                     else if (p.gt.1.0d0) then
                        if (p.gt.3.0d0) then
                           q = p - 4.0d0
                           f1 = ((((((((((-2.62453564772299d-11*q+
     +                          3.24031041623823d-10)
     +                          *q-3.614965656163d-09)
     +                          *q+3.760256799971d-08)
     +                          *q-3.553558319675d-07)
     +                          *q+3.022556449731d-06)
     +                          *q-2.290098979647d-05)
     +                          *q+1.526537461148d-04)
     +                          *q-8.81947375894379d-04)
     +                          *q+4.33207949514611d-03)
     +                          *q-1.75257821619926d-02)
     +                          *q + 5.28406320615584d-02
                           ww1 = (p+p)*f1 + dexp(-p)
                        else
                           q = p - 2.0d0
                           f1 = ((((((((((-1.61702782425558d-10*q+
     +                          1.96215250865776d-09)
     +                          *q-2.14234468198419d-08)
     +                          *q+2.17216556336318d-07)
     +                          *q-1.98850171329371d-06)
     +                          *q+1.62429321438911d-05)
     +                          *q-1.16740298039895d-04)
     +                          *q+7.24888732052332d-04)
     +                          *q-3.79490003707156d-03)
     +                          *q+1.61723488664661d-02)
     +                          *q-5.29428148329736d-02)
     +                          *q + 1.15702180856167d-01
                           ww1 = (p+p)*f1 + dexp(-p)
                        end if
                     else if (p.gt.3.0d-07) then
                        f1 = ((((((((-8.36313918003957d-08*p+
     +                       1.21222603512827d-06)
     +                       *p-1.15662609053481d-05)
     +                       *p+9.25197374512647d-05)
     +                       *p-6.40994113129432d-04)
     +                       *p+3.78787044215009d-03)
     +                       *p-1.85185172458485d-02)
     +                       *p+7.14285713298222d-02)
     +                       *p-1.99999999997023d-01)
     +                       *p + 3.33333333333318d-01
                        ww1 = (p+p)*f1 + dexp(-p)
                     else
                        ww1 = 1.0d0 - p/3.0d0
                     end if
                     sum = sum + dij(nn)*ww1*expe*dsqrt(ab)
                  end if
                  nn = nn + 16
 20            continue
               gout = gout + d2*sum
            end if
 30      continue
 40   continue
      g(1) = gout*pi252*qq4
      return
      end
      subroutine shells_3c(gout,iao_tag,icd_tag,
     &                     nelec,ish,jsh,ksh,lsh,iexch)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(../m4/common/sizes)
      parameter (mxp2 = mxprms * mxprms)
INCLUDE(../m4/common/infoa)
INCLUDE(common/dft_mbasis)
c
INCLUDE(common/dft_root)
      common/junkx/cxyz(3,5625),aaa(21*mxp2),ijaaa(225)
INCLUDE(common/dft_shlnos)
INCLUDE(common/dft_indez)
INCLUDE(common/dft_shlinf)
INCLUDE(common/dft_misc)
INCLUDE(common/dft_flips)
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
         ngtk = kgt(iexch)
         ngtl = lgt(iexch)
c
c     ----- kshell
c
         k = katom(icd_tag,ksh)
         pk = c(1,k)
         qk = c(2,k)
         rk = c(3,k)
         k1 = kstart(icd_tag,ksh)
         k2 = k1 + kng(icd_tag,ksh) - 1
         lkt = ktype(icd_tag,ksh)
         mink = kmin(icd_tag,ksh)
         maxk = kmax(icd_tag,ksh)
         lock = kloc(icd_tag,ksh) - mink
         ngc = 0
         do 20 k = k1 , k2
            ngc = ngc + 1
            cgg(ngc) = ex_m(2,k)
            csc(ngc) = cs(2,k)
            cpc(ngc) = cp(2,k)
            cdc(ngc) = cd(2,k)
            cfc(ngc) = cf(2,k)
            cgc(ngc) = cg(2,k)
 20      continue
c
c     ----- lshell
c
         l = k
         pl = c(1,l)
         ql = c(2,l)
         rl = c(3,l)
         l1 = 1
         l2 = 1
         llt = 1
         minl = 1
         maxl = 1
         locl = 0
         ngd = 1
         dg(ngd) = 0.0d0
         csd(ngd) = 1.0d0
         cpd(ngd) = 1.0d0
         cdd(ngd) = 1.0d0
         cfd(ngd) = 1.0d0
         cgd(ngd) = 1.0d0
 30      continue
         nroots = (lit+ljt+lkt+llt-4)/2 + 1
c        rrk = ((pk-pl)**2+(qk-ql)**2+(rk-rl)**2)
         rrk=0.0d0
c
c     ----- prepare indices for pairs of (k,l) functions
c
         okanl = .false.
         kl = 0
         lmax = maxl
         do 50 k = mink , maxk
           kl = kl + 1
           klx(kl) = kx(k)
           kly(kl) = ky(k)
           klz(kl) = kz(k)
           klgt(kl) = ngtk*(k-mink) 
 50      continue
         max = kl
         do 60 i = 1 , ij
c           if (oident) max = i
            ik(i) = max
 60      continue
         ijkl = ij*kl
c        if (oident) ijkl = ij*(ij+1)/2
c
c     zero integral storage
c
         do 65 i = 1,ij
         ngij = ijgt(i)
         do 65 k = 1,ik(i)
 65      gout(ngij + klgt(k)) = 0.0d0
         return
      else
c        oianj = ish.eq.jsh
         oianj = .false.
         ngti = igt(iexch)
         ngtj = jgt(iexch)
c
c     ----- ishell
c
         i = katom(iao_tag,ish)
         pi = c(1,i)
         qi = c(2,i)
         ri = c(3,i)
         i1 = kstart(iao_tag,ish)
         i2 = i1 + kng(iao_tag,ish) - 1
         lit = ktype(iao_tag,ish)
         mini = kmin(iao_tag,ish)
         maxi = kmax(iao_tag,ish)
         loci = kloc(iao_tag,ish) - mini
c        write(6,*) 'ISHELL',ish
c        write(6,*) 'mini',mini
c        write(6,*) 'maxi',maxi
c        write(6,*) 'loci',loci
         nga = 0
         do 70 i = i1 , i2
            nga = nga + 1
            ag(nga) = ex_m(iao_tag,i)
            csa(nga) = cs(iao_tag,i)
            cpa(nga) = cp(iao_tag,i)
            cda(nga) = cd(iao_tag,i)
            cfa(nga) = cf(iao_tag,i)
            cga(nga) = cg(iao_tag,i)
 70      continue
c
c     ----- jshell
c
         j = katom(iao_tag,jsh)
         pj = c(1,j)
         qj = c(2,j)
         rj = c(3,j)
         j1 = kstart(iao_tag,jsh)
         j2 = j1 + kng(iao_tag,jsh) - 1
         ljt = ktype(iao_tag,jsh)
         minj = kmin(iao_tag,jsh)
         maxj = kmax(iao_tag,jsh)
         locj = kloc(iao_tag,jsh) - minj
c        write(6,*) 'JSHELL',jsh
c        write(6,*) 'minj',minj
c        write(6,*) 'maxj',maxj
c        write(6,*) 'locj',locj
         ngb = 0
         do 80 j = j1 , j2
            ngb = ngb + 1
            bg(ngb) = ex_m(iao_tag,j)
            csb(ngb) = cs(iao_tag,j)
            cpb(ngb) = cp(iao_tag,j)
            cdb(ngb) = cd(iao_tag,j)
            cfb(ngb) = cf(iao_tag,j)
            cgb(ngb) = cg(iao_tag,j)
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
      subroutine spdint_3c(gout)
c
c     ----- form integrals over functions -----
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
INCLUDE(../m4/common/sizes)
      parameter (mxp2 = mxprms * mxprms)
INCLUDE(common/dft_root)
INCLUDE(common/dft_indez)
      common/junkx/pin(5625),qin(5625),rin(5625),aaa(21*mxp2),ijaaa(225)
INCLUDE(common/dft_shlnos)
INCLUDE(common/dft_denss)
      dimension gout(*)
c
      go to (20,50,80,110,140,170,200,230,260) , nroots
 20   do 40 i = 1 , ij
         d1 = dij(i)
         nx = ijx(i)
         ny = ijy(i)
         nz = ijz(i)
         n1 = ijgt(i)
         max = ik(i)
         do 30 k = 1 , max
            mx = nx + klx(k)
            my = ny + kly(k)
            mz = nz + klz(k)
            n = n1 + klgt(k)
            gout(n) = (pin(mx)*qin(my)*rin(mz))*d1*dkl(k) + gout(n)
 30      continue
 40   continue
      return
 50   do 70 i = 1 , ij
         d1 = dij(i)
         nx = ijx(i)
         ny = ijy(i)
         nz = ijz(i)
         n1 = ijgt(i)
         max = ik(i)
         do 60 k = 1 , max
            mx = nx + klx(k)
            my = ny + kly(k)
            mz = nz + klz(k)
            n = n1 + klgt(k)
            gout(n) = (pin(mx)*qin(my)*rin(mz)+pin(mx+625)*qin(my+625)
     +                *rin(mz+625))*d1*dkl(k) + gout(n)
 60      continue
 70   continue
      return
 80   do 100 i = 1 , ij
         d1 = dij(i)
         nx = ijx(i)
         ny = ijy(i)
         nz = ijz(i)
         n1 = ijgt(i)
         max = ik(i)
         do 90 k = 1 , max
            mx = nx + klx(k)
            my = ny + kly(k)
            mz = nz + klz(k)
            n = n1 + klgt(k)
            gout(n) = (pin(mx)*qin(my)*rin(mz)+
     +                 pin(mx+625)*qin(my+625)*rin(mz+625)+
     +                 pin(mx+1250)*qin(my+1250)*rin(mz+1250))
     +                *d1*dkl(k) + gout(n)
 90      continue
 100  continue
      return
 110  do 130 i = 1 , ij
         d1 = dij(i)
         nx = ijx(i)
         ny = ijy(i)
         nz = ijz(i)
         n1 = ijgt(i)
         max = ik(i)
         do 120 k = 1 , max
            mx = nx + klx(k)
            my = ny + kly(k)
            mz = nz + klz(k)
            n = n1 + klgt(k)
            gout(n) = (pin(mx)*qin(my)*rin(mz)+
     +                 pin(mx+625)*qin(my+625)*rin(mz+625)+
     +                 pin(mx+1250)*qin(my+1250)*rin(mz+1250)+
     +                 pin(mx+1875)*qin(my+1875)*rin(mz+1875))
     +                *d1*dkl(k) + gout(n)
 120     continue
 130  continue
      return
 140  do 160 i = 1 , ij
         d1 = dij(i)
         nx = ijx(i)
         ny = ijy(i)
         nz = ijz(i)
         n1 = ijgt(i)
         max = ik(i)
         do 150 k = 1 , max
            mx = nx + klx(k)
            my = ny + kly(k)
            mz = nz + klz(k)
            n = n1 + klgt(k)
            gout(n) = (pin(mx)*qin(my)*rin(mz)+
     +                 pin(mx+625)*qin(my+625)*rin(mz+625)+
     +                 pin(mx+1250)*qin(my+1250)*rin(mz+1250)+
     +                 pin(mx+1875)*qin(my+1875)*rin(mz+1875)+
     +                 pin(mx+2500)*qin(my+2500)*rin(mz+2500))
     +                *d1*dkl(k) + gout(n)
 150     continue
 160  continue
      return
 170  do 190 i = 1 , ij
         d1 = dij(i)
         nx = ijx(i)
         ny = ijy(i)
         nz = ijz(i)
         n1 = ijgt(i)
         max = ik(i)
         do 180 k = 1 , max
            mx = nx + klx(k)
            my = ny + kly(k)
            mz = nz + klz(k)
            n = n1 + klgt(k)
            gout(n) = (pin(mx)*qin(my)*rin(mz)+
     +                 pin(mx+625)*qin(my+625)*rin(mz+625)+
     +                 pin(mx+1250)*qin(my+1250)*rin(mz+1250)+
     +                 pin(mx+1875)*qin(my+1875)*rin(mz+1875)+
     +                 pin(mx+2500)*qin(my+2500)*rin(mz+2500)+
     +                 pin(mx+3125)*qin(my+3125)*rin(mz+3125))
     +                *d1*dkl(k) + gout(n)
 180     continue
 190  continue
      return
 200  do 220 i = 1 , ij
         d1 = dij(i)
         nx = ijx(i)
         ny = ijy(i)
         nz = ijz(i)
         n1 = ijgt(i)
         max = ik(i)
         do 210 k = 1 , max
            mx = nx + klx(k)
            my = ny + kly(k)
            mz = nz + klz(k)
            n = n1 + klgt(k)
            gout(n) = (pin(mx)*qin(my)*rin(mz)+
     +                 pin(mx+625)*qin(my+625)*rin(mz+625)+
     +                 pin(mx+1250)*qin(my+1250)*rin(mz+1250)+
     +                 pin(mx+1875)*qin(my+1875)*rin(mz+1875)+
     +                 pin(mx+2500)*qin(my+2500)*rin(mz+2500)+
     +                 pin(mx+3125)*qin(my+3125)*rin(mz+3125)+
     +                 pin(mx+3750)*qin(my+3750)*rin(mz+3750))
     +                *d1*dkl(k) + gout(n)
 210     continue
 220  continue
      return
 230  do 250 i = 1 , ij
         d1 = dij(i)
         nx = ijx(i)
         ny = ijy(i)
         nz = ijz(i)
         n1 = ijgt(i)
         max = ik(i)
         do 240 k = 1 , max
            mx = nx + klx(k)
            my = ny + kly(k)
            mz = nz + klz(k)
            n = n1 + klgt(k)
            gout(n) = (pin(mx)*qin(my)*rin(mz)+
     +                 pin(mx+625)*qin(my+625)*rin(mz+625)+
     +                 pin(mx+1250)*qin(my+1250)*rin(mz+1250)+
     +                 pin(mx+1875)*qin(my+1875)*rin(mz+1875)+
     +                 pin(mx+2500)*qin(my+2500)*rin(mz+2500)+
     +                 pin(mx+3125)*qin(my+3125)*rin(mz+3125)+
     +                 pin(mx+3750)*qin(my+3750)*rin(mz+3750)+
     +                 pin(mx+4375)*qin(my+4375)*rin(mz+4375))
     +                *d1*dkl(k)+gout(n)
 240     continue
 250  continue
      return
c
 260  do 280 i = 1,ij
         d1 = dij(i)
         nx = ijx(i)
         ny = ijy(i)
         nz = ijz(i)
         n1 = ijgt(i)
         max = ik(i)
         do 270 k = 1,max
            mx = nx+klx(k)
            my = ny+kly(k)
            mz = nz+klz(k)
            n = n1+klgt(k)
            gout(n) = (pin(mx)*qin(my)*rin(mz)+
     +                 pin(mx+625)*qin(my+625)*rin(mz+625)+
     +                 pin(mx+1250)*qin(my+1250)*rin(mz+1250)+
     +                 pin(mx+1875)*qin(my+1875)*rin(mz+1875)+
     +                 pin(mx+2500)*qin(my+2500)*rin(mz+2500)+
     +                 pin(mx+3125)*qin(my+3125)*rin(mz+3125)+
     +                 pin(mx+3750)*qin(my+3750)*rin(mz+3750)+
     +                 pin(mx+4375)*qin(my+4375)*rin(mz+4375)+
     +                 pin(mx+5000)*qin(my+5000)*rin(mz+5000))
     +                *d1*dkl(k)+gout(n)
 270     continue
 280  continue
      return
      end
      subroutine sskl_3c(gout)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(../m4/common/sizes)
      parameter (mxp2 = mxprms * mxprms)
INCLUDE(../m4/common/restar)
c
INCLUDE(common/dft_shlinf)
      common/junkx/pint(5625),qint(5625),rint(5625),
     + aa(mxp2),r(mxp2),p1(mxp2),q1(mxp2),r1(mxp2),dd(16*mxp2),
     + ijd(225)
INCLUDE(common/dft_indez)
INCLUDE(common/dft_shlnos)
INCLUDE(common/dft_root)
INCLUDE(common/dft_shlt)
INCLUDE(common/dft_misc)
INCLUDE(common/dft_setint)
INCLUDE(common/dft_denss)
      dimension gout(*)
      data kln1,kln2 /5,1/
      data sqrt3 /1.73205080756888d0/
      data sqrt5 /2.23606797749979d0/
      data sqrt7 /2.64575131106459d0/
      data pi252 /34.986836655250d0/
      data dzero,pt5,done /0.0d0,0.5d0,1.0d0/
      factor = pi252*qq4
      onorm = normf.ne.1 .or. normp.ne.1
      onorm = .false.
      if (lkt.lt.llt) then
         nk = llt - 1
         nl = lkt - 1
         kl1 = kln2
         kl2 = kln1
         pd = pl
         qd = ql
         rd = rl
         dxkl = pl - pk
         dykl = ql - qk
         dzkl = rl - rk
      else
c
c     ----- select expansion centre for -xyz- integrals -----
c
         nk = lkt - 1
         nl = llt - 1
         kl1 = kln1
         kl2 = kln2
         pd = pk
         qd = qk
         rd = rk
         dxkl = pk - pl
         dykl = qk - ql
         dzkl = rk - rl
      end if
      mmax = nk + nl
      max = mmax + 1
      do 20 k = 1 , max
         n = k - 1
         if (n.le.nk) kn(k) = kl1*n
         if (n.gt.nk) kn(k) = kl1*nk + kl2*(n-nk)
 20   continue
      lgmax = ngd
c
c     ----- k primitive
c
      do 370 kg = 1 , ngc
         ak = cgg(kg)
         brrk = ak*rrk
         akxk = ak*pk
         akyk = ak*qk
         akzk = ak*rk
         csk = csc(kg)*factor
         cpk = cpc(kg)*factor
         cdk = cdc(kg)*factor
         cfk = cfc(kg)*factor
         cgk = cgc(kg)*factor
         if (okanl) lgmax = kg
c
c     ----- l primitive
c
         do 360 lg = 1 , lgmax
            al = dg(lg)
            b = ak + al
            b1 = done/b
            bbrrk = al*brrk*b1
            if (bbrrk.le.tol) then
               csl = csd(lg)
               cpl = cpd(lg)
               cdl = cdd(lg)
               cfl = cfd(lg)
               cgl = cgd(lg)
               pb = (akxk+al*pl)*b1
               qb = (akyk+al*ql)*b1
               rb = (akzk+al*rl)*b1
               bxbd = b*(pb-pd)
               bybd = b*(qb-qd)
               bzbd = b*(rb-rd)
               odoub = okanl .and. kg.gt.lg
c
c     ----- density factor
c
               n = 0
               max = maxl
               do 200 k = mink , maxk
                  go to (30,40,100,100,
     +                   50,100,100,60,100,100,
     +                   70,100,100,80,100,100,100,100,100,90,
     +                   92,100,100,94,100,100,100,100,100,96,
     +                  100,100,98,100,100), k
c
 30               dum1 = csk*b1
                  go to 100
 40               dum1 = cpk*b1
                  go to 100
 50               dum1 = cdk*b1
                  go to 100
 60               if (onorm) dum1 = dum1*sqrt3
                  go to 100
 70               dum1 = cfk*b1
                  go to 100
 80               if (onorm) dum1 = dum1*sqrt5
                  go to 100
 90               if (onorm) dum1 = dum1*sqrt3
                  go to 100
 92               dum1 = cgk*b1
                  go to 100
 94               if (onorm) dum1 = dum1*sqrt7
                  go to 100
 96               if (onorm) dum1 = dum1*sqrt5/sqrt3
                  go to 100
 98               if (onorm) dum1 = dum1*sqrt3
 100              if (okanl) max = k
                  do 190 l = minl , max
                     go to (110,120,180,180,
     +                      130,180,180,140,180,180,
     +                      150,180,180,160,180,180,180,180,180,170,
     +                      182,180,180,184,180,180,180,180,180,186,
     +                      180,180,188,180,180), l
 110                 dum2 = dum1*csl
                     if (odoub) then
                        if (k.gt.1) then
                           dum2 = dum2 + csk*cpl*b1
                        else
                           dum2 = dum2 + dum2
                        end if
                     end if
                     go to 180
 120                 dum2 = dum1*cpl
                     if (odoub) dum2 = dum2 + dum2
                     go to 180
 130                 dum2 = dum1*cdl
                     if (odoub) dum2 = dum2 + dum2
                     go to 180
 140                 if (onorm) dum2 = dum2*sqrt3
                     go to 180
 150                 dum2 = dum1*cfl
                     if (odoub) dum2 = dum2 + dum2
                     go to 180
 160                 if (onorm) dum2 = dum2*sqrt5
                     go to 180
 170                 if (onorm) dum2 = dum2*sqrt3
                     go to 180
 182                 dum2 = dum1*cgl
                     if (odoub) dum2 = dum2+dum2
                     go to 180
 184                 if (onorm) dum2 = dum2 * sqrt7
                     go to 180
 186                 if (onorm) dum2 = dum2 * sqrt5/sqrt3
                     go to 180
 188                 if (onorm) dum2 = dum2 * sqrt3
 180                 n = n + 1
                     dkl(n) = dum2
 190              continue
 200           continue
               nn = 0
c
c     ----- pair of i,j primitives
c
               do 350 n = 1 , nij
                  dum = bbrrk + r(n)
                  if (dum.gt.tol) then
                     nn = nn + 16
                     go to 350
                  else
                     a = aa(n)
                     ab = a*b
                     aandb = a + b
                     expe = dd(1+nn)*dexp(-dum)/dsqrt(aandb)
                     rho = ab/aandb
                     pa = p1(n)
                     qa = q1(n)
                     ra = r1(n)
                     pp = rho*((pa-pb)**2+(qa-qb)**2+(ra-rb)**2)
                     axad = a*(pa-pd)
                     ayad = a*(qa-qd)
                     azad = a*(ra-rd)
                     c1x = bxbd + axad
                     c2x = a*bxbd
                     c1y = bybd + ayad
                     c2y = a*bybd
                     c1z = bzbd + azad
                     c2z = a*bzbd
                     if (nroots.le.3) then
                        call rt123_dft
                     else if (nroots.eq.4) then
                        call roots4_dft
                     else
                        call roots5_dft
                     end if
c
c     ----- roots and weights for quadrature
c
                     mm = 0
                     do 260 m = 1 , nroots
c
c     compute two-electron  integrals for each root
c
                        u2 = u(m)*rho
                        f00 = expe*w(m)
                        dum = done/(ab+u2*aandb)
                        bp01 = (a+u2)*pt5*dum
                        pcp00 = (u2*c1x+c2x)*dum
                        qcp00 = (u2*c1y+c2y)*dum
                        rcp00 = (u2*c1z+c2z)*dum
c
c     ----- i(0,0) -----
c
                        i1 = 1 + mm
                        pint(i1) = done
                        qint(i1) = done
                        rint(i1) = f00
                        if (mmax.eq.0) then
                           mm = mm + 625
                           go to 260
                        else
c
c     ----- i(0,1) -----
c
                           k2 = kn(2)
                           i3 = i1 + k2
                           pint(i3) = pcp00
                           qint(i3) = qcp00
                           rint(i3) = rcp00*f00
                           if (mmax.le.1) then
                              mm = mm + 625
                              go to 260
                           else
c
c     ----- i(0,kkk) -----
c
                              cp01 = dzero
                              i3 = i1
                              i4 = i1 + k2
                              do 210 kkk = 2 , mmax
                                 cp01 = cp01 + bp01
                                 i5 = i1 + kn(kkk+1)
                                 pint(i5) = cp01*pint(i3)
     +                              + pcp00*pint(i4)
                                 qint(i5) = cp01*qint(i3)
     +                              + qcp00*qint(i4)
                                 rint(i5) = cp01*rint(i3)
     +                              + rcp00*rint(i4)
                                 i3 = i4
                                 i4 = i5
 210                          continue
                              if (nl.eq.0) then
                                 mm = mm + 625
                                 go to 260
                              else
c
c     ----- i(0,0,kkk,lll) -----
c
                                 i5 = kn(mmax+1)
                                 min = nk
                              end if
                           end if
                        end if
 220                    kkk = mmax
                        i3 = i1 + i5
 230                    i4 = i1 + kn(kkk)
                        pint(i3) = pint(i3) + dxkl*pint(i4)
                        qint(i3) = qint(i3) + dykl*qint(i4)
                        rint(i3) = rint(i3) + dzkl*rint(i4)
                        i3 = i4
                        kkk = kkk - 1
                        if (kkk.gt.min) go to 230
                        min = min + 1
                        if (min.lt.mmax) go to 220
                        if (nk.ne.0) then
                           i3 = i1 + kl2
                           do 250 lll = 1 , nl
                              i4 = i3
                              do 240 kkk = 1 , nk
                                 pint(i4) = pint(i4+kl1-kl2)
     +                              + dxkl*pint(i4-kl2)
                                 qint(i4) = qint(i4+kl1-kl2)
     +                              + dykl*qint(i4-kl2)
                                 rint(i4) = rint(i4+kl1-kl2)
     +                              + dzkl*rint(i4-kl2)
                                 i4 = i4 + kl1
 240                          continue
                              i3 = i3 + kl2
 250                       continue
                        end if
                        mm = mm + 625
 260                 continue
c
c     ----- form (i,j//k,l) integrals over functions
c
                     go to (270,290,310,330,335) , nroots
                  end if
 270              do 280 kkk = 1 , kl
                     mx = 1 + klx(kkk)
                     my = 1 + kly(kkk)
                     mz = 1 + klz(kkk)
                     nnn = 1 + klgt(kkk)
                     gout(nnn) = (pint(mx)*qint(my)*rint(mz))
     +                         *dkl(kkk) + gout(nnn)
 280              continue
                  nn = nn + 16
                  go to 350
 290              do 300 kkk = 1 , kl
                     mx = 1 + klx(kkk)
                     my = 1 + kly(kkk)
                     mz = 1 + klz(kkk)
                     nnn = 1 + klgt(kkk)
                     gout(nnn) = (pint(mx)*qint(my)*rint(mz)+
     +                       pint(mx+625)*qint(my+625)*rint(mz+625))
     +                       *dkl(kkk) + gout(nnn)
 300              continue
                  nn = nn + 16
                  go to 350
 310              do 320 kkk = 1 , kl
                     mx = 1 + klx(kkk)
                     my = 1 + kly(kkk)
                     mz = 1 + klz(kkk)
                     nnn = 1 + klgt(kkk)
                     gout(nnn) = (pint(mx)*qint(my)*rint(mz)+
     +                       pint(mx+625)*qint(my+625)*rint(mz+625)+
     +                       pint(mx+1250)*qint(my+1250)*rint(mz+1250))
     +                       *dkl(kkk) + gout(nnn)
 320              continue
                  nn = nn + 16
                  go to 350
 330              do 340 kkk = 1 , kl
                     mx = 1 + klx(kkk)
                     my = 1 + kly(kkk)
                     mz = 1 + klz(kkk)
                     nnn = 1 + klgt(kkk)
                     gout(nnn) = (pint(mx)*qint(my)*rint(mz)+
     +                       pint(mx+625)*qint(my+625)*rint(mz+625)+
     +                       pint(mx+1250)*qint(my+1250)*rint(mz+1250)+
     +                       pint(mx+1875)*qint(my+1875)*rint(mz+1875))
     +                         *dkl(kkk) + gout(nnn)
 340              continue
                  nn = nn + 16
                  go to 350
 335              do 345 kkk = 1 , kl
                     mx = 1 + klx(kkk)
                     my = 1 + kly(kkk)
                     mz = 1 + klz(kkk)
                     nnn = 1 + klgt(kkk)
                     gout(nnn) = (pint(mx)*qint(my)*rint(mz)+
     +                       pint(mx+625)*qint(my+625)*rint(mz+625)+
     +                       pint(mx+1250)*qint(my+1250)*rint(mz+1250)+
     +                       pint(mx+1875)*qint(my+1875)*rint(mz+1875)+
     +                       pint(mx+2500)*qint(my+2500)*rint(mz+2500))
     +                         *dkl(kkk) + gout(nnn)
 345              continue
                  nn = nn + 16
 350           continue
            end if
 360     continue
 370  continue
      return
      end
      subroutine ssprm_3c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(../m4/common/sizes)
      parameter (mxp2 = mxprms * mxprms)
INCLUDE(common/dft_shlt)
INCLUDE(common/dft_shlinf)
INCLUDE(common/dft_misc)
INCLUDE(common/dft_shlnos)
      common/junkx/cxyz(3,5625),
     *a(mxp2),r(mxp2),p1(mxp2),q1(mxp2),r1(mxp2),dij(16*mxp2),
     + ijd(225)
      data done /1.0d0/
c
c     ----- i primitive
c
      nij = 0
      jbmax = ngb

      do 30 ia = 1 , nga
         ai = ag(ia)
         arri = ai*rri
         axi = ai*pi
         ayi = ai*qi
         azi = ai*ri
         csi = csa(ia)
         if (oianj) jbmax = ia
c
c     ----- j primitive
c
         do 20 jb = 1 , jbmax
            aj = bg(jb)
            aa = ai + aj
            aa1 = done/aa
            dum = aj*arri*aa1
            if (dum.le.tol) then
               nn = 1 + nij*16
               nij = nij + 1
               r(nij) = dum
               a(nij) = aa
               p1(nij) = (axi+aj*pj)*aa1
               q1(nij) = (ayi+aj*qj)*aa1
               r1(nij) = (azi+aj*rj)*aa1
               dij(nn) = csi*csb(jb)*aa1
               if (oianj .and. (ia.gt.jb)) dij(nn) = dij(nn) + dij(nn)
            end if
 20      continue
 30   continue
      return
      end
      subroutine xyzint_3c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
INCLUDE(common/dft_setint)
      common/junkx/pint(5625),qint(5625),rint(5625)
c
      data dzero,done /0.0d0,1.0d0/
c
      on0 = nmax.eq.0
      on1 = nmax.le.1
      om0 = mmax.eq.0
      om1 = mmax.le.1
c
c     ----- in(0,0) -----
c
      i1 = in(1)
      pint(i1) = done
      qint(i1) = done
      rint(i1) = f00
      if (on0 .and. om0) return
      if (.not.(on0)) then
c
c     ----- in(1,0) -----
c
         i2 = in(2)
         pint(i2) = pc00
         qint(i2) = qc00
         rint(i2) = rc00*f00
      end if
      if (.not.(om0)) then
c
c     ----- in(0,1) -----
c
         k2 = kn(2)
         i3 = i1 + k2
         pint(i3) = pcp00
         qint(i3) = qcp00
         rint(i3) = rcp00*f00
         if (.not.(on0)) then
c
c     ----- in(1,1) -----
c
            i3 = i2 + k2
            cp10 = b00
            pint(i3) = pcp00*pint(i2) + cp10
            qint(i3) = qcp00*qint(i2) + cp10
            rint(i3) = rcp00*rint(i2) + cp10*f00
         end if
      end if
      if (.not.(on1)) then
         c10 = dzero
         i3 = i1
         i4 = i2
         do 20 n = 2 , nmax
            c10 = c10 + b10
c
c     ----- in(n,0) -----
c
            i5 = in(n+1)
            pint(i5) = c10*pint(i3) + pc00*pint(i4)
            qint(i5) = c10*qint(i3) + qc00*qint(i4)
            rint(i5) = c10*rint(i3) + rc00*rint(i4)
            if (.not.(om0)) then
               cp10 = cp10 + b00
c
c     ----- in(n,1) -----
c
               i3 = i5 + k2
               pint(i3) = pcp00*pint(i5) + cp10*pint(i4)
               qint(i3) = qcp00*qint(i5) + cp10*qint(i4)
               rint(i3) = rcp00*rint(i5) + cp10*rint(i4)
            end if
            i3 = i4
            i4 = i5
 20      continue
      end if
      if (.not.(om1)) then
         cp01 = dzero
         c01 = b00
         i3 = i1
         i4 = i1 + k2
         do 30 m = 2 , mmax
            cp01 = cp01 + bp01
c
c     ----- in(0,m) -----
c
            i5 = i1 + kn(m+1)
            pint(i5) = cp01*pint(i3) + pcp00*pint(i4)
            qint(i5) = cp01*qint(i3) + qcp00*qint(i4)
            rint(i5) = cp01*rint(i3) + rcp00*rint(i4)
            if (.not.(on0)) then
               c01 = c01 + b00
c
c     ----- in(1,m) -----
c
               i3 = i2 + kn(m+1)
               pint(i3) = pc00*pint(i5) + c01*pint(i4)
               qint(i3) = qc00*qint(i5) + c01*qint(i4)
               rint(i3) = rc00*rint(i5) + c01*rint(i4)
            end if
            i3 = i4
            i4 = i5
 30      continue
      end if
      if (.not.(on1 .or. om1)) then
c
c     ----- in(n,m) -----
c
         c01 = b00
         k3 = k2
         do 50 m = 2 , mmax
            k4 = kn(m+1)
            c01 = c01 + b00
            i3 = i1
            i4 = i2
            c10 = b10
            do 40 n = 2 , nmax
               i5 = in(n+1)
               pint(i5+k4) = c10*pint(i3+k4) + pc00*pint(i4+k4)
     +                       + c01*pint(i4+k3)
               qint(i5+k4) = c10*qint(i3+k4) + qc00*qint(i4+k4)
     +                       + c01*qint(i4+k3)
               rint(i5+k4) = c10*rint(i3+k4) + rc00*rint(i4+k4)
     +                       + c01*rint(i4+k3)
               c10 = c10 + b10
               i3 = i4
               i4 = i5
 40         continue
            k3 = k4
 50      continue
      end if
      if (nj.eq.0) go to 110
c
c     ----- in(mi,mj,m) -----
c
      m = 0
      i5 = in(nmax+1)
 60   min = ni
      km = kn(m+1)
 70   n = nmax
      i3 = i5 + km
 80   i4 = in(n) + km
      pint(i3) = pint(i3) + dxij*pint(i4)
      qint(i3) = qint(i3) + dyij*qint(i4)
      rint(i3) = rint(i3) + dzij*rint(i4)
      i3 = i4
      n = n - 1
      if (n.gt.min) go to 80
      min = min + 1
      if (min.lt.nmax) go to 70
      if (ni.ne.0) then
         i3 = ij2 + i1 + km
         do 100 mj = 1 , nj
            i4 = i3
            do 90 mi = 1 , ni
               pint(i4) = pint(i4+ij1-ij2) + dxij*pint(i4-ij2)
               qint(i4) = qint(i4+ij1-ij2) + dyij*qint(i4-ij2)
               rint(i4) = rint(i4+ij1-ij2) + dzij*rint(i4-ij2)
               i4 = i4 + ij1
 90         continue
            i3 = i3 + ij2
 100     continue
      end if
      m = m + 1
      if (m.le.mmax) go to 60
 110  if (nl.eq.0) go to 170
c
c     ----- in(mi,mj,mk,ml) -----
c
      i5 = kn(mmax+1)
      ia = i1
      mi = 0
 120  mj = 0
      ib = ia
      min = nk
 130  m = mmax
      i3 = ib + i5
 140  i4 = ib + kn(m)
      pint(i3) = pint(i3) + dxkl*pint(i4)
      qint(i3) = qint(i3) + dykl*qint(i4)
      rint(i3) = rint(i3) + dzkl*rint(i4)
      i3 = i4
      m = m - 1
      if (m.gt.min) go to 140
      min = min + 1
      if (min.lt.mmax) go to 130
      if (nk.ne.0) then
         i3 = ib + kl2
         do 160 ml = 1 , nl
            i4 = i3
            do 150 mk = 1 , nk
               pint(i4) = pint(i4+kl1-kl2) + dxkl*pint(i4-kl2)
               qint(i4) = qint(i4+kl1-kl2) + dykl*qint(i4-kl2)
               rint(i4) = rint(i4+kl1-kl2) + dzkl*rint(i4-kl2)
               i4 = i4 + kl1
 150        continue
            i3 = i3 + kl2
 160     continue
      end if
      mj = mj + 1
      ib = ib + ij2
      if (mj.le.nj) then
         min = nk
         go to 130
      else
         mi = mi + 1
         ia = ia + ij1
         if (mi.le.ni) go to 120
      end if
 170  return
      end
      subroutine mat_form_3c(gout,mini,maxi,minj,maxj,mink,maxk,
     &                       ao_basfn,cd_basfn,imc,
     &                       eri_scr)
C *****************************************************************
C *Description:                                                   *
C *Form two centre repulsion integral matrix using array pointers *
C *****************************************************************
      implicit none
C *****************************************************************
C *Declarations                                                   *
C *                                                              *
C *In variables                                                  *
      REAL gout(*)
      integer mini,maxi,minj,maxj,mink,maxk
      integer ao_basfn,cd_basfn
INCLUDE(common/dft_indez)
C *Out variables                                                 *
      REAL eri_scr(225,*)
      integer imc
C *Local variables                                               *
      integer ltyi,ltyj,ltyk
      integer ijn,knn,nn,ann,imm
C *End declarations                                              *
C ****************************************************************
      ijn=0
      ann=0
      do ltyi=mini,maxi
        do ltyj=minj,maxj
          ijn=ijn+1
          knn=0
          ann=ann+1
          imm=imc
          do ltyk=mink,maxk
            knn=knn+1
            imm=imm+1
            nn=ijgt(ijn)+klgt(knn)
            eri_scr(ann,imm)=gout(nn)
          enddo
        enddo
      enddo
      return
      end
      subroutine mat_save_3c(gout,mini,maxi,minj,maxj,mink,maxk,
     &                       te3c_int,ite3c_int)
C *****************************************************************
C *Description:                                                   *
C *Form two centre repulsion integral matrix using array pointers *
C *****************************************************************
      implicit none
C *****************************************************************
C *Declarations                                                   *
C *                                                              *
C *In variables                                                  *
      REAL gout(*)
      integer mini,maxi,minj,maxj,mink,maxk
INCLUDE(common/dft_indez)
C *Out variables                                                 *
      integer ite3c_int
      REAL te3c_int(*)
C *Local variables                                               *
      integer ltyi,ltyj,ltyk
      integer ijn,knn,nn,ann,imm
C *End declarations                                              *
C ****************************************************************
      ijn=0
      ann=0
      imm=ite3c_int
      do ltyi=mini,maxi
        do ltyj=minj,maxj
          ijn=ijn+1
          knn=0
          do ltyk=mink,maxk
            knn=knn+1
            imm=imm+1
            nn=ijgt(ijn)+klgt(knn)
            te3c_int(imm)=gout(nn)
          enddo
        enddo
      enddo
      ite3c_int = imm
      return
      end
      subroutine dvec_fill(locij,adens,
     &                     ii,jj,mini,maxi,minj,maxj,
     &                     dvec_scr)
C ***************************************************************
      implicit none
C ***************************************************************
C *Declarations				      *
      REAL adens(*)
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/nshel)
      integer locij
      integer mini,maxi,minj,maxj,ii,jj
      integer ibas_num,jbas_num
      REAL dvec_scr(*)
      integer lbai,lbaj,ijbas_num,inum,jnum
      integer dpos,ipos,jpos,maxa,mina
C ***************************************************************
      dpos=0
      ibas_num=0
      jbas_num=0
      ipos=kloc(ii)
      jpos=kloc(jj)
      inum=(kmax(ii)-kmin(ii))
      jnum=(kmax(jj)-kmin(jj))
      do lbai=0,inum
        ibas_num=ipos+(lbai)
        do lbaj=0,jnum
          dpos=dpos+1
          jbas_num=jpos+(lbaj)
          maxa=max(ibas_num,jbas_num)
          mina=min(ibas_num,jbas_num)
          ijbas_num=(((maxa*maxa)-maxa)/2)+mina
          dvec_scr(dpos)=adens(ijbas_num)
c         write(6,*) 'Locations:',ipos,jpos,ibas_num,jbas_num
c         write(6,*) 'DVEC:',ipos,jpos,ijbas_num,dvec_scr(dpos)
        enddo
      enddo
      locij=locij+ibas_num
      return
      end 
      subroutine km_fill(locij,ii,jj,kma,dvec_scr)
C *************************************************************
      implicit none
C *************************************************************
C *Declarations                                               *
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/nshel)
      REAL dvec_scr(*)
      integer ii,jj,locij
      integer ibas_num,jbas_num
      REAL kma(*)
      integer lbai,lbaj,inum,jnum
      integer dpos,ipos,jpos,maxa,mina,ijbas_num
C *************************************************************
      dpos=0
      ibas_num=0
      jbas_num=0
      ipos=kloc(ii)
      jpos=kloc(jj)
      inum=kmax(ii)-kmin(ii)
      jnum=kmax(jj)-kmin(jj)
      do lbai=0,inum
        ibas_num=ipos+lbai
        do lbaj=0,jnum
          dpos=dpos+1
          jbas_num=jpos+lbaj
          maxa=max(ibas_num,jbas_num)
          mina=min(ibas_num,jbas_num)
          ijbas_num=(((maxa*maxa)-maxa)/2)+mina
          kma(ijbas_num)=dvec_scr(dpos)
        enddo
      enddo
      locij=dpos
      return
      end 

      subroutine aux_find(tag)
      implicit none
INCLUDE(../m4/common/sizes)
INCLUDE(common/dft_mbasis)
INCLUDE(common/dft_auxg)
      integer tag
      REAL exad_find
      integer lshl,lprm
      integer pfirst,qfirst
      REAL ex1,ex2,ex_min 
      do lshl=1,nshell(tag)
        pfirst=kstart(tag,lshl)
        ex1=ex_m(tag,pfirst)
        qfirst=pfirst
        do lprm=1,kng(tag,lshl)
          ex2=ex_m(tag,qfirst)
          ex_min=min(ex1,ex2)
          qfirst=qfirst+1
        enddo
        auxg(lshl)=ex_min
      enddo
      return
      end 
      REAL function ovtest_3c(shl1,shl2,rsquar)
      implicit none
INCLUDE(common/dft_auxg)
      integer shl1,shl2
      REAL rsquar
      REAL exad1,exad2
      REAL abi,p

      exad1=auxg(shl1)
      exad2=auxg(shl2)
    
      abi    = 1.0d0/(exad1+exad2)
      p      = exad1*exad2*abi
      ovtest_3c = log(p*abi*4.0d0)*0.75d0
      ovtest_3c = ovtest_3c-rsquar*p
      return
      end     
