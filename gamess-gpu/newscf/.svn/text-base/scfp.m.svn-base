c
c temporary SCF driver for test SCF routine
c
      subroutine scfp(q)
c
      use newscf_routines
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/cslosc)
INCLUDE(../m4/common/field)
INCLUDE(../m4/common/scrf)
INCLUDE(../m4/common/gvalue)
INCLUDE(../m4/common/tran)
INCLUDE(../m4/common/saveco)
INCLUDE(../m4/common/harmon)
INCLUDE(../m4/common/runlab)
INCLUDE(../m4/common/scra7)
INCLUDE(../m4/common/symtry)
INCLUDE(../m4/common/scfopt)
INCLUDE(../m4/common/iofile)
INCLUDE(../m4/common/dump3)
INCLUDE(../m4/common/funct)
INCLUDE(../m4/common/infoa)
INCLUDE(../m4/common/restar)
INCLUDE(../m4/common/restrj)
INCLUDE(../m4/common/filel)
INCLUDE(../m4/common/phycon)
      common/diisd/ sta(210),cca(20),ra(20),scalea(20),iposa(20),
     + nstora,mpa,odiisa,
     + nspaca,stb(210),ccb(20),rb(20),scaleb(20),iposb(20),nstorb,
     + mpb,odiisb,nspacb,nsss(30),igvbo(maxorb),igvb1(maxorb),igsp(4),
     + ekk(63),intci(150)
      common/bufb/ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720)
      common/scra/iso(mxshel,48)
INCLUDE(../m4/common/zorac)
INCLUDE(../m4/common/statis)
INCLUDE(../m4/common/timeperiods)
      common/craypk/ijkl(1360)
INCLUDE(../m4/common/segm)
INCLUDE(../m4/common/dnfnw)
INCLUDE(../m4/common/xfield)
      common/blkin/corev(512),array(10)
_IF(parallel)
INCLUDE(../m4/common/nodeio)
INCLUDE(../m4/common/parcntl)
_ENDIF

INCLUDE(../m4/common/newscf)
      logical mpassi
      common/integmp/mpassi

      common/restri/jjfile(63),lds(508),isect(508),ldsect(508)
      character *8 title,guess
      common/restrz/title(12),guess
      dimension zcas(2),q(*),o1e(6)
      integer loop

c
      data zrhf /'rhf'/,zcas   /'casscf','mcscf'/
      data zuhf /'uhf'/
      data m1,m5,m16/1,5,16/
      data zgvb,zgrhf /'gvb','grhf'/
      data dzero /0.0d0/
c
      nav = lenwrd()
      if (nt.gt.1) then
       call rdedx(ptr,nw196(1),ibl196(1),idaf)
       call rdedx(dtr,nw196(2),ibl196(2),idaf)
       call rdedx(ftr,nw196(3),ibl196(3),idaf)
       call rdedx(gtr,nw196(4),ibl196(4),idaf)
       call readi(iso,nw196(5)*nav,ibl196(5),idaf)
      endif
      call cpuwal(begin,ebegin)
      call start_time_period(TP_SCF)
      if(nprint.ne.-5)write(iwr,180)
 180  format(//1x,129('-'))
c
c     ----- set convergence criteria -----
c
      iextin = 4
      exttol = 1.0d-3
      dmptol = 1.0d-4
      vshtol = 0.4d0
      if(nconv.le.0)nconv=5
      acurcy=10.0d0**(-nconv)
      nstora= 0
      nstorb= 0
      accdi2=acurcy*acurcy
c
c     clear buffer for integral output
c
      call setsto(1360,0,ijkl)
      odiisa=.false.
      odiisb=.false.
      if(zscftp.eq.zcas(1).or.zscftp.eq.zcas(2)  )go to 130

c$$$      if (nt.gt.1 .or. idadap .ne. 1) then
c$$$       call caserr('disable symmetry and adaption to use newscf')
c$$$      endif

c$$$      if (idadap .ne. 1) then
c$$$       call caserr('disable symmetry and adaption to use newscf')
c$$$      endif


c$$$      if (oharm) then
c$$$       call caserr('disable harmonic to use newscf')
c$$$      endif

c
c     ----- read in transformation matrices for s,p,d,f basis functions.
c
      if (nt.gt.1) then
       call rdedx(ptr,nw196(1),ibl196(1),idaf)
       call rdedx(dtr,nw196(2),ibl196(2),idaf)
       call rdedx(ftr,nw196(3),ibl196(3),idaf)
       call rdedx(gtr,nw196(4),ibl196(4),idaf)
       call readi(iso,nw196(5)*nav,ibl196(5),idaf)
      endif
c     minimize triangle usage?
      otriang = mpassi.and. .not.ofield.and. .not.oscrf
      l2=num*(num+1)/2
      if(otriang) then
       ilen = l2 + l2
       i10 = igmem_alloc(ilen)
       i20 = i10 + l2
       last = i20 + l2
      else
       len2=l2+l2
       ilen = 3*l2
       if (ofield) then
         ilen = ilen + 3*l2
       endif
       i10 = igmem_alloc(ilen)      
       i20=i10+l2
       i30=i20+l2
       i40=i30+l2
       if (ofield) then
        i50  = i40 + l2
        i60  = i50 + l2
        last = i60 + l2
       else
        last=i40
       endif
      endif
      length=last-i10

      if(length .ne. ilen)then
         write(6,*)length,ilen
         call caserr('mem size error')
      endif
c
c ----- if scrf then restore dipole integrals
c       and write to ed7                       ----
c
      if(oscrf) call dipmat(zscftp,q)
c
c ----- first restore s,t,f from section 192 and store on ed7
c
c ----- are field values to be introduced
c
      if(ofield) then
c
       do loop =1,6
        o1e(loop) = .true.
       enddo
       call getmat(q(i10),q(i20 ),q(i30),q(i40),q(i50),q(i60),
     * array,num,o1e,ionsec)
c ----- field calculation
       do 220 loop=1,l2
       q(i30+loop-1) = q(i30+loop-1) - fieldx * q(i40+loop-1)
     *                               - fieldy * q(i50+loop-1)
     *                               - fieldz * q(i60+loop-1)
 220   continue
       if(nprint.ne.-5) then
        write(iwr,230)fieldx,fieldy,fieldz
 230    format(/5x,38('*')/
     *          5x,'one-electron hamiltonian modified with'/
     *          5x,'x *',f11.7/
     *          5x,'y *',f11.7/
     *          5x,'z *',f11.7/
     *          5x,38('*') /)
       endif
      else
c ------ scrf calculation
       if(oscrf) then
c ------ convert angstrom to a.u.in g aradix,y,z
        aradix=aradix/toang(1)
        aradiy=aradiy/toang(1)
        aradiz=aradiz/toang(1)
c ------ calculate g factor (if sphere then only need use x)
        gx = (2.0d0*(dielec-1.0d0))/
     +                ((2.0d0*dielec + 1.0d0)*(aradix**3))
c ----- check if charged, and alter g accordingly
c ----- +ve ions
        if(ich.gt.0)then
                gx = gx*ne
                gx = gx/(ne+ich)
        endif
c ----- -ve ions
        if(ich.lt.0)then
                gx = gx*(ne+ich)
                gx = gx/ne
        endif
c
        gy = gx
        gz = gx
c
        if(nprint.ne.-5) then
        write(iwr,*)'    **********************************************'
        write(iwr,231)gx
        write(iwr,232)aradix
        write(iwr,*)'    **********************************************'
 231    format(/5x,'*** f factor set at            ',f11.7,' ***')
 232    format(/5x,'*** solvent cavity set at      ',f7.3,'au   ***',/)
        endif
       endif
c
_IF(parallel)
      if(ipiomode .eq. IO_NZ_S) then
       oswed3(4) = .false.
       oswed3(8) = .false.
      endif

      omaster = .true.
      if(ipiomode .eq. IO_NZ_S .and. .not. opg_root())
     &     omaster = .false.
       if(omaster) then
_ENDIF
        if(otriang) then
c
c      do this in 3 passes to reduce overall memory usage
c      to two triangles (including that for tranp)
c
         do mpass= 1,3
          do loop = 1,6
            o1e(loop) = .false.
          enddo
          o1e(mpass) = .true.
          call getmat(q(i10),q(i10),q(i10),q(i10),q(i10),q(i10),
     *    array,num,o1e,ionsec)
          if (mpass.eq.1) then
           call wrt3(q(i10),l2,ibl7s,num8)
          else if (mpass.eq.2) then
           call wrt3(q(i10),l2,ibl7t,num8)
          else
           call wrt3(q(i10),l2,ibl7f,num8)
          endif
         enddo
        else
         do loop =1,3
          o1e(loop) = .true.
          o1e(loop+3) = .false.
         enddo
         call getmat(q(i10),q(i20 ),q(i30),q(i10),q(i10),q(i10),
     *   array,num,o1e,ionsec)
        endif
_IF(parallel)
       endif
_ENDIF
      endif
c
_IF(parallel)
      if(omaster) then
_ENDIF
       if(otriang) then
c      restore s matrix
        call rdedx(q(i10),l2,ibl7s,num8)
       else
        call wrt3(q(i10),l2,ibl7s,num8)
        call wrt3(q(i20),l2,ibl7t,num8)
        call wrt3(q(i30),l2,ibl7f,num8)
       endif
c
c ----- transform the s-matrix
c
       call tranp(q(i10),q(i20))
       call wrt3(q(i20),l2,ibl7st,num8)
_IF(parallel)
      endif
      oswed3(4) = .true.
      oswed3(8) = .true.
_ENDIF
c
c ----- reset core
c
      call gmem_free(i10)
c
c     ----- establish buffers for the two electron integral file(s)
c
 130  continue

c     opank = .false.
c     if (nopk .eq. (-1)) opank = .true.
c     if ((zscftp .ne.zrhf .or. na .ne. nb) .and. nopk .ne. 1) opank =
c    +     .true.

      ogrhf = .false.

c
c     ----- execute scf procedure -----
c
      if (zscftp .eq. zgrhf) ogrhf = .true.
      sz = dzero
      s2 = dzero
c
      if(.not.odscf) then
       lfile=m2file
       do 141 i=1,lfile
       lotape(i)=m2tape(i)
       liblk(i)  =m2blk(i)
 141   llblk(i)  =liblk(i)-m2last(i)
      endif
c
c        atomic density startup for open and closed shells
c        modified for correct restart
c
c$$$      if (guess.eq.'atdens'.or.guess.eq.'atoms' .or.
c$$$     +    (zguess.eq.'mosaved'. and. (guess.eq.'atoms'.or.
c$$$     +      guess.eq.'atdens') ) ) then
c$$$c
c$$$c         do one cycle scf with the density matrix from denat
c$$$c
c$$$         call denscf(q,zscftp)
c$$$
c$$$         if (irest.ne.0) go to 160
c$$$c
c$$$c        set zguess to 'anything' so denscf will be called but once
c$$$c
c$$$         zguess = 'anything'
c$$$         guess  = 'anything'
c$$$      end if
c
c
c VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
c
      if (onews) then

         if (zscftp .eq. zrhf ) then
            ouhf = .false.
         else if (zscftp .eq. zuhf )then
            ouhf=.true.
         else
            call caserr(
     &           'newscf only available for (closed shell) rhf and uhf')
         endif

chvd     call newscf(q,ouhf,odscf,en,ehf,etot,sz,s2,ek,vir)

         call newscf_prepare(q,ouhf,odscf,en,ehf,etot,sz,s2,ek,vir)
         goto 160
      endif
c
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c

      if(zscftp . eq. zcas(1) ) go to 150
      if(zscftp . eq. zcas(2) ) go to 250
      if (zscftp .eq. zgvb ) go to 140
      if (zscftp .eq. zgrhf) go to 145

      if (zscftp .eq. zuhf )go to 120
      if (zscftp .ne. zuhf .and. na .ne. nb) go to 100
c
c     ----- closed shell restricted hartree fock -----
c

      if (zscftp.eq.zrhf) then
        if (odscf) then
              if (oscrf) then
                 call dscrf(q)
              else
                 if(odnew) then
                  call drcdrv(q,'scf')
                  goto 160
                 else
                  call drhfcl(q,c,czan)
                 endif
              endif
        else
c
c  first determine total amount of free memory available
         lwor = igmem_max_memory()
c  now determine core required given that lwor is available
c
cdft
           lscdf = 0
           if(osc)lscdf = l2
c  memory for grid - @ make this switchable
           if (osc.or.odenfu) then
            lscdf = lscdf + 50*194*nat
            write(6,*)'lscdf',lscdf
           endif

            nex = 0
            if(ocryst)nex = 1
            if(lwor.gt.((35+nex)*l2+num+lscdf)) then
               len=(31+nex)*l2+num+lscdf
               omem=.true.
            else
               len=len2+2*num+lscdf + nex*l2
               omem=.false.
            endif

_IF1()c         if(lwor.gt.(35*l2+num)) then
_IF1()c         len=31*l2+num
_IF1()c         omem=.true.
_IF1()c         else
_IF1()c         len=len2+2*num
_IF1()c         omem=.false.
_IF1()c         endif
_IF(unicos)
         nss=max(2,min(10,(lwor-len)/len2))
_ENDIF
_IF(convex)
         nss=max(2,min(18,(lwor-len)/len2))
_ENDIF
_IF(titan)
         nss=max(2,min(5,(lwor-len)/len2))
_ENDIF
_IFN(unicos,convex,titan)
      nss=2
_ENDIF
         lword=nss*len2+len
c 
         if (oso .and. ozora) call scf_so(q)
c spin orbit zora calculation sf
c
c     now allocate required core 
c 
c        ibase = igmem_alloc(lword)
c ------ force scrf calc through 'd' routine not 'm'
         odvsg=.true.
         if(oscrf) odvsg=.false.
c ------ carry on original code

              if(omem.and.odvsg.and..not.odisc) then
                 ibase = igmem_alloc(lword)
                 call rhfclm(q,q(ibase),c,czan,lword,nss)
              else
                 ibase = igmem_alloc(lword)
                 call rhfcld(q,q(ibase),c,czan,lword,nss)
              endif
         call gmem_free(ibase)
        endif
      endif
      go to 160
c
c     ----- open shell restricted hartree fock -----
c           e.r.davidson rohf not available in this
c           version.  ocbse type rohf in gvb.
c
  100 call caserr('davidson rohf not available')
c
c     ----- open shell unrestricted hartree fock -----
c
  120 if (odscf) then
         if(odnew) then
             call caserr('direct-UHF not available')
         else
             call duhfop(sz,s2,q)
         endif
      else
         call uhfop(sz,s2,q)
      endif
      go to 160
c
c     ----- generalized valence bond calculations -----
c
  140 if (odscf) then
         if(odnew) then
             call caserr('direct-GVB not available')
         else
             call drhfgvb(q,ogrhf)
         endif
      else
         call rhfgvb(q,ogrhf)
      endif
      go to 160
c
c     ----- generalized hartree fock
c
  145 if (odscf) then
         if(odnew) then
             call caserr('direct-GVB not available')
         else
             call drhfgvb(q,ogrhf)
         endif
      else
         call rhfgvb(q,ogrhf)
      endif
      go to 160
c
c     ----- complete active space scf calculation
c
  150 call casscf(q)
      call timana(9)
      go to 170
c
c     ----- 2nd order mcscf calculation
c
  250 call multi(q)
      call timana(10)
      go to 170
c
  160 call timana(5)
  170 enrgy = etot
c
c     ----- save data
c
      array(1) = en
      array(2) = ehf
      array(3) = etot
      array(4) = sz
      array(5) = s2
      array(6) = ek
      array(7) = vir
      do loop  = 8, 10
         array(loop) = 0.0d0
      enddo
c
      call secput(isect(494),m16,m1,iblk16)
      call wrt3(array,10,iblk16,idaf)
c
      call end_time_period(TP_SCF)
c
      return
      end

      logical function divide_and_conquer_requested()
INCLUDE(../m4/common/parcntl)
      divide_and_conquer_requested = (ipdiagmode.eq.IDIAG_PDSYEVD)
      return
      end

            
