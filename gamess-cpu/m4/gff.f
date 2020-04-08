c 
c  $Author: jmht $
c  $Date: 2010-08-10 16:49:47 +0200 (Tue, 10 Aug 2010) $
c  $Locker:  $
c  $Revision: 6176 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/gff.m,v $
c  $State: Exp $
c  
c ******************************************************
c ******************************************************
c             =   gff     =
c ******************************************************
c ******************************************************
      subroutine gf(core)
      implicit real*8  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c ...
c ... timing statistics
c ...
      real*8 begin, ebegin, timsec, walsec, tstart, estart
      common/statis/begin,ebegin,timsec(50),walsec(50),
     + tstart,estart
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      character*10 charwall
      dimension core(*)
      begin=cpulft(1)
      write(iwr,1100)
 1100 format(/1x,104('=')/
     *40x,28('*')/
     *40x,'ionization spectra gf module'/
     *40x,28('*')//)
c
c     allocate available memory
c
      i10 = igmem_alloc_all(lword)
c
      if(lword.lt.1)call caserr(
     *'insufficient memory for gf module')
      call initgf(core(i10),core(i10),lword)
      at=cpulft(1)
      write(iwr,1300) at ,charwall()
 1300 format(/1x,'end of gf module at',f8.2,' seconds',a10,' wall'/
     */1x,104('=')/)
c     free memory
      call gmem_free(i10)
c
      call clredx
      call timana(15)
      return
      end
      subroutine gfctl(v,wmeg,polst,wmeg_1,wmeg_2,wmeg_3,
     + pol_1,pol_2,pol_3, b_0,b_1,b_2,f_0,
     + polst_m,wmeg_m,sig_di,sig_dii,
     + indvec,jndvec,kndvec,lndvec)
      implicit real*8  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
      common/miscop/epsil(maxorb),wmega,conver,nocc,nunocc,nstart,ntot,
     *              nsym,imul(8,8),idorb(maxorb),
     *              nsecd,insecd(maxorb),nthd,inthd(maxorb),m,mm
      common /junk/ sigd,sigdi,sigdii,ders,dersi,dersii
      common /blkin/  a1,a2,a3,a5,c1,c2,c4,c6,d1,d3,d5,d6,dc1,dc2,dc4,
     *               dc6,dd1,dd3,dd5,dd6
      dimension v(*),wmeg(*),polst(*)
      dimension wmeg_1(*),wmeg_2(*),wmeg_3(*),pol_1(*),pol_2(*)
      dimension pol_3(*), b_0(*),b_1(*),b_2(*),f_0(*)
      dimension polst_m(*),wmeg_m(*),sig_di(*),sig_dii(*)
      dimension indvec(*),jndvec(*),kndvec(*),lndvec(*)
      data thr   /0.03d0/
      data z2,z3 /'second','third'/
      call dcopy(ntot,epsil,1,wmeg,1)
      thresh=thr/conver
      thrash=0.01d0/conver
      att=cpulft(1)
      if (nsecd.gt.0) then
      write(iwr,10010) z2,att
      do 100 kk=1,nsecd
      m=insecd(kk)
      mm=idorb(m)
      do 80 n=1,7
      wmega=wmeg(m)
      call sig2d(v,indvec,jndvec,kndvec,lndvec)
      pa=sigdi*conver
      pb=sigdii*conver
      pc=sigd*conver
      if(oprint(54))write (iwr,16) pa,pb,pc
      if(oprint(54))write (iwr,17) dersi,dersii,ders
      pols=1.0d0/(1.0d0-ders)
      wmegb=pols*(epsil(m)+sigd-wmega*ders)
      dif= dabs(wmegb-wmega)
      wmeg(m)=wmegb
      a=epsil(m)*conver
      b=wmeg(m)*conver
      if(oprint(54))write (iwr,4) m,a,b,pols
      if (dif.lt.thresh) go to 90
   80 continue
   90 continue
      wmega=wmeg(m)
      wmeg_m(m)=wmega
      call sig2d(v,indvec,jndvec,kndvec,lndvec)
      polst(m)=1.0d0/(1.0d0-ders)
  100 continue
      write (iwr,2)
      write (iwr,3)
      do 130 ii=1,nsecd
      i=insecd(ii)
      a=epsil(i)*conver
      b=wmeg_m(i)*conver
      write (iwr,4) i,a,b,polst(i)
  130 continue
c
      endif
c
      if(nthd.le.0) return
c
      write(iwr,1)
      att=cpulft(1)
      write(iwr,10010) z3,att
      do 700 kk=1,nthd
      m=inthd(kk)
      mm=idorb(m)
      call grapha(v,wmeg,indvec,jndvec,kndvec,lndvec)
      wmega=epsil(m)*0.92d0
      do 250 n=1,5
      if (n.gt.1) wmega=wmegm
      call sig2d(v,indvec,jndvec,kndvec,lndvec)
      call gc1c2(v,indvec,jndvec,kndvec,lndvec)
      call gc4c6(v,indvec,jndvec,kndvec,lndvec)
      call grphda(v,indvec,jndvec,kndvec,lndvec)
      call grphdb(v,indvec,jndvec,kndvec,lndvec)
      pa=sigdi*conver
      pb=sigdii*conver
      pc=sigd*conver
      if(oprint(54))write (iwr,16) pa,pb,pc
      if(oprint(54))write (iwr,17) dersi,dersii,ders
      pa=a1*conver
      pb=a2*conver
      pc=a3*conver
      pd=a5*conver
      if(oprint(54))write (iwr,11) pa,pb,pc,pd
      pa=c1*conver
      pb=c2*conver
      if(oprint(54))write (iwr,12) pa,pb,dc1,dc2
      pa=c4*conver
      pb=c6*conver
      if(oprint(54))write (iwr,13) pa,pb,dc4,dc6
      pa=d1*conver
      pb=d3*conver
      if(oprint(54))write (iwr,14) pa,pb,dd1,dd3
      pa=d5*conver
      pb=d6*conver
      if(oprint(54))write (iwr,15) pa,pb,dd5,dd6
      tder=ders+dc1+dc2+dc4+dc6+dd1+dd3+dd5+dd6
      tsig=sigd+a1+a2+a3+a5+c1+c2+c4+c6+d1+d3+d5+d6
      pols=1.0d0/(1.0d0-tder)
      wmegb=pols*(epsil(m)+tsig-wmega*tder)
      dif= dabs(wmegb-wmega)
      wmegm=wmegb
      a=epsil(m)*conver
      b=wmegm*conver
      if(oprint(54))write (iwr,4) m,a,b,pols
      if (dif.lt.thresh) go to 260
  250 continue
  260 continue
      wmega=wmegm
      if (dif.gt.thrash) then
      call sig2d(v,indvec,jndvec,kndvec,lndvec)
      call gc1c2(v,indvec,jndvec,kndvec,lndvec)
      call gc4c6(v,indvec,jndvec,kndvec,lndvec)
      call grphda(v,indvec,jndvec,kndvec,lndvec)
      call grphdb(v,indvec,jndvec,kndvec,lndvec)
      endif
      tder=ders+dc1+dc2+dc4+dc6+dd1+dd3+dd5+dd6
      polst_m(m)=1.0d0/(1.0d0-tder)
      wmeg_m(m)=wmegm
      pa=sigdi*conver
      pb=sigdii*conver
      pc=sigd*conver
      if(oprint(54))write (iwr,16) pa,pb,pc
      if(oprint(54))write (iwr,17) dersi,dersii,ders
      pa=a1*conver
      pb=a2*conver
      pc=a3*conver
      pd=a5*conver
      if(oprint(54))write (iwr,11) pa,pb,pc,pd
      pa=c1*conver
      pb=c2*conver
      if(oprint(54))write (iwr,12) pa,pb,dc1,dc2
      pa=c4*conver
      pb=c6*conver
      if(oprint(54))write (iwr,13) pa,pb,dc4,dc6
      pa=d1*conver
      pb=d3*conver
      if(oprint(54))write (iwr,14) pa,pb,dd1,dd3
      pa=d5*conver
      pb=d6*conver
      if(oprint(54))write (iwr,15) pa,pb,dd5,dd6
      b=-(c2+c4+d3+d5)/sigd
      b_0(m)=b
      s3=a1+a2+a3+a5+c1+c2+c4+c6+d1+d3+d5+d6
      if(oprint(54))write (iwr,24) b
      db=-(dc2+dc4+dd3+dd5+b*ders)/sigd
      tsig=sigd+s3/(1.0d0+b)
      tder=ders+(dc1+dc2+dc4+dc6+dd1+dd3+dd5+dd6)/(1.0d0+b)-
     *db*s3/(1.0d0+b)**2
      pols=1.0d0/(1.0d0-tder)
      wmegb=pols*(epsil(m)+tsig-wmega*tder)
      pol_1(m)=pols
      wmeg_1(m)=wmegb*conver
      wmega=wmegm
      b1=-(c2+d3)/sigdii
      sig_dii(m) = sigdii
      b_1(m)=b1
      sig_di(m) = sigdi
      b2=-(c4+d5)/sigdi
      b_2(m)=b2
      if(oprint(54))write (iwr,25) b1,b2
      pa=1.0d0/(1.0d0+b1)-1.0d0
      pb=1.0d0/(1.0d0+b2)-1.0d0
      db1=(dc2+dd3)/sigdii-(c2+d3)*dersii/sigdii**2
      db2=(dc4+dd5)/sigdi-(c4+d5)*dersi/sigdi**2
      tsig=sigd+s3+(c2+d3+c1+d1)*pa+(c4+d5+c6+d6)*pb
      tder=ders+dc1+dc2+dc4+dc6+dd1+dd3+dd5+dd6
      tder=tder+(dc2+dd3+dd1+dc1)*pa-(c2+d3+d1+c1+a3)*db1/(1.0d0+b1)**2
      tder=tder+(dc4+dd5+dc6+dd6)*pb-(c4+d5+d6+c6+a5)*db2/(1.0d0+b2)**2
      pols=1.0d0/(1.0d0-tder)
      wmegb=pols*(epsil(m)+tsig-wmega*tder)
      wmeg_2(m)=wmegb*conver
      pol2=pols
      pol_2(m)=pol2
      scd=c1+c2+c4+c6+d1+d3+d5+d6
      dscd=dc1+dc2+dc4+dc6+dd1+dd3+dd5+dd6
      f=(b2*(c4+c6+d5+d6)+b1*(c1+c2+d1+d3))/scd
      f_0(m)=f
      if(oprint(54))write (iwr,27) f
      df=(db2*(c4+c6+d5+d6)+b2*(dc4+dc6+dd5+dd6)+db1*(c1+c2+d1+d3)+
     *   b1*(dc1+dc2+dd1+dd3))/scd-f*dscd/scd
      tder=ders+dscd/(1.0d0+f)-s3*df/(1.0d0+f)**2
      tsig=sigd+s3/(1.0d0+f)
      pols=1.0d0/(1.0d0-tder)
      wmegc=pols*(epsil(m)+tsig-wmega*tder)
      wmeg_3(m) =wmegc*conver
      pol_3(m)=pols
c
 700  continue
c
c     now access results from other nodes for printing
c
c
c     final printing of results
c
      do 7000 kk=1,nthd
      m=inthd(kk)
      mm=idorb(m)
      write (iwr,1)
      write (iwr,9)
      write (iwr,3)
      b=wmeg_m(m)*conver
      a=epsil(m)*conver
      write (iwr,4) m,a,b,polst_m(m)
      write (iwr,23)
c
      a=-epsil(m)*conver
      wmeg1=-wmeg_1(m)
      wmeg2=-wmeg_2(m)
      wmeg3=-wmeg_3(m)
      pol1=pol_1(m)
      pol2=pol_2(m)
      pol3=pol_3(m)
      b =b_0(m)
      b1=b_1(m)
      b2=b_2(m)
      sigdi = sig_di(m)
      sigdii= sig_dii(m)
      f=f_0(m)
      write (iwr,26)
      write (iwr,3)
      write (iwr,4) m,a,wmeg1,pol1
      write (iwr,28)
      write (iwr,3)
      write (iwr,4) m,a,wmeg2,pol2
      write (iwr,29)
      write (iwr,3)
      write (iwr,4) m,a,wmeg3,pol3
      write (iwr,10)
      if(dabs(b).gt.0.85d0.and.
     *   dabs(b1).gt.0.85d0.and.
     *   dabs(b2).gt.0.85d0.and. 
     *   dabs(f).gt.0.85d0) go to 500
      go to 505
  500 write (iwr,6) m,a,b,wmeg1,pol1,b1,b2,wmeg2,pol2,f,wmeg3,pol3
      write (iwr,7)
      return
  505 if (wmeg1.gt.15.0d0.and. dabs(b).le.0.85d0) go to 510
      go to 515
  510 if ((sigdi+sigdii)*conver.lt.0.6d0.and. dabs(f).le.0.85d0)
     * go to 560
      write (iwr,3)
      write (iwr,4) m,a,wmeg1,pol1
      go to 7000
  515 if(wmeg1.gt.15.0d0.and.dabs(b).gt.0.85d0.and.
     *   wmeg2.gt.15.0d0.and.dabs(b1).le.0.85d0.and.
     *   dabs(b2).le.0.85d0) go to 520
      go to 530
  520 if ((sigdi+sigdii)*conver.lt.0.6d0.and. 
     *     dabs(b1).le.0.85d0.and.dabs(b2).le.0.85d0) go to 560
      write (iwr,3)
      write (iwr,4) m,a,wmeg2,pol2
      go to 7000
  530 if(wmeg1.gt.15.0d0.and.dabs(b).gt.0.85d0.and.
     *   wmeg2.le.15.0d0.and.dabs(b1).le.0.85d0.and.
     *   dabs(b2).le.0.85d0) go to 550
      if(wmeg1.gt.15.0d0.and.dabs(b).gt.0.85d0.and.
     *   wmeg3.gt.15.0d0.and.dabs(f).le.0.85d0) go to 540
      go to 550
  540 write (iwr,3)
      write (iwr,4) m,a,wmeg3,pol3
      go to 7000
  550 if ((sigdi+sigdii)*conver.lt.0.6d0) go to 560
      go to 570
  560 write (iwr,3)
      write (iwr,4) m,a,wmeg3,pol3
      go to 7000
  570 write (iwr,3)
      write (iwr,4) m,a,wmeg2,pol2
 7000 continue
c
      return
 1    format(/1x,104('*')/)
10010 format(/1x,104('*')
     *//5x,'commence  ',a8,'order perturbation calculation ',
     *            'at',f8.2,' secs'/5x,64('=')/)
    2 format ('   second order perturbation theory, diagonal form,',
     * 'final results'/3x,62('='))
    3 format ('   orbital   hf-eigenvalue     corrected value    pole',
     *'strength'/3x,60('='))
    4 format (4x,i3,7x,f12.5,8x,f12.5,6x,f8.5/)
c   5 format ('   second order perturbation theory, complete form,final
c    *results'/)
    6 format (//2x,i3,11(1x,f10.5))
    7 format (//'   contact author with these results')
c   8 format (/10(2x,f10.6))
    9 format(/'   third order perturbation theory, diagonal form,',
     *'final  results'/1x,63('-')/)
   10 format (///'   renormalization of interaction, final results'///)
   11 format ('   a1  =  ',f12.6,'   a2  =  ',f12.6,'   a3   =  ',f12.6,
     *'   a5   =  ',f12.6)
   12 format ('   c1  =  ',f12.6,'   c2  =  ',f12.6,'   dc1  =  ',f12.6,
     *'   dc2  =  ',f12.6)
   13 format ('   c4  =  ',f12.6,'   c6  =  ',f12.6,'   dc4  =  ',f12.6,
     *'   dc6  =  ',f12.6)
   14 format ('   d1  =  ',f12.6,'   d3  =  ',f12.6,'   dd1  =  ',f12.6,
     *'   dd3  =  ',f12.6)
   15 format ('   d5  =  ',f12.6,'   d6  =  ',f12.6,'   dd5  =  ',f12.6,
     *'   dd6  =  ',f12.6)
   16 format ('   sigdi  =  ',f12.6,'   sigdii  =  ',f12.6,'   sigd  ='
     *,f12.6)
   17 format ('   dersi  =  ',f12.6,'   dersii  =  ',f12.6,'   ders  ='
     *,f12.6)
c  22 format ('   matrix of collected eigenvectors'///)
   23 format ('   renormalization of interaction, methods involving ',
     *'2nd and 3rd order diagrams'//)
   24 format ('   method a : b =  ',f10.6//)
   25 format ('   method b : b1 =  ',f10.6,'   b2 =  ',f10.6//)
   26 format(/'   final results for method a'/)
   27 format ('   method c : f =  ',f10.6//)
   28 format(/'   final results for method b'/)
   29 format(/'   final results for method c'/)
      end
      subroutine initgf(q,iq,lword)
      implicit real*8  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      parameter (mxorb1=maxorb+1)
      common/miscop/epsil(maxorb),wmega,conver,nocc,nunocc,nstart,ntot,
     +       nsym,imul(8,8),idorb(maxorb),
     +       nsecd,insecd(maxorb),nthd,inthd(maxorb),m,mm,
     +       mapie(maxorb),nmv
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
      dimension iworka(maxorb),iworkb(maxorb)
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
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common/junk/pop(maxorb),potn,core,ncolo(3),ncore,
     *mapcie(maxorb),map2(maxorb),nact,mapaie(maxorb),mapaei(maxorb)
     *,iqsec,nacta(maxorb),nactb(maxorb),nactc(5),isecor,
     *norb2,nseco(maxorb),norb3,nthird(maxorb),ntoto,iwork(maxorb),
     *nscrap(maxorb),nscrp,
     *evalue(maxorb),eocc(mxorb1),nbas,newb,ncol,ieval,ipop,ispp
     *,nirr,mult(8,8),isymao(maxorb),isymmo(maxorb),nsp
c
      integer lfile, lotape, liblk, llblk
      common /filel/ lfile,lotape(20),liblk(20),llblk(20)
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
      common/junkc/zjob,zdate,ztime,zprog,ztype,zspace(14),ztext(10)
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
      dimension q(*),iq(*)
      data yblank,xstar,thresh,two,ydomo,yvmo
     *     /' ','*',1.0d-5,2.0d0,'domo','vmo'/
      data z2,z3 /'second','third'/
      data xdash /'-'/
      data m0,m8/0,8/
      data m51,m29/51,29/
c
      nav = lenwrd()
      lwor4=lword*nav
      conver=27.2116d0
      call secget(isect(470),1005,iblka)
      call readi(nacta,mach(14)*nav,iblka,idaf)
      call readis(norb2,mach(16)*nav,idaf)
      call secget(isecor,1004,iblka)
      call rdedx(pop,mach(15),iblka,idaf)
      write(iwr,5004) yed(idaf),ibl3d,isecor
 5004 format(/5x,'dumpfile resides on ',a4,' at block ',i6//
     *          5x,'core hamiltonian to be retrieved from section ',i4/)
      write(iwr,2)nact,ncore,potn,core
    2 format(5x,'header block information :'/
     *       5x,'number of active orbitals ',i4/
     *       5x,'number of core   orbitals ',i4/
     *       5x,'nuclear repulsion energy ',e21.12/
     *       5x,'core energy              ',e21.12/)
      if(nact.gt.1) go to 10005
80000 call caserr('invalid number of active orbitals')
80010 call caserr('parameter error in gf preprocessor')
80020 call caserr(
     *'2nd and 3rd order calculations supressed')
80030 call caserr('insufficient memory for gf module')
10005 do 82003 i=1,nact
      j=mapaie(i)
82003 mapaei(j)=i
      lfile=m6file
      do i=1,lfile
       lotape(i)=m6tape(i)
       liblk  (i)=m6blk(i)
       llblk  (i)=liblk(i)-m6last(i)
      enddo
      ntot=ntoto
      if(ntot.ne.0)go to 5005
      ntot=nact
      do i=1,ntot
       iworka(i)=i
       mapie(i)=i
      enddo
      go to 8460
5005  continue
      call igthr(ntot,mapaie,iworka,iwork)
      call rddir(ntot,iworka)
      write(iwr,15002)
      write(iwr,9) (iworka(i),i=1,ntot)
      if(ntot.eq.1) go to 80000
      if(ntot.gt.nact.or.iworka(ntot).gt.nact) go to 80000
      call setsto(nact,0,mapie)
      do i=1,ntot
       mapie(iworka(i))=i
      enddo
 8460 nsecd=norb2
      if(nsecd)8462,8461,8463
8463  call rdba(nact,nsecd,insecd,nseco,mapaei,iworka,mapie)
      write(iwr,7) z2,nsecd
      write(iwr,9) (nseco(i),i=1,nsecd)
      go to 8462
 8461 write(iwr,8) z2
    7 format(//5x,a8,'order perturbation calculation requested ',
     *                'for the following',i4,'  orbitals :')
    8 format(/5x,a8,'order perturbation calculation suppressed'/)
 8462 nthd=norb3
      if(nthd)60006,8464,8465
 8465 call rdba(nact,nthd,inthd,nthird,mapaei,iworka,mapie)
      write(iwr,7) z3,nthd
      write(iwr,9) (nthird(i),i=1,nthd)
      go to 60006
 8464 write(iwr,8) z3
      if(nsecd.eq.0) go to 80020
60006 call secget(isect(490),m51,iblka)
      osym=.true.
      call readi(nirr,mach(13)*nav,iblka,idaf)
      nsym=nirr
      call icopy(64,mult,1,imul,1)
      call setsto(ntot,m0,idorb)
      call secget(iqsec,3,iblka)
      call rdchr(zjob,m29,iblka,idaf)
      call reads(evalue,mach(8),idaf)
      write(iwr,89005) iqsec,ztype,ztime,zdate,zjob,ztext,nbas,ncol
89005 format(/5x,'scf mo specifications restored from section ',i4,
     *          ' of dumpfile'//
     *       5x,'header block information :'/
     *       5x,a7,' vectors created at ',
     *          a8,' on ',a8,' in the job ',a8/
     *       5x,'with the title :',5x,10a8/
     *       5x,'no. of gtos        ',i4/
     *       5x,'no. of scf vectors ',i4//)
      call symvec(q,isymao,isymmo,nbas,ncol,iblka)
      ic=0
      do 44 i=1,nact
      ida=mapie(i)
      if(ida.eq.0)go to 44
      ic=ic+1
      if(idorb(ida).ne.0)go to 80010
      idorb(ida)=isymmo(mapaie(ida))
 44   continue
      if(ic.ne.ntot)go to 80010
      nocc=0
      nunocc=0
      do 89006 i=1,ntot
      jscf=mapaie(iworka(i))
      epsil(i)=evalue(jscf)
      qio=eocc(jscf)
      if( dabs(qio-two).le.thresh) nocc=nocc+1
      if(qio.le.thresh) nunocc=nunocc+1
89006 continue
      if(nocc+nunocc.ne.ntot) go to 80010
      if(nunocc.le.0) go to 80010
      call vsmul(epsil,1,conver,evalue,1,ntot)
      nstart=nocc+1
      if(nsecd.gt.0.or.nthd.gt.0) go to 89100
      nsecd=nocc
      nthd=nocc
      do i=1,nocc
       insecd(i)=i
       inthd(i)=i
      enddo
89100 nmv=0
      do 90 i=1,ntot
      ii=idorb(i)
      do 90 j=1,i
      jj=idorb(j)
      ijm=imul(ii,jj)
      do 90 k=1,i
      kk=idorb(k)
      le=k
      if (i.eq.k) le=j
      do 85 l=1,le
      ll=idorb(l)
      if (ijm.eq.imul(kk,ll)) nmv=nmv+1
   85 continue
   90 continue
    3 format(  5x,'case :',5x,10a8
     *       //5x,'number of active orbitals         ',i4/
     *        /5x,'number of doubly occupied orbitals',i4/
     *        /5x,'number of virtual orbitals        ',i4)
    5 format (/50x,'scf orbital energies (au)'/50x,25('*'))
    6 format(/10x,8f14.7)
   11 format (/50x,'scf orbital energies (ev)'/50x,25('*'))
    9 format(/20(3x,i3))
   13 format(5x,'no symmetry is assumed for the present case'//)
   14 format(5x,'input orbital specifications :'//
     *  5x,'orbital   tran4    scf    symmetry   occup   2-nd order',
     *     '    3-rd order'/
     *  5x,'  no.     label   label                      calculation',
     *     '   calculation')
  144 format(4x,120a1)
   16 format(5x,i4,6x,i4,4x,i4,7x,i1,8x,a4,8x,a8,6x,a8)
   17 format (//5x,i7,'   integrals are to be read in and processed')
   24 format(//5x,'insufficient storage available, increase by',
     *      i10,' words'//)
   25 format(/5x,'main core not used:',2x,i10,' words'//)
      write(iwr,144) (xstar,i=1,93)
      write(iwr,3)ztitle,ntot,nocc,nunocc
      write(iwr,144) (xstar,i=1,93)
15002 format(/5x,'the following mo''s included in active list :'/)
      write(iwr,5)
      write(iwr,6) (epsil(i),i=1,ntot)
      write(iwr,11)
      write(iwr,6) (evalue(i),i=1,ntot)
      if(.not.osym) write(iwr,13)
      write(iwr,144) (xdash,i=1,72)
      write(iwr,14)
      write(iwr,144) (xdash,i=1,72)
      isat2=1
      isat3=1
      do 20002 i=1,ntot
      jtr=iworka(i)
      jscf=mapaie(jtr)
      ytype=yvmo
      ycalc2=yblank
      ycalc3=yblank
      if(i.le.nocc) ytype=ydomo
      nsy=idorb(i)
      if(i.ne.insecd(isat2)) go to 20011
      isat2=isat2+1
      ycalc2=xstar
20011 if(i.ne.inthd(isat3)) go to 20001
      isat3=isat3+1
      ycalc3=xstar
20001 write(iwr,16) i,jtr,jscf,nsy,ytype,ycalc2,ycalc3
20002 continue
      write(iwr,144) (xdash,i=1,72)
      ixx=ntot*(ntot+1)/2
      isigm=nmv+1
      ipol=isigm+ ixx
      i1 = ipol + ntot
      i2 = i1   + ntot
      i3 = i2   + ntot
      i4 = i3   + ntot
      i5 = i4   + ntot
      i6 = i5   + ntot
      i7 = i6   + ntot
      i8 = i7   + ntot
      i9 = i8   + ntot
      i10= i9   + ntot
      i11= i10  + ntot
      i12= i11  + ntot
      i13= i12  + ntot
      i14= i13  + ntot
c
c     iyy=ixx*(ntot+ntot+1)/3
      iyy=0
      do 60009 i=1,ntot
      do 60009 j=1,i
      do 60009 k=1,i
60009 iyy=iyy+1
      iln=ixx
      if(iln.le.3) iln=iln+1
      lenr=(isigm+15*ntot+iln-1)*nav
      ju1=lenr+1
      ju2=ju1+ntot
      ju3=ju2+ntot
      ju4=ju3+ixx
      iz1=ju4+iyy -1
      if(lwor4.ge.iz1) go to 96
      need=(iz1-lwor4)/nav+1
      write(iwr,24) need
      go to 80030
   96 nfree=(lwor4-iz1)/nav
      write(iwr,25) nfree
      write(iwr,1)
 1    format(/1x,'transformed integral files'/1x,26('*'))
      call filprn(m6file,m6blk,m6last,m6tape)
      write(iwr,17) nmv
      call inhelp(iq(ju1),iq(ju2),iq(ju3),iq(ju4),ntot,idorb,imul,m8)
c     call linkin(iq(ju4),iq(ju3),iq(ju2),iq(ju1))
      if(ntot.ne.nact) go to 23000
      call rdv01(q,iq(ju4),iq(ju3),iq(ju2),iq(ju1) )
      go to 25000
23000 do 23001 i=1,ntot
      if(iworka(i).ne.i) go to 23002
23001 continue
      call rdv02(q,iq(ju4),iq(ju3),iq(ju2),iq(ju1) )
      go to 25000
23002 if(iworka(ntot).ne.nact) go to 23003
      call rdv03(q,iq(ju4),iq(ju3),iq(ju2),iq(ju1) )
      go to 25000
23003 call rdv04(q,iworka(ntot),iq(ju4),iq(ju3),iq(ju2),iq(ju1) )
25000 call gfctl(q(1),q(isigm),q(ipol),q(i1),q(i2),q(i3),
     +           q(i4),q(i5),q(i6),q(i7),q(i8),q(i9),
     +           q(i10),q(i11),q(i12),q(i13),q(i14),
     +           iq(ju4),iq(ju3),iq(ju2),iq(ju1))
      return
      end
      subroutine rddir(n,jvec)
      implicit real*8  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension jvec(*)
      nm1=n-1
      do 23 i=1,nm1
      ip1=i+1
      j=jvec(i)
      do 22 k=ip1,n
      if(jvec(k).gt.j) go to 22
      if(jvec(k).eq.j) call caserr(
     *'duplicate orbital index in active set')
      it=jvec(k)
      jvec(k)=j
      jvec(i)=it
      j=it
   22 continue
   23 continue
      return
      end
      subroutine rdba(nact,nr,istor,iw1,mapaei,iwa,mapie)
      implicit real*8  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension istor(*),iw1(*),mapaei(*),iwa(*),mapie(*)
      do 30 i=1,nr
 30   iw1(i)=mapaei(iw1(i))
      nm1=nr-1
      do 23 i=1,nm1
      ip1=i+1
      j=iw1(i)
      do 22 k=ip1,nr
      if(j-iw1(k)) 22,19,21
   19 call caserr('orbital doubly defined')
   21 it=iw1(k)
      iw1(k)=j
      iw1(i)=it
      j=it
   22 continue
   23 continue
      if(nr.gt.nact.or.iw1(nr).gt.nact) call caserr(
     *'invalid orbital parameters specified')
      ii=0
      do 7005 i=1,nr
      j=mapie(iw1(i))
      if(j.le.0) go to 7005
      ii=ii+1
      istor(ii)=j
 7005 continue
      if(nr.eq.ii) return
      nr=ii
      do 7006 i=1,nr
 7006 iw1(i)=iwa(istor(i))
      return
      end
      subroutine rdv01(v,indvec,jndvec,kndvec,lndvec)
      implicit real*8  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension v(*)
      dimension indvec(*),jndvec(*),kndvec(*),lndvec(*)
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
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
c
      integer lfile, lotape, liblk, llblk
      common /filel/ lfile,lotape(20),liblk(20),llblk(20)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common/craypk/i205(1360)
      common/miscop/epsil(maxorb),wmega,conver,nocc,nunocc,nstart,ntot,
     *              nsym,imul(8,8),idorb(maxorb),
     *              nsecd,insecd(maxorb),nthd,inthd(maxorb),m,mm
     *             ,mapie(maxorb),nmv
      common/blkin/gin(510),numi
      indv (i,j,k,l) =  lndvec(l) +
     *                  indvec(k+jndvec(j+kndvec(i)))
c
      ic=0
      call vclr(v,1,nmv)
      itot=0
      do 8 ifile=1,lfile
      iunit=lotape(ifile)
      call search(liblk(ifile),iunit)
      call find(iunit)
      lbl=llblk(ifile)
    6 lbl=lbl+1
      call get(gin(1),mw)
      if(mw.eq.0) go to 8
      itot=itot+numi
      if(lbl.ne.0) call find(iunit)
      int4=1
      call unpack(gin(num2e+1),lab816,i205,numlab)
      do 5 num=1,numi
      j=i205(int4  )
      i=i205(int4+1)
      l=i205(int4+2)
      k=i205(int4+3)
      ijm=imul(idorb(i),idorb(j))
      klm=imul(idorb(k),idorb(l))
      if(ijm.ne.klm) go to 5
      ic=ic+1
      v(indv(i,j,k,l))=gin(num)
    5 int4=int4+4
      if(lbl.ne.0) go to 6
    8 continue
      write(iwr,999) itot
  999 format(/5x,i7,'  2-electron integrals scanned'/)
      if(ic.ne.nmv) write(iwr,1000) ic
 1000 format(/5x,'*** warning ***'/
     *         5x,'only ',i7,'  symmetry allowed integrals found on',
     *            ' final mainfile'/
     *         5x,'remaining integrals assumed zero'//)
      return
      end
      subroutine rdv02(v,indvec,jndvec,kndvec,lndvec)
      implicit real*8  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension v(*)
      dimension indvec(*),jndvec(*),kndvec(*),lndvec(*)
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
      integer lfile, lotape, liblk, llblk
      common /filel/ lfile,lotape(20),liblk(20),llblk(20)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      common/craypk/i205(1360)
      common/blkin/gin(510),numi
      common/miscop/epsil(maxorb),wmega,conver,nocc,nunocc,nstart,ntot,
     *              nsym,imul(8,8),idorb(maxorb),
     *              nsecd,insecd(maxorb),nthd,inthd(maxorb),m,mm
     *             ,mapie(maxorb),nmv
      indv (i,j,k,l) =  lndvec(l) +
     *                  indvec(k+jndvec(j+kndvec(i)))
c
      ic=0
      call vclr(v,1,nmv)
      itot=0
      do 8 ifile=1,lfile
      iunit=lotape(ifile)
      call search(liblk(ifile),iunit)
      call find(iunit)
      lbl=llblk(ifile)
    6 lbl=lbl+1
      call get(gin(1),mw)
      if(mw.eq.0) go to 8
      itot=itot+numi
      if(lbl.ne.0) call find(iunit)
      int4=1
      call unpack(gin(num2e+1),lab816,i205,numlab)
      do 5 num=1,numi
      i=i205(int4+1)
      if(i.le.ntot) go to 10
      if(ic.lt.nmv) go to 7
      go to 12
 10   continue
      j=i205(int4  )
      l=i205(int4+2)
      k=i205(int4+3)
      ijm=imul(idorb(i),idorb(j))
      klm=imul(idorb(k),idorb(l))
      if(ijm.ne.klm) go to 5
      ic=ic+1
      v(indv(i,j,k,l))=gin(num)
    5 int4=int4+4
    7 if(lbl.ne.0) go to 6
    8 continue
   12 write(iwr,999) itot
  999 format(/5x,i7,'  2-electron integrals scanned'/)
      if(ic.ne.nmv) write(iwr,1000) ic
 1000 format(/5x,'*** warning ***'/
     *         5x,'only ',i7,'  symmetry allowed integrals found on',
     *            ' tran4 final mainfile'/
     *         5x,'remaining integrals assumed zero'/)
      return
      end
      subroutine rdv03(v,indvec,jndvec,kndvec,lndvec)
      implicit real*8  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension v(*)
      dimension indvec(*),jndvec(*),kndvec(*),lndvec(*)
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
      integer lfile, lotape, liblk, llblk
      common /filel/ lfile,lotape(20),liblk(20),llblk(20)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      common/craypk/i205(1360)
      common/blkin/gin(510),numi
      common/miscop/epsil(maxorb),wmega,conver,nocc,nunocc,nstart,ntot,
     *              nsym,imul(8,8),idorb(maxorb),
     *              nsecd,insecd(maxorb),nthd,inthd(maxorb),m,mm
     *             ,mapie(maxorb),nmv
      indv (i,j,k,l) =  lndvec(l) +
     *                  indvec(k+jndvec(j+kndvec(i)))
c
      ic=0
      call vclr(v,1,nmv)
      itot=0
      do 8 ifile=1,lfile
      iunit=lotape(ifile)
      call search(liblk(ifile),iunit)
      call find(iunit)
      lbl=llblk(ifile)
    6 lbl=lbl+1
      call get(gin(1),mw)
      if(mw.eq.0) go to 8
      itot=itot+numi
      if(lbl.ne.0) call find(iunit)
      int4=1
      call unpack(gin(num2e+1),lab816,i205,numlab)
      do 5 num=1,numi
      i=mapie(i205(int4+1))
      if(i.eq.0)go to 5
      j=mapie(i205(int4  ))
      if(j.eq.0) go to 5
      k=mapie(i205(int4+3))
      if(k.eq.0) go to 5
      l=mapie(i205(int4+2))
      if(l.eq.0) go to 5
      ijm=imul(idorb(i),idorb(j))
      klm=imul(idorb(k),idorb(l))
      if(ijm.ne.klm) go to 5
      ic=ic+1
      v(indv(i,j,k,l))=gin(num)
    5 int4=int4+4
      if(lbl.ne.0) go to 6
    8 continue
      write(iwr,999) itot
  999 format(//5x,i7,'  2-electron integrals scanned'//)
      if(ic.ne.nmv) write(iwr,1000) ic
 1000 format(//5x,'*** warning ***'/
     *         5x,'only ',i7,'  symmetry allowed integrals found on',
     *            ' tran4 final mainfile'/
     *         5x,'remaining integrals assumed zero'//)
      return
      end
      subroutine rdv04(v,imax,indvec,jndvec,kndvec,lndvec)
      implicit real*8  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension v(*)
      dimension indvec(*),jndvec(*),kndvec(*),lndvec(*)
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
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
c
      integer lfile, lotape, liblk, llblk
      common /filel/ lfile,lotape(20),liblk(20),llblk(20)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common/craypk/i205(1360)
      common/blkin/gin(510),numi
      common/miscop/epsil(maxorb),wmega,conver,nocc,nunocc,nstart,ntot,
     *              nsym,imul(8,8),idorb(maxorb),
     *              nsecd,insecd(maxorb),nthd,inthd(maxorb),m,mm
     *             ,mapie(maxorb),nmv
      indv (i,j,k,l) =  lndvec(l) +
     *                  indvec(k+jndvec(j+kndvec(i)))
c
      ic=0
      call setsto(1360,0,i205)
      call vclr(v,1,nmv)
      itot=0
      do 8 ifile=1,lfile
      iunit=lotape(ifile)
      call search(liblk(ifile),iunit)
      call find(iunit)
      lbl=llblk(ifile)
    6 lbl=lbl+1
      call get(gin(1),mw)
      if(mw.eq.0) go to 8
      itot=itot+numi
      if(lbl.ne.0) call find(iunit)
      int4=1
      call unpack(gin(num2e+1),lab816,i205,numlab)
      do 5 num=1,numi
      i=i205(int4+1)
      if(i.le.imax) go to 10
      if(ic.ge.nmv) go to 12
      go to 7
   10 i=mapie(i)
      if(i.eq.0) go to 5
      j=mapie(i205(int4  ))
      if(j.eq.0) go to 5
      k=mapie(i205(int4+3))
      if(k.eq.0) go to 5
      l=mapie(i205(int4+2))
      if(l.eq.0) go to 5
      ijm=imul(idorb(i),idorb(j))
      klm=imul(idorb(k),idorb(l))
      if(ijm.ne.klm) go to 5
      ic=ic+1
      v(indv(i,j,k,l))=gin(num)
    5 int4=int4+4
    7 if(lbl.ne.0) go to 6
    8 continue
   12 write(iwr,999) itot
  999 format(/5x,i7,'  2-electron integrals scanned'/)
      if(ic.ne.nmv) write(iwr,1000) ic
 1000 format(/5x,'*** warning ***'/
     *         5x,'only ',i7,'  symmetry allowed integrals found on',
     *            ' tran4 final mainfile'/
     *         5x,'remaining integrals assumed zero'/)
      return
      end
      subroutine inhelp (lndvec,kndvec,jndvec,indvec,ntot,idorb,imul
     * ,ndim)
      implicit real*8  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension indvec(*),jndvec(*),kndvec(*),lndvec(*)
      dimension idorb(*),imul(ndim,*)
      isum=0
      ind=0
      do 20 i=1,ntot
      ii=idorb(i)
      do 20 j=1,i
      jj=idorb(j)
      ijm=imul(ii,jj)
      do 20 k=1,i
      kk=idorb(k)
      ind=ind+1
      indvec(ind)=isum
      le=k
      if(i.eq.k) le=j
      do 15 l=1,le
      ll=idorb(l)
      if(ijm.eq.imul(kk,ll)) isum=isum+1
   15 continue
   20 continue
      isum=0
      isum2=0
      ind=0
      do 60 i=1,ntot
      kndvec(i)=isum
      isum=isum+i
      lsym=idorb(i)
      isum3=0
      do 50 j=1,i
      ind=ind+1
      jndvec(ind)=isum2
      isum2=isum2+i
      if(idorb(j).eq.lsym) isum3=isum3+1
   50 continue
      lndvec(i)=isum3
   60 continue
      return
      end
      subroutine sig2d(v,indvec,jndvec,kndvec,lndvec)
      implicit real*8  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      common/miscop/epsil(maxorb),wmega,conver,nocc,nunocc,nstart,ntot,
     *              nsym,imul(8,8),idorb(maxorb),
     *              nsecd,insecd(maxorb),nthd,inthd(maxorb),m,mm
      common /junk/ sigd,sigdi,sigdii,ders,dersi,dersii
      dimension v(*)
      dimension indvec(*),jndvec(*),kndvec(*),lndvec(*)
      indv (i,j,k,l) =  lndvec(l) +
     *                  indvec(k+jndvec(j+kndvec(i)))
      sigd=0.0d0
      sigdi=0.0d0
      sigdii=0.0d0
      ders=0.0d0
      dersi=0.0d0
      dersii=0.0d0
      do 200 i=1,nocc
      ii=idorb(i)
      epsa=wmega-epsil(i)
      ijm=imul(ii,mm)
      ii3=m
      ii4=i
      if (m.ge.i) go to 10
      ii3=i
      ii4=m
   10 continue
      do 190 j=1,nocc
      jj=idorb(j)
      epsb=epsa-epsil(j)
      jj3=m
      jj4=j
      if (m.ge.j) go to 20
      jj3=j
      jj4=m
   20 continue
      do 180 ia=nstart,ntot
      iaa=idorb(ia)
      if (ijm.ne.imul(jj,iaa)) go to 180
      epsc=epsb+epsil(ia)
      rz=1.0d0/epsc
      i1=ia
      i2=j
      i3=ii3
      i4=ii4
      if (i1.gt.i3) go to 40
      if (i1.eq.i3) go to 30
      it1=i1
      it2=i2
      i1=i3
      i2=i4
      i3=it1
      i4=it2
      go to 40
   30 if (i2.ge.i4) go to 40
      it1=i2
      i2=i4
      i4=it1
   40 continue
      ind1=indv(i1,i2,i3,i4)
      j3=jj3
      j4=jj4
      j1=ia
      j2=i
      if (j1.gt.j3) go to 70
      if (j1.eq.j3) go to 60
      it1=j1
      it2=j2
      j1=j3
      j2=j4
      j3=it1
      j4=it2
      go to 70
   60 if (j2.ge.j4) go to 70
      it1=j2
      j2=j4
      j4=it1
   70 continue
      ind2=indv(j1,j2,j3,j4)
      pa=v(ind1)
      pb=v(ind2)
      pi=pa*(pa+pa-pb)*rz
      sigdi=sigdi+pi
      dersi=dersi-pi*rz
  180 continue
  190 continue
  200 continue
      do 400 ia=nstart,ntot
      iaa=idorb(ia)
      ijm=imul(mm,iaa)
      epsa=wmega-epsil(ia)
      ii1=ia
      ii2=m
      if (ia.ge.m) go to 210
      ii1=m
      ii2=ia
  210 continue
      do 390 ib=nstart,ntot
      ibb=idorb(ib)
      epsb=epsa-epsil(ib)
      jj3=ib
      jj4=m
      if (ib.ge.m) go to 220
      jj3=m
      jj4=ib
  220 continue
      do 380 i=1,nocc
      ii=idorb(i)
      if (imul(ii,ibb).ne.ijm) go to 380
      epsc=epsb+epsil(i)
      rz=1.0d0/epsc
      i1=ii1
      i2=ii2
      i3=ib
      i4=i
      if (i1.gt.i3) go to 250
      if (i1.eq.i3) go to 240
      it1=i1
      it2=i2
      i1=i3
      i2=i4
      i3=it1
      i4=it2
      go to 250
  240 if (i2.ge.i4) go to 250
      it1=i2
      i2=i4
      i4=it1
  250 continue
      ind1=indv(i1,i2,i3,i4)
      j3=jj3
      j4=jj4
      j1=ia
      j2=i
      if (j1.gt.j3) go to 280
      if (j1.eq.j3) go to 270
      it1=j1
      it2=j2
      j1=j3
      j2=j4
      j3=it1
      j4=it2
      go to 280
  270 if (j2.ge.j4) go to 280
      it1=j2
      j2=j4
      j4=it1
  280 continue
      ind2=indv(j1,j2,j3,j4)
      pa=v(ind1)
      pb=v(ind2)
      pi=pa*(pa+pa-pb)*rz
      sigdii=sigdii+pi
      dersii=dersii-pi*rz
  380 continue
  390 continue
  400 continue
      sigd=sigdi+sigdii
      ders=dersi+dersii
      return
      end
      subroutine grapha(v,a,indvec,jndvec,kndvec,lndvec)
      implicit real*8  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      common/miscop/epsil(maxorb),wmega,conver,nocc,nunocc,nstart,ntot,
     *              nsym,imul(8,8),idorb(maxorb),
     *              nsecd,insecd(maxorb),nthd,inthd(maxorb),m,mm
      common /blkin/  a1,a2,a3,a5,c1,c2,c4,c6,d1,d3,d5,d6,dc1,dc2,dc4,
     *               dc6,dd1,dd3,dd5,dd6
      dimension a(*),v(*)
      dimension indvec(*),jndvec(*),kndvec(*),lndvec(*)
      ind(m,n)= max(m,n)*(max(m,n)-1)/2 + min(m,n)
      indv (i,j,k,l) =  lndvec(l) +
     *                  indvec(k+jndvec(j+kndvec(i)))
      a1=0.0d0
      a2=0.0d0
      a3=0.0d0
      a5=0.0d0
      ix=0
      do 100 i=1,ntot
      ii=idorb(i)
      do 90 j=1,i
      ix=ix+1
      a(ix)=0.0d0
      if (imul(ii,idorb(j)).ne.1) go to 90
      if(m.gt.j) go to 50
      pa=v(indv(i,j,m,m))
      pb=v(indv(i,m,j,m))
      go to 80
   50 if(i-m) 55,56,57
   55 pa=v(indv(m,m,i,j))
      pb=v(indv(m,i,m,j))
      go to 80
   56 pa=v(indv(m,m,i,j))
      pb=pa
      go to 80
   57 pa=v(indv(i,j,m,m))
      pb=v(indv(i,m,m,j))
   80 a(ix)=pa+pa-pb
   90 continue
  100 continue
      do 600 i=1,nocc
      ii=idorb(i)
      epsa=epsil(i)
      do 590 ia=nstart,ntot
      iaa=idorb(ia)
      ijm=imul(ii,iaa)
      epsb=epsa-epsil(ia)
      jj1=ia
      jj2=i
      do 580 ib=nstart,ntot
      ibb=idorb(ib)
      epsib=epsil(ib)
      epsc=epsb-epsib
      ii3=ib
      ii4=i
      kk3=ib
      kk4=i
      do 570 j=1,nocc
      jj=idorb(j)
      if (ijm.ne.imul(jj,ibb)) go to 570
      epsij=epsil(j)
      epsd=epsc+epsij
      epse=-epsb-epsij
      ll1=ia
      ll2=j
      i1=ia
      i2=j
      i3=ii3
      i4=ii4
      j1=jj1
      j2=jj2
      j3=ib
      j4=j
      if (i1.gt.i3) go to 170
      if (i1.eq.i3) go to 150
      it1=i1
      it2=i2
      i1=i3
      i2=i4
      i3=it1
      i4=it2
      it1=j1
      it2=j2
      j1=j3
      j2=j4
      j3=it1
      j4=it2
      go to 170
  150 if (i2.ge.i4) go to 160
      it1=i2
      i2=i4
      i4=it1
      go to 170
  160 it1=j2
      j2=j4
      j4=it1
  170 continue
      ind1=indv(i1,i2,i3,i4)
      ind2=indv(j1,j2,j3,j4)
      pa=v(ind1)
      pb=v(ind2)
      pi=(pa+pa-pb)/epsd
      do 250 k=1,nocc
      kk=idorb(k)
      if (imul(jj,kk).ne.1) go to 250
      if (imul(kk,ibb).ne.ijm) go to 250
      epsf=epsc+epsil(k)
      aa=a(ind(k,j))
      k1=ia
      k2=k
      k3=kk3
      k4=kk4
      if (k1.gt.k3) go to 210
      if (k1.eq.k3) go to 200
      it1=k1
      it2=k2
      k1=k3
      k2=k4
      k3=it1
      k4=it2
      go to 210
  200 if (k2.ge.k4) go to 210
      it1=k2
      k2=k4
      k4=it1
  210 continue
      ind1=indv(k1,k2,k3,k4)
      a1=a1+pi*aa*v(ind1)/epsf
  250 continue
      do 300 ic=nstart,ntot
      icc=idorb(ic)
      if (imul(ibb,icc).ne.1) go to 300
      if (imul(jj,icc).ne.ijm) go to 300
      epsg=epse+epsil(ic)
      aa=a(ind(ib,ic))
      l1=ll1
      l2=ll2
      l3=ic
      l4=i
      if (l1.gt.l3) go to 280
      if (l1.eq.l3) go to 270
      it1=l1
      it2=l2
      l1=l3
      l2=l4
      l3=it1
      l4=it2
      go to 280
  270 if (l2.ge.l4) go to 280
      it1=l2
      l2=l4
      l4=it1
  280 continue
      ind1=indv(l1,l2,l3,l4)
      a2=a2+pi*aa*v(ind1)/epsg
  300 continue
      do 350 ic=nstart,ntot
      icc=idorb(ic)
      if (imul(jj,icc).ne.1) go to 350
      if (imul(ibb,icc).ne.ijm) go to 350
      epsf=epsij-epsil(ic)
      aa=a(ind(ic,j))
      k1=ia
      k2=ic
      if (ia.ge.ic) go to 320
      k1=ic
      k2=ia
  320 continue
      k3=kk3
      k4=kk4
      if (k1.gt.k3) go to 340
      if (k1.eq.k3) go to 330
      it1=k1
      it2=k2
      k1=k3
      k2=k4
      k3=it1
      k4=it2
      go to 340
  330 if (k2.ge.k4) go to 340
      it1=k2
      k2=k4
      k4=it1
  340 continue
      ind1=indv(k1,k2,k3,k4)
      a3=a3+pi*aa*v(ind1)/epsf
  350 continue
      do 400 k=1,nocc
      kk=idorb(k)
      if (imul(kk,ibb).ne.1) go to 400
      if (imul(kk,jj).ne.ijm) go to 400
      epsg=epsib-epsil(k)
      aa=a(ind(ib,k))
      l1=ll1
      l2=ll2
      l3=i
      l4=k
      if (i.ge.k) go to 360
      l3=k
      l4=i
  360 continue
      if (l1.gt.l3) go to 380
      if (l1.eq.l3) go to 370
      it1=l1
      it2=l2
      l1=l3
      l2=l4
      l3=it1
      l4=it2
      go to 380
  370 if (l2.ge.l4) go to 380
      it1=l2
      l2=l4
      l4=it1
  380 continue
      ind1=indv(l1,l2,l3,l4)
      a5=a5+pi*aa*v(ind1)/epsg
  400 continue
  570 continue
  580 continue
  590 continue
  600 continue
      a1=-a1
      a2=-a2
      a3=a3+a3
      a5=a5+a5
      return
      end
      subroutine gc1c2(v,indvec,jndvec,kndvec,lndvec)
      implicit real*8  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      common/miscop/epsil(maxorb),wmega,conver,nocc,nunocc,nstart,ntot,
     *              nsym,imul(8,8),idorb(maxorb),
     *              nsecd,insecd(maxorb),nthd,inthd(maxorb),m,mm
      common /blkin/  a1,a2,a3,a5,c1,c2,c4,c6,d1,d3,d5,d6,dc1,dc2,dc4,
     *               dc6,dd1,dd3,dd5,dd6
      dimension v(*)
      dimension indvec(*),jndvec(*),kndvec(*),lndvec(*)
      indv (i,j,k,l) =  lndvec(l) +
     *                  indvec(k+jndvec(j+kndvec(i)))
      c1=0.0d0
      c2=0.0d0
      dc1=0.0d0
      dc2=0.0d0
      do 900 i=1,nocc
      ii=idorb(i)
      iim=imul(mm,ii)
      epsoi=wmega+epsil(i)
      do 890 ia=nstart,ntot
      iaa=idorb(ia)
      epsia=epsil(ia)
      ii1=ia
      ii2=m
      if (ia.ge.m) go to 10
      ii1=m
      ii2=ia
   10 continue
      jj3=ia
      jj4=i
      do 880 ib=nstart,ntot
      ibb=idorb(ib)
      iab=imul(iaa,ibb)
      if (iim.ne.iab) go to 880
      epsab=-epsia-epsil(ib)
      epsa=epsoi+epsab
      rz=1.0d0/epsa
      i1=ii1
      i2=ii2
      i3=ib
      i4=i
      if (i1.gt.i3) go to 40
      if (i1.eq.i3) go to 30
      it1=i1
      it2=i2
      i1=i3
      i2=i4
      i3=it1
      i4=it2
      go to 40
   30 if (i2.ge.i4) go to 40
      it1=i2
      i2=i4
      i4=it1
   40 continue
      j3=jj3
      j4=jj4
      j1=ib
      j2=m
      if (ib.ge.m) go to 50
      j1=m
      j2=ib
   50 continue
      if (j1.gt.j3) go to 80
      if (j1.eq.j3) go to 70
      it1=j1
      it2=j2
      j1=j3
      j2=j4
      j3=it1
      j4=it2
      go to 80
   70 if (j2.ge.j4) go to 80
      it1=j2
      j2=j4
      j4=it1
   80 continue
      ind1=indv(i1,i2,i3,i4)
      pa=v(ind1)
      ind1=indv(j1,j2,j3,j4)
      pb=v(ind1)
      pi=(pa+pa-pb)*rz
      pii=pi*rz
      do 400 ic=nstart,ntot
      icc=idorb(ic)
      epsic=epsil(ic)
      kk1=ia
      kk2=ic
      if (ia.ge.ic) go to 110
      kk1=ic
      kk2=ia
  110 continue
      ll1=ic
      ll2=m
      if (ic.ge.m) go to 120
      ll1=m
      ll2=ic
  120 continue
      do 390 id=nstart,ntot
      idd=idorb(id)
      icd=imul(icc,idd)
      if (iim.ne.icd) go to 390
      epscd=-epsic-epsil(id)
      epsb=epsoi+epscd
      am=epsa+epsb
      rzb=1.0d0/epsb
      k1=kk1
      k2=kk2
      k3=ib
      k4=id
      if (ib.ge.id) go to 130
      k3=id
      k4=ib
  130 continue
      if (k1.gt.k3) go to 160
      if (k1.eq.k3) go to 150
      it1=k1
      it2=k2
      k1=k3
      k2=k4
      k3=it1
      k4=it2
      go to 160
  150 if (k2.ge.k4) go to 160
      it1=k2
      k2=k4
      k4=it1
  160 continue
      ind1=indv(k1,k2,k3,k4)
      qa=v(ind1)
      l1=ll1
      l2=ll2
      l3=id
      l4=i
      if (l1.gt.l3) go to 200
      if (l1.eq.l3) go to 180
      it1=l1
      it2=l2
      l1=l3
      l2=l4
      l3=it1
      l4=it2
      go to 200
  180 if (l2.ge.l4) go to 200
      it1=l2
      l2=l4
      l4=it1
  200 continue
      ind1=indv(l1,l2,l3,l4)
      qb=v(ind1)
      q=qa*qb*rzb
      c1=c1+pi*q
      dc1=dc1+pii*q*rzb*am
  390 continue
  400 continue
      do 700 j=1,nocc
      jj=idorb(j)
      epsb=epsab+epsil(j)
      kk1=ia
      kk2=j
      ll1=m
      ll2=j
      if (m.ge.j) go to 510
      ll1=j
      ll2=m
  510 continue
      do 690 k=1,nocc
      kk=idorb(k)
      icd=imul(jj,kk)
      if (iim.ne.icd) go to 690
      epsc=epsb+epsil(k)
      k1=kk1
      k2=kk2
      k3=ib
      k4=k
      if (k1.gt.k3) go to 540
      if (k1.eq.k3) go to 530
      it1=k1
      it2=k2
      k1=k3
      k2=k4
      k3=it1
      k4=it2
      go to 540
  530 if (k2.ge.k4) go to 540
      it1=k2
      k2=k4
      k4=it1
  540 continue
      ind1=indv(k1,k2,k3,k4)
      qa=v(ind1)
      l1=ll1
      l2=ll2
      l3=i
      l4=k
      if (i.ge.k) go to 550
      l3=k
      l4=i
  550 continue
      if (l1.gt.l3) go to 580
      if (l1.eq.l3) go to 570
      it1=l1
      it2=l2
      l1=l3
      l2=l4
      l3=it1
      l4=it2
      go to 580
  570 if (l2.ge.l4) go to 580
      it1=l2
      l2=l4
      l4=it1
  580 continue
      ind1=indv(l1,l2,l3,l4)
      q=qa*v(ind1)/epsc
      c2=c2+pi*q
      dc2=dc2+pii*q
  690 continue
  700 continue
  880 continue
  890 continue
  900 continue
      c2=c2+c2
      dc1=-dc1
      dc2=-dc2-dc2
      return
      end
      subroutine gc4c6(v,indvec,jndvec,kndvec,lndvec)
      implicit real*8  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      common/miscop/epsil(maxorb),wmega,conver,nocc,nunocc,nstart,ntot,
     *              nsym,imul(8,8),idorb(maxorb),
     *              nsecd,insecd(maxorb),nthd,inthd(maxorb),m,mm
      common /blkin/  a1,a2,a3,a5,c1,c2,c4,c6,d1,d3,d5,d6,dc1,dc2,dc4,
     *               dc6,dd1,dd3,dd5,dd6
      dimension v(*)
      dimension indvec(*),jndvec(*),kndvec(*),lndvec(*)
      indv (i,j,k,l) =  lndvec(l) +
     *                  indvec(k+jndvec(j+kndvec(i)))
      c4=0.0d0
      c6=0.0d0
      dc4=0.0d0
      dc6=0.0d0
      do 400 ia=nstart,ntot
      iaa=idorb(ia)
      iam=imul(mm,iaa)
      epsa=wmega+epsil(ia)
      do 390 i=1,nocc
      ii=idorb(i)
      epsii=epsil(i)
      ii3=m
      ii4=i
      if (m.ge.i) go to 10
      ii3=i
      ii4=m
   10 continue
      jj1=ia
      jj2=i
      do 380 j=1,nocc
      jj=idorb(j)
      iij=imul(ii,jj)
      if (iam.ne.iij) go to 380
      epsab=epsii+epsil(j)
      epsb=epsa-epsab
      rz=1.0d0/epsb
      i3=ii3
      i4=ii4
      i1=ia
      i2=j
      if (i1.gt.i3) go to 40
      if (i1.eq.i3) go to 30
      it1=i1
      it2=i2
      i1=i3
      i2=i4
      i3=it1
      i4=it2
      go to 40
   30 if (i2.ge.i4) go to 40
      it1=i2
      i2=i4
      i4=it1
   40 continue
      ind1=indv(i1,i2,i3,i4)
      pa=v(ind1)
      j1=jj1
      j2=jj2
      j3=m
      j4=j
      if (m.ge.j) go to 50
      j3=j
      j4=m
   50 continue
      if (j1.gt.j3) go to 80
      if (j1.eq.j3) go to 70
      it1=j1
      it2=j2
      j1=j3
      j2=j4
      j3=it1
      j4=it2
      go to 80
   70 if (j2.ge.j4) go to 80
      it1=j2
      j2=j4
      j4=it1
   80 ind1=indv(j1,j2,j3,j4)
      pb=v(ind1)
      pi=(pa+pa-pb)*rz
      pii=pi*rz
      do 200 ib=nstart,ntot
      ibb=idorb(ib)
      epsc=epsab-epsil(ib)
      kk1=ib
      kk2=i
      ll3=ib
      ll4=m
      if (ib.ge.m) go to 110
      ll3=m
      ll4=ib
  110 continue
      do 190 ic=nstart,ntot
      icc=idorb(ic)
      if (iam.ne.imul(ibb,icc)) go to 190
      epsd=epsc-epsil(ic)
      k3=ic
      k4=j
      k1=kk1
      k2=kk2
      if (k1.gt.k3) go to 140
      if (k1.eq.k3) go to 130
      it1=k1
      it2=k2
      k1=k3
      k2=k4
      k3=it1
      k4=it2
      go to 140
  130 if (k2.ge.k4) go to 140
      it1=k2
      k2=k4
      k4=it1
  140 continue
      ind1=indv(k1,k2,k3,k4)
      qa=v(ind1)
      l3=ll3
      l4=ll4
      l1=ia
      l2=ic
      if (ia.ge.ic) go to 150
      l1=ic
      l2=ia
  150 continue
      if (l1.gt.l3) go to 180
      if (l1.eq.l3) go to 170
      it1=l1
      it2=l2
      l1=l3
      l2=l4
      l3=it1
      l4=it2
      go to 180
  170 if (l2.ge.l4) go to 180
      it1=l2
      l2=l4
      l4=it1
  180 continue
      ind1=indv(l1,l2,l3,l4)
      q=qa*v(ind1)/epsd
      c4=c4+pi*q
      dc4=dc4+pii*q
  190 continue
  200 continue
      do 300 k=1,nocc
      kk=idorb(k)
      epse=epsa-epsil(k)
      mm3=m
      mm4=k
      if (m.ge.k) go to 210
      mm3=k
      mm4=m
  210 continue
      nn1=i
      nn2=k
      if (i.ge.k) go to 220
      nn1=k
      nn2=i
  220 continue
      do 290 l=1,nocc
      ll=idorb(l)
      if (iam.ne.imul(kk,ll)) go to 290
      epsf=epse-epsil(l)
      rzb=1.0d0/epsf
      m1=ia
      m2=l
      m3=mm3
      m4=mm4
      if (m1.gt.m3) go to 240
      if (m1.eq.m3) go to 230
      it1=m1
      it2=m2
      m1=m3
      m2=m4
      m3=it1
      m4=it2
      go to 240
  230 if (m2.ge.m4) go to 240
      it1=m2
      m2=m4
      m4=it1
  240 continue
      ind1=indv(m1,m2,m2,m4)
      ind1=indv(m1,m2,m3,m4)
      ra=v(ind1)
      n1=nn1
      n2=nn2
      n3=j
      n4=l
      if (j.ge.l) go to 250
      n3=l
      n4=j
  250 continue
      if (n1.gt.n3) go to 270
      if (n1.eq.n3) go to 260
      it1=n1
      it2=n2
      n1=n3
      n2=n4
      n3=it1
      n4=it2
      go to 270
  260 if (n2.ge.n4) go to 270
      it1=n2
      n2=n4
      n4=it1
  270 continue
      ind1=indv(n1,n2,n3,n4)
      rb=v(ind1)
      ri=ra*rb*rzb
      am=epsb+epsf
      rii=ri*rzb
      c6=c6+pi*ri
      dc6=dc6+pii*rii*am
  290 continue
  300 continue
  380 continue
  390 continue
  400 continue
      c4=c4+c4
      c6=-c6
      dc4=-dc4-dc4
      return
      end
      subroutine grphda (v,indvec,jndvec,kndvec,lndvec)
      implicit real*8  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      common/miscop/epsil(maxorb),wmega,conver,nocc,nunocc,nstart,ntot,
     *              nsym,imul(8,8),idorb(maxorb),
     *              nsecd,insecd(maxorb),nthd,inthd(maxorb),m,mm
      common /blkin/  a1,a2,a3,a5,c1,c2,c4,c6,d1,d3,d5,d6,dc1,dc2,dc4,
     *               dc6,dd1,dd3,dd5,dd6
      dimension v(*)
      dimension indvec(*),jndvec(*),kndvec(*),lndvec(*)
      indv (i,j,k,l) =  lndvec(l) +
     *                  indvec(k+jndvec(j+kndvec(i)))
      d1=0.0d0
      d3=0.0d0
      dd1=0.0d0
      dd3=0.0d0
      do 900 ia=nstart,ntot
      iaa=idorb(ia)
      ijm=imul(iaa,mm)
      epsia=epsil(ia)
      epsoi=wmega-epsia
      jj1=ia
      jj2=m
      mm1=ia
      mm2=m
      if (ia.ge.m) go to 10
      jj1=m
      jj2=ia
      mm1=m
      mm2=ia
   10 continue
      do 890 ib=nstart,ntot
      ibb=idorb(ib)
      epsib=epsil(ib)
      epsa=epsoi-epsib
      ii3=ib
      ii4=m
      if (ib.ge.m) go to 20
      ii3=m
      ii4=ib
   20 continue
      do 880 i=1,nocc
      ii=idorb(i)
      ikm=imul(ii,ibb)
      if (ijm.ne.ikm) go to 880
      epsii=epsil(i)
      epsb=epsa+epsii
      rz=1.0d0/epsb
      epscd=epsii-epsib
      i1=ia
      i2=i
      j3=ib
      j4=i
      ll1=ib
      ll2=i
      i3=ii3
      i4=ii4
      if (i1.gt.i3) go to 40
      if (i1.eq.i3) go to 30
      it1=i1
      it2=i2
      i1=i3
      i2=i4
      i3=it1
      i4=it2
      go to 40
   30 if (i2.ge.i4) go to 40
      it1=i2
      i2=i4
      i4=it1
   40 continue
      ind1=indv(i1,i2,i3,i4)
      pa=v(ind1)
      j1=jj1
      j2=jj2
      if (j1.gt.j3) go to 60
      if (j1.eq.j3) go to 50
      it1=j1
      it2=j2
      j1=j3
      j2=j4
      j3=it1
      j4=it2
      go to 60
   50 if (j2.ge.j4) go to 60
      it1=j2
      j2=j4
      j4=it1
   60 continue
      ind1=indv(j1,j2,j3,j4)
      pb=v(ind1)
      pai=(pa-pb-pb)*rz
      paii=pai*rz
      pbi=(pb-pa-pa)*rz
      pbii=pbi*rz
      do 850 j=1,nocc
      jj=idorb(j)
      epsij=epsil(j)
      epsc=epsoi+epsij
      epsf=epscd+epsij
      kk1=ia
      kk2=j
      nn3=i
      nn4=j
      if (i.ge.j) go to 110
      nn3=j
      nn4=i
  110 continue
      irr1=ib
      irr2=j
      iss3=m
      iss4=j
      if (m.ge.j) go to 120
      iss3=j
      iss4=m
  120 continue
      do 840 ic=nstart,ntot
      icc=idorb(ic)
      ijc=imul(jj,icc)
      if (ijm.ne.ijc) go to 840
      epsic=epsil(ic)
      epsg=epsc-epsic
      rza1=1.0d0/epsg
      epsh=epsf-epsic
      rzb1=1.0d0/epsh
      am=epsb+epsg
      k1=kk1
      k2=kk2
      k3=ic
      k4=m
      if (ic.ge.m) go to 130
      k3=m
      k4=ic
  130 continue
      if (k1.gt.k3) go to 150
      if (k1.eq.k3) go to 140
      it1=k1
      it2=k2
      k1=k3
      k2=k4
      k3=it1
      k4=it2
      go to 150
  140 if (k2.ge.k4) go to 150
      it1=k2
      k2=k4
      k4=it1
  150 continue
      ind1=indv(k1,k2,k3,k4)
      qa=v(ind1)
      l1=ll1
      l2=ll2
      l3=ic
      l4=j
      if (l1.gt.l3) go to 170
      if (l1.eq.l3) go to 160
      it1=l1
      it2=l2
      l1=l3
      l2=l4
      l3=it1
      l4=it2
      go to 170
  160 if (l2.ge.l4) go to 170
      it1=l2
      l2=l4
      l4=it1
  170 continue
      ind1=indv(l1,l2,l3,l4)
      qb=v(ind1)
      m1=mm1
      m2=mm2
      m3=ic
      m4=j
      if (m1.gt.m3) go to 190
      if (m1.eq.m3) go to 180
      it1=m1
      it2=m2
      m1=m3
      m2=m4
      m3=it1
      m4=it2
      go to 190
  180 if (m2.ge.m4) go to 190
      it1=m2
      m2=m4
      m4=it1
  190 continue
      ind1=indv(m1,m2,m3,m4)
      qc=v(ind1)
      n1=ib
      n2=ic
      if (ib.ge.ic) go to 200
      n1=ic
      n2=ib
  200 continue
      n3=nn3
      n4=nn4
      ind1=indv(n1,n2,n3,n4)
      qd=v(ind1)
      ir1=irr1
      ir2=irr2
      ir3=ic
      ir4=i
      if (ir1.gt.ir3) go to 220
      if (ir1.eq.ir3) go to 210
      it1=ir1
      it2=ir2
      ir1=ir3
      ir2=ir4
      ir3=it1
      ir4=it2
      go to 220
  210 if (ir2.ge.ir4) go to 220
      it1=ir2
      ir2=ir4
      ir4=it1
  220 continue
      ind1=indv(ir1,ir2,ir3,ir4)
      qe=v(ind1)
      is3=iss3
      is4=iss4
      is1=ia
      is2=ic
      if (ia.ge.ic) go to 230
      is1=ic
      is2=ia
  230 continue
      if (is1.ge.is3) go to 250
      it1=is1
      it2=is2
      is1=is3
      is2=is4
      is3=it1
      is4=it2
  250 continue
      ind1=indv(is1,is2,is3,is4)
      qf=v(ind1)
      qai=(qa*qb+qc*(qd-qb-qb))*rza1
      pmul=rza1*am
      qaii=qai*pmul
      qbi=(qc*qe+qb*(qf-qc-qc))*rzb1
      qci=qa*qd*rza1
      qcii=qci*pmul
      qdi=qe*qf*rzb1
      d1=d1+pai*qai+pbi*qci
      d3=d3+pai*qbi+pbi*qdi
      dd1=dd1+paii*qaii+pbii*qcii
      dd3=dd3+paii*qbi+pbii*qdi
  840 continue
  850 continue
  880 continue
  890 continue
  900 continue
      d3=d3+d3
      dd1=-dd1
      dd3=-dd3-dd3
      return
      end
      subroutine grphdb(v,indvec,jndvec,kndvec,lndvec)
      implicit real*8  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      common/miscop/epsil(maxorb),wmega,conver,nocc,nunocc,nstart,ntot,
     *              nsym,imul(8,8),idorb(maxorb),
     *              nsecd,insecd(maxorb),nthd,inthd(maxorb),m,mm
      common /blkin/  a1,a2,a3,a5,c1,c2,c4,c6,d1,d3,d5,d6,dc1,dc2,dc4,
     *               dc6,dd1,dd3,dd5,dd6
      dimension v(*)
      dimension indvec(*),jndvec(*),kndvec(*),lndvec(*)
      indv (i,j,k,l) =  lndvec(l) +
     *                  indvec(k+jndvec(j+kndvec(i)))
      d5=0.0d0
      d6=0.0d0
      dd5=0.0d0
      dd6=0.0d0
      do 900 j=1,nocc
      jj=idorb(j)
      ijm=imul(mm,jj)
      epsij=epsil(j)
      epsoj=wmega-epsij
      jj3=m
      jj4=j
      mm3=m
      mm4=j
      if (m.ge.j) go to 10
      jj3=j
      jj4=m
      mm3=j
      mm4=m
   10 continue
      do 890 i=1,nocc
      ii=idorb(i)
      epsii=epsil(i)
      epsa=epsoj-epsii
      ii3=m
      ii4=i
      if (m.ge.i) go to 20
      ii3=i
      ii4=m
   20 continue
      do 880 ia=nstart,ntot
      iaa=idorb(ia)
      ikm=imul(ii,iaa)
      if (ijm.ne.ikm) go to 880
      epsia=epsil(ia)
      epsb=epsa+epsia
      rz=1.0d0/epsb
      epsc=epsii-epsia
      i1=ia
      i2=j
      i3=ii3
      i4=ii4
      if (i1.gt.i3) go to 40
      if (i1.eq.i3) go to 30
      it1=i1
      it2=i2
      i1=i3
      i2=i4
      i3=it1
      i4=it2
      go to 40
   30 if (i2.ge.i4) go to 40
      it1=i2
      i2=i4
      i4=it1
   40 continue
      ind1=indv(i1,i2,i3,i4)
      pa=v(ind1)
      j1=ia
      j2=i
      j3=jj3
      j4=jj4
      if (j1.gt.j3) go to 60
      if (j1.eq.j3) go to 50
      it1=j1
      it2=j2
      j1=j3
      j2=j4
      j3=it1
      j4=it2
      go to 60
   50 if (j2.ge.j4) go to 60
      it1=j2
      j2=j4
      j4=it1
   60 continue
      ind1=indv(j1,j2,j3,j4)
      pb=v(ind1)
      ll1=ia
      ll2=i
      pai=(pa-pb-pb)*rz
      paii=pai*rz
      pbi=(pb-pa-pa)*rz
      pbii=pbi*rz
      do 850 k=1,nocc
      kk=idorb(k)
      epsik=epsil(k)
      epse=epsc+epsik
      epsh=epsoj-epsik
      kk3=m
      kk4=k
      if (m.ge.k) go to 70
      kk3=k
      kk4=m
   70 continue
      nn3=i
      nn4=k
      if (i.ge.k) go to 80
      nn3=k
      nn4=i
   80 continue
      iss3=j
      iss4=k
      if (j.ge.k) go to 90
      iss3=k
      iss4=j
   90 continue
      irr1=ia
      irr2=k
      do 840 ib=nstart,ntot
      ibb=idorb(ib)
      ikb=imul(kk,ibb)
      if (ijm.ne.ikb) go to 840
      epsib=epsil(ib)
      epsi=epsh+epsib
      rza1=1.0d0/epsi
      epsj=epse-epsib
      rzb1=1.0d0/epsj
      am=epsb+epsi
      k1=ib
      k2=j
      k3=kk3
      k4=kk4
      if (k1.gt.k3) go to 150
      if (k1.eq.k3) go to 140
      it1=k1
      it2=k2
      k1=k3
      k2=k4
      k3=it1
      k4=it2
      go to 150
  140 if (k2.ge.k4) go to 150
      it1=k2
      k2=k4
      k4=it1
  150 continue
      ind1=indv(k1,k2,k3,k4)
      qa=v(ind1)
      l1=ll1
      l2=ll2
      l3=ib
      l4=k
      if (l1.gt.l3) go to 170
      if (l1.eq.l3) go to 160
      it1=l1
      it2=l2
      l1=l3
      l2=l4
      l3=it1
      l4=it2
      go to 170
  160 if (l2.ge.l4) go to 170
      it1=l2
      l2=l4
      l4=it1
  170 continue
      ind1=indv(l1,l2,l3,l4)
      qb=v(ind1)
      m1=ib
      m2=k
      m3=mm3
      m4=mm4
      if (m1.gt.m3) go to 190
      if (m1.eq.m3) go to 180
      it1=m1
      it2=m2
      m1=m3
      m2=m4
      m3=it1
      m4=it2
      go to 190
  180 if (m2.ge.m4) go to 190
      it1=m2
      m2=m4
      m4=it1
  190 continue
      ind1=indv(m1,m2,m3,m4)
      qc=v(ind1)
      n1=ia
      n2=ib
      n3=nn3
      n4=nn4
      if (ia.ge.ib) go to 200
      n1=ib
      n2=ia
  200 continue
      ind1=indv(n1,n2,n3,n4)
      qd=v(ind1)
      ir1=irr1
      ir2=irr2
      ir3=ib
      ir4=i
      if (ir1.gt.ir3) go to 240
      if (ir1.eq.ir3) go to 230
      it1=ir1
      it2=ir2
      ir1=ir3
      ir2=ir4
      ir3=it1
      ir4=it2
      go to 240
  230 if (ir2.ge.ir4)    go to 240
      it1=ir2
      ir2=ir4
      ir4=it1
  240 continue
      ind1=indv(ir1,ir2,ir3,ir4)
      qe=v(ind1)
      is1=ib
      is2=m
      if (ib.ge.m) go to 250
      is1=m
      is2=ib
  250 continue
      is3=iss3
      is4=iss4
      ind1=indv(is1,is2,is3,is4)
      qf=v(ind1)
      qai=(qa*qb+ qc*(qd-qb-qb))*rza1
      pmul=rza1*am
      qaii=qai*pmul
      qbi=(qc*qe+qb*(qf-qc-qc))*rzb1
      qci=qa*qd*rza1
      qcii=qci*pmul
      qdi=qe*qf*rzb1
      d6=d6+pai*qai+pbi*qci
      d5=d5+pai*qbi+pbi*qdi
      dd6=dd6+paii*qaii+pbii*qcii
      dd5=dd5+paii*qbi+pbii*qdi
  840 continue
  850 continue
  880 continue
  890 continue
  900 continue
      d5=d5+d5
      d6=-d6
      dd5=-dd5-dd5
      return
      end
c
c symvec now moved to util1.m
c
      subroutine inip
      implicit real*8  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
c     mxtoken - the maximum number of tokens on an input line
c               for space separated tokens: 132/2 plus 1 to allow
c               for look ahead (132 is the line-length defined in
c               common/workc)
c
      integer mxtoken
      parameter(mxtoken=67)
      integer jrec, jump, istrt, inumb, iwidth, nend, nstart
      integer nline, noline, jwidth, nerr
      logical oswit, oterm, oflush
      common/work/jrec,jump,istrt(mxtoken),inumb(mxtoken),iwidth,
     + nend(mxtoken),nstart(mxtoken), nline,noline,jwidth,nerr,
     + oswit,oterm,oflush
c
      common/junkc/ylab(26),ztype(831),zoptio(17),ztagg(100)
      common/data1/vlist(400),np(206),np1(maxorb,2),norbt(6),
     * norbta(maxorb,2),
     * norb2,nseco(maxorb),norb3,nthird(maxorb),norbg,norbga(maxorb),
     * nscrp,range(2),alimit,blimit,isymm,nx,lsym,nsymg(8),iuse(8,6),
     * idorb(100),ncont,nmax(10)
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
      dimension yopt(5),ygf(3),ytda(9),ysymm(2),norbaa(2)
      dimension itest(10),iprod(10),ytest(3)
c
      equivalence(norba,norb2),(norbaa(1),nseco(1))
c
      data ysymm/'a','b'/
      data zgf,ztda/'gf','tda'/
      data ygf/'seco','thir','acti'/
      data ytda/'dege','rang','limi','band','acti','on',
     *          'off','symm','nmax'/
      data yopt/'psa','adia','quar','sele','diag'/
      data oyes,ono/.true.,.false./
      data yoff,ytest/'off',' ','prin','test'/
      data zend/'end'/
c
      data itest/
     *30,200,20,30,4000,3200,200,20,6,200/
      data iprod/
     *50,600,20,30,176000,60000,3000,100,6,600/
c
      if(zruntp.eq.ztda)go to 1000
      if(zruntp.ne.zgf )
     *call caserr('faulty directive ordering')
c
      call inpa4(ytext)
      if(ytext.eq.ytest(2))oprint(54)=oyes
40030 call input
      call inpa4(ytext)
      i=locatc(ygf,3,ytext)+1
      go to (40090,40040,40050,40060),i
40040 call inpa4(ytext)
      if(ytext.ne.yoff)go to 40041
      norb2=0
      go to 40030
40041 jrec=jrec-1
      call active(norb2,nseco,ono,ono)
      go to 40030
40050 call inpa4(ytext)
      if(ytext.ne.yoff)go to 40051
      norb3=0
      go to 40030
40051 jrec=jrec-1
      call active(norb3,nthird,ono,ono)
      go to 40030
40060 call active(norbg,norbga,ono,oyes)
      go to 40030
c
c     tda input
c
 1000 do 6000 i=1,10
 6000 nmax(i)=iprod(i)
 6040 call inpa4(ytext)
      i=locatc(ytest,3,ytext)
      go to (6010,6020,6030),i
 6020 oprint(55)=oyes
      go to 6040
 6030 do 6050 i=1,10
 6050 nmax(i)=itest(i)
      go to 6040
 6010 zoptio(1)='selectio'
      zoptio(2)='closed  '
      zoptio(3)='open    '
      zoptio(4)='correl  '
      zoptio(5)='closed  '
      zoptio(6)='selectio'
      zoptio(7)='degenera'
      zoptio(8)='no degen'
      zoptio(9)='no quart'
      zoptio(10)='quart   '
      zoptio(11)='psa     '
      zoptio(12)='psa     '
      zoptio(13)='adiagram'
      zoptio(14)='no diago'
      zoptio(15)='        '
      zoptio(16)='        '
      zoptio(17)='        '
c
      ncont=0
      nx=1
      alimit=-10000.0d0
      blimit= 10000.0d0
      range(1)=0.0d0
      range(2)=0.0d0
      do 1 i=1,100
 1    ztagg(i)=' '
      isymm=1
      call setsto(100,0,idorb)
      lsym=0
      do 2 i=1,8
      nsymg(i)=0
      do 2 j=1,6
 2    iuse(i,j)=0
c
50000 call input
      call inpa4(ytext)
      k=locatc(ytda,9,ytext)+1
      go to (40080,100,200,300,400,500,600,700,800,900),k
 100  zoptio(8)='degenera'
      call input
      lsym=jump
      if(lsym.gt.0)go to 104
 105  call caserr('invalid degeneracy data in tda input')
 104  do 101 i=1,lsym
 101  call inpi(nsymg(i))
      do 102 i=1,lsym
      jend=nsymg(i)
      call input
      do 103 j=1,jend
 103  call inpi(iuse(i,j))
 102  continue
      j=0
 109  call input
      do 110 i=1,jump
      call inpa(ztext)
      if(ztext.eq.zend)go to 111
       jrec=jrec-1
      j=j+1
110   call inpi(idorb(j))
      go to 109
 111  j=0
 106  call input
      do 107 i=1,jump
      call inpa(ztext)
      if(ztext.eq.zend)go to 108
      j=j+1
 107  ztagg(j)=ztext
      go to 106
 108  if(j.le.0.or.j.gt.100)go to 105
      go to 50000
 200  call inpf(range(1))
      call inpf(range(2))
      go to 50000
 300  call inpf(alimit)
      call inpf(blimit)
      go to 50000
 400  call active(norba,norbaa,ono,ono)
      go to 50000
 500  call active(norbg,norbga,ono,oyes)
      go to 50000
 600   call inpa4(ytext)
       j=locatc(yopt,5,ytext)+ 1
       go to (50000,600,460,470,600,480),j
 460   call inpi(nx)
       go to 600
 470   zoptio(9)=zoptio(10)
       go to 600
 480   zoptio(14)='diagonal'
       goto 600
 700   call inpa4(ytext)
       j=locatc(yopt,4,ytext)+ 1
       go to (50000,350,360,700,380,390),j
 350   zoptio(12)='no psa'
       go to 700
 360   zoptio(13)='no adiag'
       go to 700
 380   zoptio(1)='no selec'
       go to 700
 390   zoptio(14)='no diago'
       goto 700
 800  call inpa4(ytext)
      isymm=locatc(ysymm,2,ytext)
      go to 50000
 900  do 991 i = 1,10
         call inpi(kkk)
         if(kkk.ge.0)nmax(i)=kkk
         write(iwr,*)'set nmax',i,kkk
991   continue
      goto 50000
40080 continue
c
c check nmax settings for overflow
c
      overfl=.false.
      if(nmax(7)+nmax(2).gt.mxtda1.or.
     &     nmax(7)+nmax(10).gt.mxtda1)then
         write(iwr,*)
     &        'dimension overflow - increase mxtda1 and recompile'
         overfl=.true.
      endif
      if(nmax(1).gt.mxtda2)then
         write(iwr,*)
     &        'dimension overflow - increase mxtda2 and recompile'
         overfl=.true.
      endif
      if(nmax(3).gt.mxtda3)then
         write(iwr,*)
     &        'dimension overflow - increase mxtda3 and recompile'
         overfl=.true.
      endif
      if(nmax(2).gt.mxtda4.or.
     &     nmax(10).gt.mxtda4)then
         write(iwr,*)
     &        'dimension overflow - increase mxtda4 and recompile'
         overfl=.true.
      endif
      if(overfl)call caserr('tda dimensioning problem')
40090 jrec=jrec-1
      return
      end
      subroutine ver_gff(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/gff.m,v $
     +     "/
      data revision /"$Revision: 6176 $"/
      data date /"$Date: 2010-08-10 16:49:47 +0200 (Tue, 10 Aug 2010) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
