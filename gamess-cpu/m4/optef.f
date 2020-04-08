**==checkd.f
      subroutine checkd(d,stepmx,dd,skal,nvar,iw)
c
      implicit real*8  (a-h,p-w) , integer  (i-n) , logical (o)
      implicit character*8 (z) , character*4 (y) , character*1 (x)
      dimension d(*)
c ...
c                   get current stepsize
c ...
      skal = 1.0d0
      dd = dsqrt(ddot(nvar,d,1,d,1))
      if (dd.gt.stepmx) then
         skal = stepmx/dd
         write (iw,6010) skal
         call dvskal(d,d,skal,nvar)
         dd = stepmx
      end if
      write (iw,6020) dd
      return
 6010 format (' calculated step too large. step scaled by ',f9.6)
 6020 format (' step taken. stepsize is ',f9.6)
      end
**==convef.f
      subroutine convef(neg,oconv,tvec,g,rx,d,iw)
c
      implicit real*8  (a-h,p-w) , integer  (i-n) , logical (o)
      implicit character*8 (z) , character*4 (y) , character*1 (x)
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
      dimension tvec(*),g(*),rx(*),d(*)
c
      common/jorgs/ maxj,irecal,opowel,obfgs,obfgx,onrfo,ocut,
     1 oprnt,odump
      common/optef/ jnksp,mode,es,maxstp,ihess,is,negreq,
     1 stpmax,dd,ggrms,gmax
     2 ,drms,dmax,iupdat,maxhes
c
      integer ianz, iz, lbl, lalpha, lbeta, nz, nvar
      real*8 bl, alpha, beta
      common /czmat/ ianz(maxnz),iz(maxnz,4),bl(maxnz),
     +               alpha(maxnz),beta(maxnz),lbl(maxnz),
     +               lalpha(maxnz),lbeta(maxnz),nz,nvar
c
c
      real*8 en, accd
      integer ig, ierec, igrec, iterat, jpoint
      common /cntl2/ en(200),ig(200),ierec(200),igrec(200),
     +               accd(6,200),iterat,jpoint
c
      logical o_var_high,odumdum
      common /csubst/ values(maxvar),intvec(maxvar),fpvec(maxvar),
     1                cmin10(maxvar),cmin20(maxvar),o_var_high,odumdum
      data yno /'no  '/, yes /'yes '/
c ...
c             print out current rx, current gradients, displacements
c             and newx
c ...
      write (iw,6010)
      write (iw,6020)
      do 20 i = 1 , nvar
         tvec(i) = rx(i) + d(i)
         write (iw,6030) i , rx(i) , g(i) , d(i) , tvec(i)
 20   continue
c ...
c              find maximum and rms gradients and displacements
c ...
      fmax = dabs(g(1))
      dxmax = dabs(d(1))
      if (nvar.gt.1) then
         do 30 i = 2 , nvar
            if (dabs(g(i)).gt.fmax) fmax = dabs(g(i))
            if (dabs(d(i)).gt.dxmax) dxmax = dabs(d(i))
 30      continue
      end if
      frms = dsqrt(ddot(nvar,g,1,g,1)/nvar)
      dxrms = dsqrt(dd*dd/nvar)
c ...
c                         check convergence
c ...
      write (iw,6040)
      ycon = yno
      if (fmax.lt.gmax) ycon = yes
      write (iw,6050) fmax , gmax , ycon
      ycon = yno
      if (frms.lt.ggrms) ycon = yes
      write (iw,6060) frms , ggrms , ycon
      ycon = yno
      if (dxmax.lt.dmax) ycon = yes
      write (iw,6070) dxmax , dmax , ycon
      ycon = yno
      if (dxrms.lt.drms) ycon = yes
      write (iw,6080) dxrms , drms , ycon
c ...
c                set convergence flag and test for maximum iteratio
c                exceeded
c ...
      if (fmax.gt.gmax .or. frms.gt.ggrms .or. dxmax.gt.dmax .or.
     +    dxrms.gt.drms) then
         oconv = .false.
         if (iterat.gt.maxstp) then
            write (iw,6090)
            oconv = .true.
         end if
         return
      end if
c ...
c                 converged. print final solution
c ...
      oconv = .true.
      write (iw,6100)
      if (neg.ne.negreq) write (iw,6110)
      return
 6010 format (1x,' current parameter values (internal coordinates)'/)
 6020 format (4x,'i',8x,'x',9x,'gradient',4x,'displacement',7x,
     +        'newx'/' ',3x,'=',8x,'=',9x,'========',4x,'============',
     +        7x,'====')
 6030 format (1x,i4,4(2x,f10.6,2x))
 6040 format (9x,'item',15x,'value',5x,'threshold',2x,'converged?'/' ',
     +        8x,'====',15x,'=====',5x,'=========',2x,'==========')
 6050 format (' maximum force',12x,f8.6,5x,f8.6,5x,a4)
 6060 format (' rms     force',12x,f8.6,5x,f8.6,5x,a4)
 6070 format (' maximum displacement',5x,f8.6,5x,f8.6,5x,a4)
 6080 format (' rms     displacement',5x,f8.6,5x,f8.6,5x,a4,/)
 6090 format (//' ***********************************'/
     +        ' **     optimization stopped      **'/
     +        ' **  maximum iterations exceeded  **'/
     +        ' ***********************************'//)
 6100 format (//' *************************************************'/
     +        ' **  convergence criteria apparently satisfied  **'/
     +        ' *************************************************'//)
 6110 format (1x,'!! warning: the hessian has the wrong number of',
     +        ' negative eigenvalues !!'/)
      end
**==dmatdv.f
      subroutine dmatdv(a,b,c,nvar)
c
      implicit real*8  (a-h,p-w) , integer  (i-n) , logical (o)
      implicit character*8 (z) , character*4 (y) , character*1 (x)
      dimension a(*),b(nvar,nvar),c(*)
      data dzero/0.0d0/
      do 30 i = 1 , nvar
         a(i) = dzero
         do 20 j = 1 , nvar
            a(i) = a(i) + b(i,j)*c(j)
 20      continue
 30   continue
      return
      end
**==dvplv.f
      subroutine dvplv(a,b,c,n)
c
      implicit real*8  (a-h,p-w) , integer  (i-n) , logical (o)
      implicit character*8 (z) , character*4 (y) , character*1 (x)
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
      dimension a(*),b(*),c(*)
      nb2 = n/2
      do 20 i = 1 , n - 1 , 2
         a(i) = b(i) + c(i)
         a(i+1) = b(i+1) + c(i+1)
 20   continue
      if (nb2*2.ne.n) a(n) = b(n) + c(n)
      return
      end
**==dvskal.f
      subroutine dvskal(a,b,scale,n)
c
      implicit real*8  (a-h,p-w) , integer  (i-n) , logical (o)
      implicit character*8 (z) , character*4 (y) , character*1 (x)
      dimension a(*),b(*)
      nb2 = n/2
      do 20 i = 1 , n - 1 , 2
         a(i) = b(i)*scale
         a(i+1) = b(i+1)*scale
 20   continue
      if (nb2*2.ne.n) a(n) = b(n)*scale
      return
      end
**==estime.f
      subroutine estime(eigval,fx,dd,skal,chnge,nvar)
c
      implicit real*8  (a-h,p-w) , integer  (i-n) , logical (o)
      implicit character*8 (z) , character*4 (y) , character*1 (x)
      real*8  lamda,lamda0
      dimension eigval(*),fx(*)
      common/lambda/lamda,lamda0,it
c ...
c         estimate the likely energy change due to the step d
c ...
      chnge = 0.0d0
      do 20 i = 1 , nvar
         dummy = lamda
         if (i.eq.it) dummy = lamda0
         temp = fx(i)*fx(i)*(dummy-eigval(i)/2.0d0)
         temp = temp/((dummy-eigval(i))*(dummy-eigval(i)))
         chnge = chnge + temp
 20   continue
      chnge = chnge/(1+dd*dd)
c ...
c     scale the energy, in case the step was scaled
c     (a scaled energy change is likely to be an overestimate)
c ...
      chnge = chnge*skal
      return
      end
**==fndneg.f
      subroutine fndneg(eigval,neg,nvar)
c
c
      implicit real*8  (a-h,o-z)
      dimension   eigval(*)
      data dzero/0.0d0/
      do 20 i = 1 , nvar
         if (eigval(i).gt.dzero) go to 30
 20   continue
 30   neg = i - 1
      return
      end
**==formd.f
      subroutine formd(u,eigval,fx,nvar,d,vmode,iw)
c
      implicit real*8  (a-h,p-w) , integer  (i-n) , logical (o)
      implicit character*8 (z) , character*4 (y) , character*1 (x)
      real*8  lamda,lamda0,lamda1,lamda2
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
      dimension u(nvar,nvar),eigval(*),fx(*),d(*),vmode(*)
      common/jorgs/ maxj,irecal,opowel,obfgs,obfgx,onrfo,ocut,
     1 oprnt,odump
      common/optef/ jnksp,mode,es,maxstp,ihess,is,negreq,
     1 stpmax,dd,ggrms,gmax
     2 ,drms,dmax,iupdat,maxhes
      common/lambda/lamda,lamda0,it
      logical o_var_high,odumdum
      common /csubst/ values(maxvar),intvec(maxvar),fpvec(maxvar),
     1                cmin10(maxvar),cmin20(maxvar),o_var_high,odumdum
c
      data dzero/0.0d0/, half/0.5d0/, toll/1.0d-8/
      data step/0.05d0/, big/1.0d+3/, maxit/999/
c ...
c
c  ts  search   forms a step by p-rfo that tries to maximize
c               along the direction of a choosen hessian mode
c               and minimize along all other modes.
c  min search   forms a step by simple rfo that attempts to
c               minimize along all hessian modes.
c
c ...
      numit = 0
      it = 0
      if (negreq.ne.0) then
c ...
c  (a) maximization along one of the hessian modes
c ...
         if (mode.ne.0) then
            call ovlap(u,newmod,nvar,vmode,iw)
c ...
c  on return from overlap, newmod is the mode along which
c  the energy is to be maximized
c ...
            if (newmod.ne.mode) write (iw,6010) mode , newmod
            mode = newmod
            if (oprnt) write (iw,6020) mode
            it = mode
c ...
c  if the mode being followed is now the lowest mode,
c  then switch off mode following
c ...
            if (mode.eq.1) then
               mode = 0
               write (iw,6030)
            end if
         else
            if (oprnt) write (iw,6040)
            it = 1
         end if
         lamda0 = eigval(it) + dsqrt(eigval(it)*eigval(it)+4.0d0*fx(it)
     +            *fx(it))
         lamda0 = half*lamda0
         if (oprnt) write (iw,6050) lamda0
         if (nvar.eq.1) go to 40
      end if
c
c  (b) minimization along all other modes
c
      jt = 1 + it
      if (jt.gt.2) jt = 1
c
      if (oprnt .and. negreq.eq.1) write (iw,6060)
      if (oprnt .and. negreq.eq.0) write (iw,6070)
c
c  solve iteratively for lamda
c  initial guess for lamda is zero except note that
c  lamda should be less than eigval(jt)
c
      lamda = dzero
      if (eigval(jt).lt.dzero) then
         lamda = eigval(jt) - step
         lamda1 = eigval(jt)
         lamda2 = -big
      end if
c
 20   numit = numit + 1
      temp = dzero
      do 30 i = 1 , nvar
         if (i.ne.it) then
            temp = temp + (fx(i)*fx(i))/(lamda-eigval(i))
         end if
 30   continue
      if (numit.gt.(maxit-20)) write (iw,6080) numit , lamda , temp
      if (odump) write (iw,6080) numit , lamda , temp
c
c  check for convergence of lamda
c
      if (dabs(lamda-temp).lt.toll) then
c ...
c         at this point we have an acceptable value for lamda
c         make final check on acceptability
c ...
         if (oprnt) write (iw,6050) lamda
         if (lamda.gt.eigval(jt)) then
            call caserr(' error in determining lamda in formd ')
            go to 70
         else if (eigval(jt).gt.dzero .and. lamda.gt.dzero) then
c ....
c    error section
c ....
            call caserr(' error in determining lamda in formd ')
            go to 70
         end if
      else
c
c  check for maximum iterations exceeded
c
         if (numit.gt.maxit) go to 70
c
c  (a) simple iterative scheme
c      (used if eigval(jt) > zero)
c
         if (eigval(jt).gt.dzero) then
            lamda = temp
c
         else
c
c  (b) cautious bracketing scheme
c      (used if eigval(jt) < zero)
c
            if (temp.lt.lamda) lamda1 = lamda
            if (temp.gt.lamda) lamda2 = lamda
            if (lamda2.gt.-big) lamda = half*(lamda1+lamda2)
            if (lamda2.eq.-big) lamda = lamda - step
         end if
         go to 20
      end if
c ...
c  calculate the step
c ...
 40   call vclr(d,1,nvar)
      do 60 i = 1 , nvar
         temp = fx(i)/(lamda-eigval(i))
         if (i.eq.it) temp = fx(i)/(lamda0-eigval(i))
         do 50 j = 1 , nvar
            d(j) = d(j) + temp*u(j,i)
 50      continue
 60   continue
c
      return
 70   call caserr('unable to determine lamda in formd')
 6010 format (' warning!! mode switching. was following mode ',i3,
     +        ' now following mode ',i3,/)
 6020 format (' searching for lamda that maximizes along mode ',i3,/)
 6030 format (' ============================ ',
     +        /' mode following switched off. ',
     +        /' ============================ ')
 6040 format (' searching for lamda that maximizes along the',
     +        ' lowest mode',/)
 6050 format (' value taken    lamda= ',f12.8,/)
 6060 format (' searching for lamda that minimizes along all',
     +        ' other modes',/)
 6070 format (' searching for lamda that minimizes along all modes'/)
 6080 format (' in iterative cycle: ',i3,'  lamda= ',f12.8,' temp= ',
     +        f12.8/)
      end
**==formnr.f
      subroutine formnr(u,eigval,fx,d,nvar)
c
      implicit real*8  (a-h,p-w) , integer  (i-n) , logical (o)
      implicit character*8 (z) , character*4 (y) , character*1 (x)
      dimension u(nvar,nvar),eigval(*),fx(*),d(*)
c ...
c     take the newton-raphson step
c     (i.e. set lamda to zero)
c ...
      call vclr(d,1,nvar)
      do 30 i = 1 , nvar
         temp = -fx(i)/eigval(i)
         do 20 j = 1 , nvar
            d(j) = d(j) + temp*u(j,i)
 20      continue
 30   continue
      return
      end
**==magchk.f
      subroutine magchk(eigval,nvar,iw)
c
      implicit real*8  (a-h,p-w) , integer  (i-n) , logical (o)
      implicit character*8 (z) , character*4 (y) , character*1 (x)
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
      dimension eigval(*)
      common/jorgs/ maxj,irecal,opowel,obfgs,obfgx,onrfo,ocut,
     1 oprnt,odump
      common/optef/ jnksp,mode,es,maxstp,ihess,is,negreq,
     1 stpmax,dd,ggrms,gmax
     2 ,drms,dmax,iupdat,maxhes
c
      real*8 acc, stol, stepmx, tolmax, tolstp, tanstp
      real*8 accin, eigmin, eigmax, grms
      integer jm, mnls, nls, isadle, ifcm, iblfcm, isecfcm
      integer nentry, iterch, nftotl, ismax, lintyp, lqstyp
      logical ofcm

      common /cntl1/ acc,stol,stepmx,tolmax,tolstp,tanstp,
     +               jm,mnls,nls,isadle,imnter,idumcntl1,
     +               accin(6),nftotl,ismax,lintyp,lqstyp,
     +               eigmin,eigmax,grms(2),nentry,iterch,
     +               ifcm,iblfcm,isecfcm,ofcm
c
      logical o_var_high,odumdum
      common /csubst/ values(maxvar),intvec(maxvar),fpvec(maxvar),
     1                cmin10(maxvar),cmin20(maxvar),o_var_high,odumdum
c
      data dzero/0.0d0/
c ...
c      first check if the magnitudes of the lowest eigenvalues
c      exceed eigmin
c ...
      do 20 i = 1 , nvar
         if (eigval(i).gt.eigmin) go to 30
         if (dabs(eigval(i)).lt.eigmin) then
            if (eigval(i).lt.dzero) eigval(i) = -eigmin
            if (eigval(i).ge.dzero) eigval(i) = eigmin
            if (oprnt) write (iw,6010) i , eigval(i)
         end if
 20   continue
c ...
c                 now check if the magnitude of the highest eigenvalues
c                 exceed eigmax
c ...
 30   do 40 i = 1 , nvar
         ii = nvar - i + 1
         if (eigval(ii).lt.eigmax) go to 50
         if (dabs(eigval(ii)).gt.eigmax) then
            if (eigval(ii).lt.dzero) eigval(ii) = -eigmax
            if (eigval(ii).ge.dzero) eigval(ii) = eigmax
            if (oprnt) write (iw,6020) ii , eigval(ii)
         end if
 40   continue
 50   return
 6010 format (' warning! eigenvalue ',i2,' too small. replaced',' by ',
     +        f12.6)
 6020 format (' warning! eigenvalue ',i2,' too large. replaced',' by ',
     +        f12.6)
      end
**==matprn.f
      subroutine matprn(vec,nvar,eig,zvar,ofd,iw)
c
c ...
c     print out matrix of eigenvectors along with eigenvalues.
c ...
      implicit real*8  (a-h,p-w) , integer  (i-n) , logical (o)
      implicit character*8 (z) , character*4 (y) , character*1 (x)
      dimension vec(nvar,nvar),eig(nvar),zvar(nvar)
      ohead = .true.
      is = 1
      it = nvar
      if (it.gt.7) it = 7
 20   do 30 i = 1 , nvar
         if (ohead) then
            write (iw,6010) (j,j=is,it)
            write (iw,6020) (eig(j),j=is,it)
            write (iw,6030)
            ohead = .false.
         end if
         if (ofd) then
            write (iw,6040) i , zvar(i) , (vec(i,j),j=is,it)
         else
            write (iw,6050) i , (vec(i,j),j=is,it)
         end if
 30   continue
      ohead = .true.
      if (it.lt.nvar) then
         is = is + 7
         it = is + 6
         if (it.gt.nvar) it = nvar
         go to 20
      end if
      return
 6010 format (21x,7(i3,13x))
 6020 format (1x,'eigenvalues:-',3x,7(f11.7,5x))
 6030 format (/)
 6040 format (3x,i3,2x,a8,1x,7(f11.7,5x))
 6050 format (3x,i3,11x,7(f11.7,5x))
      end
**==optjs.f
      subroutine optjs(hes   ,hesi10,hesi20,hesi40,hesi50,hesiix,
     &                 hesiig,hessp1,hessp2,hessp3,hessp4,hessp5,hessp6,
     &                 heson1,heson2,nprino)
c  **************************************************************
c  *                                                            *
c  *  driver for eigenvector following transition state search  *
c  *  based on "on finding transition states" by cerjan and     *
c  *  miller (j.chem.phys. 75 (1981) 2800);   "walking on       *
c  *  potential energy surfaces" by simons, jorgenson, taylor   *
c  *  and ozmet (j.phys.chem. 87 (1983) 2745)  and  "search     *
c  *  for stationary points on surfaces" by banerjee, adams,    *
c  *  simons and shepard (j.phys.chem. 89 (1985) 52).           *
c  *                                                            *
c  **************************************************************
c  mode following:  mode following turned on via the intvec array, which
c  is input with the z-matrix variables :adding 4 to the intvec
c  entry for a particular variable will cause a transition state search
c  to follow the hessian mode with the largest magnitude component for
c  that variable: adding 10 to intvec entry for the kth variable will
c  follow the kth hessian mode. only one intvec entry should be modified
c  in this way  i.e. only one mode should be followed at a time.
c ...
      implicit real*8  (a-h,p-w) , integer  (i-n) , logical (o)
      implicit character*8 (z) , character*4 (y) , character*1 (x)
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
      parameter ( mvar10 = maxvar * 10 )
      dimension hes(*),hesi10(*),hesi20(*),hesi40(*),hesi50(*)
      dimension hesiix(*),hesiig(*),hessp1(*),hessp2(*),hessp3(*)
      dimension hessp4(*),hessp5(*),hessp6(*),heson1(*),heson2(*)
      dimension accu(6)
      common/jorgs/ maxj,irecal,opowel,obfgs,obfgx,onrfo,ocut,
     1 oprnt,odump
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer ibl3s, ibl3t, ibl3f, ibl3x, ibl3y, ibl3z
      integer ibl3rs, ibl3qa, ibl3pa, ibl3ea, ibl3qb, ibl3pb
      integer ibl3eb, ibl3g, ibl3hs, ionsec, ibl3op, ions2
      integer isecqm
      common /dump3/ ibl3s, ibl3t, ibl3f, ibl3x, ibl3y, ibl3z,
     +               ibl3rs,ibl3qa,ibl3pa,ibl3ea,ibl3qb,ibl3pb,
     +               ibl3eb,ibl3g,ibl3hs,ionsec,ibl3op,ions2,
     +               isecqm
c
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
      common/optef/ jnksp,mode,es,maxstp,ihess,is,negreq,
     1 stpmax,dd,ggrms,gmax
     2 ,drms,dmax,iupdat,maxhes
      logical o_var_high,odumdum
      common /csubst/ values(maxvar),intvec(maxvar),fpvec(maxvar),
     1                cmin10(maxvar),cmin20(maxvar),o_var_high,odumdum
      character *8 zvar
      common/csubch/zvar(maxvar)
c
      integer ianz, iz, lbl, lalpha, lbeta, nz, nvar
      real*8 bl, alpha, beta
      common /czmat/ ianz(maxnz),iz(maxnz,4),bl(maxnz),
     +               alpha(maxnz),beta(maxnz),lbl(maxnz),
     +               lalpha(maxnz),lbeta(maxnz),nz,nvar
c
c
      real*8 acc, stol, stepmx, tolmax, tolstp, tanstp
      real*8 accin, eigmin, eigmax, grms
      integer jm, mnls, nls, isadle, ifcm, iblfcm, isecfcm
      integer nentry, iterch, nftotl, ismax, lintyp, lqstyp
      logical ofcm

      common /cntl1/ acc,stol,stepmx,tolmax,tolstp,tanstp,
     +               jm,mnls,nls,isadle,imnter,idumcntl1,
     +               accin(6),nftotl,ismax,lintyp,lqstyp,
     +               eigmin,eigmax,grms(2),nentry,iterch,
     +               ifcm,iblfcm,isecfcm,ofcm
c
c
      real*8 en, accd
      integer ig, ierec, igrec, iterat, jpoint
      common /cntl2/ en(200),ig(200),ierec(200),igrec(200),
     +               accd(6,200),iterat,jpoint
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
      common /miscop/ cmin1(maxvar),cmin2(maxvar),
     1 coefb(maxvar),coefc(maxvar),coefd(maxvar),coefe(maxvar),
     2 coefsp(mvar10),spppp(31),fspa(3*maxat),fspaa(maxvar),nt,ia(maxat)
c
      real*8 enrgy, egrad
      common/gms_funct/enrgy,egrad(3*maxat)
c
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
      equivalence ( accu(1) , acc )
c
      data zprop /'prop'/
      data dzero/0.0d0/
c
      write (iwr,6230)
      ihess = -1
      ioutp = iwr
      negreq = 1
      nvar2 = nvar*nvar
      nvar60 = nvar*200
      nvarb = lensec(nvar60)
      lenh = lensec(mach(1)) + lensec(mach(10)) + nvarb + nvarb
c ...
c     check the intvec array.
c     the ic array is read in with the variables and the z-matrix
c     it can control the setting up of the initial hessian (see below)
c     and also determine if hessian mode following is switched on.
c
c     hessian mode following. follow mode mode.
c ...
      mode = 0
      do 20 i = 1 , nvar
         ici = intvec(i)
         if (ici.gt.9) then
            intvec(i) = ici - 10
            mode = i + nvar
         else if (ici.gt.3) then
            intvec(i) = ici - 4
            mode = i
         else if (ici.lt.0 .or. ici.gt.3) then
            write (iwr,6050) i
            call caserr(' optjs : bad ic entry for variable. ')
         end if
 20   continue
c ...
c      set optimization control parameters.
c ...
      if (isadle.ne.5) negreq = 0
      maxstp = maxj
      if (opowel) then
         iupdat = 0
      else if (obfgs) then
         iupdat = 1
      else
         iupdat = 2
      end if
      stpmax = stepmx
      if (ocut) then
         stpmax = 1.0d0
         eigmin = dzero
         eigmax = 1.0d12
      end if
      if (irecal.eq.0) irecal = 99999
c ...
c                 calculate remaining convergence factors
c ...
      gmax = stol*0.25d0
c ...
c                 maximum gradient.
c ...
      ggrms = stol*0.1666666666d0
c ....
c                 rms force
c ....
      dmax = stol
c ...
c                 maximum displacement.
c ...
      drms = stol*0.666666666d0
c ...
c                 rms displacement convergence.
c ...
      write (iwr,6060) gmax , drms , dmax , maxstp , stpmax , irecal
c ...
c     check hessian mode following if this is switched on.
c ...
      if (mode.gt.0) then
         if (mode.gt.nvar) then
            ihop = mode - nvar
            write (iwr,6080) ihop
         else
            write (iwr,6090) mode
         end if
         if (mode.gt.2*nvar) call caserr('mode following out of range')
      end if
      if (opowel) write (iwr,6030)
      if (obfgs) write (iwr,6010)
      if (obfgx) write (iwr,6020)
      nt = nat
c ...
c     jpoint  -- counts number of points in optimisation
c                ( calls to newpt )
c                ( including points in the evaluation of the hessian ).
c     iteraj  -- counts points in the jorgensen-simons optimisation
c                ( excludes points in the evaluation of the hessian ).
c     itera   -- ?
c ...
      iteraj = 1
      ifail = 0
      imnter = 0
c ...
c     hessian is still in hesi10, diagonalized hessian is in hes(i20)
c     print the hessian, its eigenvectors and eigenvalues if
c     requested. vectors are in hesi20, eigenvalues in hessp6)
c ...
 30   if (odump) then
         write (iwr,6100)
         ofd = .false.
         call matprn(hesi10,nvar,hessp6,zvar,ofd,iwr)
      end if
      if (oprnt) then
         write (iwr,6110)
         ofd = .true.
         call matprn(hesi20,nvar,hessp6,zvar,ofd,iwr)
      else
         write (iwr,6120)
         write (iwr,6130) (hessp6(l),l=1,nvar)
      end if
c ...
c      check for maximum and minimum allowed magnitudes
c      for the hessian eigenvalues
c ...
      call magchk(hessp6,nvar,iwr)
c ...
c     calculate neg, the number of negative eigenvalues
c ...
      call fndneg(hessp6,neg,nvar)
c ...
c     form the fx vector
c     ( the component of f along the local eigendirections )
c ...
      do 40 i = 1 , nvar
         hessp5(i) = ddot(nvar,hesi20(nvar*(i-1)+1),1,hesiig,1)
 40   continue
c ...
c          if we are in the "right" region of the energy surface
c          i.e. the hessian has the required number of negative eigenval
c          and onrflg ( onrfo ) is on, take the newton-raphson step
c ...
      if (onrfo .and. neg.eq.negreq) then
         write (iwr,6140)
         call formnr(hesi20,hessp6,hessp5,hessp2,nvar)
      else
c ...
c              take the p-rfo step for a ts search and
c              the simple rfo step for a minimum search
c ...
         if (neg.eq.negreq .and. negreq.eq.1) write (iwr,6150)
         if (neg.eq.negreq .and. negreq.eq.0) write (iwr,6160)
         if (neg.ne.negreq .and. negreq.eq.1) write (iwr,6170)
         if (neg.ne.negreq .and. negreq.eq.0) write (iwr,6180)
         call formd(hesi20,hessp6,hessp5,nvar,hessp2,hessp3,iwr)
      end if
c ...
c                        we now have a new step in hessp2
c                        check that the stepsize does not exceed stpmax
c                        if so, scale.
c ...
      call checkd(hessp2,stpmax,dd,skal,nvar,iwr)
c ...
c                        predict the change in energy
c ...
      call estime(hessp6,hessp5,dd,skal,chnge,nvar)
c ...
c            test for convergence and print out current parameter values
c ...
      call convef(neg,oconv,hessp4,hesiig,hesiix,hessp2,iwr)
      write (iwr,6190) chnge
c ...
c                         if converged then exit
c ...
      if (oconv) write (iwr,6200)
 50   if (oconv) then
c ...
c         finish up optimisation run by call to calcfg wtih iflag = 8 ,
c ...
         call calcfg(8,nvar,hesiix,enrgy,hesiig,hesi40,hesi50,ioutp,hes,
     +               heson1,heson2)
         nprint = nprino
         if (ifail.ne.0) then
            if (.not.oprint(44)) go to 60
         end if
         write (iwr,6040)
         call intr(hes)
         call newpt(nvar,hesiix,enrgy,0,ioutp,hes)
         call optend(hes,nprint)
 60      zruntp = zprop
         write (iwr,6240)
         return
      end if
c ...
c                form the new hesiix ( internal coordinates) vector.
c ...
      call dvplv(hesiix,hesiix,hessp2,nvar)
c ...
c                save the current gradients.
c ...
      call dcopy(nvar,hesiig,1,hessp1,1)
      iteraj = iteraj + 1
      if (iteraj.gt.maxstp) then
         write (iwr,6210)
         oconv = .true.
         ifail = 0
         go to 50
      end if
      ihess = ihess + 1
c ...
c                         s.c.f. calculation.
c ...
      call calcfg(0,nvar,hesiix,enrgy,hesiig,hesi40,hesi50,ioutp,hes,
     +            heson1,heson2)
c ...
c                         gradients calculation.
c ...
      call calcfg(3,nvar,hesiix,enrgy,hesiig,hesi40,hesi50,ioutp,hes,
     +            heson1,heson2)
c ...
c ...
c                        check what is to be done with the hessian.
c ...
      if (((ihess-irecal).eq.0) .or. (iteraj.eq.1) .and. (ifcm.le.0))
     +    then
c ...
c         on 2nd. or greater iteration of optimisation and the hessian
c         is to be recalculated.                               *******
c                  ************
c                                  or
c         on 1st. iteration and hessian is to be calculated numerically.
c ...
         write (iwr,6220)
         call star(hesi10,hesi40,hesi50,hesiix,enrgy,hesiig,nvar,ioutp,
     +             hes,heson1,heson2)
         ihess = 0
c ...
c             numerically estimated hessian is in hesi10
c             update the hessian (unless first point).
c ...
      else
c ...
c
c         on 1st. iteration and the hessian has been taken off a foreign
c         dumpfile.
c                                 or
c         on 2nd. or greater iteration of optimisation and the hessian
c         need not be recalculted this iteration. therefore  only update
c         hessian.
c                                 or
c         on 2nd. or greater iteration of optimisation and the hessian
c         is not to be recalculted during this optimisation run, merely
c         updated each iteration.
c ...
         if (iteraj.gt.1) call updhes(hesi10,hessp6,hessp4,hesiig,
     +                                hessp1,hessp2,iwr)
      end if
c ...
c          write out hessian to dumpfile.
c ...
      call wrt3(hesi10,nvar2,lenh+ibl3op,idaf)
      write (iwr,6070) iteraj
c ...
c                         diagonalize hessian.
c ...
      ireord = 0
      ineg = -1
      call diaeig(hesi10,hesi20,nvar,hessp6,eigmin,eigmax,ineg,ireord,
     +            iwr)
      go to 30
c
 6010 format (/' update of hessian using bfgs method')
 6020 format (/' update of hessian using modified bfgs method')
 6030 format (/' update of hessian using powell method')
 6040 format (/,40x,22('*')/40x,'geometry of last point'/40x,22('*')/)
 6050 format (1x,'call from initef:  bad ic entry for variable ',i3)
 6060 format (/11x,53('*')/10x,
     +        ' **  o p t i m i z a t i o n   p a r a m e t e r s  **'/
     +        10x,' **  maximum gradient                 :  ',f10.6,
     +        ' **'/10x,' **  rms displacement convergence     :  ',
     +        f10.6,' **'/10x,
     +        ' **  maximum displacement             :  ',f10.6,
     +        ' **'/10x,' **  maximum no. of iterations        :  ',i10,
     +        ' **'/10x,' **  maximum allowed step length      :  ',
     +        f10.6,' **'/10x,
     +        ' **  steps before hessian reevaluated :  ',i10,' **'/11x,
     +        53('*')//)
 6070 format (/1x,
     +        '********************************************************'
     +        ,/1x,'* jorgenson & simons optimisation iteration number',
     +        i4,' *',/1x,
     +        '********************************************************'
     +        ,//)
 6080 format (' hessian mode following switched on',/' following mode ',
     +        i3,/)
 6090 format (' hessian mode following switched on',
     +        /' following mode with largest magnitude component',
     +        ' for variable ',i3,/)
 6100 format (/' force constant matrix and eigenvalues',
     +        /' =====================================',/)
 6110 format (/' hessian eigenvectors and eigenvalues',/1x,
     +        '====================================',/)
 6120 format (' eigenvalues of the hessian',1x,
     +        '=========================='/)
 6130 format (6x,6f12.6,/)
 6140 format (' ',/,' hessian has required local structure'//
     +        ' taking newton-raphson step',/)
 6150 format (1x,'transition state search: taking p-rfo step',/)
 6160 format (1x,'minimum search. taking simple rfo step',/)
 6170 format (1x,'hessian does not have the desired local structure'/
     +        ' taking p-rfo step',/)
 6180 format (1x,'hessian does not have the desired local structure'/
     +        ' taking simple rfo step',/)
 6190 format (1x,'predicted change in energy  ',e20.10/)
 6200 format (' ***************************',
     +        /' * optimisation converged. *',
     +        /' ***************************'/)
 6210 format (' ',' ************* too many steps in optimisation. ',/)
 6220 format (1x,
     +        '*******************************************************',
     +        '**********************',/
     +       ' * numerically reevaluating the hessian before continuing'
     +       ,' with optimisation. *',
     +       /' ******************************************************',
     +       '***********************')
 6230 format (' ',/,' ',
     +        '********************************************************'
     +        ,/' ',
     +        '* jorgenson & simons optimisation code now in control  *'
     +        ,/' ',
     +        '********************************************************'
     +        ,//)
 6240 format (' ',/,' ',
     +        '********************************************************'
     +        ,/' ',
     +        '* jorgenson & simons optimisation control now ended    *'
     +        ,/' ',
     +        '********************************************************'
     +        ,//)
      end
**==ovlap.f
      subroutine ovlap(u,newmod,nvar,vmode,iw)
c
      implicit real*8  (a-h,p-w) , integer  (i-n) , logical (o)
      implicit character*8 (z) , character*4 (y) , character*1 (x)
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
      dimension u(nvar,nvar),vmode(*)
      common/jorgs/ maxj,irecal,opowel,obfgs,obfgx,onrfo,ocut,
     1 oprnt,odump
      common/optef/ jnksp,mode,es,maxstp,ihess,is,negreq,
     1 stpmax,dd,ggrms,gmax
     2 ,drms,dmax,iupdat,maxhes
c
      real*8 en, accd
      integer ig, ierec, igrec, iterat, jpoint
      common /cntl2/ en(200),ig(200),ierec(200),igrec(200),
     +               accd(6,200),iterat,jpoint
c
      logical o_var_high,odumdum
      common /csubst/ values(maxvar),intvec(maxvar),fpvec(maxvar),
     1                cmin10(maxvar),cmin20(maxvar),o_var_high,odumdum
c ...
c           on the first step simply determine which mode to follow
c ...
      if (iterat.eq.1) then
c ...
c  (a) following a given mode
c ...
         if (mode.gt.nvar) then
            it = mode - nvar
            go to 30
         end if
c ...
c  (b) find the hessian eigenvector with the largest magnitude
c  component in the position given by mode (read in via ic array)
c ...
         it = 1
         tovlp = dabs(u(mode,1))
         do 20 i = 2 , nvar
            if (dabs(u(mode,i)).gt.tovlp) then
               tovlp = dabs(u(mode,i))
               it = i
            end if
 20      continue
 30      mode = it
         if (oprnt) write (iw,6020) mode
      else
c
c  on subsequent steps determine which hessian eigenvector has
c  the greatest overlap with the mode we are following
c
         it = 1
         tovlp = ddot(nvar,u(1,1),1,vmode,1)
         tovlp = dabs(tovlp)
         do 40 i = 2 , nvar
            rovlp = ddot(nvar,u(1,i),1,vmode,1)
            rovlp = dabs(rovlp)
            if (rovlp.gt.tovlp) then
               tovlp = rovlp
               it = i
            end if
 40      continue
         if (oprnt) write (iw,6010) tovlp
      end if
c
c  save the eigenvector in vmode
c
      do 50 i = 1 , nvar
         vmode(i) = u(i,it)
 50   continue
      newmod = it
      return
 6010 format (' overlap of current mode with previous mode is ',f12.6)
 6020 format (' hessian mode following switched on'/' following mode ',
     +        i3)
      end
**==updhes.f
      subroutine updhes(hess,svec,tvec,g,roldf,d,iw)
c
      implicit real*8  (a-h,p-w) , integer  (i-n) , logical (o)
      implicit character*8 (z) , character*4 (y) , character*1 (x)
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
      common/jorgs/ maxj,irecal,opowel,obfgs,obfgx,onrfo,ocut,
     1 oprnt,odump
      common/optef/ jnksp,mode,es,maxstp,ihess,is,negreq,
     1 stpmax,dd,ggrms,gmax
     2 ,drms,dmax,iupdat,maxhes
c
      integer ianz, iz, lbl, lalpha, lbeta, nz, nvar
      real*8 bl, alpha, beta
      common /czmat/ ianz(maxnz),iz(maxnz,4),bl(maxnz),
     +               alpha(maxnz),beta(maxnz),lbl(maxnz),
     +               lalpha(maxnz),lbeta(maxnz),nz,nvar
c
      dimension tvec(*),svec(*),hess(nvar,nvar),g(*),d(*),roldf(*)
      data dzero/0.0d0/
c ...
c  updating of the hessian
c  depends on current gradients, old gradients and the
c  correction vector used on the last cycle
c  svec & tvecare for temporary storage
c
c  2 updating procedures are possible
c  (i)   the powell update
c        this preserves the symmetric character of the hessian
c        whilst allowing its eigenvalue structure to change.
c        it is the default update for a transition state search
c  (ii)  the bfgs update
c        this update has the important characteristic of retaining
c        positive definiteness (note: this is not rigorously
c        guaranteed, but can be checked for by the program).
c        it is the default update for a minimum search
c ...
      call dmatdv(tvec,hess,d,nvar)
      if (iupdat.eq.0) then
c ...
c   (i) powell update
c ...
         if (oprnt) write (iw,6010)
         do 20 i = 1 , nvar
            tvec(i) = g(i) - roldf(i) - tvec(i)
 20      continue
         dds = dd*dd
         ddtd = ddot(nvar,tvec,1,d,1)
         ddtd = ddtd/dds
         do 40 i = 1 , nvar
            do 30 j = 1 , i
               temp = tvec(i)*d(j) + d(i)*tvec(j) - d(i)*ddtd*d(j)
               hess(i,j) = hess(i,j) + temp/dds
               hess(j,i) = hess(i,j)
 30         continue
 40      continue
      else
c ...
c  (ii) bfgs update
c ...
         if (oprnt) write (iw,6040)
         do 50 i = 1 , nvar
            svec(i) = g(i) - roldf(i)
 50      continue
         dds = ddot(nvar,svec,1,d,1)
c ...
c  if dds is negative, retention of positive definiteness is not
c  guaranteed. print a warning and skip the update this time if
c  requested
c ...
         if (dds.lt.dzero) then
            write (iw,6020)
            if (iupdat.eq.2) then
               write (iw,6030)
               return
            end if
         end if
         ddtd = ddot(nvar,d,1,tvec,1)
         do 70 i = 1 , nvar
            do 60 j = 1 , i
               temp = (svec(i)*svec(j))/dds - (tvec(i)*tvec(j))/ddtd
               hess(i,j) = hess(i,j) + temp
               hess(j,i) = hess(i,j)
 60         continue
 70      continue
      end if
      return
 6010 format (/' hessian updated using powell update')
 6020 format (/' warning! hereditary positive definiteness endangered')
 6030 format (/' update skipped this cycle')
 6040 format (/' hessian updated using bfgs update')
      end
c ******************************************************
c ******************************************************
c             =   irc    =
c ******************************************************
c ******************************************************
      subroutine gldia(ldvect,ndim,h,wrk,eig,vector,ierr,iwrk,iw)
c
c     ----- general routine to diagonalize a matrix -----
c     if kdiag = 0, use a routine from the vector library,
c                   if available (see the subroutine 'gldiag'
c                   in vector.src), or evvrsp otherwise
c              = 1, use evvrsp
c              = 2, use giveis
c              = 3, use jacobi
c
c           ldvect = row dimension of vector
c           ndim   = dimension (order) of matrix to be solved
c           h      = matrix to be diagonalized
c           wrk    = n*8 w.p. real words of scratch space
c           eig    = eigenvectors (output)
c           vector = eigenvectors (output)
c           ierr   = error flag (output)
c           iwrk   = n integer words of scratch space
c
      implicit real*8 (a-h,o-z)
c
      dimension vector(ldvect,*), wrk(ndim,8), eig(ndim), iwrk(ndim)
     *         ,h(*)
c
c
      ierr = 0
      call jacdia(iw,h,vector,eig,iwrk,wrk,ldvect)
      return
      end
      subroutine jacdia(iw, a, vec, eig, ia, big, n)
c
c     ----- routine for diagonalization -----
c              of symmetric matrix a  in triangular form
c
      implicit real*8 (a-h,o-z)
      dimension a(*), vec(n,*), eig(*), ia(*), big(*)
c
c...   common for overruling local criteria (e.g. dependency, diag thresholds)
      integer  i_depen_check, i_global_diag, n_depen_check
      logical  o_depen_print
      common/overrule/ i_depen_check,i_global_diag, o_depen_print, 
     1                 n_depen_check
c
c     i_depen_check : criterion for throwing out vectors
c     n_depen_check : number  of vectors to be thrown out 
c     o_depen_print : if true +> extra print 

      parameter (zero=0.0d+00, one=1.0d+00)
c
      conv = 1.0d-14
      if (i_global_diag.ne.-999) conv = 10.0d0**(-i_global_diag*1.0d0)
c  
      maxit = 2048
      do 20 i = 1, n
        do 10 j = 1, n
          vec(j,i) = zero
   10   continue
        vec(i,i) = one
   20 continue
      call jcdg(iw, vec, eig, a, n, conv, maxit, convg, nit, big, ia)
      if (nit .ge. maxit) write (iw,30) nit
      return
   30 format (/' jacobi diagonalization did not converge in ', i4,
     1       ' iterations')
      end
      subroutine jcdg(iw, vc, eps, h, n, conv8, mxniw, convg, ni, big,
     1           jbig)
      implicit real*8 (a-h,o-z)
      dimension vc(n,n), eps(*), h(*)
      dimension big(*), jbig(*)
c
      parameter (dmin=1.d-4)
      parameter (root2=0.707106781186547531d0)
      parameter (thr1=8796093022208.00000d0)
      parameter (thr2=0.390625000000000000d-02)
      parameter (thr3=0.119209289550781250d-06)
      parameter (thr4=0.244140625000000000d-03)
      parameter
     * (zero=0.0d+00, one=1.0d+00, half=0.5d+00, thalvs=1.5d+00,
     *  elevei=1.375d+00, quart=0.25d+00, thonei=3.875d+00,
     *  ten=10.0d+00, vsmall=1.0d-37)
c
c
      msweep = 0
   10 ni = 1
      msweep = msweep + 1
      if (msweep .gt. 5) then
        write (iw,300) msweep
 300    format(' ***** jcdg: msweep = ',i6)
        call caserr('msweep exceeded in jcdg')
      end if
      j = 0
      do 20 i = 1, n
        j = j + i
   20 eps(i) = h(j)
      if (n .eq. 1) return
      iear = 1
      big(1) = zero
      jbig(1) = 0
      do 40 i = 2, n
        big(i) = zero
        jbig(i) = 0
        j = i - 1
        do 30 k = 1, j
          if (dabs(big(i)) .ge. dabs(h(iear + k))) go to 30
          big(i) = h(iear + k)
          jbig(i) = k
   30   continue
   40 iear = iear + i
      do 230 ni = 1, mxniw
        sd = dabs(eps(1))
        do 50 j = 2, n
          if (sd .gt. dabs(eps(j))) sd = dabs(eps(j))
   50   continue
        sd = max(sd,dmin)
        t = zero
        do 60 i = 1, n
          if (t .ge. dabs(big(i))) go to 60
          t = dabs(big(i))
          ib = i
          ia = jbig(i)
   60   continue
        if (t .le. conv8*sd) go to 240
        sx = root2
        cx = root2
        if (dabs(eps(ia) - eps(ib)) .le. vsmall) go to 100
        t2x2 = big(ib) / (eps(ia) - eps(ib))
        abst = dabs(t2x2)
        if (abst .ge. thr1) go to 100
        if (abst .gt. thr2) go to 90
        if (abst .gt. thr3) go to 70
        cx = one
        sx = t2x2
        go to 100
   70   t2x25 = abst * abst
        if (abst .gt. thr4) go to 80
        cx = one-half * t2x25
        sx = t2x2 * (one-thalvs*t2x25)
        go to 100
   80   cx = one+t2x25 * (t2x25*elevei-half)
        sx = t2x2 * (one+t2x25*(t2x25*thonei-thalvs))
        go to 100
   90   t = quart / dsqrt(quart+abst*abst)
        cx = dsqrt(half+t)
        sx = dsqrt(half-t)
        sx = dsign(sx,t2x2)
  100   iear = (ia*(ia - 1)) / 2 + 1
        iebr = (ib*(ib - 1)) / 2 + 1
        do 210 ir = 1, n
          t = h(iear) * sx
          h(iear) = h(iear) * cx + h(iebr) * sx
          h(iebr) = t - h(iebr) * cx
          if (ir - ia) 200, 110, 120
  110     tt = h(iebr)
          ieaa = iear
          ieab = iebr
          h(iebr) = big(ib)
          iear = iear + ir - 1
          if (jbig(ir) .eq. 0) go to 200
          go to 180
  120     t = h(iear)
          it = ia
          t = zero
          it = 0
          iear = iear + ir - 1
          if (ir - ib) 160, 130, 140
  130     h(ieaa) = h(ieaa) * cx + h(ieab) * sx
          eps(ia) = h(ieaa)
          h(ieab) = tt * cx + h(iebr) * sx
          h(iebr) = tt * sx - h(iebr) * cx
          eps(ib) = h(iebr)
          iebr = iebr + ir - 1
          go to 180
  140     if (dabs(t) .ge. dabs(h(iebr))) go to 150
          t = h(iebr)
          it = ib
  150     iebr = iebr + ir - 1
  160     if (dabs(t) .lt. dabs(big(ir))) go to 170
          big(ir) = t
          jbig(ir) = it
          go to 200
  170     if (ia .ne. jbig(ir) .and. ib .ne. jbig(ir)) go to 200
  180     k = iear - ir - ia + 2
          big(ir) = zero
          ir1 = ir - 1
          do 190 i = 1, ir1
            if (dabs(big(ir)) .ge. dabs(h(k))) go to 190
            big(ir) = h(k)
            jbig(ir) = i
  190     k = k + 1
  200     iear = iear + 1
  210   iebr = iebr + 1
        do 220 i = 1, n
          t = vc(i,ib) * sx
          vc(i,ib) = vc(i,ia) * sx - vc(i,ib) * cx
  220   vc(i,ia) = vc(i,ia) * cx + t
  230 continue
      ni = mxniw + 1
  240 convg = t
      mkl = 1
      do 260 mki = 2, n
        mkl = mkl + mki
        mkj1 = mkl - mki + 1
        mkj2 = mkl - 1
        do 250 mkj = mkj1, mkj2
          if (dabs(h(mkj)) .gt. conv8*ten) then
c-mk        write(iw,fmt='(105x,'' i j h(j)'',2i4,1p,d9.1)')
c-mk  *         mki,mkj,h(mkj)
            go to 10
          end if
  250   continue
  260 continue
      do 290 i = 1, n
        jj = i
        do 270 j = i, n
          if (eps(j) .lt. eps(jj)) jj = j
  270   continue
        if (jj .eq. i) go to 290
        t = eps(jj)
        eps(jj) = eps(i)
        eps(i) = t
        do 280 j = 1, n
          t = vc(j,jj)
          vc(j,jj) = vc(j,i)
  280   vc(j,i) = t
  290 continue
      return
      end
      subroutine iniirc(nc2,core,hess,a,ia,scr,vec,e,iw)
      implicit real*8 (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
c
      parameter (szero=0.0d0,sone=1.0d0,nnam=23,mxa3=maxat*3,
     *           two=2.0d+00, small=1.0d-05,fact=2.642461d+07)
c
c ...
c ... timing statistics
c ...
      real*8 begin, ebegin, timsec, walsec, tstart, estart
      common/statis/begin,ebegin,timsec(50),walsec(50),
     + tstart,estart
c
c
      integer nmaxly, ntotly, icurly, isecly
      logical oprintm
      common /segm/ nmaxly,ntotly,icurly,isecly(100),
     +              oprintm
c
      common/junkc /ylab(26),ztype(831),zoptio(13),ztagg(100)
      common/data1 /vlist(400),np(206),np1(maxorb,2),norbt(6),
     *              norbta(maxorb,2),norb2,nseco(maxorb),norb3,
     *              nthird(maxorb),norbg,norbga(maxorb),nscrp,
     *              range(2),alimit,blimit,isymm,nnx,lsym,nsymg(8),
     *              iuse(8,6),idorb(100),ncont,nmax(10)
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
c from March 97, amass is used only for user-input
c mass values
c czanr contains a copy of czan made in basis.m after reorganisation and
c before pseudo potential; so real charges in real order 
c beware, infob is explicit in direct.m and  dirrpa.m
c
      real*8 czann, czanr, cat, amass, cnew
      integer nonsym, map80
      logical ozmat
      common/infob/czanr(maxat),czann(maxat),cat(3,maxat),amass(maxat),
     +             cnew(maxat3,3),nonsym,map80(maxat),ozmat

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
      integer numdu, iblkdu, maxb, iblkla
      logical orevis
      real*8 apos
      common/sector/numdu,iblkdu,orevis(2),
     + apos(508),maxb,iblkla
c
      common /irc1 /ea,eb,ec,sab,sbc,stotal,ganorm,gbnorm,gcnorm,
     *              qa(mxa3),qb(mxa3),qc(mxa3),ga(mxa3),gb(mxa3),
     *              gc(mxa3),qa2(mxa3),qa1(mxa3),qa0(mxa3),ga2(mxa3),
     *              ga1(mxa3),ga0(mxa3),qcp(mxa3),qcc(mxa3),gcp(mxa3),
     *              cmode(mxa3),ebis(10),cbis(10),gbis,rtrms(mxa3),
     *              ck0(mxa3),ck1(mxa3),ck2(mxa3),ck3(mxa3),
     *              ipirc,npirc
      common /irc2 /oforw,osadd,ostab,otsen,ostbs,olimn,oline,oquad,
     *              oampc,orung,ovtst,ndummy,evib,delta,elbow,stride,
     *              npoint,nextpt,freq,ypac,yyy,ncoord,nhopt,nampc,
     *              ofirst,nn
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
      dimension hess(ncoord,ncoord),a(nc2),e(ncoord),
     *          ia(ncoord),scr(ncoord,8),vec(ncoord,ncoord) 
c    *          cmode(ncoord)
      dimension core(*)
c
      odebug=.false.
c     odebug=.true.
      ihess=ncoord*ncoord
c
c ---- this routine generates an initial point with a nonzero
c      gradient vector by following the imaginary normal mode
c      a short distance away from the saddle point 
c
      write(iw,9000)
c
c ---- store saddle point coords in qa, with szero gradient 
c
      stotal = szero
      write(iw,9002) 
         ii=0
      do 100 i=1,nat
         qa(ii+1) = c(1,i)
         qa(ii+2) = c(2,i)
         qa(ii+3) = c(3,i)
         ga(ii+1) = szero
         ga(ii+2) = szero
         ga(ii+3) = szero
         write(iw,9003) zaname(i),qa(ii+1),qa(ii+2),qa(ii+3)
         ii=ii+3
         ga(i) = szero
  100 continue
      if(nhopt.eq.0) go to 200
c
c ---- obtain imaginary mode from the input force field 
c
c ---- read force constant matrix, which is returned
c      in cartesian coordinates 
c
      call secloc(isect(46),oexist,isec46)
      if(.not.oexist) 
     *   call caserr('no force matrix found on dumpfile')
      write(iw,9110) 
      call rdedx(hess,ihess,isec46,numdu)
      if(odebug) then
         write(iw,9004)
         do i=1,ncoord
            write(iw,'(10(f7.4,2x))') (hess(i,j),j=1,ncoord)
         enddo
      endif
c
c     ----- mass weight the hessian, store in symmetric form -----
c
      ij = 0
      do 130 i=1,ncoord
         rtrmsi = rtrms(i)
         do 120 j=1,i
            ij = ij + 1
            a(ij) = rtrmsi * hess(i,j) * rtrms(j)
  120    continue
  130 continue
c
c     ----- obtain normal modes and frequencies -----
c
      ierr = 0
c
c ---- use different gldia version in order to aviod problems
c      has to be cleaned up !!!!
c
c     call gldiag(ncoord,ncoord,ncoord,a,scr,e,vec,ia,ierr)
      call gldia(ncoord,ncoord,a,scr,e,vec,ierr,ia,iw)
      if(odebug) write(iw,9005) (e(i),i=1,ncoord)
      eone = e(1)
      if(eone.ge.szero) then
         write(iw,9020) eone
         call caserr('no imaginary mode found')
      end if
      freq = dsqrt(-fact*eone)
c
c ---- convert the imaginary mode from mass weighted
c      to ordinary cartesian coordinates 
c
      do 150 i=1,ncoord
         cmode(i) = rtrms(i) * vec(i,1)
  150 continue
      go to 300
c
c ---- normalize the input normal mode so its length
c      in mass weighted cartesian space is unity -----
c
  200 continue
      eone = -freq*freq/fact
      if(eone.ge.szero) 
     + call caserr('vibrational energy for imaginary mode > zero')
      vnorm = szero
      do 220 i=1,ncoord
         v = cmode(i)
         v = v/rtrms(i)
         vnorm = vnorm + v * v
  220 continue
      vnorm = dsqrt(vnorm)
      if(vnorm.lt.small) call caserr('norm of step too small')
      vnorm = sone/vnorm
      if(dabs(sone-vnorm).gt.small) write(iw,9100) vnorm
      do 240 i=1,ncoord
         cmode(i) = vnorm * cmode(i)
  240 continue
c
c ---- at this point, cmode is in cartesian coordinates,
c      with unit length in mass-weighted cartesian space.
c      eone is the corresponding negative curvature. 
c
c ---- pick direction along the imaginary mode,
c      and the magnitude of the harmonic step 
c
  300 continue
      step = dsqrt(-two * evib/eone)
      imax = 1
      vmax = dabs(cmode(1))
      do 320 i=2,ncoord
         v = dabs(cmode(i))
         if(v.le.vmax) go to 320
            imax = i
            vmax = v
  320 continue
      if(     oforw  .and.  cmode(imax).lt.szero) step = -step
      if(.not.oforw  .and.  cmode(imax).gt.szero) step = -step
c
c     ----- form a symmetrized step along the normal mode -----
c
      do 340 i=1,ncoord
         cmode(i) = step * cmode(i)
  340 continue
      call symdr(cmode,core)
      write(iw,9040) step,freq
      ii=1
      do 1220 i=1,nat
         write(iw,9060) (cmode(k),k=ii,ii+2)
         ii=ii+3
1220  continue
c
c ---- generate cartesian coordinates of the point that
c      is displaced from the saddle point 
c
      do 400 i=1,ncoord
         qb(i) = qa(i) + cmode(i)
  400 continue
      dist = szero
      call ircdst(qa,qb,rtrms,ncoord,dist)
      write(iw,9080) dist
      sab = dist
      stotal = stotal + dist
      return
c
 
 9000 format(//10x,'**** jumping off saddle point along the imaginary',
     *             ' normal mode'//)
 9002 format(//,10x,'    cartesian coordiantes of saddle point',/
     *       10x,'----------------------------------------------'/
     *       10x,' atom         x           y           z       '/
     *       10x,'----------------------------------------------')
 9003 format(10x,a8,2x,3(f10.6,2x))
 9004 format(//10x,'hessian restored from dumpfile',/
     *         10x,'------------------------------')
 9005 format(//10x,'eigen values of cartesian force matrix',/
     *         10x,'--------------------------------------',/
     *         10x,(f7.4,2x))
 9020 format(//10x,'**** the input hessian does not',
     *          ' have an imaginary mode, e(1)=',f16.6)
 9040 format(//10x,'initial displacement=',f6.2,' along mode with ',
     *           'frequency=',f12.3,' i cm**-1'///
     *       10x,'the components of the displacement are (bohr) -'/,
     *       10x,'----------------------------------------------')
 9060 format(10x,3(f11.6,4x))
 9080 format(//10x,'the distance from the saddle point to',
     *          ' the first irc point is',f10.5,' dsqrt(amu)-bohr)'//)
 9100 format(//10x,'*** input normal mode normalized by factor of',
     *f10.4)
 9110 format(//10x,'*** force matrix restored from dumpfile ***'//) 
      end
      subroutine iniirm(qx,iw)
      implicit real*8 (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      common /irc2 /oforw,osadd,ostab,otsen,ostbs,olimn,oline,oquad,
     *              oampc,orung,ovtst,ndummy,evib,delta,elbow,stride,
     *              npoint,nextpt,freq,ypac,yyy,ncoord,nhopt,nampc,
     *              ofirst,nn
c
      dimension qx(*)
c
c ---- grow memory for iniirc
c
      nc2 = (ncoord*ncoord + ncoord)/2
      nc3 = ncoord * ncoord
c
c ---- setup core partition for irc calculation
c
      i10 =   0
      i20 = i10 + nc3
      i30 = i20 + nc2
      i40 = i30 + ncoord
      i50 = i40 + ncoord * 8
      i60 = i50 + nc3
      need= i60 + ncoord
c
      i10 = igmem_alloc(need)
      i20 = i10 + nc3
      i30 = i20 + nc2
      i40 = i30 + ncoord
      i50 = i40 + ncoord * 8
      i60 = i50 + nc3
c
c
      write(iw,1009) need
1009  format(10x,'***** ',i7,' words of fast memory allocated for irc')
c
      call iniirc(nc2,qx,qx(i10),qx(i20),qx(i30),qx(i40),
     +  qx(i50),qx(i60),iw)
c
c ---- free memory used by irc step
c
      call gmem_free(i10)
c
      return
      end
      subroutine irc(q,iw)
      implicit real*8 (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
c
      parameter(szero=0.0d0,sone=1.0d0,nnam=23,mxa3=maxat*3)
c
c ...
c ... timing statistics
c ...
      real*8 begin, ebegin, timsec, walsec, tstart, estart
      common/statis/begin,ebegin,timsec(50),walsec(50),
     + tstart,estart
c
c
      integer nmaxly, ntotly, icurly, isecly
      logical oprintm
      common /segm/ nmaxly,ntotly,icurly,isecly(100),
     +              oprintm
c
      common/junkc /ylab(26),ztype(831),zoptio(13),ztagg(100)
      common/data1 /vlist(400),np(206),np1(maxorb,2),norbt(6),
     *              norbta(maxorb,2),norb2,nseco(maxorb),norb3,
     *              nthird(maxorb),norbg,norbga(maxorb),nscrp,
     *              range(2),alimit,blimit,isymm,nnx,lsym,nsymg(8),
     *              iuse(8,6),idorb(100),ncont,nmax(10)
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
c from March 97, amass is used only for user-input
c mass values
c czanr contains a copy of czan made in basis.m after reorganisation and
c before pseudo potential; so real charges in real order 
c beware, infob is explicit in direct.m and  dirrpa.m
c
      real*8 czann, czanr, cat, amass, cnew
      integer nonsym, map80
      logical ozmat
      common/infob/czanr(maxat),czann(maxat),cat(3,maxat),amass(maxat),
     +             cnew(maxat3,3),nonsym,map80(maxat),ozmat

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
      integer numdu, iblkdu, maxb, iblkla
      logical orevis
      real*8 apos
      common/sector/numdu,iblkdu,orevis(2),
     + apos(508),maxb,iblkla
      common /irc1 /ea,eb,ec,sab,sbc,stotal,ganorm,gbnorm,gcnorm,
     *              qa(mxa3),qb(mxa3),qc(mxa3),ga(mxa3),gb(mxa3),
     *              gc(mxa3),qa2(mxa3),qa1(mxa3),qa0(mxa3),ga2(mxa3),
     *              ga1(mxa3),ga0(mxa3),qcp(mxa3),qcc(mxa3),gcp(mxa3),
     *              bis(mxa3),ebis(10),cbis(10),gbis,rtrms(mxa3),
     *              ck0(mxa3),ck1(mxa3),ck2(mxa3),ck3(mxa3),
     *              ipirc,npirc
      common /irc2 /oforw,osadd,ostab,otsen,ostbs,olimn,oline,oquad,
     *              oampc,orung,ovtst,ndummy,evib,delta,elbow,stride,
     *              npoint,nextpt,freq,ypac,yyy,ncoord,nhopt,nampc,
     *              ofirst,nn
      character*10 charwall
      dimension q(*)
      begin=cpulft(1)
      write(iw,1100)
 1100 format(/1x,104('=')/
     *27x,36('*')/
     *27x,'intrinsic reaction coordinate module'/
     *27x,36('*')//)
      ncoord=3*nat
c
c ---- print atomic masses used, then get triples
c ---- of invert square roots of the masses 
c
      write(iw,9200)
      i0=0
c     nzmat=0
      do 100 i=1,nat
ccc
ccc   use atomic mass_get (NB original reference to 
ccc   infob array amass was incorrect if symmetry present
ccc
ccc         temp = sone/sqrt(amass(i))
ccc         write(iw,9220) zaname(i),amass(i),temp

         wght =  amass_get(1,i)
         temp = sone/sqrt(wght)
         write(iw,9220) zaname(i),wght,temp
         rtrms(i0+1)=temp
         rtrms(i0+2)=temp
         rtrms(i0+3)=temp
         i0=i0+3
  100 continue
c
c ---- move coords of input group to current point qb -----
c
      ii=0
      do 110 i=1,nat
         qb(ii+1) = c(1,i)
         qb(ii+2) = c(2,i)
         qb(ii+3) = c(3,i)
         ii=ii+3
  110 continue
c
c ----- generate a first point on the irc just off
c       the saddle point (bis holds cmode)
c
      if(osadd) then
         cmodln = dsqrt(ddot(ncoord,bis,1,bis,1))
         nhopt=1
         if(cmodln.gt.0.001) nhopt=0
         if(nhopt.eq.0) then
            write(iw,9300) freq
            write(iw,9320) (bis(i),i=1,ncoord)
         else
            write(iw,9340)
         endif
         if(otsen) then
            e=szero
c
c ---- go and calculate the energy
c
            call scfgvb(q)
c
c ---- this has to be modified in order to run casscf etc
c
            e=etot
            if(e.eq.szero) 
     *        call caserr(' zero energy calculated')
            ea = e
         end if
c
c ---- allocate fast memory and generate an inital point
c      with nozero gradient vector by following the imaginary
c      normal mode a short distance away from the saddle point
c
         call iniirm(q,iw)
         nextpt=1
      end if
c
c     ----- generate additional points on the irc -----

       npirc = nextpt
  200 continue
      call ircstp(q,iw)
      if(npirc-nextpt.lt.npoint) go to 200
      write(iw,9120) npoint
 9120 format(//10x,'**********************************',/,
     *         10x,'max number of points reached :',i4,/,
     *         10x,'**********************************'//)
      nextpt = npirc
c
c ---- print/punch predicted coordinates of next irc point 
c      the predicted next point is now in qb, qa holds the
c      last irc point where energy/gradient was computed.
c
      write(iw,9400)
      i0 = 0
      do 420 i=1,nat
         cznuc = czan(i)
         write(iw,9410) zaname(i),cznuc,(qb(i0+j),j=1,3)
         i0 = i0 + 3
  420 continue
      ii=0
      do 430 i=1,nat
         c(1,i) = qb(ii+1)
         c(2,i) = qb(ii+2)
         c(3,i) = qb(ii+3)
         ii=ii+3
  430 continue
      en = enucf(nat,czan,c)
      write(iw,9420) en
      call intr(q)
c     if(nzmat.gt.0) then
c         call bandbi
c         call pzandg(0,0)
c     end if
c
      write(iw,9430) ypac,stride,stotal,nextpt
c
      call ircdmp
c 
      if(ostab) write(iw,9440) delta,elbow
c 
      if(ypac.eq.'quad') write(iw,9450) sab
c
 9200 format(/10x,'atomic masses and rtrms used (in amu) are'/
     *        10x,'-----------------------------------------',/)
 9220 format(10x,a8,2x,f10.4,2x,f10.4)
c 
 9300 format(/10x,'jump away from saddle point along mode with',
     *       ' frequency',f11.5,' input cmode=')
 9320 format(10x,3f15.8,5x,3f15.8)
 9340 format(//10x,'jump away from saddle point will be found from',
     *       ' the input hessian matrix')
c
 9400 format(//10x,'**** irc restart information ',
     *           ' written to dumpfile   ****'//
     *       10x,'coordinates of the predicted point qb'/,
     *       10x,'-------------------------------------')
 9410 format(10x,a8,2x,f5.1,3f20.10)
 9420 format(//10x,
     *     'the nuclear repulsion at the next irc point is',f20.10//)
 9430 format(//10x,'short info'/
     +         10x,'----------'/
     +         10x,'pace  =',a8/
     +         10x,'stride=',f10.6/
     +         10x,'stotal=',f10.6/
     +         10x,'nextpt=',i8/
     +         10x,'saddle=.f.')
 9440 format(10x,'stablz=.t.'/
     *       10x,'delta =',f10.6/
     *       10x,'elbow =',f10.6)
 9450 format(10x,'sab   =',f10.6)
 1300 format(//1x,'end of irc module at',f8.2,' seconds',a10,' wall'/
     */1x,104('=')//)
      end=cpulft(1)
      write(iw,1300) end ,charwall()
      call clredx
      call timana(28)
      return
      end
      subroutine ircdat(iw)
      implicit real*8 (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c ============================
c = irc input driver
c ============================
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
c
      parameter(szero=0.0d0,sone=1.0d0,nnam=23,mxa3=maxat*3)
c
      common/junkc/ylab(26),ztype(831),zoptio(13),ztagg(100)
      common/data1/vlist(400),np(206),np1(maxorb,2),norbt(6),
     *             norbta(maxorb,2),norb2,nseco(maxorb),norb3,
     *             nthird(maxorb),norbg,norbga(maxorb),nscrp,
     *             range(2),alimit,blimit,isymm,nnx,lsym,nsymg(8),
     *             iuse(8,6),idorb(100),ncont,nmax(10)
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
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
      common /irc1 /ea,eb,ec,sab,sbc,stotal,ganorm,gbnorm,gcnorm,
     *              qa(mxa3),qb(mxa3),qc(mxa3),ga(mxa3),gb(mxa3),
     *              gc(mxa3),qa2(mxa3),qa1(mxa3),qa0(mxa3),ga2(mxa3),
     *              ga1(mxa3),ga0(mxa3),qcp(mxa3),qcc(mxa3),gcp(mxa3),
     *              bis(mxa3),ebis(10),cbis(10),gbis,rtrms(mxa3),
     *              ck0(mxa3),ck1(mxa3),ck2(mxa3),ck3(mxa3),
     *              ipirc,npirc
      common /irc2 /oforw,osadd,ostab,otsen,ostbs,olimn,oline,oquad,
     *              oampc,orung,ovtst,ndummy,evib,delta,elbow,stride,
     *              npoint,nextpt,freq,ypac,yyy,ncoord,nhopt,nampc,
     *              ofirst,nn
c
      character*80 pace
c
      dimension zircip(15),zpace(10),zpace1(10),zpace2(10),zpace3(10)
      dimension kqnam(nnam)
c
      data ircmax/12/
      data ipace / 5/
      data ipace1,ipace2,ipace3/8,2,5/
      data zirc/'irc'/
      data zend/'end'/
      data zircip/
     *     'pace    ','stride  ','npoints ','nextpnt ','s-total ',
     *     'saddle  ','tsenergy','forward ','enervib ','frequen ',
     *     'c-mode  ','vtst    ','        ','        ','        '/
c
      data zpace /
     *     'linear  ','quad-grd','ampc4   ','runge-ku','        ',
     *   5*'        ' /
      data zpace1/
     *     'imk-stab','delta   ','elbow   ','readqb  ','        ',
     *     'eb      ','gbnorm  ','gb      ','        ','        '/
      data zpace2/
     *     's-ab    ','grad-vec',8*'        '/
      data zpace3/
     *     'grd0    ','grd1    ','grd2    ','grda    ','grdb    ',
     *   5*'        '/
      data kqnam/0,0,0,3,3,3,3,3,3,-3,1,1,5,0,3,-3,0,3,3,-3,
     *           -3,-3,-3/
c     data zblank/10*'        '/
c
c ---- defaults for irc
c
      oline=.true.
      oquad=.false.
      oampc=.false.
      orung=.false.
      ostab=.false.
      osadd=.false.
      otsen=.false.
      oforw=.true.
      ostbs=.false.
      ovtst=.false.
      ofirst=.true.
      evib=0.0005d0
      delta=0.025d0
      elbow=175.0d0
      stride=0.15d0
      npoint=1
      nextpt=0
c
c ---- fortran channel for punching results for vtst
c
      ipirc=50
c
      pace='linear gradient following selected (eulers method)'
      ypac=ytrunc(zpace(1))
c
      if(zruntp.ne.zirc )
     *call caserr('faulty directive ordering')
c
      top=cpulft(1)
      write(iw,444)top
 444  format(//1x,104('=')//
     *' **** irc-input processor called at',f9.2,' secs')
1001   call input
      call inpa(ztest)
      if(ztest.eq.zend) go to 1000
      ytest=ytrunc(ztest)
      k=locats(zircip,ircmax,ytest)
      if(k.eq.0) goto 999
      goto(10,20,30,40,50,60,70,80,90,100,110,120,999,999,999)k
c
c ---- input integration method
c
10    call inpa(ztest)
      ytest=ytrunc(ztest)
      l=locats(zpace,ipace,ytest)
      if ( l.eq.0 ) 
     *   call caserr('invalid integration method') 
         go to ( 101,102,103,104) l
c
c ---- linear gradient following selected (euler's method) default
c     
101   call inpa(ztest)
      pace='linear gradient following selected (eulers method)'
      ypac=ytrunc(zpace(1))
      oline=.true.
      if(ztest.eq.zend) go to 1001
      ytest=ytrunc(ztest)
      l=locats(zpace1,ipace1,ytest)
      if(l.eq.0)
     *  call caserr('invalid input for integration scheme: linear')
      if (l.eq.1) then
         ostab=.true.
         goto  101
      elseif (l.eq.2) then
         call inpf(delta)
         go to 101
      elseif (l.eq.3) then
         call inpf(elbow)
         go to 101
      elseif (l.eq.4) then
         call inpf(eb)
         call inpf(gbnorm)
         call inpf(gb)
         oread=.true.
         goto 101
      endif
c
c ---- quadratic gradient following selected
c
102   oquad=.true.
      pace='quadratic gradient following selected'
      ypac=ytrunc(zpace(2))
      call inpa(ztest)
      if(ztest.eq.zend) go to 1001
      ytest=ytrunc(ztest)
      l=locats(zpace2,ipace2,ytest)
      if(l.eq.0)
     * call caserr('invalid input for integration scheme: quadratic')
      if(l.eq.1) then
        call inpf(sab)
        go to 102
      elseif (l.eq.2) then
        call inpf(ga)
        goto 102
      endif
c
c ---- fourth order adams-moulton variable step predictor-corrector
c
103   oampc=.true.
      pace=
     * 'fourth order adams-moulton variable step predictor-corrector'
      ypac=ytrunc(zpace(3))
      call inpa(ztest)
      if(ztest.eq.zend) go to 1001
      ytest=ytrunc(ztest)
      l=locats(zpace3,ipace3,ytest)
      if(l.eq.0)
     *  call caserr('invalid input for integration scheme: ampc4')
      if (l.eq.1) then
         call inpf(ga0)
         go to 103
      elseif (l.eq.2) then
         call inpf(ga1)
         go to 103
      elseif (l.eq.3) then
         call inpf(ga2)
         go to 103
      endif
c
c ---- 4th order runge-kutta variable step method selected
c
104   orung=.true.
      pace='4th order runge-kutta variable step method'
      ypac=ytrunc(zpace(4))
      go to 1001

c
c ---- stride is used to calculate the step taken
c
20    call inpf(stride)
      go to 1001
c
c ---- npoint number of points to be located in this run
c
30    call inpi(npoint)
      go to 1001
c
c ---- nextpt the number of the next point to be computed
c
40    call inpi(nextpt)
      go to 1001
c
c ---- stotal total distance along the reac.path to the next irc point
c
50    call inpf(stotal)
      go to 1001
c
c ---- saddle point ? or poin on the irc                 
c
60    osadd=.true.
      go to 1001
c 
c ---- tsenergy is energy evaluated at transistion state ?
c
70    otsen=.true. 
      go to 1001
c
c ---- forwrd direction to proceed away from the ts
c
80    oforw=.false.
      go to 1001
c
c --- evib decrease in energy when following the imag. normal mode away from ts
c
90    call inpf(evib)
      go to 1001
c
c --- frequency magnitude of the imagniary freq. given in cm**-1
c 
100   call inpf(freq)
      go to 1001
c
c --- cmode array of the components of the normal mode whose frequency is
c           imaginary, in cartesian coordinates 
110   call caserr('cmode not yet implemented ')
      go to 1001
c
c --- vtst open punch file for further use of truhlar's polyrate program
c
120   ovtst=.true.
      open(unit=ipirc,file='vtst.dat',status='unknown')
      go to 1001
c
999   call caserr('invalid input option in irc input driver') 
c
1000  continue
c
      do 1900 i=1,nnam
         if(kqnam(i).eq.-3) kqnam(i) = 10*mxa3 + 3
 1900 continue
c
      ncoord = 3 * nat
      ea = szero
      eb = szero
      ec = szero
      sab = szero
      sbc = szero
      stotal = szero
      ganorm = szero
      gbnorm = szero
      gcnorm = szero
      gbis = szero
      do 1330 i=1,mxa3
         qa2(i) = szero
         qa1(i) = szero
         qa0(i) = szero
         qa(i) = szero
         qb(i) = szero
         qc(i) = szero
         ga2(i) = szero
         ga1(i) = szero
         ga0(i) = szero
         ga(i) = szero
         gb(i) = szero
         gc(i) = szero
         bis(i) = szero
         rtrms(i) = szero
 1330 continue
      do 1440 i=1,10
         cbis(i) = szero
         ebis(i) = szero
 1440 continue
c
      freq = szero
c                      bis is used as temp storage of cmode
      do 1650 i=1,ncoord
         bis(i) = szero
 1650 continue
c
      jret = 0
      evib = dabs(evib)
      stride = dabs(stride)
c
      if(ostab.and.ypac.ne.ytrunc(zpace(1))) 
     *call caserr('stabilization is available for linear pace only')
      if(ypac.eq.ytrunc(zpace(4))) nampc = 1
      if(jret.eq.2) call caserr('jret eq 2')
c
      write(iw,9000) pace,stride,npoint,nextpt,stotal,sab
      write(iw,9020) osadd,otsen,oforw,evib,freq
      write(iw,9040) ostab,elbow,delta,oread,eb,gbnorm,ovtst
c
      top=cpulft(1)
      write(iw,445)top
c
 445  format(//1x,
     *' **** irc-input processor finished at',f9.2,' secs',//
     *104('=')//)
 9000 format(/10x,29(1h*),/
     *        10x,'intrinsic reaction coordinate',/
     *        10x,'     input parameters'/
     *        10x,29(1h*),//
     *        10x,'integration method :',a80/
     *        10x,50('-'),/
     *        10x,'stride between points       :',f8.5,
     *                                          ' bohr*sqrt(amu)'/
     *        10x,'number of irc points        :',i8/
     *        10x,'number of next point        :',i8/
     *        10x,'total distance to the irc   :',f8.4,
     *                                          ' bohr*sqrt(amu)'/
     *        10x,'distance to previous point  :',f8.4,
     *                                          ' bohr*sqrt(amu)')
 9020 format(/10x,'saddle point given          :',l8/
     *        10x,'energy of transition state  :',l8/
     *        10x,'direction away form ts      :',l8/
     *        10x,'energy decrease away from ts:',f8.5,' hartree'/
     *        10x,'amplitude of imag. frequency:',f8.1,' cm**-1')
 9040 format(/10x,'ihsida/morokuma/komornicki  :',l8/
     *        10x,'collinearity threshold      :',f8.2,' degrees'/
     *        10x,'initial step size           :',f8.5,' bohr'//
     *        10x,'read qb                     :',l8/
     *        10x,'energy b                    :',f20.10,' hartrees'/
     *        10x,'gradient norm b             :',f20.10,
     *                             ' hartrees/sqrt(amu)-bohr'//
     *        10x,'punching results for vtst   :',l8,/,10x,29(1h*)/)
      return
      end
      subroutine ircdmp
      implicit real*8 (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
c
      parameter(mxa3=maxat*3)
c
      common /irc1 /ea,eb,ec,sab,sbc,stotal,ganorm,gbnorm,gcnorm,
     *              qa(mxa3),qb(mxa3),qc(mxa3),ga(mxa3),gb(mxa3),
     *              gc(mxa3),qa2(mxa3),qa1(mxa3),qa0(mxa3),ga2(mxa3),
     *              ga1(mxa3),ga0(mxa3),qcp(mxa3),qcc(mxa3),gcp(mxa3),
     *              bis(mxa3),ebis(10),cbis(10),gbis,rtrms(mxa3),
     *              ck0(mxa3),ck1(mxa3),ck2(mxa3),ck3(mxa3),
     *              ipirc,npirc
      common /irc2 /oforw,osadd,ostab,otsen,ostbs,olimn,oline,oquad,
     *              oampc,orung,ovtst,ndummy,evib,delta,elbow,stride,
     *              npoint,nextpt,freq,ypac,yyy,ncoord,nhopt,nampc,
     *              ofirst,nn
      integer numdu, iblkdu, maxb, iblkla
      logical orevis
      real*8 apos
      common/sector/numdu,iblkdu,orevis(2),
     + apos(508),maxb,iblkla
c 
      dimension dummy(mxa3),ia(50)
c
      data ircsec,irctyp/250,22/
c 
      odebug=.false.
c
      do 110 i=1,22
         ia(i)=ncoord*(i-1)+10
 110  continue
      ia(23)=ia(22)+10
      ia(24)=ia(23)+10
c
      dummy(1)=ea     
      dummy(2)=eb
      dummy(3)=ec    
      dummy(4)=sab
      dummy(5)=sbc
      dummy(6)=stotal
      dummy(7)=ganorm
      dummy(8)=gbnorm
      dummy(9)=gcnorm
      dummy(10)=gbis
      do 100 i=1,ncoord
          dummy(ia( 1)+i)=   qa(i)
          dummy(ia( 2)+i)=   qb(i)
          dummy(ia( 3)+i)=   qc(i)
          dummy(ia( 4)+i)=   ga(i)
          dummy(ia( 5)+i)=   gb(i)
          dummy(ia( 6)+i)=   gc(i)
          dummy(ia( 7)+i)=  qa2(i)
          dummy(ia( 8)+i)=  qa1(i)
          dummy(ia( 9)+i)=  qa0(i)
          dummy(ia(10)+i)=  ga2(i)
          dummy(ia(11)+i)=  ga1(i)
          dummy(ia(12)+i)=  ga0(i)
          dummy(ia(13)+i)=  qcp(i)
          dummy(ia(14)+i)=  qcc(i)
          dummy(ia(15)+i)=  qcp(i)
          dummy(ia(16)+i)=  bis(i)
          dummy(ia(17)+i)=rtrms(i)
          dummy(ia(18)+i)=  ck0(i)
          dummy(ia(19)+i)=  ck1(i)
          dummy(ia(20)+i)=  ck2(i)
          dummy(ia(21)+i)=  ck3(i)
 100  continue
      do 120 i=1,10
         dummy(ia(22)+i)=  ebis(i)
         dummy(ia(23)+i)=  cbis(i)
 120  continue 
      if(odebug) then
         do i=1,ia(24)
            write(6,*) 'dummy(',i,')=',dummy(i)
         enddo 
      endif
      length=ia(24)/512+1
      call secput(ircsec,irctyp,length,iblk)
      call wrt3(dummy ,ia(24),iblk  ,numdu) 
      call revind
      return
      end
      subroutine ircdst(a,b,rtrms,ncoord,dist)
      implicit real*8 (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      parameter(szero=0.0d0)
c
c
      dimension a(ncoord),b(ncoord),rtrms(ncoord)
c
c ---- compute mass-weighted distance from point a to b 
c      rtrms contains inverse square root of masses
c
      dist = szero
      do 100 i=1,ncoord
         diff = (a(i) - b(i))/rtrms(i)
         dist = dist + diff * diff
  100 continue
      dist = dsqrt(dist)
      return
      end
      subroutine ircstp(core,iw)
      implicit real*8 (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
c
      parameter(szero=0.0d+00,sone=1.0d+00, two=2.0d+00, three=3.0d+00
     *         ,sone80=180.0d+00, gtol=1.0d-08, tol5=5.0d-06,
     *          nnam=23,mxa3=maxat*3)
c
c ...
c ... timing statistics
c ...
      real*8 begin, ebegin, timsec, walsec, tstart, estart
      common/statis/begin,ebegin,timsec(50),walsec(50),
     + tstart,estart
c
c
      integer nmaxly, ntotly, icurly, isecly
      logical oprintm
      common /segm/ nmaxly,ntotly,icurly,isecly(100),
     +              oprintm
c
      common/junkc /ylab(26),ztype(831),zoptio(13),ztagg(100)
      common/data1 /vlist(400),np(206),np1(maxorb,2),norbt(6),
     *              norbta(maxorb,2),norb2,nseco(maxorb),norb3,
     *              nthird(maxorb),norbg,norbga(maxorb),nscrp,
     *              range(2),alimit,blimit,isymm,nnx,lsym,nsymg(8),
     *              iuse(8,6),idorb(100),ncont,nmax(10)
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
c from March 97, amass is used only for user-input
c mass values
c czanr contains a copy of czan made in basis.m after reorganisation and
c before pseudo potential; so real charges in real order 
c beware, infob is explicit in direct.m and  dirrpa.m
c
      real*8 czann, czanr, cat, amass, cnew
      integer nonsym, map80
      logical ozmat
      common/infob/czanr(maxat),czann(maxat),cat(3,maxat),amass(maxat),
     +             cnew(maxat3,3),nonsym,map80(maxat),ozmat

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
      common/miscop/paa(mxa3),grd(mxa3)
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
c
      common /irc1 /ea,eb,ec,sab,sbc,stotal,ganorm,gbnorm,gcnorm,
     *              qa(mxa3),qb(mxa3),qc(mxa3),ga(mxa3),gb(mxa3),
     *              gc(mxa3),qa2(mxa3),qa1(mxa3),qa0(mxa3),ga2(mxa3),
     *              ga1(mxa3),ga0(mxa3),qcp(mxa3),qcc(mxa3),gcp(mxa3),
     *              bis(mxa3),ebis(10),cbis(10),gbis,rtrms(mxa3),
     *              ck0(mxa3),ck1(mxa3),ck2(mxa3),ck3(mxa3),
     *              ipirc,npirc
      common /irc2 /oforw,osadd,ostab,otsen,ostbs,olimn,oline,oquad,
     *              oampc,orung,ovtst,ndummy,evib,delta,elbow,stride,
     *              npoint,nextpt,freq,ypac,yyy,ncoord,nhopt,nampc,
     *              ofirst,nn
c
c
c
      real*8 enrgy, egrad
      common/gms_funct/enrgy,egrad(3*maxat)
c
c
      dimension core(*)
c     
c
c ---- jump from point qb on the irc to a point qc by
c      a linear (possibly stabilized) or quadratic pace 
c
c ---- if a linear pace with stabilization has been chosen,
c      and the path is not curved very much, we may have
c      skipped the stabilization, and already know eb, gb ---
c
      cpu=cpulft(1)
      if(cpu.ge.timlim-safety*1.2d0) then
            write(iw,9455)
            npoint=npirc-nextpt
            return
      endif
      if(ostbs) write(iw,9010) npirc,stotal
      if(ostbs) go to 300 
c
c     ----- print coordinates of point qb -----
c
      cpu=cpulft(1)
      write(iw,9000) npirc,cpu,stotal
c     write(iw,9040)
c     i0 = 0
c     do 120 i=1,nat
c        write(iw,9060) zaname(i),(qb(i0+j),j=1,3)
c        i0 = i0 + 3
c 120 continue
c
c     ----- energy at point qb -----
c
      ii=0
      do 130 i=1,nat
         c(1,i) = qb(ii+1)
         c(2,i) = qb(ii+2)
         c(3,i) = qb(ii+3)
         ii=ii+3
  130 continue
c
c ---- calculate interatomic distances
c
c     call intr(core)
      enrgy = szero
c
c ---- energy and gradient at point qb -----
c
      if(ofirst) then
         call gamgn2
         ofirst=.false.
      endif
      nprint=-5 
      call optx(core)
      enrgy=etot
      if(enrgy .eq. szero ) call caserr('zero energy evaluated')
      eb = enrgy
c
      if(eb.gt.ea) then
        write(iw,9090) eb,ea
        npirc=99999 
        return
      end if
      do 200 i=1,ncoord
         egrad(i) = szero
 200  continue
      do 220 i=1,ncoord
         gb(i) = grd(i)
         egrad(i) = grd(i)
 220  continue
c
c     ----- print internal coordinates, with gradient -----
c
c     if(nzmat.gt.0) then
c        call bandbi
c        call trang(egrad,nvar,ncoord)
c        call pzandg(egrad,1)
c     end if

c
c     ----- normalized mass weighted gradient at point qb -----
c
      gbnorm = szero
      do 240 i=1,ncoord
          g = rtrms(i) * gb(i)
          if(dabs(g).lt.gtol) g=szero
          gb(i) = g
          gbnorm = gbnorm + g * g
 240  continue
      gbnorm = dsqrt(gbnorm)
      write(iw,9100) gbnorm
      fact = sone/gbnorm
      do 250 i=1,ncoord
         gb(i) = fact * gb(i)
 250  continue
      write(iw,9120)
       i0 = 0
       do 260 i=1,nat
          write(iw,9140) zaname(i),(gb(i0+j),j=1,3)
          i0 = i0 + 3
 260  continue
  300 continue
c
      if(ovtst) then
c
c     ----- punch data from point qb on irc for possible vtst use -----
c
         ostbs=.false.
         write(ipirc,8000) npirc,stotal,eb
         i0=0
         do 310 i=1,nat
c           cznuc = czan(i) + izcore(i)
            cznuc = czan(i) 
            write(ipirc,8020) zaname,cznuc,(c(j,i),j=1,3)
  310    continue
         write(ipirc,8040) gbnorm
         i0=0
         do 320 i=1,nat
c           znuc = zan(i) + izcore(i)
            cznuc = czan(i) 
            write(ipirc,8020) zaname,cznuc,(gb(i0+j),j=1,3)
            i0=i0+3
  320    continue
      endif
c
c     ----- step to point qc by the appropriate pace -----
c
      nprint = -5
      if(ypac.eq.'rung')   then
        zpace='rk4     '
        write(iw,1039) zpace
        go to 660
      elseif(ypac.eq.'ampc') then
        zpace='ampc4   '
        write(iw,1039) zpace
        go to 550
      elseif(osadd .and.npirc.eq.1) then
        if(ypac.eq.'line') zpace='linear  '
        if(ypac.eq.'quad') zpace='quadrat '
        write(iw,1039) zpace
        go to 400
      elseif(ypac.eq.'line') then
        zpace='linear  '
        write(iw,1039) zpace
        go to 400
      elseif(ypac.eq.'quad') then
        zpace='quadrat '
        write(iw,1039) zpace
        go to 500
      else
        call caserr('invalid pace given (should not occur)')
      endif
 1039 format(//10x,'**** jumping to point qc via a ',a,' pace ***'//)
c
c
c     ----- linear gradient following -----
c
  400 continue
      a = -stride
      do 420 i=1,ncoord
         qc(i) = a * gb(i)
  420 continue
      call symdr(qc,core)
      do 440 i=1,ncoord
         qc(i) = qb(i) + qc(i)*rtrms(i)
  440 continue
      if(ostab) then
         call stbimk(core,iw)
c
c ---- if not jumping saddle point, calculate angle between
c
         if (npirc .gt. 1) then
            dotprd = szero
            do 450 i=1,ncoord
               dotprd = dotprd + ga(i)*gb(i)
  450       continue
            theta = dacos(dotprd)
            pi = dacos(-sone)
            theta = sone80*theta/pi
            theta = sone80 - theta
            write(iw,9200) theta
         end if
      end if
      go to 700
c
c ----- quadratic gradient following 
c
  500 continue
      a = -stride
      b = (stride*stride)/(two*sab)
      do 510 i=1,ncoord
         qc(i) = (a-b)*gb(i) - b*ga(i)
  510 continue
      call symdr(qc,core)
      do 520 i=1,ncoord
         qc(i) = qb(i) + qc(i)*rtrms(i)
  520 continue
      go to 700
c
c     ----- adams-moulton predictor-corrector order four
c           the information needed to start this method
c           is obtained by doing 1 step using rk4.
c
  550 continue
      imark = 0
      h = -stride
      go to (660,564) nampc
c
  564 continue
      do 578 i=1,ncoord
         qc(i) = (h/24.0d0)*(55.0d0*gb(i)-59.0d0*ga0(i)+37.0d0*ga1(i)-
     *                        9.0d0*ga2(i))
  578 continue
c
      call symdr(qc,core)
      do 580 i=1,ncoord
         qcp(i) = qb(i) + qc(i)*rtrms(i)
  580 continue
c
c ---- energy and gradient at predicted point
c
      ii=0
      do 590 i=1,nat
         c(1,i) = qcp(ii+1)
         c(2,i) = qcp(ii+2)
         c(3,i) = qcp(ii+3)
         egrad(ii+1) =  szero
         egrad(ii+2) =  szero
         egrad(ii+3) =  szero
         ii=ii+3
  590 continue
      enrgy = szero
      call optx(core)
      enrgy=etot
      do 592 i=1,ncoord
         egrad(i)=grd(i)
  592 continue
      if(enrgy .eq. szero ) call caserr('zero energy evaluated')
c
c     ----- normalized mass weighted gradient at predicted point -----
c
      gbnorm = szero
      do 610 i=1,ncoord
         g = rtrms(i) * egrad(i)
         if(dabs(g).lt.gtol) g=szero
         egrad(i) = g
         gbnorm = gbnorm + g * g
  610 continue
      gbnorm = dsqrt(gbnorm)
      fact = sone/gbnorm
      do 620 i=1,ncoord
         gcp(i) = fact * egrad(i)
  620 continue
c
c     ----- corrector
c
      do 630 i=1,ncoord
         qc(i) = (h/24.0d0)*(9.0d0*gcp(i)+19.0d0*gb(i)-5.0d0*ga0(i)
     +            +ga1(i))
  630 continue
      call symdr(qc,core)
      do 640 i=1,ncoord
         qcc(i) = qb(i) + qc(i)*rtrms(i)
  640 continue
c
c     ----- variable stepper
c
      diff = szero
      sigma = szero
      do 650 i=1,ncoord
        diff = (qcp(i) - qcc(i))/rtrms(i)
         sigma = sigma + diff*diff
  650 continue
      sigma = dsqrt(sigma)/14.0d0
      power = 1.0d0/6.0d0
      hcorr = ((tol5/3.0d0)/sigma)**power
      if((sigma.lt.tol5).or.(imark.eq.1)) go to 700
      imark = 1
      h = hcorr*h
      go to 564
c
c     ----- runge-kutta order four
c
  660 continue
      nrk4 = 1
      h = -stride
  665 continue
         go to (672,674,676,678) nrk4
c
  672    continue
         do 682 i=1,ncoord
            ck0(i) = h*gb(i)
            qc(i) = 0.5d0*ck0(i)
  682    continue
         go to 670
c
  674    continue
         do 684 i=1,ncoord
            if(ypac.eq.'ampc') ga2(i) = egrad(i)
            ck1(i) = h*egrad(i)
            qc(i) = 0.5d0*ck1(i)
  684    continue
         go to 670
c
  676    continue
         do 686 i=1,ncoord
            if(ypac.eq.'ampc') ga1(i) = egrad(i)
            ck2(i) = h*egrad(i)
            qc(i) = ck2(i)
  686    continue
         go to 670
c
  678    continue
         do 688 i=1,ncoord
            if(ypac.eq.'ampc') ga0(i) = egrad(i)
            ck3(i) = h*egrad(i)
            qc(i) = 1.0d0/6.0d0*(ck0(i)+2.0d0*ck1(i)+2.0d0*ck2(i)
     +              +ck3(i))
  688    continue
c
  670 continue
c
      nrk4 = nrk4 + 1
      call symdr(qc,core)
      do 690 i=1,ncoord
         qc(i) = qb(i) + qc(i)*rtrms(i)
  690 continue
      if(nrk4.gt.4) go to 700
c
c    ----- energy and gradient at intermediate r-k point -----
c
      ii=0
      do 750 i=1,nat
          c(1,i) = qc(ii+1)
          c(2,i) = qc(ii+2)
          c(3,i) = qc(ii+3)
          egrad(ii+1) = szero
          egrad(ii+2) = szero
          egrad(ii+3) = szero
          ii=ii+3
  750 continue
      enrgy = szero
      write(iw,1299) nrk4
1299  format(//10x,
     * 'energy and gradient at intermediate r-k point: ',i3,//)
      call optx(core)
      enrgy=etot
      do 760 i=1,ncoord
         egrad(i)=grd(i)
 760  continue
c
      if(enrgy .eq. szero ) call caserr('zero energy evaluated') 
c
c     ----- normalized mass weighted gradient at predicted point -----
c
      gbnorm = szero
      do 694 i=1,ncoord
         g = rtrms(i) * egrad(i)
         if(dabs(g).lt.gtol) g=szero
         egrad(i) = g
         gbnorm = gbnorm + g * g
  694 continue
      gbnorm = dsqrt(gbnorm)
      fact = sone/gbnorm
      do 696 i=1,ncoord
         if(nrk4.le.4) then
            egrad(i) = fact * egrad(i)
         else
            gb(i) = fact * egrad(i)
         endif
  696 continue
      if(nrk4.le.4) then
         goto 665
      else
         nrk4 = 0
      endif
c
c     ----- qc is the predicted next point on the irc -----
c
  700 continue
      if((ypac.eq.'ampc4').and.(nampc.gt.1)) then
         call ircdst(qb,qcc,rtrms,ncoord,sbc)
      else
         call ircdst(qb,qc,rtrms,ncoord,sbc)
      endif
      stotal = stotal + sbc
      write(iw,9160) sbc,stotal
      write(iw,9180)
c
c     ----- current point qb becomes old point qa,
c           next point qc becomes current point qb. -----
c
      npirc = npirc+1
      ea = eb
      eb = ec
      sab = sbc
      ganorm = gbnorm
      gbnorm = gcnorm
      do 740 i=1,ncoord
         if(ypac.eq.'ampc4') then
           if(nampc.eq.1) then
              qb(i)  = qc(i)
              gb(i)  = gc(i)
           else
              ga2(i) = ga1(i)
              qa2(i) = qa1(i)
              ga1(i) = ga0(i)
              qa1(i) = qa0(i)
              ga0(i) = gb(i)
              qa0(i) = qb(i)
              qb(i) = qcc(i)
           endif
         else
            qa(i) = qb(i)
            qb(i) = qc(i)
            ga(i) = gb(i)
            gb(i) = gc(i)
         endif
  740 continue
      if(nampc.eq.1) nampc=nampc+1
      return
c
 8000 format('***** begin irc information packet *****'/
     *    'point=',i4,' stotal=',f10.5,'  e=',f20.10/
     *       'cartesian coordinates (bohr)')
 8020 format(a8,2x,f5.1,3f20.10)
 8040 format('mass-weighted gradient - orig.norm=',1p,e13.6,
     *       ' (hartree/bohr-sqrt(amu))')
c
 9000 format(//10x,75(1h*),/,
     *       10x,'calculation for point',i4,' on the reaction path',
     *           ' starting at ',f8.2,' s'/10x,75(1h*)/
     *       10x,'at path distance stotal=',f10.5//
     *       10x,'the coordinates of point qb are (bohrs)')
 9010 format(1x,'because of the elbow criteria, the previous geometry',
     *          ' is irc point',i5/1x,'at path distance stotal=',f10.5)
c9040 format(/10x,'atom',16x,'x',14x,'y',14x,'z'/
c    *        10x,59(1h-)/)
c9060 format(10x,a,4x,3f15.10)
c9080 format(1x,'**** error, scf failed to converge')
 9090 format(//10x,59(1h*)/,
     *       10x,'current energy=',f20.10,', is higher than the'/
     *       10x,'   last energy=',f20.10/
     *       10x,'you are no longer on a steepest descent path!'/
     *       10x,'this irc should be run with smaller steps and/or',
     *          ' a better integration method.'/10x,59(1h*)//)
 9100 format(//10x,'normalized mass weighted gradient at point qb'/
     *   10x,'(original gbnorm=',f10.6,' hartree/bohr-sqrt(amu))')
 9120 format(//10x,'atom',12x,'de/dx',7x,'de/dy',7x,'de/dz'/
     *        10x,46(1h-))
 9140 format(10x,a8,2x,3(2x,f10.6))
c9150 format(//10x,'*** jumping to point qc via a ',a8,' pace ***')
 9160 format(//10x,'the distance from qb to qc (or qc'') is',
     *          f10.5,' dsqrt(amu)-bohr)'/
     *       10x,' total path distance to qc (or qc'') is',
     *          f10.5,' sqrt(amu)-bohr)')
 9180 format(//,10x,'......next point on irc found......',
     *          'renaming qc(qc'')-->qb-->qa.')
 9200 format(//10x,'the angle between gradients ga and gb is',f8.2,
     *         ' degrees')
 9455 format(//10x,'***************************************',/,
     *         10x,'not enough time to calculate next point ',/,
     *         10x,'writing restart information on dumpfile'/,
     *         10x,'***************************************'//)
      end
      subroutine qfit1v(e1,g1,e2,c0,c1,c2,step,iwr,ofit)
      implicit real*8 (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c ---- fit parabola to energy and gradient at qc,
c      and energy at step along vector bis -----
c
      c0 = e1
      c1 = g1
      c2 = ((e2-e1) - step*g1)/(step*step)
      ofit = .true.
c
c     ----- test fit for reliability -----
c
      cxmin = -c1/(c2+c2)
      if(cxmin.gt.(step+step)) ofit = .false.
      if(cxmin.lt.-step)       ofit = .false.
      if(ofit) return
c
      emin = c2 * cxmin*cxmin  + c1 * cxmin  + c0
      write(iwr,9000) cxmin,emin,c0,c1,c2
      return
c
 9000 format(10x,'parabolic fit to 2 energies and 1 gradient',
     *      ' rejected due to large xmin.'/
     *    10x,'xmin,emin=',f12.6,f16.9,2x,
     *     'c0,c1,c2=',3f16.10/
     *   10x,'explicit linear search will be performed')
      end
      subroutine seabis(core,iw)
      implicit real*8 (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
c
      parameter(szero=0.0d+00,mxa3=maxat*3)
c
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
c from March 97, amass is used only for user-input
c mass values
c czanr contains a copy of czan made in basis.m after reorganisation and
c before pseudo potential; so real charges in real order 
c beware, infob is explicit in direct.m and  dirrpa.m
c
      real*8 czann, czanr, cat, amass, cnew
      integer nonsym, map80
      logical ozmat
      common/infob/czanr(maxat),czann(maxat),cat(3,maxat),amass(maxat),
     +             cnew(maxat3,3),nonsym,map80(maxat),ozmat

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
      common/miscop/paa(mxa3),grd(mxa3)
c
      common /irc1 /ea,eb,ec,sab,sbc,stotal,ganorm,gbnorm,gcnorm,
     *              qa(mxa3),qb(mxa3),qc(mxa3),ga(mxa3),gb(mxa3),
     *              gc(mxa3),qa2(mxa3),qa1(mxa3),qa0(mxa3),ga2(mxa3),
     *              ga1(mxa3),ga0(mxa3),qcp(mxa3),qcc(mxa3),gcp(mxa3),
     *              bis(mxa3),ebis(10),cbis(10),gbis,rtrms(mxa3),
     *              ck0(mxa3),ck1(mxa3),ck2(mxa3),ck3(mxa3),
     *              ipirc,npirc
      common /irc2 /oforw,osadd,ostab,otsen,ostbs,olimn,oline,oquad,
     *              oampc,orung,ovtst,ndummy,evib,delta,elbow,stride,
     *              npoint,nextpt,freq,ypac,yyy,ncoord,nhopt,nampc,
     *              ofirst,nn
c
c
      real*8 enrgy, egrad
      common/gms_funct/enrgy,egrad(3*maxat)
c
c
c
c     ----- search along the line bis for a minimum energy.
c           this minimum will be point qc' on the irc -----
c
      np = 1
      cbis(np) = szero
      ebis(np) = ec
      step = delta
c
c     ----- begin linear search -----
c
  100 continue
      np = np + 1
      if(np.gt.10) write(iw,9000)
      if(np.gt.10) call caserr('more than 10 point in seabis')
      ii=0
      do 120 i=1,nat
         c(1,i) = qc(ii+1) + step * bis(ii+1)
         c(2,i) = qc(ii+2) + step * bis(ii+2)
         c(3,i) = qc(ii+3) + step * bis(ii+3)
         ii=ii+3
  120 continue
      write(iw,9020) np,step
      write(iw,9040)
      ii=0
      do 140 i=1,nat
         write(iw,9060) zaname(i),(c(j,i),j=1,3)
         ii=ii+3
  140 continue
c
c     ----- energy evaluation -----
c
      enrgy = szero
      call scfgvb(core)
      enrgy=etot
      if(enrgy.eq.szero) write(iw,9080)
      if(enrgy.eq.szero) call caserr('zero energy evaluated')
      ebis(np) = enrgy
      cbis(np) = step
      ofit = .false.
      c0 = szero
      c1 = szero
      c2 = szero
      if(np.eq.2) call qfit1v(ec,gbis,enrgy,c0,c1,c2,step,iw,ofit)
      if(ofit) go to 320
      if(ebis(np).gt.ebis(np-1)) go to 160
      step = step + step
      go to 100
c
  160 continue
      if(np.ge.3) go to 300
      step = -step
      enrgy = ebis(1)
      cx = cbis(1)
      ebis(1) = ebis(2)
      cbis(1) = cbis(2)
      ebis(2) = enrgy
      cbis(2) = cx
      go to 100
c
c ---- minimum energy on the line is bracketed.
c      fit parabola to last three points 
c
  300 continue
      cx1 = cbis(np-2)
      f1  = ebis(np-2)
      cx2 = cbis(np-1)
      f2  = ebis(np-1)
      cx3 = cbis(np  )
      f3  = ebis(np  )
c
      c2 = ((f1-f2)*(cx2-cx3) - (f2-f3)*(cx1-cx2))
      c2 = c2/((cx1*cx1-cx2*cx2)*(cx2-cx3) - 
     *     (cx2*cx2-cx3*cx3)*(cx1-cx2))
      c1 = (f1 - f2 + c2 * (cx2*cx2-cx1*cx1))/(cx1-cx2)
      c0 = f1 - c2*cx1*cx1 - c1*cx1
c
c ---- predict extremum of the parabola -----
c
  320 continue
      write(iw,9100) c0,c1,c2
      cxmin = -c1/(c2+c2)
      emin = c2 * cxmin*cxmin + c1 * cxmin + c0
      write(iw,9120) cxmin,emin
c
c ----- generate coordinates of point qc' (in qc) at the
c       predicted minimum of the linear search -----
c
      do 400 i=1,ncoord
         qc(i) = qc(i) + cxmin * bis(i)
  400 continue
      return
c
 9000 format(//10x,'**** error, linear search takes too',
     *          ' many points to find minimum')
 9020 format(//10x,'point',i3,' on linear search...',
     *          'step=',f9.5,' bohr')
 9040 format(//10x,'atom',16x,'x',14x,'y',14x,'z'/
     *        10x,55(1h-)/)
 9060 format(10x,a8,2x,3f15.6)
 9080 format(//10x,'**** error, scf failed to converge')
 9100 format(//10x,'the parabolic fit is -'/
     *       10x,f16.10,' + ',f16.10,' * step +',f16.10,
     *           ' * step**2 = e(step)')
 9120 format(//10x,'the minimum of the fit is at step=',
     *        f12.6,' with a predicted energy of',f16.9//
     *    10x,'stabilized point qc'' on the irc has been found'/)
      end
      subroutine stbimk(core,iw)
      implicit real*8 (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
c
      parameter(szero=0.0d+00,sone=1.0d+00, two=2.0d+00, three=3.0d+00
     *         ,sone80=180.0d+00, gtol=1.0d-08, tol5=5.0d-06,
     *          nnam=23,mxa3=maxat*3)
c
c ...
c ... timing statistics
c ...
      real*8 begin, ebegin, timsec, walsec, tstart, estart
      common/statis/begin,ebegin,timsec(50),walsec(50),
     + tstart,estart
c
c
      integer nmaxly, ntotly, icurly, isecly
      logical oprintm
      common /segm/ nmaxly,ntotly,icurly,isecly(100),
     +              oprintm
c
      common/junkc /ylab(26),ztype(831),zoptio(13),ztagg(100)
      common/data1 /vlist(400),np(206),np1(maxorb,2),norbt(6),
     *              norbta(maxorb,2),norb2,nseco(maxorb),norb3,
     *              nthird(maxorb),norbg,norbga(maxorb),nscrp,
     *              range(2),alimit,blimit,isymm,nnx,lsym,nsymg(8),
     *              iuse(8,6),idorb(100),ncont,nmax(10)
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
c from March 97, amass is used only for user-input
c mass values
c czanr contains a copy of czan made in basis.m after reorganisation and
c before pseudo potential; so real charges in real order 
c beware, infob is explicit in direct.m and  dirrpa.m
c
      real*8 czann, czanr, cat, amass, cnew
      integer nonsym, map80
      logical ozmat
      common/infob/czanr(maxat),czann(maxat),cat(3,maxat),amass(maxat),
     +             cnew(maxat3,3),nonsym,map80(maxat),ozmat

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
      common/miscop/paa(mxa3),grd(mxa3)
c
      common /irc1 /ea,eb,ec,sab,sbc,stotal,ganorm,gbnorm,gcnorm,
     *              qa(mxa3),qb(mxa3),qc(mxa3),ga(mxa3),gb(mxa3),
     *              gc(mxa3),qa2(mxa3),qa1(mxa3),qa0(mxa3),ga2(mxa3),
     *              ga1(mxa3),ga0(mxa3),qcp(mxa3),qcc(mxa3),gcp(mxa3),
     *              bis(mxa3),ebis(10),cbis(10),gbis,rtrms(mxa3),
     *              ck0(mxa3),ck1(mxa3),ck2(mxa3),ck3(mxa3),
     *              ipirc,npirc
      common /irc2 /oforw,osadd,ostab,otsen,ostbs,olimn,oline,oquad,
     *              oampc,orung,ovtst,ndummy,evib,delta,elbow,stride,
     *              npoint,nextpt,freq,ypac,yyy,ncoord,nhopt,nampc,
     *              ofirst,nn
c
c
      real*8 enrgy, egrad
      common/gms_funct/enrgy,egrad(3*maxat)
c
c
c ----- stabilize the reaction path by setting up the
c       search along the bisector as described by
c       ishida, morokuma, and komornicki ----
c
      write(iw,9000)
      i0 = 0
      do 100 i=1,nat
         write(iw,9020) zaname(i),(qc(i0+j),j=1,3)
         i0 = i0 + 3
  100 continue
c
c ---- energy and gradient at point qc -----
c      qc is the point not quite on the path, which the bisector
c      will pass through.
c
      ii=0
      do 120 i=1,nat
         c(1,i) = qc(ii+1)
         c(2,i) = qc(ii+2)
         c(3,i) = qc(ii+3)
         ii=ii+3
  120 continue
      enrgy = szero
      call optx(core)
      enrgy=etot
      if (enrgy .eq. szero) call caserr('zero energy evaluated')
      ec = enrgy
      do 400 i=1,ncoord
         egrad(i) = grd(i)
  400 continue
c
c     ----- normalized, mass weighted gradient at point qc -----
c
      gcnorm = szero
      do 420 i=1,ncoord
         g = rtrms(i) * egrad(i)
         if(dabs(g).lt.gtol) g = szero
         gcnorm = gcnorm + g * g
         gc(i) = g
  420 continue
      gcnorm = dsqrt(gcnorm)
      write(iw,9060) gcnorm
      fact = sone/gcnorm
      do 440 i=1,ncoord
         gc(i) = fact * gc(i)
  440 continue
      write(iw,9080)
      i0 = 0
      do 460 i=1,nat
         write(iw,9100) zaname(i),(gc(i0+j),j=1,3)
         i0 = i0 + 3
  460 continue
c
c ---- test gradient at point qc and qb for collinearity -----
c      if these are sufficiently parallel, the stabilization is skipped.
c      note that gb and gc are already normalized.
c
      dotprd = szero
      do 500 i=1,ncoord
         dotprd = dotprd + gc(i)*gb(i)
  500 continue
      theta = dacos(dotprd)
      pi = dacos(-sone)
      theta = sone80*theta/pi
      theta = sone80 - theta
      write(iw,9120) theta
      if(theta.lt.elbow) go to 600
c
c     ----- stabilization is not necessary.  print and exit -----
c
      next = npirc + 1
      write(iw,9140) elbow,next
      ostbs = .true.
c     if(nzmat.gt.0) then
c           call bandbi
c           call trang(egrad,nvar,ncoord)
c           call pzandg(egrad,1)
c        end if
      return
c
c ----- determine unit vector along bisector
c       of the unit vectors along -gb and -gc -----
c
  600 continue
      write(iw,9180)
      ostbs = .false.
      do 610 i=1,ncoord
         bis(i) = gb(i) - gc(i)
  610 continue
      call symdr(bis,core)
      bisnrm = szero
      do 620 i=1,ncoord
         bisnrm = bisnrm + bis(i) * bis(i)
  620 continue
      bisnrm = dsqrt(bisnrm)
      if(bisnrm.lt.gtol)
     *   call caserr('norm of gradient difference less than gtol')
      fact = sone/bisnrm
      do 640 i=1,ncoord
         bis(i) = fact * bis(i)
  640 continue
c
c ---- get component of gradient along bisector direction -----
c      eg still holds un-mass-weighted, un-normalized gradient at qc.
c
      gbis = szero
      do 700 i=1,ncoord
         gbis = gbis + egrad(i) * bis(i)
  700 continue
      write(iw,9200) gbis
      write(iw,9220)
      i0 = 0
      do 720 i=1,nat
         write(iw,9100) zaname(i),(bis(i0+j),j=1,3)
         i0 = i0 + 3
  720 continue
c
c     ----- search along the bisector for the energy minimum qc' -----
c
      call seabis(core,iw)
      return
c
 9000 format(//10x,'the coordinates of point qc (bohrs) are'/
     *        10x,'atom',16x,'x',14x,'y',14x,'z'/
     *        10x,59(1h-)/)
 9020 format(10x,a8,2x,4x,3f15.10)
c9040 format(//10x,'**** error, scf failed to converge')
 9060 format(//10x,'the normalized mass weighted gradient at point qc'/
     *        10x,'(original gcnorm=',f10.6,' hartree/bohr-sqrt(amu))')
 9080 format(/10x,'atom',12x,'de/dx',7x,'de/dy',7x,'de/dz'/
     *        10x,46(1h-))
 9100 format(10x,a8,2x,3(2x,f10.6))
 9120 format(//10x,'the angle between gradients gb and gc is',f8.2,
     *         ' degrees')
 9140 format(//10x,'**** stabilization is being skipped * * *'/
     *        10x,'      (point qc is taken to be qc''.)'/
     *        10x,'because the angle between gradients exceeds',
     *           ' elbow angle',f8.2/
     *        10x,'the preceeding point is therefore point',i4,
     *           ' on the reaction path.'/
     *        10x,'when this point is relabeled from qc'' to qb,'/
     *        10x,'its energy and gradient will be already known,'/
     *        10x,'so that we can proceed immediately ',
     *           'to look for the next irc point.')
 9180 format(///10x,'*** ishida, morokuma, komornicki stabilization',
     *              ' of the reaction path ***'///
     *       10x,'the point qc above is at step=0.0 bohr'/
     *       10x,'looking for minimum energy point qc'' along the',
     *           ' gradient bisector')
 9200 format(//10x,'the component of the gradient at qc along',
     *     ' the bisector direction is',f11.6,' hartree/bohr')
 9220 format(//10x,'the unit vector along the bisector',
     *            ' direction is (dimensionless)'//
     *            10x,'atom',12x,'x',11x,'y',11x,'z'/
     *            10x,50(1h-))
      end
      subroutine ver_optef(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/optef.m,v $
     +     "/
      data revision /"$Revision: 6176 $"/
      data date /"$Date: 2010-08-10 16:49:47 +0200 (Tue, 10 Aug 2010) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
