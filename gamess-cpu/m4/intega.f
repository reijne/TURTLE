      subroutine chek1e(ochek)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
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
      common/blkin/dxyz(4),potold(4),nato,numo
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
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
      data m2/2/,thresh/1.0d-6/
c
      ochek = .false.
      dxyz(1) = enucf(nat,czan,c)
      do 20 i = 2 , 4
         dxyz(i) = 0.0d0
 20   continue
      do 30 i = 1 , nat
         dxyz(2) = dxyz(2) + czan(i)*c(1,i)
         dxyz(3) = dxyz(3) + czan(i)*c(2,i)
         dxyz(4) = dxyz(4) + czan(i)*c(3,i)
 30   continue
      call sectst(ionsec,itest)
      if (itest.ne.0) then
         call secget(ionsec,m2,ibl192)
         call rdedx(potold,511,ibl192,idaf)
         do 40 i = 1 , 4
            if (dabs(dxyz(i)-potold(i)).gt.thresh) go to 50
 40      continue
         if (nato.eq.nat .and. numo.eq.num) then
            ochek = .true.
         end if
      end if
 50   return
      end
      subroutine chkout(ii,jj,kk,ll,fock,core)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
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
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
      logical osortp,oschw
      common/sortpdat/
     + igbf,iklbf,iiptbf,ijbase,ngbf,lengbf,osortp,oschw
c
c     dlnmxd: the maximum value of the logarithm of the reduced
c             density matrix.
c     dlnmxs: the maximum value of the logarithm of the Schwarz 
c             integrals.
c     dlntol: minus the logarithm of the direct SCF integral cutoff.
c
      real*8 dlnmxd, dlnmxs, dlntol, tolitr, deltol, delfac
      integer intcut, intmag
      integer ibl171, lentri, itrtol, iexch, iof171, m171t
      integer nhamd, nsheld
      logical odscf, oimag, odnew, opref, odelta
      common/cslosc/dlnmxd,dlnmxs,dlntol,tolitr(3),deltol,delfac,
     + odscf,
     + intcut(10),oimag,intmag(1060),
     + ibl171,lentri,itrtol(3),iexch,odnew,opref,iof171(6),m171t,
     + nhamd, nsheld,odelta
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
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
c
      logical ospbas, onocnt, ochek, opdbas, opfbas, opgbas
      integer kad, istd
      common /ijlab/ kad(mxshel),ospbas,onocnt,istd,ochek,
     +               opdbas,opfbas,opgbas
c
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer n2file, n2tape, n2blk, n2last
      integer n4file, n4tape, n4blk, n4last
      integer n6file, n6tape, n6blk, n6last
      integer n5file, n5tape, n5blk, n5last
      integer n9file, n9tape, n9blk, n9last
      integer ntfile, nttape, ntblk, ntlast
      integer n1file, n1tape, n1blk, n1last
      integer n11fil, n11tap, n11bl, n11lst
      integer n12fil, n12tap, n12bl, n12lst
      integer n13fil, n13tap, n13bl, n13lst
      integer maxbl, minbl, maxset
      common/gms_files/
     +      n2file,n2tape(20),n2blk(20),n2last(20),
     +      n4file,n4tape(20),n4blk(20),n4last(20),
     +      n6file,n6tape(20),n6blk(20),n6last(20),
     +      n5file,n5tape(20),n5blk(20),n5last(20),
     +      n9file,n9tape(20),n9blk(20),n9last(20),
     +      ntfile,nttape(20),ntblk(20),ntlast(20),
     +      n1file,n1tape(20),n1blk(20),n1last(20),
     +      n11fil,n11tap(20),n11bl(20),n11lst(20),
     +      n12fil,n12tap(20),n12bl(20),n12lst(20),
     +      n13fil,n13tap(20),n13bl(20),n13lst(20),
     +       maxbl(maxlfn),minbl(maxlfn),maxset(maxlfn)
c
c
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
c
      logical oianj, okanl, oident, omisc, oham, opdipd, omp2
      integer ipos1, ipos2
      common /misc/ oianj,okanl,oident,omisc,
     +            oham,opdipd,omp2,ipos1,ipos2
c
      common/scfopt/maxit(4),accdi(2),odiis(4),dmpcut(7),iter
c ***
      dimension core(*),fock(*)
c ***
      if (omaxb) then
         call maxout
         return
      end if
      if (ochek) return
      ochek = .true.
      t4 = cpulft(0) - safety
      if ((timlst-t4)*2.5d0.lt.t4) then
         timlst = t4
      else
         osafe = safety.le.safe
         irest = 1
         ista = ii
         jsta = jj
         ksta = kk
         lsta = ll + 1
         if (nopk.eq.1) then
            max = ksta
            if (ista.eq.ksta) max = jsta
            if (lsta.le.max) go to 20
            lsta = 1
            ksta = kk + 1
            if (ksta.le.ii) go to 20
         else
            if (lsta.le.kk) go to 20
            lsta = 1
            ksta = kk + 1
            if (ksta.le.jj) go to 20
         end if
         ksta = 1
         jsta = jj + 1
         if (jsta.gt.ii) then
            jsta = 1
            ista = ii + 1
            if (ista.gt.nshell) return
         end if
 20      if (odscf) then
            if (icount.gt.1) then
               ochek = .false.
               nrec = nrec + 1
            end if
         else
            if (osortp) then
               call zsortp(core(igbf),core(iklbf),core(iiptbf),ijbase)
               do 30 ibuf = 1 , ngbf
                  call clrgbf(ibuf,core(igbf),core(iklbf),core(iiptbf))
 30            continue
            end if
            call blocki
         end if
c
         isti = ista
         jsti = jsta
         ksti = ksta
         lsti = lsta
         lastb = iblkmp
         lastu = mfilep
         m2file = mfilep
         m2last(m2file) = iblkmp
         ist = ista
         jst = jsta
         kst = ksta
         lst = lsta
         ss = cpulft(1)
         if (osafe .or. outv) write (iwr,6030) ss , ist , jst , kst ,
     +                               lst
         if (odscf) then
            call wrt3(fock,lentri,ibl171+iof171(3),idaf)
         else if (osafe .or. outv) then
            write (iwr,6040)
            call filprn(m2file,m2blk,m2last,m2tape)
         end if
         itask(mtask) = irest
         call revise
         if (omaxb) then
            call maxout
            return
         end if
         t4 = dumtim*1.5d0
         if (osafe) then
            write (iwr,6010)
            write (iwr,6050)
            tim = timlim + 0.5d0
         else
            safety = cpulft(0)
            if (safety.gt.t4) then
               safety = safety - dumtim
            else
               safety = safe
            end if
            timlst = cpulft(0) - safety
            if (ista.lt.nshell) then
               if (outv) then
                  if (odscf) then
                     if (iter.le.0) write (iwr,6020)
                  else
                     write (iwr,6020)
                  end if
               end if
            end if
         end if
         call clredx
      end if
      return
c
 6010 format (/' *** insufficient time to complete integral evaluation'/
     +        )
 6020 format (//1x,58('-')/1x,'ist',2x,'jst',2x,'kst',2x,'lst',7x,
     +        'nrec',3x,'intloc',5x,'del(t)',5x,'time'/1x,58('-')/)
 6030 format (/1x,58('-')//' job dumped at ',f8.2,
     +        ' seconds'//' next batch ',4i5/)
 6040 format (/' status of mainfile'/1x,18('-'))
 6050 format (/10x,27('*')/10x,'*** warning ***'/10x,
     +        'this job must be restarted'/10x,27('*')/)
c
      end
      subroutine denfan(dij,csi,cpi,cdi,cfi,cgi,
     +                      csj,cpj,cdj,cfj,cgj,
     +    mini,maxi,minj,maxj,oianj,double,norm)
      implicit real*8  (a-h,o-z)
      logical oianj,double,norm
c
      real*8 pito52, pidiv4, root3, root5, root53, root7
      common /picon/ pito52,pidiv4,root3,root5,root53,root7
c
      dimension dij(*)
      max = maxj
      nn = 0
      do 170 i = mini , maxi
         go to (20,30,80,80,
     +          40,80,80,60,80,80,
     +          50,80,80,70,80,80,80,80,80,60,
     +          82,80,80,84,80,80,80,80,80,86,
     +          80,80,88,80,80),i
 20      dum1 = csi
         go to 80
 30      dum1 = cpi
         go to 80
 40      dum1 = cdi
         go to 80
 50      dum1 = cfi
         go to 80
 60      if (norm) dum1 = dum1*root3
         go to 80
 70      if (norm) dum1 = dum1*root5
         go to 80
 82      dum1 = cgi
         go to 80
 84      if (norm) dum1 = dum1*root7
         go to 80
 86      if (norm) dum1 = dum1*root5/root3
         go to 80
 88      if (norm) dum1 = dum1*root3
 80      if (oianj) max = i
         do 160 j = minj , max
            go to (90,100,150,150,
     +            110,150,150,130,150,150,
     +            120,150,150,140,150,150,150,150,150,130,
     +            142,150,150,144,150,150,150,150,150,146,
     +            150,150,148,150,150),j
 90         dum2 = dum1*csj
            if (double) then
               if (i.gt.1) then
                  dum2 = dum2 + csi*cpj
               else
                  dum2 = dum2 + dum2
               end if
            end if
            go to 150
 100        dum2 = dum1*cpj
            if (double) dum2 = dum2 + dum2
            go to 150
 110        dum2 = dum1*cdj
            if (double) dum2 = dum2 + dum2
            go to 150
 120        dum2 = dum1*cfj
            if (double) dum2 = dum2 + dum2
            go to 150
 130        if (norm) dum2 = dum2*root3
            go to 150
 140        if (norm) dum2 = dum2*root5
            go to 150
 142        dum2 = dum1*cgj
            if (double) dum2 = dum2+dum2
            go to 150
 144        if (norm) dum2 = dum2*root7
            go to 150
 146        if (norm) dum2 = dum2*root5/root3
            go to 150
 148        if (norm) dum2 = dum2*root3
 150        nn = nn + 1
            dij(nn) = dum2
 160     continue
 170  continue
      return
      end
      subroutine dipint
c
c     ----- gauss-hermite quadrature using minimum point formula -----
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      common/junk/ssgg(450),
     *pint0,qint0,rint0,pintx,qinty,rintz,t,p0,q0,r0,
     + pi,qi,ri,pj,qj,rj,ni,nj
      common/hermit/h(45)
      common/wermit/w(45)
cdrf drfneeds possibility of other centre than zero,zero,zero
c    so this common ONLY carries the cmx,cmy and cmz values
c    which should be zero on GAMESS-start! (added in mainsrf.f)
      common/c_of_m/pcm,qcm,rcm
      dimension min(6),max(6)
      data min /1,2,4,7,11,16/
      data max /1,3,6,10,15,21/
      data dzero /0.0d0/
      pint0 = dzero
      qint0 = dzero
      rint0 = dzero
      pintx = dzero
      qinty = dzero
      rintz = dzero
      npts = (ni+nj-2+1)/2 + 1
      imin = min(npts)
      imax = max(npts)
      do 140 i = imin , imax
         dum = w(i)
         px = dum
         py = dum
         pz = dum
         dum = h(i)/t
         ptx = dum + p0
         pty = dum + q0
         ptz = dum + r0
         ax = ptx - pi
         ay = pty - qi
         az = ptz - ri
         bx = ptx - pj
         by = pty - qj
         bz = ptz - rj
         go to (60,50,40,30,20) , ni
 20      px = px*ax
         py = py*ay
         pz = pz*az
 30      px = px*ax
         py = py*ay
         pz = pz*az
 40      px = px*ax
         py = py*ay
         pz = pz*az
 50      px = px*ax
         py = py*ay
         pz = pz*az
 60      go to (130,120,110,100,90,80,70) , nj
 70      px = px*bx
         py = py*by
         pz = pz*bz
 80      px = px*bx
         py = py*by
         pz = pz*bz
 90      px = px*bx
         py = py*by
         pz = pz*bz
 100     px = px*bx
         py = py*by
         pz = pz*bz
 110     px = px*bx
         py = py*by
         pz = pz*bz
 120     px = px*bx
         py = py*by
         pz = pz*bz
 130     pint0 = pint0 + px
         qint0 = qint0 + py
         rint0 = rint0 + pz
cdrf note the addiditon of 'pcm,qcm,rcm' for c_of_mass position!
         pintx = pintx + px*(ptx-pcm)
         qinty = qinty + py*(pty-qcm)
         rintz = rintz + pz*(ptz-rcm)
 140  continue
      return
      end
      subroutine dipxyz(ps,qs,rs)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
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
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
      common/junk/ssgg(450),
     + pint0,qint0,rint0,pintx,qinty,rintz,t,p0,q0,r0,
     + pi,qi,ri,pj,qj,rj,ni,nj
      common/blkin/pxyz(4),sx(225),sy(225),sz(225),
     + dij(225),pin(125),qin(125),rin(125),
     + ijx(225),ijy(225),ijz(225)
      dimension ps(*),qs(*),rs(*)
      dimension ix(35),iy(35),iz(35),jx(35),jy(35),jz(35)
      data dzero /0.0d0/, done/1.0d0/
      data sqrt3 /1.73205080756888d0/
      data sqrt5 /2.23606797749979d0/
      data sqrt7 /2.64575131106459d0/
      data rln10 /2.30258d0/
c
      data jx / 0, 1, 0, 0, 2, 0, 0, 1, 1, 0,
     +          3, 0, 0, 2, 2, 1, 0, 1, 0, 1,
     +          4, 0, 0, 3, 3, 1, 0, 1, 0, 2,
     +          2, 0, 2, 1, 1/
      data ix / 1, 6, 1, 1,11, 1, 1, 6, 6, 1,
     +         16, 1, 1,11,11, 6, 1, 6, 1, 6,
     +         21, 1, 1,16,16, 6, 1, 6, 1,11,
     +         11, 1,11, 6, 6/
      data jy / 0, 0, 1, 0, 0, 2, 0, 1, 0, 1,
     +          0, 3, 0, 1, 0, 2, 2, 0, 1, 1,
     +          0, 4, 0, 1, 0, 3, 3, 0, 1, 2,
     +          0, 2, 1, 2, 1/
      data iy / 1, 1, 6, 1, 1,11, 1, 6, 1, 6,
     +          1,16, 1, 6, 1,11,11, 1, 6, 6,
     +          1,21, 1, 6, 1,16,16, 1, 6,11,
     +          1,11, 6,11, 6/
      data jz / 0, 0, 0, 1, 0, 0, 2, 0, 1, 1,
     +          0, 0, 3, 0, 1, 0, 1, 2, 2, 1,
     +          0, 0, 4, 0, 1, 0, 1, 3, 3, 0,
     +          2, 2, 1, 1, 2/
      data iz / 1, 1, 1, 6, 1, 1,11, 1, 6, 6,
     +          1, 1,16, 1, 6, 1, 6,11,11, 6,
     +          1, 1,21, 1, 6, 1, 6,16,16, 1,
     +         11,11, 6, 6,11/
c
      tol = rln10*itol
      onorm = normf.ne.1 .or. normp.ne.1
      do 310 ii = 1 , nshell
c
c     ----- ishell
c
         i = katom(ii)
         pi = c(1,i)
         qi = c(2,i)
         ri = c(3,i)
         i1 = kstart(ii)
         i2 = i1 + kng(ii) - 1
         lit = ktype(ii)
         mini = kmin(ii)
         maxi = kmax(ii)
         loci = kloc(ii) - mini
c
c     ----- jshell
c
         do 300 jj = 1 , ii
            j = katom(jj)
            pj = c(1,j)
            qj = c(2,j)
            rj = c(3,j)
            j1 = kstart(jj)
            j2 = j1 + kng(jj) - 1
            ljt = ktype(jj)
            minj = kmin(jj)
            maxj = kmax(jj)
            locj = kloc(jj) - minj
            rr = (pi-pj)**2 + (qi-qj)**2 + (ri-rj)**2
            oiandj = ii.eq.jj
c
c     ----- prepare indices for pairs of (i,j) functions
c
            ij = 0
            max = maxj
            do 30 i = mini , maxi
               nnx = ix(i)
               nny = iy(i)
               nnz = iz(i)
               if (oiandj) max = i
               do 20 j = minj , max
                  ij = ij + 1
                  ijx(ij) = nnx + jx(j)
                  ijy(ij) = nny + jy(j)
                  ijz(ij) = nnz + jz(j)
 20            continue
 30         continue
            do 40 i = 1 , ij
               sx(i) = dzero
               sy(i) = dzero
               sz(i) = dzero
 40         continue
c
c     ----- i primitive
c
            jgmax = j2
            do 270 ig = i1 , i2
               ai = ex(ig)
               arri = ai*rr
               axi = ai*pi
               ayi = ai*qi
               azi = ai*ri
               csi = cs(ig)
               cpi = cp(ig)
               cdi = cd(ig)
               cfi = cf(ig)
               cgi = cg(ig)
c
c     ----- j primtive
c
               if (oiandj) jgmax = ig
               do 260 jg = j1 , jgmax
                  aj = ex(jg)
                  aa = ai + aj
                  dum = aj*arri/aa
                  if (dum.le.tol) then
                     fac = dexp(-dum)
                     csj = cs(jg)
                     cpj = cp(jg)
                     cdj = cd(jg)
                     cfj = cf(jg)
                     cgj = cg(jg)
                     ax = (axi+aj*pj)/aa
                     ay = (ayi+aj*qj)/aa
                     az = (azi+aj*rj)/aa
c
c     ----- density factor
c
                     odoub = oiandj .and. ig.ne.jg
                     max = maxj
                     nn = 0
                     do 220 i = mini , maxi
                        go to (50,60,120,120,
     +                      70,120,120,80,120,120,
     +                      90,120,120,100,120,120,120,120,120,110,
     +                      112,120,120,114,120,120,120,120,120,116,
     +                      120,120,118,120,120),i
 50                     dum1 = csi*fac
                        go to 120
 60                     dum1 = cpi*fac
                        go to 120
 70                     dum1 = cdi*fac
                        go to 120
 80                     if (onorm) dum1 = dum1*sqrt3
                        go to 120
 90                     dum1 = cfi*fac
                        go to 120
 100                    if (onorm) dum1 = dum1*sqrt5
                        go to 120
 110                    if (onorm) dum1 = dum1*sqrt3
                        go to 120
 112                    dum1 = cgi*fac
                        go to 120
 114                    if (onorm) dum1 = dum1*sqrt7
                        go to 120
 116                    if (onorm) dum1 = dum1*sqrt5/sqrt3
                        go to 120
 118                    if (onorm) dum1 = dum1*sqrt3
 120                    if (oiandj) max = i
                        do 210 j = minj , max
                           go to (130,140,200,200,
     +                      150,200,200,160,200,200,
     +                      170,200,200,180,200,200,200,200,200,190,
     +                      192,200,200,194,200,200,200,200,200,196,
     +                      200,200,198,200,200),j
 130                       dum2 = dum1*csj
                           if (odoub) then
                              if (i.gt.1) then
                                 dum2 = dum2 + csi*cpj*fac
                              else
                                 dum2 = dum2 + dum2
                              end if
                           end if
                           go to 200
 140                       dum2 = dum1*cpj
                           if (odoub) dum2 = dum2 + dum2
                           go to 200
 150                       dum2 = dum1*cdj
                           if (odoub) dum2 = dum2 + dum2
                           go to 200
 160                       if (onorm) dum2 = dum2*sqrt3
                           go to 200
 170                       dum2 = dum1*cfj
                           if (odoub) dum2 = dum2 + dum2
                           go to 200
 180                       if (onorm) dum2 = dum2*sqrt5
                           go to 200
 190                       if (onorm) dum2 = dum2*sqrt3
                           go to 200
 192                       dum2 = dum1*cgj
                           if (odoub) dum2 = dum2+dum2
                           go to 200
 194                       if (onorm) dum2 = dum2*sqrt7
                           go to 200
 196                       if (onorm) dum2 = dum2*sqrt5/sqrt3
                           go to 200
 198                       if (onorm) dum2 = dum2*sqrt3
 200                       nn = nn + 1
                           dij(nn) = dum2
 210                    continue
 220                 continue
c
c     ----- dipole moment integrals -----
c
                     t = dsqrt(aa)
                     tinv = done/t
                     p0 = ax
                     q0 = ay
                     r0 = az
                     in = -5
                     do 240 i = 1 , lit
                        in = in + 5
                        ni = i
                        do 230 j = 1 , ljt
                           jn = in + j
                           nj = j
                           call dipint
                           pin(jn) = pint0*tinv
                           qin(jn) = qint0*tinv
                           rin(jn) = rint0*tinv
                           pin(jn+25) = pintx*tinv
                           qin(jn+25) = qinty*tinv
                           rin(jn+25) = rintz*tinv
 230                    continue
 240                 continue
                     do 250 i = 1 , ij
                       nnx = ijx(i)
                       nny = ijy(i)
                       nnz = ijz(i)
                       sx(i) = sx(i) + dij(i)
     +                         *(pin(nnx+25)*qin(nny)*rin(nnz))
                       sy(i) = sy(i) + dij(i)
     +                         *(pin(nnx)*qin(nny+25)*rin(nnz))
                       sz(i) = sz(i) + dij(i)
     +                         *(pin(nnx)*qin(nny)*rin(nnz+25))
 250                 continue
                  end if
 260           continue
 270        continue
c
c     ----- set up dipole moment matrices -----
c
            max = maxj
            nn = 0
            do 290 i = mini , maxi
               li = loci + i
               in = (li*(li-1))/2
               if (oiandj) max = i
               do 280 j = minj , max
                  lj = locj + j
                  jn = lj + in
                  nn = nn + 1
                  ps(jn) = sx(nn)
                  qs(jn) = sy(nn)
                  rs(jn) = sz(nn)
 280           continue
 290        continue
 300     continue
 310  continue
      return
      end
      subroutine ffun(x,f,npt,dji,madd)
      implicit real*8  (a-h,o-z)
      dimension f(60),madd(20),dji(60)
      dimension xmax(20)
      data xmax/24.0d0,29.0d0,32.0d0,35.0d0,37.0d0,40.0d0,42.0d0,
     &   44.0d0,47.0d0,49.0d0,51.0d0,53.0d0,55.0d0,57.0d0,58.0d0,
     &   59.0d0,59.0d0,59.0d0,59.0d0,59.0d0/
      m = npt + npt - 1
      fact = 4.8d0 + 0.4d0*npt
      if (x.ge.fact) then
         x2 = 0.5d0/x
         xinv = x2 + x2
         f1 = dsqrt(1.5707963267949d0*x2)
         if (x.ge.xmax(m)) then
c...... very high argument
            f(1) = f1
            if (m.gt.0) then
               do 20 i = 1 , m
                  f(i+1) = f(i)*x2
                  x2 = x2 + xinv
 20            continue
               return
            end if
         else
c.....  high argument
            x3 = dexp(-x)
            if (x.gt.21.6d0) then
               f1 = f1 - x3*x2
            else if (x.gt.18.2d0) then
               f1 = f1 - x3*x2*(1.0d0-x2)
            else if (x.gt.(12.0d0+0.1d0*npt)) then
               f1 = ((1.9623264149430d-1*xinv-4.9695241464490d-1)
     +              *xinv-6.0156581186481d-5)*x3 + f1
            else if (x.gt.(9.2d0+0.2d0*npt)) then
               f1 = (((-1.8784686463512d-1*xinv+2.2991849164985d-1)
     +              *xinv-4.9893752514047d-1)*xinv-2.1916512131607d-5)
     +              *x3 + f1
            else
               f1 = f1 +
     +              ((((((4.6897511375022d-1*xinv-6.9955602298985d-1)
     +              *xinv+5.3689283271887d-1)*xinv-3.2883030418398d-1)
     +              *xinv+2.4645596956002d-1)*xinv-4.9984072848436d-1)
     +              *xinv-3.1501078774085d-6)*x3
            end if
            f(1) = f1
            if (m.gt.0) then
               x23 = x2*x3
               do 30 i = 1 , m
                  f(i+1) = f(i)*x2 - x23
                  x2 = x2 + xinv
 30            continue
            end if
         end if
c.....  low argument
      else if (x.lt.1.0d-10) then
         j = m + 1
         do 40 i = 1 , j
            f(i) = dji(i)
 40      continue
         return
      else if (x.lt.1.0d-5) then
         j = m + 1
         do 50 i = 1 , j
            f(i) = dji(i) - x*dji(i+1)
 50      continue
         return
      else
         x2 = x + x
         i = x2 + madd(m)
         j = i
         x3 = dexp(-x)
         f(j+1) = 0.0d0
         do 60 k = 1 , i
            f(j) = (f(j+1)*x2+x3)*dji(j)
            j = j - 1
 60      continue
         return
      end if
      return
      end
      subroutine final(core,fock,dens)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
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
c
c     dlnmxd: the maximum value of the logarithm of the reduced
c             density matrix.
c     dlnmxs: the maximum value of the logarithm of the Schwarz 
c             integrals.
c     dlntol: minus the logarithm of the direct SCF integral cutoff.
c
      real*8 dlnmxd, dlnmxs, dlntol, tolitr, deltol, delfac
      integer intcut, intmag
      integer ibl171, lentri, itrtol, iexch, iof171, m171t
      integer nhamd, nsheld
      logical odscf, oimag, odnew, opref, odelta
      common/cslosc/dlnmxd,dlnmxs,dlntol,tolitr(3),deltol,delfac,
     + odscf,
     + intcut(10),oimag,intmag(1060),
     + ibl171,lentri,itrtol(3),iexch,odnew,opref,iof171(6),m171t,
     + nhamd, nsheld,odelta
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
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
      common/blkin/gout(510),nword
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
c
      logical ospbas, onocnt, ochek, opdbas, opfbas, opgbas
      integer kad, istd
      common /ijlab/ kad(mxshel),ospbas,onocnt,istd,ochek,
     +               opdbas,opfbas,opgbas
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
      integer n2file, n2tape, n2blk, n2last
      integer n4file, n4tape, n4blk, n4last
      integer n6file, n6tape, n6blk, n6last
      integer n5file, n5tape, n5blk, n5last
      integer n9file, n9tape, n9blk, n9last
      integer ntfile, nttape, ntblk, ntlast
      integer n1file, n1tape, n1blk, n1last
      integer n11fil, n11tap, n11bl, n11lst
      integer n12fil, n12tap, n12bl, n12lst
      integer n13fil, n13tap, n13bl, n13lst
      integer maxbl, minbl, maxset
      common/gms_files/
     +      n2file,n2tape(20),n2blk(20),n2last(20),
     +      n4file,n4tape(20),n4blk(20),n4last(20),
     +      n6file,n6tape(20),n6blk(20),n6last(20),
     +      n5file,n5tape(20),n5blk(20),n5last(20),
     +      n9file,n9tape(20),n9blk(20),n9last(20),
     +      ntfile,nttape(20),ntblk(20),ntlast(20),
     +      n1file,n1tape(20),n1blk(20),n1last(20),
     +      n11fil,n11tap(20),n11bl(20),n11lst(20),
     +      n12fil,n12tap(20),n12bl(20),n12lst(20),
     +      n13fil,n13tap(20),n13bl(20),n13lst(20),
     +       maxbl(maxlfn),minbl(maxlfn),maxset(maxlfn)
c
c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
c
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
c
      logical oianj, okanl, oident, omisc, oham, opdipd, omp2
      integer ipos1, ipos2
      common /misc/ oianj,okanl,oident,omisc,
     +            oham,opdipd,omp2,ipos1,ipos2
c
c ***
      dimension core(*),fock(*),dens(*)
c ***
      data m0/0/
c
      if (odscf) then
         if (icount.gt.1) then
            ochek = .false.
            nrec = nrec + 1
         end if
      else
         call blocki
      end if
      if (omaxb) then
         ist = isti
         jst = jsti
         kst = ksti
         lst = lsti
         iblkmp = lastb
         irest = 1
         ss = cpulft(1)
         write (iwr,6010) ss , ist , jst , kst , lst
         m2last(m2file) = iblkmp
         itask(mtask) = irest
         call revise
         if (odscf) then
            call wrt3(fock,lentri,ibl171+iof171(3),idaf)
         else
            write (iwr,6020)
            call filprn(m2file,m2blk,m2last,m2tape)
         end if
         write (iwr,6050)
         tim = timlim + 0.5d0
         call clredx
c
      else
c
         m2last(m2file) = iblkmp
         if (.not.(odscf.or.imaxb_ic.eq.1)) then
            call search(iblkmp,mainp)
            if (onocnt) then
               call put(gout,m0,mainp)
               m2last(m2file) = m2last(m2file) + 1
               if (outv) then
                  write (iwr,6030)
                  call filprn(m2file,m2blk,m2last,m2tape)
               end if
               irest = 0
               ist = 1
               nrec = 1
            else
               ist = istd
            end if
         else if (onocnt) then
            m2last(m2file) = m2last(m2file) + 1
            if (outv) then
               write (iwr,6040) nrec
            end if
            irest = 0
            ist = 1
            nrec = 1
         else
            ist = istd
         end if
         jst = 1
         kst = 1
         lst = 1
         call timit(0)
         if (tim.lt.timlim) return
c
         if (onocnt) then
            irest = 2
         else
            irest = 1
            nindmx = 1
         end if
         itask(mtask) = irest
         call revise
         if (odscf) call wrt3(fock,lentri,ibl171+iof171(3),idaf)
         write (iwr,6050)
         tim = timlim + 0.5d0
         call clredx
      end if
      return
 6010 format (/1x,58('-')//' job dumped at ',f8.2,
     +        ' seconds'//' next batch ',4i5/)
 6020 format (/' status of mainfile'/1x,18('-'))
 6030 format (1x,58('-')/' status of mainfile'/1x,18('-'))
 6040 format (/1x,'no. of 2-electron integral records = ',i10)
 6050 format (/10x,27('*')//10x,'*** warning ***'/10x,
     +        'this job must be restarted'//10x,27('*')/)
      end
      subroutine indexa(ijx,ijy,ijz,ij,mini,maxi,
     &minj,maxj,iandj,inc1,inc2,inc3)
      implicit real*8  (a-h,o-z)
      logical iandj
      dimension ijx(225),ijy(225),ijz(225)
c
      integer ix, iy, iz
      common /inxblk/ ix(35),iy(35),iz(35)
c
      dimension jx(35),jy(35),jz(35)
c
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
      return
      end
      subroutine jandk(zscftp,core,f,fb,exch,d,db,prefac,rdmat)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
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
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
      common/scfopt/maxit(4),accdi(2),odiis(4),dmpcut(7),iter
c
c     dlnmxd: the maximum value of the logarithm of the reduced
c             density matrix.
c     dlnmxs: the maximum value of the logarithm of the Schwarz 
c             integrals.
c     dlntol: minus the logarithm of the direct SCF integral cutoff.
c
      real*8 dlnmxd, dlnmxs, dlntol, tolitr, deltol, delfac
      integer intcut, intmag
      integer ibl171, lentri, itrtol, iexch, iof171, m171t
      integer nhamd, nsheld
      logical odscf, oimag, odnew, opref, odelta
      common/cslosc/dlnmxd,dlnmxs,dlntol,tolitr(3),deltol,delfac,
     + odscf,
     + intcut(10),oimag,intmag(1060),
     + ibl171,lentri,itrtol(3),iexch,odnew,opref,iof171(6),m171t,
     + nhamd, nsheld,odelta
c
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
c ...
c ... timing statistics
c ...
      real*8 begin, ebegin, timsec, walsec, tstart, estart
      common/statis/begin,ebegin,timsec(50),walsec(50),
     + tstart,estart
c
c
      logical ospbas, onocnt, ochek, opdbas, opfbas, opgbas
      integer kad, istd
      common /ijlab/ kad(mxshel),ospbas,onocnt,istd,ochek,
     +               opdbas,opfbas,opgbas
c
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
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
      logical osortp,oschw
      common/sortpdat/
     + igbf,iklbf,iiptbf,ijbase,ngbf,lengbf,osortp,oschw
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
c
      logical omem, ofirt
      integer maxq, maxqq, iqqoff, mxblock
c
      common/maxlen/maxq,maxqq,omem(maxlfn),iqqoff(maxlfn),
     +              ofirt,mxblock
c     logical opk,opank,oskips
c     common/pkfil/opk,opank,oskips
      logical opk,opank
      common/pkfil/opk,opank
c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      integer icount_dlb, icount_slb
      common /par_dlb/ icount_dlb, icount_slb
c
      integer idpdiag, ipdiagmode, ntchnk, limchnk, iparapr,
     & idpdiis, idpmult2, idporth, ipiomode, iptest, ipinvmode,
     & ipdiagif
c
c     ipdiagmode - which diag (IDIAG_PEIGS/_PDSYEV/_PDSYEVD....)
c
c     ipdiagif   - which interface (IDIAG_NO_GAWRAP/IDIAG_GAWRAP)
c
c     idpdiag = dimension for parallel diag
c
c     ipinvmode : what matrix inversion algorithm to use
c                 either the old diag based approach or the more
c                 appropriate ScaLAPACK Cholesky factorisation can 
c                 be used (INV_CHOLESKY or INV_DIAG )
c
      logical odebugp

      common /parcntl/ idpdiag, ipdiagmode, ntchnk, limchnk, iparapr,
     & idpdiis, idpmult2, idporth, ipiomode, iptest, ipinvmode,
     & ipdiagif, odebugp
c
      logical ga_initted
      common/ gainit/ ga_initted

      integer IO_NZ,IO_NZ_S,IO_A
      parameter(IO_NZ=1,IO_NZ_S=2,IO_A=3)

      integer INV_CHOLESKY, INV_DIAG
      parameter(INV_CHOLESKY = 100)
      parameter(INV_DIAG     = INV_CHOLESKY + 1)
c
      integer IDIAG_PEIGS,   IDIAG_PDSYEV, IDIAG_PDSYEVX
      integer IDIAG_PDSYEVD, IDIAG_PDSYEVR
      parameter(IDIAG_PEIGS   = 10)
      parameter(IDIAG_PDSYEV  = IDIAG_PEIGS   + 1)
      parameter(IDIAG_PDSYEVX = IDIAG_PDSYEV  + 1)
      parameter(IDIAG_PDSYEVD = IDIAG_PDSYEVX + 1)
      parameter(IDIAG_PDSYEVR = IDIAG_PDSYEVD + 1)

      integer IDIAG_NO_GAWRAP, IDIAG_GAWRAP
      parameter(IDIAG_NO_GAWRAP=200)
      parameter(IDIAG_GAWRAP=201)

c
c     ZORA common
c
      logical ozora,oscalz,onlyp,onlyd,onlyf,onlyg,opzora,
     1        op2zora,onlys,ocontr,oint_zora,o1e_zora,osmall_zora,
     2        opre_zora,oso,onoext_z,oioraz,oatint_z,oatscf_z,
     3        oscalatz
      real*8 cspeed,critat_z
      integer nbas_zora,ibls_z,iblv_z,nwv_z,ibiso_z,nwiso_z,ibcor_z,
     1        ibsx_z,icoul_z,ibsy_z,ibsz_z,ibsmin_z,ibsplus_z,nwmat_z,
     2        ibtrin_z,ibtrout_z,ibscalf_z,nwtrat_z,ibscale_z,ibcoul_z,
     3        niter_z,int_zora,num_ext,igauge_z,nwcor_z,isexch_z,
     4        nat_z,ibdens_z,nwdens_z,ibscalao_z,ibshat_z,is_z,
     5        irest_z,numdu_z,ibldu_z,iblock_z,lshift_z,ibl7ew_z
      common/zorac/ cspeed,critat_z,ozora,oscalz,
     1              onlyp,onlyd,onlyf,onlyg,opzora,
     1              op2zora,onlys,ocontr,nbas_zora,
     2              oint_zora,o1e_zora,icoul_z,
     3              ibls_z,iblv_z,nwv_z,
     3              ibiso_z,nwiso_z,ibcor_z,nwcor_z,
     4              ibsx_z,ibsy_z,ibsz_z,ibsmin_z,ibsplus_z,nwmat_z,
     5              ibtrin_z,ibtrout_z,ibscalf_z,nwtrat_z,
     6              ibscale_z,ibcoul_z,ibdens_z,nwdens_z,osmall_zora,
     7              niter_z,opre_zora,int_zora,num_ext,isexch_z,
     8              oso,onoext_z,is_z,igauge_z,oioraz,
     9              nat_z,ibscalao_z,oatint_z,oatscf_z,ibshat_z,
     1              oscalatz,irest_z,numdu_z,ibldu_z,iblock_z,lshift_z,
     2              ibl7ew_z
c....    icoul_z : 1 : Full coulomb
c...               2 : atomic coulomb (default)
c...               3 : no coulomb (just 1-electron ints)
c...               4 : small basis coulomb
c...     isexch_z  0 : no exchange in zora correction (default)
c...               1 : exchange in zora correction (like in dft)
c...     niter_z : # iterations between coulomb matrix update
c...     int_zora :internal basis gen. : 0 none, 1, automatic, 2 by hand
c
c...     oscalz  : scaling on or off (true or false)
c...     oioraz  : change metric to perform inf. order scaling 
c...               (Dyall/EvLenthe (not quite))
c...     oscalatz: atoms only scaled zora (together with get atom)
c...     only'x' : internal basis only op to 'x'-type functions
c...     op(2)zora: print flags (op2zora+opzora: prints *everything*)
c...     ocontr  : adds l+1 block to internal basis in contracted way 
c...     nbas_zora: size of internal basis
c...     oint_zora: in internal or external mode???
c...     o1e_zora: one electron integrals have been calculated in *this*
c...               internal basis
c...     osmall_zora: small swap for projected coulomb
c...     opre_zora: pre_zora has been called
c...     oso     : spin orbit calculation 
c...     onoext_z: no external functions copied to internal basis
c...     is_z : section for orbitals used for  coulomb matrix
c...     igauge_z : add so many s-marices to h-matrix
c...                the energy should change by igauge
c...     nat_z  : indicates (if ne 0) that starting density is used
c...              to calc. coulomb; gives # it. to get atomic density
c...     critat_z : accuracy of atomic energy for atomic zora
c...     oatint_z : indicates that integrals to be calculated are atomic
c...     oatscf_z : indicates we are doing a set of atomic SCFs
c...     ibshat_z : block on num8 to store blocked s/h
c...     for disk-adress and lengths (ib.., nw.. ) see subroutine atoms2
c...     irest_z  : .ne.0 indicates atomic zora corrections from dumpfile; 1 not on restart
c...     numdu_z,ibldu_z : if ne 0 other dumpfile to read zora corrections from
c...     iblock_z is the block number on ed3, where the zora corrections (section 425) reside
c...     ibl7ew is where the energy weighted density matrix resides
      dimension core(*),f(*),fb(*),exch(*),d(*),db(*),prefac(*),rdmat(*)
      character*10 charwall
      character *8 fnm
      character *5 snm
      data fnm/'intega.m'/
      data snm/'jandk'/
      call cpuwal(begin,ebegin)
      ochek = .true.
c     oskips = .false.
c ***
      if (odscf) then
         call setsto(10,0,intcut)
         call setsto(1060,0,intmag)
         nopk = 1
      end if
      outvv = nprint.ne.-5.and.oprint(55)
c***
      if (nopk.ne.0) osortp = .false.
c     check for pure sp basis. if so, do gaussian integrals.
      call spchck
      call debut(zscftp)
      nopkr = nopk
      iofrst = iofsym
c
c *** allocate space for integral buffer (gout)
c *** integral driven    = 50625 (spdfg) ... nopk = 1
c *** supermatrix driven = 151875 (spdfg) ... nopk .ne. 1
c
c     nav = lenwrd()
      if (nopk.ne.1) then
         ilen = 151875
      else
         ilen =  50625
      endif
      i10 = igmem_alloc_inf(ilen,fnm,snm,'gout',IGMEM_DEBUG)
c
c *** allocate space for iso array
c
      if (oint_zora) then
         i20 = igmem_alloc_inf(nwiso_z,fnm,snm,'iso',IGMEM_DEBUG)
      else
         i20 = igmem_alloc_inf(nw196(5),fnm,snm,'iso',IGMEM_DEBUG)
      end if

c     schwarz inequality coulomb integral calculation ?
      if (oschw) then
c        oskips = .true.
         i30 = igmem_alloc_inf(ikyp(nshell),fnm,snm,'clints',
     +                         IGMEM_NORMAL)
         if (odscf) then
            if(irest.le.0. and. iter.le.0) 
     +      call coulmb(core(i30),core(i10))
         else
            call coulmb(core(i30),core(i10))
         end if
c        oskips = .false.
         call gmem_free_inf(i30,fnm,snm,'clints') 
      end if
c
      if (oint_zora) then
         call rdedx(core(i20),nwiso_z,ibiso_z,num8)
      else
         call rdedx(core(i20),nw196(5),ibl196(5),idaf)
      end if
c
      if (nopk.ne.1) then
c
c     supermatrix driven
c
         if (intg76.eq.0) go to 20
         if (nindmx.ne.0) go to 20
         call pkin70(core,core(i20),core(i10),nshell,outvv)
         if (omaxb .or. tim.ge.timlim) go to 30
         if (.not.(ospbas)) go to 20
         go to 30
c
c     integral driven
c
      else if (intg76.ne.0) then
         if (nindmx.eq.0) then
            if (oschw) then
             cpu1 = cpulft(1)
             call jkin7s(zscftp,core,f,fb,exch,d,db,prefac,rdmat,
     +                   core(i20),core(i10),nshell,outvv)
             if (outvv.and.odscf) then
              cpu2 = cpulft(1) - cpu1
              write (iwr,6046) cpu2
 6046         format (/' ** rotated axis integral evaluation required ',
     +        f8.2,' seconds')
             endif
            else
             call jkin7a(zscftp,core,f,fb,exch,d,db,prefac,rdmat,
     +                   core(i20),core(i10),nshell)
            end if

            if (omaxb .or. tim.ge.timlim) go to 30
            if (ospbas) go to 30
            ochek = .true.
         end if
      end if
      nindmx = 1
      onocnt = .true.

      if (oschw) then
        cpu1 = cpulft(1)
        call jkints(zscftp,core,f,fb,exch,d,db,prefac,rdmat,
     +              core(i20),core(i10),nshell,outvv)
        if (outvv.and.odscf) then
        cpu2 = cpulft(1) - cpu1
        write (iwr,6045) cpu2
 6045   format (/' ** gauss-rys integral evaluation required ',
     +             f8.2,' seconds')
       endif
      else
       call jkinta(zscftp,core,f,fb,exch,d,db,prefac,rdmat,
     +             core(i20),core(i10),nshell)
      end if
      go to 30
 20   nindmx = 1
      onocnt = .true.
      if (oschw) then
       call pkints(core,core(i20),core(i10),nshell,outvv)
      else
       call pkinta(core,core(i20),core(i10),nshell)
      end if
 30   cpu = cpulft(1)

      call gmem_free_inf(i20,fnm,snm,'iso')
      call gmem_free_inf(i10,fnm,snm,'gout')

c ***
      if (odscf) then
         if (outvv) then
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
      end if
c
      if (nprint.ne.-5) write (iwr,6040) cpu,charwall()
c
      if (omem(mainp).and.(nprint.ne.-5.or.imaxb_ic.eq.1)) then
c
c...     sort out how many integral blocks were needed for incore scf
c
         mi_bl = iblkmp + 1
         ma_bl = iblkmp + 1
         write(iwr,6050)  mi_bl,mi_bl*4096.0d0/(1024*1024), 
     +                    ma_bl,ma_bl*4096.0d0/(1024*1024) 
         call flushn(iwr)
c
c....    try only to abort root   ; might reduce the confusion
c
         if (imaxb_ic.eq.1.and.oroot()) 
     +   call caserr('incore mainfile ran out of memory')
      end if
c
      call timana(4)
      if (omem(3)) call clredx
      return
 6010 format (/' integral test counts'/1x,32('=')
     +        /' on ij shell         ',i12/
     +         ' on ijkl shells      ',i12/
     +         ' on ijkl shells & den',i12/1x,32('='))
 6020 format (' magnitudes of computed integrals '//' 2**',4x,11i8)
 6030 format (' ',8x,11i8)
 6040 format (/' end of 2-electron integral evaluation at   ',f12.2,
     +        ' seconds',a10,' wall')
 6050 format (' status of in-core Mainfile size :',/,
     +        ' minimal : ',i10,' blocks =',f12.3,' Mbyte',/,
     +        ' maximal : ',i10,' blocks =',f12.3,' Mbyte')
      end
      subroutine maxout
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
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
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
c
      logical ospbas, onocnt, ochek, opdbas, opfbas, opgbas
      integer kad, istd
      common /ijlab/ kad(mxshel),ospbas,onocnt,istd,ochek,
     +               opdbas,opfbas,opgbas
c
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer n2file, n2tape, n2blk, n2last
      integer n4file, n4tape, n4blk, n4last
      integer n6file, n6tape, n6blk, n6last
      integer n5file, n5tape, n5blk, n5last
      integer n9file, n9tape, n9blk, n9last
      integer ntfile, nttape, ntblk, ntlast
      integer n1file, n1tape, n1blk, n1last
      integer n11fil, n11tap, n11bl, n11lst
      integer n12fil, n12tap, n12bl, n12lst
      integer n13fil, n13tap, n13bl, n13lst
      integer maxbl, minbl, maxset
      common/gms_files/
     +      n2file,n2tape(20),n2blk(20),n2last(20),
     +      n4file,n4tape(20),n4blk(20),n4last(20),
     +      n6file,n6tape(20),n6blk(20),n6last(20),
     +      n5file,n5tape(20),n5blk(20),n5last(20),
     +      n9file,n9tape(20),n9blk(20),n9last(20),
     +      ntfile,nttape(20),ntblk(20),ntlast(20),
     +      n1file,n1tape(20),n1blk(20),n1last(20),
     +      n11fil,n11tap(20),n11bl(20),n11lst(20),
     +      n12fil,n12tap(20),n12bl(20),n12lst(20),
     +      n13fil,n13tap(20),n13bl(20),n13lst(20),
     +       maxbl(maxlfn),minbl(maxlfn),maxset(maxlfn)
c
c ***
      ist = isti
      jst = jsti
      kst = ksti
      lst = lsti
      iblkmp = lastb
      mfilep = lastu
      irest = 1
      m2file = mfilep
      m2last(m2file) = iblkmp
      ss = cpulft(1)
      write (iwr,6010) ss , ist , jst , kst , lst
      write (iwr,6020)
      call filprn(m2file,m2blk,m2last,m2tape)
      itask(mtask) = irest
      call revise
      write (iwr,6030)
      tim = timlim + 0.5d0
      call clredx
      return
 6010 format (/1x,58('-')//' job dumped at ',f8.2,
     +        ' seconds'//' next batch ',4i5/)
 6020 format (/' status of mainfile'/1x,18('-'))
 6030 format (///10x,27('*')///10x,'*** warning ***'/10x,
     +        'this job must be restarted'//10x,27('*')///)
      end
      subroutine pr2ei(iunit,iblock)
      implicit real*8  (a-h,p-w), integer (i-n), logical (o)
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
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
      common/craypk/i205(4,340)
      common/blkin/g(510),mword,mspace,f(4),ia(4),ja(4),ka(4),la(4)
c
      write (iwr,6040) iblock , yed(iunit)
      n = 0
      if (mword.ne.0) then
      call setsto(numlab,0,i205)
      call unpack(g(num2e+1),lab816,i205,numlab)
         do 20 m = 1 , mword
            n = n + 1
            f(n) = g(m)
      ia(n)=i205(2,m)
      ja(n)=i205(1,m)
      ka(n)=i205(4,m)
      la(n)=i205(3,m)
            if (n.ge.4) then
               write (iwr,6030) (ia(n),ja(n),ka(n),la(n),f(n),n=1,4)
               n = 0
            end if
 20      continue
         if (n.ne.0) write (iwr,6030) (ia(m),ja(m),ka(m),la(m),f(m),m=1,
     +                                n)
         write (iwr,6020) mword
      else
         write (iwr,6010)
      end if
      return
 6010 format (/3x,'endfile block')
 6020 format (/' no. of integrals =',i4)
 6030 format (4(1x,4i4,f15.8))
 6040 format (/1x,104('*')//45x,'list of block',i6,' from ',a4/45x,
     +        27('-')//4(4x,'i',3x,'j',3x,'k',3x,'l',6x,'g',8x)/1x,
     +        115('+')/)
      end
      subroutine pr2ejk(iunit,iblock)
      implicit real*8  (a-h,p-w), integer (i-n), logical (o)
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
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
      common/craypk/ij205(2,204)
      common/blkin/g(510),mword,mspace,ia(4),ja(4),fp(4),fk(4)
c
      write (iwr,6040) iblock , yed(iunit)
      n = 0
      if (mword.ne.0) then
           labss = num2ejk+num2ejk+1
      call unpack(g(labss),lab1632,ij205,numlabjk)
         do 20 m = 1 , mword
            n = n + 1
            fp(n) = g(m)
            fk(n) = g(num2ejk+m)
      ia(n)=ij205(1,m)
      ja(n)=ij205(2,m)
            if (n.ge.4) then
               write (iwr,6030) (ia(n),ja(n),fp(n),fk(n),n=1,4)
               n = 0
            end if
 20      continue
         if (n.ne.0) write (iwr,6030) (ia(m),ja(m),fp(m),fk(m),m=1,n)
         write (iwr,6020) mword
      else
         write (iwr,6010)
      end if
      return
 6010 format (/3x,'endfile block')
 6020 format (/' no. of pk-supermatrix elements =',i4)
 6030 format (4(1x,2i5,2f10.6))
 6040 format (/1x,104('*')//45x,'list of block',i6,' from ',a4/45x,
     +        27('-')//4(4x,'ij',3x,'kl',5x,'p',9x,'k',4x)/1x,123('+')/)
      end
      subroutine pr2ep(iunit,iblock)
      implicit real*8  (a-h,p-w), integer (i-n), logical (o)
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
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
      common/craypk/ij205(2,340)
      common/blkin/g(510),mword,mspace,f(4),ia(4),ja(4)
c
      write (iwr,6040) iblock , yed(iunit)
      n = 0
      if (mword.ne.0) then
      call unpack(g(num2ep+1),lab1632,ij205,numlabp)
         do 20 m = 1 , mword
            n = n + 1
            f(n) = g(m)
      ia(n)=ij205(1,m)
      ja(n)=ij205(2,m)
            if (n.ge.4) then
               write (iwr,6030) (ia(n),ja(n),f(n),n=1,4)
               n = 0
            end if
 20      continue
         if (n.ne.0) write (iwr,6030) (ia(m),ja(m),f(m),m=1,n)
         write (iwr,6020) mword
      else
         write (iwr,6010)
      end if
      return
 6010 format (/3x,'endfile block')
 6020 format (/' no. of p-supermatrix elements =',i4)
 6030 format (4(1x,2i6,f16.9))
 6040 format (/1x,104('*')//45x,'list of block',i6,' from ',a4/45x,
     +        27('-')//4(5x,'ij',4x,'kl',6x,'p',9x)/1x,115('+')/)
      end
      subroutine prt1e(q,ionsec)
c
c     ----- print 1-electron integrals resident in section
c           ionsec as a sequence of matrices   ----
c
      implicit real*8  (a-h,p-w), integer (i-n), logical (o)
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
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
      common/blkin/corev(512),potnuc(4),pint(508)
      dimension yprin(6),o1e(6),q(*)
      data yprin/'  s-','  t-','t+v-','  x-','  y-','  z-'/
c
c     process the 1e-integrals
c
      do loop =1,6
       o1e(loop) = .false.
      enddo
c
c     first restore header block
c
      call getmat(q,q,q,q,q,q,potnuc,num,o1e,ionsec)
c
      write (iwr,6010) ionsec
      write (iwr,6020) potnuc
c
c     now process each matrx in turn (S, T, T+V, X, Y, Z)
c
      do imat=1,6
        if(oprn(30+imat)) then
          o1e(imat) = .true.
          call getmat(q,q,q,q,q,q,potnuc,num,o1e,ionsec)
          write(iwr,6030)yprin(imat)
          call writel(q,num)
          o1e(imat) = .false.
        endif
      enddo
      return
 6030  format(/1x,104('*')//
     * 35x,a4,'matrix over gaussian basis set'/
     * 35x,34('-')/)
 6010 format (/1x,104('=')/45x,
     +        'list of 1-electron integrals in section',i6/45x,43('-'))
 6020 format (/' potnuc,dx,dy,dz = ',4f19.8)
      end
      subroutine roots4
c          *****   version february 16,1975   *****
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      real*8 pp, u, w, dji
      integer nroots, nroot2
      common /root/ pp,u(12),w(12),nroots,nroot2,dji(60)
c
      equivalence (u(1),rt1),(u(2),rt2),(u(3),rt3),(u(4),rt4),(u(5),rt5)
      equivalence (w(1),ww1),(w(2),ww2),(w(3),ww3),(w(4),ww4),(w(5),ww5)
      data r14,pie4/1.45303521503316d-01, 7.85398163397448d-01/
      data r24,w24/ 1.33909728812636d+00, 2.34479815323517d-01/
      data r34,w34/ 3.92696350135829d+00, 1.92704402415764d-02/
      data r44,w44/ 8.58863568901199d+00, 2.25229076750736d-04/
      if (pp.gt.15.0d+00) then
         ww1 = dsqrt(pie4/pp)
         if (pp.gt.35.0d+00) then
            if (pp.gt.53.0d+00) then
               rt1 = r14/(pp-r14)
               rt2 = r24/(pp-r24)
               rt3 = r34/(pp-r34)
               rt4 = r44/(pp-r44)
               ww4 = w44*ww1
               ww3 = w34*ww1
               ww2 = w24*ww1
               ww1 = ww1 - ww2 - ww3 - ww4
               return
            else
               g = dexp(-pp)*(pp*pp)**2
c     x=35.0 to 53.0                             nroots = 4
               rt4 = ((-2.19135070169653d-03*pp-1.19108256987623d-01)
     +               *pp-7.50238795695573d-01)*g + r44/(pp-r44)
               rt3 = ((-9.65842534508637d-04*pp-4.49822013469279d-02)
     +               *pp+6.08784033347757d-01)*g + r34/(pp-r34)
               rt2 = ((-3.62569791162153d-04*pp-9.09231717268466d-03)
     +               *pp+1.84336760556262d-01)*g + r24/(pp-r24)
               rt1 = ((-4.07557525914600d-05*pp-6.88846864931685d-04)
     +               *pp+1.74725309199384d-02)*g + r14/(pp-r14)
               ww4 = ((5.76631982000990d-06*pp-7.89187283804890d-05)
     +               *pp+3.28297971853126d-04)*g + w44*ww1
               ww3 = ((2.08294969857230d-04*pp-3.77489954837361d-03)
     +               *pp+2.09857151617436d-02)*g + w34*ww1
               ww2 = ((6.16374517326469d-04*pp-1.26711744680092d-02)
     +               *pp+8.14504890732155d-02)*g + w24*ww1
               ww1 = ww1 - ww2 - ww3 - ww4
               return
            end if
         else if (pp.gt.20.0d+00) then
c     x=20.0 to 35.0                             nroots = 4
            g = dexp(-pp)
          rt1 = ((((((-4.45711399441838d-05*pp+1.27267770241379d-03)*pp-
     +            2.36954961381262d-01)*pp+1.54330657903756d+01)
     +            *pp-5.22799159267808d+02)*pp+1.05951216669313d+04)
     +            *pp+(-2.51177235556236d+06/pp+8.72975373557709d+05)
     +            /pp-1.29194382386499d+05)*g + r14/(pp-r14)
          rt2 = (((((-7.85617372254488d-02*pp+6.35653573484868d+00)*pp-
     +            3.38296938763990d+02)*pp+1.25120495802096d+04)
     +            *pp-3.16847570511637d+05)
     +            *pp+((-1.02427466127427d+09/pp+3.70104713293016d+08)
     +            /pp-5.87119005093822d+07)/pp+5.38614211391604d+06)
     +            *g + r24/(pp-r24)
          rt3 = (((((-2.37900485051067d-01*pp+1.84122184400896d+01)*pp-
     +            1.00200731304146d+03)*pp+3.75151841595736d+04)
     +            *pp-9.50626663390130d+05)
     +            *pp+((-2.88139014651985d+09/pp+1.06625915044526d+09)
     +            /pp-1.72465289687396d+08)/pp+1.60419390230055d+07)
     +            *g + r34/(pp-r34)
          rt4 = ((((((-6.00691586407385d-04*pp-3.64479545338439d-01)*pp+
     +            1.57496131755179d+01)*pp-6.54944248734901d+02)
     +            *pp+1.70830039597097d+04)*pp-2.90517939780207d+05)
     +            *pp+(+3.49059698304732d+07/pp-1.64944522586065d+07)
     +            /pp+2.96817940164703d+06)*g + r44/(pp-r44)
            if (pp.le.25.0d+00) ww4 = (((((((2.33766206773151d-07*pp-
     +                               3.81542906607063d-05)
     +                               *pp+3.51416601267000d-03)
     +                               *pp-1.66538571864728d-01)
     +                               *pp+4.80006136831847d+00)
     +                               *pp-8.73165934223603d+01)
     +                               *pp+9.77683627474638d+02)
     +                               *pp+1.66000945117640d+04/pp-
     +                               6.14479071209961d+03)*g + w44*ww1
            if (pp.gt.25.0d+00) ww4 = ((((((5.74245945342286d-06*pp-
     +                               7.58735928102351d-05)
     +                               *pp+2.35072857922892d-04)
     +                               *pp-3.78812134013125d-03)
     +                               *pp+3.09871652785805d-01)
     +                               *pp-7.11108633061306d+00)
     +                               *pp+5.55297573149528d+01)
     +                               *g + w44*ww1
           ww3 = ((((((2.36392855180768d-04*pp-9.16785337967013d-03)*pp+
     +            4.62186525041313d-01)*pp-1.96943786006540d+01)
     +            *pp+4.99169195295559d+02)*pp-6.21419845845090d+03)
     +            *pp+((+5.21445053212414d+07/pp-1.34113464389309d+07)
     +            /pp+1.13673298305631d+06)/pp-2.81501182042707d+03)
     +            *g + w34*ww1
           ww2 = ((((((7.29841848989391d-04*pp-3.53899555749875d-02)*pp+
     +            2.07797425718513d+00)*pp-1.00464709786287d+02)
     +            *pp+3.15206108877819d+03)*pp-6.27054715090012d+04)
     +            *pp+(+1.54721246264919d+07/pp-5.26074391316381d+06)
     +            /pp+7.67135400969617d+05)*g + w24*ww1
           ww1 = ((1.9623264149430d-01/pp-4.9695241464490d-01)
     +            /pp-6.0156581186481d-05)*g + ww1 - ww2 - ww3 - ww4
            return
         else
c     x=15.0 to 20.0                             nroots = 4
            q = pp - 17.5d+00
            rt1 = (((((((((((4.36701759531398d-17*q-1.12860600219889d-16
     +            )*q-6.149849164164d-15)*q+5.820231579541d-14)
     +            *q+4.396602872143d-13)*q-1.24330365320172d-11)
     +            *q+6.71083474044549d-11)*q+2.43865205376067d-10)
     +            *q+1.67559587099969d-08)*q-9.32738632357572d-07)
     +            *q+2.39030487004977d-05)*q-4.68648206591515d-04)
     +            *q + 8.34977776583956d-03
            rt2 = (((((((((((4.98913142288158d-16*q-2.60732537093612d-16
     +            )*q-7.775156445127d-14)*q+5.766105220086d-13)
     +            *q+6.432696729600d-12)*q-1.39571683725792d-10)
     +            *q+5.95451479522191d-10)*q+2.42471442836205d-09)
     +            *q+2.47485710143120d-07)*q-1.14710398652091d-05)
     +            *q+2.71252453754519d-04)*q-4.96812745851408d-03)
     +            *q + 8.26020602026780d-02
            rt3 = (((((((((((1.91498302509009d-15*q+1.48840394311115d-14
     +            )*q-4.316925145767d-13)*q+1.186495793471d-12)
     +            *q+4.615806713055d-11)*q-5.54336148667141d-10)
     +            *q+3.48789978951367d-10)*q-2.79188977451042d-09)
     +            *q+2.09563208958551d-06)*q-6.76512715080324d-05)
     +            *q+1.32129867629062d-03)*q-2.05062147771513d-02)
     +            *q + 2.88068671894324d-01
            rt4 = (((((((((((-5.43697691672942d-15*q-
     +            1.12483395714468d-13)*q+2.826607936174d-12)
     +            *q-1.266734493280d-11)*q-4.258722866437d-10)
     +            *q+9.45486578503261d-09)*q-5.86635622821309d-08)
     +            *q-1.28835028104639d-06)*q+4.41413815691885d-05)
     +            *q-7.61738385590776d-04)*q+9.66090902985550d-03)
     +            *q-1.01410568057649d-01)*q + 9.54714798156712d-01
            ww4 = ((((((((((((-7.56882223582704d-19*q+
     +            7.53541779268175d-18)*q-1.157318032236d-16)
     +            *q+2.411195002314d-15)*q-3.601794386996d-14)
     +            *q+4.082150659615d-13)*q-4.289542980767d-12)
     +            *q+5.086829642731d-11)*q-6.35435561050807d-10)
     +            *q+6.82309323251123d-09)*q-5.63374555753167d-08)
     +            *q+3.57005361100431d-07)*q-2.40050045173721d-06)
     +            *q + 4.94171300536397d-05
            ww3 = (((((((((((-5.54451040921657d-17*q+
     +            2.68748367250999d-16)*q+1.349020069254d-14)
     +            *q-2.507452792892d-13)*q+1.944339743818d-12)
     +            *q-1.29816917658823d-11)*q+3.49977768819641d-10)
     +            *q-8.67270669346398d-09)*q+1.31381116840118d-07)
     +            *q-1.36790720600822d-06)*q+1.19210697673160d-05)
     +            *q-1.42181943986587d-04)*q + 4.12615396191829d-03
            ww2 = (((((((((((-1.86506057729700d-16*q+
     +            1.16661114435809d-15)*q+2.563712856363d-14)
     +            *q-4.498350984631d-13)*q+1.765194089338d-12)
     +            *q+9.04483676345625d-12)*q+4.98930345609785d-10)
     +            *q-2.11964170928181d-08)*q+3.98295476005614d-07)
     +            *q-5.49390160829409d-06)*q+7.74065155353262d-05)
     +            *q-1.48201933009105d-03)*q + 4.97836392625268d-02
            ww1 = ((1.9623264149430d-01/pp-4.9695241464490d-01)
     +            /pp-6.0156581186481d-05)*dexp(-pp) + ww1 - ww2 - ww3 -
     +            ww4
            return
         end if
      else if (pp.gt.5.0d+00) then
         if (pp.gt.10.0d+00) then
c     x=10.0 to 15.0                             nroots = 4
            q = pp - 12.5d+00
            rt1 = (((((((((((4.94869622744119d-17*q+8.03568805739160d-16
     +            )*q-5.599125915431d-15)*q-1.378685560217d-13)
     +            *q+7.006511663249d-13)*q+1.30391406991118d-11)
     +            *q+8.06987313467541d-11)*q-5.20644072732933d-09)
     +            *q+7.72794187755457d-08)*q-1.61512612564194d-06)
     +            *q+4.15083811185831d-05)*q-7.87855975560199d-04)
     +            *q + 1.14189319050009d-02
            rt2 = (((((((((((4.89224285522336d-16*q+1.06390248099712d-14
     +            )*q-5.446260182933d-14)*q-1.613630106295d-12)
     +            *q+3.910179118937d-12)*q+1.90712434258806d-10)
     +            *q+8.78470199094761d-10)*q-5.97332993206797d-08)
     +            *q+9.25750831481589d-07)*q-2.02362185197088d-05)
     +            *q+4.92341968336776d-04)*q-8.68438439874703d-03)
     +            *q + 1.15825965127958d-01
            rt3 = ((((((((((6.12419396208408d-14*q+1.12328861406073d-13)
     +            *q-9.051094103059d-12)*q-4.781797525341d-11)
     +            *q+1.660828868694d-09)*q+4.499058798868d-10)
     +            *q-2.519549641933d-07)*q+4.977444040180d-06)
     +            *q-1.25858350034589d-04)*q+2.70279176970044d-03)
     +            *q-3.99327850801083d-02)*q + 4.33467200855434d-01
            rt4 = (((((((((((4.63414725924048d-14*q-4.72757262693062d-14
     +            )*q-1.001926833832d-11)*q+6.074107718414d-11)
     +            *q+1.576976911942d-09)*q-2.01186401974027d-08)
     +            *q-1.84530195217118d-07)*q+5.02333087806827d-06)
     +            *q+9.66961790843006d-06)*q-1.58522208889528d-03)
     +            *q+2.80539673938339d-02)*q-2.78953904330072d-01)
     +            *q + 1.82835655238235d+00
            ww4 = (((((((((((((2.90401781000996d-18*q-
     +            4.63389683098251d-17)*q+6.274018198326d-16)
     +            *q-8.936002188168d-15)*q+1.194719074934d-13)
     +            *q-1.45501321259466d-12)*q+1.64090830181013d-11)
     +            *q-1.71987745310181d-10)*q+1.63738403295718d-09)
     +            *q-1.39237504892842d-08)*q+1.06527318142151d-07)
     +            *q-7.27634957230524d-07)*q+4.12159381310339d-06)
     +            *q-1.74648169719173d-05)*q + 8.50290130067818d-05
            ww3 = ((((((((((((-4.19569145459480d-17*q+
     +            5.94344180261644d-16)*q-1.148797566469d-14)
     +            *q+1.881303962576d-13)*q-2.413554618391d-12)
     +            *q+3.372127423047d-11)*q-4.933988617784d-10)
     +            *q+6.116545396281d-09)*q-6.69965691739299d-08)
     +            *q+7.52380085447161d-07)*q-8.08708393262321d-06)
     +            *q+6.88603417296672d-05)*q-4.67067112993427d-04)
     +            *q + 5.42313365864597d-03
            ww2 = ((((((((((-6.22272689880615d-15*q+1.04126809657554d-13
     +            )*q-6.842418230913d-13)*q+1.576841731919d-11)
     +            *q-4.203948834175d-10)*q+6.287255934781d-09)
     +            *q-8.307159819228d-08)*q+1.356478091922d-06)
     +            *q-2.08065576105639d-05)*q+2.52396730332340d-04)
     +            *q-2.94484050194539d-03)*q + 6.01396183129168d-02
            ww1 = (((-1.8784686463512d-01/pp+2.2991849164985d-01)
     +            /pp-4.9893752514047d-01)/pp-2.1916512131607d-05)
     +            *dexp(-pp) + dsqrt(pie4/pp) - ww4 - ww3 - ww2
            return
         else
            q = pp - 7.5d+00
c     x=5.0 to 10.0                              nroots = 4
            rt1 = (((((((((4.64217329776215d-15*q-6.27892383644164d-15)*
     +            q+3.462236347446d-13)*q-2.927229355350d-11)
     +            *q+5.090355371676d-10)*q-9.97272656345253d-09)
     +            *q+2.37835295639281d-07)*q-4.60301761310921d-06)
     +            *q+8.42824204233222d-05)*q-1.37983082233081d-03)
     +            *q + 1.66630865869375d-02
            rt2 = (((((((((2.93981127919047d-14*q+8.47635639065744d-13)*
     +            q-1.446314544774d-11)*q-6.149155555753d-12)
     +            *q+8.484275604612d-10)*q-6.10898827887652d-08)
     +            *q+2.39156093611106d-06)*q-5.35837089462592d-05)
     +            *q+1.00967602595557d-03)*q-1.57769317127372d-02)
     +            *q + 1.74853819464285d-01
            rt3 = ((((((((((2.93523563363000d-14*q-6.40041776667020d-14)
     +            *q-2.695740446312d-12)*q+1.027082960169d-10)
     +            *q-5.822038656780d-10)*q-3.159991002539d-08)
     +            *q+4.327249251331d-07)*q+4.856768455119d-06)
     +            *q-2.54617989427762d-04)*q+5.54843378106589d-03)
     +            *q-7.95013029486684d-02)*q + 7.20206142703162d-01
            rt4 = (((((((((((-1.62212382394553d-14*q+
     +            7.68943641360593d-13)*q+5.764015756615d-12)
     +            *q-1.380635298784d-10)*q-1.476849808675d-09)
     +            *q+1.84347052385605d-08)*q+3.34382940759405d-07)
     +            *q-1.39428366421645d-06)*q-7.50249313713996d-05)
     +            *q-6.26495899187507d-04)*q+4.69716410901162d-02)
     +            *q-6.66871297428209d-01)*q + 4.11207530217806d+00
            ww1 = ((((((((((-1.65995045235997d-15*q+6.91838935879598d-14
     +            )*q-9.131223418888d-13)*q+1.403341829454d-11)
     +            *q-3.672235069444d-10)*q+6.366962546990d-09)
     +            *q-1.039220021671d-07)*q+1.959098751715d-06)
     +            *q-3.33474893152939d-05)*q+5.72164211151013d-04)
     +            *q-1.05583210553392d-02)*q + 2.26696066029591d-01
            ww2 = ((((((((((((-3.57248951192047d-16*q+
     +            6.25708409149331d-15)*q-9.657033089714d-14)
     +            *q+1.507864898748d-12)*q-2.332522256110d-11)
     +            *q+3.428545616603d-10)*q-4.698730937661d-09)
     +            *q+6.219977635130d-08)*q-7.83008889613661d-07)
     +            *q+9.08621687041567d-06)*q-9.86368311253873d-05)
     +            *q+9.69632496710088d-04)*q-8.14594214284187d-03)
     +            *q + 8.50218447733457d-02
            ww3 = (((((((((((((1.64742458534277d-16*q-
     +            2.68512265928410d-15)*q+3.788890667676d-14)
     +            *q-5.508918529823d-13)*q+7.555896810069d-12)
     +            *q-9.69039768312637d-11)*q+1.16034263529672d-09)
     +            *q-1.28771698573873d-08)*q+1.31949431805798d-07)
     +            *q-1.23673915616005d-06)*q+1.04189803544936d-05)
     +            *q-7.79566003744742d-05)*q+5.03162624754434d-04)
     +            *q-2.55138844587555d-03)*q + 1.13250730954014d-02
            ww4 = ((((((((((((((-1.55714130075679d-17*q+
     +            2.57193722698891d-16)*q-3.626606654097d-15)
     +            *q+5.234734676175d-14)*q-7.067105402134d-13)
     +            *q+8.793512664890d-12)*q-1.006088923498d-10)
     +            *q+1.050565098393d-09)*q-9.91517881772662d-09)
     +            *q+8.35835975882941d-08)*q-6.19785782240693d-07)
     +            *q+3.95841149373135d-06)*q-2.11366761402403d-05)
     +            *q+9.00474771229507d-05)*q-2.78777909813289d-04)
     +            *q + 5.26543779837487d-04
            return
         end if
      else if (pp.gt.1.0d+00) then
c     x= 1.0 to 5.0                              nroots = 4
         q = pp - 3.0d+00
         rt1 = (((((((((-1.48570633747284d-15*q-1.33273068108777d-13)*q+
     +         4.068543696670d-12)*q-9.163164161821d-11)
     +         *q+2.046819017845d-09)*q-4.03076426299031d-08)
     +         *q+7.29407420660149d-07)*q-1.23118059980833d-05)
     +         *q+1.88796581246938d-04)*q-2.53262912046853d-03)
     +         *q + 2.51198234505021d-02
         rt2 = (((((((((1.35830583483312d-13*q-2.29772605964836d-12)*q-
     +         3.821500128045d-12)*q+6.844424214735d-10)
     +         *q-1.048063352259d-08)*q+1.50083186233363d-08)
     +         *q+3.48848942324454d-06)*q-1.08694174399193d-04)
     +         *q+2.08048885251999d-03)*q-2.91205805373793d-02)
     +         *q + 2.72276489515713d-01
         rt3 = (((((((((5.02799392850289d-13*q+1.07461812944084d-11)*q-
     +         1.482277886411d-10)*q-2.153585661215d-09)
     +         *q+3.654087802817d-08)*q+5.15929575830120d-07)
     +         *q-9.52388379435709d-06)*q-2.16552440036426d-04)
     +         *q+9.03551469568320d-03)*q-1.45505469175613d-01)
     +         *q + 1.21449092319186d+00
         rt4 = (((((((((-1.08510370291979d-12*q+6.41492397277798d-11)*q+
     +         7.542387436125d-10)*q-2.213111836647d-09)
     +         *q-1.448228963549d-07)*q-1.95670833237101d-06)
     +         *q-1.07481314670844d-05)*q+1.49335941252765d-04)
     +         *q+4.87791531990593d-02)*q-1.10559909038653d+00)
     +         *q + 8.09502028611780d+00
         ww1 = ((((((((((-4.65801912689961d-14*q+7.58669507106800d-13)*q
     +         -1.186387548048d-11)*q+1.862334710665d-10)
     +         *q-2.799399389539d-09)*q+4.148972684255d-08)
     +         *q-5.933568079600d-07)*q+8.168349266115d-06)
     +         *q-1.08989176177409d-04)*q+1.41357961729531d-03)
     +         *q-1.87588361833659d-02)*q + 2.89898651436026d-01
         ww2 = ((((((((((((-1.46345073267549d-14*q+2.25644205432182d-13)
     +         *q-3.116258693847d-12)*q+4.321908756610d-11)
     +         *q-5.673270062669d-10)*q+7.006295962960d-09)
     +         *q-8.120186517000d-08)*q+8.775294645770d-07)
     +         *q-8.77829235749024d-06)*q+8.04372147732379d-05)
     +         *q-6.64149238804153d-04)*q+4.81181506827225d-03)
     +         *q-2.88982669486183d-02)*q + 1.56247249979288d-01
         ww3 = (((((((((((((9.06812118895365d-15*q-1.40541322766087d-13)
     +         *q+1.919270015269d-12)*q-2.605135739010d-11)
     +         *q+3.299685839012d-10)*q-3.86354139348735d-09)
     +         *q+4.16265847927498d-08)*q-4.09462835471470d-07)
     +         *q+3.64018881086111d-06)*q-2.88665153269386d-05)
     +         *q+2.00515819789028d-04)*q-1.18791896897934d-03)
     +         *q+5.75223633388589d-03)*q-2.09400418772687d-02)
     +         *q + 4.85368861938873d-02
         ww4 = ((((((((((((((-9.74835552342257d-16*q+
     +         1.57857099317175d-14)*q-2.249993780112d-13)
     +         *q+3.173422008953d-12)*q-4.161159459680d-11)
     +         *q+5.021343560166d-10)*q-5.545047534808d-09)
     +         *q+5.554146993491d-08)*q-4.99048696190133d-07)
     +         *q+3.96650392371311d-06)*q-2.73816413291214d-05)
     +         *q+1.60106988333186d-04)*q-7.64560567879592d-04)
     +         *q+2.81330044426892d-03)*q-7.16227030134947d-03)
     +         *q + 9.66077262223353d-03
         return
      else if (pp.gt.3.0d-07) then
c     x=0.0 to 1.0                               nroots = 4
        rt1 = ((((((-1.95309614628539d-10*pp+5.19765728707592d-09)*pp-
     +         1.01756452250573d-07)*pp+1.72365935872131d-06)
     +         *pp-2.61203523522184d-05)*pp+3.52921308769880d-04)
     +         *pp-4.09645850658433d-03)*pp + 3.48198973061469d-02
        rt2 = (((((-1.89554881382342d-08*pp+3.07583114342365d-07)*pp+
     +         1.270981734393d-06)*pp-1.417298563884d-04)
     +         *pp+3.226979163176d-03)*pp-4.48902570678178d-02)
     +         *pp + 3.81567185080039d-01
        rt3 = ((((((1.77280535300416d-09*pp+3.36524958870615d-08)*pp-
     +         2.58341529013893d-07)*pp-1.13644895662320d-05)
     +         *pp-7.91549618884063d-05)*pp+1.03825827346828d-02)
     +         *pp-2.04389090525137d-01)*pp + 1.73730726945889d+00
        rt4 = (((((-5.61188882415248d-08*pp-2.49480733072460d-07)*pp+
     +         3.428685057114d-06)*pp+1.679007454539d-04)
     +         *pp+4.722855585715d-02)*pp-1.39368301737828d+00)
     +         *pp + 1.18463056481543d+01
        ww1 = ((((((-1.14649303201279d-08*pp+1.88015570196787d-07)*pp-
     +         2.33305875372323d-06)*pp+2.68880044371597d-05)
     +         *pp-2.94268428977387d-04)*pp+3.06548909776613d-03)
     +         *pp-3.13844305680096d-02)*pp + 3.62683783378335d-01
        ww2 = ((((((((-4.11720483772634d-09*pp+6.54963481852134d-08)*pp-
     +         7.20045285129626d-07)*pp+6.93779646721723d-06)
     +         *pp-6.05367572016373d-05)*pp+4.74241566251899d-04)
     +         *pp-3.26956188125316d-03)*pp+1.91883866626681d-02)
     +         *pp-8.98046242565811d-02)*pp + 3.13706645877886d-01
        ww3 = ((((((((-3.41688436990215d-08*pp+5.07238960340773d-07)*pp-
     +         5.01675628408220d-06)*pp+4.20363420922845d-05)
     +         *pp-3.08040221166823d-04)*pp+1.94431864731239d-03)
     +         *pp-1.02477820460278d-02)*pp+4.28670143840073d-02)
     +         *pp-1.29314370962569d-01)*pp + 2.22381034453369d-01
        ww4 = (((((((((4.99660550769508d-09*pp-7.94585963310120d-08)*pp+
     +         8.359072409485d-07)*pp-7.422369210610d-06)
     +         *pp+5.763374308160d-05)*pp-3.86645606718233d-04)
     +         *pp+2.18417516259781d-03)*pp-9.99791027771119d-03)
     +         *pp+3.48791097377370d-02)*pp-8.28299075413889d-02)
     +         *pp + 1.01228536290376d-01
         return
      else
c     x is approximately zero.                   nroots = 4
         rt1 = 3.48198973061471d-02 - 4.09645850660395d-03*pp
         rt2 = 3.81567185080042d-01 - 4.48902570656719d-02*pp
         rt3 = 1.73730726945891d+00 - 2.04389090547327d-01*pp
         rt4 = 1.18463056481549d+01 - 1.39368301742312d+00*pp
         ww1 = 3.62683783378362d-01 - 3.13844305713928d-02*pp
         ww2 = 3.13706645877886d-01 - 8.98046242557724d-02*pp
         ww3 = 2.22381034453372d-01 - 1.29314370958973d-01*pp
         ww4 = 1.01228536290376d-01 - 8.28299075414321d-02*pp
         return
      end if
      end
      subroutine roots5
c          *****   version  february 27,1975   *****
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      real*8 pp, u, w, dji
      integer nroots, nroot2
      common /root/ pp,u(12),w(12),nroots,nroot2,dji(60)
c
      equivalence (u(1),rt1),(u(2),rt2),(u(3),rt3),(u(4),rt4),(u(5),rt5)
      equivalence (w(1),ww1),(w(2),ww2),(w(3),ww3),(w(4),ww4),(w(5),ww5)
      data r15,pie4/1.17581320211778d-01, 7.85398163397448d-01/
      data r25,w25/ 1.07456201243690d+00, 2.70967405960535d-01/
      data r35,w35/ 3.08593744371754d+00, 3.82231610015404d-02/
      data r45,w45/ 6.41472973366203d+00, 1.51614186862443d-03/
      data r55,w55/ 1.18071894899717d+01, 8.62130526143657d-06/
      if (pp.gt.15.0d+00) then
         if (pp.gt.25.0d+00) then
            ww1 = dsqrt(pie4/pp)
            if (pp.le.40.0d+00) then
c     x=25.0 to 40.0                             nroots = 5
               g = dexp(-pp)
              rt1 = ((((((((-1.73363958895356d-06*pp+
     +               1.19921331441483d-04)*pp-1.59437614121125d-02)
     +               *pp+1.13467897349442d+00)*pp-4.47216460864586d+01)
     +               *pp+1.06251216612604d+03)*pp-1.52073917378512d+04)
     +               *pp+1.20662887111273d+05)*pp-4.07186366852475d+05)
     +               *g + r15/(pp-r15)
              rt2 = ((((((((-1.60102542621710d-05*pp+
     +               1.10331262112395d-03)*pp-1.50043662589017d-01)
     +               *pp+1.05563640866077d+01)*pp-4.10468817024806d+02)
     +               *pp+9.62604416506819d+03)*pp-1.35888069838270d+05)
     +               *pp+1.06107577038340d+06)*pp-3.51190792816119d+06)
     +               *g + r25/(pp-r25)
              rt3 = ((((((((-4.48880032128422d-05*pp+
     +              2.69025112122177d-03)*pp-4.01048115525954d-01)
     +               *pp+2.78360021977405d+01)*pp-1.04891729356965d+03)
     +               *pp+2.36985942687423d+04)*pp-3.19504627257548d+05)
     +               *pp+2.34879693563358d+06)*pp-7.16341568174085d+06)
     +               *g + r35/(pp-r35)
              rt4 = ((((((((-6.38526371092582d-05*pp-
     +               2.29263585792626d-03)*pp-7.65735935499627d-02)
     +               *pp+9.12692349152792d+00)*pp-2.32077034386717d+02)
     +               *pp+2.81839578728845d+02)*pp+9.59529683876419d+04)
     +               *pp-1.77638956809518d+06)*pp+1.02489759645410d+07)
     +               *g + r45/(pp-r45)
              rt5 = ((((((((-3.59049364231569d-05*pp-
     +               2.25963977930044d-02)*pp+1.12594870794668d+00)
     +               *pp-4.56752462103909d+01)*pp+1.05804526830637d+03)
     +               *pp-1.16003199605875d+04)*pp-4.07297627297272d+04)
     +               *pp+2.22215528319857d+06)*pp-1.61196455032613d+07)
     +               *g + r55/(pp-r55)
              ww5 = (((((((((-4.61100906133970d-10*pp+
     +               1.43069932644286d-07)*pp-1.63960915431080d-05)
     +               *pp+1.15791154612838d-03)*pp-5.30573476742071d-02)
     +               *pp+1.61156533367153d+00)*pp-3.23248143316007d+01)
     +               *pp+4.12007318109157d+02)*pp-3.02260070158372d+03)
     +               *pp+9.71575094154768d+03)*g + w55*ww1
              ww4 = (((((((((-2.40799435809950d-08*pp+
     +               8.12621667601546d-06)*pp-9.04491430884113d-04)
     +               *pp+6.37686375770059d-02)*pp-2.96135703135647d+00)
     +               *pp+9.15142356996330d+01)*pp-1.86971865249111d+03)
     +               *pp+2.42945528916947d+04)*pp-1.81852473229081d+05)
     +               *pp+5.96854758661427d+05)*g + w45*ww1
              ww3 = ((((((((1.83574464457207d-05*pp-1.54837969489927d-03
     +               )*pp+1.18520453711586d-01)*pp-6.69649981309161d+00)
     +               *pp+2.44789386487321d+02)*pp-5.68832664556359d+03)
     +               *pp+8.14507604229357d+04)*pp-6.55181056671474d+05)
     +               *pp+2.26410896607237d+06)*g + w35*ww1
              ww2 = ((((((((2.77778345870650d-05*pp-2.22835017655890d-03
     +               )*pp+1.61077633475573d-01)*pp-8.96743743396132d+00)
     +               *pp+3.28062687293374d+02)*pp-7.65722701219557d+03)
     +               *pp+1.10255055017664d+05)*pp-8.92528122219324d+05)
     +               *pp+3.10638627744347d+06)*g + w25*ww1
              ww1 = ww1 - 0.01962d+00*g - ww2 - ww3 - ww4 - ww5
               return
            else if (pp.gt.59.0d+00) then
c     x=59.0 to infinity                         nroots = 5
               rt1 = r15/(pp-r15)
               rt2 = r25/(pp-r25)
               rt3 = r35/(pp-r35)
               rt4 = r45/(pp-r45)
               rt5 = r55/(pp-r55)
               ww2 = w25*ww1
               ww3 = w35*ww1
               ww4 = w45*ww1
               ww5 = w55*ww1
               ww1 = ww1 - ww2 - ww3 - ww4 - ww5
               return
            else
c     x=40.0 to 59.0                             nroots = 5
               ppp = pp**3
               g = ppp*dexp(-pp)
               rt1 = (((-2.43758528330205d-02*pp+2.07301567989771d+00)
     +               *pp-6.45964225381113d+01)*pp+7.14160088655470d+02)
     +               *g + r15/(pp-r15)
               rt2 = (((-2.28861955413636d-01*pp+1.93190784733691d+01)
     +               *pp-5.99774730340912d+02)*pp+6.61844165304871d+03)
     +               *g + r25/(pp-r25)
               rt3 = (((-6.95053039285586d-01*pp+5.76874090316016d+01)
     +               *pp-1.77704143225520d+03)*pp+1.95366082947811d+04)
     +               *g + r35/(pp-r35)
               rt4 = (((-1.58072809087018d+00*pp+1.27050801091948d+02)
     +               *pp-3.86687350914280d+03)*pp+4.23024828121420d+04)
     +               *g + r45/(pp-r45)
               rt5 = (((-3.33963830405396d+00*pp+2.51830424600204d+02)
     +               *pp-7.57728527654961d+03)*pp+8.21966816595690d+04)
     +               *g + r55/(pp-r55)
               g = ppp*g
               ww5 = ((1.35482430510942d-08*pp-3.27722199212781d-07)
     +               *pp+2.41522703684296d-06)*g + w55*ww1
               ww4 = ((1.23464092261605d-06*pp-3.55224564275590d-05)
     +               *pp+3.03274662192286d-04)*g + w45*ww1
               ww3 = ((1.34547929260279d-05*pp-4.19389884772726d-04)
     +               *pp+3.87706687610809d-03)*g + w35*ww1
               ww2 = ((2.09539509123135d-05*pp-6.87646614786982d-04)
     +               *pp+6.68743788585688d-03)*g + w25*ww1
               ww1 = ww1 - ww2 - ww3 - ww4 - ww5
               return
            end if
         else if (pp.gt.20.0d+00) then
c     x=20.0 to 25.0                             nroots = 5
            q = pp - 22.5d+00
            rt1 = (((((((((-1.13927848238726d-15*q+7.39404133595713d-15)
     +            *q+1.445982921243d-13)*q-2.676703245252d-12)
     +            *q+5.823521627177d-12)*q+2.17264723874381d-10)
     +            *q+3.56242145897468d-09)*q-3.03763737404491d-07)
     +            *q+9.46859114120901d-06)*q-2.30896753853196d-04)
     +            *q + 5.24663913001114d-03
            rt2 = ((((((((((2.89872355524581d-16*q-1.22296292045864d-14)
     +            *q+6.184065097200d-14)*q+1.649846591230d-12)
     +            *q-2.729713905266d-11)*q+3.709913790650d-11)
     +            *q+2.216486288382d-09)*q+4.616160236414d-08)
     +            *q-3.32380270861364d-06)*q+9.84635072633776d-05)
     +            *q-2.30092118015697d-03)*q + 5.00845183695073d-02
            rt3 = ((((((((((1.97068646590923d-15*q-4.89419270626800d-14)
     +            *q+1.136466605916d-13)*q+7.546203883874d-12)
     +            *q-9.635646767455d-11)*q-8.295965491209d-11)
     +            *q+7.534109114453d-09)*q+2.699970652707d-07)
     +            *q-1.42982334217081d-05)*q+3.78290946669264d-04)
     +            *q-8.03133015084373d-03)*q + 1.58689469640791d-01
            rt4 = ((((((((((1.33642069941389d-14*q-1.55850612605745d-13)
     +            *q-7.522712577474d-13)*q+3.209520801187d-11)
     +            *q-2.075594313618d-10)*q-2.070575894402d-09)
     +            *q+7.323046997451d-09)*q+1.851491550417d-06)
     +            *q-6.37524802411383d-05)*q+1.36795464918785d-03)
     +            *q-2.42051126993146d-02)*q + 3.97847167557815d-01
            rt5 = ((((((((((-6.07053986130526d-14*q+1.04447493138843d-12
     +            )*q-4.286617818951d-13)*q-2.632066100073d-10)
     +            *q+4.804518986559d-09)*q-1.835675889421d-08)
     +            *q-1.068175391334d-06)*q+3.292234974141d-05)
     +            *q-5.94805357558251d-04)*q+8.29382168612791d-03)
     +            *q-9.93122509049447d-02)*q + 1.09857804755042d+00
            ww1 = (((((((((-9.10338640266542d-15*q+1.00438927627833d-13)
     +            *q+7.817349237071d-13)*q-2.547619474232d-11)
     +            *q+1.479321506529d-10)*q+1.52314028857627d-09)
     +            *q+9.20072040917242d-09)*q-2.19427111221848d-06)
     +            *q+8.65797782880311d-05)*q-2.82718629312875d-03)
     +            *q + 1.28718310443295d-01
            ww2 = (((((((((5.52380927618760d-15*q-6.43424400204124d-14)*
     +            q-2.358734508092d-13)*q+8.261326648131d-12)
     +            *q+9.229645304956d-11)*q-5.68108973828949d-09)
     +            *q+1.22477891136278d-07)*q-2.11919643127927d-06)
     +            *q+4.23605032368922d-05)*q-1.14423444576221d-03)
     +            *q + 5.06607252890186d-02
            ww3 = (((((((((3.99457454087556d-15*q-5.11826702824182d-14)*
     +            q-4.157593182747d-14)*q+4.214670817758d-12)
     +            *q+6.705582751532d-11)*q-3.36086411698418d-09)
     +            *q+6.07453633298986d-08)*q-7.40736211041247d-07)
     +            *q+8.84176371665149d-06)*q-1.72559275066834d-04)
     +            *q + 7.16639814253567d-03
            ww4 = (((((((((((-2.14649508112234d-18*q-
     +            2.45525846412281d-18)*q+6.126212599772d-16)
     +            *q-8.526651626939d-15)*q+4.826636065733d-14)
     +            *q-3.39554163649740d-13)*q+1.67070784862985d-11)
     +            *q-4.42671979311163d-10)*q+6.77368055908400d-09)
     +            *q-7.03520999708859d-08)*q+6.04993294708874d-07)
     +            *q-7.80555094280483d-06)*q + 2.85954806605017d-04
            ww5 = ((((((((((((-5.63938733073804d-21*q+
     +            6.92182516324628d-20)*q-1.586937691507d-18)
     +            *q+3.357639744582d-17)*q-4.810285046442d-16)
     +            *q+5.386312669975d-15)*q-6.117895297439d-14)
     +            *q+8.441808227634d-13)*q-1.18527596836592d-11)
     +            *q+1.36296870441445d-10)*q-1.17842611094141d-09)
     +            *q+7.80430641995926d-09)*q-5.97767417400540d-08)
     +            *q + 1.65186146094969d-06
            return
         else
c     x=15.0 to 20.0                             nroots = 5
            q = pp - 17.5d+00
            rt1 = ((((((((((1.91875764545740d-16*q+7.8357401095707d-16)*
     +            q-3.260875931644d-14)*q-1.186752035569d-13)
     +            *q+4.275180095653d-12)*q+3.357056136731d-11)
     +            *q-1.123776903884d-09)*q+1.231203269887d-08)
     +            *q-3.99851421361031d-07)*q+1.45418822817771d-05)
     +            *q-3.49912254976317d-04)*q + 6.67768703938812d-03
            rt2 = ((((((((((2.02778478673555d-15*q+1.01640716785099d-14)
     +            *q-3.385363492036d-13)*q-1.615655871159d-12)
     +            *q+4.527419140333d-11)*q+3.853670706486d-10)
     +            *q-1.184607130107d-08)*q+1.347873288827d-07)
     +            *q-4.47788241748377d-06)*q+1.54942754358273d-04)
     +            *q-3.55524254280266d-03)*q + 6.44912219301603d-02
            rt3 = ((((((((((7.79850771456444d-15*q+6.00464406395001d-14)
     +            *q-1.249779730869d-12)*q-1.020720636353d-11)
     +            *q+1.814709816693d-10)*q+1.766397336977d-09)
     +            *q-4.603559449010d-08)*q+5.863956443581d-07)
     +            *q-2.03797212506691d-05)*q+6.31405161185185d-04)
     +            *q-1.30102750145071d-02)*q + 2.10244289044705d-01
            rt4 = (((((((((((-2.92397030777912d-15*q+
     +            1.94152129078465d-14)*q+4.859447665850d-13)
     +            *q-3.217227223463d-12)*q-7.484522135512d-11)
     +            *q+7.19101516047753d-10)*q+6.88409355245582d-09)
     +            *q-1.44374545515769d-07)*q+2.74941013315834d-06)
     +            *q-1.02790452049013d-04)*q+2.59924221372643d-03)
     +            *q-4.35712368303551d-02)*q + 5.62170709585029d-01
            rt5 = (((((((((((1.17976126840060d-14*q+1.24156229350669d-13
     +            )*q-3.892741622280d-12)*q-7.755793199043d-12)
     +            *q+9.492190032313d-10)*q-4.98680128123353d-09)
     +            *q-1.81502268782664d-07)*q+2.69463269394888d-06)
     +            *q+2.50032154421640d-05)*q-1.33684303917681d-03)
     +            *q+2.29121951862538d-02)*q-2.45653725061323d-01)
     +            *q + 1.89999883453047d+00
            ww1 = ((((((((((1.74841995087592d-15*q-6.95671892641256d-16)
     +            *q-3.000659497257d-13)*q+2.021279817961d-13)
     +            *q+3.853596935400d-11)*q+1.461418533652d-10)
     +            *q-1.014517563435d-08)*q+1.132736008979d-07)
     +            *q-2.86605475073259d-06)*q+1.21958354908768d-04)
     +            *q-3.86293751153466d-03)*q + 1.45298342081522d-01
            ww2 = ((((((((((-1.11199320525573d-15*q+1.85007587796671d-15
     +            )*q+1.220613939709d-13)*q+1.275068098526d-12)
     +            *q-5.341838883262d-11)*q+6.161037256669d-10)
     +            *q-1.009147879750d-08)*q+2.907862965346d-07)
     +            *q-6.12300038720919d-06)*q+1.00104454489518d-04)
     +            *q-1.80677298502757d-03)*q + 5.78009914536630d-02
            ww3 = ((((((((((-9.49816486853687d-16*q+6.67922080354234d-15
     +            )*q+2.606163540537d-15)*q+1.983799950150d-12)
     +            *q-5.400548574357d-11)*q+6.638043374114d-10)
     +            *q-8.799518866802d-09)*q+1.791418482685d-07)
     +            *q-2.96075397351101d-06)*q+3.38028206156144d-05)
     +            *q-3.58426847857878d-04)*q + 8.39213709428516d-03
            ww4 = (((((((((((1.33829971060180d-17*q-3.44841877844140d-16
     +            )*q+4.745009557656d-15)*q-6.033814209875d-14)
     +            *q+1.049256040808d-12)*q-1.70859789556117d-11)
     +            *q+2.15219425727959d-10)*q-2.52746574206884d-09)
     +            *q+3.27761714422960d-08)*q-3.90387662925193d-07)
     +            *q+3.46340204593870d-06)*q-2.43236345136782d-05)
     +            *q + 3.54846978585226d-04
            ww5 = (((((((((((((2.69412277020887d-20*q-
     +            4.24837886165685d-19)*q+6.030500065438d-18)
     +            *q-9.069722758289d-17)*q+1.246599177672d-15)
     +            *q-1.56872999797549d-14)*q+1.87305099552692d-13)
     +            *q-2.09498886675861d-12)*q+2.11630022068394d-11)
     +            *q-1.92566242323525d-10)*q+1.62012436344069d-09)
     +            *q-1.23621614171556d-08)*q+7.72165684563049d-08)
     +            *q-3.59858901591047d-07)*q + 2.43682618601000d-06
            return
         end if
      else if (pp.gt.5.0d+00) then
         if (pp.gt.10.0d+00) then
c     x=10.0 to 15.0                             nroots = 5
            q = pp - 12.5d+00
            rt1 = ((((((((((-4.16387977337393d-17*q+7.20872997373860d-16
     +            )*q+1.395993802064d-14)*q+3.660484641252d-14)
     +            *q-4.154857548139d-12)*q+2.301379846544d-11)
     +            *q-1.033307012866d-09)*q+3.997777641049d-08)
     +            *q-9.35118186333939d-07)*q+2.38589932752937d-05)
     +            *q-5.35185183652937d-04)*q + 8.85218988709735d-03
            rt2 = ((((((((((-4.56279214732217d-16*q+6.24941647247927d-15
     +            )*q+1.737896339191d-13)*q+8.964205979517d-14)
     +            *q-3.538906780633d-11)*q+9.561341254948d-11)
     +            *q-9.772831891310d-09)*q+4.240340194620d-07)
     +            *q-1.02384302866534d-05)*q+2.57987709704822d-04)
     +            *q-5.54735977651677d-03)*q + 8.68245143991948d-02
            rt3 = ((((((((((-2.52879337929239d-15*q+2.13925810087833d-14
     +            )*q+7.884307667104d-13)*q-9.023398159510d-13)
     +            *q-5.814101544957d-11)*q-1.333480437968d-09)
     +            *q-2.217064940373d-08)*q+1.643290788086d-06)
     +            *q-4.39602147345028d-05)*q+1.08648982748911d-03)
     +            *q-2.13014521653498d-02)*q + 2.94150684465425d-01
            rt4 = ((((((((((-6.42391438038888d-15*q+5.37848223438815d-15
     +            )*q+8.960828117859d-13)*q+5.214153461337d-11)
     +            *q-1.106601744067d-10)*q-2.007890743962d-08)
     +            *q+1.543764346501d-07)*q+4.520749076914d-06)
     +            *q-1.88893338587047d-04)*q+4.73264487389288d-03)
     +            *q-7.91197893350253d-02)*q + 8.60057928514554d-01
            rt5 = (((((((((((-2.24366166957225d-14*q+
     +            4.87224967526081d-14)*q+5.587369053655d-12)
     +            *q-3.045253104617d-12)*q-1.223983883080d-09)
     +            *q-2.05603889396319d-09)*q+2.58604071603561d-07)
     +            *q+1.34240904266268d-06)*q-5.72877569731162d-05)
     +            *q-9.56275105032191d-04)*q+4.23367010370921d-02)
     +            *q-5.76800927133412d-01)*q + 3.87328263873381d+00
            ww1 = (((((((((8.98007931950169d-15*q+7.25673623859497d-14)*
     +            q+5.851494250405d-14)*q-4.234204823846d-11)
     +            *q+3.911507312679d-10)*q-9.65094802088511d-09)
     +            *q+3.42197444235714d-07)*q-7.51821178144509d-06)
     +            *q+1.94218051498662d-04)*q-5.38533819142287d-03)
     +            *q + 1.68122596736809d-01
            ww2 = ((((((((((-1.05490525395105d-15*q+1.96855386549388d-14
     +            )*q-5.500330153548d-13)*q+1.003849567976d-11)
     +            *q-1.720997242621d-10)*q+3.533277061402d-09)
     +            *q-6.389171736029d-08)*q+1.046236652393d-06)
     +            *q-1.73148206795827d-05)*q+2.57820531617185d-04)
     +            *q-3.46188265338350d-03)*q + 7.03302497508176d-02
            ww3 = (((((((((((3.60020423754545d-16*q-6.24245825017148d-15
     +            )*q+9.945311467434d-14)*q-1.749051512721d-12)
     +            *q+2.768503957853d-11)*q-4.08688551136506d-10)
     +            *q+6.04189063303610d-09)*q-8.23540111024147d-08)
     +            *q+1.01503783870262d-06)*q-1.20490761741576d-05)
     +            *q+1.26928442448148d-04)*q-1.05539461930597d-03)
     +            *q + 1.15543698537013d-02
            ww4 = (((((((((((((2.51163533058925d-18*q-
     +            4.31723745510697d-17)*q+6.557620865832d-16)
     +            *q-1.016528519495d-14)*q+1.491302084832d-13)
     +            *q-2.06638666222265d-12)*q+2.67958697789258d-11)
     +            *q-3.23322654638336d-10)*q+3.63722952167779d-09)
     +            *q-3.75484943783021d-08)*q+3.49164261987184d-07)
     +            *q-2.92658670674908d-06)*q+2.12937256719543d-05)
     +            *q-1.19434130620929d-04)*q + 6.45524336158384d-04
            ww5 = ((((((((((((((-1.29043630202811d-19*q+
     +            2.16234952241296d-18)*q-3.107631557965d-17)
     +            *q+4.570804313173d-16)*q-6.301348858104d-15)
     +            *q+8.031304476153d-14)*q-9.446196472547d-13)
     +            *q+1.018245804339d-11)*q-9.96995451348129d-11)
     +            *q+8.77489010276305d-10)*q-6.84655877575364d-09)
     +            *q+4.64460857084983d-08)*q-2.66924538268397d-07)
     +            *q+1.24621276265907d-06)*q-4.30868944351523d-06)
     +            *q + 9.94307982432868d-06
            return
         else
c     x=5.0 to 10.0                              nroots = 5
            q = pp - 7.5d+00
            rt1 = ((((((((-1.13825201010775d-14*q+1.89737681670375d-13)*
     +            q-4.81561201185876d-12)*q+1.56666512163407d-10)
     +            *q-3.73782213255083d-09)*q+9.15858355075147d-08)
     +            *q-2.13775073585629d-06)*q+4.56547356365536d-05)
     +            *q-8.68003909323740d-04)*q + 1.22703754069176d-02
            rt2 = (((((((((-3.67160504428358d-15*q+1.27876280158297d-14)
     +            *q-1.296476623788d-12)*q+1.477175434354d-11)
     +            *q+5.464102147892d-10)*q-2.42538340602723d-08)
     +            *q+8.20460740637617d-07)*q-2.20379304598661d-05)
     +            *q+4.90295372978785d-04)*q-9.14294111576119d-03)
     +            *q + 1.22590403403690d-01
            rt3 = (((((((((1.39017367502123d-14*q-6.96391385426890d-13)*
     +            q+1.176946020731d-12)*q+1.725627235645d-10)
     +            *q-3.686383856300d-09)*q+2.87495324207095d-08)
     +            *q+1.71307311000282d-06)*q-7.94273603184629d-05)
     +            *q+2.00938064965897d-03)*q-3.63329491677178d-02)
     +            *q + 4.34393683888443d-01
            rt4 = ((((((((((-1.27815158195209d-14*q+1.99910415869821d-14
     +            )*q+3.753542914426d-12)*q-2.708018219579d-11)
     +            *q-1.190574776587d-09)*q+1.106696436509d-08)
     +            *q+3.954955671326d-07)*q-4.398596059588d-06)
     +            *q-2.01087998907735d-04)*q+7.89092425542937d-03)
     +            *q-1.42056749162695d-01)*q + 1.39964149420683d+00
            rt5 = ((((((((((-1.19442341030461d-13*q-2.34074833275956d-12
     +            )*q+6.861649627426d-12)*q+6.082671496226d-10)
     +            *q+5.381160105420d-09)*q-6.253297138700d-08)
     +            *q-2.135966835050d-06)*q-2.373394341886d-05)
     +            *q+2.88711171412814d-06)*q+4.85221195290753d-02)
     +            *q-1.04346091985269d+00)*q + 7.89901551676692d+00
            ww1 = (((((((((7.95526040108997d-15*q-2.48593096128045d-13)*
     +            q+4.761246208720d-12)*q-9.535763686605d-11)
     +            *q+2.225273630974d-09)*q-4.49796778054865d-08)
     +            *q+9.17812870287386d-07)*q-1.86764236490502d-05)
     +            *q+3.76807779068053d-04)*q-8.10456360143408d-03)
     +            *q + 2.01097936411496d-01
            ww2 = (((((((((((1.25678686624734d-15*q-2.34266248891173d-14
     +            )*q+3.973252415832d-13)*q-6.830539401049d-12)
     +            *q+1.140771033372d-10)*q-1.82546185762009d-09)
     +            *q+2.77209637550134d-08)*q-4.01726946190383d-07)
     +            *q+5.48227244014763d-06)*q-6.95676245982121d-05)
     +            *q+8.05193921815776d-04)*q-8.15528438784469d-03)
     +            *q + 9.71769901268114d-02
            ww3 = ((((((((((((-8.20929494859896d-16*q+
     +            1.37356038393016d-14)*q-2.022863065220d-13)
     +            *q+3.058055403795d-12)*q-4.387890955243d-11)
     +            *q+5.923946274445d-10)*q-7.503659964159d-09)
     +            *q+8.851599803902d-08)*q-9.65561998415038d-07)
     +            *q+9.60884622778092d-06)*q-8.56551787594404d-05)
     +            *q+6.66057194311179d-04)*q-4.17753183902198d-03)
     +            *q + 2.25443826852447d-02
            ww4 = ((((((((((((((-1.08764612488790d-17*q+
     +            1.85299909689937d-16)*q-2.730195628655d-15)
     +            *q+4.127368817265d-14)*q-5.881379088074d-13)
     +            *q+7.805245193391d-12)*q-9.632707991704d-11)
     +            *q+1.099047050624d-09)*q-1.15042731790748d-08)
     +            *q+1.09415155268932d-07)*q-9.33687124875935d-07)
     +            *q+7.02338477986218d-06)*q-4.53759748787756d-05)
     +            *q+2.41722511389146d-04)*q-9.75935943447037d-04)
     +            *q + 2.57520532789644d-03
            ww5 = (((((((((((((((7.28996979748849d-19*q-
     +            1.26518146195173d-17)*q+1.886145834486d-16)
     +            *q-2.876728287383d-15)*q+4.114588668138d-14)
     +            *q-5.44436631413933d-13)*q+6.64976446790959d-12)
     +            *q-7.44560069974940d-11)*q+7.57553198166848d-10)
     +            *q-6.92956101109829d-09)*q+5.62222859033624d-08)
     +            *q-3.97500114084351d-07)*q+2.39039126138140d-06)
     +            *q-1.18023950002105d-05)*q+4.52254031046244d-05)
     +            *q-1.21113782150370d-04)*q + 1.75013126731224d-04
            return
         end if
      else if (pp.gt.1.0d+00) then
c     x=1.0 to 5.0                               nroots = 5
         q = pp - 3.0d+00
         rt1 = ((((((((-2.58163897135138d-14*q+8.14127461488273d-13)*q-
     +         2.11414838976129d-11)*q+5.09822003260014d-10)
     +         *q-1.16002134438663d-08)*q+2.46810694414540d-07)
     +         *q-4.92556826124502d-06)*q+9.02580687971053d-05)
     +         *q-1.45190025120726d-03)*q + 1.73416786387475d-02
         rt2 = (((((((((1.04525287289788d-14*q+5.44611782010773d-14)*q-
     +         4.831059411392d-12)*q+1.136643908832d-10)
     +         *q-1.104373076913d-09)*q-2.35346740649916d-08)
     +         *q+1.43772622028764d-06)*q-4.23405023015273d-05)
     +         *q+9.12034574793379d-04)*q-1.52479441718739d-02)
     +         *q + 1.76055265928744d-01
         rt3 = (((((((((-6.89693150857911d-14*q+5.92064260918861d-13)*q+
     +         1.847170956043d-11)*q-3.390752744265d-10)
     +         *q-2.995532064116d-09)*q+1.57456141058535d-07)
     +         *q-3.95859409711346d-07)*q-9.58924580919747d-05)
     +         *q+3.23551502557785d-03)*q-5.97587007636479d-02)
     +         *q + 6.46432853383057d-01
         rt4 = ((((((((-3.61293809667763d-12*q-2.70803518291085d-11)*q+
     +         8.83758848468769d-10)*q+1.59166632851267d-08)
     +         *q-1.32581997983422d-07)*q-7.60223407443995d-06)
     +         *q-7.41019244900952d-05)*q+9.81432631743423d-03)
     +         *q-2.23055570487771d-01)*q + 2.21460798080643d+00
         rt5 = (((((((((7.12332088345321d-13*q+3.16578501501894d-12)*q-
     +         8.776668218053d-11)*q-2.342817613343d-09)
     +         *q-3.496962018025d-08)*q-3.03172870136802d-07)
     +         *q+1.50511293969805d-06)*q+1.37704919387696d-04)
     +         *q+4.70723869619745d-02)*q-1.47486623003693d+00)
     +         *q + 1.35704792175847d+01
         ww1 = (((((((((1.04348658616398d-13*q-1.94147461891055d-12)*q+
     +         3.485512360993d-11)*q-6.277497362235d-10)
     +         *q+1.100758247388d-08)*q-1.88329804969573d-07)
     +         *q+3.12338120839468d-06)*q-5.04404167403568d-05)
     +         *q+8.00338056610995d-04)*q-1.30892406559521d-02)
     +         *q + 2.47383140241103d-01
         ww2 = (((((((((((3.23496149760478d-14*q-5.24314473469311d-13)*q
     +         +7.743219385056d-12)*q-1.146022750992d-10)
     +         *q+1.615238462197d-09)*q-2.15479017572233d-08)
     +         *q+2.70933462557631d-07)*q-3.18750295288531d-06)
     +         *q+3.47425221210099d-05)*q-3.45558237388223d-04)
     +         *q+3.05779768191621d-03)*q-2.29118251223003d-02)
     +         *q + 1.59834227924213d-01
         ww3 = ((((((((((((-3.42790561802876d-14*q+5.26475736681542d-13)
     +         *q-7.184330797139d-12)*q+9.763932908544d-11)
     +         *q-1.244014559219d-09)*q+1.472744068942d-08)
     +         *q-1.611749975234d-07)*q+1.616487851917d-06)
     +         *q-1.46852359124154d-05)*q+1.18900349101069d-04)
     +         *q-8.37562373221756d-04)*q+4.93752683045845d-03)
     +         *q-2.25514728915673d-02)*q + 6.95211812453929d-02
         ww4 = (((((((((((((1.04072340345039d-14*q-1.60808044529211d-13)
     +         *q+2.183534866798d-12)*q-2.939403008391d-11)
     +         *q+3.679254029085d-10)*q-4.23775673047899d-09)
     +         *q+4.46559231067006d-08)*q-4.26488836563267d-07)
     +         *q+3.64721335274973d-06)*q-2.74868382777722d-05)
     +         *q+1.78586118867488d-04)*q-9.68428981886534d-04)
     +         *q+4.16002324339929d-03)*q-1.28290192663141d-02)
     +         *q + 2.22353727685016d-02
         ww5 = ((((((((((((((-8.16770412525963d-16*q+
     +         1.31376515047977d-14)*q-1.856950818865d-13)
     +         *q+2.596836515749d-12)*q-3.372639523006d-11)
     +         *q+4.025371849467d-10)*q-4.389453269417d-09)
     +         *q+4.332753856271d-08)*q-3.82673275931962d-07)
     +         *q+2.98006900751543d-06)*q-2.00718990300052d-05)
     +         *q+1.13876001386361d-04)*q-5.23627942443563d-04)
     +         *q+1.83524565118203d-03)*q-4.37785737450783d-03)
     +         *q + 5.36963805223095d-03
         return
      else if (pp.gt.3.0d-07) then
c     x=0.0 to 1.0                               nroots = 5
        rt1 = ((((((-4.46679165328413d-11*pp+1.21879111988031d-09)*pp-
     +         2.62975022612104d-08)*pp+5.15106194905897d-07)
     +         *pp-9.27933625824749d-06)*pp+1.51794097682482d-04)
     +         *pp-2.15865967920301d-03)*pp + 2.26659266316985d-02
        rt2 = ((((((1.93117331714174d-10*pp-4.57267589660699d-09)*pp+
     +         2.48339908218932d-08)*pp+1.50716729438474d-06)
     +         *pp-6.07268757707381d-05)*pp+1.37506939145643d-03)
     +         *pp-2.20258754419939d-02)*pp + 2.31271692140905d-01
        rt3 = (((((4.84989776180094d-09*pp+1.31538893944284d-07)*pp-
     +         2.766753852879d-06)*pp-7.651163510626d-05)
     +         *pp+4.033058545972d-03)*pp-8.16520022916145d-02)
     +         *pp + 8.57346024118779d-01
        rt4 = ((((-2.48581772214623d-07*pp-4.34482635782585d-06)*pp-
     +         7.46018257987630d-07)*pp+1.01210776517279d-02)
     +         *pp-2.83193369640005d-01)*pp + 2.97353038120345d+00
        rt5 = (((((-8.92432153868554d-09*pp+1.77288899268988d-08)*pp+
     +         3.040754680666d-06)*pp+1.058229325071d-04)
     +         *pp+4.596379534985d-02)*pp-1.75382723579114d+00)
     +         *pp + 1.84151859759049d+01
        ww1 = ((((((-2.03822632771791d-09*pp+3.89110229133810d-08)*pp-
     +         5.84914787904823d-07)*pp+8.30316168666696d-06)
     +         *pp-1.13218402310546d-04)*pp+1.49128888586790d-03)
     +         *pp-1.96867576904816d-02)*pp + 2.95524224714749d-01
        ww2 = (((((((8.62848118397570d-09*pp-1.38975551148989d-07)*pp+
     +         1.602894068228d-06)*pp-1.646364300836d-05)
     +         *pp+1.538445806778d-04)*pp-1.28848868034502d-03)
     +         *pp+9.38866933338584d-03)*pp-5.61737590178812d-02)
     +         *pp + 2.69266719309991d-01
        ww3 = ((((((((-9.41953204205665d-09*pp+1.47452251067755d-07)*pp-
     +         1.57456991199322d-06)*pp+1.45098401798393d-05)
     +         *pp-1.18858834181513d-04)*pp+8.53697675984210d-04)
     +         *pp-5.22877807397165d-03)*pp+2.60854524809786d-02)
     +         *pp-9.71152726809059d-02)*pp + 2.19086362515979d-01
        ww4 = ((((((((-3.84961617022042d-08*pp+5.66595396544470d-07)*pp-
     +         5.52351805403748d-06)*pp+4.53160377546073d-05)
     +         *pp-3.22542784865557d-04)*pp+1.95682017370967d-03)
     +         *pp-9.77232537679229d-03)*pp+3.79455945268632d-02)
     +         *pp-1.02979262192227d-01)*pp + 1.49451349150573d-01
        ww5 = (((((((((4.09594812521430d-09*pp-6.47097874264417d-08)*pp+
     +         6.743541482689d-07)*pp-5.917993920224d-06)
     +         *pp+4.531969237381d-05)*pp-2.99102856679638d-04)
     +         *pp+1.65695765202643d-03)*pp-7.40671222520653d-03)
     +         *pp+2.50889946832192d-02)*pp-5.73782817487958d-02)
     +         *pp + 6.66713443086877d-02
         return
      else
c     x is approximately zero.                   nroots = 5
         rt1 = 2.26659266316985d-02 - 2.15865967920897d-03*pp
         rt2 = 2.31271692140903d-01 - 2.20258754389745d-02*pp
         rt3 = 8.57346024118836d-01 - 8.16520023025515d-02*pp
         rt4 = 2.97353038120346d+00 - 2.83193369647137d-01*pp
         rt5 = 1.84151859759051d+01 - 1.75382723579439d+00*pp
         ww1 = 2.95524224714752d-01 - 1.96867576909777d-02*pp
         ww2 = 2.69266719309995d-01 - 5.61737590184721d-02*pp
         ww3 = 2.19086362515981d-01 - 9.71152726793658d-02*pp
         ww4 = 1.49451349150580d-01 - 1.02979262193565d-01*pp
         ww5 = 6.66713443086877d-02 - 5.73782817488315d-02*pp
         return
      end if
      end
      subroutine rootss
      implicit real*8  (a-h,o-z)
      logical jump
      dimension f(60)
      common/rtdata/whi(45),wlow(45),rhi(45),rlow(45),rfac(45),
     1amps(9),ipoint(9),madd(20)
      dimension  beta(12),gamma(12),v(12),d(12)
c
c
      real*8 pp, u, w, dji
      integer nroots, nroot2
      common /root/ pp,u(12),w(12),nroots,nroot2,dji(60)
c
c
c      finds weights and points for rys integration
c      using legendre or hermite points as starting
c      estimates.
c
      data m1/1/
      data ms2,m24/-2,24/
      data eps/1.0d-12/
      data d(1)/1.0d0/
c
c
      nroot2 = nroots - 2
      mprime = nroots + nroots - 1
      mprim2 = mprime - 1
      iflag = ipoint(nroots)
      xmax = amps(nroots)
      if (pp.le.1.0d-11) then
         do 20 ipt = 1 , nroots
            u(ipt) = rlow(iflag+ipt)
            w(ipt) = wlow(iflag+ipt)
 20      continue
      else if (pp.ge.xmax) then
         top = 1.0d0/pp
         scale = dsqrt(top)
         do 30 ipt = 1 , nroots
            u(ipt) = rhi(iflag+ipt)*top
            w(ipt) = whi(iflag+ipt)*scale
 30      continue
      else
         call ffun(pp,f,nroots,dji,madd)
         a = f(2)/f(1)
         if (nroots.gt.m1) then
            a1 = a
            do 40 mm = 1 , mprim2
               f(mm+24) = f(mm+2) - a*f(mm+1)
 40         continue
            ibas = 0
            jbas = m24
            max = mprim2
            do 60 k = 2 , nroots
               cfac = f(jbas+1)
               ax = f(jbas+2)/cfac
               bet = a - ax
               beta(k-1) = bet
               a = ax
               gam = -cfac/f(ibas+1)
               gamma(k-1) = gam
               if (k.eq.nroots) go to 70
               max = max + ms2
               do 50 mm = 1 , max
                  f(ibas+mm) = bet*f(jbas+mm+1) + gam*f(ibas+mm+2)
     +                         + f(jbas+mm+2)
 50            continue
               ibas = jbas
               jbas = m24 - jbas
 60         continue
 70         gam = 0.0d0
            do 80 k = 1 , nroots
               r = rlow(iflag+k)/(rfac(iflag+k)*pp+1.0d0)
               if (pp.gt.17.0d0) r = rhi(iflag+k)/pp
               u(k) = r
               gam = gam + r
 80         continue
            epsi = eps*ax
            do 120 mm = 1 , nroots
               rrr = u(mm)
               r = rrr*ax/gam
               jump = .false.
               do 100 icount = 1 , 50
                  v(1) = r - a1
                  top = r + beta(1)
                  dd = v(1) + top
                  d(2) = dd
                  vv = v(1)*top + gamma(1)
                  if (nroot2.gt.0) then
                     v(2) = vv
                     do 90 k = 3 , nroots
                        top = r + beta(k-1)
                        dd = top*dd + gamma(k-1)*d(k-2) + vv
                        d(k) = dd
                        vv = top*vv + gamma(k-1)*v(k-2)
                        v(k) = vv
 90                  continue
                  end if
                  scale = vv/dd
                  r = r - scale
                  if (jump) go to 110
                  if (dabs(scale).lt.epsi) jump = .true.
 100           continue
               write (6,6010) pp , nroots
 110           u(mm) = r
               w(mm) = cfac/(dd*v(nroots-1))
               ax = ax - r
               gam = gam - rrr
 120        continue
         else
            u(1) = a
            w(1) = f(1)
         end if
      end if
      do 130 i = 1 , nroots
         rr = u(i)
         u(i) = rr/(1.0d0-rr)
 130  continue
      return
 6010 format (1x,'no convergence',f15.6,i6)
      end
      subroutine rt123
c             *****   version february 13,1975   *****
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      real*8 pp, u, w, dji
      integer nroots, nroot2
      common /root/ pp,u(12),w(12),nroots,nroot2,dji(60)
c
      equivalence (u(1),rt1),(u(2),rt2),(u(3),rt3),(u(4),rt4),(u(5),rt5)
      equivalence (w(1),ww1),(w(2),ww2),(w(3),ww3),(w(4),ww4),(w(5),ww5)
      data r12,pie4/2.75255128608411d-01, 7.85398163397448d-01/
      data r22,w22/ 2.72474487139158d+00, 9.17517095361369d-02/
      data r13/     1.90163509193487d-01/
      data r23,w23/ 1.78449274854325d+00, 1.77231492083829d-01/
      data r33,w33/ 5.52534374226326d+00, 5.11156880411248d-03/
      if (pp.gt.5.0d+00) then
         if (pp.le.15.0d+00) then
            g = dexp(-pp)
            if (pp.gt.10.0d+00) then
c     x = 10.0 to 15.0                 nroots=1,2, or 3
               ww1 = (((-1.8784686463512d-01/pp+2.2991849164985d-01)
     +               /pp-4.9893752514047d-01)/pp-2.1916512131607d-05)
     +               *g + dsqrt(pie4/pp)
               f1 = (ww1-g)/(pp+pp)
               if (nroots.lt.2) then
                  rt1 = f1/(ww1-f1)
                  return
               else if (nroots.eq.2) then
                 rt1 = ((((-1.01041157064226d-05*pp+1.19483054115173d-03
     +                  )*pp-6.73760231824074d-02)
     +                  *pp+1.25705571069895d+00)
     +                  *pp+(((-8.57609422987199d+03/pp+
     +                  5.91005939591842d+03)/pp-1.70807677109425d+03)
     +                /pp+2.64536689959503d+02)/pp-2.38570496490846d+01)
     +                  *g + r12/(pp-r12)
                 rt2 = (((3.39024225137123d-04*pp-9.34976436343509d-02)
     +                  *pp-4.22216483306320d+00)
     +                  *pp+(((-2.08457050986847d+03/pp-
     +                  1.04999071905664d+03)/pp+3.39891508992661d+02)
     +                /pp-1.56184800325063d+02)/pp+8.00839033297501d+00)
     +                  *g + r22/(pp-r22)
                 ww2 = ((f1-ww1)*rt1+f1)*(1.0d+00+rt2)/(rt2-rt1)
                  ww1 = ww1 - ww2
                  return
               else
                  f2 = (f1+f1+f1-g)/(pp+pp)
                  q = pp - 12.5d+00
                  rt1 = (((((((((((4.42133001283090d-16*q-
     +                  2.77189767070441d-15)*q-4.084026087887d-14)
     +                  *q+5.379885121517d-13)*q+1.882093066702d-12)
     +                  *q-8.67286219861085d-11)*q+7.11372337079797d-10)
     +                  *q-3.55578027040563d-09)*q+1.29454702851936d-07)
     +                  *q-4.14222202791434d-06)*q+8.04427643593792d-05)
     +                  *q-1.18587782909876d-03)
     +                  *q + 1.53435577063174d-02
                  rt2 = (((((((((((6.85146742119357d-15*q-
     +                  1.08257654410279d-14)*q-8.579165965128d-13)
     +                  *q+6.642452485783d-12)*q+4.798806828724d-11)
     +                  *q-1.13413908163831d-09)*q+7.08558457182751d-09)
     +                  *q-5.59678576054633d-08)*q+2.51020389884249d-06)
     +                  *q-6.63678914608681d-05)*q+1.11888323089714d-03)
     +                  *q-1.45361636398178d-02)
     +                  *q + 1.65077877454402d-01
                  rt3 = ((((((((((((3.20622388697743d-15*q-
     +                  2.73458804864628d-14)*q-3.157134329361d-13)
     +                  *q+8.654129268056d-12)*q-5.625235879301d-11)
     +                  *q-7.718080513708d-10)*q+2.064664199164d-08)
     +                  *q-1.567725007761d-07)*q-1.57938204115055d-06)
     +                  *q+6.27436306915967d-05)*q-1.01308723606946d-03)
     +                  *q+1.13901881430697d-02)*q-1.01449652899450d-01)
     +                  *q + 7.77203937334739d-01
                  go to 20
               end if
            else
c     x = 5.0 to 10.0                  nroots =1,2, or 3
             ww1 = ((((((4.6897511375022d-01/pp-6.9955602298985d-01)/pp+
     +               5.3689283271887d-01)/pp-3.2883030418398d-01)
     +               /pp+2.4645596956002d-01)/pp-4.9984072848436d-01)
     +               /pp-3.1501078774085d-06)*g + dsqrt(pie4/pp)
               f1 = (ww1-g)/(pp+pp)
               if (nroots.lt.2) then
                  rt1 = f1/(ww1-f1)
                  return
               else if (nroots.eq.2) then
                  q = pp - 7.5d+00
                  rt1 = (((((((((((((-1.43632730148572d-16*q+
     +                  2.38198922570405d-16)*q+1.358319618800d-14)
     +                  *q-7.064522786879d-14)*q-7.719300212748d-13)
     +                  *q+7.802544789997d-12)*q+6.628721099436d-11)
     +                  *q-1.775564159743d-09)*q+1.713828823990d-08)
     +                  *q-1.497500187053d-07)*q+2.283485114279d-06)
     +                  *q-3.76953869614706d-05)*q+4.74791204651451d-04)
     +                  *q-4.60448960876139d-03)
     +                  *q + 3.72458587837249d-02
                  rt2 = ((((((((((((2.48791622798900d-14*q-
     +                  1.36113510175724d-13)*q-2.224334349799d-12)
     +                  *q+4.190559455515d-11)*q-2.222722579924d-10)
     +                  *q-2.624183464275d-09)*q+6.128153450169d-08)
     +                  *q-4.383376014528d-07)*q-2.49952200232910d-06)
     +                  *q+1.03236647888320d-04)*q-1.44614664924989d-03)
     +                  *q+1.35094294917224d-02)*q-9.53478510453887d-02)
     +                  *q + 5.44765245686790d-01
                  ww2 = ((f1-ww1)*rt1+f1)*(1.0d+00+rt2)/(rt2-rt1)
                  ww1 = ww1 - ww2
                  return
               else
                  f2 = (f1+f1+f1-g)/(pp+pp)
                  q = pp - 7.5d+00
                  rt1 = (((((((((((5.74429401360115d-16*q+
     +                  7.11884203790984d-16)*q-6.736701449826d-14)
     +                  *q-6.264613873998d-13)*q+1.315418927040d-11)
     +                  *q-4.23879635610964d-11)*q+1.39032379769474d-09)
     +                  *q-4.65449552856856d-08)*q+7.34609900170759d-07)
     +                  *q-1.08656008854077d-05)*q+1.77930381549953d-04)
     +                  *q-2.39864911618015d-03)
     +                  *q + 2.39112249488821d-02
                  rt2 = (((((((((((1.13464096209120d-14*q+
     +                  6.99375313934242d-15)*q-8.595618132088d-13)
     +                  *q-5.293620408757d-12)*q-2.492175211635d-11)
     +                  *q+2.73681574882729d-09)*q-1.06656985608482d-08)
     +                  *q-4.40252529648056d-07)*q+9.68100917793911d-06)
     +                  *q-1.68211091755327d-04)*q+2.69443611274173d-03)
     +                  *q-3.23845035189063d-02)
     +                  *q + 2.75969447451882d-01
                  rt3 = ((((((((((((6.66339416996191d-15*q+
     +                  1.84955640200794d-13)*q-1.985141104444d-12)
     +                  *q-2.309293727603d-11)*q+3.917984522103d-10)
     +                  *q+1.663165279876d-09)*q-6.205591993923d-08)
     +                  *q+8.769581622041d-09)*q+8.97224398620038d-06)
     +                  *q-3.14232666170796d-05)*q-1.83917335649633d-03)
     +                  *q+3.51246831672571d-02)*q-3.22335051270860d-01)
     +                  *q + 1.73582831755430d+00
                  go to 20
               end if
            end if
         else if (pp.gt.33.0d+00) then
c     x = 33.0 to infinity             nroots=1,2, or 3
            ww1 = dsqrt(pie4/pp)
            if (nroots.lt.2) then
               rt1 = 0.5d+00/(pp-0.5d+00)
               return
            else if (nroots.eq.2) then
               if (pp.gt.40.0d+00) then
                  rt1 = r12/(pp-r12)
                  rt2 = r22/(pp-r22)
                  ww2 = w22*ww1
                  ww1 = ww1 - ww2
                  return
               else
                  g = dexp(-pp)
                  rt1 = (-8.78947307498880d-01*pp+1.09243702330261d+01)
     +                  *g + r12/(pp-r12)
                  rt2 = (-9.28903924275977d+00*pp+8.10642367843811d+01)
     +                  *g + r22/(pp-r22)
                  ww2 = (4.46857389308400d+00*pp-7.79250653461045d+01)
     +                  *g + w22*ww1
                  ww1 = ww1 - ww2
                  return
               end if
            else if (pp.gt.47.0d+00) then
               rt1 = r13/(pp-r13)
               rt2 = r23/(pp-r23)
               rt3 = r33/(pp-r33)
               ww2 = w23*ww1
               ww3 = w33*ww1
               ww1 = ww1 - ww2 - ww3
               return
            else
               g = dexp(-pp)
               rt1 = ((-7.39058467995275d+00*pp+3.21318352526305d+02)
     +               *pp-3.99433696473658d+03)*g + r13/(pp-r13)
               rt2 = ((-7.38726243906513d+01*pp+3.13569966333873d+03)
     +               *pp-3.86862867311321d+04)*g + r23/(pp-r23)
               rt3 = ((-2.63750565461336d+02*pp+1.04412168692352d+04)
     +               *pp-1.28094577915394d+05)*g + r33/(pp-r33)
               ww3 = (((1.52258947224714d-01*pp-8.30661900042651d+00)
     +               *pp+1.92977367967984d+02)*pp-1.67787926005344d+03)
     +               *g + w33*ww1
               ww2 = ((6.15072615497811d+01*pp-2.91980647450269d+03)
     +               *pp+3.80794303087338d+04)*g + w23*ww1
               ww1 = ww1 - ww2 - ww3
               return
            end if
         else
c     x = 15.0 to 33.0                 nroots=1,2, or 3
            g = dexp(-pp)
            ww1 = ((1.9623264149430d-01/pp-4.9695241464490d-01)
     +            /pp-6.0156581186481d-05)*g + dsqrt(pie4/pp)
            f1 = (ww1-g)/(pp+pp)
            if (nroots.lt.2) then
               rt1 = f1/(ww1-f1)
               return
            else if (nroots.eq.2) then
             rt1 = ((((-1.14906395546354d-06*pp+1.76003409708332d-04)*pp
     +               -1.71984023644904d-02)*pp-1.37292644149838d-01)
     +               *pp+(-4.75742064274859d+01/pp+9.21005186542857d+00)
     +               /pp-2.31080873898939d-02)*g + r12/(pp-r12)
             rt2 = (((3.64921633404158d-04*pp-9.71850973831558d-02)
     +               *pp-4.02886174850252d+00)
     +               *pp+(-1.35831002139173d+02/pp-8.66891724287962d+01)
     +               /pp+2.98011277766958d+00)*g + r22/(pp-r22)
               ww2 = ((f1-ww1)*rt1+f1)*(1.0d+00+rt2)/(rt2-rt1)
               ww1 = ww1 - ww2
               return
            else
               f2 = (f1+f1+f1-g)/(pp+pp)
               if (pp.gt.20.0d+00) then
                 rt1 = ((((-4.97561537069643d-04*pp-5.00929599665316d-02
     +                  )*pp+1.31099142238996d+00)
     +                  *pp-1.88336409225481d+01)
     +                 *pp-6.60344754467191d+02/pp+1.64931462413877d+02)
     +                  *g + r13/(pp-r13)
                 rt2 = ((((-4.48218898474906d-03*pp-5.17373211334924d-01
     +                  )*pp+1.13691058739678d+01)
     +                  *pp-1.65426392885291d+02)
     +                 *pp-6.30909125686731d+03/pp+1.52231757709236d+03)
     +                  *g + r23/(pp-r23)
                 rt3 = ((((-1.38368602394293d-02*pp-1.77293428863008d+00
     +                  )*pp+1.73639054044562d+01)
     +                  *pp-3.57615122086961d+02)
     +                 *pp-1.45734701095912d+04/pp+2.69831813951849d+03)
     +                  *g + r33/(pp-r33)
               else
                 rt1 = ((((((-2.43270989903742d-06*pp+
     +                  3.57901398988359d-04)*pp-2.34112415981143d-02)
     +                *pp+7.81425144913975d-01)*pp-1.73209218219175d+01)
     +                  *pp+2.43517435690398d+02)
     +                  *pp+(-1.97611541576986d+04/pp+
     +                  9.82441363463929d+03)/pp-2.07970687843258d+03)
     +                  *g + r13/(pp-r13)
                 rt2 = (((((-2.62627010965435d-04*pp+
     +                  3.49187925428138d-02)*pp-3.09337618731880d+00)
     +                *pp+1.07037141010778d+02)*pp-2.36659637247087d+03)
     +                  *pp+
     +                  ((-2.91669113681020d+06/pp+1.41129505262758d+06)
     +                /pp-2.91532335433779d+05)/pp+3.35202872835409d+04)
     +                  *g + r23/(pp-r23)
                 rt3 = (((((9.31856404738601d-05*pp-2.87029400759565d-02
     +                  )*pp-7.83503697918455d-01)
     +                *pp-1.84338896480695d+01)*pp+4.04996712650414d+02)
     +                  *pp+
     +                  (-1.89829509315154d+05/pp+5.11498390849158d+04)
     +                  /pp-6.88145821789955d+03)*g + r33/(pp-r33)
               end if
               go to 20
            end if
         end if
      else if (pp.gt.1.0d+00) then
         if (pp.gt.3.0d+00) then
c     x = 3.0 to 5.0                   nroots =1,2, or 3
            q = pp - 4.0d+00
            if (nroots.eq.3) then
               rt1 = (((((((1.44265709189601d-11*q-4.66622033006074d-10)
     +               *q+7.649155832025d-09)*q-1.229940017368d-07)
     +               *q+2.026002142457d-06)*q-2.87048671521677d-05)
     +               *q+3.70326938096287d-04)*q-4.21006346373634d-03)
     +               *q + 3.50898470729044d-02
               rt2 = ((((((((-2.65526039155651d-11*q+
     +               1.97549041402552d-10)*q+2.15971131403034d-09)
     +               *q-7.95045680685193d-08)*q+5.15021914287057d-07)
     +               *q+1.11788717230514d-05)*q-3.33739312603632d-04)
     +               *q+5.30601428208358d-03)*q-5.93483267268959d-02)
     +               *q + 4.31180523260239d-01
               rt3 = ((((((((-3.92833750584041d-10*q-
     +               4.16423229782280d-09)*q+4.42413039572867d-08)
     +               *q+6.40574545989551d-07)*q-3.05512456576552d-06)
     +               *q-1.05296443527943d-04)*q-6.14120969315617d-04)
     +               *q+4.89665802767005d-02)*q-6.24498381002855d-01)
     +               *q + 3.36412312243724d+00
               f2 = ((((((((((-2.36788772599074d-11*q+
     +              2.89147476459092d-10)*q-3.18111322308846d-09)
     +              *q+3.25336816562485d-08)*q-3.00873821471489d-07)
     +              *q+2.48749160874431d-06)*q-1.81353179793672d-05)
     +              *q+1.14504948737066d-04)*q-6.10614987696677d-04)
     +              *q+2.64584212770942d-03)*q-8.66415899015349d-03)
     +              *q + 1.75257821619922d-02
            else
               f1 = ((((((((((-2.62453564772299d-11*q+
     +              3.24031041623823d-10)*q-3.614965656163d-09)
     +              *q+3.760256799971d-08)*q-3.553558319675d-07)
     +              *q+3.022556449731d-06)*q-2.290098979647d-05)
     +              *q+1.526537461148d-04)*q-8.81947375894379d-04)
     +              *q+4.33207949514611d-03)*q-1.75257821619926d-02)
     +              *q + 5.28406320615584d-02
               ww1 = (pp+pp)*f1 + dexp(-pp)
               if (nroots.eq.2) then
                  rt1 = ((((((((-4.11560117487296d-12*q+
     +                  7.10910223886747d-11)*q-1.73508862390291d-09)
     +                  *q+5.93066856324744d-08)*q-9.76085576741771d-07)
     +                  *q+1.08484384385679d-05)*q-1.12608004981982d-04)
     +                  *q+1.16210907653515d-03)*q-9.89572595720351d-03)
     +                  *q + 6.12589701086408d-02
                  rt2 = (((((((((-1.80555625241001d-10*q+
     +                  5.44072475994123d-10)*q+1.603498045240d-08)
     +                  *q-1.497986283037d-07)*q-7.017002532106d-07)
     +                  *q+1.85882653064034d-05)*q-2.04685420150802d-05)
     +                  *q-2.49327728643089d-03)*q+3.56550690684281d-02)
     +                  *q-2.60417417692375d-01)
     +                  *q + 1.12155283108289d+00
                  ww2 = ((f1-ww1)*rt1+f1)*(1.0d+00+rt2)/(rt2-rt1)
                  ww1 = ww1 - ww2
                  return
               else
                  rt1 = f1/(ww1-f1)
                  return
               end if
            end if
         else
c     x = 1.0 to 3.0                   nroots=1,2, or 3
            q = pp - 2.0d+00
            if (nroots.eq.3) then
               rt1 = ((((((((1.44687969563318d-12*q+4.85300143926755d-12
     +               )*q-6.55098264095516d-10)*q+1.56592951656828d-08)
     +               *q-2.60122498274734d-07)*q+3.86118485517386d-06)
     +               *q-5.13430986707889d-05)*q+6.03194524398109d-04)
     +               *q-6.11219349825090d-03)*q + 4.52578254679079d-02
               rt2 = (((((((6.95964248788138d-10*q-5.35281831445517d-09)
     +               *q-6.745205954533d-08)*q+1.502366784525d-06)
     +               *q+9.923326947376d-07)*q-3.89147469249594d-04)
     +               *q+7.51549330892401d-03)*q-8.48778120363400d-02)
     +               *q + 5.73928229597613d-01
               rt3 = ((((((((-2.81496588401439d-10*q+
     +               3.61058041895031d-09)*q+4.53631789436255d-08)
     +               *q-1.40971837780847d-07)*q-6.05865557561067d-06)
     +               *q-5.15964042227127d-05)*q+3.34761560498171d-05)
     +               *q+5.04871005319119d-02)*q-8.24708946991557d-01)
     +               *q + 4.81234667357205d+00
               f2 = ((((((((((-1.48044231072140d-10*q+
     +              1.78157031325097d-09)*q-1.92514145088973d-08)
     +              *q+1.92804632038796d-07)*q-1.73806555021045d-06)
     +              *q+1.39195169625425d-05)*q-9.74574633246452d-05)
     +              *q+5.83701488646511d-04)*q-2.89955494844975d-03)
     +              *q+1.13847001113810d-02)*q-3.23446977320647d-02)
     +              *q + 5.29428148329709d-02
            else
               f1 = ((((((((((-1.61702782425558d-10*q+
     +              1.96215250865776d-09)*q-2.14234468198419d-08)
     +              *q+2.17216556336318d-07)*q-1.98850171329371d-06)
     +              *q+1.62429321438911d-05)*q-1.16740298039895d-04)
     +              *q+7.24888732052332d-04)*q-3.79490003707156d-03)
     +              *q+1.61723488664661d-02)*q-5.29428148329736d-02)
     +              *q + 1.15702180856167d-01
               ww1 = (pp+pp)*f1 + dexp(-pp)
               if (nroots.eq.2) then
                  rt1 = (((((((((-6.36859636616415d-12*q+
     +                  8.47417064776270d-11)*q-5.152207846962d-10)
     +                  *q-3.846389873308d-10)*q+8.472253388380d-08)
     +                  *q-1.85306035634293d-06)*q+2.47191693238413d-05)
     +                  *q-2.49018321709815d-04)*q+2.19173220020161d-03)
     +                  *q-1.63329339286794d-02)
     +                  *q + 8.68085688285261d-02
                  rt2 = (((((((((1.45331350488343d-10*q+
     +                  2.07111465297976d-09)*q-1.878920917404d-08)
     +                  *q-1.725838516261d-07)*q+2.247389642339d-06)
     +                  *q+9.76783813082564d-06)*q-1.93160765581969d-04)
     +                  *q-1.58064140671893d-03)*q+4.85928174507904d-02)
     +                  *q-4.30761584997596d-01)
     +                  *q + 1.80400974537950d+00
                  ww2 = ((f1-ww1)*rt1+f1)*(1.0d+00+rt2)/(rt2-rt1)
                  ww1 = ww1 - ww2
                  return
               else
                  rt1 = f1/(ww1-f1)
                  return
               end if
            end if
         end if
      else if (pp.gt.3.0d-07) then
c     x = 0.0 to 1.0                   nroots=1,2, or 3
         if (nroots.eq.3) then
          rt1 = ((((((-5.10186691538870d-10*pp+2.40134415703450d-08)*pp-
     +            5.01081057744427d-07)*pp+7.58291285499256d-06)
     +            *pp-9.55085533670919d-05)*pp+1.02893039315878d-03)
     +            *pp-9.28875764374337d-03)*pp + 6.03769246832810d-02
          rt2 = ((((((-1.29646524960555d-08*pp+7.74602292865683d-08)*pp+
     +            1.56022811158727d-06)*pp-1.58051990661661d-05)
     +            *pp-3.30447806384059d-04)*pp+9.74266885190267d-03)
     +            *pp-1.19511285526388d-01)*pp + 7.76823355931033d-01
          rt3 = ((((((-9.28536484109606d-09*pp-3.02786290067014d-07)*pp-
     +            2.50734477064200d-06)*pp-7.32728109752881d-06)
     +            *pp+2.44217481700129d-04)*pp+4.94758452357327d-02)
     +            *pp-1.02504611065774d+00)*pp + 6.66279971938553d+00
          f2 = ((((((((-7.60911486098850d-08*pp+1.09552870123182d-06)*pp
     +           -1.03463270693454d-05)*pp+8.16324851790106d-05)
     +           *pp-5.55526624875562d-04)*pp+3.20512054753924d-03)
     +           *pp-1.51515139838540d-02)*pp+5.55555554649585d-02)
     +           *pp-1.42857142854412d-01)*pp + 1.99999999999986d-01
         else
          f1 = ((((((((-8.36313918003957d-08*pp+1.21222603512827d-06)*pp
     +           -1.15662609053481d-05)*pp+9.25197374512647d-05)
     +           *pp-6.40994113129432d-04)*pp+3.78787044215009d-03)
     +           *pp-1.85185172458485d-02)*pp+7.14285713298222d-02)
     +           *pp-1.99999999997023d-01)*pp + 3.33333333333318d-01
            ww1 = (pp+pp)*f1 + dexp(-pp)
            if (nroots.eq.2) then
              rt1 = (((((((-2.35234358048491d-09*pp+2.49173650389842d-08
     +               )*pp-4.558315364581d-08)*pp-2.447252174587d-06)
     +               *pp+4.743292959463d-05)*pp-5.33184749432408d-04)
     +               *pp+4.44654947116579d-03)*pp-2.90430236084697d-02)
     +               *pp + 1.30693606237085d-01
              rt2 = (((((((-2.47404902329170d-08*pp+2.36809910635906d-07
     +               )*pp+1.835367736310d-06)*pp-2.066168802076d-05)
     +               *pp-1.345693393936d-04)*pp-5.88154362858038d-05)
     +               *pp+5.32735082098139d-02)*pp-6.37623643056745d-01)
     +               *pp + 2.86930639376289d+00
               ww2 = ((f1-ww1)*rt1+f1)*(1.0d+00+rt2)/(rt2-rt1)
               ww1 = ww1 - ww2
               return
            else
               rt1 = f1/(ww1-f1)
               return
            end if
         end if
c     x is approximately zero.         nroots=1,2, or 3
      else if (nroots.lt.2) then
         rt1 = 0.5d+00 - pp/5.0d+00
         ww1 = 1.0d+00 - pp/3.0d+00
         return
      else if (nroots.eq.2) then
         rt1 = 1.30693606237085d-01 - 2.90430236082028d-02*pp
         rt2 = 2.86930639376291d+00 - 6.37623643058102d-01*pp
         ww1 = 6.52145154862545d-01 - 1.22713621927067d-01*pp
         ww2 = 3.47854845137453d-01 - 2.10619711404725d-01*pp
         return
      else
         rt1 = 6.03769246832797d-02 - 9.28875764357368d-03*pp
         rt2 = 7.76823355931043d-01 - 1.19511285527878d-01*pp
         rt3 = 6.66279971938567d+00 - 1.02504611068957d+00*pp
         ww1 = 4.67913934572691d-01 - 5.64876917232519d-02*pp
         ww2 = 3.60761573048137d-01 - 1.49077186455208d-01*pp
         ww3 = 1.71324492379169d-01 - 1.27768455150979d-01*pp
         return
      end if
      g = dexp(-pp)
      f1 = ((pp+pp)*f2+g)/3.0d+00
      ww1 = (pp+pp)*f1 + g
 20   t1 = rt1/(rt1+1.0d+00)
      t2 = rt2/(rt2+1.0d+00)
      t3 = rt3/(rt3+1.0d+00)
      a2 = f2 - t1*f1
      a1 = f1 - t1*ww1
      ww3 = (a2-t2*a1)/((t3-t2)*(t3-t1))
      ww2 = (t3*a1-a2)/((t3-t2)*(t2-t1))
      ww1 = ww1 - ww2 - ww3
      return
      end
      subroutine sec192(s,f,t,lenbas)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
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
      dimension s(lenbas),t(lenbas),f(lenbas)
c
      common/blkin/potnuc,dx,dy,dz,nato,numo
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
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
c
      data m2/2/
c
      lenb = lensec(lenbas)
      i = lenb*6 + 1
      call secput(ionsec,m2,i,ib)
      nato = nat
      numo = num
c     header block 
      call wrt3(potnuc, 511, ib, idaf)
c     S, T, T+V
      ibl3s = ib + 1
      call wrt3s(s,lenbas,idaf)
      ibl3t = ibl3s + lenb
      call wrt3s(t,lenbas,idaf)
      ibl3f = ibl3t + lenb
      call wrt3s(f,lenbas,idaf)
c
      return
      end
      subroutine secdip(p,q,r,lenbas)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
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
      common/blkin/potnuc(4),nato,numo
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
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
c
      dimension p(lenbas),q(lenbas),r(lenbas)
      data m2/2/
c
      call secget(ionsec,m2,ib)
      call rdedx(potnuc,511,ib,idaf)
      do ij = 1,lenbas
      p(ij) = -p(ij)
      q(ij) = -q(ij)
      r(ij) = -r(ij)
      enddo
c
      lenb = lensec(lenbas)
      ibl3x = lenb*3 + 1 +ib
      call wrt3(p,lenbas,ibl3x,idaf)
      ibl3y = ibl3x + lenb
      call wrt3s(q,lenbas,idaf)
      ibl3z = ibl3y + lenb
      call wrt3s(r,lenbas,idaf)
      return
      end
      subroutine spchck
c
c     ----- check for pure sp basis -----
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
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
      logical ospbas, onocnt, ochek, opdbas, opfbas, opgbas
      integer kad, istd
      common /ijlab/ kad(mxshel),ospbas,onocnt,istd,ochek,
     +               opdbas,opfbas,opgbas
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
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
      ospbas = .true.
      opdbas = .false.
      opfbas = .false.
      opgbas = .false.
      onocnt = .true.
      istd = 1
      isptype = -1
      do 20 i = 1 , nshell
         kad(i) = 0
         if (ktype(i).gt.2) then
            if (onocnt) istd = i
            ospbas = .false.
            if (ktype(i).eq.3) opdbas = .true.
            if (ktype(i).eq.4) opfbas = .true.
            if (ktype(i).eq.5) opgbas = .true.
            onocnt = .false.
            kad(i) = -1
         end if
 20   continue
      if (intg76.eq.0) then
         onocnt = .true.
         istd = 1
         do 30 i = 1 , nshell
            kad(i) = -1
 30      continue
      end if
      if (nprint.eq.-5) return
c     if(.not.ospbas) write(iwr,110) ii
c110  format(/1x,'intg76 switch reset -- orbital ',i3)
      return
      end
      subroutine stvstv(ncall,q,iso,nshels)
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
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
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
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
      integer nmaxly, ntotly, icurly, isecly
      logical oprintm
      common /segm/ nmaxly,ntotly,icurly,isecly(100),
     +              oprintm
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
c
      integer ibl7s, ibl7t, ibl7f, ibl7st,ibl7tt,ibl7ft,ibl3qs
      integer ibl7qa,ibl7pa,ibl7ea,ibl7qb,ibl7pb,ibl7eb
      integer ibl7x, ibl7y, ibl7z
      integer ibl7so,ibl7o, ibl7la,ibl7ec,ibl7ha
      common /scra7/ ibl7s, ibl7t, ibl7f, ibl7st,ibl7tt,ibl7ft,ibl3qs,
     +               ibl7qa,ibl7pa,ibl7ea,ibl7qb,ibl7pb,ibl7eb,
     +               ibl7x, ibl7y, ibl7z,
     +               ibl7so,ibl7o(8), ibl7la, ibl7ec, ibl7ha
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
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
caleko
      logical vector
      common/nottwi/obeen,obeen2,obeen3,obeen4
caleko
      common/junk/s(225),g(225),
     *pint,qint,rint,t,p0,q0,r0,pi,qi,ri,pj,qj,rj,ni,nj,
     *tol,ii,jj,lit,ljt,mini,minj,maxi,maxj,iandk
c
      real*8 pp, u, w, dji
      integer nroots, nroot2
      common /root/ pp,u(12),w(12),nroots,nroot2,dji(60)
c
      common/blkin/dxyz(4),gg(225),ft(225),fx(225),dij(225),
     + pin(125),qin(125),rin(125),
     + ijx(225),ijy(225),ijz(225)
c mechanics
c
      real*8 tchg
      integer ifree, kform
      logical opress,otemp,omapp,omodel,ocharg,oaminc,oamhes,osubs
      common /modj/ kform,opress,otemp,omapp,omodel,ocharg,
     +              oaminc,oamhes,tchg(maxat),osubs,ifree
c
      common/g80nb/natmod,nnbg,mapmod(maxat),ipnonb(maxat+1),
     +             inonb(maxat*6)
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
c ***** omit specified charges from attraction terms ***
c---------
c  omtchg - omit 1e- integrals involving the specified charges
c
c   integrals of the form  < bf(a) | q(b) | bf (a) > are deleted
c
c   where a and b are centres indicated as follows
c
c  ztagcc - list of centre names a
c  ztagcq - for each a in ztagcc, the list of centres b
c
c  omtslf - omit self-energy of point charge array (applied to all
c           atoms with centre label of form bq* )
c      CHECK
c      parameter (mxombq=20,mxomat=50)
c---------

c
c obqf  - special treatment of forces .. skip force contributions
c         on non-bq centres
c         for use in qm/mm optimisations of surrounding
c         assumption is that bq's in this case don't have basis
c         functions or ecps
c
c omitted charge specifications
c

      integer mxombq, mxomat
      parameter (mxombq=50,mxomat=20)
      character*8 ztagcc, ztagcq
      logical omtchg, omtslf, obqf
      common /chgcc/ ztagcc(mxomat),ztagcq(mxomat,mxombq),omtchg,
     &     omtslf,obqf
c crystal field
      integer isecfx 
      logical ocryst
      common/xfield/isecfx,ocryst
c
c     ZORA common
c
      logical ozora,oscalz,onlyp,onlyd,onlyf,onlyg,opzora,
     1        op2zora,onlys,ocontr,oint_zora,o1e_zora,osmall_zora,
     2        opre_zora,oso,onoext_z,oioraz,oatint_z,oatscf_z,
     3        oscalatz
      real*8 cspeed,critat_z
      integer nbas_zora,ibls_z,iblv_z,nwv_z,ibiso_z,nwiso_z,ibcor_z,
     1        ibsx_z,icoul_z,ibsy_z,ibsz_z,ibsmin_z,ibsplus_z,nwmat_z,
     2        ibtrin_z,ibtrout_z,ibscalf_z,nwtrat_z,ibscale_z,ibcoul_z,
     3        niter_z,int_zora,num_ext,igauge_z,nwcor_z,isexch_z,
     4        nat_z,ibdens_z,nwdens_z,ibscalao_z,ibshat_z,is_z,
     5        irest_z,numdu_z,ibldu_z,iblock_z,lshift_z,ibl7ew_z
      common/zorac/ cspeed,critat_z,ozora,oscalz,
     1              onlyp,onlyd,onlyf,onlyg,opzora,
     1              op2zora,onlys,ocontr,nbas_zora,
     2              oint_zora,o1e_zora,icoul_z,
     3              ibls_z,iblv_z,nwv_z,
     3              ibiso_z,nwiso_z,ibcor_z,nwcor_z,
     4              ibsx_z,ibsy_z,ibsz_z,ibsmin_z,ibsplus_z,nwmat_z,
     5              ibtrin_z,ibtrout_z,ibscalf_z,nwtrat_z,
     6              ibscale_z,ibcoul_z,ibdens_z,nwdens_z,osmall_zora,
     7              niter_z,opre_zora,int_zora,num_ext,isexch_z,
     8              oso,onoext_z,is_z,igauge_z,oioraz,
     9              nat_z,ibscalao_z,oatint_z,oatscf_z,ibshat_z,
     1              oscalatz,irest_z,numdu_z,ibldu_z,iblock_z,lshift_z,
     2              ibl7ew_z
c....    icoul_z : 1 : Full coulomb
c...               2 : atomic coulomb (default)
c...               3 : no coulomb (just 1-electron ints)
c...               4 : small basis coulomb
c...     isexch_z  0 : no exchange in zora correction (default)
c...               1 : exchange in zora correction (like in dft)
c...     niter_z : # iterations between coulomb matrix update
c...     int_zora :internal basis gen. : 0 none, 1, automatic, 2 by hand
c
c...     oscalz  : scaling on or off (true or false)
c...     oioraz  : change metric to perform inf. order scaling 
c...               (Dyall/EvLenthe (not quite))
c...     oscalatz: atoms only scaled zora (together with get atom)
c...     only'x' : internal basis only op to 'x'-type functions
c...     op(2)zora: print flags (op2zora+opzora: prints *everything*)
c...     ocontr  : adds l+1 block to internal basis in contracted way 
c...     nbas_zora: size of internal basis
c...     oint_zora: in internal or external mode???
c...     o1e_zora: one electron integrals have been calculated in *this*
c...               internal basis
c...     osmall_zora: small swap for projected coulomb
c...     opre_zora: pre_zora has been called
c...     oso     : spin orbit calculation 
c...     onoext_z: no external functions copied to internal basis
c...     is_z : section for orbitals used for  coulomb matrix
c...     igauge_z : add so many s-marices to h-matrix
c...                the energy should change by igauge
c...     nat_z  : indicates (if ne 0) that starting density is used
c...              to calc. coulomb; gives # it. to get atomic density
c...     critat_z : accuracy of atomic energy for atomic zora
c...     oatint_z : indicates that integrals to be calculated are atomic
c...     oatscf_z : indicates we are doing a set of atomic SCFs
c...     ibshat_z : block on num8 to store blocked s/h
c...     for disk-adress and lengths (ib.., nw.. ) see subroutine atoms2
c...     irest_z  : .ne.0 indicates atomic zora corrections from dumpfile; 1 not on restart
c...     numdu_z,ibldu_z : if ne 0 other dumpfile to read zora corrections from
c...     iblock_z is the block number on ed3, where the zora corrections (section 425) reside
c...     ibl7ew is where the energy weighted density matrix resides
      character *8 title,guess
      common/restrz/title(12),guess
c
      character*8 field
      character*24 version
      common /drfopt/ field,version
      logical oreact
      integer intdrf
      common /drfopti/ oreact,intdrf
c
      dimension q(*),iso(nshels,*)
      dimension ix(35),iy(35),iz(35),jx(35),jy(35),jz(35)
      dimension m0(48)
c ******
cdrf
c     drf extension
c     ====================================================================
c         note: in hondo, o,x,y and z are real*8
          common/hefcpar/edumm(5,1000),nedumm,iefc
          common/hfldpar/fldxyz(3),ifld
         integer idafh,navh
         common/hdafile/idafh,navh,ioda(2,1000)
         common/c_of_m/pcm,qcm,rcm
c
      integer nsamp, nblock, nmoves, ncheck
      common /mcrunp1/ nsamp,nblock,nmoves,ncheck
c
      real*8 darot, dtrans, temp, onekt, excld, delmax, enmin
      common /mcrunp2/ darot,dtrans,temp,onekt,excld,delmax,enmin
c
      logical mcupdt,secstat,excess
      integer notrans, norot, imcout, iseed, imcst, iactst
      common /mcrunp3/ notrans,norot,imcout,iseed,imcst,iactst,
     +               mcupdt,secstat,excess
c
      character*16 outfor
      common /mcrunp4/ outfor
c
      real*8 ratmin, ratmax, amxrot, amxtrn, amnrot, amntrn
      common /mcrunp5/ ratmin,ratmax,amxrot,amxtrn,amnrot,amntrn
c
      real*8 gammc
      common /mcrunp6/ gammc(5)
c
c
      integer isingl, nbits
      common /machin1/ isingl,nbits
c
      integer ixddaf
      integer ndar10, ndar20, ndar31, ndar41, ndar42, ndar43, ndar44
      integer ndar45, ndar46, ndar47, ndar48, ndar49
      common /dafindx/ ixddaf(8192),ndar10,ndar20,
     +                 ndar31,ndar41,ndar42,ndar43,ndar44,ndar45,
     +                 ndar46,ndar47,ndar48,ndar49
c
      integer lrec10, lrec20, lrec31, lrec41, lrec42, lrec43
      integer lrec44, lrec45, lrec46, lrec47, lrec48, lrec49
      common /dafrec/ lrec10,lrec20,lrec31,lrec41,lrec42,lrec43,
     +                lrec44,lrec45,lrec46,lrec47,lrec48,lrec49
c
cafc
c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     this file contains the parameter settings for piet van duijnen's
c     reactionfield programme together with the sizes needed for the
c     v. duijnen - gamess interface as developed by 
c     nederkoorn and v. lenthe
c
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c     ==================================================================
c
c     parameter setting for reaction field programme van duijnen:
c
c=======================================================================
c
c     
      parameter (mxat = 10000)
      parameter (mxpol = 10000)
      parameter (mxpts = 10000)
      parameter (mxspe = 50)
      parameter (mxamb = 10)
      parameter (mxgrp = 500)
      parameter (mxnum = 1024)
      parameter (mxsh = 512)
      parameter (mxgs = 2048)
      parameter (mxnpts = 20)
      parameter (mxneq = 255)
      parameter (mcwrkmx = 10)
c
c ==================================================================
c
c     parameter setting new drf routines (drf/dimpar):
c
c===================================================================
c
      parameter (mxst = 5)
      parameter (mxgran = 10)
cmw      parameter (mxnpts = 50)
      parameter (mxgrpar = (mxgran*(mxgran+1))/2)
c
c     ==================================================================
c
c     parameter setting for v. duijnen - gamess interface:
c
c=======================================================================
c
      parameter (ncatom=500)
c
c  ncatom is maximal allowed classical atoms in zmat
c  total allowed atoms in zmat = maxnz + ncatom
c
c  maxrfa is maximal number of classical atoms in zmat + coordinates
c  total maximal atoms allowed is maxrfa + maxat
c
      parameter (maxrfa=mxpts)
c
c  maxamb is maximal number of ambigous atoms. it is set equal
c  to maximal z-cards possible in gamess input , i.e. 50
c
      parameter (maxamb=mxamb)
c
c=============================================
c for the pes-routines
c=============================================
c
      parameter (mxgr1=128, mxpes=100)
c
c=============================================
c for the drfexpc routines
c=============================================
c
      parameter (maxex=1000)
c
c=============================================
c for the drfsurf routine
c=============================================
c
      parameter (maxinv=10,maxnrc=10000)
c      parameter (maxnrc=10000)
c
c=============================================
c for the analrf routine
c=============================================
c
c
      character*8 scftyp
      common /scfopt2/ scftyp
c
c
      integer ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep
      integer neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga
      integer idipcal,iclinte,iclintd,ieffpol,isodis,iclintr
      integer iextdip,itsolv,itermax,isolsav,imomsav,idisadd
      integer ianal,ineqex,ithole,irevdis,ihbond,igrppol
      integer maxci_drf
      common/drfpar1/ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep,
     +               neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga,
     +               idipcal,iclinte,iclintd,ieffpol,isodis,iclintr,
     +               iextdip,itsolv,itermax,isolsav,imomsav,idisadd,
     +               ianal,ineqex,ithole,irevdis,ihbond,igrppol,
     +               maxci_drf
c
      integer nxtpts, npol, nexp, natint, nshint, namb, nspec, ngran
      common/drfpar2/nxtpts,npol,nexp,natint,nshint,namb,nspec,ngran
c
      integer ndima, ndimb, nbem, ndim, nexp4, npol3, nwtr, nwtc, nzfa
      integer lomga, nomga, nchd, ngrnam, ngrnam2, nzfp, nzfn
      integer nneq, nneqrf, maxneq, neqdim, nqdim, nqcls
      common/drfpar3/ndima,ndimb,nbem,ndim,nexp4,npol3,nwtr,nwtc,nzfa,
     +               lomga,nomga,nchd,ngrnam,ngrnam2,nzfp,nzfn,
     +               nneq,nneqrf,maxneq,neqdim,nqdim,nqcls
c
      integer ifldin, ifldout, idrfout, modxza, irepopt
      integer iexpza, icmexp, iadexp, nodpe, igetden, iarfcal, iqmclr
      common/drfpar4/ifldin,ifldout,idrfout,modxza,irepopt,
     +               iexpza,icmexp,iadexp,nodpe,igetden,iarfcal,
     +               iqmclr
c
      real*8 gamdrf, dstmin, dstmax, hbondl, hbondr, afact, rfact
      real*8 cvgrel, agrpe, agrpm, agrpc, scffact, acur
      real*8 conci_drf
      common/drfpar5/gamdrf,dstmin,dstmax,hbondl,hbondr,afact,rfact,
     +               cvgrel,agrpe,agrpm,agrpc,scffact,acur,conci_drf
c
      integer iind, isur1, isur2, ilwt, ilvr, illur, ilindx
      common/drfpar6/iind,isur1,isur2,ilwt,ilvr,illur,ilindx
          character*4 keyfld, keyefc, iefc, ifld
          data keyfld, keyefc /' fld',' efc'/
cdrf  ===================  end drfexts ===============================
      character*8 fnm
      character*6 snm
      data fnm,snm/"intega.m","stvstv"/
      data  m51/51/
c mechanics
      data dzero,pt5,done,two,three,five,seven,
     +     rnine,eleven  /0.0d0,0.5d0,1.0d0,
     + 2.0d0,3.0d0,5.0d0,7.0d0,9.0d0,11.0d0/
      data pi212 /1.1283791670955d0/
      data sqrt3 /1.73205080756888d0/
      data sqrt5 /2.23606797749979d0/
      data sqrt7 /2.64575131106459d0/
      data rln10 /2.30258d0/
c
      data jx / 0, 1, 0, 0, 2, 0, 0, 1, 1, 0,
     +          3, 0, 0, 2, 2, 1, 0, 1, 0, 1,
     +          4, 0, 0, 3, 3, 1, 0, 1, 0, 2,
     +          2, 0, 2, 1, 1/
      data ix / 1, 6, 1, 1,11, 1, 1, 6, 6, 1,
     *         16, 1, 1,11,11, 6, 1, 6, 1, 6,
     *         21, 1, 1,16,16, 6, 1, 6, 1,11,
     *         11, 1,11, 6, 6/
      data jy / 0, 0, 1, 0, 0, 2, 0, 1, 0, 1,
     +          0, 3, 0, 1, 0, 2, 2, 0, 1, 1,
     +          0, 4, 0, 1, 0, 3, 3, 0, 1, 2,
     +          0, 2, 1, 2, 1/
      data iy / 1, 1, 6, 1, 1,11, 1, 6, 1, 6,
     +          1,16, 1, 6, 1,11,11, 1, 6, 6,
     +          1,21, 1, 6, 1,16,16, 1, 6,11,
     +          1,11, 6,11, 6/
      data jz / 0, 0, 0, 1, 0, 0, 2, 0, 1, 1,
     +          0, 0, 3, 0, 1, 0, 1, 2, 2, 1,
     +          0, 0, 4, 0, 1, 0, 1, 3, 3, 0,
     +          2, 2, 1, 1, 2/
      data iz / 1, 1, 1, 6, 1, 1,11, 1, 6, 6,
     +          1, 1,16, 1, 6, 1, 6,11,11, 6,
     +          1, 1,21, 1, 6, 1, 6,16,16, 1,
     +         11,11, 6, 6,11/
      data maxbit /64/
      l1 = num
      l2 = (num*(num+1))/2
      l3 = num*num
      l4 = natmod*3
      l5 = natmod
      ono  = .false.
      outv = oprint(59)
      nav = lenwrd()
      if (nprint.eq.-5) outv = .false.
c
c     ----- set pointers for partitioning of core -----
c
      i10 = 0
      i20 = i10 + l2
      i30 = i20 + l2
      if(oreact) then
cdrf  i301 space for hamilonian in fldint
       i301 = i30 + l2
       last = i301 + l2
      endif
      if (lpseud.le.1) then
         i40 = i30 + l2
         last = i40 + l2 + l1
      else
         i40 = i30 + l2
         i50 = i40 + l2
         i60 = i50 + l3
         i70 = i60 + 225*num
         i80 = i70 + 225*20
         last = i80 + (nat*nt+nav-1)/nav
      end if
c
c crystal field correction
c
      if(ocryst)then
         ic10 = last
         last = ic10 + l2
      endif
c ...
c mechanics
c ...
      if (oaminc) then
         ia40 = last
         ia50 = ia40 + l4
         ia60 = ia50 + l5
         last = ia60 + l5 + 1
      end if
c
c ***** allocate memory for charge exclusion
c
      if (omtchg) then
         i110 = last
         last = i110 + nat * nat
      endif
c
      length = last - i10

      i10 = igmem_alloc_inf(length,fnm,snm,"i10",IGMEM_DEBUG)
      i20 = i10 + l2
      i30 = i20 + l2
      if (lpseud.le.1) then
         i40 = i30 + l2
         last = i40 + l2 + l1
      else
         i40 = i30 + l2
         i50 = i40 + l2
         i60 = i50 + l3
         i70 = i60 + 225*num
         i80 = i70 + 225*20
         last = i80 + (nat*nt+nav-1)/nav
      end if
c
c crystal field correction
c
      if(ocryst)then
         ic10 = last
         last = ic10 + l2
      endif
c ...
c mechanics
c ...
      if (oaminc) then
         ia40 = last
         ia50 = ia40 + l4
         ia60 = ia50 + l5
         last = ia60 + l5 + 1
      end if
c ***** allocate memory for charge exclusion
      if (omtchg) then
         i110 = last
         last = i110 + nat * nat
      endif
c
      if (ncall.eq.0 .and. nprint.eq.3) write (iwr,6010) i10 , i20 ,
     +    i30 , last

c **************
      if (omtchg) then
       call omitch(q(i110),nat,nona)
      endif
c **************
      call chek1e(o1e)
      if (oint_zora) o1e = o1e_zora
      if (oatscf_z) o1e = .false.
      if (oreact. and. .not. obeen) o1e = .false.
      if (.not.o1e .or. mcupdt) then
c...     avoid to reset if atomic zora corrections are used
c...     or if we are in the middle of atomic startup (first point)
         if (nat_z.le.0.and.guess.ne.'atoms')
     1       opre_zora = .false.
c
         if (oaminc) then
c     retrieve mechanics coordinates. xmod
            call secget(isect(472),m51,iblk172)
            call rdedx(q(ia40),l4,iblk172,idaf)
c                      retrieve nbact array. nbact
            call secget(isect(473),m51,iblk173)
            idum = (natmod-1) /nav + 1
            call rdedx(q(ia60),idum,iblk173,idaf)
c                      retrieve charges. chgmod
            call secget(isect(474),m51,iblk174)
            call rdedx(q(ia50),l5,iblk174,idaf)
         end if
c
c     ----- calculate -s- and -h0- matrices -----
c
c     - s- at x(i10)
c     -h0- at x(i20)
c
         cpu = cpulft(1)
c
         if (ncall.eq.0) call cpuwal(begin,ebegin)
         if (ncall.eq.0 .and. outv) write (iwr,6030) cpu
         tol = rln10*itol
         out = nprint.eq.3
         if (out) then
           do ii = 1,6
            oprn(30+ii) = .true.
           enddo
          else
           do ii = 1,6
            if(oprn(30+ii)) out = .true.
           enddo
         endif
         onorm = normf.ne.1 .or. normp.ne.1
         call vclr(q(i10),1,l2)
         call vclr(q(i20),1,l2)
         call vclr(q(i30),1,l2)
c
c     ----- ishell
c
         do 440 ii = 1 , nshell
            i = katom(ii)
c
c     ----- eliminate ishell -----
c
            do 450 it = 1 , nt
               id = iso(ii,it)
               if (id.gt.ii) go to 440
               m0(it) = id
450         continue
            icent = i
            pi = c(1,i)
            qi = c(2,i)
            ri = c(3,i)
            i1 = kstart(ii)
            i2 = i1 + kng(ii) - 1
            lit = ktype(ii)
            mini = kmin(ii)
            maxi = kmax(ii)
            loci = kloc(ii) - mini
c
c     ----- jshell
c
            do 430 jj = 1 , ii
               j = katom(jj)
               if (oatint_z.and.j.ne.icent) go to 430
               n2 = 0
               do 470 it = 1 , nt
                  jd = iso(jj,it)
                  if (jd.gt.ii) go to 430
                  id = m0(it)
                  if (id.lt.jd) then
                     nd = id
                     id = jd
                     jd = nd
                  end if
                  if (id.eq.ii .and. jd.gt.jj) go to 430
                  if (id.eq.ii.and.jd.eq.jj) then
                      n2 = n2 + 1
                  end if
470            continue
               q2 = dble(nt)/dble(n2)
               jcent = j
               pj = c(1,j)
               qj = c(2,j)
               rj = c(3,j)
               if (oatint_z) then
c..              atomic scf all on origin
                 pi = 0.0d0
                 qi = 0.0d0
                 ri = 0.0d0
                 pj = 0.0d0
                 qj = 0.0d0
                 rj = 0.0d0
               end if
               j1 = kstart(jj)
               j2 = j1 + kng(jj) - 1
               ljt = ktype(jj)
               minj = kmin(jj)
               maxj = kmax(jj)
               locj = kloc(jj) - minj
               nroots = (lit+ljt-2)/2 + 1
               rr = (pi-pj)**2 + (qi-qj)**2 + (ri-rj)**2
               oiandj = ii.eq.jj
c
c     ----- prepare indices for pairs of (i,j) functions
c
               ij = 0
               max = maxj
               do 30 i = mini , maxi
                  nnx = ix(i)
                  nny = iy(i)
                  nnz = iz(i)
                  if (oiandj) max = i
                  do 20 j = minj , max
                     ij = ij + 1
                     ijx(ij) = nnx + jx(j)
                     ijy(ij) = nny + jy(j)
                     ijz(ij) = nnz + jz(j)
                     if (j.le.1) then
                        ft(ij) = three
                     else if (j.le.4) then
                        ft(ij) = five
                     else if (j.le.10) then
                        ft(ij) = seven
                     else if (j.gt.20) then
                        ft(ij) = eleven
                     else
                        ft(ij) = rnine
                     end if
 20               continue
 30            continue
c
               do 40 i = 1 , ij
                  s(i) = dzero
                  gg(i) = dzero
                  g(i) = dzero
                  fx(i) = dzero
 40            continue
c
c     ----- i primitive
c
               jgmax = j2
               do 400 ig = i1 , i2
                  ai = ex(ig)
                  arri = ai*rr
                  axi = ai*pi
                  ayi = ai*qi
                  azi = ai*ri
                  csi = cs(ig)
                  cpi = cp(ig)
                  cdi = cd(ig)
                  cfi = cf(ig)
                  cgi = cg(ig)
c
c     ----- j primtive
c
                  if (oiandj) jgmax = ig
                  do 390 jg = j1 , jgmax
                     aj = ex(jg)
                     aa = ai + aj
                     aa1 = done/aa
                     dum = aj*arri*aa1
                     if (dum.le.tol) then
                        fac = dexp(-dum)
                        csj = cs(jg)
                        cpj = cp(jg)
                        cdj = cd(jg)
                        cfj = cf(jg)
                        cgj = cg(jg)
                        ax = (axi+aj*pj)*aa1
                        ay = (ayi+aj*qj)*aa1
                        az = (azi+aj*rj)*aa1
                        odoub = oiandj .and. ig.ne.jg
c
c     ----- density factor
c
                        max = maxj
                        nn = 0
                        do 220 i = mini , maxi
                           go to (50,60,120,120,
     +                            70,120,120,80,120,120,
     +                            90,120,120,100,120,120,120,120,120,
     +                            110,
     +                            112,120,120,114,120,120,120,120,120,
     +                            116,120,120,118,120,120), i
c
 50                        dum1 = csi*fac
                           go to 120
 60                        dum1 = cpi*fac
                           go to 120
 70                        dum1 = cdi*fac
                           go to 120
 80                        if (onorm) dum1 = dum1*sqrt3
                           go to 120
 90                        dum1 = cfi*fac
                           go to 120
 100                       if (onorm) dum1 = dum1*sqrt5
                           go to 120
 110                       if (onorm) dum1 = dum1*sqrt3
                           go to 120
 112                       dum1 = cgi*fac
                           go to 120
 114                       if (onorm) dum1 = dum1*sqrt7
                           go to 120
 116                       if (onorm) dum1 = dum1*sqrt5/sqrt3
                           go to 120
 118                       if (onorm) dum1 = dum1*sqrt3
 120                       if (oiandj) max = i
                           do 210 j = minj , max
                              go to (130,140,200,200,
     +                               150,200,200,160,200,200,
     +                               170,200,200,180,200,200,
     +                               200,200,200,190,
     +                               192,200,200,194,200,200,200,200,
     +                               200,196,200,200,198,200,200),j
 130                          dum2 = dum1*csj
                              if (odoub) then
                                 if (i.gt.1) then
                                    dum2 = dum2 + csi*cpj*fac
                                 else
                                    dum2 = dum2 + dum2
                                 end if
                              end if
                              go to 200
 140                          dum2 = dum1*cpj
                              if (odoub) dum2 = dum2 + dum2
                              go to 200
 150                          dum2 = dum1*cdj
                              if (odoub) dum2 = dum2 + dum2
                              go to 200
 160                          if (onorm) dum2 = dum2*sqrt3
                              go to 200
 170                          dum2 = dum1*cfj
                              if (odoub) dum2 = dum2 + dum2
                              go to 200
 180                          if (onorm) dum2 = dum2*sqrt5
                              go to 200
 190                          if (onorm) dum2 = dum2*sqrt3
                              go to 200
 192                          dum2 = dum1*cgj
                              if (odoub) dum2 = dum2 + dum2
                              go to 200
 194                          if (onorm) dum2 = dum2*sqrt7
                              go to 200
 196                          if (onorm) dum2 = dum2*sqrt5/sqrt3
                              go to 200
 198                          if (onorm) dum2 = dum2*sqrt3
 200                          nn = nn + 1
                              dij(nn) = dum2
 210                       continue
 220                    continue
c
c     ----- overlap and kinetic energy
c
                        t = dsqrt(aa1)
                        t1 = -two*aj*aj*t
                        t2 = -pt5*t
                        p0 = ax
                        q0 = ay
                        r0 = az
                        in = -5
                        do 240 i = 1 , lit
                           in = in + 5
                           ni = i
                           do 230 j = 1 , ljt
                              jn = in + j
                              nj = j
                              call stvint
                              pin(jn) = pint*t
                              qin(jn) = qint*t
                              rin(jn) = rint*t
                              nj = j + 2
                              call stvint
                              pin(jn+25) = pint*t1
                              qin(jn+25) = qint*t1
                              rin(jn+25) = rint*t1
                              nj = j - 2
                              if (nj.gt.0) then
                                 call stvint
                              else
                                 pint = dzero
                                 qint = dzero
                                 rint = dzero
                              end if
                              n = (j-1)*(j-2)
                              dum = dble(n)*t2
                              pin(jn+50) = pint*dum
                              qin(jn+50) = qint*dum
                              rin(jn+50) = rint*dum
 230                       continue
 240                    continue
                        do 250 i = 1 , ij
                         nnx = ijx(i)
                         nny = ijy(i)
                         nnz = ijz(i)
                         pyz = qin(nny)*rin(nnz)
                         dum = pyz*pin(nnx)
                         dum1 = (pin(nnx+25)+pin(nnx+50))*pyz + 
     +                    (qin(nny+25)+qin(nny+50))*pin(nnx)*rin(nnz) + 
     +                    (rin(nnz+25)+rin(nnz+50))*pin(nnx)*qin(nny)
                         s(i) = s(i) + dij(i)*dum
                         g(i) = g(i) + dij(i)*(dum*aj*ft(i)+dum1)
                         gg(i) = gg(i) + dij(i)*(dum*aj*ft(i)+dum1)
 250                    continue
c
c     ----- nuclear attraction
c
                        dum = pi212*aa1
c                       facinv = aa/(fac*pi212)
                        do 260 i = 1 , ij
                           dij(i) = dij(i)*dum
 260                    continue
                        aax = aa*ax
                        aay = aa*ay
                        aaz = aa*az

                        do 320 ic = 1 , nat

                           if(omtchg .and. icent.eq.jcent) then
c
c ****** charge exclusion
c
                              pnuc = -q(i110-1+(ic-1)*nat + icent)
                           else if ((oatint_z.or.
     1                              (oint_zora.and.icoul_z.eq.2)).and.
     1                             (ic.ne.icent.or.ic.ne.jcent)) then 
c**ZORA
                              pnuc = 0.0d0
                           else
                              pnuc = -czan(ic)
                           endif
                           if(ocryst)then
c
c ****** crystal field point charge
c
                              obq = (zaname(ic)(1:2).eq.'bq')
                           endif

                           cx = c(1,ic)
                           cy = c(2,ic)
                           cz = c(3,ic)
                           if (oatint_z) then
c....                        atomic zora all on origin
                             cx = 0.0d0
                             cy = 0.0d0
                             cz = 0.0d0
                           end if
                           pp = aa*((ax-cx)**2+(ay-cy)**2+(az-cz)**2)
                           if (nroots.le.3) call rt123
                           if (nroots.eq.4) call roots4
                           if (nroots.eq.5) call roots5
                           mm = 0
                           do 290 k = 1 , nroots
                              uu = aa*u(k)
                              ww = w(k)*pnuc
                              tt = done/(aa+uu)
                              t = dsqrt(tt)
                              p0 = (aax+uu*cx)*tt
                              q0 = (aay+uu*cy)*tt
                              r0 = (aaz+uu*cz)*tt
                              in = -5 + mm
                              do 280 i = 1 , lit
                                 in = in + 5
                                 ni = i
                                 do 270 j = 1 , ljt
                                    jn = in + j
                                    nj = j
                                    call stvint
                                    pin(jn) = pint
                                    qin(jn) = qint
                                    rin(jn) = rint*ww
 270                             continue
 280                          continue
                              mm = mm + 25
 290                       continue
                           do 310 i = 1 , ij
                            nnx = ijx(i)
                            nny = ijy(i)
                            nnz = ijz(i)
                            dum = dzero
                            mm = 0
                            do 300 k = 1 , nroots
                               dum = dum + pin(nnx+mm)*qin(nny+mm)
     +                               *rin(nnz+mm)
                               mm = mm + 25
 300                        continue
                            g(i) = g(i) + dum*dij(i)
                            if(ocryst.and.obq)
     +                           fx(i) = fx(i) + dum*dij(i)
 310                       continue
 320                    continue

                        if (oaminc) then
c ...
c    skip this code if dealing with a dummy atom.
c ...
                           ic3 = -4
                           do 380 icmod = 1 , natmod
                              ic3 = ic3 + 3
                              call skip80(icent,jcent,icmod,q(ia60),
     +                           oskmod)
                              if (.not.(oskmod)) then
                                 ia40p = ia40 + ic3
                                 pnuc = -q(ia50+icmod-1)
                                 cx = q(ia40p+1)
                                 cy = q(ia40p+2)
                                 cz = q(ia40p+3)
                                 pp = aa*((ax-cx)**2+(ay-cy)**2+(az-cz)
     +                                **2)
                                 if (nroots.le.3) call rt123
                                 if (nroots.eq.4) call roots4
                                 if (nroots.eq.5) call roots5
                                 mm = 0
                                 do 350 kmod = 1 , nroots
                                    uu = aa*u(kmod)
                                    ww = w(kmod)*pnuc
                                    tt = done/(aa+uu)
                                    t = dsqrt(tt)
                                    p0 = (aax+uu*cx)*tt
                                    q0 = (aay+uu*cy)*tt
                                    r0 = (aaz+uu*cz)*tt
                                    in = -5 + mm
                                    do 340 imod = 1 , lit
                                       in = in + 5
                                       ni = imod
                                       do 330 jmod = 1 , ljt
                                         jn = in + jmod
                                         nj = jmod
                                         call stvint
                                         pin(jn) = pint
                                         qin(jn) = qint
                                         rin(jn) = rint*ww
 330                                   continue
 340                                continue
                                    mm = mm + 25
 350                             continue
                                 do 370 imod = 1 , ij
                                  nnx = ijx(imod)
                                  nny = ijy(imod)
                                  nnz = ijz(imod)
                                  dum = 0.0d0
                                  mm = 0
                                  do 360 kmod = 1 , nroots
                                     dum = dum + pin(nnx+mm)*qin(nny+mm)
     +                                  *rin(nnz+mm)
                                     mm = mm + 25
 360                              continue
                                  g(imod) = g(imod) + dum*dij(imod)
 370                             continue
                              end if
 380                       continue
                        end if
                     end if
c ...
c ...
 390              continue
 400           continue
c
c
c     - s- at x(i10)
c     -h0- at x(i20)
c     -t-  at x(i30)
c
c
c     ----- set up overlap and h-core matrices
c
               max = maxj
               nn = 0
               do 420 i = mini , maxi
                  li = loci + i
                  in = (li*(li-1))/2
                  if (oiandj) max = i
                  do 410 j = minj , max
                     lj = locj + j
                     jn = lj + in
                     nn = nn + 1
                     q(jn-1+i10) = s(nn)*q2
                     q(jn-1+i20) = g(nn)*q2
                     q(jn-1+i30) = gg(nn)*q2
                     if(ocryst)q(jn-1+ic10) = fx(nn)*q2
 410              continue
 420           continue
 430        continue
 440     continue
c
         if (lpseud.gt.0) i20a = igmem_alloc_inf(l2,fnm,snm,"i20a",
     +                                           IGMEM_DEBUG)
         if (lpseud.eq.1) then
c
c ---- if ecp pseudopotentials required, call ecpint
c ----  keep them in i20a
c
            call vclr(q(i40),1,l2)
            call vclr(q(i20a),1,l2)
            call ecpint(q,q(i40),q(i20a),iso,nshell,ono,out)
         endif
         if (lpseud.eq.2) then
c
c ---- if nonlocal pseudopotentials required, call xpsnlc
c
            call vclr(q(i20a),1,l2)
            call xpsnlc(q(i20a),q(i50),q(i60),q(i40),q(i70),q(i80),
     +           iso,l2,nat,num,nshell,ono)
         end if
         if (lpseud.gt.0) call vadd(q(i20a),1,q(i20),1,q(i20),1,l2)
c
c     ----- symmetrize 1-e integral arrays
c
         call sym1e(q(i10),q(i40),iso,1,nshell)
         call sym1e(q(i20),q(i40),iso,1,nshell)
         call sym1e(q(i30),q(i40),iso,1,nshell)
         if (lpseud.gt.0) call gmem_free_inf(i20a,fnm,snm,"i20a")
c
         if (ocryst) 
     +    call sym1e(q(ic10),q(i40),iso,1,nshell)
c
c...   Gauge shift for testing Zora schemes
c
         if (igauge_z.ne.0) then
            if (oatint_z) call caserr('gauge+atint_z not obvious')
            write(iwr,601) igauge_z,ne,igauge_z*1.0d0
601      format(1x,/,' ****************************************'
     1            ,/,' ** adding ',i10,   ' to potential     **'
     2            ,/,' ** with   ',i10,   ' electrons        **'
     3            ,/,' ** the energy should change',f10.0 ,' **'
     4            ,/,' ****************************************')
c...      add s to h
            do i=1,l2
               q(i20+i-1) =  q(i20+i-1) +  q(i10+i-1)*igauge_z /  
     +         dble(ne)
            end do
         end if
c
c   ----- output 1-electron integrals to dumpfile
c
         if (oint_zora) then
c...     write s,v=h-t
            if (l2.ne.nwv_z) call caserr(' zora problem in stvstv')
            call wrt3(q(i10),l2,ibls_z,num8)
            do i=1,l2
               q(i30+i-1) = q(i20+i-1) - q(i30+i-1)
            end do
            call wrt3(q(i30),l2,iblv_z,num8)
            o1e_zora = .true.
            call gmem_free_inf(i10,fnm,snm,"i10")
            return
         else if (oatscf_z) then
c...        block and write 1-electron ints for atoms only scf
            call wrt3(q(i10),l2,ibshat_z,num8)  
            call wrt3(q(i20),l2,ibshat_z+lensec(l2),num8)  
            call gmem_free_inf(i10,fnm,snm,"i10")
            return
         else
            call sec192(q(i10),q(i20),q(i30),l2)
         end if
c        call sec192(q(i10),q(i20),q(i30),l2)
         if (oreact) then
          if (field .eq. ' ' .and. .not.obeen2) then
c
c -- for Dipole Preserving Charge analysis
c    in the absence of a DRF calculation
c
           isingl=lenwrd()
           nbits =maxbit/isingl
           vector=.false.
           idafh=110
           ndar=999
           call daopen(idafh,ioda,navh,ndar)
           call dawrit(idafh,ioda,q(i10),l2,12,navh)
           scftyp = zscftp  
           obeen2 = .true.
          endif
          if (field .ne. '  ' .and. .not. obeen) then
c    h0, s and t are safe on sec192 now, their space is
c    changed by fldint (via dipxyz)
c    dipxyz needs scratch-space for h0, i301 is used
c    h0 in fldint is taken and written in dafile 11
           if (idrfout .ge. 3)
     1     write (iwr,6031) 'drf-routines initialised'
c in rfin:         call inithondo(0)
           call inithondo(1)
           call dawrit(idafh,ioda,q(i10),l2,12,navh)
           call dawrit(idafh,ioda,q(i20),l2,11,navh)
           call dawrit(idafh,ioda,q(i30),l2,41,navh)
c    -written the hcore to da10 ,11 because hondo expects it there
c     taken from 'inttel.f' hondo style (375:377)
c    -efcint(hefc,h-core)  =>  efcint(q(i10),q(i20))
c    the efcint is not effective cope pot. what is it?
c           if (iefc.eq.keyefc) then
c              print *,'standv: call efcint'
c if used: activate gen/opt: /hefcpar/
c              call efcint(q(i10),q(i20))
c           endif
c    -fldint(h-core,x,y,z)    =>  fldint(q(i20),q(x),q(y),q(z))
c    optional fieldint. xfld,yfld and zfld are
c    taken from gamess /field/ which is only used by runtyp=analy
c    (maybe 'field'-directive to add always allowing fields)
c    q(i301) will contain the changed hamiltonian which
c    must be written to sec192 in order to affect anything
c           if (ifld.eq.keyfld) then
c              print *,'standv: call fldint'
c if used: activate gen/opt: /hfldpar/
c              call fldint(q(i301),q(i10),q(i20),q(i30))
c           endif
          endif
         endif
cdrf =================================================end drf extension
c
c     ----- now compute dipole moment integrals
c
cdrf
         pcm=dzero
         qcm=dzero
         rcm=dzero
cdrf
         call dipxyz(q(i10),q(i20),q(i30))
c
c     ---- load to ionsec
c
        if (oreact) then 
c
cdrf =================================================== drf extension
         if (.not.obeen3) then
c    in hondo, dipoles are at rec.53,54 and 55 resp. on the idaf file.
           call dawrit(idafh,ioda,q(i10),l2,53,navh)
           call dawrit(idafh,ioda,q(i20),l2,54,navh)
           call dawrit(idafh,ioda,q(i30),l2,55,navh)
           obeen3 = .true.
         endif
         call secdip(q(i10),q(i20),q(i30),l2)
c
c     ----- reset core memory -----
c
cahv  -- first free memory used by vacuum stvstv
c
         call gmem_free_inf(i10,fnm,snm,"i10")
c
         if (field .ne. '  ' .and. .not. obeen4 ) then
c
c    here, the external field options must be included:
c    to be called: -drfexpc
c    (calculate expansion centra if necessary)
c    and           -drfzfa
c    (calc pot. due to external charges at exp.cent)
c    (not necessary if second or later non-eq rf)
c    and           -rfcal
c    (calc necessary mats for rf-evaluation)
c
           if (.not. mcupdt) write(iwr,1011)
 1011 format(//,1x,104('-'),//)
           if (idrfout .ge. 3) 
     1       print *,"standv: call to drfexpc"
           if (.not. mcupdt) write(iwr,1021) 
 1021      format(/,' Defining expansion centres',
     1 ' for DRF embedding')
           call drfexpc(q)
cahv
	   i10 = igmem_alloc_inf(3*nexp,fnm,snm,"i10a",IGMEM_DEBUG)
           if (idrfout .ge. 3) 
     1       print *,"standv: call to drfzfa"
           if (.not. mcupdt) write(iwr,1031) 
 1031      format(/,' Calculating potential/field',
     1 ' for DRF embedding')
           call drfzfa(q(i10))
c
	   call gmem_free_inf(i10,fnm,snm,"i10a")
c
c  reset core memory
c          call setc(loadcm)
           if (idrfout .ge. 3) 
     1       print *,"standv: call to rfcal"
           if (.not. mcupdt) write(iwr,1041) 
 1041      format(/' Calculating couplings for DRF response embedding')
           call rfcal(q)
           obeen4 = .true.
          endif
         else
c
         call secdip(q(i10),q(i20),q(i30),l2)
c
         endif
cdrf ======================================================== end drf extension
c
      if(ocryst)then
c
c   ---- save crystal field part
c
         m53=53  
         write(6,*)'writing section',isecfx,m53
         len = lensec(l2)
         call secput(isecfx,m53,len,iblkcr)
         call wrt3(q(ic10),l2,iblkcr,idaf)
      endif
c
c     ---- print if requested
c
      if (out) call prt1e(q(i10),ionsec)
c
         cpu = cpulft(1)
c
         if (ncall.eq.0) then
            if (outv) write (iwr,6020) cpu
            tim = cpulft(1)
            if (tim.ge.timlim) then
               nindmx = -1
               irest = 1
            end if
            call timana(3)
         end if
      end if
c
c     ----- reset core memory -----
c
      if (oreact) then
         if (o1e .and. .not. mcupdt) then
           call gmem_free_inf(i10,fnm,snm,"i10")
         endif
      else
         call gmem_free_inf(i10,fnm,snm,"i10")
      endif
 
      return
 6010 format (1x,'core assignement'/1x,'i10, i20, i30 = ',3i10/1x,
     +        'last = ',i10)
 6020 format (/' integral evaluation complete at ',f8.2,' seconds')
 6030 format (/1x,20('*')/1x,'1-electron integrals'/1x,20('*')
     +        //' commence integral evaluation at ',f8.2,' seconds')
 6031 format (/,79('-')/,' hondo   ',a,/,79('-'))
      end
      subroutine sym1e(f,h,iso,move,nshels)
c
c     ----- symmetrize 1e-matrix in f to matrix h
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
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
      dimension f(*),h(*),iso(nshels,*)
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
      common/junk3/ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720)
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
      common/blkin/pxyz(4),t(35,35),mini,maxi,lit,minj,maxj,ljt,ntr
      dimension mi(48)
      if (nt .eq. 1) go to 600
      call vclr(h,1,nx)
c
c     ----- find a block (i,j)
c
      do 520 ii = 1,nshell
      do 140 itr = 1,nt
      ish = iso(ii,itr)
      if (ish .gt. ii) go to 520
  140 mi(itr) = ish
      lit = ktype(ii)
      mini = kmin(ii)
      maxi = kmax(ii)
      loci = kloc(ii)-mini
      do 500 jj = 1,ii
      do 200 itr = 1,nt
      jsh = iso(jj,itr)
      if (jsh .gt. ii) go to 500
      ish = mi(itr)
      if (ish .ge. jsh) go to 180
      n = ish
      ish = jsh
      jsh = n
  180 if (ish .eq. ii .and. jsh .gt. jj) go to 500
  200 continue
      ljt = ktype(jj)
      minj = kmin(jj)
      maxj = kmax(jj)
      locj = kloc(jj)-minj
      oiandj = ii .eq. jj
c
c     ----- find the equivalent blocks -----
c     ----- transfer equivalent block into t-matrix
c     ----- compute (r) t (r)
c     ----- put the result back into the (i,j) block of the h-matrix
c
      jmax = maxj
      do 300 itr = 1,nt
      ntr = itr
      kk = mi(itr)
      ll = iso(jj,itr)
      lock = kloc(kk)-kmin(kk)
      locl = kloc(ll)-kmin(ll)
      do 260 k = mini,maxi
      lck = lock+k
      if (oiandj) jmax = k
      do 260 l = minj,jmax
      kl = iky(max(lck,locl+l))+min(lck,locl+l)
      t(k,l) = f(kl)
      if (oiandj) t(l,k) = f(kl)
  260 continue
      call rhrrhr
      do 280 i = mini,maxi
      lci = iky(loci+i)+locj
      if (oiandj) jmax = i
      do 280 j = minj,jmax
      ij = lci+j
  280 h(ij) = h(ij)+t(i,j)
c
c     ----- for each block (k,l) equivalent to (i,j)
c     ----- find the transformation that maps (k,l) into (i,j)
c     ----- compute (r) t (r)
c     ----- put the result back into the (k,l) block of the h-matrix
c
  300 continue
      do 480 itr = 1,nt
      if (itr .eq. 1) go to 480
      kk = mi(itr)
      ll = iso(jj,itr)
      if (kk .ge. ll) go to 320
      k = ll
      l = kk
      go to 340
  320 k = kk
      l = ll
  340 if (k .eq. ii .and. l .eq. jj) go to 480
      ntr = itr+1
      if (ntr .gt. nt) go to 400
      do 380 it = ntr,nt
      i = mi(it)
      j = iso(jj,it)
      if (i .ge. j) go to 360
      ij = i
      i = j
      j = ij
  360 if (i .eq. k .and. j .eq. l) go to 480
  380 continue
  400 continue
      ntr = invt(itr)
      do 420 i = mini,maxi
      lci = iky(loci+i)+locj
      if (oiandj) jmax = i
      do 420 j = minj,jmax
      t(i,j) = h(lci+j)
      if (oiandj) t(j,i) = h(lci+j)
  420 continue
      call rhrrhr
      lock = kloc(kk)-kmin(kk)
      locl = kloc(ll)-kmin(ll)
      do 460 k = mini,maxi
      lck = lock+k
      if (oiandj) jmax = k
      do 460 l = minj,jmax
      kl = iky(max(lck,locl+l))+min(lck,locl+l)
  460 h(kl) = t(k,l)
  480 continue
  500 continue
  520 continue
      dum = 1.0d0/ dble(nt)
      call vsmul(h,1,dum,f,1,nx)
 600  continue
      if(move.ne.0) then
        call dcopy(nx,f,1,h,1)
      endif
      return
      end
      subroutine rhrrhr
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
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
      common/junk3/ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720)
      common/blkin/pxyz(4),t(35,35),mink,maxk,lkt,minl,maxl,llt,ntr
      dimension v(35)
c
c     ----- right multiply  t  by  r,
c           result back in  t
c
      data dzero /0.0d0/
c
      if (llt.eq.1) then
       go to 260
      else if(llt.eq.2) then
c
c     ----- p shell
c
       np = 3*(ntr-1)- 1
       do 240 k = mink,maxk
       do 220 l = 2,4
       dum = dzero
       do 200 n = 2,4
  200  dum = dum+t(k,n)*ptr(n-1,np+l)
  220  v(l) = dum
       do 240 l = 2,4
  240  t(k,l) = v(l)
c
      else if(llt.eq.3) then
c
c     ----- d shell
c
       nd = 6*(ntr-1)- 4
       do 160 k = mink,maxk
       do 140 l = 5,10
       dum = dzero
       do 120 n = 5,10
  120  dum = dum+t(k,n)*dtr(n-4,nd+l)
  140  v(l) = dum
       do 160 l = 5,10
  160  t(k,l) = v(l)
c
      elseif(llt.eq.4) then
c
c     ----- f shell
c
       nf=10*(ntr-1)-10
       do 60 k=mink,maxk
       do  80 l=11,20
       dum = dzero
       do  70 n=11,20
   70  dum=dum+t(k,n)*ftr(n-10,nf+l)
   80  v(l)=dum
       do 60 l=11,20
   60  t(k,l)=v(l)
c
       else
c
c     ----- g shell
c
      ng=15*(ntr-1)-20
      do k=mink,maxk
        do l=21,35
        dum= dzero
          do  n=21,35
           dum=dum+t(k,n)*gtr(n-20,ng+l)
          enddo
        v(l)=dum
        enddo
        do  l=21,35
          t(k,l)=v(l)
        enddo
      enddo
c
      endif
c
c     ----- left multiply  t  by r
c           result back in  t
c
260   if (lkt.eq.1) then
       return
      else if(lkt.eq.2) then
c
c     ----- p shell
c
       np = 3*(ntr-1)- 1
       do 420 l = minl,maxl
       do 400 k = 2,4
       dum = dzero
       do 380 n = 2,4
  380  dum = dum+ptr(n-1,np+k)*t(n,l)
  400  v(k) = dum
       do 420 k = 2,4
  420  t(k,l) = v(k)
c
      else if(lkt.eq.3) then
c
c     ----- d shell
c
       nd = 6*(ntr-1)-4
       do 340 l = minl,maxl
       do 320 k = 5,10
       dum = dzero
       do 300 n = 5,10
  300  dum = dum+dtr(n-4,nd+k)*t(n,l)
  320  v(k) = dum
       do 340 k = 5,10
  340  t(k,l) = v(k)
c
      else if(lkt.eq.4) then
c
c     ----- f shell
c
       nf=10*(ntr-1)-10
       do 530 l=minl,maxl
       do 520 k=11,20
       dum = dzero
       do 510 n=11,20
  510  dum=dum+ftr(n-10,nf+k)*t(n,l)
  520  v(k)=dum
       do 530 k=11,20
  530  t(k,l)=v(k)
c
      else
c
c     ----- g shell
c  
      ng=15*(ntr-1)-20
      do l=minl,maxl
        do k=21,35
        dum=dzero
          do n=21,35
            dum=dum+gtr(n-20,ng+k)*t(n,l)
          enddo
        v(k)=dum
        enddo
        do k=21,35
          t(k,l)=v(k)
        enddo
      enddo
c
      endif
c
      return
      end
      subroutine omitch(qzan,natom,nona)
      implicit none
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
c---------
c  omtchg - omit 1e- integrals involving the specified charges
c
c   integrals of the form  < bf(a) | q(b) | bf (a) > are deleted
c
c   where a and b are centres indicated as follows
c
c  ztagcc - list of centre names a
c  ztagcq - for each a in ztagcc, the list of centres b
c
c  omtslf - omit self-energy of point charge array (applied to all
c           atoms with centre label of form bq* )
c      CHECK
c      parameter (mxombq=20,mxomat=50)
c---------

c
c obqf  - special treatment of forces .. skip force contributions
c         on non-bq centres
c         for use in qm/mm optimisations of surrounding
c         assumption is that bq's in this case don't have basis
c         functions or ecps
c
c omitted charge specifications
c

      integer mxombq, mxomat
      parameter (mxombq=50,mxomat=20)
      character*8 ztagcc, ztagcq
      logical omtchg, omtslf, obqf
      common /chgcc/ ztagcc(mxomat),ztagcq(mxomat,mxombq),omtchg,
     &     omtslf,obqf
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer natom, nona
      real*8  qzan(natom,natom)
c
      integer i,i2,j,j2,icount
c
      integer mxprint
      parameter(mxprint=10)
      integer  ipr(mxprint)
c
*****
***** extract mapping for omitted charges (nona>0)
*****
      nona = 0
      do 70 i=1,nat
      do 70 j=1,nat
         qzan(i,j) = czan(j)
 70   continue
 
      do 10 i =1,mxombq
         if(ztagcc(i).ne.' ') then
            nona = nona + 1
         else
            go to 20
         endif
 10   continue
      if(nona.le.0)call caserr('invalid number of a-centres in omitc')
c
 20   continue
c
      write(iwr,6000)
6000  format(/
     + 10x,'#######################################'/
     + 10x,'Centre-Charge Interaction Specification'/
     + 10x,'#######################################'/,
     +  '   Atom       Omitted charges')
      do 30 i = 1,nona
         do i2 = 1,nat
            if(ztagcc(i) .eq. zaname(i2))then
               icount = 0
               do 40 j = 1,mxomat
                  if(ztagcq(i,j).ne.' ') then
                     do j2 = 1,nat
                        if(ztagcq(i,j) .eq. zaname(j2))then
                           icount = icount + 1
                           qzan(i2,j2) = 0.0d0
                           if(icount.le.mxprint)then
                              ipr(icount)=j2
                           endif
                        endif
                     enddo
                  else
                     go to 50
                  endif
 40            continue

 50            if(icount.eq.0) call 
     &              caserr('invalid no. of b-centres in omitc')
            
               if(icount.gt.mxprint)then
                  write(iwr,6011) i2,zaname(i2),(ipr(j),zaname(ipr(j)),
     &                 j=1,mxprint)
 6011             format(1x,i4,'(',a8,')',5x,10(1x,i4,1x,'(',a8,')'),
     &                 '...')
               else
                  write(iwr,6010) i2,zaname(i2),(ipr(j),zaname(ipr(j)),
     &                 j=1,icount)
 6010             format(1x,i4,'(',a8,')',5x,10(1x,i4,1x,'(',a8,')'))
               endif
c
            endif
         enddo
 30   continue
c
      if(nat.le.20)then
         write(iwr,6019) (zaname(i),i=1,nat)
 6019    format(9x,20a8)
         do 80 i =1,nat
            write(iwr,6020) zaname(i),(qzan(i,j),j=1,nat)
 6020       format(1x,a8, 20f8.3)
 80      continue
      else
         write(iwr,*)'Charge matrix print supressed for nat>12'
      endif
c
      call prsq(qzan,nat,nat,nat)
c
      return
      end
      subroutine stvint
c
c     ----- gauss-hermite quadrature using minimum point formula -----
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      common/junk/ssgg(450)
     *,pint,qint,rint,t,p0,q0,r0,pi,qi,ri,pj,qj,rj,ni,nj
      common/hermit/h(45)
      common/wermit/w(45)
      dimension min(7),max(7)
      data min /1,2,4,7,11,16,22/
      data max /1,3,6,10,15,21,28/
      data dzero /0.0d0/
      pint = dzero
      qint = dzero
      rint = dzero
      npts = (ni+nj-2)/2 + 1
      imin = min(npts)
      imax = max(npts)
      do 160 i = imin , imax
         dum = w(i)
         px = dum
         py = dum
         pz = dum
         dum = h(i)*t
         ptx = dum + p0
         pty = dum + q0
         ptz = dum + r0
         ax = ptx - pi
         ay = pty - qi
         az = ptz - ri
         bx = ptx - pj
         by = pty - qj
         bz = ptz - rj
         go to (80,70,60,50,40,30,20) , ni
 20      px = px*ax
         py = py*ay
         pz = pz*az
 30      px = px*ax
         py = py*ay
         pz = pz*az
 40      px = px*ax
         py = py*ay
         pz = pz*az
 50      px = px*ax
         py = py*ay
         pz = pz*az
 60      px = px*ax
         py = py*ay
         pz = pz*az
 70      px = px*ax
         py = py*ay
         pz = pz*az
 80      go to (150,140,130,120,110,100,90) , nj
 90      px = px*bx
         py = py*by
         pz = pz*bz
 100     px = px*bx
         py = py*by
         pz = pz*bz
 110     px = px*bx
         py = py*by
         pz = pz*bz
 120     px = px*bx
         py = py*by
         pz = pz*bz
 130     px = px*bx
         py = py*by
         pz = pz*bz
 140     px = px*bx
         py = py*by
         pz = pz*bz
 150     pint = pint + px
         qint = qint + py
         rint = rint + pz
 160  continue
      return
      end
      subroutine coulmb(clints,gout)
c
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
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
      dimension clints(*),gout(*)
      dimension ibpop(4,4),ibr(4,4)
c
c
      integer ib
      common /flip70/ ib(4,3)
c
c
      integer igt, jgt, kgt, lgt
      integer ibt, jbt, kbt, lbt, lenddd
      common /flips/ igt(3),jgt(3),kgt(3),lgt(3),
     +               ibt(mxprms),jbt(mxprms),kbt(mxprms),lbt(mxprms),
     +               lenddd
c
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
c
      integer ik, ijx, ijy, ijz, klx, kly, klz
      integer ijgt, klgt
      common /indez/ ijgt(225),ijx(225),ijy(225),ijz(225),
     +       ik(225),klgt(225),klx(225),kly(225),klz(225)
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
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
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
      integer ishell, jshell, kshell, lshell
      integer inew, jnew, knew, lnew
      common/shlg70/ishell,jshell,kshell,lshell,inew,jnew,knew,lnew
c
c
      real*8 pito52, pidiv4, root3, root5, root53, root7
      common /picon/ pito52,pidiv4,root3,root5,root53,root7
c
c
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
c
      logical ospbas, onocnt, ochek, opdbas, opfbas, opgbas
      integer kad, istd
      common /ijlab/ kad(mxshel),ospbas,onocnt,istd,ochek,
     +               opdbas,opfbas,opgbas
c
c
c     dlnmxd: the maximum value of the logarithm of the reduced
c             density matrix.
c     dlnmxs: the maximum value of the logarithm of the Schwarz 
c             integrals.
c     dlntol: minus the logarithm of the direct SCF integral cutoff.
c
      real*8 dlnmxd, dlnmxs, dlntol, tolitr, deltol, delfac
      integer intcut, intmag
      integer ibl171, lentri, itrtol, iexch, iof171, m171t
      integer nhamd, nsheld
      logical odscf, oimag, odnew, opref, odelta
      common/cslosc/dlnmxd,dlnmxs,dlntol,tolitr(3),deltol,delfac,
     + odscf,
     + intcut(10),oimag,intmag(1060),
     + ibl171,lentri,itrtol(3),iexch,odnew,opref,iof171(6),m171t,
     + nhamd, nsheld,odelta
c
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
c
      data m25/25/
c     data m22,mword1/22,6000/
c
      data ibpop/0,0,0,0,64,16,4,1,128,32,8,2,192,48,12,3/
      data ibr/64,16,4,1,216,36,6,1,1000,100,10,1,
     +         3375,225,15,1 /
      data done/1.0d0/,two/2.0d0/,twopt5/2.5d0/,four/4.0d0/
      data dzero/0.0d0/
c
c     ----- compute all coulomb integrals (ii,jj//ii,jj) -----
c     the largest coulomb integral from each integral block
c     is picked out and saved, so that the schwarz inequality
c     can be used later to avoid entire integral blocks.
c
c     instead of saving SQRT( (ii,jj//ii,jj) ) we now save
c     its logarithm 0.5*LOG( (ii,jj//ii,jj) ) because we need that
c     for the direct SCF.
c
      call spchck
      outi = nprint.ne.-5.and.oprint(55)
      nclint = ikyp(nshell)
      lensch = lensec(nclint)
      tim0 = cpulft(1)
c
c === save tolerances set in debut and increase accuracy
c
      tolsv = tol
      tol = 75.0d00
      cutsv = cutoff
      cutoff = 1.0d0/(10.0d0**20)
c
      qq4   = done
      nint  = 0
      call vclr(clints,1,nclint)
      if (intg76.ne.0) then
c
c     ----- loop over all shell blocks -----
c     ----- first over rotated axis integrals
c
         pidiv4 = datan(done)
         pi = four*pidiv4
         pito52 = two*pi**twopt5
         call sinset
         call filmax
         do 60 ishell = 1, nshell
            iiii = iky(ishell)
            if (kad(ishell).ge.0) then
               do 50 jshell = 1, ishell
                  ijij = iiii + jshell
                  if (kad(jshell).ge.0) then
c
c     use pople code for any pure sp integral blocks,
c     use hondo rys polynomial code for other blocks
c
                     kshell = ishell
                     lshell = jshell
                     call genr70(gout,1,.false.)
c
c     ----- pick out largest coulomb integral for this block -----
c
                     vmax = dzero
                     mini = kmin(ishell)
                     minj = kmin(jshell)
                     maxi = kmax(ishell)
                     jmax = kmax(jshell)
                     oiandj = ishell.eq.jshell
                     ibb = ib(1,1)
                     jbb = ib(2,1)
                     kbb = ib(3,1)
                     lbb = ib(4,1)
                     ijn = 0
                     do 40 i = mini, maxi
                        if (oiandj) jmax = i
                        do 30 j = minj, jmax
                           nn = ibpop(ibb,i) + ibpop(jbb,j)
     +                        + ibpop(kbb,i) + ibpop(lbb,j) + 1
                           val = gout(nn)
                           if (val.gt.dzero) nint = nint+1
                           if (val.gt.vmax) vmax = val
 30                     continue
 40                  continue
                     clints(ijij) = vmax
                  end if
 50            continue
            end if
 60      continue
         if (ospbas) go to 110
      end if
c
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
      do 300 loop = 1 , 3
        igt(loop) = ibr(1,ii)
        jgt(loop) = ibr(2,ii)
        kgt(loop) = ibr(3,ii)
        lgt(loop) = ibr(4,ii)
 300  continue
      qq4   = done
c
      do 100 ishell = 1, nshell
         iiii = iky(ishell)
         kadi = kad(ishell)
         do 90 jshell = 1, ishell
            ijij = iiii + jshell
            if ((kadi+kad(jshell)).lt.0) then
c
c     use hondo rys polynomial code for other blocks
c
               call shells(gout,1,ishell,jshell,ishell,jshell,1)
               call ijprim
               if (nij.ne.0) then 
                  call shells(gout,2,ishell,jshell,ishell,jshell,1)
                  call genral(gout)
c
c     ----- pick out largest coulomb integral for this block -----
c
                  vmax = dzero
                  mini = kmin(ishell)
                  minj = kmin(jshell)
                  maxi = kmax(ishell)
                  jmax = kmax(jshell)
                  oiandj = ishell.eq.jshell
                  ijn = 0
                  do 80 i = mini, maxi
                     if (oiandj) jmax = i
                     do 70 j = minj, jmax
                        ijn = ijn+1
                        nn = ijgt(ijn) + klgt(ijn)
                        val = gout(nn)
                        if (val.gt.dzero) nint = nint+1
                        if (val.gt.vmax) vmax = val
 70                  continue
 80               continue
                  clints(ijij) = vmax
               end if
            end if
 90      continue
 100  continue
c
 110  if (outv) then
         write(iwr,6020)
         call prtri(clints,nshell)
      end if
c
      dlnmxs = -1.0d50
      do 120 i = 1, nclint
         clints(i) = 0.5d0*dlog(dmax1(1.0d-60,clints(i)))
         dlnmxs = dmax1(dlnmxs,clints(i))
 120  continue
c
      tcoul = cpulft(1)-tim0
      if (outi) write(iwr,6030) nint, tcoul
c
      call secput(isect(421),m25,lensch,iblock)
      call wrt3(clints,nclint,iblock,idaf)
c
c === now reset tolerances
c
      tol = tolsv
      cutoff = cutsv
c
      call spchck
c
      return
 6020 format (/20x,'*****************************'/
     +         20x,'max coulomb integral in shell'/
     +         20x,'*****************************'/)
 6030 format(1x,
     + 'schwarz inequality overhead:',i10,' integrals, t=',
     *       f8.2,' seconds')
      end
      subroutine debut(zscftp)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
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
c     dlnmxd: the maximum value of the logarithm of the reduced
c             density matrix.
c     dlnmxs: the maximum value of the logarithm of the Schwarz 
c             integrals.
c     dlntol: minus the logarithm of the direct SCF integral cutoff.
c
      real*8 dlnmxd, dlnmxs, dlntol, tolitr, deltol, delfac
      integer intcut, intmag
      integer ibl171, lentri, itrtol, iexch, iof171, m171t
      integer nhamd, nsheld
      logical odscf, oimag, odnew, opref, odelta
      common/cslosc/dlnmxd,dlnmxs,dlntol,tolitr(3),deltol,delfac,
     + odscf,
     + intcut(10),oimag,intmag(1060),
     + ibl171,lentri,itrtol(3),iexch,odnew,opref,iof171(6),m171t,
     + nhamd, nsheld,odelta
c
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
      common/scfopt/maxit(4),accdi(2),odiis(4),dmpcut(7),iter
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
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
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
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer n2file, n2tape, n2blk, n2last
      integer n4file, n4tape, n4blk, n4last
      integer n6file, n6tape, n6blk, n6last
      integer n5file, n5tape, n5blk, n5last
      integer n9file, n9tape, n9blk, n9last
      integer ntfile, nttape, ntblk, ntlast
      integer n1file, n1tape, n1blk, n1last
      integer n11fil, n11tap, n11bl, n11lst
      integer n12fil, n12tap, n12bl, n12lst
      integer n13fil, n13tap, n13bl, n13lst
      integer maxbl, minbl, maxset
      common/gms_files/
     +      n2file,n2tape(20),n2blk(20),n2last(20),
     +      n4file,n4tape(20),n4blk(20),n4last(20),
     +      n6file,n6tape(20),n6blk(20),n6last(20),
     +      n5file,n5tape(20),n5blk(20),n5last(20),
     +      n9file,n9tape(20),n9blk(20),n9last(20),
     +      ntfile,nttape(20),ntblk(20),ntlast(20),
     +      n1file,n1tape(20),n1blk(20),n1last(20),
     +      n11fil,n11tap(20),n11bl(20),n11lst(20),
     +      n12fil,n12tap(20),n12bl(20),n12lst(20),
     +      n13fil,n13tap(20),n13bl(20),n13lst(20),
     +       maxbl(maxlfn),minbl(maxlfn),maxset(maxlfn)
c
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
c
      real*8 p, q, r
      common /bshell/ p(mxshel),q(mxshel),r(mxshel)
c
c     logical opk,opank,oskips
c     common/pkfil/opk,opank,oskips
      logical opk,opank
      common/pkfil/opk,opank
      character*10 charwall
      data zuhf,zgvb/'uhf','gvb'/
      data zgrhf/'grhf'/
      data done,ten,e /1.0d0,1.0d+01,2.30258d0/
c
      outv = oprint(59)
      if (nprint.eq.-5) outv = .false.
      out = nprint.eq.4
      opk = nopk.ne.1
      opank = (opk .and. zscftp.eq.zuhf) .or. (opk .and. na.ne.nb) .or.
     +        (opk .and. zscftp.eq.zgvb) .or.
     +        (opk .and. zscftp.eq.zgrhf)
      if (opk .and. nopk.eq.(-1)) opank = .true.
c
c..    for denscf adjust nopk
c
      if (opank .and. nopk.eq.0) nopk = -1
c
      if (.not.odscf) then
         if (nprint.ne.-5) then
            write (iwr,6010)
            if (opk) then
               write (iwr,6020)
               if (opank) then
                  write (iwr,6040)
               else
                  write (iwr,6050)
               end if
            else
               write (iwr,6030)
            end if
            write (iwr,6080) nintmx
         end if
      end if
      omaxb = .false.
c
c ----- setup dumping on time
c
      cutoff = dumtim*1.5d0
      safety = cpulft(0)
      if (safety.gt.cutoff) then
         safety = safety - dumtim
      else
         safety = safe
      end if
      timlst = cpulft(0) - safety
      icut0 = iabs(icut)
      if (icut0.eq.0) icut0 = 9
      cutoff = done/(ten**icut0)
      if (itol.eq.0) itol = 20
      tol = e*itol
      if (irest.ge.1) then
c
c     -----                  irest = 1                   -----
c
c     check on input from 1-e ints (nindmx = -1)
c
         if (nindmx.ne.-1) then
c
c     ----- position the integral file for a restart job -----
c
            icount = 1
            ic4 = 1
            if (m2blk(mfilep).lt.0) then
               iblkmp = n2blk(mfilep)
               m2tape(mfilep) = n2tape(mfilep)
               m2blk(mfilep) = iblkmp
               m2last(mfilep) = -1
            else
               iblkmp = m2last(mfilep)
            end if
            mainp = m2tape(mfilep)
            mblp = iblkmp - n2last(mfilep)
            if (nprint.ne.-5 .and. (.not.odscf)) write (iwr,6090)
     +          yed(mainp) , m2blk(mfilep) , ist , jst , kst , lst ,
     +          iblkmp
            do 20 i = 1 , nshell
               icc = katom(i)
               p(i) = c(1,icc)
               q(i) = c(2,icc)
               r(i) = c(3,icc)
 20         continue
            go to 50
         end if
      end if
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
      do 40 i = 1 , nshell
         icc = katom(i)
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
 50   if (.not.odscf) call search(iblkmp,mainp)
      cpu = cpulft(1)
      isti = ist
      jsti = jst
      ksti = kst
      lsti = lst
      lastb = iblkmp
      lastu = mfilep
      if (nprint.ne.-5) then
         write (iwr,6060) cpu,charwall()
         if (odscf) then
            if (iter.le.0 .and. outv) write (iwr,6070)
         else
            if (outv) write (iwr,6070)
         end if
      end if
      return
 6010 format (//1x,20('*')/1x,'2-electron integrals'/1x,20('*')/)
 6020 format (/' integrals are in a supermatrix form')
 6030 format (/
     +    ' integrals are not in a supermatrix form : p-k option is off'
     +    )
 6040 format (' generate -p- and -k- supermatrices')
 6050 format (' generate -p- supermatrix only')
 6060 format (/' commence 2-electron integral evaluation at ',f12.2,
     +        ' seconds',a10,' wall')
 6070 format (/1x,'ist',2x,'jst',2x,'kst',2x,'lst',7x,'nrec',3x,
     +        'intloc',5x,'del(t)',5x,'time'/1x,58('-'))
 6080 format (' number of integrals per block = ',i3)
 6090 format (/' integrals are on mainfile section ',a4,
     +        ' starting at block ',i6//' starting shells are : ',
     +        4i5//' restarting at block ',i6/)
      end
      subroutine genral(gout)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
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
      parameter (mxp2 = mxprms * mxprms)
c
c     dlnmxd: the maximum value of the logarithm of the reduced
c             density matrix.
c     dlnmxs: the maximum value of the logarithm of the Schwarz 
c             integrals.
c     dlntol: minus the logarithm of the direct SCF integral cutoff.
c
      real*8 dlnmxd, dlnmxs, dlntol, tolitr, deltol, delfac
      integer intcut, intmag
      integer ibl171, lentri, itrtol, iexch, iof171, m171t
      integer nhamd, nsheld
      logical odscf, oimag, odnew, opref, odelta
      common/cslosc/dlnmxd,dlnmxs,dlntol,tolitr(3),deltol,delfac,
     + odscf,
     + intcut(10),oimag,intmag(1060),
     + ibl171,lentri,itrtol(3),iexch,odnew,opref,iof171(6),m171t,
     + nhamd, nsheld,odelta
c
c
      real*8 ag, csa, cpa, cda, cfa, cga
      real*8 bg, csb, cpb, cdb, cfb, cgb
      real*8 cgg, csc, cpc, cdc, cfc, cgc
      real*8 dg, csd, cpd ,cdd, cfd, cgd
      real*8 pi, qi, ri, pj, qj, rj, rri, pk, qk, rk, pl, ql, rl, rrk
      integer nga, ngb, ngc, ngd
      common /shlinf/ ag(mxprms),csa(mxprms),cpa(mxprms),cda(mxprms),
     +                cfa(mxprms),cga(mxprms),
     +                bg(mxprms),csb(mxprms),cpb(mxprms),cdb(mxprms),
     +                cfb(mxprms),cgb(mxprms),
     +                cgg(mxprms),csc(mxprms),cpc(mxprms),cdc(mxprms),
     +                cfc(mxprms),cgc(mxprms),
     +                dg(mxprms),csd(mxprms),cpd(mxprms),cdd(mxprms),
     +                cfd(mxprms),cgd(mxprms),
     +                pi,qi,ri,pj,qj,rj,rri,pk,qk,rk,pl,ql,rl,rrk,
     +                nga,ngb,ngc,ngd
c
      common/junk/pin(5625), qin(5625), rin(5625),
     + aa(mxp2),r(mxp2),p1(mxp2),q1(mxp2),r1(mxp2),dd(16*mxp2),
     + ijd(225)
c
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
c
      real*8 pp, u, w, dji
      integer nroots, nroot2
      common /root/ pp,u(12),w(12),nroots,nroot2,dji(60)
c
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
c
      logical oianj, okanl, oident, omisc, oham, opdipd, omp2
      integer ipos1, ipos2
      common /misc/ oianj,okanl,oident,omisc,
     +            oham,opdipd,omp2,ipos1,ipos2
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
      integer in, kn, ni, nj, nk, nl, nmax, mmax
      integer ij1, ij2, kl1, kl2
      real*8 bp01, b00, b10, pcp00, pc00, qcp00, qc00, rcp00, rc00
      real*8 f00, dxij, dyij, dzij, dxkl, dykl, dzkl
      common /setint/ in(9),kn(9),ni,nj,nk,nl,nmax,mmax,
     +                bp01,b00,b10,pcp00,pc00,qcp00,qc00,rcp00,rc00,
     +                f00,dxij,dyij,dzij,dxkl,dykl,dzkl,
     +                ij1,ij2,kl1,kl2
c
c
      real*8 dij, dkl
      common /denss/ dij(225),dkl(225)
c
      dimension gout(*),in1(9)
c
      data ijn1,ijn2,kln1,kln2 /125,25,5,1/
      data sqrt3 /1.73205080756888d0/
      data sqrt5 /2.23606797749979d+00/
      data sqrt7 /2.64575131106459d+00/
      data pi252 /34.986836655250d0/
      data pt5,done /0.5d0,1.0d0/
c
      if (ijkl.eq.1) then
c
c     ----- (s,s//s,s) -----
c
         call s0000(gout)
      else if (ij.eq.1) then
c
c     ----- (s,s//k,l) -----
c
         call sskl(gout)
      else if (kl.eq.1) then
         call ijss(gout)
      else
         factor = pi252*qq4
         onorm = normf.ne.1 .or. normp.ne.1
c
c     ----- select expansion center for -xyz- integrals -----
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
            do 260 lg = 1 , lgmax
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
 70                  if (onorm) dum1 = dum1*sqrt3
                     go to 110
 80                  dum1 = cfk*b1
                     go to 110
 90                  if (onorm) dum1 = dum1*sqrt5
                     go to 110
 100                 if (onorm) dum1 = dum1*sqrt3
                     go to 110
 92                  dum1 = cgk*b1
                     go to 110
 94                  if (onorm) dum1 = dum1*sqrt7
                     go to 110
 96                  if (onorm) dum1 = dum1*sqrt5/sqrt3
                     go to 110
 98                  if (onorm) dum1 = dum1*sqrt3
 110                 if (okanl) max = k
                     do 200 l = minl , max
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
 150                    if (onorm) dum2 = dum2 * sqrt3
                        go to 190
 160                    dum2 = dum1*cfl
                        if (odoub) dum2 = dum2 + dum2
                        go to 190
 170                    if (onorm) dum2 = dum2 * sqrt5
                        go to 190
 180                    if (onorm) dum2 = dum2 * sqrt3
                        go to 190
 182                    dum2 = dum1*cgl
                        if (odoub) dum2 = dum2+dum2
                        go to 190
 184                   if (onorm) dum2 = dum2 * sqrt7
                        go to 190
 186                   if (onorm) dum2 = dum2 * sqrt5/sqrt3
                        go to 190
 188                   if (onorm) dum2 = dum2 * sqrt3
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
                        if (nroots.le.3) call rt123
                        if (nroots.eq.4) call roots4
                        if (nroots.eq.5) call roots5
                        if (nroots.ge.6) call rootss
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
                           call xyzint
                           mm = mm + 625
 240                    continue
c
c     ----- form (i,j//k,l) integrals over functions
c
                        call spdint(gout)
                     end if
                     nn = nn + 16
 250              continue
               end if
 260        continue
 270     continue
      end if
      return
      end
      subroutine ijprim
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
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
      parameter (mxp2 = mxprms * mxprms)
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
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
c
      real*8 ag, csa, cpa, cda, cfa, cga
      real*8 bg, csb, cpb, cdb, cfb, cgb
      real*8 cgg, csc, cpc, cdc, cfc, cgc
      real*8 dg, csd, cpd ,cdd, cfd, cgd
      real*8 pi, qi, ri, pj, qj, rj, rri, pk, qk, rk, pl, ql, rl, rrk
      integer nga, ngb, ngc, ngd
      common /shlinf/ ag(mxprms),csa(mxprms),cpa(mxprms),cda(mxprms),
     +                cfa(mxprms),cga(mxprms),
     +                bg(mxprms),csb(mxprms),cpb(mxprms),cdb(mxprms),
     +                cfb(mxprms),cgb(mxprms),
     +                cgg(mxprms),csc(mxprms),cpc(mxprms),cdc(mxprms),
     +                cfc(mxprms),cgc(mxprms),
     +                dg(mxprms),csd(mxprms),cpd(mxprms),cdd(mxprms),
     +                cfd(mxprms),cgd(mxprms),
     +                pi,qi,ri,pj,qj,rj,rri,pk,qk,rk,pl,ql,rl,rrk,
     +                nga,ngb,ngc,ngd
c
c
      logical oianj, okanl, oident, omisc, oham, opdipd, omp2
      integer ipos1, ipos2
      common /misc/ oianj,okanl,oident,omisc,
     +            oham,opdipd,omp2,ipos1,ipos2
c
c
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
      common/junk/cxyz(3,5625),
     +     a(mxp2),r(mxp2),p1(mxp2),q1(mxp2),r1(mxp2),dij(16*mxp2),
     +     ijd(225)
c
      data sqrt3 /1.73205080756888d0/
      data sqrt5 /2.23606797749979d+00/
      data sqrt7 /2.64575131106459d+00/
      data done /1.0d0/
c
      if (ij.eq.1) then
c
c     ----- (s,s// -----
c
         call ssprim
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
      subroutine ijss(gout)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
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
      parameter (mxp2 = mxprms * mxprms)
c
      real*8 ag, csa, cpa, cda, cfa, cga
      real*8 bg, csb, cpb, cdb, cfb, cgb
      real*8 cgg, csc, cpc, cdc, cfc, cgc
      real*8 dg, csd, cpd ,cdd, cfd, cgd
      real*8 pi, qi, ri, pj, qj, rj, rri, pk, qk, rk, pl, ql, rl, rrk
      integer nga, ngb, ngc, ngd
      common /shlinf/ ag(mxprms),csa(mxprms),cpa(mxprms),cda(mxprms),
     +                cfa(mxprms),cga(mxprms),
     +                bg(mxprms),csb(mxprms),cpb(mxprms),cdb(mxprms),
     +                cfb(mxprms),cgb(mxprms),
     +                cgg(mxprms),csc(mxprms),cpc(mxprms),cdc(mxprms),
     +                cfc(mxprms),cgc(mxprms),
     +                dg(mxprms),csd(mxprms),cpd(mxprms),cdd(mxprms),
     +                cfd(mxprms),cgd(mxprms),
     +                pi,qi,ri,pj,qj,rj,rri,pk,qk,rk,pl,ql,rl,rrk,
     +                nga,ngb,ngc,ngd
c
      common/junk/pint(5625),qint(5625),rint(5625),
     + aa(mxp2),r(mxp2),p1(mxp2),q1(mxp2),r1(mxp2),dd(16*mxp2),
     + ijd(225)
c
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
c
      real*8 pp, u, w, dji
      integer nroots, nroot2
      common /root/ pp,u(12),w(12),nroots,nroot2,dji(60)
c
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
c
      logical oianj, okanl, oident, omisc, oham, opdipd, omp2
      integer ipos1, ipos2
      common /misc/ oianj,okanl,oident,omisc,
     +            oham,opdipd,omp2,ipos1,ipos2
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
      integer ik, ijx, ijy, ijz, klx, kly, klz
      integer ijgt, klgt
      common /indez/ ijgt(225),ijx(225),ijy(225),ijz(225),
     +       ik(225),klgt(225),klx(225),kly(225),klz(225)
c
c
      integer in, kn, ni, nj, nk, nl, nmax, mmax
      integer ij1, ij2, kl1, kl2
      real*8 bp01, b00, b10, pcp00, pc00, qcp00, qc00, rcp00, rc00
      real*8 f00, dxij, dyij, dzij, dxkl, dykl, dzkl
      common /setint/ in(9),kn(9),ni,nj,nk,nl,nmax,mmax,
     +                bp01,b00,b10,pcp00,pc00,qcp00,qc00,rcp00,rc00,
     +                f00,dxij,dyij,dzij,dxkl,dykl,dzkl,
     +                ij1,ij2,kl1,kl2
c
c
      real*8 dij, dkl
      common /denss/ dij(225),dkl(225)
c
      dimension gout(*)
      data ijn1,ijn2 /125,25/
      data pi252 /34.986836655250d0/
      data dzero,pt5,done /0.0d0,0.5d0,1.0d0/
      factor = pi252*qq4
c
c     ----- select expansion center for -xyz- integrals -----
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
                        call rt123
                     else if (nroots.eq.4) then
                        call roots4
                     else
                        call roots5
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
      subroutine jkinta(zscftp,q,fock,fockb,exch,dens,densb,prefac,
     +                  rdmat,iso,gout,nshels)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *59 text
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
      parameter (mxp2 = mxprms * mxprms)
c
c     dlnmxd: the maximum value of the logarithm of the reduced
c             density matrix.
c     dlnmxs: the maximum value of the logarithm of the Schwarz 
c             integrals.
c     dlntol: minus the logarithm of the direct SCF integral cutoff.
c
      real*8 dlnmxd, dlnmxs, dlntol, tolitr, deltol, delfac
      integer intcut, intmag
      integer ibl171, lentri, itrtol, iexch, iof171, m171t
      integer nhamd, nsheld
      logical odscf, oimag, odnew, opref, odelta
      common/cslosc/dlnmxd,dlnmxs,dlntol,tolitr(3),deltol,delfac,
     + odscf,
     + intcut(10),oimag,intmag(1060),
     + ibl171,lentri,itrtol(3),iexch,odnew,opref,iof171(6),m171t,
     + nhamd, nsheld,odelta
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      common/scfopt/maxit(4),accdi(2),odiis(4),dmpcut(7),iter
c
      character*8 field
      character*24 version
      common /drfopt/ field,version
      logical oreact
      integer intdrf
      common /drfopti/ oreact,intdrf
c
c
c  control parameters for morokuma energy decomposition
c  analysis
c
      logical omorok, oifc, oife, oifh
      integer imode, ifrag, ifocc, imoro
      real*8 emorok
      common/morok1/emorok(0:4,2),imode,ifrag(maxorb),
     &     ifocc(maxorb),imoro, omorok, 
     &     oifc, oife, oifh
      character*20 fragname
      common/morokc/fragname(3)

c
c  settings
c
      integer ALL_BLOCKS
      parameter (ALL_BLOCKS=0)
      
      integer ESX
      parameter (ESX=1)

      integer ESX_PLX
      parameter (ESX_PLX=2)

      integer ESX_CTT
      parameter (ESX_CTT=3)

      integer ESX_EX
      parameter (ESX_EX=4)

      integer YES_K
      parameter (YES_K=1)

      integer NO_K
      parameter (NO_K=2)

c
c  field,intdrf:
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical odumpm, ooprnt, opun7, omax, ologf
      integer maxbrd, logf, iofsp
      common /utilc/ odumpm,ooprnt,opun7,omax,maxbrd,ologf,
     +               logf,iofsp(4)
c
c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
c
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
c
      logical ospbas, onocnt, ochek, opdbas, opfbas, opgbas
      integer kad, istd
      common /ijlab/ kad(mxshel),ospbas,onocnt,istd,ochek,
     +               opdbas,opfbas,opgbas
c
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
      common/junk/cxyz(3,5625),aaa(21*mxp2),ijaaa(225)
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
c     logical opk,opank,oskips
c     common/pkfil/opk,opank,oskips
      logical opk,opank
      common/pkfil/opk,opank
c
      integer igt, jgt, kgt, lgt
      integer ibt, jbt, kbt, lbt, lenddd
      common /flips/ igt(3),jgt(3),kgt(3),lgt(3),
     +               ibt(mxprms),jbt(mxprms),kbt(mxprms),lbt(mxprms),
     +               lenddd
c
      logical osortp,oschw
      common/sortpdat/
     + igbf,iklbf,iiptbf,ijbase,ngbf,lengbf,osortp,oschw
      integer icount_dlb, icount_slb
      common /par_dlb/ icount_dlb, icount_slb
c
c     ZORA common
c
      logical ozora,oscalz,onlyp,onlyd,onlyf,onlyg,opzora,
     1        op2zora,onlys,ocontr,oint_zora,o1e_zora,osmall_zora,
     2        opre_zora,oso,onoext_z,oioraz,oatint_z,oatscf_z,
     3        oscalatz
      real*8 cspeed,critat_z
      integer nbas_zora,ibls_z,iblv_z,nwv_z,ibiso_z,nwiso_z,ibcor_z,
     1        ibsx_z,icoul_z,ibsy_z,ibsz_z,ibsmin_z,ibsplus_z,nwmat_z,
     2        ibtrin_z,ibtrout_z,ibscalf_z,nwtrat_z,ibscale_z,ibcoul_z,
     3        niter_z,int_zora,num_ext,igauge_z,nwcor_z,isexch_z,
     4        nat_z,ibdens_z,nwdens_z,ibscalao_z,ibshat_z,is_z,
     5        irest_z,numdu_z,ibldu_z,iblock_z,lshift_z,ibl7ew_z
      common/zorac/ cspeed,critat_z,ozora,oscalz,
     1              onlyp,onlyd,onlyf,onlyg,opzora,
     1              op2zora,onlys,ocontr,nbas_zora,
     2              oint_zora,o1e_zora,icoul_z,
     3              ibls_z,iblv_z,nwv_z,
     3              ibiso_z,nwiso_z,ibcor_z,nwcor_z,
     4              ibsx_z,ibsy_z,ibsz_z,ibsmin_z,ibsplus_z,nwmat_z,
     5              ibtrin_z,ibtrout_z,ibscalf_z,nwtrat_z,
     6              ibscale_z,ibcoul_z,ibdens_z,nwdens_z,osmall_zora,
     7              niter_z,opre_zora,int_zora,num_ext,isexch_z,
     8              oso,onoext_z,is_z,igauge_z,oioraz,
     9              nat_z,ibscalao_z,oatint_z,oatscf_z,ibshat_z,
     1              oscalatz,irest_z,numdu_z,ibldu_z,iblock_z,lshift_z,
     2              ibl7ew_z
c....    icoul_z : 1 : Full coulomb
c...               2 : atomic coulomb (default)
c...               3 : no coulomb (just 1-electron ints)
c...               4 : small basis coulomb
c...     isexch_z  0 : no exchange in zora correction (default)
c...               1 : exchange in zora correction (like in dft)
c...     niter_z : # iterations between coulomb matrix update
c...     int_zora :internal basis gen. : 0 none, 1, automatic, 2 by hand
c
c...     oscalz  : scaling on or off (true or false)
c...     oioraz  : change metric to perform inf. order scaling 
c...               (Dyall/EvLenthe (not quite))
c...     oscalatz: atoms only scaled zora (together with get atom)
c...     only'x' : internal basis only op to 'x'-type functions
c...     op(2)zora: print flags (op2zora+opzora: prints *everything*)
c...     ocontr  : adds l+1 block to internal basis in contracted way 
c...     nbas_zora: size of internal basis
c...     oint_zora: in internal or external mode???
c...     o1e_zora: one electron integrals have been calculated in *this*
c...               internal basis
c...     osmall_zora: small swap for projected coulomb
c...     opre_zora: pre_zora has been called
c...     oso     : spin orbit calculation 
c...     onoext_z: no external functions copied to internal basis
c...     is_z : section for orbitals used for  coulomb matrix
c...     igauge_z : add so many s-marices to h-matrix
c...                the energy should change by igauge
c...     nat_z  : indicates (if ne 0) that starting density is used
c...              to calc. coulomb; gives # it. to get atomic density
c...     critat_z : accuracy of atomic energy for atomic zora
c...     oatint_z : indicates that integrals to be calculated are atomic
c...     oatscf_z : indicates we are doing a set of atomic SCFs
c...     ibshat_z : block on num8 to store blocked s/h
c...     for disk-adress and lengths (ib.., nw.. ) see subroutine atoms2
c...     irest_z  : .ne.0 indicates atomic zora corrections from dumpfile; 1 not on restart
c...     numdu_z,ibldu_z : if ne 0 other dumpfile to read zora corrections from
c...     iblock_z is the block number on ed3, where the zora corrections (section 425) reside
c...     ibl7ew is where the energy weighted density matrix resides
c
c
      real*8 ag, csa, cpa, cda, cfa, cga
      real*8 bg, csb, cpb, cdb, cfb, cgb
      real*8 cgg, csc, cpc, cdc, cfc, cgc
      real*8 dg, csd, cpd ,cdd, cfd, cgd
      real*8 pi, qi, ri, pj, qj, rj, rri, pk, qk, rk, pl, ql, rl, rrk
      integer nga, ngb, ngc, ngd
      common /shlinf/ ag(mxprms),csa(mxprms),cpa(mxprms),cda(mxprms),
     +                cfa(mxprms),cga(mxprms),
     +                bg(mxprms),csb(mxprms),cpb(mxprms),cdb(mxprms),
     +                cfb(mxprms),cgb(mxprms),
     +                cgg(mxprms),csc(mxprms),cpc(mxprms),cdc(mxprms),
     +                cfc(mxprms),cgc(mxprms),
     +                dg(mxprms),csd(mxprms),cpd(mxprms),cdd(mxprms),
     +                cfd(mxprms),cgd(mxprms),
     +                pi,qi,ri,pj,qj,rj,rri,pk,qk,rk,pl,ql,rl,rrk,
     +                nga,ngb,ngc,ngd
c
c
c flag specifying if we are doing a dft calc
c
      logical occpdft
      common/ccpdft/occpdft
c
c  CCPDFT API declarations
c
      logical CD_2e
      integer CD_4c2eon
      integer CD_abort
      integer CD_accuracy
      logical CD_active
      integer CD_assign_grid
      integer CD_ang_npoints_row
      integer CD_auto_ang_prune
      logical CD_check_print
      integer CD_chf_dksm_mo
      integer CD_chf_lhs_ao
      integer CD_chf_lhs_mo
      integer CD_chf_rhs_ao
      integer CD_chf_rhs_mo
      integer CD_clone_grid
      integer CD_conv_prune_on
      integer CD_create_grid
      integer CD_debug
      integer CD_defaults
      integer CD_defaults_old
      integer CD_dksm_exp_ao
      integer CD_dksm_exp_mo
      integer CD_euleron
      integer CD_gausslon
      integer CD_generation
      logical CD_gradcorr
      integer CD_gradquad
      integer CD_gridatomradius
      integer CD_gridscale
      integer CD_hess_ao
      integer CD_hess_mo
      logical CD_HF_coulomb
      logical CD_HF_coulomb_deriv
      logical CD_has_HF_exchange
      real*8    CD_has_HF_exchange_weight
      logical CD_HF_exchange
      real*8    CD_HF_exchange_weight
      logical CD_ignore_accuracy
      integer CD_init
      integer CD_inttol
      logical CD_is_rks
      logical CD_is_jfiton
      logical CD_is_jfitmem
      integer CD_jfit_clean1
      integer CD_jfit_clean2
      logical CD_jfit_incore
      integer CD_jfit_init1
      integer CD_jfit_init2
      integer CD_jfitoff
      integer CD_jfiton
      integer CD_jfitgon
      integer CD_jmulton
      integer CD_lebedevon
      integer CD_logon
      integer CD_lypon
      integer CD_memreq_chf_dksm_ao
      integer CD_memreq_chf_dksm_mo
      integer CD_memreq_chf_lhs_mo
      integer CD_memreq_chf_rhs_ao
      integer CD_memreq_chf_rhs_mo
      integer CD_memreq_energy
      integer CD_memreq_energy_ao
      integer CD_memreq_energy_mo
      integer CD_MHL_ang_prune
      integer CD_over
      integer CD_pener
      integer CD_pole
      integer CD_pruneatomradius
      integer CD_psitol
      integer CD_rad_npoints_row
      integer CD_radscale_scheme
      integer CD_request
      logical CD_request_multstate
      integer CD_reset_2e
      integer CD_rks
      integer CD_schwarz
      integer CD_screen
      integer CD_screenatomradius
      integer CD_set_2e
      integer CD_set_functional
      integer CD_set_ignore_accuracy
      integer CD_set_print_level
      integer CD_set_weight
      integer CD_sortpoints
      integer CD_energy
      integer CD_energy_ao
      integer CD_energy_mo
      integer CD_forces_ao
      integer CD_uks
      integer CD_xcfiton
      integer CD_import_geom
      integer CD_update_geom
      integer CD_warn
      integer CD_weightatomradius
c
      integer gden_init
      integer gden_energy
      integer gden_forces
c
c     declare API routines as external
c
      external CD_2e
      external CD_4c2eon
      external CD_abort
      external CD_accuracy
      external CD_active
      external CD_assign_grid
      external CD_auto_ang_prune
      external CD_check_print
      external CD_chf_dksm_ao
      external CD_chf_dksm_mo
      external CD_chf_lhs_ao
      external CD_chf_lhs_mo
      external CD_chf_rhs_ao
      external CD_chf_rhs_mo
      external CD_clone_grid
      external CD_conv_prune_on
      external CD_create_grid
      external CD_debug
      external CD_defaults
      external CD_defaults_old
      external CD_dksm_exp_ao
      external CD_dksm_exp_mo
      external CD_euleron
      external CD_gausslon
      external CD_generation
      external CD_gradcorr
      external CD_gradquad
      external CD_gridatomradius
      external CD_gridscale
      external CD_hess_ao
      external CD_hess_mo
      external CD_HF_coulomb
      external CD_HF_coulomb_deriv
      external CD_has_HF_exchange
      external CD_has_HF_exchange_weight
      external CD_HF_exchange
      external CD_HF_exchange_weight
      external CD_ignore_accuracy
      external CD_init
      external CD_inttol
      external CD_is_rks
      external CD_is_jfiton
      external CD_is_jfitmem
      external CD_jfit_clean1
      external CD_jfit_clean2
      external CD_jfit_incore
      external CD_jfit_init1
      external CD_jfit_init2
      external CD_jfitoff
      external CD_jfiton
      external CD_jfitgon
      external CD_jmulton
      external CD_lebedevon
      external CD_logon
      external CD_memreq_chf_dksm_mo
      external CD_memreq_chf_lhs_mo
      external CD_memreq_chf_rhs_mo
      external CD_memreq_energy
      external CD_memreq_energy_ao
      external CD_memreq_energy_mo
      external CD_MHL_ang_prune
      external CD_over
      external CD_pener
      external CD_pole
      external CD_psitol
      external CD_request
      external CD_request_multstate
      external CD_reset_2e
      external CD_rks
      external CD_schwarz
      external CD_screen
      external CD_screenatomradius
      external CD_set_2e
      external CD_set_functional
      external CD_set_ignore_accuracy
      external CD_set_print_level
      external CD_set_weight
      external CD_sortpoints
      external CD_energy
      external CD_energy_ao
      external CD_energy_mo
      external CD_forces_ao
      external CD_uks
      external CD_xcfiton
      external CD_import_geom
      external CD_update_geom
      external CD_warn
      external CD_weightatomradius
c
      external gden_init
      external gden_energy
      external gden_forces
c
c print control
c
      integer PRINT_NONE
      parameter(PRINT_NONE=0)

      integer PRINT_LOW
      parameter(PRINT_LOW=2)

      integer PRINT_DEFAULT
      parameter(PRINT_DEFAULT=5)

      integer PRINT_HIGH
      parameter(PRINT_HIGH=7)

      integer PRINT_ALL
      parameter(PRINT_ALL=10)
      real*8 ovtest,testlim,IL_shlove_tol
      logical IL_test4c
      dimension icountb(3)
c
      dimension fock(*),fockb(*),exch(*),dens(*),densb(*)
      dimension prefac(*),rdmat(*)
      dimension ib(4,4)
      dimension iso(nshels,*),q(*),gout(*)
      dimension mi(48),mj(48),mk(48),m0(48)
      data ib/64,16,4,1,216,36,6,1,1000,100,10,1,
     +         3375,225,15,1 /
c
c set control switches for penetration tests and fock build
c
      odft = CD_2e()
c
c exchange
c
      if (odft) then
         if(CD_HF_exchange())then
            facex = CD_HF_exchange_weight()
            oexch = .true.
         else
            facex = 0.0d0
            oexch = .false.
         endif
c multipoles
         odft_jmult = CD_request_multstate()
c coulomb
         ocoul = CD_HF_coulomb()
      else
         facex = 1.0d0
         odft_jmult = .false.
         oexch = .true.
         ocoul = .true.
      endif
c
c coulomb weights, modified later in multipole case
c
      fac1 = 1.0d0
      fac2 = 1.0d0
c
c  whether to compute coulomb terms for a shell quartet
c  [ reset later from penetration tests when in multipole mode ]
c
      odft_test = .true.
c
c initialise counters for screening statistics
c
      nmult = 0
      ncomp = 0
      icountb(1)=0
      icountb(2)=0
      icountb(3)=0
c
      testlim=IL_shlove_tol()
      if(opg_root() .and. odft_jmult)
     &     write(6,*) 'Using ',testlim,' for overlap tol'
c
czora    zora uses dft flags to build coulomb matrix
c
      if (oint_zora.or.osmall_zora) then
          ocoul = .true.
          if (isexch_z .eq. 1) then
             if (.not. oatscf_z)
     1      write(iwr,*)'exchange included in ZORA correction (intega)'
             oexch = .true.
          else 
             oexch = .false.
          end if
      end if
c
c     ----- two-electron integrals -----
c
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
c     ----- allocate core memory
c
      if(.not.odscf) then
         call rdmake(prefac)
         dlntol=-dlog(cutoff*0.1d0)
       else
         l2 = ikyp(numorb)
       endif
c
c  initialisation  rf contribution two-electron integrals
      if (field(5:) .eq. 'scf' .and. intdrf .eq. 1) call rftwin(q)
c
      time = cpulft(1)
      tim0 = time
      tim1 = time
c
      ist0 = ist
      jst0 = jst
      kst0 = kst
      lst0 = lst
c
c     ----- ishell -----
c
      if (ist0.le.nshels) then
         do 150 ii = ist0 , nshels
c
c     ----- print intermediate restart data -----
c
            kadi = kad(ii)
            dt0 = time - tim0
            dt1 = time - tim1
            tim1 = time
            if (outv) then
               if (odscf) then
                  if (iter.le.0) write (iwr,6020) ii , jst0 , kst0 , 
     +                           lst0 , nrec , dt1 , dt0
                  if (ologf) then
                     write (text,6020) ii , jst0 , kst0 , lst0 , nrec ,
     +                                 dt1 , dt0
                     call sptchk(text)
                  end if
               else
                  write (iwr,6010) ii , jst0 , kst0 , lst0 , nrec , 
     +                             icount , dt1 , dt0
               end if
            end if
c
c     ----- eliminate ishell -----
c
            do 50 it = 1 , nt
               id = iso(ii,it)
               if (id.gt.ii) go to 150
               m0(it) = id
 50         continue
            ikyii = iky(ii)
c
c     ----- jshell -----
c
            j0 = jst0
            do 140 jj = j0 , ii
               jst0 = 1
               itrij = ikyii + jj
               tolij = dlntol + prefac(itrij)
               do 70 it = 1 , nt
                  jd = iso(jj,it)
                  if (jd.gt.ii) go to 140
                  id = m0(it)
                  if (id.lt.jd) then
                     nd = id
                     id = jd
                     jd = nd
                  end if
                  if (id.eq.ii .and. jd.gt.jj) go to 140
                  mi(it) = id
                  mj(it) = jd
 70            continue
               kadij = kad(jj) + kadi
               mij = itrij
c
c     ----- get information about i-shell and j-shell -----
c
               ish = ii
               jsh = jj
               call shells(gout,1,ish,jsh,ksh,lsh,1)
               if (odft_jmult) then
                  if (ovtest(ii,jj,rri).lt.testlim) goto 140
               end if
               call ijprim
               if (nij.ne.0) then
c
c     ----- kshell -----
c
c
c store information for penetration tests relating to i and j
c shells
                  if(odft_jmult)call IL_test4a(ii,jj)
                  k0 = kst0
                  do 130 kk = k0 , ii
                     kst0 = 1
                     do 90 it = 1 , nt
                        kd = iso(kk,it)
                        if (kd.gt.ii) go to 130
                        mk(it) = kd
 90                  continue
                     kadijk = kad(kk) + kadij
                     ikykk = iky(kk)
                     if (odscf) then
                        itrik = ikyii + kk
                        mik = itrik
                        itrjk = iky(max(jj,kk)) + min(jj,kk)
                        mjk = itrjk
                        tijk = dmax1(rdmat(mij),rdmat(mik),rdmat(mjk))
                     end if
                     if(odft_jmult)call IL_test4b(kk)
c
c     ----- lshell ----
c
                     maxll = kk
                     if (kk.eq.ii) maxll = jj
                     l0 = lst0
                     do 120 ll = l0 , maxll
czora
                      if (oatint_z) then
                       if (
     1                    (katom(jj).ne.katom(ii).or.
     2                     katom(kk).ne.katom(ii).or.
     3                     katom(ll).ne.katom(ii))) go to 120
                      end if
czora
                      lst0 = 1
                      if (kadijk+kad(ll).lt.0) then
                        itrkl = ikykk+ll
                        tijkl = tolij+prefac(itrkl)
                        if (tijkl.le.0.0d0) then
                          if (odscf) intcut(2) = intcut(2)+1
                          goto 120
                        end if
                        if (odscf) then
                          mil = ikyii+ll
                          mjl = iky(max(jj,ll))+min(jj,ll)
                          mkl = itrkl
                          tijkl = tijkl+dmax1(tijk,rdmat(mil),
     +                                        rdmat(mjl),rdmat(mkl))
                          if (tijkl.le.0.0d0) then
                            intcut(3) = intcut(3)+1
                            goto 120
                          end if
                        end if
                        n4 = 0
                        do 110 it = 1 , nt
                           ld = iso(ll,it)
                           if (ld.gt.ii) go to 120
                           kd = mk(it)
                           if (kd.lt.ld) then
                              nd = kd
                              kd = ld
                              ld = nd
                           end if
                           id = mi(it)
                           jd = mj(it)
                           if (id.eq.ii .or. kd.eq.ii) then
                              if (kd.ge.id) then
                                 if (kd.ne.id .or. ld.gt.jd) then
                                    nd = id
                                    id = kd
                                    kd = nd
                                    nd = jd
                                    jd = ld
                                    ld = nd
                                 end if
                              end if
                              if (jd.ge.jj) then
                                 if (jd.gt.jj) go to 120
                                 if (kd.ge.kk) then
                                    if (kd.gt.kk) go to 120
                                    if (ld.ge.ll) then
                                       if (ld.gt.ll) go to 120
                                       n4 = n4 + 1
                                    end if
                                 end if
                              end if
                           end if
 110                    continue
                        q4 = dble(nt)/dble(n4)
c
c     ----- (ii,jj//kk,ll) -----
c
                        ish = ii
                        jsh = jj
                        ksh = kk
                        lsh = ll
                        qq4 = q4
c
c     ----- initialize gout to zero -----
c     ----- get information about ksh and lsh -----
c
                        call shells(gout,2,ish,jsh,ksh,lsh,1)
                        if (odft_jmult) then
                           if (ovtest(kk,ll,rrk).lt.testlim) goto 120
                        end if
c
c compute weights and test for relevance of this block of integrals
c also returns stats in icountb
c
                        if (odft_jmult)
     &                     odft_test=IL_test4c(ll,fac1,fac2,icountb)
c
                        if (oexch .or. odft_test) then
c
c     ----- compute two-electron integrals ----
c     ----- write them out on mainfile -----
c
                        call genral(gout)
                        if (odscf) then
                           if (zscftp.eq.'uhf')then
                              call dir_build_uhf(fock,fockb,
     +                             dens,densb,gout,
     +                             fac1,fac2, facex, ocoul, oexch)
                           else if (zscftp.eq.'gvb'.or.
     +                              zscftp.eq.'grhf') then
                              if (nsheld.le.1) then
                                 call dir_build_open(fock,exch,
     +                                dens,gout)
                              else
                                 call dir_build_open2(l2,fock,
     +                                exch,dens,gout)
                              endif
                           else
                              if (omorok) then
                              call dbuild_morok(fock,dens,gout)
                              else
                              call dbuild(fock,dens,gout,
     +                                    fac1,fac2,facex,ocoul,oexch)
                              endif
                        endif
                        else if (field(5:) .eq. 'scf' .and.
     +                           intdrf .eq. 1) then
                           call qoutdrf(q,gout)
                        else
                           call qout(gout)
                        endif
                        ncomp=ncomp+1
                     else
                        nmult=nmult+1
                     endif
c
c
c     ----- check cpu time/ maxblock condition -----
c
                     call chkout(ii,jj,kk,ll,fock,q)
                     if (omaxb .or. tim.gt.timlim) go to 160
                  end if
c
 120           continue
 130        continue
            end if
 140        continue
            time = cpulft(1)
 150     continue
         call final(q,fock,dens)
 160     continue
      end if
c
c     ----- reset core memory
c
      if(odft_jmult)then
      write(iwr,*)'shell quartets: ',ncomp,'computed ',nmult,' skipped'
      write(iwr,*)'counts',icountb
      endif
c
      return
 6010 format (i4,3i5,1x,i10,i9,f11.2,f9.2)
 6020 format (i4,3i5,1x,i10,9x,f11.2,f9.2)
      end
      subroutine jkints(zscftp,q,fock,fockb,exch,dens,densb,prefac,
     +                  rdmat,iso,gout,nshels,outvv)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *59 text
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
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
      parameter (mxp2 = mxprms * mxprms)
c
c     dlnmxd: the maximum value of the logarithm of the reduced
c             density matrix.
c     dlnmxs: the maximum value of the logarithm of the Schwarz 
c             integrals.
c     dlntol: minus the logarithm of the direct SCF integral cutoff.
c
      real*8 dlnmxd, dlnmxs, dlntol, tolitr, deltol, delfac
      integer intcut, intmag
      integer ibl171, lentri, itrtol, iexch, iof171, m171t
      integer nhamd, nsheld
      logical odscf, oimag, odnew, opref, odelta
      common/cslosc/dlnmxd,dlnmxs,dlntol,tolitr(3),deltol,delfac,
     + odscf,
     + intcut(10),oimag,intmag(1060),
     + ibl171,lentri,itrtol(3),iexch,odnew,opref,iof171(6),m171t,
     + nhamd, nsheld,odelta
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      common/scfopt/maxit(4),accdi(2),odiis(4),dmpcut(7),iter
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
      logical odumpm, ooprnt, opun7, omax, ologf
      integer maxbrd, logf, iofsp
      common /utilc/ odumpm,ooprnt,opun7,omax,maxbrd,ologf,
     +               logf,iofsp(4)
c
c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
c
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
c
      logical ospbas, onocnt, ochek, opdbas, opfbas, opgbas
      integer kad, istd
      common /ijlab/ kad(mxshel),ospbas,onocnt,istd,ochek,
     +               opdbas,opfbas,opgbas
c
c
c  control parameters for morokuma energy decomposition
c  analysis
c
      logical omorok, oifc, oife, oifh
      integer imode, ifrag, ifocc, imoro
      real*8 emorok
      common/morok1/emorok(0:4,2),imode,ifrag(maxorb),
     &     ifocc(maxorb),imoro, omorok, 
     &     oifc, oife, oifh
      character*20 fragname
      common/morokc/fragname(3)

c
c  settings
c
      integer ALL_BLOCKS
      parameter (ALL_BLOCKS=0)
      
      integer ESX
      parameter (ESX=1)

      integer ESX_PLX
      parameter (ESX_PLX=2)

      integer ESX_CTT
      parameter (ESX_CTT=3)

      integer ESX_EX
      parameter (ESX_EX=4)

      integer YES_K
      parameter (YES_K=1)

      integer NO_K
      parameter (NO_K=2)

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
      common/junk/cxyz(3,5625),aaa(21*mxp2),ijaaa(225)
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
c     logical opk,opank,oskips
c     common/pkfil/opk,opank,oskips
      logical opk,opank
      common/pkfil/opk,opank
c
      integer igt, jgt, kgt, lgt
      integer ibt, jbt, kbt, lbt, lenddd
      common /flips/ igt(3),jgt(3),kgt(3),lgt(3),
     +               ibt(mxprms),jbt(mxprms),kbt(mxprms),lbt(mxprms),
     +               lenddd
c
      logical osortp,oschw
      common/sortpdat/
     + igbf,iklbf,iiptbf,ijbase,ngbf,lengbf,osortp,oschw
      integer icount_dlb, icount_slb
      common /par_dlb/ icount_dlb, icount_slb
c
c
c flag specifying if we are doing a dft calc
c
      logical occpdft
      common/ccpdft/occpdft
c
c  CCPDFT API declarations
c
      logical CD_2e
      integer CD_4c2eon
      integer CD_abort
      integer CD_accuracy
      logical CD_active
      integer CD_assign_grid
      integer CD_ang_npoints_row
      integer CD_auto_ang_prune
      logical CD_check_print
      integer CD_chf_dksm_mo
      integer CD_chf_lhs_ao
      integer CD_chf_lhs_mo
      integer CD_chf_rhs_ao
      integer CD_chf_rhs_mo
      integer CD_clone_grid
      integer CD_conv_prune_on
      integer CD_create_grid
      integer CD_debug
      integer CD_defaults
      integer CD_defaults_old
      integer CD_dksm_exp_ao
      integer CD_dksm_exp_mo
      integer CD_euleron
      integer CD_gausslon
      integer CD_generation
      logical CD_gradcorr
      integer CD_gradquad
      integer CD_gridatomradius
      integer CD_gridscale
      integer CD_hess_ao
      integer CD_hess_mo
      logical CD_HF_coulomb
      logical CD_HF_coulomb_deriv
      logical CD_has_HF_exchange
      real*8    CD_has_HF_exchange_weight
      logical CD_HF_exchange
      real*8    CD_HF_exchange_weight
      logical CD_ignore_accuracy
      integer CD_init
      integer CD_inttol
      logical CD_is_rks
      logical CD_is_jfiton
      logical CD_is_jfitmem
      integer CD_jfit_clean1
      integer CD_jfit_clean2
      logical CD_jfit_incore
      integer CD_jfit_init1
      integer CD_jfit_init2
      integer CD_jfitoff
      integer CD_jfiton
      integer CD_jfitgon
      integer CD_jmulton
      integer CD_lebedevon
      integer CD_logon
      integer CD_lypon
      integer CD_memreq_chf_dksm_ao
      integer CD_memreq_chf_dksm_mo
      integer CD_memreq_chf_lhs_mo
      integer CD_memreq_chf_rhs_ao
      integer CD_memreq_chf_rhs_mo
      integer CD_memreq_energy
      integer CD_memreq_energy_ao
      integer CD_memreq_energy_mo
      integer CD_MHL_ang_prune
      integer CD_over
      integer CD_pener
      integer CD_pole
      integer CD_pruneatomradius
      integer CD_psitol
      integer CD_rad_npoints_row
      integer CD_radscale_scheme
      integer CD_request
      logical CD_request_multstate
      integer CD_reset_2e
      integer CD_rks
      integer CD_schwarz
      integer CD_screen
      integer CD_screenatomradius
      integer CD_set_2e
      integer CD_set_functional
      integer CD_set_ignore_accuracy
      integer CD_set_print_level
      integer CD_set_weight
      integer CD_sortpoints
      integer CD_energy
      integer CD_energy_ao
      integer CD_energy_mo
      integer CD_forces_ao
      integer CD_uks
      integer CD_xcfiton
      integer CD_import_geom
      integer CD_update_geom
      integer CD_warn
      integer CD_weightatomradius
c
      integer gden_init
      integer gden_energy
      integer gden_forces
c
c     declare API routines as external
c
      external CD_2e
      external CD_4c2eon
      external CD_abort
      external CD_accuracy
      external CD_active
      external CD_assign_grid
      external CD_auto_ang_prune
      external CD_check_print
      external CD_chf_dksm_ao
      external CD_chf_dksm_mo
      external CD_chf_lhs_ao
      external CD_chf_lhs_mo
      external CD_chf_rhs_ao
      external CD_chf_rhs_mo
      external CD_clone_grid
      external CD_conv_prune_on
      external CD_create_grid
      external CD_debug
      external CD_defaults
      external CD_defaults_old
      external CD_dksm_exp_ao
      external CD_dksm_exp_mo
      external CD_euleron
      external CD_gausslon
      external CD_generation
      external CD_gradcorr
      external CD_gradquad
      external CD_gridatomradius
      external CD_gridscale
      external CD_hess_ao
      external CD_hess_mo
      external CD_HF_coulomb
      external CD_HF_coulomb_deriv
      external CD_has_HF_exchange
      external CD_has_HF_exchange_weight
      external CD_HF_exchange
      external CD_HF_exchange_weight
      external CD_ignore_accuracy
      external CD_init
      external CD_inttol
      external CD_is_rks
      external CD_is_jfiton
      external CD_is_jfitmem
      external CD_jfit_clean1
      external CD_jfit_clean2
      external CD_jfit_incore
      external CD_jfit_init1
      external CD_jfit_init2
      external CD_jfitoff
      external CD_jfiton
      external CD_jfitgon
      external CD_jmulton
      external CD_lebedevon
      external CD_logon
      external CD_memreq_chf_dksm_mo
      external CD_memreq_chf_lhs_mo
      external CD_memreq_chf_rhs_mo
      external CD_memreq_energy
      external CD_memreq_energy_ao
      external CD_memreq_energy_mo
      external CD_MHL_ang_prune
      external CD_over
      external CD_pener
      external CD_pole
      external CD_psitol
      external CD_request
      external CD_request_multstate
      external CD_reset_2e
      external CD_rks
      external CD_schwarz
      external CD_screen
      external CD_screenatomradius
      external CD_set_2e
      external CD_set_functional
      external CD_set_ignore_accuracy
      external CD_set_print_level
      external CD_set_weight
      external CD_sortpoints
      external CD_energy
      external CD_energy_ao
      external CD_energy_mo
      external CD_forces_ao
      external CD_uks
      external CD_xcfiton
      external CD_import_geom
      external CD_update_geom
      external CD_warn
      external CD_weightatomradius
c
      external gden_init
      external gden_energy
      external gden_forces
c
c print control
c
      integer PRINT_NONE
      parameter(PRINT_NONE=0)

      integer PRINT_LOW
      parameter(PRINT_LOW=2)

      integer PRINT_DEFAULT
      parameter(PRINT_DEFAULT=5)

      integer PRINT_HIGH
      parameter(PRINT_HIGH=7)

      integer PRINT_ALL
      parameter(PRINT_ALL=10)
c
      dimension fock(*),fockb(*),exch(*),dens(*),densb(*)
      dimension prefac(*),rdmat(*)
      dimension ib(4,4)
      dimension iso(nshels,*),q(*),gout(*)
      dimension mi(48),mj(48),mk(48),m0(48)
      character *8 fnm
      character *6 snm
      data fnm/'intega.m'/
      data snm/'jkints'/
      data ib/64,16,4,1,216,36,6,1,1000,100,10,1,
     +         3375,225,15,1 /
      data m25/25/
c
c for simplicity, only implement full coulomb version
c here for the moment
c
      if (CD_2e())then
        odft = .true.

        if(CD_request_multstate())
     &    call caserr('disable schwarz inequality for multipole runs')

        if(CD_HF_exchange())then
           facex = CD_HF_exchange_weight()
           oexch = .true.
        else
           facex = 0.0d0
           oexch = .false.
        endif
      else
         odft = .false.
         facex = 1.0d0
         oexch = .true.
      endif
c
c coulomb
      ocoul = CD_HF_coulomb() .or. .not. odft
c
c coulomb weights, modified later in multipole case
c
      fac1 = 1.0d0
      fac2 = 1.0d0
c
c     ----- two-electron integrals -----
c
      ltri = ikyp(nshels)
      l2 = ikyp(numorb)
      oskipp = .false.
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
      dlncutoff = dlog(cutoff)
c
c     ----- allocate core memory
c
      ischw = igmem_alloc_inf(ltri,fnm,snm,'schwarz',IGMEM_DEBUG)
      nschwz = 0
c
c *** read in ints for schwarz inequality test
c
      call secget(isect(421),m25,iblk25)
      call rdedx(q(ischw),ltri,iblk25,idaf)

      dlnmxs = -1.0d50
      do ii = 0, ltri-1
         dlnmxs = dmax1(dlnmxs,q(ischw+ii))
      enddo
c
      time = cpulft(1)
      tim0 = time
      tim1 = time
c
      ist0 = ist
      jst0 = jst
      kst0 = kst
      lst0 = lst
c
c     ----- ishell -----
c
      if (ist0.le.nshels) then
         do 150 ii = ist0 , nshels
c
c     ----- print intermediate restart data -----
c
            kadi = kad(ii)
            dt0 = time - tim0
            dt1 = time - tim1
            tim1 = time
            if (outv) then
               if (odscf) then
                  if (iter.le.0) write (iwr,6020) ii , jst0 , kst0 , 
     +                                  lst0 , nrec , dt1 , dt0
                  if (ologf) then
                     write (text,6020) ii , jst0 , kst0 , lst0 , nrec ,
     +                                 dt1 , dt0
                     call sptchk(text)
                  end if
               else
                  write (iwr,6010) ii , jst0 , kst0 , lst0 , nrec , 
     +                             icount , dt1 , dt0
               end if
            end if
c
c     ----- eliminate ishell -----
c
            do 50 it = 1 , nt
               id = iso(ii,it)
               if (id.gt.ii) go to 150
               m0(it) = id
 50         continue
            ikyii = iky(ii)
c
c     ----- jshell -----
c
            j0 = jst0
            do 140 jj = j0 , ii
            jst0 = 1
            itrij = ikyii + jj
             do 70 it = 1 , nt
                 jd = iso(jj,it)
                 if (jd.gt.ii) go to 140
                 id = m0(it)
                 if (id.lt.jd) then
                    nd = id
                    id = jd
                    jd = nd
                 end if
                 if (id.eq.ii .and. jd.gt.jj) go to 140
                 mi(it) = id
                 mj(it) = jd
 70           continue
               ijij = itrij + ischw - 1
               if (odscf) then
                 tijkl = dlntol + q(ijij) + dlnmxs + dlnmxd
                 if (tijkl .le. 0.0d0) then
c
c                   if the (ij/ij) integral is too small we can 
c                   skip the whole kl loop.
c
                    nschwz = nschwz + itrij
                    intcut(1)=intcut(1)+1
                    go to 140
                 endif
               else
                 if (q(ijij)+dlnmxs.lt.dlncutoff) then
c                   see above...
                    nschwz = nschwz + itrij
                    go to 140
                 endif
               endif
               kadij = kad(jj) + kadi
               mij = itrij
c
c     ----- get information about i-shell and j-shell -----
c
               ish = ii
               jsh = jj
               call shells(gout,1,ish,jsh,ksh,lsh,1)
               call ijprim
               if (nij.ne.0) then
c
c     ----- kshell -----
c
                  k0 = kst0
                  do 130 kk = k0 , ii
                     kst0 = 1
                     do 90 it = 1 , nt
                        kd = iso(kk,it)
                        if (kd.gt.ii) go to 130
                        mk(it) = kd
 90                  continue
                     kadijk = kad(kk) + kadij
                     ikykk = iky(kk)
                     if (odscf) then
                        itrik = ikyii + kk
                        mik = itrik
                        itrjk = iky(max(jj,kk)) + min(jj,kk)
                        mjk = itrjk
                        tijk = dmax1(rdmat(mij),rdmat(mik),rdmat(mjk))
                     end if
c
c     ----- lshell ----
c
                     maxll = kk
                     if (kk.eq.ii) maxll = jj
                     l0 = lst0
                     do 120 ll = l0 , maxll
                        lst0 = 1
                        if (kadijk+kad(ll).lt.0) then
                           itrkl = ikykk + ll
c                          ijij = itrij + ischw - 1
                           klkl = itrkl + ischw - 1
                           test = q(ijij)+q(klkl)
                           if (odscf) then
                              mil = ikyii+ll
                              mjl = iky(max(jj,ll)) 
     +                            + min(jj,ll)
                              mkl = itrkl
                              tijkl = dlntol + test
     +                              + dmax1(tijk,rdmat(mil),
     +                                      rdmat(mjl),rdmat(mkl))
                              if (tijkl.le.0.0d0) then
                                 nschwz = nschwz + 1
                                 intcut(3)=intcut(3)+1
                                 goto 120
                              end if
                           else
                              oskipp = test.lt.dlncutoff
                              if (oskipp) then
                                 nschwz = nschwz + 1
                                 go to 120
                              end if
                           end if
                           n4 = 0
                           do 110 it = 1 , nt
                              ld = iso(ll,it)
                              if (ld.gt.ii) go to 120
                              kd = mk(it)
                              if (kd.lt.ld) then
                                 nd = kd
                                 kd = ld
                                 ld = nd
                              end if
                              id = mi(it)
                              jd = mj(it)
                              if (id.eq.ii .or. kd.eq.ii) then
                                 if (kd.ge.id) then
                                    if (kd.ne.id .or. ld.gt.jd) then
                                       nd = id
                                       id = kd
                                       kd = nd
                                       nd = jd
                                       jd = ld
                                       ld = nd
                                    end if
                                 end if
                                 if (jd.ge.jj) then
                                    if (jd.gt.jj) go to 120
                                    if (kd.ge.kk) then
                                       if (kd.gt.kk) go to 120
                                       if (ld.ge.ll) then
                                          if (ld.gt.ll) go to 120
                                          n4 = n4 + 1
                                       end if
                                    end if
                                 end if
                              end if
 110                       continue
                           q4 = dble(nt)/dble(n4)
c
c     ----- (ii,jj//kk,ll) -----
c
                           ish = ii
                           jsh = jj
                           ksh = kk
                           lsh = ll
                           qq4 = q4
c
c     ----- initialize gout to zero -----
c     ----- get information about ksh and lsh -----
c
                           call shells(gout,2,ish,jsh,ksh,lsh,1)
c
c     ----- compute two-electron integrals ----
c     ----- write them out on mainfile -----
c
                           call genral(gout)
                           if (odscf) then
                              if (zscftp.eq.'uhf') then
                              call dir_build_uhf(fock,fockb,
     +                             dens,densb,gout,
     +                             fac1,fac2, facex, ocoul, oexch)
                              else if (zscftp.eq.'gvb'.or.
     +                                 zscftp.eq.'grhf') then
                               if(nsheld.le.1) then
                                call dir_build_open(fock,exch,
     +                               dens,gout)
                               else
                                call dir_build_open2(l2,fock,
     +                               exch,dens,gout)
                               endif
                              else
                              if (omorok) then
                              call dbuild_morok(fock,dens,gout)
                              else
                              call dbuild(fock,dens,gout,
     +                                    fac1,fac2,facex,ocoul,oexch)
                              endif
                              endif
                           else
                              call qout(gout)
                           endif
c
c
c     ----- check cpu time/ maxblock condition -----
c
                           call chkout(ii,jj,kk,ll,fock,q)
                           if (omaxb .or. tim.gt.timlim) go to 160
                        end if
c
 120                 continue
 130              continue
               end if
 140        continue
            time = cpulft(1)
 150     continue
         call final(q,fock,dens)
 160     if(outvv) write(iwr,6030) nschwz
      end if
c
c     ----- reset core memory
c
      call gmem_free_inf(ischw,fnm,snm,'schwarz')

      return
 6010 format (i4,3i5,1x,i10,i9,f11.2,f9.2)
 6020 format (i4,3i5,1x,i10,9x,f11.2,f9.2)
 6030 format(/1x,'schwarz inequality test skipped ',i10,
     + ' integral blocks')
      end
      subroutine pkinta(q,iso,gout,nshels)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
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
      parameter (mxp2 = mxprms * mxprms)
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
c
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
c
      logical ospbas, onocnt, ochek, opdbas, opfbas, opgbas
      integer kad, istd
      common /ijlab/ kad(mxshel),ospbas,onocnt,istd,ochek,
     +               opdbas,opfbas,opgbas
c
      common/junk/cxyz(3,5625),aaa(21*mxp2),ijaaa(225)
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
c     logical opk,opank,oskips
c     common/pkfil/opk,opank,oskips
      logical opk,opank
      common/pkfil/opk,opank
c
      integer igt, jgt, kgt, lgt
      integer ibt, jbt, kbt, lbt, lenddd
      common /flips/ igt(3),jgt(3),kgt(3),lgt(3),
     +               ibt(mxprms),jbt(mxprms),kbt(mxprms),lbt(mxprms),
     +               lenddd
c
      logical osortp,oschw
      common/sortpdat/
     + igbf,iklbf,iiptbf,ijbase,ngbf,lengbf,osortp,oschw
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      integer icount_dlb, icount_slb
      common /par_dlb/ icount_dlb, icount_slb
c
      dimension iso(nshels,*),q(*),gout(*)
c
      character*8 field
      character*24 version
      common /drfopt/ field,version
      logical oreact
      integer intdrf
      common /drfopti/ oreact,intdrf
c
      logical odrf
c
      integer ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep
      integer neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga
      integer idipcal,iclinte,iclintd,ieffpol,isodis,iclintr
      integer iextdip,itsolv,itermax,isolsav,imomsav,idisadd
      integer ianal,ineqex,ithole,irevdis,ihbond,igrppol
      integer maxci_drf
      common/drfpar1/ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep,
     +               neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga,
     +               idipcal,iclinte,iclintd,ieffpol,isodis,iclintr,
     +               iextdip,itsolv,itermax,isolsav,imomsav,idisadd,
     +               ianal,ineqex,ithole,irevdis,ihbond,igrppol,
     +               maxci_drf
c
      integer nxtpts, npol, nexp, natint, nshint, namb, nspec, ngran
      common/drfpar2/nxtpts,npol,nexp,natint,nshint,namb,nspec,ngran
c
      integer ndima, ndimb, nbem, ndim, nexp4, npol3, nwtr, nwtc, nzfa
      integer lomga, nomga, nchd, ngrnam, ngrnam2, nzfp, nzfn
      integer nneq, nneqrf, maxneq, neqdim, nqdim, nqcls
      common/drfpar3/ndima,ndimb,nbem,ndim,nexp4,npol3,nwtr,nwtc,nzfa,
     +               lomga,nomga,nchd,ngrnam,ngrnam2,nzfp,nzfn,
     +               nneq,nneqrf,maxneq,neqdim,nqdim,nqcls
c
      integer ifldin, ifldout, idrfout, modxza, irepopt
      integer iexpza, icmexp, iadexp, nodpe, igetden, iarfcal, iqmclr
      common/drfpar4/ifldin,ifldout,idrfout,modxza,irepopt,
     +               iexpza,icmexp,iadexp,nodpe,igetden,iarfcal,
     +               iqmclr
c
      real*8 gamdrf, dstmin, dstmax, hbondl, hbondr, afact, rfact
      real*8 cvgrel, agrpe, agrpm, agrpc, scffact, acur
      real*8 conci_drf
      common/drfpar5/gamdrf,dstmin,dstmax,hbondl,hbondr,afact,rfact,
     +               cvgrel,agrpe,agrpm,agrpc,scffact,acur,conci_drf
c
      integer iind, isur1, isur2, ilwt, ilvr, illur, ilindx
      common/drfpar6/iind,isur1,isur2,ilwt,ilvr,illur,ilindx
c
      integer i10, i20, i30, i40, i50, i60, i70, i80, i90 
      integer i100, i110, i120, i130, i140, i150, i160
      integer i170, i180
      common /scmprt/ i10,i20,i30,i40,i50,i60,i70,i80,i90,
     +          i100,i110,i120,i130,i140,i150,i160,i170,i180
c
c
      integer itwoeps,ioponly,itwosur,ixnsurf,ibemout,iuniout
      common /drfbem1/ itwoeps,ioponly,itwosur,ixnsurf,ibemout,iuniout
c
      integer leveli, levelo, nbem1, nbem2, nraw
      common /drfbem2/ leveli,levelo,nbem1,nbem2,nraw
c
      character*32 solnam
      real*8 radmax, eps1, eps2, spherad, rprobe, rprobej
      real*8 solrad, swidth, sdist
      real*8 kappa
      real*8 kappa1, kappa2, kappas
      common /drfbem3/ radmax,eps1,eps2,kappa1,kappa2,spherad,rprobe,
     +                 rprobej,solnam,solrad,swidth,sdist
c
      real*8 radat, cgrav
      common /drfbem4/ radat(maxat),cgrav(3)
c
      dimension mi(48),mj(48),mk(48),m0(48)
      dimension norgin(3)
      dimension ib(4,4)
c
      data norgin/0,50625,101250/
      data ib/64,16,4,1,216,36,6,1,1000,100,10,1,
     +         3375,225,15,1 /
c
      inull = igmem_null()
      nav = lenwrd()
c     ltri = ikyp(nshels)
cafc
      odrf = (field(5:) .eq. 'scf ') .and. (intdrf .ne. 0)
cxxx
c       allocate memory for drf matrices
c
      if (odrf) then
c-----  set acur for use in 2-electron reaction-field
c       routines
c
****    acurcy not passed into this routine
****    use value hardwired in scfinp (1.0d-05)
****    acur = acurcy
        acur = 1.0d-05
c  ovlap at i110, da10, 12, l2
c  dipx at  i120, da10, 53, l2
c  dipy at  i130, da10, 54, l2
c  dipz at  i140, da10, 55, l2
c  iexp at      i150, da31, 2, l2
c  omega(s) at  i160, da31, 50, nomga
c  omega(op) at i170, da31, 51, nomga (itwoeps = 1)
c
        i110 = igmem_alloc(nx)
        i120 = igmem_alloc(nx)
        i130 = igmem_alloc(nx)
        i140 = igmem_alloc(nx)
        i150 = igmem_alloc(2*nx)
        i160 = igmem_alloc(nomga)
        if (itwoeps .eq. 1) then
          i170 = igmem_alloc(nomga)
        endif
cafc  initialisation  rf contribution two-electron integrals
        call rftwin(q)
cafc
      endif
cafc
      if (opank) osortp = .false.
      ngout = 256
      lenddd = 4
      ii = 1
      if (opdbas) then
         ngout = 1296
         ii = 2
         lenddd = 6
      end if
      if (opfbas) then
         ngout = 10000
         lenddd = 10
         ii = 3
      end if
      if (opgbas) then
         ngout = 50625
         lenddd = 15
         ii = 4
      end if
      do 30 loop = 1 , 3
        igt(loop) = ib(1,ii)
        jgt(loop) = ib(2,ii)
        kgt(loop) = ib(3,ii)
        lgt(loop) = ib(4,ii)
 30   continue
      ii = 0
      do 40 loop = 1 , lenddd
         lbt(loop) = ii
         kbt(loop) = lbt(loop)*lenddd
         jbt(loop) = kbt(loop)*lenddd
         ibt(loop) = jbt(loop)*lenddd
         ii = ii + 1
 40   continue
      if (osortp) then
c *** num*maxshellsize buffers length blocksize
         ngbf = num*6
         if (opfbas) ngbf = num*10
         if (opgbas) ngbf = num*15
c ***
c *** need mods to clrgbf to use anything other than 340 here
c     allocate memory for sorting
         lengbf = num2ep
         igbf = igmem_alloc(ngbf*lengbf)
         itmp = (ngbf*lengbf+1)/nav
         iklbf = igmem_alloc(itmp)
         itmp = (ngbf+1)/nav
         iiptbf = igmem_alloc(itmp)
c *** zero out p sort buffer counts
         call setsto(ngbf,0,q(iiptbf))
      end if
cafc  initialisation  rf contribution two-electron integrals
      if (odrf) call rftwin(q)
cafc
      time = cpulft(1)
      tim0 = time
      tim1 = time
      ist0 = ist
      jst0 = jst
      kst0 = kst
      lst0 = lst
      if (ist0.le.nshels) then
         do 300 ii = ist0 , nshels
c *** set ijbase to triangle of first basis function in shell
c *** and check integrity of counters
            if (osortp) then
               ijbase = (kloc(ii)*(kloc(ii)-1))/2
               if (icount.ne.1) call caserr('icount.ne.1 for new ii')
            end if
            kadi = kad(ii)
            dt0 = time - tim0
            dt1 = time - tim1
            tim1 = time
            if (outv) write (iwr,6010) ii , jst0 , kst0 , lst0 , nrec ,
     +                                icount , dt1 , dt0
            do 60 it = 1 , nt
               id = iso(ii,it)
               if (id.gt.ii) go to 280
               mi(it) = id
 60         continue
            j0 = jst0
            do 270 jj = j0 , ii
               jst0 = 1
               kadij = kadi + kad(jj)
               do 80 it = 1 , nt
                  id = mi(it)
                  jd = iso(jj,it)
                  mj(it) = jd
                  if (id.lt.jd) then
                     nd = id
                     id = jd
                     jd = nd
                  end if
                  if (id.lt.ii) then
                  else if (id.eq.ii) then
                     if (jd.gt.jj) go to 270
                  else
                     go to 270
                  end if
 80            continue
               k0 = kst0
               do 260 kk = k0 , jj
                  kst0 = 1
                  kadijk = kadij + kad(kk)
                  do 110 it = 1 , nt
                     id = mi(it)
                     jd = mj(it)
                     kd = iso(kk,it)
                     mk(it) = kd
 100                 if (id.lt.jd) then
                        nd = id
                        id = jd
                        jd = nd
                     end if
                     if (jd.lt.kd) then
                        nd = jd
                        jd = kd
                        kd = nd
                        go to 100
                     else if (id.lt.ii) then
                     else if (id.eq.ii) then
                        if (jd.lt.jj) then
                        else if (jd.eq.jj) then
                           if (kd.gt.kk) go to 260
                        else
                           go to 260
                        end if
                     else
                        go to 260
                     end if
 110              continue
                  l0 = lst0
                  do 250 ll = l0 , kk
                     lst0 = 1
                     if (kadijk+kad(ll).ge.0) go to 250
                     n4 = 0
                     do 150 it = 1 , nt
                        id = mi(it)
                        jd = mj(it)
                        kd = mk(it)
                        ld = iso(ll,it)
 130                    if (id.lt.jd) then
                           nd = id
                           id = jd
                           jd = nd
                        end if
 140                    if (jd.lt.kd) then
                           nd = jd
                           jd = kd
                           kd = nd
                           go to 130
                        else if (kd.lt.ld) then
                           nd = kd
                           kd = ld
                           ld = nd
                           go to 140
                        else if (id.lt.ii) then
                        else if (id.eq.ii) then
                           if (jd.lt.jj) then
                           else if (jd.eq.jj) then
                              if (kd.lt.kk) then
                              else if (kd.eq.kk) then
                                 if (ld.lt.ll) then
                                 else if (ld.eq.ll) then
                                    n4 = n4 + 1
                                    m0(n4) = it
                                 else
                                    go to 250
                                 end if
                              else
                                 go to 250
                              end if
                           else
                              go to 250
                           end if
                        else
                           go to 250
                        end if
 150                 continue
                     oskpa = jj.eq.kk
                     oskpb = (ii.eq.kk) .or. (jj.eq.ll)
                     oskpc = (ii.eq.jj) .or. (kk.eq.ll)
                     onpsym = .false.
                     if (.not.(oskpa .or. oskpb .or. oskpc)) then
                        onpsym = .true.
                        do 160 m = 1 , n4
                           it = m0(m)
                           ih = mi(it)
                           jh = mj(it)
                           if (jh.gt.ih) then
                              id = jh
                              jd = ih
                           else
                              id = ih
                              jd = jh
                           end if
                           if (.not.oskpa)
     +                         oskpa = (id.eq.ii .and. jd.eq.kk) .or.
     +                         (id.eq.jj .and. jd.eq.ll)
                           if (.not.oskpb)
     +                         oskpb = (id.eq.ii .and. jd.eq.ll) .or.
     +                         (id.eq.jj .and. jd.eq.kk)
                           if (oskpa .and. oskpb) go to 170
                           kh = mk(it)
                           if (kh.gt.ih) then
                              id = kh
                              kd = ih
                           else
                              id = ih
                              kd = kh
                           end if
                           if (.not.oskpc)
     +                         oskpc = (id.eq.ii .and. kd.eq.ll) .or.
     +                         (id.eq.jj .and. kd.eq.kk)
                           if (oskpa .and. oskpc) go to 180
                           if (oskpb .and. oskpc) go to 190
 160                    continue
                     end if
                     go to 200
 170                 oskpc = .true.
                     go to 200
 180                 oskpb = .true.
                     go to 200
 190                 oskpa = .true.
 200                 q4 = dble(nt)/dble(n4)
                     iexch = 1
                     ish = ii
                     jsh = jj
                     ksh = kk
                     lsh = ll
                     qq4 = q4
                     if (oskpa .and. onpsym) qq4 = qq4 + q4
                     if (oskpb .and. onpsym) qq4 = qq4 + q4
                     go to 230
 210                 if (.not.(oskpa)) then
                        iexch = 2
                        ish = ii
                        jsh = kk
                        ksh = jj
                        lsh = ll
                        qq4 = q4
                        if (oskpc .and. onpsym) qq4 = qq4 + q4
                        go to 230
                     end if
 220                 if (oskpb .or. oskpc) go to 240
                     iexch = 3
                     ish = ii
                     jsh = ll
                     ksh = jj
                     lsh = kk
                     qq4 = q4
 230                 norg = norgin(iexch)
                     call shells(gout(norg+1),1,ish,jsh,ksh,lsh,iexch)
                     call ijprim
                     if (nij.ne.0) then
                        call shells(gout(norg+1),2,ish,jsh,ksh,lsh,
     +                              iexch)
                        call genral(gout(norg+1))
                     else
                        call vclr(gout(norg+1),1,ngout)
                     end if
                     go to (210,220,240) , iexch
 240                 if (opk) then
                       if (odrf) then
                         call pkfilrf(ii,jj,kk,ll,oskpa,oskpb,oskpc,
     +                   onpsym,gout,q)
                       else
                         call pkfile(ii,jj,kk,ll,oskpa,oskpb,oskpc,
     +                   onpsym,gout,q)
                       endif
                     endif
                     call chkout(ii,jj,kk,ll,q(inull),q)
                     if (omaxb .or. tim.gt.timlim) go to 310
 250              continue
 260           continue
 270        continue
            time = cpulft(1)
c *** clear p sort buffers
 280        if (osortp) then
               call zsortp(q(igbf),q(iklbf),q(iiptbf),ijbase)
               do 290 ibuf = 1 , ngbf
                  call clrgbf(ibuf,q(igbf),q(iklbf),q(iiptbf))
 290           continue
               call blocki
            end if
c ***
 300     continue
         call final(q,q(inull),q(inull))
      end if
 310  continue
      if (osortp) then
        call gmem_free(iiptbf)
        call gmem_free(iklbf)
        call gmem_free(igbf)
      endif
      if (odrf) then
        if (itwoeps .eq. 1) then
          call gmem_free(i170)
        endif
        call gmem_free(i160)
        call gmem_free(i150)
        call gmem_free(i140)
        call gmem_free(i130)
        call gmem_free(i120)
        call gmem_free(i110)
      endif
      return
 6010 format (i4,3i5,1x,i10,i9,f11.2,f9.2)
      end
      subroutine pkints(q,iso,gout,nshels,outvv)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
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
      parameter (mxp2 = mxprms * mxprms)
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
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
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
c
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
c
      logical ospbas, onocnt, ochek, opdbas, opfbas, opgbas
      integer kad, istd
      common /ijlab/ kad(mxshel),ospbas,onocnt,istd,ochek,
     +               opdbas,opfbas,opgbas
c
      common/junk/cxyz(3,5625),aaa(21*mxp2),ijaaa(225)
c
      character*8 field
      character*24 version
      common /drfopt/ field,version
      logical oreact
      integer intdrf
      common /drfopti/ oreact,intdrf
c
      logical odrf
c
      integer ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep
      integer neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga
      integer idipcal,iclinte,iclintd,ieffpol,isodis,iclintr
      integer iextdip,itsolv,itermax,isolsav,imomsav,idisadd
      integer ianal,ineqex,ithole,irevdis,ihbond,igrppol
      integer maxci_drf
      common/drfpar1/ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep,
     +               neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga,
     +               idipcal,iclinte,iclintd,ieffpol,isodis,iclintr,
     +               iextdip,itsolv,itermax,isolsav,imomsav,idisadd,
     +               ianal,ineqex,ithole,irevdis,ihbond,igrppol,
     +               maxci_drf
c
      integer nxtpts, npol, nexp, natint, nshint, namb, nspec, ngran
      common/drfpar2/nxtpts,npol,nexp,natint,nshint,namb,nspec,ngran
c
      integer ndima, ndimb, nbem, ndim, nexp4, npol3, nwtr, nwtc, nzfa
      integer lomga, nomga, nchd, ngrnam, ngrnam2, nzfp, nzfn
      integer nneq, nneqrf, maxneq, neqdim, nqdim, nqcls
      common/drfpar3/ndima,ndimb,nbem,ndim,nexp4,npol3,nwtr,nwtc,nzfa,
     +               lomga,nomga,nchd,ngrnam,ngrnam2,nzfp,nzfn,
     +               nneq,nneqrf,maxneq,neqdim,nqdim,nqcls
c
      integer ifldin, ifldout, idrfout, modxza, irepopt
      integer iexpza, icmexp, iadexp, nodpe, igetden, iarfcal, iqmclr
      common/drfpar4/ifldin,ifldout,idrfout,modxza,irepopt,
     +               iexpza,icmexp,iadexp,nodpe,igetden,iarfcal,
     +               iqmclr
c
      real*8 gamdrf, dstmin, dstmax, hbondl, hbondr, afact, rfact
      real*8 cvgrel, agrpe, agrpm, agrpc, scffact, acur
      real*8 conci_drf
      common/drfpar5/gamdrf,dstmin,dstmax,hbondl,hbondr,afact,rfact,
     +               cvgrel,agrpe,agrpm,agrpc,scffact,acur,conci_drf
c
      integer iind, isur1, isur2, ilwt, ilvr, illur, ilindx
      common/drfpar6/iind,isur1,isur2,ilwt,ilvr,illur,ilindx
c
      integer i10, i20, i30, i40, i50, i60, i70, i80, i90 
      integer i100, i110, i120, i130, i140, i150, i160
      integer i170, i180
      common /scmprt/ i10,i20,i30,i40,i50,i60,i70,i80,i90,
     +          i100,i110,i120,i130,i140,i150,i160,i170,i180
c
c
      integer itwoeps,ioponly,itwosur,ixnsurf,ibemout,iuniout
      common /drfbem1/ itwoeps,ioponly,itwosur,ixnsurf,ibemout,iuniout
c
      integer leveli, levelo, nbem1, nbem2, nraw
      common /drfbem2/ leveli,levelo,nbem1,nbem2,nraw
c
      character*32 solnam
      real*8 radmax, eps1, eps2, spherad, rprobe, rprobej
      real*8 solrad, swidth, sdist
      real*8 kappa
      real*8 kappa1, kappa2, kappas
      common /drfbem3/ radmax,eps1,eps2,kappa1,kappa2,spherad,rprobe,
     +                 rprobej,solnam,solrad,swidth,sdist
c
      real*8 radat, cgrav
      common /drfbem4/ radat(maxat),cgrav(3)
c
c ***
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
c     logical opk,opank,oskips
c     common/pkfil/opk,opank,oskips
      logical opk,opank
      common/pkfil/opk,opank
c
      integer igt, jgt, kgt, lgt
      integer ibt, jbt, kbt, lbt, lenddd
      common /flips/ igt(3),jgt(3),kgt(3),lgt(3),
     +               ibt(mxprms),jbt(mxprms),kbt(mxprms),lbt(mxprms),
     +               lenddd
c
      logical osortp,oschw
      common/sortpdat/
     + igbf,iklbf,iiptbf,ijbase,ngbf,lengbf,osortp,oschw
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      integer icount_dlb, icount_slb
      common /par_dlb/ icount_dlb, icount_slb
c
      dimension iso(nshels,*),q(*),gout(*)
      dimension mi(48),mj(48),mk(48),m0(48)
      dimension norgin(3)
      dimension ib(4,4)
      data norgin/0,50625,101250/
      data ib/64,16,4,1,216,36,6,1,1000,100,10,1,
     +         3375,225,15,1 /
      data m25/25/
      inull = igmem_null()
cafc
      odrf = (field(5:) .eq. 'scf ') .and. (intdrf .ne. 0)
cxxx
c       allocate memory for drf matrices
c
      if (odrf) then
c-----  set acur for use in 2-electron reaction-field
c       routines
c
****    acurcy not passed into this routine
****    use value hardwired in scfinp (1.0d-05)
****    acur = acurcy
        acur = 1.0d-05
c  ovlap at i110, da10, 12, l2
c  dipx at  i120, da10, 53, l2
c  dipy at  i130, da10, 54, l2
c  dipz at  i140, da10, 55, l2
c  iexp at      i150, da31, 2, l2
c  omega(s) at  i160, da31, 50, nomga
c  omega(op) at i170, da31, 51, nomga (itwoeps = 1)
c
        i110 = igmem_alloc(nx)
        i120 = igmem_alloc(nx)
        i130 = igmem_alloc(nx)
        i140 = igmem_alloc(nx)
        i150 = igmem_alloc(2*nx)
        i160 = igmem_alloc(nomga)
        if (itwoeps .eq. 1) then
          i170 = igmem_alloc(nomga)
        endif
cafc  initialisation  rf contribution two-electron integrals
        call rftwin(q)
cafc
      endif
cafc
c
      nav = lenwrd()
      ltri = ikyp(nshels)
      oskipp = .false.
      if (opank) osortp = .false.
      ngout = 256
      lenddd = 4
      ii = 1
      if (opdbas) then
         ngout = 1296
         ii = 2
         lenddd = 6
      end if
      if (opfbas) then
         ngout = 10000
         lenddd = 10
         ii = 3
      end if
      if (opgbas) then
         ngout = 50625
         lenddd = 15
         ii = 4
      end if
      do 30 loop = 1 , 3
         igt(loop) = ib(1,ii)
         jgt(loop) = ib(2,ii)
         kgt(loop) = ib(3,ii)
         lgt(loop) = ib(4,ii)
 30   continue
      ii = 0
      do 40 loop = 1 , lenddd
         lbt(loop) = ii
         kbt(loop) = lbt(loop)*lenddd
         jbt(loop) = kbt(loop)*lenddd
         ibt(loop) = jbt(loop)*lenddd
         ii = ii + 1
 40   continue
      ischw = igmem_alloc(ltri)
      nschwz = 0
      if (osortp) then
c *** num*maxshellsize buffers length blocksize
         ngbf = num*6
         if (opfbas) ngbf = num*10
         if (opgbas) ngbf = num*15
c *** need mods to clrgbf to use anything other than 340 here
         lengbf = num2ep
         igbf = igmem_alloc(ngbf*lengbf)
         itmp = (ngbf*lengbf+1)/nav
         iklbf = igmem_alloc(itmp)
         itmp = (ngbf+1)/nav
         iiptbf = igmem_alloc(itmp)
      end if
c *** zero out p sort buffer counts
      if (osortp) call setsto(ngbf,0,q(iiptbf))
c *** read in ints for schwarz inequality test
      call secget(isect(421),m25,iblk25)
      call rdedx(q(ischw),ltri,iblk25,idaf)
      dlncutoff = dlog(cutoff)
c
cafc  initialisation  rf contribution two-electron integrals
c     if (odrf) call rftwin(q)
cafc
      time = cpulft(1)
      tim0 = time
      tim1 = time
      ist0 = ist
      jst0 = jst
      kst0 = kst
      lst0 = lst
      if (ist0.le.nshels) then
         do 300 ii = ist0 , nshels
c *** set ijbase to triangle of first basis function in shell
c *** and check integrity of counters
            if (osortp) then
               ijbase = (kloc(ii)*(kloc(ii)-1))/2
               if (icount.ne.1) call caserr('icount.ne.1 for new ii')
            end if
            kadi = kad(ii)
            dt0 = time - tim0
            dt1 = time - tim1
            tim1 = time
            if (outv) write (iwr,6010) ii , jst0 , kst0 , lst0 , nrec ,
     +                                icount , dt1 , dt0
            do 60 it = 1 , nt
               id = iso(ii,it)
               if (id.gt.ii) go to 280
               mi(it) = id
 60         continue
            j0 = jst0
            do 270 jj = j0 , ii
               jst0 = 1
               kadij = kadi + kad(jj)
               do 80 it = 1 , nt
                  id = mi(it)
                  jd = iso(jj,it)
                  mj(it) = jd
                  if (id.lt.jd) then
                     nd = id
                     id = jd
                     jd = nd
                  end if
                  if (id.lt.ii) then
                  else if (id.eq.ii) then
                     if (jd.gt.jj) go to 270
                  else
                     go to 270
                  end if
 80            continue
               k0 = kst0
               do 260 kk = k0 , jj
                  kst0 = 1
                  kadijk = kadij + kad(kk)
                  do 110 it = 1 , nt
                     id = mi(it)
                     jd = mj(it)
                     kd = iso(kk,it)
                     mk(it) = kd
 100                 if (id.lt.jd) then
                        nd = id
                        id = jd
                        jd = nd
                     end if
                     if (jd.lt.kd) then
                        nd = jd
                        jd = kd
                        kd = nd
                        go to 100
                     else if (id.lt.ii) then
                     else if (id.eq.ii) then
                        if (jd.lt.jj) then
                        else if (jd.eq.jj) then
                           if (kd.gt.kk) go to 260
                        else
                           go to 260
                        end if
                     else
                        go to 260
                     end if
 110              continue
                  l0 = lst0
                  do 250 ll = l0 , kk
                     lst0 = 1
                     if (kadijk+kad(ll).ge.0) go to 250
                     n4 = 0
                     do 150 it = 1 , nt
                        id = mi(it)
                        jd = mj(it)
                        kd = mk(it)
                        ld = iso(ll,it)
 130                    if (id.lt.jd) then
                           nd = id
                           id = jd
                           jd = nd
                        end if
 140                    if (jd.lt.kd) then
                           nd = jd
                           jd = kd
                           kd = nd
                           go to 130
                        else if (kd.lt.ld) then
                           nd = kd
                           kd = ld
                           ld = nd
                           go to 140
                        else if (id.lt.ii) then
                        else if (id.eq.ii) then
                           if (jd.lt.jj) then
                           else if (jd.eq.jj) then
                              if (kd.lt.kk) then
                              else if (kd.eq.kk) then
                                 if (ld.lt.ll) then
                                 else if (ld.eq.ll) then
                                    n4 = n4 + 1
                                    m0(n4) = it
                                 else
                                    go to 250
                                 end if
                              else
                                 go to 250
                              end if
                           else
                              go to 250
                           end if
                        else
                           go to 250
                        end if
 150                 continue
                     oskpa = jj.eq.kk
                     oskpb = (ii.eq.kk) .or. (jj.eq.ll)
                     oskpc = (ii.eq.jj) .or. (kk.eq.ll)
                     onpsym = .false.
                     if (.not.(oskpa .or. oskpb .or. oskpc)) then
                        onpsym = .true.
                        do 160 m = 1 , n4
                           it = m0(m)
                           ih = mi(it)
                           jh = mj(it)
                           if (jh.gt.ih) then
                              id = jh
                              jd = ih
                           else
                              id = ih
                              jd = jh
                           end if
                           if (.not.oskpa)
     +                         oskpa = (id.eq.ii .and. jd.eq.kk) .or.
     +                         (id.eq.jj .and. jd.eq.ll)
                           if (.not.oskpb)
     +                         oskpb = (id.eq.ii .and. jd.eq.ll) .or.
     +                         (id.eq.jj .and. jd.eq.kk)
                           if (oskpa .and. oskpb) go to 170
                           kh = mk(it)
                           if (kh.gt.ih) then
                              id = kh
                              kd = ih
                           else
                              id = ih
                              kd = kh
                           end if
                           if (.not.oskpc)
     +                         oskpc = (id.eq.ii .and. kd.eq.ll) .or.
     +                         (id.eq.jj .and. kd.eq.kk)
                           if (oskpa .and. oskpc) go to 180
                           if (oskpb .and. oskpc) go to 190
 160                    continue
                     end if
                     go to 200
 170                 oskpc = .true.
                     go to 200
 180                 oskpb = .true.
                     go to 200
 190                 oskpa = .true.
 200                 q4 = dble(nt)/dble(n4)
                     iexch = 1
                     ish = ii
                     jsh = jj
                     ksh = kk
                     lsh = ll
                     qq4 = q4
                     if (oskpa .and. onpsym) qq4 = qq4 + q4
                     if (oskpb .and. onpsym) qq4 = qq4 + q4
                     go to 230
 210                 if (.not.(oskpa)) then
                        iexch = 2
                        ish = ii
                        jsh = kk
                        ksh = jj
                        lsh = ll
                        qq4 = q4
                        if (oskpc .and. onpsym) qq4 = qq4 + q4
                        go to 230
                     end if
 220                 if (oskpb .or. oskpc) go to 240
                     iexch = 3
                     ish = ii
                     jsh = ll
                     ksh = jj
                     lsh = kk
                     qq4 = q4
 230                 norg = norgin(iexch)
                     ijij = iky(ish) + jsh + ischw -1
                     klkl = iky(ksh) + lsh + ischw -1
                     test = q(ijij) + q(klkl)
                     oskipp = test.lt.dlncutoff
                     if (oskipp) then
                        nschwz = nschwz + 1
                        call vclr(gout(norg+1),1,ngout)
                        go to 320
                     endif
                     call shells(gout(norg+1),1,ish,jsh,ksh,lsh,iexch)
                     call ijprim
                     if (nij.ne.0) then
                        call shells(gout(norg+1),2,ish,jsh,ksh,lsh,
     +                              iexch)
                        call genral(gout(norg+1))
                     else
                        call vclr(gout(norg+1),1,ngout)
                     end if
 320                 go to (210,220,240) , iexch
 240                 if (opk) then
                       if (odrf) then
                          call pkfilrf(ii,jj,kk,ll,oskpa,oskpb,oskpc,
     +                    onpsym,gout,q)
                       else
                          call pkfile(ii,jj,kk,ll,oskpa,oskpb,oskpc,
     +                    onpsym,gout,q)
                       endif
                     endif
                     call chkout(ii,jj,kk,ll,q(inull),q)
                     if (omaxb .or. tim.gt.timlim) go to 310
 250              continue
 260           continue
 270        continue
            time = cpulft(1)
c *** clear p sort buffers
 280        if (osortp) then
               call zsortp(q(igbf),q(iklbf),q(iiptbf),ijbase)
               do 290 ibuf = 1 , ngbf
                  call clrgbf(ibuf,q(igbf),q(iklbf),q(iiptbf))
 290           continue
               call blocki
            end if
c ***
 300     continue
         call final(q,q(inull),q(inull))
         if(outvv) write(iwr,6020) nschwz
      end if
 310  if(osortp) then
        call gmem_free(iiptbf)
        call gmem_free(iklbf)
        call gmem_free(igbf)
      endif
      call gmem_free(ischw)
      if (odrf) then
        if (itwoeps .eq. 1) then
          call gmem_free(i170)
        endif
        call gmem_free(i160)
        call gmem_free(i150)
        call gmem_free(i140)
        call gmem_free(i130)
        call gmem_free(i120)
        call gmem_free(i110)
      endif
      return
 6010 format (i4,3i5,1x,i10,i9,f11.2,f9.2)
 6020 format(/1x,'schwarz inequality test skipped ',i10,
     + ' integral blocks')
      end
      subroutine s0000(g)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
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
      parameter (mxp2 = mxprms * mxprms)
c
      real*8 ag, csa, cpa, cda, cfa, cga
      real*8 bg, csb, cpb, cdb, cfb, cgb
      real*8 cgg, csc, cpc, cdc, cfc, cgc
      real*8 dg, csd, cpd ,cdd, cfd, cgd
      real*8 pi, qi, ri, pj, qj, rj, rri, pk, qk, rk, pl, ql, rl, rrk
      integer nga, ngb, ngc, ngd
      common /shlinf/ ag(mxprms),csa(mxprms),cpa(mxprms),cda(mxprms),
     +                cfa(mxprms),cga(mxprms),
     +                bg(mxprms),csb(mxprms),cpb(mxprms),cdb(mxprms),
     +                cfb(mxprms),cgb(mxprms),
     +                cgg(mxprms),csc(mxprms),cpc(mxprms),cdc(mxprms),
     +                cfc(mxprms),cgc(mxprms),
     +                dg(mxprms),csd(mxprms),cpd(mxprms),cdd(mxprms),
     +                cfd(mxprms),cgd(mxprms),
     +                pi,qi,ri,pj,qj,rj,rri,pk,qk,rk,pl,ql,rl,rrk,
     +                nga,ngb,ngc,ngd
c
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
      common/junk/cxyz(3,5625),
     + a(mxp2),r(mxp2),p1(mxp2),q1(mxp2),r1(mxp2),dij(16*mxp2),
     + ijd(225)
c
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
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
      logical oianj, okanl, oident, omisc, oham, opdipd, omp2
      integer ipos1, ipos2
      common /misc/ oianj,okanl,oident,omisc,
     +            oham,opdipd,omp2,ipos1,ipos2
c
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
      subroutine shells(gout,nelec,ish,jsh,ksh,lsh,iexch)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
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
      parameter (mxp2 = mxprms * mxprms)
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
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
c
      real*8 pp, u, w, dji
      integer nroots, nroot2
      common /root/ pp,u(12),w(12),nroots,nroot2,dji(60)
c
      common/junk/cxyz(3,5625),aaa(21*mxp2),ijaaa(225)
c
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
c
      integer ik, ijx, ijy, ijz, klx, kly, klz
      integer ijgt, klgt
      common /indez/ ijgt(225),ijx(225),ijy(225),ijz(225),
     +       ik(225),klgt(225),klx(225),kly(225),klz(225)
c
c
      real*8 ag, csa, cpa, cda, cfa, cga
      real*8 bg, csb, cpb, cdb, cfb, cgb
      real*8 cgg, csc, cpc, cdc, cfc, cgc
      real*8 dg, csd, cpd ,cdd, cfd, cgd
      real*8 pi, qi, ri, pj, qj, rj, rri, pk, qk, rk, pl, ql, rl, rrk
      integer nga, ngb, ngc, ngd
      common /shlinf/ ag(mxprms),csa(mxprms),cpa(mxprms),cda(mxprms),
     +                cfa(mxprms),cga(mxprms),
     +                bg(mxprms),csb(mxprms),cpb(mxprms),cdb(mxprms),
     +                cfb(mxprms),cgb(mxprms),
     +                cgg(mxprms),csc(mxprms),cpc(mxprms),cdc(mxprms),
     +                cfc(mxprms),cgc(mxprms),
     +                dg(mxprms),csd(mxprms),cpd(mxprms),cdd(mxprms),
     +                cfd(mxprms),cgd(mxprms),
     +                pi,qi,ri,pj,qj,rj,rri,pk,qk,rk,pl,ql,rl,rrk,
     +                nga,ngb,ngc,ngd
c
c
      logical oianj, okanl, oident, omisc, oham, opdipd, omp2
      integer ipos1, ipos2
      common /misc/ oianj,okanl,oident,omisc,
     +            oham,opdipd,omp2,ipos1,ipos2
c
c
      integer igt, jgt, kgt, lgt
      integer ibt, jbt, kbt, lbt, lenddd
      common /flips/ igt(3),jgt(3),kgt(3),lgt(3),
     +               ibt(mxprms),jbt(mxprms),kbt(mxprms),lbt(mxprms),
     +               lenddd
c
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
         oident = ish.eq.ksh .and. jsh.eq.lsh
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
         max = kl
         do 60 i = 1 , ij
            if (oident) max = i
            ik(i) = max
 60      continue
         ijkl = ij*kl
         if (oident) ijkl = ij*(ij+1)/2
c
c     zero integral storage
c
         do 65 i = 1,ij
         ngij = ijgt(i)
         do 65 k = 1,ik(i)
  65     gout(ngij + klgt(k)) = 0.0d0
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
      subroutine spdint(gout)
c
c     ----- form integrals over functions -----
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
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
      parameter (mxp2 = mxprms * mxprms)
c
      real*8 pp, u, w, dji
      integer nroots, nroot2
      common /root/ pp,u(12),w(12),nroots,nroot2,dji(60)
c
c
      integer ik, ijx, ijy, ijz, klx, kly, klz
      integer ijgt, klgt
      common /indez/ ijgt(225),ijx(225),ijy(225),ijz(225),
     +       ik(225),klgt(225),klx(225),kly(225),klz(225)
c
      common/junk/pin(5625),qin(5625),rin(5625),aaa(21*mxp2),ijaaa(225)
c
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
c
      real*8 dij, dkl
      common /denss/ dij(225),dkl(225)
c
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
      subroutine sskl(gout)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
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
      parameter (mxp2 = mxprms * mxprms)
c
      real*8 ag, csa, cpa, cda, cfa, cga
      real*8 bg, csb, cpb, cdb, cfb, cgb
      real*8 cgg, csc, cpc, cdc, cfc, cgc
      real*8 dg, csd, cpd ,cdd, cfd, cgd
      real*8 pi, qi, ri, pj, qj, rj, rri, pk, qk, rk, pl, ql, rl, rrk
      integer nga, ngb, ngc, ngd
      common /shlinf/ ag(mxprms),csa(mxprms),cpa(mxprms),cda(mxprms),
     +                cfa(mxprms),cga(mxprms),
     +                bg(mxprms),csb(mxprms),cpb(mxprms),cdb(mxprms),
     +                cfb(mxprms),cgb(mxprms),
     +                cgg(mxprms),csc(mxprms),cpc(mxprms),cdc(mxprms),
     +                cfc(mxprms),cgc(mxprms),
     +                dg(mxprms),csd(mxprms),cpd(mxprms),cdd(mxprms),
     +                cfd(mxprms),cgd(mxprms),
     +                pi,qi,ri,pj,qj,rj,rri,pk,qk,rk,pl,ql,rl,rrk,
     +                nga,ngb,ngc,ngd
c
      common/junk/pint(5625),qint(5625),rint(5625),
     + aa(mxp2),r(mxp2),p1(mxp2),q1(mxp2),r1(mxp2),dd(16*mxp2),
     + ijd(225)
c
      integer ik, ijx, ijy, ijz, klx, kly, klz
      integer ijgt, klgt
      common /indez/ ijgt(225),ijx(225),ijy(225),ijz(225),
     +       ik(225),klgt(225),klx(225),kly(225),klz(225)
c
c
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
c
      real*8 pp, u, w, dji
      integer nroots, nroot2
      common /root/ pp,u(12),w(12),nroots,nroot2,dji(60)
c
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
c
      logical oianj, okanl, oident, omisc, oham, opdipd, omp2
      integer ipos1, ipos2
      common /misc/ oianj,okanl,oident,omisc,
     +            oham,opdipd,omp2,ipos1,ipos2
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
      integer in, kn, ni, nj, nk, nl, nmax, mmax
      integer ij1, ij2, kl1, kl2
      real*8 bp01, b00, b10, pcp00, pc00, qcp00, qc00, rcp00, rc00
      real*8 f00, dxij, dyij, dzij, dxkl, dykl, dzkl
      common /setint/ in(9),kn(9),ni,nj,nk,nl,nmax,mmax,
     +                bp01,b00,b10,pcp00,pc00,qcp00,qc00,rcp00,rc00,
     +                f00,dxij,dyij,dzij,dxkl,dykl,dzkl,
     +                ij1,ij2,kl1,kl2
c
c
      real*8 dij, dkl
      common /denss/ dij(225),dkl(225)
c
      dimension gout(*)
      data kln1,kln2 /5,1/
      data sqrt3 /1.73205080756888d0/
      data sqrt5 /2.23606797749979d0/
      data sqrt7 /2.64575131106459d0/
      data pi252 /34.986836655250d0/
      data dzero,pt5,done /0.0d0,0.5d0,1.0d0/
      factor = pi252*qq4
      onorm = normf.ne.1 .or. normp.ne.1
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
c     ----- select expansion center for -xyz- integrals -----
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
                        call rt123
                     else if (nroots.eq.4) then
                        call roots4
                     else
                        call roots5
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
      subroutine xyzint
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c
      integer in, kn, ni, nj, nk, nl, nmax, mmax
      integer ij1, ij2, kl1, kl2
      real*8 bp01, b00, b10, pcp00, pc00, qcp00, qc00, rcp00, rc00
      real*8 f00, dxij, dyij, dzij, dxkl, dykl, dzkl
      common /setint/ in(9),kn(9),ni,nj,nk,nl,nmax,mmax,
     +                bp01,b00,b10,pcp00,pc00,qcp00,qc00,rcp00,rc00,
     +                f00,dxij,dyij,dzij,dxkl,dykl,dzkl,
     +                ij1,ij2,kl1,kl2
c
      common/junk/pint(5625),qint(5625),rint(5625)
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
      subroutine ssprim
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
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
      parameter (mxp2 = mxprms * mxprms)
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
c
      real*8 ag, csa, cpa, cda, cfa, cga
      real*8 bg, csb, cpb, cdb, cfb, cgb
      real*8 cgg, csc, cpc, cdc, cfc, cgc
      real*8 dg, csd, cpd ,cdd, cfd, cgd
      real*8 pi, qi, ri, pj, qj, rj, rri, pk, qk, rk, pl, ql, rl, rrk
      integer nga, ngb, ngc, ngd
      common /shlinf/ ag(mxprms),csa(mxprms),cpa(mxprms),cda(mxprms),
     +                cfa(mxprms),cga(mxprms),
     +                bg(mxprms),csb(mxprms),cpb(mxprms),cdb(mxprms),
     +                cfb(mxprms),cgb(mxprms),
     +                cgg(mxprms),csc(mxprms),cpc(mxprms),cdc(mxprms),
     +                cfc(mxprms),cgc(mxprms),
     +                dg(mxprms),csd(mxprms),cpd(mxprms),cdd(mxprms),
     +                cfd(mxprms),cgd(mxprms),
     +                pi,qi,ri,pj,qj,rj,rri,pk,qk,rk,pl,ql,rl,rrk,
     +                nga,ngb,ngc,ngd
c
c
      logical oianj, okanl, oident, omisc, oham, opdipd, omp2
      integer ipos1, ipos2
      common /misc/ oianj,okanl,oident,omisc,
     +            oham,opdipd,omp2,ipos1,ipos2
c
c
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
      common/junk/cxyz(3,5625),
     +   a(mxp2),r(mxp2),p1(mxp2),q1(mxp2),r1(mxp2),dij(16*mxp2),
     +   ijd(225)
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
      subroutine ver_intega(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/intega.m,v $
     +     "/
      data revision /"$Revision: 5958 $"/
      data date /"$Date: 2009-05-07 13:08:36 +0200 (Thu, 07 May 2009) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
cafc
      subroutine rftwin(q)
c
      implicit real*8  (a-h,o-w),integer  (i-n)
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
      character*8 field
      character*24 version
      common /drfopt/ field,version
      logical oreact
      integer intdrf
      common /drfopti/ oreact,intdrf
c
cafc
c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     this file contains the parameter settings for piet van duijnen's
c     reactionfield programme together with the sizes needed for the
c     v. duijnen - gamess interface as developed by 
c     nederkoorn and v. lenthe
c
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c     ==================================================================
c
c     parameter setting for reaction field programme van duijnen:
c
c=======================================================================
c
c     
      parameter (mxat = 10000)
      parameter (mxpol = 10000)
      parameter (mxpts = 10000)
      parameter (mxspe = 50)
      parameter (mxamb = 10)
      parameter (mxgrp = 500)
      parameter (mxnum = 1024)
      parameter (mxsh = 512)
      parameter (mxgs = 2048)
      parameter (mxnpts = 20)
      parameter (mxneq = 255)
      parameter (mcwrkmx = 10)
c
c ==================================================================
c
c     parameter setting new drf routines (drf/dimpar):
c
c===================================================================
c
      parameter (mxst = 5)
      parameter (mxgran = 10)
cmw      parameter (mxnpts = 50)
      parameter (mxgrpar = (mxgran*(mxgran+1))/2)
c
c     ==================================================================
c
c     parameter setting for v. duijnen - gamess interface:
c
c=======================================================================
c
      parameter (ncatom=500)
c
c  ncatom is maximal allowed classical atoms in zmat
c  total allowed atoms in zmat = maxnz + ncatom
c
c  maxrfa is maximal number of classical atoms in zmat + coordinates
c  total maximal atoms allowed is maxrfa + maxat
c
      parameter (maxrfa=mxpts)
c
c  maxamb is maximal number of ambigous atoms. it is set equal
c  to maximal z-cards possible in gamess input , i.e. 50
c
      parameter (maxamb=mxamb)
c
c=============================================
c for the pes-routines
c=============================================
c
      parameter (mxgr1=128, mxpes=100)
c
c=============================================
c for the drfexpc routines
c=============================================
c
      parameter (maxex=1000)
c
c=============================================
c for the drfsurf routine
c=============================================
c
      parameter (maxinv=10,maxnrc=10000)
c      parameter (maxnrc=10000)
c
c=============================================
c for the analrf routine
c=============================================
c
c
c
      integer ir, iwr, ipnch,ipadiofil
      common /iofile/ ir,iwr,ipnch,ipadiofil(20)
c
      integer list
      common /hlistng/ list
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
      integer maxcor, maxlcm
      common /hmemory/ maxcor,maxlcm
c
c
      character*8 wfntyp
      common /wfnopt/ wfntyp
c
      character*8 runtyp
      common /runopt/ runtyp
c
      character*8 hndtyp
      common /hndopt/ hndtyp
c
c
      integer ipople, ihondo
      common /inttyp/ ipople,ihondo
c
      real*8 qint, valint
      integer jcint, out
      common /intprt/ qint(3),valint(3),jcint(16),out
c
      integer nopk, nok, nosqur
      common /intopt/ nopk,nok,nosqur
c
      integer nintmx
      common /intfil/ nintmx
c
      integer jkskip
      common /intcal/ jkskip
c
c
c
       integer idafh, navh, ioda
       common /hdafile/ idafh,navh,ioda(2,1000)
c
c
c
      integer ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep
      integer neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga
      integer idipcal,iclinte,iclintd,ieffpol,isodis,iclintr
      integer iextdip,itsolv,itermax,isolsav,imomsav,idisadd
      integer ianal,ineqex,ithole,irevdis,ihbond,igrppol
      integer maxci_drf
      common/drfpar1/ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep,
     +               neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga,
     +               idipcal,iclinte,iclintd,ieffpol,isodis,iclintr,
     +               iextdip,itsolv,itermax,isolsav,imomsav,idisadd,
     +               ianal,ineqex,ithole,irevdis,ihbond,igrppol,
     +               maxci_drf
c
      integer nxtpts, npol, nexp, natint, nshint, namb, nspec, ngran
      common/drfpar2/nxtpts,npol,nexp,natint,nshint,namb,nspec,ngran
c
      integer ndima, ndimb, nbem, ndim, nexp4, npol3, nwtr, nwtc, nzfa
      integer lomga, nomga, nchd, ngrnam, ngrnam2, nzfp, nzfn
      integer nneq, nneqrf, maxneq, neqdim, nqdim, nqcls
      common/drfpar3/ndima,ndimb,nbem,ndim,nexp4,npol3,nwtr,nwtc,nzfa,
     +               lomga,nomga,nchd,ngrnam,ngrnam2,nzfp,nzfn,
     +               nneq,nneqrf,maxneq,neqdim,nqdim,nqcls
c
      integer ifldin, ifldout, idrfout, modxza, irepopt
      integer iexpza, icmexp, iadexp, nodpe, igetden, iarfcal, iqmclr
      common/drfpar4/ifldin,ifldout,idrfout,modxza,irepopt,
     +               iexpza,icmexp,iadexp,nodpe,igetden,iarfcal,
     +               iqmclr
c
      real*8 gamdrf, dstmin, dstmax, hbondl, hbondr, afact, rfact
      real*8 cvgrel, agrpe, agrpm, agrpc, scffact, acur
      real*8 conci_drf
      common/drfpar5/gamdrf,dstmin,dstmax,hbondl,hbondr,afact,rfact,
     +               cvgrel,agrpe,agrpm,agrpc,scffact,acur,conci_drf
c
      integer iind, isur1, isur2, ilwt, ilvr, illur, ilindx
      common/drfpar6/iind,isur1,isur2,ilwt,ilvr,illur,ilindx
c
      integer idafdrf, navdrf, iodadrf
      common /drfdaf/ idafdrf,navdrf,iodadrf(2,255)
c
c
      integer i10, i20, i30, i40, i50, i60, i70, i80, i90 
      integer i100, i110, i120, i130, i140, i150, i160
      integer i170, i180
      common /scmprt/ i10,i20,i30,i40,i50,i60,i70,i80,i90,
     +          i100,i110,i120,i130,i140,i150,i160,i170,i180
c
c
c
      integer itwoeps,ioponly,itwosur,ixnsurf,ibemout,iuniout
      common /drfbem1/ itwoeps,ioponly,itwosur,ixnsurf,ibemout,iuniout
c
      integer leveli, levelo, nbem1, nbem2, nraw
      common /drfbem2/ leveli,levelo,nbem1,nbem2,nraw
c
      character*32 solnam
      real*8 radmax, eps1, eps2, spherad, rprobe, rprobej
      real*8 solrad, swidth, sdist
      real*8 kappa
      real*8 kappa1, kappa2, kappas
      common /drfbem3/ radmax,eps1,eps2,kappa1,kappa2,spherad,rprobe,
     +                 rprobej,solnam,solrad,swidth,sdist
c
      real*8 radat, cgrav
      common /drfbem4/ radat(maxat),cgrav(3)
c
c
c
      integer nsamp, nblock, nmoves, ncheck
      common /mcrunp1/ nsamp,nblock,nmoves,ncheck
c
      real*8 darot, dtrans, temp, onekt, excld, delmax, enmin
      common /mcrunp2/ darot,dtrans,temp,onekt,excld,delmax,enmin
c
      logical mcupdt,secstat,excess
      integer notrans, norot, imcout, iseed, imcst, iactst
      common /mcrunp3/ notrans,norot,imcout,iseed,imcst,iactst,
     +               mcupdt,secstat,excess
c
      character*16 outfor
      common /mcrunp4/ outfor
c
      real*8 ratmin, ratmax, amxrot, amxtrn, amnrot, amntrn
      common /mcrunp5/ ratmin,ratmax,amxrot,amxtrn,amnrot,amntrn
c
      real*8 gammc
      common /mcrunp6/ gammc(5)
c
c
c
      real*8 val
      integer i, j, k, l
      common /drfint/ i,j,k,l,val
c
c
      dimension q(*)
c
c  allocate mem
c     call setscm(i10)
c     i20 = i10 + nchd
c     i30 = i20 + nchd
c     i40 = i30 + nchd
c     i50 = i40 + nchd
c     i60 = i50 + 2*nchd
c     i70 = i60 + nomga
c
c     last = i10 + i70 + itwoeps*nomga
c     call cmem(loadcm)
c     length = last - i10
c     loc10 = loccm(q(i10))
c     need = loc10 + length
c     call setc(need)
c
c  read from da10 and da31
      call daread(idafh,ioda,q(i110),nchd,12)
      call daread(idafh,ioda,q(i120),nchd,53)
      call daread(idafh,ioda,q(i130),nchd,54)
      call daread(idafh,ioda,q(i140),nchd,55)
      call daread(idafdrf,iodadrf,q(i150),nchd,2)
      call daread(idafdrf,iodadrf,q(i160),nomga,50)
      if (itwoeps .ne. 0) call
     +  daread(idafdrf,iodadrf,q(i170),nomga,51)
      return
      end
c
      subroutine qoutdrf(q,g)
cafc
c  ordens and stores two-electron integrals including
c  the reaction field contribution; analogous to qout
c
      implicit real*8  (a-h,p-w),integer  (i-n)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y),logical (o)
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
c  i,j,k,l:
c
      real*8 val
      integer i, j, k, l
      common /drfint/ i,j,k,l,val
c
c  scffact:
c
      integer ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep
      integer neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga
      integer idipcal,iclinte,iclintd,ieffpol,isodis,iclintr
      integer iextdip,itsolv,itermax,isolsav,imomsav,idisadd
      integer ianal,ineqex,ithole,irevdis,ihbond,igrppol
      integer maxci_drf
      common/drfpar1/ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep,
     +               neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga,
     +               idipcal,iclinte,iclintd,ieffpol,isodis,iclintr,
     +               iextdip,itsolv,itermax,isolsav,imomsav,idisadd,
     +               ianal,ineqex,ithole,irevdis,ihbond,igrppol,
     +               maxci_drf
c
      integer nxtpts, npol, nexp, natint, nshint, namb, nspec, ngran
      common/drfpar2/nxtpts,npol,nexp,natint,nshint,namb,nspec,ngran
c
      integer ndima, ndimb, nbem, ndim, nexp4, npol3, nwtr, nwtc, nzfa
      integer lomga, nomga, nchd, ngrnam, ngrnam2, nzfp, nzfn
      integer nneq, nneqrf, maxneq, neqdim, nqdim, nqcls
      common/drfpar3/ndima,ndimb,nbem,ndim,nexp4,npol3,nwtr,nwtc,nzfa,
     +               lomga,nomga,nchd,ngrnam,ngrnam2,nzfp,nzfn,
     +               nneq,nneqrf,maxneq,neqdim,nqdim,nqcls
c
      integer ifldin, ifldout, idrfout, modxza, irepopt
      integer iexpza, icmexp, iadexp, nodpe, igetden, iarfcal, iqmclr
      common/drfpar4/ifldin,ifldout,idrfout,modxza,irepopt,
     +               iexpza,icmexp,iadexp,nodpe,igetden,iarfcal,
     +               iqmclr
c
      real*8 gamdrf, dstmin, dstmax, hbondl, hbondr, afact, rfact
      real*8 cvgrel, agrpe, agrpm, agrpc, scffact, acur
      real*8 conci_drf
      common/drfpar5/gamdrf,dstmin,dstmax,hbondl,hbondr,afact,rfact,
     +               cvgrel,agrpe,agrpm,agrpc,scffact,acur,conci_drf
c
      integer iind, isur1, isur2, ilwt, ilvr, illur, ilindx
      common/drfpar6/iind,isur1,isur2,ilwt,ilvr,illur,ilindx
c  itwoeps:
c
      integer itwoeps,ioponly,itwosur,ixnsurf,ibemout,iuniout
      common /drfbem1/ itwoeps,ioponly,itwosur,ixnsurf,ibemout,iuniout
c
      integer leveli, levelo, nbem1, nbem2, nraw
      common /drfbem2/ leveli,levelo,nbem1,nbem2,nraw
c
      character*32 solnam
      real*8 radmax, eps1, eps2, spherad, rprobe, rprobej
      real*8 solrad, swidth, sdist
      real*8 kappa
      real*8 kappa1, kappa2, kappas
      common /drfbem3/ radmax,eps1,eps2,kappa1,kappa2,spherad,rprobe,
     +                 rprobej,solnam,solrad,swidth,sdist
c
      real*8 radat, cgrav
      common /drfbem4/ radat(maxat),cgrav(3)
c
c  i110,i120,i130,i140,i150,i160,i170:
c
      integer i10, i20, i30, i40, i50, i60, i70, i80, i90 
      integer i100, i110, i120, i130, i140, i150, i160
      integer i170, i180
      common /scmprt/ i10,i20,i30,i40,i50,i60,i70,i80,i90,
     +          i100,i110,i120,i130,i140,i150,i160,i170,i180
c
c
c
      logical ospbas, onocnt, ochek, opdbas, opfbas, opgbas
      integer kad, istd
      common /ijlab/ kad(mxshel),ospbas,onocnt,istd,ochek,
     +               opdbas,opfbas,opgbas
c
      common/craypk/integ(680)
      common/mapper/iky(maxorb),ikyp(4,maxorb),i4096(maxorb)
      common/iofile/ ir,iw,ip,main,iblkm,idaf,iblkd,num8,iblk21,
     * iblkqa,iblkpa,iblkea,iblkqb,iblkpb,iblkeb,iblkg,iblkhs,
     * iblksp(3),nav
      common/shlt/tol,cutoff,icount,ic4,out
      common/restar/nprint,itol,icut,normf,normp,nopk,
     *irest,nrec,omxblk,ist,jst,kst,lst,
     *nintmx,nindmx,intg76
      common/shlnos/qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl
      common/misc/oianj,okanl,osame
      common/blkin/gout(510),nword
c
      integer ik, ijx, ijy, ijz, klx, kly, klz
      integer ijgt, klgt
      common /indez/ ijgt(225),ijx(225),ijy(225),ijz(225),
     +       ik(225),klgt(225),klx(225),kly(225),klz(225)
c
      dimension q(*),g(*)
c
c     ----- pack the 4 indices of two integrals into one word
c     ----- write label + integral on mainfile
c
      ijn = 0
      jmax = maxj
      do 2600 iii = mini,maxi
      if (oianj) jmax = iii
      i1 = loci + iii
      ipack = i4096(i1)
      do 2400 jjj = minj,jmax
      ijn = ijn + 1
      n1 = ijgt(ijn)
      i2 = locj + jjj
      if(i1-i2)10, 20, 20
 10   lab1 = i4096(i2) + i1
      go to 30
 20   lab1 = ipack + i2
 30   lmax = maxl
      kln = 0
      do 2200 kkk = mink,maxk
      if (okanl) lmax = kkk
      i3 = lock + kkk
      kpack = i4096(i3)
      do 2000 lll = minl,lmax
      kln = kln + 1
      if (osame .and. kln .gt. ijn) go to 2400
      nn = n1+klgt(kln)
      i4 = locl + lll
c
      i = i1
      j = i2
      k = i3
      l = i4
      val = g(nn) + scffact * rftwog(q(i110),q(i120),q(i130),
     +                               q(i140),q(i150),q(i160))
      if (itwoeps .ne. 0) then
        val = g(nn) + scffact * rftwog(q(i110),q(i120),q(i130),
     +                                 q(i140),q(i150),q(i170))
      endif
c
      if (dabs(val) .lt. cutoff) go to 2000
      if(i3-i4) 40, 50, 50
 40   lab2 = i4096(i4) + i3
      go to 60
 50   lab2 = kpack + i4
 60   gout(icount) = val
      integ(ic4  )=max(lab1,lab2)
      integ(ic4+1)=min(lab1,lab2)
      ic4=ic4+2
      icount = icount+1
      if (icount .le. nintmx) go to 2000
      call blocki
      if(omxblk)go to 261
 2000 continue
 2200 continue
 2400 continue
 2600 continue
 261  return
      end
c
      function rftwog(s,dx,dy,dz,iexpc,omega)
cafc
c  function value == two-electron contribution of reaction field
c
      implicit real*8  (a-h,o-w),integer  (i-n)
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
cafc
c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     this file contains the parameter settings for piet van duijnen's
c     reactionfield programme together with the sizes needed for the
c     v. duijnen - gamess interface as developed by 
c     nederkoorn and v. lenthe
c
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c     ==================================================================
c
c     parameter setting for reaction field programme van duijnen:
c
c=======================================================================
c
c     
      parameter (mxat = 10000)
      parameter (mxpol = 10000)
      parameter (mxpts = 10000)
      parameter (mxspe = 50)
      parameter (mxamb = 10)
      parameter (mxgrp = 500)
      parameter (mxnum = 1024)
      parameter (mxsh = 512)
      parameter (mxgs = 2048)
      parameter (mxnpts = 20)
      parameter (mxneq = 255)
      parameter (mcwrkmx = 10)
c
c ==================================================================
c
c     parameter setting new drf routines (drf/dimpar):
c
c===================================================================
c
      parameter (mxst = 5)
      parameter (mxgran = 10)
cmw      parameter (mxnpts = 50)
      parameter (mxgrpar = (mxgran*(mxgran+1))/2)
c
c     ==================================================================
c
c     parameter setting for v. duijnen - gamess interface:
c
c=======================================================================
c
      parameter (ncatom=500)
c
c  ncatom is maximal allowed classical atoms in zmat
c  total allowed atoms in zmat = maxnz + ncatom
c
c  maxrfa is maximal number of classical atoms in zmat + coordinates
c  total maximal atoms allowed is maxrfa + maxat
c
      parameter (maxrfa=mxpts)
c
c  maxamb is maximal number of ambigous atoms. it is set equal
c  to maximal z-cards possible in gamess input , i.e. 50
c
      parameter (maxamb=mxamb)
c
c=============================================
c for the pes-routines
c=============================================
c
      parameter (mxgr1=128, mxpes=100)
c
c=============================================
c for the drfexpc routines
c=============================================
c
      parameter (maxex=1000)
c
c=============================================
c for the drfsurf routine
c=============================================
c
      parameter (maxinv=10,maxnrc=10000)
c      parameter (maxnrc=10000)
c
c=============================================
c for the analrf routine
c=============================================
c
c
      dimension s(nchd),dx(nchd),dy(nchd),dz(nchd)
      dimension iexpc(nchd)
      dimension omega(nwtc,nwtc)
c
c
      integer ir, iwr, ipnch,ipadiofil
      common /iofile/ ir,iwr,ipnch,ipadiofil(20)
c
      integer list
      common /hlistng/ list
c
c
      integer ia
      common /ijpair/ ia(3*mxpts)
c
c
c
      integer ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep
      integer neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga
      integer idipcal,iclinte,iclintd,ieffpol,isodis,iclintr
      integer iextdip,itsolv,itermax,isolsav,imomsav,idisadd
      integer ianal,ineqex,ithole,irevdis,ihbond,igrppol
      integer maxci_drf
      common/drfpar1/ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep,
     +               neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga,
     +               idipcal,iclinte,iclintd,ieffpol,isodis,iclintr,
     +               iextdip,itsolv,itermax,isolsav,imomsav,idisadd,
     +               ianal,ineqex,ithole,irevdis,ihbond,igrppol,
     +               maxci_drf
c
      integer nxtpts, npol, nexp, natint, nshint, namb, nspec, ngran
      common/drfpar2/nxtpts,npol,nexp,natint,nshint,namb,nspec,ngran
c
      integer ndima, ndimb, nbem, ndim, nexp4, npol3, nwtr, nwtc, nzfa
      integer lomga, nomga, nchd, ngrnam, ngrnam2, nzfp, nzfn
      integer nneq, nneqrf, maxneq, neqdim, nqdim, nqcls
      common/drfpar3/ndima,ndimb,nbem,ndim,nexp4,npol3,nwtr,nwtc,nzfa,
     +               lomga,nomga,nchd,ngrnam,ngrnam2,nzfp,nzfn,
     +               nneq,nneqrf,maxneq,neqdim,nqdim,nqcls
c
      integer ifldin, ifldout, idrfout, modxza, irepopt
      integer iexpza, icmexp, iadexp, nodpe, igetden, iarfcal, iqmclr
      common/drfpar4/ifldin,ifldout,idrfout,modxza,irepopt,
     +               iexpza,icmexp,iadexp,nodpe,igetden,iarfcal,
     +               iqmclr
c
      real*8 gamdrf, dstmin, dstmax, hbondl, hbondr, afact, rfact
      real*8 cvgrel, agrpe, agrpm, agrpc, scffact, acur
      real*8 conci_drf
      common/drfpar5/gamdrf,dstmin,dstmax,hbondl,hbondr,afact,rfact,
     +               cvgrel,agrpe,agrpm,agrpc,scffact,acur,conci_drf
c
      integer iind, isur1, isur2, ilwt, ilvr, illur, ilindx
      common/drfpar6/iind,isur1,isur2,ilwt,ilvr,illur,ilindx
c
      real*8 p, q, omgab
      common /drfexp/ p(4),q(4),omgab(4,4)
c
c
      real*8 val
      integer i, j, k, l
      common /drfint/ i,j,k,l,val
c
c
      dimension sump(4),sumq(4)
c
      vall = 0.0d0
c
      ij = ia(max(i,j)) + min(i,j)
      ijexp = iexpc(ij)
      p(1) = dx(ij)
      p(2) = dy(ij)
      p(3) = dz(ij)
      p(4) = s(ij)
c
      kl = ia(max(k,l)) + min(k,l)
      klexp = iexpc(kl)
      q(1) = dx(kl)
      q(2) = dy(kl)
      q(3) = dz(kl)
      q(4) = s(kl)
c
      call drfoab(ijexp,klexp,nwtc,omega)
      call matvec(omgab,q,sumq,4,.false.)
c
      call drfoab(klexp,ijexp,nwtc,omega)
      call matvec(omgab,p,sump,4,.false.)
c
c-----  note: + sign because source and recipient are electrons
c
      do 50, ii= 1, 4
        vall = vall + sump(ii)*q(ii) + sumq(ii)*p(ii)
   50 continue
c
c-----  test output
c
      if (idrfout .ge. 4) write(iwr,*) i, j, k, l, vall
      if (idrfout .ge. 5) then
        call hatout(p,1,4,2,'p')
        call hatout(q,1,4,2,'q')
        call hatout(omgab,4,4,2,'omgab')
      endif
c
      rftwog = vall
c
      return
      end
      subroutine pkfilrf(ii,jj,kk,ll,oskpa,oskpb,oskpc,onpsym,
     +                   g,q)
      implicit real*8  (a-h,p-w),integer  (i-n)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y),logical (o)
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
cafc
c
      character*8 field
      character*24 version
      common /drfopt/ field,version
      logical oreact
      integer intdrf
      common /drfopti/ oreact,intdrf
c
c
      integer ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep
      integer neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga
      integer idipcal,iclinte,iclintd,ieffpol,isodis,iclintr
      integer iextdip,itsolv,itermax,isolsav,imomsav,idisadd
      integer ianal,ineqex,ithole,irevdis,ihbond,igrppol
      integer maxci_drf
      common/drfpar1/ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep,
     +               neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga,
     +               idipcal,iclinte,iclintd,ieffpol,isodis,iclintr,
     +               iextdip,itsolv,itermax,isolsav,imomsav,idisadd,
     +               ianal,ineqex,ithole,irevdis,ihbond,igrppol,
     +               maxci_drf
c
      integer nxtpts, npol, nexp, natint, nshint, namb, nspec, ngran
      common/drfpar2/nxtpts,npol,nexp,natint,nshint,namb,nspec,ngran
c
      integer ndima, ndimb, nbem, ndim, nexp4, npol3, nwtr, nwtc, nzfa
      integer lomga, nomga, nchd, ngrnam, ngrnam2, nzfp, nzfn
      integer nneq, nneqrf, maxneq, neqdim, nqdim, nqcls
      common/drfpar3/ndima,ndimb,nbem,ndim,nexp4,npol3,nwtr,nwtc,nzfa,
     +               lomga,nomga,nchd,ngrnam,ngrnam2,nzfp,nzfn,
     +               nneq,nneqrf,maxneq,neqdim,nqdim,nqcls
c
      integer ifldin, ifldout, idrfout, modxza, irepopt
      integer iexpza, icmexp, iadexp, nodpe, igetden, iarfcal, iqmclr
      common/drfpar4/ifldin,ifldout,idrfout,modxza,irepopt,
     +               iexpza,icmexp,iadexp,nodpe,igetden,iarfcal,
     +               iqmclr
c
      real*8 gamdrf, dstmin, dstmax, hbondl, hbondr, afact, rfact
      real*8 cvgrel, agrpe, agrpm, agrpc, scffact, acur
      real*8 conci_drf
      common/drfpar5/gamdrf,dstmin,dstmax,hbondl,hbondr,afact,rfact,
     +               cvgrel,agrpe,agrpm,agrpc,scffact,acur,conci_drf
c
      integer iind, isur1, isur2, ilwt, ilvr, illur, ilindx
      common/drfpar6/iind,isur1,isur2,ilwt,ilvr,illur,ilindx
c
      integer i10, i20, i30, i40, i50, i60, i70, i80, i90 
      integer i100, i110, i120, i130, i140, i150, i160
      integer i170, i180
      common /scmprt/ i10,i20,i30,i40,i50,i60,i70,i80,i90,
     +          i100,i110,i120,i130,i140,i150,i160,i170,i180
c
      common/drfint/ii1,ii2,ii3,ii4,valext
c
      integer itwoeps,ioponly,itwosur,ixnsurf,ibemout,iuniout
      common /drfbem1/ itwoeps,ioponly,itwosur,ixnsurf,ibemout,iuniout
c
      integer leveli, levelo, nbem1, nbem2, nraw
      common /drfbem2/ leveli,levelo,nbem1,nbem2,nraw
c
      character*32 solnam
      real*8 radmax, eps1, eps2, spherad, rprobe, rprobej
      real*8 solrad, swidth, sdist
      real*8 kappa
      real*8 kappa1, kappa2, kappas
      common /drfbem3/ radmax,eps1,eps2,kappa1,kappa2,spherad,rprobe,
     +                 rprobej,solnam,solrad,swidth,sdist
c
      real*8 radat, cgrav
      common /drfbem4/ radat(maxat),cgrav(3)
c
cafc
c ***
      dimension q(*),g(*)
      logical osortp,oschw
      common/sortpdat/
     + igbf,iklbf,iiptbf,ijbase,ngbf,lengbf,osortp,oschw
c ***
      common/craypk/integ(680)
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
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
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
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
c
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
c
      logical oianj, okanl, oident, omisc, oham, opdipd, omp2
      integer ipos1, ipos2
      common /misc/ oianj,okanl,oident,omisc,
     +            oham,opdipd,omp2,ipos1,ipos2
c
c
      integer ik, ijx, ijy, ijz, klx, kly, klz
      integer ijgt, klgt
      common /indez/ ijgt(225),ijx(225),ijy(225),ijz(225),
     +       ik(225),klgt(225),klx(225),kly(225),klz(225)
c
      common/blkin/gout(510),nword
c     logical opk,opank,oskips
c     common/pkfil/opk,opank,oskips
      logical opk,opank
      common/pkfil/opk,opank
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      integer igt, jgt, kgt, lgt
      integer ibt, jbt, kbt, lbt, lenddd
      common /flips/ igt(3),jgt(3),kgt(3),lgt(3),
     +               ibt(mxprms),jbt(mxprms),kbt(mxprms),lbt(mxprms),
     +               lenddd
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
cafc
      logical odrf
cafc
      data dzero,pt5 /0.0d0,0.5d0/
cafc
      odrf = (field(5:) .eq. 'scf ') .and. (intdrf .ne. 0)
cafc
      lit = ktype(ii)
      mini = kmin(ii)
      maxi = kmax(ii)
      loci = kloc(ii)-mini
      minj = kmin(jj)
      maxj = kmax(jj)
      locj = kloc(jj)-minj
      lkt = ktype(kk)
      mink = kmin(kk)
      maxk = kmax(kk)
      lock = kloc(kk)-mink
      minl = kmin(ll)
      maxl = kmax(ll)
      locl = kloc(ll)-minl
      oianj = ii .eq. jj
      okanl = kk .eq. ll
      osame = (ii .eq. kk) .and. (jj .eq. ll)
c
c     type = 1 for (ii ii ii ii)
c            2     (ii jj jj jj)
c            3     (ii ii kk kk)
c            4     (ii ii ii ll)
c            5     (ii jj kk kk)
c            6     (ii jj jj ll)
c            7     (ii ii kk ll)
c            8     (ii jj kk ll)
c
      if (ii-jj) 100, 100, 240
  100 if (jj-kk) 120, 120, 180
  120 if (kk-ll) 140, 140, 160
  140 ntyp = 1
      go to 380
  160 ntyp = 4
      go to 380
  180 if (kk-ll) 200, 200, 220
  200 ntyp = 3
      go to 380
  220 ntyp = 7
      go to 380
  240 if (jj-kk) 260, 260,320
  260 if (kk-ll) 280,280,300
  280 ntyp = 2
      go to 380
  300 ntyp = 6
      go to 380
  320 if (kk-ll) 340,340,360
  340 ntyp = 5
      go to 380
  360 ntyp = 8
  380 continue
c
c     -----
c
      ind = 1
      if (oskpa .and. onpsym) ind = 2
      if (oskpb .and. onpsym) ind = 3
      if (oskpc .and. onpsym) ind = 3
      if (oskpa .and. oskpb .and. onpsym) ind = 4
      norg1 = 1
      norg2 = 50626
      norg3 = 101251
      if (oskpa .and. .not. onpsym) norg2 = 1
      if (oskpc .and. .not. onpsym) norg3 = 50626
      if (oskpb .and. .not. onpsym) norg3 = 1
      jmax = maxj
      kmaxx = maxk
      lmax = maxl
      do 1060 i = mini,maxi
      if (oianj) jmax = i
      do 1040 j = minj,jmax
      if (jj .eq. kk) kmaxx = j
      do 1020 k = mink,kmaxx
      if (okanl) lmax = k
      do 1000 l = minl,lmax
      ia = i-mini+1
      ja = j-minj+1
      ka = k-mink+1
      la = l-minl+1
c
c     ----- calculate n1 for (i,j//k,l) -----
c
      n1 = ibt(ia)+jbt(ja)+kbt(ka)+lbt(la)+norg1
c
c     ----- calculate n2 for (i,k//j,l) -----
c
      n2 = ibt(ia)+jbt(ka)+kbt(ja)+lbt(la)+norg2
c
c     ----- calculate n3 for (i,l//j,k) -----
c
      go to (400,420,400,440,420,420,440,420),ntyp
  400 if (ia .eq. ja) go to 440
  420 n3 = ibt(ia)+jbt(la)+kbt(ja)+lbt(ka)+norg3
      go to 460
  440 n3 = ibt(ja)+jbt(ka)+kbt(ia)+lbt(la)+norg3
  460 continue
c
c     ----- form first linear combination -----
c
cxxx
cafc      jump = 1
      jump = 1
      i1 = loci+i
      i2 = locj+j
      i3 = lock+k
      i4 = locl+l
c     i1 = ii
c     i2 = jj    
c     i3 = kk    
c     i4 = ll     
c
c     ------------ drf extension
c
c      --------  p.th. van duijnen, groningen 1991 -----
c
      if (odrf) then
        ii1 = loci+i
        ii2 = locj+j
        ii3 = lock+k
        ii4 = locl+l
c 1-----
c  -----  total dielectric contribution
c
        valext = dzero
        call drftwo(q(i110),q(i120),q(i130),q(i140),
     1               q(i150),q(i160))
        if (abs(valext) .lt. cutoff) valext = dzero
        gdrf1 = scffact*valext
        gdrfo1 = gdrf1
c
        ii3 = locj + j
        ii2 = lock + k
        valext = dzero
        call drftwo(q(i110),q(i120),q(i130),q(i140),
     1               q(i150),q(i160))
        if (abs(valext) .lt. cutoff) valext = dzero
        gdrf2 = scffact*valext
        gdrfo2 = gdrf2
c
        ii4 = locj + j
        ii3 = lock + k
        ii2 = locl + l
        valext = dzero
        call drftwo(q(i110),q(i120),q(i130),q(i140),
     1               q(i150),q(i160))
        if (abs(valext) .lt. cutoff) valext = dzero
        gdrf3 = scffact*valext
        gdrfo3 = gdrf3
c
c  -----  optic dielectric contribution
c
        if (itwoeps .ne. 0) then
c   2-----
          ii2 = locj + j
          ii3 = lock + k
          ii4 = locl + l
          valext = dzero
          call drftwo(q(i110),q(i120),q(i130),q(i140),
     1               q(i150),q(i170))
          if (abs(valext) .lt. cutoff) valext = dzero
          gdrfo1 = scffact*valext
c
          ii3 = locj + j
          ii2 = lock + k
          valext = dzero
          call drftwo(q(i110),q(i120),q(i130),q(i140),
     1               q(i150),q(i170))
          if (abs(valext) .lt. cutoff) valext = dzero
          gdrfo2 = scffact*valext
c
          ii4 = locj + j
          ii3 = lock + k
          ii2 = locl + l
          valext = dzero
          call drftwo(q(i110),q(i120),q(i130),q(i140),
     1               q(i150),q(i170))
          if (abs(valext) .lt. cutoff) valext = dzero
          gdrfo3 = scffact*valext
c   2-----
        endif
c 1-----
      endif
c
c-----  end drf extension
c
      jump = 1
      if (i2 .eq. i3) jump = 2
      if ((i2 .eq. i4) .or. (i1 .eq. i3)) jump = 3
      g1 = g(n1)
      g2 = g(n2)
      g3 = g(n3)
      go to (480,500,520,540),ind
  480 valk = g2+g3
      if (odrf) valk = valk + (gdrfo2+gdrfo3)*gamdrf
      valp = (g1+g1)+(g1+g1)-valk
      if (odrf) valp = valp + (gdrf1+gdrf1) +(gdrf1+gdrf1)
      go to 560
  500 valk = g3
      valp = (g1+g1)+(g1+g1)-valk
      go to 560
  520 valk = g2
      valp = (g1+g1)+(g1+g1)-valk
      go to 560
  540 valk = dzero
      valp = (g1+g1)+(g1+g1)
  560 continue
      g(n1) = valp
      nn = n1
      go to 820
c
c     ----- form second linear combination -----
c
  580 go to (600,620,640,660),ind
  600 valk = g3+g1
      if (odrf) valk = valk + (gdrfo3+gdrfo1)*gamdrf
      valp = (g2+g2)+(g2+g2)-valk
      if (odrf) valp = valp + (gdrf2+gdrf2) +(gdrf2+gdrf2)
      go to 680
  620 valk = g1+g3
      valp = -valk
      go to 680
  640 valk = g1
      valp = (g2+g2)+(g2+g2)-valk
      go to 680
  660 valk = g1
      valp = -valk
  680 continue
      g(n2) = valp
      nn = n2
      jump = 2
      if ((i1 .eq. i2) .or. (i3 .eq. i4)) jump = 3
      n = i2
      i2 = i3
      i3 = n
      go to 820
c
c     ----- form third linear combination -----
c
  700 go to (720,740,760,780),ind
  720 valk = g1+g2
      if (odrf) valk = valk + (gdrfo1+gdrfo2)*gamdrf
      valp = (g3+g3)+(g3+g3)-valk
      if (odrf) valp = valp + (gdrf3+gdrf3) +(gdrf3+gdrf3)
      go to 800
  740 valk = g1
      valp = (g3+g3)+(g3+g3)-valk
      go to 800
  760 valk = g1+g2
      valp = -valk
      go to 800
  780 valk = g1
      valp = -valk
  800 continue
      g(n3) = valp
      nn = n3
      i2 = locl+l
      i3 = locj+j
      i4 = lock+k
      jump = 3
  820 continue
c
c     ----- store integral and indices -----
c
      if (opank) go to 880
c
c     ----- -p- supermatrix only. -----
c
      if ( dabs(valp) .lt. cutoff) go to 980
      valp0 = valp
      if (i1 .eq. i3 .and. i2 .eq. i4) valp = valp*pt5
      integ(ic4  )=iky(i1)+i2
      integ(ic4+1)=iky(i3)+i4
      ic4=ic4+2
      gout(icount) = valp
      icount = icount+1
      if (icount .le. nintmx) go to 980
c ***
      if (osortp) then
        call zsortp(q(igbf),q(iklbf),q(iiptbf),ijbase)
      else
        call blocki
      endif
c ***
      if(omaxb)go to 1061
      go to 980
  880 continue
c
c     ----- -p- and -k- supermatrices -----
c
      if (dabs(valp).lt.cutoff .and. dabs(valk).lt.cutoff) go to 980
      valp0 = valp
      valk0 = valk
      if (i1 .ne. i3 .or. i2 .ne. i4) go to 900
      valp = valp*pt5
      valk = valk*pt5
  900 continue
      integ(ic4  )=iky(i1)+i2
      integ(ic4+1)=iky(i3)+i4
      ic4 = ic4 + 2
      gout(icount) = valp
      gout(num2ejk+icount) = valk
      icount = icount+1
      if (icount .le. nintmx) go to 980
      call blocki
      if(omaxb)go to 1061
  980 go to (580,700,1000),jump
 1000 continue
 1020 continue
 1040 continue
 1060 continue
 1061 return
      end
      subroutine standv(ncall,core)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
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
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
c
      integer nmaxly, ntotly, icurly, isecly
      logical oprintm
      common /segm/ nmaxly,ntotly,icurly,isecly(100),
     +              oprintm
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
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
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
      logical ciopt, ciforc, mp2, hfgr, bfgs, ump2, lmeth2
      logical ump3, rmp3, ordmo, mp2w, loptor, ladp, lcpf
      logical lopti, lmcscf, lforce, lci, lcart, lmcdat
      logical lfdtrn, unit7, lcontr, lvcd, lgten
      logical ldenom, ignore, ldens, lset, ladapt, lsym, latmol
      logical berny, llibry, limpt, fpres, oss, ldiag, lskip
      logical opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common /restrl/ ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,
     +rmp3,ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
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
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
c
c     ZORA common
c
      logical ozora,oscalz,onlyp,onlyd,onlyf,onlyg,opzora,
     1        op2zora,onlys,ocontr,oint_zora,o1e_zora,osmall_zora,
     2        opre_zora,oso,onoext_z,oioraz,oatint_z,oatscf_z,
     3        oscalatz
      real*8 cspeed,critat_z
      integer nbas_zora,ibls_z,iblv_z,nwv_z,ibiso_z,nwiso_z,ibcor_z,
     1        ibsx_z,icoul_z,ibsy_z,ibsz_z,ibsmin_z,ibsplus_z,nwmat_z,
     2        ibtrin_z,ibtrout_z,ibscalf_z,nwtrat_z,ibscale_z,ibcoul_z,
     3        niter_z,int_zora,num_ext,igauge_z,nwcor_z,isexch_z,
     4        nat_z,ibdens_z,nwdens_z,ibscalao_z,ibshat_z,is_z,
     5        irest_z,numdu_z,ibldu_z,iblock_z,lshift_z,ibl7ew_z
      common/zorac/ cspeed,critat_z,ozora,oscalz,
     1              onlyp,onlyd,onlyf,onlyg,opzora,
     1              op2zora,onlys,ocontr,nbas_zora,
     2              oint_zora,o1e_zora,icoul_z,
     3              ibls_z,iblv_z,nwv_z,
     3              ibiso_z,nwiso_z,ibcor_z,nwcor_z,
     4              ibsx_z,ibsy_z,ibsz_z,ibsmin_z,ibsplus_z,nwmat_z,
     5              ibtrin_z,ibtrout_z,ibscalf_z,nwtrat_z,
     6              ibscale_z,ibcoul_z,ibdens_z,nwdens_z,osmall_zora,
     7              niter_z,opre_zora,int_zora,num_ext,isexch_z,
     8              oso,onoext_z,is_z,igauge_z,oioraz,
     9              nat_z,ibscalao_z,oatint_z,oatscf_z,ibshat_z,
     1              oscalatz,irest_z,numdu_z,ibldu_z,iblock_z,lshift_z,
     2              ibl7ew_z
c....    icoul_z : 1 : Full coulomb
c...               2 : atomic coulomb (default)
c...               3 : no coulomb (just 1-electron ints)
c...               4 : small basis coulomb
c...     isexch_z  0 : no exchange in zora correction (default)
c...               1 : exchange in zora correction (like in dft)
c...     niter_z : # iterations between coulomb matrix update
c...     int_zora :internal basis gen. : 0 none, 1, automatic, 2 by hand
c
c...     oscalz  : scaling on or off (true or false)
c...     oioraz  : change metric to perform inf. order scaling 
c...               (Dyall/EvLenthe (not quite))
c...     oscalatz: atoms only scaled zora (together with get atom)
c...     only'x' : internal basis only op to 'x'-type functions
c...     op(2)zora: print flags (op2zora+opzora: prints *everything*)
c...     ocontr  : adds l+1 block to internal basis in contracted way 
c...     nbas_zora: size of internal basis
c...     oint_zora: in internal or external mode???
c...     o1e_zora: one electron integrals have been calculated in *this*
c...               internal basis
c...     osmall_zora: small swap for projected coulomb
c...     opre_zora: pre_zora has been called
c...     oso     : spin orbit calculation 
c...     onoext_z: no external functions copied to internal basis
c...     is_z : section for orbitals used for  coulomb matrix
c...     igauge_z : add so many s-marices to h-matrix
c...                the energy should change by igauge
c...     nat_z  : indicates (if ne 0) that starting density is used
c...              to calc. coulomb; gives # it. to get atomic density
c...     critat_z : accuracy of atomic energy for atomic zora
c...     oatint_z : indicates that integrals to be calculated are atomic
c...     oatscf_z : indicates we are doing a set of atomic SCFs
c...     ibshat_z : block on num8 to store blocked s/h
c...     for disk-adress and lengths (ib.., nw.. ) see subroutine atoms2
c...     irest_z  : .ne.0 indicates atomic zora corrections from dumpfile; 1 not on restart
c...     numdu_z,ibldu_z : if ne 0 other dumpfile to read zora corrections from
c...     iblock_z is the block number on ed3, where the zora corrections (section 425) reside
c...     ibl7ew is where the energy weighted density matrix resides
c
c
      character*8 field
      character*24 version
      common /drfopt/ field,version
      logical oreact
      integer intdrf
      common /drfopti/ oreact,intdrf
c
c crystal field
      integer isecfx 
      logical ocryst
      common/xfield/isecfx,ocryst
c---------
c  omtchg - omit 1e- integrals involving the specified charges
c
c   integrals of the form  < bf(a) | q(b) | bf (a) > are deleted
c
c   where a and b are centres indicated as follows
c
c  ztagcc - list of centre names a
c  ztagcq - for each a in ztagcc, the list of centres b
c
c  omtslf - omit self-energy of point charge array (applied to all
c           atoms with centre label of form bq* )
c      CHECK
c      parameter (mxombq=20,mxomat=50)
c---------

c
c obqf  - special treatment of forces .. skip force contributions
c         on non-bq centres
c         for use in qm/mm optimisations of surrounding
c         assumption is that bq's in this case don't have basis
c         functions or ecps
c
c omitted charge specifications
c

      integer mxombq, mxomat
      parameter (mxombq=50,mxomat=20)
      character*8 ztagcc, ztagcq
      logical omtchg, omtslf, obqf
      common /chgcc/ ztagcc(mxomat),ztagcq(mxomat,mxombq),omtchg,
     &     omtslf,obqf
c
      real*8 tchg
      integer ifree, kform
      logical opress,otemp,omapp,omodel,ocharg,oaminc,oamhes,osubs
      common /modj/ kform,opress,otemp,omapp,omodel,ocharg,
     +              oaminc,oamhes,tchg(maxat),osubs,ifree
c
      logical mpassi
      common/integmp/mpassi
c
c
      common/junk3/ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720)
c
      character*8 fnm
      character*6 snm
      data fnm,snm/"intega.m","standv"/
c
      dimension core(*) 
c
      if (oint_zora) then
         i10 = igmem_alloc_inf(nwiso_z,fnm,snm,"i10",IGMEM_DEBUG)
      else    
         i10 = igmem_alloc_inf(nw196(5),fnm,snm,"i10",IGMEM_DEBUG)
      end if  
c
c----- read in transformation matrices for s,p,d,f and g basis functions
c
      call rdedx(ptr,nw196(1),ibl196(1),idaf)
      if (odbas) call rdedx(dtr,nw196(2),ibl196(2),idaf)
      if (ofbas) call rdedx(ftr,nw196(3),ibl196(3),idaf)
      if (ogbas) call rdedx(gtr,nw196(4),ibl196(4),idaf)
c
c----- read in symmetry array - iso
c
      if (oint_zora) then
         call rdedx(core(i10),nwiso_z,ibiso_z,num8)
      else    
         call rdedx(core(i10),nw196(5),ibl196(5),idaf)
      end if  
c
c     now decide if this is a full triangle, single pass call, 
c     or one where only 2 triangles are to to be used in each call 
c     to stvstv and secdip
c
      if (ocryst.or.(lpseud.eq.1.or.lpseud.eq.2).or.oreact.
     + or.oaminc)then
       call stvstv(ncall,core,core(i10),nshell)
      else if (mpassi) then
c     multipass one electron integrals
       o1e = .false.
       do ipass = 1,3
       call stvstv_mp(ncall,core,core(i10),nshell,ipass,o1e)
c      write(iwr,*)' IPASS, O1E = ', ipass,o1e
       if (o1e) go to 100
       enddo
      else
       call stvstv(ncall,core,core(i10),nshell)
      endif
c
100   call gmem_free_inf(i10,fnm,snm,"i10")
      return  
      end 
      subroutine stvstv_mp(ncall,q,iso,nshels,ipass,o1e)
c
c    this routine is a modified version of stvstv designed
c    to minimise the use of storage. It now computes S, T and T+V
c    in three passes (ipass), thereby using just 2 triangles
c    Note that the this requires a modified version of sym1e
c    to allocate the 2nd triangle dynamically given the workings
c    of pg_dgop that also allocates a triangle dynamically
c
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
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
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
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
      integer nmaxly, ntotly, icurly, isecly
      logical oprintm
      common /segm/ nmaxly,ntotly,icurly,isecly(100),
     +              oprintm
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
c
      integer ibl7s, ibl7t, ibl7f, ibl7st,ibl7tt,ibl7ft,ibl3qs
      integer ibl7qa,ibl7pa,ibl7ea,ibl7qb,ibl7pb,ibl7eb
      integer ibl7x, ibl7y, ibl7z
      integer ibl7so,ibl7o, ibl7la,ibl7ec,ibl7ha
      common /scra7/ ibl7s, ibl7t, ibl7f, ibl7st,ibl7tt,ibl7ft,ibl3qs,
     +               ibl7qa,ibl7pa,ibl7ea,ibl7qb,ibl7pb,ibl7eb,
     +               ibl7x, ibl7y, ibl7z,
     +               ibl7so,ibl7o(8), ibl7la, ibl7ec, ibl7ha
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
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
      common/junk/s(225),g(225),
     *pint,qint,rint,t,p0,q0,r0,pi,qi,ri,pj,qj,rj,ni,nj,
     *tol,ii,jj,lit,ljt,mini,minj,maxi,maxj,iandk
c
      real*8 pp, u, w, dji
      integer nroots, nroot2
      common /root/ pp,u(12),w(12),nroots,nroot2,dji(60)
c
      common/blkin/dxyz(4),gg(225),ft(225),fx(225),dij(225),
     + pin(125),qin(125),rin(125),
     + ijx(225),ijy(225),ijz(225)
c mechanics
c
      real*8 tchg
      integer ifree, kform
      logical opress,otemp,omapp,omodel,ocharg,oaminc,oamhes,osubs
      common /modj/ kform,opress,otemp,omapp,omodel,ocharg,
     +              oaminc,oamhes,tchg(maxat),osubs,ifree
c
      common/g80nb/natmod,nnbg,mapmod(maxat),ipnonb(maxat+1),
     +             inonb(maxat*6)
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
c ***** omit specified charges from attraction terms ***
c---------
c  omtchg - omit 1e- integrals involving the specified charges
c
c   integrals of the form  < bf(a) | q(b) | bf (a) > are deleted
c
c   where a and b are centres indicated as follows
c
c  ztagcc - list of centre names a
c  ztagcq - for each a in ztagcc, the list of centres b
c
c  omtslf - omit self-energy of point charge array (applied to all
c           atoms with centre label of form bq* )
c      CHECK
c      parameter (mxombq=20,mxomat=50)
c---------

c
c obqf  - special treatment of forces .. skip force contributions
c         on non-bq centres
c         for use in qm/mm optimisations of surrounding
c         assumption is that bq's in this case don't have basis
c         functions or ecps
c
c omitted charge specifications
c

      integer mxombq, mxomat
      parameter (mxombq=50,mxomat=20)
      character*8 ztagcc, ztagcq
      logical omtchg, omtslf, obqf
      common /chgcc/ ztagcc(mxomat),ztagcq(mxomat,mxombq),omtchg,
     &     omtslf,obqf
c crystal field
      integer isecfx 
      logical ocryst
      common/xfield/isecfx,ocryst
c
c     ZORA common
c
      logical ozora,oscalz,onlyp,onlyd,onlyf,onlyg,opzora,
     1        op2zora,onlys,ocontr,oint_zora,o1e_zora,osmall_zora,
     2        opre_zora,oso,onoext_z,oioraz,oatint_z,oatscf_z,
     3        oscalatz
      real*8 cspeed,critat_z
      integer nbas_zora,ibls_z,iblv_z,nwv_z,ibiso_z,nwiso_z,ibcor_z,
     1        ibsx_z,icoul_z,ibsy_z,ibsz_z,ibsmin_z,ibsplus_z,nwmat_z,
     2        ibtrin_z,ibtrout_z,ibscalf_z,nwtrat_z,ibscale_z,ibcoul_z,
     3        niter_z,int_zora,num_ext,igauge_z,nwcor_z,isexch_z,
     4        nat_z,ibdens_z,nwdens_z,ibscalao_z,ibshat_z,is_z,
     5        irest_z,numdu_z,ibldu_z,iblock_z,lshift_z,ibl7ew_z
      common/zorac/ cspeed,critat_z,ozora,oscalz,
     1              onlyp,onlyd,onlyf,onlyg,opzora,
     1              op2zora,onlys,ocontr,nbas_zora,
     2              oint_zora,o1e_zora,icoul_z,
     3              ibls_z,iblv_z,nwv_z,
     3              ibiso_z,nwiso_z,ibcor_z,nwcor_z,
     4              ibsx_z,ibsy_z,ibsz_z,ibsmin_z,ibsplus_z,nwmat_z,
     5              ibtrin_z,ibtrout_z,ibscalf_z,nwtrat_z,
     6              ibscale_z,ibcoul_z,ibdens_z,nwdens_z,osmall_zora,
     7              niter_z,opre_zora,int_zora,num_ext,isexch_z,
     8              oso,onoext_z,is_z,igauge_z,oioraz,
     9              nat_z,ibscalao_z,oatint_z,oatscf_z,ibshat_z,
     1              oscalatz,irest_z,numdu_z,ibldu_z,iblock_z,lshift_z,
     2              ibl7ew_z
c....    icoul_z : 1 : Full coulomb
c...               2 : atomic coulomb (default)
c...               3 : no coulomb (just 1-electron ints)
c...               4 : small basis coulomb
c...     isexch_z  0 : no exchange in zora correction (default)
c...               1 : exchange in zora correction (like in dft)
c...     niter_z : # iterations between coulomb matrix update
c...     int_zora :internal basis gen. : 0 none, 1, automatic, 2 by hand
c
c...     oscalz  : scaling on or off (true or false)
c...     oioraz  : change metric to perform inf. order scaling 
c...               (Dyall/EvLenthe (not quite))
c...     oscalatz: atoms only scaled zora (together with get atom)
c...     only'x' : internal basis only op to 'x'-type functions
c...     op(2)zora: print flags (op2zora+opzora: prints *everything*)
c...     ocontr  : adds l+1 block to internal basis in contracted way 
c...     nbas_zora: size of internal basis
c...     oint_zora: in internal or external mode???
c...     o1e_zora: one electron integrals have been calculated in *this*
c...               internal basis
c...     osmall_zora: small swap for projected coulomb
c...     opre_zora: pre_zora has been called
c...     oso     : spin orbit calculation 
c...     onoext_z: no external functions copied to internal basis
c...     is_z : section for orbitals used for  coulomb matrix
c...     igauge_z : add so many s-marices to h-matrix
c...                the energy should change by igauge
c...     nat_z  : indicates (if ne 0) that starting density is used
c...              to calc. coulomb; gives # it. to get atomic density
c...     critat_z : accuracy of atomic energy for atomic zora
c...     oatint_z : indicates that integrals to be calculated are atomic
c...     oatscf_z : indicates we are doing a set of atomic SCFs
c...     ibshat_z : block on num8 to store blocked s/h
c...     for disk-adress and lengths (ib.., nw.. ) see subroutine atoms2
c...     irest_z  : .ne.0 indicates atomic zora corrections from dumpfile; 1 not on restart
c...     numdu_z,ibldu_z : if ne 0 other dumpfile to read zora corrections from
c...     iblock_z is the block number on ed3, where the zora corrections (section 425) reside
c...     ibl7ew is where the energy weighted density matrix resides
      character *8 title,guess
      common/restrz/title(12),guess
c
      character*8 field
      character*24 version
      common /drfopt/ field,version
      logical oreact
      integer intdrf
      common /drfopti/ oreact,intdrf
c
      dimension q(*),iso(nshels,*)
      dimension ix(35),iy(35),iz(35),jx(35),jy(35),jz(35)
      dimension m0(48)
c mechanics
      data dzero,pt5,done,two,three,five,seven,
     +     rnine,eleven  /0.0d0,0.5d0,1.0d0,
     + 2.0d0,3.0d0,5.0d0,7.0d0,9.0d0,11.0d0/
      data pi212 /1.1283791670955d0/
      data sqrt3 /1.73205080756888d0/
      data sqrt5 /2.23606797749979d0/
      data sqrt7 /2.64575131106459d0/
      data rln10 /2.30258d0/
c
      data jx / 0, 1, 0, 0, 2, 0, 0, 1, 1, 0,
     +          3, 0, 0, 2, 2, 1, 0, 1, 0, 1,
     +          4, 0, 0, 3, 3, 1, 0, 1, 0, 2,
     +          2, 0, 2, 1, 1/
      data ix / 1, 6, 1, 1,11, 1, 1, 6, 6, 1,
     *         16, 1, 1,11,11, 6, 1, 6, 1, 6,
     *         21, 1, 1,16,16, 6, 1, 6, 1,11,
     *         11, 1,11, 6, 6/
      data jy / 0, 0, 1, 0, 0, 2, 0, 1, 0, 1,
     +          0, 3, 0, 1, 0, 2, 2, 0, 1, 1,
     +          0, 4, 0, 1, 0, 3, 3, 0, 1, 2,
     +          0, 2, 1, 2, 1/
      data iy / 1, 1, 6, 1, 1,11, 1, 6, 1, 6,
     +          1,16, 1, 6, 1,11,11, 1, 6, 6,
     +          1,21, 1, 6, 1,16,16, 1, 6,11,
     +          1,11, 6,11, 6/
      data jz / 0, 0, 0, 1, 0, 0, 2, 0, 1, 1,
     +          0, 0, 3, 0, 1, 0, 1, 2, 2, 1,
     +          0, 0, 4, 0, 1, 0, 1, 3, 3, 0,
     +          2, 2, 1, 1, 2/
      data iz / 1, 1, 1, 6, 1, 1,11, 1, 6, 6,
     +          1, 1,16, 1, 6, 1, 6,11,11, 6,
     +          1, 1,21, 1, 6, 1, 6,16,16, 1,
     +         11,11, 6, 6,11/
      l1 = num
      ono  = .false.
      outv = oprint(59)
      nav = lenwrd()
      if (nprint.eq.-5) outv = .false.
c
      if (ipass.eq.1) then
c     have the 1e-integrals already been calculated ..
       call chek1e(o1e)
       if (o1e) return
      endif
c
c     ----- set pointers for partitioning of core -----
c
      l2 = (num*(num+1))/2
      length = l2
      if(ipass.eq.3) then
       if (omtchg) then
        length = length + nat * nat
       endif
      endif
      i10 = igmem_alloc(length)
      if(ipass.eq.3.and.omtchg) i110 = i10 + nat * nat
c
      if (ncall.eq.0 .and. nprint.eq.3) then
       if(ipass.eq.1) then
        last = i10 + length
        write (iwr,6010) i10 , last
       endif
      endif
c
      if (ipass.eq.3.and.omtchg) then
       call omitch(q(i110),nat,nona)
      endif
c
c     ----- calculate -s- and -h0- matrices -----
c
c     - s- at x(i10) in the 1st pass
c     - t- at x(i10) on the 2nd pass
c     -h0- at x(i10) on the 3rd pass
c
         cpu = cpulft(1)
c
         if (ncall.eq.0) call cpuwal(begin,ebegin)
         if (ncall.eq.0 .and. outv) write (iwr,6030) cpu
         tol = rln10*itol
         out = nprint.eq.3
         if (out) then
           do ii = 1,6
            oprn(30+ii) = .true.
           enddo
          else
           do ii = 1,6
            if(oprn(30+ii)) out = .true.
           enddo
         endif
         onorm = normf.ne.1 .or. normp.ne.1
         call vclr(q(i10),1,l2)
c
c     ----- ishell
c
         do 440 ii = 1 , nshell
            i = katom(ii)
c
c     ----- eliminate ishell -----
c
            do 450 it = 1 , nt
               id = iso(ii,it)
               if (id.gt.ii) go to 440
               m0(it) = id
450         continue
            icent = i
            pi = c(1,i)
            qi = c(2,i)
            ri = c(3,i)
            i1 = kstart(ii)
            i2 = i1 + kng(ii) - 1
            lit = ktype(ii)
            mini = kmin(ii)
            maxi = kmax(ii)
            loci = kloc(ii) - mini
c
c     ----- jshell
c
            do 430 jj = 1 , ii
               j = katom(jj)
               n2 = 0
               do 470 it = 1 , nt
                  jd = iso(jj,it)
                  if (jd.gt.ii) go to 430
                  id = m0(it)
                  if (id.lt.jd) then
                     nd = id
                     id = jd
                     jd = nd
                  end if
                  if (id.eq.ii .and. jd.gt.jj) go to 430
                  if (id.eq.ii.and.jd.eq.jj) then
                      n2 = n2 + 1
                  end if
470            continue
               q2 = dble(nt)/dble(n2)
               jcent = j
               pj = c(1,j)
               qj = c(2,j)
               rj = c(3,j)
               j1 = kstart(jj)
               j2 = j1 + kng(jj) - 1
               ljt = ktype(jj)
               minj = kmin(jj)
               maxj = kmax(jj)
               locj = kloc(jj) - minj
               nroots = (lit+ljt-2)/2 + 1
               rr = (pi-pj)**2 + (qi-qj)**2 + (ri-rj)**2
               oiandj = ii.eq.jj
c
c     ----- prepare indices for pairs of (i,j) functions
c
               ij = 0
               max = maxj
               do 30 i = mini , maxi
                  nnx = ix(i)
                  nny = iy(i)
                  nnz = iz(i)
                  if (oiandj) max = i
                  do 20 j = minj , max
                     ij = ij + 1
                     ijx(ij) = nnx + jx(j)
                     ijy(ij) = nny + jy(j)
                     ijz(ij) = nnz + jz(j)
                     if (j.le.1) then
                        ft(ij) = three
                     else if (j.le.4) then
                        ft(ij) = five
                     else if (j.le.10) then
                        ft(ij) = seven
                     else if (j.gt.20) then
                        ft(ij) = eleven
                     else
                        ft(ij) = rnine
                     end if
 20               continue
 30            continue
c
               if (ipass.eq.1) then
                do i = 1 , ij
                 s(i) = dzero
                enddo
               else if (ipass.eq.2) then
                do i = 1 , ij
                 gg(i) = dzero
                enddo
               else
                do i = 1 , ij
                 g(i) = dzero
                enddo
               endif
c
c     ----- i primitive
c
               jgmax = j2
               do 400 ig = i1 , i2
                  ai = ex(ig)
                  arri = ai*rr
                  axi = ai*pi
                  ayi = ai*qi
                  azi = ai*ri
                  csi = cs(ig)
                  cpi = cp(ig)
                  cdi = cd(ig)
                  cfi = cf(ig)
                  cgi = cg(ig)
c
c     ----- j primtive
c
                  if (oiandj) jgmax = ig
                  do 390 jg = j1 , jgmax
                     aj = ex(jg)
                     aa = ai + aj
                     aa1 = done/aa
                     dum = aj*arri*aa1
                     if (dum.le.tol) then
                        fac = dexp(-dum)
                        csj = cs(jg)
                        cpj = cp(jg)
                        cdj = cd(jg)
                        cfj = cf(jg)
                        cgj = cg(jg)
                        ax = (axi+aj*pj)*aa1
                        ay = (ayi+aj*qj)*aa1
                        az = (azi+aj*rj)*aa1
                        odoub = oiandj .and. ig.ne.jg
c
c     ----- density factor
c
                        max = maxj
                        nn = 0
                        do 220 i = mini , maxi
                           go to (50,60,120,120,
     +                            70,120,120,80,120,120,
     +                            90,120,120,100,120,120,120,120,120,
     +                            110,
     +                            112,120,120,114,120,120,120,120,120,
     +                            116,120,120,118,120,120), i
c
 50                        dum1 = csi*fac
                           go to 120
 60                        dum1 = cpi*fac
                           go to 120
 70                        dum1 = cdi*fac
                           go to 120
 80                        if (onorm) dum1 = dum1*sqrt3
                           go to 120
 90                        dum1 = cfi*fac
                           go to 120
 100                       if (onorm) dum1 = dum1*sqrt5
                           go to 120
 110                       if (onorm) dum1 = dum1*sqrt3
                           go to 120
 112                       dum1 = cgi*fac
                           go to 120
 114                       if (onorm) dum1 = dum1*sqrt7
                           go to 120
 116                       if (onorm) dum1 = dum1*sqrt5/sqrt3
                           go to 120
 118                       if (onorm) dum1 = dum1*sqrt3
 120                       if (oiandj) max = i
                           do 210 j = minj , max
                              go to (130,140,200,200,
     +                               150,200,200,160,200,200,
     +                               170,200,200,180,200,200,
     +                               200,200,200,190,
     +                               192,200,200,194,200,200,200,200,
     +                               200,196,200,200,198,200,200),j
 130                          dum2 = dum1*csj
                              if (odoub) then
                                 if (i.gt.1) then
                                    dum2 = dum2 + csi*cpj*fac
                                 else
                                    dum2 = dum2 + dum2
                                 end if
                              end if
                              go to 200
 140                          dum2 = dum1*cpj
                              if (odoub) dum2 = dum2 + dum2
                              go to 200
 150                          dum2 = dum1*cdj
                              if (odoub) dum2 = dum2 + dum2
                              go to 200
 160                          if (onorm) dum2 = dum2*sqrt3
                              go to 200
 170                          dum2 = dum1*cfj
                              if (odoub) dum2 = dum2 + dum2
                              go to 200
 180                          if (onorm) dum2 = dum2*sqrt5
                              go to 200
 190                          if (onorm) dum2 = dum2*sqrt3
                              go to 200
 192                          dum2 = dum1*cgj
                              if (odoub) dum2 = dum2 + dum2
                              go to 200
 194                          if (onorm) dum2 = dum2*sqrt7
                              go to 200
 196                          if (onorm) dum2 = dum2*sqrt5/sqrt3
                              go to 200
 198                          if (onorm) dum2 = dum2*sqrt3
 200                          nn = nn + 1
                              dij(nn) = dum2
 210                       continue
 220                    continue
c
c     ----- overlap and kinetic energy
c
                        t = dsqrt(aa1)
                        t1 = -two*aj*aj*t
                        t2 = -pt5*t
                        p0 = ax
                        q0 = ay
                        r0 = az
                        in = -5
                        do 240 i = 1 , lit
                           in = in + 5
                           ni = i
                           do 230 j = 1 , ljt
                              jn = in + j
                              nj = j
                              call stvint
                              pin(jn) = pint*t
                              qin(jn) = qint*t
                              rin(jn) = rint*t
                              nj = j + 2
                              call stvint
                              pin(jn+25) = pint*t1
                              qin(jn+25) = qint*t1
                              rin(jn+25) = rint*t1
                              nj = j - 2
                              if (nj.gt.0) then
                                 call stvint
                              else
                                 pint = dzero
                                 qint = dzero
                                 rint = dzero
                              end if
                              n = (j-1)*(j-2)
                              dum = dble(n)*t2
                              pin(jn+50) = pint*dum
                              qin(jn+50) = qint*dum
                              rin(jn+50) = rint*dum
 230                       continue
 240                    continue
                        if(ipass.eq.1) then
                         do i = 1 , ij
                          nnx = ijx(i)
                          nny = ijy(i)
                          nnz = ijz(i)
                          pyz = qin(nny)*rin(nnz)
                          dum = pyz*pin(nnx)
                          s(i) = s(i) + dij(i)*dum
                         enddo
                        else if(ipass.eq.2) then
                         do i = 1 , ij
                          nnx = ijx(i)
                          nny = ijy(i)
                          nnz = ijz(i)
                          pyz = qin(nny)*rin(nnz)
                          dum = pyz*pin(nnx)
                          dum1 = (pin(nnx+25)+pin(nnx+50))*pyz + 
     +                    (qin(nny+25)+qin(nny+50))*pin(nnx)*rin(nnz) + 
     +                    (rin(nnz+25)+rin(nnz+50))*pin(nnx)*qin(nny)
                          gg(i) = gg(i) + dij(i)*(dum*aj*ft(i)+dum1)
                         enddo
                        else
                         do i = 1 , ij
                         nnx = ijx(i)
                         nny = ijy(i)
                         nnz = ijz(i)
                         pyz = qin(nny)*rin(nnz)
                         dum = pyz*pin(nnx)
                         dum1 = (pin(nnx+25)+pin(nnx+50))*pyz + 
     +                   (qin(nny+25)+qin(nny+50))*pin(nnx)*rin(nnz) + 
     +                   (rin(nnz+25)+rin(nnz+50))*pin(nnx)*qin(nny)
                         g(i) = g(i) + dij(i)*(dum*aj*ft(i)+dum1)
                         enddo
                        endif

c
c     ----- nuclear attraction
c
                        if(ipass.eq.3) then
                        dum = pi212*aa1
c                       facinv = aa/(fac*pi212)
                        do i = 1 , ij
                         dij(i) = dij(i)*dum
                        enddo
                        aax = aa*ax
                        aay = aa*ay
                        aaz = aa*az

                        do 320 ic = 1 , nat
                         if(omtchg .and. icent.eq.jcent) then
c ****** charge exclusion
                          pnuc = -q(i110-1+(ic-1)*nat + icent)
                         else
                          pnuc = -czan(ic)
                         endif
                           cx = c(1,ic)
                           cy = c(2,ic)
                           cz = c(3,ic)
                           pp = aa*((ax-cx)**2+(ay-cy)**2+(az-cz)**2)
                           if (nroots.le.3) call rt123
                           if (nroots.eq.4) call roots4
                           if (nroots.eq.5) call roots5
                           mm = 0
                           do k = 1 , nroots
                             uu = aa*u(k)
                             ww = w(k)*pnuc
                             tt = done/(aa+uu)
                             t = dsqrt(tt)
                             p0 = (aax+uu*cx)*tt
                             q0 = (aay+uu*cy)*tt
                             r0 = (aaz+uu*cz)*tt
                             in = -5 + mm
                             do i = 1 , lit
                                in = in + 5
                                ni = i
                                do j = 1 , ljt
                                   jn = in + j
                                   nj = j
                                   call stvint
                                   pin(jn) = pint
                                   qin(jn) = qint
                                   rin(jn) = rint*ww
                                enddo
                             enddo
                             mm = mm + 25
                           enddo
                           do i = 1 , ij
                            nnx = ijx(i)
                            nny = ijy(i)
                            nnz = ijz(i)
                            dum = dzero
                            mm = 0
                            do k = 1 , nroots
                             dum = dum + pin(nnx+mm)*qin(nny+mm)
     +                             *rin(nnz+mm)
                             mm = mm + 25
                            enddo
                            g(i) = g(i) + dum*dij(i)
                           enddo
 320                    continue
                        endif

                     end if
c ...
c ...
 390              continue
 400           continue
c
c
c     - s-  at x(i10) on 1st pass
c     - t-  at x(i10) on 2nd pass
c     - h0- at x(i10) on 3rd pass
c
c
c     ----- set up overlap and h-core matrices
c
               max = maxj
               nn = 0
               do i = mini , maxi
                  li = loci + i
                  in = (li*(li-1))/2
                  if (oiandj) max = i
                  do j = minj , max
                     lj = locj + j
                     jn = lj + in
                     nn = nn + 1
                     if(ipass.eq.1) then
                      q(jn-1+i10) = s(nn)*q2
                     else if(ipass.eq.2) then
                      q(jn-1+i10) = gg(nn)*q2
                     else
                      q(jn-1+i10) = g(nn)*q2
                     endif
                  enddo
               enddo
 430        continue
 440     continue
c
c
c     ----- symmetrize 1-e integral arrays
c
         if(nt.gt.1) then
           call sym1e_mp(q(i10),q,iso,nshell)
         endif
c
c   ----- output 1-electron integrals to dumpfile
c
        call sec192_mp(q(i10),l2,ipass)
c
c     ----- now compute dipole moment integrals
c
          call dipxyz_mp(q(i10),ipass)
c
c     ---- load to ionsec
c
      call secdip_mp(q(i10),l2,ipass)

c
c     ---- print if requested
c
      if (out) call prt1e_mp(q(i10),ionsec,ipass)
c
         cpu = cpulft(1)
c
         if (ncall.eq.0) then
            if (outv) write (iwr,6020) cpu
            tim = cpulft(1)
            if (tim.ge.timlim) then
               nindmx = -1
               irest = 1
            end if
            call timana(3)
         end if
c
c     ----- reset core memory -----
c
      call gmem_free(i10)
 
      return
 6010 format (1x,'core assignement'/1x,'i10 = ',i10/1x,
     +        'last = ',i10)
 6020 format (/' integral evaluation complete at ',f8.2,' seconds')
 6030 format (/1x,20('*')/1x,'1-electron integrals'/1x,20('*')
     +        //' commence integral evaluation at ',f8.2,' seconds')
      end
      subroutine sec192_mp(sft,lenbas,ipass)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
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
      dimension sft(lenbas)
c
      common/blkin/potnuc,dx,dy,dz,nato,numo
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
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
c
      data m2/2/
c
      lenb = lensec(lenbas)
      i = lenb*6 + 1
c     S, T, T+V as a function of ipass
      if (ipass.eq.1) then
       call secput(ionsec,m2,i,ib)
       nato = nat
       numo = num
c      header block
       call wrt3(potnuc, 511, ib, idaf)
       ibl3s = ib + 1
       call wrt3(sft,lenbas,ibl3s,idaf)
c      and write subsequent blocks to fill out section
c      including dipole ints
       do loop=1,5
        call wrt3s(sft,lenbas,idaf)
       enddo
      else if (ipass.eq.2) then
       ibl3t = ibl3s + lenb
       call wrt3(sft,lenbas,ibl3t,idaf)
      else
       ibl3f = ibl3t + lenb 
       call wrt3(sft,lenbas,ibl3f,idaf)
      endif
c
      return
      end
      subroutine sym1e_mp(f,h,iso,nshels)
c
c     ----- symmetrize 1e-matrix in f to matrix h
c     ----- matrix h is now allocated in this routine
c     ----- note that h is now allocated in this routine
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
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
      dimension f(*),h(*),iso(nshels,*)
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
      common/junk3/ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720)
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
      common/blkin/pxyz(4),t(35,35),mini,maxi,lit,minj,maxj,ljt,ntr
      dimension mi(48)
c
      i10 = igmem_alloc(nx)
      call vclr(h(i10),1,nx)
c
c     ----- find a block (i,j)
c
      do 520 ii = 1,nshell
      do 140 itr = 1,nt
      ish = iso(ii,itr)
      if (ish .gt. ii) go to 520
  140 mi(itr) = ish
      lit = ktype(ii)
      mini = kmin(ii)
      maxi = kmax(ii)
      loci = kloc(ii)-mini
      do 500 jj = 1,ii
      do 200 itr = 1,nt
      jsh = iso(jj,itr)
      if (jsh .gt. ii) go to 500
      ish = mi(itr)
      if (ish .ge. jsh) go to 180
      n = ish
      ish = jsh
      jsh = n
  180 if (ish .eq. ii .and. jsh .gt. jj) go to 500
  200 continue
      ljt = ktype(jj)
      minj = kmin(jj)
      maxj = kmax(jj)
      locj = kloc(jj)-minj
      oiandj = ii .eq. jj
c
c     ----- find the equivalent blocks -----
c     ----- transfer equivalent block into t-matrix
c     ----- compute (r) t (r)
c     ----- put the result back into the (i,j) block of the h-matrix
c
      jmax = maxj
      do 300 itr = 1,nt
      ntr = itr
      kk = mi(itr)
      ll = iso(jj,itr)
      lock = kloc(kk)-kmin(kk)
      locl = kloc(ll)-kmin(ll)
      do 260 k = mini,maxi
      lck = lock+k
      if (oiandj) jmax = k
      do 260 l = minj,jmax
      kl = iky(max(lck,locl+l))+min(lck,locl+l)
      t(k,l) = f(kl)
      if (oiandj) t(l,k) = f(kl)
  260 continue
      call rhrrhr
      do 280 i = mini,maxi
      lci = iky(loci+i)+locj
      if (oiandj) jmax = i
      do 280 j = minj,jmax
      ij = lci+j-1
  280 h(i10+ij) = h(i10+ij)+t(i,j)
c
c     ----- for each block (k,l) equivalent to (i,j)
c     ----- find the transformation that maps (k,l) into (i,j)
c     ----- compute (r) t (r)
c     ----- put the result back into the (k,l) block of the h-matrix
c
  300 continue
      do 480 itr = 1,nt
      if (itr .eq. 1) go to 480
      kk = mi(itr)
      ll = iso(jj,itr)
      if (kk .ge. ll) go to 320
      k = ll
      l = kk
      go to 340
  320 k = kk
      l = ll
  340 if (k .eq. ii .and. l .eq. jj) go to 480
      ntr = itr+1
      if (ntr .gt. nt) go to 400
      do 380 it = ntr,nt
      i = mi(it)
      j = iso(jj,it)
      if (i .ge. j) go to 360
      ij = i
      i = j
      j = ij
  360 if (i .eq. k .and. j .eq. l) go to 480
  380 continue
  400 continue
      ntr = invt(itr)
      do 420 i = mini,maxi
      lci = iky(loci+i)+locj-1
      if (oiandj) jmax = i
      do 420 j = minj,jmax
      t(i,j) = h(i10+lci+j)
      if (oiandj) t(j,i) = h(i10+lci+j)
  420 continue
      call rhrrhr
      lock = kloc(kk)-kmin(kk)
      locl = kloc(ll)-kmin(ll)
      do 460 k = mini,maxi
      lck = lock+k
      if (oiandj) jmax = k
      do 460 l = minj,jmax
      kl = iky(max(lck,locl+l))+min(lck,locl+l)
  460 h(i10-1+kl) = t(k,l)
  480 continue
  500 continue
  520 continue
      dum = 1.0d0/ dble(nt)
      call vsmul(h(i10),1,dum,f,1,nx)
c
c     ----- reset core memory -----
c
      call gmem_free(i10)
      return
      end
      subroutine dipxyz_mp(ps,ipass)
c     compute integrals in separate passes
c     so no more than 1 triangle is used
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
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
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
      common/junk/ssgg(450),
     + pint0,qint0,rint0,pintx,qinty,rintz,t,p0,q0,r0,
     + pi,qi,ri,pj,qj,rj,ni,nj
      common/blkin/pxyz(4),sx(225),sy(225),sz(225),
     + dij(225),pin(125),qin(125),rin(125),
     + ijx(225),ijy(225),ijz(225)
      dimension ps(*)
      dimension ix(35),iy(35),iz(35),jx(35),jy(35),jz(35)
      data dzero /0.0d0/, done/1.0d0/
      data sqrt3 /1.73205080756888d0/
      data sqrt5 /2.23606797749979d0/
      data sqrt7 /2.64575131106459d0/
      data rln10 /2.30258d0/
c
      data jx / 0, 1, 0, 0, 2, 0, 0, 1, 1, 0,
     +          3, 0, 0, 2, 2, 1, 0, 1, 0, 1,
     +          4, 0, 0, 3, 3, 1, 0, 1, 0, 2,
     +          2, 0, 2, 1, 1/
      data ix / 1, 6, 1, 1,11, 1, 1, 6, 6, 1,
     +         16, 1, 1,11,11, 6, 1, 6, 1, 6,
     +         21, 1, 1,16,16, 6, 1, 6, 1,11,
     +         11, 1,11, 6, 6/
      data jy / 0, 0, 1, 0, 0, 2, 0, 1, 0, 1,
     +          0, 3, 0, 1, 0, 2, 2, 0, 1, 1,
     +          0, 4, 0, 1, 0, 3, 3, 0, 1, 2,
     +          0, 2, 1, 2, 1/
      data iy / 1, 1, 6, 1, 1,11, 1, 6, 1, 6,
     +          1,16, 1, 6, 1,11,11, 1, 6, 6,
     +          1,21, 1, 6, 1,16,16, 1, 6,11,
     +          1,11, 6,11, 6/
      data jz / 0, 0, 0, 1, 0, 0, 2, 0, 1, 1,
     +          0, 0, 3, 0, 1, 0, 1, 2, 2, 1,
     +          0, 0, 4, 0, 1, 0, 1, 3, 3, 0,
     +          2, 2, 1, 1, 2/
      data iz / 1, 1, 1, 6, 1, 1,11, 1, 6, 6,
     +          1, 1,16, 1, 6, 1, 6,11,11, 6,
     +          1, 1,21, 1, 6, 1, 6,16,16, 1,
     +         11,11, 6, 6,11/
c
      tol = rln10*itol
      onorm = normf.ne.1 .or. normp.ne.1
      do 310 ii = 1 , nshell
c
c     ----- ishell
c
         i = katom(ii)
         pi = c(1,i)
         qi = c(2,i)
         ri = c(3,i)
         i1 = kstart(ii)
         i2 = i1 + kng(ii) - 1
         lit = ktype(ii)
         mini = kmin(ii)
         maxi = kmax(ii)
         loci = kloc(ii) - mini
c
c     ----- jshell
c
         do 300 jj = 1 , ii
            j = katom(jj)
            pj = c(1,j)
            qj = c(2,j)
            rj = c(3,j)
            j1 = kstart(jj)
            j2 = j1 + kng(jj) - 1
            ljt = ktype(jj)
            minj = kmin(jj)
            maxj = kmax(jj)
            locj = kloc(jj) - minj
            rr = (pi-pj)**2 + (qi-qj)**2 + (ri-rj)**2
            oiandj = ii.eq.jj
c
c     ----- prepare indices for pairs of (i,j) functions
c
            ij = 0
            max = maxj
            do 30 i = mini , maxi
               nnx = ix(i)
               nny = iy(i)
               nnz = iz(i)
               if (oiandj) max = i
               do 20 j = minj , max
                  ij = ij + 1
                  ijx(ij) = nnx + jx(j)
                  ijy(ij) = nny + jy(j)
                  ijz(ij) = nnz + jz(j)
 20            continue
 30         continue
            do 40 i = 1 , ij
               sx(i) = dzero
               sy(i) = dzero
               sz(i) = dzero
 40         continue
c
c     ----- i primitive
c
            jgmax = j2
            do 270 ig = i1 , i2
               ai = ex(ig)
               arri = ai*rr
               axi = ai*pi
               ayi = ai*qi
               azi = ai*ri
               csi = cs(ig)
               cpi = cp(ig)
               cdi = cd(ig)
               cfi = cf(ig)
               cgi = cg(ig)
c
c     ----- j primtive
c
               if (oiandj) jgmax = ig
               do 260 jg = j1 , jgmax
                  aj = ex(jg)
                  aa = ai + aj
                  dum = aj*arri/aa
                  if (dum.le.tol) then
                     fac = dexp(-dum)
                     csj = cs(jg)
                     cpj = cp(jg)
                     cdj = cd(jg)
                     cfj = cf(jg)
                     cgj = cg(jg)
                     ax = (axi+aj*pj)/aa
                     ay = (ayi+aj*qj)/aa
                     az = (azi+aj*rj)/aa
c
c     ----- density factor
c
                     odoub = oiandj .and. ig.ne.jg
                     max = maxj
                     nn = 0
                     do 220 i = mini , maxi
                        go to (50,60,120,120,
     +                      70,120,120,80,120,120,
     +                      90,120,120,100,120,120,120,120,120,110,
     +                      112,120,120,114,120,120,120,120,120,116,
     +                      120,120,118,120,120),i
 50                     dum1 = csi*fac
                        go to 120
 60                     dum1 = cpi*fac
                        go to 120
 70                     dum1 = cdi*fac
                        go to 120
 80                     if (onorm) dum1 = dum1*sqrt3
                        go to 120
 90                     dum1 = cfi*fac
                        go to 120
 100                    if (onorm) dum1 = dum1*sqrt5
                        go to 120
 110                    if (onorm) dum1 = dum1*sqrt3
                        go to 120
 112                    dum1 = cgi*fac
                        go to 120
 114                    if (onorm) dum1 = dum1*sqrt7
                        go to 120
 116                    if (onorm) dum1 = dum1*sqrt5/sqrt3
                        go to 120
 118                    if (onorm) dum1 = dum1*sqrt3
 120                    if (oiandj) max = i
                        do 210 j = minj , max
                           go to (130,140,200,200,
     +                      150,200,200,160,200,200,
     +                      170,200,200,180,200,200,200,200,200,190,
     +                      192,200,200,194,200,200,200,200,200,196,
     +                      200,200,198,200,200),j
 130                       dum2 = dum1*csj
                           if (odoub) then
                              if (i.gt.1) then
                                 dum2 = dum2 + csi*cpj*fac
                              else
                                 dum2 = dum2 + dum2
                              end if
                           end if
                           go to 200
 140                       dum2 = dum1*cpj
                           if (odoub) dum2 = dum2 + dum2
                           go to 200
 150                       dum2 = dum1*cdj
                           if (odoub) dum2 = dum2 + dum2
                           go to 200
 160                       if (onorm) dum2 = dum2*sqrt3
                           go to 200
 170                       dum2 = dum1*cfj
                           if (odoub) dum2 = dum2 + dum2
                           go to 200
 180                       if (onorm) dum2 = dum2*sqrt5
                           go to 200
 190                       if (onorm) dum2 = dum2*sqrt3
                           go to 200
 192                       dum2 = dum1*cgj
                           if (odoub) dum2 = dum2+dum2
                           go to 200
 194                       if (onorm) dum2 = dum2*sqrt7
                           go to 200
 196                       if (onorm) dum2 = dum2*sqrt5/sqrt3
                           go to 200
 198                       if (onorm) dum2 = dum2*sqrt3
 200                       nn = nn + 1
                           dij(nn) = dum2
 210                    continue
 220                 continue
c
c     ----- dipole moment integrals -----
c
                     t = dsqrt(aa)
                     tinv = done/t
                     p0 = ax
                     q0 = ay
                     r0 = az
                     in = -5
                     do 240 i = 1 , lit
                        in = in + 5
                        ni = i
                        do 230 j = 1 , ljt
                           jn = in + j
                           nj = j
                           call dipint
                           pin(jn) = pint0*tinv
                           qin(jn) = qint0*tinv
                           rin(jn) = rint0*tinv
                           pin(jn+25) = pintx*tinv
                           qin(jn+25) = qinty*tinv
                           rin(jn+25) = rintz*tinv
 230                    continue
 240                 continue
                     do 250 i = 1 , ij
                       nnx = ijx(i)
                       nny = ijy(i)
                       nnz = ijz(i)
                       sx(i) = sx(i) + dij(i)
     +                         *(pin(nnx+25)*qin(nny)*rin(nnz))
                       sy(i) = sy(i) + dij(i)
     +                         *(pin(nnx)*qin(nny+25)*rin(nnz))
                       sz(i) = sz(i) + dij(i)
     +                         *(pin(nnx)*qin(nny)*rin(nnz+25))
 250                 continue
                  end if
 260           continue
 270        continue
c
c     ----- set up dipole moment matrices -----
c
            max = maxj
            nn = 0
            if(ipass.eq.1) then
             do i = mini , maxi
                li = loci + i
                in = (li*(li-1))/2
                if (oiandj) max = i
                do j = minj , max
                   lj = locj + j
                   jn = lj + in
                   nn = nn + 1
                   ps(jn) = sx(nn)
                enddo
             enddo
            else if (ipass.eq.2) then
             do i = mini , maxi
                li = loci + i
                in = (li*(li-1))/2
                if (oiandj) max = i
                do j = minj , max
                   lj = locj + j
                   jn = lj + in
                   nn = nn + 1
                   ps(jn) = sy(nn)
                enddo
             enddo
            else
             do i = mini , maxi
                li = loci + i
                in = (li*(li-1))/2
                if (oiandj) max = i
                do j = minj , max
                   lj = locj + j
                   jn = lj + in
                   nn = nn + 1
                   ps(jn) = sz(nn)
                enddo
             enddo
            endif
 300     continue
 310  continue
      return
      end
      subroutine secdip_mp(pqr,lenbas,ipass)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
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
      common/blkin/potnuc(4),nato,numo
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
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
c
      dimension pqr(lenbas)
      data m2/2/
c
c     call secget(ionsec,m2,ib)
c     call rdedx(potnuc,511,ib,idaf)
       do ij = 1,lenbas
        pqr(ij) = -pqr(ij)
       enddo
c
      lenb = lensec(lenbas)
      if(ipass.eq.1) then
       ibl3x = ibl3s + 3 * lenb
       call wrt3(pqr,lenbas,ibl3x,idaf)
      else if(ipass.eq.2) then
       ibl3y = ibl3t + 3 * lenb
       call wrt3(pqr,lenbas,ibl3y,idaf)
      else
       ibl3z = ibl3f + 3 * lenb 
       call wrt3(pqr,lenbas,ibl3z,idaf)
      endif
c
      return
      end
      subroutine prt1e_mp(q,ionsec,ipass)
c
c     ----- print 1-electron integrals resident in section
c           ionsec as a sequence of matrices   ----
c
      implicit real*8  (a-h,p-w), integer (i-n), logical (o)
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
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
      common/blkin/corev(512),potnuc(4),pint(508)
      dimension yprin(6),o1e(6),q(*)
      data yprin/'  s-','  t-','t+v-','  x-','  y-','  z-'/
c
c     process the 1e-integrals
c
      do loop =1,6
       o1e(loop) = .false.
      enddo
c
c     first restore header block
c
      call getmat(q,q,q,q,q,q,potnuc,num,o1e,ionsec)
c
      write (iwr,6010) ionsec
      write (iwr,6020) potnuc
c
c     now process each matrix in turn (S, X, T, Y, T+V, Z)
c
      if(ipass.eq.1) then
       imat=1
       if(oprn(31)) then
        o1e(imat) = .true.
        call getmat(q,q,q,q,q,q,potnuc,num,o1e,ionsec)
         write(iwr,6030)yprin(imat)
         call writel(q,num)
        o1e(imat) = .false.
       endif
       imat=4
       if(oprn(34)) then
        o1e(imat) = .true.
        call getmat(q,q,q,q,q,q,potnuc,num,o1e,ionsec)
         write(iwr,6030)yprin(imat)
         call writel(q,num)
        o1e(imat) = .false.
       endif
      else if(ipass.eq.2) then
       imat=2
       if(oprn(32)) then
        o1e(imat) = .true.
        call getmat(q,q,q,q,q,q,potnuc,num,o1e,ionsec)
         write(iwr,6030)yprin(imat)
         call writel(q,num)
        o1e(imat) = .false.
       endif
       imat=5
       if(oprn(35)) then
        o1e(imat) = .true.
        call getmat(q,q,q,q,q,q,potnuc,num,o1e,ionsec)
         write(iwr,6030)yprin(imat)
         call writel(q,num)
        o1e(imat) = .false.
       endif
      else
       imat=3
       if(oprn(33)) then
        o1e(imat) = .true.
        call getmat(q,q,q,q,q,q,potnuc,num,o1e,ionsec)
         write(iwr,6030)yprin(imat)
         call writel(q,num)
        o1e(imat) = .false.
       endif
       imat=6
       if(oprn(36)) then
        o1e(imat) = .true.
        call getmat(q,q,q,q,q,q,potnuc,num,o1e,ionsec)
         write(iwr,6030)yprin(imat)
         call writel(q,num)
        o1e(imat) = .false.
       endif
c
      endif

      return
 6030  format(/1x,104('*')//
     * 35x,a4,'matrix over gaussian basis set'/
     * 35x,34('-')/)
 6010 format (/1x,104('=')/45x,
     +        'list of 1-electron integrals in section',i6/45x,43('-'))
 6020 format (/' potnuc,dx,dy,dz = ',4f19.8)
      end

c
c to support additional code in server
c

      subroutine amints(xl,yl,zl,dd,v,space,otran)
c
c     angular momentum integrals
c
      implicit real*8  (a-h,o-z)
      dimension xl(*),space(*),dd(*),v(*),yl(*),zl(*)
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
      logical out,norm,iandj,double
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      common/intdip/xint0,yint0,zint0,xintx,yinty,zintz,
     1t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj,origx,origy,origz
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
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
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
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
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
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
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      common/small/dij(100),
     1    xs(25),ys(25),zs(25),xd(25),yd(25),zd(25),
     2    xdip(25),ydip(25),zdip(25),
     3    slx(100),sly(100),slz(100),
     4    ijx(100),ijy(100),ijz(100)
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
      data ndim/5/
      logical otran
c
c     origin for lx,ly,lz
c
      origx = gx
      origy = gy
      origz = gz
      tol = 2.30258d0*itol
      out = odebug(1)
      norm = normf.ne.1 .or. normp.ne.1
c     ----- ishell
      do 110 ii = 1 , nshell
         iat = katom(ii)
         xi = c(1,iat)
         yi = c(2,iat)
         zi = c(3,iat)
         i1 = kstart(ii)
         i2 = i1 + kng(ii) - 1
         lit = ktype(ii)
         mini = kmin(ii)
         maxi = kmax(ii)
         loci = kloc(ii) - mini
c     ----- jshell
         do 100 jj = 1 , ii
            jat = katom(jj)
            xj = c(1,jat)
            yj = c(2,jat)
            zj = c(3,jat)
            j1 = kstart(jj)
            j2 = j1 + kng(jj) - 1
            ljt = ktype(jj)
            minj = kmin(jj)
            maxj = kmax(jj)
            locj = kloc(jj) - minj
            rr = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
c     ----- prepare indices for pairs of (i,j) functions
            iandj = ii.eq.jj
            call indexa(ijx,ijy,ijz,ij,mini,maxi,minj,maxj,iandj,ndim,
     +           1,1)
            do 20 i = 1 , ij
               slx(i) = 0.0d0
               sly(i) = 0.0d0
               slz(i) = 0.0d0
 20         continue
c     ----- i primitive
c           jgmax = j2
            do 70 ig = i1 , i2
               ai = ex(ig)
               arri = ai*rr
               axi = ai*xi
               ayi = ai*yi
               azi = ai*zi
               csi = cs(ig)
               cpi = cp(ig)
               cdi = cd(ig)
               cfi = cf(ig)
               cgi = cg(ig)
c     ----- j primtive
               do 60 jg = j1 , j2
                  aj = ex(jg)
                  aa = ai + aj
                  aainv = 1.0d0/aa
                  dum = aj*arri*aainv
                  if (dum.le.tol) then
                     fac = dexp(-dum)
                     csj = cs(jg)*fac
                     cpj = cp(jg)*fac
                     cdj = cd(jg)*fac
                     cfj = cf(jg)*fac
                     cgj = cg(jg)*fac
                     ax = (axi+aj*xj)*aainv
                     ay = (ayi+aj*yj)*aainv
                     az = (azi+aj*zj)*aainv
c     ----- density factor
                     double = .false.
                     call denfan(dij,csi,cpi,cdi,cfi,cgi,
     +                               csj,cpj,cdj,cfj,cgj,
     +                           mini,maxi,minj,maxj,iandj,double,norm)
c     ----- overlap
                     t = dsqrt(aa)
                     tinv = 1.0d0/t
                     x0 = ax
                     y0 = ay
                     z0 = az
                     ljt1 = ljt + 1
                     in = -ndim
                     do 40 i = 1 , lit
                        in = in + ndim
                        ni = i
                        do 30 j = 1 , ljt1
                           jn = in + j
                           nj = j
                           call dmsint()
                           xs(jn) = xint0*tinv
                           ys(jn) = yint0*tinv
                           zs(jn) = zint0*tinv
                           xdip(jn) = xintx*tinv
                           ydip(jn) = yinty*tinv
                           zdip(jn) = zintz*tinv
 30                     continue
 40                  continue
                     call oneld(xs,ys,zs,xd,yd,zd,aj,lit,ljt,2,ndim)
                     do 50 i = 1 , ij
                        nnx = ijx(i)
                        nny = ijy(i)
                        nnz = ijz(i)
                        dum = dij(i)
                        slx(i) = slx(i) + dum*xs(nnx)
     +                           *(ydip(nny)*zd(nnz)-yd(nny)*zdip(nnz))
                        sly(i) = sly(i) + dum*ys(nny)
     +                           *(xd(nnx)*zdip(nnz)-xdip(nnx)*zd(nnz))
                        slz(i) = slz(i) + dum*zs(nnz)
     +                           *(xdip(nnx)*yd(nny)-xd(nnx)*ydip(nny))
c
c     lx,ly,lz are calculated from dipole integrals and derivative
c     overlap integrals,
c     lz   =  (x d/dy - y d/dx )  etc
c
c
 50                  continue
                  end if
 60            continue
 70         continue
            max = maxj
            n = 0
            do 90 i = mini , maxi
               li = loci + i
               in = iky(li)
               if (iandj) max = i
               do 80 j = minj , max
                  lj = locj + j
                  jn = in + lj
                  n = n + 1
                  xl(jn) = -0.5d0*slx(n)
                  yl(jn) = -0.5d0*sly(n)
                  zl(jn) = -0.5d0*slz(n)
 80            continue
 90         continue
 100     continue
 110  continue
      if (out) then
         write (iwr,*) 'angular momentum in a.o. basis'
         write (iwr,*)
     +              '(note : these matrices are really anti-symmetric!)'
         write (iwr,*)
         call prtris(xl,num,iwr)
         call prtris(yl,num,iwr)
         call prtris(zl,num,iwr)
      end if
      if (.not.otran) return
c
c
c     get vectors and transform to mo basis
c
c
      ityp = 0
      call secget(isect(8),ityp,isecv)
      call rdedx(v,num*ncoorb,isecv+mvadd,ifild)
      call dcopy(nx,xl,1,dd,1)
      call vclr(xl,1,nx)
      call skwtr(xl,dd,v,space,iky,ncoorb,num,num)
      call dcopy(nx,yl,1,dd,1)
      call vclr(yl,1,nx)
      call skwtr(yl,dd,v,space,iky,ncoorb,num,num)
      call dcopy(nx,zl,1,dd,1)
      call vclr(zl,1,nx)
      call skwtr(zl,dd,v,space,iky,ncoorb,num,num)
      len = lensec(nx)
      lds(isect(87)) = nx
      lds(isect(88)) = nx
      lds(isect(89)) = nx
c
c
c     write to sections 87,88,89 of dumpfile
c
c
      call secput(isect(87),87,len,isec87)
      call secput(isect(88),88,len,isec88)
      call secput(isect(89),89,len,isec89)
      call wrt3(xl,nx,isec87,ifild)
      call wrt3(yl,nx,isec88,ifild)
      call wrt3(zl,nx,isec89,ifild)
      call revind
      if (.not.out) return
      write (iwr,*) 'angular momentum in m.o. basis'
      write (iwr,*) '(note : these matrices are really anti-symmetric!)'
      write (iwr,*)
      call prtris(xl,num,iwr)
      call prtris(yl,num,iwr)
      call prtris(zl,num,iwr)
      return
      end

      subroutine lmints(xl,yl,zl,dd,v,space,otran)
c
c     linear momentum integrals
c
      implicit real*8  (a-h,o-z)
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
      logical out,norm,iandj,double
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      common/intdip/xint0,yint0,zint0,xintx,yinty,zintz,
     1t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj,origx,origy,origz
c
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
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
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
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
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
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
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      common/small/dij(100),
     1    xs(25),ys(25),zs(25),xd(25),yd(25),zd(25),
     2    xdip(25),ydip(25),zdip(25),
     3    px(100),py(100),pz(100),
     4    ijx(100),ijy(100),ijz(100)
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
c
      dimension xl(*),yl(*),zl(*),v(*),dd(*),space(*)
      data ndim/5/
      logical otran
c
c
c
      origx = gx
      origy = gy
      origz = gz
      tol = 2.30258d0*itol
      out = odebug(1)
      norm = normf.ne.1 .or. normp.ne.1
c     ----- ishell
      do 110 ii = 1 , nshell
         iat = katom(ii)
         xi = c(1,iat)
         yi = c(2,iat)
         zi = c(3,iat)
         i1 = kstart(ii)
         i2 = i1 + kng(ii) - 1
         lit = ktype(ii)
         mini = kmin(ii)
         maxi = kmax(ii)
         loci = kloc(ii) - mini
c     ----- jshell
         do 100 jj = 1 , ii
            jat = katom(jj)
            xj = c(1,jat)
            yj = c(2,jat)
            zj = c(3,jat)
            j1 = kstart(jj)
            j2 = j1 + kng(jj) - 1
            ljt = ktype(jj)
            minj = kmin(jj)
            maxj = kmax(jj)
            locj = kloc(jj) - minj
            rr = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
c     ----- prepare indices for pairs of (i,j) functions
            iandj = ii.eq.jj
            call indexa(ijx,ijy,ijz,ij,mini,maxi,minj,maxj,iandj,ndim,
     +           1,1)
            do 20 i = 1 , ij
               px(i) = 0.0d0
               py(i) = 0.0d0
               pz(i) = 0.0d0
 20         continue
c     ----- i primitive
            do 70 ig = i1 , i2
               ai = ex(ig)
               arri = ai*rr
               axi = ai*xi
               ayi = ai*yi
               azi = ai*zi
               csi = cs(ig)
               cpi = cp(ig)
               cdi = cd(ig)
               cfi = cf(ig)
               cgi = cg(ig)
c     ----- j primtive
               do 60 jg = j1 , j2
                  aj = ex(jg)
                  aa = ai + aj
                  aainv = 1.0d0/aa
                  dum = aj*arri*aainv
                  if (dum.le.tol) then
                     fac = dexp(-dum)
                     csj = cs(jg)*fac
                     cpj = cp(jg)*fac
                     cdj = cd(jg)*fac
                     cfj = cf(jg)*fac
                     cgj = cg(jg)*fac
                     ax = (axi+aj*xj)*aainv
                     ay = (ayi+aj*yj)*aainv
                     az = (azi+aj*zj)*aainv
c     ----- density factor
                     double = .false.
                     call denfan(dij,csi,cpi,cdi,cfi,cgi,
     +                               csj,cpj,cdj,cfj,cgj,
     +                           mini,maxi,minj,maxj,iandj,double,norm)
c     ----- overlap
                     t = dsqrt(aa)
                     tinv = 1.0d0/t
                     x0 = ax
                     y0 = ay
                     z0 = az
                     in = -ndim
                     do 40 i = 1 , lit
                        in = in + ndim
                        ni = i
                        do 30 j = 1 , ljt + 1
                           jn = in + j
                           nj = j
                           call dmsint()
                           xs(jn) = xint0*tinv
                           ys(jn) = yint0*tinv
                           zs(jn) = zint0*tinv
 30                     continue
 40                  continue
                     call oneld(xs,ys,zs,xd,yd,zd,aj,lit,ljt,2,ndim)
                     do 50 i = 1 , ij
                        nnx = ijx(i)
                        nny = ijy(i)
                        nnz = ijz(i)
                        dum = dij(i)
                        px(i) = px(i) + dum*xd(nnx)*ys(nny)*zs(nnz)
                        py(i) = py(i) + dum*xs(nnx)*yd(nny)*zs(nnz)
                        pz(i) = pz(i) + dum*xs(nnx)*ys(nny)*zd(nnz)
c
c
 50                  continue
                  end if
 60            continue
 70         continue
            max = maxj
            n = 0
            do 90 i = mini , maxi
               li = loci + i
               if (iandj) max = i
               do 80 j = minj , max
                  lj = locj + j
                  n = n + 1
                  xl(iky(li)+lj) = px(n)
                  yl(iky(li)+lj) = py(n)
                  zl(iky(li)+lj) = pz(n)
 80            continue
 90         continue
 100     continue
 110  continue
      if (out) then
       write (iwr,*)'linear momentum in a.o. basis'
       write (iwr,*)'(note : these matrices are really anti-symmetric!)'
       write (iwr,*)
         call prtris(xl,num,iwr)
         call prtris(yl,num,iwr)
         call prtris(zl,num,iwr)
      end if
      if (.not.otran) return
c
c
c     get vectors and transform to mo basis
c
c
      ityp = 0
      call secget(isect(8),ityp,isecv)
      call rdedx(v,num*ncoorb,isecv+mvadd,ifild)
      call dcopy(nx,xl,1,dd,1)
      call vclr(xl,1,nx)
      call skwtr(xl,dd,v,space,iky,ncoorb,num,num)
      call dcopy(nx,yl,1,dd,1)
      call vclr(yl,1,nx)
      call skwtr(yl,dd,v,space,iky,ncoorb,num,num)
      call dcopy(nx,zl,1,dd,1)
      call vclr(zl,1,nx)
      call skwtr(zl,dd,v,space,iky,ncoorb,num,num)
      len = lensec(nx)
      lds(isect(84)) = nx
      lds(isect(85)) = nx
      lds(isect(86)) = nx
c
c
c     write to sections 84,85,86 of dumpfile
c
c
      call secput(isect(84),84,len,isec87)
      call secput(isect(85),85,len,isec88)
      call secput(isect(86),86,len,isec89)
      call wrt3(xl,nx,isec87,ifild)
      call wrt3(yl,nx,isec88,ifild)
      call wrt3(zl,nx,isec89,ifild)
      call revind
      if (.not.out) return
      write (iwr,*) 'linear momentum in m.o. basis'
      write (iwr,*)'(note : these matrices are really anti-symmetric!)'
      write (iwr,*)
      call prtris(xl,num,iwr)
      call prtris(yl,num,iwr)
      call prtris(zl,num,iwr)
      return
      end

      subroutine dmsint
      implicit real*8  (a-h,o-z)
      common/intdip/xint0,yint0,zint0,xintx,yinty,zintz,t,x0,y0,z0,
     1 xi,yi,zi,xj,yj,zj,ni,nj,gx,gy,gz
      common/hermit/h1,h2(2),h3(3),h4(4),h5(5),h6(6),h7(7)
      common/wermit/w1,w2(2),w3(3),w4(4),w5(5),w6(6),w7(7)
      dimension h(21),w(21),min(7),max(7)
      equivalence (h(1),h1),(w(1),w1)
      data min /1,2,4,7,11,16,22/
      data max /1,3,6,10,15,21,28/
      data zero /0.0d0/
      xint0 = zero
      yint0 = zero
      zint0 = zero
      xintx = zero
      yinty = zero
      zintz = zero
      npts = (ni+nj-2+2)/2 + 1
      imin = min(npts)
      imax = max(npts)
      do 160 i = imin , imax
         dum = w(i)
         px = dum
         py = dum
         pz = dum
         dum = h(i)/t
         ptx = dum + x0
         pty = dum + y0
         ptz = dum + z0
         ax = ptx - xi
         ay = pty - yi
         az = ptz - zi
         bx = ptx - xj
         by = pty - yj
         bz = ptz - zj
         go to (70,60,50,40,30,20,10) , ni
 10      px = px*ax
         py = py*ay
         pz = pz*az
 20      px = px*ax
         py = py*ay
         pz = pz*az
 30      px = px*ax
         py = py*ay
         pz = pz*az
 40      px = px*ax
         py = py*ay
         pz = pz*az
 50      px = px*ax
         py = py*ay
         pz = pz*az
 60      px = px*ax
         py = py*ay
         pz = pz*az
 70      go to (150,140,130,120,110,100,90,80) , nj
 80      px = px*bx
         py = py*by
         pz = pz*bz
 90      px = px*bx
         py = py*by
         pz = pz*bz
 100     px = px*bx
         py = py*by
         pz = pz*bz
 110     px = px*bx
         py = py*by
         pz = pz*bz
 120     px = px*bx
         py = py*by
         pz = pz*bz
 130     px = px*bx
         py = py*by
         pz = pz*bz
 140     px = px*bx
         py = py*by
         pz = pz*bz
 150     xint0 = xint0 + px
         yint0 = yint0 + py
         zint0 = zint0 + pz
         xintx = xintx + px*(ptx-gx)
         yinty = yinty + py*(pty-gy)
         zintz = zintz + pz*(ptz-gz)
 160  continue
      return
      end
      subroutine skwtr(h,f,c,t,ia,nvec,nbasis,ndim)
c
c    transformation routine for case when (i/ /j) = - (j/ /i)
c
      implicit real*8  (a-h,o-z)
      dimension h(*),f(*),c(*),t(*),ia(*)
      data zero /0.0d0/
      nj = 0
      ij = 0
      do 50 j = 1 , nvec
         do 30 k = 1 , nbasis
            ckj = c(k+nj)
            kk = ia(k)
            dum = zero
            kless1 = k - 1
            if (kless1.gt.0) then
               do 20 l = 1 , kless1
                  fkl = f(kk+l)
                  dum = dum + fkl*c(nj+l)
                  t(l) = t(l) - fkl*ckj
 20            continue
            end if
            t(k) = dum
 30      continue
         ni = 0
         do 40 i = 1 , j
            ij = ij + 1
            h(ij) = ddot(nbasis,c(ni+1),1,t(1),1)
            ni = ni + ndim
 40      continue
         nj = nj + ndim
 50   continue
      return
      end
