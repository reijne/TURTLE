c     deck=analb
c ******************************************************
c ******************************************************
c             =   analb =
c ******************************************************
c ******************************************************
c
**==aainer.f
      function aainer(r,s,t,l,m,n,d,dp,fnu,fn,fd)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer r,s,t,u,v,w,uvwt,rstt,uvwth
      dimension d(3),dp(4),fnu(7),fn(7),fd(7)
      rstt = r + s + t
      lmnt = l + m + n
      lmnrst = rstt + lmnt
      lh = l/2
      mh = m/2
      nh = n/2
      lhp = lh + 1
      mhp = mh + 1
      nhp = nh + 1
      aainer = 0.0d0
      do 70 ii = 1 , lhp
         i = ii - 1
         il = l - 2*i + 1
         ilr = r + il
         ilrh = (ilr-1)/2
         ilrhp = ilrh + 1
         fi = fn(ilr)*fd(il)*fd(i+1)
         do 60 jj = 1 , mhp
            j = jj - 1
            jm = m - 2*j + 1
            jms = s + jm
            jmsh = (jms-1)/2
            jmshp = jmsh + 1
            fij = fn(jms)*fd(jm)*fd(j+1)*fi
            do 50 kk = 1 , nhp
               k = kk - 1
               kn = n - 2*k + 1
               knt = t + kn
               knth = (knt-1)/2
               knthp = knth + 1
               ijkt = i + j + k
               fijk = fn(knt)*fd(kn)*fd(k+1)*dp(ijkt+1)*fij
               lmrsij = lmnrst - 2*ijkt
               do 40 iu = 1 , ilrhp
                  u = iu - 1
                  ilru = ilr - 2*u
                  fu = fd(u+1)*fd(ilru)
                  if (dabs(d(1)).gt.0.0000001d0) then
                     fu = fu*d(1)**(ilru-1)
                  else if (ilru.gt.1) then
                     go to 40
                  end if
                  do 30 iv = 1 , jmshp
                     v = iv - 1
                     jmsv = jms - 2*v
                     fuv = fu*fd(v+1)*fd(jmsv)
                     if (dabs(d(2)).gt.0.0000001d0) then
                        fuv = fuv*d(2)**(jmsv-1)
                     else if (jmsv.gt.1) then
                        go to 30
                     end if
                     do 20 iw = 1 , knthp
                        w = iw - 1
                        kntw = knt - 2*w
                        fuvw = fuv*fd(w+1)*fd(kntw)
                        if (dabs(d(3)).gt.0.0000001d0) then
                           fuvw = fuvw*d(3)**(kntw-1)
                        else if (kntw.gt.1) then
                           go to 20
                        end if
                        uvwt = u + v + w
                        uvwth = uvwt/2
                        if (2*uvwth.ne.uvwt) then
                           fuvw = -fuvw
                        end if
                        nuindx = lmrsij - uvwt
                        fuvw = fuvw*fnu(nuindx+1)*dp(uvwt+1)
                        aainer = fijk*fuvw + aainer
 20                  continue
 30               continue
 40            continue
 50         continue
 60      continue
 70   continue
c
      return
      end
**==agnuc.f
      subroutine agnuc(noc,vlist,c,dint,nocmx)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension dint(7),vlist(nocmx,4),c(3)
      common/miscop/qnuc(6)
c...
c...compute nuclear contributions to diamagnetic suceptibilities.
c...ch(xx)=1.5*(rr-xx), similarly for ch(yy) and ch(zz),
c...ch(xy)=-1.5*xy, similarly for ch(xz) and ch(yz)., also
c...ch(rr)=ch(average)=rr.
c...
      do 20 i = 1 , 6
         qnuc(i) = 0.0d0
 20   continue
c...
      do 50 i = 1 , noc
         mu = 0
         do 40 j = 1 , 3
            do 30 k = 1 , j
               mu = mu + 1
               qnuc(mu) = qnuc(mu) + vlist(i,4)*(vlist(i,j)-c(j))
     +                    *(vlist(i,k)-c(k))
 30         continue
 40      continue
 50   continue
c...
      rr = qnuc(1) + qnuc(3) + qnuc(6)
      dint(1) = 1.5d0*(rr-qnuc(1))
      dint(2) = 1.5d0*(rr-qnuc(3))
      dint(3) = 1.5d0*(rr-qnuc(6))
      dint(4) = -1.5d0*qnuc(2)
      dint(5) = -1.5d0*qnuc(4)
      dint(6) = -1.5d0*qnuc(5)
      dint(7) = rr
c...
      return
      end
**==allrm.f
      subroutine allrm(iwr)
c
c     one configuration
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
      common/junk/pqn(99),expa(99),ctran(mxprms,30),ajmn(12),czn,
     *expmn,c(15,15),itran(mxprms,30),nrow(30),nbcon(4),nbas(4),
     *nsym,nosh(4),ncsh(4),nccup(4),ibase(33),
     *popul(4),facto(30),fatt(30),thre(15),a(15,15),eivr(15,15),
     *cm2scf(15,15),cm1scf(15,15),eps(15),gmin(15,15),smin(15,15),
     *populs(4),term,energ,pot,cin,vir,thrscf,enorm,
     *red,expdif,encor,
     *nsize(10),nsavf(4),n1(6),nbvar(30),nbas1,nbast,nbvar1,
     *nblock,numvar,nrun,nitscf,nitsc1,nitsc2,nconv,nalarm,no,nstep,
     *iiii,nsht,nosht,nflag3,nbasm,mvar,iex,newbas,nit,nbpos(8)
      common/blkin/bias,dgath,scfat,dmnxp,rmned,rmxed,con1,con2,
     *nvar,ndiag,nxtrp,mxxtrp,norbit,nprint,mxvar,nextra,ifcont,irun,
     *nflag1,nflag2,nsvd,lim
c
      nalar1 = nalarm - 1
      go to (20,30,40,50,60,70,80) , nalar1
 20   write (iwr,6010)
      write (iwr,6020) thre(no) , no , dgath
      go to 90
 30   write (iwr,6030)
      write (iwr,6040) thrscf , scfat
      go to 90
 40   write (iwr,6050)
      write (iwr,6060) nitscf , mxxtrp
      go to 90
 50   write (iwr,6070) red , rmned
      go to 90
 60   write (iwr,6080) red , rmxed
      go to 90
 70   if (ifcont.eq.1) then
         write (iwr,6090) term , nbvar1 , iiii
         return
      else
         write (iwr,6100) expa(nbvar1) , nbvar1 , expmn
         go to 90
      end if
 80   if (ifcont.eq.1) then
         write (iwr,6110) term , nbvar1 , iiii , dmnxp
      else
         write (iwr,6120) nbvar1 , iiii , expdif , dmnxp
         return
      end if
 90   return
 6010 format ('0',5x,'diagonalization divergence')
 6020 format (6x,'threshhold',e9.2,1x,'for vector ',i2,1x,
     +        'greater than dgath ',e9.2)
 6030 format ('0',5x,'self consistence divergence')
 6040 format (6x,'scf threshold',e9.2,' greater than scfat',e9.2)
 6050 format ('0','     too many scf iterations ')
 6060 format (6x,'iterations number',i2,'  equal mxxtrp ',i2)
 6070 format (6x,'relative energy difference ',e9.2,
     +        ' smaller than rmned ',e9.2/)
 6080 format (6x,'relative energy difference ',e9.2,
     +        ' greater than rmxed ',e9.2/)
 6090 format (6x,'contraction coefficient ',f9.5,2x,'position(',2i2,
     +        ') greater than one ')
 6100 format (6x,'exponent zeta',f9.5,' for orbital ',i2,
     +        ' smaller than expmn ',f9.5)
 6110 format (6x,'contraction coefficient ',f9.5,2x,'position(',2i2,
     +        ') smaller than dmnxp',f9.5)
 6120 format (6x,'exponent redundancy for orbitals',i2,' and ',i2,
     +        '  zetdif ',f9.5,' smaller than dmnxp ',f9.5)
c
      end
**==anagri.f
c
c grid analysis - for reqular (2 or 3d) grids only
c
      subroutine anagri(data, nx, ny, nz, iwr)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      logical omax
      dimension zmax(2),data(*)
      common/junk/gmax(5),imax(5),jmax(5),kmax(5)
      data zmax/'maxima','minima'/
      data mm1,m0,m1,m4/-1,0,1,4/
c      write (iwr,6010)
      if(nz.eq.0)then
         call avegri(data,nx*ny,iwr)
      else
         call avegri(data,nx*ny*nz,iwr)
      endif
      omax = .true.
 20   kk = 0
      ktot = 0
      if(omax)then
         write (iwr,6020) zmax(1)
      else
         write (iwr,6020) zmax(2)
      endif
      if(nz.eq.0)then
         write (iwr,6030)
         do 50 j = 2 , ny - 1
            j4 = (j-1)*nx
            k4 = j4 + nx
            l4 = j4 - nx
            do 40 i = 2 , nx - 1
               ii = i + m1
               nt1 = i + mm1
               pow1 = data(j4+i)
               if (omax) then
                 if (pow1.lt.data(j4+ii) .or. pow1.lt.data(j4+nt1) .or.
     +           pow1.lt.data(k4+i) .or. pow1.lt.data(l4+i)) go to 40
            else if (pow1.gt.data(j4+ii) .or. pow1.gt.data(j4+nt1) .or.
     +              pow1.gt.data(k4+i) .or. pow1.gt.data(l4+i)) then
                  go to 40
               end if
               kk = kk + m1
               ktot = ktot + m1
               imax(kk) = i
               jmax(kk) = j
               gmax(kk) = pow1
               if (kk.eq.m4) then
                  write (iwr,6040)
     &                   (imax(kk),jmax(kk),gmax(kk),kk=1,4)
                  kk = 0
               end if
 40         continue
 50      continue
         if (kk.ne.m0) write (iwr,6040)
     &        (imax(i),jmax(i),gmax(i),i=1,kk)
      else
         write (iwr,6035)
         n0 = 1
         n1 = n0*nx
         n2 = n1*ny
         do 80 k = 2 , nz - 1
           do 70 j = 2 , ny - 1
              do 60 i = 2 , nx - 1
               pow1 = data(1 + (i-1)*n0 + (j-1)*n1 + (k-1)*n2)
               if (omax) then
                  if (
     &             pow1.lt.data(1 +  i   *n0 + (j-1)*n1 + (k-1)*n2).or.
     &             pow1.lt.data(1 + (i-2)*n0 + (j-1)*n1 + (k-1)*n2).or.
     &             pow1.lt.data(1 + (i-1)*n0 +  j   *n1 + (k-1)*n2).or.
     &             pow1.lt.data(1 + (i-1)*n0 + (j-2)*n1 + (k-1)*n2).or.
     &             pow1.lt.data(1 + (i-1)*n0 + (j-1)*n1 +  k   *n2).or.
     &             pow1.lt.data(1 + (i-1)*n0 + (j-1)*n1 + (k-2)*n2))
     &            goto 60

               else
                  if (
     &             pow1.gt.data(1 +  i   *n0 + (j-1)*n1 + (k-1)*n2).or.
     &             pow1.gt.data(1 + (i-2)*n0 + (j-1)*n1 + (k-1)*n2).or.
     &             pow1.gt.data(1 + (i-1)*n0 +  j   *n1 + (k-1)*n2).or.
     &             pow1.gt.data(1 + (i-1)*n0 + (j-2)*n1 + (k-1)*n2).or.
     &             pow1.gt.data(1 + (i-1)*n0 + (j-1)*n1 +  k   *n2).or.
     &             pow1.gt.data(1 + (i-1)*n0 + (j-1)*n1 + (k-2)*n2))
     &                 goto 60
               end if
               kk = kk + m1
               ktot = ktot + m1
               imax(kk) = i
               jmax(kk) = j
               kmax(kk) = k
               gmax(kk) = pow1
               if (kk.eq.4) then
                  write (iwr,6060)
     &                 (imax(kk),jmax(kk),kmax(kk),gmax(kk),kk=1,4)
                  kk = 0
               end if
 60         continue
 70      continue
 80      continue
         if (kk.ne.m0)write (iwr,6060)
     &        (imax(i),jmax(i),kmax(i),gmax(i),i=1,kk)
      endif
      if ((ktot.eq.m0).and.omax) write (iwr,6050)zmax(1)
      if ((ktot.eq.m0).and..not.omax) write (iwr,6050)zmax(2)

      if(.not.omax)return
      omax = .false.
      go to 20
 6020 format (//20x,'grid ',a8/20x,11('-')//)
 6030 format (4x,4(3x,'i',4x,'j',7x,'value',6x)/)
 6035 format (4x,4(3x,'i',4x,'j',4x,'k',4x,'value',7x)/)
 6040 format (4x,4(2i5,f12.3,5x))
 6060 format (4x,4(3i5,f12.3,5x))
 6050 format (/8x,'no ',a7,'found')
      end
**==avegri.f
c
c crude grid analysis
c
      subroutine avegri(data, np, iwr)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension data(*)
      write (iwr,6010)
      bmin=+1.0d20
      bmax=-1.0d20
      sum=0.0d0
      if(np.eq.0)then
        write(iwr,*)'grid has no points'
        return
      endif
      do 10 i=1,np
         bmin=dmin1(bmin,data(i))
         bmax=dmax1(bmax,data(i))
         sum = sum+data(i)
 10   continue
      sum = sum/dble(np)
      write (iwr,6020) np, bmin, bmax, sum
 6010 format (//38x,'grid data analysis'/38x,15('-'))
 6020 format (/1x,'grid has ',i6,', points, min =',
     &     f18.8,', max =',f18.8,', average =',f18.8/)
      end
**==analm.f
      subroutine analm(q,lword)
c
c...    mulliken-population analysis (atmol-style)
c...    see atmol mulliken analysis program
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
      common/blkcore/corev(512),charge(10)
      dimension q(lword),o1e(6)
c
      logical osmall
      osmall = .false.
c
      write (iwr,6010)
      l2 = (num*(num+1))/2
      l3 = num*num
      i10 = 1
      i20 = i10 + l3
      i30 = i20 + l2
      i40 = i30 + l3
      i50 = i40 + max(l3,nat*(nat+1)/2)
      i60 = i50 + max(l3,nat)
      if (i60+l3.gt.lword) then
c....   revert to old ways / disable mull unpair
        i50 = i40 + l3
        i60 = i50 + l3
        if (i60+l3.gt.lword)
     1  call caserr('memory overfow in analm')
        osmall = .true.
      end if
c
      do loop = 1,6
       o1e(loop) = .false.
      enddo
      o1e(1) = .true.
      call getmat(q(i20),q(i10),q(i10),q(i10),q(i10),q(i10),
     +            charge,num,o1e,ionsec)
      call mullr(q(i10),q(i20),q(i30),q(i40),q(i50),q(i60),osmall)
      top = cpulft(1)
      write (iwr,6020) top
      return
 6010 format (/40x,26('-')/40x,'population analysis module'/40x,26('-')
     +        /)
 6020 format (//1x,104('-')//1x,'end of mulliken analysis at ',f8.2,
     +        ' secs')
      end
**==anuc2.f
      subroutine anuc2(noc,vlist,c,dint,nocmx)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension dint(7),vlist(nocmx,4),c(3)
      common/miscop/qnuc(6)
c...
c...compute nuclear contributions to second moments.
c...xx,yy,zz,xy,xz,yz,rr.
c...
      do 20 i = 1 , 6
         qnuc(i) = 0.0d0
 20   continue
c...
      do 50 i = 1 , noc
         mu = 0
         do 40 j = 1 , 3
            do 30 k = 1 , j
               mu = mu + 1
               qnuc(mu) = qnuc(mu) + vlist(i,4)*(vlist(i,j)-c(j))
     +                    *(vlist(i,k)-c(k))
 30         continue
 40      continue
 50   continue
c...
      dint(1) = qnuc(1)
      dint(2) = qnuc(3)
      dint(3) = qnuc(6)
      dint(4) = qnuc(2)
      dint(5) = qnuc(4)
      dint(6) = qnuc(5)
      dint(7) = dint(1) + dint(2) + dint(3)
c...
      return
      end
**==anuc3.f
      subroutine anuc3(noc,vlist,c,dint,nocmx)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      dimension dint(10),vlist(nocmx,4),c(3)
      common/miscop/cnuc(10)
c...
c...nuclear contributions to third moments.
c...xxx,xyy,xzz,yxx,yyy,yzz,zxx,zyy,zzz,xyz.
c...
      call vclr(cnuc,1,10)
c...
      do 50 i = 1 , noc
         mu = 0
         do 40 j = 1 , 3
            do 30 k = 1 , j
               do 20 m = 1 , k
                  mu = mu + 1
                  cnuc(mu) = cnuc(mu) + vlist(i,4)*(vlist(i,j)-c(j))
     +                       *(vlist(i,k)-c(k))*(vlist(i,m)-c(m))
 20            continue
 30         continue
 40      continue
 50   continue
c...
      dint(1) = cnuc(1)
      dint(2) = cnuc(3)
      dint(3) = cnuc(8)
      dint(4) = cnuc(2)
      dint(5) = cnuc(4)
      dint(6) = cnuc(9)
      dint(7) = cnuc(5)
      dint(8) = cnuc(7)
      dint(9) = cnuc(10)
      dint(10) = cnuc(6)
c...
      return
      end
**==anuc4.f
      subroutine anuc4(non,vlist,c,dint,nocmx)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension dint(10),vlist(nocmx,4),c(3)
      common/miscop/t(15),a(3)
c...
c...compute nuclear contributions to even fourth moments.
c...xxxx,yyyy,zzzz,xxyy,xxzz,yyzz,xxrr,yyrr,zzrr,rrrr.
c...
      call vclr(t,1,15)
c...
      do 60 nc = 1 , non
         a(1) = vlist(nc,1)
         a(2) = vlist(nc,2)
         a(3) = vlist(nc,3)
         chg = vlist(nc,4)
         mu = 0
         do 50 i = 1 , 3
            do 40 j = 1 , i
               do 30 k = 1 , j
                  do 20 l = 1 , k
                     mu = mu + 1
                     t(mu) = t(mu) + chg*(a(i)-c(i))*(a(j)-c(j))
     +                       *(a(k)-c(k))*(a(l)-c(l))
 20               continue
 30            continue
 40         continue
 50      continue
 60   continue
c...
      dint(1) = t(1)
      dint(2) = t(5)
      dint(3) = t(15)
      dint(4) = t(3)
      dint(5) = t(10)
      dint(6) = t(12)
      dint(7) = t(1) + t(3) + t(10)
      dint(8) = t(3) + t(5) + t(12)
      dint(9) = t(10) + t(12) + t(15)
      dint(10) = dint(7) + dint(8) + dint(9)
c...
      return
      end
**==assia.f
      subroutine assia(glab,iwr)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *(*) glab
      character *24 zbflab
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
      common/bufb/zbflab(maxorb)
      dimension glab(*)
      write (iwr,6010)
      n = 0
      do 40 ii = 1 , nshell
         iat = katom(ii)
         mini = kmin(ii)
         maxi = kmax(ii)
         if (iat.ge.100) then
            do 20 i = mini , maxi
               n = n + 1
               glab(n) = ' '
               glab(n)(1:6) = zbflab(n)(1:6)
 20         continue
         else
            do 30 i = mini , maxi
               n = n + 1
               glab(n) = ' '
               glab(n)(1:5) = zbflab(n)(1:5)
 30         continue
         end if
 40   continue
      return
 6010 format (/40x,13('*')/40x,'atom analysis'/40x,13('*')/)
      end
**==assib.f
c****f* analysis/assib
c
c     NAME
c
c       assib - assign a basis function type to atomic orbitals
c
c     FUNCTION
c
c       This routine goes through all the shells and assigns a
c       basis function type for atomic orbitals in the shell. These
c       types correspond to the appropriate angular momentum. So
c       d-, f- and g-shells should get assigned multiple types
c       corresponding to the contaminants and the true angular momentum
c       functions. In practice this separation of contaminants is not
c       used for the f- and g-shells as problems arose identifying the
c       correct transformations to set in subroutine mulltr. So far
c       all attempts led to electrons in all angular momenta (even if
c       the wavefunction clearly occupied only one), and the number of
c       electrons was not conserved.
c       Finally a table it printed detailing the classification.
c
c     SOURCE
c
      subroutine assib(glab,iwr)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *(*) glab
      character *24 zbflab
      character *2 type(5)
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
      common/bufb/zbflab(maxorb)
      dimension glab(*)
      data type/' s',' p',' d',' f',' g'/
      write (iwr,6010)
      do i = 1 , maxorb
         glab(i) = ' '
      enddo ! i
      n = 1
      do i = 1 , nshell
         j = kmax(i) - kmin(i)
         iat = katom(i)
         ind1 = 5
         if (iat.gt.100) ind1 = 6
         if (j.le.1) then
c
c           s-shell
c
            glab(n)(1:ind1) = zbflab(n)(1:ind1)
            glab(n)(ind1+1:ind1+2) = type(1)
            n = n + 1
         else if (j.lt.3) then
c
c           p-shell
c
            do ii = 1 , 3
              glab(n)(1:ind1) = zbflab(n)(1:ind1)
              glab(n)(ind1+1:ind1+2) = type(2)
              n = n + 1
           enddo ! ii
         else if (j.eq.3) then
c
c           sp-shell
c
            glab(n)(1:ind1) = zbflab(n)(1:ind1)
            glab(n)(ind1+1:ind1+2) = type(1)
            n = n + 1
            do ii = 1 , 3
              glab(n)(1:ind1) = zbflab(n)(1:ind1)
              glab(n)(ind1+1:ind1+2) = type(2)
              n = n + 1
           enddo ! ii
         else if (j.eq.5) then
c
c           d-shell
c
            glab(n)(1:ind1) = zbflab(n)(1:ind1)
            glab(n)(ind1+1:ind1+2) = type(1)
            n = n + 1
            do ii = 1 , 5
               glab(n)(ind1+1:ind1+2) = type(3)
               glab(n)(1:ind1) = zbflab(n)(1:ind1)
               n = n + 1
            enddo ! ii
         else if (j.eq.9) then
c
c           f-shell
c
            do ii = 1 , 3
               glab(n)(1:ind1) = zbflab(n)(1:ind1)
               glab(n)(ind1+1:ind1+2) = type(2)
               n = n + 1
            enddo ! ii
            do ii = 1 , 7
               glab(n)(ind1+1:ind1+2) = type(4)
               glab(n)(1:ind1) = zbflab(n)(1:ind1)
               n = n + 1
            enddo ! ii
         else if (j.eq.14) then
c
c           g-shell
c
            glab(n)(1:ind1) = zbflab(n)(1:ind1)
            glab(n)(ind1+1:ind1+2) = type(1)
            n = n + 1
            do ii = 1 , 5
               glab(n)(1:ind1) = zbflab(n)(1:ind1)
               glab(n)(ind1+1:ind1+2) = type(3)
               n = n + 1
            enddo ! ii
            do ii = 1 , 9
               glab(n)(ind1+1:ind1+2) = type(5)
               glab(n)(1:ind1) = zbflab(n)(1:ind1)
               n = n + 1
            enddo ! ii
         end if
      enddo ! i
      return
 6010 format (/40x,26('*')/40x,'s,p,d,f,g orbital analysis'/40x,
     +             26('*')/)
      end
c******
**==atlp.f
      subroutine atlp(grid, rlg1, czeta, ctran, pop, scr, scrat, rvec,
     &     mconf, icen, popin, nq, iprin)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c calculate atomic density over a grid of points
c at present, only allow user setting of ao
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
      dimension grid(*), rlg1(*), rvec(*), czeta(*), ctran(*)
      dimension pop(*), popin(*), icen(*) , scr(*), scrat(*)
c
c
      common/blkin/bias(8),nvar(10),nflag1,nope(5),
     *ajm(330),
     *ncs(90),nos(90),npopul(90),aapop(maxorb),popa(maxorb)
c
      common/miscop/p(maxat),q(maxat),r(maxat),chg(maxat),
     *cx(maxorb),cy(maxorb),cz(maxorb),popocc(maxorb)
c
      common/craypk/nterm(maxorb),icentr(maxorb),itype(maxorb),
     *mx(maxorb),my(maxorb),mz(maxorb),nx(maxorb),ny(maxorb),nz(maxorb),
     *nbasin,norb,lenbas,ngrps,non,newbas,ncol,jcol,ngrp(3),
     *iblkg,iblka,iblkv
c
c for zaname
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
      common/gjs/noc(maxorb)
c
      dimension zconca(15),nosh(4),ncsh(4),zconda(15),nsymg(3)

c
      data mm1,m0,m1,m30/
     *-1,0,1,30/
      data nsymg/4,10,0/
      data thresh,dzero,done,two/
     *1.0d-8,0.0d0,1.0d0,2.0d0/
      data zconda/
     *'1s','2s','3s','4s',
     *'2px','2py','2pz','3px','3py','3pz',
     *'3dxy','3dxz','3dyz','3dx2-y2','3dz2'/
c error flag
      oatom = .false.
      call getnpt(npt,npx,npy,npz)
      do 20 i = 1 , maxorb
         popa(i) = dzero
 20   continue
c
c set print flag
c
      nvar(6) = iprin
      write(iwr,*)' input occs', mconf
      if(mconf.ne.0)then
         do i = 1,mconf
            write(iwr,*)icen(i),popin(i)
         enddo
      endif
c
c ... default ao occupancies
c
      jmo = m0
      do 80 i = 1 , non
         izn = nint(chg(i))
         if (izn.le.m30) then
            j = m1
 30         it = itype(j)
            if (j.ge.norb) then
               do 40 j = 1 , nsym
                  it = (j+mm1)*m30 + izn
                  ncsh(j) = ncs(it)
                  nosh(j) = nos(it)
 40            continue
               do 70 j = 1 , nsym
                  nj = ncsh(j)
                  nsh = nj + nosh(j)
                  nso = j + j + mm1
                  if (nsh.ne.m0) then
                     dum = j
                     dum = dum + dum - done
                     do 60 kk = 1 , nsh
                        if (kk.gt.nj) then
                           dumx = npopul(izn+(j+mm1)*m30)/dum
                        else
                           dumx = two
                        end if
                        do 50 jj = 1 , nso
                           jmo = jmo + m1
                           popa(jmo) = dumx
 50                     continue
 60                  continue
                  end if
 70            continue
            else
               if (i.eq.icentr(j)) then
                  nsym = it
               end if
               j = j + ngrp(it)
               go to 30
            end if
         else
            oatom = .true.
            go to 90
         end if
 80   continue
      jcol = jmo
c ...
c ... config
c ...
 90   if (mconf.ne.0) then
         kmo = 0
         do 160 ico = 1 , mconf
            kcen = icen(ico)
            if (kcen.eq.0)call caserr
     &           ('unrecognised tag present in config data')
            jmo = m0
            do 150 m = 1 , non
               izn = nint(chg(m))
               j = m1
 100           it = itype(j)
               if (j.ge.norb) then
                  do 110 j = 1 , nsym
                     it = (j+mm1)*m30 + izn
                     ncsh(j) = ncs(it)
                     nosh(j) = nos(it)
 110              continue
                  do 140 j = 1 , nsym
                     nsh = ncsh(j) + nosh(j)
                     nso = j + j + mm1
                     if (nsh.ne.m0) then
                        do 130 lll = 1 , nsh
                           do 120 ll = 1 , nso
                              jmo = jmo + m1
                              if (kcen.eq.m) then
                                 kmo = kmo + 1
                                 popa(jmo) = popin(kmo)
                              end if
 120                       continue
 130                    continue
                     end if
 140              continue
               else
                  if (m.eq.icentr(j)) then
                     nsym = it
                  end if
                  j = j + ngrp(it)
                  go to 100
               end if
 150        continue
 160     continue
      end if
c
      nbsq = ncol*norb
c
      write (iwr,6060)
c .. output atomic configurations
      write (iwr,6070)
      if (oatom)
     & call caserr('invalid atom detected in atom-difference analysis')
c
      jmo = 0
      do 240 i = 1 , non
         izn = nint(chg(i))
         j = m1
 190     it = itype(j)
         if (j.ge.norb) then
            do 200 j = 1 , nsym
               it = (j+mm1)*m30 + izn
               ncsh(j) = ncs(it)
               nosh(j) = nos(it)
 200        continue
            k = jmo + m1
            kk = m0
            l = m0
            do 230 j = 1 , nsym
               kkk = m0
               nsh = ncsh(j) + nosh(j)
               nso = j + j + mm1
               if (nsh.ne.m0) then
                  do 220 lll = 1 , nsh
                     do 210 ll = 1 , nso
                        kkk = kkk + m1
                        jmo = jmo + m1
                        kk = kk + m1
                        zconca(kk) = zconda(kkk+l)
 210                 continue
 220              continue
               end if
               l = l + nsymg(j)
 230        continue
            write (iwr,6080) zaname(i) , (zconca(j),j=1,kk)
            write (iwr,6090) (popa(j),j=k,jmo)
         else
            if (i.eq.icentr(j)) then
               nsym = it
            end if
            j = j + ngrp(it)
            go to 190
         end if
 240  continue
c
c perform molecular calculation
c
      do 305 i = 1 , ncol
         noc(i)=0
         if (dabs(pop(i)).gt.thresh) noc(i) = 1
 305  continue
c
      i0 = 1
      i1 = i0 + ncol
c flag for total density
      imo = 0
c
      do 310 i = 1,npt
         call getpt(i,ax,ay,az,grid,1,idat)
         call dplota(rvec, ctran, czeta, pop, scr(i0), scr(i1),
     &        ax, ay, az, rlg1(idat),imo)
 310  continue
c
c  scratch file for atomic calculations
      iblka = ibl7la
c
c generate atomic wavefn in rvec
c
      call atom(scrat,rvec,ctran,czeta,nq)
c
c perform molecular calculation
c
      do 405 i = 1 , ncol
         noc(i)=0
         if (dabs(popa(i)).gt.thresh) noc(i) = 1
 405  continue
c
      do 410 i = 1,npt
         call getpt(i,ax,ay,az,grid,1,idat)
         call dplota(rvec, ctran, czeta, popa, scr(i0), scr(i1),
     &        ax, ay, az, dum, imo)
         rlg1(idat) = rlg1(idat) - dum
 410  continue
c
c restore molecular vectors
c
      call rdedx(rvec,nbsq,iblkv,num8)
      return
*6030 format (/' output intermediate atom-scf results')
 6060 format (/' **** generate atom-difference density grid')
 6070 format (/5x,'atomic configurations'/5x,21('-'))
 6080 format (/' atom  ',a8//' orbital :',10x,15a7)
 6090 format (/' occupancy :',8x,15(f6.4,1x))
*6170 format (/40x,'orbital occupations  (from vectors section)'/40x,
*    +        43('-'))
*6180 format (/40x,'orbital populations (from cards)'/40x,32('-')/)
*6190 format (/4x,9f14.7)
      end
**==atnear.f
      function atnear(x, y, z)
c
c DESCRIPTION : Function which returns the distance to the closest
c               atom centre from a given point.
c
c
c PARAMETERS  : NAME               DESCRIPTION
c               ----               -----------
c               x,
c               y,                 Coordinates of point from which
c               z                  distances are calculated.
c
c
c AUTHOR      : M. A. FOX.
c
c Date  :  26/3/90
c
      implicit real*8  (a-h, o-z)
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
c find distance to first atom.
c
      small = discal(x, y, z, c(1, 1), c(2, 1), c(3, 1))
c
c compare first distance found against distances to all other atoms.
c
      do 10 i = 2, nat
       dist = discal(x, y, z, c(1, i), c(2, i), c(3, i))
       if (dist .lt. small) small = dist
   10 continue
      atnear = small
      return
      end
**==atom.f
      subroutine atom(qatom,qq,ctran,expa,nq)
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
      dimension pq(3),qatom(*),qq(*),ilif(3),ctran(*),expa(*)
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
      common/miscop/p(maxat),q(maxat),r(maxat),chg(maxat),
     *cx(maxorb),cy(maxorb),cz(maxorb),occ(maxorb)
      common/craypk/nterm(maxorb),icentr(maxorb),itype(maxorb),
     *mx(maxorb),my(maxorb),mz(maxorb),nx(maxorb),ny(maxorb),nz(maxorb),
     *nbasin,norb,lenbas,ngrps,non,newbas,ncol,jcol,ngrp(3),iblkg,iblka
      common/junk/pqn(99),expat(99),ctat(mxprms,30),ajmn(12),czn,
     *expmn,c(15,15),itran(mxprms,30),nrow(30),nbcon(4),nbas(4),nsym,
     *nosh(4),ncsh(4),nccup(4),ibase(33),
     *corb(10,9),fac(20)
      common/blkin/bias(8),nvar(10),nflag1,nope(5),
     *ajm(330),
     *ncs(90),nos(90),npopul(90),pop(maxorb),popa(maxorb)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      data pq/1.0d0,2.0d0,3.0d0/
      data m0,m1,ms1,m4,m30,m225/
     *0,1,-1,4,30,225/
      data done/1.0d0/
      data ilif/0,1,4/
      data pi/3.14159265359d0/
      iblkb = iblka
      do 110 i = 1 , non
         fac(1) = done
         bot = 4.0d0*pi
         top = done
         ss = done
         do 20 ifac = 2 , 20
            top = top + done
            ss = ss*top
            fac(ifac) = ss*bot
 20      continue
         icen = m0
         kk = m0
         czn = chg(i)
         iczn = nint(czn)
         write (iwr,6010) zaname(i)
         do 30 j = 1 , m4
            nbas(j) = m0
            nbcon(j) = m0
 30      continue
         j = m1
 40      if (j.gt.norb) then
            k = nsym
            do 50 j = 1 , nsym
               it = (j-1)*30 + iczn
               ncsh(j) = ncs(it)
               nosh(j) = nos(it)
               if ((nosh(j)+ncsh(j)).eq.m0) k = k + ms1
               nccup(j) = npopul(it)
 50         continue
            nsym = k
            iczn = (iczn-1)*11
            do 60 j = 1 , 11
               ajmn(j) = ajm(iczn+j)
 60         continue
            if (icen.gt.m30 .or. kk.gt.99) go to 230
            call atomss(qatom(1),nq,iwr)
            call wrt3(c,m225,iblka,num8)
            iblka = iblka + m1
         else
            it = itype(j)
            pqni = pq(it)
            if (i.eq.icentr(j)) then
               icen = icen + m1
               nsym = it
               ibase(icen) = j
               n1 = nterm(j)
               jj = (j+ms1)*mxprms
               do 100 k = 1 , n1
                  kk = kk + m1
                  itran(k,icen) = kk
                  mbase = jj + k
                  dum = expa(mbase)
                  dum2 = dum + dum
                  cc = (pi/dum2)**1.5d0
                  go to (90,70,80) , it
 70               cc = cc/(dum2+dum2)
                  go to 90
 80               cc = cc/(4.0d0*dum2*dum2)
 90               cc = ctran(mbase)*dsqrt(cc)
                  ctat(k,icen) = cc
                  expat(kk) = dum
                  pqn(kk) = pqni
 100           continue
               nrow(icen) = n1
               nbas(it) = nbas(it) + n1
               nbcon(it) = nbcon(it) + m1
            end if
            j = j + ngrp(it)
            go to 40
         end if
 110  continue
      inc = m0
      i = jcol*norb
      write (iwr,6020)
c
c      if ((ng8+i).gt.lword) go to 240
      call vclr(qq,1,i)
      call vclr(corb,1,90)
      if (nflag1.eq.m1) then
         do 120 i = 1 , 7
            corb(i,i) = done
 120     continue
         corb(8,8) = 0.8660254d0
         corb(9,8) = -0.8660254d0
         corb(8,9) = -0.5d0
         corb(9,9) = -0.5d0
         corb(10,9) = done
      else
         do 130 i = 1 , 4
            corb(i,i) = done
 130     continue
         corb(5,8) = done
         corb(6,5) = done
         corb(7,6) = done
         corb(8,7) = done
         corb(9,9) = done
      end if
      jmo = m0
      do 220 i = 1 , non
         iczn = nint(chg(i))
         do 140 j = 1 , m4
            nbcon(j) = m0
 140     continue
         j = m1
         icen = m0
 150     it = itype(j)
         if (j.ge.norb) then
            k = nsym
            do 160 j = 1 , nsym
               it = (j+ms1)*m30 + iczn
               ncsh(j) = ncs(it)
               nosh(j) = nos(it)
               if ((nosh(j)+ncsh(j)).eq.m0) k = k + ms1
 160        continue
            nsym = k
            call rdedx(c,m225,iblkb,num8)
            nstep = m0
            kkk = m0
            do 210 j = 1 , nsym
               mbase = ilif(j)
               nsymj = ngrp(j)
               nj = ncsh(j)
               nsh = nj + nosh(j)
               nbas1 = nbcon(j)
               nso = j + j + ms1
               do 200 kk = 1 , nsh
                  k = nstep + kk
                  do 190 jj = 1 , nso
                     mb = mbase + jj
                     jmo = jmo + m1
                     do 180 l = 1 , nbas1
                        iba = ibase(l+kkk)
                        do 170 lll = 1 , nsymj
                           kstep = iba + lll + ms1
                           qq(inc+kstep) = c(l,k)*corb(mbase+lll,mb)
 170                    continue
 180                 continue
                     inc = inc + norb
 190              continue
 200           continue
               nstep = nstep + nsh
               kkk = kkk + nbas1
 210        continue
            iblkb = iblkb + m1
         else
            if (icentr(j).eq.i) then
               nsym = it
               icen = icen + m1
               ibase(icen) = j
               nbcon(it) = nbcon(it) + m1
            end if
            j = j + ngrp(it)
            go to 150
         end if
 220  continue
      return
 230  call caserr('invalid atom in atom-difference graphical analysis')
c 240  call caserr('insufficient memory for atom-difference analysis')
      return
 6010 format (/1x,104('-')//40x,'scf calculation on atom ',a8/,40x,
     +        32('-')/)
 6020 format (/1x,104('-'))
      end
**==atomco.f
      subroutine atomco
c     one configuration
c     one electron problem
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
      common/junk/pqn(99),expa(99),ctran(mxprms,30),ajmn(12),czn,
     *expmn,c(15,15),itran(mxprms,30),nrow(30),nbcon(4),nbas(4),
     *nsym,nosh(4),ncsh(4),nccup(4),ibase(33),
     *popul(4),facto(30),fatt(30),thre(15),a(15,15),eivr(15,15),
     *cm2scf(15,15),cm1scf(15,15),eps(15),gmin(15,15),smin(15,15),
     *populs(4),term,energ,pot,cin,vir,thrscf,enorm,
     *red,expdif,encor,
     *nsize(10),nsavf(4),n1(6),nbvar(30),nbas1,nbast,nbvar1,
     *nblock,numvar,nrun,nitscf,nitsc1,nitsc2,nconv,nalarm,no,nstep,
     *iiii,nsht,nosht,nflag3,nbasm,mvar,iex,newbas,nit,nbpos(8)
      common/blkin/bias,dgath,scfat,dmnxp,rmned,rmxed,con1,con2,
     *nvar,ndiag,nxtrp,mxxtrp,norbit,nprint,mxvar,nextra,ifcont,irun,
     *nflag1,nflag2,nsvd,lim
      data dzero,done,two,half/
     *0.0d0,1.0d0,2.0d0,0.5d0/
c
c     normalization of contracted coefficients  (slater or gaussian)
c
      do 50 k = 1 , newbas
         sum = dzero
         n = nrow(k)
         do 30 i = 1 , n
            ctr = ctran(i,k)
            i1 = itran(i,k)
            np = pqn(i1)
            dump = expa(i1)
            do 20 j = 1 , i
               if (j.lt.i) then
                  j1 = itran(j,k)
                  dumq = expa(j1)
                  nq = pqn(j1)
                  npq = np + nq + 1
                  dumpq = half*(dump+dumq)
                  if (nflag1.ne.0) then
                     vpq = facto(npq-1)/dumpq**(half*npq)
                     vq = facto(nq+nq)/dumq**(nq+half)
                     vp = facto(np+np)/dump**(np+half)
                  else
                     nn = 2*np + 1
                     nm = 2*nq + 1
                     vp = fatt(nn)/dump**nn
                     vq = fatt(nm)/dumq**nm
                     vpq = fatt(npq)/dumpq**npq
                  end if
                  sum = sum + two*ctran(j,k)*ctr*vpq/dsqrt(vp*vq)
               else
                  sum = sum + ctr**2
                  go to 30
               end if
 20         continue
 30      continue
         sum = done/dsqrt(sum)
         do 40 i = 1 , n
            ctran(i,k) = ctran(i,k)*sum
 40      continue
 50   continue
c
      return
      end
**==atomex.f
      subroutine atomex(s)
c
c     one configuration
c     scf iterations and extrapolations
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
      dimension s(*)
      common/junk/pqn(99),expa(99),ctran(mxprms,30),ajmn(12),czn,
     *expmn,c(15,15),itran(mxprms,30),nrow(30),nbcon(4),nbas(4),
     *nsym,nosh(4),ncsh(4),nccup(4),ibase(33),
     *popul(4),facto(30),fatt(30),thre(15),a(15,15),eivr(15,15),
     *cm2scf(15,15),cm1scf(15,15),eps(15),gmin(15,15),smin(15,15),
     *populs(4),term,energ,pot,cin,vir,thrscf,enorm,
     *red,expdif,encor,
     *nsize(10),nsavf(4),n1(6),nbvar(30),nbas1,nbast,nbvar1,
     *nblock,numvar,nrun,nitscf,nitsc1,nitsc2,nconv,nalarm,no,nstep,
     *iiii,nsht,nosht,nflag3,nbasm,mvar,iex,newbas,nit,nbpos(8)
      common/blkin/bias,dgath,scfat,dmnxp,rmned,rmxed,con1,con2,
     *nvar,ndiag,nxtrp,mxxtrp,norbit,nprint,mxvar,nextra,ifcont,irun,
     *nflag1,nflag2,nsvd,lim
      data d8,dzero,two/
     *1.0d-8,0.0d0,2.0d0/
      data done/1.0d0/
c
      nbias = bias
      threm = dzero
      do 20 i = 1 , nsht
         if (thre(i).gt.threm) then
            threm = thre(i)
         end if
 20   continue
      if (nsvd.eq.1) threm = d8
      thrscf = threm*(two**nbias)
      if (thrscf.lt.scfat) then
         sup = dzero
         nstep = 0
         do 50 i = 1 , nsym
            nbas1 = nbas(i)
            nsh = ncsh(i) + nosh(i)
            do 40 j = 1 , nsh
               j1 = j + nstep
               do 30 m = 1 , nbas1
                  term = dabs(c(m,j1)-cm1scf(m,j1))
                  if (term.gt.sup) then
                     sup = term
                  end if
 30            continue
 40         continue
            nstep = nstep + nsh
 50      continue
         if (sup.le.thrscf) then
c    ********     convergence  sets  nconv  non zero
            nconv = 1
         else
            if (nitsc2.ge.nxtrp) then
               bias = bias + done
               nitsc2 = 0
            end if
            if (nitscf.lt.mxxtrp) then
               nitscf = nitscf + 1
               nitsc2 = nitsc2 + 1
               if (nitsc1.lt.2) then
                  nitsc1 = nitsc1 + 1
               else
                  nstep1 = 0
                  nstep2 = 0
                  do 100 i = 1 , nsym
                     nbas1 = nbas(i)
                     nsh = ncsh(i) + nosh(i)
                     do 90 j = 1 , nsh
                        j1 = j + nstep1
                        do 60 m = 1 , nbas1
                           sup = c(m,j1)
                           enorm = cm1scf(m,j1)
                           term = cm2scf(m,j1)
                           if (sup.eq.enorm) then
                              if (sup.eq.term) go to 60
                           end if
                           term = (sup-enorm)/(two*enorm-sup-term)
                           if (dabs(term).lt.two) then
                              c(m,j1) = sup + term*(sup-enorm)
                           end if
 60                     continue
                        enorm = dzero
                        do 80 m = 1 , nbas1
                           do 70 k = 1 , m
                              kp = nstep2 + (m*(m-1)/2) + k
                              term = c(m,j1)*s(kp)*c(k,j1)
                              if (m.ne.k) then
                                 term = term + term
                              end if
                              enorm = enorm + term
 70                        continue
 80                     continue
                        enorm = done/dsqrt(enorm)
                call dscal(nbas1,enorm,c(1,j1),1)
 90                  continue
                     nstep2 = nstep2 + n1(i)
                     nstep1 = nstep1 + nsh
 100              continue
                  nitsc1 = 1
               end if
            else
               nalarm = 4
            end if
         end if
      else
         nalarm = 3
      end if
      return
c
      end
**==atomss.f
      subroutine atomss(space,lword,iwr)
c
c   main program one configuration lcao atomic scf calculation
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension space(*)
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
      common/junk/pqn(99),expa(99),ctran(mxprms,30),ajmn(12),czn,
     *expmn,c(15,15),itran(mxprms,30),nrow(30),nbcon(4),nbas(4),
     *nsym,nosh(4),ncsh(4),nccup(4),ibase(33),
     *popul(4),facto(30),fatt(30),thre(15),a(15,15),eivr(15,15),
     *cm2scf(15,15),cm1scf(15,15),eps(15),gmin(15,15),smin(15,15),
     *populs(4),term,energ,pot,cin,vir,thrscf,enorm,
     *red,expdif,encor,
     *nsize(10),nsavf(4),n1(6),nbvar(30),nbas1,nbast,nbvar1,
     *nblock,numvar,nrun,nitscf,nitsc1,nitsc2,nconv,nalarm,no,nstep,
     *iiii,nsht,nosht,nflag3,nbasm,mvar,iex,newbas,nit,nbpos(8)
      common/blkin/bias,dgath,scfat,dmnxp,rmned,rmxed,con1,con2,
     *nvar,ndiag,nxtrp,mxxtrp,norbit,nprint,mxvar,nextra,ifcont,irun,
     *nflag1,nflag2,nsvd,lim
c
      data done,dm5/1.0d0,.1d-5/
      data m3/3/
c     nflag1=0   slater orbitals
c     nflag1=1   gauss orbitals
c      nflag2=0   hf energy  no corr.energy
c      nflag2=1   hf energy and cor.energy
c      nflag2=2   chf energy and orbitals
c     fatt   factorials
c       facto   double  factorials
      facto(1) = done
      facto(2) = done
      fatt(1) = done
      fatt(2) = done
      do 20 i = 3 , 16
         im1 = i - 1
         facto(i) = im1*facto(i-2)
         fatt(i) = im1*fatt(im1)
 20   continue
      do 30 i = 17 , 30
         im1 = i - 1
         facto(i) = im1*facto(i-2)
         fatt(i) = im1*fatt(im1)
 30   continue
      do 100 numb = 1 , lim
c       numb  counts  the number of complete  scf  calcs
         nalarm = 0
c       p  and  q  matrices  are not initialised  to zero
c      the   c   matrix  is initialised  in   input
         call jnput(iwr)
         mone = 0
         msq = 0
         do 40 i = 1 , nsym
            j = nbas(i)
            msq = msq + j*j
            mone = mone + j*(j+1)/2
 40      continue
         mtwo = mone*(mone+1)/2
         ip = 1
         iq = ip + mtwo
         is = iq + mtwo
         it = is + mone
         iu = it + mone
         ifo = iu + mone
         ifc = ifo + mone
         idt = ifc + mone
         idos = idt + msq
         ipcap = idos + msq
         iqcap = ipcap + mone
         itot = iqcap + mone
         call corchk(lword,itot,m3,iwr)
c      iex,numvar,mvar   used  in exponent  variation
         numvar = 0
         mvar = 0
         iex = 0
         nrun = 1
         call atomco
         call coeisg(space(is),space(it),space(iu))
         call cteisg(space(ip),space(iq))
         do 50 i = 1 , nsht
            thre(i) = dm5
 50      continue
         nitscf = 1
         nitsc1 = 1
         nitsc2 = 1
         nconv = 0
 60      call densi(space(idt),space(idos))
         if (nconv.eq.0) then
            call hamil(space(ip),space(iq),space(idt),space(idos),
     +                 space(ipcap),space(iqcap),space(is),space(it),
     +                 space(iu),space(ifc),space(ifo))
c
c     save the results of the preceding scf runs
c
            nstep2 = 0
            do 90 i = 1 , nsym
               nbas1 = nbas(i)
               nsh = ncsh(i) + nosh(i)
               do 80 j = 1 , nsh
                  j1 = j + nstep2
                  do 70 m = 1 , nbas1
                     cm2scf(m,j1) = cm1scf(m,j1)
                     cm1scf(m,j1) = c(m,j1)
 70               continue
 80            continue
               nstep2 = nstep2 + nsh
 90         continue
c
c
c**** jacobi diagonalization
            call jdaaa(space(idt),space(idos),space(is),space(ifo),
     +                 space(ifc))
c
            call atomex(space(is))
c
            if (nalarm.le.0) go to 60
            call allrm(iwr)
            call output(space(is),space(it),space(iu),space(ip),
     +                  space(iq),space(ipcap),space(iqcap),space(ifc),
     +                  space(ifo),space(idt),space(idos),iwr)
         else
c        nconv  is  set  non zero  on convergence  by  atomex
            call enerat(space(ip),space(iq),space(idt),space(idos),
     +                  space(iu),space(it))
c
            call output(space(is),space(it),space(iu),space(ip),
     +                  space(iq),space(ipcap),space(iqcap),space(ifc),
     +                  space(ifo),space(idt),space(idos),iwr)
         end if
 100  continue
c
      return
      end
**==bnuc3.f
      subroutine bnuc3(noc,vlist,c,dint,nocmx)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension dint(3),vlist(nocmx,4),c(3)
      common/miscop/cnuc(10)
c...
c...nuclear contributions to components of third moments.
c...xrr,yrr,zrr.
c...
      do 20 i = 1 , 10
         cnuc(i) = 0.0d0
 20   continue
c
      do 60 i = 1 , noc
         mu = 0
         do 50 j = 1 , 3
            do 40 k = 1 , j
               do 30 m = 1 , k
                  mu = mu + 1
                  cnuc(mu) = cnuc(mu) + vlist(i,4)*(vlist(i,j)-c(j))
     +                       *(vlist(i,k)-c(k))*(vlist(i,m)-c(m))
 30            continue
 40         continue
 50      continue
 60   continue
c...
      rrx = cnuc(1) + cnuc(3) + cnuc(8)
      rry = cnuc(2) + cnuc(4) + cnuc(9)
      rrz = cnuc(5) + cnuc(7) + cnuc(10)
      dint(1) = rrx
      dint(2) = rry
      dint(3) = rrz
c
      return
      end
**==bnuc4.f
      subroutine bnuc4(non,vlist,c,dint,nocmx)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension dint(9),vlist(nocmx,4),c(3)
      common/miscop/t(15),a(3)
c...
c...compute nuclear contributions to odd fourth moments.
c...xxxy,xyyy,xyzz,xxxz,xyyz,xzzz,xxyz,yyyz,yzzz.
c...
      call vclr(t,1,15)
c...
      do 60 nc = 1 , non
         a(1) = vlist(nc,1)
         a(2) = vlist(nc,2)
         a(3) = vlist(nc,3)
         chg = vlist(nc,4)
         mu = 0
         do 50 i = 1 , 3
            do 40 j = 1 , i
               do 30 k = 1 , j
                  do 20 l = 1 , k
                     mu = mu + 1
                     t(mu) = t(mu) + chg*(a(i)-c(i))*(a(j)-c(j))
     +                       *(a(k)-c(k))*(a(l)-c(l))
 20               continue
 30            continue
 40         continue
 50      continue
 60   continue
c...
      dint(1) = t(2)
      dint(2) = t(4)
      dint(3) = t(11)
      dint(4) = t(6)
      dint(5) = t(8)
      dint(6) = t(13)
      dint(7) = t(7)
      dint(8) = t(9)
      dint(9) = t(14)
c...
      return
      end
**==coeisg.f
      subroutine coeisg(s,t,u)
c     one configuration
c     computes the contraction matrices s,u,t, with slater or gaussian
c     orbitals as basis
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
      dimension s(*),t(*),u(*)
      common/junk/pqn(99),expa(99),ctran(mxprms,30),ajmn(12),czn,
     *expmn,c(15,15),itran(mxprms,30),nrow(30),nbcon(4),nbas(4),
     *nsym,nosh(4),ncsh(4),nccup(4),ibase(33),
     *popul(4),facto(30),fatt(30),thre(15),a(15,15),eivr(15,15),
     *cm2scf(15,15),cm1scf(15,15),eps(15),gmin(15,15),smin(15,15),
     *populs(4),term,energ,pot,cin,vir,thrscf,enorm,
     *red,expdif,encor,
     *nsize(10),nsavf(4),n1(6),nbvar(30),nbas1,nbast,nbvar1,
     *nblock,numvar,nrun,nitscf,nitsc1,nitsc2,nconv,nalarm,no,nstep,
     *iiii,nsht,nosht,nflag3,nbasm,mvar,iex,newbas,nit,nbpos(8)
      common/blkin/bias,dgath,scfat,dmnxp,rmned,rmxed,con1,con2,
     *nvar,ndiag,nxtrp,mxxtrp,norbit,nprint,mxvar,nextra,ifcont,irun,
     *nflag1,nflag2,nsvd,lim
      data dzero,half,done,two,cons/
     %0.0d0,0.5d0,1.0d0,2.0d0,1.59576912d0/
c
      nstep = 0
      kc = 1
      do 60 l = 1 , nsym
         nbas1 = nbcon(l)
         do 50 i = 1 , nbas1
            i2 = i + nstep
            kf = nrow(i2)
            do 40 j = 1 , i
               j2 = j + nstep
               lf = nrow(j2)
               sr = dzero
               ur = dzero
               tr = dzero
               do 30 k = 1 , kf
                  ctr = ctran(k,i2)
                  i1 = itran(k,i2)
                  np = pqn(i1)
                  dump = expa(i1)
                  wp = two*(np-l)/dump
                  do 20 la = 1 , lf
                     j1 = itran(la,j2)
                     dumq = expa(j1)
                     nq = pqn(j1)
                     dumpq = half*(dump+dumq)
                     npq = np + nq + 1
                     npq1 = npq - 1
                     npq2 = npq - 2
                     wq = two*(nq-l)/dumq
                     if (nflag1.ne.0) then
                        vp = facto(np+np)/dump**(np+half)
                        vq = facto(nq+nq)/dumq**(nq+half)
                        vpq = facto(np+nq)/dumpq**(half*(np+nq+1))
                        if (np+nq.le.2) then
                           vpqm2 = done
                        else
                           vpqm2 = facto(np+nq-2)
     +                             /dumpq**(half*(np+nq-1))
                        end if
                        vpq1 = facto(np+nq-1)/dumpq**(half*(np+nq))*cons
                        vpqp2 = facto(np+nq+2)/dumpq**(half*(np+nq+3))
                        term2 = vpqp2 - vpq*(wp+wq) + wp*wq*vpqm2
                     else
                        nn = 2*np + 1
                        nm = 2*nq + 1
                        vp = fatt(nn)/dump**nn
                        vq = fatt(nm)/dumq**nm
                        vpq = fatt(npq)/dumpq**npq
                        vpq1 = two*fatt(npq1)/dumpq**npq1
                        vpq2 = fatt(npq2)/dumpq**npq2
                        term2 = vpq - (wp+wq)*vpq1/two + wp*wq*vpq2
                     end if
                     term1 = done/dsqrt(vp*vq)
                     cs = ctran(la,j2)*ctr
                     ur = term1*vpq1*cs + ur
                     sr = term1*vpq*cs + sr
                     tr = half*dump*dumq*term1*term2*cs + tr
 20               continue
 30            continue
               s(kc) = sr
c  ....  the (  contracted)   nuclear attraction integrals  are in u
c  ....+ many  centre contributions  must be added here
               u(kc) = ur
               t(kc) = tr
               kc = kc + 1
 40         continue
 50      continue
         nstep = nstep + nbas1
 60   continue
c
      return
      end
**==contjk.f
      subroutine contjk(gr, index, vect, ctran, czeta, pop,
     &     scr, value, nmax, iwr)
c
c DESCRIPTION : Driver routine for wrap, loop over
c               edges of 3D elements
c
c PARAMETERS  : NAME               DESCRIPTION
c               ----               -----------
c               gr            Array for grid points on the contour surface
c               index         Number of equidensity points located.
c               vect,         Arrays used by GAMESS for calculating
c               ctran,        densities.
c               czeta
c               value         the contour height
c               nmax          the memory limit on number of pts
c
c AUTHOR      : M. A. FOX (adapted P. Sherwood)
c
c Date  :  26/3/90
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

      dimension  gr(*), ctran(*), czeta(*), vect(*), pop(*),scr(*)
      parameter (toler = 0.3d0)
      common/craypk/nbuf(9,maxorb),nbasin,norb,lenbas,ngrps,non,newbas,
     &     ncol,jcol,ngrp(3),iblkg,iblka
      common/gjs/noc(maxorb)
      data thresh/1.0d-8/
c
c - initialisation
c
      do 5 i = 1 , ncol
         noc(i)=0
         if (dabs(pop(i)).gt.thresh) noc(i) = 1
 5    continue
c
c scratch space
c
      i0 = 1
      i1 = i0 + ncol
c - loop d,e,f=1 loop moved to ngraph

      dmax = -1.0d+30
      dmin = +1.0d30

      call getnpt(n,npx,npy,npz)
      index = 0
      ipt = 0
      ozmiss = .false.
      do 10 i = 1,npz
         oymiss = .false.
         do 20 j = 1,npy
            oxmiss = .false.
            do 30 k = 1,npx
c
c check for space
c
               if(index.eq.nmax)then
                  write(iwr,*)'too many points generated in wrap'
                  return
               endif
               ipt = ipt+1
c
c get the coordinates of the current point
c
               call getpt(ipt,ax,ay,az,dum,1,idum)

               dist1 = atnear(ax, ay, az)
               if(dist1.le.toler)goto 30

               if (.not. oxmiss) then
c
c generate a point along the x plot axis
c
                  ipt2 = ipt+1
                  call getpt(ipt2,pointx,pointy,pointz,dum,1,idum)

                  dist2 = atnear(pointx, pointy, pointz)
                  if (dist2 .gt. toler)
     &                 call wrap(vect, ctran, czeta, pop,
     &                 scr(i0), scr(i1),
     &                 ax, ay, az, pointx, pointy, pointz, gr,
     &                 index, dmax, dmin, value, iwr)
               end if
               if (.not. oymiss) then
c
c generate a point along the y plot axis
c
                  ipt2 = ipt+npx
                  call getpt(ipt2,pointx,pointy,pointz,dum,1,idum)

                  dist2 = atnear(pointx, pointy, pointz)
                  if (dist2 .gt. toler)
     &                 call wrap(vect, ctran, czeta, pop,
     &                 scr(i0), scr(i1),
     &                 ax, ay, az, pointx, pointy, pointz, gr,
     &                 index, dmax, dmin, value, iwr)
               end if
               if (.not. ozmiss) then
c
c generate a point along the z plot axis
c
                  ipt2 = ipt+npx*npy
                  call getpt(ipt2,pointx,pointy,pointz,dum,1,idum)

                  dist2 = atnear(pointx, pointy, pointz)
                  if (dist2 .gt. toler)
     &                 call wrap(vect, ctran, czeta, pop,
     &                 scr(i0), scr(i1),
     &                 ax, ay, az, pointx, pointy, pointz, gr,
     &                 index, dmax, dmin, value, iwr)
               end if

               if (k .eq. npx - 1) oxmiss = .true.
 30         continue
            if (j .eq. npy - 1) oymiss = .true.
 20      continue
         if (i .eq. npz - 1) ozmiss = .true.
 10   continue
      write(iwr,*)
      write(iwr,*) 'Largest value found on grid is  ',dmax
      write(iwr,*) 'Smallest value found on grid is ',dmin
      return
      end
**==contra.f
      subroutine contra(dint,dint11,q,nocc,nbas,nt,consf)
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
      dimension dint(*),q(*),dint11(*)
      common/junk/eta(2,mxgaus),vlist(maxat,4),frocc(maxorb),
     * nfirst(4,maxorb),isky(mxgaus),nocco(maxorb)
c
      real*8 cconv, atwt, cspin
      integer lin, kin, nin, ngmx, nfmx, nomx, ncmx, ntmx, ntcol, npmx
      integer icon, ncom, maxnt,itps, maxcon
c
      common /datprp/ cconv(11),atwt(30),cspin(30),
     + lin(84),kin(84),nin(84),
     + ngmx,nfmx,nomx,ncmx,ntmx,ntcol,npmx,
     + icon(24),ncom(22),maxnt,itps,maxcon
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      integer mp, lwordp, lvec, ied3, newbas, lenb
      integer nbfnd, lenbad, lentrd, len3, ngrpm
      integer idena, nwroot, nosec
      logical noscf
      common /limy/ mp(10),lwordp,lvec,ied3,newbas,lenb,
     +              nbfnd,lenbad,lentrd,len3,ngrpm,
     +              idena,noscf,nwroot,nosec(mxroot)
c
      nocc2 = ikyp(nocc)
      do 70 mm = 1 , nt
         ii4 = mp(mm) + 2
         nu = 0
         do 60 i = 1 , nocc
            imap = ilifq(i)
            do 50 j = 1 , i
               jmap = ilifq(j)
               nu = nu + 1
               pp = 0.0d0
               if (icon(8).le.0) then
                  if ((nocco(i)+nocco(j)).lt.2 .or. i.ne.j) then
                     dint11(nu) = pp
                     go to 50
                  end if
               end if
               do 30 k = 2 , nbas
                  kk = k - 1
                  mapk = iky(k)
                  iimap = ii4 + mapk
                  qmapi = q(imap+k)
                  qmapj = q(jmap+k)
                  do 20 l = 1 , kk
                     iimap = iimap + 1
                     fac1 = qmapi*q(jmap+l)
                     fac2 = q(imap+l)*qmapj
                     pp = pp + (fac1+consf*fac2)*dint(iimap)
 20               continue
 30            continue
               do 40 k = 1 , nbas
                  mapk = ikyp(k)
                  pp = pp + q(imap+k)*q(jmap+k)*dint(ii4+mapk)
 40            continue
               dint11(nu) = pp
 50         continue
 60      continue
         call dcopy(nocc2,dint11,1,dint(ii4+1),1)
 70   continue
c
      return
      end
**==corchk.f
      subroutine corchk(lword,iparam,itype,iwr)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      dimension ztype(4)
      data ztype/'basic','grid','atomscf','plot'/
      lwor = lword - iparam
      if (lwor.lt.0) then
         lwor = -lwor
         write (iwr,6010) ztype(itype) , lwor
         call caserr('insufficient memory for requested analysis')
      end if
      write (iwr,6020) ztype(itype) , lwor
c
      return
 6010 format (/1x,'**** insufficient core assigned in ',a7,' mode'//1x,
     +        '**** additional ',i10,' words required')
 6020 format (/' main core not used in ',a7,' mode =',i6,' words')
      end
**==cteisg.f
      subroutine cteisg(p,q)
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
      dimension p(*),q(*),almn(4)
      common/junk/pqn(99),expa(99),ctran(mxprms,30),ajmn(12),czn,
     *expmn,c(15,15),itran(mxprms,30),nrow(30),nbcon(4),nbas(4),
     *nsym,nosh(4),ncsh(4),nccup(4),ibase(33),
     *popul(4),facto(30),fatt(30),thre(15),a(15,15),eivr(15,15),
     *cm2scf(15,15),cm1scf(15,15),eps(15),gmin(15,15),smin(15,15),
     *populs(4),term,energ,pot,cin,vir,thrscf,enorm,
     *red,expdif,encor,
     *nsize(10),nsavf(4),n1(6),nbvar(30),nbas1,nbast,nbvar1,
     *nblock,numvar,nrun,nitscf,nitsc1,nitsc2,nconv,nalarm,no,nstep,
     *iiii,nsht,nosht,nflag3,nbasm,mvar,iex,newbas,nit,nbpos(8)
      common/blkin/bias,dgath,scfat,dmnxp,rmned,rmxed,con1,con2,
     *nvar,ndiag,nxtrp,mxxtrp,norbit,nprint,mxvar,nextra,ifcont,irun,
     *nflag1,nflag2,nsvd,lim
      equivalence (almn(1),smin(1,1))
c
      data dzero,half,cons,two/
     %0.0d0,0.5d0,1.59576912d0,2.0d0/
      data done/1.0d0/
      data m0,m2/0,2/
c
c      teipsg computes the supermatrices p and q using a slater basis
c    set if nflag1=0 and a a gaussian basis set if nflag1=1 (contraction
c
c    this part sets up the coefficients lambda,p,q, and mu,r,s
c
      j = 0
      kmx = 0
      do 130 i = 1 , nsym
         kin = kmx + 1
         kmx = kin + nbcon(i) - 1
         do 120 k = kin , kmx
            nrowk = nrow(k)
            do 110 l = kin , k
               nrowl = nrow(l)
               mmx = 0
               do 100 im = 1 , i
                  nulo = i - im + 1
                  nuhi = i + im - 1
                  inum = 0
                  do 20 kik = nulo , nuhi , 2
                     nu = kik - nulo
                     inum = inum + 1
                     bum1 = fatt(nu+1)/fatt((nu+2)/2)**2
                     nu = kik + nulo
                     bum1 = bum1*fatt(nu-1)/fatt(nu/2)**2
                     nu = nuhi - kik
                     bum1 = bum1*fatt(nu+1)/fatt((nu+2)/2)**2
                     nu = nuhi + kik
                     bum1 = bum1/(fatt(nu-1)/fatt(nu/2)**2)
                     almn(inum) = bum1/(nu-1)
 20               continue
c
c    almn is the a(lambda,mu,nu) coeffficient
c
                  min = mmx + 1
                  mmx = min - 1 + nbcon(im)
                  if (im.eq.i) then
                     mmxp = k
                  else
                     mmxp = mmx
                  end if
                  do 90 m = min , mmxp
                     nrowm = nrow(m)
                     if (m.eq.k) then
                        nmx = l
                     else
                        nmx = m
                     end if
                     do 80 n = min , nmx
                        nrown = nrow(n)
                        j = j + 1
c
c    j is the label of the matrix elements to be calculated
c    i=lambda+1,k=p,l=q, im=mu+1,m=r,n=s
c
c   this part computes the integrals j and k and the supermatrices
c
                        if (numvar.ne.0) then
                           if (k.ne.nbvar1) then
                              if (l.ne.nbvar1) then
                                 if (m.ne.nbvar1) then
                                    if (n.eq.nbvar1) go to 80
                                 end if
                              end if
                           end if
                        end if
                        pr = dzero
                        qr = dzero
                        do 70 nnp = 1 , nrowk
                           npp = itran(nnp,k)
                           np = pqn(npp)
                           np2p1 = 2*np + 1
                           wp = expa(npp)
                           ctrp = ctran(nnp,k)
                           do 60 nnq = 1 , nrowl
                              nqq = itran(nnq,l)
                              nq = pqn(nqq)
                              nq2p1 = 2*nq + 1
                              npq = np + nq
                              wq = expa(nqq)
                              wpq = half*(wp+wq)
                              ctrq = ctran(nnq,l)
                              do 50 nnr = 1 , nrowm
                                 ctrr = ctran(nnr,m)
                                 nrr = itran(nnr,m)
                                 nr = pqn(nrr)
                                 nr2p1 = 2*nr + 1
                                 wr = expa(nrr)
                                 npr = np + nr
                                 nqr = nq + nr
                                 wpr = half*(wr+wp)
                                 wqr = half*(wq+wr)
                                 do 40 nns = 1 , nrown
                                    nss = itran(nns,n)
                                    ns = pqn(nss)
                                    ns2p1 = 2*ns + 1
                                    ws = expa(nss)
                                    nqs = nq + ns
                                    nps = np + ns
                                    wps = half*(wp+ws)
                                    wqs = half*(wq+ws)
                                    cccc = ctrp*ctrq*ctrr*ctran(nns,n)
                                    nrs = nr + ns
                                    wrs = half*(wr+ws)
                                    if (nflag1.eq.1) then
                                       vp = facto(np2p1-1)/wp**(np+half)
                                       vp = vp*facto(nq2p1-1)
     +                                    /wq**(nq+half)
                                       vp = vp*facto(nr2p1-1)
     +                                    /wr**(nr+half)
                                       vp = vp*facto(ns2p1-1)
     +                                    /ws**(ns+half)
                                       sss = cons
                                       vj = giny(m0,npq,nrs,wpq,wrs)
     +                                    + giny(m0,nrs,npq,wrs,wpq)
                                    else
                                       vp = fatt(np2p1)/wp**np2p1
                                       vp = vp*fatt(nq2p1)/wq**nq2p1
                                       vp = vp*fatt(nr2p1)/wr**nr2p1
                                       vp = vp*fatt(ns2p1)/ws**ns2p1
                                       sss = done
                                       vj = two*
     +                                    (siny(m0,npq,nrs,wpq,wrs)
     +                                    +siny(m0,nrs,npq,wrs,wpq))
                                    end if
                                    term2 = sss/dsqrt(vp)
                                    vj = vj*term2
c
c    xj is the integral j0(lambda,p,q,mu,r,s)
c
                                    sss = dzero
                                    ssu = dzero
                                    inum = 0
                                    do 30 kik = nulo , nuhi , 2
                                       inum = inum + 1
                                       nu = kik - 1
                                       if (nflag1.eq.1) then
                                         wjnu = half*
     +                                      (giny(nu,npr,nqs,wpr,wqs)
     +                                      +giny(nu,nqs,npr,wqs,wpr)
     +                                      +giny(nu,nps,nqr,wps,wqr)
     +                                      +giny(nu,nqr,nps,wqr,wps))
     +                                      *term2
                                       else
                                         wjnu = siny(nu,npr,nqs,wpr,wqs)
     +                                      + siny(nu,nqs,npr,wqs,wpr)
     +                                      + siny(nu,nps,nqr,wps,wqr)
     +                                      + siny(nu,nqr,nps,wqr,wps)
                                         wjnu = wjnu*term2
                                       end if
c
c    yjnu is the integral k(nu,lambda,p,q,mu,r,s)
c
                                       ssu = ssu + wjnu*almn(inum)
                                       if (nosht.gt.0) then
                                         if (i.lt.im) go to 30
                                         if (i.eq.im) then
                                         if (nu.lt.2) then
                                         if (i.lt.2) then
                                         akoo = ajmn(1)
                                         else if (i.eq.2) then
                                         akoo = ajmn(3)
                                         else
                                         akoo = ajmn(8)
                                         end if
                                         else if (nu.eq.2) then
                                         if (i.lt.2) go to 30
                                         if (i.eq.2) then
                                         akoo = ajmn(4)
                                         else
                                         akoo = ajmn(9)
                                         end if
                                         else
                                         akoo = ajmn(10)
                                         end if
                                         else if (nu.lt.2) then
                                         if (i.lt.2) go to 30
                                         if (i.eq.2) then
                                         akoo = ajmn(2)
                                         else
                                         akoo = ajmn(6)
                                         end if
                                         else if (nu.eq.2) then
                                         akoo = ajmn(5)
                                         else
                                         akoo = ajmn(7)
                                         end if
                                         sss = sss + akoo*wjnu
                                       end if
 30                                 continue
                                    pr = pr + (vj-half*ssu)*cccc
                                    qr = qr - half*sss*cccc
                                    if (im.eq.2) then
                                       if (i.eq.3) then
                                         if (ajmn(11).ne.dzero) then
                                         if (nflag1.eq.1) then
                                         vjnu = giny(m2,npq,nrs,wpq,wrs)
     +                                      + giny(m2,nrs,npq,wrs,wpq)
                                         else
                                         vjnu = two*
     +                                      (siny(m2,npq,nrs,wpq,wrs)
     +                                      +siny(m2,nrs,npq,wrs,wpq))
                                         end if
                                         qr = qr + ajmn(11)*vjnu*term2
                                         end if
                                       end if
                                    end if
c
c    ajmn(11) is the parameter j(d,p,2) occurring if both p and d are
c    open shells
c
 40                              continue
 50                           continue
 60                        continue
 70                     continue
                        p(j) = pr
                        q(j) = qr
 80                  continue
 90               continue
 100           continue
 110        continue
 120     continue
 130  continue
c      write   p   q    to tape  around  here
      do 140 i = 1 , nsym
         nbas(i) = nbcon(i)
 140  continue
c
      return
      end
**==ddenlp.f
      subroutine ddenlp(czeta, ctran, vect, pop, scr,
     & grid, data)
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
c calculate gradient of the density on a grid
c
      dimension czeta(*), ctran(*), vect(*), grid(*), data(*), pop(*)
      dimension scr(*)
c
      common/craypk/nbuf(9,maxorb),nbasin,norb,lenbas,ngrps,non,newbas,
     &     ncol,jcol,ngrp(3),iblkg,iblka
      common/gjs/noc(maxorb)
      data thresh/1.0d-8/
c
c - initialisation
c
      do 5 i = 1 , ncol
         noc(i)=0
         if (dabs(pop(i)).gt.thresh) noc(i) = 1
 5    continue
c
c scratch memory
c
      i0 = 1
      i1 = i0 + ncol
      i2 = i1 + ncol
      i3 = i2 + ncol
      i4 = i3 + ncol
c
      imo = 0
      call getnpt(np,npx,npy,npz)
      do 10 i = 1,np
         call getpt(i,ax,ay,az,grid,3,idat)
         call dplotd(vect, ctran, czeta, pop,
     &        scr(i0),scr(i1),scr(i2),scr(i3),scr(i4),
     &        ax, ay, az,
     &        data(idat), data(idat+1), data(idat+2),dum,imo)
 10   continue
      return
      end
**==denlp.f
      subroutine denlp(czeta, ctran, vect, pop, scr, grid, data,
     +                 core,lmem)
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
      dimension czeta(*), ctran(*), vect(*), pop(*), scr(*),
     &     grid(*), data(*)
      dimension core(*)
      common/craypk/nbuf(9,maxorb),nbasin,norb,lenbas,ngrps,non,newbas,
     &     ncol,jcol,ngrp(3),iblkg,iblka
      common/gjs/noc(maxorb)
      common/tmdata/isp2(5),otran,isp4(2)
      data thresh/1.0d-8/
c
c  transition dens. plotting has all its information in itmfile
c  and does not need anything from here with exception of
c  x,y,z and inputted task specification
c
      if(otran) call tmpl00(core,lmem)
c
c
c - initialisation
c
      do 5 i = 1 , ncol
         noc(i)=0
         if (dabs(pop(i)).gt.thresh) noc(i) = 1
 5    continue

c
      i0 = 1
      i1 = i0 + ncol
c
c flag density calc
c
      imo = 0
c
      call getnpt(np,npx,npy,npy)

      do 10 i = 1,np
         call getpt(i,ax,ay,az,grid,1,idat)
         if(otran) then
         data(idat)=tmpl1(core,ax,ay,az)
         else
         call dplota(vect, ctran, czeta, pop, scr(i0), scr(i1),
     &        ax, ay, az, data(idat), imo)
         endif
 10   continue
      return
      end
**==densi.f
      subroutine densi(dt,dos)
c
c    one configuration
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
      dimension dt(*),dos(*)
      common/junk/pqn(99),expa(99),ctran(mxprms,30),ajmn(12),czn,
     *expmn,c(15,15),itran(mxprms,30),nrow(30),nbcon(4),nbas(4),
     *nsym,nosh(4),ncsh(4),nccup(4),ibase(33),
     *popul(4),facto(30),fatt(30),thre(15),a(15,15),eivr(15,15),
     *cm2scf(15,15),cm1scf(15,15),eps(15),gmin(15,15),smin(15,15),
     *populs(4),term,energ,pot,cin,vir,thrscf,enorm,
     *red,expdif,encor,
     *nsize(10),nsavf(4),n1(6),nbvar(30),nbas1,nbast,nbvar1,
     *nblock,numvar,nrun,nitscf,nitsc1,nitsc2,nconv,nalarm,no,nstep,
     *iiii,nsht,nosht,nflag3,nbasm,mvar,iex,newbas,nit,nbpos(8)
      data dzero,two,six,ten,fourt/
     *0.0d0,2.0d0,6.0d0,10.0d0,14.0d0/
c
c    compute the density matrices
c
      populs(1) = two
      populs(2) = six
      populs(3) = ten
      populs(4) = fourt
c
      k = 1
      nstep = 0
      j1 = 0
      do 50 i = 1 , nsym
         is = 0
         ik = 0
         nbas1 = nbas(i)
         noshic = nosh(i)
         ncshic = ncsh(i)
         nocsh = noshic + ncshic
         j1 = j1 + nocsh
         if (noshic.gt.0) is = 1
         if (ncshic.gt.0) ik = 1
         do 40 m = 1 , nbas1
            do 30 n = 1 , m
               term1 = popul(i)*c(m,j1)*c(n,j1)
               term2 = dzero
               do 20 j = 1 , ncshic
                  j2 = j + nstep
                  term2 = term2 + c(m,j2)*c(n,j2)
 20            continue
               term2 = term2*populs(i)
               if (m.ne.n) then
                  term1 = term1 + term1
                  term2 = term2 + term2
               end if
               dosa = term1*is
               dt(k) = term2*ik + dosa
               dos(k) = dosa
               k = k + 1
 30         continue
 40      continue
         nstep = j1
 50   continue
c
      return
      end
**==dipnuc.f
      subroutine dipnuc(noc,vlist,c,dint,nocmx)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension dint(3),vlist(nocmx,4),c(3)
c...
c...compute nuclear contributions to dipole moments.
c...x,y,z.
c...
      do 20 i = 1 , 3
         dint(i) = 0.0d0
 20   continue
c...
      do 40 i = 1 , noc
         do 30 j = 1 , 3
            dint(j) = dint(j) + vlist(i,4)*(vlist(i,j)-c(j))
 30      continue
 40   continue
c...
      return
      end
      function discal(x1, y1, z1, x2, y2, z2)
c
c find distance from point (x1, y1, z1) to point (x2, y2, z2)
c
      implicit real*8  (a-h, o-z)
      xdist = (x2 - x1) ** 2
      ydist = (y2 - y1) ** 2
      zdist = (z2 - z1) ** 2
      discal = dsqrt(xdist + ydist + zdist)
      return
      end
**==divpt.f
      subroutine divpt (a,alpha,b,beta,ecoord,eexp,efact)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension a(3), b(3), ecoord(3)
c
      do 20 i = 1 , 3
         ecoord(i) = (alpha*a(i)+beta*b(i))/(alpha+beta)
 20   continue
      eexp = alpha + beta
      absqd = 0.0d0
      do 30 i = 1 , 3
         absqd = absqd + (b(i)-a(i))*(b(i)-a(i))
 30   continue
      efact = dexp(-alpha*beta*absqd/(alpha+beta))
c
      return
      end
**==dmolp.f
      subroutine dmolp(czeta, ctran, vect, scr, grid, data, imo)
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
c analytical derivitive of mo amplitude on a grid
c
      common/craypk/nbuf(9,maxorb),nbasin,norb,lenbas,ngrps,non,newbas,
     &     ncol,jcol,ngrp(3),iblkg,iblka
      dimension czeta(*), ctran(*), vect(*), scr(*), grid(*), data(*)
c
c scratch memory
c
      i0 = 1
      i1 = i0 + ncol
      i2 = i1 + ncol
      i3 = i2 + ncol
      i4 = i3 + ncol
      call getnpt(np,npx,npy,npz)
      do 10 i = 1,np
         call getpt(i,ax,ay,az,grid,3,idat)
         call dplotd(vect, ctran, czeta, pop,
     &        scr(i0),scr(i1),scr(i2),scr(i3),scr(i4),
     &        ax, ay, az,
     &        data(idat), data(idat+1), data(idat+2),dum,imo)
 10   continue
      return
      end
**==dmtx2.f
      subroutine dmtx2(d,v,p,ia,m,n,ndim)
      implicit real*8  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension ia(*),d(*),v(ndim,*),p(*)
      data dzero /0.0d0/
      call vclr(d,1,ia(n+1))
      if(m.eq.0) return
      ii=0
      do 120 i = 1,n
      do 130 k = 1,m
      dum = p(k)*v(i,k)
      if(dum.ne.dzero) then
      call daxpy(i,dum,v(1,k),1,d(ii+1),1)
         endif
  130 continue
  120 ii=ii+i
      return
      end
**==dplota.f
      subroutine dplota(vect, ctran, czeta, pop, v, gt,
     &     grdx, grdy, grdz, den, imo)
c
c Calculate the density at x, y and z
c The point at which this is to be calculated is given by
c grdx, grdy and grdz
c The results are returned in den
c Vect are the MO's
c ctran are the contraction coefficients for the a.o.'s
c czeta are the exponents
c
c if imo is non-zero we will just calculate a MO amplitude
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
      parameter (mcen32=3*maxat-2)
      dimension vect(*),czeta(*),ctran(*),pop(*),v(*),gt(*)
      common/bufb/rbuf(72),ibuf(56),
     +            pa(18),cxi(3),cxj(3),b27,side,
     +            d(1),pp(1),ppp(mcen32),
     +            e(1),qq(1),qqq(mcen32),
     +            f(1),rr(1),rrr(mcen32),p1(200),p2(200),p3(200)
      common/miscop/p(maxat),q(maxat),r(maxat),chg(maxat),
     +             cx(maxorb),cy(maxorb),cz(maxorb),popocc(maxorb)
      common/craypk/nterm(maxorb),icentr(maxorb),itype(maxorb),
     +             mx(maxorb),my(maxorb),mz(maxorb),
     +             nx(maxorb),ny(maxorb),nz(maxorb),
     +     nbasin,norb,lenbas,ngrps,non,newbas,
     &     ncol,jcol,ngrp(3),iblkg,iblka
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      common/gjs/noc(maxorb),noca(maxorb)
      data dzero/0.0d0/
c
c Loop over atoms calculates x, x**2 . It only works up to d's
c
      nsp=1

      do 61 ii=1,non
         pp(nsp)=grdx-p(ii)
         ppp(nsp)=pp(nsp)*pp(nsp)
         qq(nsp)=grdy-q(ii)
         qqq(nsp)=qq(nsp)*qq(nsp)
         rr(nsp)=grdz-r(ii)
         rrr(nsp)=rr(nsp)*rr(nsp)
         v(ii)=rrr(nsp)+ppp(nsp)+qqq(nsp)
 61   nsp=nsp+3
c
c Loop over atomic orbitals
c
      do 62 jj = 1,norb
         pow1 = v(icentr(jj))
c
c Loop over primitive orbitals
c
         ctj = dzero
         nsp=(jj-1)*mxprms
         do 63 ii=1,nterm(jj)
            nsp = nsp+1
            ctj = ctj + ctran(nsp) * dexp(-czeta(nsp)*pow1)
 63         continue

c
c This bit is almost too clever
c It multiplies the radial part of the distribution by the angular
c component ( in a cartesian sense ).
c mx(jj) access the correct power of x in the array d.
c
         gt(jj) = ctj * d(mx(jj)) * e(my(jj)) * f(mz(jj))
 62   continue
c
c Multiply by mo coefficients
c
      den=dzero
      if(imo.eq.0)then
         do 64 ii=1,ncol
            if(noc(ii) .ge. 0) then
               nsp=ilifq(ii)
               b=ddot(norb,gt,1,vect(nsp+1),1)
               den = b * b * pop(ii) + den
            endif
 64      continue
      else
         den = ddot(norb,gt,1,vect(ilifq(imo)+1),1)
      endif
c
c Result is in den
c
      return
      end
      subroutine dplotd(vect, ctran, czeta, pop, v, gt, gtx, gty, gtz,
     &     grdx, grdy, grdz, denx, deny, denz, den, imo)
c
c Calculate the derivative of the density wrt x, y and z
c The point at which this is to be calculated is given by
c grdx, grdy and grdz
c The results are returned in denx,deny,denz
c Vect are the MO's
c ctran are the contraction coefficients for the a.o.'s
c czeta are the exponents
c
c if imo != 0, we will just calc a MO derivative
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
      parameter (mcen32=3*maxat-2)
      dimension vect(*),czeta(*),ctran(*), pop(*), v(*),
     +     gt(*),gtx(*),gty(*),gtz(*)
      common/bufb/rbuf(72),ibuf(56),
     +            pa(18),cxi(3),cxj(3),b27,side,
     +            d(1),pp(1),ppp(mcen32),
     +            e(1),qq(1),qqq(mcen32),
     +            f(1),rr(1),rrr(mcen32),p1(200),p2(200),p3(200)
      common/miscop/p(maxat),q(maxat),r(maxat),chg(maxat),
     +             cx(maxorb),cy(maxorb),cz(maxorb),popocc(maxorb)
      common/craypk/nterm(maxorb),icentr(maxorb),itype(maxorb),
     +             mx(maxorb),my(maxorb),mz(maxorb),
     +             nx(maxorb),ny(maxorb),nz(maxorb),
     +     nbasin,norb,lenbas,ngrps,non,newbas,
     &     ncol,jcol,ngrp(3),iblkg,iblka
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      common/gjs/noc(maxorb)
      data m1,mm1,dzero/1,-1,0.0d0/
c
c Loop over atoms calculates x, x**2 . It only works up to d's
c
      nsp=1
      do 61 ii=1,non
         pp(nsp)=grdx-p(ii)
         ppp(nsp)=pp(nsp)*pp(nsp)
         qq(nsp)=grdy-q(ii)
         qqq(nsp)=qq(nsp)*qq(nsp)
         rr(nsp)=grdz-r(ii)
         rrr(nsp)=rr(nsp)*rr(nsp)
         v(ii)=rrr(nsp)+ppp(nsp)+qqq(nsp)
 61   nsp=nsp+3
c
c Loop over atomic orbitals
c
      do 62 jj = 1,norb
         ctj = dzero
         ctjx= dzero
         ctjy= dzero
         ctjz= dzero
         nsp = (jj+mm1)*mxprms
         pow1 = v(icentr(jj))
         powx = 2.0d0 * pp((icentr(jj)-1)*3+1)
         powy = 2.0d0 * qq((icentr(jj)-1)*3+1)
         powz = 2.0d0 * rr((icentr(jj)-1)*3+1)
c
c Loop over primitive orbitals
c
         do 63 ii=1,nterm(jj)
            nsp = nsp+m1
c Radial part
            ctj = ctj + ctran(nsp) * dexp(-czeta(nsp)*pow1)
c and its derivative
            ctjx= ctjx+ ctran(nsp) * dexp(-czeta(nsp)*pow1) *
     -                  ( -czeta(nsp) * powx )
            ctjy= ctjy+ ctran(nsp) * dexp(-czeta(nsp)*pow1) *
     -                  ( -czeta(nsp) * powy )
            ctjz= ctjz+ ctran(nsp) * dexp(-czeta(nsp)*pow1) *
     -                  ( -czeta(nsp) * powz )
63       continue
c
c This bit is almost too clever
c It multiplies the radial part of the distribution by the angular
c component ( in a cartesian sense ).
c mx(jj) access the correct power of x in the array d.
c
         gt(jj) = ctj * d(mx(jj)) * e(my(jj)) * f(mz(jj))
         gtx(jj) = ctjx * d(mx(jj)) * e(my(jj)) * f(mz(jj))
         gty(jj) = ctjy * d(mx(jj)) * e(my(jj)) * f(mz(jj))
         gtz(jj) = ctjz * d(mx(jj)) * e(my(jj)) * f(mz(jj))
         mmx = mod(mx(jj)-1,3)
         mmy = mod(my(jj)-1,3)
         mmz = mod(mz(jj)-1,3)
         if(mmx.gt.0)
     +       gtx(jj)=gtx(jj)+ctj*d(mx(jj)-1)*e(my(jj))*f(mz(jj))*mmx
         if(mmy.gt.0)
     +       gty(jj)=gty(jj)+ctj*d(mx(jj))*e(my(jj)-1)*f(mz(jj))*mmy
         if(mmz.gt.0)
     +       gtz(jj)=gtz(jj)+ctj*d(mx(jj))*e(my(jj))*f(mz(jj)-1)*mmz
 62   continue
c
c Multiply by mo coefficients
c
      den=dzero
      denx=dzero
      deny=dzero
      denz=dzero
c
c density loop
c
      if(imo.eq.0)then
         do 64 ii=1,norb
            if(noc(ii) .ge. 0) then
               nsp=ilifq(ii)
               b=ddot(norb,gt,1,vect(nsp+1),1)
               bx=ddot(norb,gtx,1,vect(nsp+1),1)
               by=ddot(norb,gty,1,vect(nsp+1),1)
               bz=ddot(norb,gtz,1,vect(nsp+1),1)
               den = b * b * pop(ii) + den
               denx = 2.0d0 * b * bx * pop(ii) + denx
               deny = 2.0d0 * b * by * pop(ii) + deny
               denz = 2.0d0 * b * bz * pop(ii) + denz
            endif
 64      continue
      else
c
c orbital calc
c
         nsp=ilifq(imo)
         den=ddot(norb,gt,1,vect(nsp+1),1)
         denx=ddot(norb,gtx,1,vect(nsp+1),1)
         deny=ddot(norb,gty,1,vect(nsp+1),1)
         denz=ddot(norb,gtz,1,vect(nsp+1),1)
      endif
c
c Result is in den, denx, deny, denz
c
      return
      end
**==dpogr2.f
      subroutine dpogr2(istar,ien,m225,grid,dens,ga,data)
      implicit real*8  (a-h,o-z)
      logical iandj,norm,double
      character*10 charwall
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
      real*8 pp, u, w, dji
      integer nroots, nroot2
      common /root/ pp,u(12),w(12),nroots,nroot2,dji(60)
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
      common/junk/s(225),g(225),
     + xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj,
     + tol,ii,jj,lit,ljt,mini,minj,maxi,maxj,iandk
      common/junk2/dxyz(4),ggg(225),ft(225),dij(225),
     + xin(125),yin(125),zin(125),
     + ijx(225),ijy(225),ijz(225)
      dimension dens(*), ga(*), data(*), grid(*)
      dimension ix(35),iy(35),iz(35),jx(35),jy(35),jz(35)
      data zero,one /0.0d+00,1.0d+00 /
      data pi212 /1.1283791670955d0/
      data sqrt3 /1.73205080756888d0/
      data sqrt5 /2.23606797749979d0/
      data sqrt7 /2.64575131106459d0/
      data rln10 /2.30258d+00/
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
c
      tol=rln10*itol
      norm=normf.ne.1.or.normp.ne.1
      ngrid = ien - istar + 1
      np36=ngrid * m225
c
c     ----- ishell
c
      do 9000 ii=1,nshell
      i=katom(ii)
      xi=c(1,i)
      yi=c(2,i)
      zi=c(3,i)
      i1=kstart(ii)
      i2=i1+kng(ii)-1
      lit=ktype(ii)
      mini=kmin(ii)
      maxi=kmax(ii)
      loci=kloc(ii)-mini
c
c     ----- jshell
c
      do 8000 jj=1,ii
      j=katom(jj)
      xj=c(1,j)
      yj=c(2,j)
      zj=c(3,j)
      j1=kstart(jj)
      j2=j1+kng(jj)-1
      ljt=ktype(jj)
      minj=kmin(jj)
      maxj=kmax(jj)
      locj=kloc(jj)-minj
      nroots=(lit+ljt-2)/2+1
      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
      iandj=ii.eq.jj
c
c     ----- prepare indices for pairs of (i,j) functions
c
      ij=0
      max=maxj
      do i=mini,maxi
      nnx=ix(i)
      nny=iy(i)
      nnz=iz(i)
      if(iandj) max=i
       do j=minj,max
       ij=ij+1
       ijx(ij)=nnx+jx(j)
       ijy(ij)=nny+jy(j)
       ijz(ij)=nnz+jz(j)
       enddo
      enddo
c
      call vclr(ga,1,np36)
c
c     ----- i primitive
c
      jgmax=j2
      do 7000 ig=i1,i2
      ai=ex(ig)
      arri=ai*rr
      axi=ai*xi
      ayi=ai*yi
      azi=ai*zi
      csi=cs(ig)
      cpi=cp(ig)
      cdi=cd(ig)
      cfi=cf(ig)
      cgi=cg(ig)
c
c     ----- j primitive
c
      if(iandj) jgmax=ig
      do 6000 jg=j1,jgmax
      aj=ex(jg)
      aa=ai+aj
      aa1=one/aa
      dum=aj*arri*aa1
      if(dum.gt.tol) go to 6000
      fac=dexp(-dum)
      csj=cs(jg)
      cpj=cp(jg)
      cdj=cd(jg)
      cfj=cf(jg)
      cgj=cg(jg)
      ax=(axi+aj*xj)*aa1
      ay=(ayi+aj*yj)*aa1
      az=(azi+aj*zj)*aa1
c
c     ----- density factor
c
      double=iandj.and.ig.ne.jg
      max=maxj
      nn=0
      do 170 i = mini,maxi
      go to (200,220,280,280,
     +       240,280,280,260,280,280,
     +       210,280,280,215,280,280,280,280,280,217,
     +       222,280,280,224,280,280,280,280,280,
     +       226,280,280,228,280,280), i
  200 dum1 = csi*fac
      go to 280
  220 dum1 = cpi*fac
      go to 280
  240 dum1 = cdi*fac
      go to 280
  260 if (norm) dum1 = dum1*sqrt3
      go to 280
  210 dum1=cfi*fac
      go to 280
  215 if (norm) dum1=dum1*sqrt5
      go to 280
  217 if (norm) dum1=dum1*sqrt3
      go to 280
  222 dum1=cgi*fac
      go to 280
  224 if (norm) dum1 = dum1*sqrt7
      go to 280
  226 if (norm) dum1 = dum1*sqrt5/sqrt3
      go to 280
  228 if (norm) dum1 = dum1*sqrt3
  280 if (iandj) max = i
      do 170 j = minj,max
      go to (300,340,400,400,
     +       360,400,400,380,400,400,
     +       310,400,400,315,400,400,400,400,400,317,
     +       322,400,400,324,400,400,400,400,
     +       400,326,400,400,328,400,400),j
  300 dum2 = dum1*csj
      if ( .not. double) go to 400
      if (i .gt. 1) go to 320
      dum2 = dum2+dum2
      go to 400
  320 dum2 = dum2+csi*cpj*fac
      go to 400
  340 dum2 = dum1*cpj
      if (double) dum2 = dum2+dum2
      go to 400
  360 dum2 = dum1*cdj
      if (double) dum2 = dum2+dum2
      go to 400
  380 if (norm) dum2 = dum2*sqrt3
      go to 400
  310 dum2=dum1*cfj
      if (double) dum2 = dum2+dum2
      go to 400
  315 if (norm) dum2 = dum2*sqrt5
      go to 400
  317 if (norm) dum2 = dum2*sqrt3
      go to 400
  322 dum2=dum1*cgj
      if (double) dum2 = dum2 + dum2
      go to 400
  324 if (norm) dum2 = dum2*sqrt7
      go to 400
  326 if (norm) dum2 = dum2*sqrt5/sqrt3
      go to 400
  328 if (norm) dum2 = dum2*sqrt3
  400 nn = nn+1
  170 dij(nn) = dum2
c
c     ----- nuclear attraction
c
      dum=pi212*aa1
      do i=1,ij
       dij(i)=dij(i)*dum
      enddo
      aax=aa*ax
      aay=aa*ay
      aaz=aa*az
c
c loop over points
c
      do ipt = 1,ngrid
      znuc=one
      call getpt(istar + ipt - 1,cx,cy,cz,grid,1,idat)
      pp=aa*((ax-cx)**2+(ay-cy)**2+(az-cz)**2)
      if(nroots.le.3)call rt123
      if(nroots.eq.4)call roots4
      if(nroots.eq.5)call roots5
      mm=0
       do k=1,nroots
       uu=aa*u(k)
       ww=w(k)*znuc
       tt=one/(aa+uu)
       t= dsqrt(tt)
       x0=(aax+uu*cx)*tt
       y0=(aay+uu*cy)*tt
       z0=(aaz+uu*cz)*tt
c
       in=-5+mm
        do i=1,lit
        in=in+5
        ni=i
         do j=1,ljt
         jn=in+j
         nj=j
         call stvint
         xin(jn)=xint
         yin(jn)=yint
         zin(jn)=zint*ww
         enddo
        enddo
       mm=mm+25
       enddo
c
       do i=1,ij
       nnx=ijx(i)
       nny=ijy(i)
       nnz=ijz(i)
       dum=zero
       mm=0
        do k=1,nroots
        dum=dum+xin(nnx+mm)*yin(nny+mm)*zin(nnz+mm)
        mm=mm+25
        enddo
       g(i)= dum*dij(i)
       enddo
c
       index=m225*(ipt-1)
       max=maxj
       nn2=0
       do i = mini,maxi
       if (iandj) max = i
        do j = minj, max
        index=index+1
        nn2 = nn2 +1
        ga(index) = ga(index) + g(nn2)
        enddo
       enddo
      enddo
c
 6000 continue
 7000 continue
c
c     ----- set up data
c
      do ipt = 1,ngrid
         index=m225*(ipt-1)
         vt=zero
         max = maxj
         do i=mini,maxi
            li=loci+i
            in=iky(li)
            if (iandj) max=i
            do j =minj,max
               lj= locj + j
               jn = lj + in
               index=index+1
               dw = dens(jn) * ga(index)
               if(li.ne.lj) dw= dw + dw
               vt = vt +dw
            enddo
         enddo
         data(istar+ipt-1)=data(istar+ipt-1)-vt
      enddo

 8000 continue
 9000 continue
      ss=cpulft(1)
      write(iwr,3001)ss,charwall()
 3001 format(' potential batch complete at ',f14.2,' secs ',a10,' wall')
      return
      end
**==dponc2.f
      subroutine dponc2(ngrid, grid, data, inozr)
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
      dimension grid(*), data(*)
c
      data dzero/0.0d0/
c
c ... setup nuclear contribution to potential grid
c
cTWK  The nuclear contribution can be neglected with keyword NOZR.
c     Useful e.g. when combining alpha and beta potentials
c     to avoid double counting the nuclear contribution.
c
      do i=1,ngrid
         call getpt(i,pppp,qqqq,rrrr,grid,1,idat)
         top=dzero
         if (inozr .ne. 1) then
            do ii=1,nat
               ctj=pppp-c(1,ii)
               pow1=qqqq-c(2,ii)
               pow2=rrrr-c(3,ii)
               ctj=ctj*ctj+pow1*pow1+pow2*pow2
               top=top+czan(ii)/dsqrt(ctj)
            enddo
         endif
      data(idat)=top
      enddo
c
      return
      end
      subroutine eldeng(istar,ien,grid,dens,elden,
     +                  xx0,yy0,zz0,data)
      implicit real*8 (a-h,o-z), integer   (i-n)
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
      logical iandj
      logical norm,double
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
      common/xyzden/xint,yint,zint,x0,y0,z0,xi,yi,zi,xj,yj,zj,
     +              ni,nj
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      dimension data(*),grid(*),dens(*), elden(*)
      dimension xx0(*),yy0(*),zz0(*)
      common /junk3/ dij(225),
     +               xin(5,5),yin(5,5),zin(5,5)
      dimension ijx(35),ijy(35),ijz(35)
      character*10 charwall
c
      data zero  /0.0d+00/
      data one   /1.0d+00/
      data rln10 /2.30258d+00/
      data sqrt3 /1.73205080756888d+00/
      data sqrt5 /2.23606797749979d+00/
      data sqrt7 /2.64575131106459d+00/
      data ijx    / 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,
     1              4, 1, 1, 3, 3, 2, 1, 2, 1, 2,
     2              5, 1, 1, 4, 4, 2, 1, 2, 1, 3,
     3              3, 1, 3, 2, 2/
      data ijy    / 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,
     1              1, 4, 1, 2, 1, 3, 3, 1, 2, 2,
     2              1, 5, 1, 2, 1, 4, 4, 1, 2, 3,
     3              1, 3, 2, 3, 2/
      data ijz    / 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,
     1              1, 1, 4, 1, 2, 1, 2, 3, 3, 2,
     2              1, 1, 5, 1, 2, 1, 2, 4, 4, 1,
     3              3, 3, 2, 2, 3/
c
c
      tol=rln10*itol
      norm=normf.ne.1.or.normp.ne.1
      ngrid = ien - istar + 1
      call vclr(elden,1,ngrid)
      do ipt=1,ngrid
      call getpt(istar+ipt-1,xx0(ipt),yy0(ipt),zz0(ipt),grid,1,idat)
      enddo
c
c     ----- ishell -----
c
      do 9000 ii=1,nshell
      i=katom(ii)
      xi=c(1,i)
      yi=c(2,i)
      zi=c(3,i)
      i1=kstart(ii)
      i2=i1+kng(ii)-1
      lit=ktype(ii)
      mini=kmin(ii)
      maxi=kmax(ii)
      loci=kloc(ii)-mini
c
c     ----- jshell -----
c
      do 8000 jj=1,ii
      j=katom(jj)
      xj=c(1,j)
      yj=c(2,j)
      zj=c(3,j)
      j1=kstart(jj)
      j2=j1+kng(jj)-1
      ljt=ktype(jj)
      minj=kmin(jj)
      maxj=kmax(jj)
      locj=kloc(jj)-minj
c
      iandj=ii.eq.jj
      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
c
c     ----- i primitive -----
c
      do 7000 ig=i1,i2
      ai=ex(ig)
      arri=ai*rr
      axi=ai*xi
      ayi=ai*yi
      azi=ai*zi
      csi=cs(ig)
      cpi=cp(ig)
      cdi=cd(ig)
      cfi=cf(ig)
      cgi=cg(ig)
c
c     ----- j primitive -----
c
      jgmax=j2
      if(iandj) jgmax=ig
      do 6000 jg=j1,jgmax
      aj=ex(jg)
      aa=ai+aj
      aa1=one/aa
      dum=aj*arri*aa1
      if(dum.gt.tol) go to 6000
      fac= exp(-dum)
      csj=cs(jg)
      cpj=cp(jg)
      cdj=cd(jg)
      cfj=cf(jg)
      cgj=cg(jg)
      ax=(axi+aj*xj)*aa1
      ay=(ayi+aj*yj)*aa1
      az=(azi+aj*zj)*aa1
c
c     ----- density factor -----
c
      double=iandj.and.ig.ne.jg
      ij=0
      do 360 i=mini,maxi
      go to (110,120,220,220,130,220,220,140,220,220,
     1       150,220,220,160,220,220,220,220,220,170,
     2       180,220,220,190,220,220,220,220,220,200,
     3       220,220,210,220,220),i
  110 dum1=csi*fac
      go to 220
  120 dum1=cpi*fac
      go to 220
  130 dum1=cdi*fac
      go to 220
  140 if(norm) dum1=dum1*sqrt3
      go to 220
  150 dum1=cfi*fac
      go to 220
  160 if(norm) dum1=dum1*sqrt5
      go to 220
  170 if(norm) dum1=dum1*sqrt3
      go to 220
  180 dum1=cgi*fac
      go to 220
  190 if(norm) dum1=dum1*sqrt7
      go to 220
  200 if(norm) dum1=dum1*sqrt5/sqrt3
      go to 220
  210 if(norm) dum1=dum1*sqrt3
  220 continue
c
      jmax=maxj
      if(iandj) jmax=i
      do 360 j=minj,jmax
      go to (230,250,350,350,260,350,350,270,350,350,
     1       280,350,350,290,350,350,350,350,350,300,
     2       310,350,350,320,350,350,350,350,350,330,
     3       350,350,340,350,350),j
  230 dum2=dum1*csj
      if(.not.double) go to 350
      if(i.gt.1) go to 240
      dum2=dum2+dum2
      go to 350
  240 dum2=dum2+csi*cpj*fac
      go to 350
  250 dum2=dum1*cpj
      if(double) dum2=dum2+dum2
      go to 350
  260 dum2=dum1*cdj
      if(double) dum2=dum2+dum2
      go to 350
  270 if(norm) dum2=dum2*sqrt3
      go to 350
  280 dum2=dum1*cfj
      if(double) dum2=dum2+dum2
      go to 350
  290 if(norm) dum2=dum2*sqrt5
      go to 350
  300 if(norm) dum2=dum2*sqrt3
      go to 350
  310 dum2=dum1*cgj
      if(double) dum2=dum2+dum2
      go to 350
  320 if(norm) dum2=dum2*sqrt7
      go to 350
  330 if(norm) dum2=dum2*sqrt5/sqrt3
      go to 350
  340 if(norm) dum2=dum2*sqrt3
  350 continue
c
      ij=ij+1
  360 dij(ij)=dum2
c
c     ----- electron density -----
c     ------loop over points -----
c
      do ipt=1,ngrid
c
c     call getpt(istar+ipt-1,x0,y0,z0,grid,1,idat)
c
      x0 = xx0(ipt)
      y0 = yy0(ipt)
      z0 = zz0(ipt)
      dum = aa*((x0-ax)**2+(y0-ay)**2+(z0-az)**2)
      if(dum.le.tol) then
       fac = exp(-dum)
c
c     ----- density values -----
c
       do j=1,ljt
        nj=j
        do i=1,lit
         ni=i
         call denint2
         xin(i,j)=xint
         yin(i,j)=yint
         zin(i,j)=zint
        enddo
       enddo
c
       ij=0
       max = maxj
       do i=mini,maxi
        ix=ijx(i)
        iy=ijy(i)
        iz=ijz(i)
        li=loci+i
        in=iky(li)
        if(iandj) max=i
        do j=minj,max
         lj=locj + j
         jn=lj + in
         jx=ijx(j)
         jy=ijy(j)
         jz=ijz(j)
         ij=ij+1
         eld=fac*dij(ij)*xin(ix,jx)*yin(iy,jy)*zin(iz,jz)
         den=dens(jn)
         if(li.ne.lj) den=den+den
         elden(ipt)=elden(ipt)+eld*den
        enddo
       enddo
c
      endif
      enddo
c
 6000 continue
 7000 continue
 8000 continue
 9000 continue
c
c     ----- set up data
c
      do ipt = 1,ngrid
       eld = elden(ipt)
       data(istar+ipt-1)=eld
       enddo
c
      ss=cpulft(1)
      write(iwr,3001)ss,charwall()
 3001 format(' density batch complete at ',f8.2,' secs ',a10,' wall')
      return
      end
**==enerat.f
      subroutine enerat(p,q,dt,dos,u,t)
c
c    one configuration
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
      dimension p(*),q(*),dt(*),dos(*),u(*),t(*)
      common/junk/pqn(99),expa(99),ctran(mxprms,30),ajmn(12),czn,
     *expmn,c(15,15),itran(mxprms,30),nrow(30),nbcon(4),nbas(4),
     *nsym,nosh(4),ncsh(4),nccup(4),ibase(33),
     *popul(4),facto(30),fatt(30),thre(15),a(15,15),eivr(15,15),
     *cm2scf(15,15),cm1scf(15,15),eps(15),gmin(15,15),smin(15,15),
     *populs(4),term,energ,pot,cin,vir,thrscf,enorm,
     *red,expdif,encor,
     *nsize(10),nsavf(4),n1(6),nbvar(30),nbas1,nbast,nbvar1,
     *nblock,numvar,nrun,nitscf,nitsc1,nitsc2,nconv,nalarm,no,nstep,
     *iiii,nsht,nosht,nflag3,nbasm,mvar,iex,newbas,nit,nbpos(8)
      data dzero,half/0.0d0,0.5d0/
c
      pot = dzero
      potn = dzero
      cin = dzero
c
      do i = 1 , nit
         energ = dt(i)
         do j = 1 , i
            k = (i*(i-1)/2) + j
            term = energ*p(k)*dt(j) - dos(i)*q(k)*dos(j)
            if (i.eq.j) then
               term = term*half
            end if
            pot = pot + term
         enddo
         potn = potn + u(i)*energ
         cin = cin + t(i)*energ
      enddo
      pot = pot - czn*potn
      energ = cin + pot
      vir = pot/cin
c
      return
      end
**==fldnuc.f
      subroutine fldnuc(noc,vlist,c,dint,nocmx)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      dimension dint(3),vlist(nocmx,4),c(3)
c
      do i = 1 , 3
         dint(i) = 0.0d0
      enddo
      do i = 1 , noc
         dst = (vlist(i,1)-c(1))**2 + (vlist(i,2)-c(2))
     +         **2 + (vlist(i,3)-c(3))**2
         if (dst.gt.1.0d-10) then
            do j = 1 , 3
               dint(j) = dint(j) + vlist(i,4)*(vlist(i,j)-c(j))
     +                   /(dsqrt(dst))**3
             enddo
         end if
      enddo
c
      do l = 1 , 3
         dint(l) = -dint(l)
      enddo
c
      return
      end
**==fmc.f
      subroutine fmc(mvar,varx,expmx,fmch)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      m = mvar
      var = varx
      expmx = dexp(-var)
      if (var.le.20.0d0) then
c         x less than or equal to 20.0.
         a = m
         a = a + 0.5d0
         term = 1.0d0/a
         ptlsum = term
         do 20 i = 2 , 60
            a = a + 1.0d0
            term = term*var/a
            ptlsum = ptlsum + term
            if (term/ptlsum.lt.0.00000001d0) go to 50
 20      continue
         write (iwr,6010) m , var
      else
c         x greater than 20.0.
         a = m
         vard = 1.0d0/var
         b = a + 0.5d0
         a = a - 0.5d0
         approx = 0.88622692d0*(dsqrt(vard)*vard**m)
         if (m.ne.0) then
            do 30 i = 1 , m
               b = b - 1.0d0
               approx = approx*b
 30         continue
         end if
         fimult = 0.5d0*expmx*vard
         fiprop = fimult/approx
         term = 1.0d0
         ptlsum = term
         notrms = var
         notrms = notrms + m
         do 40 i = 2 , notrms
            term = term*a*vard
            ptlsum = ptlsum + term
            if (dabs(term*fiprop/ptlsum).le.0.00000001d0) go to 60
            a = a - 1.0d0
 40      continue
         write (iwr,6010) m , var
         go to 60
      end if
 50   fmch = 0.5d0*ptlsum*expmx
      return
 60   fmch = approx - fimult*ptlsum
c
      return
 6010 format (' *** warning == no convergence for fmch: ',i6,e16.9)
      end
**==get3ix.f
      subroutine get3ix(ipt,ix,iy,iz)
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
      common/dfgrid/geom(12,mxgrid),gdata(5,mxgrid),igtype(mxgrid),
     &             igsect(mxgrid),nptot(mxgrid),npt(3,mxgrid),
     &             igdata(5,mxgrid),igstat(mxgrid),ngdata(mxgrid),
     &             ngrid
      common/sgrd/ stepx(3),stepy(3),stepz(3),igrd
c
      ip = ipt - 1
      if(igtype(igrd).eq.1)then
         iy = ip/npt(1,igrd)
         ix = ip - iy*npt(1,igrd)
         iz = 0
      else if(igtype(igrd).eq.2)then
         iz = ip/(npt(1,igrd)*npt(2,igrd))
         ir = ip - iz*npt(1,igrd)*npt(2,igrd)
         iy = ir/npt(1,igrd)
         ix = ir - iy*npt(1,igrd)
      else
         call caserr('bad get3ix call')
      endif
      return
      end
**==getgty.f
      subroutine getgty(i)
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
      common/dfgrid/geom(12,mxgrid),gdata(5,mxgrid),igtype(mxgrid),
     &             igsect(mxgrid),nptot(mxgrid),npt(3,mxgrid),
     &             igdata(5,mxgrid),igstat(mxgrid),ngdata(mxgrid),
     &             ngrid
      common/sgrd/ stepx(3),stepy(3),stepz(3),igrd
      i = igtype(igrd)
      return
      end
**==getnpt.f
      subroutine getnpt(n,nx,ny,nz)
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
      common/dfgrid/geom(12,mxgrid),gdata(5,mxgrid),igtype(mxgrid),
     &             igsect(mxgrid),nptot(mxgrid),npt(3,mxgrid),
     &             igdata(5,mxgrid),igstat(mxgrid),ngdata(mxgrid),
     &             ngrid
      common/sgrd/ stepx(3),stepy(3),stepz(3),igrd
      n = nptot(igrd)
      nx = npt(1,igrd)
      ny = npt(2,igrd)
      nz = npt(3,igrd)
      return
      end
**==getpro.f
      subroutine getpro(dint,zcentg,ipos,itpro,mop,iwr)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      dimension dint(*)
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
      common/miscop/pnuc(15),
     *mopp,ncom,len,lentri
      common/junkc/zheadt(3),zprop(5),ztitp(10),ztagp,zspip
      integer numdu, iblkdu, maxb, iblkla
      logical orevis
      real*8 apos
      common/sector/numdu,iblkdu,orevis(2),
     + apos(508),maxb,iblkla
c
      integer mp, lwordp, lvec, ied3, newbas, lenb
      integer nbfnd, lenbad, lentrd, len3, ngrpm
      integer idena, nwroot, nosec
      logical noscf
      common /limy/ mp(10),lwordp,lvec,ied3,newbas,lenb,
     +              nbfnd,lenbad,lentrd,len3,ngrpm,
     +              idena,noscf,nwroot,nosec(mxroot)
c
c
      call secget(ipos,itpro,iblk)
      write (iwr,6010) ipos , iblk
      mword1 = 20
      call rdchr(zheadt,mword1,iblk,numdu)
      mword2 = 15 + 4/lenwrd()
      call reads(pnuc,mword2,numdu)
      write (iwr,6020) zprop , ztagp , zheadt
      if (mop.ne.mopp .or. ztagp.ne.zcentg)
     +    call caserr('invalid property section retrieved')
      iblk = iblk + 2
      do 20 i = 1 , ncom
         dint(mp(i)) = pnuc(i)
         int4 = mp(i) + 3
         call rdedx(dint(int4),len,iblk,numdu)
         iblk = iblk + lentri
 20   continue
c
      return
 6010 format (/' integrals restored from dumpfile :section ',i3,
     +        ' (at block ',i6,')')
 6020 format (/1x,5a6,' (centre ',a8,')'/' created by job ',a8,
     +        ' :  at ',a8,' on ',a8//)
      end
**==getpt.f
      subroutine getpt(ipt,ax,ay,az,q,nel,idat)
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
      common/dfgrid/geom(12,mxgrid),gdata(5,mxgrid),igtype(mxgrid),
     &             igsect(mxgrid),nptot(mxgrid),npt(3,mxgrid),
     &             igdata(5,mxgrid),igstat(mxgrid),ngdata(mxgrid),
     &             ngrid
      common/sgrd/ stepx(3),stepy(3),stepz(3),igrd
c
      dimension q(*)
      ip = ipt - 1
      if(igtype(igrd).eq.1)then
         iy = ip/npt(1,igrd)
         ix = ip - iy*npt(1,igrd)
         ax = geom(1,igrd) + ix*stepx(1) + iy*stepy(1)
         ay = geom(2,igrd) + ix*stepx(2) + iy*stepy(2)
         az = geom(3,igrd) + ix*stepx(3) + iy*stepy(3)
      else if(igtype(igrd).eq.2)then
         iz = ip/(npt(1,igrd)*npt(2,igrd))
         ir = ip - iz*npt(1,igrd)*npt(2,igrd)
         iy = ir/npt(1,igrd)
         ix = ir - iy*npt(1,igrd)
         ax = geom(1,igrd) + ix*stepx(1) + iy*stepy(1) + iz*stepz(1)
         ay = geom(2,igrd) + ix*stepx(2) + iy*stepy(2) + iz*stepz(2)
         az = geom(3,igrd) + ix*stepx(3) + iy*stepy(3) + iz*stepz(3)
      else
         ax = q(ip*3+1)
         ay = q(ip*3+2)
         az = q(ip*3+3)
      endif
      idat=nel*ip+1
      return
      end
**==ggntde.f
       function ggntde (i,iski,nceni,j,iskj,ncenj)
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
      common/junk/eta(mxgaus),wnorm(mxgaus),vlist(maxat,4)
      common/miscop/aaa(80),aa(3),bb(3),e(3),a(10),b(10)
c
      do 20 m = 1 , 3
         aa(m) = vlist(nceni,m)
         bb(m) = vlist(ncenj,m)
 20   continue
      do 30 m = 1 , 10
         a(m) = 0.0d0
         b(m) = 0.0d0
 30   continue
      dumx = eta(i)
      dumy = eta(j)
      a(iski) = 1.0d0
      b(iskj) = 1.0d0
      call divpt(aa,dumx,bb,dumy,e,eexp,efact)
      call shftce(a,aa,e,iski)
      call shftce(b,bb,e,iskj)
      fexp = 2.0d0*eexp
      ggntde = 0.0d0
      ggntde = (a(1)*b(1)
     +         +(a(2)*b(2)+a(3)*b(3)+a(4)*b(4)+a(1)*(b(5)+b(6)+b(7))
     +         +b(1)*(a(5)+a(6)+a(7))
     +         +(3.0d0*(a(5)*b(5)+a(6)*b(6)+a(7)*b(7))+a(8)*b(8)+a(9)
     +         *b(9)+a(10)*b(10)+a(5)*(b(6)+b(7))+a(6)*(b(5)+b(7))+a(7)
     +         *(b(5)+b(6))+ggntde)/fexp)/fexp)
      ggntde = efact*15.7496099d0*ggntde/(fexp*dsqrt(fexp))
c
      return
      end
**==giny.f
      function giny(ny,nab,ncd,vab,vcd)
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
      common/junk/pqn(99),expa(99),ctran(mxprms,30),ajmn(12),czn,
     *expmn,c(15,15),itran(mxprms,30),nrow(30),nbcon(4),nbas(4),
     *nsym,nosh(4),ncsh(4),nccup(4),ibase(33),
     *popul(4),facto(30),fatt(30),thre(15),a(15,15),eivr(15,15),
     *cm2scf(15,15),cm1scf(15,15),eps(15),gmin(15,15),smin(15,15),
     *populs(4),term,energ,pot,cin,vir,thrscf,enorm,
     *red,expdif,encor,
     *nsize(10),nsavf(4),n1(6),nbvar(30),nbas1,nbast,nbvar1,
     *nblock,numvar,nrun,nitscf,nitsc1,nitsc2,nconv,nalarm,no,nstep,
     *iiii,nsht,nosht,nflag3,nbasm,mvar,iex,newbas,nit,nbpos(8)
c
c    calculates the integral i(ny) needed for the supermatrices p and q
c    when a gaussian basis set is used
c    *******
      data dzero,done,half/
     *0.0d0,1.0d0,0.5d0/
c
      tabcd = vab/vcd
      nabny = nab - ny - 1
      ncdny = ncd + ny
      vnab = facto(nabny)/vab**(half*(nabny+1))
      vncd = facto(ncdny)/vcd**(half*(ncdny+1))
      term1 = dzero
      do 20 i = 1 , nabny , 2
         term1 = term1 + facto(i+ncdny-1)/(facto(i)*facto(ncdny))
     +           *(tabcd/(done+tabcd))**(half*(i-1))
 20   continue
      cabcd = term1*(done+tabcd)**(-half*(ncdny+1))
      giny = vnab*vncd*cabcd
c
      return
      end
*==graph.f
      subroutine graph(q,ocopy)
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
c grid definition parameters
c
      common/dfgrid/geom(12,mxgrid),gdata(5,mxgrid),igtype(mxgrid),
     &             igsect(mxgrid),nptot(mxgrid),npt(3,mxgrid),
     &             igdata(5,mxgrid),igstat(mxgrid),ngdata(mxgrid),
     &             ngrid
c
c data calculation parameters
c
      common/dfcalc/cdata(5,mxcalc),ictype(mxcalc),icsect(mxcalc),
     &            icgrid(mxcalc),icgsec(mxcalc),ndata(mxcalc),
     &            icstat(mxcalc),icdata(5,mxcalc),ncalc
      common/dfclc2/cscal(mxcalc),iblkp(mxcalc)
c
c plot definition parameters
c
      common/dfplot/pdata(7,mxplot),iptype(mxplot),ipcalc(mxplot),
     &            ipcsec(mxplot),ipcont(mxplot),ncont(mxplot),
     &            ipdata(3,mxplot),nplot
c
c requests for restore of data from foreign dumpfiles
c
      common/dfrest/irestu(mxrest),irestb(mxrest),irests(mxrest),
     &               iresec(mxrest),nrest
c
c labels and titles
c
      common/cplot/zgridt(10,mxgrid),zgrid(mxgrid),
     &             zcalct(10,mxcalc),zcalc(mxcalc),
     &             zplott(10,mxplot),zplot(mxplot)
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
      integer mp, lwordp, lvec, ied3, newbas, lenb
      integer nbfnd, lenbad, lentrd, len3, ngrpm
      integer idena, nwroot, nosec
      logical noscf
      common /limy/ mp(10),lwordp,lvec,ied3,newbas,lenb,
     +              nbfnd,lenbad,lentrd,len3,ngrpm,
     +              idena,noscf,nwroot,nosec(mxroot)
c
      common/craypk/nbuf(9,maxorb),nbasin,norb,lenbas,ngrps,non,newb,
     &     ncol,jcol,ngrp(3),iblkg,iblka,iblkv
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common/miscop/qx(4,maxat),cx(3,maxorb),popocc(maxorb)
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

c - just to set d,e, and f
      parameter (mcen32=3*maxat-2)
      common/bufb/rbuf(72),ibuf(56),
     +            pa(18),cxi(3),cxj(3),b27,side,
     +            d(1),pp(1),ppp(mcen32),
     +            e(1),qq(1),qqq(mcen32),
     +            f(1),rr(1),rrr(mcen32),p1(200),p2(200),p3(200)
      character*10 charwall
      dimension q(*)
      data thresh/1.0d-8/
c
c driver for graphics routines
c
c  keys
c
c   igtype  1  2D
c           2  3D
c           3  sphere                       ; radial mesh
c           4  user defined                 ; iblk
c           5  interpolation          height;
c           6  isodensity             height;
c           7  atom positions
c           8  screen on scalar data value
c
c   ictype  1  density
c           2  mo                          ; mo #
c           3  atom difference             ; print, nconf, npop, iblkc
c           5  grad density
c           6  grad MO                     ; mo #
c           7  combination (previous  scale; unit,block,section
c                 + dumpfile*scale)
c           8  diff (not implemented)
c           9  vdw    
c           10 log vdw
c           11 vdw gradient
c           12 log vdw surface normal (same as 11)
c
c   iptype  1  lineprinter
c           2  contour
c           3  perspective (relief)
c
c   /route/
c   stype   1  grid generation
c           2  calculation
c           3  plot
c           4  restore grid
c           5  restore calculation
c
      write(iwr,6010)
c initialisation of ghost
      oinit = .false.
c density matrix for density and potential
      opotdm = .false.
c
c     icurgr = 0
      icurca = 0
      nav = lenwrd()
c
c     do we even need to drive off old properties modules?
c
      if (ocopy) then
c
c q(i1), czeta; q(i2) ctran (expanded out over orbitals)
       i1 = 1
       i2 = i1 + num*mxprms
c q(i3) - vectors q(i3b) - shell ordered vectors
       i3 = i2 + num*mxprms
c q(i3a) - scratch space for density grad rho and MO
c
       call graph1
       call graph2(q(i1),q(i2))
       if(newbas.ne.0)then
          call gridd(q(i3),ocopy)
       else
          write(iwr,*)
          write(iwr,*)' skipping vectors (nbas=0)'
          write(iwr,*)
          ncol = 0
          norb = 0
       endif
       i3b = i3 + ncol*norb
       i3a = i3b + ncol*norb
       if(i3a.gt.lwordp)goto 999
c
c  q(i4) - density matrix
       i4 = i3a + ncol*5
       i4a = i4 + norb*(norb + 1)/2 
       i5 = i4a + num
       if(i5.gt.lwordp)goto 999
c
      else
c
c     startup for density and potential grids
c     skip setup phase if we don't need it (allows limited
c     use of f- and g-functions for density, potentials etc)
c     this code is to minimise the amount of memory required
c     from the excessive 2 * num *mxprms + 5 * norb(norb+1/2
c     to 3 * norb(norb+1/2 )  [+ grid storage]
c     ONLY for density and potential at present
c
       i1 = 1
       i2 = i1
c     q(i3) - vectors 
       i3 = i2 
c
       iblkv = ibl7la
       norb = newbas
       call gridd(q(i3),ocopy)
       i3a = i3 + ncol*norb
       i3b = i3 
       if(i3a.gt.lwordp)goto 999
c
c  q(i4) - density matrix
       i4 = i3a + ncol*5
       i4a = i4 + norb*(norb + 1)/2 
       i5 = i4a + num
       if(i5.gt.lwordp)goto 999
c
      endif

      n = 1
      do i=1,maxat
         d(n)=1.0d0
         e(n)=1.0d0
         f(n)=1.0d0
         n=n+3
      enddo
c
c
c debugging listing
c
c      call listgr
c
c list job steps
c
      call listrt(iwr)
c
c start sequence at beginning
c
      call setstp(1)
c
c get the next task
c
 100  call nxtstp(ijob, i)
c
c  ***** ijob = 0 (completion) ********
c
      if(ijob.eq.0)then
         goto 900
c
c  ***** ijob = 1 (grid generation) ********
c
      else if(ijob.eq.1)then
         nmax = (lwordp-i5)/3
         if(igtype(i).eq.1.or.igtype(i).eq.2)then
            write(iwr,8000)' set up regular grid          ',
     &           zgrid(i),i,cpulft(1),charwall(),(zgridt(j,i),j=1,10)
c
c no work to do for regular grids, and no memory req
c
            ngdata(i)=0
         else if(igtype(i).eq.3)then
            write(iwr,8000)' generation of spherical grid  ',
     &           zgrid(i),i,cpulft(1),charwall(),(zgridt(j,i),j=1,10)
            call spgrid(q(i5),gdata(1,i),npt(1,i),
     &           geom(1,i),nptot(i),igdata(1,i),igdata(2,i),nmax,
     &           iwr)
            ngdata(i) = nptot(i)*3
         else if(igtype(i).eq.4)then
            write(iwr,8000)' setup of user defined grid    ',
     &           zgrid(i),i,cpulft(1),charwall(),(zgridt(j,i),j=1,10)
            if(nptot(i).gt.nmax)goto 999
c
c     retrieve data for user grid
c
            call rdedx(q(i5),nptot(i)*3,igdata(1,i),num8)
            ngdata(i) = nptot(i)*3
         else if(igtype(i).eq.5)then
            nmax = (lwordp-i5)/3
            write(iwr,8000)' generation of contour grid    ',
     &           zgrid(i),i,cpulft(1),charwall(),(zgridt(j,i),j=1,10)
c
c     generate by interpolation contouring of existing data grid
c     (at q(i6)) load the known data into q(i5)
c     results into q(i7) transfer to q(i5)
c
            nmax = (lwordp-i7)/3
            idata = icurca
            icube = icurgr
            if(idata.eq.0.or.icube.eq.0)call caserr
     &           ('must define dataset to contour')
            if(igtype(icube).ne.2)call caserr
     &           ('can only contour 3D data')
            if(icstat(idata).ne.2)call caserr('data not ready')
            call intrp3(q(i6),q(i7),
     &           npt(1,icube),npt(2,icube),npt(3,icube),geom(1,icube),
     &           ntot,gdata(1,i),igdata(1,i),iwr)
            call dcopy(ntot*3,q(i7),1,q(i5),1)
            nptot(i) = ntot
            ngdata(i) = nptot(i)*3
            i6 = i5 + ngdata(i)
            igstat(i)=2
            do 120 j = 1,3
               npt(j,i) = npt(j,icube)
               do 110 jj = 1,4
 110              geom((jj-1)*3+j,i)=geom((jj-1)*3+j,icube)
 120        continue
         else if(igtype(i).eq.6)then
            nmax = (lwordp-i5)/3
            write(iwr,8000)' generation of isodensity grid ',
     &           zgrid(i),i,cpulft(1),charwall(),(zgridt(j,i),j=1,10)
c     jk wrap - number of points unknown bring back
c     x,y,z, for each point in q(i5)
c
            icube = icurgr
            if(icube.eq.0.or.igtype(icube).ne.2)
     &           call caserr('need cubic grid for wrap')
            call contjk(q(i5),nptot(i),q(i3),q(i2),q(i1),popocc,
     &           q(i3a),gdata(1,i),nmax,iwr)
            write(iwr,*)nptot(i),' grid points generated by wrap'
            ngdata(i) = nptot(i)*3
            i6 = i5 + ngdata(i)
            igstat(i)=2
         else if(igtype(i).eq.7)then
            write(iwr,8000)' generation of atom coord grid',
     &           zgrid(i),i,cpulft(1),charwall(),(zgridt(j,i),j=1,10)
            do 122 iat = 1,nat
               q(i5+3*(iat - 1) + 0) = c(1,iat)
               q(i5+3*(iat - 1) + 1) = c(2,iat)
               q(i5+3*(iat - 1) + 2) = c(3,iat)
 122        continue
            nptot(i)=nat
            ngdata(i) = nptot(i)*3
         else if(igtype(i).eq.8)then
            if(icurca.eq.0)then
              call caserr('no data available for grid point selection')
            endif
c
c scratch memory - leave enough space to write
c the new grid of max length before starting scratch array
            i8 = i5 + nptot(icurgr)*3
c
c need also to leave space for the data
c
            if(i8.lt.i7)i8 = i7
            nmax = nav*(lwordp-i8)/3
            if (nmax.lt.nptot(icurgr))then
               call 
     &        caserr('no memory available for grid point selection')
            endif
            write(iwr,8000)' select points by data value',
     &           zgrid(i),i,cpulft(1),charwall(),(zgridt(j,i),j=1,10)
            call selgrd(q(i5),q(i6),nptot(icurgr),gdata(1,i),
     &           gdata(2,i), q(i8), nptot(i), igdata(1,i),iwr)
            ngdata(i) = nptot(i)*3
         else
            call caserr('unrecognised grid type in graph')
         endif
         i6 = i5 + ngdata(i)
         igstat(i)=2
         if(igsect(i).ne.0) call wgrsec(i,q(i5),iwr)
c
c set the current grid for subsequent calcs
c
         call setgrd(i,iwr)
         icurgr = i
         icurca = 0
c
c ****** ijob = 2 (data generation)
c
      else if(ijob.eq.2)then
c
c memory requirements
c    i5 - irregular grid
c    i6 - results
c    i7 - will be scratch array for difference calculations
c         and potential, and atomic scf
c
         write(iwr,8000)' calculation of data on grid   ',
     &        zcalc(i),i,cpulft(1),charwall(),(zcalct(j,i),j=1,10)
c
c  we assume that the grid parameters are in memory
c  but must retrieve any data into q(i5)
c
         if(iblkp(i).ne.-1)then
            call rdedx(q(i4a),num,iblkp(i),num8)
            itop=0
            do 124 ii = num,1,-1
               if(dabs(q(i4a+ii-1)).gt.thresh.and.itop.eq.0)itop=ii+5
 124        continue
            write(iwr,7010)
            write(iwr,7011)(ii,q(i4a+ii-1),ii=1,5*((itop+1)/5))
            opotdm = .false.
         else
            call dcopy(num,popocc,1,q(i4a),1)
            write(iwr,*)
            write(iwr,*)'MO populations taken from dumpfile'
            opotdm = .false.
         endif
c
         igrid = icurgr
         if(igrid.eq.0)call caserr('no grid defined')
         icgrid(i) = icurgr
         ndata(i)=nptot(igrid)
         nvec=1
         if(ictype(i).eq.5.or.
     &      ictype(i).eq.6.or.
     &      ictype(i).eq.11.or.
     &      ictype(i).eq.12)nvec=3
         ndata(i)=ndata(i)*nvec
         i7 = i6 + ndata(i)
         if(i7.gt.lwordp)goto 999
         if(ictype(i).eq.1)then
c ****    OLD code
c *           call denlp(q(i1),q(i2),q(i3),q(i4a),q(i3a),q(i5),q(i6),
c *    +                 q(i7),i7)
c
           if(.not.opotdm)then
cc              write(6,*) 'potdm', q(i3b), ncol, norb
              call potdm(q(i3b),q(i4a),q(i4),ncol,norb)
              opotdm = .true.
              write(iwr,'(a24,f15.2,a10)')
     1              ' density matrix done at ',cpulft(1),charwall()
            endif
c
            call getnpt(np,nnx,nny,nnx)
            nw = lwordp - i7
*****
            nsp = 4 * np
            npp = np
            npage = 1
1002        if (nsp - nw) 1000,1000,1001
1001        npage = npage + 1
            npp = npp / npage
            nsp = 4 * npp
            go to 1002
1000        i71 = i7 + npp
            i72 = i71+ npp
            i73 = i72+ npp
            last = i73 + npp
            if (last.gt.lwordp) then
             write(iwr,1003) last, lwordp
1003         format(/
     +       '*** memory required = ', i10/
     +       '*** memory available = ', i10/)
             call caserr2('insufficient memory for density')
            endif
*****
            ien = 0
c
c DEBUG     oroot = opg_root()
c           if (oroot) then
c            call prev(q(i3b),popocc,ncol,num,num)
c            call prtril(q(i4),num)
c           endif
            do iii = 1,npage
               ist = ien + 1
               if(iii.ne.npage)then
                  ien = ist + npp - 1
               else
                  ien = np
               endif
               write(iwr,3009)iii,ist,ien,ss,charwall()
 3009          format(' density batch ',3i15,' starts at ',
     1                f8.2,' secs ',a10,' wall')
               call eldeng(ist,ien,q(i5),q(i4),q(i7),q(i71),
     +         q(i72),q(i73),q(i6))
            enddo
         else if(ictype(i).eq.2)then
            call molp(q(i1),q(i2),q(i3),q(i3a),q(i5),q(i6),icdata(1,i))
         else if(ictype(i).eq.3)then
c
c  --- atom difference density
c
            iprin = icdata(1,i)
c
c  --- user defined orbital populations
c
            mconf = icdata(2,i)
            mpop  = icdata(3,i)
            iblkc = icdata(4,i)
            if(mconf.ne.0)then
               nw = (mconf-1)/nav + 1
               i7a = i7 + mpop
               i8 = i7a + nw
               call rdedx(q(i7),mpop,iblkc,num8)
               call reads(q(i7a),nw,num8)
            else
               i8 = i7
            endif
            nwd = lwordp - i8
            call atlp(q(i5),q(i6),q(i1),q(i2),q(i4a),q(i3a),q(i8),
     &           q(i3),mconf,q(i7a),q(i7),nwd,iprin)
         else if(ictype(i).eq.4)then
            call getnpt(np,nnx,nny,nnx)
            nw = lwordp - i7
c **** the factor m225 reflects max. shell type, 
c **** originally set to 225 for g shells
            m225 = nspchk()
            nsp = nw/m225
            if(nsp .eq. 0)then
               write(iwr,*)'memory shortage for potential',nw
               goto 999
            endif
            if(.not.opotdm)then
cc               write(6,*) 'potdm', q(i3b), ncol, norb
               call potdm(q(i3b),q(i4a),q(i4),ncol,norb)
               opotdm = .true.
               write(iwr,'(a19,f15.2,a10)')
     1               ' density matrix at ',cpulft(1),charwall()
            endif
            call dponc2(np,q(i5),q(i6),icdata(2,i))
            npage = 1 + np / nsp
            npp = np / npage
            ien = 0
            do iii = 1,npage
               ist = ien + 1
               if(iii.ne.npage)then
                  ien = ist + npp - 1
               else
                  ien = np
               endif
               write(iwr,*)'potential batch',iii,ist,ien
               call dpogr2(ist,ien,m225,q(i5),q(i4),q(i7),q(i6))
            enddo
         else if(ictype(i).eq.5)then
            call ddenlp(q(i1),q(i2),q(i3),q(i4a),q(i3a),q(i5),q(i6))
         else if(ictype(i).eq.6)then
            call dmolp(q(i1),q(i2),q(i3),q(i3a),q(i5),q(i6),icdata(1,i))
         else if (ictype(i).eq.7)then
c
c - comb   q(i7) for scratch
c
c  add to the current grid, make the result the current grid
c
            i8=i7+ndata(i)
            if(i8.gt.lwordp)goto 999
            ier = 0
            ic = 0
            call rcasec(icdata(3,i),ic,icdata(1,i),icdata(2,i),
     &           q(i7),.true.,ier)
            if(ndata(ic).ne.ndata(i))call caserr
     *           ('invalid comb - grids must be of equal size')
c
c basic grid checking
c
            ochk = .false.
            omism = .false.
            icg = icgrid(ic)
            ig = icgrid(i)
            if(icg.ne.0.and.ig.ne.0)then
               if(igtype(icg).ne.igtype(ig))omism = .true.
               if(.not.omism.and.igtype(ig).le.2)then
                  ochk = .true.
                  tol = 1.0d-5
                  if(dist(geom(1,icg),geom(1,ig)).gt.tol)omism = .true.
                  if(dist(geom(4,icg),geom(4,ig)).gt.tol)omism = .true.
                  if(dist(geom(7,icg),geom(7,ig)).gt.tol)omism = .true.
                  if(igtype(ig).eq.2)then
                     if(dist(geom(10,icg),geom(10,ig)).gt.tol)
     &                    omism = .true.
                  endif
               endif
            endif
            if(omism)then
               write(iwr,6040)
            else if(.not.ochk)then
               write(iwr,6045)
            endif

c perform daxpy on root node only
c other nodes hold null result

           if(opg_root())then
               call daxpy(ndata(i),cdata(1,i),q(i7),1,q(i6),1)
           else 
               call vclr(q(i6),1,ndata(i))
           endif
c           write(iwr,*)' diff contrib',zcalc(ic),cdata(1,i)
         else if (ictype(i).eq.8)then
c
c - diff   q(i7) for scratch
c
c use first dataset for reference
c at present diff is not available - comb can
c be used to add in one grid at a time
c


            i8=i7+ndata(i)
            if(i8.gt.lwordp)goto 999
            igrid = icgrid(icdata(1,i))
            icgrid(i) = igrid
            icgsec(i) = igsect(igrid)
            ndat = nptot(igrid)
            write(iwr,*)' diff grid',igrid,ndat
c
c everything from dumpfile at present
c
           call vclr(q(i6),1,ndat)
            do 130 j = 1,5
               ic = icdata(j,i)
               if(ic.ne.0)then
                  ier = 0
                  call rcasec(icsect(ic),ic,idaf,ibl3d,
     &                 q(i7),.true.,ier)
                  call daxpy(ndat,cdata(j,i),q(i7),1,q(i6),1)
                  write(iwr,*)' diff contrib',zcalc(ic),cdata(j,i)
               endif
 130        continue
         else if (ictype(i).ge.9.and.ictype(i).le.12)then
c
c  van der waals functions and surface normals
c


            olog = (ictype(i).eq.10.or.ictype(i).eq.12)
            ograd = (ictype(i).eq.11.or.ictype(i).eq.12)
            nvdw = icdata(1,i)
            iblkc = icdata(2,i)
            if(nvdw.ne.0)then
               nw = (nvdw-1)/nav + 1
               i7a = i7 + nvdw
               i8 = i7a + nw
               call rdedx(q(i7),nvdw,iblkc,num8)
               call reads(q(i7a),nw,num8)
            else
               i8 = i7
            endif
            if((lwordp - i8 - nat).le.0)goto 999
            call vdwf(q(i5),q(i6),nvdw,q(i7a),q(i7),q(i8),olog,
     +                ograd,iwr)
         else if(ictype(i).eq.-1)then
            call gridck(q(i5),iwr)
         endif
c
c  scaling
c
         if(dabs(cscal(i)-1.0d0).gt.1.0d-20)then
            write(iwr,7015)cscal(i)
            call dscal(ndata(i),cscal(i),q(i6),1)
         endif
c
c  grid analysis - only for regular grids, scalar data
c
         if(igtype(igrid).eq.1.and.nvec.eq.1)then
            call anagri(q(i6),npt(1,igrid),npt(2,igrid),0,iwr)
         else if(igtype(igrid).eq.2.and.nvec.eq.1)then
            call anagri(q(i6),npt(1,igrid),npt(2,igrid),
     &                        npt(3,igrid),iwr)
         else if(nvec.eq.1)then
            call avegri(q(i6),nptot(igrid),iwr)
         endif
c
      write(iwr,3003) cpulft(1),charwall()
 3003 format(2x,'end grid no analysis at ',f14.2,' secs ',a10,' wall')
c
c put data on dumpfile
c
         icstat(i) = 2
         if(icsect(i).ne.0)call wcasec(i,q(i6),iwr)
         icurca=i
c
c  === ijob = 3 ====== plot generation ==========
c
      else if(ijob.eq.3)then
         write(iwr,8000)' generation of plot            ',
     &        zplot(i),i,cpulft(1),charwall(),(zplott(j,i),j=1,10)
         ic = icurca
         ig = icurgr
         if(ic.eq.0)call caserr('no calculation defined')
c
c   get contour heights from ed7
c
         i8 = i7
         if(ipcont(i).ne.0)then
            i8 = i7 + ncont(i)
            if(i8.gt.lwordp)goto 999
            call rdedx(q(i7),ncont(i),ipcont(i),num8)
c            write(iwr,*)' contours',(q(i7+ii-1),ii=1,ncont(i))
         endif
         if(iptype(i).eq.1)then
            write (iwr,7020) (zplott(j,i),j=1,10)
            ax = npt(1,ig)-1
            ay = npt(2,ig)-1
            write (iwr,7030)'x  ',((geom(3+j,ig)-geom(j,ig))/ax,j=1,3),
     &           'y  ',((geom(6+j,ig)-geom(j,ig))/ay,j=1,3)
            call linep(q(i6),npt(1,ig),npt(2,ig),q(i7),ncont(i),iwr)
         else if(iptype(i).eq.2)then
c
c memory for cntour 4*n*n real*4 words
c
            if(igtype(ig).ne.1.or.npt(1,ig).ne.npt(2,ig))
     &           call caserr('unsuitable data for contour')
            i9 = i8 + 4*npt(1,ig)*npt(1,ig)/nav + 1
            if(i9.gt.lwordp)goto 999
c
            if(.not.oinit)call iniplt
            oinit = .true.
            write (iwr,7050) (zplott(j,i),j=1,10)
            ax = npt(1,ig)-1
            write (iwr,7030)'x  ',((geom(3+j,ig)-geom(j,ig))/ax,j=1,3),
     &           'y  ',((geom(6+j,ig)-geom(j,ig))/ax,j=1,3)
c.... call frame see if arguments are required
            call cntour(q(i6),q(i8),zplott(1,i),npt(1,ig),
     +           q(i7),ncont(i), ictype(ic),.false.)
         else if(iptype(i).eq.3)then
            if(igtype(ig).ne.1.or.npt(1,ig).ne.npt(2,ig))
     &           call caserr('unsuitable data for surface plot')
            if(.not.oinit)call iniplt
            oinit = .true.
            scamax=pdata(1,i)
            scamin=pdata(2,i)
            facmax=pdata(3,i)
            facmin=pdata(4,i)
            anga=pdata(5,i)
            angb=pdata(6,i)
            dist3d = pdata(7,i)
            write (iwr,7040) (zplott(j,i),j=1,10) , scamax , scamin ,
     +           anga , angb , dist3d
            dist3d = dist3d/dist(geom(1,ig),geom(4,ig))
            call arch3d(q(i6),q(i8),zplott(1,i),npt(1,ig),
     +           scamax,scamin,facmax,facmin,anga,angb,dist3d)
         else
            call caserr('requested plot option unavailable')
         endif
c
c restore grid from a foreign dumpfile
c
      else if (ijob.eq.4)then
         ier = 0
         call rgrsec(irests(i),ig,irestu(i),irestb(i),
     &        q(i5),.true.,ier)
c         write(iwr,*)' restore grid',i,ig
         ng = ngdata(ig)
         i6 = i5+ng
         if(i6.gt.lwordp)goto 999
         call setgrd(ig,iwr)
         icurgr = ig
         icurca = 0
         if(iresec(i).ne.0) call wgrsec(ig,q(i5),iwr)
      else if (ijob.eq.5)then
c
c restore data from foreign dumpfile - we will get the grid
c without its data
c
         ier = 0
         i6 = i5
         call rcasec(irests(i),ic,irestu(i),irestb(i),
     &        q(i6),.true.,ier)
c         write(iwr,*)' restore  calc',i,ic
         i7 = i6 + ndata(ic)
         if(i7.gt.lwordp)goto 999
         icurca = ic
         icurgr = icgrid(ic)
         call setgrd(icurgr,iwr)
         if(iresec(i).ne.0)call wcasec(ic,q(i6),iwr)
      else
         call caserr('unrecognised task type in graph')
      endif
      goto 100
 900  continue
      call secsum
      top = cpulft(1)
      write (iwr,6020) top,charwall()
c
      return
 999  write(iwr,6030)
      return
 6010 format(1x,104('=')//55x,'graphics module'//1x,104('='))
 6020 format (//1x,104('-')//' end of graphical analysis at ',f14.2,
     +        ' seconds ',a10,' wall')
 6030 format(//
     + 1x,'Error in memory allocation of graphical routines'/
     + 1x,'Not enough memory to continue'/)
 6040 format(1x,'** warning ** the grids for the two data sets',
     &     ' do not match',/,15x,
     &     'the resultant grid is probably invalid')

 6045 format(1x,'the correspondence of the combined grids has not ',
     &     'been checked')
 8000 format(//,1x,104('*'),/,1x,'* commence ',a30,1x,a8,
     &     ' (number ',i2,')',' at ',f11.2,' seconds',a10,' wall',
     &      3x,'*',/,1x,'* title :',10a8,14x,'*',/,
     &  1x,104('*'))
 7010 format(1x,'user defined MO occupations:',/,4x,5('  orb    occ '))
 7011 format(4x,5(i5,f8.3))
 7015 format(/,1x,'scaling data by ',f20.6)
 7020 format (/20x,'construct line printer contour plot'//20x,
     +        'title : ',10a8//' normalised plane definition',2x,
     +        '(for 1 character on line printer plot)'//9x,'x',13x,'y',
     +        13x,'z'/)
 7030 format (a3,3f14.7)
 7040 format (//20x,'construct  perspective plot'//20x,'title : ',
     +        10a8//20x,'scale plot to maxima : ',f10.4,' units'//20x,
     +        'scale plot to minima : ',f10.4,' units'///20x,
     +        'viewing angles  (',f6.1,' ,',f6.1,' ) degrees'///20x,
     +        'viewing distance ',f10.4,' bohr')
 7050 format (//20x,'construct contour plot'//20x,'title : ',
     +        10a8//' normalised plane definition',2x,
     +        '(for 1 cm on plot)'/9x,'x',13x,'y',13x,'z'/)
      end
**==nspchk.f
      function nspchk()
c
c     ----- check for s,p,d,f,g basis requirements -----
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
      odbas = .false.
      ofbas = .false.
      ogbas = .false.
      lendd = 4
      do i = 1 , nshell
         if (ktype(i).gt.2) then
            if (ktype(i).eq.3) odbas = .true.
            if (ktype(i).eq.4) ofbas = .true.
            if (ktype(i).eq.5) ogbas = .true.
         end if
      enddo
      if (odbas) lendd = 6
      if (ofbas) lendd = 10
      if (ogbas) lendd = 15
      nspchk = lendd * lendd
      return
      end
**==graph1.f
      subroutine graph1
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
      parameter (mxgr10=mxgrps*mxprms)
      dimension ngrp(3),minf(3),maxf(3),ioff(10)
      common/junk/eta(mxgr10),wnorm(mxgr10),vlist(maxat,4),
     *itype(mxgrps),itypp(mxgrps),
     *icent(mxgrps),mbase(mxgrps),ntrm(mxgrps),buff(2),ibuff(8),ktw,mon,
     *ibl(3),mgrps,ncentg(2),cgx(36),idhar,norbd,mmbase(mxgrps)
     *,ntype(maxorb)
     *,ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),cf(mxprim),cg(mxprim)
     *,abuf(maxat)
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
      common/junkc/zcomm(29),ztita(10),ztagga(maxat)
      common/miscop/
     *kstart(mxshel),katom(mxshel),ktype(mxshel),
     *kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),nshell,non3(3)
c
      integer mapsie
      common /linkan/ mapsie(maxorb)
c
      common/gjs/c1(mxprms),alp(mxprms)
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
      common/restri/jjfile(63),lds(508),isect(508),ldsect(508)
c
      data mone/1/
      data ngrp/1,3,6/
      data ioff/0,0,0,0,3,3,3,-3,-3,-3/
      data minf/1,2,5/
      data maxf/1,4,10/
      data pi32/5.56832799683170d0/
      data pi/3.14159265359d0/
      data pt5,pt75,dzero/0.5d0,0.75d0,0.0d0/
      data sqrt3/1.73205080756888d0/
c
      nav= lenwrd()
      mdump = mxgrps*(2*mxprms+6/nav) + maxat*4 + 38 + 18/nav
      call secget(isect(491),mone,iblkv)
      m110 = 10 + maxat
      m1420 = mxprim*5 + maxat

      call rdedx(ex,mxprim,iblkv,idaf)
      call reads(cs,m1420,idaf)
      call rdchrs(ztita,m110,idaf)
      call readis(kstart,mach(2)*nav,idaf)
c
      if (nat.le.0 .or. nat.gt.maxat) then
       call caserr('invalid number of atoms in analysis preprocessor')
      else if (nshell.lt.0 .or. nshell.gt.mxshel) then
       call caserr('invalid number of shells in analysis preprocessor')
      else if (num.lt.0 .or. num.gt.maxorb) then
       call caserr('invalid number of bfns in analysis preprocessor')
      else
         idhar = 0
         mon = nat
         ktw = 0
         fact = 2.0d0**0.25d0*pi**0.625d0
         do 20 i = 1 , nat
            vlist(i,1) = c(1,i)
            vlist(i,2) = c(2,i)
            vlist(i,3) = c(3,i)
            vlist(i,4) = czan(i)
 20      continue
c
c ..... 1st decide no. of s,p,d groups for reordering
c
         ns = 0
         np = 0
         nd = 0
         do 30 i = 1 , nshell
            j = kmax(i) - kmin(i)
            if (j.le.1) then
               ns = ns + 1
            else if (j.lt.3) then
               np = np + 1
            else if (j.eq.3) then
               ns = ns + 1
               np = np + 1
            else
               nd = nd + 1
            end if
 30      continue
         mgrps = ns + np + nd
         if (mgrps.lt.0 .or. mgrps.gt.mxgrps) then
            call caserr(
     +      'invalid number of groups in analysis preprocessor')
         end if
      end if
c
 40   nss = 0
      npp = 0
      ndd = 0
      ngrps = 0
      ospbas = .false.
c in case nshell=0
      kk = 0
      do 240 i = 1 , nshell
         l = katom(i)
         nt1 = kng(i)
         is = kstart(i)
         if = is + nt1 - 1
         mini = kmin(i)
         maxi = kmax(i)
         ngrps = ngrps + 1
         jj = maxi - mini
         if (jj.ne.0) then
            if (jj.lt.3) go to 130
            if (jj.eq.3) then
               ospbas = .true.
            else
               k = 0
               ndd = ndd + 1
               igrps = ns + np + ndd
               mbase(igrps) = ns + np*3 + (ndd-1)*6
               m = 3
               do 50 ig = is , if
                  k = k + 1
                  alp(k) = ex(ig)
                  ee = alp(k) + alp(k)
                  fac = pi32/(ee*dsqrt(ee))
                  fac = pt75*fac/(ee*ee)
                  c1(k) = cd(ig)/dsqrt(fac)
 50            continue
               facd = dzero
               do 70 ig = 1 , k
                  do 60 jg = 1 , ig
                     ee = alp(ig) + alp(jg)
                     fac = ee*dsqrt(ee)
                     dum = pt75*c1(ig)*c1(jg)/(ee**2*fac)
                     if (ig.ne.jg) then
                        dum = dum + dum
                     end if
                     facd = facd + dum
 60               continue
 70            continue
               do 80 ig = 1 , k
                  c1(ig) = c1(ig)/dsqrt(facd*pi32)
 80            continue
               go to 180
            end if
         end if
         m = 1
         nss = nss + 1
         igrps = nss
         mbase(igrps) = nss - 1
         k = 0
         do 90 ig = is , if
            k = k + 1
            alp(k) = ex(ig)
            ee = alp(k) + alp(k)
            fac = pi32/(ee*dsqrt(ee))
            c1(k) = cs(ig)/dsqrt(fac)
 90      continue
         facs = dzero
         do 110 ig = 1 , k
            do 100 jg = 1 , ig
               ee = alp(ig) + alp(jg)
               fac = ee*dsqrt(ee)
               dum = c1(ig)*c1(jg)/fac
               if (ig.ne.jg) then
                  dum = dum + dum
               end if
               facs = facs + dum
 100        continue
 110     continue
         do 120 ig = 1 , k
            c1(ig) = c1(ig)/dsqrt(facs*pi32)
 120     continue
         go to 180
 130     k = 0
         npp = npp + 1
         igrps = ns + npp
         mbase(igrps) = ns + (npp-1)*3
         m = 2
         do 140 ig = is , if
            k = k + 1
            alp(k) = ex(ig)
            ee = alp(k) + alp(k)
            fac = pi32/(ee*dsqrt(ee))
            fac = pt5*fac/ee
            c1(k) = cp(ig)/dsqrt(fac)
 140     continue
         facp = dzero
         do 160 ig = 1 , k
            do 150 jg = 1 , ig
               ee = alp(ig) + alp(jg) 
               fac = ee*dsqrt(ee)
               dum = pt5*c1(ig)*c1(jg)/(ee*fac)
               if (ig.ne.jg) then
                  dum = dum + dum
               end if
               facp = facp + dum
 150        continue
 160     continue
         do 170 ig = 1 , k
            c1(ig) = c1(ig)/dsqrt(facp*pi32)
 170     continue
 180     nn = ngrp(m)
         kk = (igrps-1)*mxprms
         icent(igrps) = l
         itype(igrps) = m
         ntrm(igrps) = nt1
         mmbase(igrps) = mbase(igrps)
         do 220 j = 1 , nn
            ktw = ktw + 1
            go to (190,200,210) , m
 190        ntype(ktw) = j
            go to 220
 200        ntype(ktw) = j + 1
            go to 220
 210        ntype(ktw) = j + 4
 220     continue
         factd = fact
         if (m.eq.3) factd = fact*sqrt3
         do 230 k = 1 , nt1
            kk = kk + 1
            wnorm(kk) = c1(k)*factd
            eta(kk) = alp(k)
 230     continue
         if (ospbas) then
            ngrps = ngrps + 1
            ospbas = .false.
            go to 130
         end if
 240  continue
c
      if (num.ne.ktw) then
       call caserr('parameter error detected in analysis preprocessor')
      else if (ngrps.ne.mgrps) then
       call caserr('parameter error detected in analysis preprocessor')
      else if (kk.lt.0 .or. kk.gt.mxgaus) then
       call caserr('parameter error detected in analysis preprocessor')
      else
         norbd = ktw
         mgrps = ngrps
c
c ----- now setup gamess/atmol mapping vectors
c
         ij = 0
         do 260 i = 1 , 3
            m1 = minf(i)
            m2 = maxf(i)
            do 250 j = 1 , norbd
               jj = ntype(j)
               if (jj.ge.m1 .and. jj.le.m2) then
                  ij = ij + 1
                  ik = ij + ioff(jj)
                  mapsie(ik) = j
               end if
 250        continue
 260     continue
c
         call wrt3(eta,mdump,ibl7la,num8)
c
         return
      end if
      end
**==graph2.f
      subroutine graph2(czeta,ctran)
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      parameter (mxgr10=mxprms*mxgrps)
      dimension lx(10),ly(10),lz(10),ilif(3),xs(3)
      common/miscop/p(maxat),q(maxat),r(maxat),chg(maxat),
     *cx(maxorb),cy(maxorb),cz(maxorb),pop(maxorb)
      common/craypk/nterm(maxorb),icentr(maxorb),itype(maxorb),
     *mx(maxorb),my(maxorb),mz(maxorb),nx(maxorb),ny(maxorb),nz(maxorb),
     *nbasin,norb,lenbas,ngrps,non,newb,ncol,jcol,ngrp(3),
     *iblkg,iblka,iblkv
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
      common/junk/czet(mxgr10),anorm(mxgr10),
     *ccx(maxat),ccy(maxat),ccz(maxat),czz(maxat),ityp(mxgrps),
     *itypp(mxgrps),
     *icent(mxgrps),mbase(mxgrps),ntrm(mxgrps),title(2),ititl(8),
     *nbasis,ncent,ibl(3),mgrps,ncentg(2),cgx(36),idddd,norbd,
     *mmbase(mxgrps)
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
      integer mp, lwordp, lvec, ied3, newbas, lenb
      integer nbfnd, lenbad, lentrd, len3, ngrpm
      integer idena, nwroot, nosec
      logical noscf
      common /limy/ mp(10),lwordp,lvec,ied3,newbas,lenb,
     +              nbfnd,lenbad,lentrd,len3,ngrpm,
     +              idena,noscf,nwroot,nosec(mxroot)
c
      common/gjs/pr(mxprms),alp(mxprms),c1(mxprms)
c
      dimension czeta(*),ctran(*)
c
      data m0,m1,mm1/0,1,-1/
      data ilif/0,1,4/
      data lx,ly,lz/
     *0,1,0,0,1,1,0,2,0,0,
     *0,0,1,0,1,0,1,0,2,0,
     *0,0,0,1,0,1,1,0,0,2/
      data xs/'s','p','d'/
      data pi/3.14159265359d0/
c
      ngrp(1) = 1
      ngrp(2) = 3
      ngrp(3) = 6
      nav = lenwrd()
      ndump = mxgrps*(2*mxprms+6/nav) + maxat*4 + 38 + 18/nav
      pip = 1.0d0/(2.0d0**0.25d0*pi**0.625d0)
      t = pip/dsqrt(3.0d0)
      do 20 i = 1 , 3
         pr(i) = pip
         pr(i+3) = t
 20   continue
      call rdedx(czet,ndump,ibl7la,num8)
      iblkv = iposun(num8)
      if (ncent.gt.m0 .and. ncent.le.maxat) go to 40
 30   call caserr('invalid parameters detected in analysis preprocessor'
     +            )
 40   if (mgrps.lt.m0 .or. mgrps.gt.mxgrps) go to 30
      if (nbasis.lt.m0 .or. nbasis.gt.maxorb .or. nbasis.ne.nbfnd)
     +    go to 30
      nbasin = nbfnd
      ngrps = mgrps

      non = ncent
      do 50 i = 1 , non
         p(i) = ccx(i)
         q(i) = ccy(i)
         r(i) = ccz(i)
         chg(i) = czz(i)
 50   continue
      norb = 0
      do 90 i = 1 , ngrps
         m = ityp(i)
         nn = ngrp(m)
         l = icent(i)
         nt1 = ntrm(i)
         ll = (i+mm1)*mxprms
         do 60 k = 1 , nt1
            alp(k) = czet(k+ll)
            c1(k) = anorm(k+ll)
 60      continue
         il = ilif(m)
         top = p(l)
         bot = q(l)
         param = r(l)
         ktw = l*3 - 2
         do 80 j = 1 , nn
            norb = norb + 1
            itype(norb) = m
            nterm(norb) = nt1
            icentr(norb) = l
            cx(norb) = top
            cy(norb) = bot
            cz(norb) = param
            il = il + 1
            k = lx(il)
            nx(norb) = k
            mx(norb) = k + ktw
            k = ly(il)
            ny(norb) = k
            my(norb) = k + ktw
            k = lz(il)
            nz(norb) = k
            mz(norb) = k + ktw
            kk = (norb+mm1)*mxprms
            at = pr(j)
            do 70 k = 1 , nt1
               czeta(k+kk) = alp(k)
               ctran(k+kk) = c1(k)*at
 70         continue
 80      continue
 90   continue
      if (norb.ne.nbasis) go to 30
      mbasin= norb
      lenbas = norb*(norb+1)/2
      If (nbasin.ne.mbasin) go to 30
c
      if (oprint(53)) then
         write (iwr,6010)
         do 100 i = 1 , non
            write (iwr,6020) p(i) , q(i) , r(i) , chg(i) , zaname(i)
 100     continue
         write (iwr,6030) non
         if (oprint(53)) write (iwr,6040)
      end if
      ngr = 0
      i = 1
 110  if (i.gt.norb) then
         if (oprint(53)) write (iwr,6050) ngr , norb
         if (ngr.ne.ngrps) go to 30
c
         return
      else
         n = itype(i)
         l = icentr(i)
         nn = nterm(i)
         jj = ngrp(n)
         ngr = ngr + 1
         if (oprint(53)) write (iwr,6060) ngr , xs(n) , zaname(l) , nn
         if (l.lt.m1 .or. l.gt.non) go to 30
         k = (i+mm1)*mxprms
         do 150 j = 1 , nn
            dum = czeta(k+j)
            dum2 = dum + dum
            c = (pi/dum2)**1.5d0
            go to (140,120,130) , n
 120        c = c/(dum2+dum2)
            go to 140
 130        c = c/(4.0d0*dum2*dum2)
 140        c = ctran(j+k)*dsqrt(c)
            if (oprint(53)) write (iwr,6020) c , dum
            if (dum.lt.1.0d-8) go to 30
            if (dabs(c).lt.1.0d-8) go to 30
 150     continue
         i = i + jj
         go to 110
      end if
 6010 format (//28x,'geometry'/28x,8('-')///8x,'x',15x,'y',15x,'z',13x,
     +        'charge',5x,'label'//)
 6020 format (4f16.7,1x,a8)
 6030 format (/' no. of nuclei =',i3)
 6040 format (//8x,'ordered list of groups (normalised)'/8x,35('-')/)
 6050 format (//' no. of groups=',4x,i4/' no. of gtos=',6x,i4/)
 6060 format (/' group type centre  nterm'/1x,25('-')/1x,i5,1x,a1,4x,a8,
     +        i5//11x,'ctran',12x,'zeta'/11x,21('-'))
      end
**==grdnuc.f
      subroutine grdnuc(noc,vlist,c,dint,nocmx)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      dimension dint(6),vlist(nocmx,4),c(3)
      common/miscop/gnuc(6)
c
      do 20 i = 1 , 6
         gnuc(i) = 0.0d0
 20   continue
      do 50 i = 1 , noc
         dst = (vlist(i,1)-c(1))**2 + (vlist(i,2)-c(2))
     +         **2 + (vlist(i,3)-c(3))**2
         if (dst.gt.1.0d-10) then
            mu = 0
            do 40 j = 1 , 3
               do 30 k = 1 , j
                  mu = mu + 1
                  gnuc(mu) = gnuc(mu) + vlist(i,4)*(vlist(i,j)-c(j))
     +                       *(vlist(i,k)-c(k))/(dsqrt(dst))**5
 30            continue
 40         continue
         end if
 50   continue
      dst = gnuc(1) + gnuc(3) + gnuc(6)
      dint(1) = -(3.0d0*gnuc(1)-dst)
      dint(2) = -(3.0d0*gnuc(3)-dst)
      dint(3) = -(3.0d0*gnuc(6)-dst)
      dint(4) = -3.0d0*gnuc(2)
      dint(5) = -3.0d0*gnuc(4)
      dint(6) = -3.0d0*gnuc(5)
c
      return
      end
**==gridck.f
      subroutine gridck(q,iwr)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension q(*)
      call getnpt(np,nx,ny,nz)
      write(iwr,*)' ---- grid check ----',np,nx,ny,nz
      do 10 i = 1,np
         call getpt(i,ax,ay,az,q,1,idum)
         write(iwr,100)i,ax,ay,az
 10   continue
 100  format(1x,i5,3f10.3)
      write(iwr,*)'---- end check ----'
      return
      end
**==gridd.f
      subroutine gridd (vect,ocopy)
c
c     writes vectors to ed7
c     ocopy drives specification of reordered vectors for old modules
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
      common/junk/eig(maxorb)
c
      integer mapsie
      common /linkan/ mapsie(maxorb)
c
c
      real*8 bias, ajm, ajm1, ajm2
      integer nvar, mope, ncs, nos, npopul
c
      common/datplt/bias(8),nvar(14),mope(2),ajm(110),ajm1(110),
     + ajm2(110),ncs(90),nos(90),npopul(90)
      common/blkin/bbias(8),mvar(14),nope(2),ajmm(330),nncs(270),
     *popm(maxorb),popa(maxorb)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common/miscop/q(4,maxat),cx(3,maxorb),popocc(maxorb)
      common/craypk/nterm(9,maxorb),nbasin,norb,lenbas,ngrps,non,
     &newbas,ncol,jcol,ngrp(3),iblkg,iblka,iblkv
c
      logical oming, oextg, odirec, ovcfre, osemi, ostopm, olevd
      logical osym_op,osym_dens
      integer mina, minb, mouta, moutb, lock, ibrk
      integer numg, isecg, iblk3g, mextra, nsymo, iorbsy, isunor
      real*8  gapa1, gapa2, gapb1, gapb2, scaleg, symtag
      common/atmol3/mina,minb,mouta,moutb,lock,ibrk,
     + numg(2),isecg(2),iblk3g(2),mextra,nsymo,iorbsy(100),
     + oming, oextg,
     + gapa1,gapa2,gapb1,gapb2,scaleg,symtag(9),odirec(50),
     + isunor,ovcfre,osemi,ostopm,olevd,osym_op,osym_dens
c
c
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
c
      real*8 ctran
      integer ilifc, ntran, itran
      logical otran,otri
      common/tran/ilifc(maxorb),ntran(maxorb),itran(mxorb3),
     +            ctran(mxorb3),otran,otri
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      common/junkc/zspace(29),zcont(100),z1(15),zname(15),
     * zbuff(10),zcomm(19),zzbuff(10)
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
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
      dimension vect(*)
c
      data itvec/3/
c
      nav = lenwrd()
      nw0 = 338 + 286/nav
c
      if (mina.le.0 .or. mina.gt.350)
     +    call caserr('invalid section specified for eigenvectors')
c     restore vectors for analysis
      call secget(mina,itvec,iblvec)
      call getqp(zcomm,zzbuff,eig,popocc,nbas,newbas,ncol,ieig,idiff,
     +           maxorb,iblvec)
      write (iwr,6010) mina , ibl3d , yed(idaf) , (zcomm(7-i),i=1,6) ,
     +                 zzbuff
      nbsq = ncol*nbas
      nbsq2= nbas*nbas

      call rdedx(vect,nbsq2,iblvec,idaf)

      call tdown(vect,ilifq,vect,ilifq,ncol)

      if(ocopy) then
c save vectors before reordering
       call dcopy(ncol*norb,vect,1,vect(ncol*norb+1),1)
      endif

      if (oprint(53)) then
         if (.not.otran) then
            write (iwr,6020)
            do 30 i = 1 , newbas
               n = ntran(i)
               l = ilifc(i) + 1
               write (iwr,6030) ctran(l) , itran(l) , n , i
               if (n.ne.1) then
                  do 20 j = 2 , n
                     l = l + 1
                     write (iwr,6030) ctran(l) , itran(l)
 20               continue
               end if
 30         continue
         end if
         write (iwr,6040) zcomm(5)
         call prsql(vect,ncol,nbas,nbas)
         if (idiff.gt.0) then
            write (iwr,6050) mina
            write (iwr,6060) (popocc(i),i=1,ncol)
         end if
      end if
      if (ocopy) then
       ii = 0
       do 50 i = 1 , ncol
          do 40 j = 1 , nbas
             eig(j) = vect(ii+mapsie(j))
 40       continue
          call dcopy(nbas,eig,1,vect(ii+1),1)
          ii = ii + nbas
 50    continue
      endif

      call wrt3(vect,nbsq,iblkv,num8)
      iblkg = iposun(num8)
      ibl7la = iblkg
c
c     now migrate /datplt/ common block data
c
      call dcopy(nw0,bias,1,bbias,1)
      mvar(11) = 1
c
      return
 6010 format (//1x,'vectors restored from section',i4,
     +        ' of dumpfile starting at block',i6,' of ',
     +        a4//' header block information : '/1x,27('-')
     +        /' vectors created under acct. ',a8/1x,a7,
     +        'vectors created by ',a8,' program at ',a8,' on ',a8,
     +        ' in the job ',a8/' with the title : ',10a8/)
 6020 format (/40x,12('*')/40x,'list of sabf'/40x,12('*')///11x,
     +        'transformed orbitals (normalised)'//
     +        ' coefficient old orbital nterm new orbital'//)
 6030 format (1x,f11.7,i12,i6,i12)
 6040 format (//1x,104('-')//40x,16('*')/40x,a8,' vectors'/40x,16('*')
     +        //)
 6050 format (//40x,'orbital populations  (from vectors section ',i3,
     +        ')'//)
 6060 format (/15x,7f15.7)
      end
**==hamil.f
      subroutine hamil(p,q,dt,dos,pcap,qcap,s,t,u,fc,fo)
c
c    one configuration
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
      dimension p(*),q(*),dt(*),dos(*),pcap(*),qcap(*),
     *s(*),t(*),u(*),fc(*),fo(*)
      common/junk/pqn(99),expa(99),ctran(mxprms,30),ajmn(12),czn,
     *expmn,c(15,15),itran(mxprms,30),nrow(30),nbcon(4),nbas(4),
     *nsym,nosh(4),ncsh(4),nccup(4),ibase(33),
     *popul(4),facto(30),fatt(30),thre(15),a(15,15),eivr(15,15),
     *cm2scf(15,15),cm1scf(15,15),eps(15),gmin(15,15),smin(15,15),
     *populs(4),term,energ,pot,cin,vir,thrscf,enorm,
     *red,expdif,encor,
     *nsize(10),nsavf(4),n1(6),nbvar(30),nbas1,nbast,nbvar1,
     *nblock,numvar,nrun,nitscf,nitsc1,nitsc2,nconv,nalarm,no,nstep,
     *iiii,nsht,nosht,nflag3,nbasm,mvar,iex,newbas,nit,nbpos(8)
      data dzero/0.0d0/
c
c    compute in the following order the matrices qcap,pcap,smin,qmin,
c    fc,fo
c    qcap mean the q matrix of roothaan
c    pcap mean the p matrix of roothaan
c
c    compute qcap and pcap
c
      nit = 0
      do 20 i = 1 , nsym
         nbas1 = nbas(i)
         nstep1 = nbas1*(nbas1+1)/2
         nit = nit + nstep1
         n1(i) = nstep1
 20   continue
      nstep1 = 0
      k = 1
      do 60 i = 1 , nsym
c     i=lambda
         n11 = n1(i)
         do 50 m = 1 , n11
c     m=pq
            fact1 = dzero
            fact2 = dzero
            kk = 1
            mm = m + nstep1
            nstep2 = 0
            do 40 j = 1 , nsym
c     j=mu
               n12 = n1(j)
               do 30 n = 1 , n12
c     n=rs
                  nn = n + nstep2
                  l = max(mm,nn)*(max(mm,nn)-1)/2 + min(mm,nn)
                  fact1 = fact1 + p(l)*dt(kk)
                  fact2 = fact2 + q(l)*dos(kk)
                  kk = kk + 1
 30            continue
               nstep2 = nstep2 + n12
 40         continue
            pcap(k) = fact1
            qcap(k) = fact2
            k = k + 1
 50      continue
         nstep1 = nstep1 + n11
 60   continue
c
c    computes smin and gmin
c
      nstep1 = 0
      nstep2 = 0
      do 100 i = 1 , nsym
         nsh = ncsh(i) + nosh(i)
         nbas1 = nbas(i)
         do 90 j = 1 , nsh
            j1 = nstep2 + j
            do 80 m = 1 , nbas1
               fact1 = dzero
               fact2 = dzero
               do 70 n = 1 , nbas1
                  k = nstep1 + max(m,n)*(max(m,n)-1)/2 + min(m,n)
                  fact1 = fact1 + s(k)*c(n,j1)
                  fact2 = fact2 + qcap(k)*c(n,j1)
 70            continue
               smin(m,j1) = fact1
               gmin(m,j1) = fact2
 80         continue
 90      continue
         nstep2 = nstep2 + nsh
         nstep1 = nstep1 + n1(i)
 100  continue
c
c    compute fc and fo
c
      k = 1
      nstep = 0
      do 140 i = 1 , nsym
         nosh1 = nosh(i)
         nsh = ncsh(i) + nosh1
         fact2 = popul(i)
         fact1 = populs(i) - fact2
         nbas1 = nbas(i)
         do 130 m = 1 , nbas1
            do 120 n = 1 , m
               term1 = dzero
               term2 = dzero
               do 110 j = 1 , nsh
                  j1 = j + nstep
                  if (j.ne.nsh) then
                     term1 = term1 + smin(m,j1)*gmin(n,j1) + gmin(m,j1)
     +                       *smin(n,j1)
                  else if (nosh1.ne.0) then
                     term2 = smin(m,j1)*gmin(n,j1) + smin(n,j1)
     +                       *gmin(m,j1)
                  else
                     term1 = term1 + smin(m,j1)*gmin(n,j1) + gmin(m,j1)
     +                       *smin(n,j1)
                  end if
 110           continue
               hppc = pcap(k) + t(k) - czn*u(k)
               fc(k) = hppc + term2*fact2/fact1
               fo(k) = hppc - qcap(k) + term1*populs(i)/fact1
               k = k + 1
 120        continue
 130     continue
         nstep = nstep + nsh
 140  continue
c
      return
      end
**==hdiag.f
      subroutine hdiag(h,u,p,iq,n,nbmx,iegen,nr)
c
c      ...   jacobi .... (?)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      dimension h(nbmx,nbmx),u(nbmx,nbmx),p(*),iq(*)
c...   common for overruling local criteria (e.g. dependency, diag thresholds)
      integer  i_depen_check, i_global_diag, n_depen_check
      logical  o_depen_print
      common/overrule/ i_depen_check,i_global_diag, o_depen_print, 
     1                 n_depen_check
c
c     i_depen_check : criterion for throwing out vectors
c     n_depen_check : number  of vectors to be thrown out 
c     o_depen_print : if true +> extra print 

c
      if (iegen.eq.0) then
         do 20 i = 1 , n
            call vclr(u(1,i),1,n)
            u(i,i) = 1.00d0
 20      continue
      end if
      nr = 0
      if (n.le.1) return
      nmi1 = n - 1
      do 40 i = 1 , nmi1
         p(i) = 0.0d0
         ipl1 = i + 1
         do 30 j = ipl1 , n
            if (p(i).le.dabs(h(i,j))) then
               p(i) = dabs(h(i,j))
               iq(i) = j
            end if
 30      continue
 40   continue
      rap = 5.0d-8
      hdtest = 1.0d38
 50   pmax = 0.0d0
      do 60 i = 1 , nmi1
         if (i.gt.1) then
            if (pmax.ge.p(i)) go to 60
         end if
         pmax = p(i)
         ipiv = i
         jpiv = iq(i)
 60   continue
      if (pmax.le.0.0d0) return
      if (hdtest.gt.0.0d0) then
         if (pmax.gt.hdtest) go to 80
      end if
      hdimin = dabs(h(1,1))
      do 70 i = 2 , n
         if (hdimin.gt.dabs(h(i,i))) hdimin = dabs(h(i,i))
 70   continue
      hdtest = hdimin*rap
c
      if (i_global_diag.ne.-999) hdtest = 10.0d0**(-i_global_diag*1.0d0)
c
      if (hdtest.ge.pmax) return
 80   nr = nr + 1
      tang = dsign(2.0d0,(h(ipiv,ipiv)-h(jpiv,jpiv)))*h(ipiv,jpiv)
     +   /(dabs(h(ipiv,ipiv)-h(jpiv,jpiv))
     +   +dsqrt((h(ipiv,ipiv)-h(jpiv,jpiv))**2+4.0d0*h(ipiv,jpiv)**2))
      cosine = 1.0d0/dsqrt(1.0d0+tang**2)
      sine = tang*cosine
      hii = h(ipiv,ipiv)
      h(ipiv,ipiv) = cosine**2*(hii+tang*(2.0d0*h(ipiv,jpiv)+tang*h(jpiv
     +               ,jpiv)))
      h(jpiv,jpiv) = cosine**2*(h(jpiv,jpiv)-tang*(2.0d0*h(ipiv,jpiv)-
     +               tang*hii))
      h(ipiv,jpiv) = 0.0d0
      if (h(ipiv,ipiv).lt.h(jpiv,jpiv)) then
         htemp = h(ipiv,ipiv)
         h(ipiv,ipiv) = h(jpiv,jpiv)
         h(jpiv,jpiv) = htemp
         htemp = dsign(1.0d0,-sine)*cosine
         cosine = dabs(sine)
         sine = htemp
      end if
      do 100 i = 1 , nmi1
         if (i.ne.ipiv .and. i.ne.jpiv) then
            if (iq(i).eq.ipiv .or. iq(i).eq.jpiv) then
               k = iq(i)
               htemp = h(i,k)
               h(i,k) = 0.0d0
               ipl1 = i + 1
               p(i) = 0.0d0
               do 90 j = ipl1 , n
                  if (p(i).le.dabs(h(i,j))) then
                     p(i) = dabs(h(i,j))
                     iq(i) = j
                  end if
 90            continue
               h(i,k) = htemp
            end if
         end if
 100  continue
      p(ipiv) = 0.0d0
      p(jpiv) = 0.0d0
      do 110 i = 1 , n
         if (i.lt.ipiv) then
            htemp = h(i,ipiv)
            h(i,ipiv) = cosine*htemp + sine*h(i,jpiv)
            if (p(i).lt.dabs(h(i,ipiv))) then
               p(i) = dabs(h(i,ipiv))
               iq(i) = ipiv
            end if
            h(i,jpiv) = -sine*htemp + cosine*h(i,jpiv)
            if (p(i).lt.dabs(h(i,jpiv))) then
               p(i) = dabs(h(i,jpiv))
               iq(i) = jpiv
            end if
         else if (i.ne.ipiv) then
            if (i.lt.jpiv) then
               htemp = h(ipiv,i)
               h(ipiv,i) = cosine*htemp + sine*h(i,jpiv)
               if (p(ipiv).lt.dabs(h(ipiv,i))) then
                  p(ipiv) = dabs(h(ipiv,i))
                  iq(ipiv) = i
               end if
               h(i,jpiv) = -sine*htemp + cosine*h(i,jpiv)
               if (p(i).lt.dabs(h(i,jpiv))) then
                  p(i) = dabs(h(i,jpiv))
                  iq(i) = jpiv
               end if
            else if (i.ne.jpiv) then
               htemp = h(ipiv,i)
               h(ipiv,i) = cosine*htemp + sine*h(jpiv,i)
               if (p(ipiv).lt.dabs(h(ipiv,i))) then
                  p(ipiv) = dabs(h(ipiv,i))
                  iq(ipiv) = i
               end if
               h(jpiv,i) = -sine*htemp + cosine*h(jpiv,i)
               if (p(jpiv).le.dabs(h(jpiv,i))) then
                  p(jpiv) = dabs(h(jpiv,i))
                  iq(jpiv) = i
               end if
            end if
         end if
 110  continue
      if (iegen.eq.0) then
           call drot(n,u(1,ipiv),1,u(1,jpiv),1,cosine,sine)
      end if
      go to 50
c
      end
**==hexnuc.f
      subroutine hexnuc(non,vlist,c,dint,nocmx)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension dint(10),vlist(nocmx,4),c(3)
      common/miscop/t(15),a(3)
c...
c...compute nuclear contributions to hexadecapole moments.
c...p4z=(35zzzz-30zzrr+3rrrr)/8,xxrr,yyrr,zzrr,rrrr,xyrr,xzrr,yzrr.
c...p4x=(35xxxx-30xxrr+3rrrr)/8,p4y=(35yyyy-30yyrr+3rrrr)/8,
c...
      call vclr(t,1,15)
c
      do 60 nc = 1 , non
         a(1) = vlist(nc,1)
         a(2) = vlist(nc,2)
         a(3) = vlist(nc,3)
         chg = vlist(nc,4)
         mu = 0
         do 50 i = 1 , 3
            do 40 j = 1 , i
               do 30 k = 1 , j
                  do 20 l = 1 , k
                     mu = mu + 1
                     t(mu) = t(mu) + chg*(a(i)-c(i))*(a(j)-c(j))
     +                       *(a(k)-c(k))*(a(l)-c(l))
 20               continue
 30            continue
 40         continue
 50      continue
 60   continue
c...
      pprr = t(1) + t(3) + t(10)
      qqrr = t(3) + t(5) + t(12)
      ssrr = t(10) + t(12) + t(15)
      rrrr = pprr + qqrr + ssrr
      dint(1) = (35.0d0*t(1)-30.0d0*pprr+3.0d0*rrrr)/8.0d0
      dint(2) = (35.0d0*t(5)-30.0d0*qqrr+3.0d0*rrrr)/8.0d0
      dint(3) = (35.0d0*t(15)-30.0d0*ssrr+3.0d0*rrrr)/8.0d0
      dint(4) = pprr
      dint(5) = qqrr
      dint(6) = ssrr
      dint(7) = rrrr
      dint(8) = t(2) + t(4) + t(11)
      dint(9) = t(6) + t(8) + t(13)
      dint(10) = t(7) + t(9) + t(14)
c...
      return
      end
**==intrp3.f
      subroutine intrp3(data,grid,nx,ny,nz,geom,np,value,
     +                  iflag,iwr)
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
      dimension iof1(12),jof1(12),kof1(12),iof2(12),jof2(12),kof2(12)
      dimension data(nx,ny,nz),grid(3,*),geom(3,4)
c
      data iof1/0,0,0,0,0,0,1,1,0,0,1,1/
      data jof1/0,0,1,1,0,0,0,0,0,1,1,0/
      data kof1/0,1,1,0,0,1,1,0,0,0,0,0/
c
      data iof2/1,1,1,1,0,0,1,1,0,0,1,1/
      data jof2/0,0,1,1,1,1,1,1,0,1,1,0/
      data kof2/0,1,1,0,0,1,1,0,1,1,1,1/
c
      np = 0
c
c step lengths
c
      ap = nx - 1
      stepxx = (geom(1,2) - geom(1,1))/ap
      stepyx = (geom(2,2) - geom(2,1))/ap
      stepzx = (geom(3,2) - geom(3,1))/ap
      ap = ny - 1
      stepxy = (geom(1,3) - geom(1,1))/ap
      stepyy = (geom(2,3) - geom(2,1))/ap
      stepzy = (geom(3,3) - geom(3,1))/ap
      ap = nz - 1
      stepxz = (geom(1,4) - geom(1,1))/ap
      stepyz = (geom(2,4) - geom(2,1))/ap
      stepzz = (geom(3,4) - geom(3,1))/ap
c
c      write(iwr,*)' interp3d',nx,ny,nz,value
      write(iwr,1000)value
      nwarn = 0
      do 130 k = 1,nz-1
         do 120 j = 1,ny-1
            do 110 i = 1,nx - 1
c
c determine whether a point is in the current cell
c
               icount = 0
               tx = 0.0d0
               ty = 0.0d0
               tz = 0.0d0
               do 100 ii = 1,12
                  val1 = data(i+iof1(ii),j+jof1(ii),k+kof1(ii))
                  val2 = data(i+iof2(ii),j+jof2(ii),k+kof2(ii))
                  if((val1 - value)*(val2-value).lt.0.0d0)then
                     df = val2 - val1
                     icount = icount + 1
                     if(dabs(df).lt.1.0d-20)then
                        t = 0.5d0
                     else
                        t = (value - val1) / df
                     endif

                     gx1 = (i+iof1(ii)-1)*stepxx +
     &                     (j+jof1(ii)-1)*stepxy +
     &                     (k+kof1(ii)-1)*stepxz

                     gy1 = (i+iof1(ii)-1)*stepyx +
     &                     (j+jof1(ii)-1)*stepyy +
     &                     (k+kof1(ii)-1)*stepyz

                     gz1 = (i+iof1(ii)-1)*stepzx +
     &                     (j+jof1(ii)-1)*stepzy +
     &                     (k+kof1(ii)-1)*stepzz

                     dgx = (iof2(ii)-iof1(ii))*stepxx +
     &                     (jof2(ii)-jof1(ii))*stepxy +
     &                     (kof2(ii)-kof1(ii))*stepxz

                     dgy = (iof2(ii)-iof1(ii))*stepyx +
     &                     (jof2(ii)-jof1(ii))*stepyy +
     &                     (kof2(ii)-kof1(ii))*stepyz

                     dgz = (iof2(ii)-iof1(ii))*stepzx +
     &                     (jof2(ii)-jof1(ii))*stepzy +
     &                     (kof2(ii)-kof1(ii))*stepzz
c
                     tx = tx + gx1 + dgx * t
                     ty = ty + gy1 + dgy * t
                     tz = tz + gz1 + dgz * t
                  endif
 100           continue
c
c  require warning that contouring occupies the edge cells.
c
               if(icount.ne.0)then
                  acount = icount
                  np = np + 1
                  grid(1,np) = geom(1,1) + tx / acount
                  grid(2,np) = geom(2,1) + ty / acount
                  grid(3,np) = geom(3,1) + tz / acount
                  if ( (k .eq. 1) .or. (k .eq. (nz-1)) .or.
     &                 (j .eq. 1) .or. (j .eq. (ny-1)) .or.
     &                 (i .eq. 1) .or. (i .eq. (nx-1)) ) then
                     nwarn = nwarn + 1
                  endif
               endif
 110        continue
 120     continue
 130  continue
c      
      if(nwarn .ne. 0)then
         if(iflag.eq.0)then
            write(iwr,1020)'*ERROR*',nwarn
            write(iwr,1030)
            call caserr('grid contour size error')
         else
            write(iwr,1020)'WARNING',nwarn
         endif
      endif
      write(iwr,1010)np
c
 1000 format(/,1x,'commence contouring of 3D field at level',f20.6)
 1010 format(/,1x,'contouring has generated ',i7,' surface points')
 1020 format(//,1x,a7,' -',i4,
     &     ' points generated are in edge cells ',
     &     /,1x,' check that the 3D grid you are contouring ',
     &     'is large enough')
 1030 format(/,1x,' add the keyword nocheck to the type cont ',
     &     'directive if you want the job to ',/,1x,
     &     'continue past this check.')
      return
      end
**==jacvao.f
      subroutine jacvao(aa,n,eivu,eeivr)
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
      dimension aa(*),eivu(*),eeivr(*)
      common/junk/pqn(99),expa(99),ctran(mxprms,30),ajmn(12),czn,
     *expmn,cvec(15,15),itran(mxprms,30),nrow(30),nbcon(4),nbas(4),
     *nsym,nosh(4),ncsh(4),nccup(4),ibase(33),
     *popul(4),facto(30),fatt(30),thre(15),a(15,15),eivr(15,15),
     *cm2scf(15,15),cm1scf(15,15),eps(15),gmin(15,15),smin(15,15),
     *populs(4),term,energ,pot,cin,vir,thrscf,enorm,
     *red,expdif,encor,
     *nsize(10),nsavf(4),n1(6),nbvar(30),nbas1,nbast,nbvar1,
     *nblock,numvar,nrun,nitscf,nitsc1,nitsc2,nconv,nalarm,no,nstep,
     *iiii,nsht,nosht,nflag3,nbasm,mvar,iex,newbas,nit,nbpos(8)
      data dzero,done,d55,d9,d10,sq7,quart,two/
     *0.0d0,1.0d0,.55d0,1.0d-9,1.0d-10,.707106781d0,.25d0,2.0d0/
      data half/0.5d0/
c
      do 30 i = 1 , n
         i2 = (i-1)*n
         do 20 j = 1 , n
            a(j,i) = aa(i2+j)
 20      continue
 30   continue
      if (n.le.1) then
         eeivr(1) = done
         eivu(1) = a(1,1)
         return
      else
         do 50 j = 1 , n
            do 40 i = 1 , n
               eivr(i,j) = dzero
 40         continue
            eivr(j,j) = done
 50      continue
c     find absolutely largest element of a
         atop = dzero
         do 70 i = 1 , n
            do 60 j = 1 , n
               if (atop.lt.dabs(a(i,j))) then
                  atop = dabs(a(i,j))
               end if
 60         continue
 70      continue
         if (atop.le.0) then
            return
         else
c     calculate the stopping criterion-- dstop
            avgf = dble(n*(n-1))*d55
            d = dzero
            do 90 jj = 2 , n
               do 80 ii = 2 , jj
                  s = a(ii-1,jj)/atop
                  d = s*s + d
 80            continue
 90         continue
            dstop = d9*d
c
c      calculate the threshold, thresh
            thrsh = dsqrt(d/avgf)*atop
c
c     start a sweep
c
 100        iflag = 0
            do 160 jcol = 2 , n
               jcol1 = jcol - 1
               do 150 irow = 1 , jcol1
                  aij = a(irow,jcol)
c
c    compare   the off diagonal element with thrsh
c
                  if (dabs(aij).gt.thrsh) then
                     aii = a(irow,irow)
                     ajj = a(jcol,jcol)
                     s = ajj - aii
c
c    check to see if the chosen rotation is less than the rounding error
c    if so then do not rotate
c
                     if (dabs(aij).gt.d9*dabs(s)) then
                        iflag = 1
c
c    if the routation is very close yo 45 degrees set sin and cos
c  to 1/root 2).
c
                        if (d10*dabs(aij).lt.dabs(s)) then
c
c      calculation of sin and cos for rotation that is not very close to
c            45 degrees
c
                           t = aij/s
                           s = quart/dsqrt(quart+t*t)
c
c     cas=c  sin=s
c
                           c = dsqrt(s+half)
                           s = two*t*s/c
                        else
                           s = sq7
                           c = s
                        end if
c
c   calculation of the new elements of matrix a
c
                        call drot(irow,a(1,jcol),1,a(1,irow),1,c,s)
                        i2 = irow + 2
                        if (i2.le.jcol) then
                           do 110 i = i2 , jcol
                              t = a(i-1,jcol)
                              u = a(irow,i-1)
                              a(i-1,jcol) = s*u + c*t
                              a(irow,i-1) = c*u - s*t
 110                       continue
                        end if
                        a(jcol,jcol) = s*aij + c*ajj
                        a(irow,irow) = c*a(irow,irow) - s*(c*aij-s*ajj)
                        do 120 j = jcol , n
                           t = a(irow,j)
                           u = a(jcol,j)
                           a(irow,j) = c*t - s*u
                           a(jcol,j) = s*t + c*u
 120                    continue
c
c   rotation completed
c
                          call drot(n,eivr(1,jcol),1,eivr(1,irow),1,c,s)
c
c   calculate   the new norm d and compare with dstop
c
                        s = aij/atop
                        d = d - s*s
                        if (d.lt.dstop) then
c
c    recalculate dstop and thres to discard rounding errors
c
                           d = dzero
                           do 140 jj = 2 , n
                              do 130 ii = 2 , jj
                                 s = a(ii-1,jj)/atop
                                 d = s*s + d
 130                          continue
 140                       continue
                           dstop = d9*d
                        end if
                        thrsh = dsqrt(d/avgf)*atop
                     end if
                  end if
 150           continue
 160        continue
            if (iflag.ne.0) go to 100
c
c    place eigenvectors in eivu
            do 170 j = 1 , n
               eivu(j) = a(j,j)
 170        continue
            do 190 i = 1 , n
               i2 = (i-1)*n
               do 180 j = 1 , n
                  eeivr(i2+j) = eivr(j,i)
 180           continue
 190        continue
c
            return
         end if
      end if
      end
**==jdaaa.f
      subroutine jdaaa(dt,v,s,fo,fc)
c
c    this subroutine diagonalises the hartree fock matrix fc and fo
c   symmetry by symmetry, nbas (n) is smaller or equal to 15
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
      dimension dt(*),v(*),s(*),fo(*),fc(*),w(225),epsi(15)
      common/junk/pqn(99),expa(99),ctran(mxprms,30),ajmn(12),czn,
     *expmn,c(15,15),itran(mxprms,30),nrow(30),nbcon(4),nbas(4),
     *nsym,nosh(4),ncsh(4),nccup(4),ibase(33),
     *popul(4),facto(30),fatt(30),thre(15),a(15,15),eivr(15,15),
     *cm2scf(15,15),cm1scf(15,15),eps(15),gmin(15,15),smin(15,15),
     *populs(4),term,energ,pot,cin,vir,thrscf,enorm,
     *red,expdif,encor,
     *nsize(10),nsavf(4),n1(6),nbvar(30),nbas1,nbast,nbvar1,
     *nblock,numvar,nrun,nitscf,nitsc1,nitsc2,nconv,nalarm,no,nstep,
     *iiii,nsht,nosht,nflag3,nbasm,mvar,iex,newbas,nit,nbpos(8)
      equivalence (w(1),gmin(1,1)),(epsi(1),smin(1,15))
      data done,dzero,d40/
     *1.0d0,0.0d0,1.0d39/
c
c**** move s to square form at dt
c
      m1 = 0
      m = 0
      nstep1 = 1
      do 170 n = 1 , nsym
         nbas1 = nbas(n)
         do 30 i = 1 , nbas1
            nstep = nbas1*(i-1)
            do 20 j = 1 , i
               ij = nstep + j
               ji = nbas1*(j-1) + i
               m = m + 1
               term = s(m)
               dt(ij) = term
               dt(ji) = term
 20         continue
 30      continue
         call jacvao(dt,nbas1,epsi,v)
         do 50 i = 1 , nbas1
            nstep = nbas1*(i-1)
            term = done/dsqrt(epsi(i))
            do 40 j = 1 , nbas1
               ji = nstep + j
               v(ji) = v(ji)*term
 40         continue
 50      continue
c
c**** v contains vt/sqrt(so) (so=eigenvalue matrix)
c
         noshic = nosh(n)
         ncshic = ncsh(n)
         if (noshic.eq.0) then
            ifin = 1
         else if (ncshic.eq.0) then
            ifin = 1
         else
            ifin = 2
         end if
c
c**** form dt=v*f*vt  (f is equal to fc of fo)
c
         do 160 ish = 1 , ifin
            do 90 i = 1 , nbas1
               nstep = nbas1*(i-1)
               do 80 j = 1 , i
                  ji = nstep + j
                  ik = nbas1*(j-1)
                  ij = ik + i
                  term = dzero
                  do 70 k = 1 , nbas1
                     do 60 l = 1 , nbas1
                        kl = m1 + max(k,l)*(max(k,l)-1)/2 + min(k,l)
                        lj = ik + l
                        ki = nstep + k
                        if (ifin.eq.1) then
                           if (noshic.gt.0) then
c    this  is  the  case  of  only an  open  shell of a  symmetry
                              iloop = 1
                              aa = fo(kl)
                              term = term + v(ki)*v(lj)*aa
                              go to 60
                           end if
                        else if (ish.ge.ifin) then
c   this  is  the  case of one open shell which is not the lowest
c     of  that  symmetry
                           iloop = ncshic + 1
                           aa = fo(kl)
                           term = term + v(ki)*v(lj)*aa
                           go to 60
                        end if
                        iloop = ncshic
                        aa = fc(kl)
                        term = term + v(ki)*v(lj)*aa
 60                  continue
 70               continue
                  dt(ij) = term
                  dt(ji) = term
 80            continue
 90         continue
c
c**** diagonalises dt and transformation for vectors dt=vt*u
c
            call jacvao(dt,nbas1,epsi,w)
            do 120 i = 1 , nbas1
               do 110 j = 1 , nbas1
                  term = dzero
                  nstep = nbas1*(j-1)
                  ji = nstep + i
                  do 100 k = 1 , nbas1
                     ik = nbas1*(k-1) + i
                     kj = nstep + k
                     term = term + v(ik)*w(kj)
 100              continue
                  dt(ji) = term
 110           continue
 120        continue
c
c**** order eigenvectors and eigenvalues whithin symmetry
c
            do 150 i = 1 , iloop
               if ((ish.eq.1) .and. (i.gt.ncshic) .and. (ncshic.ne.0))
     +             go to 150
               mmin = 1
               if (nbas1.ne.1) then
                  do 130 j = 2 , nbas1
                     if (epsi(mmin).gt.epsi(j)) then
                        mmin = j
                     end if
 130              continue
               end if
               if (ish.ne.1) then
                  if (i.lt.iloop) then
                     epsi(mmin) = d40
                     go to 150
                  end if
               end if
               eps(nstep1) = epsi(mmin)
               epsi(mmin) = d40
               do 140 j = 1 , nbas1
                  jj = nbas1*(mmin-1) + j
                  c(j,nstep1) = dt(jj)
 140           continue
               nstep1 = nstep1 + 1
               epsi(mmin) = d40
 150        continue
 160     continue
         m1 = n1(n) + m1
 170  continue
c
      return
      end
**==jnput.f
      subroutine jnput(iwr)
c
c     one configuration
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
      common/junk/pqn(99),expa(99),ctran(mxprms,30),ajmn(12),czn,
     *expmn,c(15,15),itran(mxprms,30),nrow(30),nbcon(4),nbas(4),
     *nsym,nosh(4),ncsh(4),nccup(4),ibase(33),
     *popul(4),facto(30),fatt(30),thre(15),a(15,15),eivr(15,15),
     *cm2scf(15,15),cm1scf(15,15),eps(15),gmin(15,15),smin(15,15),
     *populs(4),term,energ,pot,cin,vir,thrscf,enorm,
     *red,expdif,encor,
     *nsize(10),nsavf(4),n1(6),nbvar(30),nbas1,nbast,nbvar1,
     *nblock,numvar,nrun,nitscf,nitsc1,nitsc2,nconv,nalarm,no,nstep,
     *iiii,nsht,nosht,nflag3,nbasm,mvar,iex,newbas,nit,nbpos(8)
      common/blkin/bias,dgath,scfat,dmnxp,rmned,rmxed,con1,con2,
     *nvar,ndiag,nxtrp,mxxtrp,norbit,nprint,mxvar,nextra,ifcont,irun,
     *nflag1,nflag2,nsvd,lim
c
      if (nprint.ne.0) then
c
         if (nsvd.eq.1) then
            write (iwr,6010)
         else
            write (iwr,6020)
         end if
c
         write (iwr,6030) czn
         write (iwr,6040) (nbas(i),i=1,nsym)
         write (iwr,6050) (ncsh(i),i=1,nsym)
         write (iwr,6060) (nccup(i),i=1,nsym)
c
         write (iwr,6070)
         write (iwr,6080) (ajmn(i),i=1,11)
      end if
c
      nbast = 0
      nosht = 0
      ncsht = 0
      do 20 i = 1 , nsym
         nbas1 = nbas(i)
         nsavf(i) = nbas1
         nosht = nosht + nosh(i)
         ncsht = ncsht + ncsh(i)
         nbast = nbast + nbas1
 20   continue
      nsht = ncsht + nosht
      bias = 5.0d0
      expmn = expa(1)
      do 30 i = 1 , 4
         popul(i) = nccup(i)
 30   continue
      dmnxp = dabs(expa(2)-expmn)
      do 50 i = 2 , nbast
         dumx = expa(i)
         if (expmn.gt.dumx) expmn = dumx
         j = i - 1
         do 40 k = 1 , j
            dumy = dabs(dumx-expa(k))
            if (dumy.lt.dmnxp) dmnxp = dumy
 40      continue
 50   continue
c
      if (nprint.ne.0) then
         write (iwr,6090) bias , dgath
         write (iwr,6100) nxtrp , scfat
         write (iwr,6110) ndiag , expmn
         write (iwr,6120) mxxtrp , dmnxp
         write (iwr,6130) rmned , rmxed
      end if
c
      if ((expmn.lt.1.d-3) .or. (expmn.gt.99999.d0))
     +    call caserr('invalid exponent in atomic scf calculation')
      if (nprint.eq.0) go to 150
      write (iwr,6140)
c
      nbas1 = nbas(1)
      if (nsym.gt.1) then
c
         nbas2 = nbas(2)
         if (nsym.gt.2) then
c
            nbas3 = nbas(3)
            if (nsym.gt.3) then
               nbas4 = nbas(4)
               do 60 i = 1 , nbas4
                  j = nbas1 + i
                  k = j + nbas2
                  l = k + nbas3
                  write (iwr,6150) pqn(i) , expa(i) , pqn(j) , expa(j) ,
     +                             pqn(k) , expa(k) , pqn(l) , expa(l)
 60            continue
               if (nbas4.eq.nbas3) go to 140
               nba4 = nbas4 + 1
               go to 120
            else
               nba4 = 1
               go to 120
            end if
         else
            nba3 = 1
         end if
      else
         do 70 i = 1 , nbas1
            write (iwr,6150) pqn(i) , expa(i)
 70      continue
         go to 150
      end if
 80   do 90 i = nba3 , nbas2
         j = nbas1 + i
         write (iwr,6150) pqn(i) , expa(i) , pqn(j) , expa(j)
 90   continue
 100  if (nbas1.ne.nbas2) then
         nba2 = nbas2 + 1
         do 110 i = nba2 , nbas1
            write (iwr,6150) pqn(i) , expa(i)
 110     continue
      end if
      go to 150
 120  do 130 i = nba4 , nbas3
         j = nbas1 + i
         k = nbas1 + nbas2 + i
         write (iwr,6150) pqn(i) , expa(i) , pqn(j) , expa(j) , pqn(k) ,
     +                    expa(k)
 130  continue
 140  if (nbas3.eq.nbas2) go to 100
      nba3 = nbas3 + 1
      go to 80
 150  newbas = 0
      do 160 i = 1 , nsym
         newbas = newbas + nbcon(i)
 160  continue
      if (nprint.ne.0) then
         write (iwr,6160)
         do 170 i = 1 , newbas
c*****maximum of 10 old per new
            nbf = nrow(i)
            write (iwr,6170) i , (ctran(l,i),l=1,nbf)
 170     continue
      end if
      nextra = 0
      do 180 i = 1 , nsht
         call vclr(c(1,i),1,15)
 180  continue
      nstep = 1
      nbasm = 0
      do 190 i = 1 , nsym
         nbas1 = nbcon(i)
         if (nbasm.lt.nbas1) then
            nbasm = nbas1
         end if
 190  continue
c
      return
 6010 format (6x,'jacobi diagonalization')
 6020 format (6x,'single vector diagonalization')
 6030 format (6x,'charge = ',f10.6)
 6040 format (6x,'symmetry species',11x,'s',5x,'p',5x,'d',5x,'f'/6x,
     +        'number of basis functions=',4(i2,4x))
 6050 format (6x,'number of closed shells  = ',i1,5x,i1,5x,i1,5x,i1)
 6060 format (6x,'open shell occupation  =   ',i1,5x,i1,5x,i1,5x,i1)
 6070 format (6x,'vector coupling coefficients k'//)
 6080 format (6f16.8)
 6090 format (6x,'diagonalisation scf threshold bias = ',e10.3,7x,
     +        '  diagonalisation divergence = ',e10.3)
 6100 format (6x,'number of scf iterations            = ',i2,16x,
     +        'scf divergence             = ',7x,e10.3)
 6110 format (6x,'number of diagonalization iterations= ',i2,16x,
     +        'minimum exponent           = ',7x,e10.3)
 6120 format (6x,'maximum number of iterations        = ',i2,16x,
     +        'minimum exponent difference= ',7x,e10.3)
 6130 format (6x,'relative minimum energy difference  = ',e10.3,8x,
     +        'relative maximum energy difference= ',e10.3)
 6140 format (6x,
     +  'basis functions (principal quantum numbers, orbital exponents)'
     +  ,//15x,'s',28x,'p',28x,'d',28x,'f'//)
 6150 format (4(5x,f4.1,5x,f15.6))
 6160 format ('0',5x,'contraction coefficients'/' ')
 6170 format (6x,i2,10f11.4)
      end
**==linep.f
      subroutine linep(grid,nx,ny,val,ncont,iwr)
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
      dimension xplotp(26),xplotn(15),iy(200),grid(*),val(*)
      common/junk/val1(200),val2(200),cy1(200),cy2(200),
     *ilifv(200)
      common/junkc/zspace(159),zbuff(10),
     * x2(200),xplot1(200),xplot2(200),xplot3(50),xplot4(200),
     * xchar(50)
      equivalence(iy(1),cy1(1))
      data m5,mxsize/5,125/
      data xplotp/
     *'a','b','c','d','e','f','g','h','i','j','k','l',
     *'m','n','o','p','q','r','s','t','u','v','w','x','y','z'/
      data xplotn/
     *'0','1','2','3','4','5','6','7','8','9','*','-','+','=','$'/
      data m0,m1,m15,m27,dzero,two,xspace,xpoint,xpoin/
     *0,1,15,26,0.0d0,2.0d0,' ','*','.'/
      if(nx.gt.mxsize)then
         write (iwr,6090)
         return
      endif

      i = m0
      j = m0
      do 20 k = 1 , ncont
         if (val(k).le.dzero) then
            j = j + m1
            if (j.le.m15) then
               xchar(k) = xplotn(j)
            end if
         else
            i = i + m1
            if (i.le.m27) then
               xchar(k) = xplotp(i)
            end if
         end if
 20   continue
      s1 = dzero
      do 30 i = 1 , ncont
         if (val(i).lt.s1)s1 = val(i)
 30   continue
      ten = -10.0d0 + s1
      k = ncont - m1
      do 50 i = 1 , k
         l = i + m1
         do 40 m = l , ncont
            if (val(m).gt.val(i)) then
               s1 = val(m)
               val(m) = val(i)
               val(i) = s1
               x1 = xchar(m)
               xchar(m) = xchar(i)
               xchar(i) = x1
            end if
 40      continue
 50   continue
      write (iwr,6010)
      write (iwr,6020) (xchar(i) , val(i),i=1,ncont)
c to index into grid array
      do 70 i = 1 , ny
         ilifv(i) = (i-1)*nx
 70   continue
c
c set up labels
c
      k = 5
      l = 0
      do 90 i = 1 , nx
         if (k.gt.nx) go to 100
         iy(i) = k
         x2(i) = xpoin
         l = i
         k = k + m5
 90   continue
 100  write (iwr,6030) (iy(i),i=1,l)
      write (iwr,6040) (x2(i),i=1,l)

      do 110 i = 1 , nx+2
         xplot4(i) = xpoint
         xplot2(i) = xspace
         val2(i) = ten
 110  continue
      write (iwr,6050) (xplot4(i),i=2,nx+2)
      i1 = ilifv(m1)
      call dcopy(nx,grid(i1+1),1,cy2(2),1)
c fill in end points
      cy2(1) = two*cy2(2)-cy2(3)
      cy2(nx+2) = two*cy2(nx+1)-cy2(nx)
c
      do 230 i = 2 , ny+1
         i1 = ilifv(i)
         do 120 j = 1 , nx+2
            xplot1(j) = xplot2(j)
            xplot2(j) = xspace
            val1(j) = val2(j)
            val2(j) = ten
            cy1(j) = cy2(j)
 120     continue
         if(i.eq.ny+1)then
            do 121 j = 2,nx+1
               i1 = ilifv(ny)
               i2 = ilifv(ny - 1)
               cy2(j) = two*grid(i1+j-1)-grid(i2+j-1)
 121        continue
         else
         call dcopy(nx,grid(i1+1),1,cy2(2),1)
         endif
c fill in end points
         cy2(1) = two*cy2(2)-cy2(3)
         cy2(nx+2) = two*cy2(nx+1)-cy2(nx)
         do 170 j = 1 , nx + 1
            r1 = cy1(j)
            s1 = val1(j)
            k = j + m1
            r2 = cy1(k)
            s2 = val1(k)
            do 160 l = 1 , ncont
               v = val(l)
               v1 = r1 - v
               v2 = r2 - v
               if (v1.lt.0) then
                  if (v2.lt.0) go to 160
                  if (v2.eq.0) go to 140
                  go to 150
               else if (v1.ne.0) then
                  if (v2.lt.0) go to 150
                  if (v2.eq.0) go to 140
                  go to 170
               end if
 130           if (s1.lt.v) then
                  xplot1(j) = xchar(l)
                  val1(j) = v
               end if
               go to 170
 140           if (s2.lt.v) then
                  xplot1(k) = xchar(l)
                  val1(k) = v
               end if
               go to 170
 150           if (dabs(v1).gt.dabs(v2)) go to 140
               go to 130
 160        continue
 170     continue
         do 220 j = 1 , nx + 1
            r1 = cy1(j)
            r2 = cy2(j)
            s1 = val1(j)
            do 210 l = 1 , ncont
               v = val(l)
               v1 = r1 - v
               v2 = r2 - v
               if (v1.lt.0) then
                  if (v2.lt.0) go to 210
                  if (v2.eq.0) go to 190
                  go to 200
               else if (v1.ne.0) then
                  if (v2.lt.0) go to 200
                  if (v2.eq.0) go to 190
                  go to 220
               end if
 180           if (s1.lt.v) then
                  xplot1(j) = xchar(l)
                  val1(j) = v
               end if
               go to 220
 190           xplot2(j) = xchar(l)
               val2(j) = v
               go to 220
 200           if (dabs(v1).gt.dabs(v2)) go to 190
               go to 180
 210        continue
 220     continue
         xplot1(nx+2) = xpoint
         write (iwr,6060) i-1 , (xplot1(j),j=2,nx+2)
 230  continue
      write (iwr,6050) (xplot4(i),i=2,nx+1)
c
      return
 6010 format (//14x,'legend'/14x,6('-')//5x,'character',5x,'value'/5x,
     +        '---------',5x,'-----'/)
 6020 format (9x,a1,f16.7)
 6030 format ('1'/6x,25i5)
 6040 format (6x,25(4x,a1))
 6050 format (5x,'*',124a1)
 6060 format (i5,'*',124a1)
c 6070 format (6x,62a1)
c 6080 format (/5x,10a8)
 6090 format (//5x,'*** no line printer plot generated'/
     +        '*** plot size exceeds page width'//)
      end
**==listrt.f
      subroutine listrt(iwr)
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
c grid definition parameters
c
      common/dfgrid/geom(12,mxgrid),gdata(5,mxgrid),igtype(mxgrid),
     &             igsect(mxgrid),nptot(mxgrid),npt(3,mxgrid),
     &             igdata(5,mxgrid),igstat(mxgrid),ngdata(mxgrid),
     &             ngrid
c
c data calculation parameters
c
      common/dfcalc/cdata(5,mxcalc),ictype(mxcalc),icsect(mxcalc),
     &            icgrid(mxcalc),icgsec(mxcalc),ndata(mxcalc),
     &            icstat(mxcalc),icdata(5,mxcalc),ncalc
      common/dfclc2/cscal(mxcalc),iblkp(mxcalc)
c
c plot definition parameters
c
      common/dfplot/pdata(7,mxplot),iptype(mxplot),ipcalc(mxplot),
     &            ipcsec(mxplot),ipcont(mxplot),ncont(mxplot),
     &            ipdata(3,mxplot),nplot
c
c requests for restore of data from foreign dumpfiles
c
      common/dfrest/irestu(mxrest),irestb(mxrest),irests(mxrest),
     &               iresec(mxrest),nrest
c
c labels and titles
c
      common/cplot/zgridt(10,mxgrid),zgrid(mxgrid),
     &             zcalct(10,mxcalc),zcalc(mxcalc),
     &             zplott(10,mxplot),zplot(mxplot)
c the job sequence
      integer stype(mxstp), arg(mxstp)
      common/route/stype,arg,nstep,istep
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
      character *12 caltyp, ctype
      write(iwr,1000)
      do 10 i = 1,nstep
         if(stype(i).eq.1)then
            write(iwr,1010)arg(i),(zgridt(j,arg(i)),j=1,9)
            if(igtype(arg(i)).eq.1)then
               write(iwr,1020)(npt(j,arg(i)),j=1,2)
               write(iwr,1030)(geom(j,arg(i)),j=1,9)
            else if(igtype(arg(i)).eq.2)then
               write(iwr,1025)(npt(j,arg(i)),j=1,3)
               write(iwr,1030)(geom(j,arg(i)),j=1,12)
            else if(igtype(arg(i)).eq.3)then
               write(iwr,1040)gdata(1,arg(i)),(geom(j,arg(i)),j=1,3),
     &              npt(1,arg(i))
            else if(igtype(arg(i)).eq.4)then
               write(iwr,1050)
            else if(igtype(arg(i)).eq.5)then
               write(iwr,1060)gdata(1,arg(i))
            else if(igtype(arg(i)).eq.6)then
               write(iwr,1070)gdata(1,arg(i))
            else if(igtype(arg(i)).eq.7)then
               write(iwr,1080)
            endif
            if(igsect(arg(i)).ne.0)write(iwr,1100)igsect(arg(i))
         else if(stype(i).eq.2)then
            write(iwr,2010)arg(i),(zcalct(j,arg(i)),j=1,9)
            ctype = caltyp(ictype(arg(i)),icdata(1,arg(i)))
            write(iwr,2020)ctype
            if(ictype(arg(i)).eq.7)write(iwr,2030)icdata(3,arg(i)),
     &           yed(icdata(1,arg(i))),icdata(2,arg(i)),
     &           cdata(1,arg(i))
            if(dabs(cscal(arg(i))-1.0d0).gt.1.0d-20)
     &           write(iwr,2040)cscal(arg(i))
            if(icsect(arg(i)).ne.0)write(iwr,1100)icsect(arg(i))
         else if(stype(i).eq.3)then
            write(iwr,3010)arg(i),(zplott(j,arg(i)),j=1,9)
            if(iptype(arg(i)).eq.1)then
               write(iwr,3020)'line printer plot'
            else if(iptype(arg(i)).eq.2)then
               write(iwr,3020)'contour plot     '
            else if(iptype(arg(i)).eq.3)then
               write(iwr,3020)'relief plot      '
            endif
         else if(stype(i).eq.4)then
            write(iwr,4010)'grid',irests(arg(i)),yed(irestu(arg(i))),
     &           irestb(arg(i))
            if(iresec(arg(i)).ne.0)write(iwr,1100)iresec(arg(i))
         else if(stype(i).eq.5)then
            write(iwr,4010)'data',irests(arg(i)),yed(irestu(arg(i))),
     &           irestb(arg(i))
            if(iresec(arg(i)).ne.0)write(iwr,1100)iresec(arg(i))
         endif
         write(iwr,1005)
 10   continue
      write(iwr,9000)
 1000 format(/1x,104('*')/1x,'*',41x,'job steps requested ',41x,'*'
     &/1x,104('*'))
 1005 format(1x,104('*'))
 1010 format(1x,'* generate grid # ',i3,' title: ',9a8,2x,'*')
 1020 format(1x,'* grid is 2D, size ',i3,' * ',i3,75x,'*')
 1025 format(1x,'* grid is 3D, size ',i3,' * ',i3,' * ',i3,69x,'*')
 1030 format(1x,'* grid corners : origin ', 3f8.3,55x,'*',/,
     &       1x,'*',21x,'x ',3f8.3,55x,'*',/,
     &       1x,'*',21x,'y ',3f8.3,55x,'*':/,
     &       1x,'*',21x,'z ',3f8.3,55x,'*')
 1040 format(1x,'* grid is spherical, radius ',f8.4,' centred at ',
     &     3f8.3,' angular mesh ',i4,27x,'*')
 1050 format(1x,'* grid taken from user input file',94x,'*')
 1060 format(1x,'* grid on isovalue surface at level ',f10.5,81x,'*')
 1070 format(1x,'* grid on isodensity surface at level ',f10.5,79x,'*')
 1080 format(1x,'* grid at atom coordinates ',100x,'*')
 1100 format(1x,'* to be written to dumpfile section ',i4,63x,'*')
 2010 format(1x,'* calculate data # ',i3,' title: ',9a8,1x,'*')
 2020 format(1x,'* data type : ',a12,77x,'*')
 2030 format(1x,'* add data on section ',i4,' of dumpfile on ',a4,
     &     ' (block ',i3,') to current grid with scale factor',f10.3,
     &     1x,'*')
 2040 format(1x,'* apply scale factor : ',f10.5,71x,'*')
 3010 format(1x,'* generate plot # ',i3,' title: ',9a8,2x,'*')
 3020 format(1x,'* plot type : ',a17,72x,'*')
 4010 format(1x,'* restore ',a4,' from section ',i4,' of dumpfile on ',
     &     a4,' starting at block ',i3,29x,'*')
 9000 format(///)
      return
      end
**==mabat.f
      subroutine mabat(a,b,nr,nc,temp,nbmx)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension a(nbmx,*),b(nbmx,*),temp(nbmx,*)
c...
c...this subroutine takes the matrix product -
c...          c = a * b * at
c...where the matrix at is the transpose of a, a is nr by nc, b is
c...nc by nc, and c is nr by nr.
c...the matrix b is assumed symmetric thus only its upper triangular
c...half (including the diagonal) is used. the matrix c is also assumed
c...to be symmetric, and therefore only its upper triangular half is
c...calculated. this upper triangular half of c is returned in the upper
c...triangular half of b.
c...compute a * b . store that in y temporarily.
c...
      do 50 i = 1 , nr
         do 40 k = 1 , nc
            sum = 0.0d0
            do 20 j = 1 , k
               sum = sum + a(i,j)*b(j,k)
 20         continue
            k1 = k + 1
            if (k1.le.nc) then
               do 30 j = k1 , nc
                  sum = sum + a(i,j)*b(k,j)
 30            continue
            end if
            temp(i,k) = sum
 40      continue
 50   continue
c...multiply y * a to get upper half of c, storing that in b.
      do 80 k = 1 , nr
         do 70 i = 1 , k
            sum = 0.0d0
            do 60 j = 1 , nc
               sum = sum + temp(i,j)*a(k,j)
 60         continue
            b(i,k) = sum
 70      continue
 80   continue
c
      return
      end
**==molp.f
      subroutine molp(czeta, ctran, vect, scr, grid, data, imo)
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
      common/craypk/nbuf(9,maxorb),nbasin,norb,lenbas,ngrps,non,newbas,
     &     ncol,jcol,ngrp(3),iblkg,iblka
      dimension czeta(*), ctran(*), vect(*), scr(*), grid(*), data(*)
      dimension pop(2)
c
c mo amplitude on a grid
c
      i0 = 1
      i1 = i0 + ncol
      call getnpt(np,npx,npy,npz)
      do 10 i = 1,np
         call getpt(i,ax,ay,az,grid,1,idat)
c dummy pop
         call dplota(vect, ctran, czeta, pop, scr(i0),scr(i1),
     &        ax, ay, az, data(idat),imo)
 10   continue
      return
      end
**==mullin.f
      subroutine mullin(a,tran,b,scr,m,n,ndim)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension a(*),b(*),tran(ndim,*),scr(*)
      ij = 0
      do 40 j = 1 , m
         ik = 1
         scr(1) = a(1)*tran(1,j)
         do 20 i = 2 , n
            call daxpy(i-1,tran(i,j),a(ik+1),1,scr,1)
            scr(i) = ddot(i,a(ik+1),1,tran(1,j),1)
            ik = ik + i
 20      continue
         do 30 i = 1 , j
            ij = ij + 1
            b(ij) = ddot(n,tran(1,i),1,scr,1)
 30      continue
 40   continue
      return
      end
**==mullr.f
      subroutine mullr(vect,s,pmat,qmat,rmat,qqq,osmall)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      logical osmall
      character * 24 gao,glabel,gtit
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
c....    junk => junk2 here and in mulwrt (aolim uses junk) 
      common/junk2/omulla,omullo,omullg,omullp,omullu,
     +             mullac(maxorb),mulla,
     +eig(maxorb),frocc(maxorb),qdiag(maxorb),pdiag(maxorb),
     +rdiag(maxorb),omt(maxorb),ilabel(maxorb),igrp(maxorb)
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
      logical oming, oextg, odirec, ovcfre, osemi, ostopm, olevd
      logical osym_op,osym_dens
      integer mina, minb, mouta, moutb, lock, ibrk
      integer numg, isecg, iblk3g, mextra, nsymo, iorbsy, isunor
      real*8  gapa1, gapa2, gapb1, gapb2, scaleg, symtag
      common/atmol3/mina,minb,mouta,moutb,lock,ibrk,
     + numg(2),isecg(2),iblk3g(2),mextra,nsymo,iorbsy(100),
     + oming, oextg,
     + gapa1,gapa2,gapb1,gapb2,scaleg,symtag(9),odirec(50),
     + isunor,ovcfre,osemi,ostopm,olevd,osym_op,osym_dens
c
c
      real*8 ctran
      integer ilifc, ntran, itran
      logical otran,otri
      common/tran/ilifc(maxorb),ntran(maxorb),itran(mxorb3),
     +            ctran(mxorb3),otran,otri
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      integer mp, lwordp, lvec, ied3, newbas, lenb
      integer nbfnd, lenbad, lentrd, len3, ngrpm
      integer idena, nwroot, nosec
      logical noscf
      common /limy/ mp(10),lwordp,lvec,ied3,newbas,lenb,
     +              nbfnd,lenbad,lentrd,len3,ngrpm,
     +              idena,noscf,nwroot,nosec(mxroot)
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
      integer ltask, lcoord, iplot, iprop, imull, idma, inbo
      logical opass1, opass2, opass3, orege, opass4, offsym
      logical opass5, opass6, opass8, opass9, odisc, oirc
      logical orgall, orestj, opass10, opass11
      logical omrdci, ofulci, omech, omopac
      logical orpa, omclr, odirpa, occsd, occsdt, oqcisd, oqcisdt
      logical omrdcid,opark,oclunp,ofill,odiesel
c
      common /restrj/ ltask,lcoord,opass1,opass2,opass3,orege,iplot,
     + iprop,imull,opass4,offsym,opass5,opass6,opass10,idma,inbo,
     + omrdci,ofulci,omech,omopac,opass8,opass9,odisc,
     + oirc,orpa,omclr,odirpa,occsd,occsdt,oqcisd,oqcisdt,orgall,
     + omrdcid,orestj(7),opark,oclunp,ofill(6),odiesel, opass11
c
      common/bufb/gao(maxorb),glabel(maxorb),gtit(maxorb),
     * zpr(maxorb),zcomm(19),zzbuff(10),
     * zgroupp(maxorb),zmul(maxorb)
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
      dimension vect(*),s(*),pmat(*),qmat(*),rmat(*),qqq(*)
      data itvec/3/
      data zprin,zprino/'print','noprint'/
c
      nav = lenwrd()
      nwl = (2*maxorb+1)/nav
      nw1e2 = 400 + 204/nav
      nwp2 = 2*maxorb + 82 + 60/nav
      nwdma = 5 + 4*maxat + (13+maxat+nav)/nav
      iblk = 1
      if (idma.gt.0) iblk = iblk + lensec(maxat) + lensec(nwdma)
      if (iplot.gt.0) iblk = iblk + iplot*(lensec(nwp2)+lensec(140))
      if (iprop.gt.0) iblk = iblk + lensec(nw1e2) + lensec(831)
      if (local.gt.0) iblk = iblk + lensec(nwl)
      nw = (maxorb+6)
      call readi(omulla,nw,iblk,num8)
      call rdchrs(zmul,num,num8)
      iposv = mouta
      if (moutb.ne.0) iposv = moutb
      if (iposv.le.0 .or. iposv.gt.350)
     +    call caserr('invalid section specified for eigenvectors')
      call secget(iposv,itvec,iblvec)
      call getqp(zcomm,zzbuff,eig,frocc,nbas,newbas,ncol,ieig,idiff,
     +           maxorb,iblvec)
      write (iwr,6010) iposv , ibl3d , yed(ied3) , (zcomm(7-i),i=1,6) ,
     +                 zzbuff
c     nbsq = ncol*nbas (putq always write square)
      nbsq = nbas*nbas
      call rdedx(vect,nbsq,iblvec,idaf)
      call tdown(vect,ilifq,vect,ilifq,ncol)
      call mulltr(vect,s,pmat,qmat,qqq)
c
      if (omullp) then
         if (.not.otran) then
            write (iwr,6020)
            do 30 i = 1 , newbas
               n = ntran(i)
               l = ilifc(i) + 1
               write (iwr,6030) ctran(l) , itran(l) , n , i
               if (n.ne.1) then
                  do 20 j = 2 , n
                     l = l + 1
                     write (iwr,6030) ctran(l) , itran(l)
 20               continue
               end if
 30         continue
         end if
         write (iwr,6040) zcomm(5)
         call prsql(vect,ncol,nbas,nbas)
         if (idiff.gt.0) then
            write (iwr,6050) iposv
            write (iwr,6060) (frocc(i),i=1,ncol)
         end if
      end if
c
      if (omullu) then
         l1 = num
         l2 = l1*(l1+1)/2
         write(iwr,6300)
c....    check frocc vs. ne
         if (idiff.ne.1) call caserr('no occupations in mullr')
         sum = 0.0d0
         do i=1,l1
            sum = sum + frocc(i)
         end do
         if (dabs(sum-ne).gt.1.0d-10) then
            write(iwr,6307) sum, dabs(sum-ne)
6307        format(/' ** total occupation on natorb ',f8.2,
     1              ' deviates by ',e10.3,' from ne ')
            call caserr('wrong natorbs in mullr')
         end if
         if (osmall) call caserr('not enough core for mull unpair')
c....    generate density matrix
         call dmtx(pmat,vect,frocc,iky,l1,l1,l1)
c....    determine ao-s / atom (uses junk)
         call aolim
c
c...     calculate unpaired density-matrix
c...     V.N.Staroverov,E.R.Davidson,J.Am.Chem.Soc.122,186-187 (2000)
c...     V.N.Staroverov,E.R.Davidson,IJQC,77,316-323(2000)
c
c...      square d and s
c
         call squr(pmat,qmat,l1)
         call squr(s,rmat,l1)
         call mxma(qmat,1,l1,rmat,1,l1,qqq,1,l1,l1,l1,l1)
         call mxma(qqq,1,l1,qmat,1,l1,rmat,1,l1,l1,l1,l1)
         call trsqsq(rmat,qmat,l1)
         do i=1,l2
            pmat(i) = 2.0d0*pmat(i) - qmat(i)
         end do
c
c     ----- do a mulliken population analysis ----
c           calculate overlap population
c
         call vmul(pmat,1,s,1,pmat,1,l2)
c
c     ----- calculate total gross population in ao*s ----
c
         call grossc(pmat,rmat,l1)
         write (iwr,6301)
         do i = 1 , l1
            write (iwr,6302) i , gao(i) , rmat(i)
         end do
c
c     ----- compress from orbitals to atoms -----
c
         call atpop(pmat,qmat,nat)
c
c     ----- calculate total gross population on atoms -----
c
         call grossc(qmat,rmat,nat)
         write (iwr,6305)
         sum = 0.0d0
         do i = 1 , nat
            if (rmat(i).ne.0.0d0)
     1      write (iwr,6308) i , zaname(i) , czan(i) , rmat(i)
            sum = sum + rmat(i)
         end do
         write (iwr,6306) sum
c
6300    format(//7x,'&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&',
     1          /7x,'          Unpaired Electron Density ',
     2          /7x,' V.N.Staroverov,E.R.Davidson,JACS.122,186 (2000)',
     3          /7x,'&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&',
     4          /)
 6301 format (/10x,'----- total gross population in aos ------'/)
 6302 format (10x,i5,2x,a12,2x,f12.5)
 6305 format (/10x,'----- total gross population on atoms ----'/)    
 6306 format (/10x,'----- total gross population',f9.5,' ----'/)    
 6308 format (10x,i5,2x,a8,2x,f6.1,f12.5)   
         if (.not.(omulla.or.omullo.or.omullg)) return
      end if
c
      write (iwr,6090)
      call setstl(maxorb,.true.,omt)
      omullq = mulla.ne.0
      if (omullq) then
         do 40 i = 1 , mulla
            omt(mullac(i)) = .false.
 40      continue
      end if
      do 50 i = 1 , ncol
         zpr(i) = zprino
         if (.not.omt(i)) then
            zpr(i) = zprin
         end if
 50   continue
      i = 1
      j = 8
 60   if (ncol.ge.i) then
         if (j.gt.ncol) j = ncol
         write (iwr,6070) (k,k=i,j)
         write (iwr,6080) (zpr(k),k=i,j)
         i = i + 8
         j = j + 8
         go to 60
      else if (omulla) then
         call assia(glabel,iwr)
         go to 80
      end if
 70   if (omullo) then
         call assib(glabel,iwr)
         omullo = .false.
         go to 80
      end if
 160  if (.not.omullg) then
         return
      else
         do 170 i = 1,nbas
         glabel(i) = ' '
         glabel(i)(1:8) = zmul(i)
 170     continue   
         omullg = .false.
         write (iwr,6160)
      endif
 80   ngrpm = 1
      gtit(1) = glabel(1)
      ilabel(1) = 1
      do 90 i = 2 , nbas
         j = locatc(gtit,ngrpm,glabel(i))
         ilabel(i) = j
         if (j.eq.0) then
            ngrpm = ngrpm + 1
            gtit(ngrpm) = glabel(i)
            ilabel(i) = ngrpm
         end if
 90   continue
      write (iwr,6100)
      do 120 k = 1 , ngrpm
         ij = 0
         do 100 l = 1 , newbas
            if (ilabel(l).eq.k) then
               ij = ij + 1
               igrp(ij) = l
            end if
 100     continue
         i = 1
         j = 18
 110     if (ij.lt.i) then
            write (iwr,6120)
         else
            if (j.gt.ij) j = ij
            write (iwr,6110) gtit(k) , (igrp(m),m=i,j)
            i = i + 18
            j = j + 18
            go to 110
         end if
 120  continue
      lengrp = ikyp(ngrpm)
      call vclr(pdiag,1,newbas)
      call vclr(pmat,1,lengrp)
      if(omullq) then
        call vclr(rdiag,1,newbas)
        call vclr(rmat,1,lengrp)
      endif
      do 150 i = 1 , ncol
         if (.not.omt(i)) write (iwr,6130) i , (zzbuff(j),j=1,8)
         call vclr(qmat,1,lengrp)
         m = ilifq(i)
         do 140 j = 1 , newbas
            r = vect(m+j)
            qqq(j) = r
            ij = ilabel(j)
            ikyij = iky(ij)
            ikyj = iky(j)
            ijdiag = ij + ikyij
            dum = s(j+ikyj)*r*r
            qdiag(j) = dum
            qmat(ijdiag) = qmat(ijdiag) + dum
            l = j - 1
            if (l.gt.0) then
               do 130 k = 1 , l
                  dum = s(k+ikyj)*qqq(k)*r
                  qdiag(k) = qdiag(k) + dum
                  qdiag(j) = qdiag(j) + dum
                  ik = ilabel(k)
                  if (ij.lt.ik) then
                     ik = iky(ik) + ij
                     qmat(ik) = qmat(ik) + dum
                  else if (ij.eq.ik) then
                     qmat(ijdiag) = qmat(ijdiag) + dum + dum
                  else
                     ik = ik + ikyij
                     qmat(ik) = qmat(ik) + dum
                  end if
 130           continue
            end if
 140     continue
         dum = frocc(i)
         if (.not.omt(i)) then
            call mulwrt(qmat,qdiag,iwr)
            if(omullq) then
              call daxpy(lengrp,dum,qmat,1,rmat,1)
              call daxpy(newbas,dum,qdiag,1,rdiag,1)
            endif
         end if
         call daxpy(lengrp,dum,qmat,1,pmat,1)
         call daxpy(newbas,dum,qdiag,1,pdiag,1)
 150  continue
      if(omullq) then
        write (iwr,6170) (zzbuff(i),i=1,9)
        call mulwrt(rmat,rdiag,iwr)
      endif
      write (iwr,6140) (zzbuff(i),i=1,9)
      call mulwrt(pmat,pdiag,iwr)
      write (iwr,6150)
      go to 70
 6010 format (/1x,'vectors restored from section',i4,
     +        ' of dumpfile starting at block',i6,' of ',
     +        a4//' header block information : '/1x,27('-')
     +        /' vectors created under acct. ',a8/1x,a7,
     +        'vectors created by ',a8,' gamess module  at ',a8,' on ',
     +        a8,' in the job ',a8/' with the title : ',10a8/)
 6020 format (/40x,12('*')/40x,'list of sabf'/40x,12('*')//11x,
     +        'symmetry adapted orbitals'//
     +        ' coefficient old orbital nterm new orbital'/1x,41('-'))
 6030 format (1x,f11.7,i12,i6,i12)
 6040 format (//1x,104('-')//40x,16('*')/40x,a8,' vectors'/40x,16('*')
     +        //)
 6050 format (//40x,'orbital populations  (from vectors section ',i3,
     +        ')'//)
 6060 format (/15x,7f15.7)
 6070 format (//3x,8i14)
 6080 format (/8x,8(6x,a8))
 6090 format (/40x,23('-')/40x,'orbital analysis output'/40x,23('-'))
 6100 format (/15x,'basis functions are assigned to groups as follows :'
     +        //15x,'group',15x,'functions'//)
 6110 format (15x,a12,5x,18i5)
 6120 format (/)
 6130 format (/1x,104('*')//25x,'analysis of molecular orbital',i4,
     +        ' of ',8a8)
 6140 format (/1x,104('*')//25x,'total analysis of (all MOs)',9a8)
 6150 format (//1x,104('*'))
 6160 format (/40x,20('*')/40x,'input group analysis'/40x,20('*')/)
 6170 format (/1x,104('*')//25x,'total analysis (over printed MOs) of ',
     +         9a8)
      end
**==mullsq.f
      subroutine mullsq(a,tran,scr,m,n,ndim)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension a(ndim,*),tran(ndim,*),scr(*)
      do 30 j = 1 , m
         do 20 i = 1 , n
            scr(i) = ddot(n,tran(1,i),1,a(1,j),1)
 20      continue
         call dcopy(n,scr,1,a(1,j),1)
 30   continue
      return
      end
**==mulltr.f
      subroutine mulltr(q,s,scra,temp,qqq)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *24 zlabel
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
      common/bufb/zlabel(maxorb)
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
      dimension dtr(3,3),dtrinv(3,3),zdlab(3)
      dimension ftr(9,9),ftrinv(9,9)
      dimension gtr(15,15),gtrinv(15,15)
      character *24 zflab(9), zglab(15)
      dimension q(*),s(*),scra(*),temp(*),lench(3),qqq(*)
      dimension lenchf(9), lenchg(15)
c     dtrinv chosen such that ddot(dtr(1,i),dtrinv(1,i)).eq.1 and
c     i.ne.j ddot(dtrinv(1,i),dtrinv(1,j)).eq.0.
      data dtr/
     + 0.4472135954999579d0, 0.4472135954999579d0, 0.4472135954999579d0,
     + 0.8660254037844386d0,-0.8660254037844386d0, 0.0000000000000000d0,
     +-0.5000000000000000d0,-0.5000000000000000d0, 1.0000000000000000d0/
      data dtrinv/
     + 0.7453559924999299d0, 0.7453559924999299d0, 0.7453559924999299d0,
     + 0.5773502691896258d0,-0.5773502691896258d0, 0.0000000000000000d0,
     +-0.3333333333333333d0,-0.3333333333333333d0, 0.6666666666666667d0/
      data zdlab/'rr','xx-yy','zz-rr'/
      data lench/2,5,5/
c
      data ftr/
     + 0.452369896045359d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.000000000000000d0, 0.415088297188186d0,
     + 0.000000000000000d0, 0.415088297188186d0, 0.000000000000000d0,
c
     + 0.000000000000000d0, 0.452369896045359d0, 0.000000000000000d0,
     + 0.415088297188186d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.000000000000000d0, 0.415088297188186d0,
c
     + 0.000000000000000d0, 0.000000000000000d0, 0.452369896045359d0,
     + 0.000000000000000d0, 0.415088297188186d0, 0.000000000000000d0,
     + 0.415088297188186d0, 0.000000000000000d0, 0.000000000000000d0,
c
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.866025403784439d0, 0.000000000000000d0,
     +-0.866025403784439d0, 0.000000000000000d0, 0.000000000000000d0,
c
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.866025403784439d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.000000000000000d0,-0.866025403784439d0,
c
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.000000000000000d0, 0.866025403784439d0,
     + 0.000000000000000d0,-0.866025403784439d0, 0.000000000000000d0,
c
     + 1.106315011975948d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.000000000000000d0,-0.602839938334960d0,
     + 0.000000000000000d0,-0.602839938334960d0, 0.000000000000000d0,
c
     + 0.000000000000000d0, 1.106315011975948d0, 0.000000000000000d0,
     +-0.602839938334960d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.000000000000000d0,-0.602839938334960d0,
c
     + 0.000000000000000d0, 0.000000000000000d0, 1.106315011975948d0,
     + 0.000000000000000d0,-0.602839938334960d0, 0.000000000000000d0,
     +-0.602839938334960d0, 0.000000000000000d0, 0.000000000000000d0/
c
      data ftrinv/
     + 0.823636155716326d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.000000000000000d0, 0.755757030623969d0,
     + 0.000000000000000d0, 0.755757030623969d0, 0.000000000000000d0,
c
     + 0.000000000000000d0, 0.823636155716326d0, 0.000000000000000d0,
     + 0.755757030623969d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.000000000000000d0, 0.755757030623969d0,
c
     + 0.000000000000000d0, 0.000000000000000d0, 0.823636155716326d0,
     + 0.000000000000000d0, 0.755757030623969d0, 0.000000000000000d0,
     + 0.755757030623969d0, 0.000000000000000d0, 0.000000000000000d0,
c
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.577350269189626d0, 0.000000000000000d0,
     +-0.577350269189626d0, 0.000000000000000d0, 0.000000000000000d0,
c
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.577350269189626d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.000000000000000d0,-0.577350269189626d0,
c
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.000000000000000d0, 0.577350269189626d0,
     + 0.000000000000000d0,-0.577350269189626d0, 0.000000000000000d0,
c
     + 0.567118579308447d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.000000000000000d0,-0.309027470185270d0,
     + 0.000000000000000d0,-0.309027470185270d0, 0.000000000000000d0,
c
     + 0.000000000000000d0, 0.567118579308447d0, 0.000000000000000d0,
     +-0.309027470185270d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.000000000000000d0,-0.309027470185270d0,
c
     + 0.000000000000000d0, 0.000000000000000d0, 0.567118579308447d0,
     + 0.000000000000000d0,-0.309027470185270d0, 0.000000000000000d0,
     +-0.309027470185270d0, 0.000000000000000d0, 0.000000000000000d0/
c
      data zflab/"(xx+yy+zz)x","(xx+yy+zz)y","(xx+yy+zz)z",
     +           "(xx-yy)z","(xx-zz)y","(yy-zz)x",
     +           "(xx-yy-zz)x","(yy-xx-zz)y","(zz-xx-yy)z"/
      data lenchf/11,11,11,8,8,8,11,11,11/
c
      data gtr/
c     s contaminant (xx+yy+zz)(xx+yy+zz)
     + 0.226400900893753d0, 0.226400900893753d0, 0.226400900893753d0,
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.284569229358613d0, 0.284569229358613d0, 0.284569229358613d0,
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
c     d contaminants (xx+yy+zz)xy
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.422577127364258d0, 0.000000000000000d0, 0.422577127364258d0,
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.000000000000000d0, 0.377964473009227d0,
c     d contaminants (xx+yy+zz)xz
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.422577127364258d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.422577127364258d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.377964473009227d0, 0.000000000000000d0,
c     d contaminants (xx+yy+zz)yz
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.422577127364258d0, 0.000000000000000d0, 0.422577127364258d0,
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.377964473009227d0, 0.000000000000000d0, 0.000000000000000d0,
c     d contaminants (xx+yy+zz)(2zz-xx-yy)
     +-0.265652889512170d0,-0.265652889512170d0, 0.531305779024340d0,
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     +-0.531305779024340d0, 0.265652889512170d0, 0.265652889512170d0,
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
c     d contaminants (xx+yy+zz)(xx-yy)
     + 0.460124301812560d0,-0.460124301812560d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.460124301812560d0,-0.460124301812560d0,
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
c     g function (2zz-xx-yy)(2zz-xx-yy)
     +-0.059946136315707d0,-0.059946136315707d0, 1.082029681713666d0, 
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.684507202397193d0,-0.724987842512736d0,-0.724987842512736d0,
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
c     g function (2zz-xx-yy)(xx-yy)
     +-0.719341569456438d0, 0.719341569456438d0, 0.000000000000000d0, 
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.887854507721722d0,-0.887854507721722d0,
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
c     g function (xx-yy)(xx-yy)
     + 0.901153502522778d0, 0.901153502522778d0, 0.402908470755982d0, 
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     +-0.994792935325224d0,-0.379829034742677d0,-0.379829034742677d0,
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
c     g function (xxxy-yyyx)
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     + 1.118033988749660d0, 0.000000000000000d0,-1.118033988749660d0,
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
c     g function (xxxz-zzzx)
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.000000000000000d0, 1.118033988749660d0, 0.000000000000000d0,
     + 0.000000000000000d0,-1.118033988749660d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
c     g function (yyyz-zzzy)
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.000000000000000d0, 00000000000000000d0, 0.000000000000000d0,
     + 1.118033988749660d0, 00000000000000000d0,-1.118033988749660d0,
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
c     g function (xxxy-zzxy+yyyx)
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.487950036474267d0, 0.000000000000000d0, 0.487950036474267d0,
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.000000000000000d0,-1.091089451179962d0,
c     g function (xxxz-yyxz+zzzx)
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.487950036474267d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.487950036474267d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.000000000000000d0,-1.091089451179962d0, 0.000000000000000d0,
c     g function (yyyz-xxyz+zzzy)
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     + 0.487950036474267d0, 0.000000000000000d0, 0.487950036474267d0, 
     + 0.000000000000000d0, 0.000000000000000d0, 0.000000000000000d0,
     +-1.091089451179962d0, 0.000000000000000d0, 0.000000000000000d0/
c
      data gtrinv/
c 1
     +           0.570694728763102d0,
     +           0.570694728763102d0,
     +           0.570694728763102d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.717321170198664d0,
     +           0.717321170198664d0,
     +           0.717321170198664d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
c 2
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.845154254728527d0,
     +           0.000000000000000d0,
     +           0.845154254728527d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.755928946018463d0,
c 3
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.845154254728527d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.845154254728527d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.755928946018463d0,
     +           0.000000000000000d0,
c 4
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.845154254728527d0,
     +           0.000000000000000d0,
     +           0.845154254728527d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.755928946018463d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
c 5
     +          -0.313692553792150d0,
     +          -0.313692553792150d0,
     +           0.627385107584300d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +          -0.627385107584300d0,
     +           0.313692553792150d0,
     +           0.313692553792150d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
c 6
     +           0.543331441124037d0,
     +          -0.543331441124037d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.543331441124037d0,
     +          -0.543331441124037d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
c 7
     +          -0.022220871049042d0,
     +          -0.022220871049042d0,
     +           0.401087434592432d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.253733555011834d0,
     +          -0.268738943837110d0,
     +          -0.268738943837110d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
c 8
     +          -0.275454039513396d0,
     +           0.275454039513396d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.339981895995427d0,
     +          -0.339981895995427d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
c 9
     +           0.294048377621490d0,
     +           0.294048377621490d0,
     +           0.131469923630195d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +          -0.324603131300939d0,
     +          -0.123939052699624d0,
     +          -0.123939052699623d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
c 10
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.447213595500052d0,
     +           0.000000000000000d0,
     +          -0.447213595500052d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
c 11
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.447213595500052d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +          -0.447213595500052d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
c 12
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.447213595500052d0,
     +           0.000000000000000d0,
     +          -0.447213595500052d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
c 13
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.292770021884567d0,
     +           0.000000000000000d0,
     +           0.292770021884567d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +          -0.654653670707980d0,
c 14
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.292770021884567d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.292770021884567d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +          -0.654653670707980d0,
     +           0.000000000000000d0,
c 15
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.292770021884567d0,
     +           0.000000000000000d0,
     +           0.292770021884567d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0,
     +          -0.654653670707980d0,
     +           0.000000000000000d0,
     +           0.000000000000000d0/
c
      data zglab/"rrrr","rrxy","rrxz",
     +           "rryz","rr(2zz-xx-yy)",
     +           "rr(xx-yy)","zzzz+xxyy-xxzz-yyzz",
     +           "xxxx-yyyy+xxzz-yyzz","rrrr-xxyy-xxzz-yyzz",
     +           "xxxy-yyyx","xxxz-zzzx","yyyz-zzzy",
     +           "xxxy-zzxy+yyyx","xxxz-yyxz+zzzx","yyyz-xxyz+zzzy"/
      data lenchg/4,4,4,4,13,7,19,19,19,9,9,9,14,14,14/
c
c     Fudge the g-function inverse transformation to conserve the
c     number of electrons, for some reason the conservation of the 
c     electrons depends on the choice of transformation, simply
c     ensuring the transformation vectors are orthonormal is somehow
c     not enough. Why?
c
      fudge = dsqrt(2.0d0/2.06706349206384)
      call dscal(15,fudge,gtrinv(1,5),1)
      fudge = dsqrt(2.0d0/2.06706349206385)
      call dscal(15,fudge,gtrinv(1,6),1)
      fudge = dsqrt(2.0d0/2.01629460008033)
      call dscal(15,fudge,gtrinv(1,7),1)
      fudge = dsqrt(2.0d0/2.00705236542181)
      call dscal(15,fudge,gtrinv(1,8),1)
      fudge = dsqrt(2.0d0/2.00914346107270)
      call dscal(15,fudge,gtrinv(1,9),1)
c
c     Enough fudging
c
      l2 = num*(num+1)/2
      l3 = num*num
c
c     re-set zbflab
c
      call setlab
c
      do 20 i = 1 , num
         zlabel(i) = zbflab(i)
 20   continue
      odbas = .false.
      ofbas = .false.
      ogbas = .false.
      call vclr(scra,1,l3)
      idiag = 1
      do 30 i = 1 , num
         scra(idiag) = 1.0d0
         idiag = idiag + num + 1
 30   continue
      do 70 i = 1 , nshell
         if (ktype(i).eq.3) then
            odbas = .true.
            ifbf = kloc(i)
            icol = ilifq(ifbf)
            do 50 j = 1 , 3
               irow = ifbf
               do 40 k = 1 , 3
                  scra(icol+irow) = dtr(k,j)
                  irow = irow + 1
 40            continue
               icol = icol + num
 50         continue
            ibfj = ifbf
            do 60 j = 1 , 3
               if (nat.ge.100) then
                 len = lench(j) + 6
                 zlabel(ibfj)(7:len) = zdlab(j)(1:lench(j))
               else
                 len = lench(j) + 5
                 zlabel(ibfj)(6:len) = zdlab(j)(1:lench(j))
               endif
               ibfj = ibfj + 1
 60         continue
         else if (ktype(i).eq.4) then
            ofbas = .true.
            ifbf = kloc(i)
            icol = ilifq(ifbf)
            do j = 1 , 9
               irow = ifbf
               do k = 1 , 9
                  scra(icol+irow) = ftr(k,j)
                  irow = irow + 1
               enddo ! k
               icol = icol + num
            enddo ! j
            ibfj = ifbf
            do j = 1 , 9
               if (nat.ge.100) then
                 len = lenchf(j) + 6
                 zlabel(ibfj)(7:len) = zflab(j)(1:lenchf(j))
               else
                 len = lenchf(j) + 5
                 zlabel(ibfj)(6:len) = zflab(j)(1:lenchf(j))
               endif
               ibfj = ibfj + 1
            enddo ! j
         else if (ktype(i).eq.5) then
            ogbas = .true.
            ifbf = kloc(i)
            icol = ilifq(ifbf)
            do j = 1 , 15
               irow = ifbf
               do k = 1 , 15
                  scra(icol+irow) = gtr(k,j)
                  irow = irow + 1
               enddo ! k
               icol = icol + num
            enddo ! j
            ibfj = ifbf
            do j = 1 , 15
               if (nat.ge.100) then
                 len = lenchg(j) + 6
                 zlabel(ibfj)(7:len) = zglab(j)(1:lenchg(j))
               else
                 len = lenchg(j) + 5
                 zlabel(ibfj)(6:len) = zglab(j)(1:lenchg(j))
               endif
               ibfj = ibfj + 1
            enddo ! j
         end if
 70   continue
      if (odbas.or.ofbas.or.ogbas) then
c
c        Note: the locodi if statement does not work if f- or 
c        g-functions are to be treated.
c
c        locodi = iky(ifbf+1) + ifbf
c        if (s(locodi).gt.0.001d0) then
            call mullin(s,scra,temp,qqq,num,num,num)
            call dcopy(l2,temp,1,s,1)
            if (nprint.eq.5) call prtril(s,num)
            do 100 i = 1 , nshell
               if (ktype(i).eq.3) then
                  ifbf = kloc(i)
                  icol = ilifq(ifbf)
                  do 90 j = 1 , 3
                     irow = ifbf
                     do 80 k = 1 , 3
                        scra(icol+irow) = dtrinv(k,j)
                        irow = irow + 1
 80                  continue
                     icol = icol + num
 90               continue
               else if (ktype(i).eq.4) then
                  ifbf = kloc(i)
                  icol = ilifq(ifbf)
                  do j = 1 , 9
                     irow = ifbf
                     do k = 1 , 9
                        scra(icol+irow) = ftrinv(k,j)
                        irow = irow + 1
                     enddo ! k
                     icol = icol + num
                  enddo ! j
               else if (ktype(i).eq.5) then
                  ifbf = kloc(i)
                  icol = ilifq(ifbf)
                  do j = 1 , 15
                     irow = ifbf
                     do k = 1 , 15
                        scra(icol+irow) = gtrinv(k,j)
                        irow = irow + 1
                     enddo ! k
                     icol = icol + num
                  enddo ! j
               end if
 100        continue
            call mullsq(q,scra,qqq,num,num,num)
c        end if
      end if
      return
      end
**==mulwrt.f
      subroutine mulwrt(am,ad,iwr)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character * 24 zlabel,glabel,gtit
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
      parameter ( mxorb5 = maxorb*5)
       dimension am(*),ad(*)
c
      integer mp, lwordp, lvec, ied3, newbas, lenb
      integer nbfnd, lenbad, lentrd, len3, ngrpm
      integer idena, nwroot, nosec
      logical noscf
      common /limy/ mp(10),lwordp,lvec,ied3,newbas,lenb,
     +              nbfnd,lenbad,lentrd,len3,ngrpm,
     +              idena,noscf,nwroot,nosec(mxroot)
c
      common/bufb/zlabel(maxorb),glabel(maxorb),gtit(maxorb)
      common/junk2/omull(5),mullac(maxorb),mulla,eig(mxorb5),
     * omt(maxorb),ilab(maxorb)
      ind(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)
      write (iwr,6010)
      m = 1
      n = 4
 20   if (newbas.lt.m) then
         write (iwr,6030)
         i = 1
         j = 4
 30      if (ngrpm.lt.i) then
            write (iwr,6040)
            call vclr(eig,1,ngrpm)
            do 40 i = 1 , newbas
               j = ilab(i)
               eig(j) = eig(j) + ad(i)
 40         continue
            i = 1
            j = 4
 50         if (ngrpm.lt.i) then
               return
            else
               if (j.gt.ngrpm) j = ngrpm
               write (iwr,6050) (gtit(k),k=i,j)
               write (iwr,6060) (eig(k),k=i,j)
               write (iwr,6070)
               i = i + 4
               j = j + 4
               go to 50
            end if
         else
            if (j.gt.ngrpm) j = ngrpm
            write (iwr,6050) (gtit(k),k=i,j)
            write (iwr,6070)
            do 60 k = 1 , ngrpm
               write (iwr,6080) gtit(k) , (am(ind(k,l)),l=i,j)
 60         continue
            i = i + 4
            j = j + 4
            go to 30
         end if
      else
         if (n.gt.newbas) n = newbas
         write (iwr,6090) (zlabel(k),k=m,n)
         write (iwr,6020) (ad(k),k=m,n)
         m = m + 4
         n = n + 4
         go to 20
      end if
 6010 format (///50x,'basis function populations'/50x,26('*'))
 6020 format (14x,4(f10.7,16x))
 6030 format (///50x,'group population matrix'/50x,23('*'))
 6040 format (///50x,'group populations'/50x,17('*'))
 6050 format (/14x,4(2x,a12))
 6060 format (/14x,4(f10.7,4x))
 6070 format (/)
 6080 format (2x,a12,4(f10.7,4x))
 6090 format (/14x,4(2x,a24))
      end
**==mwt.f
      subroutine mwt(w,nr,nc,nmx)
c
c...       transpose a matrix
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension w(nmx,nmx)
      do 30 j = 1 , nr
         do 20 i = j , nc
            tem = w(j,i)
            w(j,i) = w(i,j)
            w(i,j) = tem
 20      continue
 30   continue
c
      return
      end
**==nstuf.f
      function nstuf()
      common/stuf/nwrd
      nstuf = 1 + (nwrd - 1)/lenwrd()
      return
      end
**==nxtstp.f
      subroutine nxtstp(jtype,jarg)
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
c the job sequence
      integer stype(mxstp), arg(mxstp)
      common/route/stype,arg,nstep,istep
c
      if(istep.gt.nstep)then
         jtype = 0
      else
         jtype = stype(istep)
         jarg = arg(istep)
      endif
      istep = istep+1
      return
      end
**==octnuc.f
      subroutine octnuc(noc,vlist,c,dint,nocmx)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension dint(10),vlist(nocmx,4),c(3)
      common/miscop/cnuc(10)
c...
c...nuclear contributions to the octopole moments.
c...xxx,xyy,xzz,yxx,yyy,yzz,zxx,zyy,zzz,xyz.
c...
      call vclr(cnuc,1,10)
c...
      do 50 i = 1 , noc
         mu = 0
         do 40 j = 1 , 3
            do 30 k = 1 , j
               do 20 m = 1 , k
                  mu = mu + 1
                  cnuc(mu) = cnuc(mu) + vlist(i,4)*(vlist(i,j)-c(j))
     +                       *(vlist(i,k)-c(k))*(vlist(i,m)-c(m))
 20            continue
 30         continue
 40      continue
 50   continue
c...
      rrx = cnuc(1) + cnuc(3) + cnuc(8)
      rry = cnuc(2) + cnuc(4) + cnuc(9)
      rrz = cnuc(5) + cnuc(7) + cnuc(10)
      dint(1) = (5.0d0*cnuc(1)-3.0d0*rrx)/2.0d0
      dint(2) = (5.0d0*cnuc(3)-rrx)/2.0d0
      dint(3) = (5.0d0*cnuc(8)-rrx)/2.0d0
      dint(4) = (5.0d0*cnuc(2)-rry)/2.0d0
      dint(5) = (5.0d0*cnuc(4)-3.0d0*rry)/2.0d0
      dint(6) = (5.0d0*cnuc(9)-rry)/2.0d0
      dint(7) = (5.0d0*cnuc(5)-rrz)/2.0d0
      dint(8) = (5.0d0*cnuc(7)-rrz)/2.0d0
      dint(9) = (5.0d0*cnuc(10)-3.0d0*rrz)/2.0d0
      dint(10) = 5.0d0*cnuc(6)/2.0d0
c...
      return
      end
**==onuc.f
      subroutine  onuc(dint)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension dint(*)
      common/prpspc/cc(3),mop,icent,nt
c
      do 20 i = 1 , nt
         dint(i) = 0.0d0
 20   continue
c
      return
      end
**==opaa0.f
      subroutine opaa0(l,m,n,ga,v,nt,d)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension v(10),d(3),fnu(5),fn(5),fd(5),dp(3)
      data pi/3.14159265d0/
      data fn/1.0d0,1.0d0,2.0d0,6.0d0,24.0d0/
      data fd/1.0d0,1.0d0,.5d0,.166666667d0,0.041666667d0/
c
      nt = 1
      itot = l + m + n
      itoth = itot/2
      pre = (2.0d0*pi/ga)*fn(l+1)*fn(m+1)*fn(n+1)
      if (2*itoth.ne.itot) then
         pre = -pre
      end if
      del = .25d0/ga
      p = ga*(d(1)**2+d(2)**2+d(3)**2)
      pp = 2.0d0*p
      call fmc(4,p,expmx,fmch)
      fnu(5) = fmch
      fnu(4) = (expmx+pp*fnu(5))/7.0d0
      fnu(3) = (expmx+pp*fnu(4))/5.0d0
      fnu(2) = (expmx+pp*fnu(3))/3.0d0
      fnu(1) = expmx + pp*fnu(2)
      dp(1) = 1.0d0
      dp(2) = del
      dp(3) = del**2
      v(1) = pre*aainer(0,0,0,l,m,n,d,dp,fnu,fn,fd)
      do 20 i = 1 , nt
         v(i) = -v(i)
 20   continue
c
      return
      end
**==opaa1.f
      subroutine opaa1(l,m,n,ga,v,nt,d)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension v(10),d(3),fd(7),fn(7)
      common/miscop/fnu(7),dp(4)
      data pi/3.14159265d0/
      data fn/1.0d0,1.0d0,2.0d0,6.0d0,24.0d0,120.0d0,720.d0/
      data fd/1.0d0,1.0d0,.5d0,
     1 0.166666667d0,.416666667d-1,.833333333d-2,.138888889d-2/
c
      nt = 3
      itot = l + m + n
      itoth = itot/2
      pre = 4.0d0*pi*fn(l+1)*fn(m+1)*fn(n+1)
      if (2*itoth.ne.itot) then
         pre = -pre
      end if
      del = .25d0/ga
      p = ga*(d(1)**2+d(2)**2+d(3)**2)
      pp = 2.0d0*p
      call fmc(6,p,expmx,fmch)
      fnu(7) = fmch
      fnu(6) = (expmx+pp*fnu(7))/11.0d0
      fnu(5) = (expmx+pp*fnu(6))/9.0d0
      fnu(4) = (expmx+pp*fnu(5))/7.0d0
      fnu(3) = (expmx+pp*fnu(4))/5.0d0
      fnu(2) = (expmx+pp*fnu(3))/3.0d0
      fnu(1) = expmx + pp*fnu(2)
      dp(1) = 1.0d0
      dp(2) = del
      dp(3) = del**2
      dp(4) = del*dp(3)
      v(1) = pre*aainer(1,0,0,l,m,n,d,dp,fnu,fn,fd)
      v(2) = pre*aainer(0,1,0,l,m,n,d,dp,fnu,fn,fd)
      v(3) = pre*aainer(0,0,1,l,m,n,d,dp,fnu,fn,fd)
c
      return
      end
**==opaa2.f
      subroutine opaa2(l,m,n,ga,v,nt,d)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension v(10),d(3),fd(7),fn(7)
      common/miscop/fnu(7),dp(4)
      data pi/3.14159265d0/,fn/1.0d0,1.0d0,2.0d0,6.0d0,24.0d0,120.0d0,72
     10.0d0/,fd/1.0d0,1.0d0,.5d0,
     1.166666667d0,.041666667d0,.008333333d0,.001388889d0/
c
      itot = l + m + n
      itoth = itot/2
      pre = (8.0d0*pi*ga)*fn(l+1)*fn(m+1)*fn(n+1)
      if (2*itoth.ne.itot) then
         pre = -pre
      end if
      del = .25d0/ga
      p = ga*(d(1)**2+d(2)**2+d(3)**2)
      pp = 2.0d0*p
      call fmc(6,p,expmx,fmch)
      fnu(7) = fmch
      fnu(6) = (expmx+pp*fnu(7))/11.0d0
      fnu(5) = (expmx+pp*fnu(6))/9.0d0
      fnu(4) = (expmx+pp*fnu(5))/7.0d0
      fnu(3) = (expmx+pp*fnu(4))/5.0d0
      fnu(2) = (expmx+pp*fnu(3))/3.0d0
      fnu(1) = expmx + pp*fnu(2)
      dp(1) = 1.0d0
      dp(2) = del
      dp(3) = del**2
      dp(4) = del*dp(3)
      v(2) = pre*aainer(0,2,0,l,m,n,d,dp,fnu,fn,fd)
      v(1) = pre*aainer(2,0,0,l,m,n,d,dp,fnu,fn,fd)
      v(3) = pre*aainer(0,0,2,l,m,n,d,dp,fnu,fn,fd)
      v(4) = pre*aainer(1,1,0,l,m,n,d,dp,fnu,fn,fd)
      v(5) = pre*aainer(1,0,1,l,m,n,d,dp,fnu,fn,fd)
      v(6) = pre*aainer(0,1,1,l,m,n,d,dp,fnu,fn,fd)
      call opac3(l,m,n,ga,dp(1),nt,d)
      dp(1) = 1.33333333d0*pi*dp(1)
      v(1) = v(1) + dp(1)
      v(2) = v(2) + dp(1)
      v(3) = v(3) + dp(1)
      nt = 6
c
      return
      end
**==opaa3.f
      subroutine opaa3(l,m,n,ga,v,nt,d)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension vx(3),d(3),v(10)
c
      nt = 10
      call opaa1(l,m,n,ga,vx,nx,d)
      v(1) = vx(1)*d(1)
      v(2) = vx(2)*d(2)
      v(3) = vx(3)*d(3)
      v(4) = vx(2)*d(1)
      v(5) = vx(3)*d(1)
      v(6) = vx(3)*d(2)
      call opaa1(l+1,m,n,ga,vx,nx,d)
      v(1) = v(1) + vx(1)
      v(4) = v(4) + vx(2)
      v(5) = v(5) + vx(3)
      call opaa1(l,m+1,n,ga,vx,nx,d)
      v(2) = v(2) + vx(2)
      v(6) = v(6) + vx(3)
      call opaa1(l,m,n+1,ga,vx,nx,d)
      v(3) = v(3) + vx(3)
      do 20 i = 1 , 6
         v(i) = -v(i)
 20   continue
      v(7) = v(1) + v(2) + v(3)
      v(8) = v(1)
      v(9) = v(2)
      v(10) = v(3)
      v(1) = 1.5d0*(v(9)+v(10))
      v(2) = 1.5d0*(v(8)+v(10))
      v(3) = 1.5d0*(v(8)+v(9))
      v(4) = -1.5d0*v(4)
      v(5) = -1.5d0*v(5)
      v(6) = -1.5d0*v(6)
c
      return
      end
**==opab1.f
      subroutine opab1(l,m,n,ga,v,nt,d)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension v(10),d(3)
      common/miscop/g(4)
c...
c...compute integrals for dipole moment.
c...x,y,z.
c...
      nt = 3
      g(1) = slap(l,m,n,ga)
      g(2) = slap(l+1,m,n,ga)
      g(3) = slap(l,m+1,n,ga)
      g(4) = slap(l,m,n+1,ga)
      v(1) = g(2) + d(1)*g(1)
      v(2) = g(3) + d(2)*g(1)
      v(3) = g(4) + d(3)*g(1)
c...
      do 20 i = 1 , nt
         v(i) = -v(i)
 20   continue
c...
      return
      end
**==opab10.f
      subroutine opab10(l,m,n,ga,v,nt,dd)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension v(10),dd(3)
      common/miscop/g(35),d(35)
c...
c...evaluate integrals for odd fourth moments.
c...xxxy,xyyy,xyzz,xzxx,xzyy,xzzz,yzxx,yzyy,yzzz.
c...
      nt = 9
      g(1) = slap(l,m,n,ga)
      g(2) = slap(l+1,m,n,ga)
      g(3) = slap(l,m+1,n,ga)
      g(4) = slap(l,m,n+1,ga)
      g(5) = slap(l+2,m,n,ga)
      g(6) = slap(l,m+2,n,ga)
      g(7) = slap(l,m,n+2,ga)
      g(8) = slap(l+1,m+1,n,ga)
      g(9) = slap(l+1,m,n+1,ga)
      g(10) = slap(l,m+1,n+1,ga)
      g(11) = slap(l+3,m,n,ga)
      g(12) = slap(l,m+3,n,ga)
      g(13) = slap(l,m,n+3,ga)
      g(14) = slap(l+2,m+1,n,ga)
      g(15) = slap(l+2,m,n+1,ga)
      g(16) = slap(l+1,m+2,n,ga)
      g(17) = slap(l,m+2,n+1,ga)
      g(18) = slap(l+1,m,n+2,ga)
      g(19) = slap(l,m+1,n+2,ga)
      g(20) = slap(l+1,m+1,n+1,ga)
      g(21) = slap(l+4,m,n,ga)
      g(22) = slap(l,m+4,n,ga)
      g(23) = slap(l,m,n+4,ga)
      g(24) = slap(l+3,m+1,n,ga)
      g(25) = slap(l+3,m,n+1,ga)
      g(26) = slap(l+1,m+3,n,ga)
      g(27) = slap(l,m+3,n+1,ga)
      g(28) = slap(l+1,m,n+3,ga)
      g(29) = slap(l,m+1,n+3,ga)
      g(30) = slap(l+2,m+2,n,ga)
      g(31) = slap(l+2,m,n+2,ga)
      g(32) = slap(l,m+2,n+2,ga)
      g(33) = slap(l+2,m+1,n+1,ga)
      g(34) = slap(l+1,m+2,n+1,ga)
      g(35) = slap(l+1,m+1,n+2,ga)
      d(2) = dd(1)
      d(3) = dd(2)
      d(4) = dd(3)
      d(5) = dd(1)**2
      d(6) = dd(2)**2
      d(7) = dd(3)**2
      d(8) = d(2)*d(3)
      d(9) = d(2)*d(4)
      d(10) = d(3)*d(4)
      d(11) = d(5)*d(2)
      d(12) = d(6)*d(3)
      d(13) = d(7)*d(4)
      d(14) = d(5)*d(3)
      d(15) = d(5)*d(4)
      d(16) = d(6)*d(2)
      d(17) = d(6)*d(4)
      d(18) = d(7)*d(2)
      d(19) = d(7)*d(3)
      d(20) = d(8)*d(4)
      d(21) = d(5)*d(5)
      d(22) = d(6)*d(6)
      d(23) = d(7)*d(7)
      d(24) = d(5)*d(8)
      d(25) = d(5)*d(9)
      d(26) = d(6)*d(8)
      d(27) = d(6)*d(10)
      d(28) = d(7)*d(9)
      d(29) = d(7)*d(10)
      d(30) = d(5)*d(6)
      d(31) = d(5)*d(7)
      d(32) = d(6)*d(7)
      d(33) = d(5)*d(10)
      d(34) = d(6)*d(9)
      d(35) = d(7)*d(8)
c...
      v(1) = g(24) + 3.0d0*d(2)*g(14) + 3.0d0*d(5)*g(8) + d(11)*g(3)
     +       + d(3)*g(11) + 3.0d0*d(8)*g(5) + 3.0d0*d(14)*g(2) + d(24)
     +       *g(1)
      v(2) = g(26) + 3.0d0*d(3)*g(16) + 3.0d0*d(6)*g(8) + d(12)*g(2)
     +       + d(2)*g(12) + 3.0d0*d(8)*g(6) + 3.0d0*d(16)*g(3) + d(26)
     +       *g(1)
      v(3) = g(35) + d(35)*g(1) + 2.0d0*(d(4)*g(20)+d(20)*g(4)) + d(7)
     +       *g(8) + d(8)*g(7) + 2.0d0*(d(10)*g(9)+d(9)*g(10)) + d(3)
     +       *g(18) + d(18)*g(3) + d(2)*g(19) + d(19)*g(2)
      v(4) = g(25) + d(25)*g(1) + d(4)*g(11) + d(11)*g(4)
     +       + 3.0d0*(d(2)*g(15)+d(15)*g(2)+d(9)*g(5)+d(5)*g(9))
      v(5) = g(34) + d(34)*g(1) + d(2)*g(17) + d(17)*g(2) + d(4)*g(16)
     +       + d(16)*g(4) + d(9)*g(6) + d(6)*g(9)
     +       + 2.0d0*(d(3)*g(20)+d(20)*g(3)+d(8)*g(10)+d(10)*g(8))
      v(6) = g(28) + d(28)*g(1) + d(13)*g(2) + d(2)*g(13)
     +       + 3.0d0*(d(4)*g(18)+d(18)*g(4)+d(7)*g(9)+d(9)*g(7))
      v(7) = g(33) + d(33)*g(1) + d(3)*g(15) + d(15)*g(3) + d(4)*g(14)
     +       + d(14)*g(4) + d(10)*g(5) + d(5)*g(10)
     +       + 2.0d0*(d(2)*g(20)+d(20)*g(2)+d(8)*g(9)+d(9)*g(8))
      v(8) = g(27) + d(27)*g(1) + d(12)*g(4) + d(4)*g(12)
     +       + 3.0d0*(d(3)*g(17)+d(17)*g(3)+d(10)*g(6)+d(6)*g(10))
      v(9) = g(29) + d(29)*g(1) + d(3)*g(13) + d(13)*g(3)
     +       + 3.0d0*(d(4)*g(19)+d(19)*g(4)+d(7)*g(10)+d(10)*g(7))
c...
      do 20 i = 1 , nt
         v(i) = -v(i)
 20   continue
c...
      return
      end
**==opab2.f
      subroutine opab2(l,m,n,ga,v,nt,d)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension v(10),d(3)
      common/miscop/g(10)
c...
c...evaluate integrals for quadrupole moments.
c...qm(xx)=(3*xx-rr)/2, similarly for qm(yy) and qm(zz).
c...qm(xy)=1.5*xy, similarly for qm(xz) and qm(yz).
c...also  rr=xx+yy+zz.
c...
      nt = 7
      g(1) = slap(l,m,n,ga)
      g(2) = slap(l+1,m,n,ga)
      g(3) = slap(l,m+1,n,ga)
      g(4) = slap(l,m,n+1,ga)
      g(5) = slap(l+2,m,n,ga)
      g(6) = slap(l,m+2,n,ga)
      g(7) = slap(l,m,n+2,ga)
      g(8) = slap(l+1,m+1,n,ga)
      g(9) = slap(l+1,m,n+1,ga)
      g(10) = slap(l,m+1,n+1,ga)
      v(1) = g(5) + 2.0d0*d(1)*g(2) + d(1)**2*g(1)
      v(2) = g(6) + 2.0d0*d(2)*g(3) + d(2)**2*g(1)
      v(3) = g(7) + 2.0d0*d(3)*g(4) + d(3)**2*g(1)
      v(4) = (g(8)+d(1)*g(3)+d(2)*g(2)+d(1)*d(2)*g(1))/2.0d0
      v(5) = (g(9)+d(1)*g(4)+d(3)*g(2)+d(1)*d(3)*g(1))/2.0d0
      v(6) = (g(10)+d(2)*g(4)+d(3)*g(3)+d(2)*d(3)*g(1))/2.0d0
c...
      v(7) = v(1) + v(2) + v(3)
      v(1) = (3.0d0*v(1)-v(7))/2.0d0
      v(2) = (3.0d0*v(2)-v(7))/2.0d0
      v(3) = (3.0d0*v(3)-v(7))/2.0d0
      v(4) = 3.0d0*v(4)
      v(5) = 3.0d0*v(5)
      v(6) = 3.0d0*v(6)
c...
      do 20 i = 1 , nt
         v(i) = -v(i)
 20   continue
c...
      return
      end
**==opab3.f
      subroutine opab3(l,m,n,ga,v,nt,dd)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension v(10),dd(3)
      common/miscop/g(20),d(20)
c...
c...evaluate octopole moment integrals.
c...xxx,xyy,xzz,yxx,yyy,yzz,zxx,zyy,zzz,xyz.
c...
      nt = 10
      g(1) = slap(l,m,n,ga)
      g(2) = slap(l+1,m,n,ga)
      g(3) = slap(l,m+1,n,ga)
      g(4) = slap(l,m,n+1,ga)
      g(5) = slap(l+2,m,n,ga)
      g(6) = slap(l,m+2,n,ga)
      g(7) = slap(l,m,n+2,ga)
      g(8) = slap(l+1,m+1,n,ga)
      g(9) = slap(l+1,m,n+1,ga)
      g(10) = slap(l,m+1,n+1,ga)
      g(11) = slap(l+3,m,n,ga)
      g(12) = slap(l,m+3,n,ga)
      g(13) = slap(l,m,n+3,ga)
      g(14) = slap(l+2,m+1,n,ga)
      g(15) = slap(l+2,m,n+1,ga)
      g(16) = slap(l+1,m+2,n,ga)
      g(17) = slap(l,m+2,n+1,ga)
      g(18) = slap(l+1,m,n+2,ga)
      g(19) = slap(l,m+1,n+2,ga)
      g(20) = slap(l+1,m+1,n+1,ga)
      d(2) = dd(1)
      d(3) = dd(2)
      d(4) = dd(3)
      d(5) = dd(1)**2
      d(6) = dd(2)**2
      d(7) = dd(3)**2
      d(8) = d(2)*d(3)
      d(9) = d(2)*d(4)
      d(10) = d(3)*d(4)
      d(11) = d(5)*d(2)
      d(12) = d(6)*d(3)
      d(13) = d(7)*d(4)
      d(14) = d(5)*d(3)
      d(15) = d(5)*d(4)
      d(16) = d(6)*d(2)
      d(17) = d(6)*d(4)
      d(18) = d(7)*d(2)
      d(19) = d(7)*d(3)
      d(20) = d(8)*d(4)
c...
      v(1) = g(11) + 3.0d0*d(2)*g(5) + 3.0d0*d(5)*g(2) + d(11)*g(1)
      v(2) = g(16) + 2.0d0*d(3)*g(8) + d(2)*g(6) + d(6)*g(2)
     +       + 2.0d0*d(8)*g(3) + d(16)*g(1)
      v(3) = g(18) + 2.0d0*d(4)*g(9) + d(2)*g(7) + d(7)*g(2)
     +       + 2.0d0*d(9)*g(4) + d(18)*g(1)
      v(4) = g(14) + 2.0d0*d(2)*g(8) + d(3)*g(5) + d(5)*g(3)
     +       + 2.0d0*d(8)*g(2) + d(14)*g(1)
      v(5) = g(12) + 3.0d0*d(3)*g(6) + 3.0d0*d(6)*g(3) + d(12)*g(1)
      v(6) = g(19) + 2.0d0*d(4)*g(10) + d(3)*g(7) + d(7)*g(3)
     +       + 2.0d0*d(10)*g(4) + d(19)*g(1)
      v(7) = g(15) + 2.0d0*d(2)*g(9) + d(4)*g(5) + d(5)*g(4)
     +       + 2.0d0*d(9)*g(2) + d(15)*g(1)
      v(8) = g(17) + 2.0d0*d(3)*g(10) + d(4)*g(6) + d(6)*g(4)
     +       + 2.0d0*d(10)*g(3) + d(17)*g(1)
      v(9) = g(13) + 3.0d0*d(4)*g(7) + 3.0d0*d(7)*g(4) + d(13)*g(1)
      v(10) = g(20) + d(2)*g(10) + d(3)*g(9) + d(4)*g(8) + d(10)*g(2)
     +        + d(9)*g(3) + d(8)*g(4) + d(20)*g(1)
c...
      do 20 i = 1 , nt
         v(i) = -v(i)
 20   continue
c...
      rrx = v(1) + v(2) + v(3)
      rry = v(4) + v(5) + v(6)
      rrz = v(7) + v(8) + v(9)
      v(1) = (5.0d0*v(1)-3.0d0*rrx)/2.0d0
      v(2) = (5.0d0*v(2)-rrx)/2.0d0
      v(3) = (5.0d0*v(3)-rrx)/2.0d0
      v(4) = (5.0d0*v(4)-rry)/2.0d0
      v(5) = (5.0d0*v(5)-3.0d0*rry)/2.0d0
      v(6) = (5.0d0*v(6)-rry)/2.0d0
      v(7) = (5.0d0*v(7)-rrz)/2.0d0
      v(8) = (5.0d0*v(8)-rrz)/2.0d0
      v(9) = (5.0d0*v(9)-3.0d0*rrz)/2.0d0
      v(10) = 5.0d0*v(10)/2.0d0
c...
      return
      end
**==opab4.f
      subroutine opab4(l,m,n,ga,v,nt,d)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension v(10),d(3)
      common/miscop/g(10)
c...
c...compute integrals for diamagnetic suceptibilities.
c...
      nt = 7
      g(1) = slap(l,m,n,ga)
      g(2) = slap(l+1,m,n,ga)
      g(3) = slap(l,m+1,n,ga)
      g(4) = slap(l,m,n+1,ga)
      g(5) = slap(l+2,m,n,ga)
      g(6) = slap(l,m+2,n,ga)
      g(7) = slap(l,m,n+2,ga)
      g(8) = slap(l+1,m+1,n,ga)
      g(9) = slap(l+1,m,n+1,ga)
      g(10) = slap(l,m+1,n+1,ga)
      v(1) = g(5) + 2.0d0*d(1)*g(2) + d(1)**2*g(1)
      v(2) = g(6) + 2.0d0*d(2)*g(3) + d(2)**2*g(1)
      v(3) = g(7) + 2.0d0*d(3)*g(4) + d(3)**2*g(1)
      v(4) = g(8) + d(1)*g(3) + d(2)*g(2) + d(1)*d(2)*g(1)
      v(5) = g(9) + d(1)*g(4) + d(3)*g(2) + d(1)*d(3)*g(1)
      v(6) = g(10) + d(2)*g(4) + d(3)*g(3) + d(2)*d(3)*g(1)
c...
      do 20 i = 1 , nt
         v(i) = -v(i)
 20   continue
c...
      v(7) = v(1) + v(2) + v(3)
      v(1) = 1.5d0*(v(7)-v(1))
      v(2) = 1.5d0*(v(7)-v(2))
      v(3) = 1.5d0*(v(7)-v(3))
      v(4) = -1.5d0*v(4)
      v(5) = -1.5d0*v(5)
      v(6) = -1.5d0*v(6)
c...
      return
      end
**==opab5.f
      subroutine opab5(l,m,n,ga,v,nt,dd)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension v(10),dd(3)
      common/miscop/g(35),d(35)
c...
c...evaluate integrals for hexadecapole moments.
c...p4x=(35xxxx-30xxrr+3rrrr)/8,p4y=(35yyyy-30yyrr+3rrrr)/8,
c...p4z=(35zzzz-30zzrr+3rrrr)/8,xxrr,yyrr,zzrr,rrrr,xyrr,xzrr,yzrr.
c...
      nt = 10
      g(1) = slap(l,m,n,ga)
      g(2) = slap(l+1,m,n,ga)
      g(3) = slap(l,m+1,n,ga)
      g(4) = slap(l,m,n+1,ga)
      g(5) = slap(l+2,m,n,ga)
      g(6) = slap(l,m+2,n,ga)
      g(7) = slap(l,m,n+2,ga)
      g(8) = slap(l+1,m+1,n,ga)
      g(9) = slap(l+1,m,n+1,ga)
      g(10) = slap(l,m+1,n+1,ga)
      g(11) = slap(l+3,m,n,ga)
      g(12) = slap(l,m+3,n,ga)
      g(13) = slap(l,m,n+3,ga)
      g(14) = slap(l+2,m+1,n,ga)
      g(15) = slap(l+2,m,n+1,ga)
      g(16) = slap(l+1,m+2,n,ga)
      g(17) = slap(l,m+2,n+1,ga)
      g(18) = slap(l+1,m,n+2,ga)
      g(19) = slap(l,m+1,n+2,ga)
      g(20) = slap(l+1,m+1,n+1,ga)
      g(21) = slap(l+4,m,n,ga)
      g(22) = slap(l,m+4,n,ga)
      g(23) = slap(l,m,n+4,ga)
      g(24) = slap(l+3,m+1,n,ga)
      g(25) = slap(l+3,m,n+1,ga)
      g(26) = slap(l+1,m+3,n,ga)
      g(27) = slap(l,m+3,n+1,ga)
      g(28) = slap(l+1,m,n+3,ga)
      g(29) = slap(l,m+1,n+3,ga)
      g(30) = slap(l+2,m+2,n,ga)
      g(31) = slap(l+2,m,n+2,ga)
      g(32) = slap(l,m+2,n+2,ga)
      g(33) = slap(l+2,m+1,n+1,ga)
      g(34) = slap(l+1,m+2,n+1,ga)
      g(35) = slap(l+1,m+1,n+2,ga)
      d(2) = dd(1)
      d(3) = dd(2)
      d(4) = dd(3)
      d(5) = dd(1)**2
      d(6) = dd(2)**2
      d(7) = dd(3)**2
      d(8) = d(2)*d(3)
      d(9) = d(2)*d(4)
      d(10) = d(3)*d(4)
      d(11) = d(5)*d(2)
      d(12) = d(6)*d(3)
      d(13) = d(7)*d(4)
      d(14) = d(5)*d(3)
      d(15) = d(5)*d(4)
      d(16) = d(6)*d(2)
      d(17) = d(6)*d(4)
      d(18) = d(7)*d(2)
      d(19) = d(7)*d(3)
      d(20) = d(8)*d(4)
      d(21) = d(5)*d(5)
      d(22) = d(6)*d(6)
      d(23) = d(7)*d(7)
      d(24) = d(5)*d(8)
      d(25) = d(5)*d(9)
      d(26) = d(6)*d(8)
      d(27) = d(6)*d(10)
      d(28) = d(7)*d(9)
      d(29) = d(7)*d(10)
      d(30) = d(5)*d(6)
      d(31) = d(5)*d(7)
      d(32) = d(6)*d(7)
      d(33) = d(5)*d(10)
      d(34) = d(6)*d(9)
      d(35) = d(7)*d(8)
c...
      v(1) = g(21) + 4.0d0*d(2)*g(11) + 6.0d0*d(5)*g(5) + 4.0d0*d(11)
     +       *g(2) + d(21)*g(1)
      v(2) = g(22) + 4.0d0*d(3)*g(12) + 6.0d0*d(6)*g(6) + 4.0d0*d(12)
     +       *g(3) + d(22)*g(1)
      v(3) = g(23) + 4.0d0*d(4)*g(13) + 6.0d0*d(7)*g(7) + 4.0d0*d(13)
     +       *g(4) + d(23)*g(1)
      v(4) = g(30) + 2.0d0*d(3)*g(14) + d(6)*g(5) + 2.0d0*d(2)*g(16)
     +       + 4.0d0*d(8)*g(8) + 2.0d0*d(16)*g(2) + d(5)*g(6)
     +       + 2.0d0*d(14)*g(3) + d(30)*g(1)
      v(5) = g(31) + 2.0d0*d(4)*g(15) + d(7)*g(5) + 2.0d0*d(2)*g(18)
     +       + 4.0d0*d(9)*g(9) + 2.0d0*d(18)*g(2) + d(5)*g(7)
     +       + 2.0d0*d(15)*g(4) + d(31)*g(1)
      v(6) = g(32) + 2.0d0*d(4)*g(17) + d(7)*g(6) + 2.0d0*d(3)*g(19)
     +       + 4.0d0*d(10)*g(10) + 2.0d0*d(19)*g(3) + d(6)*g(7)
     +       + 2.0d0*d(17)*g(4) + d(32)*g(1)
      pprr = v(1) + v(4) + v(5)
      qqrr = v(4) + v(2) + v(6)
      ssrr = v(5) + v(6) + v(3)
      rrrr = pprr + qqrr + ssrr
      pppp = (35.0d0*v(1)-30.0d0*pprr+3.0d0*rrrr)/8.0d0
      qqqq = (35.0d0*v(2)-30.0d0*qqrr+3.0d0*rrrr)/8.0d0
      ssss = (35.0d0*v(3)-30.0d0*ssrr+3.0d0*rrrr)/8.0d0
      v(1) = g(24) + 3.0d0*d(2)*g(14) + 3.0d0*d(5)*g(8) + d(11)*g(3)
     +       + d(3)*g(11) + 3.0d0*d(8)*g(5) + 3.0d0*d(14)*g(2) + d(24)
     +       *g(1)
      v(2) = g(26) + 3.0d0*d(3)*g(16) + 3.0d0*d(6)*g(8) + d(12)*g(2)
     +       + d(2)*g(12) + 3.0d0*d(8)*g(6) + 3.0d0*d(16)*g(3) + d(26)
     +       *g(1)
      v(3) = g(35) + d(35)*g(1) + 2.0d0*(d(4)*g(20)+d(20)*g(4)) + d(7)
     +       *g(8) + d(8)*g(7) + 2.0d0*(d(10)*g(9)+d(9)*g(10)) + d(3)
     +       *g(18) + d(18)*g(3) + d(2)*g(19) + d(19)*g(2)
      v(4) = g(25) + d(25)*g(1) + d(4)*g(11) + d(11)*g(4)
     +       + 3.0d0*(d(2)*g(15)+d(15)*g(2)+d(9)*g(5)+d(5)*g(9))
      v(5) = g(34) + d(34)*g(1) + d(2)*g(17) + d(17)*g(2) + d(4)*g(16)
     +       + d(16)*g(4) + d(9)*g(6) + d(6)*g(9)
     +       + 2.0d0*(d(3)*g(20)+d(20)*g(3)+d(8)*g(10)+d(10)*g(8))
      v(6) = g(28) + d(28)*g(1) + d(13)*g(2) + d(2)*g(13)
     +       + 3.0d0*(d(4)*g(18)+d(18)*g(4)+d(7)*g(9)+d(9)*g(7))
      v(7) = g(33) + d(33)*g(1) + d(3)*g(15) + d(15)*g(3) + d(4)*g(14)
     +       + d(14)*g(4) + d(10)*g(5) + d(5)*g(10)
     +       + 2.0d0*(d(2)*g(20)+d(20)*g(2)+d(8)*g(9)+d(9)*g(8))
      v(8) = g(27) + d(27)*g(1) + d(12)*g(4) + d(4)*g(12)
     +       + 3.0d0*(d(3)*g(17)+d(17)*g(3)+d(10)*g(6)+d(6)*g(10))
      v(9) = g(29) + d(29)*g(1) + d(3)*g(13) + d(13)*g(3)
     +       + 3.0d0*(d(4)*g(19)+d(19)*g(4)+d(7)*g(10)+d(10)*g(7))
      pqrr = v(1) + v(2) + v(3)
      psrr = v(4) + v(5) + v(6)
      qsrr = v(7) + v(8) + v(9)
      v(1) = pppp
      v(2) = qqqq
      v(3) = ssss
      v(4) = pprr
      v(5) = qqrr
      v(6) = ssrr
      v(7) = rrrr
      v(8) = pqrr
      v(9) = psrr
      v(10) = qsrr
c...
      do 20 i = 1 , nt
         v(i) = -v(i)
 20   continue
c...
      return
      end
**==opab6.f
      subroutine opab6(l,m,n,ga,v,nt,d)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension v(10),d(3)
      common/miscop/g(10)
c...
c...compute integrals for second moments.
c...xx,yy,zz,xy,xz,yz,rr.
c...
      nt = 7
      g(1) = slap(l,m,n,ga)
      g(2) = slap(l+1,m,n,ga)
      g(3) = slap(l,m+1,n,ga)
      g(4) = slap(l,m,n+1,ga)
      g(5) = slap(l+2,m,n,ga)
      g(6) = slap(l,m+2,n,ga)
      g(7) = slap(l,m,n+2,ga)
      g(8) = slap(l+1,m+1,n,ga)
      g(9) = slap(l+1,m,n+1,ga)
      g(10) = slap(l,m+1,n+1,ga)
      v(1) = g(5) + 2.0d0*d(1)*g(2) + d(1)**2*g(1)
      v(2) = g(6) + 2.0d0*d(2)*g(3) + d(2)**2*g(1)
      v(3) = g(7) + 2.0d0*d(3)*g(4) + d(3)**2*g(1)
      v(4) = g(8) + d(1)*g(3) + d(2)*g(2) + d(1)*d(2)*g(1)
      v(5) = g(9) + d(1)*g(4) + d(3)*g(2) + d(1)*d(3)*g(1)
      v(6) = g(10) + d(2)*g(4) + d(3)*g(3) + d(2)*d(3)*g(1)
      v(7) = v(1) + v(2) + v(3)
c...
      do 20 i = 1 , nt
         v(i) = -v(i)
 20   continue
c...
      return
      end
**==opab7.f
      subroutine opab7(l,m,n,ga,v,nt,dd)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension v(10),dd(3)
      common/miscop/g(20),d(20)
c...
c...evaluate third moment integrals.
c...xxx,xyy,xzz,yxx,yyy,yzz,zxx,zyy,zzz,xyz.
c...
      nt = 10
      g(1) = slap(l,m,n,ga)
      g(2) = slap(l+1,m,n,ga)
      g(3) = slap(l,m+1,n,ga)
      g(4) = slap(l,m,n+1,ga)
      g(5) = slap(l+2,m,n,ga)
      g(6) = slap(l,m+2,n,ga)
      g(7) = slap(l,m,n+2,ga)
      g(8) = slap(l+1,m+1,n,ga)
      g(9) = slap(l+1,m,n+1,ga)
      g(10) = slap(l,m+1,n+1,ga)
      g(11) = slap(l+3,m,n,ga)
      g(12) = slap(l,m+3,n,ga)
      g(13) = slap(l,m,n+3,ga)
      g(14) = slap(l+2,m+1,n,ga)
      g(15) = slap(l+2,m,n+1,ga)
      g(16) = slap(l+1,m+2,n,ga)
      g(17) = slap(l,m+2,n+1,ga)
      g(18) = slap(l+1,m,n+2,ga)
      g(19) = slap(l,m+1,n+2,ga)
      g(20) = slap(l+1,m+1,n+1,ga)
      d(2) = dd(1)
      d(3) = dd(2)
      d(4) = dd(3)
      d(5) = dd(1)**2
      d(6) = dd(2)**2
      d(7) = dd(3)**2
      d(8) = d(2)*d(3)
      d(9) = d(2)*d(4)
      d(10) = d(3)*d(4)
      d(11) = d(5)*d(2)
      d(12) = d(6)*d(3)
      d(13) = d(7)*d(4)
      d(14) = d(5)*d(3)
      d(15) = d(5)*d(4)
      d(16) = d(6)*d(2)
      d(17) = d(6)*d(4)
      d(18) = d(7)*d(2)
      d(19) = d(7)*d(3)
      d(20) = d(8)*d(4)
c...
      v(1) = g(11) + 3.0d0*d(2)*g(5) + 3.0d0*d(5)*g(2) + d(11)*g(1)
      v(2) = g(16) + 2.0d0*d(3)*g(8) + d(2)*g(6) + d(6)*g(2)
     +       + 2.0d0*d(8)*g(3) + d(16)*g(1)
      v(3) = g(18) + 2.0d0*d(4)*g(9) + d(2)*g(7) + d(7)*g(2)
     +       + 2.0d0*d(9)*g(4) + d(18)*g(1)
      v(4) = g(14) + 2.0d0*d(2)*g(8) + d(3)*g(5) + d(5)*g(3)
     +       + 2.0d0*d(8)*g(2) + d(14)*g(1)
      v(5) = g(12) + 3.0d0*d(3)*g(6) + 3.0d0*d(6)*g(3) + d(12)*g(1)
      v(6) = g(19) + 2.0d0*d(4)*g(10) + d(3)*g(7) + d(7)*g(3)
     +       + 2.0d0*d(10)*g(4) + d(19)*g(1)
      v(7) = g(15) + 2.0d0*d(2)*g(9) + d(4)*g(5) + d(5)*g(4)
     +       + 2.0d0*d(9)*g(2) + d(15)*g(1)
      v(8) = g(17) + 2.0d0*d(3)*g(10) + d(4)*g(6) + d(6)*g(4)
     +       + 2.0d0*d(10)*g(3) + d(17)*g(1)
      v(9) = g(13) + 3.0d0*d(4)*g(7) + 3.0d0*d(7)*g(4) + d(13)*g(1)
      v(10) = g(20) + d(2)*g(10) + d(3)*g(9) + d(4)*g(8) + d(10)*g(2)
     +        + d(9)*g(3) + d(8)*g(4) + d(20)*g(1)
c...
      do 20 i = 1 , nt
         v(i) = -v(i)
 20   continue
c...
      return
      end
**==opab8.f
      subroutine opab8(l,m,n,ga,v,nt,dd)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension v(10),dd(3)
      common/miscop/g(20),d(20)
c...
c...evaluate certain sums of third moment integrals.
c...rrx,rry,rrz.
c...
      nt = 3
      g(1) = slap(l,m,n,ga)
      g(2) = slap(l+1,m,n,ga)
      g(3) = slap(l,m+1,n,ga)
      g(4) = slap(l,m,n+1,ga)
      g(5) = slap(l+2,m,n,ga)
      g(6) = slap(l,m+2,n,ga)
      g(7) = slap(l,m,n+2,ga)
      g(8) = slap(l+1,m+1,n,ga)
      g(9) = slap(l+1,m,n+1,ga)
      g(10) = slap(l,m+1,n+1,ga)
      g(11) = slap(l+3,m,n,ga)
      g(12) = slap(l,m+3,n,ga)
      g(13) = slap(l,m,n+3,ga)
      g(14) = slap(l+2,m+1,n,ga)
      g(15) = slap(l+2,m,n+1,ga)
      g(16) = slap(l+1,m+2,n,ga)
      g(17) = slap(l,m+2,n+1,ga)
      g(18) = slap(l+1,m,n+2,ga)
      g(19) = slap(l,m+1,n+2,ga)
      d(2) = dd(1)
      d(3) = dd(2)
      d(4) = dd(3)
      d(5) = dd(1)**2
      d(6) = dd(2)**2
      d(7) = dd(3)**2
      d(8) = d(2)*d(3)
      d(9) = d(2)*d(4)
      d(10) = d(3)*d(4)
      d(11) = d(5)*d(2)
      d(12) = d(6)*d(3)
      d(13) = d(7)*d(4)
      d(14) = d(5)*d(3)
      d(15) = d(5)*d(4)
      d(16) = d(6)*d(2)
      d(17) = d(6)*d(4)
      d(18) = d(7)*d(2)
      d(19) = d(7)*d(3)
c...
      v(1) = g(11) + 3.0d0*d(2)*g(5) + 3.0d0*d(5)*g(2) + d(11)*g(1)
      v(2) = g(16) + 2.0d0*d(3)*g(8) + d(2)*g(6) + d(6)*g(2)
     +       + 2.0d0*d(8)*g(3) + d(16)*g(1)
      v(3) = g(18) + 2.0d0*d(4)*g(9) + d(2)*g(7) + d(7)*g(2)
     +       + 2.0d0*d(9)*g(4) + d(18)*g(1)
      v(4) = g(14) + 2.0d0*d(2)*g(8) + d(3)*g(5) + d(5)*g(3)
     +       + 2.0d0*d(8)*g(2) + d(14)*g(1)
      v(5) = g(12) + 3.0d0*d(3)*g(6) + 3.0d0*d(6)*g(3) + d(12)*g(1)
      v(6) = g(19) + 2.0d0*d(4)*g(10) + d(3)*g(7) + d(7)*g(3)
     +       + 2.0d0*d(10)*g(4) + d(19)*g(1)
      v(7) = g(15) + 2.0d0*d(2)*g(9) + d(4)*g(5) + d(5)*g(4)
     +       + 2.0d0*d(9)*g(2) + d(15)*g(1)
      v(8) = g(17) + 2.0d0*d(3)*g(10) + d(4)*g(6) + d(6)*g(4)
     +       + 2.0d0*d(10)*g(3) + d(17)*g(1)
      v(9) = g(13) + 3.0d0*d(4)*g(7) + 3.0d0*d(7)*g(4) + d(13)*g(1)
c...
      do 20 i = 1 , 9
         v(i) = -v(i)
 20   continue
c...
      v(1) = v(1) + v(2) + v(3)
      v(2) = v(4) + v(5) + v(6)
      v(3) = v(7) + v(8) + v(9)
c...
      return
      end
**==opab9.f
      subroutine opab9(l,m,n,ga,v,nt,dd)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension v(10),dd(3)
      common/miscop/g(35),d(35)
c...
c...evaluate integrals for even fourth moments.
c...xxxx,yyyy,zzzz,xxyy,xxzz,yyzz,rrxx,rryy,rrzz,rrrr.
c...
      nt = 10
      g(1) = slap(l,m,n,ga)
      g(2) = slap(l+1,m,n,ga)
      g(3) = slap(l,m+1,n,ga)
      g(4) = slap(l,m,n+1,ga)
      g(5) = slap(l+2,m,n,ga)
      g(6) = slap(l,m+2,n,ga)
      g(7) = slap(l,m,n+2,ga)
      g(8) = slap(l+1,m+1,n,ga)
      g(9) = slap(l+1,m,n+1,ga)
      g(10) = slap(l,m+1,n+1,ga)
      g(11) = slap(l+3,m,n,ga)
      g(12) = slap(l,m+3,n,ga)
      g(13) = slap(l,m,n+3,ga)
      g(14) = slap(l+2,m+1,n,ga)
      g(15) = slap(l+2,m,n+1,ga)
      g(16) = slap(l+1,m+2,n,ga)
      g(17) = slap(l,m+2,n+1,ga)
      g(18) = slap(l+1,m,n+2,ga)
      g(19) = slap(l,m+1,n+2,ga)
      g(21) = slap(l+4,m,n,ga)
      g(22) = slap(l,m+4,n,ga)
      g(23) = slap(l,m,n+4,ga)
      g(30) = slap(l+2,m+2,n,ga)
      g(31) = slap(l+2,m,n+2,ga)
      g(32) = slap(l,m+2,n+2,ga)
      d(2) = dd(1)
      d(3) = dd(2)
      d(4) = dd(3)
      d(5) = dd(1)**2
      d(6) = dd(2)**2
      d(7) = dd(3)**2
      d(8) = d(2)*d(3)
      d(9) = d(2)*d(4)
      d(10) = d(3)*d(4)
      d(11) = d(5)*d(2)
      d(12) = d(6)*d(3)
      d(13) = d(7)*d(4)
      d(14) = d(5)*d(3)
      d(15) = d(5)*d(4)
      d(16) = d(6)*d(2)
      d(17) = d(6)*d(4)
      d(18) = d(7)*d(2)
      d(19) = d(7)*d(3)
      d(21) = d(5)*d(5)
      d(22) = d(6)*d(6)
      d(23) = d(7)*d(7)
      d(30) = d(5)*d(6)
      d(31) = d(5)*d(7)
      d(32) = d(6)*d(7)
c...
      v(1) = g(21) + 4.0d0*d(2)*g(11) + 6.0d0*d(5)*g(5) + 4.0d0*d(11)
     +       *g(2) + d(21)*g(1)
      v(2) = g(22) + 4.0d0*d(3)*g(12) + 6.0d0*d(6)*g(6) + 4.0d0*d(12)
     +       *g(3) + d(22)*g(1)
      v(3) = g(23) + 4.0d0*d(4)*g(13) + 6.0d0*d(7)*g(7) + 4.0d0*d(13)
     +       *g(4) + d(23)*g(1)
      v(4) = g(30) + 2.0d0*d(3)*g(14) + d(6)*g(5) + 2.0d0*d(2)*g(16)
     +       + 4.0d0*d(8)*g(8) + 2.0d0*d(16)*g(2) + d(5)*g(6)
     +       + 2.0d0*d(14)*g(3) + d(30)*g(1)
      v(5) = g(31) + 2.0d0*d(4)*g(15) + d(7)*g(5) + 2.0d0*d(2)*g(18)
     +       + 4.0d0*d(9)*g(9) + 2.0d0*d(18)*g(2) + d(5)*g(7)
     +       + 2.0d0*d(15)*g(4) + d(31)*g(1)
      v(6) = g(32) + 2.0d0*d(4)*g(17) + d(7)*g(6) + 2.0d0*d(3)*g(19)
     +       + 4.0d0*d(10)*g(10) + 2.0d0*d(19)*g(3) + d(6)*g(7)
     +       + 2.0d0*d(17)*g(4) + d(32)*g(1)
      v(7) = v(1) + v(4) + v(5)
      v(8) = v(4) + v(2) + v(6)
      v(9) = v(5) + v(6) + v(3)
      v(10) = v(7) + v(8) + v(9)
c...
      do 20 i = 1 , nt
         v(i) = -v(i)
 20   continue
c...
      return
      end
**==opac1.f
      subroutine opac1(l,m,n,ga,v,nt,d)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension v(10),d(3),dftr(6)
      common/miscop/dd(3)
      data pi/3.14159265d0/
      data dftr/1.0d0,1.0d0,3.0d0,15.0d0,105.0d0,945.0d0/
c
      nt = 3
      do 20 i = 1 , 3
         dd(i) = -d(i)
         v(i) = 0.0d0
 20   continue
      g = .5d0/ga
      lh = l/2
      mh = m/2
      nh = n/2
      pre = 2.0d0*g*pi
      if (2*nh.ne.n) then
c
         if (2*mh.ne.m) go to 30
         if (2*lh.ne.l) go to 30
      else if (2*mh.ne.m) then
c
         if (2*lh.ne.l) go to 30
         v(2) = 1.0d0
         if (m.ne.0) then
            v(2) = dd(2)**m
         end if
         v(2) = pre*v(2)*dexp(-ga*dd(2)**2)*g**(lh+nh)*dftr(lh+1)
     +          *dftr(nh+1)
         return
      else
         v(1) = 1.0d0
         if (l.ne.0) then
            v(1) = dd(1)**l
         end if
         v(1) = pre*v(1)*dexp(-ga*dd(1)**2)*g**(mh+nh)*dftr(mh+1)
     +          *dftr(nh+1)
         if (2*lh.ne.l) go to 30
         v(2) = 1.0d0
         if (m.ne.0) then
            v(2) = dd(2)**m
         end if
         v(2) = pre*v(2)*dexp(-ga*dd(2)**2)*g**(lh+nh)
     +          *dftr(lh+1)*dftr(nh+1)
      end if
      v(3) = 1.0d0
      if (n.ne.0) then
         v(3) = dd(3)**n
      end if
      v(3) = pre*v(3)*dexp(-ga*dd(3)**2)*g**(lh+mh)
     +       *dftr(lh+1)*dftr(mh+1)
      return
 30   return
c
      end
**==opac2.f
      subroutine opac2(l,m,n,ga,v,nt,d)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension v(10),d(3),dftr(5)
      common/miscop/dd(3),ds(3)
      data sqrpi/1.77245385d0/,dftr/1.0d0,1.0d0,3.0d0,15.0d0,105.0d0/
c
      nt = 3
      do 20 i = 1 , 3
         dd(i) = -d(i)
         ds(i) = d(i)**2
         v(i) = 0.0d0
 20   continue
      g = .5d0/ga
      pre = sqrpi*dsqrt(2.0d0*g)
      lh = l/2
      mh = m/2
      nh = n/2
      if (2*nh.eq.n) then
         v(1) = 1.0d0
         if (l.ne.0) then
            v(1) = dd(1)**l
         end if
         if (m.ne.0) then
            v(1) = v(1)*dd(2)**m
         end if
         v(1) = pre*g**nh*v(1)*dexp(-ga*(ds(1)+ds(2)))*dftr(nh+1)
      end if
      if (2*mh.eq.m) then
         v(2) = 1.0d0
         if (l.ne.0) then
            v(2) = dd(1)**l
         end if
         if (n.ne.0) then
            v(2) = v(2)*dd(3)**n
         end if
         v(2) = pre*g**mh*v(2)*dexp(-ga*(ds(1)+ds(3)))*dftr(mh+1)
      end if
      if (2*lh.eq.l) then
         v(3) = 1.0d0
         if (m.ne.0) then
            v(3) = dd(2)**m
         end if
         if (n.ne.0) then
            v(3) = v(3)*dd(3)**n
         end if
         v(3) = pre*g**lh*v(3)*dexp(-ga*(ds(2)+ds(3)))*dftr(lh+1)
      end if
c
      return
      end
**==opac3.f
      subroutine opac3(l,m,n,ga,v,nt,d)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension v(10),d(3)
      common/miscop/dd(3)
c
      nt = 1
      ddsq = 0.0d0
      do 20 i = 1 , 3
         dd(i) = -d(i)
         ddsq = ddsq + d(i)**2
 20   continue
      v(1) = 1.0d0
      if (l.ne.0) then
         v(1) = dd(1)**l
      end if
      if (m.ne.0) then
         v(1) = v(1)*dd(2)**m
      end if
      if (n.ne.0) then
         v(1) = v(1)*dd(3)**n
      end if
      v(1) = v(1)*dexp(-ga*ddsq)
c
      return
      end
**==opad1.f
      subroutine opad1(l,m,n,ga,v,nt)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension v(10)
c
      nt = 1
      v(1) = slap(l,m,n,ga)
c
      return
      end
**==output.f
      subroutine output(s,t,u,p,q,pcap,qcap,
     *fc,fo,dt,dos,iwr)
c
c    one configuration
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
      dimension s(*),t(*),u(*),p(*),q(*),pcap(*),qcap(*)
      dimension fc(*),fo(*),dt(*),dos(*)
      common/junk/pqn(99),expa(99),ctran(mxprms,30),ajmn(12),czn,
     *expmn,c(15,15),itran(mxprms,30),nrow(30),nbcon(4),nbas(4),
     *nsym,nosh(4),ncsh(4),nccup(4),ibase(33),
     *popul(4),facto(30),fatt(30),thre(15),a(15,15),eivr(15,15),
     *cm2scf(15,15),cm1scf(15,15),eps(15),gmin(15,15),smin(15,15),
     *populs(4),term,energ,pot,cin,vir,thrscf,enorm,
     *red,expdif,encor,
     *nsize(10),nsavf(4),n1(6),nbvar(30),nbas1,nbast,nbvar1,
     *nblock,numvar,nrun,nitscf,nitsc1,nitsc2,nconv,nalarm,no,nstep,
     *iiii,nsht,nosht,nflag3,nbasm,mvar,iex,newbas,nit,nbpos(8)
      common/blkin/bias,dgath,scfat,dmnxp,rmned,rmxed,con1,con2,
     *nvar,ndiag,nxtrp,mxxtrp,norbit,nprint,mxvar,nextra,ifcont,irun,
     *nflag1,nflag2,nsvd,lim
c
      if (nalarm.le.0) then
         if (nprint.eq.0) go to 60
      end if
      write (iwr,6010) (s(i),i=1,nit)
      write (iwr,6030) (u(i),i=1,nit)
      write (iwr,6040) (t(i),i=1,nit)
      nn12 = nit*(nit+1)/2
      write (iwr,6050) (p(i),i=1,nn12)
      write (iwr,6060)
      k = 0
      do 40 i = 1 , nn12
         if (q(i).eq.0) go to 30
         k = k + 1
         p(k) = q(i)
         if (k.ne.10) go to 30
 20      write (iwr,6020) (p(m),m=1,k)
         k = 0
         go to 40
 30      if (i.eq.nn12) then
            if (k.gt.0) go to 20
         end if
 40   continue
c
      write (iwr,6070) (s(i),i=1,nit)
      do 50 i = 1 , nit
         u(i) = -czn*u(i) + t(i)
 50   continue
      write (iwr,6080) (u(i),i=1,nit)
      write (iwr,6090) (pcap(i),i=1,nit)
      write (iwr,6100) (qcap(i),i=1,nit)
      write (iwr,6110) (fc(i),i=1,nit)
      write (iwr,6120) (fo(i),i=1,nit)
      write (iwr,6130) (dt(i),i=1,nit)
      write (iwr,6140) (dos(i),i=1,nit)
c
      write (iwr,6150)
      write (iwr,6160)
      if (nsvd.eq.1) then
         write (iwr,6170)
      else
         write (iwr,6180)
      end if
      write (iwr,6190) czn , con1 , con2
      if (ifcont.ne.0) then
         write (iwr,6210) (nsavf(i),i=1,nsym)
         write (iwr,6200) (nbas(i),i=1,nsym)
         go to 70
      end if
 60   write (iwr,6210) (nbas(i),i=1,nsym)
 70   write (iwr,6220) (ncsh(i),i=1,nsym)
      write (iwr,6230) (popul(i),i=1,nsym)
      write (iwr,6240)
      nbas1 = nsavf(1)
      if (nsym.gt.1) then
         nbas2 = nsavf(2)
         if (nsym.gt.2) then
c
            nbas3 = nsavf(3)
            if (nsym.gt.3) then
               nbas4 = nsavf(4)
               do 80 i = 1 , nbas4
                  j = nbas1 + i
                  k = j + nbas2
                  l = k + nbas3
                  write (iwr,6250) pqn(i) , expa(i) , pqn(j) , expa(j) ,
     +                             pqn(k) , expa(k) , pqn(l) , expa(l)
 80            continue
               if (nbas4.eq.nbas3) go to 160
               nba4 = nbas4 + 1
               go to 140
            else
               nba4 = 1
               go to 140
            end if
         else
            nba3 = 1
         end if
      else
         do 90 i = 1 , nbas1
            write (iwr,6250) pqn(i) , expa(i)
 90      continue
         go to 170
      end if
 100  do 110 i = nba3 , nbas2
         j = nbas1 + i
         write (iwr,6250) pqn(i) , expa(i) , pqn(j) , expa(j)
 110  continue
 120  if (nbas1.ne.nbas2) then
         nba2 = nbas2 + 1
         do 130 i = nba2 , nbas1
            write (iwr,6250) pqn(i) , expa(i)
 130     continue
      end if
      go to 170
 140  do 150 i = nba4 , nbas3
         j = nbas1 + i
         k = nbas1 + nbas2 + i
         write (iwr,6250) pqn(i) , expa(i) , pqn(j) , expa(j) , pqn(k) ,
     +                    expa(k)
 150  continue
 160  if (nbas3.eq.nbas2) go to 120
      nba3 = nbas3 + 1
      go to 100
 170  if (ifcont.ne.0) then
         write (iwr,6260)
         do 180 i = 1 , newbas
            n = nrow(i)
            write (iwr,6270) i , (ctran(k,i),k=1,n)
 180     continue
      end if
      write (iwr,6280) nitscf
      if (nflag2.lt.1) then
         write (iwr,6290)
         write (iwr,6300) energ , pot , cin , vir
      else if (nflag2.eq.1) then
         write (iwr,6310)
         write (iwr,6300) energ , pot , cin , vir , encor
      else
         write (iwr,6320)
         write (iwr,6300) energ , pot , cin , vir
      end if
      write (iwr,6330)
      write (iwr,6340) (eps(i),i=1,nsht)
      write (iwr,6350)
      do 190 i = 1 , nbasm
c         nbasm  is the re  quired  number  in place  of  nsht
         write (iwr,6360) (c(i,j),j=1,nsht)
 190  continue
      write (iwr,6370) (thre(j),j=1,nsht)
      write (iwr,6380) thrscf
c
      return
 6010 format (//1x,104('-')//5x,
     +'matrices and supermatrices of integrals between basis functions'
     + //5x,'s matrix'//(1x,10f12.6))
 6020 format (6x,10f12.6)
 6030 format (6x,'u matrix'/(1x,10f12.6))
 6040 format (6x,'t matrix'/(1x,10f12.6))
 6050 format (6x,'p supermatrix'/(1x,10f12.6))
 6060 format (6x,'q supermatrix'/)
 6070 format (//1x,104('-')//5x,'final matrices'//5x,
     +        's matrix'/(1x,10f12.6))
 6080 format (5x,'h matrix'/(1x,10f12.6))
 6090 format (5x,'p matrix'/(1x,10f12.6))
 6100 format (5x,'q matrix'/(1x,10f12.6))
 6110 format (5x,'closed shell f matrix '/(1x,10f12.6))
 6120 format (5x,'open shell f matrix '/(1x,10f12.6))
 6130 format (5x,'total density matrix '/(1x,10f12.6))
 6140 format (5x,'open shell density matrix '/(1x,10f12.6))
 6150 format (//1x,104('-')////5x,'final results'/5x,13('-')//)
 6160 format (5x,'type gauss orbitals')
 6170 format (5x,'jacobi diagonalization')
 6180 format (5x,'single vector diagonalization')
 6190 format (5x,'charge =',f10.6//6x,'cut off =',f10.6,3x,f10.6//,
     +        6x,'symmetry species',11x,'s',5x,'p',5x,'d')
 6200 format (6x,'contracted bas.functions =',i2,3(4x,i2))
 6210 format (6x,'number of basis functions=',i2,3(4x,i2))
 6220 format (6x,'number of closed shells =  ',4(i1,4x))
 6230 format (6x,'open shell occupation =  ',4(f4.1,2x))
 6240 format (/5x,
     + 'basis functions (principal quantum number, orbital exponents)'
     +   //15x,'s',28x,'p',28x,'d',28x,'f'/)
 6250 format (4(5x,f4.1,5x,f15.6))
 6260 format (/5x,'contraction coefficients (normalized)'/' ')
 6270 format (6x,i2,8f14.6)
 6280 format (/6x,'final scf results obtained at computation no ',i2)
 6290 format (6x,'total hf energy',4x,'potential energy',4x,
     +        'kinetic energy',4x,'virial theorem')
 6300 format (4x,5(e19.10))
 6310 format (6x,'total hf energy',4x,'potential energy',4x,
     +        'kinetic energy',4x,'virial theorem',4x,'cor. energy')
 6320 format (6x,'total chf energy',3x,'potential energy',4x,
     +        'kinetic energy',4x,'virial theorem')
 6330 format (/6x,
     + 'orbital energies,eigenvectors,threshold diagonalization and scf'
     + )
c      changing  nsht  to the sum (up to nsym)  of nbas(i) gives
c         the  virtual orbitals  as  well
 6340 format (2x,8(4x,f10.5))
 6350 format (' ')
 6360 format (2x,8f14.6)
 6370 format (2x,8(4x,e10.4))
 6380 format (5x,e11.4)
      end
**==paxis.f
      subroutine paxis(dint,zspi,cspi,nocco,norb,iwr)
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
      dimension dint(*),nocco(*),map(6)
c
      integer mp, lwordp, lvec, ied3, newbas, lenb
      integer nbfnd, lenbad, lentrd, len3, ngrpm
      integer idena, nwroot, nosec
      logical noscf
      common /limy/ mp(10),lwordp,lvec,ied3,newbas,lenb,
     +              nbfnd,lenbad,lentrd,len3,ngrpm,
     +              idena,noscf,nwroot,nosec(mxroot)
c
c
      real*8 cconv, atwt, cspin
      integer lin, kin, nin, ngmx, nfmx, nomx, ncmx, ntmx, ntcol, npmx
      integer icon, ncom, maxnt,itps, maxcon
c
      common /datprp/ cconv(11),atwt(30),cspin(30),
     + lin(84),kin(84),nin(84),
     + ngmx,nfmx,nomx,ncmx,ntmx,ntcol,npmx,
     + icon(24),ncom(22),maxnt,itps,maxcon
c
      common/junkc/ztype(435),zhead(5,22),zunit(3,11),ztit(220),
     * zspin(30),zsta,zhz,zgauss
      common/prpspc/cod(3),mm,ilnmcm
      common/miscop/p(3,3),prnct(3,3),prnc(3,3),q(3,3),buf(6),
     1u(3,3),iq(3,3)
      data m0,m1,mm1,m2,map/0,1,-1,2,1,4,5,2,6,3/
      data m3,nr/3,3/
      data dzero/0.0d0/
c
      ij = m0
      do 30 i = 1 , 3
         do 20 j = i , 3
            ij = ij + m1
            prnct(i,j) = dint(mp(map(ij))+m2)
 20      continue
 30   continue
      call hdiag(prnct,u,p,iq,3,3,0,nr)
      ij = m0
      do 50 i = 1 , 3
         do 40 j = i , 3
            ij = ij + m1
            dint(mp(map(ij))+m2) = prnct(i,j)
 40      continue
 50   continue
      call mwt(u,m3,m3,m3)
      do 120 n = 1 , 2
         nn = n + mm1
         ij = m0
         do 70 i = 1 , 3
            do 60 j = i , 3
               ij = ij + m1
               prnc(i,j) = dint(mp(map(ij))+nn)
 60         continue
 70      continue
         call mabat(u,prnc,m3,m3,q,m3)
         ij = m0
         do 90 i = 1 , 3
            do 80 j = i , 3
               ij = ij + m1
               dint(mp(map(ij))+nn) = prnc(i,j)
 80         continue
 90      continue
         if (n.eq.1) then
            do 100 m = 1 , 6
               buf(m) = dint(mp(m))
 100        continue
            write (iwr,6080) buf
         else
            do 110 m = 1 , 6
               buf(m) = dint(mp(m)+m1)
 110        continue
            write (iwr,6090) buf
         end if
 120  continue
      do 130 m = 1 , 6
         buf(m) = dint(mp(m)+m2)
 130  continue
      write (iwr,6100) buf
      if (mm.eq.4) then
         pp = dabs(buf(1))
         qq = dabs(buf(2))
         rr = dabs(buf(3))
         if (pp.lt.qq) then
            temp = qq
            qq = pp
            pp = temp
         end if
         if (qq.lt.rr) then
            temp = rr
            rr = qq
            qq = temp
            if (pp.lt.qq) then
               temp = qq
               qq = pp
               pp = temp
            end if
         end if
         if (pp.ne.dzero) then
            asym = (qq-rr)/pp
            write (iwr,6010) asym
         end if
      end if
      if (mm.ne.20) then
         cc = cconv(mm)
         do 140 i = 1 , 6
            buf(i) = buf(i)*cc
 140     continue
         write (iwr,6020) (zunit(mmm,mm),mmm=1,3) , (buf(j),j=1,6)
      else
         bot = 1600.48d0*cspi*0.119366207d0
c ... 3/8pi  =  .119366207
         do 150 i = 1 , 6
            buf(i) = buf(i)*bot
 150     continue
         write (iwr,6030) zspi , cspi
         write (iwr,6040) zhz , (buf(j),j=1,6)
         do 160 i = 1 , 6
            buf(i) = buf(i)*.356835d0
 160     continue
         write (iwr,6040) zgauss , (buf(i),i=1,6)
      end if
      write (iwr,6070)
      n = 3
      do 230 i = 1 , norb
         do 220 j = 1 , i
            if (icon(8).le.0) then
               if (i.ne.j) then
                  n = n + m1
                  go to 220
               else if ((nocco(i)+nocco(j)).ne.2) then
                  n = n + m1
                  go to 220
               end if
            else if ((nocco(i)+nocco(j)).lt.2) then
               n = n + m1
               go to 220
            end if
            ij = m0
            do 180 k = 1 , 3
               do 170 l = k , 3
                  ij = ij + m1
                  prnc(k,l) = dint(mp(map(ij))+n)
 170           continue
 180        continue
            call mabat(u,prnc,m3,m3,q,m3)
            ij = m0
            do 200 k = 1 , 3
               do 190 l = k , 3
                  ij = ij + m1
                  dint(mp(map(ij))+n) = prnc(k,l)
 190           continue
 200        continue
            do 210 k = 1 , 6
               buf(k) = dint(n+mp(k))
 210        continue
            write (iwr,6110) i , j , buf
            n = n + m1
 220     continue
 230  continue
      write (iwr,6060)
      do 240 n = 1 , 3
         pp = dint(mp(n)+m2)
         write (iwr,6050) pp , (u(n,j),j=1,3)
 240  continue
c
      return
 6010 format (/1x,'asymmetry parameter =',f10.3)
 6020 format (/50x,'in ',3a8//1x,'total',4x,10f12.5)
 6030 format (//2x,'isotope',6x,'nuclear'/15x,'moment (u/2i)'//5x,a4,4x,
     +        f15.7/)
 6040 format (/50x,'in ',a5//1x,'total',4x,10f12.5)
c
 6050 format (f12.6,5x,3(f12.6,3x))
 6060 format (//1x,'rotation matrix'//3x,'eigenvalue',23x,
     +        'eigenvectors'/)
 6070 format (/40x,'integrals over molecular orbitals'//)
 6080 format (1x,'nuclear',2x,10f12.5)
 6090 format (1x,'electron',1x,10f12.5)
 6100 format (1x,'total',4x,10f12.5)
 6110 format (2i5,10f12.5)
c
      end
**==potdm.f
      subroutine potdm(vect,pop,dens,ncol,norb)
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      dimension vect(*),pop(*),dens(*)
c
c only for hfscf at present
c
      call dmtx2(dens,vect,pop,iky,ncol,norb,norb)
      return
      end
**==potnuc.f
      subroutine potnuc(noc,vlist,c,dint,nocmx)
c...     nuclear repulsion
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension dint(*),vlist(nocmx,4),c(3)
c
      pnuc = 0.0d0
      do 20 i = 1 , noc
         dst = (vlist(i,1)-c(1))**2 + (vlist(i,2)-c(2))
     +         **2 + (vlist(i,3)-c(3))**2
         if (dst.gt.1.0d-10) then
            pnuc = pnuc + vlist(i,4)/dsqrt(dst)
         end if
 20   continue
      dint(1) = pnuc
c
      return
      end
**==prnt.f
      subroutine prnt(dint,zspi,cspi,norb,iwr)
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
      dimension dint(*),zglob(4),zprin(4)
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
c
      real*8 cconv, atwt, cspin
      integer lin, kin, nin, ngmx, nfmx, nomx, ncmx, ntmx, ntcol, npmx
      integer icon, ncom, maxnt,itps, maxcon
c
      common /datprp/ cconv(11),atwt(30),cspin(30),
     + lin(84),kin(84),nin(84),
     + ngmx,nfmx,nomx,ncmx,ntmx,ntcol,npmx,
     + icon(24),ncom(22),maxnt,itps,maxcon
c
      common/junk/eta(2,mxgaus),vlist(maxat,4),frocc(maxorb),
     *nfirst(4,maxorb),isky(mxgaus),nocco(maxorb),
     *noc,non
      common/junkc/ztype(435),zhead(5,22),zunit(3,11),ztitle(10,22)
c
      integer mp, lwordp, lvec, ied3, newbas, lenb
      integer nbfnd, lenbad, lentrd, len3, ngrpm
      integer idena, nwroot, nosec
      logical noscf
      common /limy/ mp(10),lwordp,lvec,ied3,newbas,lenb,
     +              nbfnd,lenbad,lentrd,len3,ngrpm,
     +              idena,noscf,nwroot,nosec(mxroot)
c
      common/miscop/buf(10),buf1(10),buf2(10),buf3(10)
      common/prpspc/c(3),mm,ilnmcm,nt
      data m1,m2/1,2/
      data zdash/'------'/
      data zglob,zprin/
     *'global c','oordinat','es',' ',
     *'principa','l axis c','oordinat','es'/
c...print out nuclear and electronic contributions and total value.
      write (iwr,6070) zglob
      write (iwr,6080) (ztitle(i,mm),i=1,nt)
      nt2 = nt + nt
      write (iwr,6090) (zdash,i=1,nt2)
      do 20 i = 1 , nt
         ii4 = mp(i)
         buf(i) = dint(ii4)
         buf1(i) = dint(ii4+m1)
         buf2(i) = dint(ii4+m2)
 20   continue
      write (iwr,6100) (buf(i),i=1,nt)
      write (iwr,6110) (buf1(i),i=1,nt)
      write (iwr,6120) (buf2(i),i=1,nt)
      if (mm.ne.19) then
         if (mm.eq.3) then
            if (ilnmcm.le.non) then
               do 30 i = 1 , nt
                  buf3(i) = buf2(i)*vlist(ilnmcm,4)
 30            continue
            end if
         end if
         if (mm.le.11) then
            if (mm.ne.4 .and. mm.ne.6) then
               cc = cconv(mm)
               do 40 i = 1 , nt
                  buf2(i) = buf2(i)*cc
 40            continue
               write (iwr,6010) (zunit(i,mm),i=1,3) , (buf2(j),j=1,nt)
               if (mm.eq.3) then
                  if (ilnmcm.le.non) then
                     write (iwr,6020) (buf3(i),i=1,nt)
                  end if
               end if
            end if
         end if
      else
         bot = 1600.48d0*cspi*buf2(1)
         hot = bot*.356835d0
         write (iwr,6030) zspi , cspi , buf2(1) , bot , hot
      end if
      if (icon(8).gt.0) then
c...
c...print integrals over molecular orbitals.
c...
c     reset format flag for high precision on print
         otemp = oprint(20)
         oprint(20) = .true.
         write (iwr,6130)
         do 50 loop = 1 , nt
            write (iwr,6040) ztitle(loop,mm)
            call writel(dint(mp(loop)+3),norb)
 50      continue
         oprint(20) = otemp
      else
         write (iwr,6050) (ztitle(i,mm),i=1,nt)
         write (iwr,6090) (zdash,i=1,nt2)
         mu = 3
         do 80 i = 1 , norb
            do 70 j = 1 , i
               if ((nocco(i)+nocco(j)).ge.2 .and. i.eq.j) then
                  do 60 k = 1 , nt
                     ii4 = mp(k)
                     buf(k) = dint(ii4+mu)
 60               continue
                  write (iwr,6140) i , j , (buf(k),k=1,nt)
               end if
               mu = mu + 1
 70         continue
 80      continue
      end if
      if (mm.eq.4 .or. mm.eq.6 .or. mm.eq.8 .or. mm.eq.20) then
c...
c...perform principal axis transformation.
c...
         write (iwr,6060)
         write (iwr,6070) zprin
         write (iwr,6080) (ztitle(i,mm),i=1,6)
         write (iwr,6090) (zdash,i=1,nt2)
         call paxis(dint,zspi,cspi,nocco,norb,iwr)
         write (iwr,6060)
c
         return
      else
         write (iwr,6060)
         return
      end if
 6010 format (//50x,'in ',3a8//1x,'total',4x,10f12.5)
 6020 format (/1x,'force at nucleus (a.u.)'//10x,3f12.7)
 6030 format (//2x,'isotope',6x,'nuclear',9x,'spin density',3x,
     +        'hyperfine splitting'/15x,'moment (u/2i)',25x,'mc/s',5x,
     +        'gauss'//5x,a4,4x,f15.7,2x,f15.7,2x,2f10.2)
 6040 format (/1x,104('=')//40x,'component : ',a6)
 6050 format (//6x,'mos',1x,10(6x,a6))
 6060 format (/)
c
 6070 format (/1x,34('=')//' expectation values in atomic units'/' --- '
     +        ,4a8/1x,34('=')//)
 6080 format (/10x,10(6x,a6))
 6090 format (10x,20a6)
 6100 format (1x,'nuclear',2x,10f12.5)
 6110 format (1x,'electron',1x,10f12.5)
 6120 format (1x,'total',4x,10f12.5)
 6130 format (/40x,'integrals over molecular orbitals (a.u.)'/40x,
     +        '########################################'/)
 6140 format (2i5,10f12.5)
c
      end
**==prod.f
      subroutine prod(cc,isky1,dd,isky2,s,mrange)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension cc(*),dd(*),s(*),c(10),d(10),ss(35)
      common/miscop/aaa(70),
     *c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,
     *d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,
     *s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,
     *s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,
     *s21,s22,s23,s24,s25,s26,s27,s28,s29,s30,
     *s31,s32,s33,s34,s35
      equivalence(c(1),c1),(d(1),d1),(ss(1),s1)
      data m10,m20,m35/10,20,35/
c
      id = 0
      call vclr(ss,1,35)
      call dcopy(m10,cc(1),1,c(1),1)
      call dcopy(m10,dd(1),1,d(1),1)
      s1 = c1*d1
      s2 = c2*d1 + c1*d2
      s3 = c3*d1 + c1*d3
      s4 = c4*d1 + c1*d4
      s5 = c2*d2
      s6 = c3*d3
      s7 = c4*d4
      s8 = c3*d2 + c2*d3
      s9 = c4*d2 + c2*d4
      s10 = c4*d3 + c3*d4
      mrange = m10
      if (isky1.gt.4) then
         id = id + 1
         do 20 i = 5 , 10
            ss(i) = ss(i) + c(i)*d1
 20      continue
         s11 = c5*d2
         s12 = c6*d3
         s13 = c7*d4
         s14 = c5*d3 + c8*d2
         s15 = c5*d4 + c9*d2
         s16 = c6*d2 + c8*d3
         s17 = c6*d4 + c10*d3
         s18 = c7*d2 + c9*d4
         s19 = c7*d3 + c10*d4
         s20 = c8*d4 + c9*d3 + c10*d2
         mrange = m20
      end if
      if (isky2.gt.4) then
         id = id + 1
         do 30 i = 5 , 10
            ss(i) = ss(i) + d(i)*c1
 30      continue
         s11 = s11 + c2*d5
         s12 = s12 + c3*d6
         s13 = s13 + c4*d7
         s14 = s14 + c3*d5 + c2*d8
         s15 = s15 + c4*d5 + c2*d9
         s16 = s16 + c2*d6 + c3*d8
         s17 = s17 + c4*d6 + c3*d10
         s18 = s18 + c2*d7 + c4*d9
         s19 = s19 + c3*d7 + c4*d10
         s20 = s20 + c4*d8 + c3*d9 + c2*d10
         mrange = m20
      end if
      if (id.ge.2) then
         s21 = c5*d5
         s22 = c6*d6
         s23 = c7*d7
         s24 = c8*d5 + c5*d8
         s25 = c9*d5 + c5*d9
         s26 = c8*d6 + c6*d8
         s27 = c10*d6 + c6*d10
         s28 = c9*d7 + c7*d9
         s29 = c10*d7 + c7*d10
         s30 = c6*d5 + c5*d6 + c8*d8
         s31 = c7*d5 + c5*d7 + c9*d9
         s32 = c7*d6 + c6*d7 + c10*d10
         s33 = c9*d8 + c8*d9 + c10*d5 + c5*d10
         s34 = c9*d6 + c6*d9 + c10*d8 + c8*d10
         s35 = c8*d7 + c7*d8 + c10*d9 + c9*d10
         mrange = m35
      end if
      call dcopy(m35,ss(1),1,s(1),1)
c
      return
      end
**==propdr.f
      subroutine propdr(iskyx,iskyy,iskyz,q,dint11,dint)
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
      dimension dint(*),q(*),dint11(*),iskyx(*),iskyy(*),iskyz(*)
      dimension jskx(10),jsky(10),jskz(10),ccc(3),bbb(80)
      dimension vvale(12),e(3)
c
      real*8 cconv, atwt, cspin
      integer lin, kin, nin, ngmx, nfmx, nomx, ncmx, ntmx, ntcol, npmx
      integer icon, ncom, maxnt,itps, maxcon
c
      common /datprp/ cconv(11),atwt(30),cspin(30),
     + lin(84),kin(84),nin(84),
     + ngmx,nfmx,nomx,ncmx,ntmx,ntcol,npmx,
     + icon(24),ncom(22),maxnt,itps,maxcon
c
      common/junk/eta(mxgaus),wnorm(mxgaus),vlist(maxat,4),
     *frocc(maxorb),
     *nfirst(maxorb),nlast(maxorb),nctr(maxorb),ntype(maxorb),
     *isky(mxgaus),nocco(maxorb),
     *noc,non,maxtyp,ngaus,nbfns,newb,
     *ncol,ngrps,idhar,maxcom
      common/junkc/ztype(35),zlist(300),znumb(100),zhead(5,22),
     * zunit(33),zbuff(10,22),zspin(30),zsta,zhz,zgauss
      common/reso/iresc(100),jresc(100),istat(100),ipos(100)
      common/sver/la,lb,lc,ntt,gama,v(10),d(3)
c
      integer mp, lwordp, lvec, ied3, newbas, lenb
      integer nbfnd, lenbad, lentrd, len3, ngrpm
      integer idena, nwroot, nosec
      logical noscf
      common /limy/ mp(10),lwordp,lvec,ied3,newbas,lenb,
     +              nbfnd,lenbad,lentrd,len3,ngrpm,
     +              idena,noscf,nwroot,nosec(mxroot)
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
      integer mapsie
      common /linkan/ mapsie(maxorb)
c
c
      real*8 ctran
      integer ilifc, ntran, itran
      logical otran,otri
      common/tran/ilifc(maxorb),ntran(maxorb),itran(mxorb3),
     +            ctran(mxorb3),otran,otri
c
      common/miscop/ext(125),a(3),b(3),ee(3),aa(10),aa1(10),
     *bb(10),bb1(10),bb2(10),bb3(10),dd(35),
     *dd1(35),dd2(35),dd3(35),vval(10),vual(3),
     *vval1(3),vval2(3),vval3(3),dimt(10),
     *cbuf(10)
      common/prpspc/centx,centy,centz,mop,icent,nt
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
      common/restri/jjfile(63),lds(508),isect(508),ldsect(508)
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
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
c
      character*10 charwall
      equivalence (centx,ccc(1)),(bb(1),bbb(1)),(vvale(1),vual(1))
c
      data mm1,m0,m1,m20,m30/
     *-1,0,1,20,30/
      data m10/10/
      data jskx,jsky,jskz/
     *0,1,0,0,2,0,0,3,4,0,
     *0,0,1,0,0,3,0,2,0,4,
     *0,0,0,1,0,0,4,0,2,3/
c
      ind(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)
c
      nav = lenwrd()
c
      timest = cpulft(1)
      npmx2 = (4*npmx)/nav
      timax = 0.0d0
      write (iwr,6010) timest
      do 720 idr = 1 , npmx
         istatu = istat(idr)
         if (istatu.le.0) go to 720
         mop = iresc(idr)
         consf = 1.0d0
         if (mop.gt.m20) consf = -consf
         itpro = mop + itps
         icent = jresc(idr)
         iposi = ipos(idr)
         zcentg = zaname(icent)
         centx = vlist(icent,1)
         centy = vlist(icent,2)
         centz = vlist(icent,3)
         centc = vlist(icent,4)
         i = centc + .005d0
         if (i.gt.m0 .and. i.le.m30) then
            cspi = cspin(i)
            zspi = zspin(i)
         else
            cspi = 0.0d0
            zspi = zsta
         end if
         write (iwr,6020) (zhead(ki,mop),ki=1,5) , zcentg , ccc
         nt = ncom(mop)
         if (nt.le.maxcom) then
            lenp = lentrd*nt + m1
            do 20 jdr = 1 , npmx
               l = ipos(jdr)
               if (iresc(jdr).eq.mop .and. jresc(jdr).eq.icent .and.
     +             l.gt.m0) go to 30
 20         continue
c
            go to (40,50,60,70,80,90,100,110,120,130,140,150,160,170,
     +             180,180,180,180,180,180,180,180) , mop
            go to 40
         else
            write (iwr,6030) (zhead(ki,mop),ki=1,5) , zcentg
            go to 720
         end if
 30      call getpro(dint,zcentg,l,itpro,mop,iwr)
         go to 650
 40      call potnuc(non,vlist,ccc,dimt,maxat)
         go to 190
 50      call sldnuc(non,vlist,ccc,dimt,maxat)
         go to 190
 60      call fldnuc(non,vlist,ccc,dimt,maxat)
         go to 190
 70      call grdnuc(non,vlist,ccc,dimt,maxat)
         go to 190
 80      call dipnuc(non,vlist,ccc,dimt,maxat)
         go to 190
 90      call qudnuc(non,vlist,ccc,dimt,maxat)
         go to 190
 100     call agnuc(non,vlist,ccc,dimt,maxat)
         go to 190
 110     call anuc2(non,vlist,ccc,dimt,maxat)
         go to 190
 120     call octnuc(non,vlist,ccc,dimt,maxat)
         go to 190
 130     call anuc3(non,vlist,ccc,dimt,maxat)
         go to 190
 140     call bnuc3(non,vlist,ccc,dimt,maxat)
         go to 190
 150     call hexnuc(non,vlist,ccc,dimt,maxat)
         go to 190
 160     call anuc4(non,vlist,ccc,dimt,maxat)
         go to 190
 170     call bnuc4(non,vlist,ccc,dimt,maxat)
         go to 190
 180     call onuc(dimt)
c
 190     do 200 i = 1 , nt
            dint(mp(i)) = dimt(i)
 200     continue
         if (mop.gt.m20) then
c ... additional code for l and l/rrr integrals
c ... the integrals are evaluated
c ...      <i/op/j>
c ... where i=1,nbfns   ;   j=1,i
            do 340 i = 1 , nbfns
               iiii = mapsie(i)
               is = nfirst(iiii)
               if = nlast(iiii)
               m = nctr(iiii)
               do 210 j = 1 , 3
                  a(j) = vlist(m,j)
 210           continue
               j = ntype(iiii)
               ki = jskx(j)
               kj = jsky(j)
               ii = jskz(j)
               do 220 j = is , if
                  iskyx(j) = ki
                  iskyy(j) = kj
                  iskyz(j) = ii
 220           continue
               do 330 j = 1 , i
                  jjjj = mapsie(j)
                  do 230 m = 1 , 3
                     cbuf(m) = 0.0d0
 230              continue
                  js = nfirst(jjjj)
                  jf = nlast(jjjj)
                  m = nctr(jjjj)
                  nu = 2 + ind(iiii,jjjj)
                  do 240 jj = 1 , 3
                     b(jj) = vlist(m,jj)
 240              continue
                  do 310 ii = is , if
                     isky1 = isky(ii)
                     call vclr(aa1,1,10)
                     aa1(isky1) = wnorm(ii)
                     dumx = eta(ii)
                     do 300 jj = js , jf
                        wj = wnorm(jj)
                        dumy = eta(jj)
                        isky2 = isky(jj)
                        do 250 m = 1 , 3
                           e(m) = b(m) - ccc(m)
 250                    continue
                        wx = lin(isky2)
                        wy = kin(isky2)
                        wz = nin(isky2)
                        isky3 = iskyx(jj)
                        isky4 = iskyy(jj)
                        isky5 = iskyz(jj)
                        do 260 m = 1 , 40
                           bbb(m) = 0.0d0
 260                    continue
                        bb(isky2) = wj*dumy*2.0d0
                        if (isky3.gt.0) bb1(isky3) = wj*wx
                        if (isky4.gt.0) bb2(isky4) = wj*wy
                        if (isky5.gt.0) bb3(isky5) = wj*wz
                        call divpt(a,dumx,b,dumy,ee,gama,efact)
                        call dcopy(m10,aa1(1),1,aa(1),1)
                        call shftce(aa,a,ee,isky1)
                        call shftce(bb,b,ee,isky2)
                        call shftce(bb1,b,ee,isky3)
                        call shftce(bb2,b,ee,isky4)
                        call shftce(bb3,b,ee,isky5)
                        call prod(aa,isky1,bb,isky2,dd,mrange)
                        call prod(aa,isky1,bb1,isky3,dd1,mrange)
                        call prod(aa,isky1,bb2,isky4,dd2,mrange)
                        call prod(aa,isky1,bb3,isky5,dd3,mrange)
                        call vclr(vvale,1,12)
                        do 270 ki = 1 , 3
                           d(ki) = ee(ki) - ccc(ki)
 270                    continue
                        do 290 ki = 1 , 35
                           la = lin(ki)
                           lb = kin(ki)
                           lc = nin(ki)
                           dumy = dabs(dd(ki)) + dabs(dd1(ki))
     +                            + dabs(dd2(ki)) + dabs(dd3(ki))
                           if (dumy.ge.1.d-8) then
                              if (mop.le.21) then
                                 call opab1(la,lb,lc,gama,v,ntt,d)
                              else
                                 call opaa1(la,lb,lc,gama,v,ntt,d)
                              end if
                              do 280 kj = 1 , nt
                                 dumy = v(kj)
                                 vual(kj) = dd(ki)*dumy + vual(kj)
                                 vval1(kj) = dd1(ki)*dumy + vval1(kj)
                                 vval2(kj) = dd2(ki)*dumy + vval2(kj)
                                 vval3(kj) = dd3(ki)*dumy + vval3(kj)
 280                          continue
                           end if
 290                    continue
                        cbuf(1) = cbuf(1)
     +                            + efact*(e(3)*vual(2)-e(2)*vual(3)
     +                            +vval3(2)-vval2(3))
                        cbuf(2) = cbuf(2)
     +                            + efact*(e(1)*vual(3)-e(3)*vual(1)
     +                            +vval1(3)-vval3(1))
                        cbuf(3) = cbuf(3)
     +                            + efact*(e(2)*vual(1)-e(1)*vual(2)
     +                            +vval2(1)-vval1(2))
 300                 continue
 310              continue
                  do 320 ki = 1 , nt
                     dumy = cbuf(ki)
                     if (mop.eq.22) dumy = -dumy
                     dint(nu+mp(ki)) = dumy
 320              continue
 330           continue
 340        continue
         else
c
c ----- loops over internal ordered basis functions
c ----- integrals stored over extranal functions ie
c ----- identical ordering for atmol, revised ordering for gamess
c
            do 640 i = 1 , nbfns
               iiii = mapsie(i)
               is = nfirst(iiii)
               if = nlast(iiii)
               m = nctr(iiii)
               do 350 j = 1 , 3
                  a(j) = vlist(m,j)
 350           continue
               do 630 j = 1 , i
                  jjjj = mapsie(j)
                  js = nfirst(jjjj)
                  jf = nlast(jjjj)
                  m = nctr(jjjj)
                  do 360 jj = 1 , 3
                     b(jj) = vlist(m,jj)
 360              continue
                  nu = 2 + ind(jjjj,iiii)
                  call vclr(cbuf,1,10)
                  do 610 ii = is , if
                     isky1 = isky(ii)
                     call vclr(aa1,1,10)
                     aa1(isky1) = wnorm(ii)
                     dumx = eta(ii)
                     do 600 jj = js , jf
                        isky2 = isky(jj)
                        call vclr(bb,1,10)
                        bb(isky2) = wnorm(jj)
                        dumy = eta(jj)
                        call dcopy(m10,aa1(1),1,aa(1),1)
                        call divpt(a,dumx,b,dumy,ee,gama,efact)
                        call shftce(aa,a,ee,isky1)
                        call shftce(bb,b,ee,isky2)
                        call prod(aa,isky1,bb,isky2,dd,mrange)
                        call vclr(vval,1,10)
                        do 370 ki = 1 , 3
                           d(ki) = ee(ki) - ccc(ki)
 370                    continue
                        do 580 ki = 1 , mrange
                           if (dabs(dd(ki)).le.1.0d-8) go to 580
                           la = lin(ki)
                           lb = kin(ki)
                           lc = nin(ki)
c
                           go to (380,390,400,410,420,430,440,450,460,
     +                            470,480,490,500,510,520,530,540,550,
     +                            550,410) , mop
 380                       call opaa0(la,lb,lc,gama,v,ntt,d)
                           go to 560
 390                       call opaa3(la,lb,lc,gama,v,ntt,d)
                           go to 560
 400                       call opaa1(la,lb,lc,gama,v,ntt,d)
                           go to 560
 410                       call opaa2(la,lb,lc,gama,v,ntt,d)
                           go to 560
 420                       call opab1(la,lb,lc,gama,v,ntt,d)
                           go to 560
 430                       call opab2(la,lb,lc,gama,v,ntt,d)
                           go to 560
 440                       call opab4(la,lb,lc,gama,v,ntt,d)
                           go to 560
 450                       call opab6(la,lb,lc,gama,v,ntt,d)
                           go to 560
 460                       call opab3(la,lb,lc,gama,v,ntt,d)
                           go to 560
 470                       call opab7(la,lb,lc,gama,v,ntt,d)
                           go to 560
 480                       call opab8(la,lb,lc,gama,v,ntt,d)
                           go to 560
 490                       call opab5(la,lb,lc,gama,v,ntt,d)
                           go to 560
 500                       call opab9(la,lb,lc,gama,v,ntt,d)
                           go to 560
 510                       call opab10(la,lb,lc,gama,v,ntt,d)
                           go to 560
 520                       call opad1(la,lb,lc,gama,v,ntt)
                           go to 560
 530                       call opac1(la,lb,lc,gama,v,ntt,d)
                           go to 560
 540                       call opac2(la,lb,lc,gama,v,ntt,d)
                           go to 560
 550                       call opac3(la,lb,lc,gama,v,ntt,d)
c
 560                       do 570 kj = 1 , nt
                              vval(kj) = dd(ki)*v(kj) + vval(kj)
 570                       continue
 580                    continue
                        do 590 ki = 1 , nt
                           cbuf(ki) = cbuf(ki) + efact*vval(ki)
 590                    continue
 600                 continue
 610              continue
                  do 620 ki = 1 , nt
                     dint(nu+mp(ki)) = cbuf(ki)
 620              continue
 630           continue
 640        continue
         end if
 650     if (iposi.lt.0) then
            iposi = -iposi
            call secput(iposi,itpro,lenp,iblpro)
            write (iwr,6040) iposi , iblpro
            call putpro(dint,zhead,zbuff,zcentg,centx,centy,centz,
     +                  centc,zspi,cspi,mop,nt,iblpro)
            call revind
         end if
         ipos(idr) = iposi
         istat(idr) = mm1
         if (icon(9).gt.0) then
            write (iwr,6050) (zbuff(ki,mop),ki=1,nt)
            write (iwr,6060)
            mu = 3
c
c ----- loop over external basis functions
c
            do 680 ki = 1 , nbfnd
               do 670 kj = 1 , ki
                  do 660 k = 1 , nt
                     aa(k) = dint(mp(k)+mu)
 660              continue
                  mu = mu + 1
                  write (iwr,6070) ki , kj , (aa(k),k=1,nt)
 670           continue
 680        continue
         end if
         do 690 ki = 1 , nt
            ii4 = mp(ki) + 3
            call tranpr(dint(ii4),dint11,consf,newb,ilifc,ntran,
     +                  itran,ctran,otran,lenb)
 690     continue
         call contra(dint,dint11,q,ncol,newb,nt,consf)
c .. compute electronic+total contribution
         do 710 i = 1 , nt
            wx = 0.0d0
            ii4 = mp(i)
            mu = 2
            do 700 j = 1 , ncol
               mu = mu + j
               if (frocc(j).ne.0.0d0) then
                  wx = wx + frocc(j)*dint(ii4+mu)
               end if
 700        continue
            dint(ii4+1) = wx
            dint(ii4+2) = wx + dint(ii4)
 710     continue
         call prnt(dint(1),zspi,cspi,ncol,iwr)
         top = cpulft(1)
         call secput(isect(468),itps,m1,ibl468)
         call wrt3i(iresc,npmx2*nav,ibl468,idaf)
         call revind
         if ((top-timest).ge.timax) then
            timax = top - timest
         end if
         if (timlim-timax*1.2d0.lt.top) go to 730
         timest = top
         write (iwr,6080) top,charwall()
 720  continue
      irest = 0
c
      write (iwr,6090)
c
      return
c
 730  irest = 7
      call texit(0,irest)
      tim = timlim + 0.5d0
      return
 6010 format (/1x,104('*')//' commence property evaluation at ',f9.2,
     +        ' secs.')
 6020 format (/1x,104('-')//50x,30('*')/50x,5a6/50x,30('*')///50x,
     +        'calculate with reference to the point ',5x,a8//50x,
     +        'coordinates'/50x,11('-')//45x,'x',11x,'y',11x,'z'/39x,
     +        3f12.6)
 6030 format (//10x,'----- insufficient core allocated to handle'/10x,
     +        '      ',5a6,' at point ',a8//10x,
     +        '-----  resubmit with increased core'/10x,
     +        '      as specified'//)
 6040 format (/' integrals output to dumpfile :section ',i3,' (block ',
     +        i6,')')
 6050 format (//40x,'integrals over gtos (a.u.)'//5x,'gtos ',10(6x,a6))
 6060 format (/)
 6070 format (2i5,10(1x,f11.4))
 6080 format (/' end of property evaluation at ',f9.2,' seconds',
     +        a10,' wall'/)
 6090 format (//1x,104('*')/)
      end
**==propg.f
      subroutine propg
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
      dimension ngrp(3),minf(3),maxf(3)
      dimension c1(mxprms),alp(mxprms)
      common/restri/jjfile(63),lds(508),isect(508),ldsect(508)
c
      integer ltask, lcoord, iplot, iprop, imull, idma, inbo
      logical opass1, opass2, opass3, orege, opass4, offsym
      logical opass5, opass6, opass8, opass9, odisc, oirc
      logical orgall, orestj, opass10, opass11
      logical omrdci, ofulci, omech, omopac
      logical orpa, omclr, odirpa, occsd, occsdt, oqcisd, oqcisdt
      logical omrdcid,opark,oclunp,ofill,odiesel
c
      common /restrj/ ltask,lcoord,opass1,opass2,opass3,orege,iplot,
     + iprop,imull,opass4,offsym,opass5,opass6,opass10,idma,inbo,
     + omrdci,ofulci,omech,omopac,opass8,opass9,odisc,
     + oirc,orpa,omclr,odirpa,occsd,occsdt,oqcisd,oqcisdt,orgall,
     + omrdcid,orestj(7),opark,oclunp,ofill(6),odiesel, opass11
c
      common/junk/eta(mxgaus),wnorm(mxgaus),vlist(maxat,4),
     +frocc(maxorb),nfirst(maxorb),nlast(maxorb),nctr(maxorb),
     +ntype(maxorb),isky(mxgaus),nocco(maxorb),noc,mon,maxtyp,kk,ktw,
     +newb,ncol,mgrps,idhar,maxcom,
     +kstart(mxshel),katom(mxshel),ktype(mxshel),
     +kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),nshell,non3(3)
      common/junkc/ztita(10),ztitb(maxat)
      common/miscop/
     *ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),cf(mxprim),cg(mxprim)
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
      integer mapsie
      common /linkan/ mapsie(maxorb)
c
c
      integer mp, lwordp, lvec, ied3, newbas, lenb
      integer nbfnd, lenbad, lentrd, len3, ngrpm
      integer idena, nwroot, nosec
      logical noscf
      common /limy/ mp(10),lwordp,lvec,ied3,newbas,lenb,
     +              nbfnd,lenbad,lentrd,len3,ngrpm,
     +              idena,noscf,nwroot,nosec(mxroot)
c
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
      data mone/1/
      data ngrp/1,3,6/
      data minf/1,2,5/
      data maxf/1,4,10/
      data pi32/5.56832799683170d0/
      data pt5,pt75,dzero/0.5d0,0.75d0,0.0d0/
      data sqrt3/1.73205080756888d0/
c
c ----- retrieve basis set from section 491
c
      nav = lenwrd()
c
      call secget(isect(491),mone,iblkv)
      m1420 = 5*mxprim + maxat
      m110 = 10 + maxat
      call rdedx(ex,mxprim,iblkv,idaf)
      call reads(cs,m1420,idaf)
      call rdchrs(ztita,m110,idaf)
      call readis(kstart,mach(2)*nav,idaf)
c
c ----- position ed7 for retrieval of analysis data 
c
      nwp2 = 2*maxorb + 82 + 60/nav
      nwdma = 5 + 4*maxat + (13+maxat+nav)/nav
      i = 1
      if (idma.gt.0) i = i + lensec(maxat) + lensec(nwdma)
      if (iplot.gt.0) i = i + iplot*(lensec(nwp2)+lensec(140))
      call search(i,num8)
c
 20   if (nat.le.0 .or. nat.gt.maxat) then
         call caserr(
     +       'invalid parameters detected in analysis preprocessor')
      else if (nshell.le.0 .or. nshell.gt.mxshel) then
         call caserr(
     +       'invalid parameters detected in analysis preprocessor')
      else if (num.le.0 .or. num.gt.maxorb) then
         call caserr(
     +       'invalid parameters detected in analysis preprocessor')
      else
         idhar = 0
         noc = nat
         mon = noc
         kk = 0
         maxtyp = 4
         kls = 0
         ktw = 0
         do 30 i = 1 , nat
            vlist(i,1) = c(1,i)
            vlist(i,2) = c(2,i)
            vlist(i,3) = c(3,i)
            vlist(i,4) = czan(i)
 30      continue
c
         ngrps = 0
         ospbas = .false.
         do 240 i = 1 , nshell
            l = katom(i)
            nt1 = kng(i)
            is = kstart(i)
            if = is + nt1 - 1
            mini = kmin(i)
            maxi = kmax(i)
            ngrps = ngrps + 1
            jj = maxi - mini
            if (jj.ne.0) then
               if (jj.lt.3) go to 120
               if (jj.eq.3) then
                  ospbas = .true.
               else
                  k = 0
                  m = 3
                  do 40 ig = is , if
                     k = k + 1
                     alp(k) = ex(ig)
                     ee = alp(k) + alp(k)
                     fac = pi32/(ee*dsqrt(ee))
                     fac = pt75*fac/(ee*ee)
                     c1(k) = cd(ig)/dsqrt(fac)
 40               continue
                  facd = dzero
                  do 60 ig = 1 , k
                     do 50 jg = 1 , ig
                        ee = alp(ig) + alp(jg)
                        fac = ee*dsqrt(ee)
                        dum = pt75*c1(ig)*c1(jg)/(ee**2*fac)
                        if (ig.ne.jg) then
                           dum = dum + dum
                        end if
                        facd = facd + dum
 50                  continue
 60               continue
                  do 70 ig = 1 , k
                     c1(ig) = c1(ig)/dsqrt(facd*pi32)
 70               continue
                  go to 170
               end if
            end if
            m = 1
            k = 0
            do 80 ig = is , if
               k = k + 1
               alp(k) = ex(ig)
               ee = alp(k) + alp(k)
               fac = pi32/(ee*dsqrt(ee))
               c1(k) = cs(ig)/dsqrt(fac)
 80         continue
            facs = dzero
            do 100 ig = 1 , k
               do 90 jg = 1 , ig
                  ee = alp(ig) + alp(jg)
                  fac = ee*dsqrt(ee)
                  dum = c1(ig)*c1(jg)/fac
                  if (ig.ne.jg) then
                     dum = dum + dum
                  end if
                  facs = facs + dum
 90            continue
 100        continue
            do 110 ig = 1 , k
               c1(ig) = c1(ig)/dsqrt(facs*pi32)
 110        continue
            go to 170
 120        k = 0
            m = 2
            do 130 ig = is , if
               k = k + 1
               alp(k) = ex(ig)
               ee = alp(k) + alp(k)
               fac = pi32/(ee*dsqrt(ee))
               fac = pt5*fac/ee
               c1(k) = cp(ig)/dsqrt(fac)
 130        continue
            facp = dzero
            do 150 ig = 1 , k
               do 140 jg = 1 , ig
                  ee = alp(ig) + alp(jg)
                  fac = ee*dsqrt(ee)
                  dum = pt5*c1(ig)*c1(jg)/(ee*fac)
                  if (ig.ne.jg) then
                     dum = dum + dum
                  end if
                  facp = facp + dum
 140           continue
 150        continue
            do 160 ig = 1 , k
               c1(ig) = c1(ig)/dsqrt(facp*pi32)
 160        continue
 170        nn = ngrp(m)
            do 230 j = 1 , nn
               ktw = ktw + 1
               kls = kls + 1
               nfirst(ktw) = kls
               kls = kls + nt1 - 1
               nlast(ktw) = kls
               nctr(ktw) = l
c
               go to (180,190,200) , m
 180           ntype(ktw) = j
               go to 210
 190           ntype(ktw) = j + 1
               go to 210
 200           ntype(ktw) = j + 4
               maxtyp = 10
 210           dum = 1.0d0
               if (j.ge.4) dum = sqrt3
               do 220 k = 1 , nt1
                  kk = kk + 1
                  isky(kk) = ntype(ktw)
                  wnorm(kk) = c1(k)*dum
                  eta(kk) = alp(k)
 220           continue
 230        continue
            if (ospbas) then
               ngrps = ngrps + 1
               ospbas = .false.
               go to 120
            end if
 240     continue
c
         if (num.ne.ktw) then
            call caserr(
     +       'invalid parameters detected in analysis preprocessor')
         else if (ngrps.le.0 .or. ngrps.gt.mxgrps) then
            call caserr(
     +       'invalid parameters detected in analysis preprocessor')
         else if (kk.le.0 .or. kk.gt.mxgaus) then
            call caserr(
     +       'invalid parameters detected in analysis preprocessor')
         else
            nbfnd = ktw
            len3 = nx
            newbas = num
            lenb = len3 + 3
            lenbad = lenb
            mgrps = ngrps
            lentrd = lensec(len3)
            if (newbas.ne.nbfnd) then
               call caserr(
     +       'invalid parameters detected in analysis preprocessor')
            else
c
c ----- now setup gamess/atmol mapping vectors
c
               ij = 0
               do 260 i = 1 , 3
                  m1 = minf(i)
                  m2 = maxf(i)
                  do 250 j = 1 , newbas
                     jj = ntype(j)
                     if (jj.ge.m1 .and. jj.le.m2) then
                        ij = ij + 1
                        mapsie(ij) = j
                     end if
 250              continue
 260           continue
c
               return
            end if
         end if
      end if
      end
**==propin.f
      subroutine propin(q)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *3 char3i
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
      dimension x(3),ngrp(3),q(*)
c
      integer mp, lwordp, lvec, ied3, newbas, lenb
      integer nbfnd, lenbad, lentrd, len3, ngrpm
      integer idena, nwroot, nosec
      logical noscf
      common /limy/ mp(10),lwordp,lvec,ied3,newbas,lenb,
     +              nbfnd,lenbad,lentrd,len3,ngrpm,
     +              idena,noscf,nwroot,nosec(mxroot)
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
      common/restri/jjfile(63),lds(508),isect(508),ldsect(508)
      common/miscop/clist(100,3),cenmax(100),
     * newpro,nocp,mcen,mspace,itemp(100),jtemp(100),ktemp(100),
     * ltemp(maxat),cenmas(maxat),eig(maxorb),eig1(mxorb1),jeig(6)
c
      real*8 cconv, atwt, cspin
      integer lin, kin, nin, ngmx, nfmx, nomx, ncmx, ntmx, ntcol, npmx
      integer icon, ncom, maxnt,itps, maxcon
c
      common /datprp/ cconv(11),atwt(30),cspin(30),
     + lin(84),kin(84),nin(84),
     + ngmx,nfmx,nomx,ncmx,ntmx,ntcol,npmx,
     + icon(24),ncom(22),maxnt,itps,maxcon
c
      common/junk/eta(mxgaus),wnorm(mxgaus),vlist(maxat,4),
     *frocc(maxorb),
     *nfirst(maxorb),nlast(maxorb),nctr(maxorb),ntype(maxorb)
     *,isky(mxgaus),nocco(maxorb),noc,non,maxtyp,ngaus,nbfns,newb,
     *ncol,ngrps,idddd,maxcom
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
      common/junkc/ztype(35),zlist(100),zcentg(100),ztagc(100),
     * znumb(100),zhead(5,22),zunit(33),ztitb(10,22),
     * zspin(33),zcomm(19),zzbuff(10)
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
      real*8 ctran
      integer ilifc, ntran, itran
      logical otran,otri
      common/tran/ilifc(maxorb),ntran(maxorb),itran(mxorb3),
     +            ctran(mxorb3),otran,otri
c
      common/reso/iresc(100),jresc(100),istat(100),iposp(100)
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
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
      logical oming, oextg, odirec, ovcfre, osemi, ostopm, olevd
      logical osym_op,osym_dens
      integer mina, minb, mouta, moutb, lock, ibrk
      integer numg, isecg, iblk3g, mextra, nsymo, iorbsy, isunor
      real*8  gapa1, gapa2, gapb1, gapb2, scaleg, symtag
      common/atmol3/mina,minb,mouta,moutb,lock,ibrk,
     + numg(2),isecg(2),iblk3g(2),mextra,nsymo,iorbsy(100),
     + oming, oextg,
     + gapa1,gapa2,gapb1,gapb2,scaleg,symtag(9),odirec(50),
     + isunor,ovcfre,osemi,ostopm,olevd,osym_op,osym_dens
c
c
      data zspace,x/' ','s','p','d'/
      data itvec,m350/3,350/
      data zcomp,zretr/'compute','retrieve'/
      data dzero,thrsh/0.0d0,1.0d-7/
      data mm1,m0,m1,m2/-1,0,1,2/
      data m3,ngrp/3,1,3,6/
c
      do 10 i = 1 , non
         ltemp(i) = m0
 10   continue
      do 20 i = 1 , npmx
         iposp(i) = m0
         ktemp(i) = m0
         istat(i) = mm1
 20   continue
      do 30 i = 1 , maxorb
         nocco(i) = m0
         frocc(i) = dzero
 30   continue
      nav = lenwrd()
      npmx2 = npmx*4/nav
      non1 = non + m1
      noc = non1
c
c     print flags
c
      if (oprint(52)) icon(8) = 1
      if (oprint(51)) icon(9) = 1
      if (oprint(53)) icon(7) = 1
c
c
c ----- retrieve remainder of =analysis= data from ed7
c
      nw = 400 + 204/nav
      call reads(clist,nw,num8)
      nw = 831
      call rdchrs(ztype,nw,num8)
c
c     check valid data
c
      if (mcen.ge.0 .and. mcen.le.ncmx) then
         if (newpro.ge.0 .and. newpro.le.npmx) then
            if (nocp.ge.0 .and. nocp.le.ncmx) then
c
c ...
c ... nucl
c ...
               if (mcen.ne.0) then

                  do 40 ii = 1 , mcen
                     i = locatc(zaname,non,zcentg(ii))
                     if (i.eq.0) call caserr(
     +               'nuclidic mass specified for unrecognised centre')
                     cenmas(i) = cenmax(ii)
                     ltemp(i) = 1
 40               continue
               end if
c ...
c ... centres
c ...
               if (nocp.ne.0) then
                  do 50 i = 1 , nocp
                     noc = noc + 1
                     if (noc.gt.maxat) call caserr(
     +       'maximum number of centres exceeded in property analysis')
                     vlist(noc,1) = clist(i,1)
                     vlist(noc,2) = clist(i,2)
                     vlist(noc,3) = clist(i,3)
                     vlist(noc,4) = dzero
                     zs = zlist(i)
                     if (zs.eq.zspace) zs = char3i(noc)
                     zaname(noc) = zs
 50               continue
               end if
c ...
c ... restore/vectors--mouta
c ...
               if (mouta.le.m0 .or. mouta.gt.m350) then
                  call caserr('retrieved eigenvectors are invalid')
               else
                  call secget(mouta,itvec,iblkv)
                  j = iblkv
                  call getqp(zcomm,zzbuff,eig,frocc,nbas,newb,ncol,
     +                       ieig,idiff,nomx,iblkv)
                  if (idiff.eq.1.d0) go to 70
                  call caserr(
     +                'no occupation nos. found in vectors section')
               end if
            end if
         end if
      end if
 60   call caserr('invalid parameters detected in property routines')
 70   do 80 i = 1 , ncol
         if (dabs(frocc(i)).gt.thrsh) nocco(i) = m1
 80   continue
      write (iwr,6010) mouta , ibl3d , yed(ied3) , (zcomm(7-i),i=1,6) ,
     +                 zzbuff
c ...
c ... find
c ...
      if (irest.eq.7) then
         write (iwr,6020)
         call secget(isect(468),itps,ibl468)
         call readi(iresc,npmx2*nav,ibl468,idaf)
         irest = 0
      end if
c
c ... enter
c ...

c
c determine how many isotope patterns have been input using
c the weights directive
c
      nmv = mass_numvec()
      if(nmv  .ne. 1)then
        write(iwr,6021)nmv
 6021   format(1x,/,'Properties will be computed for ',i4,
     x       ' different isotopic substitution patterns :'/)

         call mass_print

      endif

      do imv = 1, nmv

      write (iwr,6030)
      ss = dzero
      tw = dzero
      rat = dzero
      sss = dzero
      do 90 i = 1 , non
         if (ltemp(i).le.0) then
            pow1 = amass_get(imv,i)
         else
c
c explicit input from nucl directive, overrides
c weights
c
            pow1 = cenmas(i)
         end if
         ss = ss + vlist(i,1)*pow1
         tw = tw + vlist(i,2)*pow1
         sss = sss + vlist(i,3)*pow1
         write (iwr,6040) (vlist(i,j),j=1,4) , zaname(i) , pow1
         rat = rat + pow1
 90   continue
      write (iwr,6050) non
      write (iwr,6060)
      if (non1.gt.maxat) go to 60
      zaname(non1) = char3i(non1)
      vlist(non1,4) = dzero
      vlist(non1,1) = ss/rat
      vlist(non1,2) = tw/rat
      vlist(non1,3) = sss/rat
      write (iwr,6080) (vlist(non1,j),j=1,3) , zaname(non1)
      if (noc.ne.non1) then
         write (iwr,6070)
         non2 = non1 + m1
         do 100 i = non2 , noc
            write (iwr,6080) (vlist(i,j),j=1,3) , zaname(i)
 100     continue
      end if
      if (oprint(53)) write (iwr,6090)
      ngr = 0
      i = 1
 110  if (i.gt.nbfns) then
         write (iwr,6100) ngr , nbfns , ngaus
         if (ngr.ne.ngrps) go to 60
         if (idddd.eq.m1) write (iwr,6110) nbfnd
         if (newpro.gt.0) then
            do 170 i = 1 , newpro
               ip = itemp(i)
               jp = locatc(zaname,noc,ztagc(i))
               if (jp.eq.0) call caserr(
     +             'unrecognised centre specified in property directive'
     +             )
               kp = jtemp(i)
               if (kp.ne.m0) then
c ... is any 'old' property overwritten on dumpfile
                  do 120 j = 1 , npmx
                     jj = iposp(j)
                     if (jj.lt.0) then
                        jj = -jj
                     else if (jj.eq.0) then
                        go to 120
                     end if
                     if (jj.eq.kp) then
                        if (iresc(j).ne.ip .or. jresc(j).ne.jp) then
                           iposp(j) = m0
                        end if
                     end if
 120              continue
               end if
c ... is property already on list
               k = m1
               do 130 j = 1 , npmx
                  kk = iposp(j)
                  if (kk.eq.m0 .and. istat(j).eq.mm1) go to 130
                  if (iresc(j).ne.ip .or. jresc(j).ne.jp) go to 130
                  if (kk.lt.0) then
                     kk = -kk
                  else if (kk.eq.0) then
                     if (kp.ne.m0) then
                        iposp(j) = -kp
                     end if
                     go to 170
                  end if
                  if (kk.eq.kp) go to 160
                  k = j
                  go to 140
 130           continue
c ... insert property into list
 140           do 150 j = k , npmx
                  if (iposp(j).eq.m0 .and. istat(j).eq.mm1) then
                     iposp(j) = -kp
                     iresc(j) = ip
                     jresc(j) = jp
                     go to 160
                  end if
 150           continue
               call caserr('size of property list exceeded')
 160           istat(j) = m1
 170        continue
         end if
      call icopy(npmx2*nav,iresc,1,itemp,1)
         call secput(isect(468),itps,m1,ibl468)
         write (iwr,6120)
         ij = m0
         do 210 i = 1 , npmx
            is = ktemp(i)
            if (is.le.0) go to 210
            ip = ltemp(i)
            mop = itemp(i)
            ityp = mop + itps
            if (mop.le.m0 .or. mop.gt.maxnt)
     +          call caserr('invalid property index specified')
            j = jtemp(i)
            zss = zaname(j)
            lenp = lentrd*ncom(mop) + 2
            ij = ij + m1
            do 180 jdr = 1 , npmx
               if (mop.eq.itemp(jdr) .and. j.eq.jtemp(jdr) .and.
     +             ltemp(jdr).gt.m0) go to 190
 180        continue
            write (iwr,6130) (zhead(ki,mop),ki=1,5) , zss , zcomp
            go to 200
 190        call secget(ltemp(jdr),ityp,ibl468)
            write (iwr,6140) (zhead(ki,mop),ki=1,5) , zss , zretr ,
     +                       ltemp(jdr) , ibl468
 200        if (ip.lt.0) then
               ip = -ip
               call secput(ip,ityp,lenp,ibl468)
               write (iwr,6150) ip , ibl468
            else if (ip.eq.0) then
               write (iwr,6160)
            end if
            ltemp(i) = ip
            ktemp(i) = mm1
 210     continue
         call secini(ibl3d,idaf)
         if (ij.le.0) then
            write (iwr,6170)
            return
         else
            ns1 = m0
            do 220 j = 1 , npmx
               if (istat(j).ge.m1) then
                  kk = ncom(iresc(j))
                  if (kk.gt.ns1) ns1 = kk
               end if
 220        continue
c ... ns1 is max.no. of components(i.e. triangles) requested
            maxcom = ns1
            int4 = ns1*lenb
c           nbsq = ncol*nbas   last vector is there but zeroed
            nbsq = nbas*nbas
c
            i00 = 1
            i01 = i00 + mxgaus/nav
            i02 = i01 + mxgaus/nav
            i10 = i02 + mxgaus/nav
            i20 = i10 + nbsq
            i30 = i20 + lenb
            need = i30 + int4
            if (oprint(53)) write (iwr,6180) lwordp , need
            if (lwordp.lt.need) then
c ... insufficient core to handle components
               if (lwordp.lt.i30) then
                  call caserr(
     +          'insufficient memory available for requested analysis')
               else
                  iparam = lwordp - i30
                  maxcom = iparam/lenb
                  if (maxcom.le.m0) then
                     call caserr(
     +          'insufficient memory available for requested analysis')
                  else
c ... maxcom is max. no. of triangles that may be handled
                     int4 = (ns1-maxcom)*lenb
                     write (iwr,6190) maxcom , int4
                     go to 240
                  end if
               end if
            else if (lwordp.eq.need) then
               go to 240
            else
               int4 = (lwordp-need)
               if (int4.ne.0) then
                  write (iwr,6200) int4
               end if
               go to 240
            end if
         end if
      else
         n = ntype(i)
         l = nctr(i)
         ia = nfirst(i)
         ib = nlast(i)
         nn = ib - ia + 1
         if (n.lt.1) go to 60
         if (n.eq.1) then
            jj = ngrp(1)
            jjj = m1
         else if (n.le.4) then
            jj = ngrp(2)
            jjj = m2
         else
            if (n.gt.10) go to 60
            jj = ngrp(3)
            jjj = m3
         end if
      end if
      ngr = ngr + 1
      if (oprint(53)) write (iwr,6210) ngr, x(jjj), zaname(l), nn
      if (l.lt.m1 .or. l.gt.non) go to 60
      do 230 j = ia , ib
         dum = eta(j)
         iski = isky(j)
         c = wnorm(j)*dsqrt(ggntde(j,iski,l,j,iski,l))
         if (oprint(53)) write (iwr,6220) c , dum
         if (dum.lt.1.0d-8) go to 60
         if (dabs(c).lt.1.0d-8) go to 60
 230  continue
      i = i + jj
      go to 110
 240  mp(1) = 1
      do 250 i = 1 , 9
         mp(i+1) = mp(i) + lenb
 250  continue
      call rdedx(q(i10),nbsq,iblkv,idaf)
      call tdown(q(i20),ilifq,q(i10),ilifq,ncol)
      if (oprint(53)) then
         if (.not.otran) then
            write (iwr,6230)
            do 270 i = 1 , newb
               n = ntran(i)
               l = ilifc(i) + m1
               write (iwr,6240) ctran(l) , itran(l) , n , i
               if (n.ne.m1) then
                  do 260 j = 2 , n
                     l = l + m1
                     write (iwr,6240) ctran(l) , itran(l)
 260              continue
               end if
 270        continue
         end if
         write (iwr,6250) zcomm(5)
         call prsql(q(i20),ncol,nbas,nbas)
         write (iwr,6260) mouta
         write (iwr,6270) (frocc(i),i=1,ncol)
      end if
      call secput(isect(468),itps,m1,ibl468)
      call wrt3i(iresc,npmx2*nav,ibl468,idaf)
      call revind
c
      call propdr(q(i00),q(i01),q(i02),q(i10),q(i20),q(i30))
c

c
c end of loop over mass substitutions
c
      enddo

      return
 6010 format (1x,'vectors restored from section',i4,
     +        ' of dumpfile starting at block',i6,' of ',
     +        a4//' header block information : '/
     +        ' vectors created under acct. ',a8/1x,a7,
     +        'vectors created by ',a8,' program at ',a8,' on ',a8,
     +        ' in the job ',a8/' with the title : ',10a8/)
 6020 format (/20x,'restart computation of properties'//)
 6030 format (/1x,104('-')//40x,18('*')/40x,'molecular geometry'/40x,
     +        18('*')//8x,'x',15x,'y',15x,'z',13x,'charge',5x,'label',
     +        4x,'nuclidic mass'/)
 6040 format (4f16.7,1x,a8,f14.7)
 6050 format (' number of nuclei  ',i3)
 6060 format (/40x,14('*')/40x,'centre of mass'/40x,14('*')/)
 6070 format (/40x,18('*')/40x,'additional centres'/40x,18('*')//8x,'x',
     +        15x,'y',15x,'z',24x,'label'/)
 6080 format (3f16.7,17x,a8)
 6090 format (/1x,104('-')//40x,19('*')/40x,'molecular basis set'/40x,
     +        19('*')///8x,'ordered list of groups (normalised)'/)
 6100 format (/' total number of groups ',4x,
     +        i4/' total number of gtos ',6x,
     +        i4/' total number of primitives ',i4/)
 6110 format (/' no. of spherical harmonic gtos=',i4)
 6120 format (/1x,104('-')//40x,25('*')/40x,
     +        'list of active properties'/40x,25('*')//10x,'property',
     +        14x,'centre',5x,'integrals'/1x,51('=')/)
 6130 format (/1x,5a6,1x,a8,3x,a8)
 6140 format (/1x,5a6,1x,a8,3x,a8,5x,'from dumpfile :section ',i3,
     +        ' (at block ',i6,')')
 6150 format (/56x,'route to dumpfile :section ',i3,' (at block ',i6,
     +        ')')
 6160 format (/56x,'no routing of integrals specified')
 6170 format (/40x,'**** no active properties-terminate analysis')
 6180 format (/1x,104('-')//20x,'core store analysis'/20x,19('-')//20x,
     +        'main core available=',i7,' words'//20x,
     +        'main core required =',i7,' words'/)
 6190 format (/15x,'**** insufficient core has been assigned to compute'
     +        /'     all properties'//15x,'**** only those with ',i2,
     +        ' or less components will be'/'     evaluated'//15x,
     +        '**** submit  startup job,modifying properties'/
     +        '     directive) assigning'/15x,'     additional ',i7,
     +        'words of core'//)
 6200 format (/15x,'**** step has requested ',i7,
     +        ' words than it will use'/)
 6210 format (/' group type centre  nterm'/1x,24('-')/1x,i5,1x,a1,4x,a8,
     +        i5//31x,'ctran',12x,'zeta'/31x,21('-')/)
 6220 format (20x,2f16.7)
 6230 format (//1x,104('-')//40x,12('*')/40x,'list of sabf'/40x,12('*')
     +        //11x,'transformed orbitals (normalised)'//
     +        ' coefficient old orbital nterm new orbital'//)
 6240 format (1x,f11.7,i12,i6,i12)
 6250 format (//1x,104('-')//40x,16('*')/40x,a8,' vectors'/40x,16('*')
     +        //)
 6260 format (//40x,'orbital populations  (from vectors section ',i3,
     +        ')'//)
 6270 format (/15x,7f15.7)
      end
**==putpro.f
      subroutine putpro(dint,zhead,ztit,ztagg,pcen,qcen,rcen,
     *chgcen,zspi,cspi,mop,nt,iblk)
c...
c...    write property integrals to dumpfile
c...
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
      dimension dint(*),zhead(5,*),ztit(10,*)
c
      integer mp, lwordp, lvec, ied3, newbas, lenb
      integer nbfnd, lenbad, lentrd, len3, ngrpm
      integer idena, nwroot, nosec
      logical noscf
      common /limy/ mp(10),lwordp,lvec,ied3,newbas,lenb,
     +              nbfnd,lenbad,lentrd,len3,ngrpm,
     +              idena,noscf,nwroot,nosec(mxroot)
c
      integer numdu, iblkdu, maxb, iblkla
      logical orevis
      real*8 apos
      common/sector/numdu,iblkdu,orevis(2),
     + apos(508),maxb,iblkla
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
      common/junkc/zheadp(8),ztitp(10),ztagp,zspip
      common/miscop/pnuc(10),cx,cy,cz,chgp,cspip,
     *mopp,ncom,lenp,lentrp
c
      do 20 i = 1 , 5
         zheadp(i+3) = zhead(i,mop)
 20   continue
      do 30 i = 1 , 10
         ztitp(i) = ztit(i,mop)
 30   continue
      zheadp(1) = zcom(1)
      zheadp(2) = zcom(3)
      zheadp(3) = zcom(2)
      mword1 = 20
      mword2 = 15 + 4/lenwrd()
      cx = pcen
      cy = qcen
      cz = rcen
      chgp = chgcen
      ztagp = ztagg
      zspip = zspi
      cspip = cspi
      mopp = mop
      ncom = nt
      lenp = len3
      lentrp = lentrd
      do 40 i = 1 , nt
         pnuc(i) = dint(mp(i))
 40   continue
c
      call wrtc(zheadp,mword1,iblk,numdu)
      call wrt3s(pnuc,mword2,numdu)
      iblk = iblk + 2
      do 50 i = 1 , nt
         itest = mp(i) + 3
         call wrt3(dint(itest),len3,iblk,numdu)
         iblk = iblk + lentrd
 50   continue
c
      return
      end
**==qudnuc.f
      subroutine qudnuc(noc,vlist,c,dint,nocmx)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension dint(7),vlist(nocmx,4),c(3)
      common/miscop/qnuc(6)
c...
c...compute nuclear contributions to quadrupole moments.
c...qm(xx)=(3*xx-rr)/2, similarly for qm(yy) and qm(zz).
c...qm(xy)=1.5*xy, similarly for qm(xz) and qm(yz).
c...also  rr=xx+yy+zz.
c...
      do 20 i = 1 , 6
         qnuc(i) = 0.0d0
 20   continue
c...
      do 50 i = 1 , noc
         mu = 0
         do 40 j = 1 , 3
            do 30 k = 1 , j
               mu = mu + 1
               qnuc(mu) = qnuc(mu) + vlist(i,4)*(vlist(i,j)-c(j))
     +                    *(vlist(i,k)-c(k))
 30         continue
 40      continue
 50   continue
c...
      dint(1) = qnuc(1)
      dint(2) = qnuc(3)
      dint(3) = qnuc(6)
      dint(4) = 3.0d0*qnuc(2)/2.0d0
      dint(5) = 3.0d0*qnuc(4)/2.0d0
      dint(6) = 3.0d0*qnuc(5)/2.0d0
      dint(7) = dint(1) + dint(2) + dint(3)
      dint(1) = (3.0d0*dint(1)-dint(7))/2.0d0
      dint(2) = (3.0d0*dint(2)-dint(7))/2.0d0
      dint(3) = (3.0d0*dint(3)-dint(7))/2.0d0
c...
      return
      end
**==rcasec.f
      subroutine rcasec(isec,icalc,iunit,nsp,cq,odata,ierr)
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
c  retrieve a grid data definition from dumpfile, put data in cq
c  the grid will be retrieved, with data into gq if not in memory
c
c
c grid definition parameters
c
      common/dfgrid/geom(12,mxgrid),gdata(5,mxgrid),igtype(mxgrid),
     &             igsect(mxgrid),nptot(mxgrid),npt(3,mxgrid),
     &             igdata(5,mxgrid),igstat(mxgrid),ngdata(mxgrid),
     &             ngrid
c
c data calculation parameters
c
      common/dfcalc/cdata(5,mxcalc),ictype(mxcalc),icsect(mxcalc),
     &            icgrid(mxcalc),icgsec(mxcalc),ndata(mxcalc),
     &            icstat(mxcalc),icdata(5,mxcalc),ncalc
      common/dfclc2/cscal(mxcalc),iblkp(mxcalc)
c
c labels and titles
c
      common/cplot/zgridt(10,mxgrid),zgrid(mxgrid),
     &             zcalct(10,mxcalc),zcalc(mxcalc),
     &             zplott(10,mxplot),zplot(mxplot)
      dimension cq(*)
      dimension geom1(12), igdat1(5), gdat1(5), npt1(3)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common/psh/crap(100)
      integer numdu, iblkdu, maxb, iblkla
      logical orevis
      real*8 apos
      common/sector/numdu,iblkdu,orevis(2),
     + apos(508),maxb,iblkla
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
      data itcalc/51/
c
      nav = lenwrd()
c
      numd1 = numdu
      iblkd1 = iblkdu

      ofor = (iunit.ne.idaf.or.nsp.ne.ibl3d)
      oreset = (iunit.ne.numdu.or.nsp.ne.iblkdu)

      if(oreset)then
         call filec(yed(iunit),nsp,iunit,irep)
         if(irep.ne.0)then
            if(ierr.ne.0)then
               ierr = 3
               return
            else
               call caserr('file not allocated')
            endif
         endif
         call secini(nsp,iunit)
      endif

      if(.not.ofor)then
         icalc = 0
         do 10 i = 1,ncalc
            if(icsect(i).eq.isec)icalc = i
 10      continue
      endif

      if(icalc.eq.0)then
         ncalc = ncalc + 1
         icalc = ncalc
      endif
c
      call secget(isec,itcalc,iblk)
      write (iwr,1000) isec , iblk, yed(numdu)
      if(iblk.eq.0)then
         if(ierr.ne.0)then
            ierr = 1
            goto 999
         else
            call  caserr
     &           ('attempt to retrieve undefined grid data section')
         endif
      endif
      nw = 23 + 1 + (20-1)/nav
      call rdedx(crap,nw,iblk,numdu)
      call stuf1
      call unstfi(crap,ictype(icalc),1)
      call unstfi(crap,ndata(icalc),1)
      call unstfi(crap,icdata(1,icalc),5)
      call unstff(crap,cdata(1,icalc),5)
      call unstff(crap,cscal(icalc),1)
      call unstff(crap,geom1,12)
      call unstfi(crap,igtyp1,1)
      call unstfi(crap,nptot1,1)
      call unstfi(crap,npt1  ,3)
      call unstfi(crap,igdat1,5)
      call unstff(crap,gdat1,5)
      call unstfi(crap,igstt1,1)
      call unstfi(crap,igsct1,1)
      call unstfi(crap,ngdat1,1)
      nw1 = nstuf()
      if(nw.ne.nw1)then
         write(iwr,*)nw,nw1
         call caserr('nw error in rcasec')
      endif
c  do we have the grid stored
      ig = 0
      if(.not.ofor)then
         do 20 i = 1,ngrid
            if ((igsect(i).ne.0).and.(igsect(i).eq.igsct1)) ig = i
 20      continue
      endif
      if(ig.eq.0)then
         ngrid = ngrid + 1
         ig = ngrid
         do 25 i=1,12
            geom(i,ig) = geom1(i)
 25      continue
         do 26 i = 1,3
            npt(i,ig) = npt1(i)
 26      continue
         do 27 i=1,5
            gdata(i,ig)=gdat1(i)
            igdata(i,ig)=igdat1(i)
 27      continue
         igtype(ig) = igtyp1
         nptot(ig) = nptot1
         if(.not.ofor)then
            igsect(ig) = igsct1
         else
            igsect(ig) = 0
         endif
         ngdata(ig) = ngdat1
         igstat(ig) = 1
      endif
      icgrid(icalc) = ig
c @ need to account for restart with icstat
      icstat(icalc) = 1
      nwc = ndata(icalc)
c
      call rdchrs(zcalct(1,icalc),10,numdu)
      if(odata)then
         call reads(cq,nwc,numdu)
         icstat(icalc) = 2
      endif
      if(.not.ofor)icsect(icalc) = isec
      ierr=0
 999  if(oreset)call secini(iblkd1,numd1)
      return
 1000 format (//' restore grid values from section',i4,
     +        ' of dumpfile starting at block',i6,' of ',a4/)
      end
**==rgrsec.f
      subroutine rgrsec(isec,igrid,iunit,nsp,gq,odata,ierr)
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
c  retrieve a grid definition from dumpfile, put irreg data in gq
c
c
c grid definition parameters
c
      common/dfgrid/geom(12,mxgrid),gdata(5,mxgrid),igtype(mxgrid),
     &             igsect(mxgrid),nptot(mxgrid),npt(3,mxgrid),
     &             igdata(5,mxgrid),igstat(mxgrid),ngdata(mxgrid),
     &             ngrid
c
c labels and titles
c
      common/cplot/zgridt(10,mxgrid),zgrid(mxgrid),
     &             zcalct(10,mxcalc),zcalc(mxcalc),
     &             zplott(10,mxplot),zplot(mxplot)
      integer numdu, iblkdu, maxb, iblkla
      logical orevis
      real*8 apos
      common/sector/numdu,iblkdu,orevis(2),
     + apos(508),maxb,iblkla
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common/psh/crap(100)
      data itgrid/50/
c
      nav = lenwrd()
c
      numd1 = numdu
      iblkd1 = iblkdu
c
      ofor = (iunit.ne.idaf.or.nsp.ne.ibl3d)
      oreset = (iunit.ne.numdu.or.nsp.ne.iblkdu)
c
c
      if(iunit.ne.numdu)then
         call filec(yed(iunit),nsp,iunit,irep)
         if(irep.ne.0)then
            if(ierr.ne.0)then
               ierr = 3
               return
            else
               call caserr('file not allocated')
            endif
         endif
         call secini(nsp,iunit)
      endif
c
      if(.not.ofor)then
         igrid = 0
         do 10 i = 1,ngrid
            if(igsect(i).eq.isec)igrid = i
 10      continue
      endif
c
      if(igrid.eq.0)then
         ngrid = ngrid + 1
         igrid = ngrid
      endif
c
c
c read from a grid dumpfile section - irregular data in gq
c
      call secget(isec,itgrid,iblk)
      write (iwr,1000) isec , iblk, yed(numdu)
c
      if(iblk.eq.0)then
         if(ierr.ne.0)then
            ierr = 1
            goto 999
         else
            call  caserr
     &           ('attempt to retrieve undefined grid section')
         endif
      endif
c
      nw = 17  + 1 + (12-1)/nav
c
      call rdedx(crap,nw,iblk,numdu)
c
      call stuf1
      call unstff(crap,geom(1,igrid),12)
      call unstfi(crap,igtype(igrid),1)
      call unstfi(crap,nptot(igrid),1)
      call unstfi(crap,npt(1,igrid),3)
      call unstfi(crap,igdata(1,igrid),5)
      call unstff(crap,gdata(1,igrid),5)
      call unstfi(crap,igstat(igrid),1)
      call unstfi(crap,ngdata(igrid),1)
      nw1 = nstuf()
c
      if(nw.ne.nw1)then
         write(iwr,*)nw,nw1
         call caserr('nw error in rgrsec')
      endif
c
      igstat(igrid)=1
      nwg = ngdata(igrid)
      call rdchrs(zgridt(1,igrid),10,numdu)
      if(odata)then
         if(nwg.ne.0)call reads(gq,nwg,numdu)
         igstat(igrid) = 2
      endif
c
c save the section
c
      if(.not.ofor)then
         igsect(igrid)=isec
      endif
c
      ierr=0
 999  if(oreset)call secini(iblkd1,numd1)
      return
 1000 format (//' restore grid from section',i4,
     +        ' of dumpfile starting at block',i6,' of ',a4/)
      end
**==setgrd.f
      subroutine setgrd(i,iwr)
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
      common/dfgrid/geom(12,mxgrid),gdata(5,mxgrid),igtype(mxgrid),
     &             igsect(mxgrid),nptot(mxgrid),npt(3,mxgrid),
     &             igdata(5,mxgrid),igstat(mxgrid),ngdata(mxgrid),
     &             ngrid
c
c labels and titles
c
      common/cplot/zgridt(10,mxgrid),zgrid(mxgrid),
     &             zcalct(10,mxcalc),zcalc(mxcalc),
     &             zplott(10,mxplot),zplot(mxplot)
      common/sgrd/ stepx(3),stepy(3),stepz(3),igrd
c
c
      do 10 j = 1,3
         if(igtype(i).eq.1.or.igtype(i).eq.2)then
            stepx(j)= (geom(3+j,i) - geom(j,i))/(npt(1,i)-1)
            stepy(j)= (geom(6+j,i) - geom(j,i))/(npt(2,i)-1)
         endif
         if(igtype(i).eq.2)
     &        stepz(j)= (geom(9+j,i) - geom(j,i))/(npt(3,i)-1)
 10   continue
      igrd = i
      write(iwr,100)i,(zgridt(j,i),j=1,6)
 100  format(1x,'current grid is number ',i3,' :',10a8)
      return
      end
**==setstp.f
      subroutine setstp(i)
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
c setstp - where to start the route
c
c the job sequence
      integer stype(mxstp), arg(mxstp)
      common/route/stype,arg,nstep,istep

      istep = i
      return
      end
**==shftce.f
      subroutine shftce(c,a,e,isky)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension c(4),a(3),e(3),cfs(10)
      common/miscop/aaa(70),
     *c1,c2,c3,c4,c5,c6,c7,c8,c9,c10
      equivalence(c1,cfs(1))
      data m10/10/
c
      call dcopy(m10,c(1),1,cfs(1),1)
      aex = e(1) - a(1)
      aey = e(2) - a(2)
      aez = e(3) - a(3)
      c(1) = c1 + aex*c2 + aey*c3 + aez*c4
      if (isky.gt.4) then
         c(1) = c(1) + aex*aex*c5 + aey*aey*c6 + aez*aez*c7 +
     +          aex*aey*c8 + aex*aez*c9 + aey*aez*c10
         c(2) = c2 + 2.0d0*aex*c5 + aey*c8 + aez*c9
         c(3) = c3 + 2.0d0*aey*c6 + aez*c10 + aex*c8
         c(4) = c4 + 2.0d0*aez*c7 + aex*c9 + aey*c10
      end if
c
      return
      end
**==siny.f
      function siny(ny,nab,ncd,pab,pcd)
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
      common/junk/pqn(99),expa(99),ctran(mxprms,30),ajmn(12),czn,
     *expmn,c(15,15),itran(mxprms,30),nrow(30),nbcon(4),nbas(4),
     *nsym,nosh(4),ncsh(4),nccup(4),ibase(33),
     *popul(4),facto(30),fatt(30),thre(15),a(15,15),eivr(15,15),
     *cm2scf(15,15),cm1scf(15,15),eps(15),gmin(15,15),smin(15,15),
     *populs(4),term,energ,pot,cin,vir,thrscf,enorm,
     *red,expdif,encor,
     *nsize(10),nsavf(4),n1(6),nbvar(30),nbas1,nbast,nbvar1,
     *nblock,numvar,nrun,nitscf,nitsc1,nitsc2,nconv,nalarm,no,nstep,
     *iiii,nsht,nosht,nflag3,nbasm,mvar,iex,newbas,nit,nbpos(8)
      data dzero,done/0.0d0,1.0d0/
c
c    calculates the integral i(ny) needed for the supermatrices p and q
c    when a slater   basis set is used
c    *******
      tabcd = pab/pcd
      nabny = nab - ny - 1
      ncdny = ncd + ny
      vnab = fatt(nabny+1)/pab**(nabny+1)
      vncd = fatt(ncdny+1)/pcd**(ncdny+1)
      term1 = dzero
      nalf = nabny + 1
      do 20 i = 1 , nalf
         term1 = term1 + fatt(nabny+ncdny+2)*tabcd**(i-1)
     +           /(fatt(i)*fatt(nabny+ncdny+3-i))
 20   continue
      cabcd = term1*(tabcd+done)**(-nabny-ncdny-1)
      siny = vnab*vncd*cabcd
c
      return
      end
**==slap.f
       function slap (l,m,n,gama)
c...    originally function olap
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension dftr(7)
      data pi/3.14159265d0/
      data dftr/1.0d0,1.0d0,3.0d0,15.0d0,105.0d0,945.0d0,10395.0d0/
      lh = l/2
      mh = m/2
      nh = n/2
      if (2*lh.eq.l) then
         if (2*mh.eq.m) then
            if (2*nh.eq.n) then
               slap = (dsqrt(pi/gama))**3*(.5d0/gama)**(lh+mh+nh)
     +                *dftr(lh+1)*dftr(mh+1)*dftr(nh+1)
               return
            end if
         end if
      end if
      slap = 0.0d0
c
      return
      end
**==sldnuc.f
      subroutine sldnuc(noc,vlist,c,dint,nocmx)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension dint(10),vlist(nocmx,4),c(3)
      common/miscop/pnuc(6)
c
      do 20 i = 1 , 6
         pnuc(i) = 0.0d0
 20   continue
      do 50 i = 1 , noc
         dst = (vlist(i,1)-c(1))**2 + (vlist(i,2)-c(2))
     +         **2 + (vlist(i,3)-c(3))**2
         if (dst.ne.0.0d0) then
            mu = 0
            do 40 j = 1 , 3
               do 30 k = 1 , j
                  mu = mu + 1
                  pnuc(mu) = pnuc(mu) + vlist(i,4)*(vlist(i,j)-c(j))
     +                       *(vlist(i,k)-c(k))/(dsqrt(dst))**3
 30            continue
 40         continue
         end if
 50   continue
c
      dint(1) = 1.5d0*(pnuc(3)+pnuc(6))
      dint(2) = 1.5d0*(pnuc(1)+pnuc(6))
      dint(3) = 1.5d0*(pnuc(1)+pnuc(3))
      dint(4) = -1.5d0*pnuc(2)
      dint(5) = -1.5d0*pnuc(4)
      dint(6) = -1.5d0*pnuc(5)
      dint(7) = pnuc(1) + pnuc(3) + pnuc(6)
      dint(8) = pnuc(1)
      dint(9) = pnuc(3)
      dint(10) = pnuc(6)
c
      return
      end
**==spgrid.f
      subroutine spgrid(grid, radius, ngrid, orig, index, 
     &                  itype, idat2, nmax, iwr)
c
c DESCRIPTION : This routine generates a spherical grid around a
c               molecule.  Grid points are calculated using spherical
c               coordinates which are then converted into x, y and z
c               coordinates.
c
c PARAMETERS  : NAME               DESCRIPTION
c               ----               -----------
c
c               grid        grid array
c               radius      Radius of sphere in Bohr's .
c               ngrid       Number of grid points on equator and also
c                           number of lattitudes.
c               orig        sphere's origin.
c               index       Number of grid points generated.
c
c AUTHOR      : M. A. FOX.
c
c Date        : 2/2/90
c
      implicit real*8 (a-h,o-z)
      integer point, bottom
      logical stag, ax
      dimension grid(*), orig(*)
      parameter (pi = 3.141592654d0, twopi = 6.283185307d0)
      data aaa,ccc,rand/54891.d0,65536.d0,290351.d0/
c 
      write(iwr,1000)(orig(i),i=1,3),radius,ngrid
      if(itype.eq.1)then
         write(iwr,1010)
      else if (itype.eq.2)then
         write(iwr,1020)idat2
      else
         call caserr('bad type integer in spgrid')
      endif
c
      index = 0
      ii = 0
      if(itype.eq.1)then
c
c calculate the longitudinal incremental angle.
c
         vrtinc = pi / (ngrid - 1)
         bottom = ngrid - 1
         do 10 lattde = 0, bottom
            phi = vrtinc * lattde
c
c calculate number of points on current lattitude
c
            numpts = idint((ngrid - 1) * dabs(dsin(phi))) + 1
            step = twopi / numpts
c
c generate initial random starting point to ensure no relation between
c initial starting points on different lattitudes
c
            randum=(aaa/ccc)*(rand/ccc)-idint((aaa/ccc)*(rand/ccc))
            rand=ccc*(ccc*randum)
            start = randum * twopi
            if(numpts*bottom.gt.nmax)
     &       call caserr('insufficient memory for sphere grid')
            do 20 point = 1, numpts
               index = index + 1
               theta = point * step + start
              grid(ii + 1) = radius * dsin(phi) * dcos(theta) + orig(1)
              grid(ii + 2) = radius * dsin(phi) * dsin(theta) + orig(2)
               grid(ii + 3) = radius * dcos(phi) + orig(3)
               ii = ii+3
 20         continue
 10      continue
      else if(itype.eq.2)then
c
c  axial symmetry preserving spherical grid
c
c axial symmetry number
         nsym = idat2
c
c calculate the longitudinal incremental angle.
c
         vrtinc = pi / (ngrid - 1)
         bottom = ngrid - 1
         stag=.false.
         ngrid = ngrid / nsym
         if (ngrid.eq.0)ngrid = 1
         asym = twopi/dble(nsym)
         do 110 lattde = 0, bottom
            ax=(lattde.eq.0.or.lattde.eq.bottom)
            phi = vrtinc * lattde
c
c calculate number of points on segment of current lattitude
c
         numpts = idint(0.5d0 + dble(ngrid - 1) * dabs(dsin(phi))) +1
            write(6,*)'layer',lattde,'numpts',numpts
            step = twopi / (nsym*numpts)
            start = pi*0.5d0
            if (stag) start = start + asym/2.0d0
            nloop=numpts
            if(.not.ax)nloop=numpts*nsym
            do 120 point = 1, nloop
               index = index + 1
               theta = point * step + start
              grid(ii + 1) = radius * dsin(phi) * dcos(theta) + orig(1)
              grid(ii + 2) = radius * dsin(phi) * dsin(theta) + orig(2)
               grid(ii + 3) = radius * dcos(phi) + orig(3)
               ii = ii+3
 120         continue
 110      continue
      endif
      write(iwr,1030)index
 1000 format(1x,/,1x,'generating spherical grid with origin (',
     &     f12.5,',',f12.5,',',f12.5,')',/,1x,
     &     'radius ',f12.5,/,1x,
     &     'number of longitude and lattitude increments',i4)
 1010 format(1x,'generate grid with random starting lattitude values')
 1020 format(1x,'generate grid to preserve axial symmetry ',i3)
 1030 format(1x,'total number of points generated ',i6,/)
      end
**==stuf1.f
      subroutine stuf1
      common/stuf/nwrd
      nwrd=0
      return
      end
**==stuff.f
      subroutine stuff(iq,ia,n)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension iq(*),ia(*)
      common/stuf/nwrd
      nav = lenwrd()
      do 10 i = 1,n*nav
 10      iq(nwrd+i)=ia(i)
      nwrd = nwrd + n*nav
      return
      end
**==stufi.f
      subroutine stufi(iq,ia,n)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension iq(*),ia(*)
      common/stuf/nwrd
      do 10 i = 1,n
 10      iq(nwrd+i)=ia(i)
      nwrd = nwrd + n
      return
      end
**==tmpl00.f
      subroutine tmpl00(core,lmem)
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
      dimension core(*)
c
c  store some large temporaries in dynamic memory 
c  (NOT in a common block (mxprim.gt.maxorb)
c
c     common/scrtch/tmao(mxprim*(mxprim+1)/2),tmp1(mxprim*maxorb),
c    + tmp2(mxprim*mxprim),vectp(mxprim*maxorb)
c
c
      integer mp, lwordp, lvec, ied3, newbas, lenb
      integer nbfnd, lenbad, lentrd, len3, ngrpm
      integer idena, nwroot, nosec
      logical noscf
      common /limy/ mp(10),lwordp,lvec,ied3,newbas,lenb,
     +              nbfnd,lenbad,lentrd,len3,ngrpm,
     +              idena,noscf,nwroot,nosec(mxroot)
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
      l1 = mxprim*(mxprim+1)/2
      l2 = mxprim*maxorb
      l3 = mxprim*mxprim
      l4 = maxorb*maxorb
c .. tmao
      i10 = 1
c .. tmp1
      i20 = i10 + l1
c .. tmp2
      i30 = i20 + l2
c .. vectp
      i40 = i30 + l3
c
c  .. vect and tmat (originally in /tmdata/
c
      i50 = i40 + l2
      i60 = i50 + l4
      last = i60 + l4
c
      ltot = lmem-1+last
       if(ltot.gt.lwordp) then
       write(iwr,10)ltot,lwordp
 10    format(/
     + 1x,'error in memory allocation of graphical routines'/
     + 1x,'not enough memory to continue'/
     + 1x,'TM plotting requires  ',1x,i9,' real words'/
     + 1x,'but has been allocated',1x,i9,' real words'/)
        call caserr('insufficient memory')
      endif
      call tmpl0(core(i10),core(i20),core(i30),core(i40),
     +           core(i50),core(i60),
     + l1,l2,l3,l2,l4,l4)
c
      return
      end
      function maketm(naos)
      implicit real*8 (a-h,p-z), integer  (i-n), logical (o)
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
      parameter (maxtrm = 10)
      common/junk2/idum(6),
     +wc(maxat),wx(maxat),wy(maxat),wz(maxat),
     +new(maxorb),icen(maxorb),ii(maxorb),ij(maxorb),ik(maxorb),
     +ic(maxorb),ccoef(maxorb,maxtrm),zeta(maxorb,maxtrm)
      common/blkin/x0(3,mxprim),zet(mxprim),inew(3,mxprim),
     +ipoint(maxorb),jpoint(mxprim),jcount(mxprim)
      dimension xx(3),in(3)
c
c  do for all atomic orbitals
c
      ntotal=0
      do 10, iii=1,naos
      i=new(iii)
      in(1)=ii(i)
      in(2)=ij(i)
      in(3)=ik(i)
      iatom=icen(i)
      xx(1)=wx(iatom)
      xx(2)=wy(iatom)
      xx(3)=wz(iatom)
      ncon=ic(i)
      ipoint(iii)=ntotal+1
      do 20,j=1,ncon
      ntotal=ntotal+1
      if(ntotal.gt.mxprim) 
     +  call caserr('maketm: too many primitive gaussians')
      jpoint(ntotal)=iii
      jcount(ntotal)=j
      zet(ntotal)=zeta(i,j)
      do 30,k=1,3
      x0(k,ntotal)=xx(k)
      inew(k,ntotal)=in(k)
30    continue
20    continue
c
10    continue
c
      maketm=ntotal
      return
      end
      function tmaov(iao,x,y,z,onorm)
      implicit real*8 (a-h,p-z), integer  (i-n), logical (o)
      real*8 koef
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
      parameter (maxtrm = 10)
      common/junk2/idum(6),
     +wc(maxat),wx(maxat),wy(maxat),wz(maxat),
     +new(maxorb),icen(maxorb),ii(maxorb),ij(maxorb),ik(maxorb),
     +ic(maxorb),ccoef(maxorb,maxtrm),zeta(maxorb,maxtrm)
c
c  calculates value of an contracted gaussian no. iao at position x,y,z
c  if ccoef refer to unnormalised primitives and make the 
c  contracted gaussian c  normalised (usual case) onorm 
c  should be .false. and the returned value will
c  be of normalised contracted AO
c
c     (2i-1)!!, i=0,5
      dimension dfact(0:5)
      data dfact /1.0d0,1.0d0,3.0d0,15.0d0,105.0d0,945.0d0/
c     (2/pi)**1.5 =
      data koef /0.507949087473927818336960627d0/
c
c
c
c  in tm.f the output of basis indirected through the new permutation
c  corresponds to the printout of scf vectors after scf run and
c  this corresponds to the lcao coefficients as they are saved in the
c  mrdci tran. moment module
c
c  onorm=.true. is necessary, if contract. coefficients refer to 
c  normalised  primitives, false for unnormalised primitives. 
c  The normalisation of whole contracted gaussian and the 
c  recalculation of contraction coefficients is a question of an 
c  another routine.
c
      n=new(iao)
c
      ncontr=ic(n)
      i=icen(n)
      dx= x-wx(i)
      dy= y-wy(i)
      dz= z-wz(i)
      dr2= dx*dx+dy*dy+dz*dz
      ix=ii(n)
      iy=ij(n)
      iz=ik(n)
      ie=ix+iy+iz
c
c make also limit of 0**0
      if(dx.ne.0.0d0.or.ix.ne.0)then
        tempx=dx**ix
      else
        tempx=1.0d0
      endif
      if(dy.ne.0.0d0.or.iy.ne.0)then
        tempy=dy**iy
      else
        tempy=1.0d0
      endif
      if(dz.ne.0.0d0.or.iz.ne.0)then
        tempz=dz**iz
      else
        tempz=1.0d0
      endif
c
c constant for all primitives
      if(onorm) then
          com= tempx*tempy*tempz * 
     +         dsqrt(koef/(dfact(ix)*dfact(iy)*dfact(iz)))
      else
          com= tempx*tempy*tempz
      endif
c
c  loop for primitives
c
      sum=0.0d0
      if(onorm) then
          do 10,i=1,ncontr
          zet=zeta(n,i)
          prim= dsqrt(zet*dsqrt(zet)*(4.0d0*zet)**ie) * dexp(-zet*dr2)
          sum= sum+ ccoef(n,i)*prim
10        continue
      else
          do 20,i=1,ncontr
          sum= sum+ ccoef(n,i)*dexp(-dr2*zeta(n,i))
20        continue
      endif
c
      tmaov=sum*com
      return
      end
      function tmgaus(zet,ix)
      implicit real*8 (a-h,p-z), integer  (i-n), logical (o)
      real*8 koef
c
c  calculates value of the integral of one dimens. gaussian and is 
c  called from calculations of overlap integrals
c
c     (2i-1)!!, i=0,10
      dimension dfact(0:10)
      data dfact /1.0d0,1.0d0,3.0d0,15.0d0,105.0d0,945.0d0,
     +            10395.0d0,135135.0d0,
     +     2027025.0d0,34459425.0d0,654729075.0d0/
c     (2/pi)**0.5 =
      data koef /0.797884560802865383521975673d0/
c
      if(and(ix,1).eq.0) then
         i=ix/2
         tmgaus= dfact(i)/(koef*dsqrt(0.5d0*zet)*(2.0d0*zet)**i)
      else
c  odd function
         tmgaus= 0.0d0
      endif
      return
      end
      subroutine tminte(integ,naos,tmao)
c naos is number of primitives
      implicit real*8 (a-h,p-z), integer  (i-n), logical (o)
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
      dimension tmao(naos*(naos+1)/2)
      common/blkin/x0(3,mxprim),zet(mxprim),ii(3,mxprim)
c
c  multiplies triangular transition matrix by one dimensional integrals
c  which are constant in all points of the grid
c
      index=0
      do 30,i=1,naos
      xi=x0(integ,i)
      iii=ii(integ,i)
      zi=zet(i)
      do 30,j=1,i
      index=index+1
      tdum = tmprin(xi,x0(integ,j),zi,zet(j),iii,ii(integ,j))
      tmao(index)=tmao(index)*tdum
30    continue
c
      return
      end
      subroutine tmnbas
      implicit real*8 (a-h,p-z), integer  (i-n), logical (o)
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
      parameter (maxtrm = 10)
      common/junk2/nstate,naos,imo,icore,melec,nuclei,
     +wc(maxat),wx(maxat),wy(maxat),wz(maxat),
     +new(maxorb),icen(maxorb),ii(maxorb),ij(maxorb),ik(maxorb),
     +ic(maxorb),ccoef(maxorb,maxtrm),zeta(maxorb,maxtrm)
c
c  calculates value of an contracted gaussian no. iao at position x,y,z
c
c  do for all atomic orbitals
c
      do 10, iii=1,naos
      i=new(iii)
c
c  first convert contr. coefficients for unnormalized primitives
c
      nprim=ic(i)
      do 20,j=1,nprim
      ccoef(i,j)=ccoef(i,j)*tmnorm(.false.,zeta(i,j),ii(i),ij(i),ik(i))
20    continue
c
c  then normalize the whole contracted gaussian
c  in the case of d,f normalize as dxx, fxxx - see also atoms.f
c  pozor, jeste zkontrolovat!!!!!
c
      ie=ii(i)+ij(i)+ik(i)
      sum=0.0d0
      do 30,j=1,nprim
      cj=ccoef(i,j)
      zj=zeta(i,j)
      do 30,k=1,j
c  zeta must be divided by two - tmnorm doesn't calculate integral from
c  a gaussian but from a square of it!!!; ie must not be divided!
      val=cj*ccoef(i,k)*tmnorm(.true.,(zj+zeta(i,k))*0.5d0,ie,0,0)
      if(j.ne.k) val=val+val
      sum=sum+val
30    continue
      sum=dsqrt(1.0d0/sum)
      do 40,j=1,nprim
      ccoef(i,j)=ccoef(i,j)*sum
40    continue
c
10    continue
c
      return
      end
      function tmnorm(orecsq,zet,ix,iy,iz)
      implicit real*8 (a-h,p-z), integer  (i-n), logical (o)
      real*8 koef
c
c  calculates value of the normalization coefficient of a prim. GTO 
c  (false), or square of its reciprococal value = integral over 
c  unnorm. GTO (true)
c
c     (2i-1)!!, i=0,5
      dimension dfact(0:5)
      data dfact /1.0d0,1.0d0,3.0d0,15.0d0,105.0d0,945.0d0/
c     (2/pi)**1.5 =
      data koef /0.507949087473927818336960627d0/
c
c
      ie=ix+iy+iz
      if(orecsq) then
         tmnorm= dfact(ix)*dfact(iy)*dfact(iz)/
     +          (koef*zet*dsqrt(zet)*(4.0d0*zet)**ie)
      else
         tmnorm=dsqrt(koef*zet*dsqrt(zet)*(4.0d0*zet)**ie/
     +          (dfact(ix)*dfact(iy)*dfact(iz)))
      endif
      return
      end
      subroutine tmpl0(tmao,tmp1,tmp2,vectp,vect,tmat,
     +   ldim1,ldim2,ldim3,ldim4,ldim5,ldim6)
c
      implicit real*8    (a-h,p-w),integer   (i-n),logical    (o)
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
      parameter (maxtrm = 10)
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
      dimension tmao(ldim1),tmp1(ldim2),tmp2(ldim3),vectp(ldim4)
      dimension tmat(ldim6),vect(ldim5)
c
c
c  plane definition for check of the integration coordinate
c from old version
c     common/junk/dummy(maxorb),plane(26)
c
c  here are positions of nuclei to check if the molecule is the same as
c  in tran. mom. calcn. Basis specified here is irrelevant, because
c  the basis of tran. matrix is stored in itmfil. Nuclei are important
c  because they are later projected to the picture (in utility module)
c  and to avoid misunderstandings in plane positioning. The charge of 
c  atoms and molecule is not checked (czan is not ok initialized).
c  
c   this remains
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
c  information from mrdci transition moment module
      common/junk2/nstate,iorbs,imo,icore,melec,nuclei,
     +wc(maxat),wx(maxat),wy(maxat),wz(maxat),
     +new(maxorb),icen(maxorb),ii(maxorb),ij(maxorb),ik(maxorb),
     +ic(maxorb),ccoef(maxorb,maxtrm),zeta(maxorb,maxtrm)
c  from input and communicating to tmpl1 (norbs)
      common/tmdata/iwanted,imultcoor,isymcoor,ianticoor,integcoor
     +,otran,norbs,msptm
c
c
c
      character*1 axis
      character*10 charwall
      dimension axis(3)
      data axis/'x','y','z'/
c
      call tmrdtm(vect,tmat,.true.,0)
      if(iwanted.le.0)then
         write(iwr,*)' tmpl0: specified state is invalid'
         write(iwr,*)' tmpl0: or no state was specified at all!'
       call caserr('transition matrix plot requested for unknown state')
      endif
      if(nstate.lt.iwanted)then
        write(iwr,*)' tmpl0: specified state is higher than number of'
        write(iwr,*)' tmpl0: states calculated in tran. moment section!'
       call caserr('transition matrix plot requested for unknown state')
      endif
c
c  now check geometry
      write(iwr,*)
      write(iwr,*)
      write(iwr,*)' Checking positions of nuclei with'
      write(iwr,*)' respect to the transition moment job.'
      eps=1.0d-7
      if(nat.ne.nuclei)call caserr(
     +  'inconsistent no. of atoms with respect to the moment job')
      do 100,i=1,nuclei
      if(dabs(wx(i)-c(1,i))+dabs(wy(i)-c(2,i))
     +     +dabs(wz(i)-c(3,i)).gt.eps)then
         write(iwr,*)' tmpl0: difference at nuclei no. ',i
         call caserr('inconsistency detected in new and old geometry')
      endif
100   continue
      write(iwr,*)' OK'
c
      call tmrdtm(vect,tmat,.false.,iwanted)
c
c ta geometrie na zacatku jobu je v a.u.
c scf MO jsou lcao vzhledem k normalizovanym kontrahovanym gaussianum
c na zacatku jobu se prepocitavaji kontr. koeficienty na nenorm.
c primit. gaussiany a pak se cely kontr. gaussian normalizuje.
c
      write(iwr,*)
      write(iwr,*)' Transforming transition density matrix to ao basis'
      tt0=cpulft(1)
c
c  first recalculate contraction from normalized primitives to 
c  unnormalised ones and normalize the whole contracted gaussian
c because to this normalised contr. the scf vectors refer
      call tmnbas
c
c  if integration is requested, we will calculate in basis of 
c  primitive GTOs and then make the transf. itsself
c
      if(integcoor.eq.0)then
       norbs=iorbs
      call tmtftm(norbs,imo,vect,tmat,tmao,tmp1,tmp2)
      else
       norbs=maketm(iorbs)
       call tmtvec(imo,iorbs,norbs,vect,vectp)
       call tmtftm(norbs,imo,vectp,tmat,tmao,tmp1,tmp2)
       call tminte(integcoor,norbs,tmao)
      endif
c
      write(iwr,1)cpulft(1)-tt0,charwall()
1     format(/' Transformation done in',f7.3,' econds at',a10,' wall'/)
      write(iwr,*)
     + ' Grid calculation for transition density matrix is running.'
      if(isymcoor.eq.0)then
      write(iwr,*)
     + ' No symmetrisation of the transition matrix was requested.'
      else
      write(iwr,99)axis(isymcoor)
99    format(
     + ' Symmetrisation of the transition matrix with respect to ',a1,
     + ' coordinate was requested.')
      endif
      if(ianticoor.eq.0)then
      write(iwr,*)
     +  ' No antisymmetrisation of the transition matrix was requested.'
      else
      write(iwr,991)axis(ianticoor)
991   format(
     + ' Antisymmetrisation of the transition matrix with respect to ',
     + a1,' coordinate was requested.')
      endif
      if(imultcoor.eq.0)then
      write(iwr,*)
     + ' No multiplication of the transition matrix with coordinate',
     + ' was requested.'
      else
      write(iwr,992)axis(imultcoor)
992   format(' Multiplication of the transition matrix with ',a1,
     + ' coordinate was requested.')
      endif
      if(integcoor.eq.0)then
      write(iwr,*)
     + ' No integration of the transition matrix was requested.'
      else
      write(iwr,993)axis(integcoor)
993   format(' Integration of the transition matrix over ',a1,
     + ' coordinate was requested.')
      endif
      write(iwr,*)
c
c  pozor jeste na integraci spolu s (anti)symetrizaci !!
c
      if(integcoor.eq.imultcoor.and.integcoor.ne.0) call caserr(
     + 'integration not possible for same coordinate as multiplication')
      if(isymcoor.ne.0.and.ianticoor.ne.0) call caserr(
     + 'simultaneous antisymetrisation and symetrisation not available')
c
c
c
c  the following check not implemented yet in this version
c
c
c      if(integcoor.ne.0)then
cc
cc set up axes of plane and checks direction of integration
cc
c      do 90 i=1,3
cc  d = c - a
c         dvec(i)=plane(6+i)-plane(i)
cc  e = b - a
c         evec(i)=plane(3+i)-plane(i)
c  90     continue
cc calculate orthogonal set of axis
c      call crprod(dvec,evec,fvec)
cc  e,g are in plane, f is perpendicular
c      oflag=.false.
c      do 91,i=1,3
c      if(i.ne.integcoor .and. dabs(fvec(i)).gt.1d-9) oflag=.true.
c91    continue
c      if(oflag)write(iwr,*)
c     +  ' Warning - direction of integration is not '
c      if(oflag)write(iwr,*)' perpendicular to the specified plane!'
cc
c      endif
cc
      return
      end
      function tmpl1(tmao,pppp,qqqq,rrrr)
c
      implicit real*8    (a-h,p-w),integer   (i-n),logical    (o)
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
c  store some large temporaries in a common block (mxprim.gt.maxorb)
c     common/scrtch/tmao(mxprim*(mxprim+1)/2)
      dimension tmao(*)
      common/tmdata/iwanted,imultcoor,isymcoor,ianticoor,integcoor,
     +              otran,norbs,msptm
c
c
       grid0=tmval(pppp,qqqq,rrrr,norbs,tmao,integcoor)
c
       if(ianticoor.eq.1)
     + grid0=0.5d0*(grid0-tmval(-pppp,qqqq,rrrr,norbs,tmao,integcoor))
       if(ianticoor.eq.2)
     + grid0=0.5d0*(grid0-tmval(pppp,-qqqq,rrrr,norbs,tmao,integcoor))
       if(ianticoor.eq.3)
     + grid0=0.5d0*(grid0-tmval(pppp,qqqq,-rrrr,norbs,tmao,integcoor))
c
       if(isymcoor.eq.1)
     + grid0=0.5d0*(grid0+tmval(-pppp,qqqq,rrrr,norbs,tmao,integcoor))
       if(isymcoor.eq.2)
     + grid0=0.5d0*(grid0+tmval(pppp,-qqqq,rrrr,norbs,tmao,integcoor))
       if(isymcoor.eq.3)
     + grid0=0.5d0*(grid0+tmval(pppp,qqqq,-rrrr,norbs,tmao,integcoor))
c
       if(imultcoor.eq.1)grid0=grid0*pppp
       if(imultcoor.eq.2)grid0=grid0*qqqq
       if(imultcoor.eq.3)grid0=grid0*rrrr
       tmpl1=grid0
c
      return
      end
      function tmprin(xa,xb,za,zb,la,lb)
      implicit real*8 (a-h,o-z), integer  (i-n)
c
c  calculates value of the integral over product of two generally
c  differently centered one dimensional gaussians with different l's 
c  and zeta's
c
c     g*g
      parameter(maxn=8)
      dimension ibinom(0:maxn,0:maxn)
      data ibinom/1,0,0,0,0,0,0,0,0,
     +            1,1,0,0,0,0,0,0,0,
     +            1,2,1,0,0,0,0,0,0,
     +            1,3,3,1,0,0,0,0,0,
     +            1,4,6,4,1,0,0,0,0,
     +            1,5,10,10,5,1,0,0,0,
     +            1,6,15,20,15,6,1,0,0,
     +            1,7,21,35,35,21,7,1,0,
     +            1,8,28,56,70,56,28,8,1/
c
c
      z=za+zb
      dx=xa-xb
      l=la+lb
      if(dabs(dx).le.1.0d-9)then
c  identical center
      tmprin=tmgaus(z,l)
      return
      else
c  different centers
      dx2=dx*dx
      xp=(za*xa+zb*xb)/z
c  factor from conversion of exponential parts
      expcoef=dexp(-dx2*za*zb/z)
c
c  make decomposition over l=0...la+lb
c
      sum=0.0d0
      do 10,k=0,l,2
c
c  calculate decomposition coefficient
c
      factor=0.0d0
      pa=xp-xa
      pb=xp-xb
      ih=min(k,la)
      il=max(0,k-lb)
      do 20,i=il,ih
      j=k-i
      factor=factor+ pa**(la-i)*pb**(lb-j)*ibinom(i,la)*ibinom(j,lb)
20    continue
c
      sum=sum+factor*tmgaus(z,k)
10    continue
c
      tmprin=sum*expcoef
      return
      endif
      end
      subroutine tmrdtm(vect,tmat,oreadbas,istate)
      implicit real*8 (a-h,p-z), integer  (i-n), logical (o)
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
      parameter (maxtrm = 10)
      common/junk2/nstate,iorbs,imo,icore,melec,nuclei,
     +wc(maxat),wx(maxat),wy(maxat),wz(maxat),
     +new(maxorb),icen(maxorb),ii(maxorb),ij(maxorb),ik(maxorb),
     +ic(maxorb),ccoef(maxorb,maxtrm),zeta(maxorb,maxtrm)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common/ftape/lscr(20),itmfil
      dimension vect(*),tmat(*)
c
c
      if(oreadbas) then
c
c read important information for tmatrix conversion to cartesian 
c representation
c
      rewind(itmfil)
c  basic parameters
      read(itmfil) nstate,iorbs,imo,icore,melec,nuclei
c  geometry
      read(itmfil) (wc(j),j=1,nuclei)
      read(itmfil) (wx(j),j=1,nuclei)
      read(itmfil) (wy(j),j=1,nuclei)
      read(itmfil) (wz(j),j=1,nuclei)
c
c  basis
c     new is permutation of 1..iorbs
      read(itmfil) (new(j),j=1,iorbs)
      onew=.false.
      do i=1,iorbs
      if(new(i).ne.i) onew=.true.
      enddo
      if(onew) then
      write(iwr,*)' tmrdtm: array new is NOT identical permutation'
      write(iwr,32)(new(i),i=1,iorbs)
32    format(20i5)
      else
      write(iwr,*)' tmrdtm: array new is identical permutation'
      endif
c centers
      read(itmfil) (icen(j),j=1,iorbs)
c x,y,z exponents
      read(itmfil) (ii(j),j=1,iorbs)
      read(itmfil) (ij(j),j=1,iorbs)
      read(itmfil) (ik(j),j=1,iorbs)
c number of primitives
      read(itmfil) (ic(j),j=1,iorbs)
c contraction coefficients and exponents
      do 8000, j=1,iorbs
      do 8001, k=1,ic(j)
      read(itmfil) ccoef(j,k),zeta(j,k)
8001  continue
8000  continue
c
c
c
c  read vect(ao,mo) coefficients - coefficients of transformed mos to 
c  aos. core and discarded orbitals are skipped and the rest is sorted 
c  according to irreps
c
      ibase=0
      do 8004,i=1,imo
      read(itmfil)(vect(ibase+j),j=1,iorbs)
      ibase=ibase+iorbs
8004  continue
c
c
      endif
c
c
c read tmat; the istate-th one since the last call of tmrdtm
c
c read in columns
c
      do 8003,i=1,istate
      jpoint=0
      do 8003,j=1,imo
      read(itmfil) (tmat(k+jpoint),k=1,imo)
      jpoint=jpoint+imo
8003  continue
c
      return
      end
      subroutine tmtftm(naos,nmos,vect,tmat,tmao,tmp1,tmp2)
      implicit real*8 (a-h,p-z), integer  (i-n), logical (o)
      dimension vect(naos,nmos),tmat(nmos,nmos),tmao(naos*(naos+1)/2),
     +  tmp1(nmos,naos),tmp2(naos,naos)
c
c  transforms tmat into ao basis and converts it to the symetrized 
c  triangular form: vect can be here coef to normalised contracted 
c  or unnormalised primitive gaussians
c
      do 10,j=1,naos
      do 10,i=1,nmos
      val=0.0d0
      do 20,k=1,nmos
20    val=val+tmat(i,k)*vect(j,k)
10    tmp1(i,j)=val
      call mxm(vect,naos,tmp1,nmos,tmp2,naos)
c
c  here convert to the triangular
c  because we need only the diagonal (in cart. coords. - gamma(r|r)) 
c  element and therefore
c  it is not necessary to keep the whole matrix in ao representation
c  but it was necessary in mo representation
c
      index=0
      do 30,i=1,naos
      do 30,j=1,i
      index=index+1
      if(i.eq.j) then
        tmao(index)=tmp2(j,i)
      else
        tmao(index)=tmp2(i,j)+tmp2(j,i)
      endif
30    continue
c
      return
      end
      subroutine tmtvec(nmos,naos,npri,vect,vectp)
      implicit real*8 (a-h,p-z), integer  (i-n), logical (o)
      dimension vect(naos,nmos),vectp(npri,nmos)
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
      parameter (maxtrm = 10)
      common/junk2/idum(6),
     +wc(maxat),wx(maxat),wy(maxat),wz(maxat),
     +new(maxorb),icen(maxorb),ii(maxorb),ij(maxorb),ik(maxorb),
     +ic(maxorb),ccoef(maxorb,maxtrm)
      common/blkin/x0(3,mxprim),zet(mxprim),
     + inew(3,mxprim),ipoint(maxorb),jpoint(mxprim),jcount(mxprim)
c
c
c  transforms scf vectors into primitive basis
c  vectp will refer to unnormalised primitives
c  vect refered to normalised contracted aos
c
c
      do 10,i=1,nmos
      do 10,j=1,npri
      jbase=jpoint(j)
      vectp(j,i)=vect(jbase,i)*ccoef(new(jbase),jcount(j))
10    continue
      return
      end
      function tmval(x,y,z,naos,tmao,integ)
c
c  calculates value of an diagonal element of the first order
c  reduced transition density matrix in coordinate representation
c  at point x,y,z
c  optionally integrated over integ coordinate (now must be different 
c  from the polarisation coordinate)
c  if integ.ne.0 then 
c  tmao is the tran dens. mat. in prim. GTO basis symetrised 
c  triangular and already multiplied by the constant 
c (in all points of grid) one dimensional integrals
c
      implicit real*8 (a-h,p-z), integer  (i-n), logical (o)
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
      dimension tmao(naos*(naos+1)/2)
c
      dimension valuesao(mxprim),xx(3)
c
c  branch according to integ
c
      if(integ.eq.0) then
c
c  first calculate aos over contracted gaussians, false because 
c  ccoefs refer to unnormalised primitives
c  renormalisation of basis was called before - these values will be 
c  of normalised contracted gaussians
c
      do 10, i=1,naos
      valuesao(i)=tmaov(i,x,y,z,.false.)
10    continue
      else
c
c  first calculate parts of aos over remaining 2 dimensions of 
c  primitive gaussians
c
      xx(1)=x
      xx(2)=y
      xx(3)=z
      do 80,i=1,naos
      valuesao(i)=tmval2(i,xx,integ)
80    continue
      endif
c
c
c and then calc tmatrix in the coordinate representation
c the matrix is in triangular form
c
      sum=0.0d0
      index=0
      do 50,i=1,naos
      val=0.0d0
      do 40,j=1,i
      index=index+1
      val=val+tmao(index)*valuesao(j)
40    continue
      sum=sum+val*valuesao(i)
50    continue
c
      tmval=sum
      return
c
      end
      function tmval2(n,xx,integ)
c
c  calculates the two dim. after integration remaining part of  
c  primit. GTO value of x corresponding to integ is ignored
c
      implicit real*8 (a-h,p-z), integer  (i-n), logical (o)
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
      dimension xx(3)
      common/blkin/x0(3,mxprim),zeta(mxprim),ii(3,mxprim)
c
c
      com=1.0d0
      dr2=0.0d0
      do 10,i=1,3
      if(i.ne.integ) then
          dx=xx(i)-x0(i,n)
          ixi=ii(i,n)
          if(dx.ne.0.0d0 .or. ixi.ne.0) then
            com=com* dx**ixi
          endif
          dr2=dr2+dx*dx
      endif
10    continue
c
      tmval2=com*dexp(-dr2*zeta(n))
      return
      end
**==tranpr.f
      subroutine tranpr(dint,dint11,consf,newbas,
     *ilifc,nterm,it,ct,otran,lenbas)
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
      dimension dint(*),dint11(*)
      dimension ilifc(*),nterm(*),it(*),ct(*)
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      common/miscop/ct1(mxorb3),it1(mxorb3)
      data mm3/-3/
c
      if (.not.otran) then
         m = 0
         len3 = lenbas + mm3
         do 60 i = 1 , newbas
            nt1 = nterm(i)
            ii = ilifc(i)
            do 20 k = 1 , nt1
               l = k + ii
               ct1(k) = ct(l)
               it1(k) = it(l)
 20         continue
            do 50 j = 1 , i
               m = m + 1
               dum = 0.0d0
               nt2 = nterm(j)
               jj = ilifc(j)
               do 40 k = 1 , nt2
                  kk = k + jj
                  l = it(kk)
                  do 30 ii = 1 , nt1
                     ll = it1(ii)
                     if (ll.ge.l) then
                        top = dint(l+iky(ll))
                     else
                        top = consf*dint(ll+iky(l))
                     end if
                     dum = dum + ct1(ii)*ct(kk)*top
 30               continue
 40            continue
               dint11(m) = dum
 50         continue
 60      continue
         call dcopy(len3,dint11(1),1,dint(1),1)
      end if
c
      return
      end
**==unstff.f
      subroutine unstff(iq,ia,n)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension iq(*),ia(*)
      common/stuf/nwrd
      nav = lenwrd()
      do 10 i = 1,n*nav
 10      ia(i)=iq(nwrd+i)
      nwrd = nwrd + n*nav
      return
      end
**==unstfi.f
      subroutine unstfi(iq,ia,n)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension iq(*),ia(*)
      common/stuf/nwrd
      do 10 i = 1,n
 10      ia(i)=iq(nwrd+i)
      nwrd = nwrd + n
      return
      end
      subroutine vdwf(grid, data, nvdw, iuvdw, ruvdw, rad, olog, 
     +                ograd, iwr)
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
      real*8 radvdw
      common/vdwrad/radvdw(105)
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
      dimension grid(*), data(*), rad(*), iuvdw(*), ruvdw(*), v(3)

      if(.not.ograd)then
         nel=1
         if(.not.olog)then
            write(iwr,1000)
         else
            write(iwr,1005)
         endif
      else
         nel=3
c         if(.not.olog)then
            write(iwr,1002)
c         else
c            write(iwr,1007)
c         endif
      endif
         write(iwr,1008)
 1000 format(/1x,'van der waals function calculation: ',
     &     ' f = max(rvdw-r)')
 1002 format(/1x,'generation of vdw normals ')
c     &     ' f = max(rvdw-r)')
 1005 format(/1x,'log van der waals function calculation: ',
     &     ' f = max(-ln(r/rvdw))')
c1007 format(/1x,'grad log van der waals function calculation: ')
c     &     ' f = max(-ln(r/rvdw))')
 1008 format(/1x,'table of van der waals radii',//,12x,
     &     'at #  tag        z    rvdw (au)',/12x,30('-'))
      do 100 i = 1,nat
         izi=nint(czan(i))
         if(izi.le.0)goto 100
         ouser = .false.
         do 50 j = 1, nvdw
            if(zaname(i).eq.zaname(iuvdw(j)))then
               ouser = .true.
               rad(i) = ruvdw(j)
               write(iwr,1010)i,zaname(i),izi,rad(i),'(user defined)'
            endif
 50      continue
         if(.not.ouser)then
            rad(i) = radvdw(izi)
            if(radvdw(izi).eq.0.0d0)then
               write(iwr,1010)i,ztag(i),izi,rad(i),'** undefined **'
            else
               write(iwr,1010)i,ztag(i),izi,rad(i),' '
            endif
         endif
 100  continue
 1010 format(10x,i3,5x,a8,1x,i3,f12.6,2x,a15)
      call getnpt(np,npx,npy,npz)
      do 10 ipt = 1,np
         call getpt(ipt,ax,ay,az,grid,nel,idat)
         f0 = -1.0d20
         do 200 i = 1,nat
            if(izi.le.0)goto 200
            r=dsqrt((ax-c(1,i))**2 +(ay-c(2,i))**2 + (az-c(3,i))**2)
            if(olog)then
               r = max(r,1.0d-10)
               f = -dlog(r/rad(i))
            else
               f = rad(i) - r
            endif
            if(f.ge.f0)then
               f0 = f
               inear = i
            endif
 200     continue
         if(ograd)then
            v(1) = ax - c(1,inear)
            v(2) = ay - c(2,inear)
            v(3) = az - c(3,inear)
            call normal(v)
            data(idat+0) = v(1)
            data(idat+1) = v(2)
            data(idat+2) = v(3)
         else
            data(idat) = f0
         endif
 10   continue
      return
      end
      subroutine wcasec(icalc,cq,iwr)
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
c  routines to put calculation results onto dumpfile
c  (we do not force the user to put the grid there as well
c  but we save the section if present
c
c grid definition parameters
c
      common/dfgrid/geom(12,mxgrid),gdata(5,mxgrid),igtype(mxgrid),
     &             igsect(mxgrid),nptot(mxgrid),npt(3,mxgrid),
     &             igdata(5,mxgrid),igstat(mxgrid),ngdata(mxgrid),
     &             ngrid
c
c data calculation parameters
c
      common/dfcalc/cdata(5,mxcalc),ictype(mxcalc),icsect(mxcalc),
     &            icgrid(mxcalc),icgsec(mxcalc),ndata(mxcalc),
     &            icstat(mxcalc),icdata(5,mxcalc),ncalc
      common/dfclc2/cscal(mxcalc),iblkp(mxcalc)
c
c labels and titles
c
      common/cplot/zgridt(10,mxgrid),zgrid(mxgrid),
     &             zcalct(10,mxcalc),zcalc(mxcalc),
     &             zplott(10,mxplot),zplot(mxplot)
      dimension cq(*)
      common/psh/crap(100)
      integer numdu, iblkdu, maxb, iblkla
      logical orevis
      real*8 apos
      common/sector/numdu,iblkdu,orevis(2),
     + apos(508),maxb,iblkla
      data itcalc/51/
c
c create a grid data dumpfile section - irregular data in cq
c
c @ check icgsec
      if(icstat(icalc).ne.2)
     &     call caserr('attempt to write grid def. before calculation')

      igrid = icgrid(icalc)
      call stuf1
      call stufi(crap,ictype(icalc),1)
      call stufi(crap,ndata(icalc),1)
      call stufi(crap,icdata(1,icalc),5)
      call stuff(crap,cdata(1,icalc),5)
      call stuff(crap,cscal(icalc),1)
      call stuff(crap,geom(1,igrid),12)
      call stufi(crap,igtype(igrid),1)
      call stufi(crap,nptot(igrid),1)
      call stufi(crap,npt(1,igrid),3)
      call stufi(crap,igdata(1,igrid),5)
      call stuff(crap,gdata(1,igrid),5)
      call stufi(crap,igstat(igrid),1)
      call stufi(crap,igsect(igrid),1)
      call stufi(crap,ngdata(igrid),1)
      nw = nstuf()
      nwc = ndata(icalc)

      lengrd = 1 + lensec(nw) + lensec(nwc) + lensec(10)

      call secput(icsect(icalc),itcalc,lengrd,iblk)
      write (iwr,1000) icsect(icalc) , iblk

      call wrt3(crap,nw,iblk,numdu)
      call wrtcs(zcalct(1,icalc),10,numdu)
      call wrt3s(cq,nwc,numdu)
      call clredx
      call revind
      return
 1000 format (/' route grid data to dumpfile :section ',i3,
     &' (at block ',i6,')'/)
      end
      subroutine wgrsec(igrid,gq,iwr)
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
c grid definition parameters
c
      common/dfgrid/geom(12,mxgrid),gdata(5,mxgrid),igtype(mxgrid),
     &             igsect(mxgrid),nptot(mxgrid),npt(3,mxgrid),
     &             igdata(5,mxgrid),igstat(mxgrid),ngdata(mxgrid),
     &             ngrid
c
c labels and titles
c
      common/cplot/zgridt(10,mxgrid),zgrid(mxgrid),
     &             zcalct(10,mxcalc),zcalc(mxcalc),
     &             zplott(10,mxplot),zplot(mxplot)
      integer numdu, iblkdu, maxb, iblkla
      logical orevis
      real*8 apos
      common/sector/numdu,iblkdu,orevis(2),
     + apos(508),maxb,iblkla
      common/psh/crap(100)
      data itgrid/50/
      if(igstat(igrid).ne.2)call
     &  caserr('trying to write grid to dumpfile before calculation')
c
c create a grid dumpfile section - irregular data in gq
c
      call stuf1
      call stuff(crap,geom(1,igrid),12)
      call stufi(crap,igtype(igrid),1)
      call stufi(crap,nptot(igrid),1)
      call stufi(crap,npt(1,igrid),3)
      call stufi(crap,igdata(1,igrid),5)
      call stuff(crap,gdata(1,igrid),5)
      call stufi(crap,igstat(igrid),1)
      call stufi(crap,ngdata(igrid),1)
      nw = nstuf()

      nwg = ngdata(igrid)
      lengrd = 1 + lensec(nw) + lensec(10)
      if(nwg.ne.0)lengrd = lengrd + lensec(nwg)

      call secput(igsect(igrid),itgrid,lengrd,iblk)
      write (iwr,1000) igsect(igrid) , iblk

      call wrt3(crap,nw,iblk,numdu)
      call wrtcs(zgridt(1,igrid),10,numdu)
      if(nwg.ne.0)call wrt3s(gq,nwg,numdu)

      call clredx
      call revind
      return
 1000 format (/' route grid definition to dumpfile :section ',i3
     &     ,' (at block ',i6, ')'/)
      end
**==wrap.f
      subroutine wrap(vect, ctran, czeta, pop, s1, s2,
     +     startx, starty, startz,
     +     pointx, pointy, pointz, gr, index, bigg, small,value, iwr)
      implicit real*8  (a-h,p-w),integer   (i-n)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c modified wrap - single contour, no gradient
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
      dimension vect(*), ctran(*), czeta(*), pop(*), gr(3,*),
     &     s1(*), s2(*)
c
c Accuracy for finding contours
      accry = 0.01d0
c factor to overshoot interpolations
      del = 0.02d0
c
c at present density iso only
c
      imo = 0
      call dplota(vect, ctran, czeta, pop, s1, s2,
     + startx, starty, startz, stden,imo)
      bigg = max (bigg,stden)
      small = min(small,stden)
c
c     Take a unit step
c
      tryx0 = pointx
      tryy0 = pointy
      tryz0 = pointz
      call dplota(vect, ctran, czeta, pop, s1, s2,
     &     tryx0, tryy0, tryz0, trydn0,imo)
      bigg = max (bigg,trydn0)
      small = min(small,trydn0)

      tryden = trydn0
      tryx=tryx0
      tryy=tryy0
      tryz=tryz0
      maxstp = 15
c
c     Initialise the bisection
c
      ihi = 0
      ilo = 0
      if(stden.gt.value .and. tryden .lt.value ) then
         denmax = stden
         stpmax = 0.0d0
         denmin = tryden
         stpmin = 1.0d0
      else if ( stden .lt. value .and. tryden .gt. value ) then
         denmax = tryden
         stpmax = 1.0d0
         denmin = stden
         stpmin = 0.0d0
      else
c     Initialisation failed try the next contour
         goto 600
      endif
c
       nstep = 0

 100   continue
c
c     test on convergence
c
       convx = (stpmax - stpmin) * ( pointx - startx )
       convy = (stpmax - stpmin) * ( pointy - starty )
       convz = (stpmax - stpmin) * ( pointz - startz )
       conv = convx * convx + convy * convy + convz * convz
       if ( conv .lt. accry*accry) goto 200
*
*     Avoid an infinite loop
*
       nstep = nstep + 1
       if(nstep .gt. maxstp ) then
          write(iwr,6000)
          goto 300
       end if
c
c     No convergence yet, interpolate
c

       if(denmax-denmin.lt.1.0d-30)then
          t = 0.5d0
       else
          t = (denmax-value)/(denmax-denmin)
       endif

       if(ilo.gt.1)then
          t = t + del
       else if (ihi.gt.1)then
          t = t - del
       endif

c       write(iwr,*)nstep,denmin,denmax,stpmin,stpmax,t
       stpnew = t*stpmin + (1.0d0-t)*stpmax

       tryx = startx + stpnew * ( pointx - startx )
       tryy = starty + stpnew * ( pointy - starty )
       tryz = startz + stpnew * ( pointz - startz )
c
       call dplota(vect, ctran, czeta, pop, s1, s2,
     &      tryx, tryy, tryz, tryden,imo)

       if( tryden .gt. value) then
          ihi = 0
          ilo = ilo+1
          denmax = tryden
          stpmax = stpnew
       else
          ilo = 0
          ihi = ihi+1
          denmin = tryden
          stpmin = stpnew
       endif
       goto 100
c
c     Converged or too many iterations
c
 200   continue
       index = index+1
       gr(1,index)=tryx
       gr(2,index)=tryy
       gr(3,index)=tryz
 300   continue
 600   continue
       return
 6000  format(/1x,'too many steps taken, convergence failed'/)
       end
c
c  generate a grid by point selection from another grid, based on
c  data values.
c
      subroutine selgrd(grid, data, ndata, vlow, vhigh, imap, ntot,
     &     iflag,iw)
      implicit real*8 (a-h,p-w),integer   (i-n)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      dimension grid(3,*), data(*), imap(*)
c
      write(iw,1000)vlow,vhigh
c
      n = 0
      do 10 i = 1, ndata
         if(data(i).ge.vlow .and. data(i).le.vhigh)then
            n = n + 1
            imap(n)=i
         endif
 10   continue
c
      do 30 i = 1, n
         do 20 j = 1, 3
            call getpt(imap(i),ax,ay,az,grid,1,idat)
            grid(1,i)=ax
            grid(2,i)=ay
            grid(3,i)=az
 20      continue
 30   continue
c      
      ntot = n
c
c  is there an edge cell in the list ( only apply tests to 3d grids)
c
      call getgty(it)
      if(it.eq.2)then
         nwarn=0
         call getnpt(np,nx,ny,nz)
         do 40 j = 1,n
            call get3ix(imap(j),ii,ij,ik)
            if(ii.eq.1 .or. ij.eq.1 .or. ik.eq.1 .or.
     &         ii.eq.nx .or. ij.eq.ny .or. ik.eq. nz)
     &         nwarn=nwarn+1
 40      continue
         if(nwarn .ne. 0)then
            if(iflag.eq.0)then
               write(iw,1020)'*ERROR*',nwarn
               write(iw,1030)
               call caserr('grid select size error')
            else
               write(iw,1020)'WARNING',nwarn
            endif
         endif
      endif
c
      write(iw,1010)ntot,ndata
c
 1000 format(/,1x,'select grid points by selection data values between',
     &     g14.6, ' and ',g14.6,' from current dataset')
 1010 format(/,1x,i7,' points were selected by data value out of a ',
     &     'a total of ',i7)
 1020 format(//,1x,a7,' -',i4,
     &     ' points selected are in edge cells ',
     &     /,1x,' check that the 3D grid you started with ',
     &     'is large enough')
 1030 format(/,1x,' add the keyword nocheck to the type select ',
     &     'directive if you want the job to ',/,1x,
     &     'continue past this check.')
      return
      end
      subroutine ver_analb(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/analb.m,v $
     +     "/
      data revision /
     +     "$Revision: 6299 $"
     +     /
      data date /
     +     "$Date: 2014-10-06 16:33:14 +0200 (Mon, 06 Oct 2014) $"
     +     /
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
