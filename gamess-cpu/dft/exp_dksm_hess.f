c***********************************************************************
c
c  Exp_dksm_hess.m
c
c  This file contains the subroutines that calculate the explicit 
c  derivatives with respect to the nuclear coordinates needed for a
c  Hessian calculation. The quantities calculated are the explicit
c  1st derivative of the Kohn-Sham matrix and the explicit 2nd 
c  derivative of the functional.
c
c***********************************************************************
c
      subroutine ver_dft_exp_dksm_hess(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/dft/exp_dksm_hess.m,v $
     +     "/
      data revision /
     +     "$Revision: 5774 $"
     +      /
      data date /
     +     "$Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $"
     +     /
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine den_pert_exp_ao(rkstyp_sw,gradcorr_sw,ao_tag,
     &     nao,npts,npert,ngridcentres,bfn_val,bfng_val,bfn_hess,
     &     adens,bdens,drhodb,dgrhodb,mxp)
      implicit none
c
c     Computes the derivatives of rho, grho, and gamma with respect
c     to a nuclear coordinate. Only the terms stemming from 
c     differentiating the explicit functions of the geometry are 
c     included. 
c
c     Parameters:
c
C *Parameter module for CCP1/DFT            
C *Memory information		
C *----------------
C * max_block 	- maximum number of blocks which can be allocate
      integer max_block
      parameter(max_block=20)
C *							
C *Basis sets information                              
C *----------------------				
C *max_tag 	- maximum number of basis sets which can be alloc
C *max_atype	- maximum number of centre types	
C *max_gtype	- maximum number of grid types (at least every element 
C                 has different grid type)
C *max_grids    - maximum number of grids (different terms may have
C                 different grids, i.e. CPKS equations use a different
C                 grid than the one used for the KS-matrix)
C *max_shel	- maximum number of shells on centre
C *max_prm	- maximum number of primitives for any given centre
C *maxL		- maximum angular momentum allowed	
C *max_func	- maximum number of basis functions for any centre
      integer max_tag,max_atype,max_shel,max_prm,max_ang
      integer max_gtype,max_grids
      parameter(max_tag=3,max_atype=30,max_shel=500,max_prm=5000)
      parameter(max_ang=5,max_gtype=10,max_grids=2)
C *								
C *Geometry information					
C *--------------------					
C *max_atom	- maximum number of atoms in system
      integer max_atom
      parameter(max_atom=750)
C *
C *Accuracy information				
C ---------------------			
C *global_accuracy 	- global accuracy	
      real*8  global_accuracy
      parameter(global_accuracy=1.0d-14)

C *Grid information
c *----------------
c *maxradzn - The maximum number of radial zones. 
      integer maxgpt,maxrad,maxfpt,maxang,maxradzn,maxtablerows
      parameter(maxgpt=2900,maxrad=50,maxang=302,maxfpt=100)
      parameter(maxradzn=35,maxtablerows=7)
c
c     The parameters defined here label the 1st, 2nd and 3rd derivative
c     components of the basis functions.
c
c     2nd derivatives:
c
      integer hxx, hyy, hzz, hxy, hxz, hyz
c
      parameter (hxx = 1)
      parameter (hyy = 2)
      parameter (hzz = 3)
      parameter (hxy = 4)
      parameter (hxz = 5)
      parameter (hyz = 6)
c
      integer hcc(3,3)
      data hcc/hxx,hxy,hxz,
     &         hxy,hyy,hyz,
     &         hxz,hyz,hzz/
c
c     3rd derivatives:
c
      integer txxx, txxy, txxz, txyy, txyz
      integer txzz, tyyy, tyyz, tyzz, tzzz
c
      parameter (txxx =  1)
      parameter (txxy =  2)
      parameter (txxz =  3)
      parameter (txyy =  4)
      parameter (txyz =  5)
      parameter (txzz =  6)
      parameter (tyyy =  7)
      parameter (tyyz =  8)
      parameter (tyzz =  9)
      parameter (tzzz = 10)
c
      integer tc3(3,3,3)
      data tc3/txxx,txxy,txxz,
     &         txxy,txyy,txyz,
     &         txxz,txyz,txzz,
     &         txxy,txyy,txyz,
     &         txyy,tyyy,tyyz,
     &         txyz,tyyz,tyzz,
     &         txxz,txyz,txzz,
     &         txyz,tyyz,tyzz,
     &         txzz,tyzz,tzzz/
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
c     In variables:
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
C *Basis sets include file
c
      real*8  alpha
      real*8  cont_coeff
      integer num_bset,bset_tags
      common/basis_sets/alpha(max_tag,max_atype,max_prm),
     &                  cont_coeff(max_tag,max_atype,max_prm,max_ang),
     &                  num_bset,bset_tags(max_tag)

      integer Ashl, Aprm, Abfn, totshl, totprm, totbfn, size_shlA
      integer size_basA, size_primA, maxi_shlA, maxi_basA, maxi_primA
      integer num_types, atom_tag , num_shl, atm_typ, nprim, angmom
      integer hybrid, pstart

      common/basis_size_info/Ashl(max_tag,max_atype),
     &                       Aprm(max_tag,max_atype),
     &                       Abfn(max_tag,max_atype),
     &                       totshl(max_tag),
     &                       totprm(max_tag),
     &                       totbfn(max_tag),
     &                       size_shlA,maxi_shlA,
     &                       size_basA,maxi_basA,
     &                       size_primA,maxi_primA
C
C Descriptions
C 
C  num_types		-	number of types of atoms for a given basis
C  atom_tag		-	list of atom type id numbers
C
      common/basis_cent_info/num_types(max_tag),
     &                       atom_tag(max_tag,max_atom)
C
C Descriptions
C
C  num_shl		-	number of shells for given atom type
C  atm_typ		-	atomic number of atom type
C  nprim		-  	number of primitives 
C  angmom		-	angular momentum
C  hybrid		-	level of hybrid shell. Is same as angmom if
C				shell is not a hybrid. For hybrid shells is
C				always less than angmom e.g. for an sp shell
C				angmom=2, hybrid=1
C  pstart		-	start of exponents and contraction coeffs 
C				contained in include basis.hf77, for a given
C				shell
C
      common/basis_cont_info/num_shl(max_tag,max_atype),
     &                       atm_typ(max_tag,max_atype),
     &                       nprim(max_tag,max_atype,max_shel),
     &                       angmom(max_tag,max_atype,max_shel),
     &                       hybrid(max_tag,max_atype,max_shel),
     &                       pstart(max_tag,max_atype,max_shel)

c
      logical rkstyp_sw ! .true. if closed shell
      logical gradcorr_sw
      integer ao_tag
      integer nao, npts, npert, ngridcentres, mxp
      real*8 bfn_val(mxp,nao)
      real*8 bfng_val(mxp,nao,3)
      real*8 bfn_hess(mxp,nao,6)
      real*8 adens(nao*(nao+1)/2)
      real*8 bdens(nao*(nao+1)/2)
c
c     Out variables:
c
      real*8 drhodb(mxp,2,npert), dgrhodb(mxp,2,3,npert)
c
c     Local variables:
c
      integer bstart, bend, latm, lpert, atmt, mu, nu, iki
      integer ipt
      integer k, l
c
c     Code:
c
c     For the time being I don't want to worry about using translational
c     invariance.
      if (npert.ne.3*ngridcentres) call caserr('not implemented yet 1')
c     if (gradcorr_sw) call caserr('not implemented yet 2')
c
      call aclear_dp(drhodb,mxp*2*npert,0.0d0)
      if (gradcorr_sw) then
         call aclear_dp(dgrhodb,mxp*2*3*npert,0.0d0)
      endif
c
      if (rkstyp_sw) then
         bend = 0
         do 110 latm=1,ngridcentres
            lpert=3*(latm-1)
            atmt = atom_tag(ao_tag,latm)
            bstart=bend+1
            bend=bend+Abfn(ao_tag,atmt)
c           if (gwt_sw.and.latm.eq.iatm) goto 110
            do mu=bstart,bend
               iki=iky(mu)
               do nu=1, mu
                  do k=1,3
                     do ipt=1,npts
                        drhodb(ipt,1,lpert+k) = drhodb(ipt,1,lpert+k)
     &                  - adens(iki+nu)*bfng_val(ipt,mu,k)*
     &                    bfn_val(ipt,nu)
                     enddo ! ipt
                  enddo ! k
               enddo ! nu
               do nu=mu+1, nao
                  iki=iky(nu)
                  do k=1,3
                     do ipt=1,npts
                        drhodb(ipt,1,lpert+k) = drhodb(ipt,1,lpert+k)
     &                  - adens(iki+mu)*bfng_val(ipt,mu,k)*
     &                    bfn_val(ipt,nu)
                      enddo ! ipt
                  enddo ! k
               enddo ! nu
            enddo ! mu
 110     continue ! latm
         if (gradcorr_sw) then
            bend = 0
            do 120 latm=1,ngridcentres
               lpert=3*(latm-1)
               atmt = atom_tag(ao_tag,latm)
               bstart=bend+1
               bend=bend+Abfn(ao_tag,atmt)
c              if (gwt_sw.and.latm.eq.iatm) goto 120
cDEBUG
c        write(*,*)'latm  =',latm
c        write(*,*)'atmt  =',atmt
c        write(*,*)'bstart=',bstart
c        write(*,*)'bend  =',bend
cDEBUG
               do mu=bstart,bend
                  iki=iky(mu)
                  do nu=1, mu
                     do ipt=1,npts
                        do l=1,3
                           do k=1,3
                              dgrhodb(ipt,1,k,lpert+l) 
     &                        = dgrhodb(ipt,1,k,lpert+l)
     &                        - adens(iki+nu)*bfng_val(ipt,mu,l)*
     &                          bfng_val(ipt,nu,k)
                           enddo ! k
                        enddo ! l
c
                        dgrhodb(ipt,1,1,lpert+1) 
     &                  = dgrhodb(ipt,1,1,lpert+1)
     &                  - adens(iki+nu)*bfn_hess(ipt,mu,hxx)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,1,1,lpert+2) 
     &                  = dgrhodb(ipt,1,1,lpert+2)
     &                  - adens(iki+nu)*bfn_hess(ipt,mu,hxy)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,1,1,lpert+3) 
     &                  = dgrhodb(ipt,1,1,lpert+3)
     &                  - adens(iki+nu)*bfn_hess(ipt,mu,hxz)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,1,2,lpert+1) 
     &                  = dgrhodb(ipt,1,2,lpert+1)
     &                  - adens(iki+nu)*bfn_hess(ipt,mu,hxy)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,1,2,lpert+2) 
     &                  = dgrhodb(ipt,1,2,lpert+2)
     &                  - adens(iki+nu)*bfn_hess(ipt,mu,hyy)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,1,2,lpert+3) 
     &                  = dgrhodb(ipt,1,2,lpert+3)
     &                  - adens(iki+nu)*bfn_hess(ipt,mu,hyz)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,1,3,lpert+1) 
     &                  = dgrhodb(ipt,1,3,lpert+1)
     &                  - adens(iki+nu)*bfn_hess(ipt,mu,hxz)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,1,3,lpert+2) 
     &                  = dgrhodb(ipt,1,3,lpert+2)
     &                  - adens(iki+nu)*bfn_hess(ipt,mu,hyz)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,1,3,lpert+3) 
     &                  = dgrhodb(ipt,1,3,lpert+3)
     &                  - adens(iki+nu)*bfn_hess(ipt,mu,hzz)*
     &                    bfn_val(ipt,nu)
                     enddo ! ipt
                  enddo ! nu
                  do nu=mu+1, nao
                     iki=iky(nu)
                     do ipt=1,npts
                        do l=1,3
                           do k=1,3
                              dgrhodb(ipt,1,k,lpert+l) 
     &                        = dgrhodb(ipt,1,k,lpert+l)
     &                        - adens(iki+mu)*bfng_val(ipt,mu,l)*
     &                          bfng_val(ipt,nu,k)
                           enddo ! k
                        enddo ! l
c
                        dgrhodb(ipt,1,1,lpert+1) 
     &                  = dgrhodb(ipt,1,1,lpert+1)
     &                  - adens(iki+mu)*bfn_hess(ipt,mu,hxx)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,1,1,lpert+2) 
     &                  = dgrhodb(ipt,1,1,lpert+2)
     &                  - adens(iki+mu)*bfn_hess(ipt,mu,hxy)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,1,1,lpert+3) 
     &                  = dgrhodb(ipt,1,1,lpert+3)
     &                  - adens(iki+mu)*bfn_hess(ipt,mu,hxz)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,1,2,lpert+1) 
     &                  = dgrhodb(ipt,1,2,lpert+1)
     &                  - adens(iki+mu)*bfn_hess(ipt,mu,hxy)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,1,2,lpert+2) 
     &                  = dgrhodb(ipt,1,2,lpert+2)
     &                  - adens(iki+mu)*bfn_hess(ipt,mu,hyy)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,1,2,lpert+3) 
     &                  = dgrhodb(ipt,1,2,lpert+3)
     &                  - adens(iki+mu)*bfn_hess(ipt,mu,hyz)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,1,3,lpert+1) 
     &                  = dgrhodb(ipt,1,3,lpert+1)
     &                  - adens(iki+mu)*bfn_hess(ipt,mu,hxz)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,1,3,lpert+2) 
     &                  = dgrhodb(ipt,1,3,lpert+2)
     &                  - adens(iki+mu)*bfn_hess(ipt,mu,hyz)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,1,3,lpert+3) 
     &                  = dgrhodb(ipt,1,3,lpert+3)
     &                  - adens(iki+mu)*bfn_hess(ipt,mu,hzz)*
     &                    bfn_val(ipt,nu)
                     enddo ! ipt
                  enddo ! nu
               enddo ! mu
 120        continue ! latm
         endif ! gradcorr_sw
      else ! rkstyp_sw
         bend = 0
         do 130 latm=1,ngridcentres
            lpert=3*(latm-1)
            atmt = atom_tag(ao_tag,latm)
            bstart=bend+1
            bend=bend+Abfn(ao_tag,atmt)
c           if (gwt_sw.and.latm.eq.iatm) goto 130
cDEBUG
c        write(*,*)'latm  =',latm
c        write(*,*)'atmt  =',atmt
c        write(*,*)'bstart=',bstart
c        write(*,*)'bend  =',bend
cDEBUG
            do mu=bstart,bend
               iki=iky(mu)
               do nu=1, mu
                  do ipt=1,npts
                     drhodb(ipt,1,lpert+1) = drhodb(ipt,1,lpert+1)
     &               - 2*adens(iki+nu)*bfng_val(ipt,mu,1)*
     &                 bfn_val(ipt,nu)
                     drhodb(ipt,1,lpert+2) = drhodb(ipt,1,lpert+2)
     &               - 2*adens(iki+nu)*bfng_val(ipt,mu,2)*
     &                 bfn_val(ipt,nu)
                     drhodb(ipt,1,lpert+3) = drhodb(ipt,1,lpert+3)
     &               - 2*adens(iki+nu)*bfng_val(ipt,mu,3)*
     &                 bfn_val(ipt,nu)
                  enddo ! ipt
               enddo ! nu
               do nu=mu+1, nao
                  iki=iky(nu)
                  do ipt=1,npts
                     drhodb(ipt,1,lpert+1) = drhodb(ipt,1,lpert+1)
     &               - 2*adens(iki+mu)*bfng_val(ipt,mu,1)*
     &                 bfn_val(ipt,nu)
                     drhodb(ipt,1,lpert+2) = drhodb(ipt,1,lpert+2)
     &               - 2*adens(iki+mu)*bfng_val(ipt,mu,2)*
     &                 bfn_val(ipt,nu)
                     drhodb(ipt,1,lpert+3) = drhodb(ipt,1,lpert+3)
     &               - 2*adens(iki+mu)*bfng_val(ipt,mu,3)*
     &                 bfn_val(ipt,nu)
                  enddo ! ipt
               enddo ! nu
            enddo ! mu
c
            do mu=bstart,bend
               iki=iky(mu)
               do nu=1, mu
                  do ipt=1,npts
                     drhodb(ipt,2,lpert+1) = drhodb(ipt,2,lpert+1)
     &               - 2*bdens(iki+nu)*bfng_val(ipt,mu,1)*
     &                 bfn_val(ipt,nu)
                     drhodb(ipt,2,lpert+2) = drhodb(ipt,2,lpert+2)
     &               - 2*bdens(iki+nu)*bfng_val(ipt,mu,2)*
     &                 bfn_val(ipt,nu)
                     drhodb(ipt,2,lpert+3) = drhodb(ipt,2,lpert+3)
     &               - 2*bdens(iki+nu)*bfng_val(ipt,mu,3)*
     &                 bfn_val(ipt,nu)
                  enddo ! ipt
               enddo ! nu
            enddo ! mu
            do mu=bstart,bend
               do nu=mu+1, nao
                  iki=iky(nu)
                  do ipt=1,npts
                     drhodb(ipt,2,lpert+1) = drhodb(ipt,2,lpert+1)
     &               - 2*bdens(iki+mu)*bfng_val(ipt,mu,1)*
     &                 bfn_val(ipt,nu)
                     drhodb(ipt,2,lpert+2) = drhodb(ipt,2,lpert+2)
     &               - 2*bdens(iki+mu)*bfng_val(ipt,mu,2)*
     &                 bfn_val(ipt,nu)
                     drhodb(ipt,2,lpert+3) = drhodb(ipt,2,lpert+3)
     &               - 2*bdens(iki+mu)*bfng_val(ipt,mu,3)*
     &                 bfn_val(ipt,nu)
                  enddo ! ipt
               enddo ! nu
            enddo  ! mu
 130     continue ! latm
         if (gradcorr_sw) then
            bend = 0
            do 140 latm=1,ngridcentres
               lpert=3*(latm-1)
               atmt = atom_tag(ao_tag,latm)
               bstart=bend+1
               bend=bend+Abfn(ao_tag,atmt)
c              if (gwt_sw.and.latm.eq.iatm) goto 140
cDEBUG
c        write(*,*)'latm  =',latm
c        write(*,*)'atmt  =',atmt
c        write(*,*)'bstart=',bstart
c        write(*,*)'bend  =',bend
cDEBUG
               do mu=bstart,bend
                  iki=iky(mu)
                  do nu=1, mu
                     do ipt=1,npts
                        dgrhodb(ipt,1,1,lpert+1) 
     &                  = dgrhodb(ipt,1,1,lpert+1)
     &                  - 2*adens(iki+nu)*bfng_val(ipt,mu,1)*
     &                    bfng_val(ipt,nu,1)
                        dgrhodb(ipt,1,1,lpert+2) 
     &                  = dgrhodb(ipt,1,1,lpert+2)
     &                  - 2*adens(iki+nu)*bfng_val(ipt,mu,2)*
     &                    bfng_val(ipt,nu,1)
                        dgrhodb(ipt,1,1,lpert+3) 
     &                  = dgrhodb(ipt,1,1,lpert+3)
     &                  - 2*adens(iki+nu)*bfng_val(ipt,mu,3)*
     &                    bfng_val(ipt,nu,1)
                        dgrhodb(ipt,1,2,lpert+1) 
     &                  = dgrhodb(ipt,1,2,lpert+1)
     &                  - 2*adens(iki+nu)*bfng_val(ipt,mu,1)*
     &                    bfng_val(ipt,nu,2)
                        dgrhodb(ipt,1,2,lpert+2) 
     &                  = dgrhodb(ipt,1,2,lpert+2)
     &                  - 2*adens(iki+nu)*bfng_val(ipt,mu,2)*
     &                    bfng_val(ipt,nu,2)
                        dgrhodb(ipt,1,2,lpert+3) 
     &                  = dgrhodb(ipt,1,2,lpert+3)
     &                  - 2*adens(iki+nu)*bfng_val(ipt,mu,3)*
     &                    bfng_val(ipt,nu,2)
                        dgrhodb(ipt,1,3,lpert+1) 
     &                  = dgrhodb(ipt,1,3,lpert+1)
     &                  - 2*adens(iki+nu)*bfng_val(ipt,mu,1)*
     &                    bfng_val(ipt,nu,3)
                        dgrhodb(ipt,1,3,lpert+2) 
     &                  = dgrhodb(ipt,1,3,lpert+2)
     &                  - 2*adens(iki+nu)*bfng_val(ipt,mu,2)*
     &                    bfng_val(ipt,nu,3)
                        dgrhodb(ipt,1,3,lpert+3) 
     &                  = dgrhodb(ipt,1,3,lpert+3)
     &                  - 2*adens(iki+nu)*bfng_val(ipt,mu,3)*
     &                    bfng_val(ipt,nu,3)
c
                        dgrhodb(ipt,1,1,lpert+1) 
     &                  = dgrhodb(ipt,1,1,lpert+1)
     &                  - 2*adens(iki+nu)*bfn_hess(ipt,mu,hxx)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,1,1,lpert+2) 
     &                  = dgrhodb(ipt,1,1,lpert+2)
     &                  - 2*adens(iki+nu)*bfn_hess(ipt,mu,hxy)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,1,1,lpert+3) 
     &                  = dgrhodb(ipt,1,1,lpert+3)
     &                  - 2*adens(iki+nu)*bfn_hess(ipt,mu,hxz)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,1,2,lpert+1) 
     &                  = dgrhodb(ipt,1,2,lpert+1)
     &                  - 2*adens(iki+nu)*bfn_hess(ipt,mu,hxy)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,1,2,lpert+2) 
     &                  = dgrhodb(ipt,1,2,lpert+2)
     &                  - 2*adens(iki+nu)*bfn_hess(ipt,mu,hyy)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,1,2,lpert+3) 
     &                  = dgrhodb(ipt,1,2,lpert+3)
     &                  - 2*adens(iki+nu)*bfn_hess(ipt,mu,hyz)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,1,3,lpert+1) 
     &                  = dgrhodb(ipt,1,3,lpert+1)
     &                  - 2*adens(iki+nu)*bfn_hess(ipt,mu,hxz)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,1,3,lpert+2) 
     &                  = dgrhodb(ipt,1,3,lpert+2)
     &                  - 2*adens(iki+nu)*bfn_hess(ipt,mu,hyz)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,1,3,lpert+3) 
     &                  = dgrhodb(ipt,1,3,lpert+3)
     &                  - 2*adens(iki+nu)*bfn_hess(ipt,mu,hzz)*
     &                    bfn_val(ipt,nu)
                     enddo
                  enddo
                  do nu=mu+1, nao
                     iki=iky(nu)
                     do ipt=1,npts
                        dgrhodb(ipt,1,1,lpert+1) 
     &                  = dgrhodb(ipt,1,1,lpert+1)
     &                  - 2*adens(iki+mu)*bfng_val(ipt,mu,1)*
     &                    bfng_val(ipt,nu,1)
                        dgrhodb(ipt,1,1,lpert+2) 
     &                  = dgrhodb(ipt,1,1,lpert+2)
     &                  - 2*adens(iki+mu)*bfng_val(ipt,mu,2)*
     &                    bfng_val(ipt,nu,1)
                        dgrhodb(ipt,1,1,lpert+3) 
     &                  = dgrhodb(ipt,1,1,lpert+3)
     &                  - 2*adens(iki+mu)*bfng_val(ipt,mu,3)*
     &                    bfng_val(ipt,nu,1)
                        dgrhodb(ipt,1,2,lpert+1) 
     &                  = dgrhodb(ipt,1,2,lpert+1)
     &                  - 2*adens(iki+mu)*bfng_val(ipt,mu,1)*
     &                    bfng_val(ipt,nu,2)
                        dgrhodb(ipt,1,2,lpert+2) 
     &                  = dgrhodb(ipt,1,2,lpert+2)
     &                  - 2*adens(iki+mu)*bfng_val(ipt,mu,2)*
     &                    bfng_val(ipt,nu,2)
                        dgrhodb(ipt,1,2,lpert+3) 
     &                  = dgrhodb(ipt,1,2,lpert+3)
     &                  - 2*adens(iki+mu)*bfng_val(ipt,mu,3)*
     &                    bfng_val(ipt,nu,2)
                        dgrhodb(ipt,1,3,lpert+1) 
     &                  = dgrhodb(ipt,1,3,lpert+1)
     &                  - 2*adens(iki+mu)*bfng_val(ipt,mu,1)*
     &                    bfng_val(ipt,nu,3)
                        dgrhodb(ipt,1,3,lpert+2) 
     &                  = dgrhodb(ipt,1,3,lpert+2)
     &                  - 2*adens(iki+mu)*bfng_val(ipt,mu,2)*
     &                    bfng_val(ipt,nu,3)
                        dgrhodb(ipt,1,3,lpert+3) 
     &                  = dgrhodb(ipt,1,3,lpert+3)
     &                  - 2*adens(iki+mu)*bfng_val(ipt,mu,3)*
     &                    bfng_val(ipt,nu,3)
c
                        dgrhodb(ipt,1,1,lpert+1) 
     &                  = dgrhodb(ipt,1,1,lpert+1)
     &                  - 2*adens(iki+mu)*bfn_hess(ipt,mu,hxx)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,1,1,lpert+2) 
     &                  = dgrhodb(ipt,1,1,lpert+2)
     &                  - 2*adens(iki+mu)*bfn_hess(ipt,mu,hxy)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,1,1,lpert+3) 
     &                  = dgrhodb(ipt,1,1,lpert+3)
     &                  - 2*adens(iki+mu)*bfn_hess(ipt,mu,hxz)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,1,2,lpert+1) 
     &                  = dgrhodb(ipt,1,2,lpert+1)
     &                  - 2*adens(iki+mu)*bfn_hess(ipt,mu,hxy)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,1,2,lpert+2) 
     &                  = dgrhodb(ipt,1,2,lpert+2)
     &                  - 2*adens(iki+mu)*bfn_hess(ipt,mu,hyy)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,1,2,lpert+3) 
     &                  = dgrhodb(ipt,1,2,lpert+3)
     &                  - 2*adens(iki+mu)*bfn_hess(ipt,mu,hyz)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,1,3,lpert+1) 
     &                  = dgrhodb(ipt,1,3,lpert+1)
     &                  - 2*adens(iki+mu)*bfn_hess(ipt,mu,hxz)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,1,3,lpert+2) 
     &                  = dgrhodb(ipt,1,3,lpert+2)
     &                  - 2*adens(iki+mu)*bfn_hess(ipt,mu,hyz)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,1,3,lpert+3) 
     &                  = dgrhodb(ipt,1,3,lpert+3)
     &                  - 2*adens(iki+mu)*bfn_hess(ipt,mu,hzz)*
     &                    bfn_val(ipt,nu)
                     enddo
                  enddo
               enddo
c
               do mu=bstart,bend
                  iki=iky(mu)
                  do nu=1, mu
                     do ipt=1,npts
                        dgrhodb(ipt,2,1,lpert+1) 
     &                  = dgrhodb(ipt,2,1,lpert+1)
     &                  - 2*bdens(iki+nu)*bfng_val(ipt,mu,1)*
     &                    bfng_val(ipt,nu,1)
                        dgrhodb(ipt,2,1,lpert+2) 
     &                  = dgrhodb(ipt,2,1,lpert+2)
     &                  - 2*bdens(iki+nu)*bfng_val(ipt,mu,2)*
     &                    bfng_val(ipt,nu,1)
                        dgrhodb(ipt,2,1,lpert+3) 
     &                  = dgrhodb(ipt,2,1,lpert+3)
     &                  - 2*bdens(iki+nu)*bfng_val(ipt,mu,3)*
     &                    bfng_val(ipt,nu,1)
                        dgrhodb(ipt,2,2,lpert+1) 
     &                  = dgrhodb(ipt,2,2,lpert+1)
     &                  - 2*bdens(iki+nu)*bfng_val(ipt,mu,1)*
     &                    bfng_val(ipt,nu,2)
                        dgrhodb(ipt,2,2,lpert+2) 
     &                  = dgrhodb(ipt,2,2,lpert+2)
     &                  - 2*bdens(iki+nu)*bfng_val(ipt,mu,2)*
     &                    bfng_val(ipt,nu,2)
                        dgrhodb(ipt,2,2,lpert+3) 
     &                  = dgrhodb(ipt,2,2,lpert+3)
     &                  - 2*bdens(iki+nu)*bfng_val(ipt,mu,3)*
     &                    bfng_val(ipt,nu,2)
                        dgrhodb(ipt,2,3,lpert+1) 
     &                  = dgrhodb(ipt,2,3,lpert+1)
     &                  - 2*bdens(iki+nu)*bfng_val(ipt,mu,1)*
     &                    bfng_val(ipt,nu,3)
                        dgrhodb(ipt,2,3,lpert+2) 
     &                  = dgrhodb(ipt,2,3,lpert+2)
     &                  - 2*bdens(iki+nu)*bfng_val(ipt,mu,2)*
     &                    bfng_val(ipt,nu,3)
                        dgrhodb(ipt,2,3,lpert+3) 
     &                  = dgrhodb(ipt,2,3,lpert+3)
     &                  - 2*bdens(iki+nu)*bfng_val(ipt,mu,3)*
     &                    bfng_val(ipt,nu,3)
c
                        dgrhodb(ipt,2,1,lpert+1) 
     &                  = dgrhodb(ipt,2,1,lpert+1)
     &                  - 2*bdens(iki+nu)*bfn_hess(ipt,mu,hxx)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,2,1,lpert+2) 
     &                  = dgrhodb(ipt,2,1,lpert+2)
     &                  - 2*bdens(iki+nu)*bfn_hess(ipt,mu,hxy)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,2,1,lpert+3) 
     &                  = dgrhodb(ipt,2,1,lpert+3)
     &                  - 2*bdens(iki+nu)*bfn_hess(ipt,mu,hxz)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,2,2,lpert+1) 
     &                  = dgrhodb(ipt,2,2,lpert+1)
     &                  - 2*bdens(iki+nu)*bfn_hess(ipt,mu,hxy)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,2,2,lpert+2) 
     &                  = dgrhodb(ipt,2,2,lpert+2)
     &                  - 2*bdens(iki+nu)*bfn_hess(ipt,mu,hyy)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,2,2,lpert+3) 
     &                  = dgrhodb(ipt,2,2,lpert+3)
     &                  - 2*bdens(iki+nu)*bfn_hess(ipt,mu,hyz)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,2,3,lpert+1) 
     &                  = dgrhodb(ipt,2,3,lpert+1)
     &                  - 2*bdens(iki+nu)*bfn_hess(ipt,mu,hxz)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,2,3,lpert+2) 
     &                  = dgrhodb(ipt,2,3,lpert+2)
     &                  - 2*bdens(iki+nu)*bfn_hess(ipt,mu,hyz)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,2,3,lpert+3) 
     &                  = dgrhodb(ipt,2,3,lpert+3)
     &                  - 2*bdens(iki+nu)*bfn_hess(ipt,mu,hzz)*
     &                    bfn_val(ipt,nu)
                     enddo
                  enddo
                  do nu=mu+1, nao
                     iki=iky(nu)
                     do ipt=1,npts
                        dgrhodb(ipt,2,1,lpert+1) 
     &                  = dgrhodb(ipt,2,1,lpert+1)
     &                  - 2*bdens(iki+mu)*bfng_val(ipt,mu,1)*
     &                    bfng_val(ipt,nu,1)
                        dgrhodb(ipt,2,1,lpert+2) 
     &                  = dgrhodb(ipt,2,1,lpert+2)
     &                  - 2*bdens(iki+mu)*bfng_val(ipt,mu,2)*
     &                    bfng_val(ipt,nu,1)
                        dgrhodb(ipt,2,1,lpert+3) 
     &                  = dgrhodb(ipt,2,1,lpert+3)
     &                  - 2*bdens(iki+mu)*bfng_val(ipt,mu,3)*
     &                    bfng_val(ipt,nu,1)
                        dgrhodb(ipt,2,2,lpert+1) 
     &                  = dgrhodb(ipt,2,2,lpert+1)
     &                  - 2*bdens(iki+mu)*bfng_val(ipt,mu,1)*
     &                    bfng_val(ipt,nu,2)
                        dgrhodb(ipt,2,2,lpert+2) 
     &                  = dgrhodb(ipt,2,2,lpert+2)
     &                  - 2*bdens(iki+mu)*bfng_val(ipt,mu,2)*
     &                    bfng_val(ipt,nu,2)
                        dgrhodb(ipt,2,2,lpert+3) 
     &                  = dgrhodb(ipt,2,2,lpert+3)
     &                  - 2*bdens(iki+mu)*bfng_val(ipt,mu,3)*
     &                    bfng_val(ipt,nu,2)
                        dgrhodb(ipt,2,3,lpert+1) 
     &                  = dgrhodb(ipt,2,3,lpert+1)
     &                  - 2*bdens(iki+mu)*bfng_val(ipt,mu,1)*
     &                    bfng_val(ipt,nu,3)
                        dgrhodb(ipt,2,3,lpert+2) 
     &                  = dgrhodb(ipt,2,3,lpert+2)
     &                  - 2*bdens(iki+mu)*bfng_val(ipt,mu,2)*
     &                    bfng_val(ipt,nu,3)
                        dgrhodb(ipt,2,3,lpert+3) 
     &                  = dgrhodb(ipt,2,3,lpert+3)
     &                  - 2*bdens(iki+mu)*bfng_val(ipt,mu,3)*
     &                    bfng_val(ipt,nu,3)
c
                        dgrhodb(ipt,2,1,lpert+1) 
     &                  = dgrhodb(ipt,2,1,lpert+1)
     &                  - 2*bdens(iki+mu)*bfn_hess(ipt,mu,hxx)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,2,1,lpert+2) 
     &                  = dgrhodb(ipt,2,1,lpert+2)
     &                  - 2*bdens(iki+mu)*bfn_hess(ipt,mu,hxy)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,2,1,lpert+3) 
     &                  = dgrhodb(ipt,2,1,lpert+3)
     &                  - 2*bdens(iki+mu)*bfn_hess(ipt,mu,hxz)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,2,2,lpert+1) 
     &                  = dgrhodb(ipt,2,2,lpert+1)
     &                  - 2*bdens(iki+mu)*bfn_hess(ipt,mu,hxy)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,2,2,lpert+2) 
     &                  = dgrhodb(ipt,2,2,lpert+2)
     &                  - 2*bdens(iki+mu)*bfn_hess(ipt,mu,hyy)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,2,2,lpert+3) 
     &                  = dgrhodb(ipt,2,2,lpert+3)
     &                  - 2*bdens(iki+mu)*bfn_hess(ipt,mu,hyz)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,2,3,lpert+1) 
     &                  = dgrhodb(ipt,2,3,lpert+1)
     &                  - 2*bdens(iki+mu)*bfn_hess(ipt,mu,hxz)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,2,3,lpert+2) 
     &                  = dgrhodb(ipt,2,3,lpert+2)
     &                  - 2*bdens(iki+mu)*bfn_hess(ipt,mu,hyz)*
     &                    bfn_val(ipt,nu)
                        dgrhodb(ipt,2,3,lpert+3) 
     &                  = dgrhodb(ipt,2,3,lpert+3)
     &                  - 2*bdens(iki+mu)*bfn_hess(ipt,mu,hzz)*
     &                    bfn_val(ipt,nu)
                     enddo
                  enddo
               enddo
 140        continue
         endif ! gradcorr_sw
      endif ! rkstype_sw

c     do 100 latm=1,ngridcentres
c        lpert=3*(latm-1)
c        atmt = atom_tag(ao_tag,latm)
c        bstart=bend+1
c        bend=bend+Abfn(ao_tag,atmt)
c        if (gwt_sw.and.latm.eq.iatm) goto 100
cDEBUG
c        write(*,*)'latm  =',latm
c        write(*,*)'atmt  =',atmt
c        write(*,*)'bstart=',bstart
c        write(*,*)'bend  =',bend
cDEBUG
c        do mu=bstart,bend
c           iki=iky(mu)
c           do nu=1, mu
c              do ipt=1,npts
c                 drhodb(ipt,1,lpert+1) = drhodb(ipt,1,lpert+1)
c    &            - 2*adens(iki+nu)*bfng_val(ipt,mu,1)*
c    &              bfn_val(ipt,nu)
c                 drhodb(ipt,1,lpert+2) = drhodb(ipt,1,lpert+2)
c    &            - 2*adens(iki+nu)*bfng_val(ipt,mu,2)*
c    &              bfn_val(ipt,nu)
c                 drhodb(ipt,1,lpert+3) = drhodb(ipt,1,lpert+3)
c    &            - 2*adens(iki+nu)*bfng_val(ipt,mu,3)*
c    &              bfn_val(ipt,nu)
c              enddo
c           enddo
c           do nu=mu+1, nao
c              iki=iky(nu)
c              do ipt=1,npts
c                 drhodb(ipt,1,lpert+1) = drhodb(ipt,1,lpert+1)
c    &            - 2*adens(iki+mu)*bfng_val(ipt,mu,1)*
c    &              bfn_val(ipt,nu)
c                 drhodb(ipt,1,lpert+2) = drhodb(ipt,1,lpert+2)
c    &            - 2*adens(iki+mu)*bfng_val(ipt,mu,2)*
c    &              bfn_val(ipt,nu)
c                 drhodb(ipt,1,lpert+3) = drhodb(ipt,1,lpert+3)
c    &            - 2*adens(iki+mu)*bfng_val(ipt,mu,3)*
c    &              bfn_val(ipt,nu)
c              enddo
c           enddo
c        enddo
c
c        if (gradcorr_sw) then
c           do mu=bstart,bend
c              iki=iky(mu)
c              do nu=1, mu
c                 do ipt=1,npts
c                    dgrhodb(ipt,1,1,lpert+1) = dgrhodb(ipt,1,1,lpert+1)
c    &               - 2*adens(iki+nu)*bfng_val(ipt,mu,1)*
c    &                 bfng_val(ipt,nu,1)
c                    dgrhodb(ipt,1,1,lpert+2) = dgrhodb(ipt,1,1,lpert+2)
c    &               - 2*adens(iki+nu)*bfng_val(ipt,mu,2)*
c    &                 bfng_val(ipt,nu,1)
c                    dgrhodb(ipt,1,1,lpert+3) = dgrhodb(ipt,1,1,lpert+3)
c    &               - 2*adens(iki+nu)*bfng_val(ipt,mu,3)*
c    &                 bfng_val(ipt,nu,1)
c                    dgrhodb(ipt,1,2,lpert+1) = dgrhodb(ipt,1,2,lpert+1)
c    &               - 2*adens(iki+nu)*bfng_val(ipt,mu,1)*
c    &                 bfng_val(ipt,nu,2)
c                    dgrhodb(ipt,1,2,lpert+2) = dgrhodb(ipt,1,2,lpert+2)
c    &               - 2*adens(iki+nu)*bfng_val(ipt,mu,2)*
c    &                 bfng_val(ipt,nu,2)
c                    dgrhodb(ipt,1,2,lpert+3) = dgrhodb(ipt,1,2,lpert+3)
c    &               - 2*adens(iki+nu)*bfng_val(ipt,mu,3)*
c    &                 bfng_val(ipt,nu,2)
c                    dgrhodb(ipt,1,3,lpert+1) = dgrhodb(ipt,1,3,lpert+1)
c    &               - 2*adens(iki+nu)*bfng_val(ipt,mu,1)*
c    &                 bfng_val(ipt,nu,3)
c                    dgrhodb(ipt,1,3,lpert+2) = dgrhodb(ipt,1,3,lpert+2)
c    &               - 2*adens(iki+nu)*bfng_val(ipt,mu,2)*
c    &                 bfng_val(ipt,nu,3)
c                    dgrhodb(ipt,1,3,lpert+3) = dgrhodb(ipt,1,3,lpert+3)
c    &               - 2*adens(iki+nu)*bfng_val(ipt,mu,3)*
c    &                 bfng_val(ipt,nu,3)
c
c                    dgrhodb(ipt,1,1,lpert+1) = dgrhodb(ipt,1,1,lpert+1)
c    &               - 2*adens(iki+nu)*bfn_hess(ipt,mu,hxx)*
c    &                 bfn_val(ipt,nu)
c                    dgrhodb(ipt,1,1,lpert+2) = dgrhodb(ipt,1,1,lpert+2)
c    &               - 2*adens(iki+nu)*bfn_hess(ipt,mu,hxy)*
c    &                 bfn_val(ipt,nu)
c                    dgrhodb(ipt,1,1,lpert+3) = dgrhodb(ipt,1,1,lpert+3)
c    &               - 2*adens(iki+nu)*bfn_hess(ipt,mu,hxz)*
c    &                 bfn_val(ipt,nu)
c                    dgrhodb(ipt,1,2,lpert+1) = dgrhodb(ipt,1,2,lpert+1)
c    &               - 2*adens(iki+nu)*bfn_hess(ipt,mu,hxy)*
c    &                 bfn_val(ipt,nu)
c                    dgrhodb(ipt,1,2,lpert+2) = dgrhodb(ipt,1,2,lpert+2)
c    &               - 2*adens(iki+nu)*bfn_hess(ipt,mu,hyy)*
c    &                 bfn_val(ipt,nu)
c                    dgrhodb(ipt,1,2,lpert+3) = dgrhodb(ipt,1,2,lpert+3)
c    &               - 2*adens(iki+nu)*bfn_hess(ipt,mu,hyz)*
c    &                 bfn_val(ipt,nu)
c                    dgrhodb(ipt,1,3,lpert+1) = dgrhodb(ipt,1,3,lpert+1)
c    &               - 2*adens(iki+nu)*bfn_hess(ipt,mu,hxz)*
c    &                 bfn_val(ipt,nu)
c                    dgrhodb(ipt,1,3,lpert+2) = dgrhodb(ipt,1,3,lpert+2)
c    &               - 2*adens(iki+nu)*bfn_hess(ipt,mu,hyz)*
c    &                 bfn_val(ipt,nu)
c                    dgrhodb(ipt,1,3,lpert+3) = dgrhodb(ipt,1,3,lpert+3)
c    &               - 2*adens(iki+nu)*bfn_hess(ipt,mu,hzz)*
c    &                 bfn_val(ipt,nu)
c                 enddo
c              enddo
c              do nu=mu+1, nao
c                 iki=iky(nu)
c                 do ipt=1,npts
c                    dgrhodb(ipt,1,1,lpert+1) = dgrhodb(ipt,1,1,lpert+1)
c    &               - 2*adens(iki+mu)*bfng_val(ipt,mu,1)*
c    &                 bfng_val(ipt,nu,1)
c                    dgrhodb(ipt,1,1,lpert+2) = dgrhodb(ipt,1,1,lpert+2)
c    &               - 2*adens(iki+mu)*bfng_val(ipt,mu,2)*
c    &                 bfng_val(ipt,nu,1)
c                    dgrhodb(ipt,1,1,lpert+3) = dgrhodb(ipt,1,1,lpert+3)
c    &               - 2*adens(iki+mu)*bfng_val(ipt,mu,3)*
c    &                 bfng_val(ipt,nu,1)
c                    dgrhodb(ipt,1,2,lpert+1) = dgrhodb(ipt,1,2,lpert+1)
c    &               - 2*adens(iki+mu)*bfng_val(ipt,mu,1)*
c    &                 bfng_val(ipt,nu,2)
c                    dgrhodb(ipt,1,2,lpert+2) = dgrhodb(ipt,1,2,lpert+2)
c    &               - 2*adens(iki+mu)*bfng_val(ipt,mu,2)*
c    &                 bfng_val(ipt,nu,2)
c                    dgrhodb(ipt,1,2,lpert+3) = dgrhodb(ipt,1,2,lpert+3)
c    &               - 2*adens(iki+mu)*bfng_val(ipt,mu,3)*
c    &                 bfng_val(ipt,nu,2)
c                    dgrhodb(ipt,1,3,lpert+1) = dgrhodb(ipt,1,3,lpert+1)
c    &               - 2*adens(iki+mu)*bfng_val(ipt,mu,1)*
c    &                 bfng_val(ipt,nu,3)
c                    dgrhodb(ipt,1,3,lpert+2) = dgrhodb(ipt,1,3,lpert+2)
c    &               - 2*adens(iki+mu)*bfng_val(ipt,mu,2)*
c    &                 bfng_val(ipt,nu,3)
c                    dgrhodb(ipt,1,3,lpert+3) = dgrhodb(ipt,1,3,lpert+3)
c    &               - 2*adens(iki+mu)*bfng_val(ipt,mu,3)*
c    &                 bfng_val(ipt,nu,3)
c
c                    dgrhodb(ipt,1,1,lpert+1) = dgrhodb(ipt,1,1,lpert+1)
c    &               - 2*adens(iki+mu)*bfn_hess(ipt,mu,hxx)*
c    &                 bfn_val(ipt,nu)
c                    dgrhodb(ipt,1,1,lpert+2) = dgrhodb(ipt,1,1,lpert+2)
c    &               - 2*adens(iki+mu)*bfn_hess(ipt,mu,hxy)*
c    &                 bfn_val(ipt,nu)
c                    dgrhodb(ipt,1,1,lpert+3) = dgrhodb(ipt,1,1,lpert+3)
c    &               - 2*adens(iki+mu)*bfn_hess(ipt,mu,hxz)*
c    &                 bfn_val(ipt,nu)
c                    dgrhodb(ipt,1,2,lpert+1) = dgrhodb(ipt,1,2,lpert+1)
c    &               - 2*adens(iki+mu)*bfn_hess(ipt,mu,hxy)*
c    &                 bfn_val(ipt,nu)
c                    dgrhodb(ipt,1,2,lpert+2) = dgrhodb(ipt,1,2,lpert+2)
c    &               - 2*adens(iki+mu)*bfn_hess(ipt,mu,hyy)*
c    &                 bfn_val(ipt,nu)
c                    dgrhodb(ipt,1,2,lpert+3) = dgrhodb(ipt,1,2,lpert+3)
c    &               - 2*adens(iki+mu)*bfn_hess(ipt,mu,hyz)*
c    &                 bfn_val(ipt,nu)
c                    dgrhodb(ipt,1,3,lpert+1) = dgrhodb(ipt,1,3,lpert+1)
c    &               - 2*adens(iki+mu)*bfn_hess(ipt,mu,hxz)*
c    &                 bfn_val(ipt,nu)
c                    dgrhodb(ipt,1,3,lpert+2) = dgrhodb(ipt,1,3,lpert+2)
c    &               - 2*adens(iki+mu)*bfn_hess(ipt,mu,hyz)*
c    &                 bfn_val(ipt,nu)
c                    dgrhodb(ipt,1,3,lpert+3) = dgrhodb(ipt,1,3,lpert+3)
c    &               - 2*adens(iki+mu)*bfn_hess(ipt,mu,hzz)*
c    &                 bfn_val(ipt,nu)
c                 enddo
c              enddo
c           enddo
c        endif
c
c        if (.not.rkstyp_sw) then
c           do mu=bstart,bend
c              iki=iky(mu)
c              do nu=1, mu
c                 do ipt=1,npts
c                    drhodb(ipt,2,lpert+1) = drhodb(ipt,2,lpert+1)
c    &               - 2*bdens(iki+nu)*bfng_val(ipt,mu,1)*
c    &                 bfn_val(ipt,nu)
c                    drhodb(ipt,2,lpert+2) = drhodb(ipt,2,lpert+2)
c    &               - 2*bdens(iki+nu)*bfng_val(ipt,mu,2)*
c    &                 bfn_val(ipt,nu)
c                    drhodb(ipt,2,lpert+3) = drhodb(ipt,2,lpert+3)
c    &               - 2*bdens(iki+nu)*bfng_val(ipt,mu,3)*
c    &                 bfn_val(ipt,nu)
c                 enddo
c              enddo
c           enddo
c           do mu=bstart,bend
c              do nu=mu+1, nao
c                 iki=iky(nu)
c                 do ipt=1,npts
c                    drhodb(ipt,2,lpert+1) = drhodb(ipt,2,lpert+1)
c    &               - 2*bdens(iki+mu)*bfng_val(ipt,mu,1)*
c    &                 bfn_val(ipt,nu)
c                    drhodb(ipt,2,lpert+2) = drhodb(ipt,2,lpert+2)
c    &               - 2*bdens(iki+mu)*bfng_val(ipt,mu,2)*
c    &                 bfn_val(ipt,nu)
c                    drhodb(ipt,2,lpert+3) = drhodb(ipt,2,lpert+3)
c    &               - 2*bdens(iki+mu)*bfng_val(ipt,mu,3)*
c    &                 bfn_val(ipt,nu)
c                 enddo
c              enddo
c           enddo
c
c           if (gradcorr_sw) then
c              do mu=bstart,bend
c                 iki=iky(mu)
c                 do nu=1, mu
c                    do ipt=1,npts
c                       dgrhodb(ipt,2,1,lpert+1) 
c    &                  = dgrhodb(ipt,2,1,lpert+1)
c    &                  - 2*bdens(iki+nu)*bfng_val(ipt,mu,1)*
c    &                    bfng_val(ipt,nu,1)
c                       dgrhodb(ipt,2,1,lpert+2) 
c    &                  = dgrhodb(ipt,2,1,lpert+2)
c    &                  - 2*bdens(iki+nu)*bfng_val(ipt,mu,2)*
c    &                    bfng_val(ipt,nu,1)
c                       dgrhodb(ipt,2,1,lpert+3) 
c    &                  = dgrhodb(ipt,2,1,lpert+3)
c    &                  - 2*bdens(iki+nu)*bfng_val(ipt,mu,3)*
c    &                    bfng_val(ipt,nu,1)
c                       dgrhodb(ipt,2,2,lpert+1) 
c    &                  = dgrhodb(ipt,2,2,lpert+1)
c    &                  - 2*bdens(iki+nu)*bfng_val(ipt,mu,1)*
c    &                    bfng_val(ipt,nu,2)
c                       dgrhodb(ipt,2,2,lpert+2) 
c    &                  = dgrhodb(ipt,2,2,lpert+2)
c    &                  - 2*bdens(iki+nu)*bfng_val(ipt,mu,2)*
c    &                    bfng_val(ipt,nu,2)
c                       dgrhodb(ipt,2,2,lpert+3) 
c    &                  = dgrhodb(ipt,2,2,lpert+3)
c    &                  - 2*bdens(iki+nu)*bfng_val(ipt,mu,3)*
c    &                    bfng_val(ipt,nu,2)
c                       dgrhodb(ipt,2,3,lpert+1) 
c    &                  = dgrhodb(ipt,2,3,lpert+1)
c    &                  - 2*bdens(iki+nu)*bfng_val(ipt,mu,1)*
c    &                    bfng_val(ipt,nu,3)
c                       dgrhodb(ipt,2,3,lpert+2) 
c    &                  = dgrhodb(ipt,2,3,lpert+2)
c    &                  - 2*bdens(iki+nu)*bfng_val(ipt,mu,2)*
c    &                    bfng_val(ipt,nu,3)
c                       dgrhodb(ipt,2,3,lpert+3) 
c    &                  = dgrhodb(ipt,2,3,lpert+3)
c    &                  - 2*bdens(iki+nu)*bfng_val(ipt,mu,3)*
c    &                    bfng_val(ipt,nu,3)
c
c                       dgrhodb(ipt,2,1,lpert+1) 
c    &                  = dgrhodb(ipt,2,1,lpert+1)
c    &                  - 2*bdens(iki+nu)*bfn_hess(ipt,mu,hxx)*
c    &                    bfn_val(ipt,nu)
c                       dgrhodb(ipt,2,1,lpert+2) 
c    &                  = dgrhodb(ipt,2,1,lpert+2)
c    &                  - 2*bdens(iki+nu)*bfn_hess(ipt,mu,hxy)*
c    &                    bfn_val(ipt,nu)
c                       dgrhodb(ipt,2,1,lpert+3) 
c    &                  = dgrhodb(ipt,2,1,lpert+3)
c    &                  - 2*bdens(iki+nu)*bfn_hess(ipt,mu,hxz)*
c    &                    bfn_val(ipt,nu)
c                       dgrhodb(ipt,2,2,lpert+1) 
c    &                  = dgrhodb(ipt,2,2,lpert+1)
c    &                  - 2*bdens(iki+nu)*bfn_hess(ipt,mu,hxy)*
c    &                    bfn_val(ipt,nu)
c                       dgrhodb(ipt,2,2,lpert+2) 
c    &                  = dgrhodb(ipt,2,2,lpert+2)
c    &                  - 2*bdens(iki+nu)*bfn_hess(ipt,mu,hyy)*
c    &                    bfn_val(ipt,nu)
c                       dgrhodb(ipt,2,2,lpert+3) 
c    &                  = dgrhodb(ipt,2,2,lpert+3)
c    &                  - 2*bdens(iki+nu)*bfn_hess(ipt,mu,hyz)*
c    &                    bfn_val(ipt,nu)
c                       dgrhodb(ipt,2,3,lpert+1) 
c    &                  = dgrhodb(ipt,2,3,lpert+1)
c    &                  - 2*bdens(iki+nu)*bfn_hess(ipt,mu,hxz)*
c    &                    bfn_val(ipt,nu)
c                       dgrhodb(ipt,2,3,lpert+2) 
c    &                  = dgrhodb(ipt,2,3,lpert+2)
c    &                  - 2*bdens(iki+nu)*bfn_hess(ipt,mu,hyz)*
c    &                    bfn_val(ipt,nu)
c                       dgrhodb(ipt,2,3,lpert+3) 
c    &                  = dgrhodb(ipt,2,3,lpert+3)
c    &                  - 2*bdens(iki+nu)*bfn_hess(ipt,mu,hzz)*
c    &                    bfn_val(ipt,nu)
c                    enddo
c                 enddo
c                 do nu=mu+1, nao
c                    iki=iky(nu)
c                    do ipt=1,npts
c                       dgrhodb(ipt,2,1,lpert+1) 
c    &                  = dgrhodb(ipt,2,1,lpert+1)
c    &                  - 2*bdens(iki+mu)*bfng_val(ipt,mu,1)*
c    &                    bfng_val(ipt,nu,1)
c                       dgrhodb(ipt,2,1,lpert+2) 
c    &                  = dgrhodb(ipt,2,1,lpert+2)
c    &                  - 2*bdens(iki+mu)*bfng_val(ipt,mu,2)*
c    &                    bfng_val(ipt,nu,1)
c                       dgrhodb(ipt,2,1,lpert+3) 
c    &                  = dgrhodb(ipt,2,1,lpert+3)
c    &                  - 2*bdens(iki+mu)*bfng_val(ipt,mu,3)*
c    &                    bfng_val(ipt,nu,1)
c                       dgrhodb(ipt,2,2,lpert+1) 
c    &                  = dgrhodb(ipt,2,2,lpert+1)
c    &                  - 2*bdens(iki+mu)*bfng_val(ipt,mu,1)*
c    &                    bfng_val(ipt,nu,2)
c                       dgrhodb(ipt,2,2,lpert+2) 
c    &                  = dgrhodb(ipt,2,2,lpert+2)
c    &                  - 2*bdens(iki+mu)*bfng_val(ipt,mu,2)*
c    &                    bfng_val(ipt,nu,2)
c                       dgrhodb(ipt,2,2,lpert+3) 
c    &                  = dgrhodb(ipt,2,2,lpert+3)
c    &                  - 2*bdens(iki+mu)*bfng_val(ipt,mu,3)*
c    &                    bfng_val(ipt,nu,2)
c                       dgrhodb(ipt,2,3,lpert+1) 
c    &                  = dgrhodb(ipt,2,3,lpert+1)
c    &                  - 2*bdens(iki+mu)*bfng_val(ipt,mu,1)*
c    &                    bfng_val(ipt,nu,3)
c                       dgrhodb(ipt,2,3,lpert+2) 
c    &                  = dgrhodb(ipt,2,3,lpert+2)
c    &                  - 2*bdens(iki+mu)*bfng_val(ipt,mu,2)*
c    &                    bfng_val(ipt,nu,3)
c                       dgrhodb(ipt,2,3,lpert+3) 
c    &                  = dgrhodb(ipt,2,3,lpert+3)
c    &                  - 2*bdens(iki+mu)*bfng_val(ipt,mu,3)*
c    &                    bfng_val(ipt,nu,3)
c
c                       dgrhodb(ipt,2,1,lpert+1) 
c    &                  = dgrhodb(ipt,2,1,lpert+1)
c    &                  - 2*bdens(iki+mu)*bfn_hess(ipt,mu,hxx)*
c    &                    bfn_val(ipt,nu)
c                       dgrhodb(ipt,2,1,lpert+2) 
c    &                  = dgrhodb(ipt,2,1,lpert+2)
c    &                  - 2*bdens(iki+mu)*bfn_hess(ipt,mu,hxy)*
c    &                    bfn_val(ipt,nu)
c                       dgrhodb(ipt,2,1,lpert+3) 
c    &                  = dgrhodb(ipt,2,1,lpert+3)
c    &                  - 2*bdens(iki+mu)*bfn_hess(ipt,mu,hxz)*
c    &                    bfn_val(ipt,nu)
c                       dgrhodb(ipt,2,2,lpert+1) 
c    &                  = dgrhodb(ipt,2,2,lpert+1)
c    &                  - 2*bdens(iki+mu)*bfn_hess(ipt,mu,hxy)*
c    &                    bfn_val(ipt,nu)
c                       dgrhodb(ipt,2,2,lpert+2) 
c    &                  = dgrhodb(ipt,2,2,lpert+2)
c    &                  - 2*bdens(iki+mu)*bfn_hess(ipt,mu,hyy)*
c    &                    bfn_val(ipt,nu)
c                       dgrhodb(ipt,2,2,lpert+3) 
c    &                  = dgrhodb(ipt,2,2,lpert+3)
c    &                  - 2*bdens(iki+mu)*bfn_hess(ipt,mu,hyz)*
c    &                    bfn_val(ipt,nu)
c                       dgrhodb(ipt,2,3,lpert+1) 
c    &                  = dgrhodb(ipt,2,3,lpert+1)
c    &                  - 2*bdens(iki+mu)*bfn_hess(ipt,mu,hxz)*
c    &                    bfn_val(ipt,nu)
c                       dgrhodb(ipt,2,3,lpert+2) 
c    &                  = dgrhodb(ipt,2,3,lpert+2)
c    &                  - 2*bdens(iki+mu)*bfn_hess(ipt,mu,hyz)*
c    &                    bfn_val(ipt,nu)
c                       dgrhodb(ipt,2,3,lpert+3) 
c    &                  = dgrhodb(ipt,2,3,lpert+3)
c    &                  - 2*bdens(iki+mu)*bfn_hess(ipt,mu,hzz)*
c    &                    bfn_val(ipt,nu)
c                    enddo
c                 enddo
c              enddo
c           endif
c        endif
c100  continue
c
cDEBUG
c     write(*,*)'*** drhodb end'
cDEBUG
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine den_pert_exp_ao_scr(rkstyp_sw,gradcorr_sw,
     &     active_bfn_list,active_bfn_indx,
     &     n_active_bfn,n_active_atm,
     &     nao,npts,bfn_val,bfng_val,bfn_hess,
     &     adens,bdens,drhodb,dgrhodb,mxp,dentol)
      implicit none
c
c     Computes the derivatives of rho, grho, and gamma with respect
c     to a nuclear coordinate. Only the terms stemming from 
c     differentiating the explicit functions of the geometry are 
c     included. 
c
c     Parameters:
c
C *Parameter module for CCP1/DFT            
C *Memory information		
C *----------------
C * max_block 	- maximum number of blocks which can be allocate
      integer max_block
      parameter(max_block=20)
C *							
C *Basis sets information                              
C *----------------------				
C *max_tag 	- maximum number of basis sets which can be alloc
C *max_atype	- maximum number of centre types	
C *max_gtype	- maximum number of grid types (at least every element 
C                 has different grid type)
C *max_grids    - maximum number of grids (different terms may have
C                 different grids, i.e. CPKS equations use a different
C                 grid than the one used for the KS-matrix)
C *max_shel	- maximum number of shells on centre
C *max_prm	- maximum number of primitives for any given centre
C *maxL		- maximum angular momentum allowed	
C *max_func	- maximum number of basis functions for any centre
      integer max_tag,max_atype,max_shel,max_prm,max_ang
      integer max_gtype,max_grids
      parameter(max_tag=3,max_atype=30,max_shel=500,max_prm=5000)
      parameter(max_ang=5,max_gtype=10,max_grids=2)
C *								
C *Geometry information					
C *--------------------					
C *max_atom	- maximum number of atoms in system
      integer max_atom
      parameter(max_atom=750)
C *
C *Accuracy information				
C ---------------------			
C *global_accuracy 	- global accuracy	
      real*8  global_accuracy
      parameter(global_accuracy=1.0d-14)

C *Grid information
c *----------------
c *maxradzn - The maximum number of radial zones. 
      integer maxgpt,maxrad,maxfpt,maxang,maxradzn,maxtablerows
      parameter(maxgpt=2900,maxrad=50,maxang=302,maxfpt=100)
      parameter(maxradzn=35,maxtablerows=7)
c
c     The parameters defined here label the 1st, 2nd and 3rd derivative
c     components of the basis functions.
c
c     2nd derivatives:
c
      integer hxx, hyy, hzz, hxy, hxz, hyz
c
      parameter (hxx = 1)
      parameter (hyy = 2)
      parameter (hzz = 3)
      parameter (hxy = 4)
      parameter (hxz = 5)
      parameter (hyz = 6)
c
      integer hcc(3,3)
      data hcc/hxx,hxy,hxz,
     &         hxy,hyy,hyz,
     &         hxz,hyz,hzz/
c
c     3rd derivatives:
c
      integer txxx, txxy, txxz, txyy, txyz
      integer txzz, tyyy, tyyz, tyzz, tzzz
c
      parameter (txxx =  1)
      parameter (txxy =  2)
      parameter (txxz =  3)
      parameter (txyy =  4)
      parameter (txyz =  5)
      parameter (txzz =  6)
      parameter (tyyy =  7)
      parameter (tyyz =  8)
      parameter (tyzz =  9)
      parameter (tzzz = 10)
c
      integer tc3(3,3,3)
      data tc3/txxx,txxy,txxz,
     &         txxy,txyy,txyz,
     &         txxz,txyz,txzz,
     &         txxy,txyy,txyz,
     &         txyy,tyyy,tyyz,
     &         txyz,tyyz,tyzz,
     &         txxz,txyz,txzz,
     &         txyz,tyyz,tyzz,
     &         txzz,tyzz,tzzz/
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
c     In variables:
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
C *Basis sets include file
c
      real*8  alpha
      real*8  cont_coeff
      integer num_bset,bset_tags
      common/basis_sets/alpha(max_tag,max_atype,max_prm),
     &                  cont_coeff(max_tag,max_atype,max_prm,max_ang),
     &                  num_bset,bset_tags(max_tag)

      integer Ashl, Aprm, Abfn, totshl, totprm, totbfn, size_shlA
      integer size_basA, size_primA, maxi_shlA, maxi_basA, maxi_primA
      integer num_types, atom_tag , num_shl, atm_typ, nprim, angmom
      integer hybrid, pstart

      common/basis_size_info/Ashl(max_tag,max_atype),
     &                       Aprm(max_tag,max_atype),
     &                       Abfn(max_tag,max_atype),
     &                       totshl(max_tag),
     &                       totprm(max_tag),
     &                       totbfn(max_tag),
     &                       size_shlA,maxi_shlA,
     &                       size_basA,maxi_basA,
     &                       size_primA,maxi_primA
C
C Descriptions
C 
C  num_types		-	number of types of atoms for a given basis
C  atom_tag		-	list of atom type id numbers
C
      common/basis_cent_info/num_types(max_tag),
     &                       atom_tag(max_tag,max_atom)
C
C Descriptions
C
C  num_shl		-	number of shells for given atom type
C  atm_typ		-	atomic number of atom type
C  nprim		-  	number of primitives 
C  angmom		-	angular momentum
C  hybrid		-	level of hybrid shell. Is same as angmom if
C				shell is not a hybrid. For hybrid shells is
C				always less than angmom e.g. for an sp shell
C				angmom=2, hybrid=1
C  pstart		-	start of exponents and contraction coeffs 
C				contained in include basis.hf77, for a given
C				shell
C
      common/basis_cont_info/num_shl(max_tag,max_atype),
     &                       atm_typ(max_tag,max_atype),
     &                       nprim(max_tag,max_atype,max_shel),
     &                       angmom(max_tag,max_atype,max_shel),
     &                       hybrid(max_tag,max_atype,max_shel),
     &                       pstart(max_tag,max_atype,max_shel)

c
      logical rkstyp_sw ! .true. if closed shell
      logical gradcorr_sw
      integer n_active_bfn
      integer n_active_atm
      integer active_bfn_list(n_active_bfn)
      integer active_bfn_indx(n_active_bfn)
      integer nao, npts, mxp
      real*8 bfn_val(mxp,nao)
      real*8 bfng_val(mxp,nao,3)
      real*8 bfn_hess(mxp,nao,6)
      real*8 adens(nao*(nao+1)/2)
      real*8 bdens(nao*(nao+1)/2)
      real*8 dentol
c
c     Out variables:
c
      real*8 drhodb(mxp,2,3*n_active_atm)
      real*8 dgrhodb(mxp,2,3,3*n_active_atm)
c
c     Local variables:
c
      integer lpert, mu, nu, iki, im, in
      integer ipt
      integer k, l
      real*8 dd
c
c     Code:
c
c     For the time being I don't want to worry about using translational
c     invariance.
c     if (npert.ne.3*ngridcentres) call caserr('not implemented yet 1')
c     if (gradcorr_sw) call caserr('not implemented yet 2')
c
      call aclear_dp(drhodb,mxp*2*3*n_active_atm,0.0d0)
      if (gradcorr_sw) then
         call aclear_dp(dgrhodb,mxp*2*3*3*n_active_atm,0.0d0)
      endif
c
      if (rkstyp_sw) then
         do im = 1, n_active_bfn
            lpert = 3*(active_bfn_indx(im)-1)
            mu = active_bfn_list(im)
            iki = iky(mu)
            do in = 1, im
               nu = active_bfn_list(in)
               dd = adens(iki+nu)
               if (abs(dd).gt.dentol) then
                  do k=1,3
                     do ipt=1,npts
                        drhodb(ipt,1,lpert+k) = drhodb(ipt,1,lpert+k)
     &                  - dd*bfng_val(ipt,im,k)*bfn_val(ipt,in)
                     enddo ! ipt
                  enddo ! k
               endif
            enddo ! in
            do in = im+1, n_active_bfn
               nu = active_bfn_list(in)
               iki = iky(nu)
               dd = adens(iki+mu)
               if (abs(dd).gt.dentol) then
                  do k=1,3
                     do ipt=1,npts
                        drhodb(ipt,1,lpert+k) = drhodb(ipt,1,lpert+k)
     &                  - dd*bfng_val(ipt,im,k)*bfn_val(ipt,in)
                     enddo ! ipt
                  enddo ! k
               endif
            enddo ! in
         enddo ! im

         if (gradcorr_sw) then
            do im = 1, n_active_bfn
               lpert = 3*(active_bfn_indx(im)-1)
               mu = active_bfn_list(im)
               iki = iky(mu)
               do in = 1, im
                  nu = active_bfn_list(in)
                  dd = adens(iki+nu)
                  if (abs(dd).gt.dentol) then
                     do l=1,3
                        do k=1,3
                           do ipt=1,npts
                              dgrhodb(ipt,1,k,lpert+l) 
     &                        = dgrhodb(ipt,1,k,lpert+l)
     &                        - dd*bfng_val(ipt,im,l)*
     &                          bfng_val(ipt,in,k)
     &                        - dd*bfn_hess(ipt,im,hcc(k,l))*
     &                          bfn_val(ipt,in)
                           enddo ! ipt
                        enddo ! k
                     enddo ! l
                  endif
               enddo ! in
               do in = im+1, n_active_bfn
                  nu = active_bfn_list(in)
                  iki = iky(nu)
                  dd = adens(iki+mu)
                  if (abs(dd).gt.dentol) then
                     do l=1,3
                        do k=1,3
                           do ipt=1,npts
                              dgrhodb(ipt,1,k,lpert+l) 
     &                        = dgrhodb(ipt,1,k,lpert+l)
     &                        - dd*bfng_val(ipt,im,l)*
     &                          bfng_val(ipt,in,k)
     &                        - dd*bfn_hess(ipt,im,hcc(k,l))*
     &                          bfn_val(ipt,in)
                           enddo ! ipt
                        enddo ! k
                     enddo ! l
                  endif
               enddo ! in
            enddo ! im
         endif ! gradcorr_sw
      else ! rkstyp_sw
         do im = 1, n_active_bfn
            lpert = 3*(active_bfn_indx(im)-1)
            mu = active_bfn_list(im)
            iki = iky(mu)
            do in = 1, im
               nu = active_bfn_list(in)
               dd = 2*adens(iki+nu)
               if (abs(dd).gt.dentol) then
                  do k=1,3
                     do ipt=1,npts
                        drhodb(ipt,1,lpert+k) = drhodb(ipt,1,lpert+k)
     &                  - dd*bfng_val(ipt,im,k)*
     &                    bfn_val(ipt,in)
                     enddo ! ipt
                  enddo
               endif
               dd = 2*bdens(iki+nu)
               if (abs(dd).gt.dentol) then
                  do k=1,3
                     do ipt=1,npts
                        drhodb(ipt,2,lpert+k) = drhodb(ipt,2,lpert+k)
     &                  - dd*bfng_val(ipt,im,k)*
     &                    bfn_val(ipt,in)
                     enddo ! ipt
                  enddo ! k
               endif
            enddo ! in
            do in = im+1, n_active_bfn
               nu = active_bfn_list(in)
               iki = iky(nu)
               dd = 2*adens(iki+mu)
               if (abs(dd).gt.dentol) then
                  do k=1,3
                     do ipt=1,npts
                        drhodb(ipt,1,lpert+k) = drhodb(ipt,1,lpert+k)
     &                  - dd*bfng_val(ipt,im,k)*
     &                    bfn_val(ipt,in)
                     enddo ! ipt
                  enddo
               endif
               dd = 2*adens(iki+mu)
               if (abs(dd).gt.dentol) then
                  do k=1,3
                     do ipt=1,npts
                        drhodb(ipt,2,lpert+k) = drhodb(ipt,2,lpert+k)
     &                  - dd*bfng_val(ipt,im,k)*
     &                    bfn_val(ipt,in)
                     enddo ! ipt
                  enddo ! k
               endif
            enddo ! in
         enddo ! im

         if (gradcorr_sw) then
            do im = 1, n_active_bfn
               lpert = 3*(active_bfn_indx(im)-1)
               mu = active_bfn_list(im)
               iki = iky(mu)
               do in = 1, im
                  nu = active_bfn_list(in)
                  dd = 2*adens(iki+nu)
                  if (abs(dd).gt.dentol) then
                     do l=1,3
                        do k=1,3
                           do ipt=1,npts
                              dgrhodb(ipt,1,k,lpert+l) 
     &                        = dgrhodb(ipt,1,k,lpert+l)
     &                        - dd*bfng_val(ipt,im,l)*
     &                          bfng_val(ipt,in,k)
     &                        - dd*bfn_hess(ipt,im,hcc(k,l))*
     &                          bfn_val(ipt,in)
                           enddo ! ipt
                        enddo
                     enddo
                  endif
                  dd = 2*bdens(iki+nu)
                  if (abs(dd).gt.dentol) then
                     do l=1,3
                        do k=1,3
                           do ipt=1,npts
                              dgrhodb(ipt,2,k,lpert+l) 
     &                        = dgrhodb(ipt,2,k,lpert+l)
     &                        - dd*bfng_val(ipt,im,l)*
     &                          bfng_val(ipt,in,k)
     &                        - dd*bfn_hess(ipt,im,hcc(k,l))*
     &                          bfn_val(ipt,in)
                           enddo ! ipt
                        enddo ! k
                     enddo ! l
                  endif
               enddo ! in
               do in = im+1, n_active_bfn
                  nu = active_bfn_list(in)
                  iki = iky(nu)
                  dd = 2*adens(iki+mu)
                  if (abs(dd).gt.dentol) then
                     do l=1,3
                        do k=1,3
                           do ipt=1,npts
                              dgrhodb(ipt,1,k,lpert+l) 
     &                        = dgrhodb(ipt,1,k,lpert+l)
     &                        - dd*bfng_val(ipt,im,l)*
     &                          bfng_val(ipt,in,k)
     &                        - dd*bfn_hess(ipt,im,hcc(k,l))*
     &                          bfn_val(ipt,in)
                           enddo ! ipt
                        enddo
                     enddo
                  endif
                  dd = 2*bdens(iki+mu)
                  if (abs(dd).gt.dentol) then
                     do l=1,3
                        do k=1,3
                           do ipt=1,npts
                              dgrhodb(ipt,2,k,lpert+l) 
     &                        = dgrhodb(ipt,2,k,lpert+l)
     &                        - dd*bfng_val(ipt,im,l)*
     &                          bfng_val(ipt,in,k)
     &                        - dd*bfn_hess(ipt,im,hcc(k,l))*
     &                          bfn_val(ipt,in)
                           enddo ! ipt
                        enddo ! k
                     enddo ! l
                  endif
               enddo ! in
            enddo ! im
         endif ! gradcorr_sw
      endif ! rkstype_sw

      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine dksm_exp_dft_ao(rkstyp_sw,gradcorr_sw,gwt_sw,
     &     gwt_avail_sw,nao,npts,npert,ngridcentres,iatom,wt,gwt,
     &     bfn_val,bfng_val,bfn_hess,grho,drhodb,dgrhodb,
     &     xc_vpt,xc_dvpt,xc_hpt,xc_dhpt,
     &     qa,qb,mxp)
      implicit none
c
c     Calculates the DFT contributions to the derivative Fock matrices.
c
c     Parameters:
c
C *Parameter module for CCP1/DFT            
C *Memory information		
C *----------------
C * max_block 	- maximum number of blocks which can be allocate
      integer max_block
      parameter(max_block=20)
C *							
C *Basis sets information                              
C *----------------------				
C *max_tag 	- maximum number of basis sets which can be alloc
C *max_atype	- maximum number of centre types	
C *max_gtype	- maximum number of grid types (at least every element 
C                 has different grid type)
C *max_grids    - maximum number of grids (different terms may have
C                 different grids, i.e. CPKS equations use a different
C                 grid than the one used for the KS-matrix)
C *max_shel	- maximum number of shells on centre
C *max_prm	- maximum number of primitives for any given centre
C *maxL		- maximum angular momentum allowed	
C *max_func	- maximum number of basis functions for any centre
      integer max_tag,max_atype,max_shel,max_prm,max_ang
      integer max_gtype,max_grids
      parameter(max_tag=3,max_atype=30,max_shel=500,max_prm=5000)
      parameter(max_ang=5,max_gtype=10,max_grids=2)
C *								
C *Geometry information					
C *--------------------					
C *max_atom	- maximum number of atoms in system
      integer max_atom
      parameter(max_atom=750)
C *
C *Accuracy information				
C ---------------------			
C *global_accuracy 	- global accuracy	
      real*8  global_accuracy
      parameter(global_accuracy=1.0d-14)

C *Grid information
c *----------------
c *maxradzn - The maximum number of radial zones. 
      integer maxgpt,maxrad,maxfpt,maxang,maxradzn,maxtablerows
      parameter(maxgpt=2900,maxrad=50,maxang=302,maxfpt=100)
      parameter(maxradzn=35,maxtablerows=7)
c
c     The parameters defined here label the 1st and 2nd derivative
c     components of the density functional. Because a functional
c     may be purely local or of a GGA type the components are grouped
c     into derivative with respect to the density only and others.
c
c     1st derivatives:
c
      integer ira, irb
      integer igaa, igab, igbb
      integer ita, itb
c
c     2nd derivatives:
c
      integer irara, irarb, irbrb
      integer iragaa, iragab, iragbb, irbgaa, irbgab, irbgbb
      integer igaagaa, igaagab, igaagbb, igabgab, igabgbb, igbbgbb
      integer irata, iratb, irbta, irbtb
      integer igaata, igaatb, igabta, igabtb, igbbta, igbbtb
      integer itata, itatb, itbtb
c
      parameter (ira     = 1)
      parameter (irb     = 2)
c
      parameter (igaa    = 1)
      parameter (igab    = 2)
      parameter (igbb    = 3)
c
      parameter (ita     = 1)
      parameter (itb     = 2)
c
c
      parameter (irara   = 1)
      parameter (irarb   = 2)
      parameter (irbrb   = 3)
c
      parameter (iragaa  = 1)
      parameter (iragab  = 3)
      parameter (iragbb  = 4)
      parameter (irbgaa  = 5)
      parameter (irbgab  = 6)
      parameter (irbgbb  = 7)
      parameter (igaagaa = 2)
      parameter (igaagab = 8)
      parameter (igaagbb = 9)
      parameter (igabgab = 10)
      parameter (igabgbb = 11)
      parameter (igbbgbb = 12)
c
      parameter (irata   = 1)
      parameter (iratb   = 3)
      parameter (irbta   = 4)
      parameter (irbtb   = 5)
      parameter (itata   = 2)
      parameter (itatb   = 6)
      parameter (itbtb   = 7)
c
      parameter (igaata  = 1)
      parameter (igaatb  = 2)
      parameter (igabta  = 3)
      parameter (igabtb  = 4)
      parameter (igbbta  = 5)
      parameter (igbbtb  = 6)
c
c     The parameters defined here label the 1st, 2nd and 3rd derivative
c     components of the basis functions.
c
c     2nd derivatives:
c
      integer hxx, hyy, hzz, hxy, hxz, hyz
c
      parameter (hxx = 1)
      parameter (hyy = 2)
      parameter (hzz = 3)
      parameter (hxy = 4)
      parameter (hxz = 5)
      parameter (hyz = 6)
c
      integer hcc(3,3)
      data hcc/hxx,hxy,hxz,
     &         hxy,hyy,hyz,
     &         hxz,hyz,hzz/
c
c     3rd derivatives:
c
      integer txxx, txxy, txxz, txyy, txyz
      integer txzz, tyyy, tyyz, tyzz, tzzz
c
      parameter (txxx =  1)
      parameter (txxy =  2)
      parameter (txxz =  3)
      parameter (txyy =  4)
      parameter (txyz =  5)
      parameter (txzz =  6)
      parameter (tyyy =  7)
      parameter (tyyz =  8)
      parameter (tyzz =  9)
      parameter (tzzz = 10)
c
      integer tc3(3,3,3)
      data tc3/txxx,txxy,txxz,
     &         txxy,txyy,txyz,
     &         txxz,txyz,txzz,
     &         txxy,txyy,txyz,
     &         txyy,tyyy,tyyz,
     &         txyz,tyyz,tyzz,
     &         txxz,txyz,txzz,
     &         txyz,tyyz,tyzz,
     &         txzz,tyzz,tzzz/
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
c     In variables:
c
c
C *Basis sets include file
c
      real*8  alpha
      real*8  cont_coeff
      integer num_bset,bset_tags
      common/basis_sets/alpha(max_tag,max_atype,max_prm),
     &                  cont_coeff(max_tag,max_atype,max_prm,max_ang),
     &                  num_bset,bset_tags(max_tag)

      integer Ashl, Aprm, Abfn, totshl, totprm, totbfn, size_shlA
      integer size_basA, size_primA, maxi_shlA, maxi_basA, maxi_primA
      integer num_types, atom_tag , num_shl, atm_typ, nprim, angmom
      integer hybrid, pstart

      common/basis_size_info/Ashl(max_tag,max_atype),
     &                       Aprm(max_tag,max_atype),
     &                       Abfn(max_tag,max_atype),
     &                       totshl(max_tag),
     &                       totprm(max_tag),
     &                       totbfn(max_tag),
     &                       size_shlA,maxi_shlA,
     &                       size_basA,maxi_basA,
     &                       size_primA,maxi_primA
C
C Descriptions
C 
C  num_types		-	number of types of atoms for a given basis
C  atom_tag		-	list of atom type id numbers
C
      common/basis_cent_info/num_types(max_tag),
     &                       atom_tag(max_tag,max_atom)
C
C Descriptions
C
C  num_shl		-	number of shells for given atom type
C  atm_typ		-	atomic number of atom type
C  nprim		-  	number of primitives 
C  angmom		-	angular momentum
C  hybrid		-	level of hybrid shell. Is same as angmom if
C				shell is not a hybrid. For hybrid shells is
C				always less than angmom e.g. for an sp shell
C				angmom=2, hybrid=1
C  pstart		-	start of exponents and contraction coeffs 
C				contained in include basis.hf77, for a given
C				shell
C
      common/basis_cont_info/num_shl(max_tag,max_atype),
     &                       atm_typ(max_tag,max_atype),
     &                       nprim(max_tag,max_atype,max_shel),
     &                       angmom(max_tag,max_atype,max_shel),
     &                       hybrid(max_tag,max_atype,max_shel),
     &                       pstart(max_tag,max_atype,max_shel)

c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      logical rkstyp_sw   ! .true. if closed shell
      logical gradcorr_sw
      logical gwt_sw       ! .true. if gradients of quadrature used
      logical gwt_avail_sw ! .true. if gradients of quadrature present
      integer nao, npts, npert, ngridcentres, mxp
      integer iatom
      real*8 wt(mxp)
      real*8 gwt(3,mxp,*)
      real*8 bfn_val(mxp,nao)
      real*8 bfng_val(mxp,nao,3)
      real*8 bfn_hess(mxp,nao,6)
      real*8 grho(mxp,2,3)
      real*8 drhodb(mxp,2,npert)    ! explicit derivatives only
      real*8 dgrhodb(mxp,2,3,npert) ! explicit derivatives only
      real*8 xc_vpt(mxp,2)
      real*8 xc_dvpt(mxp,3)
      real*8 xc_hpt(mxp,3)
      real*8 xc_dhpt(mxp,12)
c
c     Out variables:
c
      real*8 qa((nao+1)*nao/2,npert)
      real*8 qb((nao+1)*nao/2,npert)
c
c     Local variables:
c
      integer latm, ipert, jpert, i, j, k, ipt, n
      integer bstart, bend, mu, nu
      real*8 scr1(3), scr2(3)
      real*8 t1, t2, t3, t
      integer atmt, iki
c
      real*8 gagmn, gbgmn
      real*8 gadbxgmn, gadbygmn, gadbzgmn
      real*8 gbdbxgmn, gbdbygmn, gbdbzgmn
      real*8 gaadbx, gabdbx, gbbdbx
      real*8 gaadby, gabdby, gbbdby
      real*8 gaadbz, gbbdbz, gabdbz
      real*8 gagmndbx, gagmndby, gagmndbz
      real*8 gbgmndbx, gbgmndby, gbgmndbz
c
      real*8 gadbgmn(3)
      real*8 gaadb(3)
c     real*8 gagmndb(3)
c
      integer ao_tag
      parameter(ao_tag=1)
c
c     Code:
c
c     For the time being I don't want to worry about using translational
c     invariance.
      if (npert.ne.3*ngridcentres) call caserr('not implemented yet')
c
cDEBUG
c     if (gwt_sw) then
c     write(*,*)'*** WE should never get here 1'
c     call caserr('problem 1')
c     endif
c     write(*,*)'*** rhs_dft grho   =',grho(3,1,1)
c     write(*,*)'*** rhs_dft drhodb =',drhodb(3,1,1)
c     write(*,*)'*** rhs_dft dgrhodb=',dgrhodb(3,1,1,1)
cDEBUG
      if (rkstyp_sw) then
c
c        the gradient of the weights contributions
c
         if (gwt_sw.and.gwt_avail_sw) then
            do latm=1,ngridcentres
               ipert=3*(latm-1)
               n=0
               do i=1,nao
                  do j=1,i
                     n=n+1
                     do ipt=1,npts
                        do k = 1, 3
                           qa(n,ipert+k)=qa(n,ipert+k) +
     &                       gwt(k,ipt,latm)*xc_vpt(ipt,ira)*
     &                       bfn_val(ipt,i)*bfn_val(ipt,j)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            if (gradcorr_sw) then
               do latm=1,ngridcentres
                  ipert=3*(latm-1)
                  n=0
                  do i=1,nao
                     do j=1,i
                        n=n+1
                        do ipt=1,npts
                           t1 = xc_dvpt(ipt,igaa)*grho(ipt,1,1)
                           t2 = xc_dvpt(ipt,igaa)*grho(ipt,1,2)
                           t3 = xc_dvpt(ipt,igaa)*grho(ipt,1,3)
                           t  = (bfng_val(ipt,i,1)*t1
     &                        +  bfng_val(ipt,i,2)*t2
     &                        +  bfng_val(ipt,i,3)*t3)*bfn_val(ipt,j)
     &                        + (bfng_val(ipt,j,1)*t1
     &                        +  bfng_val(ipt,j,2)*t2
     &                        +  bfng_val(ipt,j,3)*t3)*bfn_val(ipt,i)
                           do k = 1, 3
                              qa(n,ipert+k)=qa(n,ipert+k) +
     &                          0.5d0*gwt(k,ipt,latm)*t
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            endif ! gradcorr_sw
         endif ! gwt_sw
c
         jpert = 3*(iatom-1)
         bend = 0
         do 90 latm=1,ngridcentres
            ipert=3*(latm-1)
            atmt = atom_tag(ao_tag,latm)
            bstart=bend+1
            bend=bend+Abfn(ao_tag,atmt)
            if (latm.eq.iatom.and.gwt_sw) goto 90
            n=0
            do mu=1,nao
               do nu=1,mu
                  n=n+1
                  call aclear_dp(scr1,3,0.0d0)
                  do k=1,3
                     do ipt=1,npts
                        scr1(k)=scr1(k) 
     &                  +wt(ipt)*xc_hpt(ipt,irara)*bfn_val(ipt,mu)*
     &                   bfn_val(ipt,nu)*drhodb(ipt,1,ipert+k)
                     enddo ! ipt
                  enddo ! k
                  qa(n,ipert+1)=qa(n,ipert+1)+scr1(1)
                  qa(n,ipert+2)=qa(n,ipert+2)+scr1(2)
                  qa(n,ipert+3)=qa(n,ipert+3)+scr1(3)
                  if (gwt_sw) then
                    qa(n,jpert+1)=qa(n,jpert+1)-scr1(1)
                    qa(n,jpert+2)=qa(n,jpert+2)-scr1(2)
                    qa(n,jpert+3)=qa(n,jpert+3)-scr1(3)
                  endif ! gwt_sw
 
                  if (gradcorr_sw) then
                     call aclear_dp(scr1,3,0.0d0)
                     do ipt=1,npts
                        gagmn = grho(ipt,1,1)*(
     &                            bfng_val(ipt,mu,1)*bfn_val(ipt,nu)
     &                          + bfng_val(ipt,nu,1)*bfn_val(ipt,mu))
     &                        + grho(ipt,1,2)*(
     &                            bfng_val(ipt,mu,2)*bfn_val(ipt,nu)
     &                          + bfng_val(ipt,nu,2)*bfn_val(ipt,mu))
     &                        + grho(ipt,1,3)*(
     &                            bfng_val(ipt,mu,3)*bfn_val(ipt,nu)
     &                          + bfng_val(ipt,nu,3)*bfn_val(ipt,mu))
                        gagmn = 0.5d0*gagmn
                        do k=1,3
                           gadbgmn(k) 
     &                           = dgrhodb(ipt,1,1,ipert+k)*(
     &                               bfng_val(ipt,mu,1)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,1)*bfn_val(ipt,mu))
     &                           + dgrhodb(ipt,1,2,ipert+k)*(
     &                               bfng_val(ipt,mu,2)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,2)*bfn_val(ipt,mu))
     &                           + dgrhodb(ipt,1,3,ipert+k)*(
     &                               bfng_val(ipt,mu,3)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,3)*bfn_val(ipt,mu))
                           gaadb(k) = (
     &                           dgrhodb(ipt,1,1,ipert+k)*grho(ipt,1,1)
     &                         + dgrhodb(ipt,1,2,ipert+k)*grho(ipt,1,2)
     &                         + dgrhodb(ipt,1,3,ipert+k)*grho(ipt,1,3))
                        enddo ! k
                        do k=1,3
                        scr1(k)=scr1(k) 
     &                  +2*wt(ipt)*xc_dhpt(ipt,iragaa)*gagmn*
     &                     drhodb(ipt,1,ipert+k)
c    &                  +  wt(ipt)*xc_dhpt(ipt,iragab)*gbgmn*
c    &                     drhodb(ipt,1,ipert+k)
c    &                  +2*wt(ipt)*xc_dhpt(ipt,irbgaa)*gagmn*
c    &                     drhodb(ipt,2,ipert+k)
c    &                  +  wt(ipt)*xc_dhpt(ipt,irbgab)*gbgmn*
c    &                     drhodb(ipt,2,ipert+k)
     &                  +  wt(ipt)*xc_dhpt(ipt,iragaa)*gaadb(k)*
     &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
c    &                  +  wt(ipt)*xc_dhpt(ipt,iragab)*gabdbx*
c    &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
c    &                  +  wt(ipt)*xc_dhpt(ipt,iragbb)*gbbdbx*
c    &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
     &                  +  wt(ipt)*xc_dhpt(ipt,igaagaa)*gagmn*gaadb(k)
c    &                  +  wt(ipt)*xc_dhpt(ipt,igaagab)*gbgmn*gaadbx
c    &                  +2*wt(ipt)*xc_dhpt(ipt,igaagab)*gagmn*gabdbx
c    &                  +  wt(ipt)*xc_dhpt(ipt,igabgab)*gbgmn*gabdbx
c    &                  +2*wt(ipt)*xc_dhpt(ipt,igaagbb)*gagmn*gbbdbx
c    &                  +  wt(ipt)*xc_dhpt(ipt,igabgbb)*gbgmn*gbbdbx
     &                  +  wt(ipt)*xc_dvpt(ipt,igaa)*gadbgmn(k)
c    &                  +  wt(ipt)*xc_dvpt(ipt,igab)*gbdbxgmn
                        enddo ! k
                     enddo ! ipt
                     qa(n,ipert+1)=qa(n,ipert+1)+scr1(1)
                     qa(n,ipert+2)=qa(n,ipert+2)+scr1(2)
                     qa(n,ipert+3)=qa(n,ipert+3)+scr1(3)
                     if (gwt_sw) then
                        qa(n,jpert+1)=qa(n,jpert+1)-scr1(1)
                        qa(n,jpert+2)=qa(n,jpert+2)-scr1(2)
                        qa(n,jpert+3)=qa(n,jpert+3)-scr1(3)
                     endif ! gwt_sw
                  endif ! gradcorr_sw
               enddo ! nu
            enddo ! mu
c
            do mu=bstart,bend
               iki=iky(mu)
               do nu=1, mu
c
                  call aclear_dp(scr1,3,0.0d0)
                  do k=1,3
                     do ipt=1,npts
                        scr1(k)=scr1(k)
     &                  -wt(ipt)*xc_vpt(ipt,ira)*bfng_val(ipt,mu,k)*
     &                   bfn_val(ipt,nu)
                     enddo ! ipt
                  enddo ! k
                  if (gradcorr_sw) then
                     do ipt=1,npts
                        gagmndbx = grho(ipt,1,1)*(
     &                           - bfng_val(ipt,mu,1)*bfng_val(ipt,nu,1)
     &                           - bfn_hess(ipt,mu,hxx)*bfn_val(ipt,nu))
     &                           + grho(ipt,1,2)*(
     &                           - bfng_val(ipt,mu,1)*bfng_val(ipt,nu,2)
     &                           - bfn_hess(ipt,mu,hxy)*bfn_val(ipt,nu))
     &                           + grho(ipt,1,3)*(
     &                           - bfng_val(ipt,mu,1)*bfng_val(ipt,nu,3)
     &                           - bfn_hess(ipt,mu,hxz)*bfn_val(ipt,nu))
                        gagmndby = grho(ipt,1,1)*(
     &                           - bfng_val(ipt,mu,2)*bfng_val(ipt,nu,1)
     &                           - bfn_hess(ipt,mu,hxy)*bfn_val(ipt,nu))
     &                           + grho(ipt,1,2)*(
     &                           - bfng_val(ipt,mu,2)*bfng_val(ipt,nu,2)
     &                           - bfn_hess(ipt,mu,hyy)*bfn_val(ipt,nu))
     &                           + grho(ipt,1,3)*(
     &                           - bfng_val(ipt,mu,2)*bfng_val(ipt,nu,3)
     &                           - bfn_hess(ipt,mu,hyz)*bfn_val(ipt,nu))
                        gagmndbz = grho(ipt,1,1)*(
     &                           - bfng_val(ipt,mu,3)*bfng_val(ipt,nu,1)
     &                           - bfn_hess(ipt,mu,hxz)*bfn_val(ipt,nu))
     &                           + grho(ipt,1,2)*(
     &                           - bfng_val(ipt,mu,3)*bfng_val(ipt,nu,2)
     &                           - bfn_hess(ipt,mu,hyz)*bfn_val(ipt,nu))
     &                           + grho(ipt,1,3)*(
     &                           - bfng_val(ipt,mu,3)*bfng_val(ipt,nu,3)
     &                           - bfn_hess(ipt,mu,hzz)*bfn_val(ipt,nu))
                        gagmndbx = 0.5d0*gagmndbx
                        gagmndby = 0.5d0*gagmndby
                        gagmndbz = 0.5d0*gagmndbz
                        scr1(1) = scr1(1)
     &                          +   wt(ipt)*xc_dvpt(ipt,igaa)*gagmndbx
                        scr1(2) = scr1(2)
     &                          +   wt(ipt)*xc_dvpt(ipt,igaa)*gagmndby
                        scr1(3) = scr1(3)
     &                          +   wt(ipt)*xc_dvpt(ipt,igaa)*gagmndbz
                     enddo ! ipt
                  endif ! gradcorr_sw
                  qa(iki+nu,ipert+1)=qa(iki+nu,ipert+1)+scr1(1)
                  qa(iki+nu,ipert+2)=qa(iki+nu,ipert+2)+scr1(2)
                  qa(iki+nu,ipert+3)=qa(iki+nu,ipert+3)+scr1(3)
                  if (gwt_sw) then
                     qa(iki+nu,jpert+1)=qa(iki+nu,jpert+1)-scr1(1)
                     qa(iki+nu,jpert+2)=qa(iki+nu,jpert+2)-scr1(2)
                     qa(iki+nu,jpert+3)=qa(iki+nu,jpert+3)-scr1(3)
                  endif ! gwt_sw
               enddo ! nu
               do nu=mu, nao
                  iki = iky(nu)
c
                  call aclear_dp(scr1,3,0.0d0)
c                 call aclear_dp(scr2,3,0.0d0)
                  do k=1,3
                     do ipt=1,npts
                        scr1(k)=scr1(k)
     &                  -wt(ipt)*xc_vpt(ipt,ira)*bfng_val(ipt,mu,k)*
     &                   bfn_val(ipt,nu)
                     enddo ! ipt
                  enddo ! k
                  if (gradcorr_sw) then
                     do ipt=1,npts
                        gagmndbx = grho(ipt,1,1)*(
     &                           - bfng_val(ipt,mu,1)*bfng_val(ipt,nu,1)
     &                           - bfn_hess(ipt,mu,hxx)*bfn_val(ipt,nu))
     &                           + grho(ipt,1,2)*(
     &                           - bfng_val(ipt,mu,1)*bfng_val(ipt,nu,2)
     &                           - bfn_hess(ipt,mu,hxy)*bfn_val(ipt,nu))
     &                           + grho(ipt,1,3)*(
     &                           - bfng_val(ipt,mu,1)*bfng_val(ipt,nu,3)
     &                           - bfn_hess(ipt,mu,hxz)*bfn_val(ipt,nu))
                        gagmndby = grho(ipt,1,1)*(
     &                           - bfng_val(ipt,mu,2)*bfng_val(ipt,nu,1)
     &                           - bfn_hess(ipt,mu,hxy)*bfn_val(ipt,nu))
     &                           + grho(ipt,1,2)*(
     &                           - bfng_val(ipt,mu,2)*bfng_val(ipt,nu,2)
     &                           - bfn_hess(ipt,mu,hyy)*bfn_val(ipt,nu))
     &                           + grho(ipt,1,3)*(
     &                           - bfng_val(ipt,mu,2)*bfng_val(ipt,nu,3)
     &                           - bfn_hess(ipt,mu,hyz)*bfn_val(ipt,nu))
                        gagmndbz = grho(ipt,1,1)*(
     &                           - bfng_val(ipt,mu,3)*bfng_val(ipt,nu,1)
     &                           - bfn_hess(ipt,mu,hxz)*bfn_val(ipt,nu))
     &                           + grho(ipt,1,2)*(
     &                           - bfng_val(ipt,mu,3)*bfng_val(ipt,nu,2)
     &                           - bfn_hess(ipt,mu,hyz)*bfn_val(ipt,nu))
     &                           + grho(ipt,1,3)*(
     &                           - bfng_val(ipt,mu,3)*bfng_val(ipt,nu,3)
     &                           - bfn_hess(ipt,mu,hzz)*bfn_val(ipt,nu))
                        gagmndbx = 0.5d0*gagmndbx
                        gagmndby = 0.5d0*gagmndby
                        gagmndbz = 0.5d0*gagmndbz
                        scr1(1) = scr1(1)
     &                          +   wt(ipt)*xc_dvpt(ipt,igaa)*gagmndbx
                        scr1(2) = scr1(2)
     &                          +   wt(ipt)*xc_dvpt(ipt,igaa)*gagmndby
                        scr1(3) = scr1(3)
     &                          +   wt(ipt)*xc_dvpt(ipt,igaa)*gagmndbz
                     enddo ! ipt
                  endif ! gradcorr_sw
                  qa(iki+mu,ipert+1)=qa(iki+mu,ipert+1)+scr1(1)
                  qa(iki+mu,ipert+2)=qa(iki+mu,ipert+2)+scr1(2)
                  qa(iki+mu,ipert+3)=qa(iki+mu,ipert+3)+scr1(3)
                  if (gwt_sw) then
                     qa(iki+mu,jpert+1)=qa(iki+mu,jpert+1)-scr1(1)
                     qa(iki+mu,jpert+2)=qa(iki+mu,jpert+2)-scr1(2)
                     qa(iki+mu,jpert+3)=qa(iki+mu,jpert+3)-scr1(3)
                  endif ! gwt_sw
               enddo ! nu
            enddo ! mu
  90     continue
c
      else ! rkstyp_sw
c
c        the gradient of the weights contributions
c
         if (gwt_sw.and.gwt_avail_sw) then
            do latm=1,ngridcentres
               ipert=3*(latm-1)
               n=0
               do i=1,nao
                  do j=1,i
                     n=n+1
                     do ipt=1,npts
                        do k = 1, 3
                           qa(n,ipert+k)=qa(n,ipert+k) +
     &                       gwt(k,ipt,latm)*xc_vpt(ipt,ira)*
     &                       bfn_val(ipt,i)*bfn_val(ipt,j)
                        enddo
                     enddo
                     do ipt=1,npts
                        do k = 1, 3
                           qb(n,ipert+k)=qb(n,ipert+k) +
     &                       gwt(k,ipt,latm)*xc_vpt(ipt,irb)*
     &                       bfn_val(ipt,i)*bfn_val(ipt,j)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            if (gradcorr_sw) then
               do latm=1,ngridcentres
                  ipert=3*(latm-1)
                  n=0
                  do i=1,nao
                     do j=1,i
                        n=n+1
                        do ipt=1,npts
                           t1 = (2.0d0*xc_dvpt(ipt,igaa)*grho(ipt,1,1)
     &                        +        xc_dvpt(ipt,igab)*grho(ipt,2,1))
                           t2 = (2.0d0*xc_dvpt(ipt,igaa)*grho(ipt,1,2)
     &                        +        xc_dvpt(ipt,igab)*grho(ipt,2,2))
                           t3 = (2.0d0*xc_dvpt(ipt,igaa)*grho(ipt,1,3)
     &                        +        xc_dvpt(ipt,igab)*grho(ipt,2,3))
                           t  = (bfng_val(ipt,i,1)*t1
     &                        +  bfng_val(ipt,i,2)*t2
     &                        +  bfng_val(ipt,i,3)*t3)*bfn_val(ipt,j)
     &                        + (bfng_val(ipt,j,1)*t1
     &                        +  bfng_val(ipt,j,2)*t2
     &                        +  bfng_val(ipt,j,3)*t3)*bfn_val(ipt,i)
                           do k = 1, 3
                              qa(n,ipert+k)=qa(n,ipert+k) +
     &                          gwt(k,ipt,latm)*t
                           enddo
                        enddo
                     enddo
                  enddo
                  n=0
                  do i=1,nao
                     do j=1,i
                        n=n+1
                        do ipt=1,npts
                           t1 = (2.0d0*xc_dvpt(ipt,igbb)*grho(ipt,2,1)
     &                        +        xc_dvpt(ipt,igab)*grho(ipt,1,1))
                           t2 = (2.0d0*xc_dvpt(ipt,igbb)*grho(ipt,2,2)
     &                        +        xc_dvpt(ipt,igab)*grho(ipt,1,2))
                           t3 = (2.0d0*xc_dvpt(ipt,igbb)*grho(ipt,2,3)
     &                        +        xc_dvpt(ipt,igab)*grho(ipt,1,3))
                           t  = (bfng_val(ipt,i,1)*t1
     &                        +  bfng_val(ipt,i,2)*t2
     &                        +  bfng_val(ipt,i,3)*t3)*bfn_val(ipt,j)
     &                        + (bfng_val(ipt,j,1)*t1
     &                        +  bfng_val(ipt,j,2)*t2
     &                        +  bfng_val(ipt,j,3)*t3)*bfn_val(ipt,i)
                           do k = 1, 3
                              qb(n,ipert+k)=qb(n,ipert+k) +
     &                          gwt(k,ipt,latm)*t
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            endif ! gradcorr_sw
         endif ! gwt_sw
c
         jpert = 3*(iatom-1)
         bend = 0
         do 100 latm=1,ngridcentres
            ipert=3*(latm-1)
            atmt = atom_tag(ao_tag,latm)
            bstart=bend+1
            bend=bend+Abfn(ao_tag,atmt)
            if (latm.eq.iatom.and.gwt_sw) goto 100
            n=0
            do mu=1,nao
               do nu=1,mu
                  n=n+1
                  call aclear_dp(scr1,3,0.0d0)
                  do ipt=1,npts
                     scr1(1)=scr1(1) 
     &               +wt(ipt)*xc_hpt(ipt,irara)*bfn_val(ipt,mu)*
     &                bfn_val(ipt,nu)*drhodb(ipt,1,ipert+1)
     &               +wt(ipt)*xc_hpt(ipt,irarb)*bfn_val(ipt,mu)*
     &                bfn_val(ipt,nu)*drhodb(ipt,2,ipert+1)
                     scr1(2)=scr1(2) 
     &               +wt(ipt)*xc_hpt(ipt,irara)*bfn_val(ipt,mu)*
     &                bfn_val(ipt,nu)*drhodb(ipt,1,ipert+2)
     &               +wt(ipt)*xc_hpt(ipt,irarb)*bfn_val(ipt,mu)*
     &                bfn_val(ipt,nu)*drhodb(ipt,2,ipert+2)
                     scr1(3)=scr1(3) 
     &               +wt(ipt)*xc_hpt(ipt,irara)*bfn_val(ipt,mu)*
     &                bfn_val(ipt,nu)*drhodb(ipt,1,ipert+3)
     &               +wt(ipt)*xc_hpt(ipt,irarb)*bfn_val(ipt,mu)*
     &                bfn_val(ipt,nu)*drhodb(ipt,2,ipert+3)
                  enddo
                  qa(n,ipert+1)=qa(n,ipert+1)+scr1(1)
                  qa(n,ipert+2)=qa(n,ipert+2)+scr1(2)
                  qa(n,ipert+3)=qa(n,ipert+3)+scr1(3)
                  if (gwt_sw) then
                    qa(n,jpert+1)=qa(n,jpert+1)-scr1(1)
                    qa(n,jpert+2)=qa(n,jpert+2)-scr1(2)
                    qa(n,jpert+3)=qa(n,jpert+3)-scr1(3)
                  endif
 
                  call aclear_dp(scr1,3,0.0d0)
                  do ipt=1,npts
                     scr1(1)=scr1(1) 
     &               +wt(ipt)*xc_hpt(ipt,irbrb)*bfn_val(ipt,mu)*
     &                bfn_val(ipt,nu)*drhodb(ipt,2,ipert+1)
     &               +wt(ipt)*xc_hpt(ipt,irarb)*bfn_val(ipt,mu)*
     &                bfn_val(ipt,nu)*drhodb(ipt,1,ipert+1)
                     scr1(2)=scr1(2) 
     &               +wt(ipt)*xc_hpt(ipt,irbrb)*bfn_val(ipt,mu)*
     &                bfn_val(ipt,nu)*drhodb(ipt,2,ipert+2)
     &               +wt(ipt)*xc_hpt(ipt,irarb)*bfn_val(ipt,mu)*
     &                bfn_val(ipt,nu)*drhodb(ipt,1,ipert+2)
                     scr1(3)=scr1(3) 
     &               +wt(ipt)*xc_hpt(ipt,irbrb)*bfn_val(ipt,mu)*
     &                bfn_val(ipt,nu)*drhodb(ipt,2,ipert+3)
     &               +wt(ipt)*xc_hpt(ipt,irarb)*bfn_val(ipt,mu)*
     &                bfn_val(ipt,nu)*drhodb(ipt,1,ipert+3)
                  enddo
                  qb(n,ipert+1)=qb(n,ipert+1)+scr1(1)
                  qb(n,ipert+2)=qb(n,ipert+2)+scr1(2)
                  qb(n,ipert+3)=qb(n,ipert+3)+scr1(3)
                  if (gwt_sw) then
                     qb(n,jpert+1)=qb(n,jpert+1)-scr1(1)
                     qb(n,jpert+2)=qb(n,jpert+2)-scr1(2)
                     qb(n,jpert+3)=qb(n,jpert+3)-scr1(3)
                  endif

                  if (gradcorr_sw) then
                     call aclear_dp(scr1,3,0.0d0)
                     call aclear_dp(scr2,3,0.0d0)
                     do ipt=1,npts
                        gagmn = grho(ipt,1,1)*(
     &                            bfng_val(ipt,mu,1)*bfn_val(ipt,nu)
     &                          + bfng_val(ipt,nu,1)*bfn_val(ipt,mu))
     &                        + grho(ipt,1,2)*(
     &                            bfng_val(ipt,mu,2)*bfn_val(ipt,nu)
     &                          + bfng_val(ipt,nu,2)*bfn_val(ipt,mu))
     &                        + grho(ipt,1,3)*(
     &                            bfng_val(ipt,mu,3)*bfn_val(ipt,nu)
     &                          + bfng_val(ipt,nu,3)*bfn_val(ipt,mu))
                        gbgmn = grho(ipt,2,1)*(
     &                            bfng_val(ipt,mu,1)*bfn_val(ipt,nu)
     &                          + bfng_val(ipt,nu,1)*bfn_val(ipt,mu))
     &                        + grho(ipt,2,2)*(
     &                            bfng_val(ipt,mu,2)*bfn_val(ipt,nu)
     &                          + bfng_val(ipt,nu,2)*bfn_val(ipt,mu))
     &                        + grho(ipt,2,3)*(
     &                            bfng_val(ipt,mu,3)*bfn_val(ipt,nu)
     &                          + bfng_val(ipt,nu,3)*bfn_val(ipt,mu))
                        gadbxgmn = dgrhodb(ipt,1,1,ipert+1)*(
     &                               bfng_val(ipt,mu,1)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,1)*bfn_val(ipt,mu))
     &                           + dgrhodb(ipt,1,2,ipert+1)*(
     &                               bfng_val(ipt,mu,2)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,2)*bfn_val(ipt,mu))
     &                           + dgrhodb(ipt,1,3,ipert+1)*(
     &                               bfng_val(ipt,mu,3)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,3)*bfn_val(ipt,mu))
                        gbdbxgmn = dgrhodb(ipt,2,1,ipert+1)*(
     &                               bfng_val(ipt,mu,1)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,1)*bfn_val(ipt,mu))
     &                           + dgrhodb(ipt,2,2,ipert+1)*(
     &                               bfng_val(ipt,mu,2)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,2)*bfn_val(ipt,mu))
     &                           + dgrhodb(ipt,2,3,ipert+1)*(
     &                               bfng_val(ipt,mu,3)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,3)*bfn_val(ipt,mu))
                        gadbygmn = dgrhodb(ipt,1,1,ipert+2)*(
     &                               bfng_val(ipt,mu,1)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,1)*bfn_val(ipt,mu))
     &                           + dgrhodb(ipt,1,2,ipert+2)*(
     &                               bfng_val(ipt,mu,2)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,2)*bfn_val(ipt,mu))
     &                           + dgrhodb(ipt,1,3,ipert+2)*(
     &                               bfng_val(ipt,mu,3)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,3)*bfn_val(ipt,mu))
                        gbdbygmn = dgrhodb(ipt,2,1,ipert+2)*(
     &                               bfng_val(ipt,mu,1)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,1)*bfn_val(ipt,mu))
     &                           + dgrhodb(ipt,2,2,ipert+2)*(
     &                               bfng_val(ipt,mu,2)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,2)*bfn_val(ipt,mu))
     &                           + dgrhodb(ipt,2,3,ipert+2)*(
     &                               bfng_val(ipt,mu,3)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,3)*bfn_val(ipt,mu))
                        gadbzgmn = dgrhodb(ipt,1,1,ipert+3)*(
     &                               bfng_val(ipt,mu,1)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,1)*bfn_val(ipt,mu))
     &                           + dgrhodb(ipt,1,2,ipert+3)*(
     &                               bfng_val(ipt,mu,2)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,2)*bfn_val(ipt,mu))
     &                           + dgrhodb(ipt,1,3,ipert+3)*(
     &                               bfng_val(ipt,mu,3)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,3)*bfn_val(ipt,mu))
                        gbdbzgmn = dgrhodb(ipt,2,1,ipert+3)*(
     &                               bfng_val(ipt,mu,1)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,1)*bfn_val(ipt,mu))
     &                           + dgrhodb(ipt,2,2,ipert+3)*(
     &                               bfng_val(ipt,mu,2)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,2)*bfn_val(ipt,mu))
     &                           + dgrhodb(ipt,2,3,ipert+3)*(
     &                               bfng_val(ipt,mu,3)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,3)*bfn_val(ipt,mu))
                        gaadbx = 2*(
     &                           dgrhodb(ipt,1,1,ipert+1)*grho(ipt,1,1)
     &                         + dgrhodb(ipt,1,2,ipert+1)*grho(ipt,1,2)
     &                         + dgrhodb(ipt,1,3,ipert+1)*grho(ipt,1,3))
                        gbbdbx = 2*(
     &                           dgrhodb(ipt,2,1,ipert+1)*grho(ipt,2,1)
     &                         + dgrhodb(ipt,2,2,ipert+1)*grho(ipt,2,2)
     &                         + dgrhodb(ipt,2,3,ipert+1)*grho(ipt,2,3))
                        gabdbx = dgrhodb(ipt,1,1,ipert+1)*grho(ipt,2,1)
     &                         + dgrhodb(ipt,1,2,ipert+1)*grho(ipt,2,2)
     &                         + dgrhodb(ipt,1,3,ipert+1)*grho(ipt,2,3)
     &                         + dgrhodb(ipt,2,1,ipert+1)*grho(ipt,1,1)
     &                         + dgrhodb(ipt,2,2,ipert+1)*grho(ipt,1,2)
     &                         + dgrhodb(ipt,2,3,ipert+1)*grho(ipt,1,3)
                        gaadby = 2*(
     &                           dgrhodb(ipt,1,1,ipert+2)*grho(ipt,1,1)
     &                         + dgrhodb(ipt,1,2,ipert+2)*grho(ipt,1,2)
     &                         + dgrhodb(ipt,1,3,ipert+2)*grho(ipt,1,3))
                        gbbdby = 2*(
     &                           dgrhodb(ipt,2,1,ipert+2)*grho(ipt,2,1)
     &                         + dgrhodb(ipt,2,2,ipert+2)*grho(ipt,2,2)
     &                         + dgrhodb(ipt,2,3,ipert+2)*grho(ipt,2,3))
                        gabdby = dgrhodb(ipt,1,1,ipert+2)*grho(ipt,2,1)
     &                         + dgrhodb(ipt,1,2,ipert+2)*grho(ipt,2,2)
     &                         + dgrhodb(ipt,1,3,ipert+2)*grho(ipt,2,3)
     &                         + dgrhodb(ipt,2,1,ipert+2)*grho(ipt,1,1)
     &                         + dgrhodb(ipt,2,2,ipert+2)*grho(ipt,1,2)
     &                         + dgrhodb(ipt,2,3,ipert+2)*grho(ipt,1,3)
                        gaadbz = 2*(
     &                           dgrhodb(ipt,1,1,ipert+3)*grho(ipt,1,1)
     &                         + dgrhodb(ipt,1,2,ipert+3)*grho(ipt,1,2)
     &                         + dgrhodb(ipt,1,3,ipert+3)*grho(ipt,1,3))
                        gbbdbz = 2*(
     &                           dgrhodb(ipt,2,1,ipert+3)*grho(ipt,2,1)
     &                         + dgrhodb(ipt,2,2,ipert+3)*grho(ipt,2,2)
     &                         + dgrhodb(ipt,2,3,ipert+3)*grho(ipt,2,3))
                        gabdbz = dgrhodb(ipt,1,1,ipert+3)*grho(ipt,2,1)
     &                         + dgrhodb(ipt,1,2,ipert+3)*grho(ipt,2,2)
     &                         + dgrhodb(ipt,1,3,ipert+3)*grho(ipt,2,3)
     &                         + dgrhodb(ipt,2,1,ipert+3)*grho(ipt,1,1)
     &                         + dgrhodb(ipt,2,2,ipert+3)*grho(ipt,1,2)
     &                         + dgrhodb(ipt,2,3,ipert+3)*grho(ipt,1,3)
                        scr1(1)=scr1(1) 
     &                  +2*wt(ipt)*xc_dhpt(ipt,iragaa)*gagmn*
     &                     drhodb(ipt,1,ipert+1)
     &                  +  wt(ipt)*xc_dhpt(ipt,iragab)*gbgmn*
     &                     drhodb(ipt,1,ipert+1)
     &                  +2*wt(ipt)*xc_dhpt(ipt,irbgaa)*gagmn*
     &                     drhodb(ipt,2,ipert+1)
     &                  +  wt(ipt)*xc_dhpt(ipt,irbgab)*gbgmn*
     &                     drhodb(ipt,2,ipert+1)
     &                  +  wt(ipt)*xc_dhpt(ipt,iragaa)*gaadbx*
     &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
     &                  +  wt(ipt)*xc_dhpt(ipt,iragab)*gabdbx*
     &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
     &                  +  wt(ipt)*xc_dhpt(ipt,iragbb)*gbbdbx*
     &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
     &                  +2*wt(ipt)*xc_dhpt(ipt,igaagaa)*gagmn*gaadbx
     &                  +  wt(ipt)*xc_dhpt(ipt,igaagab)*gbgmn*gaadbx
     &                  +2*wt(ipt)*xc_dhpt(ipt,igaagab)*gagmn*gabdbx
     &                  +  wt(ipt)*xc_dhpt(ipt,igabgab)*gbgmn*gabdbx
     &                  +2*wt(ipt)*xc_dhpt(ipt,igaagbb)*gagmn*gbbdbx
     &                  +  wt(ipt)*xc_dhpt(ipt,igabgbb)*gbgmn*gbbdbx
     &                  +2*wt(ipt)*xc_dvpt(ipt,igaa)*gadbxgmn
     &                  +  wt(ipt)*xc_dvpt(ipt,igab)*gbdbxgmn
                        scr1(2)=scr1(2) 
     &                  +2*wt(ipt)*xc_dhpt(ipt,iragaa)*gagmn*
     &                     drhodb(ipt,1,ipert+2)
     &                  +  wt(ipt)*xc_dhpt(ipt,iragab)*gbgmn*
     &                     drhodb(ipt,1,ipert+2)
     &                  +2*wt(ipt)*xc_dhpt(ipt,irbgaa)*gagmn*
     &                     drhodb(ipt,2,ipert+2)
     &                  +  wt(ipt)*xc_dhpt(ipt,irbgab)*gbgmn*
     &                     drhodb(ipt,2,ipert+2)
     &                  +  wt(ipt)*xc_dhpt(ipt,iragaa)*gaadby*
     &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
     &                  +  wt(ipt)*xc_dhpt(ipt,iragab)*gabdby*
     &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
     &                  +  wt(ipt)*xc_dhpt(ipt,iragbb)*gbbdby*
     &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
     &                  +2*wt(ipt)*xc_dhpt(ipt,igaagaa)*gagmn*gaadby
     &                  +  wt(ipt)*xc_dhpt(ipt,igaagab)*gbgmn*gaadby
     &                  +2*wt(ipt)*xc_dhpt(ipt,igaagab)*gagmn*gabdby
     &                  +  wt(ipt)*xc_dhpt(ipt,igabgab)*gbgmn*gabdby
     &                  +2*wt(ipt)*xc_dhpt(ipt,igaagbb)*gagmn*gbbdby
     &                  +  wt(ipt)*xc_dhpt(ipt,igabgbb)*gbgmn*gbbdby
     &                  +2*wt(ipt)*xc_dvpt(ipt,igaa)*gadbygmn
     &                  +  wt(ipt)*xc_dvpt(ipt,igab)*gbdbygmn
                        scr1(3)=scr1(3) 
     &                  +2*wt(ipt)*xc_dhpt(ipt,iragaa)*gagmn*
     &                     drhodb(ipt,1,ipert+3)
     &                  +  wt(ipt)*xc_dhpt(ipt,iragab)*gbgmn*
     &                     drhodb(ipt,1,ipert+3)
     &                  +2*wt(ipt)*xc_dhpt(ipt,irbgaa)*gagmn*
     &                     drhodb(ipt,2,ipert+3)
     &                  +  wt(ipt)*xc_dhpt(ipt,irbgab)*gbgmn*
     &                     drhodb(ipt,2,ipert+3)
     &                  +  wt(ipt)*xc_dhpt(ipt,iragaa)*gaadbz*
     &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
     &                  +  wt(ipt)*xc_dhpt(ipt,iragab)*gabdbz*
     &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
     &                  +  wt(ipt)*xc_dhpt(ipt,iragbb)*gbbdbz*
     &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
     &                  +2*wt(ipt)*xc_dhpt(ipt,igaagaa)*gagmn*gaadbz
     &                  +  wt(ipt)*xc_dhpt(ipt,igaagab)*gbgmn*gaadbz
     &                  +2*wt(ipt)*xc_dhpt(ipt,igaagab)*gagmn*gabdbz
     &                  +  wt(ipt)*xc_dhpt(ipt,igabgab)*gbgmn*gabdbz
     &                  +2*wt(ipt)*xc_dhpt(ipt,igaagbb)*gagmn*gbbdbz
     &                  +  wt(ipt)*xc_dhpt(ipt,igabgbb)*gbgmn*gbbdbz
     &                  +2*wt(ipt)*xc_dvpt(ipt,igaa)*gadbzgmn
     &                  +  wt(ipt)*xc_dvpt(ipt,igab)*gbdbzgmn
                        scr2(1)=scr2(1) 
     &                  +2*wt(ipt)*xc_dhpt(ipt,iragbb)*gbgmn*
     &                     drhodb(ipt,1,ipert+1)
     &                  +  wt(ipt)*xc_dhpt(ipt,iragab)*gagmn*
     &                     drhodb(ipt,1,ipert+1)
     &                  +2*wt(ipt)*xc_dhpt(ipt,irbgbb)*gbgmn*
     &                     drhodb(ipt,2,ipert+1)
     &                  +  wt(ipt)*xc_dhpt(ipt,irbgab)*gagmn*
     &                     drhodb(ipt,2,ipert+1)
     &                  +  wt(ipt)*xc_dhpt(ipt,irbgaa)*gaadbx*
     &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
     &                  +  wt(ipt)*xc_dhpt(ipt,irbgab)*gabdbx*
     &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
     &                  +  wt(ipt)*xc_dhpt(ipt,irbgbb)*gbbdbx*
     &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
     &                  +2*wt(ipt)*xc_dhpt(ipt,igaagbb)*gbgmn*gaadbx
     &                  +  wt(ipt)*xc_dhpt(ipt,igaagab)*gagmn*gaadbx
     &                  +2*wt(ipt)*xc_dhpt(ipt,igabgbb)*gbgmn*gabdbx
     &                  +  wt(ipt)*xc_dhpt(ipt,igabgab)*gagmn*gabdbx
     &                  +2*wt(ipt)*xc_dhpt(ipt,igbbgbb)*gbgmn*gbbdbx
     &                  +  wt(ipt)*xc_dhpt(ipt,igabgbb)*gagmn*gbbdbx
     &                  +2*wt(ipt)*xc_dvpt(ipt,igbb)*gbdbxgmn
     &                  +  wt(ipt)*xc_dvpt(ipt,igab)*gadbxgmn
                        scr2(2)=scr2(2) 
     &                  +2*wt(ipt)*xc_dhpt(ipt,iragbb)*gbgmn*
     &                     drhodb(ipt,1,ipert+2)
     &                  +  wt(ipt)*xc_dhpt(ipt,iragab)*gagmn*
     &                     drhodb(ipt,1,ipert+2)
     &                  +2*wt(ipt)*xc_dhpt(ipt,irbgbb)*gbgmn*
     &                     drhodb(ipt,2,ipert+2)
     &                  +  wt(ipt)*xc_dhpt(ipt,irbgab)*gagmn*
     &                     drhodb(ipt,2,ipert+2)
     &                  +  wt(ipt)*xc_dhpt(ipt,irbgaa)*gaadby*
     &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
     &                  +  wt(ipt)*xc_dhpt(ipt,irbgab)*gabdby*
     &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
     &                  +  wt(ipt)*xc_dhpt(ipt,irbgbb)*gbbdby*
     &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
     &                  +2*wt(ipt)*xc_dhpt(ipt,igaagbb)*gbgmn*gaadby
     &                  +  wt(ipt)*xc_dhpt(ipt,igaagab)*gagmn*gaadby
     &                  +2*wt(ipt)*xc_dhpt(ipt,igabgbb)*gbgmn*gabdby
     &                  +  wt(ipt)*xc_dhpt(ipt,igabgab)*gagmn*gabdby
     &                  +2*wt(ipt)*xc_dhpt(ipt,igbbgbb)*gbgmn*gbbdby
     &                  +  wt(ipt)*xc_dhpt(ipt,igabgbb)*gagmn*gbbdby
     &                  +2*wt(ipt)*xc_dvpt(ipt,igbb)*gbdbygmn
     &                  +  wt(ipt)*xc_dvpt(ipt,igab)*gadbygmn
                        scr2(3)=scr2(3) 
     &                  +2*wt(ipt)*xc_dhpt(ipt,iragbb)*gbgmn*
     &                     drhodb(ipt,1,ipert+3)
     &                  +  wt(ipt)*xc_dhpt(ipt,iragab)*gagmn*
     &                     drhodb(ipt,1,ipert+3)
     &                  +2*wt(ipt)*xc_dhpt(ipt,irbgbb)*gbgmn*
     &                     drhodb(ipt,2,ipert+3)
     &                  +  wt(ipt)*xc_dhpt(ipt,irbgab)*gagmn*
     &                     drhodb(ipt,2,ipert+3)
     &                  +  wt(ipt)*xc_dhpt(ipt,irbgaa)*gaadbz*
     &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
     &                  +  wt(ipt)*xc_dhpt(ipt,irbgab)*gabdbz*
     &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
     &                  +  wt(ipt)*xc_dhpt(ipt,irbgbb)*gbbdbz*
     &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
     &                  +2*wt(ipt)*xc_dhpt(ipt,igaagbb)*gbgmn*gaadbz
     &                  +  wt(ipt)*xc_dhpt(ipt,igaagab)*gagmn*gaadbz
     &                  +2*wt(ipt)*xc_dhpt(ipt,igabgbb)*gbgmn*gabdbz
     &                  +  wt(ipt)*xc_dhpt(ipt,igabgab)*gagmn*gabdbz
     &                  +2*wt(ipt)*xc_dhpt(ipt,igbbgbb)*gbgmn*gbbdbz
     &                  +  wt(ipt)*xc_dhpt(ipt,igabgbb)*gagmn*gbbdbz
     &                  +2*wt(ipt)*xc_dvpt(ipt,igbb)*gbdbzgmn
     &                  +  wt(ipt)*xc_dvpt(ipt,igab)*gadbzgmn
                     enddo
                     qa(n,ipert+1)=qa(n,ipert+1)+scr1(1)
                     qa(n,ipert+2)=qa(n,ipert+2)+scr1(2)
                     qa(n,ipert+3)=qa(n,ipert+3)+scr1(3)
                     qb(n,ipert+1)=qb(n,ipert+1)+scr2(1)
                     qb(n,ipert+2)=qb(n,ipert+2)+scr2(2)
                     qb(n,ipert+3)=qb(n,ipert+3)+scr2(3)
                     if (gwt_sw) then
                        qa(n,jpert+1)=qa(n,jpert+1)-scr1(1)
                        qa(n,jpert+2)=qa(n,jpert+2)-scr1(2)
                        qa(n,jpert+3)=qa(n,jpert+3)-scr1(3)
                        qb(n,jpert+1)=qb(n,jpert+1)-scr2(1)
                        qb(n,jpert+2)=qb(n,jpert+2)-scr2(2)
                        qb(n,jpert+3)=qb(n,jpert+3)-scr2(3)
                     endif
                  endif
               enddo
            enddo
c
            if (latm.eq.iatom.and.gwt_sw) goto 100
            do mu=bstart,bend
               iki=iky(mu)
               do nu=1, mu
c
                  call aclear_dp(scr1,3,0.0d0)
                  call aclear_dp(scr2,3,0.0d0)
                  do ipt=1,npts
                     scr1(1)=scr1(1)
     &               -wt(ipt)*xc_vpt(ipt,ira)*bfng_val(ipt,mu,1)*
     &                bfn_val(ipt,nu)
                     scr1(2)=scr1(2)
     &               -wt(ipt)*xc_vpt(ipt,ira)*bfng_val(ipt,mu,2)*
     &                bfn_val(ipt,nu)
                     scr1(3)=scr1(3)
     &               -wt(ipt)*xc_vpt(ipt,ira)*bfng_val(ipt,mu,3)*
     &                bfn_val(ipt,nu)
                     scr2(1)=scr2(1)
     &               -wt(ipt)*xc_vpt(ipt,irb)*bfng_val(ipt,mu,1)*
     &                bfn_val(ipt,nu)
                     scr2(2)=scr2(2)
     &               -wt(ipt)*xc_vpt(ipt,irb)*bfng_val(ipt,mu,2)*
     &                bfn_val(ipt,nu)
                     scr2(3)=scr2(3)
     &               -wt(ipt)*xc_vpt(ipt,irb)*bfng_val(ipt,mu,3)*
     &                bfn_val(ipt,nu)
                  enddo
                  if (gradcorr_sw) then
                     do ipt=1,npts
                        gagmndbx = grho(ipt,1,1)*(
     &                           - bfng_val(ipt,mu,1)*bfng_val(ipt,nu,1)
     &                           - bfn_hess(ipt,mu,hxx)*bfn_val(ipt,nu))
     &                           + grho(ipt,1,2)*(
     &                           - bfng_val(ipt,mu,1)*bfng_val(ipt,nu,2)
     &                           - bfn_hess(ipt,mu,hxy)*bfn_val(ipt,nu))
     &                           + grho(ipt,1,3)*(
     &                           - bfng_val(ipt,mu,1)*bfng_val(ipt,nu,3)
     &                           - bfn_hess(ipt,mu,hxz)*bfn_val(ipt,nu))
                        gbgmndbx = grho(ipt,2,1)*(
     &                           - bfng_val(ipt,mu,1)*bfng_val(ipt,nu,1)
     &                           - bfn_hess(ipt,mu,hxx)*bfn_val(ipt,nu))
     &                           + grho(ipt,2,2)*(
     &                           - bfng_val(ipt,mu,1)*bfng_val(ipt,nu,2)
     &                           - bfn_hess(ipt,mu,hxy)*bfn_val(ipt,nu))
     &                           + grho(ipt,2,3)*(
     &                           - bfng_val(ipt,mu,1)*bfng_val(ipt,nu,3)
     &                           - bfn_hess(ipt,mu,hxz)*bfn_val(ipt,nu))
                        gagmndby = grho(ipt,1,1)*(
     &                           - bfng_val(ipt,mu,2)*bfng_val(ipt,nu,1)
     &                           - bfn_hess(ipt,mu,hxy)*bfn_val(ipt,nu))
     &                           + grho(ipt,1,2)*(
     &                           - bfng_val(ipt,mu,2)*bfng_val(ipt,nu,2)
     &                           - bfn_hess(ipt,mu,hyy)*bfn_val(ipt,nu))
     &                           + grho(ipt,1,3)*(
     &                           - bfng_val(ipt,mu,2)*bfng_val(ipt,nu,3)
     &                           - bfn_hess(ipt,mu,hyz)*bfn_val(ipt,nu))
                        gbgmndby = grho(ipt,2,1)*(
     &                           - bfng_val(ipt,mu,2)*bfng_val(ipt,nu,1)
     &                           - bfn_hess(ipt,mu,hxy)*bfn_val(ipt,nu))
     &                           + grho(ipt,2,2)*(
     &                           - bfng_val(ipt,mu,2)*bfng_val(ipt,nu,2)
     &                           - bfn_hess(ipt,mu,hyy)*bfn_val(ipt,nu))
     &                           + grho(ipt,2,3)*(
     &                           - bfng_val(ipt,mu,2)*bfng_val(ipt,nu,3)
     &                           - bfn_hess(ipt,mu,hyz)*bfn_val(ipt,nu))
                        gagmndbz = grho(ipt,1,1)*(
     &                           - bfng_val(ipt,mu,3)*bfng_val(ipt,nu,1)
     &                           - bfn_hess(ipt,mu,hxz)*bfn_val(ipt,nu))
     &                           + grho(ipt,1,2)*(
     &                           - bfng_val(ipt,mu,3)*bfng_val(ipt,nu,2)
     &                           - bfn_hess(ipt,mu,hyz)*bfn_val(ipt,nu))
     &                           + grho(ipt,1,3)*(
     &                           - bfng_val(ipt,mu,3)*bfng_val(ipt,nu,3)
     &                           - bfn_hess(ipt,mu,hzz)*bfn_val(ipt,nu))
                        gbgmndbz = grho(ipt,2,1)*(
     &                           - bfng_val(ipt,mu,3)*bfng_val(ipt,nu,1)
     &                           - bfn_hess(ipt,mu,hxz)*bfn_val(ipt,nu))
     &                           + grho(ipt,2,2)*(
     &                           - bfng_val(ipt,mu,3)*bfng_val(ipt,nu,2)
     &                           - bfn_hess(ipt,mu,hyz)*bfn_val(ipt,nu))
     &                           + grho(ipt,2,3)*(
     &                           - bfng_val(ipt,mu,3)*bfng_val(ipt,nu,3)
     &                           - bfn_hess(ipt,mu,hzz)*bfn_val(ipt,nu))
                        scr1(1) = scr1(1)
     &                          + 2*wt(ipt)*xc_dvpt(ipt,igaa)*gagmndbx
     &                          +   wt(ipt)*xc_dvpt(ipt,igab)*gbgmndbx
                        scr1(2) = scr1(2)
     &                          + 2*wt(ipt)*xc_dvpt(ipt,igaa)*gagmndby
     &                          +   wt(ipt)*xc_dvpt(ipt,igab)*gbgmndby
                        scr1(3) = scr1(3)
     &                          + 2*wt(ipt)*xc_dvpt(ipt,igaa)*gagmndbz
     &                          +   wt(ipt)*xc_dvpt(ipt,igab)*gbgmndbz
                        scr2(1) = scr2(1)
     &                          + 2*wt(ipt)*xc_dvpt(ipt,igbb)*gbgmndbx
     &                          +   wt(ipt)*xc_dvpt(ipt,igab)*gagmndbx
                        scr2(2) = scr2(2)
     &                          + 2*wt(ipt)*xc_dvpt(ipt,igbb)*gbgmndby
     &                          +   wt(ipt)*xc_dvpt(ipt,igab)*gagmndby
                        scr2(3) = scr2(3)
     &                          + 2*wt(ipt)*xc_dvpt(ipt,igbb)*gbgmndbz
     &                          +   wt(ipt)*xc_dvpt(ipt,igab)*gagmndbz
                     enddo
                  endif
                  qa(iki+nu,ipert+1)=qa(iki+nu,ipert+1)+scr1(1)
                  qa(iki+nu,ipert+2)=qa(iki+nu,ipert+2)+scr1(2)
                  qa(iki+nu,ipert+3)=qa(iki+nu,ipert+3)+scr1(3)
                  if (gwt_sw) then
                     qa(iki+nu,jpert+1)=qa(iki+nu,jpert+1)-scr1(1)
                     qa(iki+nu,jpert+2)=qa(iki+nu,jpert+2)-scr1(2)
                     qa(iki+nu,jpert+3)=qa(iki+nu,jpert+3)-scr1(3)
                  endif
                  qb(iki+nu,ipert+1)=qb(iki+nu,ipert+1)+scr2(1)
                  qb(iki+nu,ipert+2)=qb(iki+nu,ipert+2)+scr2(2)
                  qb(iki+nu,ipert+3)=qb(iki+nu,ipert+3)+scr2(3)
                  if (gwt_sw) then
                     qb(iki+nu,jpert+1)=qb(iki+nu,jpert+1)-scr2(1)
                     qb(iki+nu,jpert+2)=qb(iki+nu,jpert+2)-scr2(2)
                     qb(iki+nu,jpert+3)=qb(iki+nu,jpert+3)-scr2(3)
                  endif
               enddo
               do nu=mu, nao
                  iki = iky(nu)
c
                  call aclear_dp(scr1,3,0.0d0)
                  call aclear_dp(scr2,3,0.0d0)
                  do ipt=1,npts
                     scr1(1)=scr1(1)
     &               -wt(ipt)*xc_vpt(ipt,ira)*bfng_val(ipt,mu,1)*
     &                bfn_val(ipt,nu)
                     scr1(2)=scr1(2)
     &               -wt(ipt)*xc_vpt(ipt,ira)*bfng_val(ipt,mu,2)*
     &                bfn_val(ipt,nu)
                     scr1(3)=scr1(3)
     &               -wt(ipt)*xc_vpt(ipt,ira)*bfng_val(ipt,mu,3)*
     &                bfn_val(ipt,nu)
                     scr2(1)=scr2(1)
     &               -wt(ipt)*xc_vpt(ipt,irb)*bfng_val(ipt,mu,1)*
     &                bfn_val(ipt,nu)
                     scr2(2)=scr2(2)
     &               -wt(ipt)*xc_vpt(ipt,irb)*bfng_val(ipt,mu,2)*
     &                bfn_val(ipt,nu)
                     scr2(3)=scr2(3)
     &               -wt(ipt)*xc_vpt(ipt,irb)*bfng_val(ipt,mu,3)*
     &                bfn_val(ipt,nu)
                  enddo
                  if (gradcorr_sw) then
                     do ipt=1,npts
                        gagmndbx = grho(ipt,1,1)*(
     &                           - bfng_val(ipt,mu,1)*bfng_val(ipt,nu,1)
     &                           - bfn_hess(ipt,mu,hxx)*bfn_val(ipt,nu))
     &                           + grho(ipt,1,2)*(
     &                           - bfng_val(ipt,mu,1)*bfng_val(ipt,nu,2)
     &                           - bfn_hess(ipt,mu,hxy)*bfn_val(ipt,nu))
     &                           + grho(ipt,1,3)*(
     &                           - bfng_val(ipt,mu,1)*bfng_val(ipt,nu,3)
     &                           - bfn_hess(ipt,mu,hxz)*bfn_val(ipt,nu))
                        gbgmndbx = grho(ipt,2,1)*(
     &                           - bfng_val(ipt,mu,1)*bfng_val(ipt,nu,1)
     &                           - bfn_hess(ipt,mu,hxx)*bfn_val(ipt,nu))
     &                           + grho(ipt,2,2)*(
     &                           - bfng_val(ipt,mu,1)*bfng_val(ipt,nu,2)
     &                           - bfn_hess(ipt,mu,hxy)*bfn_val(ipt,nu))
     &                           + grho(ipt,2,3)*(
     &                           - bfng_val(ipt,mu,1)*bfng_val(ipt,nu,3)
     &                           - bfn_hess(ipt,mu,hxz)*bfn_val(ipt,nu))
                        gagmndby = grho(ipt,1,1)*(
     &                           - bfng_val(ipt,mu,2)*bfng_val(ipt,nu,1)
     &                           - bfn_hess(ipt,mu,hxy)*bfn_val(ipt,nu))
     &                           + grho(ipt,1,2)*(
     &                           - bfng_val(ipt,mu,2)*bfng_val(ipt,nu,2)
     &                           - bfn_hess(ipt,mu,hyy)*bfn_val(ipt,nu))
     &                           + grho(ipt,1,3)*(
     &                           - bfng_val(ipt,mu,2)*bfng_val(ipt,nu,3)
     &                           - bfn_hess(ipt,mu,hyz)*bfn_val(ipt,nu))
                        gbgmndby = grho(ipt,2,1)*(
     &                           - bfng_val(ipt,mu,2)*bfng_val(ipt,nu,1)
     &                           - bfn_hess(ipt,mu,hxy)*bfn_val(ipt,nu))
     &                           + grho(ipt,2,2)*(
     &                           - bfng_val(ipt,mu,2)*bfng_val(ipt,nu,2)
     &                           - bfn_hess(ipt,mu,hyy)*bfn_val(ipt,nu))
     &                           + grho(ipt,2,3)*(
     &                           - bfng_val(ipt,mu,2)*bfng_val(ipt,nu,3)
     &                           - bfn_hess(ipt,mu,hyz)*bfn_val(ipt,nu))
                        gagmndbz = grho(ipt,1,1)*(
     &                           - bfng_val(ipt,mu,3)*bfng_val(ipt,nu,1)
     &                           - bfn_hess(ipt,mu,hxz)*bfn_val(ipt,nu))
     &                           + grho(ipt,1,2)*(
     &                           - bfng_val(ipt,mu,3)*bfng_val(ipt,nu,2)
     &                           - bfn_hess(ipt,mu,hyz)*bfn_val(ipt,nu))
     &                           + grho(ipt,1,3)*(
     &                           - bfng_val(ipt,mu,3)*bfng_val(ipt,nu,3)
     &                           - bfn_hess(ipt,mu,hzz)*bfn_val(ipt,nu))
                        gbgmndbz = grho(ipt,2,1)*(
     &                           - bfng_val(ipt,mu,3)*bfng_val(ipt,nu,1)
     &                           - bfn_hess(ipt,mu,hxz)*bfn_val(ipt,nu))
     &                           + grho(ipt,2,2)*(
     &                           - bfng_val(ipt,mu,3)*bfng_val(ipt,nu,2)
     &                           - bfn_hess(ipt,mu,hyz)*bfn_val(ipt,nu))
     &                           + grho(ipt,2,3)*(
     &                           - bfng_val(ipt,mu,3)*bfng_val(ipt,nu,3)
     &                           - bfn_hess(ipt,mu,hzz)*bfn_val(ipt,nu))
                        scr1(1) = scr1(1)
     &                          + 2*wt(ipt)*xc_dvpt(ipt,igaa)*gagmndbx
     &                          +   wt(ipt)*xc_dvpt(ipt,igab)*gbgmndbx
                        scr1(2) = scr1(2)
     &                          + 2*wt(ipt)*xc_dvpt(ipt,igaa)*gagmndby
     &                          +   wt(ipt)*xc_dvpt(ipt,igab)*gbgmndby
                        scr1(3) = scr1(3)
     &                          + 2*wt(ipt)*xc_dvpt(ipt,igaa)*gagmndbz
     &                          +   wt(ipt)*xc_dvpt(ipt,igab)*gbgmndbz
                        scr2(1) = scr2(1)
     &                          + 2*wt(ipt)*xc_dvpt(ipt,igbb)*gbgmndbx
     &                          +   wt(ipt)*xc_dvpt(ipt,igab)*gagmndbx
                        scr2(2) = scr2(2)
     &                          + 2*wt(ipt)*xc_dvpt(ipt,igbb)*gbgmndby
     &                          +   wt(ipt)*xc_dvpt(ipt,igab)*gagmndby
                        scr2(3) = scr2(3)
     &                          + 2*wt(ipt)*xc_dvpt(ipt,igbb)*gbgmndbz
     &                          +   wt(ipt)*xc_dvpt(ipt,igab)*gagmndbz
                     enddo
                  endif
                  qa(iki+mu,ipert+1)=qa(iki+mu,ipert+1)+scr1(1)
                  qa(iki+mu,ipert+2)=qa(iki+mu,ipert+2)+scr1(2)
                  qa(iki+mu,ipert+3)=qa(iki+mu,ipert+3)+scr1(3)
                  if (gwt_sw) then
                     qa(iki+mu,jpert+1)=qa(iki+mu,jpert+1)-scr1(1)
                     qa(iki+mu,jpert+2)=qa(iki+mu,jpert+2)-scr1(2)
                     qa(iki+mu,jpert+3)=qa(iki+mu,jpert+3)-scr1(3)
                  endif
                  qb(iki+mu,ipert+1)=qb(iki+mu,ipert+1)+scr2(1)
                  qb(iki+mu,ipert+2)=qb(iki+mu,ipert+2)+scr2(2)
                  qb(iki+mu,ipert+3)=qb(iki+mu,ipert+3)+scr2(3)
                  if (gwt_sw) then
                     qb(iki+mu,jpert+1)=qb(iki+mu,jpert+1)-scr2(1)
                     qb(iki+mu,jpert+2)=qb(iki+mu,jpert+2)-scr2(2)
                     qb(iki+mu,jpert+3)=qb(iki+mu,jpert+3)-scr2(3)
                  endif
               enddo
            enddo
 100     continue
      endif ! rkstyp_sw
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine dksm_exp_dft_ao_scr(rkstyp_sw,gradcorr_sw,
     &     gwt_sw,gwt_avail_sw,
     &     active_bfn_list,active_bfn_indx,active_bfn_atms,
     &     n_active_bfn,n_active_atm,
     &     near_atom_list,num_near_atoms,
     &     nao,npts,npert,iatom,wt,gwt,bfn_val,bfng_val,
     &     bfn_hess,grho,drhodb,dgrhodb,
     &     xc_vpt,xc_dvpt,xc_hpt,xc_dhpt,
     &     qa,qb,mxp)
      implicit none
c
c     Calculates the DFT contributions to the derivative Fock matrices.
c
c     Parameters:
c
C *Parameter module for CCP1/DFT            
C *Memory information		
C *----------------
C * max_block 	- maximum number of blocks which can be allocate
      integer max_block
      parameter(max_block=20)
C *							
C *Basis sets information                              
C *----------------------				
C *max_tag 	- maximum number of basis sets which can be alloc
C *max_atype	- maximum number of centre types	
C *max_gtype	- maximum number of grid types (at least every element 
C                 has different grid type)
C *max_grids    - maximum number of grids (different terms may have
C                 different grids, i.e. CPKS equations use a different
C                 grid than the one used for the KS-matrix)
C *max_shel	- maximum number of shells on centre
C *max_prm	- maximum number of primitives for any given centre
C *maxL		- maximum angular momentum allowed	
C *max_func	- maximum number of basis functions for any centre
      integer max_tag,max_atype,max_shel,max_prm,max_ang
      integer max_gtype,max_grids
      parameter(max_tag=3,max_atype=30,max_shel=500,max_prm=5000)
      parameter(max_ang=5,max_gtype=10,max_grids=2)
C *								
C *Geometry information					
C *--------------------					
C *max_atom	- maximum number of atoms in system
      integer max_atom
      parameter(max_atom=750)
C *
C *Accuracy information				
C ---------------------			
C *global_accuracy 	- global accuracy	
      real*8  global_accuracy
      parameter(global_accuracy=1.0d-14)

C *Grid information
c *----------------
c *maxradzn - The maximum number of radial zones. 
      integer maxgpt,maxrad,maxfpt,maxang,maxradzn,maxtablerows
      parameter(maxgpt=2900,maxrad=50,maxang=302,maxfpt=100)
      parameter(maxradzn=35,maxtablerows=7)
c
c     The parameters defined here label the 1st and 2nd derivative
c     components of the density functional. Because a functional
c     may be purely local or of a GGA type the components are grouped
c     into derivative with respect to the density only and others.
c
c     1st derivatives:
c
      integer ira, irb
      integer igaa, igab, igbb
      integer ita, itb
c
c     2nd derivatives:
c
      integer irara, irarb, irbrb
      integer iragaa, iragab, iragbb, irbgaa, irbgab, irbgbb
      integer igaagaa, igaagab, igaagbb, igabgab, igabgbb, igbbgbb
      integer irata, iratb, irbta, irbtb
      integer igaata, igaatb, igabta, igabtb, igbbta, igbbtb
      integer itata, itatb, itbtb
c
      parameter (ira     = 1)
      parameter (irb     = 2)
c
      parameter (igaa    = 1)
      parameter (igab    = 2)
      parameter (igbb    = 3)
c
      parameter (ita     = 1)
      parameter (itb     = 2)
c
c
      parameter (irara   = 1)
      parameter (irarb   = 2)
      parameter (irbrb   = 3)
c
      parameter (iragaa  = 1)
      parameter (iragab  = 3)
      parameter (iragbb  = 4)
      parameter (irbgaa  = 5)
      parameter (irbgab  = 6)
      parameter (irbgbb  = 7)
      parameter (igaagaa = 2)
      parameter (igaagab = 8)
      parameter (igaagbb = 9)
      parameter (igabgab = 10)
      parameter (igabgbb = 11)
      parameter (igbbgbb = 12)
c
      parameter (irata   = 1)
      parameter (iratb   = 3)
      parameter (irbta   = 4)
      parameter (irbtb   = 5)
      parameter (itata   = 2)
      parameter (itatb   = 6)
      parameter (itbtb   = 7)
c
      parameter (igaata  = 1)
      parameter (igaatb  = 2)
      parameter (igabta  = 3)
      parameter (igabtb  = 4)
      parameter (igbbta  = 5)
      parameter (igbbtb  = 6)
c
c     The parameters defined here label the 1st, 2nd and 3rd derivative
c     components of the basis functions.
c
c     2nd derivatives:
c
      integer hxx, hyy, hzz, hxy, hxz, hyz
c
      parameter (hxx = 1)
      parameter (hyy = 2)
      parameter (hzz = 3)
      parameter (hxy = 4)
      parameter (hxz = 5)
      parameter (hyz = 6)
c
      integer hcc(3,3)
      data hcc/hxx,hxy,hxz,
     &         hxy,hyy,hyz,
     &         hxz,hyz,hzz/
c
c     3rd derivatives:
c
      integer txxx, txxy, txxz, txyy, txyz
      integer txzz, tyyy, tyyz, tyzz, tzzz
c
      parameter (txxx =  1)
      parameter (txxy =  2)
      parameter (txxz =  3)
      parameter (txyy =  4)
      parameter (txyz =  5)
      parameter (txzz =  6)
      parameter (tyyy =  7)
      parameter (tyyz =  8)
      parameter (tyzz =  9)
      parameter (tzzz = 10)
c
      integer tc3(3,3,3)
      data tc3/txxx,txxy,txxz,
     &         txxy,txyy,txyz,
     &         txxz,txyz,txzz,
     &         txxy,txyy,txyz,
     &         txyy,tyyy,tyyz,
     &         txyz,tyyz,tyzz,
     &         txxz,txyz,txzz,
     &         txyz,tyyz,tyzz,
     &         txzz,tyzz,tzzz/
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
c     In variables:
c
c
C *Basis sets include file
c
      real*8  alpha
      real*8  cont_coeff
      integer num_bset,bset_tags
      common/basis_sets/alpha(max_tag,max_atype,max_prm),
     &                  cont_coeff(max_tag,max_atype,max_prm,max_ang),
     &                  num_bset,bset_tags(max_tag)

      integer Ashl, Aprm, Abfn, totshl, totprm, totbfn, size_shlA
      integer size_basA, size_primA, maxi_shlA, maxi_basA, maxi_primA
      integer num_types, atom_tag , num_shl, atm_typ, nprim, angmom
      integer hybrid, pstart

      common/basis_size_info/Ashl(max_tag,max_atype),
     &                       Aprm(max_tag,max_atype),
     &                       Abfn(max_tag,max_atype),
     &                       totshl(max_tag),
     &                       totprm(max_tag),
     &                       totbfn(max_tag),
     &                       size_shlA,maxi_shlA,
     &                       size_basA,maxi_basA,
     &                       size_primA,maxi_primA
C
C Descriptions
C 
C  num_types		-	number of types of atoms for a given basis
C  atom_tag		-	list of atom type id numbers
C
      common/basis_cent_info/num_types(max_tag),
     &                       atom_tag(max_tag,max_atom)
C
C Descriptions
C
C  num_shl		-	number of shells for given atom type
C  atm_typ		-	atomic number of atom type
C  nprim		-  	number of primitives 
C  angmom		-	angular momentum
C  hybrid		-	level of hybrid shell. Is same as angmom if
C				shell is not a hybrid. For hybrid shells is
C				always less than angmom e.g. for an sp shell
C				angmom=2, hybrid=1
C  pstart		-	start of exponents and contraction coeffs 
C				contained in include basis.hf77, for a given
C				shell
C
      common/basis_cont_info/num_shl(max_tag,max_atype),
     &                       atm_typ(max_tag,max_atype),
     &                       nprim(max_tag,max_atype,max_shel),
     &                       angmom(max_tag,max_atype,max_shel),
     &                       hybrid(max_tag,max_atype,max_shel),
     &                       pstart(max_tag,max_atype,max_shel)

c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      logical rkstyp_sw   ! .true. if closed shell
      logical gradcorr_sw
      logical gwt_sw       ! .true. if gradients of quadrature used
      logical gwt_avail_sw ! .true. if gradients of quadrature present
      integer n_active_bfn
      integer n_active_atm
      integer active_bfn_list(n_active_bfn)
      integer active_bfn_indx(n_active_bfn)
      integer active_bfn_atms(n_active_atm)
      integer num_near_atoms
      integer near_atom_list(num_near_atoms)
      integer nao, npts, npert, mxp
      integer iatom
      real*8 wt(mxp)
      real*8 gwt(3,mxp,*)
      real*8 bfn_val(mxp,nao)
      real*8 bfng_val(mxp,nao,3)
      real*8 bfn_hess(mxp,nao,6)
      real*8 grho(mxp,2,3)
      real*8 drhodb(mxp,2,3*n_active_atm)    ! explicit derivatives only
      real*8 dgrhodb(mxp,2,3,3*n_active_atm) ! explicit derivatives only
      real*8 xc_vpt(mxp,2)
      real*8 xc_dvpt(mxp,3)
      real*8 xc_hpt(mxp,3)
      real*8 xc_dhpt(mxp,12)
c
c     Out variables:
c
      real*8 qa((nao+1)*nao/2,npert)
      real*8 qb((nao+1)*nao/2,npert)
c
c     Local variables:
c
      integer latm, ipert, jpert, lpert, i, j, k, ipt, n, nn
c     integer bstart
      integer bend, mu, nu
      real*8 scr1(3), scr2(3)
      real*8 t1, t2, t3, t
      integer iki
c
      real*8 gagmn, gbgmn
      real*8 gadbxgmn, gadbygmn, gadbzgmn
      real*8 gbdbxgmn, gbdbygmn, gbdbzgmn
      real*8 gaadbx, gabdbx, gbbdbx
      real*8 gaadby, gabdby, gbbdby
      real*8 gaadbz, gbbdbz, gabdbz
      real*8 gagmndbx, gagmndby, gagmndbz
      real*8 gbgmndbx, gbgmndby, gbgmndbz
c
      real*8 gadbgmn(3)
      real*8 gaadb(3)
c     real*8 gagmndb(3)
c
      integer ao_tag
      parameter(ao_tag=1)
c
c     Code:
c
c     For the time being I don't want to worry about using translational
c     invariance.
c     if (npert.ne.3*ngridcentres) call caserr('not implemented yet')
c
cDEBUG
c     if (gwt_sw) then
c     write(*,*)'*** WE should never get here 1'
c     call caserr('problem 1')
c     endif
c     write(*,*)'*** rhs_dft grho   =',grho(3,1,1)
c     write(*,*)'*** rhs_dft drhodb =',drhodb(3,1,1)
c     write(*,*)'*** rhs_dft dgrhodb=',dgrhodb(3,1,1,1)
c     write(*,*)'*** dksm_exp_dft_ao_scr ',loc(rkstyp_sw)
cDEBUG
      if (rkstyp_sw) then
c
c        the gradient of the weights contributions
c
         if (gwt_sw.and.gwt_avail_sw) then
            do latm=1,num_near_atoms
               ipert=3*(near_atom_list(latm)-1)
               do i=1,n_active_bfn
                  nn = iky(active_bfn_list(i))
                  do j=1,i
                     n=nn+active_bfn_list(j)
                     do ipt=1,npts
                        do k = 1, 3
                           qa(n,ipert+k)=qa(n,ipert+k) +
     &                       gwt(k,ipt,latm)*xc_vpt(ipt,ira)*
     &                       bfn_val(ipt,i)*bfn_val(ipt,j)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            if (gradcorr_sw) then
               do latm=1,num_near_atoms
                  ipert=3*(near_atom_list(latm)-1)
                  do i=1,n_active_bfn
                     nn = iky(active_bfn_list(i))
                     do j=1,i
                        n=nn+active_bfn_list(j)
                        do ipt=1,npts
                           t1 = xc_dvpt(ipt,igaa)*grho(ipt,1,1)
                           t2 = xc_dvpt(ipt,igaa)*grho(ipt,1,2)
                           t3 = xc_dvpt(ipt,igaa)*grho(ipt,1,3)
                           t  = (bfng_val(ipt,i,1)*t1
     &                        +  bfng_val(ipt,i,2)*t2
     &                        +  bfng_val(ipt,i,3)*t3)*bfn_val(ipt,j)
     &                        + (bfng_val(ipt,j,1)*t1
     &                        +  bfng_val(ipt,j,2)*t2
     &                        +  bfng_val(ipt,j,3)*t3)*bfn_val(ipt,i)
                           do k = 1, 3
                              qa(n,ipert+k)=qa(n,ipert+k) +
     &                          0.5d0*gwt(k,ipt,latm)*t
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            endif ! gradcorr_sw
         endif ! gwt_sw
c
         jpert = 3*(iatom-1)
         bend = 0
         do 100 latm=1,n_active_atm
            lpert=3*(latm-1)
            ipert=3*(active_bfn_atms(latm)-1)
            if (ipert.eq.jpert.and.gwt_sw) goto 100
            do mu=1,n_active_bfn
               nn=iky(active_bfn_list(mu))
               do nu=1,mu
                  n=nn+active_bfn_list(nu)
                  call aclear_dp(scr1,3,0.0d0)
                  do k=1,3
                     do ipt=1,npts
                        scr1(k)=scr1(k) 
     &                  +wt(ipt)*xc_hpt(ipt,irara)*bfn_val(ipt,mu)*
     &                   bfn_val(ipt,nu)*drhodb(ipt,1,lpert+k)
                     enddo ! ipt
                  enddo ! k
                  qa(n,ipert+1)=qa(n,ipert+1)+scr1(1)
                  qa(n,ipert+2)=qa(n,ipert+2)+scr1(2)
                  qa(n,ipert+3)=qa(n,ipert+3)+scr1(3)
                  if (gwt_sw) then
                    qa(n,jpert+1)=qa(n,jpert+1)-scr1(1)
                    qa(n,jpert+2)=qa(n,jpert+2)-scr1(2)
                    qa(n,jpert+3)=qa(n,jpert+3)-scr1(3)
                  endif ! gwt_sw
 
                  if (gradcorr_sw) then
                     call aclear_dp(scr1,3,0.0d0)
                     do ipt=1,npts
                        gagmn = grho(ipt,1,1)*(
     &                            bfng_val(ipt,mu,1)*bfn_val(ipt,nu)
     &                          + bfng_val(ipt,nu,1)*bfn_val(ipt,mu))
     &                        + grho(ipt,1,2)*(
     &                            bfng_val(ipt,mu,2)*bfn_val(ipt,nu)
     &                          + bfng_val(ipt,nu,2)*bfn_val(ipt,mu))
     &                        + grho(ipt,1,3)*(
     &                            bfng_val(ipt,mu,3)*bfn_val(ipt,nu)
     &                          + bfng_val(ipt,nu,3)*bfn_val(ipt,mu))
                        gagmn = 0.5d0*gagmn
                        do k=1,3
                           gadbgmn(k) 
     &                           = dgrhodb(ipt,1,1,lpert+k)*(
     &                               bfng_val(ipt,mu,1)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,1)*bfn_val(ipt,mu))
     &                           + dgrhodb(ipt,1,2,lpert+k)*(
     &                               bfng_val(ipt,mu,2)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,2)*bfn_val(ipt,mu))
     &                           + dgrhodb(ipt,1,3,lpert+k)*(
     &                               bfng_val(ipt,mu,3)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,3)*bfn_val(ipt,mu))
                           gaadb(k) = 2*(
     &                           dgrhodb(ipt,1,1,lpert+k)*grho(ipt,1,1)
     &                         + dgrhodb(ipt,1,2,lpert+k)*grho(ipt,1,2)
     &                         + dgrhodb(ipt,1,3,lpert+k)*grho(ipt,1,3))
                           gaadb(k) = 0.5d0*gaadb(k)
                        enddo ! k
                        do k=1,3
                        scr1(k)=scr1(k) 
     &                  +2*wt(ipt)*xc_dhpt(ipt,iragaa)*gagmn*
     &                     drhodb(ipt,1,lpert+k)
     &                  +  wt(ipt)*xc_dhpt(ipt,iragaa)*gaadb(k)*
     &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
     &                  +  wt(ipt)*xc_dhpt(ipt,igaagaa)*gagmn*gaadb(k)
     &                  +  wt(ipt)*xc_dvpt(ipt,igaa)*gadbgmn(k)
                        enddo ! k
                     enddo ! ipt
                     qa(n,ipert+1)=qa(n,ipert+1)+scr1(1)
                     qa(n,ipert+2)=qa(n,ipert+2)+scr1(2)
                     qa(n,ipert+3)=qa(n,ipert+3)+scr1(3)
                     if (gwt_sw) then
                        qa(n,jpert+1)=qa(n,jpert+1)-scr1(1)
                        qa(n,jpert+2)=qa(n,jpert+2)-scr1(2)
                        qa(n,jpert+3)=qa(n,jpert+3)-scr1(3)
                     endif ! gwt_sw
                  endif ! gradcorr_sw
               enddo ! nu
            enddo ! mu
 100     continue
c
         do 110 mu=1,n_active_bfn
c           lpert=3*(active_bfn_indx(mu)-1)
            ipert=3*(active_bfn_atms(active_bfn_indx(mu))-1)
            if (ipert.eq.jpert.and.gwt_sw) goto 110
            iki=iky(active_bfn_list(mu))
            do nu=1, mu
               n=iki+active_bfn_list(nu)
c
               call aclear_dp(scr1,3,0.0d0)
               do k=1,3
                  do ipt=1,npts
                     scr1(k)=scr1(k)
     &               -wt(ipt)*xc_vpt(ipt,ira)*bfng_val(ipt,mu,k)*
     &                bfn_val(ipt,nu)
                  enddo ! ipt
               enddo ! k
               if (gradcorr_sw) then
                  do ipt=1,npts
                     gagmndbx = grho(ipt,1,1)*(
     &                        - bfng_val(ipt,mu,1)*bfng_val(ipt,nu,1)
     &                        - bfn_hess(ipt,mu,hxx)*bfn_val(ipt,nu))
     &                        + grho(ipt,1,2)*(
     &                        - bfng_val(ipt,mu,1)*bfng_val(ipt,nu,2)
     &                        - bfn_hess(ipt,mu,hxy)*bfn_val(ipt,nu))
     &                        + grho(ipt,1,3)*(
     &                        - bfng_val(ipt,mu,1)*bfng_val(ipt,nu,3)
     &                        - bfn_hess(ipt,mu,hxz)*bfn_val(ipt,nu))
                     gagmndby = grho(ipt,1,1)*(
     &                        - bfng_val(ipt,mu,2)*bfng_val(ipt,nu,1)
     &                        - bfn_hess(ipt,mu,hxy)*bfn_val(ipt,nu))
     &                        + grho(ipt,1,2)*(
     &                        - bfng_val(ipt,mu,2)*bfng_val(ipt,nu,2)
     &                        - bfn_hess(ipt,mu,hyy)*bfn_val(ipt,nu))
     &                        + grho(ipt,1,3)*(
     &                        - bfng_val(ipt,mu,2)*bfng_val(ipt,nu,3)
     &                        - bfn_hess(ipt,mu,hyz)*bfn_val(ipt,nu))
                     gagmndbz = grho(ipt,1,1)*(
     &                        - bfng_val(ipt,mu,3)*bfng_val(ipt,nu,1)
     &                        - bfn_hess(ipt,mu,hxz)*bfn_val(ipt,nu))
     &                        + grho(ipt,1,2)*(
     &                        - bfng_val(ipt,mu,3)*bfng_val(ipt,nu,2)
     &                        - bfn_hess(ipt,mu,hyz)*bfn_val(ipt,nu))
     &                        + grho(ipt,1,3)*(
     &                        - bfng_val(ipt,mu,3)*bfng_val(ipt,nu,3)
     &                        - bfn_hess(ipt,mu,hzz)*bfn_val(ipt,nu))
                     gagmndbx = 0.5d0*gagmndbx
                     gagmndby = 0.5d0*gagmndby
                     gagmndbz = 0.5d0*gagmndbz
                     scr1(1) = scr1(1)
     &                       +   wt(ipt)*xc_dvpt(ipt,igaa)*gagmndbx
                     scr1(2) = scr1(2)
     &                       +   wt(ipt)*xc_dvpt(ipt,igaa)*gagmndby
                     scr1(3) = scr1(3)
     &                       +   wt(ipt)*xc_dvpt(ipt,igaa)*gagmndbz
                  enddo ! ipt
               endif ! gradcorr_sw
               qa(n,ipert+1)=qa(n,ipert+1)+scr1(1)
               qa(n,ipert+2)=qa(n,ipert+2)+scr1(2)
               qa(n,ipert+3)=qa(n,ipert+3)+scr1(3)
               if (gwt_sw) then
                  qa(n,jpert+1)=qa(n,jpert+1)-scr1(1)
                  qa(n,jpert+2)=qa(n,jpert+2)-scr1(2)
                  qa(n,jpert+3)=qa(n,jpert+3)-scr1(3)
               endif ! gwt_sw
            enddo ! nu
            do nu=mu, n_active_bfn
               iki = iky(active_bfn_list(nu))
               n   = iki + active_bfn_list(mu)
               call aclear_dp(scr1,3,0.0d0)
               do k=1,3
                  do ipt=1,npts
                     scr1(k)=scr1(k)
     &               -wt(ipt)*xc_vpt(ipt,ira)*bfng_val(ipt,mu,k)*
     &                bfn_val(ipt,nu)
                  enddo ! ipt
               enddo ! k
               if (gradcorr_sw) then
                  do ipt=1,npts
                     gagmndbx = grho(ipt,1,1)*(
     &                        - bfng_val(ipt,mu,1)*bfng_val(ipt,nu,1)
     &                        - bfn_hess(ipt,mu,hxx)*bfn_val(ipt,nu))
     &                        + grho(ipt,1,2)*(
     &                        - bfng_val(ipt,mu,1)*bfng_val(ipt,nu,2)
     &                        - bfn_hess(ipt,mu,hxy)*bfn_val(ipt,nu))
     &                        + grho(ipt,1,3)*(
     &                        - bfng_val(ipt,mu,1)*bfng_val(ipt,nu,3)
     &                        - bfn_hess(ipt,mu,hxz)*bfn_val(ipt,nu))
                     gagmndby = grho(ipt,1,1)*(
     &                        - bfng_val(ipt,mu,2)*bfng_val(ipt,nu,1)
     &                        - bfn_hess(ipt,mu,hxy)*bfn_val(ipt,nu))
     &                        + grho(ipt,1,2)*(
     &                        - bfng_val(ipt,mu,2)*bfng_val(ipt,nu,2)
     &                        - bfn_hess(ipt,mu,hyy)*bfn_val(ipt,nu))
     &                        + grho(ipt,1,3)*(
     &                        - bfng_val(ipt,mu,2)*bfng_val(ipt,nu,3)
     &                        - bfn_hess(ipt,mu,hyz)*bfn_val(ipt,nu))
                     gagmndbz = grho(ipt,1,1)*(
     &                        - bfng_val(ipt,mu,3)*bfng_val(ipt,nu,1)
     &                        - bfn_hess(ipt,mu,hxz)*bfn_val(ipt,nu))
     &                        + grho(ipt,1,2)*(
     &                        - bfng_val(ipt,mu,3)*bfng_val(ipt,nu,2)
     &                        - bfn_hess(ipt,mu,hyz)*bfn_val(ipt,nu))
     &                        + grho(ipt,1,3)*(
     &                        - bfng_val(ipt,mu,3)*bfng_val(ipt,nu,3)
     &                        - bfn_hess(ipt,mu,hzz)*bfn_val(ipt,nu))
                     gagmndbx = 0.5d0*gagmndbx
                     gagmndby = 0.5d0*gagmndby
                     gagmndbz = 0.5d0*gagmndbz
                     scr1(1) = scr1(1)
     &                       +   wt(ipt)*xc_dvpt(ipt,igaa)*gagmndbx
                     scr1(2) = scr1(2)
     &                       +   wt(ipt)*xc_dvpt(ipt,igaa)*gagmndby
                     scr1(3) = scr1(3)
     &                       +   wt(ipt)*xc_dvpt(ipt,igaa)*gagmndbz
                  enddo ! ipt
               endif ! gradcorr_sw
               qa(n,ipert+1)=qa(n,ipert+1)+scr1(1)
               qa(n,ipert+2)=qa(n,ipert+2)+scr1(2)
               qa(n,ipert+3)=qa(n,ipert+3)+scr1(3)
               if (gwt_sw) then
                  qa(n,jpert+1)=qa(n,jpert+1)-scr1(1)
                  qa(n,jpert+2)=qa(n,jpert+2)-scr1(2)
                  qa(n,jpert+3)=qa(n,jpert+3)-scr1(3)
               endif ! gwt_sw
            enddo ! nu
 110     continue ! mu
c
      else ! rkstyp_sw
c
c        the gradient of the weights contributions
c
         if (gwt_sw.and.gwt_avail_sw) then
            do latm=1,num_near_atoms
               ipert=3*(near_atom_list(latm)-1)
               do i=1,n_active_bfn
                  nn = iky(active_bfn_list(i))
                  do j=1,i
                     n=nn+active_bfn_list(j)
                     do ipt=1,npts
                        do k = 1, 3
                           qa(n,ipert+k)=qa(n,ipert+k) +
     &                       gwt(k,ipt,latm)*xc_vpt(ipt,ira)*
     &                       bfn_val(ipt,i)*bfn_val(ipt,j)
                        enddo
                     enddo
                     do ipt=1,npts
                        do k = 1, 3
                           qb(n,ipert+k)=qb(n,ipert+k) +
     &                       gwt(k,ipt,latm)*xc_vpt(ipt,irb)*
     &                       bfn_val(ipt,i)*bfn_val(ipt,j)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            if (gradcorr_sw) then
               do latm=1,num_near_atoms
                  ipert=3*(near_atom_list(latm)-1)
                  do i=1,n_active_bfn
                     nn = iky(active_bfn_list(i))
                     do j=1,i
                        n=nn+active_bfn_list(j)
                        do ipt=1,npts
                           t1 = (2.0d0*xc_dvpt(ipt,igaa)*grho(ipt,1,1)
     &                        +        xc_dvpt(ipt,igab)*grho(ipt,2,1))
                           t2 = (2.0d0*xc_dvpt(ipt,igaa)*grho(ipt,1,2)
     &                        +        xc_dvpt(ipt,igab)*grho(ipt,2,2))
                           t3 = (2.0d0*xc_dvpt(ipt,igaa)*grho(ipt,1,3)
     &                        +        xc_dvpt(ipt,igab)*grho(ipt,2,3))
                           t  = (bfng_val(ipt,i,1)*t1
     &                        +  bfng_val(ipt,i,2)*t2
     &                        +  bfng_val(ipt,i,3)*t3)*bfn_val(ipt,j)
     &                        + (bfng_val(ipt,j,1)*t1
     &                        +  bfng_val(ipt,j,2)*t2
     &                        +  bfng_val(ipt,j,3)*t3)*bfn_val(ipt,i)
                           do k = 1, 3
                              qa(n,ipert+k)=qa(n,ipert+k) +
     &                          gwt(k,ipt,latm)*t
                           enddo
                        enddo
                     enddo
                  enddo
                  do i=1,n_active_bfn
                     nn = iky(active_bfn_list(i))
                     do j=1,i
                        n=nn+active_bfn_list(j)
                        do ipt=1,npts
                           t1 = (2.0d0*xc_dvpt(ipt,igbb)*grho(ipt,2,1)
     &                        +        xc_dvpt(ipt,igab)*grho(ipt,1,1))
                           t2 = (2.0d0*xc_dvpt(ipt,igbb)*grho(ipt,2,2)
     &                        +        xc_dvpt(ipt,igab)*grho(ipt,1,2))
                           t3 = (2.0d0*xc_dvpt(ipt,igbb)*grho(ipt,2,3)
     &                        +        xc_dvpt(ipt,igab)*grho(ipt,1,3))
                           t  = (bfng_val(ipt,i,1)*t1
     &                        +  bfng_val(ipt,i,2)*t2
     &                        +  bfng_val(ipt,i,3)*t3)*bfn_val(ipt,j)
     &                        + (bfng_val(ipt,j,1)*t1
     &                        +  bfng_val(ipt,j,2)*t2
     &                        +  bfng_val(ipt,j,3)*t3)*bfn_val(ipt,i)
                           do k = 1, 3
                              qb(n,ipert+k)=qb(n,ipert+k) +
     &                          gwt(k,ipt,latm)*t
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            endif ! gradcorr_sw
         endif ! gwt_sw
c
         jpert = 3*(iatom-1)
         bend = 0
         do 200 latm=1,n_active_atm
            lpert=3*(latm-1)
            ipert=3*(active_bfn_atms(latm)-1)
            if (ipert.eq.jpert.and.gwt_sw) goto 200
            do mu=1,n_active_bfn
               nn=iky(active_bfn_list(mu))
               do nu=1,mu
                  n=nn+active_bfn_list(nu)
                  call aclear_dp(scr1,3,0.0d0)
                  do ipt=1,npts
                     scr1(1)=scr1(1) 
     &               +wt(ipt)*xc_hpt(ipt,irara)*bfn_val(ipt,mu)*
     &                bfn_val(ipt,nu)*drhodb(ipt,1,lpert+1)
     &               +wt(ipt)*xc_hpt(ipt,irarb)*bfn_val(ipt,mu)*
     &                bfn_val(ipt,nu)*drhodb(ipt,2,lpert+1)
                     scr1(2)=scr1(2) 
     &               +wt(ipt)*xc_hpt(ipt,irara)*bfn_val(ipt,mu)*
     &                bfn_val(ipt,nu)*drhodb(ipt,1,lpert+2)
     &               +wt(ipt)*xc_hpt(ipt,irarb)*bfn_val(ipt,mu)*
     &                bfn_val(ipt,nu)*drhodb(ipt,2,lpert+2)
                     scr1(3)=scr1(3) 
     &               +wt(ipt)*xc_hpt(ipt,irara)*bfn_val(ipt,mu)*
     &                bfn_val(ipt,nu)*drhodb(ipt,1,lpert+3)
     &               +wt(ipt)*xc_hpt(ipt,irarb)*bfn_val(ipt,mu)*
     &                bfn_val(ipt,nu)*drhodb(ipt,2,lpert+3)
                  enddo
                  qa(n,ipert+1)=qa(n,ipert+1)+scr1(1)
                  qa(n,ipert+2)=qa(n,ipert+2)+scr1(2)
                  qa(n,ipert+3)=qa(n,ipert+3)+scr1(3)
                  if (gwt_sw) then
                    qa(n,jpert+1)=qa(n,jpert+1)-scr1(1)
                    qa(n,jpert+2)=qa(n,jpert+2)-scr1(2)
                    qa(n,jpert+3)=qa(n,jpert+3)-scr1(3)
                  endif
 
                  call aclear_dp(scr1,3,0.0d0)
                  do ipt=1,npts
                     scr1(1)=scr1(1) 
     &               +wt(ipt)*xc_hpt(ipt,irbrb)*bfn_val(ipt,mu)*
     &                bfn_val(ipt,nu)*drhodb(ipt,2,lpert+1)
     &               +wt(ipt)*xc_hpt(ipt,irarb)*bfn_val(ipt,mu)*
     &                bfn_val(ipt,nu)*drhodb(ipt,1,lpert+1)
                     scr1(2)=scr1(2) 
     &               +wt(ipt)*xc_hpt(ipt,irbrb)*bfn_val(ipt,mu)*
     &                bfn_val(ipt,nu)*drhodb(ipt,2,lpert+2)
     &               +wt(ipt)*xc_hpt(ipt,irarb)*bfn_val(ipt,mu)*
     &                bfn_val(ipt,nu)*drhodb(ipt,1,lpert+2)
                     scr1(3)=scr1(3) 
     &               +wt(ipt)*xc_hpt(ipt,irbrb)*bfn_val(ipt,mu)*
     &                bfn_val(ipt,nu)*drhodb(ipt,2,lpert+3)
     &               +wt(ipt)*xc_hpt(ipt,irarb)*bfn_val(ipt,mu)*
     &                bfn_val(ipt,nu)*drhodb(ipt,1,lpert+3)
                  enddo
                  qb(n,ipert+1)=qb(n,ipert+1)+scr1(1)
                  qb(n,ipert+2)=qb(n,ipert+2)+scr1(2)
                  qb(n,ipert+3)=qb(n,ipert+3)+scr1(3)
                  if (gwt_sw) then
                     qb(n,jpert+1)=qb(n,jpert+1)-scr1(1)
                     qb(n,jpert+2)=qb(n,jpert+2)-scr1(2)
                     qb(n,jpert+3)=qb(n,jpert+3)-scr1(3)
                  endif

                  if (gradcorr_sw) then
                     call aclear_dp(scr1,3,0.0d0)
                     call aclear_dp(scr2,3,0.0d0)
                     do ipt=1,npts
                        gagmn = grho(ipt,1,1)*(
     &                            bfng_val(ipt,mu,1)*bfn_val(ipt,nu)
     &                          + bfng_val(ipt,nu,1)*bfn_val(ipt,mu))
     &                        + grho(ipt,1,2)*(
     &                            bfng_val(ipt,mu,2)*bfn_val(ipt,nu)
     &                          + bfng_val(ipt,nu,2)*bfn_val(ipt,mu))
     &                        + grho(ipt,1,3)*(
     &                            bfng_val(ipt,mu,3)*bfn_val(ipt,nu)
     &                          + bfng_val(ipt,nu,3)*bfn_val(ipt,mu))
                        gbgmn = grho(ipt,2,1)*(
     &                            bfng_val(ipt,mu,1)*bfn_val(ipt,nu)
     &                          + bfng_val(ipt,nu,1)*bfn_val(ipt,mu))
     &                        + grho(ipt,2,2)*(
     &                            bfng_val(ipt,mu,2)*bfn_val(ipt,nu)
     &                          + bfng_val(ipt,nu,2)*bfn_val(ipt,mu))
     &                        + grho(ipt,2,3)*(
     &                            bfng_val(ipt,mu,3)*bfn_val(ipt,nu)
     &                          + bfng_val(ipt,nu,3)*bfn_val(ipt,mu))
                        gadbxgmn = dgrhodb(ipt,1,1,lpert+1)*(
     &                               bfng_val(ipt,mu,1)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,1)*bfn_val(ipt,mu))
     &                           + dgrhodb(ipt,1,2,ipert+1)*(
     &                               bfng_val(ipt,mu,2)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,2)*bfn_val(ipt,mu))
     &                           + dgrhodb(ipt,1,3,lpert+1)*(
     &                               bfng_val(ipt,mu,3)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,3)*bfn_val(ipt,mu))
                        gbdbxgmn = dgrhodb(ipt,2,1,lpert+1)*(
     &                               bfng_val(ipt,mu,1)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,1)*bfn_val(ipt,mu))
     &                           + dgrhodb(ipt,2,2,lpert+1)*(
     &                               bfng_val(ipt,mu,2)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,2)*bfn_val(ipt,mu))
     &                           + dgrhodb(ipt,2,3,lpert+1)*(
     &                               bfng_val(ipt,mu,3)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,3)*bfn_val(ipt,mu))
                        gadbygmn = dgrhodb(ipt,1,1,lpert+2)*(
     &                               bfng_val(ipt,mu,1)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,1)*bfn_val(ipt,mu))
     &                           + dgrhodb(ipt,1,2,lpert+2)*(
     &                               bfng_val(ipt,mu,2)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,2)*bfn_val(ipt,mu))
     &                           + dgrhodb(ipt,1,3,lpert+2)*(
     &                               bfng_val(ipt,mu,3)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,3)*bfn_val(ipt,mu))
                        gbdbygmn = dgrhodb(ipt,2,1,lpert+2)*(
     &                               bfng_val(ipt,mu,1)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,1)*bfn_val(ipt,mu))
     &                           + dgrhodb(ipt,2,2,lpert+2)*(
     &                               bfng_val(ipt,mu,2)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,2)*bfn_val(ipt,mu))
     &                           + dgrhodb(ipt,2,3,lpert+2)*(
     &                               bfng_val(ipt,mu,3)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,3)*bfn_val(ipt,mu))
                        gadbzgmn = dgrhodb(ipt,1,1,lpert+3)*(
     &                               bfng_val(ipt,mu,1)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,1)*bfn_val(ipt,mu))
     &                           + dgrhodb(ipt,1,2,lpert+3)*(
     &                               bfng_val(ipt,mu,2)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,2)*bfn_val(ipt,mu))
     &                           + dgrhodb(ipt,1,3,lpert+3)*(
     &                               bfng_val(ipt,mu,3)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,3)*bfn_val(ipt,mu))
                        gbdbzgmn = dgrhodb(ipt,2,1,lpert+3)*(
     &                               bfng_val(ipt,mu,1)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,1)*bfn_val(ipt,mu))
     &                           + dgrhodb(ipt,2,2,lpert+3)*(
     &                               bfng_val(ipt,mu,2)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,2)*bfn_val(ipt,mu))
     &                           + dgrhodb(ipt,2,3,lpert+3)*(
     &                               bfng_val(ipt,mu,3)*bfn_val(ipt,nu)
     &                             + bfng_val(ipt,nu,3)*bfn_val(ipt,mu))
                        gaadbx = 2*(
     &                           dgrhodb(ipt,1,1,lpert+1)*grho(ipt,1,1)
     &                         + dgrhodb(ipt,1,2,lpert+1)*grho(ipt,1,2)
     &                         + dgrhodb(ipt,1,3,lpert+1)*grho(ipt,1,3))
                        gbbdbx = 2*(
     &                           dgrhodb(ipt,2,1,lpert+1)*grho(ipt,2,1)
     &                         + dgrhodb(ipt,2,2,lpert+1)*grho(ipt,2,2)
     &                         + dgrhodb(ipt,2,3,lpert+1)*grho(ipt,2,3))
                        gabdbx = dgrhodb(ipt,1,1,lpert+1)*grho(ipt,2,1)
     &                         + dgrhodb(ipt,1,2,lpert+1)*grho(ipt,2,2)
     &                         + dgrhodb(ipt,1,3,lpert+1)*grho(ipt,2,3)
     &                         + dgrhodb(ipt,2,1,lpert+1)*grho(ipt,1,1)
     &                         + dgrhodb(ipt,2,2,lpert+1)*grho(ipt,1,2)
     &                         + dgrhodb(ipt,2,3,lpert+1)*grho(ipt,1,3)
                        gaadby = 2*(
     &                           dgrhodb(ipt,1,1,lpert+2)*grho(ipt,1,1)
     &                         + dgrhodb(ipt,1,2,lpert+2)*grho(ipt,1,2)
     &                         + dgrhodb(ipt,1,3,lpert+2)*grho(ipt,1,3))
                        gbbdby = 2*(
     &                           dgrhodb(ipt,2,1,lpert+2)*grho(ipt,2,1)
     &                         + dgrhodb(ipt,2,2,lpert+2)*grho(ipt,2,2)
     &                         + dgrhodb(ipt,2,3,lpert+2)*grho(ipt,2,3))
                        gabdby = dgrhodb(ipt,1,1,lpert+2)*grho(ipt,2,1)
     &                         + dgrhodb(ipt,1,2,lpert+2)*grho(ipt,2,2)
     &                         + dgrhodb(ipt,1,3,lpert+2)*grho(ipt,2,3)
     &                         + dgrhodb(ipt,2,1,lpert+2)*grho(ipt,1,1)
     &                         + dgrhodb(ipt,2,2,lpert+2)*grho(ipt,1,2)
     &                         + dgrhodb(ipt,2,3,lpert+2)*grho(ipt,1,3)
                        gaadbz = 2*(
     &                           dgrhodb(ipt,1,1,lpert+3)*grho(ipt,1,1)
     &                         + dgrhodb(ipt,1,2,lpert+3)*grho(ipt,1,2)
     &                         + dgrhodb(ipt,1,3,lpert+3)*grho(ipt,1,3))
                        gbbdbz = 2*(
     &                           dgrhodb(ipt,2,1,lpert+3)*grho(ipt,2,1)
     &                         + dgrhodb(ipt,2,2,lpert+3)*grho(ipt,2,2)
     &                         + dgrhodb(ipt,2,3,lpert+3)*grho(ipt,2,3))
                        gabdbz = dgrhodb(ipt,1,1,lpert+3)*grho(ipt,2,1)
     &                         + dgrhodb(ipt,1,2,lpert+3)*grho(ipt,2,2)
     &                         + dgrhodb(ipt,1,3,lpert+3)*grho(ipt,2,3)
     &                         + dgrhodb(ipt,2,1,lpert+3)*grho(ipt,1,1)
     &                         + dgrhodb(ipt,2,2,lpert+3)*grho(ipt,1,2)
     &                         + dgrhodb(ipt,2,3,lpert+3)*grho(ipt,1,3)
                        scr1(1)=scr1(1) 
     &                  +2*wt(ipt)*xc_dhpt(ipt,iragaa)*gagmn*
     &                     drhodb(ipt,1,lpert+1)
     &                  +  wt(ipt)*xc_dhpt(ipt,iragab)*gbgmn*
     &                     drhodb(ipt,1,lpert+1)
     &                  +2*wt(ipt)*xc_dhpt(ipt,irbgaa)*gagmn*
     &                     drhodb(ipt,2,lpert+1)
     &                  +  wt(ipt)*xc_dhpt(ipt,irbgab)*gbgmn*
     &                     drhodb(ipt,2,lpert+1)
     &                  +  wt(ipt)*xc_dhpt(ipt,iragaa)*gaadbx*
     &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
     &                  +  wt(ipt)*xc_dhpt(ipt,iragab)*gabdbx*
     &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
     &                  +  wt(ipt)*xc_dhpt(ipt,iragbb)*gbbdbx*
     &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
     &                  +2*wt(ipt)*xc_dhpt(ipt,igaagaa)*gagmn*gaadbx
     &                  +  wt(ipt)*xc_dhpt(ipt,igaagab)*gbgmn*gaadbx
     &                  +2*wt(ipt)*xc_dhpt(ipt,igaagab)*gagmn*gabdbx
     &                  +  wt(ipt)*xc_dhpt(ipt,igabgab)*gbgmn*gabdbx
     &                  +2*wt(ipt)*xc_dhpt(ipt,igaagbb)*gagmn*gbbdbx
     &                  +  wt(ipt)*xc_dhpt(ipt,igabgbb)*gbgmn*gbbdbx
     &                  +2*wt(ipt)*xc_dvpt(ipt,igaa)*gadbxgmn
     &                  +  wt(ipt)*xc_dvpt(ipt,igab)*gbdbxgmn
                        scr1(2)=scr1(2) 
     &                  +2*wt(ipt)*xc_dhpt(ipt,iragaa)*gagmn*
     &                     drhodb(ipt,1,lpert+2)
     &                  +  wt(ipt)*xc_dhpt(ipt,iragab)*gbgmn*
     &                     drhodb(ipt,1,lpert+2)
     &                  +2*wt(ipt)*xc_dhpt(ipt,irbgaa)*gagmn*
     &                     drhodb(ipt,2,lpert+2)
     &                  +  wt(ipt)*xc_dhpt(ipt,irbgab)*gbgmn*
     &                     drhodb(ipt,2,lpert+2)
     &                  +  wt(ipt)*xc_dhpt(ipt,iragaa)*gaadby*
     &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
     &                  +  wt(ipt)*xc_dhpt(ipt,iragab)*gabdby*
     &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
     &                  +  wt(ipt)*xc_dhpt(ipt,iragbb)*gbbdby*
     &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
     &                  +2*wt(ipt)*xc_dhpt(ipt,igaagaa)*gagmn*gaadby
     &                  +  wt(ipt)*xc_dhpt(ipt,igaagab)*gbgmn*gaadby
     &                  +2*wt(ipt)*xc_dhpt(ipt,igaagab)*gagmn*gabdby
     &                  +  wt(ipt)*xc_dhpt(ipt,igabgab)*gbgmn*gabdby
     &                  +2*wt(ipt)*xc_dhpt(ipt,igaagbb)*gagmn*gbbdby
     &                  +  wt(ipt)*xc_dhpt(ipt,igabgbb)*gbgmn*gbbdby
     &                  +2*wt(ipt)*xc_dvpt(ipt,igaa)*gadbygmn
     &                  +  wt(ipt)*xc_dvpt(ipt,igab)*gbdbygmn
                        scr1(3)=scr1(3) 
     &                  +2*wt(ipt)*xc_dhpt(ipt,iragaa)*gagmn*
     &                     drhodb(ipt,1,lpert+3)
     &                  +  wt(ipt)*xc_dhpt(ipt,iragab)*gbgmn*
     &                     drhodb(ipt,1,lpert+3)
     &                  +2*wt(ipt)*xc_dhpt(ipt,irbgaa)*gagmn*
     &                     drhodb(ipt,2,lpert+3)
     &                  +  wt(ipt)*xc_dhpt(ipt,irbgab)*gbgmn*
     &                     drhodb(ipt,2,lpert+3)
     &                  +  wt(ipt)*xc_dhpt(ipt,iragaa)*gaadbz*
     &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
     &                  +  wt(ipt)*xc_dhpt(ipt,iragab)*gabdbz*
     &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
     &                  +  wt(ipt)*xc_dhpt(ipt,iragbb)*gbbdbz*
     &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
     &                  +2*wt(ipt)*xc_dhpt(ipt,igaagaa)*gagmn*gaadbz
     &                  +  wt(ipt)*xc_dhpt(ipt,igaagab)*gbgmn*gaadbz
     &                  +2*wt(ipt)*xc_dhpt(ipt,igaagab)*gagmn*gabdbz
     &                  +  wt(ipt)*xc_dhpt(ipt,igabgab)*gbgmn*gabdbz
     &                  +2*wt(ipt)*xc_dhpt(ipt,igaagbb)*gagmn*gbbdbz
     &                  +  wt(ipt)*xc_dhpt(ipt,igabgbb)*gbgmn*gbbdbz
     &                  +2*wt(ipt)*xc_dvpt(ipt,igaa)*gadbzgmn
     &                  +  wt(ipt)*xc_dvpt(ipt,igab)*gbdbzgmn
                        scr2(1)=scr2(1) 
     &                  +2*wt(ipt)*xc_dhpt(ipt,iragbb)*gbgmn*
     &                     drhodb(ipt,1,lpert+1)
     &                  +  wt(ipt)*xc_dhpt(ipt,iragab)*gagmn*
     &                     drhodb(ipt,1,lpert+1)
     &                  +2*wt(ipt)*xc_dhpt(ipt,irbgbb)*gbgmn*
     &                     drhodb(ipt,2,lpert+1)
     &                  +  wt(ipt)*xc_dhpt(ipt,irbgab)*gagmn*
     &                     drhodb(ipt,2,lpert+1)
     &                  +  wt(ipt)*xc_dhpt(ipt,irbgaa)*gaadbx*
     &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
     &                  +  wt(ipt)*xc_dhpt(ipt,irbgab)*gabdbx*
     &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
     &                  +  wt(ipt)*xc_dhpt(ipt,irbgbb)*gbbdbx*
     &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
     &                  +2*wt(ipt)*xc_dhpt(ipt,igaagbb)*gbgmn*gaadbx
     &                  +  wt(ipt)*xc_dhpt(ipt,igaagab)*gagmn*gaadbx
     &                  +2*wt(ipt)*xc_dhpt(ipt,igabgbb)*gbgmn*gabdbx
     &                  +  wt(ipt)*xc_dhpt(ipt,igabgab)*gagmn*gabdbx
     &                  +2*wt(ipt)*xc_dhpt(ipt,igbbgbb)*gbgmn*gbbdbx
     &                  +  wt(ipt)*xc_dhpt(ipt,igabgbb)*gagmn*gbbdbx
     &                  +2*wt(ipt)*xc_dvpt(ipt,igbb)*gbdbxgmn
     &                  +  wt(ipt)*xc_dvpt(ipt,igab)*gadbxgmn
                        scr2(2)=scr2(2) 
     &                  +2*wt(ipt)*xc_dhpt(ipt,iragbb)*gbgmn*
     &                     drhodb(ipt,1,lpert+2)
     &                  +  wt(ipt)*xc_dhpt(ipt,iragab)*gagmn*
     &                     drhodb(ipt,1,lpert+2)
     &                  +2*wt(ipt)*xc_dhpt(ipt,irbgbb)*gbgmn*
     &                     drhodb(ipt,2,lpert+2)
     &                  +  wt(ipt)*xc_dhpt(ipt,irbgab)*gagmn*
     &                     drhodb(ipt,2,lpert+2)
     &                  +  wt(ipt)*xc_dhpt(ipt,irbgaa)*gaadby*
     &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
     &                  +  wt(ipt)*xc_dhpt(ipt,irbgab)*gabdby*
     &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
     &                  +  wt(ipt)*xc_dhpt(ipt,irbgbb)*gbbdby*
     &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
     &                  +2*wt(ipt)*xc_dhpt(ipt,igaagbb)*gbgmn*gaadby
     &                  +  wt(ipt)*xc_dhpt(ipt,igaagab)*gagmn*gaadby
     &                  +2*wt(ipt)*xc_dhpt(ipt,igabgbb)*gbgmn*gabdby
     &                  +  wt(ipt)*xc_dhpt(ipt,igabgab)*gagmn*gabdby
     &                  +2*wt(ipt)*xc_dhpt(ipt,igbbgbb)*gbgmn*gbbdby
     &                  +  wt(ipt)*xc_dhpt(ipt,igabgbb)*gagmn*gbbdby
     &                  +2*wt(ipt)*xc_dvpt(ipt,igbb)*gbdbygmn
     &                  +  wt(ipt)*xc_dvpt(ipt,igab)*gadbygmn
                        scr2(3)=scr2(3) 
     &                  +2*wt(ipt)*xc_dhpt(ipt,iragbb)*gbgmn*
     &                     drhodb(ipt,1,lpert+3)
     &                  +  wt(ipt)*xc_dhpt(ipt,iragab)*gagmn*
     &                     drhodb(ipt,1,lpert+3)
     &                  +2*wt(ipt)*xc_dhpt(ipt,irbgbb)*gbgmn*
     &                     drhodb(ipt,2,lpert+3)
     &                  +  wt(ipt)*xc_dhpt(ipt,irbgab)*gagmn*
     &                     drhodb(ipt,2,ipert+3)
     &                  +  wt(ipt)*xc_dhpt(ipt,irbgaa)*gaadbz*
     &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
     &                  +  wt(ipt)*xc_dhpt(ipt,irbgab)*gabdbz*
     &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
     &                  +  wt(ipt)*xc_dhpt(ipt,irbgbb)*gbbdbz*
     &                     bfn_val(ipt,mu)*bfn_val(ipt,nu)
     &                  +2*wt(ipt)*xc_dhpt(ipt,igaagbb)*gbgmn*gaadbz
     &                  +  wt(ipt)*xc_dhpt(ipt,igaagab)*gagmn*gaadbz
     &                  +2*wt(ipt)*xc_dhpt(ipt,igabgbb)*gbgmn*gabdbz
     &                  +  wt(ipt)*xc_dhpt(ipt,igabgab)*gagmn*gabdbz
     &                  +2*wt(ipt)*xc_dhpt(ipt,igbbgbb)*gbgmn*gbbdbz
     &                  +  wt(ipt)*xc_dhpt(ipt,igabgbb)*gagmn*gbbdbz
     &                  +2*wt(ipt)*xc_dvpt(ipt,igbb)*gbdbzgmn
     &                  +  wt(ipt)*xc_dvpt(ipt,igab)*gadbzgmn
                     enddo
                     qa(n,ipert+1)=qa(n,ipert+1)+scr1(1)
                     qa(n,ipert+2)=qa(n,ipert+2)+scr1(2)
                     qa(n,ipert+3)=qa(n,ipert+3)+scr1(3)
                     qb(n,ipert+1)=qb(n,ipert+1)+scr2(1)
                     qb(n,ipert+2)=qb(n,ipert+2)+scr2(2)
                     qb(n,ipert+3)=qb(n,ipert+3)+scr2(3)
                     if (gwt_sw) then
                        qa(n,jpert+1)=qa(n,jpert+1)-scr1(1)
                        qa(n,jpert+2)=qa(n,jpert+2)-scr1(2)
                        qa(n,jpert+3)=qa(n,jpert+3)-scr1(3)
                        qb(n,jpert+1)=qb(n,jpert+1)-scr2(1)
                        qb(n,jpert+2)=qb(n,jpert+2)-scr2(2)
                        qb(n,jpert+3)=qb(n,jpert+3)-scr2(3)
                     endif
                  endif
               enddo
            enddo
 200     continue
c
         do 210 mu=1,n_active_bfn
            ipert=3*(active_bfn_atms(active_bfn_indx(mu))-1)
            if (ipert.eq.jpert.and.gwt_sw) goto 210
            iki=iky(active_bfn_list(mu))
            do nu=1, mu
               n=iki+active_bfn_list(nu)
c
               call aclear_dp(scr1,3,0.0d0)
               call aclear_dp(scr2,3,0.0d0)
               do ipt=1,npts
                  scr1(1)=scr1(1)
     &            -wt(ipt)*xc_vpt(ipt,ira)*bfng_val(ipt,mu,1)*
     &             bfn_val(ipt,nu)
                  scr1(2)=scr1(2)
     &            -wt(ipt)*xc_vpt(ipt,ira)*bfng_val(ipt,mu,2)*
     &             bfn_val(ipt,nu)
                  scr1(3)=scr1(3)
     &            -wt(ipt)*xc_vpt(ipt,ira)*bfng_val(ipt,mu,3)*
     &             bfn_val(ipt,nu)
                  scr2(1)=scr2(1)
     &            -wt(ipt)*xc_vpt(ipt,irb)*bfng_val(ipt,mu,1)*
     &             bfn_val(ipt,nu)
                  scr2(2)=scr2(2)
     &            -wt(ipt)*xc_vpt(ipt,irb)*bfng_val(ipt,mu,2)*
     &             bfn_val(ipt,nu)
                  scr2(3)=scr2(3)
     &            -wt(ipt)*xc_vpt(ipt,irb)*bfng_val(ipt,mu,3)*
     &             bfn_val(ipt,nu)
               enddo
               if (gradcorr_sw) then
                  do ipt=1,npts
                     gagmndbx = grho(ipt,1,1)*(
     &                        - bfng_val(ipt,mu,1)*bfng_val(ipt,nu,1)
     &                        - bfn_hess(ipt,mu,hxx)*bfn_val(ipt,nu))
     &                        + grho(ipt,1,2)*(
     &                        - bfng_val(ipt,mu,1)*bfng_val(ipt,nu,2)
     &                        - bfn_hess(ipt,mu,hxy)*bfn_val(ipt,nu))
     &                        + grho(ipt,1,3)*(
     &                        - bfng_val(ipt,mu,1)*bfng_val(ipt,nu,3)
     &                        - bfn_hess(ipt,mu,hxz)*bfn_val(ipt,nu))
                     gbgmndbx = grho(ipt,2,1)*(
     &                        - bfng_val(ipt,mu,1)*bfng_val(ipt,nu,1)
     &                        - bfn_hess(ipt,mu,hxx)*bfn_val(ipt,nu))
     &                        + grho(ipt,2,2)*(
     &                        - bfng_val(ipt,mu,1)*bfng_val(ipt,nu,2)
     &                        - bfn_hess(ipt,mu,hxy)*bfn_val(ipt,nu))
     &                        + grho(ipt,2,3)*(
     &                        - bfng_val(ipt,mu,1)*bfng_val(ipt,nu,3)
     &                        - bfn_hess(ipt,mu,hxz)*bfn_val(ipt,nu))
                     gagmndby = grho(ipt,1,1)*(
     &                        - bfng_val(ipt,mu,2)*bfng_val(ipt,nu,1)
     &                        - bfn_hess(ipt,mu,hxy)*bfn_val(ipt,nu))
     &                        + grho(ipt,1,2)*(
     &                        - bfng_val(ipt,mu,2)*bfng_val(ipt,nu,2)
     &                        - bfn_hess(ipt,mu,hyy)*bfn_val(ipt,nu))
     &                        + grho(ipt,1,3)*(
     &                        - bfng_val(ipt,mu,2)*bfng_val(ipt,nu,3)
     &                        - bfn_hess(ipt,mu,hyz)*bfn_val(ipt,nu))
                     gbgmndby = grho(ipt,2,1)*(
     &                        - bfng_val(ipt,mu,2)*bfng_val(ipt,nu,1)
     &                        - bfn_hess(ipt,mu,hxy)*bfn_val(ipt,nu))
     &                        + grho(ipt,2,2)*(
     &                        - bfng_val(ipt,mu,2)*bfng_val(ipt,nu,2)
     &                        - bfn_hess(ipt,mu,hyy)*bfn_val(ipt,nu))
     &                        + grho(ipt,2,3)*(
     &                        - bfng_val(ipt,mu,2)*bfng_val(ipt,nu,3)
     &                        - bfn_hess(ipt,mu,hyz)*bfn_val(ipt,nu))
                     gagmndbz = grho(ipt,1,1)*(
     &                        - bfng_val(ipt,mu,3)*bfng_val(ipt,nu,1)
     &                        - bfn_hess(ipt,mu,hxz)*bfn_val(ipt,nu))
     &                        + grho(ipt,1,2)*(
     &                        - bfng_val(ipt,mu,3)*bfng_val(ipt,nu,2)
     &                        - bfn_hess(ipt,mu,hyz)*bfn_val(ipt,nu))
     &                        + grho(ipt,1,3)*(
     &                        - bfng_val(ipt,mu,3)*bfng_val(ipt,nu,3)
     &                        - bfn_hess(ipt,mu,hzz)*bfn_val(ipt,nu))
                     gbgmndbz = grho(ipt,2,1)*(
     &                        - bfng_val(ipt,mu,3)*bfng_val(ipt,nu,1)
     &                        - bfn_hess(ipt,mu,hxz)*bfn_val(ipt,nu))
     &                        + grho(ipt,2,2)*(
     &                        - bfng_val(ipt,mu,3)*bfng_val(ipt,nu,2)
     &                        - bfn_hess(ipt,mu,hyz)*bfn_val(ipt,nu))
     &                        + grho(ipt,2,3)*(
     &                        - bfng_val(ipt,mu,3)*bfng_val(ipt,nu,3)
     &                        - bfn_hess(ipt,mu,hzz)*bfn_val(ipt,nu))
                     scr1(1) = scr1(1)
     &                       + 2*wt(ipt)*xc_dvpt(ipt,igaa)*gagmndbx
     &                       +   wt(ipt)*xc_dvpt(ipt,igab)*gbgmndbx
                     scr1(2) = scr1(2)
     &                       + 2*wt(ipt)*xc_dvpt(ipt,igaa)*gagmndby
     &                       +   wt(ipt)*xc_dvpt(ipt,igab)*gbgmndby
                     scr1(3) = scr1(3)
     &                       + 2*wt(ipt)*xc_dvpt(ipt,igaa)*gagmndbz
     &                       +   wt(ipt)*xc_dvpt(ipt,igab)*gbgmndbz
                     scr2(1) = scr2(1)
     &                       + 2*wt(ipt)*xc_dvpt(ipt,igbb)*gbgmndbx
     &                       +   wt(ipt)*xc_dvpt(ipt,igab)*gagmndbx
                     scr2(2) = scr2(2)
     &                       + 2*wt(ipt)*xc_dvpt(ipt,igbb)*gbgmndby
     &                       +   wt(ipt)*xc_dvpt(ipt,igab)*gagmndby
                     scr2(3) = scr2(3)
     &                       + 2*wt(ipt)*xc_dvpt(ipt,igbb)*gbgmndbz
     &                       +   wt(ipt)*xc_dvpt(ipt,igab)*gagmndbz
                  enddo
               endif
               qa(n,ipert+1)=qa(n,ipert+1)+scr1(1)
               qa(n,ipert+2)=qa(n,ipert+2)+scr1(2)
               qa(n,ipert+3)=qa(n,ipert+3)+scr1(3)
               qb(n,ipert+1)=qb(n,ipert+1)+scr2(1)
               qb(n,ipert+2)=qb(n,ipert+2)+scr2(2)
               qb(n,ipert+3)=qb(n,ipert+3)+scr2(3)
               if (gwt_sw) then
                  qa(n,jpert+1)=qa(n,jpert+1)-scr1(1)
                  qa(n,jpert+2)=qa(n,jpert+2)-scr1(2)
                  qa(n,jpert+3)=qa(n,jpert+3)-scr1(3)
                  qb(n,jpert+1)=qb(n,jpert+1)-scr2(1)
                  qb(n,jpert+2)=qb(n,jpert+2)-scr2(2)
                  qb(n,jpert+3)=qb(n,jpert+3)-scr2(3)
               endif
            enddo ! nu
            do nu=mu, n_active_bfn
               iki = iky(active_bfn_list(nu))
               n   = iki + active_bfn_list(mu)
c
               call aclear_dp(scr1,3,0.0d0)
               call aclear_dp(scr2,3,0.0d0)
               do ipt=1,npts
                  scr1(1)=scr1(1)
     &            -wt(ipt)*xc_vpt(ipt,ira)*bfng_val(ipt,mu,1)*
     &             bfn_val(ipt,nu)
                  scr1(2)=scr1(2)
     &            -wt(ipt)*xc_vpt(ipt,ira)*bfng_val(ipt,mu,2)*
     &             bfn_val(ipt,nu)
                  scr1(3)=scr1(3)
     &            -wt(ipt)*xc_vpt(ipt,ira)*bfng_val(ipt,mu,3)*
     &             bfn_val(ipt,nu)
                  scr2(1)=scr2(1)
     &            -wt(ipt)*xc_vpt(ipt,irb)*bfng_val(ipt,mu,1)*
     &             bfn_val(ipt,nu)
                  scr2(2)=scr2(2)
     &            -wt(ipt)*xc_vpt(ipt,irb)*bfng_val(ipt,mu,2)*
     &             bfn_val(ipt,nu)
                  scr2(3)=scr2(3)
     &            -wt(ipt)*xc_vpt(ipt,irb)*bfng_val(ipt,mu,3)*
     &             bfn_val(ipt,nu)
               enddo
               if (gradcorr_sw) then
                  do ipt=1,npts
                     gagmndbx = grho(ipt,1,1)*(
     &                        - bfng_val(ipt,mu,1)*bfng_val(ipt,nu,1)
     &                        - bfn_hess(ipt,mu,hxx)*bfn_val(ipt,nu))
     &                        + grho(ipt,1,2)*(
     &                        - bfng_val(ipt,mu,1)*bfng_val(ipt,nu,2)
     &                        - bfn_hess(ipt,mu,hxy)*bfn_val(ipt,nu))
     &                        + grho(ipt,1,3)*(
     &                        - bfng_val(ipt,mu,1)*bfng_val(ipt,nu,3)
     &                        - bfn_hess(ipt,mu,hxz)*bfn_val(ipt,nu))
                     gbgmndbx = grho(ipt,2,1)*(
     &                        - bfng_val(ipt,mu,1)*bfng_val(ipt,nu,1)
     &                        - bfn_hess(ipt,mu,hxx)*bfn_val(ipt,nu))
     &                        + grho(ipt,2,2)*(
     &                        - bfng_val(ipt,mu,1)*bfng_val(ipt,nu,2)
     &                        - bfn_hess(ipt,mu,hxy)*bfn_val(ipt,nu))
     &                        + grho(ipt,2,3)*(
     &                        - bfng_val(ipt,mu,1)*bfng_val(ipt,nu,3)
     &                        - bfn_hess(ipt,mu,hxz)*bfn_val(ipt,nu))
                     gagmndby = grho(ipt,1,1)*(
     &                        - bfng_val(ipt,mu,2)*bfng_val(ipt,nu,1)
     &                        - bfn_hess(ipt,mu,hxy)*bfn_val(ipt,nu))
     &                        + grho(ipt,1,2)*(
     &                        - bfng_val(ipt,mu,2)*bfng_val(ipt,nu,2)
     &                        - bfn_hess(ipt,mu,hyy)*bfn_val(ipt,nu))
     &                        + grho(ipt,1,3)*(
     &                        - bfng_val(ipt,mu,2)*bfng_val(ipt,nu,3)
     &                        - bfn_hess(ipt,mu,hyz)*bfn_val(ipt,nu))
                     gbgmndby = grho(ipt,2,1)*(
     &                        - bfng_val(ipt,mu,2)*bfng_val(ipt,nu,1)
     &                        - bfn_hess(ipt,mu,hxy)*bfn_val(ipt,nu))
     &                        + grho(ipt,2,2)*(
     &                        - bfng_val(ipt,mu,2)*bfng_val(ipt,nu,2)
     &                        - bfn_hess(ipt,mu,hyy)*bfn_val(ipt,nu))
     &                        + grho(ipt,2,3)*(
     &                        - bfng_val(ipt,mu,2)*bfng_val(ipt,nu,3)
     &                        - bfn_hess(ipt,mu,hyz)*bfn_val(ipt,nu))
                     gagmndbz = grho(ipt,1,1)*(
     &                        - bfng_val(ipt,mu,3)*bfng_val(ipt,nu,1)
     &                        - bfn_hess(ipt,mu,hxz)*bfn_val(ipt,nu))
     &                        + grho(ipt,1,2)*(
     &                        - bfng_val(ipt,mu,3)*bfng_val(ipt,nu,2)
     &                        - bfn_hess(ipt,mu,hyz)*bfn_val(ipt,nu))
     &                        + grho(ipt,1,3)*(
     &                        - bfng_val(ipt,mu,3)*bfng_val(ipt,nu,3)
     &                        - bfn_hess(ipt,mu,hzz)*bfn_val(ipt,nu))
                     gbgmndbz = grho(ipt,2,1)*(
     &                        - bfng_val(ipt,mu,3)*bfng_val(ipt,nu,1)
     &                        - bfn_hess(ipt,mu,hxz)*bfn_val(ipt,nu))
     &                        + grho(ipt,2,2)*(
     &                        - bfng_val(ipt,mu,3)*bfng_val(ipt,nu,2)
     &                        - bfn_hess(ipt,mu,hyz)*bfn_val(ipt,nu))
     &                        + grho(ipt,2,3)*(
     &                        - bfng_val(ipt,mu,3)*bfng_val(ipt,nu,3)
     &                        - bfn_hess(ipt,mu,hzz)*bfn_val(ipt,nu))
                     scr1(1) = scr1(1)
     &                       + 2*wt(ipt)*xc_dvpt(ipt,igaa)*gagmndbx
     &                       +   wt(ipt)*xc_dvpt(ipt,igab)*gbgmndbx
                     scr1(2) = scr1(2)
     &                       + 2*wt(ipt)*xc_dvpt(ipt,igaa)*gagmndby
     &                       +   wt(ipt)*xc_dvpt(ipt,igab)*gbgmndby
                     scr1(3) = scr1(3)
     &                       + 2*wt(ipt)*xc_dvpt(ipt,igaa)*gagmndbz
     &                       +   wt(ipt)*xc_dvpt(ipt,igab)*gbgmndbz
                     scr2(1) = scr2(1)
     &                       + 2*wt(ipt)*xc_dvpt(ipt,igbb)*gbgmndbx
     &                       +   wt(ipt)*xc_dvpt(ipt,igab)*gagmndbx
                     scr2(2) = scr2(2)
     &                       + 2*wt(ipt)*xc_dvpt(ipt,igbb)*gbgmndby
     &                       +   wt(ipt)*xc_dvpt(ipt,igab)*gagmndby
                     scr2(3) = scr2(3)
     &                       + 2*wt(ipt)*xc_dvpt(ipt,igbb)*gbgmndbz
     &                       +   wt(ipt)*xc_dvpt(ipt,igab)*gagmndbz
                  enddo
               endif
               qa(n,ipert+1)=qa(n,ipert+1)+scr1(1)
               qa(n,ipert+2)=qa(n,ipert+2)+scr1(2)
               qa(n,ipert+3)=qa(n,ipert+3)+scr1(3)
               qb(n,ipert+1)=qb(n,ipert+1)+scr2(1)
               qb(n,ipert+2)=qb(n,ipert+2)+scr2(2)
               qb(n,ipert+3)=qb(n,ipert+3)+scr2(3)
               if (gwt_sw) then
                  qa(n,jpert+1)=qa(n,jpert+1)-scr1(1)
                  qa(n,jpert+2)=qa(n,jpert+2)-scr1(2)
                  qa(n,jpert+3)=qa(n,jpert+3)-scr1(3)
                  qb(n,jpert+1)=qb(n,jpert+1)-scr2(1)
                  qb(n,jpert+2)=qb(n,jpert+2)-scr2(2)
                  qb(n,jpert+3)=qb(n,jpert+3)-scr2(3)
               endif
            enddo
 210     continue
      endif ! rkstyp_sw
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine rks_dksm_exp_mo(npts,nao,nvec,naocc,npert,ngridcent,
     &           mxp,
     &           wt,xc_vpt,xc_hpt,amo_val,bfng_val,
     &           avec,dmo,drho,tmo,qa)
      implicit none
c
c     Compute the explicit derivative contribution to the Kohn-Sham
c     matrix in MO-basis for the local density closed shell case.
c
c     The tactics are to construct to first and second derivatives
c     of the density with respect to the nuclear coordinates on the
c     fly first and subsequently consume them. Derivates of the MO's
c     are generated per orbital so as to minimise memory requirements.
c
c     Inputs
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
C *Parameter module for CCP1/DFT            
C *Memory information		
C *----------------
C * max_block 	- maximum number of blocks which can be allocate
      integer max_block
      parameter(max_block=20)
C *							
C *Basis sets information                              
C *----------------------				
C *max_tag 	- maximum number of basis sets which can be alloc
C *max_atype	- maximum number of centre types	
C *max_gtype	- maximum number of grid types (at least every element 
C                 has different grid type)
C *max_grids    - maximum number of grids (different terms may have
C                 different grids, i.e. CPKS equations use a different
C                 grid than the one used for the KS-matrix)
C *max_shel	- maximum number of shells on centre
C *max_prm	- maximum number of primitives for any given centre
C *maxL		- maximum angular momentum allowed	
C *max_func	- maximum number of basis functions for any centre
      integer max_tag,max_atype,max_shel,max_prm,max_ang
      integer max_gtype,max_grids
      parameter(max_tag=3,max_atype=30,max_shel=500,max_prm=5000)
      parameter(max_ang=5,max_gtype=10,max_grids=2)
C *								
C *Geometry information					
C *--------------------					
C *max_atom	- maximum number of atoms in system
      integer max_atom
      parameter(max_atom=750)
C *
C *Accuracy information				
C ---------------------			
C *global_accuracy 	- global accuracy	
      real*8  global_accuracy
      parameter(global_accuracy=1.0d-14)

C *Grid information
c *----------------
c *maxradzn - The maximum number of radial zones. 
      integer maxgpt,maxrad,maxfpt,maxang,maxradzn,maxtablerows
      parameter(maxgpt=2900,maxrad=50,maxang=302,maxfpt=100)
      parameter(maxradzn=35,maxtablerows=7)
c
c     The parameters defined here label the 1st, 2nd and 3rd derivative
c     components of the basis functions.
c
c     2nd derivatives:
c
      integer hxx, hyy, hzz, hxy, hxz, hyz
c
      parameter (hxx = 1)
      parameter (hyy = 2)
      parameter (hzz = 3)
      parameter (hxy = 4)
      parameter (hxz = 5)
      parameter (hyz = 6)
c
      integer hcc(3,3)
      data hcc/hxx,hxy,hxz,
     &         hxy,hyy,hyz,
     &         hxz,hyz,hzz/
c
c     3rd derivatives:
c
      integer txxx, txxy, txxz, txyy, txyz
      integer txzz, tyyy, tyyz, tyzz, tzzz
c
      parameter (txxx =  1)
      parameter (txxy =  2)
      parameter (txxz =  3)
      parameter (txyy =  4)
      parameter (txyz =  5)
      parameter (txzz =  6)
      parameter (tyyy =  7)
      parameter (tyyz =  8)
      parameter (tyzz =  9)
      parameter (tzzz = 10)
c
      integer tc3(3,3,3)
      data tc3/txxx,txxy,txxz,
     &         txxy,txyy,txyz,
     &         txxz,txyz,txzz,
     &         txxy,txyy,txyz,
     &         txyy,tyyy,tyyz,
     &         txyz,tyyz,tyzz,
     &         txxz,txyz,txzz,
     &         txyz,tyyz,tyzz,
     &         txzz,tyzz,tzzz/
c
c     The parameters defined here label the 1st and 2nd derivative
c     components of the density functional. Because a functional
c     may be purely local or of a GGA type the components are grouped
c     into derivative with respect to the density only and others.
c
c     1st derivatives:
c
      integer ira, irb
      integer igaa, igab, igbb
      integer ita, itb
c
c     2nd derivatives:
c
      integer irara, irarb, irbrb
      integer iragaa, iragab, iragbb, irbgaa, irbgab, irbgbb
      integer igaagaa, igaagab, igaagbb, igabgab, igabgbb, igbbgbb
      integer irata, iratb, irbta, irbtb
      integer igaata, igaatb, igabta, igabtb, igbbta, igbbtb
      integer itata, itatb, itbtb
c
      parameter (ira     = 1)
      parameter (irb     = 2)
c
      parameter (igaa    = 1)
      parameter (igab    = 2)
      parameter (igbb    = 3)
c
      parameter (ita     = 1)
      parameter (itb     = 2)
c
c
      parameter (irara   = 1)
      parameter (irarb   = 2)
      parameter (irbrb   = 3)
c
      parameter (iragaa  = 1)
      parameter (iragab  = 3)
      parameter (iragbb  = 4)
      parameter (irbgaa  = 5)
      parameter (irbgab  = 6)
      parameter (irbgbb  = 7)
      parameter (igaagaa = 2)
      parameter (igaagab = 8)
      parameter (igaagbb = 9)
      parameter (igabgab = 10)
      parameter (igabgbb = 11)
      parameter (igbbgbb = 12)
c
      parameter (irata   = 1)
      parameter (iratb   = 3)
      parameter (irbta   = 4)
      parameter (irbtb   = 5)
      parameter (itata   = 2)
      parameter (itatb   = 6)
      parameter (itbtb   = 7)
c
      parameter (igaata  = 1)
      parameter (igaatb  = 2)
      parameter (igabta  = 3)
      parameter (igabtb  = 4)
      parameter (igbbta  = 5)
      parameter (igbbtb  = 6)
c
      integer npts
      integer nao
      integer nvec
      integer naocc 
      integer npert
      integer ngridcent
      integer mxp 
      real*8 wt(npts)
      real*8 xc_vpt(mxp,2)
      real*8 xc_hpt(mxp,3)
      real*8 amo_val(npts,nvec)
      real*8 bfng_val(mxp,nao,3)
      real*8 avec(nao,nvec)
c
c     Outputs
c
      real*8 qa(nvec*(nvec+1)/2,npert)
c
c     Workspace
c
      real*8 dmo(npts,npert)
c     real*8 ddmo(npts,6,ngridcent)
      real*8 drho(npts,npert)
      real*8 tmo(npts,nvec)
c     real*8 ddrho(npts,npert*(npert+1)/2)
c
c     Local variables
c
      integer i, j
      integer n, n1
      integer ipert
      integer ipt
c     integer iatm
c     integer ic, jc
      integer iki, ikj
c     real*8 tmp
      integer ao_tag
      parameter(ao_tag=1)
c
c     Code
c
      call aclear_dp(drho,npts*npert,0.0d0)
      do ipt=1,npts
         dmo(ipt,1)=wt(ipt)*xc_vpt(ipt,ira)
      enddo
      do i=1,nvec
         do ipt=1,npts
            tmo(ipt,i)=dmo(ipt,1)*amo_val(ipt,i)
         enddo
      enddo
c
c     Construct the derivatives of the density (times 0.5)
c
      iki=0
      do i=1,naocc
         call dmo_val(npts,nao,nvec,npert,ngridcent,mxp,i,
     &                bfng_val,avec,dmo)
         do ipert=1,npert
            do ipt=1,npts
               drho(ipt,ipert)=drho(ipt,ipert)
     &                        +amo_val(ipt,i)*dmo(ipt,ipert)
            enddo
         enddo
c
         iki=iki+i-1
         do ipert=1,npert
            do j=1,i-1
               do ipt=1,npts
                  qa(iki+j,ipert)=qa(iki+j,ipert)
     &            +dmo(ipt,ipert)*tmo(ipt,j)
               enddo ! ipt
            enddo ! j
            do ipt=1,npts
               qa(iki+i,ipert)=qa(iki+i,ipert)
     &         +2*dmo(ipt,ipert)*tmo(ipt,i)
            enddo ! ipt
            ikj=iki
            do j=i+1,nvec
               ikj=ikj+j-1
               do ipt=1,npts
                  qa(ikj+i,ipert)=qa(ikj+i,ipert)
     &            +dmo(ipt,ipert)*tmo(ipt,j)
               enddo
            enddo ! j
         enddo ! ipert
      enddo
c
      n=(naocc-1)*naocc/2
      do i=naocc+1,nvec
         call dmo_val(npts,nao,nvec,npert,ngridcent,mxp,i,
     &                bfng_val,avec,dmo)
         n=n+i-1
         do ipert=1,npert
            do j=1,naocc
               do ipt=1,npts
                  qa(n+j,ipert)=qa(n+j,ipert)
     &            +dmo(ipt,ipert)*tmo(ipt,j)
               enddo ! ipt
            enddo ! j
         enddo ! ipert
      enddo
c
c     Add the contributions to the Kohn-Sham matrix
c
      do ipt=1,npts 
         dmo(ipt,1)=2*wt(ipt)*xc_hpt(ipt,irara)
      enddo
      do ipert=1,npert
         do ipt=1,npts 
            dmo(ipt,2)=dmo(ipt,1)*drho(ipt,ipert)
         enddo
         n=0
         do j=1,naocc
            n1=n
            n=n+j-1
            do ipt=1,npts
               dmo(ipt,3)=dmo(ipt,2)*amo_val(ipt,j)
            enddo
            do i=j,nvec
               n1=n1+i-1
               do ipt=1,npts
                  qa(n1+j,ipert)=qa(n1+j,ipert)
     &            +dmo(ipt,3)*amo_val(ipt,i)
               enddo ! ipt
            enddo ! j
         enddo ! i
      enddo ! ipert
c 
      end
c
c-----------------------------------------------------------------------
c
      subroutine rks_dksm_exp_mo_gga(npts,nao,nvec,naocc,npert,
     &           ngridcent,mxp,wt,xc_vpt,xc_dvpt,xc_hpt,xc_dhpt,
     &           amo_val,amo_grad,bfng_val,bfn_hess,avec,grho,dmo,dgmo,
     &           drho,dgrho,gmorho,grdgmo,t1,t2,t3,t4,tptmo,qa)
      implicit none
c
c     Compute the explicit derivative contribution to the Kohn-Sham 
c     matrix in MO-basis for the gradient corrected closed shell case.
c
c     The tactics are to construct to first and second derivatives
c     of the density with respect to the nuclear coordinates on the
c     fly first and subsequently consume them. Derivates of the MO's
c     are generated per orbital so as to minimise memory requirements.
c
c     Inputs
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
C *Parameter module for CCP1/DFT            
C *Memory information		
C *----------------
C * max_block 	- maximum number of blocks which can be allocate
      integer max_block
      parameter(max_block=20)
C *							
C *Basis sets information                              
C *----------------------				
C *max_tag 	- maximum number of basis sets which can be alloc
C *max_atype	- maximum number of centre types	
C *max_gtype	- maximum number of grid types (at least every element 
C                 has different grid type)
C *max_grids    - maximum number of grids (different terms may have
C                 different grids, i.e. CPKS equations use a different
C                 grid than the one used for the KS-matrix)
C *max_shel	- maximum number of shells on centre
C *max_prm	- maximum number of primitives for any given centre
C *maxL		- maximum angular momentum allowed	
C *max_func	- maximum number of basis functions for any centre
      integer max_tag,max_atype,max_shel,max_prm,max_ang
      integer max_gtype,max_grids
      parameter(max_tag=3,max_atype=30,max_shel=500,max_prm=5000)
      parameter(max_ang=5,max_gtype=10,max_grids=2)
C *								
C *Geometry information					
C *--------------------					
C *max_atom	- maximum number of atoms in system
      integer max_atom
      parameter(max_atom=750)
C *
C *Accuracy information				
C ---------------------			
C *global_accuracy 	- global accuracy	
      real*8  global_accuracy
      parameter(global_accuracy=1.0d-14)

C *Grid information
c *----------------
c *maxradzn - The maximum number of radial zones. 
      integer maxgpt,maxrad,maxfpt,maxang,maxradzn,maxtablerows
      parameter(maxgpt=2900,maxrad=50,maxang=302,maxfpt=100)
      parameter(maxradzn=35,maxtablerows=7)
c
c     The parameters defined here label the 1st, 2nd and 3rd derivative
c     components of the basis functions.
c
c     2nd derivatives:
c
      integer hxx, hyy, hzz, hxy, hxz, hyz
c
      parameter (hxx = 1)
      parameter (hyy = 2)
      parameter (hzz = 3)
      parameter (hxy = 4)
      parameter (hxz = 5)
      parameter (hyz = 6)
c
      integer hcc(3,3)
      data hcc/hxx,hxy,hxz,
     &         hxy,hyy,hyz,
     &         hxz,hyz,hzz/
c
c     3rd derivatives:
c
      integer txxx, txxy, txxz, txyy, txyz
      integer txzz, tyyy, tyyz, tyzz, tzzz
c
      parameter (txxx =  1)
      parameter (txxy =  2)
      parameter (txxz =  3)
      parameter (txyy =  4)
      parameter (txyz =  5)
      parameter (txzz =  6)
      parameter (tyyy =  7)
      parameter (tyyz =  8)
      parameter (tyzz =  9)
      parameter (tzzz = 10)
c
      integer tc3(3,3,3)
      data tc3/txxx,txxy,txxz,
     &         txxy,txyy,txyz,
     &         txxz,txyz,txzz,
     &         txxy,txyy,txyz,
     &         txyy,tyyy,tyyz,
     &         txyz,tyyz,tyzz,
     &         txxz,txyz,txzz,
     &         txyz,tyyz,tyzz,
     &         txzz,tyzz,tzzz/
c
c     The parameters defined here label the 1st and 2nd derivative
c     components of the density functional. Because a functional
c     may be purely local or of a GGA type the components are grouped
c     into derivative with respect to the density only and others.
c
c     1st derivatives:
c
      integer ira, irb
      integer igaa, igab, igbb
      integer ita, itb
c
c     2nd derivatives:
c
      integer irara, irarb, irbrb
      integer iragaa, iragab, iragbb, irbgaa, irbgab, irbgbb
      integer igaagaa, igaagab, igaagbb, igabgab, igabgbb, igbbgbb
      integer irata, iratb, irbta, irbtb
      integer igaata, igaatb, igabta, igabtb, igbbta, igbbtb
      integer itata, itatb, itbtb
c
      parameter (ira     = 1)
      parameter (irb     = 2)
c
      parameter (igaa    = 1)
      parameter (igab    = 2)
      parameter (igbb    = 3)
c
      parameter (ita     = 1)
      parameter (itb     = 2)
c
c
      parameter (irara   = 1)
      parameter (irarb   = 2)
      parameter (irbrb   = 3)
c
      parameter (iragaa  = 1)
      parameter (iragab  = 3)
      parameter (iragbb  = 4)
      parameter (irbgaa  = 5)
      parameter (irbgab  = 6)
      parameter (irbgbb  = 7)
      parameter (igaagaa = 2)
      parameter (igaagab = 8)
      parameter (igaagbb = 9)
      parameter (igabgab = 10)
      parameter (igabgbb = 11)
      parameter (igbbgbb = 12)
c
      parameter (irata   = 1)
      parameter (iratb   = 3)
      parameter (irbta   = 4)
      parameter (irbtb   = 5)
      parameter (itata   = 2)
      parameter (itatb   = 6)
      parameter (itbtb   = 7)
c
      parameter (igaata  = 1)
      parameter (igaatb  = 2)
      parameter (igabta  = 3)
      parameter (igabtb  = 4)
      parameter (igbbta  = 5)
      parameter (igbbtb  = 6)
c
      integer npts
      integer nao
      integer nvec
      integer naocc 
      integer npert
      integer ngridcent
      integer mxp 
      real*8 wt(npts)
      real*8 xc_vpt(mxp,2)
      real*8 xc_dvpt(mxp,3)
      real*8 xc_hpt(mxp,3)
      real*8 xc_dhpt(mxp,12)
      real*8 amo_val(npts,nvec)
      real*8 amo_grad(npts,nvec,3)
      real*8 bfng_val(mxp,nao,3)
      real*8 bfn_hess(mxp,nao,6)
c     real*8 bfn_3rd(mxp,nao,10)
      real*8 avec(nao,nvec)
      real*8 grho(mxp,2,3)
c
c     Outputs
c
      real*8 qa(nao*(nao+1)/2,npert)
c
c     Workspace
c
      real*8 dmo(npts,npert)
      real*8 dgmo(npts,npert,3)
c     real*8 ddmo(npts,6,ngridcent)
c     real*8 ddgmo(npts,10,ngridcent)
      real*8 drho(npts,npert)
c     real*8 ddrho(npts,npert*(npert+1)/2)
      real*8 dgrho(npts,npert,3)
c     real*8 ddgrho(npts,npert*(npert+1)/2,3)
      real*8 gmorho(npts,nvec)
      real*8 grdgmo(npts)
      real*8 t1(npts)
      real*8 t2(npts)
      real*8 t3(npts)
      real*8 t4(npts)
      real*8 tptmo(npts,nvec)
c
c     Local variables
c
      integer i, j
c     integer n, n1
      integer ipert
      integer ipt
c     integer iatm
      integer ic
      integer iki, ikj
      real*8 tmp
      integer ao_tag
      parameter(ao_tag=1)
c
c     Code
c
      call aclear_dp(drho,npts*npert,0.0d0)
      call aclear_dp(dgrho,npts*npert*3,0.0d0)
c
c     Construct dot product of grad rho and grad mo
c
      do i=1,nvec
         do ipt=1,npts
            gmorho(ipt,i)=0.5d0*
     &      (grho(ipt,1,1)*amo_grad(ipt,i,1)
     &      +grho(ipt,1,2)*amo_grad(ipt,i,2)
     &      +grho(ipt,1,3)*amo_grad(ipt,i,3))
         enddo ! ipt
      enddo ! i
      do ipt=1,npts
         t1(ipt)=wt(ipt)*xc_dvpt(ipt,igaa)
         t2(ipt)=wt(ipt)*xc_vpt(ipt,ira)
      enddo
      do i=1,nvec
         do ipt=1,npts
            tptmo(ipt,i)=t2(ipt)*amo_val(ipt,i)
     &                  +t1(ipt)*gmorho(ipt,i)
         enddo
      enddo
c
c     Construct the derivatives of the density (times 0.5)
c
      iki=0
      do i=1,naocc
         call dmo_val(npts,nao,nvec,npert,ngridcent,mxp,i,
     &                bfng_val,avec,dmo)
         do ipert=1,npert
            do ipt=1,npts
               drho(ipt,ipert)=drho(ipt,ipert)
     &                        +amo_val(ipt,i)*dmo(ipt,ipert)
            enddo ! ipt
         enddo ! ipert
         do ic=1,3
            do ipert=1,npert
               do ipt=1,npts
                  dgrho(ipt,ipert,ic)=dgrho(ipt,ipert,ic)
     &            +amo_grad(ipt,i,ic)*dmo(ipt,ipert)
               enddo ! ipt
            enddo ! ipert
         enddo ! ic

         iki=iki+i-1
         do ipert=1,npert
            do j=1,i-1
               do ipt=1,npts
                  qa(iki+j,ipert)=qa(iki+j,ipert)
c    &            +wt(ipt)*xc_vpt(ipt,ira)*dmo(ipt,ipert)*amo_val(ipt,j)
c    &            +wt(ipt)*xc_dvpt(ipt,igaa)*
c    &             dmo(ipt,ipert)*gmorho(ipt,j)
     &            +dmo(ipt,ipert)*tptmo(ipt,j)
               enddo ! ipt
            enddo ! j
            tmp=0.0d0
            do ipt=1,npts
               tmp=tmp
c    &         +2*wt(ipt)*xc_vpt(ipt,ira)*dmo(ipt,ipert)*amo_val(ipt,i)
c    &         +2*wt(ipt)*xc_dvpt(ipt,igaa)*dmo(ipt,ipert)*gmorho(ipt,i)
     &         +dmo(ipt,ipert)*tptmo(ipt,i)
            enddo ! ipt
            qa(iki+i,ipert)=qa(iki+i,ipert)+2.0d0*tmp
            ikj=iki
            do j=i+1,nvec
               ikj=ikj+j-1
               do ipt=1,npts
                  qa(ikj+i,ipert)=qa(ikj+i,ipert)
c    &            +wt(ipt)*xc_vpt(ipt,ira)*dmo(ipt,ipert)*amo_val(ipt,j)
c    &            +wt(ipt)*xc_dvpt(ipt,igaa)*
c    &             dmo(ipt,ipert)*gmorho(ipt,j)
     &            +dmo(ipt,ipert)*tptmo(ipt,j)
               enddo
            enddo ! j
         enddo ! ipert
      enddo ! i
c
      iki=naocc*(naocc-1)/2
      do i=naocc+1,nvec
         call dmo_val(npts,nao,nvec,npert,ngridcent,mxp,i,
     &                bfng_val,avec,dmo)
         iki=iki+i-1
         do ipert=1,npert
            do j=1,naocc
               do ipt=1,npts
                  qa(iki+j,ipert)=qa(iki+j,ipert)
c    &            +wt(ipt)*xc_vpt(ipt,ira)*dmo(ipt,ipert)*amo_val(ipt,j)
c    &            +wt(ipt)*xc_dvpt(ipt,igaa)*
c    &             dmo(ipt,ipert)*gmorho(ipt,j)
     &            +dmo(ipt,ipert)*tptmo(ipt,j)
               enddo ! ipt
            enddo ! j
         enddo ! ipert
      enddo
c
      do ipt=1,npts
         t2(ipt)=4*wt(ipt)*xc_dhpt(ipt,iragaa)
         t3(ipt)=2*wt(ipt)*xc_hpt(ipt,irara)
      enddo
      do i=1,nvec
         do ipt=1,npts
            tptmo(ipt,i)=t3(ipt)*amo_val(ipt,i)
     &                  +t2(ipt)*gmorho(ipt,i)
         enddo
      enddo
      do ipert=1,npert
         iki=0
         do i=1,nvec
c           iki=iky(i)
            iki=iki+i-1
            do j=1,min(i,naocc)
               do ipt=1,npts
                  qa(iki+j,ipert)=qa(iki+j,ipert)
c    &            +2*wt(ipt)*xc_hpt(ipt,irara)*drho(ipt,ipert)*
c    &             amo_val(ipt,i)*amo_val(ipt,j)
c    &            +4*wt(ipt)*xc_dhpt(ipt,iragaa)*drho(ipt,ipert)*
c    &              gmorho(ipt,i)*amo_val(ipt,j)
     &            +tptmo(ipt,i)*amo_val(ipt,j)*drho(ipt,ipert)
               enddo ! ipt
            enddo ! j
         enddo ! i
      enddo ! ipert
      do i=1,nvec
         do ipt=1,npts
            tptmo(ipt,i)=t2(ipt)*amo_val(ipt,i)
         enddo
      enddo
      do ipert=1,npert
         iki=0
         do i=1,nvec
c           iki=iky(i)
            iki=iki+i-1
            do j=1,min(i,naocc)
               do ipt=1,npts
                  qa(iki+j,ipert)=qa(iki+j,ipert)
c    &            +4*wt(ipt)*xc_dhpt(ipt,iragaa)*drho(ipt,ipert)*
c    &              gmorho(ipt,j)*amo_val(ipt,i)
     &            +tptmo(ipt,i)*gmorho(ipt,j)*drho(ipt,ipert)
               enddo ! ipt
            enddo ! j
         enddo ! i
      enddo ! ipert
c
      do i=1,nvec
         do ipt=1,npts
            tptmo(ipt,i)=t1(ipt)*amo_val(ipt,i)
         enddo
      enddo
      iki=0
      do i=1,naocc
         call dgmo_val(npts,nao,nvec,npert,ngridcent,mxp,i,
     &                 bfn_hess,avec,dgmo)
         do ic=1,3
            do ipert=1,npert
               do ipt=1,npts
                  dgrho(ipt,ipert,ic)=dgrho(ipt,ipert,ic)
     &            +amo_val(ipt,i)*dgmo(ipt,ipert,ic)
               enddo ! ipt
            enddo ! ipert
         enddo ! ic
c
         iki=iki+i-1
         do ipert=1,npert
            do ipt=1,npts
               grdgmo(ipt) = grho(ipt,1,1)*dgmo(ipt,ipert,1)
     &                     + grho(ipt,1,2)*dgmo(ipt,ipert,2)
     &                     + grho(ipt,1,3)*dgmo(ipt,ipert,3)
               grdgmo(ipt) = 0.5d0*grdgmo(ipt)
            enddo
            do j=1,i-1
               do ipt=1,npts
                  qa(iki+j,ipert)=qa(iki+j,ipert)
c    &            +wt(ipt)*xc_dvpt(ipt,igaa)*grdgmo(ipt)*amo_val(ipt,j)
     &            +tptmo(ipt,j)*grdgmo(ipt)
               enddo ! ipt
            enddo ! j
            tmp=0.0d0
            do ipt=1,npts
               tmp=tmp
c    &         +2*wt(ipt)*xc_dvpt(ipt,igaa)*grdgmo(ipt)*amo_val(ipt,i)
     &            +tptmo(ipt,i)*grdgmo(ipt)
            enddo ! ipt
            qa(iki+i,ipert)=qa(iki+i,ipert)+2.0d0*tmp
            ikj=iki
            do j=i+1,nvec
               ikj=ikj+j-1
               do ipt=1,npts
                  qa(ikj+i,ipert)=qa(ikj+i,ipert)
c    &            +wt(ipt)*xc_dvpt(ipt,igaa)*grdgmo(ipt)*amo_val(ipt,j)
     &            +tptmo(ipt,j)*grdgmo(ipt)
               enddo ! ipt
            enddo ! j
         enddo ! npert
      enddo ! i
c
      iki=naocc*(naocc-1)/2
      do i=naocc+1,nvec
         call dgmo_val(npts,nao,nvec,npert,ngridcent,mxp,i,
     &                 bfn_hess,avec,dgmo)
c        iki=iky(i)
         iki=iki+i-1
         do ipert=1,npert
            do ipt=1,npts
               grdgmo(ipt) = grho(ipt,1,1)*dgmo(ipt,ipert,1)
     &                     + grho(ipt,1,2)*dgmo(ipt,ipert,2)
     &                     + grho(ipt,1,3)*dgmo(ipt,ipert,3)
               grdgmo(ipt) = 0.5d0*grdgmo(ipt)
            enddo
            do j=1,naocc
               do ipt=1,npts
                  qa(iki+j,ipert)=qa(iki+j,ipert)
c    &            +wt(ipt)*xc_dvpt(ipt,igaa)*grdgmo(ipt)*amo_val(ipt,j)
     &            +tptmo(ipt,j)*grdgmo(ipt)
               enddo ! ipt
            enddo ! j
         enddo ! npert
      enddo ! i
c
      do ipt=1,npts
         t3(ipt)=4*wt(ipt)*xc_dhpt(ipt,igaagaa)
      enddo
      do i=1,nvec
         do ipt=1,npts
            tptmo(ipt,i)=t2(ipt)*amo_val(ipt,i)
     &                  +t3(ipt)*gmorho(ipt,i)
         enddo
      enddo
      do ipert=1,npert
         do ipt=1,npts
            grdgmo(ipt) = grho(ipt,1,1)*dgrho(ipt,ipert,1)
     &                  + grho(ipt,1,2)*dgrho(ipt,ipert,2)
     &                  + grho(ipt,1,3)*dgrho(ipt,ipert,3)
            grdgmo(ipt) = 0.5d0*grdgmo(ipt)
         enddo ! ipt
         iki=0
         do i=1,nvec
c           iki=iky(i)
            iki=iki+i-1
            do j=1,min(i,naocc)
               do ipt=1,npts
                  qa(iki+j,ipert)=qa(iki+j,ipert)
c    &            +4*wt(ipt)*xc_dhpt(ipt,iragaa)*grdgmo(ipt)*
c    &             amo_val(ipt,i)*amo_val(ipt,j)
c    &            +4*wt(ipt)*xc_dhpt(ipt,igaagaa)*grdgmo(ipt)*
c    &             gmorho(ipt,i)*amo_val(ipt,j)
     &            +tptmo(ipt,i)*amo_val(ipt,j)*grdgmo(ipt)
               enddo ! ipt
            enddo ! j
         enddo ! i
      enddo ! ipert
      do i=1,nvec
         do ipt=1,npts
            tptmo(ipt,i)=t3(ipt)*amo_val(ipt,i)
         enddo
      enddo
      do ipert=1,npert
         do ipt=1,npts
            grdgmo(ipt) = grho(ipt,1,1)*dgrho(ipt,ipert,1)
     &                  + grho(ipt,1,2)*dgrho(ipt,ipert,2)
     &                  + grho(ipt,1,3)*dgrho(ipt,ipert,3)
            grdgmo(ipt) = 0.5d0*grdgmo(ipt)
         enddo ! ipt
         iki=0
         do i=1,nvec
c           iki=iky(i)
            iki=iki+i-1
            do j=1,min(i,naocc)
               do ipt=1,npts
                  qa(iki+j,ipert)=qa(iki+j,ipert)
c    &            +4*wt(ipt)*xc_dhpt(ipt,igaagaa)*grdgmo(ipt)*
c    &             gmorho(ipt,j)*amo_val(ipt,i)
     &            +tptmo(ipt,i)*gmorho(ipt,j)*grdgmo(ipt)
               enddo ! ipt
            enddo ! j
         enddo ! i
      enddo ! ipert
c
      do ipt=1,npts
         t1(ipt)=2*t1(ipt)
      enddo
      do i=1,nvec
         do ipt=1,npts
            tptmo(ipt,i)=t1(ipt)*amo_val(ipt,i)
         enddo
      enddo
      do ipert=1,npert
         do i=1,nvec 
            do ipt=1,npts
               gmorho(ipt,i) = dgrho(ipt,ipert,1)*amo_grad(ipt,i,1)
     &                       + dgrho(ipt,ipert,2)*amo_grad(ipt,i,2)
     &                       + dgrho(ipt,ipert,3)*amo_grad(ipt,i,3)
            enddo ! ipt
         enddo ! i
         iki=0
         do i=1,nvec
c           iki=iky(i)
            iki=iki+i-1
            do j=1,min(i,naocc)
               do ipt=1,npts
                  qa(iki+j,ipert)=qa(iki+j,ipert)
c    &            +wt(ipt)*xc_dvpt(ipt,igaa)*
c    &             (gmorho(ipt,i)*amo_val(ipt,j)+
c    &              gmorho(ipt,j)*amo_val(ipt,i))
     &            +gmorho(ipt,i)*tptmo(ipt,j)
     &            +gmorho(ipt,j)*tptmo(ipt,i)
               enddo ! ipt
            enddo ! j
         enddo ! i
      enddo ! ipert
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine rks_dksm_exp_mo_gwt(npts,nao,nvec,naocc,npert,
     &           ngridcent,mxp,gwt_avail_sw,latm,
     &           wt,gwt,xc_vpt,xc_hpt,amo_val,bfng_val,
     &           avec,dmo,drho,tmo,qa)
      implicit none
c
c     Compute the explicit derivative contribution to the Kohn-Sham
c     matrix in MO-basis for the local density closed shell case.
c
c     The tactics are to construct to first and second derivatives
c     of the density with respect to the nuclear coordinates on the
c     fly first and subsequently consume them. Derivates of the MO's
c     are generated per orbital so as to minimise memory requirements.
c
c     Inputs
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
C *Parameter module for CCP1/DFT            
C *Memory information		
C *----------------
C * max_block 	- maximum number of blocks which can be allocate
      integer max_block
      parameter(max_block=20)
C *							
C *Basis sets information                              
C *----------------------				
C *max_tag 	- maximum number of basis sets which can be alloc
C *max_atype	- maximum number of centre types	
C *max_gtype	- maximum number of grid types (at least every element 
C                 has different grid type)
C *max_grids    - maximum number of grids (different terms may have
C                 different grids, i.e. CPKS equations use a different
C                 grid than the one used for the KS-matrix)
C *max_shel	- maximum number of shells on centre
C *max_prm	- maximum number of primitives for any given centre
C *maxL		- maximum angular momentum allowed	
C *max_func	- maximum number of basis functions for any centre
      integer max_tag,max_atype,max_shel,max_prm,max_ang
      integer max_gtype,max_grids
      parameter(max_tag=3,max_atype=30,max_shel=500,max_prm=5000)
      parameter(max_ang=5,max_gtype=10,max_grids=2)
C *								
C *Geometry information					
C *--------------------					
C *max_atom	- maximum number of atoms in system
      integer max_atom
      parameter(max_atom=750)
C *
C *Accuracy information				
C ---------------------			
C *global_accuracy 	- global accuracy	
      real*8  global_accuracy
      parameter(global_accuracy=1.0d-14)

C *Grid information
c *----------------
c *maxradzn - The maximum number of radial zones. 
      integer maxgpt,maxrad,maxfpt,maxang,maxradzn,maxtablerows
      parameter(maxgpt=2900,maxrad=50,maxang=302,maxfpt=100)
      parameter(maxradzn=35,maxtablerows=7)
c
c     The parameters defined here label the 1st, 2nd and 3rd derivative
c     components of the basis functions.
c
c     2nd derivatives:
c
      integer hxx, hyy, hzz, hxy, hxz, hyz
c
      parameter (hxx = 1)
      parameter (hyy = 2)
      parameter (hzz = 3)
      parameter (hxy = 4)
      parameter (hxz = 5)
      parameter (hyz = 6)
c
      integer hcc(3,3)
      data hcc/hxx,hxy,hxz,
     &         hxy,hyy,hyz,
     &         hxz,hyz,hzz/
c
c     3rd derivatives:
c
      integer txxx, txxy, txxz, txyy, txyz
      integer txzz, tyyy, tyyz, tyzz, tzzz
c
      parameter (txxx =  1)
      parameter (txxy =  2)
      parameter (txxz =  3)
      parameter (txyy =  4)
      parameter (txyz =  5)
      parameter (txzz =  6)
      parameter (tyyy =  7)
      parameter (tyyz =  8)
      parameter (tyzz =  9)
      parameter (tzzz = 10)
c
      integer tc3(3,3,3)
      data tc3/txxx,txxy,txxz,
     &         txxy,txyy,txyz,
     &         txxz,txyz,txzz,
     &         txxy,txyy,txyz,
     &         txyy,tyyy,tyyz,
     &         txyz,tyyz,tyzz,
     &         txxz,txyz,txzz,
     &         txyz,tyyz,tyzz,
     &         txzz,tyzz,tzzz/
c
c     The parameters defined here label the 1st and 2nd derivative
c     components of the density functional. Because a functional
c     may be purely local or of a GGA type the components are grouped
c     into derivative with respect to the density only and others.
c
c     1st derivatives:
c
      integer ira, irb
      integer igaa, igab, igbb
      integer ita, itb
c
c     2nd derivatives:
c
      integer irara, irarb, irbrb
      integer iragaa, iragab, iragbb, irbgaa, irbgab, irbgbb
      integer igaagaa, igaagab, igaagbb, igabgab, igabgbb, igbbgbb
      integer irata, iratb, irbta, irbtb
      integer igaata, igaatb, igabta, igabtb, igbbta, igbbtb
      integer itata, itatb, itbtb
c
      parameter (ira     = 1)
      parameter (irb     = 2)
c
      parameter (igaa    = 1)
      parameter (igab    = 2)
      parameter (igbb    = 3)
c
      parameter (ita     = 1)
      parameter (itb     = 2)
c
c
      parameter (irara   = 1)
      parameter (irarb   = 2)
      parameter (irbrb   = 3)
c
      parameter (iragaa  = 1)
      parameter (iragab  = 3)
      parameter (iragbb  = 4)
      parameter (irbgaa  = 5)
      parameter (irbgab  = 6)
      parameter (irbgbb  = 7)
      parameter (igaagaa = 2)
      parameter (igaagab = 8)
      parameter (igaagbb = 9)
      parameter (igabgab = 10)
      parameter (igabgbb = 11)
      parameter (igbbgbb = 12)
c
      parameter (irata   = 1)
      parameter (iratb   = 3)
      parameter (irbta   = 4)
      parameter (irbtb   = 5)
      parameter (itata   = 2)
      parameter (itatb   = 6)
      parameter (itbtb   = 7)
c
      parameter (igaata  = 1)
      parameter (igaatb  = 2)
      parameter (igabta  = 3)
      parameter (igabtb  = 4)
      parameter (igbbta  = 5)
      parameter (igbbtb  = 6)
c
      integer npts
      integer nao
      integer nvec
      integer naocc 
      integer npert
      integer ngridcent
      integer mxp 
      logical gwt_avail_sw
      integer latm
      real*8 wt(npts)
      real*8 gwt(3,mxp,ngridcent)
      real*8 xc_vpt(mxp,2)
      real*8 xc_hpt(mxp,3)
      real*8 amo_val(npts,nvec)
      real*8 bfng_val(mxp,nao,3)
c     real*8 bfn_hess(mxp,nao,6)
      real*8 avec(nao,nvec)
c
c     Outputs
c
      real*8 qa(nvec*(nvec+1)/2,npert)
c
c     Workspace
c
      real*8 dmo(npts,npert)
c     real*8 ddmo(npts,6,ngridcent)
      real*8 drho(npts,npert)
      real*8 tmo(npts,nvec)
c     real*8 ddrho(npts,npert*(npert+1)/2)
c
c     Local variables
c
      integer i, j
      integer n, n1
      integer ipert, jpert
      integer ipt
      integer iatm
      integer ic
      integer iki, ikj
      real*8 tmp
      integer ao_tag
      parameter(ao_tag=1)
c
c     Code:
c
c     First the contributions of the derivatives of the weights
c
      if (gwt_avail_sw) then
         do iatm=1,ngridcent
            do ic=1,3
               ipert=3*(iatm-1)+ic
               do j=1,nvec
                  ikj=iky(j)
                  do i=1,min(j,naocc)
                     do ipt=1,npts
                        qa(ikj+i,ipert)=qa(ikj+i,ipert)
     &                  +gwt(ic,ipt,iatm)*xc_vpt(ipt,ira)
     &                  *amo_val(ipt,i)*amo_val(ipt,j)
                     enddo ! ipt
                  enddo ! i
               enddo ! j
            enddo ! ic
         enddo ! iatm
      endif
c
c     Now the derivatives of the Fock operator
c
      call aclear_dp(drho,npts*npert,0.0d0)
      do ipt=1,npts
         dmo(ipt,1)=wt(ipt)*xc_vpt(ipt,ira)
      enddo
      do i=1,nvec
         do ipt=1,npts
            tmo(ipt,i)=dmo(ipt,1)*amo_val(ipt,i)
         enddo
      enddo
c
c     Construct the derivatives of the density (times 0.5)
c
      iki=0
      do i=1,naocc
         call dmo_val(npts,nao,nvec,npert,ngridcent,mxp,i,
     &                bfng_val,avec,dmo)
         do ipert=1,npert
            do ipt=1,npts
               drho(ipt,ipert)=drho(ipt,ipert)
     &                        +amo_val(ipt,i)*dmo(ipt,ipert)
            enddo
         enddo
c
         iki=iki+i-1
         do iatm=1,ngridcent
            if (iatm.ne.latm) then
               do ic=1,3
                  ipert=3*(iatm-1)+ic
                  jpert=3*(latm-1)+ic
                  do j=1,i-1
                     tmp=0.0d0
                     do ipt=1,npts
                        tmp=tmp+dmo(ipt,ipert)*tmo(ipt,j)
                     enddo ! ipt
                     qa(iki+j,ipert)=qa(iki+j,ipert)+tmp
                     qa(iki+j,jpert)=qa(iki+j,jpert)-tmp
                  enddo ! j
                  tmp=0.0d0
                  do ipt=1,npts
                     tmp=tmp+dmo(ipt,ipert)*tmo(ipt,i)
                  enddo ! ipt
                  tmp=2.0d0*tmp
                  qa(iki+i,ipert)=qa(iki+i,ipert)+tmp
                  qa(iki+i,jpert)=qa(iki+i,jpert)-tmp
                  ikj=iki
                  do j=i+1,nvec
                     ikj=ikj+j-1
                     tmp=0.0d0
                     do ipt=1,npts
                        tmp=tmp+dmo(ipt,ipert)*tmo(ipt,j)
                     enddo
                     qa(ikj+i,ipert)=qa(ikj+i,ipert)+tmp
                     qa(ikj+i,jpert)=qa(ikj+i,jpert)-tmp
                  enddo ! j
               enddo ! ic
            endif
         enddo ! iatm
      enddo
c
      n=(naocc-1)*naocc/2
      do i=naocc+1,nvec
         call dmo_val(npts,nao,nvec,npert,ngridcent,mxp,i,
     &                bfng_val,avec,dmo)
         n=n+i-1
         do iatm=1,ngridcent
            if (iatm.ne.latm) then
               do ic=1,3
                  ipert=3*(iatm-1)+ic
                  jpert=3*(latm-1)+ic
                  do j=1,naocc
                     tmp=0.0d0
                     do ipt=1,npts
                        tmp=tmp+dmo(ipt,ipert)*tmo(ipt,j)
                     enddo ! ipt
                     qa(n+j,ipert)=qa(n+j,ipert)+tmp
                     qa(n+j,jpert)=qa(n+j,jpert)-tmp
                  enddo ! j
               enddo ! ic
            endif 
         enddo ! iatm
      enddo
c
c     Add the contributions to the Kohn-Sham matrix
c
      do ipt=1,npts 
         dmo(ipt,1)=2*wt(ipt)*xc_hpt(ipt,irara)
      enddo
      do iatm=1,ngridcent
         if (iatm.ne.latm) then
            do ic=1,3
               ipert=3*(iatm-1)+ic
               jpert=3*(latm-1)+ic
               do ipt=1,npts 
                  dmo(ipt,2)=dmo(ipt,1)*drho(ipt,ipert)
               enddo
               n=0
               do j=1,naocc
                  n1=n
                  n=n+j-1
                  do ipt=1,npts
                     dmo(ipt,3)=dmo(ipt,2)*amo_val(ipt,j)
                  enddo
                  do i=j,nvec
                     n1=n1+i-1
                     tmp=0.0d0
                     do ipt=1,npts
                        tmp=tmp+dmo(ipt,3)*amo_val(ipt,i)
                     enddo ! ipt
                     qa(n1+j,ipert)=qa(n1+j,ipert)+tmp
                     qa(n1+j,jpert)=qa(n1+j,jpert)-tmp
                  enddo ! i
               enddo ! j
            enddo ! ic
         endif
      enddo ! iatm
c 
      end
c
c-----------------------------------------------------------------------
c
      subroutine rks_dksm_exp_mo_gga_gwt(npts,nao,nvec,naocc,npert,
     &           ngridcent,mxp,gwt_avail_sw,latm,wt,gwt,xc_vpt,xc_dvpt,
     &           xc_hpt,xc_dhpt,amo_val,amo_grad,bfng_val,bfn_hess,
     &           avec,grho,dmo,dgmo,drho,dgrho,gmorho,grdgmo,t1,t2,t3,
     &           t4,tptmo,qa)
      implicit none
c
c     Compute the explicit derivative contribution to the Kohn-Sham 
c     matrix in MO-basis for the gradient corrected closed shell case.
c
c     The tactics are to construct to first and second derivatives
c     of the density with respect to the nuclear coordinates on the
c     fly first and subsequently consume them. Derivates of the MO's
c     are generated per orbital so as to minimise memory requirements.
c
c     Inputs
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
C *Parameter module for CCP1/DFT            
C *Memory information		
C *----------------
C * max_block 	- maximum number of blocks which can be allocate
      integer max_block
      parameter(max_block=20)
C *							
C *Basis sets information                              
C *----------------------				
C *max_tag 	- maximum number of basis sets which can be alloc
C *max_atype	- maximum number of centre types	
C *max_gtype	- maximum number of grid types (at least every element 
C                 has different grid type)
C *max_grids    - maximum number of grids (different terms may have
C                 different grids, i.e. CPKS equations use a different
C                 grid than the one used for the KS-matrix)
C *max_shel	- maximum number of shells on centre
C *max_prm	- maximum number of primitives for any given centre
C *maxL		- maximum angular momentum allowed	
C *max_func	- maximum number of basis functions for any centre
      integer max_tag,max_atype,max_shel,max_prm,max_ang
      integer max_gtype,max_grids
      parameter(max_tag=3,max_atype=30,max_shel=500,max_prm=5000)
      parameter(max_ang=5,max_gtype=10,max_grids=2)
C *								
C *Geometry information					
C *--------------------					
C *max_atom	- maximum number of atoms in system
      integer max_atom
      parameter(max_atom=750)
C *
C *Accuracy information				
C ---------------------			
C *global_accuracy 	- global accuracy	
      real*8  global_accuracy
      parameter(global_accuracy=1.0d-14)

C *Grid information
c *----------------
c *maxradzn - The maximum number of radial zones. 
      integer maxgpt,maxrad,maxfpt,maxang,maxradzn,maxtablerows
      parameter(maxgpt=2900,maxrad=50,maxang=302,maxfpt=100)
      parameter(maxradzn=35,maxtablerows=7)
c
c     The parameters defined here label the 1st, 2nd and 3rd derivative
c     components of the basis functions.
c
c     2nd derivatives:
c
      integer hxx, hyy, hzz, hxy, hxz, hyz
c
      parameter (hxx = 1)
      parameter (hyy = 2)
      parameter (hzz = 3)
      parameter (hxy = 4)
      parameter (hxz = 5)
      parameter (hyz = 6)
c
      integer hcc(3,3)
      data hcc/hxx,hxy,hxz,
     &         hxy,hyy,hyz,
     &         hxz,hyz,hzz/
c
c     3rd derivatives:
c
      integer txxx, txxy, txxz, txyy, txyz
      integer txzz, tyyy, tyyz, tyzz, tzzz
c
      parameter (txxx =  1)
      parameter (txxy =  2)
      parameter (txxz =  3)
      parameter (txyy =  4)
      parameter (txyz =  5)
      parameter (txzz =  6)
      parameter (tyyy =  7)
      parameter (tyyz =  8)
      parameter (tyzz =  9)
      parameter (tzzz = 10)
c
      integer tc3(3,3,3)
      data tc3/txxx,txxy,txxz,
     &         txxy,txyy,txyz,
     &         txxz,txyz,txzz,
     &         txxy,txyy,txyz,
     &         txyy,tyyy,tyyz,
     &         txyz,tyyz,tyzz,
     &         txxz,txyz,txzz,
     &         txyz,tyyz,tyzz,
     &         txzz,tyzz,tzzz/
c
c     The parameters defined here label the 1st and 2nd derivative
c     components of the density functional. Because a functional
c     may be purely local or of a GGA type the components are grouped
c     into derivative with respect to the density only and others.
c
c     1st derivatives:
c
      integer ira, irb
      integer igaa, igab, igbb
      integer ita, itb
c
c     2nd derivatives:
c
      integer irara, irarb, irbrb
      integer iragaa, iragab, iragbb, irbgaa, irbgab, irbgbb
      integer igaagaa, igaagab, igaagbb, igabgab, igabgbb, igbbgbb
      integer irata, iratb, irbta, irbtb
      integer igaata, igaatb, igabta, igabtb, igbbta, igbbtb
      integer itata, itatb, itbtb
c
      parameter (ira     = 1)
      parameter (irb     = 2)
c
      parameter (igaa    = 1)
      parameter (igab    = 2)
      parameter (igbb    = 3)
c
      parameter (ita     = 1)
      parameter (itb     = 2)
c
c
      parameter (irara   = 1)
      parameter (irarb   = 2)
      parameter (irbrb   = 3)
c
      parameter (iragaa  = 1)
      parameter (iragab  = 3)
      parameter (iragbb  = 4)
      parameter (irbgaa  = 5)
      parameter (irbgab  = 6)
      parameter (irbgbb  = 7)
      parameter (igaagaa = 2)
      parameter (igaagab = 8)
      parameter (igaagbb = 9)
      parameter (igabgab = 10)
      parameter (igabgbb = 11)
      parameter (igbbgbb = 12)
c
      parameter (irata   = 1)
      parameter (iratb   = 3)
      parameter (irbta   = 4)
      parameter (irbtb   = 5)
      parameter (itata   = 2)
      parameter (itatb   = 6)
      parameter (itbtb   = 7)
c
      parameter (igaata  = 1)
      parameter (igaatb  = 2)
      parameter (igabta  = 3)
      parameter (igabtb  = 4)
      parameter (igbbta  = 5)
      parameter (igbbtb  = 6)
c
      integer npts
      integer nao
      integer nvec
      integer naocc 
      integer npert
      integer ngridcent
      integer mxp 
      logical gwt_avail_sw
      integer latm
      real*8 wt(npts)
      real*8 gwt(3,mxp,ngridcent)
      real*8 xc_vpt(mxp,2)
      real*8 xc_dvpt(mxp,3)
      real*8 xc_hpt(mxp,3)
      real*8 xc_dhpt(mxp,12)
      real*8 amo_val(npts,nvec)
      real*8 amo_grad(npts,nvec,3)
      real*8 bfng_val(mxp,nao,3)
      real*8 bfn_hess(mxp,nao,6)
c     real*8 bfn_3rd(mxp,nao,10)
      real*8 avec(nao,nvec)
      real*8 grho(mxp,2,3)
c
c     Outputs
c
      real*8 qa(nao*(nao+1)/2,npert)
c
c     Workspace
c
      real*8 dmo(npts,npert)
      real*8 dgmo(npts,npert,3)
c     real*8 ddmo(npts,6,ngridcent)
c     real*8 ddgmo(npts,10,ngridcent)
      real*8 drho(npts,npert)
c     real*8 ddrho(npts,npert*(npert+1)/2)
      real*8 dgrho(npts,npert,3)
c     real*8 ddgrho(npts,npert*(npert+1)/2,3)
      real*8 gmorho(npts,nvec)
      real*8 grdgmo(npts)
      real*8 t1(npts)
      real*8 t2(npts)
      real*8 t3(npts)
      real*8 t4(npts)
      real*8 tptmo(npts,nvec)
c
c     Local variables
c
      integer i, j
c     integer n, n1
      integer ipert, jpert
      integer ipt
      integer iatm
      integer ic
      integer iki, ikj
      real*8 tmp
      integer ao_tag
      parameter(ao_tag=1)
c
c     Code
c
c     Construct dot product of grad rho and grad mo
c
      do i=1,nvec
         do ipt=1,npts
            gmorho(ipt,i)=
     &      (grho(ipt,1,1)*amo_grad(ipt,i,1)
     &      +grho(ipt,1,2)*amo_grad(ipt,i,2)
     &      +grho(ipt,1,3)*amo_grad(ipt,i,3))
            gmorho(ipt,i)=0.5d0*gmorho(ipt,i)
         enddo ! ipt
      enddo ! i
c
c     First add the contributions of the derivatives of the weights
c
      if (gwt_avail_sw) then
         do iatm=1,ngridcent
            do ic=1,3
               ipert=3*(iatm-1)+ic
               do j=1,nvec
                  ikj=iky(j)
                  do i=1,min(j,naocc)
                     do ipt=1,npts
                        qa(ikj+i,ipert)=qa(ikj+i,ipert)
     &                  +gwt(ic,ipt,iatm)*(
     &                   xc_vpt(ipt,ira)*amo_val(ipt,i)*amo_val(ipt,j)
     &                  +xc_dvpt(ipt,igaa)*gmorho(ipt,i)*amo_val(ipt,j)
     &                  +xc_dvpt(ipt,igaa)*amo_val(ipt,i)*gmorho(ipt,j))
                     enddo ! ipt
                  enddo ! i
               enddo ! j
            enddo ! ic
         enddo ! iatm
      endif
c
c     Now the derivatives of the Fock operator
c
      call aclear_dp(drho,npts*npert,0.0d0)
      call aclear_dp(dgrho,npts*npert*3,0.0d0)
c
      do ipt=1,npts
         t1(ipt)=wt(ipt)*xc_dvpt(ipt,igaa)
         t2(ipt)=wt(ipt)*xc_vpt(ipt,ira)
      enddo
      do i=1,nvec
         do ipt=1,npts
            tptmo(ipt,i)=t2(ipt)*amo_val(ipt,i)
     &                  +t1(ipt)*gmorho(ipt,i)
         enddo
      enddo
c
c     Construct the derivatives of the density (times 0.5)
c
      iki=0
      do i=1,naocc
         call dmo_val(npts,nao,nvec,npert,ngridcent,mxp,i,
     &                bfng_val,avec,dmo)
         do ipert=1,npert
            do ipt=1,npts
               drho(ipt,ipert)=drho(ipt,ipert)
     &                        +amo_val(ipt,i)*dmo(ipt,ipert)
            enddo ! ipt
         enddo ! ipert
         do ic=1,3
            do ipert=1,npert
               do ipt=1,npts
                  dgrho(ipt,ipert,ic)=dgrho(ipt,ipert,ic)
     &            +amo_grad(ipt,i,ic)*dmo(ipt,ipert)
               enddo ! ipt
            enddo ! ipert
         enddo ! ic

         iki=iki+i-1
         do iatm=1,ngridcent
            if (iatm.ne.latm) then
               do ic=1,3
                  ipert=3*(iatm-1)+ic
                  jpert=3*(latm-1)+ic
                  do j=1,i-1
                     tmp=0.0d0
                     do ipt=1,npts
                        tmp=tmp+dmo(ipt,ipert)*tptmo(ipt,j)
                     enddo ! ipt
                     qa(iki+j,ipert)=qa(iki+j,ipert)+tmp
                     qa(iki+j,jpert)=qa(iki+j,jpert)-tmp
                  enddo ! j
                  tmp=0.0d0
                  do ipt=1,npts
                     tmp=tmp+dmo(ipt,ipert)*tptmo(ipt,i)
                  enddo ! ipt
                  tmp=2.0d0*tmp
                  qa(iki+i,ipert)=qa(iki+i,ipert)+tmp
                  qa(iki+i,jpert)=qa(iki+i,jpert)-tmp
                  ikj=iki
                  do j=i+1,nvec
                     ikj=ikj+j-1
                     tmp=0.0d0
                     do ipt=1,npts
                        tmp=tmp+dmo(ipt,ipert)*tptmo(ipt,j)
                     enddo
                     qa(ikj+i,ipert)=qa(ikj+i,ipert)+tmp
                     qa(ikj+i,jpert)=qa(ikj+i,jpert)-tmp
                  enddo ! j
               enddo ! ic
            endif
         enddo ! iatm
      enddo ! i
c
      iki=naocc*(naocc-1)/2
      do i=naocc+1,nvec
         call dmo_val(npts,nao,nvec,npert,ngridcent,mxp,i,
     &                bfng_val,avec,dmo)
         iki=iki+i-1
         do iatm=1,ngridcent
            if (iatm.ne.latm) then
               do ic=1,3
                  ipert=3*(iatm-1)+ic
                  jpert=3*(latm-1)+ic
                  do j=1,naocc
                     tmp=0.0d0
                     do ipt=1,npts
                        tmp=tmp+dmo(ipt,ipert)*tptmo(ipt,j)
                     enddo ! ipt
                     qa(iki+j,ipert)=qa(iki+j,ipert)+tmp
                     qa(iki+j,jpert)=qa(iki+j,jpert)-tmp
                  enddo ! j
               enddo ! ic
            endif
         enddo ! iatm
      enddo
c
      do ipt=1,npts
         t2(ipt)=4*wt(ipt)*xc_dhpt(ipt,iragaa)
         t3(ipt)=2*wt(ipt)*xc_hpt(ipt,irara)
      enddo
      do i=1,nvec
         do ipt=1,npts
            tptmo(ipt,i)=t3(ipt)*amo_val(ipt,i)
     &                  +t2(ipt)*gmorho(ipt,i)
         enddo
      enddo
      do iatm=1,ngridcent
         if (iatm.ne.latm) then
            do ic=1,3
               ipert=3*(iatm-1)+ic
               jpert=3*(latm-1)+ic
               iki=0
               do i=1,nvec
                  iki=iki+i-1
                  do j=1,min(i,naocc)
                     tmp=0.0d0
                     do ipt=1,npts
                        tmp=tmp+tptmo(ipt,i)*amo_val(ipt,j)*
     &                          drho(ipt,ipert)
                     enddo ! ipt
                     qa(iki+j,ipert)=qa(iki+j,ipert)+tmp
                     qa(iki+j,jpert)=qa(iki+j,jpert)-tmp
                  enddo ! j
               enddo ! i
            enddo ! ic
         endif
      enddo ! iatm
      do i=1,nvec
         do ipt=1,npts
            tptmo(ipt,i)=t2(ipt)*amo_val(ipt,i)
         enddo
      enddo
      do iatm=1,ngridcent
         if (iatm.ne.latm) then
            do ic=1,3
               ipert=3*(iatm-1)+ic
               jpert=3*(latm-1)+ic
               iki=0
               do i=1,nvec
                  iki=iki+i-1
                  do j=1,min(i,naocc)
                     tmp=0.0d0
                     do ipt=1,npts
                        tmp=tmp+tptmo(ipt,i)*gmorho(ipt,j)*
     &                          drho(ipt,ipert)
                     enddo ! ipt
                     qa(iki+j,ipert)=qa(iki+j,ipert)+tmp
                     qa(iki+j,jpert)=qa(iki+j,jpert)-tmp
                  enddo ! j
               enddo ! i
            enddo ! ic
         endif
      enddo ! iatm
c
      do i=1,nvec
         do ipt=1,npts
            tptmo(ipt,i)=t1(ipt)*amo_val(ipt,i)
         enddo
      enddo
      iki=0
      do i=1,naocc
         call dgmo_val(npts,nao,nvec,npert,ngridcent,mxp,i,
     &                 bfn_hess,avec,dgmo)
         do ic=1,3
            do ipert=1,npert
               do ipt=1,npts
                  dgrho(ipt,ipert,ic)=dgrho(ipt,ipert,ic)
     &            +amo_val(ipt,i)*dgmo(ipt,ipert,ic)
               enddo ! ipt
            enddo ! ipert
         enddo ! ic
c
         iki=iki+i-1
         do iatm=1,ngridcent
            if (iatm.ne.latm) then
               do ic=1,3
                  ipert=3*(iatm-1)+ic
                  jpert=3*(latm-1)+ic
                  do ipt=1,npts
                     grdgmo(ipt) = grho(ipt,1,1)*dgmo(ipt,ipert,1)
     &                           + grho(ipt,1,2)*dgmo(ipt,ipert,2)
     &                           + grho(ipt,1,3)*dgmo(ipt,ipert,3)
                     grdgmo(ipt) = 0.5d0*grdgmo(ipt)
                  enddo
                  do j=1,i-1
                     tmp=0.0d0
                     do ipt=1,npts
                        tmp=tmp+tptmo(ipt,j)*grdgmo(ipt)
                     enddo ! ipt
                     qa(iki+j,ipert)=qa(iki+j,ipert)+tmp
                     qa(iki+j,jpert)=qa(iki+j,jpert)-tmp
                  enddo ! j
                  tmp=0.0d0
                  do ipt=1,npts
                     tmp=tmp+tptmo(ipt,i)*grdgmo(ipt)
                  enddo ! ipt
                  tmp=2.0d0*tmp
                  qa(iki+i,ipert)=qa(iki+i,ipert)+tmp
                  qa(iki+i,jpert)=qa(iki+i,jpert)-tmp
                  ikj=iki
                  do j=i+1,nvec
                     ikj=ikj+j-1
                     tmp=0.0d0
                     do ipt=1,npts
                        tmp=tmp+tptmo(ipt,j)*grdgmo(ipt)
                     enddo ! ipt
                     qa(ikj+i,ipert)=qa(ikj+i,ipert)+tmp
                     qa(ikj+i,jpert)=qa(ikj+i,jpert)-tmp
                  enddo ! j
               enddo ! ic
            endif
         enddo ! iatm
      enddo ! i
c
      iki=naocc*(naocc-1)/2
      do i=naocc+1,nvec
         call dgmo_val(npts,nao,nvec,npert,ngridcent,mxp,i,
     &                 bfn_hess,avec,dgmo)
         iki=iki+i-1
         do iatm=1,ngridcent
            if (iatm.ne.latm) then
               do ic=1,3
                  ipert=3*(iatm-1)+ic
                  jpert=3*(latm-1)+ic
                  do ipt=1,npts
                     grdgmo(ipt) = grho(ipt,1,1)*dgmo(ipt,ipert,1)
     &                           + grho(ipt,1,2)*dgmo(ipt,ipert,2)
     &                           + grho(ipt,1,3)*dgmo(ipt,ipert,3)
                     grdgmo(ipt) = 0.5d0*grdgmo(ipt)
                  enddo
                  do j=1,naocc
                     tmp=0.0d0
                     do ipt=1,npts
                        tmp=tmp+tptmo(ipt,j)*grdgmo(ipt)
                     enddo ! ipt
                     qa(iki+j,ipert)=qa(iki+j,ipert)+tmp
                     qa(iki+j,jpert)=qa(iki+j,jpert)-tmp
                  enddo ! j
               enddo ! ic
            endif
         enddo ! iatm
      enddo ! i
c
      do ipt=1,npts
         t3(ipt)=4*wt(ipt)*xc_dhpt(ipt,igaagaa)
      enddo
      do i=1,nvec
         do ipt=1,npts
            tptmo(ipt,i)=t2(ipt)*amo_val(ipt,i)
     &                  +t3(ipt)*gmorho(ipt,i)
         enddo
      enddo
      do iatm=1,ngridcent
         if (iatm.ne.latm) then
            do ic=1,3
               ipert=3*(iatm-1)+ic
               jpert=3*(latm-1)+ic
               do ipt=1,npts
                  grdgmo(ipt) = grho(ipt,1,1)*dgrho(ipt,ipert,1)
     &                        + grho(ipt,1,2)*dgrho(ipt,ipert,2)
     &                        + grho(ipt,1,3)*dgrho(ipt,ipert,3)
                  grdgmo(ipt) = 0.5d0*grdgmo(ipt)
               enddo ! ipt
               iki=0
               do i=1,nvec
                  iki=iki+i-1
                  do j=1,min(i,naocc)
                     tmp=0.0d0
                     do ipt=1,npts
                        tmp=tmp+tptmo(ipt,i)*amo_val(ipt,j)*grdgmo(ipt)
                     enddo ! ipt
                     qa(iki+j,ipert)=qa(iki+j,ipert)+tmp
                     qa(iki+j,jpert)=qa(iki+j,jpert)-tmp
                  enddo ! j
               enddo ! i
            enddo ! ic
         endif
      enddo ! iatm
      do i=1,nvec
         do ipt=1,npts
            tptmo(ipt,i)=t3(ipt)*amo_val(ipt,i)
         enddo
      enddo
      do iatm=1,ngridcent
         if (iatm.ne.latm) then
            do ic=1,3
               ipert=3*(iatm-1)+ic
               jpert=3*(latm-1)+ic
               do ipt=1,npts
                  grdgmo(ipt) = grho(ipt,1,1)*dgrho(ipt,ipert,1)
     &                        + grho(ipt,1,2)*dgrho(ipt,ipert,2)
     &                        + grho(ipt,1,3)*dgrho(ipt,ipert,3)
                  grdgmo(ipt) = 0.5d0*grdgmo(ipt)
               enddo ! ipt
               iki=0
               do i=1,nvec
                  iki=iki+i-1
                  do j=1,min(i,naocc)
                     tmp=0.0d0
                     do ipt=1,npts
                        tmp=tmp+tptmo(ipt,i)*gmorho(ipt,j)*grdgmo(ipt)
                     enddo ! ipt
                     qa(iki+j,ipert)=qa(iki+j,ipert)+tmp
                     qa(iki+j,jpert)=qa(iki+j,jpert)-tmp
                  enddo ! j
               enddo ! i
            enddo ! ic
         endif
      enddo ! iatm
c
      do ipt=1,npts
         t1(ipt)=2*t1(ipt)
      enddo
      do i=1,nvec
         do ipt=1,npts
            tptmo(ipt,i)=t1(ipt)*amo_val(ipt,i)
         enddo
      enddo
      do iatm=1,ngridcent
         if (iatm.ne.latm) then
            do ic=1,3
               ipert=3*(iatm-1)+ic
               jpert=3*(latm-1)+ic
               do i=1,nvec 
                  do ipt=1,npts
                     gmorho(ipt,i)=dgrho(ipt,ipert,1)*amo_grad(ipt,i,1)
     &                            +dgrho(ipt,ipert,2)*amo_grad(ipt,i,2)
     &                            +dgrho(ipt,ipert,3)*amo_grad(ipt,i,3)
                  enddo ! ipt
               enddo ! i
               iki=0
               do i=1,nvec
                  iki=iki+i-1
                  do j=1,min(i,naocc)
                     tmp=0.0d0
                     do ipt=1,npts
                        tmp=tmp+gmorho(ipt,i)*tptmo(ipt,j)
     &                         +gmorho(ipt,j)*tptmo(ipt,i)
                     enddo ! ipt
                     qa(iki+j,ipert)=qa(iki+j,ipert)+tmp
                     qa(iki+j,jpert)=qa(iki+j,jpert)-tmp
                  enddo ! j
               enddo ! i
            enddo ! ic
         endif
      enddo ! iatm
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine rks_dksm_exp_mo_scr(npts,nao,nvec,naocc,npert,mxp,
     &           active_bfn_list,active_bfn_indx,active_bfn_atms,
     &           n_active_bfn, n_active_atm,
     &           wt,xc_vpt,xc_hpt,amo_val,bfng_val,
     &           avec,dmo,drho,tmo,qa)
      implicit none
c
c     Compute the explicit derivative contribution to the Kohn-Sham
c     matrix in MO-basis for the local density closed shell case.
c
c     The tactics are to construct to first and second derivatives
c     of the density with respect to the nuclear coordinates on the
c     fly first and subsequently consume them. Derivates of the MO's
c     are generated per orbital so as to minimise memory requirements.
c
c     Inputs
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
C *Parameter module for CCP1/DFT            
C *Memory information		
C *----------------
C * max_block 	- maximum number of blocks which can be allocate
      integer max_block
      parameter(max_block=20)
C *							
C *Basis sets information                              
C *----------------------				
C *max_tag 	- maximum number of basis sets which can be alloc
C *max_atype	- maximum number of centre types	
C *max_gtype	- maximum number of grid types (at least every element 
C                 has different grid type)
C *max_grids    - maximum number of grids (different terms may have
C                 different grids, i.e. CPKS equations use a different
C                 grid than the one used for the KS-matrix)
C *max_shel	- maximum number of shells on centre
C *max_prm	- maximum number of primitives for any given centre
C *maxL		- maximum angular momentum allowed	
C *max_func	- maximum number of basis functions for any centre
      integer max_tag,max_atype,max_shel,max_prm,max_ang
      integer max_gtype,max_grids
      parameter(max_tag=3,max_atype=30,max_shel=500,max_prm=5000)
      parameter(max_ang=5,max_gtype=10,max_grids=2)
C *								
C *Geometry information					
C *--------------------					
C *max_atom	- maximum number of atoms in system
      integer max_atom
      parameter(max_atom=750)
C *
C *Accuracy information				
C ---------------------			
C *global_accuracy 	- global accuracy	
      real*8  global_accuracy
      parameter(global_accuracy=1.0d-14)

C *Grid information
c *----------------
c *maxradzn - The maximum number of radial zones. 
      integer maxgpt,maxrad,maxfpt,maxang,maxradzn,maxtablerows
      parameter(maxgpt=2900,maxrad=50,maxang=302,maxfpt=100)
      parameter(maxradzn=35,maxtablerows=7)
c
c     The parameters defined here label the 1st, 2nd and 3rd derivative
c     components of the basis functions.
c
c     2nd derivatives:
c
      integer hxx, hyy, hzz, hxy, hxz, hyz
c
      parameter (hxx = 1)
      parameter (hyy = 2)
      parameter (hzz = 3)
      parameter (hxy = 4)
      parameter (hxz = 5)
      parameter (hyz = 6)
c
      integer hcc(3,3)
      data hcc/hxx,hxy,hxz,
     &         hxy,hyy,hyz,
     &         hxz,hyz,hzz/
c
c     3rd derivatives:
c
      integer txxx, txxy, txxz, txyy, txyz
      integer txzz, tyyy, tyyz, tyzz, tzzz
c
      parameter (txxx =  1)
      parameter (txxy =  2)
      parameter (txxz =  3)
      parameter (txyy =  4)
      parameter (txyz =  5)
      parameter (txzz =  6)
      parameter (tyyy =  7)
      parameter (tyyz =  8)
      parameter (tyzz =  9)
      parameter (tzzz = 10)
c
      integer tc3(3,3,3)
      data tc3/txxx,txxy,txxz,
     &         txxy,txyy,txyz,
     &         txxz,txyz,txzz,
     &         txxy,txyy,txyz,
     &         txyy,tyyy,tyyz,
     &         txyz,tyyz,tyzz,
     &         txxz,txyz,txzz,
     &         txyz,tyyz,tyzz,
     &         txzz,tyzz,tzzz/
c
c     The parameters defined here label the 1st and 2nd derivative
c     components of the density functional. Because a functional
c     may be purely local or of a GGA type the components are grouped
c     into derivative with respect to the density only and others.
c
c     1st derivatives:
c
      integer ira, irb
      integer igaa, igab, igbb
      integer ita, itb
c
c     2nd derivatives:
c
      integer irara, irarb, irbrb
      integer iragaa, iragab, iragbb, irbgaa, irbgab, irbgbb
      integer igaagaa, igaagab, igaagbb, igabgab, igabgbb, igbbgbb
      integer irata, iratb, irbta, irbtb
      integer igaata, igaatb, igabta, igabtb, igbbta, igbbtb
      integer itata, itatb, itbtb
c
      parameter (ira     = 1)
      parameter (irb     = 2)
c
      parameter (igaa    = 1)
      parameter (igab    = 2)
      parameter (igbb    = 3)
c
      parameter (ita     = 1)
      parameter (itb     = 2)
c
c
      parameter (irara   = 1)
      parameter (irarb   = 2)
      parameter (irbrb   = 3)
c
      parameter (iragaa  = 1)
      parameter (iragab  = 3)
      parameter (iragbb  = 4)
      parameter (irbgaa  = 5)
      parameter (irbgab  = 6)
      parameter (irbgbb  = 7)
      parameter (igaagaa = 2)
      parameter (igaagab = 8)
      parameter (igaagbb = 9)
      parameter (igabgab = 10)
      parameter (igabgbb = 11)
      parameter (igbbgbb = 12)
c
      parameter (irata   = 1)
      parameter (iratb   = 3)
      parameter (irbta   = 4)
      parameter (irbtb   = 5)
      parameter (itata   = 2)
      parameter (itatb   = 6)
      parameter (itbtb   = 7)
c
      parameter (igaata  = 1)
      parameter (igaatb  = 2)
      parameter (igabta  = 3)
      parameter (igabtb  = 4)
      parameter (igbbta  = 5)
      parameter (igbbtb  = 6)
c
      integer npts
      integer nao
      integer nvec
      integer naocc 
      integer npert
      integer mxp 
      integer n_active_bfn, n_active_atm
      integer active_bfn_list(n_active_bfn)
      integer active_bfn_indx(n_active_bfn)
      integer active_bfn_atms(n_active_atm)
      real*8 wt(npts)
      real*8 xc_vpt(mxp,2)
      real*8 xc_hpt(mxp,3)
      real*8 amo_val(npts,nvec)
      real*8 bfng_val(mxp,nao,3)
c     real*8 bfn_hess(mxp,nao,6)
      real*8 avec(nao,nvec)
c
c     Outputs
c
      real*8 qa(nvec*(nvec+1)/2,npert)
c
c     Workspace
c
      real*8 dmo(npts,npert)
c     real*8 ddmo(npts,6,ngridcent)
      real*8 drho(npts,npert)
      real*8 tmo(npts,nvec)
c
c     Local variables
c
      integer i, j
      integer n, n1
      integer ipert, iprt
      integer ipt
      integer iatm
      integer ic
      integer iki, ikj
      integer ao_tag
      parameter(ao_tag=1)
c
c     Code
c
      call aclear_dp(drho,npts*npert,0.0d0)
      do ipt=1,npts
         dmo(ipt,1)=wt(ipt)*xc_vpt(ipt,ira)
      enddo
      do i=1,nvec
         do ipt=1,npts
            tmo(ipt,i)=dmo(ipt,1)*amo_val(ipt,i)
         enddo
      enddo
c
c     Construct the derivatives of the density (times 0.5)
c
      iki=0
      do i=1,naocc
         call dmo_val_scr(npts,nao,nvec,npert,mxp,i,
     &                active_bfn_list,active_bfn_indx,n_active_bfn,
     &                bfng_val,avec,dmo)
         do ipert=1,3*n_active_atm
            do ipt=1,npts
               drho(ipt,ipert)=drho(ipt,ipert)
     &                        +amo_val(ipt,i)*dmo(ipt,ipert)
            enddo
         enddo
c
         iki=iki+i-1
         do iatm=1,n_active_atm
            do ic=1,3
               ipert=3*(active_bfn_atms(iatm)-1)+ic
               iprt =3*(iatm-1)+ic
               do j=1,i-1
                  do ipt=1,npts
                     qa(iki+j,ipert)=qa(iki+j,ipert)
     &               +dmo(ipt,iprt)*tmo(ipt,j)
                  enddo ! ipt
               enddo ! j
               do ipt=1,npts
                  qa(iki+i,ipert)=qa(iki+i,ipert)
     &            +2*dmo(ipt,iprt)*tmo(ipt,i)
               enddo ! ipt
               ikj=iki
               do j=i+1,nvec
                  ikj=ikj+j-1
                  do ipt=1,npts
                     qa(ikj+i,ipert)=qa(ikj+i,ipert)
     &               +dmo(ipt,iprt)*tmo(ipt,j)
                  enddo ! ipt
               enddo ! j
            enddo ! ic
         enddo ! iatm
      enddo
c
      n=(naocc-1)*naocc/2
      do i=naocc+1,nvec
         call dmo_val_scr(npts,nao,nvec,npert,mxp,i,
     &                active_bfn_list,active_bfn_indx,n_active_bfn,
     &                bfng_val,avec,dmo)
         n=n+i-1
         do iatm=1,n_active_atm
            do ic=1,3
               ipert=3*(active_bfn_atms(iatm)-1)+ic
               iprt =3*(iatm-1)+ic
               do j=1,naocc
                  do ipt=1,npts
                     qa(n+j,ipert)=qa(n+j,ipert)
     &               +dmo(ipt,iprt)*tmo(ipt,j)
                  enddo ! ipt
               enddo ! j
            enddo ! ic
         enddo ! iatm
      enddo ! i
c
c     Add the contributions to the Kohn-Sham matrix
c
      do ipt=1,npts 
         dmo(ipt,1)=2*wt(ipt)*xc_hpt(ipt,irara)
      enddo
      do iatm=1,n_active_atm
         do ic=1,3
            ipert=3*(active_bfn_atms(iatm)-1)+ic
            iprt =3*(iatm-1)+ic
            do ipt=1,npts 
               dmo(ipt,2)=dmo(ipt,1)*drho(ipt,iprt)
            enddo
            n=0
            do j=1,naocc
               n1=n
               n=n+j-1
               do ipt=1,npts
                  dmo(ipt,3)=dmo(ipt,2)*amo_val(ipt,j)
               enddo
               do i=j,nvec
                  n1=n1+i-1
                  do ipt=1,npts
                     qa(n1+j,ipert)=qa(n1+j,ipert)
     &               +dmo(ipt,3)*amo_val(ipt,i)
                  enddo ! ipt
               enddo ! i
            enddo ! j
         enddo ! ic
      enddo ! iatm
c 
      end
c
c-----------------------------------------------------------------------
c
      subroutine rks_dksm_exp_mo_gga_scr(npts,nao,nvec,naocc,npert,
     &           mxp,active_bfn_list,active_bfn_indx,
     &           active_bfn_atms,n_active_bfn, n_active_atm,
     &           wt,xc_vpt,xc_dvpt,xc_hpt,xc_dhpt,amo_val,amo_grad,
     &           bfng_val,bfn_hess,avec,grho,dmo,dgmo,drho,dgrho,gmorho,
     &           grdgmo,t1,t2,t3,t4,tptmo,qa)
      implicit none
c
c     Compute the explicit derivative contribution to the Kohn-Sham 
c     matrix in MO-basis for the gradient corrected closed shell case.
c
c     The tactics are to construct to first and second derivatives
c     of the density with respect to the nuclear coordinates on the
c     fly first and subsequently consume them. Derivates of the MO's
c     are generated per orbital so as to minimise memory requirements.
c
c     Inputs
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
C *Parameter module for CCP1/DFT            
C *Memory information		
C *----------------
C * max_block 	- maximum number of blocks which can be allocate
      integer max_block
      parameter(max_block=20)
C *							
C *Basis sets information                              
C *----------------------				
C *max_tag 	- maximum number of basis sets which can be alloc
C *max_atype	- maximum number of centre types	
C *max_gtype	- maximum number of grid types (at least every element 
C                 has different grid type)
C *max_grids    - maximum number of grids (different terms may have
C                 different grids, i.e. CPKS equations use a different
C                 grid than the one used for the KS-matrix)
C *max_shel	- maximum number of shells on centre
C *max_prm	- maximum number of primitives for any given centre
C *maxL		- maximum angular momentum allowed	
C *max_func	- maximum number of basis functions for any centre
      integer max_tag,max_atype,max_shel,max_prm,max_ang
      integer max_gtype,max_grids
      parameter(max_tag=3,max_atype=30,max_shel=500,max_prm=5000)
      parameter(max_ang=5,max_gtype=10,max_grids=2)
C *								
C *Geometry information					
C *--------------------					
C *max_atom	- maximum number of atoms in system
      integer max_atom
      parameter(max_atom=750)
C *
C *Accuracy information				
C ---------------------			
C *global_accuracy 	- global accuracy	
      real*8  global_accuracy
      parameter(global_accuracy=1.0d-14)

C *Grid information
c *----------------
c *maxradzn - The maximum number of radial zones. 
      integer maxgpt,maxrad,maxfpt,maxang,maxradzn,maxtablerows
      parameter(maxgpt=2900,maxrad=50,maxang=302,maxfpt=100)
      parameter(maxradzn=35,maxtablerows=7)
c
c     The parameters defined here label the 1st, 2nd and 3rd derivative
c     components of the basis functions.
c
c     2nd derivatives:
c
      integer hxx, hyy, hzz, hxy, hxz, hyz
c
      parameter (hxx = 1)
      parameter (hyy = 2)
      parameter (hzz = 3)
      parameter (hxy = 4)
      parameter (hxz = 5)
      parameter (hyz = 6)
c
      integer hcc(3,3)
      data hcc/hxx,hxy,hxz,
     &         hxy,hyy,hyz,
     &         hxz,hyz,hzz/
c
c     3rd derivatives:
c
      integer txxx, txxy, txxz, txyy, txyz
      integer txzz, tyyy, tyyz, tyzz, tzzz
c
      parameter (txxx =  1)
      parameter (txxy =  2)
      parameter (txxz =  3)
      parameter (txyy =  4)
      parameter (txyz =  5)
      parameter (txzz =  6)
      parameter (tyyy =  7)
      parameter (tyyz =  8)
      parameter (tyzz =  9)
      parameter (tzzz = 10)
c
      integer tc3(3,3,3)
      data tc3/txxx,txxy,txxz,
     &         txxy,txyy,txyz,
     &         txxz,txyz,txzz,
     &         txxy,txyy,txyz,
     &         txyy,tyyy,tyyz,
     &         txyz,tyyz,tyzz,
     &         txxz,txyz,txzz,
     &         txyz,tyyz,tyzz,
     &         txzz,tyzz,tzzz/
c
c     The parameters defined here label the 1st and 2nd derivative
c     components of the density functional. Because a functional
c     may be purely local or of a GGA type the components are grouped
c     into derivative with respect to the density only and others.
c
c     1st derivatives:
c
      integer ira, irb
      integer igaa, igab, igbb
      integer ita, itb
c
c     2nd derivatives:
c
      integer irara, irarb, irbrb
      integer iragaa, iragab, iragbb, irbgaa, irbgab, irbgbb
      integer igaagaa, igaagab, igaagbb, igabgab, igabgbb, igbbgbb
      integer irata, iratb, irbta, irbtb
      integer igaata, igaatb, igabta, igabtb, igbbta, igbbtb
      integer itata, itatb, itbtb
c
      parameter (ira     = 1)
      parameter (irb     = 2)
c
      parameter (igaa    = 1)
      parameter (igab    = 2)
      parameter (igbb    = 3)
c
      parameter (ita     = 1)
      parameter (itb     = 2)
c
c
      parameter (irara   = 1)
      parameter (irarb   = 2)
      parameter (irbrb   = 3)
c
      parameter (iragaa  = 1)
      parameter (iragab  = 3)
      parameter (iragbb  = 4)
      parameter (irbgaa  = 5)
      parameter (irbgab  = 6)
      parameter (irbgbb  = 7)
      parameter (igaagaa = 2)
      parameter (igaagab = 8)
      parameter (igaagbb = 9)
      parameter (igabgab = 10)
      parameter (igabgbb = 11)
      parameter (igbbgbb = 12)
c
      parameter (irata   = 1)
      parameter (iratb   = 3)
      parameter (irbta   = 4)
      parameter (irbtb   = 5)
      parameter (itata   = 2)
      parameter (itatb   = 6)
      parameter (itbtb   = 7)
c
      parameter (igaata  = 1)
      parameter (igaatb  = 2)
      parameter (igabta  = 3)
      parameter (igabtb  = 4)
      parameter (igbbta  = 5)
      parameter (igbbtb  = 6)
c
      integer npts
      integer nao
      integer nvec
      integer naocc 
      integer npert
      integer mxp 
      integer n_active_bfn, n_active_atm
      integer active_bfn_list(n_active_bfn)
      integer active_bfn_indx(n_active_bfn)
      integer active_bfn_atms(n_active_atm)
      real*8 wt(npts)
      real*8 xc_vpt(mxp,2)
      real*8 xc_dvpt(mxp,3)
      real*8 xc_hpt(mxp,3)
      real*8 xc_dhpt(mxp,12)
      real*8 amo_val(npts,nvec)
      real*8 amo_grad(npts,nvec,3)
      real*8 bfng_val(mxp,nao,3)
      real*8 bfn_hess(mxp,nao,6)
c     real*8 bfn_3rd(mxp,nao,10)
      real*8 avec(nao,nvec)
      real*8 grho(mxp,2,3)
c
c     Outputs
c
      real*8 qa(nao*(nao+1)/2,npert)
c
c     Workspace
c
      real*8 dmo(npts,npert)
      real*8 dgmo(npts,npert,3)
c     real*8 ddmo(npts,6,ngridcent)
c     real*8 ddgmo(npts,10,ngridcent)
      real*8 drho(npts,npert)
c     real*8 ddrho(npts,npert*(npert+1)/2)
      real*8 dgrho(npts,npert,3)
c     real*8 ddgrho(npts,npert*(npert+1)/2,3)
      real*8 gmorho(npts,nvec)
      real*8 grdgmo(npts)
      real*8 t1(npts)
      real*8 t2(npts)
      real*8 t3(npts)
      real*8 t4(npts)
      real*8 tptmo(npts,nvec)
c
c     Local variables
c
      integer i, j
c     integer n
      integer ipert, iprt
      integer ipt
      integer iatm
      integer ic
      integer iki, ikj
      real*8 tmp
      integer ao_tag
      parameter(ao_tag=1)
c
c     Code
c
      call aclear_dp(drho,npts*npert,0.0d0)
      call aclear_dp(dgrho,npts*npert*3,0.0d0)
c
c     Construct dot product of grad rho and grad mo
c
      do i=1,nvec
         do ipt=1,npts
            gmorho(ipt,i)=
     &      (grho(ipt,1,1)*amo_grad(ipt,i,1)
     &      +grho(ipt,1,2)*amo_grad(ipt,i,2)
     &      +grho(ipt,1,3)*amo_grad(ipt,i,3))
            gmorho(ipt,i)=0.5d0*gmorho(ipt,i)
         enddo ! ipt
      enddo ! i
      do ipt=1,npts
         t1(ipt)=wt(ipt)*xc_dvpt(ipt,igaa)
         t2(ipt)=wt(ipt)*xc_vpt(ipt,ira)
      enddo
      do i=1,nvec
         do ipt=1,npts
            tptmo(ipt,i)=t2(ipt)*amo_val(ipt,i)
     &                  +t1(ipt)*gmorho(ipt,i)
         enddo
      enddo
c
c     Construct the derivatives of the density (times 0.5)
c
      iki=0
      do i=1,naocc
         call dmo_val_scr(npts,nao,nvec,npert,mxp,i,
     &                active_bfn_list,active_bfn_indx,n_active_bfn,
     &                bfng_val,avec,dmo)
         do ipert=1,3*n_active_atm
            do ipt=1,npts
               drho(ipt,ipert)=drho(ipt,ipert)
     &                        +amo_val(ipt,i)*dmo(ipt,ipert)
            enddo ! ipt
         enddo ! ipert
         do ic=1,3
            do ipert=1,3*n_active_atm
               do ipt=1,npts
                  dgrho(ipt,ipert,ic)=dgrho(ipt,ipert,ic)
     &            +amo_grad(ipt,i,ic)*dmo(ipt,ipert)
               enddo ! ipt
            enddo ! ipert
         enddo ! ic

         iki=iki+i-1
         do iatm=1,n_active_atm
            do ic=1,3
               ipert=3*(active_bfn_atms(iatm)-1)+ic
               iprt =3*(iatm-1)+ic
               do j=1,i-1
                  do ipt=1,npts
                     qa(iki+j,ipert)=qa(iki+j,ipert)
     &               +dmo(ipt,iprt)*tptmo(ipt,j)
                  enddo ! ipt
               enddo ! j
               tmp=0.0d0
               do ipt=1,npts
                  tmp=tmp
     &            +dmo(ipt,iprt)*tptmo(ipt,i)
               enddo ! ipt
               qa(iki+i,ipert)=qa(iki+i,ipert)+2*tmp
               ikj=iki
               do j=i+1,nvec
                  ikj=ikj+j-1
                  do ipt=1,npts
                     qa(ikj+i,ipert)=qa(ikj+i,ipert)
     &               +dmo(ipt,iprt)*tptmo(ipt,j)
                  enddo ! ipt
               enddo ! j
            enddo ! ic
         enddo ! iatm
      enddo ! i
c
      iki=naocc*(naocc-1)/2
      do i=naocc+1,nvec
         call dmo_val_scr(npts,nao,nvec,npert,mxp,i,
     &                active_bfn_list,active_bfn_indx,n_active_bfn,
     &                bfng_val,avec,dmo)
         iki=iki+i-1
         do iatm=1,n_active_atm
            do ic=1,3
               ipert=3*(active_bfn_atms(iatm)-1)+ic
               iprt =3*(iatm-1)+ic
               do j=1,naocc
                  do ipt=1,npts
                     qa(iki+j,ipert)=qa(iki+j,ipert)
     &               +dmo(ipt,iprt)*tptmo(ipt,j)
                  enddo ! ipt
               enddo ! ic
            enddo ! j
         enddo ! iatm
      enddo ! i
c
      do ipt=1,npts
         t2(ipt)=4*wt(ipt)*xc_dhpt(ipt,iragaa)
         t3(ipt)=2*wt(ipt)*xc_hpt(ipt,irara)
      enddo
      do i=1,nvec
         do ipt=1,npts
            tptmo(ipt,i)=t3(ipt)*amo_val(ipt,i)
     &                  +t2(ipt)*gmorho(ipt,i)
         enddo
      enddo
      do iatm=1,n_active_atm
         do ic=1,3
            ipert=3*(active_bfn_atms(iatm)-1)+ic
            iprt =3*(iatm-1)+ic
            iki=0
            do i=1,nvec
               iki=iki+i-1
               do j=1,min(i,naocc)
                  do ipt=1,npts
                     qa(iki+j,ipert)=qa(iki+j,ipert)
     &               +tptmo(ipt,i)*amo_val(ipt,j)*drho(ipt,iprt)
                  enddo ! ipt
               enddo ! j
            enddo ! i
         enddo ! ic
      enddo ! iatm
      do i=1,nvec
         do ipt=1,npts
            tptmo(ipt,i)=t2(ipt)*amo_val(ipt,i)
         enddo
      enddo
      do iatm=1,n_active_atm
         do ic=1,3
            ipert=3*(active_bfn_atms(iatm)-1)+ic
            iprt =3*(iatm-1)+ic
            iki=0
            do i=1,nvec
               iki=iki+i-1
               do j=1,min(i,naocc)
                  do ipt=1,npts
                     qa(iki+j,ipert)=qa(iki+j,ipert)
     &               +tptmo(ipt,i)*gmorho(ipt,j)*drho(ipt,iprt)
                  enddo ! ipt
               enddo ! j
            enddo ! i
         enddo ! ic
      enddo ! iatm
c
      do i=1,nvec
         do ipt=1,npts
            tptmo(ipt,i)=t1(ipt)*amo_val(ipt,i)
         enddo
      enddo
      iki=0
      do i=1,naocc
         call dgmo_val_scr(npts,nao,nvec,npert,mxp,i,
     &                 active_bfn_list,active_bfn_indx,n_active_bfn,
     &                 bfn_hess,avec,dgmo)
         do ic=1,3
            do ipert=1,3*n_active_atm
               do ipt=1,npts
                  dgrho(ipt,ipert,ic)=dgrho(ipt,ipert,ic)
     &            +amo_val(ipt,i)*dgmo(ipt,ipert,ic)
               enddo ! ipt
            enddo ! ipert
         enddo ! ic
c
         iki=iki+i-1
         do iatm=1,n_active_atm
            do ic=1,3
               ipert=3*(active_bfn_atms(iatm)-1)+ic
               iprt =3*(iatm-1)+ic
               do ipt=1,npts
                  grdgmo(ipt) = grho(ipt,1,1)*dgmo(ipt,iprt,1)
     &                        + grho(ipt,1,2)*dgmo(ipt,iprt,2)
     &                        + grho(ipt,1,3)*dgmo(ipt,iprt,3)
                  grdgmo(ipt) = 0.5d0*grdgmo(ipt)
               enddo
               do j=1,i-1
                  do ipt=1,npts
                     qa(iki+j,ipert)=qa(iki+j,ipert)
     &               +tptmo(ipt,j)*grdgmo(ipt)
                  enddo ! ipt
               enddo ! j
               tmp=0.0d0
               do ipt=1,npts
                  tmp=tmp
     &            +tptmo(ipt,i)*grdgmo(ipt)
               enddo ! ipt
               qa(iki+i,ipert)=qa(iki+i,ipert)+2*tmp
               ikj=iki
               do j=i+1,nvec
                  ikj=ikj+j-1
                  do ipt=1,npts
                     qa(ikj+i,ipert)=qa(ikj+i,ipert)
     &               +tptmo(ipt,j)*grdgmo(ipt)
                  enddo ! ipt
               enddo ! j
            enddo ! ic
         enddo ! iatm
      enddo ! i
c
      iki=naocc*(naocc-1)/2
      do i=naocc+1,nvec
         call dgmo_val_scr(npts,nao,nvec,npert,mxp,i,
     &                 active_bfn_list,active_bfn_indx,n_active_bfn,
     &                 bfn_hess,avec,dgmo)
         iki=iki+i-1
         do iatm=1,n_active_atm
            do ic=1,3
               ipert=3*(active_bfn_atms(iatm)-1)+ic
               iprt =3*(iatm-1)+ic
               do ipt=1,npts
                  grdgmo(ipt) = grho(ipt,1,1)*dgmo(ipt,iprt,1)
     &                        + grho(ipt,1,2)*dgmo(ipt,iprt,2)
     &                        + grho(ipt,1,3)*dgmo(ipt,iprt,3)
                  grdgmo(ipt) = 0.5d0*grdgmo(ipt)
               enddo
               do j=1,naocc
                  do ipt=1,npts
                     qa(iki+j,ipert)=qa(iki+j,ipert)
     &               +tptmo(ipt,j)*grdgmo(ipt)
                  enddo ! ipt
               enddo ! j
            enddo ! ic
         enddo ! iatm
      enddo ! i
c
      do ipt=1,npts
         t3(ipt)=4*wt(ipt)*xc_dhpt(ipt,igaagaa)
      enddo
      do i=1,nvec
         do ipt=1,npts
            tptmo(ipt,i)=t2(ipt)*amo_val(ipt,i)
     &                  +t3(ipt)*gmorho(ipt,i)
         enddo
      enddo
      do iatm=1,n_active_atm
         do ic=1,3
            ipert=3*(active_bfn_atms(iatm)-1)+ic
            iprt =3*(iatm-1)+ic
            do ipt=1,npts
               grdgmo(ipt) = grho(ipt,1,1)*dgrho(ipt,iprt,1)
     &                     + grho(ipt,1,2)*dgrho(ipt,iprt,2)
     &                     + grho(ipt,1,3)*dgrho(ipt,iprt,3)
               grdgmo(ipt) = 0.5d0*grdgmo(ipt)
            enddo ! ipt
            iki=0
            do i=1,nvec
               iki=iki+i-1
               do j=1,min(i,naocc)
                  do ipt=1,npts
                     qa(iki+j,ipert)=qa(iki+j,ipert)
     &               +tptmo(ipt,i)*amo_val(ipt,j)*grdgmo(ipt)
                  enddo ! ipt
               enddo ! j
            enddo ! i
         enddo ! ic
      enddo ! iatm
      do i=1,nvec
         do ipt=1,npts
            tptmo(ipt,i)=t3(ipt)*amo_val(ipt,i)
         enddo
      enddo
      do iatm=1,n_active_atm
         do ic=1,3
            ipert=3*(active_bfn_atms(iatm)-1)+ic
            iprt =3*(iatm-1)+ic
            do ipt=1,npts
               grdgmo(ipt) = grho(ipt,1,1)*dgrho(ipt,iprt,1)
     &                     + grho(ipt,1,2)*dgrho(ipt,iprt,2)
     &                     + grho(ipt,1,3)*dgrho(ipt,iprt,3)
               grdgmo(ipt) = 0.5d0*grdgmo(ipt)
            enddo ! ipt
            iki=0
            do i=1,nvec
               iki=iki+i-1
               do j=1,min(i,naocc)
                  do ipt=1,npts
                     qa(iki+j,ipert)=qa(iki+j,ipert)
     &               +tptmo(ipt,i)*gmorho(ipt,j)*grdgmo(ipt)
                  enddo ! ipt
               enddo ! j
            enddo ! i
         enddo ! ic
      enddo ! iatm
c
      do ipt=1,npts
         t1(ipt)=2*t1(ipt)
      enddo
      do i=1,nvec
         do ipt=1,npts
            tptmo(ipt,i)=t1(ipt)*amo_val(ipt,i)
         enddo
      enddo
      do iatm=1,n_active_atm
         do ic=1,3
            ipert=3*(active_bfn_atms(iatm)-1)+ic
            iprt =3*(iatm-1)+ic
            do i=1,nvec 
               do ipt=1,npts
                  gmorho(ipt,i) = dgrho(ipt,iprt,1)*amo_grad(ipt,i,1)
     &                          + dgrho(ipt,iprt,2)*amo_grad(ipt,i,2)
     &                          + dgrho(ipt,iprt,3)*amo_grad(ipt,i,3)
               enddo ! ipt
            enddo ! i
            iki=0
            do i=1,nvec
               iki=iki+i-1
               do j=1,min(i,naocc)
                  do ipt=1,npts
                     qa(iki+j,ipert)=qa(iki+j,ipert)
     &                +gmorho(ipt,i)*tptmo(ipt,j)
     &                +gmorho(ipt,j)*tptmo(ipt,i)
                  enddo ! ipt
               enddo ! j
            enddo ! i
         enddo ! ic
      enddo ! iatm
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine rks_dksm_exp_mo_gwt_scr(npts,nao,nvec,naocc,npert,mxp,
     &           active_bfn_list,active_bfn_indx,active_bfn_atms,
     &           n_active_bfn, n_active_atm,
     &           num_near,near_atom,
     &           gwt_avail_sw,latm,
     &           wt,gwt,xc_vpt,xc_hpt,amo_val,bfng_val,
     &           avec,dmo,drho,tmo,qa)
      implicit none
c
c     Compute the explicit derivative contribution to the Kohn-Sham
c     matrix in MO-basis for the local density closed shell case.
c
c     The tactics are to construct to first and second derivatives
c     of the density with respect to the nuclear coordinates on the
c     fly first and subsequently consume them. Derivates of the MO's
c     are generated per orbital so as to minimise memory requirements.
c
c     Inputs
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
C *Parameter module for CCP1/DFT            
C *Memory information		
C *----------------
C * max_block 	- maximum number of blocks which can be allocate
      integer max_block
      parameter(max_block=20)
C *							
C *Basis sets information                              
C *----------------------				
C *max_tag 	- maximum number of basis sets which can be alloc
C *max_atype	- maximum number of centre types	
C *max_gtype	- maximum number of grid types (at least every element 
C                 has different grid type)
C *max_grids    - maximum number of grids (different terms may have
C                 different grids, i.e. CPKS equations use a different
C                 grid than the one used for the KS-matrix)
C *max_shel	- maximum number of shells on centre
C *max_prm	- maximum number of primitives for any given centre
C *maxL		- maximum angular momentum allowed	
C *max_func	- maximum number of basis functions for any centre
      integer max_tag,max_atype,max_shel,max_prm,max_ang
      integer max_gtype,max_grids
      parameter(max_tag=3,max_atype=30,max_shel=500,max_prm=5000)
      parameter(max_ang=5,max_gtype=10,max_grids=2)
C *								
C *Geometry information					
C *--------------------					
C *max_atom	- maximum number of atoms in system
      integer max_atom
      parameter(max_atom=750)
C *
C *Accuracy information				
C ---------------------			
C *global_accuracy 	- global accuracy	
      real*8  global_accuracy
      parameter(global_accuracy=1.0d-14)

C *Grid information
c *----------------
c *maxradzn - The maximum number of radial zones. 
      integer maxgpt,maxrad,maxfpt,maxang,maxradzn,maxtablerows
      parameter(maxgpt=2900,maxrad=50,maxang=302,maxfpt=100)
      parameter(maxradzn=35,maxtablerows=7)
c
c     The parameters defined here label the 1st, 2nd and 3rd derivative
c     components of the basis functions.
c
c     2nd derivatives:
c
      integer hxx, hyy, hzz, hxy, hxz, hyz
c
      parameter (hxx = 1)
      parameter (hyy = 2)
      parameter (hzz = 3)
      parameter (hxy = 4)
      parameter (hxz = 5)
      parameter (hyz = 6)
c
      integer hcc(3,3)
      data hcc/hxx,hxy,hxz,
     &         hxy,hyy,hyz,
     &         hxz,hyz,hzz/
c
c     3rd derivatives:
c
      integer txxx, txxy, txxz, txyy, txyz
      integer txzz, tyyy, tyyz, tyzz, tzzz
c
      parameter (txxx =  1)
      parameter (txxy =  2)
      parameter (txxz =  3)
      parameter (txyy =  4)
      parameter (txyz =  5)
      parameter (txzz =  6)
      parameter (tyyy =  7)
      parameter (tyyz =  8)
      parameter (tyzz =  9)
      parameter (tzzz = 10)
c
      integer tc3(3,3,3)
      data tc3/txxx,txxy,txxz,
     &         txxy,txyy,txyz,
     &         txxz,txyz,txzz,
     &         txxy,txyy,txyz,
     &         txyy,tyyy,tyyz,
     &         txyz,tyyz,tyzz,
     &         txxz,txyz,txzz,
     &         txyz,tyyz,tyzz,
     &         txzz,tyzz,tzzz/
c
c     The parameters defined here label the 1st and 2nd derivative
c     components of the density functional. Because a functional
c     may be purely local or of a GGA type the components are grouped
c     into derivative with respect to the density only and others.
c
c     1st derivatives:
c
      integer ira, irb
      integer igaa, igab, igbb
      integer ita, itb
c
c     2nd derivatives:
c
      integer irara, irarb, irbrb
      integer iragaa, iragab, iragbb, irbgaa, irbgab, irbgbb
      integer igaagaa, igaagab, igaagbb, igabgab, igabgbb, igbbgbb
      integer irata, iratb, irbta, irbtb
      integer igaata, igaatb, igabta, igabtb, igbbta, igbbtb
      integer itata, itatb, itbtb
c
      parameter (ira     = 1)
      parameter (irb     = 2)
c
      parameter (igaa    = 1)
      parameter (igab    = 2)
      parameter (igbb    = 3)
c
      parameter (ita     = 1)
      parameter (itb     = 2)
c
c
      parameter (irara   = 1)
      parameter (irarb   = 2)
      parameter (irbrb   = 3)
c
      parameter (iragaa  = 1)
      parameter (iragab  = 3)
      parameter (iragbb  = 4)
      parameter (irbgaa  = 5)
      parameter (irbgab  = 6)
      parameter (irbgbb  = 7)
      parameter (igaagaa = 2)
      parameter (igaagab = 8)
      parameter (igaagbb = 9)
      parameter (igabgab = 10)
      parameter (igabgbb = 11)
      parameter (igbbgbb = 12)
c
      parameter (irata   = 1)
      parameter (iratb   = 3)
      parameter (irbta   = 4)
      parameter (irbtb   = 5)
      parameter (itata   = 2)
      parameter (itatb   = 6)
      parameter (itbtb   = 7)
c
      parameter (igaata  = 1)
      parameter (igaatb  = 2)
      parameter (igabta  = 3)
      parameter (igabtb  = 4)
      parameter (igbbta  = 5)
      parameter (igbbtb  = 6)
c
      integer npts
      integer nao
      integer nvec
      integer naocc 
      integer npert
      integer mxp 
      integer n_active_bfn, n_active_atm
      integer active_bfn_list(n_active_bfn)
      integer active_bfn_indx(n_active_bfn)
      integer active_bfn_atms(n_active_atm)
      integer num_near, near_atom(num_near)
      logical gwt_avail_sw
      integer latm
      real*8 wt(npts)
      real*8 gwt(3,mxp,num_near)
      real*8 xc_vpt(mxp,2)
      real*8 xc_hpt(mxp,3)
      real*8 amo_val(npts,nvec)
      real*8 bfng_val(mxp,nao,3)
c     real*8 bfn_hess(mxp,nao,6)
      real*8 avec(nao,nvec)
c
c     Outputs
c
      real*8 qa(nvec*(nvec+1)/2,npert)
c
c     Workspace
c
      real*8 dmo(npts,npert)
c     real*8 ddmo(npts,6,ngridcent)
      real*8 drho(npts,npert)
      real*8 tmo(npts,nvec)
c     real*8 ddrho(npts,npert*(npert+1)/2)
c
c     Local variables
c
      integer i, j
      integer n, n1
      integer ipert, jpert, iprt
      integer ipt
      integer iatm
      integer ic
      integer iki, ikj
      integer iab, lab
      real*8 tmp
      integer ao_tag
      parameter(ao_tag=1)
c
c     Code:
c
c     First the contributions of the derivatives of the weights
c
      if (gwt_avail_sw) then
         do iatm=1,num_near
            do ic=1,3
               ipert=3*(near_atom(iatm)-1)+ic
               do j=1,nvec
                  ikj=iky(j)
                  do i=1,min(j,naocc)
                     do ipt=1,npts
                        qa(ikj+i,ipert)=qa(ikj+i,ipert)
     &                  +gwt(ic,ipt,iatm)*xc_vpt(ipt,ira)
     &                  *amo_val(ipt,i)*amo_val(ipt,j)
                     enddo ! ipt
                  enddo ! i
               enddo ! j
            enddo ! ic
         enddo ! iatm
      endif
c
c     Now the derivatives of the Fock operator
c
      call aclear_dp(drho,npts*npert,0.0d0)
      do ipt=1,npts
         dmo(ipt,1)=wt(ipt)*xc_vpt(ipt,ira)
      enddo
      do i=1,nvec
         do ipt=1,npts
            tmo(ipt,i)=dmo(ipt,1)*amo_val(ipt,i)
         enddo
      enddo
c
c     Construct the derivatives of the density (times 0.5)
c
      iki=0
      lab=3*(latm-1)
      do i=1,naocc
         call dmo_val_scr(npts,nao,nvec,npert,mxp,i,
     &                active_bfn_list,active_bfn_indx,n_active_bfn,
     &                bfng_val,avec,dmo)
         do ipert=1,3*n_active_atm
            do ipt=1,npts
               drho(ipt,ipert)=drho(ipt,ipert)
     &                        +amo_val(ipt,i)*dmo(ipt,ipert)
            enddo
         enddo
c
         iki=iki+i-1
         do iatm=1,n_active_atm
            iab=3*(active_bfn_atms(iatm)-1)
            if (iab.ne.lab) then
               do ic=1,3
                  ipert=iab+ic
                  jpert=lab+ic
                  iprt =3*(iatm-1)+ic
                  do j=1,i-1
                     tmp=0.0d0
                     do ipt=1,npts
                        tmp=tmp+dmo(ipt,iprt)*tmo(ipt,j)
                     enddo ! ipt
                     qa(iki+j,ipert)=qa(iki+j,ipert)+tmp
                     qa(iki+j,jpert)=qa(iki+j,jpert)-tmp
                  enddo ! j
                  tmp=0.0d0
                  do ipt=1,npts
                     tmp=tmp+dmo(ipt,iprt)*tmo(ipt,i)
                  enddo ! ipt
                  tmp=2.0d0*tmp
                  qa(iki+i,ipert)=qa(iki+i,ipert)+tmp
                  qa(iki+i,jpert)=qa(iki+i,jpert)-tmp
                  ikj=iki
                  do j=i+1,nvec
                     ikj=ikj+j-1
                     tmp=0.0d0
                     do ipt=1,npts
                        tmp=tmp+dmo(ipt,iprt)*tmo(ipt,j)
                     enddo
                     qa(ikj+i,ipert)=qa(ikj+i,ipert)+tmp
                     qa(ikj+i,jpert)=qa(ikj+i,jpert)-tmp
                  enddo ! j
               enddo ! ic
            endif
         enddo ! iatm
      enddo
c
      n=(naocc-1)*naocc/2
      do i=naocc+1,nvec
         call dmo_val_scr(npts,nao,nvec,npert,mxp,i,
     &                active_bfn_list,active_bfn_indx,n_active_bfn,
     &                bfng_val,avec,dmo)
         n=n+i-1
         do iatm=1,n_active_atm
            iab=3*(active_bfn_atms(iatm)-1)
            if (iab.ne.lab) then
               do ic=1,3
                  ipert=iab+ic
                  jpert=lab+ic
                  iprt =3*(iatm-1)+ic
                  do j=1,naocc
                     tmp=0.0d0
                     do ipt=1,npts
                        tmp=tmp+dmo(ipt,iprt)*tmo(ipt,j)
                     enddo ! ipt
                     qa(n+j,ipert)=qa(n+j,ipert)+tmp
                     qa(n+j,jpert)=qa(n+j,jpert)-tmp
                  enddo ! j
               enddo ! ic
            endif 
         enddo ! iatm
      enddo
c
c     Add the contributions to the Kohn-Sham matrix
c
      do ipt=1,npts 
         dmo(ipt,1)=2*wt(ipt)*xc_hpt(ipt,irara)
      enddo
      do iatm=1,n_active_atm
         iab=3*(active_bfn_atms(iatm)-1)
         if (iab.ne.lab) then
            do ic=1,3
               ipert=iab+ic
               jpert=lab+ic
               iprt =3*(iatm-1)+ic
               do ipt=1,npts 
                  dmo(ipt,2)=dmo(ipt,1)*drho(ipt,iprt)
               enddo
               n=0
               do j=1,naocc
                  n1=n
                  n=n+j-1
                  do ipt=1,npts
                     dmo(ipt,3)=dmo(ipt,2)*amo_val(ipt,j)
                  enddo
                  do i=j,nvec
                     n1=n1+i-1
                     tmp=0.0d0
                     do ipt=1,npts
                        tmp=tmp+dmo(ipt,3)*amo_val(ipt,i)
                     enddo ! ipt
                     qa(n1+j,ipert)=qa(n1+j,ipert)+tmp
                     qa(n1+j,jpert)=qa(n1+j,jpert)-tmp
                  enddo ! i
               enddo ! j
            enddo ! ic
         endif
      enddo ! iatm
c 
      end
c
c-----------------------------------------------------------------------
c
      subroutine rks_dksm_exp_mo_gga_gwt_scr(npts,nao,nvec,naocc,npert,
     &           mxp,
     &           active_bfn_list,active_bfn_indx,active_bfn_atms,
     &           n_active_bfn, n_active_atm,
     &           num_near,near_atom,
     &           gwt_avail_sw,latm,
     &           wt,gwt,xc_vpt,xc_dvpt,xc_hpt,xc_dhpt,amo_val,amo_grad,
     &           bfng_val,bfn_hess,
     &           avec,grho,dmo,dgmo,drho,dgrho,gmorho,grdgmo,
     &           t1,t2,t3,t4,tptmo,qa)
      implicit none
c
c     Compute the explicit derivative contribution to the Kohn-Sham 
c     matrix in MO-basis for the gradient corrected closed shell case.
c
c     The tactics are to construct to first and second derivatives
c     of the density with respect to the nuclear coordinates on the
c     fly first and subsequently consume them. Derivates of the MO's
c     are generated per orbital so as to minimise memory requirements.
c
c     Inputs
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
C *Parameter module for CCP1/DFT            
C *Memory information		
C *----------------
C * max_block 	- maximum number of blocks which can be allocate
      integer max_block
      parameter(max_block=20)
C *							
C *Basis sets information                              
C *----------------------				
C *max_tag 	- maximum number of basis sets which can be alloc
C *max_atype	- maximum number of centre types	
C *max_gtype	- maximum number of grid types (at least every element 
C                 has different grid type)
C *max_grids    - maximum number of grids (different terms may have
C                 different grids, i.e. CPKS equations use a different
C                 grid than the one used for the KS-matrix)
C *max_shel	- maximum number of shells on centre
C *max_prm	- maximum number of primitives for any given centre
C *maxL		- maximum angular momentum allowed	
C *max_func	- maximum number of basis functions for any centre
      integer max_tag,max_atype,max_shel,max_prm,max_ang
      integer max_gtype,max_grids
      parameter(max_tag=3,max_atype=30,max_shel=500,max_prm=5000)
      parameter(max_ang=5,max_gtype=10,max_grids=2)
C *								
C *Geometry information					
C *--------------------					
C *max_atom	- maximum number of atoms in system
      integer max_atom
      parameter(max_atom=750)
C *
C *Accuracy information				
C ---------------------			
C *global_accuracy 	- global accuracy	
      real*8  global_accuracy
      parameter(global_accuracy=1.0d-14)

C *Grid information
c *----------------
c *maxradzn - The maximum number of radial zones. 
      integer maxgpt,maxrad,maxfpt,maxang,maxradzn,maxtablerows
      parameter(maxgpt=2900,maxrad=50,maxang=302,maxfpt=100)
      parameter(maxradzn=35,maxtablerows=7)
c
c     The parameters defined here label the 1st, 2nd and 3rd derivative
c     components of the basis functions.
c
c     2nd derivatives:
c
      integer hxx, hyy, hzz, hxy, hxz, hyz
c
      parameter (hxx = 1)
      parameter (hyy = 2)
      parameter (hzz = 3)
      parameter (hxy = 4)
      parameter (hxz = 5)
      parameter (hyz = 6)
c
      integer hcc(3,3)
      data hcc/hxx,hxy,hxz,
     &         hxy,hyy,hyz,
     &         hxz,hyz,hzz/
c
c     3rd derivatives:
c
      integer txxx, txxy, txxz, txyy, txyz
      integer txzz, tyyy, tyyz, tyzz, tzzz
c
      parameter (txxx =  1)
      parameter (txxy =  2)
      parameter (txxz =  3)
      parameter (txyy =  4)
      parameter (txyz =  5)
      parameter (txzz =  6)
      parameter (tyyy =  7)
      parameter (tyyz =  8)
      parameter (tyzz =  9)
      parameter (tzzz = 10)
c
      integer tc3(3,3,3)
      data tc3/txxx,txxy,txxz,
     &         txxy,txyy,txyz,
     &         txxz,txyz,txzz,
     &         txxy,txyy,txyz,
     &         txyy,tyyy,tyyz,
     &         txyz,tyyz,tyzz,
     &         txxz,txyz,txzz,
     &         txyz,tyyz,tyzz,
     &         txzz,tyzz,tzzz/
c
c     The parameters defined here label the 1st and 2nd derivative
c     components of the density functional. Because a functional
c     may be purely local or of a GGA type the components are grouped
c     into derivative with respect to the density only and others.
c
c     1st derivatives:
c
      integer ira, irb
      integer igaa, igab, igbb
      integer ita, itb
c
c     2nd derivatives:
c
      integer irara, irarb, irbrb
      integer iragaa, iragab, iragbb, irbgaa, irbgab, irbgbb
      integer igaagaa, igaagab, igaagbb, igabgab, igabgbb, igbbgbb
      integer irata, iratb, irbta, irbtb
      integer igaata, igaatb, igabta, igabtb, igbbta, igbbtb
      integer itata, itatb, itbtb
c
      parameter (ira     = 1)
      parameter (irb     = 2)
c
      parameter (igaa    = 1)
      parameter (igab    = 2)
      parameter (igbb    = 3)
c
      parameter (ita     = 1)
      parameter (itb     = 2)
c
c
      parameter (irara   = 1)
      parameter (irarb   = 2)
      parameter (irbrb   = 3)
c
      parameter (iragaa  = 1)
      parameter (iragab  = 3)
      parameter (iragbb  = 4)
      parameter (irbgaa  = 5)
      parameter (irbgab  = 6)
      parameter (irbgbb  = 7)
      parameter (igaagaa = 2)
      parameter (igaagab = 8)
      parameter (igaagbb = 9)
      parameter (igabgab = 10)
      parameter (igabgbb = 11)
      parameter (igbbgbb = 12)
c
      parameter (irata   = 1)
      parameter (iratb   = 3)
      parameter (irbta   = 4)
      parameter (irbtb   = 5)
      parameter (itata   = 2)
      parameter (itatb   = 6)
      parameter (itbtb   = 7)
c
      parameter (igaata  = 1)
      parameter (igaatb  = 2)
      parameter (igabta  = 3)
      parameter (igabtb  = 4)
      parameter (igbbta  = 5)
      parameter (igbbtb  = 6)
c
      integer npts
      integer nao
      integer nvec
      integer naocc 
      integer npert
      integer mxp 
      integer n_active_bfn, n_active_atm
      integer active_bfn_list(n_active_bfn)
      integer active_bfn_indx(n_active_bfn)
      integer active_bfn_atms(n_active_atm)
      integer num_near, near_atom(num_near)
      logical gwt_avail_sw
      integer latm
      real*8 wt(npts)
      real*8 gwt(3,mxp,num_near)
      real*8 xc_vpt(mxp,2)
      real*8 xc_dvpt(mxp,3)
      real*8 xc_hpt(mxp,3)
      real*8 xc_dhpt(mxp,12)
      real*8 amo_val(npts,nvec)
      real*8 amo_grad(npts,nvec,3)
      real*8 bfng_val(mxp,nao,3)
      real*8 bfn_hess(mxp,nao,6)
c     real*8 bfn_3rd(mxp,nao,10)
      real*8 avec(nao,nvec)
      real*8 grho(mxp,2,3)
c
c     Outputs
c
      real*8 qa(nao*(nao+1)/2,npert)
c
c     Workspace
c
      real*8 dmo(npts,npert)
      real*8 dgmo(npts,npert,3)
c     real*8 ddmo(npts,6,ngridcent)
c     real*8 ddgmo(npts,10,ngridcent)
      real*8 drho(npts,npert)
c     real*8 ddrho(npts,npert*(npert+1)/2)
      real*8 dgrho(npts,npert,3)
c     real*8 ddgrho(npts,npert*(npert+1)/2,3)
      real*8 gmorho(npts,nvec)
      real*8 grdgmo(npts)
      real*8 t1(npts)
      real*8 t2(npts)
      real*8 t3(npts)
      real*8 t4(npts)
      real*8 tptmo(npts,nvec)
c
c     Local variables
c
      integer i, j
c     integer n, n1
      integer ipert, jpert, iprt
      integer ipt
      integer iatm
      integer ic
      integer iki, ikj
      integer iab, lab
      real*8 tmp
      integer ao_tag
      parameter(ao_tag=1)
c
c     Code
c
c     Construct dot product of grad rho and grad mo
c
      do i=1,nvec
         do ipt=1,npts
            gmorho(ipt,i)=
     &      (grho(ipt,1,1)*amo_grad(ipt,i,1)
     &      +grho(ipt,1,2)*amo_grad(ipt,i,2)
     &      +grho(ipt,1,3)*amo_grad(ipt,i,3))
            gmorho(ipt,i)=0.5d0*gmorho(ipt,i)
         enddo ! ipt
      enddo ! i
c
c     First add the contributions of the derivatives of the weights
c
      if (gwt_avail_sw) then
         do iatm=1,num_near
            do ic=1,3
               ipert=3*(near_atom(iatm)-1)+ic
               do j=1,nvec
                  ikj=iky(j)
                  do i=1,min(j,naocc)
                     do ipt=1,npts
                        qa(ikj+i,ipert)=qa(ikj+i,ipert)
     &                  +gwt(ic,ipt,iatm)*(
     &                   xc_vpt(ipt,ira)*amo_val(ipt,i)*amo_val(ipt,j)
     &                  +xc_dvpt(ipt,igaa)*gmorho(ipt,i)*amo_val(ipt,j)
     &                  +xc_dvpt(ipt,igaa)*amo_val(ipt,i)*gmorho(ipt,j))
                     enddo ! ipt
                  enddo ! i
               enddo ! j
            enddo ! ic
         enddo ! iatm
      endif
c
c     Now the derivatives of the Fock operator
c
      call aclear_dp(drho,npts*npert,0.0d0)
      call aclear_dp(dgrho,npts*npert*3,0.0d0)
c
      do ipt=1,npts
         t1(ipt)=wt(ipt)*xc_dvpt(ipt,igaa)
         t2(ipt)=wt(ipt)*xc_vpt(ipt,ira)
      enddo
      do i=1,nvec
         do ipt=1,npts
            tptmo(ipt,i)=t2(ipt)*amo_val(ipt,i)
     &                  +t1(ipt)*gmorho(ipt,i)
         enddo
      enddo
c
c     Construct the derivatives of the density (times 0.5)
c
      iki=0
      lab=3*(latm-1)
      do i=1,naocc
         call dmo_val_scr(npts,nao,nvec,npert,mxp,i,
     &                active_bfn_list,active_bfn_indx,n_active_bfn,
     &                bfng_val,avec,dmo)
         do ipert=1,3*n_active_atm
            do ipt=1,npts
               drho(ipt,ipert)=drho(ipt,ipert)
     &                        +amo_val(ipt,i)*dmo(ipt,ipert)
            enddo ! ipt
         enddo ! ipert
         do ic=1,3
            do ipert=1,3*n_active_atm
               do ipt=1,npts
                  dgrho(ipt,ipert,ic)=dgrho(ipt,ipert,ic)
     &            +amo_grad(ipt,i,ic)*dmo(ipt,ipert)
               enddo ! ipt
            enddo ! ipert
         enddo ! ic

         iki=iki+i-1
         do iatm=1,n_active_atm
            iab=3*(active_bfn_atms(iatm)-1)
            if (iab.ne.lab) then
               do ic=1,3
                  ipert=iab+ic
                  jpert=lab+ic
                  iprt =3*(iatm-1)+ic
                  do j=1,i-1
                     tmp=0.0d0
                     do ipt=1,npts
                        tmp=tmp+dmo(ipt,iprt)*tptmo(ipt,j)
                     enddo ! ipt
                     qa(iki+j,ipert)=qa(iki+j,ipert)+tmp
                     qa(iki+j,jpert)=qa(iki+j,jpert)-tmp
                  enddo ! j
                  tmp=0.0d0
                  do ipt=1,npts
                     tmp=tmp+dmo(ipt,iprt)*tptmo(ipt,i)
                  enddo ! ipt
                  tmp=2.0d0*tmp
                  qa(iki+i,ipert)=qa(iki+i,ipert)+tmp
                  qa(iki+i,jpert)=qa(iki+i,jpert)-tmp
                  ikj=iki
                  do j=i+1,nvec
                     ikj=ikj+j-1
                     tmp=0.0d0
                     do ipt=1,npts
                        tmp=tmp+dmo(ipt,iprt)*tptmo(ipt,j)
                     enddo
                     qa(ikj+i,ipert)=qa(ikj+i,ipert)+tmp
                     qa(ikj+i,jpert)=qa(ikj+i,jpert)-tmp
                  enddo ! j
               enddo ! ic
            endif
         enddo ! iatm
      enddo ! i
c
      iki=naocc*(naocc-1)/2
      do i=naocc+1,nvec
         call dmo_val_scr(npts,nao,nvec,npert,mxp,i,
     &                active_bfn_list,active_bfn_indx,n_active_bfn,
     &                bfng_val,avec,dmo)
         iki=iki+i-1
         do iatm=1,n_active_atm
            iab=3*(active_bfn_atms(iatm)-1)
            if (iab.ne.lab) then
               do ic=1,3
                  ipert=iab+ic
                  jpert=lab+ic
                  iprt =3*(iatm-1)+ic
                  do j=1,naocc
                     tmp=0.0d0
                     do ipt=1,npts
                        tmp=tmp+dmo(ipt,iprt)*tptmo(ipt,j)
                     enddo ! ipt
                     qa(iki+j,ipert)=qa(iki+j,ipert)+tmp
                     qa(iki+j,jpert)=qa(iki+j,jpert)-tmp
                  enddo ! j
               enddo ! ic
            endif
         enddo ! iatm
      enddo
c
      do ipt=1,npts
         t2(ipt)=4*wt(ipt)*xc_dhpt(ipt,iragaa)
         t3(ipt)=2*wt(ipt)*xc_hpt(ipt,irara)
      enddo
      do i=1,nvec
         do ipt=1,npts
            tptmo(ipt,i)=t3(ipt)*amo_val(ipt,i)
     &                  +t2(ipt)*gmorho(ipt,i)
         enddo
      enddo
      do iatm=1,n_active_atm
         iab=3*(active_bfn_atms(iatm)-1)
         if (iab.ne.lab) then
            do ic=1,3
               ipert=iab+ic
               jpert=lab+ic
               iprt =3*(iatm-1)+ic
               iki=0
               do i=1,nvec
                  iki=iki+i-1
                  do j=1,min(i,naocc)
                     tmp=0.0d0
                     do ipt=1,npts
                        tmp=tmp+tptmo(ipt,i)*amo_val(ipt,j)*
     &                          drho(ipt,iprt)
                     enddo ! ipt
                     qa(iki+j,ipert)=qa(iki+j,ipert)+tmp
                     qa(iki+j,jpert)=qa(iki+j,jpert)-tmp
                  enddo ! j
               enddo ! i
            enddo ! ic
         endif
      enddo ! iatm
      do i=1,nvec
         do ipt=1,npts
            tptmo(ipt,i)=t2(ipt)*amo_val(ipt,i)
         enddo
      enddo
      do iatm=1,n_active_atm
         iab=3*(active_bfn_atms(iatm)-1)
         if (iab.ne.lab) then
            do ic=1,3
               ipert=iab+ic
               jpert=lab+ic
               iprt =3*(iatm-1)+ic
               iki=0
               do i=1,nvec
                  iki=iki+i-1
                  do j=1,min(i,naocc)
                     tmp=0.0d0
                     do ipt=1,npts
                        tmp=tmp+tptmo(ipt,i)*gmorho(ipt,j)*
     &                          drho(ipt,iprt)
                     enddo ! ipt
                     qa(iki+j,ipert)=qa(iki+j,ipert)+tmp
                     qa(iki+j,jpert)=qa(iki+j,jpert)-tmp
                  enddo ! j
               enddo ! i
            enddo ! ic
         endif
      enddo ! iatm
c
      do i=1,nvec
         do ipt=1,npts
            tptmo(ipt,i)=t1(ipt)*amo_val(ipt,i)
         enddo
      enddo
      iki=0
      do i=1,naocc
         call dgmo_val_scr(npts,nao,nvec,npert,mxp,i,
     &                 active_bfn_list,active_bfn_indx,n_active_bfn,
     &                 bfn_hess,avec,dgmo)
         do ic=1,3
            do ipert=1,3*n_active_atm
               do ipt=1,npts
                  dgrho(ipt,ipert,ic)=dgrho(ipt,ipert,ic)
     &            +amo_val(ipt,i)*dgmo(ipt,ipert,ic)
               enddo ! ipt
            enddo ! ipert
         enddo ! ic
c
         iki=iki+i-1
         do iatm=1,n_active_atm
            iab=3*(active_bfn_atms(iatm)-1)
            if (iab.ne.lab) then
               do ic=1,3
                  ipert=iab+ic
                  jpert=lab+ic
                  iprt =3*(iatm-1)+ic
                  do ipt=1,npts
                     grdgmo(ipt) = grho(ipt,1,1)*dgmo(ipt,iprt,1)
     &                           + grho(ipt,1,2)*dgmo(ipt,iprt,2)
     &                           + grho(ipt,1,3)*dgmo(ipt,iprt,3)
                     grdgmo(ipt) = 0.5d0*grdgmo(ipt)
                  enddo
                  do j=1,i-1
                     tmp=0.0d0
                     do ipt=1,npts
                        tmp=tmp+tptmo(ipt,j)*grdgmo(ipt)
                     enddo ! ipt
                     qa(iki+j,ipert)=qa(iki+j,ipert)+tmp
                     qa(iki+j,jpert)=qa(iki+j,jpert)-tmp
                  enddo ! j
                  tmp=0.0d0
                  do ipt=1,npts
                     tmp=tmp+tptmo(ipt,i)*grdgmo(ipt)
                  enddo ! ipt
                  tmp=2.0d0*tmp
                  qa(iki+i,ipert)=qa(iki+i,ipert)+tmp
                  qa(iki+i,jpert)=qa(iki+i,jpert)-tmp
                  ikj=iki
                  do j=i+1,nvec
                     ikj=ikj+j-1
                     tmp=0.0d0
                     do ipt=1,npts
                        tmp=tmp+tptmo(ipt,j)*grdgmo(ipt)
                     enddo ! ipt
                     qa(ikj+i,ipert)=qa(ikj+i,ipert)+tmp
                     qa(ikj+i,jpert)=qa(ikj+i,jpert)-tmp
                  enddo ! j
               enddo ! ic
            endif
         enddo ! iatm
      enddo ! i
c
      iki=naocc*(naocc-1)/2
      do i=naocc+1,nvec
         call dgmo_val_scr(npts,nao,nvec,npert,mxp,i,
     &                 active_bfn_list,active_bfn_indx,n_active_bfn,
     &                 bfn_hess,avec,dgmo)
         iki=iki+i-1
         do iatm=1,n_active_atm
            iab=3*(active_bfn_atms(iatm)-1)
            if (iab.ne.lab) then
               do ic=1,3
                  ipert=iab+ic
                  jpert=lab+ic
                  iprt =3*(iatm-1)+ic
                  do ipt=1,npts
                     grdgmo(ipt) = grho(ipt,1,1)*dgmo(ipt,iprt,1)
     &                           + grho(ipt,1,2)*dgmo(ipt,iprt,2)
     &                           + grho(ipt,1,3)*dgmo(ipt,iprt,3)
                     grdgmo(ipt) = 0.5d0*grdgmo(ipt)
                  enddo
                  do j=1,naocc
                     tmp=0.0d0
                     do ipt=1,npts
                        tmp=tmp+tptmo(ipt,j)*grdgmo(ipt)
                     enddo ! ipt
                     qa(iki+j,ipert)=qa(iki+j,ipert)+tmp
                     qa(iki+j,jpert)=qa(iki+j,jpert)-tmp
                  enddo ! j
               enddo ! ic
            endif
         enddo ! iatm
      enddo ! i
c
      do ipt=1,npts
         t3(ipt)=4*wt(ipt)*xc_dhpt(ipt,igaagaa)
      enddo
      do i=1,nvec
         do ipt=1,npts
            tptmo(ipt,i)=t2(ipt)*amo_val(ipt,i)
     &                  +t3(ipt)*gmorho(ipt,i)
         enddo
      enddo
      do iatm=1,n_active_atm
         iab=3*(active_bfn_atms(iatm)-1)
         if (iab.ne.lab) then
            do ic=1,3
               ipert=iab+ic
               jpert=lab+ic
               iprt =3*(iatm-1)+ic
               do ipt=1,npts
                  grdgmo(ipt) = grho(ipt,1,1)*dgrho(ipt,iprt,1)
     &                        + grho(ipt,1,2)*dgrho(ipt,iprt,2)
     &                        + grho(ipt,1,3)*dgrho(ipt,iprt,3)
                  grdgmo(ipt) = 0.5d0*grdgmo(ipt)
               enddo ! ipt
               iki=0
               do i=1,nvec
                  iki=iki+i-1
                  do j=1,min(i,naocc)
                     tmp=0.0d0
                     do ipt=1,npts
                        tmp=tmp+tptmo(ipt,i)*amo_val(ipt,j)*grdgmo(ipt)
                     enddo ! ipt
                     qa(iki+j,ipert)=qa(iki+j,ipert)+tmp
                     qa(iki+j,jpert)=qa(iki+j,jpert)-tmp
                  enddo ! j
               enddo ! i
            enddo ! ic
         endif
      enddo ! iatm
      do i=1,nvec
         do ipt=1,npts
            tptmo(ipt,i)=t3(ipt)*amo_val(ipt,i)
         enddo
      enddo
      do iatm=1,n_active_atm
         iab=3*(active_bfn_atms(iatm)-1)
         if (iab.ne.lab) then
            do ic=1,3
               ipert=iab+ic
               jpert=lab+ic
               iprt =3*(iatm-1)+ic
               do ipt=1,npts
                  grdgmo(ipt) = grho(ipt,1,1)*dgrho(ipt,iprt,1)
     &                        + grho(ipt,1,2)*dgrho(ipt,iprt,2)
     &                        + grho(ipt,1,3)*dgrho(ipt,iprt,3)
                  grdgmo(ipt) = 0.5d0*grdgmo(ipt)
               enddo ! ipt
               iki=0
               do i=1,nvec
                  iki=iki+i-1
                  do j=1,min(i,naocc)
                     tmp=0.0d0
                     do ipt=1,npts
                        tmp=tmp+tptmo(ipt,i)*gmorho(ipt,j)*grdgmo(ipt)
                     enddo ! ipt
                     qa(iki+j,ipert)=qa(iki+j,ipert)+tmp
                     qa(iki+j,jpert)=qa(iki+j,jpert)-tmp
                  enddo ! j
               enddo ! i
            enddo ! ic
         endif
      enddo ! iatm
c
      do ipt=1,npts
         t1(ipt)=2*t1(ipt)
      enddo
      do i=1,nvec
         do ipt=1,npts
            tptmo(ipt,i)=t1(ipt)*amo_val(ipt,i)
         enddo
      enddo
      do iatm=1,n_active_atm
         iab=3*(active_bfn_atms(iatm)-1)
         if (iab.ne.lab) then
            do ic=1,3
               ipert=iab+ic
               jpert=lab+ic
               iprt =3*(iatm-1)+ic
               do i=1,nvec 
                  do ipt=1,npts
                     gmorho(ipt,i)=dgrho(ipt,iprt,1)*amo_grad(ipt,i,1)
     &                            +dgrho(ipt,iprt,2)*amo_grad(ipt,i,2)
     &                            +dgrho(ipt,iprt,3)*amo_grad(ipt,i,3)
                  enddo ! ipt
               enddo ! i
               iki=0
               do i=1,nvec
                  iki=iki+i-1
                  do j=1,min(i,naocc)
                     tmp=0.0d0
                     do ipt=1,npts
                        tmp=tmp+gmorho(ipt,i)*tptmo(ipt,j)
     &                         +gmorho(ipt,j)*tptmo(ipt,i)
                     enddo ! ipt
                     qa(iki+j,ipert)=qa(iki+j,ipert)+tmp
                     qa(iki+j,jpert)=qa(iki+j,jpert)-tmp
                  enddo ! j
               enddo ! i
            enddo ! ic
         endif
      enddo ! iatm
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine dmo_val(npts,nao,nvec,npert,ngridcentres,mxp,ith,
     &           bfng_val,avec,dmo)
      implicit none
c
c     Construct the derivative of the i-th MO with respect to the 
c     nuclear coordinates.
c
c     Inputs
c
C *Parameter module for CCP1/DFT            
C *Memory information		
C *----------------
C * max_block 	- maximum number of blocks which can be allocate
      integer max_block
      parameter(max_block=20)
C *							
C *Basis sets information                              
C *----------------------				
C *max_tag 	- maximum number of basis sets which can be alloc
C *max_atype	- maximum number of centre types	
C *max_gtype	- maximum number of grid types (at least every element 
C                 has different grid type)
C *max_grids    - maximum number of grids (different terms may have
C                 different grids, i.e. CPKS equations use a different
C                 grid than the one used for the KS-matrix)
C *max_shel	- maximum number of shells on centre
C *max_prm	- maximum number of primitives for any given centre
C *maxL		- maximum angular momentum allowed	
C *max_func	- maximum number of basis functions for any centre
      integer max_tag,max_atype,max_shel,max_prm,max_ang
      integer max_gtype,max_grids
      parameter(max_tag=3,max_atype=30,max_shel=500,max_prm=5000)
      parameter(max_ang=5,max_gtype=10,max_grids=2)
C *								
C *Geometry information					
C *--------------------					
C *max_atom	- maximum number of atoms in system
      integer max_atom
      parameter(max_atom=750)
C *
C *Accuracy information				
C ---------------------			
C *global_accuracy 	- global accuracy	
      real*8  global_accuracy
      parameter(global_accuracy=1.0d-14)

C *Grid information
c *----------------
c *maxradzn - The maximum number of radial zones. 
      integer maxgpt,maxrad,maxfpt,maxang,maxradzn,maxtablerows
      parameter(maxgpt=2900,maxrad=50,maxang=302,maxfpt=100)
      parameter(maxradzn=35,maxtablerows=7)
c
C *Basis sets include file
c
      real*8  alpha
      real*8  cont_coeff
      integer num_bset,bset_tags
      common/basis_sets/alpha(max_tag,max_atype,max_prm),
     &                  cont_coeff(max_tag,max_atype,max_prm,max_ang),
     &                  num_bset,bset_tags(max_tag)

      integer Ashl, Aprm, Abfn, totshl, totprm, totbfn, size_shlA
      integer size_basA, size_primA, maxi_shlA, maxi_basA, maxi_primA
      integer num_types, atom_tag , num_shl, atm_typ, nprim, angmom
      integer hybrid, pstart

      common/basis_size_info/Ashl(max_tag,max_atype),
     &                       Aprm(max_tag,max_atype),
     &                       Abfn(max_tag,max_atype),
     &                       totshl(max_tag),
     &                       totprm(max_tag),
     &                       totbfn(max_tag),
     &                       size_shlA,maxi_shlA,
     &                       size_basA,maxi_basA,
     &                       size_primA,maxi_primA
C
C Descriptions
C 
C  num_types		-	number of types of atoms for a given basis
C  atom_tag		-	list of atom type id numbers
C
      common/basis_cent_info/num_types(max_tag),
     &                       atom_tag(max_tag,max_atom)
C
C Descriptions
C
C  num_shl		-	number of shells for given atom type
C  atm_typ		-	atomic number of atom type
C  nprim		-  	number of primitives 
C  angmom		-	angular momentum
C  hybrid		-	level of hybrid shell. Is same as angmom if
C				shell is not a hybrid. For hybrid shells is
C				always less than angmom e.g. for an sp shell
C				angmom=2, hybrid=1
C  pstart		-	start of exponents and contraction coeffs 
C				contained in include basis.hf77, for a given
C				shell
C
      common/basis_cont_info/num_shl(max_tag,max_atype),
     &                       atm_typ(max_tag,max_atype),
     &                       nprim(max_tag,max_atype,max_shel),
     &                       angmom(max_tag,max_atype,max_shel),
     &                       hybrid(max_tag,max_atype,max_shel),
     &                       pstart(max_tag,max_atype,max_shel)

c
      integer npts
      integer nao
      integer nvec
      integer npert
      integer ngridcentres
      integer mxp
      integer ith
      real*8 bfng_val(mxp,nao,3)
      real*8 avec(nao,nvec)
c
c     Outputs
c
      real*8 dmo(npts,npert)
c
c     Local variables
c
      integer bstart
      integer bend
      integer latm
      integer lpert
      integer k
      integer atmt
      integer ao_tag
      parameter(ao_tag=1)
      integer mu
      integer ipt
c
c     Code
c
      call aclear_dp(dmo,npts*npert,0.0d0)
c
      bend = 0
      do latm=1,ngridcentres
         lpert=3*(latm-1)
         atmt=atom_tag(ao_tag,latm)
         bstart=bend+1
         bend=bend+Abfn(ao_tag,atmt)
         do k=1,3
            do mu=bstart,bend
               do ipt=1,npts
                  dmo(ipt,lpert+k)=dmo(ipt,lpert+k)
     &            - avec(mu,ith)*bfng_val(ipt,mu,k)
               enddo ! ipt
            enddo ! mu
         enddo ! k
      enddo ! latm
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine dgmo_val(npts,nao,nvec,npert,ngridcentres,mxp,ith,
     &           bfn_hess,avec,dgmo)
      implicit none
c
c     Construct the derivative of the gradient of the i-th MO with 
c     respect to the nuclear coordinates.
c
c     Inputs
c
C *Parameter module for CCP1/DFT            
C *Memory information		
C *----------------
C * max_block 	- maximum number of blocks which can be allocate
      integer max_block
      parameter(max_block=20)
C *							
C *Basis sets information                              
C *----------------------				
C *max_tag 	- maximum number of basis sets which can be alloc
C *max_atype	- maximum number of centre types	
C *max_gtype	- maximum number of grid types (at least every element 
C                 has different grid type)
C *max_grids    - maximum number of grids (different terms may have
C                 different grids, i.e. CPKS equations use a different
C                 grid than the one used for the KS-matrix)
C *max_shel	- maximum number of shells on centre
C *max_prm	- maximum number of primitives for any given centre
C *maxL		- maximum angular momentum allowed	
C *max_func	- maximum number of basis functions for any centre
      integer max_tag,max_atype,max_shel,max_prm,max_ang
      integer max_gtype,max_grids
      parameter(max_tag=3,max_atype=30,max_shel=500,max_prm=5000)
      parameter(max_ang=5,max_gtype=10,max_grids=2)
C *								
C *Geometry information					
C *--------------------					
C *max_atom	- maximum number of atoms in system
      integer max_atom
      parameter(max_atom=750)
C *
C *Accuracy information				
C ---------------------			
C *global_accuracy 	- global accuracy	
      real*8  global_accuracy
      parameter(global_accuracy=1.0d-14)

C *Grid information
c *----------------
c *maxradzn - The maximum number of radial zones. 
      integer maxgpt,maxrad,maxfpt,maxang,maxradzn,maxtablerows
      parameter(maxgpt=2900,maxrad=50,maxang=302,maxfpt=100)
      parameter(maxradzn=35,maxtablerows=7)
c
C *Basis sets include file
c
      real*8  alpha
      real*8  cont_coeff
      integer num_bset,bset_tags
      common/basis_sets/alpha(max_tag,max_atype,max_prm),
     &                  cont_coeff(max_tag,max_atype,max_prm,max_ang),
     &                  num_bset,bset_tags(max_tag)

      integer Ashl, Aprm, Abfn, totshl, totprm, totbfn, size_shlA
      integer size_basA, size_primA, maxi_shlA, maxi_basA, maxi_primA
      integer num_types, atom_tag , num_shl, atm_typ, nprim, angmom
      integer hybrid, pstart

      common/basis_size_info/Ashl(max_tag,max_atype),
     &                       Aprm(max_tag,max_atype),
     &                       Abfn(max_tag,max_atype),
     &                       totshl(max_tag),
     &                       totprm(max_tag),
     &                       totbfn(max_tag),
     &                       size_shlA,maxi_shlA,
     &                       size_basA,maxi_basA,
     &                       size_primA,maxi_primA
C
C Descriptions
C 
C  num_types		-	number of types of atoms for a given basis
C  atom_tag		-	list of atom type id numbers
C
      common/basis_cent_info/num_types(max_tag),
     &                       atom_tag(max_tag,max_atom)
C
C Descriptions
C
C  num_shl		-	number of shells for given atom type
C  atm_typ		-	atomic number of atom type
C  nprim		-  	number of primitives 
C  angmom		-	angular momentum
C  hybrid		-	level of hybrid shell. Is same as angmom if
C				shell is not a hybrid. For hybrid shells is
C				always less than angmom e.g. for an sp shell
C				angmom=2, hybrid=1
C  pstart		-	start of exponents and contraction coeffs 
C				contained in include basis.hf77, for a given
C				shell
C
      common/basis_cont_info/num_shl(max_tag,max_atype),
     &                       atm_typ(max_tag,max_atype),
     &                       nprim(max_tag,max_atype,max_shel),
     &                       angmom(max_tag,max_atype,max_shel),
     &                       hybrid(max_tag,max_atype,max_shel),
     &                       pstart(max_tag,max_atype,max_shel)

c
c     The parameters defined here label the 1st, 2nd and 3rd derivative
c     components of the basis functions.
c
c     2nd derivatives:
c
      integer hxx, hyy, hzz, hxy, hxz, hyz
c
      parameter (hxx = 1)
      parameter (hyy = 2)
      parameter (hzz = 3)
      parameter (hxy = 4)
      parameter (hxz = 5)
      parameter (hyz = 6)
c
      integer hcc(3,3)
      data hcc/hxx,hxy,hxz,
     &         hxy,hyy,hyz,
     &         hxz,hyz,hzz/
c
c     3rd derivatives:
c
      integer txxx, txxy, txxz, txyy, txyz
      integer txzz, tyyy, tyyz, tyzz, tzzz
c
      parameter (txxx =  1)
      parameter (txxy =  2)
      parameter (txxz =  3)
      parameter (txyy =  4)
      parameter (txyz =  5)
      parameter (txzz =  6)
      parameter (tyyy =  7)
      parameter (tyyz =  8)
      parameter (tyzz =  9)
      parameter (tzzz = 10)
c
      integer tc3(3,3,3)
      data tc3/txxx,txxy,txxz,
     &         txxy,txyy,txyz,
     &         txxz,txyz,txzz,
     &         txxy,txyy,txyz,
     &         txyy,tyyy,tyyz,
     &         txyz,tyyz,tyzz,
     &         txxz,txyz,txzz,
     &         txyz,tyyz,tyzz,
     &         txzz,tyzz,tzzz/
c
      integer npts
      integer nao
      integer nvec
      integer npert
      integer ngridcentres
      integer mxp
      integer ith
      real*8 bfn_hess(mxp,nao,6)
      real*8 avec(nao,nvec)
c
c     Outputs
c
      real*8 dgmo(npts,npert,3)
c
c     Local variables
c
      integer bstart
      integer bend
      integer latm
      integer lpert
      integer k,l
      integer atmt
      integer ao_tag
      parameter(ao_tag=1)
      integer mu
      integer ipt
c
c     Code
c
      call aclear_dp(dgmo,npts*npert*3,0.0d0)
c
      bend = 0
      do latm=1,ngridcentres
         lpert=3*(latm-1)
         atmt=atom_tag(ao_tag,latm)
         bstart=bend+1
         bend=bend+Abfn(ao_tag,atmt)
         do l=1,3
            do mu=bstart,bend
               do k=1,3
                  do ipt=1,npts
                     dgmo(ipt,lpert+k,l)=dgmo(ipt,lpert+k,l)
     &               - avec(mu,ith)*bfn_hess(ipt,mu,hcc(k,l))
                  enddo ! ipt
               enddo ! k
            enddo ! mu
         enddo ! l
      enddo ! latm
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine dmo_val_scr(npts,nao,nvec,npert,mxp,ith,
     &           active_bfn_list,active_bfn_indx,
     &           n_active_bfn,
     &           bfng_val,avec,dmo)
      implicit none
c
c     Construct the derivative of the i-th MO with respect to the 
c     nuclear coordinates.
c
c     The derivative of the MO is returned in dmo in a condensed
c     format. I.e. only the non-zero derivatives are stored.
c
c     Inputs
c
C *Parameter module for CCP1/DFT            
C *Memory information		
C *----------------
C * max_block 	- maximum number of blocks which can be allocate
      integer max_block
      parameter(max_block=20)
C *							
C *Basis sets information                              
C *----------------------				
C *max_tag 	- maximum number of basis sets which can be alloc
C *max_atype	- maximum number of centre types	
C *max_gtype	- maximum number of grid types (at least every element 
C                 has different grid type)
C *max_grids    - maximum number of grids (different terms may have
C                 different grids, i.e. CPKS equations use a different
C                 grid than the one used for the KS-matrix)
C *max_shel	- maximum number of shells on centre
C *max_prm	- maximum number of primitives for any given centre
C *maxL		- maximum angular momentum allowed	
C *max_func	- maximum number of basis functions for any centre
      integer max_tag,max_atype,max_shel,max_prm,max_ang
      integer max_gtype,max_grids
      parameter(max_tag=3,max_atype=30,max_shel=500,max_prm=5000)
      parameter(max_ang=5,max_gtype=10,max_grids=2)
C *								
C *Geometry information					
C *--------------------					
C *max_atom	- maximum number of atoms in system
      integer max_atom
      parameter(max_atom=750)
C *
C *Accuracy information				
C ---------------------			
C *global_accuracy 	- global accuracy	
      real*8  global_accuracy
      parameter(global_accuracy=1.0d-14)

C *Grid information
c *----------------
c *maxradzn - The maximum number of radial zones. 
      integer maxgpt,maxrad,maxfpt,maxang,maxradzn,maxtablerows
      parameter(maxgpt=2900,maxrad=50,maxang=302,maxfpt=100)
      parameter(maxradzn=35,maxtablerows=7)
c
C *Basis sets include file
c
      real*8  alpha
      real*8  cont_coeff
      integer num_bset,bset_tags
      common/basis_sets/alpha(max_tag,max_atype,max_prm),
     &                  cont_coeff(max_tag,max_atype,max_prm,max_ang),
     &                  num_bset,bset_tags(max_tag)

      integer Ashl, Aprm, Abfn, totshl, totprm, totbfn, size_shlA
      integer size_basA, size_primA, maxi_shlA, maxi_basA, maxi_primA
      integer num_types, atom_tag , num_shl, atm_typ, nprim, angmom
      integer hybrid, pstart

      common/basis_size_info/Ashl(max_tag,max_atype),
     &                       Aprm(max_tag,max_atype),
     &                       Abfn(max_tag,max_atype),
     &                       totshl(max_tag),
     &                       totprm(max_tag),
     &                       totbfn(max_tag),
     &                       size_shlA,maxi_shlA,
     &                       size_basA,maxi_basA,
     &                       size_primA,maxi_primA
C
C Descriptions
C 
C  num_types		-	number of types of atoms for a given basis
C  atom_tag		-	list of atom type id numbers
C
      common/basis_cent_info/num_types(max_tag),
     &                       atom_tag(max_tag,max_atom)
C
C Descriptions
C
C  num_shl		-	number of shells for given atom type
C  atm_typ		-	atomic number of atom type
C  nprim		-  	number of primitives 
C  angmom		-	angular momentum
C  hybrid		-	level of hybrid shell. Is same as angmom if
C				shell is not a hybrid. For hybrid shells is
C				always less than angmom e.g. for an sp shell
C				angmom=2, hybrid=1
C  pstart		-	start of exponents and contraction coeffs 
C				contained in include basis.hf77, for a given
C				shell
C
      common/basis_cont_info/num_shl(max_tag,max_atype),
     &                       atm_typ(max_tag,max_atype),
     &                       nprim(max_tag,max_atype,max_shel),
     &                       angmom(max_tag,max_atype,max_shel),
     &                       hybrid(max_tag,max_atype,max_shel),
     &                       pstart(max_tag,max_atype,max_shel)

c
      integer npts
      integer nao
      integer nvec
      integer npert
      integer mxp
      integer ith
      integer n_active_bfn
      integer active_bfn_list(n_active_bfn)
      integer active_bfn_indx(n_active_bfn)
      real*8 bfng_val(mxp,nao,3)
      real*8 avec(nao,nvec)
c
c     Outputs
c
      real*8 dmo(npts,npert)
c
c     Local variables
c
      integer lpert
      integer k
      integer ao_tag
      parameter(ao_tag=1)
      integer ipt
      integer ibfn
c
c     Code
c
      call aclear_dp(dmo,npts*npert,0.0d0)
c
      do ibfn=1,n_active_bfn
         lpert=3*(active_bfn_indx(ibfn)-1)
         do k=1,3
            do ipt=1,npts
               dmo(ipt,lpert+k)=dmo(ipt,lpert+k)
     &         - avec(active_bfn_list(ibfn),ith)*bfng_val(ipt,ibfn,k)
            enddo ! ipt
         enddo ! k
      enddo ! ibfn
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine dgmo_val_scr(npts,nao,nvec,npert,mxp,ith,
     &           active_bfn_list,active_bfn_indx,
     &           n_active_bfn,
     &           bfn_hess,avec,dgmo)
      implicit none
c
c     Construct the derivative of the gradient of the i-th MO with 
c     respect to the nuclear coordinates.
c
c     The derivative of the gradient of the MO is returned in dgmo in a 
c     condensed format. I.e. only the non-zero derivatives are stored.
c
c     Inputs
c
C *Parameter module for CCP1/DFT            
C *Memory information		
C *----------------
C * max_block 	- maximum number of blocks which can be allocate
      integer max_block
      parameter(max_block=20)
C *							
C *Basis sets information                              
C *----------------------				
C *max_tag 	- maximum number of basis sets which can be alloc
C *max_atype	- maximum number of centre types	
C *max_gtype	- maximum number of grid types (at least every element 
C                 has different grid type)
C *max_grids    - maximum number of grids (different terms may have
C                 different grids, i.e. CPKS equations use a different
C                 grid than the one used for the KS-matrix)
C *max_shel	- maximum number of shells on centre
C *max_prm	- maximum number of primitives for any given centre
C *maxL		- maximum angular momentum allowed	
C *max_func	- maximum number of basis functions for any centre
      integer max_tag,max_atype,max_shel,max_prm,max_ang
      integer max_gtype,max_grids
      parameter(max_tag=3,max_atype=30,max_shel=500,max_prm=5000)
      parameter(max_ang=5,max_gtype=10,max_grids=2)
C *								
C *Geometry information					
C *--------------------					
C *max_atom	- maximum number of atoms in system
      integer max_atom
      parameter(max_atom=750)
C *
C *Accuracy information				
C ---------------------			
C *global_accuracy 	- global accuracy	
      real*8  global_accuracy
      parameter(global_accuracy=1.0d-14)

C *Grid information
c *----------------
c *maxradzn - The maximum number of radial zones. 
      integer maxgpt,maxrad,maxfpt,maxang,maxradzn,maxtablerows
      parameter(maxgpt=2900,maxrad=50,maxang=302,maxfpt=100)
      parameter(maxradzn=35,maxtablerows=7)
c
C *Basis sets include file
c
      real*8  alpha
      real*8  cont_coeff
      integer num_bset,bset_tags
      common/basis_sets/alpha(max_tag,max_atype,max_prm),
     &                  cont_coeff(max_tag,max_atype,max_prm,max_ang),
     &                  num_bset,bset_tags(max_tag)

      integer Ashl, Aprm, Abfn, totshl, totprm, totbfn, size_shlA
      integer size_basA, size_primA, maxi_shlA, maxi_basA, maxi_primA
      integer num_types, atom_tag , num_shl, atm_typ, nprim, angmom
      integer hybrid, pstart

      common/basis_size_info/Ashl(max_tag,max_atype),
     &                       Aprm(max_tag,max_atype),
     &                       Abfn(max_tag,max_atype),
     &                       totshl(max_tag),
     &                       totprm(max_tag),
     &                       totbfn(max_tag),
     &                       size_shlA,maxi_shlA,
     &                       size_basA,maxi_basA,
     &                       size_primA,maxi_primA
C
C Descriptions
C 
C  num_types		-	number of types of atoms for a given basis
C  atom_tag		-	list of atom type id numbers
C
      common/basis_cent_info/num_types(max_tag),
     &                       atom_tag(max_tag,max_atom)
C
C Descriptions
C
C  num_shl		-	number of shells for given atom type
C  atm_typ		-	atomic number of atom type
C  nprim		-  	number of primitives 
C  angmom		-	angular momentum
C  hybrid		-	level of hybrid shell. Is same as angmom if
C				shell is not a hybrid. For hybrid shells is
C				always less than angmom e.g. for an sp shell
C				angmom=2, hybrid=1
C  pstart		-	start of exponents and contraction coeffs 
C				contained in include basis.hf77, for a given
C				shell
C
      common/basis_cont_info/num_shl(max_tag,max_atype),
     &                       atm_typ(max_tag,max_atype),
     &                       nprim(max_tag,max_atype,max_shel),
     &                       angmom(max_tag,max_atype,max_shel),
     &                       hybrid(max_tag,max_atype,max_shel),
     &                       pstart(max_tag,max_atype,max_shel)

c
c     The parameters defined here label the 1st, 2nd and 3rd derivative
c     components of the basis functions.
c
c     2nd derivatives:
c
      integer hxx, hyy, hzz, hxy, hxz, hyz
c
      parameter (hxx = 1)
      parameter (hyy = 2)
      parameter (hzz = 3)
      parameter (hxy = 4)
      parameter (hxz = 5)
      parameter (hyz = 6)
c
      integer hcc(3,3)
      data hcc/hxx,hxy,hxz,
     &         hxy,hyy,hyz,
     &         hxz,hyz,hzz/
c
c     3rd derivatives:
c
      integer txxx, txxy, txxz, txyy, txyz
      integer txzz, tyyy, tyyz, tyzz, tzzz
c
      parameter (txxx =  1)
      parameter (txxy =  2)
      parameter (txxz =  3)
      parameter (txyy =  4)
      parameter (txyz =  5)
      parameter (txzz =  6)
      parameter (tyyy =  7)
      parameter (tyyz =  8)
      parameter (tyzz =  9)
      parameter (tzzz = 10)
c
      integer tc3(3,3,3)
      data tc3/txxx,txxy,txxz,
     &         txxy,txyy,txyz,
     &         txxz,txyz,txzz,
     &         txxy,txyy,txyz,
     &         txyy,tyyy,tyyz,
     &         txyz,tyyz,tyzz,
     &         txxz,txyz,txzz,
     &         txyz,tyyz,tyzz,
     &         txzz,tyzz,tzzz/
c
      integer npts
      integer nao
      integer nvec
      integer npert
c     integer ngridcentres
      integer mxp
      integer ith
      integer n_active_bfn
      integer active_bfn_list(n_active_bfn)
      integer active_bfn_indx(n_active_bfn)
c     integer active_bfn_atms(n_active_atm)
      real*8 bfn_hess(mxp,nao,6)
      real*8 avec(nao,nvec)
c
c     Outputs
c
      real*8 dgmo(npts,npert,3)
c
c     Local variables
c
      integer lpert
      integer k,l
      integer ao_tag
      parameter(ao_tag=1)
      integer ipt
      integer ibfn
c
c     Code
c
      call aclear_dp(dgmo,npts*npert*3,0.0d0)
c
      do ibfn=1,n_active_bfn
         lpert=3*(active_bfn_indx(ibfn)-1)
         do l=1,3
            do k=1,3
               do ipt=1,npts
                  dgmo(ipt,lpert+k,l)=dgmo(ipt,lpert+k,l)
     &            - avec(active_bfn_list(ibfn),ith)*
     &              bfn_hess(ipt,ibfn,hcc(k,l))
               enddo ! ipt
            enddo ! k
         enddo ! l
      enddo ! latm
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine hess_dft_ao(rkstyp_sw,gradcorr_sw,gwt_sw,gwt_avail_sw,
     &     ao_tag,nao,npert,ngridcentres,npts,
     &     wt,gwt,ggwt,grho,drho,dgrho,ddrhoa,ddrhob,ddgrhoa,ddgrhob,
     &     bfn_val,bfng_val,bfn_hess,bfn_3rd,adens,bdens,
     &     xc_ept,xc_vpt,xc_dvpt,xc_hpt,xc_dhpt,iatm,hess,
     &     mxp)
      implicit none
c
c     Calculates the DFT contributions to the 2nd derivative of the 
c     energy. Only the explicit geometry dependencies are being 
c     differentiated.
c
c     Parameters:
c
C *Parameter module for CCP1/DFT            
C *Memory information		
C *----------------
C * max_block 	- maximum number of blocks which can be allocate
      integer max_block
      parameter(max_block=20)
C *							
C *Basis sets information                              
C *----------------------				
C *max_tag 	- maximum number of basis sets which can be alloc
C *max_atype	- maximum number of centre types	
C *max_gtype	- maximum number of grid types (at least every element 
C                 has different grid type)
C *max_grids    - maximum number of grids (different terms may have
C                 different grids, i.e. CPKS equations use a different
C                 grid than the one used for the KS-matrix)
C *max_shel	- maximum number of shells on centre
C *max_prm	- maximum number of primitives for any given centre
C *maxL		- maximum angular momentum allowed	
C *max_func	- maximum number of basis functions for any centre
      integer max_tag,max_atype,max_shel,max_prm,max_ang
      integer max_gtype,max_grids
      parameter(max_tag=3,max_atype=30,max_shel=500,max_prm=5000)
      parameter(max_ang=5,max_gtype=10,max_grids=2)
C *								
C *Geometry information					
C *--------------------					
C *max_atom	- maximum number of atoms in system
      integer max_atom
      parameter(max_atom=750)
C *
C *Accuracy information				
C ---------------------			
C *global_accuracy 	- global accuracy	
      real*8  global_accuracy
      parameter(global_accuracy=1.0d-14)

C *Grid information
c *----------------
c *maxradzn - The maximum number of radial zones. 
      integer maxgpt,maxrad,maxfpt,maxang,maxradzn,maxtablerows
      parameter(maxgpt=2900,maxrad=50,maxang=302,maxfpt=100)
      parameter(maxradzn=35,maxtablerows=7)
c
c     The parameters defined here label the 1st, 2nd and 3rd derivative
c     components of the basis functions.
c
c     2nd derivatives:
c
      integer hxx, hyy, hzz, hxy, hxz, hyz
c
      parameter (hxx = 1)
      parameter (hyy = 2)
      parameter (hzz = 3)
      parameter (hxy = 4)
      parameter (hxz = 5)
      parameter (hyz = 6)
c
      integer hcc(3,3)
      data hcc/hxx,hxy,hxz,
     &         hxy,hyy,hyz,
     &         hxz,hyz,hzz/
c
c     3rd derivatives:
c
      integer txxx, txxy, txxz, txyy, txyz
      integer txzz, tyyy, tyyz, tyzz, tzzz
c
      parameter (txxx =  1)
      parameter (txxy =  2)
      parameter (txxz =  3)
      parameter (txyy =  4)
      parameter (txyz =  5)
      parameter (txzz =  6)
      parameter (tyyy =  7)
      parameter (tyyz =  8)
      parameter (tyzz =  9)
      parameter (tzzz = 10)
c
      integer tc3(3,3,3)
      data tc3/txxx,txxy,txxz,
     &         txxy,txyy,txyz,
     &         txxz,txyz,txzz,
     &         txxy,txyy,txyz,
     &         txyy,tyyy,tyyz,
     &         txyz,tyyz,tyzz,
     &         txxz,txyz,txzz,
     &         txyz,tyyz,tyzz,
     &         txzz,tyzz,tzzz/
c
C *Basis sets include file
c
      real*8  alpha
      real*8  cont_coeff
      integer num_bset,bset_tags
      common/basis_sets/alpha(max_tag,max_atype,max_prm),
     &                  cont_coeff(max_tag,max_atype,max_prm,max_ang),
     &                  num_bset,bset_tags(max_tag)

      integer Ashl, Aprm, Abfn, totshl, totprm, totbfn, size_shlA
      integer size_basA, size_primA, maxi_shlA, maxi_basA, maxi_primA
      integer num_types, atom_tag , num_shl, atm_typ, nprim, angmom
      integer hybrid, pstart

      common/basis_size_info/Ashl(max_tag,max_atype),
     &                       Aprm(max_tag,max_atype),
     &                       Abfn(max_tag,max_atype),
     &                       totshl(max_tag),
     &                       totprm(max_tag),
     &                       totbfn(max_tag),
     &                       size_shlA,maxi_shlA,
     &                       size_basA,maxi_basA,
     &                       size_primA,maxi_primA
C
C Descriptions
C 
C  num_types		-	number of types of atoms for a given basis
C  atom_tag		-	list of atom type id numbers
C
      common/basis_cent_info/num_types(max_tag),
     &                       atom_tag(max_tag,max_atom)
C
C Descriptions
C
C  num_shl		-	number of shells for given atom type
C  atm_typ		-	atomic number of atom type
C  nprim		-  	number of primitives 
C  angmom		-	angular momentum
C  hybrid		-	level of hybrid shell. Is same as angmom if
C				shell is not a hybrid. For hybrid shells is
C				always less than angmom e.g. for an sp shell
C				angmom=2, hybrid=1
C  pstart		-	start of exponents and contraction coeffs 
C				contained in include basis.hf77, for a given
C				shell
C
      common/basis_cont_info/num_shl(max_tag,max_atype),
     &                       atm_typ(max_tag,max_atype),
     &                       nprim(max_tag,max_atype,max_shel),
     &                       angmom(max_tag,max_atype,max_shel),
     &                       hybrid(max_tag,max_atype,max_shel),
     &                       pstart(max_tag,max_atype,max_shel)

c
c     The parameters defined here label the 1st and 2nd derivative
c     components of the density functional. Because a functional
c     may be purely local or of a GGA type the components are grouped
c     into derivative with respect to the density only and others.
c
c     1st derivatives:
c
      integer ira, irb
      integer igaa, igab, igbb
      integer ita, itb
c
c     2nd derivatives:
c
      integer irara, irarb, irbrb
      integer iragaa, iragab, iragbb, irbgaa, irbgab, irbgbb
      integer igaagaa, igaagab, igaagbb, igabgab, igabgbb, igbbgbb
      integer irata, iratb, irbta, irbtb
      integer igaata, igaatb, igabta, igabtb, igbbta, igbbtb
      integer itata, itatb, itbtb
c
      parameter (ira     = 1)
      parameter (irb     = 2)
c
      parameter (igaa    = 1)
      parameter (igab    = 2)
      parameter (igbb    = 3)
c
      parameter (ita     = 1)
      parameter (itb     = 2)
c
c
      parameter (irara   = 1)
      parameter (irarb   = 2)
      parameter (irbrb   = 3)
c
      parameter (iragaa  = 1)
      parameter (iragab  = 3)
      parameter (iragbb  = 4)
      parameter (irbgaa  = 5)
      parameter (irbgab  = 6)
      parameter (irbgbb  = 7)
      parameter (igaagaa = 2)
      parameter (igaagab = 8)
      parameter (igaagbb = 9)
      parameter (igabgab = 10)
      parameter (igabgbb = 11)
      parameter (igbbgbb = 12)
c
      parameter (irata   = 1)
      parameter (iratb   = 3)
      parameter (irbta   = 4)
      parameter (irbtb   = 5)
      parameter (itata   = 2)
      parameter (itatb   = 6)
      parameter (itbtb   = 7)
c
      parameter (igaata  = 1)
      parameter (igaatb  = 2)
      parameter (igabta  = 3)
      parameter (igabtb  = 4)
      parameter (igbbta  = 5)
      parameter (igbbtb  = 6)
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
c     In variables:
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      logical rkstyp_sw
      logical gradcorr_sw
      logical gwt_sw
      logical gwt_avail_sw
      integer ao_tag
      integer nao
      integer npert
      integer ngridcentres
      integer npts
      integer iatm
      integer mxp
      real*8 wt(mxp)
      real*8 gwt(3,mxp,ngridcentres)
      real*8 ggwt(mxp,3*ngridcentres,3*ngridcentres)
      real*8 grho(mxp,2,3)
      real*8 drho(mxp,2,npert)     ! explicit derivatives only
      real*8 dgrho(mxp,2,3,npert)  ! explicit derivatives only
      real*8 ddrhoa(npts,9)        ! explicit derivatives only
      real*8 ddgrhoa(npts,3,9)     ! explicit derivatives only
      real*8 ddrhob(npts,9)        ! explicit derivatives only
      real*8 ddgrhob(npts,3,9)     ! explicit derivatives only
      real*8 bfn_val(mxp,nao)
      real*8 bfng_val(mxp,nao,3)
      real*8 bfn_hess(mxp,nao,6)
      real*8 bfn_3rd(mxp,nao,10)
      real*8 adens(nao*(nao+1)/2)
      real*8 bdens(nao*(nao+1)/2)
      real*8 xc_ept(mxp)
      real*8 xc_vpt(mxp,2)
      real*8 xc_dvpt(mxp,3)
      real*8 xc_hpt(mxp,3)
      real*8 xc_dhpt(mxp,12)
c
c     Out variables
c
      real*8 hess(npert,npert)
c
c     Local variables
c
      integer katm, latm, kpert, lpert, ipert
      integer kc, kx, lc, mc
      integer mu, nu
      integer klxx, ipt, iki
      integer bend, aend, batmt, aatmt, bstart, astart
      real*8 scr(9)
      real*8 gaada, gaadb
      real*8 gabda, gabdb
      real*8 gbbda, gbbdb
      real*8 gaadab, gabdab, gbbdab
c
c     Code:
c
      ipert = 3*(iatm-1)
      if (rkstyp_sw) then
c***  BEGIN OF CLOSED SHELL
         if (gwt_sw.and.gwt_avail_sw) then
cc
cc           First the terms involving the 2nd derivative of the weights
cc
c            do katm = 1, ngridcentres
c               kpert = 3*(katm-1)
cc              do latm = 1, katm-1
c               do latm = 1, ngridcentres
c                  lpert = 3*(latm-1)
c                  call aclear_dp(scr,9,0.0d0)
c                  do kc=1,3
c                     kx=3*(kc-1)
c                     do lc=1,3
cc                       klxx = iky(max(kpert+kc,lpert+lc))
cc    &                       + min(kpert+kc,lpert+lc)
cc                       klxx = iky(kpert+kc) + lpert+lc
c                        do ipt=1,npts
c                           scr(kx+lc)=scr(kx+lc)
c     &                     + ggwt(ipt,kpert+kc,lpert+lc)
c     &                       *xc_ept(ipt)
c                        enddo
c                     enddo
c                  enddo
c                  do kc=1,3
c                     kx=3*(kc-1)
c                     do lc=1,3
c                        hess(kpert+kc,lpert+lc) 
c     &                  = hess(kpert+kc,lpert+lc) + scr(kx+lc)
cc                       hess(lpert+lc,kpert+kc) 
cc    &                  = hess(lpert+lc,kpert+kc) + scr(kx+lc)
c                     enddo
c                  enddo
c               enddo
cc              call aclear_dp(scr,6,0.0d0)
cc              do kc = 1, 3
cc                 kx = kc*(kc-1)/2
cc                 do lc = 1, kc
cc                    klxx = iky(kpert+kc)+kpert+lc
cc                    do ipt = 1, npts
cc                       scr(kx+lc) = scr(kx+lc)
cc    &                  + ggwt(ipt,kpert+kc,kpert+lc)*xc_ept(ipt)
cc                    enddo
cc                 enddo
cc              enddo
cc              do kc = 1, 3
cc                 kx = kc*(kc-1)/2
cc                 do lc = 1, kc-1
cc                    hess(kpert+kc,kpert+lc) 
cc    &               = hess(kpert+kc,kpert+lc) + scr(kx+lc)
cc                    hess(kpert+lc,kpert+kc) 
cc    &               = hess(kpert+lc,kpert+kc) + scr(kx+lc)
cc                 enddo
cc                 hess(kpert+kc,kpert+kc) 
cc    &            = hess(kpert+kc,kpert+kc) + scr(kx+kc)
cc              enddo
c            enddo
c
c           Now the terms involving the 1st derivative of the weights
c
            do 100 katm = 1, ngridcentres
               if (katm.eq.iatm) goto 100
               kpert = 3*(katm-1)
               do latm = 1, ngridcentres
                  call aclear_dp(scr,9,0.0d0)
                  do kc = 1, 3
                     kx = 3*(kc-1)
                     do lc = 1, 3
                        do ipt = 1, npts
                           scr(kx+lc) = scr(kx+lc)
     &                     + 2*gwt(lc,ipt,latm)*
     &                       xc_vpt(ipt,1)*drho(ipt,1,kpert+kc)
                        enddo
                     enddo
                  enddo
                  if (gradcorr_sw) then
                     do kc = 1, 3
                        kx = 3*(kc-1)
                        do lc = 1, 3
                           do ipt = 1, npts
                              gaada 
     &                        = (dgrho(ipt,1,1,kpert+kc)*grho(ipt,1,1)
     &                         + dgrho(ipt,1,2,kpert+kc)*grho(ipt,1,2)
     &                         + dgrho(ipt,1,3,kpert+kc)*grho(ipt,1,3))
                              scr(kx+lc) = scr(kx+lc)
     &                        +gwt(lc,ipt,latm)*xc_dvpt(ipt,igaa)*gaada
                           enddo
                        enddo
                     enddo
                  endif
                  lpert = 3*(latm-1)
                  do kc = 1, 3
                     kx = 3*(kc-1)
                     do lc = 1, 3
                        hess(kpert+kc,lpert+lc)
     &                  = hess(kpert+kc,lpert+lc) + scr(kx+lc)
                        hess(ipert+kc,lpert+lc)
     &                  = hess(ipert+kc,lpert+lc) - scr(kx+lc)
                     enddo
                  enddo
               enddo
 100        continue
            do katm = 1, ngridcentres
               kpert = 3*(katm-1)
               do 110 latm = 1, ngridcentres
                  if (latm.eq.iatm) goto 110
                  lpert = 3*(latm-1)
                  call aclear_dp(scr,9,0.0d0)
                  do kc = 1, 3
                     kx = 3*(kc-1)
                     do lc = 1, 3
                        do ipt = 1, npts
                           scr(kx+lc) = scr(kx+lc)
     &                     + 2*gwt(kc,ipt,katm)*
     &                       xc_vpt(ipt,1)*drho(ipt,1,lpert+lc)
                        enddo
                     enddo
                  enddo
                  if (gradcorr_sw) then
                     do kc = 1, 3
                        kx = 3*(kc-1)
                        do lc = 1, 3
                           do ipt = 1, npts
                              gaada 
     &                        = (dgrho(ipt,1,1,lpert+lc)*grho(ipt,1,1)
     &                         + dgrho(ipt,1,2,lpert+lc)*grho(ipt,1,2)
     &                         + dgrho(ipt,1,3,lpert+lc)*grho(ipt,1,3))
                              scr(kx+lc) = scr(kx+lc)
     &                        +gwt(kc,ipt,katm)*xc_dvpt(ipt,igaa)*gaada
                           enddo
                        enddo
                     enddo
                  endif
                  do kc = 1, 3
                     kx = 3*(kc-1)
                     do lc = 1, 3
                        hess(kpert+kc,lpert+lc)
     &                  = hess(kpert+kc,lpert+lc) + scr(kx+lc)
                        hess(kpert+kc,ipert+lc)
     &                  = hess(kpert+kc,ipert+lc) - scr(kx+lc)
                     enddo
                  enddo
 110           continue
            enddo
         endif
c
         bend = 0
         do 120 katm = 1, ngridcentres
            batmt = atom_tag(ao_tag,katm)
            bstart= bend+1
            bend  = bend+Abfn(ao_tag,batmt)
            if (gwt_sw.and.katm.eq.iatm) goto 120
            kpert = 3*(katm-1)
            aend  = 0
            do 130 latm = 1, katm-1
               aatmt = atom_tag(ao_tag,latm)
               astart= aend+1
               aend  = aend+Abfn(ao_tag,aatmt)
               if (gwt_sw.and.latm.eq.iatm) goto 130
               lpert = 3*(latm-1)
c
               call aclear_dp(ddrhoa,npts*9,0.0d0)
               do mu=bstart,bend
                  do nu=astart,aend
                     iki = iky(mu)+nu
                     do kc = 1, 3
                        kx = 3*(kc-1)
                        do lc = 1, 3
                           klxx = kx+lc
                           do ipt=1,npts
                              ddrhoa(ipt,klxx) = ddrhoa(ipt,klxx)
     &                        + adens(iki)*bfng_val(ipt,mu,kc)*
     &                          bfng_val(ipt,nu,lc)
                           enddo ! ipt
                        enddo ! lc
                     enddo ! kc
                  enddo ! nu
               enddo ! mu
               if (gradcorr_sw) then
                  call aclear_dp(ddgrhoa,npts*3*9,0.0d0)
                  do mu=bstart,bend
                     do nu=astart,aend
                        iki = iky(mu)+nu
                        do kc = 1, 3
                           kx = 3*(kc-1)
                           do lc = 1, 3
                              klxx = kx+lc
                              do mc=1,3
                                 do ipt=1,npts
                                    ddgrhoa(ipt,mc,klxx)
     &                              = ddgrhoa(ipt,mc,klxx)
     &                              + adens(iki)*
     &                                bfn_hess(ipt,mu,hcc(mc,kc))*
     &                                bfng_val(ipt,nu,lc)
     &                              + adens(iki)*
     &                                bfng_val(ipt,mu,kc)*
     &                                bfn_hess(ipt,nu,hcc(mc,lc))
                                 enddo ! ipt
                              enddo ! mc
                           enddo ! lc
                        enddo ! kc
                     enddo ! nu
                  enddo ! mu
               endif ! gradcorr_sw
c
               call aclear_dp(scr,9,0.0d0)
               do kc=1,3
                  kx=3*(kc-1)
                  do lc=1,3
                     klxx = kx+lc
                     do ipt=1,npts
                        scr(kx+lc)=scr(kx+lc)
     &                  + 2*wt(ipt)*xc_hpt(ipt,irara)*
     &                    drho(ipt,1,kpert+kc)*drho(ipt,1,lpert+lc)
     &                  + 2*wt(ipt)*xc_vpt(ipt,ira)*ddrhoa(ipt,klxx)
                     enddo
                     if (gradcorr_sw) then
                        do ipt=1,npts
                           gaada 
     &                     =   (dgrho(ipt,1,1,kpert+kc)*grho(ipt,1,1)
     &                        + dgrho(ipt,1,2,kpert+kc)*grho(ipt,1,2)
     &                        + dgrho(ipt,1,3,kpert+kc)*grho(ipt,1,3))
                           gaadb 
     &                     =   (dgrho(ipt,1,1,lpert+lc)*grho(ipt,1,1)
     &                        + dgrho(ipt,1,2,lpert+lc)*grho(ipt,1,2)
     &                        + dgrho(ipt,1,3,lpert+lc)*grho(ipt,1,3))
                           gaadab 
     &                     =   (ddgrhoa(ipt,1,klxx)*grho(ipt,1,1)
     &                         +ddgrhoa(ipt,2,klxx)*grho(ipt,1,2)
     &                         +ddgrhoa(ipt,3,klxx)*grho(ipt,1,3))
     &                     + 2*(dgrho(ipt,1,1,kpert+kc)*
     &                          dgrho(ipt,1,1,lpert+lc)
     &                         +dgrho(ipt,1,2,kpert+kc)*
     &                          dgrho(ipt,1,2,lpert+lc)
     &                         +dgrho(ipt,1,3,kpert+kc)*
     &                          dgrho(ipt,1,3,lpert+lc))
                           scr(kx+lc)=scr(kx+lc)
     &                     + 2*wt(ipt)*xc_dhpt(ipt,iragaa)*
     &                       (drho(ipt,1,kpert+kc)*gaadb
     &                       +drho(ipt,1,lpert+lc)*gaada)
     &                     + wt(ipt)*xc_dhpt(ipt,igaagaa)*
     &                       gaada*gaadb
     &                     + wt(ipt)*xc_dvpt(ipt,igaa)*gaadab
                        enddo
                     endif
                  enddo
               enddo
               do kc=1,3
                  kx=3*(kc-1)
                  do lc=1,3
                     hess(kpert+kc,lpert+lc) = hess(kpert+kc,lpert+lc)
     &                                       + scr(kx+lc)
                     hess(lpert+lc,kpert+kc) = hess(lpert+lc,kpert+kc)
     &                                       + scr(kx+lc)
                     if (gwt_sw) then
                        hess(kpert+kc,ipert+lc) 
     &                  = hess(kpert+kc,ipert+lc) - scr(kx+lc)
                        hess(ipert+lc,kpert+kc) 
     &                  = hess(ipert+lc,kpert+kc) - scr(kx+lc)
c
                        hess(ipert+kc,lpert+lc) 
     &                  = hess(ipert+kc,lpert+lc) - scr(kx+lc)
                        hess(lpert+lc,ipert+kc) 
     &                  = hess(lpert+lc,ipert+kc) - scr(kx+lc)
c
                        hess(ipert+kc,ipert+lc) 
     &                  = hess(ipert+kc,ipert+lc) + scr(kx+lc)
                        hess(ipert+lc,ipert+kc) 
     &                  = hess(ipert+lc,ipert+kc) + scr(kx+lc)
                     endif
                  enddo
               enddo
 130        continue
c
            call aclear_dp(ddrhoa,npts*6,0.0d0)
c           2nd derivative of 1 basis function
            do mu=bstart,bend
               iki=iky(mu)
               do nu=1, mu
                  do kc=1,3
                     do lc=1,kc
                        klxx = kc*(kc-1)/2+lc
                        do ipt=1,npts
                           ddrhoa(ipt,klxx) = ddrhoa(ipt,klxx)
     &                     + adens(iki+nu)*bfn_hess(ipt,mu,hcc(lc,kc))*
     &                       bfn_val(ipt,nu)
                        enddo ! ipt
                     enddo ! lc
                  enddo ! kc
               enddo ! nu
               do nu=mu+1, nao
                  iki=iky(nu)
                  do kc=1,3
                     do lc=1,kc
                        klxx = kc*(kc-1)/2+lc
                        do ipt=1,npts
                           ddrhoa(ipt,klxx) = ddrhoa(ipt,klxx)
     &                     + adens(iki+mu)*bfn_hess(ipt,mu,hcc(lc,kc))*
     &                       bfn_val(ipt,nu)
                        enddo ! ipt
                     enddo ! lc
                  enddo ! kc
               enddo ! nu
            enddo ! mu
c           1st derivative of 2 basis functions
            do mu=bstart,bend
               do nu=bstart,mu
                  iki = iky(mu)+nu
                  do kc = 1, 3
                     do lc = 1, kc
                        klxx = kc*(kc-1)/2+lc
                        do ipt=1,npts
                           ddrhoa(ipt,klxx) = ddrhoa(ipt,klxx)
     &                     + adens(iki)*bfng_val(ipt,mu,kc)*
     &                       bfng_val(ipt,nu,lc)
                        enddo ! ipt
                     enddo ! lc
                  enddo ! kc
               enddo ! nu
               do nu=mu+1,bend
                  iki = iky(nu)+mu
                  do kc = 1, 3
                     do lc = 1, kc
                        klxx = kc*(kc-1)/2+lc
                        do ipt=1,npts
                           ddrhoa(ipt,klxx) = ddrhoa(ipt,klxx)
     &                     + adens(iki)*bfng_val(ipt,mu,kc)*
     &                       bfng_val(ipt,nu,lc)
                        enddo ! ipt
                     enddo ! lc
                  enddo ! kc
               enddo ! nu
            enddo ! mu
            if (gradcorr_sw) then
               call aclear_dp(ddgrhoa,npts*3*6,0.0d0)
c              2nd derivative of 1 basis function
               do mu=bstart,bend
                  iki=iky(mu)
                  do nu=1, mu
                     do kc=1,3
                        do lc=1,kc
                           klxx = kc*(kc-1)/2+lc
                           do mc=1,3
                              do ipt=1,npts
                                 ddgrhoa(ipt,mc,klxx)
     &                           = ddgrhoa(ipt,mc,klxx)
     &                           + adens(iki+nu)*(
     &                              bfn_3rd(ipt,mu,tc3(mc,lc,kc))
     &                             *bfn_val(ipt,nu)+
     &                              bfn_hess(ipt,mu,hcc(lc,kc))
     &                             *bfng_val(ipt,nu,mc))
                              enddo ! ipt
                           enddo ! mc
                        enddo ! lc
                     enddo ! kc
                  enddo ! nu
                  do nu=mu+1, nao
                     iki=iky(nu)
                     do kc=1,3
                        do lc=1,kc
                           klxx = kc*(kc-1)/2+lc
                           do mc=1,3
                              do ipt=1,npts
                                 ddgrhoa(ipt,mc,klxx)
     &                           = ddgrhoa(ipt,mc,klxx)
     &                           + adens(iki+mu)*(
     &                              bfn_3rd(ipt,mu,tc3(mc,lc,kc))
     &                             *bfn_val(ipt,nu)+
     &                              bfn_hess(ipt,mu,hcc(lc,kc))
     &                             *bfng_val(ipt,nu,mc))
                              enddo ! ipt
                           enddo ! mc
                        enddo ! lc
                     enddo ! kc
                  enddo ! nu
               enddo ! mu
c              1st derivatives of 2 basis functions
               do mu=bstart,bend
                  do nu=bstart, mu
                     iki = iky(mu)+nu
                     do kc = 1, 3
                        do lc = 1, kc
                           klxx = kc*(kc-1)/2+lc
                           do mc=1,3
                              do ipt=1,npts
                                 ddgrhoa(ipt,mc,klxx)
     &                           = ddgrhoa(ipt,mc,klxx)
     &                           + adens(iki)*
     &                             bfn_hess(ipt,mu,hcc(mc,kc))*
     &                             bfng_val(ipt,nu,lc)
     &                           + adens(iki)*
     &                             bfng_val(ipt,mu,kc)*
     &                             bfn_hess(ipt,nu,hcc(mc,lc))
                              enddo ! ipt
                           enddo ! m
                        enddo ! lc
                     enddo ! kc
                  enddo ! nu
                  do nu=mu+1, bend
                     iki = iky(nu)+mu
                     do kc = 1, 3
                        do lc = 1, kc
                           klxx = kc*(kc-1)/2+lc
                           do mc=1,3
                              do ipt=1,npts
                                 ddgrhoa(ipt,mc,klxx)
     &                           = ddgrhoa(ipt,mc,klxx)
     &                           + adens(iki)*
     &                             bfn_hess(ipt,mu,hcc(mc,kc))*
     &                             bfng_val(ipt,nu,lc)
     &                           + adens(iki)*
     &                             bfng_val(ipt,mu,kc)*
     &                             bfn_hess(ipt,nu,hcc(mc,lc))
                              enddo ! ipt
                           enddo ! mc
                        enddo ! lc
                     enddo ! kc
                  enddo ! nu
               enddo ! mu
            endif ! gradcorr_sw
c
            call aclear_dp(scr,6,0.0d0)
            do kc=1,3
               kx=kc*(kc-1)/2
               do lc=1,kc
                  klxx = kx+lc
                  do ipt=1,npts
                     scr(klxx)=scr(klxx)
     &               + 2*wt(ipt)*xc_hpt(ipt,irara)*
     &                 drho(ipt,1,kpert+kc)*drho(ipt,1,kpert+lc)
     &               + 2*wt(ipt)*xc_vpt(ipt,ira)*ddrhoa(ipt,klxx)
                  enddo
                  if (gradcorr_sw) then
                     do ipt=1,npts
                        gaada 
     &                  =   (dgrho(ipt,1,1,kpert+kc)*grho(ipt,1,1)
     &                     + dgrho(ipt,1,2,kpert+kc)*grho(ipt,1,2)
     &                     + dgrho(ipt,1,3,kpert+kc)*grho(ipt,1,3))
                        gaadb 
     &                  =   (dgrho(ipt,1,1,kpert+lc)*grho(ipt,1,1)
     &                     + dgrho(ipt,1,2,kpert+lc)*grho(ipt,1,2)
     &                     + dgrho(ipt,1,3,kpert+lc)*grho(ipt,1,3))
                        gaadab 
     &                  =   (ddgrhoa(ipt,1,klxx)*grho(ipt,1,1)
     &                      +ddgrhoa(ipt,2,klxx)*grho(ipt,1,2)
     &                      +ddgrhoa(ipt,3,klxx)*grho(ipt,1,3))
     &                  + 2*(dgrho(ipt,1,1,kpert+kc)*
     &                       dgrho(ipt,1,1,kpert+lc)
     &                      +dgrho(ipt,1,2,kpert+kc)*
     &                       dgrho(ipt,1,2,kpert+lc)
     &                      +dgrho(ipt,1,3,kpert+kc)*
     &                       dgrho(ipt,1,3,kpert+lc))
                        scr(klxx)=scr(klxx)
     &                  + 2*wt(ipt)*xc_dhpt(ipt,iragaa)*
     &                    (drho(ipt,1,kpert+kc)*gaadb
     &                    +drho(ipt,1,kpert+lc)*gaada)
     &                  + wt(ipt)*xc_dhpt(ipt,igaagaa)*
     &                    gaada*gaadb
     &                  + wt(ipt)*xc_dvpt(ipt,igaa)*gaadab
                     enddo
                  endif
               enddo
            enddo
            do kc=1,3
               kx=kc*(kc-1)/2
               do lc=1,kc-1
                  hess(kpert+kc,kpert+lc) = hess(kpert+kc,kpert+lc)
     &                                    + scr(kx+lc)
                  hess(kpert+lc,kpert+kc) = hess(kpert+lc,kpert+kc)
     &                                    + scr(kx+lc)
                  if (gwt_sw) then
                     hess(kpert+kc,ipert+lc) 
     &               = hess(kpert+kc,ipert+lc) - scr(kx+lc)
                     hess(ipert+lc,kpert+kc) 
     &               = hess(ipert+lc,kpert+kc) - scr(kx+lc)
c
                     hess(ipert+kc,kpert+lc) 
     &               = hess(ipert+kc,kpert+lc) - scr(kx+lc)
                     hess(kpert+lc,ipert+kc) 
     &               = hess(kpert+lc,ipert+kc) - scr(kx+lc)
c
                     hess(ipert+kc,ipert+lc) 
     &               = hess(ipert+kc,ipert+lc) + scr(kx+lc)
                     hess(ipert+lc,ipert+kc) 
     &               = hess(ipert+lc,ipert+kc) + scr(kx+lc)
                  endif
               enddo
               hess(kpert+kc,kpert+kc) = hess(kpert+kc,kpert+kc)
     &                                 + scr(kx+kc)
               if (gwt_sw) then
                  hess(kpert+kc,ipert+kc) = hess(kpert+kc,ipert+kc)
     &                                    - scr(kx+kc)
                  hess(ipert+kc,kpert+kc) = hess(ipert+kc,kpert+kc)
     &                                    - scr(kx+kc)
                  hess(ipert+kc,ipert+kc) = hess(ipert+kc,ipert+kc)
     &                                    + scr(kx+kc)
               endif
            enddo
 120     continue

c***  END OF CLOSED SHELL
      else ! rkstyp_sw
         if (gwt_sw.and.gwt_avail_sw) then
cc
cc           First the terms involving the 2nd derivative of the weights
cc
c            do katm = 1, ngridcentres
c               kpert = 3*(katm-1)
cc              do latm = 1, katm-1
c               do latm = 1, ngridcentres
c                  lpert = 3*(latm-1)
c                  call aclear_dp(scr,9,0.0d0)
c                  do kc=1,3
c                     kx=3*(kc-1)
c                     do lc=1,3
cc                       klxx = iky(max(kpert+kc,lpert+lc))
cc    &                       + min(kpert+kc,lpert+lc)
cc                       klxx = iky(kpert+kc) + lpert+lc
c                        do ipt=1,npts
c                           scr(kx+lc)=scr(kx+lc)
c     &                     + ggwt(ipt,kpert+kc,lpert+lc)
c     &                       *xc_ept(ipt)
c                        enddo
c                     enddo
c                  enddo
c                  do kc=1,3
c                     kx=3*(kc-1)
c                     do lc=1,3
c                        hess(kpert+kc,lpert+lc) 
c     &                  = hess(kpert+kc,lpert+lc) + scr(kx+lc)
cc                       hess(lpert+lc,kpert+kc) 
cc    &                  = hess(lpert+lc,kpert+kc) + scr(kx+lc)
c                     enddo
c                  enddo
c               enddo
cc              call aclear_dp(scr,6,0.0d0)
cc              do kc = 1, 3
cc                 kx = kc*(kc-1)/2
cc                 do lc = 1, kc
cc                    klxx = iky(kpert+kc)+kpert+lc
cc                    do ipt = 1, npts
cc                       scr(kx+lc) = scr(kx+lc)
cc    &                  + ggwt(ipt,kpert+kc,kpert+lc)*xc_ept(ipt)
cc                    enddo
cc                 enddo
cc              enddo
cc              do kc = 1, 3
cc                 kx = kc*(kc-1)/2
cc                 do lc = 1, kc-1
cc                    hess(kpert+kc,kpert+lc) 
cc    &               = hess(kpert+kc,kpert+lc) + scr(kx+lc)
cc                    hess(kpert+lc,kpert+kc) 
cc    &               = hess(kpert+lc,kpert+kc) + scr(kx+lc)
cc                 enddo
cc                 hess(kpert+kc,kpert+kc) 
cc    &            = hess(kpert+kc,kpert+kc) + scr(kx+kc)
cc              enddo
c            enddo
c
c           Now the terms involving the 1st derivative of the weights
c
            do 200 katm = 1, ngridcentres
               if (katm.eq.iatm) goto 200
               kpert = 3*(katm-1)
               do latm = 1, ngridcentres
                  call aclear_dp(scr,9,0.0d0)
                  do kc = 1, 3
                     kx = 3*(kc-1)
                     do lc = 1, 3
                        do ipt = 1, npts
                           scr(kx+lc) = scr(kx+lc)
     &                     + gwt(lc,ipt,latm)*
     &                       (xc_vpt(ipt,1)*drho(ipt,1,kpert+kc)
     &                       +xc_vpt(ipt,2)*drho(ipt,2,kpert+kc))
                        enddo
                     enddo
                  enddo
                  lpert = 3*(latm-1)
                  do kc = 1, 3
                     kx = 3*(kc-1)
                     do lc = 1, 3
                        hess(kpert+kc,lpert+lc)
     &                  = hess(kpert+kc,lpert+lc) + scr(kx+lc)
                        hess(ipert+kc,lpert+lc)
     &                  = hess(ipert+kc,lpert+lc) - scr(kx+lc)
                     enddo
                  enddo
               enddo
 200        continue
            do katm = 1, ngridcentres
               kpert = 3*(katm-1)
               do 210 latm = 1, ngridcentres
                  if (latm.eq.iatm) goto 210
                  lpert = 3*(latm-1)
                  call aclear_dp(scr,9,0.0d0)
                  do kc = 1, 3
                     kx = 3*(kc-1)
                     do lc = 1, 3
                        do ipt = 1, npts
                           scr(kx+lc) = scr(kx+lc)
     &                     + gwt(kc,ipt,katm)*
     &                       (xc_vpt(ipt,1)*drho(ipt,1,lpert+lc)
     &                       +xc_vpt(ipt,2)*drho(ipt,2,lpert+lc))
                        enddo
                     enddo
                  enddo
                  do kc = 1, 3
                     kx = 3*(kc-1)
                     do lc = 1, 3
                        hess(kpert+kc,lpert+lc)
     &                  = hess(kpert+kc,lpert+lc) + scr(kx+lc)
                        hess(kpert+kc,ipert+lc)
     &                  = hess(kpert+kc,ipert+lc) - scr(kx+lc)
                     enddo
                  enddo
 210           continue
            enddo
         endif

         bend = 0
         do 220 katm = 1, ngridcentres
            batmt = atom_tag(ao_tag,katm)
            bstart= bend+1
            bend  = bend+Abfn(ao_tag,batmt)
            if (gwt_sw.and.katm.eq.iatm) goto 220
            kpert = 3*(katm-1)
            aend  = 0
            do 230 latm = 1, katm-1
               aatmt = atom_tag(ao_tag,latm)
               astart= aend+1
               aend  = aend+Abfn(ao_tag,aatmt)
               if (gwt_sw.and.latm.eq.iatm) goto 230
               lpert = 3*(latm-1)
c
               call aclear_dp(ddrhoa,npts*9,0.0d0)
               call aclear_dp(ddrhob,npts*9,0.0d0)
               do mu=bstart,bend
                  do nu=astart,aend
                     iki = iky(mu)+nu
                     do kc = 1, 3
                        kx = 3*(kc-1)
                        do lc = 1, 3
                           klxx = kx+lc
                           do ipt=1,npts
                              ddrhoa(ipt,klxx) = ddrhoa(ipt,klxx)
     &                        + adens(iki)*bfng_val(ipt,mu,kc)*
     &                          bfng_val(ipt,nu,lc)
                              ddrhob(ipt,klxx) = ddrhob(ipt,klxx)
     &                        + bdens(iki)*bfng_val(ipt,mu,kc)*
     &                          bfng_val(ipt,nu,lc)
                           enddo ! ipt
                        enddo ! lc
                     enddo ! kc
                  enddo ! nu
               enddo ! mu
               if (gradcorr_sw) then
                  call aclear_dp(ddgrhoa,npts*3*9,0.0d0)
                  call aclear_dp(ddgrhob,npts*3*9,0.0d0)
                  do mu=bstart,bend
                     do nu=astart,aend
                        iki = iky(mu)+nu
                        do kc = 1, 3
                           kx = 3*(kc-1)
                           do lc = 1, 3
                              klxx = kx+lc
                              do mc=1,3
                                 do ipt=1,npts
                                    ddgrhoa(ipt,mc,klxx)
     &                              = ddgrhoa(ipt,mc,klxx)
     &                              + adens(iki)*
     &                                bfn_hess(ipt,mu,hcc(mc,kc))*
     &                                bfng_val(ipt,nu,lc)
     &                              + adens(iki)*
     &                                bfng_val(ipt,mu,kc)*
     &                                bfn_hess(ipt,nu,hcc(mc,lc))
                                    ddgrhob(ipt,mc,klxx)
     &                              = ddgrhob(ipt,mc,klxx)
     &                              + bdens(iki)*
     &                                bfn_hess(ipt,mu,hcc(mc,kc))*
     &                                bfng_val(ipt,nu,lc)
     &                              + bdens(iki)*
     &                                bfng_val(ipt,mu,kc)*
     &                                bfn_hess(ipt,nu,hcc(mc,lc))
                                 enddo ! ipt
                              enddo ! mc
                           enddo ! lc
                        enddo ! kc
                     enddo ! nu
                  enddo ! mu
               endif ! gradcorr_sw
c
               call aclear_dp(scr,9,0.0d0)
               do kc=1,3
                  kx=3*(kc-1)
                  do lc=1,3
                     klxx = kx+lc
                     do ipt=1,npts
                        scr(kx+lc)=scr(kx+lc)
     &                  + wt(ipt)*xc_hpt(ipt,irara)*
     &                    drho(ipt,1,kpert+kc)*drho(ipt,1,lpert+lc)
     &                  + wt(ipt)*xc_hpt(ipt,irbrb)*
     &                    drho(ipt,2,kpert+kc)*drho(ipt,2,lpert+lc)
     &                  + wt(ipt)*xc_hpt(ipt,irarb)*
     &                    (drho(ipt,1,kpert+kc)*drho(ipt,2,lpert+lc)+
     &                     drho(ipt,2,kpert+kc)*drho(ipt,1,lpert+lc))
     &                  + wt(ipt)*xc_vpt(ipt,ira)*ddrhoa(ipt,klxx)
     &                  + wt(ipt)*xc_vpt(ipt,irb)*ddrhob(ipt,klxx)
                     enddo
                     if (gradcorr_sw) then
                        do ipt=1,npts
                           gaada 
     &                     = 2*(dgrho(ipt,1,1,kpert+kc)*grho(ipt,1,1)
     &                        + dgrho(ipt,1,2,kpert+kc)*grho(ipt,1,2)
     &                        + dgrho(ipt,1,3,kpert+kc)*grho(ipt,1,3))
                           gbbda 
     &                     = 2*(dgrho(ipt,2,1,kpert+kc)*grho(ipt,2,1)
     &                        + dgrho(ipt,2,2,kpert+kc)*grho(ipt,2,2)
     &                        + dgrho(ipt,2,3,kpert+kc)*grho(ipt,2,3))
                           gabda 
     &                     =    dgrho(ipt,1,1,kpert+kc)*grho(ipt,2,1)
     &                        + dgrho(ipt,1,2,kpert+kc)*grho(ipt,2,2)
     &                        + dgrho(ipt,1,3,kpert+kc)*grho(ipt,2,3)
     &                        + dgrho(ipt,2,1,kpert+kc)*grho(ipt,1,1)
     &                        + dgrho(ipt,2,2,kpert+kc)*grho(ipt,1,2)
     &                        + dgrho(ipt,2,3,kpert+kc)*grho(ipt,1,3)
                           gaadb 
     &                     = 2*(dgrho(ipt,1,1,lpert+lc)*grho(ipt,1,1)
     &                        + dgrho(ipt,1,2,lpert+lc)*grho(ipt,1,2)
     &                        + dgrho(ipt,1,3,lpert+lc)*grho(ipt,1,3))
                           gbbdb 
     &                     = 2*(dgrho(ipt,2,1,lpert+lc)*grho(ipt,2,1)
     &                        + dgrho(ipt,2,2,lpert+lc)*grho(ipt,2,2)
     &                        + dgrho(ipt,2,3,lpert+lc)*grho(ipt,2,3))
                           gabdb 
     &                     =    dgrho(ipt,1,1,lpert+lc)*grho(ipt,2,1)
     &                        + dgrho(ipt,1,2,lpert+lc)*grho(ipt,2,2)
     &                        + dgrho(ipt,1,3,lpert+lc)*grho(ipt,2,3)
     &                        + dgrho(ipt,2,1,lpert+lc)*grho(ipt,1,1)
     &                        + dgrho(ipt,2,2,lpert+lc)*grho(ipt,1,2)
     &                        + dgrho(ipt,2,3,lpert+lc)*grho(ipt,1,3)
                           gaadab 
     &                     = 2*(ddgrhoa(ipt,1,klxx)*grho(ipt,1,1)
     &                         +ddgrhoa(ipt,2,klxx)*grho(ipt,1,2)
     &                         +ddgrhoa(ipt,3,klxx)*grho(ipt,1,3)
     &                         +dgrho(ipt,1,1,kpert+kc)*
     &                          dgrho(ipt,1,1,lpert+lc)
     &                         +dgrho(ipt,1,2,kpert+kc)*
     &                          dgrho(ipt,1,2,lpert+lc)
     &                         +dgrho(ipt,1,3,kpert+kc)*
     &                          dgrho(ipt,1,3,lpert+lc))
                           gbbdab 
     &                     = 2*(ddgrhob(ipt,1,klxx)*grho(ipt,2,1)
     &                         +ddgrhob(ipt,2,klxx)*grho(ipt,2,2)
     &                         +ddgrhob(ipt,3,klxx)*grho(ipt,2,3)
     &                         +dgrho(ipt,2,1,kpert+kc)*
     &                          dgrho(ipt,2,1,lpert+lc)
     &                         +dgrho(ipt,2,2,kpert+kc)*
     &                          dgrho(ipt,2,2,lpert+lc)
     &                         +dgrho(ipt,2,3,kpert+kc)*
     &                          dgrho(ipt,2,3,lpert+lc))
                           gabdab 
     &                     =   (ddgrhoa(ipt,1,klxx)*grho(ipt,2,1)
     &                         +ddgrhoa(ipt,2,klxx)*grho(ipt,2,2)
     &                         +ddgrhoa(ipt,3,klxx)*grho(ipt,2,3)
     &                         +ddgrhob(ipt,1,klxx)*grho(ipt,1,1)
     &                         +ddgrhob(ipt,2,klxx)*grho(ipt,1,2)
     &                         +ddgrhob(ipt,3,klxx)*grho(ipt,1,3)
     &                         +dgrho(ipt,1,1,kpert+kc)*
     &                          dgrho(ipt,2,1,lpert+lc)
     &                         +dgrho(ipt,1,2,kpert+kc)*
     &                          dgrho(ipt,2,2,lpert+lc)
     &                         +dgrho(ipt,1,3,kpert+kc)*
     &                          dgrho(ipt,2,3,lpert+lc)
     &                         +dgrho(ipt,2,1,kpert+kc)*
     &                          dgrho(ipt,1,1,lpert+lc)
     &                         +dgrho(ipt,2,2,kpert+kc)*
     &                          dgrho(ipt,1,2,lpert+lc)
     &                         +dgrho(ipt,2,3,kpert+kc)*
     &                          dgrho(ipt,1,3,lpert+lc))
                           scr(kx+lc)=scr(kx+lc)
     &                     + wt(ipt)*xc_dhpt(ipt,iragaa)*
     &                       (drho(ipt,1,kpert+kc)*gaadb
     &                       +drho(ipt,1,lpert+lc)*gaada)
     &                     + wt(ipt)*xc_dhpt(ipt,iragab)*
     &                       (drho(ipt,1,kpert+kc)*gabdb
     &                       +drho(ipt,1,lpert+lc)*gabda)
     &                     + wt(ipt)*xc_dhpt(ipt,iragbb)*
     &                       (drho(ipt,1,kpert+kc)*gbbdb
     &                       +drho(ipt,1,lpert+lc)*gbbda)
     &                     + wt(ipt)*xc_dhpt(ipt,irbgaa)*
     &                       (drho(ipt,2,kpert+kc)*gaadb
     &                       +drho(ipt,2,lpert+lc)*gaada)
     &                     + wt(ipt)*xc_dhpt(ipt,irbgab)*
     &                       (drho(ipt,2,kpert+kc)*gabdb
     &                       +drho(ipt,2,lpert+lc)*gabda)
     &                     + wt(ipt)*xc_dhpt(ipt,irbgbb)*
     &                       (drho(ipt,2,kpert+kc)*gbbdb
     &                       +drho(ipt,2,lpert+lc)*gbbda)
     &                     + wt(ipt)*xc_dhpt(ipt,igaagaa)*
     &                       gaada*gaadb
     &                     + wt(ipt)*xc_dhpt(ipt,igaagab)*
     &                       gaada*gabdb
     &                     + wt(ipt)*xc_dhpt(ipt,igaagbb)*
     &                       gaada*gbbdb
     &                     + wt(ipt)*xc_dhpt(ipt,igaagab)*
     &                       gabda*gaadb
     &                     + wt(ipt)*xc_dhpt(ipt,igabgab)*
     &                       gabda*gabdb
     &                     + wt(ipt)*xc_dhpt(ipt,igabgbb)*
     &                       gabda*gbbdb
     &                     + wt(ipt)*xc_dhpt(ipt,igaagbb)*
     &                       gbbda*gaadb
     &                     + wt(ipt)*xc_dhpt(ipt,igabgbb)*
     &                       gbbda*gabdb
     &                     + wt(ipt)*xc_dhpt(ipt,igbbgbb)*
     &                       gbbda*gbbdb
     &                     + wt(ipt)*xc_dvpt(ipt,igaa)*gaadab
     &                     + wt(ipt)*xc_dvpt(ipt,igab)*gabdab
     &                     + wt(ipt)*xc_dvpt(ipt,igbb)*gbbdab
                        enddo
                     endif
                  enddo
               enddo
               do kc=1,3
                  kx=3*(kc-1)
                  do lc=1,3
                     hess(kpert+kc,lpert+lc) = hess(kpert+kc,lpert+lc)
     &                                       + scr(kx+lc)
                     hess(lpert+lc,kpert+kc) = hess(lpert+lc,kpert+kc)
     &                                       + scr(kx+lc)
                     if (gwt_sw) then
                        hess(kpert+kc,ipert+lc) 
     &                  = hess(kpert+kc,ipert+lc) - scr(kx+lc)
                        hess(ipert+lc,kpert+kc) 
     &                  = hess(ipert+lc,kpert+kc) - scr(kx+lc)
c
                        hess(ipert+kc,lpert+lc) 
     &                  = hess(ipert+kc,lpert+lc) - scr(kx+lc)
                        hess(lpert+lc,ipert+kc) 
     &                  = hess(lpert+lc,ipert+kc) - scr(kx+lc)
c
                        hess(ipert+kc,ipert+lc) 
     &                  = hess(ipert+kc,ipert+lc) + scr(kx+lc)
                        hess(ipert+lc,ipert+kc) 
     &                  = hess(ipert+lc,ipert+kc) + scr(kx+lc)
                     endif
                  enddo
               enddo
 230        continue
c
            call aclear_dp(ddrhoa,npts*6,0.0d0)
            call aclear_dp(ddrhob,npts*6,0.0d0)
c           2nd derivative of 1 basis function
            do mu=bstart,bend
               iki=iky(mu)
               do nu=1, mu
                  do kc=1,3
                     do lc=1,kc
                        klxx = kc*(kc-1)/2+lc
                        do ipt=1,npts
                           ddrhoa(ipt,klxx) = ddrhoa(ipt,klxx)
     &                     + adens(iki+nu)*bfn_hess(ipt,mu,hcc(lc,kc))*
     &                       bfn_val(ipt,nu)
                           ddrhob(ipt,klxx) = ddrhob(ipt,klxx)
     &                     + bdens(iki+nu)*bfn_hess(ipt,mu,hcc(lc,kc))*
     &                       bfn_val(ipt,nu)
                        enddo ! ipt
                     enddo ! l
                  enddo ! k
               enddo ! nu
               do nu=mu+1, nao
                  iki=iky(nu)
                  do kc=1,3
                     do lc=1,kc
                        klxx = kc*(kc-1)/2+lc
                        do ipt=1,npts
                           ddrhoa(ipt,klxx) = ddrhoa(ipt,klxx)
     &                     + adens(iki+mu)*bfn_hess(ipt,mu,hcc(lc,kc))*
     &                       bfn_val(ipt,nu)
                           ddrhob(ipt,klxx) = ddrhob(ipt,klxx)
     &                     + bdens(iki+mu)*bfn_hess(ipt,mu,hcc(lc,kc))*
     &                       bfn_val(ipt,nu)
                        enddo ! ipt
                     enddo ! lc
                  enddo ! kc
               enddo ! nu
            enddo ! mu
c           1st derivative of 2 basis functions
            do mu=bstart,bend
               do nu=bstart,mu
                  iki = iky(mu)+nu
                  do kc = 1, 3
                     do lc = 1, kc
                        klxx = kc*(kc-1)/2+lc
                        do ipt=1,npts
                           ddrhoa(ipt,klxx) = ddrhoa(ipt,klxx)
     &                     + adens(iki)*bfng_val(ipt,mu,kc)*
     &                       bfng_val(ipt,nu,lc)
                           ddrhob(ipt,klxx) = ddrhob(ipt,klxx)
     &                     + bdens(iki)*bfng_val(ipt,mu,kc)*
     &                       bfng_val(ipt,nu,lc)
                        enddo ! ipt
                     enddo ! lc
                  enddo ! kc
               enddo ! nu
               do nu=mu+1,bend
                  iki = iky(nu)+mu
                  do kc = 1, 3
                     do lc = 1, kc
                        klxx = kc*(kc-1)/2+lc
                        do ipt=1,npts
                           ddrhoa(ipt,klxx) = ddrhoa(ipt,klxx)
     &                     + adens(iki)*bfng_val(ipt,mu,kc)*
     &                       bfng_val(ipt,nu,lc)
                           ddrhob(ipt,klxx) = ddrhob(ipt,klxx)
     &                     + bdens(iki)*bfng_val(ipt,mu,kc)*
     &                       bfng_val(ipt,nu,lc)
                        enddo ! ipt
                     enddo ! lc
                  enddo ! kc
               enddo ! nu
            enddo ! mu
            if (gradcorr_sw) then
               call aclear_dp(ddgrhoa,npts*3*6,0.0d0)
               call aclear_dp(ddgrhob,npts*3*6,0.0d0)
c              2nd derivative of 1 basis function
               do mu=bstart,bend
                  iki=iky(mu)
                  do nu=1, mu
                     do kc=1,3
                        do lc=1,kc
                           klxx = kc*(kc-1)/2+lc
                           do mc=1,3
                              do ipt=1,npts
                                 ddgrhoa(ipt,mc,klxx)
     &                           = ddgrhoa(ipt,mc,klxx)
     &                           + adens(iki+nu)*(
     &                              bfn_3rd(ipt,mu,tc3(mc,lc,kc))
     &                             *bfn_val(ipt,nu)+
     &                              bfn_hess(ipt,mu,hcc(lc,kc))
     &                             *bfng_val(ipt,nu,mc))
                                 ddgrhob(ipt,mc,klxx)
     &                           = ddgrhob(ipt,mc,klxx)
     &                           + bdens(iki+nu)*(
     &                              bfn_3rd(ipt,mu,tc3(mc,lc,kc))
     &                             *bfn_val(ipt,nu)+
     &                              bfn_hess(ipt,mu,hcc(lc,kc))
     &                             *bfng_val(ipt,nu,mc))
                              enddo ! ipt
                           enddo ! mc
                        enddo ! lc
                     enddo ! kc
                  enddo ! nu
                  do nu=mu+1, nao
                     iki=iky(nu)
                     do kc=1,3
                        do lc=1,kc
                           klxx = kc*(kc-1)/2+lc
                           do mc=1,3
                              do ipt=1,npts
                                 ddgrhoa(ipt,mc,klxx)
     &                           = ddgrhoa(ipt,mc,klxx)
     &                           + adens(iki+mu)*(
     &                              bfn_3rd(ipt,mu,tc3(mc,lc,kc))
     &                             *bfn_val(ipt,nu)+
     &                              bfn_hess(ipt,mu,hcc(lc,kc))
     &                             *bfng_val(ipt,nu,mc))
                                 ddgrhob(ipt,mc,klxx)
     &                           = ddgrhob(ipt,mc,klxx)
     &                           + bdens(iki+mu)*(
     &                              bfn_3rd(ipt,mu,tc3(mc,lc,kc))
     &                             *bfn_val(ipt,nu)+
     &                              bfn_hess(ipt,mu,hcc(lc,kc))
     &                             *bfng_val(ipt,nu,mc))
                              enddo ! ipt
                           enddo ! mc
                        enddo ! lc
                     enddo ! kc
                  enddo ! nu
               enddo ! mu
c              1st derivatives of 2 basis functions
               do mu=bstart,bend
                  do nu=bstart, mu
                     iki = iky(mu)+nu
                     do kc = 1, 3
                        do lc = 1, kc
                           klxx = kc*(kc-1)/2+lc
                           do mc=1,3
                              do ipt=1,npts
                                 ddgrhoa(ipt,mc,klxx)
     &                           = ddgrhoa(ipt,mc,klxx)
     &                           + adens(iki)*
     &                             bfn_hess(ipt,mu,hcc(mc,kc))*
     &                             bfng_val(ipt,nu,lc)
     &                           + adens(iki)*
     &                             bfng_val(ipt,mu,kc)*
     &                             bfn_hess(ipt,nu,hcc(mc,lc))
                                 ddgrhob(ipt,mc,klxx)
     &                           = ddgrhob(ipt,mc,klxx)
     &                           + bdens(iki)*
     &                             bfn_hess(ipt,mu,hcc(mc,kc))*
     &                             bfng_val(ipt,nu,lc)
     &                           + bdens(iki)*
     &                             bfng_val(ipt,mu,kc)*
     &                             bfn_hess(ipt,nu,hcc(mc,lc))
                              enddo ! ipt
                           enddo ! mc
                        enddo ! lc
                     enddo ! kc
                  enddo ! nu
                  do nu=mu+1, bend
                     iki = iky(nu)+mu
                     do kc = 1, 3
                        do lc = 1, kc
                           klxx = kc*(kc-1)/2+lc
                           do mc=1,3
                              do ipt=1,npts
                                 ddgrhoa(ipt,mc,klxx)
     &                           = ddgrhoa(ipt,mc,klxx)
     &                           + adens(iki)*
     &                             bfn_hess(ipt,mu,hcc(mc,kc))*
     &                             bfng_val(ipt,nu,lc)
     &                           + adens(iki)*
     &                             bfng_val(ipt,mu,kc)*
     &                             bfn_hess(ipt,nu,hcc(mc,lc))
                                 ddgrhob(ipt,mc,klxx)
     &                           = ddgrhob(ipt,mc,klxx)
     &                           + bdens(iki)*
     &                             bfn_hess(ipt,mu,hcc(mc,kc))*
     &                             bfng_val(ipt,nu,lc)
     &                           + bdens(iki)*
     &                             bfng_val(ipt,mu,kc)*
     &                             bfn_hess(ipt,nu,hcc(mc,lc))
                              enddo ! ipt
                           enddo ! mc
                        enddo ! lc
                     enddo ! kc
                  enddo ! nu
               enddo ! mu
            endif ! gradcorr_sw
c
            call aclear_dp(scr,6,0.0d0)
            do kc=1,3
               kx=kc*(kc-1)/2
               do lc=1,kc
                  klxx = kx+lc
                  do ipt=1,npts
                     scr(klxx)=scr(klxx)
     &               + wt(ipt)*xc_hpt(ipt,irara)*
     &                 drho(ipt,1,kpert+kc)*drho(ipt,1,kpert+lc)
     &               + wt(ipt)*xc_hpt(ipt,irbrb)*
     &                 drho(ipt,2,kpert+kc)*drho(ipt,2,kpert+lc)
     &               + wt(ipt)*xc_hpt(ipt,irarb)*
     &                 (drho(ipt,1,kpert+kc)*drho(ipt,2,kpert+lc)+
     &                  drho(ipt,2,kpert+kc)*drho(ipt,1,kpert+lc))
     &               + wt(ipt)*xc_vpt(ipt,ira)*ddrhoa(ipt,klxx)
     &               + wt(ipt)*xc_vpt(ipt,irb)*ddrhob(ipt,klxx)
                  enddo
                  if (gradcorr_sw) then
                     do ipt=1,npts
                        gaada 
     &                  = 2*(dgrho(ipt,1,1,kpert+kc)*grho(ipt,1,1)
     &                     + dgrho(ipt,1,2,kpert+kc)*grho(ipt,1,2)
     &                     + dgrho(ipt,1,3,kpert+kc)*grho(ipt,1,3))
                        gbbda 
     &                  = 2*(dgrho(ipt,2,1,kpert+kc)*grho(ipt,2,1)
     &                     + dgrho(ipt,2,2,kpert+kc)*grho(ipt,2,2)
     &                     + dgrho(ipt,2,3,kpert+kc)*grho(ipt,2,3))
                        gabda 
     &                  =    dgrho(ipt,1,1,kpert+kc)*grho(ipt,2,1)
     &                     + dgrho(ipt,1,2,kpert+kc)*grho(ipt,2,2)
     &                     + dgrho(ipt,1,3,kpert+kc)*grho(ipt,2,3)
     &                     + dgrho(ipt,2,1,kpert+kc)*grho(ipt,1,1)
     &                     + dgrho(ipt,2,2,kpert+kc)*grho(ipt,1,2)
     &                     + dgrho(ipt,2,3,kpert+kc)*grho(ipt,1,3)
                        gaadb 
     &                  = 2*(dgrho(ipt,1,1,kpert+lc)*grho(ipt,1,1)
     &                     + dgrho(ipt,1,2,kpert+lc)*grho(ipt,1,2)
     &                     + dgrho(ipt,1,3,kpert+lc)*grho(ipt,1,3))
                        gbbdb 
     &                  = 2*(dgrho(ipt,2,1,kpert+lc)*grho(ipt,2,1)
     &                     + dgrho(ipt,2,2,kpert+lc)*grho(ipt,2,2)
     &                     + dgrho(ipt,2,3,kpert+lc)*grho(ipt,2,3))
                        gabdb 
     &                  =    dgrho(ipt,1,1,kpert+lc)*grho(ipt,2,1)
     &                     + dgrho(ipt,1,2,kpert+lc)*grho(ipt,2,2)
     &                     + dgrho(ipt,1,3,kpert+lc)*grho(ipt,2,3)
     &                     + dgrho(ipt,2,1,kpert+lc)*grho(ipt,1,1)
     &                     + dgrho(ipt,2,2,kpert+lc)*grho(ipt,1,2)
     &                     + dgrho(ipt,2,3,kpert+lc)*grho(ipt,1,3)
                        gaadab 
     &                  = 2*(ddgrhoa(ipt,1,klxx)*grho(ipt,1,1)
     &                      +ddgrhoa(ipt,2,klxx)*grho(ipt,1,2)
     &                      +ddgrhoa(ipt,3,klxx)*grho(ipt,1,3)
     &                      +dgrho(ipt,1,1,kpert+kc)*
     &                       dgrho(ipt,1,1,kpert+lc)
     &                      +dgrho(ipt,1,2,kpert+kc)*
     &                       dgrho(ipt,1,2,kpert+lc)
     &                      +dgrho(ipt,1,3,kpert+kc)*
     &                       dgrho(ipt,1,3,kpert+lc))
                        gbbdab 
     &                  = 2*(ddgrhob(ipt,1,klxx)*grho(ipt,2,1)
     &                      +ddgrhob(ipt,2,klxx)*grho(ipt,2,2)
     &                      +ddgrhob(ipt,3,klxx)*grho(ipt,2,3)
     &                      +dgrho(ipt,2,1,kpert+kc)*
     &                       dgrho(ipt,2,1,kpert+lc)
     &                      +dgrho(ipt,2,2,kpert+kc)*
     &                       dgrho(ipt,2,2,kpert+lc)
     &                      +dgrho(ipt,2,3,kpert+kc)*
     &                       dgrho(ipt,2,3,kpert+lc))
                        gabdab 
     &                  =   (ddgrhoa(ipt,1,klxx)*grho(ipt,2,1)
     &                      +ddgrhoa(ipt,2,klxx)*grho(ipt,2,2)
     &                      +ddgrhoa(ipt,3,klxx)*grho(ipt,2,3)
     &                      +ddgrhob(ipt,1,klxx)*grho(ipt,1,1)
     &                      +ddgrhob(ipt,2,klxx)*grho(ipt,1,2)
     &                      +ddgrhob(ipt,3,klxx)*grho(ipt,1,3)
     &                      +dgrho(ipt,1,1,kpert+kc)*
     &                       dgrho(ipt,2,1,kpert+lc)
     &                      +dgrho(ipt,1,2,kpert+kc)*
     &                       dgrho(ipt,2,2,kpert+lc)
     &                      +dgrho(ipt,1,3,kpert+kc)*
     &                       dgrho(ipt,2,3,kpert+lc)
     &                      +dgrho(ipt,2,1,kpert+kc)*
     &                       dgrho(ipt,1,1,kpert+lc)
     &                      +dgrho(ipt,2,2,kpert+kc)*
     &                       dgrho(ipt,1,2,kpert+lc)
     &                      +dgrho(ipt,2,3,kpert+kc)*
     &                       dgrho(ipt,1,3,kpert+lc))
                        scr(klxx)=scr(klxx)
     &                  + wt(ipt)*xc_dhpt(ipt,iragaa)*
     &                    (drho(ipt,1,kpert+kc)*gaadb
     &                    +drho(ipt,1,kpert+lc)*gaada)
     &                  + wt(ipt)*xc_dhpt(ipt,iragab)*
     &                    (drho(ipt,1,kpert+kc)*gabdb
     &                    +drho(ipt,1,kpert+lc)*gabda)
     &                  + wt(ipt)*xc_dhpt(ipt,iragbb)*
     &                    (drho(ipt,1,kpert+kc)*gbbdb
     &                    +drho(ipt,1,kpert+lc)*gbbda)
     &                  + wt(ipt)*xc_dhpt(ipt,irbgaa)*
     &                    (drho(ipt,2,kpert+kc)*gaadb
     &                    +drho(ipt,2,kpert+lc)*gaada)
     &                  + wt(ipt)*xc_dhpt(ipt,irbgab)*
     &                    (drho(ipt,2,kpert+kc)*gabdb
     &                    +drho(ipt,2,kpert+lc)*gabda)
     &                  + wt(ipt)*xc_dhpt(ipt,irbgbb)*
     &                    (drho(ipt,2,kpert+kc)*gbbdb
     &                    +drho(ipt,2,kpert+lc)*gbbda)
     &                  + wt(ipt)*xc_dhpt(ipt,igaagaa)*
     &                    gaada*gaadb
     &                  + wt(ipt)*xc_dhpt(ipt,igaagab)*
     &                    gaada*gabdb
     &                  + wt(ipt)*xc_dhpt(ipt,igaagbb)*
     &                    gaada*gbbdb
     &                  + wt(ipt)*xc_dhpt(ipt,igaagab)*
     &                    gabda*gaadb
     &                  + wt(ipt)*xc_dhpt(ipt,igabgab)*
     &                    gabda*gabdb
     &                  + wt(ipt)*xc_dhpt(ipt,igabgbb)*
     &                    gabda*gbbdb
     &                  + wt(ipt)*xc_dhpt(ipt,igaagbb)*
     &                    gbbda*gaadb
     &                  + wt(ipt)*xc_dhpt(ipt,igabgbb)*
     &                    gbbda*gabdb
     &                  + wt(ipt)*xc_dhpt(ipt,igbbgbb)*
     &                    gbbda*gbbdb
     &                  + wt(ipt)*xc_dvpt(ipt,igaa)*gaadab
     &                  + wt(ipt)*xc_dvpt(ipt,igab)*gabdab
     &                  + wt(ipt)*xc_dvpt(ipt,igbb)*gbbdab
                     enddo
                  endif
               enddo
            enddo
            do kc=1,3
               kx=kc*(kc-1)/2
               do lc=1,kc-1
                  hess(kpert+kc,kpert+lc) = hess(kpert+kc,kpert+lc)
     &                                    + scr(kx+lc)
                  hess(kpert+lc,kpert+kc) = hess(kpert+lc,kpert+kc)
     &                                    + scr(kx+lc)
                  if (gwt_sw) then
                     hess(kpert+kc,ipert+lc) 
     &               = hess(kpert+kc,ipert+lc) - scr(kx+lc)
                     hess(ipert+lc,kpert+kc) 
     &               = hess(ipert+lc,kpert+kc) - scr(kx+lc)
c
                     hess(ipert+kc,kpert+lc) 
     &               = hess(ipert+kc,kpert+lc) - scr(kx+lc)
                     hess(kpert+lc,ipert+kc) 
     &               = hess(kpert+lc,ipert+kc) - scr(kx+lc)
c
                     hess(ipert+kc,ipert+lc) 
     &               = hess(ipert+kc,ipert+lc) + scr(kx+lc)
                     hess(ipert+lc,ipert+kc) 
     &               = hess(ipert+lc,ipert+kc) + scr(kx+lc)
                  endif
               enddo
               hess(kpert+kc,kpert+kc) = hess(kpert+kc,kpert+kc)
     &                                 + scr(kx+kc)
               if (gwt_sw) then
                  hess(kpert+kc,ipert+kc) = hess(kpert+kc,ipert+kc)
     &                                    - scr(kx+kc)
                  hess(ipert+kc,kpert+kc) = hess(ipert+kc,kpert+kc)
     &                                    - scr(kx+kc)
                  hess(ipert+kc,ipert+kc) = hess(ipert+kc,ipert+kc)
     &                                    + scr(kx+kc)
               endif
            enddo
 220     continue
      endif
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine hess_dft_ao_scr(rkstyp_sw,gradcorr_sw,gwt_sw,
     &     gwt_avail_sw,active_bfn_list,active_bfn_indx,active_bfn_atms,
     &     n_active_bfn,n_active_atm,
     &     near_atom_list,num_near_atoms,
     &     nao,npert,npts,
     &     wt,gwt,ggwt,grho,drho,dgrho,ddrhoa,ddrhob,ddgrhoa,ddgrhob,
     &     bfn_val,bfng_val,bfn_hess,bfn_3rd,adens,bdens,
     &     xc_ept,xc_vpt,xc_dvpt,xc_hpt,xc_dhpt,iatm,hess,
     &     mxp,dentol,istart_bfn)
      implicit none
c
c     Calculates the DFT contributions to the 2nd derivative of the 
c     energy. Only the explicit geometry dependencies are being 
c     differentiated.
c
c     Parameters:
c
C *Parameter module for CCP1/DFT            
C *Memory information		
C *----------------
C * max_block 	- maximum number of blocks which can be allocate
      integer max_block
      parameter(max_block=20)
C *							
C *Basis sets information                              
C *----------------------				
C *max_tag 	- maximum number of basis sets which can be alloc
C *max_atype	- maximum number of centre types	
C *max_gtype	- maximum number of grid types (at least every element 
C                 has different grid type)
C *max_grids    - maximum number of grids (different terms may have
C                 different grids, i.e. CPKS equations use a different
C                 grid than the one used for the KS-matrix)
C *max_shel	- maximum number of shells on centre
C *max_prm	- maximum number of primitives for any given centre
C *maxL		- maximum angular momentum allowed	
C *max_func	- maximum number of basis functions for any centre
      integer max_tag,max_atype,max_shel,max_prm,max_ang
      integer max_gtype,max_grids
      parameter(max_tag=3,max_atype=30,max_shel=500,max_prm=5000)
      parameter(max_ang=5,max_gtype=10,max_grids=2)
C *								
C *Geometry information					
C *--------------------					
C *max_atom	- maximum number of atoms in system
      integer max_atom
      parameter(max_atom=750)
C *
C *Accuracy information				
C ---------------------			
C *global_accuracy 	- global accuracy	
      real*8  global_accuracy
      parameter(global_accuracy=1.0d-14)

C *Grid information
c *----------------
c *maxradzn - The maximum number of radial zones. 
      integer maxgpt,maxrad,maxfpt,maxang,maxradzn,maxtablerows
      parameter(maxgpt=2900,maxrad=50,maxang=302,maxfpt=100)
      parameter(maxradzn=35,maxtablerows=7)
c
c     The parameters defined here label the 1st, 2nd and 3rd derivative
c     components of the basis functions.
c
c     2nd derivatives:
c
      integer hxx, hyy, hzz, hxy, hxz, hyz
c
      parameter (hxx = 1)
      parameter (hyy = 2)
      parameter (hzz = 3)
      parameter (hxy = 4)
      parameter (hxz = 5)
      parameter (hyz = 6)
c
      integer hcc(3,3)
      data hcc/hxx,hxy,hxz,
     &         hxy,hyy,hyz,
     &         hxz,hyz,hzz/
c
c     3rd derivatives:
c
      integer txxx, txxy, txxz, txyy, txyz
      integer txzz, tyyy, tyyz, tyzz, tzzz
c
      parameter (txxx =  1)
      parameter (txxy =  2)
      parameter (txxz =  3)
      parameter (txyy =  4)
      parameter (txyz =  5)
      parameter (txzz =  6)
      parameter (tyyy =  7)
      parameter (tyyz =  8)
      parameter (tyzz =  9)
      parameter (tzzz = 10)
c
      integer tc3(3,3,3)
      data tc3/txxx,txxy,txxz,
     &         txxy,txyy,txyz,
     &         txxz,txyz,txzz,
     &         txxy,txyy,txyz,
     &         txyy,tyyy,tyyz,
     &         txyz,tyyz,tyzz,
     &         txxz,txyz,txzz,
     &         txyz,tyyz,tyzz,
     &         txzz,tyzz,tzzz/
c
c     The parameters defined here label the 1st and 2nd derivative
c     components of the density functional. Because a functional
c     may be purely local or of a GGA type the components are grouped
c     into derivative with respect to the density only and others.
c
c     1st derivatives:
c
      integer ira, irb
      integer igaa, igab, igbb
      integer ita, itb
c
c     2nd derivatives:
c
      integer irara, irarb, irbrb
      integer iragaa, iragab, iragbb, irbgaa, irbgab, irbgbb
      integer igaagaa, igaagab, igaagbb, igabgab, igabgbb, igbbgbb
      integer irata, iratb, irbta, irbtb
      integer igaata, igaatb, igabta, igabtb, igbbta, igbbtb
      integer itata, itatb, itbtb
c
      parameter (ira     = 1)
      parameter (irb     = 2)
c
      parameter (igaa    = 1)
      parameter (igab    = 2)
      parameter (igbb    = 3)
c
      parameter (ita     = 1)
      parameter (itb     = 2)
c
c
      parameter (irara   = 1)
      parameter (irarb   = 2)
      parameter (irbrb   = 3)
c
      parameter (iragaa  = 1)
      parameter (iragab  = 3)
      parameter (iragbb  = 4)
      parameter (irbgaa  = 5)
      parameter (irbgab  = 6)
      parameter (irbgbb  = 7)
      parameter (igaagaa = 2)
      parameter (igaagab = 8)
      parameter (igaagbb = 9)
      parameter (igabgab = 10)
      parameter (igabgbb = 11)
      parameter (igbbgbb = 12)
c
      parameter (irata   = 1)
      parameter (iratb   = 3)
      parameter (irbta   = 4)
      parameter (irbtb   = 5)
      parameter (itata   = 2)
      parameter (itatb   = 6)
      parameter (itbtb   = 7)
c
      parameter (igaata  = 1)
      parameter (igaatb  = 2)
      parameter (igabta  = 3)
      parameter (igabtb  = 4)
      parameter (igbbta  = 5)
      parameter (igbbtb  = 6)
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
c     In variables:
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      logical rkstyp_sw
      logical gradcorr_sw
      logical gwt_sw
      logical gwt_avail_sw
      integer n_active_bfn
      integer n_active_atm
      integer active_bfn_list(n_active_bfn)
      integer active_bfn_indx(n_active_bfn)
      integer active_bfn_atms(n_active_atm)
      integer num_near_atoms
      integer near_atom_list(num_near_atoms)
      integer nao
      integer npert
      integer npts
      integer iatm
      integer mxp
      real*8 wt(mxp)
      real*8 gwt(3,mxp,num_near_atoms)
      real*8 ggwt(mxp,3*num_near_atoms,3*num_near_atoms)
      real*8 grho(mxp,2,3)
      real*8 bfn_val(mxp,nao)
      real*8 bfng_val(mxp,nao,3)
      real*8 bfn_hess(mxp,nao,6)
      real*8 bfn_3rd(mxp,nao,10)
      real*8 adens(nao*(nao+1)/2)
      real*8 bdens(nao*(nao+1)/2)
      real*8 dentol
c
c     explicit derivatives only (((
c
      real*8 drho(mxp,2,3*n_active_atm)
      real*8 dgrho(mxp,2,3,3*n_active_atm)
      real*8 ddrhoa(npts,9)
      real*8 ddrhob(npts,9)
      real*8 ddgrhoa(npts,3,9)
      real*8 ddgrhob(npts,3,9)
c
c     ))) explicit derivatives only 
c
      real*8 xc_ept(mxp)
      real*8 xc_vpt(mxp,2)
      real*8 xc_dvpt(mxp,3)
      real*8 xc_hpt(mxp,3)
      real*8 xc_dhpt(mxp,12)
c
c     Out variables
c
      real*8 hess(npert,npert)
c
c     Local variables
c
      integer katm, latm, kpert, lpert, ipert, kprt, lprt
      integer kc, kx, lc, m, inu, imu, nu, mu, iki
      integer klxx, ipt
      real*8 scr(9)
      real*8 gaada, gaadb
      real*8 gabda, gabdb
      real*8 gbbda, gbbdb
      real*8 gaadab, gabdab, gbbdab
      real*8 dd
c
c     istart_bfn: For every atom where the basis functions start.
c                 Assumption is that all basis functions associated
c                 with 1 atom are consecutively located on the
c                 active_bfn_list.
c                 Beware: not every active atom necessarily has any
c                 active basis functions. (See below also.)
c
      integer istart_bfn(n_active_atm+1)
      real*8 a
c
c     Code:
c
c     if (gwt_sw) then
c     if (n_active_atm.ne.num_near_atom) then
c        call caserr('hess_dft_scr: n_active_atm.ne.num_near_atom')
c     endif
c     do lx=1,num_near_atom
c        if (near_atom_list(lx).ne.active_bfn_atms(lx) then
c        call caserr('hess_dft_scr: near_atom_list.ne.active_bfn_atms')
c        endif
c     enddo
c     endif
c
      ipert = 3*(iatm-1)
c
      do m=1,n_active_atm+1
         istart_bfn(m) = 0
      enddo
      do m=n_active_bfn,1,-1
         istart_bfn(active_bfn_indx(m))=m
      enddo
      istart_bfn(n_active_atm+1)=n_active_bfn+1
c
c     BEWARE: The list of active_bfn_atms is set up based on the radius
c     of the current shell of grid points. The list of active_bfn_list
c     is setup based on the grid points in the current batch. The
c     latter criterion is more selective!!!
c     THEREFORE: The active_bfn_atms may contain atoms that have NO
c     active basis functions for the current batch.
c
      do m=n_active_atm,1,-1
         if (istart_bfn(m).eq.0) istart_bfn(m) = istart_bfn(m+1)
      enddo
c
      if (rkstyp_sw) then
c***  BEGIN OF CLOSED SHELL
         if (gwt_sw.and.gwt_avail_sw) then
cc
cc           First the terms involving the 2nd derivative of the weights
cc
c            do katm = 1, num_near_atoms
c               kprt = 3*(katm-1)
c               kpert = 3*(near_atom_list(katm)-1)
c               do latm = 1, num_near_atoms
c                  lprt = 3*(latm-1)
c                  lpert = 3*(near_atom_list(latm)-1)
c                  call aclear_dp(scr,9,0.0d0)
c                  do kc=1,3
c                     kx=3*(kc-1)
c                     do lc=1,3
c                        do ipt=1,npts
c                           scr(kx+lc)=scr(kx+lc)
c     &                     + ggwt(ipt,kprt+kc,lprt+lc)
c     &                       *xc_ept(ipt)
c                        enddo
c                     enddo
c                  enddo
c                  do kc=1,3
c                     kx=3*(kc-1)
c                     do lc=1,3
c                        hess(kpert+kc,lpert+lc) 
c     &                  = hess(kpert+kc,lpert+lc) + scr(kx+lc)
c                     enddo
c                  enddo
c               enddo
c            enddo
c
c           Now the terms involving the 1st derivative of the weights
c
            do 100 katm = 1, n_active_atm
               kpert = 3*(active_bfn_atms(katm)-1)
               if (kpert.eq.ipert) goto 100
               kprt = 3*(katm-1)
               do latm = 1, num_near_atoms
                  call aclear_dp(scr,9,0.0d0)
                  do kc = 1, 3
                     kx = 3*(kc-1)
                     do lc = 1, 3
                        do ipt = 1, npts
                           scr(kx+lc) = scr(kx+lc)
     &                     + 2*gwt(lc,ipt,latm)*
     &                       xc_vpt(ipt,1)*drho(ipt,1,kprt+kc)
                        enddo
                     enddo
                  enddo
                  if (gradcorr_sw) then
                     do kc = 1, 3
                        kx = 3*(kc-1)
                        do lc = 1, 3
                           do ipt = 1, npts
                              gaada
     &                        = (dgrho(ipt,1,1,kprt+kc)*grho(ipt,1,1)
     &                         + dgrho(ipt,1,2,kprt+kc)*grho(ipt,1,2)
     &                         + dgrho(ipt,1,3,kprt+kc)*grho(ipt,1,3))
                              scr(kx+lc) = scr(kx+lc)
     &                        +gwt(lc,ipt,latm)*xc_dvpt(ipt,igaa)*gaada
                           enddo
                        enddo
                     enddo
                  endif
                  lpert = 3*(near_atom_list(latm)-1)
                  do kc = 1, 3
                     kx = 3*(kc-1)
                     do lc = 1, 3
                        hess(kpert+kc,lpert+lc)
     &                  = hess(kpert+kc,lpert+lc) + scr(kx+lc)
                        hess(ipert+kc,lpert+lc)
     &                  = hess(ipert+kc,lpert+lc) - scr(kx+lc)
                     enddo
                  enddo
               enddo
 100        continue
            do katm = 1, num_near_atoms
               kprt = 3*(katm-1)
               kpert = 3*(near_atom_list(katm)-1)
               do 110 latm = 1, n_active_atm
                  lpert = 3*(active_bfn_atms(latm)-1)
                  if (lpert.eq.ipert) goto 110
                  lprt = 3*(latm-1)
                  call aclear_dp(scr,9,0.0d0)
                  do kc = 1, 3
                     kx = 3*(kc-1)
                     do lc = 1, 3
                        do ipt = 1, npts
                           scr(kx+lc) = scr(kx+lc)
     &                     + 2*gwt(kc,ipt,katm)*
     &                       xc_vpt(ipt,1)*drho(ipt,1,lprt+lc)
                        enddo
                     enddo
                  enddo
                  if (gradcorr_sw) then
                     do kc = 1, 3
                        kx = 3*(kc-1)
                        do lc = 1, 3
                           do ipt = 1, npts
                              gaada
     &                        = (dgrho(ipt,1,1,lprt+lc)*grho(ipt,1,1)
     &                         + dgrho(ipt,1,2,lprt+lc)*grho(ipt,1,2)
     &                         + dgrho(ipt,1,3,lprt+lc)*grho(ipt,1,3))
                              scr(kx+lc) = scr(kx+lc)
     &                        +gwt(kc,ipt,katm)*xc_dvpt(ipt,igaa)*gaada
                           enddo
                        enddo
                     enddo
                  endif
                  do kc = 1, 3
                     kx = 3*(kc-1)
                     do lc = 1, 3
                        hess(kpert+kc,lpert+lc)
     &                  = hess(kpert+kc,lpert+lc) + scr(kx+lc)
                        hess(kpert+kc,ipert+lc)
     &                  = hess(kpert+kc,ipert+lc) - scr(kx+lc)
                     enddo
                  enddo
 110           continue
            enddo
         endif

         do 120 katm = 1, n_active_atm
            kpert = 3*(active_bfn_atms(katm)-1)
            if (gwt_sw.and.kpert.eq.ipert) goto 120
            kprt = 3*(katm-1)
            do 130 latm = 1, katm-1
               lpert = 3*(active_bfn_atms(latm)-1)
               if (gwt_sw.and.lpert.eq.ipert) goto 130
               lprt = 3*(latm-1)
c
               call aclear_dp(ddrhoa,npts*9,0.0d0)
               do imu=istart_bfn(katm),istart_bfn(katm+1)-1
                  mu=active_bfn_list(imu)
                  do inu=istart_bfn(latm),istart_bfn(latm+1)-1
                     nu=active_bfn_list(inu)
                     iki = iky(mu)+nu
                     dd=adens(iki)
                     if (abs(dd).gt.dentol) then
                        do kc = 1, 3
                           kx = 3*(kc-1)
                           do lc = 1, 3
                              klxx = kx+lc
                              do ipt=1,npts
                                 ddrhoa(ipt,klxx) = ddrhoa(ipt,klxx)
     &                           + dd*bfng_val(ipt,imu,kc)*
     &                             bfng_val(ipt,inu,lc)
                              enddo ! ipt
                           enddo ! lc
                        enddo ! kc
                     endif
                  enddo ! nu
               enddo ! mu
c
               if (gradcorr_sw) then
                  call aclear_dp(ddgrhoa,npts*3*9,0.0d0)
                  do imu=istart_bfn(katm),istart_bfn(katm+1)-1
                     mu=active_bfn_list(imu)
                     do inu=istart_bfn(latm),istart_bfn(latm+1)-1
                        nu=active_bfn_list(inu)
                        iki = iky(mu)+nu
                        dd=adens(iki)
                        if (abs(dd).gt.dentol) then
                           do kc = 1, 3
                              kx = 3*(kc-1)
                              do lc = 1, 3
                                 klxx = kx+lc
                                 do m=1,3
                                    do ipt=1,npts
                                       ddgrhoa(ipt,m,klxx)
     &                                 = ddgrhoa(ipt,m,klxx)
     &                                 + dd*
     &                                   bfn_hess(ipt,imu,hcc(m,kc))*
     &                                   bfng_val(ipt,inu,lc)
     &                                 + dd*
     &                                   bfng_val(ipt,imu,kc)*
     &                                   bfn_hess(ipt,inu,hcc(m,lc))
                                    enddo ! ipt
                                 enddo ! m
                              enddo ! lc
                           enddo ! kc
                        endif
                     enddo ! nu
                  enddo ! mu
               endif ! gradcorr_sw
c
               call aclear_dp(scr,9,0.0d0)
               do kc=1,3
                  kx=3*(kc-1)
                  do lc=1,3
                     klxx = kx+lc
                     do ipt=1,npts
                        scr(klxx)=scr(klxx)
     &                  + 2*wt(ipt)*xc_hpt(ipt,irara)*
     &                    drho(ipt,1,kprt+kc)*drho(ipt,1,lprt+lc)
     &                  + 2*wt(ipt)*xc_vpt(ipt,ira)*ddrhoa(ipt,klxx)
                     enddo
                     if (gradcorr_sw) then
                        do ipt=1,npts
                           gaada 
     &                     =   (dgrho(ipt,1,1,kprt+kc)*grho(ipt,1,1)
     &                        + dgrho(ipt,1,2,kprt+kc)*grho(ipt,1,2)
     &                        + dgrho(ipt,1,3,kprt+kc)*grho(ipt,1,3))
                           gaadb 
     &                     =   (dgrho(ipt,1,1,lprt+lc)*grho(ipt,1,1)
     &                        + dgrho(ipt,1,2,lprt+lc)*grho(ipt,1,2)
     &                        + dgrho(ipt,1,3,lprt+lc)*grho(ipt,1,3))
                           gaadab 
     &                     =   (ddgrhoa(ipt,1,klxx)*grho(ipt,1,1)
     &                         +ddgrhoa(ipt,2,klxx)*grho(ipt,1,2)
     &                         +ddgrhoa(ipt,3,klxx)*grho(ipt,1,3))
     &                     + 2*(dgrho(ipt,1,1,kprt+kc)*
     &                          dgrho(ipt,1,1,lprt+lc)
     &                         +dgrho(ipt,1,2,kprt+kc)*
     &                          dgrho(ipt,1,2,lprt+lc)
     &                         +dgrho(ipt,1,3,kprt+kc)*
     &                          dgrho(ipt,1,3,lprt+lc))
                           scr(klxx)=scr(klxx)
     &                     + 2*wt(ipt)*xc_dhpt(ipt,iragaa)*
     &                       (drho(ipt,1,kprt+kc)*gaadb
     &                       +drho(ipt,1,lprt+lc)*gaada)
     &                     + wt(ipt)*xc_dhpt(ipt,igaagaa)*
     &                       gaada*gaadb
     &                     + wt(ipt)*xc_dvpt(ipt,igaa)*gaadab
                        enddo
                     endif
                  enddo
               enddo
               do kc=1,3
                  kx=3*(kc-1)
                  do lc=1,3
                     hess(kpert+kc,lpert+lc) = hess(kpert+kc,lpert+lc)
     &                                       + scr(kx+lc)
                     hess(lpert+lc,kpert+kc) = hess(lpert+lc,kpert+kc)
     &                                       + scr(kx+lc)
                     if (gwt_sw) then
                        hess(kpert+kc,ipert+lc) 
     &                  = hess(kpert+kc,ipert+lc) - scr(kx+lc)
                        hess(ipert+lc,kpert+kc) 
     &                  = hess(ipert+lc,kpert+kc) - scr(kx+lc)
c
                        hess(ipert+kc,lpert+lc) 
     &                  = hess(ipert+kc,lpert+lc) - scr(kx+lc)
                        hess(lpert+lc,ipert+kc) 
     &                  = hess(lpert+lc,ipert+kc) - scr(kx+lc)
c
                        hess(ipert+kc,ipert+lc) 
     &                  = hess(ipert+kc,ipert+lc) + scr(kx+lc)
                        hess(ipert+lc,ipert+kc) 
     &                  = hess(ipert+lc,ipert+kc) + scr(kx+lc)
                     endif
                  enddo
               enddo
 130        continue
c
            call aclear_dp(ddrhoa,npts*6,0.0d0)
c           2nd derivatives of 1 basis function
            do imu=istart_bfn(katm),istart_bfn(katm+1)-1
               mu=active_bfn_list(imu)
               iki=iky(mu)
               do inu=1,imu
                  nu=active_bfn_list(inu)
                  dd=adens(iki+nu)
                  if (abs(dd).gt.dentol) then
                     do kc=1,3
                        kx=kc*(kc-1)/2
                        do lc=1,kc
                           klxx = kx+lc
                           do ipt=1,npts
                              ddrhoa(ipt,klxx) = 
     &                        ddrhoa(ipt,klxx)
     &                        + dd*bfn_hess(ipt,imu,hcc(lc,kc))*
     &                          bfn_val(ipt,inu)
                           enddo ! ipt
                        enddo ! lc
                     enddo ! kc
                  endif
               enddo ! nu
               do inu=imu+1, n_active_bfn
                  nu=active_bfn_list(inu)
                  iki=iky(nu)
                  dd=adens(iki+mu)
                  if (abs(dd).gt.dentol) then
                     do kc=1,3
                        kx=kc*(kc-1)/2
                        do lc=1,kc
                           klxx = kx+lc
                           do ipt=1,npts
                              ddrhoa(ipt,klxx) = 
     &                        ddrhoa(ipt,klxx)
     &                        + dd*bfn_hess(ipt,imu,hcc(lc,kc))*
     &                          bfn_val(ipt,inu)
                           enddo ! ipt
                        enddo ! lc
                     enddo ! kc
                  endif
               enddo ! nu
            enddo ! mu
c           1st derivative of 2 basis functions
            do imu=istart_bfn(katm),istart_bfn(katm+1)-1
               mu=active_bfn_list(imu)
               do inu=istart_bfn(katm),imu
                  nu=active_bfn_list(inu)
                  iki = iky(mu)+nu
                  dd=adens(iki)
                  if (abs(dd).gt.dentol) then
                     do kc = 1, 3
                        kx = kc*(kc-1)/2
                        do lc = 1, kc
                           klxx = kx+lc
                           do ipt=1,npts
                              ddrhoa(ipt,klxx) = ddrhoa(ipt,klxx)
     &                        + dd*bfng_val(ipt,imu,kc)*
     &                          bfng_val(ipt,inu,lc)
                           enddo ! ipt
                        enddo ! lc
                     enddo ! kc
                  endif
               enddo ! nu
               do inu=imu+1,istart_bfn(katm+1)-1
                  nu=active_bfn_list(inu)
                  iki = iky(nu)+mu
                  dd=adens(iki)
                  if (abs(dd).gt.dentol) then
                     do kc = 1, 3
                        kx = kc*(kc-1)/2
                        do lc = 1, kc
                           klxx = kx+lc
                           do ipt=1,npts
                              ddrhoa(ipt,klxx) = ddrhoa(ipt,klxx)
     &                        + dd*bfng_val(ipt,imu,kc)*
     &                          bfng_val(ipt,inu,lc)
                           enddo ! ipt
                        enddo ! lc
                     enddo ! kc
                  endif
               enddo ! nu
            enddo ! mu
c
            if (gradcorr_sw) then
               call aclear_dp(ddgrhoa,npts*3*9,0.0d0)
c              2nd derivative of 1 basis function
               do imu=istart_bfn(katm),istart_bfn(katm+1)-1
                  mu=active_bfn_list(imu)
                  iki=iky(mu)
                  do inu=1, imu
                     nu=active_bfn_list(inu)
                     dd=adens(iki+nu)
                     if (abs(dd).gt.dentol) then
                        do kc=1,3
                           kx=kc*(kc-1)/2
                           do lc=1,kc
                              klxx=kx+lc
                              do m=1,3
                                 do ipt=1,npts
                                    ddgrhoa(ipt,m,klxx)
     &                              = ddgrhoa(ipt,m,klxx)
     &                              + dd*(
     &                                 bfn_3rd(ipt,imu,tc3(m,lc,kc))
     &                                *bfn_val(ipt,inu)+
     &                                 bfn_hess(ipt,imu,hcc(lc,kc))
     &                                *bfng_val(ipt,inu,m))
                                 enddo ! ipt
                              enddo ! m
                           enddo ! lc
                        enddo ! kc
                     endif
                  enddo ! inu
                  do inu=imu+1, n_active_bfn
                     nu=active_bfn_list(inu)
                     iki=iky(nu)
                     dd=adens(iki+mu)
                     if (abs(dd).gt.dentol) then
                        do kc=1,3
                           kx=kc*(kc-1)/2
                           do lc=1,kc
                              klxx=kx+lc
                              do m=1,3
                                 do ipt=1,npts
                                    ddgrhoa(ipt,m,klxx)
     &                              = ddgrhoa(ipt,m,klxx)
     &                              + dd*(
     &                                 bfn_3rd(ipt,imu,tc3(m,lc,kc))
     &                                *bfn_val(ipt,inu)+
     &                                 bfn_hess(ipt,imu,hcc(lc,kc))
     &                                *bfng_val(ipt,inu,m))
                                 enddo ! ipt
                              enddo ! m
                           enddo ! lc
                        enddo ! kc
                     endif
                  enddo ! nu
               enddo ! mu
c              1st derivative of 2 basis functions
               do imu=istart_bfn(katm),istart_bfn(katm+1)-1
                  mu=active_bfn_list(imu)
                  do inu=istart_bfn(katm),imu
                     nu=active_bfn_list(inu)
                     iki = iky(mu)+nu
                     dd=adens(iki)
                     if (abs(dd).gt.dentol) then
                        do kc = 1, 3
                           kx = kc*(kc-1)/2
                           do lc = 1, kc
                              klxx = kx+lc
                              do m=1,3
                                 do ipt=1,npts
                                    ddgrhoa(ipt,m,klxx)
     &                              = ddgrhoa(ipt,m,klxx)
     &                              + dd*
     &                                bfn_hess(ipt,imu,hcc(m,kc))*
     &                                bfng_val(ipt,inu,lc)
     &                              + dd*
     &                                bfng_val(ipt,imu,kc)*
     &                                bfn_hess(ipt,inu,hcc(m,lc))
                                 enddo ! ipt
                              enddo ! m
                           enddo ! lc
                        enddo ! kc
                     endif
                  enddo ! nu
                  do inu=imu+1,istart_bfn(katm+1)-1
                     nu=active_bfn_list(inu)
                     iki = iky(nu)+mu
                     dd=adens(iki)
                     if (abs(dd).gt.dentol) then
                        do kc = 1, 3
                           kx = kc*(kc-1)/2
                           do lc = 1, kc
                              klxx = kx+lc
                              do m=1,3
                                 do ipt=1,npts
                                    ddgrhoa(ipt,m,klxx)
     &                              = ddgrhoa(ipt,m,klxx)
     &                              + dd*
     &                                bfn_hess(ipt,imu,hcc(m,kc))*
     &                                bfng_val(ipt,inu,lc)
     &                              + dd*
     &                                bfng_val(ipt,imu,kc)*
     &                                bfn_hess(ipt,inu,hcc(m,lc))
                                 enddo ! ipt
                              enddo ! m
                           enddo ! lc
                        enddo ! kc
                     endif
                  enddo ! nu
               enddo ! mu
            endif
c
            call aclear_dp(scr,6,0.0d0)
            do kc=1,3
               kx=kc*(kc-1)/2
               do lc=1,kc
                  klxx = kx+lc
                  do ipt=1,npts
                     scr(klxx)=scr(klxx)
     &               + 2*wt(ipt)*xc_hpt(ipt,irara)*
     &                 drho(ipt,1,kprt+kc)*drho(ipt,1,kprt+lc)
     &               + 2*wt(ipt)*xc_vpt(ipt,ira)*ddrhoa(ipt,klxx)
                  enddo
                  if (gradcorr_sw) then
                     do ipt=1,npts
                        gaada 
     &                  =   (dgrho(ipt,1,1,kprt+kc)*grho(ipt,1,1)
     &                     + dgrho(ipt,1,2,kprt+kc)*grho(ipt,1,2)
     &                     + dgrho(ipt,1,3,kprt+kc)*grho(ipt,1,3))
                        gaadb 
     &                  =   (dgrho(ipt,1,1,kprt+lc)*grho(ipt,1,1)
     &                     + dgrho(ipt,1,2,kprt+lc)*grho(ipt,1,2)
     &                     + dgrho(ipt,1,3,kprt+lc)*grho(ipt,1,3))
                        gaadab 
     &                  =   (ddgrhoa(ipt,1,klxx)*grho(ipt,1,1)
     &                      +ddgrhoa(ipt,2,klxx)*grho(ipt,1,2)
     &                      +ddgrhoa(ipt,3,klxx)*grho(ipt,1,3))
     &                  + 2*(dgrho(ipt,1,1,kprt+kc)*
     &                       dgrho(ipt,1,1,kprt+lc)
     &                      +dgrho(ipt,1,2,kprt+kc)*
     &                       dgrho(ipt,1,2,kprt+lc)
     &                      +dgrho(ipt,1,3,kprt+kc)*
     &                       dgrho(ipt,1,3,kprt+lc))
                        scr(klxx)=scr(klxx)
     &                  + 2*wt(ipt)*xc_dhpt(ipt,iragaa)*
     &                    (drho(ipt,1,kprt+kc)*gaadb
     &                    +drho(ipt,1,kprt+lc)*gaada)
     &                  + wt(ipt)*xc_dhpt(ipt,igaagaa)*
     &                    gaada*gaadb
     &                  + wt(ipt)*xc_dvpt(ipt,igaa)*gaadab
                     enddo
                  endif
               enddo
            enddo
            do kc=1,3
               kx=kc*(kc-1)/2
               do lc=1,kc-1
                  hess(kpert+kc,kpert+lc) = hess(kpert+kc,kpert+lc)
     &                                    + scr(kx+lc)
                  hess(kpert+lc,kpert+kc) = hess(kpert+lc,kpert+kc)
     &                                    + scr(kx+lc)
                  if (gwt_sw) then
                     hess(kpert+kc,ipert+lc) 
     &               = hess(kpert+kc,ipert+lc) - scr(kx+lc)
                     hess(ipert+lc,kpert+kc) 
     &               = hess(ipert+lc,kpert+kc) - scr(kx+lc)
c
                     hess(ipert+kc,kpert+lc) 
     &               = hess(ipert+kc,kpert+lc) - scr(kx+lc)
                     hess(kpert+lc,ipert+kc) 
     &               = hess(kpert+lc,ipert+kc) - scr(kx+lc)
c
                     hess(ipert+kc,ipert+lc) 
     &               = hess(ipert+kc,ipert+lc) + scr(kx+lc)
                     hess(ipert+lc,ipert+kc) 
     &               = hess(ipert+lc,ipert+kc) + scr(kx+lc)
                  endif
               enddo
               hess(kpert+kc,kpert+kc) = hess(kpert+kc,kpert+kc)
     &                                 + scr(kx+kc)
               if (gwt_sw) then
                  hess(kpert+kc,ipert+kc) = hess(kpert+kc,ipert+kc)
     &                                    - scr(kx+kc)
                  hess(ipert+kc,kpert+kc) = hess(ipert+kc,kpert+kc)
     &                                    - scr(kx+kc)
                  hess(ipert+kc,ipert+kc) = hess(ipert+kc,ipert+kc)
     &                                    + scr(kx+kc)
               endif
            enddo
 120     continue

c***  END OF CLOSED SHELL
      else ! rkstyp_sw
         if (gwt_sw.and.gwt_avail_sw) then
cc
cc           First the terms involving the 2nd derivative of the weights
cc
c            do katm = 1, num_near_atoms
c               kprt = 3*(katm-1)
c               kpert = 3*(near_atom_list(katm)-1)
c               do latm = 1, num_near_atoms
c                  lprt = 3*(latm-1)
c                  lpert = 3*(near_atom_list(latm)-1)
c                  call aclear_dp(scr,9,0.0d0)
c                  do kc=1,3
c                     kx=3*(kc-1)
c                     do lc=1,3
c                        do ipt=1,npts
c                           scr(kx+lc)=scr(kx+lc)
c     &                     + ggwt(ipt,kprt+kc,lprt+lc)
c     &                       *xc_ept(ipt)
c                        enddo
c                     enddo
c                  enddo
c                  do kc=1,3
c                     kx=3*(kc-1)
c                     do lc=1,3
c                        hess(kpert+kc,lpert+lc) 
c     &                  = hess(kpert+kc,lpert+lc) + scr(kx+lc)
c                     enddo
c                  enddo
c               enddo
c            enddo
c
c           Now the terms involving the 1st derivative of the weights
c
            do 200 katm = 1, n_active_atm
               kpert = 3*(active_bfn_atms(katm)-1)
               if (kpert.eq.ipert) goto 200
               kprt = 3*(katm-1)
               do latm = 1, num_near_atoms
                  call aclear_dp(scr,9,0.0d0)
                  do kc = 1, 3
                     kx = 3*(kc-1)
                     do lc = 1, 3
                        do ipt = 1, npts
                           scr(kx+lc) = scr(kx+lc)
     &                     + gwt(lc,ipt,latm)*
     &                       (xc_vpt(ipt,1)*drho(ipt,1,kprt+kc)
     &                       +xc_vpt(ipt,2)*drho(ipt,2,kprt+kc))
                        enddo
                     enddo
                  enddo
c                 lprt = 3*(latm-1)
                  lpert = 3*(near_atom_list(latm)-1)
                  do kc = 1, 3
                     kx = 3*(kc-1)
                     do lc = 1, 3
                        hess(kpert+kc,lpert+lc)
     &                  = hess(kpert+kc,lpert+lc) + scr(kx+lc)
                        hess(ipert+kc,lpert+lc)
     &                  = hess(ipert+kc,lpert+lc) - scr(kx+lc)
                     enddo
                  enddo
               enddo
 200        continue
            do katm = 1, num_near_atoms
c              kprt = 3*(katm-1)
               kpert = 3*(near_atom_list(katm)-1)
               do 210 latm = 1, n_active_atm
                  lpert = 3*(active_bfn_atms(latm)-1)
                  if (lpert.eq.ipert) goto 210
                  lprt = 3*(latm-1)
                  call aclear_dp(scr,9,0.0d0)
                  do kc = 1, 3
                     kx = 3*(kc-1)
                     do lc = 1, 3
                        do ipt = 1, npts
                           scr(kx+lc) = scr(kx+lc)
     &                     + gwt(kc,ipt,katm)*
     &                       (xc_vpt(ipt,1)*drho(ipt,1,lprt+lc)
     &                       +xc_vpt(ipt,2)*drho(ipt,2,lprt+lc))
                        enddo
                     enddo
                  enddo
                  do kc = 1, 3
                     kx = 3*(kc-1)
                     do lc = 1, 3
                        hess(kpert+kc,lpert+lc)
     &                  = hess(kpert+kc,lpert+lc) + scr(kx+lc)
                        hess(kpert+kc,ipert+lc)
     &                  = hess(kpert+kc,ipert+lc) - scr(kx+lc)
                     enddo
                  enddo
 210           continue
            enddo
         endif

         do 220 katm = 1, n_active_atm
            kpert = 3*(active_bfn_atms(katm)-1)
            if (gwt_sw.and.kpert.eq.ipert) goto 220
            kprt = 3*(katm-1)
            do 230 latm = 1, katm-1
               lpert = 3*(active_bfn_atms(latm)-1)
               if (gwt_sw.and.lpert.eq.ipert) goto 230
               lprt = 3*(latm-1)
c
               call aclear_dp(ddrhoa,npts*9,0.0d0)
               call aclear_dp(ddrhob,npts*9,0.0d0)
               do imu=istart_bfn(katm),istart_bfn(katm+1)-1
                  mu=active_bfn_list(imu)
                  do inu=istart_bfn(latm),istart_bfn(latm+1)-1
                     nu=active_bfn_list(inu)
                     iki = iky(mu)+nu
                     dd=adens(iki)
                     if (abs(dd).gt.dentol) then
                        do kc = 1, 3
                           kx = 3*(kc-1)
                           do lc = 1, 3
                              klxx = kx+lc
                              do ipt=1,npts
                                 ddrhoa(ipt,klxx) = ddrhoa(ipt,klxx)
     &                           + dd*bfng_val(ipt,imu,kc)*
     &                             bfng_val(ipt,inu,lc)
                              enddo ! ipt
                           enddo ! lc
                        enddo ! kc
                     endif
                     dd=bdens(iki)
                     if (abs(dd).gt.dentol) then
                        do kc = 1, 3
                           kx = 3*(kc-1)
                           do lc = 1, 3
                              klxx = kx+lc
                              do ipt=1,npts
                                 ddrhob(ipt,klxx) = ddrhob(ipt,klxx)
     &                           + dd*bfng_val(ipt,imu,kc)*
     &                             bfng_val(ipt,inu,lc)
                              enddo ! ipt
                           enddo ! lc
                        enddo ! kc
                     endif
                  enddo ! nu
               enddo ! mu
c
               if (gradcorr_sw) then
                  call aclear_dp(ddgrhoa,npts*3*9,0.0d0)
                  call aclear_dp(ddgrhob,npts*3*9,0.0d0)
                  do imu=istart_bfn(katm),istart_bfn(katm+1)-1
                     mu=active_bfn_list(imu)
                     do inu=istart_bfn(latm),istart_bfn(latm+1)-1
                        nu=active_bfn_list(inu)
                        iki = iky(mu)+nu
                        dd=adens(iki)
                        if (abs(dd).gt.dentol) then
                           do kc = 1, 3
                              kx = 3*(kc-1)
                              do lc = 1, 3
                                 klxx = kx+lc
                                 do m=1,3
                                    do ipt=1,npts
                                       ddgrhoa(ipt,m,klxx)
     &                                 = ddgrhoa(ipt,m,klxx)
     &                                 + dd*
     &                                   bfn_hess(ipt,imu,hcc(m,kc))*
     &                                   bfng_val(ipt,inu,lc)
     &                                 + dd*
     &                                   bfng_val(ipt,imu,kc)*
     &                                   bfn_hess(ipt,inu,hcc(m,lc))
                                    enddo ! ipt
                                 enddo ! m
                              enddo ! lc
                           enddo ! kc
                        endif
                        dd=bdens(iki)
                        if (abs(dd).gt.dentol) then
                           do kc = 1, 3
                              kx = 3*(kc-1)
                              do lc = 1, 3
                                 klxx = kx+lc
                                 do m=1,3
                                    do ipt=1,npts
                                       ddgrhob(ipt,m,klxx)
     &                                 = ddgrhob(ipt,m,klxx)
     &                                 + dd*
     &                                   bfn_hess(ipt,imu,hcc(m,kc))*
     &                                   bfng_val(ipt,inu,lc)
     &                                 + dd*
     &                                   bfng_val(ipt,imu,kc)*
     &                                   bfn_hess(ipt,inu,hcc(m,lc))
                                    enddo ! ipt
                                 enddo ! m
                              enddo ! lc
                           enddo ! kc
                        endif
                     enddo ! nu
                  enddo ! mu
               endif ! gradcorr_sw
c
               call aclear_dp(scr,9,0.0d0)
               do kc=1,3
                  kx=3*(kc-1)
                  do lc=1,3
                     klxx = kx+lc
                     do ipt=1,npts
                        scr(klxx)=scr(klxx)
     &                  + wt(ipt)*xc_hpt(ipt,irara)*
     &                    drho(ipt,1,kprt+kc)*drho(ipt,1,lprt+lc)
     &                  + wt(ipt)*xc_hpt(ipt,irbrb)*
     &                    drho(ipt,2,kprt+kc)*drho(ipt,2,lprt+lc)
     &                  + wt(ipt)*xc_hpt(ipt,irarb)*
     &                    (drho(ipt,1,kprt+kc)*drho(ipt,2,lprt+lc)+
     &                     drho(ipt,2,kprt+kc)*drho(ipt,1,lprt+lc))
     &                  + wt(ipt)*xc_vpt(ipt,ira)*ddrhoa(ipt,klxx)
     &                  + wt(ipt)*xc_vpt(ipt,irb)*ddrhob(ipt,klxx)
                     enddo
                     if (gradcorr_sw) then
                        do ipt=1,npts
                           gaada 
     &                     = 2*(dgrho(ipt,1,1,kprt+kc)*grho(ipt,1,1)
     &                        + dgrho(ipt,1,2,kprt+kc)*grho(ipt,1,2)
     &                        + dgrho(ipt,1,3,kprt+kc)*grho(ipt,1,3))
                           gbbda 
     &                     = 2*(dgrho(ipt,2,1,kprt+kc)*grho(ipt,2,1)
     &                        + dgrho(ipt,2,2,kprt+kc)*grho(ipt,2,2)
     &                        + dgrho(ipt,2,3,kprt+kc)*grho(ipt,2,3))
                           gabda 
     &                     =    dgrho(ipt,1,1,kprt+kc)*grho(ipt,2,1)
     &                        + dgrho(ipt,1,2,kprt+kc)*grho(ipt,2,2)
     &                        + dgrho(ipt,1,3,kprt+kc)*grho(ipt,2,3)
     &                        + dgrho(ipt,2,1,kprt+kc)*grho(ipt,1,1)
     &                        + dgrho(ipt,2,2,kprt+kc)*grho(ipt,1,2)
     &                        + dgrho(ipt,2,3,kprt+kc)*grho(ipt,1,3)
                           gaadb
     &                     = 2*(dgrho(ipt,1,1,lprt+lc)*grho(ipt,1,1)
     &                        + dgrho(ipt,1,2,lprt+lc)*grho(ipt,1,2)
     &                        + dgrho(ipt,1,3,lprt+lc)*grho(ipt,1,3))
                           gbbdb 
     &                     = 2*(dgrho(ipt,2,1,lprt+lc)*grho(ipt,2,1)
     &                        + dgrho(ipt,2,2,lprt+lc)*grho(ipt,2,2)
     &                        + dgrho(ipt,2,3,lprt+lc)*grho(ipt,2,3))
                           gabdb
     &                     =    dgrho(ipt,1,1,lprt+lc)*grho(ipt,2,1)
     &                        + dgrho(ipt,1,2,lprt+lc)*grho(ipt,2,2)
     &                        + dgrho(ipt,1,3,lprt+lc)*grho(ipt,2,3)
     &                        + dgrho(ipt,2,1,lprt+lc)*grho(ipt,1,1)
     &                        + dgrho(ipt,2,2,lprt+lc)*grho(ipt,1,2)
     &                        + dgrho(ipt,2,3,lprt+lc)*grho(ipt,1,3)
                           gaadab 
     &                     = 2*(ddgrhoa(ipt,1,klxx)*grho(ipt,1,1)
     &                         +ddgrhoa(ipt,2,klxx)*grho(ipt,1,2)
     &                         +ddgrhoa(ipt,3,klxx)*grho(ipt,1,3)
     &                         +dgrho(ipt,1,1,kprt+kc)*
     &                          dgrho(ipt,1,1,lprt+lc)
     &                         +dgrho(ipt,1,2,kprt+kc)*
     &                          dgrho(ipt,1,2,lprt+lc)
     &                         +dgrho(ipt,1,3,kprt+kc)*
     &                          dgrho(ipt,1,3,lprt+lc))
                           gbbdab 
     &                     = 2*(ddgrhob(ipt,1,klxx)*grho(ipt,2,1)
     &                         +ddgrhob(ipt,2,klxx)*grho(ipt,2,2)
     &                         +ddgrhob(ipt,3,klxx)*grho(ipt,2,3)
     &                         +dgrho(ipt,2,1,kprt+kc)*
     &                          dgrho(ipt,2,1,lprt+lc)
     &                         +dgrho(ipt,2,2,kprt+kc)*
     &                          dgrho(ipt,2,2,lprt+lc)
     &                         +dgrho(ipt,2,3,kprt+kc)*
     &                          dgrho(ipt,2,3,lprt+lc))
                           gabdab 
     &                     =   (ddgrhoa(ipt,1,klxx)*grho(ipt,2,1)
     &                         +ddgrhoa(ipt,2,klxx)*grho(ipt,2,2)
     &                         +ddgrhoa(ipt,3,klxx)*grho(ipt,2,3)
     &                         +ddgrhob(ipt,1,klxx)*grho(ipt,1,1)
     &                         +ddgrhob(ipt,2,klxx)*grho(ipt,1,2)
     &                         +ddgrhob(ipt,3,klxx)*grho(ipt,1,3)
     &                         +dgrho(ipt,1,1,kprt+kc)*
     &                          dgrho(ipt,2,1,lprt+lc)
     &                         +dgrho(ipt,1,2,kprt+kc)*
     &                          dgrho(ipt,2,2,lprt+lc)
     &                         +dgrho(ipt,1,3,kprt+kc)*
     &                          dgrho(ipt,2,3,lprt+lc)
     &                         +dgrho(ipt,2,1,kprt+kc)*
     &                          dgrho(ipt,1,1,lprt+lc)
     &                         +dgrho(ipt,2,2,kprt+kc)*
     &                          dgrho(ipt,1,2,lprt+lc)
     &                         +dgrho(ipt,2,3,kprt+kc)*
     &                          dgrho(ipt,1,3,lprt+lc))
                           scr(kx+lc)=scr(kx+lc)
     &                     + wt(ipt)*xc_dhpt(ipt,iragaa)*
     &                       (drho(ipt,1,kprt+kc)*gaadb
     &                       +drho(ipt,1,lprt+lc)*gaada)
     &                     + wt(ipt)*xc_dhpt(ipt,iragab)*
     &                       (drho(ipt,1,kprt+kc)*gabdb
     &                       +drho(ipt,1,lprt+lc)*gabda)
     &                     + wt(ipt)*xc_dhpt(ipt,iragbb)*
     &                       (drho(ipt,1,kprt+kc)*gbbdb
     &                       +drho(ipt,1,lprt+lc)*gbbda)
     &                     + wt(ipt)*xc_dhpt(ipt,irbgaa)*
     &                       (drho(ipt,2,kprt+kc)*gaadb
     &                       +drho(ipt,2,lprt+lc)*gaada)
     &                     + wt(ipt)*xc_dhpt(ipt,irbgab)*
     &                       (drho(ipt,2,kprt+kc)*gabdb
     &                       +drho(ipt,2,lprt+lc)*gabda)
     &                     + wt(ipt)*xc_dhpt(ipt,irbgbb)*
     &                       (drho(ipt,2,kprt+kc)*gbbdb
     &                       +drho(ipt,2,lprt+lc)*gbbda)
     &                     + wt(ipt)*xc_dhpt(ipt,igaagaa)*
     &                       gaada*gaadb
     &                     + wt(ipt)*xc_dhpt(ipt,igaagab)*
     &                       gaada*gabdb
     &                     + wt(ipt)*xc_dhpt(ipt,igaagbb)*
     &                       gaada*gbbdb
     &                     + wt(ipt)*xc_dhpt(ipt,igaagab)*
     &                       gabda*gaadb
     &                     + wt(ipt)*xc_dhpt(ipt,igabgab)*
     &                       gabda*gabdb
     &                     + wt(ipt)*xc_dhpt(ipt,igabgbb)*
     &                       gabda*gbbdb
     &                     + wt(ipt)*xc_dhpt(ipt,igaagbb)*
     &                       gbbda*gaadb
     &                     + wt(ipt)*xc_dhpt(ipt,igabgbb)*
     &                       gbbda*gabdb
     &                     + wt(ipt)*xc_dhpt(ipt,igbbgbb)*
     &                       gbbda*gbbdb
     &                     + wt(ipt)*xc_dvpt(ipt,igaa)*gaadab
     &                     + wt(ipt)*xc_dvpt(ipt,igab)*gabdab
     &                     + wt(ipt)*xc_dvpt(ipt,igbb)*gbbdab
                        enddo
                     endif
                  enddo
               enddo
               do kc=1,3
                  kx=3*(kc-1)
                  do lc=1,3
                     hess(kpert+kc,lpert+lc) = hess(kpert+kc,lpert+lc)
     &                                       + scr(kx+lc)
                     hess(lpert+lc,kpert+kc) = hess(lpert+lc,kpert+kc)
     &                                       + scr(kx+lc)
                     if (gwt_sw) then
                        hess(kpert+kc,ipert+lc) 
     &                  = hess(kpert+kc,ipert+lc) - scr(kx+lc)
                        hess(ipert+lc,kpert+kc) 
     &                  = hess(ipert+lc,kpert+kc) - scr(kx+lc)
c
                        hess(ipert+kc,lpert+lc) 
     &                  = hess(ipert+kc,lpert+lc) - scr(kx+lc)
                        hess(lpert+lc,ipert+kc) 
     &                  = hess(lpert+lc,ipert+kc) - scr(kx+lc)
c
                        hess(ipert+kc,ipert+lc) 
     &                  = hess(ipert+kc,ipert+lc) + scr(kx+lc)
                        hess(ipert+lc,ipert+kc) 
     &                  = hess(ipert+lc,ipert+kc) + scr(kx+lc)
                     endif
                  enddo
               enddo
 230        continue
c
            call aclear_dp(ddrhoa,npts*6,0.0d0)
            call aclear_dp(ddrhob,npts*6,0.0d0)
c           2nd derivatives of 1 basis function
            do imu=istart_bfn(katm),istart_bfn(katm+1)-1
               mu=active_bfn_list(imu)
               iki=iky(mu)
               do inu=1,imu
                  nu=active_bfn_list(inu)
                  dd=adens(iki+nu)
                  if (abs(dd).gt.dentol) then
                     do kc=1,3
                        kx=kc*(kc-1)/2
                        do lc=1,kc
                           klxx = kx+lc
                           do ipt=1,npts
                              ddrhoa(ipt,klxx) =
     &                        ddrhoa(ipt,klxx)
     &                        + dd*bfn_hess(ipt,imu,hcc(lc,kc))*
     &                          bfn_val(ipt,inu)
                           enddo ! ipt
                        enddo ! lc
                     enddo ! kc
                  endif
                  dd=bdens(iki+nu)
                  if (abs(dd).gt.dentol) then
                     do kc=1,3
                        kx=kc*(kc-1)/2
                        do lc=1,kc
                           klxx = kx+lc
                           do ipt=1,npts
                              ddrhob(ipt,klxx) =
     &                        ddrhob(ipt,klxx)
     &                        + dd*bfn_hess(ipt,imu,hcc(lc,kc))*
     &                          bfn_val(ipt,inu)
                           enddo ! ipt
                        enddo ! lc
                     enddo ! kc
                  endif
               enddo ! nu
               do inu=imu+1, n_active_bfn
                  nu=active_bfn_list(inu)
                  iki=iky(nu)
                  dd=adens(iki+mu)
                  if (abs(dd).gt.dentol) then
                     do kc=1,3
                        kx=kc*(kc-1)/2
                        do lc=1,kc
                           klxx = kx+lc
                           do ipt=1,npts
                              ddrhoa(ipt,klxx) =
     &                        ddrhoa(ipt,klxx)
     &                        + dd*bfn_hess(ipt,imu,hcc(lc,kc))*
     &                          bfn_val(ipt,inu)
                           enddo ! ipt
                        enddo ! lc
                     enddo ! kc
                  endif
                  dd=bdens(iki+mu)
                  if (abs(dd).gt.dentol) then
                     do kc=1,3
                        kx=kc*(kc-1)/2
                        do lc=1,kc
                           klxx = kx+lc
                           do ipt=1,npts
                              ddrhob(ipt,klxx) =
     &                        ddrhob(ipt,klxx)
     &                        + dd*bfn_hess(ipt,imu,hcc(lc,kc))*
     &                          bfn_val(ipt,inu)
                           enddo ! ipt
                        enddo ! lc
                     enddo ! kc
                  endif
               enddo ! nu
            enddo ! mu
c           1st derivative of 2 basis functions
            do imu=istart_bfn(katm),istart_bfn(katm+1)-1
               mu=active_bfn_list(imu)
               do inu=istart_bfn(katm),imu
                  nu=active_bfn_list(inu)
                  iki = iky(mu)+nu
                  dd=adens(iki)
                  if (abs(dd).gt.dentol) then
                     do kc = 1, 3
                        kx = kc*(kc-1)/2
                        do lc = 1, kc
                           klxx = kx+lc
                           do ipt=1,npts
                              ddrhoa(ipt,klxx) = ddrhoa(ipt,klxx)
     &                        + dd*bfng_val(ipt,imu,kc)*
     &                          bfng_val(ipt,inu,lc)
                           enddo ! ipt
                        enddo ! lc
                     enddo ! kc
                  endif
                  dd=bdens(iki)
                  if (abs(dd).gt.dentol) then
                     do kc = 1, 3
                        kx = kc*(kc-1)/2
                        do lc = 1, kc
                           klxx = kx+lc
                           do ipt=1,npts
                              ddrhob(ipt,klxx) = ddrhob(ipt,klxx)
     &                        + dd*bfng_val(ipt,imu,kc)*
     &                          bfng_val(ipt,inu,lc)
                           enddo ! ipt
                        enddo ! lc
                     enddo ! kc
                  endif
               enddo ! nu
               do inu=imu+1,istart_bfn(katm+1)-1
                  nu=active_bfn_list(inu)
                  iki = iky(nu)+mu
                  dd=adens(iki)
                  if (abs(dd).gt.dentol) then
                     do kc = 1, 3
                        kx = kc*(kc-1)/2
                        do lc = 1, kc
                           klxx = kx+lc
                           do ipt=1,npts
                              ddrhoa(ipt,klxx) = ddrhoa(ipt,klxx)
     &                        + dd*bfng_val(ipt,imu,kc)*
     &                          bfng_val(ipt,inu,lc)
                           enddo ! ipt
                        enddo ! lc
                     enddo ! kc
                  endif
                  dd=bdens(iki)
                  if (abs(dd).gt.dentol) then
                     do kc = 1, 3
                        kx = kc*(kc-1)/2
                        do lc = 1, kc
                           klxx = kx+lc
                           do ipt=1,npts
                              ddrhob(ipt,klxx) = ddrhob(ipt,klxx)
     &                        + dd*bfng_val(ipt,imu,kc)*
     &                          bfng_val(ipt,inu,lc)
                           enddo ! ipt
                        enddo ! lc
                     enddo ! kc
                  endif
               enddo ! nu
            enddo ! mu
c
            if (gradcorr_sw) then
               call aclear_dp(ddgrhoa,npts*3*9,0.0d0)
               call aclear_dp(ddgrhob,npts*3*9,0.0d0)
c              2nd derivative of 1 basis function
               do imu=istart_bfn(katm),istart_bfn(katm+1)-1
                  mu=active_bfn_list(imu)
                  iki=iky(mu)
                  do inu=1, imu
                     nu=active_bfn_list(inu)
                     dd=adens(iki+nu)
                     if (abs(dd).gt.dentol) then
                        do kc=1,3
                           kx=kc*(kc-1)/2
                           do lc=1,kc
                              klxx=kx+lc
                              do m=1,3
                                 do ipt=1,npts
                                    ddgrhoa(ipt,m,klxx)
     &                              = ddgrhoa(ipt,m,klxx)
     &                              + dd*(
     &                                 bfn_3rd(ipt,imu,tc3(m,lc,kc))
     &                                *bfn_val(ipt,inu)+
     &                                 bfn_hess(ipt,imu,hcc(lc,kc))
     &                                *bfng_val(ipt,inu,m))
                                 enddo ! ipt
                              enddo ! m
                           enddo ! lc
                        enddo ! kc
                     endif
                     dd=bdens(iki+nu)
                     if (abs(dd).gt.dentol) then
                        do kc=1,3
                           kx=kc*(kc-1)/2
                           do lc=1,kc
                              klxx=kx+lc
                              do m=1,3
                                 do ipt=1,npts
                                    ddgrhob(ipt,m,klxx)
     &                              = ddgrhob(ipt,m,klxx)
     &                              + dd*(
     &                                 bfn_3rd(ipt,imu,tc3(m,lc,kc))
     &                                *bfn_val(ipt,inu)+
     &                                 bfn_hess(ipt,imu,hcc(lc,kc))
     &                                *bfng_val(ipt,inu,m))
                                 enddo ! ipt
                              enddo ! m
                           enddo ! lc
                        enddo ! kc
                     endif
                  enddo ! inu
                  do inu=imu+1, n_active_bfn
                     nu=active_bfn_list(inu)
                     iki=iky(nu)
                     dd=adens(iki+mu)
                     if (abs(dd).gt.dentol) then
                        do kc=1,3
                           kx=kc*(kc-1)/2
                           do lc=1,kc
                              klxx=kx+lc
                              do m=1,3
                                 do ipt=1,npts
                                    ddgrhoa(ipt,m,klxx)
     &                              = ddgrhoa(ipt,m,klxx)
     &                              + dd*(
     &                                 bfn_3rd(ipt,imu,tc3(m,lc,kc))
     &                                *bfn_val(ipt,inu)+
     &                                 bfn_hess(ipt,imu,hcc(lc,kc))
     &                                *bfng_val(ipt,inu,m))
                                 enddo ! ipt
                              enddo ! m
                           enddo ! lc
                        enddo ! kc
                     endif
                     dd=bdens(iki+mu)
                     if (abs(dd).gt.dentol) then
                        do kc=1,3
                           kx=kc*(kc-1)/2
                           do lc=1,kc
                              klxx=kx+lc
                              do m=1,3
                                 do ipt=1,npts
                                    ddgrhob(ipt,m,klxx)
     &                              = ddgrhob(ipt,m,klxx)
     &                              + dd*(
     &                                 bfn_3rd(ipt,imu,tc3(m,lc,kc))
     &                                *bfn_val(ipt,inu)+
     &                                 bfn_hess(ipt,imu,hcc(lc,kc))
     &                                *bfng_val(ipt,inu,m))
                                 enddo ! ipt
                              enddo ! m
                           enddo ! lc
                        enddo ! kc
                     endif
                  enddo ! nu
               enddo ! mu
c              1st derivative of 2 basis functions
               do imu=istart_bfn(katm),istart_bfn(katm+1)-1
                  mu=active_bfn_list(imu)
                  do inu=istart_bfn(katm),imu
                     nu=active_bfn_list(inu)
                     iki = iky(mu)+nu
                     dd=adens(iki)
                     if (abs(dd).gt.dentol) then
                        do kc = 1, 3
                           kx = kc*(kc-1)/2
                           do lc = 1, kc
                              klxx = kx+lc
                              do m=1,3
                                 do ipt=1,npts
                                    ddgrhoa(ipt,m,klxx)
     &                              = ddgrhoa(ipt,m,klxx)
     &                              + dd*
     &                                bfn_hess(ipt,imu,hcc(m,kc))*
     &                                bfng_val(ipt,inu,lc)
     &                              + dd*
     &                                bfng_val(ipt,imu,kc)*
     &                                bfn_hess(ipt,inu,hcc(m,lc))
                                 enddo ! ipt
                              enddo ! m
                           enddo ! lc
                        enddo ! kc
                     endif
                     dd=bdens(iki)
                     if (abs(dd).gt.dentol) then
                        do kc = 1, 3
                           kx = kc*(kc-1)/2
                           do lc = 1, kc
                              klxx = kx+lc
                              do m=1,3
                                 do ipt=1,npts
                                    ddgrhob(ipt,m,klxx)
     &                              = ddgrhob(ipt,m,klxx)
     &                              + dd*
     &                                bfn_hess(ipt,imu,hcc(m,kc))*
     &                                bfng_val(ipt,inu,lc)
     &                              + dd*
     &                                bfng_val(ipt,imu,kc)*
     &                                bfn_hess(ipt,inu,hcc(m,lc))
                                 enddo ! ipt
                              enddo ! m
                           enddo ! lc
                        enddo ! kc
                     endif
                  enddo ! nu
                  do inu=imu+1,istart_bfn(katm+1)-1
                     nu=active_bfn_list(inu)
                     iki = iky(nu)+mu
                     dd=adens(iki)
                     if (abs(dd).gt.dentol) then
                        do kc = 1, 3
                           kx = kc*(kc-1)/2
                           do lc = 1, kc
                              klxx = kx+lc
                              do m=1,3
                                 do ipt=1,npts
                                    ddgrhoa(ipt,m,klxx)
     &                              = ddgrhoa(ipt,m,klxx)
     &                              + dd*
     &                                bfn_hess(ipt,imu,hcc(m,kc))*
     &                                bfng_val(ipt,inu,lc)
     &                              + dd*
     &                                bfng_val(ipt,imu,kc)*
     &                                bfn_hess(ipt,inu,hcc(m,lc))
                                 enddo ! ipt
                              enddo ! m
                           enddo ! lc
                        enddo ! kc
                     endif
                     dd=bdens(iki)
                     if (abs(dd).gt.dentol) then
                        do kc = 1, 3
                           kx = kc*(kc-1)/2
                           do lc = 1, kc
                              klxx = kx+lc
                              do m=1,3
                                 do ipt=1,npts
                                    ddgrhob(ipt,m,klxx)
     &                              = ddgrhob(ipt,m,klxx)
     &                              + dd*
     &                                bfn_hess(ipt,imu,hcc(m,kc))*
     &                                bfng_val(ipt,inu,lc)
     &                              + dd*
     &                                bfng_val(ipt,imu,kc)*
     &                                bfn_hess(ipt,inu,hcc(m,lc))
                                 enddo ! ipt
                              enddo ! m
                           enddo ! lc
                        enddo ! kc
                     endif
                  enddo ! nu
               enddo ! mu
            endif
c
            call aclear_dp(scr,6,0.0d0)
            do kc=1,3
               kx=kc*(kc-1)/2
               do lc=1,kc
                  klxx = kx+lc
                  do ipt=1,npts
                     scr(klxx)=scr(klxx)
     &               + wt(ipt)*xc_hpt(ipt,irara)*
     &                 drho(ipt,1,kprt+kc)*drho(ipt,1,kprt+lc)
     &               + wt(ipt)*xc_hpt(ipt,irbrb)*
     &                 drho(ipt,2,kprt+kc)*drho(ipt,2,kprt+lc)
     &               + wt(ipt)*xc_hpt(ipt,irarb)*
     &                 (drho(ipt,1,kprt+kc)*drho(ipt,2,kprt+lc)+
     &                  drho(ipt,2,kprt+kc)*drho(ipt,1,kprt+lc))
     &               + wt(ipt)*xc_vpt(ipt,ira)*ddrhoa(ipt,klxx)
     &               + wt(ipt)*xc_vpt(ipt,irb)*ddrhob(ipt,klxx)
                  enddo
                  if (gradcorr_sw) then
                     do ipt=1,npts
                        gaada 
     &                  = 2*(dgrho(ipt,1,1,kprt+kc)*grho(ipt,1,1)
     &                     + dgrho(ipt,1,2,kprt+kc)*grho(ipt,1,2)
     &                     + dgrho(ipt,1,3,kprt+kc)*grho(ipt,1,3))
                        gbbda 
     &                  = 2*(dgrho(ipt,2,1,kprt+kc)*grho(ipt,2,1)
     &                     + dgrho(ipt,2,2,kprt+kc)*grho(ipt,2,2)
     &                     + dgrho(ipt,2,3,kprt+kc)*grho(ipt,2,3))
                        gabda 
     &                  =    dgrho(ipt,1,1,kprt+kc)*grho(ipt,2,1)
     &                     + dgrho(ipt,1,2,kprt+kc)*grho(ipt,2,2)
     &                     + dgrho(ipt,1,3,kprt+kc)*grho(ipt,2,3)
     &                     + dgrho(ipt,2,1,kprt+kc)*grho(ipt,1,1)
     &                     + dgrho(ipt,2,2,kprt+kc)*grho(ipt,1,2)
     &                     + dgrho(ipt,2,3,kprt+kc)*grho(ipt,1,3)
                        gaadb 
     &                  = 2*(dgrho(ipt,1,1,kprt+lc)*grho(ipt,1,1)
     &                     + dgrho(ipt,1,2,kprt+lc)*grho(ipt,1,2)
     &                     + dgrho(ipt,1,3,kprt+lc)*grho(ipt,1,3))
                        gbbdb 
     &                  = 2*(dgrho(ipt,2,1,kprt+lc)*grho(ipt,2,1)
     &                     + dgrho(ipt,2,2,kprt+lc)*grho(ipt,2,2)
     &                     + dgrho(ipt,2,3,kprt+lc)*grho(ipt,2,3))
                        gabdb 
     &                  =    dgrho(ipt,1,1,kprt+lc)*grho(ipt,2,1)
     &                     + dgrho(ipt,1,2,kprt+lc)*grho(ipt,2,2)
     &                     + dgrho(ipt,1,3,kprt+lc)*grho(ipt,2,3)
     &                     + dgrho(ipt,2,1,kprt+lc)*grho(ipt,1,1)
     &                     + dgrho(ipt,2,2,kprt+lc)*grho(ipt,1,2)
     &                     + dgrho(ipt,2,3,kprt+lc)*grho(ipt,1,3)
                        gaadab 
     &                  = 2*(ddgrhoa(ipt,1,klxx)*grho(ipt,1,1)
     &                      +ddgrhoa(ipt,2,klxx)*grho(ipt,1,2)
     &                      +ddgrhoa(ipt,3,klxx)*grho(ipt,1,3)
     &                      +dgrho(ipt,1,1,kprt+kc)*
     &                       dgrho(ipt,1,1,kprt+lc)
     &                      +dgrho(ipt,1,2,kprt+kc)*
     &                       dgrho(ipt,1,2,kprt+lc)
     &                      +dgrho(ipt,1,3,kprt+kc)*
     &                       dgrho(ipt,1,3,kprt+lc))
                        gbbdab 
     &                  = 2*(ddgrhob(ipt,1,klxx)*grho(ipt,2,1)
     &                      +ddgrhob(ipt,2,klxx)*grho(ipt,2,2)
     &                      +ddgrhob(ipt,3,klxx)*grho(ipt,2,3)
     &                      +dgrho(ipt,2,1,kprt+kc)*
     &                       dgrho(ipt,2,1,kprt+lc)
     &                      +dgrho(ipt,2,2,kprt+kc)*
     &                       dgrho(ipt,2,2,kprt+lc)
     &                      +dgrho(ipt,2,3,kprt+kc)*
     &                       dgrho(ipt,2,3,kprt+lc))
                        gabdab 
     &                  =   (ddgrhoa(ipt,1,klxx)*grho(ipt,2,1)
     &                      +ddgrhoa(ipt,2,klxx)*grho(ipt,2,2)
     &                      +ddgrhoa(ipt,3,klxx)*grho(ipt,2,3)
     &                      +ddgrhob(ipt,1,klxx)*grho(ipt,1,1)
     &                      +ddgrhob(ipt,2,klxx)*grho(ipt,1,2)
     &                      +ddgrhob(ipt,3,klxx)*grho(ipt,1,3)
     &                      +dgrho(ipt,1,1,kprt+kc)*
     &                       dgrho(ipt,2,1,kprt+lc)
     &                      +dgrho(ipt,1,2,kprt+kc)*
     &                       dgrho(ipt,2,2,kprt+lc)
     &                      +dgrho(ipt,1,3,kprt+kc)*
     &                       dgrho(ipt,2,3,kprt+lc)
     &                      +dgrho(ipt,2,1,kprt+kc)*
     &                       dgrho(ipt,1,1,kprt+lc)
     &                      +dgrho(ipt,2,2,kprt+kc)*
     &                       dgrho(ipt,1,2,kprt+lc)
     &                      +dgrho(ipt,2,3,kprt+kc)*
     &                       dgrho(ipt,1,3,kprt+lc))
                        scr(klxx)=scr(klxx)
     &                  + wt(ipt)*xc_dhpt(ipt,iragaa)*
     &                    (drho(ipt,1,kprt+kc)*gaadb
     &                    +drho(ipt,1,kprt+lc)*gaada)
     &                  + wt(ipt)*xc_dhpt(ipt,iragab)*
     &                    (drho(ipt,1,kprt+kc)*gabdb
     &                    +drho(ipt,1,kprt+lc)*gabda)
     &                  + wt(ipt)*xc_dhpt(ipt,iragbb)*
     &                    (drho(ipt,1,kprt+kc)*gbbdb
     &                    +drho(ipt,1,kprt+lc)*gbbda)
     &                  + wt(ipt)*xc_dhpt(ipt,irbgaa)*
     &                    (drho(ipt,2,kprt+kc)*gaadb
     &                    +drho(ipt,2,kprt+lc)*gaada)
     &                  + wt(ipt)*xc_dhpt(ipt,irbgab)*
     &                    (drho(ipt,2,kprt+kc)*gabdb
     &                    +drho(ipt,2,kprt+lc)*gabda)
     &                  + wt(ipt)*xc_dhpt(ipt,irbgbb)*
     &                    (drho(ipt,2,kprt+kc)*gbbdb
     &                    +drho(ipt,2,kprt+lc)*gbbda)
     &                  + wt(ipt)*xc_dhpt(ipt,igaagaa)*
     &                    gaada*gaadb
     &                  + wt(ipt)*xc_dhpt(ipt,igaagab)*
     &                    gaada*gabdb
     &                  + wt(ipt)*xc_dhpt(ipt,igaagbb)*
     &                    gaada*gbbdb
     &                  + wt(ipt)*xc_dhpt(ipt,igaagab)*
     &                    gabda*gaadb
     &                  + wt(ipt)*xc_dhpt(ipt,igabgab)*
     &                    gabda*gabdb
     &                  + wt(ipt)*xc_dhpt(ipt,igabgbb)*
     &                    gabda*gbbdb
     &                  + wt(ipt)*xc_dhpt(ipt,igaagbb)*
     &                    gbbda*gaadb
     &                  + wt(ipt)*xc_dhpt(ipt,igabgbb)*
     &                    gbbda*gabdb
     &                  + wt(ipt)*xc_dhpt(ipt,igbbgbb)*
     &                    gbbda*gbbdb
     &                  + wt(ipt)*xc_dvpt(ipt,igaa)*gaadab
     &                  + wt(ipt)*xc_dvpt(ipt,igab)*gabdab
     &                  + wt(ipt)*xc_dvpt(ipt,igbb)*gbbdab
                     enddo
                  endif
               enddo
            enddo
            do kc=1,3
               kx=kc*(kc-1)/2
               do lc=1,kc-1
                  hess(kpert+kc,kpert+lc) = hess(kpert+kc,kpert+lc)
     &                                    + scr(kx+lc)
                  hess(kpert+lc,kpert+kc) = hess(kpert+lc,kpert+kc)
     &                                    + scr(kx+lc)
                  if (gwt_sw) then
                     hess(kpert+kc,ipert+lc) 
     &               = hess(kpert+kc,ipert+lc) - scr(kx+lc)
                     hess(ipert+lc,kpert+kc) 
     &               = hess(ipert+lc,kpert+kc) - scr(kx+lc)
c
                     hess(ipert+kc,kpert+lc) 
     &               = hess(ipert+kc,kpert+lc) - scr(kx+lc)
                     hess(kpert+lc,ipert+kc) 
     &               = hess(kpert+lc,ipert+kc) - scr(kx+lc)
c
                     hess(ipert+kc,ipert+lc) 
     &               = hess(ipert+kc,ipert+lc) + scr(kx+lc)
                     hess(ipert+lc,ipert+kc) 
     &               = hess(ipert+lc,ipert+kc) + scr(kx+lc)
                  endif
               enddo
               hess(kpert+kc,kpert+kc) = hess(kpert+kc,kpert+kc)
     &                                 + scr(kx+kc)
               if (gwt_sw) then
                  hess(kpert+kc,ipert+kc) = hess(kpert+kc,ipert+kc)
     &                                    - scr(kx+kc)
                  hess(ipert+kc,kpert+kc) = hess(ipert+kc,kpert+kc)
     &                                    - scr(kx+kc)
                  hess(ipert+kc,ipert+kc) = hess(ipert+kc,ipert+kc)
     &                                    + scr(kx+kc)
               endif
            enddo
 220     continue
      endif
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine rks_hess_dft_mo(gradcorr_sw,gwt_sw,gwt_avail_sw,
     &                           npts,naocc,nao,ngridcent,npert,
     &                           mxp,ao_tag,iatm,
     &                           wt,gwt,g2wt,
     &                           bfng_val,bfn_hess,bfn_3rd,
     &                           xc_ept,xc_vpt,xc_dvpt,xc_hpt,xc_dhpt,
     &                           avec,amo_val,amo_grad,grho,
     &                           drho,dgrho,dgaa,ddrhoa,ddgrhoa,
     &                           bdmo,admo,bddmo,bdgmo,adgmo,bddgmo,t1,
     &                           hess)
      implicit none
c
c     Calculate the DFT contribution to the Hessian in the MO-basis. 
c     The 2nd derivatives of the density and its gradient are calculated
c     on the fly for each atom pair. This way massive memory 
c     requirements can be avoided. 
c
c     Inputs
c
C *Parameter module for CCP1/DFT            
C *Memory information		
C *----------------
C * max_block 	- maximum number of blocks which can be allocate
      integer max_block
      parameter(max_block=20)
C *							
C *Basis sets information                              
C *----------------------				
C *max_tag 	- maximum number of basis sets which can be alloc
C *max_atype	- maximum number of centre types	
C *max_gtype	- maximum number of grid types (at least every element 
C                 has different grid type)
C *max_grids    - maximum number of grids (different terms may have
C                 different grids, i.e. CPKS equations use a different
C                 grid than the one used for the KS-matrix)
C *max_shel	- maximum number of shells on centre
C *max_prm	- maximum number of primitives for any given centre
C *maxL		- maximum angular momentum allowed	
C *max_func	- maximum number of basis functions for any centre
      integer max_tag,max_atype,max_shel,max_prm,max_ang
      integer max_gtype,max_grids
      parameter(max_tag=3,max_atype=30,max_shel=500,max_prm=5000)
      parameter(max_ang=5,max_gtype=10,max_grids=2)
C *								
C *Geometry information					
C *--------------------					
C *max_atom	- maximum number of atoms in system
      integer max_atom
      parameter(max_atom=750)
C *
C *Accuracy information				
C ---------------------			
C *global_accuracy 	- global accuracy	
      real*8  global_accuracy
      parameter(global_accuracy=1.0d-14)

C *Grid information
c *----------------
c *maxradzn - The maximum number of radial zones. 
      integer maxgpt,maxrad,maxfpt,maxang,maxradzn,maxtablerows
      parameter(maxgpt=2900,maxrad=50,maxang=302,maxfpt=100)
      parameter(maxradzn=35,maxtablerows=7)
c
c     The parameters defined here label the 1st, 2nd and 3rd derivative
c     components of the basis functions.
c
c     2nd derivatives:
c
      integer hxx, hyy, hzz, hxy, hxz, hyz
c
      parameter (hxx = 1)
      parameter (hyy = 2)
      parameter (hzz = 3)
      parameter (hxy = 4)
      parameter (hxz = 5)
      parameter (hyz = 6)
c
      integer hcc(3,3)
      data hcc/hxx,hxy,hxz,
     &         hxy,hyy,hyz,
     &         hxz,hyz,hzz/
c
c     3rd derivatives:
c
      integer txxx, txxy, txxz, txyy, txyz
      integer txzz, tyyy, tyyz, tyzz, tzzz
c
      parameter (txxx =  1)
      parameter (txxy =  2)
      parameter (txxz =  3)
      parameter (txyy =  4)
      parameter (txyz =  5)
      parameter (txzz =  6)
      parameter (tyyy =  7)
      parameter (tyyz =  8)
      parameter (tyzz =  9)
      parameter (tzzz = 10)
c
      integer tc3(3,3,3)
      data tc3/txxx,txxy,txxz,
     &         txxy,txyy,txyz,
     &         txxz,txyz,txzz,
     &         txxy,txyy,txyz,
     &         txyy,tyyy,tyyz,
     &         txyz,tyyz,tyzz,
     &         txxz,txyz,txzz,
     &         txyz,tyyz,tyzz,
     &         txzz,tyzz,tzzz/
c
c     The parameters defined here label the 1st and 2nd derivative
c     components of the density functional. Because a functional
c     may be purely local or of a GGA type the components are grouped
c     into derivative with respect to the density only and others.
c
c     1st derivatives:
c
      integer ira, irb
      integer igaa, igab, igbb
      integer ita, itb
c
c     2nd derivatives:
c
      integer irara, irarb, irbrb
      integer iragaa, iragab, iragbb, irbgaa, irbgab, irbgbb
      integer igaagaa, igaagab, igaagbb, igabgab, igabgbb, igbbgbb
      integer irata, iratb, irbta, irbtb
      integer igaata, igaatb, igabta, igabtb, igbbta, igbbtb
      integer itata, itatb, itbtb
c
      parameter (ira     = 1)
      parameter (irb     = 2)
c
      parameter (igaa    = 1)
      parameter (igab    = 2)
      parameter (igbb    = 3)
c
      parameter (ita     = 1)
      parameter (itb     = 2)
c
c
      parameter (irara   = 1)
      parameter (irarb   = 2)
      parameter (irbrb   = 3)
c
      parameter (iragaa  = 1)
      parameter (iragab  = 3)
      parameter (iragbb  = 4)
      parameter (irbgaa  = 5)
      parameter (irbgab  = 6)
      parameter (irbgbb  = 7)
      parameter (igaagaa = 2)
      parameter (igaagab = 8)
      parameter (igaagbb = 9)
      parameter (igabgab = 10)
      parameter (igabgbb = 11)
      parameter (igbbgbb = 12)
c
      parameter (irata   = 1)
      parameter (iratb   = 3)
      parameter (irbta   = 4)
      parameter (irbtb   = 5)
      parameter (itata   = 2)
      parameter (itatb   = 6)
      parameter (itbtb   = 7)
c
      parameter (igaata  = 1)
      parameter (igaatb  = 2)
      parameter (igabta  = 3)
      parameter (igabtb  = 4)
      parameter (igbbta  = 5)
      parameter (igbbtb  = 6)
c
C *Basis sets include file
c
      real*8  alpha
      real*8  cont_coeff
      integer num_bset,bset_tags
      common/basis_sets/alpha(max_tag,max_atype,max_prm),
     &                  cont_coeff(max_tag,max_atype,max_prm,max_ang),
     &                  num_bset,bset_tags(max_tag)

      integer Ashl, Aprm, Abfn, totshl, totprm, totbfn, size_shlA
      integer size_basA, size_primA, maxi_shlA, maxi_basA, maxi_primA
      integer num_types, atom_tag , num_shl, atm_typ, nprim, angmom
      integer hybrid, pstart

      common/basis_size_info/Ashl(max_tag,max_atype),
     &                       Aprm(max_tag,max_atype),
     &                       Abfn(max_tag,max_atype),
     &                       totshl(max_tag),
     &                       totprm(max_tag),
     &                       totbfn(max_tag),
     &                       size_shlA,maxi_shlA,
     &                       size_basA,maxi_basA,
     &                       size_primA,maxi_primA
C
C Descriptions
C 
C  num_types		-	number of types of atoms for a given basis
C  atom_tag		-	list of atom type id numbers
C
      common/basis_cent_info/num_types(max_tag),
     &                       atom_tag(max_tag,max_atom)
C
C Descriptions
C
C  num_shl		-	number of shells for given atom type
C  atm_typ		-	atomic number of atom type
C  nprim		-  	number of primitives 
C  angmom		-	angular momentum
C  hybrid		-	level of hybrid shell. Is same as angmom if
C				shell is not a hybrid. For hybrid shells is
C				always less than angmom e.g. for an sp shell
C				angmom=2, hybrid=1
C  pstart		-	start of exponents and contraction coeffs 
C				contained in include basis.hf77, for a given
C				shell
C
      common/basis_cont_info/num_shl(max_tag,max_atype),
     &                       atm_typ(max_tag,max_atype),
     &                       nprim(max_tag,max_atype,max_shel),
     &                       angmom(max_tag,max_atype,max_shel),
     &                       hybrid(max_tag,max_atype,max_shel),
     &                       pstart(max_tag,max_atype,max_shel)

c
      logical gradcorr_sw
      logical gwt_sw
      logical gwt_avail_sw
c
      integer npts
      integer naocc
      integer nao
      integer ngridcent
      integer npert
      integer mxp
      integer ao_tag
      integer iatm
c
      real*8 wt(mxp)
      real*8 gwt(3,mxp,ngridcent)
      real*8 g2wt(mxp,npert,npert)
c
      real*8 bfn_val(mxp,nao)
      real*8 bfng_val(mxp,nao,3)
      real*8 bfn_hess(mxp,nao,6)
      real*8 bfn_3rd(mxp,nao,10)
c
      real*8 xc_ept(mxp)
      real*8 xc_vpt(mxp,2)
      real*8 xc_dvpt(mxp,3)
      real*8 xc_hpt(mxp,3)
      real*8 xc_dhpt(mxp,12)
c
      real*8 avec(nao,naocc)
      real*8 amo_val(npts,naocc)
      real*8 amo_grad(npts,naocc,3)
      real*8 grho(mxp,2,3)
c
c     Outputs
c
      real*8 hess(npert,npert)
c
c     Workspave
c
      real*8 drho(npts,npert)
      real*8 dgrho(npts,3,npert)
      real*8 dgaa(npts,npert)
      real*8 ddrhoa(npts,9)
      real*8 ddgrhoa(npts,3,9)
      real*8 bdmo(npts,3)
      real*8 admo(npts,3)
      real*8 bddmo(npts,6)
      real*8 bdgmo(npts,6)
      real*8 adgmo(npts,6)
      real*8 bddgmo(npts,10)
      real*8 t1(npts)
c
c     Local
c
      integer ipt
      integer katm, latm
      integer kpert, lpert, ipert
      integer kc, lc, mc
      integer klxx, kx
      integer mu, nu
      integer ith
      integer batmt, bstart, bend
      integer aatmt, astart, aend
c
      real*8 scr(9)
      real*8 gaadab
c
c     Code
c
cDEBUG
c     call aclear_dp(xc_vpt(1,ira),npts,0.0d0)
c     call aclear_dp(xc_hpt(1,irara),npts,0.0d0)
c     call aclear_dp(xc_dvpt(1,igaa),npts,0.0d0)
c     call aclear_dp(xc_dhpt(1,iragaa),npts,0.0d0)
c     call aclear_dp(xc_dhpt(1,igaagaa),npts,0.0d0)
cDEBUG
      ipert=3*(iatm-1)
c
c      if (gwt_sw.and.gwt_avail_sw) then
c         do kpert=1,npert
c            do lpert=1,npert
c               do ipt=1,npts
c                  hess(lpert,kpert)=hess(lpert,kpert)
c     &            +g2wt(ipt,lpert,kpert)*xc_ept(ipt)
c               enddo
c            enddo ! lpert
c         enddo ! kpert
c      endif
c
c     evaluate drho and dgrho
c
      call aclear_dp(drho,npts*npert,0.0d0)
      if (gradcorr_sw) then
         call aclear_dp(dgrho,3*npts*npert,0.0d0)
      endif
      bend = 0
      do katm=1,ngridcent
         kpert=3*(katm-1)
         batmt=atom_tag(ao_tag,katm)
         bstart=bend+1
         bend=bend+Abfn(ao_tag,batmt)
         do kc=1,3
            do ith=1,naocc
               call aclear_dp(t1,npts,0.0d0)
               do mu=bstart,bend
                  do ipt=1,npts
                     t1(ipt)=t1(ipt)
     &               - avec(mu,ith)*bfng_val(ipt,mu,kc)
                  enddo ! ipt
               enddo ! mu
c
c              t1 = dmo(ipt,kpert+k,ith)
c
               do ipt=1,npts 
                  drho(ipt,kpert+kc)=drho(ipt,kpert+kc)
     &            + amo_val(ipt,ith)*t1(ipt)
               enddo ! ipt
c
               if (gradcorr_sw) then
                  do lc=1,3
                     do ipt=1,npts
                        dgrho(ipt,lc,kpert+kc)=dgrho(ipt,lc,kpert+kc)
     &                  + amo_grad(ipt,ith,lc)*t1(ipt)
                     enddo ! ipt
                  enddo ! lc
               endif ! gradcorr_sw
            enddo ! ith
         enddo ! kc
      enddo ! katm
c
      if (gradcorr_sw) then
         bend = 0
         do katm=1,ngridcent
            kpert=3*(katm-1)
            batmt=atom_tag(ao_tag,katm)
            bstart=bend+1
            bend=bend+Abfn(ao_tag,batmt)
            do kc=1,3
               do lc=1,3
                  do ith=1,naocc
                     call aclear_dp(t1,npts,0.0d0)
                     do mu=bstart,bend
                        do ipt=1,npts
                           t1(ipt)=t1(ipt)
     &                     - avec(mu,ith)*bfn_hess(ipt,mu,hcc(lc,kc))
                        enddo ! ipt
                     enddo ! mu
c
c                    t1 = dgmo(ipt,lc,kpert+kc,ith)
c
                     do ipt=1,npts
                        dgrho(ipt,lc,kpert+kc)=dgrho(ipt,lc,kpert+kc)
     &                  + amo_val(ipt,ith)*t1(ipt)
                     enddo ! ipt
                  enddo ! ith
               enddo ! lc
            enddo ! kc
         enddo ! katm
c
         call aclear_dp(dgaa,npts*npert,0.0d0)
         do katm=1,ngridcent
            kpert=3*(katm-1)
            do kc=1,3
               do lc=1,3
                  do ipt=1,npts
                     dgaa(ipt,kpert+kc)=dgaa(ipt,kpert+kc)
     &               + dgrho(ipt,lc,kpert+kc)*grho(ipt,1,lc)
                  enddo ! ipt
               enddo ! lc
            enddo ! kc
         enddo ! katm
      endif ! gradcorr_sw
c
c     Now do the terms with the gradients of the weights
c
      if (gwt_sw.and.gwt_avail_sw) then
         do 100 katm=1,ngridcent
            if (katm.eq.iatm) goto 100
            kpert=3*(katm-1)
            do 110 latm=1,ngridcent
c              if (latm.eq.iatm) goto 110
               lpert=3*(latm-1)
               call aclear_dp(scr,9,0.0d0)
               do kc=1,3
                  kx=3*(kc-1)
                  do lc=1,3
                     klxx=kx+lc
                     do ipt=1,npts
                        scr(klxx)=scr(klxx)
     &                  + 4*gwt(lc,ipt,latm)*xc_vpt(ipt,ira)*
     &                    drho(ipt,kpert+kc)
                     enddo ! ipt
c                    scr(klxx)=4*scr(klxx)
                  enddo ! lc
               enddo ! kc
               if (gradcorr_sw) then
                  do kc=1,3
                     kx=3*(kc-1)
                     do lc=1,3
                        klxx=kx+lc
                        do ipt=1,npts
                           scr(klxx)=scr(klxx)
     &                     + 2*gwt(lc,ipt,latm)*xc_dvpt(ipt,igaa)*
     &                       dgaa(ipt,kpert+kc)
                        enddo ! ipt
                     enddo ! lc
                  enddo ! kc
               endif ! gradcorr_sw
               do kc=1,3
                  kx=3*(kc-1)
                  do lc=1,3
                     klxx=kx+lc
                     hess(lpert+lc,kpert+kc)
     &               =hess(lpert+lc,kpert+kc)+scr(klxx)
                     hess(lpert+lc,ipert+kc)
     &               =hess(lpert+lc,ipert+kc)-scr(klxx)
                     hess(kpert+kc,lpert+lc)
     &               =hess(kpert+kc,lpert+lc)+scr(klxx)
                     hess(ipert+kc,lpert+lc)
     &               =hess(ipert+kc,lpert+lc)-scr(klxx)
                  enddo ! lc
               enddo ! kc
 110        continue ! latm
 100     continue ! katm
      endif
c
      bend=0
      do 120 katm=1,ngridcent
         batmt  = atom_tag(ao_tag,katm)
         bstart = bend+1
         bend   = bend+Abfn(ao_tag,batmt)
         if (gwt_sw.and.katm.eq.iatm) goto 120
         kpert  = 3*(katm-1)
         aend=0
         do 130 latm=1,katm-1
            aatmt  = atom_tag(ao_tag,latm)
            astart = aend+1
            aend   = aend+Abfn(ao_tag,aatmt)
            if (gwt_sw.and.latm.eq.iatm) goto 130
            lpert  = 3*(latm-1)
c
            call aclear_dp(ddrhoa,npts*9,0.0d0)
            if (gradcorr_sw) then
               call aclear_dp(ddgrhoa,npts*3*9,0.0d0)
            endif
            do ith=1,naocc
               call aclear_dp(admo,npts*3,0.0d0)
               call aclear_dp(bdmo,npts*3,0.0d0)
               do kc=1,3
                  do mu=bstart,bend
                     do ipt=1,npts
                        bdmo(ipt,kc)=bdmo(ipt,kc)
     &                  -avec(mu,ith)*bfng_val(ipt,mu,kc)
                     enddo ! ipt
                  enddo ! mu
               enddo ! kc
               do kc=1,3
                  do mu=astart,aend
                     do ipt=1,npts
                        admo(ipt,kc)=admo(ipt,kc)
     &                  -avec(mu,ith)*bfng_val(ipt,mu,kc)
                     enddo ! ipt
                  enddo ! mu
               enddo ! kc
               do kc=1,3
                  kx=3*(kc-1)
                  do lc=1,3
                     klxx=kx+lc
                     do ipt=1,npts
                        ddrhoa(ipt,klxx)=ddrhoa(ipt,klxx)
     &                  +bdmo(ipt,kc)*admo(ipt,lc)
                     enddo
                  enddo ! lc
               enddo ! kc
               if (gradcorr_sw) then
                  call aclear_dp(adgmo,npts*6,0.0d0)
                  call aclear_dp(bdgmo,npts*6,0.0d0)
                  do kx=1,6
                     do mu=bstart,bend
                        do ipt=1,npts
                           bdgmo(ipt,kx)=bdgmo(ipt,kx)
     &                     -avec(mu,ith)*bfn_hess(ipt,mu,kx)
                        enddo ! ipt
                     enddo ! mu
                  enddo ! kx
                  do kx=1,6
                     do mu=astart,aend
                        do ipt=1,npts
                           adgmo(ipt,kx)=adgmo(ipt,kx)
     &                     -avec(mu,ith)*bfn_hess(ipt,mu,kx)
                        enddo ! ipt
                     enddo ! mu
                  enddo ! kx
                  do kc=1,3
                     kx=3*(kc-1)
                     do lc=1,3
                        klxx=kx+lc
                        do mc=1,3
                           do ipt=1,npts
                              ddgrhoa(ipt,mc,klxx)=ddgrhoa(ipt,mc,klxx)
     &                        +bdgmo(ipt,hcc(mc,kc))*admo(ipt,lc)
     &                        +bdmo(ipt,kc)*adgmo(ipt,hcc(mc,lc))
                           enddo
                        enddo ! mc
                     enddo ! lc
                  enddo ! kc
               endif
            enddo ! ith
c
            call aclear_dp(scr,9,0.0d0)
            do kc=1,3
               kx=3*(kc-1)
               do lc=1,3
                  klxx = kx+lc
                  do ipt=1,npts
                     scr(klxx)=scr(klxx)
     &               + 8*wt(ipt)*xc_hpt(ipt,irara)*
     &                 drho(ipt,kpert+kc)*drho(ipt,lpert+lc)
     &               + 4*wt(ipt)*xc_vpt(ipt,ira)*ddrhoa(ipt,klxx)
                  enddo
                  if (gradcorr_sw) then
                     do ipt=1,npts
                        gaadab  
     &                  =   (ddgrhoa(ipt,1,klxx)*grho(ipt,1,1)
     &                      +ddgrhoa(ipt,2,klxx)*grho(ipt,1,2)
     &                      +ddgrhoa(ipt,3,klxx)*grho(ipt,1,3))
     &                  + 4*(dgrho(ipt,1,kpert+kc)*
     &                       dgrho(ipt,1,lpert+lc)
     &                      +dgrho(ipt,2,kpert+kc)*
     &                       dgrho(ipt,2,lpert+lc)
     &                      +dgrho(ipt,3,kpert+kc)*
     &                       dgrho(ipt,3,lpert+lc))
                        scr(kx+lc)=scr(kx+lc)
     &                  + 8*wt(ipt)*xc_dhpt(ipt,iragaa)*
     &                    (drho(ipt,kpert+kc)*dgaa(ipt,lpert+lc) 
     &                    +drho(ipt,lpert+lc)*dgaa(ipt,kpert+kc))
     &                  + 4*wt(ipt)*xc_dhpt(ipt,igaagaa)*
     &                    dgaa(ipt,lpert+lc)*dgaa(ipt,kpert+kc)
     &                  + 2*wt(ipt)*xc_dvpt(ipt,igaa)*gaadab
                     enddo
                  endif
               enddo
            enddo
            do kc=1,3
               kx=3*(kc-1)
               do lc=1,3
                  hess(kpert+kc,lpert+lc) = hess(kpert+kc,lpert+lc)
     &                                    + scr(kx+lc)
                  hess(lpert+lc,kpert+kc) = hess(lpert+lc,kpert+kc)
     &                                    + scr(kx+lc)
                  if (gwt_sw) then
                     hess(kpert+kc,ipert+lc) 
     &               = hess(kpert+kc,ipert+lc) - scr(kx+lc)
                     hess(ipert+lc,kpert+kc) 
     &               = hess(ipert+lc,kpert+kc) - scr(kx+lc)
c
                     hess(ipert+kc,lpert+lc) 
     &               = hess(ipert+kc,lpert+lc) - scr(kx+lc)
                     hess(lpert+lc,ipert+kc) 
     &               = hess(lpert+lc,ipert+kc) - scr(kx+lc)
c
                     hess(ipert+kc,ipert+lc) 
     &               = hess(ipert+kc,ipert+lc) + scr(kx+lc)
                     hess(ipert+lc,ipert+kc) 
     &               = hess(ipert+lc,ipert+kc) + scr(kx+lc)
                  endif
               enddo
            enddo
 130     continue
c
         call aclear_dp(ddrhoa,npts*6,0.0d0)
c        2nd derivative of 1 basis function
         do ith=1,naocc
            call aclear_dp(bddmo,npts*6,0.0d0)
            do kc=1,3
               kx=kc*(kc-1)/2
               do lc=1,kc
                  klxx=kx+lc
                  do mu=bstart,bend
                     do ipt=1,npts
                        bddmo(ipt,klxx)=bddmo(ipt,klxx)
     &                  +avec(mu,ith)*bfn_hess(ipt,mu,klxx)
                     enddo ! ipt
                  enddo ! mu
               enddo ! lc
            enddo ! kc
            do kc=1,3
               kx=kc*(kc-1)/2
               do lc=1,kc
                  klxx=kx+lc
                  do ipt=1,npts
                     ddrhoa(ipt,klxx)=ddrhoa(ipt,klxx)
     &               +amo_val(ipt,ith)*bddmo(ipt,hcc(lc,kc))
                  enddo ! ipt
               enddo ! lc
            enddo ! kc
         enddo ! ith
c        1st derivative of 2 basis functions
         do ith=1,naocc
            call aclear_dp(bdmo,npts*3,0.0d0)
            do kc=1,3
               do mu=bstart,bend
                  do ipt=1,npts
                     bdmo(ipt,kc)=bdmo(ipt,kc)
     &               -avec(mu,ith)*bfng_val(ipt,mu,kc)
                  enddo ! ipt
               enddo ! mu
            enddo ! kc
            do kc=1,3
               kx=kc*(kc-1)/2
               do lc=1,kc
                  klxx=kx+lc
                  do ipt=1,npts
                     ddrhoa(ipt,klxx)=ddrhoa(ipt,klxx)
     &               +bdmo(ipt,kc)*bdmo(ipt,lc)
                  enddo ! ipt
               enddo ! lc
            enddo ! kc
         enddo ! ith
c
         if (gradcorr_sw) then
            call aclear_dp(ddgrhoa,npts*3*6,0.0d0)
            do ith=1,naocc
c              2nd derivative of 1 basis function
               call aclear_dp(bddmo,npts*6,0.0d0)
               call aclear_dp(bddgmo,npts*10,0.0d0)
               do kx=1,6
                  do mu=bstart,bend
                     do ipt=1,npts
                        bddmo(ipt,kx)=bddmo(ipt,kx)
     &                  +avec(mu,ith)*bfn_hess(ipt,mu,kx)
                     enddo ! ipt
                  enddo ! mu
               enddo ! kx
               do kx=1,10
                  do mu=bstart,bend
                     do ipt=1,npts
                        bddgmo(ipt,kx)=bddgmo(ipt,kx)
     &                  +avec(mu,ith)*bfn_3rd(ipt,mu,kx)
                     enddo ! ipt
                  enddo ! mu
               enddo ! kx
               do kc=1,3
                  kx=kc*(kc-1)/2
                  do lc=1,kc
                     klxx=kx+lc
                     do mc=1,3
                        do ipt=1,npts
                           ddgrhoa(ipt,mc,klxx)=ddgrhoa(ipt,mc,klxx)
     &                     + amo_val(ipt,ith)*bddgmo(ipt,tc3(mc,lc,kc))
     &                     + amo_grad(ipt,ith,mc)*bddmo(ipt,hcc(lc,kc))
                        enddo 
                     enddo ! mc
                  enddo ! lc
               enddo ! kc
c              1st derivative of 2 basis functions
c              ... nasty ... use the fact that -*- = +
               call aclear_dp(bdmo,npts*3,0.0d0)
               do kc=1,3
                  do mu=bstart,bend
                     do ipt=1,npts
                        bdmo(ipt,kc)=bdmo(ipt,kc)
     &                  +avec(mu,ith)*bfng_val(ipt,mu,kc)
                     enddo ! ipt
                  enddo ! mu
               enddo ! kc
               do kc=1,3
                  kx=kc*(kc-1)/2
                  do lc=1,kc
                     klxx=kx+lc
                     do mc=1,3
                        do ipt=1,npts
                           ddgrhoa(ipt,mc,klxx)=ddgrhoa(ipt,mc,klxx)
     &                     + bdmo(ipt,lc)*bddmo(ipt,hcc(mc,kc))
     &                     + bddmo(ipt,hcc(mc,lc))*bdmo(ipt,kc)
                        enddo ! ipt
                     enddo ! mc
                  enddo ! lc
               enddo ! kc
            enddo ! ith
         endif
c
         call aclear_dp(scr,6,0.0d0)
         do kc=1,3
            kx=kc*(kc-1)/2
            do lc=1,kc
               klxx = kx+lc
               do ipt=1,npts
                  scr(klxx)=scr(klxx)
     &            + 8*wt(ipt)*xc_hpt(ipt,irara)*
     &              drho(ipt,kpert+kc)*drho(ipt,kpert+lc)
     &            + 4*wt(ipt)*xc_vpt(ipt,ira)*ddrhoa(ipt,klxx)
               enddo
               if (gradcorr_sw) then
                  do ipt=1,npts
c                    gaada = dgaa(ipt,kpert+kc)
c    &               =   (dgrho(ipt,1,1,kpert+kc)*grho(ipt,1,1)
c    &                  + dgrho(ipt,1,2,kpert+kc)*grho(ipt,1,2)
c    &                  + dgrho(ipt,1,3,kpert+kc)*grho(ipt,1,3))
c                    gaadb = dgaa(ipt,kpert+lc)
c    &               =   (dgrho(ipt,1,1,kpert+lc)*grho(ipt,1,1)
c    &                  + dgrho(ipt,1,2,kpert+lc)*grho(ipt,1,2)
c    &                  + dgrho(ipt,1,3,kpert+lc)*grho(ipt,1,3))
                     gaadab 
     &               =   (ddgrhoa(ipt,1,klxx)*grho(ipt,1,1)
     &                   +ddgrhoa(ipt,2,klxx)*grho(ipt,1,2)
     &                   +ddgrhoa(ipt,3,klxx)*grho(ipt,1,3))
     &               + 4*(dgrho(ipt,1,kpert+kc)*
     &                    dgrho(ipt,1,kpert+lc)
     &                   +dgrho(ipt,2,kpert+kc)*
     &                    dgrho(ipt,2,kpert+lc)
     &                   +dgrho(ipt,3,kpert+kc)*
     &                    dgrho(ipt,3,kpert+lc))
                     scr(klxx)=scr(klxx)
     &               + 8*wt(ipt)*xc_dhpt(ipt,iragaa)*
     &                 (drho(ipt,kpert+kc)*dgaa(ipt,kpert+lc)
     &                 +drho(ipt,kpert+lc)*dgaa(ipt,kpert+kc))
     &               + 4*wt(ipt)*xc_dhpt(ipt,igaagaa)*
     &                 dgaa(ipt,kpert+kc)*dgaa(ipt,kpert+lc)
     &               + 2*wt(ipt)*xc_dvpt(ipt,igaa)*gaadab
                  enddo
               endif
            enddo
         enddo
         do kc=1,3
            kx=kc*(kc-1)/2
            do lc=1,kc-1
               hess(kpert+kc,kpert+lc) = hess(kpert+kc,kpert+lc)
     &                                 + scr(kx+lc)
               hess(kpert+lc,kpert+kc) = hess(kpert+lc,kpert+kc)
     &                                 + scr(kx+lc)
               if (gwt_sw) then
                  hess(kpert+kc,ipert+lc) 
     &            = hess(kpert+kc,ipert+lc) - scr(kx+lc)
                  hess(ipert+lc,kpert+kc) 
     &            = hess(ipert+lc,kpert+kc) - scr(kx+lc)
c
                  hess(ipert+kc,kpert+lc) 
     &            = hess(ipert+kc,kpert+lc) - scr(kx+lc)
                  hess(kpert+lc,ipert+kc) 
     &            = hess(kpert+lc,ipert+kc) - scr(kx+lc)
c
                  hess(ipert+kc,ipert+lc) 
     &            = hess(ipert+kc,ipert+lc) + scr(kx+lc)
                  hess(ipert+lc,ipert+kc) 
     &            = hess(ipert+lc,ipert+kc) + scr(kx+lc)
               endif
            enddo
            hess(kpert+kc,kpert+kc) = hess(kpert+kc,kpert+kc)
     &                              + scr(kx+kc)
            if (gwt_sw) then
               hess(kpert+kc,ipert+kc) = hess(kpert+kc,ipert+kc)
     &                                 - scr(kx+kc)
               hess(ipert+kc,kpert+kc) = hess(ipert+kc,kpert+kc)
     &                                 - scr(kx+kc)
               hess(ipert+kc,ipert+kc) = hess(ipert+kc,ipert+kc)
     &                                 + scr(kx+kc)
            endif
         enddo
 120  continue
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine rks_hess_dft_mo_scr(
     &           gradcorr_sw,gwt_sw,gwt_avail_sw,
     &           active_bfn_list,active_bfn_indx,active_bfn_atms,
     &           n_active_bfn,n_active_atm,
     &           near_atom_list,num_near_atoms,
     &           npts,naocc,nao,ngridcent,npert,
     &           mxp,ao_tag,iatm,
     &           wt,gwt,g2wt,
     &           bfng_val,bfn_hess,bfn_3rd,
     &           xc_ept,xc_vpt,xc_dvpt,xc_hpt,xc_dhpt,
     &           avec,amo_val,amo_grad,grho,
     &           drho,dgrho,dgaa,ddrhoa,ddgrhoa,
     &           bdmo,admo,bddmo,bdgmo,adgmo,bddgmo,t1,istart_bfn,
     &           hess)
      implicit none
c
c     Calculate the DFT contribution to the Hessian in the MO-basis. 
c     The 2nd derivatives of the density and its gradient are calculated
c     on the fly for each atom pair. This way massive memory 
c     requirements can be avoided. 
c
c     Inputs
c
C *Parameter module for CCP1/DFT            
C *Memory information		
C *----------------
C * max_block 	- maximum number of blocks which can be allocate
      integer max_block
      parameter(max_block=20)
C *							
C *Basis sets information                              
C *----------------------				
C *max_tag 	- maximum number of basis sets which can be alloc
C *max_atype	- maximum number of centre types	
C *max_gtype	- maximum number of grid types (at least every element 
C                 has different grid type)
C *max_grids    - maximum number of grids (different terms may have
C                 different grids, i.e. CPKS equations use a different
C                 grid than the one used for the KS-matrix)
C *max_shel	- maximum number of shells on centre
C *max_prm	- maximum number of primitives for any given centre
C *maxL		- maximum angular momentum allowed	
C *max_func	- maximum number of basis functions for any centre
      integer max_tag,max_atype,max_shel,max_prm,max_ang
      integer max_gtype,max_grids
      parameter(max_tag=3,max_atype=30,max_shel=500,max_prm=5000)
      parameter(max_ang=5,max_gtype=10,max_grids=2)
C *								
C *Geometry information					
C *--------------------					
C *max_atom	- maximum number of atoms in system
      integer max_atom
      parameter(max_atom=750)
C *
C *Accuracy information				
C ---------------------			
C *global_accuracy 	- global accuracy	
      real*8  global_accuracy
      parameter(global_accuracy=1.0d-14)

C *Grid information
c *----------------
c *maxradzn - The maximum number of radial zones. 
      integer maxgpt,maxrad,maxfpt,maxang,maxradzn,maxtablerows
      parameter(maxgpt=2900,maxrad=50,maxang=302,maxfpt=100)
      parameter(maxradzn=35,maxtablerows=7)
c
c     The parameters defined here label the 1st, 2nd and 3rd derivative
c     components of the basis functions.
c
c     2nd derivatives:
c
      integer hxx, hyy, hzz, hxy, hxz, hyz
c
      parameter (hxx = 1)
      parameter (hyy = 2)
      parameter (hzz = 3)
      parameter (hxy = 4)
      parameter (hxz = 5)
      parameter (hyz = 6)
c
      integer hcc(3,3)
      data hcc/hxx,hxy,hxz,
     &         hxy,hyy,hyz,
     &         hxz,hyz,hzz/
c
c     3rd derivatives:
c
      integer txxx, txxy, txxz, txyy, txyz
      integer txzz, tyyy, tyyz, tyzz, tzzz
c
      parameter (txxx =  1)
      parameter (txxy =  2)
      parameter (txxz =  3)
      parameter (txyy =  4)
      parameter (txyz =  5)
      parameter (txzz =  6)
      parameter (tyyy =  7)
      parameter (tyyz =  8)
      parameter (tyzz =  9)
      parameter (tzzz = 10)
c
      integer tc3(3,3,3)
      data tc3/txxx,txxy,txxz,
     &         txxy,txyy,txyz,
     &         txxz,txyz,txzz,
     &         txxy,txyy,txyz,
     &         txyy,tyyy,tyyz,
     &         txyz,tyyz,tyzz,
     &         txxz,txyz,txzz,
     &         txyz,tyyz,tyzz,
     &         txzz,tyzz,tzzz/
c
c     The parameters defined here label the 1st and 2nd derivative
c     components of the density functional. Because a functional
c     may be purely local or of a GGA type the components are grouped
c     into derivative with respect to the density only and others.
c
c     1st derivatives:
c
      integer ira, irb
      integer igaa, igab, igbb
      integer ita, itb
c
c     2nd derivatives:
c
      integer irara, irarb, irbrb
      integer iragaa, iragab, iragbb, irbgaa, irbgab, irbgbb
      integer igaagaa, igaagab, igaagbb, igabgab, igabgbb, igbbgbb
      integer irata, iratb, irbta, irbtb
      integer igaata, igaatb, igabta, igabtb, igbbta, igbbtb
      integer itata, itatb, itbtb
c
      parameter (ira     = 1)
      parameter (irb     = 2)
c
      parameter (igaa    = 1)
      parameter (igab    = 2)
      parameter (igbb    = 3)
c
      parameter (ita     = 1)
      parameter (itb     = 2)
c
c
      parameter (irara   = 1)
      parameter (irarb   = 2)
      parameter (irbrb   = 3)
c
      parameter (iragaa  = 1)
      parameter (iragab  = 3)
      parameter (iragbb  = 4)
      parameter (irbgaa  = 5)
      parameter (irbgab  = 6)
      parameter (irbgbb  = 7)
      parameter (igaagaa = 2)
      parameter (igaagab = 8)
      parameter (igaagbb = 9)
      parameter (igabgab = 10)
      parameter (igabgbb = 11)
      parameter (igbbgbb = 12)
c
      parameter (irata   = 1)
      parameter (iratb   = 3)
      parameter (irbta   = 4)
      parameter (irbtb   = 5)
      parameter (itata   = 2)
      parameter (itatb   = 6)
      parameter (itbtb   = 7)
c
      parameter (igaata  = 1)
      parameter (igaatb  = 2)
      parameter (igabta  = 3)
      parameter (igabtb  = 4)
      parameter (igbbta  = 5)
      parameter (igbbtb  = 6)
c
C *Basis sets include file
c
      real*8  alpha
      real*8  cont_coeff
      integer num_bset,bset_tags
      common/basis_sets/alpha(max_tag,max_atype,max_prm),
     &                  cont_coeff(max_tag,max_atype,max_prm,max_ang),
     &                  num_bset,bset_tags(max_tag)

      integer Ashl, Aprm, Abfn, totshl, totprm, totbfn, size_shlA
      integer size_basA, size_primA, maxi_shlA, maxi_basA, maxi_primA
      integer num_types, atom_tag , num_shl, atm_typ, nprim, angmom
      integer hybrid, pstart

      common/basis_size_info/Ashl(max_tag,max_atype),
     &                       Aprm(max_tag,max_atype),
     &                       Abfn(max_tag,max_atype),
     &                       totshl(max_tag),
     &                       totprm(max_tag),
     &                       totbfn(max_tag),
     &                       size_shlA,maxi_shlA,
     &                       size_basA,maxi_basA,
     &                       size_primA,maxi_primA
C
C Descriptions
C 
C  num_types		-	number of types of atoms for a given basis
C  atom_tag		-	list of atom type id numbers
C
      common/basis_cent_info/num_types(max_tag),
     &                       atom_tag(max_tag,max_atom)
C
C Descriptions
C
C  num_shl		-	number of shells for given atom type
C  atm_typ		-	atomic number of atom type
C  nprim		-  	number of primitives 
C  angmom		-	angular momentum
C  hybrid		-	level of hybrid shell. Is same as angmom if
C				shell is not a hybrid. For hybrid shells is
C				always less than angmom e.g. for an sp shell
C				angmom=2, hybrid=1
C  pstart		-	start of exponents and contraction coeffs 
C				contained in include basis.hf77, for a given
C				shell
C
      common/basis_cont_info/num_shl(max_tag,max_atype),
     &                       atm_typ(max_tag,max_atype),
     &                       nprim(max_tag,max_atype,max_shel),
     &                       angmom(max_tag,max_atype,max_shel),
     &                       hybrid(max_tag,max_atype,max_shel),
     &                       pstart(max_tag,max_atype,max_shel)

c
      logical gradcorr_sw
      logical gwt_sw
      logical gwt_avail_sw
c
      integer n_active_bfn
      integer n_active_atm
      integer active_bfn_list(n_active_bfn)
      integer active_bfn_indx(n_active_bfn)
      integer active_bfn_atms(n_active_atm)
      integer num_near_atoms
      integer near_atom_list(num_near_atoms)
c
      integer npts
      integer naocc
      integer nao
      integer ngridcent
      integer npert
      integer mxp
      integer ao_tag
      integer iatm
c
      real*8 wt(mxp)
      real*8 gwt(3,mxp,num_near_atoms)
      real*8 g2wt(mxp,3*num_near_atoms,3*num_near_atoms)
c
      real*8 bfn_val(mxp,nao)
      real*8 bfng_val(mxp,nao,3)
      real*8 bfn_hess(mxp,nao,6)
      real*8 bfn_3rd(mxp,nao,10)
c
      real*8 xc_ept(mxp)
      real*8 xc_vpt(mxp,2)
      real*8 xc_dvpt(mxp,3)
      real*8 xc_hpt(mxp,3)
      real*8 xc_dhpt(mxp,12)
c
      real*8 avec(nao,naocc)
      real*8 amo_val(npts,naocc)
      real*8 amo_grad(npts,naocc,3)
      real*8 grho(mxp,2,3)
c
c     Outputs
c
      real*8 hess(npert,npert)
c
c     Workspave
c
      real*8 drho(npts,3*n_active_atm)
      real*8 dgrho(npts,3,3*n_active_atm)
      real*8 dgaa(npts,3*n_active_atm)
      real*8 ddrhoa(npts,9)
      real*8 ddgrhoa(npts,3,9)
      real*8 bdmo(npts,3)
      real*8 admo(npts,3)
      real*8 bddmo(npts,6)
      real*8 bdgmo(npts,6)
      real*8 adgmo(npts,6)
      real*8 bddgmo(npts,10)
      real*8 t1(npts)
c
c     istart_bfn: For every atom where the basis functions start.
c                 Assumption is that all basis functions associated
c                 with 1 atom are consecutively located on the
c                 active_bfn_list.
c                 Beware: not every active atom necessarily has any
c                 active basis functions. (See below also.)
c
      integer istart_bfn(n_active_atm+1)
c
c     Local
c
      integer ipt
      integer katm, latm
      integer kpert, lpert, ipert
      integer kprt,  lprt
      integer kc, lc, mc
      integer klxx, kx
      integer mu, imu
      integer ith
      integer batmt, bstart, bend
      integer aatmt, astart, aend
c
      real*8 scr(9)
      real*8 gaadab
c
c     Code
c
      do katm=1,n_active_atm+1
         istart_bfn(katm) = 0
      enddo
      do ith=n_active_bfn,1,-1
         istart_bfn(active_bfn_indx(ith))=ith
      enddo
      istart_bfn(n_active_atm+1)=n_active_bfn+1
      do katm=n_active_atm,1,-1
         if (istart_bfn(katm).eq.0) 
     &       istart_bfn(katm) = istart_bfn(katm+1)
      enddo
c
      ipert=3*(iatm-1)
c
c     if (gwt_sw.and.gwt_avail_sw) then
c           do katm = 1, num_near_atoms
c              kprt = 3*(katm-1)
c              kpert = 3*(near_atom_list(katm)-1)
c              do latm = 1, num_near_atoms
c                 lprt = 3*(latm-1)
c                 lpert = 3*(near_atom_list(latm)-1)
c                 call aclear_dp(scr,9,0.0d0)
c                 do kc=1,3
c                    kx=3*(kc-1)
c                    do lc=1,3
c                       do ipt=1,npts
c                          scr(kx+lc)=scr(kx+lc)
c    &                     + g2wt(ipt,kprt+kc,lprt+lc)
c    &                       *xc_ept(ipt)
c                       enddo
c                    enddo
c                 enddo
c                 do kc=1,3
c                    kx=3*(kc-1)
c                    do lc=1,3
c                       hess(kpert+kc,lpert+lc)
c    &                  = hess(kpert+kc,lpert+lc) + scr(kx+lc)
c                    enddo
c                 enddo
c              enddo
c           enddo
c     endif
c
c     evaluate drho and dgrho
c
      call aclear_dp(drho,npts*3*n_active_atm,0.0d0)
      if (gradcorr_sw) then
         call aclear_dp(dgrho,3*npts*3*n_active_atm,0.0d0)
      endif
      do katm=1,n_active_atm
         kprt=3*(katm-1)
         bstart=istart_bfn(katm)
         bend  =istart_bfn(katm+1)-1
         do kc=1,3
            do ith=1,naocc
               call aclear_dp(t1,npts,0.0d0)
               do imu=bstart,bend
                  mu=active_bfn_list(imu)
                  do ipt=1,npts
                     t1(ipt)=t1(ipt)
     &               - avec(mu,ith)*bfng_val(ipt,imu,kc)
                  enddo ! ipt
               enddo ! imu
c
c              t1 = dmo(ipt,kprt+k,ith)
c
               do ipt=1,npts 
                  drho(ipt,kprt+kc)=drho(ipt,kprt+kc)
     &            + amo_val(ipt,ith)*t1(ipt)
               enddo ! ipt
c
               if (gradcorr_sw) then
                  do lc=1,3
                     do ipt=1,npts
                        dgrho(ipt,lc,kprt+kc)=dgrho(ipt,lc,kprt+kc)
     &                  + amo_grad(ipt,ith,lc)*t1(ipt)
                     enddo ! ipt
                  enddo ! lc
               endif ! gradcorr_sw
            enddo ! ith
         enddo ! kc
      enddo ! katm
c
      if (gradcorr_sw) then
         do katm=1,n_active_atm
            kprt=3*(katm-1)
            bstart=istart_bfn(katm)
            bend  =istart_bfn(katm+1)-1
            do kc=1,3
               do lc=1,3
                  do ith=1,naocc
                     call aclear_dp(t1,npts,0.0d0)
                     do imu=bstart,bend
                        mu=active_bfn_list(imu)
                        do ipt=1,npts
                           t1(ipt)=t1(ipt)
     &                     - avec(mu,ith)*bfn_hess(ipt,imu,hcc(lc,kc))
                        enddo ! ipt
                     enddo ! imu
c
c                    t1 = dgmo(ipt,lc,kprt+kc,ith)
c
                     do ipt=1,npts
                        dgrho(ipt,lc,kprt+kc)=dgrho(ipt,lc,kprt+kc)
     &                  + amo_val(ipt,ith)*t1(ipt)
                     enddo ! ipt
                  enddo ! ith
               enddo ! lc
            enddo ! kc
         enddo ! katm
c
         call aclear_dp(dgaa,npts*3*n_active_atm,0.0d0)
         do katm=1,n_active_atm
            kprt=3*(katm-1)
            do kc=1,3
               do lc=1,3
                  do ipt=1,npts
                     dgaa(ipt,kprt+kc)=dgaa(ipt,kprt+kc)
     &               + dgrho(ipt,lc,kprt+kc)*grho(ipt,1,lc)
                  enddo ! ipt
               enddo ! lc
            enddo ! kc
         enddo ! katm
      endif ! gradcorr_sw
c
c     Now do the terms with the gradients of the weights
c
      if (gwt_sw.and.gwt_avail_sw) then
         do 100 katm=1,n_active_atm
            kpert=3*(active_bfn_atms(katm)-1)
            if (kpert.eq.ipert) goto 100
            kprt=3*(katm-1)
            do 110 latm=1,num_near_atoms
               lpert=3*(near_atom_list(latm)-1)
c              lprt=3*(latm-1)
               call aclear_dp(scr,9,0.0d0)
               do kc=1,3
                  kx=3*(kc-1)
                  do lc=1,3
                     klxx=kx+lc
                     do ipt=1,npts
                        scr(klxx)=scr(klxx)
     &                  + 4*gwt(lc,ipt,latm)*xc_vpt(ipt,ira)*
     &                    drho(ipt,kprt+kc)
                     enddo ! ipt
                  enddo ! lc
               enddo ! kc
               if (gradcorr_sw) then
                  do kc=1,3
                     kx=3*(kc-1)
                     do lc=1,3
                        klxx=kx+lc
                        do ipt=1,npts
                           scr(klxx)=scr(klxx)
     &                     + 2*gwt(lc,ipt,latm)*xc_dvpt(ipt,igaa)*
     &                       dgaa(ipt,kprt+kc)
                        enddo ! ipt
                     enddo ! lc
                  enddo ! kc
               endif ! gradcorr_sw
               do kc=1,3
                  kx=3*(kc-1)
                  do lc=1,3
                     klxx=kx+lc
                     hess(lpert+lc,kpert+kc)
     &               =hess(lpert+lc,kpert+kc)+scr(klxx)
                     hess(lpert+lc,ipert+kc)
     &               =hess(lpert+lc,ipert+kc)-scr(klxx)
                     hess(kpert+kc,lpert+lc)
     &               =hess(kpert+kc,lpert+lc)+scr(klxx)
                     hess(ipert+kc,lpert+lc)
     &               =hess(ipert+kc,lpert+lc)-scr(klxx)
                  enddo ! lc
               enddo ! kc
 110        continue ! latm
 100     continue ! katm
      endif
c
      do 120 katm=1,n_active_atm
         kpert  = 3*(active_bfn_atms(katm)-1)
         if (gwt_sw.and.kpert.eq.ipert) goto 120
         bstart = istart_bfn(katm)
         bend   = istart_bfn(katm+1)-1
         kprt   = 3*(katm-1)
         do 130 latm=1,katm-1
            lpert  = 3*(active_bfn_atms(latm)-1)
            if (gwt_sw.and.lpert.eq.ipert) goto 130
            astart = istart_bfn(latm)
            aend   = istart_bfn(latm+1)-1
            lprt   = 3*(latm-1)
c
            call aclear_dp(ddrhoa,npts*9,0.0d0)
            if (gradcorr_sw) then
               call aclear_dp(ddgrhoa,npts*3*9,0.0d0)
            endif
            do ith=1,naocc
               call aclear_dp(admo,npts*3,0.0d0)
               call aclear_dp(bdmo,npts*3,0.0d0)
               do kc=1,3
                  do imu=bstart,bend
                     mu=active_bfn_list(imu)
                     do ipt=1,npts
                        bdmo(ipt,kc)=bdmo(ipt,kc)
     &                  -avec(mu,ith)*bfng_val(ipt,imu,kc)
                     enddo ! ipt
                  enddo ! imu
               enddo ! kc
               do kc=1,3
                  do imu=astart,aend
                     mu=active_bfn_list(imu)
                     do ipt=1,npts
                        admo(ipt,kc)=admo(ipt,kc)
     &                  -avec(mu,ith)*bfng_val(ipt,imu,kc)
                     enddo ! ipt
                  enddo ! imu
               enddo ! kc
               do kc=1,3
                  kx=3*(kc-1)
                  do lc=1,3
                     klxx=kx+lc
                     do ipt=1,npts
                        ddrhoa(ipt,klxx)=ddrhoa(ipt,klxx)
     &                  +bdmo(ipt,kc)*admo(ipt,lc)
                     enddo
                  enddo ! lc
               enddo ! kc
               if (gradcorr_sw) then
                  call aclear_dp(adgmo,npts*6,0.0d0)
                  call aclear_dp(bdgmo,npts*6,0.0d0)
                  do kx=1,6
                     do imu=bstart,bend
                        mu=active_bfn_list(imu)
                        do ipt=1,npts
                           bdgmo(ipt,kx)=bdgmo(ipt,kx)
     &                     -avec(mu,ith)*bfn_hess(ipt,imu,kx)
                        enddo ! ipt
                     enddo ! imu
                  enddo ! kx
                  do kx=1,6
                     do imu=astart,aend
                        mu=active_bfn_list(imu)
                        do ipt=1,npts
                           adgmo(ipt,kx)=adgmo(ipt,kx)
     &                     -avec(mu,ith)*bfn_hess(ipt,imu,kx)
                        enddo ! ipt
                     enddo ! imu
                  enddo ! kx
                  do kc=1,3
                     kx=3*(kc-1)
                     do lc=1,3
                        klxx=kx+lc
                        do mc=1,3
                           do ipt=1,npts
                              ddgrhoa(ipt,mc,klxx)=ddgrhoa(ipt,mc,klxx)
     &                        +bdgmo(ipt,hcc(mc,kc))*admo(ipt,lc)
     &                        +bdmo(ipt,kc)*adgmo(ipt,hcc(mc,lc))
                           enddo
                        enddo ! mc
                     enddo ! lc
                  enddo ! kc
               endif
            enddo ! ith
c
            call aclear_dp(scr,9,0.0d0)
            do kc=1,3
               kx=3*(kc-1)
               do lc=1,3
                  klxx = kx+lc
                  do ipt=1,npts
                     scr(klxx)=scr(klxx)
     &               + 8*wt(ipt)*xc_hpt(ipt,irara)*
     &                 drho(ipt,kprt+kc)*drho(ipt,lprt+lc)
     &               + 4*wt(ipt)*xc_vpt(ipt,ira)*ddrhoa(ipt,klxx)
                  enddo
                  if (gradcorr_sw) then
                     do ipt=1,npts
                        gaadab  
     &                  =   (ddgrhoa(ipt,1,klxx)*grho(ipt,1,1)
     &                      +ddgrhoa(ipt,2,klxx)*grho(ipt,1,2)
     &                      +ddgrhoa(ipt,3,klxx)*grho(ipt,1,3))
     &                  + 4*(dgrho(ipt,1,kprt+kc)*
     &                       dgrho(ipt,1,lprt+lc)
     &                      +dgrho(ipt,2,kprt+kc)*
     &                       dgrho(ipt,2,lprt+lc)
     &                      +dgrho(ipt,3,kprt+kc)*
     &                       dgrho(ipt,3,lprt+lc))
                        scr(kx+lc)=scr(kx+lc)
     &                  + 8*wt(ipt)*xc_dhpt(ipt,iragaa)*
     &                    (drho(ipt,kprt+kc)*dgaa(ipt,lprt+lc) 
     &                    +drho(ipt,lprt+lc)*dgaa(ipt,kprt+kc))
     &                  + 4*wt(ipt)*xc_dhpt(ipt,igaagaa)*
     &                    dgaa(ipt,lprt+lc)*dgaa(ipt,kprt+kc)
     &                  + 2*wt(ipt)*xc_dvpt(ipt,igaa)*gaadab
                     enddo
                  endif
               enddo
            enddo
            do kc=1,3
               kx=3*(kc-1)
               do lc=1,3
                  hess(kpert+kc,lpert+lc) = hess(kpert+kc,lpert+lc)
     &                                    + scr(kx+lc)
                  hess(lpert+lc,kpert+kc) = hess(lpert+lc,kpert+kc)
     &                                    + scr(kx+lc)
                  if (gwt_sw) then
                     hess(kpert+kc,ipert+lc) 
     &               = hess(kpert+kc,ipert+lc) - scr(kx+lc)
                     hess(ipert+lc,kpert+kc) 
     &               = hess(ipert+lc,kpert+kc) - scr(kx+lc)
c
                     hess(ipert+kc,lpert+lc) 
     &               = hess(ipert+kc,lpert+lc) - scr(kx+lc)
                     hess(lpert+lc,ipert+kc) 
     &               = hess(lpert+lc,ipert+kc) - scr(kx+lc)
c
                     hess(ipert+kc,ipert+lc) 
     &               = hess(ipert+kc,ipert+lc) + scr(kx+lc)
                     hess(ipert+lc,ipert+kc) 
     &               = hess(ipert+lc,ipert+kc) + scr(kx+lc)
                  endif
               enddo
            enddo
 130     continue
c
         call aclear_dp(ddrhoa,npts*6,0.0d0)
c        2nd derivative of 1 basis function
         do ith=1,naocc
            call aclear_dp(bddmo,npts*6,0.0d0)
            do kc=1,3
               kx=kc*(kc-1)/2
               do lc=1,kc
                  klxx=kx+lc
                  do imu=bstart,bend
                     mu=active_bfn_list(imu)
                     do ipt=1,npts
                        bddmo(ipt,klxx)=bddmo(ipt,klxx)
     &                  +avec(mu,ith)*bfn_hess(ipt,imu,klxx)
                     enddo ! ipt
                  enddo ! imu
               enddo ! lc
            enddo ! kc
            do kc=1,3
               kx=kc*(kc-1)/2
               do lc=1,kc
                  klxx=kx+lc
                  do ipt=1,npts
                     ddrhoa(ipt,klxx)=ddrhoa(ipt,klxx)
     &               +amo_val(ipt,ith)*bddmo(ipt,hcc(lc,kc))
                  enddo ! ipt
               enddo ! lc
            enddo ! kc
         enddo ! ith
c        1st derivative of 2 basis functions
         do ith=1,naocc
            call aclear_dp(bdmo,npts*3,0.0d0)
            do kc=1,3
               do imu=bstart,bend
                  mu=active_bfn_list(imu)
                  do ipt=1,npts
                     bdmo(ipt,kc)=bdmo(ipt,kc)
     &               -avec(mu,ith)*bfng_val(ipt,imu,kc)
                  enddo ! ipt
               enddo ! imu
            enddo ! kc
            do kc=1,3
               kx=kc*(kc-1)/2
               do lc=1,kc
                  klxx=kx+lc
                  do ipt=1,npts
                     ddrhoa(ipt,klxx)=ddrhoa(ipt,klxx)
     &               +bdmo(ipt,kc)*bdmo(ipt,lc)
                  enddo ! ipt
               enddo ! lc
            enddo ! kc
         enddo ! ith
c
         if (gradcorr_sw) then
            call aclear_dp(ddgrhoa,npts*3*6,0.0d0)
            do ith=1,naocc
c              2nd derivative of 1 basis function
               call aclear_dp(bddmo,npts*6,0.0d0)
               call aclear_dp(bddgmo,npts*10,0.0d0)
               do kx=1,6
                  do imu=bstart,bend
                     mu=active_bfn_list(imu)
                     do ipt=1,npts
                        bddmo(ipt,kx)=bddmo(ipt,kx)
     &                  +avec(mu,ith)*bfn_hess(ipt,imu,kx)
                     enddo ! ipt
                  enddo ! imu
               enddo ! kx
               do kx=1,10
                  do imu=bstart,bend
                     mu=active_bfn_list(imu)
                     do ipt=1,npts
                        bddgmo(ipt,kx)=bddgmo(ipt,kx)
     &                  +avec(mu,ith)*bfn_3rd(ipt,imu,kx)
                     enddo ! ipt
                  enddo ! imu
               enddo ! kx
               do kc=1,3
                  kx=kc*(kc-1)/2
                  do lc=1,kc
                     klxx=kx+lc
                     do mc=1,3
                        do ipt=1,npts
                           ddgrhoa(ipt,mc,klxx)=ddgrhoa(ipt,mc,klxx)
     &                     + amo_val(ipt,ith)*bddgmo(ipt,tc3(mc,lc,kc))
     &                     + amo_grad(ipt,ith,mc)*bddmo(ipt,hcc(lc,kc))
                        enddo 
                     enddo ! mc
                  enddo ! lc
               enddo ! kc
c              1st derivative of 2 basis functions
c              ... nasty ... use the fact that -*- = +
               call aclear_dp(bdmo,npts*3,0.0d0)
               do kc=1,3
                  do imu=bstart,bend
                     mu=active_bfn_list(imu)
                     do ipt=1,npts
                        bdmo(ipt,kc)=bdmo(ipt,kc)
     &                  +avec(mu,ith)*bfng_val(ipt,imu,kc)
                     enddo ! ipt
                  enddo ! imu
               enddo ! kc
               do kc=1,3
                  kx=kc*(kc-1)/2
                  do lc=1,kc
                     klxx=kx+lc
                     do mc=1,3
                        do ipt=1,npts
                           ddgrhoa(ipt,mc,klxx)=ddgrhoa(ipt,mc,klxx)
     &                     + bdmo(ipt,lc)*bddmo(ipt,hcc(mc,kc))
     &                     + bddmo(ipt,hcc(mc,lc))*bdmo(ipt,kc)
                        enddo ! ipt
                     enddo ! mc
                  enddo ! lc
               enddo ! kc
            enddo ! ith
         endif
c
         call aclear_dp(scr,6,0.0d0)
         do kc=1,3
            kx=kc*(kc-1)/2
            do lc=1,kc
               klxx = kx+lc
               do ipt=1,npts
                  scr(klxx)=scr(klxx)
     &            + 8*wt(ipt)*xc_hpt(ipt,irara)*
     &              drho(ipt,kprt+kc)*drho(ipt,kprt+lc)
     &            + 4*wt(ipt)*xc_vpt(ipt,ira)*ddrhoa(ipt,klxx)
               enddo ! ipt
               if (gradcorr_sw) then
                  do ipt=1,npts
c                    gaada = dgaa(ipt,kpert+kc)
c    &               =   (dgrho(ipt,1,1,kpert+kc)*grho(ipt,1,1)
c    &                  + dgrho(ipt,1,2,kpert+kc)*grho(ipt,1,2)
c    &                  + dgrho(ipt,1,3,kpert+kc)*grho(ipt,1,3))
c                    gaadb = dgaa(ipt,kpert+lc)
c    &               =   (dgrho(ipt,1,1,kpert+lc)*grho(ipt,1,1)
c    &                  + dgrho(ipt,1,2,kpert+lc)*grho(ipt,1,2)
c    &                  + dgrho(ipt,1,3,kpert+lc)*grho(ipt,1,3))
                     gaadab 
     &               =   (ddgrhoa(ipt,1,klxx)*grho(ipt,1,1)
     &                   +ddgrhoa(ipt,2,klxx)*grho(ipt,1,2)
     &                   +ddgrhoa(ipt,3,klxx)*grho(ipt,1,3))
     &               + 4*(dgrho(ipt,1,kprt+kc)*
     &                    dgrho(ipt,1,kprt+lc)
     &                   +dgrho(ipt,2,kprt+kc)*
     &                    dgrho(ipt,2,kprt+lc)
     &                   +dgrho(ipt,3,kprt+kc)*
     &                    dgrho(ipt,3,kprt+lc))
                     scr(klxx)=scr(klxx)
     &               + 8*wt(ipt)*xc_dhpt(ipt,iragaa)*
     &                 (drho(ipt,kprt+kc)*dgaa(ipt,kprt+lc)
     &                 +drho(ipt,kprt+lc)*dgaa(ipt,kprt+kc))
     &               + 4*wt(ipt)*xc_dhpt(ipt,igaagaa)*
     &                 dgaa(ipt,kprt+kc)*dgaa(ipt,kprt+lc)
     &               + 2*wt(ipt)*xc_dvpt(ipt,igaa)*gaadab
                  enddo ! ipt
               endif
            enddo ! lc
         enddo ! kc
         do kc=1,3
            kx=kc*(kc-1)/2
            do lc=1,kc-1
               hess(kpert+kc,kpert+lc) = hess(kpert+kc,kpert+lc)
     &                                 + scr(kx+lc)
               hess(kpert+lc,kpert+kc) = hess(kpert+lc,kpert+kc)
     &                                 + scr(kx+lc)
               if (gwt_sw) then
                  hess(kpert+kc,ipert+lc) 
     &            = hess(kpert+kc,ipert+lc) - scr(kx+lc)
                  hess(ipert+lc,kpert+kc) 
     &            = hess(ipert+lc,kpert+kc) - scr(kx+lc)
c
                  hess(ipert+kc,kpert+lc) 
     &            = hess(ipert+kc,kpert+lc) - scr(kx+lc)
                  hess(kpert+lc,ipert+kc) 
     &            = hess(kpert+lc,ipert+kc) - scr(kx+lc)
c
                  hess(ipert+kc,ipert+lc) 
     &            = hess(ipert+kc,ipert+lc) + scr(kx+lc)
                  hess(ipert+lc,ipert+kc) 
     &            = hess(ipert+lc,ipert+kc) + scr(kx+lc)
               endif
            enddo ! lc
            hess(kpert+kc,kpert+kc) = hess(kpert+kc,kpert+kc)
     &                              + scr(kx+kc)
            if (gwt_sw) then
               hess(kpert+kc,ipert+kc) = hess(kpert+kc,ipert+kc)
     &                                 - scr(kx+kc)
               hess(ipert+kc,kpert+kc) = hess(ipert+kc,kpert+kc)
     &                                 - scr(kx+kc)
               hess(ipert+kc,ipert+kc) = hess(ipert+kc,ipert+kc)
     &                                 + scr(kx+kc)
            endif
         enddo ! kc
 120  continue
c
      return
      end
